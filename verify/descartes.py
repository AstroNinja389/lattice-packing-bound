#!/usr/bin/env python3
"""[proven] Exact, independent Descartes certificate verifier for even dimensions."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

from flint import fmpq, fmpq_poly, fmpz


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

RATIONAL = re.compile(r"-?[0-9]+(?:/[1-9][0-9]*)?\Z", re.ASCII)
DECIMAL = re.compile(
    r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z",
    re.ASCII,
)
FIELDS = {"n", "K", "tau0", "E", "c", "Gam", "claims"}
MAX_DEPTH = 256
MAX_NODES = 200000
EXP_TERMS = 512


class Rejected(Exception):
    """A failed exact condition or an exhausted conservative proof budget."""


def require(condition, message):
    if not condition:
        raise Rejected(message)


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "duplicate JSON key: " + key)
        result[key] = value
    return result


def invalid_constant(value):
    raise Rejected("non-finite JSON constant: " + value)


def load_json(text):
    try:
        return json.loads(
            text, object_pairs_hook=unique_object, parse_constant=invalid_constant
        )
    except (ValueError, RecursionError) as error:
        raise Rejected("invalid JSON: " + str(error)) from error


def rational(value, name):
    require(type(value) is str and RATIONAL.fullmatch(value),
            name + " must be a rational string")
    pieces = value.split("/")
    return fmpq(int(pieces[0]), int(pieces[1]) if len(pieces) == 2 else 1)


def decimal(value, name):
    require(type(value) is str and DECIMAL.fullmatch(value),
            name + " must be a decimal string")
    pieces = re.split("[eE]", value)
    exponent = int(pieces[1]) if len(pieces) == 2 else 0
    mantissa = pieces[0]
    if "." in mantissa:
        whole, fractional = mantissa.split(".")
        exponent -= len(fractional)
        digits = whole + fractional
    else:
        digits = mantissa
    number = int(digits)
    require(number > 0, name + " must be positive")
    if exponent >= 0:
        return fmpq(number * 10**exponent)
    return fmpq(number, 10**(-exponent))


def validate(data):
    require(type(data) is dict and set(data) == FIELDS, "wrong certificate fields")
    n, degree = data["n"], data["K"]
    require(type(n) is int and n >= 2, "n must be an integer at least 2")
    require(type(degree) is int and degree >= 1, "K must be an integer at least 1")
    require(n % 2 == 0, "odd n is unsupported by this even-dimensional verifier")
    tau0 = rational(data["tau0"], "tau0")
    require(tau0 > 0, "tau0 must be positive")
    vectors = {}
    for key in ("E", "c", "Gam"):
        require(type(data[key]) is list, key + " must be a JSON array")
        vectors[key] = [rational(x, key + "[" + str(i) + "]")
                        for i, x in enumerate(data[key])]
    edges, coefficients, gamma = vectors["E"], vectors["c"], vectors["Gam"]
    require(len(edges) >= 3, "E must have at least three edges")
    bins = len(edges) - 1
    require(len(coefficients) == degree + 1, "c must have K+1 entries")
    require(len(gamma) == bins, "Gam must have one entry per bin")
    require(edges[0] == 1, "E[0] must equal 1")
    require(all(a < b for a, b in zip(edges, edges[1:])),
            "E must be strictly increasing")
    require(fmpq(2) in edges, "E must contain 2")
    split = edges.index(fmpq(2))
    require(0 < split < bins, "2 must be an interior edge")
    require(edges[-1] < 4, "last edge must be below 4")
    require(n * (1 - 2 / edges[-1])**2 < 1, "last long cap denominator is not positive")
    require(all(x >= 0 for x in gamma), "Gam must be nonnegative")
    claims = data["claims"]
    require(type(claims) is dict and set(claims) == {"density", "hermite"},
            "wrong claims fields")
    density = decimal(claims["density"], "claimed density")
    hermite = decimal(claims["hermite"], "claimed hermite")
    return n, degree, tau0, edges, coefficients, gamma, split, density, hermite


def laguerre_pair(n, coefficients):
    alpha = fmpq(n, 2) - 1
    previous = fmpq_poly([0])
    current = fmpq_poly([1])
    p, q = fmpq_poly([0]), fmpq_poly([0])
    for k, coefficient in enumerate(coefficients):
        p += coefficient * current
        q += ((-1)**k * coefficient) * current
        if k + 1 < len(coefficients):
            following = (fmpq_poly([2*k + alpha + 1, -1]) * current
                         - (k + alpha) * previous) / (k + 1)
            previous, current = current, following
    return p, q


def primitive(values):
    divisor = fmpz(0)
    for value in values:
        divisor = divisor.gcd(value)
    if divisor > 1:
        return [value // divisor for value in values]
    return values


def integer_coefficients(poly, degree):
    values = [poly[i] for i in range(degree + 1)]
    denominator = fmpz(1)
    for value in values:
        other = value.denominator
        denominator = denominator // denominator.gcd(other) * other
    return primitive([value.numerator * (denominator // value.denominator)
                      for value in values])


def variations(values):
    previous, result = 0, 0
    for value in values:
        sign = 1 if value > 0 else -1 if value < 0 else 0
        if sign:
            if previous and previous != sign:
                result += 1
            previous = sign
    return result


def bisect_coefficients(values, binomial):
    # [proven] These are positive multiples of the two child interval transforms.
    degree = len(values) - 1
    left, right = [], []
    for k in range(degree + 1):
        left.append(sum((values[i] * binomial[degree-i][k-i]
                         for i in range(k + 1)), fmpz(0)) * (1 << (degree-k)))
        right.append(sum((values[i] * binomial[i][k]
                          for i in range(k, degree + 1)), fmpz(0)) * (1 << k))
    return primitive(left), primitive(right)


def prove_transformed(values, name, max_depth=MAX_DEPTH, max_nodes=MAX_NODES):
    degree = len(values) - 1
    require(values[0] > 0 and values[-1] > 0,
            name + ": a finite or compactified endpoint is not strictly positive")
    binomial = [[math.comb(i, j) for j in range(i + 1)]
                for i in range(degree + 1)]
    pending = [(values, 0)]
    nodes = 0
    while pending:
        current, depth = pending.pop()
        nodes += 1
        require(nodes <= max_nodes, name + ": Descartes node budget exhausted")
        if variations(current) == 0:
            continue
        require(depth < max_depth, name + ": Descartes depth budget exhausted")
        require(sum(current, fmpz(0)) > 0,
                name + ": a bisection point is not strictly positive")
        left, right = bisect_coefficients(current, binomial)
        pending.append((right, depth + 1))
        pending.append((left, depth + 1))


def interval_transform(poly, lower, upper):
    degree = poly.degree()
    numerator = fmpq_poly([lower, upper])
    denominator = fmpq_poly([1, 1])
    transformed = fmpq_poly([poly[degree]])
    denominator_power = fmpq_poly([1])
    for k in range(degree - 1, -1, -1):
        denominator_power *= denominator
        transformed = transformed * numerator + poly[k] * denominator_power
    return integer_coefficients(transformed, degree)


def positive_interval(poly, lower, upper, name, **budgets):
    require(lower < upper, name + ": invalid interval")
    require(poly.degree() >= 0, name + ": zero polynomial")
    require(poly(lower) > 0 and poly(upper) > 0,
            name + ": an interval endpoint is not strictly positive")
    prove_transformed(interval_transform(poly, lower, upper), name, **budgets)


def positive_half_line(poly, lower, name, **budgets):
    degree = poly.degree()
    require(degree >= 0, name + ": zero polynomial")
    require(poly(lower) > 0 and poly[degree] > 0,
            name + ": nonpositive start or leading coefficient")
    shifted = fmpq_poly([poly[degree]])
    linear = fmpq_poly([lower, 1])
    for k in range(degree - 1, -1, -1):
        shifted = shifted * linear + poly[k]
    # [proven] f(lower+x) is the Descartes transform of
    # (1-u)^degree f(lower+u/(1-u)) on the compact interval [0,1].
    prove_transformed(integer_coefficients(shifted, degree), name, **budgets)


def exponential_upper(x):
    require(x >= 0, "exponential argument must be nonnegative")
    term, partial = fmpq(1), fmpq(1)
    for k in range(1, EXP_TERMS + 1):
        term = term * x / k
        partial += term
    return 1 / partial


def floor_rational(value):
    return value.numerator // value.denominator


def long_cap(n, edge):
    cosine = 1 - 2 / edge
    denominator = 1 - n * cosine**2
    require(denominator > 0, "invalid long cap denominator")
    return floor_rational(n * (1 - cosine**2) / denominator)


def short_cap(n, short_edge, long_edge):
    if short_edge + long_edge <= 4:
        return fmpz(0)
    cosine = 1 - 2 / long_edge
    u = (short_edge + long_edge - 4)**2 / (4 * short_edge * long_edge)
    require(u < 1, "invalid projected cap denominator")
    projected = (cosine + u) / (1 - u)
    bound = long_cap(n, long_edge) - 1
    if (n - 1) * projected**2 < 1:
        projected_bound = floor_rational(
            (n - 1) * (1 - projected**2) / (1 - (n - 1) * projected**2)
        )
        bound = min(bound, projected_bound)
    require(bound >= 0, "negative short cap bound")
    return bound


def kappa_bound(n, tau0, edges, gamma, split):
    bins = len(gamma)
    weights = [gamma[b] * exponential_upper(tau0 * edges[b] / 2)
               for b in range(bins)] + [fmpq(0)]
    for b in range(bins - 2, split - 1, -1):
        weights[b] = max(weights[b], weights[b + 1])
    differences = [weights[b] - weights[b + 1] for b in range(split, bins)]
    empty = sum((differences[b-split] * long_cap(n, edges[b+1])
                 for b in range(split, bins)), fmpq(0))
    candidates = [empty]
    for j in range(split):
        candidates.append(weights[j] + sum(
            (differences[b-split] * short_cap(n, edges[j+1], edges[b+1])
             for b in range(split, bins)), fmpq(0)))
    return max(candidates)


def atan_enclosure(reciprocal):
    x = fmpq(1, reciprocal)
    partial, power = fmpq(0), x
    upper = None
    for k in range(80):
        partial += ((-1)**k * power) / (2*k + 1)
        if k == 78:
            upper = partial
        power *= x*x
    lower = partial
    require(0 < lower < upper, "invalid alternating arctangent enclosure")
    return lower, upper


def pi_lower_bound():
    lower5, upper5 = atan_enclosure(5)
    lower239, upper239 = atan_enclosure(239)
    lower = 16 * lower5 - 4 * upper239
    upper = 16 * upper5 - 4 * lower239
    # [proven] The rational tangent computation is the algebra in Machin's identity.
    tangent2 = 2*fmpq(1, 5)/(1-fmpq(1, 5)**2)
    tangent4 = 2*tangent2/(1-tangent2**2)
    require((tangent4-fmpq(1, 239))/(1+tangent4*fmpq(1, 239)) == 1,
            "Machin identity algebra failed")
    require(3 < lower < upper < 4, "invalid Machin enclosure")
    return lower


def verify(data):
    n, degree, tau0, edges, coefficients, gamma, split, density, hermite = validate(data)
    p, q = laguerre_pair(n, coefficients)
    positive_half_line(q, fmpq(0), "q on [0,infinity)")
    positive_half_line(-p, tau0*edges[-1], "-p on the tail")
    for b, magnitude in enumerate(gamma):
        lower, upper = tau0*edges[b], tau0*edges[b+1]
        taylor = fmpq_poly([1-lower/2+lower**2/8, fmpq(1, 2)-lower/4, fmpq(1, 8)])
        positive_interval(magnitude*taylor - p, lower, upper,
                          "bin " + str(b) + " majorant")
    kappa = kappa_bound(n, tau0, edges, gamma, split)
    rho = (p(0) + 2*((1 << n)-1)*kappa) / q(0)
    require(rho > 0, "rho must be positive")
    pi_lower = pi_lower_bound()
    half_dimension = n // 2
    density_upper = (tau0**half_dimension * rho
                     / (8**half_dimension * math.factorial(half_dimension)))
    require(density >= density_upper, "claimed density is below the rational upper bound")
    require(hermite**n >= (tau0/(2*pi_lower))**n * rho**2,
            "claimed hermite is below the rational sufficient bound")
    return "PASS n={} K={} density={} hermite={} [proven]".format(
        n, degree, data["claims"]["density"], data["claims"]["hermite"]
    )


def self_test():
    positive_interval(fmpq_poly([fmpq(101, 100), -2, 1]), fmpq(0), fmpq(2),
                      "positive quadratic")
    positive_half_line(fmpq_poly([2, -2, 1]), fmpq(0), "positive half-line")
    positive_half_line(fmpq_poly([7]), fmpq(100), "constant half-line")
    checks = 3

    def expect_rejection(call):
        nonlocal checks
        try:
            call()
        except Rejected:
            checks += 1
        else:
            raise AssertionError("self-test expected rejection")

    expect_rejection(lambda: positive_interval(fmpq_poly([1, -2, 1]), fmpq(0),
                                               fmpq(2), "rational zero"))
    expect_rejection(lambda: positive_interval(fmpq_poly([4, 0, -4, 0, 1]), fmpq(1),
                                               fmpq(2), "irrational zero", max_depth=12))
    expect_rejection(lambda: positive_half_line(fmpq_poly([9999, -200, 1]), fmpq(0),
                                                "negative pocket"))
    expect_rejection(lambda: positive_half_line(fmpq_poly([0, 1]), fmpq(0),
                                                "zero endpoint"))
    expect_rejection(lambda: load_json('{"n":2,"n":2}'))
    expect_rejection(lambda: load_json('{"n":NaN}'))
    expect_rejection(lambda: rational("1/01", "test rational"))
    expect_rejection(lambda: decimal("0e100", "test decimal"))
    require(decimal("+.125e+1", "test decimal") == fmpq(5, 4), "decimal self-test")
    checks += 1
    certificate = {"n": 2, "K": 1, "tau0": "4", "E": ["1", "2", "3"],
                   "c": ["1", "1/2"], "Gam": ["0", "0"],
                   "claims": {"density": "2", "hermite": "2"}}
    require(verify(certificate).startswith("PASS "), "acceptance self-test")
    checks += 1
    for change in ({"n": 3}, {"n": True}, {"tau0": "0"}, {"extra": "1"},
                   {"claims": {"density": "1.49", "hermite": "2"}},
                   {"claims": {"density": "2", "hermite": "1.9"}}):
        expect_rejection(lambda change=change: verify(dict(certificate, **change)))
    print("SELF-TEST PASS checks={} [computed]".format(checks))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", nargs="?", type=Path)
    parser.add_argument("--self-test", action="store_true", help="run exact internal checks")
    args = parser.parse_args(argv)
    if args.self_test:
        if args.certificate is not None:
            parser.error("--self-test does not take a certificate path")
        try:
            self_test()
            return 0
        except Exception as error:
            print("SELF-TEST REJECT: {} [computed]".format(error), file=sys.stderr)
            return 1
    if args.certificate is None:
        parser.error("a certificate path is required")
    try:
        result = verify(load_json(args.certificate.read_text(encoding="utf-8")))
    except (Rejected, OSError, UnicodeError, ValueError, TypeError, OverflowError,
            ZeroDivisionError, RecursionError, MemoryError) as error:
        print("REJECT: {} [computed]".format(error), file=sys.stderr)
        return 1
    print(result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
