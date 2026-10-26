#!/usr/bin/env python3
"""[proven] Standalone exact integer Taylor certificate verifier."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import comb, factorial, gcd, lcm
from pathlib import Path
import re
import sys


EXP_TERMS = 512
MAX_DEPTH = 40
MAX_NODES = 200_000
RATIONAL = re.compile(r"-?[0-9]+(?:/[1-9][0-9]*)?\Z")
DECIMAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")
FIELDS = {"n", "K", "tau0", "E", "c", "Gam", "claims"}


class Rejected(Exception):
    pass


def require(condition, message):
    if not condition:
        raise Rejected(message)


def rational(value, name):
    require(type(value) is str and RATIONAL.fullmatch(value) is not None,
            f"{name}: invalid rational string")
    numerator, slash, denominator = value.partition("/")
    return Fraction(int(numerator), int(denominator) if slash else 1)


def decimal(value, name):
    require(type(value) is str and DECIMAL.fullmatch(value) is not None,
            f"{name}: invalid decimal string")
    result = Fraction(value)
    require(result > 0, f"{name}: claim must be positive")
    return result


def unique_object(items):
    result = {}
    for key, value in items:
        require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def forbidden_json_number(value):
    raise Rejected("JSON noninteger number is forbidden")


def parse_certificate(source):
    obj = json.loads(source, object_pairs_hook=unique_object,
                     parse_constant=forbidden_json_number,
                     parse_float=forbidden_json_number)
    require(type(obj) is dict and set(obj) == FIELDS,
            "certificate must have exactly the seven required fields")
    n, degree = obj["n"], obj["K"]
    require(type(n) is int and n >= 2, "n must be an integer at least 2")
    require(type(degree) is int and degree >= 1,
            "K must be an integer at least 1")
    tau0 = rational(obj["tau0"], "tau0")
    require(tau0 > 0, "tau0 must be positive")
    for key in ("E", "c", "Gam"):
        require(type(obj[key]) is list, f"{key} must be an array")
    edges = [rational(x, f"E[{i}]") for i, x in enumerate(obj["E"])]
    coefficients = [rational(x, f"c[{i}]") for i, x in enumerate(obj["c"])]
    gamma = [rational(x, f"Gam[{i}]") for i, x in enumerate(obj["Gam"])]
    require(len(edges) >= 3, "E must have at least three edges")
    require(edges[0] == 1, "E[0] must equal 1")
    require(all(a < b for a, b in zip(edges, edges[1:])),
            "E must be strictly increasing")
    require(Fraction(2) in edges, "E must include 2")
    split = edges.index(Fraction(2))
    bins = len(edges) - 1
    require(0 < split < bins, "the edge 2 must be internal")
    require(edges[-1] < 4, "last edge must be less than 4")
    require(n * (1 - 2 / edges[-1]) ** 2 < 1,
            "last edge violates the strict cap condition")
    require(len(coefficients) == degree + 1, "c must have K+1 coefficients")
    require(len(gamma) == bins and all(x >= 0 for x in gamma),
            "Gam must have one nonnegative value per bin")
    claims = obj["claims"]
    require(type(claims) is dict and set(claims) == {"density", "hermite"},
            "claims must have exactly density and hermite")
    density = decimal(claims["density"], "density")
    hermite = decimal(claims["hermite"], "hermite")
    return n, degree, tau0, edges, coefficients, gamma, split, density, hermite, claims


def trim(polynomial):
    result = list(polynomial)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result or [0]


def primitive(polynomial):
    result = trim(polynomial)
    common = 0
    for value in result:
        common = gcd(common, abs(value))
    return [x // common for x in result] if common > 1 else result


def integer_polynomial(polynomial):
    denominator = 1
    for value in polynomial:
        denominator = lcm(denominator, Fraction(value).denominator)
    return primitive([int(Fraction(x) * denominator) for x in polynomial])


def evaluate(polynomial, argument):
    result = 0
    for coefficient in reversed(polynomial):
        result = result * argument + coefficient
    return result


def laguerre_pair(n, coefficients):
    alpha = Fraction(n, 2) - 1
    p = [Fraction(0)] * len(coefficients)
    q = p.copy()
    previous = [Fraction(0)]
    current = [Fraction(1)]
    for k, weight in enumerate(coefficients):
        for j, coefficient in enumerate(current):
            p[j] += weight * coefficient
            q[j] += (-1) ** k * weight * coefficient
        if k + 1 == len(coefficients):
            break
        following = [Fraction(0)] * (len(current) + 1)
        for j, coefficient in enumerate(current):
            following[j] += (2 * k + alpha + 1) * coefficient
            following[j + 1] -= coefficient
        for j, coefficient in enumerate(previous):
            following[j] -= (k + alpha) * coefficient
        previous, current = current, [x / (k + 1) for x in following]
    return trim(p), trim(q)


def affine_integer(polynomial, origin, width):
    """[proven] Return a positive integer multiple of P(origin+width*x)."""
    origin, width = Fraction(origin), Fraction(width)
    denominator = lcm(origin.denominator, width.denominator)
    a, b = int(origin * denominator), int(width * denominator)
    result = [polynomial[-1]]
    scale = denominator
    for coefficient in reversed(polynomial[:-1]):
        following = [0] * (len(result) + 1)
        for j, value in enumerate(result):
            following[j] += a * value
            following[j + 1] += b * value
        following[0] += coefficient * scale
        scale *= denominator
        result = following
    return primitive(result)


def midpoint_expansion(polynomial):
    """[proven] Return 2**d * P((1+z)/2), in ascending powers of z."""
    degree = len(polynomial) - 1
    result = [0] * (degree + 1)
    for i, coefficient in enumerate(polynomial):
        weighted = coefficient << (degree - i)
        for j in range(i + 1):
            result[j] += weighted * comb(i, j)
    return result


def certify_interval(polynomial, left, right, strict, max_depth=MAX_DEPTH,
                     max_nodes=MAX_NODES):
    """[proven] Accept only an exact Taylor lower-bound covering proof."""
    require(left <= right, "invalid proof interval")
    if left == right:
        value = evaluate(polynomial, left)
        require(value > 0 if strict else value >= 0,
                "sign condition fails on a singleton")
        return
    stack = [(affine_integer(polynomial, left, right - left), 0)]
    visited = 0
    while stack:
        current, depth = stack.pop()
        visited += 1
        require(visited <= max_nodes, "Taylor proof exhausted its node budget")
        endpoints = (current[0], sum(current))
        require(all(x > 0 if strict else x >= 0 for x in endpoints),
                "polynomial violates the required endpoint sign")
        lower = current[0] + sum(min(0, x) for x in current[1:])
        if lower > 0 if strict else lower >= 0:
            continue
        centered = midpoint_expansion(current)
        require(centered[0] > 0 if strict else centered[0] >= 0,
                "polynomial violates the required midpoint sign")
        lower = centered[0]
        for power, coefficient in enumerate(centered[1:], 1):
            lower += -abs(coefficient) if power % 2 else min(0, coefficient)
        if lower > 0 if strict else lower >= 0:
            continue
        require(depth < max_depth,
                "Taylor proof exhausted its depth budget; sign is unresolved")
        degree = len(current) - 1
        left_child = [x << (degree - i) for i, x in enumerate(current)]
        stack.append((primitive(centered), depth + 1))
        stack.append((primitive(left_child), depth + 1))


def ceiling(value):
    value = Fraction(value)
    return -((-value.numerator) // value.denominator)


def certify_positive_tail(polynomial, start, max_depth=MAX_DEPTH):
    """[proven] Cover a finite prefix and an exact nonnegative Taylor tail."""
    polynomial = primitive(polynomial)
    require(polynomial[-1] > 0, "tail has nonpositive leading coefficient")
    shifted = affine_integer(polynomial, start, 1)
    require(shifted[0] > 0, "tail is not strictly positive at its left endpoint")
    if all(x >= 0 for x in shifted[1:]):
        return
    lower_size = max((abs(x) for x in polynomial[:-1]), default=0)
    root_bound = 1 + ceiling(Fraction(lower_size, polynomial[-1]))
    boundary = max(1, ceiling(start))
    limit = max(boundary, root_bound)
    while True:
        shifted = affine_integer(polynomial, boundary, 1)
        if shifted[0] > 0 and all(x >= 0 for x in shifted[1:]):
            certify_interval(polynomial, start, Fraction(boundary), True, max_depth)
            return
        require(boundary < limit,
                "positive Taylor tail could not be proved at the root bound")
        boundary = min(2 * boundary, limit)


def exponential_upper_reciprocal(argument):
    require(argument >= 0, "exponential argument must be nonnegative")
    term = Fraction(1)
    total = term
    for k in range(1, EXP_TERMS + 1):
        term *= argument / k
        total += term
    return 1 / total


def long_cap(n, squared_edge):
    cosine = 1 - 2 / squared_edge
    denominator = 1 - n * cosine ** 2
    require(denominator > 0, "long cap denominator must be positive")
    ratio = n * (1 - cosine ** 2) / denominator
    return ratio.numerator // ratio.denominator


def short_cap(n, short_edge, long_edge):
    if short_edge + long_edge <= 4:
        return 0
    cosine = 1 - 2 / long_edge
    u = (short_edge + long_edge - 4) ** 2 / (4 * short_edge * long_edge)
    require(u < 1, "short cap denominator must be positive")
    reduced = (cosine + u) / (1 - u)
    result = long_cap(n, long_edge) - 1
    if (n - 1) * reduced ** 2 < 1:
        ratio = (n - 1) * (1 - reduced ** 2) / (1 - (n - 1) * reduced ** 2)
        result = min(result, ratio.numerator // ratio.denominator)
    return result


def arctangent_partial(inverse, terms):
    argument = Fraction(1, inverse)
    power = argument
    total = Fraction(0)
    for k in range(terms):
        total += power / (2 * k + 1) if k % 2 == 0 else -power / (2 * k + 1)
        power *= argument ** 2
    return total


def pi_lower_bound():
    """[proven] Machin plus alternating sums proves the returned lower bound."""
    tangent = Fraction(1, 5)
    for _ in range(2):
        tangent = 2 * tangent / (1 - tangent ** 2)
    beta_tangent = Fraction(1, 239)
    require((tangent - beta_tangent) / (1 + tangent * beta_tangent) == 1,
            "Machin rational identity check failed")
    lower = 16 * arctangent_partial(5, 80) - 4 * arctangent_partial(239, 79)
    require(lower > 0, "proved pi lower bound must be positive")
    return lower


def verify(source, max_depth=MAX_DEPTH):
    n, degree, tau0, edges, coefficients, gamma, split, density, hermite, claims = parse_certificate(source)
    p, q = laguerre_pair(n, coefficients)
    tau = [tau0 * edge for edge in edges]
    certify_positive_tail(integer_polynomial(q), Fraction(0), max_depth)
    certify_positive_tail(integer_polynomial([-x for x in p]), tau[-1], max_depth)
    bins = len(gamma)
    for b in range(bins):
        gap = [-x for x in p] + [Fraction(0)] * max(0, 3 - len(p))
        gap[0] += gamma[b] * (1 - tau[b] / 2 + tau[b] ** 2 / 8)
        gap[1] += gamma[b] * (Fraction(1, 2) - tau[b] / 4)
        gap[2] += gamma[b] / 8
        certify_interval(integer_polynomial(gap), tau[b], tau[b + 1], False, max_depth)
    g = [gamma[b] * exponential_upper_reciprocal(tau[b] / 2) for b in range(bins)]
    g.append(Fraction(0))
    for b in range(bins - 2, split - 1, -1):
        g[b] = max(g[b], g[b + 1])
    empty = sum((g[b] - g[b + 1]) * long_cap(n, edges[b + 1])
                for b in range(split, bins))
    kappas = [empty]
    for j in range(split):
        kappas.append(g[j] + sum((g[b] - g[b + 1]) * short_cap(n, edges[j + 1], edges[b + 1])
                                  for b in range(split, bins)))
    kappa = max(kappas)
    rho = (p[0] + 2 * (2 ** n - 1) * kappa) / q[0]
    require(rho > 0, "rho must be positive")
    pi_lower = pi_lower_bound()
    half_dimension = n // 2
    if n % 2 == 0:
        bound = tau0 ** half_dimension * rho / (8 ** half_dimension * factorial(half_dimension))
        require(density >= bound, "claimed density is below the proved rational bound")
    else:
        odd_product = 1
        for factor in range(1, n + 1, 2):
            odd_product *= factor
        prefactor = 2 * tau0 ** half_dimension * rho / (4 ** half_dimension * odd_product)
        require(density ** 2 >= prefactor ** 2 * tau0 / (8 * pi_lower),
                "claimed density fails the exact squared inequality")
    require(hermite ** n >= (tau0 / (2 * pi_lower)) ** n * rho ** 2,
            "claimed Hermite bound fails the exact power inequality")
    return f"PASS [proven] n={n} K={degree} density={claims['density']} hermite={claims['hermite']}"


def self_test():
    count = 0

    def check(condition):
        nonlocal count
        require(condition, "internal self-test failed")
        count += 1

    def fails(action):
        nonlocal count
        try:
            action()
        except (Rejected, ValueError):
            count += 1
            return
        raise Rejected("internal rejection self-test failed")

    p, q = laguerre_pair(3, [Fraction(0), Fraction(0), Fraction(1)])
    check(p == q == [Fraction(15, 8), Fraction(-5, 2), Fraction(1, 2)])
    check(affine_integer([1, -2, 3], Fraction(1, 3), Fraction(2, 5)) == [25, 0, 18])
    check(decimal("+.125e+2", "test") == Fraction(25, 2))
    for polynomial, left, right, strict in (
        ([1370, -7400, 10000], Fraction(0), Fraction(1), True),
        ([0, 1, -1], Fraction(0), Fraction(1), False),
        ([1, -4, 4], Fraction(0), Fraction(1), False),
        ([7], Fraction(-2), Fraction(8), True),
    ):
        certify_interval(polynomial, left, right, strict)
        count += 1
    fails(lambda: certify_interval([3, -16, 16], Fraction(0), Fraction(1), True))
    fails(lambda: certify_interval([1, -4, 4], Fraction(0), Fraction(1), True))
    fails(lambda: certify_interval([1370, -7400, 10000], Fraction(0), Fraction(1), True, 0))
    fails(lambda: certify_interval([1, 0, -4, 0, 4], Fraction(0), Fraction(1), False, 0))
    certify_positive_tail([1, -2, 2], Fraction(0))
    count += 1
    fails(lambda: certify_positive_tail([1, -1], Fraction(0)))
    fails(lambda: certify_positive_tail([0, 1], Fraction(0)))
    check(Fraction(3) < pi_lower_bound() < Fraction(22, 7))
    check(exponential_upper_reciprocal(Fraction(0)) == 1)
    base = {"n": 2, "K": 1, "tau0": "2", "E": ["1", "2", "3"],
            "c": ["2", "1"], "Gam": ["1", "0"],
            "claims": {"density": "2", "hermite": "2"}}
    check(verify(json.dumps(base)).startswith("PASS [proven] n=2 K=1"))
    odd = dict(base, n=3, c=["3", "1"], Gam=["3", "1"],
               claims={"density": "10000", "hermite": "10000"})
    check(verify(json.dumps(odd)).startswith("PASS [proven] n=3 K=1"))
    for key, value in (("n", True), ("K", 0), ("tau0", "1/0"),
                       ("E", ["1", "2", "2"]), ("Gam", ["1", "-1"]),
                       ("claims", {"density": "0", "hermite": "2"}),
                       ("claims", {"density": "1e-20", "hermite": "2"}),
                       ("claims", {"density": "2", "hermite": "1e-20"}),
                       ("claims", {"density": "-2", "hermite": "2"}),
                       ("claims", {"density": 2, "hermite": "2"}),
                       ("c", ["2", "1", "0"])):
        damaged = dict(base)
        damaged[key] = value
        fails(lambda damaged=damaged: verify(json.dumps(damaged)))
    fails(lambda: parse_certificate('{"n":2,"n":3}'))
    fails(lambda: parse_certificate('{"claims":{"density":"1","density":"2"}}'))
    fails(lambda: parse_certificate('{"n":NaN}'))
    fails(lambda: parse_certificate('{"n":Infinity}'))
    fails(lambda: parse_certificate('{"n":2.0}'))
    fails(lambda: parse_certificate(json.dumps(dict(base, extra="1"))))
    return count


def main(argv=None):
    parser = argparse.ArgumentParser(description="Exact integer Taylor certificate verification")
    parser.add_argument("certificate", nargs="?", type=Path)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--max-depth", type=int, default=MAX_DEPTH)
    args = parser.parse_args(argv)
    if not 0 <= args.max_depth <= MAX_DEPTH:
        parser.error("--max-depth must lie between 0 and 40")
    if args.self_test == (args.certificate is not None):
        parser.error("provide either a certificate path or --self-test")
    try:
        if args.self_test:
            print(f"[computed] SELF-TEST PASS checks={self_test()}")
        else:
            print(verify(args.certificate.read_text(encoding="utf-8"), args.max_depth))
        return 0
    except (Rejected, ValueError, OSError, ArithmeticError, MemoryError, RecursionError) as exc:
        print(f"[computed] REJECT: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
