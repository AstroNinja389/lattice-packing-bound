#!/usr/bin/env python3
"""[proven] Conservative certificate verification with exact Bernstein bounds."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from math import comb, factorial, lcm
import json
from pathlib import Path
import re
import sys


MAX_DEPTH = 40
EXP_TERMS = 512
ATAN_TERMS = 80
RATIONAL_TOKEN = re.compile(r"-?[0-9]+(?:/[1-9][0-9]*)?\Z", re.ASCII)
DECIMAL_TOKEN = re.compile(
    r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z",
    re.ASCII,
)
ZERO = Fraction(0)
ONE = Fraction(1)


class Reject(ValueError):
    """A malformed certificate or an obligation that was not proved."""


@dataclass(frozen=True)
class Certificate:
    n: int
    K: int
    tau0: Fraction
    E: tuple[Fraction, ...]
    c: tuple[Fraction, ...]
    Gam: tuple[Fraction, ...]
    J: int
    density: Fraction
    hermite: Fraction
    density_text: str
    hermite_text: str


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Reject(message)


def unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        require(key not in result, "duplicate JSON object key")
        result[key] = value
    return result


def forbidden_number(_token: str) -> object:
    raise Reject("JSON decimal numbers, NaN, and Infinity are forbidden")


def rational(value: object, name: str) -> Fraction:
    require(isinstance(value, str), name + " must be a rational string")
    require(RATIONAL_TOKEN.fullmatch(value) is not None,
            name + " has invalid rational syntax")
    return Fraction(value)


def decimal(value: object, name: str) -> Fraction:
    require(isinstance(value, str), name + " must be a decimal string")
    require(DECIMAL_TOKEN.fullmatch(value) is not None,
            name + " has invalid decimal syntax")
    # [proven] Fraction parses finite decimal and scientific notation exactly.
    result = Fraction(value)
    require(result > 0, name + " must be positive")
    return result


def rational_array(value: object, name: str) -> tuple[Fraction, ...]:
    require(type(value) is list, name + " must be a JSON array")
    return tuple(rational(item, name + " entry") for item in value)


def parse_certificate(text: str) -> Certificate:
    try:
        obj = json.loads(text, object_pairs_hook=unique_object,
                         parse_float=forbidden_number,
                         parse_constant=forbidden_number)
    except json.JSONDecodeError as error:
        raise Reject("invalid JSON") from error
    require(type(obj) is dict, "certificate must be a JSON object")
    require(set(obj) == {"n", "K", "tau0", "E", "c", "Gam", "claims"},
            "certificate fields do not match the contract")
    n, K = obj["n"], obj["K"]
    require(type(n) is int and n >= 2, "n must be an integer at least 2")
    require(type(K) is int and K >= 1, "K must be an integer at least 1")
    tau0 = rational(obj["tau0"], "tau0")
    E = rational_array(obj["E"], "E")
    c = rational_array(obj["c"], "c")
    Gam = rational_array(obj["Gam"], "Gam")
    require(tau0 > 0, "tau0 must be positive")
    require(len(E) >= 3, "E must contain at least three squared edges")
    require(E[0] == 1, "the first squared edge must be 1")
    require(all(left < right for left, right in zip(E, E[1:])),
            "squared edges must be strictly increasing")
    require(Fraction(2) in E, "squared edges must contain 2")
    J = E.index(Fraction(2))
    M = len(E) - 1
    require(0 < J < M, "squared edge 2 must be an interior edge")
    require(E[-1] < 4, "the last squared edge must be less than 4")
    require(n * (1 - 2 / E[-1]) ** 2 < 1,
            "the last squared edge violates the cap denominator condition")
    require(len(c) == K + 1, "c must have K+1 entries")
    require(len(Gam) == M, "Gam must have one entry per bin")
    require(all(value >= 0 for value in Gam), "Gam must be nonnegative")
    claims = obj["claims"]
    require(type(claims) is dict and set(claims) == {"density", "hermite"},
            "claims must contain exactly density and hermite")
    density = decimal(claims["density"], "claimed density")
    hermite = decimal(claims["hermite"], "claimed hermite")
    return Certificate(n, K, tau0, E, c, Gam, J, density, hermite,
                       claims["density"], claims["hermite"])


def trim(poly: list[Fraction]) -> list[Fraction]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly or [ZERO]


def poly_add(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    result = [ZERO] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] += value
    for index, value in enumerate(right):
        result[index] += value
    return trim(result)


def poly_scale(poly: list[Fraction], scalar: Fraction) -> list[Fraction]:
    return trim([scalar * value for value in poly])


def multiply_linear(poly: list[Fraction], a: Fraction,
                    b: Fraction) -> list[Fraction]:
    result = [ZERO] * (len(poly) + 1)
    for index, value in enumerate(poly):
        result[index] += a * value
        result[index + 1] += b * value
    return trim(result)


def compose_affine(poly: list[Fraction], a: Fraction,
                   b: Fraction) -> list[Fraction]:
    """[proven] Return coefficients of poly(a+b*x), in ascending powers."""
    result = [ZERO]
    for value in reversed(poly):
        result = multiply_linear(result, a, b)
        result[0] += value
    return trim(result)


def laguerre_pair(cert: Certificate) -> tuple[list[Fraction], list[Fraction]]:
    alpha = Fraction(cert.n, 2) - 1
    previous, current = [ZERO], [ONE]
    p, q = [ZERO], [ZERO]
    for k, coefficient in enumerate(cert.c):
        p = poly_add(p, poly_scale(current, coefficient))
        q = poly_add(q, poly_scale(current, coefficient if k % 2 == 0
                                  else -coefficient))
        if k < cert.K:
            numerator = poly_add(
                multiply_linear(current, 2 * k + alpha + 1, -ONE),
                poly_scale(previous, -(k + alpha)),
            )
            following = poly_scale(numerator, Fraction(1, k + 1))
            previous, current = current, following
    return p, q


def interval_bernstein(poly: list[Fraction], left: Fraction,
                       right: Fraction) -> list[Fraction]:
    """[proven] Convert poly(left+(right-left)*x) to the Bernstein basis."""
    require(left < right, "finite sign interval must have positive length")
    power = compose_affine(trim(list(poly)), left, right - left)
    degree = len(power) - 1
    return [sum((power[j] * Fraction(comb(i, j), comb(degree, j))
                 for j in range(i + 1)), ZERO)
            for i in range(degree + 1)]


def halfline_bernstein(poly: list[Fraction], start: Fraction) -> list[Fraction]:
    """[proven] Return Bernstein coefficients after t=start+x/(1-x)."""
    power = compose_affine(trim(list(poly)), start, ONE)
    degree = len(power) - 1
    return [value / comb(degree, index) for index, value in enumerate(power)]


def bisect_bernstein(coefficients: list[Fraction]) -> tuple[list[Fraction],
                                                         list[Fraction]]:
    """[proven] Exact de Casteljau subdivision at x=1/2."""
    row = list(coefficients)
    left, right = [row[0]], [row[-1]]
    while len(row) > 1:
        row = [(row[index] + row[index + 1]) / 2
               for index in range(len(row) - 1)]
        left.append(row[0])
        right.append(row[-1])
    right.reverse()
    return left, right


def integer_scale(values) -> list[int]:
    """Multiply rational coefficients by one positive common denominator."""
    values = [Fraction(value) for value in values]
    common = lcm(*(value.denominator for value in values))
    return [value.numerator * (common // value.denominator) for value in values]


def integer_affine(poly, origin: Fraction, width: Fraction) -> list[int]:
    """Coefficients of D*Q^d*p(origin+width*x), with positive integers D,Q."""
    values = integer_scale(trim(list(poly)))
    degree = len(values) - 1
    Q = lcm(origin.denominator, width.denominator)
    A, B = int(origin * Q), int(width * Q)
    result = [values[-1]]
    power = Q
    for index in range(degree - 1, -1, -1):
        old = result
        result = [A * old[0]]
        result.extend(A * old[j] + B * old[j - 1] for j in range(1, len(old)))
        result.append(B * old[-1])
        result[0] += values[index] * power
        power *= Q
    return result


def integer_interval_bernstein(poly, left: Fraction, right: Fraction) -> list[int]:
    """Bernstein coefficients multiplied by the positive factor D*Q^d*d!."""
    require(left < right, "finite sign interval must have positive length")
    power = integer_affine(poly, left, right - left)
    degree = len(power) - 1
    weighted = [value * factorial(j) * factorial(degree - j)
                for j, value in enumerate(power)]
    return [sum(comb(i, j) * weighted[j] for j in range(i + 1))
            for i in range(degree + 1)]


def integer_halfline_bernstein(poly, start: Fraction) -> list[int]:
    """Integer scaling of the Bernstein coefficients for the half-line map."""
    power = integer_affine(poly, start, ONE)
    degree = len(power) - 1
    return [value * factorial(j) * factorial(degree - j)
            for j, value in enumerate(power)]


def integer_bisection(values: list[int]) -> tuple[list[int], list[int]]:
    """Both children's exact de Casteljau coefficients multiplied by 2^d."""
    degree = len(values) - 1
    row = list(values)
    left, right = [row[0] << degree], [row[-1] << degree]
    for level in range(1, degree + 1):
        row = [row[j] + row[j + 1] for j in range(len(row) - 1)]
        left.append(row[0] << (degree - level))
        right.append(row[-1] << (degree - level))
    right.reverse()
    return left, right


def prove_bernstein(coefficients: list[Fraction], *, strict: bool,
                    label: str, max_depth: int = MAX_DEPTH) -> None:
    """[proven] Return only if all leaves certify the requested sign."""
    require(bool(coefficients), "empty Bernstein coefficient list")
    require(type(max_depth) is int and max_depth >= 0,
            "invalid Bernstein depth limit")
    pending = [(integer_scale(coefficients), 0)]
    while pending:
        values, depth = pending.pop()
        if strict:
            require(values[0] > 0 and values[-1] > 0,
                    label + ": a required strictly positive endpoint fails")
        else:
            require(values[0] >= 0 and values[-1] >= 0,
                    label + ": a required nonnegative endpoint fails")
        # [proven] Nonnegative coefficients and positive endpoints give strict
        # positivity: every Bernstein basis function is positive in (0,1).
        if all(value >= 0 for value in values):
            continue
        require(depth < max_depth,
                label + ": Bernstein depth limit exhausted without proof")
        left, right = integer_bisection(values)
        pending.append((right, depth + 1))
        pending.append((left, depth + 1))


def prove_polynomial_obligations(cert: Certificate, p: list[Fraction],
                                 q: list[Fraction]) -> None:
    tau = [cert.tau0 * edge for edge in cert.E]
    prove_bernstein(integer_halfline_bernstein(q, ZERO), strict=True,
                    label="q on the nonnegative half-line")
    prove_bernstein(integer_halfline_bernstein(poly_scale(p, -ONE), tau[-1]),
                    strict=True, label="-p on the final half-line")
    for b, gamma in enumerate(cert.Gam):
        start = tau[b]
        envelope = [gamma * (1 - start / 2 + start * start / 8),
                    gamma * (Fraction(1, 2) - start / 4), gamma / 8]
        difference = poly_add(envelope, poly_scale(p, -ONE))
        prove_bernstein(integer_interval_bernstein(difference, start, tau[b + 1]),
                        strict=False, label="envelope bin " + str(b))


def exp_partial_sum(x: Fraction) -> Fraction:
    """[proven] Return sum(x**k/k!, k=0..512), evaluated exactly."""
    require(x >= 0, "exponential bound requires a nonnegative argument")
    result = ONE
    for k in range(EXP_TERMS, 0, -1):
        result = 1 + x * result / k
    return result


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def long_cap(n: int, S: Fraction) -> int:
    C = 1 - 2 / S
    denominator = 1 - n * C * C
    require(S > 2 and denominator > 0, "invalid long cap argument")
    return floor_fraction(n * (1 - C * C) / denominator)


def short_cap(n: int, T: Fraction, S: Fraction) -> int:
    require(T > 0 and T <= 2 and S > 2, "invalid short cap argument")
    if T + S <= 4:
        return 0
    U = (T + S - 4) ** 2 / (4 * T * S)
    require(U < 1, "invalid short cap projection denominator")
    C = 1 - 2 / S
    projected = (C + U) / (1 - U)
    result = long_cap(n, S) - 1
    denominator = 1 - (n - 1) * projected * projected
    if denominator > 0:
        result = min(result, floor_fraction(
            (n - 1) * (1 - projected * projected) / denominator))
    return result


def compute_rho(cert: Certificate, p: list[Fraction],
                q: list[Fraction]) -> Fraction:
    M, J = len(cert.Gam), cert.J
    g = [gamma / exp_partial_sum(cert.tau0 * cert.E[b] / 2)
         for b, gamma in enumerate(cert.Gam)] + [ZERO]
    for b in range(M - 1, J - 1, -1):
        g[b] = max(g[b], g[b + 1])
    empty = sum(((g[b] - g[b + 1]) * long_cap(cert.n, cert.E[b + 1])
                 for b in range(J, M)), ZERO)
    candidates = [empty]
    for j in range(J):
        candidates.append(g[j] + sum(
            ((g[b] - g[b + 1]) * short_cap(cert.n, cert.E[j + 1],
                                         cert.E[b + 1])
             for b in range(J, M)), ZERO))
    require(q[0] > 0, "q(0) must be positive")
    rho = (p[0] + 2 * (2 ** cert.n - 1) * max(candidates)) / q[0]
    require(rho > 0, "rho must be positive")
    return rho


def atan_bounds(x: Fraction) -> tuple[Fraction, Fraction]:
    """[proven] For 0<x<1, the 80-term sum is lower, 79-term sum upper."""
    require(0 < x < 1, "arctangent argument must lie in (0,1)")
    power, total, last = x, ZERO, ZERO
    for j in range(ATAN_TERMS):
        last = power / (2 * j + 1)
        if j % 2:
            last = -last
        total += last
        power *= x * x
    require(ATAN_TERMS % 2 == 0 and last < 0,
            "arctangent lower sum must have an even number of terms")
    return total, total - last


def pi_lower_bound() -> Fraction:
    """[proven] Machin's identity and alternating remainders give this bound."""
    x, y = Fraction(1, 5), Fraction(1, 239)
    twice = 2 * x / (1 - x * x)
    four_times = 2 * twice / (1 - twice * twice)
    require((four_times - y) / (1 + four_times * y) == 1,
            "Machin tangent identity failed")
    lower_x, _upper_x = atan_bounds(x)
    _lower_y, upper_y = atan_bounds(y)
    result = 16 * lower_x - 4 * upper_y
    require(result > 0, "the proven pi lower bound must be positive")
    return result


def check_claims(cert: Certificate, rho: Fraction) -> None:
    pi_lower = pi_lower_bound()
    m = cert.n // 2
    if cert.n % 2 == 0:
        density_bound = cert.tau0 ** m * rho / (8 ** m * factorial(m))
        require(cert.density >= density_bound,
                "claimed density is below the verified rational bound")
    else:
        odd_product = 1
        for integer in range(1, cert.n + 1, 2):
            odd_product *= integer
        R = 2 * cert.tau0 ** m * rho / (4 ** m * odd_product)
        require(cert.density ** 2 >= R ** 2 * cert.tau0 / (8 * pi_lower),
                "claimed density fails the verified squared bound")
    require(cert.hermite ** cert.n >=
            (cert.tau0 / (2 * pi_lower)) ** cert.n * rho ** 2,
            "claimed hermite fails the verified power bound")


def verify(cert: Certificate) -> None:
    p, q = laguerre_pair(cert)
    prove_polynomial_obligations(cert, p, q)
    rho = compute_rho(cert, p, q)
    check_claims(cert, rho)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Verify an exact certificate using Bernstein subdivision.")
    parser.add_argument("certificate", type=Path, help="certificate JSON path")
    args = parser.parse_args(argv)
    # [proven] No verdict is obtained by truncating an integer's decimal input.
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    try:
        cert = parse_certificate(args.certificate.read_text(encoding="utf-8"))
        verify(cert)
    except (Reject, OSError, UnicodeError, ValueError, OverflowError,
            RecursionError, MemoryError) as error:
        print("[computed] REJECT: " + str(error), file=sys.stderr)
        return 1
    print(f"PASS n={cert.n} K={cert.K} density={cert.density_text} "
          f"hermite={cert.hermite_text} [computed]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
