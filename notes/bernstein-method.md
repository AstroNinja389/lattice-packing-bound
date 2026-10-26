# Exact Bernstein verification

[proven] `verify/bernstein.py` uses only integers and `fractions.Fraction` for every arithmetic operation that determines acceptance. It independently parses the JSON, constructs Laguerre polynomials, proves the polynomial inequalities, computes the cap and exponential bounds, and compares the claimed decimals. A failed proof is a rejection; it does not assert that the underlying inequality is false.

## Input and output

[computed] Run `python verify/bernstein.py certificate.json`. A successful run prints `PASS n=... K=... density=... hermite=... [computed]` and exits 0. Rejection prints a diagnostic to stderr and exits 1. Command-line usage errors exit 2.

[proven] The parser rejects duplicate keys at every JSON nesting level, unexpected or missing fields, JSON noninteger numeric literals, NaN, Infinity, boolean dimensions, invalid array lengths, and rational strings outside `-?[0-9]+(/[1-9][0-9]*)?`. It checks every edge, dimension, and sign restriction in the verification contract before constructing a polynomial.

[computed] Claimed decimals use the grammar `[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?` and must have strictly positive exact rational value. Thus `+1`, `.5`, `1.`, and `1e-3` are supported spellings, while negative and zero values are rejected.

## Polynomial construction

[proven] Polynomials are stored as ascending rational power coefficients, with trailing zero coefficients removed and the zero polynomial represented by one zero coefficient. The recurrence starts with `L[-1]=0`, `L[0]=1` and uses `(k+1)L[k+1]=(2k+alpha+1-t)L[k]-(k+alpha)L[k-1]`, where `alpha=n/2-1`. Each term contributes `c[k] L[k]` to `p` and `(-1)^k c[k] L[k]` to `q`.

## Finite intervals and subdivision

[proven] The acceptance path clears one positive common denominator `D` from the power coefficients. If the affine substitution is `(A+B*x)/Q`, with positive integer `Q`, integer Horner evaluation gives the coefficients of `D*Q^d*f((A+B*x)/Q)`. Multiplication by `d!` clears the Bernstein denominators: `d!/binom(d,j)=j!*(d-j)!`. Thus the integer arrays are one positive multiple of the original Bernstein arrays. Midpoint subdivision uses integer pairwise sums; its level-`i` boundary entries are shifted by `d-i` bits, so both child arrays are exactly `2^d` times their rational de Casteljau counterparts. All sign decisions are unchanged. The original rational routines remain available for direct independent comparison, but are not imported from another verifier.

[proven] For a polynomial `f` on `[a,b]`, first form `f(a+(b-a)x)=sum A[j]x^j`. For its actual degree `d`, its degree-`d` Bernstein coefficients are `beta[i]=sum(j=0..i) A[j]*binom(i,j)/binom(d,j)`. This follows from `x^j=sum(i=j..d) binom(i,j)/binom(d,j) * binom(d,i)x^i(1-x)^(d-i)`.

[proven] The Bernstein basis functions are nonnegative on `[0,1]` and sum to one. Nonnegative coefficients therefore prove a nonnegative polynomial. If both endpoint coefficients are strictly positive, the same condition proves strict positivity throughout the closed interval, since each basis function is positive in the interior. Constant polynomials use the same criterion with a single coefficient.

[proven] If a coefficient is negative, midpoint de Casteljau subdivision computes the two child Bernstein coefficient arrays using only repeated pairwise averaging. Each child represents the same polynomial on its respective half interval. The algorithm checks endpoint signs at every node, accepts a node only from its Bernstein coefficient proof, and requires proofs for both children. Finding a negative endpoint rejects the requested nonnegative claim; finding a nonpositive endpoint rejects a strict claim.

[computed] Subdivision begins at depth 0 and permits nodes through depth 40. An uncertified depth-40 node causes rejection. No sampled sign is ever sufficient for acceptance. The traversal visits the left child first and has deterministic output. Nonnegative finite-bin polynomials with endpoint zeros, including the identically zero polynomial, can pass. Interior zeros that do not occur at a subdivision boundary may cause conservative rejection.

## Half-lines

[proven] To prove `f(t)>0` on `[a,infinity)`, let `d` be the actual degree and translate `f(a+y)=sum h[j]y^j`. Substitute `t=a+x/(1-x)` and multiply by `(1-x)^d`. The resulting polynomial is `F(x)=sum h[j]x^j(1-x)^(d-j)`, so its Bernstein coefficients are exactly `h[j]/binom(d,j)`. For `0<=x<1`, the multiplying factor is positive and the substitution covers exactly the original half-line.

[proven] Proving `F>0` on the entire closed interval `[0,1]` therefore suffices. Its value at 1 is the nonzero leading coefficient of `f`; any polynomial strictly positive on the entire half-line must have positive leading coefficient, including the degree-zero case. Thus using a closed transformed interval introduces no additional mathematical obstruction, although the finite subdivision budget can still reject a valid claim. This method is applied to `q` with `a=0` and to `-p` with `a=tau[M]`.

## Exponentials, caps, and final comparisons

[proven] At every nonnegative argument `x`, the Taylor series has nonnegative terms, hence `exp(-x)<=1/sum(k=0..512)x^k/k!`. The sum is evaluated exactly by backward Horner recurrence. The code uses the contract's suffix maxima only for long bins and then evaluates the stated `A`, `B`, `kappa`, and `rho` formulas with exact rational floor operations. Cap denominator checks remain explicit.

[proven] For `0<x<1`, integration of the finite geometric expansion of `1/(1+x^2)` shows that the 80-term alternating arctangent sum `L(x)` lies below `atan(x)` and the 79-term sum `U(x)` lies above it. Their terms have indices `0..79` and `0..78`, respectively. In particular, `pi_lower=16*L(1/5)-4*U(1/239)` is a strict rational lower bound once Machin's identity is established.

[proven] To establish that identity, let `theta=4 atan(1/5)-atan(1/239)`. Double-angle tangent formulas give `tan(2 atan(1/5))=5/12`, `tan(4 atan(1/5))=120/119`, and the tangent subtraction formula gives `tan(theta)=1`; the last equality is also checked using rational arithmetic at runtime. The inequalities `x-x^3/3<atan(x)<x` give `0<theta<4/5`. Also `pi/2=2 atan(1)>1`, since the integral of `1/(1+t^2)` from 0 to 1 exceeds `1/2`. Thus theta is in the principal tangent interval and equals `pi/4`. Consequently `pi=16 atan(1/5)-4 atan(1/239)`.

[proven] The verifier checks `pi_lower>0`, evaluates the contract's even-dimensional rational density bound directly, and uses the prescribed squared bound for odd dimensions. The positive Hermite claim is checked by its `n`th power against `(tau0/(2*pi_lower))^n*rho^2`. All decimal claims are exact fractions. Replacing pi by the proven smaller positive rational increases the tested upper bounds, so every successful comparison is sufficient for the stated contract.

## Scope of acceptance

[proven] Acceptance certifies the inequalities and formulas in `notes/verification-contract.md`. It does not independently derive the external geometric theorem that motivates those formulas. Exact Bernstein tests are sufficient certificates; rejection at depth exhaustion is deliberately inconclusive about the mathematical sign.

## Reproducible edge-case checks

[computed] The implementation checks completed on 2026-10-23 covered 243 exact arithmetic, parser, sign, and coefficient-identity assertions, followed by seven CLI assertions. The CLI accepted both the even and odd examples below with exit 0, rejected a Hermite claim reduced to 1 with exit 1, and returned exit 2 without a positional argument. The coefficient checks compared direct polynomial evaluation with both Bernstein transformations and both halves of subdivision in degrees 0 through 8 using rational arguments.

[computed] This smaller repeatable check covers constants, endpoint zeros, interior zeros, depth exhaustion, and complete even and odd certificates. Run it in the repository directory with PowerShell:

```powershell
@'
from fractions import Fraction as F
import json
import runpy

v = runpy.run_path('verify/bernstein.py')
passed = 0

def check_sign(power, strict, expected, halfline=False, depth=40):
    global passed
    power = list(map(F, power))
    coefficients = (v['halfline_bernstein'](power, F(0)) if halfline
                    else v['interval_bernstein'](power, F(0), F(1)))
    try:
        v['prove_bernstein'](coefficients, strict=strict,
                             label='edge case', max_depth=depth)
        accepted = True
    except v['Reject']:
        accepted = False
    assert accepted == expected
    passed += 1

check_sign([1], True, True)
check_sign([1], True, True, halfline=True)
check_sign([0], False, True)
check_sign([0], True, False)
check_sign([-1], False, False)
check_sign([0, 1, -1], False, True)
check_sign([0, 1, -1], True, False)
check_sign([F(1, 4), -1, 1], False, True)
check_sign([F(1, 4), -1, 1], True, False)
check_sign([F(1, 9), F(-2, 3), 1], False, False)
check_sign([F(26, 100), -1, 1], True, True)
check_sign([F(26, 100), -1, 1], True, False, depth=0)
check_sign([1, -2, 2], True, True, halfline=True)
check_sign([1, -2, 1], True, False, halfline=True)
check_sign([1, -1], True, False, halfline=True)

example = {'n': 2, 'K': 1, 'tau0': '2', 'E': ['1', '2', '3'],
           'c': ['2', '1'], 'Gam': ['1', '0'],
           'claims': {'density': '2', 'hermite': '2'}}
v['verify'](v['parse_certificate'](json.dumps(example)))
passed += 1
example.update(n=3, Gam=['3/2', '0'],
               claims={'density': '3', 'hermite': '4'})
v['verify'](v['parse_certificate'](json.dumps(example)))
passed += 1
print('[computed] PASS:', passed, 'repeatable checks')
'@ | python -
```

[computed] Expected and observed output: `[computed] PASS: 17 repeatable checks`.
