# Integer Taylor verification

[proven] `verify/taylor.py` uses Python standard-library integer and rational arithmetic for every acceptance decision. Its polynomial construction, sign proof, cap arithmetic, exponential estimate, and decimal comparison are self-contained.

[computed] Run `python verify/taylor.py certificate.json`; acceptance prints `PASS [proven] n=... K=... density=... hermite=...`. Exit codes are 0 for acceptance, 1 for rejection, and 2 for command-line usage errors. `--self-test` exercises exact arithmetic, both parities, malformed inputs, sign failures, and conservative depth exhaustion. `--max-depth` defaults to 40 and can only decrease that limit. Each bounded interval also has a 200,000-node limit.

## Input and polynomial construction

[proven] JSON object keys are checked for duplication at every nesting level; JSON decimal numbers and nonfinite constants are rejected. Only `n` and `K` may be JSON integers, and booleans are rejected for both. The top-level object and claims object must have exactly their specified fields. Rational strings match `-?[0-9]+(/[1-9][0-9]*)?`. Decimal claims match `[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?` and must have strictly positive exact values. Arrays, lengths, edge ordering, internal edge 2, nonnegative Gamma, and the strict last-edge cap condition are checked before the proof.

[proven] The Laguerre recurrence in the contract constructs rational coefficient lists for p and q. Every sign-test polynomial is multiplied by a positive common denominator, then divided by the positive greatest common divisor of its integer coefficients. These transformations preserve its sign. No floating-point evaluation participates in a verdict.

## Bounded intervals

[proven] Write P(a+(b-a)x) with x in [0,1]. If a=A/D and b-a=B/D, the homogeneous Horner recurrence computes D^d P((A+Bx)/D) using integers only. Its coefficients h_i therefore represent a positive multiple of the original polynomial on the entire interval.

[proven] The first lower bound is h_0 + sum(min(h_i,0), i>=1), since each x^i lies in [0,1]. A second bound is centered: compute S(z)=2^d H((1+z)/2)=sum(s_i z^i), where z lies in [-1,1]. Then S(z) is at least s_0 minus the sum of absolute values of odd coefficients plus the sum of negative even coefficients. This follows term by term from z^i in [-1,1] for odd i and in [0,1] for positive even i. A strictly positive bound certifies a strict sign condition; a nonnegative bound certifies a weak condition.

[proven] If neither bound succeeds, bisect. The left child has coefficients 2^(d-i)h_i; the right child has coefficients s_i. Positive common factors can be removed from either child. Covering both closed children proves the parent interval. Endpoint or midpoint values can refute the requested sign but never suffice to accept an interval. Exhausting either budget rejects an unresolved sign. Some valid nonnegative polynomials with interior zeros can therefore be rejected; acceptance always has a covering proof.

## Infinite intervals

[proven] A positive constant term and nonnegative remaining coefficients of P(R+u) prove P(t)>0 for every t>=R. The verifier first tries the specified starting point, then integer R values increasing by doubling. A finite prefix [start,R] is proved by the bounded-interval algorithm.

[proven] The search has a finite exact stopping point. For P(t)=a_d t^d+...+a_0 with a_d>0, put A=max(|a_i|, i<d) and R_0=1+ceil(A/a_d). For t>=R_0, t>=1 and, when A>0, the geometric-series estimate gives P(t)>=a_d t^d-A(t^d-1)/(t-1)>0. When A=0 the conclusion follows from the positive monomial. Each nonzero derivative has positive leading coefficient and lower-to-leading absolute coefficient ratios at most A/a_d, because each lower coefficient is multiplied by a falling-factorial factor no larger than the leading coefficient's factor. The same estimate proves every derivative positive at R_0. Consequently the Taylor coefficients of P(R+u) are positive for every R>=R_0; an explicit coefficient check is still required. Degree-zero polynomials are handled directly. The chosen limit is also at least ceil(start).

[proven] Applying this procedure to q proves q>0 on [0,infinity); applying it to -p proves p<0 on the required tail. The slab checks use Gamma times the quadratic in the contract minus p and require nonnegativity on each closed interval.

## Exponentials, caps, and final comparisons

[proven] For x>=0, the degree-512 exponential partial sum is positive and at most exp(x), so its reciprocal is at least exp(-x). Each sum is formed by the exact recurrence term_k=term_(k-1)*x/k. Gamma is nonnegative, preserving the required direction. The long-bin suffix maximum and terminal zero are applied exactly as specified; short-bin values are unchanged. Cap floors, kappa, and rho use exact rationals. A nonpositive rho rejects the certificate.

[proven] For 0<x<1, the alternating arctangent sum with 80 terms (indices 0 through 79) is a lower bound, and the sum with 79 terms (indices 0 through 78) is an upper bound. The verifier uses pi_lower=16*lower_80(1/5)-4*upper_79(1/239). This direction accounts for the negative coefficient of the second arctangent.

[proven] Machin's identity is justified independently of decimal constants. Let a=atan(1/5) and b=atan(1/239). Twice applying tan(2x)=2tan(x)/(1-tan(x)^2) gives tan(4a)=120/119. Thus tan(4a-b)=(120/119-1/239)/(1+(120/119)(1/239))=1, an identity also checked with rationals in the program. The integral bounds x/(1+x^2)<atan(x)<x give 0<4a-b<4/5. Moreover atan(1)>1/2, so pi/2>1. Hence 4a-b lies in (0,pi/2), where tangent is injective, and equals pi/4. This proves pi=16a-4b and the asserted rational lower bound.

[proven] The even-dimensional density comparison uses its exact rational expression without squaring. In odd dimensions, positivity of both sides allows the specified squared inequality using pi_lower. The positive Hermite claim is checked after raising it to n, against (tau0/(2*pi_lower))^n*rho^2. The pi lower bound enlarges these right-hand sides, so each accepted comparison is sufficient for the corresponding bound in the contract.

## Recorded verification

[computed] The 35 built-in self-test checks pass. An additional deterministic set of 139 constructive checks passed: 60 products of squared rational linear factors plus 1/100 were certified on a bounded interval and on an infinite tail, and 19 quadratic polynomials with strictly negative interior minima were rejected. Seven subprocess checks confirmed full certificate-file acceptance and rejection and the actual child-process exit codes 0, 1, and 2. These tests check the implementation; the covering and tail arguments above justify acceptance.
