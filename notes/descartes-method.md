# Exact Descartes verifier

[proven] `verify/descartes.py` constructs the two Laguerre polynomials independently from the recurrence in the certificate contract, using `python-flint` rational polynomial arithmetic. It accepts only even dimensions. It validates every field and rejects duplicate JSON keys and non-finite JSON constants before using any certificate data. Positive decimal claims use the grammar `[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?` and are converted to fractions using integer powers of ten.

## Compact interval proof

[proven] For a degree-d polynomial f and a<b, form

    Q(z) = (1+z)^d f((a+b z)/(1+z)).

[proven] The substitution sends z in [0,infinity) onto [a,b). The coefficients of Q are rational; clearing denominators by a positive common multiple and removing a positive common divisor does not change any sign. The constant and highest coefficient are positive multiples of f(a) and f(b). Both are required to be strictly positive.

[proven] Descartes' rule of signs states that the number of positive real zeros, counted with multiplicity, does not exceed the number of coefficient sign variations after zeros are omitted. In particular, zero sign variations together with positive endpoint coefficients proves strict positivity throughout the interval. This special case also follows immediately because every nonzero coefficient is positive. No acceptance is based on point samples.

[proven] When sign variations remain, let m=(a+b)/2. The two exact integer coefficient transformations are

    L(w) = (2+w)^d Q(w/(2+w)),
    R(w) = Q(1+2w).

[proven] L and R are 2^d times the corresponding Descartes transforms for [a,m] and [m,b]. Their coefficients are

    L[k] = 2^(d-k) sum(i=0..k) binom(d-i,k-i) Q[i],
    R[k] = 2^k     sum(i=k..d) binom(i,k) Q[i].

[proven] The shared midpoint has the sign of sum(Q[i]); a nonpositive sum rejects the strict test. Otherwise each child is processed with the same rule. Acceptance requires every leaf to have zero variations. Every transformation uses integer arithmetic. Exceeding depth 256 or 200000 visited nodes for an individual sign proof rejects the certificate conservatively.

## Half-line proof

[proven] For f on [a,infinity), use the polynomial compactification

    h(u) = (1-u)^d f(a+u/(1-u)),  0 <= u < 1.

[proven] The expression extends polynomially to u=1, where h(1) is the leading coefficient of f. The initial Descartes transform for h on [0,1] is exactly f(a+z). Thus the same integer subdivision proof handles the entire half-line, including an exact endpoint test at infinity. The constant-polynomial case is included. No finite cutoff or approximate root bound is used. The requirements are q>0 on its full half-line, -p>0 on its tail, and strict positivity of every bin majorant. The latter is stronger than the contract's nonnegative requirement and can conservatively reject a valid boundary certificate.

## Rational bound calculations

[proven] For x>=0, the exponential series has positive terms, so exp(x)>=sum(k=0..512) x^k/k! and therefore exp(-x)<=1/sum(k=0..512) x^k/k!. The implementation constructs that sum exactly. It computes each cap using the contract's rational expressions and integer floor, takes the suffix maximum only for long-bin weights, and then constructs every kappa candidate and rho independently.

[proven] For 0<x<1, the alternating arctangent series gives a lower bound from 80 terms (indices 0 through 79) and an upper bound from 79 terms (indices 0 through 78). Write these as l(x) and u(x). Machin's identity gives the certified rational lower bound

    pi_lower = 16*l(1/5) - 4*u(1/239).

[proven] For completeness, setting a=atan(1/5), b=atan(1/239) gives tan(2a)=5/12, tan(4a)=120/119, and tan(4a-b)=1. The angle is positive since a>5/26 and b<1/239, and it is below 4/5<pi/2. The inequalities for arctangent follow by integrating 1/(1+t^2); pi/2=2*atan(1)>1 follows by the same integral. Thus 4a-b=pi/4. The implementation checks the rational tangent identity and checks that the resulting enclosure lies between 3 and 4.

[proven] For even n=2m the verifier computes D=tau0^m*rho/(8^m*m!) and checks the claimed positive density is at least D. It checks the claimed positive Hermite number H using H^n >= (tau0/(2*pi_lower))^n*rho^2. These operations and comparisons are exact. Odd dimensions are explicitly rejected with exit status 1.

## Interface and checks

[computed] Running `descartes.py certificate.json` prints deterministic `PASS n=... K=... density=... hermite=... [proven]` fields only after every proof and comparison has succeeded. Rejection returns 1; command-line misuse returns 2. The optional `--self-test` mode checks accepted strict polynomials, excluded zeros, a negative half-line pocket, budget exhaustion, malformed data, rejected claims, and a complete accepted certificate.
