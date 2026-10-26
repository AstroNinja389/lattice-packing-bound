# Lattice sphere-packing certificates from classes modulo two

This repository gives exact polynomial certificates for lattice packing density and Hermite bounds obtained from Theorem 1 of the paper. Each certificate is checked independently by three methods in every dimension, and by a fourth method in even dimensions. The numerical search supplies candidates; it is not part of the proof.

Licence selection is pending; no licence is granted by this checkout.

## What a certificate proves

Let a lattice in dimension `n` have minimum length `r`, and let `rho` denote its reciprocal covolume. For a real Schwartz function `h` with nonnegative Fourier transform and positive `hat h(0)`, Theorem 1 proves

```
rho <= (h(0) + 2*(2^n-1)*kappa(G)) / hat h(0).
```

Here `h(x) <= G(|x|/r)` from radius `r` to `e_M*r`, and `h <= 0` beyond. The step weight has partition points `1=e_0<...<e_M<2`, including `e_J=sqrt(2)`, and nonnegative bin values, nonincreasing for `b>=J`. Its last partition point satisfies `n*(1-2/e_M^2)^2<1`. The geometric caps and budget are

```
c = 1 - 2/s^2
B(s) = floor(n*(1-c^2)/(1-n*c^2))
A(a,s) = 0                                  if a^2+s^2 <= 4
c1 = (a^2+s^2-4)/(2*a*s)                    otherwise
cprime = (c+c1^2)/(1-c1^2)
A(a,s) = B(s)-1
```

In the last case, when `(n-1)*cprime^2<1`, replace `A` by its minimum with `floor((n-1)*(1-cprime^2)/(1-(n-1)*cprime^2))`. With `g_M=0`,

```
kappa_empty = sum(b=J..M-1, (g_b-g_(b+1))*B(e_(b+1)))
kappa_j = g_j + sum(b=J..M-1, (g_b-g_(b+1))*A(e_(j+1),e_(b+1)))
kappa = max(kappa_empty, kappa_0, ..., kappa_(J-1)).
```

The paper proves these formulas, including their half-open endpoint conventions.

## Certificate format and exact conditions

One JSON file contains exactly `n`, `K`, `tau0`, `E`, `c`, `Gam`, and `claims`. Dimensions and degrees are JSON integers. All other input numbers are integer or fraction strings, such as `"3"` or `"7/5"`; claims are strictly positive decimal strings. `E` contains **squared** partition points. `c` has `K+1` entries and `Gam` has one entry per bin. Duplicate keys and invalid numeric forms are rejected.

Set `alpha=n/2-1`, `tau_b=tau0*E_b`, and

```
p(t) = sum(k=0..K, c_k*L_k^alpha(t))
q(t) = sum(k=0..K, (-1)^k*c_k*L_k^alpha(t)).
```

Every verifier checks the partition-point structure, `tau0>0`, `Gam_b>=0`, and

1. `q(t)>0` for all `t>=0`.
2. `Gam_b*(1+u+u^2/2)-p(t)>=0` on every closed bin, where `u=(t-tau_b)/2`.
3. `p(t)<0` for all `t>=tau_M`.

Some verification methods require strict positivity in condition 2; this stronger condition suffices. An exhausted subdivision budget rejects instead of accepting an unresolved sign.

For `x=tau_b/2`, replace `exp(-x)` by `1/sum(j=0..512,x^j/j!)`. Multiply by `Gam_b` and take backward running maxima for long bins. Each verifier independently computes these rational weights, the caps, the budget and `rho_+`. It checks the decimal claims against

```
density <= tau0^(n/2)*rho_+ / (8^(n/2)*Gamma(n/2+1))
Hermite <= tau0/(2*pi)*rho_+^(2/n).
```

For odd dimensions, the density comparison is squared after checking that the claimed decimal is positive. The Hermite comparison is raised to the `n`th power. A rational lower bound for pi is derived from Machin's identity and alternating arctangent sums. No floating-point arithmetic determines a verdict.

More explicitly, for `n=2*m+1` the density expression is
`2*tau0^m*rho_+/(4^m*(2*m+1)!!) * sqrt(tau0/(8*pi))`.
The rational lower bound for pi uses 80 terms of the alternating series for
`arctan(1/5)` and 79 terms for `arctan(1/239)` in
`16*arctan(1/5)-4*arctan(1/239)`.

## Verification

Requirements: Python 3.10 or later for the standard-library verifiers; PARI/GP 2.19.0 for root counting; python-flint 0.9.0 for the even-dimensional Descartes verifier. The tested Python version is 3.14.6. Replace `python` below with the executable for the intended installation.

```
python verify/run.py --verifier taylor
python verify/run.py --verifier bernstein
python verify/run.py --verifier pari --gp gp
python verify/run.py --verifier descartes
```

The driver checks every applicable manifest entry, validates file hashes and metadata, and compares output with the committed reference log, ignoring timing fields. Descartes applies only to even dimensions. Empty runs cannot pass. Exit status is `0` for successful checks, `1` for rejection or incomplete verification, and `2` for command-line usage errors.

Add `--jobs 2` to check two certificates concurrently. The default is one process;
values from 1 through 16 are accepted. Output remains in manifest order.

To check a modified copy without comparing its output with a reference log:

```
python verify/run.py --verifier bernstein --file modified.json --no-reference
```

The individual verifiers also accept a certificate path directly. PARI is always invoked with `-q -f`; its driver validates every integer and fraction before constructing input. The GP program independently rejects non-rational arguments. These checks prevent certificate strings from becoming executable code.

Reference logs are never overwritten by ordinary runs. To explicitly record a complete manifest run:

```
python verify/run.py --verifier taylor --full-run --record-reference
```

The four verifiers share no implementation code. Their mathematical algorithms are exact real-root counting, integer Taylor bounds with interval bisection, Bernstein coefficients with subdivision depth at most 40, and Descartes sign variations. The certificate format and geometric theorem are their common specification.

## Upper-bound results

<!-- RESULTS_BEGIN -->
The manifest contains 83 certificates in dimensions 9 through 91. The paper reports improvements in 67 dimensions: 11â€“23, 25â€“40, 43â€“57, 62â€“75, 83â€“91.

The comparisons are rigorous published three-point bounds through dimension 16 and numerical linear-programming values from dimension 17 onward. Dimensions 9 and 10 are included for comparison with the three-point bounds, although stronger lattice bounds are available: the [dimension-9 preprint](https://arxiv.org/abs/2508.20719) reports the exact optimum in dimension 9, and Mordell's inequality applied to that result gives a stronger bound in dimension 10. Dimension 24 is a control near the known optimum. Every positive margin is included; the four smallest are listed below. These exact fractions express percentage gains relative to the stated comparison decimals.

| Dimension | Gain (%) | Exact gain (%) |
| ---: | ---: | --- |
| 62 | 0.0117098 | `5272553708663909626274429300/450269625537086639096262744293` |
| 91 | 0.0166780 | `884518112685501987118054241100/53034941781126855019871180542411` |
| 75 | 0.0176948 | `77745549783113201672042272700/4393698715497831132016720422727` |
| 83 | 0.0859357 | `93687099345589990246783970100/1090200460993455899902467839701` |

In dimension 16, the certificate gives density at most `0.0187461206` using Theorem 1, with polynomial degree 81.

Full density and Hermite bounds, comparison status, truncation checks and exact gains are in [the result data](data/results.json) and [the paper](paper/main.pdf).
<!-- RESULTS_END -->

## Repository contents

The repository contains 83 upper-bound certificates in dimensions 9 through 91, the manifest, four independent exact verification methods, and their reference logs. It also includes the mathematical verification notes, the result data for those 83 certificates, and the self-contained paper source and PDF. Each certificate records its density and Hermite bounds in its claims field.

The certificates were obtained by numerical linear programming followed by exact verification. The numerical search is not part of the proof. Exact verification requires no network access or numerical optimizer.
