# Exact certificate contract

Specification for the independent exact verifiers.

Each certificate is a JSON object with exactly n, K, tau0, E, c, Gam, claims. n and K are JSON integers (not booleans); n >= 2, K >= 1. Other numeric data are strings matching `-?[0-9]+(/[1-9][0-9]*)?`. Claims have exactly density and hermite, positive decimal strings, optionally scientific notation. Reject duplicate JSON keys, NaN, Infinity, invalid types, extra fields, zero denominators, and negative or zero claims. Exit 0 on acceptance, 1 on rejection, 2 for command-line usage errors.

E contains squared edges: E[0]=1, strictly increasing, one E[J]=2 with 0<J<M, E[M]<4, n*(1-2/E[M])^2<1. There are M+1 edges, M nonnegative Gam values, K+1 coefficients c. tau0>0. All arithmetic used for a verdict is exact rational or integer arithmetic.

Set alpha=n/2-1, p=sum c[k] L_k^alpha(t), q=sum (-1)^k c[k] L_k^alpha(t), tau[b]=tau0 E[b]. Verify q strictly positive on [0,infinity), p strictly negative on [tau[M],infinity), and Gam[b]*(1+(t-tau[b])/2+(t-tau[b])^2/8)-p(t)>=0 on [tau[b],tau[b+1]]. Stronger strict interval positivity is acceptable but may reject valid boundary cases. The Laguerre recurrence is (k+1)L[k+1]=(2k+alpha+1-t)L[k]-(k+alpha)L[k-1].

For deterministic identical rational bounds in every implementation, use N=512 and exp(-x)<=1/sum_{k=0}^N x^k/k!. Define raw g[b]=Gam[b]/sum_{k=0}^N (tau[b]/2)^k/k!. Replace g[b], b>=J, by max(raw g[b],...,raw g[M-1]); g[M]=0. Short-bin g values are not changed.

For a long squared edge S, C=1-2/S and B(S)=floor(n*(1-C^2)/(1-n*C^2)). For a short squared edge T<=2: A(T,S)=0 if T+S<=4. Otherwise U=(T+S-4)^2/(4*T*S), Cprime=(C+U)/(1-U), A=B(S)-1; if (n-1)*Cprime^2<1 replace A by min(A, floor((n-1)*(1-Cprime^2)/(1-(n-1)*Cprime^2))).

kappa_empty=sum_{b=J}^{M-1}(g[b]-g[b+1])*B(E[b+1]); kappa_j=g[j]+sum_{b=J}^{M-1}(g[b]-g[b+1])*A(E[j+1],E[b+1]). kappa=max(kappa_empty,kappa_0,...,kappa_{J-1}); rho=(p(0)+2*(2^n-1)*kappa)/q(0). Require rho>0.

Prove a rational pi lower bound within each verifier using Machin's identity pi=16 atan(1/5)-4 atan(1/239) and alternating Taylor bounds: use the 80-term lower sum for atan(1/5) and the 79-term upper sum for atan(1/239). Never trust an unproved decimal constant. Even n=2m: density bound D=tau0^m*rho/(8^m*m!). Odd n=2m+1: D=R*sqrt(tau0/(8*pi)) where R=2*tau0^m*rho/(4^m*(1*3*...*(2m+1))). Check positive claimed density d by d^2>=R^2*tau0/(8*pi_lower). Check positive claimed Hermite H by H^n >= (tau0/(2*pi_lower))^n*rho^2. Even density is checked without squaring. These are sufficient inequalities, not assertions that rational overestimation is optimal.

Every implementation independently parses, builds polynomials, proves signs, computes caps, exponential bounds and claimed decimals. No cross-imports. Each supports positional certificate path and produces deterministic PASS fields n,K,density,hermite; optional diagnostics go to stderr. Default no acceptance from sampled signs. Conservative rejection on depth exhaustion is required.

PARI input: Python driver emits only validated integer or fraction tokens, never certificate strings without validation. GP itself also validates input types. Run gp -q -f. Set parisizemax before read(verify.gp), in driver-generated input. Require explicit machine-readable completion marker as well as successful process status. No success inferred from a partial run.
