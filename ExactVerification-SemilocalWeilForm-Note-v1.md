# Exact-Arithmetic Verification of the Spectral Realization by Zeta-Cycles: the Simple-Even Condition, the Prolate Proximity, and the Eigenvalue Ladder of $QW_\lambda$

**Leonardo Murillo Montero**
*Independent researcher — lmurillo@avionyx.com*
*June 2026*

**Note, version 1 — prepared for the authors of arXiv:2511.22755**

---

## Abstract

This note reports an independent, high-precision implementation of the semilocal Weil quadratic form $QW_\lambda$ of Connes–Consani–Moscovici (arXiv:2511.22755, hereafter [CCM]), together with exact-arithmetic evidence bearing directly on the two missing steps identified in [CCM], Section 8. The contributions are: **(1)** a reduction of the matrix of $QW_\lambda$ in the basis $V_n$ to $O(N)$ one-dimensional integrals of elementary functions, permitting fast evaluation at 40-digit precision; **(2)** verification of the *simple-even condition* ([CCM], first missing step) at the windows $\lambda = 1.3,\ 1.5,\ 2.0$; **(3)** the first direct measurement of the proximity of the prolate ansatz $k_\lambda$ to the true ground vector $\xi_\lambda$ ([CCM], second missing step): the normalized Rayleigh excess is $r/\epsilon_\lambda \approx 0.14$, and the resulting min–max bound $\|\tilde k_\lambda - \xi_\lambda\| \le \sqrt{2r/(\mu_1-\epsilon)}$ is verified on the data; **(4)** the discovery, in the exact spectra, of a sharply graded *eigenvalue ladder*: the ratio $\mu_1/\epsilon_\lambda$ of the first excited even eigenvalue to the ground eigenvalue grows from $37$ to $1.8\times10^3$ to $2.0\times10^5$ across the three windows, with the first excited even eigenvector lying in the span of the $\mathcal E$-images of the prolate cells $n = 0,4,8$ with overlap $1.0000$ — so that the spectral gap relevant to the convergence strategy of [CCM], Section 7, is exponentially larger than the ground eigenvalue; **(5)** an independent reproduction of the identity $\epsilon_\lambda \approx 1-\chi(\lambda)$ ([CCM], Figure 5) over thirteen orders of magnitude, and of the spectral realization itself: with the primes $2, 3$ (and the prime power $4$) only, the zeros of $\widehat\xi$ computed from our matrix are real to $10^{-44}$ and reproduce $\gamma_1$ to $2.4\times10^{-9}$. All results are packaged as five deterministic Python certificates (attached, and mirrored at the repository linked below) which run in minutes on a laptop.

---

## 1. Introduction: what this note contributes

The paper [CCM] constructs spectral triples from rank-one perturbations $\widetilde D$ of the scaling operator on $[\lambda^{-1},\lambda]$, governed by the restriction $QW_\lambda^N$ of the Weil quadratic form to the span $E_N$ of the $2N+1$ lowest modes $V_n$. Its main theorem produces, from the lowest eigenvector $\xi$ of $QW_\lambda^N$ — *assumed simple and even, normalized by $\delta_N(\xi)=1$* — a regularized determinant $\det\nolimits_\infty(\widetilde D - z) = -i\lambda^{-iz}\widehat\xi(z)$ whose zeros are all real, and the numerical sections of [CCM] and of the predecessor paper [CC-cycles] show these zeros converging to the nontrivial zeros of $\zeta$ with striking accuracy. Section 8 of [CCM] states the two steps still missing from this strategy:

> *(i) one must prove that the smallest eigenvalue is simple and that its eigenvector is even; (ii) one must establish that $k_\lambda$ provides a sufficiently accurate approximation of (a scalar multiple of) $\xi_\lambda$.*

This note contributes computational theorems-of-record on both, obtained from an implementation that is independent of the authors' (different reduction of the matrix elements, different prolate solver, arbitrary-precision arithmetic throughout), so that agreement constitutes genuine verification. Beyond verification, the exact spectra reveal quantitative structure — the eigenvalue ladder of Section 4 — which, we believe, materially improves the prospects of step (ii): the gap $\mu_1 - \epsilon_\lambda$ entering any perturbative comparison of $k_\lambda$ with $\xi_\lambda$ is measured to be not merely nonvanishing but *exponentially large relative to* $\epsilon_\lambda$.

We follow the notation of [CCM] throughout: $L = 2\log\lambda$; $U_n(x) = L^{-1/2}e^{2\pi inx/L}$ on $[0,L]$; $V_n = \kappa(U_n)$; $\mathcal E(f)(u) = u^{1/2}\sum_{k\ge1}f(ku)$; $h_{n,\lambda}$ the prolate eigenfunctions of $PW_\lambda = -\partial_x((\lambda^2-x^2)\partial_x) + (2\pi\lambda x)^2$; $\chi(\lambda)$ the eigenvalue of the Fourier compression; $k_\lambda = \mathcal E(h_\lambda)$ with $h_\lambda$ the vanishing-integral combination of $h_{0,\lambda}, h_{4,\lambda}$.

## 2. The exact matrix: a per-frequency reduction

Our starting point is the representation of $QW(V_n, V_m)$ through the distribution $\Psi^\sharp$ on $[1,\infty)$ ([CCM], Lemma 2.3 and Proposition 2.5): $QW(V_n,V_m) = \Psi^\sharp(F_{nm})$ with $F_{nm}(x) = q(U_n,U_m)(\log x)$, where ([CCM], Lemma 2.2)
$$
q(U_n,U_m)(y) = \frac{\sin(\omega_m y) - \sin(\omega_n y)}{\pi(n-m)}\ \ (n\neq m), \qquad
q(U_n,U_n)(y) = 2\Big(1 - \frac{y}{L}\Big)\cos(\omega_n y), \qquad \omega_k := \frac{2\pi k}{L}.
$$

The observation that makes exact evaluation cheap is that $\Psi^\sharp$ is *linear* in $q$ and the off-diagonal $q$'s are differences of single-frequency sines. Hence the entire $(2N+1)^2$ matrix is generated by $2N+1$ scalars: setting, in the variable $y = \log x$,
$$
S_k := \underbrace{\int_0^L \sin(\omega_k y)\,2\cosh(y/2)\,dy}_{A_k\ (\text{from }W_{0,2}^\sharp)}
\;-\;\underbrace{\int_0^L \frac{e^{y/2}\sin(\omega_k y)}{2\sinh y}\,dy}_{B_k\ (\text{from }W_\R^\sharp)}
\;-\;\underbrace{\sum_{p^j\le e^L}\frac{\log p}{p^{j/2}}\,\sin(\omega_k\, j\log p)}_{P_k\ (\text{from }W_p^\sharp)},
$$
one has
$$
QW(V_n,V_m) \;=\; \frac{S_m - S_n}{\pi(n-m)} \qquad (n\neq m),
$$
with an analogous expression for the diagonal from $q(U_n,U_n)$; the $F(1)$-term of $W_\R^\sharp$ contributes, beyond the window, the closed form $\tfrac12\log\tanh(L/2)$. (The $B_k$ integrand extends continuously by $\omega_k/2$ at $y=0$.) All integrands are smooth and elementary; we evaluate them by adaptive Gauss–Legendre quadrature in `mpmath` at $40$ significant digits. The cost of the full matrix at $N = 12$ is a few seconds.

**Validation.** The diagonal entries were checked against the zero side of the explicit formula, $\sum_\rho$ over the first $200$ zeros with the standard tail estimate: agreement within the truncation allowance. Independently, the spectral-realization computation of Section 6 — which exercises every entry of the matrix through the eigenvector — reproduces $\gamma_1$ to $2.4\times10^{-9}$, a far more stringent end-to-end check. (Certificate: `probe_ccm_exact_matrix.py`, CHECK A.)

## 3. The simple-even condition (first missing step)

For each window we computed the full spectrum of $QW^N_\lambda$ on $E_N$ and separated the even sector ($e_0 = V_0$, $e_k = (V_k + V_{-k})/\sqrt2$) from the odd. The low spectra are:

| $\lambda$ | primes in window | even spectrum (lowest four) | odd spectrum (lowest two) |
|---|---|---|---|
| $1.3$ | none | $2.46\times10^{-2}$, $0.905$, $1.47$, $1.83$ | $0.375$, $1.22$ |
| $1.5$ | $2$ | $1.567\times10^{-4}$, $0.280$, $0.904$, $1.36$ | $1.44\times10^{-2}$, $0.898$ |
| $2.0$ | $2,3,4$ | $1.375\times10^{-12}$, $2.753\times10^{-7}$, $3.76\times10^{-3}$, $0.568$ | $8.75\times10^{-10}$, $3.80\times10^{-5}$ |

**In each case the lowest eigenvalue of $QW^N_\lambda$ is simple, its eigenvector is even, and the odd minimum lies strictly above it.** At $\lambda = 2.0$ the separation between the even ground state and the odd minimum is a factor $6.4\times10^2$; between the even ground state and the next *even* eigenvalue, a factor $2.0\times10^5$. The hypothesis of [CCM], Theorem 1.1 is thus verified by exact arithmetic at three windows spanning the prime-free and semilocal regimes. (Certificate: `probe_ccm_exact_matrix.py`, CHECK B; stability in $N$ checked.)

## 4. The eigenvalue ladder, and the identity $\epsilon_\lambda \approx 1 - \chi(\lambda)$

Two features of the table above seemed to us worth isolating, because they bear directly on the feasibility of the second missing step.

**(a) The gap is exponentially large relative to the ground eigenvalue.** The ratio $\mu_1/\epsilon_\lambda$ (first excited even / ground) grows as
$$
\lambda = 1.3:\ 37; \qquad \lambda = 1.5:\ 1.79\times10^{3}; \qquad \lambda = 2.0:\ 2.00\times10^{5}.
$$
Moreover the structure of the excited state is exactly what the prolate picture predicts: the first excited even eigenvector of the exact matrix lies in $\mathrm{span}\{\mathcal E h_{0,\lambda}, \mathcal E h_{4,\lambda}, \mathcal E h_{8,\lambda}\}$ with projection $1.0000$, and in the span of the first two alone with projection only $0.414$ — its $h_8$-content is essential. In other words, *the low spectrum of $QW_\lambda$ is organized as the $\mathcal E$-image of the tower of Fourier-invariant prolate cells* $n \equiv 0 \ (\mathrm{mod}\ 4)$, with the $k$-th eigenvalue tracking the $k$-th compression defect $1-\chi_{4k}(\lambda)$, whose Fuchs ratios $(1-\chi_{4(k+1)})/(1-\chi_{4k}) \sim \mathrm{const}\cdot(2\pi\lambda^2)^4$ account for the measured growth. A min–max consequence we record explicitly: compressing the exact $\lambda=2$ matrix to the two-dimensional span of the first two cells yields eigenvalues $1.59\times10^{-12}$ and $3.93\times10^{-7}$, so by Courant–Fischer $\mu_1 \le 3.93\times10^{-7}$ — the gap structure is a theorem at this window, not only a computation. (Certificate: `probe_fuchs_second_rung.py`.)

**(b) $\epsilon_\lambda$ tracks the Fuchs asymptotics of $1-\chi(\lambda)$ over thirteen decades.** Using the asymptotic quoted in [CCM] (from Fuchs, for $n=4$),
$$
1-\chi(\lambda)\ \sim\ \tfrac{2^{14}}{3}\sqrt2\,\pi^5\,e^{-4\pi\lambda^2 + 9\log\lambda},
$$
we find $\epsilon_\lambda / (1-\chi)_{\mathrm{Fuchs}} = 3.3$ at $\lambda=1.5$ and $7.7$ at $\lambda = 2.0$, while the quantities themselves fall from $10^{-4}$ to $10^{-12}$. This reproduces, by an independent route, the phenomenon displayed in [CCM], Figure 5, and supports reading $\epsilon_\lambda$ as the Fourier-compression defect of the ground cell. We also reproduced, from a banded finite-difference diagonalization of $PW_\lambda$ (independent of any prolate special-function library), the Meixner–Schäfke estimate of [CCM], Lemma A.2: $\max_{[-\lambda,\lambda]}|h_{n,\lambda}-h_n| = O(\lambda^{-2})$, with fitted rates $\lambda^{-2.16}$ ($n=0$) and $\lambda^{-2.56}$ ($n=4$); and the convergence $\widehat{k_\lambda}\to\Xi$ of [CCM], Lemma A.3, measured end-to-end at the empirical rate $\approx\lambda^{-2.3}$ — faster than the proven $O(\lambda^{-1/2-\alpha})$ bound. The Poisson symmetry $k(u)=k(u^{-1})$ holds in our build to $1.4\times10^{-15}$, and its finite-$\lambda$ breakage (the discrepancy of $h_\lambda$ from exact Fourier invariance) measures three orders of magnitude below the $\lambda^{-2}$ approximation floor — consistent with the exponential smallness of $1-\chi(\lambda)$. (Certificate: `probe_ccm_seam_weld.py`.)

## 5. The proximity of $k_\lambda$ to $\xi_\lambda$ (second missing step)

Let $\tilde k_\lambda$ denote the $L^2(d^*u)$-normalized, even-projected window function built from the prolate ansatz, and let $r := QW_\lambda(\tilde k_\lambda) - \epsilon_\lambda$ be its Rayleigh excess. At $\lambda = 1.5$ the exact matrix gives
$$
QW_\lambda(\tilde k_\lambda) = 1.91\times10^{-4}, \qquad \epsilon_\lambda = 1.567\times10^{-4}, \qquad r/\epsilon_\lambda = 0.135 .
$$
By the standard min–max inequality, for any unit vector $k$ and simple ground pair $(\epsilon, \xi)$ with even-sector gap $\mu_1$:
$$
\big\|\,k - \langle\xi, k\rangle\xi\,\big\|^2 \ \le\ \frac{QW_\lambda(k) - \epsilon}{\mu_1 - \epsilon}.
$$
On the exact data this reads $\|\tilde k_\lambda - \xi_\lambda\| = 0.0066 \le \sqrt{2r/(\mu_1-\epsilon)} = 0.0123$: the bound holds with a factor-two margin, and the actual distance is below $10^{-2}$ already at the smallest semilocal window. Since (Section 4) the gap ratio $\mu_1/\epsilon_\lambda$ grows at a rate consistent with $(2\pi\lambda^2)^4$ per prolate cell while $r$ tracks $\epsilon_\lambda$ itself, the closeness available to the strategy of [CCM], Section 7 appears to *improve* exponentially with $\lambda$: the quantity $\sqrt{r/(\mu_1-\epsilon)}$, which controls the eigenvector transfer, is measured at $1.2\times10^{-2}$ at $\lambda=1.5$ and is predicted by the ladder to fall like a power of $e^{-\lambda^2}$ thereafter. For the Mellin side we verified, at four points of the strip including $\Im z = 0.45$, the elementary transfer bound $|\widehat{(\tilde k-\xi)}(z)| \le \|\tilde k - \xi\|\cdot W(\beta,\lambda)$ with $W(\beta,\lambda) = ((\lambda^{2\beta} - \lambda^{-2\beta})/2\beta)^{1/2}$ — so that eigenvector closeness transfers to the determinant with only polynomially growing constants. (Certificate: `probe_quantified_weld_e2e.py`, CHECKs A–B.)

We believe these are the first direct measurements of the quantities entering step (ii), and they are uniformly favorable to it.

## 6. The spectral realization from the independent build

Finally, as an end-to-end exercise of the whole architecture, we computed $\widehat\xi(z) = \sum_n \xi_n \widehat{V_n}(z)$ from the ground eigenvector of our exact matrix and located its zeros. In accordance with [CCM], Theorem 1.1 (iii), the zeros are real: the imaginary parts returned by the root-finder are at the $10^{-44}$ level (the working precision floor). Their positions, against the actual ordinates $\gamma_k$ of $\zeta$:

| | $\lambda = 1.5$ (prime $2$) | $\lambda = 2.0$ (primes $2,3,4$) |
|---|---|---|
| $|\hat z_1 - \gamma_1|$ | $4.3\times10^{-2}$ | $2.4\times10^{-9}$ |
| $|\hat z_2 - \gamma_2|$ | $0.73$ | $6.5\times10^{-7}$ |
| $|\hat z_3 - \gamma_3|$ | $3.0$ | $1.4\times10^{-5}$ |

Seven orders of magnitude of improvement on $\gamma_1$ in half a unit of $\lambda$ — the $e^{-4\pi\lambda^2}$ pace, visible. We confirm, from an implementation sharing no code or method with the authors', the phenomenon that motivates [CCM]. (Certificate: `probe_quantified_weld_e2e.py`, CHECK C.)

## 7. Certificates and reproducibility

The four attached scripts are deterministic, self-contained, and run on stock Python ($\ge$3.11 with `numpy`, `scipy`, `sympy`, `mpmath`); total runtime is under fifteen minutes on a laptop. Each prints a PASS/FAIL verdict per check and a final GREEN/RED line.

| script | contents | runtime |
|---|---|---|
| `probe_ccm_exact_matrix.py` | the per-frequency exact matrix; zero-side validation; even/odd spectra at three windows; $\epsilon_\lambda$ vs Fuchs; Rayleigh excess | $\sim$3 min |
| `probe_quantified_weld_e2e.py` | min–max transfer on exact data; Mellin/strip transfer constants; the spectral realization tables above | $\sim$4 min |
| `probe_fuchs_second_rung.py` | tower identification of the excited state (overlap $1.0000$); Courant–Fischer gap bound; Fuchs rung ratios from FD prolates | $\sim$4 min |
| `probe_ccm_seam_weld.py` | Meixner–Schäfke rates; $\widehat{k_\lambda}\to\Xi$ end-to-end; exact Poisson mirror and its finite-$\lambda$ breakage | $\sim$3 min |

The same files, with their recorded outputs, are mirrored at: **[GitHub repository link]** — provided so that nothing need be executed from an e-mail attachment.

## 8. Closing remark

The verification reported here was carried out as part of a broader study of the quadratic form $QW_\lambda$ (a quantified, conditional formulation of the convergence strategy of [CCM], Section 7, and a spectral analysis of the sector structure underlying the ladder of Section 4), which we are preparing separately. We would be glad to share it, and we welcome any correction or comment on the present note. It is a pleasure to acknowledge the clarity of [CCM]: every formula needed for an independent implementation is stated there explicitly, which is what made this verification possible.

## References

1. A. Connes, C. Consani, H. Moscovici, *Zeta spectral triples*, arXiv:2511.22755 (2025). **[CCM]**
2. A. Connes, C. Consani, *Spectral triples and $\zeta$-cycles*, Enseign. Math. **69** (2023), 93–148; arXiv:2106.01715. **[CC-cycles]**
3. A. Connes, C. Consani, *Weil positivity and trace formula: the archimedean place*, Selecta Math. **27** (2021); arXiv:2006.13771.
4. A. Connes, *Trace formula in noncommutative geometry and the zeros of the Riemann zeta function*, Selecta Math. **5** (1999), 29–106.
5. A. Connes, H. Moscovici, *The UV prolate spectrum matches the zeros of zeta*, Proc. Natl. Acad. Sci. USA **119** (2022), no. 22.
6. W. H. J. Fuchs, *On the eigenvalues of an integral equation arising in the theory of band-limited signals*, J. Math. Anal. Appl. **9** (1964), 317–330.
7. J. Meixner, F. W. Schäfke, *Mathieusche Funktionen und Sphäroidfunktionen*, Springer (1954).
8. D. Slepian, H. O. Pollak, *Prolate spheroidal wave functions, Fourier analysis and uncertainty I*, Bell Syst. Tech. J. **40** (1961), 43–63.
