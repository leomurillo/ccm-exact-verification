"""probe_ccm_exact_matrix.py — certificate for P-infinity-E SS21 (the exact CCM matrix).

Settles SS20's open experimental questions by building the EXACT matrix of the semilocal
Weil form QW_lambda in the CCM basis V_n = kappa(U_n) (arXiv 2511.22755v1, Prop toadd +
Lemma polarize0 + Lemma wsharp), in mpmath high precision:

  QW(V_n, V_m) = Psi#(F_nm),  F_nm(x) = q(U_n, U_m)(log x),
  Psi# = W02# - WR# - sum_p Wp#,
  q(U_n,U_m)(y) = [sin(w_m y) - sin(w_n y)]/(pi(n-m))  (n != m),  w_k = 2 pi k / L,
  q(U_n,U_n)(y) = 2(1 - y/L) cos(w_n y),    L = 2 log lambda,
  W02#(q) = int_0^L q(y) 2 cosh(y/2) dy,
  WR#(q)  = (1/2)(log 4pi + gamma) q(0) + int_0^L [e^{y/2} q(y) - q(0)]/(2 sinh y) dy
            + q(0) * (1/2) log tanh(L/2)          [the F(1)-tail beyond L, closed form],
  Wp#(q)  = sum_{p^j <= e^L} (log p) p^{-j/2} q(j log p).

The form is LINEAR in q, so all off-diagonal entries come from 2N+1 per-frequency
integrals S_k = A_k - B_k - P_k:  QW[n,m] = (S_m - S_n)/(pi(n-m)).

CHECK A  validation: QW[0,0] and QW[1,1] against the zero side (200 zeros + tail estimate).
CHECK B  THE VERDICT (cliff or ladder): even/odd low spectra of QW_lambda^N at
         lambda = 1.3 (pure archimedean), 1.5 (prime 2), 2.0 (primes 2,3,4):
         is the even ground state simple? is the gap O(1) (cliff) or ~lambda^2*eps (ladder)?
CHECK C  eps_lambda vs the Fuchs scale 1-chi_4 ~ (2^14/3) sqrt2 pi^5 e^{-4 pi l^2 + 9 log l}
         (their Figure fpro1 identity, reproduced at two lambdas).
CHECK D  r(lambda) at lambda = 1.5: Rayleigh excess of the FD-prolate guess k_lambda in the
         EXACT matrix: the weld's crux ratio r/eps measured for the first time.

Deterministic, local, ~2-4 min.
"""
import sys
import numpy as np
from scipy.linalg import eig_banded
import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 40

def primes_upto(x):
    out = []
    n = 2
    while n <= x:
        if all(n % p for p in out):
            out.append(n)
        n += 1
    return out

def build_QW(lam, N):
    """Exact QW matrix on V_n, |n| <= N, plus even/odd sector matrices."""
    L = 2 * mp.log(lam)
    w = lambda k: 2 * mp.pi * k / L
    pj = []                                          # (j log p, log p * p^{-j/2})
    for p in primes_upto(float(mp.exp(L)) + 1e-12):
        j = 1
        while j * mp.log(p) <= L + mp.mpf("1e-30"):
            pj.append((j * mp.log(p), mp.log(p) * mp.power(p, -mp.mpf(j) / 2)))
            j += 1
    pts = mp.linspace(0, L, 16)
    # per-frequency functionals for the sin-parts (off-diagonal building blocks)
    S = {}
    for k in range(-N, N + 1):
        wk = w(k)
        A = mp.quad(lambda y: mp.sin(wk * y) * 2 * mp.cosh(y / 2), pts)
        B = mp.quad(lambda y: mp.exp(y / 2) * mp.sin(wk * y) / (2 * mp.sinh(y)) if y != 0 else wk / 2, pts)
        P = mp.fsum(cw * mp.sin(wk * yp) for yp, cw in pj)
        S[k] = A - B - P
    # diagonal entries
    cgamma = mp.euler
    tailL = mp.log(mp.tanh(L / 2)) / 2               # -int_L^inf dy/(2 sinh y) = +tail term sign
    D = {}
    for k in range(0, N + 1):
        wk = w(k)
        qd = lambda y: 2 * (1 - y / L) * mp.cos(wk * y)
        A = mp.quad(lambda y: qd(y) * 2 * mp.cosh(y / 2), pts)
        Bint = mp.quad(lambda y: (mp.exp(y / 2) * qd(y) - 2) / (2 * mp.sinh(y)) if y != 0
                       else (mp.mpf(1) / 2 - 1 / L), pts)
        WR = (mp.log(4 * mp.pi) + cgamma) / 2 * 2 + Bint + 2 * tailL
        P = mp.fsum(cw * (2 * (1 - yp / L) * mp.cos(wk * yp)) for yp, cw in pj)
        D[k] = A - WR - P
    # assemble full matrix, indices n = -N..N
    dim = 2 * N + 1
    QW = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            n, m = i - N, j - N
            QW[i, j] = D[abs(n)] if n == m else (S[m] - S[n]) / (mp.pi * (n - m))
    # even sector: e_0 = V_0, e_k = (V_k + V_{-k})/sqrt2 ; odd: o_k = (V_k - V_{-k})/sqrt2
    Ev = mp.matrix(N + 1, N + 1)
    Od = mp.matrix(N, N)
    for k in range(N + 1):
        for l in range(N + 1):
            a, b = QW[N + k, N + l], QW[N + k, N - l]
            if k == 0 and l == 0:
                Ev[k, l] = a
            elif k == 0 or l == 0:
                Ev[k, l] = mp.sqrt(2) * a
            else:
                Ev[k, l] = a + b
    for k in range(1, N + 1):
        for l in range(1, N + 1):
            Od[k - 1, l - 1] = QW[N + k, N + l] - QW[N + k, N - l]
    return QW, Ev, Od, L

def low_spec(Mt, count=6):
    E = sorted([mp.re(e) for e in mp.eigsy(Mt, eigvals_only=True)])
    return E[:count]

# ---------------- CHECK A: validation against the zero side ----------------
print("== CHECK A: exact matrix vs zero side (lambda = 1.5, 200 zeros) ==")
lam0 = mp.mpf("1.5")
QW0, Ev0, Od0, L0 = build_QW(lam0, 8)
GAM = [mp.im(mp.zetazero(k)) for k in range(1, 201)]

def zero_side_diag(k, L):
    wk = 2 * mp.pi * k / L
    tot = mp.mpf(0)
    for g in GAM:
        qhat = mp.quad(lambda y: 2 * (1 - y / L) * mp.cos(wk * y) * mp.cos(g * y),
                       mp.linspace(0, L, 12))
        tot += 2 * 2 * qhat              # both +-gamma, and qhat = 2 int q cos
    return tot

okA = True
for k in (0, 1):
    zs = zero_side_diag(k, L0)
    ex = QW0[8 + k, 8 + k]
    gmax = GAM[-1]
    tail = 8 / (mp.pi * L0) * (mp.log(gmax / (2 * mp.pi)) + 1) / gmax
    e = abs(ex - zs)
    okA &= e < 3 * tail
    print(f"   QW[V_{k},V_{k}]: exact = {mp.nstr(ex, 8)}   zero side = {mp.nstr(zs, 8)}   "
          f"|diff| = {mp.nstr(e, 2)} (tail allowance {mp.nstr(3*tail, 2)})")
print("   " + ("PASS — exact matrix certified against the zeros" if okA else "FAIL"))

# ---------------- CHECK B: the verdict ----------------
print("\n== CHECK B: low spectra — cliff or ladder ==")
results = {}
for lam_s, N in (("1.3", 10), ("1.5", 10), ("2.0", 12)):
    lam = mp.mpf(lam_s)
    QWm, Ev, Od, L = build_QW(lam, N)
    ev = low_spec(Ev, 6)
    od = low_spec(Od, 3)
    results[lam_s] = (ev, od)
    print(f"   lambda = {lam_s} (L = {mp.nstr(L,4)}, primes <= {float(lam)**2:.2f}):")
    print("      even: " + "  ".join(mp.nstr(e, 6) for e in ev))
    print("      odd : " + "  ".join(mp.nstr(e, 6) for e in od))
    eps, mu1 = ev[0], ev[1]
    simple_even = eps < od[0] and mu1 / max(eps, mp.mpf("1e-38")) > 10
    print(f"      eps = {mp.nstr(eps, 4)}  (even, {'SIMPLE' if simple_even else 'check'});  "
          f"gap mu1 - eps = {mp.nstr(mu1 - eps, 4)};  mu1/eps = {mp.nstr(mu1/max(eps, mp.mpf('1e-38')), 4)}")
okB = True

# ---------------- CHECK C: eps vs the Fuchs scale ----------------
print("\n== CHECK C: eps_lambda vs Fuchs 1-chi_4 scale ==")
for lam_s in ("1.5", "2.0"):
    lam = mp.mpf(lam_s)
    fuchs = mp.mpf(2)**14 / 3 * mp.sqrt(2) * mp.pi**5 * mp.exp(-4 * mp.pi * lam**2 + 9 * mp.log(lam))
    eps = results[lam_s][0][0]
    print(f"   lambda = {lam_s}: eps = {mp.nstr(eps, 4)}   Fuchs(1-chi_4) = {mp.nstr(fuchs, 4)}   "
          f"ratio = {mp.nstr(eps / fuchs, 4)}")

# ---------------- CHECK D: the weld crux r/eps at lambda = 1.5 ----------------
print("\n== CHECK D: Rayleigh excess of the prolate guess in the EXACT matrix ==")
lam_f = 1.5
Lf = 2 * np.log(lam_f)
M = 4001
x = np.linspace(-lam_f, lam_f, M); dx = x[1] - x[0]
xm = (x[:-1] + x[1:]) / 2
a = lam_f**2 - xm**2
diag = np.zeros(M); diag[:-1] += a / dx**2; diag[1:] += a / dx**2
diag += (2 * np.pi * lam_f * x)**2
bands = np.vstack([diag, np.append(-a / dx**2, 0.0)])
wv, vv = eig_banded(bands, lower=True, select='i', select_range=(0, 5))
h0l = vv[:, 0] / np.sqrt(dx); h4l = vv[:, 4] / np.sqrt(dx)
if h0l[M // 2] < 0: h0l = -h0l
if h4l[M // 2] < 0: h4l = -h4l
I0, I4 = np.trapezoid(h0l, x), np.trapezoid(h4l, x)
hl = h4l - (I4 / I0) * h0l
xg = np.linspace(0.0, Lf, 2001)
u = np.exp(xg - Lf / 2)
kl = np.zeros_like(u)
for n in range(1, int(lam_f / u.min()) + 2):
    msk = n * u <= lam_f
    kl[msk] += np.interp(n * u[msk], x[x >= 0], hl[x >= 0])
kl *= np.sqrt(u)
N8 = 8
c = []
for n in range(-N8, N8 + 1):
    ph = np.exp(-2j * np.pi * n * xg / Lf) / np.sqrt(Lf)
    c.append(np.trapezoid(kl * ph, xg))
c = mp.matrix([mp.mpc(z) for z in c])
num = mp.fsum(mp.conj(c[i]) * QW0[i, j] * c[j] for i in range(17) for j in range(17))
den = mp.fsum(mp.conj(c[i]) * c[i] for i in range(17))
R = mp.re(num / den)
eps15 = results["1.5"][0][0]
r = R - eps15
print(f"   Rayleigh(k_guess) = {mp.nstr(R, 6)}   eps = {mp.nstr(eps15, 6)}   "
      f"r = R - eps = {mp.nstr(r, 4)}   r/eps = {mp.nstr(r / eps15, 4)}")
print("   (FD-guess quality bounds r from above; r_true <= measured r)")

allok = okA
print("\nRESULT: " + ("GREEN — exact CCM matrix built and certified; see CHECK B for the "
                      "cliff/ladder verdict and CHECK D for the weld crux" if allok else "RED"))
