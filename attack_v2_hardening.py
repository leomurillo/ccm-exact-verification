"""attack_v2_hardening.py — adversarial audit of the four headline numbers of the
ExactVerification note (v1 -> v2 hardening). Referee-style attacks:

ATTACK 1  N-convergence + quadrature-degree + dps stability of the spectra:
          eps, mu_1 (even), odd minimum at lambda = 1.3/1.5/2.0 under
          N -> N+2 -> N+4, quadrature nodes 16 -> 32, dps 40 -> 60.
          Freezes (or kills) the table in note Section 3 and the ladder ratios 37/1787/2e5
          and the Fuchs ratios 3.3/7.7.
ATTACK 2  the overlap "1.0000": print 1 - P^2 to 10 digits; vary FD grid (4001 -> 8001),
          E-grid (3001 -> 6001), N (12 -> 14); control: random 3-dim spans in the
          13-dim even sector (expected ~ sqrt(3/13) = 0.48). Also the 0.4142 ~ sqrt2 - 1
          curiosity on the first-cell projection.
ATTACK 3  zeros: root-finder basin (4 distinct starts per zero), reality at dps 40 vs 60
          (if |Im z| is a precision floor it must FALL with dps), N = 12 vs 14 stability
          of z_1..z_3 at lambda = 2.
ATTACK 4  r/eps at lambda = 1.5: FD contamination budget by grid refinement
          (M = 2001/4001/8001, E-grid 3001/6001): the drift bounds the FD error;
          consistent triple (QW(k~), eps, r/eps) for the note's Section 5.

Deterministic, local. ~15-25 min.
"""
import sys
import numpy as np
from scipy.linalg import eig_banded
import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")

def primes_upto(x):
    out = []; n = 2
    while n <= x:
        if all(n % p for p in out):
            out.append(n)
        n += 1
    return out

def build_QW(lam, N, npts=16):
    """Exact QW matrix on V_n, |n| <= N, + even/odd sectors. Quadrature subdivision npts."""
    L = 2 * mp.log(lam)
    w = lambda k: 2 * mp.pi * k / L
    pj = []
    for p in primes_upto(float(mp.exp(L)) + 1e-12):
        j = 1
        while j * mp.log(p) <= L + mp.mpf("1e-30"):
            pj.append((j * mp.log(p), mp.log(p) * mp.power(p, -mp.mpf(j) / 2)))
            j += 1
    pts = mp.linspace(0, L, npts)
    S = {}
    for k in range(-N, N + 1):
        wk = w(k)
        A = mp.quad(lambda y: mp.sin(wk * y) * 2 * mp.cosh(y / 2), pts)
        B = mp.quad(lambda y: mp.exp(y / 2) * mp.sin(wk * y) / (2 * mp.sinh(y)) if y != 0 else wk / 2, pts)
        P = mp.fsum(cw * mp.sin(wk * yp) for yp, cw in pj)
        S[k] = A - B - P
    tailL = mp.log(mp.tanh(L / 2)) / 2
    D = {}
    for k in range(0, N + 1):
        wk = w(k)
        qd = lambda y: 2 * (1 - y / L) * mp.cos(wk * y)
        A = mp.quad(lambda y: qd(y) * 2 * mp.cosh(y / 2), pts)
        Bint = mp.quad(lambda y: (mp.exp(y / 2) * qd(y) - 2) / (2 * mp.sinh(y)) if y != 0
                       else (mp.mpf(1) / 2 - 1 / L), pts)
        WR = (mp.log(4 * mp.pi) + mp.euler) + Bint + 2 * tailL
        P = mp.fsum(cw * (2 * (1 - yp / L) * mp.cos(wk * yp)) for yp, cw in pj)
        D[k] = A - WR - P
    dim = 2 * N + 1
    QW = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            n, m = i - N, j - N
            QW[i, j] = D[abs(n)] if n == m else (S[m] - S[n]) / (mp.pi * (n - m))
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

def low_even_odd(Ev, Od, ne=4, no=2):
    ee = sorted([mp.re(e) for e in mp.eigsy(Ev, eigvals_only=True)])[:ne]
    oo = sorted([mp.re(e) for e in mp.eigsy(Od, eigvals_only=True)])[:no]
    return ee, oo

def fd_prolates(lam_f, M, idx=(0, 4, 8)):
    x = np.linspace(-lam_f, lam_f, M); dx = x[1] - x[0]
    xm = (x[:-1] + x[1:]) / 2
    a = lam_f**2 - xm**2
    diag = np.zeros(M); diag[:-1] += a / dx**2; diag[1:] += a / dx**2
    diag += (2 * np.pi * lam_f * x)**2
    bands = np.vstack([diag, np.append(-a / dx**2, 0.0)])
    _, vv = eig_banded(bands, lower=True, select='i', select_range=(0, max(idx)))
    out = {}
    for n in idx:
        h = vv[:, n] / np.sqrt(dx)
        if h[M // 2] < 0: h = -h
        out[n] = h
    return x, out

def Emap_coeffs(h, x, lam_f, L, N, Mg):
    """E(h) on the window, Fourier coefficients on U_n, n = -N..N."""
    Lf = float(L)
    xg = np.linspace(0.0, Lf, Mg)
    u = np.exp(xg - Lf / 2)
    kl = np.zeros_like(u)
    for n in range(1, int(lam_f / u.min()) + 2):
        msk = n * u <= lam_f
        kl[msk] += np.interp(n * u[msk], x[x >= 0], h[x >= 0])
    kl *= np.sqrt(u)
    c = []
    for n in range(-N, N + 1):
        ph = np.exp(-2j * np.pi * n * xg / Lf) / np.sqrt(Lf)
        c.append(complex(np.trapezoid(kl * ph, xg)))
    return c

def to_even(c, N):
    """Full U_n coefficients -> even-sector coordinates (e_0, e_1..e_N)."""
    v = [mp.mpc(c[N])]
    for k in range(1, N + 1):
        v.append((mp.mpc(c[N + k]) + mp.mpc(c[N - k])) / mp.sqrt(2))
    return mp.matrix(v)

def gs_proj(vecs, target):
    """orthonormalize vecs (list of mp column vectors), return |proj of target|."""
    basis = []
    for v in vecs:
        w = v.copy()
        for b in basis:
            ov = mp.fsum(mp.conj(b[i]) * w[i] for i in range(len(w)))
            w = mp.matrix([w[i] - ov * b[i] for i in range(len(w))])
        nrm = mp.sqrt(mp.re(mp.fsum(mp.conj(w[i]) * w[i] for i in range(len(w)))))
        if nrm > mp.mpf("1e-25"):
            basis.append(w / nrm)
    p2 = mp.mpf(0)
    for b in basis:
        ov = mp.fsum(mp.conj(b[i]) * target[i] for i in range(len(target)))
        p2 += abs(ov) ** 2
    tn = mp.re(mp.fsum(mp.conj(target[i]) * target[i] for i in range(len(target))))
    return mp.sqrt(p2 / tn)

def Vhat(n, z, L):
    wn = 2 * mp.pi * n / L
    if abs(wn - z) < mp.mpf("1e-25"):
        return mp.exp(1j * z * L / 2) * mp.sqrt(L)
    return mp.exp(1j * z * L / 2) / mp.sqrt(L) * (mp.exp(1j * (wn - z) * L) - 1) / (1j * (wn - z))

# =====================================================================
print("==================== ATTACK 1: spectra stability ====================")
mp.mp.dps = 40
configs = [
    ("1.3", 10, 16, 40, "baseline"),
    ("1.3", 12, 16, 40, "N+2"),
    ("1.5", 10, 16, 40, "baseline"),
    ("1.5", 12, 16, 40, "N+2"),
    ("1.5", 14, 16, 40, "N+4"),
    ("1.5", 10, 32, 40, "quad x2"),
    ("1.5", 10, 16, 60, "dps 60"),
    ("2.0", 12, 16, 40, "baseline"),
    ("2.0", 14, 16, 40, "N+2"),
    ("2.0", 16, 16, 40, "N+4"),
    ("2.0", 12, 32, 40, "quad x2"),
    ("2.0", 12, 16, 60, "dps 60"),
]
store = {}
for lam_s, N, npts, dps, tag in configs:
    mp.mp.dps = dps
    QWm, Ev, Od, L = build_QW(mp.mpf(lam_s), N, npts)
    ee, oo = low_even_odd(Ev, Od)
    store[(lam_s, tag)] = (ee, oo, QWm, Ev, L, N)
    mu_ratio = ee[1] / ee[0]
    print(f"  lam={lam_s} N={N} pts={npts} dps={dps} [{tag}]")
    print(f"    even: {'  '.join(mp.nstr(e, 9) for e in ee)}")
    print(f"    odd : {'  '.join(mp.nstr(o, 9) for o in oo)}")
    print(f"    eps = {mp.nstr(ee[0], 9)}   mu1/eps = {mp.nstr(mu_ratio, 6)}")
mp.mp.dps = 40
print("\n  --- drift summary (relative |delta eps| and |delta mu1/eps| vs baseline) ---")
for lam_s in ("1.3", "1.5", "2.0"):
    base = store[(lam_s, "baseline")]
    e0, m0 = base[0][0], base[0][1] / base[0][0]
    for tag in ("N+2", "N+4", "quad x2", "dps 60"):
        if (lam_s, tag) in store:
            ee = store[(lam_s, tag)][0]
            de = abs(ee[0] - e0) / e0
            dm = abs(ee[1] / ee[0] - m0) / m0
            print(f"    lam={lam_s} [{tag}]: |d eps|/eps = {mp.nstr(de, 3)}   |d(mu1/eps)|/ = {mp.nstr(dm, 3)}")
print("\n  --- Fuchs ratios under variation ---")
for lam_s in ("1.5", "2.0"):
    lam = mp.mpf(lam_s)
    fuchs = mp.mpf(2)**14 / 3 * mp.sqrt(2) * mp.pi**5 * mp.exp(-4 * mp.pi * lam**2 + 9 * mp.log(lam))
    for tag in ("baseline", "N+2", "N+4", "quad x2", "dps 60"):
        if (lam_s, tag) in store:
            eps = store[(lam_s, tag)][0][0]
            print(f"    lam={lam_s} [{tag}]: eps/Fuchs = {mp.nstr(eps / fuchs, 5)}")

# =====================================================================
print("\n==================== ATTACK 2: the overlap 1.0000 ====================")
mp.mp.dps = 40
for (Ntag, Mfd, Mg) in ((("2.0", "baseline"), 4001, 3001),
                        (("2.0", "baseline"), 8001, 6001),
                        (("2.0", "N+2"), 4001, 3001)):
    ee, oo, QWm, Ev, L, N = store[Ntag]
    E, Q = mp.eigsy(Ev)
    order = sorted(range(N + 1), key=lambda i: mp.re(E[i]))
    v1 = mp.matrix([Q[k, order[1]] for k in range(N + 1)])     # first excited even
    lam_f = 2.0
    x, hs = fd_prolates(lam_f, Mfd, idx=(0, 4, 8))
    cells = []
    for n in (0, 4, 8):
        c = Emap_coeffs(hs[n], x, lam_f, L, N, Mg)
        cells.append(to_even(c, N))
    P3 = gs_proj(cells, v1)
    P2 = gs_proj(cells[:2], v1)
    print(f"  lam=2.0 [{Ntag[1]}] FDgrid={Mfd} Egrid={Mg}:")
    print(f"    |proj v1 onto span(Eh0,Eh4,Eh8)| = {mp.nstr(P3, 10)}    1 - P = {mp.nstr(1 - P3, 4)}")
    print(f"    |proj v1 onto span(Eh0,Eh4)|     = {mp.nstr(P2, 10)}")
    print(f"    (sqrt2 - 1 = {mp.nstr(mp.sqrt(2) - 1, 10)} for comparison)")
# random-span control at baseline
ee, oo, QWm, Ev, L, N = store[("2.0", "baseline")]
E, Q = mp.eigsy(Ev)
order = sorted(range(N + 1), key=lambda i: mp.re(E[i]))
v1 = mp.matrix([Q[k, order[1]] for k in range(N + 1)])
rng = np.random.default_rng(20260611)
projs = []
for t in range(20):
    vecs = [mp.matrix([mp.mpf(float(z)) for z in rng.standard_normal(N + 1)]) for _ in range(3)]
    projs.append(gs_proj(vecs, v1))
print(f"  control: 20 RANDOM 3-dim spans in the {N+1}-dim even sector:")
print(f"    max |proj| = {mp.nstr(max(projs), 4)}   mean = {mp.nstr(mp.fsum(projs)/20, 4)}   "
      f"(sqrt(3/{N+1}) = {mp.nstr(mp.sqrt(mp.mpf(3)/(N+1)), 4)})")

# =====================================================================
print("\n==================== ATTACK 3: zeros — basin + dps + N ====================")
GAM = [mp.im(mp.zetazero(k)) for k in (1, 2, 3)]
for (tag, dps) in (("baseline", 40), ("N+2", 40), ("baseline", 60)):
    mp.mp.dps = dps
    if dps == 60:
        QWm, Ev, Od, L = build_QW(mp.mpf("2.0"), 12, 16)
        N = 12
    else:
        ee, oo, QWm, Ev, L, N = store[("2.0", tag)]
    E, Q = mp.eigsy(QWm)
    o = sorted(range(2 * N + 1), key=lambda i: mp.re(E[i]))
    xv = mp.matrix([Q[k, o[0]] for k in range(2 * N + 1)])
    xihat = lambda z: mp.fsum(xv[N + n] * Vhat(n, z, L) for n in range(-N, N + 1))
    print(f"  lam=2.0 [{tag}] dps={dps}:")
    for k, g in enumerate(GAM):
        roots = []
        for s in (mp.mpc(g, 0.01), mp.mpc(g + 0.3, 0.0), mp.mpc(g - 0.3, -0.05), mp.mpc(g, 0.5)):
            try:
                zr = mp.findroot(xihat, s)
                roots.append(zr)
            except Exception:
                roots.append(None)
        good = [z for z in roots if z is not None]
        re_spread = max(abs(mp.re(a) - mp.re(b)) for a in good for b in good) if len(good) > 1 else mp.mpf(0)
        im_max = max(abs(mp.im(z)) for z in good)
        print(f"    gamma_{k+1}: {len(good)}/4 starts converged; Re spread = {mp.nstr(re_spread, 3)}; "
              f"max|Im| = {mp.nstr(im_max, 3)}; |Re z - gamma| = {mp.nstr(abs(mp.re(good[0]) - g), 4)}")
mp.mp.dps = 40

# =====================================================================
print("\n==================== ATTACK 4: r/eps FD contamination budget ====================")
ee, oo, QWm, Ev, L15, N15 = store[("1.5", "baseline")]
E, Q = mp.eigsy(QWm)
order = sorted(range(2 * N15 + 1), key=lambda i: mp.re(E[i]))
eps = mp.re(E[order[0]])
mu1 = None
for idx in order[1:]:
    v = Q[:, idx]
    if all(abs(v[N15 + k] - v[N15 - k]) < mp.mpf("1e-10") for k in range(1, N15 + 1)):
        mu1 = mp.re(E[idx]); break
for (Mfd, Mg) in ((2001, 3001), (4001, 3001), (8001, 3001), (4001, 6001), (8001, 6001)):
    lam_f = 1.5
    x, hs = fd_prolates(lam_f, Mfd, idx=(0, 4))
    I0 = np.trapezoid(hs[0], x); I4 = np.trapezoid(hs[4], x)
    hl = hs[4] - (I4 / I0) * hs[0]
    c = Emap_coeffs(hl, x, lam_f, L15, N15, Mg)
    ce = to_even(c, N15)
    nrm = mp.sqrt(mp.re(mp.fsum(mp.conj(ce[i]) * ce[i] for i in range(N15 + 1))))
    ce = ce / nrm
    R = mp.re(mp.fsum(mp.conj(ce[i]) * Ev[i, j] * ce[j] for i in range(N15 + 1) for j in range(N15 + 1)))
    r = R - eps
    print(f"  FDgrid={Mfd} Egrid={Mg}: QW(k~) = {mp.nstr(R, 8)}   r = {mp.nstr(r, 6)}   "
          f"r/eps = {mp.nstr(r / eps, 6)}   sqrt(2r/gap) = {mp.nstr(mp.sqrt(2 * r / (mu1 - eps)), 4)}")
print(f"  eps = {mp.nstr(eps, 8)}   mu1(even) = {mp.nstr(mu1, 8)}")
print("\nDONE — attack run complete; freeze or fix the numbers per the drifts above.")
