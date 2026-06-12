# AUDIT.md — adversarial hardening audit of the ExactVerification note (v1 → v2)

**Date:** 2026-06-11/12. **Mode:** paranoid-first GPD; the adversarial stance aimed at our own work.
**Scope:** every claimed contribution, every [CCM] reference number, every headline number of
`ExactVerification-SemilocalWeilForm-Note-v1.md`, before transmission to Connes–Consani–Moscovici.
**Evidence base:** (a) clean re-runs of all four certificates (logs `rerun_*.log`, all GREEN);
(b) the stability/attack run `attack_v2_hardening.py` → `attack_v2_hardening.log`;
(c) hand-trace of the full theorem numbering of `mc2arXiv.tex` (arXiv:2511.22755v1);
(d) full-source reads of arXiv:2602.04022 (Letter), 2511.23257 (CvS), 2602.15941 (Jacobian),
2606.06604 (Absolute Geometry); (e) systematic web/citation sweep for third-party verification work.

---

## 1. Redundancy audit (Mission 1): per-contribution verdicts

| # | Claimed contribution (v1 abstract) | Verdict | Evidence / required reframe |
|---|---|---|---|
| 1 | O(N) per-frequency reduction of the QW_λ matrix | **NOVEL** (implementation technique) | No prior reduction published; Groskin implements the CvS Galerkin form, different basis/object. Keep. |
| 2 | Simple-even verified at λ = 1.3, 1.5, 2.0 | **NOVEL — sharpened** | The Letter (2602.04022) restates simple-even as *open* ("The analogue of this property is known for the prolate wave operator"). No third party computed semilocal QW_λ^N spectra at λ ≤ 2; Groskin (2605.20224) works at c = λ² ≥ 13 with the even sector **imposed by basis** (parity never tested); CC's own data point is λ² = 11 ([CC-cycles] §2.5). Our claim survives as: *first at small λ, first with parity verified rather than imposed, first with full even+odd spectra*. v2 cites Groskin + [CC-cycles]. |
| 3 | First direct measurement of k_λ-to-ξ_λ proximity (Rayleigh excess r/ε) | **NOVEL (clean first) — but reframe "first"** | No measurement of r, r/ε, or the min–max transfer exists anywhere (web sweep: nothing). BUT [CCM] §8 indication (3) and [CC-cycles] §3 already report *qualitative* numerical evidence of proximity. v2 wording: the proximity is theirs; the *quantitative* Rayleigh excess + verified min–max transfer is ours. |
| 4 | Eigenvalue ladder discovery (μ₁/ε = 37 / 1.8e3 / 2e5; excited state = tower cell) | **NOVEL as measurement — reframe "discovery"** | [CCM] §8 indication (3) says proximity "extends to the higher eigenfunctions" — the tower is qualitatively known to the authors. No ladder *ratios*, no gap measurement, no Courant–Fischer bound exists in print. Nearest third-party: Śliwiński (2601.12133), 7-digit D_log^(λ,N) spectra vs zeros — different operator. v2: "quantifies the structure indicated in [CCM, §8 (3)]", cites Śliwiński. |
| 5 | Reproduction of ε_λ ≈ 1−χ(λ) over 13 decades | **KNOWN phenomenon, correctly framed as reproduction** | Their Figure 4 (label `fpro1`; v1 miscited as "Figure 5") + [CC-cycles]. Independent reproduction has verification value; keep framing. |
| 6 | Spectral realization (zeros real to 1e-44; γ₁ to 2.3e-9 with primes 2,3,4) | **PARTIALLY KNOWN — regime novel** | Phenomenon is [CCM] §6 / [CC-cycles]; Groskin reproduced it at c ≥ 13 to 1e-55…1e-334. **Nobody has touched the x < 13 / primes-{2,3} regime.** v2 cites Groskin and frames ours as the minimal-Euler-product regime + independence. |
| 7 | "Independent implementation" (blanket) | **OVERCLAIM — fix** | Groskin's arXiv:2605.20224 (May 2026, PyPI `connes-cvs`) claims verbatim "the first independent public implementation of the Connes–van Suijlekom Galerkin matrix". v2 says "independent of the authors' and of the other implementations we are aware of", citing Groskin + Śliwiński. Not citing them risks looking unaware (both cite [CCM]; Groskin cites the Letter). |

**Letter audit (does 2602.04022 supersede the missing-steps framing?):** No. It restates both steps as open
(§"Remaining steps"), proves the substrip convergence k̂_λ → Ξ **with explicit rate** c·λ^(−1/2−α)(1−2α)^(−1)
(the rate is in the Letter; [CCM] Lemma 7.3 states the convergence without a rate — v2 re-attributes), and
contains no QW_λ spectra, no Rayleigh/min–max language, no ladder. Its one numerical table is the
50-zero error table for the primes ≤ 13 Euler-product construction (different object).
**2511.23257 (CvS):** published CMP **406:312 (2025)**; purely theoretical ("no data has been created or
analysed"); it is the engine behind [CCM] Thm 1.1 — v2 adds the citation.
**2602.15941 / 2606.06604:** geometric track (Jacobian of Spec Z-bar; F₁-curve / Scholze heuristic /
Fargues–Fontaine); no bearing on either missing step; no numerics. Not cited in the note (out of scope).

**Third-party landscape (full citation sweep, 2026-06-11):** citing 2511.22755: Groskin 2605.20224
(numerical), Śliwiński 2601.12133 (numerical), Suzuki 2606.09096 (theoretical), the Letter (self).
GitHub: Siche/Architect-Resonance (targets simple-even at large λ via Krein–Rutman; unvetted, not arXiv —
noted here, not cited in v2), JeWaVe stub (no code), Waschtl904 prolate-Gram TeX, research-line
even-dominance atlas (λ ∈ [100, 1.3e6], different object). Nothing else found (searches logged in session).

## 2. Reference-number audit (Mission 2): every [CCM] citation hand-traced

CCM uses one counter per section shared across theorem/lemma/prop/cor/defn/remark; **there is no appendix**;
4 figures total. Corrections (v1 → v2):

| v1 citation | Actual item in arXiv:2511.22755v1 | Status |
|---|---|---|
| Theorem 1.1 (+ part (iii)) | Theorem 1.1 = label `finmainintro`; has (i)–(iii); "assumed simple … assumed even, normalized by δ_N(ξ)=1"; det_reg(D̃−z) = −iλ^(−iz)·ξ̂(z) | ✓ CORRECT (verbatim verified) |
| Lemma 2.2 (the q(U_n,U_m) formulas) | **Lemma 2.3** = `polarize0` (Lemma 2.2 = `qtrans` is translation invariance) | ✗ FIXED |
| Lemma 2.3 + Proposition 2.5 (Ψ♯ representation) | **Lemma 3.1** = `wsharp` (Ψ♯ = W₀,₂♯ − W_R♯ − ΣW_p♯) + **Proposition 3.2** = `toadd` (QW(κf,κg) = Ψ♯(F)) | ✗ FIXED |
| Section 8 | Section 8 "The missing steps" = `missing` | ✓ CORRECT |
| Figure 5 | **Figure 4** = `fpro1` (4 captions in the paper; fpro1 is the last, in §8) | ✗ FIXED |
| Lemma A.2 (Meixner–Schäfke O(λ⁻²)) | **Lemma 7.2** = `meixnerlem` | ✗ FIXED |
| Lemma A.3 (k̂_λ → Ξ on substrips) | **Lemma 7.3** = `hermfact1`; stated **without** rate — explicit rate c·λ^(−1/2−α)(1−2α)^(−1) is in the Letter §6.5 | ✗ FIXED + rate re-attributed |
| §1 blockquote of the two missing steps | v1 was a **paraphrase** presented as quote (invented "(i)/(ii)" markers; "approximation **of**" vs their "**to**"; dropped the Theorem-3.6 clause) | ✗ FIXED — v2 quotes verbatim |
| (uncited) existence of smallest eigenvalue | Theorem 3.6 = `thmsmallest` | added where the verbatim quote needs it |

Notation spot-checks all verbatim-correct: PW_λ = −∂_x((λ²−x²)∂_x)+(2πλx)²; 𝓔(f)(u) = u^(1/2)Σ_{n≥1}f(nu);
χ(λ) = Fourier-compression eigenvalue; k_λ = 𝓔(h_λ), h_λ the vanishing-integral combination of h_{0,λ}, h_{4,λ};
Fuchs n=4 asymptotic (2^14/3)√2·π⁵·e^(−4πλ²+9logλ) as quoted by CCM from Fuchs Thm 1.
References block verified against CCM's own bibliography: [CC-cycles] Enseign. Math. 69 (2023) no. 1–2, 93–148 ✓;
[Co-zeta] Selecta Math. (N.S.) 5 (1999) no. 1, 29–106 ✓; [CM] PNAS 119 (2022) no. 22 ✓; Fuchs JMAA 9 (1964)
317–330 ✓; Slepian–Pollak BSTJ 40 (1961) 43–63 ✓; Meixner–Schäfke Springer 1954 ✓. Added in v2: CvS CMP
406:312 (2025); the Letter arXiv:2602.04022; CC Weil-positivity Selecta 27 (2021) completed with no. 4, Paper No. 77.

## 3. Adversarial audit of headline numbers (Mission 3)

All four certificates re-ran GREEN from clean state (2026-06-11, logs kept). Attack runs:
`attack_v2_hardening.py` (ATTACK 1–4). Per-number verdicts:

**(a) Overlap "1.0000" (excited state = tower cell) — SURVIVES, with digits.**
Even-sector projection of the first excited even eigenvector onto span{𝓔h₀, 𝓔h₄, 𝓔h₈} at λ=2:
**0.999997 (deficit 3.187e-6)**; FD grid 4001→8001 and 𝓔-grid 3001→6001: deficit 3.184e-6 (stable);
N=12→14: deficit 2.824e-6 (stable). Control: 20 random 3-dim subspaces of the 13-dim even sector give
mean projection 0.399, max 0.760 — the cell-span value is five nines above chance. NOT a projection-rank
artifact. **Convention finding (fired falsifier F4):** v1's companion number "0.414" is the *full-space,
raw-cell* convention; the even-sector value is **0.593**. Both ≪ 1 (h₈-content essential — conclusion
unchanged); v2 states the convention and gives both.

**(b) γ₁ to 2.4e-9 / zeros real to 1e-44 — SURVIVES, strengthened; one rounding fix.**
Root-finder basin: 4/4 distinct starts converge to the same root for each of ẑ₁, ẑ₂, ẑ₃ (Re spread ≤ 3e-44).
Precision-floor test: at dps=60 the imaginary parts fall 1e-44 → **1e-66** while Re is unchanged to 46
digits — the reality claim is a genuine precision floor, machine-verified. N-dependence: |ẑ₁−γ₁| = 2.346e-9
at N=12, 1.889e-9 at N=14 (zeros drift mildly with truncation; table is N-labelled in v2).
Rounding fix: 2.346e-9 → quote **2.3e-9** (v1's "2.4e-9" was mis-rounded).

**(c) Ladder 37/1787/2e5 and ε vs Fuchs 3.3/7.7 — STRUCTURE SURVIVES; digits required N-labels (fired falsifier F3).**
Quadrature 16→32 panels: zero drift (10 digits). dps 40→60: drift ≤ 1e-29. **But N-refinement at λ=2.0:**
ε = 1.375e-12 (N=12) → 1.108e-12 (N=14) → 1.056e-12 (N=16): **19–23% drift — ε is NOT N-converged at the
quoted truncation**; μ₁/ε = 2.00e5 → 2.45e5 → 2.33e5; odd min 8.75e-10 → 6.47e-10. At λ=1.5: ε drifts 2.6%
(N+2), 4.5% (N+4); μ₁/ε = 1787 → 1850. At λ=1.3: 0.3%. Fuchs ratios: λ=1.5: 3.28 → 3.13; λ=2.0: 7.68 →
6.19 → 5.90. **Verdict:** the simple-even ordering, the ladder structure, and order-unity Fuchs tracking
hold at every truncation tested (and the finite sections are the actual objects of [CCM] Thm 1.1); the
4-digit values are properties of QW_λ^N at stated N, not of the λ-limit. v2 adds N to every table and a
drift-disclosure paragraph; "7.7" is requoted as "7.7 at N=12, decreasing to 5.9 at N=16 — order-unity
tracking is the stable statement".

**(d) r/ε = 0.135 — SURVIVES; FD budget computed; v1 internal inconsistency fixed (fired falsifier F2).**
FD grid 2001/4001/8001 × 𝓔-grid 3001/6001: r/ε = 0.1335 / 0.1355 / 0.1367 / 0.1355 / 0.1366 — total spread
2.4%, drifting *up* with refinement (FD-limit ≈ 0.137). Each QW(k̃) is the exact Rayleigh value of an
explicitly constructed vector (hence individually a rigorous upper bound for ε); FD uncertainty affects only
the identification of k̃ with the ideal prolate ansatz, at the 2% level. **v1 bug:** §5 quoted the
inconsistent triple QW(k̃) = 1.91e-4 with ε = 1.567e-4 and r/ε = 0.135 (1.91e-4 matches neither the
even-projected run, 1.779e-4, nor the unprojected one, 1.953e-4). v2 quotes the consistent even-sector
triple: QW(k̃) = 1.779e-4, r = 2.12e-5, r/ε = 0.135 (N=10 window), with the grid-refinement budget.
Also fixed: v1's "√(r/(μ₁−ε)) … 1.2e-2" — that value is √(2r/(μ₁−ε)).

**(e) Smaller fixes found by cross-checking:** abstract said "five certificates", §7 said four (now five,
including the stability audit); §7 listed `sympy` (no probe imports it); "three orders below the λ⁻² floor"
→ "2.5–3 orders" (measured ratios 9.7e-4–5.0e-3); diagonal validation claim tightened to the two entries
actually checked; "primes in window" column header → "prime powers".

### Fired falsifiers (retained per GPD; recorded in README.md)
- **F1 (Mission 2):** "internal [CCM] labels inferred from TeX labels are correct" — FIRED: no appendix
  exists (A.2/A.3 → 7.2/7.3); "Prop 2.5" → Prop 3.2 + Lemma 3.1; "Lemma 2.2" → Lemma 2.3; "Figure 5" → Figure 4;
  the §1 "quote" was a paraphrase.
- **F2 (Mission 3d):** "the §5 triple is internally consistent" — FIRED (1.91e-4 ≠ ε(1+r/ε) under any run).
- **F3 (Mission 3c):** "ε_N is converged at the quoted N" — FIRED at λ=2.0 (19–23% drift to N+4); structure robust.
- **F4 (Mission 3a):** "0.414 is the even-sector first-cell projection" — FIRED (convention-dependent: 0.593
  even-sector vs 0.414 full-space raw-cell).
- **F5 (Mission 1):** "this is the first independent implementation" — FIRED as blanket claim (Groskin May 2026);
  survives in the four sharpened forms above.

## 4. Mission 4 — elevation fruit (ranked; default conservative)

Their open questions swept: [CCM] §8 (two steps + three indications), Letter §"Remaining steps" + §6.5
(no numbered "Open Question 1" exists in the Letter — the memory shorthand was imprecise), CC-2021 Thm 3
(semilocal extension unproven), [CC-cycles] §§2.5/3, CM-PNAS. Cross-referenced against the certified
inventory (P∞-E §§11–26, P∞-F, N-series).

1. **Radical-membership + exact tail identity for k_λ** (P∞-E §25, Lemmas E.22.1/E.22.2 [P, certified]):
   *their* step (ii) mechanism — the two defining conditions of h (∫h = 0, evenness) kill exactly the two
   pole terms of the Weil functional ⟹ global 𝓔-image lies in the radical of the global form ⟹ windowed
   Rayleigh numerator = Weil form of the below-window tail (linear in 1−χ₄, certified). Effort: 2–3 pages.
   **Recommendation: one-sentence trailer in v2 §8 (done); full statement in a SECOND note after contact.**
2. **Quantified weld theorem** (P∞-E §22, Thm E.20 [P as conditional]): (H1)+(H2)+(H3) ⟹ RH with explicit
   poly transfer constants — the precise conditional form of their §7 strategy; the note already verifies its
   inequalities on data. Effort: small. **Recommendation: second note (same package as #1).**
3. **The named gap bridge (B2)** (P∞-E §23, Thm E.21(b)): μ₁(QW_λ) ≥ (1−χ₈)/poly(λ) is the exact missing
   inequality connecting their indication (1) to step (ii), with c⁴ room measured. Effort: zero (one remark).
   **Recommendation: kept implicit in v2 §4 (Courant–Fischer paragraph); state explicitly in second note.**
4. **Off-tower sector = zero-representer ladder; flat floor 2θ′(γ₁)** (P∞-F §§4–7 [measured + partial proofs]):
   the semilocal quantitative face of CC-2021's Sonin positivity. Deeper program material; needs CC-2021
   vocabulary alignment. **Recommendation: HOLD (third communication / paper).**
5. **N14 "Gibbs/relative-entropy reason for archimedean positivity":** EXCLUDED — our own 2026-06-08
   novelty audit found it textbook (Amari–Nagaoka; Garbouge–Nielsen) and the Connes quote it answers
   mis-deployed (he supplies his own operator-theoretic reason, §7.2 Sonin space). Do not send.

**Verification note stays a verification note.** Items 1–3 are flagged in the §8 closing remark only.

## 5. Number freeze (v2 canonical values)

| Quantity | Frozen value | Source |
|---|---|---|
| ε(λ=1.3, N=10) | 2.4634e-2 (N+2 drift 0.3%) | attack log |
| ε(λ=1.5, N=10) | 1.56653e-4 (N+2: −2.6%, N+4: −4.5%) | attack log |
| ε(λ=2.0, N=12) | 1.3754e-12 (N+2: −19%, N+4: −23%) | attack log |
| μ₁/ε ladder | 36.7 (λ=1.3, N=10) / 1787 (λ=1.5, N=10; 1850 at N=14) / 2.00e5 (λ=2.0, N=12; 2.3–2.5e5 at N=14,16) | attack log |
| odd minima | 0.3746 / 1.4407e-2 / 8.749e-10 (at stated N; λ=2 N+4: 6.47e-10) | attack log |
| Fuchs ratios | 3.28 (λ=1.5, N=10 → 3.13 at N=14); 7.68 (λ=2.0, N=12 → 5.90 at N=16) | attack log |
| overlap deficit | 1 − P = 3.2e-6 (grid- and N-stable); first-cell: 0.593 even-sector / 0.414 full-space | attack log |
| r/ε (λ=1.5, N=10) | 0.135 at baseline grids; 0.133–0.137 across refinements (FD-limit ≈ 0.137) | attack log |
| weld triple | QW(k̃) = 1.779e-4, r = 2.12e-5, ‖k̃−ξ‖ = 0.0066 ≤ √(2r/gap) = 0.0123 | rerun_weld.log |
| zeros (λ=2.0, N=12) | |ẑ−γ| = 2.3e-9 / 6.5e-7 / 1.4e-5; |Im| ≤ 8e-44 at dps 40, ≤ 3e-63 at dps 60; basin 4/4 | attack log |
| quadrature / dps sensitivity | 0 (10 digits) / ≤ 1e-29 | attack log |
| Poisson mirror | 1.37e-15; h-identity 1.94e-16; duality scalar 4.000000 | rerun_seam.log |
| M–S rates | λ^−2.16 (n=0), λ^−2.56 (n=4); k̂→Ξ end-to-end ≈ λ^−2.3 | rerun_seam.log |
