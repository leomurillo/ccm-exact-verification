# Exact-arithmetic verification of the CCM semilocal Weil form $QW_\lambda$

Independent implementation and exact-arithmetic verification of the semilocal Weil quadratic form
of Connes–Consani–Moscovici, *Zeta spectral triples*, arXiv:2511.22755 — bearing on the two missing
steps of their Section 8 (the simple-even condition; the prolate-ansatz proximity).

**Note:** `ExactVerification-SemilocalWeilForm-Note-v2.md` (v1 retained for the record).
**Audit trail:** `AUDIT.md` — per-claim adversarial audit performed before transmission.

## Certificates

All scripts are deterministic and self-contained; stock Python ≥ 3.11 with `numpy`, `scipy`, `mpmath`.

| script | contents | runtime |
|---|---|---|
| `probe_ccm_exact_matrix.py` | per-frequency exact matrix; zero-side validation; even/odd spectra at λ = 1.3, 1.5, 2.0; ε vs Fuchs; Rayleigh excess | ~3 min |
| `probe_quantified_weld_e2e.py` | min–max transfer; Mellin/strip transfer constants; spectral realization (zeros vs γ₁, γ₂, γ₃) | ~4 min |
| `probe_fuchs_second_rung.py` | tower identification of the excited state; Courant–Fischer gap bound; Fuchs rung ratios | ~4 min |
| `probe_ccm_seam_weld.py` | Meixner–Schäfke rates; k̂_λ → Ξ end-to-end; exact Poisson mirror and finite-λ breakage | ~3 min |
| `attack_v2_hardening.py` | stability audit: N-refinement, quadrature doubling, 40 vs 60 digits, overlap deficit + random-subspace control, root-finder basin, FD budget for r/ε | ~15 min |

Recorded outputs: `rerun_*.log`, `attack_v2_hardening.log`.

## Failures retained (fired falsifiers from the pre-send audit)

Per our working method, falsifiers that fired are recorded, not erased:

1. **Reference-number falsifier (fired).** v1 cited [CCM] item numbers partly inferred from TeX labels.
   Hand-tracing the actual per-section counters showed: there is **no appendix** in 2511.22755
   (v1's "Lemma A.2/A.3" are Lemmas **7.2/7.3**); the Ψ♯-representation is Lemma **3.1** + Proposition
   **3.2** (not "Lemma 2.3 and Prop 2.5"); the q(U_n,U_m) formulas are Lemma **2.3** (not 2.2); the
   ε ≈ 1−χ graph is Figure **4** (not 5); and v1's Section-8 "quote" was a paraphrase. All fixed in v2,
   which quotes Section 8 verbatim.
2. **Internal-consistency falsifier (fired).** v1 §5 quoted QW(k̃) = 1.91e-4 together with r/ε = 0.135 —
   numbers from two different runs (even-projected: QW = 1.779e-4, r/ε = 0.135; unprojected: 1.953e-4,
   0.247). v2 quotes one consistent even-sector triple.
3. **N-convergence falsifier (fired).** The assumption "ε is N-converged at the quoted truncation" is
   false at λ = 2.0: ε drifts −19% (N=14) and −23% (N=16); the Fuchs ratio 7.7 → 5.9. The simple-even
   ordering and the ladder structure hold at every truncation; v2 labels every number with its N and
   discloses the drifts.
4. **Convention falsifier (fired).** The "first-cell projection 0.414" of v1 is the full-space raw-cell
   convention; the even-sector value is 0.593. Both ≪ 1 (the h₈-content conclusion is unchanged);
   v2 states the convention.
5. **Priority falsifier (fired).** "First independent implementation" as a blanket claim is false:
   Groskin (arXiv:2605.20224) independently implemented the Connes–van Suijlekom Galerkin matrix at
   c = λ² ≥ 13 (even sector imposed); Śliwiński (arXiv:2601.12133) computed low-precision D_log spectra.
   v2 cites both and claims only what survives: first at the semilocal windows λ ≤ 2, first with parity
   verified rather than imposed, first measurement of the Rayleigh excess r/ε and of the eigenvalue
   ladder, first reproduction of the spectral realization in the primes-{2,3} regime.

What survived the audit unweakened: the simple-even verification at all three windows and all tested
truncations; the reality of the zeros (imaginary parts fall 10⁻⁴⁴ → 10⁻⁶⁶ when precision is raised
40 → 60 digits — a genuine precision floor; root-finder basin 4/4); γ₁ to 2.3×10⁻⁹ from primes
{2,3,4}; the tower identification (squared-projection deficit 3.2×10⁻⁶, five nines above the
random-subspace control); r/ε = 0.135 (0.133–0.137 across grid refinements); quadrature and precision
sensitivity zero at the quoted digits.

## Reproducing

```
pip install numpy scipy mpmath
python probe_ccm_exact_matrix.py
python probe_quantified_weld_e2e.py
python probe_fuchs_second_rung.py
python probe_ccm_seam_weld.py
python attack_v2_hardening.py
```

Each script is a standalone certificate: it prints its checks and a final GREEN/RED verdict, and is
deterministic (no randomness beyond seeded controls). Total runtime ≈ 30 minutes on a laptop.
Recorded reference outputs are in `rerun_*.log` and `attack_v2_hardening.log`.

## License

- Text (the notes, `AUDIT.md`, this README): [CC BY 4.0](LICENSE-CC-BY-4.0.txt)
- Code (`probe_*.py`, `attack_v2_hardening.py`): [MIT](LICENSE-MIT.txt)

## Links

- Target paper: [arXiv:2511.22755](https://arxiv.org/abs/2511.22755) (Connes–Consani–Moscovici, *Zeta spectral triples*)
- Archived release: Zenodo DOI — *to be filled at publication*
