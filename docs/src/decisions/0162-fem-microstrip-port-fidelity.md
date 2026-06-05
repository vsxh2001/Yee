# ADR-0162: Higher-fidelity FEM microstrip port — power-wave normalization + multi-mode (the filter-S21 floor)

**Status:** Accepted (track kickoff; maintainer-funded 2026-06-05)
**Date:** 2026-06-05
**Related:** ADR-0147 (#1 goal: a mask-clearing full-wave filter S21), ADR-0153/0154 (FEM driven-sweep +
the N1/N2/N3 numerical-eigenmode port), ADR-0159 (B2: dimensioning is a minor lever → port-bound),
[[fem-driven-sweep-s21-viable]], [[project-filter-design-final-goal]].

---

## Context

The FEM driven-sweep is validated for ε_eff (B4, 0.61 % of Hammerstad-Jensen) and coupling-k (K1/K2),
and the numerical-eigenmode microstrip port (ADR-0154) is matched + high-transmission on a **straight
thru** (`|S21|≈0.778`, `|S11|≈0.087`). But the 3-pole edge-coupled **filter** S21 floors at ≈ −27 dB
(N3) / −21.5 dB (B2), missing the Chebyshev mask. ADR-0159 (B2) showed gap-dimensioning lifts it only
+5.8 dB → the dominant floor is **port fidelity**, not dimensioning. The maintainer funded a
higher-fidelity-port track (2026-06-05).

**Research + diagnosis (two read-only agents, research-first — see the source list):** the floor has a
specific, testable root cause.

1. **The S21 extraction uses an E-field L² normalization, not power-wave normalization.**
   `OpenBoundarySolver::extract_s_qp` / `extract_s11` (`crates/yee-fem/src/open_boundary.rs` ~1853–2124)
   compute `S_qp = ⟨E_FEM,t, e_mode_q⟩_port / M_qq − a_inc`, with `M_qq = Σ_face w_g (e_mode·e_mode)` —
   the real E-field self-overlap `∫|e_mode|²`, NOT the modal power `Re ∫(e_mode × h_mode*)·ẑ`. The
   `PortDefinition` carries only `modal_e_t` + a scalar `beta_mode` — **no modal H-field, no per-point
   wave impedance**.

2. **Smoking gun:** even on the matched, lossless straight thru, `|S11|²+|S21|² ≈ 0.087²+0.778² ≈ 0.61`
   — ~39 % of incident power is unaccounted-for. A correct power-conserving extraction gives `≈ 1` for a
   lossless 2-port. So the extraction is **not power-unitary**.

3. **Why a thru is "OK" but a filter floors.** Microstrip is **inhomogeneous** (air + dielectric), so the
   modal wave impedance varies across the cross-section; the E-only norm mis-weights the power. On a
   uniform mode-matched thru the S-parameter *ratio* survives (≈0.78), but on the filter (a) each coupling
   gap launches higher-order / evanescent modes carrying real power that the **single-mode** projection
   silently discards, and (b) the inhomogeneous-Z mis-weighting compounds. Result: the extracted
   fundamental-mode |S21| collapses.

4. **Caveat (must be separated in the de-risk):** the thru's 39 % deficit may be *partly* the known ~30 %
   numerical Q-floor (ADR-0156/K3: `|S11|@res` coupling-invariant ≈0.84). The de-risk B1 must distinguish
   an **extraction-normalization artifact** (fixable by power-normalization) from **real numerical loss**
   (a different, harder problem).

## Decision

Pursue the port-fidelity track in **de-risk-first bricks**, each machine-checkable. The diagnosis suggests
the first fix may be far cheaper than "multi-week" (a power-normalization correction), so B1 measures
before B2+ commits.

| # | Brick | Gate / decision |
|---|-------|-----------------|
| **B1** | **De-risk (decisive, cheap — on the THRU, not the heavy filter).** (a) Print `\|S11\|²+\|S21\|²` for the thru + (the existing) filter sweep. (b) Implement a **power-normalized re-extraction** of the SAME solved thru field: obtain the modal H (from the cross-section eigensolver if it exposes `h_t`, else the quasi-TEM relation `h_t = ẑ×e_t · Y(x,y)` with the per-point modal admittance, or `h = ∇×E/(−jωμ)`), normalize by `Re ∫(e×h*)·ẑ`. | If the thru's `\|S11\|²+\|S21\|²` rises `0.61 → ~1` under power-normalization ⇒ **the L² normalization is the bug** ⇒ proceed to B2. If it stays ≈0.61 ⇒ the deficit is **real numerical loss** (the K3 Q-floor) ⇒ NO-GO for the normalization fix; re-scope to the Q-floor (honest documented outcome). Decisive + cheap (line solve, not the 40-min filter). |
| B2 | **Power-wave normalization** (if B1 implicates the norm). Add `modal_h_t` (or a per-point wave impedance) to `PortDefinition`; replace the L² `M_qq` in `extract_s_qp`/`extract_s11` with `Re ∫(e×h*)·ẑ`; keep the excitation reciprocal. | The thru gate: `\|S11\|²+\|S21\|²→~1`, `\|S21\|`/ε_eff stay within their N2/B4 tolerances (no regression). Then re-run the **filter** (heavy, boxed): does the in-band peak lift toward the mask? Record honestly. |
| B3 | **Multi-mode port** (if B2 leaves a floor — power leaking into unmodeled modes at the coupling gaps). Add higher-order cross-section modes as `a_inc=0` absorbing projectors (`PortDefinition` already supports `Vec<PortMode>`) + sum per-mode power in the extraction. | The filter in-band peak lifts further / clears the mask; per-mode power accounting closes `\|S11\|²+\|S21\|²→1` on the filter. |
| B4 | **Hardening / fallback** (only if needed): thru-line de-embed (scikit-rf `Network.inv`/`**` style — a systematic launch-error removal, NOT a floor fix) and/or a PML-backed absorbing port for residual higher-order leakage. | Filter S21 robust to the launch / reference-plane; documented. |

**Start B1** — the cheapest experiment that determines the entire track's direction (and whether it is
even a normalization problem vs the numerical Q-floor). Misfire-split: an agent writes the code, the
orchestrator runs the (cheap) thru de-risk boxed.

## Consequences

- If B1+B2 confirm + fix the normalization, a mask-clearing (or much-closer) full-wave filter S21 may be
  reachable **without** the full multi-week effort — a large win for the ADR-0147 #1 goal.
- If B1 shows the deficit is the numerical Q-floor, that is an honest NO-GO for the port-normalization
  fix and re-scopes the track (the Q-floor is the real wall) — documented, not faked.
- Scope: `crates/yee-fem/src/open_boundary.rs` (the extraction + `PortDefinition`),
  `crates/yee-fem/src/microstrip_port_numerical.rs` (the port builder, modal H), possibly
  `crates/yee-mom/src/ports.rs` (expose the cross-section modal H), + the line/filter gates.
- **Not in scope / do NOT reopen:** the dimensioning lever (ADR-0159, minor), FDTD-resonant (ADR-0108),
  mom-002/003 (ADR-0064), fem-eig-006.

## Method references (research-first — implement, don't reinvent)
- **Modal power normalization / wave-port S-extraction:** Jin, *The Finite Element Method in
  Electromagnetics* (wave-port chapter); COMSOL RF "S-Parameter Calculations" (power-flow normalization,
  conjugate-mode overlap, frequency-dependent unless TEM); arXiv 2407.21766 (modal power
  `κ_m = ∫(e_m×h_m*)·ẑ`, coefficient `α_i = ∫(E_tot×h_i*)·ẑ / κ_i`); Palace boundaries
  (`S_ij = ∫E·E_inc/∫E_inc·E_inc − δ`, valid only because its wave-port mode is unit-incident-power
  normalized first — the step Yee skipped).
- **V/I de-embedded microstrip port (the canonical OSS recipe):** openEMS `AddMSLPort.m` +
  `calcTLPort.m` (3 voltage planes + 2 current loops → differential `Z0`/`β` à la Gwarek, IEEE MGWL
  6(5):187 1996; telegrapher de-embed `U′=u·cos βΔ − j i Z0 sin βΔ`; incident/reflected split). Sheen-Ali-
  Abouzahra-Kong, IEEE T-MTT 38(7):849 1990 (the V/I recipe the Yee comments already cite).
- **De-embedding:** scikit-rf `skrf.calibration` (TRL, IEEEP370 2x-thru, `SplitPi/Tee`) — hardening, not
  the primary fix (won't lift an extraction/multi-mode floor; leaves mismatch ripple).
- **Multi-mode / absorbing ports:** arXiv 2407.21766 (number of modes per port); HFSS/COMSOL numeric-port
  multi-mode; the repo's own ADR-0049/0070 absorbing-complement (note: first-order barely moved |S11| on
  WR-90 — treat PML-backed as a refinement, not the lead).

## References (code)
- Extraction + `PortDefinition`: `crates/yee-fem/src/open_boundary.rs` (`extract_s_qp` ~1853,
  `extract_s11` ~1970, `M_qq` self-inner ~1927–1953, `PortDefinition` ~556–725, RHS
  `scatter_port_face_gauss` ~2410).
- Numerical port builder: `crates/yee-fem/src/microstrip_port_numerical.rs` (`single_mode` return ~334).
- Cross-section mode (modal H source): `crates/yee-mom/src/ports.rs` (`e_tangential_at` ~689).
- Gates: thru `crates/yee-fem/tests/microstrip_eeff.rs` (N2/B4); filter
  `crates/yee-fem/tests/microstrip_filter_s21.rs` (N3/B2).
