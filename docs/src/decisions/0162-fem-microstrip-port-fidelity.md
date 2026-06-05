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

## De-risk outcome (B1 + B1.5 — DECISIVE; branch `feature/fem-port-power-norm`, B1 `c2f4f97`, B1.5 `cb84ab0`)

Self-verified boxed on the matched straight thru (lossless FR-4, PEC metal — physically lossless):

| measurement | result |
|---|---|
| `\|S11\|²+\|S21\|²`, current **E-only L²** extraction | **0.6145** (reproduces the smoking gun) |
| `\|S11\|²+\|S21\|²`, **quasi-TEM √ε_r power-norm** (B1) | 0.6725 (lifts only +0.06 — the approximate modal-H *shape* barely helps) |
| **P_out/P_in, true-field Poynting flux** `½Re∫(E×H*)·n̂`, H=∇×E/(−jωμ) (B1.5) | **0.9982** |
| `\|S21\|²/(1−\|S11\|²)` (same field, via the S-formula) | 0.6109 |

**Conclusion — the B1-only NO-GO is OVERTURNED:** the solved FEM field **conserves energy** (transmits
0.998 of the incident power port-to-port — lossless to 0.2 %, the per-port fluxes near-equal-and-opposite
`[−1.714e-10, +1.711e-10] W`, the lossless-2-port signature). The ~39 % deficit is therefore a **pure
EXTRACTION-normalization artifact**, NOT real solver/ABC/numerical loss. The K3 Q-floor is **not** what
floors a *non-resonant* thru. **⇒ the filter-S21 floor is an EXTRACTION problem (the cheap, tractable
side of the fork), not a solver wall — the track is SALVAGEABLE.**

B1.5's true H came from a new `pub(crate) element::tet_whitney_e_and_curl` (reuses the assembly's exact
`∇×N_α = 2∇λ_i×∇λ_j` Whitney-1 curl — no re-derivation) + `OpenBoundarySolver::poynting_flux_audit`.

**This reframes B2** (supersedes the table's B2): NOT "productionize the quasi-TEM √ε_r power-norm" (B1
showed it barely helps — its modal-H *shape* is approximate), but a **flux-calibrated extraction** —
unit-incident-power-normalize the modal projection against the **true** modal power flux (the
Palace/COMSOL recipe Yee skipped: normalize the modal field so `½Re∫(e_m×h_m*)·ẑ = 1` first, with an
accurate modal H — e.g. from a reference thru-solve's `∇×E`, which B1.5's evaluator now provides — then
`a_m = ∫(E_FEM×h_m*)·ẑ`). Target: the thru `|S21| → ~1` (matching the 0.998 flux) + `|S11|²+|S21|² → 1`,
with ε_eff/β unchanged (N2/B4 stay green). Then B3' applies the corrected extraction to the **filter**
(heavy boxed) + a flux audit on the filter passband → does the in-band |S21| lift toward the mask? (If
the filter field also transmits in-band, the −27 dB was largely an extraction artifact and the mask may
be reachable; if the filter field genuinely reflects in-band, that residual is real.)

## Final outcome (B2'→B3'' — SHIPPED, merge `bb9608a`; the filter floor is COUPLING, not extraction)

**B2' (the extraction fix — REAL WIN):** a power-correct **E+H modal wave-port extraction** —
`a_fwd/a_bwd = ½(projE ± projH)`, `projE = ∫(E_FEM×h_m)·ẑ/κ`, `projH = ∫(e_m×H_FEM)·ẑ/κ`,
`κ = ∫(e_m×h_m)·ẑ` (un-conjugated, the fwd/bwd split via the H sign flip; Jin/COMSOL recipe), with the
true modal H from `∇×E/(−jωμ)`. On the matched thru it recovers `|S21| 0.778 → 1.0001`,
`|S11|²+|S21|² 0.6145 → 1.0037`, ε_eff 0.66 % (β unchanged). The old E-only L² normalization
under-counted a *transmitting* port by ~40 %; this is a genuine S-parameter-correctness fix.

**B3'→B3'' (the filter test — the honest reversal):** the in-situ modal reference (B3') gave a *spurious*
in-band peak −0.86 dB (+26 dB) but an **UNPHYSICAL** curve (`|S21|→2.96`, `|S|²sum→12.7`) — standing-wave
contamination of the in-situ reference at the reflective filter feeds; the honest gate (passivity) caught
it. A **clean modal basis** from a matched reference line + a port-FACE E&H projection (B3'') makes the
curve PHYSICAL (`|S21|≤0.05`, in-band `|S|²sum 0.79–0.93`) and reveals the **true** filter: in-band peak
**−26.14 dB (only +1.24 dB over N3)**, in-band **`|S11| ≈ 0.85–0.91`** — the filter **GENUINELY REFLECTS
in-band** (transmits ~0.2 %), mask MISS by ~33.7 dB (asymmetry PASS +2.17 dB).

**DECISION / honest conclusion:** the −27 dB filter floor is **REAL coupling-bound reflection** (the
resonators do not realize the Chebyshev in-band match), **NOT an extraction or port-fidelity artifact**.
The B3' +26 dB was an unphysical contamination false-positive (recorded as such; no fake claimed). The
funded port-fidelity hypothesis is **refuted as the filter floor**. **What this track delivered:** (a) a
correct power-conserving E+H wave-port extraction (`power_modal_extract`; thru `|S21|→1`, reusable —
though NOT yet wired into the production `extract_s_qp`, a deferred follow-on); (b) reusable diagnostics
(`poynting_flux_audit` + `element::tet_whitney_e_and_curl`); (c) the **definitive diagnosis** — the filter
floor is COUPLING (ADR-0159 territory; dimensioning lifted only +5.8 dB), corroborated by N3's independent
E-only −27 dB and the physical in-band `|S|²sum`. The reviewer (APPROVE-WITH-FIXES, no P0/P1) confirmed the
extraction math, that the reflection conclusion is sound (not an x-shift artifact — symmetric mirror feeds
cancel in S21), and the gate honesty. **A mask-clearing full-wave filter S21 needs the COUPLING addressed
(why the FEM filter mismatches in-band) — a different track from port-fidelity; the port/extraction is now
exonerated as the floor.**

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
