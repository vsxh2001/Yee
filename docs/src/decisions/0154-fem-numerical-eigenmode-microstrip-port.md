# ADR-0154: FEM numerical-eigenmode microstrip port (the ADR-0153 B4 fork — higher-fidelity port)

**Status:** Accepted
**Date:** 2026-06-03
**Related:** ADR-0153 (the 7-brick FEM driven-sweep track this forks from — specifically the **B4
GO/fork decision point** and the **B7 v1-port ceiling**), ADR-0147 (the EM-sim-wall finding that
ranked FEM driven-sweep #1), ADR-0064 (the planar-MoM port obstruction, **provably non-binding on
3-D Whitney-1** — see ADR-0153 §Context), ADRs 0050–0061 (the shipped `yee-mom`
`NumericalCrossSection::with_quasi_tem` cross-section eigensolver this port reuses),
[[fem-driven-sweep-s21-viable]], [[step5-mixed-solver-dielectric-underweight]].

---

## Context

ADR-0153 drove the FEM driven-sweep microstrip-S21 track to its two milestones:

- **B4 (the make-or-break):** straight-microstrip ε_eff to **0.61 % of Hammerstad-Jensen** —
  8× inside the 5 % ADR-0147 GO gate, reviewer-validated as **real physics** (a non-circularity
  probe fed the port a *wrong* β and still measured ε_eff tracking the FEM propagation, not a
  match-by-construction).
- **B7 (the payoff):** a 3-pole Chebyshev microstrip filter S21 — a recognizable bandpass with the
  **correct geometric asymmetry** (depth(1.6 GHz) > depth(2.4 GHz), margin +1.47 dB) — but the
  absolute level **floors at ~−42 dB**, missing the strict Cheb mask by ~42 dB.

B7's ceiling was reviewer-confirmed as **port fidelity, not a bug**: the v1 analytic flat-E_z
`modal_e_t` only ~9 % overlaps the true FEM eigenmode, so each port radiates ~−21 dB of the
incident mode (|S21|≈0.089 per port, squared across two ports ⇒ the ~−42 dB filter floor). The
direct LU fits; the mesh and scaling are fine. The single missing piece is a **higher-fidelity
modal shape** at the port face. ADR-0153 §Consequences named the fork explicitly:
*"fail → mesh-refine vs bridge `yee-mom::NumericalCrossSection::with_quasi_tem` vs TL-de-embed
Z₀."*

**A decisive de-risk probe ran the bridge before any build commit** (the maintainer's "scope the
options first" gate). The probe mirrored the B4 straight-line setup **exactly** — same
mesh / box / interior-PEC / two-length extraction, same analytic-HJ β so the **only** changed
variable is the modal *shape* — and swapped the v1 `modal_e_t` for one sampling `yee-mom`'s shipped
quasi-TEM cross-section eigenmode (`NumericalCrossSection::with_quasi_tem` →
`e_tangential_at(x, z)`). Measured, then **independently re-verified in the box**:

| Quantity | v1 analytic flat-E_z | numerical eigenmode | |
|---|---|---|---|
| **\|S21\|** (L2) | 0.0890 | **0.7781** | **8.74× lift — floor BROKEN** |
| **\|S11\|** (L2) | 0.573 (−4.8 dB) | **0.087 (−21 dB)** | port now genuinely **matched** |
| ε_eff (phase) | 0.61 % vs HJ | 0.61 % vs HJ | unchanged ⇒ β/phase isolated; only the shape moved |

The collapse of |S11| (not merely a |S21| bump) is the signature of a **high-fidelity modal
shape** — a coincidental amplitude gain would not also match the port. The ~9 % overlap floor was
the analytic *shape*, **not an intrinsic limit**. The cross-lane dependency `yee-fem → yee-mom` is
**acyclic** (`yee-mom` depends only on yee-core/yee-mesh/yee-io (+ optional yee-cuda); no edge back
to `yee-fem`). The maintainer reviewed the measurement and picked **"Full Option-1 now."**

## Decision

Productionize the numerical-eigenmode microstrip port and re-grade the filter, as an **ordered
3-brick decomposition** (N-bricks, continuing ADR-0153's lettering), each with a
**machine-checkable gate**. The probe code (de-risk branch `feature/fem-port-numerical-probe`,
`c102d16`) is the **seed** — its cross-section builder, one-shot Arc-shared eigensolve, and frame
map are validated and promote directly from the test into `yee-fem/src`.

| # | Brick | Gate (machine-checkable) | Risk | Deps |
|---|-------|--------------------------|------|------|
| **N1** | `microstrip_port_numerical` production API in `yee-fem/src` — promote `yee-mom` dev-dep → real dep; public constructor builds the (x, substrate-normal) cross-section, solves `with_quasi_tem` once (Arc-shared across the two port closures), returns a `PortDefinition` with the numerical `modal_e_t` + analytic-HJ β | builds clean (clippy `-D warnings`, fmt); unit test: `modal_e_t` finite/nonzero + **E_z-dominant in the gap, decaying in air**; β matches `yee_layout::eps_eff` | **eng** (de-risked by probe) | ADR-0153 B1–B4 |
| **N2** | Validated straight-line |S21| **gate** with the numerical port (`fem_line_eeff_numerical_001`) — promote the probe MEASUREMENT into a PASS/FAIL gate with a **|S21| lower-bound tripwire** (the probe was non-failing by design) | `|S21| ≥ 0.6` **AND** `|S11| ≤ 0.2` **AND** `ε_eff` within 5 % of HJ; `#[ignore]`'d + `--release` gate job | **eng** (de-risked: probe = 0.778 / 0.087 / 0.61 %) | N1 |
| **N3** | Re-grade B7's 3-pole Chebyshev filter S21 with the numerical port — swap the port in `microstrip_filter_s21.rs`, re-measure vs the strict Cheb mask | **HONEST** gate: assert the measured lift + whether |S21| clears (or by how much it approaches) the `oracle_grade` mask, **AND** the asymmetry discriminator (depth(1.6) > depth(2.4)) still fires | **research-open** (the real test: 2 good ports + 3 resonators) | N1, N2 |

**Critical path** = N1 → N2 → N3. **N1+N2 are one tight, fully-de-risked increment** (the line case
the probe already proved); they ship together. **N3 is the payoff and the one remaining
research-open question**: with two high-fidelity ports, does the *filter* clear the Cheb mask, or
does a second-order effect (resonator coupling, mesh at the gaps) cap it short? N3's gate stays
**honest** — it asserts the measured truth (the lift, the mask margin, the asymmetry), never a
match-by-construction.

## Consequences

**N1+N2 achievable in ~2-4 days** (promotion of validated probe code + a tripwire gate), **high**
confidence. **N3 is ~1-1.5 weeks** with **moderate** confidence on clearing the *strict* mask:
the probe proves the *port* is now matched and high-transmission on a straight line, which is the
dominant term, but a 3-pole filter adds resonator-coupling and gap-mesh sensitivity the line does
not exercise. If N3 lands a clean in-mask S21, the FEM driven-sweep track delivers its original
goal (a validated full-board microstrip-filter EM result, ADR-0147's #1 path). If N3 lifts the
floor dramatically but stops short of the strict mask, that is itself an **honest, documented
result** (a real bandpass with a quantified mask margin) and the remaining gap is a finer-mesh /
resonator-coupling follow-on, **not** a port-fidelity wall.

**β source (N1 decision):** the probe isolated the shape by keeping analytic-HJ β. The eigensolve
*also* yields `mode.beta` (its own ε_eff was 2.91 vs HJ 3.17 — an ~8 % box-truncation gap on the
6 mm cross-section, normal for the loose-tolerance eigensolver). N1 takes **β from analytic HJ**
(validated to 0.61 % on the driven line) and uses the numerical mode **only for the shape** — the
non-circularity probe (ADR-0153 B4) showed the measured driven ε_eff is robust to a mistuned port
β regardless, so this is both defensible and the more accurate phase.

**Resource discipline (inherited from ADR-0153, standing):** every heavy cargo invocation runs
through the bounded Docker box (`scripts/yee-box.sh`, ≤14 g / 3 cpu cgroup) so a build/solve can
never OOM the host; heavy solver tests are `#[ignore]`'d + run only in dedicated `--release` gate
jobs (never the debug workspace test). **Honesty gates:** the reviewer enforces `gate_is_real`
(no tautology / match-by-construction); N2's tripwire is a real lower bound on a measured quantity;
N3's mask grade asserts the measured truth. No EM result merges until its gate genuinely passes.

**Not in scope / do NOT reopen:** the ADR-0133 FDTD cavity wall, the ADR-0064 planar-MoM port
(non-binding here), the `fem-eig-006` *eigen* modal-projection wave-port (distinct from this
*driven* track). The full Sommerfeld-integral tail and `faer::matrix_free::bicgstab` scaling
(ADR-0153 B5b) remain deferred — N1–N3 fit the existing per-ω faer complex LU in the box.

---

## References
- De-risk probe: branch `feature/fem-port-numerical-probe` (`c102d16`, base `93d3c49`),
  `crates/yee-fem/tests/microstrip_eeff.rs::fem_line_eeff_001_numerical_port` (the seed for N1/N2)
  + `crates/yee-fem/Cargo.toml` (`yee-mom` dep).
- The reused eigensolver: `yee_mom::ports::NumericalCrossSection::with_quasi_tem`
  (`crates/yee-mom/src/ports/...`, ADRs 0050–0061), validated by
  `crates/yee-mom/tests/eigensolver_microstrip_quasi_tem.rs`.
- The v1 port + B4/B7 baseline: `crates/yee-fem/src/microstrip_port.rs`,
  `crates/yee-fem/tests/microstrip_eeff.rs` (B4), `crates/yee-fem/tests/microstrip_filter_s21.rs`
  (B7, the floored filter this re-grades).
- Spec: `docs/superpowers/specs/2026-06-03-fem-numerical-eigenmode-port-design.md`;
  plan: `docs/superpowers/plans/2026-06-03-fem-numerical-eigenmode-port.md`.
