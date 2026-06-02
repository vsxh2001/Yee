# ADR-0153: FEM driven-sweep microstrip-S21 EM-verification track (7-brick decomposition)

**Status:** Accepted
**Date:** 2026-06-02
**Related:** ADR-0147 (the EM-sim-wall finding that ranked FEM driven-sweep **#1**, blocker-free),
ADR-0064 (the planar-MoM port obstruction that provably does **not** bind 3-D FEM), ADR-0133
(the FDTD cavity wall — a *different* wall), ADR-0108 (the microstrip z-stack), the
`fem-eig-006` memory (the *eigen* modal-projection port — distinct from this *driven* track),
[[fem-driven-sweep-s21-viable]], [[emwall-oracle-validation]]

---

## Context

ADR-0147 decomposed the full-board EM-sim wall and ranked **FEM freq-domain driven-sweep #1**:
the only direction judged genuinely blocker-free (no cavity-Q wall, no single-cell aperture-port
floor) on the real 3-D microstrip, at ~2-3 weeks of engineering plus an LU-scaling risk. The
maintainer picked this direction ("continue 1 until done").

A **resource-bounded scoping workflow** (`fem-em-blocker`, 8 agents: 5 read-only de-risk +
synthesis + build + adversarial review; only the read-only phase fanned out, all cargo serial +
boxed at 6 g / 2 cpu, host monitored — never OOM'd) confirmed the path and the three reusable
primitives already on `main`:

1. **Interior-PEC assembly** — `assemble_complex_with_pec_edges` (`assembly.rs:442`) already reads
   an arbitrary global-edge `HashSet` and is already on the driven path (`open_boundary.rs:1297`).
2. **Per-tet ε_r threading** — `assembly.rs:484` `db.eps_at(tag, ω)`, proven by
   `dispersive_solve.rs` centroid-repaint.
3. **Port contract = plain closures** — `beta_mode` / `modal_e_t` are `Box<dyn Fn>`
   (`open_boundary.rs:590/595`), so a quasi-TEM microstrip port injects an **analytic** HJ β + an
   E_z-dominant modal field with **zero `yee-mom` dependency** — collapsing the cross-lane risk
   three scope reports flagged.

Decisively, **the ADR-0064 obstruction does not bind here.** ADR-0064 is planar-MoM-specific: the
in-plane RWG surface-current basis cannot represent the microstrip quasi-TEM mode's dominant
substrate-normal `E_z`. `yee-fem` is first-order **Nédélec/Whitney-1 on tets with true 3-D edge
DoFs** — vertical edges carry `E_z` as a first-class field. The scatter projection
(`open_boundary.rs:2390`) drops the propagation component and projects the in-(y,z)-face components
onto edges including the vertical ones. This is exactly why FEM can do what planar MoM provably
cannot.

## Decision

Pursue the microstrip-S21 verification as an **ordered 7-brick decomposition**, each with a
**machine-checkable gate**. Earliest bricks de-risk the open questions (quasi-TEM port modal
closure, solver scaling) on the **smallest geometry**. The **ADR-0147 milestone — straight-line
ε_eff within 5 % of Hammerstad-Jensen (B4)** — is the explicit **GO / research-fork decision
point, BEFORE any filter**. The final gate is a 3-pole Chebyshev microstrip-filter S21 graded vs
`yee_filter::ladder_s21`, **including the geometric-asymmetry discriminator** (the 1.6 GHz notch
deeper than the 2.4 GHz notch — the honest check that the curve is real physics, not a fitted
artifact).

| # | Brick | Gate (machine-checkable) | Risk | Deps |
|---|-------|--------------------------|------|------|
| **B1** | Wire interior-PEC edges into `OpenBoundarySolver` + geometric edge-picker | interior-DoF count drops by **exactly** `|E|`; idempotent re-tag; **SHIPPED in this merge** | eng | — |
| B2 | Layered straight-microstrip tet mesh (substrate+air+trace+ground, per-tet ε_r) | tet count + substrate/air tag proportions; trace/ground pickers non-empty; `eps_at(1)=4.4` | eng | B1 |
| B3 | Quasi-TEM port: analytic HJ β(ω) + E_z-dominant `modal_e_t` closures (no `yee-mom` dep) | β matches `yee_layout::eps_eff` to <1e-9; mode E_z-dominant in gap; finite nonzero self-inner-product | **research-open** | B1,B2 |
| **B4** | **Straight-microstrip ε_eff within 5 % of Hammerstad-Jensen — the ADR-0147 milestone** | `eps_eff_fem = (βc/ω)²` within 5 % of HJ (relax→≤15 % FDTD floor if coarse-mesh tight); low `|S11|` | **research-open (highest)** | B1,B2,B3 |
| B5 | Symbolic-factorization reuse across the sweep (banked scaling; B5b bicgstab deferred) | refactored `sweep_matrix` bit-identical (`<1e-12`) to pre-refactor on the WR-90 thru fixture | eng | — |
| B6 | Grading harness: port `oracle_reference`/`oracle_grade` to `main` + wire `ladder_s21` | reference grid passband ripple ≤0.5 dB, −3 dB at ~1.9/2.1 GHz; grader asymmetry discriminator fires | eng | — |
| B7 | 3-pole Chebyshev microstrip-filter geometry + S21 graded vs ladder (incl. asymmetry) | `|S21|(f)` within `oracle_grade` mask **AND** depth(1.6 GHz) > depth(2.4 GHz) | research-open | B4,B5,B6 |

**Critical path** = B1 → B3 → B4 → B7 (the port-physics spine, one continuous context). **Parallel
periphery** = B2 (mesh), B5 (scaling), B6 (grading). The `ladder_s21` grading **primitive** is on
`main` (`lumped.rs:202`) but the `oracle_reference`/`oracle_grade` **examples are not** (they live
only on a worktree/branch) — B6 ports them.

## Consequences

**Achievable in ~2.5-3 weeks** to a first *graded* filter S21, **high** confidence on the
engineering bricks (B1/B2/B5/B6), **moderate** on the make-or-break physics (B4). **Risk
concentrates entirely in B4**: whether the injected quasi-TEM modal *shape* (not just β, which is
robust/phase-derived) is faithful enough — the air/substrate `E_z` discontinuity and the
cross-section→port-face frame map have never fed a `yee-fem` port. **B4 GO/fork is explicit**: pass
→ proceed to filter; fail → mesh-refine vs bridge `yee-mom::NumericalCrossSection::with_quasi_tem`
vs TL-de-embed Z₀ from line currents (days-to-a-week each). **Scaling does not block the first
curve**: a single line (B4) and a small coarse 3-pole filter (B7) fit the existing per-ω faer
complex LU in a 6 g box (~40-70 k tets ceiling); `faer::matrix_free::bicgstab` is the deferred
B5b escalation (near-resonance convergence the only un-de-risked unknown).

**Resource discipline (standing for this track):** every heavy cargo invocation runs through the
bounded Docker box (`scripts/yee-box.sh`, 6 g / 2 cpu cgroup) so a build/solve can never OOM the
host; parallel fan-out is restricted to read-only agents. **Honesty gates:** the reviewer enforces
`gate_is_real` (no tautology / match-by-construction — the #1 failure mode); B4 is the real physics
check on a closed-form answer; the B7 asymmetry discriminator guards against a fitted artifact. No
EM result merges until its gate genuinely passes.

**Not in scope / do NOT reopen:** the ADR-0133 FDTD cavity wall, the ADR-0064 planar-MoM port,
the `fem-eig-006` *eigen* modal-projection wave-port (this is the *driven*-sweep track, distinct).

---

## References
- Brick-1 code: `crates/yee-fem/src/open_boundary.rs`
  (`with_interior_pec_edges` / `interior_edges_matching`),
  `crates/yee-fem/tests/open_boundary_interior_pec.rs`.
- Driven-sweep baseline: `OpenBoundarySolver::sweep_matrix` (`open_boundary.rs:1580`), validated
  WR-90 TE₁₀ thru in `crates/yee-fem/tests/open_boundary_sweep_matrix.rs`.
- Scoping workflow: `fem-em-blocker` (run `wf_5fccf1e8-72a`).
- Spec: `docs/superpowers/specs/2026-06-02-fem-em-microstrip-s21-track-design.md`;
  plan: `docs/superpowers/plans/2026-06-02-fem-em-microstrip-s21-track.md`.
