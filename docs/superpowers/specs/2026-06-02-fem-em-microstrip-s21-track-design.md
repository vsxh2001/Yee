# FEM driven-sweep microstrip-S21 EM-verification track — Design

**ADR:** ADR-0153 · **Date:** 2026-06-02 · **Status:** Accepted

## Problem

ADR-0147 ranked the FEM freq-domain driven-sweep **#1** (blocker-free) for breaking the full-board
EM-sim wall, but stopped at a *finding*. We need a concrete, validated path from `yee-fem`'s
existing WR-90 TE₁₀ waveguide driven-sweep thru to a **real microstrip-filter S21 graded vs the
analytic `ladder_s21` reference** (3-pole Chebyshev, 0.5 dB ripple, 2 GHz, 10 % FBW), with the
geometric-asymmetry discriminator (1.6 GHz notch deeper than 2.4 GHz) as the honesty check.

## Why this is tractable (the three primitives + the basis argument)

Confirmed on `main` by an 8-agent resource-bounded scoping workflow:

1. **Interior-PEC assembly** — `assemble_complex_with_pec_edges` (`assembly.rs:442`) already
   eliminates an arbitrary global-edge set and is already on the driven path
   (`open_boundary.rs:1297`). A floating trace conductor is just an interior-edge set.
2. **Per-tet ε_r** — already threaded (`assembly.rs:484`, proven by `dispersive_solve.rs`).
3. **Port = closures** — `beta_mode`/`modal_e_t` are `Box<dyn Fn>`; a quasi-TEM port is analytic,
   **no `yee-mom` dependency**.

**ADR-0064 does not bind.** It is planar-MoM-specific (in-plane RWG basis cannot represent
substrate-normal `E_z`). `yee-fem` uses Whitney-1 edge elements on tets — full 3-D `E`, `E_z`
first-class. This is the crux that makes FEM succeed where planar MoM provably cannot.

## Architecture: 7 bricks

The decomposition, gates, risk classes, dependencies, parallelism, and the B4 GO/fork decision
point are specified in **ADR-0153** (table + consequences) and detailed step-by-step in the
companion plan. Summary of the spine: **B1 interior-PEC wiring → B3 quasi-TEM port → B4 ε_eff<5 %
HJ (the milestone, GO/fork) → B7 filter S21 vs ladder**; periphery **B2 mesh / B5 scaling / B6
grading** run parallel.

## Components & boundaries

- **`yee-fem` interior-PEC seam** (B1): `with_interior_pec_edges` + `interior_edges_matching` —
  the only API the mesh/port bricks consume; keeps `yee-fem` off `yee-mom`.
- **Layered-mesh helper** (B2, `microstrip_mesh.rs`): `cavity_uniform` + centroid ε_r repaint +
  world-coordinate trace/ground predicates feeding B1's picker.
- **Quasi-TEM port** (B3, `microstrip_port.rs`): analytic HJ β + E_z-dominant `modal_e_t`,
  injected via the existing `single_mode` port API.
- **Grading harness** (B6, `yee-filter/examples`): `oracle_reference` (builds the ladder ref) +
  `oracle_grade` (grades extracted S21 + the asymmetry discriminator).

## Testing / validation gates

Each brick's gate is machine-checkable (see ADR-0153 table). The two heavy gates (B4 ε_eff, B7
filter S21) run `--release` inside `scripts/yee-box.sh` (6 g / 2 cpu) — minutes, not hours, on a
small coarse mesh that fits the per-ω faer LU. The honesty bar: the reviewer must confirm
`gate_is_real` (no match-by-construction); B4 validates against a closed-form answer; B7's
asymmetry discriminator rejects a fitted curve.

## Resource & integration constraints

Boxed serial cargo for all heavy work; read-only-only parallel fan-out; ≤4 concurrent agents.
Lanes: `yee-fem/`, `yee-mesh/`, `yee-filter/examples/`. The dispatcher merges (agents never touch
`main`); the gate must genuinely pass before any merge; no EM result merges unvalidated.

## Out of scope

The ADR-0133 FDTD cavity wall; the ADR-0064 planar-MoM port; the `fem-eig-006` *eigen*
modal-projection port (this is the *driven* track). B5b iterative/AMG solver (deferred until a
filter exceeds the direct-LU ceiling).
