# FEM microstrip-S21 track — Implementation Plan

**Spec:** `2026-06-02-fem-em-microstrip-s21-track-design.md` · **ADR:** ADR-0153

Driven brick-by-brick as validated increments (one or more per loop tick): each brick →
worktree off current `main` → agent on the brick's lane → **dispatcher re-verifies the gate
(boxed)** → adversarial code-reviewer (`gate_is_real`, never self-review) → fix P0/P1 →
merge `--no-ff` → cleanup. Heavy cargo **always** boxed: `YEE_BOX_DIR=$(pwd) YEE_BOX_MEM=6g
YEE_BOX_CPUS=2 scripts/yee-box.sh cargo …`.

## Bricks

### B1 — interior-PEC wiring + geometric edge-picker  ✅ SHIPPED (this ADR's merge)
- **Lane:** `crates/yee-fem/src/open_boundary.rs`, `crates/yee-fem/tests/open_boundary_interior_pec.rs`
- **Built:** `with_interior_pec_edges(impl IntoIterator<Item=usize>)` (set-union into
  `pec_global_edges`, idempotent) + `interior_edges_matching(Fn(Vector3,Vector3)->bool)`
  (rebuilds the canonical `EdgeKey`/`LOCAL_EDGES` map, returns ascending deduped global IDs).
  Reuses `assemble_complex_with_pec_edges` verbatim — no solve-path change.
- **Gate (PASSED, boxed):** `cargo test -p yee-fem --test open_boundary_interior_pec` — interior-DoF
  count drops by exactly `|E|`; idempotent re-tag. 2/2 in 0.00 s.

### B2 — layered straight-microstrip mesh   (eng · parallel w/ B5,B6)
- **Lane:** `crates/yee-fem/src/microstrip_mesh.rs`, `crates/yee-fem/tests/layered_microstrip_mesh.rs`
- `layered_microstrip_mesh(W,H,L,h,trace_w,nx,ny,nz)` → `(TetMesh3D, MaterialDatabase, ground_pred,
  trace_pred)`; `cavity_uniform` lattice with a cell boundary on `z=h` (ADR-0108 z-stack); centroid
  repaint `z<h→FR-4 tag 1, else air tag 0`; predicates feed B1's picker.
- **Gate:** tet count + substrate/air tag proportions; trace/ground pickers non-empty; `eps_at(1,ω).re==4.4`. No solve.

### B3 — quasi-TEM microstrip port closures   (research-open · spine)
- **Lane:** `crates/yee-fem/src/microstrip_port.rs`, `crates/yee-fem/tests/microstrip_port_closures.rs`
- `beta_mode(ω) = (ω/c)·√eps_eff(w,h,εr)` via `yee_layout::eps_eff`; `modal_e_t` = analytic
  E_z-dominant transverse field (parallel-plate-like in the gap, decaying into air). Inject via
  `single_mode`. Two `WavePort` faces at the line z-ends.
- **Gate:** β matches `eps_eff` to <1e-9; E_z-dominant nonzero in gap, decays above trace; finite
  nonzero modal self-inner-product.

### B4 — straight-microstrip ε_eff < 5 % HJ — THE MILESTONE / GO-fork   (research-open, highest risk · spine)
- **Lane:** `crates/yee-fem/tests/microstrip_eeff.rs`, `crates/yee-fem/src/microstrip_port.rs`
- End-to-end driven sweep (B1 trace+ground, B2 mesh, B3 port): `sweep_matrix`, extract β from port
  phase progression, `eps_eff_fem=(βc/ω)²`, compare HJ. Small mesh (~10-40 k tets) fits the 6 g box.
- **Gate (boxed `--release`):** `YEE_BOX_DIR=. scripts/yee-box.sh cargo test -p yee-fem --release
  --test microstrip_eeff -- --ignored fem_line_eeff_001 --nocapture` — within 5 % HJ (relax ≤15 %
  FDTD floor only if coarse-mesh-tight); low `|S11|`.
- **DECISION:** pass → B7; fail → fork (mesh-refine / `yee-mom::NumericalCrossSection` bridge /
  TL-de-embed). Surface the fork to the maintainer honestly — do not fake a pass.

### B5 — symbolic-factorization reuse   (eng · parallel)
- **Lane:** `crates/yee-fem/src/open_boundary.rs`
- Split `sp_lu()` (`open_boundary.rs:1602`) into `SymbolicLu::try_new` (once) + numeric per ω.
- **Gate:** refactored `sweep_matrix` bit-identical (`<1e-12`) to pre-refactor on the WR-90 thru.
- B5b (bicgstab + shifted-LU precond) **deferred** — only if B4/B7 exceed the ~40-70 k-tet ceiling.

### B6 — grading harness to main   (eng · parallel)
- **Lane:** `crates/yee-filter/examples/{oracle_reference,oracle_grade}.rs`
- Port from the worktree/branch to `main`'s `examples/`, retargeted at `main`'s `ladder_s21`
  (`lumped.rs:202`). `oracle_grade` includes the asymmetry discriminator.
- **Gate:** reference passband ripple ≤0.5 dB, −3 dB at ~1.9/2.1 GHz; grader fires on a synthetic asymmetric input.

### B7 — 3-pole filter geometry + S21 vs ladder   (research-open · spine)
- **Lane:** `crates/yee-fem/{src/microstrip_mesh.rs, tests/microstrip_filter_s21.rs}`
- Extend B2 to a coarse 3-pole coupled-resonator geometry; drive via `sweep_matrix` + B3 ports;
  de-embed reference plane; extract S21; grade vs B6 ladder.
- **Gate (boxed `--release`):** `|S21|(f)` within `oracle_grade` mask **AND** depth(1.6 GHz) >
  depth(2.4 GHz).

## Team / order
Spine B1→B3→B4→B7 = one strong physics context. Periphery B2/B5/B6 parallel. ≤4 concurrent agents.
B4 is the make-or-break — staff it best, decide GO/fork there before any filter.

## Escape hatch
If B3/B4 modal closure proves research-grade stuck after the fork options, **surface honestly to
the maintainer with the finding** (do not fake a pass, do not weaken the gate, do not reopen the
ADR-0064/0133 walls). A genuine NO-GO at B4 is a valid, documented outcome.
