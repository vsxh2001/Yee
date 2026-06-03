# FEM numerical-eigenmode microstrip port — implementation plan

**Spec:** [2026-06-03-fem-numerical-eigenmode-port-design.md](../specs/2026-06-03-fem-numerical-eigenmode-port-design.md)
**ADR:** [ADR-0154](../../src/decisions/0154-fem-numerical-eigenmode-microstrip-port.md)
**Seed:** branch `feature/fem-port-numerical-probe` (`c102d16`) — validated de-risk probe.

## Brick N1 — `microstrip_port_numerical` production API

1. `crates/yee-fem/Cargo.toml`: move `yee-mom = { workspace = true }` from `[dev-dependencies]`
   to `[dependencies]`; drop the probe comment, add a one-line production rationale (acyclic).
2. New `crates/yee-fem/src/microstrip_port_numerical.rs`:
   - `pub struct MicrostripPortGeom { trace_w, sub_h, eps_r, box_w, box_h }` (+ doc).
   - `fn microstrip_cross_section(...) -> yee_mesh::TriMesh2D` — promote `probe_cross_section` from
     the probe test verbatim (strip-as-hole; FR-4 tag 1 below `sub_h`, air tag 0 above).
   - `pub fn microstrip_port_numerical(geom: &MicrostripPortGeom, f_hz: f64)
     -> Result<PortDefinition, FemError>` — build the cross-section at internal default density
     (nx=20, ny=12), `NumericalCrossSection::new(mesh, eps, mu).with_quasi_tem().solve(f)` once,
     `Arc`-wrap, return `PortDefinition::single_mode(analytic-HJ β, numerical modal_e_t)`. Frame
     map: `e_tangential_at(p.x, p.z)` → `Vector3::new(et[0], 0.0, et[1])`.
   - Re-export from `crates/yee-fem/src/lib.rs` (`pub use microstrip_port_numerical::*`).
3. Fast unit test (non-ignored, in the new module or `tests/`): `modal_e_t` finite + nonzero +
   E_z-dominant in the gap, decaying in air (mirror the v1 `modal_e_t_is_ez_dominant_*` test); β
   matches `yee_layout::eps_eff`.
4. Verify (boxed): `cargo test -p yee-fem --lib microstrip_port_numerical` (or the unit test name);
   `cargo clippy -p yee-fem --all-targets -- -D warnings`; `cargo fmt --check`. All exit 0.

## Brick N2 — straight-line |S21| gate

5. `crates/yee-fem/tests/microstrip_eeff.rs`: replace the probe's inline numerical-port helpers
   (`probe_cross_section`, `solve_numerical_mode`, `numerical_port`, `solve_line_numerical`) with
   calls to the new `src` API. Rename the test `fem_line_eeff_001_numerical_port` →
   `fem_line_eeff_numerical_001` and turn it into a real gate.
6. Assertions: `|S21|(L2) >= 0.6`, `|S11|(L2) <= 0.2`, `ε_eff` within 5 % of HJ. Keep `#[ignore]`'d
   + the multi-minute-solve docstring so debug skips it.
7. Verify (boxed, release): `cargo test -p yee-fem --release --test microstrip_eeff -- --ignored
   fem_line_eeff_numerical_001 --nocapture` → exit 0, prints |S21|≈0.778.

## Brick N3 — filter S21 re-grade (separate increment, after N1/N2 merge)

8. `crates/yee-fem/tests/microstrip_filter_s21.rs`: swap port construction to
   `microstrip_port_numerical`; re-run the driven sweep.
9. Honest gate: assert measured lift over the −42 dB v1 floor + asymmetry (depth(1.6) > depth(2.4));
   record the strict-mask margin (clear → assert pass; short → assert measured margin honestly).
10. Verify (boxed, release): `cargo test -p yee-fem --release --test microstrip_filter_s21 --
    --ignored fem_filter_s21_vs_ladder --nocapture` → exit 0.

## CI

- N1 unit test runs in the default debug workspace test (fast, non-ignored).
- N2/N3 release gates run in the existing `fem-eigen` `--release` gate job (`.github/workflows/
  ci.yml`); confirm the job's test list includes them. `libfontconfig1-dev` already installed there.

## Dispatch

- **N1+N2 = one agent**, one worktree off `c102d16` (the probe seed — already has the dep + code),
  lane `crates/yee-fem/src/**`, `crates/yee-fem/tests/microstrip_eeff.rs`, `crates/yee-fem/Cargo.toml`.
- N3 = a later agent after N1/N2 merge (depends on the merged `src` API).
- Reviewer (never self-review) after each; I verify boxed + merge `--no-ff`.
