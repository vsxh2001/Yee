# FS.0a implementation plan — auto-mesh rulebook + convergence loop

**Spec:** `docs/superpowers/specs/2026-07-07-fs0-automesh-design.md`
**ADR:** ADR-0204. **Lane:** `crates/yee-engine/**` (+ this doc lane).

1. `yee_engine::automesh` module: `auto_dx` (λ/20-in-dielectric, h/3,
   min_feature/2; clamp [1 µm, 1 mm]) + `min_feature_m` (AABB widths and
   axis-aligned inter-box gaps). Unit tests: each rule binds on a scenario
   built to make it the binding constraint; clamp test.
2. `converge_two_port(layout, reference, opts, freqs, tol, max_passes)`
   riding `yee_engine::board::two_port_board_job`: per pass rescale
   n_steps, margin_cells, air_above_cells, npml to hold the physical
   window/margins/absorber constant; dx → dx/√2; stop on max per-bin
   **linear** Δ|S21| ≤ tol; report `converged` honestly.
3. Expose `TwoPortBoardOptions::npml` (default 10 — no behavioral change
   for existing callers; `for_band` sets it).
4. Gate `crates/yee-engine/tests/board_automesh.rs`
   (`engine-automesh-001`, `#[ignore]`'d, release): S.6 stub-notch layout,
   no hand dx; asserts notch ≤ 5 % of TL theory, depth ≤ −20 dB,
   `converged == true`. Picked up by the blanket yee-engine CI gates step
   (matches no skip pattern).
5. Verification: `cargo fmt --check --all`; clippy `-D warnings`;
   `cargo test -p yee-engine --lib`; the gate in release; studio workspace
   `cargo check --tests` (consumes the widened options struct).
