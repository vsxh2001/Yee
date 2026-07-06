# Plan — S.8 / F1.3.0 engine-verified stepped-impedance LPF

**Spec:** `docs/superpowers/specs/2026-07-06-s8-engine-filter-verify-lpf-design.md`

1. `crates/yee-filter/Cargo.toml`: dev-deps `yee-engine` (path) + `yee-voxel` (path).
2. Gate `crates/yee-filter/tests/engine_lpf_verify.rs` (`#[ignore]`, release):
   - Synthesize: `prototype(Butterworth, 5)` →
     `dimension_stepped_impedance_layout(…, 2 GHz, 50 Ω, 120 Ω, 20 Ω, FR-4)`.
   - Reference layout: Z₀-width through line over the same port-to-port extent,
     same bbox → identical grid; same two ports (drive + passive load), same probes
     (P2 near load, P1 on the input feed).
   - Two engine jobs (S.5 materials + dt); |S21| via `transmission_db`, |S11| via
     `reflection_db`; asserts per spec; print measured vs ideal table.
3. CI: add `cargo test -p yee-filter --release --test engine_lpf_verify -- --ignored
   --nocapture` to `compute-engine-gates`.
4. Verify: fmt, clippy `-D warnings`, fast build, release gate run (~2×80 s), record
   measured numbers; iterate tolerances only with justification.
5. Ship: ADR-0185, ENGINE-STUDIO-ROADMAP S.8 row + footer, FILTER-DESIGN-ROADMAP
   F1.3.0 note, SUMMARY.md, commit + push.
