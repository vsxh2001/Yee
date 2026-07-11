# FS.5a — Monte-Carlo yield analysis (plan)

**Spec:** `docs/superpowers/specs/2026-07-11-fs5a-yield-mc-design.md`

1. `crates/yee-surrogate/src/yield_mc.rs`: `splitmix64` + Box-Muller normals,
   `ToleranceSpec`, `YieldEstimate` (Wilson 95 % CI), `yield_estimate`.
   Register `pub mod yield_mc` + re-exports in `lib.rs`.
2. Unit tests in-module: RNG determinism, normal-moment sanity, Wilson CI
   endpoints (0 %, 100 % yield stay in [0, 1]).
3. `crates/yee-surrogate/tests/yield_mc.rs`: gates `yield-mc-001` (Φ(z)
   bracket), `yield-mc-002` (determinism), `surrogate-yield-001` (GP vs
   brute-force on the patch-resonance closed form).
4. ADR-0211; roadmap FS.5 row → FS.5a shipped; CI already runs
   `cargo test --workspace` so no workflow edit needed (tests are instant).
5. Verification: `cargo test -p yee-surrogate` exit 0;
   `cargo clippy --workspace --all-targets -- -D warnings` + `cargo fmt --check --all` exit 0.
