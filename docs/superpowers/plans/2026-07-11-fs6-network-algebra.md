# FS.6.0 — Network algebra walking skeleton (plan)

**Spec:** `docs/superpowers/specs/2026-07-11-fs6-network-algebra-design.md`

1. `Error::Network(String)` variant in `yee-io/src/lib.rs`.
2. `crates/yee-io/src/network.rs`: `s_to_t`, `t_to_s`, `cascade`,
   `deembed_left`, `cascade_files`; module docs derive the T convention.
   Register `pub mod network` + re-exports.
3. `crates/yee-io/tests/network_algebra.rs`: gate `net-001` (7 cases per
   the spec).
4. ADR-0212; roadmap FS.6 row → FS.6.0 shipped. CI: covered by the
   workspace test job (instant, non-ignored).
5. Verify: `cargo test -p yee-io` exit 0; clippy floor + fmt check exit 0.
