# Filter Phase F1.2.0 — closed-form edge-coupled dimensional synthesis — Plan

**Spec:** `2026-05-30-filter-f1-2-0-dimensional-synthesis-design.md` · **ADR:** ADR-0097

## Lane
`crates/yee-filter/**` ONLY (new `src/dimension.rs`, `lib.rs` re-export,
`tests/`, `Cargo.toml` to add the `yee-layout` dep). Do NOT edit `yee-layout` or
any other crate — consume `yee-layout`'s existing public API only. Keep
`yee-filter` WASM-safe (the new `yee-layout` dep is `serde`-only; no native dep,
no FDTD). Out of lane → finding, do NOT fix.

## Base
New worktree off current `main` (base SHA pinned in the brief). Branch
`feature/filter-f1-2-0-dimensional-synthesis`.

## Pattern files
- `crates/yee-filter/src/extract.rs` — house style for a self-contained
  `yee-filter` numeric module + its doc/test shape.
- `crates/yee-filter/src/lib.rs` — the `CouplingMatrix { m, qe_in, qe_out }` /
  `FilterProject { spec, prototype, coupling, topology }` / `FilterSpec` types +
  the `pub use extract::…` re-export pattern to mirror.
- `crates/yee-layout/src/lib.rs` + `src/coupled.rs` — the API to consume:
  `microstrip_width`, `eps_eff`, `Substrate`, `EdgeCoupledParams`,
  `edge_coupled_bpf`, `coupled_microstrip`, `coupling_coefficient`. READ
  `EdgeCoupledParams`/`EdgeCoupledSection` fields before mapping into them.
- `crates/yee-synth/src/lib.rs` ~line 218 — confirm `k_{i,i+1}=FBW/√(g_i g_{i+1})`
  so the `target_k = FBW·m_{i,i+1}` equality is correct (cross-check in a comment).

## Steps
1. `Cargo.toml`: add `yee-layout = { workspace = true }` to `yee-filter`
   `[dependencies]`. Confirm no dependency cycle (`yee-layout` must NOT depend on
   `yee-filter` — it does not at base).
2. `src/dimension.rs`: `EdgeCoupledDimensions`, `DimError`,
   `dimension_edge_coupled`, `dimension_edge_coupled_layout` per the spec.
   Bisection helper for the monotonic gap inversion (relative tol ≤ 1e-4, capped
   iterations, `GapNotBracketed` if `target_k` is unreachable in the bracket).
3. `src/lib.rs`: `pub mod dimension;` + `pub use dimension::{…};`. Document every
   public item (`missing_docs = warn`).
4. `tests/dim_001_inversion_roundtrip.rs`, `tests/dim_002_sanity.rs`,
   `tests/dim_003_layout_serde.rs` per DoD 4–6.

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-filter --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-filter --jobs 2
```
Pure math — sub-second. Do NOT run `cargo test --workspace`, FDTD, mom-001.

## Escape hatch
Blocked > 15 min — the `coupled_microstrip` coupling range cannot bracket the
fixture's `target_k` for any sane gap (the N=5 0.5 dB / FBW 0.10 couplings fall
outside what FR-4 W/h≈1 edge-coupled lines can realize) → STOP and surface: the
computed `target_k` values, the achievable `coupling_coefficient` range over the
gap bracket at the chosen width, and propose either a different substrate/width
or a wider bracket. Do NOT silently clamp gaps or weaken `dim-001`'s < 1 %. Do
NOT invent a `qe`→gap formula (that is F1.2.1) — leave I/O feed as the documented
placeholder. Do NOT edit `yee-layout`.

## Done when
DoD 1–6 pass; `git diff --stat <base>..HEAD` shows only `crates/yee-filter/**`
(+ `Cargo.lock`) + the 3 committed docs; `yee-filter` has no FDTD/native dep
(still WASM-safe); `yee-layout` untouched.
