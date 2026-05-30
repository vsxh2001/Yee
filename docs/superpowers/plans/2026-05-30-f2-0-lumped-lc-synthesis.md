# Filter Phase F2.0 — lumped LC ladder synthesis — Plan

**Spec:** `2026-05-30-f2-0-lumped-lc-synthesis-design.md` · **ADR:** ADR-0111

## Lane
`crates/yee-filter/**` ONLY (new `src/lumped.rs`, `lib.rs` re-export, `tests/`).
Consume `yee-synth` (already a dep) for the prototype g-values. Out of lane →
finding. WASM-safe: pure `f64` + serde, no native dep.

## Base
New worktree off current `main` (re-fetch first — cloud-race). Branch
`feature/filter-f2-0-lumped-lc`.

## Pattern files (MIRROR)
- `crates/yee-filter/src/dimension.rs` — module-doc style, `Result`/`*Error`
  pattern, serde structs (`EdgeCoupledDimensions`), `lib.rs` re-exports, the
  cheb-N=5 fixture usage. Mirror its shape for `lumped.rs`.
- `crates/yee-filter/tests/dim_001_inversion_roundtrip.rs` — the gate idiom +
  the exact committed cheb_bpf fixture values (f0=2e9, fbw=0.10, z0=50, ripple
  0.5 dB, N=5, stopband [(2.4e9,40)]). Clone for `lumped_001`.
- `crates/yee-filter/src/lib.rs` — `FilterProject`, `prototype`/g-values access,
  `synthesize`, `SpecMask`, the ideal-response/mask-verdict logic (or see
  `crates/yee-cli/src/filter.rs` for how the mask verdict is computed — reuse the
  same comparison for the ladder S21 mask check).

## Steps
1. `src/lumped.rs`: `LcBranch`, `LcResonator`, `LumpedLadder`, `synthesize_lumped`
   (the LPF→BPF transform from the spec), `LumpedError`, and a private
   `ladder_s21(&LumpedLadder, f_hz) -> num_complex::Complex<f64>` ABCD cascade
   (source Z0 → per-resonator ABCD → load Z0; S21 from ABCD). Document all public
   items. (num-complex is already in the workspace — check yee-filter deps; if
   absent, do the 2×2 complex ABCD by hand with (re,im) tuples to avoid a new dep,
   OR add num-complex if another workspace crate already uses it.)
2. `lib.rs`: re-export the public items.
3. `tests/lumped_001.rs` per DoD 3 (len==N, each resonator tuned L·C·ω0²≈1, ladder
   |S21| meets the spec mask, values physical).

## Verify (exit 0; pure-math, fast — NO FDTD, NO container needed)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-filter --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-filter --jobs 2
```
Sub-second. Do NOT run `cargo test --workspace`/FDTD/mom-001.

## Escape hatch
Blocked > 25 min (the prototype g-values aren't exposed on `FilterProject` as
expected; the band-pass transform's series/shunt ordering disagrees with the
synth convention; the ladder S21 doesn't meet the mask → transform bug) → STOP +
surface the `FilterProject`/prototype shape + the computed element values + the
S21 at band edges. Do NOT weaken the mask gate; do NOT add heavy deps.

## Done when
DoD 1–3 pass; `lumped_001` green (ladder S21 meets the mask); `git diff --stat
<base>..HEAD` = only `crates/yee-filter/**` (+ maybe Cargo.lock if num-complex
added); WASM-safe preserved.
