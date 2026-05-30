# Filter Phase F1.2.2 — hairpin dimensional synthesis — Plan

**Spec:** `2026-05-30-f1-2-2-hairpin-dimensional-synthesis-design.md` · **ADR:** ADR-0109

## Lane
`crates/yee-filter/**` (primary). `crates/yee-layout/**` ONLY if extending
`hairpin_bpf`/`HairpinParams` to per-section gaps (option (a)) — keep that change
minimal + backward-compatible, and keep the existing hairpin geometry gate green.
Out of lane → finding.

## Base
New worktree off current `main` (re-fetch first — cloud-race lesson). Branch
`feature/filter-f1-2-2-hairpin-dim`.

## Pattern files (MIRROR)
- `crates/yee-filter/src/dimension.rs` — `EdgeCoupledDimensions`,
  `dimension_edge_coupled` (~line 202), `dimension_edge_coupled_layout` (~263),
  the gap-bisection + the `GAP_MIN_M`/`GAP_MAX_M`/`GAP_REL_TOL`/`GAP_MAX_ITERS`
  consts, `DimError`, and the `target_k = fbw·m[i][i+1]` derivation in the module
  doc. Mirror ALL of it for hairpin.
- `crates/yee-filter/tests/dim_001_inversion_roundtrip.rs` — clone to
  `hairpin_dim_001` (same cheb N=5 / FR-4 fixture + <1% round-trip).
- `crates/yee-layout/src/lib.rs` — `hairpin_bpf` / `HairpinParams` (fields: n,
  arm_length_m, line_width_m, fold_spacing_m, coupling_gap_m, tap_offset_m,
  feed_width_m, feed_length_m), `microstrip_width`, `eps_eff`,
  `coupled_microstrip`, `coupling_coefficient`.

## Steps
1. `dimension.rs`: add `HairpinDimensions` + `dimension_hairpin` +
   `dimension_hairpin_layout` per the spec. Reuse `microstrip_width`, `eps_eff`,
   the gap-bisection helper (factor it shared with edge-coupled if clean, else
   duplicate with a comment). `arm_length = c/(4·f0·√ε_eff)`.
2. Layout: reuse `hairpin_bpf`; resolve the single-`coupling_gap_m` vs per-section
   gaps per the spec (prefer a minimal backward-compatible per-section extension
   in yee-layout; else representative-gap + documented limitation — surface which).
3. `tests/hairpin_dim_001.rs` per DoD 3.

## Verify (exit 0; pure-math, fast — no FDTD, no container needed)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-filter --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-filter --jobs 2
# if yee-layout touched:
nice -n 19 cargo clippy -p yee-layout --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-layout --jobs 2
```
Pure-math — sub-second. Do NOT run `cargo test --workspace`/FDTD/mom-001. (The
bounded container `scripts/yee-box.sh` exists for heavy work but is NOT needed
here — this is closed-form.)

## Escape hatch
Blocked > 20 min (the hairpin arm-length/coupling geometry doesn't map cleanly to
`HairpinParams`, or per-section gaps force a large yee-layout change) → STOP and
surface the exact mismatch + the HairpinParams shape. Do NOT weaken the <1%
round-trip gate; do NOT change `dimension_edge_coupled`.

## Done when
DoD 1–4 pass; `git diff --stat <base>..HEAD` = `crates/yee-filter/**` (+ minimal
`crates/yee-layout/**` if extended) + the 3 committed docs; WASM-safety preserved
(no native dep on yee-filter/yee-layout default paths).
