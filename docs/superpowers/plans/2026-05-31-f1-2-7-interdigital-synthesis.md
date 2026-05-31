# F1.2.7 — Interdigital dimensional synthesis — Plan

**Spec:** `2026-05-31-f1-2-7-interdigital-synthesis-design.md` · **ADR:** ADR-0148

## Lane
`crates/yee-filter/src/dimension.rs` (new `dimension_interdigital` + `InterdigitalDimensions`),
`crates/yee-filter/src/lib.rs` (re-export the two new public items),
`crates/yee-filter/tests/dim_interdigital_001.rs` (new gate). NO other crates; NO studio;
NO layout. Out of lane → finding, not a fix.

## Base / worktree
New worktree off `main` (re-fetch first; main 5a1de74). Branch `feature/f1-2-7-interdigital`.

## Pattern files (READ FIRST — edit ONLY in the worktree, never the main checkout)
- `crates/yee-filter/src/dimension.rs` — MIRROR `dimension_combline` (signature, the
  `solve_gap` gap loop, `target_k = FBW·m`, the `Topology::CoupledResonator` / `N<2` error
  guards) and `ComblineDimensions`/`HairpinDimensions` (the struct shape). Reuse `solve_gap`,
  `microstrip_width`, `eps_eff`, `coupling_coefficient` verbatim. The interdigital function
  takes NO θ0 (θ = π/2 fixed); `resonator_length_m = (π/2)/β(f0) = λg/4` (the hairpin
  `arm_length_m = λg/4` computation is the exact length formula to copy). NO loading cap.
- `crates/yee-filter/tests/dim_combline_001.rs` — MIRROR all three gate parts (see the spec
  DoD). The H&L Qe/M benchmark (gate 1) is copied verbatim (`spec_5pole_cheb_01db`,
  `qe_m12_m23`, the published constants). Gate 2 is the combline resonance check with θ = π/2
  and the cap term REMOVED (`B(f) = −(1/Z0)·cot((π/2)·f/f0)`). Gate 3 mirrors
  `dim_combline_001_dims_solved_and_positive` minus the loading-cap assertions, plus
  `resonator_length_m == (π/2)/β(f0)`.
- `crates/yee-filter/src/lib.rs` — how `dimension_combline` / `ComblineDimensions` are
  re-exported (add `dimension_interdigital` / `InterdigitalDimensions` the same way).

## Steps
1. `InterdigitalDimensions { line_width_m, resonator_length_m, gaps_m, target_k }` (derive
   `Debug, Clone, PartialEq, Serialize, Deserialize` like the siblings) + full doc comment
   (cite H&L §5; note θ = π/2 / λg/4 / no cap is the interdigital distinction).
2. `dimension_interdigital(project, substrate) -> Result<InterdigitalDimensions, DimError>`:
   topology + `N<2` guards; `line_width_m = microstrip_width(z0, h, εr)`; `e_eff =
   eps_eff(w, h, εr)`; `β0 = 2π f0 √e_eff / c`; `resonator_length_m = (π/2)/β0`; gaps via the
   same `solve_gap` loop as combline (`target_k[i] = FBW·m[i][i+1]`). NO cap.
3. Re-export both from `lib.rs`.
4. `tests/dim_interdigital_001.rs` — the three gate fns (spec DoD). Keep the H&L benchmark
   IDENTICAL to combline's (it is the shared, published synthesis core); gate 2 with θ = π/2
   no-cap; gate 3 with the λg/4 length equality + no-cap struct.

## Verify (run FROM THE WORKTREE; expected EXIT 0; quote output)
- `cargo test -p yee-filter --test dim_interdigital_001 -- --nocapture` → the H&L Qe/M lines
  print, all three tests pass ("test result: ok"). Quote the H&L benchmark line.
- `cargo test -p yee-filter` (no regression in the existing dim/combline/hairpin gates).
- `cargo clippy -p yee-filter --all-targets -- -D warnings` ; `cargo fmt --check -p yee-filter`.
- `cargo check --workspace`.
- `git -C /home/hadassi/Code/Yee status --porcelain crates/` EMPTY (main untouched).

Commit (in the worktree): `yee-filter: interdigital dimensional synthesis (F1.2.7, ADR-0148)`
+ the Co-Authored-By trailer.

## Escape hatch
The gate MUST be the real non-tautological one (H&L published Qe/M + the θ=π/2 λg/4 resonance
+ length=λg/4), NOT a self-consistency check against the engine's own output. Do NOT add a
loading cap (interdigital has none). Reuse `solve_gap` — do NOT write a new coupling solver.
If `solve_gap`'s bracket can't reach a `target_k` for the 5-pole FBW=0.10 fixture (it does for
combline, so it should here — same widths), surface it rather than widening the bracket
silently. NEVER edit the main checkout or another crate. Blocked > 30 min → stop + surface.

## Done when
`dimension_interdigital` + `InterdigitalDimensions` exist + are re-exported; `dim_interdigital_001`
(3 parts) passes incl. the H&L published Qe/M and the no-cap λg/4 resonance; existing gates
unregressed; clippy/fmt/check clean; diff = `crates/yee-filter/**` only. Then I (dispatcher)
verify + adversarial code-review + merge. (Layout F1.2.8 + studio lighting follow.)
