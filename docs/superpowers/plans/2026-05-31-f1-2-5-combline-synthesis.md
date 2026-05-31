# Filter F1.2.5 — Combline dimensional synthesis — Plan

**Spec:** `2026-05-31-f1-2-5-combline-synthesis-design.md` · **ADR:** ADR-0144

## Lane
`crates/yee-filter/**` (engine + gate). NO yee-layout/studio edits. Out of lane → finding.

## Base / worktree
New worktree off `main` (re-fetch first). Branch `feature/f1-2-5-combline`.

## Pattern files (READ FIRST)
- `crates/yee-filter/src/dimension.rs` — `dimension_hairpin` / `HairpinDimensions`
  (MIRROR exactly: the `yee_layout` imports, `microstrip_width`, `eps_eff`, `solve_gap`,
  `target_k = fbw·m[i][i+1]`, `DimError`, `C = 299_792_458`, the `#![warn(missing_docs)]`
  doc density). `dimension_edge_coupled` for the resonator-length idiom (λg/2 there;
  combline is `θ0/β`).
- `crates/yee-synth/src/lib.rs` — `prototype(Approximation::Chebyshev{ripple_db}, n)`
  → g-values; `crates/yee-filter/src/lib.rs` — `synthesize`, `CouplingMatrix`
  (`m`, `qe_in`, `qe_out`). The gate uses these for the H&L benchmark.
- `crates/yee-filter/tests/dim_stepped_001.rs` / `hairpin_dim_001` — the gate-test style.
- The spec §Method + §DoD (the H&L eq 5.46 numbers + the resonance root-find) + ADR-0144.

## Steps
1. `ComblineDimensions` + `dimension_combline(project, theta0_rad, substrate)` per the
   spec: validate `theta0_rad ∈ (0, π/2)` (else `DimError`); `line_width = microstrip_width(z0,…)`;
   `ε_eff = eps_eff(width,…)`; `β(f0) = 2π·f0·√ε_eff/c`; `resonator_length = θ0/β(f0)`;
   `loading_cap = cot(θ0)/(2π·f0·z0)`; gaps from `solve_gap(fbw·m[i][i+1])` (reuse).
   Re-export from the crate root. Document all public items.
2. Gate `crates/yee-filter/tests/dim_combline_001.rs` (spec §DoD):
   - **(1) H&L eq 5.46:** `synthesize` a 5-pole `Chebyshev{0.1}` BP at FBW=0.1; assert the
     `CouplingMatrix` `qe_in` ≈ 11.468 and the adjacent `m[i][i+1]·FBW` ≈ [0.07975 (M12),
     0.06077 (M23)] within ±1% (tighten to ±1e-3 if g-values reproduce H&L). Second point
     FBW=0.15 → Qe≈7.645, M12≈0.11962, M23≈0.09115. (Read the actual `m`/`qe` field names
     from `CouplingMatrix`; M = the normalized `m[i][i+1]` × FBW, or whatever the studio's
     `target_k` derivation uses — match it.)
   - **(2) resonance consistency:** from `dimension_combline` output build
     `B(f) = −(1/z0)·cot(θ0·f/f0) + 2π·f·C_L`, root-find `B(f)=0` (bisection/scan) over
     [0.5·f0, 1.5·f0], assert the root ≈ f0 within ±1%. Independent of the cap emit formula.
   - **(3)** gaps solved (Ok), `θ0≥π/2` → `DimError`, dims/`C_L` positive + finite.
3. (If the `CouplingMatrix` does not directly expose `qe`/`m` in the form needed, compute
   the H&L Qe/M from `project.prototype.g` + FBW directly — `Qe=g0·g1/FBW`,
   `M=FBW/√(g_i g_{i+1})` — and ALSO assert they match `synthesize`'s coupling output, so
   the gate ties the published numbers to the real synthesis.)

## Verify (run these; expected EXIT 0; quote output)
- `cargo test -p yee-filter --test dim_combline_001` — quote "test result: ok" + the
  synthesized Qe/M values + the resonance root.
- `cargo test -p yee-filter` (full crate, no regressions).
- `cargo clippy -p yee-filter --all-targets -- -D warnings` ; `cargo fmt --check`.
- `cargo check --workspace`.
(yee-filter is light/pure — host fine; NO Docker box, NO FDTD.)

Commit: `yee-filter: combline dimensional synthesis + H&L §5.2.5 gate (F1.2.5, ADR-0144)`
+ the Co-Authored-By trailer.

## Escape hatch
If the synthesized Qe/M do NOT match H&L eq 5.46 within ±1% (e.g. the crate's Chebyshev
g-values or the coupling convention differ), SURFACE it with the actual values — do NOT
loosen the tolerance to force a pass or hardcode the answer. If the resonance root is not
at f0, the cap/length formula has a bug — surface it (do not widen tol). Do NOT gate on
`C_L == cot(θ0)/(2π f0 z0)` (tautological — the research flagged this). NEVER edit
yee-synth/physics to force the benchmark. Blocked > 30 min → surface.

## Done when
`dimension_combline` + the non-vacuous `dim_combline_001` gate (H&L eq 5.46 published
Qe/M + the first-principles resonance check) are green; clippy/fmt/check clean; diff =
`crates/yee-filter/**`. Then I (dispatcher) verify + adversarial-review + merge. Studio
lighting is a follow-on.
