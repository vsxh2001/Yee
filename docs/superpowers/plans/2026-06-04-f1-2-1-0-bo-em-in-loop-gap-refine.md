# F1.2.1.0 — single-gap surrogate-BO EM-in-loop refinement — implementation plan

**Spec:** [2026-06-04-f1-2-1-0-bo-em-in-loop-gap-refine-design.md](../specs/2026-06-04-f1-2-1-0-bo-em-in-loop-gap-refine-design.md)
**ADR:** [ADR-0157](../../src/decisions/0157-f1-2-1-0-bo-em-in-loop-gap-refine.md)
**Fork from:** `main` (`e6b4043` or later — has `coupled_resonator_k`, ADR-0155).

## Brick F1.2.1.0

1. `crates/yee-validation/Cargo.toml`: ensure `yee-surrogate` is a dep (it deps `yee-filter` +
   `yee-fem` already; add `yee-surrogate` if missing). Do NOT touch `yee-filter`'s deps.
2. New `crates/yee-validation/tests/bo_coupling_001.rs`:
   - A fixture: build a small edge-coupled Chebyshev `FilterProject` (reuse the synth/`yee-filter`
     fixtures used by `dim-001`/`coupled-001`), `dimension_edge_coupled(project, &substrate)` →
     pick one inter-resonator gap index `i`; `seed = gaps_m[i]`, `target_k = target_k[i]`.
   - `fn refine_gap_objective(x_norm) -> f64`: unscale `x_norm∈[0,1]` → gap `g` in
     `[g_lo,g_hi]=[clamp(seed·0.5),clamp(seed·1.5)]`; build a `CoupledResonatorGeom` at `gap_s=g`
     (fixture W/h/ε_r/f0); `yee_fem::coupled_resonator_k(geom, N_PTS)` → `k_fem`; return
     `(k_fem − target_k).abs()`.
   - `yee_surrogate::minimize(objective, vec![(0.0,1.0)], BoConfig{ n_initial:3, n_iters:9, seed,.. })`;
     unscale `x_best` → `refined_gap`; one final `coupled_resonator_k` at `refined_gap` for
     `k_fem_refined` (or reuse the best history eval).
   - Gate asserts (ADR §Gate): seed-off ≥10 %, strict improve, refined <8 %. Print seed/target/
     refined gaps + k_fem's + the BO history. `#[ignore]`'d + multi-minute docstring.
   - Fast unit test (non-ignored, no FEM): normalization round-trip + fixture seed/target finite.
3. CI: add a `bo-coupling-001` step to the `fem-eigen-gate` (or a sibling) `--release` job in
   `.github/workflows/ci.yml` (mirror `fem-coupling-001/002`). `libfontconfig1-dev` already there.

## Dispatch (the misfire-lesson SPLIT)

- Agent WRITES the code only: implements the test + the objective + the fast unit test, verifies it
  COMPILES (boxed `cargo clippy -p yee-validation --all-targets -- -D warnings` + `cargo fmt --check`
  + the fast unit test), commits. **Does NOT run the heavy ~57 min gate.**
- Orchestrator (me) RUNS the heavy gate via `Bash(run_in_background)` boxed
  (`cargo test -p yee-validation --release --test bo_coupling_001 -- --ignored --nocapture`), reads
  the result: seed/refined k_fem, did BO converge <8 %, the history.
- Then: verify myself → code-reviewer (gate honesty: non-circular, seed-off real, no weakened
  tolerance; the [0,1] normalization; the objective wiring) → fix P0/P1 → merge `--no-ff`.

## Verification (orchestrator)
- compile (agent): clippy `-D warnings` + fmt + fast unit → 0.
- heavy gate (me, boxed): `bo_coupling_001 --ignored --release` → 0, BO converged <8 %, strict
  improvement, seed-off ≥10 %.

## Escape hatch
- If BO cannot reach `target_k` in-bracket (the EM-achievable k range doesn't span `target_k` for the
  fixture) → surface as an honest finding (re-pick the fixture/bracket, or report the reachability
  limit); do NOT widen the 8 % tolerance or fake convergence.
- If a per-eval re-mesh makes the loop far exceed ~57 min → surface (reduce N_PTS / evals, or note
  the cost).
