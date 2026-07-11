# FS.5b.0 — Aggressive space mapping walking skeleton (plan)

**Spec:** `docs/superpowers/specs/2026-07-11-fs5b-space-mapping-design.md`

1. `crates/yee-surrogate/src/spacemap.rs`: `ExtractConfig` + `extract`
   (Gauss–Newton, FD Jacobian, step-halving), `SpaceMapConfig` +
   `SpaceMapResult` + `space_map` (Broyden ASM). Register + re-export.
2. In-module unit tests: extraction identity, fine=coarse 1-eval
   convergence, Broyden step sanity.
3. `crates/yee-surrogate/tests/spacemap.rs`: gate `surrogate-sm-001`
   (patch two-mode warp, ASM ≤ 5 fine evals to ≤ 0.1 %, BO same budget
   ≥ 5× worse — measured numbers pinned).
4. ADR-0213; roadmap FS.5 row → FS.5b.0 shipped (FS.5b.1 = EM-fine on the
   R.4 scenario, queued).
5. Verify: `cargo test -p yee-surrogate` + clippy floor + fmt.
