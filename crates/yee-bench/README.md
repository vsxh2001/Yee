# yee-bench

Criterion benchmarks for Yee's hot paths.

## Usage

```sh
cargo bench -p yee-bench
```

Run a single bench binary:

```sh
cargo bench -p yee-bench --bench mom_solve
cargo bench -p yee-bench --bench fdtd_step
cargo bench -p yee-bench --bench gmres_vs_direct
```

Compile (without running) all benches:

```sh
cargo bench -p yee-bench --no-run
```

## Benches

| Binary            | What it measures                                              |
|-------------------|---------------------------------------------------------------|
| `mom_solve`       | Full `PlanarMoM::run` on a thin-cylinder dipole (`n_axial = 8`, `n_around = 8`, 128 triangles), single-frequency sweep near the half-wave resonance. |
| `fdtd_step`       | One `update_h + update_e` pair on a 50³ vacuum grid via `WalkingSkeletonSolver::step`. |
| `gmres_vs_direct` | `PartialPivLu` vs `gmres_jacobi` on a 128×128 diagonally-dominant complex system. |

## Expected ballparks

These are not regression gates — they're rough sanity checks for a typical
dev laptop (one core, release profile, no GPU). If you see numbers
dramatically off from these, suspect a build configuration issue
(e.g. accidentally running a debug build, `RUSTFLAGS` overriding the
workspace `[profile.release]` settings) before you suspect a real
regression.

- `mom_solve_dipole_8x8_single_freq` — full sweep completes in **well under
  30 s** total (Criterion's default ~100 samples included).
- `fdtd_step_50cubed_vacuum` — **~50 ms per call** for a single Yee step on
  125 000 cells.
- `direct_lu_128` vs `gmres_jacobi_128` — both within **one order of
  magnitude** of each other. The direct LU is normally the faster of the two
  at this size; GMRES wins asymptotically once `n` is large enough that the
  `O(n³)` factorisation dominates.

## When to add a bench

Add a new bench binary when you are about to optimise something. The bench
that pre-dates the optimisation is the only credible artifact that the
optimisation actually moved the needle. Conversely, do not add benches
"just in case" — every bench is wall-time on every `cargo bench` run, and
the criterion baseline data has to be re-recorded after any meaningful
change to the surrounding code.

## ML benches

- `gp_fit` — GP fit + fit_ml across n = 10, 25, 50 (+ 100 for plain fit).
- `bo_step` — full 25-eval BO minimize on the 1-D deceptive objective.

## FDTD specialty benches

- `plane_wave_step` — TF/SF 1-D incident-grid advance (one step).
- `lumped_resistor` — LumpedRlcPort `correct_e` at a single cell.
