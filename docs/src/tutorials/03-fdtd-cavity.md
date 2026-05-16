# Tutorial 3 — FDTD cavity resonance

Yee's third top-level solver is a 3D FDTD on the staggered grid that
gave the project its name. This tutorial steps off the planar-MoM rail
and shows you how to drive the time-domain solver directly from a tiny
Rust binary. You will build a 50x50x50 vacuum grid, inject a Gaussian
pulse on `E_z` at the centre, and run the same source through two
different outer-boundary treatments — a hard perfect electric conductor
(reflecting cavity) and a convolutional perfectly-matched layer
(open-domain absorber). Comparing what happens on each lets you *see*
why CPML matters.

## Goal

Step the FDTD walking skeleton in two configurations side by side, log
`|E_z|^2` at a probe cell, and watch the ratio of the two traces. The
PEC run will keep ringing — energy bounces off the walls and feeds back
into the probe. The CPML run will decay smoothly toward the noise
floor as the pulse leaves the domain. This is the qualitative pattern
every FDTD user needs to internalise before doing anything more
ambitious.

## Prerequisites

- Rust 1.88+ and a workspace checkout. No Gmsh, no Python, no CUDA.
- Five minutes. The grid is small; a release build runs in under a
  second on a modest laptop.

## What's in the box

The walking-skeleton FDTD lives in `crates/yee-fdtd/`. The public
surface is small:

- `YeeGrid::vacuum(nx, ny, nz, dx)` builds a uniform vacuum grid and
  sets `dt` to 0.9x the Courant limit.
- `WalkingSkeletonSolver::new(grid)` wraps a grid with **hard PEC**
  outer faces (perfect reflector).
- `WalkingSkeletonSolver::with_cpml(grid, CpmlParams::for_grid(&grid, 10))`
  wraps the same grid with a 10-cell CPML absorber on every face.
- `step()` advances one Yee leapfrog.
  `step_with_source(i, j, k, t0, sigma)` does the same plus a
  Gaussian-in-time pulse on `E_z(i, j, k)`.

No GPU kernels, no dispersive materials, no NTFF. That is intentional —
the walking skeleton is exactly the surface area downstream crates
integrate against while the high-performance kernels are written.

## Code

Drop the following into `examples/fdtd-cavity/src/main.rs` (or any
binary crate in the workspace; the dependencies are `yee-fdtd` and
`yee-core`). It builds two solvers from identical grids, steps both
for 600 iterations, and prints the ratio of `|E_z|^2` at a probe cell
near the corner.

```rust
use yee_fdtd::{
    cpml::CpmlParams, FdtdSolver, WalkingSkeletonSolver, YeeGrid,
};

fn main() {
    // 50^3 cubic-cell vacuum grid, 1 mm cells.
    let nx = 50;
    let dx = 1.0e-3;
    let grid_pec  = YeeGrid::vacuum(nx, nx, nx, dx);
    let grid_cpml = YeeGrid::vacuum(nx, nx, nx, dx);

    // Same Courant-stable dt on both grids.
    let dt = grid_pec.dt;
    let sigma = 20.0 * dt;
    let t0    = 4.0 * sigma;

    // Source at centre, probe at a corner-ish cell.
    let (si, sj, sk) = (nx / 2, nx / 2, nx / 2);
    let (pi, pj, pk) = (5, 5, 5);

    let mut pec  = WalkingSkeletonSolver::new(grid_pec);
    let mut cpml = WalkingSkeletonSolver::with_cpml(
        grid_cpml,
        CpmlParams::for_grid(pec.grid(), 10),
    );

    let n_steps = 600;
    for n in 0..n_steps {
        pec .step_with_source(si, sj, sk, t0, sigma);
        cpml.step_with_source(si, sj, sk, t0, sigma);
        if n % 50 == 0 {
            let e_pec  = pec .grid().ez[(pi, pj, pk)];
            let e_cpml = cpml.grid().ez[(pi, pj, pk)];
            let ratio  = (e_pec / e_cpml.abs().max(1e-30)).abs();
            println!(
                "step {n:4}: |Ez|_pec = {e_pec: .3e}  |Ez|_cpml = {e_cpml: .3e}  ratio = {ratio:.2e}"
            );
        }
    }
}
```

Build and run:

```bash
cargo run --release -p fdtd-cavity
```

(If you would rather not create a new crate, paste the body of `main`
into an integration test under `crates/yee-fdtd/tests/` and run
`cargo test -p yee-fdtd --release -- --nocapture` — the pattern is
identical.)

## What you should see

The pulse leaves the source cell at the centre, propagates outward at
*c*, and reaches the probe cell after roughly
`sqrt(3) * (nx / 2 - 5) * dx / c` seconds — a few dozen Yee steps. From
that point on:

- **PEC run.** The wavefront hits the outer face and reflects back.
  Energy never leaves the domain, so `|E_z|` at the probe keeps
  bouncing in a bounded oscillation. Run for a few thousand steps and
  you will resolve discrete cavity eigenmodes — the 50x50x50 box is a
  rectangular resonator.

- **CPML run.** The wavefront enters the 10-cell PML, gets attenuated
  by the polynomial-graded conductivity profile (Roden & Gedney 2000;
  `crates/yee-fdtd/src/cpml.rs`), and exits the domain with a
  reflection coefficient several orders of magnitude smaller than the
  pulse itself. `|E_z|` at the probe rises with the outgoing wave and
  then decays toward the floor.

The `crates/yee-fdtd/tests/cpml_reflection.rs` regression test gates
the CPML implementation at **>= 30 dB** of late-time reflection
reduction versus the PEC reference. The qualitative behaviour is what
FDTD theory promises: one bounces, the other absorbs.

## Tuning knobs

- **`CpmlParams::for_grid(grid, npml)`** — bump `npml` to 12 or 15 for
  more isolation at the cost of a thicker absorber.
- **Pulse width `sigma`** — keep `sigma >> dt`. `sigma = 20 * dt` is
  generous and keeps numerical dispersion well-behaved.
- **Probe placement** — keep probes outside the outer `npml` cells so
  you measure physical fields, not PML state.

## Next

The walking skeleton is enough to develop downstream tooling against.
The full theory page at [`theory/fdtd.md`](../theory/fdtd.md) covers
Yee staggering, leapfrog ordering, Courant stability, CFS-PML, and
what Phase 2.1+ adds (per-cell materials, dispersion, NTFF, conformal
cells, CUDA kernels). Read that next.
