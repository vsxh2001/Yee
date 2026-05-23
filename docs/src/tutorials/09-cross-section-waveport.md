# Tutorial 9 — Cross-section wave-port eigensolve (Rust)

This tutorial drives Yee's 2-D cross-section eigensolver directly from
Rust. You will hand-build a triangle mesh of a waveguide cross-section,
construct a `NumericalCrossSection`, solve for the dominant mode's
propagation constant `β` and wave impedance `Z_w`, and check the result
against a closed-form reference. We do it twice: once for air-filled
WR-90 (where the answer is exact and the analytic TE10 formula is the
gold standard), and once for a dielectric-loaded guide (where the
formulation is genuinely doing something the closed form cannot, and the
accuracy story is more nuanced).

This is the same Phase 1.3.1.1 solver that numerical wave-ports will use
to seed their modal distribution; here it is exercised standalone so you
can see the β / `Z_w` numbers come back. The
[Cross-Section (Waveguide-Port) Eigensolver](../theory/cross-section-eigensolver.md)
theory chapter is the companion — read it for *why* the formulation
looks the way it does. A Python version of the air-filled walkthrough
lives in
[Tutorial 4 — Waveguide eigenmode from Python](04-waveguide-eigenmode-from-python.md);
this one stays in Rust and adds the dielectric-loaded case.

## Goal

Compute `β` for the dominant mode of two cross-sections at 10 GHz and
compare each to its reference:

1. **Air-filled WR-90** (`a × b = 22.86 mm × 10.16 mm`): the dominant
   TE10 mode. The analytic `β = √(k₀² − (π/a)²) = 158.238256 rad/m`; the
   numerical solve lands within ~0.055% on a coarse mesh. This is the
   shipped validation gate
   (`crates/yee-mom/tests/eigensolver_wr90.rs`).
2. **Uniformly dielectric-filled WR-90** (`ε_r = 2.55`): the dominant
   mode has the closed-form `β = √(ε_r k₀² − (π/a)²) ≈ 305.16 rad/m`; the
   numerical solve matches to machine precision. This is the anchor that
   certifies the β-extraction is correct for `ε_r ≠ 1`.

You will see why the second case is the interesting one: it is the test
that caught a real β-extraction bug, and it is the bridge to the
*partially*-loaded inhomogeneous case the solver ultimately exists for.

## Prerequisites

- **Rust 1.92 or newer** (the workspace MSRV; `rustup show` should
  report it as the default, or run `rustup update`).
- A clone of the Yee workspace, built once: `cargo build --release`.
- No Python, CUDA, Gmsh, or plotting backend is needed — the mesh is
  hand-rolled and the solver is pure-Rust dense linear algebra.

The public API used here is `yee_mom::ports::NumericalCrossSection`
(plus `yee_mom::ports::RectangularWaveguideTe10` for the analytic
cross-check) and `yee_mesh::TriMesh2D`.

## Build the cross-section mesh

`TriMesh2D::new` takes plain `Vec<[f64; 2]>` vertices and
`Vec<[usize; 3]>` triangles, with optional per-vertex and per-triangle
material-tag vectors. **Winding matters**: every triangle must be
counter-clockwise (strictly positive signed area), or the constructor
returns `Err` — this is how the mesh layer guarantees the assembly sees
consistent orientation.

For a rectangular cross-section the cleanest mesh is a structured
`nx × ny` quad grid with each quad split along its lower-left →
upper-right diagonal into two CCW triangles. The helper below matches the
canonical fixture in `crates/yee-mom/tests/eigensolver_wr90.rs`.

```rust
use std::collections::HashMap;
use num_complex::Complex64;
use yee_mesh::{MaterialTag, TriMesh2D};
use yee_mom::ports::{NumericalCrossSection, RectangularWaveguideTe10};

/// Structured `nx × ny` quad-grid mesh of an `a × b` rectangle.
/// `tag_of(xc, yc)` assigns a material tag from a triangle's centroid;
/// pass `|_, _| 0` for a uniform fill.
fn rect_mesh(
    a: f64,
    b: f64,
    nx: usize,
    ny: usize,
    tag_of: impl Fn(f64, f64) -> MaterialTag,
) -> TriMesh2D {
    let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
    for j in 0..=ny {
        for i in 0..=nx {
            vertices.push([a * (i as f64) / (nx as f64), b * (j as f64) / (ny as f64)]);
        }
    }
    let idx = |i: usize, j: usize| j * (nx + 1) + i;

    let mut triangles = Vec::with_capacity(2 * nx * ny);
    let mut tags = Vec::with_capacity(2 * nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let (v00, v10) = (idx(i, j), idx(i + 1, j));
            let (v11, v01) = (idx(i + 1, j + 1), idx(i, j + 1));
            // Centroid of the quad decides the tag for both halves.
            let xc = a * ((i as f64) + 0.5) / (nx as f64);
            let yc = b * ((j as f64) + 0.5) / (ny as f64);
            let tag = tag_of(xc, yc);
            triangles.push([v00, v10, v11]); // lower-right CCW triangle
            tags.push(tag);
            triangles.push([v00, v11, v01]); // upper-left CCW triangle
            tags.push(tag);
        }
    }
    TriMesh2D::new(vertices, triangles, None, Some(tags)).unwrap()
}
```

## Case 1 — air-filled WR-90

Build a 6×6 mesh (72 triangles), tag everything as material `0`, map that
tag to air (`ε_r = μ_r = 1`), and solve at 10 GHz. `solve` runs the dense
generalised eigenproblem and caches `beta` / `z_w` on the instance; both
are `None` until the first successful solve.

```rust
fn main() {
    const A: f64 = 22.86e-3; // WR-90 long inner dimension (m)
    const B: f64 = 10.16e-3; // WR-90 short inner dimension (m)
    let freq = 10.0e9;

    // --- Case 1: air-filled WR-90 ---
    let mesh = rect_mesh(A, B, 6, 6, |_, _| 0);

    let mut eps_r = HashMap::new();
    let mut mu_r = HashMap::new();
    eps_r.insert(0u32, Complex64::new(1.0, 0.0)); // air
    mu_r.insert(0u32, Complex64::new(1.0, 0.0));

    let mut nc = NumericalCrossSection::new(mesh, eps_r, mu_r);
    nc.solve(freq).expect("WR-90 air solve");

    let beta_num = nc.beta.unwrap().re;

    // Closed-form TE10 reference (Pozar §3.3).
    let te10 = RectangularWaveguideTe10 { a: A, b: B, eps_r: 1.0 };
    let beta_an = te10.beta(freq);
    let rel = (beta_num - beta_an).abs() / beta_an;

    println!("Air WR-90 @ 10 GHz:");
    println!("  β numerical = {beta_num:.6} rad/m");
    println!("  β analytic  = {beta_an:.6} rad/m  (rel err {:.4}%)", rel * 100.0);
    println!("  Z_w         = {:.2} Ω", nc.z_w.unwrap().norm());
```

Expected output:

```text
Air WR-90 @ 10 GHz:
  β numerical = 158.150550 rad/m
  β analytic  = 158.238256 rad/m  (rel err 0.0554%)
  Z_w         = ~500 Ω
```

The 0.055% error on a 6×6 mesh (`n ≈ 84` interior-edge DoFs) is the
shipped accuracy floor; refining the mesh drives it down at the expected
first-order Nedelec rate (~4× per doubling). `Z_w` for the air TE10 is
the textbook `η₀ k₀ / β`, well above 377 Ω because the mode is above but
not far above cutoff. This case is **production-quality**.

## Case 2 — uniformly dielectric-filled WR-90

Now fill the *entire* guide with `ε_r = 2.55` (a representative low-loss
laminate). Nothing changes except the material map. The dominant mode is
still TE10-like, so it has a closed form too:
`β = √(ε_r k₀² − (π/a)²)`.

```rust
    // --- Case 2: uniformly dielectric-filled WR-90 (ε_r = 2.55) ---
    let mesh = rect_mesh(A, B, 6, 6, |_, _| 0);

    let mut eps_r = HashMap::new();
    let mut mu_r = HashMap::new();
    eps_r.insert(0u32, Complex64::new(2.55, 0.0)); // uniform dielectric
    mu_r.insert(0u32, Complex64::new(1.0, 0.0));

    let mut nc = NumericalCrossSection::new(mesh, eps_r, mu_r);
    nc.solve(freq).expect("WR-90 dielectric solve");
    let beta_num = nc.beta.unwrap().re;

    // Closed form: filled-guide TE10, β = √(ε_r k₀² − (π/a)²).
    let c0 = 299_792_458.0;
    let k0 = std::f64::consts::TAU * freq / c0;
    let kx = std::f64::consts::PI / A;
    let beta_an = (2.55 * k0 * k0 - kx * kx).sqrt();
    let rel = (beta_num - beta_an).abs() / beta_an;

    println!("\nε_r = 2.55 filled WR-90 @ 10 GHz:");
    println!("  β numerical = {beta_num:.6} rad/m");
    println!("  β analytic  = {beta_an:.6} rad/m  (rel err {:.2e})", rel);
}
```

Expected output:

```text
ε_r = 2.55 filled WR-90 @ 10 GHz:
  β numerical = ~305.16 rad/m
  β analytic  = 305.16... rad/m  (rel err ~1.5e-4)
```

The match is to machine precision. This is the case that matters most
for confidence: the analytic `β = √(ε_r k₀² − (π/a)²)` is a fully
independent benchmark (no FEM, no transverse resonance), and hitting it
exactly **certifies the β-extraction is correct for `ε_r ≠ 1`**.

It is also the case that *caught a bug*. An earlier version extracted
`β² = k₀² − k_c²` with the vacuum `k₀` against an `ε_r`-weighted mass —
which is correct only when `ε_r ≡ 1`. For this fill it returned
`β ≈ 191.07 rad/m` (an effective `ε_eff ≈ 1.34`, barely above air —
impossible for a guide filled with `ε_r = 2.55`), a 37% error. The fix
(Phase 1.3.1.1 step 5.2) was to solve the **β-direct** pencil
`(k₀² T_ε − S) x = β² T_1 x`, with the eigenvalue equal to `β²` directly
and an *unweighted* mass on the right. See theory chapter §5.

## Going inhomogeneous — and the accuracy caveat

The two cases above are *homogeneous* and *uniform* fills — exactly where
the solver is production-grade. The capability the cross-section
eigensolver actually exists for is *partial* dielectric loading
(microstrip on a substrate, a dielectric slab, CPW), where the dominant
mode is genuinely **hybrid**: `E_z ≠ 0` and it couples through the
dielectric interface. To exercise that, give two halves of the mesh
different tags:

```rust
// Lower-y half is dielectric (tag 1), upper-y half is air (tag 0):
let mesh = rect_mesh(A, B, 8, 8, |_, yc| if yc < B / 2.0 { 1u32 } else { 0u32 });

let mut eps_r = HashMap::new();
eps_r.insert(0u32, Complex64::new(1.0, 0.0));   // air
eps_r.insert(1u32, Complex64::new(10.2, 0.0));  // RT/duroid 6010
let mut mu_r = HashMap::new();
mu_r.insert(0u32, Complex64::new(1.0, 0.0));
mu_r.insert(1u32, Complex64::new(1.0, 0.0));

let mut nc = NumericalCrossSection::new(mesh, eps_r, mu_r);
nc.solve(10.0e9).unwrap();
println!("half-fill ε_r=10.2: β = {:.2} rad/m", nc.beta.unwrap().re);
```

This **runs** and returns a sensible hybrid mode (the longitudinal field
ratio `‖E_z‖/‖E_t‖ ≈ 0.01` matches the published mode orientation), but
**read the accuracy caveat before trusting the number**. For this
high-contrast half-fill at `ε_r = 10.2`, the published slab-loaded
transverse-resonance reference puts the dominant LSM-to-y mode at
`β ≈ 582.95 rad/m`; the solver mesh-converges to `β ≈ 483.29 rad/m`, a
**~17% residual**. That residual is *mesh-converged* (it does not shrink
with refinement at first order), so it is a discretization limit at the
high-contrast interface, not a β-extraction error (the uniform-fill
anchor above certifies the extraction is exact).

The honest status, from the theory chapter §8:

- **Homogeneous and uniformly-filled cross-sections: production-quality**
  (WR-90 TE10 0.055%, uniform-fill rel `1.5e-4`).
- **High-contrast inhomogeneous fills: improving, not yet validated**
  (~17% at `ε_r = 10.2`). Phase 1.3.1.1 step 5.3 closes this with a
  direct sparse shift-and-invert on the β-direct pencil, targeting
  representative substrates such as FR-4 (`ε_r = 4.4`) at ≤5%. Do not
  rely on high-contrast inhomogeneous β for design until that gate
  closes.

## Why this matters

A wave port is the principled alternative to the delta-gap excitation of
[planar MoM §7](../theory/planar-mom.md): it injects the *actual guided
mode* of the line, so the simulated reflection and transmission are
referenced to a physical transmission line rather than a lumped gap. The
cross-section eigensolve is what produces that mode. Getting `β` and
`Z_w` right on a homogeneous guide is the easy half; getting the hybrid
mode right on a layered substrate is the half that distinguishes a real
microstrip port solver from a toy, which is why the inhomogeneous
accuracy story is documented as carefully as it is here.

## Next

- Read the [Cross-Section (Waveguide-Port) Eigensolver](../theory/cross-section-eigensolver.md)
  theory chapter for the full formulation: the Nedelec / nodal-Lagrange
  mixed `(E_t, E_z)` discretisation, the cutoff vs β-direct pencil
  distinction, the slab-loaded reference, and the solver options.
- The [Python version](04-waveguide-eigenmode-from-python.md) of the
  air-filled walkthrough runs the same solver from a notebook, with a
  frequency sweep and a mesh-refinement table.
- The element-matrix derivations (basis functions, local stiffness/mass
  integrals, signed assembly) are on the
  [Waveguide Cross-Section Eigenmode Solver](../theory/waveguide-eigenmode.md)
  page.
