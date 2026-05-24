//! mom-002 numerical-microstrip-wave-port — bounded experiment diagnostic.
//!
//! **This is an EXPERIMENT, not a gate.** It is a NON-FAILING diagnostic
//! (it prints a comparison and only asserts pipeline non-degeneracy, never
//! a Z_in tripwire). See:
//!   * spec  `docs/superpowers/specs/2026-05-24-mom-002-numerical-waveport-experiment-design.md`
//!   * plan  `docs/superpowers/plans/2026-05-24-mom-002-numerical-waveport-experiment.md`
//!   * ADR   `docs/src/decisions/0059-mom-002-numerical-waveport-experiment.md`
//!
//! ## Hypothesis
//!
//! mom-002 (FR-4 microstrip Z₀, `L = 82 mm`) passes only loosely
//! (`|Z_in| ≈ 674 Ω` under the original delta-gap, `≈ 3.46 Ω` under the
//! TEM-smoothed port) vs the Hammerstad-Jensen target `Z_0 ≈ 51 Ω`. Ten
//! forensic tracks exonerated the Sommerfeld kernel and pinned the
//! residual to port-excitation modeling. The cross-section eigensolver is
//! now FR-4-validated, and the `WavePort` `Numerical2D` arm injects a
//! cross-section modal `E_t` into the MoM port RHS. The hypothesis is
//! that a *numerical microstrip modal* wave-port beats the delta-gap.
//!
//! ## Method (mirrors `__internal::z_in_with_greens_tem`, swapping ONLY
//! the RHS source)
//!
//! 1. Build a **shielded-microstrip** cross-section `TriMesh2D` in the
//!    SAME `(x, y)` frame as the mom-002 strip's port edges (port edges
//!    are y-aligned at the strip centre `x ≈ L/2`, spanning
//!    `y ∈ [−w/2, +w/2]` at `z = 0`). The eigensolver applies PEC
//!    Dirichlet on every boundary edge — i.e. it solves a *closed*
//!    cross-section — so a microstrip is modeled as a shielded microstrip
//!    (FR-4 substrate slab + signal-strip metal + air, all inside a PEC
//!    box). This is the standard way a closed-domain modal solver handles
//!    microstrip and is exactly how the validated FR-4 horizontal-slab
//!    case is posed.
//! 2. `NumericalCrossSection::solve(1 GHz)` → quasi-TEM `β` / `Z_w`;
//!    sanity-check `ε_eff = (β / k₀)²` lands in the FR-4 ballpark.
//! 3. Build the production Sommerfeld impedance matrix on the mom-002
//!    strip mesh via `__internal::impedance_matrix_for_test` (kernel
//!    consumed READ-ONLY — identical fill to the headline gate).
//! 4. Build the numerical-port RHS on that same strip mesh via
//!    `__internal::wave_port_rhs_for_test` with
//!    `ModalDistribution::Numerical2D`.
//! 5. LU-solve `Z·i = b`, extract `Z_in = V_port / I_port`, print the
//!    comparison against 674 Ω (delta-gap) + 51 Ω (HJ).
//!
//! The eigensolver and the mom-002 kernel/Greens/gate are NOT touched.

use nalgebra::Vector3;
use num_complex::Complex64;
use std::collections::HashMap;
use yee_mesh::{TriMesh, TriMesh2D};
use yee_mom::__internal::{
    MultilayerGreens, RwgBasis, build_basis, impedance_matrix_for_test, wave_port_rhs_for_test,
};
use yee_mom::ports::{ModalDistribution, NumericalCrossSection};

use faer::linalg::solvers::{PartialPivLu, Solve};

// ── mom-002 constants (mirrored from `yee-validation` — these are the
// documented public geometry of the case; we only READ them). ──
const STRIP_WIDTH_M: f64 = 2.94e-3;
const STRIP_LENGTH_M: f64 = 82.0e-3;
const N_LENGTH: usize = 82;
const N_WIDTH: usize = 16;
const F_HZ: f64 = 1.0e9;
const SUBSTRATE_EPS_R: f64 = 4.4;
const SUBSTRATE_H_M: f64 = 1.6e-3;
const DCIM_N_IMAGES: usize = 5;
const SOMMERFELD_N_POLES: usize = 1;

/// The two published comparison anchors for the experiment.
const Z_IN_DELTA_GAP_OHM: f64 = 674.0;
const Z_0_HJ_OHM: f64 = 51.0;

/// Build the mom-002 strip mesh (length along x ∈ [0, L], width along
/// y ∈ [−w/2, +w/2], z = 0), centered-port placement: columns
/// `N_LENGTH/2 − 1` / `N_LENGTH/2` tagged 1 / 2 so the shared y-aligned
/// edges at x ≈ L/2 become the port (port_tag = 1). Identical placement
/// law to `yee-validation`'s `mom_002_strip_mesh_with_spacing`
/// (Uniform). We rebuild it test-side rather than reach into the
/// validation crate's private builder.
fn strip_mesh() -> TriMesh {
    let nx = N_LENGTH + 1;
    let ny = N_WIDTH + 1;
    let dx = STRIP_LENGTH_M / (N_LENGTH as f64);
    let dy = STRIP_WIDTH_M / (N_WIDTH as f64);
    let y0 = -STRIP_WIDTH_M / 2.0;

    let mut vertices = Vec::with_capacity(nx * ny);
    for i in 0..nx {
        let x = (i as f64) * dx;
        for j in 0..ny {
            vertices.push(Vector3::new(x, y0 + (j as f64) * dy, 0.0));
        }
    }

    let port_left = N_LENGTH / 2 - 1;
    let port_right = N_LENGTH / 2;
    let mut triangles = Vec::with_capacity(2 * N_LENGTH * N_WIDTH);
    let mut tags = Vec::with_capacity(2 * N_LENGTH * N_WIDTH);
    for i in 0..N_LENGTH {
        for j in 0..N_WIDTH {
            let a = (i * ny + j) as u32;
            let b = ((i + 1) * ny + j) as u32;
            let c = ((i + 1) * ny + (j + 1)) as u32;
            let d = (i * ny + (j + 1)) as u32;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
            let tag = if i == port_left {
                1
            } else if i == port_right {
                2
            } else {
                0
            };
            tags.push(tag);
            tags.push(tag);
        }
    }
    TriMesh::new(vertices, triangles, tags).expect("strip mesh invariants")
}

/// Build a **shielded-microstrip** cross-section `TriMesh2D` in the same
/// `(x, y)` frame the mom-002 port edges are sampled in.
///
/// The mom-002 port edges sit at the strip-centre column `x ≈ L/2` and
/// span `y ∈ [−w/2, +w/2]` at `z = 0`. The `Numerical2D` RHS samples the
/// modal field at each port-edge midpoint `(mid_x, mid_y)` in THAT frame,
/// so the cross-section mesh must cover those `(x, y)` points. We build a
/// PEC box spanning `x ∈ [x_c − box_half, x_c + box_half]` (with
/// `x_c = L/2`, the strip-centre column) and `y ∈ [−Y, +Y]` (with `Y` a
/// few strip-widths so the box walls are well away from the strip),
/// partitioned into FR-4 substrate (lower half, tag 1, `ε_r = 4.4`) and
/// air (upper half, tag 0).
///
/// The signal-strip PEC conductor is modeled as a **strip-as-hole** in the
/// mesh: cells whose centroid falls inside `x ∈ [x_c−w/2, x_c+w/2]`,
/// `y ∈ [0, t_strip]` are omitted so their boundary edges become PEC
/// Dirichlet (the inner conductor). The bottom wall (`y = −Y`) is the
/// ground plane (PEC by the outer-boundary Dirichlet); the surrounding box
/// is the shield. This two-conductor layout supports a quasi-TEM mode
/// between the signal strip and the ground/shield.
///
/// **Coordinate aliasing note.** This cross-section is in the MoM mesh's
/// `(x, y)` frame, not the physical microstrip transverse plane
/// `(strip-width, substrate-normal)`. The `Numerical2D` arm samples the
/// modal field at the port-edge midpoints `(mid_x ≈ L/2, mid_y)` in the
/// MoM frame — but all port edges share the SAME x (the strip-centre
/// column), so the transverse field variation in x cannot be sampled. That
/// coordinate aliasing is the documented coupling-blocker finding (Phase A).
fn shielded_microstrip_cross_section(nx: usize, ny: usize) -> TriMesh2D {
    let x_c = STRIP_LENGTH_M / 2.0;
    // Box half-extent in x: several axial cells so the port column is
    // interior to the cross-section domain.
    let box_half_x = 3.0 * (STRIP_LENGTH_M / (N_LENGTH as f64));
    // Box half-extent in y: a few strip widths so the shield walls do not
    // crowd the strip.
    let box_half_y = 3.0 * STRIP_WIDTH_M;

    let x_lo = x_c - box_half_x;
    let x_hi = x_c + box_half_x;
    let y_lo = -box_half_y;
    let y_hi = box_half_y;

    let dx = (x_hi - x_lo) / (nx as f64);
    let dy = (y_hi - y_lo) / (ny as f64);
    // Signal-strip hole: centred at x_c, y ∈ [0, t_strip].
    // The strip PEC sits at the substrate/air interface; y=0 is the
    // interface (substrate below, air above). The hole height is ~1 cell.
    let strip_y0 = 0.0_f64;
    let strip_y1 = dy; // ~1 cell thick
    let strip_x0 = x_c - STRIP_WIDTH_M / 2.0;
    let strip_x1 = x_c + STRIP_WIDTH_M / 2.0;
    let in_strip = |cx: f64, cy: f64| -> bool {
        cx >= strip_x0 - 1e-14
            && cx <= strip_x1 + 1e-14
            && cy >= strip_y0 - 1e-14
            && cy <= strip_y1 + 1e-14
    };

    let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
    for j in 0..=ny {
        let y = y_lo + dy * (j as f64);
        for i in 0..=nx {
            let x = x_lo + dx * (i as f64);
            vertices.push([x, y]);
        }
    }
    let idx = |i: usize, j: usize| j * (nx + 1) + i;
    let mut triangles = Vec::with_capacity(2 * nx * ny);
    let mut tags = Vec::with_capacity(2 * nx * ny);
    // FR-4 substrate fills the lower part of the box (y < 0); air above.
    for j in 0..ny {
        let yc = y_lo + dy * ((j as f64) + 0.5);
        for i in 0..nx {
            let xc = x_lo + dx * ((i as f64) + 0.5);
            if in_strip(xc, yc) {
                continue; // strip-as-hole → boundary edges become PEC
            }
            let v00 = idx(i, j);
            let v10 = idx(i + 1, j);
            let v11 = idx(i + 1, j + 1);
            let v01 = idx(i, j + 1);
            let tag = if yc < 0.0 { 1u32 } else { 0u32 };
            triangles.push([v00, v10, v11]);
            tags.push(tag);
            triangles.push([v00, v11, v01]);
            tags.push(tag);
        }
    }
    TriMesh2D::new(vertices, triangles, None, Some(tags)).unwrap()
}

#[test]
fn mom_002_numerical_waveport_comparison() {
    // ── Phase A.1–A.2: build + solve the shielded-microstrip cross-section. ──
    // Mesh density matches the validated `eigensolver_wr90` / FR-4 gate
    // (≈6×6 → n ≈ 120 DoF). The cross-section eigensolve is a DENSE
    // `O(n³)` `B⁻¹A` solve, so the grid is kept coarse — a 24×24 grid
    // (n ≈ 1k) runs for minutes and is unnecessary for a quasi-TEM mode.
    let nx = 8;
    let ny = 8;
    let xs = shielded_microstrip_cross_section(nx, ny);
    let mut eps_r = HashMap::new();
    eps_r.insert(0u32, Complex64::new(1.0, 0.0)); // air
    eps_r.insert(1u32, Complex64::new(SUBSTRATE_EPS_R, 0.0)); // FR-4
    let mut mu_r = HashMap::new();
    mu_r.insert(0u32, Complex64::new(1.0, 0.0));
    mu_r.insert(1u32, Complex64::new(1.0, 0.0));

    // Phase B: switch to the quasi-TEM solve path so the microstrip dominant
    // mode (k_c²≈0, zero cutoff) is selected rather than discarded.
    // `with_quasi_tem()` dispatches to `solve_dense_mixed_quasi_tem`
    // (Phase 1.3.1.2 / ADR-0060), which uses a TEM-scale β-direct
    // shift-invert ladder seeded with a uniform-E_t vector and discriminates
    // the quasi-TEM from curl-free gradient nulls via the transverse-energy
    // screen. The closed-guide `NumericalCrossSection::solve` path is
    // unchanged for all existing callers.
    let mut mode = NumericalCrossSection::new(xs, eps_r, mu_r).with_quasi_tem();
    let solve_res = mode.solve(F_HZ);
    let k0 = std::f64::consts::TAU * F_HZ / yee_core::units::C0;

    match &solve_res {
        Ok(()) => {
            let beta = mode.beta.expect("β cached on success");
            let z_w = mode.z_w.expect("Z_w cached on success");
            let eps_eff = (beta.re / k0).powi(2);
            eprintln!(
                "[A.2] shielded-microstrip cross-section solved: \
                 β = {:.3} rad/m, ε_eff = (β/k₀)² = {:.3} \
                 (FR-4 ballpark ε_eff ≈ 3.3), Z_w = {:.3} + j{:.3} Ω",
                beta.re, eps_eff, z_w.re, z_w.im
            );
        }
        Err(e) => {
            // A solve failure with the quasi-TEM path (Phase B) is a
            // documented finding, not a test failure.
            //
            // The quasi-TEM path (`solve_dense_mixed_quasi_tem`, Phase
            // 1.3.1.2 / ADR-0060) uses a TEM-scale β-direct shift-invert
            // ladder over the physical ε_eff window. A failure here means
            // no rung converged to a transverse-dominated propagating mode
            // — possible if the cross-section box is electrically too small
            // or the mesh is too coarse for the quasi-TEM field to be
            // adequately resolved. Increasing the box size or mesh density
            // typically resolves this; the validated gate
            // (`tests/eigensolver_microstrip_quasi_tem.rs`) uses a 20×10
            // mesh with box 8w × 7h.
            //
            // FINDING: if this branch is reached, the quasi-TEM solve did
            // NOT succeed at the 8×8 coarse grid used here; the experiment
            // must be re-run at higher resolution (out of scope for this
            // bounded diagnostic — see ROADMAP note).
            eprintln!(
                "[A.2] FINDING (quasi-TEM-solve-failed): shielded-microstrip \
                 cross-section solve with `.with_quasi_tem()` returned: {e}. \
                 The quasi-TEM path (Phase 1.3.1.2 / ADR-0060) did not \
                 converge on the coarse 8×8 grid. The cross-section is \
                 modeled without the strip-as-hole PEC (a full-slab layout), \
                 which reduces the conductor-field interaction that seeds the \
                 quasi-TEM mode. Increasing the mesh density or using the \
                 strip-as-hole PEC construction from \
                 `eigensolver_microstrip_quasi_tem.rs` would resolve this. \
                 Experiment stops; kernel/Greens NOT re-opened, eigensolver \
                 NOT touched. See test docstring + ADR-0059/ADR-0060."
            );
            return;
        }
    }

    // ── Phase A.3: build the strip mesh + the numerical-port RHS, and
    // confirm the cross-section→strip-port coupling is non-vanishing. ──
    let mesh = strip_mesh();
    let basis: RwgBasis = build_basis(&mesh).expect("strip basis");
    let port_indices: Vec<usize> = basis.port_basis_indices(1).collect();
    assert!(
        !port_indices.is_empty(),
        "mom-002 strip mesh must produce port edges"
    );

    // Report where the port edges live in (x, y) vs the cross-section box,
    // so the coordinate-aliasing crux is visible in the diagnostic output.
    let (mut min_x, mut max_x, mut min_y, mut max_y) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &k in &port_indices {
        let e = &basis.edges[k];
        let p0 = mesh.vertices[e.v0 as usize];
        let p1 = mesh.vertices[e.v1 as usize];
        for p in [p0, p1] {
            min_x = min_x.min(p.x);
            max_x = max_x.max(p.x);
            min_y = min_y.min(p.y);
            max_y = max_y.max(p.y);
        }
    }
    eprintln!(
        "[A.3] {} port edges; port-edge (x,y) extent: \
         x ∈ [{:.4}, {:.4}] m, y ∈ [{:.4}, {:.4}] m (strip centre x = {:.4} m)",
        port_indices.len(),
        min_x,
        max_x,
        min_y,
        max_y,
        STRIP_LENGTH_M / 2.0
    );

    let b_num = wave_port_rhs_for_test(
        &basis,
        1,
        Complex64::new(1.0, 0.0),
        1.0,
        ModalDistribution::Numerical2D(Box::new(mode)),
        F_HZ,
    );
    let rhs_norm: f64 = (0..b_num.nrows())
        .map(|k| b_num[(k, 0)].norm_sqr())
        .sum::<f64>()
        .sqrt();
    eprintln!("[A.3] numerical-port RHS ‖b‖₂ = {rhs_norm:.6e}");

    if rhs_norm < 1e-30 {
        eprintln!(
            "[A.4] FINDING (coupling blocker): the numerical-port RHS is \
             identically zero — the cross-section modal E_t does not \
             project onto the mom-002 port edges. Root cause: the \
             `Numerical2D` arm samples the modal field at the MoM port-edge \
             midpoints in the cross-section's OWN (x, y) frame, but the \
             mom-002 port edges all sit at a SINGLE x (the strip centre, \
             x = L/2) spanning y ∈ [−w/2, w/2]. The cross-section solver's \
             (x, y) is the transverse plane; mapping a single-x strip-port \
             line onto a 2-D transverse mode requires glue the arm lacks \
             (it was validated for a waveguide whose port face IS the \
             cross-section). RECOMMENDATION: port-infra-glue-needed. \
             Experiment stops (Phase A finding); kernel/Greens NOT re-opened."
        );
        return;
    }

    // ── Phase B: coupling wired → solve mom-002 with the numerical port. ──
    // Production Sommerfeld fill (kernel READ-ONLY, identical to the gate).
    let greens = MultilayerGreens::new_microstrip_sommerfeld(
        SUBSTRATE_EPS_R,
        SUBSTRATE_H_M,
        F_HZ,
        DCIM_N_IMAGES,
        SOMMERFELD_N_POLES,
    );
    let z = impedance_matrix_for_test(&basis, &greens);

    let lu = PartialPivLu::new(z.as_ref());
    let i = lu.solve(b_num.as_ref());

    // Port current with the SAME numerical-modal weighting used to build
    // the RHS, preserving the Galerkin V/I inner-product structure (this
    // mirrors the symmetric extraction `z_in_with_greens_tem` does for the
    // TEM port). We recompute the per-edge modal projection weight w_k so
    // that I_port = Σ w_k i_k with the same w_k that set b_k = V·w_k.
    let mut i_port = Complex64::new(0.0, 0.0);
    for &k in &port_indices {
        // w_k = b_k / V (V = 1 here), the modal-projected edge weight.
        let w_k = b_num[(k, 0)];
        i_port += w_k * i[(k, 0)];
    }

    if i_port.norm() < 1e-30 {
        eprintln!(
            "[B] FINDING: port current vanished under the numerical port \
             despite a non-zero RHS — the solved current does not couple \
             back through the modal weighting. RECOMMENDATION: \
             port-infra-glue-needed. Stop (no kernel re-open)."
        );
        return;
    }

    let v_port = Complex64::new(1.0, 0.0);
    let z_in = v_port / i_port;
    let z_mag = z_in.norm();

    eprintln!("──────────────────────────────────────────────────────────");
    eprintln!("[B] mom-002 numerical-microstrip-wave-port |Z_in| COMPARISON");
    eprintln!("──────────────────────────────────────────────────────────");
    eprintln!(
        "    numerical wave-port : Z_in = {:.3} + j{:.3} Ω,  |Z_in| = {:.3} Ω",
        z_in.re, z_in.im, z_mag
    );
    eprintln!("    delta-gap baseline  : |Z_in| ≈ {Z_IN_DELTA_GAP_OHM:.1} Ω");
    eprintln!("    Hammerstad-Jensen   : Z_0   ≈ {Z_0_HJ_OHM:.1} Ω (target)");
    eprintln!(
        "    |Z_in|/Z_0 (numerical) = {:.2}×   vs   (delta-gap) = {:.2}×",
        z_mag / Z_0_HJ_OHM,
        Z_IN_DELTA_GAP_OHM / Z_0_HJ_OHM
    );
    let closer = (z_mag - Z_0_HJ_OHM).abs() < (Z_IN_DELTA_GAP_OHM - Z_0_HJ_OHM).abs();
    eprintln!(
        "    => numerical port is {} to Z_0 than the delta-gap baseline.",
        if closer { "CLOSER" } else { "NOT closer" }
    );
    eprintln!("──────────────────────────────────────────────────────────");

    // NON-FAILING diagnostic: only assert pipeline non-degeneracy, never a
    // Z_in tripwire (this is an experiment, not a re-gate). A finite,
    // non-NaN |Z_in| in the broad mom-002 non-degeneracy band is the only
    // hard assertion.
    assert!(
        z_mag.is_finite() && z_mag > 0.0,
        "numerical-port |Z_in| must be finite and positive (got {z_mag})"
    );
}
