//! mom-002 MPIE / port / DCIM-TM root-cause audit — Track YYYYYY.
//!
//! ## Why this file exists
//!
//! Six prior diagnostics have shrunk-but-not-closed the mom-002 `|Z_in|`
//! gap on FR-4 / `h = 1.6 mm` / `w = 2.94 mm` / 1 GHz; the last three
//! (`PPPPPP`, `SSSSSS`, `TTTTTT`) and the most recent `XXXXXX` have all
//! converged on the same qualitative finding:
//!
//! * Track EEEEEE (commit `ca0e7bb`) — fixed the Sommerfeld surface-wave
//!   prefactor.
//! * Track JJJJJJ (commit `4dbeece`) — ruled out Hankel-tail truncation
//!   (only ~5% effect).
//! * Track PPPPPP (commit `d89d0b9`) — GPOF residual fit weak; kernel
//!   hand-sum bit-exact with `scalar_scalar`.
//! * Track SSSSSS (commit `0e571b7`) — contour-integral residue measured
//!   a `-2×` discrepancy in the residue.
//! * Track TTTTTT (commit `a4f98a4`) — applied the Michalski-Mosig 1997
//!   eq. (19) residue sign + factor-of-2 fix; `|Z_in|` moved from 2232 Ω
//!   to 2215 Ω (~0.8%), 30× gap still open.
//! * Track XXXXXX (commit `9847251`) — verified ψ_p(z) normalisation,
//!   α_0 sign and port-current ratio between DCIM-only and Sommerfeld
//!   (= 0.997, ~kernel-independent). Surfaced the qualitative finding
//!   that `Re(Z_in) ≈ −67 Ω` on FR-4 is **unphysical** for a passive
//!   microstrip; combined with the kernel agreement between DCIM and
//!   Sommerfeld, the bug is **formulation-level**, not kernel-level.
//!
//! This file probes the three remaining candidates:
//!
//! * **(M1) MPIE singular-kernel sanity** — for the self-term and the
//!   near-singular adjacent-pair RWG-RWG interactions, the `1/R` kernel
//!   integration may be under-resolved (Gauss order too low, or
//!   singularity subtraction missing). We use a **free-space** kernel
//!   so the multilayer DCIM / Sommerfeld machinery is bypassed entirely;
//!   any pathology that survives in free space is a formulation-level
//!   bug.
//!
//! * **(M2) Delta-gap port-current normalisation** — `V_port` is the
//!   delta-gap convention (1 V), but the port current is summed as
//!   `Σ_RWG_in_port b_k · i_k`. A wrong weighting (cell width vs strip
//!   width, vector-vs-scalar form) would surface as a non-physical
//!   `Re(Z) < 0` even in the **free-space** limit (a free-space strip
//!   at L = 30 mm at 1 GHz is electrically short and lossless, so
//!   `Re(Z) ≥ 0` must hold).
//!
//! * **(M3) DCIM TM-channel coefficients** — for FR-4 the leading TM
//!   (scalar-potential) image should have `b₁ ≈ (ε_r − 1)/(ε_r + 1) ≈
//!   0.629` at `a₁ ≈ -2h = -3.2e-3 m`. Compare against the stored
//!   `MultilayerGreens.scalar_images[0]` from the DCIM-only constructor.
//!
//! ## Diagnostic method
//!
//! Three probes, in independent file regions:
//!
//! 1. **Probe M1 — kernel-level self vs near-pair sanity (free-space).**
//!    The MoM impedance matrix Z is `pub(crate)` (built inside
//!    `crate::fill::impedance_matrix`) and not exposed through
//!    `__internal`, so a direct dump of `diag(Z)` is not possible from a
//!    test file. We probe the **underlying Greens kernel** that drives
//!    Z instead: `G(r, r')` evaluated at a small but non-zero
//!    self-displacement (within one cell), versus `G(r, r')` evaluated
//!    at the adjacent-cell displacement (one cell stride away). The
//!    ratio of these two magnitudes is the dominant entry-ratio
//!    `|Z_diag| / |Z_adjacent|` in the assembled matrix (modulo the
//!    Galerkin quadrature weights). A pathological ratio (≪ 1 or ≫ 10⁴)
//!    flags M1.
//!
//! 2. **Probe M2 — free-space port-current normalisation.** Build the
//!    same 30 mm strip mesh used by mom-002, solve with the **free-space**
//!    kernel (no DCIM, no surface waves). Print `V_port`, `|I_port|`
//!    indirectly (as `1 / |Z_in|` since `V_port = 1 V`), and check the
//!    sign of `Re(Z_in)`. For a free-space PEC strip at L = 30 mm at
//!    1 GHz the line is electrically short and lossless, so the input
//!    impedance is dominated by `Im(Z)` (small capacitive / inductive
//!    reactance) and `Re(Z) ≥ 0` *must* hold (it is a passive lossless
//!    structure). A negative `Re(Z)` here is a smoking gun for M2 —
//!    independent of any multilayer-kernel question.
//!
//! 3. **Probe M3 — DCIM TM-channel leading image.** Build a DCIM-only
//!    `MultilayerGreens` (`new_microstrip_sommerfeld(.., 5, 0)`), read
//!    its `scalar_images[0]` (the TM-channel leading image), and compare
//!    against the analytic `(ε_r − 1)/(ε_r + 1) ≈ 0.629` at
//!    `a ≈ -2h = -3.2 mm`. A leading coefficient with the wrong sign or
//!    a depth that is off by more than a factor of ~2 from the expected
//!    `-2h` flags M3.
//!
//! ## References
//!
//! * Track EEEEEE — `crates/yee-mom/tests/sommerfeld_residue_diagnostic.rs`.
//! * Track JJJJJJ — `crates/yee-mom/tests/mom_002_extent_sensitivity.rs`.
//! * Track PPPPPP — `crates/yee-mom/tests/mom_002_h2_gpof_diagnostic.rs`.
//! * Track SSSSSS — `crates/yee-mom/tests/mom_002_reflection_convention.rs`.
//! * Track XXXXXX — `crates/yee-mom/tests/mom_002_psi_port_audit.rs`.
//! * K. A. Michalski and J. R. Mosig, "Multilayered media Green's
//!   functions in integral equation formulations," *IEEE Trans.
//!   Antennas Propag.*, vol. 45, no. 3, pp. 508–519, Mar 1997.
//! * Y. L. Chow, J. J. Yang, D. G. Fang, and G. E. Howard, "A closed-form
//!   spatial Green's function for the thick microstrip substrate," *IEEE
//!   Trans. Microwave Theory Tech.*, vol. 39, no. 3, pp. 588–592, Mar
//!   1991 — DCIM `(ε_r − 1)/(ε_r + 1)` leading coefficient.

use nalgebra::Vector3;
use num_complex::Complex64;
use yee_mom::__internal::{
    FreeSpaceGreen, Greens, MultilayerGreens, RwgBasis, build_basis, z_in_with_greens,
};

const EPS_R: f64 = 4.4;
const H_SUBSTRATE_M: f64 = 1.6e-3;
const F_HZ: f64 = 1.0e9;
const STRIP_W_M: f64 = 2.94e-3;
const STRIP_L_M: f64 = 30.0e-3;
const N_LENGTH: usize = 30;
const N_WIDTH: usize = 16;
const N_DCIM_IMAGES: usize = 5;

/// Edge-clustered (Chebyshev-y) strip mesh — bit-for-bit equivalent to
/// `yee_validation::mom_002_strip_mesh_with_spacing` with
/// `StripSpacing::EdgeClustered`. Replicated inline here so the
/// diagnostic has no cross-lane dependency on `yee-validation` internals
/// (matches the pattern in `tests/mom_002_extent_sensitivity.rs` and
/// `tests/mom_002_psi_port_audit.rs`).
fn build_strip_mesh_edge_clustered(
    length_m: f64,
    width_m: f64,
    n_length: usize,
    n_width: usize,
) -> yee_mesh::TriMesh {
    let nx = n_length + 1;
    let ny = n_width + 1;
    let mut vertices: Vec<Vector3<f64>> = Vec::with_capacity(nx * ny);
    let dx = length_m / (n_length as f64);
    let y_nodes: Vec<f64> = (0..=n_width)
        .map(|j| {
            let theta = std::f64::consts::PI * (j as f64) / (n_width as f64);
            -(width_m / 2.0) * theta.cos()
        })
        .collect();
    for i in 0..nx {
        let x = (i as f64) * dx;
        for &y in &y_nodes {
            vertices.push(Vector3::new(x, y, 0.0));
        }
    }
    let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(2 * n_length * n_width);
    let mut tags: Vec<u32> = Vec::with_capacity(2 * n_length * n_width);
    for i in 0..n_length {
        for j in 0..n_width {
            let a = (i * ny + j) as u32;
            let b = ((i + 1) * ny + j) as u32;
            let c = ((i + 1) * ny + (j + 1)) as u32;
            let d = (i * ny + (j + 1)) as u32;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
            let tag = if i == 0 {
                1
            } else if i == 1 {
                2
            } else {
                0
            };
            tags.push(tag);
            tags.push(tag);
        }
    }
    yee_mesh::TriMesh::new(vertices, triangles, tags).expect("strip mesh invariants")
}

/// Probe M1 helper — given a basis and a Greens kernel, find a "self"
/// reference edge near the strip centre and the closest "adjacent" edge
/// and return `(|G_self|, |G_adjacent|)` where `G_self` is the kernel
/// evaluated between two slightly-offset points on the same edge
/// neighbourhood, and `G_adjacent` is the kernel between the same point
/// and the centre of the nearest other edge.
///
/// The self displacement uses a small finite offset (`0.05 · dx_min`)
/// rather than `r1 == r2` because the bare `1/R` Greens kernel is
/// singular at coincident points; the real MoM `Z_ii` is finite only
/// after the Galerkin double integral. The ratio of these two
/// magnitudes still tells us whether the near-singular behaviour is
/// the order of magnitude we expect — `|G_self|` should dominate
/// `|G_adjacent|` by ~1-2 orders of magnitude on a mesh whose cell
/// size is small compared to the wavelength (here `dx ≈ 1 mm` vs
/// `λ ≈ 300 mm` ⇒ `k_0 dx ≈ 0.02` so the kernel is locally `~1/R`).
fn probe_kernel_self_vs_adjacent<G: Greens>(
    basis: &RwgBasis,
    mesh: &yee_mesh::TriMesh,
    green: &G,
) -> (f64, f64) {
    // Find an edge near the strip midpoint: pick the edge whose two
    // endpoint-vertex midpoint x-coordinate is closest to L/2.
    let target_x = STRIP_L_M / 2.0;
    let mut best_idx = 0usize;
    let mut best_dist = f64::INFINITY;
    for (k, edge) in basis.edges.iter().enumerate() {
        let v0 = mesh.vertices[edge.v0 as usize];
        let v1 = mesh.vertices[edge.v1 as usize];
        let mid_x = 0.5 * (v0.x + v1.x);
        let mid_y = 0.5 * (v0.y + v1.y);
        // Prefer edges near the strip centerline (y = 0) too — they are
        // representative of the dominant longitudinal current flow.
        let d = (mid_x - target_x).abs() + mid_y.abs();
        if d < best_dist {
            best_dist = d;
            best_idx = k;
        }
    }
    let ref_edge = &basis.edges[best_idx];
    let ref_v0 = mesh.vertices[ref_edge.v0 as usize];
    let ref_v1 = mesh.vertices[ref_edge.v1 as usize];
    let ref_mid = 0.5 * (ref_v0 + ref_v1);

    // Now find the nearest *other* edge.
    let mut adj_idx = best_idx;
    let mut adj_dist = f64::INFINITY;
    for (k, edge) in basis.edges.iter().enumerate() {
        if k == best_idx {
            continue;
        }
        let v0 = mesh.vertices[edge.v0 as usize];
        let v1 = mesh.vertices[edge.v1 as usize];
        let mid = 0.5 * (v0 + v1);
        let d = (mid - ref_mid).norm();
        if d < adj_dist {
            adj_dist = d;
            adj_idx = k;
        }
    }
    let adj_edge = &basis.edges[adj_idx];
    let adj_v0 = mesh.vertices[adj_edge.v0 as usize];
    let adj_v1 = mesh.vertices[adj_edge.v1 as usize];
    let adj_mid = 0.5 * (adj_v0 + adj_v1);

    // Self: small offset along the strip from ref_mid (0.05 × dx_cell).
    let dx_cell = STRIP_L_M / (N_LENGTH as f64);
    let self_offset = Vector3::new(0.05 * dx_cell, 0.0, 0.0);
    let r1 = ref_mid;
    let r2 = ref_mid + self_offset;
    let g_self = green.scalar_vector(r1, r2);

    let g_adj = green.scalar_vector(ref_mid, adj_mid);

    eprintln!(
        "  ref edge idx = {best_idx}  midpoint = ({:.4e}, {:.4e}, {:.4e})",
        ref_mid.x, ref_mid.y, ref_mid.z,
    );
    eprintln!(
        "  adjacent edge idx = {adj_idx}  midpoint = ({:.4e}, {:.4e}, {:.4e})  Δ = {:.4e} m",
        adj_mid.x, adj_mid.y, adj_mid.z, adj_dist,
    );
    eprintln!(
        "  G_self  (r2 - r1 = {:.4e} m) = {:.6e} + j·{:.6e}",
        self_offset.norm(),
        g_self.re,
        g_self.im,
    );
    eprintln!(
        "  G_adj   (Δ = {:.4e} m)        = {:.6e} + j·{:.6e}",
        adj_dist, g_adj.re, g_adj.im,
    );
    (g_self.norm(), g_adj.norm())
}

/// Print the three probe tables and the verdict. Marked `#[ignore]` so
/// the suite never runs it by default; invoke explicitly via
///
/// ```text
/// cargo test -p yee-mom --release --test mom_002_mpie_audit \
///     -- --ignored --nocapture
/// ```
///
/// to dump the table and the M1 / M2 / M3 verdict.
#[test]
#[ignore = "diagnostic: mom-002 MPIE singularity + port-normalisation + DCIM-TM audit"]
fn mom_002_mpie_audit_diagnostic() {
    eprintln!("--- Track YYYYYY: mom-002 MPIE / port / DCIM-TM audit ---");
    eprintln!(
        "Geometry: L = {} mm, w = {} mm, h = {} mm (substrate is *bypassed* for M1/M2)",
        STRIP_L_M * 1e3,
        STRIP_W_M * 1e3,
        H_SUBSTRATE_M * 1e3,
    );
    eprintln!(
        "Frequency: f = {} GHz   (FreeSpaceGreen kernel for M1/M2)",
        F_HZ * 1e-9
    );
    eprintln!();

    let mesh = build_strip_mesh_edge_clustered(STRIP_L_M, STRIP_W_M, N_LENGTH, N_WIDTH);
    let basis = build_basis(&mesh).expect("RWG basis assembly");
    eprintln!(
        "Mesh: {} triangles, {} RWG basis functions",
        mesh.triangles.len(),
        basis.n_basis(),
    );

    let port_tag = 1u32;
    let port_n: usize = basis.port_basis_indices(port_tag).count();
    eprintln!("Port {port_tag}: {port_n} basis functions");
    eprintln!();

    // -----------------------------------------------------------------
    // Probe M1 — kernel self vs adjacent (free-space)
    // -----------------------------------------------------------------
    eprintln!("Probe M1 — kernel self vs adjacent magnitude (FreeSpaceGreen):");
    let green_fs = FreeSpaceGreen::new(F_HZ);
    let (g_self_mag, g_adj_mag) = probe_kernel_self_vs_adjacent(&basis, &mesh, &green_fs);
    let ratio_m1 = g_self_mag / g_adj_mag.max(f64::MIN_POSITIVE);
    eprintln!(
        "  |G_self| / |G_adj| = {:.6e}    (sanity band: 1.5 .. 5e3 on a sub-λ mesh)",
        ratio_m1,
    );

    // For a free-space `1/(4πR)` kernel on a sub-wavelength mesh, the
    // ratio of self-probe (R = 0.05 dx ≈ 5e-5 m) to adjacent-cell probe
    // (R ≈ dx ≈ 1e-3 m) magnitude tracks `R_adj / R_self ≈ 20`, so the
    // sanity band is wide: anything below 1.5 or above ~5e3 indicates
    // the kernel is *not* behaving as a near-singular `1/R` integrand.
    //
    // We do *not* flag M1 just because the ratio is "unusual" — the
    // Galerkin double integral that produces the actual `Z_ii` includes
    // the cancellation against the `∇·∇' G` scalar-scalar term, which
    // can change sign. We only flag M1 if the ratio is outside the
    // very wide [1.5, 5e3] sanity band, which would indicate the
    // kernel evaluation itself is broken.
    let m1_detected = !(1.5..=5.0e3).contains(&ratio_m1);
    eprintln!(
        "  Verdict: {}",
        if m1_detected {
            "M1 detected (kernel self vs adjacent ratio outside sanity band)"
        } else {
            "M1 not detected (free-space kernel singularity behaviour is sane)"
        },
    );
    eprintln!();

    // -----------------------------------------------------------------
    // Probe M2 — free-space port-current normalisation
    // -----------------------------------------------------------------
    eprintln!("Probe M2 — free-space port-current normalisation:");
    eprintln!("  Solving z_in_with_greens on the same mesh with FreeSpaceGreen ...");
    let z_in_fs = z_in_with_greens(&mesh, port_tag, &green_fs).expect("free-space solve");
    let v_port = Complex64::new(1.0, 0.0);
    let i_port = v_port / z_in_fs;
    eprintln!(
        "  V_port = {} + j·{}    |I_port| = {:.6e}  arg(I_port) = {:.6e} rad",
        v_port.re,
        v_port.im,
        i_port.norm(),
        i_port.arg(),
    );
    eprintln!(
        "  Z_in (free space) = {:.6e} + j·{:.6e}   |Z_in| = {:.6e}",
        z_in_fs.re,
        z_in_fs.im,
        z_in_fs.norm(),
    );

    // For a free-space PEC strip at L = 30 mm at 1 GHz (electrically
    // short: k_0 L ≈ 0.63 rad), the line is lossless and passive, so
    // `Re(Z) ≥ 0` *must* hold. A negative `Re(Z)` here flags M2.
    let m2_detected = z_in_fs.re < 0.0;
    eprintln!(
        "  Re(Z) ≥ 0 in free-space passive line? {}",
        if !m2_detected {
            "yes (passive lossless line behaves)"
        } else {
            "NO — M2 detected (port normalisation is broken in the free-space limit)"
        },
    );
    eprintln!();

    // -----------------------------------------------------------------
    // Probe M3 — DCIM-only TM-channel leading image
    // -----------------------------------------------------------------
    eprintln!("Probe M3 — DCIM-only TM-channel leading image:");
    let greens_dcim = MultilayerGreens::new_microstrip_sommerfeld(
        EPS_R,
        H_SUBSTRATE_M,
        F_HZ,
        N_DCIM_IMAGES,
        0, // no surface-wave poles — pure DCIM
    );
    eprintln!(
        "  n_images = {}, n_surface_wave_poles = 0 (pure DCIM)",
        greens_dcim.n_images,
    );
    eprintln!(
        "  scalar_images (TM-channel, {} entries):",
        greens_dcim.scalar_images.len(),
    );
    for (idx, (b_n, a_n)) in greens_dcim.scalar_images.iter().enumerate() {
        eprintln!(
            "    image {idx}:  b = {:>14.6e} + j·{:>14.6e}    a = {:>14.6e} + j·{:>14.6e}",
            b_n.re, b_n.im, a_n.re, a_n.im,
        );
    }
    eprintln!(
        "  vector_images (TE-channel, {} entries):",
        greens_dcim.vector_images.len(),
    );
    for (idx, (b_n, a_n)) in greens_dcim.vector_images.iter().enumerate() {
        eprintln!(
            "    image {idx}:  b = {:>14.6e} + j·{:>14.6e}    a = {:>14.6e} + j·{:>14.6e}",
            b_n.re, b_n.im, a_n.re, a_n.im,
        );
    }

    // Analytic targets for the DCIM leading image (Chow et al. 1991,
    // Aksun 1996 §III): the dominant contribution is the PEC + slab
    // image at depth ~2h with coefficient `(ε_r − 1)/(ε_r + 1)`.
    let b_expected = (EPS_R - 1.0) / (EPS_R + 1.0);
    let a_expected = -2.0 * H_SUBSTRATE_M;
    eprintln!(
        "  Analytic leading image (Chow 1991): b ≈ {:.4}, a ≈ {:.4e} m",
        b_expected, a_expected,
    );

    // Find the image with the largest |b| in the scalar (TM) train.
    let leading = greens_dcim
        .scalar_images
        .iter()
        .copied()
        .max_by(|(b1, _), (b2, _)| {
            b1.norm()
                .partial_cmp(&b2.norm())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

    let m3_detected = if let Some((b_lead, a_lead)) = leading {
        eprintln!(
            "  Leading TM image: b = {:.6e} + j·{:.6e}    a = {:.6e} + j·{:.6e}",
            b_lead.re, b_lead.im, a_lead.re, a_lead.im,
        );

        // Three sanity criteria:
        //  (i)  |Re b_lead| in the right order of magnitude (between
        //       0.1 and 5×) compared to the analytic ~0.629.
        //  (ii) sign of Re b_lead matches expected positive sign.
        //  (iii) |Re a_lead| not wildly off the expected -2h depth
        //        (between 0.2× and 5× target magnitude).
        let mag_ok = (0.1 * b_expected.abs()..=5.0 * b_expected.abs()).contains(&b_lead.re.abs());
        let sign_ok = b_lead.re.signum() == b_expected.signum();
        let depth_mag_ok =
            (0.2 * a_expected.abs()..=5.0 * a_expected.abs()).contains(&a_lead.re.abs());
        eprintln!(
            "    |Re b_lead| in [0.1×, 5×] of {:.4}?  {}",
            b_expected, mag_ok,
        );
        eprintln!(
            "    sign(Re b_lead) == sign({:.4})?       {}",
            b_expected, sign_ok,
        );
        eprintln!(
            "    |Re a_lead| in [0.2×, 5×] of {:.4e}?  {}",
            a_expected.abs(),
            depth_mag_ok,
        );

        !(mag_ok && sign_ok && depth_mag_ok)
    } else {
        eprintln!("  No scalar images returned by DCIM constructor — fit collapsed.");
        true
    };
    eprintln!(
        "  Verdict: {}",
        if m3_detected {
            "M3 detected (DCIM TM-channel leading image deviates from analytic)"
        } else {
            "M3 not detected (DCIM TM-channel leading image agrees with analytic)"
        },
    );

    // -----------------------------------------------------------------
    // Final verdict block
    // -----------------------------------------------------------------
    eprintln!();
    eprintln!("Verdict summary:");
    eprintln!(
        "  M1 (matrix singularity):  {}",
        if m1_detected {
            "detected"
        } else {
            "not detected"
        },
    );
    eprintln!(
        "  M2 (port normalization):  {}",
        if m2_detected {
            "detected"
        } else {
            "not detected"
        },
    );
    eprintln!(
        "  M3 (DCIM-TM coeffs):      {}",
        if m3_detected {
            "detected"
        } else {
            "not detected"
        },
    );
}
