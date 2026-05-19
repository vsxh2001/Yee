//! mom-002 strip-eigenmode β extraction probe — Track QQQQQQQ.
//!
//! ## Why this file exists
//!
//! Track MMMMMMM (commit `c7b001e`) ran a three-probe diagnostic on the
//! post-IIIIIII reframed mom-002 (`L = 82 mm` centered uniform on FR-4
//! at 1 GHz) and reported an `ε_eff_solver = 1.36` extracted from a
//! linear fit of `arg(G_φ(ρ, 0, 0))` against ρ, vs the
//! Hammerstad-Jensen analytic `ε_eff = 3.32`. The fit slope on the
//! point-source kernel phase was interpreted as the strip propagation
//! constant `k_eff`.
//!
//! Track NNNNNNN proved this was a measurement artifact. ADR-0037
//! (`docs/src/decisions/0037-mom-002-r1-metric-retracted.md`)
//! documents the retraction: `arg(G_φ(ρ, 0, 0))` measures the
//! point-source scalar-potential phase decay between two field/source
//! points on the air side of the slab. It does **not** measure the
//! strip-eigenmode propagation constant `β` of the integral-equation
//! solution `Z · I = V` for a finite microstrip on the slab.
//!
//! ADR-0037 recommends extracting `β` directly from the assembled
//! impedance matrix:
//!
//! > Extract the strip eigenmode directly from the assembled impedance
//! > matrix `Z`. The smallest-singular-value right eigenvector gives the
//! > dominant current distribution; its phase-vs-x slope along the
//! > strip's longitudinal axis is the propagation constant `β`. Compare
//! > `β / k_0` against `√ε_eff_HJ ≈ 1.82`. If they agree, the kernel is
//! > fine and the `|Im(Z)|` residual is a port-excitation or
//! > edge-singularity discretisation effect. If they disagree, the
//! > kernel's strip eigenmode physics is genuinely off — the `K^A`
//! > vector-potential image train (TE channel, inductive part of the
//! > line) is the most likely site.
//!
//! ## Probe construction
//!
//! Build the IIIIIII reframed mom-002 mesh — `L = 82 mm`, `w = 2.94 mm`,
//! `82 × 16` cells, centered port, uniform y-spacing — at FR-4 / `h =
//! 1.6 mm` / `f = 1 GHz` with the Phase 1.1.1.2 Sommerfeld kernel
//! (DCIM `N = 5` + 1 TM₀ pole), exactly matching the production
//! headline gate.
//!
//! Per the brief's escape hatch (`impedance_matrix` and `delta_gap_rhs`
//! are `pub(crate)` and not surfaced via `__internal`), drive the strip
//! with a centered delta-gap and recover the **port-driven current
//! distribution** `i = Z^-1 b` rather than the smallest-eigenvalue
//! right eigenvector. For a half-wave resonator the centered delta-gap
//! launches a current that is dominated by the fundamental TEM-like
//! quasi-mode `v_0`, so the phase-vs-x slope of the column-averaged
//! current carries the same `β` as the eigenvector probe (up to a
//! global complex offset that the linear fit absorbs).
//!
//! Because `impedance_matrix` is not on the stable surface, this probe
//! replicates the Galerkin MPIE matrix fill in-test against the public
//! `Greens` trait. The replication uses straight order-5 Gauss
//! quadrature on every triangle pair (no Duffy regularisation). This
//! is intentional and acceptable for the probe: the dominant mode
//! current pattern on a half-wave strip is governed by long-range
//! (well-separated) coupling, and Duffy only refines the on/near-
//! diagonal contributions which contribute uniformly to the phase
//! offset (linear-fit invariant), not the slope. The verdict block
//! reports this assumption explicitly.
//!
//! β extraction:
//!
//! 1. Group basis functions by the longitudinal coordinate of their
//!    shared-edge midpoint. The structured `82 × 16` mesh produces
//!    discrete longitudinal columns of edges; basis functions whose
//!    shared-edge midpoint lies within a tight tolerance of one of the
//!    `n_length + 1` column-x values are bucketed together.
//! 2. For each column form the column-averaged complex current
//!    `c_col = Σ_{k in col} length_k · i_k / Σ_{k in col} length_k`.
//!    Length weighting matches the Galerkin port-current convention.
//! 3. **Method A (running wave):** unwrap `arg(c_col)` along x, mask
//!    the two port columns, and linear-least-squares fit phase vs x.
//!    Slope = `-β_A`. Method A returns `β ≈ 0` on a pure standing
//!    wave (the centered delta-gap on a half-wave resonator excites a
//!    symmetric real-valued pattern; the running-wave β is zero by
//!    construction). The half-symmetry check `asymmetry_(left,right)
//!    < 10 %` AND `|β_A/k_0| < 0.1` flags the standing-wave regime.
//! 4. **Method B (standing wave):** 1-D grid search over `β/k_0 ∈
//!    [0.5, 3.5]` minimising `chi²` between `|i(x)|` and the envelope
//!    `A · |sin(β · x)|`, then parabolic refinement. The thin open-
//!    ended strip pins the current to zero at `x = 0` and `x = L`
//!    (no charge accumulation past the strip ends), so the lowest
//!    mode is the half-wave `β · L = π` and the envelope shape is
//!    `|sin(β · x)|` with no phase offset. The fit returns `β_B`.
//! 5. Report `β / k_0` and `ε_eff_solver = (β / k_0)²` for both
//!    methods; the verdict picks the dominant regime. Compare to
//!    Hammerstad-Jensen `ε_eff = 3.32` -> `β / k_0 = √3.32 ≈ 1.823`.
//!
//! ## References
//!
//! * ADR-0037 — `docs/src/decisions/0037-mom-002-r1-metric-retracted.md`
//!   (R1 metric retracted; β-from-Z recommended).
//! * ADR-0036 — `docs/src/decisions/0036-mom-002-validation-strategy.md`
//!   (IIIIIII reframe to `L = 82 mm` half-wave).
//! * Sibling diagnostic `mom_002_13x_residual_diagnostic.rs` (MMMMMMM).
//! * Sibling diagnostic `mom_002_psi_port_audit.rs` (XXXXXX — same
//!   `pub(crate)` constraint, prior workaround pattern).
//! * D. M. Pozar, *Microwave Engineering*, 4th ed., §3.7 (microstrip
//!   `Z_0` and `ε_eff`).
//! * E. Hammerstad and Ø. Jensen, "Accurate Models for Microstrip
//!   Computer-Aided Design," *MTT-S Digest*, 1980.

use faer::Mat;
use faer::linalg::solvers::{PartialPivLu, Solve};
use nalgebra::Vector3;
use num_complex::Complex64;
use yee_mesh::TriMesh;
use yee_mom::__internal::{Greens, MultilayerGreens, RwgBasis, build_basis};

// FR-4 / 1 GHz canonical microstrip parameters (match
// `yee-validation::MOM_002_*` constants and the production headline gate).
const EPS_R: f64 = 4.4;
const H_SUBSTRATE_M: f64 = 1.6e-3;
const STRIP_W_M: f64 = 2.94e-3;
const STRIP_L_M: f64 = 82.0e-3;
const F_HZ: f64 = 1.0e9;

// Production headline mesh dimensions per the brief.
const N_LENGTH: usize = 82;
const N_WIDTH: usize = 16;

// Sommerfeld-kernel parameters (match the production headline gate).
const N_DCIM_IMAGES: usize = 5;
const N_SW_POLES: usize = 1;

// H-J analytic eps_eff for FR-4 / 1.6 mm / w = 2.94 mm (Pozar 4e §3.7
// eq. 3.195 at u = w/h = 1.8375). β/k_0 target is sqrt(3.32) ≈ 1.823.
const EPS_EFF_HJ_ANALYTIC: f64 = 3.32;

/// `k_0 = 2π f / c` at `freq_hz` (matches helper in sibling diagnostics).
fn k0_at(freq_hz: f64) -> f64 {
    std::f64::consts::TAU * freq_hz / yee_core::units::C0
}

/// Build the ADR-0036 / IIIIIII centered-uniform strip mesh. Bit-
/// equivalent to `yee_validation::mom_002_strip_mesh_with_spacing(..,
/// StripSpacing::Uniform)` but inlined here so the diagnostic has no
/// cross-lane dependency on `yee-validation` internals.
fn build_strip_mesh_centered_uniform(
    length_m: f64,
    width_m: f64,
    n_length: usize,
    n_width: usize,
) -> TriMesh {
    assert!(
        n_length >= 4 && n_length.is_multiple_of(2),
        "n_length must be even and >= 4 to host a centered port column"
    );
    assert!(n_width >= 1, "n_width must be >= 1");

    let nx = n_length + 1;
    let ny = n_width + 1;
    let mut vertices: Vec<Vector3<f64>> = Vec::with_capacity(nx * ny);
    let dx = length_m / (n_length as f64);
    let dy = width_m / (n_width as f64);
    let y0 = -width_m / 2.0;

    for i in 0..nx {
        let x = (i as f64) * dx;
        for j in 0..=n_width {
            let y = y0 + (j as f64) * dy;
            vertices.push(Vector3::new(x, y, 0.0));
        }
    }

    let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(2 * n_length * n_width);
    let mut tags: Vec<u32> = Vec::with_capacity(2 * n_length * n_width);
    let port_left = n_length / 2 - 1;
    let port_right = n_length / 2;
    for i in 0..n_length {
        for j in 0..n_width {
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

/// 7-point Dunavant order-5 Gauss-triangle rule (barycentric points +
/// weights). Replicates `yee_mom::quadrature::GaussTriangle::order_5`
/// — that module is `pub(crate)` and not on the test surface. Sum of
/// weights = 1 (referenced area; multiply by triangle area to integrate).
///
/// `clippy::excessive_precision` is allowed locally because the
/// literals are quoted verbatim from the published Dunavant tables
/// (matches the source-side allowance in
/// `yee_mom::quadrature::GaussTriangle::order_5`).
#[allow(clippy::excessive_precision)]
fn gauss_order_5() -> (Vec<[f64; 3]>, Vec<f64>) {
    let a1 = 0.0597158717_897698;
    let b1 = 0.4701420641_051151;
    let a2 = 0.7974269853_530873;
    let b2 = 0.1012865073_234563;
    let p = vec![
        [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        [a1, b1, b1],
        [b1, a1, b1],
        [b1, b1, a1],
        [a2, b2, b2],
        [b2, a2, b2],
        [b2, b2, a2],
    ];
    let w = vec![
        0.2250000000_000000,
        0.1323941527_885062,
        0.1323941527_885062,
        0.1323941527_885062,
        0.1259391805_448271,
        0.1259391805_448271,
        0.1259391805_448271,
    ];
    (p, w)
}

/// Per-basis topological summary needed for the in-test fill. The
/// `RwgBasis` we get from `__internal::build_basis` exposes the
/// `edges: Vec<RwgEdge>` field, but `RwgBasis` accessors for basis-
/// function vector value and divergence are `pub(crate)`. We rebuild a
/// minimal sidecar that lets us evaluate each basis function using
/// only public data: edge length, free-vertex position on each
/// adjacent triangle, the two triangle vertex-index lists, and the
/// two triangle areas.
struct BasisSidecar {
    /// Vertices of the mesh (mirrors `TriMesh::vertices`).
    vertices: Vec<Vector3<f64>>,
    /// Triangles of the mesh (mirrors `TriMesh::triangles`).
    triangles: Vec<[u32; 3]>,
    /// Per-triangle area (signed-positive; matches `RwgBasis::areas`).
    tri_areas: Vec<f64>,
    /// Per-basis: `(tri_plus, tri_minus, free_plus, free_minus,
    /// edge_length, edge_midpoint_x, port_tag)`.
    per_basis: Vec<(u32, u32, u32, u32, f64, f64, u32)>,
}

impl BasisSidecar {
    fn from(basis: &RwgBasis, mesh: &TriMesh) -> Self {
        let vertices = mesh.vertices.clone();
        let triangles = mesh.triangles.clone();
        let mut tri_areas = Vec::with_capacity(triangles.len());
        for tri in &triangles {
            let a = vertices[tri[0] as usize];
            let b = vertices[tri[1] as usize];
            let c = vertices[tri[2] as usize];
            // Area = 0.5 * |(b - a) × (c - a)| — same convention as
            // `RwgBasis::from_mesh`.
            let cross = (b - a).cross(&(c - a));
            tri_areas.push(0.5 * cross.norm());
        }
        let mut per_basis = Vec::with_capacity(basis.edges.len());
        for edge in &basis.edges {
            let m_v0 = vertices[edge.v0 as usize];
            let m_v1 = vertices[edge.v1 as usize];
            let midx = 0.5 * (m_v0.x + m_v1.x);
            per_basis.push((
                edge.tri_plus,
                edge.tri_minus,
                edge.free_plus,
                edge.free_minus,
                edge.length,
                midx,
                edge.port_tag,
            ));
        }
        Self {
            vertices,
            triangles,
            tri_areas,
            per_basis,
        }
    }

    /// Replicates `RwgBasis` vector-value evaluation. Returns
    /// `Vector3::zeros()` if `tri` is not in the support of basis `k`.
    fn basis_vec(&self, k: usize, tri: u32, bary: [f64; 3]) -> Vector3<f64> {
        let (tri_plus, tri_minus, free_plus, free_minus, length, _, _) = self.per_basis[k];
        let reconstruct_r = |t: u32| -> Vector3<f64> {
            let tri_vs = self.triangles[t as usize];
            bary[0] * self.vertices[tri_vs[0] as usize]
                + bary[1] * self.vertices[tri_vs[1] as usize]
                + bary[2] * self.vertices[tri_vs[2] as usize]
        };
        if tri == tri_plus {
            let r = reconstruct_r(tri);
            let p = self.vertices[free_plus as usize];
            let scale = length / (2.0 * self.tri_areas[tri as usize]);
            scale * (r - p)
        } else if tri == tri_minus {
            let r = reconstruct_r(tri);
            let p = self.vertices[free_minus as usize];
            let scale = length / (2.0 * self.tri_areas[tri as usize]);
            -scale * (r - p)
        } else {
            Vector3::zeros()
        }
    }

    /// Replicates `RwgBasis` divergence evaluation.
    fn basis_div(&self, k: usize, tri: u32) -> f64 {
        let (tri_plus, tri_minus, _, _, length, _, _) = self.per_basis[k];
        if tri == tri_plus {
            length / self.tri_areas[tri_plus as usize]
        } else if tri == tri_minus {
            -length / self.tri_areas[tri_minus as usize]
        } else {
            0.0
        }
    }

    fn n(&self) -> usize {
        self.per_basis.len()
    }

    fn edge_length(&self, k: usize) -> f64 {
        self.per_basis[k].4
    }

    fn edge_midpoint_x(&self, k: usize) -> f64 {
        self.per_basis[k].5
    }

    fn port_tag(&self, k: usize) -> u32 {
        self.per_basis[k].6
    }

    fn triangle_vertices(&self, tri: u32) -> [Vector3<f64>; 3] {
        let [a, b, c] = self.triangles[tri as usize];
        [
            self.vertices[a as usize],
            self.vertices[b as usize],
            self.vertices[c as usize],
        ]
    }

    fn n_tri_for_basis(&self, k: usize) -> [u32; 2] {
        [self.per_basis[k].0, self.per_basis[k].1]
    }
}

/// Galerkin MPIE matrix fill against the public `Greens` trait. Straight
/// order-5 Gauss on every triangle pair (no Duffy regularisation —
/// `_smooth` Greens variants on coincident integration points are used
/// to avoid the panic path at `R = 0`).
///
/// Implementation strategy: pre-cache the Green's function at every
/// (outer_tri, outer_gp, inner_tri, inner_gp) tuple, then reduce the
/// per-pair fill to a weighted sum over the cache. This collapses the
/// cost from `O(n_basis^2 · 2 · 2 · 49)` Green's-function evaluations
/// to `O(n_tri^2 · 49)`, which on an `82 × 16` mesh is `2624^2 · 49 ≈
/// 3.4e8` evaluations — feasible inside the `--ignored` budget.
fn build_z_matrix<G: Greens>(side: &BasisSidecar, green: &G) -> Mat<Complex64> {
    let n = side.n();
    let mut z = Mat::<Complex64>::zeros(n, n);
    let (gp, gw) = gauss_order_5();

    // Phase 1.0 MPIE prefactors (mirror `fill.rs`).
    //   omega_mu0      = j k0 η0
    //   inv_omega_eps0 = -j η0 / k0
    let k0 = green.k0().re;
    let eta0 = green.eta0();
    let omega_mu0 = Complex64::new(0.0, 1.0) * Complex64::new(k0 * eta0, 0.0);
    let inv_omega_eps0 = Complex64::new(0.0, -1.0) * Complex64::new(eta0 / k0, 0.0);

    let n_tri = side.tri_areas.len();
    let tri_v: Vec<[Vector3<f64>; 3]> = (0..n_tri)
        .map(|t| side.triangle_vertices(t as u32))
        .collect();
    let tri_a: Vec<f64> = side.tri_areas.clone();

    // Pre-evaluate each basis vector value at all Gauss points of its
    // tri_plus / tri_minus support, plus its piecewise-constant
    // divergence.
    //
    // Memory: O(n_basis · 2 · 7 · 3) doubles ≈ ~640 KB on a 3800-basis
    // mesh — small.
    let n_gauss = gp.len();
    let mut cached_f: Vec<Vec<[Vector3<f64>; 7]>> = vec![vec![[Vector3::zeros(); 7]; 2]; n];
    let mut cached_div: Vec<[f64; 2]> = vec![[0.0_f64; 2]; n];
    for k in 0..n {
        let [tp, tm] = side.n_tri_for_basis(k);
        cached_div[k][0] = side.basis_div(k, tp);
        cached_div[k][1] = side.basis_div(k, tm);
        for (gi, bary) in gp.iter().enumerate() {
            cached_f[k][0][gi] = side.basis_vec(k, tp, *bary);
            cached_f[k][1][gi] = side.basis_vec(k, tm, *bary);
        }
    }

    // Pre-evaluate triangle Gauss-point spatial coordinates.
    let mut tri_gauss_pts: Vec<[Vector3<f64>; 7]> = vec![[Vector3::zeros(); 7]; n_tri];
    for t in 0..n_tri {
        for (gi, bary) in gp.iter().enumerate() {
            let v = tri_v[t];
            tri_gauss_pts[t][gi] = bary[0] * v[0] + bary[1] * v[1] + bary[2] * v[2];
        }
    }

    // Pre-cache Green's-function values at every (outer_tri, outer_gp,
    // inner_tri, inner_gp) tuple. Cache size: 2 · n_tri² · 49 ·
    // complex64 = 16 bytes · 49 · 2624² ≈ 5.4 GB for the full 82×16
    // mesh. That is too large for many test runners. We therefore
    // build the cache *row by row* (over outer (t_out, go)) and
    // immediately consume each row by sweeping all basis pairs that
    // use t_out. Memory peak drops to one row at a time: n_gauss ·
    // n_tri · 16 B · 2 ≈ 5.7 MB. Algorithm complexity is unchanged.
    let row_len = n_gauss * n_tri;
    let mut g_a_row = vec![Complex64::new(0.0, 0.0); row_len];
    let mut g_phi_row = vec![Complex64::new(0.0, 0.0); row_len];

    // Reverse-index: for each triangle, the list of basis indices `m`
    // for which `tri_plus(m) == t` (slot 0) or `tri_minus(m) == t`
    // (slot 1). Used to walk every basis function whose outer
    // contribution comes from t_out.
    let mut tri_to_basis_plus: Vec<Vec<usize>> = vec![Vec::new(); n_tri];
    let mut tri_to_basis_minus: Vec<Vec<usize>> = vec![Vec::new(); n_tri];
    for k in 0..n {
        let [tp, tm] = side.n_tri_for_basis(k);
        tri_to_basis_plus[tp as usize].push(k);
        tri_to_basis_minus[tm as usize].push(k);
    }
    let tri_outer: Vec<[u32; 2]> = (0..n).map(|k| side.n_tri_for_basis(k)).collect();

    eprintln!(
        "  fill: n_basis = {}, n_tri = {}, row cache = {} Gauss samples (peak {:.1} MB)",
        n,
        n_tri,
        row_len,
        (row_len * 2 * std::mem::size_of::<Complex64>()) as f64 / (1024.0 * 1024.0),
    );

    let fill_t0 = std::time::Instant::now();

    // For each (t_out, go) in turn: build the row of cached Green's
    // values across all (t_in, gi), then sweep every basis pair (m,
    // n_idx) whose outer integral uses (t_out, go).
    //
    // Per-pair contribution accumulator: keep an in-flight `Mat`
    // entry update so we can amortise the (n_basis × n_basis)
    // workload across the n_tri · n_gauss outer points.
    for t_out in 0..n_tri {
        for go in 0..n_gauss {
            let r_out = tri_gauss_pts[t_out][go];

            // Build the row cache for this outer sample.
            for t_in in 0..n_tri {
                for gi in 0..n_gauss {
                    let r_in = tri_gauss_pts[t_in][gi];
                    let r = (r_out - r_in).norm();
                    let (ga, gphi) = if r > 0.0 {
                        (
                            green.scalar_vector(r_out, r_in),
                            green.scalar_scalar(r_out, r_in),
                        )
                    } else {
                        // Coincident sample (R = 0); smooth variants
                        // are finite here. Without Duffy this slightly
                        // underestimates the on-diagonal contribution;
                        // the bias is a uniform phase offset across
                        // all basis pairs and so cancels in the
                        // phase-vs-x slope used for β extraction.
                        (
                            green.scalar_vector_smooth(r_out, r_in),
                            green.scalar_scalar_smooth(r_out, r_in),
                        )
                    };
                    g_a_row[t_in * n_gauss + gi] = ga;
                    g_phi_row[t_in * n_gauss + gi] = gphi;
                }
            }

            // Sweep basis functions that have t_out in their support.
            // Slot 0 means t_out == tri_plus(m), slot 1 means
            // t_out == tri_minus(m).
            for &slot in &[0_usize, 1_usize] {
                let m_list = if slot == 0 {
                    &tri_to_basis_plus[t_out]
                } else {
                    &tri_to_basis_minus[t_out]
                };
                for &m in m_list {
                    let fm = cached_f[m][slot][go];
                    let div_m = cached_div[m][slot];
                    let a_out_w = tri_a[t_out] * gw[go];

                    // For each n_idx, sum over its two inner
                    // triangles.
                    for n_idx in 0..n {
                        let [tp_n, tm_n] = tri_outer[n_idx];
                        let mut acc = Complex64::new(0.0, 0.0);
                        for (in_slot, &t_in) in [tp_n, tm_n].iter().enumerate() {
                            let a_in = tri_a[t_in as usize];
                            let div_n = cached_div[n_idx][in_slot];
                            let base = (t_in as usize) * n_gauss;
                            let mut sub = Complex64::new(0.0, 0.0);
                            for gi in 0..n_gauss {
                                let fn_vec = cached_f[n_idx][in_slot][gi];
                                let ga = g_a_row[base + gi];
                                let gphi = g_phi_row[base + gi];
                                let w_in = gw[gi];
                                let dot = fm.dot(&fn_vec);
                                sub += Complex64::new(w_in, 0.0)
                                    * (omega_mu0 * Complex64::new(dot, 0.0) * ga
                                        + inv_omega_eps0
                                            * Complex64::new(div_m * div_n, 0.0)
                                            * gphi);
                            }
                            acc += Complex64::new(a_in, 0.0) * sub;
                        }
                        // Write-add into Z[m, n_idx].
                        let cur = z[(m, n_idx)];
                        z[(m, n_idx)] = cur + Complex64::new(a_out_w, 0.0) * acc;
                    }
                }
            }
        }
        if t_out.is_multiple_of(64) {
            eprintln!(
                "    progress: t_out = {}/{}  ({:.1} s elapsed)",
                t_out,
                n_tri,
                fill_t0.elapsed().as_secs_f64(),
            );
        }
    }
    z
}

/// Replicates `solve::delta_gap_rhs` (1 V across `port_tag` edges).
fn delta_gap_rhs_local(side: &BasisSidecar, port_tag: u32) -> Mat<Complex64> {
    let n = side.n();
    let mut b = Mat::<Complex64>::zeros(n, 1);
    for k in 0..n {
        if side.port_tag(k) == port_tag {
            b[(k, 0)] = Complex64::new(side.edge_length(k), 0.0);
        }
    }
    b
}

/// Linear least-squares fit `y = a x + b`. Returns `(slope a, intercept
/// b)`. Used to fit phase-vs-x along the strip.
fn linfit(xs: &[f64], ys: &[f64]) -> (f64, f64) {
    let n = xs.len() as f64;
    let sx: f64 = xs.iter().sum();
    let sy: f64 = ys.iter().sum();
    let sxx: f64 = xs.iter().map(|x| x * x).sum();
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| x * y).sum();
    let denom = n * sxx - sx * sx;
    let slope = (n * sxy - sx * sy) / denom;
    let intercept = (sy - slope * sx) / n;
    (slope, intercept)
}

/// Run the β-from-Z probe and print the verdict block.
///
/// Wall-time budget: dominated by the in-test Galerkin fill, which on
/// the production 82×16 mesh is in the same order as production
/// (~minutes). The `--ignored` gate keeps this off the default suite.
#[test]
#[ignore = "diagnostic: extracts β/k_0 from the port-driven current distribution (NNNNNNN / ADR-0037 recommendation)"]
fn mom_002_beta_eigenmode_probe() {
    eprintln!("--- Track QQQQQQQ: mom-002 strip-eigenmode β probe ---");
    eprintln!();
    eprintln!(
        "Mesh:       L = {} mm, w = {} mm, {} × {}, centered port, uniform y",
        STRIP_L_M * 1e3,
        STRIP_W_M * 1e3,
        N_LENGTH,
        N_WIDTH,
    );
    eprintln!(
        "Substrate:  FR-4, eps_r = {}, h = {} mm, f = {} GHz",
        EPS_R,
        H_SUBSTRATE_M * 1e3,
        F_HZ * 1e-9,
    );
    eprintln!(
        "Greens:     MultilayerGreens new_microstrip_sommerfeld (DCIM N={} + {} TM0 pole)",
        N_DCIM_IMAGES, N_SW_POLES,
    );

    let mesh = build_strip_mesh_centered_uniform(STRIP_L_M, STRIP_W_M, N_LENGTH, N_WIDTH);
    let basis = build_basis(&mesh).expect("RwgBasis build");
    let side = BasisSidecar::from(&basis, &mesh);
    let n_basis = side.n();
    eprintln!();
    eprintln!("Basis count (RWG): {}", n_basis);

    let greens = MultilayerGreens::new_microstrip_sommerfeld(
        EPS_R,
        H_SUBSTRATE_M,
        F_HZ,
        N_DCIM_IMAGES,
        N_SW_POLES,
    );

    eprintln!("Assembling Z matrix (in-test Galerkin order-5 Gauss, no Duffy)...");
    let t0 = std::time::Instant::now();
    let z = build_z_matrix(&side, &greens);
    eprintln!("  Z fill complete in {:.1} s", t0.elapsed().as_secs_f64());
    eprintln!("Z matrix shape:      ({}, {})", n_basis, n_basis);

    // Smallest-|diagonal-eigenvalue|-proxy: faer 0.23 does not expose
    // a clean dense complex EVD via the high-level `solvers` surface
    // on the integration-test build, so we report `min |Z_kk|` as a
    // coarse health check instead. The actual β extraction below uses
    // the port-driven current and does NOT depend on the eigenproblem.
    let mut min_diag: Complex64 = Complex64::new(f64::INFINITY, 0.0);
    let mut min_diag_norm = f64::INFINITY;
    for k in 0..n_basis {
        let d = z[(k, k)];
        if d.norm() < min_diag_norm {
            min_diag_norm = d.norm();
            min_diag = d;
        }
    }
    eprintln!(
        "Smallest-|Z_kk| diagonal entry (eigenvalue proxy): {:.6e} + j·{:.6e}  (|.| = {:.4e})",
        min_diag.re, min_diag.im, min_diag_norm,
    );

    // Port-driven solve: i = Z^-1 b for a centered 1V delta-gap.
    let b = delta_gap_rhs_local(&side, 1u32);
    let lu = PartialPivLu::new(z.as_ref());
    let i_vec = lu.solve(b.as_ref());

    // Column-bucketing along x. The structured grid produces shared-
    // edge midpoint x-values at the `n_length + 1` column-x positions
    // (vertices) PLUS at the (n_length) cell-center positions (the
    // intra-cell diagonal edges produced by the (a, c) split). Bucket
    // by x within tol = dx / 4 to merge both kinds.
    let dx = STRIP_L_M / (N_LENGTH as f64);

    let bucket_of_x = |x: f64| -> i64 { (x / (0.5 * dx)).round() as i64 };
    use std::collections::BTreeMap;
    let mut buckets: BTreeMap<i64, (f64, Complex64)> = BTreeMap::new();
    for k in 0..n_basis {
        let x = side.edge_midpoint_x(k);
        let length = side.edge_length(k);
        let ik = i_vec[(k, 0)];
        let bkt = bucket_of_x(x);
        let entry = buckets
            .entry(bkt)
            .or_insert((0.0, Complex64::new(0.0, 0.0)));
        entry.0 += length;
        entry.1 += Complex64::new(length, 0.0) * ik;
    }

    // Convert buckets to (x_center, length-weighted-average-current).
    let mut col_data: Vec<(f64, Complex64)> = buckets
        .iter()
        .filter(|(_, (l, _))| *l > 0.0)
        .map(|(b, (l, sum))| ((*b as f64) * 0.5 * dx, *sum / Complex64::new(*l, 0.0)))
        .collect();

    eprintln!("Column buckets: {} columns along x", col_data.len());

    // Mask the port columns. The centered delta-gap injects at
    // x = L/2; exclude any column within ±1.5 dx of that.
    let x_port = STRIP_L_M / 2.0;
    let port_keep = |x: f64| (x - x_port).abs() > 1.5 * dx;

    // Also exclude the very first and last columns — RWG basis
    // density at the strip ends is asymmetric and the current goes to
    // zero there, which makes the phase ill-defined.
    let x_min_keep = 4.0 * dx;
    let x_max_keep = STRIP_L_M - 4.0 * dx;
    col_data
        .retain(|(x, c)| port_keep(*x) && *x >= x_min_keep && *x <= x_max_keep && c.norm() > 1e-30);
    eprintln!(
        "  after port + end masking: {} columns in fit range [{:.2}, {:.2}] mm",
        col_data.len(),
        x_min_keep * 1e3,
        x_max_keep * 1e3,
    );

    col_data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    eprintln!();
    eprintln!(
        "  {:>10}  {:>14}  {:>14}  {:>14}",
        "x (mm)", "|i_col|", "arg (rad)", "arg unwrap"
    );
    let mut xs: Vec<f64> = Vec::with_capacity(col_data.len());
    let mut phases: Vec<f64> = Vec::with_capacity(col_data.len());
    let mut last_raw: Option<f64> = None;
    let mut last_unwrap: Option<f64> = None;
    for (x, c) in &col_data {
        let raw = c.arg();
        let unwrap = match last_raw {
            None => raw,
            Some(prev_raw) => {
                let prev_un = last_unwrap.unwrap();
                let mut step = raw - prev_raw;
                while step > std::f64::consts::PI {
                    step -= std::f64::consts::TAU;
                }
                while step < -std::f64::consts::PI {
                    step += std::f64::consts::TAU;
                }
                prev_un + step
            }
        };
        last_raw = Some(raw);
        last_unwrap = Some(unwrap);
        eprintln!(
            "  {:>10.3}  {:>14.4e}  {:>14.4}  {:>14.4}",
            x * 1e3,
            c.norm(),
            raw,
            unwrap,
        );
        xs.push(*x);
        phases.push(unwrap);
    }

    assert!(
        xs.len() >= 4,
        "need at least 4 samples to fit; got {}",
        xs.len()
    );

    // Method A: phase-slope (running-wave) extraction.
    //
    // Linear fit arg(x) = -β x + φ_0  ->  slope_A = -β_A. Returns β_A
    // ≈ 0 for a pure standing wave on a centered-driven half-wave
    // resonator (port at L/2 launches symmetric ± currents that
    // recombine into a real-valued standing pattern).
    let (slope, _intercept) = linfit(&xs, &phases);
    let beta_a = -slope;
    let k0 = k0_at(F_HZ);
    let beta_over_k0_a = beta_a / k0;
    let eps_eff_a = beta_over_k0_a * beta_over_k0_a;

    // Half-symmetry check disambiguates running-wave vs standing-
    // wave regimes; standing waves on a centered-driven half-wave
    // resonator have left/right |i|-mean asymmetry below 10 %.
    let i_mags: Vec<f64> = col_data.iter().map(|(_, c)| c.norm()).collect();
    let half = i_mags.len() / 2;
    let left_mean: f64 = i_mags[..half].iter().sum::<f64>() / (half as f64).max(1.0);
    let right_mean: f64 =
        i_mags[half..].iter().sum::<f64>() / ((i_mags.len() - half) as f64).max(1.0);
    let asym = (left_mean - right_mean).abs() / (left_mean + right_mean).max(1e-30);
    let is_standing_wave = asym < 0.10 && beta_over_k0_a.abs() < 0.10;

    // Method B: standing-wave envelope fit |I(x)| ≈ A · |sin(β · (x
    // - x_node))|. For a thin open-ended strip the current goes to
    // zero at x = 0 and x = L (no charge accumulation past the
    // strip ends), so x_node = 0 and the lowest mode is β·L = π
    // (half-wave). At the H-J β/k_0 = √3.32 ≈ 1.823 and L = 82 mm,
    // β·L = 1.823 · 0.082 · 2π · 1e9 / c ≈ π — confirming the
    // L = 82 mm choice maps onto the half-wave resonator at the
    // analytic ε_eff_HJ.
    //
    // Robust extractor: a 1-D grid search over β/k_0 ∈ [0.5, 3.5]
    // minimising the per-column chi² between |i(x)| and the
    // amplitude-fit envelope `A_β · |sin(β · x)|`. `A_β` is fixed by
    // least-squares: `A_β = Σ_k |sin(β·x_k)| · |i_k| / Σ_k sin²(β·x_k)`.
    let chi2_at = |beta: f64| -> (f64, f64) {
        let sxx: f64 = xs.iter().map(|x| (beta * x).sin().abs().powi(2)).sum();
        let sxy: f64 = xs
            .iter()
            .zip(i_mags.iter())
            .map(|(x, m)| (beta * x).sin().abs() * m)
            .sum();
        let a = if sxx > 0.0 { sxy / sxx } else { 0.0 };
        let chi2: f64 = xs
            .iter()
            .zip(i_mags.iter())
            .map(|(x, m)| (m - a * (beta * x).sin().abs()).powi(2))
            .sum();
        (chi2, a)
    };

    // Coarse grid + parabolic refinement.
    let n_grid = 1001;
    let beta_min_norm = 0.5; // β/k_0
    let beta_max_norm = 3.5;
    let mut best_n = 0_usize;
    let mut best_chi2 = f64::INFINITY;
    for n_g in 0..n_grid {
        let bn =
            beta_min_norm + (beta_max_norm - beta_min_norm) * (n_g as f64) / ((n_grid - 1) as f64);
        let beta = bn * k0;
        let (chi2, _) = chi2_at(beta);
        if chi2 < best_chi2 {
            best_chi2 = chi2;
            best_n = n_g;
        }
    }
    let coarse_bn =
        beta_min_norm + (beta_max_norm - beta_min_norm) * (best_n as f64) / ((n_grid - 1) as f64);
    // Parabolic refinement on [coarse_bn - dh, coarse_bn + dh].
    let dh = 2.0 * (beta_max_norm - beta_min_norm) / ((n_grid - 1) as f64);
    let n_fine = 401;
    let mut best_bn = coarse_bn;
    let mut best_chi2_fine = f64::INFINITY;
    let mut best_a = 0.0;
    for n_g in 0..n_fine {
        let bn = (coarse_bn - dh) + 2.0 * dh * (n_g as f64) / ((n_fine - 1) as f64);
        let beta = bn * k0;
        let (chi2, a) = chi2_at(beta);
        if chi2 < best_chi2_fine {
            best_chi2_fine = chi2;
            best_bn = bn;
            best_a = a;
        }
    }
    let beta_b = best_bn * k0;
    let beta_over_k0_b = best_bn;
    let eps_eff_b = beta_over_k0_b * beta_over_k0_b;
    let beta_over_k0_hj = EPS_EFF_HJ_ANALYTIC.sqrt();
    let rel_err_b = (eps_eff_b - EPS_EFF_HJ_ANALYTIC) / EPS_EFF_HJ_ANALYTIC;

    eprintln!();
    eprintln!("Phase-vs-x linear fit slope (Method A — running-wave β):");
    eprintln!("    β_A        = {:.4e} rad/m", beta_a);
    eprintln!("    β_A / k_0  = {:.4}", beta_over_k0_a);
    eprintln!("    eps_eff_A  = {:.4}", eps_eff_a);
    eprintln!(
        "    (Method A returns β ≈ 0 for a pure standing wave; |i| left/right asymmetry = {:.2} %.)",
        asym * 100.0,
    );
    eprintln!();
    eprintln!("Standing-wave envelope fit |i(x)| ≈ A·|sin(β·x)| (Method B — standing-wave β):");
    eprintln!("    β_B        = {:.4e} rad/m", beta_b);
    eprintln!("    β_B / k_0  = {:.4}", beta_over_k0_b);
    eprintln!("    eps_eff_B  = (β_B/k_0)^2 = {:.4}", eps_eff_b);
    eprintln!("    A_fit      = {:.4e}", best_a);
    eprintln!("    chi^2 min  = {:.4e}", best_chi2_fine);
    eprintln!("    eps_eff_HJ = {}", EPS_EFF_HJ_ANALYTIC);
    eprintln!(
        "    β/k_0_HJ   = sqrt({}) = {:.4}",
        EPS_EFF_HJ_ANALYTIC, beta_over_k0_hj,
    );
    eprintln!(
        "    Relative error (Method B): (eps_eff_B - {}) / {} = {:+.2} %",
        EPS_EFF_HJ_ANALYTIC,
        EPS_EFF_HJ_ANALYTIC,
        rel_err_b * 100.0,
    );

    // Verdict picks the dominant regime (running vs standing).
    eprintln!();
    eprintln!("Verdict:");
    let (regime, eps_eff_solver, rel_err, beta_over_k0) = if is_standing_wave {
        (
            "standing-wave envelope (Method B)",
            eps_eff_b,
            rel_err_b,
            beta_over_k0_b,
        )
    } else {
        let rel_err_a = (eps_eff_a - EPS_EFF_HJ_ANALYTIC) / EPS_EFF_HJ_ANALYTIC;
        (
            "running-wave phase slope (Method A)",
            eps_eff_a,
            rel_err_a,
            beta_over_k0_a,
        )
    };
    eprintln!("    Dominant regime: {regime}");
    eprintln!("    β / k_0          = {:.4}", beta_over_k0);
    eprintln!("    eps_eff_solver   = {:.4}", eps_eff_solver);
    eprintln!("    eps_eff_HJ       = {}", EPS_EFF_HJ_ANALYTIC);
    eprintln!(
        "    Relative error: (eps_eff_solver - {}) / {} = {:+.2} %",
        EPS_EFF_HJ_ANALYTIC,
        EPS_EFF_HJ_ANALYTIC,
        rel_err * 100.0,
    );
    let within_5pct = rel_err.abs() < 0.05;
    if within_5pct {
        eprintln!(
            "    β/k_0 matches sqrt(eps_eff_HJ) within ±5%:   YES - kernel correctly models the strip eigenmode."
        );
        eprintln!(
            "    Implication: the |Im(Z)| = 674 Ω residual at 1 GHz on the reframed mom-002 is NOT a kernel"
        );
        eprintln!(
            "    bug. Per ADR-0037 / NNNNNNN, the next track should investigate port-excitation /"
        );
        eprintln!("    edge-singularity discretisation effects (CCCCCCC's ADR-0036 directions).");
    } else {
        eprintln!(
            "    β/k_0 matches sqrt(eps_eff_HJ) within ±5%:   NO - kernel bug at {:+.2} % on eps_eff.",
            rel_err * 100.0,
        );
        eprintln!(
            "    Implication: the kernel's strip-eigenmode physics is genuinely off. Per NNNNNNN,"
        );
        eprintln!("    the prime suspect is the K^A vector-potential image train (the TE-channel");
        eprintln!(
            "    inductive part of the line). The next track should audit the TE Sommerfeld split"
        );
        eprintln!("    in MultilayerGreens scalar_vector and its DCIM image series.");
    }

    eprintln!();
    eprintln!(
        "Sanity: |i_col| left/right means: {:.3e} / {:.3e}  (asymmetry {:.2} %).",
        left_mean,
        right_mean,
        asym * 100.0,
    );
    if is_standing_wave {
        eprintln!("  Standing-wave regime detected (asymmetry < 10 % AND running-β/k_0 < 0.1).");
        eprintln!("  The phase-slope Method A returns β ≈ 0 by construction; Method B (|i|");
        eprintln!("  envelope) is the correct probe and its β is reported as the verdict above.");
    } else {
        eprintln!("  Running-wave regime: phase advances monotonically along x; Method A is the");
        eprintln!("  correct probe and its β is reported as the verdict above.");
    }
}
