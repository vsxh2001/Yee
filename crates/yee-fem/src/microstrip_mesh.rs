//! Layered straight-microstrip tetrahedral mesh helper for the FEM-EM
//! driven-sweep track (FEM-EM brick 2, ADR-0153).
//!
//! [`layered_microstrip_mesh`] builds a Kuhn-6 tetrahedralised box
//! ([`yee_mesh::TetMesh3D::cavity_uniform`]), repaints the per-tet
//! material tags so the substrate slab is FR-4 and the volume above it is
//! air, and returns two edge-endpoint predicates — `ground_pred` and
//! `trace_pred` — shaped to feed FEM-EM brick 1's
//! [`crate::OpenBoundarySolver::interior_edges_matching`] /
//! [`crate::OpenBoundarySolver::with_interior_pec_edges`] interior-PEC
//! picker. The output is the geometry input for the eventual
//! [`crate::OpenBoundarySolver::sweep_matrix`] driven S-parameter sweep;
//! this helper performs **no solve**.
//!
//! ## Coordinate / z-stack convention (ADR-0108)
//!
//! The substrate-normal axis is **z**; the propagation axis is **y**; the
//! cross-section width is **x**. The z-stack matches `yee-voxel`'s
//! `voxelize_microstrip` (ADR-0108) exactly — there is **no air series
//! gap** between the ground plane and the dielectric:
//!
//! ```text
//!   z = 0          ground plane (k = 0)            — ground_pred footprint
//!   z ∈ [0, sub_h] substrate slab, ε_r = FR-4      — k = 0 .. n_sub cells
//!   z = sub_h      trace (k_top = n_sub)           — trace_pred footprint
//!   z ∈ [sub_h, box_h] air, ε_r = 1                — k = n_sub .. nz cells
//! ```
//!
//! A mesh plane is forced to land **exactly** on `z = sub_h` so the
//! substrate/air interface is a cell boundary, not a mid-cell — the
//! helper returns [`Error::Invalid`] if the requested `(sub_h, box_h,
//! nz)` triple does not place an integer number of z-cells inside the
//! substrate (mirroring `cavity_uniform`'s own validation style).
//!
//! ## Material tags
//!
//! Per-tet `ε_r` is assigned by **centroid**: a tet whose centroid
//! satisfies `centroid.z < sub_h` is tagged [`FR4_TAG`] (`1`) and resolves
//! to FR-4 (`ε_∞ = 4.4`); every other tet is tagged [`AIR_TAG`] (`0`) and
//! resolves to air (`ε_∞ = 1`). This is the exact centroid-repaint pattern
//! used by `crates/yee-fem/tests/dispersive_solve.rs`. Per-tet `ε_r` is
//! already threaded through assembly via
//! [`crate::MaterialDatabase::eps_at`]; this helper only sets the mesh
//! tags and builds the database.
//!
//! ## Predicate shape
//!
//! Both returned closures take the **two endpoint world coordinates** of a
//! mesh edge (`Fn(Vector3<f64>, Vector3<f64>) -> bool`) — the exact
//! signature [`crate::OpenBoundarySolver::interior_edges_matching`]
//! consumes — and return `true` iff **both** endpoints lie on the relevant
//! conductor footprint:
//!
//! - `ground_pred` — both endpoints on the `z ≈ 0` ground plane (full
//!   box footprint).
//! - `trace_pred` — both endpoints on the `z ≈ sub_h` plane **and** within
//!   the trace's `x`-window `[x_lo, x_hi]` (the strip runs the full `y`
//!   extent, so `y` is unconstrained).
//!
//! The trace `x`-window is snapped to lattice lines (the trace edges fall
//! on `x = i · dx` grid columns), so the picker returns the in-window edges
//! exactly. See [`layered_microstrip_mesh`] for the centring rule.

use nalgebra::Vector3;
use yee_mesh::{Error, MaterialTag, TetMesh3D};

use crate::material::{Material, MaterialDatabase};

/// Material tag applied to substrate (FR-4) tetrahedra — those whose
/// centroid lies below the substrate/air interface (`centroid.z < sub_h`).
pub const FR4_TAG: MaterialTag = 1;

/// Material tag applied to air tetrahedra — those whose centroid lies on
/// or above the substrate/air interface (`centroid.z >= sub_h`). Matches
/// the [`Material::default`] free-space tag convention.
pub const AIR_TAG: MaterialTag = 0;

/// FR-4 high-frequency relative permittivity used for the substrate slab.
pub const FR4_EPS_R: f64 = 4.4;

/// Build a layered straight-microstrip tetrahedral volume mesh plus the
/// per-tet material database and the ground/trace interior-PEC edge
/// predicates (FEM-EM brick 2, ADR-0153).
///
/// The box spans `[0, box_w] × [0, line_len] × [0, box_h]` in `(x, y, z)`,
/// tetrahedralised by [`TetMesh3D::cavity_uniform`] into
/// `nx · ny · nz · 6` Kuhn tets. The substrate-normal axis is `z`; the
/// propagation axis (down the line) is `y`; the cross-section width is
/// `x`. A mesh plane is forced onto `z = sub_h` (see Errors).
///
/// # Arguments
///
/// * `box_w` — total box width along `x` (m).
/// * `box_h` — total box height along `z` (m); spans substrate + air.
/// * `line_len` — line length along `y` (m), the propagation direction.
/// * `sub_h` — substrate thickness (m); the dielectric fills `z ∈ [0,
///   sub_h]`. Must be an integer number of z-cells and strictly inside the
///   box (`0 < sub_h < box_h`).
/// * `trace_w` — trace width along `x` (m); the trace strip is centred on
///   `x = box_w / 2` and snapped to the nearest lattice column boundaries.
///   Must be `> 0` and `<= box_w`.
/// * `nx`, `ny`, `nz` — Kuhn-brick subdivisions along `x`, `y`, `z`. All
///   `>= 1`.
///
/// # Returns
///
/// `(mesh, material_db, ground_pred, trace_pred)` where:
///
/// * `mesh` — the tetrahedralised box with `tetrahedron_material` repainted
///   so substrate tets carry [`FR4_TAG`] and air tets carry [`AIR_TAG`].
/// * `material_db` — a [`MaterialDatabase`] with [`FR4_TAG`] → FR-4
///   (`ε_∞ = 4.4`, `μ_r = 1`, no poles) and [`AIR_TAG`] → air
///   (`ε_∞ = 1`).
/// * `ground_pred` — `Fn(Vector3<f64>, Vector3<f64>) -> bool` selecting
///   edges with both endpoints on `z ≈ 0`.
/// * `trace_pred` — `Fn(Vector3<f64>, Vector3<f64>) -> bool` selecting
///   edges with both endpoints on `z ≈ sub_h` and within the trace
///   `x`-window.
///
/// Feed the predicates through
/// [`crate::OpenBoundarySolver::interior_edges_matching`] then
/// [`crate::OpenBoundarySolver::with_interior_pec_edges`] to tag the ground
/// plane and trace as PEC.
///
/// # Errors
///
/// Returns [`Error::Invalid`] if:
///
/// * any extent (`box_w`, `box_h`, `line_len`, `sub_h`, `trace_w`) is not
///   strictly positive and finite, or any of `nx`, `ny`, `nz` is zero
///   (these are also re-checked by the inner `cavity_uniform`);
/// * `sub_h >= box_h` (the substrate must leave room for air above);
/// * `trace_w > box_w`;
/// * the substrate/air interface `z = sub_h` does not land on a mesh
///   plane — i.e. `sub_h` is not an integer multiple of `dz = box_h / nz`
///   (checked to a relative tolerance of `1e-9 · dz`). Choose `nz` so
///   `sub_h · nz / box_h` is integral.
///
/// # Examples
///
/// ```
/// use yee_fem::microstrip_mesh::{layered_microstrip_mesh, FR4_TAG, AIR_TAG};
/// // 4 mm × 4 mm box, 10 mm line, 1 mm substrate, 0.5 mm trace.
/// // box_h = 4 mm, nz = 8 → dz = 0.5 mm, so z = 1 mm is the 2nd plane.
/// let (mesh, db, ground, trace) =
///     layered_microstrip_mesh(4e-3, 4e-3, 10e-3, 1e-3, 0.5e-3, 4, 10, 8).unwrap();
/// assert_eq!(mesh.n_tets(), 4 * 10 * 8 * 6);
/// // Substrate is the bottom 1/4 of the box height (1 mm of 4 mm).
/// let n_fr4 = mesh
///     .tetrahedron_material
///     .iter()
///     .filter(|&&t| t == FR4_TAG)
///     .count();
/// assert_eq!(n_fr4, mesh.n_tets() / 4);
/// # let _ = (db, ground, trace);
/// ```
// The eight geometry/subdivision arguments are the irreducible inputs of
// a layered straight-microstrip box (three box extents, substrate height,
// trace width, three subdivisions); a parameter struct would only move the
// same eight fields one call away, so the lint is allowed at the boundary.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::type_complexity)]
pub fn layered_microstrip_mesh(
    box_w: f64,
    box_h: f64,
    line_len: f64,
    sub_h: f64,
    trace_w: f64,
    nx: usize,
    ny: usize,
    nz: usize,
) -> Result<
    (
        TetMesh3D,
        MaterialDatabase,
        impl Fn(Vector3<f64>, Vector3<f64>) -> bool,
        impl Fn(Vector3<f64>, Vector3<f64>) -> bool,
    ),
    Error,
> {
    // Extent / subdivision validation (mirrors cavity_uniform's style so
    // the error messages are consistent; cavity_uniform re-checks too).
    for (name, v) in [
        ("box_w", box_w),
        ("box_h", box_h),
        ("line_len", line_len),
        ("sub_h", sub_h),
        ("trace_w", trace_w),
    ] {
        if !(v.is_finite() && v > 0.0) {
            return Err(Error::Invalid(format!(
                "layered_microstrip_mesh: {name} must be a positive finite extent, got {v}"
            )));
        }
    }
    if nx == 0 || ny == 0 || nz == 0 {
        return Err(Error::Invalid(format!(
            "layered_microstrip_mesh: nx, ny, nz must all be >= 1, got ({nx}, {ny}, {nz})"
        )));
    }
    if sub_h >= box_h {
        return Err(Error::Invalid(format!(
            "layered_microstrip_mesh: sub_h ({sub_h}) must be strictly less than box_h \
             ({box_h}) to leave room for air above the substrate"
        )));
    }
    if trace_w > box_w {
        return Err(Error::Invalid(format!(
            "layered_microstrip_mesh: trace_w ({trace_w}) must be <= box_w ({box_w})"
        )));
    }

    // The substrate/air interface must land on a mesh plane: z = sub_h
    // must be an integer multiple of dz = box_h / nz (ADR-0108: the
    // interface is a cell boundary, never mid-cell).
    let dz = box_h / nz as f64;
    let n_sub_f = sub_h / dz;
    let n_sub = n_sub_f.round();
    if (n_sub_f - n_sub).abs() > 1e-9 * (1.0 + n_sub.abs()) {
        return Err(Error::Invalid(format!(
            "layered_microstrip_mesh: substrate/air interface z = sub_h ({sub_h}) does not \
             land on a mesh plane — sub_h / dz = {n_sub_f} is not integral (dz = box_h / nz = \
             {dz}). Choose nz so sub_h * nz / box_h is an integer."
        )));
    }

    // Box: x ∈ [0, box_w] (width), y ∈ [0, line_len] (propagation),
    // z ∈ [0, box_h] (substrate-normal). Kuhn-6 per brick.
    let mut mesh = TetMesh3D::cavity_uniform(box_w, line_len, box_h, nx, ny, nz)?;

    // Per-tet ε_r repaint by centroid (mirrors
    // `crates/yee-fem/tests/dispersive_solve.rs` centroid-repaint): a tet
    // whose centroid is below the substrate/air interface is FR-4 (tag 1),
    // else air (tag 0).
    for (tet_idx, tet) in mesh.tetrahedra.iter().enumerate() {
        let centroid_z = 0.25
            * (mesh.vertices[tet[0]].z
                + mesh.vertices[tet[1]].z
                + mesh.vertices[tet[2]].z
                + mesh.vertices[tet[3]].z);
        mesh.tetrahedron_material[tet_idx] = if centroid_z < sub_h { FR4_TAG } else { AIR_TAG };
    }

    let material_db = MaterialDatabase::new()
        .with_material(
            FR4_TAG,
            Material {
                eps_inf: FR4_EPS_R,
                mu_r: 1.0,
                poles: vec![],
            },
        )
        .with_material(AIR_TAG, Material::default());

    // Trace x-window, centred on the box and snapped to lattice columns so
    // the trace edges fall on x = i·dx grid lines. With dx = box_w / nx,
    // the half-width in columns is round((trace_w / 2) / dx), clamped to
    // [1, ...] so the window is never empty (at least one column-pair
    // straddling the centre) but never wider than the box.
    let dx = box_w / nx as f64;
    let centre_col_f = (box_w / 2.0) / dx;
    // Snap the trace centre to the nearest lattice column index, then take
    // a symmetric half-width in whole columns.
    let half_cols = ((trace_w / 2.0) / dx).round().max(1.0);
    let lo_col = (centre_col_f - half_cols).round().max(0.0);
    let hi_col = (centre_col_f + half_cols).round().min(nx as f64);
    let x_lo = lo_col * dx;
    let x_hi = hi_col * dx;

    // Geometric tolerances: a small fraction of the smallest in-plane cell
    // so an endpoint exactly on a lattice line passes but the next line
    // over does not.
    let z_tol = 1e-6 * box_h;
    let x_tol = 1e-6 * box_w;

    let ground_pred =
        move |a: Vector3<f64>, b: Vector3<f64>| -> bool { a.z.abs() < z_tol && b.z.abs() < z_tol };

    let trace_pred = move |a: Vector3<f64>, b: Vector3<f64>| -> bool {
        let on_top = (a.z - sub_h).abs() < z_tol && (b.z - sub_h).abs() < z_tol;
        let in_x = a.x >= x_lo - x_tol
            && a.x <= x_hi + x_tol
            && b.x >= x_lo - x_tol
            && b.x <= x_hi + x_tol;
        on_top && in_x
    };

    Ok((mesh, material_db, ground_pred, trace_pred))
}
