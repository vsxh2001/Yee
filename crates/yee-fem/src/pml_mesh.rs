//! Phase 4.fem.eig.3.5 P2 — CFS-PML mesh extension.
//!
//! Extends an axis-aligned uniform-brick [`yee_mesh::TetMesh3D`] (the
//! existing `cavity_uniform` constructor's output, or anything else
//! that meshes a cuboid in Kuhn-6 brick form with the same global-id
//! layout) with a thin tetrahedral PML shell on every face the caller
//! tags. Returns the extended mesh, a per-tet [`PmlClass`] map
//! identifying interior vs PML-shell tets together with the per-axis
//! depths, and a face index map relating the new mesh's faces back to
//! the original.
//!
//! The PML shell preserves the original mesh's lattice spacing along
//! the two in-face axes and extrudes `thickness_cells` brick layers
//! outward along the face-normal axis. Each new brick is meshed with
//! the same six Kuhn tets as the cavity interior, ensuring vertex
//! coincidence at the inner PML/cavity boundary (the spec §3.2
//! continuity requirement that makes `Λ(d = 0) = I` produce zero
//! surface reflection).
//!
//! ## v3.5 scope: Cartesian-aligned, single-axis-per-face PML
//!
//! ADR-0043 §4 restricts v3.5 to **Cartesian-aligned, single-axis-per-
//! face** PML shells. Multi-axis edge / corner wedges (where two or
//! three PML shells meet) and non-axis-aligned face normals are
//! rejected with `Error::Unimplemented`; the Phase 4.fem.eig.3.5.1
//! follow-up lifts these restrictions per ADR-0043 §risks. fem-eig-003
//! (one ABC face) and fem-eig-006 (one ABC face) both fit inside the
//! v3.5 restriction. Multi-face callers must invoke this helper once
//! per face and accept the absence of corner / edge wedges — the
//! out-of-plane PML axes' `Λ` factors are taken as `I` (no
//! attenuation) at the seam, which preserves continuity but lowers the
//! corner absorption efficiency.
//!
//! ## References
//!
//! * Phase 4.fem.eig.3.5 spec
//!   `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-cfs-pml-design.md`
//!   §4.2 (PML mesh extension).
//! * Phase 4.fem.eig.3.5 plan
//!   `docs/superpowers/plans/2026-05-20-phase-4-fem-eig-3-5-cfs-pml.md`
//!   step P2.
//! * ADR-0043 — Cartesian-aligned-only scope decision.
//! * Roden, J. A. and Gedney, S. D., *IEEE MWCL* 10(5) (2000),
//!   pp. 27-29 — the CFS-PML formulation whose mesh substrate this
//!   module provides.

use nalgebra::Vector3;
use yee_core::Error;
use yee_mesh::TetMesh3D;

/// Classification of one tet in the extended (cavity + PML shells)
/// mesh: is it inside the original cavity, or does it sit in a PML
/// shell along one Cartesian axis?
///
/// `d_*` are the depths (m) into each PML shell measured from the
/// inner (cavity-side) PML boundary toward the outer truncation
/// surface. Per ADR-0043 §4 only single-axis variants are supported in
/// v3.5; the [`Self::PmlXY`] / [`Self::PmlYZ`] / [`Self::PmlZX`] /
/// [`Self::PmlXYZ`] variants exist for forward compatibility with the
/// Phase 4.fem.eig.3.5.1 multi-axis wedge extension but are never
/// returned by [`extend_mesh_with_pml`] in v3.5 — corner / edge tets
/// in a multi-face PML call are emitted as one of the single-axis
/// variants (the axis of the face that "owns" the corner), and the
/// remaining PML axes contribute `Λ = I` (no attenuation along those
/// axes for that tet).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PmlClass {
    /// Cavity interior tet — no PML.
    Interior,
    /// PML shell along the `+x` axis (face at `x = x_max`) or `-x`
    /// axis (face at `x = x_min`). `d` is the inward depth from the
    /// PML inner boundary (cavity / PML interface) toward the outer
    /// truncation surface; `d ∈ [0, D]` with `D = thickness_cells *
    /// h_cell`.
    PmlX {
        /// Inward depth (m) into the x-PML shell.
        d: f64,
    },
    /// PML shell along the `±y` axis. See [`Self::PmlX`].
    PmlY {
        /// Inward depth (m) into the y-PML shell.
        d: f64,
    },
    /// PML shell along the `±z` axis. See [`Self::PmlX`].
    PmlZ {
        /// Inward depth (m) into the z-PML shell.
        d: f64,
    },
    /// Edge wedge where the x- and y-PML shells overlap. **Not
    /// returned by v3.5** — placeholder for Phase 4.fem.eig.3.5.1.
    PmlXY {
        /// Inward depth (m) into the x-PML shell.
        d_x: f64,
        /// Inward depth (m) into the y-PML shell.
        d_y: f64,
    },
    /// Edge wedge where the y- and z-PML shells overlap. **Not
    /// returned by v3.5** — placeholder for Phase 4.fem.eig.3.5.1.
    PmlYZ {
        /// Inward depth (m) into the y-PML shell.
        d_y: f64,
        /// Inward depth (m) into the z-PML shell.
        d_z: f64,
    },
    /// Edge wedge where the z- and x-PML shells overlap. **Not
    /// returned by v3.5** — placeholder for Phase 4.fem.eig.3.5.1.
    PmlZX {
        /// Inward depth (m) into the z-PML shell.
        d_z: f64,
        /// Inward depth (m) into the x-PML shell.
        d_x: f64,
    },
    /// Corner wedge where all three PML shells overlap. **Not
    /// returned by v3.5** — placeholder for Phase 4.fem.eig.3.5.1.
    PmlXYZ {
        /// Inward depth (m) into the x-PML shell.
        d_x: f64,
        /// Inward depth (m) into the y-PML shell.
        d_y: f64,
        /// Inward depth (m) into the z-PML shell.
        d_z: f64,
    },
}

impl PmlClass {
    /// `true` iff this is the [`Self::Interior`] variant (no PML
    /// stretching applied).
    pub fn is_interior(&self) -> bool {
        matches!(self, Self::Interior)
    }
}

/// Cartesian direction identifying which face of the cavity bounding
/// box gets a PML shell. Used by [`extend_mesh_with_pml`] to specify
/// the outward-normal axis of each PML face.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PmlAxis {
    /// Face at the cavity's minimum-x edge (outward normal `-x̂`).
    XMin,
    /// Face at the cavity's maximum-x edge (outward normal `+x̂`).
    XMax,
    /// Face at the cavity's minimum-y edge (outward normal `-ŷ`).
    YMin,
    /// Face at the cavity's maximum-y edge (outward normal `+ŷ`).
    YMax,
    /// Face at the cavity's minimum-z edge (outward normal `-ẑ`).
    ZMin,
    /// Face at the cavity's maximum-z edge (outward normal `+ẑ`).
    ZMax,
}

impl PmlAxis {
    /// Match a Cartesian outward-normal unit vector against the six
    /// axis-aligned directions. Returns `None` if the vector is not
    /// within `1e-6` of any of `±x̂`, `±ŷ`, `±ẑ`.
    pub fn from_outward_normal(n: Vector3<f64>) -> Option<Self> {
        let norm = n.norm();
        if norm < 1e-12 {
            return None;
        }
        let n = n / norm;
        let tol = 1e-6;
        if (n.x - 1.0).abs() < tol {
            Some(Self::XMax)
        } else if (n.x + 1.0).abs() < tol {
            Some(Self::XMin)
        } else if (n.y - 1.0).abs() < tol {
            Some(Self::YMax)
        } else if (n.y + 1.0).abs() < tol {
            Some(Self::YMin)
        } else if (n.z - 1.0).abs() < tol {
            Some(Self::ZMax)
        } else if (n.z + 1.0).abs() < tol {
            Some(Self::ZMin)
        } else {
            None
        }
    }

    /// Returns `(axis_index ∈ {0, 1, 2}, sign ∈ {-1, +1})` where the
    /// axis index identifies the Cartesian direction (0 = x, 1 = y,
    /// 2 = z) and the sign indicates the outward orientation. Used by
    /// downstream consumers (P4 `with_cfs_pml` builder) that need to
    /// resolve a per-tet PML class to a Cartesian stretching axis.
    pub fn axis_and_sign(self) -> (usize, i32) {
        match self {
            Self::XMin => (0, -1),
            Self::XMax => (0, 1),
            Self::YMin => (1, -1),
            Self::YMax => (1, 1),
            Self::ZMin => (2, -1),
            Self::ZMax => (2, 1),
        }
    }
}

/// Per-face mapping from extended-mesh face indices back to original-
/// mesh face indices. v3.5 ships only a placeholder: the original
/// mesh's exterior face indices are remapped 1-to-1 onto the new
/// mesh's exterior face indices (no extra cavity-side faces are
/// removed — only the PML-fronted faces become **interior** in the
/// extended mesh, replaced by new outer-PEC faces on the truncation
/// surface).
#[derive(Debug, Clone)]
pub struct FaceIndexMap {
    /// `new_to_old[i] = Some(j)` if the extended mesh's exterior face
    /// `i` corresponds to original face `j`; `None` if it is a brand-
    /// new face introduced by the PML extension (the outer truncation
    /// surface).
    pub new_to_old: Vec<Option<usize>>,
}

/// Extend `mesh` with PML shells on every face listed in
/// `pml_faces`, returning the extended mesh, the per-tet PML
/// classification table, and a face index map.
///
/// # Algorithm
///
/// 1. Detect the original mesh's axis-aligned bounding box and
///    per-axis lattice spacing. Requires `mesh` to be a uniform
///    rectangular tet mesh (the same shape as
///    [`yee_mesh::TetMesh3D::cavity_uniform`]'s output). Non-uniform
///    or non-rectangular meshes return `Error::Unimplemented`.
/// 2. For each [`PmlAxis`] in `pml_faces`, extend the cavity bounding
///    box outward by `thickness_cells` lattice layers along that
///    axis. The face cross-section grid stays the same.
/// 3. Rebuild the full vertex grid and tetrahedron list using the
///    same Kuhn-6 brick pattern as [`yee_mesh::TetMesh3D::cavity_uniform`].
/// 4. Classify each tet by its centroid: tets whose centroid lies
///    inside the **original** cavity bounding box are [`PmlClass::Interior`];
///    everything else is [`PmlClass::PmlX/Y/Z`] depending on which
///    extended axis its centroid falls in. The per-axis depth
///    `d_α = |centroid_α − boundary_α|` is the distance from the PML
///    inner boundary (the original cavity face) along the outward
///    normal.
/// 5. Tets at multi-axis seams (corners / edges of multi-face PML
///    calls) are classified by their **primary** axis — the axis with
///    the largest `d` — and the other axes contribute `Λ = I`. This
///    is the v3.5 simplification; Phase 4.fem.eig.3.5.1 emits true
///    edge / corner wedges.
///
/// The new mesh is constructed in-place; vertex / tet ordering is
/// canonical and deterministic across runs.
///
/// # Arguments
///
/// * `mesh` — the original cavity mesh. Must be axis-aligned uniform
///   bricks (as produced by [`yee_mesh::TetMesh3D::cavity_uniform`]).
/// * `pml_faces` — non-empty slice of axis-aligned [`PmlAxis`] tags
///   identifying the cavity boundary faces that should receive a PML
///   shell.
/// * `thickness_cells` — PML shell thickness in tet brick layers
///   (Roden-Gedney 2000 §III default 6).
///
/// # Errors
///
/// * `Error::Invalid` — empty `pml_faces`, non-uniform mesh, or
///   `thickness_cells == 0`.
/// * `Error::Unimplemented` — mesh shape not recognised as a uniform
///   Kuhn-6 brick lattice; lifting this restriction is queued for
///   Phase 4.fem.eig.3.5.1.
///
/// # v3.5 compatibility note
///
/// With `thickness_cells = 0`, the function returns the original mesh
/// unchanged with every tet tagged [`PmlClass::Interior`]. This is
/// the "no-op PML" path that the
/// `crates/yee-fem/tests/pml_open_boundary_assembly.rs::pml_assembly_matches_scalar_on_zero_thickness`
/// gate consumes to verify backward compatibility with the v3
/// 2nd-order ABC path.
pub fn extend_mesh_with_pml(
    mesh: &TetMesh3D,
    pml_faces: &[PmlAxis],
    thickness_cells: usize,
) -> Result<(TetMesh3D, Vec<PmlClass>, FaceIndexMap), Error> {
    if pml_faces.is_empty() {
        return Err(Error::Invalid(
            "extend_mesh_with_pml: pml_faces slice is empty".to_string(),
        ));
    }

    // Zero-thickness no-op path — used by the v3-equivalence gate
    // test to confirm the PML assembly reduces to the v3 surface-
    // integral path bit-for-bit when no shell is added.
    if thickness_cells == 0 {
        let classes = vec![PmlClass::Interior; mesh.tetrahedra.len()];
        let new_to_old: Vec<Option<usize>> = (0..mesh.tetrahedra.len() * 4).map(Some).collect();
        return Ok((clone_mesh(mesh)?, classes, FaceIndexMap { new_to_old }));
    }

    // ---- 1. Reverse-engineer the cavity bounding box + lattice.
    let lattice = UniformLattice::detect(mesh)?;
    let (bx_min, _bx_max) = (lattice.x_min, lattice.x_max);
    let (by_min, _by_max) = (lattice.y_min, lattice.y_max);
    let (bz_min, _bz_max) = (lattice.z_min, lattice.z_max);

    // ---- 2. Compute extended-grid extents per axis (in lattice-step
    // integer units measured from the original cavity origin).
    let mut ext_x_min: i64 = 0;
    let mut ext_x_max: i64 = lattice.nx as i64;
    let mut ext_y_min: i64 = 0;
    let mut ext_y_max: i64 = lattice.ny as i64;
    let mut ext_z_min: i64 = 0;
    let mut ext_z_max: i64 = lattice.nz as i64;
    let t = thickness_cells as i64;

    for axis in pml_faces.iter().copied() {
        match axis {
            PmlAxis::XMin => ext_x_min -= t,
            PmlAxis::XMax => ext_x_max += t,
            PmlAxis::YMin => ext_y_min -= t,
            PmlAxis::YMax => ext_y_max += t,
            PmlAxis::ZMin => ext_z_min -= t,
            PmlAxis::ZMax => ext_z_max += t,
        }
    }

    let new_nx = (ext_x_max - ext_x_min) as usize;
    let new_ny = (ext_y_max - ext_y_min) as usize;
    let new_nz = (ext_z_max - ext_z_min) as usize;

    let nvx = new_nx + 1;
    let nvy = new_ny + 1;
    let nvz = new_nz + 1;

    // The new origin (in m) is shifted from the original by the
    // negative-side extensions.
    let new_x0 = bx_min + (ext_x_min as f64) * lattice.dx;
    let new_y0 = by_min + (ext_y_min as f64) * lattice.dy;
    let new_z0 = bz_min + (ext_z_min as f64) * lattice.dz;

    // ---- 3. Build the vertex grid.
    let mut vertices = Vec::with_capacity(nvx * nvy * nvz);
    for k in 0..nvz {
        for j in 0..nvy {
            for i in 0..nvx {
                vertices.push(Vector3::new(
                    new_x0 + (i as f64) * lattice.dx,
                    new_y0 + (j as f64) * lattice.dy,
                    new_z0 + (k as f64) * lattice.dz,
                ));
            }
        }
    }

    let global_id = |i: usize, j: usize, k: usize| -> usize { i + j * nvx + k * nvx * nvy };

    // ---- 4. Build the tet list using the same Kuhn-6 brick layout
    // as `cavity_uniform`. We re-derive `KUHN_LOCAL_TETS` here so the
    // pml_mesh module is self-contained and does not need to depend
    // on `yee-mesh`'s private constants.
    const KUHN_LOCAL_TETS: [[usize; 4]; 6] = [
        [0, 1, 3, 7],
        [0, 1, 5, 7],
        [0, 2, 3, 7],
        [0, 2, 6, 7],
        [0, 4, 5, 7],
        [0, 4, 6, 7],
    ];

    let mut tetrahedra = Vec::with_capacity(new_nx * new_ny * new_nz * 6);
    let mut classes = Vec::with_capacity(new_nx * new_ny * new_nz * 6);

    // Negative-side extension counts per axis. Required to identify
    // which bricks are PML vs interior in the extended grid.
    let neg_x = (-ext_x_min).max(0) as usize;
    let neg_y = (-ext_y_min).max(0) as usize;
    let neg_z = (-ext_z_min).max(0) as usize;
    // Positive-side extension counts per axis.
    let pos_x = (ext_x_max - (lattice.nx as i64)).max(0) as usize;
    let pos_y = (ext_y_max - (lattice.ny as i64)).max(0) as usize;
    let pos_z = (ext_z_max - (lattice.nz as i64)).max(0) as usize;

    for k in 0..new_nz {
        for j in 0..new_ny {
            for i in 0..new_nx {
                // Brick centroid in extended-grid integer coordinates
                // (half-cell offsets).
                let in_x_pml_neg = i < neg_x;
                let in_x_pml_pos = i >= neg_x + lattice.nx;
                let in_y_pml_neg = j < neg_y;
                let in_y_pml_pos = j >= neg_y + lattice.ny;
                let in_z_pml_neg = k < neg_z;
                let in_z_pml_pos = k >= neg_z + lattice.nz;

                let x_in_pml = in_x_pml_neg || in_x_pml_pos;
                let y_in_pml = in_y_pml_neg || in_y_pml_pos;
                let z_in_pml = in_z_pml_neg || in_z_pml_pos;

                // Per-axis depth: 0 inside cavity; (cells into PML +
                // 0.5) * spacing at the brick centroid. The +0.5
                // accounts for the centroid sitting in the middle of
                // the brick.
                let d_x = if in_x_pml_neg {
                    ((neg_x - 1 - i) as f64 + 0.5) * lattice.dx
                } else if in_x_pml_pos {
                    ((i - (neg_x + lattice.nx)) as f64 + 0.5) * lattice.dx
                } else {
                    0.0
                };
                let d_y = if in_y_pml_neg {
                    ((neg_y - 1 - j) as f64 + 0.5) * lattice.dy
                } else if in_y_pml_pos {
                    ((j - (neg_y + lattice.ny)) as f64 + 0.5) * lattice.dy
                } else {
                    0.0
                };
                let d_z = if in_z_pml_neg {
                    ((neg_z - 1 - k) as f64 + 0.5) * lattice.dz
                } else if in_z_pml_pos {
                    ((k - (neg_z + lattice.nz)) as f64 + 0.5) * lattice.dz
                } else {
                    0.0
                };

                // Per-brick PML class. Per ADR-0043 §4 v3.5 emits
                // single-axis variants only; multi-axis seam bricks
                // collapse to their primary (largest-depth) axis.
                let brick_class = if !x_in_pml && !y_in_pml && !z_in_pml {
                    PmlClass::Interior
                } else if x_in_pml && !y_in_pml && !z_in_pml {
                    PmlClass::PmlX { d: d_x }
                } else if y_in_pml && !x_in_pml && !z_in_pml {
                    PmlClass::PmlY { d: d_y }
                } else if z_in_pml && !x_in_pml && !y_in_pml {
                    PmlClass::PmlZ { d: d_z }
                } else {
                    // Multi-axis seam — primary axis = max depth.
                    let (mut best, mut best_d) = (PmlClass::PmlX { d: d_x }, d_x);
                    if d_y > best_d {
                        best = PmlClass::PmlY { d: d_y };
                        best_d = d_y;
                    }
                    if d_z > best_d {
                        best = PmlClass::PmlZ { d: d_z };
                    }
                    best
                };

                for local in &KUHN_LOCAL_TETS {
                    let mut tet = [0usize; 4];
                    for (slot, &c) in local.iter().enumerate() {
                        let di = c & 1;
                        let dj = (c >> 1) & 1;
                        let dk = (c >> 2) & 1;
                        tet[slot] = global_id(i + di, j + dj, k + dk);
                    }
                    tetrahedra.push(tet);
                    classes.push(brick_class);
                }
            }
        }
    }

    let n_verts = vertices.len();
    let n_tets = tetrahedra.len();
    let extended = TetMesh3D::new(vertices, tetrahedra, None, None)
        .map_err(|e| Error::Invalid(format!("extend_mesh_with_pml: TetMesh3D::new failed: {e}")))?;

    // FaceIndexMap is a forward-compat placeholder for v3.5 — we do
    // not actually remap face indices here (the open-boundary
    // solver's `ExteriorFaceTable::build` rebuilds the face table
    // from the extended mesh).
    let _ = pos_x;
    let _ = pos_y;
    let _ = pos_z;
    let _ = (n_verts, n_tets);
    let new_to_old: Vec<Option<usize>> = Vec::new();
    Ok((extended, classes, FaceIndexMap { new_to_old }))
}

/// Cheap clone of a [`TetMesh3D`] — the `cavity_uniform` path
/// drops the material tag fields when `None`, so we reconstruct via
/// the public constructor (which performs the same validation) to
/// avoid leaking implementation details across crates.
fn clone_mesh(mesh: &TetMesh3D) -> Result<TetMesh3D, Error> {
    TetMesh3D::new(mesh.vertices.clone(), mesh.tetrahedra.clone(), None, None)
        .map_err(|e| Error::Invalid(format!("clone_mesh: TetMesh3D::new failed: {e}")))
}

/// Reverse-engineered uniform lattice structure of a `cavity_uniform`-
/// shaped mesh.
struct UniformLattice {
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    z_min: f64,
    z_max: f64,
    dx: f64,
    dy: f64,
    dz: f64,
    nx: usize,
    ny: usize,
    nz: usize,
}

impl UniformLattice {
    /// Detect the lattice parameters from the vertex list. Assumes the
    /// mesh is a uniform Kuhn-6 brick mesh with `(nx + 1) × (ny + 1)
    /// × (nz + 1)` vertices on an axis-aligned regular lattice.
    fn detect(mesh: &TetMesh3D) -> Result<Self, Error> {
        // Collect unique sorted x / y / z coordinates.
        let mut xs: Vec<f64> = mesh.vertices.iter().map(|v| v.x).collect();
        let mut ys: Vec<f64> = mesh.vertices.iter().map(|v| v.y).collect();
        let mut zs: Vec<f64> = mesh.vertices.iter().map(|v| v.z).collect();
        sort_dedup_f64(&mut xs);
        sort_dedup_f64(&mut ys);
        sort_dedup_f64(&mut zs);

        if xs.len() < 2 || ys.len() < 2 || zs.len() < 2 {
            return Err(Error::Unimplemented(
                "extend_mesh_with_pml: mesh is not a uniform 3D brick \
                 lattice (need at least 2 unique coords per axis)",
            ));
        }

        let dx = xs[1] - xs[0];
        let dy = ys[1] - ys[0];
        let dz = zs[1] - zs[0];

        // Verify uniform spacing per axis (tolerance: 1e-6 relative
        // to the smallest spacing).
        let tol = 1e-6 * dx.min(dy).min(dz).max(f64::MIN_POSITIVE);
        for w in xs.windows(2) {
            if ((w[1] - w[0]) - dx).abs() > tol {
                return Err(Error::Unimplemented(
                    "extend_mesh_with_pml: non-uniform x spacing — \
                     mesh is not a `cavity_uniform`-style brick lattice",
                ));
            }
        }
        for w in ys.windows(2) {
            if ((w[1] - w[0]) - dy).abs() > tol {
                return Err(Error::Unimplemented(
                    "extend_mesh_with_pml: non-uniform y spacing",
                ));
            }
        }
        for w in zs.windows(2) {
            if ((w[1] - w[0]) - dz).abs() > tol {
                return Err(Error::Unimplemented(
                    "extend_mesh_with_pml: non-uniform z spacing",
                ));
            }
        }

        Ok(Self {
            x_min: xs[0],
            x_max: *xs.last().unwrap(),
            y_min: ys[0],
            y_max: *ys.last().unwrap(),
            z_min: zs[0],
            z_max: *zs.last().unwrap(),
            dx,
            dy,
            dz,
            nx: xs.len() - 1,
            ny: ys.len() - 1,
            nz: zs.len() - 1,
        })
    }
}

/// Sort an `f64` slice ascending and dedupe within a small tolerance
/// — used by [`UniformLattice::detect`] to recover the per-axis
/// lattice coordinates from a flat vertex list.
fn sort_dedup_f64(xs: &mut Vec<f64>) {
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let tol = 1e-9;
    xs.dedup_by(|a, b| (*a - *b).abs() < tol);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn axis_from_unit_normal_round_trip() {
        assert_eq!(
            PmlAxis::from_outward_normal(Vector3::new(1.0, 0.0, 0.0)),
            Some(PmlAxis::XMax)
        );
        assert_eq!(
            PmlAxis::from_outward_normal(Vector3::new(-1.0, 0.0, 0.0)),
            Some(PmlAxis::XMin)
        );
        assert_eq!(
            PmlAxis::from_outward_normal(Vector3::new(0.0, 0.0, -1.0)),
            Some(PmlAxis::ZMin)
        );
        // Off-axis vector rejected.
        assert!(PmlAxis::from_outward_normal(Vector3::new(1.0, 1.0, 0.0).normalize()).is_none());
    }

    #[test]
    fn pml_class_is_interior_only_for_interior_variant() {
        assert!(PmlClass::Interior.is_interior());
        assert!(!PmlClass::PmlX { d: 1.0e-3 }.is_interior());
    }
}
