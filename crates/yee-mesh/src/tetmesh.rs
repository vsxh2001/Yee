//! Three-dimensional tetrahedral mesh — the input contract for the
//! Phase 4 FEM eigenmode solver in `yee-fem`.
//!
//! Mirror of [`crate::TriMesh2D`] one dimension up: vertices in
//! `(x, y, z)`, tetrahedra (4 vertex indices each), per-vertex and
//! per-tet [`crate::MaterialTag`]s for `eps_r` / `mu_r` lookup during
//! Nedelec element-matrix assembly.
//!
//! **Status:** Phase 4.fem.eig.0 step T2. Construction enforces positive
//! signed volume on every tet (auto-reorienting silently if a caller
//! hands us a tet with reversed orientation — easier than rejecting
//! outright, and consistent with first-order Nedelec which has no
//! orientation-of-the-cell preference). The element-matrix assembly
//! lands in step T3+ of the SSSSS plan
//! (`docs/superpowers/plans/2026-05-18-phase-4-fem-eigenmode.md`).
//!
//! The uniform-cavity constructor `cavity_uniform(a, b, d, nx, ny, nz)`
//! using a Kuhn 6-tet decomposition is **deferred to step T6** of the
//! plan; it is a separate brick of work with its own validation
//! (cell-count, total-volume, boundary-edge closure) and is dispatched
//! as a follow-up track.

use crate::{Error, MaterialTag, Result};
use nalgebra::Vector3;

/// A three-dimensional tetrahedral mesh.
///
/// **Field invariants** (enforced by [`Self::new`]):
///
/// * `vertices.len() >= 4` (need at least one tet's worth of vertices).
/// * `tetrahedra.len() >= 1`.
/// * Every tet index is `< vertices.len()`.
/// * Every tet has strictly positive signed volume `V_e =
///   (1/6) · (v₁ − v₀) · ((v₂ − v₀) × (v₃ − v₀)) > ε`. Tets whose
///   computed volume is negative are auto-reoriented by swapping
///   `v₂ ↔ v₃` (the silent-fix path documented above). Tets whose
///   volume is exactly zero or non-finite are rejected outright.
/// * Material-tag lengths match `vertices` and `tetrahedra` respectively
///   when `Some(_)`; `None` defaults to `vec![0; n]`.
#[derive(Debug, Default, Clone)]
pub struct TetMesh3D {
    /// Vertices in world coordinates `(x, y, z)`, metres.
    pub vertices: Vec<Vector3<f64>>,
    /// Tetrahedron vertex indices, 4 per tet. Indices reference
    /// [`Self::vertices`]. After [`Self::new`], every tet has positive
    /// signed volume (clockwise tets were silently re-oriented by
    /// swapping `v₂ ↔ v₃`).
    pub tetrahedra: Vec<[usize; 4]>,
    /// Per-vertex material tag. Length matches [`Self::vertices`]. Used
    /// by the FEM eigensolver to flag vertices sitting on material
    /// interfaces / PEC boundaries.
    pub vertex_material: Vec<MaterialTag>,
    /// Per-tet material tag. Length matches [`Self::tetrahedra`]. The
    /// FEM eigensolver looks up `eps_r` / `mu_r` from this tag during
    /// element-matrix assembly.
    pub tetrahedron_material: Vec<MaterialTag>,
}

/// Minimum signed volume below which a tet is treated as degenerate.
///
/// Chosen at `1e-18` so that double-precision round-off on cavities
/// with extents of order `1e-2 m` (the WR-90 cavity is 22.86 × 10.16 ×
/// 30 mm) cannot trigger a false rejection: a single Kuhn tet in that
/// cavity has volume of order `1e-7 m³`, eleven orders of magnitude
/// above the threshold.
const MIN_SIGNED_VOLUME: f64 = 1.0e-18;

impl TetMesh3D {
    /// Build a `TetMesh3D` after validating its invariants.
    ///
    /// See the type-level docs for the full invariant list. The most
    /// important guarantee for downstream consumers (`yee-fem`
    /// assembly) is that **every tet in the returned mesh has
    /// strictly positive signed volume** — even if the caller handed
    /// us a tet with `(v₀, v₁, v₃, v₂)` ordering, the swap is performed
    /// silently and the caller does not need to worry about
    /// orientation.
    ///
    /// Material tags default to `0` for every vertex and every tet if
    /// the caller passes `None`. Pass `Some(vec)` to supply explicit
    /// tags; lengths must match the corresponding primary array.
    ///
    /// Returns [`Error::Invalid`] with a descriptive message on any
    /// invariant violation that cannot be silently fixed (degenerate
    /// volume, out-of-range index, material-length mismatch).
    pub fn new(
        vertices: Vec<Vector3<f64>>,
        tetrahedra: Vec<[usize; 4]>,
        vertex_material: Option<Vec<MaterialTag>>,
        tetrahedron_material: Option<Vec<MaterialTag>>,
    ) -> Result<Self> {
        if vertices.len() < 4 {
            return Err(Error::Invalid(format!(
                "TetMesh3D requires >= 4 vertices, got {}",
                vertices.len()
            )));
        }
        if tetrahedra.is_empty() {
            return Err(Error::Invalid(
                "TetMesh3D requires >= 1 tetrahedron".to_string(),
            ));
        }

        // Validate index bounds first so the signed-volume computation
        // below can assume valid indices.
        for (i, tet) in tetrahedra.iter().enumerate() {
            for &idx in tet {
                if idx >= vertices.len() {
                    return Err(Error::Invalid(format!(
                        "tetrahedron {i} references vertex index {idx} >= vertex count {}",
                        vertices.len()
                    )));
                }
            }
        }

        // Auto-reorient any tet whose signed volume is negative; reject
        // any tet whose signed volume is too small or non-finite to
        // recover (degenerate / collinear vertices).
        let mut tetrahedra = tetrahedra;
        for (i, tet) in tetrahedra.iter_mut().enumerate() {
            let v = signed_volume_from_indices(&vertices, tet);
            if !v.is_finite() {
                return Err(Error::Invalid(format!(
                    "tetrahedron {i} has non-finite signed volume (degenerate vertices)"
                )));
            }
            if v.abs() < MIN_SIGNED_VOLUME {
                return Err(Error::Invalid(format!(
                    "tetrahedron {i} has near-zero signed volume {v:.3e} (degenerate or coplanar vertices)"
                )));
            }
            if v < 0.0 {
                // Silently re-orient by swapping v2 <-> v3. The
                // resulting signed volume is the negation of the
                // original, which is now strictly positive.
                tet.swap(2, 3);
            }
        }

        let vertex_material = match vertex_material {
            Some(v) => {
                if v.len() != vertices.len() {
                    return Err(Error::Invalid(format!(
                        "vertex_material ({}) must match vertices ({})",
                        v.len(),
                        vertices.len()
                    )));
                }
                v
            }
            None => vec![0; vertices.len()],
        };
        let tetrahedron_material = match tetrahedron_material {
            Some(v) => {
                if v.len() != tetrahedra.len() {
                    return Err(Error::Invalid(format!(
                        "tetrahedron_material ({}) must match tetrahedra ({})",
                        v.len(),
                        tetrahedra.len()
                    )));
                }
                v
            }
            None => vec![0; tetrahedra.len()],
        };

        Ok(Self {
            vertices,
            tetrahedra,
            vertex_material,
            tetrahedron_material,
        })
    }

    /// Number of tetrahedra.
    pub fn n_tets(&self) -> usize {
        self.tetrahedra.len()
    }

    /// Number of vertices.
    pub fn n_verts(&self) -> usize {
        self.vertices.len()
    }

    /// Signed volume of tetrahedron `i`, computed as
    /// `(1/6) · (v₁ − v₀) · ((v₂ − v₀) × (v₃ − v₀))`. After
    /// [`Self::new`] every tet has positive signed volume, so this
    /// returns a positive number for any `i < self.n_tets()`.
    ///
    /// # Panics
    /// Panics if `i >= self.n_tets()`. Callers driven by the
    /// eigensolver iterate over `0..n_tets()` so the bounds check is
    /// upheld by construction.
    pub fn signed_volume(&self, i: usize) -> f64 {
        signed_volume_from_indices(&self.vertices, &self.tetrahedra[i])
    }

    /// Centroid of tetrahedron `i` — arithmetic mean of its four
    /// vertices.
    ///
    /// # Panics
    /// Panics if `i >= self.n_tets()`. Same rationale as
    /// [`Self::signed_volume`].
    pub fn centroid(&self, i: usize) -> Vector3<f64> {
        let tet = &self.tetrahedra[i];
        let v0 = &self.vertices[tet[0]];
        let v1 = &self.vertices[tet[1]];
        let v2 = &self.vertices[tet[2]];
        let v3 = &self.vertices[tet[3]];
        (v0 + v1 + v2 + v3) / 4.0
    }
}

/// Compute the signed volume `(1/6) · (v₁ − v₀) · ((v₂ − v₀) × (v₃ − v₀))`
/// of the tet given by the four vertex indices into `vertices`.
///
/// Free function (not a method) so the constructor can call it during
/// validation, before `self` exists. Indices are assumed to be in
/// bounds; the caller validates that.
fn signed_volume_from_indices(vertices: &[Vector3<f64>], tet: &[usize; 4]) -> f64 {
    let v0 = &vertices[tet[0]];
    let v1 = &vertices[tet[1]];
    let v2 = &vertices[tet[2]];
    let v3 = &vertices[tet[3]];
    let e1 = v1 - v0;
    let e2 = v2 - v0;
    let e3 = v3 - v0;
    e1.dot(&e2.cross(&e3)) / 6.0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference unit tet `[(0,0,0), (1,0,0), (0,1,0), (0,0,1)]` has
    /// known signed volume `1/6`.
    fn reference_unit_tet() -> TetMesh3D {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let tetrahedra = vec![[0, 1, 2, 3]];
        TetMesh3D::new(vertices, tetrahedra, None, None).unwrap()
    }

    #[test]
    fn construct_reference_tet_volume_one_sixth() {
        let m = reference_unit_tet();
        assert_eq!(m.n_verts(), 4);
        assert_eq!(m.n_tets(), 1);
        assert!((m.signed_volume(0) - 1.0 / 6.0).abs() < 1e-15);
        // Default material tags are zero.
        assert_eq!(m.vertex_material, vec![0u32; 4]);
        assert_eq!(m.tetrahedron_material, vec![0u32; 1]);
    }

    #[test]
    fn construct_enforces_positive_signed_volume_via_reorientation() {
        // Same geometry as the reference tet but with v2 and v3 swapped,
        // giving a negative signed volume. `new` must silently swap
        // them back so the stored tet has positive signed volume.
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0), // v3 in canonical order
            Vector3::new(0.0, 1.0, 0.0), // v2 in canonical order
        ];
        let tetrahedra = vec![[0, 1, 2, 3]];
        let m = TetMesh3D::new(vertices, tetrahedra, None, None).unwrap();
        assert!(
            m.signed_volume(0) > 0.0,
            "post-construction signed volume must be positive, got {}",
            m.signed_volume(0)
        );
        assert!((m.signed_volume(0) - 1.0 / 6.0).abs() < 1e-15);
        // The reorientation swaps positions 2 and 3 of the stored tet.
        assert_eq!(m.tetrahedra[0], [0, 1, 3, 2]);
    }

    #[test]
    fn out_of_range_tet_index_rejected() {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let tetrahedra = vec![[0, 1, 2, 9]];
        let err = TetMesh3D::new(vertices, tetrahedra, None, None).unwrap_err();
        match err {
            Error::Invalid(msg) => assert!(msg.contains("vertex index 9"), "got: {msg}"),
            _ => panic!("expected Invalid"),
        }
    }

    #[test]
    fn degenerate_coplanar_tet_rejected() {
        // Four coplanar points (all at z = 0) → signed volume 0 → rejected.
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
        ];
        let tetrahedra = vec![[0, 1, 2, 3]];
        let err = TetMesh3D::new(vertices, tetrahedra, None, None).unwrap_err();
        match err {
            Error::Invalid(msg) => {
                assert!(msg.contains("near-zero signed volume"), "got: {msg}");
            }
            _ => panic!("expected Invalid"),
        }
    }

    #[test]
    fn too_few_vertices_rejected() {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let tetrahedra = vec![[0, 1, 2, 0]]; // degenerate but caught earlier
        let err = TetMesh3D::new(vertices, tetrahedra, None, None).unwrap_err();
        match err {
            Error::Invalid(msg) => assert!(msg.contains(">= 4 vertices"), "got: {msg}"),
            _ => panic!("expected Invalid"),
        }
    }

    #[test]
    fn material_length_mismatch_rejected() {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let tetrahedra = vec![[0, 1, 2, 3]];
        let err = TetMesh3D::new(
            vertices,
            tetrahedra,
            Some(vec![1, 2]), // length 2, but 4 vertices
            None,
        )
        .unwrap_err();
        match err {
            Error::Invalid(msg) => assert!(msg.contains("vertex_material"), "got: {msg}"),
            _ => panic!("expected Invalid"),
        }
    }

    #[test]
    fn explicit_material_tags_preserved() {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
        ];
        // Two tets sharing face (0, 1, 2). The second tet has v4 above the plane.
        let tetrahedra = vec![[0, 1, 2, 3], [0, 1, 2, 4]];
        let m = TetMesh3D::new(
            vertices,
            tetrahedra,
            Some(vec![1, 2, 3, 4, 5]),
            Some(vec![10, 20]),
        )
        .unwrap();
        assert_eq!(m.vertex_material, vec![1, 2, 3, 4, 5]);
        assert_eq!(m.tetrahedron_material, vec![10, 20]);
        // Both stored tets have positive signed volume (reorientation
        // handled silently for whichever orientation the caller used).
        assert!(m.signed_volume(0) > 0.0);
        assert!(m.signed_volume(1) > 0.0);
    }

    #[test]
    fn centroid_of_reference_tet_at_known_position() {
        // Centroid of (0,0,0), (1,0,0), (0,1,0), (0,0,1) is (1/4, 1/4, 1/4).
        let m = reference_unit_tet();
        let c = m.centroid(0);
        assert!((c.x - 0.25).abs() < 1e-15);
        assert!((c.y - 0.25).abs() < 1e-15);
        assert!((c.z - 0.25).abs() < 1e-15);
    }

    /// Five-tet decomposition of the unit cube `[0,1]³`. This is the
    /// well-known canonical 5-tet split that exhibits 4 corner tets +
    /// 1 central tet; the total volume is `1` and each tet has
    /// volume `1/6` or `1/3`. Used by the cube-volume-sum test below.
    fn unit_cube_five_tet() -> TetMesh3D {
        // Cube vertex layout (z = 0 bottom face, z = 1 top face):
        //   0 = (0,0,0)   1 = (1,0,0)   2 = (1,1,0)   3 = (0,1,0)
        //   4 = (0,0,1)   5 = (1,0,1)   6 = (1,1,1)   7 = (0,1,1)
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
        ];
        // Standard 5-tet split of a cube (one central tet + four
        // corner tets). The central tet is (1, 4, 6, 3); the four
        // corner tets share one face with it each. Orientation is
        // handled silently by `new`, so the absolute volumes here add
        // to 1.
        let tetrahedra = vec![
            [0, 1, 3, 4], // corner near v0
            [1, 2, 3, 6], // corner near v2
            [1, 5, 4, 6], // corner near v5
            [3, 4, 7, 6], // corner near v7
            [1, 3, 4, 6], // central
        ];
        TetMesh3D::new(vertices, tetrahedra, None, None).unwrap()
    }

    #[test]
    fn unit_cube_five_tet_total_volume_is_one() {
        let m = unit_cube_five_tet();
        assert_eq!(m.n_tets(), 5);
        let total: f64 = (0..m.n_tets()).map(|i| m.signed_volume(i)).sum();
        assert!(
            (total - 1.0).abs() < 1e-12,
            "five-tet cube total volume = {total}, expected 1.0 within 1e-12"
        );
        // Every stored tet has positive signed volume after construction.
        for i in 0..m.n_tets() {
            assert!(m.signed_volume(i) > 0.0, "tet {i} has non-positive volume");
        }
    }
}
