//! RWG basis on triangle meshes. Phase 1.0 free-space MoM.
//!
//! This module implements the Rao–Wilton–Glisson (RWG) vector basis functions
//! used by the planar / free-space MoM solver. Each interior edge of the input
//! [`TriMesh`] becomes one RWG basis function whose support is the pair of
//! triangles sharing that edge.
//!
//! Reference: S. M. Rao, D. R. Wilton, A. W. Glisson, "Electromagnetic
//! scattering by surfaces of arbitrary shape," *IEEE Trans. Antennas Propag.*,
//! vol. 30, no. 3, pp. 409–418, May 1982.
//!
//! On `tri_plus`, the basis evaluates to
//! `+ length / (2 * area_plus)  * (r - p_free_plus)` and its surface
//! divergence to `+ length / area_plus`. On `tri_minus`, both quantities
//! flip sign. Outside the two-triangle support the basis is identically zero.
//!
//! Dead-code allowances: this module lands in Phase 1.0 Task 3+4 ahead of
//! the impedance fill that consumes it (Tasks 5–8). Clippy with
//! `-D warnings` flags every yet-unused symbol as an error, so the struct
//! and its associated items are explicitly tagged `#[allow(dead_code)]`
//! at the module boundary. The allow will be removed implicitly once
//! `fill.rs` references these items.
#![allow(dead_code)]

use nalgebra::Vector3;
use std::collections::BTreeMap;
use yee_mesh::TriMesh;

/// One Rao–Wilton–Glisson edge in the triangulation.
///
/// `v0`, `v1` are sorted (`v0 < v1`) so the same shared edge is keyed exactly
/// once when enumerating. `tri_plus`/`tri_minus` index into [`TriMesh::triangles`];
/// `free_plus`/`free_minus` are the *vertex* indices (into
/// [`TriMesh::vertices`]) of the free (non-shared) vertex on each adjacent
/// triangle. `port_tag` is non-zero only when the two adjacent triangles
/// carry DIFFERENT non-zero mesh tags (the edge sits on the boundary
/// between two tagged regions, the delta-gap port convention).
///
/// Note: `tri_plus` is the first triangle in mesh-iteration order that
/// touches this shared edge — it carries no geometric meaning. The
/// `eval`/`div` sign conventions ensure the basis function is well-defined
/// regardless of which side is called `+`.
// NOTE: visibility is `pub` (not `pub(crate)`) only to allow the
// doc-hidden `lib::__internal` test-helper surface to re-export `RwgEdge`
// and `RwgBasis` for diagnostic integration tests. The type is **not** part
// of the stable API — depending on it from downstream crates is unsupported.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RwgEdge {
    /// Lower-indexed shared-edge vertex (vertex index into the mesh).
    pub v0: u32,
    /// Higher-indexed shared-edge vertex (vertex index into the mesh).
    pub v1: u32,
    /// Triangle index of the "+" adjacent triangle. See struct-level note:
    /// the `+`/`-` labeling is enumeration-order, not geometric.
    pub tri_plus: u32,
    /// Triangle index of the "-" adjacent triangle. See struct-level note:
    /// the `+`/`-` labeling is enumeration-order, not geometric.
    pub tri_minus: u32,
    /// Free (non-shared) vertex index of the "+" triangle.
    pub free_plus: u32,
    /// Free (non-shared) vertex index of the "-" triangle.
    pub free_minus: u32,
    /// Edge length (meters).
    pub length: f64,
    /// Non-zero when the two adjacent triangles carry DIFFERENT non-zero
    /// tags (the edge straddles the boundary between two tagged regions —
    /// convention for delta-gap port placement). Set to `0` otherwise.
    pub port_tag: u32,
}

/// Geometric and topological summary of an [`RwgEdge`] enumeration on a
/// [`TriMesh`]: per-triangle centroids/normals/areas plus the deduplicated
/// edge list.
// NOTE: visibility is `pub` only for the doc-hidden `__internal` test-helper
// surface. See the `RwgEdge` comment for the stability stance.
#[derive(Debug, Clone)]
pub struct RwgBasis {
    /// The underlying mesh. Owned so basis indices remain stable.
    pub(crate) mesh: TriMesh,
    /// Deduplicated interior edges. One entry per RWG basis function.
    pub edges: Vec<RwgEdge>,
    /// Triangle centroids `(a + b + c) / 3`.
    pub(crate) centroids: Vec<Vector3<f64>>,
    /// Triangle outward normals (right-hand rule from the vertex winding).
    pub(crate) normals: Vec<Vector3<f64>>,
    /// Triangle areas (meters squared). Always strictly positive.
    pub(crate) areas: Vec<f64>,
}

impl RwgBasis {
    /// Build an [`RwgBasis`] by computing per-triangle geometry and
    /// enumerating the interior edges of `mesh`.
    ///
    /// Returns [`yee_core::Error::Invalid`] in two cases:
    /// 1. Any triangle has non-positive area (degenerate / collinear).
    /// 2. Any edge is shared by 3 or more triangles (non-manifold mesh).
    ///
    /// Boundary edges (touched by exactly one triangle) are silently dropped
    /// — they do not carry an RWG basis function.
    pub fn from_mesh(mesh: TriMesh) -> Result<Self, yee_core::Error> {
        let n_tris = mesh.n_tris();

        let mut centroids = Vec::with_capacity(n_tris);
        let mut normals = Vec::with_capacity(n_tris);
        let mut areas = Vec::with_capacity(n_tris);

        for (t, tri) in mesh.triangles.iter().enumerate() {
            let a = mesh.vertices[tri[0] as usize];
            let b = mesh.vertices[tri[1] as usize];
            let c = mesh.vertices[tri[2] as usize];
            let centroid = (a + b + c) / 3.0;
            let cross = (b - a).cross(&(c - a));
            let cross_norm = cross.norm();
            let area = 0.5 * cross_norm;
            // Non-positive (or NaN) area means the triangle is degenerate.
            // The basis functions divide by area, so refusing degenerate
            // input here prevents NaN/Inf propagation downstream. We check
            // via `partial_cmp` to robustly catch NaN — `area <= 0.0` would
            // miss it because every NaN comparison is `false`.
            if !matches!(
                area.partial_cmp(&0.0),
                Some(std::cmp::Ordering::Greater)
            ) {
                return Err(yee_core::Error::Invalid(format!(
                    "triangle {t} has non-positive area {area}"
                )));
            }
            centroids.push(centroid);
            normals.push(cross / cross_norm);
            areas.push(area);
        }

        // Edge key: sorted (v_lo, v_hi). Value: list of (triangle_idx, free_vertex_idx).
        //
        // `BTreeMap` (not `HashMap`) so the resulting `edges` vector is
        // ordered by sorted edge key and is therefore *deterministic*
        // across runs and machines. Basis-function indices are derived
        // from this iteration order, and downstream matrix-entry
        // assertions need a stable mapping.
        let mut edge_map: BTreeMap<(u32, u32), Vec<(u32, u32)>> = BTreeMap::new();
        for (t_idx, tri) in mesh.triangles.iter().enumerate() {
            let t = t_idx as u32;
            // Each of the three edges of the triangle is keyed by sorted
            // (v_lo, v_hi); the remaining vertex is the free vertex for that
            // edge on this triangle.
            for &(i, j, k) in &[(0usize, 1usize, 2usize), (1, 2, 0), (2, 0, 1)] {
                let vi = tri[i];
                let vj = tri[j];
                let vk = tri[k];
                let key = if vi < vj { (vi, vj) } else { (vj, vi) };
                edge_map.entry(key).or_default().push((t, vk));
            }
        }

        let mut edges = Vec::new();
        for ((v_lo, v_hi), adj) in edge_map {
            match adj.len() {
                1 => {
                    // Boundary edge — no RWG basis function lives here.
                }
                2 => {
                    let (tri_plus, free_plus) = adj[0];
                    let (tri_minus, free_minus) = adj[1];
                    let p0 = mesh.vertices[v_lo as usize];
                    let p1 = mesh.vertices[v_hi as usize];
                    let length = (p0 - p1).norm();
                    let tag_plus = mesh.tags[tri_plus as usize];
                    let tag_minus = mesh.tags[tri_minus as usize];
                    // Delta-gap port convention: an edge is a port edge iff
                    // its two adjacent triangles carry DIFFERENT non-zero
                    // tags (the edge straddles the boundary between two
                    // tagged regions). The earlier "same non-zero tag" rule
                    // mis-tagged every interior edge of a uniformly tagged
                    // region as a port edge — see the dipole-gate diagnosis
                    // in the Phase 1.0 Task 11 report. Both-zero, one-zero,
                    // and equal non-zero pairs all yield `port_tag = 0`.
                    let port_tag = match (tag_plus, tag_minus) {
                        (a, b) if a != 0 && b != 0 && a != b => a.min(b),
                        _ => 0,
                    };
                    edges.push(RwgEdge {
                        v0: v_lo,
                        v1: v_hi,
                        tri_plus,
                        tri_minus,
                        free_plus,
                        free_minus,
                        length,
                        port_tag,
                    });
                }
                n => {
                    return Err(yee_core::Error::Invalid(format!(
                        "non-manifold edge ({v_lo}, {v_hi}) shared by {n} triangles"
                    )));
                }
            }
        }

        Ok(Self {
            mesh,
            edges,
            centroids,
            normals,
            areas,
        })
    }

    /// Number of RWG basis functions (one per interior edge).
    pub fn n_basis(&self) -> usize {
        self.edges.len()
    }

    /// Number of triangles in the underlying mesh.
    pub(crate) fn n_tris(&self) -> usize {
        self.mesh.n_tris()
    }

    /// Iterate the basis-function indices whose edge carries the given
    /// non-zero `port_tag`. With `port_tag == 0` no basis is selected (the
    /// "no port" sentinel is intentionally not iterable).
    pub fn port_basis_indices(
        &self,
        port_tag: u32,
    ) -> impl Iterator<Item = usize> + '_ {
        // `port_tag == 0` means "not a port" inside `RwgEdge`, so iterating
        // it would spuriously return every boundary-mismatched edge — guard
        // it out instead of trusting the caller.
        self.edges
            .iter()
            .enumerate()
            .filter(move |(_, e)| port_tag != 0 && e.port_tag == port_tag)
            .map(|(i, _)| i)
    }

    /// Evaluate the `k`-th RWG basis function at barycentric coordinates
    /// `bary` on triangle `tri`.
    ///
    /// Returns `Vector3::zeros()` if `tri` is neither `tri_plus` nor
    /// `tri_minus` for basis `k` (i.e. `tri` lies outside the basis
    /// function's two-triangle support). This includes sentinel/invalid
    /// triangle indices.
    ///
    /// # Panics
    ///
    /// Panics if `k >= self.n_basis()`.
    pub(crate) fn eval(&self, k: usize, tri: u32, bary: [f64; 3]) -> Vector3<f64> {
        let edge = &self.edges[k];
        // Reconstruct the spatial point from barycentric coordinates of the
        // triangle that was supplied.
        let reconstruct_r = |t: u32| -> Vector3<f64> {
            let tri_vs = self.mesh.triangles[t as usize];
            bary[0] * self.mesh.vertices[tri_vs[0] as usize]
                + bary[1] * self.mesh.vertices[tri_vs[1] as usize]
                + bary[2] * self.mesh.vertices[tri_vs[2] as usize]
        };
        if tri == edge.tri_plus {
            let r = reconstruct_r(tri);
            let p = self.mesh.vertices[edge.free_plus as usize];
            let scale = edge.length / (2.0 * self.areas[tri as usize]);
            scale * (r - p)
        } else if tri == edge.tri_minus {
            let r = reconstruct_r(tri);
            let p = self.mesh.vertices[edge.free_minus as usize];
            let scale = edge.length / (2.0 * self.areas[tri as usize]);
            -scale * (r - p)
        } else {
            Vector3::zeros()
        }
    }

    /// Surface divergence of the `k`-th RWG basis function on triangle `tri`.
    ///
    /// The RWG divergence is piecewise constant: `+length/area_plus` on
    /// `tri_plus`, `-length/area_minus` on `tri_minus`, and `0.0` elsewhere.
    /// Returns `0.0` if `tri` is neither `tri_plus` nor `tri_minus` for
    /// basis `k` (i.e. `tri` lies outside the basis function's
    /// two-triangle support). This includes sentinel/invalid triangle
    /// indices.
    ///
    /// # Panics
    ///
    /// Panics if `k >= self.n_basis()`.
    pub(crate) fn div(&self, k: usize, tri: u32) -> f64 {
        let edge = &self.edges[k];
        if tri == edge.tri_plus {
            edge.length / self.areas[edge.tri_plus as usize]
        } else if tri == edge.tri_minus {
            -edge.length / self.areas[edge.tri_minus as usize]
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    /// 2-triangle unit-square mesh used by several tests below.
    fn two_tri_mesh() -> TriMesh {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let triangles = vec![[0u32, 1, 2], [0u32, 2, 3]];
        let tags = vec![0u32, 0u32];
        TriMesh::new(vertices, triangles, tags).unwrap()
    }

    #[test]
    fn edge_count_two_tri_mesh() {
        // Two triangles sharing one diagonal edge: 5 mesh edges total,
        // 4 boundary + 1 interior, hence exactly one RWG basis function.
        let basis = RwgBasis::from_mesh(two_tri_mesh()).expect("valid mesh");
        assert_eq!(basis.n_basis(), 1, "expected exactly one RWG edge");
        assert_eq!(basis.n_tris(), 2);
        assert_eq!(basis.areas.len(), 2);
        assert_eq!(basis.centroids.len(), 2);
        assert_eq!(basis.normals.len(), 2);
        // Each half of the unit square has area 1/2.
        for a in &basis.areas {
            assert!((a - 0.5).abs() < 1e-12);
        }
        // The shared edge is the diagonal (0, 2) of length sqrt(2).
        let edge = basis.edges[0];
        assert_eq!((edge.v0, edge.v1), (0, 2));
        assert!((edge.length - 2.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn rejects_non_positive_area() {
        // Three collinear vertices give a degenerate triangle with zero
        // area — `from_mesh` must reject it rather than divide by zero
        // downstream.
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
        ];
        let triangles = vec![[0u32, 1, 2]];
        let tags = vec![0u32];
        let mesh = TriMesh::new(vertices, triangles, tags).unwrap();
        let err = RwgBasis::from_mesh(mesh).expect_err("collinear triangle must be rejected");
        match err {
            yee_core::Error::Invalid(_) => {}
            other => panic!("expected Error::Invalid, got {other:?}"),
        }
    }

    #[test]
    fn divergence_sign_and_magnitude() {
        // div on tri_plus is +length/area_plus, on tri_minus -length/area_minus.
        // For the symmetric 2-triangle square (equal areas) the two divergences
        // are exact opposites and sum to zero.
        let basis = RwgBasis::from_mesh(two_tri_mesh()).expect("valid mesh");
        let edge = basis.edges[0];
        let dp = basis.div(0, edge.tri_plus);
        let dm = basis.div(0, edge.tri_minus);
        assert!(dp > 0.0, "div on tri_plus must be positive, got {dp}");
        assert!(dm < 0.0, "div on tri_minus must be negative, got {dm}");
        // Robust same-magnitude check: opposite signs and (here) equal
        // magnitudes — both halves of the square have the same area.
        assert!(
            (dp + dm).abs() < 1e-12,
            "equal-area divergences must cancel, got {dp} + {dm}"
        );
        // Off-support divergence is identically zero.
        assert_eq!(basis.div(0, 99), 0.0);
    }

    #[test]
    fn eval_zero_outside_support() {
        // A triangle index that does not appear in the basis's two-triangle
        // support must yield the zero vector regardless of barycentrics.
        let basis = RwgBasis::from_mesh(two_tri_mesh()).expect("valid mesh");
        let v = basis.eval(0, 99, [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]);
        assert_eq!(v, Vector3::zeros());
    }

    #[test]
    fn eval_vanishes_at_free_vertex() {
        // RWG basis is (length / 2A) * (r - p_free). Evaluated at the free
        // vertex itself, r == p_free, so the field is exactly zero. We
        // realise that point by setting the barycentric coordinate of the
        // free vertex to 1 within the triangle that contains it.
        let basis = RwgBasis::from_mesh(two_tri_mesh()).expect("valid mesh");
        let edge = basis.edges[0];
        let tri = edge.tri_plus;
        let tri_vs = basis.mesh.triangles[tri as usize];
        let mut bary = [0.0f64; 3];
        let mut found = false;
        for (i, &vi) in tri_vs.iter().enumerate() {
            if vi == edge.free_plus {
                bary[i] = 1.0;
                found = true;
                break;
            }
        }
        assert!(
            found,
            "free_plus must be a vertex of tri_plus by construction"
        );
        let v = basis.eval(0, tri, bary);
        assert!(
            v.norm() < 1e-12,
            "RWG basis must vanish at the free vertex, got {v:?}"
        );
    }

    #[test]
    fn port_tag_marks_boundary_between_tagged_regions() {
        // 4-tri mesh: two adjacent cells with different non-zero tags.
        // Tags [1, 1, 2, 2] → edge between cell 1 (tag 1) and cell 2 (tag 2)
        // is the port edge; other edges have either same tag or untagged.
        let vertices = vec![
            nalgebra::Vector3::new(0.0, 0.0, 0.0),
            nalgebra::Vector3::new(1.0, 0.0, 0.0),
            nalgebra::Vector3::new(1.0, 1.0, 0.0),
            nalgebra::Vector3::new(0.0, 1.0, 0.0),
            nalgebra::Vector3::new(2.0, 0.0, 0.0),
            nalgebra::Vector3::new(2.0, 1.0, 0.0),
        ];
        let triangles = vec![
            [0u32, 1, 2],
            [0u32, 2, 3], // left cell, tag 1
            [1u32, 4, 5],
            [1u32, 5, 2], // right cell, tag 2
        ];
        let tags = vec![1u32, 1, 2, 2];
        let mesh = yee_mesh::TriMesh::new(vertices, triangles, tags).unwrap();
        let basis = RwgBasis::from_mesh(mesh).unwrap();
        let port_count = basis.port_basis_indices(1).count();
        assert!(
            port_count >= 1,
            "expected at least one boundary edge tagged as port"
        );
        // Within-cell edges (same-tag pair) should NOT be port.
        let within_tagged = basis
            .edges
            .iter()
            .filter(|e| {
                let t0 = basis.mesh.tags[e.tri_plus as usize];
                let t1 = basis.mesh.tags[e.tri_minus as usize];
                t0 == t1 && t0 != 0
            })
            .count();
        let port_with_same_tag = basis
            .edges
            .iter()
            .filter(|e| {
                let t0 = basis.mesh.tags[e.tri_plus as usize];
                let t1 = basis.mesh.tags[e.tri_minus as usize];
                t0 == t1 && t0 != 0 && e.port_tag != 0
            })
            .count();
        assert_eq!(
            port_with_same_tag, 0,
            "same-tag edges must not be port edges (found {port_with_same_tag} of {within_tagged})"
        );
    }
}
