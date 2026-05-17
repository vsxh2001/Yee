//! `TriMesh2D` adaptor: edge enumeration, boundary detection, and
//! per-triangle local→global edge orientation map.
//!
//! The Nedelec edge basis requires a **canonical orientation per global
//! edge**, with each triangle reporting whether its local edge direction
//! agrees (`+1`) or opposes (`−1`) the canonical direction. Without this
//! sign convention the assembled curl-curl block has wrong off-diagonal
//! entries and the dominant eigenvalue is meaningless. This module owns
//! the bookkeeping.

use std::collections::HashMap;
use yee_mesh::TriMesh2D;

/// Canonical orientation of a global edge: `from < to` in vertex-index
/// order. A triangle's local edge `(a, b)` with `a < b` matches the
/// canonical orientation (sign `+1`); `a > b` opposes it (sign `−1`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct EdgeKey {
    /// Smaller vertex index.
    pub from: usize,
    /// Larger vertex index.
    pub to: usize,
}

impl EdgeKey {
    /// Build a canonical key from an unordered pair `(a, b)`.
    fn new(a: usize, b: usize) -> Self {
        if a < b {
            Self { from: a, to: b }
        } else {
            Self { from: b, to: a }
        }
    }
}

/// Per-triangle edge connectivity into the global edge table.
///
/// Local edge numbering follows Jin §8.5 convention: edge `e` lies
/// opposite local vertex `e`, traversed in CCW order around the triangle.
/// For a CCW triangle with local vertices `(v0, v1, v2)`:
///
/// * local edge 0 connects `v1 → v2` (opposite v0)
/// * local edge 1 connects `v2 → v0` (opposite v1)
/// * local edge 2 connects `v0 → v1` (opposite v2)
///
/// `global_edge[e]` is the global edge index; `sign[e]` is `+1` if the
/// local direction (above) matches the canonical `from < to` orientation
/// of that global edge, else `−1`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct TriEdgeConnectivity {
    /// Global edge index for local edges 0, 1, 2.
    pub global_edge: [usize; 3],
    /// Orientation sign for each local edge (`+1` or `−1`).
    pub sign: [f64; 3],
}

/// Edge table for a [`TriMesh2D`].
///
/// Holds the global edge list, per-edge boundary flag, and the
/// per-triangle local→global edge map with orientation signs. Built
/// once per mesh by [`EdgeTable::build`] and consumed read-only by
/// [`super::assembly`].
#[derive(Debug, Clone)]
pub(crate) struct EdgeTable {
    /// All distinct edges (canonical orientation `from < to`).
    pub edges: Vec<EdgeKey>,
    /// `true` iff the edge is on the mesh boundary (lies in exactly
    /// one triangle). PEC walls in the eigensolver are tagged via
    /// these boundary edges.
    pub is_boundary: Vec<bool>,
    /// Per-triangle local→global edge map and orientation.
    pub tri_edges: Vec<TriEdgeConnectivity>,
}

impl EdgeTable {
    /// Walk the mesh and build the edge table. `O(n_tris)` with a
    /// `HashMap` for canonical-key dedup.
    pub fn build(mesh: &TriMesh2D) -> Self {
        let mut edge_map: HashMap<EdgeKey, usize> = HashMap::new();
        let mut edges: Vec<EdgeKey> = Vec::new();
        // edge_count[i] is the number of triangles incident on global
        // edge i; the boundary predicate is `count == 1`.
        let mut edge_count: Vec<usize> = Vec::new();
        let mut tri_edges: Vec<TriEdgeConnectivity> = Vec::with_capacity(mesh.n_tris());

        // Local edge endpoints per Jin §8.5 (edge e opposite local vertex e,
        // traversed CCW around the triangle).
        let local_edges: [[usize; 2]; 3] = [[1, 2], [2, 0], [0, 1]];

        for tri in &mesh.triangles {
            let mut global_edge = [0usize; 3];
            let mut sign = [0.0f64; 3];
            for (e, &[la, lb]) in local_edges.iter().enumerate() {
                let a = tri[la];
                let b = tri[lb];
                let key = EdgeKey::new(a, b);
                let idx = *edge_map.entry(key).or_insert_with(|| {
                    edges.push(key);
                    edge_count.push(0);
                    edges.len() - 1
                });
                edge_count[idx] += 1;
                global_edge[e] = idx;
                // Local direction is `a -> b`. Canonical is `from -> to`
                // with `from < to`. Sign is `+1` iff `a < b`.
                sign[e] = if a < b { 1.0 } else { -1.0 };
            }
            tri_edges.push(TriEdgeConnectivity { global_edge, sign });
        }

        let is_boundary: Vec<bool> = edge_count.iter().map(|&c| c == 1).collect();

        Self {
            edges,
            is_boundary,
            tri_edges,
        }
    }

    /// Number of global edges (boundary + interior).
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }

    /// Number of interior edges — these are the DoFs of the Nedelec
    /// transverse-`E_t` block after PEC Dirichlet boundary elimination.
    /// Used by unit tests (`n_interior_edges == 1` for the canonical
    /// two-triangle unit-square fixture).
    #[allow(dead_code)] // unit-test helper; not used by the assembly path
    pub fn n_interior_edges(&self) -> usize {
        self.is_boundary.iter().filter(|&&b| !b).count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_tri_unit_square() -> TriMesh2D {
        // Two CCW triangles sharing the diagonal v0-v2.
        // v3 (0,1) ── v2 (1,1)
        //  │ ╲          │
        //  │   ╲   tri1 │
        //  │     ╲      │
        //  │ tri0  ╲    │
        // v0 (0,0) ── v1 (1,0)
        TriMesh2D::new(
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            vec![[0, 1, 2], [0, 2, 3]],
            None,
            None,
        )
        .unwrap()
    }

    #[test]
    fn edge_enumeration_two_tri_mesh() {
        // Two triangles sharing one diagonal → 4 boundary + 1 interior = 5 edges.
        let mesh = two_tri_unit_square();
        let table = EdgeTable::build(&mesh);
        assert_eq!(table.n_edges(), 5);
        assert_eq!(table.n_interior_edges(), 1);
        // The shared diagonal connects vertex 0 and vertex 2.
        let diag = EdgeKey::new(0, 2);
        let diag_idx = table.edges.iter().position(|&e| e == diag).unwrap();
        assert!(!table.is_boundary[diag_idx]);
    }

    #[test]
    fn orientation_signs_consistent_on_shared_edge() {
        // The shared edge is traversed in opposite local directions by
        // the two triangles, so their orientation signs must differ.
        let mesh = two_tri_unit_square();
        let table = EdgeTable::build(&mesh);
        // tri0 = [0,1,2]: local edges (1->2, 2->0, 0->1)
        //   edge "2->0" is the diagonal, traversed as 2→0, canonical 0→2,
        //   so sign = -1.
        // tri1 = [0,2,3]: local edges (2->3, 3->0, 0->2)
        //   edge "0->2" is the diagonal, traversed as 0→2, canonical 0→2,
        //   so sign = +1.
        let tri0 = &table.tri_edges[0];
        let tri1 = &table.tri_edges[1];
        let diag = EdgeKey::new(0, 2);
        let diag_idx = table.edges.iter().position(|&e| e == diag).unwrap();
        // Find which local edge of each triangle maps to the diagonal.
        let local0 = tri0
            .global_edge
            .iter()
            .position(|&g| g == diag_idx)
            .unwrap();
        let local1 = tri1
            .global_edge
            .iter()
            .position(|&g| g == diag_idx)
            .unwrap();
        // Signs on the same global edge must be opposites — that is the
        // condition that makes the Nedelec curl assembly come out
        // single-valued.
        assert!((tri0.sign[local0] * tri1.sign[local1] + 1.0).abs() < 1e-15);
    }

    #[test]
    fn boundary_flag_on_perimeter() {
        let mesh = two_tri_unit_square();
        let table = EdgeTable::build(&mesh);
        // Every perimeter edge (4 of them) must be flagged boundary.
        let perimeter = [
            EdgeKey::new(0, 1),
            EdgeKey::new(1, 2),
            EdgeKey::new(2, 3),
            EdgeKey::new(0, 3),
        ];
        for e in perimeter {
            let idx = table.edges.iter().position(|&x| x == e).unwrap();
            assert!(table.is_boundary[idx], "edge {e:?} should be boundary");
        }
    }
}
