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
use yee_mesh::{MaterialTag, TriMesh2D};

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

    /// Flag, per mesh vertex, whether the vertex lies on the PEC
    /// boundary (i.e. is an endpoint of at least one boundary edge).
    ///
    /// The mixed `(E_t, E_z)` formulation
    /// ([`super::assembly::assemble_mixed`]) imposes homogeneous
    /// Dirichlet on the longitudinal nodal `E_z` DoFs at PEC walls,
    /// exactly mirroring the boundary-edge elimination applied to the
    /// transverse `E_t` DoFs. A boundary vertex is therefore dropped
    /// from the interior-vertex DoF set the same way a boundary edge is
    /// dropped from the interior-edge DoF set.
    ///
    /// `n_verts` is the mesh's vertex count
    /// ([`yee_mesh::TriMesh2D::n_verts`]); the returned vector has that
    /// length with `true` marking boundary vertices.
    pub fn boundary_vertices(&self, n_verts: usize) -> Vec<bool> {
        let mut is_bnd = vec![false; n_verts];
        for (gid, edge) in self.edges.iter().enumerate() {
            if self.is_boundary[gid] {
                is_bnd[edge.from] = true;
                is_bnd[edge.to] = true;
            }
        }
        is_bnd
    }
}

/// Build the geometrically-graded `y`-grid lines for a horizontal-slab
/// cross-section spanning `[0, b]`, clustered toward the interior
/// dielectric interface `y = d1` (Phase 1.3.1.1 step 5.4).
///
/// The slab geometry stratifies the guide into a dielectric layer
/// `0 ≤ y ≤ d1` and an air layer `d1 ≤ y ≤ b`. The dominant inhomogeneous
/// mode concentrates a sharp field peak at the dielectric/air interface,
/// which uniform first-order Nedelec/nodal elements under-resolve. This
/// helper concentrates element rows there by giving each layer a
/// **geometric** node distribution whose cells shrink toward `d1`:
///
/// * Lower layer `[0, d1]`: `ny_lo` intervals of sizes `h, h/r, h/r², …`
///   from the PEC wall (`y = 0`) up to the interface — finest at `d1`.
/// * Upper layer `[d1, b]`: `ny_hi` intervals of sizes `…, h'/r², h'/r, h'`
///   from the interface up to the PEC wall (`y = b`) — finest at `d1`.
///
/// `ratio` is the geometric grading factor `r ≥ 1` (consecutive cells
/// differ by `r` toward the interface). `r = 1` recovers a uniform layer;
/// larger `r` clusters more aggressively. A node is placed **exactly** at
/// `y = d1` (it is the shared layer boundary), so a mesh built on these
/// lines keeps the dielectric/air material partition sharp — no element
/// straddles the interface.
///
/// Returns the strictly-increasing `y`-coordinates `0 = y₀ < y₁ < … < b`
/// (length `ny_lo + ny_hi + 1`), with `d1` present exactly once at index
/// `ny_lo`.
///
/// # Panics (debug)
/// Debug-asserts `0 < d1 < b`, `ny_lo ≥ 1`, `ny_hi ≥ 1`, and `ratio ≥ 1`.
#[allow(dead_code)] // step-5.4 convergence-study helper; consumed by mesh unit tests + the mirror in tests/eigensolver_inhomogeneous.rs
pub(crate) fn graded_y_lines(b: f64, d1: f64, ny_lo: usize, ny_hi: usize, ratio: f64) -> Vec<f64> {
    debug_assert!(
        d1 > 0.0 && d1 < b,
        "interface d1={d1} must lie in (0, b={b})"
    );
    debug_assert!(ny_lo >= 1 && ny_hi >= 1, "each layer needs ≥1 interval");
    debug_assert!(ratio >= 1.0, "grading ratio r={ratio} must be ≥ 1");

    // Geometric cell sizes within a layer of thickness `thickness` split
    // into `n` intervals, finest at the interface end. With ratio `r`, the
    // interval nearest the interface has size `s`, the next `s·r`, … so the
    // sizes (interface→wall) are `s·r^k`, summing to a geometric series
    // `s·(r^n − 1)/(r − 1) = thickness`. For `r = 1` the sizes are uniform.
    let layer_cell_sizes = |thickness: f64, n: usize| -> Vec<f64> {
        if (ratio - 1.0).abs() < 1e-12 {
            return vec![thickness / (n as f64); n];
        }
        let geom_sum: f64 = (0..n).map(|k| ratio.powi(k as i32)).sum();
        let s = thickness / geom_sum; // smallest cell, at the interface
        (0..n).map(|k| s * ratio.powi(k as i32)).collect()
    };

    let mut ys = Vec::with_capacity(ny_lo + ny_hi + 1);

    // Lower layer [0, d1]: walk from the wall (y=0) up to the interface.
    // Cell sizes interface→wall are `s·r^k`; from the wall they appear in
    // reverse (largest first), so the wall cell is the coarsest and the
    // interface cell the finest — exactly the desired clustering at d1.
    ys.push(0.0);
    let lo_sizes = layer_cell_sizes(d1, ny_lo); // index 0 = interface cell (finest)
    let mut y = 0.0;
    for k in (0..ny_lo).rev() {
        y += lo_sizes[k];
        ys.push(y);
    }
    // Snap the interface node to exactly d1 (kill geometric-sum roundoff so
    // the material-tag predicate `yc < d1` stays on the right side).
    let iface_idx = ys.len() - 1;
    ys[iface_idx] = d1;

    // Upper layer [d1, b]: walk from the interface up to the wall (y=b).
    // Interface→wall sizes are `s'·r^k`, so the first (interface) cell is
    // finest and the wall cell coarsest — clustering at d1 from above.
    let hi_sizes = layer_cell_sizes(b - d1, ny_hi); // index 0 = interface cell (finest)
    let mut y = d1;
    for (k, &h) in hi_sizes.iter().enumerate() {
        y += h;
        if k + 1 == ny_hi {
            ys.push(b); // snap the top wall to exactly b
        } else {
            ys.push(y);
        }
    }

    ys
}

/// Build a structured WR-90-style horizontal-slab [`TriMesh2D`] with `nx`
/// uniform columns in `x` and a **geometrically interface-graded** row
/// distribution in `y` (Phase 1.3.1.1 step 5.4).
///
/// The cross-section spans `[0, a] × [0, b]`. The dielectric layer
/// `0 ≤ y ≤ d1` is tagged material `1`; the air layer `d1 ≤ y ≤ b` is
/// tagged material `0`. The `y`-grid lines come from [`graded_y_lines`]
/// (clustered toward the interface `y = d1`, finest cells there, a node
/// placed exactly at `d1`), while `x` stays uniform (the dominant mode
/// varies slowly in `x`). Each rectangular cell splits along its
/// `(low-x, low-y) → (high-x, high-y)` diagonal into two CCW triangles —
/// the same split the uniform builders use — so the edge-table /
/// assembly path is unchanged.
///
/// Because a grid line sits exactly at `d1`, every cell lies entirely in
/// one layer; the triangle material is decided by the cell's `y`-midpoint
/// (`< d1` ⇒ dielectric tag 1, else air tag 0), and **no element straddles
/// the interface** — the material partition stays sharp under grading.
///
/// `ratio` is the geometric grading factor (`1` = uniform rows; larger =
/// finer at the interface). `ny_lo` / `ny_hi` are the row counts of the
/// dielectric / air layers.
///
/// This builder is **additive**: the uniform `horizontal_slab_mesh`
/// fixtures in the test/solve modules keep their builders and values; this
/// one only adds interface clustering for the step-5.4 convergence study.
#[allow(dead_code)] // step-5.4 convergence-study helper; consumed by mesh unit tests + the mirror in tests/eigensolver_inhomogeneous.rs
pub(crate) fn horizontal_slab_graded_mesh(
    a: f64,
    b: f64,
    d1: f64,
    nx: usize,
    ny_lo: usize,
    ny_hi: usize,
    ratio: f64,
) -> TriMesh2D {
    let ys = graded_y_lines(b, d1, ny_lo, ny_hi, ratio);
    let ny = ys.len() - 1;

    let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
    for &yj in &ys {
        for i in 0..=nx {
            vertices.push([a * (i as f64) / (nx as f64), yj]);
        }
    }
    let idx = |i: usize, j: usize| j * (nx + 1) + i;
    let mut triangles = Vec::with_capacity(2 * nx * ny);
    let mut tags: Vec<MaterialTag> = Vec::with_capacity(2 * nx * ny);
    for j in 0..ny {
        // Cell y-midpoint decides the layer; the interface node at d1 means
        // no cell straddles, so this is exact.
        let yc = 0.5 * (ys[j] + ys[j + 1]);
        let tag = if yc < d1 { 1u32 } else { 0u32 };
        for i in 0..nx {
            let v00 = idx(i, j);
            let v10 = idx(i + 1, j);
            let v11 = idx(i + 1, j + 1);
            let v01 = idx(i, j + 1);
            triangles.push([v00, v10, v11]);
            tags.push(tag);
            triangles.push([v00, v11, v01]);
            tags.push(tag);
        }
    }
    TriMesh2D::new(vertices, triangles, None, Some(tags)).unwrap()
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

    // ── Phase 1.3.1.1 step 5.4 — interface-graded horizontal-slab mesh ──

    const WR90_A: f64 = 22.86e-3;
    const WR90_B: f64 = 10.16e-3;

    #[test]
    fn graded_y_lines_endpoints_and_interface_exact() {
        // The grid must span exactly [0, b], be strictly increasing, contain
        // d1 EXACTLY (so the material partition stays sharp), and have the
        // documented node count `ny_lo + ny_hi + 1`.
        let (d1, ny_lo, ny_hi, r) = (WR90_B / 2.0, 6usize, 6usize, 1.5);
        let ys = graded_y_lines(WR90_B, d1, ny_lo, ny_hi, r);
        assert_eq!(ys.len(), ny_lo + ny_hi + 1, "node count");
        assert_eq!(ys[0], 0.0, "first node at the y=0 wall");
        assert_eq!(*ys.last().unwrap(), WR90_B, "last node at the y=b wall");
        // Strictly increasing.
        for w in ys.windows(2) {
            assert!(w[1] > w[0], "y-lines must strictly increase: {w:?}");
        }
        // d1 present exactly once, at the layer boundary index ny_lo.
        assert_eq!(ys[ny_lo], d1, "interface node must be exactly d1");
        let n_at_d1 = ys.iter().filter(|&&y| (y - d1).abs() < 1e-15).count();
        assert_eq!(n_at_d1, 1, "d1 must appear exactly once");
    }

    #[test]
    fn graded_y_lines_cluster_toward_interface() {
        // With ratio > 1 the cells must shrink toward the interface: the
        // cell touching d1 (from either side) is the finest in its layer,
        // and the wall cell the coarsest.
        let (d1, ny_lo, ny_hi, r) = (WR90_B / 2.0, 6usize, 6usize, 1.6);
        let ys = graded_y_lines(WR90_B, d1, ny_lo, ny_hi, r);

        // Lower-layer cell sizes (indices 0..ny_lo). Last one abuts d1.
        let lo: Vec<f64> = (0..ny_lo).map(|j| ys[j + 1] - ys[j]).collect();
        assert!(
            lo[ny_lo - 1] < lo[0],
            "lower-layer interface cell {} must be finer than the wall cell {}",
            lo[ny_lo - 1],
            lo[0]
        );
        // Geometric ratio between consecutive lower-layer cells ≈ r.
        for w in lo.windows(2) {
            let q = w[0] / w[1]; // coarser (wall side) / finer (interface side)
            assert!(
                (q - r).abs() < 1e-9,
                "lower-layer cells must grade by r={r}, got ratio {q}"
            );
        }

        // Upper-layer cell sizes (indices ny_lo..ny). First one abuts d1.
        let hi: Vec<f64> = (ny_lo..ys.len() - 1).map(|j| ys[j + 1] - ys[j]).collect();
        assert!(
            hi[0] < hi[hi.len() - 1],
            "upper-layer interface cell {} must be finer than the wall cell {}",
            hi[0],
            hi[hi.len() - 1]
        );
    }

    #[test]
    fn graded_y_lines_uniform_when_ratio_one() {
        // ratio = 1 must recover a uniform grid in each layer (additive: the
        // graded builder degenerates exactly to the uniform spacing).
        let (d1, ny_lo, ny_hi) = (WR90_B / 2.0, 4usize, 4usize);
        let ys = graded_y_lines(WR90_B, d1, ny_lo, ny_hi, 1.0);
        let h_lo = d1 / ny_lo as f64;
        let h_hi = (WR90_B - d1) / ny_hi as f64;
        for j in 0..ny_lo {
            assert!(
                (ys[j + 1] - ys[j] - h_lo).abs() < 1e-12,
                "lower layer uniform"
            );
        }
        for j in ny_lo..ys.len() - 1 {
            assert!(
                (ys[j + 1] - ys[j] - h_hi).abs() < 1e-12,
                "upper layer uniform"
            );
        }
    }

    #[test]
    fn graded_mesh_tags_sharp_and_no_straddle() {
        // The graded mesh must tag the dielectric layer (y < d1) as 1 and air
        // as 0, with NO triangle straddling the interface — every triangle's
        // three vertices lie on one side of (or on) y = d1 consistent with
        // its tag. This is the §5(c) "grading must not distort the interface
        // tag boundary" guard.
        let d1 = WR90_B / 2.0;
        let mesh = horizontal_slab_graded_mesh(WR90_A, WR90_B, d1, 4, 6, 6, 1.5);
        assert_eq!(
            mesh.triangle_material.len(),
            mesh.triangles.len(),
            "one tag per triangle"
        );
        let mut saw_diel = false;
        let mut saw_air = false;
        for (t, tri) in mesh.triangles.iter().enumerate() {
            let tag = mesh.triangle_material[t];
            let ys: Vec<f64> = tri.iter().map(|&v| mesh.vertices[v][1]).collect();
            let ymin = ys.iter().cloned().fold(f64::INFINITY, f64::min);
            let ymax = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            if tag == 1 {
                saw_diel = true;
                // Dielectric triangle: entirely at or below d1.
                assert!(
                    ymax <= d1 + 1e-12,
                    "dielectric triangle {t} straddles d1 (ymax={ymax}, d1={d1})"
                );
            } else {
                saw_air = true;
                // Air triangle: entirely at or above d1.
                assert!(
                    ymin >= d1 - 1e-12,
                    "air triangle {t} straddles d1 (ymin={ymin}, d1={d1})"
                );
            }
        }
        assert!(saw_diel && saw_air, "mesh must have both layers populated");
    }

    #[test]
    fn graded_mesh_node_and_edge_counts() {
        // Sanity: a graded nx × (ny_lo+ny_hi) slab has the structured-grid
        // vertex count and builds a valid edge table (boundary edges form the
        // rectangle perimeter). Mirrors the uniform-grid bookkeeping.
        let (nx, ny_lo, ny_hi) = (4usize, 5usize, 5usize);
        let ny = ny_lo + ny_hi;
        let mesh = horizontal_slab_graded_mesh(WR90_A, WR90_B, WR90_B / 2.0, nx, ny_lo, ny_hi, 1.4);
        assert_eq!(mesh.n_verts(), (nx + 1) * (ny + 1));
        assert_eq!(mesh.n_tris(), 2 * nx * ny);
        let table = EdgeTable::build(&mesh);
        // A structured nx×ny tri grid has 2·nx·ny+nx+ny boundary+interior
        // edges; the perimeter (boundary) count is 2(nx+ny).
        let n_boundary = table.is_boundary.iter().filter(|&&x| x).count();
        assert_eq!(n_boundary, 2 * (nx + ny), "perimeter boundary-edge count");
    }
}
