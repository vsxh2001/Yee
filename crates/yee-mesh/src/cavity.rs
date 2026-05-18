//! Hand-rolled rectangular-cavity tetrahedralization for the Phase 4
//! FEM eigenmode solver.
//!
//! Implements the canonical Kuhn 6-tet decomposition (Kuhn 1960; see
//! <https://en.wikipedia.org/wiki/Tetrahedron#Subdivision_into_orthoschemes>)
//! of an axis-aligned brick. A regular `(nx+1) × (ny+1) × (nz+1)` grid
//! of vertices is laid over `[0, a] × [0, b] × [0, d]`, each brick is
//! split into the same six Kuhn tets sharing the body diagonal from the
//! brick's `(0,0,0)` corner to its `(1,1,1)` corner, and every tet is
//! tagged with the default [`crate::MaterialTag`] (currently `0` —
//! there is no `MaterialTag::AIR` variant since `MaterialTag` is a bare
//! `u32` alias, not an enum; see the module-level finding below).
//!
//! Properties guaranteed by the Kuhn decomposition for first-order
//! Nedelec edge elements:
//!
//! * **No slivers.** The six tets are congruent orthoschemes; each has
//!   the same volume `(1/6) · (a/nx) · (b/ny) · (d/nz)`. Aspect ratio
//!   is bounded by the brick aspect ratio, so a near-cubic brick gives
//!   near-equilateral tets in the conditioning sense relevant to the
//!   Nedelec stiffness/mass spectrum.
//! * **Consistent orientation across brick boundaries.** Each pair of
//!   adjacent bricks splits its shared face into the same two
//!   triangles, so edge DoFs glue cleanly without remeshing the
//!   interface (a property the spec calls out in §5 / §11).
//! * **Positive signed volume on every tet** — although the Kuhn
//!   table written below is sign-flipped on three of the six tets,
//!   [`crate::TetMesh3D::new`] silently re-orients them (T2 contract),
//!   so the returned mesh has uniformly positive signed volumes.
//!
//! **Out-of-lane finding:** `crate::MaterialTag` is currently a bare
//! `u32` alias (see `crates/yee-mesh/src/lib.rs`). There is no
//! `MaterialTag::AIR` enum variant — the v0 plan refers to one
//! aspirationally. This constructor tags every tet with `0` (the
//! default tag used by `TetMesh3D::new`'s `None` path). Introducing a
//! proper enum is a separate, cross-cutting change (`yee-mom` and
//! `yee-py` both consume the `u32` keying via `HashMap<MaterialTag,
//! Complex64>` for `eps_r` / `mu_r` lookup), so it is surfaced as a
//! finding rather than fixed in this lane.

use crate::{Error, MaterialTag, Result, TetMesh3D};
use nalgebra::Vector3;

/// Default material tag applied to every tet emitted by
/// [`TetMesh3D::cavity_uniform`]. Currently `0` (free space / air);
/// see the module-level out-of-lane finding for why this is not a
/// named `MaterialTag::AIR` variant.
const DEFAULT_AIR_TAG: MaterialTag = 0;

/// Kuhn 6-tet decomposition of a unit cube, indexed by the local
/// corner numbering documented in [`TetMesh3D::cavity_uniform`].
///
/// Each row lists the four corners of one Kuhn tet as
/// `[c0, c1, c2, c3]`. All six tets share the body diagonal `c0 = 0`
/// → `c3 = 7` (i.e. corners `(0,0,0)` and `(1,1,1)` of the brick),
/// and the two intermediate corners trace one of the 3! = 6 monotone
/// edge paths from corner 0 to corner 7.
///
/// Three of the six tets below have negative signed volume for the
/// canonical ordering of brick corners; [`TetMesh3D::new`] silently
/// reorients them (T2 contract).
const KUHN_LOCAL_TETS: [[usize; 4]; 6] = [
    // path x → y → z : 0 → 1 → 3 → 7
    [0, 1, 3, 7],
    // path x → z → y : 0 → 1 → 5 → 7
    [0, 1, 5, 7],
    // path y → x → z : 0 → 2 → 3 → 7
    [0, 2, 3, 7],
    // path y → z → x : 0 → 2 → 6 → 7
    [0, 2, 6, 7],
    // path z → x → y : 0 → 4 → 5 → 7
    [0, 4, 5, 7],
    // path z → y → x : 0 → 4 → 6 → 7
    [0, 4, 6, 7],
];

impl TetMesh3D {
    /// Build a uniform rectangular-cavity tet mesh on
    /// `[0, a] × [0, b] × [0, d]` with `nx × ny × nz` axis-aligned
    /// bricks, each brick decomposed into six Kuhn tetrahedra.
    ///
    /// The returned mesh has:
    /// * `(nx + 1) · (ny + 1) · (nz + 1)` vertices on a regular grid;
    /// * `nx · ny · nz · 6` tetrahedra (Kuhn decomposition per brick);
    /// * `0` tag (the default / air tag — see module-level finding)
    ///   on every vertex and every tet;
    /// * strictly positive signed volume on every tet (T2 re-orients
    ///   silently as needed).
    ///
    /// **Local corner indexing within each brick** — used by the
    /// `KUHN_LOCAL_TETS` table:
    ///
    /// ```text
    ///   bit 0 (1) → x-offset (0 → 0,  1 → 1)
    ///   bit 1 (2) → y-offset (0 → 0,  1 → 1)
    ///   bit 2 (4) → z-offset (0 → 0,  1 → 1)
    ///
    ///   0 = (0,0,0)   1 = (1,0,0)   2 = (0,1,0)   3 = (1,1,0)
    ///   4 = (0,0,1)   5 = (1,0,1)   6 = (0,1,1)   7 = (1,1,1)
    /// ```
    ///
    /// # Arguments
    ///
    /// * `a`, `b`, `d` — cavity extents along x, y, z (metres).
    /// * `nx`, `ny`, `nz` — number of bricks along each axis. All
    ///   three must be `>= 1`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Invalid`] if any of `a`, `b`, `d` is not
    /// strictly positive and finite, or if any of `nx`, `ny`, `nz`
    /// is zero. The downstream [`TetMesh3D::new`] call may also
    /// surface a validation failure if a tet ends up degenerate
    /// (should not happen for any caller-supplied positive extents
    /// and non-zero subdivisions — guarded for completeness).
    ///
    /// # Examples
    ///
    /// ```
    /// use yee_mesh::TetMesh3D;
    /// // WR-90-based cavity: 22.86 × 10.16 × 30 mm, 4 × 2 × 4 bricks.
    /// let mesh = TetMesh3D::cavity_uniform(0.02286, 0.01016, 0.030, 4, 2, 4).unwrap();
    /// assert_eq!(mesh.n_tets(), 4 * 2 * 4 * 6);
    /// assert_eq!(mesh.n_verts(), 5 * 3 * 5);
    /// ```
    pub fn cavity_uniform(
        a: f64,
        b: f64,
        d: f64,
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> Result<TetMesh3D> {
        if !(a.is_finite() && a > 0.0) {
            return Err(Error::Invalid(format!(
                "cavity_uniform: a must be a positive finite extent, got {a}"
            )));
        }
        if !(b.is_finite() && b > 0.0) {
            return Err(Error::Invalid(format!(
                "cavity_uniform: b must be a positive finite extent, got {b}"
            )));
        }
        if !(d.is_finite() && d > 0.0) {
            return Err(Error::Invalid(format!(
                "cavity_uniform: d must be a positive finite extent, got {d}"
            )));
        }
        if nx == 0 || ny == 0 || nz == 0 {
            return Err(Error::Invalid(format!(
                "cavity_uniform: nx, ny, nz must all be >= 1, got ({nx}, {ny}, {nz})"
            )));
        }

        let dx = a / nx as f64;
        let dy = b / ny as f64;
        let dz = d / nz as f64;

        let nvx = nx + 1;
        let nvy = ny + 1;
        let nvz = nz + 1;

        // Build the regular vertex grid. Index layout:
        //   global_id(i, j, k) = i + j * nvx + k * nvx * nvy
        // matches the local-corner bit layout in `KUHN_LOCAL_TETS`
        // (bit 0 → x, bit 1 → y, bit 2 → z), so a local corner
        // `c ∈ 0..8` of the brick anchored at `(i, j, k)` resolves to
        // `global_id(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1))`.
        let mut vertices = Vec::with_capacity(nvx * nvy * nvz);
        for k in 0..nvz {
            for j in 0..nvy {
                for i in 0..nvx {
                    vertices.push(Vector3::new(i as f64 * dx, j as f64 * dy, k as f64 * dz));
                }
            }
        }

        let global_id = |i: usize, j: usize, k: usize| -> usize { i + j * nvx + k * nvx * nvy };

        // For each brick (i, j, k), emit the six Kuhn tets. We resolve
        // each local-corner index through the bit layout described
        // above, then push the four global vertex indices.
        let mut tetrahedra = Vec::with_capacity(nx * ny * nz * 6);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    for local in &KUHN_LOCAL_TETS {
                        let mut tet = [0usize; 4];
                        for (slot, &c) in local.iter().enumerate() {
                            let di = c & 1;
                            let dj = (c >> 1) & 1;
                            let dk = (c >> 2) & 1;
                            tet[slot] = global_id(i + di, j + dj, k + dk);
                        }
                        tetrahedra.push(tet);
                    }
                }
            }
        }

        let n_verts = vertices.len();
        let n_tets = tetrahedra.len();
        let vertex_material = vec![DEFAULT_AIR_TAG; n_verts];
        let tetrahedron_material = vec![DEFAULT_AIR_TAG; n_tets];

        TetMesh3D::new(
            vertices,
            tetrahedra,
            Some(vertex_material),
            Some(tetrahedron_material),
        )
    }
}
