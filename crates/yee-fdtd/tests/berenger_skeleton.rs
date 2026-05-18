//! Phase 2.fdtd.7.x B1 — Berenger Huygens-surface skeleton tests.
//!
//! Pins the type surface added for the Berenger 2006 fine → coarse closure:
//!
//! - [`BerengerHuygensFace`] — the six face identifiers.
//! - [`assign_edge_to_face`] — the lower-numbered-axis-wins edge-ownership
//!   tie-break for the 12 cuboid edges.
//! - [`SubgridRegion::face_index_table`] — per-face coarse-cell enumeration.
//! - [`SubgridRegion::inject_equivalent_currents_to_coarse`] — the B1 stub
//!   (current re-radiation lands in B2).
//!
//! Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-x-berenger-huygens-design.md`.
//! Plan: `docs/superpowers/plans/2026-05-19-phase-2-fdtd-7-x-berenger-huygens.md` step B1.
//! ADR:  `docs/src/decisions/0035-berenger-huygens-subgridding.md`.

use std::collections::HashMap;

use yee_fdtd::{BerengerHuygensFace, SubgridRegion, YeeGrid, assign_edge_to_face};

const N: usize = 16;
const DX: f64 = 1.0e-3;

fn parent_grid() -> YeeGrid {
    YeeGrid::vacuum(N, N, N, DX)
}

/// `lo = (2, 2, 2)`, `hi = (5, 5, 5)` → 3 coarse cells per axis →
/// 3 × 3 = 9 cells per Huygens face. Verifies the per-face cell counts
/// and the array ordering `[XLow, XHigh, YLow, YHigh, ZLow, ZHigh]`.
#[test]
fn face_index_table_cell_counts() {
    let parent = parent_grid();
    let region = SubgridRegion::new(&parent, (2, 2, 2), (5, 5, 5))
        .expect("valid subgrid bounds should construct");

    let table = region.face_index_table();
    for (face, cells) in BerengerHuygensFace::all().iter().zip(table.iter()) {
        assert_eq!(
            cells.len(),
            9,
            "face {face:?} should have 3 × 3 = 9 cells on a (3, 3, 3) coarse subgrid"
        );
    }

    // XLow lives at i = lo.0 − 1 = 1; XHigh at i = hi.0 = 5; etc.
    let [x_low, x_high, y_low, y_high, z_low, z_high] = &table;
    for &(i, _, _) in x_low {
        assert_eq!(i, 1, "XLow cells live at i = lo.0 - 1 = 1");
    }
    for &(i, _, _) in x_high {
        assert_eq!(i, 5, "XHigh cells live at i = hi.0 = 5");
    }
    for &(_, j, _) in y_low {
        assert_eq!(j, 1, "YLow cells live at j = lo.1 - 1 = 1");
    }
    for &(_, j, _) in y_high {
        assert_eq!(j, 5, "YHigh cells live at j = hi.1 = 5");
    }
    for &(_, _, k) in z_low {
        assert_eq!(k, 1, "ZLow cells live at k = lo.2 - 1 = 1");
    }
    for &(_, _, k) in z_high {
        assert_eq!(k, 5, "ZHigh cells live at k = hi.2 = 5");
    }
}

/// Spot-check the spec's three example edge assignments. The rule is
/// "lower-numbered axis wins" with `X = 0`, `Y = 1`, `Z = 2`.
#[test]
fn assign_edge_to_face_lower_axis_wins() {
    use BerengerHuygensFace::*;

    assert_eq!(assign_edge_to_face(XLow, YLow), XLow);
    assert_eq!(assign_edge_to_face(YLow, ZLow), YLow);
    assert_eq!(assign_edge_to_face(XHigh, ZLow), XHigh);

    // Commutativity: argument order does not change the winner.
    assert_eq!(assign_edge_to_face(YLow, XLow), XLow);
    assert_eq!(assign_edge_to_face(ZLow, YLow), YLow);
    assert_eq!(assign_edge_to_face(ZLow, XHigh), XHigh);
}

/// The B1 stub is a no-op: calling it leaves the parent grid's E and H
/// arrays bit-exact equal to their pre-call values. The actual current
/// accumulation lands in B2.
#[test]
fn inject_equivalent_currents_to_coarse_is_currently_noop() {
    let mut parent = parent_grid();

    // Seed the parent with a non-trivial field pattern so a stub that
    // accidentally zeros / overwrites a cell would be caught.
    for i in 0..parent.ex.shape()[0] {
        for j in 0..parent.ex.shape()[1] {
            for k in 0..parent.ex.shape()[2] {
                parent.ex[(i, j, k)] = 1.0 + (i + j + k) as f64;
            }
        }
    }
    for i in 0..parent.ey.shape()[0] {
        for j in 0..parent.ey.shape()[1] {
            for k in 0..parent.ey.shape()[2] {
                parent.ey[(i, j, k)] = 2.0 + (i + 2 * j + k) as f64;
            }
        }
    }
    for i in 0..parent.ez.shape()[0] {
        for j in 0..parent.ez.shape()[1] {
            for k in 0..parent.ez.shape()[2] {
                parent.ez[(i, j, k)] = 3.0 + (i + j + 2 * k) as f64;
            }
        }
    }
    for i in 0..parent.hx.shape()[0] {
        for j in 0..parent.hx.shape()[1] {
            for k in 0..parent.hx.shape()[2] {
                parent.hx[(i, j, k)] = -1.0 - (i + j + k) as f64;
            }
        }
    }
    for i in 0..parent.hy.shape()[0] {
        for j in 0..parent.hy.shape()[1] {
            for k in 0..parent.hy.shape()[2] {
                parent.hy[(i, j, k)] = -2.0 - (i + 2 * j + k) as f64;
            }
        }
    }
    for i in 0..parent.hz.shape()[0] {
        for j in 0..parent.hz.shape()[1] {
            for k in 0..parent.hz.shape()[2] {
                parent.hz[(i, j, k)] = -3.0 - (i + j + 2 * k) as f64;
            }
        }
    }

    let region = SubgridRegion::new(&parent, (2, 2, 2), (5, 5, 5))
        .expect("valid subgrid bounds should construct");

    let ex_before = parent.ex.clone();
    let ey_before = parent.ey.clone();
    let ez_before = parent.ez.clone();
    let hx_before = parent.hx.clone();
    let hy_before = parent.hy.clone();
    let hz_before = parent.hz.clone();

    region.inject_equivalent_currents_to_coarse(&mut parent);

    // Bit-exact equality — the B1 stub must touch nothing.
    assert_eq!(parent.ex, ex_before, "stub must not modify coarse E_x");
    assert_eq!(parent.ey, ey_before, "stub must not modify coarse E_y");
    assert_eq!(parent.ez, ez_before, "stub must not modify coarse E_z");
    assert_eq!(parent.hx, hx_before, "stub must not modify coarse H_x");
    assert_eq!(parent.hy, hy_before, "stub must not modify coarse H_y");
    assert_eq!(parent.hz, hz_before, "stub must not modify coarse H_z");
}

/// Enumerate all 12 cuboid edges (each shared by exactly two faces) and
/// verify [`assign_edge_to_face`] picks exactly one owner per edge. The
/// 12 edges are the pairs `{face_a, face_b}` whose axes differ; pairs
/// of opposite faces along the same axis (e.g. `XLow`–`XHigh`) do not
/// share an edge and are excluded.
#[test]
fn no_edge_double_counted() {
    use BerengerHuygensFace::*;

    let faces = BerengerHuygensFace::all();
    let mut edges: Vec<(BerengerHuygensFace, BerengerHuygensFace)> = Vec::new();
    for (i, &face_a) in faces.iter().enumerate() {
        for &face_b in faces.iter().skip(i + 1) {
            if face_a.axis() != face_b.axis() {
                edges.push((face_a, face_b));
            }
        }
    }

    // Sanity: 6 faces × 4 adjacent neighbours / 2 (each edge counted
    // once) = 12 cuboid edges.
    assert_eq!(edges.len(), 12, "cuboid has exactly 12 edges");

    // Every edge resolves to exactly one of its two faces.
    for &(face_a, face_b) in &edges {
        let owner = assign_edge_to_face(face_a, face_b);
        assert!(
            owner == face_a || owner == face_b,
            "edge ({face_a:?}, {face_b:?}) must be owned by one of its two faces, got {owner:?}"
        );
    }

    // Tally per-face ownership counts: each face owns the edges it
    // shares with higher-axis faces. With axes X=0 / Y=1 / Z=2 and
    // each face having 2 same-axis-orientation neighbours per other
    // axis, the counts come out as XLow=4, XHigh=4, YLow=2, YHigh=2,
    // ZLow=0, ZHigh=0 — and the sum is 12, matching the edge count.
    let mut counts: HashMap<BerengerHuygensFace, usize> = HashMap::new();
    for &(face_a, face_b) in &edges {
        *counts
            .entry(assign_edge_to_face(face_a, face_b))
            .or_insert(0) += 1;
    }
    assert_eq!(counts.get(&XLow).copied().unwrap_or(0), 4);
    assert_eq!(counts.get(&XHigh).copied().unwrap_or(0), 4);
    assert_eq!(counts.get(&YLow).copied().unwrap_or(0), 2);
    assert_eq!(counts.get(&YHigh).copied().unwrap_or(0), 2);
    assert_eq!(counts.get(&ZLow).copied().unwrap_or(0), 0);
    assert_eq!(counts.get(&ZHigh).copied().unwrap_or(0), 0);

    // Total ownership = 12 edges, no double-count.
    let total: usize = counts.values().sum();
    assert_eq!(total, 12, "every edge must be owned by exactly one face");
}
