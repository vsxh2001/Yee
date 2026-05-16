//! Hand-coded thin-cylinder mesh generator for the dipole validation.

use nalgebra::Vector3;
use yee_mesh::TriMesh;

/// Triangulates the lateral surface of a cylinder (no end caps).
///
/// The cylinder's axis is along `z`, centred at the origin. `length_m` is
/// the total length; `radius_m` is the cylinder radius. `n_axial` is the
/// number of axial segments (rings of triangles between adjacent z-cuts);
/// `n_around` is the number of segments around the circumference.
///
/// Two triangles are produced per `(axial × around)` cell, so the total
/// triangle count is `2 * n_axial * n_around`.
///
/// The two axial rings flanking the centre plane carry DIFFERENT non-zero
/// tags: the ring immediately below `z = 0` is tagged `1`, the ring
/// immediately above is tagged `2`, all others `0`. The single
/// circumferential edge ring at `z = 0` therefore lies on the boundary
/// between the two tagged regions and is picked up as the delta-gap port
/// by `RwgBasis::from_mesh`'s "different non-zero tags" port convention.
pub fn thin_cylinder(length_m: f64, radius_m: f64, n_axial: usize, n_around: usize) -> TriMesh {
    assert!(
        n_axial >= 2 && n_axial.is_multiple_of(2),
        "n_axial must be even and >= 2"
    );
    assert!(n_around >= 3, "n_around must be >= 3");

    let mut vertices: Vec<Vector3<f64>> = Vec::with_capacity((n_axial + 1) * n_around);
    let dz = length_m / (n_axial as f64);
    let z0 = -length_m / 2.0;
    let dtheta = std::f64::consts::TAU / (n_around as f64);

    for i in 0..=n_axial {
        let z = z0 + (i as f64) * dz;
        for j in 0..n_around {
            let theta = (j as f64) * dtheta;
            vertices.push(Vector3::new(
                radius_m * theta.cos(),
                radius_m * theta.sin(),
                z,
            ));
        }
    }

    let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(2 * n_axial * n_around);
    let mut tags: Vec<u32> = Vec::with_capacity(2 * n_axial * n_around);
    let central_ring = n_axial / 2;

    for i in 0..n_axial {
        for j in 0..n_around {
            let j_next = (j + 1) % n_around;
            let a = (i * n_around + j) as u32;
            let b = (i * n_around + j_next) as u32;
            let c = ((i + 1) * n_around + j_next) as u32;
            let d = ((i + 1) * n_around + j) as u32;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
            let tag = if i == central_ring - 1 {
                1
            } else if i == central_ring {
                2
            } else {
                0
            };
            tags.push(tag);
            tags.push(tag);
        }
    }

    TriMesh::new(vertices, triangles, tags).expect("cylinder mesh invariants")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangle_count_matches_formula() {
        let mesh = thin_cylinder(1.0, 0.005, 24, 24);
        assert_eq!(mesh.n_tris(), 2 * 24 * 24);
    }

    #[test]
    fn vertex_count_matches_formula() {
        let mesh = thin_cylinder(1.0, 0.005, 24, 24);
        assert_eq!(mesh.vertices.len(), 25 * 24);
    }

    #[test]
    fn central_ring_tag_counts() {
        let mesh = thin_cylinder(1.0, 0.005, 24, 24);
        let tagged_1 = mesh.tags.iter().filter(|&&t| t == 1).count();
        let tagged_2 = mesh.tags.iter().filter(|&&t| t == 2).count();
        // Each ring has n_around cells × 2 triangles = 48 triangles.
        assert_eq!(tagged_1, 2 * 24);
        assert_eq!(tagged_2, 2 * 24);
    }
}
