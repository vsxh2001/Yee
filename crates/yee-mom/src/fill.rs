//! MoM impedance matrix assembly.
//!
//! Implements the mixed-potential integral equation (MPIE) form for
//! free-space PEC surfaces:
//!
//! ```text
//! Z_{mn} = j ω μ₀ ⟨f_m, A_n⟩ + (1 / (j ω ε₀)) ⟨∇·f_m, φ_n⟩
//! ```
//!
//! Evaluated as nested quadrature over triangle pairs. Inner integration
//! switches to Duffy when triangles share a vertex / edge / face; otherwise
//! straight Gauss order 5 in both inner and outer.
//!
//! Reference: Gibson, *The Method of Moments in Electromagnetics* (2nd ed.,
//! 2014), Ch. 7.
#![allow(dead_code)]
// solve.rs (Task 9) consumes impedance_matrix.

use crate::basis::RwgBasis;
use crate::greens::FreeSpaceGreen;
use crate::quadrature::{DuffyTopology, DuffyTransform, GaussTriangle};
use faer::Mat;
use nalgebra::Vector3;
use num_complex::Complex64;
use rayon::prelude::*;

pub(crate) fn impedance_matrix(basis: &RwgBasis, green: &FreeSpaceGreen) -> Mat<Complex64> {
    let n = basis.n_basis();

    // Row-parallel fill. faer::Mat is not Sync across rows, so build
    // Vec<Vec<_>> in parallel then copy into the matrix.
    let rows: Vec<Vec<Complex64>> = (0..n)
        .into_par_iter()
        .map(|m| {
            let gauss = GaussTriangle::order_5();
            let mut row = vec![Complex64::new(0.0, 0.0); n];
            for (nidx, cell) in row.iter_mut().enumerate() {
                *cell = matrix_element(basis, green, &gauss, m, nidx);
            }
            row
        })
        .collect();

    let mut z = Mat::<Complex64>::zeros(n, n);
    for m in 0..n {
        for nidx in 0..n {
            z[(m, nidx)] = rows[m][nidx];
        }
    }
    z
}

fn matrix_element(
    basis: &RwgBasis,
    green: &FreeSpaceGreen,
    gauss: &GaussTriangle,
    m: usize,
    n: usize,
) -> Complex64 {
    let em = &basis.edges[m];
    let en = &basis.edges[n];

    let mut z_mn = Complex64::new(0.0, 0.0);
    for &t_outer in &[em.tri_plus, em.tri_minus] {
        for &t_inner in &[en.tri_plus, en.tri_minus] {
            z_mn += pair_contribution(basis, green, gauss, m, n, t_outer, t_inner);
        }
    }
    z_mn
}

fn pair_contribution(
    basis: &RwgBasis,
    green: &FreeSpaceGreen,
    gauss: &GaussTriangle,
    m: usize,
    n: usize,
    t_outer: u32,
    t_inner: u32,
) -> Complex64 {
    let outer_v = triangle_vertices(basis, t_outer);
    let inner_v = triangle_vertices(basis, t_inner);
    let outer_area = basis.areas[t_outer as usize];
    let inner_area = basis.areas[t_inner as usize];

    let div_m = basis.div(m, t_outer);
    let div_n = basis.div(n, t_inner);

    // k0 = ω/c real for free space. η0 stored on Green struct.
    let k0 = green.k0.re;
    let omega_mu0 = Complex64::new(0.0, 1.0) * Complex64::new(k0 * green.eta0, 0.0); // j k0 η0
    let inv_omega_eps0 = Complex64::new(0.0, -1.0) * Complex64::new(green.eta0 / k0, 0.0); // -j η0/k0

    let topology = topology_of(basis, t_outer, t_inner);

    // Duffy regularizes 1/R through its Jacobian — the integrand must remain
    // the FULL Green's function G(R), never `scalar_smooth`. Using
    // `scalar_smooth` (= G − 1/(4πR)) inside the Duffy path would double-
    // subtract the 1/R term and systematically bias every singular and
    // near-singular pair contribution.
    let integrand_duffy = |r_outer: Vector3<f64>, r_inner: Vector3<f64>| -> Complex64 {
        let fm = basis_value_at_point(basis, m, t_outer, r_outer, &outer_v);
        let fn_vec = basis_value_at_point(basis, n, t_inner, r_inner, &inner_v);
        let r = (r_outer - r_inner).norm();
        // `green.scalar` panics at r == 0. Dunavant order-5 has no vertex
        // point, so a Duffy sub-triangle Gauss point coinciding bit-exactly
        // with the outer anchor `r_outer` cannot occur in practice — but
        // the bit-exact guard is cheap and removes the panic risk. At r == 0
        // the Duffy Jacobian also vanishes, so the analytic limit
        // −j k0 / (4π) of G as R → 0 (matching `scalar_smooth`) is a safe
        // value here; it does not double-subtract because the Jacobian is
        // zero on the same point.
        let g = if r > 0.0 {
            green.scalar(r_outer, r_inner)
        } else {
            Complex64::new(0.0, -green.k0.re / (4.0 * std::f64::consts::PI))
        };
        omega_mu0 * Complex64::new(fm.dot(&fn_vec), 0.0) * g
            + inv_omega_eps0 * Complex64::new(div_m * div_n, 0.0) * g
    };

    // Well-separated pairs never have r near zero; full G is always defined.
    let integrand_gauss = |r_outer: Vector3<f64>, r_inner: Vector3<f64>| -> Complex64 {
        let fm = basis_value_at_point(basis, m, t_outer, r_outer, &outer_v);
        let fn_vec = basis_value_at_point(basis, n, t_inner, r_inner, &inner_v);
        let g = green.scalar(r_outer, r_inner);
        omega_mu0 * Complex64::new(fm.dot(&fn_vec), 0.0) * g
            + inv_omega_eps0 * Complex64::new(div_m * div_n, 0.0) * g
    };

    match topology {
        Some(t) => {
            let duffy = DuffyTransform {
                topology: t,
                outer_vertices: outer_v,
                inner_vertices: inner_v,
            };
            duffy.integrate(5, integrand_duffy)
        }
        None => {
            // Well-separated pair: straight nested Gauss.
            let mut acc = Complex64::new(0.0, 0.0);
            for (p_out, w_out) in gauss.points.iter().zip(gauss.weights.iter()) {
                let r_outer = bary_to_point(&outer_v, *p_out);
                for (p_in, w_in) in gauss.points.iter().zip(gauss.weights.iter()) {
                    let r_inner = bary_to_point(&inner_v, *p_in);
                    let val = integrand_gauss(r_outer, r_inner);
                    acc += Complex64::new(*w_out * *w_in * outer_area * inner_area, 0.0) * val;
                }
            }
            acc
        }
    }
}

fn triangle_vertices(basis: &RwgBasis, tri: u32) -> [Vector3<f64>; 3] {
    let [a, b, c] = basis.mesh.triangles[tri as usize];
    [
        basis.mesh.vertices[a as usize],
        basis.mesh.vertices[b as usize],
        basis.mesh.vertices[c as usize],
    ]
}

// PERF/DRY(phase-1.1): consolidate with quadrature::bary_to_point.
fn bary_to_point(v: &[Vector3<f64>; 3], bary: [f64; 3]) -> Vector3<f64> {
    bary[0] * v[0] + bary[1] * v[1] + bary[2] * v[2]
}

/// Reconstruct barycentric coordinates for `r` projected into the plane of
/// `tri_v`, then evaluate basis `k`. This is more work than carrying the
/// barycentric coord through, but isolates the API at Phase 1.0.
// PERF(phase-1.1): hoist barycentric pre-computation.
fn basis_value_at_point(
    basis: &RwgBasis,
    k: usize,
    tri: u32,
    r: Vector3<f64>,
    tri_v: &[Vector3<f64>; 3],
) -> Vector3<f64> {
    let v0 = tri_v[0];
    let e1 = tri_v[1] - v0;
    let e2 = tri_v[2] - v0;
    let d = r - v0;
    let g11 = e1.dot(&e1);
    let g12 = e1.dot(&e2);
    let g22 = e2.dot(&e2);
    let rhs1 = e1.dot(&d);
    let rhs2 = e2.dot(&d);
    let det = g11 * g22 - g12 * g12;
    let b1 = (g22 * rhs1 - g12 * rhs2) / det;
    let b2 = (-g12 * rhs1 + g11 * rhs2) / det;
    let b0 = 1.0 - b1 - b2;
    basis.eval(k, tri, [b0, b1, b2])
}

/// Classify a triangle pair by how many vertices they share.
fn topology_of(basis: &RwgBasis, t1: u32, t2: u32) -> Option<DuffyTopology> {
    if t1 == t2 {
        return Some(DuffyTopology::SameTriangle);
    }
    let [a1, b1, c1] = basis.mesh.triangles[t1 as usize];
    let [a2, b2, c2] = basis.mesh.triangles[t2 as usize];
    let set1 = [a1, b1, c1];
    let set2 = [a2, b2, c2];
    let shared = set1.iter().filter(|v| set2.contains(v)).count();
    match shared {
        0 => None,
        1 => Some(DuffyTopology::SharedVertex),
        2 => Some(DuffyTopology::SharedEdge),
        _ => Some(DuffyTopology::SameTriangle),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::RwgBasis;
    use nalgebra::Vector3;
    use yee_mesh::TriMesh;

    fn two_tri_mesh() -> TriMesh {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.1, 0.0, 0.0),
            Vector3::new(0.1, 0.1, 0.0),
            Vector3::new(0.0, 0.1, 0.0),
        ];
        let triangles = vec![[0u32, 1, 2], [0u32, 2, 3]];
        let tags = vec![0u32, 0u32];
        TriMesh::new(vertices, triangles, tags).unwrap()
    }

    #[test]
    fn two_rwg_matrix_is_finite_and_symmetric() {
        let basis = RwgBasis::from_mesh(two_tri_mesh()).unwrap();
        let green = FreeSpaceGreen::new(1.0e9);
        let z = impedance_matrix(&basis, &green);
        let n = basis.n_basis();
        assert!(n >= 1);
        for m in 0..n {
            for nidx in 0..n {
                let a = z[(m, nidx)];
                let b = z[(nidx, m)];
                assert!(a.re.is_finite() && a.im.is_finite());
                // Reciprocal MoM: Z is symmetric (NOT Hermitian). Tightened
                // to 1e-9 per plan spec — catches genuine reciprocity
                // violations that a looser tolerance would hide.
                assert!(
                    (a - b).norm() < 1.0e-9 * a.norm().max(1.0),
                    "asymmetry Z[{m},{nidx}]={a} vs Z[{nidx},{m}]={b}"
                );
            }
        }
    }
}
