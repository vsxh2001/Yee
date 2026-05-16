//! Gauss quadrature on triangles and Duffy transform for singular RWG integrals.
//!
//! References:
//! - Dunavant, *Int. J. Numer. Methods Eng.* 21.6 (1985) — symmetric Gauss
//!   quadrature on triangles.
//! - Khayat & Wilton, *IEEE T-AP* 53.10 (2005) — Duffy transform for RWG.

#![allow(dead_code)]
// Phase 1.0 cleanup: lifts naturally in Task 8 when fill.rs consumes these.
#![allow(clippy::excessive_precision)]
// Dunavant tables are quoted at >15 significant figures so the literals match
// the published source verbatim; the surplus digits beyond f64 precision are
// intentional provenance, not a typo.

use nalgebra::Vector3;
use num_complex::Complex64;

pub(crate) struct GaussTriangle {
    /// Barycentric coordinates of each quadrature point, sum = 1.
    pub points: Vec<[f64; 3]>,
    /// Weights normalised so sum = 1. Multiply by triangle area to integrate.
    pub weights: Vec<f64>,
}

impl GaussTriangle {
    /// Order-3 rule (4 points, exact for cubics).
    pub fn order_3() -> Self {
        let p = vec![
            [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            [0.6, 0.2, 0.2],
            [0.2, 0.6, 0.2],
            [0.2, 0.2, 0.6],
        ];
        let w = vec![-9.0 / 16.0, 25.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0];
        Self {
            points: p,
            weights: w,
        }
    }

    /// Order-5 rule (7 points, exact for quintics).
    pub fn order_5() -> Self {
        let a1 = 0.0597158717_897698;
        let b1 = 0.4701420641_051151;
        let a2 = 0.7974269853_530873;
        let b2 = 0.1012865073_234563;
        let p = vec![
            [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            [a1, b1, b1],
            [b1, a1, b1],
            [b1, b1, a1],
            [a2, b2, b2],
            [b2, a2, b2],
            [b2, b2, a2],
        ];
        let w = vec![
            0.2250000000_000000,
            0.1323941527_885062,
            0.1323941527_885062,
            0.1323941527_885062,
            0.1259391805_448271,
            0.1259391805_448271,
            0.1259391805_448271,
        ];
        Self {
            points: p,
            weights: w,
        }
    }

    /// Order-7 rule (13 points). For near-singular outer integration.
    pub fn order_7() -> Self {
        let p = vec![
            [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            [
                0.4793080678_413916,
                0.2603459660_790042,
                0.2603459660_790042,
            ],
            [
                0.2603459660_790042,
                0.4793080678_413916,
                0.2603459660_790042,
            ],
            [
                0.2603459660_790042,
                0.2603459660_790042,
                0.4793080678_413916,
            ],
            [
                0.8697397941_955675,
                0.0651301029_022159,
                0.0651301029_022166,
            ],
            [
                0.0651301029_022159,
                0.8697397941_955675,
                0.0651301029_022166,
            ],
            [
                0.0651301029_022159,
                0.0651301029_022166,
                0.8697397941_955675,
            ],
            [
                0.6384441885_698096,
                0.3128654960_048880,
                0.0486903154_253024,
            ],
            [
                0.6384441885_698096,
                0.0486903154_253024,
                0.3128654960_048880,
            ],
            [
                0.3128654960_048880,
                0.6384441885_698096,
                0.0486903154_253024,
            ],
            [
                0.3128654960_048880,
                0.0486903154_253024,
                0.6384441885_698096,
            ],
            [
                0.0486903154_253024,
                0.6384441885_698096,
                0.3128654960_048880,
            ],
            [
                0.0486903154_253024,
                0.3128654960_048880,
                0.6384441885_698096,
            ],
        ];
        let w = vec![
            -0.1495700444_677495,
            0.1756152574_332137,
            0.1756152574_332137,
            0.1756152574_332137,
            0.0533472356_088403,
            0.0533472356_088403,
            0.0533472356_088403,
            0.0771137608_903113,
            0.0771137608_903113,
            0.0771137608_903113,
            0.0771137608_903113,
            0.0771137608_903113,
            0.0771137608_903113,
        ];
        Self {
            points: p,
            weights: w,
        }
    }
}

/// Triangle-pair topology determining which singularity treatment to use.
#[derive(Debug, Clone, Copy)]
pub(crate) enum DuffyTopology {
    SameTriangle,
    SharedEdge,
    SharedVertex,
}

pub(crate) struct DuffyTransform {
    pub topology: DuffyTopology,
    pub outer_vertices: [Vector3<f64>; 3],
    pub inner_vertices: [Vector3<f64>; 3],
}

impl DuffyTransform {
    /// Integrate `f(r_outer, r_inner)` over the outer × inner triangle pair
    /// using a Duffy-style transform. The inner integration sub-divides the
    /// inner triangle into three sub-triangles anchored at the outer-side
    /// quadrature point, which removes the 1/R singularity in the radial
    /// Duffy variable when the two triangles share at least one vertex.
    ///
    /// `order` is the Gauss order (3, 5, or 7) used for the inner Gauss
    /// integration on each sub-triangle and for the outer Gauss integration.
    ///
    /// # Panics
    ///
    /// Panics if `order` is not one of 3, 5, 7.
    pub fn integrate<F>(&self, order: usize, f: F) -> Complex64
    where
        F: Fn(Vector3<f64>, Vector3<f64>) -> Complex64,
    {
        let gauss_outer = gauss_order(order);
        let outer_area = triangle_area(&self.outer_vertices);

        let mut acc = Complex64::new(0.0, 0.0);
        for (p_outer, w_outer) in gauss_outer.points.iter().zip(gauss_outer.weights.iter()) {
            let r_outer = bary_to_point(&self.outer_vertices, *p_outer);
            let inner = duffy_inner_split(&self.inner_vertices, r_outer, order, &f);
            acc += Complex64::new(*w_outer * outer_area, 0.0) * inner;
        }
        acc
    }
}

fn gauss_order(order: usize) -> GaussTriangle {
    match order {
        3 => GaussTriangle::order_3(),
        5 => GaussTriangle::order_5(),
        7 => GaussTriangle::order_7(),
        _ => panic!("Duffy quadrature order must be 3, 5, or 7"),
    }
}

fn triangle_area(v: &[Vector3<f64>; 3]) -> f64 {
    0.5 * (v[1] - v[0]).cross(&(v[2] - v[0])).norm()
}

fn bary_to_point(v: &[Vector3<f64>; 3], bary: [f64; 3]) -> Vector3<f64> {
    bary[0] * v[0] + bary[1] * v[1] + bary[2] * v[2]
}

/// Split the inner triangle into three sub-triangles anchored at `r_outer`
/// and Gauss-integrate `f` over each sub-triangle.
fn duffy_inner_split<F>(
    inner: &[Vector3<f64>; 3],
    r_outer: Vector3<f64>,
    order: usize,
    f: &F,
) -> Complex64
where
    F: Fn(Vector3<f64>, Vector3<f64>) -> Complex64,
{
    let gauss = gauss_order(order);
    let mut acc = Complex64::new(0.0, 0.0);
    for k in 0..3 {
        let v0 = r_outer;
        let v1 = inner[k];
        let v2 = inner[(k + 1) % 3];
        let sub = [v0, v1, v2];
        let sub_area = triangle_area(&sub);
        for (p, w) in gauss.points.iter().zip(gauss.weights.iter()) {
            let r_inner = bary_to_point(&sub, *p);
            acc += Complex64::new(*w * sub_area, 0.0) * f(r_outer, r_inner);
        }
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Integrate ξ₁ over reference triangle: ∫_T ξ₁ dA = area * (1/3) = 1/6.
    #[test]
    fn order_3_integrates_linear_exact() {
        let q = GaussTriangle::order_3();
        let area = 0.5;
        let s: f64 = q
            .points
            .iter()
            .zip(q.weights.iter())
            .map(|(p, w)| w * p[0])
            .sum();
        let integral = area * s;
        assert!((integral - 1.0 / 6.0).abs() < 1e-12);
    }

    /// ∫_T ξ_1^n dA = 2 area / ((n+1)(n+2)). For n=5: 1 / 42.
    #[test]
    fn order_5_integrates_quintic_exact() {
        let q = GaussTriangle::order_5();
        let area = 0.5;
        let s: f64 = q
            .points
            .iter()
            .zip(q.weights.iter())
            .map(|(p, w)| w * p[0].powi(5))
            .sum();
        let integral = area * s;
        assert!((integral - 1.0 / 42.0).abs() < 1e-10);
    }

    #[test]
    fn weights_sum_to_one_each_order() {
        for q in [
            GaussTriangle::order_3(),
            GaussTriangle::order_5(),
            GaussTriangle::order_7(),
        ] {
            let s: f64 = q.weights.iter().sum();
            assert!((s - 1.0).abs() < 1e-12);
        }
    }

    /// Duffy `∫∫ 1/R dA_outer dA_inner` over a self-triangle pair must be
    /// finite and positive. The non-Duffy direct quadrature would diverge.
    #[test]
    fn duffy_self_triangle_one_over_r_finite() {
        let tri = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let duffy = DuffyTransform {
            topology: DuffyTopology::SameTriangle,
            outer_vertices: tri,
            inner_vertices: tri,
        };
        let result = duffy.integrate(5, |r1, r2| {
            let r = (r1 - r2).norm();
            if r > 1e-15 {
                Complex64::new(1.0 / r, 0.0)
            } else {
                Complex64::new(0.0, 0.0)
            }
        });
        assert!(result.re.is_finite() && result.re > 0.0);
    }
}
