//! Flat six-component field state shared by every backend.

use crate::spec::{FdtdSpec, idx3, len3};

/// The six staggered field components as flat row-major FP64 buffers.
///
/// Layout is `ndarray`'s default C order (`(i * dim_j + j) * dim_k + k`) with
/// the staggered shapes from [`FdtdSpec`], so a `YeeGrid` component maps onto
/// the matching buffer via `as_slice()` with no reindexing, and the GPU
/// backend uses the identical linearization in WGSL.
#[derive(Debug, Clone, PartialEq)]
pub struct Fields {
    /// `E_x`, shape `[nx, ny+1, nz+1]`.
    pub ex: Vec<f64>,
    /// `E_y`, shape `[nx+1, ny, nz+1]`.
    pub ey: Vec<f64>,
    /// `E_z`, shape `[nx+1, ny+1, nz]`.
    pub ez: Vec<f64>,
    /// `H_x`, shape `[nx+1, ny, nz]`.
    pub hx: Vec<f64>,
    /// `H_y`, shape `[nx, ny+1, nz]`.
    pub hy: Vec<f64>,
    /// `H_z`, shape `[nx, ny, nz+1]`.
    pub hz: Vec<f64>,
}

impl Fields {
    /// All-zero fields for `spec`.
    pub fn zero(spec: &FdtdSpec) -> Self {
        Self {
            ex: vec![0.0; len3(spec.ex_dims())],
            ey: vec![0.0; len3(spec.ey_dims())],
            ez: vec![0.0; len3(spec.ez_dims())],
            hx: vec![0.0; len3(spec.hx_dims())],
            hy: vec![0.0; len3(spec.hy_dims())],
            hz: vec![0.0; len3(spec.hz_dims())],
        }
    }

    /// Zero fields with a unit-amplitude Gaussian ball injected into `E_z`,
    /// centred on primary cell `center` with standard deviation `sigma_cells`
    /// (in cells). The standard initial condition for the E.0 parity gates.
    ///
    /// # Panics
    ///
    /// Panics if `center` lies outside the `E_z` array or `sigma_cells` is
    /// non-positive.
    pub fn with_gaussian_ez(
        spec: &FdtdSpec,
        center: (usize, usize, usize),
        sigma_cells: f64,
    ) -> Self {
        assert!(sigma_cells > 0.0, "sigma_cells must be positive");
        let dims = spec.ez_dims();
        let (ci, cj, ck) = center;
        assert!(
            ci < dims.0 && cj < dims.1 && ck < dims.2,
            "center must lie inside the E_z array"
        );

        let mut fields = Self::zero(spec);
        let two_sigma_sq = 2.0 * sigma_cells * sigma_cells;
        for i in 0..dims.0 {
            for j in 0..dims.1 {
                for k in 0..dims.2 {
                    let di = i as f64 - ci as f64;
                    let dj = j as f64 - cj as f64;
                    let dk = k as f64 - ck as f64;
                    let r_sq = di * di + dj * dj + dk * dk;
                    fields.ez[idx3(dims, i, j, k)] = (-r_sq / two_sigma_sq).exp();
                }
            }
        }
        fields
    }
}
