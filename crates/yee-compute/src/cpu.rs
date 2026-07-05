//! Rayon-parallel FP64 CPU backend.
//!
//! The six kernels are line-for-line ports of the uniform lossless arm of
//! `yee_fdtd::update::{update_h, update_e}` on flat buffers, parallelized by
//! slabbing the outermost `i` index of the *target* component. Every cell's
//! half-step update is independent (H reads only E and vice versa), so the
//! parallel result is required to be **bit-exact** against the scalar
//! reference — gate `compute-001` asserts max |Δ| == 0.0, not a tolerance.
//! Do not reorder the per-cell arithmetic.

use rayon::prelude::*;

use yee_core::units::{EPS0, MU0};

use crate::fields::Fields;
use crate::spec::{FdtdSpec, idx3, len3};

/// Multi-threaded FP64 FDTD stepper (uniform lossless vacuum, PEC box).
#[derive(Debug, Clone)]
pub struct CpuFdtd {
    spec: FdtdSpec,
    fields: Fields,
}

impl CpuFdtd {
    /// Build a stepper from a spec and an initial field state.
    ///
    /// # Panics
    ///
    /// Panics if any buffer length disagrees with the spec's staggered shapes.
    pub fn new(spec: FdtdSpec, fields: Fields) -> Self {
        assert_eq!(fields.ex.len(), len3(spec.ex_dims()), "ex length mismatch");
        assert_eq!(fields.ey.len(), len3(spec.ey_dims()), "ey length mismatch");
        assert_eq!(fields.ez.len(), len3(spec.ez_dims()), "ez length mismatch");
        assert_eq!(fields.hx.len(), len3(spec.hx_dims()), "hx length mismatch");
        assert_eq!(fields.hy.len(), len3(spec.hy_dims()), "hy length mismatch");
        assert_eq!(fields.hz.len(), len3(spec.hz_dims()), "hz length mismatch");
        Self { spec, fields }
    }

    /// The problem description this stepper was built from.
    pub fn spec(&self) -> &FdtdSpec {
        &self.spec
    }

    /// Advance the state by `n` leapfrog steps (H half-step, then E).
    pub fn step_n(&mut self, n: usize) {
        for _ in 0..n {
            self.update_h();
            self.update_e();
        }
    }

    /// Borrow the current field state.
    pub fn fields(&self) -> &Fields {
        &self.fields
    }

    fn update_h(&mut self) {
        let s = self.spec;
        let coeff = s.dt / (MU0 * s.mu_r);
        let eyd = s.ey_dims();
        let ezd = s.ez_dims();
        let exd = s.ex_dims();
        let hxd = s.hx_dims();
        let hyd = s.hy_dims();
        let hzd = s.hz_dims();
        let Fields {
            ex,
            ey,
            ez,
            hx,
            hy,
            hz,
        } = &mut self.fields;

        // ---- H_x: shape [nx+1, ny, nz], full extent ----
        hx.par_chunks_mut(hxd.1 * hxd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                for j in 0..s.ny {
                    for k in 0..s.nz {
                        let dey_dz = (ey[idx3(eyd, i, j, k + 1)] - ey[idx3(eyd, i, j, k)]) / s.dz;
                        let dez_dy = (ez[idx3(ezd, i, j + 1, k)] - ez[idx3(ezd, i, j, k)]) / s.dy;
                        slab[j * hxd.2 + k] += coeff * (dey_dz - dez_dy);
                    }
                }
            });

        // ---- H_y: shape [nx, ny+1, nz], full extent ----
        hy.par_chunks_mut(hyd.1 * hyd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                for j in 0..=s.ny {
                    for k in 0..s.nz {
                        let dez_dx = (ez[idx3(ezd, i + 1, j, k)] - ez[idx3(ezd, i, j, k)]) / s.dx;
                        let dex_dz = (ex[idx3(exd, i, j, k + 1)] - ex[idx3(exd, i, j, k)]) / s.dz;
                        slab[j * hyd.2 + k] += coeff * (dez_dx - dex_dz);
                    }
                }
            });

        // ---- H_z: shape [nx, ny, nz+1], full extent ----
        hz.par_chunks_mut(hzd.1 * hzd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                for j in 0..s.ny {
                    for k in 0..=s.nz {
                        let dex_dy = (ex[idx3(exd, i, j + 1, k)] - ex[idx3(exd, i, j, k)]) / s.dy;
                        let dey_dx = (ey[idx3(eyd, i + 1, j, k)] - ey[idx3(eyd, i, j, k)]) / s.dx;
                        slab[j * hzd.2 + k] += coeff * (dex_dy - dey_dx);
                    }
                }
            });
    }

    fn update_e(&mut self) {
        let s = self.spec;
        let coeff = s.dt / (EPS0 * s.eps_r);
        let exd = s.ex_dims();
        let eyd = s.ey_dims();
        let ezd = s.ez_dims();
        let hxd = s.hx_dims();
        let hyd = s.hy_dims();
        let hzd = s.hz_dims();
        let Fields {
            ex,
            ey,
            ez,
            hx,
            hy,
            hz,
        } = &mut self.fields;

        // ---- E_x: shape [nx, ny+1, nz+1]; interior j ∈ [1, ny), k ∈ [1, nz) ----
        // Outer tangential faces are the PEC box and stay untouched.
        ex.par_chunks_mut(exd.1 * exd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                for j in 1..s.ny {
                    for k in 1..s.nz {
                        let dhz_dy = (hz[idx3(hzd, i, j, k)] - hz[idx3(hzd, i, j - 1, k)]) / s.dy;
                        let dhy_dz = (hy[idx3(hyd, i, j, k)] - hy[idx3(hyd, i, j, k - 1)]) / s.dz;
                        slab[j * exd.2 + k] += coeff * (dhz_dy - dhy_dz);
                    }
                }
            });

        // ---- E_y: shape [nx+1, ny, nz+1]; interior i ∈ [1, nx), k ∈ [1, nz) ----
        ey.par_chunks_mut(eyd.1 * eyd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                if i == 0 || i >= s.nx {
                    return;
                }
                for j in 0..s.ny {
                    for k in 1..s.nz {
                        let dhx_dz = (hx[idx3(hxd, i, j, k)] - hx[idx3(hxd, i, j, k - 1)]) / s.dz;
                        let dhz_dx = (hz[idx3(hzd, i, j, k)] - hz[idx3(hzd, i - 1, j, k)]) / s.dx;
                        slab[j * eyd.2 + k] += coeff * (dhx_dz - dhz_dx);
                    }
                }
            });

        // ---- E_z: shape [nx+1, ny+1, nz]; interior i ∈ [1, nx), j ∈ [1, ny) ----
        ez.par_chunks_mut(ezd.1 * ezd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                if i == 0 || i >= s.nx {
                    return;
                }
                for j in 1..s.ny {
                    for k in 0..s.nz {
                        let dhy_dx = (hy[idx3(hyd, i, j, k)] - hy[idx3(hyd, i - 1, j, k)]) / s.dx;
                        let dhx_dy = (hx[idx3(hxd, i, j, k)] - hx[idx3(hxd, i, j - 1, k)]) / s.dy;
                        slab[j * ezd.2 + k] += coeff * (dhy_dx - dhx_dy);
                    }
                }
            });
    }
}
