//! Rayon-parallel FP64 CPU backend.
//!
//! The kernels are line-for-line ports of `yee_fdtd::update::{update_h,
//! update_e}` (including the per-cell ε_r / μ_r / σ arms) on flat buffers,
//! parallelized by slabbing the outermost `i` index of the *target*
//! component. Every cell's half-step update is independent (H reads only E
//! and vice versa), so the parallel result is required to be **bit-exact**
//! against the scalar reference — gates `compute-001` (uniform vacuum) and
//! `compute-003` (heterogeneous + CPML + masks) assert max |Δ| == 0.0, not a
//! tolerance. Do not reorder the per-cell arithmetic.
//!
//! The full-step orchestration mirrors `yee_fdtd::WalkingSkeletonSolver`:
//! `update_h` → boundary-H (CPML or legacy PEC clamp) → optional source →
//! `update_e` → boundary-E → interior PEC mask → clock.

use rayon::prelude::*;

use yee_core::units::{EPS0, MU0};

use crate::cpml::CpuCpmlState;
use crate::fields::Fields;
use crate::materials::{Boundary, Materials};
use crate::spec::{FdtdSpec, idx3, len3};

/// Lossy CA/CB coefficients for the Yee E-update (Taflove §3.7), identical
/// to `yee_fdtd::update::ca_cb`.
#[inline]
fn ca_cb(eps_r: f64, sigma: f64, dt: f64) -> (f64, f64) {
    let denom = 2.0 * EPS0 * eps_r + sigma * dt;
    let ca = (2.0 * EPS0 * eps_r - sigma * dt) / denom;
    let cb = 2.0 * dt / denom;
    (ca, cb)
}

/// Multi-threaded FP64 FDTD stepper.
#[derive(Debug, Clone)]
pub struct CpuFdtd {
    spec: FdtdSpec,
    fields: Fields,
    materials: Materials,
    cpml: Option<CpuCpmlState>,
    pec_box: bool,
    step: u64,
}

impl CpuFdtd {
    /// Build a uniform-vacuum stepper with no boundary phase — the raw E.0
    /// kernel semantics (outer tangential E faces are never written).
    ///
    /// # Panics
    ///
    /// Panics if any buffer length disagrees with the spec's staggered shapes.
    pub fn new(spec: FdtdSpec, fields: Fields) -> Self {
        Self::with_config(spec, fields, Materials::default(), Boundary::None)
    }

    /// Build a stepper with per-cell materials / masks and an outer-boundary
    /// treatment (E.1).
    ///
    /// # Panics
    ///
    /// Panics if any field, material, or mask buffer length disagrees with
    /// the spec's staggered shapes.
    pub fn with_config(
        spec: FdtdSpec,
        fields: Fields,
        materials: Materials,
        boundary: Boundary,
    ) -> Self {
        assert_eq!(fields.ex.len(), len3(spec.ex_dims()), "ex length mismatch");
        assert_eq!(fields.ey.len(), len3(spec.ey_dims()), "ey length mismatch");
        assert_eq!(fields.ez.len(), len3(spec.ez_dims()), "ez length mismatch");
        assert_eq!(fields.hx.len(), len3(spec.hx_dims()), "hx length mismatch");
        assert_eq!(fields.hy.len(), len3(spec.hy_dims()), "hy length mismatch");
        assert_eq!(fields.hz.len(), len3(spec.hz_dims()), "hz length mismatch");
        materials.validate(&spec);
        let (cpml, pec_box) = match boundary {
            Boundary::None => (None, false),
            Boundary::PecBox => (None, true),
            Boundary::Cpml(config) => (Some(CpuCpmlState::new(&spec, config)), false),
        };
        Self {
            spec,
            fields,
            materials,
            cpml,
            pec_box,
            step: 0,
        }
    }

    /// The problem description this stepper was built from.
    pub fn spec(&self) -> &FdtdSpec {
        &self.spec
    }

    /// Simulation time at the start of the next step (seconds).
    pub fn current_time(&self) -> f64 {
        self.step as f64 * self.spec.dt
    }

    /// Advance the state by `n` full leapfrog steps (no source).
    pub fn step_n(&mut self, n: usize) {
        for _ in 0..n {
            self.update_h();
            self.boundary_h();
            self.update_e();
            self.boundary_e();
            self.step += 1;
        }
    }

    /// One full step injecting a Gaussian-in-time soft pulse on `E_z` at
    /// `source`, sampled at the current simulation time — identical timing
    /// and amplitude to `WalkingSkeletonSolver::step_with_source` /
    /// `sources::gaussian_pulse_ez`.
    ///
    /// # Panics
    ///
    /// Panics if `sigma` is non-positive or `source` is out of bounds.
    pub fn step_with_gaussian_ez(&mut self, source: (usize, usize, usize), t0: f64, sigma: f64) {
        assert!(
            sigma > 0.0 && sigma.is_finite(),
            "gaussian sigma must be positive and finite"
        );
        let t = self.current_time();
        self.update_h();
        self.boundary_h();
        let arg = (t - t0) / sigma;
        let amplitude = (-arg * arg).exp();
        let ezd = self.spec.ez_dims();
        self.fields.ez[idx3(ezd, source.0, source.1, source.2)] += amplitude;
        self.update_e();
        self.boundary_e();
        self.step += 1;
    }

    /// Borrow the current field state.
    pub fn fields(&self) -> &Fields {
        &self.fields
    }

    /// Outer-boundary phase after the H half-step: CPML auxiliary update,
    /// or the legacy PEC clamp in [`Boundary::PecBox`] mode.
    fn boundary_h(&mut self) {
        if let Some(cpml) = self.cpml.as_mut() {
            cpml.update_h(
                &self.spec,
                &mut self.fields,
                self.materials.mu_r_cells.as_deref(),
            );
        } else if self.pec_box {
            apply_pec_box(&self.spec, &mut self.fields);
        }
    }

    /// Outer-boundary phase after the E half-step, then the interior PEC
    /// mask (the mask is the final word for the step, as in the reference).
    fn boundary_e(&mut self) {
        if let Some(cpml) = self.cpml.as_mut() {
            cpml.update_e(
                &self.spec,
                &mut self.fields,
                self.materials.eps_r_cells.as_deref(),
            );
        } else if self.pec_box {
            apply_pec_box(&self.spec, &mut self.fields);
        }
        self.apply_pec_mask();
    }

    /// Clamp masked E cells to zero (no-op when no masks are attached),
    /// mirroring `YeeGrid::apply_pec_mask`.
    fn apply_pec_mask(&mut self) {
        for (field, mask) in [
            (&mut self.fields.ex, &self.materials.pec_mask_ex),
            (&mut self.fields.ey, &self.materials.pec_mask_ey),
            (&mut self.fields.ez, &self.materials.pec_mask_ez),
        ] {
            if let Some(mask) = mask {
                for (e, &m) in field.iter_mut().zip(mask.iter()) {
                    if m {
                        *e = 0.0;
                    }
                }
            }
        }
    }

    fn update_h(&mut self) {
        let s = self.spec;
        let celld = (s.nx + 1, s.ny + 1, s.nz + 1);
        let coeff_scalar = s.dt / (MU0 * s.mu_r);
        let mu_r_cells = self.materials.mu_r_cells.as_deref();
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

        // ---- H_x: shape [nx+1, ny, nz], full extent ----
        hx.par_chunks_mut(hxd.1 * hxd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                for j in 0..s.ny {
                    for k in 0..s.nz {
                        let dey_dz = (ey[idx3(eyd, i, j, k + 1)] - ey[idx3(eyd, i, j, k)]) / s.dz;
                        let dez_dy = (ez[idx3(ezd, i, j + 1, k)] - ez[idx3(ezd, i, j, k)]) / s.dy;
                        let coeff = match mu_r_cells {
                            None => coeff_scalar,
                            Some(m) => s.dt / (MU0 * m[idx3(celld, i, j, k)]),
                        };
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
                        let coeff = match mu_r_cells {
                            None => coeff_scalar,
                            Some(m) => s.dt / (MU0 * m[idx3(celld, i, j, k)]),
                        };
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
                        let coeff = match mu_r_cells {
                            None => coeff_scalar,
                            Some(m) => s.dt / (MU0 * m[idx3(celld, i, j, k)]),
                        };
                        slab[j * hzd.2 + k] += coeff * (dex_dy - dey_dx);
                    }
                }
            });
    }

    fn update_e(&mut self) {
        let s = self.spec;
        let celld = (s.nx + 1, s.ny + 1, s.nz + 1);
        let coeff_scalar = s.dt / (EPS0 * s.eps_r);
        let eps_r_cells = self.materials.eps_r_cells.as_deref();
        let sigma_cells = self.materials.sigma_cells.as_deref();
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

        // Per-cell E update body, shared by the three components. `e` is the
        // target cell, `curl_h` its curl term — the match structure and
        // operation order are the reference's, verbatim.
        let update_cell = |e: &mut f64, curl_h: f64, cell: usize| match sigma_cells {
            None => {
                let coeff = match eps_r_cells {
                    None => coeff_scalar,
                    Some(eps) => s.dt / (EPS0 * eps[cell]),
                };
                *e += coeff * curl_h;
            }
            Some(sig) => {
                let eps_r = match eps_r_cells {
                    None => s.eps_r,
                    Some(eps) => eps[cell],
                };
                let (ca, cb) = ca_cb(eps_r, sig[cell], s.dt);
                *e = ca * *e + cb * curl_h;
            }
        };

        // ---- E_x: shape [nx, ny+1, nz+1]; interior j ∈ [1, ny), k ∈ [1, nz) ----
        // Outer tangential faces are managed by the boundary phase.
        ex.par_chunks_mut(exd.1 * exd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                for j in 1..s.ny {
                    for k in 1..s.nz {
                        let dhz_dy = (hz[idx3(hzd, i, j, k)] - hz[idx3(hzd, i, j - 1, k)]) / s.dy;
                        let dhy_dz = (hy[idx3(hyd, i, j, k)] - hy[idx3(hyd, i, j, k - 1)]) / s.dz;
                        update_cell(
                            &mut slab[j * exd.2 + k],
                            dhz_dy - dhy_dz,
                            idx3(celld, i, j, k),
                        );
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
                        update_cell(
                            &mut slab[j * eyd.2 + k],
                            dhx_dz - dhz_dx,
                            idx3(celld, i, j, k),
                        );
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
                        update_cell(
                            &mut slab[j * ezd.2 + k],
                            dhy_dx - dhx_dy,
                            idx3(celld, i, j, k),
                        );
                    }
                }
            });
    }
}

/// Zero the tangential E field on all six outer faces (legacy reflecting
/// PEC box), mirroring `yee_fdtd::boundary::apply_pec`.
pub(crate) fn apply_pec_box(s: &FdtdSpec, fields: &mut Fields) {
    let exd = s.ex_dims();
    let eyd = s.ey_dims();
    let ezd = s.ez_dims();

    // x = 0 and x = nx faces: clamp E_y and E_z.
    for j in 0..s.ny {
        for k in 0..=s.nz {
            fields.ey[idx3(eyd, 0, j, k)] = 0.0;
            fields.ey[idx3(eyd, s.nx, j, k)] = 0.0;
        }
    }
    for j in 0..=s.ny {
        for k in 0..s.nz {
            fields.ez[idx3(ezd, 0, j, k)] = 0.0;
            fields.ez[idx3(ezd, s.nx, j, k)] = 0.0;
        }
    }

    // y = 0 and y = ny faces: clamp E_x and E_z.
    for i in 0..s.nx {
        for k in 0..=s.nz {
            fields.ex[idx3(exd, i, 0, k)] = 0.0;
            fields.ex[idx3(exd, i, s.ny, k)] = 0.0;
        }
    }
    for i in 0..=s.nx {
        for k in 0..s.nz {
            fields.ez[idx3(ezd, i, 0, k)] = 0.0;
            fields.ez[idx3(ezd, i, s.ny, k)] = 0.0;
        }
    }

    // z = 0 and z = nz faces: clamp E_x and E_y.
    for i in 0..s.nx {
        for j in 0..=s.ny {
            fields.ex[idx3(exd, i, j, 0)] = 0.0;
            fields.ex[idx3(exd, i, j, s.nz)] = 0.0;
        }
    }
    for i in 0..=s.nx {
        for j in 0..s.ny {
            fields.ey[idx3(eyd, i, j, 0)] = 0.0;
            fields.ey[idx3(eyd, i, j, s.nz)] = 0.0;
        }
    }
}
