//! Flat-buffer CPML state for the CPU backend (E.1).
//!
//! Line-for-line port of `yee_fdtd::cpml::CpmlState` (Roden & Gedney 2000)
//! onto the flat row-major buffers of [`crate::Fields`]. Per-cell arithmetic,
//! branch structure, and branch *order* are identical to the reference —
//! gate `compute-003` asserts bit-exact agreement, so do not "clean up" the
//! math here. Parallelization slabs the outermost `i` index of the written
//! component; the two ψ arrays touched by a pass share that component's
//! shape, so field and ψ slabs zip together and stay disjoint.

use rayon::prelude::*;

use yee_core::units::{EPS0, MU0};

use crate::fields::Fields;
use crate::materials::CpmlConfig;
use crate::spec::{FdtdSpec, idx3, len3};

/// One profile triple `(b, c, κ)`, each of length `npml`.
pub(crate) type ProfileTriple = (Vec<f64>, Vec<f64>, Vec<f64>);

/// σ, κ, α at depth fraction `rho_over_d` (R&G eq. 17 grading).
fn grading_sample(config: &CpmlConfig, rho_over_d: f64) -> (f64, f64, f64) {
    let rho_m = rho_over_d.powi(config.m);
    let sigma = config.sigma_max * rho_m;
    let kappa = 1.0 + (config.kappa_max - 1.0) * rho_m;
    let alpha = config.alpha_max * (1.0 - rho_over_d);
    (sigma, kappa, alpha)
}

/// `(b, c)` from `(σ, κ, α)` and `dt` (R&G eq. 25).
fn finalize_coeffs(sigma: f64, kappa: f64, alpha: f64, dt: f64) -> (f64, f64) {
    let exponent = -(sigma / kappa + alpha) * dt / EPS0;
    let b = exponent.exp();
    let denom = sigma * kappa + kappa * kappa * alpha;
    let c = if denom.abs() > 1.0e-30 {
        sigma * (b - 1.0) / denom
    } else {
        0.0
    };
    (b, c)
}

/// E-grid and H-grid `(b, c, κ)` profiles: E graded at `(d+1)/npml`,
/// H at `(d+0.5)/npml` (half-cell shift). Shared with the GPU upload path.
pub(crate) fn make_profiles(config: &CpmlConfig, dt: f64) -> (ProfileTriple, ProfileTriple) {
    let n = config.npml;
    let mut b_e = vec![1.0; n];
    let mut c_e = vec![0.0; n];
    let mut kappa_e = vec![1.0; n];
    let mut b_h = vec![1.0; n];
    let mut c_h = vec![0.0; n];
    let mut kappa_h = vec![1.0; n];
    for d in 0..n {
        let rho_e = (d as f64 + 1.0) / (n as f64);
        let (sigma, kappa, alpha) = grading_sample(config, rho_e);
        let (b, c) = finalize_coeffs(sigma, kappa, alpha, dt);
        b_e[d] = b;
        c_e[d] = c;
        kappa_e[d] = kappa;

        let rho_h = (d as f64 + 0.5) / (n as f64);
        let (sigma, kappa, alpha) = grading_sample(config, rho_h);
        let (b, c) = finalize_coeffs(sigma, kappa, alpha, dt);
        b_h[d] = b;
        c_h[d] = c;
        kappa_h[d] = kappa;
    }
    ((b_e, c_e, kappa_e), (b_h, c_h, kappa_h))
}

/// PML profile index for absolute index `i` on an axis of length `n`, or
/// `None` outside the PML / on a disabled axis. Mirrors
/// `CpmlState::pml_depth` (the "side" flag there is unused and dropped).
#[inline]
fn pml_depth(faces: [[bool; 2]; 3], npml: usize, axis: usize, i: usize, n: usize) -> Option<usize> {
    if i < npml {
        if !faces[axis][0] {
            return None;
        }
        Some(npml - 1 - i)
    } else if i >= n.saturating_sub(npml) && n >= npml {
        if !faces[axis][1] {
            return None;
        }
        let depth = i - (n - npml);
        if depth < npml { Some(depth) } else { None }
    } else {
        None
    }
}

/// CPML auxiliary state on flat buffers.
#[derive(Debug, Clone)]
pub(crate) struct CpuCpmlState {
    /// ψ_E arrays, order (xy, xz, yx, yz, zx, zy); shapes: E_x, E_x, E_y,
    /// E_y, E_z, E_z.
    psi_e: [Vec<f64>; 6],
    /// ψ_H arrays, same ordering; shapes: H_x, H_x, H_y, H_y, H_z, H_z.
    psi_h: [Vec<f64>; 6],
    b: Vec<f64>,
    c: Vec<f64>,
    kappa: Vec<f64>,
    b_h: Vec<f64>,
    c_h: Vec<f64>,
    kappa_h: Vec<f64>,
    npml: usize,
    faces: [[bool; 2]; 3],
}

impl CpuCpmlState {
    pub(crate) fn new(spec: &FdtdSpec, config: CpmlConfig) -> Self {
        let psi_e = [
            vec![0.0; len3(spec.ex_dims())],
            vec![0.0; len3(spec.ex_dims())],
            vec![0.0; len3(spec.ey_dims())],
            vec![0.0; len3(spec.ey_dims())],
            vec![0.0; len3(spec.ez_dims())],
            vec![0.0; len3(spec.ez_dims())],
        ];
        let psi_h = [
            vec![0.0; len3(spec.hx_dims())],
            vec![0.0; len3(spec.hx_dims())],
            vec![0.0; len3(spec.hy_dims())],
            vec![0.0; len3(spec.hy_dims())],
            vec![0.0; len3(spec.hz_dims())],
            vec![0.0; len3(spec.hz_dims())],
        ];
        let ((b, c, kappa), (b_h, c_h, kappa_h)) = make_profiles(&config, spec.dt);
        Self {
            psi_e,
            psi_h,
            b,
            c,
            kappa,
            b_h,
            c_h,
            kappa_h,
            npml: config.npml,
            faces: config.faces,
        }
    }

    /// CPML correction for the E half-step (call after the bulk E update).
    pub(crate) fn update_e(
        &mut self,
        s: &FdtdSpec,
        fields: &mut Fields,
        eps_r_cells: Option<&[f64]>,
    ) {
        let celld = (s.nx + 1, s.ny + 1, s.nz + 1);
        let coeff_scalar = s.dt / (EPS0 * s.eps_r);
        let exd = s.ex_dims();
        let eyd = s.ey_dims();
        let ezd = s.ez_dims();
        let hxd = s.hx_dims();
        let hyd = s.hy_dims();
        let hzd = s.hz_dims();
        let [p_exy, p_exz, p_eyx, p_eyz, p_ezx, p_ezy] = &mut self.psi_e;
        let (b, c, kappa) = (&self.b, &self.c, &self.kappa);
        let (npml, faces) = (self.npml, self.faces);
        let Fields {
            ex,
            ey,
            ez,
            hx,
            hy,
            hz,
        } = fields;

        // ---- E_x: i ∈ [0, nx), j ∈ [1, ny), k ∈ [1, nz) ----
        let sz = exd.1 * exd.2;
        ex.par_chunks_mut(sz)
            .zip(p_exy.par_chunks_mut(sz))
            .zip(p_exz.par_chunks_mut(sz))
            .enumerate()
            .for_each(|(i, ((ex_s, pxy_s), pxz_s))| {
                for j in 1..s.ny {
                    let dep_y = pml_depth(faces, npml, 1, j, s.ny + 1);
                    for k in 1..s.nz {
                        let dep_z = pml_depth(faces, npml, 2, k, s.nz + 1);
                        if dep_y.is_none() && dep_z.is_none() {
                            continue;
                        }
                        let dhz_dy = (hz[idx3(hzd, i, j, k)] - hz[idx3(hzd, i, j - 1, k)]) / s.dy;
                        let dhy_dz = (hy[idx3(hyd, i, j, k)] - hy[idx3(hyd, i, j, k - 1)]) / s.dz;
                        let coeff = match eps_r_cells {
                            None => coeff_scalar,
                            Some(e) => s.dt / (EPS0 * e[idx3(celld, i, j, k)]),
                        };
                        let off = j * exd.2 + k;
                        if let Some(d) = dep_y {
                            let p = b[d] * pxy_s[off] + c[d] * dhz_dy;
                            pxy_s[off] = p;
                            ex_s[off] += coeff * (p - (1.0 - 1.0 / kappa[d]) * dhz_dy);
                        }
                        if let Some(d) = dep_z {
                            let p = b[d] * pxz_s[off] + c[d] * dhy_dz;
                            pxz_s[off] = p;
                            ex_s[off] -= coeff * (p - (1.0 - 1.0 / kappa[d]) * dhy_dz);
                        }
                    }
                }
            });

        // ---- E_y: i ∈ [1, nx), j ∈ [0, ny), k ∈ [1, nz) ----
        let sz = eyd.1 * eyd.2;
        ey.par_chunks_mut(sz)
            .zip(p_eyx.par_chunks_mut(sz))
            .zip(p_eyz.par_chunks_mut(sz))
            .enumerate()
            .for_each(|(i, ((ey_s, pyx_s), pyz_s))| {
                if i == 0 || i >= s.nx {
                    return;
                }
                let dep_x = pml_depth(faces, npml, 0, i, s.nx + 1);
                for j in 0..s.ny {
                    for k in 1..s.nz {
                        let dep_z = pml_depth(faces, npml, 2, k, s.nz + 1);
                        if dep_x.is_none() && dep_z.is_none() {
                            continue;
                        }
                        let dhx_dz = (hx[idx3(hxd, i, j, k)] - hx[idx3(hxd, i, j, k - 1)]) / s.dz;
                        let dhz_dx = (hz[idx3(hzd, i, j, k)] - hz[idx3(hzd, i - 1, j, k)]) / s.dx;
                        let coeff = match eps_r_cells {
                            None => coeff_scalar,
                            Some(e) => s.dt / (EPS0 * e[idx3(celld, i, j, k)]),
                        };
                        let off = j * eyd.2 + k;
                        if let Some(d) = dep_z {
                            let p = b[d] * pyz_s[off] + c[d] * dhx_dz;
                            pyz_s[off] = p;
                            ey_s[off] += coeff * (p - (1.0 - 1.0 / kappa[d]) * dhx_dz);
                        }
                        if let Some(d) = dep_x {
                            let p = b[d] * pyx_s[off] + c[d] * dhz_dx;
                            pyx_s[off] = p;
                            ey_s[off] -= coeff * (p - (1.0 - 1.0 / kappa[d]) * dhz_dx);
                        }
                    }
                }
            });

        // ---- E_z: i ∈ [1, nx), j ∈ [1, ny), k ∈ [0, nz) ----
        let sz = ezd.1 * ezd.2;
        ez.par_chunks_mut(sz)
            .zip(p_ezx.par_chunks_mut(sz))
            .zip(p_ezy.par_chunks_mut(sz))
            .enumerate()
            .for_each(|(i, ((ez_s, pzx_s), pzy_s))| {
                if i == 0 || i >= s.nx {
                    return;
                }
                let dep_x = pml_depth(faces, npml, 0, i, s.nx + 1);
                for j in 1..s.ny {
                    let dep_y = pml_depth(faces, npml, 1, j, s.ny + 1);
                    if dep_x.is_none() && dep_y.is_none() {
                        continue;
                    }
                    for k in 0..s.nz {
                        let dhy_dx = (hy[idx3(hyd, i, j, k)] - hy[idx3(hyd, i - 1, j, k)]) / s.dx;
                        let dhx_dy = (hx[idx3(hxd, i, j, k)] - hx[idx3(hxd, i, j - 1, k)]) / s.dy;
                        let coeff = match eps_r_cells {
                            None => coeff_scalar,
                            Some(e) => s.dt / (EPS0 * e[idx3(celld, i, j, k)]),
                        };
                        let off = j * ezd.2 + k;
                        if let Some(d) = dep_x {
                            let p = b[d] * pzx_s[off] + c[d] * dhy_dx;
                            pzx_s[off] = p;
                            ez_s[off] += coeff * (p - (1.0 - 1.0 / kappa[d]) * dhy_dx);
                        }
                        if let Some(d) = dep_y {
                            let p = b[d] * pzy_s[off] + c[d] * dhx_dy;
                            pzy_s[off] = p;
                            ez_s[off] -= coeff * (p - (1.0 - 1.0 / kappa[d]) * dhx_dy);
                        }
                    }
                }
            });
    }

    /// CPML correction for the H half-step (call after the bulk H update).
    pub(crate) fn update_h(
        &mut self,
        s: &FdtdSpec,
        fields: &mut Fields,
        mu_r_cells: Option<&[f64]>,
    ) {
        let celld = (s.nx + 1, s.ny + 1, s.nz + 1);
        let coeff_scalar = s.dt / (MU0 * s.mu_r);
        let exd = s.ex_dims();
        let eyd = s.ey_dims();
        let ezd = s.ez_dims();
        let hxd = s.hx_dims();
        let hyd = s.hy_dims();
        let hzd = s.hz_dims();
        let [p_hxy, p_hxz, p_hyx, p_hyz, p_hzx, p_hzy] = &mut self.psi_h;
        let (b_h, c_h, kappa_h) = (&self.b_h, &self.c_h, &self.kappa_h);
        let (npml, faces) = (self.npml, self.faces);
        let Fields {
            ex,
            ey,
            ez,
            hx,
            hy,
            hz,
        } = fields;

        // ---- H_x: i ∈ [0, nx], j ∈ [0, ny), k ∈ [0, nz) ----
        let sz = hxd.1 * hxd.2;
        hx.par_chunks_mut(sz)
            .zip(p_hxy.par_chunks_mut(sz))
            .zip(p_hxz.par_chunks_mut(sz))
            .enumerate()
            .for_each(|(i, ((hx_s, pxy_s), pxz_s))| {
                for j in 0..s.ny {
                    let dep_y = pml_depth(faces, npml, 1, j, s.ny);
                    for k in 0..s.nz {
                        let dep_z = pml_depth(faces, npml, 2, k, s.nz);
                        if dep_y.is_none() && dep_z.is_none() {
                            continue;
                        }
                        let dey_dz = (ey[idx3(eyd, i, j, k + 1)] - ey[idx3(eyd, i, j, k)]) / s.dz;
                        let dez_dy = (ez[idx3(ezd, i, j + 1, k)] - ez[idx3(ezd, i, j, k)]) / s.dy;
                        let coeff = match mu_r_cells {
                            None => coeff_scalar,
                            Some(m) => s.dt / (MU0 * m[idx3(celld, i, j, k)]),
                        };
                        let off = j * hxd.2 + k;
                        if let Some(d) = dep_z {
                            let p = b_h[d] * pxz_s[off] + c_h[d] * dey_dz;
                            pxz_s[off] = p;
                            hx_s[off] += coeff * (p - (1.0 - 1.0 / kappa_h[d]) * dey_dz);
                        }
                        if let Some(d) = dep_y {
                            let p = b_h[d] * pxy_s[off] + c_h[d] * dez_dy;
                            pxy_s[off] = p;
                            hx_s[off] -= coeff * (p - (1.0 - 1.0 / kappa_h[d]) * dez_dy);
                        }
                    }
                }
            });

        // ---- H_y: i ∈ [0, nx), j ∈ [0, ny], k ∈ [0, nz) ----
        let sz = hyd.1 * hyd.2;
        hy.par_chunks_mut(sz)
            .zip(p_hyx.par_chunks_mut(sz))
            .zip(p_hyz.par_chunks_mut(sz))
            .enumerate()
            .for_each(|(i, ((hy_s, pyx_s), pyz_s))| {
                let dep_x = pml_depth(faces, npml, 0, i, s.nx);
                for j in 0..=s.ny {
                    for k in 0..s.nz {
                        let dep_z = pml_depth(faces, npml, 2, k, s.nz);
                        if dep_x.is_none() && dep_z.is_none() {
                            continue;
                        }
                        let dez_dx = (ez[idx3(ezd, i + 1, j, k)] - ez[idx3(ezd, i, j, k)]) / s.dx;
                        let dex_dz = (ex[idx3(exd, i, j, k + 1)] - ex[idx3(exd, i, j, k)]) / s.dz;
                        let coeff = match mu_r_cells {
                            None => coeff_scalar,
                            Some(m) => s.dt / (MU0 * m[idx3(celld, i, j, k)]),
                        };
                        let off = j * hyd.2 + k;
                        if let Some(d) = dep_x {
                            let p = b_h[d] * pyx_s[off] + c_h[d] * dez_dx;
                            pyx_s[off] = p;
                            hy_s[off] += coeff * (p - (1.0 - 1.0 / kappa_h[d]) * dez_dx);
                        }
                        if let Some(d) = dep_z {
                            let p = b_h[d] * pyz_s[off] + c_h[d] * dex_dz;
                            pyz_s[off] = p;
                            hy_s[off] -= coeff * (p - (1.0 - 1.0 / kappa_h[d]) * dex_dz);
                        }
                    }
                }
            });

        // ---- H_z: i ∈ [0, nx), j ∈ [0, ny), k ∈ [0, nz] ----
        let sz = hzd.1 * hzd.2;
        hz.par_chunks_mut(sz)
            .zip(p_hzx.par_chunks_mut(sz))
            .zip(p_hzy.par_chunks_mut(sz))
            .enumerate()
            .for_each(|(i, ((hz_s, pzx_s), pzy_s))| {
                let dep_x = pml_depth(faces, npml, 0, i, s.nx);
                for j in 0..s.ny {
                    let dep_y = pml_depth(faces, npml, 1, j, s.ny);
                    if dep_x.is_none() && dep_y.is_none() {
                        continue;
                    }
                    for k in 0..=s.nz {
                        let dex_dy = (ex[idx3(exd, i, j + 1, k)] - ex[idx3(exd, i, j, k)]) / s.dy;
                        let dey_dx = (ey[idx3(eyd, i + 1, j, k)] - ey[idx3(eyd, i, j, k)]) / s.dx;
                        let coeff = match mu_r_cells {
                            None => coeff_scalar,
                            Some(m) => s.dt / (MU0 * m[idx3(celld, i, j, k)]),
                        };
                        let off = j * hzd.2 + k;
                        if let Some(d) = dep_y {
                            let p = b_h[d] * pzy_s[off] + c_h[d] * dex_dy;
                            pzy_s[off] = p;
                            hz_s[off] += coeff * (p - (1.0 - 1.0 / kappa_h[d]) * dex_dy);
                        }
                        if let Some(d) = dep_x {
                            let p = b_h[d] * pzx_s[off] + c_h[d] * dey_dx;
                            pzx_s[off] = p;
                            hz_s[off] -= coeff * (p - (1.0 - 1.0 / kappa_h[d]) * dey_dx);
                        }
                    }
                }
            });
    }
}
