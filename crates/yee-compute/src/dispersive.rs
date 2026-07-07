//! ADE dispersive materials (Drude / Lorentz / Debye) on the CPU backend
//! (E.5c). Verbatim port of `yee_fdtd::dispersive` — per-cell arithmetic and
//! op order identical to the reference `ade_step`; gate `compute-011`
//! asserts bit-exact agreement, so do not restructure the math.
//!
//! Storage follows the reference walking skeleton: nine full-grid
//! `[nx+1, ny+1, nz+1]` auxiliary arrays (`jp_*`, `jp_*_prev`, `p_*`),
//! each component addressed by its own `(i, j, k)`.

use rayon::prelude::*;

use yee_core::units::EPS0;

use crate::fields::Fields;
use crate::spec::{FdtdSpec, idx3};

/// Single-pole dispersion models, mirroring `yee_fdtd::material::Material`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DispersiveMaterial {
    /// Non-dispersive vacuum cell (`E += (Δt/ε₀)·curl H`).
    Vacuum,
    /// Drude metal: `ε(ω) = ε∞ − ω_p²/(ω² + jγω)`.
    Drude {
        /// High-frequency permittivity ε∞.
        eps_inf: f64,
        /// Plasma frequency ω_p (rad/s).
        omega_p: f64,
        /// Collision rate γ (1/s).
        gamma: f64,
    },
    /// Lorentz pole: `ε(ω) = ε∞ + Δε·ω₀²/(ω₀² + 2jδω − ω²)`.
    Lorentz {
        /// High-frequency permittivity ε∞.
        eps_inf: f64,
        /// Pole strength Δε.
        delta_eps: f64,
        /// Resonance ω₀ (rad/s).
        omega_0: f64,
        /// Damping δ (1/s).
        delta: f64,
    },
    /// Debye relaxation: `ε(ω) = ε∞ + Δε/(1 + jωτ)`.
    Debye {
        /// High-frequency permittivity ε∞.
        eps_inf: f64,
        /// Relaxation strength Δε.
        delta_eps: f64,
        /// Relaxation time τ (s).
        tau: f64,
    },
}

/// Per-cell dispersive material map, `[nx+1, ny+1, nz+1]` row-major
/// (the `yee_fdtd::MaterialMap` convention: all three E components in a
/// primary cell see the same material).
#[derive(Debug, Clone)]
pub struct DispersiveMap {
    /// Flat per-cell materials.
    pub cells: Vec<DispersiveMaterial>,
}

impl DispersiveMap {
    /// All-vacuum map for `spec`.
    pub fn vacuum(spec: &FdtdSpec) -> Self {
        let n = (spec.nx + 1) * (spec.ny + 1) * (spec.nz + 1);
        Self {
            cells: vec![DispersiveMaterial::Vacuum; n],
        }
    }

    /// Assign `mat` on the half-open index box
    /// `[i0, i1) × [j0, j1) × [k0, k1)`.
    #[allow(clippy::too_many_arguments)]
    pub fn set_box(
        &mut self,
        spec: &FdtdSpec,
        i0: usize,
        i1: usize,
        j0: usize,
        j1: usize,
        k0: usize,
        k1: usize,
        mat: DispersiveMaterial,
    ) {
        let celld = (spec.nx + 1, spec.ny + 1, spec.nz + 1);
        for i in i0..i1.min(celld.0) {
            for j in j0..j1.min(celld.1) {
                for k in k0..k1.min(celld.2) {
                    self.cells[idx3(celld, i, j, k)] = mat;
                }
            }
        }
    }

    pub(crate) fn validate(&self, spec: &FdtdSpec) {
        let n = (spec.nx + 1) * (spec.ny + 1) * (spec.nz + 1);
        assert_eq!(self.cells.len(), n, "dispersive map length mismatch");
    }
}

/// Auxiliary ADE state on flat buffers (reference `DispersiveState` layout).
#[derive(Debug, Clone)]
pub(crate) struct CpuDispersiveState {
    jp: [Vec<f64>; 3],
    jp_prev: [Vec<f64>; 3],
    p: [Vec<f64>; 3],
}

impl CpuDispersiveState {
    pub(crate) fn new(spec: &FdtdSpec) -> Self {
        let n = (spec.nx + 1) * (spec.ny + 1) * (spec.nz + 1);
        let z = || vec![0.0; n];
        Self {
            jp: [z(), z(), z()],
            jp_prev: [z(), z(), z()],
            p: [z(), z(), z()],
        }
    }

    /// Fused ADE + E update, replacing the standard E half-step — the
    /// reference `update_e_with_dispersion` on flat buffers, slab-parallel
    /// over `i` (each component's aux arrays are written at its own
    /// `(i, j, k)`, so field and aux slabs zip disjointly).
    pub(crate) fn update_e(&mut self, s: &FdtdSpec, fields: &mut Fields, map: &DispersiveMap) {
        let celld = (s.nx + 1, s.ny + 1, s.nz + 1);
        let cell_sz = celld.1 * celld.2;
        let exd = s.ex_dims();
        let eyd = s.ey_dims();
        let ezd = s.ez_dims();
        let hxd = s.hx_dims();
        let hyd = s.hy_dims();
        let hzd = s.hz_dims();
        let mats = &map.cells;
        let [jp_x, jp_y, jp_z] = &mut self.jp;
        let [jpp_x, jpp_y, jpp_z] = &mut self.jp_prev;
        let [p_x, p_y, p_z] = &mut self.p;
        let Fields {
            ex,
            ey,
            ez,
            hx,
            hy,
            hz,
        } = fields;

        // ---- E_x: interior j ∈ [1, ny), k ∈ [1, nz) ----
        ex.par_chunks_mut(exd.1 * exd.2)
            .zip(jp_x.par_chunks_mut(cell_sz))
            .zip(jpp_x.par_chunks_mut(cell_sz))
            .zip(p_x.par_chunks_mut(cell_sz))
            .enumerate()
            .for_each(|(i, (((ex_s, jp_s), jpp_s), p_s))| {
                for j in 1..s.ny {
                    for k in 1..s.nz {
                        let dhz_dy = (hz[idx3(hzd, i, j, k)] - hz[idx3(hzd, i, j - 1, k)]) / s.dy;
                        let dhy_dz = (hy[idx3(hyd, i, j, k)] - hy[idx3(hyd, i, j, k - 1)]) / s.dz;
                        let curl_h = dhz_dy - dhy_dz;
                        let cell = j * celld.2 + k;
                        let e_off = j * exd.2 + k;
                        let (e_new, jp_new, jp_prev_new, p_new) = ade_step(
                            mats[i * cell_sz + cell],
                            ex_s[e_off],
                            curl_h,
                            jp_s[cell],
                            jpp_s[cell],
                            p_s[cell],
                            s.dt,
                        );
                        ex_s[e_off] = e_new;
                        jpp_s[cell] = jp_prev_new;
                        jp_s[cell] = jp_new;
                        p_s[cell] = p_new;
                    }
                }
            });

        // ---- E_y: interior i ∈ [1, nx), k ∈ [1, nz) ----
        ey.par_chunks_mut(eyd.1 * eyd.2)
            .zip(jp_y.par_chunks_mut(cell_sz))
            .zip(jpp_y.par_chunks_mut(cell_sz))
            .zip(p_y.par_chunks_mut(cell_sz))
            .enumerate()
            .for_each(|(i, (((ey_s, jp_s), jpp_s), p_s))| {
                if i == 0 || i >= s.nx {
                    return;
                }
                for j in 0..s.ny {
                    for k in 1..s.nz {
                        let dhx_dz = (hx[idx3(hxd, i, j, k)] - hx[idx3(hxd, i, j, k - 1)]) / s.dz;
                        let dhz_dx = (hz[idx3(hzd, i, j, k)] - hz[idx3(hzd, i - 1, j, k)]) / s.dx;
                        let curl_h = dhx_dz - dhz_dx;
                        let cell = j * celld.2 + k;
                        let e_off = j * eyd.2 + k;
                        let (e_new, jp_new, jp_prev_new, p_new) = ade_step(
                            mats[i * cell_sz + cell],
                            ey_s[e_off],
                            curl_h,
                            jp_s[cell],
                            jpp_s[cell],
                            p_s[cell],
                            s.dt,
                        );
                        ey_s[e_off] = e_new;
                        jpp_s[cell] = jp_prev_new;
                        jp_s[cell] = jp_new;
                        p_s[cell] = p_new;
                    }
                }
            });

        // ---- E_z: interior i ∈ [1, nx), j ∈ [1, ny) ----
        ez.par_chunks_mut(ezd.1 * ezd.2)
            .zip(jp_z.par_chunks_mut(cell_sz))
            .zip(jpp_z.par_chunks_mut(cell_sz))
            .zip(p_z.par_chunks_mut(cell_sz))
            .enumerate()
            .for_each(|(i, (((ez_s, jp_s), jpp_s), p_s))| {
                if i == 0 || i >= s.nx {
                    return;
                }
                for j in 1..s.ny {
                    for k in 0..s.nz {
                        let dhy_dx = (hy[idx3(hyd, i, j, k)] - hy[idx3(hyd, i - 1, j, k)]) / s.dx;
                        let dhx_dy = (hx[idx3(hxd, i, j, k)] - hx[idx3(hxd, i, j - 1, k)]) / s.dy;
                        let curl_h = dhy_dx - dhx_dy;
                        let cell = j * celld.2 + k;
                        let e_off = j * ezd.2 + k;
                        let (e_new, jp_new, jp_prev_new, p_new) = ade_step(
                            mats[i * cell_sz + cell],
                            ez_s[e_off],
                            curl_h,
                            jp_s[cell],
                            jpp_s[cell],
                            p_s[cell],
                            s.dt,
                        );
                        ez_s[e_off] = e_new;
                        jpp_s[cell] = jp_prev_new;
                        jp_s[cell] = jp_new;
                        p_s[cell] = p_new;
                    }
                }
            });
    }
}

/// Per-cell ADE step — verbatim `yee_fdtd::dispersive::ade_step`.
#[inline]
fn ade_step(
    mat: DispersiveMaterial,
    e_old: f64,
    curl_h: f64,
    jp_old: f64,
    jp_prev_old: f64,
    p_old: f64,
    dt: f64,
) -> (f64, f64, f64, f64) {
    match mat {
        DispersiveMaterial::Vacuum => {
            let e_new = e_old + (dt / EPS0) * curl_h;
            (e_new, jp_old, jp_prev_old, p_old)
        }
        DispersiveMaterial::Drude {
            eps_inf,
            omega_p,
            gamma,
        } => {
            let half = 0.5 * gamma * dt;
            let alpha = (1.0 - half) / (1.0 + half);
            let beta = (EPS0 * omega_p * omega_p * dt) / (1.0 + half);
            let jp_new = alpha * jp_old + beta * e_old;
            let j_avg = 0.5 * (jp_new + jp_old);
            let e_new = e_old + (dt / (eps_inf * EPS0)) * (curl_h - j_avg);
            (e_new, jp_new, jp_old, p_old)
        }
        DispersiveMaterial::Lorentz {
            eps_inf,
            delta_eps,
            omega_0,
            delta,
        } => {
            let denom = 1.0 + delta * dt;
            let alpha_l = (2.0 - omega_0 * omega_0 * dt * dt) / denom;
            let beta_l = (delta * dt - 1.0) / denom;
            let gamma_l = (EPS0 * delta_eps * omega_0 * omega_0 * dt * dt) / denom;
            let p_n = jp_old;
            let p_nm1 = jp_prev_old;
            let p_np1 = alpha_l * p_n + beta_l * p_nm1 + gamma_l * e_old;
            let e_new = e_old + (dt / (eps_inf * EPS0)) * curl_h - (p_np1 - p_n) / (eps_inf * EPS0);
            (e_new, p_np1, p_n, p_old)
        }
        DispersiveMaterial::Debye {
            eps_inf,
            delta_eps,
            tau,
        } => {
            let p_new = p_old + (dt / tau) * (EPS0 * delta_eps * e_old - p_old);
            let e_new =
                e_old + (dt / (eps_inf * EPS0)) * curl_h - (p_new - p_old) / (eps_inf * EPS0);
            (e_new, jp_old, jp_prev_old, p_new)
        }
    }
}

/// Unified ADE coefficients for the GPU path (tolerance-gated, so the
/// algebraic flattening below is allowed; the CPU path stays verbatim):
///
/// ```text
/// aux1' = c0·aux1 + c1·aux2 + c2·E
/// aux2' = aux1
/// E'    = E + ce·curl_H + q·(aux1' + s·aux1)
/// ```
///
/// The `(aux1' + s·aux1)` grouping is load-bearing: Lorentz/Debye subtract
/// two nearly-equal polarizations (`s = −1`), and computing the difference
/// *before* scaling preserves the cancellation in FP32 (a flattened
/// `d_new·aux1' + d_old·aux1` loses ~10 bits and fails the parity gate).
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct AdeCoeffs {
    pub ce: f64,
    pub c0: f64,
    pub c1: f64,
    pub c2: f64,
    pub q: f64,
    pub s: f64,
}

pub(crate) fn ade_coeffs(mat: DispersiveMaterial, dt: f64) -> AdeCoeffs {
    match mat {
        DispersiveMaterial::Vacuum => AdeCoeffs {
            ce: dt / EPS0,
            ..Default::default()
        },
        DispersiveMaterial::Drude {
            eps_inf,
            omega_p,
            gamma,
        } => {
            let half = 0.5 * gamma * dt;
            let ce = dt / (eps_inf * EPS0);
            AdeCoeffs {
                ce,
                c0: (1.0 - half) / (1.0 + half),
                c1: 0.0,
                c2: (EPS0 * omega_p * omega_p * dt) / (1.0 + half),
                q: -0.5 * ce,
                s: 1.0,
            }
        }
        DispersiveMaterial::Lorentz {
            eps_inf,
            delta_eps,
            omega_0,
            delta,
        } => {
            let denom = 1.0 + delta * dt;
            let k = 1.0 / (eps_inf * EPS0);
            AdeCoeffs {
                ce: dt * k,
                c0: (2.0 - omega_0 * omega_0 * dt * dt) / denom,
                c1: (delta * dt - 1.0) / denom,
                c2: (EPS0 * delta_eps * omega_0 * omega_0 * dt * dt) / denom,
                q: -k,
                s: -1.0,
            }
        }
        DispersiveMaterial::Debye {
            eps_inf,
            delta_eps,
            tau,
        } => {
            let k = 1.0 / (eps_inf * EPS0);
            AdeCoeffs {
                ce: dt * k,
                c0: 1.0 - dt / tau,
                c1: 0.0,
                c2: (dt / tau) * EPS0 * delta_eps,
                q: -k,
                s: -1.0,
            }
        }
    }
}
