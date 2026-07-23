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
use crate::dispersive::{CpuDispersiveState, DispersiveMap};
use crate::drive::{Drive, EComponent, HComponent};
use crate::fields::Fields;
use crate::materials::{Boundary, Materials};
use crate::spec::{FdtdSpec, GradedSpacings, SpacingArrays, idx3, len3};

/// Lossy CA/CB coefficients for the Yee E-update (Taflove §3.7), identical
/// to `yee_fdtd::update::ca_cb`.
#[inline]
fn ca_cb(eps_r: f64, sigma: f64, dt: f64) -> (f64, f64) {
    let denom = 2.0 * EPS0 * eps_r + sigma * dt;
    let ca = (2.0 * EPS0 * eps_r - sigma * dt) / denom;
    let cb = 2.0 * dt / denom;
    (ca, cb)
}

/// Fixed-length row view into a flat field buffer (E.4 row-sliced kernels).
#[inline]
fn row(buf: &[f64], lo: usize, len: usize) -> &[f64] {
    &buf[lo..lo + len]
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
    dispersive: Option<(DispersiveMap, CpuDispersiveState)>,
    drive: Drive,
    /// `e_z_prev` per resistive port (`LumpedRlcPort::e_z_prev` equivalent).
    port_state: Vec<f64>,
    /// Cached modal terminal voltage `V_T^n` per aperture port
    /// (`LumpedRlcPort::e_z_prev` in the aperture sense).
    aperture_state: Vec<f64>,
    /// Recorded probe series, one inner `Vec` per [`Drive::probes`] entry.
    probe_series: Vec<Vec<f64>>,
    /// Recorded H-probe series (FS.4.2a), one inner `Vec` per
    /// [`Drive::h_probes`] entry.
    h_probe_series: Vec<Vec<f64>>,
    /// Per-step `(v_src, v_terminal, i_branch)` per aperture port (empty
    /// inner `Vec` when `AperturePort::record` is false) — FS.2a,
    /// ADR-0207.
    aperture_records: Vec<Vec<(f64, f64, f64)>>,
    /// Per-axis primal/dual spacings the kernels divide by (FS.0b.0). The
    /// uniform constructors fill these with the scalar `spec.dx/dy/dz`, so
    /// the uniform path is bit-exact by construction (gate `compute-018`).
    spacings: SpacingArrays,
    /// True once [`CpuFdtd::set_spacings`] attached a graded grid.
    graded: bool,
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
        Self::with_drive(spec, fields, materials, boundary, Drive::default())
    }

    /// Build a driven stepper (E.2): soft sources injected between the H and
    /// E half-steps, resistive ports applied after the E boundary phase, and
    /// probes recorded once per step — the exact `WalkingSkeletonSolver` +
    /// `LumpedRlcPort` orchestration.
    ///
    /// # Panics
    ///
    /// Panics on any shape mismatch or out-of-bounds drive cell.
    pub fn with_drive(
        spec: FdtdSpec,
        fields: Fields,
        materials: Materials,
        boundary: Boundary,
        drive: Drive,
    ) -> Self {
        assert_eq!(fields.ex.len(), len3(spec.ex_dims()), "ex length mismatch");
        assert_eq!(fields.ey.len(), len3(spec.ey_dims()), "ey length mismatch");
        assert_eq!(fields.ez.len(), len3(spec.ez_dims()), "ez length mismatch");
        assert_eq!(fields.hx.len(), len3(spec.hx_dims()), "hx length mismatch");
        assert_eq!(fields.hy.len(), len3(spec.hy_dims()), "hy length mismatch");
        assert_eq!(fields.hz.len(), len3(spec.hz_dims()), "hz length mismatch");
        materials.validate(&spec);
        drive.validate(&spec);
        let (cpml, pec_box) = match boundary {
            Boundary::None => (None, false),
            Boundary::PecBox => (None, true),
            Boundary::Cpml(config) => (Some(CpuCpmlState::new(&spec, config)), false),
        };
        let port_state = vec![0.0; drive.ports.len()];
        let aperture_state = vec![0.0; drive.aperture_ports.len()];
        let probe_series = vec![Vec::new(); drive.probes.len()];
        let h_probe_series = vec![Vec::new(); drive.h_probes.len()];
        let aperture_records = vec![Vec::new(); drive.aperture_ports.len()];
        let spacings = SpacingArrays::uniform(&spec);
        Self {
            spec,
            fields,
            materials,
            cpml,
            pec_box,
            step: 0,
            dispersive: None,
            drive,
            port_state,
            aperture_state,
            probe_series,
            h_probe_series,
            aperture_records,
            spacings,
            graded: false,
        }
    }

    /// Attach per-axis nonuniform primal spacings (FS.0b.0, ADR-0208): the
    /// H updates divide by the primal cell width at the H sample, the E
    /// updates by the dual spacing (average of the two adjacent primal
    /// cells; the single adjacent primal at the domain edges). Constant
    /// arrays equal to the scalar `spec.dx/dy/dz` are bit-exact against the
    /// uniform path (gate `compute-018`). Call before stepping.
    ///
    /// # Panics
    ///
    /// Panics on invalid spacings (wrong lengths, non-positive widths),
    /// grading inside an absorbing CPML layer, or when a dispersive map is
    /// attached (the ADE E-step is uniform-only in FS.0b.0). Untrusted
    /// callers should pre-flight with [`GradedSpacings::validate`] /
    /// [`GradedSpacings::validate_cpml_layers`].
    pub fn set_spacings(&mut self, graded: &GradedSpacings) {
        if let Err(e) = graded.validate(&self.spec) {
            panic!("set_spacings: {e}");
        }
        if let Some(cpml) = &self.cpml
            && let Err(e) = graded.validate_cpml_layers(cpml.npml, cpml.faces)
        {
            panic!("set_spacings: {e}");
        }
        assert!(
            self.dispersive.is_none(),
            "set_spacings: dispersive ADE materials are uniform-grid only (FS.0b.0)"
        );
        assert!(
            self.spec.dt <= graded.courant_limit(),
            "set_spacings: dt = {} s exceeds the graded Courant limit {} s \
             (use the minimum spacing per axis)",
            self.spec.dt,
            graded.courant_limit()
        );
        self.spacings = SpacingArrays::graded(graded);
        self.graded = true;
    }

    /// Attach an ADE dispersive-material map (E.5c). Replaces the standard
    /// E half-step with the fused ADE update (`yee_fdtd::dispersive`
    /// semantics: the non-dispersive per-cell ε/σ arms do not apply — the
    /// vacuum ADE arm assumes ε_r = 1 like the reference).
    ///
    /// # Panics
    ///
    /// Panics if per-cell ε/σ maps are also attached, or on a length
    /// mismatch.
    pub fn set_dispersive(&mut self, map: DispersiveMap) {
        assert!(
            self.materials.eps_r_cells.is_none() && self.materials.sigma_cells.is_none(),
            "dispersive map is exclusive with per-cell eps/sigma (reference semantics)"
        );
        assert!(
            !self.graded,
            "dispersive ADE materials are uniform-grid only (FS.0b.0)"
        );
        map.validate(&self.spec);
        let state = CpuDispersiveState::new(&self.spec);
        self.dispersive = Some((map, state));
    }

    /// The problem description this stepper was built from.
    pub fn spec(&self) -> &FdtdSpec {
        &self.spec
    }

    /// Simulation time at the start of the next step (seconds).
    pub fn current_time(&self) -> f64 {
        self.step as f64 * self.spec.dt
    }

    /// Advance the state by `n` full leapfrog steps, applying the drive (if
    /// any): H → boundary-H → soft sources → E → boundary-E (incl. masks) →
    /// resistive ports → probe recording → clock.
    pub fn step_n(&mut self, n: usize) {
        use yee_core::units::EPS0;
        let dt = self.spec.dt;
        for _ in 0..n {
            let n_step = self.step as usize;
            self.update_h();
            self.boundary_h();
            if !self.drive.soft_sources.is_empty() {
                for src in &self.drive.soft_sources {
                    let flat = src.component.flat(&self.spec, src.cell);
                    let field = match src.component {
                        EComponent::Ex => &mut self.fields.ex,
                        EComponent::Ey => &mut self.fields.ey,
                        EComponent::Ez => &mut self.fields.ez,
                    };
                    field[flat] += src.waveform.value(n_step, dt);
                }
            }
            match self.dispersive.take() {
                None => self.update_e(),
                Some((map, mut state)) => {
                    state.update_e(&self.spec, &mut self.fields, &map);
                    self.dispersive = Some((map, state));
                }
            }
            self.boundary_e();
            if !self.drive.ports.is_empty() {
                // Verbatim `LumpedRlcPort::{correct_e, update_pure_resistor}`,
                // with the local cell sizes at the port cell (dual transverse
                // spacings × primal dz): bit-equal to `s.dx*s.dy` / `s.dz` on
                // a uniform grid (compute-018).
                let s = self.spec;
                for (port, e_z_prev) in self.drive.ports.iter().zip(&mut self.port_state) {
                    let (ci, cj, ck) = port.cell;
                    let area = self.spacings.x.dual[ci] * self.spacings.y.dual[cj];
                    let dz_c = self.spacings.z.primal[ck];
                    let flat = EComponent::Ez.flat(&s, port.cell);
                    let e1_star = self.fields.ez[flat];
                    let e0 = *e_z_prev;
                    let v_src = port.waveform.value(n_step, dt);
                    let alpha = dt * dz_c / (2.0 * EPS0 * port.resistance * area);
                    let gamma = dt / (EPS0 * port.resistance * area);
                    let e1 = (e1_star - alpha * e0 + gamma * v_src) / (1.0 + alpha);
                    self.fields.ez[flat] = e1;
                    *e_z_prev = e1;
                }
            }
            if !self.drive.aperture_ports.is_empty() {
                // Verbatim `LumpedRlcPort::correct_e_aperture` (pure-R arm):
                // modal V = ∫E_z·dz averaged over the width columns, aggregate
                // branch current with β = dt·h/(2·ε₀·A), sheet-current
                // back-action referenced to the physical area A. The height
                // integral uses each cell's primal dz (bit-equal to `s.dz` on
                // a uniform grid, compute-018).
                let s = self.spec;
                let dzp = &self.spacings.z.primal;
                for ((port, v_prev), record) in self
                    .drive
                    .aperture_ports
                    .iter()
                    .zip(&mut self.aperture_state)
                    .zip(&mut self.aperture_records)
                {
                    let v_src = port.waveform.value(n_step, dt);
                    let n_col = port.n_columns as f64;
                    let mut v_sum = 0.0;
                    for &cell in &port.cells {
                        v_sum += self.fields.ez[EComponent::Ez.flat(&s, cell)] * dzp[cell.2];
                    }
                    let v_term_star = v_sum / n_col;
                    let v_term_mid = 0.5 * (v_term_star + *v_prev);
                    let beta = dt * port.height / (2.0 * EPS0 * port.area);
                    let i_branch = if port.resistance.is_infinite() {
                        0.0
                    } else {
                        (v_term_mid - v_src) / (port.resistance + beta)
                    };
                    if port.record {
                        record.push((v_src, v_term_mid, i_branch));
                    }
                    let back = (dt / (EPS0 * port.area)) * i_branch;
                    for &cell in &port.cells {
                        self.fields.ez[EComponent::Ez.flat(&s, cell)] -= back;
                    }
                    let mut v_sum_post = 0.0;
                    for &cell in &port.cells {
                        v_sum_post += self.fields.ez[EComponent::Ez.flat(&s, cell)] * dzp[cell.2];
                    }
                    *v_prev = v_sum_post / n_col;
                }
            }
            self.step += 1;
            if !self.drive.probes.is_empty() {
                for (probe, series) in self.drive.probes.iter().zip(&mut self.probe_series) {
                    let flat = probe.component.flat(&self.spec, probe.cell);
                    let field = match probe.component {
                        EComponent::Ex => &self.fields.ex,
                        EComponent::Ey => &self.fields.ey,
                        EComponent::Ez => &self.fields.ez,
                    };
                    series.push(field[flat]);
                }
            }
            if !self.drive.h_probes.is_empty() {
                // Same step iteration as the E-probe recording above, so the
                // two series share an index — but the H state read here was
                // written by `update_h` at the TOP of this iteration
                // (t = (n+½)·Δt), a half-step behind the E sample just taken
                // above (t = (n+1)·Δt). See `HProbe`'s doc comment.
                for (probe, series) in self.drive.h_probes.iter().zip(&mut self.h_probe_series) {
                    let flat = probe.component.flat(&self.spec, probe.cell);
                    let field = match probe.component {
                        HComponent::Hx => &self.fields.hx,
                        HComponent::Hy => &self.fields.hy,
                        HComponent::Hz => &self.fields.hz,
                    };
                    series.push(field[flat]);
                }
            }
        }
    }

    /// Recorded probe series (one inner slice per [`Drive::probes`] entry,
    /// one sample per completed step).
    pub fn probe_series(&self) -> &[Vec<f64>] {
        &self.probe_series
    }

    /// Recorded H-probe series (FS.4.2a; one inner slice per
    /// [`Drive::h_probes`] entry, one sample per completed step — see
    /// `HProbe`'s doc comment for the half-step timing relative to
    /// [`CpuFdtd::probe_series`]).
    pub fn h_probe_series(&self) -> &[Vec<f64>] {
        &self.h_probe_series
    }

    /// Per-step `(v_src, v_terminal, i_branch)` records, one inner slice
    /// per [`Drive::aperture_ports`] entry (empty unless the port set
    /// [`AperturePort::record`]). Sign convention: positive `i_branch`
    /// flows FROM the aperture INTO the branch (`i = (v_term − v_src) /
    /// (R + β)`), so a passive load reads `v_term·i ≥ 0` and the EMF's
    /// supplied power is `−v_src·i`.
    pub fn aperture_records(&self) -> &[Vec<(f64, f64, f64)>] {
        &self.aperture_records
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
                &self.spacings,
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
                &self.spacings,
            );
        } else if self.pec_box {
            apply_pec_box(&self.spec, &mut self.fields);
        }
        self.apply_pec_mask();
    }

    /// Apply the interior conductor model to masked E edges, mirroring
    /// `YeeGrid::apply_pec_mask`: plain PEC clamps `E_tan = 0`; with a
    /// resistive sheet attached (R.0b, ADR-0202) the masked `E_x`/`E_y`
    /// edges instead carry the Leontovich sheet relation
    /// `E_tan = R_s·K`, `K = ẑ × (H_above − H_below)` — planar (z-normal)
    /// trace loss at the design frequency. `E_z` masks (vias) always stay
    /// PEC, and `R_s = 0`/`None` reproduces the PEC clamp bit-exactly.
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
        let Some(r_s) = self.materials.sheet_r_ohm else {
            return;
        };
        if r_s == 0.0 {
            return;
        }
        let s = self.spec;
        // E_x(i+1/2, j, k): K_x = H_y(k-1/2) - H_y(k+1/2).
        if let Some(mask) = &self.materials.pec_mask_ex {
            let exd = s.ex_dims();
            let hyd = s.hy_dims();
            for i in 0..exd.0 {
                for j in 0..exd.1 {
                    for k in 1..exd.2 - 1 {
                        if mask[idx3(exd, i, j, k)] && j < hyd.1 {
                            self.fields.ex[idx3(exd, i, j, k)] = r_s
                                * (self.fields.hy[idx3(hyd, i, j, k - 1)]
                                    - self.fields.hy[idx3(hyd, i, j, k)]);
                        }
                    }
                }
            }
        }
        // E_y(i, j+1/2, k): K_y = H_x(k+1/2) - H_x(k-1/2).
        if let Some(mask) = &self.materials.pec_mask_ey {
            let eyd = s.ey_dims();
            let hxd = s.hx_dims();
            for i in 0..eyd.0 {
                for j in 0..eyd.1 {
                    for k in 1..eyd.2 - 1 {
                        if mask[idx3(eyd, i, j, k)] && i < hxd.0 {
                            self.fields.ey[idx3(eyd, i, j, k)] = r_s
                                * (self.fields.hx[idx3(hxd, i, j, k)]
                                    - self.fields.hx[idx3(hxd, i, j, k - 1)]);
                        }
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
        // H samples sit at cell centres: curl-E differences span one PRIMAL
        // cell (FS.0b.0). Uniform fill == the old scalar divisors bit-exact.
        let (dxp, dyp, dzp) = (
            self.spacings.x.primal.as_slice(),
            self.spacings.y.primal.as_slice(),
            self.spacings.z.primal.as_slice(),
        );
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

        // Row-sliced inner loops (E.4): fixed-length row slices let LLVM
        // elide bounds checks and vectorize; the material branch is hoisted
        // to row level. Per-cell arithmetic is unchanged — bit-exactness is
        // enforced by compute-001/003/007.

        // ---- H_x: shape [nx+1, ny, nz], full extent ----
        hx.par_chunks_mut(hxd.1 * hxd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                for j in 0..s.ny {
                    let hx_row = &mut slab[j * hxd.2..(j + 1) * hxd.2];
                    let ey_row = row(ey, idx3(eyd, i, j, 0), s.nz + 1);
                    let ez_row0 = row(ez, idx3(ezd, i, j, 0), s.nz);
                    let ez_row1 = row(ez, idx3(ezd, i, j + 1, 0), s.nz);
                    let dy_j = dyp[j];
                    match mu_r_cells {
                        None => {
                            for k in 0..s.nz {
                                let dey_dz = (ey_row[k + 1] - ey_row[k]) / dzp[k];
                                let dez_dy = (ez_row1[k] - ez_row0[k]) / dy_j;
                                hx_row[k] += coeff_scalar * (dey_dz - dez_dy);
                            }
                        }
                        Some(m) => {
                            let m_row = row(m, idx3(celld, i, j, 0), s.nz);
                            for k in 0..s.nz {
                                let dey_dz = (ey_row[k + 1] - ey_row[k]) / dzp[k];
                                let dez_dy = (ez_row1[k] - ez_row0[k]) / dy_j;
                                let coeff = s.dt / (MU0 * m_row[k]);
                                hx_row[k] += coeff * (dey_dz - dez_dy);
                            }
                        }
                    }
                }
            });

        // ---- H_y: shape [nx, ny+1, nz], full extent ----
        hy.par_chunks_mut(hyd.1 * hyd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                let dx_i = dxp[i];
                for j in 0..=s.ny {
                    let hy_row = &mut slab[j * hyd.2..(j + 1) * hyd.2];
                    let ez_row0 = row(ez, idx3(ezd, i, j, 0), s.nz);
                    let ez_row1 = row(ez, idx3(ezd, i + 1, j, 0), s.nz);
                    let ex_row = row(ex, idx3(exd, i, j, 0), s.nz + 1);
                    match mu_r_cells {
                        None => {
                            for k in 0..s.nz {
                                let dez_dx = (ez_row1[k] - ez_row0[k]) / dx_i;
                                let dex_dz = (ex_row[k + 1] - ex_row[k]) / dzp[k];
                                hy_row[k] += coeff_scalar * (dez_dx - dex_dz);
                            }
                        }
                        Some(m) => {
                            let m_row = row(m, idx3(celld, i, j, 0), s.nz);
                            for k in 0..s.nz {
                                let dez_dx = (ez_row1[k] - ez_row0[k]) / dx_i;
                                let dex_dz = (ex_row[k + 1] - ex_row[k]) / dzp[k];
                                let coeff = s.dt / (MU0 * m_row[k]);
                                hy_row[k] += coeff * (dez_dx - dex_dz);
                            }
                        }
                    }
                }
            });

        // ---- H_z: shape [nx, ny, nz+1], full extent ----
        hz.par_chunks_mut(hzd.1 * hzd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                let dx_i = dxp[i];
                for j in 0..s.ny {
                    let hz_row = &mut slab[j * hzd.2..(j + 1) * hzd.2];
                    let ex_row0 = row(ex, idx3(exd, i, j, 0), s.nz + 1);
                    let ex_row1 = row(ex, idx3(exd, i, j + 1, 0), s.nz + 1);
                    let ey_row0 = row(ey, idx3(eyd, i, j, 0), s.nz + 1);
                    let ey_row1 = row(ey, idx3(eyd, i + 1, j, 0), s.nz + 1);
                    let dy_j = dyp[j];
                    match mu_r_cells {
                        None => {
                            for k in 0..=s.nz {
                                let dex_dy = (ex_row1[k] - ex_row0[k]) / dy_j;
                                let dey_dx = (ey_row1[k] - ey_row0[k]) / dx_i;
                                hz_row[k] += coeff_scalar * (dex_dy - dey_dx);
                            }
                        }
                        Some(m) => {
                            let m_row = row(m, idx3(celld, i, j, 0), s.nz + 1);
                            for k in 0..=s.nz {
                                let dex_dy = (ex_row1[k] - ex_row0[k]) / dy_j;
                                let dey_dx = (ey_row1[k] - ey_row0[k]) / dx_i;
                                let coeff = s.dt / (MU0 * m_row[k]);
                                hz_row[k] += coeff * (dex_dy - dey_dx);
                            }
                        }
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
        // E nodes sit between H samples at cell centres: curl-H differences
        // span the DUAL spacing (FS.0b.0). Uniform fill == scalar divisors.
        let (dxd, dyd, dzd) = (
            self.spacings.x.dual.as_slice(),
            self.spacings.y.dual.as_slice(),
            self.spacings.z.dual.as_slice(),
        );
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

        // Row-sliced inner loops with the (σ, ε) material match hoisted to
        // row level (E.4). Per-cell arithmetic and operation order are the
        // reference's, verbatim — the four arms below are the flattening of
        // the reference's nested per-cell match.
        macro_rules! e_rows {
            ($e_row:ident, $curl:ident, $cell_lo:expr, $range:expr) => {{
                match (sigma_cells, eps_r_cells) {
                    (None, None) => {
                        for k in $range {
                            $e_row[k] += coeff_scalar * $curl(k);
                        }
                    }
                    (None, Some(eps)) => {
                        let eps_row = row(eps, $cell_lo, celld.2);
                        for k in $range {
                            let coeff = s.dt / (EPS0 * eps_row[k]);
                            $e_row[k] += coeff * $curl(k);
                        }
                    }
                    (Some(sig), None) => {
                        let sig_row = row(sig, $cell_lo, celld.2);
                        for k in $range {
                            let (ca, cb) = ca_cb(s.eps_r, sig_row[k], s.dt);
                            $e_row[k] = ca * $e_row[k] + cb * $curl(k);
                        }
                    }
                    (Some(sig), Some(eps)) => {
                        let sig_row = row(sig, $cell_lo, celld.2);
                        let eps_row = row(eps, $cell_lo, celld.2);
                        for k in $range {
                            let (ca, cb) = ca_cb(eps_row[k], sig_row[k], s.dt);
                            $e_row[k] = ca * $e_row[k] + cb * $curl(k);
                        }
                    }
                }
            }};
        }

        // ---- E_x: shape [nx, ny+1, nz+1]; interior j ∈ [1, ny), k ∈ [1, nz) ----
        // Outer tangential faces are managed by the boundary phase.
        ex.par_chunks_mut(exd.1 * exd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                for j in 1..s.ny {
                    let ex_row = &mut slab[j * exd.2..(j + 1) * exd.2];
                    let hz_row0 = row(hz, idx3(hzd, i, j - 1, 0), s.nz + 1);
                    let hz_row1 = row(hz, idx3(hzd, i, j, 0), s.nz + 1);
                    let hy_row = row(hy, idx3(hyd, i, j, 0), s.nz);
                    let dy_j = dyd[j];
                    let curl = |k: usize| {
                        let dhz_dy = (hz_row1[k] - hz_row0[k]) / dy_j;
                        let dhy_dz = (hy_row[k] - hy_row[k - 1]) / dzd[k];
                        dhz_dy - dhy_dz
                    };
                    e_rows!(ex_row, curl, idx3(celld, i, j, 0), 1..s.nz);
                }
            });

        // ---- E_y: shape [nx+1, ny, nz+1]; interior i ∈ [1, nx), k ∈ [1, nz) ----
        ey.par_chunks_mut(eyd.1 * eyd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                if i == 0 || i >= s.nx {
                    return;
                }
                let dx_i = dxd[i];
                for j in 0..s.ny {
                    let ey_row = &mut slab[j * eyd.2..(j + 1) * eyd.2];
                    let hx_row = row(hx, idx3(hxd, i, j, 0), s.nz);
                    let hz_row0 = row(hz, idx3(hzd, i - 1, j, 0), s.nz + 1);
                    let hz_row1 = row(hz, idx3(hzd, i, j, 0), s.nz + 1);
                    let curl = |k: usize| {
                        let dhx_dz = (hx_row[k] - hx_row[k - 1]) / dzd[k];
                        let dhz_dx = (hz_row1[k] - hz_row0[k]) / dx_i;
                        dhx_dz - dhz_dx
                    };
                    e_rows!(ey_row, curl, idx3(celld, i, j, 0), 1..s.nz);
                }
            });

        // ---- E_z: shape [nx+1, ny+1, nz]; interior i ∈ [1, nx), j ∈ [1, ny) ----
        ez.par_chunks_mut(ezd.1 * ezd.2)
            .enumerate()
            .for_each(|(i, slab)| {
                if i == 0 || i >= s.nx {
                    return;
                }
                let dx_i = dxd[i];
                for j in 1..s.ny {
                    let ez_row = &mut slab[j * ezd.2..(j + 1) * ezd.2];
                    let hy_row0 = row(hy, idx3(hyd, i - 1, j, 0), s.nz);
                    let hy_row1 = row(hy, idx3(hyd, i, j, 0), s.nz);
                    let hx_row0 = row(hx, idx3(hxd, i, j - 1, 0), s.nz);
                    let hx_row1 = row(hx, idx3(hxd, i, j, 0), s.nz);
                    let dy_j = dyd[j];
                    let curl = |k: usize| {
                        let dhy_dx = (hy_row1[k] - hy_row0[k]) / dx_i;
                        let dhx_dy = (hx_row1[k] - hx_row0[k]) / dy_j;
                        dhy_dx - dhx_dy
                    };
                    e_rows!(ez_row, curl, idx3(celld, i, j, 0), 0..s.nz);
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
