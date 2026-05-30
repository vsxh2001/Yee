//! cap-cw-001: a CW (single-frequency steady-state) diagnostic of the aperture
//! lumped-port CAPACITOR arm (Phase 2.fdtd.6.10, ADR-0125/0126 follow-on).
//!
//! # Why this test exists
//!
//! ADR-0126 wired F2.3's band-pass onto the aperture port ([`LumpedRlcPort::aperture`])
//! and found that, while the shunt L‖C tanks now *load* the line (the O(dx²)
//! inductor collapse of ADR-0124 is gone), the band-pass still doesn't form and —
//! decisively — a *longer FDTD record makes the in-band loss WORSE* (12.5→22.1 dB
//! at 24k→60k steps). For a linear, passive element a pulse→DFT transfer function
//! must NOT degrade with window length, so ADR-0126 flagged the shunt-tank
//! **capacitor** as reading a "deepening near-short over time" — a possible
//! steady-state / DC-windup defect in the cap arm rather than merely an unsettled
//! pulse.
//!
//! # Two probes, one verdict
//!
//! The decisive question — *does the cap ARM have a steady-state defect?* — is
//! cleanest when the arm is driven by a **controlled** terminal voltage, free of
//! line standing-wave / de-embed noise. So this test runs two probes:
//!
//! 1. **ISOLATED cap arm under CW (asserted).** The aperture port's update
//!    ([`LumpedRlcPort::correct_e_aperture`]) is exercised on a single field cell
//!    whose pre-correction value is forced to a clean CW history, in the SAME
//!    closed loop the real driver uses (the port back-action modifies the cell,
//!    which feeds the next step's terminal voltage). We read the port's OWN
//!    `(V_T, I, V_C)` each step and form `Z = V_T(ω)/I(ω)` over sliding windows.
//!    This is exactly the Piket-May/Taflove cap update with NO line in the loop,
//!    so it isolates the arm's steady-state behaviour. Does `|Z|` CONVERGE to the
//!    analytic `1/(ωC)` (a fixed `−jX`) or DRIFT toward a short? Does `V_C` settle
//!    to a bounded oscillation or WIND UP? This block is asserted.
//!
//! 2. **On-guide field-mediated probe (recorded, NOT asserted).** The same port
//!    on a PEC parallel-plate guide under a soft CW source sheet. The realized
//!    `Z = V_T/I` here is corrupted by the short-line standing-wave de-embed (the
//!    exact fragility `aperture_port_001.rs` documents — `Z₀` itself swings and
//!    the cap branch current is 90° from `V_T`), so it is logged for context but
//!    not used for the verdict.
//!
//! # Verdict logic (ADR-0126)
//!
//! - **MEASUREMENT-limit:** the isolated cap presents a STABLE ≈`1/(jωC)` and
//!   `V_C` is bounded → the F2.3 failure is the PULSE drive (DFT of an unsettled
//!   transient on a noisy short line), NOT the cap arm. The fix is downstream (a
//!   CW drive + clean de-embed in F2.3); the cap arm is correct.
//! - **CAP-UPDATE BUG:** the isolated cap's `|Z|` drifts toward a short / `V_C`
//!   winds up → the cap arm has a steady-state defect to fix in `lumped.rs`.

use std::f64::consts::PI;

use yee_core::units::EPS0;
use yee_fdtd::{
    ApertureSpec, LumpedRlcPort, SourceWaveform, WalkingSkeletonSolver, YeeGrid, boundary, update,
};

/// CW drive frequency (Hz). 2 GHz matches the F2.3 band centre.
const F0: f64 = 2.0e9;
/// Cell size (m) for the FDTD timestep / aperture references.
const DX: f64 = 1.0e-3;

/// A minimal complex number so the bench has no extra deps.
#[derive(Clone, Copy, Debug)]
struct Cplx {
    re: f64,
    im: f64,
}
impl Cplx {
    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    fn div(self, o: Cplx) -> Cplx {
        let d = o.re * o.re + o.im * o.im;
        Cplx::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
}

/// One complex DFT bin of `series[a..b]`, `Σ v·e^{-jωt}` with absolute step
/// indices (so different windows share a common phase reference).
fn dft_bin(series: &[f64], a: usize, b: usize, omega: f64, dt: f64) -> Cplx {
    let mut re = 0.0_f64;
    let mut im = 0.0_f64;
    for (n, &v) in series[a..b].iter().enumerate() {
        let ph = omega * ((a + n) as f64) * dt;
        re += v * ph.cos();
        im -= v * ph.sin();
    }
    Cplx::new(re, im)
}

/// Sliding-window `Z = V_T/I` at `omega`: `win`-wide windows striding by `win`
/// from `start`. Returns `(window_end_step, Z, |Z|)` per stride.
fn sliding_z(
    v_t: &[f64],
    i_branch: &[f64],
    dt: f64,
    omega: f64,
    start: usize,
    win: usize,
) -> Vec<(usize, Cplx, f64)> {
    let mut out = Vec::new();
    let n = v_t.len();
    let mut a = start;
    while a + win <= n {
        let b = a + win;
        let vf = dft_bin(v_t, a, b, omega, dt);
        let iff = dft_bin(i_branch, a, b, omega, dt);
        let z = vf.div(iff);
        out.push((b, z, z.abs()));
        a += win;
    }
    out
}

/// Analytic capacitor reactance `−1/(ωC)` (Ω).
fn x_cap(c: f64, omega: f64) -> f64 {
    -1.0 / (omega * c)
}

// =====================================================================
// Probe 1 — ISOLATED cap arm under CW (the decisive, asserted probe).
// =====================================================================

/// A degenerate single-cell aperture (`A = dx²`, `h = dz = dx`) so the aperture
/// references collapse onto a clean per-cell capacitor — the cap arm in isolation.
fn one_cell_spec() -> ApertureSpec {
    ApertureSpec {
        cells: vec![(0, 0, 0)],
        n_columns: 1,
        area: DX * DX,
        height: DX,
    }
}

/// Recorded traces from one isolated-arm CW run.
struct ArmRun {
    /// Port modal terminal voltage `V_T[n]` (V).
    v_t: Vec<f64>,
    /// Aggregate aperture branch current `I[n]` (A).
    i_branch: Vec<f64>,
    /// Capacitor-voltage state `V_C[n]` (V). Empty if no capacitive port.
    v_c: Vec<f64>,
    /// The run's FDTD time step (s) — used to phase the sliding-window DFT.
    dt: f64,
}

/// Drive the isolated cap arm with a Hann-ramped CW terminal voltage in the SAME
/// closed loop the real driver uses: each step the field cell's pre-correction
/// value is `E_prev + soft_CW`, then `correct_e_aperture` solves the branch
/// current and applies the back-action (which `E_prev` carries to the next step).
///
/// This is the Piket-May/Taflove aperture cap update with no transmission line in
/// the feedback path — so the realized `(V_T, I, V_C)` reflect the ARM alone.
fn run_isolated_arm(r: f64, l: f64, c: f64, f0: f64, n_steps: usize, ramp_steps: usize) -> ArmRun {
    // A tiny 1×1×1-usable vacuum grid; only cell (0,0,0)'s E_z is touched.
    let grid = YeeGrid::vacuum(3, 3, 3, DX);
    let dt = grid.dt;
    let mut solver = WalkingSkeletonSolver::new(grid);

    let mut port = LumpedRlcPort::aperture(one_cell_spec(), r, l, c, SourceWaveform::None);
    let src = SourceWaveform::HannSine {
        v0: 1.0,
        frequency: f0,
        ramp_steps,
    };

    let mut run = ArmRun {
        v_t: Vec::with_capacity(n_steps),
        i_branch: Vec::with_capacity(n_steps),
        v_c: Vec::with_capacity(n_steps),
        dt,
    };
    // The cell's E_z carries between steps (the closed loop). Each step we ADD a
    // soft CW increment (like a soft source) so the drive is a clean single tone
    // and the port's back-action genuinely feeds back.
    for n in 0..n_steps {
        let (grid, _) = solver.grid_and_cpml_mut();
        grid.ez[(0, 0, 0)] += src.value(n, dt);
        port.correct_e_aperture(grid, n, dt);
        run.v_t.push(port.last_terminal_voltage());
        run.i_branch.push(port.inductor_current());
        run.v_c.push(port.capacitor_voltage());
        solver.advance_clock();
    }
    run
}

/// Peak `|V_C|` over `[a, b)` (steady-state envelope tracker).
fn vc_envelope(v_c: &[f64], a: usize, b: usize) -> f64 {
    v_c[a..b].iter().fold(0.0_f64, |m, &v| m.max(v.abs()))
}

// =====================================================================
// Probe 2 — on-guide field-mediated de-embed (recorded, NOT asserted).
// =====================================================================

const NX: usize = 200;
const NT: usize = 6;
const SRC_I: usize = 20;
const PORT_I: usize = 100;

fn guide_aperture_cells() -> Vec<(usize, usize, usize)> {
    let mut v = Vec::new();
    for j in 1..NT {
        for k in 0..NT {
            v.push((PORT_I, j, k));
        }
    }
    v
}

fn guide_spec() -> ApertureSpec {
    let w = (NT - 1) as f64 * DX;
    let height = NT as f64 * DX;
    ApertureSpec {
        cells: guide_aperture_cells(),
        n_columns: NT - 1,
        area: w * height,
        height,
    }
}

/// One PEC FDTD step with a Hann-ramped full-width CW `E_z` source sheet, then
/// the aperture-load correction(s).
fn step_guide(
    solver: &mut WalkingSkeletonSolver,
    ports: &mut [LumpedRlcPort],
    n_step: usize,
    dt: f64,
    src: &SourceWaveform,
) {
    {
        let (grid, _) = solver.grid_and_cpml_mut();
        update::update_h(grid);
        #[allow(deprecated)]
        boundary::apply_pec(grid);
    }
    {
        let (grid, _) = solver.grid_and_cpml_mut();
        let s = src.value(n_step, dt);
        for j in 1..NT {
            for k in 0..NT {
                grid.ez[(SRC_I, j, k)] += s;
            }
        }
    }
    {
        let (grid, _) = solver.grid_and_cpml_mut();
        update::update_e(grid);
        #[allow(deprecated)]
        boundary::apply_pec(grid);
    }
    let (grid, _) = solver.grid_and_cpml_mut();
    for p in ports.iter_mut() {
        p.correct_e_aperture(grid, n_step, dt);
    }
    solver.advance_clock();
}

/// On-guide run: record the (shared) modal `V_T`, the aggregate branch current,
/// and (if a capacitive port exists) its `V_C`.
fn run_guide(loads: &[(f64, f64, f64)], n_steps: usize, ramp_steps: usize) -> ArmRun {
    let grid = YeeGrid::vacuum(NX, NT, NT, DX);
    let dt = grid.dt;
    let mut solver = WalkingSkeletonSolver::new(grid);
    let spec = guide_spec();
    let mut ports: Vec<LumpedRlcPort> = loads
        .iter()
        .map(|&(r, l, c)| LumpedRlcPort::aperture(spec.clone(), r, l, c, SourceWaveform::None))
        .collect();
    let cap_idx = loads.iter().position(|&(_, _, c)| c.is_finite());
    let src = SourceWaveform::HannSine {
        v0: 1.0,
        frequency: F0,
        ramp_steps,
    };
    let mut run = ArmRun {
        v_t: Vec::with_capacity(n_steps),
        i_branch: Vec::with_capacity(n_steps),
        v_c: Vec::with_capacity(n_steps),
        dt,
    };
    for n in 0..n_steps {
        step_guide(&mut solver, &mut ports, n, dt, &src);
        run.v_t.push(ports[0].last_terminal_voltage());
        run.i_branch
            .push(ports.iter().map(|p| p.inductor_current()).sum());
        if let Some(ci) = cap_idx {
            run.v_c.push(ports[ci].capacitor_voltage());
        }
    }
    run
}

#[test]
#[ignore = "slow: CW steady-state aperture-cap diagnostic (Phase 2.fdtd.6.10)"]
fn cap_cw_001() {
    let dt = YeeGrid::vacuum(3, 3, 3, DX).dt;
    let omega = 2.0 * PI * F0;
    let steps_per_cycle = (1.0 / (F0 * dt)).round() as usize;

    // |X_C| ≈ 100 Ω at f0 (same scale as aperture_port_001); the L‖C tank then
    // resonates AT f0 (ωL = 1/(ωC) ⇒ ω² = 1/LC).
    let x_ref = 100.0;
    let c_react = 1.0 / (omega * x_ref);
    let l_react = x_ref / omega;

    let ramp_steps = 30 * steps_per_cycle;
    let n_steps = 220 * steps_per_cycle;
    let settle = ramp_steps + 60 * steps_per_cycle; // probe well after the ramp
    let win = steps_per_cycle * 8; // 8-cycle DFT window → clean single-tone bin

    // β for each aperture (isolated single-cell vs the wide guide aperture).
    let beta_cell = dt * one_cell_spec().height / (2.0 * EPS0 * one_cell_spec().area);
    let gspec = guide_spec();
    let beta_guide = dt * gspec.height / (2.0 * EPS0 * gspec.area);

    eprintln!(
        "Phase 2.fdtd.6.10 — CW capacitor steady-state diagnostic (ADR-0126 follow-on)
  dt={dt:.3e} s, f0={f0} GHz, {spc} steps/cycle, n_steps={ns} (~{cyc} cyc)
  Hann ramp {rs} steps (~30 cyc), settle@{st}, win={win} ({wc} cyc)
  C={c:.3e} F (1/(ωC)={xc:.1} Ω @ f0)   L={l:.3e} H (ωL={xl:.1} Ω @ f0)
  parallel-resonant f = 1/(2π√LC) = {fr:.3} GHz (== f0)
  β_isolated(single-cell)={bc:.2} Ω   β_guide(wide aperture)={bg:.2} Ω\n",
        dt = dt,
        f0 = F0 / 1e9,
        spc = steps_per_cycle,
        ns = n_steps,
        cyc = n_steps / steps_per_cycle,
        rs = ramp_steps,
        st = settle,
        win = win,
        wc = win / steps_per_cycle,
        c = c_react,
        xc = x_ref,
        l = l_react,
        xl = x_ref,
        fr = 1.0 / (2.0 * PI * (l_react * c_react).sqrt()) / 1e9,
        bc = beta_cell,
        bg = beta_guide,
    );

    // ================================================================
    // PROBE 1 (DECISIVE): isolated cap arm under CW.
    // ================================================================
    eprintln!(
        "==== PROBE 1 (DECISIVE): ISOLATED cap arm under CW @ {f0} GHz ====",
        f0 = F0 / 1e9
    );
    // Tiny series R (constructor requires R>0; R=0 rejected). Matches the
    // reactive_deembed_001 / aperture_port_001 convention.
    let cap_run = run_isolated_arm(1.0e-6, 0.0, c_react, F0, n_steps, ramp_steps);
    assert!(
        cap_run
            .v_t
            .iter()
            .chain(cap_run.i_branch.iter())
            .chain(cap_run.v_c.iter())
            .all(|x| x.is_finite()),
        "ISOLATED CAP: a CW trace went non-finite (instability)"
    );

    let x_an = x_cap(c_react, omega);
    let z_an = Cplx::new(beta_cell, x_an); // realized Z = β + 1/(jωC)
    let cap_z = sliding_z(
        &cap_run.v_t,
        &cap_run.i_branch,
        cap_run.dt,
        omega,
        settle,
        win,
    );
    eprintln!(
        "  sliding Z=V_T/I (analytic realized Z = {br:.2} − j{xc:.2} Ω, |Z|={za:.2}):",
        br = beta_cell,
        xc = -x_an,
        za = z_an.abs()
    );
    let mut cap_abs = Vec::new();
    let mut cap_im = Vec::new();
    let mut cap_re = Vec::new();
    for &(end, z, za) in &cap_z {
        eprintln!(
            "    end={end:6} ({cyc:3} cyc) | Z={zr:9.3}{zi:+9.3}j  |Z|={za:9.3}",
            cyc = end / steps_per_cycle,
            zr = z.re,
            zi = z.im,
            za = za,
        );
        cap_abs.push(za);
        cap_im.push(z.im);
        cap_re.push(z.re);
    }
    // V_C envelope window-over-window.
    eprintln!("  V_C envelope (peak |V_C| per window):");
    let mut vc_env = Vec::new();
    {
        let mut a = settle;
        while a + win <= cap_run.v_c.len() {
            let env = vc_envelope(&cap_run.v_c, a, a + win);
            eprintln!(
                "    end={end:6} ({cyc:3} cyc) | peak|V_C|={env:.5e}",
                end = a + win,
                cyc = (a + win) / steps_per_cycle,
                env = env,
            );
            vc_env.push(env);
            a += win;
        }
    }

    // ---- isolated cap metrics ----
    let z_first = *cap_abs.first().expect("≥1 settled window");
    let z_last = *cap_abs.last().expect("≥1 settled window");
    let z_drift = (z_last - z_first).abs() / z_first.max(1.0);
    let z_mean = cap_abs.iter().sum::<f64>() / cap_abs.len() as f64;
    let z_acc = (z_mean - z_an.abs()).abs() / z_an.abs();
    let im_all_cap = cap_im.iter().all(|&x| x < 0.0);
    let vc_first = *vc_env.first().expect("≥1 V_C window");
    let vc_last = *vc_env.last().expect("≥1 V_C window");
    let vc_growth = (vc_last - vc_first) / vc_first.max(1e-30);

    eprintln!(
        "\n  ISOLATED CAP METRICS:
    |Z| first={zf:.3} Ω  last={zl:.3} Ω  (analytic |Z|={zan:.3} Ω)
    |Z| drift = {zd:.4}   |Z| mean={zm:.3} Ω  accuracy={za:.4}
    reactance all −jX (capacitive)? {sgn}
    V_C envelope first={vf:.4e} last={vl:.4e} growth={vg:+.4}",
        zf = z_first,
        zl = z_last,
        zan = z_an.abs(),
        zd = z_drift,
        zm = z_mean,
        za = z_acc,
        sgn = im_all_cap,
        vf = vc_first,
        vl = vc_last,
        vg = vc_growth,
    );

    // L‖C tank (isolated, two ports on the same single-cell aperture).
    eprintln!(
        "\n==== PROBE 1b: ISOLATED shunt L‖C tank under CW @ {f0} GHz ====",
        f0 = F0 / 1e9
    );
    // Two ports threading the SAME isolated aperture cell = parallel L‖C.
    let tank_run = {
        let grid = YeeGrid::vacuum(3, 3, 3, DX);
        let dt = grid.dt;
        let mut solver = WalkingSkeletonSolver::new(grid);
        let spec = one_cell_spec();
        let mut ports = [
            LumpedRlcPort::aperture(
                spec.clone(),
                1.0e-6,
                l_react,
                f64::INFINITY,
                SourceWaveform::None,
            ),
            LumpedRlcPort::aperture(spec, 1.0e-6, 0.0, c_react, SourceWaveform::None),
        ];
        let src = SourceWaveform::HannSine {
            v0: 1.0,
            frequency: F0,
            ramp_steps,
        };
        let mut run = ArmRun {
            v_t: Vec::with_capacity(n_steps),
            i_branch: Vec::with_capacity(n_steps),
            v_c: Vec::new(),
            dt,
        };
        for n in 0..n_steps {
            let (grid, _) = solver.grid_and_cpml_mut();
            grid.ez[(0, 0, 0)] += src.value(n, dt);
            for p in ports.iter_mut() {
                p.correct_e_aperture(grid, n, dt);
            }
            run.v_t.push(ports[0].last_terminal_voltage());
            run.i_branch
                .push(ports.iter().map(|p| p.inductor_current()).sum());
            solver.advance_clock();
        }
        run
    };
    assert!(
        tank_run
            .v_t
            .iter()
            .chain(tank_run.i_branch.iter())
            .all(|x| x.is_finite()),
        "ISOLATED L‖C TANK: a CW trace went non-finite (instability)"
    );
    let tank_z = sliding_z(
        &tank_run.v_t,
        &tank_run.i_branch,
        tank_run.dt,
        omega,
        settle,
        win,
    );
    eprintln!(
        "  sliding Z=V_T/I_aggregate (parallel resonance ⇒ |Z| ≫ single-arm |X|≈{x_ref:.0} Ω):"
    );
    let mut tank_abs = Vec::new();
    for &(end, z, za) in &tank_z {
        eprintln!(
            "    end={end:6} ({cyc:3} cyc) | Z={zr:11.2}{zi:+11.2}j  |Z|={za:11.2}",
            cyc = end / steps_per_cycle,
            zr = z.re,
            zi = z.im,
            za = za,
        );
        tank_abs.push(za);
    }
    let tank_mean = tank_abs.iter().sum::<f64>() / tank_abs.len() as f64;
    eprintln!(
        "\n  ISOLATED TANK METRIC: mean |Z|={tm:.2} Ω (single-arm |X|≈{x_ref:.0} Ω)",
        tm = tank_mean
    );

    // ================================================================
    // PROBE 2 (RECORDED, not asserted): on-guide field-mediated de-embed.
    // ================================================================
    eprintln!("\n==== PROBE 2 (RECORDED): on-guide cap de-embed (KNOWN-FRAGILE short-line) ====");
    let guide_cap = run_guide(&[(1.0e-6, 0.0, c_react)], n_steps, ramp_steps);
    let gz = sliding_z(
        &guide_cap.v_t,
        &guide_cap.i_branch,
        guide_cap.dt,
        omega,
        settle,
        win,
    );
    let g_abs: Vec<f64> = gz.iter().map(|&(_, _, za)| za).collect();
    let g_mean = if g_abs.is_empty() {
        f64::NAN
    } else {
        g_abs.iter().sum::<f64>() / g_abs.len() as f64
    };
    let g_min = g_abs.iter().cloned().fold(f64::INFINITY, f64::min);
    let g_max = g_abs.iter().cloned().fold(0.0_f64, f64::max);
    eprintln!(
        "  on-guide |Z| over settled windows: mean={gm:.1} Ω  min={gn:.1} Ω  max={gx:.1} Ω
  (the short-line standing-wave de-embed scatters |Z| wildly — recorded only, NOT a cap verdict)",
        gm = g_mean,
        gn = g_min,
        gx = g_max,
    );

    // ================================================================
    // VERDICT
    // ================================================================
    let cap_stable = z_drift < 0.10 && z_acc < 0.10 && im_all_cap && vc_growth.abs() < 0.10;
    let tank_resonates = tank_mean > x_ref;

    eprintln!("\n======= VERDICT (cap-cw-001, ADR-0126 follow-on) =======");
    if cap_stable {
        eprintln!(
            "  CAP (isolated, decisive): STABLE ≈1/(jωC) under CW —
       |Z| no drift (drift {zd:.3} < 0.10), accuracy {za:.3} < 0.10, −jX, V_C bounded (growth {vg:+.3}).
  ==> MEASUREMENT-LIMIT: the cap ARM IS FINE. F2.3's 'longer→worse' is the PULSE
      drive (DFT of an unsettled transient on a noisy short line), NOT a cap bug.
      Fix is downstream (a CW drive + clean de-embed in F2.3); do NOT change the cap arm.",
            zd = z_drift,
            za = z_acc,
            vg = vc_growth,
        );
    } else {
        eprintln!(
            "  CAP (isolated, decisive): does NOT present a stable 1/(jωC) —
       drift {zd:.3}, accuracy {za:.3}, −jX {sgn}, V_C growth {vg:+.3}.
  ==> CAP-UPDATE BUG: a steady-state defect in the cap arm (drift-to-short / windup).",
            zd = z_drift,
            za = z_acc,
            sgn = im_all_cap,
            vg = vc_growth,
        );
    }
    eprintln!(
        "  L‖C TANK (isolated): mean |Z|={tm:.1} Ω vs single-arm |X|={xr:.0} Ω ⇒ resonance {res}",
        tm = tank_mean,
        xr = x_ref,
        res = if tank_resonates {
            "PRESENT (|Z| > |X|)"
        } else {
            "ABSENT (|Z| ≤ |X|)"
        }
    );
    eprintln!("========================================================");

    // ---- ASSERTIONS (isolated arm only — the honest, line-noise-free probe) ----
    assert!(
        z_first.is_finite() && z_last.is_finite() && z_mean.is_finite() && tank_mean.is_finite(),
        "CW metric non-finite (NaN/instability)"
    );
    // The cap reactance SIGN must be capacitive (−jX) in every settled window.
    assert!(
        im_all_cap,
        "ISOLATED CAP realized reactance is not −jX in every settled window — \
         the cap arm is dead or sign-inverted"
    );
    // DECISIVE: |Z| must NOT drift over the CW record (a steady-state cap is a
    // FIXED reactance; the F2.3 'deepening short' would show as |Z|→0 here).
    assert!(
        z_drift < 0.10,
        "ISOLATED CAP |Z| DRIFTS across the CW record (drift {z_drift:.3} ≥ 0.10): \
         first={z_first:.3} Ω, last={z_last:.3} Ω — a deepening short / non-steady cap \
         (the CAP-UPDATE BUG ADR-0126 suspected)."
    );
    // DECISIVE: V_C must be a bounded steady-state oscillation (no windup).
    assert!(
        vc_growth.abs() < 0.10,
        "ISOLATED CAP V_C WINDS UP over the CW record (envelope growth {vc_growth:+.3}): \
         first={vc_first:.4e}, last={vc_last:.4e} — the cap state is not a bounded \
         steady-state oscillation (CAP-UPDATE windup bug)."
    );
    // The realized |Z| ≈ analytic (β + 1/(jωC)) under CW — correct steady-state
    // reactance. Tight, because the isolated arm has no line de-embed noise.
    assert!(
        z_acc < 0.10,
        "ISOLATED CAP realized |Z|={z_mean:.3} Ω != analytic |Z|={zan:.3} Ω under CW \
         (accuracy {z_acc:.3} ≥ 0.10) — the steady-state reactance is off.",
        zan = z_an.abs()
    );
    // The L‖C tank must RESONATE under CW: |Z| at f0 EXCEEDS the single-arm
    // reactance (parallel-resonance signature). A frozen-short cap would pin |Z|
    // low and make this impossible.
    assert!(
        tank_resonates,
        "ISOLATED L‖C TANK does NOT resonate under CW: mean |Z|={tank_mean:.2} Ω ≤ single-arm \
         |X|={x_ref:.0} Ω. A parallel tank at f0=1/(2π√LC) must present |Z| > |X|."
    );
}
