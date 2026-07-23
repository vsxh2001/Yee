//! Gate `engine-stripline-z0-001` (FS.4.2a, ADR-0225): a symmetric
//! stripline's characteristic impedance, measured from a time-gated V/I
//! ratio of the FDTD forward wave, against the exact conformal-mapping
//! closed form (Cohn 1954).
//!
//! Fixture: identical idiom to `engine-stripline-eeff-001`
//! (`stripline_eeff.rs`, FS.4.0) — `voxelize_stackup` symmetric stripline,
//! resistive-port drive, hard-PEC box, time-gated single-bin DFT. This
//! gate goes through `yee_compute::{CpuFdtd, Drive}` directly (not
//! `yee_engine::JobSpec::submit`): `JobSpec`/`ProbeSpec` only carry E
//! probes today, and every `JobSpec` construction workspace-wide is a
//! full struct literal (no `Default`), so widening it would force
//! unrelated edits across 29 call sites for a one-gate need — the same
//! "least churn" reasoning Task 1 used to keep `HProbe` a parallel
//! `Drive` field instead of widening `Probe`/`EComponent` (see
//! `.superpowers/sdd/fs42a-task-1-report.md`).
//!
//! # Measurement method
//!
//! **V(t)**: column of `Ez` probes from the ground plane (`k=0`) up to
//! (excluding) the trace plane (`k=k_trace`), summed × `dz`, at a plane
//! 2.5 guided wavelengths downstream of the port. Symmetric stripline ⇒
//! the lower-half ground-to-trace integral equals the line voltage once
//! the launch transient has settled into the guided quasi-TEM mode (the
//! port itself excites only the upper-half `Ez` cell, `k=k_trace`, but by
//! structural symmetry the propagating mode is symmetric about the trace
//! far from the port — same reasoning `engine-stripline-eeff-001` relies
//! on for its single below-trace phase probe).
//!
//! **I(t)**: a rectangular Ampère loop in the `(y, z)` cross-section at
//! the *same* `i` index, one `dz` tall (straddling the trace plane
//! exactly: `Hy` at `k_trace−1` below, `Hy` at `k_trace` above) and
//! spanning the trace width plus a guard margin (`Hz` at the two side
//! legs). This loop is not an approximation of the FDTD curl — it *is*
//! the curl: `Ex(i, j, k_trace)` sits at the electrical centre of exactly
//! this `Hy`/`Hz` rectangle (Ampère-law surface bounded by the loop is
//! pierced only by the trace's own PEC surface current), and summing
//! adjacent unit loops telescopes the *interior* `Hz` terms away, leaving
//! only the two outermost side samples:
//!
//! ```text
//! I = Δy·[Σ_j Hy(j, k_trace−1) − Σ_j Hy(j, k_trace)]
//!   + Δz·[Hz(j_hi, k_trace) − Hz(j_lo−1, k_trace)]
//! ```
//! (CCW loop viewed with `+x` out of the page, `y` right, `z` up — thumb
//! along `+x` for a CCW circulation, matching current flowing toward the
//! load.) Because Ampère's law is exact and the enclosed current is
//! conserved (`∇·(J + ∂D/∂t) = 0` identically), the guard margin is not
//! load-bearing for correctness — any loop that fully encircles the trace
//! without crossing it or another conductor measures the same `I`. The
//! margin is chosen only to sit comfortably clear of the trace edge and
//! the box's side walls.
//!
//! # Staggering
//!
//! - **Time**: `yee-compute` records the co-indexed E/H probe pair at
//!   *different* times — sample index `m`'s `Ez` is `t=(m+1)·Δt` (recorded
//!   right after the E half-step, after `step` increments), its `Hy`/`Hz`
//!   is `t=(m+½)·Δt` (written by `update_h` at the top of the *same*
//!   iteration — see `HProbe`'s doc comment in `yee-compute/src/drive.rs`,
//!   landed by Task 1). The single-bin DFT below uses each series' own
//!   true sample time in its phasor exponent, so this is handled exactly
//!   (no interpolation, no residual error).
//! - **Space**: `Ez`'s `i` index is an integer-`x` node; `Hy`/`Hz`'s `i`
//!   index (same numeric value, `i_probe`) is the half-`x` node a half
//!   cell downstream. Unlike the time offset this is not corrected —
//!   quantified instead: at `f0`/`λ_g` below, the induced phase error is
//!   `β·Δx/2 = π·dx/λ_g ≈ 0.019 rad` (≈ 1.1°), a `cos` magnitude error
//!   `< 0.02 %` — three orders under this gate's tolerance.
//!
//! # Closed form
//!
//! `k = tanh(πw/2b)`, `k′ = sech(πw/2b) = √(1−k²)`,
//! `Z₀ = (η₀/4√ε_r)·K(k′)/K(k)`, `K` via the ~10-line AGM iteration
//! (Abramowitz & Stegun 17.6). **Note on the `k`/`k′` labelling**: the
//! FS.4.2a design spec
//! (`docs/superpowers/specs/2026-07-23-fs4-2a-stripline-z0-design.md`)
//! writes `k = sech(πw/2b), k′ = tanh(πw/2b)` — swapped from the
//! standard/Wikipedia "Stripline" convention used here. Plugging the
//! spec's own labels into `Z₀ = K(k′)/K(k)` makes `Z₀` *increase* with
//! `w/b`, which contradicts basic transmission-line physics (wider trace
//! ⇒ more capacitance per length ⇒ *lower* `Z₀`, `Z₀ → ∞` as `w → 0`) and
//! disagrees with the Wheeler/Pozar fit by tens of percent instead of the
//! expected ≲1 % (caught by `wheeler_fit_cross_checks_the_exact_form`
//! below, *before* running any FDTD — see the Task 2 report for the
//! numeric check). This file uses the convention verified to give the
//! physically correct trend and to match the Wheeler fit to <0.1% at this
//! fixture's `w/b`.
//!
//! `#[ignore]`'d (multi-minute release run):
//!
//! ```bash
//! cargo test -p yee-engine --release --test stripline_z0 -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_compute::{
    Boundary, CpuFdtd, Drive, EComponent, FdtdSpec, Fields, HComponent, HProbe, Materials, Probe,
    ResistivePort, Waveform,
};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Stackup, Substrate};
use yee_voxel::{VoxelOptions, voxelize_stackup};

const EPS_R: f64 = 2.2;
/// Ground-to-ground spacing b — 16 cells at `DX_M` (ADR-0215/0221 lesson:
/// confined lidded stripline modes need b >= ~16 cells).
const B_M: f64 = 3.2e-3;
/// Trace width — 13 cells at `DX_M`; chosen (not guessed) by solving
/// `z0_stripline_exact(EPS_R, w, B_M) = 50 Ω` for `w` and rounding to the
/// nearest whole cell (13 cells -> Z0_exact ~= 50.7 Ω, see the Task 2
/// report for the bisection).
const W_M: f64 = 2.6e-3;
const DX_M: f64 = 0.2e-3;
const F0_HZ: f64 = 6.0e9;
const C0_M_S: f64 = 299_792_458.0;
/// Box-mode hygiene (same TE10-cutoff discipline as
/// `engine-stripline-eeff-001`): box width = W_M + 2*MARGIN_CELLS*DX_M =
/// 10.6 mm -> f_c = c/(2*w_box*sqrt(EPS_R)) ~= 9.5 GHz, safely above the
/// drive band's top (~F0_HZ + 2*bw ~= 8.4 GHz at bw = 0.4*F0_HZ below).
const MARGIN_CELLS: usize = 20;
/// Ampère-loop guard beyond the trace's y-extent (not load-bearing for
/// correctness — see the module doc; picked to clear the trace edge by
/// about the b/pi transverse decay scale, ADR-0215).
const GUARD_CELLS: usize = 5;
const PORT_R_OHM: f64 = 50.0;

/// Arithmetic-geometric mean (Abramowitz & Stegun 17.6): `K(k) = (pi/2) /
/// AGM(1, sqrt(1-k^2))`.
fn agm(a0: f64, b0: f64) -> f64 {
    let (mut a, mut b) = (a0, b0);
    for _ in 0..60 {
        let a_next = 0.5 * (a + b);
        let b_next = (a * b).sqrt();
        if (a_next - b_next).abs() <= 1e-15 * a_next.max(1e-300) {
            return a_next;
        }
        a = a_next;
        b = b_next;
    }
    a
}

/// Exact conformal-mapping Z0 of a zero-thickness symmetric stripline —
/// see the module doc's "Closed form" section for the k/k' convention
/// note.
fn z0_stripline_exact(eps_r: f64, w_m: f64, b_m: f64) -> f64 {
    const ETA0: f64 = 376.730_313_668;
    let x = PI * w_m / (2.0 * b_m);
    let k = x.tanh();
    let kp = 1.0 / x.cosh(); // sech(x) = sqrt(1 - k^2)
    // Z0 = (eta0/4 sqrt(er)) * K(k')/K(k) = (eta0/4 sqrt(er)) * AGM(1,k')/AGM(1,k)
    ETA0 / (4.0 * eps_r.sqrt()) * agm(1.0, kp) / agm(1.0, k)
}

/// Wheeler/Cohn closed-form fit (Pozar §3.8; valid `W/b >= 0.35`, zero
/// strip thickness) — a cross-check on `z0_stripline_exact`'s AGM/K
/// implementation, not itself the pinned reference.
fn z0_stripline_wheeler_fit(eps_r: f64, w_m: f64, b_m: f64) -> f64 {
    30.0 * PI / eps_r.sqrt() / (w_m / b_m + 0.441)
}

#[test]
fn wheeler_fit_cross_checks_the_exact_form() {
    // Cheap, always-on: validates the AGM/elliptic-K implementation
    // itself against an independent closed form, with no FDTD involved.
    let z_exact = z0_stripline_exact(EPS_R, W_M, B_M);
    let z_fit = z0_stripline_wheeler_fit(EPS_R, W_M, B_M);
    let err = (z_fit - z_exact).abs() / z_exact;
    eprintln!(
        "engine-stripline-z0-001: closed-form cross-check: exact = {z_exact:.4} Ohm, \
         Wheeler fit = {z_fit:.4} Ohm (err {:.4} %)",
        err * 100.0
    );
    assert!(
        err < 0.01,
        "exact-form/Wheeler-fit cross-check exceeds 1 % ({:.4} %) — the AGM/k \
         convention is probably wrong, see the module doc's k/k' note",
        err * 100.0
    );
    // Physical sanity the module doc argues for: Z0 must fall as w/b
    // grows (wider trace -> more capacitance -> lower Z0).
    assert!(
        z0_stripline_exact(EPS_R, W_M * 1.5, B_M) < z_exact,
        "Z0 must decrease as the trace widens"
    );
}

#[test]
#[ignore = "slow: multi-minute release FDTD; engine-stripline-z0-001 gate (FS.4.2a) — run with --release --ignored"]
fn stripline_z0_matches_the_exact_closed_form() {
    let z0_exact = z0_stripline_exact(EPS_R, W_M, B_M);

    let lam_g = C0_M_S / (F0_HZ * EPS_R.sqrt());
    let l_m = 8.0 * lam_g;
    let traces = vec![Polygon::rect(0.0, 0.0, l_m, W_M)];
    let bbox = BBox::from_polygons(&traces);
    let layout = Layout {
        substrate: Substrate {
            // Unused by the stackup path; kept for the Layout contract.
            eps_r: EPS_R,
            height_m: B_M,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        },
        traces,
        ports: vec![PortRef {
            at: Point2::new(0.5e-3, W_M / 2.0),
            width_m: W_M,
            ref_impedance_ohm: PORT_R_OHM,
        }],
        bbox,
    };

    let stack = Stackup::symmetric_stripline(EPS_R, B_M);
    let model = voxelize_stackup(
        &layout,
        &stack,
        0,
        &VoxelOptions {
            dx_m: DX_M,
            xy_margin_cells: MARGIN_CELLS,
            air_above_cells: 0, // lidded: ignored
        },
    );
    let (nx, ny, nz) = model.dims;
    let (_i_drive, j_strip, k_trace) = model.port_cells[0];
    let dt = model.grid.dt;
    let dx = model.dx_m;
    eprintln!(
        "engine-stripline-z0-001: grid {nx}x{ny}x{nz}, trace at k = {k_trace} (b = {nz} \
         cells), w/b = {:.4}, Z0_exact = {z0_exact:.4} Ohm, L = {:.1} mm",
        W_M / B_M,
        l_m * 1e3
    );

    let mut fdtd_spec = FdtdSpec::vacuum(nx, ny, nz, DX_M);
    fdtd_spec.dt = dt;

    let materials = Materials {
        eps_r_cells: model
            .grid
            .eps_r_cells
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        pec_mask_ex: model
            .grid
            .pec_mask_ex
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        pec_mask_ey: model
            .grid
            .pec_mask_ey
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        ..Materials::default()
    };

    // Measurement plane: 2.5 guided wavelengths downstream of the port,
    // past the launch transient (same hygiene as
    // engine-stripline-eeff-001).
    let x_probe = 2.5 * lam_g;
    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| (((xp - x0) / dx).round() as isize).clamp(0, nx as isize - 1) as usize;
    let i_probe = i_for(x_probe);

    // Ampère-loop y-span: the trace's own y-extent (via the same
    // physical-coordinate mapping the voxelizer used) plus GUARD_CELLS.
    let y0 = layout.bbox.min.y - MARGIN_CELLS as f64 * dx;
    let j_for = |yp: f64| (((yp - y0) / dx).round() as isize).clamp(0, ny as isize) as usize;
    let j_trace_lo = j_for(0.0);
    let j_trace_hi = j_for(W_M);
    let j_lo = j_trace_lo.saturating_sub(GUARD_CELLS);
    let j_hi = (j_trace_hi + GUARD_CELLS).min(ny - 1);
    assert!(
        j_lo >= 1,
        "loop's left Hz leg (j_lo-1) underflows: j_lo = {j_lo}"
    );
    assert!(
        j_hi < ny,
        "loop's right Hz leg (j_hi) is out of the Hz y-range: j_hi = {j_hi}, ny = {ny}"
    );
    let n_j = j_hi - j_lo + 1;

    // Time gate: stop before the far-end reflection reaches the probe
    // plane.
    let v_p_ref = C0_M_S / EPS_R.sqrt();
    let x_drive = 0.5e-3;
    let t_refl = ((l_m - x_drive) + (l_m - x_probe)) / v_p_ref;
    let gate_steps = (0.9 * t_refl / dt) as usize;
    let n_steps = gate_steps + 200;

    let bw = 0.4 * F0_HZ;
    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * bw)) / dt).ceil() as usize;

    let mut drive = Drive::default();
    drive.ports.push(ResistivePort {
        cell: model.port_cells[0],
        resistance: PORT_R_OHM,
        waveform: Waveform::GaussianPulse {
            v0: 1.0,
            f0: F0_HZ,
            bw,
            t0_steps,
        },
    });
    // V(t): ground (k=0) up to (excluding) the trace (k=k_trace).
    for k in 0..k_trace {
        drive.probes.push(Probe {
            component: EComponent::Ez,
            cell: (i_probe, j_strip, k),
        });
    }
    // I(t): bottom leg (Hy just below the trace), then top leg (Hy just
    // above), then the two Hz side legs — this fixed layout is what the
    // DFT loop below slices back apart.
    for j in j_lo..=j_hi {
        drive.h_probes.push(HProbe {
            component: HComponent::Hy,
            cell: (i_probe, j, k_trace - 1),
        });
    }
    for j in j_lo..=j_hi {
        drive.h_probes.push(HProbe {
            component: HComponent::Hy,
            cell: (i_probe, j, k_trace),
        });
    }
    drive.h_probes.push(HProbe {
        component: HComponent::Hz,
        cell: (i_probe, j_lo - 1, k_trace),
    });
    drive.h_probes.push(HProbe {
        component: HComponent::Hz,
        cell: (i_probe, j_hi, k_trace),
    });

    let fields = Fields::zero(&fdtd_spec);
    let mut engine = CpuFdtd::with_drive(fdtd_spec, fields, materials, Boundary::PecBox, drive);
    engine.step_n(n_steps);

    // Time-gated single-bin DFT at f0, each series phase-referenced at
    // its own true sample time (see the module doc's "Staggering" note):
    // Ez samples at t=(m+1)*dt, Hy/Hz samples at t=(m+1/2)*dt.
    let omega = 2.0 * PI * F0_HZ;
    let gate = gate_steps.min(n_steps);
    let ez = engine.probe_series();
    let hp = engine.h_probe_series();

    let mut v_acc = [0.0_f64; 2];
    let mut i_acc = [0.0_f64; 2];
    for m in 0..gate {
        let v_t: f64 = (0..k_trace).map(|k| ez[k][m]).sum::<f64>() * dx;
        let t_v = (m + 1) as f64 * dt;
        let (sv, cv) = (omega * t_v).sin_cos();
        v_acc[0] += v_t * cv;
        v_acc[1] -= v_t * sv;

        let hy_bot: f64 = (0..n_j).map(|jj| hp[jj][m]).sum();
        let hy_top: f64 = (0..n_j).map(|jj| hp[n_j + jj][m]).sum();
        let hz_left = hp[2 * n_j][m];
        let hz_right = hp[2 * n_j + 1][m];
        let i_t = dx * (hy_bot - hy_top) + dx * (hz_right - hz_left);
        let t_i = (m as f64 + 0.5) * dt;
        let (si, ci) = (omega * t_i).sin_cos();
        i_acc[0] += i_t * ci;
        i_acc[1] -= i_t * si;
    }
    let v_mag = (v_acc[0] * v_acc[0] + v_acc[1] * v_acc[1]).sqrt();
    let i_mag = (i_acc[0] * i_acc[0] + i_acc[1] * i_acc[1]).sqrt();

    // Guard against a silent-zero probe wiring bug turning the ratio into
    // an accidental 0/0 or a near-zero/near-zero coincidence. V and I are
    // raw gated-DFT phasor sums (not amplitude-normalized), so their
    // absolute scales differ by design — not comparable to each other,
    // only each to a "this is clearly not zero" floor. Measured
    // (2026-07-23, this fixture): |V| ~= 2.33, |I| ~= 0.0466; the floor
    // below sits three orders under the smaller of the two.
    assert!(
        v_mag.is_finite() && v_mag > 1e-3,
        "V phasor magnitude is not non-trivial: {v_mag}"
    );
    assert!(
        i_mag.is_finite() && i_mag > 1e-3,
        "I phasor magnitude is not non-trivial: {i_mag}"
    );

    let z0_meas = v_mag / i_mag;
    let rel_err = (z0_meas - z0_exact).abs() / z0_exact;
    eprintln!(
        "  Z0_meas = {z0_meas:.4} Ohm vs exact {z0_exact:.4} Ohm -> err {:.3} % \
         (|V| = {v_mag:.3e}, |I| = {i_mag:.3e}, {n_steps} steps, gate {gate_steps})",
        rel_err * 100.0
    );
    assert!(
        rel_err <= 0.05,
        "engine-stripline-z0-001 FAILED: Z0_meas = {z0_meas:.4} Ohm vs exact \
         {z0_exact:.4} Ohm (err {:.3} % > 5 %)",
        rel_err * 100.0
    );
}
