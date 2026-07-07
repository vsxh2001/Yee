//! Gate `compute-017` (R.0b, ADR-0202): the **resistive-sheet conductor
//! boundary** on masked `E_x`/`E_y` edges.
//!
//! 1. **Dissipation with the right sign**: in a driven PEC box containing
//!    a masked sheet, the late-time field energy with `R_s > 0` must fall
//!    strictly below the plain-PEC-mask run — a wrong sheet sign
//!    self-amplifies and fails this structurally.
//! 2. **PEC degeneracy**: `R_s = 0` (and `None`) reproduce the PEC-mask
//!    run **bit-exactly** — the sheet is a pure extension.
//!
//! Fast (a 24³ box, 400 steps, three runs); runs in the default test pass.

use yee_compute::{Boundary, CpuFdtd, FdtdSpec, Fields, Materials};

const N: usize = 24;
const STEPS: usize = 400;

fn spec() -> FdtdSpec {
    FdtdSpec::vacuum(N, N, N, 1.0e-3)
}

/// A horizontal full-septum sheet: E_x and E_y masked across the whole
/// box cross-section at mid-height, so every bounce crosses the sheet
/// (nothing diffracts around it).
fn sheet_materials(sheet_r_ohm: Option<f64>) -> Materials {
    let s = spec();
    let exd = s.ex_dims();
    let eyd = s.ey_dims();
    let mut ex = vec![false; exd.0 * exd.1 * exd.2];
    let mut ey = vec![false; eyd.0 * eyd.1 * eyd.2];
    let k = N / 2 + 3;
    for i in 0..exd.0 {
        for j in 0..exd.1 {
            ex[(i * exd.1 + j) * exd.2 + k] = true;
        }
    }
    for i in 0..eyd.0 {
        for j in 0..eyd.1 {
            ey[(i * eyd.1 + j) * eyd.2 + k] = true;
        }
    }
    Materials {
        pec_mask_ex: Some(ex),
        pec_mask_ey: Some(ey),
        sheet_r_ohm,
        ..Materials::default()
    }
}

fn run(sheet_r_ohm: Option<f64>) -> Fields {
    let s = spec();
    let mut stepper = CpuFdtd::with_config(
        s,
        Fields::zero(&s),
        sheet_materials(sheet_r_ohm),
        Boundary::PecBox,
    );
    let c = N / 2;
    for n in 0..STEPS {
        // Gaussian Ez pulse a few cells below the sheet.
        stepper.step_with_gaussian_ez((c, c, c - 4), 40.0 * s.dt, 12.0 * s.dt);
        let _ = n;
    }
    stepper.fields().clone()
}

/// Ringing (propagating-wave) energy proxy: ‖H‖². The soft Gaussian
/// source leaves a static E remnant that dominates ‖E‖² identically in
/// both runs and that no sheet can dissipate (a static field carries no
/// surface current) — the magnetic energy isolates the wave the sheet
/// actually acts on.
fn ringing_energy(f: &Fields) -> f64 {
    f.hx.iter()
        .chain(f.hy.iter())
        .chain(f.hz.iter())
        .map(|v| v * v)
        .sum()
}

#[test]
fn sheet_dissipates_and_r_zero_degenerates_to_pec() {
    let pec = run(None);
    let r_zero = run(Some(0.0));
    // eta0/2 ~ 188 ohm: near the matched-absorber point for a sheet,
    // maximizing per-pass absorption — the strongest test of the sign.
    let lossy = run(Some(188.0));

    // 2. R_s = 0 is the PEC mask, bit-exactly.
    assert_eq!(pec.ex, r_zero.ex, "R_s = 0 must be bit-exact PEC (ex)");
    assert_eq!(pec.hz, r_zero.hz, "R_s = 0 must be bit-exact PEC (hz)");

    // 1. The sheet drains the ringing wave energy from the closed box
    //    (measured 0.56x at 400 steps; gated at 0.75x for margin).
    let e_pec = ringing_energy(&pec);
    let e_lossy = ringing_energy(&lossy);
    assert!(
        e_lossy < 0.75 * e_pec,
        "sheet must dissipate the ringing energy: lossy {e_lossy:.3e} vs pec {e_pec:.3e}"
    );
    assert!(e_lossy > 0.0, "fields must not vanish entirely");
}
