//! `filt-001` — synthesized-response spec-mask gate.
//!
//! A Chebyshev 0.5 dB bandpass synthesized at an adequate order must satisfy
//! its own mask over a swept band (PASS): passband ripple ≤ 0.5 dB, in-band
//! return loss ≥ spec, and stopband rejection ≥ the required minimum at every
//! mask point. A deliberately-too-low order is the negative control and must
//! FAIL the same mask.

use yee_filter::{Approximation, FilterSpec, Response, SpecMask, check_mask, synthesize};

/// Centre 2 GHz, 10% fractional bandwidth → band edges 1.9 / 2.1 GHz.
const F0: f64 = 2.0e9;
const FBW: f64 = 0.10;

fn mask() -> SpecMask {
    SpecMask {
        passband_ripple_db: 0.5,
        return_loss_db: 9.0,
        // 2.4 GHz maps to Ω ≈ 3.67 under the bandpass transform; an order-5
        // 0.5 dB Chebyshev gives strong rejection there, an order-2 does not.
        stopband: vec![(2.4e9, 40.0)],
    }
}

fn spec_with_order(order: usize) -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: F0,
        fbw: FBW,
        order: Some(order),
        z0_ohm: 50.0,
        mask: mask(),
    }
}

/// A swept band: dense across the passband (for ripple/RL) plus the f0 anchor.
fn swept_band() -> Vec<f64> {
    let mut f = Vec::new();
    let lo = 1.85e9;
    let hi = 2.15e9;
    let steps = 601;
    for i in 0..steps {
        f.push(lo + (hi - lo) * (i as f64) / ((steps - 1) as f64));
    }
    f
}

#[test]
fn adequate_order_passes() {
    let spec = spec_with_order(5);
    let proj = synthesize(&spec);
    let report = check_mask(&proj, &swept_band());
    assert!(
        report.pass,
        "expected PASS at order 5; ripple={:.3} dB, RL={:.3} dB, failures={:?}",
        report.worst_passband_ripple_db, report.worst_return_loss_db, report.failures
    );
    // Sanity on the graded quantities.
    assert!(
        report.worst_passband_ripple_db <= 0.5 + 1e-9,
        "passband ripple {:.3} dB should be within 0.5 dB",
        report.worst_passband_ripple_db
    );
    assert!(
        report.worst_return_loss_db >= 9.0 - 1e-9,
        "in-band return loss {:.3} dB should be >= 9 dB",
        report.worst_return_loss_db
    );
    // The single stopband point should be met with margin.
    let (_f, achieved, required, met) = report.stopband[0];
    assert!(
        met && achieved >= required,
        "stopband not met at order 5: achieved {achieved:.3} dB < required {required:.3} dB"
    );
}

#[test]
fn too_low_order_fails() {
    let spec = spec_with_order(2);
    let proj = synthesize(&spec);
    let report = check_mask(&proj, &swept_band());
    assert!(
        !report.pass,
        "expected FAIL at order 2 (negative control); report={report:?}"
    );
    // The failure must be the stopband rejection, not the passband.
    let (_f, achieved, required, met) = report.stopband[0];
    assert!(
        !met && achieved < required,
        "order-2 stopband should be under-rejected: achieved {achieved:.3} dB vs required {required:.3} dB"
    );
}

/// `synthesize` with `order: None` estimates the order from the stopband mask
/// and the result must PASS its own mask.
#[test]
fn estimated_order_passes() {
    let spec = FilterSpec {
        order: None,
        ..spec_with_order(0)
    };
    let proj = synthesize(&spec);
    assert!(
        proj.prototype.order() >= 1,
        "estimated order should be >= 1, got {}",
        proj.prototype.order()
    );
    let report = check_mask(&proj, &swept_band());
    assert!(
        report.pass,
        "estimated-order design should pass its own mask; order={}, failures={:?}",
        proj.prototype.order(),
        report.failures
    );
}
