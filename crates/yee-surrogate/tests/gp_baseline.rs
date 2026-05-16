//! gp-001 — Gaussian-process baseline validation.
//!
//! Generates 25 training points on y = sin(x), x ∈ [0, 2π], holds out 10
//! interleaved test points, and asserts that:
//!
//! 1. GP mean RMS error on the held-out set < 0.05.
//! 2. GP RMS error beats the NearestNeighbor baseline RMS error on the same
//!    held-out set.
//! 3. GP posterior variance at a training point is within 1e-3 of sigma_n^2.
//! 4. GP posterior variance at a midpoint between two training samples is
//!    strictly greater than sigma_n^2 (uncertainty grows away from data).
//!
//! Loose tolerances by design: this is the walking-skeleton gate for the
//! GP backend, not a tight regression test.

use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;
use yee_surrogate::{Dataset, GaussianProcess, GpSurrogate, NearestNeighbor, Sample, Surrogate};

const N_TRAIN: usize = 25;
const N_TEST: usize = 10;
const LENGTH_SCALE: f64 = 0.5;
const SIGMA_F: f64 = 1.0;
const SIGMA_N: f64 = 1e-4;

fn linspace(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    if n == 1 {
        return vec![lo];
    }
    let step = (hi - lo) / ((n - 1) as f64);
    (0..n).map(|i| lo + step * (i as f64)).collect()
}

fn rms(errs: &[f64]) -> f64 {
    (errs.iter().map(|e| e * e).sum::<f64>() / (errs.len() as f64)).sqrt()
}

#[test]
fn gp_beats_nn_on_sin_one_d() {
    // Training grid: 25 equally spaced points on [0, 2π].
    let x_train: Vec<f64> = linspace(0.0, 2.0 * std::f64::consts::PI, N_TRAIN);
    let y_train: Vec<f64> = x_train.iter().map(|x| x.sin()).collect();

    // Test grid: 10 points shifted halfway between consecutive training
    // samples (so we are genuinely interpolating, not sitting on training x).
    let step = (2.0 * std::f64::consts::PI) / ((N_TRAIN - 1) as f64);
    let x_test: Vec<f64> = (0..N_TEST)
        .map(|i| 0.5 * step + (i as f64) * (2.0 * std::f64::consts::PI - step) / (N_TEST as f64))
        .collect();
    let y_test: Vec<f64> = x_test.iter().map(|x| x.sin()).collect();

    // Fit the GP directly (scalar API).
    let x_mat = DMatrix::from_column_slice(N_TRAIN, 1, &x_train);
    let y_vec = DVector::from_row_slice(&y_train);
    let gp = GaussianProcess::fit(x_mat, y_vec, LENGTH_SCALE, SIGMA_F, SIGMA_N)
        .expect("GP fit on sin(x) should succeed");

    // Fit the NearestNeighbor baseline via the Surrogate trait on the same data.
    let mut ds = Dataset::new();
    for (xi, yi) in x_train.iter().zip(y_train.iter()) {
        ds.push(Sample {
            params: vec![*xi],
            output: vec![Complex64::new(*yi, 0.0)],
        });
    }
    let mut nn = NearestNeighbor::new();
    nn.train(&ds).expect("NN train should succeed");

    // Also exercise the GpSurrogate wrapper for trait coverage.
    let mut gp_surr = GpSurrogate::with_hyperparams(LENGTH_SCALE, SIGMA_F, SIGMA_N);
    gp_surr
        .train(&ds)
        .expect("GpSurrogate train should succeed");

    // Compute RMS for each.
    let mut gp_errs = Vec::with_capacity(N_TEST);
    let mut gp_surr_errs = Vec::with_capacity(N_TEST);
    let mut nn_errs = Vec::with_capacity(N_TEST);
    for (xi, yi) in x_test.iter().zip(y_test.iter()) {
        let xs = DVector::from_row_slice(&[*xi]);
        let mu = gp.predict_mean(&xs);
        gp_errs.push(mu - yi);

        let mu_surr = gp_surr.predict(&[*xi]).expect("GpSurrogate predict")[0].re;
        gp_surr_errs.push(mu_surr - yi);

        let nn_pred = nn.predict(&[*xi]).expect("NN predict")[0].re;
        nn_errs.push(nn_pred - yi);
    }
    let gp_rms = rms(&gp_errs);
    let gp_surr_rms = rms(&gp_surr_errs);
    let nn_rms = rms(&nn_errs);
    println!(
        "gp_baseline: GP RMS = {gp_rms:.6e}, GpSurrogate RMS = {gp_surr_rms:.6e}, NN RMS = {nn_rms:.6e}, GP/NN ratio = {:.3}",
        gp_rms / nn_rms
    );

    // (1) GP must hit absolute accuracy.
    assert!(gp_rms < 0.05, "GP RMS = {gp_rms} exceeds 0.05 budget");
    // The wrapper must reproduce the inherent path.
    assert!(
        (gp_rms - gp_surr_rms).abs() < 1e-12,
        "GpSurrogate ({gp_surr_rms}) disagrees with GaussianProcess ({gp_rms})"
    );
    // (2) GP must beat NN baseline.
    assert!(
        gp_rms < nn_rms,
        "GP RMS = {gp_rms} did not beat NN RMS = {nn_rms}"
    );

    // (3) Variance at a training point ≈ sigma_n^2.
    let x_train_pt = DVector::from_row_slice(&[x_train[5]]);
    let (_, var_train) = gp.predict(&x_train_pt);
    let sn2 = SIGMA_N * SIGMA_N;
    assert!(
        (var_train - sn2).abs() < 1e-3,
        "variance at training point = {var_train}, expected ≈ {sn2}"
    );

    // (4) Variance grows away from training data: midpoint between two
    // adjacent training samples carries more posterior uncertainty than the
    // training point itself.
    let mid = 0.5 * (x_train[5] + x_train[6]);
    let x_mid = DVector::from_row_slice(&[mid]);
    let (_, var_mid) = gp.predict(&x_mid);
    assert!(
        var_mid > sn2,
        "variance at midpoint = {var_mid} did not exceed sigma_n^2 = {sn2}"
    );
}
