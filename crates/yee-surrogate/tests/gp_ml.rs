//! gp-ml — Marginal-likelihood hyperparameter optimization validation.
//!
//! Verifies that `GaussianProcess::fit_ml` recovers tighter posterior fit on
//! `y = sin(x)` than a deliberately bad hand-tuned hyperparameter set. The
//! same `(length_scale, sigma_f, sigma_n)` seed feeds both paths so the only
//! difference is whether `fit_ml` actually does its job.

use yee_surrogate::{GaussianProcess, MlFitConfig};

#[test]
fn fit_ml_beats_default_hyperparams_on_sin_one_d() {
    // Same setup as gp_baseline: y = sin(x) on x ∈ [0, 2π], 25 training, 10 test.
    let n_train = 25;
    let mut x_train = nalgebra::DMatrix::<f64>::zeros(n_train, 1);
    let mut y_train = nalgebra::DVector::<f64>::zeros(n_train);
    for i in 0..n_train {
        let xi = (i as f64) * 2.0 * std::f64::consts::PI / (n_train as f64 - 1.0);
        x_train[(i, 0)] = xi;
        y_train[i] = xi.sin();
    }
    let mut x_test = nalgebra::DMatrix::<f64>::zeros(10, 1);
    let mut y_test = nalgebra::DVector::<f64>::zeros(10);
    for i in 0..10 {
        let xi = 0.1 + (i as f64) * 0.6; // off-grid
        x_test[(i, 0)] = xi;
        y_test[i] = xi.sin();
    }

    // Hand-tuned baseline (deliberately bad starting hyperparams).
    let bad_gp = GaussianProcess::fit(x_train.clone(), y_train.clone(), 5.0, 0.1, 1e-2).unwrap();
    let bad_rms = test_rms(&bad_gp, &x_test, &y_test);

    // ML-optimized fit from the same starting point.
    let ml_cfg = MlFitConfig {
        initial_length_scale: 5.0,
        initial_sigma_f: 0.1,
        initial_sigma_n: 1e-2,
        ..Default::default()
    };
    let ml_gp = GaussianProcess::fit_ml(x_train, y_train, ml_cfg).unwrap();
    let ml_rms = test_rms(&ml_gp, &x_test, &y_test);

    println!(
        "gp_ml: bad_rms = {bad_rms:.6e}, ml_rms = {ml_rms:.6e}, \
         optimized (ℓ, σ_f, σ_n) = ({}, {}, {}), log_marginal_likelihood = {}",
        ml_gp.length_scale(),
        ml_gp.sigma_f(),
        ml_gp.sigma_n(),
        ml_gp.log_marginal_likelihood(),
    );

    assert!(
        ml_rms < bad_rms,
        "fit_ml ({ml_rms}) should beat bad hand-tuned ({bad_rms})"
    );
    assert!(ml_rms < 0.01, "fit_ml RMS should be tight, got {ml_rms}");
}

fn test_rms(
    gp: &GaussianProcess,
    x_test: &nalgebra::DMatrix<f64>,
    y_test: &nalgebra::DVector<f64>,
) -> f64 {
    let mut s = 0.0;
    for i in 0..x_test.nrows() {
        let xi = x_test.row(i).transpose();
        let y_pred = gp.predict_mean(&xi.into_owned());
        s += (y_pred - y_test[i]).powi(2);
    }
    (s / (x_test.nrows() as f64)).sqrt()
}
