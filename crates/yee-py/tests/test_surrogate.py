"""Tests for yee.GaussianProcess (PyO3 wrapper around yee-surrogate)."""

import numpy as np

from yee import GaussianProcess


def test_gp_fit_predict_sin():
    x_train = np.linspace(0.0, 2 * np.pi, 25).reshape(-1, 1)
    y_train = np.sin(x_train).flatten()
    gp = GaussianProcess.fit(
        x_train, y_train, length_scale=0.5, sigma_f=1.0, sigma_n=1e-4
    )

    # Predict mean at an interior point.
    x_star = np.array([1.7])
    mean = gp.predict_mean(x_star)
    assert abs(mean - np.sin(1.7)) < 0.05

    # Predict (mean, variance).
    m2, v2 = gp.predict(x_star)
    assert m2 == mean
    assert v2 >= 0.0


def test_gp_fit_ml_beats_default():
    x_train = np.linspace(0.0, 2 * np.pi, 25).reshape(-1, 1)
    y_train = np.sin(x_train).flatten()
    bad = GaussianProcess.fit(
        x_train, y_train, length_scale=5.0, sigma_f=0.1, sigma_n=1e-2
    )
    ml = GaussianProcess.fit_ml(
        x_train,
        y_train,
        initial_length_scale=5.0,
        initial_sigma_f=0.1,
        initial_sigma_n=1e-2,
    )

    x_test = np.array([0.3, 1.2, 2.4, 3.7, 5.1])
    bad_rms = float(
        np.sqrt(
            np.mean(
                [
                    (bad.predict_mean(np.array([xi])) - np.sin(xi)) ** 2
                    for xi in x_test
                ]
            )
        )
    )
    ml_rms = float(
        np.sqrt(
            np.mean(
                [
                    (ml.predict_mean(np.array([xi])) - np.sin(xi)) ** 2
                    for xi in x_test
                ]
            )
        )
    )

    assert ml_rms < bad_rms
    assert ml_rms < 0.01
