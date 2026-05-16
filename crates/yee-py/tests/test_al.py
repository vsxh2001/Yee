"""Tests for yee.active_learn (PyO3 wrapper around yee-surrogate::al)."""

import math

import numpy as np

from yee import AlConfig, AlResult, GaussianProcess, active_learn


def test_al_recovers_sin():
    cfg = AlConfig(n_initial=5, n_iters=20, seed=42)
    result = active_learn(lambda x: math.sin(x[0]), [(0.0, 2 * math.pi)], cfg)

    assert isinstance(result, AlResult)
    assert result.history_x.shape == (25, 1)
    assert result.history_y.shape == (25,)
    assert result.history_x.dtype == np.float64
    assert result.history_y.dtype == np.float64

    gp = result.final_gp()
    assert isinstance(gp, GaussianProcess)

    # GP should predict near sin(theta) at unseen points after 25 evals
    # along a single period. Tolerance is loose but well within reach for
    # 25 well-chosen samples of a smooth scalar function.
    for x in [0.3, 1.2, 2.5, 4.1, 5.7]:
        pred = gp.predict_mean(np.array([x]))
        assert abs(pred - math.sin(x)) < 0.01, (
            f"GP prediction at x={x}: {pred:.4f} vs sin(x)={math.sin(x):.4f}"
        )


def test_al_default_config():
    result = active_learn(lambda x: float(x[0] ** 2), [(-1.0, 1.0)])
    # Default cfg: n_initial=5, n_iters=20 -> 25 evals.
    assert result.history_x.shape == (25, 1)
    assert result.history_y.shape == (25,)


def test_al_config_repr_and_accessors():
    cfg = AlConfig(n_initial=4, n_iters=8, n_candidates=256, seed=13)
    assert cfg.n_initial == 4
    assert cfg.n_iters == 8
    assert cfg.n_candidates == 256
    assert cfg.seed == 13
    r = repr(cfg)
    assert "AlConfig" in r
    assert "n_initial=4" in r


def test_al_2d_quadratic():
    def quad(x):
        return float((x[0] - 0.5) ** 2 + (x[1] + 0.3) ** 2)

    cfg = AlConfig(n_initial=6, n_iters=10, seed=99)
    result = active_learn(quad, [(-1.0, 1.0), (-1.0, 1.0)], cfg)
    assert result.history_x.shape == (16, 2)
    assert result.history_y.shape == (16,)

    gp = result.final_gp()
    # GP should at least approximately track the parabola at its minimum.
    pred_min = gp.predict_mean(np.array([0.5, -0.3]))
    assert pred_min < 0.2
