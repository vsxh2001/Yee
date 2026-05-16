"""Tests for yee.bo_minimize (PyO3 wrapper around yee-surrogate::bo).

The 1-D deceptive objective `(x - 3)^2 + sin(5x)` has its global minimum
near `(x = 3.422, y = -0.807)` on `[0, 6]` and a local minimum near
`x = 4`; with the documented BoConfig defaults the EI maximizer should
clear `y < 0` within ~25 evaluations.
"""

import math

import numpy as np

from yee import BoConfig, BoResult, bo_minimize


def deceptive_1d(x):
    return (x[0] - 3.0) ** 2 + math.sin(5.0 * x[0])


def test_bo_minimize_deceptive_1d():
    cfg = BoConfig(n_initial=5, n_iters=20, n_candidates=1024, xi=0.01, seed=0xC0FFEE)
    result = bo_minimize(deceptive_1d, [(0.0, 6.0)], cfg)

    assert isinstance(result, BoResult)

    # Global minimum near (x = 3.422, y = -0.807).
    assert result.y_best < 0.0, (
        f"BO should clear y<0 within budget, got {result.y_best}"
    )
    assert abs(result.x_best[0] - 3.422) < 0.1, (
        f"x_best far from minimum, got {result.x_best[0]}"
    )

    # n_initial + n_iters = 25 chronological evaluations.
    assert result.history_x.shape == (25, 1)
    assert result.history_y.shape == (25,)
    assert result.history_x.dtype == np.float64
    assert result.history_y.dtype == np.float64

    # y_best matches the minimum of history_y.
    assert result.y_best == float(np.min(result.history_y))


def test_bo_minimize_default_config():
    result = bo_minimize(deceptive_1d, [(0.0, 6.0)])
    # Default config: n_initial=5, n_iters=20 -> 25 evals.
    assert result.history_y.shape == (25,)
    assert result.history_x.shape == (25, 1)


def test_bo_minimize_2d_rosenbrock():
    def rosen(x):
        return (1.0 - x[0]) ** 2 + 100.0 * (x[1] - x[0] ** 2) ** 2

    cfg = BoConfig(n_initial=10, n_iters=40, seed=42)
    result = bo_minimize(rosen, [(-2.0, 2.0), (-2.0, 2.0)], cfg)

    # 2-D Rosenbrock has min at (1, 1) with y = 0. BO with budget = 50
    # should reach y < 5.0 even though the valley is curved; this is a
    # signal-not-precision test.
    assert result.y_best < 5.0
    assert result.history_x.shape == (50, 2)
    assert result.history_y.shape == (50,)


def test_bo_config_repr_and_accessors():
    cfg = BoConfig(n_initial=7, n_iters=11, n_candidates=128, xi=0.05, seed=7)
    assert cfg.n_initial == 7
    assert cfg.n_iters == 11
    assert cfg.n_candidates == 128
    assert cfg.xi == 0.05
    assert cfg.seed == 7
    r = repr(cfg)
    assert "BoConfig" in r
    assert "n_initial=7" in r
