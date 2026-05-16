"""Tests for `yee.FreqRange`."""

import numpy as np
import pytest

import yee


def test_freqrange_constructs():
    f = yee.FreqRange(1.0e9, 2.0e9, 21)
    assert f.start_hz == 1.0e9
    assert f.stop_hz == 2.0e9
    assert f.n_points == 21


def test_freqrange_iter_endpoints_exact():
    f = yee.FreqRange(1.0e9, 2.0e9, 3)
    pts = f.iter()
    assert pts.shape == (3,)
    assert pts.dtype == np.float64
    # Endpoint pinning is part of the Rust contract — must round-trip bit-for-bit.
    assert pts[0] == 1.0e9
    assert pts[-1] == 2.0e9


def test_freqrange_rejects_inverted_band():
    with pytest.raises(ValueError):
        yee.FreqRange(2.0e9, 1.0e9, 5)


def test_freqrange_rejects_zero_points():
    with pytest.raises(ValueError):
        yee.FreqRange(1.0e9, 2.0e9, 0)


def test_freqrange_rejects_non_finite():
    with pytest.raises(ValueError):
        yee.FreqRange(float("nan"), 2.0e9, 5)
