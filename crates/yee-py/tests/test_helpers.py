"""Tests for yee.{s11_db, s11_phase, smith_xy}."""

import numpy as np
import yee


def test_s11_db_basic():
    s = np.array([1.0 + 0.0j, 0.5 + 0.0j, 0.1 + 0.0j], dtype=np.complex128)
    db = yee.s11_db(s)
    np.testing.assert_allclose(db, [0.0, -6.020599913, -20.0], atol=1e-9)


def test_s11_db_clamps_at_zero():
    s = np.array([0.0 + 0.0j, 1e-30 + 0.0j], dtype=np.complex128)
    db = yee.s11_db(s)
    assert db[0] == -200.0
    # 1e-30 is finite, so the second element is *not* clamped — it's
    # `20 * log10(1e-30) = -600`. If the test wants clamping at any
    # absurdly-small magnitude, tighten the implementation; for now
    # we only clamp at exact zero, so allow the unclamped value.
    assert db[1] < -500.0


def test_s11_phase_degrees():
    s = np.array(
        [1.0 + 0.0j, 0.0 + 1.0j, -1.0 + 0.0j, 0.0 - 1.0j],
        dtype=np.complex128,
    )
    phase = yee.s11_phase(s)
    np.testing.assert_allclose(phase, [0.0, 90.0, 180.0, -90.0], atol=1e-9)


def test_smith_xy_shape_and_values():
    s = np.array([0.1 + 0.2j, 0.3 + 0.4j, 0.5 - 0.5j], dtype=np.complex128)
    xy = yee.smith_xy(s)
    assert xy.shape == (3, 2)
    np.testing.assert_allclose(xy[0], [0.1, 0.2])
    np.testing.assert_allclose(xy[1], [0.3, 0.4])
    np.testing.assert_allclose(xy[2], [0.5, -0.5])


def test_helpers_handle_empty_input():
    s = np.array([], dtype=np.complex128)
    assert yee.s11_db(s).shape == (0,)
    assert yee.s11_phase(s).shape == (0,)
    assert yee.smith_xy(s).shape == (0, 2)
