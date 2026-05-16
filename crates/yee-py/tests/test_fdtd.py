"""Tests for yee.FdtdDriver (PyO3 wrapper around yee-fdtd::FdtdDriver).

The end-to-end short-dipole simulation runs ~5-10s on a release build;
it is not gated as `slow` because the cost is small enough to run on
every `pytest tests/` invocation.
"""

import numpy as np

from yee import FdtdDriver, FdtdDriverConfig, FdtdRadiationPattern


def test_short_dipole_radiates_sin_theta():
    """A short z-directed dipole should produce |E_theta| ~ sin(theta).

    Asserts:
      - normalized magnitude at theta = 90° is ~1.0,
      - nulls at theta = 0° and 180° are small.
    """
    cfg = FdtdDriverConfig(
        n_steps=800,
        dipole_center_cells=(30, 30, 30),
        dipole_length_cells=5,
        source_freq_hz=1.0e9,
        ntff_surface_pad_cells=4,
        cpml_thickness_cells=10,
    )
    driver = FdtdDriver(60, 60, 60, 5.0e-3, cfg)
    pattern = driver.run()

    assert isinstance(pattern, FdtdRadiationPattern)

    theta = np.asarray(pattern.theta_deg)
    e = np.asarray(pattern.e_theta_phi0)

    # Sanity: 5° steps from 0 to 180 inclusive.
    assert theta.shape == (37,)
    assert e.shape == (37,)
    assert theta[0] == 0.0
    assert theta[-1] == 180.0
    assert theta.dtype == np.float64
    assert e.dtype == np.float64

    # Peak should be at theta ~ 90° and very close to 1 by construction
    # (the Rust driver normalizes to max = 1).
    i90 = int(np.argmin(np.abs(theta - 90.0)))
    assert abs(e[i90] - 1.0) < 0.05, (
        f"peak should be ~1 at 90°, got {e[i90]}"
    )

    # Nulls at 0 and 180.
    i0 = int(np.argmin(np.abs(theta)))
    i180 = int(np.argmin(np.abs(theta - 180.0)))
    assert e[i0] < 0.15, f"expected null near theta=0, got {e[i0]}"
    assert e[i180] < 0.15, f"expected null near theta=180, got {e[i180]}"


def test_config_repr_smoke():
    """Construct a config with default ntff/cpml args; verify accessors."""
    cfg = FdtdDriverConfig(
        n_steps=100,
        dipole_center_cells=(10, 10, 10),
        dipole_length_cells=3,
        source_freq_hz=1e9,
    )
    assert cfg is not None
    assert cfg.n_steps == 100
    assert cfg.dipole_center_cells == (10, 10, 10)
    assert cfg.dipole_length_cells == 3
    assert cfg.source_freq_hz == 1e9
    # Defaults applied.
    assert cfg.ntff_surface_pad_cells == 4
    assert cfg.cpml_thickness_cells == 10
    # __repr__ smoke.
    r = repr(cfg)
    assert "FdtdDriverConfig" in r
    assert "n_steps=100" in r


def test_driver_run_returns_normalized_pattern_small_grid():
    """Tiny grid + few steps: exercise the API end-to-end fast."""
    cfg = FdtdDriverConfig(
        n_steps=5,
        dipole_center_cells=(20, 20, 20),
        dipole_length_cells=3,
        source_freq_hz=1.0e9,
        ntff_surface_pad_cells=2,
        cpml_thickness_cells=8,
    )
    driver = FdtdDriver(40, 40, 40, 5.0e-3, cfg)
    pattern = driver.run()
    e = np.asarray(pattern.e_theta_phi0)
    # Either max is exactly 1.0 (normalized) or the field is identically
    # zero (degenerate at 5 steps). Either is acceptable behaviour for
    # this smoke test.
    if np.any(e > 0.0):
        assert abs(e.max() - 1.0) < 1e-12
    # A single driver can be .run() multiple times.
    pattern2 = driver.run()
    np.testing.assert_array_equal(
        np.asarray(pattern.theta_deg),
        np.asarray(pattern2.theta_deg),
    )
