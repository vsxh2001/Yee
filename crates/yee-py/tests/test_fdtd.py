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


# ---------------------------------------------------------------------------
# Phase 2.fdtd.py.0 — fdtd-202 lossy-cavity Q-factor
# ---------------------------------------------------------------------------


def test_cavity_q_default_passes():
    """fdtd-202: default σ₀=2.96e-3 S/m gives Q≈20, rel_err<5%."""
    from yee import CavityQResult, run_cavity_q

    result = run_cavity_q()
    assert isinstance(result, CavityQResult)
    assert result.passed, (
        f"fdtd-202 gate failed: q_measured={result.q_measured:.4f} "
        f"q_analytic={result.q_analytic:.4f} rel_err={result.rel_err:.4f}"
    )
    assert abs(result.q_analytic - 20.0) < 1.0, (
        f"analytic Q should be ~20, got {result.q_analytic}"
    )


def test_cavity_q_probe_array_shape():
    """probe_array() returns numpy float64 array of length n_ring."""
    import numpy as np

    from yee import run_cavity_q

    result = run_cavity_q(n_ring=100)
    arr = result.probe_array()
    assert isinstance(arr, np.ndarray)
    assert arr.shape == (100,)
    assert arr.dtype == np.float64


def test_cavity_q_repr_smoke():
    """__repr__ contains CavityQResult and q_measured."""
    from yee import run_cavity_q

    r = repr(run_cavity_q())
    assert "CavityQResult" in r
    assert "q_measured" in r


# ---------------------------------------------------------------------------
# Phase 2.fdtd.py.1 — fdtd-201 cavity resonance frequency
# ---------------------------------------------------------------------------


def test_cavity_resonance_passes_gate():
    """fdtd-201: DFT scan extracts TE₁₀₁ frequency within ±2.5% of analytic."""
    from yee import CavityResonanceResult, run_cavity_resonance

    result = run_cavity_resonance(n_steps=30_000)
    assert isinstance(result, CavityResonanceResult)
    assert result.passed, (
        f"fdtd-201 gate failed: f_extracted={result.f_extracted_hz:.6e} Hz, "
        f"f_analytic={result.f_analytic_hz:.6e} Hz, rel_err={result.rel_err:.4e}"
    )
    assert result.rel_err < 0.025, (
        f"rel_err={result.rel_err:.4e} exceeds 2.5% tolerance"
    )


def test_cavity_resonance_fields_and_probe_array():
    """Field types, f_analytic_hz > 0, probe_array shape and dtype."""
    import numpy as np

    from yee import run_cavity_resonance

    result = run_cavity_resonance(n_steps=15_000)
    assert result.f_analytic_hz > 0.0, (
        f"f_analytic_hz should be positive, got {result.f_analytic_hz}"
    )
    assert result.f_extracted_hz > 0.0, (
        f"f_extracted_hz should be positive, got {result.f_extracted_hz}"
    )
    arr = result.probe_array()
    assert isinstance(arr, np.ndarray)
    assert arr.shape == (15_000,)
    assert arr.dtype == np.float64


def test_cavity_resonance_repr_smoke():
    """__repr__ contains 'CavityResonanceResult' and 'f_extracted_hz'."""
    from yee import run_cavity_resonance

    r = repr(run_cavity_resonance(n_steps=100))
    assert "CavityResonanceResult" in r
    assert "f_extracted_hz" in r


# ---------------------------------------------------------------------------
# Phase 2.fdtd.py.2 — fdtd-203 short-dipole radiation-pattern gate
# ---------------------------------------------------------------------------


def test_run_dipole_pattern_smoke():
    """Smoke test: 5-step run verifies API plumbing, not physics."""
    from yee import DipolePatternResult, run_dipole_pattern

    r = run_dipole_pattern(n_steps=5)
    assert isinstance(r, DipolePatternResult)
    assert hasattr(r, "passed")
    assert hasattr(r, "e_theta_0")
    assert hasattr(r, "e_theta_45")
    assert hasattr(r, "e_theta_90")
    assert hasattr(r, "e_theta_135")
    assert hasattr(r, "e_theta_180")
    # Array accessors
    theta = r.theta_deg_array()
    e = r.e_theta_array()
    assert theta.shape == (37,)
    assert e.shape == (37,)
    assert theta[0] == 0.0
    assert theta[-1] == 180.0


def test_run_dipole_pattern_repr_smoke():
    """__repr__ contains DipolePatternResult and e_theta_0."""
    from yee import run_dipole_pattern

    r = repr(run_dipole_pattern(n_steps=5))
    assert "DipolePatternResult" in r
    assert "e_theta_0" in r


# ---------------------------------------------------------------------------
# Phase 2.fdtd.py.3 — fdtd-204 TF/SF Fresnel-transmission gate
# ---------------------------------------------------------------------------


def test_run_fresnel_tfsf_smoke():
    """Smoke test: 5-step run verifies API plumbing, not physics."""
    from yee import FresnelTfsfResult, run_fresnel_tfsf

    r = run_fresnel_tfsf(n_steps=5)
    assert isinstance(r, FresnelTfsfResult)
    assert hasattr(r, "t_measured")
    assert hasattr(r, "t_analytic")
    assert hasattr(r, "rel_err")
    assert hasattr(r, "passed")
    assert r.t_analytic > 0.0
    assert r.t_analytic < 1.0


def test_run_fresnel_tfsf_repr_smoke():
    """__repr__ contains FresnelTfsfResult and t_measured."""
    from yee import run_fresnel_tfsf

    r = repr(run_fresnel_tfsf(n_steps=5))
    assert "FresnelTfsfResult" in r
    assert "t_measured" in r


# =============================================================================
# Phase 2.fdtd.py.4 — cpml-001 / ntff-001 / dispersive-001 Python drivers
# =============================================================================


def test_run_cpml_reflection_returns_result():
    """run_cpml_reflection() returns a CpmlReflectionResult with expected fields."""
    import yee

    r = yee.run_cpml_reflection()
    assert isinstance(r, yee.CpmlReflectionResult)
    assert hasattr(r, "reduction_db")
    assert hasattr(r, "passed")
    assert isinstance(r.reduction_db, float)
    assert r.reduction_db > 0.0


def test_run_cpml_reflection_passes_gate():
    """cpml-001 gate: CPML reduces reflection by ≥ 30 dB vs PEC."""
    import yee

    r = yee.run_cpml_reflection()
    assert r.passed, f"cpml-001 gate failed: {r.reduction_db:.2f} dB < 30 dB"


def test_run_cpml_reflection_repr_smoke():
    """__repr__ contains 'CpmlReflectionResult' and 'reduction_db'."""
    import yee

    r = repr(yee.run_cpml_reflection())
    assert "CpmlReflectionResult" in r
    assert "reduction_db" in r


def test_run_ntff_broadside_returns_result():
    """run_ntff_broadside() returns an NtffResult with expected fields."""
    import yee

    r = yee.run_ntff_broadside()
    assert isinstance(r, yee.NtffResult)
    assert hasattr(r, "ratio_db")
    assert hasattr(r, "passed")
    assert isinstance(r.ratio_db, float)
    assert r.ratio_db > 0.0


def test_run_ntff_broadside_passes_gate():
    """ntff-001 gate: NTFF broadside/endfire ratio ≥ 20 dB."""
    import yee

    r = yee.run_ntff_broadside()
    assert r.passed, f"ntff-001 gate failed: {r.ratio_db:.2f} dB < 20 dB"


def test_run_ntff_broadside_repr_smoke():
    """__repr__ contains 'NtffResult' and 'ratio_db'."""
    import yee

    r = repr(yee.run_ntff_broadside())
    assert "NtffResult" in r
    assert "ratio_db" in r


def test_run_dispersive_drude_returns_result():
    """run_dispersive_drude() returns a DispersiveDrudeResult with expected fields."""
    import yee

    r = yee.run_dispersive_drude()
    assert isinstance(r, yee.DispersiveDrudeResult)
    assert hasattr(r, "gamma_measured")
    assert hasattr(r, "gamma_analytic")
    assert hasattr(r, "rel_err")
    assert hasattr(r, "passed")
    assert r.gamma_analytic > 0.0
    assert r.rel_err >= 0.0


def test_run_dispersive_drude_passes_gate():
    """dispersive-001 gate: Drude-slab Fresnel reflection rel_err ≤ 20%."""
    import yee

    r = yee.run_dispersive_drude()
    assert r.passed, (
        f"dispersive-001 gate failed: rel_err={r.rel_err:.2%} > 20% "
        f"(measured={r.gamma_measured:.4f}, analytic={r.gamma_analytic:.4f})"
    )


def test_run_dispersive_drude_repr_smoke():
    """__repr__ contains 'DispersiveDrudeResult' and 'gamma_measured'."""
    import yee

    r = repr(yee.run_dispersive_drude())
    assert "DispersiveDrudeResult" in r
    assert "gamma_measured" in r
