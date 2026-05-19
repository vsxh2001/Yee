"""Tests for `yee.fem.solve_cavity_dispersive`.

Mirrors the Rust-side Phase 4.fem.eig.1 D5 integration tests
(`crates/yee-fem/tests/dispersive_newton.rs`) on the Python binding:

* ``test_free_space_air_matches_solve_cavity`` — drive
  :func:`yee.fem.solve_cavity_dispersive` with a single bulk material
  ``{"tag": 0, "eps_inf": 1.0, "mu_r": 1.0, "poles": []}`` (free-space
  air) on the WR-90-based rectangular cavity. The Newton tracker must
  converge to the analytic TE_{101} resonance ≈ 8.244 GHz within ±0.3 %
  on Re(f) and |Im(f)| < 1 MHz.
* ``test_lossy_drude_returns_complex_frequency`` — same geometry, but
  the bulk material is a single-pole Drude oscillator. The converged
  frequency must be complex (Im ≠ 0), with ``converged is True`` and
  ``iterations <= 8``.

The cavity dimensions in this test mirror the D5 fixture's WR-90 geometry
(a, b, d) = (22.86, 10.16, 30) mm. Pozar eq. 6.42 gives
``f_TE101 ≈ 8.244 GHz`` for that geometry.

For ``cavity_uniform``, every tet is tagged with ``BULK_TAG = 0`` per
the `yee_fem::dispersive::BULK_TAG` convention; the standard input is
therefore a single-element materials list.
"""

import math

import pytest

import yee.fem  # noqa: F401 — submodule import sanity


# WR-90 cavity dimensions used by `crates/yee-fem/tests/dispersive_newton.rs`.
CAVITY_A_M = 0.022_86
CAVITY_B_M = 0.010_16
CAVITY_D_M = 0.030

# Mesh density — matches the Rust D5 fixture so wall-time stays bounded.
NX, NY, NZ = 8, 6, 10

# Speed of light in vacuum (m/s). Matches `yee_core::units::C0` exactly.
C0 = 299_792_458.0


def _f_te101_hz() -> float:
    """Analytic Pozar §6.3 TE_{101} resonance for the air-filled cavity."""
    return 0.5 * C0 * math.sqrt((1.0 / CAVITY_A_M) ** 2 + (1.0 / CAVITY_D_M) ** 2)


def test_free_space_air_matches_solve_cavity() -> None:
    """Air-filled WR-90 cavity — converged complex frequency matches the
    analytic Pozar §6.3 TE_{101} within ±0.3 % on Re(f) and |Im(f)| < 1 MHz.
    """
    f_te101 = _f_te101_hz()

    # Warm-start at 90 % of the analytic resonance — same off-target
    # warm-start as the Rust D5 free-space gate test.
    warm_start_hz = 0.9 * f_te101

    # Free-space material database — air at the bulk tag `0`.
    materials = [
        {
            "tag": 0,
            "eps_inf": 1.0,
            "mu_r": 1.0,
            "poles": [],
        },
    ]

    result = yee.fem.solve_cavity_dispersive(
        CAVITY_A_M,
        CAVITY_B_M,
        CAVITY_D_M,
        NX,
        NY,
        NZ,
        materials,
        warm_start_hz,
        max_iter=8,
        tol=1e-6,
    )

    assert set(result.keys()) >= {
        "frequency_hz",
        "k_complex",
        "iterations",
        "converged",
    }
    assert result["converged"] is True
    assert result["iterations"] <= 8

    freq = result["frequency_hz"]
    assert isinstance(freq, complex)
    rel_err = abs(freq.real - f_te101) / f_te101
    assert rel_err <= 3e-3, (
        f"free-space converged Re(f) = {freq.real:.6e} Hz vs analytic "
        f"{f_te101:.6e} Hz — rel err {rel_err:.4%} exceeds ±0.3 %"
    )

    # Air is non-dispersive — |Im(f)| should be vanishing.
    assert abs(freq.imag) < 1.0e6, (
        f"free-space converged Im(f) = {freq.imag:.6e} Hz should be |Im| < 1 MHz"
    )

    # Composed k = ω · √(μ₀ε₀ε) must be real-positive for air.
    k = result["k_complex"]
    assert isinstance(k, complex)
    assert k.real > 0.0
    assert abs(k.imag) < 1e-6 * abs(k.real)


def test_lossy_drude_returns_complex_frequency() -> None:
    """Same WR-90 cavity, single-pole Drude bulk filler.

    Drude parameters: ``ε_∞ = 3.78``, ``ω_p = 2π · 0.4 GHz``,
    ``γ = 2π · 2 GHz``. The plasma frequency is far below the warm-start
    air resonance (~8.24 GHz), so the Newton tracker stays inside the
    monotone-convergence basin (spec §11 risk #6). The converged ω
    picks up a non-zero imaginary part from the Drude loss term.
    """
    f_te101 = _f_te101_hz()
    warm_start_hz = f_te101  # warm-start at the analytic air resonance

    materials = [
        {
            "tag": 0,
            "eps_inf": 3.78,
            "mu_r": 1.0,
            "poles": [
                {
                    "kind": "drude",
                    "omega_p": 2.0 * math.pi * 0.4e9,
                    "gamma": 2.0 * math.pi * 2.0e9,
                },
            ],
        },
    ]

    result = yee.fem.solve_cavity_dispersive(
        CAVITY_A_M,
        CAVITY_B_M,
        CAVITY_D_M,
        NX,
        NY,
        NZ,
        materials,
        warm_start_hz,
        max_iter=8,
        tol=1e-6,
    )

    assert result["converged"] is True, (
        f"Drude Newton did not converge: result = {result}"
    )
    assert result["iterations"] <= 8

    freq = result["frequency_hz"]
    assert isinstance(freq, complex)
    # Real part must be positive and below the warm-start (ε_inf > 1
    # lowers the resonance from the air case).
    assert freq.real > 0.0
    assert freq.real < f_te101, (
        f"Drude converged Re(f) = {freq.real:.4e} Hz should sit below the "
        f"air-resonance warm-start {f_te101:.4e} Hz (ε_∞ > 1 lowers ω)"
    )
    # Im(f) must be non-zero — the Drude pole introduces loss.
    assert abs(freq.imag) > 0.0, (
        f"Drude converged Im(f) = {freq.imag:.4e} Hz should be non-zero"
    )

    # Composed k must also have a non-zero imaginary part.
    k = result["k_complex"]
    assert isinstance(k, complex)
    assert abs(k.imag) > 0.0


def test_invalid_pole_kind_raises() -> None:
    """An unknown pole kind must surface as a ValueError, not a panic."""
    materials = [
        {
            "tag": 0,
            "eps_inf": 1.0,
            "mu_r": 1.0,
            "poles": [{"kind": "definitely-not-a-pole-kind", "omega_p": 1.0, "gamma": 1.0}],
        },
    ]
    with pytest.raises(ValueError):
        yee.fem.solve_cavity_dispersive(
            CAVITY_A_M,
            CAVITY_B_M,
            CAVITY_D_M,
            NX,
            NY,
            NZ,
            materials,
            _f_te101_hz(),
        )
