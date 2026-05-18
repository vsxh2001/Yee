"""Tests for `yee.fem.solve_cavity`.

Mirrors the Rust-side `fem-eig-001` validation gate
(`crates/yee-validation/tests/fem_eig_001_rectangular_cavity.rs`) on the
Python binding: build a WR-90-based rectangular cavity, solve for the
ten lowest TE/TM modes, and assert the lowest mode lands within ±0.3 %
of the Pozar §6.3 analytic TE_{101} frequency for the
`(a, b, d) = (22.86, 10.16, 30) mm` geometry.

The cavity dimensions in this test are *not* the spec's literal
`(a, b, d) = (22.86, 10.16, 20) mm` (which would give 9.660 GHz) — Track
QQQQQQ's fem-eig-001 driver landed with `d = 30 mm`, for which Pozar
eq. 6.42 gives `TE_{101} ≈ 8.244 GHz`. We track the as-shipped Rust gate
exactly so the Python binding is a thin wrapper over the same solve.
"""

import pytest

import yee.fem  # noqa: F401 — submodule import sanity


def test_solve_cavity_te101_within_0_3_percent() -> None:
    """fem-eig-001: WR-90 cavity TE_{101} ≈ 8.244 GHz analytic.

    Pozar eq. 6.42 with `(a, b, d) = (22.86, 10.16, 30) mm`:
    `f_TE101 = (c / (2π)) · sqrt((π/a)² + (π/d)²) ≈ 8.2439 GHz`.
    """
    freqs, _modes = yee.fem.solve_cavity(0.02286, 0.01016, 0.030, 12, 9, 15)
    assert len(freqs) == 10
    assert freqs == sorted(freqs)
    expected_te101 = 8.2439e9
    rel_err = abs(freqs[0] - expected_te101) / expected_te101
    assert rel_err <= 0.003, f"TE_101 error {rel_err:.4%} exceeds +/-0.3% gate"


def test_solve_cavity_mode_count() -> None:
    """`num_eigs` controls both the frequency-list length and the
    mode-coefficient column count."""
    freqs, modes = yee.fem.solve_cavity(0.02286, 0.01016, 0.030, 8, 6, 10, num_eigs=5)
    assert len(freqs) == 5
    assert modes.shape[1] == 5


def test_solve_cavity_invalid_dims_raises() -> None:
    """Zero / non-positive extents must surface as a Python exception
    (mapped from `yee_mesh::Error::Invalid` through `yee_mesh_to_py`)."""
    with pytest.raises(Exception):
        yee.fem.solve_cavity(0.0, 0.01016, 0.030, 8, 6, 10)
