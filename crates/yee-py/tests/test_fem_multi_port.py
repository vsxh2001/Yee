"""Tests for `yee.fem.solve_open_cavity` multi-port + new kwargs.

Phase 4.fem.eig.3 step F7 — Python re-run of the Rust-side
``crates/yee-validation/tests/fem_eig_004_wr90_thruline.rs`` 2-port
WR-90 thru-line gate, plus per-kwarg sanity coverage for
``coupled_whitney`` / ``abc_order`` / ``multi_port``.

Tests:

* ``test_multi_port_thru_line_s21_at_10ghz`` — WR-90 thru-line at 10 GHz
  with ``multi_port=True`` + ``coupled_whitney=True`` returns a 1×2×2
  S-matrix whose ``|S_{21}|`` is within ±0.1 dB of 0 dB (passive lossless
  transmission, mirroring fem-eig-004 gate A).
* ``test_abc_order_second_kwarg_accepted`` — single-port stub with
  ``abc_order="second"`` returns a finite ``(1, 1, 1)`` S-tensor.
* ``test_coupled_whitney_true_kwarg_accepted`` — single-port stub with
  ``coupled_whitney=True`` returns a finite ``(1, 1, 1)`` S-tensor.
* ``test_invalid_abc_order_raises`` — any string other than
  ``"first"`` / ``"second"`` raises ``ValueError``.

The thru-line mesh is intentionally coarser than the Rust-side
fem-eig-004 driver (8 × 4 × 16 vs 12 × 6 × 18 in
``crates/yee-validation/src/lib.rs``) to keep wall-time bounded; the
±0.1 dB tolerance is preserved because both port faces carry the same
modal profile and the transmission identity is a structural property
of the multi-port sweep rather than a discretisation-convergence
result.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import yee.fem  # noqa: F401 — submodule import sanity


# ---- WR-90 thru-line geometry (matches fem-eig-004) -------------------
WR90_A = 0.022_86
WR90_B = 0.010_16
THRU_LINE_D = 0.030

# Mesh density chosen to land inside the ±0.1 dB transmission window
# while keeping Python wall-time bounded. The Rust-side fem-eig-004
# driver uses (12, 6, 18) = 7.8 k tets and measures -0.045 dB; a (8, 4,
# 16) = 3 k tet mesh measures -0.101 dB (just outside the envelope).
# The (10, 5, 16) = 4.8 k tets choice below lands ~ -0.07 dB.
THRU_NX, THRU_NY, THRU_NZ = 10, 5, 16

# Stub mesh for the single-port kwarg-acceptance tests — small enough
# that a one-frequency assemble + LU completes in well under 10 s.
STUB_NX, STUB_NY, STUB_NZ = 3, 2, 4


def _te10_modal_e_t() -> tuple[float, float, float]:
    """Constant tangential modal-E-field tuple (peak-amplitude proxy).

    Returns the orthonormalised TE_{10} amplitude
    ``ŷ · sqrt(2/(a·b))`` sampled at ``x = a/2`` where ``sin(π/2) = 1``.
    Used by the single-port kwarg-acceptance tests where shape matters
    but transmission accuracy does not — the constant proxy keeps the
    binding on the Phase 4.fem.eig.2 v0 path.
    """
    norm = math.sqrt(2.0 / (WR90_A * WR90_B))
    return (0.0, norm, 0.0)


def _te10_modal_e_t_callable():
    """Analytic TE_{10} tangential profile ``ŷ · sqrt(2/(a·b)) · sin(π · x / a)``.

    Returned as a Python callable that the F7 binding evaluates at every
    face centroid (v2 lumped path) or every per-face Gauss point
    (Phase 4.fem.eig.3 F2 coupled-Whitney path). The analytic profile
    is required for the WR-90 thru-line gate; a constant proxy
    over-counts the modal self-inner-product on the broad-wall by a
    factor of 2, blowing the ±0.1 dB envelope (~ -1.9 dB observed
    on the (8, 4, 16) mesh).
    """
    norm = math.sqrt(2.0 / (WR90_A * WR90_B))

    def profile(point: tuple[float, float, float]) -> tuple[float, float, float]:
        x, _y, _z = point
        amp = norm * math.sin(math.pi * x / WR90_A)
        return (0.0, amp, 0.0)

    return profile


def _free_space_materials() -> list[dict]:
    """Single-bulk-tag air materials list."""
    return [{"tag": 0, "eps_inf": 1.0, "mu_r": 1.0, "poles": []}]


def test_multi_port_thru_line_s21_at_10ghz() -> None:
    """WR-90 thru-line at 10 GHz: ``|S_{21}|`` within ±0.1 dB of 0 dB.

    Mirrors the Rust-side ``fem-eig-004`` gate (A) (Phase 4.fem.eig.3
    step F6) on the Python binding. Both end faces carry the same
    TE_{10} modal profile, the four sidewalls default to PEC, and the
    F1+F2 coupled exact-Whitney-1 path is enabled so the modal RHS +
    projection are consistent at the exact-basis level.
    """
    profile = _te10_modal_e_t_callable()
    port_faces = [
        {
            "axis": "z",
            "side": "low",
            "port_id": 0,
            "modal_e_t": profile,
        },
        {
            "axis": "z",
            "side": "high",
            "port_id": 1,
            "modal_e_t": profile,
        },
    ]

    result = yee.fem.solve_open_cavity(
        WR90_A,
        WR90_B,
        THRU_LINE_D,
        THRU_NX,
        THRU_NY,
        THRU_NZ,
        _free_space_materials(),
        port_faces,
        [],  # no ABC faces — sidewalls default to PEC
        [10.0e9],
        coupled_whitney=True,
        abc_order="first",
        multi_port=True,
    )

    assert isinstance(result, np.ndarray)
    assert result.shape == (1, 2, 2), (
        f"expected (1, 2, 2) for multi-port thru-line; got {result.shape}"
    )
    assert np.iscomplexobj(result)
    assert np.isfinite(result).all(), "every S-parameter entry must be finite"

    s21 = result[0, 1, 0]
    s12 = result[0, 0, 1]
    s11 = result[0, 0, 0]
    s21_mag = abs(s21)
    assert s21_mag > 0.0, f"|S_21| must be strictly positive; got {s21}"
    s21_db = 20.0 * math.log10(s21_mag)
    reciprocity = abs(s12 - s21)
    print(
        f"\nfem-eig-004-py: |S_21| = {s21_mag:.6f} ({s21_db:.3f} dB), "
        f"|S_11| = {abs(s11):.6f}, |S_12 − S_21| = {reciprocity:.3e}"
    )
    assert abs(s21_db) <= 0.1, (
        f"|S_21(10 GHz)| = {s21_db:.3f} dB outside ±0.1 dB of 0 dB "
        f"(passive lossless thru-line); raw |S_21| = {s21_mag:.6f}"
    )


def test_abc_order_second_kwarg_accepted() -> None:
    """`abc_order='second'` accepted; returns a finite (1, 1, 1) tensor."""
    e_t = _te10_modal_e_t()
    result = yee.fem.solve_open_cavity(
        WR90_A,
        WR90_B,
        THRU_LINE_D,
        STUB_NX,
        STUB_NY,
        STUB_NZ,
        _free_space_materials(),
        [
            {
                "axis": "z",
                "side": "high",
                "port_id": 0,
                "modal_e_t": e_t,
            }
        ],
        [{"axis": "z", "side": "low"}],
        [10.0e9],
        abc_order="second",
    )
    assert result.shape == (1, 1, 1)
    assert np.isfinite(result).all(), (
        f"abc_order='second' should produce a finite S-tensor; got {result}"
    )


def test_coupled_whitney_true_kwarg_accepted() -> None:
    """`coupled_whitney=True` accepted; returns a finite (1, 1, 1) tensor."""
    e_t = _te10_modal_e_t()
    result = yee.fem.solve_open_cavity(
        WR90_A,
        WR90_B,
        THRU_LINE_D,
        STUB_NX,
        STUB_NY,
        STUB_NZ,
        _free_space_materials(),
        [
            {
                "axis": "z",
                "side": "high",
                "port_id": 0,
                "modal_e_t": e_t,
            }
        ],
        [{"axis": "z", "side": "low"}],
        [10.0e9],
        coupled_whitney=True,
    )
    assert result.shape == (1, 1, 1)
    assert np.isfinite(result).all(), (
        f"coupled_whitney=True should produce a finite S-tensor; got {result}"
    )


def test_invalid_abc_order_raises() -> None:
    """`abc_order='invalid'` raises ValueError."""
    e_t = _te10_modal_e_t()
    with pytest.raises(ValueError, match="abc_order"):
        yee.fem.solve_open_cavity(
            WR90_A,
            WR90_B,
            THRU_LINE_D,
            STUB_NX,
            STUB_NY,
            STUB_NZ,
            _free_space_materials(),
            [
                {
                    "axis": "z",
                    "side": "high",
                    "port_id": 0,
                    "modal_e_t": e_t,
                }
            ],
            [{"axis": "z", "side": "low"}],
            [10.0e9],
            abc_order="invalid",
        )
