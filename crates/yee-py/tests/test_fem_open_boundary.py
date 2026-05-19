"""Tests for `yee.fem.solve_open_cavity` (Phase 4.fem.eig.2 step E6).

Mirrors the Rust-side ``crates/yee-fem/tests/open_boundary_sweep.rs``
gates on the Python binding — runs a short frequency sweep against a
small WR-90-style stub, asserts shape + passivity + ABC monotonicity:

* ``test_open_cavity_sweep_returns_correct_shape`` — a 10-point sweep
  on a single TE_{10} wave-port + ABC opposite face produces a
  ``(10, 1, 1)`` complex ndarray.
* ``test_all_pec_except_one_port_returns_full_reflection`` — when every
  exterior face other than the +z wave-port is PEC, the swept
  ``|S_11|`` saturates near 1.0 (no absorber → full reflection).
* ``test_abc_smaller_than_pec`` — the same geometry once with the +z
  face PEC (no driver) is invalid; the comparison instead runs the
  driver face fixed as a wave-port and toggles the *opposite* face
  between PEC (fully reflective stub) and ABC (absorbing stub),
  asserting ``|S_11_abc| ≤ |S_11_pec|`` pointwise across the band.

The mesh is intentionally coarse (``nx = 3, ny = 2, nz = 4`` per
``crates/yee-fem/tests/open_boundary_sweep.rs``'s ``wr90_stub_mesh``)
so the per-frequency assemble + sparse-LU stays well below 30 s wall
time.
"""

from __future__ import annotations

import math

import numpy as np

import yee.fem  # noqa: F401 — submodule import sanity


# WR-90 stub dimensions matching ``crates/yee-fem/tests/open_boundary_sweep.rs``.
WR90_A = 0.022_86
WR90_B = 0.010_16
STUB_D = 0.030

NX, NY, NZ = 3, 2, 4

# Sweep band: 8-12 GHz (TE_{10} pass-band on WR-90), 10 uniform points.
F_MIN_HZ = 8.0e9
F_MAX_HZ = 12.0e9
N_FREQ = 10

C0 = 299_792_458.0


def _omegas_hz() -> list[float]:
    """Uniform sweep over [F_MIN_HZ, F_MAX_HZ] with N_FREQ points."""
    if N_FREQ == 1:
        return [F_MIN_HZ]
    return [
        F_MIN_HZ + (F_MAX_HZ - F_MIN_HZ) * (k / (N_FREQ - 1)) for k in range(N_FREQ)
    ]


def _te10_modal_e_t() -> tuple[float, float, float]:
    """Constant tangential modal E-field on the WR-90 wave-port.

    Phase 4.fem.eig.2 v0 binding samples ``modal_e_t`` at the face
    centroid only (constant across the face). We pick the orthonormalised
    TE_{10} amplitude at the broad-wall midpoint ``x = a/2``:
    ``e_mode = ŷ · sqrt(2/(a·b)) · sin(π · (a/2) / a) = ŷ · sqrt(2/(a·b))``.
    """
    norm = math.sqrt(2.0 / (WR90_A * WR90_B))
    return (0.0, norm, 0.0)


def _free_space_materials() -> list[dict]:
    """Single-bulk-tag air materials list."""
    return [{"tag": 0, "eps_inf": 1.0, "mu_r": 1.0, "poles": []}]


def test_open_cavity_sweep_returns_correct_shape() -> None:
    """A 10-point sweep × 1 port returns ``(10, 1, 1)`` complex ndarray."""
    omegas = _omegas_hz()
    e_t = _te10_modal_e_t()

    result = yee.fem.solve_open_cavity(
        WR90_A,
        WR90_B,
        STUB_D,
        NX,
        NY,
        NZ,
        _free_space_materials(),
        [
            {
                "axis": "z",
                "side": "high",
                "port_id": 0,
                "modal_e_t": e_t,
            }
        ],
        [
            {"axis": "z", "side": "low"},
        ],
        omegas,
    )

    assert isinstance(result, np.ndarray), "expected numpy.ndarray return"
    assert result.shape == (N_FREQ, 1, 1), (
        f"expected (n_omegas, n_ports, n_ports) = ({N_FREQ}, 1, 1), got {result.shape}"
    )
    assert np.iscomplexobj(result), "expected complex dtype"
    # No NaNs / Infs.
    assert np.isfinite(result).all(), "every S-parameter entry must be finite"


def test_all_pec_except_one_port_returns_full_reflection() -> None:
    """All exterior faces PEC except a single +z wave-port → ``|S_11| ≈ 1``.

    With no ABC face on the opposite end the stub is a fully-reflective
    closed cavity (modulo the wave-port radiating one mode out). The
    coarse mesh + face-centroid quadrature places ``|S_11|`` close to
    but not exactly 1.0; the gate is ``|S_11| ≥ 0.5`` per the published
    smoke-test floor in ``crates/yee-fem/tests/open_boundary_sweep.rs``
    (see ``s11_magnitude_bounded``).
    """
    omegas = _omegas_hz()
    e_t = _te10_modal_e_t()

    result = yee.fem.solve_open_cavity(
        WR90_A,
        WR90_B,
        STUB_D,
        NX,
        NY,
        NZ,
        _free_space_materials(),
        [
            {
                "axis": "z",
                "side": "high",
                "port_id": 0,
                "modal_e_t": e_t,
            }
        ],
        [],  # no ABC — every other face defaults to PEC
        omegas,
    )

    assert result.shape == (N_FREQ, 1, 1)
    mags = np.abs(result[:, 0, 0])
    # PEC stub is highly reflective. The exact saturation depends on
    # the coarse-mesh discretisation; the v0 walking-skeleton bound is
    # |S_11| >= 0.5 (Phase 4.fem.eig.2 E4 smoke test).
    assert (mags >= 0.5).all(), (
        f"all-PEC stub should be highly reflective; got |S_11| = {mags}"
    )
    # Passivity: |S_11| <= 1 + small numerical margin.
    assert (mags <= 1.5).all(), (
        f"passive structure cannot amplify by 50%+; got |S_11| = {mags}"
    )


def test_abc_smaller_than_pec() -> None:
    """ABC-terminated stub absorbs at least as much as the PEC-closed stub.

    Same geometry, same +z TE_{10} wave-port driver, only the back wall
    (-z) toggles between PEC (closed reflective stub) and ABC (absorbing
    terminator). Per passivity the ABC variant cannot produce a higher
    |S_11| than the PEC variant at any swept frequency.
    """
    omegas = _omegas_hz()
    e_t = _te10_modal_e_t()
    materials = _free_space_materials()
    port_faces = [
        {
            "axis": "z",
            "side": "high",
            "port_id": 0,
            "modal_e_t": e_t,
        }
    ]

    result_pec = yee.fem.solve_open_cavity(
        WR90_A,
        WR90_B,
        STUB_D,
        NX,
        NY,
        NZ,
        materials,
        port_faces,
        [],  # back wall defaults to PEC
        omegas,
    )
    result_abc = yee.fem.solve_open_cavity(
        WR90_A,
        WR90_B,
        STUB_D,
        NX,
        NY,
        NZ,
        materials,
        port_faces,
        [{"axis": "z", "side": "low"}],
        omegas,
    )

    assert result_pec.shape == result_abc.shape == (N_FREQ, 1, 1)
    s11_pec = np.abs(result_pec[:, 0, 0])
    s11_abc = np.abs(result_abc[:, 0, 0])
    # Pointwise passivity comparison. A small numerical tolerance
    # absorbs face-centroid-quadrature round-off; the band-averaged
    # difference must remain non-positive (ABC <= PEC).
    tol = 1e-6
    diff = s11_abc - s11_pec
    assert (diff <= tol).all(), (
        f"ABC variant must not exceed PEC variant in |S_11|; "
        f"got diff = {diff}, s11_abc = {s11_abc}, s11_pec = {s11_pec}"
    )
