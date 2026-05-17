"""Tests for `yee.eigensolver.TriMesh2D` / `yee.eigensolver.NumericalCrossSection`.

Mirrors the Rust-side WR-90 TE10 validation gate
(`crates/yee-mom/tests/eigensolver_wr90.rs`) on the Python binding: build
a 6×6-quad-grid mesh of the WR-90 cross-section, solve at 10 GHz, and
assert the cached `β` agrees with the analytic value within 1 %.
"""

import math

import pytest

import yee
from yee.eigensolver import NumericalCrossSection, TriMesh2D


# -- helpers --------------------------------------------------------------


def _rectangular_mesh(a: float, b: float, nx: int, ny: int) -> TriMesh2D:
    """Structured `nx × ny` quad-grid mesh of `[0, a] × [0, b]`.

    Each quad is split along the `(low-x, low-y) → (high-x, high-y)`
    diagonal into two CCW triangles. All triangles share material tag 0.
    Mirrors `rectangular_mesh` in `crates/yee-mom/tests/eigensolver_wr90.rs`.
    """
    vertices = [
        (a * i / nx, b * j / ny)
        for j in range(ny + 1)
        for i in range(nx + 1)
    ]

    def idx(i: int, j: int) -> int:
        return j * (nx + 1) + i

    triangles = []
    for j in range(ny):
        for i in range(nx):
            v00 = idx(i, j)
            v10 = idx(i + 1, j)
            v11 = idx(i + 1, j + 1)
            v01 = idx(i, j + 1)
            triangles.append((v00, v10, v11))
            triangles.append((v00, v11, v01))
    return TriMesh2D(vertices, triangles)


def _analytic_te10_beta(a: float, freq_hz: float, eps_r: float = 1.0) -> float:
    """Closed-form TE10 phase constant `β = sqrt(k² − (π/a)²)`."""
    c0 = 299_792_458.0
    k = 2.0 * math.pi * freq_hz * math.sqrt(eps_r) / c0
    kc = math.pi / a
    if k <= kc:
        return float("nan")
    return math.sqrt(k * k - kc * kc)


# -- tests ----------------------------------------------------------------


def test_trimesh2d_construction():
    """Build a 4-vertex 2-triangle unit-square mesh; sanity-check accessors."""
    vertices = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    triangles = [(0, 1, 2), (0, 2, 3)]
    mesh = TriMesh2D(vertices, triangles)
    assert mesh.n_verts() == 4
    assert mesh.n_tris() == 2
    assert mesh.area(0) > 0.0
    assert mesh.area(1) > 0.0
    # Both triangles partition a unit square ⇒ total area = 1.
    assert mesh.area(0) + mesh.area(1) == pytest.approx(1.0, abs=1e-12)
    cx, cy = mesh.centroid(0)
    # Centroid of (0,0), (1,0), (1,1) is (2/3, 1/3).
    assert cx == pytest.approx(2.0 / 3.0, abs=1e-12)
    assert cy == pytest.approx(1.0 / 3.0, abs=1e-12)


def test_trimesh2d_with_explicit_material_tags():
    """Explicit `vertex_material` / `triangle_material` round-trip cleanly."""
    vertices = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    triangles = [(0, 1, 2), (0, 2, 3)]
    mesh = TriMesh2D(
        vertices,
        triangles,
        vertex_material=[1, 2, 3, 4],
        triangle_material=[10, 20],
    )
    assert mesh.n_verts() == 4
    assert mesh.n_tris() == 2


def test_trimesh2d_rejects_clockwise():
    """CW winding ⇒ negative signed area ⇒ Rust returns Error::Invalid ⇒ ValueError."""
    vertices = [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0)]
    triangles = [(0, 1, 2)]
    with pytest.raises(ValueError):
        TriMesh2D(vertices, triangles)


def test_numerical_cross_section_pre_solve_caches_none():
    """`beta` / `z_w` start `None` before any solve."""
    vertices = [(0.0, 0.0), (1.0e-3, 0.0), (1.0e-3, 1.0e-3), (0.0, 1.0e-3)]
    triangles = [(0, 1, 2), (0, 2, 3)]
    mesh = TriMesh2D(vertices, triangles)
    nc = NumericalCrossSection(mesh, {0: complex(1, 0)}, {0: complex(1, 0)})
    assert nc.beta is None
    assert nc.z_w is None


def test_numerical_cross_section_solve_wr90():
    """WR-90 TE10 at 10 GHz: numerical β must be within 1 % of analytic ≈ 158.24 rad/m."""
    a = 22.86e-3
    b = 10.16e-3
    freq_hz = 10.0e9

    mesh = _rectangular_mesh(a, b, 6, 6)
    # 6 × 6 quads = 72 triangles, 49 vertices.
    assert mesh.n_tris() == 72
    assert mesh.n_verts() == 49

    nc = NumericalCrossSection(
        mesh,
        eps_r={0: complex(1, 0)},
        mu_r={0: complex(1, 0)},
    )
    nc.solve(freq_hz)

    beta = nc.beta
    assert beta is not None, "β should be cached after a successful solve"
    # WR-90 is lossless and air-filled ⇒ β is real (Im ≈ 0 to numerical noise).
    assert math.isfinite(beta.real)
    assert abs(beta.imag) < 1e-6 * abs(beta.real)

    analytic = _analytic_te10_beta(a, freq_hz)
    rel_err = abs(beta.real - analytic) / analytic
    assert rel_err < 0.01, (
        f"WR-90 TE10 β: numerical {beta.real:.6f} rad/m vs "
        f"analytic {analytic:.6f} rad/m (rel err {rel_err:.4f}); want < 1 %"
    )

    zw = nc.z_w
    assert zw is not None, "Z_w should be cached after a successful solve"
    # η₀ ≈ 376.73 Ω; TE10 Z_w at 10 GHz on WR-90 ≈ 500 Ω.
    assert 100.0 < abs(zw) < 1000.0


def test_numerical_cross_section_solve_real_only_materials():
    """Real-valued material specs (`{0: 1.0}`) should be accepted as `1 + 0j`."""
    a = 22.86e-3
    b = 10.16e-3
    mesh = _rectangular_mesh(a, b, 6, 6)
    nc = NumericalCrossSection(mesh, {0: 1.0}, {0: 1.0})
    nc.solve(10.0e9)
    assert nc.beta is not None


def test_import_from_yee_eigensolver():
    """`from yee.eigensolver import NumericalCrossSection` must work
    (the Track LLLL sys.modules registration covers this submodule too)."""
    from yee.eigensolver import NumericalCrossSection as Nc2
    from yee.eigensolver import TriMesh2D as Tm2
    assert Nc2 is NumericalCrossSection
    assert Tm2 is TriMesh2D
    # And `import yee.eigensolver` (vs attribute access) should also resolve:
    import yee.eigensolver as eig
    assert eig is yee.eigensolver
