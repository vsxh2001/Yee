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


def _horizontal_slab_mesh(a: float, b: float, nx: int, ny: int) -> TriMesh2D:
    """`nx × ny` quad-grid mesh whose lower-`y` half carries material tag 1.

    Same structured triangulation as `_rectangular_mesh`, but each
    triangle is tagged by its quad's `y`-midpoint: tag ``1`` (dielectric)
    for ``y < b/2``, tag ``0`` (air) above. This is the dielectric-loaded
    WR-90 cross-section the Rust ``fr4_loaded_beta_matches_reference``
    gate uses (the §4 inhomogeneous closure at FR-4); a horizontal slab's
    dominant mode is genuinely hybrid, so the longitudinal `E_z` block is
    load-bearing.
    """
    vertices = [
        (a * i / nx, b * j / ny)
        for j in range(ny + 1)
        for i in range(nx + 1)
    ]

    def idx(i: int, j: int) -> int:
        return j * (nx + 1) + i

    triangles = []
    triangle_material = []
    for j in range(ny):
        y_mid = b * (j + 0.5) / ny
        tag = 1 if y_mid < b / 2.0 else 0
        for i in range(nx):
            v00 = idx(i, j)
            v10 = idx(i + 1, j)
            v11 = idx(i + 1, j + 1)
            v01 = idx(i, j + 1)
            triangles.append((v00, v10, v11))
            triangle_material.append(tag)
            triangles.append((v00, v11, v01))
            triangle_material.append(tag)
    return TriMesh2D(vertices, triangles, triangle_material=triangle_material)


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


def test_numerical_cross_section_solve_inhomogeneous_fr4():
    """Dielectric-loaded WR-90 (lower half ε_r = 4.4): β obeys the monotonic
    empty/full bracket, and the new `mode_profile_ez` getter is populated.

    The cross-section eigensolver computes β/Z_w on an INHOMOGENEOUS
    (dielectric-loaded) guide. Loading the lower-`y` half with FR-4
    (ε_r = 4.4, two material tags) slows the wave, so the loaded β must
    sit strictly inside the rigorous monotonic bracket

        β_air  <  β_loaded  <  β_fully-filled

    (β increases monotonically with dielectric fill fraction; the empty
    and fully-filled TE10 closed forms are the bracket because k_c = π/a
    is fixed by the PEC walls). This mirrors the Rust-side
    ``dod_v2_prime_loaded_beta_bracket_and_regression`` /
    ``fr4_loaded_beta_matches_reference`` gates in
    ``crates/yee-mom/tests/eigensolver_inhomogeneous.rs``. We assert the
    mesh-independent physics bracket rather than a pinned regression
    value so the test is robust to mesh density / formulation tweaks.
    """
    a = 22.86e-3
    b = 10.16e-3
    freq_hz = 10.0e9
    eps_fill = 4.4  # FR-4 substrate

    # 8×8-quad horizontal-slab mesh: lower-y half tag 1 (FR-4), upper tag 0
    # (air). Matches the Rust FR-4 gate geometry/density.
    mesh = _horizontal_slab_mesh(a, b, 8, 8)
    nc = NumericalCrossSection(
        mesh,
        eps_r={0: complex(1.0, 0.0), 1: complex(eps_fill, 0.0)},
        mu_r={0: complex(1.0, 0.0), 1: complex(1.0, 0.0)},
    )
    nc.solve(freq_hz)

    beta = nc.beta
    assert beta is not None, "β should be cached after a successful solve"
    # Lossless dielectric ⇒ β is real to numerical noise.
    assert math.isfinite(beta.real)
    assert abs(beta.imag) < 1e-6 * abs(beta.real)

    beta_air = _analytic_te10_beta(a, freq_hz, eps_r=1.0)
    beta_full = _analytic_te10_beta(a, freq_hz, eps_r=eps_fill)
    assert beta.real > beta_air, (
        f"loaded β {beta.real:.4f} must exceed air β {beta_air:.4f} "
        "(dielectric slows the wave)"
    )
    assert beta.real < beta_full, (
        f"loaded β {beta.real:.4f} must be below fully-filled β "
        f"{beta_full:.4f} (partial fill is bracketed by empty/full)"
    )

    # Z_w is finite and positive-real-dominated on the lossless guide.
    zw = nc.z_w
    assert zw is not None, "Z_w should be cached after a successful solve"
    assert math.isfinite(zw.real) and math.isfinite(zw.imag)
    assert zw.real > 0.0, "Z_w must be positive-real-dominated"
    assert abs(zw.imag) < 1e-6 * abs(zw.real), "lossless guide → Z_w ~ real"

    # The new longitudinal-E_z getter is populated post-solve, one complex
    # amplitude per global mesh vertex.
    ez = nc.mode_profile_ez
    assert ez is not None, "mode_profile_ez should be cached after solve"
    assert len(ez) == mesh.n_verts()
    assert all(isinstance(c, complex) for c in ez)
    assert all(math.isfinite(c.real) and math.isfinite(c.imag) for c in ez)


def test_mode_profile_ez_none_before_solve():
    """`mode_profile_ez` starts `None` before any solve (mirrors `beta`/`z_w`)."""
    a = 22.86e-3
    b = 10.16e-3
    mesh = _horizontal_slab_mesh(a, b, 4, 4)
    nc = NumericalCrossSection(
        mesh,
        {0: complex(1.0, 0.0), 1: complex(4.4, 0.0)},
        {0: complex(1.0, 0.0), 1: complex(1.0, 0.0)},
    )
    assert nc.mode_profile_ez is None


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
