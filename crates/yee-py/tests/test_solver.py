"""Tests for `yee.PlanarMoM`.

Note: the underlying `yee_mom::PlanarMoM::run` is a Phase 0 stub that returns
`Error::Unimplemented` for every input, regardless of mesh tags. These tests
therefore verify the binding pipeline (constructor + run + error mapping)
against current upstream behavior. Once the Phase 1.0 mom-dipole physics work
lands and the solver returns real `SParameters` for valid port topologies,
test 1 should be tightened to assert on the returned data (n_ports == 1,
shape `(F, 1, 1)`, dtype `complex128`). Test 2 should then be updated to
match the actual error variant raised for port-less meshes.
"""

import numpy as np
import pytest

import yee


def test_planar_mom_runs_and_returns_sparams():
    """Wire the full pipeline: TriMesh + FreqRange -> PlanarMoM.run.

    With the Phase 0 stub solver, this raises RuntimeError("unimplemented: ...").
    The assertion captures that current contract; tighten when the solver lands.
    """
    v = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [0.1, 0.1, 0.0],
            [0.0, 0.1, 0.0],
        ],
        dtype=np.float64,
    )
    t = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.uint32)
    g = np.array([1, 2], dtype=np.uint32)
    mesh = yee.TriMesh(v, t, g)
    freq = yee.FreqRange(1.0e9, 1.5e9, 3)
    solver = yee.PlanarMoM()

    try:
        s = solver.run(mesh, freq)
    except RuntimeError as exc:
        # Phase 0 contract: stub solver raises mapped Unimplemented.
        assert "unimplemented" in str(exc), f"unexpected message: {exc}"
        return

    # Once the solver is implemented, this is the success path the spec
    # describes (n_ports == 1 for a single-port planar geometry).
    assert s.n_ports == 1
    assert s.freq_hz.shape == (3,)
    assert s.data.shape == (3, 1, 1)
    assert s.data.dtype == np.complex128


def test_planar_mom_returns_numerical_error_without_port():
    """A mesh whose tags do not define a port must surface as a RuntimeError.

    Per spec the production solver raises `Error::Numerical` here, which the
    binding maps to `RuntimeError("numerical: ...")`. Until that lands the
    stub solver returns `Error::Unimplemented`, mapped to
    `RuntimeError("unimplemented: ...")`. Either is acceptable today; the
    test asserts the error is a RuntimeError whose message matches one of
    those two expected prefixes.
    """
    v = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [0.1, 0.1, 0.0],
            [0.0, 0.1, 0.0],
        ],
        dtype=np.float64,
    )
    t = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.uint32)
    g = np.array([0, 0], dtype=np.uint32)  # no port: identical tags
    mesh = yee.TriMesh(v, t, g)
    freq = yee.FreqRange(1.0e9, 1.5e9, 3)
    solver = yee.PlanarMoM()

    with pytest.raises(RuntimeError, match=r"(numerical|unimplemented)"):
        solver.run(mesh, freq)
