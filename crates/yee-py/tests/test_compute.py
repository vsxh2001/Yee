"""Tests for the `yee.compute` engine bindings (ADR-0178).

These gate the *binding* layer — shapes, dtypes, validation, determinism,
backend selection — not the physics: the engine itself is gated in Rust
(compute-001..010) against bit-exact references and analytic results.
"""

import numpy as np
import pytest

from yee.compute import FdtdSim


def make_driven_sim():
    sim = FdtdSim(16, 12, 14, 1e-3)
    dt = sim.dt
    sim.set_boundary("pec")
    sim.add_gaussian_source("ey", (4, 6, 7), 8.0 * dt, 3.0 * dt)
    sim.add_resistive_port((10, 6, 7), 50.0, 1.0, 5.0e9, 4.0e9, 6)
    sim.add_probe("ez", (12, 6, 7))
    return sim


def test_basic_run_shapes_and_propagation():
    sim = make_driven_sim()
    res = sim.run(50, backend="cpu")
    assert res.backend == "cpu"
    assert res.n_probes == 1

    probe = res.probe(0)
    assert probe.shape == (50,)
    assert probe.dtype == np.float64
    assert np.all(np.isfinite(probe))
    assert np.any(probe != 0.0), "probe never saw the drive"

    # Staggered field shapes.
    assert res.field("ex").shape == (16, 13, 15)
    assert res.field("ey").shape == (17, 12, 15)
    assert res.field("ez").shape == (17, 13, 14)
    assert res.field("hx").shape == (17, 12, 14)
    assert res.field("hy").shape == (16, 13, 14)
    assert res.field("hz").shape == (16, 12, 15)
    assert np.any(res.field("hx") != 0.0), "H never energized"


def test_runs_are_deterministic_and_repeatable():
    sim = make_driven_sim()
    a = sim.run(40, backend="cpu")
    b = sim.run(40, backend="cpu")
    assert np.array_equal(a.probe(0), b.probe(0))
    assert np.array_equal(a.field("ez"), b.field("ez"))


def test_materials_and_cpml():
    sim = FdtdSim(12, 12, 12, 1e-3)
    eps = np.ones((13, 13, 13), dtype=np.float64)
    eps[:, :, 3:6] = 4.3
    sim.set_eps_r_cells(eps)
    mask = np.zeros((12, 13, 13), dtype=bool)
    mask[:, :5, 6] = True
    sim.set_pec_mask("ex", mask)
    sim.set_boundary("cpml", npml=4)
    dt = sim.dt
    sim.add_gaussian_source("ez", (6, 6, 6), 8.0 * dt, 3.0 * dt)
    sim.add_probe("ez", (8, 6, 6))
    res = sim.run(60, backend="cpu")
    assert np.all(np.isfinite(res.probe(0)))
    assert np.any(res.probe(0) != 0.0)
    # Masked E_x cells are clamped to zero.
    ex = res.field("ex")
    assert np.all(ex[:, :5, 6] == 0.0)


def test_shape_validation_rejects_wrong_arrays():
    sim = FdtdSim(8, 8, 8, 1e-3)
    with pytest.raises(ValueError):
        sim.set_eps_r_cells(np.ones((8, 8, 8)))  # must be (9, 9, 9)
    with pytest.raises(ValueError):
        sim.set_pec_mask("ez", np.zeros((9, 9, 9), dtype=bool))  # must be (9, 9, 8)
    with pytest.raises(ValueError):
        sim.set_boundary("reflecting")
    with pytest.raises(ValueError):
        sim.add_probe("hx", (0, 0, 0))  # probes are E-only in the binding


def test_gpu_backend_matches_cpu_or_skips():
    sim = make_driven_sim()
    try:
        gpu = sim.run(50, backend="gpu")
    except RuntimeError as e:
        assert "adapter" in str(e).lower()
        pytest.skip("no wgpu adapter on this machine")
    cpu = sim.run(50, backend="cpu")
    assert gpu.backend == "gpu"
    # FP32 GPU vs FP64 CPU: loose elementwise agreement on the probe.
    scale = np.max(np.abs(cpu.probe(0)))
    assert scale > 0.0
    assert np.max(np.abs(gpu.probe(0) - cpu.probe(0))) < 1e-3 * scale


def test_auto_backend_always_runs():
    sim = make_driven_sim()
    res = sim.run(30, backend="auto")
    assert res.backend in ("cpu", "gpu")
    assert np.any(res.probe(0) != 0.0)
