"""Tests for yee.nsga2_minimize (PyO3 wrapper around yee-surrogate::nsga2).

The ZDT1 benchmark has a known convex Pareto front
    f2 = 1 - sqrt(f1)
on f1 in [0, 1]. With a population of 100 over 100 generations NSGA-II
recovers a well-spread non-dominated set.
"""

import numpy as np

from yee import Nsga2Config, Nsga2Result, nsga2_minimize


def zdt1(x):
    d = len(x)
    f1 = float(x[0])
    g = 1.0 + 9.0 / (d - 1) * float(np.sum(x[1:]))
    f2 = g * (1.0 - np.sqrt(f1 / g))
    return [f1, f2]


def test_nsga2_zdt1():
    d = 30
    cfg = Nsga2Config(population_size=100, n_generations=100, seed=42)
    result = nsga2_minimize(zdt1, [(0.0, 1.0)] * d, 2, cfg)

    assert isinstance(result, Nsga2Result)
    assert result.population.shape == (100, d)
    assert result.objectives.shape == (100, 2)
    assert result.population.dtype == np.float64
    assert result.objectives.dtype == np.float64

    # At least 30 non-dominated solutions after 100 generations is a
    # generous lower bound (well-tuned runs usually approach the full N).
    pareto_indices = result.pareto_indices
    assert pareto_indices.dtype == np.int64
    assert len(pareto_indices) >= 30, (
        f"too few Pareto-optimal solutions: {len(pareto_indices)}"
    )

    # Every Pareto index is a valid row.
    assert pareto_indices.min() >= 0
    assert pareto_indices.max() < 100


def test_nsga2_zdt1_default_mutation_probability():
    d = 5
    # mutation_probability=None should resolve to 1/d inside the wrapper.
    cfg = Nsga2Config(population_size=20, n_generations=10, seed=7)
    result = nsga2_minimize(zdt1, [(0.0, 1.0)] * d, 2, cfg)
    assert result.population.shape == (20, d)
    assert result.objectives.shape == (20, 2)


def test_nsga2_config_repr_and_accessors():
    cfg = Nsga2Config(
        population_size=50,
        n_generations=30,
        crossover_eta=15.0,
        mutation_eta=25.0,
        mutation_probability=0.1,
        seed=11,
    )
    assert cfg.population_size == 50
    assert cfg.n_generations == 30
    assert cfg.crossover_eta == 15.0
    assert cfg.mutation_eta == 25.0
    assert cfg.mutation_probability == 0.1
    assert cfg.seed == 11
    r = repr(cfg)
    assert "Nsga2Config" in r
    assert "population_size=50" in r


def test_nsga2_default_config():
    d = 2
    # Two independent quadratic objectives — quick smoke test.
    def two_objs(x):
        return [(x[0] - 0.2) ** 2, (x[1] - 0.8) ** 2]

    result = nsga2_minimize(two_objs, [(0.0, 1.0)] * d, 2)
    # Default cfg: population_size=100, n_generations=100.
    assert result.population.shape == (100, 2)
    assert result.objectives.shape == (100, 2)


def test_nsga2_wrong_objective_length_raises():
    import pytest

    def bad_objective(x):
        return [float(x[0])]  # returns 1 value but we ask for 2

    with pytest.raises(ValueError, match="objective returned 1 values, expected 2"):
        nsga2_minimize(
            bad_objective,
            [(0.0, 1.0)] * 3,
            2,
            Nsga2Config(population_size=10, n_generations=2),
        )
