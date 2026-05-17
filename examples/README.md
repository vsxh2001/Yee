# Yee Examples

End-to-end examples that exercise the full stack. Every example is **scripted** (no GUI-only artifacts) and **runs in CI** (Phase 1+).

Planned Phase 0/1 lineup:

| Example | Phase | What it demonstrates |
|---------|-------|----------------------|
| `microstrip-line` | 0 | Gmsh import → mesh count; smoke test for the meshing pipeline. |
| `half-wave-dipole` | 0 | Free-space dipole MoM solve on CPU, Z ≈ 73 + j42 Ω. |
| `patch-2g4` | 0 → 1 | 2.4 GHz rectangular patch on FR-4. Starts as a stub; grows into the full solver demo. |
| `hairpin-bpf-5pole` | 1 | Swanson 5-pole hairpin band-pass filter on RT/Duroid 6006. |
| `wilkinson-2ghz` | 1 | Wilkinson divider, ±0.5 dB vs Pozar reference. |
| `branch-line-hybrid` | 1 | 90° hybrid amplitude / phase balance. |
| `inset-patch-ro4003c` | 1 | Inset-fed patch on RO4003C; figure-for-figure paper match. |
| `cavity-q-factor` | 2 | FDTD rectangular cavity TE/TM modes vs analytical. |
| `pyramidal-horn` | 2 | FDTD horn radiation pattern vs measured. |

Each example lives in its own directory with:

```
examples/<name>/
├── README.md         # one-page description + reference + expected output
├── Cargo.toml        # standalone or via workspace example entry
├── src/main.rs       # the example
├── reference/        # baseline numbers / plots
└── plots/            # produced output (gitignored, regenerated)
```

## Running

```bash
cargo run --release --example half-wave-dipole
```

Phase 1+ examples will additionally have a Python sibling in `examples/<name>/python/` showing the Jupyter / PyO3 workflow.

## Python notebooks

Standalone Jupyter notebooks that exercise the `yee` Python bindings end-to-end live under [`examples/python/`](python/).

| Notebook | What it demonstrates |
|----------|----------------------|
| [`python/bo_monopole.ipynb`](python/bo_monopole.ipynb) | Bayesian optimization (`yee.bo_minimize`) on a synthetic monopole-length VSWR objective. Closed-form mock objective; converges to L ≈ λ/4 in 30 evaluations. |
| [`python/nsga2_pareto.ipynb`](python/nsga2_pareto.ipynb) | Multi-objective Pareto front via NSGA-II (`yee.nsga2_minimize`) on Schaffer N1; recovers the analytic front x ∈ [0, 2]. |
| [`python/al_dipole_surrogate.ipynb`](python/al_dipole_surrogate.ipynb) | Active learning (`yee.active_learn`) building a sin(x) GP surrogate; compares AL vs random sampling RMSE. |
| [`python/surrogate_workflow.ipynb`](python/surrogate_workflow.ipynb) | Full AL → GP → BO pipeline on a synthetic monopole. |
| [`python/touchstone_workflow.ipynb`](python/touchstone_workflow.ipynb) | Load `.s1p` via `yee.touchstone.read`; plot `S11` dB + Smith + `Z_in`. |

To run:

```bash
# one-time setup
uv venv .venv && source .venv/bin/activate
uv pip install maturin pytest numpy matplotlib jupyter
(cd crates/yee-py && maturin develop --release)

# launch the notebook
jupyter notebook examples/python/bo_monopole.ipynb
```
