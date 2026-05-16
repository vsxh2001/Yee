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
