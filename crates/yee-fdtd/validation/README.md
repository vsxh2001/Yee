# yee-fdtd — Validation

Phase 2 deliverable. Phase 0/1: no live cases. Below is the planned suite.

## Cases — Phase 2

| ID | Case | Reference | Tolerance |
|----|------|-----------|-----------|
| `fdtd-201` | Rectangular cavity TE/TM Q-factor | Analytical | ±0.5% |
| `fdtd-202` | Pyramidal horn antenna pattern | Measured / Balanis | ±1 dB main beam |
| `fdtd-203` | Dipole over dielectric half-space NTFF | Sommerfeld reference | analytic match |
| `fdtd-204` | Cross-validation vs openEMS | openEMS on identical grid | numerical-noise level |
| `fdtd-205` | Microstrip transient TDR | FFT(yee-mom Sxx) | ±2% |
| `fdtd-206` | Drude-metal plasmonic dipole | Maier / textbook | ±5% resonance |
| `fdtd-207` | Multi-pole Debye human-tissue benchmark | Gabriel database | ±5% absorption |

## Running

Will require GPU + CUDA toolkit + `yee-cuda` feature `cuda`.

```bash
cargo test -p yee-fdtd --release --features cuda
```

## Cross-tool validation

For every case where openEMS or gprMax can run the same geometry, we publish side-by-side results in `validation/results/` so a reader can verify our work without trusting our numbers.
