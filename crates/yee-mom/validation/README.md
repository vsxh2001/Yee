# yee-mom — Validation

Every solver feature is held against a canonical published benchmark before it ships. No exceptions.

## Canonical references

- Pozar, *Microwave Engineering* (4th ed.) — closed-form microstrip / patch / Wilkinson
- Hong & Lancaster, *Microstrip Filters for RF / Microwave Applications* (2nd ed.) — coupled-line, hairpin BPFs
- Swanson & Hoefer, *Microwave Circuit Modeling Using Electromagnetic Field Simulation* — Sonnet-validated reference cases
- IEEE AP-S and MTT-S transactions for any paper-specific case

## Cases — Phase 0

| ID | Case | Tolerance |
|----|------|-----------|
| `mom-001` | Half-wave dipole, L=1m, a=5mm cylinder, delta-gap, NEC-4 reference Z ≈ 87 + j41 Ω | ±5% on Re(Z), ±10% on Im(Z) |
| `mom-002` | 50 Ω microstrip line Z₀ on FR-4 (h=1.6, εr=4.4) | ±3% vs TX-LINE / Hammerstad-Jensen |
| `mom-003` | 2.4 GHz rectangular patch on FR-4 (29.2×38.0 mm) | resonance ±2%; \|S11\| < −10 dB |

## Cases — Phase 1

| ID | Case | Tolerance |
|----|------|-----------|
| `mom-101` | Swanson 5-pole hairpin BPF (RT/Duroid 6006, εr=6.15, h=1.27, ~2 GHz) | ±1 dB to 4 GHz; resonances ±0.5% |
| `mom-102` | Parallel-coupled-line BPF (Hong & Lancaster Ch.5) | ±1 dB |
| `mom-103` | Wilkinson divider @ 2 GHz | ±0.5 dB vs Pozar |
| `mom-104` | Branch-line 90° hybrid | amplitude/phase balance |
| `mom-105` | Inset-fed patch on RO4003C | figure-for-figure match |
| `mom-106` | Cross-validation vs openEMS (microstrip + patch) | ±3% at resonance |

## Running

```bash
# Phase 0 (CPU)
cargo test -p yee-mom --release

# Phase 1 (GPU)
cargo test -p yee-mom --release --features cuda -- --include-ignored
```

Results land in `validation/results/` (gitignored) and are regenerated nightly in CI.

## Plot artifacts

Each validation case emits an S-parameter PNG via `plotters` for human review. CI publishes these to GitHub Pages so trend regressions are visible at a glance.
