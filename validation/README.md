# Yee — Workspace-level Validation

This directory is the **single source of truth** for end-to-end, cross-solver validation. Per-crate validation lives in `crates/<name>/validation/` and asserts crate-local correctness; this directory is for cases that exercise the whole pipeline (mesh → solve → I/O → comparison vs reference).

Rule, enforced by CI:

> **No solver feature ships without a published-benchmark validation case here and a CI run that regenerates results nightly.**

## Layout

```
validation/
├── README.md                 # this file
├── cases/
│   ├── phase-0/              # foundation: dipole, microstrip Z₀, basic patch
│   ├── phase-1/              # MoM beachhead: hairpin BPF, Wilkinson, inset patch
│   ├── phase-2/              # FDTD: cavity Q, horn, NTFF, openEMS cross-check
│   ├── phase-3/              # surrogate accuracy, NL-to-design
│   └── phase-4/              # FEM, MLFMA-scale problems
├── fixtures/                 # checked-in geometry + reference data (small)
└── results/                  # nightly outputs (gitignored)
```

## Case template

Every case is a self-contained Rust binary or Python script that:

1. Loads its geometry from `fixtures/` or builds it programmatically.
2. Runs the relevant solver via the public API.
3. Loads or computes the reference (closed-form / published-paper / cross-tool).
4. Asserts every metric against the documented tolerance.
5. Emits a plot (PNG/SVG via `plotters`) into `results/<case-id>/`.
6. Exits 0 on success; non-zero with a diagnostic on failure.

## Cross-tool validation

For every case where **openEMS** or **gprMax** can run the same geometry, we publish side-by-side numbers in `results/<case-id>/cross-tool.md`. Reproducible without trusting Yee's output.

## Phase 0 gates

| ID | Case | Tolerance |
|----|------|-----------|
| `v0-001` | Half-wave dipole impedance | ±5% vs Z ≈ 73 + j42 Ω |
| `v0-002` | 50 Ω microstrip Z₀ on FR-4 | ±3% vs TX-LINE / Hammerstad-Jensen |
| `v0-003` | 2.4 GHz rectangular patch resonance | ±2% vs published |

## Phase 1 gates

| ID | Case | Tolerance |
|----|------|-----------|
| `v1-001` | Swanson 5-pole hairpin BPF | ±1 dB to 4 GHz vs Sonnet ref |
| `v1-002` | Parallel-coupled-line BPF | ±1 dB |
| `v1-003` | Wilkinson divider 2 GHz | ±0.5 dB vs Pozar |
| `v1-004` | Branch-line 90° hybrid | amplitude + phase balance |
| `v1-005` | Inset-fed patch on RO4003C | figure-for-figure match |
| `v1-006` | Cross-validation vs openEMS, every microstrip/patch | ±3% at resonance |

Subsequent phases follow the same structure; see per-crate `validation/README.md` for the granular case list.
