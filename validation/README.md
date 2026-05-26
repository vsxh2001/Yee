# Yee — Workspace-level Validation

This directory is the **single source of truth** for end-to-end, cross-solver
validation. Per-crate validation lives in `crates/<name>/validation/` and
asserts crate-local correctness; this directory is for cases that exercise the
whole pipeline (mesh → solve → I/O → comparison vs reference).

Rule, enforced by CI:

> **No solver feature ships without a published-benchmark validation case and a
> CI run that regenerates results nightly.**

## Registered cases (`yee validate all`)

The live registry is `yee_validation::Report::run_all` in
`crates/yee-validation/src/lib.rs`. The table below reflects the **actual**
registry as of the current codebase; run `yee validate all --json` for the
ground truth at any commit.

### Method-of-Moments (MoM)

| ID | Case | Reference | Tolerance | Status |
|----|------|-----------|-----------|--------|
| `mom-001` | Half-wave dipole impedance, 24×176 cylinder mesh | NEC-4 finite-radius `Z ≈ 87 + j41 Ω` | Re ±5 %, Im ±10 % | **PASS** |
| `mom-002` | 50 Ω microstrip Z₀ on FR-4 (L = 82 mm) | Hammerstad-Jensen `Z₀ ≈ 51 Ω` | tripwire `\|Z_in\| ≤ 100 kΩ` | **PASS** (loose) |
| `mom-003` | 2.4 GHz rectangular patch on FR-4 | published resonance | loose tolerance | **PASS** (loose) |

> **IMPORTANT (CLAUDE.md §4 / ADR-0005):** the mom-001 reference is the
> NEC-4 *finite-radius* value `87 + j41 Ω` only. Do NOT cite the Balanis
> wire-limit approximation (zero-radius, ~20 % lower resistance) — see
> CLAUDE.md §4 and ADR-0005 for the rationale.

### FDTD

| ID | Case | Reference | Tolerance | Status |
|----|------|-----------|-----------|--------|
| `cpml-001` | CPML reflection ≥ 30 dB vs PEC | FDTD self-reference | ≥ 30 dB | SKIP (driver deferred) |
| `ntff-001` | NTFF broadside/endfire null ≥ 20 dB | FDTD self-reference | ≥ 20 dB | SKIP (driver deferred) |
| `dispersive-001` | Drude slab Fresnel reflection within 20 % | Fresnel analytic | ≤ 20 % | SKIP (driver deferred) |

FDTD sub-gate tests (cpml_reflection, ntff_dipole, drude_slab) live in
`crates/yee-fdtd/tests/` and are exercised by `cargo test -p yee-fdtd`.
The `fdtd-201` (TE₁₀₁) and `fdtd-201.x` (TE₂₀₁) cavity resonance gates
ship as `#[ignore]`'d tests in `crates/yee-fdtd`; aggregator drivers are
queued for a future yee-fdtd lane slice.

### FEM Eigenmode

| ID | Case | Reference | Tolerance | Status |
|----|------|-----------|-----------|--------|
| `fem-eig-001` | WR-90 cavity TE₁₀₁ eigenmode | Pozar §6.3 analytic | ±0.3 % (mode 1), RMS ±1 % (modes 1-10) | **PASS** |
| `fem-eig-002` | Lossy SiO₂ dispersive cavity TE₁₀₁ | closed-form Drude dispersion | Re ±0.5 %, Im ±5 % | **PASS** |
| `fem-eig-003` | WR-90 stub + CFS-PML absorption floor | Pozar §3.3 analytic S11 | ≤ −40 dB | SKIP (wall-time ~31 min) |
| `fem-eig-004` | WR-90 two-port thru-line `\|S21\|` / `\|S11\|` / reciprocity | lossless thru analytic | `\|S21\|` ±0.1 dB, `\|S11\|` < −20 dB, reciprocity < 1e-3 | **PASS** |
| `fem-eig-005` | WR-90 T-junction 3-port passivity + reciprocity | unitary S-matrix identity | passivity ≤ 1+ε, reciprocity < 1e-3 | **PASS** |
| `fem-eig-006` | High-aspect WR-90 wave-port `\|S11\|` < 0.1 | matched wave-port | `\|S11\|` < 0.1 | SKIP (gate open, queued Phase 4.fem.eig.3.5.6) |

Run `yee validate fem` to execute the FEM suite (fast gates only; 003 and
006 are Skipped in the default path).

## Cross-tool validation

For every case where **openEMS** or **gprMax** can run the same geometry, we
publish side-by-side numbers in `results/<case-id>/cross-tool.md`.
Reproducible without trusting Yee's output.

## Forthcoming gates (road-mapped, not yet registered)

| ID | Case | Phase |
|----|------|-------|
| `v1-001` | Swanson 5-pole hairpin BPF vs Sonnet | Phase 1 |
| `v1-002` | Parallel-coupled-line BPF | Phase 1 |
| `v1-003` | Wilkinson divider 2 GHz vs Pozar | Phase 1 |
| `v1-004` | Branch-line 90° hybrid | Phase 1 |
| `v1-005` | Inset-fed patch on RO4003C | Phase 1 |
| `fdtd-201` | Rectangular cavity TE₁₀₁ resonance | Phase 2 |
| `fdtd-201.x` | Rectangular cavity TE₂₀₁ resonance | Phase 2 |
