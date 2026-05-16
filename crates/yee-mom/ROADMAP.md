# yee-mom — Roadmap

## Phase 0 (months 0–6)
- [ ] `Solver` trait impl skeleton returning `Unimplemented`
- [ ] `SParameters` n-port container with Touchstone round-trip (`yee-io`)
- [ ] CPU reference fill for half-wave dipole (free-space Green's function, point matching)
- [ ] CPU dense LU via `faer`
- [ ] Validation case `mom-001` (half-wave dipole) green
- [ ] Stub GPU port behind `cuda` feature

## Phase 1 — beachhead (months 6–18)
- [ ] **Green's functions:** spectral-domain MPIE, DCIM extraction, rational fit, direct Sommerfeld fallback; adaptive switching
- [ ] **Basis functions:** RWG on triangles; rooftop on rectangles; consistent dual-basis testing
- [ ] **Ports:** delta-gap; edge; microstrip wave port with cross-section modal solve; CPW wave port
- [ ] **De-embedding:** TRL + SOLT, reference plane shift, port renormalization
- [ ] **Loss models:** finite conductor σ; Hammerstad-Jensen / Groiss / Huray roughness; dielectric tan δ
- [ ] **GPU matrix fill:** kernel per RWG-pair batch; cuBLAS aggregation
- [ ] **GPU dense LU:** cuSOLVER `Zgetrf` + `Zgetrs`; iterative refinement
- [ ] **Iterative path:** GMRES on GPU with block-diagonal preconditioner (n ≥ 50k)
- [ ] **Multi-GPU dense LU:** cuSOLVERMg behind `multi-gpu` feature
- [ ] **Validation:** every Phase 1 case below passes nightly

## Phase 4 — beyond beachhead
- [ ] MLFMA (multilevel fast multipole) for n ≥ 100k
- [ ] ACA (adaptive cross approximation) for off-diagonal blocks
- [ ] H-matrix compression
- [ ] Adjoint MoM for inverse design

## Validation gates — Phase 0
| ID | Case | Tolerance |
|----|------|-----------|
| mom-001 | Half-wave dipole impedance | ±5% vs Z ≈ 73 + j42 Ω |
| mom-002 | 50 Ω microstrip Z₀ on FR-4 | ±3% vs TX-LINE |
| mom-003 | 2.4 GHz rectangular patch resonance | ±2% vs published |

## Validation gates — Phase 1
| ID | Case | Tolerance |
|----|------|-----------|
| mom-101 | Swanson 5-pole hairpin BPF | ±1 dB vs Sonnet; resonances ±0.5% |
| mom-102 | Parallel-coupled-line BPF (Hong & Lancaster Ch.5) | ±1 dB |
| mom-103 | Wilkinson divider @ 2 GHz | ±0.5 dB vs Pozar |
| mom-104 | Branch-line 90° hybrid | amplitude + phase balance verified |
| mom-105 | Inset-fed patch on RO4003C | figure-for-figure published-paper match |
| mom-106 | Cross-validation vs openEMS, all microstrip/patch | ±3% at resonance |

## Risks
- DCIM accuracy across wide bands is finicky → ship multiple Green's-function evaluators + adaptive switch.
- Dense LU at n ≥ 50k overflows a single 80 GB H100 → iterative GMRES is the n ≥ 50k path; MLFMA deferred to Phase 4.
- Sonnet reference data is paywalled in places → cite peer-reviewed paper figures, not Sonnet-internal numbers.
