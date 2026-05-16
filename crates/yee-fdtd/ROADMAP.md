# yee-fdtd — Roadmap

## Phase 0 (months 0–6)
- [ ] Crate skeleton with `Error::NotYet` stub — keeps workspace green
- [ ] Public API sketch (`Grid`, `Material`, `Source`, `Probe`) drafted in `lib.rs` doc comments

## Phase 1 (months 6–18)
- [ ] Public API frozen at the trait level — no impls yet
- [ ] CPU 1D toy solver in `examples/` for documentation/teaching purposes only

## Phase 2 — production (months 18–30)
- [ ] **Grid:** 3D Yee staggered grid on CUDA; memory-bandwidth-optimized E/H kernels
- [ ] **Mixed precision:** FP32 fields + FP64 accumulators where stability requires
- [ ] **Boundaries:** CPML (Roden & Gedney) on all six faces with polynomial grading
- [ ] **Materials:** Drude / Lorentz / Debye / multi-pole Debye via ADE or PLRC
- [ ] **Sources:** plane-wave injection, Gaussian/modulated-Gaussian, modal waveguide
- [ ] **Ports:** lumped R/L/C/RLC, waveguide modal
- [ ] **NTFF:** full 3D pattern, gain, directivity, axial ratio
- [ ] **Subgridding:** Berenger / Xiao-Liu (after the rest of the solver is locked)
- [ ] **Conformal:** Dey-Mittra fractional cells; staircase fallback
- [ ] **I/O:** STEP/IGES via `opencascade-rs`; KiCad PCB; voxelizer via `yee-mesh`
- [ ] **Multi-GPU:** NCCL boundary exchange; feature-gated initially
- [ ] **Debug streaming:** live volume-data → `rerun` viewer

## Phase 4
- [ ] Time-domain FEM hybrid
- [ ] Anisotropic / nonlinear / time-varying materials

## Validation gates — Phase 2
| ID | Case | Tolerance |
|----|------|-----------|
| fdtd-201 | Rectangular cavity TE/TM Q-factor | ±0.5% vs analytical |
| fdtd-202 | Pyramidal horn antenna pattern | ±1 dB in main beam |
| fdtd-203 | Dipole over dielectric half-space NTFF | matches Sommerfeld reference |
| fdtd-204 | Cross-validation vs openEMS | numerical-noise level (grid/PML driven) |
| fdtd-205 | Microstrip transient TDR | matches FFT of `yee-mom` frequency response |

## Risks
- Memory bandwidth is the bottleneck. Beating openEMS requires hand-tuned shared-memory tiling — we will benchmark openly.
- Subgridding stability is famously fragile. **Ship without it first.**
- Multi-GPU domain decomposition adds complexity. Single-GPU first; multi-GPU behind `multi-gpu` feature with "experimental" label.
