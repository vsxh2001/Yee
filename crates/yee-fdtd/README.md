# yee-fdtd

> 3D FDTD on the Yee staggered grid. **Phase 2 deliverable.**

This crate is a **stub during Phase 0 and Phase 1** so the workspace builds end-to-end and downstream crates can already depend on its public API surface. All public functions return `Error::NotYet` until Phase 2 begins (~month 18).

## Why FDTD when MoM is the beachhead?

MoM is the right solver for planar PCB structures because most of the geometry is flat copper on a stack of dielectrics. FDTD takes over where MoM is the wrong tool: full 3D antenna radiation patterns, transient signal integrity, dispersive materials (Drude / Lorentz / Debye), and EMC-style enclosure problems. The two solvers complement each other; the Yee grid that names this project lives in this crate.

## Scope (Phase 2)

- 3D Yee staggered grid on CUDA with memory-bandwidth-optimized E/H update kernels
- Mixed precision (FP32 fields, FP64 accumulators where needed)
- Multi-GPU domain decomposition via NCCL boundary exchange
- CPML (Roden & Gedney) absorbing boundaries on all six faces
- Dispersive materials: Drude, Lorentz, Debye, multi-pole Debye via ADE or PLRC
- NTFF (near-to-far-field) transformation — full 3D pattern, gain, directivity, axial ratio
- Lumped-element ports (R/L/C/RLC), waveguide ports with modal sources, plane-wave sources
- Subgridding (Berenger / Xiao-Liu) — ships *after* the rest of the solver is locked
- Conformal (Dey-Mittra) or staircase fallback for non-aligned geometry
- KiCad PCB / STEP / IGES ingestion via `yee-mesh` + voxelization
- GPU-resident time-stepping with live `rerun` streaming

## Validation (excerpt)

- Rectangular cavity TE/TM modes within ±0.5% of analytical
- Pyramidal horn pattern within ±1 dB in main beam
- Dipole over dielectric half-space NTFF vs Sommerfeld
- Cross-validation against openEMS on identical geometries
- Microstrip transient TDR matches FFT of MoM frequency response

## Roadmap

See [`ROADMAP.md`](ROADMAP.md).
