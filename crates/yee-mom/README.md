# yee-mom

> **Planar Method of Moments solver — the v1 beachhead of Yee.**

This is the crate where Yee earns its right to exist. The wedge is: **no open-source GPU planar MoM exists today.** Commercial Sonnet, Keysight Momentum, AWR AXIEM own this regime ($3k–$60k/seat). We take it.

## Theory of operation (Phase 1)

Mixed-potential integral equation (MPIE) formulation on planar layered substrates. Sommerfeld-type spectral-domain Green's functions accelerated by DCIM (Discrete Complex Image Method) and rational-function fitting. RWG / rooftop basis functions on triangular meshes from `yee-mesh`. Dense complex linear system filled on GPU (`cudarc::cublas` for aggregation, `cusolverDnZgetrf` for LU). Iterative GMRES path with block-diagonal preconditioner for n ≥ 50k.

References: Harrington, *Field Computation by Moment Methods* (1968 / 1993); Michalski & Mosig, "Multilayered media Green's functions in integral equation formulations," *IEEE T-AP* 45.3 (1997); Aksun, "A robust approach for the derivation of closed-form Green's functions," *IEEE T-MTT* 44.5 (1996); Pozar, *Microwave Engineering* (4th ed) §3.

## Scope

### Phase 0
- Crate skeleton, `Solver` impl returning `Unimplemented`
- `SParameters` container
- CPU dense LU via `faer` as the reference solver
- Lossless single-layer PEC test on a half-wave dipole

### Phase 1 — production planar MoM
- Multilayer dielectric stack-up via spectral-domain Green's functions
- RWG and rooftop basis functions
- Lumped (delta-gap, edge) ports; wave ports for microstrip / CPW with mode extraction
- TRL / SOLT de-embedding
- Surface roughness: Hammerstad-Jensen, Groiss, Huray small-sphere
- GPU matrix fill (one block per RWG-pair batch)
- GPU dense LU via cuSOLVER `Zgetrf` / `Zgetrs`
- Iterative GMRES on GPU for large n
- Multi-GPU via cuSOLVERMg (feature-gated)

### Phase 4 (not Phase 2/3)
- MLFMA / ACA / H-matrix compression for n ≥ 100k

## Validation (excerpt — see `validation/`)

- Half-wave dipole impedance Z ≈ 73 + j42 Ω, ±5%
- 50 Ω microstrip line on FR-4, ±3% vs TX-LINE
- 2.4 GHz rectangular patch on FR-4, ±2% resonance vs published
- **Swanson 5-pole hairpin BPF** (RT/Duroid 6006): ±1 dB vs Sonnet up to 4 GHz, ±0.5% resonance
- Wilkinson divider at 2 GHz: ±0.5 dB vs Pozar closed-form
- Cross-validate against openEMS on every microstrip / patch case (±3% at resonance)

## Roadmap

See [`ROADMAP.md`](ROADMAP.md).
