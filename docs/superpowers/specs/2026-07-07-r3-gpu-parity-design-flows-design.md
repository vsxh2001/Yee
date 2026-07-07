# R.3 — GPU parity for the design flows (aperture ports + per-face CPML)

**Date:** 2026-07-07
**Track:** RF-TOOL-ROADMAP R.3
**Related:** ADR-0187 (aperture ports, CPU), ADR-0192 (per-face CPML, CPU),
ADR-0177 (E.2 drive layer + GPU drive buffers), ADR-0178 (E.3 FP32 precision
policy), ADR-0195 (R.2).

## Problem

Every design-flow scenario (S.8–S.12 filters, A.0–A.3 antennas, R.0–R.2 board
physics) runs through **aperture ports** and, for antennas, **per-face CPML** —
both of which `GpuFdtd::with_drive` rejects with `ComputeError::Unsupported`.
`BackendChoice::Auto` silently falls back to CPU, so the design loops never
touch the GPU: the engine's 20× dGPU target is unreachable for exactly the
workloads the tool exists for.

## Scope

1. **Per-face CPML in WGSL.** The shader's `pml_depth` gates on a per-axis
   3-bit mask; the CPU reference (`cpml.rs::pml_depth`) gates min/max faces
   independently. Replace `Params.axes_mask` with a 6-bit `faces_mask`
   (bit `2·axis + side`, side 0 = min, 1 = max), built host-side from
   `CpmlConfig.faces`. The WGSL `pml_depth` checks the min bit in its
   `i < npml` branch and the max bit in its `i ≥ n − npml` branch. Profile
   arrays are depth-indexed and face-agnostic — no other shader change.
   Remove the `faces_are_axis_symmetric` rejection.

2. **Aperture ports in WGSL.** Verbatim port of the CPU loop (itself the
   verbatim `LumpedRlcPort::correct_e_aperture` pure-R arm): per port,
   `V*_T = (Σ_cells E_z·dz)/n_col`, `V_mid = ½(V*_T + V_prev)`,
   `I = (V_mid − V_src)·g` with `g = 1/(R + β)`, `β = dt·h/(2ε₀A)` (g = 0
   for an open port), sheet back-action `E_z −= (dt/(ε₀A))·I` on every cell,
   then `V_prev = (Σ_cells E_z·dz)/n_col` **re-summed** (matching the CPU's
   explicit resum, not the algebraic shortcut). One invocation per port —
   ports have O(10–100) cells, and a serial per-port loop keeps the
   semi-implicit update atomically consistent without cross-workgroup
   synchronization. Port cell sets are disjoint by construction (one port
   per i-plane), so port-parallel invocations never race.

   Buffer layout (append-only, so every existing `drv_*`/`dd_*` accessor is
   untouched):
   - `drv_idx` gains, after the probe offsets:
     `[n_ap, n_cells ×n_ap, cells_start ×n_ap, cell field-offsets ...]`
     (starts relative to the cells base).
   - `drv_data` gains, after the probe region:
     `[v_prev ×n_ap, vcoef ×n_ap, g ×n_ap, back ×n_ap, v_src (max_steps × n_ap)]`
     with `vcoef = dz/n_col`, `back = dt/(ε₀A)`, all f64-precomputed then
     narrowed, exactly like the existing port constants.
   - New entry point `apply_aperture_ports`, dispatched after `apply_ports`
     and before `record_probes` — the CPU ordering.

3. **Out of scope (documented):** the engine's protocol NTFF stays CPU-only —
   it steps a scratch grid one step at a time through the host adapter; the
   GPU path for far fields is yee-compute's own on-GPU DFT accumulator
   (E.5b), and wiring the engine's NtffSpec to it is its own increment.
   Conductor loss (R.0b) unchanged.

## Validation gates (both self-skipping without an adapter, real on llvmpipe
and the GPU nightly)

- **compute-015** `gpu_aperture_parity.rs`: a miniature S.10 board — per-cell
  ε_r substrate, PEC trace/ground masks, CPML-xy, a driven and a matched
  aperture port, E_z probes — stepped N steps on CPU (FP64) and GPU (FP32);
  probe series compared family-relative (the E.3 idiom): rel L2 < 1e-4,
  rel L∞ < 1e-3 vs the probe-family norms.
- **compute-016** `gpu_perface_cpml_parity.rs`: an open-top box
  (`faces = [[t,t],[t,t],[f,t]]`, PEC ground at z-min via the PecBox-side
  host zeroing being absent — the asymmetric case the shader could not
  express), Gaussian E_z ball, full-field family-relative comparison after N
  steps (the compute-005 idiom) plus evidence the z-max face actually
  absorbs (H-family norm well below a PEC-box run).

CI: both ride the existing `cargo test -p yee-compute --release --
--include-ignored` step (llvmpipe) and the GPU nightly.

## Consequences

`BackendChoice::Gpu`/`Auto` runs the actual design flows on the GPU; nightly
perf numbers for the design scenarios become measurable (the 20× target).
The FP32 drift policy (ADR-0178) applies unchanged: design loops that need
FP64 stay on CPU; the GPU is for the wide sweeps.
