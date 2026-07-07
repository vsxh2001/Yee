# ADR-0196: R.3 — GPU parity for the design flows (aperture ports + per-face CPML)

**Status:** Accepted
**Date:** 2026-07-07
**Related:** RF-TOOL-ROADMAP R.3, ADR-0187 (aperture ports), ADR-0192 (per-face
CPML), ADR-0177 (E.2 drive buffers), ADR-0178 (E.3 FP32 precision policy).
**Spec:** `docs/superpowers/specs/2026-07-07-r3-gpu-parity-design-flows-design.md`

## Decision

The two `ComputeError::Unsupported` rejections that kept every design-flow
scenario (S.8–S.12 filters, A.0–A.3 antennas, R.0–R.2 board physics) off the GPU
are gone; the WGSL backend now runs both features, certified differentially
against the bit-exact CPU path.

- **Per-face CPML.** `Params.axes_mask` (3 bits) became `faces_mask` (6 bits,
  `bit 2·axis + side`), built host-side from `CpmlConfig.faces`; the shader's
  `pml_depth` gates its min-face branch on bit `2·axis` and its max-face branch
  on bit `2·axis + 1`. Profiles are depth-indexed and face-agnostic — nothing
  else changed. `with_axes` still sets both faces, so every pre-R.3 upload is
  bit-identical.
- **Aperture ports.** New `apply_aperture_ports` entry point — a verbatim port
  of the CPU update (itself bit-exact vs `LumpedRlcPort::correct_e_aperture`,
  gate compute-014): modal `V*_T = (dz/n_col)·Σ E_z`, semi-implicit midpoint
  against the cached `V_prev`, branch current `I = (V_mid − V_src)·g` with
  `g = 1/(R + β)` precomputed host-side (0 for the open-port arm), sheet
  back-action on every cell, then `V_prev` explicitly re-summed. **One
  invocation per port**, serial loops over its O(10–100) cells: ports own
  disjoint cell sets, so port-parallel invocations never race, and the serial
  loop needs no cross-workgroup synchronization. Buffer layout is append-only
  (cell table at the end of `drv_idx`; `v_prev`/`vcoef`/`g`/`back` constants and
  the v_src series at the end of `drv_data`), so every existing accessor offset
  is untouched. Dispatched between `apply_ports` and `record_probes` — the
  reference order.

The engine's GPU job path needed **zero changes**: it already forwarded aperture
ports and face masks and only ever hit the rejections; `BackendChoice::Gpu` now
runs the design flows, and `Auto` stops falling back.

**Out of scope:** the engine-protocol NTFF stays CPU-only (it drives the
validated `yee_fdtd::NtffState` through the host adapter one step at a time;
wiring `NtffSpec` to yee-compute's on-GPU DFT accumulator (E.5b) is its own
increment). Conductor loss (R.0b) unchanged.

## Gates (llvmpipe, 2026-07-07; both self-skip without an adapter)

- **compute-015** (`gpu_aperture_parity.rs`): miniature S.10 board — FR-4
  substrate, PEC ground + trace masks, CPML-xy, driven + matched aperture
  ports, 400 steps, three E_z probe series vs CPU. Measured family-relative
  **L2 ≤ 1.2e−7, L∞ ≤ 1.6e−7** (gates 1e−3 / 5e−3) — ulp-class, i.e. the
  kernel is algorithmically identical and only FP32 rounding separates them.
- **compute-016** (`gpu_perface_cpml_parity.rs`): open top over a reflective
  floor (`faces = [[t,t],[t,t],[f,t]]`), 150 steps, full-field comparison.
  Measured **E ≤ 2.1e−7, H ≤ 6.7e−5 L2 / 1.7e−4 L∞** (gates 1e−4 / 1e−3);
  absorption evidence: post-run ‖H‖₂ = **2.5e−5 vs 1.0e−2** in a PEC box —
  the five enabled faces absorb, the disabled floor reflects.

## Consequences

Filter and antenna design loops can now run their sweeps on
`BackendChoice::Gpu`; the GPU nightly can time the actual design scenarios
(the 20× dGPU target). The FP32 policy (ADR-0178) is unchanged: single-run
certification stays on the FP64 CPU path; the GPU is for wide sweeps and
interactive iteration. Next: R.4 (BPF end-to-end + surrogate BO).
