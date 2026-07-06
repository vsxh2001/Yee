# ADR-0176: E.1 — CPML + per-cell materials on `yee-compute` (both backends)

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0175 (engine direction; E.0 walking skeleton), `ENGINE-STUDIO-ROADMAP.md` phase
E.1, spec `docs/superpowers/specs/2026-07-06-e1-cpml-materials-design.md`, plan
`docs/superpowers/plans/2026-07-06-e1-cpml-materials.md`.

---

## Context

E.0 shipped `yee-compute` with uniform lossless vacuum only. Real workloads — starting with the
filter pipeline's `fdtd-line-eeff-001` (E.2) — need absorbing boundaries, dielectric substrates,
loss, and PEC geometry. All four exist in `yee-fdtd` (Roden–Gedney 2000 CPML, per-cell
ε_r/μ_r/σ with the Taflove §3.7 CA/CB arm, interior PEC masks, legacy outer PEC clamp); E.1
ports them to both `yee-compute` backends under the established gating discipline: CPU bit-exact
against the reference, GPU tolerance-gated against the CPU.

## Decision

1. **Public surface:** `Materials` (optional per-cell ε_r/μ_r/σ maps `[nx+1,ny+1,nz+1]` +
   per-component PEC masks, `YeeGrid` conventions) and `Boundary { None, PecBox, Cpml(CpmlConfig) }`,
   consumed by `with_config` constructors on `CpuFdtd`/`GpuFdtd`/`FdtdEngine`. `Boundary::None`
   preserves E.0 raw-kernel semantics. The CPU backend also gets `step_with_gaussian_ez` (the
   exact `sources::gaussian_pulse_ez` soft source + `step_with_source` ordering) as E.1 test
   plumbing; engine-level sources/ports remain E.2.
2. **CPU:** flat-buffer, rayon-slab ports with per-cell arithmetic, branch structure, and branch
   order identical to the reference (`update.rs` arms, `cpml.rs` passes, `apply_pec` clamp,
   `apply_pec_mask` final word). ψ arrays share the written component's shape, so CPML passes
   slab-zip (field, ψ_a, ψ_b) and stay bit-exact.
3. **GPU:** **arena-buffer refactor** — materials, masks, 12 ψ arrays, and profiles as separate
   bindings would blow WebGPU's default 8-storage-buffers-per-stage limit, so the backend packs
   five arenas (fields; ca/cb/ce_cpml/ch coefficient maps; ψ; CPML profiles; masks) + one
   uniform. Materials are never branched on in WGSL: the host materializes the four coefficient
   maps in f64 (ca = 1 for lossless cells makes `e = 1·e + cb·curl` reproduce the plain add
   exactly; `ce_cpml = Δt/(ε₀ε_r)` mirrors the reference CPML pass ignoring σ) and narrows
   once. Bulk + CPML are fused into the six update kernels (algebraically identical to the
   reference's two passes); three clamp kernels apply masks after the E half; the PEC box is a
   host-side invariant (faces zeroed at upload, never written by any kernel).

## Gates (all green)

- **compute-003** (`tests/cpu_e1_reference_parity.rs`): heterogeneous scenario — dielectric
  slab ε_r = 4.3, lossy block σ = 0.5 (CA/CB arm), μ_r = 2 band, PEC sheet with slot (masks),
  driven Gaussian source, 30 steps, 24×20×22 — **bit-exact** (max |Δ| == 0.0, all six
  components) vs `WalkingSkeletonSolver` in **both** CPML and PEC-box modes.
- **compute-004** (`tests/cpml_reflection.rs`): the `yee-fdtd` reflection methodology on
  `CpuFdtd` (50³, npml = 10, 300 steps) — measured **69.3 dB** reduction vs the 30 dB target;
  runs in ~2 s debug (the rayon backend vs the reference's ~60 s single-thread release budget).
- **compute-005** (`tests/gpu_e1_parity.rs`): GPU vs CPU on CPML + all materials + mask,
  100 steps — family-rel L2 ≈ 2–3e-7 (E) / 2–3e-6 (H) on Mesa llvmpipe, tolerances 1e-4/1e-3;
  absorption evidence on the **H-family norm only** (the initial E_z ball's electrostatic
  residual is curl-free and unabsorbable — total-field norms barely drop by design): CPML holds
  **210× less** ‖H‖ than the PEC box after 100 steps. Self-skips without an adapter; runs on
  the GPU nightly.
- compute-001/002 (E.0) unchanged and green — the arena refactor did not perturb the vacuum
  path (compute-002 numbers identical to E.0's).

## Consequences

- `FdtdEngine` variants are boxed (`Box<CpuFdtd>` / `Box<GpuFdtd>`) — the CPU stepper now
  carries material maps + CPML state by value.
- The GPU binding layout is fixed at 5 storage + 1 uniform, inside WebGPU browser limits — the
  future wasm path keeps working without limit negotiation.
- The reference CPML's σ-blindness in its E-pass coefficient (`Δt/(ε₀ε_r)`, not CB) is
  mirrored, not "fixed" — parity with `yee-fdtd` outranks local judgement; revisit upstream if
  ever warranted.
- E.2 (sources/ports + driver parity, `fdtd-line-eeff-001` on the engine) is unblocked.
