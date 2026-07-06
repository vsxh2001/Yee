# ADR-0179: engine tickets closed (E.4/E.5b/E.5c) + S.0 job API + Tauri/React studio skeleton

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0175..0178, `ENGINE-STUDIO-ROADMAP.md`. This is the studio-kickoff ADR the
roadmap's S.2 row required.

---

## 1. Engine tickets closed

**E.4 (kernel optimization).** `update_h`/`update_e` inner loops rewritten over fixed-length
row slices with the material branch hoisted to row level (bounds-check elision +
vectorization); per-cell arithmetic untouched — compute-001/003/007 stay bit-exact. Measured
(4-core container, 64³): single-thread 11.1 → 8.1 ms (−27 %), 4-thread 4.90 ms ≈ the scalar
reference (memory-bandwidth-bound). Further numbers come from real hardware via the nightly.

**E.5c (dispersive ADE).** Drude/Lorentz/Debye on both backends. CPU: verbatim `ade_step`
port, slab-parallel; gate `compute-011` **bit-exact** vs `yee_fdtd::dispersive` on a four-arm
scenario. GPU: unified ADE form with f64-materialized coefficient maps folded into the
existing coeff/psi arenas (8-binding limit still holds). Gate `compute-012` is
**differential**: ADE-pair GPU↔CPU error ≤ 20× the standard-pair error on the identical
scenario (measured ≤ 6× — the f32 aux recursions amplify rounding; the all-vacuum ADE path is
bit-identical to the standard path, ruling out indexing bugs) + a drift-class absolute
backstop. Bring-up surfaced and documented a normalization pathology: the scenario's ‖H‖
stays ~η₀ below ‖E‖ (absorbers eat the wave, the electrostatic ball persists), so naive
family-relative bounds read ~2e-4 for *any* FP32 path, standard included.

**E.5b (GPU-side NTFF).** The full-field DFT phasor at `f_probe` is accumulated **on the
GPU** every step (an `accumulate_dft` kernel over the field arena; phase from the on-GPU step
counter — no per-step readback), stored at the tail of the psi arena, read back once, and fed
to the reference `NtffState` through two synthetic samples (its accumulation is linear:
`e^{-jωt}=1` picks up Ê_re, `e^{-jωt}=+j` picks up j·Ê_im — exactly `Σ F·e^{-jωt}·Δt`). Gate
`compute-013`: the dipole scenario GPU-resident — broadside/endfire **315.4 dB** (analytic
sin θ null) and broadside magnitude matches the CPU host-adapter path to **2.9e-7** relative.

## 2. S.0 — `yee-engine` job API

New workspace crate: serde protocol (`JobSpec` / `JobEvent` / `JobResult`) + in-process
executor. `submit(spec) → JobHandle` spawns a worker thread that runs the job in ~2 % chunks,
streams `Progress`, honors cooperative `cancel()`, and finishes with `Done`/`Error`. Backend
selection (`cpu` / `gpu` / `auto`-with-fallback) mirrors the bindings. S.0 scope: driven
vacuum jobs (any boundary, soft sources, resistive ports, probes); materials/dispersion plumb
through in a later slice. Gated by 4 unit tests + a doctest (serde round-trip, progress
streaming, cancellation, auto-backend).

## 3. S.2 walking skeleton — Tauri 2 + React studio (`studio/`)

- **Stack:** Tauri 2 (Rust host, system webview) + React 18 + TypeScript + Vite. Frontend
  bundle: **147.6 kB JS / 47.9 kB gzipped** — the small-footprint alternative to Electron the
  direction called for (ADR-0175).
- **Architecture:** the webview calls one Tauri command, `run_job(spec: JobSpec)`, which runs
  `yee_engine::submit` on a blocking task, re-emits engine `Progress` as `job://progress`
  events, and resolves with the `JobResult`. The frontend is transport-agnostic by
  construction — `yee-server` (S.1) will carry the identical serde types over WebSocket.
- **UI (walking skeleton):** grid/steps/backend controls → run → live progress bar → probe
  time series plotted as a dependency-free inline SVG. Charting/3D (S.3) comes later.
- **`studio/src-tauri` is deliberately outside the root workspace** (`[workspace]` detach):
  the webkit2gtk dependency tree must not weigh down workspace-wide builds/clippy/CI.
  Verified in-container: `cargo check` green against libwebkit2gtk-4.1-dev, `npm run build`
  green; interactive `tauri dev` needs a display and stays a local-machine step.
- The Dioxus studio freeze from ADR-0175 remains in force; the S.4 parity audit decides
  retirement.

## 4. Consequences

- Every engine roadmap ticket through E.5 is now closed with a measured gate; the engine
  track is **complete** pending real-hardware perf numbers (nightly) and future physics.
- The studio track has its foundation: protocol (S.0) + shell (S.2 skeleton). Next: S.1
  `yee-server`, S.3 visualization, S.4 parity audit.
- `yee-py` (ADR-0178) and the studio consume the same engine — one validated core, three
  surfaces (Python, desktop, future web).
