# FS.2 — far-field products: gain, efficiency, full-sphere export

**Date:** 2026-07-08
**Track:** FULL-SUITE-ROADMAP FS.2, opened as FS.1a/FS.1b closed. The
NTFF path (E.5/E.5b, A.2, engine-antenna-003/006/008) returns **relative**
|E| per direction — every pattern gate so far asserts shape ratios. What
commercial deliverables need is **absolute**: gain in dBi, radiation
efficiency, and a full-sphere export.

## The missing physics: input-power normalization

Gain = 4π·U(θ,φ) / P_in with U = r²|E|²/(2η₀). The NTFF already gives
|E| at unit reference distance; the missing quantity is **P_in at the
port**. The aperture port already computes its modal terminal voltage
`V_T` every step (`aperture_state`, the LumpedRlcPort correction); the
branch current through the source resistance is `(V_emf − V_T)/R`. So the
port can record per-step `(v_t, i_branch)` and the engine can integrate
`P_acc = Σ v·i·dt / T_pulse`-style accepted energy — no new field probes,
no modal-normalization guesswork.

## Decomposition

- **FS.2a (walking skeleton): port records on the protocol.**
  `AperturePortSpec::record: bool` (serde default false) → per-step
  `(v_t, i)` series in `JobResult.port_records`. CPU first; GPU keeps
  `Unsupported` rejection for recording ports until the readback buffer
  lands (own increment, the R.3 idiom). Gate `engine-power-001`: on the
  R.0 through-line, accepted energy at port A ≈ delivered energy at
  port B's load + (tiny) CPML leakage — energy bookkeeping closes to a
  measured-then-pinned tolerance on a lossless line.
- **FS.2b: gain in dBi.** `farfield::gain_dbi(|E|, r_ref, p_in)`;
  gate `engine-gain-001`: the A.1 patch's broadside gain lands in the
  textbook 5–8 dBi window for a thin FR-4 patch, and the FS.1b 2×1 array
  reads ~2.5–3 dBi above the single element (the array-gain identity, a
  differential assert that cancels most modeling bias).
- **FS.2c: radiation efficiency + full-sphere export.** Efficiency = 1
  sanity on the lossless stack (measured-then-pinned band); with R.0
  tan δ and R.0b sheet loss enabled the efficiency drops accordingly
  (qualitative direction gate). Full-sphere (θ, φ) raster export to CSV
  via yee-plotters/yee-io; byte-checked artifact.

## Risks

- The NTFF's |E| reference distance/normalization convention must be
  audited once against the validated `yee_fdtd::NtffState` docs before
  FS.2b quotes absolute dBi (the sin θ dipole gate certified *shape*).
- GPU port recording is deferred, not skipped — the design-loop flows
  run CPU today.
