# ADR-0183: S.6 S-parameters on the engine — two-run transmission over the job protocol

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0182 (S.5 materials on the protocol — the enabler), ADR-0108 (the
validated microstrip FDTD stack), ADR-0181 (the convergence path this serves).
**Spec:** `docs/superpowers/specs/2026-07-06-s6-engine-sparams-design.md`

## Context

The filter verify (F1.3: simulated response vs spec mask) needs |S21|(f) out of engine
jobs. S.5 put voxelized geometry on the job protocol; nothing yet turned probe series
into spectra in Rust (the studio's TS `dftMagnitude` is display math, not a reusable API).

## Decision

1. **New module `yee_engine::sparams`** — pure post-processing over `JobResult` probe
   series: `single_bin_dft` (the same correlation the ε_eff gates use, as a function) and
   `transmission_db` (20·log₁₀ magnitude ratio per frequency). Usable identically from
   tests, the studio, Python, and WS clients. **No protocol or `yee-compute` changes.**
2. **Method: two-run transmission ratio** (Sheen et al. 1990, adapted to lumped ports).
   Reference run = bare feed line; DUT run = line + device; both are ordinary S.5
   `JobSpec`s. The far end is terminated by a **passive resistive port** (`v0 = 0` leaves
   the pure-resistor arm of the validated lumped-port update — a 50 Ω load). Feed-line
   loss, launch discontinuity, probe coupling, and (to first order) termination mismatch
   divide out in the ratio.

## Gate

**engine-sparams-001** (`yee-engine/tests/sparams_stub_notch.rs`, `#[ignore]`, release
CI via the existing S.5 CI step): a **λ/4 open-circuited stub** on the S.5-certified
FR-4 stack — the textbook bandstop (Pozar): |S21| notches where the stub is a quarter
wave. The stub is sized `L_s = λ_g/4 − ΔL` (Hammerstad open-end correction), so closed
forms alone predict the notch at 5 GHz. Measured over the job protocol (2×~79 s solves,
~1.7 M cells × 9000 steps): **notch at 4.850 GHz, −36.8 dB deep — 3.0 % from the
transmission-line-theory prediction** (gate: ±15 %, ≥ 8 dB). Known walking-skeleton
artifact, bounded rather than hidden: out-of-band the stub still partially reflects and
those reflections exist only in the DUT run, so the single-probe ratio carries
standing-wave ripple of either sign at the band edges (measured +8.7 dB @3 GHz /
+5.2 dB @7 GHz; gate bounds |ripple| ≤ 12 dB). Cleaning that up is what the S11 /
incident-reflected separation follow-on is for. Fast unit gates: `single_bin_dft`
recovers a known sinusoid's amplitude/phase; `transmission_db` of a half-scaled copy is
−6.02 dB.

## Consequences

The full filter-verify chain now exists end-to-end on the engine: layout → voxelize →
`JobSpec` → CPU/GPU FDTD → probe series → |S21|(f) — every step protocol-visible, so the
studio, `yee-server` clients, and Python can compute device responses. Follow-ons: S11
(incident/reflected separation at the drive plane), spec-mask overlay + pass/fail (F1.3
proper), and a reusable layout→JobSpec bridge once a product surface consumes it.
