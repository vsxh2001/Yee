# ADR-0174: CLI default `.s2p` phase parity — use the complex coupling-matrix S (T11)

**Status:** Accepted (the ADR-0172 T9 follow-on)
**Date:** 2026-06-06
**Related:** ADR-0172 (T9 — `coupling_matrix_s_params`; the studio distributed `.s2p` got real phase, the CLI
was the noted follow-on), ADR-0171 (T8 — `.s2p`), ADR-0161 (`write_s2p` + `lossless_s_pair`),
[[project-filter-design-final-goal]].

---

## Context

`yee filter synth`'s default `.s2p` (no `--q-unloaded`) builds its S-parameters as
`ideal_response(&proj, &freqs).map(lossless_s_pair)` (`yee_cli::filter::run_synth`) — `|S21|` from the
characteristic-function magnitude with `S11` placed in lossless quadrature, i.e. **flat S21 phase**. T9 gave the
studio's distributed `.s2p` real phase via `coupling_matrix_s_params`, but the CLI default path still emits flat
phase — a CLI↔studio inconsistency (the noted T9 follow-on). The `--q-unloaded` path already emits the complex
lumped finite-Q pair (`ladder_s_params_lossy`) and is correct.

## Decision

In `run_synth`, the **default** (`q_unloaded == None`) `.s2p`/`--plot` branch uses
`coupling_matrix_s_params(&proj.coupling, &freqs, proj.spec.f0_hz, proj.spec.fbw)` (the T9 complex
coupling-matrix response — real phase) instead of `ideal_response().map(lossless_s_pair)`. The resulting
`(S11, S21)` is lossless/passive by construction (`|S11|²+|S21|² = 1`), so `write_s2p`'s contract + a `yee_io`
re-import still hold. The `--q-unloaded` (lumped finite-Q) branch is UNCHANGED. Update the `run_synth` doc +
the `.s2p` comment string (the default response is now "complex coupling-matrix S", not "lossless closed-form
`ideal_response`"). If `lossless_s_pair` becomes dead after the swap, remove it (else keep).

**Gate** (`crates/yee-cli/tests/`, extend the existing CLI `.s2p` test or add one; pure-compute/fast): the
default `.s2p` carries **non-flat S21 phase** (arg(S21) varies across the band — distinguishing it from the
old flat-phase output), **round-trips** through `yee_io::touchstone` parse, and is **passive**
(`|S11|²+|S21|² ≈ 1`). Non-circular (the phase comes from the coupling-matrix solve; the test asserts the
emitted file's properties).

## Consequences

- The CLI default `.s2p` carries physical phase — CLI↔studio parity (both now use `coupling_matrix_s_params`
  for the ideal/distributed response), closing the ADR-0172 follow-on.
- Scope T11: `crates/yee-cli/src/filter.rs` + `crates/yee-cli/tests/`. Pure-compute. This ADR is the record.
- **Not in scope:** the `--plot` magnitude (unchanged visually — `|S21|` agrees with `ideal_response` per the
  ADR-0172 gate); the lumped finite-Q path (already complex); changing `ideal_response` (kept for the mask).

## Outcome (T11 — SHIPPED, merge `c681a76`)

Done as specified (+1/−1 import, default branch one-liner, doc/comment updates, `lossless_s_pair` removed as
dead, +179 gate). `run_synth`'s default (`q_unloaded == None`) `.s2p`/`--plot` branch now calls
`coupling_matrix_s_params(&proj.coupling, &freqs, spec.f0_hz, spec.fbw)` (complex `(S11,S21)`, real phase,
lossless/passive by construction) instead of `ideal_response().map(lossless_s_pair)`. The `--q-unloaded`
finite-Q lumped path (`ladder_s_params_lossy`) is byte-for-byte UNCHANGED. The now-dead `lossless_s_pair`
helper was removed; the `run_synth` + `write_s2p` + `main.rs` `Filter`-subcommand doc comments were corrected to
describe the new default.

**Gate `cli-s2p-phase` (`crates/yee-cli/tests/cli_s2p_phase.rs`, non-`#[ignore]`'d, NON-circular):** runs the
default `yee filter synth` (3-pole 0.5 dB Cheb, f0 = 2 GHz, FBW = 0.10) with **no `--q-unloaded`**, reads the
written `.s2p` back through `yee_io::touchstone::read`, and asserts on the EMITTED file — (1) round-trips as a
2-port with a frequency grid; (2) **non-flat S21 phase** (`arg(S21)` span > 0.5 rad — the old flat-phase
default had `arg(S21) ≡ 0`, so this strictly discriminates T11's complex response); (3) **passive** everywhere
(`|S11|²+|S21|² ≤ 1+ε`) and **lossless** (`≈ 1`) at midband. Non-circular: the phase comes from the
coupling-matrix linear solve; the test only inspects the read-back artifact.

`cargo test -p yee-cli` green (incl. the new gate + the existing `cli_finite_q_s2p` / `cli_plot_touchstone`);
`cargo clippy -p yee-cli --all-targets -- -D warnings` clean; `cargo fmt -p yee-cli --check` clean. Reviewer
APPROVE after one P1 fix (a stale `main.rs:227` doc comment that still described `ideal_response` as the default —
corrected in-lane); reviewer confirmed argument order/types, the unchanged `--q-unloaded` arm, full removal of
`lossless_s_pair`, the gate's non-circularity + correct Touchstone row indexing (S11 = `row[0]`, S21 = `row[2]`).
**CLI↔studio parity achieved** — both surfaces now emit the same complex coupling-matrix S; the ADR-0172 T9
follow-on is closed.

## References
- `yee_filter::coupling_matrix_s_params` (ADR-0172); `yee_cli::filter::{run_synth, write_s2p}` (`lossless_s_pair`
  removed in T11); `yee_io::touchstone` (round-trip).
