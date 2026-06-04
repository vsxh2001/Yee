# ADR-0161: CLI `--q-unloaded` — export the realistic (finite-Q) lumped response as Touchstone

**Status:** Accepted
**Date:** 2026-06-04
**Related:** ADR-0160 (finite-Q lumped response `ladder_s21_lossy` + Cohn gate), ADR-0158 (CLI lumped-board
export), [[project-filter-design-final-goal]], CLAUDE.md §4 (Touchstone is the project's primary external
interface), `yee_cli::filter::run_synth`.

---

## Context

`yee filter synth <spec> [--output out.s2p]` writes the synthesized filter's S-parameters as a Touchstone
`.s2p` — but via `yee_filter::ideal_response` (the **lossless** closed-form response). So the only
machine-readable response a user can export is the *ideal* curve (≈0 dB midband insertion loss). With
ADR-0160 the library can now compute the **realistic finite-Q response** (`ladder_s21_lossy`, validated
against Cohn's dissipation-loss formula to 2.2 %), but it is **not user-reachable** — there is no CLI path
to emit the realistic insertion-loss curve.

Touchstone is the project's primary external interface (CLAUDE.md §4), and the realistic response (what a
built filter actually measures) is a core deliverable of a filter-design tool. This closes the
library→CLI→Touchstone gap for the finite-Q response, mirroring how ADR-0158 made the lumped *board*
user-reachable.

## Decision

Add `--q-unloaded <Q>` (alias for the per-resonator unloaded quality factor) to `yee filter synth`:

- **Unset (default):** unchanged — the `.s2p` sweep uses `ideal_response` (lossless). Backward-compatible;
  every existing invocation and the existing gates are byte-identical.
- **Set to a finite `Q > 0`:** the `.s2p` (and the optional `--plot`) sweep uses the **finite-Q lumped
  realization's** response — `synthesize_lumped(&proj)` → `ladder_s21_lossy(&ladder, f, Q)` per frequency
  — so the exported Touchstone carries the realistic insertion-loss / rounded-corner curve. The CLI prints
  the realized midband insertion loss for visibility. The response model under `--q-unloaded` is the
  lumped-LC ladder (the same realization `--lumped` exports as a board), so a `--lumped --q-unloaded 100`
  run yields a coherent board + its realistic response.

This is a flag on the existing `run_synth` path (no new subcommand); the lossless and finite-Q paths share
the same `write_s2p` / sweep plumbing.

## Validation — gate `cli-finite-q-s2p` (fast, pure-compute, non-`#[ignore]`'d)

Run `run_synth` with `q_unloaded = Some(100)` for the standard 3-pole 0.5 dB Cheb spec, **read the written
`.s2p` back via `yee_io::touchstone::read`**, and assert:

1. the midband (`f0`) `|S21|` insertion loss ≈ Cohn's `4.343·Σg/(Q_u·FBW)` (≈ 1.86 dB) within the
   narrowband tolerance (≤ 15 %) — the realistic loss is present and correct (non-circular: Cohn from
   Σg/Q_u/FBW; the `.s2p` from the independent CLI sweep + round-trip);
2. the finite-Q `.s2p` **byte-differs** from the default (ideal) `.s2p` for the same spec — proving the
   `--q-unloaded` branch is actually taken (the ADR-0158 byte-diff pattern);
3. the default (no `--q-unloaded`) `.s2p` is byte-identical to before (regression guard via the existing
   path).

## Consequences

- The realistic filter response is user-reachable as Touchstone — the app/notebooks/CI can consume the
  finite-Q curve, not just the ideal one. Plotting ideal-vs-finite-Q in the GUI stays
  framework-verdict-gated (out of scope).
- Scope: `crates/yee-cli/src/filter.rs` (the flag + response routing) + the clap `Synth` command def +
  `crates/yee-cli/tests/` (the gate). Small; this ADR is the design record (no separate spec).
- **Not in scope:** distinct `Q_L`/`Q_C` (ADR-0160 single `Q_u` suffices), a lossy coupling-matrix
  `ideal_response` variant (the lumped ladder is the realization model), GUI plotting.

## References
- CLI synth path: `yee_cli::filter::run_synth` (`crates/yee-cli/src/filter.rs`), `write_s2p`,
  `yee_io::touchstone`.
- Response model: `yee_filter::ladder_s21_lossy` / `synthesize_lumped` (ADR-0160).
- Benchmark: Cohn dissipation-loss `4.343·Σg/(Q_u·FBW)` (ADR-0160 references).
