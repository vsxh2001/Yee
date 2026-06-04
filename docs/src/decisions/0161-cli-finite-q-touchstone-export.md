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

- **Unset (default):** unchanged — the `.s2p` sweep uses `ideal_response` (lossless). The reflection `S11`
  is the **true lossless reflection** for a lossless reciprocal symmetric 2-port, written in quadrature
  with the (real) transmission: `S11 = j·√(1 − |S21|²)`. This makes the exported S-matrix both energy-
  conserving (`|S11|² + |S21|² = 1`) and *passive* (`σ_max = 1`), so it round-trips through
  `yee_io::touchstone::read`, which checks passivity. Backward-compatible; every existing invocation and the
  existing gates are byte-identical.
- **Set to a finite `Q > 0`:** the `.s2p` (and the optional `--plot`) sweep uses the **finite-Q lumped
  realization's TRUE lossy 2-port** — `synthesize_lumped(&proj)` → `ladder_s_params_lossy(&ladder, f, Q)`
  per frequency, which returns `(S11, S21)` from one total ABCD:
  `S21 = 2/Δ`, `S11 = (A + B/Z0 − C·Z0 − D)/Δ`, `Δ = A + B/Z0 + C·Z0 + D`. So the exported Touchstone
  carries the realistic insertion-loss / rounded-corner curve **and the true absorptive reflection**: a
  dissipative network has `|S11|² + |S21|² < 1` (the deficit is power absorbed in the resonator losses), and
  it is still passive so it round-trips. `S11` is **not** a lossless `√(1 − |S21|²)` placeholder — that would
  mis-attribute the dissipative insertion loss to reflection (`|S11| ≈ 0.58` at midband) and falsely claim
  the filter is lossless. The CLI prints the realized midband insertion loss and the midband
  `|S11|² + |S21|²` (absorbed fraction) for visibility. The response model under `--q-unloaded` is the
  lumped-LC ladder (the same realization `--lumped` exports as a board), so a `--lumped --q-unloaded 100`
  run yields a coherent board + its realistic response.

This is a flag on the existing `run_synth` path (no new subcommand); both paths share the same `write_s2p` /
sweep plumbing. `write_s2p` takes the per-frequency `(S11, S21)` pairs and writes them verbatim (it no longer
derives `S11` internally) — the caller owns the physics: quadrature lossless `S11` for the ideal path, the
true ABCD `S11` for the finite-Q path. `ladder_s21_lossy` is refactored to return `ladder_s_params_lossy(…).1`
so there is a single ABCD implementation and the `S21` magnitude (the Cohn-validated quantity) is unchanged.

## Validation — gate `cli-finite-q-s2p` (fast, pure-compute, non-`#[ignore]`'d)

Run `run_synth` with `q_unloaded = Some(100)` for the standard 3-pole 0.5 dB Cheb spec, **read the written
`.s2p` back via `yee_io::touchstone::read`**, and assert:

1. the midband (`f0`) `|S21|` insertion loss ≈ Cohn's `4.343·Σg/(Q_u·FBW)` (≈ 1.86 dB) within the
   narrowband tolerance (≤ 15 %) — the realistic loss is present and correct (non-circular: Cohn from
   Σg/Q_u/FBW; the `.s2p` from the independent CLI sweep + round-trip);
2. **absorption present** — at midband the finite-Q `.s2p` has `|S11|² + |S21|² < 0.999`, i.e. it is a true
   lossy 2-port and *not* the lossless placeholder (which forces `≡ 1`). This is the assertion that catches
   a fictitious `S11 = √(1 − |S21|²)`. (Measured: `|S11|² + |S21|² ≈ 0.665`, ~33 % absorbed at `Q = 100`.)
   Conversely the ideal `.s2p` has `|S11|² + |S21|² ≈ 1` (true lossless reflection);
3. the finite-Q `.s2p` **byte-differs** from the default (ideal) `.s2p` for the same spec — proving the
   `--q-unloaded` branch is actually taken (the ADR-0158 byte-diff pattern);
4. the default (no `--q-unloaded`) `.s2p` is byte-identical to the prior (passive, quadrature-`S11`) ideal
   output (regression guard via the existing path).

## Consequences

- The realistic filter response is user-reachable as Touchstone — the app/notebooks/CI can consume the
  finite-Q curve, not just the ideal one. Plotting ideal-vs-finite-Q in the GUI stays
  framework-verdict-gated (out of scope).
- Scope: `crates/yee-cli/src/filter.rs` (the flag + response routing) + the clap `Synth` command def +
  `crates/yee-cli/tests/` (the gate) + `crates/yee-filter/src/lumped.rs` (the `ladder_s_params_lossy`
  `(S11, S21)` helper + lib re-export, so the CLI can read the true lossy reflection without re-implementing
  the ABCD). Small; this ADR is the design record (no separate spec).
- **Not in scope:** distinct `Q_L`/`Q_C` (ADR-0160 single `Q_u` suffices), a lossy coupling-matrix
  `ideal_response` variant (the lumped ladder is the realization model), GUI plotting.

## References
- CLI synth path: `yee_cli::filter::run_synth` (`crates/yee-cli/src/filter.rs`), `write_s2p`,
  `yee_io::touchstone`.
- Response model: `yee_filter::ladder_s_params_lossy` (true lossy `(S11, S21)`) / `ladder_s21_lossy`
  (delegates to `.1`) / `synthesize_lumped` (ADR-0160).
- S-from-ABCD for equal real `Z0` terminations: Pozar, *Microwave Engineering* Table 4.2.
- Benchmark: Cohn dissipation-loss `4.343·Σg/(Q_u·FBW)` (ADR-0160 references).
