# ADR-0065 — GUI multi-trace S-parameter overlay (egui_plot)

**Status:** Accepted
**Date:** 2026-05-24
**Context Phase:** frontend (yee-gui) — follow-on to ADR-0063

## Context

The static plotters + CLI gained multi-trace / off-diagonal S-parameter
overlay (ADR-0063), but the desktop GUI still renders a single trace
(`plots.rs:70`), so a multi-port device cannot show S21 or an S11+S21
overlay in the app. With the MoM-port accuracy thread closed as
ill-posed (ADR-0064) and the FDTD energy-balance / FEM-port tracks
deferred quagmires, a clean frontend capability gain is the best
near-term breadth increment.

## Decision

Bring the GUI to parity: overlay multiple labelled S-parameter traces in
the egui_plot dB panel with a legend, mirroring the static multi-trace
work. Keep the **series-building logic pure + unit-tested**, separate from
the egui rendering (which is only smoke-tested). Lane = `crates/yee-gui/**`
only; no new dependency (egui_plot 0.35 `Legend`/`Line` suffice).

## Rationale

(1) **Capability gain, not just polish** — the desktop app currently
cannot visualise off-diagonal / multi-port S-parameters at all; this is
the natural completion of the ADR-0063 multi-trace arc (CLI/static → GUI).

(2) **Clean + dispatchable** — additive, single-crate, no new dep, mirrors
an already-shipped extraction idiom; the pure series-builder gives a real
unit-test contract even though egui rendering itself is not unit-testable.

(3) **Breadth** — a different crate (yee-gui / egui_plot live rendering)
from the static-plotters work, and it avoids the deferred quagmires
(mom-002/003 port intrinsic per ADR-0064; FDTD Q6; FEM real-port).

(4) **Honest validation** — a frontend feature, not a §4 solver gate; the
DoD's tested surface is the data-selection/series logic, with the egui
draw layer smoke-tested. No brittle headless-render test is attempted.

## Consequences

* `yee-gui` gains a pure series-building helper (unit-tested) + an
  egui_plot overlay + a small entry-selection control; the single-trace
  default + the Smith / 3D-viewport panels are unchanged.
* No change outside `crates/yee-gui/**`; no new dependency.
* Smith-chart multi-trace + the constant-R/X arc family remain documented
  follow-ons.

## References

* `crates/yee-gui/src/plots.rs` (`:70` single-trace egui_plot dB; Smith
  `:86`), `crates/yee-gui/src/app.rs` (loaded `File` + plot tabs).
* ADR-0063 (the static multi-trace this completes), `crates/yee-cli/src/plot.rs`
  (the entry-extraction idiom), ADR-0064 (why the MoM-accuracy tracks are
  closed, motivating the breadth rotation).
* Spec + plan (2026-05-24).
