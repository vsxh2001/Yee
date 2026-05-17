# ADR-0018: yee-gui ships a Validation panel that runs the aggregator on a worker thread

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

ADR-0008 established `yee-validation` as the canonical aggregator
emitting JSON case records and PNG plots; ADR-0016 then wired it
through Python so notebooks could ingest the same data. Until Phase
1.gui.6, `yee-gui` — the egui desktop shell — had **no in-app
visibility into validation at all**: a user running the GUI had to
shell out to `yee validate` from a separate terminal, then re-open
the GUI to look at the resulting PNGs as static images. That breaks
the feedback loop the GUI exists to support.

Phase 1.gui.6 fixes this with a **Validation panel** inside the
existing egui_dock layout. The design constraints:

- The aggregator wall time is ~10 minutes (CLAUDE.md §4;
  `mom-001` alone is 7–8 min in release). **The UI cannot block
  on this run.** An egui application that stops repainting for
  10 minutes is indistinguishable from a hung process; the OS
  will offer to kill it, and the user will accept.
- The current `yee-gui` does not depend on `tokio` or any other
  async runtime. ADR-0004 (egui pinning) and ADR-0011 (the
  toolchain bump to Rust 1.92 / egui 0.34 / wgpu 29) keep the
  dependency surface minimal on purpose; pulling in `tokio`
  solely to dispatch one long-running validation job is
  disproportionate.
- The aggregator is **already** structured as a synchronous
  `Report::run_all` call that returns when done; `yee-py`
  (ADR-0016) wraps this exactly the same way. The natural shape
  is "spawn it on a worker thread, poll for completion." This
  matches the standard egui pattern from the egui demo app's
  long-running-job example.
- The user-visible feedback needed is: which cases passed,
  which failed, which were skipped, and how long each took. Per-
  case PNG drill-down (clicking a Fail row to open the PNG)
  is *nice* but not core to the feedback loop — knowing whether
  `mom-001` regressed is the question; staring at the plot is
  the follow-up.

Two structural responses were considered:

1. **Add `tokio` to `yee-gui` and run the aggregator on a
   `tokio::task::spawn_blocking` task, with the egui event loop
   awoken via a `tokio::sync::mpsc` channel.** Rejected. The
   payoff is identical to plain `std::thread::spawn` for a
   single long-running job, and the cost is an entire async
   runtime in the GUI's dependency tree. Phase 1.gui.6 should
   not be the first crate to take that hit.
2. **`std::thread::spawn` + `std::sync::mpsc` + egui's
   `request_repaint_after`.** Accepted. See decision below.

## Decision

`yee-gui` Phase 1.gui.6 adds a **Validation** tab to the existing
egui_dock layout (introduced in earlier Phase 1.gui.x work). The
tab contains:

- A **"Run validation"** button. While idle, clicking it spawns
  the aggregator. While a run is in flight, the button is
  disabled and replaced with a spinner + "Running… (N / M cases
  done)" status line.
- A **results table** with sortable columns: **Case**,
  **Status**, **Wall time (s)**. Rows are colored by status:
  green for Pass, red for Fail, yellow for Skipped. The table
  is populated incrementally as case results arrive on the
  channel — a user watching the panel sees `mom-001` show up
  green after ~8 minutes, then the remaining cases tick through
  in seconds.

**Concurrency model.**

- Clicking "Run validation" spawns the aggregator on
  `std::thread::spawn`. The worker thread calls
  `yee_validation::Report::run_all` and emits each
  `ValidationCase` on a `std::sync::mpsc::Sender<CaseEvent>` as
  it completes. A final `CaseEvent::Done(Report)` is sent when
  the aggregator returns.
- The UI thread holds the `Receiver<CaseEvent>`. On each frame,
  it drains the receiver (non-blocking `try_recv` loop),
  appending each new case to its results table.
- To avoid pegging the UI at 60 fps while waiting for slow
  cases (the GUI does not need to redraw 600 times during an
  8-minute `mom-001` wait), the panel calls
  `ctx.request_repaint_after(Duration::from_millis(250))`. This
  schedules a wake-up 250 ms in the future regardless of input;
  combined with the channel drain, it gives a responsive UI at
  a tiny CPU cost.
- **No async runtime is added.** No `tokio`, no `async-std`,
  no `futures`. Pure stdlib `std::thread` + `std::sync::mpsc`.

**Per-case PNG drill-down is explicitly out of scope** for Phase
1.gui.6. Double-clicking a Fail row does not open the plot. That
behaviour is **Phase 1.gui.7**, which adds the
plot-viewer panel and the row → PNG path wiring; doing it inside
1.gui.6 would conflate "make the run visible" with "make the
run interactive."

## Consequences

**What becomes easier:**

- **One-click validation invocation** from the running GUI: the
  user does not have to shell out, does not have to remember the
  CLI subcommand, and does not lose their place in the GUI to
  go check whether the latest geometry change still passes the
  benchmarks.
- **The UI does not freeze for 10 minutes.** The worker-thread
  + mpsc pattern keeps egui drawing at its normal rate, and
  `request_repaint_after(250 ms)` keeps the spinner /
  case-counter alive during the long idle stretches without
  burning CPU.
- **No new dependency surface.** `tokio` is not pulled in.
  ADR-0011's careful dependency pinning (Rust 1.92, egui 0.34,
  wgpu 29) is not disturbed. The Validation panel is purely
  additive in `yee-gui` and depends only on `yee-validation`
  (which the workspace already builds).
- **The panel is incremental.** A user can watch
  `cpml_reflection` and `tfsf_slab` come back Pass within the
  first few seconds and start forming a hypothesis about a
  failing run while `mom-001` is still grinding away.

**What becomes harder:**

- **No per-case plot drill-down today.** A user looking at a
  red `mom-002` row in the table cannot click through to the
  PNG; they have to find the file on disk under the report's
  output directory. **Phase 1.gui.7** closes that gap with a
  plot-viewer pane and a double-click → open handler on the row.
- A second concurrent validation run is intentionally blocked
  (the "Run validation" button is disabled while a run is in
  flight). Allowing two parallel `Report::run_all` invocations
  would require either two output directories or a coordinated
  lock, and neither is justifiable for a 10-minute job that the
  user almost certainly does not want to start twice anyway.
- If the user closes the GUI while a validation run is in
  flight, the worker thread is detached and continues until the
  aggregator's natural completion (the OS reclaims the process
  on exit). The thread does not check a shutdown flag because
  `Report::run_all` is itself uninterruptible; designing
  cooperative cancellation here would push back into
  `yee-validation` and is out of scope.

**What's now closed off:**

- An async runtime (`tokio`, `async-std`) in `yee-gui`. The
  one-job-per-click pattern does not need it, and ADR-0011's
  dependency minimalism is a deliberate constraint.
- Bundling per-case plot drill-down into Phase 1.gui.6. It
  has its own phase number specifically so the brief stays
  small and reviewable.

## References

- `crates/yee-gui/src/panels/validation.rs` — the Validation
  tab: button, status line, sortable results table, worker-
  thread spawn + mpsc drain.
- `crates/yee-gui/src/panels/mod.rs` — egui_dock layout entry
  registering the Validation tab.
- ADR-0008 — validation aggregator JSON + PNG contract; this
  panel is a GUI consumer of that contract.
- ADR-0016 — `yee-py` wraps the same `Report::run_all` call;
  the GUI panel and the Python binding sit on opposite sides
  of the same Rust function.
- ADR-0011 — toolchain bump to Rust 1.92 / egui 0.34 / wgpu
  29; the egui APIs used here (`request_repaint_after`,
  `egui_dock 0.19`) are the post-bump versions.
- ADR-0004 — egui pinning rationale; the no-async-runtime
  position taken here is in the same spirit.
- CLAUDE.md §4 — `mom-001` wall-time (~7–8 min in release)
  is the load-bearing fact behind the worker-thread choice.
- Phase 1.gui.7 (queued) — per-case plot drill-down: row
  double-click → open PNG in a plot-viewer pane.
