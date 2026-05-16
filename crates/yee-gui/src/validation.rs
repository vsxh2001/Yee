//! Validation-aggregator panel hosted in the yee-gui dock.
//!
//! The panel exposes a single button that dispatches
//! [`yee_validation::Report::run_all`] onto a background thread and
//! polls for completion via a `std::sync::mpsc` channel. While the
//! aggregator is running the UI thread stays responsive — the button
//! disables itself, the status text flips to "Running…", and the
//! frame is invalidated every ~250 ms so completion is picked up
//! promptly without spin-polling.
//!
//! The result-table view lands in the next commit in this series;
//! this commit limits itself to the runner + state machine so the
//! diff is small and reviewable.
//!
//! Lane: `crates/yee-gui/**` only. The aggregator is invoked through
//! its public API; no validation-side changes are made here.

use eframe::egui;
use std::sync::mpsc::{Receiver, channel};
use std::thread;
use yee_validation::{Report, Status};

/// Stateful UI element rendering the validation aggregator panel.
///
/// Construct via [`ValidationPanel::default`] and forward each frame
/// to [`ValidationPanel::ui`] from inside the dock tab viewer.
pub struct ValidationPanel {
    state: ValidationState,
}

/// Internal state machine for the panel.
///
/// The transitions are linear: `Idle → Running → Done`, with the
/// `Done` → `Running` edge reachable by re-clicking the run button.
enum ValidationState {
    /// No aggregator run has been triggered yet (or the panel was
    /// just constructed).
    Idle,
    /// A background thread is executing [`Report::run_all`]; the
    /// channel receiver will deliver the final report.
    Running {
        /// Receiver half of the cross-thread completion channel.
        rx: Receiver<Report>,
    },
    /// A previous run has completed; the report is cached for
    /// rendering until the user re-runs.
    Done(Report),
}

impl Default for ValidationPanel {
    fn default() -> Self {
        Self {
            state: ValidationState::Idle,
        }
    }
}

impl ValidationPanel {
    /// Render the validation panel into `ui`.
    ///
    /// Safe to call every frame; idempotent when nothing has been
    /// triggered. Starts a background thread on button click and
    /// transitions through [`ValidationState`] as the worker
    /// progresses. The full result table is rendered by the
    /// subsequent commit in this series; for now we just confirm a
    /// report was received.
    pub fn ui(&mut self, ui: &mut egui::Ui) {
        // Top row: run/re-run button + live status text.
        ui.horizontal(|ui| {
            let running = matches!(self.state, ValidationState::Running { .. });
            let label = match &self.state {
                ValidationState::Idle => "Run validation (~10 min)",
                ValidationState::Running { .. } => "Running…",
                ValidationState::Done(_) => "Re-run validation",
            };
            if ui.add_enabled(!running, egui::Button::new(label)).clicked() {
                let (tx, rx) = channel();
                thread::spawn(move || {
                    let report = Report::run_all();
                    // Receiver may have been dropped if the user closed
                    // the panel; that's fine, the report is just
                    // discarded.
                    let _ = tx.send(report);
                });
                self.state = ValidationState::Running { rx };
            }

            // Lightweight status indicator alongside the button.
            match &self.state {
                ValidationState::Idle => {
                    ui.label("Idle — click to run the full aggregator.");
                }
                ValidationState::Running { .. } => {
                    ui.spinner();
                    ui.label("Aggregator running on a background thread.");
                }
                ValidationState::Done(report) => {
                    let n = report.cases.len();
                    let failed = report
                        .cases
                        .iter()
                        .filter(|c| c.status == Status::Failed)
                        .count();
                    ui.label(format!("Done — {n} cases, {failed} failed."));
                }
            }
        });

        // Non-blocking poll for the background worker.
        let mut completed: Option<Report> = None;
        if let ValidationState::Running { rx } = &self.state {
            match rx.try_recv() {
                Ok(report) => completed = Some(report),
                Err(_) => {
                    // Keep the UI responsive without spin-polling: ask
                    // egui to repaint in 250 ms.
                    ui.ctx()
                        .request_repaint_after(std::time::Duration::from_millis(250));
                }
            }
        }
        if let Some(report) = completed {
            self.state = ValidationState::Done(report);
        }

        // Minimal result view; the full sortable table arrives in the
        // next commit.
        if let ValidationState::Done(report) = &self.state {
            ui.separator();
            ui.label(format!("Total cases: {}", report.cases.len()));
            ui.label(format!("Has failures: {}", report.has_failures()));
            let _ = status_style(Status::Passed);
        }
    }
}

/// Map a [`Status`] to its display string + accent colour. Returned
/// colour is consumed by the table view in the next commit.
fn status_style(s: Status) -> (&'static str, egui::Color32) {
    match s {
        Status::Passed => ("Passed", egui::Color32::from_rgb(60, 180, 60)),
        Status::Failed => ("Failed", egui::Color32::from_rgb(200, 60, 60)),
        Status::Skipped => ("Skipped", egui::Color32::from_rgb(160, 160, 160)),
    }
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_panel_is_idle() {
        let panel = ValidationPanel::default();
        assert!(matches!(panel.state, ValidationState::Idle));
    }

    #[test]
    fn status_styling_distinct_colors() {
        let (_, pass) = status_style(Status::Passed);
        let (_, fail) = status_style(Status::Failed);
        let (_, skip) = status_style(Status::Skipped);
        assert_ne!(pass, fail);
        assert_ne!(pass, skip);
        assert_ne!(fail, skip);
    }
}
