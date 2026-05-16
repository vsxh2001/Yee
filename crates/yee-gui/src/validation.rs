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
//! Once a [`Report`] is in hand it is rendered as a four-column
//! striped grid (`Case`, `Status`, `Wall time (s)`, `Notes`). Status
//! cells are colour-coded green / red / grey for Passed / Failed /
//! Skipped respectively. Clicking the header of either of the
//! sortable columns (`Case` or `Wall time (s)`) cycles the row order
//! between ascending and descending; clicking a third time clears the
//! sort and falls back to registration order.
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
    sort: SortState,
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

/// Column-and-direction sort key for the result table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
enum SortState {
    /// Registration order from the aggregator.
    #[default]
    None,
    /// Sort by `id` ascending.
    IdAsc,
    /// Sort by `id` descending.
    IdDesc,
    /// Sort by `wall_time_seconds` ascending.
    WallAsc,
    /// Sort by `wall_time_seconds` descending.
    WallDesc,
}

impl Default for ValidationPanel {
    fn default() -> Self {
        Self {
            state: ValidationState::Idle,
            sort: SortState::default(),
        }
    }
}

impl ValidationPanel {
    /// Render the validation panel into `ui`.
    ///
    /// Safe to call every frame; idempotent when nothing has been
    /// triggered. Starts a background thread on button click and
    /// transitions through [`ValidationState`] as the worker
    /// progresses.
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

        // Render the result table once we have one.
        if let ValidationState::Done(report) = &self.state {
            ui.separator();
            ui.horizontal(|ui| {
                ui.label(format!("Total cases: {}", report.cases.len()));
                ui.label(format!("Has failures: {}", report.has_failures()));
                if !matches!(self.sort, SortState::None)
                    && ui.button("Reset sort").clicked()
                {
                    self.sort = SortState::None;
                }
            });

            // Build an index permutation so we render a view over
            // `report.cases` without copying CaseResult values.
            let mut order: Vec<usize> = (0..report.cases.len()).collect();
            match self.sort {
                SortState::None => {}
                SortState::IdAsc => order.sort_by(|&a, &b| {
                    report.cases[a].id.cmp(&report.cases[b].id)
                }),
                SortState::IdDesc => order.sort_by(|&a, &b| {
                    report.cases[b].id.cmp(&report.cases[a].id)
                }),
                SortState::WallAsc => order.sort_by(|&a, &b| {
                    report.cases[a]
                        .wall_time_seconds
                        .partial_cmp(&report.cases[b].wall_time_seconds)
                        .unwrap_or(std::cmp::Ordering::Equal)
                }),
                SortState::WallDesc => order.sort_by(|&a, &b| {
                    report.cases[b]
                        .wall_time_seconds
                        .partial_cmp(&report.cases[a].wall_time_seconds)
                        .unwrap_or(std::cmp::Ordering::Equal)
                }),
            }

            egui::ScrollArea::vertical().show(ui, |ui| {
                egui::Grid::new("validation_grid")
                    .num_columns(4)
                    .striped(true)
                    .show(ui, |ui| {
                        // Header row. `Case` and `Wall time (s)` are
                        // click-to-sort.
                        if ui
                            .button(sort_header_label("Case", self.sort, SortDim::Id))
                            .clicked()
                        {
                            self.sort = next_sort(self.sort, SortDim::Id);
                        }
                        ui.strong("Status");
                        if ui
                            .button(sort_header_label(
                                "Wall time (s)",
                                self.sort,
                                SortDim::Wall,
                            ))
                            .clicked()
                        {
                            self.sort = next_sort(self.sort, SortDim::Wall);
                        }
                        ui.strong("Notes");
                        ui.end_row();

                        for &idx in &order {
                            let case = &report.cases[idx];
                            ui.label(&case.id);
                            let (text, color) = status_style(case.status);
                            ui.colored_label(color, text);
                            ui.label(format!("{:.3}", case.wall_time_seconds));
                            // Notes can be long — wrap into the cell.
                            ui.add(egui::Label::new(&case.notes).wrap());
                            ui.end_row();
                        }
                    });
            });
        }
    }
}

/// Which column a sort transition applies to.
#[derive(Debug, Clone, Copy)]
enum SortDim {
    Id,
    Wall,
}

/// Render a header button label including a small ASCII arrow to
/// indicate the current sort direction.
fn sort_header_label(base: &str, sort: SortState, dim: SortDim) -> String {
    let arrow = match (sort, dim) {
        (SortState::IdAsc, SortDim::Id) => " ▲",
        (SortState::IdDesc, SortDim::Id) => " ▼",
        (SortState::WallAsc, SortDim::Wall) => " ▲",
        (SortState::WallDesc, SortDim::Wall) => " ▼",
        _ => "",
    };
    format!("{base}{arrow}")
}

/// Cycle the sort state for a given dimension: `None → Asc → Desc → None`.
/// Clicking a different column resets to `Asc` on that column.
fn next_sort(current: SortState, dim: SortDim) -> SortState {
    match (current, dim) {
        (SortState::IdAsc, SortDim::Id) => SortState::IdDesc,
        (SortState::IdDesc, SortDim::Id) => SortState::None,
        (SortState::WallAsc, SortDim::Wall) => SortState::WallDesc,
        (SortState::WallDesc, SortDim::Wall) => SortState::None,
        (_, SortDim::Id) => SortState::IdAsc,
        (_, SortDim::Wall) => SortState::WallAsc,
    }
}

/// Map a [`Status`] to its display string + accent colour.
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
        assert_eq!(panel.sort, SortState::None);
    }

    #[test]
    fn id_sort_cycle() {
        let mut s = SortState::None;
        s = next_sort(s, SortDim::Id);
        assert_eq!(s, SortState::IdAsc);
        s = next_sort(s, SortDim::Id);
        assert_eq!(s, SortState::IdDesc);
        s = next_sort(s, SortDim::Id);
        assert_eq!(s, SortState::None);
    }

    #[test]
    fn wall_sort_cycle() {
        let mut s = SortState::None;
        s = next_sort(s, SortDim::Wall);
        assert_eq!(s, SortState::WallAsc);
        s = next_sort(s, SortDim::Wall);
        assert_eq!(s, SortState::WallDesc);
        s = next_sort(s, SortDim::Wall);
        assert_eq!(s, SortState::None);
    }

    #[test]
    fn switching_column_resets_to_asc() {
        let s = SortState::IdDesc;
        assert_eq!(next_sort(s, SortDim::Wall), SortState::WallAsc);
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
