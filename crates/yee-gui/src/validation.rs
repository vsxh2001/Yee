//! Validation-aggregator panel hosted in the yee-gui dock.
//!
//! Phase 1.gui.validation.0 skeleton: this lands the `ValidationPanel`
//! type with its idle-state UI (a single button placeholder). The
//! background-thread runner and the result table are added in
//! follow-up commits in this series.
//!
//! Lane: `crates/yee-gui/**` only. The aggregator is invoked through
//! its public API (`yee_validation::Report::run_all`); no
//! validation-side changes are made here.

use eframe::egui;
use yee_validation::{Report, Status};

/// Stateful UI element rendering the validation aggregator panel.
///
/// Construct via [`ValidationPanel::default`] and forward each frame
/// to [`ValidationPanel::ui`] from inside the dock tab viewer.
pub struct ValidationPanel {
    state: ValidationState,
}

/// Internal state machine for the panel.
enum ValidationState {
    /// No aggregator run has been triggered yet.
    Idle,
    /// A previous run has completed; the report is cached for
    /// rendering until the user re-runs. Populated in commit 2 of
    /// this series (background-thread runner).
    #[allow(dead_code)]
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
    /// In this skeleton commit the button is wired but does nothing
    /// when clicked; the runner lands in the next commit. Calling
    /// [`status_style`] keeps the import live so the next commit's
    /// patch is a pure addition.
    pub fn ui(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            let label = match &self.state {
                ValidationState::Idle => "Run validation (~10 min)",
                ValidationState::Done(_) => "Re-run validation",
            };
            let _clicked = ui.button(label).clicked();
            ui.label("Idle — click to run the full aggregator.");
        });

        // Render placeholder table header so the panel is not empty
        // before any run. The full table arrives in commit 3 of this
        // series.
        if let ValidationState::Done(report) = &self.state {
            ui.separator();
            ui.label(format!("Total cases: {}", report.cases.len()));
            let _ = status_style(Status::Passed);
        }
    }
}

/// Map a [`Status`] to its display string + accent colour. Returned
/// colour is used by the table view in commit 3 of this series.
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
