//! Application state and top-level UI for the Yee studio shell.
//!
//! Phase 1.gui.0 commit-1 is intentionally a "blank window" — just enough
//! eframe scaffolding to prove the workspace crate builds, runs, and hosts
//! an `egui::Context`. The following commits in this phase add the S11 /
//! Smith plotting tabs (via `egui_dock`) and the `--file` CLI flag.

/// Top-level application state. Empty for the bare skeleton commit; future
/// commits populate this with the loaded Touchstone file and a `DockState`.
pub struct YeeApp {}

impl YeeApp {
    /// Build a fresh app.
    pub fn new() -> Self {
        Self {}
    }
}

impl Default for YeeApp {
    fn default() -> Self {
        Self::new()
    }
}

impl eframe::App for YeeApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::top("menu_bar").show(ctx, |ui| {
            egui::MenuBar::new().ui(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Quit").clicked() {
                        ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                    }
                });
            });
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.centered_and_justified(|ui| {
                ui.label("Yee Studio — Phase 1.gui.0 skeleton");
            });
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_app_constructs() {
        let _ = YeeApp::new();
    }
}
