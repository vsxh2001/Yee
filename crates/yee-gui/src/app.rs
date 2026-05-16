//! Application state, dock layout, and top-level UI for the Yee studio shell.
//!
//! The shell hosts two tabs inside an `egui_dock::DockArea`:
//!
//! - `S11Db` — `20·log10|S11|` vs frequency
//! - `Smith` — `S11` trajectory on a Smith-chart canvas (unit circle reference)
//!
//! A left side panel exposes loaded-file metadata; the menu bar provides
//! `File → Quit`. File opening is driven by a `--file` CLI flag at startup
//! (Phase 1.gui.0 keeps the GUI free of `rfd`-based pickers — see README).

use crate::plots::{show_s11_db_plot, show_smith_chart};
use egui_dock::{DockArea, DockState, NodeIndex, Style};
use yee_io::touchstone::{self, File as TsFile};

/// Tabs hosted in the central dock area.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TabKind {
    /// `20·log10|S11|` line plot vs frequency.
    S11Db,
    /// Smith-chart visualisation of `S11` in the complex plane.
    Smith,
}

impl TabKind {
    fn title(self) -> &'static str {
        match self {
            TabKind::S11Db => "S11 magnitude (dB)",
            TabKind::Smith => "Smith chart",
        }
    }
}

/// Top-level application state.
pub struct YeeApp {
    /// The currently loaded Touchstone file, if any.
    file: Option<TsFile>,
    /// Most recent load error, surfaced as a status banner.
    load_error: Option<String>,
    /// Tab layout managed by `egui_dock`.
    dock: DockState<TabKind>,
}

impl YeeApp {
    /// Build a fresh app, optionally pre-loading a Touchstone file from `path`.
    pub fn new(initial_file: Option<std::path::PathBuf>) -> Self {
        // Default layout: two tabs side-by-side in one node.
        let mut dock = DockState::new(vec![TabKind::S11Db]);
        let surface = dock.main_surface_mut();
        surface.split_right(NodeIndex::root(), 0.5, vec![TabKind::Smith]);

        let mut app = Self {
            file: None,
            load_error: None,
            dock,
        };
        if let Some(path) = initial_file {
            app.load_file(&path);
        }
        app
    }

    fn load_file(&mut self, path: &std::path::Path) {
        match touchstone::read(path) {
            Ok(f) => {
                tracing::info!(
                    "loaded {} ({} ports, {} samples)",
                    path.display(),
                    f.n_ports,
                    f.freq_hz.len()
                );
                self.file = Some(f);
                self.load_error = None;
            }
            Err(e) => {
                let msg = format!("{}: {e}", path.display());
                tracing::error!("{msg}");
                self.file = None;
                self.load_error = Some(msg);
            }
        }
    }

    /// Render the left-hand metadata panel.
    fn metadata_panel(&self, ui: &mut egui::Ui) {
        ui.heading("File metadata");
        ui.separator();
        match &self.file {
            None => {
                ui.label("No file loaded.");
                ui.label("");
                ui.label("Open a .s1p file to begin:");
                ui.code("cargo run -p yee-gui --release -- \\\n    --file path/to/dipole.s1p");
            }
            Some(f) => {
                ui.label(format!("n_ports: {}", f.n_ports));
                ui.label(format!("z0 (Ω): {}", f.z0));
                ui.label(format!("format: {:?}", f.format));
                ui.label(format!("freq unit (on disk): {:?}", f.freq_unit));
                ui.label(format!("samples: {}", f.freq_hz.len()));
                if let (Some(&lo), Some(&hi)) = (f.freq_hz.first(), f.freq_hz.last()) {
                    ui.label(format!("f range: {:.3} → {:.3} GHz", lo * 1e-9, hi * 1e-9));
                }
                if !f.comments.is_empty() {
                    ui.separator();
                    ui.label("Comments:");
                    for c in &f.comments {
                        ui.label(format!("! {}", c.trim_end()));
                    }
                }
            }
        }
        if let Some(err) = &self.load_error {
            ui.separator();
            ui.colored_label(egui::Color32::LIGHT_RED, "Load error:");
            ui.label(err);
        }
    }
}

impl eframe::App for YeeApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Menu bar.
        egui::TopBottomPanel::top("menu_bar").show(ctx, |ui| {
            egui::MenuBar::new().ui(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Open .s1p… (use --file CLI flag)").clicked() {
                        // Phase 1.gui.0 keeps the file picker out of scope; the
                        // menu entry is still surfaced so the workflow is
                        // discoverable.
                        ui.close();
                    }
                    ui.separator();
                    if ui.button("Quit").clicked() {
                        ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                    }
                });
            });
        });

        // Metadata side panel.
        egui::SidePanel::left("metadata").show(ctx, |ui| {
            self.metadata_panel(ui);
        });

        // Central dock area with the two plot tabs.
        egui::CentralPanel::default().show(ctx, |ui| {
            let mut viewer = TabViewer {
                file: self.file.as_ref(),
            };
            DockArea::new(&mut self.dock)
                .style(Style::from_egui(ui.style().as_ref()))
                .show_inside(ui, &mut viewer);
        });
    }
}

/// `egui_dock` tab viewer. Borrows the loaded file (if any) so plots can
/// render either real data or a placeholder.
struct TabViewer<'a> {
    file: Option<&'a TsFile>,
}

impl<'a> egui_dock::TabViewer for TabViewer<'a> {
    type Tab = TabKind;

    fn title(&mut self, tab: &mut Self::Tab) -> egui::WidgetText {
        tab.title().into()
    }

    fn ui(&mut self, ui: &mut egui::Ui, tab: &mut Self::Tab) {
        match self.file {
            None => {
                ui.centered_and_justified(|ui| {
                    ui.label("Open a .s1p file to begin.");
                });
            }
            Some(f) => {
                // S11 lives at row-major slot 0 for any port count.
                let s11: Vec<num_complex::Complex64> = f.data.iter().map(|m| m[0]).collect();
                match tab {
                    TabKind::S11Db => show_s11_db_plot(ui, &f.freq_hz, &s11),
                    TabKind::Smith => show_smith_chart(ui, &s11),
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_app_without_file_has_no_load_error() {
        let app = YeeApp::new(None);
        assert!(app.file.is_none());
        assert!(app.load_error.is_none());
    }

    #[test]
    fn loading_missing_file_records_error() {
        let app = YeeApp::new(Some(std::path::PathBuf::from(
            "/nonexistent/path/to/missing.s1p",
        )));
        assert!(app.file.is_none());
        assert!(
            app.load_error.is_some(),
            "expected load_error to be set for a missing file"
        );
    }
}
