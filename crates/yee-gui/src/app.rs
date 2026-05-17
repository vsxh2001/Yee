//! Application state, dock layout, and top-level UI for the Yee studio shell.
//!
//! The shell hosts four tabs inside an `egui_dock::DockArea`:
//!
//! - `S11Db`       — `20·log10|S11|` vs frequency
//! - `Smith`       — `S11` trajectory on a Smith-chart canvas (unit circle reference)
//! - `Mesh3D`      — wgpu-rendered 3D triangle mesh (Phase 1.gui.1)
//! - `Validation`  — yee-validation aggregator runner + sortable result table
//!
//! A left side panel exposes loaded-file metadata + the 3D-viewport controls
//! (wireframe toggle, camera readout). The menu bar provides `File → Quit`.
//! File opening is driven by a `--file` CLI flag at startup (the GUI keeps
//! `rfd`-based pickers out of scope through Phase 1.gui.1 — see README).

use crate::plots::{show_s11_db_plot, show_smith_chart};
use crate::validation::ValidationPanel;
use crate::viewport::{MeshCallback, ViewportState, thin_cylinder};
use egui_dock::{DockArea, DockState, NodeIndex, Style};
use yee_io::touchstone::{self, File as TsFile};

/// Tabs hosted in the central dock area.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TabKind {
    /// `20·log10|S11|` line plot vs frequency.
    S11Db,
    /// Smith-chart visualisation of `S11` in the complex plane.
    Smith,
    /// Wgpu-rendered 3D triangle mesh viewport.
    Mesh3D,
    /// yee-validation aggregator runner + sortable result table.
    Validation,
}

impl TabKind {
    fn title(self) -> &'static str {
        match self {
            TabKind::S11Db => "S11 magnitude (dB)",
            TabKind::Smith => "Smith chart",
            TabKind::Mesh3D => "Mesh 3D",
            TabKind::Validation => "Validation",
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
    /// 3D viewport state (orbit camera + active mesh + wireframe toggle).
    viewport_state: ViewportState,
    /// Validation aggregator panel (background runner + sortable table).
    validation_panel: ValidationPanel,
}

impl YeeApp {
    /// Build a fresh app, optionally pre-loading a Touchstone file from `path`.
    ///
    /// The dock layout starts with the two plot tabs at the top of the central
    /// region and the Mesh 3D tab split off below. This keeps the Phase 1.gui.0
    /// plot workflow visually identical while making the new viewport
    /// immediately discoverable. The Validation tab is stacked on top of the
    /// Mesh 3D tab so it shares the lower pane without consuming an extra
    /// split — switching between them is a single click.
    pub fn new(initial_file: Option<std::path::PathBuf>) -> Self {
        let mut dock = DockState::new(vec![TabKind::S11Db]);
        {
            let surface = dock.main_surface_mut();
            // Plot pair lives in the root node.
            surface.split_right(NodeIndex::root(), 0.5, vec![TabKind::Smith]);
            // Drop the 3D viewport and the Validation panel below the plots
            // so all four tabs are visible at once in the default layout
            // (Mesh 3D / Validation stack as siblings).
            surface.split_below(
                NodeIndex::root(),
                0.55,
                vec![TabKind::Mesh3D, TabKind::Validation],
            );
        }

        // Walking-skeleton mesh: a 0.5 m × 5 mm thin cylinder, tessellated
        // finely enough for the lighting to read as smooth.
        let mesh = thin_cylinder(0.5, 0.005, 16, 32);
        let viewport_state = ViewportState::new(mesh);

        let mut app = Self {
            file: None,
            load_error: None,
            dock,
            viewport_state,
            validation_panel: ValidationPanel::default(),
        };
        if let Some(path) = initial_file {
            app.load_touchstone(&path);
        }
        app
    }

    /// Parse a Touchstone file from `path` via [`yee_io::touchstone::read`]
    /// and store the result on the app.
    ///
    /// On success the loaded file replaces any previously held one and
    /// `load_error` is cleared; the existing S11 and Smith-chart tabs read
    /// from `self.file` on every frame, so a repaint will pick up the new
    /// data automatically. On failure the file slot is cleared and the
    /// error message is surfaced in the side panel as a red banner.
    pub fn load_touchstone(&mut self, path: &std::path::Path) {
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

    /// Render the left-hand metadata + viewport-controls panel.
    fn metadata_panel(&mut self, ui: &mut egui::Ui) {
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

        ui.separator();
        ui.heading("Viewport");
        ui.checkbox(&mut self.viewport_state.wireframe, "Wireframe overlay");
        ui.label(format!(
            "yaw / pitch: {:>5.1}° / {:>5.1}°",
            self.viewport_state.camera_yaw_deg, self.viewport_state.camera_pitch_deg
        ));
        ui.label(format!(
            "distance: {:.3} m",
            self.viewport_state.camera_dist
        ));
    }
}

impl eframe::App for YeeApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();

        // Menu bar.
        egui::Panel::top("menu_bar").show_inside(ui, |ui| {
            egui::MenuBar::new().ui(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Open .s1p… (use --file CLI flag)").clicked() {
                        // The file picker stays out of scope through Phase
                        // 1.gui.1; the menu entry is surfaced so the workflow
                        // is discoverable.
                        ui.close();
                    }
                    ui.separator();
                    if ui.button("Quit").clicked() {
                        ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                    }
                });
            });
        });

        // Metadata + viewport-controls side panel.
        egui::Panel::left("metadata").show_inside(ui, |ui| {
            self.metadata_panel(ui);
        });

        // Central dock area with the four tabs.
        egui::CentralPanel::default().show_inside(ui, |ui| {
            let mut viewer = TabViewer {
                file: self.file.as_ref(),
                viewport_state: &mut self.viewport_state,
                validation_panel: &mut self.validation_panel,
            };
            DockArea::new(&mut self.dock)
                .style(Style::from_egui(ui.style().as_ref()))
                .show_inside(ui, &mut viewer);
        });
    }
}

/// `egui_dock` tab viewer. Borrows the loaded file (for the plot tabs), the
/// viewport state (for the Mesh 3D tab) and the validation panel (for the
/// Validation tab) so each tab can render its content.
struct TabViewer<'a> {
    file: Option<&'a TsFile>,
    viewport_state: &'a mut ViewportState,
    validation_panel: &'a mut ValidationPanel,
}

impl<'a> egui_dock::TabViewer for TabViewer<'a> {
    type Tab = TabKind;

    fn title(&mut self, tab: &mut Self::Tab) -> egui::WidgetText {
        tab.title().into()
    }

    fn ui(&mut self, ui: &mut egui::Ui, tab: &mut Self::Tab) {
        match tab {
            TabKind::S11Db | TabKind::Smith => match self.file {
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
                        TabKind::Mesh3D | TabKind::Validation => unreachable!(),
                    }
                }
            },
            TabKind::Mesh3D => {
                show_mesh_viewport(ui, self.viewport_state);
            }
            TabKind::Validation => {
                self.validation_panel.ui(ui);
            }
        }
    }
}

/// Render the wgpu-backed 3D mesh viewport into `ui`.
///
/// The function reserves the full available rect, wires drag-to-orbit and
/// scroll-to-zoom interactions, then schedules an [`egui_wgpu`] paint
/// callback that draws the mesh into that rect.
///
/// Interaction:
/// - Drag (any mouse button): orbit camera (0.5° per pixel; pitch clamped to
///   ±89° to avoid look-at-up degeneracy).
/// - Mouse wheel while hovered: exponential zoom, clamped to `[0.01, 1e4]`
///   metres so the camera can't slip through the origin or escape to
///   infinity.
fn show_mesh_viewport(ui: &mut egui::Ui, state: &mut ViewportState) {
    let available = ui.available_size_before_wrap();
    let (rect, response) = ui.allocate_exact_size(available, egui::Sense::click_and_drag());

    // Orbit interaction: drag to rotate.
    if response.dragged() {
        let drag = response.drag_delta();
        state.camera_yaw_deg += drag.x * 0.5;
        state.camera_pitch_deg = (state.camera_pitch_deg - drag.y * 0.5).clamp(-89.0, 89.0);
    }

    // Zoom: scroll wheel.
    if response.hovered() {
        let scroll = ui.input(|i| i.smooth_scroll_delta.y);
        if scroll != 0.0 {
            // Exponential zoom so the perceived rate is constant in
            // log-distance space.
            let factor = (-scroll * 0.005).exp();
            state.camera_dist = (state.camera_dist * factor).clamp(0.01, 1.0e4);
        }
    }

    // Build the MVP for this frame.
    let aspect = if rect.height() > 0.0 {
        rect.width() / rect.height()
    } else {
        1.0
    };
    let mvp = state.view_proj(aspect);
    let camera_pos = state.camera_position();

    let callback = egui_wgpu::Callback::new_paint_callback(
        rect,
        MeshCallback {
            mesh: state.mesh.clone(),
            mvp,
            camera_pos,
            wireframe: state.wireframe,
        },
    );
    ui.painter().add(callback);
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

    #[test]
    fn new_app_has_default_viewport_state() {
        let app = YeeApp::new(None);
        // Default mesh is non-empty.
        assert!(!app.viewport_state.mesh.vertices.is_empty());
        assert!(!app.viewport_state.mesh.indices.is_empty());
        // Camera distance follows the bbox-based default.
        assert!(app.viewport_state.camera_dist > 0.0);
        assert!(!app.viewport_state.wireframe);
    }
}
