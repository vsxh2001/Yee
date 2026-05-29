//! `StudioApp` — the `eframe` shell for the Yee Filter Studio (App.0).
//!
//! Hosts three regions over the headless [`StudioState`] logic layer:
//!
//! - a left **spec-editor** [`egui::SidePanel`] (f0, FBW, order, ripple,
//!   return-loss, stopband points, approximation),
//! - a central **synthesis** [`egui::CentralPanel`] (g-values, coupling matrix
//!   grid, external Q, coloured PASS/FAIL verdict + notes), and
//! - an [`egui_plot::Plot`] of `|S21|` (dB) vs frequency (GHz) with each
//!   spec-mask region shaded on its forbidden side, and
//! - a **physical-dimensions** section (F1.2.0; ADR-0099): editable substrate
//!   `ε_r` / `h` (mm) inputs plus the synthesized microstrip line width,
//!   resonator length, and per-section gaps (or the [`yee_filter::DimError`]
//!   string when the coupling is unrealizable on the chosen substrate).
//!
//! Every edit in the side panel — and every substrate edit in the dimensions
//! section — calls [`StudioState::recompute`] so the central panel, the plot,
//! and the dimensions stay live.

use eframe::egui;
use egui::Color32;
use egui_plot::{Legend, Line, Plot, PlotPoints, Polygon};

use crate::{MaskRegionView, StudioState};
use yee_synth::Approximation;

/// Top-level application state: the headless [`StudioState`] plus a little UI
/// scratch (the editable ripple value when the approximation is Chebyshev).
pub struct StudioApp {
    /// The headless logic state (spec + all derived fields).
    pub state: StudioState,
}

impl StudioApp {
    /// Build a [`StudioApp`] wrapping the given [`StudioState`].
    pub fn new(state: StudioState) -> Self {
        Self { state }
    }
}

impl eframe::App for StudioApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let mut dirty = false;

        egui::Panel::left("spec_editor")
            .resizable(true)
            .default_size(300.0)
            .show_inside(ui, |ui| {
                ui.heading("Filter spec");
                ui.separator();

                // f0 in GHz (stored as Hz).
                let mut f0_ghz = self.state.spec.f0_hz * 1.0e-9;
                ui.horizontal(|ui| {
                    ui.label("f0 (GHz)");
                    if ui
                        .add(
                            egui::DragValue::new(&mut f0_ghz)
                                .speed(0.01)
                                .range(0.01..=1000.0),
                        )
                        .changed()
                    {
                        self.state.spec.f0_hz = f0_ghz * 1.0e9;
                        dirty = true;
                    }
                });

                // Fractional bandwidth.
                ui.horizontal(|ui| {
                    ui.label("FBW");
                    if ui
                        .add(
                            egui::DragValue::new(&mut self.state.spec.fbw)
                                .speed(0.001)
                                .range(0.001..=2.0),
                        )
                        .changed()
                    {
                        dirty = true;
                    }
                });

                // Order (explicit).
                let mut order = self.state.spec.order.unwrap_or(5);
                ui.horizontal(|ui| {
                    ui.label("Order N");
                    if ui
                        .add(egui::DragValue::new(&mut order).speed(0.1).range(1..=20))
                        .changed()
                    {
                        self.state.spec.order = Some(order);
                        dirty = true;
                    }
                });

                ui.separator();
                ui.label("Mask");

                // Passband ripple (dB).
                ui.horizontal(|ui| {
                    ui.label("Ripple (dB)");
                    if ui
                        .add(
                            egui::DragValue::new(&mut self.state.spec.mask.passband_ripple_db)
                                .speed(0.01)
                                .range(0.001..=10.0),
                        )
                        .changed()
                    {
                        dirty = true;
                    }
                });

                // Return loss (dB).
                ui.horizontal(|ui| {
                    ui.label("Return loss (dB)");
                    if ui
                        .add(
                            egui::DragValue::new(&mut self.state.spec.mask.return_loss_db)
                                .speed(0.1)
                                .range(0.1..=60.0),
                        )
                        .changed()
                    {
                        dirty = true;
                    }
                });

                ui.separator();
                ui.label("Stopband points (GHz, dB)");
                let mut remove: Option<usize> = None;
                for (i, point) in self.state.spec.mask.stopband.iter_mut().enumerate() {
                    ui.horizontal(|ui| {
                        let mut f_ghz = point.0 * 1.0e-9;
                        if ui
                            .add(
                                egui::DragValue::new(&mut f_ghz)
                                    .speed(0.01)
                                    .range(0.01..=1000.0)
                                    .prefix("f="),
                            )
                            .changed()
                        {
                            point.0 = f_ghz * 1.0e9;
                            dirty = true;
                        }
                        if ui
                            .add(
                                egui::DragValue::new(&mut point.1)
                                    .speed(0.5)
                                    .range(0.0..=120.0)
                                    .prefix("rej="),
                            )
                            .changed()
                        {
                            dirty = true;
                        }
                        if ui.small_button("x").clicked() {
                            remove = Some(i);
                        }
                    });
                }
                if let Some(i) = remove {
                    self.state.spec.mask.stopband.remove(i);
                    dirty = true;
                }
                if ui.button("+ stopband point").clicked() {
                    self.state.spec.mask.stopband.push((2.4e9, 40.0));
                    dirty = true;
                }

                ui.separator();
                ui.label("Approximation");
                // Chebyshev ripple shares the passband-ripple value for a
                // satisfiable default; toggling only swaps the response shape.
                let is_cheby = matches!(
                    self.state.spec.approximation,
                    Approximation::Chebyshev { .. }
                );
                let mut selected_cheby = is_cheby;
                egui::ComboBox::from_id_salt("approx")
                    .selected_text(if is_cheby { "Chebyshev" } else { "Butterworth" })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut selected_cheby, false, "Butterworth");
                        ui.selectable_value(&mut selected_cheby, true, "Chebyshev");
                    });
                if selected_cheby != is_cheby {
                    self.state.spec.approximation = if selected_cheby {
                        Approximation::Chebyshev {
                            ripple_db: self.state.spec.mask.passband_ripple_db,
                        }
                    } else {
                        Approximation::Butterworth
                    };
                    dirty = true;
                }
            });

        if dirty {
            self.state.recompute();
        }

        egui::CentralPanel::default().show_inside(ui, |ui| {
            ui.heading("Yee Filter Studio");

            // ---- mask verdict --------------------------------------------
            ui.horizontal(|ui| {
                ui.label("Mask verdict:");
                if self.state.mask_pass {
                    ui.colored_label(Color32::from_rgb(40, 180, 80), "PASS");
                } else {
                    ui.colored_label(Color32::from_rgb(220, 60, 60), "FAIL");
                }
            });
            for note in &self.state.mask_notes {
                ui.label(note);
            }

            ui.separator();

            // ---- prototype g-values --------------------------------------
            let proto = &self.state.project.prototype;
            ui.label(format!("Prototype (order N={})", proto.order()));
            let g_str: Vec<String> = proto
                .g
                .iter()
                .enumerate()
                .map(|(i, gi)| format!("g{i}={gi:.4}"))
                .collect();
            ui.label(g_str.join("   "));

            // ---- external Q ----------------------------------------------
            let coupling = &self.state.project.coupling;
            ui.label(format!(
                "Qe_in = {:.4}   Qe_out = {:.4}",
                coupling.qe_in, coupling.qe_out
            ));

            // ---- coupling matrix grid ------------------------------------
            ui.separator();
            ui.label("Coupling matrix M (normalized)");
            egui::Grid::new("coupling_matrix")
                .striped(true)
                .show(ui, |ui| {
                    for row in &coupling.m {
                        for v in row {
                            ui.label(format!("{v:+.4}"));
                        }
                        ui.end_row();
                    }
                });

            // ---- |S21| vs spec-mask plot ---------------------------------
            ui.separator();
            ui.label("|S21| (dB) vs spec mask");
            show_response_plot(ui, &self.state);

            // ---- physical dimensions (F1.2.0; ADR-0099) ------------------
            ui.separator();
            show_dimensions(ui, &mut self.state);

            // ---- top-view layout preview (F1.2.0; ADR-0101) --------------
            ui.separator();
            ui.label("Layout (top view, mm)");
            show_layout(ui, &self.state);
        });
    }
}

/// Render the "Physical dimensions" section: editable substrate `ε_r` / `h`
/// (mm) inputs and the synthesized microstrip dimensions (F1.2.0; ADR-0099).
///
/// Editing `ε_r` or `h` updates [`StudioState::eps_r`] / [`StudioState::h_m`]
/// and re-derives via [`StudioState::recompute`] so the read-out stays live.
/// Shows the line width, resonator length, and per-section gaps in millimetres
/// on `Ok`, or the [`yee_filter::DimError`] string on `Err`.
fn show_dimensions(ui: &mut egui::Ui, state: &mut StudioState) {
    /// Metres → millimetres for display.
    const MM: f64 = 1.0e3;

    ui.heading("Physical dimensions");
    ui.label("Edge-coupled microstrip (F1.2.0)");

    let mut dirty = false;

    // Substrate ε_r.
    ui.horizontal(|ui| {
        ui.label("εr");
        if ui
            .add(
                egui::DragValue::new(&mut state.eps_r)
                    .speed(0.01)
                    .range(1.0..=20.0),
            )
            .changed()
        {
            dirty = true;
        }
    });

    // Substrate height h, edited in millimetres (stored in metres).
    ui.horizontal(|ui| {
        ui.label("h (mm)");
        let mut h_mm = state.h_m * MM;
        if ui
            .add(
                egui::DragValue::new(&mut h_mm)
                    .speed(0.01)
                    .range(0.01..=20.0),
            )
            .changed()
        {
            state.h_m = h_mm / MM;
            dirty = true;
        }
    });

    if dirty {
        state.recompute();
    }

    // Read-out (or the error string when the coupling is unrealizable).
    match &state.dims {
        Ok(dims) => {
            ui.label(format!("line width: {:.4} mm", dims.line_width_m * MM));
            ui.label(format!(
                "resonator length: {:.4} mm",
                dims.resonator_length_m * MM
            ));
            for (i, &gap) in dims.gaps_m.iter().enumerate() {
                ui.label(format!("gap {}–{}: {:.4} mm", i + 1, i + 2, gap * MM));
            }
        }
        Err(msg) => {
            ui.colored_label(
                Color32::from_rgb(220, 60, 60),
                format!("no dimensions: {msg}"),
            );
        }
    }
}

/// Render the "Layout" top-view canvas (F1.2.0; ADR-0101): each
/// [`yee_layout::Layout`] trace polygon drawn as a filled
/// [`egui_plot::Polygon`] in millimetres, with equal data aspect so the
/// geometry is not distorted. Shows the [`yee_filter::DimError`] string on
/// `Err`.
///
/// `Layout.traces` vertices are in metres; they are scaled to millimetres
/// (`×1e3`) for the plot. The `StudioState` itself stays egui-free — only this
/// `app` module (behind `desktop`/`web`) touches `egui_plot`.
fn show_layout(ui: &mut egui::Ui, state: &StudioState) {
    /// Metres → millimetres for the canvas coordinates.
    const MM: f64 = 1.0e3;

    match &state.layout {
        Ok(layout) => {
            let metal_fill = Color32::from_rgba_unmultiplied(200, 140, 40, 160);
            Plot::new("studio_layout_plot")
                .data_aspect(1.0)
                .x_axis_label("x (mm)")
                .y_axis_label("y (mm)")
                .legend(Legend::default())
                .show(ui, |plot_ui| {
                    // One distinct-named polygon per trace so each metal
                    // footprint gets its own egui_plot item Id (egui_plot
                    // derives the Id from the name — shared names hover as one).
                    for (i, poly) in layout.traces.iter().enumerate() {
                        let pts: Vec<[f64; 2]> =
                            poly.verts.iter().map(|p| [p.x * MM, p.y * MM]).collect();
                        plot_ui.polygon(
                            Polygon::new(format!("trace {i}"), PlotPoints::from(pts))
                                .fill_color(metal_fill),
                        );
                    }
                });
        }
        Err(msg) => {
            ui.colored_label(Color32::from_rgb(220, 60, 60), format!("no layout: {msg}"));
        }
    }
}

/// Draw the `|S21|` (dB) trace over the sweep, with each spec-mask region
/// shaded as a translucent box on its forbidden side.
fn show_response_plot(ui: &mut egui::Ui, state: &StudioState) {
    // Plot-vertical extent for the shaded boxes: take the data range padded so
    // the floor/ceiling boxes reach the visible edges.
    let mut y_min = -120.0_f64;
    let mut y_max = 5.0_f64;
    for &db in &state.s21_db {
        if db.is_finite() {
            y_min = y_min.min(db - 10.0);
            y_max = y_max.max(db + 5.0);
        }
    }

    let trace: Vec<[f64; 2]> = state
        .freqs_hz
        .iter()
        .zip(state.s21_db.iter())
        .map(|(&f, &db)| [f * 1.0e-9, db])
        .collect();

    let floor_fill = Color32::from_rgba_unmultiplied(220, 60, 60, 40);
    let ceil_fill = Color32::from_rgba_unmultiplied(60, 90, 220, 40);

    Plot::new("studio_s21_plot")
        .x_axis_label("Frequency (GHz)")
        .y_axis_label("|S21| (dB)")
        .legend(Legend::default())
        .show(ui, |plot_ui| {
            // Shade each forbidden region. Names carry the loop index so each
            // polygon gets a distinct egui Id (egui_plot derives the item Id
            // from the name — shared names make all co-named boxes hover as one
            // once there is more than one stopband point).
            for (i, region) in state.mask_regions.iter().enumerate() {
                let box_pts = forbidden_box(region, y_min, y_max);
                let (name, fill) = if region.floor {
                    (format!("passband floor {i}"), floor_fill)
                } else {
                    (format!("stopband ceiling {i}"), ceil_fill)
                };
                plot_ui.polygon(
                    Polygon::new(name, PlotPoints::from(box_pts))
                        .fill_color(fill)
                        .stroke(egui::Stroke::NONE),
                );
            }
            // |S21| trace on top.
            plot_ui.line(Line::new("|S21| (dB)", PlotPoints::from(trace)));
        });
}

/// Corner points (GHz, dB) of the translucent box covering a region's
/// forbidden side: below `limit_db` for a floor, above it for a ceiling.
fn forbidden_box(region: &MaskRegionView, y_min: f64, y_max: f64) -> Vec<[f64; 2]> {
    let x_lo = region.f_lo_hz * 1.0e-9;
    let x_hi = region.f_hi_hz * 1.0e-9;
    let (y_lo, y_hi) = if region.floor {
        // Floor: forbidden area is below the limit.
        (y_min, region.limit_db)
    } else {
        // Ceiling: forbidden area is above the limit.
        (region.limit_db, y_max)
    };
    vec![[x_lo, y_lo], [x_hi, y_lo], [x_hi, y_hi], [x_lo, y_hi]]
}
