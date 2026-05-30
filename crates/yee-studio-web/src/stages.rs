//! Stage definitions + per-stage renderers for the Shell A canvas.
//!
//! Two stages are **real** (driven by the live [`crate::engine`]):
//! [`synthesis_stage`] and [`layout_stage`]. The rest ([`spec_stage`],
//! [`technique_stage`], [`verify_stage`], [`export_stage`]) are styled-but-
//! static stubs that prove the shell, per the POC scope.

use dioxus::prelude::*;

use crate::engine::Designed;
use crate::svg::{board_svg, response_plot};

/// The six product stages (the left rail order).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stage {
    /// Spec: f0 / bandwidth / order / ripple / mask / Z0.
    Spec,
    /// Technique: topology gallery + medium + substrate library.
    Technique,
    /// Synthesis: g-values, Qe, coupling matrix, ideal response vs mask.
    Synthesis,
    /// Layout + Materials: board top-view, stackup, resonator table.
    Layout,
    /// Verify (EM): FDTD realized response (later).
    Verify,
    /// Export: Gerber / KiCad / Touchstone / STEP + parameter sheet.
    Export,
}

impl Stage {
    /// All stages in rail order.
    pub const ALL: [Stage; 6] = [
        Stage::Spec,
        Stage::Technique,
        Stage::Synthesis,
        Stage::Layout,
        Stage::Verify,
        Stage::Export,
    ];

    /// The short rail label.
    pub fn label(self) -> &'static str {
        match self {
            Stage::Spec => "Spec",
            Stage::Technique => "Technique",
            Stage::Synthesis => "Synthesis",
            Stage::Layout => "Layout",
            Stage::Verify => "Verify",
            Stage::Export => "Export",
        }
    }

    /// The rail glyph (matches the brainstorm mockups).
    pub fn icon(self) -> &'static str {
        match self {
            Stage::Spec => "◈",
            Stage::Technique => "▦",
            Stage::Synthesis => "∿",
            Stage::Layout => "▤",
            Stage::Verify => "◎",
            Stage::Export => "⤓",
        }
    }
}

// ===========================================================================
// REAL stage — Synthesis
// ===========================================================================

/// Render the Synthesis stage from live engine output: g-values, external Q,
/// the coupling matrix `M` as a styled grid, the swept ideal `|S21|`/`|S11|` vs
/// the shaded spec mask (inline SVG), and the real PASS/FAIL verdict.
pub fn synthesis_stage(designed: ReadOnlySignal<Designed>) -> Element {
    let d = designed.read();
    let n = d.order();

    // Coupling-matrix cells (real M).
    let m = &d.coupling.m;
    let plot = response_plot(&d.sweep, &d.mask_bands);

    // Mask report lines.
    let report = &d.report;
    let stopband_rows: Vec<(String, String, String, bool)> = report
        .stopband
        .iter()
        .map(|&(f, achieved, required, met)| {
            (
                format!("{:.3} GHz", f / 1e9),
                format!("{achieved:.1} dB"),
                format!("{required:.1} dB"),
                met,
            )
        })
        .collect();

    rsx! {
        div { class: "canvas-head",
            h1 { "Synthesis" }
            p { class: "sub", "Lowpass prototype, external Q, and the all-pole coupling matrix — graded against the spec mask. All values are live engine output." }
        }

        // ---- response plot ------------------------------------------------
        div { class: "card", style: "margin-bottom:16px",
            h2 { class: "card-title",
                "Ideal response vs spec mask"
                span { class: "k", "closed-form prototype · |S21|, |S11|" }
            }
            div { class: "plot", dangerous_inner_html: "{plot}" }
            div { class: "legend",
                span { span { class: "swatch", style: "background:#2dd4bf" } "|S21| (transmission)" }
                span { span { class: "swatch", style: "background:#6b7480" } "|S11| (reflection)" }
                span { span { class: "swatch", style: "background:#e35d6a" } "forbidden (mask)" }
            }
        }

        div { class: "row",
            // ---- coupling matrix --------------------------------------------
            div { class: "card", style: "flex:1.2",
                h2 { class: "card-title",
                    "Coupling matrix M"
                    span { class: "k", "normalized · {n}×{n} · synchronous" }
                }
                div {
                    class: "matrix",
                    style: "grid-template-columns: repeat({n}, auto)",
                    for (i, row) in m.iter().enumerate() {
                        for (j, v) in row.iter().enumerate() {
                            {
                                let cls = if i == j {
                                    "cell diag"
                                } else if j == i + 1 || i == j + 1 {
                                    "cell coupled"
                                } else {
                                    "cell"
                                };
                                rsx! {
                                    div { key: "{i}-{j}", class: "{cls}", "{v:+.4}" }
                                }
                            }
                        }
                    }
                }
            }

            // ---- prototype + Qe ---------------------------------------------
            div { class: "card", style: "flex:1",
                h2 { class: "card-title", "Prototype + external Q" }
                p { class: "lab", "external Q" }
                div { class: "stats",
                    div { class: "stat",
                        div { class: "v", "{d.coupling.qe_in:.3}" }
                        div { class: "l", "Qe (input)" }
                    }
                    div { class: "stat",
                        div { class: "v", "{d.coupling.qe_out:.3}" }
                        div { class: "l", "Qe (output)" }
                    }
                }
                p { class: "lab", style: "margin-top:16px", "g-values" }
                div { class: "gvals",
                    for (i, g) in d.g_values.iter().enumerate() {
                        div { key: "g{i}", class: "gval", b { "g{i} " } "{g:.4}" }
                    }
                }
            }
        }

        // ---- mask verdict -------------------------------------------------
        div { class: "card", style: "margin-top:16px",
            h2 { class: "card-title",
                "Spec-mask verdict"
                if report.pass {
                    span { class: "chip pass", style: "margin-left:auto", "PASS" }
                } else {
                    span { class: "chip fail", style: "margin-left:auto", "FAIL" }
                }
            }
            div { class: "stats", style: "margin-bottom:12px",
                div { class: "stat",
                    div { class: "v", "{report.worst_passband_ripple_db:.3} dB" }
                    div { class: "l", "passband ripple (spec {d.spec.mask.passband_ripple_db:.2})" }
                }
                div { class: "stat",
                    div { class: "v", "{report.worst_return_loss_db:.2} dB" }
                    div { class: "l", "in-band return loss (spec {d.spec.mask.return_loss_db:.2})" }
                }
            }
            table {
                thead {
                    tr {
                        th { "stopband point" }
                        th { "achieved" }
                        th { "required" }
                        th { "" }
                    }
                }
                tbody {
                    for (i, (f, achieved, required, met)) in stopband_rows.iter().enumerate() {
                        tr { key: "sb{i}",
                            td { class: "mono", "{f}" }
                            td { class: "mono", "{achieved}" }
                            td { class: "mono", "{required}" }
                            td {
                                if *met {
                                    span { class: "ok-mark", "✓ met" }
                                } else {
                                    span { style: "color:#e35d6a", "✗ under" }
                                }
                            }
                        }
                    }
                }
            }
        }

        div { class: "note honest",
            "Honest note: this is the ideal closed-form prototype response. The full-wave "
            "realized response (metal thickness, loss, dispersion) is the EM-verify stage "
            "(later) — the studio never hides the ideal-vs-realized gap."
        }
    }
}

// ===========================================================================
// REAL stage — Layout + Materials
// ===========================================================================

/// Render the Layout + Materials stage from live engine output: the dimensioned
/// board top-view (inline SVG), the material stackup cross-section, and the
/// per-resonator components table (`W`, `L`, gap, `Z0e`/`Z0o`, `εeff`, realized
/// `k`).
pub fn layout_stage(designed: ReadOnlySignal<Designed>) -> Element {
    let d = designed.read();
    let board = board_svg(&d.layout);
    let sub = &d.layout.substrate;
    let (bw, bh) = d.board_size_mm;

    rsx! {
        div { class: "canvas-head",
            h1 { "Layout + Materials" }
            p { class: "sub", "Dimensioned edge-coupled board, the material stackup that feeds the even/odd models, and the per-resonator geometry — all from the live dimensional synthesis." }
        }

        div { class: "row",
            // ---- board top view ---------------------------------------------
            div { class: "card", style: "flex:1.5",
                h2 { class: "card-title",
                    "Board · top view"
                    span { class: "k", "F.Cu · substrate · ports · {bw:.1} × {bh:.1} mm" }
                }
                div { class: "board-frame", dangerous_inner_html: "{board}" }
                div { class: "legend-row",
                    span { class: "sw-cu", "● copper" }
                    span { class: "sw-sub", "● substrate" }
                    span { style: "color:#2dd4bf", "◯ port" }
                }
            }

            // ---- material stackup -------------------------------------------
            div { class: "card", style: "flex:0 0 260px",
                h2 { class: "card-title", "Material stackup" }
                div { class: "stack-layer",
                    span { class: "lbl", "F.Cu" }
                    span { class: "swatch", style: "background:#e6b24d;height:8px", "{sub.metal_thickness_m*1e6:.0} µm" }
                }
                div { class: "stack-layer",
                    span { class: "lbl", "substrate" }
                    span { class: "swatch", style: "background:#3f9e72;height:42px", "FR-4 · εr {sub.eps_r:.1} · {sub.height_m*1e3:.2} mm" }
                }
                div { class: "stack-layer",
                    span { class: "lbl", "GND" }
                    span { class: "swatch", style: "background:#e6b24d;height:8px", "{sub.metal_thickness_m*1e6:.0} µm" }
                }
                div { style: "margin-top:14px",
                    div { class: "editrow", span { "Substrate" } span { class: "pill-sel", "FR-4" } }
                    div { class: "editrow", span { "εr" } span { class: "v", "{sub.eps_r:.2}" } }
                    div { class: "editrow", span { "height h" } span { class: "v", "{sub.height_m*1e3:.2} mm" } }
                    div { class: "editrow", span { "Cu thickness" } span { class: "v", "{sub.metal_thickness_m*1e6:.0} µm" } }
                    div { class: "editrow", span { "loss tan δ" } span { class: "v", "{sub.loss_tangent:.3}" } }
                    div { class: "editrow", span { "line εeff" } span { class: "v", "{d.line_eps_eff:.3}" } }
                }
                p { class: "lab", style: "margin-top:10px;text-transform:none;letter-spacing:0",
                    "Library: FR-4 · Rogers RO4350B · alumina · …"
                }
            }
        }

        // ---- resonator components table -----------------------------------
        div { class: "card", style: "margin-top:16px",
            h2 { class: "card-title",
                "Components · resonators"
                span { class: "k", "edge-coupled ½λ · W / L / gap → Z0e/Z0o · εeff · realized k" }
            }
            table {
                thead {
                    tr {
                        th { "id" }
                        th { "W (mm)" }
                        th { "length (mm)" }
                        th { "gap→next (mm)" }
                        th { "Z0e / Z0o (Ω)" }
                        th { "εeff (e / o)" }
                        th { "target k" }
                        th { "realized k" }
                    }
                }
                tbody {
                    for r in d.resonators.iter() {
                        tr { key: "R{r.id}",
                            td { class: "mono", "R{r.id}" }
                            td { class: "mono", "{r.width_mm:.3}" }
                            td { class: "mono", "{r.length_mm:.2}" }
                            td { class: "mono",
                                match r.gap_to_next_mm { Some(g) => format!("{g:.3}"), None => "—".into() }
                            }
                            td { class: "mono",
                                match (r.z0e_ohm, r.z0o_ohm) {
                                    (Some(e), Some(o)) => format!("{e:.1} / {o:.1}"),
                                    _ => "—".into(),
                                }
                            }
                            td { class: "mono",
                                match (r.eps_eff_e, r.eps_eff_o) {
                                    (Some(e), Some(o)) => format!("{e:.2} / {o:.2}"),
                                    _ => "—".into(),
                                }
                            }
                            td { class: "mono",
                                match r.target_k { Some(k) => format!("{k:.4}"), None => "—".into() }
                            }
                            td { class: "mono",
                                match r.realized_k {
                                    Some(k) => rsx! { "{k:.4} " span { class: "ok-mark", "✓" } },
                                    None => rsx! { "—" },
                                }
                            }
                        }
                    }
                }
            }
            p { class: "note honest",
                "Materials (εr, h) and strips (W, gap) drive the Kirschning-Jansen even/odd "
                "Z and εeff per row; the realized k is recovered from the solved gap and "
                "matches the synthesized target."
            }
        }
    }
}

// ===========================================================================
// Styled-but-static stubs — Spec / Technique / Verify / Export
// ===========================================================================

/// Spec stage stub: the design-intent fields, statically rendered from the
/// demo spec.
pub fn spec_stage(designed: ReadOnlySignal<Designed>) -> Element {
    let d = designed.read();
    let s = &d.spec;
    let stop = s
        .mask
        .stopband
        .iter()
        .map(|(f, r)| format!("{:.2} GHz ≥ {r:.0} dB", f / 1e9))
        .collect::<Vec<_>>()
        .join(", ");
    rsx! {
        div { class: "canvas-head",
            h1 { "Spec" }
            p { class: "sub", "The design intent the synthesis consumes. (POC: read-only; full editing is App.D.1.)" }
        }
        div { class: "row",
            div { class: "card", style: "flex:1",
                h2 { class: "card-title", "Requirements" }
                div { class: "fields",
                    div { class: "field", span { class: "name", "Response" } span { class: "val", "Bandpass" } }
                    div { class: "field", span { class: "name", "Approximation" } span { class: "val", "Chebyshev 0.5 dB" } }
                    div { class: "field", span { class: "name", "Centre f0" } span { class: "val", "{s.f0_hz/1e9:.3} GHz" } }
                    div { class: "field", span { class: "name", "Fractional bandwidth" } span { class: "val", "{s.fbw*100.0:.0}%" } }
                    div { class: "field", span { class: "name", "Order N" } span { class: "val", "{d.order()}" } }
                    div { class: "field", span { class: "name", "System Z0" } span { class: "val", "{s.z0_ohm:.0} Ω" } }
                }
            }
            div { class: "card", style: "flex:1",
                h2 { class: "card-title", "Spec mask" }
                div { class: "fields",
                    div { class: "field", span { class: "name", "Passband ripple" } span { class: "val", "≤ {s.mask.passband_ripple_db:.2} dB" } }
                    div { class: "field", span { class: "name", "Return loss" } span { class: "val", "≥ {s.mask.return_loss_db:.0} dB" } }
                    div { class: "field", span { class: "name", "Stopband" } span { class: "val", "{stop}" } }
                }
                div { class: "note", "Live realizability check + editable fields land in App.D.1." }
            }
        }
    }
}

/// Technique stage stub: the topology gallery (edge-coupled selected; the rest
/// greyed "Soon"). Static glyphs, matching the brainstorm mockup.
pub fn technique_stage() -> Element {
    // (name, desc, glyph svg-inner, available, selected)
    let cards: [(&str, &str, &str, bool, bool); 6] = [
        (
            "Edge-coupled",
            "½λ parallel resonators · F1.2.0",
            r##"<svg viewBox="0 0 120 54"><g stroke="#e6b24d" stroke-width="4" fill="none"><line x1="10" y1="16" x2="55" y2="16"/><line x1="35" y1="30" x2="80" y2="30"/><line x1="60" y1="16" x2="105" y2="16"/><line x1="35" y1="44" x2="80" y2="44"/></g></svg>"##,
            true,
            true,
        ),
        (
            "Hairpin",
            "U-folded ½λ · compact · F1.2.2",
            r##"<svg viewBox="0 0 120 54"><g stroke="#e6b24d" stroke-width="4" fill="none"><path d="M14,44 L14,14 L30,14 L30,44"/><path d="M44,44 L44,14 L60,14 L60,44"/><path d="M74,44 L74,14 L90,14 L90,44"/></g></svg>"##,
            true,
            false,
        ),
        (
            "Combline",
            "grounded ¼λ + via",
            r##"<svg viewBox="0 0 120 54"><g stroke="#6b7480" stroke-width="4" fill="none"><line x1="22" y1="10" x2="22" y2="40"/><line x1="46" y1="14" x2="46" y2="44"/><line x1="70" y1="10" x2="70" y2="40"/><line x1="94" y1="14" x2="94" y2="44"/></g></svg>"##,
            false,
            false,
        ),
        (
            "Interdigital",
            "interleaved grounded fingers",
            r##"<svg viewBox="0 0 120 54"><g stroke="#6b7480" stroke-width="4" fill="none"><path d="M16,44 L16,12"/><path d="M34,10 L34,42"/><path d="M52,44 L52,12"/><path d="M70,10 L70,42"/><path d="M88,44 L88,12"/></g></svg>"##,
            false,
            false,
        ),
        (
            "Lumped LC",
            "discrete L/C ladder + BOM",
            r##"<svg viewBox="0 0 120 54"><g stroke="#6b7480" stroke-width="3" fill="none"><path d="M8,30 h18 m0,-8 v16 m6,-16 v16 m6,-8 h14"/><path d="M52,30 q8,-16 16,0"/><path d="M76,30 h16 m0,-8 v16 m6,-16 v16"/></g></svg>"##,
            false,
            false,
        ),
        (
            "Stepped-impedance",
            "hi/lo-Z stub sections",
            r##"<svg viewBox="0 0 120 54"><g stroke="#6b7480" stroke-width="4" fill="none"><line x1="10" y1="30" x2="40" y2="30"/><rect x="40" y="20" width="16" height="20"/><line x1="56" y1="30" x2="72" y2="30"/><rect x="72" y="22" width="10" height="16"/><line x1="82" y1="30" x2="110" y2="30"/></g></svg>"##,
            false,
            false,
        ),
    ];
    rsx! {
        div { class: "canvas-head",
            h1 { "Technique" }
            p { class: "sub", "Each topology realizes the same coupling matrix in a different geometry. Available ones are validated end-to-end; the rest are roadmap placeholders. (POC: static gallery.)" }
        }
        div { class: "tgrid",
            for (name, desc, glyph, avail, sel) in cards {
                {
                    let cls = if sel { "tcard sel" } else if avail { "tcard" } else { "tcard soon" };
                    let badge_cls = if avail { "badge" } else { "badge soon" };
                    let badge = if avail { "Available" } else { "Soon" };
                    rsx! {
                        div { key: "{name}", class: "{cls}",
                            span { class: "{badge_cls}", "{badge}" }
                            div { class: "glyph", dangerous_inner_html: "{glyph}" }
                            div { class: "name", "{name}" }
                            div { class: "desc", "{desc}" }
                        }
                    }
                }
            }
        }
    }
}

/// Verify (EM) stage stub: the ideal-vs-realized story, marked "later".
pub fn verify_stage() -> Element {
    rsx! {
        div { class: "canvas-head",
            h1 { "Verify (EM)" }
            p { class: "sub", "Built-in full-wave FDTD verification — the differentiator the calculators leave out." }
        }
        div { class: "card",
            h2 { class: "card-title",
                "FDTD realized response"
                span { class: "chip fail", style: "margin-left:auto;background:#1b2027;border-color:#2a313b;color:#8b95a1", "Soon · FDTD in-loop" }
            }
            div { class: "grid-3",
                div { class: "stat", div { class: "v", "—" } div { class: "l", "insertion loss @ f0" } }
                div { class: "stat", div { class: "v", "—" } div { class: "l", "in-band return loss" } }
                div { class: "stat", div { class: "v", "—" } div { class: "l", "rejection @ 2.4 GHz" } }
            }
            div { class: "note honest",
                "This stage will overlay the as-built FDTD response on the ideal prototype + "
                "mask, closing the ideal-vs-realized gap, then tune / auto-optimize. It rides "
                "on the F1.1b / F1.2.1 FDTD-in-loop work (App.D.5)."
            }
        }
    }
}

/// Export stage stub: the manufacturable-output buttons + a parameter-sheet
/// teaser.
pub fn export_stage(designed: ReadOnlySignal<Designed>) -> Element {
    let d = designed.read();
    let (bw, bh) = d.board_size_mm;
    rsx! {
        div { class: "canvas-head",
            h1 { "Export" }
            p { class: "sub", "Manufacturable files + the final parameter sheet — the end state. (POC: stubbed; wiring is App.D.4.)" }
        }
        div { class: "card",
            h2 { class: "card-title", "Design summary" }
            div { class: "fields",
                div { class: "field", span { class: "name", "Topology" } span { class: "val", "edge-coupled ½λ · N={d.order()}" } }
                div { class: "field", span { class: "name", "Response" } span { class: "val", "Chebyshev 0.5 dB" } }
                div { class: "field", span { class: "name", "f0 / FBW" } span { class: "val", "{d.spec.f0_hz/1e9:.3} GHz / {d.spec.fbw*100.0:.0}%" } }
                div { class: "field", span { class: "name", "Substrate" } span { class: "val", "FR-4 · εr {d.layout.substrate.eps_r:.1} · h {d.layout.substrate.height_m*1e3:.2} mm" } }
                div { class: "field", span { class: "name", "Board" } span { class: "val", "{bw:.1} × {bh:.1} mm" } }
            }
            div { class: "export-row",
                span { class: "btn", "⤓ Gerber" }
                span { class: "btn", "⤓ KiCad .kicad_pcb" }
                span { class: "btn", "⤓ Touchstone .s2p" }
                span { class: "btn", "⤓ STEP" }
            }
            div { class: "note", "Each exporter already exists in `yee-export` / `yee-io`; App.D.4 wires the download buttons to them." }
        }
    }
}
