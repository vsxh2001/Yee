//! Stage definitions + per-stage renderers for the Shell A canvas.
//!
//! Two stages are **real** (driven by the live [`crate::engine`]):
//! [`synthesis_stage`] and [`layout_stage`]. The rest ([`spec_stage`],
//! [`technique_stage`], [`verify_stage`], [`export_stage`]) are styled-but-
//! static stubs that prove the shell, per the POC scope.

use dioxus::prelude::*;

use crate::engine::{BomView, Designed, LumpedDesigned, YieldView};
use crate::svg::{board_svg, lumped_board_svg, response_plot};

/// The realization technique the downstream stages render for.
///
/// Selecting [`Topology::LumpedLc`] on the Technique stage routes Synthesis /
/// Components / Tolerance / Layout to the lumped-LC renderers and swaps the rail
/// for the lumped flow; [`Topology::EdgeCoupled`] keeps the distributed flow.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Topology {
    /// Distributed edge-coupled microstrip (the POC's original flow).
    EdgeCoupled,
    /// Lumped-element LC ladder (ADR-0120: synth → BOM → tolerance → board).
    LumpedLc,
}

/// The product stages (the left rail order). Two stages —
/// [`Stage::Components`] and [`Stage::Tolerance`] — exist only in the lumped
/// flow; the rail filters by the active [`Topology`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stage {
    /// Spec: f0 / bandwidth / order / ripple / mask / Z0.
    Spec,
    /// Technique: topology gallery + medium + substrate library.
    Technique,
    /// Synthesis: coupling matrix (distributed) or LC ladder (lumped) vs mask.
    Synthesis,
    /// Components + BOM: E-series part selection (lumped flow only).
    Components,
    /// Tolerance / yield: Monte-Carlo yield analysis (lumped flow only).
    Tolerance,
    /// Layout + Materials: board top-view, stackup / footprints, tables.
    Layout,
    /// Verify (EM): FDTD realized response (later).
    Verify,
    /// Export: Gerber / KiCad / Touchstone / STEP + parameter sheet.
    Export,
}

impl Stage {
    /// The distributed-flow rail (the POC's original six).
    const DISTRIBUTED: [Stage; 6] = [
        Stage::Spec,
        Stage::Technique,
        Stage::Synthesis,
        Stage::Layout,
        Stage::Verify,
        Stage::Export,
    ];

    /// The lumped-flow rail (ADR-0120): adds Components + Tolerance.
    const LUMPED: [Stage; 7] = [
        Stage::Spec,
        Stage::Technique,
        Stage::Synthesis,
        Stage::Components,
        Stage::Tolerance,
        Stage::Layout,
        Stage::Export,
    ];

    /// The rail order for the active topology.
    pub fn rail(topology: Topology) -> &'static [Stage] {
        match topology {
            Topology::EdgeCoupled => &Stage::DISTRIBUTED,
            Topology::LumpedLc => &Stage::LUMPED,
        }
    }

    /// The short rail label.
    pub fn label(self) -> &'static str {
        match self {
            Stage::Spec => "Spec",
            Stage::Technique => "Technique",
            Stage::Synthesis => "Synthesis",
            Stage::Components => "Components",
            Stage::Tolerance => "Tolerance",
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
            Stage::Components => "▥",
            Stage::Tolerance => "⤬",
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

/// One Technique gallery card's static descriptor (the topology it selects, or
/// `None` for the greyed roadmap placeholders).
struct TechCard {
    /// Display name.
    name: &'static str,
    /// One-line description.
    desc: &'static str,
    /// Inline `<svg>` glyph body.
    glyph: &'static str,
    /// The topology this card selects when live (`None` → greyed "Soon").
    selects: Option<Topology>,
}

/// Technique stage: the topology gallery. **Edge-coupled** and **Lumped LC** are
/// live and selectable (clicking a live card routes the downstream stages via
/// `topology`); the rest stay greyed "Soon". The currently selected topology is
/// highlighted.
pub fn technique_stage(mut topology: Signal<Topology>, mut active: Signal<Stage>) -> Element {
    let cards: [TechCard; 6] = [
        TechCard {
            name: "Edge-coupled",
            desc: "½λ parallel resonators · F1.2.0",
            glyph: r##"<svg viewBox="0 0 120 54"><g stroke="#e6b24d" stroke-width="4" fill="none"><line x1="10" y1="16" x2="55" y2="16"/><line x1="35" y1="30" x2="80" y2="30"/><line x1="60" y1="16" x2="105" y2="16"/><line x1="35" y1="44" x2="80" y2="44"/></g></svg>"##,
            selects: Some(Topology::EdgeCoupled),
        },
        TechCard {
            name: "Lumped LC",
            desc: "discrete L/C ladder · BOM · tolerance · F2.0–F2.4",
            glyph: r##"<svg viewBox="0 0 120 54"><g stroke="#2dd4bf" stroke-width="3" fill="none"><path d="M8,30 h18 m0,-8 v16 m6,-16 v16 m6,-8 h14"/><path d="M52,30 q8,-16 16,0"/><path d="M76,30 h16 m0,-8 v16 m6,-16 v16"/></g></svg>"##,
            selects: Some(Topology::LumpedLc),
        },
        TechCard {
            name: "Hairpin",
            desc: "U-folded ½λ · compact · F1.2.2",
            glyph: r##"<svg viewBox="0 0 120 54"><g stroke="#6b7480" stroke-width="4" fill="none"><path d="M14,44 L14,14 L30,14 L30,44"/><path d="M44,44 L44,14 L60,14 L60,44"/><path d="M74,44 L74,14 L90,14 L90,44"/></g></svg>"##,
            selects: None,
        },
        TechCard {
            name: "Combline",
            desc: "grounded ¼λ + via",
            glyph: r##"<svg viewBox="0 0 120 54"><g stroke="#6b7480" stroke-width="4" fill="none"><line x1="22" y1="10" x2="22" y2="40"/><line x1="46" y1="14" x2="46" y2="44"/><line x1="70" y1="10" x2="70" y2="40"/><line x1="94" y1="14" x2="94" y2="44"/></g></svg>"##,
            selects: None,
        },
        TechCard {
            name: "Interdigital",
            desc: "interleaved grounded fingers",
            glyph: r##"<svg viewBox="0 0 120 54"><g stroke="#6b7480" stroke-width="4" fill="none"><path d="M16,44 L16,12"/><path d="M34,10 L34,42"/><path d="M52,44 L52,12"/><path d="M70,10 L70,42"/><path d="M88,44 L88,12"/></g></svg>"##,
            selects: None,
        },
        TechCard {
            name: "Stepped-impedance",
            desc: "hi/lo-Z stub sections",
            glyph: r##"<svg viewBox="0 0 120 54"><g stroke="#6b7480" stroke-width="4" fill="none"><line x1="10" y1="30" x2="40" y2="30"/><rect x="40" y="20" width="16" height="20"/><line x1="56" y1="30" x2="72" y2="30"/><rect x="72" y="22" width="10" height="16"/><line x1="82" y1="30" x2="110" y2="30"/></g></svg>"##,
            selects: None,
        },
    ];
    let cur = topology();
    rsx! {
        div { class: "canvas-head",
            h1 { "Technique" }
            p { class: "sub", "Each topology realizes the same prototype in a different way. Edge-coupled (distributed) and Lumped LC (discrete parts + BOM + tolerance) are live — click to route the flow; the rest are roadmap placeholders." }
        }
        div { class: "tgrid",
            for card in cards {
                {
                    let avail = card.selects.is_some();
                    let sel = card.selects == Some(cur);
                    let cls = if sel { "tcard sel" } else if avail { "tcard" } else { "tcard soon" };
                    let badge_cls = if avail { "badge" } else { "badge soon" };
                    let badge = if sel { "Selected" } else if avail { "Available" } else { "Soon" };
                    let target = card.selects;
                    rsx! {
                        div {
                            key: "{card.name}",
                            class: "{cls}",
                            onclick: move |_| {
                                if let Some(t) = target
                                    && topology() != t
                                {
                                    topology.set(t);
                                    // The new flow may not contain the current
                                    // stage — land on Synthesis.
                                    active.set(Stage::Synthesis);
                                }
                            },
                            span { class: "{badge_cls}", "{badge}" }
                            div { class: "glyph", dangerous_inner_html: "{card.glyph}" }
                            div { class: "name", "{card.name}" }
                            div { class: "desc", "{card.desc}" }
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

// ===========================================================================
// REAL lumped-LC stages (App.D.1L) — Synthesis / Components / Tolerance / Layout
// ===========================================================================

/// Lumped Synthesis stage: the LC ladder resonator table (index, series/shunt,
/// L [nH], C [pF]) + the **ideal** `ladder_s21` |S21| vs the spec mask (inline
/// SVG, reusing the response plot) + a PASS/FAIL chip from the realized verdict.
pub fn lumped_synthesis_stage(designed: ReadOnlySignal<LumpedDesigned>) -> Element {
    let d = designed.read();
    let plot = response_plot(&d.sweep, &d.mask_bands);
    let v = &d.verdict;
    rsx! {
        div { class: "canvas-head",
            h1 { "Synthesis · Lumped LC" }
            p { class: "sub", "The lowpass prototype mapped to a band-pass LC ladder (shunt-first, alternating). Every resonator is tuned to f0 (L·C·ω0² = 1). All values are live engine output." }
        }

        // ---- ideal ladder response vs mask --------------------------------
        div { class: "card", style: "margin-bottom:16px",
            h2 { class: "card-title",
                "Ideal ladder response vs spec mask"
                span { class: "k", "ABCD cascade · |S21|, |S11|" }
            }
            div { class: "plot", dangerous_inner_html: "{plot}" }
            div { class: "legend",
                span { span { class: "swatch", style: "background:#2dd4bf" } "|S21| (transmission)" }
                span { span { class: "swatch", style: "background:#6b7480" } "|S11| (reflection)" }
                span { span { class: "swatch", style: "background:#e35d6a" } "forbidden (mask)" }
            }
        }

        div { class: "row",
            // ---- LC ladder table -------------------------------------------
            div { class: "card", style: "flex:1.4",
                h2 { class: "card-title",
                    "LC ladder"
                    span { class: "k", "N={d.order()} · shunt-first · tuned to f0" }
                }
                table {
                    thead {
                        tr {
                            th { "index" }
                            th { "branch" }
                            th { "L (nH)" }
                            th { "C (pF)" }
                        }
                    }
                    tbody {
                        for r in d.resonators.iter() {
                            tr { key: "lc{r.index}",
                                td { class: "mono", "{r.index}" }
                                td {
                                    if r.is_series {
                                        span { class: "pill-sel", "series" }
                                    } else {
                                        span { class: "pill-sel", style: "background:#11302a;border-color:#1f5138;color:#2dd4bf", "shunt" }
                                    }
                                }
                                td { class: "mono", "{r.l_nh:.3}" }
                                td { class: "mono", "{r.c_pf:.3}" }
                            }
                        }
                    }
                }
            }

            // ---- realized verdict ------------------------------------------
            div { class: "card", style: "flex:1",
                h2 { class: "card-title",
                    "Realized verdict"
                    if v.pass {
                        span { class: "chip pass", style: "margin-left:auto", "PASS" }
                    } else {
                        span { class: "chip fail", style: "margin-left:auto", "FAIL" }
                    }
                }
                div { class: "stats",
                    div { class: "stat",
                        div { class: "v", "{v.worst_passband_ripple_db:.3} dB" }
                        div { class: "l", "passband ripple" }
                    }
                    div { class: "stat",
                        div { class: "v", "{v.worst_return_loss_db:.2} dB" }
                        div { class: "l", "in-band return loss" }
                    }
                }
                div { class: "stats", style: "margin-top:12px",
                    div { class: "stat",
                        div { class: "v", "{v.worst_stopband_rej_db:.1} dB" }
                        div { class: "l", "stopband rejection" }
                    }
                }
            }
        }

        div { class: "note honest",
            "Honest note: this is the IDEAL LC prototype response (lossless components, "
            "exact values). Realized E-series parts, tolerance spread, and parasitics "
            "follow in Components / Tolerance; full-wave EM-verify is a later stage."
        }
    }
}

/// Lumped Components + BOM stage: an E24/E96 toggle and the grouped BOM table
/// (ref kind, ideal value, chosen E-series value, deviation %, tolerance %,
/// qty), deviation colour-coded against the series tolerance.
pub fn lumped_components_stage(
    designed: ReadOnlySignal<LumpedDesigned>,
    mut series_e96: Signal<bool>,
) -> Element {
    let d = designed.read();
    let use_e96 = series_e96();
    let bom: &BomView = if use_e96 { &d.bom_e96 } else { &d.bom_e24 };
    rsx! {
        div { class: "canvas-head",
            h1 { "Components + BOM" }
            p { class: "sub", "Each ideal L/C snapped to the nearest standard IEC 60063 preferred value (log-nearest), then grouped into a purchasable bill of materials. Switch series to trade tolerance for part count." }
        }

        div { class: "card", style: "margin-bottom:16px",
            h2 { class: "card-title",
                "E-series"
                span { class: "k", "log-nearest preferred value · ±tolerance" }
            }
            div { class: "seg",
                button {
                    class: if !use_e96 { "seg-btn on" } else { "seg-btn" },
                    onclick: move |_| series_e96.set(false),
                    "E24  (±5%)"
                }
                button {
                    class: if use_e96 { "seg-btn on" } else { "seg-btn" },
                    onclick: move |_| series_e96.set(true),
                    "E96  (±1%)"
                }
            }
            div { class: "stats", style: "margin-top:14px",
                div { class: "stat",
                    div { class: "v", "±{bom.tolerance_pct:.0}%" }
                    div { class: "l", "{bom.series_name} tolerance" }
                }
                div { class: "stat",
                    div { class: "v", "{bom.total_parts}" }
                    div { class: "l", "total parts" }
                }
                div { class: "stat",
                    div { class: "v", "{bom.worst_deviation_pct:.2}%" }
                    div { class: "l", "worst deviation" }
                }
            }
        }

        div { class: "card",
            h2 { class: "card-title",
                "Bill of materials"
                span { class: "k", "{bom.series_name} · grouped by (kind, value)" }
            }
            table {
                thead {
                    tr {
                        th { "ref" }
                        th { "kind" }
                        th { "ideal" }
                        th { "chosen ({bom.series_name})" }
                        th { "deviation" }
                        th { "tol" }
                        th { "qty" }
                    }
                }
                tbody {
                    for (i, r) in bom.rows.iter().enumerate() {
                        {
                            // Colour-code deviation: green if well inside tolerance,
                            // amber if approaching it, red if at/over.
                            let mag = r.deviation_pct.abs();
                            let dev_col = if mag <= 0.5 * r.tolerance_pct {
                                "#2dd4bf"
                            } else if mag <= r.tolerance_pct {
                                "#e6b24d"
                            } else {
                                "#e35d6a"
                            };
                            rsx! {
                                tr { key: "bom{i}",
                                    td { class: "mono", "{r.ref_kind}" }
                                    td {
                                        if r.is_inductor {
                                            span { class: "pill-sel", style: "background:#11302a;border-color:#1f5138;color:#2dd4bf", "inductor" }
                                        } else {
                                            span { class: "pill-sel", "capacitor" }
                                        }
                                    }
                                    td { class: "mono", "{r.ideal_disp}" }
                                    td { class: "mono", "{r.chosen_disp}" }
                                    td { class: "mono", style: "color:{dev_col}", "{r.deviation_pct:+.2}%" }
                                    td { class: "mono", "±{r.tolerance_pct:.0}%" }
                                    td { class: "mono", "{r.qty}" }
                                }
                            }
                        }
                    }
                }
            }
            p { class: "note honest",
                "Deviation is the chosen-vs-ideal error of the standard value itself; the "
                "±tolerance column is the part's manufacturing spread. Both feed the yield "
                "analysis on the Tolerance stage."
            }
        }
    }
}

/// Lumped Tolerance / yield stage: side-by-side E24 vs E96 Monte-Carlo yield
/// cards (yield %, worst-case return loss + stopband rejection) plus the honest
/// narrowband-yield note when yield is low.
pub fn lumped_tolerance_stage(designed: ReadOnlySignal<LumpedDesigned>) -> Element {
    let d = designed.read();
    let trials = d.yield_trials;
    let lowest = d.yield_e24.yield_pct.min(d.yield_e96.yield_pct);
    rsx! {
        div { class: "canvas-head",
            h1 { "Tolerance / yield" }
            p { class: "sub", "A seeded Monte-Carlo over {trials} trials: each part perturbed uniformly within its tolerance, the ladder rebuilt and re-graded against the spec mask. Yield is the fraction that still passes." }
        }
        div { class: "row",
            {
                let v24 = d.yield_e24;
                let v96 = d.yield_e96;
                rsx! {
                    yield_card { v: v24 }
                    yield_card { v: v96 }
                }
            }
        }
        if lowest < 90.0 {
            div { class: "note honest",
                "Honest note: a narrow-band lumped band-pass is yield-sensitive — small part "
                "errors shift the band edges and collapse return loss. To raise yield, use a "
                "tighter series (E96 ±1%) or relax the spec mask (wider band / lower return-loss "
                "target). This is the trade the calculators hide."
            }
        } else {
            div { class: "note",
                "Yield is healthy at this spec/tolerance — the realized parts meet the mask "
                "across the great majority of trials."
            }
        }
    }
}

/// One Monte-Carlo yield card (E24 or E96).
#[component]
fn yield_card(v: YieldView) -> Element {
    let yield_col = if v.yield_pct >= 90.0 {
        "#2dd4bf"
    } else if v.yield_pct >= 60.0 {
        "#e6b24d"
    } else {
        "#e35d6a"
    };
    rsx! {
        div { class: "card", style: "flex:1",
            h2 { class: "card-title",
                "{v.series_name} yield"
                span { class: "k", "±{v.tolerance_pct:.0}% parts" }
            }
            div { class: "stat", style: "margin-bottom:14px",
                div { class: "v", style: "font-size:38px;color:{yield_col}", "{v.yield_pct:.1}%" }
                div { class: "l", "pass the spec mask" }
            }
            div { class: "stats",
                div { class: "stat",
                    div { class: "v", "{v.worst_rl_db:.2} dB" }
                    div { class: "l", "worst-case return loss" }
                }
                div { class: "stat",
                    div { class: "v", "{v.worst_rej_db:.1} dB" }
                    div { class: "l", "worst-case rejection" }
                }
            }
        }
    }
}

/// Lumped Layout stage: the dimensioned SMD board (footprints to scale, pads,
/// signal trace, ground rail — inline SVG) + the placement/footprint table.
pub fn lumped_layout_stage(designed: ReadOnlySignal<LumpedDesigned>) -> Element {
    let d = designed.read();
    let board = lumped_board_svg(&d.board);
    let sub = &d.board.layout.substrate;
    let (bw, bh) = d.board_size_mm;
    rsx! {
        div { class: "canvas-head",
            h1 { "Layout · Lumped board" }
            p { class: "sub", "The LC ladder placed as SMD footprints on a Z0 microstrip with a ground rail — series parts in-line, shunt parts on stubs to ground. Footprints are drawn to scale from the live board generator." }
        }
        div { class: "row",
            // ---- board top view --------------------------------------------
            div { class: "card", style: "flex:1.6",
                h2 { class: "card-title",
                    "Board · top view"
                    span { class: "k", "0603 SMD · signal · ground · {bw:.1} × {bh:.1} mm" }
                }
                div { class: "board-frame", dangerous_inner_html: "{board}" }
                div { class: "legend-row",
                    span { class: "sw-cu", "● copper (pads / line / rail)" }
                    span { class: "sw-sub", "● substrate" }
                    span { style: "color:#2dd4bf", "◯ port" }
                }
            }
            // ---- stackup ----------------------------------------------------
            div { class: "card", style: "flex:0 0 240px",
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
                    div { class: "editrow", span { "footprint" } span { class: "v", "0603" } }
                }
            }
        }

        // ---- placement table ----------------------------------------------
        div { class: "card", style: "margin-top:16px",
            h2 { class: "card-title",
                "Placement"
                span { class: "k", "ref-des → footprint → branch → centre (board frame)" }
            }
            table {
                thead {
                    tr {
                        th { "ref-des" }
                        th { "footprint" }
                        th { "branch" }
                        th { "x (mm)" }
                        th { "y (mm)" }
                    }
                }
                tbody {
                    for (i, p) in d.placements.iter().enumerate() {
                        tr { key: "pl{i}",
                            td { class: "mono", "{p.ref_des}" }
                            td { class: "mono", "{p.footprint}" }
                            td {
                                if p.kind == "series" {
                                    span { class: "pill-sel", "series" }
                                } else {
                                    span { class: "pill-sel", style: "background:#11302a;border-color:#1f5138;color:#2dd4bf", "shunt" }
                                }
                            }
                            td { class: "mono", "{p.cx_mm:.2}" }
                            td { class: "mono", "{p.cy_mm:.2}" }
                        }
                    }
                }
            }
            p { class: "note honest",
                "The same Layout feeds the shipped Gerber / KiCad exporters (Export stage). "
                "Component values are keyed by ref-des from the BOM; routing / matched-meander "
                "and a parasitic-aware land library are documented follow-ons (F2.2b)."
            }
        }
    }
}
