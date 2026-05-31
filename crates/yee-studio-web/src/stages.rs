//! Stage definitions + per-stage renderers for the Shell A canvas.
//!
//! Almost every stage is **real** (driven by the live [`crate::engine`]):
//! [`spec_stage`] (a live editable form), [`technique_stage`] (the topology
//! gallery — Edge-coupled + Lumped LC live, the rest honest "Soon"),
//! [`synthesis_stage`] / [`layout_stage`] (distributed), the lumped quartet
//! ([`lumped_synthesis_stage`], [`lumped_components_stage`],
//! [`lumped_tolerance_stage`], [`lumped_layout_stage`]), [`export_stage`]
//! (a real parameter sheet + Gerber/KiCad/BOM downloads), and [`verify_stage`]
//! (the active flow's real circuit-level mask metrics — App.2.4 / ADR-0141 —
//! honest that full-wave EM of the board is a separate native step).

use dioxus::prelude::*;

use yee_filter::{
    Approximation, FilterSpec, RealizationTechnique, Response, SpecMask, TechniqueRecommendation,
    recommend_technique,
};

use crate::engine::{
    BomView, Designed, LumpedDesigned, SteppedLowpassDesigned, TechniqueComparison, VerifyLevel,
    YieldView, compare_techniques, verify_view,
};
use crate::svg::{board_svg, lumped_board_svg, response_plot};

/// The realization technique the downstream stages render for.
///
/// Selecting [`Topology::LumpedLc`] on the Technique stage routes Synthesis /
/// Components / Tolerance / Layout to the lumped-LC renderers and swaps the rail
/// for the lumped flow; [`Topology::EdgeCoupled`] and [`Topology::Hairpin`] both
/// keep the distributed flow (same six stages, same coupled-resonator synthesis)
/// and differ only in the geometry derived for the Layout / Export stages.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Topology {
    /// Distributed edge-coupled microstrip (the POC's original flow).
    EdgeCoupled,
    /// Distributed hairpin (U-folded ½λ) microstrip (ADR-0138, F1.2.2). Shares
    /// the edge-coupled coupled-resonator synthesis; only the realized geometry
    /// differs (U-folded λ/4 arms vs straight λ/2 lines).
    Hairpin,
    /// Lumped-element LC ladder (ADR-0120: synth → BOM → tolerance → board).
    LumpedLc,
    /// Distributed **stepped-impedance low-pass** (ADR-0139, F1.2.3): alternating
    /// high-Z / low-Z microstrip sections realizing a low-pass prototype. The
    /// first **low-pass** flow — its Spec form labels f0 as the cutoff and hides
    /// the fractional bandwidth, and Synthesis / Layout route to the stepped
    /// renderers (the distributed five-stage rail, no Components / Tolerance).
    SteppedImpedance,
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
    /// Verify: the active flow's real circuit-level mask metrics (App.2.4).
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

    /// The lumped-flow rail (ADR-0120): adds Components + Tolerance. Verify
    /// (ADR-0141) shows the realized LC ladder graded vs the mask — the lumped
    /// flow's strongest verification story — so it belongs in this rail too.
    const LUMPED: [Stage; 8] = [
        Stage::Spec,
        Stage::Technique,
        Stage::Synthesis,
        Stage::Components,
        Stage::Tolerance,
        Stage::Layout,
        Stage::Verify,
        Stage::Export,
    ];

    /// The rail order for the active topology.
    pub fn rail(topology: Topology) -> &'static [Stage] {
        match topology {
            // Stepped-impedance is distributed: the same five-stage rail (Spec /
            // Technique / Synthesis / Layout / Export — no Components / Tolerance).
            Topology::EdgeCoupled | Topology::Hairpin | Topology::SteppedImpedance => {
                &Stage::DISTRIBUTED
            }
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
    // The layout is `None` when the live spec over-couples beyond what FR-4 can
    // realize; show an honest note instead of an empty board.
    let Some(layout) = d.layout.as_ref() else {
        let msg = d
            .dim_error
            .clone()
            .unwrap_or_else(|| "geometry not realizable on FR-4".into());
        let topo = d.topology_name();
        return rsx! {
            div { class: "canvas-head",
                h1 { "Layout + Materials" }
                p { class: "sub", "Dimensioned {topo} board — derived live from the synthesized coupling matrix." }
            }
            div { class: "card",
                h2 { class: "card-title",
                    "Geometry not realizable"
                    span { class: "chip fail", style: "margin-left:auto", "FR-4" }
                }
                div { class: "note honest",
                    "The current spec over-couples the resonators beyond what an "
                    "FR-4 microstrip can physically realize (gaps would go non-positive). The "
                    "synthesized prototype + ideal response on the Synthesis stage stay valid — "
                    "raise the order, narrow the fractional bandwidth, or switch to the Lumped LC "
                    "technique. Engine note: "
                    span { class: "mono", "{msg}" }
                }
            }
        };
    };
    let board = board_svg(layout);
    let sub = &layout.substrate;
    let (bw, bh) = d.board_size_mm;
    let topo = d.topology_name();
    let length_label = d.length_label();

    rsx! {
        div { class: "canvas-head",
            h1 { "Layout + Materials" }
            p { class: "sub", "Dimensioned {topo} board, the material stackup that feeds the even/odd models, and the per-resonator geometry — all from the live dimensional synthesis." }
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
                span { class: "k", "{topo} · W / L / gap → Z0e/Z0o · εeff · realized k" }
            }
            table {
                thead {
                    tr {
                        th { "id" }
                        th { "W (mm)" }
                        th { "{length_label}" }
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

/// Spec stage: a **live editable** design-intent form. Every control writes
/// back into the shared [`yee_filter::FilterSpec`] signal; the whole studio (synthesis,
/// dimensioning, BOM, yield, board) re-derives on every edit. A live PASS/FAIL
/// chip + realizability note close the loop without leaving the stage.
///
/// `designed`/`lumped`/`stepped` are read-only views of the current
/// re-derivation, used for the live verdict + the realizability hints.
///
/// `topology` makes the form **low-pass-aware**: when the stepped-impedance
/// (low-pass) flow is active, the centre frequency is labelled **Cutoff**, the
/// fractional-bandwidth field is hidden (low-pass has no FBW), and the live
/// verdict is graded against the low-pass spec mask. Every other topology keeps
/// the band-pass form.
pub fn spec_stage(
    mut spec: Signal<yee_filter::FilterSpec>,
    topology: ReadOnlySignal<Topology>,
    designed: ReadOnlySignal<Designed>,
    lumped: ReadOnlySignal<Option<LumpedDesigned>>,
    stepped: ReadOnlySignal<SteppedLowpassDesigned>,
) -> Element {
    let s = spec.read().clone();
    let d = designed.read();
    let lowpass = topology() == Topology::SteppedImpedance;
    let chebyshev = matches!(s.approximation, yee_filter::Approximation::Chebyshev { .. });
    let ripple = match s.approximation {
        yee_filter::Approximation::Chebyshev { ripple_db } => ripple_db,
        yee_filter::Approximation::Butterworth => 0.5,
    };
    // First stopband point (the form edits a single point; the engine supports
    // a vector and the verdict grades every point).
    let (stop_f_ghz, stop_rej_db) = s
        .mask
        .stopband
        .first()
        .map(|&(f, r)| (f / 1e9, r))
        .unwrap_or((2.4, 40.0));

    // Live realizability: the lumped flow may be unrealizable; the distributed
    // flow surfaces a dim_error. Show both honestly.
    let lumped_ok = lumped.read().is_some();
    // The active flow's geometry-error + verdict: low-pass reads from the stepped
    // design, every band-pass flow from the distributed `Designed`.
    let (active_dim_err, active_topo, verdict_pass, achieved_rl) = if lowpass {
        let st = stepped.read();
        (
            st.dim_error.clone(),
            "Stepped-impedance",
            st.pass,
            st.worst_return_loss_db,
        )
    } else {
        (
            d.dim_error.clone(),
            d.topology_name(),
            d.report.pass,
            d.report.worst_return_loss_db,
        )
    };
    let freq_label = if lowpass {
        "Cutoff f_c (GHz)"
    } else {
        "Centre f0 (GHz)"
    };

    rsx! {
        div { class: "canvas-head",
            h1 { "Spec" }
            p { class: "sub", "The design intent the synthesis consumes — edit any field and the whole studio (synthesis, components, tolerance, board) re-derives live." }
        }
        div { class: "row",
            // ---- requirements ------------------------------------------------
            div { class: "card", style: "flex:1",
                h2 { class: "card-title", "Requirements" }
                div { class: "fields",
                    div { class: "field",
                        span { class: "name", "Response" }
                        span { class: "val", if lowpass { "Lowpass" } else { "Bandpass" } }
                    }
                    // approximation toggle
                    div { class: "field",
                        span { class: "name", "Approximation" }
                        div { class: "seg",
                            button {
                                class: if chebyshev { "seg-btn on" } else { "seg-btn" },
                                onclick: move |_| {
                                    spec.write().approximation =
                                        yee_filter::Approximation::Chebyshev { ripple_db: 0.5 };
                                },
                                "Chebyshev"
                            }
                            button {
                                class: if !chebyshev { "seg-btn on" } else { "seg-btn" },
                                onclick: move |_| {
                                    spec.write().approximation =
                                        yee_filter::Approximation::Butterworth;
                                },
                                "Butterworth"
                            }
                        }
                    }
                    {num_field(
                        freq_label, s.f0_hz / 1e9, 0.001, 0.1, 100.0,
                        move |v| spec.write().f0_hz = v * 1e9,
                    )}
                    // Fractional bandwidth is a band-pass-only concept; the
                    // stepped-impedance low-pass flow hides it.
                    if !lowpass {
                        {num_field(
                            "Fractional bandwidth (%)", s.fbw * 100.0, 0.1, 0.5, 80.0,
                            move |v| spec.write().fbw = (v / 100.0).max(1e-4),
                        )}
                    }
                    {int_field(
                        "Order N", s.order.unwrap_or(5), 1, 11,
                        move |n| spec.write().order = Some(n),
                    )}
                    {num_field(
                        "System Z0 (Ω)", s.z0_ohm, 1.0, 10.0, 200.0,
                        move |v| spec.write().z0_ohm = v,
                    )}
                }
            }
            // ---- spec mask ---------------------------------------------------
            div { class: "card", style: "flex:1",
                h2 { class: "card-title",
                    "Spec mask"
                    if verdict_pass {
                        span { class: "chip pass", style: "margin-left:auto", "PASS" }
                    } else {
                        span { class: "chip fail", style: "margin-left:auto", "FAIL" }
                    }
                }
                div { class: "fields",
                    if chebyshev {
                        {num_field(
                            "Passband ripple (dB)", ripple, 0.01, 0.01, 3.0,
                            move |v| {
                                let v = v.max(1e-3);
                                spec.write().approximation =
                                    yee_filter::Approximation::Chebyshev { ripple_db: v };
                                spec.write().mask.passband_ripple_db = v;
                            },
                        )}
                    } else {
                        div { class: "field",
                            span { class: "name", "Passband ripple" }
                            span { class: "val", "— (maximally flat)" }
                        }
                    }
                    {num_field(
                        "Return loss (dB)", s.mask.return_loss_db, 0.1, 1.0, 30.0,
                        move |v| spec.write().mask.return_loss_db = v,
                    )}
                    {num_field(
                        "Stopband f (GHz)", stop_f_ghz, 0.001, 0.1, 100.0,
                        move |v| set_stopband(&mut spec, Some(v * 1e9), None),
                    )}
                    {num_field(
                        "Stopband rejection (dB)", stop_rej_db, 0.5, 5.0, 90.0,
                        move |v| set_stopband(&mut spec, None, Some(v)),
                    )}
                }
                // live realizability
                div { class: "stats", style: "margin-top:12px",
                    div { class: "stat",
                        div { class: "v", "{achieved_rl:.2} dB" }
                        div { class: "l", "achieved in-band RL" }
                    }
                    if !lowpass {
                        div { class: "stat",
                            div {
                                class: "v",
                                style: if lumped_ok { "color:#2dd4bf" } else { "color:#e35d6a" },
                                if lumped_ok { "yes" } else { "no" }
                            }
                            div { class: "l", "lumped realizable" }
                        }
                    }
                }
                if let Some(err) = active_dim_err {
                    div { class: "note honest",
                        "{active_topo} geometry note: " span { class: "mono", "{err}" }
                        if lowpass {
                            " — adjust the cutoff / order."
                        } else {
                            " — the Lumped LC technique may still realize this spec."
                        }
                    }
                } else if lowpass {
                    div { class: "note",
                        "The stepped-impedance low-pass realizes this spec. Edits flow straight into "
                        "synthesis — watch the PASS/FAIL chip and the downstream stages update live."
                    }
                } else {
                    div { class: "note",
                        "Both techniques realize this spec. Edits flow straight into synthesis — "
                        "watch the PASS/FAIL chip and the downstream stages update live."
                    }
                }
            }
        }
    }
}

/// A labelled numeric input row that parses on every keystroke and calls
/// `on_set(parsed)` with the clamped value. Unparseable / out-of-range input is
/// ignored (the field keeps the last good value on re-render).
fn num_field(
    label: &'static str,
    value: f64,
    step: f64,
    min: f64,
    max: f64,
    mut on_set: impl FnMut(f64) + 'static,
) -> Element {
    rsx! {
        div { class: "field",
            span { class: "name", "{label}" }
            input {
                class: "spec-input",
                r#type: "number",
                step: "{step}",
                min: "{min}",
                max: "{max}",
                value: "{value}",
                oninput: move |e| {
                    if let Ok(v) = e.value().parse::<f64>()
                        && v.is_finite()
                        && (min..=max).contains(&v)
                    {
                        on_set(v);
                    }
                },
            }
        }
    }
}

/// A labelled integer input row (clamped to `[min, max]`).
fn int_field(
    label: &'static str,
    value: usize,
    min: usize,
    max: usize,
    mut on_set: impl FnMut(usize) + 'static,
) -> Element {
    rsx! {
        div { class: "field",
            span { class: "name", "{label}" }
            input {
                class: "spec-input",
                r#type: "number",
                step: "1",
                min: "{min}",
                max: "{max}",
                value: "{value}",
                oninput: move |e| {
                    if let Ok(v) = e.value().parse::<usize>()
                        && (min..=max).contains(&v)
                    {
                        on_set(v);
                    }
                },
            }
        }
    }
}

/// Update the (single) stopband point in the spec, overriding the frequency
/// and/or rejection while keeping the other component.
fn set_stopband(spec: &mut Signal<yee_filter::FilterSpec>, f_hz: Option<f64>, rej_db: Option<f64>) {
    let mut w = spec.write();
    let (cur_f, cur_r) = w.mask.stopband.first().copied().unwrap_or((2.4e9, 40.0));
    w.mask.stopband = vec![(f_hz.unwrap_or(cur_f), rej_db.unwrap_or(cur_r))];
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

/// Whether a [`RealizationTechnique`] is realizable by the studio today
/// ("live"), and — when not — the nearest live technique to fall back to.
///
/// This is the **UI-side** live/Soon knowledge (the `yee-filter` engine stays
/// pure-domain). Edge-coupled and Lumped LC are the two live flows
/// ([`Topology`]); the four distributed-resonator techniques are roadmap
/// placeholders whose nearest live realization is the one whose flow can stand
/// in (edge-coupled for the distributed resonators).
fn technique_status(t: RealizationTechnique) -> TechStatus {
    match t {
        RealizationTechnique::EdgeCoupled => TechStatus::Live(Topology::EdgeCoupled),
        RealizationTechnique::Hairpin => TechStatus::Live(Topology::Hairpin),
        RealizationTechnique::LumpedLc => TechStatus::Live(Topology::LumpedLc),
        // Remaining coupled-resonator distributed techniques: the live
        // edge-coupled flow is the nearest stand-in.
        RealizationTechnique::Combline | RealizationTechnique::Interdigital => {
            TechStatus::Soon(Topology::EdgeCoupled)
        }
        // The stepped-impedance low-pass flow is live (ADR-0139): it routes
        // straight into the real F1.2.3 dimensioner + the low-pass response.
        RealizationTechnique::SteppedImpedance => TechStatus::Live(Topology::SteppedImpedance),
    }
}

/// The studio's live/Soon status for a recommended technique.
#[derive(Clone, Copy, PartialEq, Eq)]
enum TechStatus {
    /// Buildable today; routes straight into this [`Topology`]'s flow.
    Live(Topology),
    /// Roadmapped; the recommendation is shown honestly and this [`Topology`]
    /// is offered as the nearest live stand-in.
    Soon(Topology),
}

/// The short flow label for a [`Topology`] (used in the "nearest live" offer).
fn topology_label(t: Topology) -> &'static str {
    match t {
        Topology::EdgeCoupled => "Edge-coupled",
        Topology::Hairpin => "Hairpin",
        Topology::LumpedLc => "Lumped LC",
        Topology::SteppedImpedance => "Stepped-impedance",
    }
}

/// Plain-language word for a [`Response`] class (used in the honest note when
/// the recommendation is for a response the live synthesis flow cannot yet build).
fn response_word(r: Response) -> &'static str {
    match r {
        Response::Lowpass => "low-pass",
        Response::Highpass => "high-pass",
        Response::Bandpass => "band-pass",
        Response::Bandstop => "band-stop",
    }
}

/// The response class a topology's synthesis flow builds: the stepped-impedance
/// flow is **low-pass**; every other live flow is band-pass. Selecting a
/// technique writes this into the shared spec so the Spec form + the engine
/// re-derive in the right response domain.
fn topology_response(t: Topology) -> Response {
    match t {
        Topology::SteppedImpedance => Response::Lowpass,
        Topology::EdgeCoupled | Topology::Hairpin | Topology::LumpedLc => Response::Bandpass,
    }
}

/// Switch the spec's response to match the selected topology's flow (low-pass
/// for stepped-impedance, band-pass otherwise) so the Spec form + the engine
/// re-derive in the right domain when a technique is picked.
fn set_response_for(mut spec: Signal<FilterSpec>, t: Topology) {
    let r = topology_response(t);
    if spec.read().response != r {
        spec.write().response = r;
    }
}

/// Route into a topology's flow: set the topology signal, align the spec's
/// response class to the flow, and land on Spec (so the user can review / refine
/// the seeded spec before synthesizing).
fn route_into(
    mut topology: Signal<Topology>,
    mut active: Signal<Stage>,
    spec: Signal<FilterSpec>,
    t: Topology,
) {
    set_response_for(spec, t);
    topology.set(t);
    active.set(Stage::Spec);
}

/// Seed the editable [`FilterSpec`] from the guided form inputs, preserving the
/// existing approximation / order / Z0 / mask shape (only the response,
/// frequency, and fractional bandwidth come from the form). When the form has a
/// stopband target it is written as the single mask stopband point.
fn seed_spec_from_form(
    mut spec: Signal<FilterSpec>,
    response: Response,
    f0_hz: f64,
    fbw: f64,
    stopband_ghz: Option<f64>,
    stopband_rej_db: f64,
) {
    let mut w = spec.write();
    w.response = response;
    w.f0_hz = f0_hz;
    w.fbw = fbw.max(1e-4);
    if let Some(f_ghz) = stopband_ghz {
        w.mask.stopband = vec![(f_ghz * 1e9, stopband_rej_db)];
    }
}

/// The guided "recommend-a-technique" panel (App.2.0, ADR-0136) atop the expert
/// gallery — the studio's dual-UI entry. A small form (response, centre/cutoff
/// frequency, fractional bandwidth, optional stopband target) drives the pure
/// [`recommend_technique`] engine; the result block highlights the primary
/// technique, shows the rationale, and lists ranked alternatives. **Live**
/// techniques (edge-coupled, lumped) get a "Use this" that seeds the spec and
/// routes into the flow; **Soon** techniques are labelled honestly and offer the
/// nearest live alternative to proceed with.
#[component]
fn guided_panel(
    topology: Signal<Topology>,
    active: Signal<Stage>,
    spec: Signal<FilterSpec>,
) -> Element {
    // Form state (local to the panel; seeds the real spec only on "Use this").
    // 0 = Lowpass, 1 = Highpass, 2 = Bandpass, 3 = Bandstop.
    let mut response_idx = use_signal(|| 2usize);
    let mut f0_ghz = use_signal(|| 2.4f64);
    let mut fbw_pct = use_signal(|| 5.0f64);
    let mut use_stopband = use_signal(|| false);
    let mut stop_ghz = use_signal(|| 4.0f64);
    let mut stop_rej_db = use_signal(|| 30.0f64);
    // The latest recommendation (None until the user clicks Recommend).
    let mut rec = use_signal(|| None::<TechniqueRecommendation>);

    let response_of = |idx: usize| match idx {
        0 => Response::Lowpass,
        1 => Response::Highpass,
        3 => Response::Bandstop,
        _ => Response::Bandpass,
    };
    let cur_response = response_of(response_idx());
    let is_band = matches!(cur_response, Response::Bandpass | Response::Bandstop);
    // Lowpass/highpass call f0 the "cutoff"; band filters call it "centre".
    let freq_label = if is_band {
        "Centre f0 (GHz)"
    } else {
        "Cutoff (GHz)"
    };

    rsx! {
        div { class: "card guided", style: "margin-bottom:18px",
            h2 { class: "card-title",
                "Guided · recommend a technique"
                span { class: "k", "tell me the requirement → topology + rationale" }
            }
            p { class: "guided-lead",
                "New to filter realization? Describe the requirement and the studio recommends a "
                "technique with a plain-language rationale and ranked alternatives — then routes you "
                "into the flow. (The expert gallery below stays available.)"
            }

            // ---- the requirement form ----------------------------------------
            div { class: "fields guided-form",
                div { class: "field",
                    span { class: "name", "Response" }
                    select {
                        class: "spec-input",
                        style: "width:160px;text-align:left",
                        value: "{response_idx}",
                        onchange: move |e| {
                            if let Ok(v) = e.value().parse::<usize>() {
                                response_idx.set(v);
                            }
                        },
                        option { value: "2", "Bandpass" }
                        option { value: "3", "Bandstop" }
                        option { value: "0", "Lowpass" }
                        option { value: "1", "Highpass" }
                    }
                }
                {num_field(
                    freq_label, f0_ghz(), 0.001, 0.001, 100.0,
                    move |v| f0_ghz.set(v),
                )}
                if is_band {
                    {num_field(
                        "Fractional bandwidth (%)", fbw_pct(), 0.1, 0.01, 80.0,
                        move |v| fbw_pct.set(v),
                    )}
                }
                // optional stopband target
                div { class: "field",
                    span { class: "name", "Stopband target" }
                    button {
                        class: if use_stopband() { "seg-btn on" } else { "seg-btn" },
                        style: "border:1px solid var(--border-3);border-radius:5px",
                        onclick: move |_| {
                            let n = !use_stopband();
                            use_stopband.set(n);
                        },
                        if use_stopband() { "on" } else { "off (optional)" }
                    }
                }
                if use_stopband() {
                    {num_field(
                        "Stopband f (GHz)", stop_ghz(), 0.001, 0.001, 100.0,
                        move |v| stop_ghz.set(v),
                    )}
                    {num_field(
                        "Stopband rejection (dB)", stop_rej_db(), 0.5, 5.0, 90.0,
                        move |v| stop_rej_db.set(v),
                    )}
                }
            }

            div { class: "export-row",
                button {
                    class: "btn dl",
                    onclick: move |_| {
                        let s = FilterSpec {
                            response: response_of(response_idx()),
                            approximation: Approximation::Chebyshev { ripple_db: 0.5 },
                            f0_hz: f0_ghz() * 1e9,
                            fbw: (fbw_pct() / 100.0).max(1e-4),
                            order: Some(5),
                            z0_ohm: 50.0,
                            mask: SpecMask {
                                passband_ripple_db: 0.5,
                                return_loss_db: 10.0,
                                stopband: if use_stopband() {
                                    vec![(stop_ghz() * 1e9, stop_rej_db())]
                                } else {
                                    vec![]
                                },
                            },
                        };
                        rec.set(Some(recommend_technique(&s)));
                    },
                    "✦ Recommend a technique"
                }
            }

            // ---- the recommendation result -----------------------------------
            if let Some(r) = rec() {
                {render_recommendation(&r, topology, active, spec, cur_response, f0_ghz(), fbw_pct(),
                    if use_stopband() { Some(stop_ghz()) } else { None }, stop_rej_db())}
            }
        }
    }
}

/// Render a [`TechniqueRecommendation`]: the highlighted primary (with its
/// live/Soon status + a route-into action), the rationale, and the ranked
/// alternatives.
#[allow(clippy::too_many_arguments)]
fn render_recommendation(
    r: &TechniqueRecommendation,
    topology: Signal<Topology>,
    active: Signal<Stage>,
    spec: Signal<FilterSpec>,
    response: Response,
    f0_ghz: f64,
    fbw_pct: f64,
    stopband_ghz: Option<f64>,
    stop_rej_db: f64,
) -> Element {
    let primary = r.primary;
    let status = technique_status(primary);
    let seed = move || {
        seed_spec_from_form(
            spec,
            response,
            f0_ghz * 1e9,
            fbw_pct / 100.0,
            stopband_ghz,
            stop_rej_db,
        );
    };
    // Alternatives, cloned for the move into rsx.
    let alts: Vec<(RealizationTechnique, String)> = r.alternatives.clone();

    rsx! {
        div { class: "rec-result",
            // ---- primary -------------------------------------------------
            div { class: "rec-primary",
                div { class: "rec-head",
                    span { class: "rec-badge", "Recommended" }
                    span { class: "rec-name", "{primary.name()}" }
                    match status {
                        TechStatus::Live(_) => rsx! { span { class: "chip pass", "live" } },
                        TechStatus::Soon(_) => rsx! { span { class: "chip fail", "soon" } },
                    }
                }
                p { class: "rec-rationale", "{r.rationale}" }
                // The studio has two live synthesis domains: the coupled-resonator
                // band-pass flow (edge-coupled / hairpin / lumped) and the
                // stepped-impedance LOW-PASS flow (ADR-0139). A Live technique
                // routes when its flow builds the recommended response — band-pass
                // for the resonator flows, low-pass for stepped-impedance.
                // High-pass / band-stop have no live flow yet, so the
                // recommendation is shown as honest guidance instead.
                {
                    // Route only when SOME live flow builds the recommended
                    // response: the primary itself (Live + matching response), or
                    // a roadmapped (Soon) technique whose nearest-live stand-in
                    // builds that response. Keyed on the stand-in's flow response,
                    // not the recommendation's class, so it generalizes to the
                    // low-pass flow (and any future response) without a band-pass
                    // assumption — high-pass / band-stop, which have no live flow,
                    // fall through to the honest note.
                    let live_match = matches!(status, TechStatus::Live(t) if topology_response(t) == response);
                    let soon_standin = matches!(status, TechStatus::Soon(t) if topology_response(t) == response);
                    if live_match {
                        let TechStatus::Live(t) = status else { unreachable!() };
                        rsx! {
                            button {
                                class: "btn dl rec-use",
                                onclick: move |_| {
                                    seed();
                                    route_into(topology, active, spec, t);
                                },
                                "Use this → seed the spec + open the {topology_label(t)} flow"
                            }
                        }
                    } else if soon_standin {
                        // The primary is roadmapped, but its nearest-live stand-in
                        // builds the recommended response: offer it.
                        let TechStatus::Soon(t) = status else { unreachable!() };
                        rsx! {
                            div { class: "note honest",
                                "This technique is roadmapped (not yet buildable in the studio). "
                                "The nearest live realization is "
                                b { "{topology_label(t)}" } " — proceed with it:"
                            }
                            button {
                                class: "btn dl rec-use",
                                onclick: move |_| {
                                    seed();
                                    route_into(topology, active, spec, t);
                                },
                                "Proceed with {topology_label(t)} (seed the spec + open the flow)"
                            }
                        }
                    } else {
                        rsx! {
                            div { class: "note honest",
                                "The studio's interactive synthesis flow currently builds "
                                b { "band-pass" } " and " b { "stepped-impedance low-pass" } " filters; "
                                "{response_word(response)} synthesis is roadmapped. The recommendation "
                                "above is the technique to target when that flow lands."
                            }
                        }
                    }
                }
            }

            // ---- alternatives --------------------------------------------
            if !alts.is_empty() {
                div { class: "rec-alts",
                    p { class: "lab", "Alternatives" }
                    for (i, (alt, note)) in alts.into_iter().enumerate() {
                        {
                            let alt_status = technique_status(alt);
                            let live = matches!(alt_status, TechStatus::Live(_));
                            rsx! {
                                div { key: "alt{i}", class: "rec-alt",
                                    div { class: "rec-alt-head",
                                        span { class: "rec-alt-name", "{alt.name()}" }
                                        if live {
                                            span { class: "chip pass", "live" }
                                        } else {
                                            span { class: "chip fail", "soon" }
                                        }
                                    }
                                    p { class: "rec-alt-note", "{note}" }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/// The live [`Topology`] a [`RealizationTechnique`] resolves to, regardless of
/// whether it is live or roadmapped (both [`TechStatus`] arms carry one): used by
/// the Compare panel to mark the recommended row and to route "Use this" into the
/// right flow. Every technique `compare_techniques` emits a row for is itself
/// live ([`TechStatus::Live`]); this also resolves a roadmapped *recommendation*
/// to its nearest-live stand-in so the recommended row can still be marked.
fn technique_topology(t: RealizationTechnique) -> Topology {
    match technique_status(t) {
        TechStatus::Live(topo) | TechStatus::Soon(topo) => topo,
    }
}

/// The Compare panel (App.2.5, ADR-0142) on the Technique stage: for the current
/// spec, a side-by-side table over every **live** technique that realizes the
/// spec's response class (the pure [`compare_techniques`] helper). Each row shows
/// the technique, its board size, a PASS / FAIL / "not realizable" chip, and the
/// key graded metrics (worst ripple / return loss / stopband rejection) — all
/// real engine output. The recommended technique (from [`recommend_technique`],
/// resolved to its live [`Topology`]) is marked, and each realizable row offers a
/// "Use this" that routes into that technique's flow ([`route_into`]). A single
/// row (low-pass) still renders; an empty set (high-pass) shows an honest note.
#[component]
fn compare_panel(
    topology: Signal<Topology>,
    active: Signal<Stage>,
    spec: Signal<FilterSpec>,
) -> Element {
    let cur_spec = spec();
    let rows = compare_techniques(&cur_spec);
    // The recommended technique, resolved to the live topology its flow routes
    // into, so the matching Compare row can be marked.
    let recommended_topo = technique_topology(recommend_technique(&cur_spec).primary);

    rsx! {
        div { class: "card compare", style: "margin-top:18px",
            h2 { class: "card-title",
                "Compare · all techniques for this spec"
                span { class: "k", "side-by-side, real engine output" }
            }
            p { class: "guided-lead",
                "Every live technique that realizes this response class, synthesized for the "
                "current spec and graded side-by-side — board size, verdict, and key metrics. "
                "The recommended technique is marked; use any realizable one to route into its flow."
            }
            if rows.is_empty() {
                div { class: "note honest",
                    "No live technique realizes a "
                    b { "{response_word(cur_spec.response)}" }
                    " response yet — that synthesis flow is roadmapped. The guided recommender "
                    "above names the technique to target when it lands."
                }
            } else {
                table { class: "compare-table",
                    thead {
                        tr {
                            th { "Technique" }
                            th { "Board (W×H mm)" }
                            th { "Verdict" }
                            th { "Worst ripple" }
                            th { "Worst RL" }
                            th { "Stopband rej." }
                            th { "" }
                        }
                    }
                    tbody {
                        for row in rows.iter().copied() {
                            {compare_row(row, recommended_topo, topology, active, spec)}
                        }
                    }
                }
            }
        }
    }
}

/// One Compare-table row: the technique's metrics + a recommended marker + a
/// "Use this" route-into action for realizable techniques.
fn compare_row(
    row: TechniqueComparison,
    recommended_topo: Topology,
    topology: Signal<Topology>,
    active: Signal<Stage>,
    spec: Signal<FilterSpec>,
) -> Element {
    let row_topo = technique_topology(row.technique);
    let is_recommended = row_topo == recommended_topo;
    let board_disp = if row.realizable {
        format!("{:.1} × {:.1}", row.board_w_mm, row.board_h_mm)
    } else {
        "—".to_string()
    };
    let ripple_disp = if row.realizable {
        format!("{:.3} dB", row.worst_passband_ripple_db)
    } else {
        "—".to_string()
    };
    let rl_disp = if row.realizable {
        format!("{:.2} dB", row.worst_return_loss_db)
    } else {
        "—".to_string()
    };
    let rej_disp = match row.worst_stopband_rej_db {
        Some(rej) => format!("{rej:.2} dB"),
        None => "—".to_string(),
    };
    let row_cls = if is_recommended { "compare-rec" } else { "" };

    rsx! {
        tr { key: "{row.technique.name()}", class: "{row_cls}",
            td {
                span { class: "rec-alt-name", "{row.technique.name()}" }
                if is_recommended {
                    span { class: "chip pass", style: "margin-left:8px", "recommended" }
                }
            }
            td { "{board_disp}" }
            td {
                match row.pass {
                    Some(true) => rsx! { span { class: "chip pass", "PASS" } },
                    Some(false) => rsx! { span { class: "chip fail", "FAIL" } },
                    None => rsx! { span { class: "chip muted", "not realizable" } },
                }
            }
            td { "{ripple_disp}" }
            td { "{rl_disp}" }
            td { "{rej_disp}" }
            td {
                if row.realizable {
                    button {
                        class: "btn dl",
                        onclick: move |_| route_into(topology, active, spec, row_topo),
                        "Use this"
                    }
                }
            }
        }
    }
}

/// Technique stage: a **guided** recommender panel atop the **expert** topology
/// gallery (the dual-UI entry, App.2.0), with a **Compare** panel (App.2.5) below
/// that synthesizes every live technique for the current spec side-by-side.
/// **Edge-coupled** and **Lumped LC** are live and selectable (clicking a live
/// card routes the downstream stages via `topology`); the rest stay greyed
/// "Soon". The currently selected topology is highlighted.
pub fn technique_stage(
    mut topology: Signal<Topology>,
    mut active: Signal<Stage>,
    spec: Signal<FilterSpec>,
) -> Element {
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
            glyph: r##"<svg viewBox="0 0 120 54"><g stroke="#2dd4bf" stroke-width="4" fill="none"><path d="M14,44 L14,14 L30,14 L30,44"/><path d="M44,44 L44,14 L60,14 L60,44"/><path d="M74,44 L74,14 L90,14 L90,44"/></g></svg>"##,
            selects: Some(Topology::Hairpin),
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
            desc: "hi/lo-Z low-pass sections · F1.2.3",
            glyph: r##"<svg viewBox="0 0 120 54"><g stroke="#e6b24d" stroke-width="4" fill="none"><line x1="10" y1="30" x2="40" y2="30"/><rect x="40" y="20" width="16" height="20"/><line x1="56" y1="30" x2="72" y2="30"/><rect x="72" y="22" width="10" height="16"/><line x1="82" y1="30" x2="110" y2="30"/></g></svg>"##,
            selects: Some(Topology::SteppedImpedance),
        },
    ];
    let cur = topology();
    rsx! {
        div { class: "canvas-head",
            h1 { "Technique" }
            p { class: "sub", "Two ways in: the guided recommender (tell it the requirement, it picks a topology) or the expert gallery below. Each topology realizes the same prototype differently — Edge-coupled (distributed) and Lumped LC (discrete parts + BOM + tolerance) are live; the rest are roadmap placeholders." }
        }

        // ---- guided dual-UI entry (App.2.0) -------------------------------
        guided_panel { topology, active, spec }

        h2 { class: "card-title", style: "margin:18px 0 10px",
            "Expert gallery"
            span { class: "k", "pick a topology directly" }
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
                                    // Align the spec's response class to the flow
                                    // (low-pass for stepped-impedance, band-pass
                                    // otherwise) so the engine re-derives correctly.
                                    set_response_for(spec, t);
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

        // ---- compare-techniques panel (App.2.5) ---------------------------
        compare_panel { topology, active, spec }
    }
}

/// Verify stage (App.2.4, ADR-0141): the active flow's **real** circuit-level
/// verification metrics + an honest statement of what level was verified.
///
/// Mirrors the App.2.3 TopBar pattern — the pure [`verify_view`] helper pulls
/// the active [`Topology`]'s already-computed graded metrics (`MaskReport` /
/// `MaskVerdict` / the stepped low-pass fields). The three stat cards (worst
/// passband ripple, worst in-band return loss, worst stopband rejection) show
/// real values; a metric reads "—" only when it is genuinely absent (no
/// stopband point in the mask). A PASS / FAIL / "not realizable" chip and a
/// level label state whether the realized LC ladder or the synthesized ideal
/// response was graded. The honest note frames full-wave EM of the physical
/// board as a separate native step (the deferred ADR-0133 frontier), not run in
/// the browser — no fabricated EM numbers.
pub fn verify_stage(
    topology: ReadOnlySignal<Topology>,
    designed: ReadOnlySignal<Designed>,
    lumped: ReadOnlySignal<Option<LumpedDesigned>>,
    stepped: ReadOnlySignal<SteppedLowpassDesigned>,
) -> Element {
    // Bind each signal guard to a named local so the borrows passed to
    // `verify_view` live for the whole call (the App.2.3 TopBar idiom).
    let designed_ref = designed.read();
    let lumped_ref = lumped.read();
    let stepped_ref = stepped.read();
    let view = verify_view(topology(), &designed_ref, lumped_ref.as_ref(), &stepped_ref);

    let level_label = match view.level {
        VerifyLevel::RealizedLadder => "Realized LC ladder, graded vs the mask",
        VerifyLevel::SynthesizedIdeal => "Synthesized ideal response vs the mask",
    };
    let ripple_disp = format!("{:.3} dB", view.worst_passband_ripple_db);
    let rl_disp = format!("{:.2} dB", view.worst_return_loss_db);
    let rej_disp = match view.worst_stopband_rej_db {
        Some(rej) => format!("{rej:.2} dB"),
        None => "—".to_string(),
    };

    rsx! {
        div { class: "canvas-head",
            h1 { "Verify" }
            p { class: "sub", "The active design graded against its spec mask — the real metrics the engine already computes, with an honest statement of what level was verified." }
        }
        div { class: "card",
            h2 { class: "card-title",
                "{level_label}"
                match view.pass {
                    Some(true) => rsx! {
                        span { class: "chip pass", style: "margin-left:auto", "PASS" }
                    },
                    Some(false) => rsx! {
                        span { class: "chip fail", style: "margin-left:auto", "FAIL" }
                    },
                    None => rsx! {
                        span { class: "chip muted", style: "margin-left:auto", "not realizable" }
                    },
                }
            }
            div { class: "grid-3",
                div { class: "stat", div { class: "v", "{ripple_disp}" } div { class: "l", "worst passband ripple" } }
                div { class: "stat", div { class: "v", "{rl_disp}" } div { class: "l", "worst in-band return loss" } }
                div { class: "stat", div { class: "v", "{rej_disp}" } div { class: "l", "worst stopband rejection" } }
            }
            div { class: "note honest",
                "These are circuit / synthesis-level metrics: the lumped flow grades its "
                "realized LC ladder, the distributed and low-pass flows grade the synthesized "
                "ideal response. Full-wave EM verification of the physical board (metal "
                "thickness, loss, dispersion, coupling) is a separate native step — the "
                "deferred research frontier — and is not run in the browser. The studio never "
                "hides the ideal-vs-realized gap."
            }
        }
    }
}

/// Trigger a client-side file download of `contents` as `filename` with the
/// given MIME `mime` type. Builds a `Blob` URL in the browser and clicks a
/// synthetic anchor, then revokes the URL. WASM-safe (`dioxus::document::eval`
/// + a small JS shim; no native dep).
fn download_file(filename: &str, mime: &str, contents: &str) {
    // Pass the payload as JSON to the eval so the (possibly multi-line) file
    // body never has to be escaped into the script source by hand.
    let payload = serde_json::json!({
        "name": filename,
        "mime": mime,
        "body": contents,
    });
    let eval = document::eval(
        r#"
        const { name, mime, body } = await dioxus.recv();
        const blob = new Blob([body], { type: mime });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = name;
        document.body.appendChild(a);
        a.click();
        a.remove();
        URL.revokeObjectURL(url);
        "#,
    );
    // Best-effort: a failed send (no JS runtime, e.g. SSR) is a no-op.
    let _ = eval.send(payload);
}

/// A download button wired to [`download_file`]. `make` is invoked on click to
/// produce `(filename, mime, contents)` lazily (so the file is only emitted
/// when the user actually asks for it).
#[component]
fn download_btn(label: String, make: EventHandler<()>) -> Element {
    rsx! {
        button {
            class: "btn dl",
            onclick: move |_| make.call(()),
            "⤓ {label}"
        }
    }
}

/// Export stage: a real, topology-aware parameter sheet + working client-side
/// downloads. The distributed flow emits Gerber (F.Cu + Edge.Cuts) and a KiCad
/// `.kicad_pcb` from the **real** dimensioned layout via `yee-export`, plus a
/// parameter sheet; the lumped flow emits a BOM CSV + a netlist-style parameter
/// sheet, plus Gerber/KiCad from the placed SMD board. Each file is generated
/// on demand from live engine output — nothing is canned.
pub fn export_stage(
    topology: ReadOnlySignal<Topology>,
    designed: ReadOnlySignal<Designed>,
    lumped: ReadOnlySignal<Option<LumpedDesigned>>,
    stepped: ReadOnlySignal<SteppedLowpassDesigned>,
) -> Element {
    match topology() {
        Topology::LumpedLc => export_lumped(lumped),
        Topology::SteppedImpedance => export_stepped(stepped),
        Topology::EdgeCoupled | Topology::Hairpin => export_distributed(designed),
    }
}

/// Distributed (edge-coupled) export panel: real Gerber/KiCad from the layout +
/// a parameter sheet.
fn export_distributed(designed: ReadOnlySignal<Designed>) -> Element {
    let d = designed.read();
    let (bw, bh) = d.board_size_mm;
    let approx = approx_label(&d.spec.approximation);
    let realizable = d.layout.is_some();
    let topo = d.topology_name();

    rsx! {
        div { class: "canvas-head",
            h1 { "Export" }
            p { class: "sub", "The final parameter sheet + manufacturable files, generated live from the dimensioned {topo} layout — Gerber and KiCad are written client-side by the shipped `yee-export` emitters." }
        }
        div { class: "card",
            h2 { class: "card-title",
                "Design summary"
                span { class: "k", "{topo} microstrip" }
            }
            div { class: "fields",
                div { class: "field", span { class: "name", "Topology" } span { class: "val", "{topo} · N={d.order()}" } }
                div { class: "field", span { class: "name", "Approximation" } span { class: "val", "{approx}" } }
                div { class: "field", span { class: "name", "f0 / FBW" } span { class: "val", "{d.spec.f0_hz/1e9:.3} GHz / {d.spec.fbw*100.0:.0}%" } }
                div { class: "field", span { class: "name", "System Z0" } span { class: "val", "{d.spec.z0_ohm:.0} Ω" } }
                div { class: "field",
                    span { class: "name", "Substrate" }
                    if let Some(layout) = d.layout.as_ref() {
                        span { class: "val", "FR-4 · εr {layout.substrate.eps_r:.1} · h {layout.substrate.height_m*1e3:.2} mm" }
                    } else {
                        span { class: "val", "FR-4 (geometry not realizable)" }
                    }
                }
                div { class: "field", span { class: "name", "Board" } span { class: "val", "{bw:.1} × {bh:.1} mm" } }
                div { class: "field",
                    span { class: "name", "Spec verdict" }
                    if d.report.pass {
                        span { class: "val", style: "color:#2dd4bf", "PASS" }
                    } else {
                        span { class: "val", style: "color:#e35d6a", "FAIL" }
                    }
                }
            }
            if realizable {
                div { class: "export-row",
                    download_btn {
                        label: "Gerber F.Cu",
                        make: move |_| {
                            if let Some(layout) = designed.read().layout.as_ref() {
                                let g = yee_export::layout_to_gerber(layout, &Default::default());
                                download_file("filter-F_Cu.gbr", "application/vnd.gerber", &g);
                            }
                        },
                    }
                    download_btn {
                        label: "Gerber Edge.Cuts",
                        make: move |_| {
                            if let Some(layout) = designed.read().layout.as_ref() {
                                let g = yee_export::layout_to_gerber_outline(layout, &Default::default());
                                download_file("filter-Edge_Cuts.gbr", "application/vnd.gerber", &g);
                            }
                        },
                    }
                    download_btn {
                        label: "KiCad .kicad_pcb",
                        make: move |_| {
                            if let Some(layout) = designed.read().layout.as_ref() {
                                let k = yee_export::layout_to_kicad_pcb(layout, &Default::default());
                                download_file("filter.kicad_pcb", "application/octet-stream", &k);
                            }
                        },
                    }
                    download_btn {
                        label: "Parameter sheet",
                        make: move |_| {
                            let sheet = distributed_param_sheet(&designed.read());
                            download_file("filter-parameters.txt", "text/plain", &sheet);
                        },
                    }
                }
                div { class: "note honest",
                    "Gerber + KiCad are written by the shipped `yee-export` emitters from the same "
                    "`Layout` the board view draws — single copper layer + Edge.Cuts outline. "
                    "Drill / soldermask / silkscreen and a Touchstone .s2p (post EM-verify) are "
                    "documented follow-ons."
                }
            } else {
                div { class: "note honest",
                    "Geometry is not realizable on FR-4 for the current spec, so the board exporters "
                    "are unavailable — adjust the spec (Spec stage) or switch to the Lumped LC "
                    "technique. The parameter sheet (synthesis-only) is still available:"
                }
                div { class: "export-row",
                    download_btn {
                        label: "Parameter sheet",
                        make: move |_| {
                            let sheet = distributed_param_sheet(&designed.read());
                            download_file("filter-parameters.txt", "text/plain", &sheet);
                        },
                    }
                }
            }
        }
    }
}

/// Lumped-LC export panel: a BOM CSV + a netlist-style parameter sheet + real
/// Gerber/KiCad from the placed SMD board.
fn export_lumped(lumped: ReadOnlySignal<Option<LumpedDesigned>>) -> Element {
    let guard = lumped.read();
    let Some(d) = guard.as_ref() else {
        return rsx! {
            div { class: "canvas-head",
                h1 { "Export" }
                p { class: "sub", "Manufacturable files + the final parameter sheet." }
            }
            div { class: "card",
                h2 { class: "card-title", "Ladder not realizable" }
                div { class: "note honest",
                    "The current spec does not map to a realizable band-pass LC ladder — adjust "
                    "the order / fractional bandwidth on the Spec stage."
                }
            }
        };
    };
    let (bw, bh) = d.board_size_mm;
    let parts = d.bom_e24.total_parts;

    rsx! {
        div { class: "canvas-head",
            h1 { "Export" }
            p { class: "sub", "The final parameter sheet + manufacturable files for the lumped LC realization — a BOM CSV, a netlist-style parameter sheet, and Gerber/KiCad of the placed SMD board, all generated live." }
        }
        div { class: "card",
            h2 { class: "card-title",
                "Design summary"
                span { class: "k", "lumped LC ladder · SMD" }
            }
            div { class: "fields",
                div { class: "field", span { class: "name", "Topology" } span { class: "val", "lumped LC ladder · N={d.order()}" } }
                div { class: "field", span { class: "name", "Components" } span { class: "val", "{parts} parts · 0603 SMD" } }
                div { class: "field", span { class: "name", "E24 yield" } span { class: "val", "{d.yield_e24.yield_pct:.1}%" } }
                div { class: "field", span { class: "name", "E96 yield" } span { class: "val", "{d.yield_e96.yield_pct:.1}%" } }
                div { class: "field", span { class: "name", "Board" } span { class: "val", "{bw:.1} × {bh:.1} mm" } }
                div { class: "field",
                    span { class: "name", "Spec verdict" }
                    if d.verdict.pass {
                        span { class: "val", style: "color:#2dd4bf", "PASS" }
                    } else {
                        span { class: "val", style: "color:#e35d6a", "FAIL" }
                    }
                }
            }
            div { class: "export-row",
                download_btn {
                    label: "BOM (E24, CSV)",
                    make: move |_| {
                        if let Some(d) = lumped.read().as_ref() {
                            let csv = bom_csv(&d.bom_e24);
                            download_file("filter-bom-e24.csv", "text/csv", &csv);
                        }
                    },
                }
                download_btn {
                    label: "BOM (E96, CSV)",
                    make: move |_| {
                        if let Some(d) = lumped.read().as_ref() {
                            let csv = bom_csv(&d.bom_e96);
                            download_file("filter-bom-e96.csv", "text/csv", &csv);
                        }
                    },
                }
                download_btn {
                    label: "Gerber F.Cu",
                    make: move |_| {
                        if let Some(d) = lumped.read().as_ref() {
                            let g = yee_export::layout_to_gerber(&d.board.layout, &Default::default());
                            download_file("filter-lumped-F_Cu.gbr", "application/vnd.gerber", &g);
                        }
                    },
                }
                download_btn {
                    label: "Gerber Edge.Cuts",
                    make: move |_| {
                        if let Some(d) = lumped.read().as_ref() {
                            let g = yee_export::layout_to_gerber_outline(&d.board.layout, &Default::default());
                            download_file("filter-lumped-Edge_Cuts.gbr", "application/vnd.gerber", &g);
                        }
                    },
                }
                download_btn {
                    label: "KiCad .kicad_pcb",
                    make: move |_| {
                        if let Some(d) = lumped.read().as_ref() {
                            let k = yee_export::layout_to_kicad_pcb(&d.board.layout, &Default::default());
                            download_file("filter-lumped.kicad_pcb", "application/octet-stream", &k);
                        }
                    },
                }
                download_btn {
                    label: "Parameter sheet",
                    make: move |_| {
                        if let Some(d) = lumped.read().as_ref() {
                            let sheet = lumped_param_sheet(d);
                            download_file("filter-lumped-parameters.txt", "text/plain", &sheet);
                        }
                    },
                }
            }
            div { class: "note honest",
                "The BOM CSV is the grouped E-series selection (the Components stage); the Gerber "
                "(F.Cu + Edge.Cuts) + KiCad come from the placed SMD `Layout` (the Layout stage) via "
                "the shipped `yee-export` emitters. Footprint pad geometry + a parasitic-aware land "
                "library are documented follow-ons (F2.2b)."
            }
        }
    }
}

/// Human-readable approximation label (e.g. `"Chebyshev 0.5 dB"`).
fn approx_label(a: &yee_filter::Approximation) -> String {
    match a {
        yee_filter::Approximation::Chebyshev { ripple_db } => {
            format!("Chebyshev {ripple_db:.2} dB")
        }
        yee_filter::Approximation::Butterworth => "Butterworth".to_string(),
    }
}

/// Build the distributed (edge-coupled) parameter sheet: the spec, the
/// synthesized prototype + coupling matrix, and the realized per-resonator
/// geometry — all live engine values.
fn distributed_param_sheet(d: &Designed) -> String {
    let mut s = String::new();
    s.push_str(&format!(
        "# Yee Filter Studio — {} microstrip parameter sheet\n\n",
        d.topology_name()
    ));
    s.push_str("## Specification\n");
    s.push_str("response          : Bandpass\n");
    s.push_str(&format!(
        "approximation     : {}\n",
        approx_label(&d.spec.approximation)
    ));
    s.push_str(&format!("order N           : {}\n", d.order()));
    s.push_str(&format!(
        "f0                : {:.6} GHz\n",
        d.spec.f0_hz / 1e9
    ));
    s.push_str(&format!(
        "fractional BW     : {:.3} %\n",
        d.spec.fbw * 100.0
    ));
    s.push_str(&format!("system Z0         : {:.1} ohm\n", d.spec.z0_ohm));
    s.push_str(&format!(
        "spec mask         : ripple <= {:.3} dB, RL >= {:.1} dB\n",
        d.spec.mask.passband_ripple_db, d.spec.mask.return_loss_db
    ));
    for (f, r) in &d.spec.mask.stopband {
        s.push_str(&format!(
            "  stopband        : {:.4} GHz >= {:.1} dB\n",
            f / 1e9,
            r
        ));
    }
    s.push_str(&format!(
        "\nverdict           : {} (worst RL {:.2} dB, worst ripple {:.3} dB)\n",
        if d.report.pass { "PASS" } else { "FAIL" },
        d.report.worst_return_loss_db,
        d.report.worst_passband_ripple_db
    ));

    s.push_str("\n## Prototype g-values\n");
    for (i, g) in d.g_values.iter().enumerate() {
        s.push_str(&format!("g{i:<2} = {g:.6}\n"));
    }
    s.push_str(&format!(
        "\nQe(in) = {:.4}   Qe(out) = {:.4}\n",
        d.coupling.qe_in, d.coupling.qe_out
    ));

    s.push_str("\n## Coupling matrix M\n");
    for row in &d.coupling.m {
        let cells: Vec<String> = row.iter().map(|v| format!("{v:+.5}")).collect();
        s.push_str(&format!("{}\n", cells.join("  ")));
    }

    if !d.resonators.is_empty() {
        let (bw, bh) = d.board_size_mm;
        s.push_str(&format!(
            "\n## Realized geometry (FR-4, board {bw:.2} x {bh:.2} mm)\n"
        ));
        s.push_str("id  W(mm)   L(mm)   gap(mm)  Z0e(ohm) Z0o(ohm)  k_target k_real\n");
        for r in &d.resonators {
            s.push_str(&format!(
                "R{:<2} {:<7.3} {:<7.2} {:<8} {:<8} {:<9} {:<8} {}\n",
                r.id,
                r.width_mm,
                r.length_mm,
                r.gap_to_next_mm
                    .map(|g| format!("{g:.3}"))
                    .unwrap_or_else(|| "-".into()),
                r.z0e_ohm
                    .map(|v| format!("{v:.1}"))
                    .unwrap_or_else(|| "-".into()),
                r.z0o_ohm
                    .map(|v| format!("{v:.1}"))
                    .unwrap_or_else(|| "-".into()),
                r.target_k
                    .map(|v| format!("{v:.4}"))
                    .unwrap_or_else(|| "-".into()),
                r.realized_k
                    .map(|v| format!("{v:.4}"))
                    .unwrap_or_else(|| "-".into()),
            ));
        }
    } else if let Some(err) = &d.dim_error {
        s.push_str(&format!(
            "\n## Realized geometry\nNOT REALIZABLE on FR-4: {err}\n"
        ));
    }
    s
}

/// Build the lumped-LC parameter sheet: the spec, the ideal LC ladder
/// (SPICE-ish netlist), and the yield summary — all live engine values.
fn lumped_param_sheet(d: &LumpedDesigned) -> String {
    let mut s = String::new();
    s.push_str("# Yee Filter Studio — lumped LC parameter sheet\n\n");
    s.push_str("## Ideal LC ladder (tuned to f0)\n");
    s.push_str("index  branch  L(nH)     C(pF)\n");
    for r in &d.resonators {
        s.push_str(&format!(
            "{:<6} {:<7} {:<9.4} {:<.4}\n",
            r.index,
            if r.is_series { "series" } else { "shunt" },
            r.l_nh,
            r.c_pf
        ));
    }
    s.push_str(&format!(
        "\nverdict : {} (worst RL {:.2} dB, worst rejection {:.1} dB)\n",
        if d.verdict.pass { "PASS" } else { "FAIL" },
        d.verdict.worst_return_loss_db,
        d.verdict.worst_stopband_rej_db
    ));
    s.push_str("\n## Tolerance / yield (Monte-Carlo)\n");
    s.push_str(&format!(
        "E24 (+/-{:.0}%) : {:.1}% yield over {} trials\n",
        d.yield_e24.tolerance_pct, d.yield_e24.yield_pct, d.yield_trials
    ));
    s.push_str(&format!(
        "E96 (+/-{:.0}%) : {:.1}% yield over {} trials\n",
        d.yield_e96.tolerance_pct, d.yield_e96.yield_pct, d.yield_trials
    ));
    let (bw, bh) = d.board_size_mm;
    s.push_str(&format!(
        "\n## Board\n0603 SMD, {bw:.2} x {bh:.2} mm, {} placements\n",
        d.placements.len()
    ));
    s
}

/// Render a [`BomView`] as CSV (ref kind, value, deviation, tolerance, qty).
fn bom_csv(bom: &BomView) -> String {
    let mut s = String::new();
    s.push_str(&format!(
        "# {} BOM (+/-{:.0}%)\n",
        bom.series_name, bom.tolerance_pct
    ));
    s.push_str("kind,ideal,chosen,deviation_pct,tolerance_pct,qty\n");
    for r in &bom.rows {
        s.push_str(&format!(
            "{},{},{},{:.3},{:.1},{}\n",
            r.ref_kind, r.ideal_disp, r.chosen_disp, r.deviation_pct, r.tolerance_pct, r.qty
        ));
    }
    s.push_str(&format!("# total parts: {}\n", bom.total_parts));
    s
}

// ===========================================================================
// REAL lumped-LC stages (App.D.1L) — Synthesis / Components / Tolerance / Layout
// ===========================================================================

/// A "ladder not realizable" placeholder for the lumped stages when the current
/// spec does not map to a band-pass LC ladder.
fn lumped_unrealizable(title: &str) -> Element {
    rsx! {
        div { class: "canvas-head",
            h1 { "{title}" }
            p { class: "sub", "Lumped LC realization of the current spec." }
        }
        div { class: "card",
            h2 { class: "card-title",
                "Ladder not realizable"
                span { class: "chip fail", style: "margin-left:auto", "lumped LC" }
            }
            div { class: "note honest",
                "The current spec does not map to a realizable band-pass LC ladder (e.g. a "
                "degenerate fractional bandwidth or order). Adjust the order / fractional "
                "bandwidth on the Spec stage, or use the distributed Edge-coupled technique."
            }
        }
    }
}

/// Lumped Synthesis stage: the LC ladder resonator table (index, series/shunt,
/// L in nH, C in pF) + the **ideal** `ladder_s21` |S21| vs the spec mask (inline
/// SVG, reusing the response plot) + a PASS/FAIL chip from the realized verdict.
pub fn lumped_synthesis_stage(designed: ReadOnlySignal<Option<LumpedDesigned>>) -> Element {
    let guard = designed.read();
    let Some(d) = guard.as_ref() else {
        return lumped_unrealizable("Synthesis · Lumped LC");
    };
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
    designed: ReadOnlySignal<Option<LumpedDesigned>>,
    mut series_e96: Signal<bool>,
) -> Element {
    let guard = designed.read();
    let Some(d) = guard.as_ref() else {
        return lumped_unrealizable("Components + BOM");
    };
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
pub fn lumped_tolerance_stage(designed: ReadOnlySignal<Option<LumpedDesigned>>) -> Element {
    let guard = designed.read();
    let Some(d) = guard.as_ref() else {
        return lumped_unrealizable("Tolerance / yield");
    };
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
pub fn lumped_layout_stage(designed: ReadOnlySignal<Option<LumpedDesigned>>) -> Element {
    let guard = designed.read();
    let Some(d) = guard.as_ref() else {
        return lumped_unrealizable("Layout · Lumped board");
    };
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

// ===========================================================================
// REAL stage — stepped-impedance low-pass (App.2.2, ADR-0139)
// ===========================================================================

/// Honest placeholder when the stepped-impedance geometry is not realizable on
/// FR-4 for the current spec (mirrors [`lumped_unrealizable`]).
fn stepped_unrealizable(title: &str, msg: &str) -> Element {
    let msg = msg.to_string();
    rsx! {
        div { class: "canvas-head",
            h1 { "{title}" }
            p { class: "sub", "Stepped-impedance low-pass realization of the current spec." }
        }
        div { class: "card",
            h2 { class: "card-title",
                "Geometry not realizable"
                span { class: "chip fail", style: "margin-left:auto", "FR-4" }
            }
            div { class: "note honest",
                "The current low-pass spec does not dimension onto an FR-4 microstrip "
                "(a section width / length goes non-physical). Adjust the cutoff or order on "
                "the Spec stage. Engine note: " span { class: "mono", "{msg}" }
            }
        }
    }
}

/// Stepped-impedance low-pass **Synthesis** stage: the prototype g-values + the
/// swept low-pass `|S21|`/`|S11|` vs the shaded low-pass mask + PASS/FAIL, plus
/// the realized stepped-section table (impedance, electrical length βl, width,
/// length per section). Mirrors [`lumped_synthesis_stage`]. All values are live
/// engine output (the F1.2.3 dimensioner + the App.2.2 low-pass response).
pub fn stepped_synthesis_stage(designed: ReadOnlySignal<SteppedLowpassDesigned>) -> Element {
    let d = designed.read();
    let plot = response_plot(&d.sweep, &d.mask_bands);
    let n = d.order;
    let fc_ghz = d.cutoff_hz() / 1e9;
    let stopband_rows: Vec<(String, String, String, bool)> = d
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
            h1 { "Synthesis · Stepped-impedance low-pass" }
            p { class: "sub", "The low-pass prototype mapped to alternating high-Z / low-Z microstrip sections (Pozar §8.6, shunt-capacitor low-Z first). The swept |S21| is the closed-form low-pass response at Ω = f / f_c. All values are live engine output." }
        }

        // ---- ideal low-pass response vs mask ------------------------------
        div { class: "card", style: "margin-bottom:16px",
            h2 { class: "card-title",
                "Ideal low-pass response vs spec mask"
                span { class: "k", "closed-form |S21| at Ω = f / f_c · f_c = {fc_ghz:.3} GHz" }
                if d.pass {
                    span { class: "chip pass", style: "margin-left:auto", "PASS" }
                } else {
                    span { class: "chip fail", style: "margin-left:auto", "FAIL" }
                }
            }
            div { class: "plot", dangerous_inner_html: "{plot}" }
            div { class: "legend",
                span { span { class: "swatch", style: "background:#2dd4bf" } "|S21| (transmission)" }
                span { span { class: "swatch", style: "background:#6b7480" } "|S11| (reflection)" }
                span { span { class: "swatch", style: "background:#e35d6a" } "forbidden (mask)" }
            }
        }

        div { class: "row",
            // ---- stepped-section table -------------------------------------
            div { class: "card", style: "flex:1.6",
                h2 { class: "card-title",
                    "Line sections"
                    span { class: "k", "N={n} · low-Z first · source → load" }
                }
                table {
                    thead {
                        tr {
                            th { "#" }
                            th { "line" }
                            th { "Z (Ω)" }
                            th { "βl (°)" }
                            th { "W (mm)" }
                            th { "L (mm)" }
                        }
                    }
                    tbody {
                        for s in d.sections.iter() {
                            tr { key: "sec{s.index}",
                                td { class: "mono", "{s.index}" }
                                td {
                                    if s.high_z {
                                        span { class: "pill-sel", "high-Z (series L)" }
                                    } else {
                                        span { class: "pill-sel", style: "background:#11302a;border-color:#1f5138;color:#2dd4bf", "low-Z (shunt C)" }
                                    }
                                }
                                td { class: "mono", "{s.z_ohm:.1}" }
                                td { class: "mono", "{s.betal_deg:.2}" }
                                td { class: "mono", "{s.width_mm:.3}" }
                                td { class: "mono", "{s.length_mm:.2}" }
                            }
                        }
                    }
                }
            }

            // ---- prototype + impedance pair --------------------------------
            div { class: "card", style: "flex:1",
                h2 { class: "card-title", "Prototype + impedances" }
                div { class: "stats",
                    div { class: "stat",
                        div { class: "v", "{d.z_high():.0} Ω" }
                        div { class: "l", "high-Z line" }
                    }
                    div { class: "stat",
                        div { class: "v", "{d.z_low():.0} Ω" }
                        div { class: "l", "low-Z line" }
                    }
                }
                p { class: "lab", style: "margin-top:16px", "g-values" }
                div { class: "gvals",
                    for (i, g) in d.g_values.iter().enumerate() {
                        div { key: "g{i}", class: "gval", b { "g{i} " } "{g:.4}" }
                    }
                }
                div { class: "stats", style: "margin-top:16px",
                    div { class: "stat",
                        div { class: "v", "{d.worst_passband_ripple_db:.3} dB" }
                        div { class: "l", "passband ripple" }
                    }
                    div { class: "stat",
                        div { class: "v", "{d.worst_return_loss_db:.2} dB" }
                        div { class: "l", "in-band return loss" }
                    }
                }
            }
        }

        // ---- stopband verdict ---------------------------------------------
        if !stopband_rows.is_empty() {
            div { class: "card", style: "margin-top:16px",
                h2 { class: "card-title", "Stopband rejection" }
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
        }

        div { class: "note honest",
            "Honest note: this is the IDEAL closed-form low-pass prototype response and the "
            "first-order stepped-impedance dimensions (Pozar §8.6). Junction discontinuities, "
            "step reactances, and the full-wave realized response are EM-verify (later) — the "
            "studio never hides the ideal-vs-realized gap."
        }
    }
}

/// Stepped-impedance low-pass **Layout** stage: the dimensioned in-line board
/// top-view (the generic `Layout` SVG — reuses [`board_svg`]) + the material
/// stackup + the per-section geometry table. Mirrors [`lumped_layout_stage`] /
/// [`layout_stage`].
pub fn stepped_layout_stage(designed: ReadOnlySignal<SteppedLowpassDesigned>) -> Element {
    let d = designed.read();
    let Some(layout) = d.layout.as_ref() else {
        let msg = d
            .dim_error
            .clone()
            .unwrap_or_else(|| "geometry not realizable on FR-4".into());
        return stepped_unrealizable("Layout · Stepped-impedance board", &msg);
    };
    let board = board_svg(layout);
    let sub = &layout.substrate;
    let (bw, bh) = d.board_size_mm;

    rsx! {
        div { class: "canvas-head",
            h1 { "Layout + Materials" }
            p { class: "sub", "The stepped-impedance low-pass line — alternating high-Z / low-Z sections laid end to end with Z0 feeds — and the material stackup, all from the live F1.2.3 dimensional synthesis." }
        }

        div { class: "row",
            // ---- board top view ---------------------------------------------
            div { class: "card", style: "flex:1.6",
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
                    div { class: "editrow", span { "high / low Z" } span { class: "v", "{d.z_high():.0} / {d.z_low():.0} Ω" } }
                }
            }
        }

        // ---- per-section geometry table -----------------------------------
        div { class: "card", style: "margin-top:16px",
            h2 { class: "card-title",
                "Sections · geometry"
                span { class: "k", "source → load · Z / βl / W / L per section" }
            }
            table {
                thead {
                    tr {
                        th { "#" }
                        th { "line" }
                        th { "Z (Ω)" }
                        th { "βl (°)" }
                        th { "W (mm)" }
                        th { "L (mm)" }
                    }
                }
                tbody {
                    for s in d.sections.iter() {
                        tr { key: "lsec{s.index}",
                            td { class: "mono", "{s.index}" }
                            td {
                                if s.high_z {
                                    span { class: "pill-sel", "high-Z" }
                                } else {
                                    span { class: "pill-sel", style: "background:#11302a;border-color:#1f5138;color:#2dd4bf", "low-Z" }
                                }
                            }
                            td { class: "mono", "{s.z_ohm:.1}" }
                            td { class: "mono", "{s.betal_deg:.2}" }
                            td { class: "mono", "{s.width_mm:.3}" }
                            td { class: "mono", "{s.length_mm:.2}" }
                        }
                    }
                }
            }
            p { class: "note honest",
                "Section widths come from Hammerstad-Jensen synthesis at each impedance; the "
                "lengths from βl and the guided wavelength at f_c. The same Layout feeds the "
                "shipped Gerber / KiCad exporters (Export stage)."
            }
        }
    }
}

/// Build the stepped-impedance low-pass parameter sheet: the spec, the prototype
/// g-values, and the realized per-section geometry — all live engine values.
fn stepped_param_sheet(d: &SteppedLowpassDesigned) -> String {
    let mut s = String::new();
    s.push_str("# Yee Filter Studio — stepped-impedance low-pass parameter sheet\n\n");
    s.push_str("## Specification\n");
    s.push_str("response          : Lowpass\n");
    s.push_str(&format!(
        "approximation     : {}\n",
        approx_label(&d.spec.approximation)
    ));
    s.push_str(&format!("order N           : {}\n", d.order));
    s.push_str(&format!(
        "cutoff f_c        : {:.6} GHz\n",
        d.cutoff_hz() / 1e9
    ));
    s.push_str(&format!("system Z0         : {:.1} ohm\n", d.spec.z0_ohm));
    s.push_str(&format!(
        "high / low Z      : {:.1} / {:.1} ohm\n",
        d.z_high(),
        d.z_low()
    ));
    s.push_str(&format!(
        "verdict           : {} (worst RL {:.2} dB, worst ripple {:.3} dB)\n",
        if d.pass { "PASS" } else { "FAIL" },
        d.worst_return_loss_db,
        d.worst_passband_ripple_db
    ));
    for (f, achieved, required, met) in &d.stopband {
        s.push_str(&format!(
            "  stopband        : {:.4} GHz -> {:.1} dB (need {:.1} dB) {}\n",
            f / 1e9,
            achieved,
            required,
            if *met { "MET" } else { "UNDER" }
        ));
    }
    s.push_str("\n## Prototype g-values\n");
    for (i, g) in d.g_values.iter().enumerate() {
        s.push_str(&format!("g{i:<2} = {g:.6}\n"));
    }
    if !d.sections.is_empty() {
        let (bw, bh) = d.board_size_mm;
        s.push_str(&format!(
            "\n## Realized sections (FR-4, board {bw:.2} x {bh:.2} mm, source -> load)\n"
        ));
        s.push_str("#   line    Z(ohm)  betal(deg)  W(mm)    L(mm)\n");
        for sec in &d.sections {
            s.push_str(&format!(
                "{:<3} {:<7} {:<7.1} {:<11.2} {:<8.3} {:<.2}\n",
                sec.index,
                if sec.high_z { "high-Z" } else { "low-Z" },
                sec.z_ohm,
                sec.betal_deg,
                sec.width_mm,
                sec.length_mm
            ));
        }
    } else if let Some(err) = &d.dim_error {
        s.push_str(&format!(
            "\n## Realized geometry\nNOT REALIZABLE on FR-4: {err}\n"
        ));
    }
    s
}

/// Stepped-impedance low-pass **Export** panel: real Gerber / KiCad from the
/// dimensioned `Layout` (the shipped generic `yee-export` emitters) + a
/// parameter sheet. Mirrors [`export_distributed`].
fn export_stepped(designed: ReadOnlySignal<SteppedLowpassDesigned>) -> Element {
    let d = designed.read();
    let (bw, bh) = d.board_size_mm;
    let approx = approx_label(&d.spec.approximation);
    let realizable = d.layout.is_some();

    rsx! {
        div { class: "canvas-head",
            h1 { "Export" }
            p { class: "sub", "The parameter sheet + manufacturable files, generated live from the dimensioned stepped-impedance low-pass layout — Gerber and KiCad are written client-side by the shipped `yee-export` emitters." }
        }
        div { class: "card",
            h2 { class: "card-title",
                "Design summary"
                span { class: "k", "stepped-impedance low-pass microstrip" }
            }
            div { class: "fields",
                div { class: "field", span { class: "name", "Topology" } span { class: "val", "Stepped-impedance · N={d.order}" } }
                div { class: "field", span { class: "name", "Approximation" } span { class: "val", "{approx}" } }
                div { class: "field", span { class: "name", "Cutoff f_c" } span { class: "val", "{d.cutoff_hz()/1e9:.3} GHz" } }
                div { class: "field", span { class: "name", "System Z0" } span { class: "val", "{d.spec.z0_ohm:.0} Ω" } }
                div { class: "field", span { class: "name", "High / low Z" } span { class: "val", "{d.z_high():.0} / {d.z_low():.0} Ω" } }
                div { class: "field", span { class: "name", "Board" } span { class: "val", "{bw:.1} × {bh:.1} mm" } }
                div { class: "field",
                    span { class: "name", "Spec verdict" }
                    if d.pass {
                        span { class: "val", style: "color:#2dd4bf", "PASS" }
                    } else {
                        span { class: "val", style: "color:#e35d6a", "FAIL" }
                    }
                }
            }
            if realizable {
                div { class: "export-row",
                    download_btn {
                        label: "Gerber F.Cu",
                        make: move |_| {
                            if let Some(layout) = designed.read().layout.as_ref() {
                                let g = yee_export::layout_to_gerber(layout, &Default::default());
                                download_file("lowpass-F_Cu.gbr", "application/vnd.gerber", &g);
                            }
                        },
                    }
                    download_btn {
                        label: "Gerber Edge.Cuts",
                        make: move |_| {
                            if let Some(layout) = designed.read().layout.as_ref() {
                                let g = yee_export::layout_to_gerber_outline(layout, &Default::default());
                                download_file("lowpass-Edge_Cuts.gbr", "application/vnd.gerber", &g);
                            }
                        },
                    }
                    download_btn {
                        label: "KiCad .kicad_pcb",
                        make: move |_| {
                            if let Some(layout) = designed.read().layout.as_ref() {
                                let k = yee_export::layout_to_kicad_pcb(layout, &Default::default());
                                download_file("lowpass.kicad_pcb", "application/octet-stream", &k);
                            }
                        },
                    }
                    download_btn {
                        label: "Parameter sheet",
                        make: move |_| {
                            let sheet = stepped_param_sheet(&designed.read());
                            download_file("lowpass-parameters.txt", "text/plain", &sheet);
                        },
                    }
                }
                div { class: "note honest",
                    "Gerber + KiCad are written by the shipped `yee-export` emitters from the same "
                    "`Layout` the board view draws — single copper layer + Edge.Cuts outline. "
                    "Drill / soldermask / silkscreen and a Touchstone .s2p (post EM-verify) are "
                    "documented follow-ons."
                }
            } else {
                div { class: "note honest",
                    "Geometry is not realizable on FR-4 for the current spec, so the board exporters "
                    "are unavailable — adjust the cutoff / order on the Spec stage. The parameter "
                    "sheet (synthesis-only) is still available:"
                }
                div { class: "export-row",
                    download_btn {
                        label: "Parameter sheet",
                        make: move |_| {
                            let sheet = stepped_param_sheet(&designed.read());
                            download_file("lowpass-parameters.txt", "text/plain", &sheet);
                        },
                    }
                }
            }
        }
    }
}
