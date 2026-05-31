//! # yee-studio-web
//!
//! Yee Filter Studio — the pure-Rust, web-first **Dioxus** filter-design studio
//! (App.D.2, ADR-0130; grew out of the App.D.0 POC, ADR-0110). This is *the*
//! studio: it replaces the retired eframe `yee-studio` view as the goal's
//! polished-UI component.
//!
//! A pure-Rust view layer in [Dioxus] (real DOM + CSS) that drives the live
//! `yee-synth` / `yee-filter` / `yee-layout` / `yee-export` engines entirely
//! client-side in WASM (no native / FDTD / wgpu in the dep graph — ADR-0089).
//!
//! It renders **Shell A** — a top bar + a stage rail over a central canvas. The
//! editable [`yee_filter::FilterSpec`] is the single source of
//! truth: the **Spec** stage edits it and the whole studio re-derives live. The
//! **Technique** stage offers two live topologies that route the flow:
//!
//! - **Edge-coupled (distributed)** — the real coupling matrix `M`, g-values /
//!   external Q, the swept ideal `|S21|`/`|S11|` vs the shaded spec mask +
//!   PASS/FAIL (**Synthesis**), and the real dimensioned board top-view +
//!   material stackup + per-resonator table (**Layout + Materials**).
//! - **Lumped LC** — the App.D.1L flow (ADR-0120), all from the shipped F2.x
//!   engine: the synthesized LC **ladder** + ideal `ladder_s21` vs mask
//!   (**Synthesis**), E24/E96 **component selection + BOM** (**Components**),
//!   the Monte-Carlo **yield** analysis (**Tolerance**), and the placed SMD
//!   **board** SVG + placement table (**Layout**).
//!
//! **Spec** and **Export** are real: Spec is a live editable form driving
//! synthesis; Export emits a parameter sheet, a BOM CSV (lumped), and Gerber /
//! KiCad files from the real layout via the shipped `yee-export` emitters,
//! downloaded client-side. **Verify** shows the active flow's real circuit-level
//! mask metrics (App.2.4 / ADR-0141), honest that full-wave EM of the physical
//! board is a separate native step. The remaining distributed topologies
//! (Combline, …) are honestly labelled "Soon".
//!
//! [Dioxus]: https://dioxuslabs.com/

mod engine;
mod stages;
mod svg;

use dioxus::prelude::*;

use engine::{
    Designed, LumpedDesigned, SteppedLowpassDesigned, demo_spec, design_demo, design_demo_from,
    design_lumped, design_lumped_from, design_stepped, design_stepped_from, topbar_view,
};
use stages::{Stage, Topology};
use yee_filter::FilterSpec;

/// The design-system stylesheet (tokens + base components). Embedded via the
/// `asset!` macro so the web build bundles it and `dx`/trunk fingerprints it.
const STUDIO_CSS: Asset = asset!("/assets/studio.css");

fn main() {
    dioxus::launch(App);
}

/// Root component: holds the editable [`FilterSpec`] as the single source of
/// truth, re-derives both engine pipelines whenever it changes, holds the
/// active-stage signal, and lays out Shell A (top bar + rail + canvas).
///
/// The Spec stage edits `spec`; a [`use_effect`] re-runs `design_*_from(spec)`
/// into `designed` / `lumped` so every stage updates live. `lumped` is an
/// `Option` because some specs are not realizable as a band-pass LC ladder.
#[component]
fn App() -> Element {
    // The editable design intent. Everything else derives from it.
    let spec = use_signal(demo_spec);

    // The selected realization technique (drives the rail + the distributed
    // geometry derivation). Declared before the re-derivation effect so the
    // `designed` memo can depend on it: switching technique re-dimensions the
    // board for the new topology.
    let topology = use_signal(|| Topology::EdgeCoupled);

    // Derived engine output, recomputed reactively on every spec OR topology
    // edit. Seeded from the demo spec (edge-coupled); the `use_effect`
    // re-derives on every subsequent edit. The hairpin and edge-coupled flows
    // share synthesis / response / verdict and differ only in the dimensioned
    // geometry, so re-running `design_demo_from` on a topology change is what
    // swaps the board between the two distributed realizations.
    let mut designed = use_signal(design_demo);
    let mut lumped = use_signal(|| Some(design_lumped()));
    // The stepped-impedance low-pass design (ADR-0139): re-derived on every spec
    // edit, mirroring `lumped`. It always succeeds (the dimensioner degrades to
    // a `dim_error` rather than failing), so it is a plain value, not an Option.
    let mut stepped = use_signal(design_stepped);
    use_effect(move || {
        let s: FilterSpec = spec();
        let t: Topology = topology();
        designed.set(design_demo_from(s.clone(), t));
        lumped.set(design_lumped_from(s.clone()).ok());
        stepped.set(design_stepped_from(s));
    });

    let active = use_signal(|| Stage::Synthesis);
    // E24 (false) / E96 (true) toggle for the lumped Components + BOM stage.
    let series_e96 = use_signal(|| false);

    rsx! {
        document::Stylesheet { href: STUDIO_CSS }
        div { class: "app",
            TopBar { topology, designed, lumped, stepped }
            div { class: "body",
                Rail { active, topology }
                main { class: "canvas",
                    StageCanvas {
                        stage: active(),
                        topology,
                        active,
                        spec,
                        designed,
                        lumped,
                        stepped,
                        series_e96,
                    }
                }
            }
        }
    }
}

/// Top bar: app brand + the **active flow's** spec summary chip + its live
/// PASS/FAIL verdict chip (App.2.3, ADR-0140).
///
/// The summary + verdict are computed by the pure [`topbar_view`] helper, which
/// dispatches on the active [`Topology`]: band-pass (edge-coupled / hairpin) →
/// the distributed verdict; lumped → the lumped ladder verdict; stepped-impedance
/// → the low-pass cutoff + verdict. When the active flow's design is not
/// realizable (`None` verdict — e.g. an unrealizable lumped ladder) a muted
/// "geometry not realizable" chip is shown instead of PASS/FAIL.
#[component]
fn TopBar(
    topology: ReadOnlySignal<Topology>,
    designed: ReadOnlySignal<Designed>,
    lumped: ReadOnlySignal<Option<LumpedDesigned>>,
    stepped: ReadOnlySignal<SteppedLowpassDesigned>,
) -> Element {
    // Bind each signal guard to a named local so the borrows passed to
    // `topbar_view` live for the whole call (no reliance on argument-position
    // temporary-lifetime extension; robust to later refactors).
    let designed_ref = designed.read();
    let lumped_ref = lumped.read();
    let stepped_ref = stepped.read();
    let (summary, verdict) =
        topbar_view(topology(), &designed_ref, lumped_ref.as_ref(), &stepped_ref);

    rsx! {
        header { class: "topbar",
            span { class: "dot" }
            span { class: "brand", "Yee Filter Studio" }
            span { class: "spec-chip", "{summary}" }
            span { class: "spacer" }
            match verdict {
                Some(true) => rsx! {
                    span { class: "chip pass", span { class: "dot-sm" } "SPEC MET" }
                },
                Some(false) => rsx! {
                    span { class: "chip fail", span { class: "dot-sm" } "SPEC FAIL" }
                },
                None => rsx! {
                    span { class: "chip muted", span { class: "dot-sm" } "geometry not realizable" }
                },
            }
        }
    }
}

/// The left stage rail: the active topology's stages with icon + label;
/// clicking switches the central canvas. The rail swaps when the technique
/// changes (the lumped flow adds Components + Tolerance).
#[component]
fn Rail(active: Signal<Stage>, topology: ReadOnlySignal<Topology>) -> Element {
    rsx! {
        nav { class: "rail",
            for stage in Stage::rail(topology()).iter().copied() {
                {
                    let on = active() == stage;
                    rsx! {
                        button {
                            key: "{stage:?}",
                            class: if on { "item on" } else { "item" },
                            onclick: move |_| active.set(stage),
                            span { class: "ic", "{stage.icon()}" }
                            span { class: "lab", "{stage.label()}" }
                        }
                    }
                }
            }
        }
    }
}

/// Dispatch the active stage to its renderer, routing Synthesis / Components /
/// Tolerance / Layout to the lumped renderers for the lumped topology, and
/// Synthesis / Layout to the stepped-impedance low-pass renderers for the
/// `SteppedImpedance` topology (ADR-0139).
#[component]
fn StageCanvas(
    stage: Stage,
    topology: Signal<Topology>,
    active: Signal<Stage>,
    spec: Signal<FilterSpec>,
    designed: ReadOnlySignal<Designed>,
    lumped: ReadOnlySignal<Option<LumpedDesigned>>,
    stepped: ReadOnlySignal<SteppedLowpassDesigned>,
    series_e96: Signal<bool>,
) -> Element {
    let lumped_flow = topology() == Topology::LumpedLc;
    let stepped_flow = topology() == Topology::SteppedImpedance;
    match stage {
        Stage::Spec => stages::spec_stage(spec, topology.into(), designed, lumped, stepped),
        Stage::Technique => stages::technique_stage(topology, active, spec),
        Stage::Synthesis if stepped_flow => stages::stepped_synthesis_stage(stepped),
        Stage::Synthesis if lumped_flow => stages::lumped_synthesis_stage(lumped),
        Stage::Synthesis => stages::synthesis_stage(designed),
        Stage::Components => stages::lumped_components_stage(lumped, series_e96),
        Stage::Tolerance => stages::lumped_tolerance_stage(lumped),
        Stage::Layout if stepped_flow => stages::stepped_layout_stage(stepped),
        Stage::Layout if lumped_flow => stages::lumped_layout_stage(lumped),
        Stage::Layout => stages::layout_stage(designed),
        Stage::Verify => stages::verify_stage(topology.into(), designed, lumped, stepped),
        Stage::Export => stages::export_stage(topology.into(), designed, lumped, stepped),
    }
}
