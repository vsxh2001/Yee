//! # yee-studio-web
//!
//! Yee Filter Studio — **App.D.0 Dioxus (web) proof-of-concept** (ADR-0110).
//!
//! A pure-Rust, web-first rebuild of the filter-design studio's view layer in
//! [Dioxus]. This POC validates three things before the full build is
//! committed: (a) Dioxus delivers SaaS-class polish (real DOM + CSS, the
//! ceiling egui can't reach), (b) the Rust→WASM **engine bridge** works, and
//! (c) it builds + serves on web.
//!
//! It renders **Shell A** — a top bar + a stage rail over a central canvas —
//! and drives the live `yee-synth` / `yee-filter` / `yee-layout` engine on the
//! committed Chebyshev N=5 fixture (see [`engine`]). The **Technique** stage
//! offers two live topologies that route the downstream flow:
//!
//! - **Edge-coupled (distributed)** — the POC's original flow: the real
//!   coupling matrix `M`, g-values / external Q, the swept ideal `|S21|`/`|S11|`
//!   vs the shaded spec mask + PASS/FAIL (**Synthesis**), and the real
//!   dimensioned board top-view + material stackup + per-resonator table
//!   (**Layout + Materials**).
//! - **Lumped LC** — the App.D.1L flow (ADR-0120), all from the shipped F2.x
//!   engine: the synthesized LC **ladder** + ideal `ladder_s21` vs mask
//!   (**Synthesis**), E24/E96 **component selection + BOM** (**Components**),
//!   the Monte-Carlo **yield** analysis (**Tolerance**), and the placed SMD
//!   **board** SVG + placement table (**Layout**).
//!
//! Spec / Verify / Export are styled-but-static stubs (prove the shell, not
//! full interactivity). `StudioState`, the engine crates, and the eframe
//! `yee-studio` are all untouched — this crate is additive.
//!
//! [Dioxus]: https://dioxuslabs.com/

mod engine;
mod stages;
mod svg;

use dioxus::prelude::*;

use engine::{Designed, LumpedDesigned, design_demo, design_lumped};
use stages::{Stage, Topology};

/// The design-system stylesheet (tokens + base components). Embedded via the
/// `asset!` macro so the web build bundles it and `dx`/trunk fingerprints it.
const STUDIO_CSS: Asset = asset!("/assets/studio.css");

fn main() {
    dioxus::launch(App);
}

/// Root component: builds the live design once, holds the active-stage signal,
/// and lays out Shell A (top bar + rail + canvas).
#[component]
fn App() -> Element {
    // Run both real engine pipelines once at startup; every stage reads from
    // the one its active topology selects.
    let designed = use_signal(design_demo);
    let lumped = use_signal(design_lumped);
    let topology = use_signal(|| Topology::EdgeCoupled);
    let active = use_signal(|| Stage::Synthesis);
    // E24 (false) / E96 (true) toggle for the lumped Components + BOM stage.
    let series_e96 = use_signal(|| false);

    rsx! {
        document::Stylesheet { href: STUDIO_CSS }
        div { class: "app",
            TopBar { designed }
            div { class: "body",
                Rail { active, topology }
                main { class: "canvas",
                    StageCanvas {
                        stage: active(),
                        topology,
                        active,
                        designed,
                        lumped,
                        series_e96,
                    }
                }
            }
        }
    }
}

/// Top bar: app brand + the spec summary chip + the live PASS/FAIL verdict chip.
#[component]
fn TopBar(designed: ReadOnlySignal<Designed>) -> Element {
    let d = designed.read();
    let spec = &d.spec;
    let approx = match spec.approximation {
        yee_filter::Approximation::Chebyshev { ripple_db } => {
            format!("Chebyshev {ripple_db:.1} dB")
        }
        yee_filter::Approximation::Butterworth => "Butterworth".to_string(),
    };
    let summary = format!(
        "· {approx} · N={} · {:.2} GHz · {:.0}%",
        d.order(),
        spec.f0_hz / 1e9,
        spec.fbw * 100.0
    );
    let pass = d.report.pass;

    rsx! {
        header { class: "topbar",
            span { class: "dot" }
            span { class: "brand", "Yee Filter Studio" }
            span { class: "spec-chip", "{summary}" }
            span { class: "spacer" }
            if pass {
                span { class: "chip pass", span { class: "dot-sm" } "SPEC MET" }
            } else {
                span { class: "chip fail", span { class: "dot-sm" } "SPEC FAIL" }
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
/// Tolerance / Layout to the lumped renderers when the lumped topology is
/// selected.
#[component]
fn StageCanvas(
    stage: Stage,
    topology: Signal<Topology>,
    active: Signal<Stage>,
    designed: ReadOnlySignal<Designed>,
    lumped: ReadOnlySignal<LumpedDesigned>,
    series_e96: Signal<bool>,
) -> Element {
    let lumped_flow = topology() == Topology::LumpedLc;
    match stage {
        Stage::Spec => stages::spec_stage(designed),
        Stage::Technique => stages::technique_stage(topology, active),
        Stage::Synthesis if lumped_flow => stages::lumped_synthesis_stage(lumped),
        Stage::Synthesis => stages::synthesis_stage(designed),
        Stage::Components => stages::lumped_components_stage(lumped, series_e96),
        Stage::Tolerance => stages::lumped_tolerance_stage(lumped),
        Stage::Layout if lumped_flow => stages::lumped_layout_stage(lumped),
        Stage::Layout => stages::layout_stage(designed),
        Stage::Verify => stages::verify_stage(),
        Stage::Export => stages::export_stage(designed),
    }
}
