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
//! It renders **Shell A** — a top bar + a six-stage rail (Spec, Technique,
//! Synthesis, Layout, Verify, Export) over a central canvas — and drives the
//! live `yee-synth` / `yee-filter` / `yee-layout` engine on the committed
//! Chebyshev N=5 fixture (see [`engine`]). Two stages render **real engine
//! output**:
//!
//! - **Synthesis** — the real coupling matrix `M`, g-values, external Q, the
//!   swept ideal `|S21|` / `|S11|` (inline SVG) versus the shaded spec mask,
//!   and the real PASS/FAIL verdict.
//! - **Layout + Materials** — the real dimensioned board top-view (inline SVG
//!   from `dimension_edge_coupled_layout`), the material stackup, and the
//!   per-resonator table (`W`, `L`, gap, `Z0e`/`Z0o`, `εeff`, realized `k`).
//!
//! Spec / Technique / Export are styled-but-static stubs (prove the shell, not
//! full interactivity). `StudioState`, the engine crates, and the eframe
//! `yee-studio` are all untouched — this crate is additive.
//!
//! [Dioxus]: https://dioxuslabs.com/

mod engine;
mod stages;
mod svg;

use dioxus::prelude::*;

use engine::{Designed, design_demo};
use stages::Stage;

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
    // Run the real engine pipeline once at startup; every stage reads from it.
    let designed = use_signal(design_demo);
    let active = use_signal(|| Stage::Synthesis);

    rsx! {
        document::Stylesheet { href: STUDIO_CSS }
        div { class: "app",
            TopBar { designed }
            div { class: "body",
                Rail { active }
                main { class: "canvas",
                    StageCanvas { stage: active(), designed }
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

/// The left stage rail: six stages with icon + label; clicking switches the
/// central canvas.
#[component]
fn Rail(active: Signal<Stage>) -> Element {
    rsx! {
        nav { class: "rail",
            for stage in Stage::ALL {
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

/// Dispatch the active stage to its renderer.
#[component]
fn StageCanvas(stage: Stage, designed: ReadOnlySignal<Designed>) -> Element {
    match stage {
        Stage::Spec => stages::spec_stage(designed),
        Stage::Technique => stages::technique_stage(),
        Stage::Synthesis => stages::synthesis_stage(designed),
        Stage::Layout => stages::layout_stage(designed),
        Stage::Verify => stages::verify_stage(),
        Stage::Export => stages::export_stage(designed),
    }
}
