//! `yee-studio` binary entry — the App.0 filter-design desktop skeleton.
//!
//! A thin `eframe::run_native` wrapper that seeds a [`StudioState`] from a
//! default *satisfiable* Chebyshev 0.5 dB N=5 bandpass filter (f0 = 2 GHz,
//! FBW = 0.10, return loss 9 dB, stopband (2.4 GHz, 40 dB)) — the committed
//! `cheb_bpf.toml`-equivalent values — and hands it to [`StudioApp`]. See
//! [`yee_studio`] for the headless logic layer and [`yee_studio::app`] for the
//! UI.
//!
//! The windowed binary is build-only in CI (no display); the logic is gated by
//! the headless `studio_state_recompute_*` tests in the library crate.

use eframe::egui;

use yee_filter::{Approximation, FilterSpec, Response, SpecMask};
use yee_studio::{StudioState, app::StudioApp};

/// The default satisfiable Chebyshev 0.5 dB N=5 bandpass spec the app opens to.
fn default_spec() -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: 2.0e9,
        fbw: 0.10,
        order: Some(5),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 9.0,
            stopband: vec![(2.4e9, 40.0)],
        },
    }
}

fn main() -> eframe::Result<()> {
    let state = StudioState::from_spec(default_spec());

    let native_options = eframe::NativeOptions {
        renderer: eframe::Renderer::Wgpu,
        viewport: egui::ViewportBuilder::default()
            .with_title("Yee Filter Studio (App.0)")
            .with_inner_size([1280.0, 800.0]),
        ..Default::default()
    };

    eframe::run_native(
        "Yee Filter Studio",
        native_options,
        Box::new(move |_cc| Ok(Box::new(StudioApp::new(state.clone())))),
    )
}
