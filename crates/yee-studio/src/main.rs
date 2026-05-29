//! `yee-studio` binary entry — the App.0 filter-design studio shell.
//!
//! The same eframe [`StudioApp`] runs **natively** (desktop, `eframe::run_native`)
//! and in the **browser** (web, `eframe::WebRunner`), one codebase per ADR-0089.
//! Both seed a [`StudioState`] from a default *satisfiable* Chebyshev 0.5 dB N=5
//! bandpass filter (f0 = 2 GHz, FBW = 0.10, return loss 9 dB, stopband
//! (2.4 GHz, 40 dB)) — the committed `cheb_bpf.toml`-equivalent values. See
//! [`yee_studio`] for the headless logic layer and [`yee_studio::app`] for the
//! UI.
//!
//! Exactly one entry compiles per feature × target:
//! - `all(feature = "desktop", not(target_arch = "wasm32"))` → native
//!   `eframe::run_native` window;
//! - `all(feature = "web", target_arch = "wasm32")` → a `#[wasm_bindgen(start)]`
//!   `eframe::WebRunner` browser entry (App.1.2a; ADR-0096);
//! - otherwise (e.g. `--no-default-features`) → a no-GUI `println!` stub so the
//!   `[[bin]]` target still links while the WASM-safe [`StudioState`] logic
//!   compiles eframe-free (App.1.0; ADR-0092).
//!
//! The windowed binary is build-only in CI (no display); the logic is gated by
//! the headless `studio_state_recompute_*` tests in the library crate.

#[cfg(all(feature = "desktop", not(target_arch = "wasm32")))]
use eframe::egui;

// `default_spec` + `StudioApp`/`StudioState` are shared by the native and web
// entries, so they are gated on either GUI feature.
#[cfg(any(feature = "desktop", feature = "web"))]
use yee_filter::{Approximation, FilterSpec, Response, SpecMask};
#[cfg(any(feature = "desktop", feature = "web"))]
use yee_studio::{StudioState, app::StudioApp};

/// The default satisfiable Chebyshev 0.5 dB N=5 bandpass spec the app opens to.
///
/// Shared by the native (`run_native`) and web (`WebRunner`) entries so both
/// platforms boot to the same design.
#[cfg(any(feature = "desktop", feature = "web"))]
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

/// Native desktop entry: open the eframe window via `run_native`.
#[cfg(all(feature = "desktop", not(target_arch = "wasm32")))]
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

/// Browser entry (App.1.2a; ADR-0096): mount the same [`StudioApp`] onto the
/// `<canvas id="the_canvas_id">` element via `eframe::WebRunner`.
///
/// eframe 0.34's [`eframe::WebRunner::start`] takes a
/// `web_sys::HtmlCanvasElement` (not a string id), so we fetch the canvas from
/// the DOM first. The `index.html` that provides that canvas + the `trunk`
/// bundle are App.1.2b. If WebGPU bindings ever need it, build with
/// `RUSTFLAGS='--cfg=web_sys_unstable_apis'` — not required for this
/// compile-only gate (wgpu 29 falls back to WebGL2).
#[cfg(all(feature = "web", target_arch = "wasm32"))]
#[wasm_bindgen::prelude::wasm_bindgen(start)]
pub fn start_web() {
    use wasm_bindgen::JsCast as _;

    console_error_panic_hook::set_once();

    let web_options = eframe::WebOptions::default();

    wasm_bindgen_futures::spawn_local(async {
        let document = web_sys::window()
            .expect("no global `window` exists")
            .document()
            .expect("should have a document on window");
        let canvas = document
            .get_element_by_id("the_canvas_id")
            .expect("missing element with id `the_canvas_id`")
            .dyn_into::<web_sys::HtmlCanvasElement>()
            .expect("`the_canvas_id` is not a <canvas>");

        eframe::WebRunner::new()
            .start(
                canvas,
                web_options,
                Box::new(|_cc| {
                    Ok(Box::new(StudioApp::new(StudioState::from_spec(
                        default_spec(),
                    ))))
                }),
            )
            .await
            .expect("yee-studio web start failed");
    });
}

/// Required `main` symbol for the wasm32 `[[bin]]` target. The real browser
/// entry is [`start_web`] (`#[wasm_bindgen(start)]`, invoked by the JS loader);
/// `cargo build --bin` for `wasm32-unknown-unknown` still demands a `fn main`,
/// so this is an intentional no-op.
#[cfg(all(feature = "web", target_arch = "wasm32"))]
fn main() {}

/// Stub entry for builds with neither a usable native nor web GUI path (e.g.
/// `--no-default-features`, or `--features web` on a non-wasm target): prints a
/// notice so the `[[bin]] yee-studio` target still links while the WASM-safe
/// [`yee_studio::StudioState`] flow logic compiles eframe-free (App.1.0;
/// ADR-0092).
#[cfg(not(any(
    all(feature = "desktop", not(target_arch = "wasm32")),
    all(feature = "web", target_arch = "wasm32")
)))]
fn main() {
    println!(
        "yee-studio built without a GUI entry (no `desktop` on native / no `web` on wasm32); use the library's StudioState API."
    );
}
