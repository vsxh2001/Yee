//! `yee-gui` — the Phase 1.gui.0 walking-skeleton studio shell.
//!
//! See [`crate::app`] for the application state and tab layout, and
//! [`crate::plots`] for the math helpers behind the dB and Smith plots.
//!
//! The next commit in this phase wires up the `--file <path>` CLI flag so
//! the GUI can preload a Touchstone file at startup.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod app;
mod plots;

use crate::app::YeeApp;

fn main() -> anyhow::Result<()> {
    let _ = tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .try_init();

    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("Yee Studio (Phase 1.gui.0)")
            .with_inner_size([1280.0, 800.0]),
        ..Default::default()
    };

    eframe::run_native(
        "yee-gui",
        native_options,
        Box::new(|_cc| Ok(Box::new(YeeApp::new(None)))),
    )
    .map_err(|e| anyhow::anyhow!("eframe::run_native failed: {e}"))
}
