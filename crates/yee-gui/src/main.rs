//! `yee-gui` — the Phase 1.gui.0 walking-skeleton studio shell.
//!
//! See [`crate::app`] for the application state. This entry point initialises
//! `tracing` and hands off to `eframe::run_native`.
//!
//! Future commits in this phase will add the S11 / Smith plotting tabs and
//! the `--file` CLI flag.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod app;

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
        Box::new(|_cc| Ok(Box::new(YeeApp::new()))),
    )
    .map_err(|e| anyhow::anyhow!("eframe::run_native failed: {e}"))
}
