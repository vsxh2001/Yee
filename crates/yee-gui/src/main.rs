//! `yee-gui` — the Phase 1.gui.0 walking-skeleton studio shell.
//!
//! See [`crate::app`] for the application state and tab layout, and
//! [`crate::plots`] for the math helpers behind the dB and Smith plots.
//!
//! Architectural notes:
//!
//! - The shell is intentionally minimal: a menu bar, a metadata side panel,
//!   and two `egui_plot` tabs hosted by `egui_dock`.
//! - File ingestion is driven by a `--file <path>` CLI flag at startup; the
//!   menu entry for `Open .s1p…` is surfaced for discoverability but defers
//!   to that workflow. A real file picker (`rfd`) arrives in Phase 1.gui.1.
//! - Future phases will add multi-port plots, a wgpu 3D viewport, and a live
//!   solver hookup — none of that lives here yet.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use std::path::PathBuf;

mod app;
mod plots;

use crate::app::YeeApp;

fn main() -> anyhow::Result<()> {
    // Initialise tracing so library logs (yee-io etc.) reach stderr.
    let _ = tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .try_init();

    let initial = parse_cli_args();

    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("Yee Studio (Phase 1.gui.0)")
            .with_inner_size([1280.0, 800.0]),
        ..Default::default()
    };

    eframe::run_native(
        "yee-gui",
        native_options,
        Box::new(move |_cc| Ok(Box::new(YeeApp::new(initial.clone())))),
    )
    .map_err(|e| anyhow::anyhow!("eframe::run_native failed: {e}"))
}

/// Minimal hand-rolled CLI parser so we can avoid pulling `clap` into the GUI
/// crate just for one optional flag.
///
/// Supported forms:
/// - `yee-gui --file path/to/foo.s1p`
/// - `yee-gui --file=path/to/foo.s1p`
fn parse_cli_args() -> Option<PathBuf> {
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        if let Some(rest) = arg.strip_prefix("--file=") {
            return Some(PathBuf::from(rest));
        }
        if arg == "--file" {
            return args.next().map(PathBuf::from);
        }
        if arg == "--help" || arg == "-h" {
            eprintln!(
                "yee-gui — Phase 1.gui.0 studio shell\n\n\
                 Usage: yee-gui [--file <path-to.s1p>]\n\n\
                 With no --file flag, the GUI opens to a placeholder; use the\n\
                 flag to preload a Touchstone .s1p file at startup."
            );
            std::process::exit(0);
        }
    }
    None
}
