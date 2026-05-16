//! `yee-gui` library entry: exposes the public modules used by the
//! `yee-gui` binary and by smoke / integration tests under
//! `tests/`.
//!
//! The shell hosts a menu bar, a metadata + viewport-controls side
//! panel, two `egui_plot` tabs, a wgpu-backed 3D viewport tab, and a
//! validation aggregator tab — all hosted by `egui_dock`.
//!
//! See [`app`] for the application state and tab layout,
//! [`plots`] for the math helpers behind the dB and Smith plots,
//! [`viewport`] for the wgpu-backed 3D mesh viewport, and
//! [`validation`] for the validation aggregator panel.
//!
//! The binary entry point lives in `src/main.rs` and re-uses these
//! modules through the library crate so integration tests can
//! construct types like [`validation::ValidationPanel`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod app;
pub mod plots;
pub mod validation;
pub mod viewport;
