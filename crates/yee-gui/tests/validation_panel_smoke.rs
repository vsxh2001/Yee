//! Smoke test that `ValidationPanel` constructs and that the type
//! is re-exported from the `yee-gui` library crate.
//!
//! Doesn't render — egui needs a context — so this test only
//! verifies the type wires up. The interactive UI surface is
//! exercised by the binary, which is built but not run in CI.

#[test]
fn validation_panel_defaults() {
    let panel = yee_gui::validation::ValidationPanel::default();
    // Drop without panicking — the type is constructible and
    // reachable from the library crate.
    drop(panel);
}
