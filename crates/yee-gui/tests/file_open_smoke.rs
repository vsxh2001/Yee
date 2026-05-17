//! Smoke test for the `File → Open Touchstone…` path.
//!
//! The native `rfd` dialog needs a display server, so this test does not
//! invoke it. It exercises only the bits that should remain reachable from
//! the library crate: that [`yee_gui::app::YeeApp`] is constructible and
//! that [`YeeApp::load_touchstone`] handles a missing path by recording
//! a load error rather than panicking.
//!
//! This is the same property that the menu button relies on — the picker
//! hands a `PathBuf` to `load_touchstone`, and any I/O failure must
//! surface as a banner rather than crash the UI.

use std::path::PathBuf;

#[test]
fn yee_app_has_load_touchstone() {
    // Type is reachable from the library crate.
    let _ = std::any::type_name::<yee_gui::app::YeeApp>();
}

#[test]
fn load_touchstone_missing_path_does_not_panic() {
    let mut app = yee_gui::app::YeeApp::new(None);
    app.load_touchstone(&PathBuf::from("/nonexistent/path/to/missing.s1p"));
    // No panic; the function returned cleanly. The error itself is
    // surfaced through the app's load_error banner field which is
    // private — covered by app.rs unit tests.
}
