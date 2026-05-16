//! Integration tests for the safe `Session` skeleton.
//!
//! Without the `gmsh` feature every method returns `Error::NotEnabled`.

use yee_mesh::{Error, Session};

#[test]
fn session_new_without_feature_returns_not_enabled() {
    let result = Session::new();
    assert!(
        matches!(result, Err(Error::NotEnabled)),
        "default build should report NotEnabled, got {result:?}"
    );
}
