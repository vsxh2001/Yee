//! Integration tests for the safe `Session` skeleton.
//!
//! Without the `gmsh` feature every method returns `Error::NotEnabled`.
//! With the feature on, method bodies are `todo!()` in Phase 0, so this
//! test only runs in the default (no-feature) build.
//!
//! Coverage note: `Session::new` is the only public constructor and it
//! returns `Err(Error::NotEnabled)` in the no-feature build, so we cannot
//! obtain a `Session` value here to exercise the `import_step`, `mesh`,
//! and `tris` no-feature paths from a black-box integration test. Those
//! bodies are trivial — the same `Err(Error::NotEnabled)` block as `new`
//! under `#[cfg(not(feature = "gmsh"))]` — and are covered by inspection.
//! Phase 1 will add real end-to-end coverage once `Session::new` can
//! succeed.

#![cfg(not(feature = "gmsh"))]

use yee_mesh::{Error, Session};

#[test]
fn session_new_without_feature_returns_not_enabled() {
    let result = Session::new();
    assert!(
        matches!(result, Err(Error::NotEnabled)),
        "default build should report NotEnabled, got {result:?}"
    );
}
