//! Shim binary that exposes the fixtures module so its inner tests run
//! under `cargo test -p yee-mom --test fixtures`.

#[path = "fixtures/mod.rs"]
mod fixtures;
