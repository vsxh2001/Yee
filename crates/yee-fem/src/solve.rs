//! Sparse generalised eigensolver — placeholder.
//!
//! This module will own:
//!
//! * The `SparseEigen` trait that abstracts the `K e = k² M e` solve
//!   (load-bearing decision per Phase 4 spec §8 — keeps the library
//!   choice swappable in one PR).
//! * `LobpcgEigen`, the pure-Rust shift-invert LOBPCG implementation
//!   that consumes `nalgebra_sparse::CsrMatrix<f64>` and produces an
//!   `EigenpairList`. Shift-invert uses `faer` sparse LU as the inner
//!   `(K − σM)^{-1}` preconditioner.
//!
//! Delivered by step **T5** of the SSSSS plan
//! (`docs/superpowers/plans/2026-05-18-phase-4-fem-eigenmode.md`). The
//! escape hatch — if `lobpcg`'s public API has regressed between spec
//! freeze and implementation — is a hand-rolled deflated inverse-power
//! iteration on `faer` sparse LU; the trait keeps the swap downstream-
//! invisible.
//!
//! The Phase 4 design spec at
//! `docs/superpowers/specs/2026-05-18-phase-4-fem-eigenmode-design.md`
//! locks the public surface (`SparseEigen`, `LobpcgEigen`,
//! `EigenpairList`). Nothing is exposed here yet — this scaffold module
//! exists so the crate compiles end-to-end while T5 is in flight.
