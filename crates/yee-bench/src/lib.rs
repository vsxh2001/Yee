//! # yee-bench
//!
//! Criterion benchmarks for Yee's hot paths.
//!
//! This crate intentionally contains no library code — it exists only to host
//! the `benches/*.rs` binaries that exercise [`yee_mom`], [`yee_fdtd`], and
//! the iterative solvers in [`yee_mom::iterative`]. Run them with:
//!
//! ```text
//! cargo bench -p yee-bench
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]
