# ADR-0002: Pin Rust MSRV to 1.88

**Status:** Accepted
**Date:** 2026-05-16
**Deciders:** Yee maintainers

## Context

Yee declares a Minimum Supported Rust Version (MSRV) so contributors,
distro packagers, and CI all build against the same compiler. The MSRV
gets locked through two complementary mechanisms:

- `rust-toolchain.toml` at the workspace root, which pins the toolchain
  for `rustup` users to a specific version of the official `stable`
  channel.
- `rust-version = "..."` in the workspace `Cargo.toml`, which causes
  `cargo` to refuse to build with an older compiler regardless of which
  toolchain happens to be on `PATH`.

The initial Phase 0 walking-skeleton work pinned MSRV to **1.85** per
`TECH_STACK.md`. That choice was driven by the lowest common denominator
of the two MSRV-sensitive direct dependencies at the time:

- `maturin 1.10`, which was the wheel-build front-end and itself
  required `rustc >= 1.84` for some of its build-script logic.
- `pyo3 0.28`, with an MSRV of 1.85 declared in its
  `Cargo.toml :: package.rust-version`.

During Phase 1 pre-flight (the workspace pass that resolves the
dependency graph against the current registry before any spec gets
locked) two transitive constraints fired:

- `nalgebra 0.34.2` raised its MSRV to **1.87** to make use of the
  `#[diagnostic::on_unimplemented]` stable attribute and to drop legacy
  const-generic workarounds. `nalgebra` is a non-negotiable dependency
  of `yee-core` (`Vector3<f64>`, `Complex<f64>`, dense matrix shapes
  shared between the MoM assembly and the cuSOLVER bridge).
- `libloading 0.9` raised its MSRV to **1.88** to use stable
  `core::ffi::c_str` and the `unsafe` extern-block lint. `libloading`
  is pulled transitively by `cudarc` for `dlopen`-style CUDA driver
  loading and by `gmsh-sys` for the optional run-time-loaded Gmsh
  shared library.

Three responses were on the table:

1. Refuse to upgrade the dependencies, stay on 1.85, hand-roll
   work-arounds for everything that wanted 1.88. The maintenance cost
   of forked dependencies was judged unacceptable for a project whose
   value proposition includes "Rust ecosystem leverage".
2. Cap MSRV at the minimum the transitive graph requires (1.88) and
   move on. Reproducible builds, well-documented, ecosystem-aligned.
3. Bump aggressively to the current stable (1.92 at the time of this
   ADR) so we can later use `egui 0.34` and other newer Rust-2024
   features. Considered and rejected — see ADR-0004 — because
   `rust-toolchain.toml = 1.92` would cut off contributors on Debian
   stable and on slower-moving Linux distros that lag the Rust release
   cadence by 9–12 months.

Option 2 is the chosen compromise: track the transitive MSRV floor, but
do not chase the bleeding edge.

## Decision

Pin the Yee workspace MSRV to **Rust 1.88.0**. Concretely:

- `rust-toolchain.toml`:

  ```toml
  [toolchain]
  channel = "1.88.0"
  components = ["rustfmt", "clippy"]
  profile = "minimal"
  ```

- Workspace `Cargo.toml`:

  ```toml
  [workspace.package]
  rust-version = "1.88"
  ```

- CI matrix (`.github/workflows/ci.yml`) builds against `1.88.0`
  exactly, plus `stable` (currently 1.92) as an early-warning canary.
  Only the 1.88.0 job is required to merge.

- Any new dependency whose declared MSRV exceeds 1.88 must either pin
  an older compatible version in `Cargo.toml` or wait for a coordinated
  MSRV bump.

## Consequences

**What becomes easier:**

- Reproducible builds across all maintainer hosts and CI runners.
  `rust-toolchain.toml` is checked in, so `cargo build` after a fresh
  `git clone` triggers `rustup` to fetch the exact compiler.
- A predictable MSRV lets us write idiomatic 2024-edition Rust without
  fear, including `let-else`, `async fn` in traits, the 2024-edition
  prelude updates, and stable `#[diagnostic::*]` attributes.
- Distro packagers can target 1.88 without needing nightly or the
  current stable. Debian trixie (Rust 1.81 in `unstable` as of May
  2026) still cannot package Yee, but Fedora 41 (Rust 1.85+) is close.

**What becomes harder:**

- Some downstream Rust ecosystem dependencies are now unreachable on
  this toolchain:
  - `egui 0.34` (requires `rustc 1.92` for new `f64::midpoint` and
    `Vec::extract_if` stabilisations) — see ADR-0004.
  - `wgpu 27` (uses `naga 26` which itself requires 1.91).
  - `rust-analyzer >= 0.4.2200` — affects IDE tooling only, not the
    build.
- A future MSRV bump (planned for Phase 1.gui.3, see ADR-0004) is a
  coordinated change: bump `rust-toolchain.toml`, bump
  `workspace.package.rust-version`, bump the CI matrix, run the full
  test suite, audit any docs that mention "1.88" (this ADR, the
  README, `TECH_STACK.md`).
- Contributors using Homebrew on macOS who routinely run `brew upgrade
  rust` will find their system Rust ahead of the project's toolchain;
  the workspace falls back to the `rust-toolchain.toml` pin
  automatically, but the mental model is non-obvious.

**What's now closed off:**

- Using 1.85-specific features (e.g. older `#[doc(cfg)]` syntax) is
  not closed off — they continue to compile under 1.88. The MSRV bump
  is forwards-only.
- Going below 1.88 would re-impose the pre-`libloading 0.9` and
  pre-`nalgebra 0.34.2` versions, which are not patch-supported by
  their upstream maintainers.

## References

- `rust-toolchain.toml` (workspace root).
- Workspace `Cargo.toml`, `[workspace.package] rust-version = "1.88"`.
- `TECH_STACK.md` — original 1.85 selection and rationale.
- ADR-0004 — egui pin to 0.32 series, downstream consequence of 1.88
  MSRV.
- Rust Edition Guide, Rust 1.85 release notes (Rust 2024 edition):
  <https://doc.rust-lang.org/edition-guide/rust-2024/>
- `libloading 0.9` changelog noting MSRV bump:
  <https://github.com/nagisa/rust_libloading/releases/tag/0.9.0>
- `nalgebra 0.34.2` changelog noting MSRV 1.87:
  <https://github.com/dimforge/nalgebra/releases/tag/v0.34.2>
