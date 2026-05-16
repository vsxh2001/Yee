# Contributing to Yee

Thank you for your interest in Yee. This document is a stub during **Phase 0 — Foundation** and will be expanded as the project matures.

## Status

The project is in early bootstrapping. Public APIs, crate boundaries, and validation conventions are still in flux. Expect breaking changes between commits.

Until the first tagged release, the most useful contributions are:

- **Validation cases.** Published-benchmark reproductions for planar MoM (microstrip patches, hairpin filters, Wilkinson dividers, branch-line hybrids). See [ROADMAP.md](ROADMAP.md) for the list.
- **Toolchain issues.** Reports on builds against CUDA 12.4 / 13.0, Rust 1.85+, on Linux / Windows / WSL2.
- **Dependency drift.** Updates to `cudarc`, `pyo3`, `egui`, `opencascade-rs`, `wgpu`. See [TECH_STACK.md](TECH_STACK.md) for pinned versions and watch-outs.
- **Design feedback.** Open a Discussion before opening a PR for anything beyond a typo or a single-file fix.

## Prerequisites

- **Rust** 1.85+ via [rustup](https://rustup.rs/) (the repo pins via `rust-toolchain.toml`).
- **CUDA Toolkit** 12.4+ (13.0 also supported) installed at the system level, with `nvcc` on `PATH`.
- **Python** 3.10+ for the optional Python bindings and validation scripts.
- **Gmsh** 4.15+ available either as the official SDK download or system package (used at build time for FFI bindings).
- A working NVIDIA driver (compute capability ≥ 7.0 recommended).

## Building

```bash
git clone https://github.com/yee-em/yee
cd yee
cargo build --release --features cuda
cargo test --workspace
```

Python wheels (Phase 1):

```bash
pip install maturin
maturin develop --release --features cuda
```

## Coding conventions

- Format with `cargo fmt` and lint with `cargo clippy --workspace --all-targets -- -D warnings`. CI enforces both.
- Keep public APIs documented; every public item gets at least a one-line `///` doc.
- New solver features must ship with a validation case in `validation/` and a citation to the canonical reference.
- CUDA kernels live in `crates/yee-cuda/kernels/` as `.cu` source and are JIT-compiled via NVRTC. No `.ptx` checked in.
- Errors: `thiserror` in library crates, `anyhow` in binaries.

## Filing issues

- **Bugs.** Include the exact command, Rust version (`rustc -V`), CUDA version (`nvcc --version`), GPU model (`nvidia-smi`), and a minimal reproduction.
- **Feature requests.** Link to the canonical reference (paper, textbook, commercial-tool feature) where possible. Tag with the phase you think it belongs to per [ROADMAP.md](ROADMAP.md).
- **Validation discrepancies.** Attach the input geometry, the produced Touchstone, and the reference result. These are the highest-signal reports.

## Pull requests

1. Open a Discussion or issue first for anything beyond a trivial fix.
2. Branch from `main`. Keep PRs focused — one solver concern, one PR.
3. CI must pass on the matrix (Linux + Windows, CUDA 12.4 + 13.0, Rust 1.85 + stable + beta).
4. By contributing, you agree your work is licensed under **GPL v3.0 or later**, matching the project license.

## License

Yee is licensed under [GPL v3.0 or later](LICENSE). All contributions are accepted under the same license.

## Community

- GitHub Discussions — design conversations and Q&A
- Discord (link in [README.md](README.md)) — informal chat once it is set up
- Monthly office hours — to be scheduled once usage justifies it

We will expand this document with a code-of-conduct, a triage rota, and a release process as the project grows.
