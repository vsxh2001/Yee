# Phase 1.gui.3 — Rust 1.92 + egui 0.34 + wgpu 26 toolchain bump

**Status:** Draft  
**Owner:** TBD  
**Phase:** 1.gui.3  
**Depends on:** Phase 1.gui.0, 1.gui.1, 1.gui.2 (shipped)  
**Blocks:** future egui-ecosystem upgrades; nothing critical

## Assumption being challenged

CLAUDE.md §10 records that egui is pinned to 0.32 (not the TECH_STACK target 0.34) because 0.34 requires rustc 1.92, and rust-toolchain.toml pins 1.88. We then accepted the downstream pins: egui_plot 0.33, egui_dock 0.17, egui-wgpu 0.32, wgpu 25 (downgraded from 26).

The assumption: **bumping the toolchain is a coordinated multi-crate change that warrants its own phase.**

By 2026-05 (current date) rustc 1.92 has been stable for multiple releases. The cost of staying on the older ecosystem (missed bug fixes, missed wgpu features, growing distance from upstream) starts to outweigh the bump cost. Time to actually do it.

## Scope

In:
- `rust-toolchain.toml` 1.88 → 1.92
- `Cargo.toml` workspace `rust-version` 1.88 → 1.92
- `[workspace.dependencies]`: egui 0.32 → 0.34, eframe 0.32 → 0.34, egui_plot 0.33 → 0.34, egui_dock 0.17 → 0.18, egui-wgpu 0.32 → 0.34, wgpu 25 → 26
- Any API breaks in `crates/yee-gui/src/{app,plots,viewport}.rs`
- CI matrix `dtolnay/rust-toolchain` toolchain string

Out:
- New GUI features
- Re-theming
- Plot regressions

## Definition of done

1. `rust-toolchain.toml` reads `1.92`. `Cargo.toml` workspace `rust-version` reads `1.92`.
2. All five GUI deps bumped per scope. `Cargo.lock` updated by `cargo update -p <crate>` only on the affected crates (don't sweep transitive).
3. `cargo check --workspace --no-default-features` green.
4. `cargo clippy --workspace --all-targets --no-default-features -- -D warnings` green.
5. `cargo test --workspace --no-default-features` green.
6. `cargo fmt --check --all` green.
7. `.github/workflows/ci.yml` and `gpu-nightly.yml` toolchain strings updated to `"1.92"`.
8. CLAUDE.md §10 entries describing the 1.88 pin and the wgpu 25 downgrade — **deleted**, not edited. Remove the rationale: the pin is gone.

## Verification

```bash
cargo check --workspace --no-default-features
cargo clippy --workspace --all-targets --no-default-features -- -D warnings
cargo test --workspace --no-default-features
cargo fmt --check --all
```

GUI smoke: `cargo run -p yee-gui --release --no-default-features -- --file examples/touchstone/dipole.s1p` opens the window, shows the S11 dB plot + Smith chart, no panic.

## Known API breaks expected (to confirm during the bump)

- egui 0.34 renamed `egui::Context::input` access patterns; check viewport code.
- wgpu 26 stabilized `Device::create_render_pipeline` validation paths; egui-wgpu adapts internally but custom wgpu callbacks may need adjustment.
- egui_dock 0.18 changed `DockState::iter_tabs` signature.

If any of these surface, fix in lane (`crates/yee-gui/**`) and note in commit body.

## Escape hatch

If `cargo check` after the version bump fails with a non-trivial API break (>30 min to chase), revert the bump and surface the specific upstream change as a follow-up ticket. The 1.88/0.32 pin is not catastrophic; we can afford one more cycle.

## Reference

- Rust release notes 1.89 → 1.92
- egui changelog https://github.com/emilk/egui/blob/master/CHANGELOG.md (0.32 → 0.34)
- wgpu changelog https://github.com/gfx-rs/wgpu/blob/trunk/CHANGELOG.md (25 → 26)
