# ADR-0004: Pin egui / eframe / egui_plot / egui_dock to the 0.32 series

**Status:** Accepted (with documented forward path to 0.34)
**Date:** 2026-05-16
**Deciders:** Yee maintainers

## Context

`yee-gui` is the desktop application shell for Yee. It uses the
immediate-mode GUI stack from <https://github.com/emilk/egui>:

- `egui` — the core immediate-mode UI library.
- `eframe` — the windowing and event-loop wrapper around `egui` plus a
  WGPU or glow renderer.
- `egui_plot` — 2D line/scatter plotting, used for the S-parameter
  view.
- `egui_dock` — dockable panel manager, used for the multi-panel layout
  (3D viewport, S-params, log, tree).
- `egui_wgpu` (pulled transitively by `eframe`) — the WGPU back-end.
  Yee uses this rather than the glow back-end because the same WGPU
  device is shared with the custom 3D viewport renderer.

`TECH_STACK.md` (the workspace's living dependency document, written
before Phase 1.gui.0 implementation began) targeted **egui 0.34**, on
the basis that 0.34 was the current series at the time of writing.

The first walking-skeleton implementation pass for Phase 1.gui.0 ran
into a hard build-time block:

```
error: rustc 1.92.0 is required for egui 0.34
       (uses `Vec::extract_if` stabilised in 1.92)
note: workspace MSRV is 1.88 per rust-toolchain.toml
```

Tracing the version graph:

- `egui 0.34.x` requires `rustc >= 1.92` for `Vec::extract_if`,
  `f64::midpoint`, and the `std::sync::OnceLock` use in `Context`.
- `egui_plot 0.34.x` requires `egui 0.34`.
- `egui_dock 0.18.x` requires `egui 0.34`.
- `eframe 0.34.x` pulls `egui_wgpu 0.34` which pulls `wgpu 27`.

Two paths forward:

1. **Bump rust-toolchain.toml to 1.92.** This was rejected on two
   grounds. First, 1.88 was deliberately chosen (ADR-0002) to keep the
   project usable on slower-moving Linux distros and on contributors
   who lag the Rust release cadence; an immediate jump to 1.92 would
   moot ADR-0002 within a single phase. Second, 1.92 was released
   2026-04-23, only weeks before the Phase 1.gui.0 ship target, and
   the maintainers did not want to be the first project on the block
   forcing 1.92 on every contributor and CI runner.
2. **Downgrade egui to a 0.32 series version that supports 1.88.**
   Tested with a scratch branch; all of `egui 0.32`, `egui_plot 0.33`,
   `egui_dock 0.17`, `eframe 0.32` build cleanly on `rustc 1.88.0`.
   The API differences are syntactic and well-documented in the
   upstream changelogs.

Path 2 is the lower-risk choice for shipping a working Phase 1.gui.0.

## Decision

Pin the egui ecosystem to the 0.32 compatibility series in workspace
`Cargo.toml`, with the consistent transitive WGPU version it expects:

```toml
[workspace.dependencies]
egui       = "0.32"
eframe     = { version = "0.32", default-features = false, features = ["wgpu", "default_fonts"] }
egui_plot  = "0.33"   # 0.33 is the egui_plot release that pairs with egui 0.32
egui_dock  = "0.17"   # 0.17 is the egui_dock release that pairs with egui 0.32
wgpu       = "25"     # egui_wgpu 0.32 expects wgpu 25; do not jump to wgpu 26
```

The pin is **range-pinned to the 0.32 minor series** (`"0.32"` in
Cargo, which by SemVer rules means `>=0.32.0, <0.33.0`). Patch updates
within the series are permitted via `cargo update`.

The path forward is documented as a follow-up:

- **Phase 1.gui.3** is the planned bump. It will simultaneously raise
  `rust-toolchain.toml` to whatever Rust stable is current at that
  phase (currently 1.92, likely 1.94 by the time the phase runs),
  bump `egui` back to 0.34 (or its successor), bump `wgpu` to 26 / 27,
  and re-test the dockable panel layout for layout drift.

## Consequences

**What becomes easier:**

- The Phase 1.gui.0 walking skeleton ships against Rust 1.88, keeping
  ADR-0002 intact. No coordinated toolchain bump required mid-phase.
- The pinned WGPU 25 in the GUI matches the WGPU 25 used by the
  bare-metal 3D viewport renderer in `yee-gui::viewport`, removing a
  whole category of "two different `wgpu::Device` versions linked in
  the same binary" build errors.

**What becomes harder:**

- The egui 0.32 API has a small but real set of differences from 0.34
  that show up in the codebase and that contributors who learn egui
  from the upstream README will trip over. The most visible ones:
  - **Menus.** 0.34: `egui::menu::bar(ui, |ui| { ... })`. 0.32:
    `MenuBar::new().ui(ui, |ui| { ... })`.
  - **Plot lines.** 0.34: `Line::new(points)` with no name argument;
    name set via `.name("S11")`. 0.33 (paired with egui 0.32):
    `Line::new("S11", points)` — name is the first positional
    argument.
  - **`Context::request_repaint_after`** signature drift: in 0.32 it
    takes `Duration` directly; in 0.34 it takes a struct.
- Plot styling (axis labels, legend placement) uses an older builder
  pattern; some properties newly added in 0.33/0.34 are unavailable.
- New egui features that land between 0.32 and the eventual 1.gui.3
  bump (improved touch input, native menu integration on macOS, image
  loaders) are unavailable until the bump.

**What's now closed off:**

- Adding any new GUI dependency that itself requires egui 0.34. Such
  a dependency must wait for Phase 1.gui.3 or be re-implemented against
  egui 0.32.
- Adopting WGPU 26 features (subgroup operations, mesh shaders) in
  the GUI viewport renderer until 1.gui.3 takes the coordinated jump.

## References

- `Cargo.toml` workspace dependencies block — egui pin.
- `TECH_STACK.md` — original 0.34 target.
- ADR-0002 — Rust MSRV 1.88.
- egui changelog 0.32 → 0.33 → 0.34:
  <https://github.com/emilk/egui/blob/master/CHANGELOG.md>
- egui_plot 0.33 release notes — `Line::new(name, points)` signature.
- egui_dock 0.17 ↔ egui 0.32 compatibility note in upstream README.
- `yee-gui/src/app.rs` — concentrates the egui API surface that would
  change under a 0.34 bump.
