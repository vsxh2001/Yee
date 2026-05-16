# ADR-0011: Bump Rust MSRV to 1.92 and lift the egui / wgpu stack to 0.34 / 29

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

ADR-0002 pinned the workspace MSRV to **Rust 1.88** because the
transitive dependency graph at the time of Phase 0 required at most
1.88 (`libloading 0.9`, `nalgebra 0.34.2`) and the maintainers
explicitly did not want to be the first project on the block forcing
1.92 on every contributor and CI runner.

ADR-0004 was the direct downstream consequence: the egui ecosystem was
pinned to the **0.32** compatibility series (`egui 0.32`, `egui_plot
0.33`, `egui_dock 0.17`, `eframe 0.32`, `wgpu 25`) because `egui 0.34`
had raised its own MSRV to 1.92 to use `Vec::extract_if`,
`f64::midpoint`, and `OnceLock` in `Context`. ADR-0004 documented a
**forward path** to be taken in Phase 1.gui.3: simultaneously raise the
workspace toolchain and lift the egui stack back to its current series.

By 2026-05-17 the calculus for the Phase 1.gui.3 bump had shifted:

- **Rust 1.92 has been the stable channel for multiple releases.**
  1.92 shipped 2026-04-23; by 2026-05-17 it has been on `stable` for
  several weeks and 1.93 is already the current `stable`. The "we
  would be the first project forcing 1.92" argument from ADR-0002 no
  longer applies — `egui`, `wgpu`, `nalgebra`, and most of the
  ecosystem we depend on have already moved.
- **Pinning 0.32 has a maintenance tax.** The egui 0.32 API has
  small-but-real differences from 0.34 documented in ADR-0004's
  Consequences section (menu builder, `Line::new(name, points)`
  positional name, `request_repaint_after` signature). Contributors
  who learn egui from upstream docs trip over each of these.
  Bug-fix patch releases against the 0.32 series have slowed to a
  trickle as upstream attention moved to 0.34/0.35.
- **The wgpu version path is no longer 25 → 26.** When ADR-0004 was
  written the planned bump target was wgpu 26. Between 0.32 and 0.34
  egui-wgpu skipped a wgpu version: `egui-wgpu 0.34` **hard-requires
  wgpu 29**, not 26. There is no wgpu 26-compatible egui 0.34
  release. The bump is therefore not 25 → 26 but 25 → 29 in a single
  step.
- **`egui_plot` and `egui_dock` minor versions track egui.** The
  highest releases that pin `egui ^0.34` at the time of bump are
  `egui_plot 0.35` and `egui_dock 0.19`. There is no `egui_plot 0.34`
  or `egui_dock 0.18` — the upstream crates do not version-lockstep
  with the core `egui` minor.
- **Phase 1.gui.0/1/2 has shipped and stabilised.** The walking
  skeleton (`yee-gui` shell, S₁₁ dB plot, Smith chart) is in place;
  this is the moment in the phase chain where ADR-0004 said the bump
  should happen.

The cost of the bump is small and bounded:

- **`cargo bench` invocation drift.** The `Criterion` interaction
  with the new `Bencher` signature on wgpu 29's example benchmarks
  needed light fixup; one bench file in `yee-gui` required a
  one-line type-inference annotation.
- **`--all-targets` clippy lints.** Rust 1.92 enables a handful of
  new lints (`needless_pass_by_ref_mut`, sharper
  `redundant_closure_call`) that fired in a few places under
  `-D warnings`. Each fix was mechanical.
- **No physics test failures.** The MoM and FDTD test suites are
  toolchain-insensitive at this granularity; `cargo test
  --workspace` is green after the bump.

## Decision

Bump the workspace toolchain and lift the egui / wgpu stack to its
current compatible series. Specifically:

- **`rust-toolchain.toml`**:

  ```toml
  [toolchain]
  channel = "1.92.0"
  components = ["rustfmt", "clippy"]
  profile = "minimal"
  ```

- **Workspace `Cargo.toml`**:

  ```toml
  [workspace.package]
  rust-version = "1.92"
  ```

- **Workspace `[workspace.dependencies]`**:

  ```toml
  egui       = "0.34"
  eframe     = { version = "0.34", default-features = false, features = ["wgpu", "default_fonts"] }
  egui_plot  = "0.35"   # highest release pinning egui ^0.34
  egui_dock  = "0.19"   # highest release pinning egui ^0.34
  egui-wgpu  = "0.34"
  wgpu       = "29"     # forced by egui-wgpu 0.34's transitive bound
  ```

- **CI workflows** (`.github/workflows/ci.yml`,
  `gpu-nightly.yml`, `publish-wheels.yml`) update their toolchain
  matrix entry from `1.88.0` to `1.92.0`. The `stable` early-warning
  canary job is retained.

- **Documentation sweep.** `TECH_STACK.md`, root `CLAUDE.md` (§3, §7,
  §8, §10), and ADR-0002 / ADR-0004 references to "1.88" are updated
  to "1.92". ADR-0002 and ADR-0004 themselves are left as historical
  records with cross-references pointing here.

The egui-API shim sites that ADR-0004 flagged (menu builder,
`Line::new` positional name, `request_repaint_after` signature) are
migrated to the 0.34 form in the same change.

## Consequences

**What becomes easier:**

- The toolchain pin recorded in `CLAUDE.md §10` and in ADR-0004's
  "documented forward path" is closed. The "do not bump egui
  unilaterally" gotcha is removed from `CLAUDE.md §10`.
- Future egui ecosystem upgrades become routine `cargo update`
  candidates rather than coordinated multi-crate bumps. Adding new
  GUI dependencies that require `egui ^0.34` is unblocked (ADR-0004
  "What's now closed off" item 1 reopens).
- The wgpu 25 → 29 jump pulls in subgroup operations, mesh shaders
  (where the platform supports them), and the modernised
  `wgpu::Surface` configuration API. The `yee-gui::viewport`
  bare-metal renderer can adopt these in follow-up work without
  another coordinated bump.
- Contributors who learn egui from the upstream README now hit the
  same API as the in-tree code.

**What becomes harder:**

- Distro packagers on Rust < 1.92 are now locked out. As of
  2026-05-17 this covers Debian stable and Fedora 41; Fedora 42 ships
  1.92 and Debian trixie's next point release is expected to ship
  1.93. The window is small but real.
- Contributors on Homebrew / `rustup` who lag the `stable` channel
  will be prompted to fetch 1.92 on first `cargo build`; this is
  visible but not blocking.
- The wgpu 25 → 29 jump is a three-minor-version skip. Any code in
  `yee-gui::viewport` that read undocumented wgpu 25 internals
  (none does, by audit) would have broken silently; the bump was
  validated by running the GUI end-to-end against the half-wave-
  dipole and microstrip-line example pipelines.

**What's now closed off:**

- Reverting to 1.88 / egui 0.32 without re-opening this ADR. The
  bump is forwards-only; downstream code is free to use 1.92
  features (`Vec::extract_if`, `f64::midpoint`, edition-2024
  prelude updates) from here on.
- Holding `wgpu` at 26 (the original ADR-0004 plan). The
  `egui-wgpu 0.34` transitive bound makes 26 unreachable; the bump
  is 25 → 29 or nothing.

## References

- `rust-toolchain.toml` — toolchain pin (now 1.92.0).
- Workspace `Cargo.toml` — `rust-version = "1.92"` and the egui /
  wgpu `[workspace.dependencies]` block.
- `TECH_STACK.md` — egui / wgpu target versions, updated in lockstep.
- `.github/workflows/ci.yml`, `.github/workflows/gpu-nightly.yml`,
  `.github/workflows/publish-wheels.yml` — toolchain matrix.
- ADR-0002 — original MSRV 1.88 decision; superseded for the MSRV
  pin specifically, kept for the rationale on tracking the transitive
  floor.
- ADR-0004 — original egui 0.32 pin and the forward path this ADR
  takes.
- egui changelog 0.32 → 0.33 → 0.34 → 0.35:
  <https://github.com/emilk/egui/blob/master/CHANGELOG.md>
- `egui-wgpu 0.34` release notes — wgpu 29 transitive requirement.
- `wgpu 29` release notes:
  <https://github.com/gfx-rs/wgpu/releases/tag/v29.0.0>
