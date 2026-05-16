# yee-gui

Phase 1.gui.0 walking-skeleton studio shell for the Yee electromagnetic
simulation workspace.

This is a deliberately minimal `eframe` + `egui_dock` + `egui_plot` desktop
shell that loads a Touchstone `.s1p` file and visualises its `S11` reflection
coefficient in two side-by-side tabs:

- **S11 magnitude (dB)** — `20·log10|S11|` vs frequency (GHz on X).
- **Smith chart** — `S11` trajectory in the complex plane, drawn on a
  data-aspect 1:1 canvas with the unit circle for reference.

A left side panel surfaces the loaded file's metadata: port count,
reference impedance, on-disk format, frequency unit, sample count, and
frequency span. Preserved Touchstone comments are echoed below the
metadata block.

Phase 1.gui.1+ will add multi-port plots, a wgpu 3D viewport, a real
file picker, and the live solver hookup. None of those are wired in yet.

## Build & run

From the workspace root:

```bash
cargo run -p yee-gui --release
```

With no `--file` flag the GUI opens to a friendly "Open a .s1p file to
begin" placeholder.

To pre-load a Touchstone file at startup:

```bash
cargo run -p yee-gui --release -- --file path/to/dipole.s1p
```

Either `--file <path>` or `--file=<path>` is accepted. A sample 1-port
fixture lives at
`crates/yee-io/validation/fixtures/touchstone/1port.s1p`:

```bash
cargo run -p yee-gui --release -- \
    --file crates/yee-io/validation/fixtures/touchstone/1port.s1p
```

## Why a CLI flag instead of a file picker

Phase 1.gui.0 explicitly defers the `rfd` file-picker integration to keep
the dependency surface small and the build fast. The `File → Open .s1p…`
menu entry is rendered for discoverability but currently points the user
at the `--file` flag. A real picker arrives in Phase 1.gui.1.

## Logging

The binary initialises `tracing-subscriber` with `EnvFilter`:

```bash
RUST_LOG=info cargo run -p yee-gui --release -- --file foo.s1p
RUST_LOG=debug cargo run -p yee-gui --release
```
