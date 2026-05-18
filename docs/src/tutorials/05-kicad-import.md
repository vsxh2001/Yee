# Tutorial 5 — KiCad PCB import to MoM

This tutorial walks through Yee's Phase 1.mesh.1 `.kicad_pcb` importer
from a small Rust driver. You will hand-author a KiCad PCB containing
one signal trace and a ground plane, parse it with
`yee_mesh::KiCadBoard`, inspect the extracted geometry, and stage the
data for a follow-on MoM solve.

## Goal

Read a hand-authored `.kicad_pcb` (an open-ended microstrip stub on
FR4), print the parsed thickness / layers / segments / zones, and
sanity-check that the geometry matches what KiCad wrote out. End state:
a printout you can paste into a bug report and a clear understanding of
which parts of the KiCad file format the Phase 1.mesh.1 parser supports.

## Prerequisites

- A built Yee workspace (`cargo build --release` from the repo root).
- A `.kicad_pcb` source file. KiCad 7+ is the reference version, but
  the parser only cares about the s-expression subset covered below —
  older KiCad files that share that subset also parse. You do **not**
  need the KiCad GUI; the file is plain text and a minimal sample is
  inline below.
- The KiCad SDK is **not** required. `yee-mesh` parses the
  `.kicad_pcb` s-expression directly with a small in-tree tokenizer.
- The `gmsh` feature on `yee-mesh` is **not** required; KiCad parsing
  is independent of the Gmsh FFI.

## What `yee-mesh` actually parses from `.kicad_pcb`

Per `crates/yee-mesh/src/kicad.rs`, the Phase 1.mesh.1 walking-skeleton
scope is deliberately narrow. The parser extracts:

- The top-level `(kicad_pcb ...)` wrapper (anything else is rejected
  with `KiCadError::Unsupported`).
- `(general (thickness ...))` → total board thickness in millimetres.
- `(layers (N "Name" type) ...)` → ordered list of `LayerInfo
  { ordinal, name, kind }`.
- `(segment (start x y) (end x y) (width w) (layer "F.Cu") ...)` →
  copper trace `Segment` records (preserving the layer name verbatim).
- `(zone ... (layer "B.Cu") (polygon (pts (xy x y) ...)))` → copper
  zone fills, captured as a `Zone { layer, polygon }`.

Everything else is **silently skipped**: footprints, vias, drills,
graphical primitives (`gr_line` / `gr_arc` / `gr_circle`), silkscreen,
solder-mask, fabrication outlines, 3-D models, and net metadata. Full
footprint / via / arc parsing is the Phase 1.mesh.2 follow-up.

KiCad's native unit in the s-expression is **millimetres**. The parser
preserves that — every numeric field on `Segment` / `Zone` is in mm.
Converting to SI metres (`* 1e-3`) is the caller's job before feeding
`yee-mom`.

## Workflow

Three stages. Today all three are driven from a small Rust binary; CLI
wiring is on the Phase 1.mesh.2 follow-up list (see
[Limitations](#limitations)).

### Stage 1 — parse the `.kicad_pcb`

`yee_mesh::KiCadBoard::read` takes a path and returns the parsed
struct. From a Rust driver:

```rust
use std::path::Path;
use yee_mesh::KiCadBoard;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let board = KiCadBoard::read(Path::new("examples/stub.kicad_pcb"))?;
    println!("thickness: {:.3} mm", board.thickness_mm);
    println!("layers   : {}", board.layers.len());
    for l in &board.layers {
        println!("  {:2} {:<8} {}", l.ordinal, l.name, l.kind);
    }
    println!("segments : {}", board.segments.len());
    for s in &board.segments {
        println!(
            "  ({:.3},{:.3}) -> ({:.3},{:.3})  w={:.3}  layer={}",
            s.start.0, s.start.1, s.end.0, s.end.1, s.width_mm, s.layer
        );
    }
    println!("zones    : {}", board.zones.len());
    for z in &board.zones {
        println!("  layer={}  polygon vertices={}", z.layer, z.polygon.len());
    }
    Ok(())
}
```

Save the snippet under `examples/kicad-inspect/src/main.rs` in a
throwaway crate that depends on `yee-mesh` (path dependency on the
workspace root), then run it:

```bash
cargo run --release --bin kicad-inspect -- examples/stub.kicad_pcb
```

### Stage 2 — convert to a MoM mesh

`KiCadBoard` is **not** a `TriMesh` — it is the *raw geometry* the next
stage triangulates. Two things have to happen before a MoM solve:

1. **Filter by layer.** Pick one signal layer (typically `"F.Cu"`) for
   the trace. The ground plane lives on `"B.Cu"`. Anything in between
   (`"In1.Cu"`, `"In2.Cu"`, ...) is a stack-up layer that the current
   single-layer planar MoM solver cannot consume.
2. **Triangulate.** Each `Segment` is a rectangular strip
   (`width_mm` × segment length); a `Zone` is a polygon outline. The
   strip is decomposed into RWG-compatible triangles in the caller —
   today by hand, in Phase 1.mesh.2 via Gmsh's OCC kernel through the
   existing `yee_mesh::Session` API.

```rust
let mm_to_m = 1.0e-3;
let f_cu: Vec<_> = board
    .segments
    .iter()
    .filter(|s| s.layer == "F.Cu")
    .map(|s| {
        // Each segment becomes a rectangular strip in SI metres.
        // Hand the four corners to your downstream meshing code.
        (
            (s.start.0 * mm_to_m, s.start.1 * mm_to_m),
            (s.end.0 * mm_to_m, s.end.1 * mm_to_m),
            s.width_mm * mm_to_m,
        )
    })
    .collect();
println!("F.Cu segments (in metres): {}", f_cu.len());
```

Extend `kicad-inspect/src/main.rs` with the filter above and rerun:

```bash
cargo run --release --bin kicad-inspect -- examples/stub.kicad_pcb
```

### Stage 3 — solve and plot

Once you have a `TriMesh` for the trace and a `TriMesh` for the ground
plane, the rest of the pipeline is identical to
[Tutorial 2](02-dipole-from-python.md) (Python) or
[Tutorial 1](01-microstrip-line.md) (Rust): build a `PlanarMoM`, attach
a delta-gap or wave port at the trace end, sweep frequency, and write a
Touchstone file with `yee_io::touchstone::write`. The resulting `.s1p`
plots with:

```bash
yee plot stub.s1p --format smith --output stub-smith.png
yee plot stub.s1p --format db    --output stub-db.png
```

`yee plot` is wired up against `yee-plotters` today and works against
any Touchstone file regardless of which solver produced it.

## A concrete worked example

Below is a hand-authored 10 mm open-ended microstrip stub: a single
`F.Cu` trace 1.6 mm wide running from `(10, 50)` to `(20, 50)` over a
ground plane that fills the whole 100 × 100 mm board on `B.Cu`. The
board is 1.6 mm thick FR4. Save the listing as `stub.kicad_pcb`:

```text
(kicad_pcb (version 20221018) (generator pcbnew)
  (general (thickness 1.6))
  (layers
    (0  "F.Cu" signal)
    (31 "B.Cu" signal)
  )
  (segment (start 10 50) (end 20 50) (width 1.6) (layer "F.Cu") (net 1))
  (zone (net 0) (net_name "GND") (layer "B.Cu")
    (polygon (pts (xy 0 0) (xy 100 0) (xy 100 100) (xy 0 100)))
  )
)
```

Run the Stage 1 driver against it. Expected output:

```text
thickness: 1.600 mm
layers   : 2
   0 F.Cu     signal
  31 B.Cu     signal
segments : 1
  (10.000,50.000) -> (20.000,50.000)  w=1.600  layer=F.Cu
zones    : 1
  layer=B.Cu  polygon vertices=4
```

That is the entire imported geometry today. The trace is 10 mm long,
1.6 mm wide; the ground plane is a 100 × 100 mm rectangle on the back
side of a 1.6 mm FR4 stack-up.

### Sanity check against the open-stub formula

For a lossless open-ended stub of length `l` at angular frequency `ω`
with phase constant `β = ω √(ε_eff) / c`, the textbook input impedance
is

```text
Z_in ≈ −j · Z₀ · cot(β · l)
```

For a 1.6 mm-wide trace on 1.6 mm FR4 (`ε_r ≈ 4.4`), the standard
microstrip closed-form (Pozar §3.8) gives `Z₀ ≈ 50 Ω` and
`ε_eff ≈ 3.2`. At `f = 2.4 GHz`:

```text
β  = 2π · 2.4e9 · √3.2 / 2.998e8 ≈ 90.0 rad/m
l  = 10e-3 m
βl ≈ 0.900 rad   →   cot(βl) ≈ 0.794
|Z_in| ≈ 50 · 0.794 ≈ 40 Ω      (capacitive — open stub below λ/4)
```

A MoM solve fed by this imported geometry should land within a few
ohms of that magnitude once the multilayer Green's function and a
correctly-defined microstrip port are wired in. **It will not match
today** — see [Limitations](#limitations) for what is still loose.

## Inspecting the imported mesh

There is no `yee mesh show` subcommand at base SHA; the inspection path
is the Stage-1 driver above. Alternatives: pipe the parsed listing
into `yee-plotters` after triangulating, or use the desktop GUI
(`cargo run -p yee-gui --release`) to load the resulting Touchstone
file once Stage 3 has run. The GUI does not yet visualize raw imported
geometry (Phase 1.gui.4).

## Common pitfalls

- **Layer name conventions.** KiCad's `"F.Cu"` (front copper) and
  `"B.Cu"` (back copper) are the two-layer-board defaults; internal
  layers are `"In1.Cu"`, `"In2.Cu"`. The parser preserves the name
  verbatim, so filter on `segment.layer == "F.Cu"` (case-sensitive,
  literal dot).
- **Units.** KiCad files ship in millimetres. `yee-mesh` *does not*
  silently convert to SI metres — every `f64` on `Segment` / `Zone` is
  in mm. Multiply by `1.0e-3` before handing the data to `yee-mom`.
- **Multi-net boards.** `(net N)` / `(net_name ...)` is parsed *and
  ignored* in Phase 1.mesh.1. To pick one net, filter the `segments`
  slice in the caller — no `KiCadBoard::filter_by_net` yet.
- **Footprints, vias, arcs.** Copper drawn inside `(footprint ...)`,
  any `(via ...)` form, and arcs are **silently skipped**. Move copper
  to top-level `segment` / `zone` forms or wait for Phase 1.mesh.2.
- **The `yee mesh` CLI subcommand.** At base SHA `905ca8b` it only
  constructs a Gmsh `Session`; it does not yet route `.kicad_pcb`
  inputs to the `KiCadBoard` parser. Drive the parser from a Rust
  binary (Stage-1 snippet) until that wiring lands.

## Limitations

- **Single-layer planar MoM only at base SHA.** The current solver
  treats the trace as a 2-D PEC sheet over a `MicrostripDcim`
  (one-image DCIM placeholder) Green's function. Full 3-D PCB
  stack-up — multiple `In*.Cu` layers, dielectric thicknesses, vias —
  is the FDTD path (Phase 2) and the future Phase 4 FEM solver.
- **Loose tolerance on mom-002-class problems** until Phase 1.1.1.2
  Sommerfeld pole extraction lands. The placeholder `MultilayerGreens`
  exercises the I/O contract end-to-end but does not validate a 50 Ω
  microstrip Z₀ to better than ±20 %. Do not over-trust a single
  number out of this pipeline yet.
- **No automated KiCad → MoM driver.** The three stages above are
  manual today; Phase 1.mesh.2 stitches them together behind a
  single CLI invocation.
- **Phase 1.mesh.1 parser scope is deliberate.** See the inline doc at
  the top of `crates/yee-mesh/src/kicad.rs` for the canonical scope
  list.

## Next steps

- Read [Multilayer Green's Function — DCIM and Surface-Wave
  Poles](../theory/multilayer-greens.md) for the theory behind the
  planar microstrip solver this tutorial feeds.
- Skim `docs/superpowers/specs/2026-05-17-phase-1-1-1-2-sommerfeld-pole-extraction-design.md`
  for the spec that closes the mom-002 / mom-003 accuracy gap.
- Once Phase 1.mesh.2 lands, this tutorial will be rewritten to drive
  `yee mesh import` → `yee run` → `yee plot` from one shell session.
