# ADR-0229: FS.3.2c — imported-board-vs-native-twin measurement gate

**Date:** 2026-07-24 · **Status:** accepted · **Track:** FS.3 (`FULL-SUITE-ROADMAP.md`)
**Spec:** `docs/superpowers/specs/2026-07-24-fs3-2c-import-twin-design.md`
**Plan:** `docs/superpowers/plans/2026-07-24-fs3-2c-import-twin.md`
**Predecessors:** FS.3.0 Gerber import subset (ADR-0209, `gerber-rt-001` byte-identical
round-trip); FS.3.1 `gerber_to_outline` + studio import (`gerber-rt-002` corner-exact,
`studio-import-e2e-001`); FS.3.2a geometry generality (ADR-0217, `engine-miter-001`);
FS.3.2b arcs/flashes (ADR-0220, `gerber-rt-003`).

## Context

Round-trip byte-identity proves the file layer: `export∘import∘export` reproduces
the writer's own bytes. It does not prove that an imported board *behaves* like the
board that generated it — the file could be byte-perfect and still lose a corner, a
polygon-winding convention, or a units conversion, any of which a full-wave solve
would expose but a byte comparison would not. FS.3's remaining validation line item
is exactly this: "an imported reference board measures within tolerance of its
native-built twin." This is the gate that closes it.

## Decision

### 1. Outline→Layout twin path: reused, not rebuilt

`yee_export::gerber_to_layout(gerber, substrate, ports) -> Result<Layout,
GerberImportError>` **already exists** (shipped under FS.3.1) and is already the
helper the spec asked for: it parses copper regions via `gerber_to_polygons` and
wraps them in a `Layout` with caller-supplied `Substrate` and `Vec<PortRef>` —
Gerber carries neither. The studio's `import_gerber_impl`
(`studio/src-tauri/src/import.rs`) already calls this exact function. No new helper
was written; the gate imports and calls `yee_export::gerber_to_layout` directly, so
the gate and the studio's live import path share the identical helper by
construction, not by duplication. This is the documented API contract: **Gerber
carries geometry only — stackup (substrate ε_r/height) and ports are always
user-supplied at import time**, in the gate exactly as in the studio's
`ImportPanel`.

### 2. Gate `engine-import-twin-001`

New file: `crates/yee-engine/tests/import_twin.rs`.

**Native twin generator** (`native_stub_layout()`): rebuilds the exact trace
geometry `sparams_stub_notch.rs`'s `stub_job(true)` builds — a 3λ_g FR-4
microstrip line with a Hammerstad-open-end-corrected λ/4 open stub, the S.6
scenario named by the spec (chosen over `stackup_via.rs`'s FS.4.1 fixture, the
other `stub_notch` grep hit, because the spec names `sparams_stub_notch`
explicitly) — same constants (`EPS_R`, `H_M`, `W_M`, `F0_HZ`, `DX_M`), same polygon
math, so the "native" side of the twin is provably the same board the shipped
gate already validates.

**Measurement path**: rather than copy `sparams_stub_notch.rs`'s hand-rolled
`PortSpec`/`Pec`-boundary `JobSpec` construction, the gate reuses the R.5b library
builder (`yee_engine::board::two_port_board_job` + `reference_through_line` — the
same builder `engine_miter.rs` uses) so both the native `Layout` and the reimported
twin `Layout` are voxelized and measured through **one shared code path**, called
**independently** on each `Layout` (no grid/model sharing between the two runs —
each side derives its own grid from its own `Layout`, per the spec's "same fixture
options derived from each Layout independently" requirement).
`sparams::forward_transfer` splits forward/backward waves at each probe triple and
ratios DUT/reference — the ADR-0204-sanctioned launch-normalized double ratio,
immune to the single-ratio launch-inequality artifact `sparams_stub_notch.rs`'s
older `sparams::transmission_db` approach does not correct for.

**Gate structure**:
1. Build `native` Layout.
2. `layout_to_gerber(&native, ...)` → Gerber bytes → `gerber_to_layout(...)` → `twin` Layout.
3. **Structural assert**: polygon count, vertex count, and per-vertex coordinates
   within `0.5e-9` m (the `gerber-rt-001` house tolerance — half the 4.6
   fixed-point quantum, 1 nm — not raw `==`: a λ_g/4-derived stub length is
   generically irrational and not itself nanometre-aligned, so the reimported
   vertex is the nearest *representable* nanometre, not the pre-export float
   bit-for-bit).
4. **Measured assert**: `s21_lin(&native, ...)` and `s21_lin(&twin, ...)` each run
   a DUT + through-line-reference pair (4 FDTD solves total, uniform 0.3 mm grid,
   9000 steps — the `sparams_stub_notch` grid/step budget), then compare the two
   `|S21|` curves over a 3.0–6.2 GHz raster (65 bins).

## Root-cause reasoning for the expected (and measured) zero delta

`two_port_board_job`'s voxelizer rounds every coordinate to the nearest cell of a
0.3 mm grid — five to six orders of magnitude coarser than the ≤0.5 nm import
quantization. A sub-nanometre vertex perturbation divided by 0.3 mm is a
~2×10⁻⁶-cell fractional shift, far below any rounding boundary, so native and twin
voxelize to bit-identical grids (same PEC masks, same `eps_r_cells`, same
port/probe cell indices). The FDTD update itself has no RNG and no reduction-order
sensitivity (explicit per-cell stencil), so the two runs' probe series are then
also bit-identical. This reasoning was written into the test **before** running it
(module doc comment + assert failure message); the measured result confirmed it
exactly. This is the "root-cause a nonzero delta" requirement, inverted: the zero
delta is the *predicted*, not merely observed, outcome, and the comment records
why — satisfying the spec's honesty bar without there being an actual delta to
explain.

## Measured result

```
$ cargo test -p yee-engine --release --test import_twin -- --ignored --nocapture
engine-import-twin-001: native twin
engine-import-twin-001: imported twin
engine-import-twin-001: notch at 5.100 GHz | native -32.59 dB, imported twin -32.59 dB
  | max |Δ|S21|| across the band = 0.000e0 (linear)
test imported_board_measures_the_same_notch_as_its_native_twin ... ok
test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 293.30s
```

Bit-identical (`max |Δ|S21|| = 0.0` exactly, not merely "close") across all 65
frequency bins. Structural assert passed — all reimported vertices within the
0.5 nm tolerance. Notch depth −32.59 dB clears the ≥8 dB regression floor with
wide margin. Runtime 293.30 s (~4.9 min release), inside the spec's ~5 min target
(the `engine-miter-001` precedent for a 4-solve gate is 621 s on a graded fixture;
this uniform 0.3 mm fixture came in under that).

`#[ignore]`'d; confirmed picked up by blanket CI (and correctly skipped by
default) via `cargo test -p yee-engine` (debug, non-release).

## Tolerances pinned

Structural: `0.5e-9` m per-vertex (the `gerber-rt-001` house tolerance, unchanged,
reused rather than invented). Measured: bit-identical `|S21|` expected and
observed — no tolerance was needed because the delta is exactly zero, not merely
small. Notch-depth regression floor: ≥ 8 dB (both twins measured −32.59 dB, far
inside).

## Bit-exactness / regression discipline (unmodified gates, this commit)

The binding gate command — `cargo test -p yee-compute --release --test
graded_uniform_bitexact --test gpu_graded_parity --test gpu_cpu_parity --
--include-ignored` — stayed green (5/5, real GPU adapter `NVIDIA GeForce RTX 5060
Ti`, not SKIPPED). `cargo test -p yee-export --release` (all `gerber-rt-001/002/003`
+ arcs/flashes + `kicad-001/002` gates) stayed green, unmodified. The native stub
gate this work twins, `sparams_stub_notch.rs`
(`cargo test -p yee-engine --release --test sparams_stub_notch -- --ignored`), ran
unmodified and green: notch 4.850 GHz (−36.8 dB) vs 5.0 GHz theory, err 3.00 %,
110.41 s — no edits to `crates/yee-export/src/**` (the import/writer code was
reused verbatim) and no edits to `studio/**`. Workspace clippy (default +
`--no-default-features` on `yee-compute`) and `cargo fmt --check --all` clean
before every commit; `missing_docs` clean.

## Verdict

**GO, bit-identical.** The Gerber import chain — export → bytes → import →
`Layout` → voxelize → measure — reproduces the natively-built board's full-wave
S-parameters exactly, not approximately, for a vertex-exact rectangular-trace
fixture. This closes FS.3's remaining validation line item for the Gerber path.
The zero delta is not a coincidence of a lucky fixture: it follows from the
import's nanometre-scale quantization sitting five to six orders of magnitude
below the 0.3 mm voxelization grid, a relationship documented before the gate
was ever run and confirmed unchanged by the measurement.

## What remains (queued, not attempted here)

Per the spec's non-goals: **DXF import is the only remaining item on the FS.3
roadmap row** — arcs/flash boards as a twin fixture (rectangles prove the chain;
an arc-twin would need a tessellation-tolerance-vs-vertex-exactness argument that
adds nothing the structural assert here doesn't already cover) and studio wiring
(the studio already has import + echo, ADR-0209/FS.3.1c) are explicitly out of
scope and not queued as follow-ons.
