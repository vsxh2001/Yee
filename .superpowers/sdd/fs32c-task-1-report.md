# FS.3.2c Task 1 report — import-twin measurement gate

**Branch:** `feature/fs3.2c-import-twin`
**Spec:** `docs/superpowers/specs/2026-07-24-fs3-2c-import-twin-design.md`
**Plan:** `docs/superpowers/plans/2026-07-24-fs3-2c-import-twin.md` — Task 1 only.

## What was read first

- `crates/yee-export/src/import.rs` / `lib.rs` — the Gerber writer
  (`layout_to_gerber`) and importer (`gerber_to_polygons`,
  `gerber_to_outline`, **`gerber_to_layout`**).
- `studio/src-tauri/src/import.rs` (read-only, out of lane) — the studio's
  `import_gerber_impl` command.
- `crates/yee-engine/tests/sparams_stub_notch.rs` (the S.6 fixture named
  by the spec) and `crates/yee-engine/tests/stackup_via.rs` (the other
  `stub_notch` grep hit, FS.4.1 — not used; the spec names
  `sparams_stub_notch` explicitly).
- `crates/yee-engine/src/board.rs` (R.5b) — `two_port_board_job` +
  `reference_through_line`, the library builder that replaced the
  gates' copy-pasted voxelize→JobSpec pattern.
- `crates/yee-engine/src/sparams.rs` and `crates/yee-engine/tests/engine_miter.rs`
  (the closest existing consumer of `two_port_board_job` for a
  notch-shaped DUT-vs-reference measurement) to confirm the
  ADR-0204-sanctioned launch-normalized double ratio
  (`sparams::forward_transfer` between the two probe triples, ratioed
  DUT/reference) is the right measurement, not the older single-probe
  `sparams::transmission_db` `sparams_stub_notch.rs` itself uses (that
  fixture predates `board.rs` and hand-rolls lumped `PortSpec` ports +
  a `Pec` boundary).

## Deliverable 1 — outline→Layout twin path: REUSED, not rebuilt

`yee_export::import::gerber_to_layout(gerber, substrate, ports) -> Result<Layout, GerberImportError>`
**already exists** and is **already** the helper the spec asked for
("a documented helper that rebuilds a `Layout` from an imported outline +
user-supplied stackup/port metadata"). It:

- parses copper regions via `gerber_to_polygons`,
- wraps them in a `Layout` with the caller-supplied `Substrate` and
  `Vec<PortRef>` (Gerber carries neither),
- returns `GerberImportError::NoCopper` if the file has no copper.

The studio's `import_gerber_impl` (`studio/src-tauri/src/import.rs`)
already calls this exact function. No new helper was written — Task 1's
gate imports and calls `yee_export::gerber_to_layout` directly, so the
gate and the studio's live import path share the identical helper by
construction (not by duplication).

## Deliverable 2 — gate `engine-import-twin-001`

New file: `crates/yee-engine/tests/import_twin.rs`.

**Native twin generator**: `native_stub_layout()` rebuilds the exact
trace geometry `sparams_stub_notch::stub_job(true)` builds (3λ_g FR-4
microstrip line + Hammerstad-open-end-corrected λ/4 stub) — same
constants (`EPS_R`, `H_M`, `W_M`, `F0_HZ`, `DX_M`), same polygon math.
Chosen over `stackup_via.rs`'s fixture per the spec's explicit naming.

**Measurement path**: rather than copy-pasting `sparams_stub_notch.rs`'s
hand-rolled `PortSpec`/`Pec`-boundary JobSpec construction, the gate
reuses the R.5b library builder (`yee_engine::board::two_port_board_job`
+ `reference_through_line`) — the same builder `engine_miter.rs` uses —
so both the native Layout and the reimported twin Layout are voxelized
and measured through **one shared code path**, called independently on
each Layout (no grid/model sharing). `sparams::forward_transfer` splits
forward/backward waves at each probe triple and ratios DUT/reference —
the ADR-0204-sanctioned double ratio, immune to the single-ratio launch
artifact `sparams_stub_notch.rs`'s older approach doesn't correct for.

**Gate structure**:
1. Build `native` Layout.
2. `layout_to_gerber(&native, ...)` → Gerber bytes → `gerber_to_layout(...)` → `twin` Layout.
3. **Structural assert**: polygon count, vertex count, and per-vertex
   coordinates within `0.5e-9` m (the `gerber-rt-001` house tolerance —
   half the 4.6 fixed-point quantum, 1 nm — not raw `==`: a λ_g/4-derived
   stub length is generically irrational and not itself nanometre-aligned,
   so the reimported vertex is the nearest *representable* nanometre, not
   the pre-export float bit-for-bit).
4. **Measured assert**: `s21_lin(&native, ...)` and `s21_lin(&twin, ...)`
   each run a DUT + through-line-reference pair (4 FDTD solves total,
   uniform 0.3 mm grid, 9000 steps — the `sparams_stub_notch` grid/step
   budget), then compare the two |S21| curves over a 3.0–6.2 GHz raster
   (65 bins).

## Root-cause reasoning for the expected (and measured) zero delta

`two_port_board_job`'s voxelizer rounds every coordinate to the nearest
cell of a 0.3 mm grid — five to six orders of magnitude coarser than the
≤0.5 nm import quantization. A sub-nanometre vertex perturbation divided
by 0.3 mm is a ~2×10⁻⁶-cell fractional shift, far below any rounding
boundary, so native and twin voxelize to bit-identical grids (same PEC
masks, same `eps_r_cells`, same port/probe cell indices). The FDTD update
itself has no RNG and no reduction-order sensitivity (explicit per-cell
stencil), so the two runs' probe series are then also bit-identical. This
was written into the test **before** running it (see the module doc
comment and the assert's failure message), and the measured result
confirmed it exactly — this is the "root-caused before pinning" case
requested for a *nonzero* delta, inverted: the zero delta is the
predicted, not merely observed, outcome, and the comment records why.

## Measured results (real, `--release --ignored`, this session)

```
$ cargo test -p yee-engine --release --test import_twin -- --ignored --nocapture
engine-import-twin-001: native twin
engine-import-twin-001: imported twin
engine-import-twin-001: notch at 5.100 GHz | native -32.59 dB, imported twin -32.59 dB
  | max |Δ|S21|| across the band = 0.000e0 (linear)
test imported_board_measures_the_same_notch_as_its_native_twin ... ok
test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 293.30s
```

Bit-identical (`max |Δ|S21|| = 0.0` exactly, not merely "close") across
all 65 frequency bins. Structural assert passed (all vertices within the
0.5 nm tolerance — in fact this fixture's coordinates print identical to
visual precision). Notch depth −32.59 dB clears the ≥8 dB regression
floor with wide margin. Runtime 293.30 s (~4.9 min release), inside the
~5 min target (the `engine-miter-001` precedent for a 4-solve gate is
621 s on a graded fixture; this uniform 0.3 mm fixture came in under
that).

Total wall time including the 4 solves: 293.30 s. `#[ignore]`'d;
confirmed picked up (and correctly skipped by default) via
`cargo test -p yee-engine` (debug, non-release): reports
`ignored, slow: 4 release FDTD solves...`.

## Files touched

- `crates/yee-engine/tests/import_twin.rs` — new gate.
- `crates/yee-engine/Cargo.toml` — added `yee-export` to
  `[dev-dependencies]` (test-only; no cyclic dependency — `yee-export`
  does not depend on `yee-engine`).
- `Cargo.lock` — regenerated by the above.

No edits to `crates/yee-export/src/**` (the import/writer code was
reused verbatim) and no edits to `studio/**` (read-only per lane).

## Verification run (all green, this session, real output)

1. `cargo clippy --workspace --all-targets -- -D warnings` → clean.
2. `cargo clippy -p yee-compute --all-targets --no-default-features -- -D warnings` → clean.
3. `cargo fmt --check --all` → clean (one file needed `cargo fmt -p yee-engine`
   before this, applied).
4. `cargo doc -p yee-engine --no-deps` → no missing_docs warnings.
5. `cargo test -p yee-compute --release --test graded_uniform_bitexact --test gpu_graded_parity --test gpu_cpu_parity -- --include-ignored`
   → all pass (bitexact 2/2, gpu_graded_parity 2/2, gpu_cpu_parity 1/1).
6. `cargo test -p yee-export --release` → all pass (gerber-rt-001/002/003/004,
   arcs/flashes, kicad-001/002, doctest).
7. `cargo test -p yee-engine --release --test sparams_stub_notch -- --ignored --nocapture`
   (the native stub gate this work twins) → unmodified, green: notch
   4.850 GHz (−36.8 dB) vs 5.0 GHz theory, err 3.00%, 110.41 s.
8. New gate itself: `cargo test -p yee-engine --release --test import_twin -- --ignored --nocapture`
   → green, 293.30 s (see above).
9. `cargo test -p yee-engine` (debug, default) → `import_twin` correctly
   `ignored` — blanket CI pickup confirmed, no default-run cost added.

No assertion was weakened anywhere. No GPU evidence was required for
this Task (no compute-path change; the bit-exact/GPU suite above was run
unmodified as the global-constraint regression check, all green on the
real NVIDIA GeForce RTX 5060 Ti path already wired into those tests).

## head_before / head_after

- head_before: `92e8dea5fc37fef979615dcc2d6966f43e59bb05`
- head_after: recorded after commit below.

## Task 2 (out of scope for this report)

Not started — Task 1 only, per the orchestrator's instruction. ADR-0229
+ roadmap row are Task 2's job.
