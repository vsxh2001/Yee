# FS.3.0 implementation plan — Gerber region-subset import

**Spec:** `docs/superpowers/specs/2026-07-08-fs3-gerber-import-design.md`
**Lane:** `crates/yee-export/**` (+ docs). **ADR:** 0209.

1. `yee_export::import`: `gerber_to_polygons` (region subset, modal
   coordinates, exact `fixed46 → n·1e-9 m`), `gerber_to_layout`
   (caller-supplied substrate/ports), `GerberImportError` with one
   variant per explicit rejection.
2. Gate `gerber-rt-001` (`tests/gerber_roundtrip_import.rs`, instant,
   non-ignored → runs in the workspace lint-test CI job automatically):
   vertex-exact import + `export∘import∘export` byte-identical on the
   inset patch, quasi-Yagi, 2×1 array, and hairpin BPF generator
   layouts; rejection paths (imperial, polarity, arcs, stroked draws,
   draw-before-move, unclosed region).
3. Verification: fmt + clippy floor + `cargo test -p yee-export`.
