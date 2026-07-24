# FS.3.3 implementation plan — DXF import

Spec: `docs/superpowers/specs/2026-07-24-fs3-3-dxf-import-design.md`
Branch: `feature/fs3.3-dxf-import` (from `main` @ `ba8a53f`+)
Lane: `crates/yee-export/**`; Task 2 docs (`docs/src/decisions/0230-*.md`,
`docs/src/SUMMARY.md`, `FULL-SUITE-ROADMAP.md`).

## Global constraints (bind every task)

- After every functional commit, green + unmodified: `cargo test -p yee-export
  --release` (all existing gerber-rt/kicad gates) and `cargo test -p yee-engine
  --release --test import_twin` NOT required (5-min FDTD; run only if you touch
  shared outline types — then it IS required).
- Never weaken any assertion. No new dependencies (STOP-and-surface if you
  believe one is warranted). Typed rejections, each tested.
- clippy `-D warnings` workspace + `-p yee-compute --no-default-features`;
  `cargo fmt --check --all`; missing_docs; explicit-path staging; reports
  committed alongside work (`.superpowers/sdd/fs33-task-N-report.md`).
- Commit style: `yee-export:`/`docs:` prefix, ≤72-char imperative subject,
  why-body, trailers exactly:
Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01FFo41x449XDGJ7Xyds4L7M

## Task 1 — `import_dxf` + gates (yee-export)

- Read first: `crates/yee-export/src/import.rs` (Gerber parser structure, error
  enum style, outline types, arc tessellation helper + its 1 µm contract,
  fixture conventions in the gerber-rt tests).
- Implement the spec's subset: group-code pair scanner → entity iterator →
  LWPOLYLINE/POLYLINE-VERTEX closed chains → outline polygons; `$INSUNITS`
  (mm=4, inch=1 — reject others by name, incl. missing/0 unless the spec-chosen
  documented default applies); bulge → arc tessellation at the pinned 1 µm
  chord tolerance (reuse the Gerber helper if shape-compatible; else mirror,
  and say which in the report); layer filter option.
- Rejection matrix tests: open polyline, CIRCLE, ARC, ELLIPSE, SPLINE, TEXT,
  INSERT, nonzero-Z, bad units — one typed error each.
- Gate `dxf-rt-001` (test file `dxf_rt.rs` mirroring `gerber_rt.rs` naming):
  hand-authored DXF of the S.6 stub trace geometry → `dxf_to_outline` →
  vertex-exact match vs the native generator's polygons (0.5 nm tolerance);
  plus a bulge fixture with sagitta ≤ 1 µm assertions (CW+CCW).
- Verify per global constraints. One or two commits.

## Task 2 — ADR-0230 + roadmap row

- `docs/src/decisions/0230-fs33-dxf-import.md` per 0229's structure: subset
  boundary + rejection matrix, units decision, tessellation contract, gate
  numbers, **FS.3 COMPLETE** statement (import side: Gerber + DXF; export:
  Gerber/KiCad — DXF export explicitly not a goal), queued follow-on (studio
  DXF wiring). SUMMARY.md line after 0229; FS.3 row final update in
  FULL-SUITE-ROADMAP.md.
