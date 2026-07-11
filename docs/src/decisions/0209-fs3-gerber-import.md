# ADR-0209: FS.3.0 — Gerber import: the writer dialect, byte-stable

**Status:** Accepted
**Date:** 2026-07-08
**Related:** ADR-0198 (the byte-checked exporter this round-trips),
FULL-SUITE-ROADMAP FS.3.
**Spec:** `docs/superpowers/specs/2026-07-08-fs3-gerber-import-design.md`

## Decision

`yee_export::import::gerber_to_polygons` parses exactly the RS-274X
region-fill subset our writer emits (format spec, `%MOMM%`, comments,
aperture bookkeeping, `G36/G37` contours from modal `X…Y…D01/D02` words)
into `Polygon`s with **exact** coordinates (`fixed46 → n·1e-9 m`), and
`gerber_to_layout` wraps them with caller-supplied substrate/ports
(Gerber carries neither — the studio import flow will prompt).
Everything outside the subset — inches, polarity, arcs, macro apertures,
stroked draws — is rejected with a named `GerberImportError` variant,
never silently mis-parsed.

## Gate `gerber-rt-001` (instant, non-ignored) — GREEN first run

On four real generator layouts (inset patch, quasi-Yagi, 2×1 array,
hairpin BPF): `import(export(L))` reproduces every polygon vertex-exactly,
and `export(import(export(L)))` is **byte-identical** to `export(L)` —
the strongest cheap proof of losslessness over the dialect, in the
ADR-0198 artifact philosophy. Rejection paths unit-pinned.

## FS.3.1a — outline import: SHIPPED

`gerber_to_outline` parses the Edge.Cuts stroked-path dialect
(`layout_to_gerber_outline`'s output): one move + draw chain, closing
vertex dropped; regions in an outline file rejected (a profile is a cut
path, not copper), and the copper importer keeps rejecting stroked draws
— the two dialects stay strictly apart. Gate `gerber-rt-002`: corners
equal bbox ± margin exactly on all four generator layouts.

## FS.3.1b — the studio import command: SHIPPED

`import_gerber` (command) / `import_gerber_impl` (pure core): copper
Gerber + optional outline + user-supplied substrate/ports → trace count,
bbox, SVG preview, outline corners, the layout as JSON (verify-flow
ready), and — the trust primitive — a **byte-provable echo**: the
response re-exports what was understood, so the UI can show
"round-trip: byte-identical" before the user runs a verify on an
imported board. `GerberImportError::NoCopper` added at the source (the
empty-file case previously panicked in `BBox::from_polygons`). Gate
`studio-import-e2e-001` (instant, in CI): echo byte-identical on the
A.1 patch export, outline enclosing the bbox, layout JSON deserializes;
error paths (no ports, imperial, no copper).

## Queued (FS.3.1c+)

The React import panel; D03 flashes, arcs, DXF; the FULL-SUITE gate (an
imported reference board measures within tolerance of its native twin).

## FS.3.1c addendum (2026-07-11): studio import panel

`ImportPanel` in `studio/src/App.tsx`: file pickers + paste area for the
copper (and optional Edge.Cuts) Gerber, stackup + single-port fields
(Gerber carries neither), one `import_gerber` call. The response renders
the parsed-layout SVG preview, polygon/bbox stats, layout-JSON and
copper-echo exports, and an **echo badge** — the UI face of
`studio-import-e2e-001`: green iff `echoIsLossless(input, echo)`, which
is strict byte equality (a trailing-newline difference is NOT lossless;
"byte-provable" means bytes, not semantics — pinned in the vitest gate).
DOM gates in `studio/src/import.test.tsx` (form shape; Import action
gated on non-empty copper). Flashes/arcs/DXF remain FS.3.2+.
