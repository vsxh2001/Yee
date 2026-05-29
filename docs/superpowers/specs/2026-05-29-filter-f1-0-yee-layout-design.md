# Filter Phase F1.0 — `yee-layout` parametric geometry — Design Spec

**Phase:** F1.0 · **ADR:** ADR-0086 · **Date:** 2026-05-29 · **Status:** Accepted

## Goal
A pure-geometry crate that turns explicit physical dimensions into a parametric
microstrip coupled-resonator `Layout` (meshable/exportable later), with a
geometry-only gate. No EM, no new external dependency. The planar walking
skeleton for FILTER-DESIGN-ROADMAP Phase F1.

## API (yee-layout)
```rust
pub struct Substrate { pub eps_r: f64, pub height_m: f64,
    pub loss_tangent: f64, pub metal_thickness_m: f64 }        // serde
pub struct Point2 { pub x: f64, pub y: f64 }                   // metres, serde
pub struct Polygon { pub verts: Vec<Point2> }                  // top-metal footprint
pub struct PortRef { pub at: Point2, pub width_m: f64, pub ref_impedance_ohm: f64 }
pub struct BBox { pub min: Point2, pub max: Point2 }
pub struct Layout { pub substrate: Substrate, pub traces: Vec<Polygon>,
    pub ports: Vec<PortRef>, pub bbox: BBox }                  // serde
impl Layout { pub fn to_svg(&self) -> String; }                // dependency-free top view

// Hammerstad-Jensen microstrip synthesis (Pozar §3.8)
pub fn microstrip_width(z0_ohm: f64, eps_r: f64, h_m: f64) -> f64;  // -> W (m)
pub fn eps_eff(w_m: f64, h_m: f64, eps_r: f64) -> f64;

pub struct EdgeCoupledParams { pub substrate: Substrate, pub sections: Vec<EdgeCoupledSection>,
    pub feed_width_m: f64, pub feed_length_m: f64 }
pub struct EdgeCoupledSection { pub length_m: f64, pub width_m: f64, pub gap_m: f64 }
pub fn edge_coupled_bpf(p: &EdgeCoupledParams) -> Layout;

pub struct HairpinParams { pub substrate: Substrate, pub n: usize, pub arm_length_m: f64,
    pub line_width_m: f64, pub fold_spacing_m: f64, pub coupling_gap_m: f64,
    pub tap_offset_m: f64, pub feed_width_m: f64, pub feed_length_m: f64 }
pub fn hairpin_bpf(p: &HairpinParams) -> Layout;
```
`#![forbid(unsafe_code)]` + `#![warn(missing_docs)]` (manifest `[lints.rust]`
form, matching yee-synth/yee-core). Dep: `serde` only.

### Hammerstad-Jensen synthesis (microstrip_width)
```
A = z0/60·√((εr+1)/2) + (εr−1)/(εr+1)·(0.23 + 0.11/εr)
B = 377π/(2·z0·√εr)
W/h = 8·e^A/(e^{2A}−2)                                   if that ratio < 2
W/h = (2/π)[B−1−ln(2B−1) + (εr−1)/(2εr)·(ln(B−1)+0.39−0.61/εr)]  otherwise
eps_eff = (εr+1)/2 + (εr−1)/2 · (1+12·h/W)^(−1/2)
```

### Generators (geometry placement)
- `edge_coupled_bpf`: lay N coupled half-wave strips along x, each
  `length_m × width_m`, adjacent strips offset in y by `width + gap`, alternating
  the half-overlap stagger of an edge-coupled section; feed lines of
  `feed_width_m × feed_length_m` at the two ends; ports at the outer feed-line
  ends (`ref_impedance_ohm = 50`). Compute `bbox` from all polygons.
- `hairpin_bpf`: N U-folded resonators (two arms `arm_length_m` + a bend),
  spaced by `coupling_gap_m`, tapped-fed at `tap_offset_m`. Same port/bbox rules.

## DoD (machine-checkable)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-layout --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-layout` exit 0 (fast; pure geometry).
4. `geo-001` (`crates/yee-layout/tests/geo_001_edge_coupled.rs`): build a 2-section
   edge-coupled layout with known dims; assert `traces.len()` and `ports.len()==2`
   are as expected; `bbox` width/height within 1% of hand-computed extents; every
   `Polygon` has ≥4 verts and positive signed area (non-degenerate); the `Layout`
   round-trips through serde JSON unchanged; `to_svg()` returns non-empty SVG
   containing `<svg` and `</svg>`.
5. `geo-002` (`crates/yee-layout/tests/geo_002_hammerstad.rs`): `microstrip_width(
   50.0, 4.4, 1.6e-3)` within ±5% of `3.0e-3` m; `eps_eff(that_w, 1.6e-3, 4.4)`
   within ±5% of `3.3`. (Published HJ reference for FR-4 50 Ω.)
6. `cargo run`-able example or a `hairpin_bpf` test that emits an `.svg` artifact
   (optional: written under the crate's `target`/tmp, not committed).

## Out of scope
Any EM/FDTD meshing; coupling-matrix→dims synthesis (F1.2); KiCad/Gerber export
(F1.4); GUI. Geometry generation + HJ helpers only.
