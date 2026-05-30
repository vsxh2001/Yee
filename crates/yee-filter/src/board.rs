//! Lumped-LC **PCB board generator** (Filter Phase F2.2).
//!
//! Places a synthesized [`crate::LumpedLadder`]'s `L`/`C` resonators as **SMD
//! footprints (two copper pads each) on a board** and returns a
//! [`yee_layout::Layout`] (so the already-shipped Gerber / KiCad export renders
//! it) plus a [`Placement`] list (ref-des → footprint → position) for BOM
//! cross-reference (F2.1). Pure-geometry, WASM-safe, **no FDTD** — this is the
//! distributed track's [`crate::dimension_edge_coupled_layout`] analogue for the
//! lumped track, realizing the lumped-LC goal "to the pcb level."
//!
//! # Placement (walking skeleton)
//!
//! A horizontal **signal microstrip** runs left→right at the spec-`Z0` width
//! ([`yee_layout::microstrip_width`]); a **ground rail** (copper) runs along the
//! bottom edge at `y = 0`. Each [`crate::LcResonator`] gets its own x-slot
//! (left→right, ladder order) holding **two** footprints — an inductor `L`
//! then a capacitor `C` (ref-des `L1, C1, L2, C2, …`):
//!
//! - **Series branch** (series L–C): the two footprints sit **in-line on the
//!   signal line**, their pads strung along `x` bridging gaps in the through
//!   path. Both footprints are centred on the signal-line centreline `y_sig`.
//! - **Shunt branch** (parallel L–C): the two footprints sit on a short **stub**
//!   dropping from the signal line down to the ground rail, side by side in `x`,
//!   each with its pads stacked along `y` so the component bridges line→rail.
//!
//! Every pad is an axis-aligned copper [`yee_layout::Polygon::rect`] added to
//! [`yee_layout::Layout::traces`]; the signal-line segments and the ground rail
//! are copper rects too. The slot pitch is chosen so **no two pads overlap**
//! (the `lumped_pcb_001` gate's rect-disjoint check), with documented clearances.
//!
//! # Out of scope
//!
//! KiCad-native `(footprint)`/pad objects + courtyards/3D (F2.2b), solder
//! mask / silkscreen, auto-routing / impedance-matched meander, the FDTD sim
//! (F2.3), and the UI. Component **values** come from the F2.1 BOM, cross-
//! referenced by the [`Placement`] ref-des.

use serde::{Deserialize, Serialize};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, microstrip_width};

use crate::{LcBranch, LumpedLadder};

/// Standard SMD chip footprint (imperial code), selecting an IPC-7351 land
/// pattern via [`Footprint::pad`].
///
/// The three sizes span the common discrete L/C range; `Smd0603` is the default
/// the spec calls out.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Footprint {
    /// 0402 (1.0 × 0.5 mm body) chip.
    Smd0402,
    /// 0603 (1.6 × 0.8 mm body) chip — the F2.2 default.
    Smd0603,
    /// 0805 (2.0 × 1.25 mm body) chip.
    Smd0805,
}

/// The land-pattern geometry of one SMD [`Footprint`]: a two-pad layout.
///
/// All lengths are in metres. A footprint is two identical rectangular pads
/// whose centres are [`pitch_m`](PadSpec::pitch_m) apart along the **signal
/// direction**; each pad is [`pad_w_m`](PadSpec::pad_w_m) wide along that signal
/// direction and [`pad_len_m`](PadSpec::pad_len_m) long across it.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PadSpec {
    /// Pad extent **across** the signal direction (the longer land dimension),
    /// metres.
    pub pad_len_m: f64,
    /// Pad extent **along** the signal direction, metres.
    pub pad_w_m: f64,
    /// Centre-to-centre spacing of the footprint's two pads along the signal
    /// direction, metres.
    pub pitch_m: f64,
}

impl Footprint {
    /// The IPC-7351 (nominal "N", density level B) land pattern for this chip.
    ///
    /// Values are representative nominal SMD chip-resistor/capacitor lands from
    /// the IPC-7351B generator (e.g. as tabulated by common EDA libraries):
    /// each footprint is two pads of `pad_len_m × pad_w_m` whose centres are
    /// `pitch_m` apart. `pad_w_m` is the pad dimension *along* the signal path
    /// (the toe-to-heel direction); `pad_len_m` is *across* it. The exact land
    /// numbers vary by density level and house rules; these are sane, standard
    /// nominal values adequate for the F2.2 walking-skeleton placement (the gate
    /// validates non-overlap and topology, not land-pattern certification).
    pub fn pad(&self) -> PadSpec {
        match self {
            // 0402: ~0.6 mm pad pitch between centres, ~0.5×0.6 mm pads.
            Footprint::Smd0402 => PadSpec {
                pad_len_m: 0.60e-3,
                pad_w_m: 0.50e-3,
                pitch_m: 0.95e-3,
            },
            // 0603: ~1.6 mm pad pitch, ~0.8×0.9 mm pads.
            Footprint::Smd0603 => PadSpec {
                pad_len_m: 0.90e-3,
                pad_w_m: 0.80e-3,
                pitch_m: 1.60e-3,
            },
            // 0805: ~2.0 mm pad pitch, ~1.2×1.3 mm pads.
            Footprint::Smd0805 => PadSpec {
                pad_len_m: 1.30e-3,
                pad_w_m: 1.20e-3,
                pitch_m: 2.00e-3,
            },
        }
    }
}

/// Which ladder branch a placed component realizes.
///
/// Mirrors [`crate::LcBranch`] but is the board-domain enum recorded on every
/// [`Placement`] so a BOM/pick-and-place consumer needn't reach back into the
/// ladder.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BranchKind {
    /// Component sits in-line on the signal path (a series-branch L or C).
    Series,
    /// Component sits on a stub from the signal line to the ground rail (a
    /// shunt-branch L or C).
    Shunt,
}

/// One placed component: its ref-des, footprint, branch role, and board centre.
///
/// `center_m` is the geometric centre of the footprint's two-pad land, in
/// metres, in the board frame used by the returned [`Layout`]. Ref-des follow
/// `L1, C1, L2, C2, …` in ladder order (the inductor then the capacitor of each
/// resonator).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Placement {
    /// Reference designator, e.g. `"L1"` / `"C1"`.
    pub ref_des: String,
    /// The SMD land pattern used for this component.
    pub footprint: Footprint,
    /// Whether the component is series (in-line) or shunt (stub-to-ground).
    pub kind: BranchKind,
    /// Footprint centre `(x, y)`, metres, in the board frame.
    pub center_m: (f64, f64),
}

/// A lumped-LC board: the renderable [`Layout`] plus the [`Placement`] list.
///
/// `layout` carries the signal line, ground rail, and every component pad as
/// copper polygons (feedable straight into `layout_to_gerber` /
/// `layout_to_kicad_pcb`); `placements` is the ordered ref-des → footprint →
/// position list for BOM cross-reference. Produced by [`lumped_board`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LumpedBoard {
    /// The board geometry (substrate, copper traces/pads/rail, ports, bbox).
    pub layout: Layout,
    /// One entry per placed footprint, in ladder order (`L1, C1, L2, C2, …`).
    pub placements: Vec<Placement>,
}

/// Edge / inter-component clearance as a multiple of the footprint pad span.
///
/// A generous fixed clearance keeps the walking-skeleton placement trivially
/// non-overlapping with margin to spare; tuning it for density is F2.2b.
const CLEARANCE_FRAC: f64 = 0.6;

/// Place a synthesized [`LumpedLadder`] as SMD footprints on a board.
///
/// Lays a `Z0`-width signal microstrip along `x` with a ground rail along the
/// bottom, then walks the resonators left→right (ladder order). Each resonator
/// emits an **L** footprint then a **C** footprint (ref-des `L1, C1, L2, C2,
/// …`): **series** resonators in-line on the signal line, **shunt** resonators
/// on a stub from the line to the ground rail. See the [module docs](self) for
/// the placement rules.
///
/// Returns a [`LumpedBoard`] whose [`Layout`] holds the signal-line segments,
/// the ground rail, and every component pad as copper [`Polygon`]s, with ports
/// at the two ends of the signal line and a bounding box over all of it; and a
/// [`Placement`] per footprint. The component **values** are not encoded in the
/// geometry — they come from the F2.1 BOM keyed by ref-des.
///
/// Pure-geometry and deterministic; the x-slot pitch is derived from the
/// footprint's pad span plus a [`CLEARANCE_FRAC`] clearance so no two pads
/// overlap.
pub fn lumped_board(
    ladder: &LumpedLadder,
    substrate: &Substrate,
    footprint: Footprint,
) -> LumpedBoard {
    let pad = footprint.pad();
    let n = ladder.resonators.len();

    // Signal microstrip width from the spec Z0 (Hammerstad-Jensen).
    let w_line = microstrip_width(ladder.z0_ohm, substrate.eps_r, substrate.height_m);

    // --- Vertical layout (y) ------------------------------------------------
    // Ground rail along the bottom; signal line above it with room for a shunt
    // stub between. A shunt footprint stacks its two pads along y over a span of
    // `pitch_m + pad_w_m`; leave clearance above and below.
    let clearance = CLEARANCE_FRAC * (pad.pitch_m + pad.pad_w_m).max(w_line);
    let rail_h = w_line.max(pad.pad_len_m); // ground rail thickness
    let shunt_span_y = pad.pitch_m + pad.pad_w_m; // y-extent of a stacked footprint
    // Signal-line centreline. The shunt stub must CONNECT the signal line to the
    // ground rail: its two stacked pads span exactly from the rail top
    // (`rail_h`) up to the signal-line bottom (`y_sig - w_line/2`). So that span
    // equals `shunt_span_y`, giving `y_sig = rail_h + shunt_span_y + w_line/2`.
    // (No y-clearance in the stub — the pads are meant to abut line and rail, so
    // the shunt element's ground terminal actually reaches the ground plane.)
    let y_sig = rail_h + shunt_span_y + w_line / 2.0;

    // --- Horizontal layout (x) ----------------------------------------------
    // Each resonator owns a slot wide enough for two side-by-side footprints
    // (each spanning `pitch_m + pad_w_m` along its pad axis) with clearances.
    // Use the larger of the series-span (pads along x) and the shunt side-by-
    // side width so a slot fits either branch.
    let footprint_span = pad.pitch_m + pad.pad_w_m; // series: pad-axis span
    let shunt_pair_w = 2.0 * pad.pad_len_m + clearance; // shunt: two lands side by side in x
    let slot_w = (2.0 * footprint_span + clearance).max(shunt_pair_w + clearance) + clearance;

    let margin = slot_w; // leading/trailing signal-line lead-in
    let board_w = margin + (n as f64) * slot_w + margin;

    let mut traces: Vec<Polygon> = Vec::new();
    let mut placements: Vec<Placement> = Vec::with_capacity(2 * n);

    // Ground rail: full-width copper rect along the bottom.
    traces.push(Polygon::rect(0.0, 0.0, board_w, rail_h));

    // Signal line: emit as segments so series footprints can bridge gaps. We
    // record the gap intervals (along x) the series footprints occupy, then lay
    // the line as the complementary copper segments. Shunt branches keep the
    // line continuous (the stub drops off it).
    let mut series_gaps: Vec<(f64, f64)> = Vec::new();

    for (i, res) in ladder.resonators.iter().enumerate() {
        let slot_x0 = margin + (i as f64) * slot_w;
        let slot_cx = slot_x0 + slot_w / 2.0;
        let resnum = i + 1; // 1-based ref-des

        match res.branch {
            LcBranch::Series => {
                // L footprint then C footprint, in-line along x, centred on
                // y_sig. Two footprints sit symmetrically about the slot centre.
                let half = footprint_span / 2.0;
                let l_cx = slot_cx - (footprint_span / 2.0 + clearance / 2.0);
                let c_cx = slot_cx + (footprint_span / 2.0 + clearance / 2.0);
                for (cx, ref_des, _is_l) in [
                    (l_cx, format!("L{resnum}"), true),
                    (c_cx, format!("C{resnum}"), false),
                ] {
                    // Two pads along x, separated by pitch, centred at (cx, y_sig).
                    for s in [-1.0, 1.0] {
                        let pad_cx = cx + s * pad.pitch_m / 2.0;
                        traces.push(pad_rect_centered(pad_cx, y_sig, pad.pad_w_m, pad.pad_len_m));
                    }
                    placements.push(Placement {
                        ref_des,
                        footprint,
                        kind: BranchKind::Series,
                        center_m: (cx, y_sig),
                    });
                }
                // The two in-line footprints break the signal line over their
                // combined x-span; record that gap so the line skips it.
                let g0 = l_cx - half;
                let g1 = c_cx + half;
                series_gaps.push((g0, g1));
            }
            LcBranch::Shunt => {
                // L and C footprints side by side in x; each stacks its two pads
                // along y, bridging the signal line down to the ground rail. The
                // footprint centre sits midway in the stub span [rail_h ..
                // y_sig - w_line/2], so the bottom pad's lower edge touches the
                // rail top and the top pad's upper edge touches the line bottom —
                // electrically connecting the shunt element to ground.
                let stub_cy = rail_h + shunt_span_y / 2.0;
                let l_cx = slot_cx - (pad.pad_len_m / 2.0 + clearance / 2.0);
                let c_cx = slot_cx + (pad.pad_len_m / 2.0 + clearance / 2.0);
                for (cx, ref_des) in [(l_cx, format!("L{resnum}")), (c_cx, format!("C{resnum}"))] {
                    // Two pads stacked along y, separated by pitch.
                    for s in [-1.0, 1.0] {
                        let pad_cy = stub_cy + s * pad.pitch_m / 2.0;
                        traces.push(pad_rect_centered(cx, pad_cy, pad.pad_len_m, pad.pad_w_m));
                    }
                    placements.push(Placement {
                        ref_des,
                        footprint,
                        kind: BranchKind::Shunt,
                        center_m: (cx, stub_cy),
                    });
                }
            }
        }
    }

    // Lay the signal line as copper segments skipping the series gaps. Sort the
    // gaps and walk left→right; the centreline rect has height `w_line` centred
    // on `y_sig`.
    series_gaps.sort_by(|a, b| a.0.total_cmp(&b.0));
    let line_y0 = y_sig - w_line / 2.0;
    let mut cursor = 0.0;
    for (g0, g1) in &series_gaps {
        if *g0 > cursor {
            traces.push(Polygon::rect(cursor, line_y0, g0 - cursor, w_line));
        }
        cursor = cursor.max(*g1);
    }
    if cursor < board_w {
        traces.push(Polygon::rect(cursor, line_y0, board_w - cursor, w_line));
    }

    let bbox = BBox::from_polygons(&traces);

    let ports = vec![
        PortRef {
            at: Point2::new(0.0, y_sig),
            width_m: w_line,
            ref_impedance_ohm: ladder.z0_ohm,
        },
        PortRef {
            at: Point2::new(board_w, y_sig),
            width_m: w_line,
            ref_impedance_ohm: ladder.z0_ohm,
        },
    ];

    let layout = Layout {
        substrate: *substrate,
        traces,
        ports,
        bbox,
    };

    LumpedBoard { layout, placements }
}

/// An axis-aligned copper pad rectangle centred at `(cx, cy)` with full extents
/// `(w, h)` along `x` / `y`, metres.
fn pad_rect_centered(cx: f64, cy: f64, w: f64, h: f64) -> Polygon {
    Polygon::rect(cx - w / 2.0, cy - h / 2.0, w, h)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fr4() -> Substrate {
        Substrate {
            eps_r: 4.4,
            height_m: 1.6e-3,
            loss_tangent: 0.02,
            metal_thickness_m: 35e-6,
        }
    }

    #[test]
    fn pad_specs_are_two_pad_nonoverlapping_lands() {
        for fp in [Footprint::Smd0402, Footprint::Smd0603, Footprint::Smd0805] {
            let p = fp.pad();
            // Pads of one footprint must not overlap along the signal axis.
            assert!(
                p.pitch_m > p.pad_w_m,
                "{fp:?} pads overlap along signal axis"
            );
            assert!(p.pad_len_m > 0.0 && p.pad_w_m > 0.0);
        }
    }

    #[test]
    fn series_centered_on_line_shunt_drops_below() {
        let ladder = LumpedLadder {
            f0_hz: 2e9,
            fbw: 0.1,
            z0_ohm: 50.0,
            resonators: vec![
                crate::LcResonator {
                    branch: LcBranch::Shunt,
                    l_henry: 1e-9,
                    c_farad: 1e-12,
                },
                crate::LcResonator {
                    branch: LcBranch::Series,
                    l_henry: 1e-9,
                    c_farad: 1e-12,
                },
            ],
        };
        let board = lumped_board(&ladder, &fr4(), Footprint::Smd0603);
        assert_eq!(board.placements.len(), 4);
        let series_y: Vec<f64> = board
            .placements
            .iter()
            .filter(|p| p.kind == BranchKind::Series)
            .map(|p| p.center_m.1)
            .collect();
        let shunt_y: Vec<f64> = board
            .placements
            .iter()
            .filter(|p| p.kind == BranchKind::Shunt)
            .map(|p| p.center_m.1)
            .collect();
        // Series sit on the signal line; shunt sit below it (toward y=0).
        for &sy in &series_y {
            for &hy in &shunt_y {
                assert!(hy < sy, "shunt {hy} should be below series {sy}");
            }
        }
    }
}
