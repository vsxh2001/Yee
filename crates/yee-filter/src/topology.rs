//! Orderable-**board-topology auto-selector** (JLCPCB production track,
//! ADR-0167 brick **T3**).
//!
//! Two lumped band-pass topologies now reach the JLCPCB BOM/CPL path
//! ([`crate::jlcpcb_export`]), each orderable in a **complementary** spec regime:
//!
//! - the alternating series/shunt [`crate::LumpedLadder`]
//!   ([`synthesize_lumped`](crate::synthesize_lumped) → [`lumped_board`](crate::lumped_board)
//!   → [`join_placed_parts`](crate::join_placed_parts)) — orderable for **wideband**
//!   specs (ADR-0164: 1 GHz/70 % → zero blanks), blanks narrow-band (its *series*
//!   resonators want sub-pF/sub-nH parts); and
//! - the **top-C-coupled** [`crate::TopCNetwork`]
//!   ([`synthesize_top_c_coupled`](crate::synthesize_top_c_coupled) →
//!   [`top_c_board`](crate::top_c_board) → [`join_top_c_parts`](crate::join_top_c_parts))
//!   — orderable for the **sub-GHz / moderate-band** corner the ladder blanks in
//!   (ADR-0166: 0.5 GHz/20 % → zero blanks), blanks GHz-narrow (sub-pF coupling
//!   caps).
//!
//! The user-facing goal is "give a spec, get an orderable board" — the user
//! should not have to know which topology their spec needs. This module is the
//! **brain**: [`synthesize_orderable`] tries the topologies in a deterministic,
//! honest order and returns the one that yields a fully-orderable board, or an
//! honest "neither lumped topology is fully orderable → the distributed/planar
//! track" (it never fabricates an orderable board).
//!
//! Pure-compute, deterministic, **WASM-safe** (no I/O, network, threads, time,
//! RNG, or `unsafe`) — the same constraint as the rest of `yee-filter`.
//!
//! # `BoardTopology` vs [`crate::Topology`] (name distinction)
//!
//! This module's [`BoardTopology`] is the **board-realization** choice — *which
//! manufacturable lumped board* (alternating ladder vs top-C-coupled) yields the
//! orderable JLCPCB upload set. It is **deliberately distinct** from the
//! synthesis-realization [`crate::Topology`] enum (`{ CoupledResonator }`,
//! [`FilterProject::topology`](crate::FilterProject::topology)), which records the
//! *coupling-matrix* realization the synthesis produced. The two answer different
//! questions (synthesis form vs orderable board form), so they are separate types.
//!
//! # T4 follow-on (not in scope here)
//!
//! Wiring [`synthesize_orderable`] into `yee filter synth` (auto-route + report
//! the chosen topology) and the studio export stage is the ADR-0167 T4 follow-on
//! — two other lanes (`yee-cli`, `yee-studio-web`), a separate brick.

use yee_layout::Substrate;

use crate::{
    ESeries, FilterProject, Footprint, LumpedBoard, LumpedError, PlacedPart, join_placed_parts,
    join_top_c_parts, lumped_board, synthesize_lumped, synthesize_top_c_coupled, top_c_board,
};

/// The FR-4 reference substrate (εr 4.4, h 1.6 mm) — the project's reference
/// board, matching the `jlcpcb-orderable-001` / `top-c-board-001` gates.
///
/// The board geometry only affects the signal-line / pad placement, not the
/// component **values** the orderability decision turns on (those come from the
/// synthesized network), so a fixed reference substrate is the right default for
/// the value-driven topology choice; a per-call substrate is a trivial later
/// addition if board geometry ever needs to vary.
fn reference_substrate() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    }
}

/// The E-series the orderability check snaps component values to before
/// autopicking an LCSC part.
///
/// `E24` matches the series [`join_top_c_parts`](crate::join_top_c_parts) uses
/// internally, so the alternating-ladder and top-C arms snap on the **same**
/// grid and the two topologies' orderability is compared like-for-like.
const SELECTOR_SERIES: ESeries = ESeries::E24;

/// Which manufacturable **lumped board** topology [`synthesize_orderable`] chose.
///
/// This is the **board-realization** choice (see the [module docs](self#boardtopology-vs-cratetopology-name-distinction));
/// it is intentionally separate from the synthesis-realization
/// [`crate::Topology`] enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoardTopology {
    /// The alternating series/shunt [`crate::LumpedLadder`] (the conventional /
    /// simplest topology; orderable for **wideband** specs).
    AlternatingLadder,
    /// The **top-C-coupled** [`crate::TopCNetwork`] (`N` shunt LC resonators +
    /// `N+1` series coupling caps; orderable for the **sub-GHz / moderate-band**
    /// corner the ladder blanks in).
    TopCCoupled,
}

/// The selected orderable board: which [`BoardTopology`], its placed
/// [`LumpedBoard`], the joined [`PlacedPart`]s (ready for
/// [`jlcpcb_bom_csv`](crate::jlcpcb_bom_csv) / [`jlcpcb_cpl_csv`](crate::jlcpcb_cpl_csv)),
/// and whether every part resolved to a real LCSC #.
///
/// Produced by [`synthesize_orderable`]. When `fully_orderable` is `false` this
/// is still the honest best lumped board (the fewer-blanks topology) — its
/// `parts` carry the real blank set ([`PlacedPart::lcsc`] `== None`), an honest
/// "no lumped topology is fully orderable for this spec; the distributed/planar
/// track is the path" — never a fabricated orderable board.
#[derive(Debug, Clone, PartialEq)]
pub struct OrderableBoard {
    /// Which manufacturable board topology was chosen.
    pub topology: BoardTopology,
    /// The placed board (renderable [`yee_layout::Layout`] + [`crate::Placement`]
    /// list) for the chosen topology.
    pub board: LumpedBoard,
    /// The joined parts (placement → value → autopicked LCSC part), ready for the
    /// JLCPCB BOM/CPL emitters. A part with `lcsc == None` is an honest blank.
    pub parts: Vec<PlacedPart>,
    /// `true` iff **every** part resolved to a real LCSC part (zero blanks).
    pub fully_orderable: bool,
}

/// Count the parts that did **not** resolve to a real LCSC part (the blanks).
///
/// A part is a "blank" iff its [`PlacedPart::lcsc`] is `None` — no in-table Basic
/// LCSC part matched its kind/footprint/value. Zero blanks ⇔ a fully-orderable
/// board.
fn blanks(parts: &[PlacedPart]) -> usize {
    parts.iter().filter(|p| p.lcsc.is_none()).count()
}

/// Select the lumped board topology that yields an orderable JLCPCB board for
/// `project`, or honestly report that neither lumped topology can.
///
/// Deterministic policy (honest — never fabricates an orderable board):
///
/// 1. Try the **alternating ladder** ([`synthesize_lumped`](crate::synthesize_lumped)
///    → [`lumped_board`](crate::lumped_board) → [`join_placed_parts`](crate::join_placed_parts)).
///    If every part resolves to a real LCSC # → return
///    [`BoardTopology::AlternatingLadder`] with `fully_orderable = true`.
/// 2. Else try **top-C** ([`synthesize_top_c_coupled`](crate::synthesize_top_c_coupled)
///    → [`top_c_board`](crate::top_c_board) → [`join_top_c_parts`](crate::join_top_c_parts)),
///    at the **same** resolved order `n = project.prototype.order()` the ladder
///    uses (so the two topologies are compared like-for-like). If every part
///    resolves → return [`BoardTopology::TopCCoupled`] with
///    `fully_orderable = true`.
/// 3. Else (neither fully orderable) → return the topology with the **fewer
///    blanks** (the alternating ladder wins a tie, being the conventional /
///    simplest topology) with `fully_orderable = false` — an honest "no lumped
///    topology is fully orderable for this spec; the distributed/planar track is
///    the path." The returned `parts` carry the real blank set.
///
/// The alternating ladder is tried first (the conventional topology) so wideband
/// specs keep their existing board; top-C is the fallback that rescues the
/// narrow-band specs the ladder blanks on.
///
/// Both topologies are placed on the FR-4 [`reference_substrate`] and snapped to
/// the same [`SELECTOR_SERIES`] E-series, so the orderability comparison is
/// like-for-like. Pure-compute / WASM-safe. Use [`synthesize_orderable_on`] to
/// place the chosen board on a caller-supplied substrate (the topology decision
/// is identical — see that function's docs).
///
/// # Errors
///
/// Returns the [`LumpedError`] from [`synthesize_lumped`](crate::synthesize_lumped)
/// if the alternating ladder cannot be synthesized at all
/// ([`LumpedError::UnsupportedResponse`](crate::LumpedError::UnsupportedResponse)
/// for a non-band-pass response, or
/// [`LumpedError::OrderTooSmall`](crate::LumpedError::OrderTooSmall) for order
/// `N < 1`) — the same precondition the rest of the lumped track requires. (The
/// top-C synthesis shares those preconditions, so the ladder's error is the right
/// gate for both.)
pub fn synthesize_orderable(
    project: &FilterProject,
    footprint: Footprint,
) -> Result<OrderableBoard, LumpedError> {
    synthesize_orderable_on(project, &reference_substrate(), footprint)
}

/// Select the orderable lumped board topology for `project`, placing the chosen
/// board on the **caller-supplied** `substrate`.
///
/// Identical routing policy to [`synthesize_orderable`] (try the alternating
/// ladder, fall back to top-C, else the honest fewer-blanks board) — see that
/// function for the full policy and the `# Errors` contract. The only difference
/// is that both candidate boards are laid out on the given `substrate` rather
/// than the FR-4 [`reference_substrate`].
///
/// # Why the substrate does not change the topology decision
///
/// The topology **decision** and the BOM **orderability** turn only on the
/// component *values*, which come from the LC synthesis
/// ([`synthesize_lumped`](crate::synthesize_lumped) /
/// [`synthesize_top_c_coupled`](crate::synthesize_top_c_coupled)) and the E-series
/// snap ([`SELECTOR_SERIES`]) — neither depends on the substrate. The substrate
/// only affects the board **geometry** (the `Z0`-width signal-line trace from
/// [`microstrip_width`](crate::board) and therefore the pad/placement coordinates
/// the CPL reports). So `synthesize_orderable_on(p, &reference_substrate(), f)`
/// chooses the **same** [`BoardTopology`] and the **same** `fully_orderable` as
/// `synthesize_orderable(p, f)`; only `board.layout` / `board.placements` (the
/// geometry/CPL coords) differ. This is what lets the CLI honor the user's
/// `--eps-r`/`--h-mm` without changing which topology is orderable.
///
/// Pure-compute, deterministic, **WASM-safe** (no I/O, network, threads, time,
/// RNG, or `unsafe`).
///
/// # Errors
///
/// Same as [`synthesize_orderable`] — bubbles up the
/// [`LumpedError`](crate::LumpedError) from
/// [`synthesize_lumped`](crate::synthesize_lumped) when the spec is not a
/// realizable band-pass of order `N >= 1`.
pub fn synthesize_orderable_on(
    project: &FilterProject,
    substrate: &Substrate,
    footprint: Footprint,
) -> Result<OrderableBoard, LumpedError> {
    // ---- (1) alternating ladder -------------------------------------------
    // synthesize_lumped enforces the band-pass + order>=1 preconditions; bubble
    // its error up so a bad spec fails the same way the rest of the track does.
    let ladder = synthesize_lumped(project)?;
    let ladder_board = lumped_board(&ladder, substrate, footprint);
    let ladder_parts = join_placed_parts(
        &ladder_board.placements,
        &ladder,
        footprint,
        SELECTOR_SERIES,
    );
    let ladder_blanks = blanks(&ladder_parts);
    if ladder_blanks == 0 {
        return Ok(OrderableBoard {
            topology: BoardTopology::AlternatingLadder,
            board: ladder_board,
            parts: ladder_parts,
            fully_orderable: true,
        });
    }

    // ---- (2) top-C fallback (same resolved order as the ladder) -----------
    // Use the order the ladder already resolved so the two topologies cover the
    // same prototype — comparable, not a different N.
    let n = project.prototype.order();
    let net = synthesize_top_c_coupled(
        project.spec.approximation,
        n,
        project.spec.f0_hz,
        project.spec.fbw,
        project.spec.z0_ohm,
    );
    let top_c_brd = top_c_board(&net, substrate, footprint);
    let top_c_parts = join_top_c_parts(&top_c_brd.placements, &net);
    let top_c_blanks = blanks(&top_c_parts);
    if top_c_blanks == 0 {
        return Ok(OrderableBoard {
            topology: BoardTopology::TopCCoupled,
            board: top_c_brd,
            parts: top_c_parts,
            fully_orderable: true,
        });
    }

    // ---- (3) neither fully orderable → fewer-blanks (ladder on a tie) ------
    // Honest "distributed/planar track" pointer: return the closest-to-orderable
    // lumped board with its real blank set, fully_orderable = false.
    if ladder_blanks <= top_c_blanks {
        Ok(OrderableBoard {
            topology: BoardTopology::AlternatingLadder,
            board: ladder_board,
            parts: ladder_parts,
            fully_orderable: false,
        })
    } else {
        Ok(OrderableBoard {
            topology: BoardTopology::TopCCoupled,
            board: top_c_brd,
            parts: top_c_parts,
            fully_orderable: false,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Approximation, FilterSpec, Response, SpecMask, synthesize};

    /// A band-pass spec for the unit tests.
    fn bp_spec(f0_hz: f64, fbw: f64) -> FilterSpec {
        FilterSpec {
            response: Response::Bandpass,
            approximation: Approximation::Chebyshev { ripple_db: 0.5 },
            f0_hz,
            fbw,
            order: Some(3),
            z0_ohm: 50.0,
            mask: SpecMask {
                passband_ripple_db: 0.5,
                return_loss_db: 9.0,
                stopband: vec![(f0_hz * 1.5, 30.0)],
            },
        }
    }

    #[test]
    fn blanks_counts_none_lcsc() {
        // blanks() must count exactly the parts whose lcsc is None.
        let spec = bp_spec(2.0e9, 0.10); // narrow GHz → known to blank
        let proj = synthesize(&spec);
        let res = synthesize_orderable(&proj, Footprint::Smd0402).unwrap();
        let n_none = res.parts.iter().filter(|p| p.lcsc.is_none()).count();
        assert_eq!(blanks(&res.parts), n_none);
        // fully_orderable is the (blanks == 0) predicate, kept consistent.
        assert_eq!(res.fully_orderable, blanks(&res.parts) == 0);
    }

    #[test]
    fn wideband_picks_alternating_ladder() {
        // The ADR-0164 wideband fixture is fully orderable as an alternating
        // ladder, so the selector must pick it (arm 1, never reaching top-C).
        let spec = bp_spec(1.0e9, 0.70);
        let proj = synthesize(&spec);
        let res = synthesize_orderable(&proj, Footprint::Smd0402).unwrap();
        assert_eq!(res.topology, BoardTopology::AlternatingLadder);
        assert!(res.fully_orderable);
        assert_eq!(blanks(&res.parts), 0);
    }

    #[test]
    fn non_bandpass_errors() {
        // A non-band-pass spec cannot be realized by either lumped topology; the
        // selector bubbles up the ladder's UnsupportedResponse.
        let mut spec = bp_spec(1.0e9, 0.50);
        spec.response = Response::Lowpass;
        let proj = synthesize(&spec);
        assert_eq!(
            synthesize_orderable(&proj, Footprint::Smd0402),
            Err(LumpedError::UnsupportedResponse)
        );
    }

    #[test]
    fn top_c_order_matches_resolved_prototype_order() {
        // When top-C is reached, it must be synthesized at the SAME order the
        // ladder resolved (project.prototype.order()), so the two are comparable.
        // For a spec where neither is fully orderable, the returned board still
        // reflects that resolved order.
        let spec = bp_spec(2.0e9, 0.05); // GHz-narrow → neither orderable
        let proj = synthesize(&spec);
        let n = proj.prototype.order();
        let res = synthesize_orderable(&proj, Footprint::Smd0402).unwrap();
        // Whichever topology won the fewer-blanks tie, it carries the resolved
        // order's component count: ladder → 2N placements, top-C → 3N+1.
        let placements = res.board.placements.len();
        assert!(
            placements == 2 * n || placements == 3 * n + 1,
            "placements {placements} not consistent with resolved order N={n}"
        );
    }
}
