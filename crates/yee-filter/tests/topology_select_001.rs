//! topology-select-001 (JLCPCB production track, ADR-0167 brick **T3** — the
//! orderable-topology auto-selector): for a given spec,
//! [`yee_filter::synthesize_orderable`] returns the lumped board topology that
//! yields a **fully-orderable** JLCPCB board (or an honest "neither lumped
//! topology can → distributed/planar"), choosing between the alternating
//! [`LumpedLadder`](yee_filter::LumpedLadder) and the
//! [`TopCNetwork`](yee_filter::TopCNetwork).
//!
//! # What this proves (and why it is non-circular)
//!
//! Two complementary lumped topologies reach the JLCPCB BOM path, each orderable
//! in a different spec regime (`jlcpcb-orderable-001` — wideband ladder;
//! `top-c-board-001` — sub-GHz/moderate-band top-C). The selector routes a spec
//! to the right one. **Every assertion here runs the REAL selector**
//! ([`synthesize_orderable`]) and then **independently recomputes** the actual
//! orderability of *both* topologies from the lower-level synthesis/board/join
//! paths ([`synthesize_lumped`] → [`lumped_board`] → [`join_placed_parts`] and
//! [`synthesize_top_c_coupled`] → [`top_c_board`] → [`join_top_c_parts`]),
//! counting the parts whose [`PlacedPart::lcsc`](yee_filter::PlacedPart) is
//! `None`. The selector's returned `topology` / `fully_orderable` must agree with
//! that independent ground truth — a pick can never be "right by construction".
//!
//! # The three discriminating cases (all Chebyshev 0.5 dB, N = 3, Z0 = 50 Ω)
//!
//! 1. **Wideband — 1.0 GHz / 70 % / 0402** (the ADR-0164 fixture): the
//!    alternating ladder is fully orderable, so the selector picks
//!    [`BoardTopology::AlternatingLadder`] with `fully_orderable == true` (arm 1
//!    short-circuits — top-C is never reached).
//!
//! 2. **The discriminating cell — 0.5 GHz / 20 % / 0402** (THE LOAD-BEARING
//!    PROOF, found by an empirical (f0, FBW, footprint) grid sweep — see below):
//!    the alternating ladder **blanks 4 of its 6 parts** but the top-C network is
//!    **fully orderable (0 blanks)**, so the selector's top-C fallback rescues a
//!    real spec the ladder cannot make orderable →
//!    [`BoardTopology::TopCCoupled`] with `fully_orderable == true`. (This is the
//!    same envelope cell `top-c-board-001` pins as orderable, here reached via the
//!    auto-selector rather than by hand.)
//!
//!    The empirically-found discriminating cell, and why the ladder blanks it:
//!    ```text
//!    ladder (4/6 blank): L1=L3 2.0 nH (< inductor floor), L2 91 nH (> ceiling),
//!                        C2 1.2 pF (< cap floor); only C1=C3 51 pF resolve.
//!    top-C  (0/10):      L1..L3 16 nH → C27143; C1..C3 3.3/4.3 pF; coupling
//!                        Cc 2.4/1.0 pF → C1559/C1550. Every part orderable.
//!    ```
//!    The grid swept f0 ∈ {0.3,0.5,0.7,1.0} GHz × FBW ∈ {10,15,20,30}% ×
//!    {0402,0603}; several cells discriminate (0.3/20, 0.3/30, 0.5/20, 0.5/30,
//!    0.7/30) — 0.5 GHz / 20 % / 0402 is chosen as the cleanest (largest ladder-
//!    blank margin: 4/6) and matches the ADR-0166 envelope.
//!
//! 3. **GHz-narrow — 2.0 GHz / 5 % / 0402**: NEITHER lumped topology is fully
//!    orderable (ladder blanks 4/6, top-C blanks 4/10 — both have sub-floor
//!    elements at this narrow GHz band), so `fully_orderable == false` and the
//!    selector returns the **fewer-blanks** topology (ladder on the count tie).
//!    The honest "distributed/planar track is the path" — with a **real, non-empty
//!    blank set**, not a fabricated orderable board.
//!
//! Pure-compute, deterministic, fast (no FDTD, no `#[ignore]`). **Do NOT weaken**:
//! case 2 is the load-bearing proof that the top-C fallback strictly rescues a
//! spec the ladder blanks; if a future table/synthesis edit makes the ladder
//! orderable there (or top-C blank), this gate must FAIL (that is the point) — do
//! not relax it and do not invent C-numbers.

use yee_filter::{
    Approximation, BoardTopology, CompKind, ESeries, FilterProject, FilterSpec, Footprint,
    PlacedPart, Response, SpecMask, join_placed_parts, join_top_c_parts, lumped_board, synthesize,
    synthesize_lumped, synthesize_orderable, synthesize_top_c_coupled, top_c_board,
};
use yee_layout::Substrate;

/// FR-4 substrate (εr 4.4, h 1.6 mm) — the project's reference board, matching
/// `jlcpcb-orderable-001` / `top-c-board-001` and the selector's internal default.
fn fr4() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    }
}

/// A band-pass spec (Chebyshev 0.5 dB, N = 3, Z0 = 50 Ω) at the given centre and
/// fractional bandwidth.
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

/// Count the parts that did NOT resolve to a real LCSC part (the blanks).
fn blanks(parts: &[PlacedPart]) -> usize {
    parts.iter().filter(|p| p.lcsc.is_none()).count()
}

/// The blank count of the **alternating ladder** for this project — recomputed
/// INDEPENDENTLY of the selector (the ground truth the selector is checked
/// against). Uses the same E24 series and FR-4 board the selector uses.
fn ladder_blanks(proj: &FilterProject, footprint: Footprint) -> (usize, usize) {
    let ladder = synthesize_lumped(proj).expect("band-pass N>=1 ladder synthesizes");
    let board = lumped_board(&ladder, &fr4(), footprint);
    let parts = join_placed_parts(&board.placements, &ladder, footprint, ESeries::E24);
    (blanks(&parts), parts.len())
}

/// The blank count of the **top-C network** for this project — recomputed
/// INDEPENDENTLY of the selector, at the SAME resolved order the ladder uses.
fn top_c_blanks(proj: &FilterProject, footprint: Footprint) -> (usize, usize) {
    let n = proj.prototype.order();
    let net = synthesize_top_c_coupled(
        proj.spec.approximation,
        n,
        proj.spec.f0_hz,
        proj.spec.fbw,
        proj.spec.z0_ohm,
    );
    let board = top_c_board(&net, &fr4(), footprint);
    let parts = join_top_c_parts(&board.placements, &net);
    (blanks(&parts), parts.len())
}

/// Human label for a value (pF / nH).
fn label(kind: CompKind, v: f64) -> String {
    match kind {
        CompKind::Capacitor => format!("{:.3} pF", v * 1e12),
        CompKind::Inductor => format!("{:.3} nH", v * 1e9),
    }
}

/// Print a selector result's parts (ref-des → value → LCSC #).
fn print_parts(tag: &str, parts: &[PlacedPart]) {
    println!(
        "  {tag} parts ({} total, {} blank):",
        parts.len(),
        blanks(parts)
    );
    for p in parts {
        let lcsc = p.lcsc.map(|x| x.lcsc).unwrap_or("(NO BASIC PART)");
        println!(
            "    {:<5} {:?} {} -> {}",
            p.ref_des,
            p.kind,
            label(p.kind, p.chosen_value),
            lcsc
        );
    }
}

/// `fully_orderable` must be exactly the (every part resolved) predicate, and the
/// returned board's placements must match the returned parts one-for-one.
fn assert_internally_consistent(res: &yee_filter::OrderableBoard) {
    assert_eq!(
        res.fully_orderable,
        blanks(&res.parts) == 0,
        "fully_orderable must be exactly (zero blanks)"
    );
    assert_eq!(
        res.parts.len(),
        res.board.placements.len(),
        "every placement must value-join into a part (none dropped)"
    );
    // Every part with an LCSC # must carry a well-formed C-number of the right
    // kind (a resolved pick is genuinely real, not a placeholder).
    for p in &res.parts {
        if let Some(part) = &p.lcsc {
            assert_eq!(part.kind, p.kind, "{} picked kind must match", p.ref_des);
            assert!(
                part.lcsc.starts_with('C')
                    && part.lcsc.len() > 1
                    && part.lcsc[1..].chars().all(|c| c.is_ascii_digit()),
                "{} resolved to a non-C-number {:?}",
                p.ref_des,
                part.lcsc
            );
        }
    }
}

#[test]
fn topology_select_001() {
    let footprint = Footprint::Smd0402;

    println!("\n========================================================================");
    println!("topology-select-001 — orderable board-topology auto-selector (T3, ADR-0167)");
    println!("========================================================================");

    // =====================================================================
    // CASE 1 — WIDEBAND (1.0 GHz / 70 %): the alternating ladder is fully
    // orderable, so the selector picks it (arm 1 short-circuits before top-C).
    // =====================================================================
    println!(
        "\n--- CASE 1: WIDEBAND 1.0 GHz / 70% / 0402 → expect AlternatingLadder, orderable ---"
    );
    let wide = synthesize(&bp_spec(1.0e9, 0.70));
    let wide_res = synthesize_orderable(&wide, footprint).expect("wideband BPF must synthesize");

    // Independent ground truth: the ladder really IS fully orderable here.
    let (wide_ladder_blank, wide_ladder_tot) = ladder_blanks(&wide, footprint);
    println!("  [independent] alternating ladder blanks = {wide_ladder_blank}/{wide_ladder_tot}");
    print_parts("selector chose", &wide_res.parts);
    assert_eq!(
        wide_ladder_blank, 0,
        "ground truth: the wideband ladder must itself be fully orderable"
    );
    // The selector must agree: AlternatingLadder + fully_orderable.
    assert_eq!(
        wide_res.topology,
        BoardTopology::AlternatingLadder,
        "wideband spec must route to the alternating ladder"
    );
    assert!(
        wide_res.fully_orderable,
        "wideband spec must yield a fully-orderable board (zero blanks)"
    );
    assert_eq!(
        blanks(&wide_res.parts),
        0,
        "wideband selector parts must have ZERO blanks"
    );
    assert_internally_consistent(&wide_res);
    // The selector's parts must equal the independent ladder join (it really took
    // the ladder arm, not some other path).
    assert_eq!(
        wide_res.parts.len(),
        wide_ladder_tot,
        "selector parts count must equal the independent ladder join (2N parts)"
    );
    println!("  PASS: wideband → AlternatingLadder, fully orderable (0 blanks).");

    // =====================================================================
    // CASE 2 — THE DISCRIMINATING CELL (0.5 GHz / 20 %): the ladder BLANKS but
    // top-C is fully orderable. The load-bearing proof the fallback rescues a
    // spec the ladder can't make orderable.
    // =====================================================================
    println!(
        "\n--- CASE 2: DISCRIMINATING 0.5 GHz / 20% / 0402 → expect TopCCoupled, orderable ---"
    );
    let disc = synthesize(&bp_spec(0.5e9, 0.20));
    let disc_res = synthesize_orderable(&disc, footprint).expect("discriminating BPF synthesizes");

    // Independent ground truth: ladder BLANKS (> 0), top-C is fully orderable (0).
    let (disc_ladder_blank, disc_ladder_tot) = ladder_blanks(&disc, footprint);
    let (disc_top_c_blank, disc_top_c_tot) = top_c_blanks(&disc, footprint);
    println!(
        "  [independent] alternating ladder blanks = {disc_ladder_blank}/{disc_ladder_tot} \
         (the ladder CANNOT make this orderable)"
    );
    println!(
        "  [independent] top-C blanks = {disc_top_c_blank}/{disc_top_c_tot} \
         (top-C rescues it — fully orderable)"
    );
    print_parts("selector chose", &disc_res.parts);

    // THE LOAD-BEARING GROUND TRUTH: the two topologies genuinely DISAGREE here.
    assert!(
        disc_ladder_blank > 0,
        "DISCRIMINATING PROOF: the alternating ladder MUST blank at this cell \
         (else the fallback is not exercised) — got {disc_ladder_blank} blanks"
    );
    assert_eq!(
        disc_top_c_blank, 0,
        "DISCRIMINATING PROOF: top-C MUST be fully orderable at this cell \
         (else the fallback does not rescue it) — got {disc_top_c_blank} blanks"
    );
    // The selector must therefore route to top-C and report it orderable.
    assert_eq!(
        disc_res.topology,
        BoardTopology::TopCCoupled,
        "the discriminating cell must route to the top-C fallback (ladder blanks here)"
    );
    assert!(
        disc_res.fully_orderable,
        "the top-C fallback must yield a fully-orderable board at the discriminating cell"
    );
    assert_eq!(
        blanks(&disc_res.parts),
        0,
        "the chosen top-C board must have ZERO blanks"
    );
    assert_internally_consistent(&disc_res);
    // It really took the top-C arm: 3N+1 parts (vs the ladder's 2N).
    assert_eq!(
        disc_res.parts.len(),
        disc_top_c_tot,
        "selector parts count must equal the independent top-C join (3N+1 parts)"
    );
    assert_ne!(
        disc_top_c_tot, disc_ladder_tot,
        "top-C (3N+1) and ladder (2N) part counts differ — confirms the arm taken"
    );
    println!(
        "  PASS: discriminating 0.5 GHz/20% → ladder blanks {disc_ladder_blank}/{disc_ladder_tot}, \
         top-C 0/{disc_top_c_tot}; selector chose TopCCoupled (fully orderable). \
         The fallback rescued a real spec the ladder can't."
    );

    // =====================================================================
    // CASE 3 — GHz-NARROW (2.0 GHz / 5 %): NEITHER lumped topology is fully
    // orderable. Honest fully_orderable == false; the fewer-blanks topology with
    // a REAL, non-empty blank set.
    // =====================================================================
    println!("\n--- CASE 3: GHz-NARROW 2.0 GHz / 5% / 0402 → expect fully_orderable == false ---");
    let narrow = synthesize(&bp_spec(2.0e9, 0.05));
    let narrow_res = synthesize_orderable(&narrow, footprint).expect("narrow BPF synthesizes");

    // Independent ground truth: BOTH topologies blank here.
    let (narrow_ladder_blank, narrow_ladder_tot) = ladder_blanks(&narrow, footprint);
    let (narrow_top_c_blank, narrow_top_c_tot) = top_c_blanks(&narrow, footprint);
    println!(
        "  [independent] alternating ladder blanks = {narrow_ladder_blank}/{narrow_ladder_tot}"
    );
    println!("  [independent] top-C blanks = {narrow_top_c_blank}/{narrow_top_c_tot}");
    print_parts("selector chose", &narrow_res.parts);
    assert!(
        narrow_ladder_blank > 0 && narrow_top_c_blank > 0,
        "ground truth: NEITHER topology is fully orderable at the GHz-narrow cell \
         (ladder {narrow_ladder_blank}, top-C {narrow_top_c_blank} blanks)"
    );

    // The selector must report NOT fully orderable.
    assert!(
        !narrow_res.fully_orderable,
        "GHz-narrow spec must NOT be fully orderable (neither lumped topology resolves)"
    );
    // The blank set is REAL and non-empty (not a fabricated orderable board).
    let narrow_chosen_blanks = blanks(&narrow_res.parts);
    assert!(
        narrow_chosen_blanks > 0,
        "the not-orderable board must carry a REAL, non-empty blank set (got {narrow_chosen_blanks})"
    );
    // It is the FEWER-blanks topology (ladder wins the count tie). Cross-check the
    // returned topology against the independently-computed fewer-blanks rule.
    let expected_topology = if narrow_ladder_blank <= narrow_top_c_blank {
        BoardTopology::AlternatingLadder
    } else {
        BoardTopology::TopCCoupled
    };
    assert_eq!(
        narrow_res.topology, expected_topology,
        "GHz-narrow: the chosen topology must be the fewer-blanks one (ladder on a tie)"
    );
    // And the chosen board's blanks must be the minimum of the two.
    assert_eq!(
        narrow_chosen_blanks,
        narrow_ladder_blank.min(narrow_top_c_blank),
        "the chosen not-orderable board must have the fewer-blanks count"
    );
    assert_internally_consistent(&narrow_res);
    println!(
        "  PASS: GHz-narrow 2.0 GHz/5% → ladder {narrow_ladder_blank}/{narrow_ladder_tot}, \
         top-C {narrow_top_c_blank}/{narrow_top_c_tot}; selector chose {:?}, fully_orderable=false, \
         {narrow_chosen_blanks} real blank(s). Distributed/planar track is the path.",
        narrow_res.topology
    );

    println!("\n========================================================================");
    println!(
        "topology-select-001 PASS: wideband→ladder(orderable), 0.5GHz/20%→top-C(orderable, \
         rescued the ladder's blanks), 2GHz/5%→neither(honest distributed pointer)."
    );
    println!("========================================================================\n");
}
