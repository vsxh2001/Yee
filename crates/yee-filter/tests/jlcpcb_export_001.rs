//! jlcpcb-export-001 (JLCPCB production track, ADR-0164 bricks **J2** + **J3**):
//! the placed lumped board → a valid JLCPCB **BOM CSV** + **CPL/centroid CSV**.
//!
//! **End-to-end + non-circular.** The component values come from the project
//! synthesis path (`synthesize` → `synthesize_lumped`); the placements from
//! [`lumped_board`]; the LCSC C-numbers from the independently-curated
//! [`yee_filter::LCSC_PARTS`] table (J1). The gate joins them by ref-des index
//! and checks the two CSVs are well-formed, schema-exact, and **designator-
//! consistent** (every placed component appears in both the CPL and the BOM, none
//! orphaned), and that the honest coverage holes (the narrow-band BPF series
//! resonator's sub-pF cap, below the 1 pF Basic floor) are **emitted with a blank
//! LCSC # + a flagged comment, not dropped**.
//!
//! For the standard 3-pole 0.5 dB Chebyshev BPF (f0 = 2 GHz, FBW = 0.10,
//! Z0 = 50 Ω) on **0603** footprints. Pure-compute, deterministic, fast (no
//! FDTD, no `#[ignore]`). Do NOT weaken: an honest blank-LCSC row for a
//! physically-unrealizable value is the correct outcome — do not invent a
//! C-number to fill it, and do not drop the line.

use yee_filter::{
    Approximation, BOM_HEADER, CPL_HEADER, ESeries, FilterSpec, Footprint, Response, SpecMask,
    jlcpcb_bom_csv, jlcpcb_cpl_csv, jlcpcb_files, join_placed_parts, lumped_board, synthesize,
    synthesize_lumped,
};
use yee_layout::Substrate;

/// Chebyshev 0.5 dB **N=3** bandpass spec (the standard demo filter; clone of the
/// J1 `jlcpcb_autopick_001` fixture).
fn fixture() -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: 2.0e9,
        fbw: 0.10,
        order: Some(3),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 9.0,
            stopband: vec![(2.4e9, 30.0)],
        },
    }
}

/// FR-4 substrate (εr 4.4, h 1.6 mm) — the project's reference board.
fn fr4() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    }
}

/// Is `lcsc` a well-formed LCSC C-number (`^C\d+$`)?
fn is_cnumber(lcsc: &str) -> bool {
    lcsc.starts_with('C') && lcsc.len() > 1 && lcsc[1..].chars().all(|c| c.is_ascii_digit())
}

/// Split a CSV line into fields, honoring double-quote escaping (RFC 4180): a
/// quoted field may contain commas; `""` inside a quoted field is a literal `"`.
fn parse_csv_line(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '"' => {
                if in_quotes && chars.peek() == Some(&'"') {
                    cur.push('"'); // escaped quote
                    chars.next();
                } else {
                    in_quotes = !in_quotes;
                }
            }
            ',' if !in_quotes => {
                fields.push(std::mem::take(&mut cur));
            }
            other => cur.push(other),
        }
    }
    fields.push(cur);
    fields
}

#[test]
fn jlcpcb_export_001() {
    let footprint = Footprint::Smd0603;
    let series = ESeries::E24;
    let spec = fixture();
    let proj = synthesize(&spec);
    assert_eq!(proj.prototype.order(), 3, "fixture is order N=3");

    let ladder = synthesize_lumped(&proj).expect("N=3 bandpass fixture should synthesize");
    let board = lumped_board(&ladder, &fr4(), footprint);
    let placements = &board.placements;
    // N=3 → 2·N = 6 placements (L1,C1,L2,C2,L3,C3).
    assert_eq!(placements.len(), 6, "N=3 → 2·N = 6 placed components");

    // The joined per-placement records (the crux: ref-des-index → resonator value
    // → autopick).
    let parts = join_placed_parts(placements, &ladder, footprint, series);
    assert_eq!(
        parts.len(),
        placements.len(),
        "every placement must value-join (none skipped)"
    );

    let files = jlcpcb_files(placements, &ladder, footprint, series);
    let bom_csv = jlcpcb_bom_csv(&parts);
    // jlcpcb_files must produce the same CSVs as the direct calls.
    assert_eq!(files.bom_csv, bom_csv, "jlcpcb_files BOM == direct BOM");
    assert_eq!(
        files.cpl_csv,
        jlcpcb_cpl_csv(placements),
        "jlcpcb_files CPL == direct CPL"
    );

    println!("\n=== JLCPCB BOM CSV (3-pole 0.5dB Cheb BPF, 2 GHz, 10% FBW, 0603) ===");
    println!("{bom_csv}");
    println!("\n=== JLCPCB CPL CSV ===");
    println!("{}", files.cpl_csv);

    // ---------------------------------------------------------------------
    // (A) BOM CSV schema
    // ---------------------------------------------------------------------
    let bom_lines: Vec<&str> = bom_csv.lines().collect();
    assert_eq!(
        bom_lines[0], BOM_HEADER,
        "BOM header must be exactly `{BOM_HEADER}`"
    );
    assert_eq!(
        bom_lines[0], "Comment,Designator,Footprint,LCSC Part #",
        "BOM header literal"
    );
    let bom_data = &bom_lines[1..];
    assert!(!bom_data.is_empty(), "BOM must have >= 1 data row");

    // Track the designators seen in the BOM and the realizable/unrealizable split.
    let mut bom_designators: Vec<String> = Vec::new();
    let mut realizable_rows = 0usize;
    let mut blank_rows = 0usize;
    for row in bom_data {
        let f = parse_csv_line(row);
        assert_eq!(f.len(), 4, "every BOM row has 4 fields: {row:?} → {f:?}");
        let (comment, designator, fp, lcsc) = (&f[0], &f[1], &f[2], &f[3]);
        assert_eq!(fp, "0603", "footprint column is the JLCPCB land name");
        assert!(!comment.is_empty(), "Comment must be non-empty");

        // Collect this row's designators (a row may group several).
        let row_des: Vec<String> = designator.split(',').map(|s| s.to_string()).collect();
        assert!(
            !row_des.is_empty() && !designator.is_empty(),
            "row has designators"
        );
        bom_designators.extend(row_des);

        if lcsc.is_empty() {
            // Unrealizable line: blank LCSC #, comment must flag it.
            blank_rows += 1;
            assert!(
                comment.contains("NO BASIC PART"),
                "a blank-LCSC row's Comment must flag the unfillable line: {comment:?}"
            );
        } else {
            // Realizable line: a well-formed C-number.
            realizable_rows += 1;
            assert!(
                is_cnumber(lcsc),
                "realizable LCSC # must be `C\\d+`: {lcsc:?}"
            );
            assert!(
                !comment.contains("NO BASIC PART"),
                "a realizable row must not be flagged: {comment:?}"
            );
        }
    }

    // The 3-pole BPF: the shunt resonator (realizable) gives the well-formed parts;
    // the series resonator's sub-pF cap is the honest blank-LCSC line.
    assert!(
        realizable_rows >= 1,
        "the realizable shunt L/C must produce >= 1 well-formed C-number row"
    );
    assert!(
        blank_rows >= 1,
        "the unrealizable sub-pF series cap must be PRESENT as a blank-LCSC row (not dropped)"
    );

    // ---------------------------------------------------------------------
    // (B) CPL CSV schema
    // ---------------------------------------------------------------------
    let cpl_lines: Vec<&str> = files.cpl_csv.lines().collect();
    assert_eq!(
        cpl_lines[0], CPL_HEADER,
        "CPL header must be exactly `{CPL_HEADER}`"
    );
    assert_eq!(
        cpl_lines[0], "Designator,Mid X,Mid Y,Layer,Rotation",
        "CPL header literal"
    );
    let cpl_data = &cpl_lines[1..];
    assert_eq!(
        cpl_data.len(),
        placements.len(),
        "one CPL row per placement"
    );

    // Board extent for the coordinate-bounds check.
    let bb = board.layout.bbox;
    let mut cpl_designators: Vec<String> = Vec::new();
    for row in cpl_data {
        let f = parse_csv_line(row);
        assert_eq!(f.len(), 5, "every CPL row has 5 fields: {row:?} → {f:?}");
        let (designator, mid_x, mid_y, layer, rotation) = (&f[0], &f[1], &f[2], &f[3], &f[4]);
        cpl_designators.push(designator.clone());
        assert_eq!(layer, "Top", "lumped board is single top layer");
        assert_eq!(rotation, "0", "axis-aligned placement → rotation 0");

        // Coordinates: `<num>mm`, finite, within the board bbox (± a footprint
        // margin, since center_m is the footprint centre and the bbox is over
        // pad corners — a centre is strictly inside the pad extent).
        let x_mm = mid_x.strip_suffix("mm").expect("Mid X has mm suffix");
        let y_mm = mid_y.strip_suffix("mm").expect("Mid Y has mm suffix");
        let x_m: f64 = x_mm.parse::<f64>().expect("Mid X numeric") * 1e-3;
        let y_m: f64 = y_mm.parse::<f64>().expect("Mid Y numeric") * 1e-3;
        assert!(x_m.is_finite() && y_m.is_finite(), "coords finite");
        // Centre lies within the board bounding box (centres are inside the pad
        // rectangles, which are inside the bbox).
        assert!(
            x_m >= bb.min.x - 1e-9 && x_m <= bb.max.x + 1e-9,
            "Mid X {x_m} out of board extent [{}, {}]",
            bb.min.x,
            bb.max.x
        );
        assert!(
            y_m >= bb.min.y - 1e-9 && y_m <= bb.max.y + 1e-9,
            "Mid Y {y_m} out of board extent [{}, {}]",
            bb.min.y,
            bb.max.y
        );
    }

    // ---------------------------------------------------------------------
    // (C) Designator consistency: the set of ref-des in the CPL == the union of
    //     ref-des across all BOM rows. Every placed part appears in both; none
    //     orphaned.
    // ---------------------------------------------------------------------
    let sorted_unique = |mut v: Vec<String>| -> Vec<String> {
        v.sort();
        v.dedup();
        v
    };
    let bom_set = sorted_unique(bom_designators.clone());
    let cpl_set = sorted_unique(cpl_designators.clone());
    // No duplicate designators within either file.
    assert_eq!(
        bom_set.len(),
        bom_designators.len(),
        "no designator appears twice across BOM rows"
    );
    assert_eq!(
        cpl_set.len(),
        cpl_designators.len(),
        "no designator appears twice in the CPL"
    );
    assert_eq!(
        bom_set, cpl_set,
        "BOM designators must equal CPL designators (every placed part in both, none orphaned)"
    );
    // And they equal the placement ref-des set (the ground truth).
    let placement_set = sorted_unique(placements.iter().map(|p| p.ref_des.clone()).collect());
    assert_eq!(
        cpl_set, placement_set,
        "CPL designators must equal the placement ref-des set"
    );

    println!(
        "\ncoverage: {realizable_rows} realizable BOM row(s) (well-formed C-number), \
         {blank_rows} honest blank-LCSC row(s); {} placements consistent across BOM+CPL",
        placements.len()
    );
}

/// Designator-grouping: two components that resolve to the SAME LCSC part group
/// into ONE BOM row listing both ref-des — exercised on the full 3-pole board
/// (the symmetric ladder makes the shunt resonators 1 & 3 identical, so their L
/// (and C) share a part and group).
#[test]
fn jlcpcb_export_001_grouping_on_board() {
    let footprint = Footprint::Smd0603;
    let proj = synthesize(&fixture());
    let ladder = synthesize_lumped(&proj).expect("synth");
    let board = lumped_board(&ladder, &fr4(), footprint);
    let parts = join_placed_parts(&board.placements, &ladder, footprint, ESeries::E24);
    let bom = jlcpcb_bom_csv(&parts);

    // Find a realizable BOM row whose Designator groups >= 2 ref-des. The N=3
    // shunt-first ladder is Shunt(1),Series(2),Shunt(3); resonators 1 & 3 are
    // identical by Chebyshev symmetry (g1 = g3), so L1==L3 and C1==C3 resolve to
    // the same part and MUST group.
    let mut found_grouped = false;
    for row in bom.lines().skip(1) {
        let f = parse_csv_line(row);
        let (designator, lcsc) = (&f[1], &f[3]);
        let count = designator.split(',').count();
        if !lcsc.is_empty() && count >= 2 {
            found_grouped = true;
            // The grouped designators are the symmetric twins (e.g. C1,C3).
            assert!(
                designator.contains(','),
                "grouped row Designator joins multiple ref-des: {designator:?}"
            );
        }
    }
    assert!(
        found_grouped,
        "the symmetric 3-pole ladder must group its identical shunt parts into a multi-ref-des BOM row"
    );
}
