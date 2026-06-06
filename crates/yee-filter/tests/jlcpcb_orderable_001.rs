//! jlcpcb-orderable-001 (JLCPCB production track, ADR-0164 brick **J5** — the
//! orderability capstone): a *lumped-appropriate* filter spec yields a
//! **FULLY-ORDERABLE** JLCPCB upload set — every part autopicks to a real LCSC
//! **Basic** C-number (**ZERO** blank rows), with a schema-valid BOM + CPL.
//!
//! # What this proves (and the honest bound)
//!
//! The earlier J1/J2/J3 gates ([`jlcpcb_autopick_001`](../../tests/jlcpcb_autopick_001.rs),
//! [`jlcpcb_export_001`](../../tests/jlcpcb_export_001.rs)) run the **standard
//! NARROW-band** demo filter (2 GHz, **10 %** FBW) and correctly show ~half its
//! BOM is unrealizable: a narrow-band band-pass shrinks the *series*-resonator
//! capacitance to **sub-pF** (below the 1 pF Basic floor) and the *shunt*
//! inductor below the ~2.7 nH Basic floor, so those rows are emitted as **honest
//! blank-LCSC** lines (`… (NO BASIC PART)`), not faked.
//!
//! This gate proves the **complementary** claim: for the regime where a lumped
//! ladder *is* the right realization — a **WIDEBAND** lumped BPF — the pipeline
//! produces a **complete, orderable board with NO blanks**. The realizability
//! math (see [`yee_filter::jlcpcb`] coverage note + the low-pass→band-pass
//! transform in [`yee_filter::lumped`]):
//!
//! - series-branch resonator: `L = g·Z0/(ω0·Δ)`, `C = Δ/(ω0·Z0·g)`
//! - shunt-branch  resonator: `L = Z0·Δ/(ω0·g)`,  `C = g/(ω0·Z0·Δ)`
//!
//! As the fractional bandwidth `Δ` grows, the series cap (`∝ Δ`) and the shunt
//! inductor (`∝ Δ`) both rise off their respective floors while the series
//! inductor (`∝ 1/Δ`) and shunt cap (`∝ 1/Δ`) fall off their ceilings — so a
//! sufficiently **wideband** BPF lands EVERY L and C inside the autopick table's
//! realizable decades (caps 1–330 pF, inductors 2.7–68 nH).
//!
//! **The binding constraint is the *shunt inductor* floor** (the table's lowest
//! inductor: 2.7 nH on 0402, 4.7 nH on 0603), not the series cap as the J1 docs'
//! narrow-band note emphasizes — at wide FBW the series cap clears 1 pF easily,
//! but the smallest shunt L (`Z0·Δ/(ω0·g1)`) is what must reach the inductor
//! floor *with a stocked part within the autopick band*. This is why the chosen
//! fixture uses **0402** (its 2.7/3.9 nH parts catch the ~3.6 nH shunt L; the
//! 0603 grid starts at 4.7 nH and would blank that value) — documented, found by
//! a spec sweep, not assumed.
//!
//! # The chosen fully-orderable fixture
//!
//! **Chebyshev 0.5 dB, N = 3, f0 = 1.0 GHz, FBW = 0.70, Z0 = 50 Ω, 0402.** A
//! physically-sensible L-band wideband BPF on standard 0402 RF parts. Every one
//! of its synthesized L/C values resolves to a real curated JLCPCB Basic part
//! already in [`yee_filter::LCSC_PARTS`] (NO table extension was needed):
//!
//! ```text
//! shunt  L ≈ 3.6 nH  → C14033 (3.9 nH 0402, ~8%)   C ≈ 7.5 pF → C1576 (6.8 pF, ~10%)
//! series L ≈ 12 nH   → C24563 (12 nH 0402, exact)  C ≈ 2.0 pF → C1558 (2.0 pF, exact)
//! ```
//!
//! The symmetric ladder makes resonators 1 & 3 identical (Chebyshev g1 = g3), so
//! L1/L3 and C1/C3 group into shared part rows.
//!
//! # Caveat (ADR-0164)
//!
//! The realizable regime for a *complete* lumped board is the **wideband** lumped
//! BPF. A NARROW-band GHz BPF (the 10 % demo) has sub-floor series/shunt elements
//! and needs the **distributed/planar** track (microstrip/combline/hairpin — the
//! `dimension` module) or a **top-C-coupled** lumped topology (series coupling
//! caps instead of series LC resonators) to be fully orderable; that is a
//! documented follow-on, not a gap this gate hides.
//!
//! Pure-compute, deterministic, fast (no FDTD, no `#[ignore]`). **Do NOT weaken**:
//! the headline is ZERO blank LCSC #s — do not relax that to "≥ half" and do not
//! invent C-numbers. If a future table edit drops a part the fixture needs, this
//! gate must fail (that is the point).

use yee_filter::{
    Approximation, BOM_HEADER, Bom, CPL_HEADER, CompKind, ESeries, FilterSpec, Footprint,
    LCSC_PARTS, LcscPart, Response, SpecMask, autopick, autopick_within, jlcpcb_bom_csv,
    jlcpcb_cpl_csv, jlcpcb_files, join_placed_parts, lumped_board, select_components, synthesize,
    synthesize_lumped,
};
use yee_layout::Substrate;

/// The fully-orderable **wideband** lumped BPF: Chebyshev 0.5 dB, N = 3,
/// f0 = 1.0 GHz, FBW = 0.70, Z0 = 50 Ω. Every synthesized L/C lands in the
/// autopick table's realizable decades on **0402** (verified by the gate, found
/// by a spec sweep — see the module docs for why 0402 / this FBW).
fn orderable_fixture() -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: 1.0e9,
        fbw: 0.70,
        order: Some(3),
        z0_ohm: 50.0,
        // The mask is not the subject of this gate (orderability is); a wideband
        // BPF's stopband sits far out, so we record only an out-of-band point. The
        // realized-response mask is exercised by `lumped_001`; here the headline is
        // a COMPLETE orderable BOM, not mask compliance.
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 9.0,
            stopband: vec![(3.0e9, 30.0)],
        },
    }
}

/// The standard **narrow-band** demo BPF (the J1/J2/J3 fixture): 2 GHz, 10 % FBW.
/// Used ONLY for the one-line contrast assertion (it HAS blanks; ours has none).
fn narrow_fixture() -> FilterSpec {
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

/// Is this exact `LcscPart` present in the bundled table (real, not synthesized)?
fn in_table(part: &LcscPart) -> bool {
    LCSC_PARTS.iter().any(|p| p.lcsc == part.lcsc)
}

/// Human label for a value (pF / nH).
fn label(kind: CompKind, v: f64) -> String {
    match kind {
        CompKind::Capacitor => format!("{:.3} pF", v * 1e12),
        CompKind::Inductor => format!("{:.3} nH", v * 1e9),
    }
}

/// Split a CSV line into fields, honoring RFC-4180 double-quote escaping (a
/// quoted field may contain commas; `""` inside quotes is a literal `"`).
fn parse_csv_line(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '"' => {
                if in_quotes && chars.peek() == Some(&'"') {
                    cur.push('"');
                    chars.next();
                } else {
                    in_quotes = !in_quotes;
                }
            }
            ',' if !in_quotes => fields.push(std::mem::take(&mut cur)),
            other => cur.push(other),
        }
    }
    fields.push(cur);
    fields
}

#[test]
fn jlcpcb_orderable_001() {
    // 0402: the chosen footprint whose 2.7/3.9 nH parts catch the wideband shunt
    // inductor (see module docs — the binding constraint is the inductor floor).
    let footprint = Footprint::Smd0402;
    let series = ESeries::E24;
    let spec = orderable_fixture();
    let proj = synthesize(&spec);
    assert_eq!(proj.prototype.order(), 3, "fixture is order N=3");

    let ladder = synthesize_lumped(&proj).expect("wideband N=3 BPF must synthesize");

    // -------------------------------------------------------------------------
    // (0) Print the full spec + the synthesized ladder (the design under proof).
    // -------------------------------------------------------------------------
    println!("\n========================================================================");
    println!("jlcpcb-orderable-001 — fully-orderable JLCPCB board (J5, ADR-0164)");
    println!("========================================================================");
    println!(
        "SPEC: Chebyshev {:?} | f0 = {:.3} GHz | FBW = {:.0}% | N = {} | Z0 = {:.0} Ω | {:?}",
        spec.approximation,
        spec.f0_hz / 1e9,
        spec.fbw * 100.0,
        proj.prototype.order(),
        spec.z0_ohm,
        footprint,
    );
    println!("synthesized LC ladder (shunt-first, each tuned to ω0):");
    for (i, r) in ladder.resonators.iter().enumerate() {
        println!(
            "  res[{i}] {:?}: L = {:.3} nH, C = {:.3} pF",
            r.branch,
            r.l_henry * 1e9,
            r.c_farad * 1e12
        );
    }

    // -------------------------------------------------------------------------
    // (1) HEADLINE: the grouped BOM has ZERO blank LCSC #s — every distinct part
    //     resolves to a real, in-table, Basic LCSC C-number, value within
    //     tolerance, footprint matching.
    // -------------------------------------------------------------------------
    let bom: Bom = select_components(&ladder, series);
    assert!(!bom.lines.is_empty(), "BOM must have lines");

    // The autopicker's realistic tolerance band as a log-distance — used to
    // independently re-verify each pick is genuinely close to the synthesized
    // value (so a pick cannot be "right by construction").
    let max_log = (1.0 + yee_filter::DEFAULT_TOLERANCE_PCT / 100.0).log10();

    println!("\nBOM → LCSC (the orderable mapping):");
    let mut covered = 0usize;
    for line in &bom.lines {
        let pick = autopick(line, footprint).unwrap_or_else(|| {
            panic!(
                "ZERO-BLANK VIOLATION: {:?} {} did not resolve to any LCSC part \
                 (a fully-orderable spec must have NO blanks — see module docs)",
                line.kind,
                label(line.kind, line.chosen_value),
            )
        });
        covered += 1;
        println!(
            "  {:?} {} (chosen) → {} = {} (basic={}) ×{}",
            line.kind,
            label(line.kind, line.chosen_value),
            pick.lcsc,
            label(pick.kind, pick.value),
            pick.basic,
            line.qty,
        );

        // Correctness of every pick: kind, footprint, well-formed real C-number,
        // value genuinely close, Basic.
        assert_eq!(pick.kind, line.kind, "picked kind must match the line");
        assert_eq!(pick.footprint, footprint, "picked footprint must match");
        assert!(is_cnumber(pick.lcsc), "{} is not a C-number", pick.lcsc);
        assert!(in_table(&pick), "{} not present in LCSC_PARTS", pick.lcsc);
        let dist = (pick.value.log10() - line.chosen_value.log10()).abs();
        assert!(
            dist <= max_log + 1e-12,
            "{} value {} is {:.4} in log10 from chosen {} (band {:.4})",
            pick.lcsc,
            label(pick.kind, pick.value),
            dist,
            label(line.kind, line.chosen_value),
            max_log,
        );
        // The default and explicit-tolerance paths must agree (pick is stable).
        assert_eq!(
            autopick_within(line, footprint, yee_filter::DEFAULT_TOLERANCE_PCT),
            Some(pick),
            "default-tolerance path must be stable",
        );
        // Basic preference: a fully-orderable (free-assembly) board needs every
        // part to be a JLCPCB Basic part.
        assert!(
            pick.basic,
            "{} must be a JLCPCB Basic part (free assembly) for a fully-orderable board",
            pick.lcsc
        );
    }
    assert_eq!(
        covered,
        bom.lines.len(),
        "EVERY BOM line must resolve — zero blanks (got {covered}/{})",
        bom.lines.len()
    );
    println!(
        "coverage: {covered}/{} distinct BOM lines resolved — ZERO blanks (FULLY ORDERABLE)",
        bom.lines.len()
    );

    // -------------------------------------------------------------------------
    // (2) The full upload set: place the board, join, render BOM + CPL CSVs, and
    //     assert the BOM CSV has ZERO blank-LCSC rows (the orderable-board proof
    //     end-to-end, through the same J2/J3 export the studio/CLI ship).
    // -------------------------------------------------------------------------
    let board = lumped_board(&ladder, &fr4(), footprint);
    let placements = &board.placements;
    assert_eq!(placements.len(), 6, "N=3 → 2·N = 6 placed components");

    let parts = join_placed_parts(placements, &ladder, footprint, series);
    assert_eq!(
        parts.len(),
        placements.len(),
        "every placement must value-join (none skipped)"
    );
    // Every placed component resolved (the join's view of zero-blank).
    assert!(
        parts.iter().all(|p| p.lcsc.is_some()),
        "every placed component must autopick a real LCSC part (zero blanks on the board)"
    );

    let files = jlcpcb_files(placements, &ladder, footprint, series);
    let bom_csv = jlcpcb_bom_csv(&parts);
    assert_eq!(files.bom_csv, bom_csv, "jlcpcb_files BOM == direct BOM");
    assert_eq!(
        files.cpl_csv,
        jlcpcb_cpl_csv(placements),
        "jlcpcb_files CPL == direct CPL"
    );

    println!("\n=== JLCPCB BOM CSV (the complete orderable upload — no blank rows) ===");
    println!("{bom_csv}");
    println!("\n=== JLCPCB CPL CSV ===");
    println!("{}", files.cpl_csv);

    // ---- (2a) BOM CSV schema + ZERO blank rows ------------------------------
    let bom_lines: Vec<&str> = bom_csv.lines().collect();
    assert_eq!(bom_lines[0], BOM_HEADER, "BOM header exact");
    assert_eq!(
        bom_lines[0], "Comment,Designator,Footprint,LCSC Part #",
        "BOM header literal"
    );
    let bom_data = &bom_lines[1..];
    assert!(!bom_data.is_empty(), "BOM must have >= 1 data row");

    let mut bom_designators: Vec<String> = Vec::new();
    let mut realizable_rows = 0usize;
    for row in bom_data {
        let f = parse_csv_line(row);
        assert_eq!(f.len(), 4, "every BOM row has 4 fields: {row:?} → {f:?}");
        let (comment, designator, fp, lcsc) = (&f[0], &f[1], &f[2], &f[3]);
        assert_eq!(fp, "0402", "footprint column is the JLCPCB land name");
        assert!(!comment.is_empty(), "Comment must be non-empty");
        assert!(!designator.is_empty(), "Designator must be non-empty");

        // THE HEADLINE: no row may have a blank LCSC #, and none may be flagged
        // `(NO BASIC PART)`. A complete, orderable BOM.
        assert!(
            !lcsc.is_empty(),
            "BLANK LCSC # in a fully-orderable BOM row: {row:?} (zero blanks required)"
        );
        assert!(
            is_cnumber(lcsc),
            "LCSC # must be a well-formed C-number: {lcsc:?}"
        );
        assert!(
            !comment.contains("NO BASIC PART"),
            "no row may be flagged unfillable in a fully-orderable BOM: {comment:?}"
        );
        realizable_rows += 1;
        bom_designators.extend(designator.split(',').map(|s| s.to_string()));
    }
    assert!(
        realizable_rows >= 1,
        "the orderable BOM must have >= 1 realizable row"
    );
    println!(
        "BOM CSV: {realizable_rows} row(s), ALL with a real LCSC Basic part — ZERO blank rows"
    );

    // ---- (2b) CPL CSV schema: one row per placement -------------------------
    let cpl_lines: Vec<&str> = files.cpl_csv.lines().collect();
    assert_eq!(cpl_lines[0], CPL_HEADER, "CPL header exact");
    assert_eq!(
        cpl_lines[0], "Designator,Mid X,Mid Y,Layer,Rotation",
        "CPL header literal"
    );
    let cpl_data = &cpl_lines[1..];
    assert_eq!(
        cpl_data.len(),
        placements.len(),
        "exactly one CPL row per placement"
    );

    let bb = board.layout.bbox;
    let mut cpl_designators: Vec<String> = Vec::new();
    for row in cpl_data {
        let f = parse_csv_line(row);
        assert_eq!(f.len(), 5, "every CPL row has 5 fields: {row:?} → {f:?}");
        let (designator, mid_x, mid_y, layer, rotation) = (&f[0], &f[1], &f[2], &f[3], &f[4]);
        cpl_designators.push(designator.clone());
        assert_eq!(layer, "Top", "lumped board is single top layer");
        assert_eq!(rotation, "0", "axis-aligned placement → rotation 0");
        let x_mm = mid_x.strip_suffix("mm").expect("Mid X has mm suffix");
        let y_mm = mid_y.strip_suffix("mm").expect("Mid Y has mm suffix");
        let x_m: f64 = x_mm.parse::<f64>().expect("Mid X numeric") * 1e-3;
        let y_m: f64 = y_mm.parse::<f64>().expect("Mid Y numeric") * 1e-3;
        assert!(x_m.is_finite() && y_m.is_finite(), "coords finite");
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

    // ---- (2c) Designator consistency: BOM ref-des set == CPL set == placements
    let sorted_unique = |mut v: Vec<String>| -> Vec<String> {
        v.sort();
        v.dedup();
        v
    };
    let bom_set = sorted_unique(bom_designators.clone());
    let cpl_set = sorted_unique(cpl_designators.clone());
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
    let placement_set = sorted_unique(placements.iter().map(|p| p.ref_des.clone()).collect());
    assert_eq!(
        cpl_set, placement_set,
        "CPL designators must equal the placement ref-des set"
    );

    println!(
        "consistency: {} placements; BOM designators == CPL designators == placement set",
        placements.len()
    );

    // -------------------------------------------------------------------------
    // (3) HONEST BOUND — the one-line contrast: the standard NARROW-band 10 %
    //     fixture (J1/J2/J3) HAS blank LCSC rows, while this wideband spec has
    //     NONE. This bounds the realizable regime without faking either side.
    // -------------------------------------------------------------------------
    let narrow = narrow_fixture();
    let narrow_proj = synthesize(&narrow);
    let narrow_ladder = synthesize_lumped(&narrow_proj).expect("narrow N=3 BPF synthesizes");
    let narrow_bom = select_components(&narrow_ladder, series);
    let narrow_blanks = narrow_bom
        .lines
        .iter()
        .filter(|l| autopick(l, footprint).is_none())
        .count();
    assert!(
        narrow_blanks >= 1,
        "CONTRAST: the narrow-band 10% fixture must have >= 1 unrealizable (blank) line \
         (sub-floor series cap / shunt L) — bounding the realizable regime honestly"
    );
    println!(
        "\nCONTRAST (honest regime bound): narrow-band 2 GHz / 10% FBW BPF has {narrow_blanks} \
         blank (unrealizable) BOM line(s) on {footprint:?}; this wideband 1 GHz / 70% BPF has 0. \
         Narrow-band GHz BPFs need the distributed/planar track or a top-C-coupled topology \
         (ADR-0164 caveat)."
    );
    println!("========================================================================\n");
}
