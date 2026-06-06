//! top-c-board-001 (JLCPCB production track, ADR-0166 brick **T2** — the top-C
//! analogue of [`jlcpcb-orderable-001`](../../tests/jlcpcb_orderable_001.rs)): a
//! **top-C-coupled** BPF spec in the sub-GHz / moderate-band regime the ADR-0165
//! T1 envelope identified as orderable yields a **FULLY-ORDERABLE** JLCPCB upload
//! set — every shunt L, every shunt C, **and** every series coupling cap
//! autopicks to a real LCSC **Basic** C-number (**ZERO** blank rows), with a
//! schema-valid BOM whose designators match the CPL, and every placement inside
//! the board outline.
//!
//! # What this proves (and the regime it closes)
//!
//! [`jlcpcb-orderable-001`](../../tests/jlcpcb_orderable_001.rs) showed the
//! **alternating** `LumpedLadder` is fully orderable only for a **WIDEBAND** BPF
//! (its *series* LC resonators go sub-pF/sub-nH at narrow band). The
//! **top-C-coupled** topology ([`yee_filter::top_c`]) is the complementary fix:
//! `N` freely-realizable **shunt** parallel-LC resonators coupled by `N+1`
//! **series coupling caps** — no series LC resonator to collapse. So it extends
//! the orderable regime into the **sub-GHz / moderate-band** corner the
//! alternating ladder blanks in. This gate pins that: the
//! `synthesize_top_c_coupled` → [`top_c_board`] → [`join_top_c_parts`] →
//! [`jlcpcb_bom_csv`] / [`jlcpcb_cpl_csv`] data path produces a **complete,
//! orderable board with NO blanks** for the T1-envelope cell.
//!
//! # The chosen fully-orderable fixture (the T1 envelope cell)
//!
//! **Chebyshev 0.5 dB, N = 3, f0 = 0.5 GHz, FBW = 20 %, Z0 = 50 Ω, 0402.** The
//! sub-GHz/moderate-band cell ADR-0165 T1 mapped as orderable. Its synthesized
//! network (probed; verified by this gate) and the LCSC parts every value
//! resolves to on **0402**:
//!
//! ```text
//! shunt  L = 15.92 nH (×3) → C27143 (15 nH 0402)
//! shunt  C1=C3 = 3.30 pF    → C1565 (3.3 pF 0402)   C2 = 4.44 pF → C1569 (4.7 pF 0402)
//! coupl  Cc1=Cc4 = 2.41 pF  → C1559 (2.2 pF 0402)   Cc2=Cc3 = 0.96 pF → C1550 (1.0 pF 0402)
//! ```
//!
//! The **binding constraint is the interior coupling cap** `Cc2`/`Cc3` ≈ 0.96 pF:
//! it lands right on the **1 pF** Basic floor (`C1550`), within the autopick band.
//! This is exactly the value family that makes a *narrow*-band top-C (sub-pF
//! coupling caps) blank — at this moderate 20 % FBW it just clears the floor,
//! which is why the T1 envelope stops here and a narrower band needs the
//! distributed track (ADR-0166 caveat). **0402** is required: the 0603 RF
//! inductor grid skips from 12 nH to 22 nH, so the 15.92 nH shunt L blanks there
//! (verified by a sweep — documented, not assumed); this matches why
//! `jlcpcb-orderable-001` also uses 0402.
//!
//! The symmetric Chebyshev prototype (g1 = g3) makes resonators 1 & 3 identical,
//! so L1/L3, C1/C3, and Cc1/Cc4 (and Cc2/Cc3) group into shared part rows.
//!
//! Pure-compute, deterministic, fast (no FDTD, no `#[ignore]`). **Do NOT weaken**:
//! the headline is ZERO blank LCSC #s across all three arms — do not relax it and
//! do not invent C-numbers. If a value does not resolve, this gate must FAIL
//! loudly (that is the honest T1-envelope-vs-realized-board check the ADR calls
//! for).

use yee_filter::{
    Approximation, BOM_HEADER, CPL_HEADER, CompKind, Footprint, LCSC_PARTS, LcscPart,
    jlcpcb_bom_csv, jlcpcb_cpl_csv, join_top_c_parts, synthesize_top_c_coupled, top_c_board,
};
use yee_layout::Substrate;

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
fn top_c_board_001() {
    // 0402: the footprint whose 1 pF cap floor catches the moderate-band coupling
    // caps and whose 15 nH part catches the shunt inductor (see module docs — the
    // 0603 grid blanks the 16 nH shunt L).
    let footprint = Footprint::Smd0402;

    // ---- Synthesize the T1-envelope orderable cell --------------------------
    let net = synthesize_top_c_coupled(
        Approximation::Chebyshev { ripple_db: 0.5 },
        3,
        0.5e9,
        0.20,
        50.0,
    );
    assert_eq!(net.shunt.len(), 3, "N=3 → 3 shunt resonators");
    assert_eq!(
        net.coupling_caps_farad.len(),
        4,
        "N=3 → N+1 = 4 series coupling caps"
    );

    // -------------------------------------------------------------------------
    // (0) Print the spec + the synthesized network (the design under proof).
    // -------------------------------------------------------------------------
    println!("\n========================================================================");
    println!("top-c-board-001 — fully-orderable top-C-coupled JLCPCB board (T2, ADR-0166)");
    println!("========================================================================");
    println!(
        "SPEC: Chebyshev 0.5 dB | f0 = {:.3} GHz | FBW = {:.0}% | N = {} | Z0 = {:.0} Ω | {:?}",
        net.f0_hz / 1e9,
        net.fbw * 100.0,
        net.shunt.len(),
        net.z0_ohm,
        footprint,
    );
    println!(
        "synthesized top-C network ({} shunt resonators + {} series coupling caps):",
        net.shunt.len(),
        net.coupling_caps_farad.len()
    );
    for (i, r) in net.shunt.iter().enumerate() {
        println!(
            "  shunt[{i}]: L = {:.3} nH  C = {:.3} pF",
            r.l_henry * 1e9,
            r.c_farad * 1e12
        );
    }
    for (j, c) in net.coupling_caps_farad.iter().enumerate() {
        println!("  Cc{}: {:.3} pF", j + 1, c * 1e12);
    }

    // -------------------------------------------------------------------------
    // (1) Place the board: through-arm of 2N+1 = 7 slots → 3N+1 = 10 footprints
    //     (N+1 coupling caps + N·(L+C) resonators).
    // -------------------------------------------------------------------------
    let board = top_c_board(&net, &fr4(), footprint);
    let placements = &board.placements;
    assert_eq!(
        placements.len(),
        3 * net.shunt.len() + 1,
        "N=3 → 3·N+1 = 10 placed components (4 Cc + 3·(L+C))"
    );

    // Every footprint is on the board's chosen footprint (sanity).
    for p in placements {
        assert_eq!(p.footprint, footprint, "placement {} footprint", p.ref_des);
    }

    // -------------------------------------------------------------------------
    // (2) HEADLINE: join → every placement resolves to a real, in-table, Basic
    //     LCSC part within tolerance — ZERO blanks across ALL THREE arms (shunt
    //     L, shunt C, coupling cap).
    // -------------------------------------------------------------------------
    let parts = join_top_c_parts(placements, &net);
    assert_eq!(
        parts.len(),
        placements.len(),
        "every placement must value-join (none skipped)"
    );

    // The autopicker's realistic tolerance band as a log-distance — used to
    // independently re-verify each pick is genuinely close to the chosen value
    // (so a pick cannot be "right by construction").
    let max_log = (1.0 + yee_filter::DEFAULT_TOLERANCE_PCT / 100.0).log10();

    println!("\nplacement → LCSC (the orderable mapping; all three arms):");
    let (mut n_shunt_l, mut n_shunt_c, mut n_coupling) = (0usize, 0usize, 0usize);
    for p in &parts {
        let pick = p.lcsc.unwrap_or_else(|| {
            panic!(
                "ZERO-BLANK VIOLATION: {} ({:?} {}) did not resolve to any LCSC part \
                 — a fully-orderable top-C board must have NO blanks across shunt L / \
                 shunt C / coupling cap (see module docs / ADR-0166)",
                p.ref_des,
                p.kind,
                label(p.kind, p.chosen_value),
            )
        });
        // Tally which arm each ref-des belongs to (Cc* = coupling, L* = shunt L,
        // C* = shunt C) so we prove all three arms are covered, not just two.
        if p.ref_des.starts_with("Cc") {
            n_coupling += 1;
        } else if p.ref_des.starts_with('L') {
            n_shunt_l += 1;
        } else if p.ref_des.starts_with('C') {
            n_shunt_c += 1;
        } else {
            panic!("unexpected top-C ref-des {}", p.ref_des);
        }
        println!(
            "  {:<5} {:?} {} (chosen) → {} = {} (basic={})",
            p.ref_des,
            p.kind,
            label(p.kind, p.chosen_value),
            pick.lcsc,
            label(pick.kind, pick.value),
            pick.basic,
        );

        // Correctness of every pick: kind, footprint, well-formed real C-number,
        // value genuinely close, Basic.
        assert_eq!(pick.kind, p.kind, "{} picked kind must match", p.ref_des);
        assert_eq!(pick.footprint, footprint, "{} picked footprint", p.ref_des);
        assert!(is_cnumber(pick.lcsc), "{} is not a C-number", pick.lcsc);
        assert!(in_table(&pick), "{} not present in LCSC_PARTS", pick.lcsc);
        let dist = (pick.value.log10() - p.chosen_value.log10()).abs();
        assert!(
            dist <= max_log + 1e-12,
            "{} value {} is {:.4} in log10 from chosen {} (band {:.4})",
            pick.lcsc,
            label(pick.kind, pick.value),
            dist,
            label(p.kind, p.chosen_value),
            max_log,
        );
        assert!(
            pick.basic,
            "{} must be a JLCPCB Basic part (free assembly) for a fully-orderable board",
            pick.lcsc
        );
    }
    // All three arms must be present and covered (N shunt L, N shunt C, N+1 Cc).
    assert_eq!(
        n_shunt_l,
        net.shunt.len(),
        "every shunt inductor placed+covered"
    );
    assert_eq!(
        n_shunt_c,
        net.shunt.len(),
        "every shunt capacitor placed+covered"
    );
    assert_eq!(
        n_coupling,
        net.coupling_caps_farad.len(),
        "every coupling capacitor placed+covered"
    );
    println!(
        "coverage: {}/{} placements resolved — ZERO blanks across all 3 arms \
         ({n_shunt_l} shunt L + {n_shunt_c} shunt C + {n_coupling} coupling caps) — FULLY ORDERABLE",
        parts.len(),
        placements.len()
    );

    // -------------------------------------------------------------------------
    // (3) The full upload set: render the BOM + CPL through the SAME J2/J3
    //     emitters the lumped path uses, and assert ZERO blank-LCSC rows.
    // -------------------------------------------------------------------------
    let bom_csv = jlcpcb_bom_csv(&parts);
    let cpl_csv = jlcpcb_cpl_csv(placements);

    println!("\n=== JLCPCB BOM CSV (the complete orderable upload — no blank rows) ===");
    println!("{bom_csv}");
    println!("\n=== JLCPCB CPL CSV ===");
    println!("{cpl_csv}");

    // ---- (3a) BOM CSV schema + ZERO blank rows ------------------------------
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
        // `(NO BASIC PART)`. A complete, orderable top-C board.
        assert!(
            !lcsc.is_empty(),
            "BLANK LCSC # in a fully-orderable top-C BOM row: {row:?} (zero blanks required)"
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

    // ---- (3b) CPL CSV schema: one row per placement, inside the board outline -
    let cpl_lines: Vec<&str> = cpl_csv.lines().collect();
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
        // Every placement centre lies within the board outline (the bbox over all
        // copper) — the ADR (c) check.
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

    // ---- (3c) Designator consistency: CPL designators EXACTLY match the BOM ---
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
        "CPL designators must EXACTLY match the BOM designators (every placed part in both)"
    );
    let placement_set = sorted_unique(placements.iter().map(|p| p.ref_des.clone()).collect());
    assert_eq!(
        cpl_set, placement_set,
        "CPL designators must equal the placement ref-des set"
    );

    println!(
        "consistency: {} placements; BOM designators == CPL designators == placement set; \
         all centres inside the board outline",
        placements.len()
    );
    println!("========================================================================\n");
}
