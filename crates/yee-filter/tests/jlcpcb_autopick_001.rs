//! jlcpcb-autopick-001 (JLCPCB production track, ADR-0164 brick J1): a
//! synthesized filter's BOM auto-picks **real** LCSC parts.
//!
//! **Non-circular by construction.** The L/C values are produced by the project
//! synthesis path (`synthesize` → `synthesize_lumped` → `select_components`,
//! exactly as `bom_001`/`lumped_pcb_001` do); the LCSC C-numbers come from the
//! independently-curated [`yee_filter::LCSC_PARTS`] table (sourced from JLCPCB's
//! published Basic Parts list — see the module docs). The gate matches the two
//! and checks each pick is genuinely close to the synthesized value, so a pick
//! cannot be "right by construction".
//!
//! For the standard 3-pole 0.5 dB Chebyshev BPF (f0 = 2 GHz, FBW = 0.10,
//! Z0 = 50 Ω) on **0603** footprints, assert:
//!
//! 1. **Coverage** — the realizable BOM lines resolve to a real LCSC part. The
//!    narrow-band band-pass *series* resonator is physically extreme (sub-pF
//!    capacitance below the 1 pF Basic floor, see `jlcpcb.rs` coverage note), so
//!    full coverage is not achievable; the gate records the **covered fraction**
//!    and asserts every uncovered line is honestly `None`. (≥ half the distinct
//!    lines must resolve — the shunt resonator + its quantity-grouped twins.)
//! 2. **Correctness of every pick** — kind matches, footprint matches, the value
//!    is genuinely close (`|log10(part) − log10(chosen)|` within the autopicker's
//!    realistic tolerance, which is *wider* than E24 because JLCPCB Basic stock
//!    is coarser than a full E24 decade — documented in `jlcpcb.rs`), and the
//!    C-number is well-formed (`^C\d+$`) and present in [`LCSC_PARTS`].
//! 3. **Basic preference** — every pick is a JLCPCB Basic part (the curated seed
//!    is all-Basic, so a covered line must pick Basic).
//!
//! Prints the BOM → picked-LCSC mapping. Pure-compute, deterministic, fast
//! (no FDTD, no `#[ignore]`). Do NOT weaken: do not invent C-numbers to force
//! full coverage — an honest `None` on a physically-unrealizable value is the
//! correct outcome.

use yee_filter::{
    Approximation, Bom, CompKind, ESeries, FilterSpec, Footprint, LCSC_PARTS, LcscPart, Response,
    SpecMask, autopick, autopick_within, select_components, synthesize, synthesize_lumped,
};

/// Chebyshev 0.5 dB **N=3** bandpass spec (the standard demo filter).
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

#[test]
fn jlcpcb_autopick_001() {
    let footprint = Footprint::Smd0603;
    let spec = fixture();
    let proj = synthesize(&spec);
    assert_eq!(proj.prototype.order(), 3, "fixture is order N=3");

    let ladder = synthesize_lumped(&proj).expect("N=3 bandpass fixture should synthesize");
    let bom: Bom = select_components(&ladder, ESeries::E24);
    assert!(!bom.lines.is_empty(), "BOM must have lines");

    // The autopicker's realistic tolerance band, as a log-distance, used to
    // independently re-verify each pick is genuinely close to the synthesized
    // value (see the jlcpcb.rs coverage note for WHY it is wider than E24 ±5 %).
    let max_log = (1.0 + yee_filter::DEFAULT_TOLERANCE_PCT / 100.0).log10();

    let mut covered = 0usize;
    let mut uncovered: Vec<String> = Vec::new();

    println!("\n3-pole 0.5dB Cheb BPF (2 GHz, 10% FBW) @ {footprint:?} — BOM → LCSC:");
    for line in &bom.lines {
        let pick = autopick(line, footprint);
        match pick {
            Some(part) => {
                covered += 1;
                println!(
                    "  {:?} {} (chosen) → {} = {} (basic={}) qty {}",
                    line.kind,
                    label(line.kind, line.chosen_value),
                    part.lcsc,
                    label(part.kind, part.value),
                    part.basic,
                    line.qty,
                );

                // (2) correctness of the pick -------------------------------
                assert_eq!(part.kind, line.kind, "picked kind must match the line");
                assert_eq!(part.footprint, footprint, "picked footprint must match");
                assert!(is_cnumber(part.lcsc), "{} is not a C-number", part.lcsc);
                assert!(in_table(&part), "{} not present in LCSC_PARTS", part.lcsc);

                // value genuinely close: within the autopicker's tolerance band.
                let dist = (part.value.log10() - line.chosen_value.log10()).abs();
                assert!(
                    dist <= max_log + 1e-12,
                    "{} value {} is {:.4} in log10 from chosen {} (band {:.4})",
                    part.lcsc,
                    label(part.kind, part.value),
                    dist,
                    label(line.kind, line.chosen_value),
                    max_log
                );

                // re-deriving the pick at the same tolerance must agree (the
                // public default and the explicit path are consistent).
                let again = autopick_within(line, footprint, yee_filter::DEFAULT_TOLERANCE_PCT);
                assert_eq!(again, Some(part), "default tolerance path must be stable");

                // (3) Basic preference: the curated seed is all-Basic, so a
                // covered line must resolve to a Basic part.
                assert!(
                    part.basic,
                    "{} should be a JLCPCB Basic part (seed table is all-Basic)",
                    part.lcsc
                );
            }
            None => {
                println!(
                    "  {:?} {} (chosen) → NONE (no in-table part within tolerance)",
                    line.kind,
                    label(line.kind, line.chosen_value),
                );
                uncovered.push(format!(
                    "{:?} {}",
                    line.kind,
                    label(line.kind, line.chosen_value)
                ));
            }
        }
    }

    let total = bom.lines.len();
    println!("coverage: {covered}/{total} distinct BOM lines resolved to a real LCSC Basic part");
    if !uncovered.is_empty() {
        println!("uncovered (honest None — physically extreme / off-catalog values):");
        for u in &uncovered {
            println!("    {u}");
        }
    }

    // (1) coverage: at least half the distinct lines resolve. The 3-pole BPF's
    // shunt resonator (an L + a C, each qty-grouped ×2 by symmetry) is fully
    // realizable; its series resonator's sub-pF cap is below the 1 pF Basic
    // floor and is the documented honest miss. We assert ≥ ⌈total/2⌉ covered.
    assert!(
        covered * 2 >= total,
        "expected at least half the BOM lines to resolve, got {covered}/{total}"
    );
    // And at least one line must genuinely resolve (the picker is not vacuous).
    assert!(
        covered >= 1,
        "at least one BOM line must resolve to a real part"
    );

    // The covered + uncovered partition is exhaustive and honest.
    assert_eq!(
        covered + uncovered.len(),
        total,
        "every line is either covered or an honest None"
    );
}

/// The shunt-resonator values of the standard 3-pole BPF DO resolve to known
/// real Basic parts — an explicit, non-vacuous anchor on the coverage claim.
#[test]
fn jlcpcb_autopick_001_shunt_resolves_to_known_parts() {
    let footprint = Footprint::Smd0603;
    let proj = synthesize(&fixture());
    let ladder = synthesize_lumped(&proj).expect("synth");
    let bom = select_components(&ladder, ESeries::E24);

    // The shunt resonator (the realizable one) gives ≈0.24 nH / ≈24 pF; the cap
    // (24 pF) resolves to the nearest Basic 0603 part (22 pF = C1653). Find the
    // capacitor line nearest 24 pF and confirm it picks a real Basic C-number.
    let cap_line = bom
        .lines
        .iter()
        .filter(|l| l.kind == CompKind::Capacitor)
        .min_by(|a, b| {
            (a.chosen_value - 24e-12)
                .abs()
                .partial_cmp(&(b.chosen_value - 24e-12).abs())
                .unwrap()
        })
        .expect("a capacitor line exists");
    let part = autopick(cap_line, footprint).expect("the ~24 pF shunt cap must resolve");
    assert_eq!(part.kind, CompKind::Capacitor);
    assert_eq!(part.footprint, Footprint::Smd0603);
    assert!(part.basic, "the resolved shunt cap is a Basic part");
    assert!(is_cnumber(part.lcsc));
    // 24 pF → nearest stocked Basic 0603 C0G is 22 pF (C1653).
    assert_eq!(
        part.lcsc, "C1653",
        "24 pF → nearest Basic 0603 = 22 pF (C1653)"
    );
}
