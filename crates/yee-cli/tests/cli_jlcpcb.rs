//! Gate for `yee filter synth --jlcpcb <dir>` (JLCPCB production track, ADR-0164
//! brick **J4**).
//!
//! Runs the CLI against the committed satisfiable Chebyshev N=5 bandpass fixture
//! with `--jlcpcb <tmpdir>` and asserts the directory holds the JLCPCB
//! SMT-assembly upload set:
//!
//! - `bom.csv` — exact JLCPCB header `Comment,Designator,Footprint,LCSC Part #`,
//!   at least one **realizable** row carrying a well-formed `C\d+` LCSC #, and
//!   the **unrealizable** rows present with a blank LCSC # (the honest
//!   coverage-hole flag — values are flagged, never dropped).
//! - `cpl.csv` — exact JLCPCB header `Designator,Mid X,Mid Y,Layer,Rotation`,
//!   one row per placed component, with the CPL designator set consistent with
//!   the BOM designator union (both derive from the same placement list).
//! - a structurally-valid RS-274X copper Gerber (`%FS` / `G36` / `G37` / `M02`).
//!
//! Pure closed-form geometry (no EM/FDTD), so NOT `#[ignore]`'d — mirrors
//! `cli_lumped_export` / `cli_finite_q_s2p`.

use std::collections::BTreeSet;
use std::path::PathBuf;
use std::process::Command;

/// Absolute path to the committed Chebyshev N=5 bandpass fixture (shared with
/// the planar / lumped Gerber gates).
fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/cheb_bpf.toml")
}

/// Split a CSV data line into its (RFC-4180) fields, unwrapping a single layer
/// of quoting. The only quoted field the JLCPCB CSVs emit is a comma-joined
/// `Designator` cell (e.g. `"C1,C5"`) with no inner quotes, so a one-pass
/// quote-aware split is sufficient here.
fn csv_fields(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    for c in line.chars() {
        match c {
            '"' => in_quotes = !in_quotes,
            ',' if !in_quotes => {
                fields.push(std::mem::take(&mut cur));
            }
            _ => cur.push(c),
        }
    }
    fields.push(cur);
    fields
}

/// True iff `s` is a well-formed LCSC part number: a `C` followed by ≥ 1 digit.
fn is_lcsc_part(s: &str) -> bool {
    let mut chars = s.chars();
    chars.next() == Some('C') && {
        let rest = chars.as_str();
        !rest.is_empty() && rest.chars().all(|c| c.is_ascii_digit())
    }
}

/// Collect the designator set from a comma-joined `Designator` CSV cell.
fn designators(cell: &str) -> impl Iterator<Item = String> + '_ {
    cell.split(',').map(|s| s.trim().to_string())
}

#[test]
fn cli_jlcpcb_writes_upload_set() {
    // A fresh per-test directory under the cargo target tmpdir.
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli_jlcpcb_set");
    let _ = std::fs::remove_dir_all(&dir);
    // A throwaway .s2p next to the dir so the synth run does not litter the
    // fixture directory.
    let s2p = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli_jlcpcb.s2p");
    let _ = std::fs::remove_file(&s2p);

    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["filter", "synth"])
        .arg(fixture())
        .arg("--output")
        .arg(&s2p)
        .arg("--jlcpcb")
        .arg(&dir)
        .output()
        .expect("invoke yee");
    assert!(
        output.status.success(),
        "yee filter synth --jlcpcb exited non-zero; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let bom_path = dir.join("bom.csv");
    let cpl_path = dir.join("cpl.csv");
    assert!(
        bom_path.exists(),
        "bom.csv was not written to {}",
        dir.display()
    );
    assert!(
        cpl_path.exists(),
        "cpl.csv was not written to {}",
        dir.display()
    );

    // ---- BOM CSV ---------------------------------------------------------
    let bom = std::fs::read_to_string(&bom_path).expect("read bom.csv");
    let bom_lines: Vec<&str> = bom.lines().collect();
    assert_eq!(
        bom_lines[0], "Comment,Designator,Footprint,LCSC Part #",
        "bom.csv header must be the exact JLCPCB assembly-BOM columns"
    );
    assert!(
        bom_lines.len() > 1,
        "bom.csv has only a header — no part rows"
    );

    // Each data row has exactly four fields; collect realizable (well-formed
    // C-number) vs unrealizable (blank LCSC #) rows + the BOM designator union.
    let mut realizable_lcsc: Vec<String> = Vec::new();
    let mut saw_unrealizable_blank = false;
    let mut bom_designators: BTreeSet<String> = BTreeSet::new();
    for line in &bom_lines[1..] {
        let f = csv_fields(line);
        assert_eq!(
            f.len(),
            4,
            "bom.csv row must have 4 fields (Comment,Designator,Footprint,LCSC #); got {f:?}"
        );
        for d in designators(&f[1]) {
            assert!(!d.is_empty(), "empty designator in bom.csv row: {line}");
            bom_designators.insert(d);
        }
        let lcsc = &f[3];
        if lcsc.is_empty() {
            // An unrealizable row: blank LCSC #, flagged in the Comment.
            saw_unrealizable_blank = true;
            assert!(
                f[0].contains("NO BASIC PART"),
                "blank-LCSC row must be flagged in its Comment; got {f:?}"
            );
        } else {
            assert!(
                is_lcsc_part(lcsc),
                "realizable row LCSC # {lcsc:?} is not a well-formed C\\d+ part"
            );
            realizable_lcsc.push(lcsc.clone());
        }
    }
    assert!(
        !realizable_lcsc.is_empty(),
        "bom.csv has no realizable row with a C\\d+ LCSC # — autopick resolved nothing"
    );
    // The N=5 Chebyshev fixture has sub-pF series-resonator caps / tiny nH
    // inductors with no Basic match, so the unrealizable (blank-LCSC) rows MUST
    // be present and not silently dropped.
    assert!(
        saw_unrealizable_blank,
        "bom.csv dropped the unrealizable (no-Basic-match) rows — they must be \
         emitted with a blank LCSC #, not removed"
    );

    // ---- CPL CSV ---------------------------------------------------------
    let cpl = std::fs::read_to_string(&cpl_path).expect("read cpl.csv");
    let cpl_lines: Vec<&str> = cpl.lines().collect();
    assert_eq!(
        cpl_lines[0], "Designator,Mid X,Mid Y,Layer,Rotation",
        "cpl.csv header must be the exact JLCPCB centroid columns"
    );
    assert!(
        cpl_lines.len() > 1,
        "cpl.csv has only a header — no placement rows"
    );

    let mut cpl_designators: BTreeSet<String> = BTreeSet::new();
    for line in &cpl_lines[1..] {
        let f = csv_fields(line);
        assert_eq!(
            f.len(),
            5,
            "cpl.csv row must have 5 fields (Designator,Mid X,Mid Y,Layer,Rotation); got {f:?}"
        );
        assert!(!f[0].is_empty(), "empty designator in cpl.csv row: {line}");
        cpl_designators.insert(f[0].clone());
        // Millimetre coordinates with the `mm` suffix.
        assert!(
            f[1].ends_with("mm") && f[2].ends_with("mm"),
            "cpl.csv Mid X/Y must carry the `mm` unit suffix; got {f:?}"
        );
        assert_eq!(f[3], "Top", "lumped board is single-layer top copper");
    }

    // The CPL designators must be consistent with the BOM designator union:
    // every physical placement in the CPL is accounted for in the BOM (the BOM
    // emits even unrealizable parts), and the BOM lists no designator absent
    // from the board. For a `lumped_board` the two sets are exactly equal.
    assert_eq!(
        cpl_designators, bom_designators,
        "CPL and BOM designator sets disagree — they must both derive from the \
         same placement list (CPL {cpl_designators:?} vs BOM {bom_designators:?})"
    );

    // ---- Gerber ----------------------------------------------------------
    // A structurally-valid RS-274X copper Gerber must exist in the dir. The
    // copper file is named for the spec stem (`cheb_bpf.gbr`); assert that one
    // specifically, then validate its RS-274X structure (mirrors the
    // `cli_lumped_export` Gerber-validity check).
    let gbr_path = dir.join("cheb_bpf.gbr");
    assert!(
        gbr_path.exists(),
        "copper Gerber {} was not written",
        gbr_path.display()
    );
    let gerber = std::fs::read_to_string(&gbr_path).expect("read copper Gerber");
    assert!(
        gerber.contains("%FSLAX46Y46*%"),
        "Gerber missing %FS coordinate-format header"
    );
    assert!(gerber.contains("G36*"), "Gerber missing G36 region-open");
    assert!(gerber.contains("G37*"), "Gerber missing G37 region-close");
    assert!(
        gerber.contains("M02*"),
        "Gerber missing M02 end-of-file marker"
    );

    // Clean up.
    let _ = std::fs::remove_dir_all(&dir);
    let _ = std::fs::remove_file(&s2p);
}
