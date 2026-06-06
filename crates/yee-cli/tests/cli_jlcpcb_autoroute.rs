//! Gate `cli-jlcpcb-autoroute` for `yee filter synth --jlcpcb <dir>` topology
//! auto-routing (JLCPCB production track, ADR-0168 brick **T4**).
//!
//! T4 wires the [`yee_filter::synthesize_orderable_on`] selector into the CLI:
//! `--jlcpcb` no longer hardcodes the alternating ladder — it auto-routes
//! whichever lumped topology (alternating series/shunt ladder vs top-C-coupled)
//! yields an orderable JLCPCB board for the spec on the user's substrate, reports
//! which one was chosen, and honestly flags when neither is fully orderable.
//!
//! # What this proves (and why it is non-circular)
//!
//! The gate drives the **real CLI binary** (`Command::new(CARGO_BIN_EXE_yee)`,
//! the same invocation the J4 `cli_jlcpcb` gate uses) for three specs that route
//! differently, then reads the **actually-written** `bom.csv` / `cpl.csv` and the
//! run's **reported topology on stdout** — never the in-process selector return.
//! The topology is checked two independent ways that must agree: the stdout
//! `topology:` line *and* the file-based marker (a top-C board carries `Cc*`
//! coupling-cap designators in the BOM/CPL; the alternating ladder does not).
//!
//! # The three discriminating specs (all Chebyshev 0.5 dB, N = 3, Z0 = 50 Ω, 0402)
//!
//! Each run passes `--footprint 0402` explicitly (the CLI default is 0603); the
//! 0402 routing is what these cases assert.
//!
//! 1. **0.5 GHz / 20 %** — the alternating ladder blanks here but top-C is fully
//!    orderable, so the run routes to **top-C** with a **zero-blank** `bom.csv`
//!    (every data row carries a `C\d+` LCSC #). The headline T4 win: a sub-GHz
//!    spec the old fixed-ladder path could not make orderable now emits a
//!    fully-orderable upload set automatically.
//! 2. **1.0 GHz / 70 %** (wideband) — the ladder is fully orderable, so the run
//!    routes to the **alternating ladder** with a **zero-blank** `bom.csv` (and no
//!    `Cc*` coupling caps). The conventional topology is kept where it works.
//! 3. **2.0 GHz / 5 %** (GHz-narrow) — NEITHER lumped topology is fully orderable,
//!    so the run reports the honest "neither lumped topology is fully orderable"
//!    note and `bom.csv` carries ≥ 1 real blank LCSC row (the coverage hole is
//!    flagged, never dropped).
//!
//! Pure closed-form geometry (no EM/FDTD), so NOT `#[ignore]`'d — mirrors
//! `cli_jlcpcb` / `cli_lumped_export`. **Do NOT weaken**: case 1's zero-blank
//! top-C BOM is the load-bearing proof the auto-route rescues a spec the ladder
//! blanks; if a future synthesis/table edit changes the routing, this gate must
//! FAIL (that is the point) — do not relax it and do not invent C-numbers.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Split a CSV data line into its (RFC-4180) fields, unwrapping a single layer of
/// quoting. The only quoted field the JLCPCB CSVs emit is a comma-joined
/// `Designator` cell (e.g. `"C1,C5"`) with no inner quotes, so a one-pass
/// quote-aware split is sufficient here.
fn csv_fields(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    for c in line.chars() {
        match c {
            '"' => in_quotes = !in_quotes,
            ',' if !in_quotes => fields.push(std::mem::take(&mut cur)),
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

/// The designators in a comma-joined `Designator` CSV cell.
fn designators(cell: &str) -> impl Iterator<Item = String> + '_ {
    cell.split(',').map(|s| s.trim().to_string())
}

/// A Chebyshev 0.5 dB / N=3 / Z0=50 Ω band-pass spec TOML at the given centre and
/// fractional bandwidth (the `yee filter synth` `FilterSpec` schema). A single
/// stopband row at `1.5·f0` keeps the spec well-formed; the mask verdict does not
/// gate this test (it asserts on the written JLCPCB files + reported topology).
fn spec_toml(f0_hz: f64, fbw: f64) -> String {
    let stop = f0_hz * 1.5;
    format!(
        "response = \"Bandpass\"\n\
         f0_hz = {f0_hz:e}\n\
         fbw = {fbw}\n\
         order = 3\n\
         z0_ohm = 50.0\n\
         \n\
         [approximation.Chebyshev]\n\
         ripple_db = 0.5\n\
         \n\
         [mask]\n\
         passband_ripple_db = 0.5\n\
         return_loss_db = 9.0\n\
         stopband = [[{stop:e}, 30.0]]\n"
    )
}

/// The outcome of a single `--jlcpcb` run, read back from the written artifacts +
/// the run's stdout (non-circular: no in-process selector return is consulted).
struct Outcome {
    /// The run's stdout (carries the `topology: … (auto-selected)` line + the
    /// orderability report).
    stdout: String,
    /// The BOM data rows (header stripped), each split into its 4 CSV fields.
    bom_rows: Vec<Vec<String>>,
    /// The CPL designator set (every placed component).
    cpl_designators: BTreeSet<String>,
    /// The BOM designator union (across all rows, realizable + blank).
    bom_designators: BTreeSet<String>,
}

impl Outcome {
    /// Count of BOM data rows with a **blank** LCSC # (the honest coverage holes).
    fn blank_rows(&self) -> usize {
        self.bom_rows.iter().filter(|f| f[3].is_empty()).count()
    }

    /// True iff at least one BOM row carries a well-formed `C\d+` LCSC #.
    fn has_realizable_row(&self) -> bool {
        self.bom_rows.iter().any(|f| is_lcsc_part(&f[3]))
    }

    /// True iff any designator (BOM or CPL) is a top-C coupling cap (`Cc…`) — the
    /// robust file-based marker that the top-C board was routed (the alternating
    /// ladder has no coupling caps).
    fn has_coupling_caps(&self) -> bool {
        self.bom_designators
            .iter()
            .chain(&self.cpl_designators)
            .any(|d| d.starts_with("Cc"))
    }
}

/// Run `yee filter synth <spec> --jlcpcb <dir>` for a freshly-written spec TOML
/// and read back the artifacts. Asserts the run succeeded and the BOM/CPL headers
/// are the exact JLCPCB columns, then returns the parsed [`Outcome`].
fn run_autoroute(tag: &str, f0_hz: f64, fbw: f64) -> Outcome {
    let base = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("autoroute_{tag}"));
    let _ = std::fs::remove_dir_all(&base);
    std::fs::create_dir_all(&base).expect("create per-case tmpdir");

    // Write the spec TOML and choose an out-dir + throwaway .s2p inside `base`.
    let spec_path = base.join("spec.toml");
    std::fs::write(&spec_path, spec_toml(f0_hz, fbw)).expect("write spec.toml");
    let out_dir = base.join("jlcpcb");
    let s2p = base.join("out.s2p");

    // `--footprint 0402` matches the ADR-0168 specs (the CLI default is 0603);
    // the 0402 routing is what the T4 cases below assert. `--lumped` selects the
    // lumped realization (the JLCPCB BOM/CPL describe lumped LC parts) and skips
    // the orthogonal edge-coupled *distributed* dimensioning, which would
    // otherwise exit 1 for a wideband spec that is lumped-orderable but
    // distributed-unrealizable (e.g. the 1 GHz/70 % case) before the --jlcpcb set
    // is written. `--jlcpcb` itself is independent of `--lumped`.
    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["filter", "synth"])
        .arg(&spec_path)
        .arg("--output")
        .arg(&s2p)
        .arg("--footprint")
        .arg("0402")
        .arg("--lumped")
        .arg("--jlcpcb")
        .arg(&out_dir)
        .output()
        .expect("invoke yee");
    // NB: the process **exit code** reflects the spec-MASK verdict, not the
    // JLCPCB write — `run_synth` writes the upload set BEFORE the PASS/FAIL
    // return, so a mask-FAIL spec (e.g. the wideband 1 GHz/70 % case, whose
    // stopband is unmet) exits 1 yet still emits a complete, correct upload set.
    // This gate validates the JLCPCB artifacts + the routed topology, which are
    // independent of the mask. So we do NOT require exit 0; we require only that
    // the run did not error out (empty stderr — a real failure prints there) and
    // that the JLCPCB stdout markers (only printed once the set is written) are
    // present.
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.is_empty(),
        "[{tag}] run wrote to stderr (a real error, not a mask FAIL): {stderr}"
    );
    let stdout = String::from_utf8_lossy(&output.stdout).into_owned();

    // The run must always name the auto-selected topology and confirm the BOM/CPL
    // were written (these lines print only after the upload set is emitted).
    assert!(
        stdout.contains("topology:") && stdout.contains("(auto-selected)"),
        "[{tag}] run must report the auto-selected topology; stdout:\n{stdout}"
    );
    assert!(
        stdout.contains("wrote JLCPCB BOM:") && stdout.contains("wrote JLCPCB CPL:"),
        "[{tag}] run must report writing the JLCPCB BOM + CPL; stdout:\n{stdout}"
    );

    // ---- BOM ----
    let bom = read_csv(&out_dir.join("bom.csv"), tag, "bom.csv");
    assert_eq!(
        bom[0], "Comment,Designator,Footprint,LCSC Part #",
        "[{tag}] bom.csv header must be the exact JLCPCB assembly-BOM columns"
    );
    let mut bom_rows = Vec::new();
    let mut bom_designators = BTreeSet::new();
    for line in &bom[1..] {
        let f = csv_fields(line);
        assert_eq!(
            f.len(),
            4,
            "[{tag}] bom.csv row must have 4 fields (Comment,Designator,Footprint,LCSC #); got {f:?}"
        );
        for d in designators(&f[1]) {
            assert!(
                !d.is_empty(),
                "[{tag}] empty designator in bom.csv row: {line}"
            );
            bom_designators.insert(d);
        }
        // A blank LCSC # must be flagged in the Comment; a present one must be a
        // well-formed C-number.
        if f[3].is_empty() {
            assert!(
                f[0].contains("NO BASIC PART"),
                "[{tag}] blank-LCSC row must be flagged in its Comment; got {f:?}"
            );
        } else {
            assert!(
                is_lcsc_part(&f[3]),
                "[{tag}] realizable row LCSC # {:?} is not a well-formed C\\d+ part",
                f[3]
            );
        }
        bom_rows.push(f);
    }
    assert!(
        !bom_rows.is_empty(),
        "[{tag}] bom.csv has only a header — no part rows"
    );

    // ---- CPL ----
    let cpl = read_csv(&out_dir.join("cpl.csv"), tag, "cpl.csv");
    assert_eq!(
        cpl[0], "Designator,Mid X,Mid Y,Layer,Rotation",
        "[{tag}] cpl.csv header must be the exact JLCPCB centroid columns"
    );
    let mut cpl_designators = BTreeSet::new();
    for line in &cpl[1..] {
        let f = csv_fields(line);
        assert_eq!(
            f.len(),
            5,
            "[{tag}] cpl.csv row must have 5 fields; got {f:?}"
        );
        assert!(
            !f[0].is_empty(),
            "[{tag}] empty designator in cpl.csv row: {line}"
        );
        cpl_designators.insert(f[0].clone());
    }
    assert!(
        !cpl_designators.is_empty(),
        "[{tag}] cpl.csv has only a header — no rows"
    );

    // CPL ≡ BOM designator sets (both derive from the same placement list, for
    // either topology). Checked for every case.
    assert_eq!(
        cpl_designators, bom_designators,
        "[{tag}] CPL and BOM designator sets disagree — they must both derive from \
         the same placement list (CPL {cpl_designators:?} vs BOM {bom_designators:?})"
    );

    Outcome {
        stdout,
        bom_rows,
        cpl_designators,
        bom_designators,
    }
}

/// Read a written CSV into its non-empty lines, asserting it exists and has data.
fn read_csv(path: &Path, tag: &str, what: &str) -> Vec<String> {
    assert!(
        path.exists(),
        "[{tag}] {what} was not written to {}",
        path.display()
    );
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("[{tag}] read {what}: {e}"));
    let lines: Vec<String> = text.lines().map(str::to_string).collect();
    assert!(
        lines.len() > 1,
        "[{tag}] {what} has only a header — no data rows"
    );
    lines
}

#[test]
fn cli_jlcpcb_autoroute() {
    // =====================================================================
    // CASE 1 — 0.5 GHz / 20 % → top-C, ZERO-blank BOM (the headline rescue).
    // =====================================================================
    let c1 = run_autoroute("subghz_topc", 0.5e9, 0.20);
    // Reported topology (stdout) AND the file marker (Cc coupling caps) must both
    // say top-C — two independent observations of the chosen topology.
    assert!(
        c1.stdout
            .contains("topology: top-C-coupled (capacitively-coupled) (auto-selected)"),
        "[subghz_topc] 0.5 GHz/20% must route to top-C; stdout:\n{}",
        c1.stdout
    );
    assert!(
        c1.has_coupling_caps(),
        "[subghz_topc] top-C board must carry Cc* coupling-cap designators; \
         BOM {:?} CPL {:?}",
        c1.bom_designators,
        c1.cpl_designators
    );
    // Fully orderable: ZERO blank LCSC #s in the written BOM.
    assert_eq!(
        c1.blank_rows(),
        0,
        "[subghz_topc] top-C BOM must have ZERO blank LCSC #s (fully orderable); \
         rows: {:?}",
        c1.bom_rows
    );
    assert!(
        c1.has_realizable_row(),
        "[subghz_topc] top-C BOM must carry at least one realizable C\\d+ row"
    );
    assert!(
        c1.stdout.contains("all ") && c1.stdout.contains("parts matched a JLCPCB Basic part"),
        "[subghz_topc] fully-orderable run must report all parts matched; stdout:\n{}",
        c1.stdout
    );

    // =====================================================================
    // CASE 2 — 1.0 GHz / 70 % → alternating ladder, ZERO-blank BOM.
    // =====================================================================
    let c2 = run_autoroute("wideband_ladder", 1.0e9, 0.70);
    assert!(
        c2.stdout
            .contains("topology: alternating series/shunt ladder (auto-selected)"),
        "[wideband_ladder] 1.0 GHz/70% must route to the alternating ladder; stdout:\n{}",
        c2.stdout
    );
    // The alternating ladder has NO coupling caps (the negative file marker).
    assert!(
        !c2.has_coupling_caps(),
        "[wideband_ladder] alternating ladder must NOT carry Cc* coupling caps; \
         BOM {:?}",
        c2.bom_designators
    );
    assert_eq!(
        c2.blank_rows(),
        0,
        "[wideband_ladder] wideband ladder BOM must have ZERO blank LCSC #s; rows: {:?}",
        c2.bom_rows
    );
    assert!(
        c2.has_realizable_row(),
        "[wideband_ladder] ladder BOM must carry at least one realizable C\\d+ row"
    );
    assert!(
        c2.stdout.contains("all ") && c2.stdout.contains("parts matched a JLCPCB Basic part"),
        "[wideband_ladder] fully-orderable run must report all parts matched; stdout:\n{}",
        c2.stdout
    );

    // =====================================================================
    // CASE 3 — 2.0 GHz / 5 % → NEITHER fully orderable: ≥1 real blank row.
    // =====================================================================
    let c3 = run_autoroute("ghznarrow_blank", 2.0e9, 0.05);
    // Not fully orderable: the honest note must fire and the BOM must carry the
    // real blank rows (flagged, not dropped).
    assert!(
        c3.stdout
            .contains("neither lumped topology is fully orderable"),
        "[ghznarrow_blank] GHz-narrow spec must fire the distributed/planar-track \
         note; stdout:\n{}",
        c3.stdout
    );
    assert!(
        c3.stdout.contains("NOTE:") && c3.stdout.contains("no JLCPCB Basic match"),
        "[ghznarrow_blank] not-orderable run must print the honest no-match NOTE; \
         stdout:\n{}",
        c3.stdout
    );
    assert!(
        c3.blank_rows() >= 1,
        "[ghznarrow_blank] not-orderable BOM must carry ≥1 real blank LCSC row \
         (honest, not dropped); rows: {:?}",
        c3.bom_rows
    );
    // The honest board still emits realizable parts alongside the blanks (the hole
    // is partial, never a fabricated all-orderable board).
    assert!(
        c3.has_realizable_row(),
        "[ghznarrow_blank] the honest fewer-blanks board still has realizable rows"
    );
}
