//! JLCPCB **assembly upload** CSV export — BOM + CPL/centroid (JLCPCB
//! production track, ADR-0164 bricks **J2** + **J3**).
//!
//! Turns a placed lumped-LC board ([`crate::lumped_board`] → [`Placement`]s) and
//! its synthesized component values ([`crate::LumpedLadder`]) into the two CSVs
//! JLCPCB's SMT-assembly order consumes:
//!
//! - **BOM CSV** ([`jlcpcb_bom_csv`]) — header `Comment,Designator,Footprint,LCSC
//!   Part #`, one row per **distinct** LCSC part (the ref-des that share a part
//!   are comma-joined into the `Designator` cell).
//! - **CPL / centroid CSV** ([`jlcpcb_cpl_csv`]) — header `Designator,Mid X,Mid
//!   Y,Layer,Rotation`, one row per placed component.
//!
//! Pure data / string + `f64`, **WASM-safe** (no I/O, no network, no non-WASM
//! dep), deterministic — the same constraint as the rest of `yee-filter`. The
//! studio (the deliverable) can hand these straight to the user as the JLCPCB
//! upload set, and the J4 CLI calls [`jlcpcb_files`] to write `bom.csv` +
//! `cpl.csv`.
//!
//! # The placement ↔ value join (the crux)
//!
//! A [`Placement`] is a *physical* component (ref-des `L1`/`C1`/`L2`/…, a
//! footprint, a board position) but carries **no value** — the values live on
//! the [`LumpedLadder`]'s resonators. [`lumped_board`](crate::lumped_board) emits
//! the placements in ladder order, `L<k>` then `C<k>` for resonator `k` (1-based,
//! `resonators[k-1]`): `L<k>` realizes that resonator's `l_henry`, `C<k>` its
//! `c_farad`. So the join is the **ref-des index**: parse the trailing integer
//! `k` off the ref-des, index `resonators[k-1]`, and read `l_henry` (for an `L`
//! prefix) or `c_farad` (for a `C` prefix). The raw reactance is then snapped to
//! the requested [`ESeries`] with the **same** [`ESeries::nearest`] call
//! [`select_components`](crate::select_components) makes, so the per-placement
//! chosen value and the [`autopick`] result match the grouped [`Bom`] exactly.
//! (The final BOM CSV then groups by **orderable LCSC part #** — one row per
//! distinct part, the right key for an assembly BOM — which can be coarser than
//! the [`Bom`]'s `(kind, chosen_value)` key when two nearby E-series values map to
//! the same in-table part within the autopick band.) [`PlacedPart`] is the join.
//!
//! # JLCPCB format conventions (documented; confirm against a real upload)
//!
//! JLCPCB's assembly templates are lenient but these are the conventions used:
//!
//! - **BOM header** — `Comment,Designator,Footprint,LCSC Part #` (JLCPCB's
//!   canonical assembly-BOM columns; `Comment` is the human value, `Designator`
//!   a comma-joined ref-des list, `Footprint` the land name, `LCSC Part #` the
//!   C-number). Extra columns are ignored by JLCPCB; these four are the ones it
//!   reads.
//! - **Footprint name** — the bare imperial code (`0402`, `0603`, `0805`) via
//!   [`jlcpcb_footprint_name`]. JLCPCB accepts both `0603` and its library name
//!   `C0603`/`R0603`; the bare code is the most portable and is what their BOM
//!   examples show for generic chips, so we emit the bare code.
//! - **CPL header** — `Designator,Mid X,Mid Y,Layer,Rotation` (JLCPCB's exact
//!   centroid columns).
//! - **CPL coordinates** — millimetres with a `mm` unit suffix (e.g. `3.500mm`).
//!   JLCPCB's CPL parser accepts an explicit unit suffix and a bare number (it
//!   defaults bare numbers to mm); the explicit `mm` suffix is unambiguous, so we
//!   emit it.
//! - **Layer** — `Top` (the lumped board is single-layer top copper).
//! - **Rotation** — degrees, `0` for the axis-aligned walking-skeleton placement
//!   ([`lumped_board`] places every footprint axis-aligned). Per-component
//!   orientation is a J-track follow-on; `0` is correct for the current geometry.
//!
//! Unverified detail flagged for the maintainer to confirm against a real JLCPCB
//! upload: the exact accepted **footprint string** (bare `0603` vs `C0603`) and
//! whether the `mm` suffix or a bare number is preferred — both are documented
//! here and changeable in one place ([`jlcpcb_footprint_name`] / the CPL
//! formatter) if a real upload rejects them.

use crate::{
    BomLine, CompKind, ESeries, Footprint, LcResonator, LumpedLadder, Placement, TopCNetwork,
    autopick, board::BranchKind,
};

/// One placed component joined to its value + autopicked LCSC part.
///
/// The output of the [placement ↔ value join](self#the-placement--value-join-the-crux):
/// a physical [`Placement`] (ref-des, footprint, branch, board centre) carrying
/// the resonator value it realizes (`chosen_value`, E-series-snapped), a
/// [`BomLine`]-equivalent for [`autopick`], and the resulting LCSC part (`None`
/// when no in-table Basic part matches — an honest, surfaced coverage hole, e.g.
/// a narrow-band band-pass series resonator's sub-pF capacitor).
#[derive(Debug, Clone, PartialEq)]
pub struct PlacedPart {
    /// Reference designator, e.g. `"L1"` / `"C1"`.
    pub ref_des: String,
    /// Inductor or capacitor (from the ref-des prefix / resonator field).
    pub kind: CompKind,
    /// Series (in-line) or shunt (stub-to-ground) branch role.
    pub branch: BranchKind,
    /// The chosen E-series value (H / F) realizing this component — snapped from
    /// the resonator reactance with the same [`ESeries::nearest`] as
    /// [`select_components`](crate::select_components).
    pub chosen_value: f64,
    /// The SMD land pattern.
    pub footprint: Footprint,
    /// Footprint centre `(x, y)`, metres, in the board frame.
    pub center_m: (f64, f64),
    /// The autopicked real LCSC part, or `None` if no in-table part matches.
    pub lcsc: Option<crate::LcscPart>,
}

/// Join each [`Placement`] to its resonator value and autopicked LCSC part.
///
/// The single source of truth behind both CSVs: for every placement, parse the
/// ref-des index `k`, read `resonators[k-1]`'s `l_henry`/`c_farad` (by the `L`/`C`
/// prefix), snap to `series` ([`ESeries::nearest`]), build a [`BomLine`] and
/// [`autopick`] an LCSC part. See the [module docs](self#the-placement--value-join-the-crux)
/// for why the ref-des index is the correct join key.
///
/// Placements whose ref-des does not parse to a valid `1..=N` index into the
/// ladder are skipped (they cannot be value-joined); for a board produced by
/// [`lumped_board`] every placement joins, so the output length equals
/// `placements.len()` in practice.
pub fn join_placed_parts(
    placements: &[Placement],
    ladder: &LumpedLadder,
    footprint: Footprint,
    series: ESeries,
) -> Vec<PlacedPart> {
    let tolerance_pct = series.tolerance_pct();
    let mut out = Vec::with_capacity(placements.len());

    for p in placements {
        let Some((kind, k)) = parse_ref_des(&p.ref_des) else {
            continue; // unparseable ref-des — cannot value-join; skip honestly.
        };
        // 1-based ladder index → resonator.
        let Some(LcResonator {
            l_henry, c_farad, ..
        }) = ladder.resonators.get(k - 1)
        else {
            continue; // ref-des index out of range for this ladder; skip.
        };
        // The reactance this physical component realizes (H for L, F for C).
        let ideal_value = match kind {
            CompKind::Inductor => *l_henry,
            CompKind::Capacitor => *c_farad,
        };
        // Snap with the SAME selector select_components uses, so the chosen value
        // (and the autopick + grouping) is identical to the grouped Bom.
        let chosen_value = series.nearest(ideal_value);
        let deviation_pct = if ideal_value != 0.0 {
            (chosen_value - ideal_value) / ideal_value * 100.0
        } else {
            0.0
        };
        let line = BomLine {
            kind,
            ideal_value,
            chosen_value,
            deviation_pct,
            series,
            tolerance_pct,
            qty: 1,
            esr_ohm: None,
            srf_hz: None,
        };
        let lcsc = autopick(&line, footprint);
        out.push(PlacedPart {
            ref_des: p.ref_des.clone(),
            kind,
            branch: p.kind,
            chosen_value,
            footprint,
            center_m: p.center_m,
            lcsc,
        });
    }
    out
}

/// Parse a ref-des like `"L12"` / `"C3"` into its component kind + 1-based index.
///
/// Returns `None` if the string is not a single-letter `L`/`C` prefix followed by
/// a positive integer.
fn parse_ref_des(ref_des: &str) -> Option<(CompKind, usize)> {
    let mut chars = ref_des.chars();
    let kind = match chars.next()? {
        'L' => CompKind::Inductor,
        'C' => CompKind::Capacitor,
        _ => return None,
    };
    let digits = chars.as_str();
    if digits.is_empty() || !digits.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }
    let k: usize = digits.parse().ok()?;
    if k == 0 {
        return None;
    }
    Some((kind, k))
}

/// The arm a top-C ref-des belongs to: a shunt-resonator element or a coupling
/// cap.
///
/// A [`TopCNetwork`] board ([`crate::top_c_board`]) has three ref-des families,
/// which [`parse_ref_des`] (used for the alternating ladder) cannot tell apart —
/// it would mis-read `Cc1` as a capacitor `C` with a non-digit body. So the
/// top-C join uses its own parser ([`parse_top_c_ref_des`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TopCRef {
    /// Shunt-resonator inductor `L{k}` → `shunt[k-1].l_henry` (1-based `k`).
    ShuntL(usize),
    /// Shunt-resonator capacitor `C{k}` → `shunt[k-1].c_farad` (1-based `k`).
    ShuntC(usize),
    /// Series coupling capacitor `Cc{j}` → `coupling_caps_farad[j-1]` (1-based `j`).
    CouplingCap(usize),
}

/// Parse a top-C board ref-des into its arm + 1-based index.
///
/// Recognizes `L{k}` / `C{k}` (shunt resonator `k`) and `Cc{j}` (coupling cap
/// `j`); the `Cc` arm must be checked **before** the bare `C` arm (`Cc1` starts
/// with `C`). Returns `None` for any other shape (e.g. `0`-index or junk).
fn parse_top_c_ref_des(ref_des: &str) -> Option<TopCRef> {
    // Coupling cap `Cc{digits}` — check first (it starts with 'C' like a shunt C).
    if let Some(rest) = ref_des.strip_prefix("Cc") {
        let j = parse_positive_index(rest)?;
        return Some(TopCRef::CouplingCap(j));
    }
    // Shunt resonator inductor / capacitor.
    let mut chars = ref_des.chars();
    let arm = match chars.next()? {
        'L' => TopCRef::ShuntL,
        'C' => TopCRef::ShuntC,
        _ => return None,
    };
    let k = parse_positive_index(chars.as_str())?;
    Some(arm(k))
}

/// Parse a string of ASCII digits into a strictly-positive (1-based) index, or
/// `None` if empty / non-digit / zero.
fn parse_positive_index(digits: &str) -> Option<usize> {
    if digits.is_empty() || !digits.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }
    let k: usize = digits.parse().ok()?;
    if k == 0 { None } else { Some(k) }
}

/// Join each top-C board [`Placement`] to its [`TopCNetwork`] value + autopicked
/// LCSC part (top-C-coupled board, ADR-0166 brick **T2**).
///
/// The top-C analogue of [`join_placed_parts`]: top-C has a different topology
/// (`N` shunt LC resonators + `N+1` series coupling caps, [`crate::top_c`]) and
/// three ref-des families, so it parses with [`parse_top_c_ref_des`] (NOT the
/// alternating-ladder [`parse_ref_des`], which cannot distinguish a coupling
/// cap `Cc{j}` from a shunt cap `C{k}`) and indexes the [`TopCNetwork`]:
///
/// - `L{k}` → `net.shunt[k-1].l_henry`  ([`CompKind::Inductor`])
/// - `C{k}` → `net.shunt[k-1].c_farad`  ([`CompKind::Capacitor`])
/// - `Cc{j}` → `net.coupling_caps_farad[j-1]` ([`CompKind::Capacitor`])
///
/// Each value is snapped to `series` ([`ESeries::nearest`]) and [`autopick`]ed
/// exactly as [`join_placed_parts`] does (capacitors via [`CompKind::Capacitor`],
/// inductors via [`CompKind::Inductor`]), building a [`PlacedPart`] reusing
/// [`value_comment`]. Unmatched values keep their `lcsc = None` so the BOM emits
/// the same honest blank-LCSC `(NO BASIC PART)` row as [`join_placed_parts`] —
/// never dropped or faked. The resulting `Vec<PlacedPart>` feeds the **same**
/// [`jlcpcb_bom_csv`] / [`jlcpcb_cpl_csv`] emitters unchanged.
///
/// Placements whose ref-des does not parse to a valid in-range index are skipped
/// (they cannot be value-joined); for a board produced by [`crate::top_c_board`]
/// every placement joins, so the output length equals `placements.len()`.
pub fn join_top_c_parts(placements: &[Placement], net: &TopCNetwork) -> Vec<PlacedPart> {
    // The top-C board places one footprint family (0402/0603/…); take it from the
    // placements (every placement carries the footprint top_c_board was called
    // with). Use E24 — the same E-series the lumped autopick path snaps to before
    // matching against the coarser stocked LCSC grid.
    let series = ESeries::E24;
    let tolerance_pct = series.tolerance_pct();
    let mut out = Vec::with_capacity(placements.len());

    for p in placements {
        let Some(parsed) = parse_top_c_ref_des(&p.ref_des) else {
            continue; // unparseable ref-des — cannot value-join; skip honestly.
        };
        // Resolve (kind, ideal SI value) by arm, with an in-range index check.
        let (kind, ideal_value) = match parsed {
            TopCRef::ShuntL(k) => {
                let Some(res) = net.shunt.get(k - 1) else {
                    continue; // index out of range for this network; skip.
                };
                (CompKind::Inductor, res.l_henry)
            }
            TopCRef::ShuntC(k) => {
                let Some(res) = net.shunt.get(k - 1) else {
                    continue;
                };
                (CompKind::Capacitor, res.c_farad)
            }
            TopCRef::CouplingCap(j) => {
                let Some(cc) = net.coupling_caps_farad.get(j - 1) else {
                    continue;
                };
                (CompKind::Capacitor, *cc)
            }
        };
        // Snap with the SAME selector as join_placed_parts so the chosen value
        // (and the autopick + grouping) is identical to the grouped Bom.
        let chosen_value = series.nearest(ideal_value);
        let deviation_pct = if ideal_value != 0.0 {
            (chosen_value - ideal_value) / ideal_value * 100.0
        } else {
            0.0
        };
        let line = BomLine {
            kind,
            ideal_value,
            chosen_value,
            deviation_pct,
            series,
            tolerance_pct,
            qty: 1,
            esr_ohm: None,
            srf_hz: None,
        };
        let lcsc = autopick(&line, p.footprint);
        out.push(PlacedPart {
            ref_des: p.ref_des.clone(),
            kind,
            branch: p.kind,
            chosen_value,
            footprint: p.footprint,
            center_m: p.center_m,
            lcsc,
        });
    }
    out
}

/// Human value label for a BOM `Comment` cell, e.g. `"22pF"` / `"47nH"`.
///
/// Picks a sensible engineering unit (pF/nF for capacitors, nH/µH for inductors)
/// and trims a trailing `.0` so common values read cleanly (`22pF`, not
/// `22.0pF`). Sub-unit values keep enough significant figures to be distinct
/// (e.g. a `0.15pF` series-resonator cap).
pub fn value_comment(kind: CompKind, value_si: f64) -> String {
    let (scaled, unit) = match kind {
        CompKind::Capacitor => {
            let pf = value_si * 1e12;
            if pf >= 1000.0 {
                (pf / 1000.0, "nF")
            } else {
                (pf, "pF")
            }
        }
        CompKind::Inductor => {
            let nh = value_si * 1e9;
            if nh >= 1000.0 {
                (nh / 1000.0, "uH")
            } else {
                (nh, "nH")
            }
        }
    };
    format!("{}{}", trim_number(scaled), unit)
}

/// Format a number with up to 3 significant decimals, trimming trailing zeros
/// and a trailing dot (`22.0 → "22"`, `0.15 → "0.15"`, `4.7 → "4.7"`).
fn trim_number(v: f64) -> String {
    if !v.is_finite() {
        return format!("{v}");
    }
    // 3 decimals is enough to distinguish the filter L/C grid (and sub-pF caps);
    // strip trailing zeros / dot for readability.
    let mut s = format!("{v:.3}");
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }
    s
}

/// CSV-escape one field per RFC 4180: wrap in double quotes (doubling any inner
/// quote) iff it contains a comma, a quote, a CR, or an LF.
fn csv_escape(field: &str) -> String {
    if field.contains([',', '"', '\n', '\r']) {
        let escaped = field.replace('"', "\"\"");
        format!("\"{escaped}\"")
    } else {
        field.to_string()
    }
}

/// Join a row of already-stringified fields into a CSV line (each field escaped).
fn csv_row(fields: &[&str]) -> String {
    fields
        .iter()
        .map(|f| csv_escape(f))
        .collect::<Vec<_>>()
        .join(",")
}

/// The JLCPCB footprint string for a [`Footprint`].
///
/// The bare imperial chip code (`"0402"` / `"0603"` / `"0805"`). See the
/// [module docs](self#jlcpcb-format-conventions-documented-confirm-against-a-real-upload)
/// for why the bare code (vs the `C0603` library name) is emitted.
pub fn jlcpcb_footprint_name(footprint: Footprint) -> &'static str {
    match footprint {
        Footprint::Smd0402 => "0402",
        Footprint::Smd0603 => "0603",
        Footprint::Smd0805 => "0805",
    }
}

/// The exact JLCPCB BOM CSV header line (no trailing newline).
pub const BOM_HEADER: &str = "Comment,Designator,Footprint,LCSC Part #";

/// The exact JLCPCB CPL/centroid CSV header line (no trailing newline).
pub const CPL_HEADER: &str = "Designator,Mid X,Mid Y,Layer,Rotation";

/// Render the **JLCPCB assembly BOM CSV** for a list of [joined parts](join_placed_parts).
///
/// Header `Comment,Designator,Footprint,LCSC Part #` (exact, [`BOM_HEADER`]), then
/// one row per **distinct LCSC part**: parts that resolved to the same C-number
/// are grouped into a single row whose `Designator` is their comma-joined,
/// sorted ref-des (e.g. `"C1,C3"`). `Comment` is the human value
/// ([`value_comment`]); `Footprint` is the JLCPCB land name
/// ([`jlcpcb_footprint_name`]).
///
/// **Unrealizable parts are emitted, not dropped** (honest coverage): every
/// placed component whose [`autopick`](crate::autopick) returned `None` gets its
/// own row with a **blank** `LCSC Part #` and a `Comment` flagged `… (NO BASIC
/// PART)`, so the user sees the unfillable line rather than a silently short BOM.
/// Such `None` parts are grouped by `(kind, chosen_value)` so identical
/// unrealizable values share one flagged row.
///
/// Rows are ordered: realizable parts first (by ascending LCSC C-number for
/// determinism), then the unrealizable rows (by `Comment`). Every row has exactly
/// four CSV fields; fields containing commas are quoted.
pub fn jlcpcb_bom_csv(parts: &[PlacedPart]) -> String {
    // Group realizable parts by LCSC C-number; group unrealizable (None) parts by
    // (kind, chosen-value) so identical off-catalog values share one flagged row.
    struct Group {
        comment: String,
        designators: Vec<String>,
        footprint: Footprint,
        lcsc: Option<String>,
    }
    let mut realizable: Vec<Group> = Vec::new();
    let mut unrealizable: Vec<Group> = Vec::new();

    for p in parts {
        match &p.lcsc {
            Some(part) => {
                if let Some(g) = realizable
                    .iter_mut()
                    .find(|g| g.lcsc.as_deref() == Some(part.lcsc))
                {
                    g.designators.push(p.ref_des.clone());
                } else {
                    realizable.push(Group {
                        // Use the human value of the chosen part value (== picked
                        // part's value within tolerance) for the Comment.
                        comment: value_comment(p.kind, p.chosen_value),
                        designators: vec![p.ref_des.clone()],
                        footprint: p.footprint,
                        lcsc: Some(part.lcsc.to_string()),
                    });
                }
            }
            None => {
                // Group None parts by (kind, chosen value) — identical
                // unrealizable values share one flagged row.
                let matches = |g: &&mut Group| {
                    g.lcsc.is_none()
                        && g.comment == flagged_comment(p.kind, p.chosen_value)
                        && g.footprint == p.footprint
                };
                if let Some(g) = unrealizable.iter_mut().find(matches) {
                    g.designators.push(p.ref_des.clone());
                } else {
                    unrealizable.push(Group {
                        comment: flagged_comment(p.kind, p.chosen_value),
                        designators: vec![p.ref_des.clone()],
                        footprint: p.footprint,
                        lcsc: None,
                    });
                }
            }
        }
    }

    // Deterministic ordering: realizable by C-number, then unrealizable by comment.
    realizable.sort_by(|a, b| a.lcsc.cmp(&b.lcsc));
    unrealizable.sort_by(|a, b| a.comment.cmp(&b.comment));

    let mut lines = vec![BOM_HEADER.to_string()];
    for g in realizable.into_iter().chain(unrealizable.into_iter()) {
        let mut des = g.designators;
        des.sort_by_key(|a| ref_des_sort_key(a));
        let designator = des.join(",");
        let footprint = jlcpcb_footprint_name(g.footprint);
        let lcsc = g.lcsc.unwrap_or_default(); // blank for unrealizable
        lines.push(csv_row(&[&g.comment, &designator, footprint, &lcsc]));
    }
    lines.join("\n")
}

/// A `Comment` for an unrealizable (no-LCSC) part, flagged so the user sees it.
fn flagged_comment(kind: CompKind, value_si: f64) -> String {
    format!("{} (NO BASIC PART)", value_comment(kind, value_si))
}

/// Sort key for ref-des so `L2` precedes `L10` (letter, then numeric index).
fn ref_des_sort_key(ref_des: &str) -> (char, u64, String) {
    match parse_ref_des(ref_des) {
        Some((kind, k)) => {
            let letter = match kind {
                CompKind::Inductor => 'L',
                CompKind::Capacitor => 'C',
            };
            (letter, k as u64, String::new())
        }
        // Unparseable ref-des sort last, lexically, deterministically.
        None => ('~', u64::MAX, ref_des.to_string()),
    }
}

/// Render the **JLCPCB CPL / centroid CSV** for the placement list.
///
/// Header `Designator,Mid X,Mid Y,Layer,Rotation` (exact, [`CPL_HEADER`]), then
/// one row per [`Placement`] in the given order: `Designator` = ref-des, `Mid X`
/// / `Mid Y` = the footprint centre in **millimetres** with an `mm` suffix
/// ([format note](self#jlcpcb-format-conventions-documented-confirm-against-a-real-upload)),
/// `Layer` = `Top`, `Rotation` = `0` (axis-aligned placement). Every row has five
/// CSV fields.
///
/// This takes the raw [`Placement`]s (not [`PlacedPart`]s) so the centroid file
/// lists **every** placed component regardless of whether its value resolved to
/// an LCSC part — the CPL describes physical positions, and JLCPCB needs a
/// centroid line for each board component. (Designator consistency with the BOM
/// is the gate's job; both derive from the same placement set.)
pub fn jlcpcb_cpl_csv(placements: &[Placement]) -> String {
    let mut lines = vec![CPL_HEADER.to_string()];
    for p in placements {
        let (x_m, y_m) = p.center_m;
        let mid_x = format_mm(x_m);
        let mid_y = format_mm(y_m);
        lines.push(csv_row(&[&p.ref_des, &mid_x, &mid_y, "Top", "0"]));
    }
    lines.join("\n")
}

/// Format a metre coordinate as millimetres with an `mm` suffix, 3 decimals
/// (micron resolution), e.g. `3.500mm`.
fn format_mm(v_m: f64) -> String {
    format!("{:.3}mm", v_m * 1e3)
}

/// The JLCPCB upload CSV pair: the assembly BOM + the CPL/centroid.
///
/// Returned by [`jlcpcb_files`] for the J4 CLI / studio to write as `bom.csv` and
/// `cpl.csv`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JlcpcbFiles {
    /// The JLCPCB assembly BOM CSV ([`jlcpcb_bom_csv`]).
    pub bom_csv: String,
    /// The JLCPCB CPL / centroid CSV ([`jlcpcb_cpl_csv`]).
    pub cpl_csv: String,
}

/// One-call convenience: a placed board + its ladder → the JLCPCB upload pair.
///
/// Joins the placements to their values + LCSC parts ([`join_placed_parts`]) and
/// renders both CSVs. The CPL lists every placement; the BOM lists every distinct
/// LCSC part (with honest blank-LCSC rows for unrealizable values). This is the
/// J4 entry point — the CLI synthesizes → `lumped_board` → calls this → writes
/// `bom.csv` + `cpl.csv`.
///
/// `placements` and `ladder` must come from the same synthesis (the placements
/// are value-joined to the ladder by ref-des index); pass the
/// [`lumped_board`](crate::lumped_board) output's `placements` and the
/// [`synthesize_lumped`](crate::synthesize_lumped) ladder. `footprint` must match
/// the one `lumped_board` was called with (it sets the BOM `Footprint` and the
/// autopick land), and `series` the BOM E-series.
pub fn jlcpcb_files(
    placements: &[Placement],
    ladder: &LumpedLadder,
    footprint: Footprint,
    series: ESeries,
) -> JlcpcbFiles {
    let parts = join_placed_parts(placements, ladder, footprint, series);
    JlcpcbFiles {
        bom_csv: jlcpcb_bom_csv(&parts),
        cpl_csv: jlcpcb_cpl_csv(placements),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_ref_des_basic() {
        assert_eq!(parse_ref_des("L1"), Some((CompKind::Inductor, 1)));
        assert_eq!(parse_ref_des("C12"), Some((CompKind::Capacitor, 12)));
        assert_eq!(parse_ref_des("L0"), None); // 0 is not a valid 1-based index
        assert_eq!(parse_ref_des("X3"), None); // unknown prefix
        assert_eq!(parse_ref_des("L"), None); // no digits
        assert_eq!(parse_ref_des("LC3"), None); // non-digit body
        assert_eq!(parse_ref_des(""), None);
    }

    #[test]
    fn csv_escape_quotes_commas_and_quotes() {
        assert_eq!(csv_escape("0603"), "0603");
        assert_eq!(csv_escape("C1,C3"), "\"C1,C3\""); // comma → quoted
        assert_eq!(csv_escape("a\"b"), "\"a\"\"b\""); // inner quote doubled
        assert_eq!(csv_escape("line\nbreak"), "\"line\nbreak\""); // newline → quoted
    }

    #[test]
    fn value_comment_units() {
        assert_eq!(value_comment(CompKind::Capacitor, 22e-12), "22pF");
        assert_eq!(value_comment(CompKind::Capacitor, 1.5e-12), "1.5pF");
        assert_eq!(value_comment(CompKind::Inductor, 47e-9), "47nH");
        assert_eq!(value_comment(CompKind::Capacitor, 2.2e-9), "2.2nF"); // ≥1000 pF → nF
        assert_eq!(value_comment(CompKind::Inductor, 1.0e-6), "1uH"); // ≥1000 nH → µH
    }

    #[test]
    fn footprint_names() {
        assert_eq!(jlcpcb_footprint_name(Footprint::Smd0402), "0402");
        assert_eq!(jlcpcb_footprint_name(Footprint::Smd0603), "0603");
        assert_eq!(jlcpcb_footprint_name(Footprint::Smd0805), "0805");
    }

    /// A `PlacedPart` with a resolved LCSC C-number.
    fn part(ref_des: &str, kind: CompKind, lcsc: &'static str, value: f64) -> PlacedPart {
        PlacedPart {
            ref_des: ref_des.to_string(),
            kind,
            branch: BranchKind::Shunt,
            chosen_value: value,
            footprint: Footprint::Smd0603,
            center_m: (0.0, 0.0),
            lcsc: Some(crate::LcscPart {
                lcsc,
                kind,
                value,
                footprint: Footprint::Smd0603,
                basic: true,
                tolerance_pct: Some(5.0),
                voltage: Some(50.0),
            }),
        }
    }

    /// A `PlacedPart` with NO resolved LCSC part (the honest coverage hole).
    fn none_part(ref_des: &str, kind: CompKind, value: f64) -> PlacedPart {
        PlacedPart {
            ref_des: ref_des.to_string(),
            kind,
            branch: BranchKind::Series,
            chosen_value: value,
            footprint: Footprint::Smd0603,
            center_m: (0.0, 0.0),
            lcsc: None,
        }
    }

    #[test]
    fn bom_groups_same_part_into_one_row() {
        // Two capacitors that resolved to the SAME C-number must collapse into one
        // BOM row whose Designator lists both ref-des (sorted).
        let parts = vec![
            part("C1", CompKind::Capacitor, "C1653", 22e-12),
            part("C3", CompKind::Capacitor, "C1653", 22e-12),
        ];
        let csv = jlcpcb_bom_csv(&parts);
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines[0], BOM_HEADER);
        assert_eq!(lines.len(), 2, "two same-part caps → one data row");
        assert_eq!(lines[1], "22pF,\"C1,C3\",0603,C1653");
    }

    #[test]
    fn bom_emits_unrealizable_row_not_dropped() {
        // A None part must appear with a blank LCSC # and a flagged comment.
        let parts = vec![
            part("C1", CompKind::Capacitor, "C1653", 22e-12),
            none_part("C2", CompKind::Capacitor, 0.15e-12),
        ];
        let csv = jlcpcb_bom_csv(&parts);
        let lines: Vec<&str> = csv.lines().collect();
        // header + realizable + unrealizable = 3 lines.
        assert_eq!(lines.len(), 3, "the None part must NOT be dropped");
        // Realizable row first (sorted by C-number), then the flagged blank row.
        assert_eq!(lines[1], "22pF,C1,0603,C1653");
        assert_eq!(lines[2], "0.15pF (NO BASIC PART),C2,0603,");
        // The last row's 4th field (LCSC #) is blank.
        assert!(lines[2].ends_with(','), "unrealizable LCSC # must be blank");
    }

    #[test]
    fn bom_groups_identical_unrealizable_values() {
        // Two None parts with the same (kind, value) share one flagged row.
        let parts = vec![
            none_part("C2", CompKind::Capacitor, 0.15e-12),
            none_part("C4", CompKind::Capacitor, 0.15e-12),
        ];
        let csv = jlcpcb_bom_csv(&parts);
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 2, "identical unrealizable values → one row");
        assert_eq!(lines[1], "0.15pF (NO BASIC PART),\"C2,C4\",0603,");
    }

    #[test]
    fn cpl_row_per_placement_with_mm() {
        let placements = vec![
            Placement {
                ref_des: "L1".to_string(),
                footprint: Footprint::Smd0603,
                kind: BranchKind::Shunt,
                center_m: (3.5e-3, 1.25e-3),
            },
            Placement {
                ref_des: "C1".to_string(),
                footprint: Footprint::Smd0603,
                kind: BranchKind::Shunt,
                center_m: (5.0e-3, 1.25e-3),
            },
        ];
        let csv = jlcpcb_cpl_csv(&placements);
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines[0], CPL_HEADER);
        assert_eq!(lines.len(), 3, "header + 2 placements");
        assert_eq!(lines[1], "L1,3.500mm,1.250mm,Top,0");
        assert_eq!(lines[2], "C1,5.000mm,1.250mm,Top,0");
    }

    #[test]
    fn ref_des_sort_is_numeric() {
        // L2 must precede L10 (numeric, not lexical).
        let mut v = vec!["L10".to_string(), "L2".to_string(), "C1".to_string()];
        v.sort_by_key(|a| ref_des_sort_key(a));
        assert_eq!(v, vec!["C1", "L2", "L10"]);
    }
}
