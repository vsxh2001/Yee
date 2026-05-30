//! Engine bridge — drives the live Yee filter engine for the POC.
//!
//! The POC hard-codes the committed `cheb_bpf.toml` fixture as a
//! [`FilterSpec`] (no file IO in WASM) and runs the **real** synthesis /
//! dimensioning / response-sweep / spec-mask paths through `yee-synth`,
//! `yee-filter`, and `yee-layout` — the same calls the `yee filter synth` CLI
//! makes ([`yee_cli::filter::run_synth`], mirrored here). Nothing on these
//! stages is faked: every number rendered by the Synthesis and Layout stages
//! comes out of [`Designed`].

use yee_filter::{
    Approximation, Bom, BranchKind, CompKind, CouplingMatrix, ESeries, FilterProject, FilterSpec,
    Footprint, LcBranch, LumpedBoard, LumpedLadder, MaskReport, MaskVerdict, Response, SpecMask,
    check_mask, dimension_edge_coupled, dimension_edge_coupled_layout, ideal_response, ladder_s21,
    lumped_board, mask_verdict, monte_carlo_yield, select_components, synthesize,
    synthesize_lumped,
};
use yee_layout::{BBox, CoupledMicrostrip, Layout, Substrate, coupled_microstrip, eps_eff};

/// Number of points in the response sweep (matches the CLI's `SWEEP_POINTS`).
const SWEEP_POINTS: usize = 401;
/// Sweep span as a multiple of FBW on each side of `f0` (matches the CLI).
const SPAN_MULT: f64 = 6.0;

/// The committed FR-4 substrate the CLI defaults to for dimensioning
/// (`--eps-r 4.4 --h-mm 1.6`), with a representative 35 µm copper + loss
/// tangent for the stackup display (these last two do not enter the
/// closed-form dimensioning — see [`dimension_edge_coupled`]).
pub const SUBSTRATE: Substrate = Substrate {
    eps_r: 4.4,
    height_m: 1.6e-3,
    loss_tangent: 0.02,
    metal_thickness_m: 35e-6,
};

/// The `cheb_bpf.toml` fixture, hard-coded (WASM has no file IO): a 0.5 dB
/// Chebyshev, 2 GHz centre, 10 % fractional bandwidth, order 5, with a 40 dB
/// rejection point at 2.4 GHz. Mirrors the committed CLI fixture exactly.
pub fn demo_spec() -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: 2.0e9,
        fbw: 0.10,
        order: Some(5),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 9.0,
            stopband: vec![(2.4e9, 40.0)],
        },
    }
}

/// A single swept response sample (magnitudes only — the closed-form model is
/// magnitude-only; `|S11|² = 1 − |S21|²` by losslessness).
#[derive(Clone, Copy)]
pub struct SweepPoint {
    /// Frequency, Hz.
    pub f_hz: f64,
    /// `|S21|` in dB.
    pub s21_db: f64,
    /// `|S11|` in dB.
    pub s11_db: f64,
}

/// Forbidden mask region for the response plot (mirrors the CLI's
/// `spec_mask_regions`): a passband insertion-loss floor + per-stopband
/// rejection ceilings.
#[derive(Clone, Copy)]
pub struct MaskBand {
    /// Lower frequency edge, Hz.
    pub f_lo_hz: f64,
    /// Upper frequency edge, Hz.
    pub f_hi_hz: f64,
    /// `true` = a "floor" (passband, forbidden *below* `limit_db`); `false` =
    /// a "ceiling" (stopband, forbidden *above* `limit_db`).
    pub is_floor: bool,
    /// The dB limit of the forbidden region.
    pub limit_db: f64,
}

/// Per-resonator realized geometry + electrical parameters for the Layout
/// stage's components table. `gap_to_next_mm`/`z0e_ohm`/… are `None` for the
/// last resonator (which has no successor gap).
#[derive(Clone)]
pub struct ResonatorRow {
    /// Resonator id, 1-based.
    pub id: usize,
    /// Strip width, mm.
    pub width_mm: f64,
    /// Resonator length (`λ_g/2`), mm.
    pub length_mm: f64,
    /// Edge-coupling gap to the next resonator, mm (`None` for the last).
    pub gap_to_next_mm: Option<f64>,
    /// Even-mode characteristic impedance, Ω (`None` for the last).
    pub z0e_ohm: Option<f64>,
    /// Odd-mode characteristic impedance, Ω (`None` for the last).
    pub z0o_ohm: Option<f64>,
    /// Even-mode effective permittivity (`None` for the last).
    pub eps_eff_e: Option<f64>,
    /// Odd-mode effective permittivity (`None` for the last).
    pub eps_eff_o: Option<f64>,
    /// Realized coupling `k = (Z0e − Z0o)/(Z0e + Z0o)` (`None` for the last).
    pub realized_k: Option<f64>,
    /// Target coupling `FBW · m_{i,i+1}` the gap was solved for (`None` last).
    pub target_k: Option<f64>,
}

/// Everything the POC's two real stages render — all from the live engine.
pub struct Designed {
    /// The hard-coded demo spec.
    pub spec: FilterSpec,
    /// The synthesized project (prototype, coupling matrix, topology).
    pub project: FilterProject,
    /// Prototype element values `[g0, g1, …, gN, g_{N+1}]`.
    pub g_values: Vec<f64>,
    /// The normalized coupling matrix + external Q.
    pub coupling: CouplingMatrix,
    /// The swept ideal response.
    pub sweep: Vec<SweepPoint>,
    /// The forbidden mask regions for the response plot.
    pub mask_bands: Vec<MaskBand>,
    /// The real spec-mask verdict.
    pub report: MaskReport,
    /// The dimensioned edge-coupled board layout.
    pub layout: Layout,
    /// Per-resonator realized geometry + electricals.
    pub resonators: Vec<ResonatorRow>,
    /// Single-line effective permittivity at the synthesized width.
    pub line_eps_eff: f64,
    /// Board bounding box, mm (`(width, height)`).
    pub board_size_mm: (f64, f64),
}

impl Designed {
    /// The filter order `N`.
    pub fn order(&self) -> usize {
        self.project.prototype.order()
    }
}

/// Run the full live engine pipeline on the hard-coded demo spec.
///
/// Mirrors `yee_cli::filter::run_synth`: synthesize → sweep the ideal response
/// → grade against the spec mask → dimension the edge-coupled board → derive
/// the per-resonator even/odd electricals from the solved gaps. Every value is
/// real engine output.
pub fn design_demo() -> Designed {
    let spec = demo_spec();
    let project = synthesize(&spec);
    let g_values = project.prototype.g.clone();
    let coupling = project.coupling.clone();

    // ---- swept ideal response (same sweep grid as the CLI) ----------------
    let freqs = sweep_freqs(spec.f0_hz, spec.fbw);
    let s21 = ideal_response(&project, &freqs);
    let sweep: Vec<SweepPoint> = freqs
        .iter()
        .zip(s21.iter())
        .map(|(&f, z)| {
            let s21_mag = z.norm().min(1.0);
            let s11_sq = (1.0 - s21_mag * s21_mag).max(0.0);
            SweepPoint {
                f_hz: f,
                s21_db: 20.0 * s21_mag.max(1e-12).log10(),
                s11_db: 10.0 * s11_sq.max(1e-12).log10(),
            }
        })
        .collect();

    let mask_bands = mask_bands(&spec);
    let report = check_mask(&project, &freqs);

    // ---- physical dimensioning (FR-4 substrate, edge-coupled) -------------
    let dims =
        dimension_edge_coupled(&project, &SUBSTRATE).expect("demo spec is realizable on FR-4");
    let layout = dimension_edge_coupled_layout(&project, &SUBSTRATE)
        .expect("demo spec layout is realizable on FR-4");

    let w_m = dims.line_width_m;
    let line_eps_eff = eps_eff(w_m, SUBSTRATE.height_m, SUBSTRATE.eps_r);

    let n = project.coupling.m.len();
    let resonators: Vec<ResonatorRow> = (0..n)
        .map(|i| {
            // Gaps + couplings exist for i in 0..n-1 (the i-th gap is to i+1).
            let (gap, z0e, z0o, eff_e, eff_o, rk, tk) = if i < n - 1 {
                let s = dims.gaps_m[i];
                let cm: CoupledMicrostrip =
                    coupled_microstrip(w_m, s, SUBSTRATE.height_m, SUBSTRATE.eps_r);
                let realized_k = (cm.z0e_ohm - cm.z0o_ohm) / (cm.z0e_ohm + cm.z0o_ohm);
                (
                    Some(s * 1e3),
                    Some(cm.z0e_ohm),
                    Some(cm.z0o_ohm),
                    Some(cm.eps_eff_e),
                    Some(cm.eps_eff_o),
                    Some(realized_k),
                    Some(dims.target_k[i]),
                )
            } else {
                (None, None, None, None, None, None, None)
            };
            ResonatorRow {
                id: i + 1,
                width_mm: w_m * 1e3,
                length_mm: dims.resonator_length_m * 1e3,
                gap_to_next_mm: gap,
                z0e_ohm: z0e,
                z0o_ohm: z0o,
                eps_eff_e: eff_e,
                eps_eff_o: eff_o,
                realized_k: rk,
                target_k: tk,
            }
        })
        .collect();

    let board_size_mm = board_size_mm(&layout.bbox);

    Designed {
        spec,
        project,
        g_values,
        coupling,
        sweep,
        mask_bands,
        report,
        layout,
        resonators,
        line_eps_eff,
        board_size_mm,
    }
}

/// Board size in mm from a layout bounding box.
fn board_size_mm(b: &BBox) -> (f64, f64) {
    (b.width() * 1e3, b.height() * 1e3)
}

/// Linear sweep of `SWEEP_POINTS` frequencies (mirrors the CLI).
fn sweep_freqs(f0: f64, fbw: f64) -> Vec<f64> {
    let half = SPAN_MULT * fbw / 2.0;
    let lo = (f0 * (1.0 - half)).max(f0 * 1e-3);
    let hi = f0 * (1.0 + half);
    (0..SWEEP_POINTS)
        .map(|i| lo + (hi - lo) * (i as f64) / ((SWEEP_POINTS - 1) as f64))
        .collect()
}

/// Mask forbidden regions for the plot (mirrors the CLI's `spec_mask_regions`):
/// a passband floor at `−passband_ripple_db` over `[f0·(1−fbw/2), f0·(1+fbw/2)]`
/// and a ceiling at `−reject` over a ±2 % band at each stopband point.
fn mask_bands(spec: &FilterSpec) -> Vec<MaskBand> {
    let f1 = spec.f0_hz * (1.0 - spec.fbw / 2.0);
    let f2 = spec.f0_hz * (1.0 + spec.fbw / 2.0);
    let mut bands = vec![MaskBand {
        f_lo_hz: f1,
        f_hi_hz: f2,
        is_floor: true,
        limit_db: -spec.mask.passband_ripple_db,
    }];
    for &(f_s, reject_db) in &spec.mask.stopband {
        bands.push(MaskBand {
            f_lo_hz: f_s * 0.98,
            f_hi_hz: f_s * 1.02,
            is_floor: false,
            limit_db: -reject_db,
        });
    }
    bands
}

// ===========================================================================
// Lumped-LC adapter (App.D.1L — F2.0 / F2.1 / F2.4 / F2.2)
// ===========================================================================

/// Fixed Monte-Carlo seed so the rendered yield is reproducible across reloads.
const YIELD_SEED: u64 = 0x59_45_45_5f_4c_43_30; // "YEE_LC0"
/// Number of Monte-Carlo trials for the tolerance/yield card (≈500 per ADR-0120).
const YIELD_TRIALS: usize = 500;
/// SMD footprint family used by the lumped board (the F2.2 default).
const LUMPED_FOOTPRINT: Footprint = Footprint::Smd0603;

/// One rendered resonator row of the lumped LC ladder (Synthesis stage).
#[derive(Clone, Copy)]
pub struct LumpedResonatorRow {
    /// Resonator index, 1-based (ladder order).
    pub index: usize,
    /// `true` if a series-arm resonator, `false` if a shunt-arm resonator.
    pub is_series: bool,
    /// Resonator inductance, nanohenries.
    pub l_nh: f64,
    /// Resonator capacitance, picofarads.
    pub c_pf: f64,
}

/// One rendered BOM line (Components + BOM stage), pre-formatted for display.
#[derive(Clone)]
pub struct BomRow {
    /// Reference designator class, e.g. `"L"` / `"C"`.
    pub ref_kind: &'static str,
    /// Whether this is an inductor or a capacitor.
    pub is_inductor: bool,
    /// Pretty ideal value with unit (e.g. `"2.34 nH"`).
    pub ideal_disp: String,
    /// Pretty chosen E-series value with unit.
    pub chosen_disp: String,
    /// Signed deviation of chosen from ideal, percent.
    pub deviation_pct: f64,
    /// Series tolerance, percent (±).
    pub tolerance_pct: f64,
    /// Quantity of this grouped part.
    pub qty: usize,
}

/// One rendered placement row (Layout stage), pre-formatted for display.
#[derive(Clone)]
pub struct PlacementRow {
    /// Reference designator (e.g. `"L1"`).
    pub ref_des: String,
    /// Footprint name (e.g. `"0603"`).
    pub footprint: &'static str,
    /// Branch role (`"series"` / `"shunt"`).
    pub kind: &'static str,
    /// Board-frame centre `x`, mm.
    pub cx_mm: f64,
    /// Board-frame centre `y`, mm.
    pub cy_mm: f64,
}

/// A BOM plus its derived display rows + summary, for one E-series.
pub struct BomView {
    /// Pretty series name (`"E24"` / `"E96"`).
    pub series_name: &'static str,
    /// Per-part tolerance, percent (the series tolerance).
    pub tolerance_pct: f64,
    /// Display-ready BOM rows (grouped, in first-encountered order).
    pub rows: Vec<BomRow>,
    /// Total physical part count (sum of all `qty`).
    pub total_parts: usize,
    /// Worst-case (largest magnitude) deviation across the BOM, percent.
    pub worst_deviation_pct: f64,
}

/// A yield result for one E-series, pre-extracted for the tolerance card.
#[derive(Clone, Copy, PartialEq)]
pub struct YieldView {
    /// Pretty series name (`"E24"` / `"E96"`).
    pub series_name: &'static str,
    /// Per-part tolerance, percent (±).
    pub tolerance_pct: f64,
    /// Yield as a percentage in `[0, 100]`.
    pub yield_pct: f64,
    /// Worst-case in-band return loss across all trials, dB.
    pub worst_rl_db: f64,
    /// Worst-case stopband rejection across all trials, dB.
    pub worst_rej_db: f64,
}

/// Everything the four lumped stages render — all from the live F2.x engine.
pub struct LumpedDesigned {
    /// The synthesized ideal LC ladder.
    pub ladder: LumpedLadder,
    /// Display-ready resonator rows.
    pub resonators: Vec<LumpedResonatorRow>,
    /// The swept ideal `ladder_s21` response (reuses [`SweepPoint`]).
    pub sweep: Vec<SweepPoint>,
    /// Forbidden mask regions for the response plot (shared with the
    /// distributed flow).
    pub mask_bands: Vec<MaskBand>,
    /// The realized-response spec-mask verdict on the ideal ladder.
    pub verdict: MaskVerdict,
    /// The E24 BOM view.
    pub bom_e24: BomView,
    /// The E96 BOM view.
    pub bom_e96: BomView,
    /// The E24 Monte-Carlo yield.
    pub yield_e24: YieldView,
    /// The E96 Monte-Carlo yield.
    pub yield_e96: YieldView,
    /// The placed lumped board (geometry + placements).
    pub board: LumpedBoard,
    /// Display-ready placement rows.
    pub placements: Vec<PlacementRow>,
    /// Board bounding box, mm (`(width, height)`).
    pub board_size_mm: (f64, f64),
    /// Number of trials the yield was run over.
    pub yield_trials: usize,
}

impl LumpedDesigned {
    /// The ladder order `N` (number of resonators).
    pub fn order(&self) -> usize {
        self.ladder.resonators.len()
    }
}

/// Pretty-print an inductance (henries) with an engineering unit.
fn fmt_henry(l: f64) -> String {
    if l >= 1e-6 {
        format!("{:.3} µH", l * 1e6)
    } else if l >= 1e-9 {
        format!("{:.3} nH", l * 1e9)
    } else {
        format!("{:.3} pH", l * 1e12)
    }
}

/// Pretty-print a capacitance (farads) with an engineering unit.
fn fmt_farad(c: f64) -> String {
    if c >= 1e-9 {
        format!("{:.3} nF", c * 1e9)
    } else if c >= 1e-12 {
        format!("{:.3} pF", c * 1e12)
    } else {
        format!("{:.3} fF", c * 1e15)
    }
}

/// Build a [`BomView`] for one E-series from the ladder.
fn bom_view(ladder: &LumpedLadder, series: ESeries) -> BomView {
    let bom: Bom = select_components(ladder, series);
    let series_name = match series {
        ESeries::E24 => "E24",
        ESeries::E96 => "E96",
    };
    let mut worst_deviation_pct: f64 = 0.0;
    let rows: Vec<BomRow> = bom
        .lines
        .iter()
        .map(|l| {
            worst_deviation_pct = worst_deviation_pct.max(l.deviation_pct.abs());
            let (ref_kind, is_inductor, ideal_disp, chosen_disp) = match l.kind {
                CompKind::Inductor => (
                    "L",
                    true,
                    fmt_henry(l.ideal_value),
                    fmt_henry(l.chosen_value),
                ),
                CompKind::Capacitor => (
                    "C",
                    false,
                    fmt_farad(l.ideal_value),
                    fmt_farad(l.chosen_value),
                ),
            };
            BomRow {
                ref_kind,
                is_inductor,
                ideal_disp,
                chosen_disp,
                deviation_pct: l.deviation_pct,
                tolerance_pct: l.tolerance_pct,
                qty: l.qty,
            }
        })
        .collect();
    BomView {
        series_name,
        tolerance_pct: series.tolerance_pct(),
        rows,
        total_parts: bom.total_parts(),
        worst_deviation_pct,
    }
}

/// Run the Monte-Carlo yield for one E-series and pack it for display.
fn yield_view(ladder: &LumpedLadder, series: ESeries, mask: &SpecMask) -> YieldView {
    let r = monte_carlo_yield(ladder, series, mask, YIELD_TRIALS, YIELD_SEED);
    let series_name = match series {
        ESeries::E24 => "E24",
        ESeries::E96 => "E96",
    };
    YieldView {
        series_name,
        tolerance_pct: series.tolerance_pct(),
        yield_pct: r.yield_fraction * 100.0,
        worst_rl_db: r.worst_inband_rl_db,
        worst_rej_db: r.worst_stopband_rej_db,
    }
}

/// Run the full live lumped-LC pipeline on the same demo spec the distributed
/// flow uses: synthesize → LC ladder → swept realized response → spec-mask
/// verdict → E24/E96 component selection + BOM → Monte-Carlo yield → SMD board
/// placement. Every value is real F2.x engine output.
pub fn design_lumped() -> LumpedDesigned {
    let spec = demo_spec();
    let project = synthesize(&spec);
    let ladder = synthesize_lumped(&project).expect("demo spec is a realizable band-pass ladder");

    // ---- display rows -----------------------------------------------------
    let resonators: Vec<LumpedResonatorRow> = ladder
        .resonators
        .iter()
        .enumerate()
        .map(|(i, r)| LumpedResonatorRow {
            index: i + 1,
            is_series: r.branch == LcBranch::Series,
            l_nh: r.l_henry * 1e9,
            c_pf: r.c_farad * 1e12,
        })
        .collect();

    // ---- swept realized (ABCD) response + mask verdict --------------------
    let freqs = sweep_freqs(spec.f0_hz, spec.fbw);
    let sweep: Vec<SweepPoint> = freqs
        .iter()
        .map(|&f| {
            let s21_mag = ladder_s21(&ladder, f).norm().min(1.0);
            let s11_sq = (1.0 - s21_mag * s21_mag).max(0.0);
            SweepPoint {
                f_hz: f,
                s21_db: 20.0 * s21_mag.max(1e-12).log10(),
                s11_db: 10.0 * s11_sq.max(1e-12).log10(),
            }
        })
        .collect();
    let mask_bands = mask_bands(&spec);
    let verdict = mask_verdict(&ladder, &spec.mask, spec.f0_hz, spec.fbw, &freqs, 0.0);

    // ---- E24 / E96 component selection + BOM ------------------------------
    let bom_e24 = bom_view(&ladder, ESeries::E24);
    let bom_e96 = bom_view(&ladder, ESeries::E96);

    // ---- Monte-Carlo tolerance / yield ------------------------------------
    let yield_e24 = yield_view(&ladder, ESeries::E24, &spec.mask);
    let yield_e96 = yield_view(&ladder, ESeries::E96, &spec.mask);

    // ---- SMD board placement ----------------------------------------------
    let board = lumped_board(&ladder, &SUBSTRATE, LUMPED_FOOTPRINT);
    let placements: Vec<PlacementRow> = board
        .placements
        .iter()
        .map(|p| PlacementRow {
            ref_des: p.ref_des.clone(),
            footprint: footprint_name(p.footprint),
            kind: match p.kind {
                BranchKind::Series => "series",
                BranchKind::Shunt => "shunt",
            },
            cx_mm: p.center_m.0 * 1e3,
            cy_mm: p.center_m.1 * 1e3,
        })
        .collect();
    let board_size_mm = board_size_mm(&board.layout.bbox);

    LumpedDesigned {
        ladder,
        resonators,
        sweep,
        mask_bands,
        verdict,
        bom_e24,
        bom_e96,
        yield_e24,
        yield_e96,
        board,
        placements,
        board_size_mm,
        yield_trials: YIELD_TRIALS,
    }
}

/// The `lumped_board` requires a `Footprint` import via `yee_filter`; this is
/// the human-readable family name for the placement table.
fn footprint_name(fp: Footprint) -> &'static str {
    match fp {
        Footprint::Smd0402 => "0402",
        Footprint::Smd0603 => "0603",
        Footprint::Smd0805 => "0805",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn demo_design_is_real_and_passes() {
        let d = design_demo();
        assert_eq!(d.order(), 5, "cheb_bpf fixture is order 5");
        // Synchronous coupling matrix is N×N with a zero diagonal.
        assert_eq!(d.coupling.m.len(), 5);
        assert!(d.coupling.m[0][0].abs() < 1e-12, "zero diagonal");
        // Real off-diagonal couplings are positive and symmetric.
        assert!(d.coupling.m[0][1] > 0.0);
        assert_eq!(d.resonators.len(), 5);
        // The fixture is designed to PASS its mask.
        assert!(d.report.pass, "cheb_bpf fixture should PASS its spec mask");
        // Board has positive extent.
        assert!(d.board_size_mm.0 > 0.0 && d.board_size_mm.1 > 0.0);
        // Realized k tracks the target within the bisection tolerance.
        for r in &d.resonators {
            if let (Some(rk), Some(tk)) = (r.realized_k, r.target_k) {
                assert!((rk - tk).abs() / tk < 1e-2, "realized k near target");
            }
        }
    }

    #[test]
    fn lumped_design_is_real() {
        let d = design_lumped();
        // Order-5 ladder → 5 resonators, shunt-first alternating.
        assert_eq!(d.order(), 5);
        assert_eq!(d.resonators.len(), 5);
        assert!(d.resonators[0].l_nh > 0.0 && d.resonators[0].c_pf > 0.0);
        assert!(!d.resonators[0].is_series, "shunt-first");
        assert!(d.resonators[1].is_series, "second is series");
        // Swept response is populated and finite.
        assert_eq!(d.sweep.len(), SWEEP_POINTS);
        assert!(d.sweep.iter().all(|s| s.s21_db.is_finite()));
        // Each resonator emits an L and a C → 10 BOM parts (before grouping
        // they total 10; grouping may collapse symmetric duplicates).
        assert_eq!(d.bom_e24.total_parts, 10);
        assert_eq!(d.bom_e96.total_parts, 10);
        assert!(!d.bom_e24.rows.is_empty() && !d.bom_e96.rows.is_empty());
        // E96 deviation is tighter than E24 (finer grid).
        assert!(d.bom_e96.worst_deviation_pct <= d.bom_e24.worst_deviation_pct + 1e-9);
        // Yield is a percentage; same seed → reproducible.
        assert!((0.0..=100.0).contains(&d.yield_e24.yield_pct));
        assert!((0.0..=100.0).contains(&d.yield_e96.yield_pct));
        let again = design_lumped();
        assert_eq!(d.yield_e24.yield_pct, again.yield_e24.yield_pct);
        // Board placed every component (2 per resonator).
        assert_eq!(d.placements.len(), 10);
        assert!(d.board_size_mm.0 > 0.0 && d.board_size_mm.1 > 0.0);
    }
}
