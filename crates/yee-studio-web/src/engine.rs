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
    Approximation, CouplingMatrix, FilterProject, FilterSpec, MaskReport, Response, SpecMask,
    check_mask, dimension_edge_coupled, dimension_edge_coupled_layout, ideal_response, synthesize,
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
}
