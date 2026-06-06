//! Engine bridge — drives the live Yee filter engine for the studio.
//!
//! Seeds from the committed `cheb_bpf.toml` fixture as a
//! [`yee_filter::FilterSpec`] (no file IO in WASM); the Spec stage then edits a
//! live spec that re-drives the **real** synthesis / dimensioning /
//! response-sweep / spec-mask / lumped paths through `yee-synth`, `yee-filter`,
//! and `yee-layout` — the same calls the `yee filter synth` CLI makes, mirrored
//! here. Nothing is faked: every number the stages render comes out of
//! [`Designed`] / [`LumpedDesigned`].

use num_complex::Complex64;
use yee_filter::{
    Approximation, BoardTopology, Bom, BranchKind, CompKind, CouplingMatrix, ESeries,
    FilterProject, FilterSpec, Footprint, LcBranch, LumpedBoard, LumpedLadder, MaskReport,
    MaskVerdict, OrderableBoard, RealizationTechnique, Response, SpecMask, SteppedSection,
    check_mask, dimension_combline, dimension_combline_layout, dimension_edge_coupled,
    dimension_edge_coupled_layout, dimension_hairpin, dimension_hairpin_layout,
    dimension_interdigital, dimension_interdigital_layout, dimension_stepped_impedance,
    dimension_stepped_impedance_layout, ideal_response, ideal_response_lowpass, jlcpcb_bom_csv,
    jlcpcb_cpl_csv, ladder_s_params_lossy, ladder_s21, ladder_s21_lossy, lumped_board,
    mask_verdict, monte_carlo_yield, select_components, synthesize, synthesize_lumped,
    synthesize_orderable_on,
};
use yee_io::touchstone::{File, Format, FreqUnit};
use yee_layout::{BBox, CoupledMicrostrip, Layout, Substrate, coupled_microstrip, eps_eff};

use crate::stages::Topology;

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
    /// Resonator length, mm. For edge-coupled this is the straight half-wave
    /// (`λ_g/2`); for hairpin it is the single U-folded arm length (`≈ λ_g/4`).
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
    /// The distributed realization the geometry was dimensioned for
    /// (edge-coupled, hairpin, combline, or interdigital). Drives the topology-aware Layout / Export
    /// labels; the synthesis / response / verdict fields are independent of it.
    pub topology: Topology,
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
    /// The dimensioned board layout for the active distributed [`topology`]
    /// (edge-coupled, hairpin, combline, or interdigital), or `None` when the coupling is not realizable
    /// on FR-4 (see [`dim_error`](Designed::dim_error)).
    pub layout: Option<Layout>,
    /// Per-resonator realized geometry + electricals (empty when geometry is
    /// not realizable).
    pub resonators: Vec<ResonatorRow>,
    /// Single-line effective permittivity at the synthesized width (`0.0` when
    /// geometry is not realizable).
    pub line_eps_eff: f64,
    /// Board bounding box, mm (`(width, height)`); `(0, 0)` when not realizable.
    pub board_size_mm: (f64, f64),
    /// The distributed dimensioning error string (for the active [`topology`])
    /// when the coupling could not be realized on FR-4, else `None`. The
    /// synthesis / response / verdict stay real even when this is `Some`.
    pub dim_error: Option<String>,
    /// The combline loading capacitance `C_L = cot(θ0)/(2π·f0·Z0)`, farads — the
    /// combline-distinct quantity. A single value (uniform θ0 / Z0 → the same
    /// `C_L` at every resonator's open end). `Some` only for
    /// [`Topology::Combline`] with a realizable geometry; `None` for every other
    /// topology and when the combline geometry is not realizable.
    pub combline_loading_cap_f: Option<f64>,
}

impl Designed {
    /// The filter order `N`.
    pub fn order(&self) -> usize {
        self.project.prototype.order()
    }

    /// Human-readable name of the realized distributed topology
    /// (e.g. `"edge-coupled ½λ"` / `"hairpin (U-folded ½λ)"`), for the Layout /
    /// Export headings.
    pub fn topology_name(&self) -> &'static str {
        match self.topology {
            Topology::Hairpin => "hairpin (U-folded ½λ)",
            Topology::Combline => "combline (capacitively-loaded)",
            Topology::Interdigital => "interdigital (λ/4, alt. short)",
            // The lumped + stepped-impedance flows render their own Layout /
            // Export (the band-pass `Designed` is never built for them — they use
            // `LumpedDesigned` / `SteppedLowpassDesigned`); the carried distributed
            // `Designed` is dimensioned edge-coupled.
            Topology::EdgeCoupled | Topology::LumpedLc | Topology::SteppedImpedance => {
                "edge-coupled ½λ"
            }
        }
    }

    /// The resonator-table length-column label for the realized topology: the
    /// edge-coupled straight `λ_g/2` resonator vs the hairpin U-folded `λ_g/4`
    /// arm.
    pub fn length_label(&self) -> &'static str {
        match self.topology {
            Topology::Hairpin => "arm length (mm)",
            Topology::Combline => "resonator length (mm)",
            Topology::Interdigital => "resonator length (mm)",
            Topology::EdgeCoupled | Topology::LumpedLc | Topology::SteppedImpedance => {
                "length (mm)"
            }
        }
    }
}

/// Run the full live engine pipeline on the hard-coded demo spec, realized as an
/// edge-coupled board.
///
/// Convenience wrapper around [`design_demo_from`] for the initial boot state
/// (the Spec stage edits a live spec that re-drives [`design_demo_from`]).
pub fn design_demo() -> Designed {
    design_demo_from(demo_spec(), Topology::EdgeCoupled)
}

/// Run the full live engine pipeline on an arbitrary [`yee_filter::FilterSpec`],
/// realized for the given distributed [`Topology`].
///
/// Mirrors `yee_cli::filter::run_synth`: synthesize → sweep the ideal response
/// → grade against the spec mask → dimension the board for `topology` → derive
/// the per-resonator even/odd electricals from the solved gaps. Every value is
/// real engine output.
///
/// The synthesis / response / mask verdict are **topology-independent**
/// (edge-coupled, hairpin, combline, and interdigital all realize the *same* coupled-resonator
/// band-pass prototype); only the geometry-derived fields differ. `topology`
/// selects the dimensioner: [`Topology::EdgeCoupled`] →
/// [`dimension_edge_coupled`], [`Topology::Hairpin`] → [`dimension_hairpin`].
/// [`Topology::LumpedLc`] has no distributed board, so it falls back to the
/// edge-coupled geometry (the lumped flow uses [`design_lumped_from`] for its
/// own board).
///
/// The distributed dimensioning can fail when the requested coupling is not
/// realizable on FR-4 (e.g. an over-wide bandwidth at a low order); in that
/// case the geometry-derived fields fall back to empty (no resonator rows, a
/// zero board size) while the synthesized prototype / coupling / response /
/// mask verdict remain real. The Spec form surfaces that as a "geometry not
/// realizable" note rather than panicking the whole app.
pub fn design_demo_from(spec: FilterSpec, topology: Topology) -> Designed {
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

    // ---- physical dimensioning (FR-4 substrate, per topology) -------------
    // Fallible: an unrealizable coupling (e.g. a wide FBW at a low order from
    // the Spec form) returns a `DimError`. We keep the synthesized prototype /
    // response / verdict real and degrade only the geometry-derived fields,
    // surfacing the error through `dim_error` instead of panicking the app.
    let geom = derive_geometry(&project, topology);

    Designed {
        topology,
        spec,
        project,
        g_values,
        coupling,
        sweep,
        mask_bands,
        report,
        layout: geom.layout,
        resonators: geom.resonators,
        line_eps_eff: geom.line_eps_eff,
        board_size_mm: geom.board_size_mm,
        dim_error: geom.dim_error,
        combline_loading_cap_f: geom.combline_loading_cap_f,
    }
}

/// The geometry-derived fields of [`Designed`], bundled so the fallible
/// dimensioning path returns one value (and the unrealizable case degrades
/// every field coherently).
struct Geometry {
    layout: Option<Layout>,
    resonators: Vec<ResonatorRow>,
    line_eps_eff: f64,
    board_size_mm: (f64, f64),
    dim_error: Option<String>,
    /// The combline loading cap `C_L`, farads — `Some` only for a realizable
    /// [`Topology::Combline`] geometry, `None` otherwise.
    combline_loading_cap_f: Option<f64>,
}

/// The topology-specific output of a distributed dimensioner, normalized so the
/// shared per-resonator electrical table + layout packaging in
/// [`derive_geometry`] can treat edge-coupled, hairpin, combline, and
/// interdigital uniformly. (All invert the same coupled-microstrip model; only
/// the line width / resonator length / assembled [`Layout`] / loading cap
/// differ.)
struct SolvedDistributed {
    /// Strip / arm width for the spec `Z0`, metres.
    line_width_m: f64,
    /// Per-resonator length, metres: edge-coupled straight `λ_g/2`, or the
    /// hairpin U-folded arm `≈ λ_g/4`.
    resonator_length_m: f64,
    /// Inter-resonator edge-coupling gaps, metres (length `N − 1`).
    gaps_m: Vec<f64>,
    /// The `FBW · m_{i,i+1}` coupling each gap was solved for (length `N − 1`).
    target_k: Vec<f64>,
    /// The assembled board layout.
    layout: Layout,
    /// The combline loading cap `C_L`, farads — `Some` only for the combline
    /// dimensioner, `None` for edge-coupled / hairpin (which have no loading cap).
    loading_cap_f: Option<f64>,
}

/// Dimension the synthesized project onto FR-4 for the given [`Topology`],
/// returning the geometry fields or — when the coupling is not realizable — an
/// empty geometry carrying the error string.
///
/// Both distributed topologies invert the **same** validated coupled-microstrip
/// model from the same coupling matrix; they differ only in the realized
/// resonator geometry (line width, resonator length, the assembled board
/// [`Layout`]). The per-resonator even/odd electrical table is therefore shared
/// — it is recovered from the solved inter-resonator gaps at the common line
/// width — and only the line width / resonator length / `Layout` vary.
fn derive_geometry(project: &FilterProject, topology: Topology) -> Geometry {
    // Solve the topology-specific dimensions: the line width, the per-resonator
    // length, the inter-resonator gaps + target couplings, and the assembled
    // board `Layout`. `LumpedLc` has no distributed board; it reuses the
    // edge-coupled geometry (its own board comes from `design_lumped_from`).
    let solved: Result<SolvedDistributed, String> = match topology {
        Topology::Hairpin => {
            match (
                dimension_hairpin(project, &SUBSTRATE),
                dimension_hairpin_layout(project, &SUBSTRATE),
            ) {
                (Ok(dims), Ok(layout)) => Ok(SolvedDistributed {
                    line_width_m: dims.line_width_m,
                    resonator_length_m: dims.arm_length_m,
                    gaps_m: dims.gaps_m,
                    target_k: dims.target_k,
                    layout,
                    loading_cap_f: None,
                }),
                (dims_res, layout_res) => Err(dims_res
                    .err()
                    .map(|e| e.to_string())
                    .or_else(|| layout_res.err().map(|e| e.to_string()))
                    .unwrap_or_else(|| "hairpin geometry is not realizable on FR-4".into())),
            }
        }
        // Combline: capacitively-loaded short-circuited θ0 resonators with a
        // loading cap C_L at each open end. θ0 = π/4 = λg/8 is the compact
        // default (ADR-0146). Same coupled-microstrip gap inversion as
        // edge-coupled / hairpin; the distinct quantities are the short θ0
        // resonator length and the surfaced loading cap.
        Topology::Combline => {
            let theta0 = std::f64::consts::FRAC_PI_4;
            match (
                dimension_combline(project, theta0, &SUBSTRATE),
                dimension_combline_layout(project, theta0, &SUBSTRATE),
            ) {
                (Ok(dims), Ok(layout)) => Ok(SolvedDistributed {
                    line_width_m: dims.line_width_m,
                    resonator_length_m: dims.resonator_length_m,
                    gaps_m: dims.gaps_m,
                    target_k: dims.target_k,
                    layout,
                    loading_cap_f: Some(dims.loading_cap_f),
                }),
                (dims_res, layout_res) => Err(dims_res
                    .err()
                    .map(|e| e.to_string())
                    .or_else(|| layout_res.err().map(|e| e.to_string()))
                    .unwrap_or_else(|| "combline geometry is not realizable on FR-4".into())),
            }
        }
        // Interdigital: full λg/4 lines short-circuited at alternating ends, with
        // NO loading cap (the θ = π/2 self-resonant limit of combline). Same
        // coupled-microstrip gap inversion as edge-coupled / hairpin / combline;
        // the distinct quantity is the full λg/4 resonator length (no θ0 param,
        // no cap → `loading_cap_f: None`).
        Topology::Interdigital => {
            match (
                dimension_interdigital(project, &SUBSTRATE),
                dimension_interdigital_layout(project, &SUBSTRATE),
            ) {
                (Ok(dims), Ok(layout)) => Ok(SolvedDistributed {
                    line_width_m: dims.line_width_m,
                    resonator_length_m: dims.resonator_length_m,
                    gaps_m: dims.gaps_m,
                    target_k: dims.target_k,
                    layout,
                    loading_cap_f: None,
                }),
                (dims_res, layout_res) => Err(dims_res
                    .err()
                    .map(|e| e.to_string())
                    .or_else(|| layout_res.err().map(|e| e.to_string()))
                    .unwrap_or_else(|| "interdigital geometry is not realizable on FR-4".into())),
            }
        }
        // `SteppedImpedance` has its own low-pass `SteppedLowpassDesigned` flow
        // and never builds a band-pass `Designed`; if it ever reached here it
        // would harmlessly reuse the edge-coupled geometry.
        Topology::EdgeCoupled | Topology::LumpedLc | Topology::SteppedImpedance => {
            match (
                dimension_edge_coupled(project, &SUBSTRATE),
                dimension_edge_coupled_layout(project, &SUBSTRATE),
            ) {
                (Ok(dims), Ok(layout)) => Ok(SolvedDistributed {
                    line_width_m: dims.line_width_m,
                    resonator_length_m: dims.resonator_length_m,
                    gaps_m: dims.gaps_m,
                    target_k: dims.target_k,
                    layout,
                    loading_cap_f: None,
                }),
                (dims_res, layout_res) => Err(dims_res
                    .err()
                    .map(|e| e.to_string())
                    .or_else(|| layout_res.err().map(|e| e.to_string()))
                    .unwrap_or_else(|| "edge-coupled geometry is not realizable on FR-4".into())),
            }
        }
    };

    match solved {
        Ok(SolvedDistributed {
            line_width_m: w_m,
            resonator_length_m: length_m,
            gaps_m,
            target_k,
            layout,
            loading_cap_f,
        }) => {
            let line_eps_eff = eps_eff(w_m, SUBSTRATE.height_m, SUBSTRATE.eps_r);

            let n = project.coupling.m.len();
            let resonators: Vec<ResonatorRow> = (0..n)
                .map(|i| {
                    // Gaps + couplings exist for i in 0..n-1 (gap i is to i+1).
                    let (gap, z0e, z0o, eff_e, eff_o, rk, tk) = if i < n - 1 {
                        let s = gaps_m[i];
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
                            Some(target_k[i]),
                        )
                    } else {
                        (None, None, None, None, None, None, None)
                    };
                    ResonatorRow {
                        id: i + 1,
                        width_mm: w_m * 1e3,
                        length_mm: length_m * 1e3,
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
            Geometry {
                layout: Some(layout),
                resonators,
                line_eps_eff,
                board_size_mm,
                dim_error: None,
                combline_loading_cap_f: loading_cap_f,
            }
        }
        Err(msg) => Geometry {
            layout: None,
            resonators: Vec::new(),
            line_eps_eff: 0.0,
            board_size_mm: (0.0, 0.0),
            dim_error: Some(msg),
            combline_loading_cap_f: None,
        },
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
// Touchstone .s2p export (ADR-0171, T8)
// ===========================================================================

/// Render a 2-port reciprocal, symmetric filter's `(S11, S21)` sweep to a
/// Touchstone `.s2p` string (Real/Imag, Hz), for a client-side download.
///
/// Mirrors `yee_cli::filter::write_s2p`'s [`File`] construction exactly — the
/// same studio↔CLI `.s2p` contract — but renders to a **string** via the
/// WASM-safe [`yee_io::touchstone::to_string`] instead of writing to disk
/// (there is no filesystem in the browser). The 2-port matrix is
/// `[[S11, S21], [S21, S11]]` (reciprocal `S12 = S21`, symmetric `S22 = S11`),
/// stored row-major as `[S11, S21, S21, S11]`.
///
/// The **caller owns the physics**: the supplied `(S11, S21)` pairs are written
/// verbatim and S11 is never re-derived. The pair must be passive
/// (`|S11|² + |S21|² ≤ 1`) for the file to re-parse via [`yee_io::touchstone`] —
/// the finite-Q lumped pair ([`ladder_s_params_lossy`]) is the true lossy
/// 2-port, and the ideal pair uses the lossless quadrature
/// (`S11 = j·√(1−|S21|²)`, [`lossless_s_pair`]); both are passive.
///
/// Returns the renderer's [`yee_io::Error`] display string on failure (a
/// length mismatch, a non-finite value, etc.); the Export buttons treat an
/// `Err` as a no-op (no broken download).
fn s2p_string(
    z0: f64,
    freqs: &[f64],
    s_params: &[(Complex64, Complex64)],
    comment: &str,
) -> Result<String, String> {
    let mut data = Vec::with_capacity(freqs.len());
    for &(s11, s21) in s_params {
        // Row-major n×n: [S11, S12, S21, S22]; reciprocal + symmetric.
        data.push(vec![s11, s21, s21, s11]);
    }
    let file = File {
        n_ports: 2,
        z0,
        freq_unit: FreqUnit::Hz,
        format: Format::RealImag,
        freq_hz: freqs.to_vec(),
        data,
        comments: vec![comment.to_string()],
    };
    yee_io::touchstone::to_string(&file).map_err(|e| e.to_string())
}

/// Map a lossless transfer `S21` to the `(S11, S21)` pair written for the
/// distributed / ideal `.s2p` path. The transmission is taken with `|S21| ≤ 1`
/// and S11 is the **true lossless reflection** placed in quadrature
/// (`S11 = j·√(1−|S21|²)`): for a lossless reciprocal symmetric 2-port this is
/// the exact reflection — it makes `|S11|² + |S21|² = 1` and the S-matrix
/// passive, so the rendered file re-parses. Mirrors
/// `yee_cli::filter::lossless_s_pair`.
fn lossless_s_pair(s21: Complex64) -> (Complex64, Complex64) {
    let s21_mag = s21.norm().min(1.0);
    let s11_mag = (1.0 - s21_mag * s21_mag).max(0.0).sqrt();
    (Complex64::new(0.0, s11_mag), Complex64::new(s21_mag, 0.0))
}

/// Build the **lumped finite-Q** `(S11, S21)` sweep written to the lumped
/// `.s2p`: the realistic [`ladder_s_params_lossy`] response over the design
/// sweep grid (`sweep_freqs(ladder.f0_hz, ladder.fbw)`), at per-resonator
/// unloaded Q `q_unloaded`. This is the **same model + grid** as the displayed
/// [`LumpedDesigned::sweep_finite_q`] (`ladder_s21_lossy` is just this pair's
/// `.1`), so the exported `.s2p` matches the plotted finite-Q curve. The pair is
/// the true lossy 2-port (`|S11|² + |S21|² < 1`), hence passive.
fn lumped_s2p_sweep(
    ladder: &LumpedLadder,
    q_unloaded: f64,
) -> (Vec<f64>, Vec<(Complex64, Complex64)>) {
    let freqs = sweep_freqs(ladder.f0_hz, ladder.fbw);
    let s_params: Vec<(Complex64, Complex64)> = freqs
        .iter()
        .map(|&f| ladder_s_params_lossy(ladder, f, q_unloaded))
        .collect();
    (freqs, s_params)
}

/// Build the **distributed / ideal** `(S11, S21)` sweep written to the
/// distributed `.s2p`: the closed-form [`ideal_response`] `S21` over the design
/// sweep grid, with S11 placed in lossless quadrature via [`lossless_s_pair`].
/// This is the **same model + grid** as the displayed [`Designed::sweep`], so
/// the exported `.s2p` matches the plotted ideal curve, and the pair is passive
/// (`|S11|² + |S21|² = 1`).
fn distributed_s2p_sweep(
    project: &FilterProject,
    spec: &FilterSpec,
) -> (Vec<f64>, Vec<(Complex64, Complex64)>) {
    let freqs = sweep_freqs(spec.f0_hz, spec.fbw);
    let s_params: Vec<(Complex64, Complex64)> = ideal_response(project, &freqs)
        .into_iter()
        .map(lossless_s_pair)
        .collect();
    (freqs, s_params)
}

/// Render the lumped finite-Q `.s2p` string for `d` (ADR-0171). The Export
/// stage's lumped `.s2p` button calls this; an `Err` (a renderer failure) is a
/// no-op there (no broken download). The comment self-describes the response.
pub fn lumped_s2p(d: &LumpedDesigned) -> Result<String, String> {
    let (freqs, s_params) = lumped_s2p_sweep(&d.ladder, d.q_unloaded);
    let comment = format!(
        " yee studio — finite-Q lumped response (Q_unloaded = {})",
        d.q_unloaded
    );
    s2p_string(d.ladder.z0_ohm, &freqs, &s_params, &comment)
}

/// Render the distributed / ideal `.s2p` string for `d` (ADR-0171). The Export
/// stage's distributed `.s2p` button calls this; an `Err` is a no-op there.
pub fn distributed_s2p(d: &Designed) -> Result<String, String> {
    let (freqs, s_params) = distributed_s2p_sweep(&d.project, &d.spec);
    s2p_string(
        d.spec.z0_ohm,
        &freqs,
        &s_params,
        " yee studio — ideal closed-form response",
    )
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

/// Default per-resonator unloaded quality factor for the finite-Q realistic
/// response (ADR-0163). `Q_u = 100` is a representative surface-mount chip
/// component value (chip inductors `Q ≈ 30–100`); it produces a few-dB midband
/// insertion loss and rounded passband corners on the demo 3/5-pole filters.
pub const DEFAULT_Q_UNLOADED: f64 = 100.0;

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
    /// The swept ideal (lossless) `ladder_s21` response (reuses [`SweepPoint`]).
    pub sweep: Vec<SweepPoint>,
    /// The swept **realistic finite-Q** `ladder_s21_lossy` response on the same
    /// grid as [`sweep`](Self::sweep) (ADR-0163). Each resonator carries an
    /// unloaded Q of [`q_unloaded`](Self::q_unloaded); the curve shows the
    /// midband insertion loss + rounded corners a built filter actually exhibits.
    pub sweep_finite_q: Vec<SweepPoint>,
    /// The per-resonator unloaded quality factor the finite-Q sweep was computed
    /// at (the [`DEFAULT_Q_UNLOADED`] on boot; edited by the studio's `Q_u`
    /// control).
    pub q_unloaded: f64,
    /// The realized **midband insertion loss** `IL = −20·log10(|S21_lossy(f0)|)`,
    /// dB — the headline number a VNA measures (`≈ 0` in the lossless limit;
    /// matches Cohn's `4.343·Σg/(Q_u·FBW)` near band-centre, validated below).
    pub midband_il_db: f64,
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
    /// The auto-routed orderable JLCPCB assembly set (ADR-0169 T5).
    ///
    /// This is the **JLCPCB upload set** the Export stage offers: a
    /// `(topology × footprint)` search ([`orderable_upload`]) for an orderable
    /// realization, which may differ from the displayed [`board`](Self::board) /
    /// [`ladder`](Self::ladder) when the alternating ladder is not orderable for
    /// this spec (then the top-C-coupled board is auto-routed). The design view
    /// and finite-Q response above stay ladder-based — a documented
    /// walking-skeleton limit (routing the whole lumped flow incl. the finite-Q
    /// response to top-C is a follow-on); only the JLCPCB export auto-routes.
    pub orderable: OrderableExport,
}

impl LumpedDesigned {
    /// The ladder order `N` (number of resonators).
    pub fn order(&self) -> usize {
        self.ladder.resonators.len()
    }
}

/// The auto-routed **orderable JLCPCB assembly set** for the lumped design
/// (ADR-0169 T5).
///
/// Pure data (WASM-safe): the chosen [`BoardTopology`] + footprint, the
/// orderability verdict, the placed board (for Gerber/CPL geometry), and the
/// precomputed JLCPCB BOM/CPL CSV strings — everything the Export stage needs to
/// offer a self-consistent, honestly-bounded JLCPCB upload set, all from the
/// **same** auto-routed board.
///
/// Produced by [`orderable_upload`]. `fully_orderable` is reported exactly — when
/// `false`, `n_blank` of `n_parts` parts had no JLCPCB Basic match (the honest
/// "narrow-band lumped is distributed-only" case); the BOM/CPL still carry the
/// real blank rows, never a fabricated orderable board.
#[derive(Clone)]
pub struct OrderableExport {
    /// Human-readable chosen board topology (`"alternating ladder"` /
    /// `"top-C-coupled"`).
    pub topology_label: &'static str,
    /// Human-readable chosen SMD footprint family (`"0402"` / `"0603"` /
    /// `"0805"`).
    pub footprint_label: &'static str,
    /// `true` iff **every** part resolved to a real LCSC Basic part (zero blanks).
    pub fully_orderable: bool,
    /// Total number of placed parts on the auto-routed board.
    pub n_parts: usize,
    /// Number of parts with no JLCPCB Basic match (the honest coverage holes).
    pub n_blank: usize,
    /// The auto-routed placed board (renderable [`yee_layout::Layout`] +
    /// placements) the Gerber + CPL files are emitted from.
    pub board: LumpedBoard,
    /// The precomputed JLCPCB assembly-BOM CSV (`jlcpcb_bom_csv`).
    pub bom_csv: String,
    /// The precomputed JLCPCB component-placement (CPL) CSV (`jlcpcb_cpl_csv`).
    pub cpl_csv: String,
}

/// Footprints the auto-router searches, in preference order: **0402 first** (the
/// finest RF L/C value grid + smallest board — the RF-lumped norm), then 0603,
/// then 0805. Mirrors the ADR-0169 / ADR-0168 search order.
const ORDERABLE_FOOTPRINTS: [Footprint; 3] =
    [Footprint::Smd0402, Footprint::Smd0603, Footprint::Smd0805];

/// Search `(topology × footprint)` for an orderable JLCPCB realization of
/// `project` on `substrate`, and surface it honestly (ADR-0169 T5).
///
/// For each footprint in [`ORDERABLE_FOOTPRINTS`] (0402 → 0603 → 0805) this calls
/// [`synthesize_orderable_on`], which itself picks the orderable board *topology*
/// (alternating ladder vs top-C-coupled) for that footprint. The first footprint
/// whose result is `fully_orderable` wins; if none is fully orderable, the
/// footprint with the **fewest total blanks** is returned (the first/0402 on a
/// tie, since the search keeps the earliest minimum). The returned
/// [`OrderableExport`] reflects `fully_orderable` exactly — no fabricated
/// orderability.
///
/// The chosen board's parts feed [`jlcpcb_bom_csv`] and its placements feed
/// [`jlcpcb_cpl_csv`], so the BOM, CPL, and (in the Export stage) the Gerbers all
/// come from one self-consistent board.
///
/// # Robustness
///
/// [`synthesize_orderable_on`] only errors when the spec is not a realizable
/// band-pass of order `N ≥ 1`; a footprint that errors is skipped. If *every*
/// footprint errors (a spec the lumped track cannot realize at all), this falls
/// back to a single, honest **empty** export (`n_parts == 0`,
/// `fully_orderable == false`, empty CSVs) rather than panicking — the lumped
/// engine path that calls this has already turned a non-realizable spec into a
/// user-facing error before reaching here, so this branch is defensive only.
/// Pure-compute / WASM-safe (no I/O, network, threads, time, RNG, or `unsafe`).
pub fn orderable_upload(project: &FilterProject, substrate: &Substrate) -> OrderableExport {
    // Track the fewest-blanks candidate across footprints (kept for the
    // not-fully-orderable case). `<` (strict) keeps the EARLIEST minimum, so a
    // tie resolves to the first footprint tried (0402).
    let mut best: Option<(Footprint, OrderableBoard)> = None;
    for &fp in &ORDERABLE_FOOTPRINTS {
        let Ok(ob) = synthesize_orderable_on(project, substrate, fp) else {
            // A footprint that fails to synthesize is skipped; see `# Robustness`.
            continue;
        };
        if ob.fully_orderable {
            // First fully-orderable footprint wins outright.
            return build_export(fp, ob);
        }
        let n_blank = ob.parts.iter().filter(|p| p.lcsc.is_none()).count();
        let take = match &best {
            Some((_, prev)) => n_blank < prev.parts.iter().filter(|p| p.lcsc.is_none()).count(),
            None => true,
        };
        if take {
            best = Some((fp, ob));
        }
    }
    match best {
        Some((fp, ob)) => build_export(fp, ob),
        // Defensive: no footprint synthesized (a non-realizable spec — which the
        // lumped engine path rejects before reaching here). Return an honest
        // empty export (zero parts, not orderable, empty CSVs) rather than
        // panicking. The board is an empty placement on the spec's own f0/fbw/Z0.
        None => {
            let empty = LumpedLadder {
                f0_hz: project.spec.f0_hz,
                fbw: project.spec.fbw,
                z0_ohm: project.spec.z0_ohm,
                resonators: Vec::new(),
            };
            OrderableExport {
                topology_label: "alternating ladder",
                footprint_label: footprint_name(ORDERABLE_FOOTPRINTS[0]),
                fully_orderable: false,
                n_parts: 0,
                n_blank: 0,
                board: lumped_board(&empty, substrate, ORDERABLE_FOOTPRINTS[0]),
                bom_csv: String::new(),
                cpl_csv: String::new(),
            }
        }
    }
}

/// Build an [`OrderableExport`] from a chosen footprint + [`OrderableBoard`]:
/// derive the labels, count blanks, and precompute the JLCPCB BOM/CPL CSVs.
fn build_export(fp: Footprint, ob: OrderableBoard) -> OrderableExport {
    let topology_label = match ob.topology {
        BoardTopology::AlternatingLadder => "alternating ladder",
        BoardTopology::TopCCoupled => "top-C-coupled",
    };
    let n_parts = ob.parts.len();
    let n_blank = ob.parts.iter().filter(|p| p.lcsc.is_none()).count();
    let bom_csv = jlcpcb_bom_csv(&ob.parts);
    let cpl_csv = jlcpcb_cpl_csv(&ob.board.placements);
    OrderableExport {
        topology_label,
        footprint_label: footprint_name(fp),
        fully_orderable: ob.fully_orderable,
        n_parts,
        n_blank,
        board: ob.board,
        bom_csv,
        cpl_csv,
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

/// Run the full live lumped-LC pipeline on the hard-coded demo spec at the
/// [`DEFAULT_Q_UNLOADED`].
///
/// Convenience wrapper around [`design_lumped_from`] for the initial boot
/// state (the Spec stage edits a live spec, and the `Q_u` control edits the
/// unloaded Q, both re-driving the pipeline).
pub fn design_lumped() -> LumpedDesigned {
    design_lumped_from(demo_spec(), DEFAULT_Q_UNLOADED)
        .expect("demo spec is a realizable band-pass ladder")
}

/// Run the full live lumped-LC pipeline on an arbitrary
/// [`yee_filter::FilterSpec`] at per-resonator unloaded quality factor
/// `q_unloaded`: synthesize → LC ladder → swept ideal + **finite-Q realized**
/// response → spec-mask verdict → E24/E96 component selection + BOM →
/// Monte-Carlo yield → SMD board placement. Every value is real F2.x engine
/// output.
///
/// `q_unloaded` (ADR-0163) is the single lumped unloaded Q modelling each
/// resonator's dissipation; it drives the [`sweep_finite_q`](LumpedDesigned::sweep_finite_q)
/// curve + the realized [`midband_il_db`](LumpedDesigned::midband_il_db) via
/// [`ladder_s21_lossy`]. The ideal lossless [`sweep`](LumpedDesigned::sweep) and
/// the spec-mask [`verdict`](LumpedDesigned::verdict) are independent of it.
/// Pass [`DEFAULT_Q_UNLOADED`] for the boot default.
///
/// Returns the `LumpedError` display string when the prototype cannot be
/// mapped to a realizable band-pass ladder (e.g. a degenerate FBW); the Spec
/// form surfaces that rather than panicking.
pub fn design_lumped_from(spec: FilterSpec, q_unloaded: f64) -> Result<LumpedDesigned, String> {
    let project = synthesize(&spec);
    let ladder = synthesize_lumped(&project).map_err(|e| e.to_string())?;

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

    // ---- swept ideal (lossless ABCD) response + mask verdict --------------
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

    // ---- swept realistic finite-Q response (ADR-0163) ---------------------
    // The same ABCD cascade with per-resonator unloaded-Q loss, on the same
    // grid as `sweep`. This is dissipative (|S11|² + |S21|² < 1), so |S11| is
    // the TRUE absorptive reflection, not a lossless √(1−|S21|²) placeholder.
    let sweep_finite_q: Vec<SweepPoint> = freqs
        .iter()
        .map(|&f| {
            let s21_mag = ladder_s21_lossy(&ladder, f, q_unloaded).norm().min(1.0);
            SweepPoint {
                f_hz: f,
                s21_db: 20.0 * s21_mag.max(1e-12).log10(),
                // |S11| is not plotted for the finite-Q curve (the overlay draws
                // |S21| only); record the lossy |S21| dB as the carried value.
                s11_db: 20.0 * s21_mag.max(1e-12).log10(),
            }
        })
        .collect();
    // Realized midband insertion loss: IL = −20·log10(|S21_lossy(f0)|) at the
    // band-centre f0 (the headline VNA number). Read straight from the lossy
    // response so it is consistent with the plotted curve.
    let s21_f0 = ladder_s21_lossy(&ladder, spec.f0_hz, q_unloaded)
        .norm()
        .min(1.0);
    let midband_il_db = -20.0 * s21_f0.max(1e-12).log10();

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

    // ---- auto-routed orderable JLCPCB assembly set (ADR-0169 T5) ----------
    // Search (topology × footprint) for an orderable JLCPCB realization on the
    // same FR-4 SUBSTRATE the design board uses. This is the JLCPCB upload set
    // the Export stage offers; it may differ from the displayed ladder `board`
    // when the alternating ladder is not orderable for this spec.
    let orderable = orderable_upload(&project, &SUBSTRATE);

    Ok(LumpedDesigned {
        ladder,
        resonators,
        sweep,
        sweep_finite_q,
        q_unloaded,
        midband_il_db,
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
        orderable,
    })
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

// ===========================================================================
// Stepped-impedance low-pass adapter (App.2.2, ADR-0139)
// ===========================================================================

/// High-Z line impedance for the stepped-impedance dimensioner, ohms (the
/// standard high/low pair used in Pozar §8.6 Example 8.6). A practical FR-4
/// microstrip realizes roughly 20–120 Ω, so these are the impedances the
/// alternating series-inductor (high-Z) / shunt-capacitor (low-Z) sections use.
const STEPPED_Z_HIGH: f64 = 120.0;
/// Low-Z line impedance for the stepped-impedance dimensioner, ohms.
const STEPPED_Z_LOW: f64 = 20.0;
/// The low-pass response sweep spans `[0, SWEEP_SPAN_MULT·f_c]` so the −3 dB
/// cutoff and the stopband roll-off are both visible on the plot.
const LP_SWEEP_SPAN_MULT: f64 = 3.0;

/// One rendered transmission-line section of the stepped-impedance low-pass
/// filter (Synthesis / Layout tables), pre-extracted for display.
#[derive(Clone, Copy)]
pub struct SteppedSectionRow {
    /// Section index, 1-based (source → load).
    pub index: usize,
    /// `true` for a series-inductor **high-Z** line; `false` for a
    /// shunt-capacitor **low-Z** line.
    pub high_z: bool,
    /// Characteristic impedance of the section, ohms.
    pub z_ohm: f64,
    /// Electrical length `βl`, degrees.
    pub betal_deg: f64,
    /// Physical microstrip width, mm.
    pub width_mm: f64,
    /// Physical section length, mm.
    pub length_mm: f64,
}

/// Everything the two stepped-impedance low-pass stages render — all from the
/// live engine (the F1.2.3 [`dimension_stepped_impedance`] dimensioner + the
/// App.2.2 [`ideal_response_lowpass`] response). Mirrors [`LumpedDesigned`].
pub struct SteppedLowpassDesigned {
    /// The low-pass spec this was designed from (`Response::Lowpass`; `f0_hz`
    /// reused as the cutoff `f_c`).
    pub spec: FilterSpec,
    /// The filter order `N`.
    pub order: usize,
    /// Prototype element values `[g0, g1, …, gN, g_{N+1}]`.
    pub g_values: Vec<f64>,
    /// Display-ready stepped-section rows (source → load, low-Z first).
    pub sections: Vec<SteppedSectionRow>,
    /// The swept low-pass `|S21|`/`|S11|` response (reuses [`SweepPoint`]).
    pub sweep: Vec<SweepPoint>,
    /// Forbidden low-pass mask regions for the response plot.
    pub mask_bands: Vec<MaskBand>,
    /// `true` iff the swept response satisfies the low-pass spec mask.
    pub pass: bool,
    /// Worst-case in-band insertion-loss ripple observed, dB.
    pub worst_passband_ripple_db: f64,
    /// Worst-case (smallest) in-band return loss observed, dB.
    pub worst_return_loss_db: f64,
    /// Per stopband point: `(freq_hz, achieved_rejection_db, required_db, met)`.
    pub stopband: Vec<(f64, f64, f64, bool)>,
    /// The dimensioned stepped-impedance board, or `None` when the sections are
    /// not realizable on FR-4 (then [`dim_error`](Self::dim_error) is `Some`).
    pub layout: Option<Layout>,
    /// Board bounding box, mm (`(width, height)`); `(0, 0)` when not realizable.
    pub board_size_mm: (f64, f64),
    /// The dimensioning error string when the sections could not be realized on
    /// FR-4, else `None`. The synthesis / response / verdict stay real.
    pub dim_error: Option<String>,
}

impl SteppedLowpassDesigned {
    /// The high-Z line impedance the dimensioner used, ohms.
    pub fn z_high(&self) -> f64 {
        STEPPED_Z_HIGH
    }
    /// The low-Z line impedance the dimensioner used, ohms.
    pub fn z_low(&self) -> f64 {
        STEPPED_Z_LOW
    }
    /// The cutoff frequency `f_c`, Hz (the low-pass spec reuses `f0_hz`).
    pub fn cutoff_hz(&self) -> f64 {
        self.spec.f0_hz
    }
}

/// Run the full live stepped-impedance low-pass pipeline on the hard-coded demo
/// spec, mapped to a low-pass cutoff.
///
/// Convenience wrapper around [`design_stepped_from`] for the initial boot state
/// (the Spec stage edits a live spec that re-drives the pipeline). The band-pass
/// [`demo_spec`] is reused with its response switched to [`Response::Lowpass`]
/// (its `f0_hz` becomes the cutoff, its `fbw` is irrelevant to low-pass).
pub fn design_stepped() -> SteppedLowpassDesigned {
    design_stepped_from(stepped_demo_spec())
}

/// The demo spec adapted to a low-pass stepped-impedance design: the [`demo_spec`]
/// with [`Response::Lowpass`] and a stopband well above the cutoff (so the mask
/// has a meaningful rejection target for the low-pass roll-off).
fn stepped_demo_spec() -> FilterSpec {
    let mut spec = demo_spec();
    spec.response = Response::Lowpass;
    // f0 stays as the cutoff; place the stopband at ~2× cutoff with a target the
    // order-5 roll-off comfortably meets.
    spec.mask.stopband = vec![(2.0 * spec.f0_hz, 25.0)];
    spec
}

/// Run the full live stepped-impedance low-pass pipeline on an arbitrary
/// [`FilterSpec`] interpreted as a **low-pass** design (the `f0_hz` field is the
/// cutoff `f_c`): synthesize the prototype g-values → dimension the alternating
/// high-Z / low-Z microstrip sections ([`dimension_stepped_impedance`], Pozar
/// §8.6) → sweep the low-pass `|S21|` ([`ideal_response_lowpass`]) → grade
/// against a low-pass spec mask → assemble the in-line board layout. Every value
/// is real engine output; nothing is faked.
///
/// The synthesis / response / verdict stay real even when the geometry is not
/// realizable on FR-4 (an over-short / over-long section width-synthesis edge
/// case): the geometry-derived fields degrade to empty + a `dim_error`, mirroring
/// [`design_demo_from`].
pub fn design_stepped_from(spec: FilterSpec) -> SteppedLowpassDesigned {
    // Low-pass synthesis is the bare prototype g-values — there is no band-pass
    // coupling matrix or fractional bandwidth, so go straight to
    // `yee_synth::prototype` (the same g-values the dimensioner consumes) rather
    // than the band-pass `synthesize`. The order is the explicit spec order
    // (default 5 for the boot demo); a low-pass mask cannot be band-pass-mapped
    // to estimate an order, so an explicit order is required here.
    let order = spec.order.unwrap_or(5).max(1);
    let prototype = yee_synth::prototype(spec.approximation, order);
    let g_values = prototype.g.clone();
    let f_c = spec.f0_hz;
    let approx = spec.approximation;

    // ---- physical dimensioning (FR-4, stepped-impedance) ------------------
    let (sections, layout, board_size_mm, dim_error) = match (
        dimension_stepped_impedance(
            &prototype,
            f_c,
            spec.z0_ohm,
            STEPPED_Z_HIGH,
            STEPPED_Z_LOW,
            &SUBSTRATE,
        ),
        dimension_stepped_impedance_layout(
            &prototype,
            f_c,
            spec.z0_ohm,
            STEPPED_Z_HIGH,
            STEPPED_Z_LOW,
            &SUBSTRATE,
        ),
    ) {
        (Ok(dims), Ok(layout)) => {
            let rows: Vec<SteppedSectionRow> = dims
                .sections
                .iter()
                .enumerate()
                .map(|(i, s): (usize, &SteppedSection)| SteppedSectionRow {
                    index: i + 1,
                    high_z: s.high_z,
                    z_ohm: s.z_ohm,
                    betal_deg: s.electrical_length_rad.to_degrees(),
                    width_mm: s.width_m * 1e3,
                    length_mm: s.length_m * 1e3,
                })
                .collect();
            let bs = board_size_mm(&layout.bbox);
            (rows, Some(layout), bs, None)
        }
        (dims_res, layout_res) => {
            let msg = dims_res
                .err()
                .map(|e| e.to_string())
                .or_else(|| layout_res.err().map(|e| e.to_string()))
                .unwrap_or_else(|| "stepped-impedance geometry is not realizable on FR-4".into());
            (Vec::new(), None, (0.0, 0.0), Some(msg))
        }
    };

    // ---- swept low-pass response (Ω = f / f_c) ----------------------------
    let freqs = lowpass_sweep_freqs(f_c);
    let s21 = ideal_response_lowpass(approx, order, f_c, &freqs);
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

    let mask_bands = lowpass_mask_bands(&spec);

    // ---- low-pass spec-mask verdict ---------------------------------------
    // In-band is [0, f_c]: grade the worst ripple + the worst return loss over
    // that band, then check each stopband rejection point.
    let mut min_il = f64::INFINITY;
    let mut max_il = f64::NEG_INFINITY;
    let mut worst_rl = f64::INFINITY;
    let mut saw_passband = false;
    for s in sweep.iter().filter(|s| s.f_hz > 0.0 && s.f_hz <= f_c) {
        saw_passband = true;
        let il = -s.s21_db; // insertion loss (positive = loss)
        min_il = min_il.min(il);
        max_il = max_il.max(il);
        worst_rl = worst_rl.min(-s.s11_db);
    }
    let worst_passband_ripple_db = if saw_passband {
        (max_il - min_il).max(0.0)
    } else {
        0.0
    };
    let worst_return_loss_db = if worst_rl.is_finite() { worst_rl } else { 0.0 };

    let mut pass = saw_passband
        && worst_passband_ripple_db <= spec.mask.passband_ripple_db + 1e-9
        && worst_return_loss_db + 1e-9 >= spec.mask.return_loss_db;

    let stopband: Vec<(f64, f64, f64, bool)> = spec
        .mask
        .stopband
        .iter()
        .map(|&(f_s, required_db)| {
            let s21f = ideal_response_lowpass(approx, order, f_c, &[f_s]);
            let s21_mag = s21f[0].norm();
            let rejection_db = -20.0 * s21_mag.max(1e-12).log10();
            let met = rejection_db + 1e-9 >= required_db;
            if !met {
                pass = false;
            }
            (f_s, rejection_db, required_db, met)
        })
        .collect();

    SteppedLowpassDesigned {
        spec,
        order,
        g_values,
        sections,
        sweep,
        mask_bands,
        pass,
        worst_passband_ripple_db,
        worst_return_loss_db,
        stopband,
        layout,
        board_size_mm,
        dim_error,
    }
}

/// Linear low-pass sweep over `[0, LP_SWEEP_SPAN_MULT·f_c]` with `SWEEP_POINTS`
/// samples (the first point is a small epsilon above 0 so the dB-floor is well
/// defined and the cutoff/stopband roll-off are both on-screen).
fn lowpass_sweep_freqs(f_c: f64) -> Vec<f64> {
    let hi = (LP_SWEEP_SPAN_MULT * f_c).max(f_c * 1e-3);
    let lo = f_c * 1e-3;
    (0..SWEEP_POINTS)
        .map(|i| lo + (hi - lo) * (i as f64) / ((SWEEP_POINTS - 1) as f64))
        .collect()
}

/// Forbidden low-pass mask regions for the plot: a passband insertion-loss floor
/// at `−passband_ripple_db` over `[0, f_c]`, plus a rejection ceiling at
/// `−reject` over a ±2 % band at each stopband point above the cutoff.
fn lowpass_mask_bands(spec: &FilterSpec) -> Vec<MaskBand> {
    let f_c = spec.f0_hz;
    let mut bands = vec![MaskBand {
        f_lo_hz: f_c * 1e-3,
        f_hi_hz: f_c,
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
// TopBar view (App.2.3, ADR-0140)
// ===========================================================================

/// The human-readable approximation label for the TopBar summary, e.g.
/// `"Chebyshev 0.5 dB"` / `"Butterworth"`.
fn approximation_label(approx: Approximation) -> String {
    match approx {
        Approximation::Chebyshev { ripple_db } => format!("Chebyshev {ripple_db:.1} dB"),
        Approximation::Butterworth => "Butterworth".to_string(),
    }
}

/// The TopBar summary line + PASS/FAIL verdict for the **active flow**.
///
/// Pure and host-testable: it reads no signals, taking the three flow values by
/// reference, and is dispatched on the active [`Topology`] so the top bar
/// reflects the filter the user is actually designing (the §10 honest-UI
/// principle) rather than the always-band-pass `designed` signal. Returns the
/// summary string and the verdict; a `None` verdict means the active flow's
/// design is **not realizable** (e.g. an unrealizable lumped ladder), which the
/// `TopBar` component renders as a muted "geometry not realizable" chip.
///
/// Branches:
/// - [`Topology::EdgeCoupled`] / [`Topology::Hairpin`] → the band-pass summary
///   (`· {approx} · N={order} · {f0} GHz · {fbw}%`) + `Some(designed.report.pass)`.
/// - [`Topology::LumpedLc`] → the **same** band-pass summary (the lumped flow
///   shares the band-pass spec, read from `designed.spec`) + the lumped ladder
///   verdict `lumped.map(|l| l.verdict.pass)` (`None` → not realizable).
/// - [`Topology::SteppedImpedance`] → the **low-pass** summary
///   (`· {approx} · N={order} · cutoff {f_c} GHz`, no fractional bandwidth) +
///   `Some(stepped.pass)`.
pub fn topbar_view(
    topology: Topology,
    designed: &Designed,
    lumped: Option<&LumpedDesigned>,
    stepped: &SteppedLowpassDesigned,
) -> (String, Option<bool>) {
    match topology {
        Topology::EdgeCoupled | Topology::Hairpin | Topology::Combline | Topology::Interdigital => {
            let spec = &designed.spec;
            let summary = format!(
                "· {} · N={} · {:.2} GHz · {:.0}%",
                approximation_label(spec.approximation),
                designed.order(),
                spec.f0_hz / 1e9,
                spec.fbw * 100.0
            );
            (summary, Some(designed.report.pass))
        }
        Topology::LumpedLc => {
            // The lumped flow shares the band-pass spec; the summary mirrors the
            // band-pass branch, but the verdict is the lumped ladder's own (and
            // `None` when the ladder is not realizable).
            let spec = &designed.spec;
            let summary = format!(
                "· {} · N={} · {:.2} GHz · {:.0}%",
                approximation_label(spec.approximation),
                designed.order(),
                spec.f0_hz / 1e9,
                spec.fbw * 100.0
            );
            (summary, lumped.map(|l| l.verdict.pass))
        }
        Topology::SteppedImpedance => {
            let summary = format!(
                "· {} · N={} · cutoff {:.2} GHz",
                approximation_label(stepped.spec.approximation),
                stepped.order,
                stepped.cutoff_hz() / 1e9
            );
            (summary, Some(stepped.pass))
        }
    }
}

// ===========================================================================
// Verify view (App.2.4, ADR-0141)
// ===========================================================================

/// What level of verification the [`VerifyView`] metrics represent.
///
/// The studio grades at the **circuit / synthesis** level; full-wave EM
/// verification of the physical board is a separate native step (the deferred
/// ADR-0133 research frontier), not run in the browser. This enum names which
/// of the two circuit-level checks produced the metrics so the Verify stage can
/// state it honestly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerifyLevel {
    /// The **realized** LC ladder graded vs the mask (the lumped flow — a genuine
    /// circuit-level check on the synthesized `ladder_s21` response).
    RealizedLadder,
    /// The **synthesized ideal** / coupled-resonator response graded vs the mask
    /// (the distributed band-pass + stepped low-pass flows). The physical board's
    /// full-wave EM response is a native step, not computed here.
    SynthesizedIdeal,
}

/// The active flow's real verification metrics for the Verify stage.
///
/// Every field is a value the engine already computes for the active flow's
/// graded response (`MaskReport` / `MaskVerdict` / the stepped low-pass fields);
/// nothing here is fabricated. The Verify stage renders these directly.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VerifyView {
    /// Which circuit-level check produced these metrics (realized ladder vs
    /// synthesized ideal response).
    pub level: VerifyLevel,
    /// Overall verdict: `Some(true)` = the active flow meets its mask,
    /// `Some(false)` = it fails, `None` = the active flow's design is **not
    /// realizable** (e.g. an unrealizable lumped ladder), so no verdict exists.
    pub pass: Option<bool>,
    /// Worst-case in-band passband insertion-loss ripple observed, dB.
    pub worst_passband_ripple_db: f64,
    /// Worst-case (smallest) in-band return loss observed, dB.
    pub worst_return_loss_db: f64,
    /// Worst-case (smallest) stopband rejection across the mask's stopband
    /// points, dB. `None` when the mask has no stopband points (nothing to
    /// report — rendered as "—", not a fabricated number).
    pub worst_stopband_rej_db: Option<f64>,
}

/// The verification metrics + level for the **active flow**.
///
/// Pure and host-testable: it reads no signals, taking the three flow values by
/// reference, and dispatches on the active [`Topology`] so the Verify stage
/// surfaces the filter the user is actually designing (the §10 honest-UI
/// principle). Every metric is pulled from the flow's already-computed graded
/// response — no fabricated placeholders.
///
/// Branches:
/// - [`Topology::LumpedLc`] → [`VerifyLevel::RealizedLadder`] from
///   `lumped.verdict` (a genuine realized-ladder check). When `lumped` is `None`
///   (the ladder is not realizable) the metrics are `0.0` / `None` rejection and
///   `pass` is `None`. The `verdict.worst_stopband_rej_db` is `+∞` when the mask
///   has no stopband points; that maps to `None`.
/// - [`Topology::SteppedImpedance`] → [`VerifyLevel::SynthesizedIdeal`] from the
///   stepped low-pass fields; stopband rejection is the minimum achieved over
///   `stepped.stopband` (`None` when empty).
/// - [`Topology::EdgeCoupled`] / [`Topology::Hairpin`] →
///   [`VerifyLevel::SynthesizedIdeal`] from `designed.report`; stopband rejection
///   is the minimum achieved over `report.stopband` (`None` when empty).
pub fn verify_view(
    topology: Topology,
    designed: &Designed,
    lumped: Option<&LumpedDesigned>,
    stepped: &SteppedLowpassDesigned,
) -> VerifyView {
    /// Minimum achieved rejection over a `(freq, achieved, required, met)`
    /// stopband table (`None` when the mask has no stopband points).
    fn min_achieved(stopband: &[(f64, f64, f64, bool)]) -> Option<f64> {
        stopband
            .iter()
            .map(|&(_, achieved, _, _)| achieved)
            .fold(None, |acc, a| Some(acc.map_or(a, |m: f64| m.min(a))))
    }

    match topology {
        Topology::LumpedLc => match lumped {
            Some(l) => VerifyView {
                level: VerifyLevel::RealizedLadder,
                pass: Some(l.verdict.pass),
                worst_passband_ripple_db: l.verdict.worst_passband_ripple_db,
                worst_return_loss_db: l.verdict.worst_return_loss_db,
                // `worst_stopband_rej_db` is `+∞` when the mask has no stopband
                // points — report that as "no stopband point" (`None`), not ∞.
                worst_stopband_rej_db: l
                    .verdict
                    .worst_stopband_rej_db
                    .is_finite()
                    .then_some(l.verdict.worst_stopband_rej_db),
            },
            None => VerifyView {
                level: VerifyLevel::RealizedLadder,
                pass: None,
                worst_passband_ripple_db: 0.0,
                worst_return_loss_db: 0.0,
                worst_stopband_rej_db: None,
            },
        },
        Topology::SteppedImpedance => VerifyView {
            level: VerifyLevel::SynthesizedIdeal,
            pass: Some(stepped.pass),
            worst_passband_ripple_db: stepped.worst_passband_ripple_db,
            worst_return_loss_db: stepped.worst_return_loss_db,
            worst_stopband_rej_db: min_achieved(&stepped.stopband),
        },
        Topology::EdgeCoupled | Topology::Hairpin | Topology::Combline | Topology::Interdigital => {
            VerifyView {
                level: VerifyLevel::SynthesizedIdeal,
                pass: Some(designed.report.pass),
                worst_passband_ripple_db: designed.report.worst_passband_ripple_db,
                worst_return_loss_db: designed.report.worst_return_loss_db,
                worst_stopband_rej_db: min_achieved(&designed.report.stopband),
            }
        }
    }
}

// ===========================================================================
// Compare-techniques view (App.2.5, ADR-0142)
// ===========================================================================

/// One comparable row in the Technique-stage compare view: a single realization
/// technique synthesized for the current spec, with its board size, verdict, and
/// key graded metrics pulled straight from that technique's design.
///
/// Every numeric field is **real engine output** (the same graded structs the
/// App.2.4 [`verify_view`] reads); nothing is fabricated. When a design fails to
/// dimension (an unrealizable lumped ladder, or a distributed layout that cannot
/// be realized on FR-4), [`realizable`](Self::realizable) is `false`, `pass` is
/// `None`, and the metric fields are zeroed (`worst_stopband_rej_db` is `None`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TechniqueComparison {
    /// The realization technique this row was synthesized for.
    pub technique: RealizationTechnique,
    /// `false` when the design failed to dimension (no realizable board /
    /// ladder). The synthesized response metrics may still be real, but the
    /// physical realization did not produce a board.
    pub realizable: bool,
    /// Board bounding-box width, mm (`0.0` when not realizable).
    pub board_w_mm: f64,
    /// Board bounding-box height, mm (`0.0` when not realizable).
    pub board_h_mm: f64,
    /// Overall spec-mask verdict for this technique's design: `Some(true)` =
    /// meets the mask, `Some(false)` = fails, `None` = not realizable.
    pub pass: Option<bool>,
    /// The filter order `N`.
    pub order: usize,
    /// Worst-case in-band passband insertion-loss ripple observed, dB.
    pub worst_passband_ripple_db: f64,
    /// Worst-case (smallest) in-band return loss observed, dB.
    pub worst_return_loss_db: f64,
    /// Worst-case (smallest) stopband rejection across the mask's stopband
    /// points, dB. `None` when the mask has no stopband points (rendered as
    /// "—", never a fabricated number).
    pub worst_stopband_rej_db: Option<f64>,
}

/// Synthesize **every live technique that realizes `spec`'s response class** and
/// collect a comparable metric row for each — board size, PASS/FAIL, order, and
/// the worst passband ripple / return loss / stopband rejection — pulled directly
/// from each design's existing graded struct (the same fields [`verify_view`]
/// reads; real engine output).
///
/// Pure (reads no signals): it builds the relevant designs from `spec` and
/// returns one [`TechniqueComparison`] per technique. Building the (light,
/// synchronous) designs per call is intentional — the studio re-derives live.
///
/// Keyed on [`FilterSpec::response`]:
/// - [`Response::Bandpass`] / [`Response::Bandstop`] → four rows:
///   [`RealizationTechnique::EdgeCoupled`], [`RealizationTechnique::Hairpin`],
///   and [`RealizationTechnique::Combline`] (via [`design_demo_from`]; metrics
///   from [`Designed::report`], `realizable = layout.is_some()`), and
///   [`RealizationTechnique::LumpedLc`] (via [`design_lumped_from`]; metrics from
///   [`LumpedDesigned::verdict`] on `Ok`, `realizable = false` with zeroed
///   metrics on `Err`).
/// - [`Response::Lowpass`] → one row,
///   [`RealizationTechnique::SteppedImpedance`] (via [`design_stepped_from`];
///   metrics from the stepped low-pass fields).
/// - [`Response::Highpass`] → `[]` (no live technique realizes this yet).
pub fn compare_techniques(spec: &FilterSpec) -> Vec<TechniqueComparison> {
    /// Minimum achieved rejection over a `(freq, achieved, required, met)`
    /// stopband table (`None` when the mask has no stopband points).
    fn min_achieved(stopband: &[(f64, f64, f64, bool)]) -> Option<f64> {
        stopband
            .iter()
            .map(|&(_, achieved, _, _)| achieved)
            .fold(None, |acc, a| Some(acc.map_or(a, |m: f64| m.min(a))))
    }

    /// A `TechniqueComparison` row from a distributed [`Designed`] (edge-coupled
    /// or hairpin): metrics from the synthesized-response mask report, the board
    /// from the dimensioned layout (`realizable = layout.is_some()`).
    fn from_distributed(technique: RealizationTechnique, d: &Designed) -> TechniqueComparison {
        // When the coupling cannot be dimensioned on FR-4 the row degrades per
        // the struct contract: realizable=false, pass=None, zeroed board +
        // metrics (the order stays — a synthesis property, not geometry).
        // Mirrors the lumped `Err` arm.
        if d.layout.is_none() {
            return TechniqueComparison {
                technique,
                realizable: false,
                board_w_mm: 0.0,
                board_h_mm: 0.0,
                pass: None,
                order: d.order(),
                worst_passband_ripple_db: 0.0,
                worst_return_loss_db: 0.0,
                worst_stopband_rej_db: None,
            };
        }
        let (bw, bh) = d.board_size_mm;
        TechniqueComparison {
            technique,
            realizable: true,
            board_w_mm: bw,
            board_h_mm: bh,
            pass: Some(d.report.pass),
            order: d.order(),
            worst_passband_ripple_db: d.report.worst_passband_ripple_db,
            worst_return_loss_db: d.report.worst_return_loss_db,
            worst_stopband_rej_db: min_achieved(&d.report.stopband),
        }
    }

    match spec.response {
        Response::Bandpass | Response::Bandstop => {
            let edge = design_demo_from(spec.clone(), Topology::EdgeCoupled);
            let hairpin = design_demo_from(spec.clone(), Topology::Hairpin);
            let combline = design_demo_from(spec.clone(), Topology::Combline);
            let interdigital = design_demo_from(spec.clone(), Topology::Interdigital);
            // The lumped flow has its own fallible ladder synthesis; an
            // unrealizable ladder degrades to a not-realizable row (zeroed
            // metrics, no verdict) rather than panicking.
            let lumped_row = match design_lumped_from(spec.clone(), DEFAULT_Q_UNLOADED) {
                Ok(l) => {
                    let (bw, bh) = l.board_size_mm;
                    TechniqueComparison {
                        technique: RealizationTechnique::LumpedLc,
                        realizable: true,
                        board_w_mm: bw,
                        board_h_mm: bh,
                        pass: Some(l.verdict.pass),
                        order: l.order(),
                        worst_passband_ripple_db: l.verdict.worst_passband_ripple_db,
                        worst_return_loss_db: l.verdict.worst_return_loss_db,
                        // `worst_stopband_rej_db` is `+∞` when the mask has no
                        // stopband points — map that to `None`, never ∞.
                        worst_stopband_rej_db: l
                            .verdict
                            .worst_stopband_rej_db
                            .is_finite()
                            .then_some(l.verdict.worst_stopband_rej_db),
                    }
                }
                Err(_) => TechniqueComparison {
                    technique: RealizationTechnique::LumpedLc,
                    realizable: false,
                    board_w_mm: 0.0,
                    board_h_mm: 0.0,
                    pass: None,
                    order: 0,
                    worst_passband_ripple_db: 0.0,
                    worst_return_loss_db: 0.0,
                    worst_stopband_rej_db: None,
                },
            };
            vec![
                from_distributed(RealizationTechnique::EdgeCoupled, &edge),
                from_distributed(RealizationTechnique::Hairpin, &hairpin),
                from_distributed(RealizationTechnique::Combline, &combline),
                from_distributed(RealizationTechnique::Interdigital, &interdigital),
                lumped_row,
            ]
        }
        Response::Lowpass => {
            let stepped = design_stepped_from(spec.clone());
            let (bw, bh) = stepped.board_size_mm;
            vec![TechniqueComparison {
                technique: RealizationTechnique::SteppedImpedance,
                realizable: stepped.layout.is_some(),
                board_w_mm: bw,
                board_h_mm: bh,
                pass: Some(stepped.pass),
                order: stepped.order,
                worst_passband_ripple_db: stepped.worst_passband_ripple_db,
                worst_return_loss_db: stepped.worst_return_loss_db,
                worst_stopband_rej_db: min_achieved(&stepped.stopband),
            }]
        }
        // No live technique realizes a high-pass response yet.
        Response::Highpass => vec![],
    }
}

// ===========================================================================
// Response-overlay curves (App.2.6, ADR-0143)
// ===========================================================================

/// One labelled response curve for the Compare-panel overlay.
///
/// Each curve is a real swept response on the shared `sweep_freqs` grid (the
/// same [`SweepPoint`]s the per-flow Synthesis stages plot); nothing is faked.
/// When a realization fails to synthesize (an unrealizable lumped ladder),
/// [`realizable`](Self::realizable) is `false` and [`sweep`](Self::sweep) is
/// empty — the overlay draws no polyline for it but still names it in the
/// legend.
#[derive(Clone)]
pub struct OverlayCurve {
    /// Honest curve label, e.g.
    /// `"Coupled-resonator (edge-coupled / hairpin / combline / interdigital) — ideal"`.
    pub label: String,
    /// The swept response (empty when not realizable).
    pub sweep: Vec<SweepPoint>,
    /// `false` when the realization could not be synthesized (empty sweep).
    pub realizable: bool,
}

/// The **distinct** swept `|S21|` responses to overlay for `spec`, on the shared
/// `sweep_freqs` grid — labelled truthfully (the §Honesty-constraint).
///
/// Pure (reads no signals): it builds the relevant designs from `spec` and
/// returns one [`OverlayCurve`] per genuinely-distinct response. Edge-coupled
/// and hairpin share the **same** coupled-resonator synthesis (identical
/// coupling matrix → identical ideal `|S21|`); they differ only physically
/// (board layout/size — already in the compare table), so they are a **single**
/// shared ideal curve, never two relabelled copies.
///
/// Keyed on [`FilterSpec::response`]:
/// - [`Response::Bandpass`] / [`Response::Bandstop`] → two curves: the
///   coupled-resonator ideal (the [`Designed::sweep`] from
///   [`design_demo_from`], labelled as edge-coupled / hairpin / combline — they
///   share the same coupling matrix → identical ideal `|S21|`; always
///   realizable — it is the synthesized response) and the lumped realized
///   ladder (the [`LumpedDesigned::sweep`] = `ladder_s21` from
///   [`design_lumped_from`]; on `Err` an empty/not-realizable curve so the
///   legend stays honest).
/// - [`Response::Lowpass`] → one curve: the stepped-impedance ideal
///   (the [`SteppedLowpassDesigned::sweep`] from [`design_stepped_from`]).
/// - [`Response::Highpass`] → `[]` (no live technique realizes this yet).
pub fn overlay_curves(spec: &FilterSpec) -> Vec<OverlayCurve> {
    match spec.response {
        Response::Bandpass | Response::Bandstop => {
            // The coupled-resonator ideal is topology-independent (edge-coupled,
            // hairpin, combline, and interdigital synthesize the *same* response),
            // so a single curve labelled as all four — never relabelled copies.
            let coupled = design_demo_from(spec.clone(), Topology::EdgeCoupled);
            let coupled_curve = OverlayCurve {
                label:
                    "Coupled-resonator (edge-coupled / hairpin / combline / interdigital) — ideal"
                        .to_string(),
                sweep: coupled.sweep,
                realizable: true,
            };
            // The lumped realized ladder is a genuinely distinct response
            // (`ladder_s21`); an unrealizable ladder degrades to an empty,
            // not-realizable curve (no polyline, but still named in the legend).
            let lumped_curve = match design_lumped_from(spec.clone(), DEFAULT_Q_UNLOADED) {
                Ok(l) => OverlayCurve {
                    label: "Lumped LC — realized ladder".to_string(),
                    sweep: l.sweep,
                    realizable: true,
                },
                Err(_) => OverlayCurve {
                    label: "Lumped LC — realized ladder".to_string(),
                    sweep: Vec::new(),
                    realizable: false,
                },
            };
            vec![coupled_curve, lumped_curve]
        }
        Response::Lowpass => {
            let stepped = design_stepped_from(spec.clone());
            vec![OverlayCurve {
                label: "Stepped-impedance — ideal".to_string(),
                sweep: stepped.sweep,
                realizable: true,
            }]
        }
        // No live technique realizes a high-pass response yet.
        Response::Highpass => vec![],
    }
}

/// The forbidden spec-mask regions for the Compare-panel overlay, matching the
/// spec's response class so the shaded mask lines up with the overlaid
/// responses: the band-pass [`mask_bands`] for band-pass/band-stop, the low-pass
/// [`lowpass_mask_bands`] for low-pass, and none for high-pass (which has no
/// live overlay). Mirrors how each per-flow Synthesis stage chooses its bands.
pub fn overlay_mask_bands(spec: &FilterSpec) -> Vec<MaskBand> {
    match spec.response {
        Response::Bandpass | Response::Bandstop => mask_bands(spec),
        Response::Lowpass => lowpass_mask_bands(spec),
        Response::Highpass => Vec::new(),
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
        // Board has positive extent and a realizable layout.
        assert!(d.board_size_mm.0 > 0.0 && d.board_size_mm.1 > 0.0);
        assert!(d.layout.is_some(), "demo spec dimensions onto FR-4");
        assert!(d.dim_error.is_none(), "no dimensioning error for the demo");
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

    #[test]
    fn finite_q_midband_il_matches_cohn() {
        // ADR-0163, mirroring `lumped-q-001` at the studio-engine layer: the
        // engine's realized finite-Q MIDBAND INSERTION LOSS (from
        // `ladder_s21_lossy` via `design_lumped_from`) must match Cohn's narrow-
        // band dissipation formula `IL ≈ 4.343·Σg/(Q_u·FBW)` (Cohn 1959; Hong &
        // Lancaster §3.2), and the IDEAL (lossless) midband IL must stay ≈ 0 dB.
        //
        // Standard 3-pole 0.5 dB Chebyshev band-pass, FBW = 10 %, Q_u = 100.
        let q_u = DEFAULT_Q_UNLOADED;
        let fbw = 0.10;
        let approx = Approximation::Chebyshev { ripple_db: 0.5 };
        let mut spec = demo_spec();
        spec.approximation = approx;
        spec.order = Some(3);
        spec.fbw = fbw;

        let l = design_lumped_from(spec, q_u).expect("3-pole 0.5 dB Cheb BPF is realizable");
        assert_eq!(l.order(), 3, "3-pole ladder");
        assert_eq!(l.q_unloaded, q_u, "engine carries the requested Q_u");
        // The finite-Q sweep shares the ideal sweep's grid.
        assert_eq!(l.sweep_finite_q.len(), l.sweep.len());
        assert_eq!(l.sweep_finite_q.len(), SWEEP_POINTS);

        // Cohn reference: Σg over the reactive elements g[1..=N], computed
        // INDEPENDENTLY from the prototype (non-circular — not from the sweep).
        let proto = yee_synth::prototype(approx, 3);
        let sum_g: f64 = proto.g[1..=3].iter().sum();
        let il_cohn = 4.343 * sum_g / (q_u * fbw);

        // The engine's realized midband IL (from the lossy ABCD response).
        let il_q = l.midband_il_db;

        // The ideal lossless midband IL: the |S21| dB at the sweep point closest
        // to f0 (a passband Cheb minimum is ≈ 0 dB; allow the ripple + grid slack).
        let il_ideal = -l
            .sweep
            .iter()
            .min_by(|a, b| {
                (a.f_hz - 2.0e9)
                    .abs()
                    .partial_cmp(&(b.f_hz - 2.0e9).abs())
                    .unwrap()
            })
            .map(|p| p.s21_db)
            .unwrap();

        println!(
            "finite-Q engine IL: IL_q = {il_q:.4} dB, IL_cohn = {il_cohn:.4} dB, \
             IL_ideal = {il_ideal:.4} dB (Σg = {sum_g:.4}, Q_u = {q_u}, FBW = {fbw})"
        );

        let rel_err = (il_q - il_cohn).abs() / il_cohn;
        assert!(
            rel_err <= 0.15,
            "finite-Q midband IL {il_q:.4} dB vs Cohn {il_cohn:.4} dB \
             (rel err {:.1} % > 15 %)",
            rel_err * 100.0
        );
        // The lossless midband IL is essentially zero (passband floor).
        assert!(
            il_ideal.abs() <= 0.2,
            "ideal midband IL {il_ideal:.4} dB should be ≈ 0 dB"
        );
    }

    /// Parse a rendered Touchstone string back into its `(option_line, rows)`
    /// where each row is `(f_hz, S11, S21)` decoded from the on-disk RI columns.
    /// Splits the lines by hand (no fs), so a successful decode reproducing the
    /// input is a NON-circular round-trip check. Assumes a 2-port RI Hz file
    /// (the only thing [`s2p_string`] emits): the on-disk row layout for n=2 is
    /// `f  S11.re S11.im  S21.re S21.im  S12.re S12.im  S22.re S22.im`
    /// (the Touchstone v1 swapped-off-diagonal order — S21 before S12).
    fn parse_s2p_ri(s: &str) -> (String, Vec<(f64, Complex64, Complex64)>) {
        let mut option_line = String::new();
        let mut rows = Vec::new();
        for line in s.lines() {
            let t = line.trim();
            if t.is_empty() || t.starts_with('!') {
                continue;
            }
            if let Some(opt) = t.strip_prefix('#') {
                option_line = opt.trim().to_string();
                continue;
            }
            let nums: Vec<f64> = t
                .split_whitespace()
                .map(|tok| tok.parse::<f64>().expect("data token is a float"))
                .collect();
            assert_eq!(
                nums.len(),
                9,
                "2-port RI data row has 1 freq + 2·n² = 9 floats, got {}: {t}",
                nums.len()
            );
            let f = nums[0];
            let s11 = Complex64::new(nums[1], nums[2]);
            let s21 = Complex64::new(nums[3], nums[4]); // on-disk: S21 is slot 2
            rows.push((f, s11, s21));
        }
        (option_line, rows)
    }

    #[test]
    fn s2p_export_roundtrips_lumped_and_renders_distributed() {
        // ADR-0171 (T8), the studio-engine `.s2p` gate. Build the lumped
        // finite-Q `(S11, S21)` sweep for the demo fixture, render it via the
        // real `s2p_string` (→ `yee_io::touchstone::to_string`), and assert the
        // rendered TEXT reproduces the input S-parameters — a non-circular
        // round-trip mirroring the CLI's `.s2p` contract at the studio layer.
        let d = design_lumped();
        let z0 = d.ladder.z0_ohm;
        let (freqs, s_params) = lumped_s2p_sweep(&d.ladder, d.q_unloaded);
        assert_eq!(freqs.len(), SWEEP_POINTS);
        assert_eq!(freqs.len(), s_params.len());

        let text = lumped_s2p(&d).expect("lumped .s2p renders for the demo fixture");

        // Diagnostic sample: the comment, option line, and first data row.
        for line in text.lines().take(3) {
            println!("s2p | {line}");
        }

        // (a) The option line is the expected Touchstone 2-port RI / Hz header.
        let expected_opt = format!("Hz S RI R {z0}");
        let (opt, rows) = parse_s2p_ri(&text);
        assert_eq!(
            opt, expected_opt,
            "option line should be `# Hz S RI R {z0}` (got `# {opt}`)"
        );
        // The raw header substring is present too (so a downstream tool reading
        // the literal line, not our re-tokenization, also sees it).
        assert!(
            text.contains(&format!("# {expected_opt}\n")),
            "rendered text carries the literal `# {expected_opt}` option line"
        );

        // (b) One data row per swept frequency.
        assert_eq!(
            rows.len(),
            freqs.len(),
            "data-row count equals the frequency count"
        );

        // (c) Re-parsing the RI columns of every row recovers the input sweep:
        // the frequency (Hz) and BOTH S11 and S21 within a tight float tol. This
        // is the non-circular core — the rendered string reproduces the inputs.
        for (i, ((f_in, (s11_in, s21_in)), (f_out, s11_out, s21_out))) in freqs
            .iter()
            .zip(s_params.iter())
            .zip(rows.iter())
            .enumerate()
        {
            assert!(
                (f_in - f_out).abs() <= 1.0,
                "row {i}: freq {f_out} Hz != input {f_in} Hz"
            );
            assert!(
                (s11_in - s11_out).norm() <= 1e-9,
                "row {i}: S11 {s11_out} != input {s11_in}"
            );
            assert!(
                (s21_in - s21_out).norm() <= 1e-9,
                "row {i}: S21 {s21_out} != input {s21_in}"
            );
        }

        // The finite-Q lumped pair is the TRUE lossy 2-port → strictly passive
        // (|S11|² + |S21|² < 1) at every sample, so the file re-parses (and is a
        // physically-meaningful absorptive network, not a lossless placeholder).
        for (s11, s21) in &s_params {
            assert!(
                s11.norm_sqr() + s21.norm_sqr() <= 1.0 + 1e-9,
                "finite-Q sample is passive"
            );
        }
        // And the displayed finite-Q |S21| curve uses the SAME model → the .s2p
        // matches the plot: |S21| dB from the rendered pair equals the carried
        // sweep value at band-centre.
        let idx_f0 = freqs
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                (**a - d.ladder.f0_hz)
                    .abs()
                    .partial_cmp(&(**b - d.ladder.f0_hz).abs())
                    .unwrap()
            })
            .map(|(i, _)| i)
            .unwrap();
        let s21_db_rendered = 20.0 * rows[idx_f0].2.norm().max(1e-12).log10();
        assert!(
            (s21_db_rendered - d.sweep_finite_q[idx_f0].s21_db).abs() <= 1e-9,
            ".s2p S21 matches the displayed finite-Q sweep at f0"
        );

        // The DISTRIBUTED / ideal matrix also renders without error (it uses the
        // lossless quadrature S11, so it is passive: |S11|² + |S21|² = 1).
        let dd = design_demo();
        let dist_text =
            distributed_s2p(&dd).expect("distributed .s2p renders for the demo fixture");
        let (dist_opt, dist_rows) = parse_s2p_ri(&dist_text);
        assert_eq!(dist_opt, format!("Hz S RI R {}", dd.spec.z0_ohm));
        assert_eq!(dist_rows.len(), SWEEP_POINTS);
        // Quadrature pair is exactly passive at every sample.
        for (_f, s11, s21) in &dist_rows {
            assert!(
                (s11.norm_sqr() + s21.norm_sqr() - 1.0).abs() <= 1e-6,
                "ideal quadrature sample is unitary-passive"
            );
        }
    }

    /// Flatten a layout's copper geometry into a comparable signature: the
    /// board bounding box plus every trace vertex. Two distributed realizations
    /// of the same coupling network produce *different* signatures iff the
    /// boards differ physically.
    fn layout_signature(layout: &Layout) -> (Vec<f64>, Vec<(f64, f64)>) {
        let bbox = vec![
            layout.bbox.min.x,
            layout.bbox.min.y,
            layout.bbox.max.x,
            layout.bbox.max.y,
        ];
        let mut verts: Vec<(f64, f64)> = layout
            .traces
            .iter()
            .flat_map(|t| t.verts.iter().map(|p| (p.x, p.y)))
            .collect();
        verts.sort_by(|a, b| a.partial_cmp(b).unwrap());
        (bbox, verts)
    }

    #[test]
    fn hairpin_card_routes_to_real_dimensioner() {
        // The Hairpin gallery card must drive the REAL `dimension_hairpin`
        // engine, not a stub or an edge-coupled clone. The demo spec must
        // dimension as a hairpin (Some(layout)), and that hairpin layout must
        // DIFFER from the edge-coupled layout for the same spec — otherwise the
        // card would be silently routing to the edge-coupled stand-in.
        let spec = demo_spec();
        let hairpin = design_demo_from(spec.clone(), Topology::Hairpin);
        let edge = design_demo_from(spec, Topology::EdgeCoupled);

        // Hairpin dimensions onto FR-4.
        assert!(
            hairpin.layout.is_some(),
            "demo spec dimensions as a hairpin board"
        );
        assert!(hairpin.dim_error.is_none(), "no hairpin dimensioning error");
        assert_eq!(hairpin.topology, Topology::Hairpin);
        assert_eq!(hairpin.resonators.len(), 5, "order-5 hairpin → 5 arms");
        assert!(hairpin.board_size_mm.0 > 0.0 && hairpin.board_size_mm.1 > 0.0);

        // The shared coupled-resonator synthesis is topology-INDEPENDENT: same
        // coupling matrix, same g-values, same swept response, same verdict.
        assert_eq!(
            hairpin.coupling.m, edge.coupling.m,
            "shared coupling matrix"
        );
        assert_eq!(hairpin.report.pass, edge.report.pass, "shared verdict");

        // ...but the realized boards DIFFER (the card routes to the real
        // hairpin dimensioner). A folded λ/4-arm hairpin is geometrically
        // distinct from a straight λ/2 edge-coupled board.
        let (h_layout, e_layout) = (
            hairpin.layout.as_ref().unwrap(),
            edge.layout.as_ref().unwrap(),
        );
        assert_ne!(
            layout_signature(h_layout),
            layout_signature(e_layout),
            "hairpin layout must differ from the edge-coupled layout (real, not a clone)"
        );

        // The hairpin arm length is the U-folded ≈λ_g/4, roughly half the
        // edge-coupled straight λ_g/2 resonator length — a concrete witness
        // that the hairpin dimensioner (not the edge-coupled one) produced it.
        let h_len = hairpin.resonators[0].length_mm;
        let e_len = edge.resonators[0].length_mm;
        assert!(
            h_len < e_len * 0.7,
            "hairpin arm length ({h_len:.2} mm) ≈ λ/4 is well under the edge-coupled λ/2 ({e_len:.2} mm)"
        );
    }

    #[test]
    fn combline_card_routes_to_real_dimensioner() {
        // The Combline gallery card (App.2.7, ADR-0146) must drive the REAL
        // `dimension_combline` / `dimension_combline_layout` engine — not a stub,
        // not an edge-coupled clone, not a hairpin clone. The demo spec must
        // dimension as a combline board (Some(layout)) whose layout DIFFERS from
        // BOTH the edge-coupled and hairpin boards for the same spec, while the
        // topology-independent synthesis (coupling matrix + verdict) is SHARED,
        // and the combline-distinct loading cap C_L is real (> 0, finite).
        let spec = demo_spec();
        let combline = design_demo_from(spec.clone(), Topology::Combline);
        let edge = design_demo_from(spec.clone(), Topology::EdgeCoupled);
        let hairpin = design_demo_from(spec, Topology::Hairpin);

        // Combline dimensions onto FR-4.
        assert!(
            combline.layout.is_some(),
            "demo spec dimensions as a combline board"
        );
        assert!(
            combline.dim_error.is_none(),
            "no combline dimensioning error"
        );
        assert_eq!(combline.topology, Topology::Combline);
        assert_eq!(
            combline.resonators.len(),
            5,
            "order-5 combline → 5 resonators"
        );
        assert!(combline.board_size_mm.0 > 0.0 && combline.board_size_mm.1 > 0.0);

        // The shared coupled-resonator synthesis is topology-INDEPENDENT: same
        // coupling matrix and same verdict as the edge-coupled realization (the
        // combline card does NOT change the synthesis — only the realization).
        assert_eq!(
            combline.coupling.m, edge.coupling.m,
            "combline shares the edge-coupled coupling matrix"
        );
        assert_eq!(
            combline.report.pass, edge.report.pass,
            "combline shares the edge-coupled verdict"
        );

        // ...but the realized board DIFFERS from BOTH edge-coupled AND hairpin
        // (the card routes to the real combline dimensioner — short θ0 = λg/8
        // lines, not the straight λ/2 edge-coupled lines nor the U-folded λ/4
        // hairpin arms). A stub / clone of either would make these EQUAL.
        let cb_layout = combline.layout.as_ref().unwrap();
        let e_layout = edge.layout.as_ref().unwrap();
        let h_layout = hairpin.layout.as_ref().unwrap();
        assert_ne!(
            layout_signature(cb_layout),
            layout_signature(e_layout),
            "combline layout must differ from the edge-coupled layout (real, not a clone)"
        );
        assert_ne!(
            layout_signature(cb_layout),
            layout_signature(h_layout),
            "combline layout must differ from the hairpin layout (real, not a clone)"
        );

        // The combline-distinct loading cap C_L is surfaced and physical: a single
        // value (uniform θ0 / Z0 → the same cap per resonator), strictly positive
        // and finite. A faked / missing cap fails this.
        let c_l = combline
            .combline_loading_cap_f
            .expect("combline surfaces a loading cap C_L");
        assert!(
            c_l > 0.0 && c_l.is_finite(),
            "combline loading cap C_L = {c_l} F must be > 0 and finite"
        );
        // The edge-coupled / hairpin realizations have NO loading cap.
        assert!(
            edge.combline_loading_cap_f.is_none(),
            "edge-coupled has no loading cap"
        );
        assert!(
            hairpin.combline_loading_cap_f.is_none(),
            "hairpin has no loading cap"
        );

        // The combline θ0 = π/4 resonator is a short (≈ λg/8) line — well under
        // the straight λ/2 edge-coupled resonator — a concrete witness that the
        // combline dimensioner (not the edge-coupled one) produced the geometry.
        let cb_len = combline.resonators[0].length_mm;
        let e_len = edge.resonators[0].length_mm;
        assert!(
            cb_len < e_len * 0.5,
            "combline θ0=λg/8 resonator ({cb_len:.2} mm) is well under the edge-coupled λ/2 ({e_len:.2} mm)"
        );
    }

    #[test]
    fn interdigital_card_routes_to_real_dimensioner() {
        // The Interdigital gallery card (App.2.8, ADR-0150) must drive the REAL
        // `dimension_interdigital` / `dimension_interdigital_layout` engine — not a
        // stub, not an edge-coupled / hairpin clone, and NOT a combline-with-cap.
        // The demo spec must dimension as an interdigital board (Some(layout))
        // whose layout DIFFERS from the edge-coupled, hairpin, AND combline boards
        // for the same spec, while the topology-independent synthesis (coupling
        // matrix + verdict) is SHARED, and — the interdigital-distinct point —
        // there is NO loading cap (`combline_loading_cap_f` is None) while the
        // λg/4 resonator length is surfaced (> 0).
        let spec = demo_spec();
        let interdigital = design_demo_from(spec.clone(), Topology::Interdigital);
        let edge = design_demo_from(spec.clone(), Topology::EdgeCoupled);
        let hairpin = design_demo_from(spec.clone(), Topology::Hairpin);
        let combline = design_demo_from(spec, Topology::Combline);

        // Interdigital dimensions onto FR-4.
        assert!(
            interdigital.layout.is_some(),
            "demo spec dimensions as an interdigital board"
        );
        assert!(
            interdigital.dim_error.is_none(),
            "no interdigital dimensioning error"
        );
        assert_eq!(interdigital.topology, Topology::Interdigital);
        assert_eq!(
            interdigital.resonators.len(),
            5,
            "order-5 interdigital → 5 resonators"
        );
        assert!(interdigital.board_size_mm.0 > 0.0 && interdigital.board_size_mm.1 > 0.0);

        // The shared coupled-resonator synthesis is topology-INDEPENDENT: same
        // coupling matrix and same verdict as the edge-coupled realization (the
        // interdigital card does NOT change the synthesis — only the realization).
        assert_eq!(
            interdigital.coupling.m, edge.coupling.m,
            "interdigital shares the edge-coupled coupling matrix"
        );
        assert_eq!(
            interdigital.report.pass, edge.report.pass,
            "interdigital shares the edge-coupled verdict"
        );

        // ...but the realized board DIFFERS from edge-coupled, hairpin, AND
        // combline (the card routes to the real interdigital dimensioner — full
        // λg/4 lines short-circuited at alternating ends, distinct from the
        // straight λ/2 edge-coupled lines, the U-folded λ/4 hairpin arms, and the
        // short θ0=λg/8 combline lines). A stub / clone of any of them would make
        // these EQUAL.
        let id_layout = interdigital.layout.as_ref().unwrap();
        let e_layout = edge.layout.as_ref().unwrap();
        let h_layout = hairpin.layout.as_ref().unwrap();
        let cb_layout = combline.layout.as_ref().unwrap();
        assert_ne!(
            layout_signature(id_layout),
            layout_signature(e_layout),
            "interdigital layout must differ from the edge-coupled layout (real, not a clone)"
        );
        assert_ne!(
            layout_signature(id_layout),
            layout_signature(h_layout),
            "interdigital layout must differ from the hairpin layout (real, not a clone)"
        );
        assert_ne!(
            layout_signature(id_layout),
            layout_signature(cb_layout),
            "interdigital layout must differ from the combline layout (real, not a combline clone)"
        );

        // The interdigital-distinct point: there is NO loading cap (it is the
        // θ=π/2 self-resonant limit of combline). A combline-with-cap stand-in
        // fails this; combline itself surfaces a Some(C_L).
        assert!(
            interdigital.combline_loading_cap_f.is_none(),
            "interdigital has NO loading cap (full λg/4, alternating shorts)"
        );
        assert!(
            combline.combline_loading_cap_f.is_some(),
            "combline (the cap-bearing sibling) DOES surface a loading cap"
        );

        // The surfaced λg/4 resonator length is real (> 0, finite) — the quantity
        // the interdigital resonator table shows in place of a cap.
        let id_len = interdigital.resonators[0].length_mm;
        assert!(
            id_len > 0.0 && id_len.is_finite(),
            "interdigital surfaces a positive λg/4 resonator length ({id_len:.2} mm)"
        );
        // The λg/4 interdigital line is ≈ half the straight λ/2 edge-coupled
        // resonator — a concrete witness the interdigital dimensioner (not the
        // edge-coupled one) produced the geometry.
        let e_len = edge.resonators[0].length_mm;
        assert!(
            id_len < e_len * 0.7,
            "interdigital λg/4 resonator ({id_len:.2} mm) is well under the edge-coupled λ/2 ({e_len:.2} mm)"
        );
    }

    #[test]
    fn design_from_edited_spec_drives_synthesis() {
        // Editing the order from the Spec form re-drives synthesis: an order-3
        // spec yields a 3×3 coupling matrix + 3-resonator ladder, distinct
        // from the demo's order-5.
        let mut spec = demo_spec();
        spec.order = Some(3);
        let d = design_demo_from(spec.clone(), Topology::EdgeCoupled);
        assert_eq!(d.order(), 3, "edited order flows through synthesize");
        assert_eq!(d.coupling.m.len(), 3);

        let l = design_lumped_from(spec, DEFAULT_Q_UNLOADED).expect("order-3 BPF is realizable");
        assert_eq!(l.order(), 3, "edited order flows through synthesize_lumped");
        assert_eq!(l.resonators.len(), 3);
    }

    #[test]
    fn unrealizable_spec_degrades_gracefully() {
        // A wide fractional bandwidth at a low order over-couples the edge-
        // coupled gaps beyond what FR-4 can realize: the geometry should
        // degrade (no layout / no resonator rows, a dim_error) while the
        // synthesized prototype + response stay real, rather than panicking.
        let mut spec = demo_spec();
        spec.order = Some(2);
        spec.fbw = 0.6;
        let d = design_demo_from(spec, Topology::EdgeCoupled);
        // Synthesis is still real.
        assert_eq!(d.order(), 2);
        assert_eq!(d.sweep.len(), SWEEP_POINTS);
        // Geometry degraded coherently.
        if d.dim_error.is_some() {
            assert!(d.layout.is_none());
            assert!(d.resonators.is_empty());
            assert_eq!(d.board_size_mm, (0.0, 0.0));
        } else {
            // If FR-4 happens to realize it, the layout must be present.
            assert!(d.layout.is_some());
        }
    }

    /// Nearest-sample `|S21|` in dB at the requested frequency from a sweep.
    fn s21_db_at(sweep: &[SweepPoint], f_hz: f64) -> f64 {
        sweep
            .iter()
            .min_by(|a, b| {
                (a.f_hz - f_hz)
                    .abs()
                    .partial_cmp(&(b.f_hz - f_hz).abs())
                    .unwrap()
            })
            .expect("non-empty sweep")
            .s21_db
    }

    #[test]
    fn stepped_card_routes_to_real_lowpass_engine() {
        // The SteppedImpedance gallery card must drive the REAL F1.2.3
        // `dimension_stepped_impedance` + the App.2.2 `ideal_response_lowpass`,
        // not a stub or the band-pass stand-in. A Butterworth low-pass demo spec
        // must:
        //   (a) produce ≥ order real stepped sections, low-Z first, from the
        //       real dimensioner (positive, finite, alternating widths); and
        //   (b) sweep a low-pass |S21| that is ≈ −3 dB at the cutoff and rolls
        //       off into the stopband — proving the response is genuinely
        //       low-pass (a band-pass clone would peak at f_c, not be −3 dB).
        let mut spec = stepped_demo_spec();
        spec.approximation = Approximation::Butterworth; // clean −3 dB cutoff
        spec.order = Some(5);
        let d = design_stepped_from(spec);

        // (a) Real stepped sections from the real dimensioner.
        assert_eq!(d.order, 5, "order-5 low-pass prototype");
        assert!(
            d.dim_error.is_none(),
            "demo low-pass spec dimensions on FR-4: {:?}",
            d.dim_error
        );
        assert!(d.layout.is_some(), "stepped board is realized");
        assert!(
            d.sections.len() >= d.order,
            "≥ order sections (got {} for order {})",
            d.sections.len(),
            d.order
        );
        assert_eq!(d.sections.len(), 5, "order-5 prototype → 5 line sections");
        // Low-Z first: the standard low-pass prototype starts with a shunt
        // capacitor (low-Z line), then alternates.
        assert!(
            !d.sections[0].high_z,
            "section 1 is the low-Z (shunt-C) line"
        );
        for (i, s) in d.sections.iter().enumerate() {
            assert_eq!(s.high_z, i % 2 == 1, "section {} alternation", i + 1);
            assert!(s.z_ohm > 0.0 && s.width_mm > 0.0 && s.length_mm > 0.0);
            // High-Z sections use Z_high (120 Ω), low-Z use Z_low (20 Ω).
            let expect_z = if s.high_z { d.z_high() } else { d.z_low() };
            assert!(
                (s.z_ohm - expect_z).abs() < 1e-9,
                "section impedance pairing"
            );
        }
        // Non-vacuous: not all section lengths equal (a constant synthesizer
        // would emit identical sections).
        let first_len = d.sections[0].length_mm;
        assert!(
            d.sections
                .iter()
                .any(|s| (s.length_mm - first_len).abs() > 1e-6),
            "section lengths must differ (real, not a constant stub)"
        );
        assert!(d.board_size_mm.0 > 0.0 && d.board_size_mm.1 > 0.0);

        // (b) The swept response is a real low-pass: ≈ −3 dB at the cutoff and
        // far down in the stopband — the band-pass `ideal_response` would NOT
        // pass through a half-power cutoff at f0 (it peaks there).
        assert_eq!(d.sweep.len(), SWEEP_POINTS);
        let f_c = d.cutoff_hz();
        let cutoff_db = s21_db_at(&d.sweep, f_c);
        assert!(
            (cutoff_db - (-3.0103)).abs() <= 0.3,
            "Butterworth |S21(f_c)| = {cutoff_db:.3} dB, expected ≈ −3.01 dB (real low-pass)"
        );
        // DC passband is near 0 dB; the stopband (2 f_c) is well below −3 dB.
        let dc_db = s21_db_at(&d.sweep, f_c * 1e-3);
        assert!(
            dc_db > -0.5,
            "near-DC passband |S21| ≈ 0 dB (got {dc_db:.3})"
        );
        let stop_db = s21_db_at(&d.sweep, 2.0 * f_c);
        assert!(
            stop_db < cutoff_db - 10.0,
            "stopband |S21|({stop_db:.2} dB) rolls off well past the cutoff"
        );

        // The low-pass spec carries the cutoff in f0_hz with Response::Lowpass.
        assert_eq!(d.spec.response, Response::Lowpass);
    }

    #[test]
    fn stepped_design_boots() {
        // The boot-state convenience wrapper produces a real, realizable design.
        let d = design_stepped();
        assert_eq!(d.spec.response, Response::Lowpass);
        assert!(d.order >= 1);
        assert!(!d.sections.is_empty());
        assert_eq!(d.sweep.len(), SWEEP_POINTS);
        assert!(d.sweep.iter().all(|s| s.s21_db.is_finite()));
    }

    #[test]
    fn topbar_view_is_topology_aware() {
        // The TopBar view must reflect the ACTIVE flow, not always the band-pass
        // `designed` (the App.2.3 honesty fix). Build all three flows from the
        // real engine and assert each branch's summary shape + verdict source —
        // a band-pass-only / constant `topbar_view` fails every assertion below.
        let designed = design_demo();
        let lumped = design_lumped();
        let stepped = design_stepped();

        // (a) Stepped-impedance is a LOW-PASS flow: a cutoff, no fractional
        // bandwidth `%`, and the low-pass verdict.
        let (sp_summary, sp_verdict) = topbar_view(
            Topology::SteppedImpedance,
            &designed,
            Some(&lumped),
            &stepped,
        );
        assert!(
            sp_summary.contains("cutoff"),
            "stepped summary shows a cutoff (got {sp_summary:?})"
        );
        assert!(
            !sp_summary.contains('%'),
            "stepped (low-pass) summary has no fractional bandwidth % (got {sp_summary:?})"
        );
        assert_eq!(
            sp_verdict,
            Some(stepped.pass),
            "stepped verdict is the low-pass verdict"
        );

        // (b) Edge-coupled is a BAND-PASS flow: a fractional bandwidth `%`, no
        // "cutoff", and the distributed verdict.
        let (ec_summary, ec_verdict) =
            topbar_view(Topology::EdgeCoupled, &designed, Some(&lumped), &stepped);
        assert!(
            ec_summary.contains('%'),
            "edge-coupled (band-pass) summary shows fractional bandwidth % (got {ec_summary:?})"
        );
        assert!(
            !ec_summary.contains("cutoff"),
            "band-pass summary has no cutoff (got {ec_summary:?})"
        );
        assert_eq!(
            ec_verdict,
            Some(designed.report.pass),
            "edge-coupled verdict is the distributed verdict"
        );

        // (c) Lumped shares the band-pass summary but reports its OWN ladder
        // verdict — from `lumped`, not `designed`.
        let (lc_summary, lc_verdict) =
            topbar_view(Topology::LumpedLc, &designed, Some(&lumped), &stepped);
        assert!(
            lc_summary.contains('%'),
            "lumped shares the band-pass summary (got {lc_summary:?})"
        );
        assert_eq!(
            lc_verdict,
            Some(&lumped).map(|l| l.verdict.pass),
            "lumped verdict comes from the lumped ladder"
        );
        // An unrealizable lumped ladder (None) → no verdict (not realizable).
        let (_, lc_none) = topbar_view(Topology::LumpedLc, &designed, None, &stepped);
        assert_eq!(lc_none, None, "no lumped design → not-realizable verdict");
    }

    #[test]
    fn verify_view_surfaces_real_per_flow_metrics() {
        // The Verify view must surface the ACTIVE flow's REAL graded metrics —
        // pulled straight from each flow's source struct, not a constant or "—"
        // (the App.2.4 honesty fix). Build all three flows from the real engine
        // and assert each branch's level + metrics equal the source fields. A
        // fake / constant `verify_view` fails every assertion below.
        let designed = design_demo();
        let lumped = design_lumped();
        let stepped = design_stepped();

        // (a) Lumped → the REALIZED-ladder level, metrics == `lumped.verdict`.
        let lc = verify_view(Topology::LumpedLc, &designed, Some(&lumped), &stepped);
        assert_eq!(
            lc.level,
            VerifyLevel::RealizedLadder,
            "lumped grades the realized LC ladder"
        );
        assert_eq!(lc.pass, Some(lumped.verdict.pass), "lumped verdict pass");
        assert_eq!(
            lc.worst_passband_ripple_db, lumped.verdict.worst_passband_ripple_db,
            "lumped ripple is the verdict's"
        );
        assert_eq!(
            lc.worst_return_loss_db, lumped.verdict.worst_return_loss_db,
            "lumped return loss is the verdict's"
        );
        // The demo spec has a stopband point, so the rejection is finite + present
        // and equals the verdict's worst rejection.
        assert!(lumped.verdict.worst_stopband_rej_db.is_finite());
        assert_eq!(
            lc.worst_stopband_rej_db,
            Some(lumped.verdict.worst_stopband_rej_db),
            "lumped stopband rejection is the verdict's"
        );

        // (b) Stepped-impedance → the SYNTHESIZED-ideal level, metrics == stepped.
        let sp = verify_view(
            Topology::SteppedImpedance,
            &designed,
            Some(&lumped),
            &stepped,
        );
        assert_eq!(
            sp.level,
            VerifyLevel::SynthesizedIdeal,
            "stepped grades the synthesized ideal response"
        );
        assert_eq!(sp.pass, Some(stepped.pass), "stepped verdict pass");
        assert_eq!(
            sp.worst_passband_ripple_db, stepped.worst_passband_ripple_db,
            "stepped ripple"
        );
        assert_eq!(
            sp.worst_return_loss_db, stepped.worst_return_loss_db,
            "stepped return loss"
        );
        // Stopband rejection is the minimum achieved over the stepped stopband.
        let sp_min = stepped
            .stopband
            .iter()
            .map(|&(_, a, _, _)| a)
            .fold(f64::INFINITY, f64::min);
        assert!(stepped.stopband.is_empty() == sp.worst_stopband_rej_db.is_none());
        if !stepped.stopband.is_empty() {
            assert_eq!(
                sp.worst_stopband_rej_db,
                Some(sp_min),
                "stepped rejection is the min achieved over its stopband"
            );
        }

        // (c) Edge-coupled → the SYNTHESIZED-ideal level, metrics == report.
        let ec = verify_view(Topology::EdgeCoupled, &designed, Some(&lumped), &stepped);
        assert_eq!(ec.level, VerifyLevel::SynthesizedIdeal);
        assert_eq!(ec.pass, Some(designed.report.pass), "edge-coupled verdict");
        assert_eq!(
            ec.worst_passband_ripple_db, designed.report.worst_passband_ripple_db,
            "edge-coupled ripple is the report's"
        );
        assert_eq!(
            ec.worst_return_loss_db, designed.report.worst_return_loss_db,
            "edge-coupled return loss is the report's"
        );
        let ec_min = designed
            .report
            .stopband
            .iter()
            .map(|&(_, a, _, _)| a)
            .fold(f64::INFINITY, f64::min);
        assert!(
            !designed.report.stopband.is_empty(),
            "demo has a stopband pt"
        );
        assert_eq!(
            ec.worst_stopband_rej_db,
            Some(ec_min),
            "edge-coupled rejection is the min achieved over the report stopband"
        );

        // The LEVEL differs lumped-vs-distributed (a constant `verify_view`
        // would tag both the same).
        assert_ne!(
            lc.level, ec.level,
            "lumped (realized) and distributed (synthesized) levels differ"
        );

        // An unrealizable lumped ladder (None) → no verdict (not realizable),
        // zeroed metrics, no rejection — never a fabricated number.
        let lc_none = verify_view(Topology::LumpedLc, &designed, None, &stepped);
        assert_eq!(lc_none.level, VerifyLevel::RealizedLadder);
        assert_eq!(lc_none.pass, None, "no lumped design → not realizable");
        assert_eq!(lc_none.worst_passband_ripple_db, 0.0);
        assert_eq!(lc_none.worst_return_loss_db, 0.0);
        assert_eq!(lc_none.worst_stopband_rej_db, None);
    }

    #[test]
    fn compare_techniques_surfaces_real_per_technique_metrics() {
        // `compare_techniques` must synthesize EVERY live technique for the
        // response class and surface each one's REAL graded metrics — pulled
        // straight from that technique's freshly-built design, not a constant.
        // A constant / empty `compare_techniques` fails every assertion below.

        // ---- (a) band-pass demo spec → five techniques, real metrics --------
        let rows = compare_techniques(&demo_spec());
        assert_eq!(
            rows.len(),
            5,
            "band-pass realizes edge-coupled + hairpin + combline + interdigital + lumped"
        );
        assert_eq!(rows[0].technique, RealizationTechnique::EdgeCoupled);
        assert_eq!(rows[1].technique, RealizationTechnique::Hairpin);
        assert_eq!(rows[2].technique, RealizationTechnique::Combline);
        assert_eq!(rows[3].technique, RealizationTechnique::Interdigital);
        assert_eq!(rows[4].technique, RealizationTechnique::LumpedLc);

        // Each row's metrics equal that technique's freshly-built design's
        // graded fields (equality, not constants).
        let edge = design_demo_from(demo_spec(), Topology::EdgeCoupled);
        let hairpin = design_demo_from(demo_spec(), Topology::Hairpin);
        let combline = design_demo_from(demo_spec(), Topology::Combline);
        let interdigital = design_demo_from(demo_spec(), Topology::Interdigital);
        let lumped =
            design_lumped_from(demo_spec(), DEFAULT_Q_UNLOADED).expect("demo ladder is realizable");

        let ec = &rows[0];
        assert_eq!(ec.realizable, edge.layout.is_some());
        assert_eq!(ec.board_w_mm, edge.board_size_mm.0);
        assert_eq!(ec.board_h_mm, edge.board_size_mm.1);
        assert_eq!(ec.pass, Some(edge.report.pass));
        assert_eq!(ec.order, edge.order());
        assert_eq!(
            ec.worst_passband_ripple_db,
            edge.report.worst_passband_ripple_db
        );
        assert_eq!(ec.worst_return_loss_db, edge.report.worst_return_loss_db);
        let ec_min = edge
            .report
            .stopband
            .iter()
            .map(|&(_, a, _, _)| a)
            .fold(f64::INFINITY, f64::min);
        assert!(!edge.report.stopband.is_empty(), "demo has a stopband pt");
        assert_eq!(ec.worst_stopband_rej_db, Some(ec_min));

        let hp = &rows[1];
        assert_eq!(hp.realizable, hairpin.layout.is_some());
        assert_eq!(hp.board_w_mm, hairpin.board_size_mm.0);
        assert_eq!(hp.board_h_mm, hairpin.board_size_mm.1);
        assert_eq!(hp.pass, Some(hairpin.report.pass));
        assert_eq!(hp.order, hairpin.order());
        assert_eq!(
            hp.worst_passband_ripple_db,
            hairpin.report.worst_passband_ripple_db
        );
        assert_eq!(hp.worst_return_loss_db, hairpin.report.worst_return_loss_db);

        let cb = &rows[2];
        assert_eq!(cb.realizable, combline.layout.is_some());
        assert_eq!(cb.board_w_mm, combline.board_size_mm.0);
        assert_eq!(cb.board_h_mm, combline.board_size_mm.1);
        assert_eq!(cb.pass, Some(combline.report.pass));
        assert_eq!(cb.order, combline.order());
        assert_eq!(
            cb.worst_passband_ripple_db,
            combline.report.worst_passband_ripple_db
        );
        assert_eq!(
            cb.worst_return_loss_db,
            combline.report.worst_return_loss_db
        );

        let id = &rows[3];
        assert_eq!(id.realizable, interdigital.layout.is_some());
        assert_eq!(id.board_w_mm, interdigital.board_size_mm.0);
        assert_eq!(id.board_h_mm, interdigital.board_size_mm.1);
        assert_eq!(id.pass, Some(interdigital.report.pass));
        assert_eq!(id.order, interdigital.order());
        assert_eq!(
            id.worst_passband_ripple_db,
            interdigital.report.worst_passband_ripple_db
        );
        assert_eq!(
            id.worst_return_loss_db,
            interdigital.report.worst_return_loss_db
        );

        let lc = &rows[4];
        assert!(lc.realizable, "demo ladder is realizable");
        assert_eq!(lc.board_w_mm, lumped.board_size_mm.0);
        assert_eq!(lc.board_h_mm, lumped.board_size_mm.1);
        assert_eq!(lc.pass, Some(lumped.verdict.pass));
        assert_eq!(lc.order, lumped.order());
        assert_eq!(
            lc.worst_passband_ripple_db,
            lumped.verdict.worst_passband_ripple_db
        );
        assert_eq!(lc.worst_return_loss_db, lumped.verdict.worst_return_loss_db);
        assert!(lumped.verdict.worst_stopband_rej_db.is_finite());
        assert_eq!(
            lc.worst_stopband_rej_db,
            Some(lumped.verdict.worst_stopband_rej_db)
        );

        // The rows are NOT all identical: the techniques realize the same
        // prototype on physically different boards. In particular the set of
        // board sizes is not a single value (hairpin folds smaller than the
        // straight edge-coupled lines; lumped is an SMD board entirely).
        let board_sizes = [
            (ec.board_w_mm, ec.board_h_mm),
            (hp.board_w_mm, hp.board_h_mm),
            (cb.board_w_mm, cb.board_h_mm),
            (id.board_w_mm, id.board_h_mm),
            (lc.board_w_mm, lc.board_h_mm),
        ];
        assert!(
            board_sizes.iter().any(|&b| b != board_sizes[0]),
            "the techniques must differ physically — board sizes are not all \
             identical (edge-coupled {:?} vs hairpin {:?} vs combline {:?} vs \
             interdigital {:?} vs lumped {:?})",
            board_sizes[0],
            board_sizes[1],
            board_sizes[2],
            board_sizes[3],
            board_sizes[4]
        );

        // ---- (b) low-pass spec → exactly the stepped-impedance row ----------
        let lp = compare_techniques(&stepped_demo_spec());
        assert_eq!(lp.len(), 1, "low-pass realizes exactly stepped-impedance");
        assert_eq!(lp[0].technique, RealizationTechnique::SteppedImpedance);
        let stepped = design_stepped_from(stepped_demo_spec());
        assert_eq!(lp[0].pass, Some(stepped.pass));
        assert_eq!(lp[0].order, stepped.order);
        assert_eq!(
            lp[0].worst_passband_ripple_db,
            stepped.worst_passband_ripple_db
        );
        assert_eq!(lp[0].worst_return_loss_db, stepped.worst_return_loss_db);
        assert_eq!(lp[0].board_w_mm, stepped.board_size_mm.0);
        assert_eq!(lp[0].board_h_mm, stepped.board_size_mm.1);

        // ---- (c) high-pass spec → no live technique -------------------------
        let mut hpf = demo_spec();
        hpf.response = Response::Highpass;
        assert!(
            compare_techniques(&hpf).is_empty(),
            "no live technique realizes a high-pass response yet"
        );

        // ---- (d) unrealizable distributed design → degraded row -------------
        // A wide fractional bandwidth at a low order is unrealizable on FR-4
        // (the coupling gaps cannot be bracketed), so the edge-coupled + hairpin
        // rows must degrade per the struct contract: realizable=false, pass=None,
        // zeroed board + metrics (the order stays — synthesis, not geometry).
        let mut wide = demo_spec();
        wide.order = Some(2);
        wide.fbw = 0.6;
        // Precondition: the distributed designs are genuinely unrealizable here.
        assert!(
            design_demo_from(wide.clone(), Topology::EdgeCoupled)
                .layout
                .is_none(),
            "fixture must be an unrealizable distributed design"
        );
        let wide_rows = compare_techniques(&wide);
        for t in [
            RealizationTechnique::EdgeCoupled,
            RealizationTechnique::Hairpin,
        ] {
            let row = wide_rows
                .iter()
                .find(|r| r.technique == t)
                .expect("band-pass lists the distributed techniques");
            assert!(!row.realizable, "{t:?} is unrealizable for this fixture");
            assert_eq!(row.pass, None, "{t:?}: no verdict when unrealizable");
            assert_eq!(row.board_w_mm, 0.0, "{t:?}: zeroed board width");
            assert_eq!(row.board_h_mm, 0.0, "{t:?}: zeroed board height");
            assert_eq!(
                row.worst_passband_ripple_db, 0.0,
                "{t:?}: zeroed ripple metric"
            );
            assert_eq!(row.worst_return_loss_db, 0.0, "{t:?}: zeroed RL metric");
            assert_eq!(row.worst_stopband_rej_db, None, "{t:?}: no rejection");
        }
    }

    #[test]
    fn overlay_curves_are_real_and_distinct() {
        // `overlay_curves` must return the genuinely-distinct swept responses
        // for the spec — real engine sweeps on the shared grid, NOT a constant
        // or two relabelled copies of one curve. A constant / empty
        // `overlay_curves` fails every assertion below.

        // ---- (a) band-pass demo → coupled-resonator ideal + lumped realized -
        let curves = overlay_curves(&demo_spec());
        assert_eq!(
            curves.len(),
            2,
            "band-pass overlays the coupled-resonator ideal + the lumped realized ladder"
        );
        let coupled = &curves[0];
        let lumped = &curves[1];
        assert!(coupled.realizable, "the synthesized ideal is realizable");
        assert!(lumped.realizable, "the demo ladder is realizable");
        // Honest labels: the distributed techniques are ONE shared ideal curve
        // naming edge-coupled + hairpin + combline + interdigital (they share the
        // coupled-resonator synthesis); the lumped is a distinct realized curve.
        assert!(
            coupled.label.contains("edge-coupled")
                && coupled.label.contains("hairpin")
                && coupled.label.contains("combline")
                && coupled.label.contains("interdigital"),
            "coupled-resonator curve names edge-coupled + hairpin + combline + interdigital (got {:?})",
            coupled.label
        );
        assert!(
            lumped.label.contains("Lumped"),
            "lumped curve is labelled lumped (got {:?})",
            lumped.label
        );

        // Each curve's sweep is the corresponding design's REAL sweep (equality,
        // not a synthetic copy).
        let coupled_design = design_demo_from(demo_spec(), Topology::EdgeCoupled);
        let lumped_design =
            design_lumped_from(demo_spec(), DEFAULT_Q_UNLOADED).expect("demo ladder is realizable");
        let s21s = |sw: &[SweepPoint]| sw.iter().map(|p| p.s21_db).collect::<Vec<_>>();
        let fhzs = |sw: &[SweepPoint]| sw.iter().map(|p| p.f_hz).collect::<Vec<_>>();
        assert_eq!(
            s21s(&coupled.sweep),
            s21s(&coupled_design.sweep),
            "coupled curve == the distributed design's sweep"
        );
        assert_eq!(
            s21s(&lumped.sweep),
            s21s(&lumped_design.sweep),
            "lumped curve == the lumped design's sweep"
        );

        // Both sweeps are on the SAME frequency grid (same length + matching f).
        assert_eq!(
            coupled.sweep.len(),
            lumped.sweep.len(),
            "both curves are on the same grid length"
        );
        assert_eq!(coupled.sweep.len(), SWEEP_POINTS);
        assert_eq!(
            fhzs(&coupled.sweep),
            fhzs(&lumped.sweep),
            "both curves share the sweep frequency grid"
        );

        // The two responses genuinely DIFFER (the whole point of the overlay):
        // the lumped realized |S21| ≠ the coupled-resonator ideal at ≥1 point.
        assert!(
            coupled
                .sweep
                .iter()
                .zip(lumped.sweep.iter())
                .any(|(c, l)| (c.s21_db - l.s21_db).abs() > 1e-6),
            "lumped realized |S21| differs from the coupled-resonator ideal somewhere"
        );

        // ---- (b) low-pass → exactly the stepped-impedance ideal -------------
        let lp = overlay_curves(&stepped_demo_spec());
        assert_eq!(lp.len(), 1, "low-pass overlays exactly the stepped ideal");
        assert!(lp[0].realizable);
        assert!(
            lp[0].label.to_lowercase().contains("stepped"),
            "low-pass curve is the stepped-impedance ideal (got {:?})",
            lp[0].label
        );
        let stepped_design = design_stepped_from(stepped_demo_spec());
        assert_eq!(
            s21s(&lp[0].sweep),
            s21s(&stepped_design.sweep),
            "stepped curve == the stepped design's sweep"
        );

        // ---- (c) high-pass → no live technique ------------------------------
        let mut hpf = demo_spec();
        hpf.response = Response::Highpass;
        assert!(
            overlay_curves(&hpf).is_empty(),
            "no live technique realizes a high-pass response yet"
        );
    }

    /// A Chebyshev 0.5 dB / N = 3 / Z0 = 50 Ω band-pass spec at the given centre
    /// and fractional bandwidth — the schema the `cli-jlcpcb-autoroute` gate
    /// drives, here at the studio-engine layer. A single stopband row at 1.5·f0
    /// keeps the spec well-formed (the mask verdict does not gate this test).
    fn autoroute_spec(f0_hz: f64, fbw: f64) -> FilterSpec {
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

    #[test]
    fn orderable_upload_auto_routes_topology_and_footprint() {
        // ADR-0169 T5, mirroring `cli-jlcpcb-autoroute` at the studio-engine
        // layer (NON-circular: it asserts the engine's auto-routed
        // `orderable_upload` result — topology, footprint, orderability, blank
        // count — for three specs that route THREE different ways, not a
        // relabelled single outcome). `orderable_upload` searches
        // (topology × footprint) for an orderable JLCPCB realization.

        // ---- wideband 1 GHz / 70 % → alternating ladder, fully orderable -----
        // The conventional ladder is orderable for a wideband spec, so the
        // search returns it (arm 1) without falling back to top-C.
        let wide = synthesize(&autoroute_spec(1.0e9, 0.70));
        let o_wide = orderable_upload(&wide, &SUBSTRATE);
        assert_eq!(
            o_wide.topology_label, "alternating ladder",
            "1 GHz/70% must route to the alternating ladder"
        );
        assert!(
            o_wide.fully_orderable,
            "1 GHz/70% ladder must be fully orderable (zero blanks)"
        );
        assert_eq!(o_wide.n_blank, 0, "fully orderable ⇒ zero blanks");

        // ---- 0.5 GHz / 20 % → top-C-coupled, 0402, fully orderable, 0 blank --
        // THE SHOWCASE (the T3 discriminating cell reached via the studio
        // engine): the alternating ladder blanks here, but top-C is fully
        // orderable on the finest 0402 grid — so the search rescues it.
        let sub = synthesize(&autoroute_spec(0.5e9, 0.20));
        let o_sub = orderable_upload(&sub, &SUBSTRATE);
        assert_eq!(
            o_sub.topology_label, "top-C-coupled",
            "0.5 GHz/20% must route to the top-C-coupled topology"
        );
        assert_eq!(
            o_sub.footprint_label, "0402",
            "0.5 GHz/20% top-C must be orderable on the 0402 footprint"
        );
        assert!(
            o_sub.fully_orderable,
            "0.5 GHz/20% top-C@0402 must be fully orderable"
        );
        assert_eq!(
            o_sub.n_blank, 0,
            "0.5 GHz/20% top-C@0402 must have ZERO blanks (the headline rescue)"
        );

        // ---- 2 GHz / 5 % → NEITHER fully orderable: honest blanks ------------
        // GHz-narrow: neither lumped topology fully orders on any footprint, so
        // the search returns the honest fewest-blanks board with real blanks.
        let narrow = synthesize(&autoroute_spec(2.0e9, 0.05));
        let o_narrow = orderable_upload(&narrow, &SUBSTRATE);
        assert!(
            !o_narrow.fully_orderable,
            "2 GHz/5% must NOT be fully orderable (honest coverage hole)"
        );
        assert!(
            o_narrow.n_blank > 0,
            "2 GHz/5% must carry ≥1 real blank LCSC part (flagged, not dropped)"
        );

        // ---- all three emit non-empty, multi-line BOM + CPL CSVs -------------
        for (tag, o) in [
            ("wideband", &o_wide),
            ("subghz", &o_sub),
            ("ghznarrow", &o_narrow),
        ] {
            assert!(
                o.n_parts > 0,
                "[{tag}] auto-routed board must place ≥1 part"
            );
            assert!(
                !o.bom_csv.is_empty() && o.bom_csv.lines().count() > 1,
                "[{tag}] JLCPCB BOM CSV must have a header + ≥1 data row"
            );
            assert!(
                !o.cpl_csv.is_empty() && o.cpl_csv.lines().count() > 1,
                "[{tag}] JLCPCB CPL CSV must have a header + ≥1 data row"
            );
        }
    }
}
