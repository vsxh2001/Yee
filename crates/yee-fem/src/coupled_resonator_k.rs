//! Coupled-microstrip resonator-pair coupling-coefficient `k` extraction from
//! a frequency-domain FEM `S21` driven sweep (FEM-EM brick K1, ADR-0155).
//!
//! This is the **production promotion** of the de-risk probe
//! (`spike/fem-coupled-k-probe`, `933940f`,
//! `tests/coupled_k_probe.rs`) into a reusable `src` API, mirroring how
//! ADR-0154 N1 promoted [`crate::microstrip_port_numerical`] out of its probe.
//! The probe proved (GO) that a frequency-domain FEM `sweep_matrix` of a
//! coupled-microstrip resonator PAIR shows **two resolvable transmission
//! peaks** whose frequency-split yields a coupling coefficient `k` within the
//! ≲30 % de-risk band of the analytic reference.
//!
//! ## Why a FEM driven sweep (not FDTD-resonant)
//!
//! The filter app needs to EM-extract the inter-resonator coupling `k`. The
//! FDTD resonant coupled-resonator split was abandoned (ADR-0108: *no box is
//! simultaneously high-Q and non-confining* — a small PEC box confines the
//! fringing fields that set the split, a large PEC box rings as a cavity that
//! swamps the spectrum, an open CPML box collapses the resonator Q so there are
//! no detectable peaks). The FEM frequency-domain driven sweep is inherently
//! **wall-free** — one linear solve per frequency, no time-stepping, no cavity
//! ring-up — and `k` depends only on peak *locations*, which are robust to the
//! absolute |S21| floor that limits FEM transmission (exactly as the B4 ε_eff
//! was; see [`crate::microstrip_port_numerical`]).
//!
//! ## Two definitions of k (both reported for traceability)
//!
//! There are two distinct, both-legitimate "coupling coefficients" in play,
//! and they are NOT the same number for full-length parallel coupling:
//!
//! * **`k_fem` = (f_hi² − f_lo²)/(f_hi² + f_lo²)** — the *resonator* coupling
//!   from the split of the two measured FEM transmission peaks (the synthesis-
//!   side designer's extract; this is what [`CoupledKResult::k_fem`] holds).
//! * **`k_imp` = (Z0e − Z0o)/(Z0e + Z0o)** — the coupled-*line* voltage
//!   coupling ([`yee_layout::coupling_coefficient`], Kirschning-Jansen). The
//!   ADR-0155 gate's named analytic reference.
//! * **`k_eps`** — the even/odd **ε_eff** resonant split a fixed-length λ_g/2
//!   pair physically reproduces: `f_{e,o} = c/(2L√ε_eff,{e,o})`, so
//!   `k_eps = (f_e² − f_o²)/(f_e² + f_o²)`. This is the like-for-like
//!   resonant-split reference `k_fem` tracks.
//!
//! For full-length parallel coupling `k_imp / k_eps` ranges from ≈3.3 (tight
//! gap S/W=0.3) down to ≈1.1 (weak gap S/W=2): they converge only in the weak
//! limit. The default geometry uses a **weak gap (S = 2·W)** where the two
//! definitions agree to ≈11 %, so the gate's `k_fem`-vs-`coupling_coefficient`
//! comparison is physically meaningful and the coupling is weak enough that the
//! two peaks do not smear together. [`CoupledKResult`] reports BOTH analytic
//! references.
//!
//! ## Geometry (matches the shipped FEM line/filter work — FR-4, h = 1 mm)
//!
//! ```text
//!   substrate h   = sub_h     ε_r = eps_r (FR-4 default 4.4)
//!   strip W       = trace_w   (W/h = 1, ~the B4 line)
//!   gap S         = gap_s     (S/W = 2, weak coupling: k_imp ≈ k_eps)
//!   resonator L   = λ_g/2 at f0 using single-line ε_eff
//!   box           = box_w × ~box_len × box_h  (walls ≥ 2.5·h clear in x — B4
//!                   box-loading finding; box_h = 6 mm = open-half-space air)
//!   feed coupling = 1 dy-cell end-gap tapping each resonator open end
//! ```
//!
//! Axes match [`crate::layered_microstrip_filter_mesh`] (B2/B7): `x`
//! cross-section, `y` propagation (feed-to-feed), `z` substrate-normal (ground
//! z=0, trace z=sub_h). Ports on the `y=0` / `y=box_len` end-caps.
//!
//! ## Excitation — weakly-coupled feeds
//!
//! TWO ports, each WEAKLY coupled to ONE resonator (weak coupling is essential
//! — over-coupling smears the split). A straight feed line runs from each port
//! plane up to a small **end-gap** before its resonator's open end; the gap
//! turns the feed into a weak capacitive tap. The feeds are off-centre in `x`
//! (each aligned with its own resonator), so each numerical-eigenmode wave-port
//! is RECENTRED on its feed's `x` via [`crate::microstrip_port_numerical_at`]
//! (the filter test's per-feed recentre, ADR-0154 N3). `with_coupled_whitney(true)`
//! is MANDATORY (B4 finding: the lumped-centroid port collapses the absorbing
//! block for the substrate-normal `E_z` mode).
//!
//! ## Layering: `yee-filter` peak-finding
//!
//! Peak-finding reuses the SHIPPED [`yee_filter::extract_coupling`] primitive
//! (interior local maxima → two strongest → `k = (f_hi²−f_lo²)/(f_hi²+f_lo²)`)
//! rather than re-deriving it inline. `yee-filter` is therefore a real
//! (non-dev) dependency of `yee-fem`; this is **acyclic** — `yee-filter`
//! depends only on `yee-synth` / `yee-layout` (pure, serde-only), with no edge
//! back to `yee-fem`. Reusing the validated primitive keeps the split formula
//! single-sourced (the same `extract_coupling` the F1.1b DSP gate validates).
//!
//! ## Cost — heavy; gate is `#[ignore]`'d + `--release`
//!
//! [`coupled_resonator_k`] runs a multi-minute driven SWEEP (one per-ω sparse
//! LU per frequency point; the probe was ~280 s at ~63 k tets). Callers that
//! drive it (the `fem-coupling-001` gate) must be `#[ignore]`'d so the debug
//! `cargo test --workspace` never runs them, and run only in `--release`,
//! boxed.

use std::collections::HashMap;
use std::f64::consts::PI;

use nalgebra::Vector3;
use yee_core::Error;
use yee_layout::{coupled_microstrip, coupling_coefficient, eps_eff};
use yee_mesh::TetMesh3D;

use crate::material::MaterialDatabase;
use crate::microstrip_mesh::{TraceRect, layered_microstrip_filter_mesh};
use crate::microstrip_port_numerical::{MicrostripPortGeom, microstrip_port_numerical_at};
use crate::open_boundary::{FaceKind, OpenBoundarySolver};

/// Speed of light (m/s).
const C0: f64 = 299_792_458.0;

/// PEC-shield clearance each side in `x` (m). B4: walls ≥ 2.5·h clear or they
/// load the line / pull ε_eff below the open-microstrip Hammerstad-Jensen value.
const CLEARANCE_X: f64 = 2.5e-3;
/// Straight feed-line length at each end (m), before the coupling gap.
const FEED_RUN: f64 = 6.0e-3;
/// Feed-to-resonator end-gap (m): one `dy` cell, weak capacitive tap.
const FEED_GAP: f64 = 1.0e-3;

/// Cross-section pitch (m): dx = 0.5 mm → W = 2 cells, S = 4 cells at the
/// default geometry.
const DX: f64 = 0.5e-3;
/// Propagation pitch (m): dy = 1.0 mm → resolves the 1 mm feed gap (1 cell).
const DY: f64 = 1.0e-3;
/// Substrate-normal pitch (m): dz = 0.5 mm → 2 substrate z-cells at h = 1 mm.
const DZ: f64 = 0.5e-3;

/// Geometry of a coupled-microstrip resonator pair for FEM `k` extraction.
///
/// Two identical parallel λ_g/2 strips along the propagation axis, separated by
/// the edge-to-edge gap [`gap_s`](Self::gap_s) in the cross-section axis, each
/// weakly capacitively tapped by an off-centre straight feed. The resonator
/// length is computed internally as λ_g/2 at [`f0_hz`](Self::f0_hz) using the
/// single-line Hammerstad-Jensen ε_eff of [`trace_w`](Self::trace_w) on
/// [`sub_h`](Self::sub_h) / [`eps_r`](Self::eps_r).
///
/// The default [`CoupledResonatorGeom::probe`] is the ADR-0155 de-risk-probe
/// geometry (W = 1 mm, S = 2 mm, h = 1 mm, ε_r = 4.4, f0 = 2.4 GHz, 6 mm-tall
/// open box) that measured `k_fem = 0.0481`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledResonatorGeom {
    /// Strip width `W` of each of the two equal resonators, metres.
    pub trace_w: f64,
    /// Edge-to-edge coupling gap `S` between the two resonators, metres.
    pub gap_s: f64,
    /// Substrate height `h`, metres: dielectric fills `z ∈ [0, sub_h]`.
    pub sub_h: f64,
    /// Substrate relative permittivity `ε_r` (e.g. 4.4 for FR-4).
    pub eps_r: f64,
    /// Resonator centre frequency `f0`, Hz. The strip length is set to λ_g/2 at
    /// this frequency (single-line Hammerstad-Jensen ε_eff).
    pub f0_hz: f64,
    /// Box / substrate width along `x`, metres (the cross-section the wave-port
    /// eigenmode is solved on). Must clear the trace pattern + clearance.
    pub box_w: f64,
    /// Box height along `z`, metres: substrate + air (≈6 mm = open half-space).
    pub box_h: f64,
}

impl CoupledResonatorGeom {
    /// The ADR-0155 de-risk-probe geometry: W = 1 mm, S = 2 mm, h = 1 mm,
    /// ε_r = 4.4 (FR-4), f0 = 2.4 GHz, in a 6 mm-tall open box whose width
    /// clears the two-strip pattern by [`CLEARANCE_X`] each side. This is the
    /// geometry the probe measured `k_fem = 0.0481` on (vs `k_imp = 0.0646`,
    /// `k_eps = 0.0581`).
    ///
    /// `box_w` is derived from the strip pattern (two `W`-wide strips with gap
    /// `S`, plus `CLEARANCE_X` each side); `box_h` is `sub_h + 5 mm` air snapped
    /// so `sub_h` lands on a z-plane.
    pub fn probe() -> Self {
        Self::probe_with_gap(2.0e-3)
    }

    /// The ADR-0155 probe geometry at an arbitrary coupling gap `gap_s` (metres),
    /// keeping every other probe parameter fixed (W = 1 mm, h = 1 mm, ε_r = 4.4,
    /// f0 = 2.4 GHz, 6 mm-tall open box). [`probe`](Self::probe) is exactly
    /// `probe_with_gap(2.0e-3)`.
    ///
    /// `box_w` is re-derived for the new gap (two `W`-wide strips with gap `S`,
    /// plus [`CLEARANCE_X`] each side, so the PEC walls stay ≥ 2.5·h clear and
    /// do not load the line — B4); `box_h` is `sub_h + 5 mm` air snapped so
    /// `sub_h` lands on a z-plane (the snap is gap-independent). This is the
    /// constructor the k-vs-gap monotonicity gate (`fem-coupling-002`, ADR-0155
    /// K2) sweeps, so the box-extent derivation stays single-sourced rather than
    /// duplicated per gap.
    pub fn probe_with_gap(gap_s: f64) -> Self {
        let trace_w = 1.0e-3;
        let sub_h = 1.0e-3;
        let air_h = 5.0e-3;
        // x: two strips with a gap S, plus clearance both sides.
        let box_w = CLEARANCE_X + trace_w + gap_s + trace_w + CLEARANCE_X;
        // z: snap so sub_h lands on a z-plane (1 mm / 0.5 mm = 2 cells).
        let nz = ((sub_h + air_h) / DZ).round() as usize;
        let nz_sub = (sub_h / DZ).round().max(1.0) as usize;
        let dz_exact = sub_h / nz_sub as f64;
        let box_h = dz_exact * nz as f64;
        Self {
            trace_w,
            gap_s,
            sub_h,
            eps_r: 4.4,
            f0_hz: 2.4e9,
            box_w,
            box_h,
        }
    }

    /// Resonator length λ_g/2 at [`f0_hz`](Self::f0_hz) using the single-line
    /// Hammerstad-Jensen ε_eff (the resonator is dominantly a single line; the
    /// coupling perturbs the resonance, which is exactly what the FEM sweep
    /// measures).
    pub fn resonator_length_m(&self) -> f64 {
        let eeff = eps_eff(self.trace_w, self.sub_h, self.eps_r);
        let lam_g = C0 / self.f0_hz / eeff.sqrt();
        lam_g / 2.0
    }
}

/// Resolved coupled-pair geometry in mesh world coordinates (the internal
/// resolution of a [`CoupledResonatorGeom`]: traces + box extents + the chosen
/// subdivision + the two feed x-centres).
#[derive(Clone, Debug)]
struct ResolvedGeometry {
    box_w: f64,
    box_len: f64,
    box_h: f64,
    sub_h: f64,
    eps_r: f64,
    f0_hz: f64,
    traces: Vec<TraceRect>,
    nx: usize,
    ny: usize,
    nz: usize,
    /// Resonator length L = λ_g/2 (m).
    res_l: f64,
    /// Strip width = the feed wave-port width (m).
    line_w: f64,
    /// x-centre (m) of the INPUT feed (port 0, y=0 end-cap).
    feed_xc_in: f64,
    /// x-centre (m) of the OUTPUT feed (port 1, y=box_len end-cap).
    feed_xc_out: f64,
}

impl ResolvedGeometry {
    /// Total tetrahedra (`nx·ny·nz·6`) — used by the mesh-size unit test to
    /// catch a resolution drift. `#[cfg(test)]` because the production path
    /// never needs the count (it only drives the solve).
    #[cfg(test)]
    fn total_tets(&self) -> usize {
        self.nx * self.ny * self.nz * 6
    }
}

/// Resolve a [`CoupledResonatorGeom`] into mesh world coordinates: two
/// identical parallel λ_g/2 strips along `y`, separated by gap `S` in `x`, each
/// open-open (floating). Two straight feeds run from the `y=0` / `y=box_len`
/// port planes up to a `FEED_GAP` end-gap before the near open end of resonator
/// A (input) / resonator B (output).
///
/// y-layout: `FEED_RUN | FEED_GAP | L (resonators) | FEED_GAP | FEED_RUN`.
fn resolve_geometry(geom: &CoupledResonatorGeom) -> ResolvedGeometry {
    let CoupledResonatorGeom {
        trace_w: w,
        gap_s: s,
        sub_h,
        eps_r,
        f0_hz,
        box_w,
        box_h,
    } = *geom;
    let res_l = geom.resonator_length_m();

    // x: two strips with a gap S, plus clearance both sides.
    //   resonator A: x ∈ [CLEARANCE_X,           CLEARANCE_X + W]
    //   resonator B: x ∈ [CLEARANCE_X + W + S,    CLEARANCE_X + 2W + S]
    let res_a_x0 = CLEARANCE_X;
    let res_b_x0 = CLEARANCE_X + w + s;

    // y: feeds + gaps + resonators.
    let y_res_lo = FEED_RUN + FEED_GAP;
    let y_res_hi = y_res_lo + res_l;
    let box_len = y_res_hi + FEED_GAP + FEED_RUN;

    // z subdivisions: snap so sub_h lands on a z-plane (the caller's box_h is
    // already snapped by `CoupledResonatorGeom::probe`, but re-derive nz here so
    // a hand-built geom still maps cleanly).
    let nz = (box_h / DZ).round().max(1.0) as usize;

    // Traces: two resonators + two feeds.
    let mut traces = Vec::with_capacity(4);
    // Resonator A (coupled to the input feed).
    traces.push(TraceRect::new(res_a_x0, y_res_lo, w, res_l));
    // Resonator B (coupled to the output feed).
    traces.push(TraceRect::new(res_b_x0, y_res_lo, w, res_l));
    // Input feed: y ∈ [0, FEED_RUN], aligned with resonator A. The FEED_GAP
    // separates the feed top (y = FEED_RUN) from resonator A's bottom open end
    // (y = y_res_lo = FEED_RUN + FEED_GAP).
    traces.push(TraceRect::new(res_a_x0, 0.0, w, FEED_RUN));
    // Output feed: y ∈ [box_len − FEED_RUN, box_len], aligned with resonator B.
    let out_feed_y0 = box_len - FEED_RUN;
    traces.push(TraceRect::new(res_b_x0, out_feed_y0, w, FEED_RUN));

    let nx = (box_w / DX).round().max((w / DX).ceil().max(1.0)) as usize;
    let ny = (box_len / DY).round().max(1.0) as usize;

    ResolvedGeometry {
        box_w,
        box_len,
        box_h,
        sub_h,
        eps_r,
        f0_hz,
        traces,
        nx,
        ny,
        nz,
        res_l,
        line_w: w,
        feed_xc_in: res_a_x0 + w / 2.0,
        feed_xc_out: res_b_x0 + w / 2.0,
    }
}

/// Result of a coupled-resonator FEM `k` extraction.
///
/// `k_fem` is the measured resonant-split coupling; `k_imp_ref` and `k_eps_ref`
/// are the two analytic references (see the [module docs](self)). The full
/// swept `s21_db` curve is retained for auditing (peak + valley locations).
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledKResult {
    /// Lower transmission-peak frequency (Hz), `f_lo < f_hi`.
    pub f_lo_hz: f64,
    /// Upper transmission-peak frequency (Hz), `f_hi > f_lo`.
    pub f_hi_hz: f64,
    /// Measured coupling `k_fem = (f_hi² − f_lo²)/(f_hi² + f_lo²)`.
    pub k_fem: f64,
    /// Whether two cleanly resolvable peaks were found (a finite valley
    /// strictly between them sitting below both peaks). When `false`, the peak
    /// frequencies / `k_fem` are best-effort (or NaN if fewer than two maxima
    /// exist) and should not be trusted.
    pub peaks_resolvable: bool,
    /// |S21| (dB) of the valley (minimum strictly between the two peaks). NaN
    /// if fewer than two peaks were found.
    pub valley_db: f64,
    /// |S21| (dB) of the lower-frequency peak.
    pub peak_lo_db: f64,
    /// |S21| (dB) of the upper-frequency peak.
    pub peak_hi_db: f64,
    /// Analytic coupled-line voltage coupling `k_imp = (Z0e−Z0o)/(Z0e+Z0o)`
    /// ([`yee_layout::coupling_coefficient`]) — the synthesis-side reference
    /// the ADR-0155 gate grades `k_fem` against.
    pub k_imp_ref: f64,
    /// Analytic even/odd ε_eff resonant-split reference
    /// `k_eps = (f_e²−f_o²)/(f_e²+f_o²)` — the like-for-like split `k_fem`
    /// physically tracks (reported for traceability).
    pub k_eps_ref: f64,
    /// The full swept |S21|(f) curve as `(f_GHz, dB)` pairs (for auditing the
    /// peak / valley structure).
    pub s21_db: Vec<(f64, f64)>,
}

fn exterior_face_count(mesh: &TetMesh3D) -> usize {
    let mut face_map: HashMap<[usize; 3], usize> = HashMap::new();
    const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
    for tet in &mesh.tetrahedra {
        for &[a, b, c] in TET_FACES.iter() {
            let mut key = [tet[a], tet[b], tet[c]];
            key.sort_unstable();
            *face_map.entry(key).or_insert(0) += 1;
        }
    }
    face_map.values().filter(|&&c| c == 1).count()
}

/// Classify exterior faces: ports on the `y=0` / `y=box_len` end-caps, all else
/// PEC (identical to the filter test; the interior just has a different
/// footprint).
fn classify_faces(centroids: &[Vector3<f64>], box_len: f64) -> Vec<FaceKind> {
    let tol = 1e-9;
    centroids
        .iter()
        .map(|c| {
            if c.y < tol {
                FaceKind::WavePort(0)
            } else if (c.y - box_len).abs() < tol {
                FaceKind::WavePort(1)
            } else {
                FaceKind::Pec
            }
        })
        .collect()
}

fn db(mag: f64) -> f64 {
    20.0 * mag.log10()
}

/// Extract the coupled-resonator coupling coefficient `k` from a FEM
/// frequency-domain `S21` driven sweep.
///
/// Builds the two-λ_g/2-resonator mesh of `geom`
/// ([`crate::layered_microstrip_filter_mesh`] + [`crate::TraceRect`]), attaches
/// two weakly gap-coupled numerical-eigenmode wave-port feeds (each recentred on
/// its off-centre feed via [`crate::microstrip_port_numerical_at`], with
/// `with_coupled_whitney(true)` — both mandatory, B4 finding), drives
/// `sweep_matrix` over a band straddling the even+odd resonances at `n_pts`
/// frequency points, finds the two strongest |S21| peaks (via the shipped
/// [`yee_filter::extract_coupling`]) plus the valley between them, and returns a
/// [`CoupledKResult`] with `k_fem = (f_hi²−f_lo²)/(f_hi²+f_lo²)` alongside the
/// two analytic references.
///
/// The sweep band is fixed at 2.10–2.70 GHz (centred on the probe's f0 =
/// 2.4 GHz, spanning both analytic resonances with valley + shoulders). `n_pts`
/// sets the resolution across the split; the probe used 61 (10 MHz step). A
/// `n_pts < 5` is clamped up to 5 so the peak-finder has interior points.
///
/// # Cost
///
/// HEAVY: one per-ω sparse LU per frequency point — multi-minute at the default
/// ~63 k-tet probe geometry (the probe measured ~280 s for 61 points). Callers
/// must `#[ignore]` + `--release` + box this; never run it in the debug
/// workspace test.
///
/// # Errors
///
/// Returns a [`yee_core::Error`] if the mesh build, port construction, or the
/// driven `sweep_matrix` fails (e.g. a degenerate geometry whose trace
/// footprint selects no interior edge, or a port whose cross-section eigensolve
/// fails to surface a propagating mode).
pub fn coupled_resonator_k(
    geom: &CoupledResonatorGeom,
    n_pts: usize,
) -> Result<CoupledKResult, Error> {
    let resolved = resolve_geometry(geom);

    // ---- Analytic references --------------------------------------------
    let cm = coupled_microstrip(geom.trace_w, geom.gap_s, geom.sub_h, geom.eps_r);
    let k_imp_ref = coupling_coefficient(&cm); // (Z0e−Z0o)/(Z0e+Z0o)
    // For two parallel λ_g/2 resonators of fixed length L, the even/odd modes
    // resonate at f_e = c/(2L√ε_eff,e), f_o = c/(2L√ε_eff,o). This is the split
    // k_fem should physically reproduce (the resonant-split reference), distinct
    // from the coupled-line k_imp.
    let f_even = C0 / (2.0 * resolved.res_l * cm.eps_eff_e.sqrt());
    let f_odd = C0 / (2.0 * resolved.res_l * cm.eps_eff_o.sqrt());
    let (f_e_lo, f_e_hi) = (f_even.min(f_odd), f_even.max(f_odd));
    let k_eps_ref = (f_e_hi * f_e_hi - f_e_lo * f_e_lo) / (f_e_hi * f_e_hi + f_e_lo * f_e_lo);

    // ---- Sweep band: straddle f_even / f_odd with margin ----------------
    // Centre the band on 2.4 GHz and span generously past both analytic
    // resonances so the two FEM peaks (which may shift low from mesh dispersion
    // + feed-gap loading) are captured with valley + shoulders.
    let n_pts = n_pts.max(5);
    let f_lo_band = 2.10e9;
    let f_hi_band = 2.70e9;
    let freqs_hz: Vec<f64> = (0..n_pts)
        .map(|i| f_lo_band + (f_hi_band - f_lo_band) * (i as f64) / ((n_pts - 1) as f64))
        .collect();
    let omegas: Vec<f64> = freqs_hz.iter().map(|f| 2.0 * PI * f).collect();

    // ---- Drive the sweep -------------------------------------------------
    let sweep = solve_coupled(&resolved, &omegas)?;

    // ---- Extract |S21|(f). No de-embed needed: k depends only on peak
    //      LOCATIONS, which a unit-magnitude feed-phase rotation cannot move. --
    let mut mags: Vec<f64> = Vec::with_capacity(n_pts);
    let mut s21_db: Vec<(f64, f64)> = Vec::with_capacity(n_pts);
    let mut all_finite = true;
    for (i, _omega) in omegas.iter().enumerate() {
        let mag = sweep.s[i][(1, 0)].norm();
        all_finite &= mag.is_finite();
        let f_ghz = freqs_hz[i] / 1e9;
        mags.push(mag);
        s21_db.push((f_ghz, db(mag)));
    }
    if !all_finite {
        return Err(Error::Invalid(
            "coupled_resonator_k: a non-finite |S21| in the driven sweep — the solve broke \
             (port collapsed or mesh invalid)"
                .to_string(),
        ));
    }

    // ---- Find the two transmission peaks (reuse the shipped peak-split
    //      primitive: interior local maxima → two strongest → k split). -----
    let extraction = yee_filter::extract_coupling(&freqs_hz, &mags);

    // Default to the no-resolvable-split case; fill in from the extraction.
    let mut peaks_resolvable = false;
    let mut k_fem = f64::NAN;
    let mut f_lo_hz = f64::NAN;
    let mut f_hi_hz = f64::NAN;
    let mut valley_db = f64::NAN;
    let mut peak_lo_db = f64::NAN;
    let mut peak_hi_db = f64::NAN;

    if let Some(ex) = extraction {
        f_lo_hz = ex.f_lo_hz;
        f_hi_hz = ex.f_hi_hz;
        k_fem = ex.k;

        // Peak magnitudes at the two extracted peak frequencies (look them up in
        // the swept curve by nearest frequency point).
        let mag_at = |f_hz: f64| -> f64 {
            let i = freqs_hz
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| (*a - f_hz).abs().partial_cmp(&(*b - f_hz).abs()).unwrap())
                .map(|(i, _)| i)
                .unwrap_or(0);
            mags[i]
        };
        peak_lo_db = db(mag_at(f_lo_hz));
        peak_hi_db = db(mag_at(f_hi_hz));

        // Valley = minimum |S21| strictly between the two peaks.
        let valley = freqs_hz
            .iter()
            .zip(mags.iter())
            .filter(|(f, _)| **f > f_lo_hz && **f < f_hi_hz)
            .map(|(_, m)| *m)
            .fold(f64::INFINITY, f64::min);
        valley_db = db(valley);

        // Resolvable ⇒ a finite valley below the shallower peak with a real
        // margin. The numeric margin tripwire is the gate's job; here we set the
        // boolean as "a genuine dip exists between two distinct peaks".
        let shallower_db = peak_lo_db.min(peak_hi_db);
        peaks_resolvable = valley.is_finite() && f_hi_hz > f_lo_hz && valley_db < shallower_db;
    }

    Ok(CoupledKResult {
        f_lo_hz,
        f_hi_hz,
        k_fem,
        peaks_resolvable,
        valley_db,
        peak_lo_db,
        peak_hi_db,
        k_imp_ref,
        k_eps_ref,
        s21_db,
    })
}

/// Outcome of a per-gap FEM-k coupling design-curve root-find
/// ([`correct_gap_fem_k`], ADR-0159 B1).
///
/// Reports the corrected gap, the FEM-measured coupling there, the synthesis
/// target it was driven to, how many FEM sweeps it cost, and whether the search
/// hit the requested tolerance. When `converged` is `false`, [`gap_m`](Self::gap_m)
/// / [`k_fem`](Self::k_fem) hold the **best** (closest-to-target) usable eval the
/// bisection saw before exhausting `max_evals` (or finding no usable eval — in
/// which case `k_fem` is NaN).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GapCorrection {
    /// The corrected edge-to-edge coupling gap `S` (metres) that realizes
    /// `k_fem ≈ k_target` — the best gap the bisection found.
    pub gap_m: f64,
    /// The FEM-measured resonant-split coupling [`CoupledKResult::k_fem`] at
    /// [`gap_m`](Self::gap_m). NaN if no usable eval was obtained.
    pub k_fem: f64,
    /// The synthesis-side target coupling the root-find drove toward (a fixed
    /// design constant, NOT a measurement — keeps the gate non-circular).
    pub k_target: f64,
    /// Number of [`coupled_resonator_k`] FEM sweeps spent (each midpoint eval,
    /// including any unusable / nudged ones, counts as one).
    pub n_evals: usize,
    /// `true` iff a usable eval landed within `tol_frac` of `k_target`
    /// (`|k_fem − k_target| / k_target ≤ tol_frac`) inside `max_evals`.
    pub converged: bool,
}

/// Bisection root-find of a **monotone-decreasing** scalar function toward a
/// target value, tolerant of a single noisy / unusable evaluation.
///
/// `eval(x)` returns `Some(y)` for a usable sample (finite) or `None` for an
/// unusable one (a noisy / failed evaluation — e.g. a FEM sweep whose peaks did
/// not resolve). Because `y(x)` is decreasing, the working bracket maintains
/// `y(lo) ≥ target ≥ y(hi)` conceptually: at each midpoint `mid`,
///
/// * `y(mid) > target` ⇒ we are still on the high side ⇒ move the LOW edge up
///   (`lo = mid`) to search the larger-`x` (smaller-`y`) half;
/// * `y(mid) ≤ target` ⇒ move the HIGH edge down (`hi = mid`).
///
/// The loop stops as soon as a usable sample satisfies
/// `|y − target| / |target| ≤ tol_frac` (`converged = true`), or when the eval
/// budget `max_evals` is spent (`converged = false`). It returns the best
/// (closest-to-target) usable `(x, y)` seen, the eval count, and the flag.
///
/// **Outlier robustness (why bisection, not secant):** if `eval(mid)` is `None`
/// (unusable), the bracket is left untouched and the next midpoint is *nudged*
/// a fraction of the bracket width toward `hi`, so a single bad sample neither
/// crashes nor derails the search (secant would divide by a garbage slope). The
/// nudge still consumes one eval against `max_evals`, so an all-`None` curve
/// terminates honestly with `converged = false` rather than looping.
///
/// This is the FEM-free core of [`correct_gap_fem_k`]; the fast unit test drives
/// it with a synthetic decreasing closure (no FEM solve).
fn bisect_monotone_dec(
    mut eval: impl FnMut(f64) -> Option<f64>,
    target: f64,
    mut lo: f64,
    mut hi: f64,
    tol_frac: f64,
    max_evals: usize,
) -> (f64, f64, usize, bool) {
    let mut n_evals = 0usize;
    let mut best_x = f64::NAN;
    let mut best_y = f64::NAN;
    let mut best_err = f64::INFINITY;
    // Fraction of the current bracket to nudge by when a midpoint is unusable.
    let nudge_frac = 0.25;
    // Count consecutive unusable midpoints, to nudge progressively further per
    // miss (a FIXED nudge would re-probe the SAME point on a 2nd consecutive
    // miss, stalling until the budget drains; scaling by the count explores new
    // points each time).
    let mut consecutive_unusable: usize = 0;

    while n_evals < max_evals {
        let span = hi - lo;
        let mut mid = 0.5 * (lo + hi);
        // If recent midpoints were unusable, perturb the probe toward `hi` —
        // progressively further per consecutive miss — so we explore new points
        // instead of re-sampling the same bad one.
        if consecutive_unusable > 0 {
            mid = (mid + nudge_frac * consecutive_unusable as f64 * span).min(hi);
        }

        n_evals += 1;
        match eval(mid) {
            Some(y) if y.is_finite() => {
                consecutive_unusable = 0;
                let err = (y - target).abs() / target.abs();
                if err < best_err {
                    best_err = err;
                    best_x = mid;
                    best_y = y;
                }
                if err <= tol_frac {
                    return (best_x, best_y, n_evals, true);
                }
                // Monotone-DECREASING: y too HIGH ⇒ need a LARGER x (smaller y).
                if y > target {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            _ => {
                // Unusable eval: keep the bracket, nudge the next probe further
                // (counted against the budget).
                consecutive_unusable += 1;
            }
        }
    }

    (best_x, best_y, n_evals, false)
}

/// Correct a coupling gap to realize a target coupling `k` via the FEM resonant
/// split — the EM coupling **design-curve** root-find (Hong-Lancaster full-wave
/// coupling design; ADR-0159 B1, see [[reference-em-in-loop-space-mapping]]).
///
/// The analytic dimensioner sizes coupling gaps from the impedance-k
/// `(Z0e−Z0o)/(Z0e+Z0o)`, which diverges ~37 % from the FEM-realized resonant-k
/// at the gaps a filter actually uses (k_imp ≠ k_eps; ADR-0155 K2). This drives
/// the gap onto the **measured** `K(gap)` curve instead: it bisects
/// `g ↦ coupled_resonator_k({..base, gap_s: g}, n_pts).k_fem` over
/// `[gap_lo, gap_hi]` until the FEM coupling hits `k_target` within `tol_frac`.
///
/// `K(gap)` is **monotone-decreasing** (wider gap → weaker coupling; confirmed
/// smooth + monotone by a 12-solve probe), so the bisection direction is: if the
/// measured `k_fem(mid) > k_target` the coupling is too strong ⇒ the gap is too
/// SMALL ⇒ search the larger-gap half `[mid, hi]`; otherwise search `[lo, mid]`.
/// Bisection (not secant) is deliberate — it tolerates a single noisy / outlier
/// extraction (the F1.2.1.0 run saw one) without dividing by a garbage slope.
///
/// Only `gap_s` is varied; every other field of `base` (including
/// [`box_w`](CoupledResonatorGeom::box_w)) is held fixed, so the strips move
/// inside an otherwise-fixed cross-section. (The `box_w`-from-gap re-derivation
/// lives only in [`CoupledResonatorGeom::probe_with_gap`]; the solve path reads
/// `gap_s` / `box_w` independently — confirmed.)
///
/// # Robustness
///
/// A `coupled_resonator_k` call that returns `Err`, reports
/// `!peaks_resolvable`, or yields a non-finite `k_fem` is treated as an
/// **unusable** sample: the bracket is left intact and the next midpoint is
/// nudged toward `hi` (the eval is still counted against `max_evals`), so one
/// bad sweep neither panics nor derails the search. The returned
/// [`GapCorrection`] always holds the best usable eval seen.
///
/// # Cost
///
/// HEAVY: **one [`coupled_resonator_k`] FEM driven sweep per eval** (multi-minute
/// each — the probe was ~280 s for 61 pts). At `max_evals = 6` this is ~5-6
/// sweeps. Callers must `#[ignore]` + `--release` + box the test that drives it;
/// never run it in the debug workspace test.
pub fn correct_gap_fem_k(
    base: &CoupledResonatorGeom,
    k_target: f64,
    gap_lo: f64,
    gap_hi: f64,
    tol_frac: f64,
    max_evals: usize,
    n_pts: usize,
) -> GapCorrection {
    // Each midpoint gap → one FEM sweep; unusable (Err / !resolvable / non-finite)
    // maps to `None` so the bisection core nudges past it instead of crashing.
    // A per-eval trajectory line is emitted to stderr (visible under
    // `--nocapture`) so the heavy gate can audit the root-find path; it is
    // diagnostic-only and never affects the returned `GapCorrection`.
    let mut eval_idx = 0usize;
    let eval = |g: f64| -> Option<f64> {
        eval_idx += 1;
        let geom = CoupledResonatorGeom { gap_s: g, ..*base };
        let sample = match coupled_resonator_k(&geom, n_pts) {
            Ok(res) if res.peaks_resolvable && res.k_fem.is_finite() => Some(res.k_fem),
            _ => None,
        };
        match sample {
            Some(k) => eprintln!(
                "[correct_gap_fem_k] eval #{eval_idx}: gap={:.4}mm  k_fem={k:.4}  (target={k_target:.4})",
                g * 1e3
            ),
            None => eprintln!(
                "[correct_gap_fem_k] eval #{eval_idx}: gap={:.4}mm  UNUSABLE (err / !peaks_resolvable / non-finite) — nudging",
                g * 1e3
            ),
        }
        sample
    };

    let (gap_m, k_fem, n_evals, converged) =
        bisect_monotone_dec(eval, k_target, gap_lo, gap_hi, tol_frac, max_evals);

    GapCorrection {
        gap_m,
        k_fem,
        k_target,
        n_evals,
        converged,
    }
}

/// Build the two-port driven solver for the coupled pair and run the sweep.
/// Trace (resonators + feeds) AND ground tagged interior-PEC (B1); the two
/// y-end-caps carry the numerical-eigenmode wave-port recentred per feed
/// (ADR-0154 N3) with `with_coupled_whitney(true)` (mandatory, B4).
fn solve_coupled(
    geom: &ResolvedGeometry,
    omegas: &[f64],
) -> Result<crate::open_boundary::SParametersMatrix, Error> {
    let (mesh, material_db, ground_pred, trace_pred) = layered_microstrip_filter_mesh(
        geom.box_w,
        geom.box_len,
        geom.box_h,
        geom.sub_h,
        geom.traces.clone(),
        geom.nx,
        geom.ny,
        geom.nz,
    )?;

    let n_exterior = exterior_face_count(&mesh);
    let picker = OpenBoundarySolver::new(
        &mesh,
        vec![FaceKind::Pec; n_exterior],
        Vec::new(),
        MaterialDatabase::new(),
    )?;

    let ground_edges = picker.interior_edges_matching(&ground_pred);
    let trace_edges = picker.interior_edges_matching(&trace_pred);
    if trace_edges.is_empty() {
        return Err(Error::Invalid(
            "coupled_resonator_k: trace_pred selected no interior edge on the z = sub_h trace \
             footprint — degenerate geometry"
                .to_string(),
        ));
    }
    let mut interior_pec: Vec<usize> = ground_edges;
    interior_pec.extend(trace_edges.iter().copied());
    interior_pec.sort_unstable();
    interior_pec.dedup();
    let centroids = picker.exterior_face_centroids();
    let kinds = classify_faces(&centroids, geom.box_len);
    drop(picker);

    // Numerical-eigenmode wave-port, recentred on each off-centre feed (the
    // feeds sit over their resonators, not at box_w/2). β stays analytic-HJ on
    // the feed width; only the modal SHAPE is numerical. Each face needs its own
    // call (boxed closures are not Clone).
    let port_geom = MicrostripPortGeom {
        trace_w: geom.line_w,
        sub_h: geom.sub_h,
        eps_r: geom.eps_r,
        box_w: geom.box_w,
        box_h: geom.box_h,
    };
    let port_in = microstrip_port_numerical_at(&port_geom, geom.feed_xc_in, geom.f0_hz)?;
    let port_out = microstrip_port_numerical_at(&port_geom, geom.feed_xc_out, geom.f0_hz)?;

    let solver = OpenBoundarySolver::new(&mesh, kinds, vec![port_in, port_out], material_db)?
        .with_interior_pec_edges(interior_pec.iter().copied())
        .with_coupled_whitney(true);

    solver.sweep_matrix(omegas)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The probe coupled-pair geometry is well-formed: 2 resonators + 2 feeds =
    /// 4 rectangles, the box clears the trace pattern in x, sub_h lands on a
    /// z-plane, the feed gap is exactly one dy cell, the feeds reach both port
    /// planes while the resonators stay clear of them, and the two resonators
    /// are gap-separated in x by S. FAST (no FEM solve).
    #[test]
    fn geometry_is_well_formed() {
        let geom = CoupledResonatorGeom::probe();
        let resolved = resolve_geometry(&geom);
        assert_eq!(resolved.traces.len(), 4, "2 resonators + 2 feeds");

        // Box clears the trace x-span.
        let trace_x_hi = resolved
            .traces
            .iter()
            .map(|r| r.x0 + r.w)
            .fold(0.0_f64, f64::max);
        assert!(
            resolved.box_w > trace_x_hi,
            "box_w {:.4} must clear trace x-extent {:.4}",
            resolved.box_w,
            trace_x_hi
        );

        // sub_h on a z-plane.
        let dz = resolved.box_h / resolved.nz as f64;
        let n_sub = geom.sub_h / dz;
        assert!(
            (n_sub - n_sub.round()).abs() < 1e-9,
            "sub_h must land on a z-plane (n_sub = {n_sub})"
        );

        // Feed gap is exactly one dy cell (so it is resolved, not sub-cell).
        let dy = resolved.box_len / resolved.ny as f64;
        assert!(
            (FEED_GAP / dy - 1.0).abs() < 0.05,
            "feed gap ({:.3}mm) should be ~1 dy cell ({:.3}mm) — got {:.2} cells",
            FEED_GAP * 1e3,
            dy * 1e3,
            FEED_GAP / dy,
        );

        // An input feed reaches y=0; an output feed reaches y=box_len.
        let touches_y0 = resolved.traces.iter().any(|r| r.y0.abs() < 1e-12);
        let touches_ylen = resolved
            .traces
            .iter()
            .any(|r| (r.y0 + r.l - resolved.box_len).abs() < 1e-9);
        assert!(touches_y0, "an input feed must reach the y=0 port plane");
        assert!(
            touches_ylen,
            "an output feed must reach the y=box_len port plane"
        );

        // The resonators do NOT touch the port planes (they are floating,
        // gap-coupled): no resonator rectangle starts at y=0 or ends at box_len.
        let res_a = &resolved.traces[0];
        let res_b = &resolved.traces[1];
        assert!(
            res_a.y0 > 1e-9 && (res_a.y0 + res_a.l) < resolved.box_len - 1e-9,
            "resonator A must be floating (not touching either port plane)"
        );
        assert!(
            res_b.y0 > 1e-9 && (res_b.y0 + res_b.l) < resolved.box_len - 1e-9,
            "resonator B must be floating (not touching either port plane)"
        );

        // The two resonators are gap-separated in x by S.
        let gap_x = res_b.x0 - (res_a.x0 + res_a.w);
        assert!(
            (gap_x - geom.gap_s).abs() < 1e-9,
            "resonator x-gap ({:.3}mm) must equal S ({:.3}mm)",
            gap_x * 1e3,
            geom.gap_s * 1e3
        );
    }

    /// The three analytic k references are finite + positive, and the two
    /// directly comparable ones (`k_imp` coupled-line vs `k_eps` resonant-split)
    /// agree to ~order (ratio in `[1.0, 1.3]`) at this weak gap (S/W = 2) — the
    /// premise that makes the gate's `k_fem`-vs-`coupling_coefficient`
    /// comparison meaningful. If a future edit tightens S this fires and flags
    /// that the two references have diverged. FAST (no FEM solve).
    #[test]
    fn analytic_k_references_finite_positive_and_agree() {
        let geom = CoupledResonatorGeom::probe();
        let cm = coupled_microstrip(geom.trace_w, geom.gap_s, geom.sub_h, geom.eps_r);
        let k_imp = coupling_coefficient(&cm);
        let l = geom.resonator_length_m();
        let f_e = C0 / (2.0 * l * cm.eps_eff_e.sqrt());
        let f_o = C0 / (2.0 * l * cm.eps_eff_o.sqrt());
        let (lo, hi) = (f_e.min(f_o), f_e.max(f_o));
        let k_eps = (hi * hi - lo * lo) / (hi * hi + lo * lo);

        assert!(
            k_imp.is_finite() && k_imp > 0.0,
            "k_imp must be finite + positive, got {k_imp}"
        );
        assert!(
            k_eps.is_finite() && k_eps > 0.0,
            "k_eps must be finite + positive, got {k_eps}"
        );
        assert!(
            l.is_finite() && l > 0.0,
            "resonator length must be finite + positive, got {l}"
        );

        let ratio = k_imp / k_eps;
        assert!(
            (1.0..=1.3).contains(&ratio),
            "at S/W=2 the coupled-line k_imp ({k_imp:.4}) and resonant-split k_eps ({k_eps:.4}) \
             should agree to ~15% (ratio {ratio:.2}); a larger ratio means the gap is too tight \
             for the two definitions to be compared like-for-like"
        );
    }

    /// `probe()` is exactly `probe_with_gap(2.0e-3)` (the refactor that extracted
    /// the gap-parameterized constructor must not have changed the default), and
    /// `probe_with_gap` widens `box_w` by exactly the gap change (the PEC walls
    /// stay `CLEARANCE_X` clear of the two-strip pattern at any gap) while leaving
    /// the gap-independent `box_h` snap untouched. This guards the K2 gate's
    /// per-gap geometry against a box-derivation drift. FAST (no FEM solve).
    #[test]
    fn probe_with_gap_matches_default_and_scales_box_w() {
        assert_eq!(
            CoupledResonatorGeom::probe(),
            CoupledResonatorGeom::probe_with_gap(2.0e-3),
            "probe() must equal probe_with_gap(2 mm) — the K2 refactor changed the default"
        );

        let g15 = CoupledResonatorGeom::probe_with_gap(1.5e-3);
        let g30 = CoupledResonatorGeom::probe_with_gap(3.0e-3);
        // box_w widens by exactly the gap delta (clearance + 2·W are fixed).
        assert!(
            ((g30.box_w - g15.box_w) - (3.0e-3 - 1.5e-3)).abs() < 1e-12,
            "box_w must widen by exactly the gap change: dW={:.6}mm vs dS={:.6}mm",
            (g30.box_w - g15.box_w) * 1e3,
            (3.0e-3 - 1.5e-3) * 1e3
        );
        // box_h is gap-independent (the z-snap does not see the gap).
        assert!(
            (g15.box_h - g30.box_h).abs() < 1e-12,
            "box_h must be gap-independent (it only snaps sub_h onto a z-plane)"
        );
        // Each gap still clears CLEARANCE_X on both sides of the two-strip span.
        for g in [&g15, &g30] {
            let two_strip_span = 2.0 * g.trace_w + g.gap_s;
            assert!(
                (g.box_w - (two_strip_span + 2.0 * CLEARANCE_X)).abs() < 1e-12,
                "box_w must be the two-strip span + CLEARANCE_X each side"
            );
        }
    }

    /// The bisection core finds the root of a SYNTHETIC monotone-decreasing
    /// closure WITHOUT any FEM solve — this is the FEM-free unit test of the
    /// [`correct_gap_fem_k`] root-find logic (the heavy gate validates it on the
    /// real `coupled_resonator_k` curve). The closure `k(g) = 0.09 − 20.0·g` is a
    /// clean decreasing line over `g ∈ [1e-3, 4e-3]` (k: 0.07 → 0.01); for the
    /// target k = 0.040 the analytic root is `g* = (0.09 − 0.040)/20 = 2.5 mm`.
    /// Asserts it converges to that root within `tol_frac` and in the
    /// bisection-bound `⌈log2(range/(2·tol_in_g))⌉`-ish evals. FAST (no FEM).
    #[test]
    fn bisect_monotone_dec_finds_synthetic_root() {
        let k_of = |g: f64| Some(0.09 - 20.0 * g);
        let target = 0.040;
        let (lo, hi) = (1.0e-3, 4.0e-3);
        let tol_frac = 0.02;
        let max_evals = 20;

        let (g_star, k_star, n_evals, converged) =
            bisect_monotone_dec(k_of, target, lo, hi, tol_frac, max_evals);

        assert!(converged, "synthetic bisection must converge (clean line)");
        let err = (k_star - target).abs() / target;
        assert!(
            err <= tol_frac,
            "converged k {k_star:.5} must be within {tol_frac} of target {target} (err {err:.4})"
        );
        // Analytic root g* = (0.09 − target)/20 = 2.5 mm; the converged k is
        // within tol of target, so the gap is within tol/20 of 2.5 mm.
        let g_analytic = (0.09 - target) / 20.0;
        assert!(
            (g_star - g_analytic).abs() < 1.0e-4,
            "converged gap {:.4}mm must be near the analytic root {:.4}mm",
            g_star * 1e3,
            g_analytic * 1e3,
        );
        // Bisection halves the bracket each usable step. To reach |k−target| ≤
        // tol_frac·target the gap-error must fall below tol_frac·target/20; the
        // step bound is ⌈log2((hi−lo)/(tol_frac·target/20))⌉.
        let tol_g = tol_frac * target / 20.0;
        let bound = (((hi - lo) / tol_g).log2().ceil() as usize) + 1;
        assert!(
            n_evals <= bound,
            "clean bisection took {n_evals} evals, expected ≤ {bound}"
        );
    }

    /// The bisection core tolerates a single UNUSABLE (`None`) eval without
    /// crashing or stalling: the same synthetic line, but the first probe at the
    /// bracket midpoint (2.5 mm) is poisoned to `None` once. The search must nudge
    /// past it and still converge to the analytic root. This guards the
    /// outlier-robustness contract (the F1.2.1.0 run saw one extraction outlier).
    /// FAST (no FEM).
    #[test]
    fn bisect_monotone_dec_tolerates_one_outlier() {
        let mut poisoned = false;
        let k_of = |g: f64| {
            // Poison exactly the first midpoint sample (g = 2.5 mm) once.
            if !poisoned && (g - 2.5e-3).abs() < 1e-9 {
                poisoned = true;
                return None;
            }
            Some(0.09 - 20.0 * g)
        };
        let target = 0.040;
        let (g_star, k_star, n_evals, converged) =
            bisect_monotone_dec(k_of, target, 1.0e-3, 4.0e-3, 0.02, 20);

        assert!(
            converged,
            "bisection must survive one outlier and still converge (got n_evals={n_evals})"
        );
        assert!(
            (k_star - target).abs() / target <= 0.02,
            "post-outlier converged k {k_star:.5} must be within tol of {target}"
        );
        assert!(
            (g_star - (0.09 - target) / 20.0).abs() < 1.0e-4,
            "post-outlier converged gap {:.4}mm must still be near the analytic root",
            g_star * 1e3
        );
    }

    /// An ALL-unusable curve terminates honestly: `correct_gap_fem_k`'s core
    /// returns `converged = false` with `n_evals = max_evals` (it does not loop
    /// forever) and a NaN best-k. FAST (no FEM).
    #[test]
    fn bisect_monotone_dec_all_unusable_terminates() {
        let k_of = |_g: f64| None;
        let max_evals = 6;
        let (_g, k_star, n_evals, converged) =
            bisect_monotone_dec(k_of, 0.040, 1.0e-3, 4.0e-3, 0.02, max_evals);
        assert!(
            !converged,
            "an all-unusable curve must not report converged"
        );
        assert_eq!(n_evals, max_evals, "must spend the full budget, not loop");
        assert!(k_star.is_nan(), "no usable eval ⇒ best-k is NaN");
    }

    /// Sanity: the resolved probe geometry lands near the probe's measured
    /// ~63 k tets (a refactor that blew up the mesh by an order of magnitude
    /// would be caught here). FAST.
    #[test]
    fn resolved_probe_mesh_is_probe_sized() {
        let resolved = resolve_geometry(&CoupledResonatorGeom::probe());
        let tets = resolved.total_tets();
        assert!(
            (40_000..=90_000).contains(&tets),
            "resolved probe mesh ({tets} tets) should be ~63 k (probe-sized); a large change \
             means the geometry resolution drifted"
        );
    }
}
