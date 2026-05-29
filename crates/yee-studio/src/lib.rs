//! `yee-studio` — App.0 filter-design studio state (logic layer).
//!
//! This crate holds the **egui-free, headless-testable** application state for
//! the Yee Filter Studio desktop app (ADR-0090). It wires the shipped light
//! flow — spec → synthesis → ideal response → spec-mask verdict — into a single
//! [`StudioState`] value that the `eframe` shell (`src/app.rs`) edits and
//! re-derives on every change.
//!
//! Keeping [`StudioState`] free of any `egui`/`eframe` types is deliberate: it
//! keeps the logic WASM-safe (App.1) and unit-testable without a GUI runtime.
//! Only `src/app.rs` and `src/main.rs` depend on `egui`/`eframe`.
//!
//! ## Pipeline ([`StudioState::recompute`])
//!
//! 1. [`yee_filter::synthesize`] the spec → [`yee_filter::FilterProject`].
//! 2. Build the 401-point sweep `f0·(1 ± 6·fbw/2)` (mirrors `yee-cli`).
//! 3. [`yee_filter::ideal_response`] over the sweep → complex `S21`.
//! 4. `s21_db = 20·log10(max(|S21|, 1e-12))`.
//! 5. Closed-form spec-mask regions ([`MaskRegionView`]).
//! 6. [`yee_filter::check_mask`] → `mask_pass` + human-readable notes.
//! 7. [`yee_filter::dimension_edge_coupled`] over the synthesized project and the
//!    editable substrate → physical microstrip [`EdgeCoupledDimensions`] (F1.2.0;
//!    ADR-0099), stored as a `Result<_, String>` so the error path stays
//!    egui-free and WASM-safe.

use yee_filter::{
    EdgeCoupledDimensions, FilterProject, FilterSpec, check_mask, dimension_edge_coupled,
    dimension_edge_coupled_layout, ideal_response, synthesize,
};
use yee_layout::Substrate;

/// The `eframe` shell (spec editor + synthesis panel + `|S21|`/mask plot).
///
/// Lives in its own module so this crate root stays egui-free and WASM-safe;
/// only [`app`] and the binary entry depend on `egui`/`eframe`.
///
/// Gated behind either the `desktop` (native) or `web` (wasm32 browser) Cargo
/// feature, both of which pull in `eframe`/`egui`/`egui_plot` (App.1.0/1.2a;
/// ADR-0092/0096): a `--no-default-features` build compiles [`StudioState`] with
/// **no** `eframe`/`egui`/`wgpu` in the dep graph, satisfying the ADR-0089
/// WASM-safety constraint. The `web` feature compiles this same UI for
/// `wasm32-unknown-unknown` behind the `WebRunner` entry in `src/main.rs`; the
/// `trunk` bundle + deploy remain for App.1.2b (they need the wasm toolchain).
#[cfg(any(feature = "desktop", feature = "web"))]
pub mod app;

/// Number of points in the response sweep (mirrors `yee-cli`'s `SWEEP_POINTS`).
const SWEEP_POINTS: usize = 401;
/// Sweep span as a multiple of the fractional bandwidth on each side of `f0`
/// (mirrors `yee-cli`'s `SPAN_MULT`): `f0·(1 ± SPAN_MULT·fbw/2)`.
const SPAN_MULT: f64 = 6.0;

/// A plain, egui-free view of one spec-mask forbidden region for plotting.
///
/// Mirrors the `yee-cli` F0.2 mapping: a passband `Floor` at `−passband_ripple`
/// over `[f0·(1−fbw/2), f0·(1+fbw/2)]`, and a per-stopband `Ceiling` at
/// `−reject` over a ±2 % band around each stopband point.
#[derive(Debug, Clone, PartialEq)]
pub struct MaskRegionView {
    /// Lower frequency edge of the region, Hz.
    pub f_lo_hz: f64,
    /// Upper frequency edge of the region, Hz.
    pub f_hi_hz: f64,
    /// `true` for a passband **floor** (|S21| must stay *above* `limit_db`);
    /// `false` for a stopband **ceiling** (|S21| must stay *below* `limit_db`).
    pub floor: bool,
    /// The mask limit on `|S21|`, in dB.
    pub limit_db: f64,
}

/// Editable design state plus everything derived from the current spec.
///
/// The `spec` field is the single source of truth; every other field is a
/// cached derivation produced by [`StudioState::recompute`]. The `eframe` shell
/// mutates `spec` on edits and calls `recompute` to refresh the rest.
#[derive(Debug, Clone)]
pub struct StudioState {
    /// The editable filter specification (the design intent).
    pub spec: FilterSpec,
    /// Synthesized project (prototype g-values, coupling matrix, topology).
    pub project: FilterProject,
    /// Swept frequencies, Hz (the 401-point `f0·(1 ± 6·fbw/2)` grid).
    pub freqs_hz: Vec<f64>,
    /// `|S21|` in dB over [`StudioState::freqs_hz`]
    /// (`20·log10(max(|S21|, 1e-12))`).
    pub s21_db: Vec<f64>,
    /// Spec-mask forbidden regions for the plot.
    pub mask_regions: Vec<MaskRegionView>,
    /// Overall mask verdict (`true` ⇒ PASS).
    pub mask_pass: bool,
    /// Human-readable mask notes (ripple/RL summary + per-stopband + failures).
    pub mask_notes: Vec<String>,
    /// Substrate relative permittivity `ε_r` used for dimensional synthesis
    /// (default `4.4`, FR-4). Editable; changing it re-derives [`dims`].
    ///
    /// [`dims`]: StudioState::dims
    pub eps_r: f64,
    /// Substrate height `h` (metal-to-ground spacing), metres, used for
    /// dimensional synthesis (default `1.6e-3`, FR-4). Editable; changing it
    /// re-derives [`dims`].
    ///
    /// [`dims`]: StudioState::dims
    pub h_m: f64,
    /// Physical microstrip dimensions (F1.2.0; ADR-0099) for the synthesized
    /// edge-coupled filter on the [`eps_r`]/[`h_m`] substrate, or the
    /// [`yee_filter::DimError`] display string when the coupling cannot be
    /// realized on the chosen substrate.
    ///
    /// The error is stored as a `String` (not the `DimError`) so [`StudioState`]
    /// holds no `egui`/native type and stays WASM-safe (ADR-0089/0099).
    ///
    /// [`eps_r`]: StudioState::eps_r
    /// [`h_m`]: StudioState::h_m
    pub dims: Result<EdgeCoupledDimensions, String>,
    /// Top-view planar layout (F1.2.0; ADR-0101) — the edge-coupled microstrip
    /// trace polygons for the synthesized filter on the [`eps_r`]/[`h_m`]
    /// substrate, or the [`yee_filter::DimError`] display string when the
    /// coupling cannot be realized on the chosen substrate.
    ///
    /// Derived alongside [`dims`] from the same `substrate`. The error is stored
    /// as a `String` (not the `DimError`) and [`yee_layout::Layout`] carries no
    /// `egui`/native type, so [`StudioState`] stays egui-free and WASM-safe
    /// (ADR-0089/0101). Coordinates in the polygons are in metres.
    ///
    /// [`dims`]: StudioState::dims
    /// [`eps_r`]: StudioState::eps_r
    /// [`h_m`]: StudioState::h_m
    pub layout: Result<yee_layout::Layout, String>,
}

impl StudioState {
    /// Build a [`StudioState`] from a [`FilterSpec`], running [`recompute`]
    /// immediately so every derived field is populated.
    ///
    /// [`recompute`]: StudioState::recompute
    pub fn from_spec(spec: FilterSpec) -> Self {
        let project = synthesize(&spec);
        let mut state = Self {
            spec,
            project,
            freqs_hz: Vec::new(),
            s21_db: Vec::new(),
            mask_regions: Vec::new(),
            mask_pass: false,
            mask_notes: Vec::new(),
            // FR-4 substrate defaults (ε_r = 4.4 on a 1.6 mm board).
            eps_r: 4.4,
            h_m: 1.6e-3,
            // Populated by `apply_derived` below.
            dims: Err(String::new()),
            layout: Err(String::new()),
        };
        // `project` is already synthesized above; derive the rest (no re-synth).
        state.apply_derived();
        state
    }

    /// Re-derive every cached field from the current [`StudioState::spec`].
    ///
    /// Synthesizes the project, sweeps the ideal response, computes the `|S21|`
    /// dB trace, builds the spec-mask regions, and grades the mask. Call this
    /// after any edit to `spec`.
    pub fn recompute(&mut self) {
        // Synthesis (the only place `synthesize` runs); the response + mask
        // fields derive from the fresh project via `apply_derived`.
        self.project = synthesize(&self.spec);
        self.apply_derived();
    }

    /// Re-derive the response + spec-mask fields from the current `spec` and
    /// `project` — everything *except* synthesis. Shared by [`recompute`] and
    /// [`from_spec`] so the project is synthesized exactly once per build/edit.
    ///
    /// [`recompute`]: StudioState::recompute
    /// [`from_spec`]: StudioState::from_spec
    fn apply_derived(&mut self) {
        // Sweep grid (mirrors yee-cli's 401-pt `f0·(1 ± 6·fbw/2)`).
        self.freqs_hz = sweep_freqs(self.spec.f0_hz, self.spec.fbw);

        // Ideal response → |S21| dB (floored at 1e-12).
        let s21 = ideal_response(&self.project, &self.freqs_hz);
        self.s21_db = s21
            .iter()
            .map(|z| 20.0 * z.norm().max(1e-12).log10())
            .collect();

        // Spec-mask regions.
        self.mask_regions = spec_mask_regions(&self.spec);

        // Mask verdict + notes.
        let report = check_mask(&self.project, &self.freqs_hz);
        self.mask_pass = report.pass;
        self.mask_notes = mask_notes(&self.spec, &report);

        // Physical microstrip dimensions (F1.2.0; ADR-0099). Only `eps_r` and
        // `height_m` feed the closed-form synthesis; `loss_tangent` /
        // `metal_thickness_m` are unused here, set to neutral defaults. The
        // `DimError` is mapped to its display string so `StudioState` holds no
        // non-`Result` error type and stays egui-free / WASM-safe (ADR-0089).
        let substrate = Substrate {
            eps_r: self.eps_r,
            height_m: self.h_m,
            loss_tangent: 0.0,
            metal_thickness_m: 0.0,
        };
        self.dims = dimension_edge_coupled(&self.project, &substrate).map_err(|e| e.to_string());

        // Top-view layout polygons (F1.2.0; ADR-0101) from the same substrate.
        // The `DimError` is mapped to its display string so `StudioState` holds
        // no non-`Result` error type and stays egui-free / WASM-safe.
        self.layout =
            dimension_edge_coupled_layout(&self.project, &substrate).map_err(|e| e.to_string());
    }
}

/// Linear sweep of [`SWEEP_POINTS`] frequencies centred on `f0`, spanning
/// `f0·(1 ± SPAN_MULT·fbw/2)` (clamped strictly positive). Mirrors `yee-cli`.
fn sweep_freqs(f0: f64, fbw: f64) -> Vec<f64> {
    let half = SPAN_MULT * fbw / 2.0;
    let lo = (f0 * (1.0 - half)).max(f0 * 1e-3);
    let hi = f0 * (1.0 + half);
    (0..SWEEP_POINTS)
        .map(|i| lo + (hi - lo) * (i as f64) / ((SWEEP_POINTS - 1) as f64))
        .collect()
}

/// Map a [`FilterSpec`]'s spec mask to `|S21|` forbidden regions (mirrors
/// `yee-cli`'s `spec_mask_regions`): one passband [`MaskRegionView`] floor at
/// `−passband_ripple_db` over `[f0·(1−fbw/2), f0·(1+fbw/2)]`, and one ceiling
/// per stopband point at `−reject` over a ±2 % band.
fn spec_mask_regions(spec: &FilterSpec) -> Vec<MaskRegionView> {
    let f1 = spec.f0_hz * (1.0 - spec.fbw / 2.0);
    let f2 = spec.f0_hz * (1.0 + spec.fbw / 2.0);
    let mut regions = vec![MaskRegionView {
        f_lo_hz: f1,
        f_hi_hz: f2,
        floor: true,
        limit_db: -spec.mask.passband_ripple_db,
    }];
    for &(f_s, reject_db) in &spec.mask.stopband {
        regions.push(MaskRegionView {
            f_lo_hz: f_s * 0.98,
            f_hi_hz: f_s * 1.02,
            floor: false,
            limit_db: -reject_db,
        });
    }
    regions
}

/// Build the human-readable mask notes shown in the synthesis panel.
fn mask_notes(spec: &FilterSpec, report: &yee_filter::MaskReport) -> Vec<String> {
    let mut notes = Vec::new();
    notes.push(format!(
        "passband ripple {:.3} dB (spec {:.3})",
        report.worst_passband_ripple_db, spec.mask.passband_ripple_db
    ));
    notes.push(format!(
        "in-band return loss {:.3} dB (spec {:.3})",
        report.worst_return_loss_db, spec.mask.return_loss_db
    ));
    for (f_hz, achieved, required, met) in &report.stopband {
        notes.push(format!(
            "stopband {:.4e} Hz: rejection {:.2} dB (need {:.2} dB) {}",
            f_hz,
            achieved,
            required,
            if *met { "OK" } else { "UNDER" }
        ));
    }
    for fail in &report.failures {
        notes.push(format!("FAILURE: {fail}"));
    }
    notes
}

#[cfg(test)]
mod tests {
    use super::*;
    use yee_filter::{Approximation, Response, SpecMask};

    /// The default satisfiable Chebyshev 0.5 dB N=5 BPF used by the bin and by
    /// the pass test (f0 = 2 GHz, FBW = 0.10, stopband (2.4 GHz, 40 dB)).
    ///
    /// Return loss is 9.0 dB, matching the committed `cheb_bpf.toml` fixture the
    /// spec says to reuse: a 0.5 dB-ripple Chebyshev caps in-band RL at
    /// ≈ 9.64 dB, so the DoD's prose "RL 10" is unsatisfiable for this shape —
    /// the canonical satisfiable value is 9.0 (see report).
    fn satisfiable_spec() -> FilterSpec {
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

    #[test]
    fn studio_state_recompute_pass() {
        let state = StudioState::from_spec(satisfiable_spec());

        // Mask verdict: PASS.
        assert!(state.mask_pass, "default spec should satisfy its mask");

        // Published F0 Chebyshev 0.5 dB N=5 g-values (g1..g5).
        let expected = [1.7058_f64, 1.2296, 2.5408, 1.2296, 1.7058];
        let g = &state.project.prototype.g;
        for (i, &want) in expected.iter().enumerate() {
            let got = g[i + 1]; // g[1..=5]
            assert!(
                (got - want).abs() < 1e-3,
                "g[{}] = {got}, expected {want}",
                i + 1
            );
        }

        // Parallel arrays.
        assert_eq!(state.freqs_hz.len(), 401, "401-point sweep");
        assert_eq!(
            state.s21_db.len(),
            state.freqs_hz.len(),
            "s21_db parallels freqs_hz"
        );

        // Mask regions: one Floor ≈ [1.9, 2.1] GHz @ −0.5 dB, one Ceiling
        // ≈ [2.352, 2.448] GHz @ −40 dB.
        assert_eq!(state.mask_regions.len(), 2, "one Floor + one Ceiling");

        let floor = &state.mask_regions[0];
        assert!(floor.floor, "region 0 is the passband floor");
        assert!((floor.f_lo_hz - 1.9e9).abs() < 1.0, "floor lo edge");
        assert!((floor.f_hi_hz - 2.1e9).abs() < 1.0, "floor hi edge");
        assert!((floor.limit_db - (-0.5)).abs() < 1e-9, "floor at −ripple");

        let ceil = &state.mask_regions[1];
        assert!(!ceil.floor, "region 1 is the stopband ceiling");
        assert!((ceil.f_lo_hz - 2.352e9).abs() < 1.0, "ceiling lo (−2 %)");
        assert!((ceil.f_hi_hz - 2.448e9).abs() < 1.0, "ceiling hi (+2 %)");
        assert!((ceil.limit_db - (-40.0)).abs() < 1e-9, "ceiling at −reject");
    }

    #[test]
    fn studio_state_recompute_fail() {
        // A too-low order (N=2) cannot meet the 40 dB stopband at 2.4 GHz.
        let mut spec = satisfiable_spec();
        spec.order = Some(2);
        let state = StudioState::from_spec(spec);
        assert!(
            !state.mask_pass,
            "order N=2 should fail the 40 dB stopband mask"
        );
    }

    #[test]
    fn studio_state_dims() {
        // Default Chebyshev N=5 spec + the FR-4 substrate defaults (ε_r = 4.4,
        // h = 1.6 mm) should size onto a realizable edge-coupled microstrip
        // filter (F1.2.0; ADR-0099).
        let state = StudioState::from_spec(satisfiable_spec());

        // The FR-4 defaults are wired in the constructor.
        assert!(
            (state.eps_r - 4.4).abs() < 1e-12,
            "default ε_r = 4.4 (FR-4)"
        );
        assert!(
            (state.h_m - 1.6e-3).abs() < 1e-12,
            "default h = 1.6 mm (FR-4)"
        );

        let dims = state
            .dims
            .as_ref()
            .expect("default N=5 spec must dimension onto FR-4");

        assert!(
            dims.line_width_m > 0.0,
            "line width must be positive, got {}",
            dims.line_width_m
        );
        assert!(
            dims.resonator_length_m > 0.0,
            "resonator length must be positive, got {}",
            dims.resonator_length_m
        );
        // N=5 ⇒ N−1 = 4 inter-resonator gaps, every one positive.
        assert_eq!(dims.gaps_m.len(), 4, "N=5 ⇒ 4 inter-resonator gaps");
        for (i, &gap) in dims.gaps_m.iter().enumerate() {
            assert!(gap > 0.0, "gap[{i}] must be positive, got {gap}");
        }
    }

    #[test]
    fn studio_state_layout() {
        // Default Chebyshev N=5 spec + the FR-4 substrate defaults (ε_r = 4.4,
        // h = 1.6 mm) should produce a realizable top-view layout (ADR-0101).
        let state = StudioState::from_spec(satisfiable_spec());

        let layout = state
            .layout
            .as_ref()
            .expect("default N=5 spec must produce a layout on FR-4");

        assert!(
            !layout.traces.is_empty(),
            "layout must have at least one trace polygon"
        );
        for (i, poly) in layout.traces.iter().enumerate() {
            assert!(
                poly.verts.len() >= 3,
                "trace polygon {i} must have ≥ 3 vertices, got {}",
                poly.verts.len()
            );
        }
    }
}
