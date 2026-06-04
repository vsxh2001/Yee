//! `yee filter synth` handler (Filter Phase F0 + F1.2.0).
//!
//! Parses a [`yee_filter::FilterSpec`] TOML, synthesizes the lowpass prototype
//! and all-pole coupling matrix, sweeps the closed-form ideal response, writes
//! the S-parameters as a Touchstone `.s2p` via `yee-io`, prints the spec-mask
//! verdict, and (F1.2.0) emits the physical edge-coupled microstrip dimensions
//! — optionally writing the layout SVG. Exit 0 on PASS, 1 on FAIL (or on an
//! unrealizable coupling for the chosen substrate).

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use num_complex::Complex64;
use std::process::ExitCode;

use yee_filter::{
    FilterSpec, Footprint, check_mask, dimension_edge_coupled, dimension_edge_coupled_layout,
    ideal_response, ladder_s_params_lossy, lumped_board, synthesize, synthesize_lumped,
};
use yee_io::touchstone::{File, Format, FreqUnit};
use yee_layout::Substrate;
use yee_plotters::{MaskKind, MaskRegion, PlotConfig, PlotFormat, draw_sparam_with_mask};

/// Number of points in the response sweep written to the Touchstone file.
const SWEEP_POINTS: usize = 401;
/// Sweep span as a multiple of the fractional bandwidth on each side of `f0`.
/// `f0·(1 ± SPAN·FBW/2)` keeps a wide skirt around the passband.
const SPAN_MULT: f64 = 6.0;

/// Run `yee filter synth <spec> [--output <out.s2p>] [--plot <out.png>]
/// [--eps-r <εr>] [--h-mm <h>] [--layout-svg <out.svg>] [--gerber <out.gbr>]
/// [--kicad-pcb <out.kicad_pcb>]`.
///
/// `eps_r` / `h_mm` describe the substrate used for the F1.2.0 physical
/// dimensioning (FR-4 defaults `4.4` / `1.6 mm` are supplied by the CLI). When
/// the synthesized couplings cannot be realized on that substrate the dims path
/// prints a diagnostic and returns a non-zero [`ExitCode`] — it is never
/// silently skipped.
///
/// `--layout-svg` / `--gerber` / `--kicad-pcb` all emit the **same** single
/// [`Layout`]; it is built **once** when any of them is set, so the SVG, Gerber,
/// and KiCad board can never diverge. With `lumped == false` (default) that
/// layout is the edge-coupled distributed layout (`dimension_edge_coupled_layout`);
/// with `lumped == true` (F2.2, ADR-0158) it is the lumped-LC board
/// (`synthesize_lumped` → `lumped_board`) on the SAME `--eps-r`/`--h-mm`
/// substrate, using the supplied SMD `footprint`. The synthesis printout and the
/// Touchstone `.s2p` are unaffected by `--lumped`.
///
/// `q_unloaded` (F2-Q, ADR-0161) selects the **`(S11, S21)` response** written
/// to the `.s2p` (and the optional `--plot`'s `|S21|`): `None` (default) keeps
/// the lossless closed-form `ideal_response`, with `S11` the true lossless
/// reflection in quadrature (`j·√(1−|S21|²)`); `Some(q)` with finite `q > 0`
/// sweeps the realistic **finite-Q lumped-LC TRUE lossy 2-port**
/// (`synthesize_lumped` → `ladder_s_params_lossy(&ladder, f, q)`), so the
/// exported Touchstone carries the midband insertion loss / rounded corners a
/// built filter measures **and the true absorptive reflection** (`|S11|²+|S21|²
/// < 1`), not a lossless `√(1−|S21|²)` placeholder. It is independent of
/// `--lumped` (the lumped ladder is always synthesizable for a band-pass spec)
/// and does not affect the mask verdict, which always grades the ideal `proj`
/// via `check_mask`. A non-finite or `q <= 0` value is rejected with an error +
/// [`ExitCode::FAILURE`] (never silently treated as ideal).
// The CLI options are a flat thread-through of independent `synth` flags; a
// builder/struct here would only obscure the one-to-one mapping with the clap
// `Synth` variant.
#[allow(clippy::too_many_arguments)]
pub fn run_synth(
    spec_path: &Path,
    output: Option<&Path>,
    plot: Option<&Path>,
    eps_r: f64,
    h_mm: f64,
    layout_svg: Option<&Path>,
    gerber: Option<&Path>,
    kicad_pcb: Option<&Path>,
    lumped: bool,
    footprint: Footprint,
    q_unloaded: Option<f64>,
) -> Result<ExitCode> {
    let text = std::fs::read_to_string(spec_path)
        .with_context(|| format!("failed to read filter spec {}", spec_path.display()))?;
    let spec: FilterSpec = toml::from_str(&text)
        .with_context(|| format!("failed to parse filter spec {}", spec_path.display()))?;

    let proj = synthesize(&spec);

    // ---- prototype g-values ----------------------------------------------
    let g = &proj.prototype.g;
    let n = proj.prototype.order();
    println!("Filter synthesis ({:?}, order N={n})", spec.approximation);
    println!(
        "  f0 = {:.6e} Hz   FBW = {:.4}   Z0 = {} Ohm",
        spec.f0_hz, spec.fbw, spec.z0_ohm
    );
    print!("  prototype g-values: g0={:.4}", g[0]);
    for (i, gi) in g.iter().enumerate().skip(1) {
        print!("  g{i}={gi:.4}");
    }
    println!();

    // ---- coupling matrix + external Q ------------------------------------
    println!(
        "  external Q: Qe_in={:.4}  Qe_out={:.4}",
        proj.coupling.qe_in, proj.coupling.qe_out
    );
    println!("  coupling matrix M (normalized, {n}x{n}):");
    for row in &proj.coupling.m {
        let cells: Vec<String> = row.iter().map(|v| format!("{v:+.4}")).collect();
        println!("    [ {} ]", cells.join("  "));
    }

    // ---- response sweep (ideal lossless OR finite-Q lumped) --------------
    // `--q-unloaded` (ADR-0161) routes the `(S11, S21)` response written to the
    // .s2p (and the optional --plot's |S21|): unset → the lossless closed-form
    // `ideal_response`, with S11 the true lossless reflection placed in
    // quadrature (`S11 = j·√(1−|S21|²)`, see `write_s2p`); set → the realistic
    // finite-Q lumped-LC realization, with S11 the **true absorptive**
    // reflection from the same ABCD as S21 (`ladder_s_params_lossy`), so the
    // exported Touchstone is a true lossy 2-port (`|S11|²+|S21|² < 1`), not a
    // lossless placeholder. The mask verdict below is independent — it always
    // grades the ideal `proj` via `check_mask`.
    let freqs = sweep_freqs(spec.f0_hz, spec.fbw);
    let s_params: Vec<(Complex64, Complex64)> = match q_unloaded {
        None => ideal_response(&proj, &freqs)
            .into_iter()
            .map(lossless_s_pair)
            .collect(),
        Some(q) => {
            // Reject a non-finite / non-positive Q: the user asked for a
            // finite-Q response, so silently falling back to lossless would be
            // wrong. (ladder_s_params_lossy itself treats q<=0 as lossless.)
            if !q.is_finite() || q <= 0.0 {
                anyhow::bail!(
                    "--q-unloaded must be a finite quality factor > 0 (got {q}); \
                     omit the flag for the ideal lossless response"
                );
            }
            // The finite-Q response uses the lumped-LC realization (the same
            // ladder `--lumped` exports as a board); always synthesizable for a
            // band-pass spec, independent of `--lumped`. S11 and S21 both come
            // from one ABCD, so S11 is the TRUE absorptive reflection.
            let ladder = synthesize_lumped(&proj)
                .map_err(|e| anyhow::anyhow!("failed to synthesize lumped LC ladder: {e}"))?;
            let s_params: Vec<(Complex64, Complex64)> = freqs
                .iter()
                .map(|&f| ladder_s_params_lossy(&ladder, f, q))
                .collect();
            // Realized midband insertion loss at the sweep point nearest f0.
            let idx = nearest_freq_index(&freqs, spec.f0_hz);
            let (s11_mid, s21_mid) = s_params[idx];
            let il_db = -20.0 * s21_mid.norm().max(1e-12).log10();
            let absorbed = 1.0 - (s11_mid.norm_sqr() + s21_mid.norm_sqr());
            println!(
                "  finite-Q response (Q_unloaded = {q:.4}): midband insertion loss = {il_db:.4} dB"
            );
            println!(
                "  .s2p carries the TRUE lossy 2-port (midband |S11|²+|S21|² = {:.4}, \
                 absorbed = {absorbed:.4}); not a lossless placeholder",
                s11_mid.norm_sqr() + s21_mid.norm_sqr()
            );
            s_params
        }
    };
    // Transmission magnitudes for the optional --plot.
    let s21: Vec<Complex64> = s_params.iter().map(|&(_, s21)| s21).collect();

    // ---- spec-mask verdict -----------------------------------------------
    // Grade over the same sweep so the in-band ripple/RL is well-sampled.
    let report = check_mask(&proj, &freqs);
    println!(
        "  mask: passband ripple {:.3} dB (spec {:.3}), in-band RL {:.3} dB (spec {:.3})",
        report.worst_passband_ripple_db,
        spec.mask.passband_ripple_db,
        report.worst_return_loss_db,
        spec.mask.return_loss_db,
    );
    for (f_hz, achieved, required, met) in &report.stopband {
        println!(
            "  stopband {f_hz:.4e} Hz: rejection {achieved:.2} dB (required {required:.2} dB) {}",
            if *met { "OK" } else { "UNDER" }
        );
    }
    for fail in &report.failures {
        println!("  FAILURE: {fail}");
    }

    // ---- write Touchstone .s2p -------------------------------------------
    let out_path = output
        .map(PathBuf::from)
        .unwrap_or_else(|| spec_path.with_extension("s2p"));
    // The comment self-describes which response the .s2p carries. The default
    // (ideal) string is kept byte-for-byte so an unset `--q-unloaded` run is
    // byte-identical to the pre-ADR-0161 behavior.
    let s2p_comment = match q_unloaded {
        None => " yee filter synth — ideal closed-form response".to_string(),
        Some(q) => format!(" yee filter synth — finite-Q lumped response (Q_unloaded = {q})"),
    };
    write_s2p(&out_path, spec.z0_ohm, &freqs, &s_params, &s2p_comment)
        .with_context(|| format!("failed to write Touchstone {}", out_path.display()))?;
    println!("  wrote Touchstone: {}", out_path.display());

    // ---- optional spec-mask plot -----------------------------------------
    if let Some(plot_path) = plot {
        let s21_db: Vec<f64> = s21
            .iter()
            .map(|z| 20.0 * z.norm().max(1e-12).log10())
            .collect();
        let regions = spec_mask_regions(&spec);
        let format = if plot_path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("svg"))
        {
            PlotFormat::Svg
        } else {
            PlotFormat::Png
        };
        let cfg = PlotConfig {
            title: format!("{:?} filter |S21| vs spec mask", spec.approximation),
            format,
            ..PlotConfig::default()
        };
        draw_sparam_with_mask(plot_path, &freqs, &[("S21", &s21_db)], &regions, &cfg)
            .map_err(|e| anyhow::anyhow!("failed to write plot {}: {e}", plot_path.display()))?;
        println!("  wrote plot: {}", plot_path.display());
    }

    // ---- F1.2.0 physical dimensions --------------------------------------
    // Build the substrate from the CLI εr/h (h supplied in mm). The remaining
    // Substrate fields (loss tangent, metal thickness) do not enter the
    // closed-form dimensioning; use neutral zero defaults.
    let substrate = Substrate {
        eps_r,
        height_m: h_mm * 1e-3,
        loss_tangent: 0.0,
        metal_thickness_m: 0.0,
    };
    println!("  substrate: eps_r = {eps_r:.4}   h = {h_mm:.4} mm");

    // Planar edge-coupled dimensioning + printout — ONLY when not `--lumped`.
    // The lumped-LC board (ADR-0158) does not depend on edge-coupled
    // realizability on this substrate, so under `--lumped` we must neither fail
    // here (a spec that is lumped-realizable but distributed-unrealizable would
    // otherwise spuriously exit 1) nor print the irrelevant planar dimensions.
    if !lumped {
        let dims = match dimension_edge_coupled(&proj, &substrate) {
            Ok(dims) => dims,
            Err(e) => {
                // Surface the unrealizable-coupling (or topology/order) error and
                // exit non-zero — never silently skip the dimensions.
                eprintln!("  ERROR: cannot dimension edge-coupled filter: {e}");
                println!("VERDICT: FAIL (dimensioning)");
                return Ok(ExitCode::FAILURE);
            }
        };

        println!("  physical dimensions (edge-coupled half-wave microstrip):");
        println!(
            "    line width       = {:.6e} m  ({:.4} mm)",
            dims.line_width_m,
            dims.line_width_m * 1e3
        );
        println!(
            "    resonator length = {:.6e} m  ({:.4} mm)",
            dims.resonator_length_m,
            dims.resonator_length_m * 1e3
        );
        for (i, (gap, k)) in dims.gaps_m.iter().zip(dims.target_k.iter()).enumerate() {
            println!(
                "    gap[{i}] = {:.6e} m  ({:.4} mm)   target_k = {k:.6}",
                gap,
                gap * 1e3
            );
        }
    }

    // ---- optional layout exports (SVG / Gerber / KiCad) ------------------
    // Build the export layout ONCE when any exporter is requested, so the
    // `--layout-svg`, `--gerber`, and `--kicad-pcb` outputs can never diverge.
    // `--lumped` (F2.2, ADR-0158) selects the lumped-LC board (signal line +
    // ground rail + every L/C component pad) over the planar edge-coupled
    // layout; both reuse the SAME `substrate` built above from `--eps-r`/`--h-mm`.
    if layout_svg.is_some() || gerber.is_some() || kicad_pcb.is_some() {
        let layout = if lumped {
            let ladder = synthesize_lumped(&proj)
                .map_err(|e| anyhow::anyhow!("failed to synthesize lumped LC ladder: {e}"))?;
            let board = lumped_board(&ladder, &substrate, footprint);
            let n_placements = board.placements.len();
            println!(
                "  lumped board: {} components ({} resonators), footprint {:?}",
                n_placements,
                ladder.resonators.len(),
                footprint
            );
            board.layout
        } else {
            dimension_edge_coupled_layout(&proj, &substrate)
                .map_err(|e| anyhow::anyhow!("failed to build layout: {e}"))?
        };
        if let Some(svg_path) = layout_svg {
            std::fs::write(svg_path, layout.to_svg())
                .with_context(|| format!("failed to write layout SVG {}", svg_path.display()))?;
            println!("  wrote layout SVG: {}", svg_path.display());
        }
        if let Some(gerber_path) = gerber {
            let gerber_text =
                yee_export::layout_to_gerber(&layout, &yee_export::GerberOptions::default());
            std::fs::write(gerber_path, gerber_text).with_context(|| {
                format!("failed to write layout Gerber {}", gerber_path.display())
            })?;
            println!("  wrote layout Gerber: {}", gerber_path.display());
        }
        if let Some(kicad_path) = kicad_pcb {
            let text =
                yee_export::layout_to_kicad_pcb(&layout, &yee_export::KicadPcbOptions::default());
            std::fs::write(kicad_path, text)
                .with_context(|| format!("failed to write KiCad PCB {}", kicad_path.display()))?;
            println!("  wrote KiCad PCB: {}", kicad_path.display());
        }
    }

    if report.pass {
        println!("VERDICT: PASS");
        Ok(ExitCode::SUCCESS)
    } else {
        println!("VERDICT: FAIL");
        Ok(ExitCode::FAILURE)
    }
}

/// Index of the swept frequency nearest `f0` (used to report the midband
/// insertion loss / absorption). Returns 0 for an empty / degenerate sweep.
fn nearest_freq_index(freqs: &[f64], f0: f64) -> usize {
    freqs
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            (*a - f0)
                .abs()
                .partial_cmp(&(*b - f0).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
        .unwrap_or(0)
}

/// Map a lossless transfer `S21` to the `(S11, S21)` pair written for the ideal
/// (no `--q-unloaded`) path. The transmission is taken real with `|S21| ≤ 1`,
/// and S11 is the **true lossless reflection** placed in quadrature
/// (`S11 = j·√(1−|S21|²)`): for a lossless reciprocal symmetric 2-port this is
/// the exact reflection (it makes `|S11|²+|S21|² = 1` *and* the S-matrix
/// unitary/passive). This reproduces the prior ideal-path bytes exactly.
fn lossless_s_pair(s21: Complex64) -> (Complex64, Complex64) {
    let s21_mag = s21.norm().min(1.0);
    let s11_mag = (1.0 - s21_mag * s21_mag).max(0.0).sqrt();
    (Complex64::new(0.0, s11_mag), Complex64::new(s21_mag, 0.0))
}

/// Linear sweep of `SWEEP_POINTS` frequencies centred on `f0`, spanning
/// `f0·(1 ± SPAN_MULT·fbw/2)` (clamped to be strictly positive).
fn sweep_freqs(f0: f64, fbw: f64) -> Vec<f64> {
    let half = SPAN_MULT * fbw / 2.0;
    let lo = (f0 * (1.0 - half)).max(f0 * 1e-3);
    let hi = f0 * (1.0 + half);
    (0..SWEEP_POINTS)
        .map(|i| lo + (hi - lo) * (i as f64) / ((SWEEP_POINTS - 1) as f64))
        .collect()
}

/// Build and write a 2-port Touchstone file for a reciprocal, symmetric filter
/// from the per-frequency `(S11, S21)` pairs the caller supplies. The matrix is
/// `[[S11, S21], [S21, S11]]` (reciprocal `S12 = S21`, symmetric `S22 = S11`);
/// the **caller** owns the physics, so `write_s2p` writes the given complex
/// values verbatim and never re-derives S11.
///
/// - Ideal path: S11 is the true lossless reflection in quadrature
///   (`S11 = j·√(1−|S21|²)`, [`lossless_s_pair`]) → `|S11|²+|S21|² = 1`, unitary.
/// - Finite-Q path: S11 and S21 are the **true** lossy 2-port from one ABCD
///   ([`ladder_s_params_lossy`]) → `|S11|²+|S21|² < 1` (absorption), still
///   passive (`σ_max = √(|S11|²+|S21|²) < 1`), so `yee_io::touchstone` accepts
///   it on read-back.
///
/// (A real S11 with `|S11| = √(1−|S21|²)`, the original form, gave the
/// eigenvalue `|S21|+|S11| > 1` — a passivity violation `touchstone::read`
/// rejects, and it also mis-attributed absorptive loss to reflection.)
fn write_s2p(
    path: &Path,
    z0: f64,
    freqs: &[f64],
    s_params: &[(Complex64, Complex64)],
    comment: &str,
) -> yee_io::Result<()> {
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
    yee_io::touchstone::write(path, &file)
}

/// Map a [`FilterSpec`]'s spec mask to |S21| forbidden regions for the plot:
/// a passband `Floor` at `−passband_ripple_db` over `[f0·(1−fbw/2), f0·(1+fbw/2)]`,
/// and a `Ceiling` at `−reject` over a ±2 % band around each stopband point.
fn spec_mask_regions(spec: &FilterSpec) -> Vec<MaskRegion> {
    let f1 = spec.f0_hz * (1.0 - spec.fbw / 2.0);
    let f2 = spec.f0_hz * (1.0 + spec.fbw / 2.0);
    let mut regions = vec![MaskRegion {
        f_lo_hz: f1,
        f_hi_hz: f2,
        kind: MaskKind::Floor,
        limit_db: -spec.mask.passband_ripple_db,
    }];
    for &(f_s, reject_db) in &spec.mask.stopband {
        regions.push(MaskRegion {
            f_lo_hz: f_s * 0.98,
            f_hi_hz: f_s * 1.02,
            kind: MaskKind::Ceiling,
            limit_db: -reject_db,
        });
    }
    regions
}

#[cfg(test)]
mod tests {
    use super::*;
    use yee_filter::{Approximation, Response, SpecMask};

    #[test]
    fn spec_mask_regions_passband_floor_and_stopband_ceiling() {
        let spec = FilterSpec {
            response: Response::Bandpass,
            approximation: Approximation::Chebyshev { ripple_db: 0.5 },
            f0_hz: 2.0e9,
            fbw: 0.10,
            order: Some(5),
            z0_ohm: 50.0,
            mask: SpecMask {
                passband_ripple_db: 0.5,
                return_loss_db: 10.0,
                stopband: vec![(2.4e9, 40.0)],
            },
        };
        let r = spec_mask_regions(&spec);
        assert_eq!(r.len(), 2, "one passband Floor + one stopband Ceiling");

        assert_eq!(r[0].kind, MaskKind::Floor);
        assert!((r[0].f_lo_hz - 1.9e9).abs() < 1.0, "passband lo edge");
        assert!((r[0].f_hi_hz - 2.1e9).abs() < 1.0, "passband hi edge");
        assert!((r[0].limit_db - (-0.5)).abs() < 1e-9, "floor at -ripple");

        assert_eq!(r[1].kind, MaskKind::Ceiling);
        assert!((r[1].f_lo_hz - 2.352e9).abs() < 1.0, "stopband lo (-2%)");
        assert!((r[1].f_hi_hz - 2.448e9).abs() < 1.0, "stopband hi (+2%)");
        assert!((r[1].limit_db - (-40.0)).abs() < 1e-9, "ceiling at -reject");
    }
}
