//! Inline-SVG renderers for the response plot and the board top-view.
//!
//! Both return self-contained `<svg>…</svg>` markup strings styled with the
//! design-system palette, injected into the RSX via `dangerous_inner_html`.
//! The board view re-uses the *real* [`yee_layout::Layout`] polygons; the
//! response plot the *real* swept [`SweepPoint`]s + mask bands.

use yee_filter::LumpedBoard;
use yee_layout::Layout;

use crate::engine::{MaskBand, SweepPoint};

// Design-system colours (kept in sync with assets/studio.css).
const ACCENT: &str = "#2dd4bf";
const COPPER: &str = "#e6b24d";
const COPPER_EDGE: &str = "#b87814";
const SUBSTRATE: &str = "#0a2218";
const SUBSTRATE_EDGE: &str = "#1f5138";
const FAIL: &str = "#e35d6a";
const GRID: &str = "#1c222a";
const MUTED: &str = "#6b7480";
const TEXT: &str = "#c9d1d9";

/// Render the ideal `|S21|` / `|S11|` response over the swept band, with the
/// spec-mask forbidden regions shaded and axis labels — as an inline SVG.
///
/// X axis = frequency (GHz), Y axis = magnitude (dB), clamped to `[y_min, 0]`.
pub fn response_plot(sweep: &[SweepPoint], bands: &[MaskBand]) -> String {
    // viewBox geometry.
    let w = 720.0;
    let h = 320.0;
    let ml = 52.0; // left margin (y labels)
    let mr = 16.0;
    let mt = 14.0;
    let mb = 36.0; // bottom margin (x labels)
    let pw = w - ml - mr;
    let ph = h - mt - mb;

    let f_lo = sweep.first().map(|s| s.f_hz).unwrap_or(0.0);
    let f_hi = sweep.last().map(|s| s.f_hz).unwrap_or(1.0);
    let y_min = -80.0_f64;
    let y_max = 5.0_f64;

    let fx = |f: f64| ml + (f - f_lo) / (f_hi - f_lo) * pw;
    let fy = |db: f64| {
        let c = db.clamp(y_min, y_max);
        mt + (y_max - c) / (y_max - y_min) * ph
    };

    let mut s = String::new();
    s.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {w} {h}\" preserveAspectRatio=\"xMidYMid meet\" role=\"img\" aria-label=\"ideal S-parameter response versus spec mask\">\n"
    ));

    // ---- gridlines + y axis labels (every 20 dB) --------------------------
    let mut db = 0.0;
    while db >= y_min {
        let y = fy(db);
        s.push_str(&format!(
            "<line x1=\"{ml:.1}\" y1=\"{y:.1}\" x2=\"{:.1}\" y2=\"{y:.1}\" stroke=\"{GRID}\" stroke-width=\"1\"/>\n",
            ml + pw
        ));
        s.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" fill=\"{MUTED}\" font-size=\"10\" font-family=\"monospace\" text-anchor=\"end\">{db:.0}</text>\n",
            ml - 6.0,
            y + 3.5
        ));
        db -= 20.0;
    }

    // ---- x axis labels (GHz) ----------------------------------------------
    let n_ticks = 6;
    for i in 0..=n_ticks {
        let f = f_lo + (f_hi - f_lo) * (i as f64) / (n_ticks as f64);
        let x = fx(f);
        s.push_str(&format!(
            "<line x1=\"{x:.1}\" y1=\"{mt:.1}\" x2=\"{x:.1}\" y2=\"{:.1}\" stroke=\"{GRID}\" stroke-width=\"1\"/>\n",
            mt + ph
        ));
        s.push_str(&format!(
            "<text x=\"{x:.1}\" y=\"{:.1}\" fill=\"{MUTED}\" font-size=\"10\" font-family=\"monospace\" text-anchor=\"middle\">{:.2}</text>\n",
            mt + ph + 16.0,
            f / 1e9
        ));
    }
    // axis titles
    s.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" fill=\"{MUTED}\" font-size=\"10\" text-anchor=\"middle\">frequency (GHz)</text>\n",
        ml + pw / 2.0,
        h - 4.0
    ));
    s.push_str(&format!(
        "<text x=\"12\" y=\"{:.1}\" fill=\"{MUTED}\" font-size=\"10\" text-anchor=\"middle\" transform=\"rotate(-90 12 {:.1})\">magnitude (dB)</text>\n",
        mt + ph / 2.0,
        mt + ph / 2.0
    ));

    // ---- mask forbidden regions (shaded) ----------------------------------
    for b in bands {
        let x0 = fx(b.f_lo_hz.max(f_lo));
        let x1 = fx(b.f_hi_hz.min(f_hi));
        if x1 <= x0 {
            continue;
        }
        let (y0, yh, fill) = if b.is_floor {
            // Passband floor: forbidden BELOW limit (toward y_min).
            let yt = fy(b.limit_db);
            (yt, (mt + ph) - yt, FAIL)
        } else {
            // Stopband ceiling: forbidden ABOVE limit (toward 0 dB / top).
            let yb = fy(b.limit_db);
            (mt, yb - mt, FAIL)
        };
        if yh > 0.0 {
            s.push_str(&format!(
                "<rect x=\"{x0:.1}\" y=\"{y0:.1}\" width=\"{:.1}\" height=\"{yh:.1}\" fill=\"{fill}\" fill-opacity=\"0.12\" stroke=\"{fill}\" stroke-opacity=\"0.35\" stroke-dasharray=\"3 3\"/>\n",
                x1 - x0
            ));
        }
    }

    // ---- traces -----------------------------------------------------------
    let poly = |pick: &dyn Fn(&SweepPoint) -> f64| -> String {
        sweep
            .iter()
            .map(|p| format!("{:.2},{:.2}", fx(p.f_hz), fy(pick(p))))
            .collect::<Vec<_>>()
            .join(" ")
    };
    // |S11| first (muted), then |S21| (accent) on top.
    s.push_str(&format!(
        "<polyline points=\"{}\" fill=\"none\" stroke=\"{MUTED}\" stroke-width=\"1.5\" stroke-opacity=\"0.8\"/>\n",
        poly(&|p| p.s11_db)
    ));
    s.push_str(&format!(
        "<polyline points=\"{}\" fill=\"none\" stroke=\"{ACCENT}\" stroke-width=\"2\"/>\n",
        poly(&|p| p.s21_db)
    ));

    // plot frame
    s.push_str(&format!(
        "<rect x=\"{ml:.1}\" y=\"{mt:.1}\" width=\"{pw:.1}\" height=\"{ph:.1}\" fill=\"none\" stroke=\"{GRID}\" stroke-width=\"1\"/>\n"
    ));

    s.push_str("</svg>\n");
    s
}

/// Render the board top-view from the real [`Layout`]: copper trace polygons
/// over a substrate-green board rect, scaled to fit, dimensions in mm. Inline
/// SVG with a faint grid + port markers + a board-size caption.
pub fn board_svg(layout: &Layout) -> String {
    const MM: f64 = 1.0e3;
    let margin = 1.5; // mm
    let min_x = layout.bbox.min.x * MM - margin;
    let min_y = layout.bbox.min.y * MM - margin;
    let w_mm = layout.bbox.width() * MM + 2.0 * margin;
    let h_mm = layout.bbox.height() * MM + 2.0 * margin;

    let mut s = String::new();
    s.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"{min_x:.3} {min_y:.3} {w_mm:.3} {h_mm:.3}\" preserveAspectRatio=\"xMidYMid meet\" role=\"img\" aria-label=\"dimensioned edge-coupled board top view\">\n"
    ));

    // Board (substrate) rect.
    s.push_str(&format!(
        "<rect x=\"{:.3}\" y=\"{:.3}\" width=\"{:.3}\" height=\"{:.3}\" fill=\"{SUBSTRATE}\" stroke=\"{SUBSTRATE_EDGE}\" stroke-width=\"0.15\" rx=\"0.6\"/>\n",
        layout.bbox.min.x * MM - margin,
        layout.bbox.min.y * MM - margin,
        w_mm,
        h_mm
    ));

    // Faint grid (~ every 5 mm).
    let step = 5.0;
    let mut gx = (min_x / step).ceil() * step;
    while gx < min_x + w_mm {
        s.push_str(&format!(
            "<line x1=\"{gx:.2}\" y1=\"{min_y:.2}\" x2=\"{gx:.2}\" y2=\"{:.2}\" stroke=\"#ffffff\" stroke-opacity=\"0.05\" stroke-width=\"0.1\"/>\n",
            min_y + h_mm
        ));
        gx += step;
    }
    let mut gy = (min_y / step).ceil() * step;
    while gy < min_y + h_mm {
        s.push_str(&format!(
            "<line x1=\"{min_x:.2}\" y1=\"{gy:.2}\" x2=\"{:.2}\" y2=\"{gy:.2}\" stroke=\"#ffffff\" stroke-opacity=\"0.05\" stroke-width=\"0.1\"/>\n",
            min_x + w_mm
        ));
        gy += step;
    }

    // Copper traces.
    for poly in &layout.traces {
        let pts: Vec<String> = poly
            .verts
            .iter()
            .map(|p| format!("{:.4},{:.4}", p.x * MM, p.y * MM))
            .collect();
        s.push_str(&format!(
            "<polygon points=\"{}\" fill=\"{COPPER}\" stroke=\"{COPPER_EDGE}\" stroke-width=\"0.08\"/>\n",
            pts.join(" ")
        ));
    }

    // Port markers.
    for port in &layout.ports {
        s.push_str(&format!(
            "<circle cx=\"{:.4}\" cy=\"{:.4}\" r=\"{:.4}\" fill=\"none\" stroke=\"{ACCENT}\" stroke-width=\"0.18\"/>\n",
            port.at.x * MM,
            port.at.y * MM,
            (port.width_m * MM * 0.45).max(0.6)
        ));
    }

    // Board-size caption (bottom-left, inside the margin).
    s.push_str(&format!(
        "<text x=\"{:.2}\" y=\"{:.2}\" fill=\"{TEXT}\" font-size=\"1.6\" font-family=\"monospace\">{:.1} × {:.1} mm</text>\n",
        min_x + 0.6,
        min_y + h_mm - 0.8,
        layout.bbox.width() * MM,
        layout.bbox.height() * MM
    ));

    s.push_str("</svg>\n");
    s
}

/// Render the lumped-LC board top-view from the real [`LumpedBoard`]: the
/// substrate, the ground rail + signal line + every SMD pad as copper polygons
/// (all to scale, in mm), ref-des labels at each placement centre, port markers,
/// and overall dimension callouts (board width × height). Inline SVG, design-
/// system palette.
pub fn lumped_board_svg(board: &LumpedBoard) -> String {
    const MM: f64 = 1.0e3;
    let layout = &board.layout;
    // Extra bottom margin for the width dimension line + label.
    let margin = 1.2; // mm sides/top
    let dim_pad = 2.4; // mm bottom (dimension line room)
    let min_x = layout.bbox.min.x * MM - margin;
    let min_y = layout.bbox.min.y * MM - margin;
    let bw = layout.bbox.width() * MM;
    let bh = layout.bbox.height() * MM;
    let w_mm = bw + 2.0 * margin;
    let h_mm = bh + margin + dim_pad;

    let mut s = String::new();
    s.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"{min_x:.3} {min_y:.3} {w_mm:.3} {h_mm:.3}\" preserveAspectRatio=\"xMidYMid meet\" role=\"img\" aria-label=\"dimensioned lumped-LC board top view with SMD footprints\">\n"
    ));

    // Substrate rect (board outline).
    s.push_str(&format!(
        "<rect x=\"{:.3}\" y=\"{:.3}\" width=\"{:.3}\" height=\"{:.3}\" fill=\"{SUBSTRATE}\" stroke=\"{SUBSTRATE_EDGE}\" stroke-width=\"0.12\" rx=\"0.4\"/>\n",
        layout.bbox.min.x * MM,
        layout.bbox.min.y * MM,
        bw,
        bh
    ));

    // Faint grid (~ every 2 mm — the board is small).
    let step = 2.0;
    let mut gx = (min_x / step).ceil() * step;
    while gx < min_x + w_mm {
        s.push_str(&format!(
            "<line x1=\"{gx:.2}\" y1=\"{:.2}\" x2=\"{gx:.2}\" y2=\"{:.2}\" stroke=\"#ffffff\" stroke-opacity=\"0.05\" stroke-width=\"0.06\"/>\n",
            layout.bbox.min.y * MM,
            layout.bbox.min.y * MM + bh
        ));
        gx += step;
    }
    let mut gy = (layout.bbox.min.y * MM / step).ceil() * step;
    while gy < layout.bbox.min.y * MM + bh {
        s.push_str(&format!(
            "<line x1=\"{:.2}\" y1=\"{gy:.2}\" x2=\"{:.2}\" y2=\"{gy:.2}\" stroke=\"#ffffff\" stroke-opacity=\"0.05\" stroke-width=\"0.06\"/>\n",
            min_x, min_x + w_mm
        ));
        gy += step;
    }

    // All copper (ground rail, signal-line segments, pads) — to scale.
    for poly in &layout.traces {
        let pts: Vec<String> = poly
            .verts
            .iter()
            .map(|p| format!("{:.4},{:.4}", p.x * MM, p.y * MM))
            .collect();
        s.push_str(&format!(
            "<polygon points=\"{}\" fill=\"{COPPER}\" stroke=\"{COPPER_EDGE}\" stroke-width=\"0.04\"/>\n",
            pts.join(" ")
        ));
    }

    // Ref-des labels at each placement centre (e.g. L1, C1, …).
    for p in &board.placements {
        let (cx, cy) = p.center_m;
        s.push_str(&format!(
            "<text x=\"{:.3}\" y=\"{:.3}\" fill=\"{TEXT}\" font-size=\"0.9\" font-family=\"monospace\" text-anchor=\"middle\" paint-order=\"stroke\" stroke=\"#0b0d11\" stroke-width=\"0.18\">{}</text>\n",
            cx * MM,
            cy * MM + 0.3,
            p.ref_des
        ));
    }

    // Port markers at the two signal-line ends.
    for port in &layout.ports {
        s.push_str(&format!(
            "<circle cx=\"{:.4}\" cy=\"{:.4}\" r=\"{:.4}\" fill=\"none\" stroke=\"{ACCENT}\" stroke-width=\"0.12\"/>\n",
            port.at.x * MM,
            port.at.y * MM,
            (port.width_m * MM * 0.5).max(0.4)
        ));
    }

    // Width dimension line + label along the bottom.
    let dim_y = layout.bbox.min.y * MM + bh + 1.2;
    let x0 = layout.bbox.min.x * MM;
    let x1 = layout.bbox.min.x * MM + bw;
    s.push_str(&format!(
        "<line x1=\"{x0:.3}\" y1=\"{dim_y:.3}\" x2=\"{x1:.3}\" y2=\"{dim_y:.3}\" stroke=\"{MUTED}\" stroke-width=\"0.06\"/>\n"
    ));
    for x in [x0, x1] {
        s.push_str(&format!(
            "<line x1=\"{x:.3}\" y1=\"{:.3}\" x2=\"{x:.3}\" y2=\"{:.3}\" stroke=\"{MUTED}\" stroke-width=\"0.06\"/>\n",
            dim_y - 0.5,
            dim_y + 0.5
        ));
    }
    s.push_str(&format!(
        "<text x=\"{:.3}\" y=\"{:.3}\" fill=\"{TEXT}\" font-size=\"1.0\" font-family=\"monospace\" text-anchor=\"middle\">{bw:.1} mm</text>\n",
        (x0 + x1) / 2.0,
        dim_y + 1.5
    ));

    s.push_str("</svg>\n");
    s
}
