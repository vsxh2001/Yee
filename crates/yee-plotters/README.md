# yee-plotters

Static PNG/SVG export of S-parameter plots via the [`plotters`] crate.

This is the "save the plot to a file" companion to the live `egui_plot`
views in `yee-gui`. It is consumed by validation harnesses, examples, and
CI artifact uploads.

## Quick usage

```rust
use std::path::Path;
use num_complex::Complex64;
use yee_plotters::{plot_s11_db, PlotConfig, PlotFormat};

let freq_hz: Vec<f64> = (0..21).map(|i| 1.0e9 + i as f64 * 1.0e8).collect();
let s11: Vec<Complex64> = freq_hz.iter().map(|_| Complex64::new(0.1, 0.0)).collect();
let cfg = PlotConfig { width_px: 800, height_px: 600, title: "S11".into(), format: PlotFormat::Png };
plot_s11_db(&freq_hz, &s11, Path::new("s11.png"), &cfg).unwrap();
```

## Public surface

| Function            | Plot                                                |
| ------------------- | --------------------------------------------------- |
| `plot_s11_db`       | `\|S₁₁\|` in dB vs. frequency (GHz)                 |
| `plot_s11_phase`    | `arg S₁₁` in degrees vs. frequency (GHz)            |
| `plot_smith_chart`  | `S₁₁` on the complex unit disk + reference circle   |

All three accept a `PlotConfig` and dispatch on `PlotFormat::{Png, Svg}`.

[`plotters`]: https://crates.io/crates/plotters
