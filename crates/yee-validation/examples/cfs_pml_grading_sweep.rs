//! Phase 4.fem.eig.3.5.1 R2 + 4.fem.eig.3.5.2 S2 — CFS-PML
//! grading-parameter ablation sweep.
//!
//! Runs the v3.5.1 §4 ablation grid (1 H1 + 5 H2 + 9 H3 rows) plus the
//! v3.5.2 §5 H4 grid (18 rows: `kappa_max = 2 × m ∈ {3, 4} ×
//! thickness ∈ {12, 14, 16} × alpha_grading_order ∈ {0, 1, 2}`)
//! against fem-eig-003 (WR-90 stub + ZMin CFS-PML) and fem-eig-006
//! (100:10:1 high-aspect + XMax CFS-PML), emitting CSV to stdout.
//! Implements the spec §6 stopping rule: run fem-eig-003 first per
//! row; only run fem-eig-006 if fem-eig-003 worst-case
//! `s11_max_db < -40`.
//!
//! ## Usage
//!
//! Full sweep (33 configurations; worst-case ~3 h `--release` because
//! the H4 `thickness = 16` rows grow the extended mesh ~2.7× over the
//! v3.5.1 baseline):
//!
//! ```bash
//! cargo run -p yee-validation \
//!     --example cfs_pml_grading_sweep --release \
//!     > /tmp/cfs_pml_grading_sweep.csv
//! ```
//!
//! Dry-run (single H1 baseline row only; ~3 min `--release`):
//!
//! ```bash
//! cargo run -p yee-validation \
//!     --example cfs_pml_grading_sweep --release -- --dry-run
//! ```
//!
//! CSV columns (one row per configuration; v3.5.2 appends the new
//! `alpha_grading_order` column at the tail of the v3.5.1 header so
//! pre-v3.5.2 readers can parse the prefix unchanged):
//!
//! ```text
//! hypothesis, kappa_max, m, thickness_cells,
//!   fem_eig_003_s11_min_db, fem_eig_003_s11_max_db,
//!   fem_eig_006_s11_mag,
//!   fem_eig_003_runtime_s, fem_eig_006_runtime_s,
//!   alpha_grading_order
//! ```
//!
//! On the first row where both fixtures retire (`s11_max_db < -40 dB`
//! on fem-eig-003 and `|S_11| < 0.1` on fem-eig-006), the binary emits
//! a final `WINNER,...` row tagged with the same parameters and exits.

use std::env;
use std::time::Instant;

use yee_fem::PmlConfig;
use yee_validation::{
    run_fem_eig_003_wr90_stub_abc_with_config, run_fem_eig_006_high_aspect_pml_with_config,
};

/// One row of the §4 / §5 ablation grid.
#[derive(Clone, Copy, Debug)]
struct Configuration {
    hypothesis: &'static str,
    kappa_max: f64,
    m: usize,
    thickness_cells: usize,
    /// Phase 4.fem.eig.3.5.2 S2: `α_α(d)` polynomial grading order.
    /// `0` on every pre-v3.5.2 H1/H2/H3 row (legacy constant
    /// `α_α(d) ≡ α_max`); swept across `{0, 1, 2}` in the new H4
    /// stage per spec §5.
    alpha_grading_order: usize,
}

impl Configuration {
    fn as_pml_config(&self) -> PmlConfig {
        PmlConfig {
            thickness_cells: self.thickness_cells,
            sigma_max: 0.0,
            alpha_max: 0.0,
            kappa_max: self.kappa_max,
            m: self.m,
            alpha_grading_order: self.alpha_grading_order,
        }
    }
}

/// Build the full ablation grid:
///
/// * **v3.5.1 §4**: H1 (1 row, baseline + per-axis), H2 (5 rows,
///   kappa_max sweep excluding the H1 baseline kappa_max = 5), H3
///   (9 rows, m × thickness sweep). All v3.5.1 rows record
///   `alpha_grading_order = 0` (legacy constant `α_α(d) ≡ α_max`)
///   for reproducibility against the QQQQQQQQQ R3 CSV.
/// * **v3.5.2 §5 H4**: fixes `kappa_max = 2` per the QQQQQQQQQ H2
///   finding that kappa clusters within 1 dB at v3.5.1 grid
///   resolution; sweeps `m × thickness × alpha_grading_order` =
///   `2 × 3 × 3 = 18` rows.
fn build_grid() -> Vec<Configuration> {
    let mut grid = Vec::new();
    grid.push(Configuration {
        hypothesis: "H1",
        kappa_max: 5.0,
        m: 3,
        thickness_cells: 6,
        alpha_grading_order: 0,
    });
    for &kappa in &[1.0_f64, 1.5, 2.0, 3.0, 7.0] {
        grid.push(Configuration {
            hypothesis: "H2",
            kappa_max: kappa,
            m: 3,
            thickness_cells: 6,
            alpha_grading_order: 0,
        });
    }
    for &m in &[2_usize, 3, 4] {
        for &thickness in &[6_usize, 8, 10] {
            grid.push(Configuration {
                hypothesis: "H3",
                kappa_max: 2.0,
                m,
                thickness_cells: thickness,
                alpha_grading_order: 0,
            });
        }
    }
    // Phase 4.fem.eig.3.5.2 §5 H4: alpha_grading_order + extended
    // thickness. 2 × 3 × 3 = 18 configurations.
    for &m in &[3_usize, 4] {
        for &thickness in &[12_usize, 14, 16] {
            for &alpha_grading_order in &[0_usize, 1, 2] {
                grid.push(Configuration {
                    hypothesis: "H4",
                    kappa_max: 2.0,
                    m,
                    thickness_cells: thickness,
                    alpha_grading_order,
                });
            }
        }
    }
    grid
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let dry_run = args.iter().any(|a| a == "--dry-run");

    println!(
        "hypothesis,kappa_max,m,thickness_cells,\
         fem_eig_003_s11_min_db,fem_eig_003_s11_max_db,\
         fem_eig_006_s11_mag,\
         fem_eig_003_runtime_s,fem_eig_006_runtime_s,\
         alpha_grading_order"
    );

    let grid = build_grid();
    let iter: Vec<Configuration> = if dry_run { vec![grid[0]] } else { grid };

    let mut winner: Option<Configuration> = None;
    for cfg in &iter {
        let pml_cfg = cfg.as_pml_config();

        let t003 = Instant::now();
        let r003 = run_fem_eig_003_wr90_stub_abc_with_config(pml_cfg);
        let dt003 = t003.elapsed().as_secs_f64();

        let (s11_min_db, s11_max_db) = match r003 {
            Ok(r) => (r.s11_db_min, r.s11_db_max),
            Err(e) => {
                eprintln!("fem-eig-003 driver error on {cfg:?}: {e}");
                (f64::NAN, f64::NAN)
            }
        };

        let (s11_006_mag, dt006) = if s11_max_db.is_finite() && s11_max_db < -40.0 {
            let t006 = Instant::now();
            let r006 = run_fem_eig_006_high_aspect_pml_with_config(pml_cfg);
            let dt006 = t006.elapsed().as_secs_f64();
            let mag = match r006 {
                Ok(r) => r.s11_magnitude,
                Err(e) => {
                    eprintln!("fem-eig-006 driver error on {cfg:?}: {e}");
                    f64::NAN
                }
            };
            (mag, dt006)
        } else {
            (f64::NAN, 0.0_f64)
        };

        println!(
            "{},{},{},{},{:.4},{:.4},{:.6},{:.2},{:.2},{}",
            cfg.hypothesis,
            cfg.kappa_max,
            cfg.m,
            cfg.thickness_cells,
            s11_min_db,
            s11_max_db,
            s11_006_mag,
            dt003,
            dt006,
            cfg.alpha_grading_order,
        );

        if s11_max_db.is_finite()
            && s11_max_db < -40.0
            && s11_006_mag.is_finite()
            && s11_006_mag < 0.1
        {
            winner = Some(*cfg);
            break;
        }
    }

    if let Some(cfg) = winner {
        println!(
            "WINNER,{},{},{},,,,,,{}",
            cfg.kappa_max, cfg.m, cfg.thickness_cells, cfg.alpha_grading_order,
        );
    } else {
        eprintln!(
            "Phase 4.fem.eig.3.5.2 S2: no configuration retired both fixtures; \
             see CSV for the full ablation row-by-row picture."
        );
    }
}
