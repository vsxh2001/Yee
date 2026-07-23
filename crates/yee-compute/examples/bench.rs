//! FS.7.0 reproducible CPU-vs-GPU FDTD throughput bench (ADR-0223).
//!
//! Sweeps vacuum cubic grids and reports Mcells·step/s for the rayon CPU
//! backend and the wgpu GPU backend side by side, plus an estimated
//! memory-bandwidth utilization for the GPU path.
//!
//! This is the sync-correct successor to an earlier scratch benchmark that
//! timed a bare `step_n` call: `step_n` only *submits* GPU work (wgpu
//! queues are async), so that naive timing measured submission overhead,
//! not device time (it read a several-hundred-x "speedup" that did not
//! scale with grid size — an artifact, not a result). This bench times
//! `step_n` followed by [`yee_compute::GpuFdtd::sync`] instead — `sync()`
//! exists for exactly this purpose.
//!
//! ```bash
//! cargo run -p yee-compute --release --example bench
//! cargo run -p yee-compute --release --example bench -- --json
//! ```
//!
//! Self-skips (prints `SKIPPED`, exits 0) when no wgpu adapter is present
//! (typical hosted CI) or when built `--no-default-features` (`gpu`
//! feature off) — same posture as `tests/gpu_cpu_parity.rs` (CLAUDE.md §10:
//! a green run here is not proof the GPU path works, only the GPU-hardware
//! runner is).

use std::env;
use std::time::Instant;

use yee_compute::{CpuFdtd, FdtdSpec, Fields};

#[cfg(feature = "gpu")]
use std::panic::{AssertUnwindSafe, catch_unwind};
#[cfg(feature = "gpu")]
use yee_compute::{ComputeError, GpuFdtd};

/// Grid edge lengths (cubes) swept by the bench.
const GRIDS: [usize; 6] = [64, 96, 128, 160, 192, 224];
/// Steps executed once to warm up caches/pipelines before timing.
const WARMUP: usize = 10;
/// Steps timed per repetition.
const STEPS: usize = 200;
/// Repetitions per grid per backend; the median absorbs scheduling noise.
const REPS: usize = 3;
const DX: f64 = 1e-3;

/// Estimated bytes moved per cell per FDTD step, for the reported GB/s
/// column. Model: each of the 6 field components (Ex,Ey,Ez,Hx,Hy,Hz) is
/// updated by reading its own previous value, ~3 neighbouring components
/// for the curl stencil, and one per-cell coefficient (`ca`/`cb` for E,
/// `ch` for H — the vacuum scenario carries no CPML/dispersion terms),
/// then writing the new value: 6 f32 accesses/component (1 self + 3
/// neighbours + 1 coeff + 1 write) x 6 components x 4 bytes =
/// 144 B/cell/step. This is a rough working-set estimate, not a captured
/// memory-traffic profile — treat the derived GB/s as approximate, not a
/// hardware counter reading.
const BYTES_PER_CELL_STEP: f64 = 144.0;

/// Rel-L2 tolerance for the 64³ CPU-vs-GPU sanity check (the
/// `gpu_cpu_parity` idiom, loosened slightly for a single-precision FDTD
/// pass with no CPML/dispersion terms in play).
const SANITY_REL_L2_TOL: f64 = 1e-3;

/// Mcells·step/s from an elapsed wall time.
fn mcells_per_step(n: usize, steps: usize, secs: f64) -> f64 {
    (n * n * n * steps) as f64 / secs.max(1e-12) / 1e6
}

/// L2 norm of a slice.
fn l2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// Median of a small sample (sorts in place).
fn median(mut xs: Vec<f64>) -> f64 {
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    xs[xs.len() / 2]
}

/// One grid's CPU-side result: throughput and, only for the sanity grid,
/// the settled fields.
struct CpuRow {
    mcells: f64,
    fields: Option<Fields>,
}

fn bench_cpu(n: usize, want_fields: bool) -> CpuRow {
    let spec = FdtdSpec::vacuum(n, n, n, DX);
    let init = Fields::with_gaussian_ez(&spec, (n / 2, n / 2, n / 2), 2.0);
    let mut cpu = CpuFdtd::new(spec, init);
    cpu.step_n(WARMUP);
    let secs = median(
        (0..REPS)
            .map(|_| {
                let t0 = Instant::now();
                cpu.step_n(STEPS);
                t0.elapsed().as_secs_f64()
            })
            .collect(),
    );
    CpuRow {
        mcells: mcells_per_step(n, STEPS, secs),
        fields: want_fields.then(|| cpu.fields().clone()),
    }
}

/// One grid's GPU-side outcome.
// `Ran`/`Failed` are only constructed on the `gpu`-feature path; the
// `--no-default-features` build only ever produces `Absent`.
#[cfg_attr(not(feature = "gpu"), allow(dead_code))]
enum GpuRow {
    /// Timed successfully.
    Ran {
        adapter: String,
        mcells: f64,
        readback_ms: f64,
        fields: Option<Fields>,
    },
    /// No wgpu adapter at all (or the `gpu` feature is off at compile
    /// time) — the whole bench should skip, not just this row.
    Absent(String),
    /// This grid size specifically failed to allocate/run; other sizes
    /// still run.
    Failed(String),
}

#[cfg(feature = "gpu")]
fn bench_gpu(n: usize, want_fields: bool) -> GpuRow {
    let spec = FdtdSpec::vacuum(n, n, n, DX);
    let init = Fields::with_gaussian_ez(&spec, (n / 2, n / 2, n / 2), 2.0);
    // wgpu's default uncaptured-error handler panics the process on a
    // validation error (e.g. a bind-group buffer exceeding the adapter's
    // max binding size at very large grids). Build+run each size in a
    // fresh `catch_unwind` scope so one bad size doesn't kill the sweep.
    let outcome = catch_unwind(AssertUnwindSafe(|| -> Result<_, ComputeError> {
        let mut gpu = GpuFdtd::new(spec, init)?;
        let adapter = gpu.adapter_name().to_string();
        gpu.step_n(WARMUP)?;
        gpu.sync()?;
        let mut secs = Vec::with_capacity(REPS);
        for _ in 0..REPS {
            let t0 = Instant::now();
            gpu.step_n(STEPS)?;
            gpu.sync()?;
            secs.push(t0.elapsed().as_secs_f64());
        }
        let mcells = mcells_per_step(n, STEPS, median(secs));
        let t1 = Instant::now();
        let fields = gpu.read_fields()?;
        let readback_ms = t1.elapsed().as_secs_f64() * 1e3;
        Ok((adapter, mcells, readback_ms, fields))
    }));
    match outcome {
        Ok(Ok((adapter, mcells, readback_ms, fields))) => GpuRow::Ran {
            adapter,
            mcells,
            readback_ms,
            fields: want_fields.then_some(fields),
        },
        Ok(Err(ComputeError::NoAdapter)) => GpuRow::Absent("no wgpu adapter found".into()),
        Ok(Err(e)) => GpuRow::Failed(e.to_string()),
        Err(payload) => {
            let msg = payload
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| payload.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "wgpu validation panic".to_string());
            GpuRow::Failed(msg)
        }
    }
}

#[cfg(not(feature = "gpu"))]
fn bench_gpu(_n: usize, _want_fields: bool) -> GpuRow {
    GpuRow::Absent("built with --no-default-features (gpu feature off)".into())
}

struct Row {
    n: usize,
    cpu_mcells: f64,
    gpu_mcells: Option<f64>,
    readback_ms: Option<f64>,
    gbps: Option<f64>,
    note: String,
}

fn main() {
    let json = env::args().any(|a| a == "--json");

    // Cheap presence probe on a tiny grid, separate from the swept sizes:
    // NoAdapter (or the feature being off) means every size will fail
    // identically, so skip the whole bench rather than printing N rows of
    // "n/a".
    if let GpuRow::Absent(reason) = bench_gpu(8, false) {
        println!("SKIPPED: {reason}");
        return;
    }

    let mut rows = Vec::new();
    let mut adapter = String::from("N/A");
    let mut sanity = String::from("not run (64^3 not swept)");

    for &n in &GRIDS {
        let want_fields = n == 64;
        let cpu = bench_cpu(n, want_fields);
        let gpu = bench_gpu(n, want_fields);

        let (gpu_mcells, readback_ms, gbps, note, gpu_fields) = match gpu {
            GpuRow::Ran {
                adapter: name,
                mcells,
                readback_ms,
                fields,
            } => {
                adapter = name;
                (
                    Some(mcells),
                    Some(readback_ms),
                    Some(mcells * BYTES_PER_CELL_STEP / 1e3),
                    String::new(),
                    fields,
                )
            }
            GpuRow::Absent(reason) => (None, None, None, format!("n/a: {reason}"), None),
            GpuRow::Failed(reason) => (None, None, None, format!("n/a: {reason}"), None),
        };

        if want_fields {
            sanity = match (&cpu.fields, &gpu_fields) {
                (Some(cf), Some(gf)) => {
                    let mut worst = 0.0f64;
                    for (name, a, b) in [
                        ("ex", &cf.ex, &gf.ex),
                        ("ey", &cf.ey, &gf.ey),
                        ("ez", &cf.ez, &gf.ez),
                    ] {
                        let diff: Vec<f64> = a.iter().zip(b).map(|(x, y)| x - y).collect();
                        let rel = l2(&diff) / l2(a).max(1e-300);
                        worst = worst.max(rel);
                        eprintln!("sanity {name}: rel L2 = {rel:.3e}");
                    }
                    if worst <= SANITY_REL_L2_TOL {
                        format!("PASS (worst rel L2 = {worst:.3e} <= {SANITY_REL_L2_TOL:.0e})")
                    } else {
                        format!("FAIL (worst rel L2 = {worst:.3e} > {SANITY_REL_L2_TOL:.0e})")
                    }
                }
                _ => "SKIPPED (GPU unavailable at 64^3)".to_string(),
            };
        }

        rows.push(Row {
            n,
            cpu_mcells: cpu.mcells,
            gpu_mcells,
            readback_ms,
            gbps,
            note,
        });
    }

    if json {
        print_json(&adapter, &rows, &sanity);
    } else {
        print_table(&adapter, &rows, &sanity);
    }

    if sanity.starts_with("FAIL") {
        std::process::exit(1);
    }
}

fn print_table(adapter: &str, rows: &[Row], sanity: &str) {
    println!("adapter: {adapter}");
    println!(
        "| {:<8} | {:>13} | {:>13} | {:>8} | {:>8} | {:>12} | note |",
        "grid", "CPU Mc/s", "GPU Mc/s", "speedup", "GB/s", "readback ms"
    );
    println!(
        "|----------|---------------|---------------|----------|----------|--------------|------|"
    );
    for r in rows {
        let gpu_s = r
            .gpu_mcells
            .map(|v| format!("{v:.1}"))
            .unwrap_or_else(|| "-".into());
        let speedup = r
            .gpu_mcells
            .map(|v| format!("{:.2}x", v / r.cpu_mcells))
            .unwrap_or_else(|| "-".into());
        let gbps = r
            .gbps
            .map(|v| format!("{v:.0}"))
            .unwrap_or_else(|| "-".into());
        let rb = r
            .readback_ms
            .map(|v| format!("{v:.2}"))
            .unwrap_or_else(|| "-".into());
        println!(
            "| {:<8} | {:>13.1} | {:>13} | {:>8} | {:>8} | {:>12} | {} |",
            format!("{}^3", r.n),
            r.cpu_mcells,
            gpu_s,
            speedup,
            gbps,
            rb,
            r.note
        );
    }
    println!("sanity (64^3 per-E-component rel L2 vs CPU): {sanity}");
}

fn json_escape(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

fn print_json(adapter: &str, rows: &[Row], sanity: &str) {
    println!("{{");
    println!("  \"adapter\": \"{}\",", json_escape(adapter));
    println!("  \"bytes_per_cell_step_assumption\": {BYTES_PER_CELL_STEP},");
    println!("  \"sanity\": \"{}\",", json_escape(sanity));
    println!("  \"rows\": [");
    for (i, r) in rows.iter().enumerate() {
        let comma = if i + 1 < rows.len() { "," } else { "" };
        println!(
            "    {{\"grid\": {n}, \"cpu_mcells_per_s\": {cpu:.3}, \"gpu_mcells_per_s\": {gpu}, \"gbps\": {gbps}, \"readback_ms\": {rb}, \"note\": \"{note}\"}}{comma}",
            n = r.n,
            cpu = r.cpu_mcells,
            gpu = r
                .gpu_mcells
                .map(|v| format!("{v:.3}"))
                .unwrap_or_else(|| "null".into()),
            gbps = r
                .gbps
                .map(|v| format!("{v:.1}"))
                .unwrap_or_else(|| "null".into()),
            rb = r
                .readback_ms
                .map(|v| format!("{v:.3}"))
                .unwrap_or_else(|| "null".into()),
            note = json_escape(&r.note),
        );
    }
    println!("  ]");
    println!("}}");
}
