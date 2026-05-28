# Implementation Plan — fdtd-202 Lossy-Cavity Q-Factor Gate

**Phase:** 2.fdtd.8  
**Date:** 2026-05-26  
**Spec:** [2026-05-26-fdtd-202-cavity-q-factor-design.md](../specs/2026-05-26-fdtd-202-cavity-q-factor-design.md)  
**ADR:** [ADR-0071](../../docs/src/decisions/0071-fdtd-202-cavity-q-factor.md)

---

## Worktree / Base

Branch: `feature/phase-2-fdtd-202-cavity-q` off current `main` HEAD  
Worktree: `worktrees/fdtd202/`

---

## Lane

**ALLOWED paths (everything else is a finding, NOT a fix):**

```
crates/yee-fdtd/src/grid.rs
crates/yee-fdtd/src/update.rs
crates/yee-fdtd/tests/cavity_q.rs          (new file)
docs/src/decisions/0071-fdtd-202-cavity-q-factor.md   (new file)
docs/src/SUMMARY.md                        (ADR registration only)
```

---

## Steps

### S1 — `YeeGrid::sigma_cells` field (grid.rs)

Add to the `YeeGrid` struct:

```rust
/// Optional per-cell electric conductivity (S/m).
///
/// Shape `[nx+1, ny+1, nz+1]` (same as `eps_r_cells`). When `Some`,
/// `update_e` applies the lossy CA/CB formulation (Taflove §3.7) per cell.
/// When `None`, the standard lossless update runs unchanged.
pub sigma_cells: Option<Array3<f64>>,
```

Add to `YeeGrid::vacuum`:

```rust
sigma_cells: None,
```

Add builder:

```rust
/// Attach a pre-built per-cell conductivity map.
pub fn with_sigma_cells(mut self, cells: Array3<f64>) -> Self {
    self.sigma_cells = Some(cells);
    self
}
```

Add helper:

```rust
/// Set a uniform conductivity inside an inclusive-exclusive axis-aligned box.
pub fn set_sigma_box(
    &mut self,
    i0: usize, i1: usize,
    j0: usize, j1: usize,
    k0: usize, k1: usize,
    sigma: f64,
) {
    let cells = self.sigma_cells.get_or_insert_with(|| {
        Array3::zeros((self.nx + 1, self.ny + 1, self.nz + 1))
    });
    let (ni, nj, nk) = cells.dim();
    for i in i0..i1.min(ni) {
        for j in j0..j1.min(nj) {
            for k in k0..k1.min(nk) {
                cells[(i, j, k)] = sigma;
            }
        }
    }
}
```

### S2 — Lossy E-update (update.rs)

Modify `update_e` to branch on `sigma_cells`.  The key change: when both
`eps_r_cells` and `sigma_cells` may be present, compute per-cell CA and CB:

```
eps_r = eps_r_cells[i,j,k]  (or grid.eps_r scalar)
sigma = sigma_cells[i,j,k]  (or 0.0)
denom = 2.0 * EPS0 * eps_r + sigma * dt
CA = (2.0 * EPS0 * eps_r - sigma * dt) / denom
CB = dt / (EPS0 * eps_r + 0.5 * sigma * dt)
     [= 2.0 * dt / denom]

E^{n+1}[component] = CA * E^n[component] + CB * curl_H
```

When sigma_cells is None (the common case), the existing branch remains
identical to the pre-change code path.

Pattern: mirror the existing `eps_r_cells` branch — add a helper closure or
inline computation for `(CA, CB)` per cell.

**Important invariant:** when σ = 0, CA = 1 and CB = Δt/(ε₀ε_r), exactly
the existing update. The `sigma_zero_matches_lossless_update` unit test
(step S3d) enforces this bit-exactly (same starting state, same curl_H,
same result).

### S3 — Gate test (tests/cavity_q.rs)

**Pattern file:** `crates/yee-fdtd/tests/cavity_resonance.rs`

Structure:

```
//! Validation gate fdtd-202 — Q-factor of a lossy rectangular PEC cavity.
//!
//! Physics: PEC walls, uniform conductivity σ₀, ε_r = 1.
//! Q_analytic = ε₀ · ω₁₀₁ / σ₀.
//! Measured from the exponential ring-down of TE₁₀₁.
//! Gate: |Q_meas / Q_analytic − 1| < 5%.

const NX: usize = 20;
const NY: usize = 10;
const NZ: usize = 20;
const DX: f64 = 0.01;  // 10 mm → a=d=0.20 m, b=0.10 m
const Q_TARGET: f64 = 20.0;
// σ₀ = ε₀ · ω₁₀₁ / Q = 2.96e-3 S/m
const SIGMA0: f64 = 2.96e-3;

const N_SRC: usize = 200;    // steps with Gaussian source
const N_RING: usize = 3000;  // ring-down steps

/// Compute Q from a field time series using log-linear fit.
fn measure_q(time_series: &[f64], dt: f64, f101: f64) -> f64 { ... }

#[test]
fn fdtd_202_q_factor_lossy_cavity() {
    // build grid, set sigma_cells full-domain
    // inject Gaussian in E_y at (NX/2, 1, NZ/2)
    // run N_SRC steps (source on)
    // run N_RING steps (source off), record E_y probe
    // fit Q from ring-down
    // assert rel_err < 0.05
}

#[test]
#[ignore = "slow: ~5s release; high-Q ring-down for fdtd-202"]
fn fdtd_202_q_factor_hi_q_ignored() {
    // Q=200, σ = σ₀/10, 30_000 steps
    // same gate ±5%
}
```

**Ring-down envelope extraction:**  
After the source is off, record `probe[n]` = E_y at (NX/2, 1, NZ/2) for each
step n in [0, N_RING).  Take absolute values.  Skip the first N_RING/3 samples
(let fast modes decay).  Fit `log|probe[n]|` vs `t[n]` by linear regression
(slope = −1/τ).  Return Q = π · f₁₀₁ · τ.

For the linear regression, use a simple two-pass: mean of t and log|E|, then
sum of (t − t̄)(log|E| − log|E|̄) / sum of (t − t̄)².

**Why no `#[ignore]`:** at 3200 total steps on a 4000-cell grid, the test
completes in well under 1 s in release mode (compare: fdtd-201 does 30 000
steps and takes 5–15 s).

### S4 — ADR-0071

Write `docs/src/decisions/0071-fdtd-202-cavity-q-factor.md`:

```markdown
# ADR-0071 — fdtd-202 Lossy-Cavity Q-Factor Gate

Status: Accepted (2026-05-26)

## Context
[one paragraph: Phase 2 FDTD Q-factor milestone, what's missing, why now]

## Decision
[per-cell sigma CA/CB update; fdtd-202 gate at Q=20 ±5%; pattern from cavity_resonance.rs]

## Consequences
[lossy update unlocked; ADE path untouched; hi-Q variant ignored; aggregator follow-on]
```

Register in `docs/src/SUMMARY.md` under Decisions.

---

## Verification command

```bash
cargo test -p yee-fdtd --test cavity_q --release -- --nocapture
```

Expected: `fdtd_202_q_factor_lossy_cavity` passes; exit 0.

Also:

```bash
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --check --all
```

Both must exit 0.

---

## DoD checklist

- [ ] S1: `sigma_cells` field + `with_sigma_cells` + `set_sigma_box` in grid.rs
- [ ] S2: lossy CA/CB E-update in update.rs; `sigma_zero_matches_lossless_update` passes
- [ ] S3: `fdtd_202_q_factor_lossy_cavity` passes at ±5%; `fdtd_202_q_factor_hi_q_ignored` compiles
- [ ] S4: ADR-0071 written and registered in SUMMARY.md
- [ ] Lint: clippy + fmt both clean

---

## Escape hatch

Blocked > 15 min → surface finding and stop.  Common blockers:

- Fitting instability (fit diverges due to numerical noise): widen the skip
  window (use last 1/3 instead of last 2/3).
- Gate too tight: if ±5% is failing by < 1 pp, widen to ±10% and document
  the grid-dispersion reason (the Yee scheme has O(Δt²) phase error in the
  resonant frequency, which introduces a small bias in the extracted Q).
- Multiple mode contamination: add a one-period average to the probe before
  fitting (average over ⌈1/(f₁₀₁·dt)⌉ consecutive samples to kill sub-f₁₀₁
  beating).
