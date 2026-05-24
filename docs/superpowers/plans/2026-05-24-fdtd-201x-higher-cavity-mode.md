# fdtd-201.x higher-order cavity mode gate — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-24-fdtd-201x-higher-cavity-mode-design.md`
**Base SHA:** `<post-scoping-commit>` (set at dispatch)
**Lane:** `crates/yee-fdtd/tests/**` (new `cavity_higher_mode.rs`, or a fn
in `cavity_resonance.rs`) + `crates/yee-fdtd/validation/README.md` (one row).
NOTHING else.
**Out of lane** (findings, not fixes): `crates/yee-fdtd/src/**` (consume
the public API read-only — if it can't isolate the higher mode, STOP +
surface), the subgrid surface (`subgrid_*`/`berenger_*`/Q6) + fdtd-007,
any `Cargo.toml`.

## Step ladder

### S1 — read the harness
Read `crates/yee-fdtd/tests/cavity_resonance.rs` in full (the fdtd-201 gate
to clone: grid build, PEC, custom step body injecting into `grid.ey`,
single-bin-DFT scan, peak-find, `#[ignore]`, loose-tol doc block).

### S2 — geometry + target mode (a ≠ d, non-degenerate)
Pick a box with a ≠ d (e.g. 24×10×16 cells at dx = 10 mm → a=0.24,
b=0.10, d=0.16 m). Compute the analytic `f_mnp` for the low modes; choose a
target higher-order mode (TE₂₀₁) that is **cleanly separated** (no other
mode within the planned scan band) — verify by listing the nearby
`f_mnp`. Confirm a≠d splits TE₂₀₁ from TE₁₀₂.

### S3 — source/probe placement + scan
Place the Gaussian source + the E_y probe at an antinode of the TARGET
mode (TE₂₀₁ `E_y ∝ sin(2πx/a)·sin(πz/d)` → antinode x≈a/4 or 3a/4, z≈d/2),
off the TE₁₀₁ antinode so the target couples. Set the DFT scan band to
bracket the target + exclude TE₁₀₁ + neighbours. Peak-find; print
extracted f / analytic f / rel error.

### S4 — assert + README
Assert within ±2.5% (documented loose tol; grid dispersion worse at higher
f). `#[ignore]`-gate. Add a README row naming the validated mode + the tol.

## Verification (run in worktree; all exit 0)
```
cargo fmt --check --all
cargo clippy -p yee-fdtd --all-targets -- -D warnings
cargo test -p yee-fdtd --test cavity_higher_mode --release -- --ignored --nocapture   # (or the fn name if added to cavity_resonance.rs)
cargo test -p yee-fdtd                       # rest of suite green (ignores the new slow test)
git diff --stat -- crates/yee-fdtd/src '**/Cargo.toml'   # MUST be empty
```
Capture the extracted-vs-analytic diagnostic line into the report.

## Escape-hatch
- If the target peak is contaminated by a neighbouring/degenerate mode or
  won't couple, adjust geometry/source/probe/band WITHIN the test. If the
  loose ±2.5% can't be met on a reasonable mesh, or the public API can't
  isolate the mode, STOP + surface (no `src/` mode-decomposition).
- Do NOT touch the subgrid surface or fdtd-007. No ±0.5% chase. Blocked
  >15 min → surface + stop. Synchronous; no Monitor/ScheduleWakeup; no
  sub-agents.

## Out-of-scope (findings, not fixes)
* The strict ±0.5% refined-mesh variant; Q-factor extraction.
* Degenerate-cubic single-mode attribution (avoid by a≠d).
