# FS.4.2b implementation plan — per-layer stackup loss

Spec: `docs/superpowers/specs/2026-07-23-fs4-2b-stackup-loss-design.md`
Branch: `feature/fs4.2b-stackup-loss` (from `main` @ `53df0dd`+)
Lane: `crates/yee-voxel/**`, `crates/yee-engine/**`; Task 3 docs
(`docs/src/decisions/0226-*.md`, `docs/src/SUMMARY.md`, `FULL-SUITE-ROADMAP.md`).

## Global constraints (bind every task)

- Real GPU present; GPU evidence prints `adapter 'NVIDIA GeForce RTX 5060 Ti'`.
- After every functional commit: bit-exact suite green (`cargo test -p yee-compute
  --release --test graded_uniform_bitexact --test gpu_graded_parity --test
  gpu_cpu_parity -- --include-ignored`), gate files unmodified; plus
  `cargo test -p yee-voxel --release` and the two stripline gates
  (`engine-stripline-eeff-001`, `stripline_z0`) green unmodified.
- Never weaken any assertion. Loss-off = provable no-op (spec). Honest pins only;
  α measurement > 20 % off closed form → STOP, root-cause, don't widen.
- clippy `-D warnings` workspace + `-p yee-compute --no-default-features`;
  `cargo fmt --check --all`; missing_docs clean; explicit-path staging only.
- Commit style: crate prefix, ≤72-char imperative subject, why-body, trailers:
Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01FFo41x449XDGJ7Xyds4L7M

## Task 1 — `stackup_sigma_cells` (yee-voxel)

- Read `crates/yee-voxel/src/lib.rs`: `voxelize_stackup` k-band layer fill (how ε
  bands map layers to k ranges) and `substrate_sigma_cells` (FS.2c idiom: output
  convention `MaterialsSpec::sigma_cells`, air/PEC handling). Mirror the
  conventions exactly.
- `pub fn stackup_sigma_cells(model: &MicrostripModel, stackup: &Stackup,
  f_ref_hz: f64) -> Vec<f64>`: per-cell σ = 2π·f_ref·ε₀·ε_r·tan δ of the cell's
  layer; 0 for air/metal cells. Derive each cell's layer from the same geometry
  the ε fill used (prefer re-deriving k-bands from `stackup` heights + the model's
  dz over inferring layer from the ε value — ε values can coincide across layers).
- Unit tests (pattern: `sigma_cells_map_substrate_only_at_the_pozar_value`):
  (a) two-layer stackup, distinct (ε_r, tan δ) per layer → exact expected σ in
  each band, boundary k exact; (b) all tan δ = 0 → all-zero vector (the no-op
  guarantee); (c) single-layer consistency vs `substrate_sigma_cells`.
- Verify: `cargo test -p yee-voxel --release`; workspace clippy/fmt. One commit.

## Task 2 — gate `engine-stripline-alpha-001` (yee-engine)

- Pattern file: `crates/yee-engine/tests/stripline_z0.rs` (FS.4.2a — fixture,
  gating, DFT idioms; read its report-worthy header comments). Reuse its fixture
  geometry (ε_r 2.2, b = 16 cells, w/b = 0.8125) with tan δ = 0.02 through
  `stackup_sigma_cells` at f_ref = the gate's grading frequency.
- Extraction: two V-column measurement planes (FS.4.2a Ez-column idiom) separated
  by a known run Δx of many cells, both source-far and end-far; time-gate the
  first forward pass at each plane; α from |V_B/V_A| over Δx at f_ref
  (single-pass ratio of the SAME wave at two planes is launch-normalized by
  construction — each plane sees the identical launch; cite ADR-0204 for why
  cross-run single ratios are the thing to avoid, not this).
- Closed form in-test: α_d[dB/m] = 8.686·(π f √ε_r / c)·tan δ. Assert at f_ref
  within the honest pinned tolerance (target ≤ 10 %); also assert the lossless
  control (tan δ = 0, same fixture) measures |α| near 0 (pin a small bound from
  the measured numeric floor) — the differential kills systematic gating bias.
- Also verify loss-off no-op cheaply: the tan δ = 0 run's gated V-phasor at
  plane A matches the FS.4.2a lossless expectations (non-trivial, sane) — no
  bit-level cross-test-file coupling.
- Runtime: keep ≤ ~2× the Z₀ gate (~24 s); the blanket yee-engine
  `--include-ignored` CI job picks it up automatically (confirmed idiom in
  FS.4.2a — verify, don't assume, with a grep of ci.yml).
- Verify: new gate green with real numbers (α_meas, α_exact, error, lossless
  floor); both stripline gates + bit-exact suite green. One commit.

## Task 3 — ADR-0226 + roadmap row

- `docs/src/decisions/0226-fs42b-stackup-loss.md` per 0225's structure: the σ
  mapping (constant-σ-at-f_ref model + its ∝f deviation off-reference,
  documented), extraction method, measured α vs closed form + lossless floor,
  remaining FS.4.2 scope (MoM cross-check, automesh integration). SUMMARY.md
  line after 0225; FS.4 row append in FULL-SUITE-ROADMAP.md.
