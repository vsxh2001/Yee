# FS.4.2a implementation plan — stripline Z₀ gate

Spec: `docs/superpowers/specs/2026-07-23-fs4-2a-stripline-z0-design.md`
Branch: `feature/fs4.2a-stripline-z0` (from `main` @ `9c23167`+)
Lane: `crates/yee-compute/**`, `crates/yee-engine/**`; Task 3 docs
(`docs/src/decisions/0225-*.md`, `docs/src/SUMMARY.md`, `FULL-SUITE-ROADMAP.md`).

## Global constraints (bind every task)

- Real GPU present (RTX 5060 Ti); any GPU evidence must print
  `adapter 'NVIDIA GeForce RTX 5060 Ti'`.
- After every functional commit: `cargo test -p yee-compute --release --test
  graded_uniform_bitexact --test gpu_graded_parity --test gpu_cpu_parity --
  --include-ignored` green, gate files unmodified. Never weaken any
  assertion/tolerance anywhere. Existing E-probe behavior must be provably
  unchanged (existing probe-using tests green unmodified).
- Honest measurement: a Z₀ gate is pinned only from a measured value with
  understood physics. First measurement > 10 % off the closed form → root-cause
  (staggering, loop placement, gating window), don't widen.
- clippy `-D warnings` (default AND `--no-default-features`), `cargo fmt --check
  --all`, `missing_docs` clean before each commit.
- Commit style `yee-compute:`/`yee-engine:`/`docs:` prefix, ≤72-char imperative
  subject, body explains why, ending with exactly:
Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01FFo41x449XDGJ7Xyds4L7M

## Task 1 — H-component probes (yee-compute)

- Read `crates/yee-compute/src/drive.rs` (Probe, EComponent), `cpu.rs` probe
  recording, `gpu.rs` record_probes pipeline + `shaders/fdtd.wgsl` probe kernel
  BEFORE deciding the shape. Preference order: (a) widen `Probe.component` to a
  new `FieldComponent` enum (E{x,y,z}|H{x,y,z}) if the churn is contained (the
  arena already holds H — flat-offset math exists); (b) else a parallel
  `h_probes: Vec<HProbe>` on `Drive`. Existing `EComponent`-based public API must
  keep compiling for downstream crates OR every downstream use updated in the
  same commit (grep the workspace; yee-engine/yee-fdtd consumers).
- CPU recording exact; GPU: extend the record_probes WGSL path if it is a
  contained edit (H offsets in the same fields arena), else named
  `ComputeError::Unsupported` rejection with a test. If GPU implemented: new
  parity assertion (CPU vs GPU H-probe streams on a small vacuum run,
  rel tol per FP32 idiom of `gpu_cpu_parity`) inside an existing or new test file.
- Unit test: H probe on a known analytic situation (e.g. the E.0 vacuum Gaussian:
  H stream nonzero, antisymmetric where expected) + E-probe regression (an
  existing E-probe test must pass byte-identical — cite which).
- Verify: bit-exact suite; full `cargo test -p yee-compute --release`; clippy both
  configs; fmt. One commit.

## Task 2 — gate `engine-stripline-z0-001` (yee-engine)

- Pattern files: the FS.4.0 stripline fixture/gate (find via
  `grep -rn "stripline" crates/yee-engine/` — `engine-stripline-eeff-001`) for
  board/drive/window idiom; `crates/yee-engine/tests/` release-gate idiom.
- Fixture: symmetric stripline, ε_r = 2.2 (low-dispersion, keeps TEM exact-form
  clean), w/b chosen for Z₀ ≈ 50 Ω region (w/b ≈ 0.65–0.75 at ε_r 2.2 — compute
  from the closed form, don't guess), b ≥ 16 cells (ADR-0215 lesson).
- Closed form in the test (or a small `yee-engine` helper if the gate file gets
  long): K via AGM iteration (documented, ~10 lines), exact k = sech(πw/2b),
  k′ = tanh(πw/2b), Z₀ = η₀/(4√ε_r)·K(k′)/K(k). Cross-check against the
  Wheeler/Pozar fit in a debug assert or test comment (agreement ≲1 %).
- Measurement: V = Σ Ez·Δz column (ground→trace) at a plane ≥ several cells from
  source and ends; I = ∮H·dl rectangular loop around the trace at the same plane
  (Hy top/bottom legs, Hz side legs — mind sign conventions and the ½-cell
  staggering; document the contour in the test header). Time-gate the first
  forward pass (window before end-reflection return — reuse the eeff gate's
  arrival-time bookkeeping). Z₀ from the gated ratio; report the number printed
  with `--nocapture`-friendly eprintln (house idiom).
- Assert: |Z₀_meas − Z₀_exact|/Z₀_exact within the honest pinned tolerance
  (target ≤ 5 %; pin measured + margin). Also assert V and I are non-trivial
  (guard against a silent-zero probe wiring bug making 0/0 or a tiny ratio).
- Runtime budget: keep the grid modest (this is a Z₀ gate, not an S-param sweep);
  if > ~3 min release, mark `#[ignore]` + wire into the existing yee-engine
  release-gate CI step (grep `.github/workflows/ci.yml` for the engine gates
  job and follow its pattern) — same idiom as other engine gates.
- Verify: gate green (real numbers in report); `engine-stripline-eeff-001` still
  green (shared fixture territory); workspace clippy/fmt. One commit (+ CI wiring
  commit if needed).

## Task 3 — ADR-0225 + roadmap row

- `docs/src/decisions/0225-fs42a-stripline-z0.md` (structure per 0224): probe
  design choice, extraction method incl. staggering treatment, measured Z₀ vs
  exact + vs Wheeler fit, gate tolerance rationale, what remains of FS.4.2
  (tan δ, MoM cross-check, automesh rule). SUMMARY.md line after 0224. FS.4 row
  update in FULL-SUITE-ROADMAP.md.
