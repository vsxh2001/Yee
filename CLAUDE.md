# Yee — Project Memory for Claude Code

This file is loaded by Claude Code instances starting work in the Yee repo. It captures conventions, patterns, and decisions that have accreted over the project so they don't have to be rediscovered each session. Update this file when a recurring gotcha bites for the second time.

---

## 1. Project overview

**Yee** is an open, GPU-accelerated electromagnetic simulator written in Rust, with a planar Method-of-Moments (MoM) beachhead, a complementary FDTD volumetric solver, Touchstone I/O, Python bindings, and a desktop GUI. Development is **phase-driven** per `ROADMAP.md`: Phase 0 (walking skeleton) is complete; Phase 1.0 (free-space MoM half-wave dipole) has shipped against the published NEC-4 reference; multiple Phase 1 / Phase 2 sub-projects have landed (multilayer Green's placeholder, wave-port skeleton, FDTD CPML / NTFF / dispersive ADE materials, cuSOLVER LU, Python bindings, mdBook theory chapters, egui+wgpu desktop shell with S-parameter and Smith-chart plots, static plotters export, hardware-gated GPU nightly CI). The current shipped solver accuracy floor is the MoM dipole; everything else is either a hardware-gated path, a Phase-1.x placeholder, or an FDTD building block awaiting an end-to-end driver.

**When starting work on this repo:**

1. Read this file end to end.
2. Skim `ROADMAP.md` for the phase your task lives in (core solvers). Application phases live in `FILTER-DESIGN-ROADMAP.md`; the GPU/CPU engine + web-studio track (ADR-0175) lives in `ENGINE-STUDIO-ROADMAP.md`.
3. Skim `TECH_STACK.md` if your task touches a new dependency.
4. Look for a matching spec under `docs/superpowers/specs/` and plan under `docs/superpowers/plans/`. If none exists and the task is non-trivial, write the spec before any code.
5. Decide your lane (§6) before opening any file.

---

## 2. Workspace layout

```
crates/
  yee-core/       — shared types, units, error
  yee-cuda/       — cudarc wrapper, cuSOLVER LU (Zgetrf / Zgetrs)
  yee-mesh/       — Gmsh FFI (gmsh feature)
  yee-mom/        — planar Method of Moments solver
  yee-fdtd/       — FDTD walking skeleton + CPML + NTFF + dispersive ADE materials
  yee-compute/    — GPU/CPU execution layer: rayon CPU + wgpu/WGSL compute (ADR-0175)
  yee-engine/     — transport-agnostic simulation job API over yee-compute (S.0, ADR-0179)
  yee-server/     — axum WebSocket exposure of the job API + `yee serve` (S.1, ADR-0180)
  yee-io/         — Touchstone v1.1 I/O
  yee-cli/        — yee CLI (validate / mesh / run / export / plot)
  yee-py/         — PyO3 0.28 Python bindings (abi3-py310)
  yee-gui/        — egui desktop shell + wgpu 3D viewport
  yee-plotters/   — static PNG/SVG plot export (plotters)
examples/         — 3 runnable example binaries (half-wave-dipole, microstrip-line, patch-2g4)
studio/           — Tauri 2 + React studio (own cargo workspace — NOT a root-workspace member; ADR-0179)
docs/             — mdBook (theory + tutorials) + superpowers/specs + superpowers/plans
.github/workflows/ — CI + GPU nightly + wheels + docs deploy
```

Other root files worth knowing: `ROADMAP.md`, `FILTER-DESIGN-ROADMAP.md`, `ENGINE-STUDIO-ROADMAP.md`, `TECH_STACK.md`, `CONTRIBUTING.md`, `THIRD_PARTY_LICENSES.md`, `rust-toolchain.toml` (pins 1.92), `rustfmt.toml`, `Cargo.toml` (workspace).

`crates/yee-surrogate/` has landed (Phase 3.gp.0/1 + 3.bo.0/1 + 3.al.0 shipped per `ROADMAP.md`). It is wired into the workspace `Cargo.toml` and exposed via `yee-py`'s `yee.surrogate` Python module.

---

## 3. Conventions

- **Rust 1.92+**, pinned in `rust-toolchain.toml`. Bumped from 1.88 in Phase 1.gui.3 (2026-05-17) alongside egui 0.34 / wgpu 29; do not bump casually beyond this.
- **All public items documented.** Each crate sets `#![warn(missing_docs)]`.
- **`#![forbid(unsafe_code)]` is the default.** It is relaxed only inside FFI submodules with an explicit `#[allow(unsafe_code)]` comment:
  - `yee-mesh` — Gmsh C-API FFI
  - `yee-cuda` — cudarc / cuSOLVER raw bindings
  No other crate should contain `unsafe`. If a new FFI need arises, document the reason inline.
- **Feature flags default OFF for anything requiring an external toolchain.** `cuda` (yee-cuda), `gmsh` (yee-mesh), and similar features must build green without the toolchain present; the no-feature path returns a `NotEnabled` error or a stub.
- **Walking-skeleton first.** Ship the minimal end-to-end pipe before optimizing or generalizing. This is non-negotiable for any new sub-system. Concretely: a Phase-X.0 placeholder that compiles and exercises the I/O contract beats a half-finished Phase-X.1.
- **No solver feature ships without a published-benchmark validation case** — see §4.
- **Sub-projects are decomposed before agents are dispatched.** Each non-trivial sub-project gets:
  - a spec under `docs/superpowers/specs/<date>-<name>-design.md`
  - an implementation plan under `docs/superpowers/plans/<date>-<name>.md`

  See §5 for the multi-track orchestration pattern that consumes these.

- **Lint floor:** `cargo clippy --workspace --all-targets -- -D warnings` and `cargo fmt --check --all` are CI gates. Run them locally before opening a PR.

### Commit-message style

The repo's commit log has a consistent shape; match it:

- **Crate-scoped, lowercase prefix.** `yee-fdtd: ADE update kernels for Drude/Lorentz/Debye`. `yee-mom: solve.rs takes &dyn Port instead of port_tag`. `ci: install libfontconfig1-dev for plotters build (fix CI)`. `docs: root CLAUDE.md — workspace conventions, multi-track orchestration, validation gates`.
- **Merge commits get track identifiers and phase tags:** `Merge Track T: Phase 1.5 cuSOLVER Zgetrf/Zgetrs (hardware-gated)`, `Merge Track V: Phase 2.fdtd.3 dispersive materials (Drude/Lorentz/Debye ADE)`. Track letters identify the parallel lane; phase tags map back to `ROADMAP.md`.
- **Subject line ≤ 72 chars, no trailing period, imperative mood.** Body wraps at 80 and explains _why_ not _what_.
- **Co-authored-by trailers** are added by the agent tooling when applicable; do not strip them.

---

## 4. Validation gates

> **No solver feature ships without a published-benchmark validation case in `crates/<crate>/validation/` or `crates/<crate>/tests/`.** CI gates per crate must pass before merge.

### mom-001 — half-wave dipole

- Geometry: L = 1 m, cylinder radius a = 5 mm, delta-gap excitation at the centre.
- **Reference: NEC-4 finite-radius `Z ≈ 87 + j41 Ω`.** Tolerance ±5% on Re(Z), ±10% on Im(Z).
- This is **not** the Balanis wire-limit `73 + j42 Ω` — that value is the zero-radius / infinitely-thin-wire approximation and disagrees with any real finite-radius solver by ~20% on resistance. Quote NEC-4 and only NEC-4 in commit messages, plots, and docs for mom-001.
- The gate test runs ~7-8 min wall-time on a 24×176 cylinder mesh in release mode; budget accordingly in CI matrices.

### mom-002 / mom-003

- `mom-002` (microstrip Z₀) and `mom-003` (2.4 GHz patch resonance) run with **loose tolerances**. The multilayer Sommerfeld Green's kernel has **shipped** (Phase 1.1.1.2.2); the residual is the **port excitation**, not the kernel — see §10.

### CPML reflection

- Phase-2 FDTD CPML target: **≥30 dB reduction vs PEC** for a plane wave at normal incidence. Currently gated by `crates/yee-fdtd/tests/cpml_reflection.rs`.

### Touchstone round-trip

- `crates/yee-mom/tests/touchstone_roundtrip.rs` enforces lossless write→read fidelity on `.s1p`/`.s2p`. Do not weaken this test; Touchstone is the project's primary external interface.

---

## 5. Multi-track orchestration pattern

This is unique to Yee and worth getting right. The canonical brief template lives in `docs/superpowers/specs/2026-05-16-phase-0-multi-agent-execution-design.md`; read it before dispatching anything non-trivial.

- **Worktree per substantive track.** Layout: `worktrees/<lane>/` with branch `feature/<phase>-<name>`. Disjoint sub-projects can run concurrently because their worktrees are physically separate checkouts of the same repo, each holding its own working tree state. Track letters in commit messages (e.g. `Track T`, `Track W`, `Track V`) correspond to lanes that ran in parallel.
- **Up to 5 parallel agents** has been observed feasible on this repo without coordination overhead dominating throughput. Beyond 5 the merge train backs up: review latency rises and `Cargo.lock` conflicts compound.
- **Each agent brief contains, at minimum:**
  - **WORKTREE / BASE COMMIT** — the SHA the worktree was forked from. Surface this in the report so reviewers can reconstruct the diff base.
  - **LANE** — the allowed paths (see §6). Out-of-lane edits are findings, not fixes.
  - **DoD** — concrete, verifiable, machine-checkable success criteria. Every item must be checkable by a shell command with a known exit code or by a `grep`-able artifact.
  - **PATTERN FILE** — an existing-in-repo example to imitate so the agent picks up house style without re-deriving it (e.g. point at `crates/yee-fdtd/tests/cpml_reflection.rs` when asking for a new integration test).
  - **VERIFICATION COMMAND** — exact shell command plus expected exit code. The agent must run this before declaring done.
  - **ESCAPE HATCH** — a "stop and surface" threshold. The standard form is "blocked > 15 min → surface and stop." This prevents runaway grinding on the wrong problem.
- **Branch-divergence artifacts in `git diff --stat main..HEAD` are normal** when a worktree was created against a base SHA older than current `main`. The three-way merge handles it cleanly; do not panic and do not rebase out of habit. A `git diff base...HEAD` (three dots) is usually what you actually want for review.
- **`Cargo.lock` conflicts during merges** are the most common merge hazard. The standard resolution is:
  ```bash
  git checkout --theirs Cargo.lock
  cargo check --workspace
  git add Cargo.lock
  git commit --no-edit
  ```
  Re-resolving by hand is almost always wrong: it tends to pin transitives to stale versions that no longer satisfy the new direct-dep constraints, and `cargo check` rejects the result anyway.
- **Merge order matters.** Land foundational crates (`yee-core`, `yee-cuda` backend trait) before features that depend on them. When two tracks both touch a shared crate, the second one in is responsible for the rebase / re-test.

---

## 6. The "lane" concept

Each agent brief includes a **LANE** — the set of paths the agent is permitted to touch. Out-of-lane edits should be **surfaced as a finding in the agent's report, not fixed in place.** This is what makes parallel sub-project execution safe: it lets two tracks share a base SHA without merge-conflicting each other's work.

Examples of lanes that have been used on this repo (each line is a real prior agent brief):

- `crates/yee-mom/src/**, crates/yee-mom/tests/**` — MoM physics lane (mom-001 dipole, multilayer Greens, ports)
- `crates/yee-fdtd/**` — FDTD lane (CPML, NTFF, dispersive ADE materials)
- `crates/yee-gui/**, crates/yee-plotters/**` — frontend lane (egui shell, S-parameter plots, static PNG export)
- `crates/yee-py/**, examples/**/*.py` — Python-bindings lane (PyO3 bindings, notebook helpers)
- `docs/**` — documentation lane (this file's lane; theory chapters, tutorials)
- `.github/workflows/**` — CI lane (`ci.yml`, `gpu-nightly.yml`, `publish-wheels.yml`, `docs.yml`)
- `crates/yee-cuda/**` — GPU lane (cudarc backend trait, cuSOLVER bindings, `cuda` feature gate)

If a needed change crosses a lane (e.g. a new MoM API forces a change in `yee-cli`), the agent should either:

(a) define the API at the boundary so the cross-lane consumer can be updated in a follow-up PR — this is the strong preference; or
(b) call it out explicitly in the report and stop, so the dispatcher can either widen the lane or open a separate ticket.

Silently editing out-of-lane is the failure mode this section exists to prevent.

---

## 7. Toolchain installs

Pre-flight installs called out in implementation plans. None of these are auto-detected; if a feature gate is being exercised, install the dependency first.

```bash
# Rust 1.92 (pinned in rust-toolchain.toml)
curl -sSf https://sh.rustup.rs | sh -s -- --default-toolchain 1.92

# Faster rebuilds across worktrees
cargo install sccache --locked

# mdBook for docs/
cargo install mdbook --locked

# plotters native deps (Linux)
sudo apt install libfontconfig1-dev pkg-config

# Python venv + maturin for yee-py
uv venv .venv
uv pip install maturin pytest numpy
```

Optional, feature-gated:

- **Gmsh SDK 4.15+** for the `gmsh` feature on `yee-mesh`: download from <https://gmsh.info>, set `$GMSH_SDK_ROOT` to the unpacked SDK root before `cargo build --features gmsh`.
- **CUDA Toolkit 12.4+** for the `cuda` feature on `yee-cuda`: cuSOLVER tests are hardware-gated and run only on the GPU nightly runner.

---

## 8. CI/CD layout

- `.github/workflows/ci.yml` — Rust workspace lint + test on Linux + Rust 1.92. The no-default-features lint-test job runs `cargo check`, `cargo clippy -- -D warnings`, `cargo fmt --check`, `cargo test --workspace --no-default-features` (heavy solver gates — mom-001 dipole, the FDTD gates, and the FEM-eigen / heavy-FDTD-validation gates — are `#[ignore]`'d here and run as separate `--release` gate jobs: `MoM dipole gate`, the several `FDTD … gate` jobs, `FEM eigen gates (fem-eig-001/003, release)`, and `FDTD heavy validation gates (release)`; see §10), and `cargo doc --no-deps`. Includes a `python-bindings` job that runs `maturin develop` and `pytest` against `yee-py`. Installs `libfontconfig1-dev` and `pkg-config` so the plotters-backed crates build.
- `.github/workflows/gpu-nightly.yml` — self-hosted GPU runner. **Gated by repo variable `YEE_GPU_RUNNER_ENABLED`**; the workflow no-ops if the variable is unset, so a fork without GPU hardware will not see red nightly runs. This is where `--features cuda -- --include-ignored` tests actually execute.
- `.github/workflows/publish-wheels.yml` — builds Python wheels on tag push (`v*`) via `maturin` with `manylinux_2_28`. PyPI publish step is **commented out** until the maintainer adds a `PYPI_API_TOKEN` repo secret; uncomment when releasing.
- `.github/workflows/docs.yml` — builds the mdBook with `mdbook build docs/` and deploys to GitHub Pages on every push to `main`. **Requires the repo's Pages settings to have Source: GitHub Actions**; otherwise the deploy step fails with `404 Not Found`. This is the single most common "first time setting up the repo" failure.

---

## 9. Key references

When in doubt, read these first — they answer 80% of questions about what's already decided.

- `ROADMAP.md` — phase-by-phase plan, validation milestones per phase, risks
- `TECH_STACK.md` — dependency choices and rationale (why cudarc, why faer, why egui, etc.)
- `CONTRIBUTING.md` — PR / commit / branch-naming conventions
- `docs/superpowers/specs/` — per-sub-project design specs:
  - `2026-05-16-phase-0-multi-agent-execution-design.md` — the canonical agent-brief template and multi-track pattern
  - `2026-05-16-phase-1-0-free-space-mom-dipole-design.md` — mom-001 deep-dive
  - `2026-05-16-phase-1-frontend-0-python-bindings-design.md` — yee-py shape
- `docs/superpowers/plans/` — per-sub-project step-by-step implementation plans corresponding to the specs above
- `docs/src/theory/planar-mom.md` and `docs/src/theory/fdtd.md` — theory of operation, derivations, and references
- `docs/src/tutorials/01-microstrip-line.md`, `02-dipole-from-python.md`, `03-fdtd-cavity.md` — end-to-end walkthroughs that double as smoke tests

**Why these dependencies, in one line each** (full rationale in `TECH_STACK.md`):

- **cudarc** — pure-Rust CUDA bindings, no `bindgen` build-time CUDA dep. Pre-alpha; pinned to `=0.19.x`; abstracted behind `yee-cuda::backend::Backend`.
- **faer** — pure-Rust dense LA with good performance, used as the CPU-side LU reference and a swap point if `nalgebra-lapack` or `ndarray-linalg` becomes preferable later.
- **Gmsh** (FFI) — best-in-class free mesher; the `rgmsh` crate is unmaintained, so we generate fresh `bindgen` bindings against Gmsh 4.15+.
- **PyO3 0.28** — `abi3-py310` lets one wheel work across Python 3.10+; pairs with `maturin 1.10` for `manylinux_2_28`.
- **egui 0.34 + eframe 0.34 + egui_plot 0.35 + egui_dock 0.19 + wgpu 29** — immediate-mode UI with embedded GPU viewport. egui_plot and egui_dock minor versions track the highest release that pins egui ^0.34 (0.35 / 0.19 at time of bump); wgpu landed on 29 because egui-wgpu 0.34 hard-requires it.
- **plotters** — server-side / headless plot generation for CI artifacts and notebook helpers. Requires `libfontconfig1-dev` on Linux.
- **Touchstone v1.1** — the de-facto S-parameter file format. Our `yee-io` round-trips `.s1p` through generic `.sNp`.

---

## 10. Known limitations and gotchas

These are the things that will bite you if you skip them. Update this section whenever a new gotcha shows up.

- **`cudarc` self-describes as "pre-alpha"** in its own README and has shipped breaking minor releases (notably 0.13 → 0.14). We pin to `=0.19.x`. The internal `Backend` trait in `yee-cuda` exists as the swap point if cudarc ever forces our hand — keep it that way.
- **`MultilayerGreens` (yee-mom): multi-image DCIM + Sommerfeld surface-wave kernel SHIPPED** (Phase 1.1.1.0→1.1.1.2.2; ADRs 0020/0025). Production path = `new_microstrip_sommerfeld` (N-image GPOF fit + TM₀ surface-wave pole subtraction/add-back); every `mom-002`/`mom-003` gate runs through it. The one-image DCIM (`n_images=1`/`n_poles=0`) is now only a **back-compat tripwire**, not the production path. A 10-track forensic effort exonerated the kernel to within ~1.83% of Hammerstad-Jensen ε_eff on mom-002 — **the mom-002/003 residual is the PORT excitation, not the kernel** (ADR-0036/0037; the Phase 1.3.1.2-B numerical quasi-TEM port halved the error to `|Z_in|≈378 Ω` but a cross-section→RWG frame-mapping follow-on remains, ADR-0061). The only un-shipped Greens piece is the full **Sommerfeld-integral tail** (Phase 1.1.1.3 — large/multi-week, no loose gate to validate against; **NOT** the next increment — do not re-scope a "real Greens" track). mom-002/003 stay at loose tolerances because the PORT gates accuracy; do not tighten until a principled port lands.
- **`WavePort` (yee-mom): the cross-section eigenmode solver has SHIPPED.** Phase 1.3.1.1 (closed / slab-loaded guides — FR-4 §4 gate 1.39%) and Phase 1.3.1.2 (quasi-TEM microstrip — HJ 1.2%, ADR-0060) are both in. `ModalDistribution::Numerical2D` + `NumericalCrossSection::with_quasi_tem()` (ADR-0061) inject a numerical modal field, so a microstrip wave-port and a delta-gap excitation now produce **different** results. The mom-002 numerical-port residual is **NOT closable by a frame relabel — it is ill-posed for planar MoM (ADR-0064)**: the microstrip quasi-TEM mode's dominant field is substrate-normal `E_z`, which is orthogonal to the in-plane RWG port-edge tangents and unrepresentable by the in-plane surface-current basis (the planar "port aperture" is a 1-D line, not the 2-D (y,z) face the mode lives on). The current 378 Ω works only because the diagnostic's *wrong* (x,y)-frame cross-section accidentally exposes an in-plane component; a correct (y,z) relabel drives the RHS → 0. **Do NOT re-attempt the Numerical2D microstrip frame mapping.** A true microstrip Z₀ needs a new port formulation (aperture/frill reciprocity, or TL-based Z₀ de-embedding from line currents) — a deferred multi-week track. The `Numerical2D` arm stays correct + validated for *waveguide* ports (WR-90 TE₁₀, port-face = cross-section, in-plane mode). Separately, the fem-eig real-waveguide-port chain is deprioritized (modal-projection saturated — see the fem-eig-006 memory).
- **`fem-eig-006` (Phase 4.fem.eig.3.5.6, ADR-0070, 2026-05-25): Lee-Mittra first-order absorbing-mode complement SHIPPED with escape-hatch.** `PortDefinition::absorbing_complement: bool` + `with_absorbing_complement()` builder (default `false` → backward-compat) implement the Lee-Mittra formula `K = jk₀ B_face + Σ_m j(β_m−k₀) R_m` via `assemble_port_face_block_projected_gauss_pts` (exact Whitney-1, 3-pt Gauss) + `assemble_port_face_block_projected` (centroid path). Measured |S₁₁|(40 GHz) = 0.955500 — essentially unchanged from baseline 0.955397 (0.01% change). Root cause: β₁₀ ≈ 776 and β₂₀ ≈ 554 rad/m < k₀ ≈ 838 rad/m, so j(β_m−k₀) corrections are negative-imaginary and the rank-1 projection R_m covers a small fraction of B_face on the 16×3×2 mesh. `fem_eig_006_magnitude_bounded` stays `#[ignore]`'d; tolerance `< 0.1` **not weakened**. Phase 4.fem.eig.3.5.7 (higher-order absorbing BC — Lee-Mittra §V rational-function extension) queued. **Do NOT reopen fem-eig-006 without a new higher-order BC strategy; do NOT weaken the < 0.1 tolerance.**
- **`mom-001` dipole gate test runs ~7-8 min wall-time** on a 24×176 cylinder mesh in `--release` (60-90 min in debug). As of branch `ci/mom001-release-gate`, `dipole_z_at_resonance` and `dipole_z_diagnostics` are **`#[ignore]`'d**: the default `cargo test --workspace` skips them, so the no-default-features CI **lint-test job no longer times out** (it was being **cancelled at ~2h23m** running them in debug — App.3.0 + brick-1 both went red on this before the fix). mom-001 now runs in a dedicated release-gate CI job, **`MoM dipole gate (mom-001, release)`**, via `cargo test -p yee-mom --release --test dipole -- --ignored dipole_z_at_resonance` (name-filtered to the NEC-4 assertion so it skips the manual `dipole_z_diagnostics` 6-mesh sweep and the pre-existing `dipole_full_sweep` 21-point sweep; mirrors the FDTD release-gate idiom; asserts `Z ≈ 87 + j41 Ω`). Run it locally the same way; the diagnostic is manual-only. For other crates prefer `cargo test -p <crate>`. **`fem_eig_003_wr90_stub_abc` is now `#[ignore]`'d in debug** (branch `ci/heavy-tests-release-gate`, see the dedicated bullet below) — it formerly ran non-ignored and got SIGKILL'd on memory-constrained machines (~7.6 min/solve in release, hours in debug, ~7.4 GB RSS); `--release` (boxed at 12g) avoids both the timeout and the OOM.
- **Heavy FEM + FDTD validation tests are `#[ignore]`'d and run as `--release` CI gates** (branch `ci/heavy-tests-release-gate`; same idiom as mom-001 above). They formerly ran non-ignored under the no-default-features **lint-test** job's `cargo test --workspace --no-default-features`, which **times out in debug** (it had been RED for 10+ commits). The ones moved behind `#[ignore]`: in **`yee-validation/tests/`** — all 4 `fem_eig_001_rectangular_cavity` tests (each a 9720-tet LOBPCG eigensolve) and all 5 `fem_eig_003_wr90_stub_abc` tests (each re-runs the *same* ~72k-tet sparse-LU 50-pt sweep — ~7.6 min/solve **in release**, ~7.4 GB RSS, so OOMs at 6g; box at **12g**); in **`yee-fdtd/tests/`** — `heterogeneous_substrate::percell_eps_r_produces_fresnel_reflection`, `ntff_dipole::ntff_recovers_dipole_pattern_broadside`, and **only** `subgrid_plane_wave_traversal::subgrid_plane_wave_matches_coarse_reference` (its `subgrid_step_runs_through_short_integration_window` ~0.4 s smoke stays non-ignored). The default debug `cargo test --workspace` now reports all of these as `ignored` and the offender files finish in ≈0 s — the **lint-test job no longer times out**. Two dedicated release-gate jobs run them (no `needs:`, parallel with the FDTD/mom gates):
  - **`FEM eigen gates (fem-eig-001/003, release)`** — `cargo test -p yee-validation --release --test fem_eig_001_rectangular_cavity -- --ignored` (4 pass, ~22 s solve) **then** `cargo test -p yee-validation --release --test fem_eig_003_wr90_stub_abc -- --include-ignored --test-threads=1` (all 5; ~38 min total = 5×~7.6 min/solve — **chosen over driver-only** because 5 serial solves fit the gate's time budget and the v3.5.2 PML retune drives the band to ≈[-71.5, -55.6] dB / |S11|≈1e-3 so **all 5 pass**, including the strict continuum-limit gates). **`--test-threads=1` is mandatory**: each solve allocates ~7.4 GB, so libtest's default parallelism runs all 5 at once (~37 GB) and **OOM-kills** the 12g box *and* the ~16 GB GitHub runner — bare `--include-ignored` SIGKILL'd in 38 s. Serialized, peak RSS is one solve (~7.4 GB).
  - **`FDTD heavy validation gates (release)`** — `cargo test -p yee-fdtd --release --test heterogeneous_substrate --test ntff_dipole --test subgrid_plane_wave_traversal -- --include-ignored` (all pass in release: heterogeneous ~21 s, ntff ~4 s, subgrid ~1.4 s — vs the ~360/83/30 s debug times that were the timeout driver; the subgrid strict gate now reports `rel_err 0.0000%` per its C6 passive-fine-grid resolution).
  Run any of them locally **boxed** (`YEE_BOX_DIR=$(pwd) YEE_BOX_MEM=12g YEE_BOX_CPUS=3 scripts/yee-box.sh cargo test -p <crate> --release --test <file> -- --include-ignored`). Do **not** un-`#[ignore]` these into the default `cargo test --workspace`; do not weaken their assertions.
- **`Cargo.lock` merge conflicts: take `--theirs`, then `cargo check --workspace`, then commit.** See §5. Do not hand-merge.
- **GPU-feature tests are hidden behind `--include-ignored` and the `cuda` feature.** A green `cargo test --workspace` on a non-GPU machine is **not** evidence that the GPU path works; only the `gpu-nightly.yml` workflow can certify that. The full invocation is `cargo test --workspace --features cuda -- --include-ignored`.
- **`plotters` will fail to link on a fresh Linux box** without `libfontconfig1-dev` and `pkg-config`. CI installs both; local environments often don't. If `cargo build` on `yee-plotters` (or anything that pulls it in) complains about `fontconfig`, that's the cause.
- **Phase numbering follows ROADMAP, not commit history.** `Phase 1.gui.3`, `Phase 1.1.0`, `Phase 1.3.0`, `Phase 2.fdtd.2`, etc. are meaningful identifiers — keep them in commit messages and spec filenames so future-you can grep for the relevant decision context. Track letters (Track A, Track J, Track T, …) name the parallel lane that delivered the work; they have no inherent ordering relative to each other.
- **`docs/superpowers/specs/` and `docs/superpowers/plans/` are kept in lockstep.** Every spec gets a plan; every plan references its spec. If you find one without the other, fix the gap before dispatching an agent against it.
- **Worktree CWD silently redirects git commits.** When you shell into a `worktrees/<lane>/` directory, plain `git commit` commits to *that branch*, not to `main`. All git mutations should be run with `git -C /home/user/Yee <command>` (the repo root), or explicitly `git checkout main` first. If a fix commit lands on a feature branch instead, cherry-pick it to `main` rather than re-applying by hand: `git -C /home/user/Yee cherry-pick <sha>`. Similarly, `git merge` run inside a worktree merges INTO the worktree's branch — to merge a feature branch into `main`, always run from the main worktree root.
- **Bounded Docker dev container for heavy/FDTD work (ADR n/a; main `6bcd026`).** `docker/Dockerfile.dev` (Rust 1.92 + fontconfig) + `scripts/yee-box.sh` run any cargo command in a cgroup-capped container (`--memory`/`--cpus`, default 12g/4cpu) so a heavy release build or a multi-minute gate (mom-001, fdtd-line-eeff) runs **locally without risking the host** — worst case the container is OOM-killed, the host is untouched. `YEE_BOX_DIR=worktrees/<lane> scripts/yee-box.sh cargo test -p <crate> --release -- --ignored <gate> --nocapture` (cold release build ~25s in-container; FDTD gate ~1-4 min). This **supersedes** the older "never run FDTD/mom-001 locally" guidance — use the box. Build once: `docker build -t yee-dev:1.92 -f docker/Dockerfile.dev .`.
- **`yee-voxel::voxelize_microstrip` z-stack: dielectric fills `k = 0..n_sub`, trace PEC at `k_top = n_sub` (ground at k=0).** A prior version left a **one-cell air gap between the ground plane and the dielectric** (trace at `1+n_sub`, substrate only `1..=n_sub`) — a series air capacitance that dragged the FDTD microstrip ε_eff ~20 % LOW (≈2.5 vs 3.33). This silently broke every downstream FDTD result and took a 7-iteration F1.1b.1 saga to find (ADR-0108). The `voxel_001` gate now pins the corrected z-stack; do not reintroduce the gap. The F1.1b.1 full-wave gate is **propagation-based** (`fdtd-line-eeff-001`: driven line, time-gated phase-velocity → ε_eff vs Hammerstad-Jensen ≤15%; `fdtd-line-eeff-coupled-001`: even/odd ε_eff). The **resonant-split** method was abandoned (no box is both high-Q and non-confining: closed PEC box confines ε_eff / adds cavity modes; CPML open box collapses Q → no peaks).

---

*Last updated: 2026-06-02, post ci/heavy-tests-release-gate (+ the mom-001 release-gate before it): main CI had been RED/cancelled for 10+ commits because the no-default-features lint-test ran heavy solver tests in DEBUG. All are now `#[ignore]`'d and run in dedicated `--release` CI gate jobs — mom-001 dipole; fem-eig-001 (×4) + fem-eig-003 (×5) via a FEM-eigen gate (`--test-threads=1` MANDATORY: 5×~7.4 GB solves OOM a 16 GB runner otherwise); heterogeneous_substrate / ntff_dipole / subgrid via an FDTD-heavy gate. §8 + §10 updated. (Two pre-existing P2 doc nits left for a follow-up: `heterogeneous_substrate.rs` module docstring cites the old 80³ grid/coords vs actual 120³; `fem_eig_003_strict_absorption_floor_gate`'s assert message cites a stale `[-45,-35] dB` window.) — Prior note: post Filter Phase F1.1b.1 (ADR-0108): the first full-wave EM gate in the filter pipeline shipped (`fdtd-line-eeff-001` propagation ε_eff, `afd1eff`), after fixing a `voxelize_microstrip` ground/dielectric one-cell air-gap that had been dragging ε_eff ~20% low; plus a bounded Docker dev container (`scripts/yee-box.sh`, `6bcd026`) that makes heavy/FDTD work runnable locally without crashing the host. See the two new §10 gotchas above. — Prior note: post Phase 1.gui.6 (ADR-0081): VSWR circles on the Smith chart in both `yee-gui` (egui) and `yee-plotters`. `smith_vswr_circle_points` / `smith_vswr_circle_pts` compute constant-VSWR loci at ρ=(VSWR−1)/(VSWR+1); VSWR ∈ {1.5, 2, 3, 5, 10} drawn in muted light-blue-grey before data traces. Key implementation gotcha: the loop-closing element is stashed explicitly as `first=[ρ,0]` and pushed (rather than computing `cos(TAU)` which is not exactly `1.0` in IEEE 754 — this causes an exact-equality close test to fail). Review APPROVED after P1 fixes (bit-identical closure in yee-plotters + missing render smoke test). Update this file whenever a decision is made twice or a gotcha bites twice.*
