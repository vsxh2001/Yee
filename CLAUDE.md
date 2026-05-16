# Yee — Project Memory for Claude Code

This file is loaded by Claude Code instances starting work in the Yee repo. It captures conventions, patterns, and decisions that have accreted over the project so they don't have to be rediscovered each session. Update this file when a recurring gotcha bites for the second time.

---

## 1. Project overview

**Yee** is an open, GPU-accelerated electromagnetic simulator written in Rust, with a planar Method-of-Moments (MoM) beachhead, a complementary FDTD volumetric solver, Touchstone I/O, Python bindings, and a desktop GUI. Development is **phase-driven** per `ROADMAP.md`: Phase 0 (walking skeleton) is complete; Phase 1.0 (free-space MoM half-wave dipole) has shipped against the published NEC-4 reference; multiple Phase 1 / Phase 2 sub-projects have landed (multilayer Green's placeholder, wave-port skeleton, FDTD CPML / NTFF / dispersive ADE materials, cuSOLVER LU, Python bindings, mdBook theory chapters, egui+wgpu desktop shell with S-parameter and Smith-chart plots, static plotters export, hardware-gated GPU nightly CI). The current shipped solver accuracy floor is the MoM dipole; everything else is either a hardware-gated path, a Phase-1.x placeholder, or an FDTD building block awaiting an end-to-end driver.

**When starting work on this repo:**

1. Read this file end to end.
2. Skim `ROADMAP.md` for the phase your task lives in.
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
  yee-io/         — Touchstone v1.1 I/O
  yee-cli/        — yee CLI (validate / mesh / run / export / plot)
  yee-py/         — PyO3 0.28 Python bindings (abi3-py310)
  yee-gui/        — egui desktop shell + wgpu 3D viewport
  yee-plotters/   — static PNG/SVG plot export (plotters)
examples/         — 3 runnable example binaries (half-wave-dipole, microstrip-line, patch-2g4)
docs/             — mdBook (theory + tutorials) + superpowers/specs + superpowers/plans
.github/workflows/ — CI + GPU nightly + wheels + docs deploy
```

Other root files worth knowing: `ROADMAP.md`, `TECH_STACK.md`, `CONTRIBUTING.md`, `THIRD_PARTY_LICENSES.md`, `rust-toolchain.toml` (pins 1.88), `rustfmt.toml`, `Cargo.toml` (workspace).

A `yee-surrogate` crate is on the Phase-3 roadmap but **has not landed** as of the current `main` (`33e28db`). Do not assume it exists.

---

## 3. Conventions

- **Rust 1.88+**, pinned in `rust-toolchain.toml`. Do not bump casually — see Phase 1.gui.3 note in §10.
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

- `mom-002` (microstrip Z₀) and `mom-003` (2.4 GHz patch resonance) will run with **loose tolerances** until the real multilayer Green's function lands. See §10 for the `MultilayerGreens` placeholder status.

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
# Rust 1.88 (pinned in rust-toolchain.toml)
curl -sSf https://sh.rustup.rs | sh -s -- --default-toolchain 1.88

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

- `.github/workflows/ci.yml` — Rust workspace lint + test on Linux + Rust 1.88. The default-features matrix runs `cargo check`, `cargo clippy -- -D warnings`, `cargo fmt --check`, `cargo test --workspace`, and `cargo doc --no-deps`. Includes a `python-bindings` job that runs `maturin develop` and `pytest` against `yee-py`. Installs `libfontconfig1-dev` and `pkg-config` so the plotters-backed crates build.
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
- **egui 0.32 + eframe + wgpu 25** — immediate-mode UI with embedded GPU viewport. Pinned versions are forced by the Rust 1.88 floor; see §10.
- **plotters** — server-side / headless plot generation for CI artifacts and notebook helpers. Requires `libfontconfig1-dev` on Linux.
- **Touchstone v1.1** — the de-facto S-parameter file format. Our `yee-io` round-trips `.s1p` through generic `.sNp`.

---

## 10. Known limitations and gotchas

These are the things that will bite you if you skip them. Update this section whenever a new gotcha shows up.

- **`egui` is pinned to 0.32, not the TECH_STACK target 0.34.** The 0.34 release requires `rustc 1.92+`, and the workspace pins `1.88`. Phase 1.gui.3 will bump both in lockstep; until then, do not unilaterally bump egui.
- **`wgpu` was downgraded 26 → 25 to match `egui-wgpu 0.32`'s transitive requirement.** Same Phase 1.gui.3 follow-up. The downgrade is intentional, not a `cargo update` bug.
- **`cudarc` self-describes as "pre-alpha"** in its own README and has shipped breaking minor releases (notably 0.13 → 0.14). We pin to `=0.19.x`. The internal `Backend` trait in `yee-cuda` exists as the swap point if cudarc ever forces our hand — keep it that way.
- **`MultilayerGreens` (yee-mom) is a Phase 1.1.0 placeholder.** It implements a one-image Discrete Complex Image Method (DCIM) approximation only. `mom-002` (microstrip Z₀) and `mom-003` (2.4 GHz patch) will run with **loose tolerances** against this placeholder until Phase 1.1.1 ships real Sommerfeld-integral / multi-image DCIM extraction. Do not tighten the tolerances until then.
- **`WavePort` (yee-mom Phase 1.3.0) is an API placeholder.** Its modal-distribution path currently matches the delta-gap behaviour bit-for-bit. The cross-section eigenmode solver lands in Phase 1.3.1; until then, a microstrip wave-port and a delta-gap excitation will produce identical results, which is expected.
- **`mom-001` dipole gate test runs ~7-8 min wall-time** on a 24×176 cylinder mesh in `--release`. Budget accordingly when sizing CI jobs and when iterating locally; debug builds are dramatically slower and should not be used for this test.
- **`Cargo.lock` merge conflicts: take `--theirs`, then `cargo check --workspace`, then commit.** See §5. Do not hand-merge.
- **GPU-feature tests are hidden behind `--include-ignored` and the `cuda` feature.** A green `cargo test --workspace` on a non-GPU machine is **not** evidence that the GPU path works; only the `gpu-nightly.yml` workflow can certify that. The full invocation is `cargo test --workspace --features cuda -- --include-ignored`.
- **`plotters` will fail to link on a fresh Linux box** without `libfontconfig1-dev` and `pkg-config`. CI installs both; local environments often don't. If `cargo build` on `yee-plotters` (or anything that pulls it in) complains about `fontconfig`, that's the cause.
- **Phase numbering follows ROADMAP, not commit history.** `Phase 1.gui.3`, `Phase 1.1.0`, `Phase 1.3.0`, `Phase 2.fdtd.2`, etc. are meaningful identifiers — keep them in commit messages and spec filenames so future-you can grep for the relevant decision context. Track letters (Track A, Track J, Track T, …) name the parallel lane that delivered the work; they have no inherent ordering relative to each other.
- **`docs/superpowers/specs/` and `docs/superpowers/plans/` are kept in lockstep.** Every spec gets a plan; every plan references its spec. If you find one without the other, fix the gap before dispatching an agent against it.

---

*Last updated: 2026-05-16, base SHA `33e28db` (Track T: Phase 1.5 cuSOLVER Zgetrf/Zgetrs hardware-gated). Update this file whenever a decision is made twice or a gotcha bites twice.*
