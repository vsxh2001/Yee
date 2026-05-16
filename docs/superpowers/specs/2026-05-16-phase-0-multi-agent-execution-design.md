# Phase 0 Multi-Agent Execution — Design

**Date:** 2026-05-16
**Status:** Approved by user; ready for writing-plans
**Repo base SHA at design time:** `b6e1c00` (post-scaffold)
**Project name:** Yee

This document specifies the **walking-skeleton Phase 0** delivery plan for the Yee electromagnetic-simulation workspace, executed by a small set of dispatched agents working in isolated git worktrees on a single host.

---

## 1. Scope & Success Criteria

Phase 0 is a **walking skeleton**: every pipe between workspace crates, the CLI, the documentation pipeline, and CI is connected end-to-end. No physical-simulation accuracy is required. The Phase 0 validation milestones previously listed in the top-level `ROADMAP.md` (`mom-001` dipole impedance, `mom-002` microstrip Z₀, `mom-003` 2.4 GHz patch resonance) are explicitly **moved to Phase 1** and recorded as such in this spec.

### "Phase 0 done" means all of the following pass

1. `cargo check --workspace --no-default-features` exits 0
2. `cargo test --workspace --no-default-features` exits 0
3. `cargo clippy --workspace --all-targets -- -D warnings` exits 0
4. `cargo fmt --check --all` exits 0
5. `cargo doc --workspace --no-deps` exits 0 with no warnings
6. `cargo run --bin yee -- --help` exits 0 and lists every subcommand
7. `cargo run --bin yee -- validate all` exits 0 with a stub Phase-0 report
8. `mdbook build docs/` exits 0
9. `THIRD_PARTY_LICENSES.md` present and documents Gmsh, OCCT, and the NVIDIA CUDA proprietary dynamic-link posture
10. CI workflow on `ubuntu-latest` with Rust 1.85 runs gates 1–9 and exits green

### Out of scope for Phase 0

- Physical-solve accuracy (deferred to Phase 1).
- Real CUDA kernel execution (Device enumeration and NVRTC compile helpers must exist; kernel launches are deferred to Phase 1 / Phase 2).
- The `gmsh` feature is built but its full validation is gated on a Gmsh-SDK CI runner that does not yet exist.
- Python bindings (`pyo3`), GUI (`egui`), and surrogate framework — all later phases.

---

## 2. Design Decisions Locked During Brainstorming

| # | Decision | Choice |
|---|----------|--------|
| D1 | Phase 0 success definition | Walking skeleton (Option (i) — pure-build gates, no physical accuracy) |
| D2 | Agent execution environment | This host; install missing toolchain (Rust 1.85, Gmsh SDK, mdBook, sccache); GPU present |
| D3 | Parallelism model | True parallel up to cap-of-three; per-worktree `CARGO_TARGET_DIR`; shared `sccache` |
| D4 | Worktree mechanism | Native `git worktree`; orchestrator manages create/merge/remove |
| D5 | Reviewer policy | One `feature-dev:code-reviewer` dispatch per batch, read-only, before merge |
| D6 | Batching | Approach C — three batches (A solo → B,C,D in parallel → E,F,G in parallel) |

---

## 3. Architecture (already scaffolded at `b6e1c00`)

Cargo workspace with seven member crates plus top-level `examples/`, `validation/`, `docs/`.

```
Cargo.toml                              # workspace, shared deps from TECH_STACK.md
rust-toolchain.toml                     # pins Rust 1.85
crates/
├── yee-core/      # types, traits, units, errors          (no I/O, no CUDA, no GUI)
├── yee-cuda/      # cudarc + NVRTC wrappers               (feature `cuda`)
├── yee-mesh/      # Gmsh FFI via bindgen + safe wrapper   (feature `gmsh`)
├── yee-mom/       # planar MoM solver (Phase 0: skeleton, Phase 1: real)
├── yee-fdtd/      # 3D FDTD solver (Phase 0: stub, Phase 2: real)
├── yee-io/        # Touchstone v1.1, later CAD + HDF5     (feature `opencascade`)
└── yee-cli/       # `yee` CLI binary                      (features `cuda`, `gmsh`)
examples/                               # end-to-end runnable demos (Phase 1+)
validation/                             # workspace-level validation harnesses
docs/                                   # mdBook + ADRs + source/ archive
```

Dependency edges:

```
yee-core ──┬──> yee-mesh ──┐
           ├──> yee-cuda   ├──> yee-mom ──> yee-cli
           └──> yee-io ────┘
            (yee-fdtd stub depends on yee-core + yee-cuda for forward compat)
```

---

## 4. Per-Agent Briefs

Every agent receives the shared preamble below. Agents that exceed 15 min blocked time must surface the specific blocker and stop, not hack around it.

### Shared preamble

```bash
cd <worktree>
export CARGO_TARGET_DIR="$PWD/target"
export RUSTC_WRAPPER=sccache     # no-op if sccache absent
git pull --ff-only origin main   # verify base SHA
```

### Agent A — `yee-core` impl

- **Worktree:** `worktrees/yee-core` from `main @ b6e1c00`
- **Agent type:** `feature-dev:code-architect`
- **Lane:** `crates/yee-core/**` only. Out-of-lane edits → surface as finding, do NOT fix.
- **Definition of done:**
  - `units` module with `C0`, `EPS0`, `MU0`, `ETA0`; tests verify against CODATA 2018.
  - `FreqRange::new(start_hz, stop_hz, n_points)` validates `start < stop`, `n_points >= 1`; otherwise returns `Error::Invalid`.
  - `FreqRange::iter()` yields exact endpoints; `n=1` yields `[start]`; `n=2` yields `[start, stop]`.
  - `Error` enum has variants `Invalid(String)`, `Numerical(String)`, `Unimplemented(&'static str)`, `Io(String)`.
  - Every public item carries `///` documentation, with at least one doc-test on every public function.
  - `cargo doc --no-deps -p yee-core` is warning-free.
- **Pattern files:** `crates/yee-core/src/lib.rs` scaffold; NIST CODATA 2018 reference values.
- **Verification command:**

  ```bash
  cargo test -p yee-core \
    && cargo clippy -p yee-core --all-targets -- -D warnings \
    && cargo doc --no-deps -p yee-core
  ```

### Agent B — `yee-mesh`

- **Worktree:** `worktrees/yee-mesh` from batch-1 merge SHA
- **Agent type:** `feature-dev:code-architect`
- **Lane:** `crates/yee-mesh/**`. Out-of-lane → finding.
- **Definition of done:**
  - `TriMesh` enforces invariant `n_tris == triangles.len() == tags.len()` via its constructor.
  - `build.rs` invokes `bindgen` against `$GMSH_SDK_ROOT/include/gmshc.h` only when feature `gmsh` is enabled and the env var is set; otherwise emits an empty stub.
  - Safe wrapper skeleton (`Session::new`, `Session::drop`, `import_step`, `mesh`, `tris`) returns `Error::NotEnabled` without the feature; with the feature, calls into the generated bindings.
  - Crate builds clean both with and without `--features gmsh`.
- **Pattern files:** existing `crates/yee-mesh/src/lib.rs`; bindgen book.
- **Verification command:**

  ```bash
  cargo build -p yee-mesh \
    && cargo test -p yee-mesh \
    && cargo clippy -p yee-mesh -- -D warnings
  # If Gmsh SDK installed:
  # GMSH_SDK_ROOT=/opt/gmsh-sdk cargo build -p yee-mesh --features gmsh
  ```

### Agent C — `yee-cuda`

- **Worktree:** `worktrees/yee-cuda` from batch-1 merge SHA
- **Agent type:** `feature-dev:code-architect`
- **Lane:** `crates/yee-cuda/**`, including `crates/yee-cuda/kernels/hello.cu`.
- **Definition of done:**
  - `Device { ordinal, name, compute_cap: (u8, u8), mem_total_bytes: u64 }`; fields populated via cudarc when feature `cuda` is on, returns `Error::NotEnabled` otherwise.
  - `Device::list() -> Result<Vec<Device>>` returns visible devices on a CUDA host: ≥ 1 device on a GPU host with the toolkit installed; returns an empty `Vec` (no error) on a toolkit-installed host with no visible GPU; returns `Error::NotEnabled` without the feature.
  - `nvrtc::compile(src: &str, name: &str) -> Result<Vec<u8>>` returns PTX bytes.
  - `kernels/hello.cu` is checked in with kernel signature `__global__ void add_one(float* out, const float* in, int n)`.
  - Internal `backend` trait wraps the cudarc handle so the binding can be swapped later.
  - Crate builds clean with and without `--features cuda`.
- **Pattern files:** existing `crates/yee-cuda/src/lib.rs`; `cudarc` docs.
- **Verification command:**

  ```bash
  cargo build -p yee-cuda \
    && cargo test -p yee-cuda \
    && cargo clippy -p yee-cuda -- -D warnings
  # With CUDA toolkit:
  # cargo build -p yee-cuda --features cuda
  ```

### Agent D — `yee-io` Touchstone v1.1

- **Worktree:** `worktrees/yee-io` from batch-1 merge SHA
- **Agent type:** `feature-dev:code-architect`
- **Lane:** `crates/yee-io/**`, including `crates/yee-io/validation/fixtures/touchstone/**`.
- **Definition of done:**
  - Parser handles option line `# <freq_unit> <param_type> <format> R <Z0>`; `freq_unit ∈ {Hz, kHz, MHz, GHz}`; `param_type = S` only this phase; `format ∈ {RI, MA, DB}`.
  - `.s1p`, `.s2p`, `.s3p`, `.s4p` parse and write correctly.
  - Writer output is deterministic (option line first, comments preserved at the top, space-separated columns).
  - Property test: round-trip of every fixture yields struct equality with floats compared at 1 × 10⁻¹² relative.
  - Passivity check on read: eigenvalues of `S†S ≤ 1 + 1 × 10⁻⁹`.
  - At least three published `.sNp` fixtures checked in with provenance in `fixtures/touchstone/README.md`.
  - Malformed files yield `Error::TouchstoneParse` containing line and column.
- **Pattern files:** Touchstone v1.1 spec; existing `crates/yee-io/src/lib.rs`.
- **Verification command:**

  ```bash
  cargo test -p yee-io \
    && cargo clippy -p yee-io -- -D warnings
  ```

### Agent E — `yee-mom` skeleton wiring

- **Worktree:** `worktrees/yee-mom` from batch-2 merge SHA
- **Agent type:** `feature-dev:code-architect`
- **Lane:** `crates/yee-mom/**`.
- **Definition of done:**
  - `impl Solver for PlanarMoM` returns `yee_core::Error::Unimplemented` from `run()`.
  - `SParameters` is constructable and convertible to and from `yee_io::touchstone::File`.
  - `SParameters::write_touchstone(path)` round-trips through `yee-io`.
  - `cargo test -p yee-mom` is green, including the round-trip test.
  - `PlanarMoM::default()` is constructable; `run()` returns `Unimplemented` without panic.
- **Pattern files:** existing `crates/yee-mom/src/lib.rs`; `yee_io::touchstone::File` shape from Agent D output.
- **Verification command:**

  ```bash
  cargo test -p yee-mom \
    && cargo clippy -p yee-mom -- -D warnings
  ```

### Agent F — `yee-cli` real wiring

- **Worktree:** `worktrees/yee-cli` from batch-2 merge SHA
- **Agent type:** `feature-dev:code-architect`
- **Lane:** `crates/yee-cli/**`, including `crates/yee-cli/tests/**` and `crates/yee-cli/completions/**`.
- **Definition of done:**
  - `yee validate <mom|fdtd|all>` dispatches: `mom` calls a `yee_mom::validation::report()` stub that prints planned cases; `fdtd` prints "Phase 2 deliverable"; `all` invokes both.
  - `yee mesh <path>` constructs a `yee_mesh::Session`; without the `gmsh` feature returns a `NotEnabled` message and exits with code 2.
  - `yee export <input> --format touchstone` round-trips through `yee-io`; `--format hdf5` returns `NotEnabled`.
  - `clap_complete` generates shell completions for bash, zsh, fish into `crates/yee-cli/completions/`.
  - Integration tests via `assert_cmd`: `--help`, `--version`, every subcommand `--help` exit 0; unknown subcommand exits non-zero with suggestion text.
- **Pattern files:** existing `crates/yee-cli/src/main.rs`; `clap_complete` docs.
- **Verification command:**

  ```bash
  cargo test -p yee-cli \
    && cargo run --bin yee -- --help \
    && cargo run --bin yee -- validate all
  ```

### Agent G — Cross-cutting infrastructure

- **Worktree:** `worktrees/yee-infra` from batch-2 merge SHA
- **Agent type:** `general-purpose`
- **Lane:** `.github/**`, `THIRD_PARTY_LICENSES.md`, `docs/**` (excluding `docs/source/**`), top-level `ROADMAP.md` (Phase 0 validation list reconciliation), `.editorconfig`, `rustfmt.toml`.
- **Definition of done:**
  - `.github/workflows/ci.yml` runs fmt, clippy, test, doc jobs on `ubuntu-latest` with Rust 1.85, no features.
  - `.github/workflows/gpu-nightly.yml` is a cron-scheduled stub with `if: false` guarding the actual job, plus inline comments explaining how to flip it on when a self-hosted GPU runner exists.
  - `THIRD_PARTY_LICENSES.md` documents Gmsh GPL v2+ linking exception, OCCT LGPL 2.1 exception, the NVIDIA CUDA proprietary dynamic-link posture, and groups the permissive Rust dependencies.
  - `docs/book.toml`, `docs/src/SUMMARY.md`, and `docs/src/introduction.md` exist; `mdbook build docs/` exits 0.
  - Top-level `ROADMAP.md` is reconciled: the `mom-001 .. mom-003` cases move from the Phase 0 validation list to a Phase 1 footnote; the Phase 0 list becomes the ten "Phase 0 done" gates from §1 of this spec.
  - `.editorconfig` and `rustfmt.toml` are checked in.
- **Pattern files:** `TECH_STACK.md` License-sanity-check section; mdBook docs.
- **Verification command:**

  ```bash
  cargo fmt --check --all \
    && yamllint .github/workflows/ \
    && mdbook build docs/
  ```

---

## 5. Orchestration

### One-time pre-flight (orchestrator runs once)

```bash
# Rust toolchain
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain 1.85
source "$HOME/.cargo/env"
rustup component add rustfmt clippy rust-src

# Shared compile cache
cargo install sccache --locked
export RUSTC_WRAPPER=sccache
export SCCACHE_CACHE_SIZE=10G

# bindgen system dep (Agent B)
sudo apt-get install -y libclang-dev pkg-config

# mdBook (Agent G)
cargo install mdbook --locked

# Gmsh SDK (optional, Agent B feature-gated verify)
# Download 4.15+ from https://gmsh.info/bin/Linux/
# tar xzf gmsh-*-Linux64-sdk.tgz -C /opt/
# export GMSH_SDK_ROOT=/opt/gmsh-sdk

# CUDA toolkit (optional, Agent C feature-gated verify)
# nvcc --version  # expect 12.4+ or 13.0
# nvidia-smi -L
```

### Worktree lifecycle per agent

```bash
# Pre-dispatch
git worktree add -b agent/<lane> worktrees/<lane> <base-sha>

# Agent runs; orchestrator does not poll.

# Post-return
git -C . diff --stat main..agent/<lane>     # verify lane
# Dispatch reviewer (read-only) on this worktree
# After user approves reviewer findings:
git merge --no-ff agent/<lane>
git worktree remove worktrees/<lane>
git branch -d agent/<lane>
```

### Reviewer brief template

```
Type:   feature-dev:code-reviewer  (read-only; no edits)
Worktree: worktrees/<lane>
Base:   main @ <SHA>
Branch: agent/<lane>

1. Lane check.
   git diff --stat main..HEAD
   Expected paths: <glob from §4 brief>
   Any path outside → P1 finding.

2. DoD checklist (paste from §4 brief).

3. Re-run verification command from §4 brief; confirm exit 0.

4. Code-quality smoke: any `unsafe`, library panics, undocumented public items, TODO without owner?

Report:
- P0 (must fix before merge): <list>
- P1 (should fix this sprint): <list>
- P2 (deferrable): <list>
- Out-of-lane edits: <list with verdict — legit / hack / reject>
```

### Batch flow

```
Batch 1
  Dispatch [A]                            (1 agent)
  Reviewer for A                          (1 reviewer)
  User approval → merge → cleanup

Batch 2
  Dispatch [B, C, D] in parallel          (3 agents, single message)
  Reviewer for B + C + D                  (1 reviewer, 3 lanes)
  User approval → merge → cleanup

Batch 3
  Dispatch [E, F, G] in parallel          (3 agents, single message)
  Reviewer for E + F + G                  (1 reviewer, 3 lanes)
  User approval → merge → cleanup

Total: 7 implementation dispatches + 3 reviewer dispatches = 10 agent runs.
```

### Failure responses

| Failure | Response |
|---------|----------|
| Agent exceeds budget | Salvage uncommitted work from worktree, re-dispatch with tighter scope. |
| Agent edits out of lane | Reviewer reports; orchestrator `git checkout main -- <files>`; re-verify; then merge. |
| Build fails on agent worktree | Agent surfaces specific blocker; orchestrator decides install vs scope-cut. |
| sccache disk fills | Bump `SCCACHE_CACHE_SIZE`; `sccache --stop-server && rm -rf ~/.cache/sccache`. |
| GPU contention (Phase 1+ only) | Phase 0 launches no kernels; defer serialization plan to Phase 1. |

---

## 6. Phase 0 Validation Gates

### Repo-level gates (Agent G's CI workflow runs all of these)

| Gate | Command | Expected |
|------|---------|----------|
| `v0-build` | `cargo build --workspace --no-default-features` | exit 0 |
| `v0-test`  | `cargo test --workspace --no-default-features`  | exit 0 |
| `v0-clippy`| `cargo clippy --workspace --all-targets -- -D warnings` | exit 0 |
| `v0-fmt`   | `cargo fmt --check --all` | exit 0 |
| `v0-doc`   | `cargo doc --workspace --no-deps` | exit 0, no warnings |
| `v0-cli`   | `cargo run --bin yee -- --help` | exit 0, lists every subcmd |
| `v0-cli-validate` | `cargo run --bin yee -- validate all` | exit 0 |
| `v0-mdbook`| `mdbook build docs/` | exit 0 |
| `v0-licenses` | grep Gmsh, OCCT, CUDA in `THIRD_PARTY_LICENSES.md` | non-empty matches |

### Per-crate gates

| Crate | Gate ID | What |
|-------|---------|------|
| yee-core | `c0-units` | Constants ≤ 1 × 10⁻¹² relative vs CODATA 2018 |
| yee-core | `c0-freq`  | `FreqRange::iter` endpoint exactness; invalid input → `Error::Invalid` |
| yee-core | `c0-docs`  | Every public item has `///`; each public fn has a doc-test |
| yee-mesh | `c0-mesh-build-default` | Builds without `gmsh` feature |
| yee-mesh | `c0-mesh-build-gmsh`    | Builds with `gmsh` feature when SDK present (optional in CI) |
| yee-mesh | `c0-mesh-trimesh`       | TriMesh invariant verified |
| yee-cuda | `c0-cuda-build-default` | Builds without `cuda` feature |
| yee-cuda | `c0-cuda-build-cuda`    | Builds with `cuda` feature when toolkit present (optional in CI) |
| yee-cuda | `c0-cuda-backend-trait` | Internal `backend` trait present (grep) |
| yee-io   | `c0-io-roundtrip`       | Touchstone fixtures round-trip equal to 1 × 10⁻¹² relative |
| yee-io   | `c0-io-passivity`       | `eig(S†S) ≤ 1 + 1 × 10⁻⁹` on all fixtures |
| yee-io   | `c0-io-malformed`       | Bad files → `TouchstoneParse` with line/column |
| yee-mom  | `c0-mom-unimpl`         | `PlanarMoM::run` returns `Unimplemented` cleanly |
| yee-mom  | `c0-mom-sparams-io`     | `SParameters ↔ touchstone::File` round-trip |
| yee-fdtd | `c0-fdtd-stub`          | Crate compiles; public API returns `Error::NotYet` |
| yee-cli  | `c0-cli-subcommands`    | Every subcmd `--help` exits 0; unknown cmd → non-zero + suggestion |
| yee-cli  | `c0-cli-completions`    | Generates bash / zsh / fish completions |

### Deferred to Phase 1 (was previously listed under Phase 0)

- `mom-001` half-wave dipole impedance ±5%
- `mom-002` 50 Ω microstrip Z₀ ±3%
- `mom-003` 2.4 GHz patch resonance ±2%

Agent G reconciles `ROADMAP.md` accordingly.

### Hardware-gated cases (acceptable to skip in Phase 0)

- `cuda-002` hello kernel launch on real GPU — workflow stub `if: false` until a GPU runner is wired in.
- `mesh-002` Gmsh STEP import → tri count — workflow stub until a Gmsh-SDK runner exists.

---

## 7. Risk Register

| # | Risk | Likelihood | Impact | Mitigation |
|---|------|-----------|--------|------------|
| R1 | `cudarc` pre-alpha ships a breaking minor release mid-Phase 1 | Med | High | Pin `=0.19.x`; Agent C's internal `backend` trait makes the swap a single PR. |
| R2 | Gmsh SDK install fails or drifts version | Low | Med | Build is feature-gated; default CI doesn't exercise it; install documented in CONTRIBUTING. |
| R3 | `opencascade-rs` 5–15 min cold build | Med (Phase 1) | Med | Phase 0 does not enable the `opencascade` feature; sccache absorbs cost when Phase 1 begins. |
| R4 | `bindgen` requires `libclang` system dep | Med | Low | Pre-flight install in §5; Agent B brief notes it. |
| R5 | Per-worktree target dir disk pressure | Low | Low | sccache + `SCCACHE_CACHE_SIZE=10G`; clean up after merge. |
| R6 | Reviewer overwhelmed by 3-lane batches | Med | Med | Structured DoD checklist + P0/P1/P2 triage in template. |
| R7 | ROADMAP intro vs validation-list contradiction | Confirmed | Low | Agent G reconciles. |
| R8 | GPU runner never appears → `cuda-002` perma-skipped | High | Low (Phase 0) | Workflow stub + doc; Phase 1 makes it real. |
| R9 | Agent edits out of lane | Med | Med | Every brief states "out of lane → finding, do NOT fix"; reviewer enforces. |
| R10 | Scope creep into Phase 1 territory | Med | High | DoD is the gate; reviewer flags additions. |
| R11 | rustup version mismatch | Low | Low | `rust-toolchain.toml` auto-installs 1.85 on first cargo run. |
| R12 | `manylinux_2_28` wheel build (Phase 1) | n/a | n/a | Out of scope this phase. |
| R13 | `assert_cmd` integration tests race on shared target dir | Low | Low | Per-worktree target dir resolves. |
| R14 | `mdbook` install missing on Agent G host | Low | Low | Pre-flight install line. |
| R15 | Single-dev knowledge concentration | High | Med | Per-batch reviewer + per-crate ROADMAP + this spec = audit trail. |

---

## 8. Next Step

After this spec is approved, invoke the `superpowers:writing-plans` skill to produce a detailed implementation plan with per-agent task lists, dispatch order, merge checkpoints, and explicit reviewer hand-off prompts.
