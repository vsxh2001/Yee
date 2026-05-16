# Phase 0 Walking-Skeleton Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver the Yee Phase 0 walking skeleton — every workspace pipe connected, every build/test/clippy/doc gate green — via three batches of dispatched agents working in isolated git worktrees on a single host.

**Architecture:** Native `git worktree` per agent lane. Per-worktree `CARGO_TARGET_DIR`. Shared `sccache` compile cache. Three batches: A solo → B+C+D parallel → E+F+G parallel. One `feature-dev:code-reviewer` dispatch per batch, read-only, before merge.

**Tech Stack:** Rust 1.85 workspace; cudarc 0.19; Gmsh 4.15+ via bindgen; faer / nalgebra / ndarray; clap + clap_complete; mdBook; GitHub Actions CI.

**Companion spec:** `docs/superpowers/specs/2026-05-16-phase-0-multi-agent-execution-design.md`

---

## File Structure

The workspace is already scaffolded at commit `b6e1c00`. This plan modifies existing files and creates the additional files listed below.

| Crate / area | Files this plan creates | Files this plan modifies |
|---|---|---|
| yee-core | `crates/yee-core/tests/units.rs`, `crates/yee-core/tests/freq.rs` | `crates/yee-core/src/lib.rs` |
| yee-mesh | `crates/yee-mesh/build.rs`, `crates/yee-mesh/src/session.rs` | `crates/yee-mesh/src/lib.rs`, `crates/yee-mesh/Cargo.toml` |
| yee-cuda | `crates/yee-cuda/src/backend.rs`, `crates/yee-cuda/src/nvrtc.rs`, `crates/yee-cuda/kernels/hello.cu` | `crates/yee-cuda/src/lib.rs` |
| yee-io | `crates/yee-io/src/touchstone.rs`, `crates/yee-io/tests/touchstone_roundtrip.rs`, `crates/yee-io/validation/fixtures/touchstone/{1,2}port.s*p`, `crates/yee-io/validation/fixtures/touchstone/README.md` | `crates/yee-io/src/lib.rs` |
| yee-mom | `crates/yee-mom/tests/touchstone_roundtrip.rs` | `crates/yee-mom/src/lib.rs` |
| yee-cli | `crates/yee-cli/tests/cli.rs`, `crates/yee-cli/completions/{yee.bash,_yee,yee.fish}` | `crates/yee-cli/src/main.rs` |
| infra (G) | `.github/workflows/ci.yml`, `.github/workflows/gpu-nightly.yml`, `THIRD_PARTY_LICENSES.md`, `docs/book.toml`, `docs/src/SUMMARY.md`, `docs/src/introduction.md`, `.editorconfig`, `rustfmt.toml` | `ROADMAP.md` |

Worktrees live in `worktrees/<lane>` at the repo root and are cleaned up after merge.

---

## Conventions Used by Every Agent Brief

All agent briefs in this plan start with the shared preamble:

```bash
cd <worktree>
export CARGO_TARGET_DIR="$PWD/target"
export RUSTC_WRAPPER=sccache       # no-op if sccache absent
export SCCACHE_CACHE_SIZE=10G
git pull --ff-only origin main
```

Out-of-lane edits → surface as finding, do NOT fix.
If blocked > 15 min → surface specific blocker and stop, do NOT hack around.

---

## Task 0: Pre-flight installs

**Files:**
- No file changes; host-level installs only.

- [ ] **Step 1: Install Rust 1.85 toolchain**

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain 1.85
source "$HOME/.cargo/env"
rustup component add rustfmt clippy rust-src
```

Expected: `rustc --version` prints `rustc 1.85.x`.

- [ ] **Step 2: Install sccache (shared compile cache)**

```bash
cargo install sccache --locked
echo 'export RUSTC_WRAPPER=sccache'    >> ~/.bashrc
echo 'export SCCACHE_CACHE_SIZE=10G'  >> ~/.bashrc
source ~/.bashrc
sccache --start-server
sccache --show-stats
```

Expected: stats table with `Compile requests` initially 0.

- [ ] **Step 3: Install bindgen system dep**

```bash
sudo apt-get install -y libclang-dev pkg-config
```

Expected: `dpkg -l libclang-dev` shows it installed.

- [ ] **Step 4: Install mdBook**

```bash
cargo install mdbook --locked
mdbook --version
```

Expected: `mdbook v0.4.x` or newer.

- [ ] **Step 5: (Optional) Install Gmsh SDK for Agent B verify path**

```bash
mkdir -p /opt
cd /tmp
curl -L -O https://gmsh.info/bin/Linux/gmsh-4.15.0-Linux64-sdk.tgz
tar xzf gmsh-4.15.0-Linux64-sdk.tgz -C /opt/
mv /opt/gmsh-4.15.0-Linux64-sdk /opt/gmsh-sdk
echo 'export GMSH_SDK_ROOT=/opt/gmsh-sdk' >> ~/.bashrc
source ~/.bashrc
ls $GMSH_SDK_ROOT/include/gmshc.h
```

Expected: `gmshc.h` exists. (Skip this step if Agent B's `gmsh`-feature verify is not required; Phase 0 default CI does not exercise it.)

- [ ] **Step 6: (Optional) Verify CUDA toolkit for Agent C verify path**

```bash
nvcc --version
nvidia-smi -L
```

Expected: `nvcc` version 12.4+ or 13.0; `nvidia-smi -L` lists ≥ 1 GPU.

- [ ] **Step 7: Sanity check workspace builds on this host**

```bash
cd /home/hadassi/Code/Yee
cargo check --workspace --no-default-features 2>&1 | tail -20
```

Expected: exit 0, "Finished" line. (If this fails, fix before starting Batch 1.)

---

## Task 1: Batch 1 — Agent A (`yee-core`)

**Files:**
- Worktree: `worktrees/yee-core` on branch `agent/yee-core`
- Modify: `crates/yee-core/src/lib.rs`
- Create: `crates/yee-core/tests/units.rs`, `crates/yee-core/tests/freq.rs`

- [ ] **Step 1: Record base SHA and create worktree**

```bash
cd /home/hadassi/Code/Yee
BASE_SHA=$(git rev-parse main)
echo "$BASE_SHA"   # record this; expect b6e1c00 (or later after Task 0 sanity-check fixes)
git worktree add -b agent/yee-core worktrees/yee-core "$BASE_SHA"
```

Expected: `worktrees/yee-core` directory exists; `git worktree list` shows it.

- [ ] **Step 2: Dispatch Agent A**

Use the `Agent` tool with `subagent_type: "feature-dev:code-architect"` and the following brief:

````
You are implementing the `yee-core` crate, Phase 0 walking skeleton, for the Yee electromagnetic-simulation workspace.

PREAMBLE (run first):
```bash
cd /home/hadassi/Code/Yee/worktrees/yee-core
export CARGO_TARGET_DIR="$PWD/target"
export RUSTC_WRAPPER=sccache
export SCCACHE_CACHE_SIZE=10G
```

LANE: `crates/yee-core/**` only. Out-of-lane edits → surface as finding, do NOT fix.

ESCAPE HATCH: If blocked > 15 min, surface the specific blocker and stop.

DEFINITION OF DONE:
1. `units` module exports constants `C0`, `EPS0`, `MU0`, `ETA0` populated to CODATA 2018 values. A test file `tests/units.rs` verifies each constant within 1e-12 relative tolerance of the reference value (`ETA0 = sqrt(MU0/EPS0)` checked at the same tolerance).
2. `FreqRange::new(start_hz: f64, stop_hz: f64, n_points: usize) -> Result<Self>` validates `start_hz < stop_hz` and `n_points >= 1`; otherwise returns `Error::Invalid` with a useful message.
3. `FreqRange::iter()` yields `n_points` linearly spaced frequencies. For `n=1` yields `[start_hz]`. For `n=2` yields `[start_hz, stop_hz]` (endpoints exact). For `n>=3` yields evenly spaced including both endpoints. Test file `tests/freq.rs` asserts these cases.
4. `Error` enum carries the variants `Invalid(String)`, `Numerical(String)`, `Unimplemented(&'static str)`, `Io(String)`. `#[derive(thiserror::Error, Debug)]`.
5. Every public item carries `///` documentation; every public function has at least one doc-test that compiles and passes.
6. `cargo doc --no-deps -p yee-core` is warning-free.

WORK IN TDD ORDER for each requirement:
  a) Write the failing test (`cargo test -p yee-core` shows the new test red).
  b) Implement the minimal code.
  c) Re-run tests; expect green.
  d) Commit (one commit per requirement is fine; keep messages descriptive).

VERIFICATION COMMAND (must exit 0 before you report done):
```bash
cargo test -p yee-core \
  && cargo clippy -p yee-core --all-targets -- -D warnings \
  && cargo doc --no-deps -p yee-core
```

REPORT FORMAT:
- List the commits you made (oneline).
- Paste the final verification command output (last ~20 lines).
- Note any out-of-lane discoveries as findings.
````

Expected: agent reports DoD-met, verification command exits 0, no out-of-lane edits.

- [ ] **Step 3: Lane check**

```bash
cd /home/hadassi/Code/Yee
git -C worktrees/yee-core diff --stat main..HEAD
```

Expected: every modified path matches glob `crates/yee-core/**`. Anything outside → reject before merge.

- [ ] **Step 4: Dispatch reviewer for Batch 1**

Use the `Agent` tool with `subagent_type: "feature-dev:code-reviewer"` and the following brief:

````
You are reviewing Agent A's work on the `yee-core` crate. READ-ONLY. No edits.

Worktree: /home/hadassi/Code/Yee/worktrees/yee-core
Branch:   agent/yee-core

1. Lane check.
   ```bash
   cd /home/hadassi/Code/Yee/worktrees/yee-core
   git diff --stat main..HEAD
   ```
   Expected paths: `crates/yee-core/**` only.
   Any other path → P1 finding.

2. DoD checklist — verify each item by reading the diff and/or running the tests:
   - [ ] `units` constants exist and tests verify CODATA 2018 values ≤ 1e-12 relative
   - [ ] `FreqRange::new` validates inputs; rejects `start >= stop` and `n_points < 1` with `Error::Invalid`
   - [ ] `FreqRange::iter` exact endpoint behavior for n=1, n=2, n>=3
   - [ ] `Error` has `Invalid`, `Numerical`, `Unimplemented`, `Io` variants
   - [ ] Every public item documented; every public fn has a doc-test
   - [ ] `cargo doc --no-deps -p yee-core` warning-free

3. Re-run verification:
   ```bash
   cd /home/hadassi/Code/Yee/worktrees/yee-core
   export CARGO_TARGET_DIR="$PWD/target"
   cargo test -p yee-core \
     && cargo clippy -p yee-core --all-targets -- -D warnings \
     && cargo doc --no-deps -p yee-core
   ```
   Expected: exit 0 for the chain.

4. Code-quality smoke: any `unsafe`, panics in library code, TODO without owner, undocumented public items?

REPORT in this format:
- P0 (must fix before merge): <list, or "none">
- P1 (should fix this sprint): <list, or "none">
- P2 (deferrable): <list, or "none">
- Out-of-lane edits: <list with verdict — legit / hack / reject>
- Verification chain exit: 0 / non-zero (paste tail if non-zero)
````

Expected: P0 = none; out-of-lane = none; verification exit 0.

- [ ] **Step 5: User approval gate**

Present reviewer findings to the user. Do not proceed without explicit approval.

- [ ] **Step 6: Merge and cleanup**

```bash
cd /home/hadassi/Code/Yee
git merge --no-ff agent/yee-core -m "Merge Agent A: yee-core impl (Batch 1)"
git worktree remove worktrees/yee-core
git branch -d agent/yee-core
git log --oneline -3
```

Expected: merge commit on `main`; worktree gone; branch deleted.

---

## Task 2: Batch 2 setup — create three parallel worktrees

**Files:**
- Worktrees: `worktrees/yee-mesh`, `worktrees/yee-cuda`, `worktrees/yee-io`

- [ ] **Step 1: Record post-Batch-1 base SHA**

```bash
cd /home/hadassi/Code/Yee
BASE_SHA=$(git rev-parse main)
echo "$BASE_SHA"   # record; this is the base for Batch 2
```

- [ ] **Step 2: Create three worktrees**

```bash
git worktree add -b agent/yee-mesh worktrees/yee-mesh "$BASE_SHA"
git worktree add -b agent/yee-cuda worktrees/yee-cuda "$BASE_SHA"
git worktree add -b agent/yee-io   worktrees/yee-io   "$BASE_SHA"
git worktree list
```

Expected: three worktrees listed at the three paths, all on branches `agent/yee-mesh`, `agent/yee-cuda`, `agent/yee-io`.

---

## Task 3: Dispatch Batch 2 in parallel (B + C + D)

Send a single message containing three `Agent` tool calls so the agents run in parallel. Briefs follow.

**Files (B):**
- Worktree: `worktrees/yee-mesh`
- Create: `crates/yee-mesh/build.rs`, `crates/yee-mesh/src/session.rs`
- Modify: `crates/yee-mesh/src/lib.rs`, `crates/yee-mesh/Cargo.toml`

**Files (C):**
- Worktree: `worktrees/yee-cuda`
- Create: `crates/yee-cuda/src/backend.rs`, `crates/yee-cuda/src/nvrtc.rs`, `crates/yee-cuda/kernels/hello.cu`
- Modify: `crates/yee-cuda/src/lib.rs`

**Files (D):**
- Worktree: `worktrees/yee-io`
- Create: `crates/yee-io/src/touchstone.rs`, `crates/yee-io/tests/touchstone_roundtrip.rs`, `crates/yee-io/validation/fixtures/touchstone/{1,2}port.s*p`, `crates/yee-io/validation/fixtures/touchstone/README.md`
- Modify: `crates/yee-io/src/lib.rs`

- [ ] **Step 1: Dispatch Agent B (`yee-mesh`) in parallel block**

Brief:

````
You are implementing the `yee-mesh` crate, Phase 0 walking skeleton, for the Yee workspace.

PREAMBLE:
```bash
cd /home/hadassi/Code/Yee/worktrees/yee-mesh
export CARGO_TARGET_DIR="$PWD/target"
export RUSTC_WRAPPER=sccache
```

LANE: `crates/yee-mesh/**` only.

DEFINITION OF DONE:
1. `TriMesh::new(vertices, triangles, tags) -> Result<Self>` enforces the invariant `triangles.len() == tags.len()`. Returns `Error::Invalid` with a useful message otherwise. A test asserts both the success and rejection paths.
2. `build.rs` runs `bindgen` against `$GMSH_SDK_ROOT/include/gmshc.h` only when feature `gmsh` is enabled AND `$GMSH_SDK_ROOT` is set. Otherwise it writes an empty `bindings.rs` stub to `$OUT_DIR`.
3. New module `src/session.rs` exposes a safe skeleton:
   - `pub struct Session;`
   - `impl Session { pub fn new() -> Result<Self>; pub fn import_step(&mut self, path: &Path) -> Result<()>; pub fn mesh(&mut self, dim: i32) -> Result<()>; pub fn tris(&self) -> Result<TriMesh>; }`
   - `impl Drop for Session { ... }`
   - Without feature `gmsh`, every method returns `Error::NotEnabled`.
   - With feature `gmsh`, methods call into the generated bindings (you may stub the body with `todo!()` only when behind `#[cfg(feature = "gmsh")]` AND document the TODO in `ROADMAP.md` Phase 1 line — but PREFER a real call into gmshc init/finalize/open/generate/getElements, even if that's the minimum).
4. Crate builds clean with AND without `--features gmsh` (the `gmsh` build is verified only if `$GMSH_SDK_ROOT` is set).
5. Add `bindgen = { version = "0.71", optional = true }` to `[build-dependencies]` and gate it via the `gmsh` feature. (The scaffold has `optional = true` but no feature wiring; complete the wiring.)

WORK TDD-FIRST per requirement.

VERIFICATION COMMAND (default build):
```bash
cargo build -p yee-mesh \
  && cargo test -p yee-mesh \
  && cargo clippy -p yee-mesh -- -D warnings
```
Expected exit 0.

VERIFICATION COMMAND (if Gmsh SDK installed):
```bash
GMSH_SDK_ROOT=/opt/gmsh-sdk cargo build -p yee-mesh --features gmsh
```
Expected exit 0. If $GMSH_SDK_ROOT is not set, skip this command explicitly (do not silently pass).

REPORT FORMAT: commits, verification output (last 20 lines), out-of-lane findings.
````

- [ ] **Step 2: Dispatch Agent C (`yee-cuda`) in parallel block**

Brief:

````
You are implementing the `yee-cuda` crate, Phase 0 walking skeleton.

PREAMBLE:
```bash
cd /home/hadassi/Code/Yee/worktrees/yee-cuda
export CARGO_TARGET_DIR="$PWD/target"
export RUSTC_WRAPPER=sccache
```

LANE: `crates/yee-cuda/**` only.

DEFINITION OF DONE:
1. `Device` struct exposes fields `ordinal: usize`, `name: String`, `compute_cap: (u8, u8)`, `mem_total_bytes: u64`.
2. `Device::list() -> Result<Vec<Device>>` semantics:
   - Without feature `cuda`: returns `Err(Error::NotEnabled)`.
   - With feature `cuda` AND a visible GPU: returns one populated `Device` per visible GPU.
   - With feature `cuda` AND no visible GPU on a toolkit-installed host: returns `Ok(Vec::new())` (empty, no error).
3. Module `src/nvrtc.rs` exposes `pub fn compile(src: &str, name: &str) -> Result<Vec<u8>>` returning PTX bytes. Without feature `cuda`, returns `Err(Error::NotEnabled)`.
4. Module `src/backend.rs` defines a small `trait Backend` with associated types/methods covering the calls used (e.g. `fn device_count() -> Result<usize>`; `fn device_props(i: usize) -> Result<(String, (u8,u8), u64)>`; `fn nvrtc_compile(src, name) -> Result<Vec<u8>>`). Provide an implementation `CudarcBackend` that is `#[cfg(feature = "cuda")]`. Internal lib code routes through this trait so future swaps are local.
5. `kernels/hello.cu` is checked in with this content (exact):
   ```cuda
   extern "C" __global__ void add_one(float* out, const float* in, int n) {
       int i = blockIdx.x * blockDim.x + threadIdx.x;
       if (i < n) out[i] = in[i] + 1.0f;
   }
   ```
6. Crate builds clean with AND without `--features cuda`.
7. A simple unit test asserts `Device::list()` returns `Err(Error::NotEnabled)` when default-built.

VERIFICATION COMMAND (default build):
```bash
cargo build -p yee-cuda \
  && cargo test -p yee-cuda \
  && cargo clippy -p yee-cuda -- -D warnings
```

VERIFICATION COMMAND (if CUDA toolkit installed):
```bash
cargo build -p yee-cuda --features cuda
```
If toolkit absent, skip explicitly.

REPORT FORMAT: commits, verification output, findings.
````

- [ ] **Step 3: Dispatch Agent D (`yee-io`) in parallel block**

Brief:

````
You are implementing the `yee-io` crate Touchstone v1.1 reader/writer.

PREAMBLE:
```bash
cd /home/hadassi/Code/Yee/worktrees/yee-io
export CARGO_TARGET_DIR="$PWD/target"
export RUSTC_WRAPPER=sccache
```

LANE: `crates/yee-io/**` only.

REFERENCE: Touchstone v1.1 spec — https://ibis.org/connector/touchstone_spec11.pdf

DEFINITION OF DONE:
1. Module `src/touchstone.rs` (move the existing struct out of `lib.rs` into this module; re-export from `lib.rs`).
2. Parser handles option line `# <freq_unit> <param_type> <format> R <Z0>` where:
   - `freq_unit` ∈ {Hz, kHz, MHz, GHz} (case-insensitive)
   - `param_type` = S (this phase only)
   - `format` ∈ {RI, MA, DB} (case-insensitive)
   - `Z0` defaults to 50.0 if missing
3. Reader supports `.s1p`, `.s2p`, `.s3p`, `.s4p`. Decide port count from the file extension; verify against per-frequency column count and return `Error::TouchstoneParse{line, col, msg}` on mismatch.
4. Writer emits deterministic output:
   - Top: preserved comments (lines that started with `!` in the input).
   - Then the option line, formatted exactly: `# <unit> <S> <fmt> R <Z0>` with one space between fields, `Z0` printed with `%g`.
   - Then one line per frequency, columns space-separated.
5. `Error::TouchstoneParse { line: usize, col: usize, msg: String }`. Replace the placeholder string variant currently in `lib.rs`.
6. Property test in `tests/touchstone_roundtrip.rs`: for each fixture, `read → write → read` produces a struct equal to the first read, with floats compared at 1e-12 relative.
7. Passivity check: on read, compute `eig(S† S)` per frequency; if any eigenvalue > 1 + 1e-9, return `Error::TouchstoneParse` with a useful message. Test the failure path with a crafted non-passive fixture.
8. Fixture corpus (at least three files) checked into `validation/fixtures/touchstone/`:
   - `1port.s1p` — a published or synthetic example documented in the fixture README
   - `2port.s2p` — a published example
   - `2port_db.s2p` — same data in DB format
   Add `validation/fixtures/touchstone/README.md` listing provenance per fixture (URL or "synthetic — generated from analytical model X").

WORK TDD-FIRST.

VERIFICATION COMMAND:
```bash
cargo test -p yee-io \
  && cargo clippy -p yee-io -- -D warnings
```
Expected exit 0.

REPORT FORMAT: commits, verification output, findings.
````

- [ ] **Step 4: Wait for all three to return**

When the parallel block returns, each agent should have committed its work to its branch.

- [ ] **Step 5: Lane checks for B, C, D**

```bash
cd /home/hadassi/Code/Yee
for lane in yee-mesh yee-cuda yee-io; do
  echo "=== $lane ==="
  git -C worktrees/$lane diff --stat main..HEAD
done
```

Expected per lane: only paths under `crates/<lane>/**`. Anything else → reject.

---

## Task 4: Batch 2 reviewer

**Files:**
- No file changes; reviewer is read-only.

- [ ] **Step 1: Dispatch Batch-2 reviewer**

Use `Agent` with `subagent_type: "feature-dev:code-reviewer"` and the following brief:

````
You are reviewing Batch 2 of the Yee Phase 0 walking skeleton. Three lanes: `yee-mesh`, `yee-cuda`, `yee-io`. READ-ONLY. No edits.

Worktrees:
- /home/hadassi/Code/Yee/worktrees/yee-mesh  (branch agent/yee-mesh)
- /home/hadassi/Code/Yee/worktrees/yee-cuda  (branch agent/yee-cuda)
- /home/hadassi/Code/Yee/worktrees/yee-io    (branch agent/yee-io)

FOR EACH LANE:

1. Lane check:
   ```bash
   cd /home/hadassi/Code/Yee/worktrees/<lane>
   git diff --stat main..HEAD
   ```
   Expected paths: `crates/<lane>/**` only. Anything else → P1 finding.

2. DoD checklist — verify each item from the agent brief (copy the DoD items from the original brief in the implementation plan):
   - yee-mesh DoD: TriMesh invariant + build.rs gating + Session skeleton + dual-feature build + Cargo.toml bindgen wiring.
   - yee-cuda DoD: Device fields + list() semantics + nvrtc::compile + backend trait + hello.cu content exact + dual-feature build.
   - yee-io DoD: option-line parsing + .s1p–.s4p + deterministic writer + TouchstoneParse{line,col,msg} + property test + passivity check + 3 fixtures with provenance.

3. Re-run verification per lane:
   ```bash
   cd /home/hadassi/Code/Yee/worktrees/<lane>
   export CARGO_TARGET_DIR="$PWD/target"
   cargo test -p <lane> \
     && cargo clippy -p <lane> --all-targets -- -D warnings
   ```
   Expected: exit 0 chain per lane.

4. Quality smoke: `unsafe`, panics in lib code, undocumented public items, TODO without owner, hard-coded paths beyond the SDK env var.

REPORT FORMAT:
- Per lane: P0 / P1 / P2 / Out-of-lane / Verification exit.
- Cross-lane: any unexpected public-API drift, name collisions, or coupling across lanes (should be NONE — Batch 2 lanes are independent).
````

- [ ] **Step 2: Present findings to user**

Wait for explicit approval before merging. Do not auto-merge even if all P0 = none.

---

## Task 5: Merge Batch 2

**Files:**
- Modifies `main` via three merge commits.

- [ ] **Step 1: Merge each lane in turn**

```bash
cd /home/hadassi/Code/Yee
git merge --no-ff agent/yee-mesh -m "Merge Agent B: yee-mesh impl (Batch 2)"
git merge --no-ff agent/yee-cuda -m "Merge Agent C: yee-cuda impl (Batch 2)"
git merge --no-ff agent/yee-io   -m "Merge Agent D: yee-io impl (Batch 2)"
git log --oneline -10
```

Expected: three merge commits land on `main` after Task 1's merge.

- [ ] **Step 2: Cleanup worktrees**

```bash
git worktree remove worktrees/yee-mesh
git worktree remove worktrees/yee-cuda
git worktree remove worktrees/yee-io
git branch -d agent/yee-mesh agent/yee-cuda agent/yee-io
git worktree list
```

Expected: only the main worktree listed.

- [ ] **Step 3: Workspace sanity check on `main`**

```bash
cargo check --workspace --no-default-features 2>&1 | tail -10
cargo test --workspace --no-default-features 2>&1 | tail -10
```

Expected: exit 0 for both. (If a merge introduced a conflict the agents could not see, fix here as orchestrator and commit before Batch 3 starts.)

---

## Task 6: Batch 3 setup — three more worktrees

**Files:**
- Worktrees: `worktrees/yee-mom`, `worktrees/yee-cli`, `worktrees/yee-infra`

- [ ] **Step 1: Record post-Batch-2 base SHA**

```bash
cd /home/hadassi/Code/Yee
BASE_SHA=$(git rev-parse main)
echo "$BASE_SHA"
```

- [ ] **Step 2: Create three worktrees**

```bash
git worktree add -b agent/yee-mom   worktrees/yee-mom   "$BASE_SHA"
git worktree add -b agent/yee-cli   worktrees/yee-cli   "$BASE_SHA"
git worktree add -b agent/yee-infra worktrees/yee-infra "$BASE_SHA"
git worktree list
```

Expected: three new worktrees listed.

---

## Task 7: Dispatch Batch 3 in parallel (E + F + G)

Send a single message with three `Agent` tool calls.

**Files (E):**
- Worktree: `worktrees/yee-mom`
- Modify: `crates/yee-mom/src/lib.rs`
- Create: `crates/yee-mom/tests/touchstone_roundtrip.rs`

**Files (F):**
- Worktree: `worktrees/yee-cli`
- Modify: `crates/yee-cli/src/main.rs`
- Create: `crates/yee-cli/tests/cli.rs`, `crates/yee-cli/completions/{yee.bash,_yee,yee.fish}`

**Files (G):**
- Worktree: `worktrees/yee-infra`
- Create: `.github/workflows/ci.yml`, `.github/workflows/gpu-nightly.yml`, `THIRD_PARTY_LICENSES.md`, `docs/book.toml`, `docs/src/SUMMARY.md`, `docs/src/introduction.md`, `.editorconfig`, `rustfmt.toml`
- Modify: `ROADMAP.md`

- [ ] **Step 1: Dispatch Agent E (`yee-mom` skeleton wiring)**

Brief:

````
You are implementing the `yee-mom` Phase 0 skeleton wiring.

PREAMBLE:
```bash
cd /home/hadassi/Code/Yee/worktrees/yee-mom
export CARGO_TARGET_DIR="$PWD/target"
export RUSTC_WRAPPER=sccache
```

LANE: `crates/yee-mom/**` only.

DEFINITION OF DONE:
1. `impl Solver for PlanarMoM` is real (not just compile-stub): `fn run(...) -> yee_core::Result<SParameters>` returns `Err(yee_core::Error::Unimplemented("PlanarMoM::run not implemented in phase 0"))`. A unit test asserts this exact behavior.
2. `SParameters` carries `freq_hz: Vec<f64>`, `data: Vec<Vec<Complex64>>`, `n_ports: usize`. Add:
   - `pub fn from_touchstone(file: yee_io::touchstone::File) -> Self`
   - `pub fn to_touchstone(&self, z0: f64) -> yee_io::touchstone::File`
   - `pub fn write_touchstone(&self, path: &Path, z0: f64) -> Result<()>` which goes through `yee_io::touchstone::write`
3. `tests/touchstone_roundtrip.rs` builds a small `SParameters` (n_ports=2, three frequencies, arbitrary complex data), writes to a tempfile, reads it back via `yee_io::touchstone::read`, asserts struct equality to 1e-12 relative.
4. `PlanarMoM::default()` constructs without panic.

WORK TDD-FIRST.

VERIFICATION COMMAND:
```bash
cargo test -p yee-mom \
  && cargo clippy -p yee-mom -- -D warnings
```

REPORT: commits, verification output, findings.
````

- [ ] **Step 2: Dispatch Agent F (`yee-cli` real wiring)**

Brief:

````
You are implementing the `yee-cli` Phase 0 real wiring.

PREAMBLE:
```bash
cd /home/hadassi/Code/Yee/worktrees/yee-cli
export CARGO_TARGET_DIR="$PWD/target"
export RUSTC_WRAPPER=sccache
```

LANE: `crates/yee-cli/**` only.

DEFINITION OF DONE:
1. `yee validate <mom|fdtd|all>` dispatches:
   - `mom` calls into a new function `yee_mom::validation::report() -> String` that returns a Phase-0 stub report listing planned cases (you add this fn to `yee-mom` in this lane — it is a `yee-mom` public addition that is acceptable here because it is consumed by the CLI; if the reviewer flags this as cross-lane, accept moving it to `yee-mom` in a follow-up commit before merge).
   - `fdtd` prints `Phase 2 deliverable — yee-fdtd not yet available`.
   - `all` runs both.
   - Exit 0 on success.
2. `yee mesh <path>` constructs a `yee_mesh::Session`. Without `--features gmsh` it prints `mesh feature not enabled; rebuild with --features gmsh` and exits with code 2.
3. `yee export <input> --format touchstone <output>` reads via `yee_io::touchstone::read` and writes back to `<output>` (round-trip). `--format hdf5` prints `hdf5 not yet enabled` and exits 2.
4. `clap_complete` generates shell completions on a hidden subcommand `yee completions <shell>`. Pre-generated completions for bash, zsh, fish are committed to `crates/yee-cli/completions/` (regenerated via that subcommand).
5. Integration tests via `assert_cmd` in `tests/cli.rs`:
   - `yee --help` exits 0 and stdout contains the subcommand names `validate`, `mesh`, `run`, `export`.
   - `yee --version` exits 0 and matches the workspace `version`.
   - `yee validate all` exits 0 and stdout contains "Phase 2 deliverable".
   - `yee garbage-subcmd` exits non-zero with stderr containing a suggestion.

WORK TDD-FIRST.

VERIFICATION COMMAND:
```bash
cargo test -p yee-cli \
  && cargo run --bin yee -- --help \
  && cargo run --bin yee -- validate all
```

REPORT: commits, verification output, findings.

NOTE on `yee-mom::validation::report`: if you find it cleaner to keep this in `yee-cli` as a local helper that calls into `yee-mom` symbols already public, do so; the goal is exit-0 user-visible behavior, not cross-crate purity in Phase 0. Surface the choice in your report.
````

- [ ] **Step 3: Dispatch Agent G (cross-cutting infrastructure)**

Brief:

````
You are implementing cross-cutting infrastructure for Yee Phase 0.

PREAMBLE:
```bash
cd /home/hadassi/Code/Yee/worktrees/yee-infra
export CARGO_TARGET_DIR="$PWD/target"
export RUSTC_WRAPPER=sccache
```

LANE: `.github/**`, `THIRD_PARTY_LICENSES.md`, `docs/**` (EXCLUDING `docs/source/**` and `docs/superpowers/**`), top-level `ROADMAP.md` ONLY for the Phase 0 validation-list reconciliation, `.editorconfig`, `rustfmt.toml`.

DEFINITION OF DONE:

1. `.github/workflows/ci.yml`:
   - Trigger: pull_request to main, push to main.
   - Job `lint-test` on `ubuntu-latest`, Rust 1.85 (use `dtolnay/rust-toolchain@stable` with `toolchain: "1.85"`).
   - Steps: checkout; rustup components rustfmt, clippy; `cargo fmt --check --all`; `cargo clippy --workspace --all-targets -- -D warnings`; `cargo test --workspace --no-default-features`; `cargo doc --workspace --no-deps`.
   - Cache cargo registry and target via `Swatinem/rust-cache@v2`.

2. `.github/workflows/gpu-nightly.yml`:
   - Cron schedule `0 4 * * *` UTC.
   - Single job with `runs-on: self-hosted-gpu-placeholder` and `if: false`.
   - Inline comments explain how to enable when a self-hosted GPU runner is registered (set `if: true` and update `runs-on:`).

3. `THIRD_PARTY_LICENSES.md`:
   - Sections: Gmsh (GPL v2+ with linking exception, FAQ link), OpenCASCADE Technology (LGPL 2.1 with exception, FAQ link), NVIDIA CUDA Libraries (proprietary, dynamic link, list cuBLAS / cuSOLVER / cuSPARSE / cuFFT / NCCL), Permissive Rust Dependencies (one-line table: name, license, link).
   - Cite where each appears in `Cargo.toml`.

4. mdBook scaffold:
   - `docs/book.toml` with `title = "Yee"`, `authors = ["The Yee Authors"]`, `src = "src"`.
   - `docs/src/SUMMARY.md` minimum:
     ```
     # Summary

     - [Introduction](introduction.md)
     ```
   - `docs/src/introduction.md` — one-page overview pulling from the top-level README (paraphrase, don't duplicate verbatim).
   - `mdbook build docs/` exits 0.

5. Top-level `ROADMAP.md` reconciliation:
   - Open `ROADMAP.md`.
   - Locate the Phase 0 "Validation milestones" section.
   - Move the three bullet points (half-wave dipole, microstrip line, patch antenna) out of Phase 0 into a footnote that says "Originally listed under Phase 0; reclassified as Phase 1 in the 2026-05-16 Phase 0 walking-skeleton design (`docs/superpowers/specs/2026-05-16-phase-0-multi-agent-execution-design.md`)."
   - Replace the Phase 0 Validation milestones list with the ten gates from §1 of the design spec.

6. `.editorconfig`:
   ```
   root = true

   [*]
   charset = utf-8
   end_of_line = lf
   indent_style = space
   indent_size = 4
   insert_final_newline = true
   trim_trailing_whitespace = true

   [*.{md,yml,yaml,toml}]
   indent_size = 2
   ```

7. `rustfmt.toml`:
   ```
   edition = "2024"
   max_width = 100
   use_field_init_shorthand = true
   ```

VERIFICATION COMMAND:
```bash
cargo fmt --check --all \
  && mdbook build docs/ \
  && grep -q "Gmsh" THIRD_PARTY_LICENSES.md \
  && grep -q "OpenCASCADE" THIRD_PARTY_LICENSES.md \
  && grep -q "CUDA" THIRD_PARTY_LICENSES.md
```
Expected: exit 0.

REPORT: commits, verification output, findings.
````

- [ ] **Step 4: Wait for all three to return**

- [ ] **Step 5: Lane checks for E, F, G**

```bash
cd /home/hadassi/Code/Yee
for lane in yee-mom yee-cli yee-infra; do
  echo "=== $lane ==="
  git -C worktrees/$lane diff --stat main..HEAD
done
```

Expected per lane:
- `yee-mom` → `crates/yee-mom/**`
- `yee-cli` → `crates/yee-cli/**` (note: Agent F brief explicitly allows touching `yee-mom` for `validation::report`; if the agent took that path, treat as legit cross-lane and document in merge commit)
- `yee-infra` → `.github/**`, `THIRD_PARTY_LICENSES.md`, `docs/**` (NOT `docs/source/**` or `docs/superpowers/**`), `ROADMAP.md`, `.editorconfig`, `rustfmt.toml`

---

## Task 8: Batch 3 reviewer

- [ ] **Step 1: Dispatch Batch-3 reviewer**

Use `Agent` with `subagent_type: "feature-dev:code-reviewer"` and this brief:

````
You are reviewing Batch 3 of the Yee Phase 0 walking skeleton. Three lanes: `yee-mom`, `yee-cli`, `yee-infra`. READ-ONLY.

Worktrees:
- /home/hadassi/Code/Yee/worktrees/yee-mom    (agent/yee-mom)
- /home/hadassi/Code/Yee/worktrees/yee-cli    (agent/yee-cli)
- /home/hadassi/Code/Yee/worktrees/yee-infra  (agent/yee-infra)

PER LANE:

1. Lane check via `git diff --stat main..HEAD` from the worktree.
   - yee-mom: only `crates/yee-mom/**`.
   - yee-cli: `crates/yee-cli/**` + possibly small `crates/yee-mom/src/lib.rs` addition for `validation::report` (allowed per brief).
   - yee-infra: `.github/**`, `THIRD_PARTY_LICENSES.md`, `docs/**` excluding `docs/source/**` and `docs/superpowers/**`, `ROADMAP.md`, `.editorconfig`, `rustfmt.toml`.

2. DoD checklist:
   - yee-mom: Solver impl returns Unimplemented; SParameters has from/to/write_touchstone; round-trip test passes.
   - yee-cli: subcommands dispatch correctly; `mesh` without `gmsh` exits 2 with the right message; `export touchstone` round-trips; completions for bash/zsh/fish exist; assert_cmd tests pass.
   - yee-infra: ci.yml runs fmt+clippy+test+doc; gpu-nightly.yml stubbed with `if: false`; THIRD_PARTY_LICENSES.md has Gmsh + OCCT + CUDA sections; mdBook builds; ROADMAP reconciled per spec §6 deferral; .editorconfig and rustfmt.toml content matches brief.

3. Re-run verification per lane:
   ```bash
   # yee-mom
   cd /home/hadassi/Code/Yee/worktrees/yee-mom
   export CARGO_TARGET_DIR="$PWD/target"
   cargo test -p yee-mom && cargo clippy -p yee-mom -- -D warnings

   # yee-cli
   cd /home/hadassi/Code/Yee/worktrees/yee-cli
   export CARGO_TARGET_DIR="$PWD/target"
   cargo test -p yee-cli && cargo run --bin yee -- --help && cargo run --bin yee -- validate all

   # yee-infra
   cd /home/hadassi/Code/Yee/worktrees/yee-infra
   export CARGO_TARGET_DIR="$PWD/target"
   cargo fmt --check --all && mdbook build docs/
   ```
   Expected exit 0 per chain.

4. Quality smoke: same as Batch 2.

REPORT: per-lane P0/P1/P2/Out-of-lane/Verification exit; cross-lane interactions noted.
````

- [ ] **Step 2: Present findings; user approval gate**

---

## Task 9: Merge Batch 3 and verify Phase 0 done

**Files:**
- Three merge commits on `main`.

- [ ] **Step 1: Merge each lane**

```bash
cd /home/hadassi/Code/Yee
git merge --no-ff agent/yee-mom   -m "Merge Agent E: yee-mom skeleton wiring (Batch 3)"
git merge --no-ff agent/yee-cli   -m "Merge Agent F: yee-cli real wiring (Batch 3)"
git merge --no-ff agent/yee-infra -m "Merge Agent G: cross-cutting infra (Batch 3)"
git log --oneline -10
```

Expected: three merge commits land.

- [ ] **Step 2: Cleanup worktrees**

```bash
git worktree remove worktrees/yee-mom
git worktree remove worktrees/yee-cli
git worktree remove worktrees/yee-infra
git branch -d agent/yee-mom agent/yee-cli agent/yee-infra
git worktree list
```

Expected: only main worktree listed.

- [ ] **Step 3: Run all ten Phase 0 done gates from `main`**

```bash
cd /home/hadassi/Code/Yee
cargo check --workspace --no-default-features                              # v0-build
cargo test --workspace --no-default-features                               # v0-test
cargo clippy --workspace --all-targets -- -D warnings                      # v0-clippy
cargo fmt --check --all                                                    # v0-fmt
cargo doc --workspace --no-deps                                            # v0-doc
cargo run --bin yee -- --help                                              # v0-cli
cargo run --bin yee -- validate all                                        # v0-cli-validate
mdbook build docs/                                                         # v0-mdbook
grep -q "Gmsh" THIRD_PARTY_LICENSES.md && \
  grep -q "OpenCASCADE" THIRD_PARTY_LICENSES.md && \
  grep -q "CUDA" THIRD_PARTY_LICENSES.md && echo "licenses OK"             # v0-licenses
# v0-ci is verified on the next push by GitHub Actions
```

Expected: every command exits 0. `v0-licenses` step prints `licenses OK`.

- [ ] **Step 4: Push to remote to trigger CI (v0-ci gate)**

```bash
cd /home/hadassi/Code/Yee
git push origin main
gh run watch --exit-status   # if gh CLI configured; otherwise check the Actions tab
```

Expected: CI run completes green. (If no remote yet, skip this step and flag for the user.)

- [ ] **Step 5: Tag Phase 0 completion**

```bash
cd /home/hadassi/Code/Yee
git tag -a phase-0-done -m "Phase 0 walking skeleton complete — all ten gates green"
git push origin phase-0-done   # if remote configured
git log --oneline -15
```

Expected: tag `phase-0-done` exists at current HEAD.

- [ ] **Step 6: Final commit / no-op confirmation**

If any orchestrator-side fixes were applied during Tasks 5 or 9 (merge-conflict resolutions, etc.), confirm they are already on `main`. No further commit needed.

---

## Self-Review (run after writing this plan; orchestrator reads this once)

1. **Spec coverage:** spec §1 (Phase 0 done) covered by Task 9 Step 3 (10 gates). Spec §2 (locked decisions) covered structurally by Task 0 + Tasks 1–9. Spec §3 (architecture) was scaffolded pre-plan; this plan modifies/creates the listed files. Spec §4 (per-agent briefs) covered by Tasks 1, 3, 7 dispatches; DoDs reproduced verbatim in brief blocks. Spec §5 (orchestration) covered by worktree create/merge/remove flow in every task. Spec §6 (validation gates) covered by Task 9 Step 3. Spec §7 (risks) — mitigations are embedded in briefs (pin cudarc minor, feature-gate gmsh/cuda, lane enforcement, escape hatch).

2. **Placeholder scan:** no "TBD" / "TODO without owner" / "implement later" / "add error handling" outside of Agent B's allowed `todo!()` inside `#[cfg(feature = "gmsh")]` Phase-1-deferred path, which is documented and scoped.

3. **Type consistency:** `SParameters` field set (`freq_hz`, `data`, `n_ports`) matches scaffold `crates/yee-mom/src/lib.rs` and Agent D's `touchstone::File` shape (`z0`, `n_ports`, `freq_hz`, `data`). `Error::TouchstoneParse` is upgraded from `String` to struct variant in Agent D's brief — Agent E and Agent F's briefs do not depend on its internal shape, only on its existence.

4. **Cross-lane writes called out:** Agent F's brief explicitly flags the `yee-mom::validation::report` cross-lane addition and tells the reviewer to expect it.

5. **No assumed external services:** GPU runner is stubbed `if: false`; Gmsh and CUDA verify paths are explicitly skip-able.

No issues found that require revising the plan.
