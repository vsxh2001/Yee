# ADR-0007: Host criterion benchmarks in a separate `yee-bench` crate

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

Phase 1.bench landed three criterion-driven benchmarks: the MoM
matrix-fill / direct-LU solve at small N (`mom_solve`), an FDTD
single-step over a 50³ vacuum grid (`fdtd_step`), and an
iterative-vs-direct shootout on a 128×128 Hermitian-positive-definite
system (`gmres_vs_direct`). These are the first numbers the project
will quote when someone asks "is Yee fast?", and they are the first
artefacts that will let us detect a performance regression introduced
by a refactor.

The question is **where the benchmark code lives**. There are three
conventional patterns in the Rust ecosystem:

1. **Per-crate `benches/` directories.** Each solver crate
   (`yee-mom`, `yee-fdtd`, …) carries its own `benches/*.rs` and a
   `[dev-dependencies] criterion = ...` line. This is the Cargo
   default and is what `criterion`'s docs assume.
2. **A separate `yee-bench` crate** that depends on the solver
   crates as ordinary dependencies and hosts all benchmarks under
   `crates/yee-bench/benches/`.
3. **An external repo** holding nothing but benches plus a pinned
   commit of Yee. The `criterion`-comparison.dev model used by some
   compiler-perf efforts.

Option 1 is the path of least resistance but has real costs in a
multi-crate workspace where benchmarks reach across crates:

- A `gmres_vs_direct` bench naturally compares `yee-mom`'s iterative
  solver against the `faer`-backed direct LU, plus (in a follow-up)
  the cuSOLVER path through `yee-cuda`. Hosting it inside `yee-mom`'s
  `benches/` either forces a dev-dependency from `yee-mom` to
  `yee-cuda` (which inverts the layering) or scatters the comparison
  across two crates' bench harnesses.
- `criterion` itself is a non-trivial dev-dependency tree (`plotters`,
  `tinytemplate`, `serde_json`, `walkdir`, `clap`). Pulling it into
  every solver crate's `[dev-dependencies]` inflates every
  `cargo test --workspace` invocation, including CI's, even though
  benches don't run there.
- Per-crate `benches/` makes it harder to standardise the harness:
  warm-up duration, sample size, throughput definitions, output
  directory all want to live in one place.

Option 3 (external repo) is overkill for the current scale. It buys
real isolation but at the cost of a second repo to keep in lockstep
with `yee` itself, which is friction we don't yet need.

Option 2 is the conventional answer for workspaces that want
benchmarks as a first-class citizen without coupling them to the
default test path. The pattern shows up in `tokio` (the `tokio-perf`
crate), in `polars` (`polars-bench`), and historically in
`rustc-perf`.

The benches are **opt-in by default**:

- `cargo test --workspace` does not run them (benches are a separate
  Cargo target type and only build under `cargo bench`).
- The CI default-features matrix does not run `cargo bench`. A
  performance-regression workflow is a future-phase concern, not a
  per-PR gate.
- Local invocation is the explicit `cargo bench -p yee-bench` (or
  `cargo bench -p yee-bench --bench mom_solve` for one bench).

This deliberately accepts that performance regressions can land on
`main` without CI catching them. The trade is:

- **For:** Every PR stays fast. CI runtime is dominated by `mom-001`
  (7-8 min) already; adding `cargo bench` to per-PR CI would add
  another 10+ min of wall-time and produce noisy results on shared
  GitHub runners with no stable performance baseline.
- **Against:** Regressions are caught later, by maintainers running
  benches locally before a release or when investigating a complaint.

The verdict is that the per-release / on-demand cadence is the right
one for Phase 1. A follow-up `bench.yml` workflow that runs
`cargo bench -p yee-bench` on a self-hosted runner with stable CPU
pinning can be added later if regressions become a real problem.

## Decision

Host all criterion benchmarks in a single workspace member
**`crates/yee-bench`**:

```text
crates/yee-bench/
├── Cargo.toml
├── README.md
└── benches/
    ├── mom_solve.rs
    ├── fdtd_step.rs
    └── gmres_vs_direct.rs
```

`yee-bench`'s `Cargo.toml`:

```toml
[package]
name = "yee-bench"
version = "0.1.0"
edition = "2024"
license = "GPL-3.0-or-later"
publish = false

[dependencies]
yee-core = { path = "../yee-core" }
yee-mom  = { path = "../yee-mom" }
yee-fdtd = { path = "../yee-fdtd" }

[dev-dependencies]
criterion = { version = "0.5", features = ["html_reports"] }

[[bench]]
name = "mom_solve"
harness = false

[[bench]]
name = "fdtd_step"
harness = false

[[bench]]
name = "gmres_vs_direct"
harness = false
```

The crate is `publish = false`. It exists only to host benchmarks and
is not part of any release surface.

Benches are invoked explicitly:

```sh
cargo bench -p yee-bench                 # all benches
cargo bench -p yee-bench --bench mom_solve  # one bench
```

The default CI does not run `cargo bench`. A future
`.github/workflows/bench.yml` (not yet authored) is the intended home
for an automated baseline once a stable runner is available.

## Consequences

**What becomes easier:**

- Cross-crate comparisons (`gmres_vs_direct` against direct LU
  against, eventually, cuSOLVER) live in one file with one
  dependency graph. No per-crate `[dev-dependencies]` duplication.
- The criterion dependency tree is paid for once, not N times. A
  default `cargo test --workspace` does not compile criterion at all.
- The harness configuration (sample size, warm-up, plot output) lives
  in one place. Future additions inherit a consistent shape.
- `cargo bench -p yee-bench` is a single-command answer to "give me
  the numbers".

**What becomes harder:**

- Adding a benchmark for a method in `yee-mom` requires touching a
  different crate from the one being benchmarked. The path is two
  files away instead of one. The cost is small but real.
- The bench code does not have access to `pub(crate)` internals of
  the solver crates. Benchmarks must use the same public API the rest
  of the workspace sees, which is good discipline but occasionally
  forces a `pub(crate) fn` to become `pub fn` for the sole purpose of
  being benchmarked.

**What's now closed off:**

- Per-crate `benches/*.rs` directories. Solver crates do not carry
  their own bench harness. If a contributor adds a `benches/` folder
  under `yee-mom/`, the review feedback is to move it to `yee-bench`.
- Wiring `cargo bench` into the default CI matrix. The decision to do
  so is deferred to a separate workflow with a separate runner
  policy; it is not a per-PR gate.

## References

- `crates/yee-bench/` — the crate.
- `crates/yee-bench/benches/mom_solve.rs` — dipole 8×8 single-frequency
  matrix-fill / direct-LU.
- `crates/yee-bench/benches/fdtd_step.rs` — 50³ vacuum step.
- `crates/yee-bench/benches/gmres_vs_direct.rs` — 128×128 HPD shootout.
- `criterion` crate, <https://crates.io/crates/criterion>; HTML
  reports under `target/criterion/` after a run.
- `tokio-perf` and `polars-bench` as ecosystem precedents for the
  separate-bench-crate pattern.
- ADR-0006 — the `yee-cuda::backend` indirection that lets a future
  bench compare CPU and GPU paths through the same trait.
