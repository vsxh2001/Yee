# ADR-0016: yee-py exposes the validation aggregator with the slow real-run path env-gated

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

ADR-0008 established the `yee-validation` aggregator as the canonical
sink for Yee's published-benchmark validation cases: every case
(`mom-001`, `mom-002`, `mom-003`, `cpml_reflection`, `tfsf_slab`, …)
reports through `Report::run_all`, which emits a JSON record per case
plus optional PNG plots. The aggregator is a Rust crate
(`yee-validation`) consumed by `yee-cli validate` today.

Phase 1.frontend.4 wires the same aggregator through Python, so
notebook users can ingest validation results without re-implementing
the case-runner or shelling out to `yee` from `subprocess`. The thing
that makes this a design decision rather than a mechanical PyO3 export
is the **wall-time cost** of the real aggregator run.

The dominant case in the aggregator is `mom-001`, which CLAUDE.md §4
documents at **~7–8 minutes wall-time in release mode** on the
canonical 24×176 cylinder mesh. The whole aggregator, end-to-end,
takes **roughly 10 minutes**. The CI Python-bindings job
(`ci.yml :: python-bindings`) runs `maturin develop` followed by
`pytest`. If the pytest suite calls into a real aggregator run on
every push, the Python-bindings job alone would push CI past 15
minutes and we'd be in the business of explaining to contributors why
their typo-fix PR took a quarter of an hour to come back red.

Two structural responses present themselves:

1. **Expose only a no-op or mocked `run_validation()`.** Rejected.
   Defeats the point — notebook users want the *real* report, not
   a fixture. A mocked aggregator drifts from the real one within
   weeks and stops being useful at exactly the moment someone
   leans on it.
2. **Expose the real aggregator, gate the *test* that exercises it
   on an environment variable.** Accepted. See decision below.

The same pattern is used elsewhere in the workspace for slow /
hardware-gated paths: cuSOLVER tests are gated on `--features cuda`
and `-- --include-ignored` (CLAUDE.md §10); the GPU nightly workflow
is gated on `YEE_GPU_RUNNER_ENABLED` (ADR-0006, CLAUDE.md §8). The
common shape is "the *path* is in the public surface; the *invocation
in CI* is opt-in." This ADR applies that shape to the
validation-aggregator Python binding.

## Decision

`yee-py` exposes the validation aggregator under the **`yee`**
top-level package as a free function plus two thin wrapper classes:

```python
import yee

report: yee.ValidationReport = yee.run_validation()
for case in report.cases:
    case: yee.ValidationCase
    print(case.name, case.status, case.wall_time_s)
```

The exposed surface:

- **`yee.run_validation() -> ValidationReport`** — calls
  `yee_validation::Report::run_all` on the Rust side and returns a
  `ValidationReport` wrapper. **Runs serially on a single Python
  thread** (the calling thread, with the GIL released around the
  Rust work). No async, no thread pool.
- **`yee.ValidationReport`** — read-only wrapper exposing `.cases:
  list[ValidationCase]`, `.json: str` (the same JSON the Rust
  aggregator emits), and `.pass_count`, `.fail_count`,
  `.skipped_count` for quick filtering.
- **`yee.ValidationCase`** — read-only wrapper exposing `.name:
  str`, `.status: str` (one of `"Pass"`, `"Fail"`, `"Skipped"`),
  `.wall_time_s: float`, `.metric: dict[str, float] | None`, and
  `.png_path: str | None` (relative to the report's output
  directory; `None` if the case did not emit a plot).

The pytest suite for `yee-py` has **two tiers**:

- **Unconditional smoke tests.** Verify the binding *imports*, that
  `yee.run_validation` is a callable, that `ValidationReport` and
  `ValidationCase` expose the documented attributes, and that
  attribute types are correct. These run on every CI push. They do
  **not** call `yee.run_validation()`.
- **Real-aggregator integration test, gated on
  `YEE_RUN_VALIDATION=1`.** Decorated `@pytest.mark.skipif(
  os.environ.get("YEE_RUN_VALIDATION") != "1",
  reason="set YEE_RUN_VALIDATION=1 to run the ~10-min aggregator")`.
  When the env var is set, the test calls
  `yee.run_validation()` end-to-end and asserts the canonical pass
  set (`mom-001`, `cpml_reflection`, `tfsf_slab`, …) all report
  `Pass`. CI does not set the env var by default; the maintainer
  can flip it on a release-candidate branch.

The Rust side: `run_validation` releases the GIL via
`py.allow_threads(|| Report::run_all(...))` so the ~10 minute call
does not block other Python threads, and `ValidationReport` /
`ValidationCase` are `#[pyclass(frozen)]` immutable wrappers.

## Consequences

**What becomes easier:**

- **Notebooks can introspect Pass / Fail / Skipped cases** without
  re-implementing the aggregator or shelling out to `yee`. The
  Python wrapper hands back exactly the same JSON the Rust crate
  emits, plus typed accessors for the common queries (count by
  status, filter to fails, look up wall time per case).
- **The slow path is opt-in.** Nobody pays the ~10-minute cost
  unintentionally: CI default-skips the real run, and a developer
  pulling `yee-py` into a notebook can choose whether they want the
  fast import-and-shape smoke tests or the full real aggregator
  run.
- The opt-in mechanism (`YEE_RUN_VALIDATION=1`) is the same shape
  as the existing slow-path gates in the workspace (`--features
  cuda`, `YEE_GPU_RUNNER_ENABLED`), so contributors do not have to
  learn a new convention.

**What becomes harder:**

- **The aggregator runs serially on a single Python thread.** A
  user who wants async invocation, parallel case-runs, or a
  progress callback today has to wrap the call themselves
  (typically with `concurrent.futures.ThreadPoolExecutor`). That is
  acceptable for Phase 1.frontend.4 because the dominant cost is
  `mom-001` itself, not the aggregator orchestration — parallelism
  at the Python layer would not reduce wall time. **Phase
  1.frontend.5** is the place where an `async` / progress-callback
  API is on the queue if a user actually asks for it.
- A user can set `YEE_RUN_VALIDATION=1` locally and then be
  surprised by a 10-minute pytest run. The env-var name is
  deliberately verbose to make that surprise unlikely; the test's
  skip reason spells out the wall-time cost.

**What's now closed off:**

- Exposing a mocked or stubbed aggregator from `yee-py`. The
  Python surface is the real aggregator or nothing.
- Calling the aggregator from the unconditional pytest smoke
  tier. That tier is, by design, fast enough that CI can run it on
  every push without changing the cost profile of the
  `python-bindings` job.

## References

- `crates/yee-py/src/validation.rs` — `run_validation`,
  `ValidationReport`, `ValidationCase` pyclass wrappers; GIL
  released around `Report::run_all`.
- `crates/yee-py/tests/test_validation_smoke.py` — unconditional
  import + shape smoke tests.
- `crates/yee-py/tests/test_validation_real.py` — `skipif`-gated
  real-aggregator integration test.
- `.github/workflows/ci.yml :: python-bindings` — does **not**
  set `YEE_RUN_VALIDATION`; default-skips the slow test.
- ADR-0008 — validation aggregator JSON + PNG contract; this ADR
  is the Python wrapper over that surface.
- ADR-0006 — cudarc pre-alpha pin; same "feature-gated slow
  path" shape applied to GPU tests.
- CLAUDE.md §4 — `mom-001` wall-time (~7–8 minutes in release).
  CLAUDE.md §8 — CI shape; CLAUDE.md §10 — slow-path gating
  conventions.
- Phase 1.frontend.5 (queued, not yet specced) — async / progress
  callback API if and when a user requests it.
