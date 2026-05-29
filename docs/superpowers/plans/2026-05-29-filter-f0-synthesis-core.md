# Filter Phase F0 — synthesis core — Implementation Plan

**Spec:** `docs/superpowers/specs/2026-05-29-filter-f0-synthesis-core-design.md`
**ADR:** ADR-0084
**Phase:** F0 (filter roadmap)
**Date:** 2026-05-29

---

## Lane

`crates/yee-synth/**` (new), `crates/yee-filter/**` (new), root `Cargo.toml`
(add 2 members), `crates/yee-cli/**` (new `Filter` subcommand). Out of lane
(any solver crate, `yee-validation`, GUI, docs already committed) → finding,
not fix. **`yee-validation` registration is Phase F0.1 — do NOT touch it here.**

## Base

Worktree `worktrees/filter-f0`, branch `feature/filter-f0-synthesis-core`,
base `origin/main` `3101d9d`.

## Pattern files (imitate house style)

- A small existing crate for `Cargo.toml` + `lib.rs` shape + `#![warn(missing_docs)]`
  + `#![forbid(unsafe_code)]`: `crates/yee-io/` or `crates/yee-design/`.
- Existing CLI subcommand wiring + tests: `crates/yee-cli/src/main.rs`
  (`Command` enum, `run`) and `crates/yee-cli/tests/cli_validate.rs`.
- Touchstone writing: `crates/yee-io/` public API (use it; do not re-implement).

## Steps

1. **`crates/yee-synth`** (new lib). `Cargo.toml`: deps `yee-core`, `nalgebra`
   (workspace versions), `serde` (derive) for `Approximation`. `src/lib.rs`
   with `#![forbid(unsafe_code)]`, `#![warn(missing_docs)]`. Implement spec §2:
   `Approximation`, `Prototype`, `prototype()`, `min_order()`, `CouplingDesign`,
   `coupling_design()`. Use the exact formulas in spec §2.1–2.5.
2. **`crates/yee-filter`** (new lib). Deps `yee-core`, `yee-synth`, `yee-io`,
   `serde`, `num-complex`. Implement spec §3: `FilterSpec`, `SpecMask`,
   `Response`, `CouplingMatrix`, `Topology`, `FilterProject`, `synthesize()`,
   `ideal_response()` (closed-form Chebyshev/Butterworth transfer fn on the
   bandpass-mapped Ω), `check_mask()` + `MaskReport`. Re-export
   `yee_synth::Approximation`.
3. **Root `Cargo.toml`**: add `"crates/yee-synth"` and `"crates/yee-filter"` to
   `members` (alongside the existing crates).
4. **`crates/yee-cli`**: add a `Filter` subcommand with a `Synth { spec: PathBuf,
   output: Option<PathBuf> }` variant (clap). Parse the `FilterSpec` TOML,
   `yee_filter::synthesize`, print prototype g-values + coupling matrix + Qe,
   sweep `ideal_response`, write Touchstone via `yee-io`, print `check_mask`
   verdict; `ExitCode::SUCCESS` on PASS, `FAILURE` on mask FAIL. Add a fixture
   spec (e.g. `crates/yee-cli/tests/fixtures/cheb_bpf.toml`).
5. **Tests = the gates** (spec §DoD 4–6), as `#[test]`s:
   - `crates/yee-synth/tests/synth_001_gvalues.rs` — assert Butterworth N=3/N=5
     and Chebyshev 0.5 dB N=3/N=4/N=5 and 3.0 dB N=3 g-values vs the spec's
     published numbers to ≤1e-3 (numbers are in spec §DoD 4 — use them verbatim).
   - `crates/yee-synth/tests/synth_002_coupling.rs` — Chebyshev 0.5 dB N=3,
     FBW=0.10: assert `k_12==k_23`, and `k`, `Qe_in==Qe_out` match the §2.5
     closed form recomputed in-test to ≤1e-9 (self-consistency) AND the
     hand-computed literal to ≤1e-3.
   - `crates/yee-filter/tests/filt_001_mask.rs` — synthesize a Chebyshev 0.5 dB
     bandpass at a satisfiable order; `check_mask` over a swept band → PASS
     (ripple ≤0.5 dB, RL ≥ spec, rejection ≥ mask). A too-low order → FAIL
     (negative control).
6. **CLI smoke test** `crates/yee-cli/tests/cli_filter.rs` — `yee filter synth
   <fixture>` exits 0, stdout contains the coupling matrix + "PASS", and a
   Touchstone file is written.

## Verify (all exit 0; nice -n 19, --jobs 2; NO --workspace, NO EM)

```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-synth -p yee-filter -p yee-cli --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-synth -p yee-filter -p yee-cli --jobs 2
nice -n 19 cargo run -p yee-cli --jobs 2 -- filter synth crates/yee-cli/tests/fixtures/cheb_bpf.toml
```
Do NOT run `cargo test --workspace` (pulls the ~8 min mom-001 + ~31 min
fem-eig-003). These crates are pure-math and build/test in seconds.

## Math correctness notes (get these right)

- Chebyshev even-order load is NOT 1.0: `g_{N+1}=coth²(β/4)`. Test N=4 covers it.
- `coth(x)=cosh/sinh`; `β=ln(coth(L_Ar/17.37))`. The `17.37 = 40/ln(10)` constant
  is the standard Pozar form — keep it.
- Butterworth `g_{N+1}=1` always; Chebyshev odd-order `g_{N+1}=1`.
- Bandpass map: `Ω = (1/FBW)·(ω/ω0 − ω0/ω)` then feed `Ω` to the lowpass `|S21|²`.

## Escape hatch

Blocked > 15 min — e.g. g-values miss the published table beyond 1e-3 (suspect
the β/γ constants or the even-order load), or a workspace-member/dep cycle —
STOP and surface the exact mismatch (print computed vs expected). Do NOT loosen
a gate tolerance to pass; the published numbers are ground truth.

## Done when

DoD 1–7 pass; lane respected (`git diff --stat 3101d9d..HEAD` shows only the new
crates + root Cargo.toml + yee-cli + the 4 docs already committed); `cargo build`
of the two new crates is warning-free.
