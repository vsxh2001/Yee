# FS.1c Task 3 — ADR-0228 + roadmap row

**Branch:** `feature/fs1c-thin-wire`
**head_before:** `fcb73a7b83fbd735817ae1d94bc1adb3acf20a1e` (docs: FS.1c Task 2 report)
**head_after:** `028853b24f7c0582b6c1c45b65cac3aa68a44072` — `docs: ADR-0228 — FS.1c
thin-wire (Holland-Simpson) COMPLETE`

Plan: `docs/superpowers/plans/2026-07-24-fs1c-thin-wire.md`, Task 3.
Spec: `docs/superpowers/specs/2026-07-24-fs1c-thin-wire-design.md`.
Task 1 report: `.superpowers/sdd/fs1c-task-1-report.md`.
Task 2 report: `.superpowers/sdd/fs1c-task-2-report.md`.

## What shipped (docs-only, per the Task 3 lane exception)

1. **`docs/src/decisions/0228-fs1c-thin-wire.md`** — new ADR, structured per
   ADR-0227's shape (Context / Decision (per-task subsections) / Measured
   result / Tolerances pinned / Bit-exactness discipline / Verdict / What
   remains). Contents transcribed from, not reinterpreted beyond, the Task 1
   and Task 2 reports and the shipped module docs / gate module doc:
   - the exact Holland & Simpson (1981) citation chain (+ the Liu 2003 thesis
     derivation, Umashankar/Taflove/Beker 1987, Taflove & Hagness ch. 10) —
     copied verbatim from `drive.rs`'s `ThinWire` doc, not re-derived;
   - the seam decision (`Drive`, not `Materials`) and its rationale;
   - the two translation bugs caught by the coarse/fine gate itself (the
     missing-`dz` NaN blowup, the missing open-end `I=0` drift);
   - Task 1's coarse/fine consistency number (~8.1%, pinned `< 10%`);
   - Task 2's measured `engine-thinwire-dipole-001` numbers verbatim
     (`Z(149.8962 MHz) = 92.045 + j109.464 Ω` vs NEC-4 87+j41 Ω; Re err 5.8%
     vs 10% tol; Im err 167.0% vs 190% pinned tol; resonance 128.9389 MHz vs
     143.0 MHz, err 9.8% vs 12% pinned tol) and the full root-cause chain
     (box/runtime convergence check, naive-PEC negative control sweep, the
     feed-model swap that disproved the aperture-port β term, the non-
     monotonic coarse/fine dx sweep) that pins Im(Z)/resonance at
     measured+margin rather than the aspirational 20%/5% targets;
   - an explicit FS.1c completion statement (GO, with the reactive gap
     named and root-caused, not hidden) and queued follow-ons (full
     telegrapher-coupled/charge-continuity model, orientations, junctions,
     GPU kernel, NTFF pattern gate — matching the spec's own non-goals).
2. **`docs/src/SUMMARY.md`** — one new line after the existing ADR-0227
   entry: `- [ADR-0228: FS.1c — thin-wire subcell (Holland–Simpson) + dipole
   gate vs NEC-4](decisions/0228-fs1c-thin-wire.md)`.
3. **`FULL-SUITE-ROADMAP.md`** — FS.1 row rewritten: appended the FS.1c
   SHIPPED sentence (citing ADR-0228, the subcell model one-liner, the
   coarse/fine consistency number, and the gate's measured Re(Z)/Im(Z)/
   resonance numbers with the root-cause summary), updated the gate-sketch
   column to add a `✓ (engine-thinwire-dipole-001)` checkmark, and changed
   the status column from `FS.1a + FS.1b COMPLETE, FS.1c queued` to
   `FS.1a + FS.1b + FS.1c COMPLETE` — **FS.1 is now COMPLETE** per the gate
   having shipped in Task 2.

No `crates/`, `.github/`, or other lane files touched — confirmed by
`git status --short` before staging (only the three doc paths above were
modified/created; unrelated untracked files from other concurrent agents in
`.superpowers/sdd/` were left alone, not staged).

## Verify (real output, this session)

- `git branch --show-current` → `feature/fs1c-thin-wire` (confirmed before
  the commit).
- `cargo clippy --workspace --all-targets -- -D warnings` — clean.
- `cargo clippy -p yee-compute --all-targets --no-default-features --
  -D warnings` — clean.
- `cargo fmt --check --all` — clean (no diff).
- `cargo doc -p yee-compute -p yee-engine --no-deps` — no `missing_docs`
  warnings (two pre-existing, unrelated `rustdoc::broken_intra_doc_links`/
  `private_intra_doc_links` warnings on `AperturePort::record` and
  `STEPS_PER_SUBMIT` links in `cpu.rs`/`gpu.rs` are untouched by this
  Task 3 docs-only diff — no Rust source file was modified this task).
- Bit-exact suite: `cargo test -p yee-compute --release --test
  graded_uniform_bitexact --test gpu_graded_parity --test gpu_cpu_parity --
  --include-ignored` — **5/5 pass**. GPU evidence real (re-run with
  `--nocapture`): `compute-002: running on adapter 'NVIDIA GeForce RTX 5060
  Ti'`, `compute-020: running on adapter 'NVIDIA GeForce RTX 5060 Ti'` — not
  SKIPPED.
- `cargo test -p yee-compute --release` (full default suite): every test
  binary `test result: ok`, 0 failed (26 binaries incl. `thin_wire.rs` 3/3,
  `gpu_thinwire_rejected.rs` 1/1, plus the 1 doctest).
- Staged explicit paths only (`git add docs/src/decisions/0228-fs1c-thin-wire.md
  docs/src/SUMMARY.md FULL-SUITE-ROADMAP.md`); `git status --short` after
  staging showed exactly those three (2 `M`, 1 `A`) plus other agents'
  pre-existing untracked `.superpowers/sdd/*` files, left untouched.

## Commit

`028853b` — `docs: ADR-0228 — FS.1c thin-wire (Holland-Simpson) COMPLETE`.
3 files changed, 218 insertions(+), 1 deletion(-). Trailers per the global
constraints.

## FS.1c track status

**COMPLETE.** Task 1 (CPU subcell + GPU rejection + unit tests), Task 2
(`engine-thinwire-dipole-001` gate: Re(Z) meets target, Im(Z)/resonance
honestly pinned short with a named, well-evidenced root cause — the dropped
wire charge-continuity coupling), and this Task 3 (ADR-0228 + roadmap row)
are all landed on `feature/fs1c-thin-wire`. FS.1 (antenna catalog: FS.1a +
FS.1b + FS.1c) is now COMPLETE per `FULL-SUITE-ROADMAP.md`. No tolerance was
ever widened past a measured value; the Re(Z) STOP threshold (25% off
NEC-4) was never approached (measured 5.8%). Queued follow-ons (not
attempted, per the spec's non-goals): the full telegrapher-coupled
(charge/continuity) thin-wire model, arbitrary orientations/bent wires,
wire junctions, monopole ground planes, a GPU thin-wire kernel, loaded/
insulated wires, an NTFF pattern gate for the dipole.
