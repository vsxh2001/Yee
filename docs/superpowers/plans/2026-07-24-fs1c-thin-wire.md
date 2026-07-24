# FS.1c implementation plan — Holland–Simpson thin wire + dipole gate

Spec: `docs/superpowers/specs/2026-07-24-fs1c-thin-wire-design.md`
Branch: `feature/fs1c-thin-wire` (from `main` @ `6dec5a6`+)
Lane: `crates/yee-compute/**`, `crates/yee-engine/**`; Task 3 docs
(`docs/src/decisions/0228-*.md`, `docs/src/SUMMARY.md`, `FULL-SUITE-ROADMAP.md`).

## Global constraints (bind every task)

- Real GPU present; GPU evidence prints `adapter 'NVIDIA GeForce RTX 5060 Ti'`.
- After every functional commit, green + unmodified: bit-exact suite
  (`cargo test -p yee-compute --release --test graded_uniform_bitexact --test
  gpu_graded_parity --test gpu_cpu_parity -- --include-ignored`) and
  `cargo test -p yee-compute --release` (full default crate suite).
- No-wire jobs are a provable no-op (existing results bit-identical). Never
  weaken any assertion. Honest pins; Re(Z) > 25 % off NEC-4 → STOP and
  root-cause, never widen.
- Research-first: cite the exact published formulation used (source + equation)
  in module docs; do not invent update coefficients.
- clippy `-D warnings` workspace + `-p yee-compute --no-default-features`;
  `cargo fmt --check --all`; missing_docs; explicit-path staging.
- Commit style: crate prefix, ≤72-char imperative subject, why-body, trailers:
Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01FFo41x449XDGJ7Xyds4L7M

## Task 1 — thin-wire subcell model, CPU (yee-compute)

- Research first (WebFetch/WebSearch allowed): Holland–Simpson formulation as
  summarized in Taflove §10.3 or gprMax's thin-wire docs. Record the chosen
  update equations in the module docs with the citation BEFORE implementing.
- Read `crates/yee-compute/src/{drive.rs,materials.rs,cpu.rs}` to pick the seam
  (Drive vs Materials attachment) — document the choice + why in the report.
- Implement z-axis-only `ThinWire` per spec: modified Ez update on wire edges
  (in-cell inductance), radial E shorted at wire cells, optional feed gap cell
  (the feed cell keeps the normal/driven update — delta-gap idiom).
- GPU: `ComputeError::Unsupported` named rejection when a drive/materials
  carries thin wires, with a test (pattern: aperture-port GPU rejection).
- Unit tests: (a) no-wire no-op bit-identity (run a small vacuum job with and
  without the (empty) wire list — byte-equal fields); (b) wire-present smoke:
  fields stay finite, Ez along the wire ≠ free-space run (the wire does
  something); (c) the coarse/fine resonance consistency check from the spec —
  same physical dipole at dx and dx/√2, resonant frequencies within a pinned
  few % (measure first, pin honestly), plus (if cheap) the naive one-cell-PEC
  negative control showing worse grid-dependence.
- Verify per global constraints. One or two commits.

## Task 2 — gate `engine-thinwire-dipole-001` (yee-engine)

- Pattern files: an existing antenna gate for the open-boundary + CPML idiom
  (`grep -rn "engine-antenna" crates/yee-engine/tests/`), FS.2a port-records
  for feed V/I.
- Fixture: L = 1 m, a = 5 mm z-dipole centered in an open CPML box (all faces;
  box clearance ≥ λ/4 at 143 MHz — compute it), coarse λ/20 grid via the
  automesh wavelength rule where applicable (document dx). Broadband pulse over
  ~100–200 MHz; feed V/I → Z(f).
- Assertions: resonance frequency (Im Z zero-crossing or |Z| min) within ±5 %
  of expectation; Re(Z) at resonance vs NEC-4 87 Ω within pinned tol (target
  ≤10 %); Im(Z) at the NEC-4 resonance vs +41 Ω within pinned tol (target
  ≤20 %). Print all measured values. STOP rule per spec.
- Runtime ≤ ~3 min release; `#[ignore]` + blanket CI pickup (verify --skip list).
- Verify per global constraints + all three stripline gates still green
  (`--test stripline_eeff --test stripline_z0 --test stripline_alpha`). One commit.

## Task 3 — ADR-0228 + roadmap row

- `docs/src/decisions/0228-fs1c-thin-wire.md` per 0227's structure: the chosen
  formulation + citation, seam decision, coarse/fine consistency numbers,
  measured Z vs NEC-4, FS.1 completion statement, queued follow-ons (GPU port,
  orientations, junctions, NTFF pattern gate). SUMMARY.md line after 0227;
  FS.1 row update in FULL-SUITE-ROADMAP.md.
