# FS.4.2c implementation plan — automesh stackup integration

Spec: `docs/superpowers/specs/2026-07-24-fs4-2c-automesh-stackup-design.md`
Branch: `feature/fs4.2c-automesh-stackup` (from `main` @ `489fa23`+)
Lane: `crates/yee-engine/**`; Task 2 docs (`docs/src/decisions/0227-*.md`,
`docs/src/SUMMARY.md`, `FULL-SUITE-ROADMAP.md`).

## Global constraints (bind every task)

- Real GPU present; GPU evidence prints `adapter 'NVIDIA GeForce RTX 5060 Ti'`.
- After every functional commit, green + unmodified: bit-exact suite
  (`cargo test -p yee-compute --release --test graded_uniform_bitexact --test
  gpu_graded_parity --test gpu_cpu_parity -- --include-ignored`),
  `cargo test -p yee-engine --release` (default set), and the three stripline
  gates (`stripline_eeff`, `stripline_z0`, `stripline_alpha`, each `--release
  -- --include-ignored`; budget ~3 min total — alpha alone is ~116 s).
- Never weaken any assertion. Honest pins; rulebook-grid ε_eff > 5 % off exact →
  STOP and root-cause the rule, don't widen.
- clippy `-D warnings` workspace + `-p yee-compute --no-default-features`;
  `cargo fmt --check --all`; missing_docs clean; explicit-path staging only.
- Commit style: `yee-engine:`/`docs:` prefix, ≤72-char imperative subject,
  why-body, trailers exactly:
Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01FFo41x449XDGJ7Xyds4L7M

## Task 1 — `auto_dx_stackup` + unit tests (yee-engine)

- Read `crates/yee-engine/src/automesh.rs` `auto_dx`/`min_feature_m` (doc style,
  clamp, rule structure) and `yee_layout::Stackup`/`StackupLayer`.
- Implement per the spec's five rules; doc comment mirrors `auto_dx`'s per-rule
  bullet style, cites ADR-0215 for the b/16 lid rule.
- Unit tests in automesh.rs's test module (pattern: the existing `auto_dx_*`
  tests): one per binding rule (λ/20 via a high-ε_r layer; per-layer h/3 via one
  thin layer; feature/2 via a narrow trace; b/16 via a lidded stack where it
  binds; clamp), plus single-layer-no-lid ≡ `auto_dx` exact-equality consistency.
- Verify: `cargo test -p yee-engine --release` (unit set); clippy/fmt. One commit.

## Task 2 — gate `engine-automesh-stackup-001`

- Pattern files: `crates/yee-engine/tests/stripline_eeff.rs` (fixture + ε_eff
  extraction + exact TEM reference) and the FS.0 no-hand-dx gate idiom
  (`grep -rn "automesh" crates/yee-engine/tests/`).
- New test file `automesh_stackup.rs`: build the stripline fixture's Stackup +
  Layout, call `auto_dx_stackup` (NO hand dx anywhere — derive every
  cell-denominated quantity from the returned dx in metres, the ADR-0204
  constant-physics hygiene), eprintln which rule binds + the dx, run the ε_eff
  measurement, assert vs exact TEM ε_r within the honest pinned tolerance
  (target ≤ 2 %, the eeff-gate bar; pin measured + margin). Assert the binding
  rule is the expected one (lid b/16 for this fixture — compute expectation in
  the test, don't hardcode blindly; if a different rule binds, the assert
  message must say which and why that matters).
- Runtime target ≤ ~60 s release; `#[ignore]` + blanket yee-engine
  `--include-ignored` CI job pickup (verify via the --skip list grep as FS.4.2b
  did). One commit.

## Task 3 — ADR-0227 + roadmap row

- `docs/src/decisions/0227-fs42c-automesh-stackup.md` per 0226's structure:
  rules + rationale, which rule binds on the reference fixture, measured ε_eff
  + dx, FS.4.2 remainder (MoM cross-check; graded auto_spacings stackup
  variant). SUMMARY.md line after 0226; FS.4 row append in FULL-SUITE-ROADMAP.md.
