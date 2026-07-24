# FS.3.2c implementation plan — import-twin measurement gate

Spec: `docs/superpowers/specs/2026-07-24-fs3-2c-import-twin-design.md`
Branch: `feature/fs3.2c-import-twin` (from `main` @ `3191190`+)
Lane: `crates/yee-export/**`, `crates/yee-engine/**`; Task 2 docs
(`docs/src/decisions/0229-*.md`, `docs/src/SUMMARY.md`, `FULL-SUITE-ROADMAP.md`).

## Global constraints (bind every task)

- Real GPU present; GPU evidence prints `adapter 'NVIDIA GeForce RTX 5060 Ti'`.
- After every functional commit, green + unmodified: bit-exact suite
  (`cargo test -p yee-compute --release --test graded_uniform_bitexact --test
  gpu_graded_parity --test gpu_cpu_parity -- --include-ignored`);
  `cargo test -p yee-export --release` (all gerber-rt gates);
  the native stub gate this work twins (identify it: `grep -rn "stub_notch"
  crates/yee-engine/tests/`).
- Never weaken any assertion. Twin deltas need root-cause comments before pinning.
- clippy `-D warnings` workspace + `-p yee-compute --no-default-features`;
  `cargo fmt --check --all`; missing_docs; explicit-path staging.
- Commit style: crate prefix, ≤72-char imperative subject, why-body, trailers:
Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01FFo41x449XDGJ7Xyds4L7M

## Task 1 — twin path + gate `engine-import-twin-001`

- Read first: `crates/yee-export/src/` (export writer + `import` + `gerber_to_outline`
  types), the studio `import_gerber` command (how outline+metadata → Layout is
  already framed), and the native stub fixture (`sparams_stub_notch` or the
  cheapest stub-notch gate — pick the cheapest that measures a notch).
- If an outline→Layout helper exists (studio path), reuse/lift it into the
  library seam and document; else write the minimal helper (outline polygons +
  substrate/port metadata → `Layout`).
- Gate in `crates/yee-engine/tests/import_twin.rs`:
  1. Build the native stub Layout (generator).
  2. Export → Gerber bytes → import → outline → rebuilt Layout.
  3. Structural assert: identical trace polygons (exact coordinates — the
     import is vertex-exact per gerber-rt-001; if ordering differs, normalize
     and document).
  4. Run the identical measurement on both Layouts (same fixture options
     derived from each Layout independently — no sharing of the native grid).
  5. Compare notch (freq, depth) and, if cheap, the full S-curve: bit-identical
     expected; assert what measurement shows with root-cause comment if nonzero.
- Runtime ≤ ~5 min release; `#[ignore]` + blanket CI pickup (verify --skip list).
- Verify per global constraints. One or two commits (helper, gate).

## Task 2 — ADR-0229 + roadmap row

- `docs/src/decisions/0229-fs32c-import-twin.md` per 0228's structure: the twin
  chain, structural + measured results (bit-identical or delta + cause), the
  no-stackup-in-Gerber API contract, FS.3 remainder (DXF only). SUMMARY.md line
  after 0228; FS.3 row append in FULL-SUITE-ROADMAP.md.
