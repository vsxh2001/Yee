# FS.4.1 — vias through multilayer stackups (plan)

**Spec:** `docs/superpowers/specs/2026-07-13-fs41-stackup-vias-design.md`

1. `yee-voxel`: `with_via_between(model, i, j, k_lo, k_hi)` (blind via,
   `E_z` edges `k_lo..k_hi` PEC) + `with_through_via_at_cell(model, i,
   j)` (full stack, `0..nz`); re-express `with_via_at_cell` as
   `with_via_between(model, i, j, 0, k_top)` — bit-identical, existing
   R.1 unit test and `engine-via-001` untouched.
2. Gate `voxel-stackup-002` in
   `crates/yee-voxel/tests/voxel_stackup_002.rs`: through + blind vias
   on the FS.4.0 3-layer lidded stack mask exactly the expected `E_z`
   cells (whole-mask exact set-count), neighbour columns untouched,
   back-compat delegation bit-identical.
3. Gate `engine-stackup-via-001` in
   `crates/yee-engine/tests/stackup_via.rs` (release, ignored; no
   `antenna_`/`patch_`/`inset_`/`design_loop_`/`graded_` prefix, so the
   blanket engine CI release step picks it up): symmetric stripline
   (`Stackup::symmetric_stripline`, b = 16 cells at dx = 0.2 mm), λ/4
   open stub vs the same stub shorted by a through-via, three runs, one
   grid. Run with `--nocapture` first; pin asserts from the measured
   notch/no-notch numbers. Budget ≤ ~15 min release.
4. ADR-0221 (`docs/src/decisions/0221-fs41-stackup-vias.md`, measured
   numbers pinned) + one `docs/src/SUMMARY.md` line.
5. Verify: `cargo test -p yee-voxel` && the release gate (`--ignored
   --nocapture`) && clippy `-D warnings` && `cargo fmt --check --all`.
   Roadmap FS.4 row update is out-of-lane → report finding.
