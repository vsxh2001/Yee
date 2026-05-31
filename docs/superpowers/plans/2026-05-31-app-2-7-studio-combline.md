# App.2.7 — Light the Combline technique in the studio — Plan

**Spec:** `2026-05-31-app-2-7-studio-combline-design.md` · **ADR:** ADR-0146

## Lane
`crates/yee-studio-web/src/{engine.rs, stages.rs}` (+ `main.rs` only if a `match Topology`
arm needs it). NO yee-filter edits (combline engine + layout already exist). Out of lane
→ finding.

## Base / worktree
New worktree off `main` (re-fetch first). Branch `feature/app-2-7-combline`.

## Pattern files (READ FIRST — edit ONLY in the worktree, not the main repo)
- `crates/yee-studio-web/src/engine.rs` — the `Topology::Hairpin` arm of `derive_geometry`
  (MIRROR it for Combline); `Designed` (the geometry fields); `topology_name` /
  `length_label`; the App.2.1 hairpin engine test `hairpin_card_routes_to_real_dimensioner`
  (MIRROR for the combline test). The `SolvedDistributed`/`ResonatorRow` shape.
- `crates/yee-studio-web/src/stages.rs` — `Topology` enum, `Stage::rail`, `technique_stage`
  (the Combline card — currently `selects: None`), `technique_status` (the `Combline |
  Interdigital => Soon` arm — split Combline out to Live), `topology_response`,
  `technique_topology`, `layout_stage` (the board + resonator-table render).
- `crates/yee-filter/src/dimension.rs` — `dimension_combline` /
  `dimension_combline_layout` + `ComblineDimensions` (line_width_m, theta0_rad,
  resonator_length_m, **loading_cap_f**, gaps_m, target_k). (Reference only — do NOT edit.)
- The spec §Method (the θ0=π/4 default + the C_L surfacing) + ADR-0146.

## Steps
1. `Topology::Combline` (stages.rs) + `Stage::rail(Combline) => &Stage::DISTRIBUTED`.
2. engine.rs: `derive_geometry` `Topology::Combline` arm (mirror Hairpin) — `let θ0 =
   std::f64::consts::FRAC_PI_4;` → `dimension_combline(project, θ0, &SUBSTRATE)` +
   `dimension_combline_layout(project, θ0, &SUBSTRATE)`; build the geometry fields +
   resonator rows; carry the loading cap `C_L` (a Designed field or surface via the
   resonator table). `topology_name`/`length_label` Combline arms. Keep the demo seed
   path (`design_demo`) on EdgeCoupled.
3. engine.rs: a NON-vacuous test (mirror the hairpin one) — `design_demo_from(demo_spec(),
   Combline)` → `layout` differs from edge-coupled AND hairpin; coupling/verdict shared;
   `C_L > 0` finite.
4. stages.rs: Combline card `selects: Some(Topology::Combline)`; `technique_status`
   `Combline => Live(Topology::Combline)`; `topology_response(Combline) => Bandpass`;
   `technique_topology` Combline arm; `layout_stage` shows the board + the loading-cap line.
5. main.rs: only if a `match topology` site (StageCanvas) needs the Combline arm — combline
   uses the DISTRIBUTED flow (not lumped/stepped), so it should fall in with edge-coupled/
   hairpin; confirm no `if lumped_flow`/`stepped_flow` branch mis-routes it.

## Verify (run these FROM THE WORKTREE; expected EXIT 0; quote output)
- `cargo test -p yee-studio-web` — the new non-vacuous combline routing test passes;
  existing tests unregressed. Quote "test result: ok".
- `cargo clippy -p yee-studio-web --all-targets -- -D warnings` ; `cargo fmt --check`.
- `cargo check --workspace`.
- `cd crates/yee-studio-web && dx build --platform web --release` → EXIT 0.
- `git -C /home/hadassi/Code/Yee status --porcelain crates/` is EMPTY (main untouched).

Commit (in the worktree): `yee-studio-web: light the Combline technique (App.2.7,
ADR-0146)` + the Co-Authored-By trailer.

## Escape hatch
The lit Combline card MUST route to the real `dimension_combline` / `dimension_combline_layout`
output (no stub/clone). Surface the loading cap honestly (do not fake it). If a Dioxus
borrow/signal or exhaustive-match issue blocks > 30 min, surface it. NEVER edit
yee-filter / the main repo. Blocked > 30 min → stop + surface.

## Done when
Combline is a live, routable studio technique driven by the real combline engine + layout;
the non-vacuous test (combline layout ≠ edge-coupled/hairpin, C_L>0) passes; dx build
EXIT 0; existing flows unregressed; clippy/fmt/check clean; diff = `crates/yee-studio-web/**`.
Then I (dispatcher) verify + adversarial-review + merge.
