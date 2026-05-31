# App.2.8 — Light the Interdigital technique in the studio — Plan

**Spec:** `2026-05-31-app-2-8-studio-interdigital-design.md` · **ADR:** ADR-0150

## Lane
`crates/yee-studio-web/src/{engine.rs, stages.rs}` (+ `main.rs` ONLY if a `match Topology`
arm needs it). NO `yee-filter` edits (interdigital engine + layout already exist + are gated).
Out of lane → finding.

## Base / worktree
New worktree off `main` (re-fetch first; main post-83240b5). Branch
`feature/app-2-8-interdigital`.

## Pattern files (READ FIRST — edit ONLY in the worktree, never the main repo)
- `crates/yee-studio-web/src/engine.rs` — MIRROR the `Topology::Combline` sites: the
  `derive_geometry` Combline arm (~355, calls `dimension_combline`+`_layout`, `loading_cap_f:
  Some`); `topology_name`/`length_label` (~171/188); the two distributed match groups (~1147
  `topbar_view`, ~1295 `verify_view`); `compare_techniques` (~1416/1456, the hardcoded
  combline row) + its test (~2137, `rows[2]==Combline`); `overlay_curves` (~1525). For
  interdigital: call `dimension_interdigital`+`dimension_interdigital_layout` (NO theta0),
  `loading_cap_f: None`.
- `crates/yee-studio-web/src/stages.rs` — MIRROR the `Topology::Combline` sites: the `Topology`
  enum (~55); `Stage::rail` (~115); `technique_status` (~785 Combline→Live; flip Interdigital
  ~788 from `Soon(EdgeCoupled)` to `Live(Topology::Interdigital)`); `technique_label` (~812);
  `topology_response` (~834, band-pass group); `technique_topology` (~1167); the gallery card
  `selects` (~1397); `layout_stage` distributed group (~1603).
- `docs/superpowers/specs|plans/2026-05-31-app-2-7-studio-combline*` — the combline-lighting
  precedent (this is its interdigital twin). The spec §Method lists all 12 sites + main.rs.

## Steps
1. `stages.rs`: `Topology::Interdigital` variant + `Stage::rail` DISTRIBUTED + `technique_status`
   Live + `technique_label "Interdigital"` + `topology_response` Bandpass + `technique_topology`
   + the gallery card `selects: Some(Topology::Interdigital)`.
2. `engine.rs`: the `derive_geometry` Interdigital arm (mirror Combline; `dimension_interdigital`
   + `_layout`, `loading_cap_f: None`); `topology_name`/`length_label` arms; add
   `| Topology::Interdigital` to the topbar + verify + layout distributed groups; add the
   interdigital row to `compare_techniques` + fix the index test; handle `overlay_curves`
   (interdigital shares the coupled-resonator ideal — no new distinct curve; verify the
   overlay test still holds).
3. `main.rs`: only if a `match topology` StageCanvas site is non-exhaustive — add the
   Interdigital arm (distributed flow, like Combline). Confirm no lumped/stepped branch
   mis-routes it.
4. The non-vacuous routing test (mirror App.2.7's combline test): interdigital layout ≠
   edge-coupled AND hairpin AND combline; shared synthesis/verdict; `combline_loading_cap_f`
   None; resonator length > 0.

## Verify (run FROM THE WORKTREE; expected EXIT 0; quote output)
- `cargo test -p yee-studio-web` — the new non-vacuous interdigital routing test passes;
  existing tests unregressed (incl. the updated compare index test). Quote "test result: ok".
- `cargo clippy -p yee-studio-web --all-targets -- -D warnings` ; `cargo fmt --check`.
- `cargo check --workspace`.
- `cd crates/yee-studio-web && dx build --platform web --release` → EXIT 0. (If `dx`/wasm
  target/dioxus-cli missing, install per CLAUDE.md §7 / the App.2.7 notes: `cargo install
  dioxus-cli --version "^0.6"`, `rustup target add wasm32-unknown-unknown`.)
- `git -C /home/hadassi/Code/Yee status --porcelain crates/` EMPTY (main untouched).

Commit (in the worktree): `yee-studio-web: light the Interdigital technique (App.2.8,
ADR-0150)` + the Co-Authored-By trailer.

## Escape hatch
The lit Interdigital card MUST route to the real `dimension_interdigital` /
`dimension_interdigital_layout` output (no stub/clone, no combline-with-cap). Surface the
λg/4 resonator length honestly; do NOT invent a loading cap (interdigital has none →
`combline_loading_cap_f` stays None). NO `yee-filter` / main-repo edits. If a Dioxus
borrow/signal or exhaustive-match issue blocks > 30 min, surface it. Blocked > 30 min → stop +
surface.

## Done when
Interdigital is a live, routable studio technique driven by the real interdigital engine +
layout; the non-vacuous test (interdigital layout ≠ edge-coupled/hairpin/combline, no cap,
λg/4 length) passes; `dx build` EXIT 0; existing flows unregressed; clippy/fmt/check clean;
diff = `crates/yee-studio-web/**`. Then I (dispatcher) verify + adversarial-review + merge. The
gallery is then COMPLETE.
