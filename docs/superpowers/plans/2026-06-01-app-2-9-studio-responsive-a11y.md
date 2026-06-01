# App.2.9 — Studio responsive + a11y + P1 polish — Plan

**Spec:** `2026-06-01-app-2-9-studio-responsive-a11y-design.md` · **ADR:** ADR-0151

## Lane
`crates/yee-studio-web/{assets/studio.css, src/stages.rs, src/main.rs, Dioxus.toml}` ONLY.
NO `yee-filter` / engine edits. NO aesthetic overhaul (that's App.3.0). Out of lane → finding.

## Base / worktree
New worktree off `main` (re-fetch first; main f04e914). Branch `feature/app-2-9-responsive-a11y`.

## Pattern files (READ FIRST — edit ONLY in the worktree)
- `crates/yee-studio-web/assets/studio.css` — the design system (tokens `:root`, `.rail`,
  `.canvas`, `.grid-2/3`, `.tgrid`, `.compare-table`, `.tcard`/`.tcard.sel`, `.rec-use`,
  `.export-row .btn.dl`, `--muted`). Add `@media` at the END; bump `--muted` + `.rail .item
  .lab`; style `.rec-use`; strengthen `.tcard.sel`. Keep teal/copper/substrate identity.
- `crates/yee-studio-web/src/stages.rs` — the rail rendering (`.rail .item`, the active/`on`
  state — add `aria-current`/`aria-label`/`aria-hidden`), the technique gallery card
  (`.tcard.sel` markup), the Export "Design summary" + the recommender "Use this" button
  markup, the yield display. (Dioxus rsx! — attrs go as `aria_current: "step"` etc.)
- `crates/yee-studio-web/src/main.rs` + `Dioxus.toml` — where the document root / index
  template is, for `lang="en"` (Dioxus 0.6: check `Dioxus.toml [web.app]` or a custom
  `index.html`; if neither cleanly supports it, set it on the top-level app element or surface
  the limitation).

## Steps
1. **Contrast tokens** (studio.css `:root`): pick `--muted` ≥ 4.5:1 on `#0b0d11` (compute +
   comment the ratio; `#8b95a1`≈5.0:1 is a good candidate); rail `.lab` 9px → 11px.
2. **Responsive** (studio.css, new `@media (max-width: 760px)` + `(max-width: 420px)`): rail →
   horizontal top strip (flex-row, full width, scrollable) or wraps; `.grid-2/.grid-3/.tgrid`
   → 1 col; tables → wrap in an overflow-x container; `.canvas` padding down; inputs/selects
   `max-width: 100%`; export buttons wrap. Verify NO page horizontal overflow at 390px.
3. **`.rec-use`** → accent CTA (background `--accent-bg`, border `--accent`, color `--accent`,
   hover lift — mirror `.btn.dl`). **`.tcard.sel`** → filled-accent badge / stronger inset
   glow so selection is unmistakable.
4. **A11y** (stages.rs): active rail item `aria_current: "step"`; each rail button
   `aria_label` = stage name; glyph spans `aria_hidden: "true"`. **`lang="en"`** on the root
   (main.rs / Dioxus.toml / index template).
5. **Yield-vs-PASS wording** (stages.rs Export/Verify): label the spec verdict "nominal
   (spec mask)"; when MC yield is low, show a caveat chip/line. Wording + markup only.

## Verify (FROM THE WORKTREE; `bash -c`; EXIT 0; quote output)
- `cd crates/yee-studio-web && dx build --platform web --release` → EXIT 0.
- `cargo test -p yee-studio-web` → green; `cargo clippy -p yee-studio-web --all-targets -- -D
  warnings`; `cargo fmt --check`; `cargo check --workspace`.
- **Playwright responsive re-check** (load the playwright MCP; serve the built bundle or use
  `dx serve` / a static server on the `dist`): navigate at 390px AND 1280px, screenshot both
  to `/home/hadassi/Code/Yee/.playwright-mcp/app-2-9/` — assert NO horizontal page overflow at
  390px (e.g. `browser_evaluate` `document.documentElement.scrollWidth <= window.innerWidth`)
  and the inputs are not clipped; desktop unchanged. Quote the overflow check result.
- `git -C /home/hadassi/Code/Yee status --porcelain crates/` EMPTY (main untouched).

Commit (worktree): `yee-studio-web: responsive layout + AA contrast + a11y + P1 polish
(App.2.9, ADR-0151)` + Co-Authored-By trailer.

## Escape hatch
Stay within the existing design system — do NOT restyle the aesthetic (no new fonts, palette
change, motion — that's App.3.0). Keep desktop visually unchanged. If `lang` can't be set
cleanly in the dx flow, surface it (don't hack). If the Playwright re-check can't run
(no browser), still ship the CSS/a11y changes + state the verification gap. NO engine edits.
Blocked > 40 min → surface + stop.

## Done when
Responsive (no 390px overflow, verified), `--muted` ≥ 4.5:1 + rail labels ≥ 11px, the a11y
attrs + `lang` present, `.rec-use`/`.tcard.sel` fixed, the yield-vs-PASS wording clarified;
dx build EXIT 0; tests/clippy/fmt green; desktop unchanged; diff = `crates/yee-studio-web/**`.
Then I (dispatcher) verify + adversarial-review + merge. App.3.0 redesign follows.
