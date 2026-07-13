# FS.5c — Studio exposure of yield analysis (plan)

**Spec:** `docs/superpowers/specs/2026-07-13-fs5c-studio-yield-design.md`

1. `studio/src-tauri/Cargo.toml`: add
   `yee-surrogate = { path = "../../crates/yee-surrogate" }` (mirrors the
   `yee-engine` path-dep idiom; the studio stays its own workspace).
2. `studio/src-tauri/src/yield_mc.rs`: `YieldRequest` / `YieldResponse`
   serde types + `yield_estimate_impl` — derive `L₀ = c/(2 f₀ √ε_eff)`,
   build `ToleranceSpec { nominal: [L₀, ε_r], sigma: [σ_L, σ_εr] }`, pass
   closure `|p| |f(p) − f₀| ≤ halfwidth` (non-physical draws fail), map
   `YieldEstimate` to explicit clamped Wilson bounds. In-module unit tests:
   determinism, ADR-0211-regime yield at seed 20260711, validation errors.
3. `studio/src-tauri/src/lib.rs`: `pub mod yield_mc`, thin
   `#[tauri::command] fn yield_estimate(req) -> Result<YieldResponse, String>`,
   register in `generate_handler!`.
4. `studio/src/App.tsx`: `YieldPanel` (exported for tests) with the
   ADR-0211 defaults (2.45 GHz, ε_r 4.4, σ_L 0.1 mm, σ_εr 0.05, ±40 MHz,
   n 10⁴, seed 20260711), Run button → `invoke("yield_estimate", { req })`,
   result line (yield %, CI [lo, hi] %, n_pass/n, L₀ mm); mount after
   `<ImportPanel />`.
5. `studio/src/yield.test.tsx`: gate `studio-yield-dom-001` — mock
   `@tauri-apps/api/core` with `vi.hoisted` + `vi.mock`; assert default
   form, exact invoke args (unit conversions), and rendered yield/CI.
6. Docs: ADR-0222 + one `docs/src/SUMMARY.md` line. (FULL-SUITE-ROADMAP FS.5
   row update is out of lane — surfaced in the track report.)
7. Verification (expect exit 0):
   `cd studio && npm install && npm run build && npx vitest run && cargo check --manifest-path src-tauri/Cargo.toml`.
