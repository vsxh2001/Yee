# ADR-0180: S.1 `yee-server` (WebSocket job API) + S.3 visualization walking skeleton

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0179 (S.0 protocol + Tauri shell), ADR-0175 (direction), `ENGINE-STUDIO-ROADMAP.md`.

---

## S.1 — `yee-server`

New workspace crate exposing the `yee-engine` protocol over HTTP/WebSocket with axum 0.8:

- `GET /healthz` — liveness.
- `GET /v1/jobs` (WS) — the client sends one JSON `JobSpec` text frame; the server streams
  every `JobEvent` back as JSON text frames and closes after `done`/`error`. Events are
  forwarded **live** (a `spawn_blocking` bridge re-sends the engine's std-channel events into a
  tokio channel), and a client disconnect mid-run **cancels the job cooperatively** via the new
  `JobHandle::canceller()` / `JobCanceller` (the cancel flag detaches from the handle so the
  handle itself can move into the bridge).
- `yee serve --addr 127.0.0.1:7332` in `yee-cli` wraps `serve_blocking` (fresh tokio runtime —
  the CLI stays sync). Verified live in-container: `/healthz` answers `ok`.

The wire format is byte-identical to what the Tauri studio uses in-process — one serde
protocol, now on its second transport, exactly as ADR-0089 originally called for
("heavy EM on a native yee-server the web client calls") and ADR-0179 planned.

**Gates** (`crates/yee-server/tests/ws_end_to_end.rs`, run in the workspace suite): a real
tokio-tungstenite client against an ephemeral-port server — (1) full round trip asserting
≥ 10 streamed `progress` events, the `done` payload's probe series (length, non-trivial
signal), and the requested field slice; (2) an invalid spec yields a structured `error` event.

## S.3 — visualization walking skeleton

- **Engine:** `JobSpec.slice: Option<SliceSpec { component, k }>` (serde-defaulted — old specs
  parse unchanged) returns `JobResult.slice: Option<FieldSlice { ni, nj, data }>`, the final
  z-plane of an E component, on both backends. Unit-gated (`slice_is_returned_when_requested`).
- **Studio:** two new dependency-free views fed by the job result: `SpectrumPlot` (single-bin
  DFT magnitude scan of the probe series, peak annotated) and `SliceHeatmap` (canvas heatmap of
  the E_z mid-plane, diverging blue-white-red). The analysis code (`analysis.ts`) is pure
  TypeScript, gated by vitest against a **known reference**: the DFT scan must recover a pure
  sinusoid's frequency to within one bin; the color map is gated at its extremes/centre/clamp.
  DOM-level smoke gates (vitest + jsdom + testing-library) render both components from fixture
  data — the "DOM-level smoke gates" the S.3 roadmap row called for. 7 tests green; bundle
  grew 149.6 kB → still **48.7 kB gzipped**.
- **CI:** new `studio-build` job (Node 22 + webkit2gtk): frontend typecheck + vite build +
  vitest + `cargo check` of the Tauri shell — closing the "CI wiring" follow-on from ADR-0179.
  `yee-server` tests ride the existing workspace jobs.

## Scope notes

- Per-event WS forwarding is in place; **binary frames** (for large field payloads) and
  multi-job sockets are deferred until a consumer needs them.
- three.js 3-D volumetric rendering remains the S.3 follow-on (`S.3b`); the slice heatmap +
  spectrum are the walking skeleton that proves the engine-stream → view pipe.
- S.4 (Dioxus parity audit) is now the only queued studio phase.
