# Yee Studio (Tauri 2 + React)

Desktop studio over the `yee-engine` in-process job API (ADR-0179; S.2 walking skeleton).

```bash
# prerequisites (Linux): libwebkit2gtk-4.1-dev libgtk-3-dev librsvg2-dev
npm install
npm run tauri dev      # interactive dev app (needs a display)
npm run tauri build    # release bundle
```

`src-tauri/` is deliberately **not** a member of the repo's root cargo workspace so the
webkit2gtk dependency tree never weighs down workspace-wide builds. The frontend talks to the
engine through one Tauri command (`run_job`) plus `job://progress` events — the same serde
protocol `yee-server` (S.1) will expose over WebSocket.
