# yee-cli — Roadmap

## Phase 0 (months 0–6)
- [ ] Subcommand skeleton: `validate`, `mesh`, `run`, `export`
- [ ] `tracing_subscriber` setup with `RUST_LOG` env var
- [ ] `--help` text covers every subcommand
- [ ] Smoke test: every subcommand exits 0 on canonical inputs
- [ ] Shell completions generated via `clap_complete`

## Phase 1 (months 6–18)
- [ ] `yee validate mom` actually runs the MoM validation suite + emits a report
- [ ] `yee mesh <step>` produces a Yee project file with the meshed geometry
- [ ] `yee run <project.toml>` dispatches to MoM / FDTD per the project type
- [ ] `yee export <results.h5> --format touchstone` works for n-port S-data
- [ ] `yee bench` runs the standard benchmark suite + compares vs published numbers

## Phase 2+ (months 18+)
- [ ] `yee gui` launches the egui desktop GUI
- [ ] `yee surrogate train | predict` for ML surrogate workflow (Phase 3)
- [ ] `yee chat` invokes the natural-language design surface (Phase 3)

## Validation gates
- Phase 0: `yee --help` and each `yee <cmd> --help` exits 0 on Linux + Windows.
- Phase 1: end-to-end smoke runs `yee mesh ... && yee run ... && yee export ...` on a published microstrip case without manual intervention.
