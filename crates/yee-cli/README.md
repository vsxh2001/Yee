# yee-cli

> The `yee` command-line tool.

## Subcommands

| Command | Purpose | Phase |
|---------|---------|-------|
| `yee validate [solver]` | Run validation suite (`mom`, `fdtd`, or `all`) | 0 (stub) → 1 (real) |
| `yee mesh <input>` | Mesh STEP / IGES / KiCad PCB via Gmsh | 0 (stub) → 1 |
| `yee run <project>` | Run a simulation from a TOML project file | 1 |
| `yee export <results>` | Export Touchstone / HDF5 | 0 (stub) → 1 |

## Installing (post-release)

```bash
cargo install --path crates/yee-cli --features cuda,gmsh
```

## Feature flags

| Flag | Effect |
|------|--------|
| `cuda` | Enable GPU paths in `yee-mom` and `yee-fdtd` |
| `gmsh` | Link the meshing backend |

Both can be combined freely.

## Phase 0 status

The CLI's job in Phase 0 is to **prove the pipes are connected** — Cargo workspace → solver crate → I/O crate → CLI entrypoint. Subcommands print what they intend to do and exit 0. Real behavior arrives Phase 1 as solver features land.

## Roadmap

See [`ROADMAP.md`](ROADMAP.md).
