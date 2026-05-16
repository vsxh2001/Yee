# yee-cli — Validation

## Cases — Phase 0

| ID | Description | Tolerance |
|----|-------------|-----------|
| `cli-001` | `yee --help` exit 0, contains every subcommand name | exit 0 |
| `cli-002` | `yee --version` matches workspace `version` | exact |
| `cli-003` | Each subcommand `--help` exits 0 | exit 0 |
| `cli-004` | Unknown subcommand exits with non-zero + suggestion text | error code non-zero |

## Cases — Phase 1

| ID | Description | Tolerance |
|----|-------------|-----------|
| `cli-101` | End-to-end: `yee mesh microstrip.step && yee run patch.toml && yee export results.h5` | exit 0 chain |
| `cli-102` | `yee validate mom` runs full Phase 1 MoM suite | report file produced |

## Running

```bash
cargo test -p yee-cli
```

Integration tests live under `tests/` and shell out to the built binary via `assert_cmd`.
