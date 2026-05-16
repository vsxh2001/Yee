# yee-core — Validation

Unit-level only. No physics here; just type contracts and numerical constants.

## Cases (Phase 0)

| ID | Description | Reference | Tolerance |
|----|-------------|-----------|-----------|
| `core-001` | Speed of light c₀ | CODATA 2018: 299,792,458 m/s (exact) | bit-exact |
| `core-002` | Vacuum permittivity ε₀ | CODATA 2018 | ≤ 1e-12 rel |
| `core-003` | Vacuum permeability μ₀ | CODATA 2018 | ≤ 1e-12 rel |
| `core-004` | Free-space impedance η₀ = √(μ₀/ε₀) | Closed-form | ≤ 1e-9 rel |
| `core-005` | FreqRange iteration endpoints | n=2 → [start, stop]; n=1 → [start] | bit-exact |
| `core-006` | FreqRange rejects start > stop | `Error::Invalid` returned | — |

## Running

```bash
cargo test -p yee-core
```

Doc tests count. Examples in `README.md` must compile.

## CI

Runs on every PR. Phase 0 must pass before merge.
