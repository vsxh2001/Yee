# yee-io — Validation

## Cases — Phase 0

| ID | Description | Tolerance |
|----|-------------|-----------|
| `io-001` | Read sample `.s2p` from Sonnet docs corpus | exact struct match |
| `io-002` | Write → read round-trip on each sample | bit-exact on floats to 12 sig figs |
| `io-003` | Passivity check on read (eigenvalues of S†S ≤ 1 + ε) | ε = 1e-9 |
| `io-004` | Reject malformed file with `Error::TouchstoneParse` | error message useful |
| `io-005` | Frequency unit dispatch (GHz / MHz / Hz / kHz) | round-trip equal |

## Cases — Phase 1

| ID | Description | Tolerance |
|----|-------------|-----------|
| `io-101` | STEP import of unit cube via `opencascade-rs` | surface area = 6.0 ± 1e-9 |
| `io-102` | KiCad PCB outline import → polygon comparison vs Gerber | ±10 µm |
| `io-103` | Touchstone v2 reader on Keysight sample | exact struct match |
| `io-104` | HDF5 field-array round-trip | bit-exact |

## Fixtures

`validation/fixtures/touchstone/` is checked in (small text files); CAD fixtures live in `validation/fixtures/cad/` and may include licensed sample files — see fixture-level `README` for provenance.

## Running

```bash
cargo test -p yee-io                       # Touchstone only
cargo test -p yee-io --features opencascade  # adds CAD (slow first build)
```
