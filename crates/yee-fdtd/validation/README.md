# yee-fdtd — Validation

Phase 2 deliverable. Phase 0/1: no live cases. Below is the planned suite.

## Cases — Phase 2

| ID | Case | Reference | Tolerance |
|----|------|-----------|-----------|
| `fdtd-007 (f_res, ±2%)` | Phase 2.fdtd.7 Q7 subgridded slot — `f_res` | Maloney-Smith 1993 Fig. 9 (~8.9 GHz, TBD) | `\|df\|/f_ref ≤ 0.02` (test `#[ignore]`'d per Q7 C6 escape hatch) |
| `fdtd-007 (\|S_11\|, ±1 dB)` | Phase 2.fdtd.7 Q7 subgridded slot — `\|S_11(f_res)\|` | Maloney-Smith 1993 Fig. 9 (~−22 dB, TBD) | `\|dS_11\| ≤ 1 dB` (test `#[ignore]`'d per Q7 C6 escape hatch) |
| `fdtd-007 (subgrid sanity, 0.3% / 0.3 dB)` | subgridded vs globally-uniform `dx = 0.5 mm` reference, 5 spot frequencies | internal comparator | max `\|df\|/f ≤ 0.003` AND max `\|dS_11\| ≤ 0.3 dB` (test `#[ignore]`'d per Q7 C6 escape hatch) |
| `fdtd-201` | Rectangular cavity TE/TM Q-factor | Analytical | ±0.5% |
| `fdtd-202` | Pyramidal horn antenna pattern | Measured / Balanis | ±1 dB main beam |
| `fdtd-203` | Dipole over dielectric half-space NTFF | Sommerfeld reference | analytic match |
| `fdtd-204` | Cross-validation vs openEMS | openEMS on identical grid | numerical-noise level |
| `fdtd-205` | Microstrip transient TDR | FFT(yee-mom Sxx) | ±2% |
| `fdtd-206` | Drude-metal plasmonic dipole | Maier / textbook | ±5% resonance |
| `fdtd-207` | Multi-pole Debye human-tissue benchmark | Gabriel database | ±5% absorption |

### `fdtd-007` (Phase 2.fdtd.7 Q7 — escape-hatched) — status note

Driver: `yee_validation::run_fdtd_007_maloney_smith_slot`. Tests live
at `crates/yee-validation/tests/fdtd_007_maloney_smith_slot.rs`. All
three hard gates above are `#[ignore]`'d in default CI per the Q7
escape hatch — the Phase 2.fdtd.7.y Step C6 un-ghosted-J variant
leaves the fine sub-grid effectively passive under source-on-coarse
drive, so the Maloney-Smith dielectric-loaded resonance cannot be
recovered until the F2 inward-coupling restoration (deferred from
Track DDDDDDDD) lands as Phase 2.fdtd.7.z. Two further structural
blockers compound the C6 trade-off — `YeeGrid` exposes a scalar
`eps_r` (no per-cell map for the `ε_r = 2.2` substrate) and there is
no per-cell PEC mask for the slot-in-ground-plane geometry. The
driver compiles, runs end-to-end, and returns well-formed
measurements; it is the scaffolding the Phase 2.fdtd.7.z follow-up
work will plug into.

The Maloney-Smith reference numerics are digitised from Fig. 9 to
±5 % per the Q7 escape hatch (see `FDTD_007_FRES_REF_HZ` /
`FDTD_007_S11_DB_REF` doc-comments — both flagged TBD pending
journal-figure verification).

## Running

Will require GPU + CUDA toolkit + `yee-cuda` feature `cuda`.

```bash
cargo test -p yee-fdtd --release --features cuda
```

## Cross-tool validation

For every case where openEMS or gprMax can run the same geometry, we publish side-by-side results in `validation/results/` so a reader can verify our work without trusting our numbers.
