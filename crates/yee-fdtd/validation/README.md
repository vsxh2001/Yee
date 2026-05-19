# yee-fdtd — Validation

Phase 2 deliverable. Phase 0/1: no live cases. Below is the planned suite.

## Cases — Phase 2

| ID | Case | Reference | Tolerance |
|----|------|-----------|-----------|
| `fdtd-007 (uniform-fine smoke)` | Phase 2.fdtd.7.z slot driver — well-formedness + passivity | internal driver invariants | finite `f_res ∈ [4, 14] GHz`, `S_11 ≤ 0 dB`, populated `notes` (Track UUUUUUUU rewire — **un-`#[ignore]`'d**, retired into default CI) |
| `fdtd-007 (f_res, ±2%)` | Phase 2.fdtd.7.z uniform-fine slot — `f_res` | Maloney-Smith 1993 Fig. 9 (~8.9 GHz, TBD) | `\|df\|/f_ref ≤ 0.02` (test `#[ignore]`'d — measured 5.30 GHz, see status note) |
| `fdtd-007 (\|S_11\|, ±1 dB)` | Phase 2.fdtd.7.z uniform-fine slot — `\|S_11(f_res)\|` | Maloney-Smith 1993 Fig. 9 (~−22 dB, TBD) | `\|dS_11\| ≤ 1 dB` (test `#[ignore]`'d — measured −6.4 dB, see status note) |
| `fdtd-007 (subgrid sanity, 0.3% / 0.3 dB)` | subgridded vs globally-uniform `dx = 0.5 mm` reference, 5 spot frequencies | internal comparator | max `\|df\|/f ≤ 0.003` AND max `\|dS_11\| ≤ 0.3 dB` (test `#[ignore]`'d — subgridded variant retired by Track UUUUUUUU pending F2 inward-coupling restoration) |
| `fdtd-201` | Rectangular cavity TE/TM Q-factor | Analytical | ±0.5% |
| `fdtd-202` | Pyramidal horn antenna pattern | Measured / Balanis | ±1 dB main beam |
| `fdtd-203` | Dipole over dielectric half-space NTFF | Sommerfeld reference | analytic match |
| `fdtd-204` | Cross-validation vs openEMS | openEMS on identical grid | numerical-noise level |
| `fdtd-205` | Microstrip transient TDR | FFT(yee-mom Sxx) | ±2% |
| `fdtd-206` | Drude-metal plasmonic dipole | Maier / textbook | ±5% resonance |
| `fdtd-207` | Multi-pole Debye human-tissue benchmark | Gabriel database | ±5% absorption |

### `fdtd-007` (Phase 2.fdtd.7.z Track UUUUUUUU rewire) — status note

Driver: `yee_validation::run_fdtd_007_maloney_smith_slot`. Tests live
at `crates/yee-validation/tests/fdtd_007_maloney_smith_slot.rs`.
**Track UUUUUUUU (Phase 2.fdtd.7.z) retired the two structural
blockers from the original LLLLLLLL commit** by rewiring the driver
onto the new per-cell-ε map (`YeeGrid::with_eps_r_cells`; MMMMMMMM
commit `cb6f8ed`), per-component PEC mask
(`YeeGrid::with_pec_mask_e{y,z}`; same commit), and CPML-per-cell-ε
coupling (PPPPPPPP commit `c57592f`). The 2.2-substrate slab is now a
true per-cell region, the ground plane is a per-component PEC sheet
with the slot rectangle cut out, and both half-spaces are
CPML-terminated for radiation.

**Status after rewire:**

- *Uniform-fine smoke gate (`fdtd_007_uniform_fine_smoke`)* —
  **retired from `#[ignore]` into default CI** (Track UUUUUUUU).
  Confirms the driver returns a well-formed result with a finite
  `f_res ∈ [4, 14] GHz` and a passive `S_11 ≤ 0 dB`. Wall-time ~3 s
  release / ~40 s debug on a `32 × 80 × 25` grid.
- *Maloney-Smith physics gates (`fdtd_007_fres_within_two_percent_*`,
  `fdtd_007_s11_within_one_db_*`, `fdtd_007_within_two_percent_and_one_db`)*
  — **still `#[ignore]`'d**. Measured `f_res ≈ 5.30 GHz` (|`S_11(f_res)`|
  ≈ −6.4 dB) is `|df|/f ≈ 0.40` off the digitised `8.9 GHz` reference,
  outside the `±5 %` digitisation envelope. Per the LLLLLLLL TBD
  escape hatch (`Do NOT relax to > 5 %`), the gates stay
  `#[ignore]`'d with the measured value documented in the ignore
  reason and the module docstring. Root-cause candidates: (a) the
  Phase 2.fdtd.7 spec cites Maloney & Smith 1993 IEEE T-AP 41(5),
  which is a *cylindrical monopole* paper — the Fig. 9 attribution
  itself may be wrong; (b) the measured 5.3 GHz resonance is
  consistent with a half-wave slot mode in the `ε_eff ≈ 1.6` slab
  approximation, suggesting the published `8.9 GHz` may be a higher
  mode or different geometry. Both follow-ups land as `fdtd-007.1`
  (radiation-CPML / reference-figure verification).
- *Subgridded sanity gate (`fdtd_007_subgrid_vs_uniform_sanity_check`)*
  — **still `#[ignore]`'d**. The Phase 2.fdtd.7.y Step C6 un-ghosted-J
  Berenger closure leaves the fine sub-grid effectively passive on
  source-on-coarse drive (LLLLLLLL item 3); the rewired driver no
  longer runs a subgrid-vs-uniform comparator, so this gate cannot
  execute its intended check until the F2 inward-coupling restoration
  (deferred from Track DDDDDDDD) lands.

The Maloney-Smith reference numerics remain digitised from Fig. 9 to
±5 % per the original Q7 escape hatch (see `FDTD_007_FRES_REF_HZ` /
`FDTD_007_S11_DB_REF` doc-comments — both flagged TBD pending
journal-figure verification).

## Running

Will require GPU + CUDA toolkit + `yee-cuda` feature `cuda`.

```bash
cargo test -p yee-fdtd --release --features cuda
```

## Cross-tool validation

For every case where openEMS or gprMax can run the same geometry, we publish side-by-side results in `validation/results/` so a reader can verify our work without trusting our numbers.
