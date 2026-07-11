# ADR-0212: FS.6.0/6.1 — two-port network algebra (S↔T, cascade, de-embed, renormalize)

**Date:** 2026-07-11 · **Status:** accepted · **Track:** FS.6 (`FULL-SUITE-ROADMAP.md`)
**Spec:** `docs/superpowers/specs/2026-07-11-fs6-network-algebra-design.md`

## Context

Yee round-trips Touchstone (R.2) but could not compose networks — cascade
two .s2p blocks, de-embed a fixture. Composition is table-stakes for
commercial parity and pure closed-form linear algebra: an ideal instant-gate
walking skeleton.

## Decision

`yee_io::network` (yee-io owns the S-parameter types):

1. **Chain convention `[b1; a1] = T [a2; b2]`**, i.e.
   `T = (1/s21)[[−det S, s11], [−s22, 1]]`, chosen because connecting
   port 2 of A to port 1 of B makes the cascade a plain product
   `T_A · T_B` (derivation in the module docs). `t_to_s` inverts it.
2. **Explicit singularity handling.** `s21 = 0` (isolating network) has no
   chain matrix; `t22 = 0` has no S-image; a singular fixture cannot be
   de-embedded. All are `Error::Network` (new variant), never NaN.
3. **`cascade_files` is strict in FS.6.0**: identical frequency grids
   (relative 1e-12, no interpolation), identical z₀ (renormalization is
   FS.6.1), 2-ports only. Every rejection is a named error; a silent
   resample is how composition tools corrupt data.
4. `deembed_left(fixture, measured) = t_to_s(T_f⁻¹ · T_m)` — recovers the
   DUT from `measured = fixture · DUT`. The right-side variant is
   mechanical and lands with renormalization in FS.6.1.

## Gate — `net-001` (`crates/yee-io/tests/network_algebra.rs`, instant, GREEN first run)

S↔T round-trip 1e-15 on a non-reciprocal lossy matrix; thru is a two-sided
cascade identity; matched 3 dB + 3 dB pads → 6.000 dB with phases summed;
associativity to 1e-12; de-embed recovers the DUT to 1e-12; **series
impedances cascade to their sum** (the ABCD identity
`[[1,Z₁],[0,1]]·[[1,Z₂],[0,1]] = [[1,Z₁+Z₂],[0,1]]` in S-form — non-zero
reflections throughout, so this exercises the det-S terms the matched
cases cannot); `cascade_files` happy path + all three rejection paths.

## FS.6.1 (same date): renormalization + right de-embed

- `renormalize(s, z_old, z_new)` — both ports real z₀ → real z₀′. With
  `r = (z_new − z_old)/(z_new + z_old)` the Kurokawa power-wave transform
  reduces to the Möbius form **`S′ = (S − rI)(I − rS)⁻¹`** because the
  scalar port-normalization factors cancel when every port changes
  identically. Complex reference impedances (where the port factors do
  NOT cancel and convention wars begin) are out of scope until a consumer
  needs them.
- `deembed_right(measured, fixture)` = `t_to_s(T_m · T_f⁻¹)`.
- `renormalize_file(&File, z_new)`; `cascade_files` **stays strict** about
  z₀ — the caller renormalizes explicitly first. Silent renormalization
  hides unit mistakes; an explicit two-step reads as intent.

Gate `net-002` GREEN first run: same-z₀ renormalization bit-exact
identity; series-impedance closed form at 50 Ω renormalized to 75 Ω
equals the direct 75 Ω construction to 1e-13; 50→75→50 round-trip 1e-14;
right de-embed recovers the DUT; `renormalize_file` unblocks the strict
`cascade_files` rejection and round-trips through the File layer.

## Consequences

- FS.6.2: CLI (`yee net cascade a.s2p b.s2p`) and studio exposure;
  matching-network synthesis (L-section/stub from a measured Γ), verified
  by cascading the synthesized match against the antenna's measured S and
  asserting the improved S11 — the roadmap's FS.6 full-wave gate.
