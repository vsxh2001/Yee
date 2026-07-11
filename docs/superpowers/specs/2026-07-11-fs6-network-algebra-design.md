# FS.6.0 — Network algebra walking skeleton (design)

**Date:** 2026-07-11 · **Track:** FS.6 (`FULL-SUITE-ROADMAP.md`) · **ADR:** 0212

## Problem

Every commercial suite composes measured/simulated blocks: cascade two .s2p
files, de-embed a fixture, check the result. Yee reads and writes Touchstone
(`yee-io`, R.2) but cannot combine networks. FS.6.0 is the walking skeleton:
2-port S↔T conversion, cascade, de-embed — pure closed-form linear algebra
with textbook-identity gates, instant to test.

## Non-goals (FS.6.1+)

N-port cascade, renormalization to a different z₀, matching-network
synthesis, CLI/studio exposure, mixed-frequency-grid interpolation.

## Design

New module `yee_io::network` (yee-io owns the S-parameter types; the
algebra's natural consumers are Touchstone `File`s):

- 2-port S-matrix as `[Complex64; 4]` row-major `[s11, s12, s21, s22]` —
  matching the `File::data` flattening.
- **Transfer (chain) matrix convention:** `[b1; a1] = T [a2; b2]`, i.e.
  `T = (1/s21)·[[−det S, s11], [−s22, 1]]`, chosen because port 2 of A
  feeding port 1 of B gives exactly `T_cas = T_A · T_B` (derivation in the
  module docs). Inverse: `s21 = 1/t22, s11 = t12/t22, s22 = −t21/t22,
  s12 = det T/t22`.
- `s_to_t`, `t_to_s`, `cascade(a, b)` (via T), `deembed_left(fixture,
  measured)` = `T_f⁻¹·T_m` (recover the DUT when `measured = fixture·DUT`).
- `cascade_files(&File, &File) -> Result<File>`: per-frequency cascade with
  explicit rejections (`Error::Network`): non-2-port, mismatched frequency
  grids (exact, no interpolation in FS.6.0), mismatched z₀.
- Singularity: `|s21| = 0` (a perfectly isolating network) has no T-matrix;
  reject with `Error::Network` rather than emit NaN.

## Gates — `net-001` (`crates/yee-io/tests/network_algebra.rs`, instant)

1. S↔T round-trip exact to 1e-15 on a non-trivial matrix.
2. Thru identity: `cascade(thru, X) = X = cascade(X, thru)`.
3. Attenuator composition: two matched 3 dB pads cascade to 6.02 dB with
   phases summed.
4. Associativity: `(A·B)·C = A·(B·C)` to 1e-12.
5. De-embed inverse: `deembed_left(F, cascade(F, D)) = D` to 1e-12.
6. Mismatch physics: cascading a lossless mismatched-impedance-step twice
   reproduces the closed-form Fabry–Perot reflection at zero line length.
7. `cascade_files`: happy path on synthetic 2-port `File`s + the three
   rejection paths.

## Lane

`crates/yee-io/**`, `docs/**`, `FULL-SUITE-ROADMAP.md`.
