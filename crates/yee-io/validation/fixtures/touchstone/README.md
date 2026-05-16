# Touchstone v1.1 fixtures

Small, hand-curated `.sNp` files used by the `yee-io` test suite. All fixtures
are synthetic so we can keep them in-tree under the workspace license.

Spec reference: <https://ibis.org/connector/touchstone_spec11.pdf>.

## Per-fixture provenance

### `1port.s1p`
- **Source:** synthetic.
- **Network:** a frequency-independent 75-Ω resistive load referenced to
  Z₀ = 50 Ω.
- **Computation:** S11 = (Z_L − Z₀) / (Z_L + Z₀) = (75 − 50) / (75 + 50) = 0.2.
- **Purpose:** exercises the simplest `.s1p` path (one complex per row),
  GHz unit dispatch, and the default RI format.

### `2port.s2p`
- **Source:** synthetic.
- **Network:** ideal matched 6 dB attenuator referenced to Z₀ = 50 Ω.
- **Computation:** S11 = S22 = 0; S21 = S12 = 10^(−6/20) ≈ 0.501187…
- **Purpose:** exercises the n = 2 on-disk `S11 S21 S12 S22` reordering,
  GHz dispatch, RI format. Sigma_max(S) ≈ 0.5012 < 1, so passive.

### `2port_db.s2p`
- **Source:** synthetic; same physical network as `2port.s2p`.
- **Encoding:** DB format, MHz frequency axis. S11/S22 are encoded as −200 dB
  (a finite stand-in for the ideal −∞; |S| ≈ 1e-10 round-trips bit-identically
  in DB through our codec).
- **Purpose:** exercises the DB codec and MHz unit dispatch end-to-end. Read
  →  this fixture and `2port.s2p` produce numerically equivalent `File`
  structs up to the documented `-200 dB ≠ 0` offset.

### `2port_nonpassive.s2p`
- **Source:** synthetic.
- **Network:** every S-parameter set to 2 + 0j at f = 1 GHz.
- **Purpose:** triggers the passivity check in `touchstone::read`. The
  expected result is `Error::TouchstoneParse` whose message mentions
  "passivity".
