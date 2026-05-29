# Coupled Lines & Dimensional Synthesis — Theory of Operation

This page is the theory-of-operation reference for the second half of
Yee's RF filter-design flow: turning the *abstract* coupling matrix
produced by synthesis into the *physical* edge-coupled-microstrip
dimensions (line width, resonator length, and inter-resonator gaps) that
realise it. It is the downstream companion to the
[RF Filter Synthesis](./filter-synthesis.md) page — read that one first,
because everything here consumes its coupling coefficients `k_{i,i+1}`
and external quality factors `Qe`. Same audience (an engineer reading
source code with a textbook open), same conventions (plain-text math
because the documentation build does not render LaTeX, inline citations,
source-file references in inline code). The implementing crates are
`yee-layout` (the closed-form coupled-microstrip model) and `yee-filter`
(the dimensional-synthesis inversion).

## 1. Where this sits in the flow

The synthesis page ends with a technology-agnostic coupled-resonator
network: a normalised coupling matrix `M`, a set of inter-resonator
coupling coefficients `k_{i,i+1}` at the design fractional bandwidth, and
input/output external quality factors `Qe_in` / `Qe_out`. Those numbers
describe *how strongly* adjacent resonators must couple and *how tightly*
the end resonators must couple to the feed — but say nothing about a piece
of metal on a substrate.

```text
synthesis  ─→  coupling matrix M, k_{i,i+1}, Qe   ──┐  (this page §3)
(filter-                                            │  invert the
 synthesis.md)                                       ▼  EM model
                       edge-coupled microstrip dimensions
                       (line width W, resonator length ℓ,
                        inter-resonator gaps s_i)
```

This page covers the **closed-form, narrowband** mapping (Filter Phase
F1.2.0, ADR-0097): an initial set of dimensions, derived by inverting two
already-validated quasi-static models, that seeds the later
EM-in-the-loop refinement (§5). The realisation target is the
**edge-coupled half-wave microstrip band-pass filter**: `N` parallel
half-wavelength resonators, each coupled to its neighbour through the
fringing field across a narrow gap (Hong & Lancaster, ch. 5; Pozar §8.7).
To dimension it we need two inverse maps — *target impedance* → *line
width*, and *target coupling* → *gap* — and the latter requires an
electrical model of a **pair** of coupled lines, which is §2.

## 2. Coupled microstrip even/odd modes

A single microstrip line is described by one characteristic impedance and
one effective permittivity. A *pair* of parallel strips supports two
quasi-TEM modes, and the standard decomposition (Pozar §8.7) splits any
excitation into them by symmetry:

- the **even mode** — both strips driven in phase. The symmetry plane
  between them behaves as a *magnetic wall* (open circuit). Field lines do
  not cross the plane; more of the field sits in the substrate, so the
  even-mode effective permittivity `εeff,e` is higher and the even-mode
  characteristic impedance `Z0e` is higher.
- the **odd mode** — the strips driven out of phase. The symmetry plane is
  an *electric wall* (short circuit / virtual ground). The field
  concentrates in the air-filled gap region, lowering the effective
  permittivity `εeff,o` and the characteristic impedance `Z0o`.

For physical strips `Z0e > Z0o > 0` always, and the split between them
widens as the gap closes — that split *is* the coupling. The
single-section coupler **voltage coupling coefficient** (Pozar §7.6) is

```text
k = ( Z0e − Z0o ) / ( Z0e + Z0o )
```

which tends to `0` for widely-spaced (weakly coupled) strips and grows
toward `1` as the gap closes. This is `yee_layout::coupling_coefficient`,
and it is the quantity §3 inverts.

### 2.1 The Kirschning-Jansen quasi-static model

`yee-layout` computes `Z0e`, `Z0o`, `εeff,e`, `εeff,o` for a symmetric
pair (two equal strips of width `W`, gap `S`, on a substrate of height `h`
and permittivity `εr`) with the **Kirschning-Jansen quasi-static**
closed-form model. The model takes the zero-frequency (DC) quasi-TEM
limit — full frequency dispersion is deliberately out of scope at this
fidelity — and builds on the single-line Hammerstad-Jensen `Z0(u)` and
`εeff(u)` forms (with `u = W/h` and `g = S/h`), where the even-mode
permittivity is evaluated at a coupling-modified width and the impedances
are corrected by the model's `Q1 … Q10` rational functions. Its published
accuracy is **better than ≈ 1.4 %** over `0.1 ≤ W/h ≤ 10`,
`0.1 ≤ S/h ≤ 10` against rigorous numerical reference data.

This is `yee_layout::coupled_microstrip`
(`crates/yee-layout/src/coupled.rs`), returning a `CoupledMicrostrip {
z0e_ohm, z0o_ohm, eps_eff_e, eps_eff_o }`. The implementation follows the
canonical transcription of the Kirschning-Jansen 1984 equation set in the
QUCS circuit simulator (`qucs-core` `mscoupled.cpp`, `analysQuasiStatic`,
"Kirschning" branch), reusing the crate's Hammerstad-Jensen single-line
helpers (`hj_eps_eff`, `hj_z0_air`) for self-consistency with the
single-line `microstrip_width` / `eps_eff` sizing used elsewhere in the
flow. It is pure `f64` and WASM-safe — part of the light client-side flow,
no FDTD, no native dependency (ADR-0089 / ADR-0094).

### 2.2 Validation

`coupled_microstrip` is gated (`coupled-001`,
`crates/yee-layout/tests/coupled_001_vs_published.rs`) against the worked
example in Steer, *Microwave and RF Design II: Transmission Lines*
(3rd ed.), §5.6, Example 5.6.1: an alumina substrate (`εr = 10`,
`W = h = 500 µm` so `W/h = 1`, `S = 250 µm` so `S/h = 0.5`) with published
`Z0e = 59 Ω`, `Z0o = 37 Ω`, `εeff,e = 7.28`, `εeff,o = 5.82`. The model
returns `Z0e ≈ 59.07 Ω`, `Z0o ≈ 36.96 Ω` — under 0.2 % error on the
impedances. A second gate (`coupled-002`) certifies the monotonicity
properties §3 relies on: `Z0e > Z0o > 0`, `k ∈ (0, 1)`, and `k` strictly
decreasing as the gap `S` grows.

## 3. From coupling coefficients to dimensions

With a coupled-line electrical model in hand, dimensional synthesis is
three inversions, each closed-form. This is
`yee_filter::dimension_edge_coupled`
(`crates/yee-filter/src/dimension.rs`), which takes a synthesized
`FilterProject` plus a `Substrate` and returns an
`EdgeCoupledDimensions { line_width_m, resonator_length_m, gaps_m,
target_k }`.

**Line width.** Every resonator and feed line is sized for the spec
reference impedance `Z0` (typically `50 Ω`) using the Hammerstad-Jensen
*synthesis* form (`yee_layout::microstrip_width`), which inverts
`Z0 → W` directly. The pair-coupling adjustment lives entirely in the gap,
not the width: both strips of a coupled section keep the single-line `Z0`
width.

**Resonator length.** Each resonator is a half guided wavelength at the
centre frequency `f0`:

```text
ℓ = λ_g / 2 = c / ( 2·f0·√εeff )
```

with `εeff` from `yee_layout::eps_eff` evaluated at the synthesized width
and `c = 299_792_458` m/s.

**Inter-resonator gaps.** This is the step that needs §2's coupled model.
For each adjacent resonator pair `(i, i+1)` the synthesis page's coupling
coefficient is the *target*:

```text
target_k[i] = FBW · M[i][i+1] = FBW / √( g_i · g_{i+1} ) = k_{i,i+1}
```

where the middle equality is exact because `yee-synth` builds the
normalised matrix as `M[i][i+1] = 1/√(g_i·g_{i+1})`, so multiplying its
off-diagonal by FBW reproduces the synthesized `k` vector verbatim. The
gap `s` that realises this target is the value for which the
coupled-section voltage coupling equals it:

```text
( Z0e(W, s) − Z0o(W, s) ) / ( Z0e(W, s) + Z0o(W, s) ) = target_k[i]
```

Because `k` is **strictly decreasing in the gap** `s` (the `coupled-002`
gate), this single-variable equation has a unique root, found exactly by
**bisection** over a manufacturable bracket (`5 µm` to `5 mm`): no
optimiser, no FDTD, no surrogate. The solver converges on the realised
coupling to a relative tolerance of `1e-4`; if a target falls outside the
achievable range over the bracket it is reported as a `GapNotBracketed`
error rather than silently clamped (`solve_gap`).

The mapping is **first-order and narrowband** — the standard
initial-dimensioning approximation. It assumes the per-section voltage
coupling of an isolated coupled pair equals the inter-resonator coupling
of the assembled filter, which is accurate for narrow fractional
bandwidths but neglects the loading each resonator imposes on its
neighbours and any frequency dispersion (the model is quasi-static). The
result is therefore an *initial estimate*, not a final geometry: it seeds
the EM-refinement loop of §5.

## 4. Worked example: 5-pole Chebyshev on FR-4

The committed fixture (gates `dim-001` / `dim-002`,
`crates/yee-filter/tests/`) is an order-`N = 5`, `0.5 dB`-ripple Chebyshev
band-pass filter at `f0 = 2 GHz` with fractional bandwidth `FBW = 0.10`
and `Z0 = 50 Ω`, realised on FR-4 (`εr = 4.4`, `h = 1.6 mm`).

Synthesis (the [filter-synthesis](./filter-synthesis.md) page) produces
the `N = 5` Chebyshev g-values
`[1, 1.706, 1.230, 2.541, 1.230, 1.706, 1]`, hence the four
inter-resonator coupling targets

```text
target_k = [ 0.069, 0.057, 0.057, 0.069 ]
```

The vector is mirror-symmetric, as expected for a synchronously-tuned
symmetric prototype: the two outer couplings are equal and stronger, the
two inner couplings equal and weaker.

Dimensional synthesis (§3) then gives:

| Quantity                                   | Value                                |
| ------------------------------------------ | ------------------------------------ |
| Line width `W` (HJ, `50 Ω`)                | ≈ `3.06 mm`                          |
| Effective permittivity `εeff`              | ≈ `3.33`                             |
| Resonator length `ℓ = λ_g/2` at `2 GHz`    | ≈ `41.1 mm`                          |
| Inter-resonator gaps `s_i`                 | ≈ `[ 2.80, 3.29, 3.29, 2.80 ] mm`   |

The gaps invert the couplings exactly the way §2.1 predicts: the strong
outer couplings (`k = 0.069`) need the *smaller* gaps (`2.80 mm`), the
weaker inner couplings (`k = 0.057`) the *larger* gaps (`3.29 mm`), and
the gap vector is mirror-symmetric because the coupling vector is. The
`dim-001` gate confirms the round-trip — feeding each solved gap back
through `coupled_microstrip` recovers its `target_k` to under 1 % — and
`dim-002` confirms the physical sanity (gaps positive, monotone in
coupling, width matches the HJ synthesis, length within ±2 % of
`λ_g/2`).

## 5. Limitations & next steps

This is the *closed-form, initial* half of dimensional synthesis. Several
pieces are deliberately deferred:

- **Coupler `k` vs resonator `k`.** `coupling_coefficient` returns the
  *coupler* voltage coupling `(Z0e − Z0o)/(Z0e + Z0o)`. The *resonator*
  coupling that an EM solve extracts from the two split resonant
  frequencies of a coupled pair is a related but distinct quantity,
  governed by the even/odd phase velocities (`εeff,e` / `εeff,o`); both
  even/odd sets are exposed so the FDTD coupled-resonator driver
  (F1.1b.1) can build its own resonator reference. For the narrowband
  first-order dimensioning here the coupler `k` is the right seed; the
  distinction matters once the EM loop closes.
- **External Q → feed geometry is not yet mapped.** The dimensioning here
  realises the *inter-resonator* couplings only. Turning `Qe_in` /
  `Qe_out` into an input/output feed or tap geometry is deferred to
  F1.2.1; the convenience layout builder (`dimension_edge_coupled_layout`)
  therefore uses a documented placeholder for the feed coupling rather
  than inventing a `Qe → gap` formula.
- **EM verification is the closing gate.** The narrowband approximation,
  the deferred feed coupling, and resonator mutual loading are all
  resolved by full-wave verification with the EM solver in the loop. On
  the planar track the back-end is **FDTD, not MoM** — the microstrip
  wave-port is ill-posed for planar MoM (its quasi-TEM mode's dominant
  field is substrate-normal `E_z`, orthogonal to the in-plane RWG
  surface-current basis; ADR-0064), whereas FDTD excites microstrip
  correctly and is broadband-validated. The coupled-resonator FDTD driver
  and `k`/`Qe` extraction (F1.1b.1), then surrogate-driven dimensional
  refinement (F1.2.1, using `yee-surrogate`), and finally an end-to-end
  published-filter gate (F1.3) build on the seed dimensions produced here.

## 6. References

- Kirschning, M., and Jansen, R. H. "Accurate Wide-Range Design Equations
  for the Frequency-Dependent Characteristic of Parallel Coupled
  Microstrip Lines." *IEEE Trans. Microwave Theory Tech.* 32.1 (Jan.
  1984), pp. 83–90 (corrected Nov. 1985). The even/odd quasi-static model
  `yee-layout::coupled_microstrip` implements; transcribed via the QUCS
  `qucs-core` `mscoupled.cpp` `analysQuasiStatic` "Kirschning" branch.
- Hong, J.-S., and Lancaster, M. J. *Microstrip Filters for RF/Microwave
  Applications.* Wiley, 2001. (Ch. 8 coupling coefficients and external
  Q; ch. 5 edge-coupled resonator filters — the dimensional-synthesis
  method of §3.)
- Pozar, D. M. *Microwave Engineering.* 4th ed. Wiley, 2012. (§8.7 coupled
  lines and the even/odd-mode decomposition; §7.6 the coupler coupling
  coefficient; §3.8 the Hammerstad-Jensen microstrip width / `εeff`.)
- Hammerstad, E., and Jensen, Ø. "Accurate Models for Microstrip
  Computer-Aided Design." *IEEE MTT-S Int. Microwave Symp. Digest*
  (1980), pp. 407–409. (The single-line `Z0` / `εeff` forms the
  Kirschning-Jansen coupled model and the line-width synthesis build on.)
- Steer, M. *Microwave and RF Design II: Transmission Lines.* 3rd ed.
  NC State University / LibreTexts. (§5.6, Example 5.6.1 — the
  coupled-microstrip even/odd worked example the `coupled-001` gate
  validates against.)
- `FILTER-DESIGN-ROADMAP.md` — the end-to-end filter-design plan (Filter
  Phases F1.1b / F1.2 / F1.3) this chapter tracks.
- ADR-0094 (`yee-layout::coupled_microstrip` even/odd model + coupler `k`),
  ADR-0097 (F1.2.0 closed-form edge-coupled dimensional synthesis).
