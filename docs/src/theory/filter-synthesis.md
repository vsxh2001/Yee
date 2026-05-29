# RF Filter Synthesis & the Yee Design Flow — Theory of Operation

This page is the theory-of-operation reference for Yee's RF
filter-synthesis stack, implemented in the `yee-synth`, `yee-filter`,
and `yee-layout` crates and exposed through the `yee filter synth`
CLI. Same audience as the planar-MoM and FDTD pages (an engineer
reading source code with a textbook open), same conventions
(plain-text math because the documentation build does not render
LaTeX, inline citations, source-file references in inline code).

## 1. Overview

A microwave filter is never designed directly in physical dimensions.
The classical flow goes through a sequence of abstractions, each one a
well-posed sub-problem with a published closed-form answer or a
published reference design to validate against:

```text
spec ─→ approximation ─→ lowpass     ─→ bandpass   ─→ coupling matrix
        (g-values)        prototype      transform     + external Q
                                                            │
        manufacturing  ←─ full-wave   ←─ physical    ←──────┘
        export            verify         layout
        (KiCad/Gerber,    (FDTD/FEM)      (dimensional
         STEP)                             synthesis)
```

The left-to-right top row is pure circuit math: a specification (centre
frequency, fractional bandwidth, passband ripple / return loss,
stopband rejection mask, reference impedance) is mapped to a normalised
**lowpass prototype** — a ladder of element values `g_k` — then to a
**coupling matrix** and a set of **external quality factors** that
describe an abstract coupled-resonator network. The bottom row is where
electromagnetics enters: those abstract couplings are realised as
**physical dimensions** (resonator lengths, coupling gaps, tap points),
laid out as meshable geometry, verified by a full-wave solve against the
spec mask, and exported to manufacturing files.

Yee delivers this flow as an interactive **desktop + web application**:
one `egui`/`eframe` codebase built for both native and WASM targets, with
the light synthesis/layout/plotting flow running client-side (it is pure
Rust with no native-only dependencies) and the heavy EM verification and
surrogate-driven dimensional synthesis running behind a native
`yee-server`. The strategic plan is `FILTER-DESIGN-ROADMAP.md`; the
application architecture is fixed in ADR-0089. The literature this stack
draws on is the standard microwave-filter canon: Matthaei, Young & Jones,
*Microwave Filters, Impedance-Matching Networks, and Coupling
Structures* (McGraw-Hill, 1964; reprinted Artech House, 1980); Pozar,
*Microwave Engineering*, 4th ed. (Wiley, 2012), ch. 8; and Hong &
Lancaster, *Microstrip Filters for RF/Microwave Applications* (Wiley,
2001).

The shipped pieces, one line each:

- **`yee-synth`** — pure synthesis math: lowpass-prototype g-values
  (Butterworth, Chebyshev), the lowpass→bandpass frequency transform,
  and the all-pole coupling-matrix + external-Q design. No EM, no I/O.
- **`yee-filter`** — the filter-domain data model (`FilterSpec`,
  `Prototype`, `CouplingMatrix`, `Topology`, `FilterProject`), the
  closed-form ideal response, and the `SpecMask` pass/fail gate.
- **`yee-layout`** — technology-specific parametric planar geometry:
  edge-coupled and hairpin microstrip generators with Hammerstad-Jensen
  width / ε_eff sizing (geometry only so far; the
  coupling-matrix→dimensions mapping is later-phase work).
- **`yee filter synth [--plot]`** — the CLI that reads a `FilterSpec`
  TOML, synthesizes the prototype and coupling matrix, sweeps the ideal
  response to a Touchstone file, grades it against the spec mask, and
  (with `--plot`) renders `|S21|` with the mask overlaid.

The current shipped scope is Filter Phase F0 → F1.0 (synthesis core +
parametric geometry); dimensional synthesis and full-wave verification
(§7) are the forward-looking F1.1+ roadmap.

## 2. Approximation: lowpass-prototype g-values

The first abstraction is the **lowpass prototype**: a normalised ladder
network (source impedance `g0`, reactive elements `g1 … gN`, load
termination `g_{N+1}`) with cutoff at `Ω = 1` rad/s, whose response
shape encodes the chosen approximation. An order-`N` prototype is the
vector `[g0, g1, …, gN, g_{N+1}]` of length `N + 2`. These are the
`yee_synth::prototype` outputs, validated by gate `synth-001` to ≤ 1e-6
against the published tables (Matthaei-Young-Jones Table 4.05-2;
Pozar Tables 8.3 / 8.4).

**Butterworth** (maximally flat) has the simple closed form

```text
g0 = 1
g_k = 2·sin( (2k − 1)·π / (2N) ),    k = 1 … N
g_{N+1} = 1
```

There is no passband ripple; the response is monotone and the band edge
is the 3 dB point.

**Chebyshev** (equi-ripple) trades a flat passband for a steeper skirt
at the price of `L_Ar` dB of in-band ripple. The element values follow
the standard recursion (Pozar §8.3, eq. 8.53):

```text
β   = ln( coth( L_Ar / 17.37 ) )
γ   = sinh( β / (2N) )
a_k = sin( (2k − 1)·π / (2N) ),    k = 1 … N
b_k = γ² + sin²( k·π / N ),        k = 1 … N
g1  = 2·a_1 / γ
g_k = 4·a_{k−1}·a_k / ( b_{k−1}·g_{k−1} ),    k = 2 … N
```

with the load termination

```text
g_{N+1} = 1                  (N odd)
g_{N+1} = coth²( β / 4 )      (N even)
```

The `17.37 = 40 / ln 10` constant is the standard Pozar form (it
converts the ripple from dB into the natural-log argument of `coth`).
The even-order load `g_{N+1} ≠ 1` is the well-known Chebyshev mismatch:
an even-order equi-ripple filter cannot be perfectly matched at both
ports with equal terminations, so the prototype carries an unequal load.
Both recipes live in `crates/yee-synth/src/lib.rs` (`butterworth` /
`chebyshev`).

A companion `min_order` estimates the smallest `N` meeting a required
stopband rejection `A_s` at a stopband ratio `Ω_s`: for Butterworth
`N ≥ log10((10^{A_s/10} − 1)/(10^{L_Ar/10} − 1)) / (2·log10 Ω_s)`, and
for Chebyshev `N ≥ acosh(√((10^{A_s/10} − 1)/(10^{L_Ar/10} − 1))) /
acosh(Ω_s)` (Pozar §8.3). Elliptic and Bessel prototypes are later-phase
work.

## 3. Lowpass → bandpass transform

A bandpass filter is obtained from the lowpass prototype by a frequency
substitution that maps the lowpass cutoff `Ω = ±1` onto the two
bandpass edges `ω1, ω2`. The centre frequency and fractional bandwidth
are the geometric-mean pair

```text
ω0  = √( ω1·ω2 )
FBW = ( ω2 − ω1 ) / ω0
```

and the prototype variable `Ω` is recovered from a real frequency `ω`
by the standard reactance map (Pozar §8.4):

```text
Ω = ( 1 / FBW )·( ω / ω0 − ω0 / ω )
```

`ω0` and `ω` may be in any consistent units (Hz or rad/s) because only
the ratio enters. `Ω = 0` is band centre, `Ω = ±1` are the band edges,
and `|Ω| > 1` is the stopband. This is `yee_synth::lowpass_to_bandpass`;
every downstream evaluation (ideal response, mask check, order estimate)
maps real frequencies through it before applying a lowpass formula.

## 4. Coupling matrix & external Q

The bandpass prototype is realised physically not as lumped L/C ladders
but as a chain of **coupled resonators**. The synthesis output the
physical design actually targets is therefore a set of **inter-resonator
coupling coefficients**, the **input/output external quality factors**,
and a **normalised coupling matrix** (Hong & Lancaster, ch. 8). From the
prototype g-values at fractional bandwidth `FBW`:

```text
k_{i,i+1} = FBW / √( g_i · g_{i+1} ),    i = 1 … N−1
Qe_in  = g0·g1 / FBW
Qe_out = g_N·g_{N+1} / FBW
```

The **normalised** (FBW-stripped) `N × N` coupling matrix `M` for a
synchronously-tuned all-pole filter has zero diagonal (every resonator
sits at the same centre frequency) and only nearest-neighbour
off-diagonal entries:

```text
M[i][i+1] = M[i+1][i] = 1 / √( g_i · g_{i+1} )
```

all other entries zero. The relation `k_{i,i+1} = FBW · M[i][i+1]` ties
the two: `M` is the bandwidth-independent topology, the `k` values are
its scaling at a given FBW. These are the `yee_synth::coupling_design`
outputs (`CouplingDesign { k, qe_in, qe_out, m }`), validated by gate
`synth-002` against a published coupled-resonator example. Cross-coupled
and elliptic matrices (with prescribed transmission zeros, via Cameron
synthesis) are later-phase work; the diagonal stays zero only for the
synchronous all-pole case.

## 5. Ideal response & the spec mask

With the order, approximation, and bandpass map fixed, the **ideal
response** is evaluated in closed form directly from the lowpass
transfer function, mapped through `Ω` of §3. The forward transmission
magnitude is

```text
Chebyshev:    |S21|² = 1 / ( 1 + ε²·T_N²(Ω) ),    ε = √( 10^{L_Ar/10} − 1 )
Butterworth:  |S21|² = 1 / ( 1 + Ω^{2N} )
```

where `T_N` is the Chebyshev polynomial of the first kind
(`T_N(x) = cos(N·acos x)` for `|x| ≤ 1`, `cosh(N·acosh|x|)` for
`|x| > 1`) and `ε` is the ripple constant. Reflection follows from
losslessness: `|S11|² = 1 − |S21|²`. This is `yee_filter::ideal_response`
(`crates/yee-filter/src/lib.rs`); it models magnitude only — the
zero-phase closed-form response is the *target* the dimensional synthesis
must reproduce, not yet a full-wave result. Driving S-parameters *from*
the coupling matrix (the Hong-Lancaster `[A] = [q] + pU − jM` admittance
form) is a later-phase realisation step.

The specification is graded against a `SpecMask`:

- `passband_ripple_db` — maximum allowed in-band insertion-loss ripple,
- `return_loss_db` — minimum required in-band return loss, and
- `stopband` — a list of `(frequency, minimum rejection dB)` points.

`yee_filter::check_mask` sweeps the response, classifies each frequency
as in-band (`|Ω| ≤ 1`) or stopband via the bandpass map, and returns a
`MaskReport` (overall pass/fail, worst-case passband ripple, worst-case
return loss, and per-stopband-point achieved-vs-required rejection).
This is gate `filt-001`: the synthesized Chebyshev response must meet its
own ripple / return-loss / rejection mask. The `yee filter synth --plot`
CLI renders the swept `|S21|` with the mask overlaid (the spec→visual
pipe), via the `yee-plotters` `draw_sparam_with_mask` overlay.

## 6. From abstract circuit to geometry

The coupling matrix is technology-agnostic; turning it into a
manufacturable structure is technology-specific. For the planar
(microstrip / stripline) track, `yee-layout` generates parametric
geometry from explicit physical dimensions: an **edge-coupled** BPF
(`N` parallel half-wavelength resonators plus end feed lines) and a
**hairpin** BPF (`N` U-folded resonators plus a tapped feed). Line
sizing uses the Hammerstad-Jensen closed form for microstrip width and
effective permittivity (Hammerstad & Jensen, "Accurate Models for
Microstrip Computer-Aided Design," *IEEE MTT-S Digest*, 1980; Pozar
§3.8), in `yee_layout::microstrip_width` / `eps_eff`, validated by gates
`geo-001/002/003`. The output is a top-metal-on-substrate footprint
(`Vec<Polygon>` of traces plus port references and a bounding box) with a
dependency-free SVG preview. This is the **dimensions → geometry**
direction only; the **coupling-matrix → dimensions** mapping is the
dimensional-synthesis step of §7.

## 7. Dimensional synthesis & full-wave verification (forward-looking)

Closing the loop — turning a target coupling matrix into the gaps and
lengths that actually realise it, and proving the realised structure
meets the spec — is the F1.1+ roadmap and not yet shipped. The plan:

- **Coupling / Qe extraction.** Drive a coupled resonator pair and a
  singly-loaded resonator through the EM solver and extract the realised
  `k` and `Qe` from the response (the two split resonances for `k`, the
  loaded-resonator group delay or bandwidth for `Qe`).
- **Surrogate-BO with the EM solver in the loop.** A single full-wave
  filter solve is seconds-to-minutes and dimensional synthesis needs
  many, so the optimisation is driven by `yee-surrogate` (Gaussian-
  process surrogate + Bayesian optimisation + active learning), never raw
  grid search: the surrogate proposes candidate dimensions, the EM solver
  evaluates the few that matter, and the loop converges the extracted
  couplings onto the §4 targets.
- **The planar EM back-end is FDTD, not MoM.** The microstrip wave-port
  is ill-posed for planar MoM — the quasi-TEM mode's dominant field is
  substrate-normal `E_z`, orthogonal to the in-plane RWG surface-current
  basis, so MoM microstrip S-parameters are port-limited (ADR-0064).
  FDTD excites microstrip correctly and is broadband-validated (CPML /
  NTFF / dispersive / lumped-port / skin-depth gates), so planar filter
  verification and dimensional synthesis run on FDTD. Waveguide / cavity
  filters use the FEM solver (gated on its wave-port maturation).

The headline end-to-end gate is reproducing a published 5-pole hairpin
bandpass filter — spec → synthesis → layout → FDTD S-parameters within a
±1 dB tolerance of the reference — and emitting Gerber that re-imports to
matching geometry. Manufacturing export (KiCad/Gerber for planar and
lumped, STEP for waveguide) is the final stage.

## 8. References

- Matthaei, G. L., Young, L., and Jones, E. M. T. *Microwave Filters,
  Impedance-Matching Networks, and Coupling Structures.* McGraw-Hill,
  1964; reprinted Artech House, 1980. (Table 4.05-2 prototype g-values;
  the `synth-001` reference.)
- Pozar, D. M. *Microwave Engineering.* 4th ed. Wiley, 2012. (§8.3
  prototype g-values and the Chebyshev recursion eq. 8.53; §8.4 the
  lowpass→bandpass transform; §3.8 microstrip design equations.)
- Hong, J.-S., and Lancaster, M. J. *Microstrip Filters for RF/Microwave
  Applications.* Wiley, 2001. (Ch. 8 coupling coefficients and external
  Q; chs. 5–6 edge-coupled and hairpin microstrip filters.)
- Cameron, R. J. "Advanced Coupling Matrix Synthesis Techniques for
  Microwave Filters." *IEEE Trans. Microwave Theory Tech.* 51.1 (2003),
  pp. 1–10. (Cross-coupled / elliptic synthesis; later-phase work.)
- Hammerstad, E., and Jensen, Ø. "Accurate Models for Microstrip
  Computer-Aided Design." *IEEE MTT-S Int. Microwave Symp. Digest*
  (1980), pp. 407–409. (Microstrip width / ε_eff used by `yee-layout`.)
- `FILTER-DESIGN-ROADMAP.md` — the end-to-end filter-design plan
  (stages, phases F0–F6, validation gates) this chapter tracks.
- ADR-0084 (synthesis core), ADR-0086 (`yee-layout` parametric geometry),
  ADR-0088 (`yee filter synth --plot`), ADR-0089 (desktop + web app
  architecture).
