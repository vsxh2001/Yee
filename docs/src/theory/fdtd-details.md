# FDTD Details: CPML, NTFF, TF/SF, Lumped, Dispersive — Theory of Operation

This page is the theory-of-operation reference for the specialty
features built on top of the base FDTD walking skeleton documented in
[`fdtd.md`](./fdtd.md): the convolutional PML absorbing boundary, the
near-to-far-field surface transformation, the total-field /
scattered-field plane-wave source, lumped RLC ports, and dispersive
materials via the auxiliary differential equation (ADE). Same audience
as the planar-MoM and base-FDTD pages: an engineer reading source
with a textbook open. Equations are written in KaTeX so the inline
notation can stay close to the Rust source.

## 1. Introduction

The base FDTD chapter covers vacuum E/H leapfrog updates, the Courant
condition, and a Gaussian soft source on a single cell. That walking
skeleton is correct but inert: an unmodified Yee grid bounded by PEC
walls is only useful for closed-cavity eigenmode problems. Production
FDTD requires four orthogonal additions, each of which is a small
modification on top of the base update:

- A **broadband absorbing boundary** — the convolutional PML — so
  radiating geometry on a finite grid sees free space.
- A **near-to-far-field transformation** that turns the
  near-field history on an interior surface into the far-field
  radiation pattern.
- A **plane-wave source** that injects a coherent uniform-amplitude
  incident wave into a designated total-field region without
  polluting the scattered-field region outside it.
- **Lumped elements** that bridge the gap between the field solver
  and the standard circuit primitives (resistors, inductors,
  capacitors, series-RLC stubs).
- **Dispersive materials** whose permittivity depends on frequency —
  silicon at optical frequencies, gold below its plasma frequency,
  biological tissue across MHz to GHz — modelled by integrating
  auxiliary polarization variables alongside the fields.

Each of these features has shipped as its own Phase 2 sub-project
(`2.fdtd.1` CPML, `2.fdtd.2` NTFF, `2.fdtd.3` dispersive ADE,
`2.fdtd.5` TF/SF, `2.fdtd.6` lumped RLC, `2.fdtd.4` end-to-end
driver). This chapter is the consolidated derivation; the per-feature
unit tests in `crates/yee-fdtd/tests/` are the matching validation.

## 2. Convolutional PML (CPML)

The absorbing-boundary problem on a finite grid is to make outgoing
waves disappear without reflecting a measurable fraction of their
power back into the simulation domain. The modern answer is the
**convolutional perfectly matched layer (CPML)** of Roden and
Gedney (2000), which descends from Berenger's 1994 split-field PML
through a complex-frequency-shifted (CFS) reformulation of the
stretched-coordinate construction.

The starting point is a coordinate stretch in the frequency domain:
replace each spatial derivative inside the PML region with

$$
\frac{\partial}{\partial x} \;\to\; \frac{1}{s_x(\omega)} \frac{\partial}{\partial x},
\qquad
s_x(\omega) = \kappa_x + \frac{\sigma_x}{\alpha_x + j\omega\epsilon_0}.
$$

Outside the PML $(\kappa_x, \sigma_x, \alpha_x) = (1, 0, \text{free})$
recovers the standard derivative; inside the PML, $\sigma_x$ and
$\kappa_x$ ramp polynomially from the inner face to the outer face.
The CFS term $\alpha_x$ in the denominator absorbs the low-frequency
and evanescent components that Berenger's original PML handled
poorly — that's the only mathematical difference, and the only
reason CPML works on cavity-coupled and surface-wave problems where
the older PML failed.

The discrete time-domain implementation rewrites $1/s_x$ as the
Laplace transform of an exponential plus a delta function:

$$
\frac{1}{s_x(\omega)} = \frac{1}{\kappa_x}
+ \frac{-\sigma_x/(\kappa_x^2 \epsilon_0)}{\alpha_x/\epsilon_0 + \sigma_x/(\kappa_x \epsilon_0) + j\omega}.
$$

The second term, multiplied by a spatial derivative of a field, is
the Laplace transform of a convolution against the impulse response
$a_x \exp(-b_x t)$ for cached scalars $a_x, b_x$. Discretising on
the leapfrog grid turns the convolution into a one-tap recursive
update of auxiliary state $\psi_x$ at every PML cell:

$$
\psi_{E_y, x}^{n+1/2} = b_x\, \psi_{E_y, x}^{n-1/2}
+ a_x\, \frac{E_z^{n}(i+1,\cdot) - E_z^{n}(i,\cdot)}{\Delta x}.
$$

The corrected H or E update then adds $\psi$ to the standard
finite-difference curl term. The cost is one extra multiply-add per
cell per face per PML thickness — vastly cheaper than a fully
absorbing analytical boundary.

The grading profile inside the PML is polynomial. For a PML of
thickness $L$ cells and depth coordinate $\rho \in [0, L]$,

$$
\sigma(\rho) = \sigma_{\max} \left(\frac{\rho}{L}\right)^m,
\qquad
\kappa(\rho) = 1 + (\kappa_{\max} - 1) \left(\frac{\rho}{L}\right)^m,
$$

with order $m \in [3, 4]$ in practice. The maximum conductivity is
chosen so the theoretical reflection of a normal-incident wave hitting
a perfect PML is

$$
R(0) = \exp\!\left(-\frac{2\sigma_{\max} L \Delta x}{(m+1)\,c\,\epsilon_0\,\kappa_{\max}}\right) \approx 10^{-6},
$$

which solves for $\sigma_{\max} = -(m+1) c \epsilon_0 \ln R / (2 L \Delta x)$.
Yee uses $\kappa_{\max} = 1$ (no extra stretch beyond the absorbing
term), $m = 3$, $L = 8\text{–}12$ cells, and the $R(0) = 10^{-6}$
recipe. The CPML validation gate is **≥30 dB reflection reduction
versus PEC** for a plane wave at normal incidence, enforced by
`crates/yee-fdtd/tests/cpml_reflection.rs`.

## 3. Near-to-far-field transformation (NTFF)

FDTD only stores the near field. Antenna patterns live in the far
field, at distances large compared to $\lambda$ and to the radiator's
own extent. The standard bridge is the **Stratton-Chu surface
integral** built on the field equivalence theorem: enclose every
radiating source in a closed surface $S$ inside the computational
domain, and replace the exterior field with surface equivalent
currents that radiate the same far field.

Define $\hat{n}$ as the outward unit normal of $S$. The equivalence
theorem states that the surface currents

$$
\mathbf{J}_s(\mathbf{r}') = \hat{n}(\mathbf{r}') \times \mathbf{H}(\mathbf{r}'),
\qquad
\mathbf{M}_s(\mathbf{r}') = -\hat{n}(\mathbf{r}') \times \mathbf{E}(\mathbf{r}'),
\quad \mathbf{r}' \in S,
$$

radiating into free space, reproduce the exterior $(\mathbf{E},\mathbf{H})$
exactly. The far-field components of the radiated $\mathbf{E}$ at
observation angle $(\theta, \phi)$ and range $r$ are

$$
E_\theta(\theta, \phi) = -\frac{j k}{4\pi r} e^{-j k r}\,
\bigl(\eta_0\, N_\theta + L_\phi\bigr),
\qquad
E_\phi(\theta, \phi) = -\frac{j k}{4\pi r} e^{-j k r}\,
\bigl(\eta_0\, N_\phi - L_\theta\bigr),
$$

where $\eta_0 = \sqrt{\mu_0/\epsilon_0}$, $k = \omega/c$, and the
auxiliary potentials are the radiation integrals

$$
\mathbf{N}(\theta, \phi) = \int_S \mathbf{J}_s(\mathbf{r}')\, e^{j k\, \hat{r}\cdot\mathbf{r}'}\, dS',
\qquad
\mathbf{L}(\theta, \phi) = \int_S \mathbf{M}_s(\mathbf{r}')\, e^{j k\, \hat{r}\cdot\mathbf{r}'}\, dS'.
$$

The components $N_\theta, N_\phi, L_\theta, L_\phi$ are obtained by
projecting $\mathbf{N}, \mathbf{L}$ onto the spherical-coordinate
basis at the observation direction.

The discrete implementation accumulates a running discrete Fourier
transform at the desired probe frequency $f$. At every FDTD step,
for each Yee cell face on $S$, sample $\mathbf{E}$ and $\mathbf{H}$,
form $\mathbf{J}_s$ and $\mathbf{M}_s$ by the cross products above,
and accumulate

$$
\hat{\mathbf{J}}_s(\mathbf{r}') \mathrel{+}= \mathbf{J}_s^{n}(\mathbf{r}')\, e^{j 2\pi f n \Delta t} \Delta t,
$$

and similarly for $\hat{\mathbf{M}}_s$. After the time loop completes
the radiation integrals are evaluated by summing $\hat{\mathbf{J}}_s
e^{j k\hat r\cdot\mathbf{r}'} \Delta S$ over the surface. Because the
phase factor $e^{j k\hat r \cdot \mathbf{r}'}$ depends only on the
observation direction, sweeping $(\theta, \phi)$ is a cheap
post-processing operation; the expensive part is the on-the-fly DFT
accumulation during the run.

Yee places $S$ as an axis-aligned box one or more cells inside the
inner edge of the CPML. The implementation lives in
`crates/yee-fdtd/src/ntff.rs`; the validation gate is the short-dipole
radiation pattern — for a $z$-polarized current dipole short
compared to a wavelength the analytic far field is $|E_\theta|
\propto \sin\theta$, peaking at $\theta = 90^\circ$. The Track HH
end-to-end driver test recovers this peak with $|E_\theta(90^\circ)| =
1.000$ after normalization.

## 4. Total-field / scattered-field source

Plane-wave excitation on a finite grid is a deceptively hard problem.
A naïve approach — apply a planar array of soft sources across one
face — works in steady state but corrupts every measurement that
needs to separate incident from scattered field. The clean fix is
the **total-field / scattered-field (TF/SF)** decomposition.

Partition the grid into an interior axis-aligned box $V_\text{TF}$
and its exterior $V_\text{SF}$. Inside $V_\text{TF}$ the stored
$\mathbf{E}$, $\mathbf{H}$ represent the **total** field
$\mathbf{E}_\text{tot} = \mathbf{E}_\text{inc} + \mathbf{E}_\text{scat}$.
Outside $V_\text{TF}$ they represent only the **scattered** field
$\mathbf{E}_\text{scat}$. By linearity of Maxwell's equations both
fields satisfy the same update on either side; the only place the
decomposition needs accounting is the shared face between them, where
a curl stencil reads field values from both regions.

The fix is a one-cell correction on every face cell of $V_\text{TF}$.
For $+x$ propagation with $E_z$ polarization (the Phase 2.fdtd.5
supported case), at the front face $i = i_0$:

$$
H_y^{n+1/2}(i_0 - 1, j, k) \mathrel{-}= \frac{\Delta t}{\mu_0 \Delta x}\, E_z^{\,\text{inc}, n}(i_0),
$$

$$
E_z^{n+1}(i_0, j, k) \mathrel{-}= \frac{\Delta t}{\epsilon_0 \Delta x}\, H_y^{\,\text{inc}, n+1/2}(i_0 - 1/2),
$$

and similarly with a $+$ sign at the back face $i = i_1$. The two
corrections exactly cancel the spurious incident-field contribution
each curl stencil would otherwise pick up at the TF/SF interface.

The incident-field values $E_z^\text{inc}, H_y^\text{inc}$ at the
required locations are evaluated on a **1-D auxiliary grid** that
propagates the analytic plane wave using the same 1-D Yee update
the 3D scheme uses along its propagation axis. This automatically
matches the 1-D numerical dispersion of the 3D grid at normal
incidence, so the TF/SF correction stays consistent step-by-step
rather than slowly accumulating dispersion error.

Phase 2.fdtd.5 (Track OO) ships only the $+x$-propagating
$E_z$-polarized slab geometry — the TF region must span the full
transverse extent of the grid, with the transverse faces absorbed
by the outer CPML. Side-face corrections for a fully finite 3D TF
box are deferred. The validation gate is **≥30 dB total-field /
scattered-field contrast**; the shipped tests measure 68.5 dB,
comfortably above the gate.

## 5. Lumped R / L / C / series-RLC ports

The FDTD grid speaks fields; engineering specifications speak in
impedance, port voltages, and S-parameters. The connector is the
**lumped element**: a sub-cell modification to the standard E-update
at a single Yee cell that injects a current density consistent with
the element's voltage-current relationship, making the cell behave
as if a discrete circuit element bridged the two faces of an
$E_z$-edge.

For an element oriented along $\pm z$ at cell $(i, j, k)$, the
terminal voltage and current are

$$
V = E_z(i, j, k)\, \Delta z, \qquad I = \mathbf{J}_z \cdot \Delta x\,\Delta y,
$$

with the current convention pointing along $+z$. The lumped-element
constitutive law fixes $I$ as a function of $V$ (and history); the
resulting $\mathbf{J}_z$ enters Ampère's law as
$\epsilon_0 \partial E_z/\partial t = (\nabla \times \mathbf{H})_z - \mathbf{J}_z$.

**Pure resistor.** A series resistor with optional source EMF
$V_\text{src}(t)$ obeys $V_\text{term} = V_\text{src} + R\,I$, so
$\mathbf{J}_z = (E_z \Delta z - V_\text{src})/(R\,\Delta x\,\Delta y)$.
Substituting and using the semi-implicit average
$(E_z^n + E_z^{n+1})/2$ for the resistor current (Taflove & Hagness
§15.10) gives, after the standard Yee update has produced the
prediction $E_z^{n+1,\star}$,

$$
E_z^{n+1} = \frac{E_z^{n+1,\star} - \alpha\, E_z^{n} + \gamma\, V_\text{src}}{1 + \alpha},
\quad
\alpha = \frac{\Delta t\, \Delta z}{2\,\epsilon_0\,R\, \Delta x\,\Delta y},
\quad
\gamma = \frac{\Delta t}{\epsilon_0\, R\, \Delta x\,\Delta y}.
$$

This is the **validated path** in Phase 2.fdtd.6 — the energy
dissipation gate requires a fraction $\geq 0.3\%$ globally or local
absorption $> 5\times$ versus a passive grid, and Track WW measures
$> 30{,}000\times$ local absorption at the resistor cell.

**Series RLC.** The full series-RLC element keeps two extra state
variables: the inductor current $I_L^{n+1/2}$ at the half-step
(staggered with $E_z$) and the capacitor voltage $V_C^{n}$ at the
integer step. Kirchhoff's voltage law around the series string is

$$
V_\text{term} = V_\text{src} + R\, I_L + L \frac{dI_L}{dt} + V_C,
\qquad
C \frac{dV_C}{dt} = I_L,
$$

which discretises (forward Euler on $I_L$, integer-step accumulation
on $V_C$) into

$$
I_L^{n+1/2} = I_L^{n-1/2} + \frac{\Delta t}{L}\,\bigl(E_z^{n}\Delta z - R\,I_L^{n-1/2} - V_C^{n} - V_\text{src}^{n}\bigr),
$$

$$
E_z^{n+1} = E_z^{n+1,\star} - \frac{\Delta t}{\epsilon_0\, \Delta x\,\Delta y}\, I_L^{n+1/2},
\qquad
V_C^{n+1} = V_C^{n} + \frac{\Delta t}{C}\, I_L^{n+1/2}.
$$

Phase 2.fdtd.6 ships the series-RLC path as qualitative-only — it
compiles, runs without diverging on benign inputs, and reduces to the
pure resistor in the $L = 0$, $C \to \infty$ limit. Quantitative
validation against analytic series-RLC reflection is deferred to
Phase 2.fdtd.6.1.

**Source waveforms.** Two excitation envelopes are supported:
`HannSine` (a sinusoid at frequency $f$ ramped on by a raised-cosine
window over the first `ramp_steps` time steps) and `GaussianPulse`
(a Gaussian-modulated carrier at $f_0$ with spectral FWHM $b$).
Both are evaluated at the integer simulation time $t = n \Delta t$
inside `SourceWaveform::value`.

## 6. Dispersive materials via ADE

In real media — silicon, gold at optical frequencies, biological
tissue, loss-tangent dielectrics — the relative permittivity is a
function of frequency, and a constant $\epsilon_r$ is the wrong
model. The **Auxiliary Differential Equation (ADE)** approach
(Joseph & Hagness 1991; Luebbers & Hunsberger 1992 for the multi-pole
form) introduces one polarization vector $\mathbf{P}_k$ per pole and
integrates it alongside the fields on the same leapfrog stencil.
Ampère's law becomes

$$
\epsilon_0 \epsilon_\infty\, \frac{\partial \mathbf{E}}{\partial t} = \nabla \times \mathbf{H} - \sum_k \frac{\partial \mathbf{P}_k}{\partial t},
$$

where $\epsilon_\infty$ is the high-frequency limit of the relative
permittivity. The pole equation for $\mathbf{P}_k$ is one of:

**Drude pole** (Drude metals — gold, silver, aluminium below their
plasma frequency):

$$
\frac{\partial^2 \mathbf{P}}{\partial t^2} + \gamma\, \frac{\partial \mathbf{P}}{\partial t} = \epsilon_0\, \omega_p^2\, \mathbf{E},
$$

with plasma frequency $\omega_p$ and damping $\gamma$. The
corresponding susceptibility is
$\chi(\omega) = -\omega_p^2 / (\omega^2 + j\omega\gamma)$.

**Lorentz pole** (narrow-band resonances — bound electrons in
dielectrics, molecular vibrations):

$$
\frac{\partial^2 \mathbf{P}}{\partial t^2} + \gamma\, \frac{\partial \mathbf{P}}{\partial t} + \omega_0^2\, \mathbf{P} = \epsilon_0\, \Delta\epsilon\, \omega_0^2\, \mathbf{E},
$$

with resonant frequency $\omega_0$, damping $\gamma$, and oscillator
strength $\Delta\epsilon = \epsilon_s - \epsilon_\infty$ at the
pole. The corresponding susceptibility is
$\chi(\omega) = \Delta\epsilon\, \omega_0^2 / (\omega_0^2 - \omega^2 + j\omega\gamma)$.

**Debye pole** (orientational polarization — water, biological
tissues, lossy dielectrics):

$$
\tau\, \frac{\partial \mathbf{P}}{\partial t} + \mathbf{P} = \epsilon_0\, \Delta\epsilon\, \mathbf{E},
$$

with relaxation time $\tau$. The susceptibility is
$\chi(\omega) = \Delta\epsilon / (1 + j\omega\tau)$. Debye is a
first-order ODE — strictly cheaper per pole than Drude or Lorentz
but limited to broad, smooth dispersion. Composite materials are
modelled by summing multiple poles of any mix of the three types.

Discretisation samples $\mathbf{P}_k$ at integer steps (the same as
$\mathbf{E}$) so the polarization current $\partial\mathbf{P}_k/\partial t$
fits the centred difference $(\mathbf{P}_k^{n+1} - \mathbf{P}_k^{n-1})/(2\Delta t)$.
The Yee staggering on $\mathbf{H}$ at half-integer steps is
preserved; only $\mathbf{E}$ acquires the extra pole-current term in
its update. The implementation lives in
`crates/yee-fdtd/src/dispersive.rs`; the validation gate checks the
extracted complex permittivity against the closed-form $\chi(\omega)$
above by frequency-domain post-processing of an FDTD transmission
spectrum through a thin material slab.

## 7. Putting it together: end-to-end driver

The Phase 2.fdtd.4 driver wires CPML, a dipole current source, and
NTFF into a single public entry point. The dipole is a soft current
source on $E_z$ distributed over `dipole_length_cells` adjacent
cells along $z$, ramped on by a Hann window over the first three
periods so the source does not ring the grid. After the time loop
the driver sweeps $\theta \in [0^\circ, 180^\circ]$ in $5^\circ$
steps at $\phi = 0$, returning the normalized $|E_\theta|$ pattern.

```rust,ignore
use yee_fdtd::driver::{FdtdDriver, FdtdDriverConfig};
use yee_fdtd::grid::YeeGrid;

let grid = YeeGrid::vacuum(40, 40, 40, 5.0e-3);
let cfg = FdtdDriverConfig {
    n_steps: 5000,
    dipole_center_cells: (20, 20, 20),
    dipole_length_cells: 3,
    source_freq_hz: 1.0e9,
    ntff_surface_pad_cells: 2,
    cpml_thickness_cells: 8,
};
let pattern = FdtdDriver::new(grid, cfg).run();
// pattern.theta_deg[18] == 90.0, pattern.e_theta_phi0[18] == 1.0
```

A TF/SF plane-wave slab source is constructed against the same grid
and stepped each timestep against the auxiliary 1-D incident grid:

```rust,ignore
use yee_fdtd::sources::{PlaneWaveDirection, PlaneWaveSource};

let mut pw = PlaneWaveSource::new(
    /* i0 */ 8, /* i1 */ 32,
    /* j0 */ 0, /* j1 */ ny,
    /* k0 */ 0, /* k1 */ nz,
    PlaneWaveDirection::PlusX,
    /* frequency */ 1.0e9,
    /* ramp_steps */ 64,
    grid.dx, dt, /* pad */ 4,
);
// Each timestep:
pw.step_incident_h();
pw.correct_h(&mut grid);
pw.step_incident_e();
pw.correct_e(&mut grid);
```

A lumped resistor port is one constructor call plus a single
`correct_e` invocation per step, called **after** the standard
`update::update_e` so the cell already holds the leapfrog prediction
$E_z^{n+1,\star}$ when the correction overwrites it:

```rust,ignore
use yee_fdtd::lumped::{LumpedRlcPort, SourceWaveform};

let mut port = LumpedRlcPort::pure_resistor(
    /* cell */ (i, j, k),
    /* R (Ω) */ 50.0,
    SourceWaveform::HannSine {
        v0: 1.0, frequency: 1.0e9, ramp_steps: 64,
    },
);
// Each timestep, after update_e:
port.correct_e(&mut grid, n_step, dt);
```

The order discipline matters: TF/SF corrections wrap the curl
updates ($H$-correction after $H$-update, $E$-correction after
$E$-update), while the lumped-port correction runs strictly after
the $E$-update so it sees the unmodified prediction. Mixing the
order silently de-tunes both features.

## 8. What's not in this chapter

The features above cover the production needs of Phase 2 — antenna
patterns, scattering problems, dispersive media, lumped excitation
and termination. The following are *explicitly* out of scope and
tracked as future phases:

- **Subgridding** (Phase 2.fdtd.7 deferred). Locally refining the
  Yee lattice around small features without paying the global
  $O(N^4)$ cost of a finer uniform grid. Requires careful interface
  treatment to avoid late-time instability.
- **Conformal techniques** (Dey & Mittra 1997). Partially-filled-cell
  modifications to the curl integration weights that recover one to
  two orders of magnitude in geometric accuracy at fixed cell size,
  without the cost of subgridding. Yee currently staircases.
- **Multi-GPU domain decomposition.** One rank per GPU with NCCL
  halo exchange. Phase 4 work; the single-GPU memory wall is the
  trigger.
- **Piecewise Linear Recursive Convolution (PLRC)** for dispersive
  media. Mathematically equivalent to ADE for the standard pole
  models; we ship ADE only.
- **Hexagonal or unstructured-grid FDTD variants.** The cubic Yee
  lattice with its discrete de Rham complex is the only grid we
  support.
- **High-order schemes** such as FDTD(2,4). Trade a wider spatial
  stencil for fourth-order dispersion at the cost of a tighter
  Courant condition.
- **Time-reversal / adjoint FDTD** for sensitivity analysis. A
  natural fit for the Phase 3 surrogate-driven optimizer but
  deliberately deferred until the forward solver matures.

## 9. References

- Yee, K. S. "Numerical Solution of Initial Boundary Value Problems
  Involving Maxwell's Equations in Isotropic Media." *IEEE Trans.
  Antennas Propag.* 14.3 (1966), pp. 302–307. The original Yee
  algorithm; the staggering convention this chapter builds on. Also
  referenced in [`fdtd.md`](./fdtd.md).
- Roden, J. A., and Gedney, S. D. "Convolution PML (CPML): An
  Efficient FDTD Implementation of the CFS-PML for Arbitrary Media."
  *Microwave Opt. Technol. Lett.* 27.5 (2000), pp. 334–339. The
  CPML formulation as implemented in `crates/yee-fdtd/src/cpml.rs`.
- Taflove, A., and Hagness, S. C. *Computational Electrodynamics: The
  Finite-Difference Time-Domain Method.* 3rd ed., Artech House, 2005.
  §6 and §14 (TF/SF), §8 (NTFF), §15.10 (lumped elements), §9
  (dispersive media). The textbook reference for every feature in
  this chapter; equation numbering in the source is keyed to the
  3rd edition.
- Stratton, J. A. *Electromagnetic Theory.* McGraw-Hill, 1941. The
  surface-integral form of the radiation equivalence theorem that
  the NTFF transformation implements; §8.14 gives the Stratton-Chu
  expressions in their original notation.
- Joseph, R. M., Hagness, S. C., and Taflove, A. "Direct Time
  Integration of Maxwell's Equations in Linear Dispersive Media with
  Absorption for Scattering and Propagation of Femtosecond
  Electromagnetic Pulses." *Optics Letters* 16.18 (1991),
  pp. 1412–1414. The ADE formulation for Drude and Lorentz poles.
- Luebbers, R., and Hunsberger, F. "FDTD for Nth-Order Dispersive
  Media." *IEEE Trans. Antennas Propag.* 40.11 (1992), pp. 1297–1301.
  The multi-pole ADE generalization.
- Piket-May, M., Taflove, A., and Baron, J. "FDTD Modeling of Digital
  Signal Propagation in 3-D Circuits with Passive and Active Loads."
  *IEEE Trans. Microwave Theory Tech.* 42.8 (1994), pp. 1514–1523.
  The original lumped-element FDTD reference.
