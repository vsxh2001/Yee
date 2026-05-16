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
