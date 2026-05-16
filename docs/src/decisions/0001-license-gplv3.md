# ADR-0001: Adopt GPL v3.0 or later as the project license

**Status:** Accepted
**Date:** 2026-05-16
**Deciders:** Yee maintainers

## Context

Yee is an open, GPU-accelerated electromagnetic simulator. The market it
enters is dominated by closed-source vendors selling perpetual or annual
node-locked seats at roughly USD 40,000 to USD 200,000 per seat per year:

- Dassault Systèmes SIMULIA CST Studio Suite (3D MoM / FDTD / FIT)
- ANSYS HFSS and HFSS-IE (FEM / surface integral equation)
- Sonnet em Professional (planar MoM)
- Keysight ADS Momentum and PathWave RF Synthesis (planar MoM)
- Cadence AWR AXIEM and Microwave Office (planar MoM, harmonic balance)

Open-source alternatives exist but are partial in scope:

- `openEMS` (GPL v3, FDTD only, octave/python frontend; no planar MoM,
  no surface MoM, no commercial layout import).
- `nec2c` and `xnec2` (NEC-2 derivatives, thin-wire MoM, GPL but stuck
  on a 1981 numerical formulation).
- `meep` (FDTD photonic, GPL v2+, optical regime, not RF-friendly UX).
- `scuff-em` (GPL v2, surface integral equation, research-grade UX).

Yee aspires to be a production-grade open alternative across the planar
MoM and FDTD problem space — usable from CLI, Python, and a desktop GUI,
with GPU acceleration, layout-format importers (KiCad, ODB++), and a
Touchstone-clean S-parameter export. That ambition is what makes the
license choice load-bearing rather than ceremonial.

Three license families were considered:

1. **Permissive** (MIT, Apache-2.0, BSD-3-Clause). Maximum adoption.
   Permits closed-source proprietary forks. Every dollar of value Yee
   creates can be captured and resold without sharing improvements.
2. **Weak copyleft** (LGPL, MPL-2.0). File-level or library-level
   copyleft. Permits dynamic linking from proprietary code while
   protecting Yee's own source. Comfortable middle ground but does not
   prevent a fork-and-extend strategy by a commercial vendor.
3. **Strong copyleft** (GPL v2, GPL v3, AGPL v3). Whole-program
   copyleft. Every distributed derivative must publish source under the
   same license. AGPL extends this to network use.

The transitive dependency graph (per `THIRD_PARTY_LICENSES.md`) imposes
real constraints on which of these are legally available:

- **Gmsh** — GPL v2 or later, with the standard linking exception for
  the API. Compatible with GPL v3 downstream.
- **OpenCASCADE (OCCT)** — LGPL 2.1 with the OCCT exception.
  Compatible with GPL v3.
- **NVIDIA CUDA Toolkit (cuSOLVER, cuBLAS, cuFFT)** — proprietary,
  dynamically linked, redistributed under NVIDIA's CUDA EULA.
  Compatible with GPL distribution via the system-library exception in
  GPL v3 §1.

Permissive licensing would let a future commercial vendor take Yee
verbatim, add a layout-import plug-in, and ship a closed product on top
— exactly the failure mode the project exists to avoid. AGPL v3 was
considered briefly but would entangle anyone exposing a Yee-powered
service over a network (including hobbyist web demos and CI test
servers) with source-disclosure obligations that the maintainers judged
to be a poor fit for a simulation library that runs primarily as a
local desktop / batch tool, not as a SaaS.

## Decision

Yee is licensed under the **GNU General Public License version 3.0 or
later (GPL-3.0-or-later)**, SPDX identifier `GPL-3.0-or-later`.

The root `LICENSE` file is the canonical GPL v3 text. Every source file
in the workspace carries an SPDX header. The `Cargo.toml` of every
workspace member declares `license = "GPL-3.0-or-later"`.

Third-party dependencies whose licenses are not GPL v3 compatible MUST
NOT be added to the dependency graph. Any new dependency added under
`crates/*/Cargo.toml` must be cross-checked against
`THIRD_PARTY_LICENSES.md` before merge.

## Consequences

**What becomes easier:**

- Contributions to Yee remain in the commons. Any organisation that
  forks Yee and distributes binaries (whether internally to a customer
  or publicly to the world) MUST release the corresponding source under
  GPL v3 or later. Improvements to the solver flow back.
- Yee can freely depend on other GPL v2-or-later and LGPL libraries,
  including Gmsh, OpenCASCADE, and (transitively, through `gmsh-sys`)
  the chain of mesh-handling code that already lives in the GPL galaxy.
- The license signals seriousness to academic and government users who
  need redistribution rights to reproduce results in published work.

**What becomes harder:**

- Closed-source proprietary extensions are blocked. A vendor cannot
  ship a GUI plug-in or layout-importer that links Yee statically (or
  in many readings, dynamically) without releasing the plug-in source.
- Some corporate users will not deploy GPL software internally even
  when no distribution occurs, citing legal-team policy that predates
  the GPL FAQ's mere-use clarifications. Yee will not pursue these
  users; the LGPL alternative was considered and rejected for §Context
  reasons.
- A dual-license offering (GPL plus a paid commercial license) would
  require Contributor License Agreements (CLAs) from every contributor.
  The maintainers have explicitly decided NOT to require CLAs, which
  forecloses future dual-licensing without re-soliciting consent from
  every author.

**What's now closed off:**

- Re-licensing to a permissive license retroactively would require
  consent from every contributor in `git log`. Practically irreversible.
- Adopting a Contributor License Agreement (CLA) is deferred until and
  unless the maintainers explicitly revisit this ADR.
- Any future component that wraps a proprietary closed-source library
  (e.g. a vendor-specific layout-format SDK) must be structured as a
  separately-distributed shim that the user installs themselves, not
  bundled into the GPL workspace, to avoid license incompatibility.

## References

- `LICENSE` (root of repository) — full GPL v3 text.
- `THIRD_PARTY_LICENSES.md` — audit of transitive dependency licenses.
- GNU GPL v3, Free Software Foundation, 29 June 2007:
  <https://www.gnu.org/licenses/gpl-3.0.html>
- SPDX License List, `GPL-3.0-or-later`:
  <https://spdx.org/licenses/GPL-3.0-or-later.html>
- GPL FAQ, "Is GPLv3 compatible with GPLv2?":
  <https://www.gnu.org/licenses/gpl-faq.html#v2v3Compatibility>
- NVIDIA CUDA Toolkit EULA, system-library exception path under GPL v3
  §1 "Corresponding Source ... does not include the operating system,
  or general-purpose tools".
