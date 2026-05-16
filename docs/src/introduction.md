# Introduction

Yee is an open-source electromagnetic-simulation studio aimed at the RF,
microwave, and antenna-engineering workflows that today depend on closed-source
commercial tools. It is written in Rust, accelerated on NVIDIA CUDA, and
exposes a first-class Python interface so it slots into Jupyter without
ceremony. The project takes its name from Kane S. Yee, whose 1966 paper on
solving Maxwell's equations on a staggered grid still names the data structure
behind every modern FDTD solver.

The v1 beachhead is **planar PCB antennas and filters**, solved with a
GPU-accelerated planar Method of Moments kernel — the same regime served
commercially by Sonnet Suites and the Momentum component of Keysight ADS. A
second-track 3D FDTD solver handles the cases where planar MoM is the wrong
tool: radiation, transients, and dispersive materials.

## Why Yee exists

RF and microwave engineering is dominated by tools that cost tens to hundreds
of thousands of dollars per seat per year, ship under sales-gated commercial
licenses, and never expose their internals. The open-source landscape covers
parts of the space — openEMS and gprMax for FDTD, Palace for 3D FEM, NEC-2 for
wire-only MoM — but **no production-quality open-source planar MoM solver
exists**, and **no GPU planar MoM solver exists at all**. Planar MoM is also
the right tool for the bulk of working PCB-antenna and filter design, which is
exactly where most practicing RF engineers spend their day. That gap is the
wedge Yee aims at.

A second motivation runs underneath the first: every simulation Yee performs
is a labeled data point. The project treats machine-learning surrogate models
as first-class participants in the workflow, so that every parametric sweep
trains a model, every subsequent sweep gets cheaper, and natural-language
design — "give me a 2.4 GHz inset-fed patch on RO4003C" — becomes a tractable
front end rather than a magic trick.

## What this book covers

This book is the user-facing companion to the source tree. It will grow as the
solver does. For the multi-year plan, see
[`ROADMAP.md`](https://github.com/yee-em/yee/blob/main/ROADMAP.md) in the
repository root; for the dependency rationale, see
[`TECH_STACK.md`](https://github.com/yee-em/yee/blob/main/TECH_STACK.md); for
the active contribution guide, see
[`CONTRIBUTING.md`](https://github.com/yee-em/yee/blob/main/CONTRIBUTING.md).
