# Third-Party Licenses

Yee itself is distributed under the **GNU General Public License v3.0 or later**
(see [`LICENSE`](LICENSE)). This file documents every third-party component
linked, bundled, or invoked by the Yee codebase, so that downstream redistributors
can comply with each component's terms. Every entry below cites the **feature
flag** and **crate** that pulls the component in, so it is unambiguous which
build configurations carry which obligations.

If you ship a binary built from this tree, you are responsible for honoring all
applicable terms below.

---

## Gmsh

- **License:** GNU General Public License v2 or later, **with the Gmsh linking
  exception** that explicitly permits combining Gmsh with software under any
  GPL-compatible license (including GPL v3) without making the combined work a
  derivative work under section 5 of GPL v2.
- **Authoritative reference:** <https://gmsh.info/doc/texinfo/gmsh.html#License-and-credits>
- **Where it appears in Yee:** the `yee-mesh` crate, behind the **`gmsh`**
  cargo feature. When that feature is enabled, `crates/yee-mesh/build.rs` uses
  `bindgen` to generate Rust FFI bindings against `$GMSH_SDK_ROOT/include/gmshc.h`,
  and the resulting library dynamically links the Gmsh shared object at runtime.
- **Obligations:** if you redistribute a binary that statically or dynamically
  links Gmsh, the linked-exception text plus the upstream license must travel
  with the binary. The combination is permitted; the bookkeeping is not optional.

## OpenCASCADE Technology (OCCT)

- **License:** GNU Lesser General Public License v2.1, with the OCCT exception
  (an additional permission specific to OCCT's use of certain headers).
- **Authoritative reference:** <https://dev.opencascade.org/resources/licensing>
- **Where it appears in Yee:** the `yee-io` crate, behind the **`opencascade`**
  cargo feature. It is consumed through the workspace dependency
  `opencascade = "0.2"` (a permissive Rust binding crate) plus the LGPL-licensed
  OCCT C++ libraries that the bindings load at runtime.
- **Obligations:** under LGPL v2.1, downstream binaries that link OCCT must
  allow end users to replace the OCCT libraries with their own modified
  version. Dynamic linking (the default with the upstream `opencascade-sys`
  bindings) is the path of least friction; static linking is permitted but
  requires shipping the relevant object files or relinkable archive.

## NVIDIA CUDA Libraries

- **License:** proprietary (NVIDIA Software License Agreement; NVIDIA CUDA EULA).
  Not open source.
- **Authoritative reference:** <https://docs.nvidia.com/cuda/eula/index.html>
- **Where it appears in Yee:** the `yee-cuda` crate, behind the **`cuda`** cargo
  feature, via the workspace dependency
  `cudarc` with the **`dynamic-loading`** feature enabled. CUDA libraries are
  **never statically linked** into the Yee binaries; they are resolved at
  process start by `dlopen` against the user's locally installed CUDA Toolkit.
  The Yee project does **not** redistribute any NVIDIA binary.
- **Affected libraries** (all loaded dynamically at runtime through `cudarc`'s
  feature surface):
  - cuBLAS
  - cuBLASLt
  - cuSOLVER
  - cuSPARSE
  - cuFFT
  - cuRAND
  - NCCL (NVIDIA Collective Communications Library)
- **Obligations:** because Yee does not redistribute NVIDIA binaries, the CUDA
  EULA does not flow through to Yee's source distribution. End users who
  install the CUDA Toolkit accept the EULA from NVIDIA directly. If you create
  a downstream redistribution that *does* bundle the CUDA runtime, you must
  ship the EULA alongside it and comply with its redistribution clauses.

---

## Permissive Rust Dependencies

The table below covers every direct workspace dependency declared in the top-level
`Cargo.toml` under `[workspace.dependencies]`. License strings are reproduced from
the upstream crate's `Cargo.toml` `license` field as of this writing. Where a
crate is dual-licensed (`A OR B`), Yee links it under whichever of the two we
prefer per the project's GPL-v3 posture; either license remains a valid
recipient choice.

| Crate | License (SPDX) | Upstream |
| --- | --- | --- |
| `cudarc` | MIT OR Apache-2.0 | <https://crates.io/crates/cudarc> |
| `faer` | MIT | <https://crates.io/crates/faer> |
| `nalgebra` | Apache-2.0 | <https://crates.io/crates/nalgebra> |
| `ndarray` | MIT OR Apache-2.0 | <https://crates.io/crates/ndarray> |
| `sprs` | Apache-2.0 | <https://crates.io/crates/sprs> |
| `num-complex` | MIT OR Apache-2.0 | <https://crates.io/crates/num-complex> |
| `glam` | MIT OR Apache-2.0 | <https://crates.io/crates/glam> |
| `parry3d` | Apache-2.0 | <https://crates.io/crates/parry3d> |
| `geo` | MIT OR Apache-2.0 | <https://crates.io/crates/geo> |
| `i_overlay` | MIT | <https://crates.io/crates/i_overlay> |
| `opencascade` | LGPL-2.1 | <https://crates.io/crates/opencascade> |
| `opencascade-sys` | LGPL-2.1 | <https://crates.io/crates/opencascade-sys> |
| `egui` | MIT OR Apache-2.0 | <https://crates.io/crates/egui> |
| `eframe` | MIT OR Apache-2.0 | <https://crates.io/crates/eframe> |
| `egui_plot` | MIT OR Apache-2.0 | <https://crates.io/crates/egui_plot> |
| `egui_dock` | MIT OR Apache-2.0 | <https://crates.io/crates/egui_dock> |
| `wgpu` | MIT OR Apache-2.0 | <https://crates.io/crates/wgpu> |
| `rerun` | Apache-2.0 OR MIT | <https://crates.io/crates/rerun> |
| `plotters` | MIT | <https://crates.io/crates/plotters> |
| `pyo3` | MIT OR Apache-2.0 | <https://crates.io/crates/pyo3> |
| `numpy` | BSD-2-Clause | <https://crates.io/crates/numpy> |
| `tracing` | MIT | <https://crates.io/crates/tracing> |
| `tracing-subscriber` | MIT | <https://crates.io/crates/tracing-subscriber> |
| `thiserror` | MIT OR Apache-2.0 | <https://crates.io/crates/thiserror> |
| `anyhow` | MIT OR Apache-2.0 | <https://crates.io/crates/anyhow> |
| `clap` | MIT OR Apache-2.0 | <https://crates.io/crates/clap> |

Transitive dependencies are not listed here; they inherit the same posture and
their licenses can be enumerated with `cargo deny check licenses` or
`cargo about generate` once those tools are wired into CI in a later phase.
