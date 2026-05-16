# Yee Technology Stack

Selections and rationale, current as of **early 2026**. Versions are the recommended floor; pin in `Cargo.toml` and `pyproject.toml` and migrate deliberately.

## Languages and toolchain

| Concern | Choice | Version floor | Rationale |
|---|---|---|---|
| Systems language | **Rust** | 1.85+ | Driven by `maturin` 1.10 (Rust 2024 edition) and `pyo3` 0.28 (MSRV 1.83). Fearless concurrency, zero-cost abstractions, mature CUDA story. |
| Python | **CPython** | 3.10+ | Most RF engineers use Jupyter; 3.10 is the comfortable floor for modern scientific Python in 2026. |
| CUDA Toolkit | **CUDA** | 12.4+ (13.0 supported) | Tensor-Core / Hopper / Blackwell coverage; matches `cudarc` 0.19 support matrix. |
| Build | Cargo workspace + `maturin` | Cargo 1.85, maturin 1.10 | Single source of truth; wheels via `maturin build --release`. |

---

## CUDA host bindings — `cudarc`

**Choice: `cudarc` 0.19+.**

We evaluated three candidates:

1. **`cudarc`** (coreylowman) — latest 0.19.4, ~monthly releases, **the production-grade Rust CUDA binding in 2026**. Used by Hugging Face Candle, mistral.rs, candle-vllm. Covers CUDA driver, NVRTC, **cuBLAS, cuBLASLt, cuSOLVER, cuSPARSE, cuFFT, cuRAND, cuTENSOR, NCCL, cuDNN** out of the box. Supports CUDA 11.4–13.0 selectable at build time. Three linking modes: dynamic-loading (default, no libs at build time), dynamic-linking, static-linking. License: MIT/Apache-2.0.
2. **`cust`** (part of Rust-GPU/Rust-CUDA project) — last crates.io release 0.3.2 is years old. The project was rebooted in January 2025 and is genuinely active again, but as of late 2025 the rebooted code is **git-only** — no crates.io publication. Not yet suitable for a production project.
3. **`rustacuda`** — last release 0.1.3 (2021). **Unmaintained.**

We will **standardize on `cudarc` as the single CUDA host binding for the entire project.** Where `cudarc` lacks a top-level safe wrapper (notably cuSPARSE), we write a thin internal `yee_cuda::sparse` module that wraps `cudarc::cusparse::result::*` with RAII handles. Caveat: `cudarc` self-describes as "pre-alpha" despite massive production use; we pin to exact minor version and introduce a thin internal abstraction layer so we can swap if needed.

### Kernel authoring

We will **write CUDA kernels in CUDA-C and JIT them via NVRTC** (the approach Candle uses). This is the path of least resistance in 2026.

We considered Rust-on-GPU via `rustc_codegen_nvvm` (the Rust-CUDA project's NVVM backend). As of August 2025 its toolchain was bumped from a 3-year-old nightly to `nightly-2025-06-23`, CUDA 12.x is "experimental," and CI was restored — meaningful progress, but it still depends on a pinned nightly and is not production-ready. **We will re-evaluate in late 2026.**

```toml
cudarc = { version = "0.19", default-features = false, features = [
    "std", "driver", "nvrtc",
    "cublas", "cublaslt", "cusolver", "cusparse", "cufft", "curand",
    "f16",
    "cuda-version-from-build-system",
    "dynamic-loading",
] }
```

---

## Linear algebra

| Workload | Choice | Version | Why |
|---|---|---|---|
| Small static matrices (geometry, transforms) | **`nalgebra`** | 0.34+ | The de-facto standard. ~4.3M downloads/month. Compile-time-sized types with zero heap alloc. Pair with `nalgebra-sparse` if useful. |
| Medium/large dense CPU LU/Cholesky/SVD | **`faer`** | 0.23+ | Pure-Rust, MIT, competitive with Intel MKL on GEMM and factorizations. Supports complex types via `num_complex::Complex64`. Used as host-side fallback for MoM and reference implementation. |
| Dense LU at MoM scale (10k–100k complex) on GPU | **cuSOLVER via cudarc** | — | `cusolverDnZgetrf`/`Zgetrs`. cuSOLVERMg for multi-GPU. This is the hot path. |
| N-d arrays at the Python boundary | **`ndarray` + `numpy` crate** | ndarray 0.17, numpy 0.28 | Zero-copy NumPy ↔ `ndarray::ArrayView`. Use for FDTD field arrays exposed to notebooks. |
| Sparse SpMV on CPU (FEM/FDTD operators) | **`sprs`** or `faer` sparse | sprs 0.11+ | `faer`'s sparse module is the best pure-Rust sparse direct solver in 2026; `sprs` is the workhorse for SpMV/SpGEMM. |
| Sparse SpMV / iterative solves on GPU | **cuSPARSE via cudarc** | — | We wrap the result/sys layer ourselves into a small safe API. |

### Why not just `nalgebra` for everything?

`nalgebra` is excellent for ergonomic, small, statically-sized math but is explicitly **not optimized for medium/large dense factorization**. `faer`'s own docs recommend `nalgebra` for low-dimensional/graphics work and `faer` for medium/large dense. We use both for what each does best.

### Scaling caveat

10k×10k complex double = 1.6 GB — feasible on a single GPU via cuSOLVER. **100k×100k dense is 1.6 TB** — not feasible by brute force. For n ≥ 50k we will introduce iterative GMRES with block-diagonal preconditioning in Phase 1; **MLFMA / ACA / H-matrix compression is Phase 4 work.**

---

## CUDA numerical libraries (wrapped via cudarc)

| Library | Use in Yee | `cudarc` coverage |
|---|---|---|
| **cuBLAS** | Dense complex GEMM/GEMV for MoM matrix-fill aggregation | ✅ Full safe wrapper |
| **cuBLASLt** | Mixed-precision / Tensor-Core GEMM | ✅ Full safe wrapper |
| **cuSOLVER** | Dense LU (`Dn{C,Z}getrf`), Cholesky, QR, SVD, eigen | ✅ Safe wrapper (`cudarc::cusolver::safe`) |
| **cuSOLVERMg** | Multi-GPU dense factorization | ✅ |
| **cuSPARSE** | SpMV / SpMM for FEM/FDTD iterative solvers | ⚠️ `result` + `sys` only; we wrap |
| **cuFFT** | Near-to-far-field projection, FDTD post-processing | ✅ Safe wrapper (dynamic-loading only) |
| **cuRAND** | Random sources (rough surfaces, Monte Carlo) | ✅ |
| **NCCL** | Multi-GPU collective comms for FDTD domain decomp | ✅ Safe wrapper |

No other actively-maintained, high-coverage Rust wrapper exists for these libraries; standardizing on `cudarc` is the only sane choice.

---

## Meshing — Gmsh

**Choice: Gmsh 4.15+ via in-tree `bindgen` FFI to `gmshc.h`.**

Gmsh is the de-facto open-source mesher: actively maintained (4.15.2 stable as of March 2026), licensed **GPL v2+ with a linking exception** that is explicitly compatible with GPL v3 hosts. Its C API (`gmshc.h`) has been stable across the 4.x line since 2018.

We considered the existing Rust bindings:
- **`gmsh-sys`** and **`rgmsh`** (mxxo) — last updated November 2019, target Gmsh 4.4.1, six years stale. **Unmaintained.**

**We regenerate bindings ourselves with `bindgen` against the current `gmshc.h` from the official SDK** and write a thin safe wrapper in-tree as `yee-mesh`. The C API surface is a few hundred functions; this is roughly one day of work and lets us pin the Gmsh version we link against.

Fallback for cases where the Rust binding lags: shell out to the Gmsh CLI or call the Python API via PyO3 — slower iteration but bulletproof.

---

## CAD kernel — OpenCascade

**Choice: `opencascade-rs` 0.2+ wrapping OCCT 7.9.x (with a watch on 8.0, due May 2026).**

OCCT is **LGPL 2.1 with an exception**, confirmed GPL-2+/GPL-3 compatible by the Open CASCADE FAQ. It is the only mature, open-source B-rep CAD kernel; FreeCAD and Gmsh's OCC backend both use it.

Rust bindings:
- **`opencascade-rs`** (Bschwind) — 0.2.0, **actively maintained** (last commit February 2026), uses `cxx.rs` bridging. Exposes primitives, sketches, extrude/revolve/sweep, fillet/chamfer, booleans, STEP/STL/DXF/SVG/KiCad I/O. Builds OCCT as a git submodule — expect 5–15 minute cold builds; cache aggressively in CI.
- **`truck`** (RICOS) — pure-Rust B-rep with NURBS, MIT-licensed, **promising but immature**. Lacks robust booleans for complex topologies and mature real-world STEP import. **We watch but do not bet on it.**

We use `opencascade-rs` as primary. For cases it does not expose, we drop into raw FFI via `opencascade-sys`. For STEP/IGES *import only*, the Gmsh OCC factory (`importShapes("file.step")`) is an excellent fallback that reuses Gmsh's OCCT linkage and gives us tessellated geometry directly.

---

## Geometry primitives and spatial queries

| Concern | Crate | Version | Notes |
|---|---|---|---|
| Simulation math (vectors, matrices, transforms) | **`nalgebra`** | 0.34+ | Already chosen above. |
| Graphics math (shader-side, wgpu interop) | **`glam`** | 0.30+ | What wgpu/egui/rerun use natively. |
| Spatial queries / BVH / ray casts | **`parry3d`** | latest (mid-2025+) | Apache-2.0. BVH (`Qbvh`), TriMesh, AABB, distance queries — perfect for "is this field point inside copper?" and far-field projection ray casts. |
| 2D PCB polygon ops (copper layers, Gerber, DXF) | **`geo`** + **`i_overlay`** | geo 0.30, i_overlay 4.x | `i_overlay` is faster and more robust than `geo-booleanop` for PCB-grade Boolean operations. |
| Half-edge / triangle mesh data structures | In-tree module | — | No clear community standard; we adopt half-edge structures from `truck-topology` or roll our own thin wrapper around `parry3d::shape::TriMesh`. |

---

## Python bindings

**Choice: `pyo3` 0.28+ with `abi3-py310`, built by `maturin` 1.10+, distributed as `manylinux_2_28` wheels.**

- **PyO3 0.28.3** — supports CPython 3.7–3.14 and PyPy/GraalPy. The `Bound<'py, T>` API has been stable since 0.21. **Free-threaded Python (PEP 703 / 3.13t / 3.14t) is supported** since 0.23, mature in 0.28; note that `abi3` does not yet apply to the free-threaded ABI (PEP 803 is in flight, targeting Python 3.15). MSRV is Rust 1.83.
- **`abi3-py310`** lets us ship one wheel per OS/arch that works on Python 3.10+, vastly simplifying distribution.
- **`maturin` 1.10.2** (Nov 2025) — MSRV Rust 1.85; uses Rust 2024 edition internally. Builds wheels for Win/Linux/macOS/FreeBSD plus iOS-simulator targets in recent versions. Manylinux tagged automatically; we use `manylinux_2_28` (glibc 2.28+ / RHEL 8 baseline) for 2026.
- **`numpy` crate** matched to PyO3 0.28 for zero-copy NumPy ↔ `ndarray` interop. Essential for FDTD field arrays in Jupyter.
- **`pyo3-async-runtimes`** when we expose async solver callbacks.

CI builds via `PyO3/maturin-action@v1` with `manylinux: "2_28"`; we publish with `uv publish`.

---

## GUI — `egui` + `eframe`

**Choice: `egui` 0.34+ with `eframe`, a custom `wgpu` 3D viewport painted into an egui rect.**

We evaluated five frameworks:

| Framework | Latest | Verdict |
|---|---|---|
| **`egui`** | 0.34.2+ | **Selected.** Immediate-mode, pure-Rust, wgpu/glow backends. Sponsored by Rerun. Used by Rerun, MakerPnP, many scientific tools. |
| Tauri | 2.9.6 | Forces a JS/TS frontend; awkward for a tightly-coupled solver studio with a custom GPU viewport. |
| Dioxus | 0.7.0 | Promising (especially Dioxus Native, no webview), but the native renderer is brand new. |
| iced | 0.14 | Excellent architecture, but missing mature plotting/table widgets for scientific use. Revisit after 1.0. |
| Slint | 1.16+ | Beautiful for embedded HMIs, but the `.slint` DSL and absent plotting widget make it a poor fit. |

**Why egui wins for an EM studio:**
1. **3D viewport integration.** `egui-wgpu` lets you render into a texture and embed it in any egui rect — this is exactly how Rerun's 3D viewer is built.
2. **Plotting.** `egui_plot` drops in for S-parameter, Smith-chart, time-domain, and convergence plots.
3. **Docking and panels.** `egui_dock` or `egui_tiles` (Rerun's panel system) cover tabbed tools.
4. **Precedent.** Rerun, Bevy's editor work, several open-source FEA/CFD tools — the "scientific Rust GUI" segment is overwhelmingly egui-based today.
5. **Cross-platform and web.** `eframe` runs on Win/Mac/Linux natively and compiles to WASM unchanged.

**Honest downsides:** breaking minor releases roughly quarterly (mitigated by isolating UI behind a stable internal trait); not "native-looking" (acceptable for a technical tool).

---

## 3D visualization

**Choice: roll the viewport on `wgpu` directly, embedded inside an `egui` paint callback. Use `rerun` as a complementary debugging logger.**

| Library | Latest | Role for Yee |
|---|---|---|
| **`wgpu`** | 26.x | Primary viewport; custom shaders for surface heatmaps, far-field deformed spheres, cut planes. |
| **`rerun`** | 0.31–0.32 | Optional structured logging sink for solver internals: dump E-field samples per frequency, mesh per refinement step, far-field per iteration. Engineers run the Rerun viewer alongside the simulator. Apache-2.0. |
| `bevy` | 0.17–0.18 | **Not used.** Game engine ECS is the wrong shape for a deterministic, parameter-driven scientific viewport. |
| `kiss3d` | 0.40 (Jan 2026, wgpu rewrite) | Useful for documentation examples; rendering model too limited for main viewport. |
| `three-d` | 0.18 | Reasonable middle ground; doesn't add much over wgpu+egui for our use case. |

Required capabilities — extruded copper layers, current-density heatmaps, E/H cut planes, 3D far-field patterns — all reduce to a few hundred lines of wgpu plus `glam`, with full shader control. Rerun complements: at any time the solver can `rerun.log()` arbitrary tensors, meshes, scalars on a timeline, giving us a free interactive multi-view debugger.

---

## Plotting

| Channel | Choice | Why |
|---|---|---|
| In-app live plots | **`egui_plot`** | Native immediate-mode plot widget; pans, zooms, axis-linked panels — perfect for sweeping S-params across ports. |
| Static / publication export | **`plotters`** | Pure Rust; SVG/PNG/Cairo backends; publication quality. |
| Power-user / scripting | Matplotlib / plotly in user notebooks via the PyO3 bindings | The cleanest path — don't host Python inside the GUI. |

---

## Logging, error handling, async

- **Errors:** `thiserror` for library crates, `anyhow` for the CLI/GUI binaries.
- **Logging / tracing:** `tracing` + `tracing-subscriber`. Rerun integration is additive on top.
- **Async (where used):** `tokio` for any background services (file watch, surrogate training orchestration). The solver itself is synchronous, GPU-driven.
- **Config / project files:** TOML for human-edited (`yee.toml`); HDF5 or Arrow IPC for solver output (large field arrays, time-series).

---

## Quick reference — workspace `Cargo.toml` excerpt

```toml
[workspace.package]
rust-version = "1.85"
edition      = "2024"
license      = "GPL-3.0-or-later"

[workspace.dependencies]
# CUDA (one binding for the whole project)
cudarc = { version = "0.19", default-features = false, features = [
    "std", "driver", "nvrtc",
    "cublas", "cublaslt", "cusolver", "cusparse", "cufft", "curand",
    "f16",
    "cuda-version-from-build-system",
    "dynamic-loading",
] }

# Linear algebra
faer       = "0.23"
nalgebra   = "0.34"
ndarray    = "0.17"
sprs       = "0.11"
num-complex = "0.4"

# Geometry
glam        = "0.30"
parry3d     = "0.22"
geo         = "0.30"
i_overlay   = "4"

# CAD / Meshing
opencascade     = "0.2"   # bschwind/opencascade-rs
opencascade-sys = "0.2"
# Gmsh: in-tree bindgen against gmshc.h 4.15+, not the stale rgmsh crate

# GUI / Viz
egui        = "0.34"
eframe      = "0.34"
egui_plot   = "0.34"
egui_dock   = "0.18"
wgpu        = "26"
rerun       = { version = "0.31", default-features = false, features = ["sdk"] }
plotters    = "0.3"

# Python bindings
pyo3   = { version = "0.28", features = ["abi3-py310", "extension-module"] }
numpy  = "0.28"

# Diagnostics
tracing            = "0.1"
tracing-subscriber = "0.3"
thiserror          = "1"
anyhow             = "1"
```

---

## License sanity check

Every dependency above is verified compatible with **GPL v3** distribution:

- **Gmsh** — GPL v2+ with linking exception (FAQ-confirmed GPL-3 compatible).
- **OCCT** — LGPL 2.1 with exception (FAQ-confirmed GPL-2+ compatible).
- **`opencascade-rs`** — LGPL-2.1.
- **`egui`, `eframe`, `egui_plot`, `parry`, `wgpu`, `glam`, `nalgebra`, `ndarray`, `faer`, `rerun`, `plotters`, `pyo3`, `numpy`, `cudarc`** — all permissive (MIT/Apache/BSD), no conflict.
- **NVIDIA CUDA libraries (cuBLAS, cuSOLVER, cuSPARSE, cuFFT, NCCL)** — proprietary but **dynamically linked at runtime**, which is the same posture every CUDA-using OSS project takes; we document the dependency clearly.

No license conflicts in the recommended stack.

---

## Known watch-outs

1. **`cudarc` is "pre-alpha"** per its own README; history includes breaking minor releases (notably 0.13 → 0.14). Pin to exact minor version and route all CUDA access through a thin internal abstraction.
2. **Gmsh Rust bindings (`rgmsh`/`gmsh-sys`) are abandoned.** We bake our own.
3. **`opencascade-rs` cold build is 5–15 minutes** (rebuilds OCCT). Cache CI builds aggressively (`sccache` or workspace artifact cache).
4. **egui has breaking minor releases roughly quarterly.** Isolate UI behind a stable internal trait.
5. **PyO3 ships one breaking change per minor release.** Pin and migrate deliberately.
6. **Bevy 0.x ecosystem resets every ~3 months** — another reason we are not using Bevy as the UI shell.
7. **OCCT 8.0 targets May 2026** with significant changes (C++17, BRepGraph topology); expect a wave of `opencascade-rs` work after that lands.
8. **Dense LU at n ≥ 50k overflows a single GPU.** Iterative methods in Phase 1; MLFMA/ACA in Phase 4.
9. **`rustc_codegen_nvvm` (Rust on GPU) is not production-ready.** Re-evaluate late 2026. Until then, kernels are CUDA-C JITed via NVRTC.
