# Summary

- [Introduction](introduction.md)

# Theory

- [Planar Method of Moments](theory/planar-mom.md)
- [RWG Basis, MPIE, and Phase 1.0 Quadrature](theory/mom-rwg-mpie.md)
- [Finite-Difference Time-Domain](theory/fdtd.md)
- [FDTD: CPML, NTFF, TF/SF, lumped, dispersive](theory/fdtd-details.md)
- [CUDA backend](theory/cuda-backend.md)
- [Touchstone v1.1 File Format](theory/touchstone-format.md)
- [GP Surrogates and Bayesian Optimization](theory/surrogate-gp-bo.md)
- [Multi-objective Optimization and Active Learning](theory/multi-objective-and-active-learning.md)

# Tutorials

- [Hello, microstrip](tutorials/01-microstrip-line.md)
- [Half-wave dipole from Python](tutorials/02-dipole-from-python.md)
- [FDTD cavity resonance](tutorials/03-fdtd-cavity.md)

# Decisions

- [ADR-0001: GPL v3 license](decisions/0001-license-gplv3.md)
- [ADR-0002: Rust MSRV 1.88](decisions/0002-rust-msrv-1.88.md)
- [ADR-0003: PyO3 abi3-py310](decisions/0003-pyo3-abi3-py310.md)
- [ADR-0004: egui pinned 0.32](decisions/0004-egui-pinned-to-0.32.md)
- [ADR-0005: NEC-4 vs Balanis mom-001](decisions/0005-nec4-vs-balanis-mom-001.md)
- [ADR-0006: cudarc pre-alpha pin](decisions/0006-cudarc-prealpha-pin.md)
- [ADR-0007: yee-bench criterion benches](decisions/0007-yee-bench-criterion.md)
- [ADR-0008: validation aggregator JSON + PNG](decisions/0008-validation-aggregator-json-and-png.md)
- [ADR-0009: hand-rolled Gaussian process surrogate](decisions/0009-gaussian-process-surrogate.md)
- [ADR-0010: BO Expected Improvement before NSGA-II](decisions/0010-bayesian-optimization-ei-first.md)
- [ADR-0011: Toolchain bump to Rust 1.92 / egui 0.34 / wgpu 29](decisions/0011-toolchain-bump-rust-1-92-egui-0-34.md)
- [ADR-0012: NSGA-II as a separate module from BO](decisions/0012-multi-objective-nsga2.md)
- [ADR-0013: Active learning via variance acquisition](decisions/0013-active-learning-variance-acquisition.md)
- [ADR-0014: TF/SF slab geometry, finite box deferred](decisions/0014-tfsf-slab-not-finite-box.md)
- [ADR-0015: PlanarMoM GreensSpec builder](decisions/0015-planar-mom-greens-spec-builder.md)
- [ADR-0016: yee-py validation binding (slow path gated)](decisions/0016-yee-py-validation-binding-slow-path-gated.md)
- [ADR-0017: FDTD lumped RLC port](decisions/0017-fdtd-lumped-rlc-port.md)
- [ADR-0018: yee-gui validation panel](decisions/0018-yee-gui-validation-panel.md)
