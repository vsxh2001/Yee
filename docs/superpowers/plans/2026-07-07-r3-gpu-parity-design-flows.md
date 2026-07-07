# Plan — R.3 GPU parity for the design flows

**Spec:** `docs/superpowers/specs/2026-07-07-r3-gpu-parity-design-flows-design.md`

1. **Per-face CPML**: `Params.axes_mask` → `faces_mask` (6 bits, `2·axis+side`);
   host builds it from `CpmlConfig.faces`; WGSL `pml_depth` checks min/max bits;
   drop the `faces_are_axis_symmetric` rejection. Gate `compute-016`
   (`gpu_perface_cpml_parity.rs`): open-top box vs CPU, family-relative field
   comparison (compute-005 idiom) + z-max absorption evidence.
2. **Aperture ports**: append-only `drv_idx`/`drv_data` extension (cells table +
   per-port constants + v_src series + v_prev state); WGSL `apply_aperture_ports`
   (one invocation per port, serial cell loops, explicit resum); dispatch between
   `apply_ports` and `record_probes`; drop the aperture rejection. Gate
   `compute-015` (`gpu_aperture_parity.rs`): mini S.10 board, probe series vs CPU
   family-relative (rel L2 < 1e-4, L∞ < 1e-3).
3. Verify on llvmpipe locally (`cargo test -p yee-compute --release -- --include-ignored`),
   plus fmt/clippy floor and the fast engine tests (drive/spec unchanged upstream).
4. ADR-0196; RF-TOOL-ROADMAP R.3 row + footer; SUMMARY entry. Engine NTFF stays
   CPU-only (documented out of scope). Commit + push; continue to R.4.
