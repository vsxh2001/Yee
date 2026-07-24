//! FS.1c: the wgpu GPU backend has no thin-wire kernel (walking skeleton —
//! the CPU-only `ThinWire` subcell is the E.0-scope deliverable, a GPU port
//! is a follow-on). `GpuFdtd::with_drive` must reject a drive that carries
//! any `ThinWire` with a named [`ComputeError::Unsupported`] **before**
//! touching the adapter/device at all (same posture as the existing
//! aperture-port-recording and resistive-sheet rejections in `gpu.rs`), so
//! this doesn't need real GPU hardware to exercise the rejection path
//! itself — it fails validation first.

#![cfg(feature = "gpu")]

use yee_compute::{Boundary, ComputeError, Drive, FdtdSpec, Fields, GpuFdtd, Materials, ThinWire};

const NX: usize = 12;
const NY: usize = 12;
const NZ: usize = 12;
const DX: f64 = 1e-3;

#[test]
fn gpu_rejects_a_drive_carrying_a_thin_wire() {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let fields = Fields::zero(&spec);

    let mut drive = Drive::default();
    drive.thin_wires.push(ThinWire {
        i: 6,
        j: 6,
        k_lo: 2,
        k_hi: 9,
        radius_m: 0.1e-3,
        feed_k: Some(5),
    });

    let result = GpuFdtd::with_drive(spec, fields, Materials::default(), Boundary::None, drive, 1);
    match result {
        Err(ComputeError::Unsupported(msg)) => {
            assert!(
                msg.contains("thin-wire"),
                "expected the thin-wire-specific Unsupported message, got: {msg}"
            );
        }
        other => panic!("expected ComputeError::Unsupported for a thin-wire drive, got {other:?}"),
    }
}
