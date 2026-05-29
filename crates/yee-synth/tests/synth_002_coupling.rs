//! `synth-002` — all-pole coupling coefficients + external Q gate.
//!
//! For a worked example (Chebyshev 0.5 dB, N=3, FBW=0.10) the synthesized
//! coupling design must:
//! - be **symmetric** for a synchronous all-pole filter: `k_12 == k_23`,
//!   `Qe_in == Qe_out`;
//! - match the spec §2.5 closed form recomputed independently in-test to
//!   `≤ 1e-9` (self-consistency); and
//! - reflect the matrix relation `k_{i,i+1} = FBW · M[i][i+1]`.

use yee_synth::{Approximation, coupling_design, prototype};

const FBW: f64 = 0.10;

#[test]
fn chebyshev_0p5db_n3_coupling() {
    let proto = prototype(Approximation::Chebyshev { ripple_db: 0.5 }, 3);
    let design = coupling_design(&proto, FBW);
    let g = &proto.g; // [g0, g1, g2, g3, g4]
    let n = proto.order();
    assert_eq!(n, 3);

    // ---- symmetry (synchronous, symmetric prototype) ----------------------
    assert_eq!(design.k.len(), n - 1, "k must have N-1 = 2 entries");
    assert!(
        (design.k[0] - design.k[1]).abs() < 1e-12,
        "k_12 ({}) != k_23 ({}) for a symmetric synchronous filter",
        design.k[0],
        design.k[1]
    );
    assert!(
        (design.qe_in - design.qe_out).abs() < 1e-12,
        "Qe_in ({}) != Qe_out ({}) for a symmetric synchronous filter",
        design.qe_in,
        design.qe_out
    );

    // ---- closed-form self-consistency (spec §2.5) -------------------------
    // k_{i,i+1} = FBW / sqrt(g_i g_{i+1})
    let k12 = FBW / (g[1] * g[2]).sqrt();
    let k23 = FBW / (g[2] * g[3]).sqrt();
    assert!((design.k[0] - k12).abs() < 1e-9, "k_12 mismatch vs §2.5");
    assert!((design.k[1] - k23).abs() < 1e-9, "k_23 mismatch vs §2.5");

    let qe_in = g[0] * g[1] / FBW;
    let qe_out = g[3] * g[4] / FBW;
    assert!(
        (design.qe_in - qe_in).abs() < 1e-9,
        "Qe_in mismatch vs §2.5"
    );
    assert!(
        (design.qe_out - qe_out).abs() < 1e-9,
        "Qe_out mismatch vs §2.5"
    );

    // ---- normalized coupling matrix relation: k = FBW · M -----------------
    // M[i][i+1] = 1/sqrt(g_i g_{i+1}); zero diagonal (synchronous).
    assert_eq!(design.m.len(), n);
    for row in &design.m {
        assert_eq!(row.len(), n);
    }
    for i in 0..n {
        assert_eq!(design.m[i][i], 0.0, "synchronous filter has zero diagonal");
    }
    assert!(
        (design.m[0][1] - design.m[1][0]).abs() < 1e-15,
        "coupling matrix must be symmetric"
    );
    assert!(
        (design.m[1][2] - design.m[2][1]).abs() < 1e-15,
        "coupling matrix must be symmetric"
    );
    assert!(
        (design.k[0] - FBW * design.m[0][1]).abs() < 1e-12,
        "k_12 must equal FBW · M[0][1]"
    );
    assert!(
        (design.k[1] - FBW * design.m[1][2]).abs() < 1e-12,
        "k_23 must equal FBW · M[1][2]"
    );
    // Non-adjacent entries are zero for an all-pole network.
    assert_eq!(design.m[0][2], 0.0);
    assert_eq!(design.m[2][0], 0.0);
}
