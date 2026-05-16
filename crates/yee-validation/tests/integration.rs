//! Integration test: exercise [`yee_validation::Report::run_all`]
//! end-to-end and assert the `mom-001` case passes against the NEC-4
//! finite-radius reference.
//!
//! Gated behind `#[ignore]` because the underlying 24x176 dipole
//! solve runs ~7-8 min in `--release`. Run explicitly with
//! `cargo test -p yee-validation --release -- --include-ignored`.

#[test]
#[ignore = "slow: ~8 min for mom-001"]
fn mom_001_passes_through_aggregator() {
    let report = yee_validation::Report::run_all();
    let mom_001 = report
        .cases
        .iter()
        .find(|c| c.id == "mom-001")
        .expect("mom-001 case present");
    assert!(
        matches!(mom_001.status, yee_validation::Status::Passed),
        "mom-001 failed in aggregator: {}",
        mom_001.notes
    );
}
