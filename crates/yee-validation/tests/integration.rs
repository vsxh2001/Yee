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

/// Assert mom-001 emits non-trivial plot artifacts under
/// `validation/results/`. Ignored alongside the impedance gate
/// because it invokes the full aggregator (fine-mesh ~8 min +
/// coarse-mesh plot sweep ~3.5 min on top).
#[test]
#[ignore = "slow: invokes the real aggregator (mom-001 ~8 min)"]
fn mom_001_emits_plot_artifacts() {
    let report = yee_validation::Report::run_all();
    let mom_001 = report.cases.iter().find(|c| c.id == "mom-001").unwrap();
    assert!(
        !mom_001.plot_paths.is_empty(),
        "mom-001 should emit plot artifacts"
    );
    for p in &mom_001.plot_paths {
        assert!(p.exists(), "plot path missing: {}", p.display());
        let size = std::fs::metadata(p).unwrap().len();
        assert!(
            size > 1024,
            "plot {} too small ({} bytes)",
            p.display(),
            size
        );
    }
}

/// fem-eig-001: assert the WR-90 cavity eigenmode case is registered
/// in `run_all` and reports `Passed`.
///
/// The driver runs in ~7 s release on a 12×9×15 mesh; gated `#[ignore]`
/// here because calling `run_all` also pulls in the slow mom-001 solve
/// (~8 min). Run with
/// `cargo test -p yee-validation --release -- --include-ignored`.
#[test]
#[ignore = "slow: aggregator invokes mom-001 (~8 min) + fem-eig-001 (~7 s)"]
fn fem_eig_001_registered_and_passes_through_aggregator() {
    let report = yee_validation::Report::run_all();
    let case = report
        .cases
        .iter()
        .find(|c| c.id == "fem-eig-001")
        .expect("fem-eig-001 must be registered in run_all");
    assert!(
        matches!(case.status, yee_validation::CaseStatus::Passed),
        "fem-eig-001 did not pass: {}",
        case.notes
    );
}

/// mom-002: assert the microstrip case lands in
/// [`yee_validation::Status::Passed`] when the aggregator runs all
/// cases. Ignored because the aggregator pulls in mom-001 (~8 min)
/// even though mom-002 itself is quick (~tens of seconds at the
/// 30x2 strip mesh).
#[test]
#[ignore = "slow: aggregator invokes mom-001 (~8 min) + mom-002 (~30s)"]
fn mom_002_passes_with_loose_tolerance() {
    let report = yee_validation::Report::run_all();
    let mom_002 = report
        .cases
        .iter()
        .find(|c| c.id == "mom-002")
        .expect("mom-002 case present");
    assert!(
        matches!(mom_002.status, yee_validation::Status::Passed),
        "mom-002 failed in aggregator: {}",
        mom_002.notes
    );
    assert!(
        !mom_002.plot_paths.is_empty(),
        "mom-002 should emit plot artifacts: {}",
        mom_002.notes
    );
    for p in &mom_002.plot_paths {
        assert!(p.exists(), "plot path missing: {}", p.display());
        let size = std::fs::metadata(p).unwrap().len();
        assert!(
            size > 1024,
            "plot {} too small ({} bytes)",
            p.display(),
            size
        );
    }
}
