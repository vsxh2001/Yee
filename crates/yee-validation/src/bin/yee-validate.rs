//! `yee-validate` — command-line driver for the validation aggregator.
//!
//! Runs every registered case in [`yee_validation::Report::run_all`],
//! emits the resulting [`yee_validation::Report`] as pretty-printed
//! JSON on stdout, and exits non-zero iff any case returned
//! [`yee_validation::CaseStatus::Failed`]. Skipped cases do not count
//! as failures.

use yee_validation::Report;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let report = Report::run_all();
    println!("{}", report.to_json()?);
    if report.has_failures() {
        std::process::exit(1);
    }
    Ok(())
}
