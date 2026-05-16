# yee-validation

Validation-case aggregator for the Yee electromagnetic simulator.

This crate collects the canonical validation cases from across the
workspace (`mom-001` dipole, `cpml-001` PML reflection, `ntff-001`
near-to-far-field, `dispersive-001` Drude slab, ...) into a single
runner that emits a structured `Report` consumable as JSON or
Markdown.

## Scope

**Phase 1.validation.0 walking skeleton.** The aggregator currently
ships the report **schema** and entry point. Every case body reports
`Skipped` with an explanatory message: the underlying validation
suites live in private `#[cfg(test)]` modules of `yee-mom` /
`yee-fdtd` and are not yet callable from a sibling crate. See the
crate-level docstring for the Phase 1.validation.1 unblock paths.

The full mom-001 24x176 gate (~7-8 min wall time, `cargo test
--release -p yee-mom`) and the FDTD integration tests continue to
run unchanged in their respective `cargo test` paths.

## Usage

```rust
use yee_validation::{run_all, Report};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let report: Report = run_all();

    // JSON for machine consumers (CI dashboards, regression bots).
    println!("{}", report.to_json()?);

    // Markdown for the GitHub job-summary step.
    std::fs::write("validation-report.md", report.to_markdown())?;

    Ok(())
}
```

## Report shape

```json
{
  "generated_at": "epoch+1747426812s",
  "git_sha": null,
  "cases": [
    {
      "id": "mom-001-fast",
      "description": "Half-wave dipole, NEC-4 reference ...",
      "status": "Skipped",
      "message": "yee_mom::fixtures::thin_cylinder is private; ...",
      "wall_time_seconds": 0.0
    }
  ]
}
```
