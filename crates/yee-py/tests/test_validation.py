"""yee.run_validation bindings smoke test.

Don't actually call run_validation() by default -- it would take ~8 min
for the real mom-001 aggregator run. We only smoke the imports and
verify the classes exist with the documented methods.

The real-aggregator path is gated behind YEE_RUN_VALIDATION=1.
"""

import os

import pytest

import yee


def test_validation_imports_present():
    assert hasattr(yee, "run_validation")
    assert hasattr(yee, "ValidationReport")
    assert hasattr(yee, "ValidationCase")


def test_validation_report_has_expected_methods():
    # We don't construct one -- just confirm the methods are advertised
    # on the class.
    assert hasattr(yee.ValidationReport, "cases")
    assert hasattr(yee.ValidationReport, "has_failures")
    assert hasattr(yee.ValidationReport, "to_json")


def test_validation_case_has_expected_getters():
    assert hasattr(yee.ValidationCase, "id")
    assert hasattr(yee.ValidationCase, "status")
    assert hasattr(yee.ValidationCase, "notes")
    assert hasattr(yee.ValidationCase, "wall_time_seconds")
    assert hasattr(yee.ValidationCase, "plot_paths")


@pytest.mark.skipif(
    os.environ.get("YEE_RUN_VALIDATION") != "1",
    reason="set YEE_RUN_VALIDATION=1 to invoke the real (~8 min) aggregator",
)
def test_run_validation_real():
    report = yee.run_validation()
    assert isinstance(report.to_json(), str)
    mom_001 = next((c for c in report.cases if c.id == "mom-001"), None)
    assert mom_001 is not None
