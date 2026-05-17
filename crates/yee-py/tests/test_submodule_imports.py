"""Regression test: child submodules must be importable via `from yee.X import Y`."""


def test_import_yee_touchstone_module():
    import yee.touchstone  # noqa: F401
    assert yee.touchstone is not None


def test_from_yee_touchstone_import_read():
    from yee.touchstone import read  # noqa: F401
    assert callable(read)


def test_from_yee_touchstone_import_write():
    from yee.touchstone import write  # noqa: F401
    assert callable(write)


def test_yee_touchstone_attribute_access_still_works():
    # The pre-fix workaround path used in examples/python/touchstone_workflow.ipynb.
    import yee
    assert callable(yee.touchstone.read)
    assert callable(yee.touchstone.write)
