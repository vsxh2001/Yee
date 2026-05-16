"""Tests for `yee.touchstone.read` / `yee.touchstone.write`.

The Phase 0 solver does not yet produce real `SParameters`, so the roundtrip
exercises the explicit `touchstone.write(dict)` → `touchstone.read(path)`
path with a hand-constructed file dict. This validates the binding's
dict-encoding logic end-to-end against the Rust Touchstone codec.
"""

import numpy as np
import pytest

import yee


def test_sparams_write_and_read_roundtrip(tmp_path):
    """Hand-write a Touchstone file, read it back, verify the contents."""
    n_ports = 1
    freq_hz = np.array([1.0e9, 1.25e9, 1.5e9], dtype=np.float64)
    data = np.array(
        [
            [[0.1 + 0.2j]],
            [[0.15 + 0.25j]],
            [[0.2 + 0.3j]],
        ],
        dtype=np.complex128,
    )
    file_dict = {
        "z0": 50.0,
        "freq_unit": "Hz",
        "format": "RI",
        "n_ports": n_ports,
        "freq_hz": freq_hz,
        "data": data,
        "comments": ["written by yee-py test suite"],
    }

    path = tmp_path / "test.s1p"
    yee.touchstone.write(str(path), file_dict)

    parsed = yee.touchstone.read(str(path))
    assert parsed["n_ports"] == 1
    assert parsed["z0"] == 50.0
    assert parsed["freq_unit"] == "Hz"
    assert parsed["format"] == "RI"
    assert parsed["freq_hz"].shape == (3,)
    assert parsed["data"].shape == (3, 1, 1)
    assert parsed["data"].dtype == np.complex128
    np.testing.assert_allclose(parsed["freq_hz"], freq_hz)
    np.testing.assert_allclose(parsed["data"], data)
    assert parsed["comments"] == ["written by yee-py test suite"]


def test_touchstone_read_rejects_missing_file(tmp_path):
    """Reading a non-existent path must raise IOError (mapped from yee_io::Error::Io)."""
    bogus = tmp_path / "does_not_exist.s1p"
    with pytest.raises(IOError):
        yee.touchstone.read(str(bogus))
