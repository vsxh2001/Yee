"""Tests for `yee.TriMesh`."""

import numpy as np
import pytest

import yee


def test_trimesh_constructs_from_numpy():
    v = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        dtype=np.float64,
    )
    t = np.array([[0, 1, 2]], dtype=np.uint32)
    g = np.array([1], dtype=np.uint32)
    mesh = yee.TriMesh(v, t, g)
    assert mesh.n_tris() == 1


def test_trimesh_rejects_bad_vertex_shape():
    # 2D coordinates instead of 3D — must fail with a shape ValueError.
    v = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float64)
    t = np.array([[0, 1, 2]], dtype=np.uint32)
    g = np.array([1], dtype=np.uint32)
    with pytest.raises(ValueError, match="shape"):
        yee.TriMesh(v, t, g)


def test_trimesh_rejects_bad_triangle_shape():
    # Only two indices per triangle.
    v = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        dtype=np.float64,
    )
    t = np.array([[0, 1]], dtype=np.uint32)
    g = np.array([1], dtype=np.uint32)
    with pytest.raises(ValueError, match="shape"):
        yee.TriMesh(v, t, g)


def test_trimesh_rejects_length_mismatch():
    v = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        dtype=np.float64,
    )
    t = np.array([[0, 1, 2]], dtype=np.uint32)
    g = np.array([0, 0], dtype=np.uint32)
    with pytest.raises(ValueError, match="length"):
        yee.TriMesh(v, t, g)


def test_trimesh_getters_return_arrays_with_correct_shape():
    v = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=np.float64,
    )
    t = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.uint32)
    g = np.array([1, 2], dtype=np.uint32)
    mesh = yee.TriMesh(v, t, g)

    verts = mesh.vertices
    tris = mesh.triangles
    tags = mesh.tags

    assert verts.shape == (4, 3)
    assert verts.dtype == np.float64
    np.testing.assert_array_equal(verts, v)

    assert tris.shape == (2, 3)
    assert tris.dtype == np.uint32
    np.testing.assert_array_equal(tris, t)

    assert tags.shape == (2,)
    assert tags.dtype == np.uint32
    np.testing.assert_array_equal(tags, g)
