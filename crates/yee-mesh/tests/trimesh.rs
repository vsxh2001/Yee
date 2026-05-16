//! Integration tests for `TriMesh::new`.
//!
//! These exercise the validating constructor: success when triangles and
//! tags have matching lengths, and rejection otherwise.

use nalgebra::Vector3;
use yee_mesh::{Error, TriMesh};

#[test]
fn trimesh_new_accepts_matching_triangles_and_tags() {
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(1.0, 1.0, 0.0),
    ];
    let triangles = vec![[0, 1, 2], [1, 3, 2]];
    let tags = vec![1, 1];

    let mesh = TriMesh::new(vertices.clone(), triangles.clone(), tags.clone())
        .expect("matching triangles and tags must build");
    assert_eq!(mesh.n_tris(), 2);
    assert_eq!(mesh.vertices.len(), vertices.len());
    assert_eq!(mesh.triangles, triangles);
    assert_eq!(mesh.tags, tags);
}

#[test]
fn trimesh_new_rejects_mismatched_tag_count() {
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
    ];
    let triangles = vec![[0, 1, 2]];
    let tags = vec![1, 2]; // length mismatch

    let err = TriMesh::new(vertices, triangles, tags)
        .expect_err("mismatched tag length must be rejected");
    match err {
        Error::Invalid(msg) => {
            assert!(
                msg.contains("triangles") && msg.contains("tags"),
                "error message should mention triangles and tags: {msg}"
            );
        }
        other => panic!("expected Error::Invalid, got {other:?}"),
    }
}
