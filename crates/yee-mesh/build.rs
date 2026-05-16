//! Build script for `yee-mesh`.
//!
//! When the `gmsh` feature is enabled AND `$GMSH_SDK_ROOT` is set, this
//! invokes `bindgen` against `$GMSH_SDK_ROOT/include/gmshc.h` and writes
//! `bindings.rs` into `$OUT_DIR`.
//!
//! When the feature is enabled but `$GMSH_SDK_ROOT` is unset, we emit a
//! `cargo:warning=` and still write an empty `bindings.rs` stub so that
//! type-only consumers can compile.
//!
//! Without the `gmsh` feature, this script is a no-op aside from
//! rerun-if-env hints.

fn main() {
    // Rerun whenever any env var controlling SDK discovery or version changes.
    println!("cargo:rerun-if-env-changed=GMSH_SDK_ROOT");
    println!("cargo:rerun-if-env-changed=GMSH_VERSION");
    println!("cargo:rerun-if-changed=build.rs");

    #[cfg(feature = "gmsh")]
    generate_bindings();
}

#[cfg(feature = "gmsh")]
fn generate_bindings() {
    use std::env;
    use std::fs;
    use std::path::PathBuf;

    let out_dir = PathBuf::from(
        env::var_os("OUT_DIR")
            .unwrap_or_else(|| panic!("cargo did not provide OUT_DIR to build.rs")),
    );
    let bindings_path = out_dir.join("bindings.rs");

    let Some(sdk_root) = env::var_os("GMSH_SDK_ROOT") else {
        println!(
            "cargo:warning=GMSH_SDK_ROOT is unset; yee-mesh's `gmsh` feature \
             is enabled but bindgen will be skipped. Writing empty bindings \
             stub so type-only consumers compile."
        );
        fs::write(&bindings_path, "// GMSH_SDK_ROOT unset; empty stub.\n").unwrap_or_else(|e| {
            panic!(
                "failed to write empty bindings stub to {}: {e}",
                bindings_path.display()
            )
        });
        return;
    };

    let sdk_root = PathBuf::from(sdk_root);
    let header = sdk_root.join("include").join("gmshc.h");
    println!("cargo:rerun-if-changed={}", header.display());

    let bindings = bindgen::Builder::default()
        .header(header.to_string_lossy())
        // gmshc.h is a C header; pin the standard so clang doesn't pick up
        // a host-dependent default.
        .clang_arg("--std=c99")
        .allowlist_function("gmsh.*")
        .allowlist_type("gmsh.*")
        .allowlist_var("gmsh.*")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .unwrap_or_else(|e| {
            panic!(
                "bindgen failed to generate bindings from {}: {e}",
                header.display()
            )
        });

    bindings.write_to_file(&bindings_path).unwrap_or_else(|e| {
        panic!(
            "failed to write generated bindings to {}: {e}",
            bindings_path.display()
        )
    });
}
