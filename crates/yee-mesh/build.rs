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
    // Always rerun if the env var that controls SDK discovery changes.
    println!("cargo:rerun-if-env-changed=GMSH_SDK_ROOT");
    println!("cargo:rerun-if-changed=build.rs");

    #[cfg(feature = "gmsh")]
    generate_bindings();
}

#[cfg(feature = "gmsh")]
fn generate_bindings() {
    use std::env;
    use std::fs;
    use std::path::PathBuf;

    let out_dir = PathBuf::from(env::var_os("OUT_DIR").expect("cargo sets OUT_DIR"));
    let bindings_path = out_dir.join("bindings.rs");

    let Some(sdk_root) = env::var_os("GMSH_SDK_ROOT") else {
        println!(
            "cargo:warning=GMSH_SDK_ROOT is unset; yee-mesh's `gmsh` feature \
             is enabled but bindgen will be skipped. Writing empty bindings \
             stub so type-only consumers compile."
        );
        fs::write(&bindings_path, "// GMSH_SDK_ROOT unset; empty stub.\n")
            .expect("write empty bindings stub");
        return;
    };

    let sdk_root = PathBuf::from(sdk_root);
    let header = sdk_root.join("include").join("gmshc.h");
    println!("cargo:rerun-if-changed={}", header.display());

    let bindings = bindgen::Builder::default()
        .header(header.to_string_lossy())
        .allowlist_function("gmsh.*")
        .allowlist_type("gmsh.*")
        .allowlist_var("gmsh.*")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("bindgen failed to generate gmshc bindings");

    bindings
        .write_to_file(&bindings_path)
        .expect("write bindings.rs");
}
