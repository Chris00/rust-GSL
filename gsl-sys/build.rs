use std::{env, path::PathBuf};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut inc_dirs = vec![];
    let mut link_paths = vec![];
    let mut libs = vec![];
    let gsl = pkg_config::Config::new()
        .atleast_version("2.1")
        .probe("gsl");
    if let Ok(mut gsl) = gsl {
        inc_dirs.append(&mut gsl.include_paths);
        link_paths.append(&mut gsl.link_paths);
        libs.append(&mut gsl.libs);
    }
    #[cfg(target_family = "windows")]
    {
        if libs.is_empty() {
            let mut vcpkg = vcpkg::Config::new()
                .emit_includes(false)
                .find_package("gsl")?;
            inc_dirs.append(&mut vcpkg.include_paths);
            link_paths.append(&mut vcpkg.link_paths);
            for dll in vcpkg.found_dlls {
                libs.push(dll.display().to_string());
            }
        }
    }
    if libs.is_empty() {
        libs.push("gsl".into());
        libs.push("gslcblas".into());
    }
    for l in libs {
        println!("cargo:rustc-link-lib={}", l);
    }
    for l in link_paths {
        println!("cargo:rustc-link-arg=-L{}", l.display());
    }

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    println!("cargo:rerun-if-changed=src/wrapper.h");
    let mut builder = bindgen::Builder::default()
        .header("src/wrapper.h")
        .blocklist_type("FILE")
        .blocklist_type("_IO_FILE")
        .allowlist_file(".*/gsl_.*");
    for dir in inc_dirs {
        builder = builder.clang_arg(format!("-I{}", dir.display()))
    }
    let bindings = builder.generate()?;
    bindings.write_to_file(out_path.join("bindings.rs"))?;

    Ok(())
}
