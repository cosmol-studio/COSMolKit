fn main() {
    if std::env::var_os("CARGO_FEATURE_EXTENSION").is_some() {
        pyo3_build_config::add_extension_module_link_args();
    } else {
        pyo3_build_config::add_libpython_rpath_link_args();
    }
}
