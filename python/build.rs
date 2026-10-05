fn main() {
    if std::env::var_os("CARGO_FEATURE_EXTENSION_MODULE").is_some() {
        pyo3_build_config::add_extension_module_link_args();
    } else if std::env::var_os("CARGO_FEATURE_STUBGEN").is_some()
        || std::env::var_os("CARGO_FEATURE_PYTHON_EMBED_TESTS").is_some()
    {
        // Development executables embed the selected interpreter. Use its actual
        // library directory, not a hard-coded Python version or shell workaround.
        // Release extension modules must not carry this development-only rpath.
        pyo3_build_config::add_libpython_rpath_link_args();
    }
}
