fn main() -> pyo3_stub_gen::Result<()> {
    cosmolkit_py::stub_info()?.generate()?;
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("cosmolkit.pyi");
    let text = std::fs::read_to_string(&path)?;
    let text = text.replace("__all__ = [\n",
        "__all__ = [\n    \"DrawingError\",\n    \"OperationError\",\n    \"__version__\",\n    \"_binding_profile\",\n");
    let text = format!(
        "from __future__ import annotations\n{text}\n{}",
        r#"class DrawingError(builtins.ValueError):
    domain: builtins.str
    kind: builtins.str
    # Context attributes exist only for applicable Rust variants.
    width: builtins.int
    height: builtins.int
    field: builtins.str
    actual: builtins.int
    expected: builtins.int
    row: typing.Optional[builtins.int]
    reason: builtins.str

class OperationError(builtins.ValueError):
    domain: builtins.str
    kind: builtins.str

__version__: builtins.str
_binding_profile: builtins.str
"#
    );
    let mut text = text;
    for name in [
        "SmilesError",
        "SmilesWriteError",
        "MorganReadError",
        "FingerprintError",
    ] {
        text.push_str(&format!("\nclass {name}(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n"));
        text = text.replace("__all__ = [\n", &format!("__all__ = [\n    \"{name}\",\n"));
    }
    text = text.replace("class FingerprintError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n", "class FingerprintError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    # Context fields exist only on applicable variants.\n    index: builtins.int\n    size: builtins.int\n    left: builtins.int\n    right: builtins.int\n    factor: builtins.int\n    n_bits: builtins.int\n    value: builtins.float\n    site: builtins.str\n    what: builtins.str\n    reason: builtins.str\n");
    std::fs::write(path, text)?;
    Ok(())
}
