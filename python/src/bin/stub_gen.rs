#[cfg(not(feature = "drawing-bindings"))]
use cosmolkit_py as cosmolkit;

#[cfg(not(feature = "drawing-bindings"))]
include!("unmigrated_stub_gen.rs");

#[cfg(feature = "drawing-bindings")]
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
    std::fs::write(path, text)?;
    Ok(())
}
