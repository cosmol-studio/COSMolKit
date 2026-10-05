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

class DrawingWriteError(builtins.OSError):
    domain: builtins.str
    kind: builtins.str

class OperationError(builtins.ValueError):
    domain: builtins.str
    kind: builtins.str

__version__: builtins.str
_binding_profile: builtins.str
"#
    );
    let text = text.replace("__all__ = [\n", "__all__ = [\n    \"DrawingWriteError\",\n");
    let mut text = text;
    for name in [
        "SmilesError",
        "SmilesWriteError",
        "MorganReadError",
        "FingerprintError",
        "DescriptorReadError",
        "DescriptorError",
    ] {
        text.push_str(&format!("\nclass {name}(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n"));
        text = text.replace("__all__ = [\n", &format!("__all__ = [\n    \"{name}\",\n"));
    }
    text = text.replace("class FingerprintError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n", "class FingerprintError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    # Context fields exist only on applicable variants.\n    index: builtins.int\n    size: builtins.int\n    left: builtins.int\n    right: builtins.int\n    factor: builtins.int\n    n_bits: builtins.int\n    value: builtins.float\n    site: builtins.str\n    what: builtins.str\n    reason: builtins.str\n");
    // Constants share the exact public facade iterator used by native publication.
    let constants = ::cosmolkit::Element::iter_with_dummy()
        .map(|element| {
            let name = if element == ::cosmolkit::Element::DUMMY {
                "DUMMY".to_owned()
            } else {
                element.symbol().to_ascii_uppercase()
            };
            format!("    {name}: typing.ClassVar[Element]\n")
        })
        .collect::<String>();
    let start = text
        .find("class Element:\n")
        .expect("generated Element class");
    let end = text[start..]
        .find("\n@typing.final\nclass ElementInfo:")
        .map(|offset| start + offset)
        .expect("generated ElementInfo class");
    let class = &text[start..end];
    // PyO3's eq slots use `value` and return NotImplemented for foreign values.
    // stub-gen 0.23 only emits __eq__(other)->bool for the eq class option.
    let eq = "    def __eq__(self, other: builtins.object, /) -> builtins.bool: ...";
    assert_eq!(class.matches(eq).count(), 1);
    let class = class.replacen("class Element:\n", &format!("class Element:\n{constants}"), 1)
        .replace(eq, "    def __eq__(self, value: builtins.object, /) -> builtins.bool | types.NotImplementedType: ...\n    def __ne__(self, value: builtins.object, /) -> builtins.bool | types.NotImplementedType: ...");
    text.replace_range(start..end, &class);
    text = text.replacen("import typing\n", "import typing\nimport types\n", 1);
    text = text.replace("class DescriptorError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n", "class DescriptorError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    # Context fields exist only for applicable Rust variants.\n    function: builtins.str\n    field: builtins.str\n    actual: builtins.int\n    expected: builtins.int\n    minimum: builtins.int\n    expected_rows: builtins.int\n    actual_rows: typing.Optional[builtins.int]\n    include_sulfur_phosphorus: builtins.bool\n    contribs_len: builtins.int\n    bin_prop_len: builtins.int\n    bins_len: builtins.int\n    cell: builtins.str\n    row: builtins.int\n    detail: builtins.str\n");
    std::fs::write(path, text)?;
    Ok(())
}
