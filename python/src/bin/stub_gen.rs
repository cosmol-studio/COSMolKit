fn main() -> pyo3_stub_gen::Result<()> {
    // Render in memory: a failed contract check must not replace the last stub.
    let info = cosmolkit_py::stub_info()?;
    if info.modules.len() != 1 || !info.modules.contains_key("cosmolkit") {
        return Err(
            std::io::Error::other("expected the canonical flat cosmolkit stub module").into(),
        );
    }
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("cosmolkit.pyi");
    let mut text = info.modules["cosmolkit"].format_with_config(info.config.use_type_statement);
    // PyO3 eq_int class enums expose typed variants and integer conversion,
    // not enum.Enum.name/value. Keep the generated TAU/property stubs faithful
    // to their actual native classes rather than promising nonexistent fields.
    for name in [
        "TautomerEnumerationStatus",
        "CoordinateZPolicy",
        "PropertyValueKind",
        "SdfPropertyListTarget",
        "ValenceModel",
        "AromaticityModel",
        "SanitizeStage",
        "ParityPolicy",
        "CoordinateDimension",
        "SdfCoordinateMode",
        "BioConnectionKind",
        "BioAsu",
        "BioHelixClass",
        "BioSoftwareClassification",
        "BioAssemblySpecialKind",
        "BioPdbReadStage",
        "BioMmcifReadStage",
        "UffParameterErrorKind",
    ] {
        let prefix = format!("class {name}(enum.Enum):\n");
        assert_eq!(
            text.matches(&prefix).count(),
            1,
            "generated {name} schema changed"
        );
        let start = text.find(&prefix).unwrap();
        let end = start + text[start..].find("\n\n").expect("generated enum boundary");
        let members = text[start + prefix.len()..end]
            .lines()
            .map(|line| {
                line.strip_prefix("    ")
                    .and_then(|member| member.strip_suffix(" = ..."))
                    .expect("generated PyO3 enum variant")
                    .to_owned()
            })
            .collect::<Vec<_>>();
        assert!(!members.is_empty(), "generated {name} has no variants");
        let replacement = format!(
            "class {name}:\n{}    def __int__(self) -> builtins.int: ...\n    def __eq__(self, other: builtins.object, /) -> builtins.bool | types.NotImplementedType: ...\n    def __ne__(self, other: builtins.object, /) -> builtins.bool | types.NotImplementedType: ...",
            members
                .iter()
                .map(|member| format!("    {member}: typing.ClassVar[{name}]\n"))
                .collect::<String>()
        );
        text.replace_range(start..end, &replacement);
    }
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
        "PickleError",
        "MoleculeHashError",
        "CipRankError",
        "SmilesError",
        "SmartsParseError",
        "SmilesWriteError",
        "MorganReadError",
        "AtomPairReadError",
        "TopologicalTorsionReadError",
        "TopologicalTorsionPathScoreError",
        "PatternFingerprintError",
        "FingerprintPreparationError",
        "FingerprintError",
        "TopologicalFingerprintError",
        "FingerprintJsonError",
        "DescriptorReadError",
        "DescriptorError",
        "TautomerRunError",
        "TautomerCatalogError",
        "PropertyValueError",
        "ValenceError",
        "CipDescriptorError",
        "MmffMolPropertiesError",
        "MmffOptimizationError",
        "UffOptimizationError",
        "UffParameterQueryError",
        "StereoReadError",
        "CoordinateInputError",
        "Coordinate3DReadError",
        "AlignmentError",
        "SanitizeError",
        "MatrixError",
        "ChemistryProblemError",
        "KekulizeError",
        "SdfError",
    ] {
        text.push_str(&format!("\nclass {name}(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n"));
        text = text.replace("__all__ = [\n", &format!("__all__ = [\n    \"{name}\",\n"));
    }
    for (name, fields) in [
        (
            "SdfError",
            "    # WrongGraphKind retains the exact expected and actual payload tags.\n    expected: builtins.str\n    actual: builtins.str\n",
        ),
        (
            "SanitizeError",
            "    bits: builtins.int\n    unknown_bits: builtins.int\n    stage: SanitizeStage\n",
        ),
        (
            "MatrixError",
            "    position: builtins.int\n    atom: builtins.int\n    atom_count: builtins.int\n    bond: builtins.int\n    bond_count: builtins.int\n    first_position: builtins.int\n    second_position: builtins.int\n    endpoint: builtins.str\n    order: BondOrder\n    conformer_id: builtins.int\n    dimension: builtins.int\n",
        ),
        (
            "ValenceError",
            "    atom: builtins.int\n    atomic_number: builtins.int\n    formal_charge: builtins.int\n    phase: builtins.str\n    calculated: typing.Optional[builtins.int]\n    reason: builtins.str\n    message: builtins.str\n    atom_count: builtins.int\n    neighbor_atom: builtins.int\n    bond: typing.Optional[builtins.int]\n    bond_count: builtins.int\n    begin: builtins.int\n    end: builtins.int\n    value: builtins.int\n    field: builtins.str\n    explicit_valence: builtins.int\n    physical_bonds: builtins.int\n    explicit: builtins.int\n    implicit: builtins.int\n    neighbor_hydrogens: builtins.int\n    order: BondOrder\n",
        ),
        (
            "KekulizeError",
            "    expected: builtins.int\n    actual: builtins.int\n    atom: builtins.int\n    atom_count: builtins.int\n    field: builtins.str\n    begin: builtins.int\n    end: builtins.int\n    questions: builtins.int\n    bit_width: builtins.int\n    problem_atoms: typing.List[builtins.int]\n    before: builtins.int\n    after: builtins.int\n    bond: builtins.int\n    detail: builtins.str\n",
        ),
    ] {
        let prefix = format!(
            "class {name}(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n"
        );
        text = text.replace(
            &prefix,
            &format!("{prefix}    # Fields exist only on applicable Rust variants.\n{fields}"),
        );
    }
    // Pattern publishes a ValueError subclass; project the attributes set by
    // canonical_pattern::pattern_pyerr, with variant-specific context fields.
    text = text.replace(
        "class PatternFingerprintError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n",
        "class PatternFingerprintError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    # Context fields exist only on applicable Rust variants.\n    reason: builtins.str\n    left: builtins.int\n    right: builtins.int\n",
    );
    text = text.replace("class TopologicalFingerprintError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n", "class TopologicalFingerprintError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    # Present only on applicable variants.\n    reason: builtins.str\n    field: builtins.str\n");
    // Alignment publishes one create_exception! type. Context attributes match
    // alignment_binding::alignment_pyerr and exist only on applicable variants.
    text = text.replace(
        "class AlignmentError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n",
        "class AlignmentError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    # Context fields exist only on applicable Rust variants.\n    id: builtins.int\n    index: builtins.int\n    atom_count: builtins.int\n    map_len: builtins.int\n    weight_len: builtins.int\n    message: builtins.str\n",
    );
    text = text.replace("class PickleError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n", "class PickleError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    # Context fields exist only on applicable Rust variants.\n    version: builtins.int\n    major: builtins.int\n    minor: builtins.int\n    section: builtins.int\n    expected: builtins.int\n    actual: builtins.int\n    value: builtins.int\n    type_name: builtins.str\n    count: builtins.int\n    message: builtins.str\n");
    text = text.replace("class MoleculeHashError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n", "class MoleculeHashError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    # Context fields exist only on applicable Rust variants.\n    actual: builtins.int\n    atom_count: builtins.int\n");
    text = text.replace("class CipRankError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n", "class CipRankError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    # Context fields exist only on applicable Rust variants.\n    field: builtins.str\n    actual: builtins.int\n    atom_count: builtins.int\n    atom: builtins.int\n    value: builtins.int\n    map_number: builtins.int\n    degree: builtins.int\n    maximum_supported: builtins.int\n    bond: builtins.int\n    order: builtins.int\n");
    text = text.replace("class FingerprintError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n", "class FingerprintError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    # Context fields exist only on applicable variants.\n    index: builtins.int\n    size: builtins.int\n    left: builtins.int\n    right: builtins.int\n    factor: builtins.int\n    n_bits: builtins.int\n    value: builtins.float\n    site: builtins.str\n    what: builtins.str\n    reason: builtins.str\n");
    text = text.replace("class SmartsParseError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n", "class SmartsParseError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    # Context fields exist only on applicable Rust variants.\n    position: builtins.int\n    character: builtins.str\n    context: builtins.str\n    detail: builtins.str\n    ring: builtins.int\n    atom: builtins.int\n    begin_atom: builtins.int\n    end_atom: builtins.int\n    feature: builtins.str\n    carrier: builtins.int\n");
    text = text.replace(
        "class TopologicalTorsionPathScoreError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n",
        "class TopologicalTorsionPathScoreError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    # Context fields exist only on applicable Rust variants.\n    actual: builtins.int\n    required: builtins.int\n    index: builtins.int\n    atom_count: builtins.int\n    code: builtins.int\n    subtract: builtins.int\n",
    );
    text = text.replace("class PropertyValueError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n", "class PropertyValueError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    def expected(self) -> PropertyValueKind: ...\n    def actual(self) -> PropertyValueKind: ...\n");
    // Forcefield exceptions are published through create_exception!, so they
    // have no pyclass stub metadata. Project only attributes set by the thin
    // native converters; parameter causes expose kind without domain.
    text = text.replace(
        "class UffOptimizationError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n",
        "class UffOptimizationError(builtins.ValueError):\n    domain: builtins.str\n    def kind(self) -> UffOptimizationErrorKind: ...\n    requested: typing.Optional[builtins.int]\n",
    );
    text.push_str("\nclass UffParameterError(builtins.ValueError):\n    def kind(self) -> UffParameterErrorKind: ...\n");
    text = text.replace("__all__ = [\n", "__all__ = [\n    \"UffParameterError\",\n");
    text=text.replace("class CoordinateInputError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n", "class CoordinateInputError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    # Fields exist only on applicable Rust variants.\n    dimension: builtins.str\n    expected: builtins.int\n    actual: builtins.int\n    row: builtins.int\n    columns: builtins.int\n    expected_columns: builtins.str\n    column: builtins.int\n    value: builtins.float | builtins.str\n    z: builtins.float\n    conformer_id: builtins.int\n    count: builtins.int\n    max_id: builtins.int\n    conformer: builtins.int\n    axis: builtins.str\n");
    text=text.replace("class Coordinate3DReadError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n", "class Coordinate3DReadError(builtins.ValueError):\n    domain: builtins.str\n    kind: builtins.str\n    conformer_id: builtins.int\n    count: builtins.int\n");
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
    // Dynamic Python enums use exactly the canonical vocabulary used by register().
    fn append_enum(text: &mut String, name: &str, base: &str, members: Vec<(String, String)>) {
        assert!(
            !text.contains(&format!("class {name}(")),
            "duplicate generated enum {name}"
        );
        if !text.contains("import enum\n") {
            text.insert_str(0, "import enum\n");
        }
        text.push_str(&format!("\nclass {name}({base}):\n"));
        for (member, value) in members {
            text.push_str(&format!("    {member} = {value}\n"));
        }
        *text = text.replace("__all__ = [\n", &format!("__all__ = [\n    \"{name}\",\n"));
    }
    append_enum(
        &mut text,
        "BondOrder",
        "enum.IntEnum",
        (0..=::cosmolkit::BondOrder::Zero.rdkit_code())
            .map(|code| {
                let value =
                    ::cosmolkit::BondOrder::from_rdkit_code(code).expect("declared source code");
                (value.rdkit_name().into(), value.rdkit_code().to_string())
            })
            .collect(),
    );
    append_enum(
        &mut text,
        "ChiralTag",
        "enum.IntEnum",
        (0..=::cosmolkit::ChiralTag::Octahedral.rdkit_code())
            .map(|code| {
                let value =
                    ::cosmolkit::ChiralTag::from_rdkit_code(code).expect("declared source code");
                (value.rdkit_name().into(), value.rdkit_code().to_string())
            })
            .collect(),
    );
    append_enum(
        &mut text,
        "BondDirection",
        "enum.IntEnum",
        (0..=::cosmolkit::BondDirection::Unknown.rdkit_code())
            .map(|code| {
                let value = ::cosmolkit::BondDirection::from_rdkit_code(code)
                    .expect("declared source code");
                (value.rdkit_name().into(), value.rdkit_code().to_string())
            })
            .collect(),
    );
    append_enum(
        &mut text,
        "BondStereo",
        "enum.IntEnum",
        (0..=::cosmolkit::BondStereo::AtropCcw.rdkit_code())
            .map(|code| {
                let value =
                    ::cosmolkit::BondStereo::from_rdkit_code(code).expect("declared source code");
                (value.rdkit_name().into(), value.rdkit_code().to_string())
            })
            .collect(),
    );
    append_enum(
        &mut text,
        "Hybridization",
        "enum.IntEnum",
        (0..=::cosmolkit::Hybridization::Other.rdkit_code())
            .map(|code| {
                let value = ::cosmolkit::Hybridization::from_rdkit_code(code)
                    .expect("declared source code");
                (value.rdkit_name().into(), value.rdkit_code().to_string())
            })
            .collect(),
    );
    append_enum(
        &mut text,
        "CipDescriptor",
        "builtins.str, enum.Enum",
        [
            ::cosmolkit::CipDescriptor::R,
            ::cosmolkit::CipDescriptor::S,
            ::cosmolkit::CipDescriptor::LowerR,
            ::cosmolkit::CipDescriptor::LowerS,
            ::cosmolkit::CipDescriptor::E,
            ::cosmolkit::CipDescriptor::Z,
            ::cosmolkit::CipDescriptor::LowerE,
            ::cosmolkit::CipDescriptor::LowerZ,
            ::cosmolkit::CipDescriptor::M,
            ::cosmolkit::CipDescriptor::P,
            ::cosmolkit::CipDescriptor::LowerM,
            ::cosmolkit::CipDescriptor::LowerP,
        ]
        .into_iter()
        .map(|v| (v.as_str().into(), format!("{:?}", v.as_str())))
        .collect(),
    );
    text.push_str("\nclass TetrahedralStereo(typing.NamedTuple):\n    center: builtins.int\n    ligands: builtins.list[typing.Optional[builtins.int]]\n");
    text = text.replace("__all__ = [\n", "__all__ = [\n    \"TetrahedralStereo\",\n");
    let original = "def tetrahedral_stereo(self) -> builtins.list[typing.Any]:";
    assert_eq!(text.matches(original).count(), 1);
    text = text.replace(
        original,
        "def tetrahedral_stereo(self) -> builtins.list[TetrahedralStereo]:",
    );
    let mut text = expose_bio_types(text);
    text.push_str("\nclass AtomCodeExplanationError(builtins.KeyError):\n    domain: builtins.str\n    kind: builtins.str\n    code: builtins.int\n");
    text = text.replace(
        "__all__ = [\n",
        "__all__ = [\n    \"AtomCodeExplanationError\",\n",
    );
    check_registered_python_callables(&text)?;
    std::fs::write(path, text)?;
    Ok(())
}

fn check_registered_python_callables(text: &str) -> pyo3_stub_gen::Result<()> {
    use ::cosmolkit::{BINDING_CONTRACT, BindingItem, BindingOwner};
    use pyo3::{prelude::*, types::PyModule};

    // The linked registry already applies the facade's actual cfg gates.
    // Neither Experimental status nor a missing binding exempts an enabled row.
    let entries = BINDING_CONTRACT
        .iter()
        .map(|entry| {
            serde_json::json!({
                "semantic_id": entry.semantic_id,
                "python_name": entry.python_name,
                "feature": entry.feature,
                "item": match entry.item {
                    BindingItem::Callable => "callable",
                    BindingItem::Type => "type",
                },
                "owner": match entry.owner {
                    BindingOwner::Module => "module",
                    BindingOwner::Molecule => "molecule",
                    BindingOwner::Type => "type",
                },
            })
        })
        .collect::<Vec<_>>();
    let contract = serde_json::to_string(&entries)?;
    let code = std::ffi::CString::new(include_str!(
        "../../../dev/tools/check_python_stub_contract.py"
    ))?;
    Python::initialize();
    let missing = Python::attach(|py| -> PyResult<Vec<String>> {
        let checker = PyModule::from_code(
            py,
            &code,
            c"check_python_stub_contract.py",
            c"_stub_contract",
        )?;
        checker
            .getattr("check_contract")?
            .call1((text, contract))?
            .extract()
    })?;
    if !missing.is_empty() {
        return Err(std::io::Error::other(format!(
            "Python binding contract failed: {} registered callables missing from generated stubs:\n{}\nImplement the registered Python projections; no stub was written.",
            missing.len(),
            missing.join("\n"),
        )).into());
    }
    eprintln!(
        "Python binding contract: all {} enabled registered callables present",
        BINDING_CONTRACT
            .iter()
            .filter(|entry| entry.item == BindingItem::Callable)
            .count()
    );
    Ok(())
}

// Generated declarations use the same canonical owners as native publication.
fn expose_bio_types(mut text: String) -> String {
    text = text.replacen("import builtins\n", "import builtins\nimport enum\n", 1);
    let mut definitions = String::from("\nclass ResidueCode(enum.IntEnum):\n");
    let mut index = 0;
    while let Some(info) = ::cosmolkit::residue_info_checked(index) {
        definitions.push_str(&format!("    {:?} = {}\n", info.code, info.code.as_u16()));
        index += 1;
    }
    definitions.push_str("\nclass ResidueInfoKind(enum.IntEnum):\n");
    use ::cosmolkit::ResidueInfoKind as K;
    for kind in [
        K::Unknown,
        K::Aa,
        K::Aad,
        K::Paa,
        K::Maa,
        K::Rna,
        K::Dna,
        K::Buf,
        K::Hoh,
        K::Pyr,
        K::Ket,
        K::Els,
    ] {
        definitions.push_str(&format!("    {} = {}\n", kind.name(), kind as u8));
    }
    definitions.push_str("    def name(self) -> builtins.str: ...\n");
    definitions.push_str("\nclass BioCoordinateFormat(enum.IntEnum):\n");
    for format in [
        ::cosmolkit::BioCoordinateFormat::Unknown,
        ::cosmolkit::BioCoordinateFormat::Detect,
        ::cosmolkit::BioCoordinateFormat::Pdb,
        ::cosmolkit::BioCoordinateFormat::Mmcif,
        ::cosmolkit::BioCoordinateFormat::Mmjson,
        ::cosmolkit::BioCoordinateFormat::ChemComp,
    ] {
        definitions.push_str(&format!("    {format:?} = {}\n", format as u8));
    }
    definitions.push_str("\nRESIDUE_CODE_MAP: typing.Mapping[builtins.str, ResidueCode]\nRESIDUE_INFO_KIND_MAP: typing.Mapping[builtins.str, ResidueInfoKind]\n");
    for (name, base) in [
        ("BioReadError", "builtins.ValueError"),
        ("BioStructureError", "builtins.ValueError"),
        ("BioPdbReadError", "BioReadError"),
        ("BioMmcifReadError", "BioReadError"),
        ("ProteinReadError", "BioReadError"),
        ("BioOperationError", "builtins.ValueError"),
        ("BioMmcifWriteError", "builtins.ValueError"),
        ("BioPdbWriteError", "builtins.ValueError"),
        ("BioSelectionParseError", "builtins.ValueError"),
        ("BioMoleculeError", "builtins.ValueError"),
        ("BioMoleculeConversionError", "builtins.ValueError"),
    ] {
        definitions.push_str(&format!(
            "\nclass {name}({base}):\n    domain: builtins.str\n    kind: builtins.str\n"
        ));
        if name == "BioPdbReadError" {
            definitions.push_str("    def stage(self) -> BioPdbReadStage: ...\n    def line_number(self) -> typing.Optional[builtins.int]: ...\n    def record_tag(self) -> typing.Optional[builtins.bytes]: ...\n");
        }
        if name == "BioMmcifReadError" {
            definitions.push_str("    def stage(self) -> BioMmcifReadStage: ...\n");
        }
        text = text.replace("__all__ = [\n", &format!("__all__ = [\n    \"{name}\",\n"));
    }
    definitions.push_str("\nclass ResidueCodeParseError(builtins.ValueError):\n    def __init__(self, input: builtins.str) -> None: ...\n    def input(self) -> builtins.str: ...\n");
    for name in [
        "ResidueCodeParseError",
        "ResidueCode",
        "ResidueInfoKind",
        "BioCoordinateFormat",
        "RESIDUE_CODE_MAP",
        "RESIDUE_INFO_KIND_MAP",
    ] {
        text = text.replace("__all__ = [\n", &format!("__all__ = [\n    \"{name}\",\n"));
    }
    text.push_str(&definitions);
    text
}
