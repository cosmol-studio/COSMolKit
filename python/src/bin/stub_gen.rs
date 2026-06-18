use pyo3_stub_gen::Result;
use std::fs;
use std::io;
use std::path::Path;

fn main() -> Result<()> {
    let stub = cosmolkit::stub_info()?;
    stub.generate()?;

    let pyi_path = [
        Path::new("./cosmolkit.pyi"),
        Path::new("./python/cosmolkit.pyi"),
    ]
    .into_iter()
    .find(|path| path.exists())
    .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, "cosmolkit.pyi was not generated"))?;
    let mut text = fs::read_to_string(pyi_path)?;
    let future_line = "from __future__ import annotations\n";
    if !text.contains(future_line) {
        text = format!("{future_line}{text}");
    }
    text = expose_chem_enums(text);
    text = expose_residue_enums_and_functions(text);
    text = expose_batch_validation_error(text);
    text = expose_batch_error_mode_inputs(text);
    text = expose_batch_getitem_overloads(text);
    text = expose_sdf_dataset_getitem_overloads(text);
    text = expose_module_functions(text);
    text = expose_confseq_module(text);
    fs::write(pyi_path, text)?;

    Ok(())
}

fn expose_confseq_module(mut text: String) -> String {
    let export = "    \"confseq\",\n";
    if !text.contains(export) {
        text = text.replace(
            "    \"__version__\",\n",
            &format!("    \"__version__\",\n{export}"),
        );
    }

    if !text.contains("class _ConfSeqModule(typing.Protocol):") {
        let declarations = r#"class _ConfSeqModule(typing.Protocol):
    def decode(self, in_smiles: builtins.str, confseq: builtins.str, *, optimize_with_uff: builtins.bool = ...) -> Molecule: ...
    def decode_batch(self, in_smiles_list: builtins.list[builtins.str], confseq_list: builtins.list[builtins.str], *, errors: builtins.str = ..., n_jobs: typing.Optional[builtins.int] = ..., optimize_with_uff: builtins.bool = ...) -> builtins.list[typing.Optional[Molecule]]: ...

confseq: _ConfSeqModule

"#;
        text = text.replace(
            "__version__: builtins.str\n",
            &format!("{declarations}__version__: builtins.str\n"),
        );
    }
    text
}

fn expose_batch_error_mode_inputs(text: String) -> String {
    text.replace(
        "errors: typing.Optional[typing.Any]",
        "errors: typing.Union[BatchErrorMode, builtins.str, None]",
    )
}

fn expose_batch_getitem_overloads(mut text: String) -> String {
    let overloads = r#"@typing.overload
    def __getitem__(self, index: builtins.int) -> typing.Optional[Molecule]: ...
    @typing.overload
    def __getitem__(self, index: builtins.slice) -> MoleculeBatch: ...
    @typing.overload
    def __getitem__(self, index: typing.Sequence[builtins.bool]) -> MoleculeBatch: ...
    @typing.overload
    def __getitem__(self, index: typing.Sequence[builtins.int]) -> MoleculeBatch: ..."#;
    if text.contains("def __getitem__(self, key: typing.Any) -> typing.Any: ...") {
        text = text.replace(
            "def __getitem__(self, key: typing.Any) -> typing.Any: ...",
            overloads,
        );
    } else if text
        .contains("def __getitem__(self, index: builtins.int) -> typing.Optional[Molecule]: ...")
    {
        text = text.replace(
            "def __getitem__(self, index: builtins.int) -> typing.Optional[Molecule]: ...",
            overloads,
        );
    }
    text
}

fn expose_sdf_dataset_getitem_overloads(mut text: String) -> String {
    let overloads = r#"    @typing.overload
    def __getitem__(self, index: builtins.int) -> SdfRecord: ...
    @typing.overload
    def __getitem__(self, index: builtins.slice) -> MoleculeBatch: ...
    @typing.overload
    def __getitem__(self, index: typing.Sequence[builtins.bool]) -> MoleculeBatch: ...
    @typing.overload
    def __getitem__(self, index: typing.Sequence[builtins.int]) -> MoleculeBatch: ..."#;
    text = text.replace(
        "    def __getitem__(self, key: typing.Any) -> typing.Union[SdfRecord, MoleculeBatch]: ...",
        overloads,
    );
    text
}

fn expose_module_functions(mut text: String) -> String {
    for name in [
        "__version__",
        "version",
        "mol_to_binary",
        "mol_from_binary",
        "parse_smarts",
        "uff_has_all_molecule_params",
        "uff_optimize_molecule",
        "uff_optimize_molecule_confs",
        "mmff_has_all_molecule_params",
        "mmff_optimize_molecule",
        "mmff_optimize_molecule_confs",
        "has_substruct_match",
        "get_substruct_match",
        "get_substruct_matches",
        "get_substruct_matches_with_params",
    ] {
        let export = format!("    \"{name}\",\n");
        if !text.contains(&export) {
            text = text.replace(
                "    \"SubstructMatchResult\",\n",
                &format!("    \"SubstructMatchResult\",\n{export}"),
            );
        }
    }

    if !text.contains("__version__: builtins.str") {
        let declarations = r#"__version__: builtins.str
def version() -> builtins.str: ...
def mol_to_binary(mol: Molecule) -> bytes: ...
def mol_from_binary(data: bytes) -> Molecule: ...
def parse_smarts(smarts: builtins.str) -> SmartsMolecule: ...
def uff_has_all_molecule_params(mol: Molecule) -> builtins.bool: ...
def uff_optimize_molecule(mol: Molecule, max_iters: builtins.int = ..., vdw_thresh: builtins.float = ..., conf_id: builtins.int = ..., ignore_interfrag_interactions: builtins.bool = ...) -> UffOptimizeMoleculeResult: ...
def uff_optimize_molecule_confs(mol: Molecule, num_threads: builtins.int = ..., max_iters: builtins.int = ..., vdw_thresh: builtins.float = ..., ignore_interfrag_interactions: builtins.bool = ...) -> UffOptimizeMoleculeConfsResult: ...
def mmff_has_all_molecule_params(mol: Molecule) -> builtins.bool: ...
def mmff_optimize_molecule(mol: Molecule, mmff_variant: builtins.str = ..., max_iters: builtins.int = ..., non_bonded_thresh: builtins.float = ..., conf_id: builtins.int = ..., ignore_interfrag_interactions: builtins.bool = ...) -> MmffOptimizeMoleculeResult: ...
def mmff_optimize_molecule_confs(mol: Molecule, num_threads: builtins.int = ..., max_iters: builtins.int = ..., mmff_variant: builtins.str = ..., non_bonded_thresh: builtins.float = ..., ignore_interfrag_interactions: builtins.bool = ...) -> MmffOptimizeMoleculeConfsResult: ...
def has_substruct_match(mol: Molecule, query: Molecule) -> builtins.bool: ...
def get_substruct_match(mol: Molecule, query: Molecule) -> typing.Optional[SubstructMatchResult]: ...
def get_substruct_matches(mol: Molecule, query: Molecule) -> builtins.list[SubstructMatchResult]: ...
def get_substruct_matches_with_params(mol: Molecule, query: Molecule, max_matches: builtins.int = ..., uniquify: builtins.bool = ...) -> builtins.list[SubstructMatchResult]: ...

"#;
        text = text.replace(
            "@typing.final\nclass SubstructMatchResult:",
            &format!("{declarations}@typing.final\nclass SubstructMatchResult:"),
        );
    }
    text
}

fn residue_stub_defs() -> String {
    let source = fs::read_to_string("crates/cosmolkit-core/src/bio/resinfo.rs")
        .or_else(|_| fs::read_to_string("../crates/cosmolkit-core/src/bio/resinfo.rs"))
        .expect("read Rust ResidueCode enum for Python stub generation");
    let start = source
        .find("pub enum ResidueCode {")
        .expect("find ResidueCode enum");
    let rest = &source[start..];
    let end = rest
        .find("\n}\n\nimpl ResidueCode")
        .expect("find ResidueCode enum end");
    let mut out = String::from("@typing.final\nclass ResidueCode(enum.IntEnum):\n");
    for raw_line in rest[..end].lines().skip(1) {
        let line = raw_line.trim().trim_end_matches(',');
        if line.is_empty() {
            continue;
        }
        let Some((name, value)) = line.split_once('=') else {
            continue;
        };
        out.push_str("    ");
        out.push_str(name.trim());
        out.push_str(" = ");
        out.push_str(value.trim());
        out.push('\n');
    }
    out.push_str("\n@typing.final\nclass ResidueInfoKind(enum.IntEnum):\n");
    for (name, value) in [
        ("UNKNOWN", 0),
        ("AA", 1),
        ("AAD", 2),
        ("PAA", 3),
        ("MAA", 4),
        ("RNA", 5),
        ("DNA", 6),
        ("BUF", 7),
        ("HOH", 8),
        ("PYR", 9),
        ("KET", 10),
        ("ELS", 11),
    ] {
        out.push_str(&format!("    {name} = {value}\n"));
    }
    out.push_str("\nRESIDUE_CODE_MAP: typing.Mapping[builtins.str, ResidueCode]\n");
    out.push_str("RESIDUE_INFO_KIND_MAP: typing.Mapping[builtins.str, ResidueInfoKind]\n\n");
    out
}

fn expose_residue_enums_and_functions(mut text: String) -> String {
    if !text.contains("import enum\n") {
        text = text.replace("import builtins\n", "import builtins\nimport enum\n");
    }

    for name in [
        "ResidueCode",
        "ResidueInfoKind",
        "RESIDUE_CODE_MAP",
        "RESIDUE_INFO_KIND_MAP",
        "ResidueInfo",
        "find_tabulated_residue",
        "find_tabulated_residue_idx",
        "get_residue_info",
        "residue_code_from_name",
        "expand_one_letter",
        "expand_protein_one_letter",
        "expand_one_letter_sequence",
        "expand_protein_one_letter_string",
    ] {
        let export = format!("    \"{name}\",\n");
        if !text.contains(&export) {
            text = text.replace("    \"Protein\",\n", &format!("{export}    \"Protein\",\n"));
        }
    }

    let residue_defs = residue_stub_defs();
    if !text.contains("class ResidueCode(enum.IntEnum)") {
        text = text.replace(
            "@typing.final\nclass BondOrder(enum.IntEnum):",
            &format!("{residue_defs}@typing.final\nclass BondOrder(enum.IntEnum):"),
        );
    }

    if !text.contains("def find_tabulated_residue(name: builtins.str) -> ResidueInfo") {
        let functions = r#"def find_tabulated_residue(name: builtins.str) -> ResidueInfo: ...
def find_tabulated_residue_idx(name: builtins.str) -> builtins.int: ...
def get_residue_info(idx: builtins.int) -> ResidueInfo: ...
def residue_code_from_name(name: builtins.str) -> ResidueCode: ...
def expand_one_letter(code: builtins.str, kind: ResidueInfoKind) -> typing.Optional[builtins.str]: ...
def expand_protein_one_letter(code: builtins.str) -> typing.Optional[builtins.str]: ...
def expand_one_letter_sequence(seq: builtins.str, kind: ResidueInfoKind) -> builtins.list[builtins.str]: ...
def expand_protein_one_letter_string(seq: builtins.str) -> builtins.list[builtins.str]: ...

"#;
        text = text.replace(
            "def version() -> builtins.str: ...\n",
            &format!("def version() -> builtins.str: ...\n{functions}"),
        );
    }

    text = text.replace("def expand_one_letter(code: builtins.str, kind: builtins.int) -> typing.Optional[builtins.str]: ...", "def expand_one_letter(code: builtins.str, kind: ResidueInfoKind) -> typing.Optional[builtins.str]: ...");
    text = text.replace("def expand_one_letter_sequence(seq: builtins.str, kind: builtins.int) -> builtins.list[builtins.str]: ...", "def expand_one_letter_sequence(seq: builtins.str, kind: ResidueInfoKind) -> builtins.list[builtins.str]: ...");
    text
}

fn expose_chem_enums(mut text: String) -> String {
    if !text.contains("import enum\n") {
        text = text.replace("import builtins\n", "import builtins\nimport enum\n");
    }

    for name in [
        "BondOrder",
        "BondDirection",
        "BondStereo",
        "ChiralTag",
        "BatchErrorMode",
        "BOND_ORDER_MAP",
        "BOND_DIRECTION_MAP",
        "BOND_STEREO_MAP",
        "CHIRAL_TAG_MAP",
        "BATCH_ERROR_MODE_MAP",
    ] {
        let export = format!("    \"{name}\",\n");
        if !text.contains(&export) {
            text = text.replace("    \"Bond\",\n", &format!("{export}    \"Bond\",\n"));
        }
    }

    let enum_defs = r#"@typing.final
class BondOrder(enum.IntEnum):
    UNSPECIFIED = 0
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    QUADRUPLE = 4
    QUINTUPLE = 5
    HEXTUPLE = 6
    ONEANDAHALF = 7
    TWOANDAHALF = 8
    THREEANDAHALF = 9
    FOURANDAHALF = 10
    FIVEANDAHALF = 11
    AROMATIC = 12
    IONIC = 13
    DATIVE = 14
    DATIVEONE = 15
    DATIVEL = 16
    DATIVER = 17
    HYDROGEN = 18
    THREECENTER = 19
    OTHER = 20
    ZERO = 21

@typing.final
class BondDirection(enum.IntEnum):
    NONE = 0
    BEGINWEDGE = 1
    BEGINDASH = 2
    ENDUPRIGHT = 3
    ENDDOWNRIGHT = 4
    EITHERDOUBLE = 5
    UNKNOWN = 6

@typing.final
class BondStereo(enum.IntEnum):
    NONE = 0
    ANY = 1
    Z = 2
    E = 3
    CIS = 4
    TRANS = 5
    ATROP_CW = 6
    ATROP_CCW = 7

@typing.final
class ChiralTag(enum.IntEnum):
    CHI_UNSPECIFIED = 0
    CHI_TETRAHEDRAL_CW = 1
    CHI_TETRAHEDRAL_CCW = 2
    CHI_OTHER = 3
    CHI_TETRAHEDRAL = 4
    CHI_ALLENE = 5
    CHI_SQUAREPLANAR = 6
    CHI_TRIGONALBIPYRAMIDAL = 7
    CHI_OCTAHEDRAL = 8

@typing.final
class BatchErrorMode(enum.IntEnum):
    RAISE = 1
    KEEP = 2

BOND_ORDER_MAP: typing.Mapping[builtins.str, BondOrder]
BOND_DIRECTION_MAP: typing.Mapping[builtins.str, BondDirection]
BOND_STEREO_MAP: typing.Mapping[builtins.str, BondStereo]
CHIRAL_TAG_MAP: typing.Mapping[builtins.str, ChiralTag]
BATCH_ERROR_MODE_MAP: typing.Mapping[builtins.str, BatchErrorMode]

"#;
    let existing_enum_defs = r#"@typing.final
class BondOrder(enum.IntEnum):
    UNSPECIFIED = 0
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    QUADRUPLE = 4
    QUINTUPLE = 5
    HEXTUPLE = 6
    ONEANDAHALF = 7
    TWOANDAHALF = 8
    THREEANDAHALF = 9
    FOURANDAHALF = 10
    FIVEANDAHALF = 11
    AROMATIC = 12
    IONIC = 13
    DATIVE = 14
    DATIVEONE = 15
    DATIVEL = 16
    DATIVER = 17
    HYDROGEN = 18
    THREECENTER = 19
    OTHER = 20
    ZERO = 21

@typing.final
class BondDirection(enum.IntEnum):
    NONE = 0
    BEGINWEDGE = 1
    BEGINDASH = 2
    ENDUPRIGHT = 3
    ENDDOWNRIGHT = 4
    EITHERDOUBLE = 5
    UNKNOWN = 6

@typing.final
class BondStereo(enum.IntEnum):
    NONE = 0
    ANY = 1
    Z = 2
    E = 3
    CIS = 4
    TRANS = 5
    ATROP_CW = 6
    ATROP_CCW = 7

@typing.final
class ChiralTag(enum.IntEnum):
    CHI_UNSPECIFIED = 0
    CHI_TETRAHEDRAL_CW = 1
    CHI_TETRAHEDRAL_CCW = 2
    CHI_OTHER = 3
    CHI_TETRAHEDRAL = 4
    CHI_ALLENE = 5
    CHI_SQUAREPLANAR = 6
    CHI_TRIGONALBIPYRAMIDAL = 7
    CHI_OCTAHEDRAL = 8

@typing.final
class BatchErrorMode(enum.IntEnum):
    RAISE = 1
    KEEP = 2

BOND_ORDER_MAP: typing.Mapping[builtins.str, BondOrder]
BOND_DIRECTION_MAP: typing.Mapping[builtins.str, BondDirection]
BOND_STEREO_MAP: typing.Mapping[builtins.str, BondStereo]
CHIRAL_TAG_MAP: typing.Mapping[builtins.str, ChiralTag]
BATCH_ERROR_MODE_MAP: typing.Mapping[builtins.str, BatchErrorMode]

"#;
    if text.contains(existing_enum_defs) {
        text = text.replace(existing_enum_defs, enum_defs);
    } else if !text.contains("class BondOrder(enum.IntEnum)") {
        text = text.replace(
            "@typing.final\nclass Atom:",
            &format!("{enum_defs}@typing.final\nclass Atom:"),
        );
    } else {
        text = replace_existing_chem_enum_block(text, enum_defs);
    }

    text
}

fn replace_existing_chem_enum_block(text: String, enum_defs: &str) -> String {
    let Some(start) = text.find("@typing.final\nclass BondOrder(enum.IntEnum):") else {
        return text;
    };
    let Some(after_start) = text[start..].find("@typing.final\nclass Atom:") else {
        return text;
    };
    let end = start + after_start;
    format!("{}{}{}", &text[..start], enum_defs, &text[end..])
}

fn expose_batch_validation_error(mut text: String) -> String {
    let export = "    \"BatchValidationError\",\n";
    if !text.contains(export) {
        text = text.replace(
            "    \"BatchExportReport\",\n",
            &format!("    \"BatchExportReport\",\n{export}"),
        );
    }

    let class_decl = "class BatchValidationError(builtins.ValueError):\n    def errors(self) -> builtins.list[BatchError]: ...\n\n";
    if !text.contains("class BatchValidationError") {
        if text.contains("@typing.final\nclass Bond:") {
            text = text.replace(
                "@typing.final\nclass Bond:",
                &format!("{class_decl}@typing.final\nclass Bond:"),
            );
        } else if text.contains("@typing.final\nclass BatchExportReport:") {
            text = text.replace(
                "@typing.final\nclass BatchExportReport:",
                &format!("{class_decl}@typing.final\nclass BatchExportReport:"),
            );
        }
    } else {
        text = text.replace(
            "class BatchValidationError(builtins.ValueError):\n    ...\n\n",
            class_decl,
        );
    }

    text
}
