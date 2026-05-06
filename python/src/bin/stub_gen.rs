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
    text = expose_batch_validation_error(text);
    text = expose_batch_error_mode_inputs(text);
    text = expose_batch_getitem_overloads(text);
    fs::write(pyi_path, text)?;

    Ok(())
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
        "BatchErrorType",
        "BOND_ORDER_MAP",
        "BOND_DIRECTION_MAP",
        "BOND_STEREO_MAP",
        "CHIRAL_TAG_MAP",
        "BATCH_ERROR_MODE_MAP",
        "BATCH_ERROR_TYPE_MAP",
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
    AROMATIC = 5
    DATIVE = 6
    HYDROGEN = 7

@typing.final
class BondDirection(enum.IntEnum):
    NONE = 0
    ENDUPRIGHT = 1
    ENDDOWNRIGHT = 2
    UNKNOWN = 3

@typing.final
class BondStereo(enum.IntEnum):
    STEREONONE = 0
    STEREOANY = 1
    STEREOCIS = 2
    STEREOTRANS = 3

@typing.final
class ChiralTag(enum.IntEnum):
    CHI_UNSPECIFIED = 0
    CHI_TETRAHEDRAL_CW = 1
    CHI_TETRAHEDRAL_CCW = 2
    CHI_TRIGONALBIPYRAMIDAL = 3

@typing.final
class BatchErrorMode(enum.IntEnum):
    RAISE = 1
    KEEP = 2
    SKIP = 3

@typing.final
class BatchErrorType(enum.IntEnum):
    UNKNOWN = 0
    SMILES_PARSE = 1
    SDF_READ = 2
    ADD_HYDROGENS = 3
    REMOVE_HYDROGENS = 4
    SANITIZE = 5
    KEKULIZE = 6
    COORDINATE_GENERATION = 7
    SMILES_WRITE = 8
    DISTANCE_GEOMETRY = 9
    FINGERPRINT = 10
    SVG_DRAW = 11
    PREPARED_DRAW = 12
    SDF_WRITE = 13
    IMAGE_EXPORT = 14
    IO = 15
    FILENAME = 16

BOND_ORDER_MAP: typing.Mapping[builtins.str, BondOrder]
BOND_DIRECTION_MAP: typing.Mapping[builtins.str, BondDirection]
BOND_STEREO_MAP: typing.Mapping[builtins.str, BondStereo]
CHIRAL_TAG_MAP: typing.Mapping[builtins.str, ChiralTag]
BATCH_ERROR_MODE_MAP: typing.Mapping[builtins.str, BatchErrorMode]
BATCH_ERROR_TYPE_MAP: typing.Mapping[builtins.str, BatchErrorType]

"#;
    let existing_enum_defs = r#"@typing.final
class BondOrder(enum.IntEnum):
    UNSPECIFIED = 0
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    QUADRUPLE = 4
    AROMATIC = 5
    DATIVE = 6

@typing.final
class BondDirection(enum.IntEnum):
    NONE = 0
    ENDUPRIGHT = 1
    ENDDOWNRIGHT = 2

@typing.final
class BondStereo(enum.IntEnum):
    STEREONONE = 0
    STEREOANY = 1
    STEREOCIS = 2
    STEREOTRANS = 3

@typing.final
class ChiralTag(enum.IntEnum):
    CHI_UNSPECIFIED = 0
    CHI_TETRAHEDRAL_CW = 1
    CHI_TETRAHEDRAL_CCW = 2
    CHI_TRIGONALBIPYRAMIDAL = 3

@typing.final
class BatchErrorMode(enum.IntEnum):
    RAISE = 1
    KEEP = 2
    SKIP = 3

@typing.final
class BatchErrorType(enum.IntEnum):
    UNKNOWN = 0
    SMILES_PARSE = 1
    SDF_READ = 2
    ADD_HYDROGENS = 3
    REMOVE_HYDROGENS = 4
    SANITIZE = 5
    KEKULIZE = 6
    COORDINATE_GENERATION = 7
    SMILES_WRITE = 8
    DISTANCE_GEOMETRY = 9
    FINGERPRINT = 10
    SVG_DRAW = 11
    PREPARED_DRAW = 12
    SDF_WRITE = 13
    IMAGE_EXPORT = 14
    IO = 15
    FILENAME = 16

BOND_ORDER_MAP: typing.Mapping[builtins.str, BondOrder]
BOND_DIRECTION_MAP: typing.Mapping[builtins.str, BondDirection]
BOND_STEREO_MAP: typing.Mapping[builtins.str, BondStereo]
CHIRAL_TAG_MAP: typing.Mapping[builtins.str, ChiralTag]
BATCH_ERROR_MODE_MAP: typing.Mapping[builtins.str, BatchErrorMode]
BATCH_ERROR_TYPE_MAP: typing.Mapping[builtins.str, BatchErrorType]

"#;
    if text.contains(existing_enum_defs) {
        text = text.replace(existing_enum_defs, enum_defs);
    } else if !text.contains("class BondOrder(enum.IntEnum)") {
        text = text.replace(
            "@typing.final\nclass Atom:",
            &format!("{enum_defs}@typing.final\nclass Atom:"),
        );
    }

    text
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
