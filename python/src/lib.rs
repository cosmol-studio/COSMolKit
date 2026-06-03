use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufReader, Write};
use std::path::PathBuf;

use cosmolkit_core::io::molblock::{self, SdfFormat};
use cosmolkit_core::io::sdf::SdfReader;
use cosmolkit_core::{BatchErrorMode, BatchRecord, BatchRecordError, SmilesWriteParams};
use numpy::PyArray2;
use pyo3::PyErr;
use pyo3::exceptions::{PyIndexError, PyNotImplementedError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{
    PyAny, PyBool, PyBytes, PyDict, PyIterator, PyList, PyMapping, PyMappingProxy, PySlice,
    PySliceMethods, PyType,
};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::define_stub_info_gatherer;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
#[cfg(not(feature = "stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;

fn parse_coordinate_mode(value: Option<&str>) -> PyResult<cosmolkit_core::SdfCoordinateMode> {
    let Some(value) = value else {
        return Ok(cosmolkit_core::SdfCoordinateMode::Preserve);
    };
    match value.to_ascii_lowercase().as_str() {
        "auto" => Ok(cosmolkit_core::SdfCoordinateMode::Preserve),
        "2d" => Ok(cosmolkit_core::SdfCoordinateMode::Require2D),
        "3d" => Ok(cosmolkit_core::SdfCoordinateMode::Require3D),
        _ => Err(PyValueError::new_err(format!(
            "unsupported coordinate_dim '{value}', expected one of: auto, 2d, 3d"
        ))),
    }
}

fn parse_batch_error_mode(errors: Option<&Bound<'_, PyAny>>) -> PyResult<BatchErrorMode> {
    let Some(errors) = errors else {
        return Ok(BatchErrorMode::Strict);
    };
    if let Ok(value) = errors.extract::<String>() {
        return match value.to_ascii_lowercase().as_str() {
            "raise" => Ok(BatchErrorMode::Strict),
            "keep" => Ok(BatchErrorMode::KeepErrors),
            _ => Err(PyValueError::new_err(format!(
                "unsupported errors mode '{value}', expected one of: raise, keep"
            ))),
        };
    }
    match errors.extract::<i64>()? {
        1 => Ok(BatchErrorMode::Strict),
        2 => Ok(BatchErrorMode::KeepErrors),
        value => Err(PyValueError::new_err(format!(
            "unsupported errors mode code {value}, expected BatchErrorMode.RAISE or KEEP"
        ))),
    }
}

#[allow(clippy::too_many_arguments)]
fn make_smiles_write_params(
    isomeric_smiles: bool,
    canonical: bool,
    kekule: bool,
    clean_stereo: bool,
    all_bonds_explicit: bool,
    all_hs_explicit: bool,
    include_dative_bonds: bool,
    ignore_atom_map_numbers: bool,
    rooted_at_atom: Option<usize>,
) -> SmilesWriteParams {
    SmilesWriteParams {
        do_isomeric_smiles: isomeric_smiles,
        canonical,
        do_kekule: kekule,
        clean_stereo,
        all_bonds_explicit,
        all_hydrogens_explicit: all_hs_explicit,
        include_dative_bonds,
        ignore_atom_map_numbers,
        rooted_at_atom,
        ..Default::default()
    }
}

fn validate_n_jobs(n_jobs: Option<usize>) -> PyResult<Option<usize>> {
    if matches!(n_jobs, Some(0)) {
        return Err(PyValueError::new_err("n_jobs must be >= 1"));
    }
    Ok(n_jobs)
}

fn maybe_batch_progress_bar(
    enabled: bool,
    total: usize,
    message: impl Into<String>,
) -> Option<cosmolkit_core::BatchProgressBar> {
    enabled.then(|| cosmolkit_core::batch_progress_bar(total, message))
}

fn batch_validation_pyerr(error: cosmolkit_core::BatchValidationError) -> PyErr {
    let message = error.to_string();
    Python::attach(|py| {
        let error_count = error.errors;
        let reason = error.reason.map(|value| value.to_string());
        let record_errors: Vec<PyBatchError> =
            error.record_errors.into_iter().map(Into::into).collect();
        match (|| -> PyResult<Bound<'_, PyAny>> {
            let cls = py.import("cosmolkit")?.getattr("BatchValidationError")?;
            cls.call1((message, error_count, reason, record_errors))
        })() {
            Ok(instance) => PyErr::from_value(instance),
            Err(error) => error,
        }
    })
}

fn fingerprint_pyerr(error: cosmolkit_core::FingerprintError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn svg_draw_pyerr(error: cosmolkit_core::SvgDrawError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn pdb_molecule_pyerr(error: cosmolkit_core::PdbMoleculeConversionError) -> PyErr {
    match error {
        cosmolkit_core::PdbMoleculeConversionError::Unsupported(message) => {
            PyNotImplementedError::new_err(message)
        }
        other => PyValueError::new_err(other.to_string()),
    }
}

fn fragment_pyerr(error: cosmolkit_core::fragment::FragmentError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn hash_pyerr(error: cosmolkit_core::mol_hash::HashError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn pickle_pyerr(error: cosmolkit_core::PickleError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn stereo_pyerr(error: cosmolkit_core::StereoError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn smarts_parse_pyerr(error: cosmolkit_core::SmartsParseError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn forcefield_pyerr(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn distgeom_pyerr(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

#[allow(clippy::too_many_arguments)]
fn make_morgan_fingerprint_params(
    radius: u32,
    n_bits: usize,
    include_chirality: bool,
    use_bond_types: bool,
    count_simulation: bool,
    count_bounds: Option<Vec<u32>>,
    only_nonzero_invariants: bool,
    include_redundant_environments: bool,
    from_atoms: Option<Vec<usize>>,
    ignore_atoms: Option<Vec<usize>>,
    custom_atom_invariants: Option<Vec<u32>>,
    custom_bond_invariants: Option<Vec<u32>>,
    atom_invariants_generator: Option<&str>,
    atom_invariants_include_ring_membership: bool,
    bond_invariants_generator: Option<&str>,
    bond_invariants_use_bond_types: bool,
    bond_invariants_use_chirality: bool,
    num_bits_per_feature: u32,
    collect_additional_output: bool,
) -> PyResult<cosmolkit_core::MorganFingerprintParams> {
    let atom_invariants_generator = match atom_invariants_generator
        .map(|value| value.to_ascii_lowercase())
        .as_deref()
    {
        None | Some("connectivity") | Some("morgan") => {
            cosmolkit_core::MorganAtomInvariantsGenerator::Connectivity {
                include_ring_membership: atom_invariants_include_ring_membership,
            }
        }
        Some("feature") | Some("fcfp") => cosmolkit_core::MorganAtomInvariantsGenerator::Feature,
        Some(value) => {
            return Err(PyValueError::new_err(format!(
                "unsupported atom_invariants_generator '{value}', expected one of: connectivity, morgan, feature, fcfp"
            )));
        }
    };

    let bond_invariants_generator = match bond_invariants_generator
        .map(|value| value.to_ascii_lowercase())
        .as_deref()
    {
        None => None,
        Some("morgan") | Some("default") | Some("bond") => {
            Some(cosmolkit_core::MorganBondInvariantsGenerator {
                use_bond_types: bond_invariants_use_bond_types,
                use_chirality: bond_invariants_use_chirality,
            })
        }
        Some(value) => {
            return Err(PyValueError::new_err(format!(
                "unsupported bond_invariants_generator '{value}', expected one of: morgan, default, bond"
            )));
        }
    };

    Ok(cosmolkit_core::MorganFingerprintParams {
        radius,
        n_bits,
        use_chirality: include_chirality,
        use_bond_types,
        count_simulation,
        count_bounds: count_bounds.unwrap_or_else(|| vec![1, 2, 4, 8]),
        only_nonzero_invariants,
        include_ring_membership: atom_invariants_include_ring_membership,
        include_redundant_environments,
        from_atoms,
        ignore_atoms,
        custom_atom_invariants,
        custom_bond_invariants,
        atom_invariants_generator,
        bond_invariants_generator,
        num_bits_per_feature,
        collect_additional_output,
    })
}

fn make_avalon_fingerprint_params(
    min_path: u32,
    max_path: u32,
    n_bits: usize,
    n_bits_per_hash: u32,
    use_bond_order: bool,
    use_hs: bool,
    tautomeric_fingerprint: bool,
    from_atoms: Option<Vec<usize>>,
) -> cosmolkit_core::avalon_fingerprint::AvalonFingerprintParams {
    cosmolkit_core::avalon_fingerprint::AvalonFingerprintParams {
        min_path,
        max_path,
        n_bits,
        n_bits_per_hash,
        use_bond_order,
        use_hs,
        tautomeric_fingerprint,
        from_atoms,
    }
}

fn make_topological_fingerprint_params(
    min_path: u32,
    max_path: u32,
    n_bits: usize,
    n_bits_per_hash: u32,
    use_bond_types: bool,
    from_atoms: Option<Vec<usize>>,
    ignore_atoms: Option<Vec<usize>>,
) -> cosmolkit_core::fingerprint::TopologicalFingerprintParams {
    cosmolkit_core::fingerprint::TopologicalFingerprintParams {
        min_path,
        max_path,
        n_bits,
        n_bits_per_hash,
        use_bond_types,
        from_atoms,
        ignore_atoms,
    }
}

fn reject_non_strict_sanitize(strict: Option<bool>) -> PyResult<()> {
    if matches!(strict, Some(false)) {
        return Err(PyValueError::new_err(
            "strict=False sanitization is not implemented; COSMolKit currently supports RDKit-style strict sanitization only",
        ));
    }
    Ok(())
}

fn reject_unsanitized_mol_reader(sanitize: Option<bool>) -> PyResult<()> {
    if matches!(sanitize, Some(false)) {
        return Err(PyValueError::new_err(
            "sanitize=False is not implemented for SDF/molfile readers; MolBlock parsing currently finalizes chemistry during read",
        ));
    }
    Ok(())
}

fn parse_mol2_variant(variant: &str) -> PyResult<cosmolkit_core::Mol2Type> {
    match variant.to_ascii_lowercase().as_str() {
        "corina" => Ok(cosmolkit_core::Mol2Type::Corina),
        _ => Err(PyValueError::new_err(format!(
            "unsupported MOL2 variant '{variant}', expected 'corina'"
        ))),
    }
}

fn write_batch_report(path: &str, report: &cosmolkit_core::BatchExportReport) -> PyResult<()> {
    let expanded_path = expand_user_path(path)?;
    let ext = expanded_path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("json")
        .to_ascii_lowercase();
    let content = if ext == "csv" {
        format!(
            "written,skipped,failed\n{},{},{}\n",
            report.written, report.skipped, report.failed
        )
    } else {
        format!(
            "{{\n  \"written\": {},\n  \"skipped\": {},\n  \"failed\": {}\n}}\n",
            report.written, report.skipped, report.failed
        )
    };
    fs::write(&expanded_path, content)
        .map_err(|err| PyValueError::new_err(format!("write error report failed: {err}")))
}

fn complete_batch_filenames(
    filenames: Option<Vec<Option<String>>>,
    total: usize,
    extension: &str,
) -> PyResult<Option<Vec<String>>> {
    let Some(filenames) = filenames else {
        return Ok(None);
    };
    if filenames.len() != total {
        return Err(PyValueError::new_err(format!(
            "filenames length must match batch length: expected {total}, got {}",
            filenames.len()
        )));
    }
    Ok(Some(
        filenames
            .into_iter()
            .enumerate()
            .map(|(index, filename)| filename.unwrap_or_else(|| format!("mol_{index}.{extension}")))
            .collect(),
    ))
}

fn to_python_tetrahedral_stereo(
    mol: &cosmolkit_core::Molecule,
) -> PyResult<Vec<(usize, Vec<Option<usize>>)>> {
    Ok(mol
        .tetrahedral_stereo()
        .map_err(|err| PyValueError::new_err(format!("tetrahedral_stereo failed: {err}")))?
        .into_iter()
        .map(|stereo| {
            let ligands = stereo
                .ligands
                .into_iter()
                .map(|ligand| match ligand {
                    cosmolkit_core::LigandRef::Atom(index) => Some(index.index()),
                    cosmolkit_core::LigandRef::ImplicitHydrogen => None,
                })
                .collect();
            (stereo.center.index(), ligands)
        })
        .collect())
}

fn enum_member<'py>(py: Python<'py>, enum_name: &str, code: i64) -> PyResult<Bound<'py, PyAny>> {
    let module = py.import("cosmolkit")?;
    module.getattr(enum_name)?.call1((code,))
}

fn residue_info_kind_code(kind: cosmolkit_core::ResidueInfoKind) -> i64 {
    kind as i64
}

fn residue_info_kind_name(kind: cosmolkit_core::ResidueInfoKind) -> &'static str {
    match kind {
        cosmolkit_core::ResidueInfoKind::Unknown => "UNKNOWN",
        cosmolkit_core::ResidueInfoKind::Aa => "AA",
        cosmolkit_core::ResidueInfoKind::Aad => "AAD",
        cosmolkit_core::ResidueInfoKind::Paa => "PAA",
        cosmolkit_core::ResidueInfoKind::Maa => "MAA",
        cosmolkit_core::ResidueInfoKind::Rna => "RNA",
        cosmolkit_core::ResidueInfoKind::Dna => "DNA",
        cosmolkit_core::ResidueInfoKind::Buf => "BUF",
        cosmolkit_core::ResidueInfoKind::Hoh => "HOH",
        cosmolkit_core::ResidueInfoKind::Pyr => "PYR",
        cosmolkit_core::ResidueInfoKind::Ket => "KET",
        cosmolkit_core::ResidueInfoKind::Els => "ELS",
    }
}

fn residue_info_kind_from_code(code: i64) -> PyResult<cosmolkit_core::ResidueInfoKind> {
    match code {
        0 => Ok(cosmolkit_core::ResidueInfoKind::Unknown),
        1 => Ok(cosmolkit_core::ResidueInfoKind::Aa),
        2 => Ok(cosmolkit_core::ResidueInfoKind::Aad),
        3 => Ok(cosmolkit_core::ResidueInfoKind::Paa),
        4 => Ok(cosmolkit_core::ResidueInfoKind::Maa),
        5 => Ok(cosmolkit_core::ResidueInfoKind::Rna),
        6 => Ok(cosmolkit_core::ResidueInfoKind::Dna),
        7 => Ok(cosmolkit_core::ResidueInfoKind::Buf),
        8 => Ok(cosmolkit_core::ResidueInfoKind::Hoh),
        9 => Ok(cosmolkit_core::ResidueInfoKind::Pyr),
        10 => Ok(cosmolkit_core::ResidueInfoKind::Ket),
        11 => Ok(cosmolkit_core::ResidueInfoKind::Els),
        _ => Err(PyValueError::new_err(format!(
            "unsupported ResidueInfoKind code {code}"
        ))),
    }
}

fn residue_info_sequence_pyerr(error: cosmolkit_core::ResidueInfoSequenceError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn residue_code_enum_member<'py>(
    py: Python<'py>,
    code: cosmolkit_core::ResidueCode,
) -> PyResult<Bound<'py, PyAny>> {
    enum_member(py, "ResidueCode", i64::from(code.as_u16()))
}

fn residue_info_kind_enum_member<'py>(
    py: Python<'py>,
    kind: cosmolkit_core::ResidueInfoKind,
) -> PyResult<Bound<'py, PyAny>> {
    enum_member(py, "ResidueInfoKind", residue_info_kind_code(kind))
}

fn add_int_enum(
    m: &Bound<'_, PyModule>,
    enum_name: &str,
    map_name: &str,
    members: &[(&str, i64)],
) -> PyResult<()> {
    let members_with_map_keys = members
        .iter()
        .map(|(name, value)| (*name, *value, *name))
        .collect::<Vec<_>>();
    add_int_enum_with_map_keys(m, enum_name, map_name, &members_with_map_keys)
}

fn add_int_enum_with_map_keys(
    m: &Bound<'_, PyModule>,
    enum_name: &str,
    map_name: &str,
    members: &[(&str, i64, &str)],
) -> PyResult<()> {
    let py = m.py();
    let enum_module = py.import("enum")?;
    let int_enum = enum_module.getattr("IntEnum")?;
    let member_dict = PyDict::new(py);
    for (name, value, _) in members {
        member_dict.set_item(name, value)?;
    }
    let enum_cls = int_enum.call1((enum_name, member_dict))?;
    m.add(enum_name, &enum_cls)?;

    let enum_map = PyDict::new(py);
    for (name, _, map_key) in members {
        enum_map.set_item(map_key, enum_cls.getattr(name)?)?;
    }
    let proxy = PyMappingProxy::new(py, enum_map.cast::<PyMapping>()?);
    m.add(map_name, proxy)?;
    Ok(())
}

fn add_int_enum_with_map_aliases(
    m: &Bound<'_, PyModule>,
    enum_name: &str,
    map_name: &str,
    members: &[(&str, i64)],
    map_aliases: &[(&str, &str)],
) -> PyResult<()> {
    let py = m.py();
    let enum_module = py.import("enum")?;
    let int_enum = enum_module.getattr("IntEnum")?;
    let member_dict = PyDict::new(py);
    for (name, value) in members {
        member_dict.set_item(name, value)?;
    }
    let enum_cls = int_enum.call1((enum_name, member_dict))?;
    m.add(enum_name, &enum_cls)?;

    let enum_map = PyDict::new(py);
    for (map_key, member_name) in map_aliases {
        enum_map.set_item(map_key, enum_cls.getattr(member_name)?)?;
    }
    let proxy = PyMappingProxy::new(py, enum_map.cast::<PyMapping>()?);
    m.add(map_name, proxy)?;
    Ok(())
}

fn add_public_enums(m: &Bound<'_, PyModule>) -> PyResult<()> {
    add_int_enum(
        m,
        "BondOrder",
        "BOND_ORDER_MAP",
        &[
            ("UNSPECIFIED", 0),
            ("SINGLE", 1),
            ("DOUBLE", 2),
            ("TRIPLE", 3),
            ("QUADRUPLE", 4),
            ("QUINTUPLE", 5),
            ("HEXTUPLE", 6),
            ("ONEANDAHALF", 7),
            ("TWOANDAHALF", 8),
            ("THREEANDAHALF", 9),
            ("FOURANDAHALF", 10),
            ("FIVEANDAHALF", 11),
            ("AROMATIC", 12),
            ("IONIC", 13),
            ("DATIVE", 14),
            ("DATIVEONE", 15),
            ("DATIVEL", 16),
            ("DATIVER", 17),
            ("HYDROGEN", 18),
            ("THREECENTER", 19),
            ("OTHER", 20),
            ("ZERO", 21),
        ],
    )?;
    add_int_enum(
        m,
        "BondDirection",
        "BOND_DIRECTION_MAP",
        &[
            ("NONE", 0),
            ("BEGINWEDGE", 1),
            ("BEGINDASH", 2),
            ("ENDUPRIGHT", 3),
            ("ENDDOWNRIGHT", 4),
            ("EITHERDOUBLE", 5),
            ("UNKNOWN", 6),
        ],
    )?;
    add_int_enum_with_map_aliases(
        m,
        "BondStereo",
        "BOND_STEREO_MAP",
        &[
            ("NONE", 0),
            ("ANY", 1),
            ("Z", 2),
            ("E", 3),
            ("CIS", 4),
            ("TRANS", 5),
            ("ATROP_CW", 6),
            ("ATROP_CCW", 7),
        ],
        &[
            ("NONE", "NONE"),
            ("STEREONONE", "NONE"),
            ("ANY", "ANY"),
            ("STEREOANY", "ANY"),
            ("Z", "Z"),
            ("STEREOZ", "Z"),
            ("E", "E"),
            ("STEREOE", "E"),
            ("CIS", "CIS"),
            ("TRANS", "TRANS"),
            ("ATROP_CW", "ATROP_CW"),
            ("STEREOATROPCW", "ATROP_CW"),
            ("ATROP_CCW", "ATROP_CCW"),
            ("STEREOATROPCCW", "ATROP_CCW"),
        ],
    )?;
    add_int_enum(
        m,
        "ChiralTag",
        "CHIRAL_TAG_MAP",
        &[
            ("CHI_UNSPECIFIED", 0),
            ("CHI_TETRAHEDRAL_CW", 1),
            ("CHI_TETRAHEDRAL_CCW", 2),
            ("CHI_OTHER", 3),
            ("CHI_TETRAHEDRAL", 4),
            ("CHI_ALLENE", 5),
            ("CHI_SQUAREPLANAR", 6),
            ("CHI_TRIGONALBIPYRAMIDAL", 7),
            ("CHI_OCTAHEDRAL", 8),
        ],
    )?;
    add_int_enum_with_map_keys(
        m,
        "BatchErrorMode",
        "BATCH_ERROR_MODE_MAP",
        &[("RAISE", 1, "raise"), ("KEEP", 2, "keep")],
    )?;
    let residue_code_members = cosmolkit_core::RESIDUE_INFO_TABLE
        .iter()
        .map(|info| (format!("{:?}", info.code), i64::from(info.code.as_u16())))
        .collect::<Vec<_>>();
    let residue_code_member_refs = residue_code_members
        .iter()
        .map(|(name, value)| (name.as_str(), *value))
        .collect::<Vec<_>>();
    let mut residue_code_aliases = cosmolkit_core::RESIDUE_INFO_TABLE
        .iter()
        .map(|info| (info.name.to_string(), format!("{:?}", info.code)))
        .collect::<Vec<_>>();
    residue_code_aliases.extend([
        ("TRY".to_string(), "TRP".to_string()),
        ("WAT".to_string(), "HOH".to_string()),
        ("H2O".to_string(), "HOH".to_string()),
        ("+A".to_string(), "DA".to_string()),
        ("+C".to_string(), "DC".to_string()),
        ("+G".to_string(), "DG".to_string()),
        ("+I".to_string(), "DI".to_string()),
        ("+T".to_string(), "DT".to_string()),
        ("+U".to_string(), "DU".to_string()),
        ("+N".to_string(), "DN".to_string()),
    ]);
    let residue_code_alias_refs = residue_code_aliases
        .iter()
        .map(|(key, name)| (key.as_str(), name.as_str()))
        .collect::<Vec<_>>();
    add_int_enum_with_map_aliases(
        m,
        "ResidueCode",
        "RESIDUE_CODE_MAP",
        &residue_code_member_refs,
        &residue_code_alias_refs,
    )?;
    add_int_enum_with_map_aliases(
        m,
        "ResidueInfoKind",
        "RESIDUE_INFO_KIND_MAP",
        &[
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
        ],
        &[
            ("UNKNOWN", "UNKNOWN"),
            ("AA", "AA"),
            ("AAD", "AAD"),
            ("PAA", "PAA"),
            ("MAA", "MAA"),
            ("RNA", "RNA"),
            ("DNA", "DNA"),
            ("BUF", "BUF"),
            ("HOH", "HOH"),
            ("PYR", "PYR"),
            ("KET", "KET"),
            ("ELS", "ELS"),
        ],
    )
}

fn add_batch_validation_error_class(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let globals = PyDict::new(py);
    globals.set_item("ValueError", py.get_type::<PyValueError>())?;
    let code = r#"
class BatchValidationError(ValueError):
    __module__ = "cosmolkit"

    def __init__(self, message, error_count=0, reason=None, record_errors=None):
        super().__init__(message)
        self.error_count = int(error_count)
        self.reason = reason
        self._errors = list(record_errors or [])

    def errors(self):
        return list(self._errors)
"#;
    py.import("builtins")?
        .getattr("exec")?
        .call1((code, &globals))?;
    let cls = globals
        .get_item("BatchValidationError")?
        .ok_or_else(|| PyValueError::new_err("failed to create BatchValidationError class"))?;
    m.add("BatchValidationError", cls)?;
    Ok(())
}

fn parse_sdf_format(format: Option<&str>) -> PyResult<SdfFormat> {
    match format.map(|s| s.to_ascii_lowercase()) {
        None => Ok(SdfFormat::V2000),
        Some(v) if v == "v2000" || v == "v2k" => Ok(SdfFormat::V2000),
        Some(v) if v == "v3000" || v == "v3k" => Ok(SdfFormat::V3000),
        Some(v) => Err(PyValueError::new_err(format!(
            "unsupported SDF format '{v}', expected one of: v2000, v3000"
        ))),
    }
}

fn molecule_to_sdf_record_string(
    mol: &cosmolkit_core::Molecule,
    format: SdfFormat,
    include_stereo: bool,
    kekulize: bool,
) -> Result<String, cosmolkit_core::io::molblock::MolWriteError> {
    let params = cosmolkit_core::io::molblock::MolBlockWriteParams {
        format,
        force_2d: false,
        include_stereo,
        kekulize,
        ..Default::default()
    };
    molblock::mol_to_sdf_record_with_params(mol, &params)
}

fn molecule_to_2d_sdf_record_string(
    mol: &cosmolkit_core::Molecule,
    format: SdfFormat,
    include_stereo: bool,
    kekulize: bool,
) -> Result<String, cosmolkit_core::io::molblock::MolWriteError> {
    let params = cosmolkit_core::io::molblock::MolBlockWriteParams {
        format,
        force_2d: true,
        include_stereo,
        kekulize,
        ..Default::default()
    };
    let with_coords = if mol.coordinates_2d().is_some() {
        mol.clone()
    } else {
        mol.with_2d_coordinates()?
    };
    molblock::mol_to_sdf_record_with_params(&with_coords, &params)
}

fn molecule_to_3d_sdf_record_string(
    mol: &cosmolkit_core::Molecule,
    format: SdfFormat,
    include_stereo: bool,
    kekulize: bool,
) -> Result<String, cosmolkit_core::io::molblock::MolWriteError> {
    let params = cosmolkit_core::io::molblock::MolBlockWriteParams {
        format,
        force_2d: false,
        include_stereo,
        kekulize,
        ..Default::default()
    };
    if !mol.conformers_3d().is_empty() {
        molblock::mol_to_sdf_record_with_params(mol, &params)
    } else {
        Err(
            cosmolkit_core::io::molblock::MolWriteError::UnsupportedSubset(
                "3D coordinates are required; read a 3D SDF record or add a 3D conformer before writing 3D SDF",
            ),
        )
    }
}

fn expand_user_path(path: &str) -> PyResult<PathBuf> {
    if path == "~" || path.starts_with("~/") {
        let home = std::env::var_os("HOME")
            .ok_or_else(|| PyValueError::new_err("cannot expand '~': HOME is not set"))?;
        let mut expanded = PathBuf::from(home);
        if let Some(rest) = path.strip_prefix("~/") {
            expanded.push(rest);
        }
        Ok(expanded)
    } else {
        Ok(PathBuf::from(path))
    }
}

fn atomic_number_from_element(element: &str) -> Option<u8> {
    match element {
        "H" => Some(1),
        "B" => Some(5),
        "C" => Some(6),
        "N" => Some(7),
        "O" => Some(8),
        "F" => Some(9),
        "P" => Some(15),
        "S" => Some(16),
        "Cl" => Some(17),
        "Br" => Some(35),
        "I" => Some(53),
        "Na" => Some(11),
        "K" => Some(19),
        "Li" => Some(3),
        "Mg" => Some(12),
        "Ca" => Some(20),
        "Fe" => Some(26),
        "Cu" => Some(29),
        "Zn" => Some(30),
        "Si" => Some(14),
        "Al" => Some(13),
        "*" => Some(0),
        _ => None,
    }
}

fn py_method<'py>(obj: &Bound<'py, PyAny>, method: &str) -> PyResult<Bound<'py, PyAny>> {
    obj.call_method0(method)
        .map_err(|err| PyValueError::new_err(format!("from_rdkit failed calling {method}: {err}")))
}

fn py_method_index<'py>(
    obj: &Bound<'py, PyAny>,
    method: &str,
    index: usize,
) -> PyResult<Bound<'py, PyAny>> {
    obj.call_method1(method, (index,))
        .map_err(|err| PyValueError::new_err(format!("from_rdkit failed calling {method}: {err}")))
}

fn py_method_extract<T>(obj: &Bound<'_, PyAny>, method: &str) -> PyResult<T>
where
    for<'a> T: FromPyObject<'a, 'a>,
{
    py_method(obj, method)?.extract::<T>().map_err(|_| {
        PyValueError::new_err(format!("from_rdkit failed extracting result from {method}"))
    })
}

fn py_attr_f64(obj: &Bound<'_, PyAny>, attr: &str) -> PyResult<f64> {
    let value = obj.getattr(attr).map_err(|err| {
        PyValueError::new_err(format!(
            "from_rdkit failed accessing attribute {attr}: {err}"
        ))
    })?;
    value.extract::<f64>().map_err(|err| {
        PyValueError::new_err(format!(
            "from_rdkit failed extracting attribute {attr}: {err}"
        ))
    })
}

fn py_method_str(obj: &Bound<'_, PyAny>, method: &str) -> PyResult<String> {
    let value = py_method(obj, method)?;
    Ok(value
        .str()
        .map_err(|err| {
            PyValueError::new_err(format!("from_rdkit failed stringifying {method}: {err}"))
        })?
        .to_string_lossy()
        .into_owned())
}

fn rdkit_chiral_tag_from_name(name: &str) -> PyResult<cosmolkit_core::ChiralTag> {
    match name {
        "CHI_UNSPECIFIED" => Ok(cosmolkit_core::ChiralTag::Unspecified),
        "CHI_TETRAHEDRAL_CW" => Ok(cosmolkit_core::ChiralTag::TetrahedralCw),
        "CHI_TETRAHEDRAL_CCW" => Ok(cosmolkit_core::ChiralTag::TetrahedralCcw),
        "CHI_OTHER" => Ok(cosmolkit_core::ChiralTag::Other),
        "CHI_TETRAHEDRAL" => Ok(cosmolkit_core::ChiralTag::Tetrahedral),
        "CHI_ALLENE" => Ok(cosmolkit_core::ChiralTag::Allene),
        "CHI_SQUAREPLANAR" => Ok(cosmolkit_core::ChiralTag::SquarePlanar),
        "CHI_TRIGONALBIPYRAMIDAL" => Ok(cosmolkit_core::ChiralTag::TrigonalBipyramidal),
        "CHI_OCTAHEDRAL" => Ok(cosmolkit_core::ChiralTag::Octahedral),
        _ => Err(PyValueError::new_err(format!(
            "from_rdkit unsupported atom chiral tag '{name}'"
        ))),
    }
}

fn rdkit_bond_order_from_name(name: &str) -> PyResult<cosmolkit_core::BondOrder> {
    match name {
        "UNSPECIFIED" | "ZERO" => Ok(cosmolkit_core::BondOrder::Unspecified),
        "SINGLE" => Ok(cosmolkit_core::BondOrder::Single),
        "DOUBLE" => Ok(cosmolkit_core::BondOrder::Double),
        "TRIPLE" => Ok(cosmolkit_core::BondOrder::Triple),
        "QUADRUPLE" => Ok(cosmolkit_core::BondOrder::Quadruple),
        "QUINTUPLE" => Ok(cosmolkit_core::BondOrder::Quintuple),
        "HEXTUPLE" => Ok(cosmolkit_core::BondOrder::Hextuple),
        "ONEANDAHALF" => Ok(cosmolkit_core::BondOrder::OneAndHalf),
        "TWOANDAHALF" => Ok(cosmolkit_core::BondOrder::TwoAndHalf),
        "THREEANDAHALF" => Ok(cosmolkit_core::BondOrder::ThreeAndHalf),
        "FOURANDAHALF" => Ok(cosmolkit_core::BondOrder::FourAndHalf),
        "FIVEANDAHALF" => Ok(cosmolkit_core::BondOrder::FiveAndHalf),
        "AROMATIC" => Ok(cosmolkit_core::BondOrder::Aromatic),
        "IONIC" => Ok(cosmolkit_core::BondOrder::Ionic),
        "DATIVE" => Ok(cosmolkit_core::BondOrder::Dative),
        "DATIVEONE" => Ok(cosmolkit_core::BondOrder::DativeOne),
        "DATIVEL" => Ok(cosmolkit_core::BondOrder::DativeLeft),
        "DATIVER" => Ok(cosmolkit_core::BondOrder::DativeRight),
        "HYDROGEN" => Ok(cosmolkit_core::BondOrder::Hydrogen),
        "THREECENTER" => Ok(cosmolkit_core::BondOrder::ThreeCenter),
        "OTHER" => Ok(cosmolkit_core::BondOrder::Other),
        _ => Err(PyValueError::new_err(format!(
            "from_rdkit unsupported bond type '{name}'"
        ))),
    }
}

fn rdkit_bond_direction_from_name(name: &str) -> PyResult<cosmolkit_core::BondDirection> {
    match name {
        "NONE" => Ok(cosmolkit_core::BondDirection::None),
        "BEGINWEDGE" => Ok(cosmolkit_core::BondDirection::BeginWedge),
        "BEGINDASH" => Ok(cosmolkit_core::BondDirection::BeginDash),
        "ENDUPRIGHT" => Ok(cosmolkit_core::BondDirection::EndUpRight),
        "ENDDOWNRIGHT" => Ok(cosmolkit_core::BondDirection::EndDownRight),
        "EITHERDOUBLE" => Ok(cosmolkit_core::BondDirection::EitherDouble),
        "UNKNOWN" => Ok(cosmolkit_core::BondDirection::Unknown),
        _ => Err(PyValueError::new_err(format!(
            "from_rdkit unsupported bond direction '{name}'"
        ))),
    }
}

fn rdkit_bond_stereo_from_name(name: &str) -> PyResult<cosmolkit_core::BondStereo> {
    match name {
        "NONE" | "STEREONONE" => Ok(cosmolkit_core::BondStereo::None),
        "ANY" | "STEREOANY" => Ok(cosmolkit_core::BondStereo::Any),
        "Z" | "STEREOZ" => Ok(cosmolkit_core::BondStereo::Z),
        "E" | "STEREOE" => Ok(cosmolkit_core::BondStereo::E),
        "CIS" => Ok(cosmolkit_core::BondStereo::Cis),
        "TRANS" => Ok(cosmolkit_core::BondStereo::Trans),
        "ATROP_CW" => Ok(cosmolkit_core::BondStereo::AtropCw),
        "ATROP_CCW" => Ok(cosmolkit_core::BondStereo::AtropCcw),
        _ => Err(PyValueError::new_err(format!(
            "from_rdkit unsupported bond stereo '{name}'"
        ))),
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
#[doc = r#"
A molecule value.

``Molecule`` stores atoms, bonds, stereochemistry, and optional coordinate
data. Transformation methods such as ``with_hydrogens()``,
``without_hydrogens()``, ``with_kekulized_bonds()``, and
``with_2d_coordinates()`` return new molecule values. The original molecule is
left unchanged.

Internally COSMolKit uses copy-on-write storage to share unchanged molecular
data efficiently, but the public Python contract is value semantics.

In-place methods mutate the receiver and always end with ``_``. COSMolKit
reserves the trailing underscore for this single public ``Molecule`` meaning.

Examples
--------
Create molecules with ``Molecule.from_smiles()``, transform them with value
methods such as ``with_2d_coordinates()``, then export strings, arrays, or
depiction files.
"#]
struct Molecule {
    inner: cosmolkit_core::Molecule,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "EmbedParameters", skip_from_py_object)]
#[derive(Clone)]
struct PyEmbedParameters {
    inner: cosmolkit_core::EmbedParameters,
}

impl PyEmbedParameters {
    fn from_inner(inner: cosmolkit_core::EmbedParameters) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyEmbedParameters {
    #[new]
    fn new() -> Self {
        Self::from_inner(cosmolkit_core::EmbedParameters::default())
    }

    #[staticmethod]
    fn dg() -> Self {
        Self::from_inner(cosmolkit_core::EmbedParameters::default())
    }

    #[staticmethod]
    fn kdg() -> Self {
        Self::from_inner(cosmolkit_core::EmbedParameters::kdg())
    }

    #[staticmethod]
    fn etdg() -> Self {
        Self::from_inner(cosmolkit_core::EmbedParameters::etdg())
    }

    #[staticmethod]
    fn etdg_v2() -> Self {
        Self::from_inner(cosmolkit_core::EmbedParameters::etdg_v2())
    }

    #[staticmethod]
    fn etkdg() -> Self {
        Self::from_inner(cosmolkit_core::EmbedParameters::etkdg())
    }

    #[staticmethod]
    fn etkdg_v2() -> Self {
        Self::from_inner(cosmolkit_core::EmbedParameters::etkdg_v2())
    }

    #[staticmethod]
    fn etkdg_v3() -> Self {
        Self::from_inner(cosmolkit_core::EmbedParameters::etkdg_v3())
    }

    #[staticmethod]
    fn sr_etkdg_v3() -> Self {
        Self::from_inner(cosmolkit_core::EmbedParameters::sr_etkdg_v3())
    }

    #[getter]
    fn random_seed(&self) -> i32 {
        self.inner.random_seed
    }

    #[setter]
    fn set_random_seed(&mut self, value: i32) {
        self.inner.random_seed = value;
    }

    #[getter]
    fn max_iterations(&self) -> u32 {
        self.inner.max_iterations
    }

    #[setter]
    fn set_max_iterations(&mut self, value: u32) {
        self.inner.max_iterations = value;
    }

    #[getter]
    fn num_threads(&self) -> i32 {
        self.inner.num_threads
    }

    #[setter]
    fn set_num_threads(&mut self, value: i32) {
        self.inner.num_threads = value;
    }

    #[getter]
    fn box_size_mult(&self) -> f64 {
        self.inner.box_size_mult
    }

    #[setter]
    fn set_box_size_mult(&mut self, value: f64) {
        self.inner.box_size_mult = value;
    }

    #[getter]
    fn rand_neg_eig(&self) -> bool {
        self.inner.rand_neg_eig
    }

    #[setter]
    fn set_rand_neg_eig(&mut self, value: bool) {
        self.inner.rand_neg_eig = value;
    }

    #[getter]
    fn num_zero_fail(&self) -> u32 {
        self.inner.num_zero_fail
    }

    #[setter]
    fn set_num_zero_fail(&mut self, value: u32) {
        self.inner.num_zero_fail = value;
    }

    #[getter]
    fn coord_map(&self) -> Option<BTreeMap<i32, (f64, f64, f64)>> {
        self.inner.coord_map.as_ref().map(|coord_map| {
            coord_map
                .iter()
                .map(|(&idx, point)| (idx, (point.x, point.y, point.z)))
                .collect()
        })
    }

    #[setter]
    fn set_coord_map(&mut self, value: Option<BTreeMap<i32, (f64, f64, f64)>>) {
        self.inner.coord_map = value.map(|coord_map| {
            coord_map
                .into_iter()
                .map(|(idx, (x, y, z))| (idx, cosmolkit_core::ForceFieldVec3::new(x, y, z)))
                .collect()
        });
    }

    #[getter]
    fn optimizer_force_tol(&self) -> f64 {
        self.inner.optimizer_force_tol
    }

    #[setter]
    fn set_optimizer_force_tol(&mut self, value: f64) {
        self.inner.optimizer_force_tol = value;
    }

    #[getter]
    fn ignore_smoothing_failures(&self) -> bool {
        self.inner.ignore_smoothing_failures
    }

    #[setter]
    fn set_ignore_smoothing_failures(&mut self, value: bool) {
        self.inner.ignore_smoothing_failures = value;
    }

    #[getter]
    fn prune_rms_thresh(&self) -> f64 {
        self.inner.prune_rms_thresh
    }

    #[setter]
    fn set_prune_rms_thresh(&mut self, value: f64) {
        self.inner.prune_rms_thresh = value;
    }

    #[getter]
    fn clear_confs(&self) -> bool {
        self.inner.clear_confs
    }

    #[setter]
    fn set_clear_confs(&mut self, value: bool) {
        self.inner.clear_confs = value;
    }

    #[getter]
    fn use_random_coords(&self) -> bool {
        self.inner.use_random_coords
    }

    #[setter]
    fn set_use_random_coords(&mut self, value: bool) {
        self.inner.use_random_coords = value;
    }

    #[getter]
    fn enforce_chirality(&self) -> bool {
        self.inner.enforce_chirality
    }

    #[setter]
    fn set_enforce_chirality(&mut self, value: bool) {
        self.inner.enforce_chirality = value;
    }

    #[getter]
    fn use_exp_torsion_angle_prefs(&self) -> bool {
        self.inner.use_exp_torsion_angle_prefs
    }

    #[setter]
    fn set_use_exp_torsion_angle_prefs(&mut self, value: bool) {
        self.inner.use_exp_torsion_angle_prefs = value;
    }

    #[getter]
    fn use_basic_knowledge(&self) -> bool {
        self.inner.use_basic_knowledge
    }

    #[setter]
    fn set_use_basic_knowledge(&mut self, value: bool) {
        self.inner.use_basic_knowledge = value;
    }

    #[getter]
    fn verbose(&self) -> bool {
        self.inner.verbose
    }

    #[setter]
    fn set_verbose(&mut self, value: bool) {
        self.inner.verbose = value;
    }

    #[getter]
    fn basin_thresh(&self) -> f64 {
        self.inner.basin_thresh
    }

    #[setter]
    fn set_basin_thresh(&mut self, value: f64) {
        self.inner.basin_thresh = value;
    }

    #[getter]
    fn only_heavy_atoms_for_rms(&self) -> bool {
        self.inner.only_heavy_atoms_for_rms
    }

    #[setter]
    fn set_only_heavy_atoms_for_rms(&mut self, value: bool) {
        self.inner.only_heavy_atoms_for_rms = value;
    }

    #[getter]
    fn et_version(&self) -> u32 {
        self.inner.et_version
    }

    #[setter]
    fn set_et_version(&mut self, value: u32) {
        self.inner.et_version = value;
    }

    #[getter]
    fn embed_fragments_separately(&self) -> bool {
        self.inner.embed_fragments_separately
    }

    #[setter]
    fn set_embed_fragments_separately(&mut self, value: bool) {
        self.inner.embed_fragments_separately = value;
    }

    #[getter]
    fn use_small_ring_torsions(&self) -> bool {
        self.inner.use_small_ring_torsions
    }

    #[setter]
    fn set_use_small_ring_torsions(&mut self, value: bool) {
        self.inner.use_small_ring_torsions = value;
    }

    #[getter]
    fn use_macrocycle_torsions(&self) -> bool {
        self.inner.use_macrocycle_torsions
    }

    #[setter]
    fn set_use_macrocycle_torsions(&mut self, value: bool) {
        self.inner.use_macrocycle_torsions = value;
    }

    #[getter]
    fn use_macrocycle14config(&self) -> bool {
        self.inner.use_macrocycle14config
    }

    #[setter]
    fn set_use_macrocycle14config(&mut self, value: bool) {
        self.inner.use_macrocycle14config = value;
    }

    #[getter]
    fn timeout(&self) -> u32 {
        self.inner.timeout
    }

    #[setter]
    fn set_timeout(&mut self, value: u32) {
        self.inner.timeout = value;
    }

    #[getter]
    fn cpci(&self) -> Option<BTreeMap<(u32, u32), f64>> {
        self.inner.cpci.clone()
    }

    #[setter]
    fn set_cpci(&mut self, value: Option<BTreeMap<(u32, u32), f64>>) {
        self.inner.cpci = value;
    }

    #[getter]
    fn force_trans_amides(&self) -> bool {
        self.inner.force_trans_amides
    }

    #[setter]
    fn set_force_trans_amides(&mut self, value: bool) {
        self.inner.force_trans_amides = value;
    }

    #[getter]
    fn use_symmetry_for_pruning(&self) -> bool {
        self.inner.use_symmetry_for_pruning
    }

    #[setter]
    fn set_use_symmetry_for_pruning(&mut self, value: bool) {
        self.inner.use_symmetry_for_pruning = value;
    }

    #[getter]
    fn bounds_mat_force_scaling(&self) -> f64 {
        self.inner.bounds_mat_force_scaling
    }

    #[setter]
    fn set_bounds_mat_force_scaling(&mut self, value: f64) {
        self.inner.bounds_mat_force_scaling = value;
    }

    #[getter]
    fn track_failures(&self) -> bool {
        self.inner.track_failures
    }

    #[setter]
    fn set_track_failures(&mut self, value: bool) {
        self.inner.track_failures = value;
    }

    #[getter]
    fn failures(&self) -> Vec<u32> {
        self.inner.failures.clone()
    }

    #[getter]
    fn enable_sequential_random_seeds(&self) -> bool {
        self.inner.enable_sequential_random_seeds
    }

    #[setter]
    fn set_enable_sequential_random_seeds(&mut self, value: bool) {
        self.inner.enable_sequential_random_seeds = value;
    }

    #[getter]
    fn symmetrize_conjugated_terminal_groups_for_pruning(&self) -> bool {
        self.inner.symmetrize_conjugated_terminal_groups_for_pruning
    }

    #[setter]
    fn set_symmetrize_conjugated_terminal_groups_for_pruning(&mut self, value: bool) {
        self.inner.symmetrize_conjugated_terminal_groups_for_pruning = value;
    }

    fn update_from_json(&mut self, json: &str) -> PyResult<()> {
        self.inner.update_from_json(json).map_err(|err| {
            PyValueError::new_err(format!("EmbedParameters.update_from_json failed: {err}"))
        })
    }

    fn to_json(&self) -> String {
        self.inner.to_json()
    }

    fn __repr__(&self) -> String {
        format!(
            "EmbedParameters(random_seed={}, num_threads={}, prune_rms_thresh={}, clear_confs={})",
            self.inner.random_seed,
            self.inner.num_threads,
            self.inner.prune_rms_thresh,
            self.inner.clear_confs
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
#[doc = r#"
A protein-focused structural value.

``Protein`` is the default high-level protein API. It keeps amino-acid
residues and excludes ligands, nucleic acids, and waters by default.

Use ``Protein.from_pdb()`` for PDB files, ``Protein.from_pdb_str()`` for PDB
text, ``Protein.from_mmcif()`` for mmCIF files, and
``Protein.from_mmcif_str()`` for mmCIF text.
"#]
struct Protein {
    inner: cosmolkit_core::Protein,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
struct ProteinChain {
    inner: cosmolkit_core::Protein,
    index: usize,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
struct ProteinResidue {
    inner: cosmolkit_core::Protein,
    index: usize,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "ResidueInfo", skip_from_py_object)]
#[derive(Clone)]
#[doc = r#"
Gemmi-derived tabulated residue information.

Use ``ResidueInfo.code()`` and ``ResidueInfo.kind()`` for enum matching instead
of matching raw residue-name strings.
"#]
struct PyResidueInfo {
    inner: cosmolkit_core::ResidueInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
struct ProteinAtom {
    inner: cosmolkit_core::Protein,
    index: usize,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
#[doc = r#"
Read-only atom feature record returned by ``Molecule.atoms()``.

The methods on this object expose common atom properties such as atomic number,
formal charge, aromaticity, chiral tag, hydrogen counts, and valence values.
"#]
struct Atom {
    idx: usize,
    atomic_num: usize,
    formal_charge: i8,
    chiral_tag_name: String,
    chiral_tag_code: i64,
    isotope: Option<u16>,
    atom_map_num: Option<u32>,
    is_aromatic: bool,
    explicit_hydrogens: usize,
    no_implicit: bool,
    num_radical_electrons: usize,
    degree: usize,
    explicit_valence: Option<usize>,
    implicit_hydrogens: Option<usize>,
    total_num_hs: Option<usize>,
    total_valence: Option<usize>,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
#[doc = r#"
Read-only bond feature record returned by ``Molecule.bonds()``.

The methods on this object expose atom endpoints, bond type, direction,
stereo labels, stereo atom indices, and aromaticity.
"#]
struct Bond {
    idx: usize,
    begin_atom_idx: usize,
    end_atom_idx: usize,
    bond_type_name: String,
    bond_type_code: i64,
    bond_dir_name: String,
    bond_dir_code: i64,
    stereo_name: String,
    stereo_code: i64,
    stereo_atoms: Vec<usize>,
    is_aromatic: bool,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "BatchError", skip_from_py_object)]
#[derive(Clone)]
#[doc = r#"
A per-record batch processing error.

Batch methods can keep invalid records when ``errors="keep"`` is used. In
that case, ``MoleculeBatch.errors()`` returns ``BatchError`` objects describing
the input index, operation, and message.
"#]
struct PyBatchError {
    index: usize,
    operation: String,
    message: String,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "SdfRecordMetadata", skip_from_py_object)]
#[derive(Clone)]
#[doc = r#"
Lightweight metadata for one indexed SDF record.

Metadata is available from ``SdfDataset`` without parsing the molecule graph.
"#]
struct PySdfRecordMetadata {
    inner: cosmolkit_core::SdfRecordMetadata,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "SdfRecord", skip_from_py_object)]
#[derive(Clone)]
#[doc = r#"
One parsed SDF record returned by ``SdfDataset``.

The record exposes the parsed molecule plus SDF data fields.
"#]
struct PySdfRecord {
    index: usize,
    molecule: cosmolkit_core::Molecule,
    data_fields: Vec<(String, String)>,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "SdfDataset", skip_from_py_object)]
#[derive(Clone)]
#[doc = r#"
Indexed, seekable SDF dataset.

``SdfDataset`` builds a lightweight in-memory index of record byte ranges first.
After opening, ``len(dataset)`` is cheap, ``dataset[i]`` parses only that record,
``dataset[:n]`` returns a ``MoleculeBatch``, and ``dataset.batches(size=...)``
yields bounded ``MoleculeBatch`` chunks.

Use ``MoleculeBatch.read_sdf()`` when you intentionally want the whole file in
memory. Use ``SdfDataset`` for large seekable files where random access,
metadata inspection, or chunked processing matter.
"#]
struct PySdfDataset {
    inner: cosmolkit_core::SdfDataset,
    coordinate_mode: cosmolkit_core::SdfCoordinateMode,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "SdfReader", skip_from_py_object)]
#[derive(Clone)]
#[doc = r#"
Forward-only SDF reader for one-pass workflows.

Use ``SdfReader`` for non-indexed stream-style processing. For seekable files
where random access or accurate record-count progress matters, prefer
``SdfDataset``.
"#]
struct PySdfReader {
    path: PathBuf,
    coordinate_mode: cosmolkit_core::io::sdf::SdfCoordinateMode,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "SdfDatasetIterator", skip_from_py_object)]
struct PySdfDatasetIterator {
    dataset: cosmolkit_core::SdfDataset,
    coordinate_mode: cosmolkit_core::SdfCoordinateMode,
    position: usize,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "SdfBatchIterator", skip_from_py_object)]
struct PySdfBatchIterator {
    dataset: cosmolkit_core::SdfDataset,
    coordinate_mode: cosmolkit_core::SdfCoordinateMode,
    position: usize,
    size: usize,
    errors: BatchErrorMode,
    n_jobs: Option<usize>,
    progress_bar: Option<cosmolkit_core::BatchProgressBar>,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "SdfReaderBatchIterator", skip_from_py_object)]
struct PySdfReaderBatchIterator {
    reader: cosmolkit_core::io::sdf::SdfReader<BufReader<File>>,
    index: usize,
    size: usize,
    errors: BatchErrorMode,
    n_jobs: Option<usize>,
}

impl From<BatchRecordError> for PyBatchError {
    fn from(error: BatchRecordError) -> Self {
        Self {
            index: error.index,
            operation: error.operation.to_string(),
            message: error.message,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "BatchExportReport", skip_from_py_object)]
#[derive(Clone)]
#[doc = r#"
Summary returned by batch export methods.

The report records how many inputs were processed successfully and includes
structured errors for records that could not be exported.
"#]
struct PyBatchExportReport {
    written: usize,
    skipped: usize,
    failed: usize,
    errors: Vec<PyBatchError>,
}

impl From<cosmolkit_core::BatchExportReport> for PyBatchExportReport {
    fn from(report: cosmolkit_core::BatchExportReport) -> Self {
        Self {
            written: report.written,
            skipped: report.skipped,
            failed: report.failed,
            errors: report.errors.into_iter().map(Into::into).collect(),
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
#[doc = r#"
An ordered collection of molecules for batch workflows.

``MoleculeBatch`` keeps input order and supports construction, transformation,
filtering, rendering, and SDF export across many molecules. Methods that
transform molecules return a new batch.

Parameters such as ``errors`` control invalid-record handling:

- ``"raise"`` raises an exception when any record fails.
- ``"keep"`` keeps failed records and exposes them through
  ``errors()``. Export methods write valid records and count invalid records as
  skipped in the returned report.

Examples
--------
Construct a batch with ``MoleculeBatch.from_smiles_list()``, choose an
``errors`` mode for invalid records, and use ``with_parallel_jobs()`` when the
same worker count should apply to later batch operations.
"#]
struct MoleculeBatch {
    inner: cosmolkit_core::MoleculeBatch,
}

impl MoleculeBatch {
    fn with_inner(&self, inner: cosmolkit_core::MoleculeBatch) -> Self {
        Self { inner }
    }

    fn records_as_molecules(&self) -> Vec<Option<Molecule>> {
        (0..self.inner.len())
            .map(|index| match self.inner.get(index) {
                Some(BatchRecord::Molecule(molecule)) => Some(Molecule {
                    inner: molecule.clone(),
                }),
                _ => None,
            })
            .collect()
    }

    fn normalize_index(&self, index: isize) -> PyResult<usize> {
        let len = self.inner.len() as isize;
        let index = if index < 0 { len + index } else { index };
        if index < 0 || index >= len {
            return Err(PyIndexError::new_err("MoleculeBatch index out of range"));
        }
        Ok(index as usize)
    }

    fn select_records(&self, indices: &[usize]) -> Self {
        let records = indices
            .iter()
            .filter_map(|index| self.inner.get(*index).cloned())
            .collect();
        let inner = cosmolkit_core::MoleculeBatch::new(records)
            .with_parallel_jobs(self.inner.parallel_jobs())
            .with_progress_bar(self.inner.progress_bar());
        Self { inner }
    }

    fn selected_batch_pyobject(&self, py: Python<'_>, indices: &[usize]) -> PyResult<Py<PyAny>> {
        Ok(Py::new(py, self.select_records(indices))?.into_any())
    }

    fn get_record_pyobject(&self, py: Python<'_>, index: usize) -> PyResult<Py<PyAny>> {
        match self.inner.get(index) {
            Some(BatchRecord::Molecule(molecule)) => Ok(Py::new(
                py,
                Molecule {
                    inner: molecule.clone(),
                },
            )?
            .into_any()),
            _ => Ok(py.None()),
        }
    }

    fn slice_indices(&self, slice: &Bound<'_, PySlice>) -> PyResult<Vec<usize>> {
        let indices = slice.indices(self.inner.len() as isize)?;
        let mut out = Vec::with_capacity(indices.slicelength);
        let mut index = indices.start;
        for _ in 0..indices.slicelength {
            out.push(index as usize);
            index += indices.step;
        }
        Ok(out)
    }

    fn sequence_indices(&self, key: &Bound<'_, PyAny>) -> PyResult<Vec<usize>> {
        let items = key.extract::<Vec<Py<PyAny>>>()?;
        if items.is_empty() {
            return Ok(Vec::new());
        }

        let py = key.py();
        let bool_mask = items
            .iter()
            .all(|item| item.bind(py).is_exact_instance_of::<PyBool>());
        if bool_mask {
            if items.len() != self.inner.len() {
                return Err(PyIndexError::new_err(format!(
                    "boolean mask length {} does not match MoleculeBatch length {}",
                    items.len(),
                    self.inner.len()
                )));
            }
            let mut indices = Vec::new();
            for (index, item) in items.iter().enumerate() {
                if item.bind(py).extract::<bool>()? {
                    indices.push(index);
                }
            }
            return Ok(indices);
        }

        let mut indices = Vec::with_capacity(items.len());
        for item in items {
            let item = item.bind(py);
            if item.is_exact_instance_of::<PyBool>() {
                return Err(PyTypeError::new_err(
                    "MoleculeBatch index lists must not mix bool and int values",
                ));
            }
            indices.push(self.normalize_index(item.extract::<isize>()?)?);
        }
        Ok(indices)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyBatchError {
    #[doc = r#"
Return the zero-based input index that produced the error.
"#]
    fn index(&self) -> usize {
        self.index
    }

    #[doc = r#"
Return the operation name.
"#]
    fn operation(&self) -> String {
        self.operation.clone()
    }

    #[doc = r#"
Return the human-readable error message.
"#]
    fn message(&self) -> String {
        self.message.clone()
    }

    #[doc = r#"
Return the error as key-value pairs.
"#]
    fn as_dict(&self) -> Vec<(String, String)> {
        vec![
            ("index".to_string(), self.index.to_string()),
            ("operation".to_string(), self.operation.clone()),
            ("message".to_string(), self.message.clone()),
        ]
    }

    fn __repr__(&self) -> String {
        format!(
            "BatchError(index={}, operation='{}', message='{}')",
            self.index, self.operation, self.message
        )
    }
}

fn sdf_record_py(index: usize, record: cosmolkit_core::io::sdf::SdfRecord) -> PySdfRecord {
    PySdfRecord {
        index,
        molecule: record.molecule,
        data_fields: record.data_fields,
    }
}

fn sdf_indices_from_key(len: usize, key: &Bound<'_, PyAny>) -> PyResult<Result<Vec<usize>, usize>> {
    if key.is_exact_instance_of::<PyBool>() {
        return Err(PyTypeError::new_err(
            "SdfDataset scalar boolean indices are not supported; use an integer index, slice, integer list, or boolean mask sequence",
        ));
    }
    if let Ok(index) = key.extract::<isize>() {
        let len_i = len as isize;
        let index = if index < 0 { len_i + index } else { index };
        if index < 0 || index >= len_i {
            return Err(PyIndexError::new_err("SdfDataset index out of range"));
        }
        return Ok(Err(index as usize));
    }
    if let Ok(slice) = key.cast::<PySlice>() {
        let indices = slice.indices(len as isize)?;
        let mut out = Vec::with_capacity(indices.slicelength);
        let mut index = indices.start;
        for _ in 0..indices.slicelength {
            out.push(index as usize);
            index += indices.step;
        }
        return Ok(Ok(out));
    }
    let items = key.extract::<Vec<Py<PyAny>>>()?;
    if items.is_empty() {
        return Ok(Ok(Vec::new()));
    }

    let py = key.py();
    let bool_mask = items
        .iter()
        .all(|item| item.bind(py).is_exact_instance_of::<PyBool>());
    if bool_mask {
        if items.len() != len {
            return Err(PyIndexError::new_err(format!(
                "boolean mask length {} does not match SdfDataset length {}",
                items.len(),
                len
            )));
        }
        let mut out = Vec::new();
        for (index, item) in items.iter().enumerate() {
            if item.bind(py).extract::<bool>()? {
                out.push(index);
            }
        }
        return Ok(Ok(out));
    }

    let mut out = Vec::with_capacity(items.len());
    for item in items {
        let item = item.bind(py);
        if item.is_exact_instance_of::<PyBool>() {
            return Err(PyTypeError::new_err(
                "SdfDataset index lists must not mix bool and int values",
            ));
        }
        let index = item.extract::<isize>()?;
        let len_i = len as isize;
        let index = if index < 0 { len_i + index } else { index };
        if index < 0 || index >= len_i {
            return Err(PyIndexError::new_err("SdfDataset index out of range"));
        }
        out.push(index as usize);
    }
    Ok(Ok(out))
}

fn sdf_batch_from_range(
    dataset: &cosmolkit_core::SdfDataset,
    coordinate_mode: cosmolkit_core::SdfCoordinateMode,
    start: usize,
    end: usize,
    errors: BatchErrorMode,
    progress_bar: Option<&cosmolkit_core::BatchProgressBar>,
) -> Result<cosmolkit_core::MoleculeBatch, cosmolkit_core::BatchValidationError> {
    let mut records = Vec::with_capacity(end.saturating_sub(start));
    let params = cosmolkit_core::SdfReadParams {
        coordinate_mode,
        ..Default::default()
    };
    for index in start..end {
        match dataset.record_with_params(index, params) {
            Ok(record) => records.push(BatchRecord::Molecule(record.molecule)),
            Err(error) => {
                let record_error = BatchRecordError::new(index, "read_sdf", error.to_string());
                match errors {
                    BatchErrorMode::Strict | BatchErrorMode::KeepErrors => {
                        records.push(BatchRecord::Error(record_error));
                    }
                }
            }
        }
        if let Some(progress_bar) = progress_bar {
            progress_bar.inc(1);
        }
    }
    cosmolkit_core::MoleculeBatch::from_records_with_mode(records, errors)
}

fn sdf_batch_from_indices(
    dataset: &cosmolkit_core::SdfDataset,
    coordinate_mode: cosmolkit_core::SdfCoordinateMode,
    indices: Vec<usize>,
) -> Result<cosmolkit_core::MoleculeBatch, cosmolkit_core::BatchValidationError> {
    let mut records = Vec::with_capacity(indices.len());
    let params = cosmolkit_core::SdfReadParams {
        coordinate_mode,
        ..Default::default()
    };
    for index in indices {
        match dataset.record_with_params(index, params) {
            Ok(record) => records.push(BatchRecord::Molecule(record.molecule)),
            Err(error) => records.push(BatchRecord::Error(BatchRecordError::new(
                index,
                "read_sdf",
                error.to_string(),
            ))),
        }
    }
    cosmolkit_core::MoleculeBatch::from_records_with_mode(records, BatchErrorMode::Strict)
}

fn sdf_batch_from_reader(
    reader: &mut cosmolkit_core::io::sdf::SdfReader<BufReader<File>>,
    start_index: usize,
    size: usize,
    errors: BatchErrorMode,
) -> Result<Option<(cosmolkit_core::MoleculeBatch, usize)>, cosmolkit_core::BatchValidationError> {
    let mut records = Vec::with_capacity(size);
    let mut seen = 0usize;
    let mut index = start_index;
    while seen < size {
        match reader.next_record() {
            Ok(Some(record)) => {
                records.push(BatchRecord::Molecule(record.molecule));
                seen += 1;
                index += 1;
            }
            Ok(None) => break,
            Err(error) => {
                let record_error = BatchRecordError::new(index, "read_sdf", error.to_string());
                match errors {
                    BatchErrorMode::Strict | BatchErrorMode::KeepErrors => {
                        records.push(BatchRecord::Error(record_error));
                    }
                }
                seen += 1;
                index += 1;
            }
        }
    }
    if seen == 0 {
        return Ok(None);
    }
    cosmolkit_core::MoleculeBatch::from_records_with_mode(records, errors)
        .map(|batch| Some((batch, index)))
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PySdfRecordMetadata {
    fn index(&self) -> usize {
        self.inner.index
    }

    fn byte_offset(&self) -> u64 {
        self.inner.byte_offset
    }

    fn byte_len(&self) -> u64 {
        self.inner.byte_len
    }

    fn byte_range(&self) -> (u64, u64) {
        (
            self.inner.byte_offset,
            self.inner.byte_offset + self.inner.byte_len,
        )
    }

    fn line_range(&self) -> (usize, usize) {
        (
            self.inner.line_offset,
            self.inner.line_offset + self.inner.line_len,
        )
    }

    fn title(&self) -> Option<String> {
        self.inner.title.clone()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PySdfRecord {
    fn index(&self) -> usize {
        self.index
    }

    fn title(&self) -> Option<String> {
        self.molecule.properties().name().map(ToOwned::to_owned)
    }

    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.molecule.clone(),
        }
    }

    fn data_fields(&self) -> Vec<(String, String)> {
        self.data_fields.clone()
    }

    fn data_field(&self, name: &str) -> Option<String> {
        self.data_fields
            .iter()
            .find_map(|(key, value)| (key == name).then(|| value.clone()))
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PySdfDataset {
    #[classmethod]
    #[pyo3(signature = (path, index=None, build=None, coordinate_dim="auto"))]
    fn open(
        _cls: &Bound<'_, PyType>,
        path: &str,
        index: Option<&Bound<'_, PyAny>>,
        build: Option<&str>,
        coordinate_dim: &str,
    ) -> PyResult<Self> {
        let expanded_path = expand_user_path(path)?;
        let coordinate_mode = parse_coordinate_mode(Some(coordinate_dim))?;
        if let Some(build) = build
            && !matches!(build, "auto" | "always" | "never")
        {
            return Err(PyValueError::new_err(
                "unsupported build mode, expected one of: auto, always, never",
            ));
        }
        if let Some(index) = index
            && !index.is_none()
        {
            if let Ok(value) = index.extract::<String>() {
                if value != "auto" && value != "memory" {
                    return Err(PyNotImplementedError::new_err(
                        "persistent SDF sidecar indexes are not implemented yet; use index='auto', index='memory', or None",
                    ));
                }
            }
        }
        let inner = cosmolkit_core::SdfDataset::open_with_params(
            expanded_path,
            cosmolkit_core::SdfReadParams {
                coordinate_mode,
                ..Default::default()
            },
        )
        .map_err(|err| PyValueError::new_err(format!("SdfDataset.open failed: {err}")))?;
        Ok(Self {
            inner,
            coordinate_mode,
        })
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn path(&self) -> String {
        self.inner.path().display().to_string()
    }

    fn metadata(&self, index: isize) -> PyResult<PySdfRecordMetadata> {
        let len = self.inner.len() as isize;
        let index = if index < 0 { len + index } else { index };
        if index < 0 || index >= len {
            return Err(PyIndexError::new_err("SdfDataset index out of range"));
        }
        Ok(PySdfRecordMetadata {
            inner: self
                .inner
                .metadata(index as usize)
                .ok_or_else(|| PyIndexError::new_err("SdfDataset index out of range"))?
                .clone(),
        })
    }

    #[gen_stub(override_return_type(type_repr = "typing.Union[SdfRecord, MoleculeBatch]", imports = ("typing")))]
    fn __getitem__(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match sdf_indices_from_key(self.inner.len(), key)? {
            Err(index) => {
                let record = self
                    .inner
                    .record_with_params(
                        index,
                        cosmolkit_core::SdfReadParams {
                            coordinate_mode: self.coordinate_mode,
                            ..Default::default()
                        },
                    )
                    .map_err(|err| {
                        PyValueError::new_err(format!("SdfDataset read failed: {err}"))
                    })?;
                Ok(Py::new(py, sdf_record_py(index, record))?.into_any())
            }
            Ok(indices) => {
                let inner = sdf_batch_from_indices(&self.inner, self.coordinate_mode, indices)
                    .map_err(batch_validation_pyerr)?;
                Ok(Py::new(py, MoleculeBatch { inner })?.into_any())
            }
        }
    }

    fn __iter__(&self) -> PySdfDatasetIterator {
        PySdfDatasetIterator {
            dataset: self.inner.clone(),
            coordinate_mode: self.coordinate_mode,
            position: 0,
        }
    }

    #[pyo3(signature = (size=1024, indices=None, errors=None, n_jobs=None, progress_bar=false))]
    fn batches(
        &self,
        size: usize,
        indices: Option<Vec<usize>>,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: bool,
    ) -> PyResult<PySdfBatchIterator> {
        if size == 0 {
            return Err(PyValueError::new_err("size must be >= 1"));
        }
        if indices.is_some() {
            return Err(PyNotImplementedError::new_err(
                "indices=... for SdfDataset.batches() is not implemented yet",
            ));
        }
        let errors = parse_batch_error_mode(errors)?;
        let n_jobs = validate_n_jobs(n_jobs)?;
        let progress_bar =
            maybe_batch_progress_bar(progress_bar, self.inner.len(), "read_sdf_batches");
        Ok(PySdfBatchIterator {
            dataset: self.inner.clone(),
            coordinate_mode: self.coordinate_mode,
            position: 0,
            size,
            errors,
            n_jobs,
            progress_bar,
        })
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PySdfReader {
    #[classmethod]
    #[pyo3(signature = (path, coordinate_dim="auto"))]
    fn open(_cls: &Bound<'_, PyType>, path: &str, coordinate_dim: &str) -> PyResult<Self> {
        Ok(Self {
            path: expand_user_path(path)?,
            coordinate_mode: parse_coordinate_mode(Some(coordinate_dim))?,
        })
    }

    #[pyo3(signature = (size=1024, errors=None, n_jobs=None, progress_bar=false))]
    fn batches(
        &self,
        size: usize,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: bool,
    ) -> PyResult<PySdfReaderBatchIterator> {
        if size == 0 {
            return Err(PyValueError::new_err("size must be >= 1"));
        }
        if progress_bar {
            return Err(PyNotImplementedError::new_err(
                "SdfReader.batches(progress_bar=True) cannot show an accurate total for forward-only streams; use SdfDataset.batches() for indexed progress",
            ));
        }
        let file = File::open(&self.path)
            .map_err(|err| PyValueError::new_err(format!("SdfReader.open failed: {err}")))?;
        Ok(PySdfReaderBatchIterator {
            reader: cosmolkit_core::io::sdf::SdfReader::with_coordinate_mode(
                BufReader::new(file),
                self.coordinate_mode,
            ),
            index: 0,
            size,
            errors: parse_batch_error_mode(errors)?,
            n_jobs: validate_n_jobs(n_jobs)?,
        })
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PySdfDatasetIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    #[gen_stub(override_return_type(type_repr = "SdfRecord"))]
    fn __next__(&mut self) -> PyResult<Option<PySdfRecord>> {
        if self.position >= self.dataset.len() {
            return Ok(None);
        }
        let index = self.position;
        self.position += 1;
        let record = self
            .dataset
            .record_with_params(
                index,
                cosmolkit_core::SdfReadParams {
                    coordinate_mode: self.coordinate_mode,
                    ..Default::default()
                },
            )
            .map_err(|err| PyValueError::new_err(format!("SdfDataset read failed: {err}")))?;
        Ok(Some(sdf_record_py(index, record)))
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PySdfBatchIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    #[gen_stub(override_return_type(type_repr = "MoleculeBatch"))]
    fn __next__(&mut self) -> PyResult<Option<MoleculeBatch>> {
        if self.position >= self.dataset.len() {
            if let Some(progress_bar) = self.progress_bar.take() {
                progress_bar.finish();
            }
            return Ok(None);
        }
        let start = self.position;
        let end = self.dataset.len().min(start + self.size);
        self.position = end;
        let inner = sdf_batch_from_range(
            &self.dataset,
            self.coordinate_mode,
            start,
            end,
            self.errors,
            self.progress_bar.as_ref(),
        )
        .map_err(batch_validation_pyerr)?;
        Ok(Some(MoleculeBatch {
            inner: inner.with_parallel_jobs(self.n_jobs),
        }))
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PySdfReaderBatchIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    #[gen_stub(override_return_type(type_repr = "MoleculeBatch"))]
    fn __next__(&mut self) -> PyResult<Option<MoleculeBatch>> {
        let start = self.index;
        let result = sdf_batch_from_reader(&mut self.reader, start, self.size, self.errors)
            .map_err(batch_validation_pyerr)?;
        let Some((inner, next_index)) = result else {
            return Ok(None);
        };
        self.index = next_index;
        Ok(Some(MoleculeBatch {
            inner: inner.with_parallel_jobs(self.n_jobs),
        }))
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyBatchExportReport {
    #[doc = r#"
Return the total number of records considered for export.
"#]
    fn total(&self) -> usize {
        self.written + self.skipped + self.failed
    }

    #[doc = r#"
Return the number of records exported successfully.
"#]
    fn success(&self) -> usize {
        self.written
    }

    #[doc = r#"
Return the number of records that failed during export.
"#]
    fn failed(&self) -> usize {
        self.failed
    }

    #[doc = r#"
Return structured errors for failed records.
"#]
    fn errors(&self) -> Vec<PyBatchError> {
        self.errors.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "BatchExportReport(written={}, skipped={}, failed={})",
            self.written, self.skipped, self.failed
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl MoleculeBatch {
    #[classmethod]
    #[pyo3(signature = (smiles, sanitize=None, errors=None, n_jobs=None))]
    #[doc = r#"
Create a batch from a list of SMILES strings.

Parameters
----------
smiles : list[str]
    Input SMILES strings.
sanitize : bool, optional
    Optional molecule preparation flag. COSMolKit applies the available
    preparation behavior during construction.
errors : {"raise", "keep"}, optional
    Invalid-record handling mode. The default is ``"raise"``.
n_jobs : int, optional
    Number of worker threads to use. ``None`` uses the default scheduler.

Returns
-------
MoleculeBatch
    A batch preserving the input order for valid and kept records.
"#]
    fn from_smiles_list(
        _cls: &Bound<'_, PyType>,
        smiles: Vec<String>,
        sanitize: Option<bool>,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
    ) -> PyResult<Self> {
        let sanitize = sanitize.unwrap_or(true);
        let mode = parse_batch_error_mode(errors)?;
        let inner = cosmolkit_core::MoleculeBatch::from_smiles_list_with_sanitize_and_options(
            &smiles,
            sanitize,
            mode,
            validate_n_jobs(n_jobs)?,
            None,
        )
        .map_err(batch_validation_pyerr)?;
        Ok(Self { inner })
    }

    #[classmethod]
    #[pyo3(signature = (sdf_text, errors=None, n_jobs=None, coordinate_dim="auto"))]
    #[doc = r#"
Read all molecule records from an SDF string.

Parameters
----------
sdf_text : str
    SDF text containing one or more records.
errors : {"raise", "keep"}, optional
    Invalid-record handling mode. The default is ``"raise"``.
n_jobs : int, optional
    Number of worker threads to use.
coordinate_dim : {"auto", "2d", "3d"}, optional
    Coordinate interpretation mode. ``"auto"`` preserves the molfile header.
"#]
    fn read_sdf_records_from_str(
        _cls: &Bound<'_, PyType>,
        sdf_text: &str,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        coordinate_dim: &str,
    ) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let coordinate_mode = parse_coordinate_mode(Some(coordinate_dim))?;
        let inner = cosmolkit_core::MoleculeBatch::read_sdf_records_from_str_with_options(
            sdf_text,
            coordinate_mode,
            mode,
            validate_n_jobs(n_jobs)?,
            None,
        )
        .map_err(batch_validation_pyerr)?;
        Ok(Self { inner })
    }

    #[classmethod]
    #[pyo3(signature = (path, errors=None, n_jobs=None, progress_bar=false, coordinate_dim="auto"))]
    #[doc = r#"
Read all molecule records from an SDF file into a batch.

Parameters
----------
path : str
    SDF file path.
errors : {"raise", "keep"}, optional
    Invalid-record handling mode. The default is ``"raise"``.
n_jobs : int, optional
    Number of worker threads to use for batch construction.
progress_bar : bool, optional
    Show a Rust-side progress bar while records are parsed. This builds a
    lightweight record index first so the total is known.
coordinate_dim : {"auto", "2d", "3d"}, optional
    Coordinate interpretation mode. ``"auto"`` preserves the molfile header.
"#]
    fn read_sdf(
        _cls: &Bound<'_, PyType>,
        path: &str,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: bool,
        coordinate_dim: &str,
    ) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let coordinate_mode = parse_coordinate_mode(Some(coordinate_dim))?;
        let expanded_path = expand_user_path(path)?;
        if progress_bar {
            let dataset = cosmolkit_core::SdfDataset::open_with_params(
                &expanded_path,
                cosmolkit_core::SdfReadParams {
                    coordinate_mode,
                    ..Default::default()
                },
            )
            .map_err(|err| PyValueError::new_err(format!("read_sdf index failed: {err}")))?;
            let inner = cosmolkit_core::MoleculeBatch::read_sdf_dataset_with_params_and_options(
                &dataset,
                cosmolkit_core::SdfReadParams {
                    coordinate_mode,
                    ..Default::default()
                },
                mode,
                validate_n_jobs(n_jobs)?,
                Some(true),
            )
            .map_err(batch_validation_pyerr)?;
            return Ok(Self { inner });
        }
        let file = File::open(&expanded_path).map_err(|error| {
            PyValueError::new_err(format!(
                "read_sdf open failed for {}: {error}",
                expanded_path.display()
            ))
        })?;
        let inner = cosmolkit_core::MoleculeBatch::read_sdf_records_from_reader_with_options(
            BufReader::new(file),
            coordinate_mode,
            mode,
            validate_n_jobs(n_jobs)?,
            None,
        )
        .map_err(batch_validation_pyerr)?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (n_jobs))]
    #[doc = r#"
Return a new batch configured to use this worker count by default.

Pass ``None`` to clear the batch-level default and let rayon decide. Method-level
``n_jobs`` arguments still override this setting for that one call.
"#]
    fn with_parallel_jobs(&self, n_jobs: Option<usize>) -> PyResult<Self> {
        Ok(Self {
            inner: self
                .inner
                .clone()
                .with_parallel_jobs(validate_n_jobs(n_jobs)?),
        })
    }

    #[doc = r#"
Return the batch-level default worker count, or ``None`` when unset.
"#]
    fn parallel_jobs(&self) -> Option<usize> {
        self.inner.parallel_jobs()
    }

    #[pyo3(signature = (progress_bar))]
    #[doc = r#"
Return a new batch configured to show Rust-side progress bars by default.

Pass ``None`` to clear the batch-level default. Method-level ``progress_bar``
arguments still override this setting for that one call.
"#]
    fn with_progress_bar(&self, progress_bar: Option<bool>) -> Self {
        Self {
            inner: self.inner.clone().with_progress_bar(progress_bar),
        }
    }

    #[doc = r#"
Return the batch-level progress-bar default, or ``None`` when unset.
"#]
    fn progress_bar(&self) -> Option<bool> {
        self.inner.progress_bar()
    }

    #[pyo3(signature = (errors=None, n_jobs=None, progress_bar=None))]
    #[doc = r#"
Return a new batch with explicit hydrogens added to each valid molecule.
"#]
    fn with_hydrogens(
        &self,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let inner = self
            .inner
            .with_hydrogens_with_options(mode, validate_n_jobs(n_jobs)?, progress_bar)
            .map_err(batch_validation_pyerr)?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (errors=None, n_jobs=None, progress_bar=None))]
    #[doc = r#"
Return a new batch with explicit hydrogens removed from each valid molecule.
"#]
    fn without_hydrogens(
        &self,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let inner = self
            .inner
            .without_hydrogens_with_options(mode, validate_n_jobs(n_jobs)?, progress_bar)
            .map_err(batch_validation_pyerr)?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (strict=None, errors=None, n_jobs=None, progress_bar=None))]
    #[doc = r#"
Return a sanitized batch.

Parameters
----------
strict : bool, optional
    Optional strictness flag for available validation steps.
errors : {"raise", "keep"}, optional
    Invalid-record handling mode.
n_jobs : int, optional
    Number of worker threads to use.
"#]
    fn sanitize(
        &self,
        strict: Option<bool>,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        reject_non_strict_sanitize(strict)?;
        let mode = parse_batch_error_mode(errors)?;
        let inner = self
            .inner
            .sanitize_with_options(mode, validate_n_jobs(n_jobs)?, progress_bar)
            .map_err(batch_validation_pyerr)?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (clear_aromatic_flags=None, errors=None, n_jobs=None, progress_bar=None))]
    #[doc = r#"
Return a new batch with aromatic bonds converted to an explicit Kekule form.
"#]
    fn with_kekulized_bonds(
        &self,
        clear_aromatic_flags: Option<bool>,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let clear_aromatic_flags = clear_aromatic_flags.unwrap_or(true);
        let mode = parse_batch_error_mode(errors)?;
        let inner = self
            .inner
            .with_kekulized_bonds_with_options(
                clear_aromatic_flags,
                mode,
                validate_n_jobs(n_jobs)?,
                progress_bar,
            )
            .map_err(batch_validation_pyerr)?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (errors=None, n_jobs=None, progress_bar=None))]
    #[doc = r#"
Return a new batch with 2D coordinates computed for each valid molecule.
"#]
    fn with_2d_coordinates(
        &self,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let inner = self
            .inner
            .with_2d_coordinates_with_options(mode, validate_n_jobs(n_jobs)?, progress_bar)
            .map_err(batch_validation_pyerr)?;
        Ok(Self { inner })
    }

    #[doc = r#"
Return a batch containing only valid molecules.
"#]
    fn filter_valid(&self) -> Self {
        self.with_inner(self.inner.filter_valid())
    }

    #[doc = r#"
Return a boolean mask indicating which records are valid.
"#]
    fn valid_mask(&self) -> Vec<bool> {
        self.inner.valid_mask()
    }

    #[doc = r#"
Return a boolean mask indicating which records are invalid.
"#]
    fn invalid_mask(&self) -> Vec<bool> {
        self.inner
            .valid_mask()
            .into_iter()
            .map(|valid| !valid)
            .collect()
    }

    #[doc = r#"
Return structured errors collected for invalid records.
"#]
    fn errors(&self) -> Vec<PyBatchError> {
        self.inner
            .errors()
            .into_iter()
            .map(PyBatchError::from)
            .collect()
    }

    #[doc = r#"
Return the number of valid records.
"#]
    fn valid_count(&self) -> usize {
        self.inner.valid_count()
    }

    #[doc = r#"
Return the number of invalid records.
"#]
    fn invalid_count(&self) -> usize {
        self.inner.invalid_count()
    }

    #[pyo3(signature = (
        isomeric_smiles=true,
        canonical=true,
        kekule=false,
        clean_stereo=true,
        all_bonds_explicit=false,
        all_hs_explicit=false,
        include_dative_bonds=true,
        ignore_atom_map_numbers=false,
        rooted_at_atom=None,
        n_jobs=None,
        progress_bar=None
    ))]
    #[doc = r#"
Return one SMILES string per record.

Invalid records are returned as ``None`` when they are kept in the batch.

Parameters
----------
isomeric_smiles : bool, default True
    Include stereochemical and isotopic information when available.
canonical : bool, default True
    Return canonical SMILES when enabled.
kekule : bool, default False
    Write aromatic systems in Kekule form.
clean_stereo : bool, default True
    Normalize stereo output where possible.
all_bonds_explicit : bool, default False
    Write explicit bond symbols.
all_hs_explicit : bool, default False
    Write explicit hydrogens.
include_dative_bonds : bool, default True
    Include dative bond notation.
ignore_atom_map_numbers : bool, default False
    Omit atom map numbers from canonical decisions.
rooted_at_atom : int, optional
    Start traversal from a selected atom index.
n_jobs : int, optional
    Number of worker threads to use.
"#]
    #[allow(clippy::too_many_arguments)]
    fn to_smiles_list(
        &self,
        isomeric_smiles: bool,
        canonical: bool,
        kekule: bool,
        clean_stereo: bool,
        all_bonds_explicit: bool,
        all_hs_explicit: bool,
        include_dative_bonds: bool,
        ignore_atom_map_numbers: bool,
        rooted_at_atom: Option<usize>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<String>>> {
        let params = make_smiles_write_params(
            isomeric_smiles,
            canonical,
            kekule,
            clean_stereo,
            all_bonds_explicit,
            all_hs_explicit,
            include_dative_bonds,
            ignore_atom_map_numbers,
            rooted_at_atom,
        );
        self.inner
            .to_smiles_optional_list_with_params_and_options(
                &params,
                validate_n_jobs(n_jobs)?,
                progress_bar,
            )
            .map_err(batch_validation_pyerr)
    }

    #[pyo3(signature = (n_jobs=None, progress_bar=None))]
    #[doc = r#"
Return distance-geometry bounds matrices for all valid records.
"#]
    #[gen_stub(override_return_type(type_repr = "builtins.list[typing.Optional[numpy.ndarray[typing.Any, typing.Any]]]", imports = ("builtins", "typing", "numpy")))]
    fn dg_bounds_matrix_list<'py>(
        &self,
        py: Python<'py>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Bound<'py, PyList>> {
        let values = self
            .inner
            .dg_bounds_matrix_list_with_options(validate_n_jobs(n_jobs)?, progress_bar)
            .map_err(batch_validation_pyerr)?;
        let out = PyList::empty(py);
        for value in values {
            if let Some(matrix) = value {
                out.append(PyArray2::from_vec2(py, &matrix).map_err(|err| {
                    PyValueError::new_err(format!(
                        "MoleculeBatch.dg_bounds_matrix_list failed: {err}"
                    ))
                })?)?;
            } else {
                out.append(py.None())?;
            }
        }
        Ok(out)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        radius=2,
        n_bits=2048,
        include_chirality=false,
        use_bond_types=true,
        count_simulation=false,
        count_bounds=None,
        only_nonzero_invariants=false,
        include_redundant_environments=false,
        from_atoms=None,
        ignore_atoms=None,
        custom_atom_invariants=None,
        custom_bond_invariants=None,
        atom_invariants_generator=None,
        atom_invariants_include_ring_membership=true,
        bond_invariants_generator=None,
        bond_invariants_use_bond_types=true,
        bond_invariants_use_chirality=false,
        num_bits_per_feature=1,
        n_jobs=None,
        progress_bar=None
    ))]
    #[doc = r#"
Return Morgan fingerprints for valid batch records.

Invalid records are returned as ``None`` in their original positions.
"#]
    fn fingerprint_morgan_list(
        &self,
        radius: u32,
        n_bits: usize,
        include_chirality: bool,
        use_bond_types: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        only_nonzero_invariants: bool,
        include_redundant_environments: bool,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        atom_invariants_generator: Option<&str>,
        atom_invariants_include_ring_membership: bool,
        bond_invariants_generator: Option<&str>,
        bond_invariants_use_bond_types: bool,
        bond_invariants_use_chirality: bool,
        num_bits_per_feature: u32,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<Fingerprint>>> {
        let params = make_morgan_fingerprint_params(
            radius,
            n_bits,
            include_chirality,
            use_bond_types,
            count_simulation,
            count_bounds,
            only_nonzero_invariants,
            include_redundant_environments,
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            atom_invariants_generator,
            atom_invariants_include_ring_membership,
            bond_invariants_generator,
            bond_invariants_use_bond_types,
            bond_invariants_use_chirality,
            num_bits_per_feature,
            false,
        )?;
        self.inner
            .morgan_fingerprint_list_with_options(&params, validate_n_jobs(n_jobs)?, progress_bar)
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| Fingerprint { inner }))
                    .collect()
            })
            .map_err(batch_validation_pyerr)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        radius=2,
        n_bits=2048,
        include_chirality=false,
        use_bond_types=true,
        count_simulation=false,
        count_bounds=None,
        only_nonzero_invariants=false,
        include_redundant_environments=false,
        from_atoms=None,
        ignore_atoms=None,
        custom_atom_invariants=None,
        custom_bond_invariants=None,
        atom_invariants_generator=None,
        atom_invariants_include_ring_membership=true,
        bond_invariants_generator=None,
        bond_invariants_use_bond_types=true,
        bond_invariants_use_chirality=false,
        num_bits_per_feature=1,
        n_jobs=None,
        progress_bar=None
    ))]
    #[doc = r#"
Return Morgan fingerprints and additional output for valid batch records.

Invalid records are returned as ``None`` in their original positions.
"#]
    fn fingerprint_morgan_with_output_list(
        &self,
        radius: u32,
        n_bits: usize,
        include_chirality: bool,
        use_bond_types: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        only_nonzero_invariants: bool,
        include_redundant_environments: bool,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        atom_invariants_generator: Option<&str>,
        atom_invariants_include_ring_membership: bool,
        bond_invariants_generator: Option<&str>,
        bond_invariants_use_bond_types: bool,
        bond_invariants_use_chirality: bool,
        num_bits_per_feature: u32,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<MorganFingerprintResult>>> {
        let params = make_morgan_fingerprint_params(
            radius,
            n_bits,
            include_chirality,
            use_bond_types,
            count_simulation,
            count_bounds,
            only_nonzero_invariants,
            include_redundant_environments,
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            atom_invariants_generator,
            atom_invariants_include_ring_membership,
            bond_invariants_generator,
            bond_invariants_use_bond_types,
            bond_invariants_use_chirality,
            num_bits_per_feature,
            true,
        )?;
        self.inner
            .morgan_fingerprint_with_output_list_with_options(
                &params,
                validate_n_jobs(n_jobs)?,
                progress_bar,
            )
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(Into::into))
                    .collect()
            })
            .map_err(batch_validation_pyerr)
    }

    #[pyo3(signature = (width=300, height=300, n_jobs=None, progress_bar=None))]
    #[doc = r#"
Render each valid molecule to an SVG string.
"#]
    fn to_svg_list(
        &self,
        width: u32,
        height: u32,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<String>>> {
        self.inner
            .to_svg_list_with_options(width, height, validate_n_jobs(n_jobs)?, progress_bar)
            .map_err(batch_validation_pyerr)
    }

    #[pyo3(signature = (out_dir, format=None, size=None, n_jobs=None, errors=None, report_path=None, filenames=None, progress_bar=None))]
    #[doc = r#"
Write molecule depictions to a directory.

Parameters
----------
out_dir : str
    Output directory.
format : {"png", "svg"}, optional
    Image format. The default is ``"png"``.
size : tuple[int, int], optional
    Output image size as ``(width, height)``.
n_jobs : int, optional
    Number of worker threads to use.
errors : {"raise", "keep"}, optional
    Export error handling mode.
report_path : str, optional
    Write a JSON or CSV error report.
filenames : list[str | None], optional
    Per-record output filenames. Names are relative to ``out_dir``; missing
    extensions are filled from ``format``.

Returns
-------
BatchExportReport
    Export summary.
"#]
    fn to_images(
        &self,
        out_dir: &str,
        format: Option<&str>,
        size: Option<(u32, u32)>,
        n_jobs: Option<usize>,
        errors: Option<&Bound<'_, PyAny>>,
        report_path: Option<&str>,
        filenames: Option<Vec<Option<String>>>,
        progress_bar: Option<bool>,
    ) -> PyResult<PyBatchExportReport> {
        let mode = parse_batch_error_mode(errors)?;
        let image_format = format.unwrap_or("png").to_string();
        let (width, height) = size.unwrap_or((300, 300));
        let out_dir = expand_user_path(out_dir)?;
        let filenames = complete_batch_filenames(filenames, self.inner.len(), &image_format)?;
        let report = self
            .inner
            .write_images_with_options(
                out_dir.as_path(),
                &image_format,
                width,
                height,
                mode,
                filenames.as_deref(),
                validate_n_jobs(n_jobs)?,
                progress_bar,
            )
            .map_err(batch_validation_pyerr)?;
        if let Some(path) = report_path {
            write_batch_report(path, &report)?;
        }
        Ok(report.into())
    }

    #[pyo3(signature = (path, format=None, errors=None, n_jobs=None, report_path=None, progress_bar=None))]
    #[doc = r#"
Write valid molecules to an SDF file.

Parameters
----------
path : str
    Output SDF path.
format : {"auto", "v2000", "v3000"}, optional
    SDF output format.
errors : {"raise", "keep"}, optional
    Export error handling mode.
n_jobs : int, optional
    Number of worker threads to use.
report_path : str, optional
    Write a JSON or CSV error report.
"#]
    fn to_sdf(
        &self,
        path: &str,
        format: Option<&str>,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        report_path: Option<&str>,
        progress_bar: Option<bool>,
    ) -> PyResult<PyBatchExportReport> {
        let mode = parse_batch_error_mode(errors)?;
        let sdf_format = parse_sdf_format(format)?;
        let path = expand_user_path(path)?;
        let report = self
            .inner
            .write_sdf_with_options(
                path.as_path(),
                sdf_format,
                mode,
                validate_n_jobs(n_jobs)?,
                progress_bar,
            )
            .map_err(batch_validation_pyerr)?;
        if let Some(report_path) = report_path {
            write_batch_report(report_path, &report)?;
        }
        Ok(report.into())
    }

    #[pyo3(signature = (out_dir, format=None, errors=None, n_jobs=None, report_path=None, filenames=None, progress_bar=None))]
    #[doc = r#"
Write each valid molecule to its own SDF file in a directory.

Parameters
----------
out_dir : str
    Output directory.
format : {"auto", "v2000", "v3000"}, optional
    SDF output format.
errors : {"raise", "keep"}, optional
    Export error handling mode.
n_jobs : int, optional
    Number of worker threads to use.
report_path : str, optional
    Write a JSON or CSV error report.
filenames : list[str | None], optional
    Per-record output filenames. Names are relative to ``out_dir``; missing
    extensions are written as ``.sdf``.
"#]
    fn to_sdf_files(
        &self,
        out_dir: &str,
        format: Option<&str>,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        report_path: Option<&str>,
        filenames: Option<Vec<Option<String>>>,
        progress_bar: Option<bool>,
    ) -> PyResult<PyBatchExportReport> {
        let mode = parse_batch_error_mode(errors)?;
        let sdf_format = parse_sdf_format(format)?;
        let out_dir = expand_user_path(out_dir)?;
        let filenames = complete_batch_filenames(filenames, self.inner.len(), "sdf")?;
        let report = self
            .inner
            .write_sdf_files_with_options(
                out_dir.as_path(),
                sdf_format,
                mode,
                filenames.as_deref(),
                validate_n_jobs(n_jobs)?,
                progress_bar,
            )
            .map_err(batch_validation_pyerr)?;
        if let Some(report_path) = report_path {
            write_batch_report(report_path, &report)?;
        }
        Ok(report.into())
    }

    #[doc = r#"
Return batch records as a Python list.

Valid records become ``Molecule`` objects and invalid records become ``None``.
"#]
    fn to_list(&self) -> Vec<Option<Molecule>> {
        self.records_as_molecules()
    }

    fn __getitem__(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        if let Ok(slice) = key.cast::<PySlice>() {
            let indices = self.slice_indices(slice)?;
            return self.selected_batch_pyobject(py, &indices);
        }
        if !key.is_exact_instance_of::<PyBool>() {
            if let Ok(index) = key.extract::<isize>() {
                return self.get_record_pyobject(py, self.normalize_index(index)?);
            }
        }
        match self.sequence_indices(key) {
            Ok(indices) => return self.selected_batch_pyobject(py, &indices),
            Err(error) => {
                if key.extract::<Vec<Py<PyAny>>>().is_ok() {
                    return Err(error);
                }
            }
        }
        Err(PyTypeError::new_err(
            "MoleculeBatch indices must be integers, slices, integer lists, or boolean masks",
        ))
    }

    #[gen_stub(override_return_type(type_repr = "typing.Iterator[typing.Optional[Molecule]]", imports = ("typing")))]
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let list = PyList::new(py, self.records_as_molecules())?;
        Ok(PyIterator::from_object(list.as_any())?.into_any())
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "MoleculeBatch(n={}, valid={}, invalid={}, parallel_jobs={:?}, progress_bar={:?})",
            self.inner.len(),
            self.inner.valid_count(),
            self.inner.invalid_count(),
            self.inner.parallel_jobs(),
            self.inner.progress_bar()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl Molecule {
    #[classmethod]
    #[pyo3(signature = (smiles, sanitize=None))]
    #[doc = r#"
Create a molecule from a SMILES string.

Parameters
----------
smiles : str
    Input SMILES string.
sanitize : bool, optional
    Optional molecule preparation flag. COSMolKit applies the available
    preparation behavior during construction.

Returns
-------
Molecule
    Parsed molecule.

Examples
--------
Use ``Molecule.from_smiles("CCO")`` to create a molecule and
``mol.to_smiles()`` to write it back.
"#]
    fn from_smiles(
        _cls: &Bound<'_, PyType>,
        smiles: &str,
        sanitize: Option<bool>,
    ) -> PyResult<Self> {
        let mol =
            cosmolkit_core::Molecule::from_smiles_with_sanitize(smiles, sanitize.unwrap_or(true))
                .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self { inner: mol })
    }

    #[classmethod]
    #[pyo3(signature = (text, *, sanitize=true, remove_hs=true, flavor=0, proximity_bonding=true))]
    #[doc = r#"
Create a molecule from a PDB block.

This follows the COSMolKit core PDB molecule conversion profile, which is
designed to match RDKit ``Chem.MolFromPDBBlock`` for modeled molecule state.
Structural parsing is handled by COSMolKit's structure reader before molecule
conversion.

Parameters
----------
text : str
    PDB block text.
sanitize : bool
    Whether to sanitize after PDB molecule construction.
remove_hs : bool
    Whether sanitization should remove hydrogens.
flavor : int
    RDKit-compatible PDB parser flavor bit mask.
proximity_bonding : bool
    Whether to add proximity bonds using RDKit's PDB proximity-bond algorithm.

Returns
-------
Molecule
    Parsed molecule.
"#]
    fn from_pdb_block(
        _cls: &Bound<'_, PyType>,
        text: &str,
        sanitize: bool,
        remove_hs: bool,
        flavor: u32,
        proximity_bonding: bool,
    ) -> PyResult<Self> {
        let profile = cosmolkit_core::RdkitPdbMolProfile {
            sanitize,
            remove_hs,
            flavor,
            proximity_bonding,
        };
        let inner = cosmolkit_core::Molecule::from_pdb_block_with_options(text, profile)
            .map_err(pdb_molecule_pyerr)?;
        Ok(Self { inner })
    }

    #[classmethod]
    #[pyo3(signature = (text, *, sanitize=true, remove_hs=true, flavor=0, proximity_bonding=true))]
    #[doc = r#"
Create a molecule from an mmCIF block.

This uses COSMolKit's mmCIF structure reader, then applies the same
RDKit-compatible molecule conversion profile used by ``Molecule.from_pdb_block``.
RDKit does not provide a direct ``Chem.MolFromMMCIFBlock`` oracle; this API is a
COSMolKit mmCIF structural reader layered into the RDKit-compatible PDB molecule
conversion state.

Parameters
----------
text : str
    mmCIF block text.
sanitize : bool
    Whether to sanitize after molecule construction.
remove_hs : bool
    Whether sanitization should remove hydrogens.
flavor : int
    RDKit-compatible PDB parser flavor bit mask applied during molecule
    conversion.
proximity_bonding : bool
    Whether to add proximity bonds using RDKit's PDB proximity-bond algorithm.

Returns
-------
Molecule
    Parsed molecule.
"#]
    fn from_mmcif_block(
        _cls: &Bound<'_, PyType>,
        text: &str,
        sanitize: bool,
        remove_hs: bool,
        flavor: u32,
        proximity_bonding: bool,
    ) -> PyResult<Self> {
        let profile = cosmolkit_core::RdkitPdbMolProfile {
            sanitize,
            remove_hs,
            flavor,
            proximity_bonding,
        };
        let inner = cosmolkit_core::Molecule::from_mmcif_block_with_options(text, profile)
            .map_err(pdb_molecule_pyerr)?;
        Ok(Self { inner })
    }

    #[classmethod]
    #[doc = r#"
Create a molecule from an XYZ block.

XYZ contains atom identities and Cartesian coordinates only. This follows
COSMolKit core's RDKit-aligned ``MolFromXYZBlock`` behavior: atoms and one 3D
conformer are parsed, and bonds are not inferred.

The returned molecule is coordinate-only. Topology-dependent operations such as
adding hydrogens or ETKDG conformer generation require a trusted bond graph.

Parameters
----------
text : str
    XYZ block text.

Returns
-------
Molecule
    Parsed molecule with zero bonds and a 3D conformer when the atom count is
    nonzero.
"#]
    fn from_xyz_block(_cls: &Bound<'_, PyType>, text: &str) -> PyResult<Self> {
        let inner = cosmolkit_core::Molecule::from_xyz_block(text)
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self { inner })
    }

    #[classmethod]
    #[pyo3(signature = (rdmol, sanitize=None))]
    #[doc = r#"
Create a molecule from an RDKit molecule object.

Parameters
----------
rdmol : object
    An object exposing RDKit's molecule API.
sanitize : bool, optional
    Optional molecule preparation flag.

Returns
-------
Molecule
    COSMolKit molecule copied from the input object.
"#]
    fn from_rdkit(
        _cls: &Bound<'_, PyType>,
        rdmol: &Bound<'_, PyAny>,
        sanitize: Option<bool>,
    ) -> PyResult<Self> {
        rdmol.py().import("rdkit.Chem").map_err(|err| {
            PyValueError::new_err(format!(
                "Molecule.from_rdkit requires rdkit to be installed and importable: {err}"
            ))
        })?;
        let atom_count: usize = py_method_extract(rdmol, "GetNumAtoms")?;
        let bond_count: usize = py_method_extract(rdmol, "GetNumBonds")?;
        let mut builder = cosmolkit_core::MoleculeBuilder::new();

        for idx in 0..atom_count {
            let atom = py_method_index(rdmol, "GetAtomWithIdx", idx)?;
            let atomic_num_raw: usize = py_method_extract(&atom, "GetAtomicNum")?;
            let atomic_num = u8::try_from(atomic_num_raw).map_err(|_| {
                PyValueError::new_err(format!(
                    "from_rdkit atom {idx} atomic number out of u8 range: {atomic_num_raw}"
                ))
            })?;
            let formal_charge_raw: i32 = py_method_extract(&atom, "GetFormalCharge")?;
            let formal_charge = i8::try_from(formal_charge_raw).map_err(|_| {
                PyValueError::new_err(format!(
                    "from_rdkit atom {idx} formal charge out of i8 range: {formal_charge_raw}"
                ))
            })?;
            let explicit_h_raw: usize = py_method_extract(&atom, "GetNumExplicitHs")?;
            let explicit_hydrogens = u8::try_from(explicit_h_raw).map_err(|_| {
                PyValueError::new_err(format!(
                    "from_rdkit atom {idx} explicit H count out of u8 range: {explicit_h_raw}"
                ))
            })?;
            let radical_raw: usize = py_method_extract(&atom, "GetNumRadicalElectrons")?;
            let num_radical_electrons = u8::try_from(radical_raw).map_err(|_| {
                PyValueError::new_err(format!(
                    "from_rdkit atom {idx} radical electron count out of u8 range: {radical_raw}"
                ))
            })?;
            let atom_map_raw: u32 = py_method_extract(&atom, "GetAtomMapNum")?;
            let isotope_raw: u16 = py_method_extract(&atom, "GetIsotope")?;
            let chiral_tag = rdkit_chiral_tag_from_name(&py_method_str(&atom, "GetChiralTag")?)?;
            let mut spec = cosmolkit_core::AtomSpec::new(
                cosmolkit_core::Element::from_atomic_number(atomic_num).ok_or_else(|| {
                    PyValueError::new_err(format!(
                        "from_rdkit atom {idx} unsupported atomic number {atomic_num}"
                    ))
                })?,
            )
            .with_formal_charge(formal_charge)
            .with_explicit_hydrogens(explicit_hydrogens)
            .with_no_implicit(py_method_extract(&atom, "GetNoImplicit")?)
            .with_chiral_tag(chiral_tag)
            .with_aromatic(py_method_extract(&atom, "GetIsAromatic")?)
            .with_radical_electrons(num_radical_electrons);
            if isotope_raw != 0 {
                spec = spec.with_isotope(isotope_raw);
            }
            if atom_map_raw != 0 {
                spec = spec.with_atom_map(atom_map_raw);
            }
            builder.add_atom(spec);
        }

        for idx in 0..bond_count {
            let bond = py_method_index(rdmol, "GetBondWithIdx", idx)?;
            let begin_atom: usize = py_method_extract(&bond, "GetBeginAtomIdx")?;
            let end_atom: usize = py_method_extract(&bond, "GetEndAtomIdx")?;
            if begin_atom >= atom_count || end_atom >= atom_count {
                return Err(PyValueError::new_err(format!(
                    "from_rdkit bond {idx} atom index out of range: {begin_atom}-{end_atom}"
                )));
            }
            let is_aromatic: bool = py_method_extract(&bond, "GetIsAromatic")?;
            let order = rdkit_bond_order_from_name(&py_method_str(&bond, "GetBondType")?)?;
            let direction = rdkit_bond_direction_from_name(&py_method_str(&bond, "GetBondDir")?)?;
            let stereo = rdkit_bond_stereo_from_name(&py_method_str(&bond, "GetStereo")?)?;
            let stereo_atoms: Vec<usize> =
                py_method(&bond, "GetStereoAtoms")?
                    .extract()
                    .map_err(|err| {
                        PyValueError::new_err(format!(
                            "from_rdkit failed extracting result from GetStereoAtoms: {err}"
                        ))
                    })?;
            let begin = cosmolkit_core::AtomId::new(begin_atom);
            let end = cosmolkit_core::AtomId::new(end_atom);
            let mut spec = cosmolkit_core::BondSpec::new(begin, end, order)
                .with_aromatic(is_aromatic)
                .with_direction(direction)
                .with_stereo(stereo);
            if stereo_atoms.len() == 2 {
                spec = spec.with_stereo_atoms(
                    cosmolkit_core::AtomId::new(stereo_atoms[0]),
                    cosmolkit_core::AtomId::new(stereo_atoms[1]),
                );
            }
            builder.add_bond(spec).map_err(|err| {
                PyValueError::new_err(format!("from_rdkit failed adding bond {idx}: {err}"))
            })?;
        }

        let conformers = py_method(rdmol, "GetConformers")?;
        let conformer_iter = PyIterator::from_object(&conformers)?;
        for conformer in conformer_iter {
            let conformer = conformer.map_err(|err| {
                PyValueError::new_err(format!(
                    "from_rdkit failed iterating result from GetConformers: {err}"
                ))
            })?;
            if !py_method_extract::<bool>(&conformer, "Is3D")? {
                continue;
            }
            let mut coords = Vec::with_capacity(atom_count);
            for atom_idx in 0..atom_count {
                let pos = py_method_index(&conformer, "GetAtomPosition", atom_idx)?;
                coords.push([
                    py_attr_f64(&pos, "x")?,
                    py_attr_f64(&pos, "y")?,
                    py_attr_f64(&pos, "z")?,
                ]);
            }
            builder.add_3d_conformer(coords).map_err(|err| {
                PyValueError::new_err(format!("from_rdkit failed adding conformer: {err}"))
            })?;
        }
        let mol = builder
            .build()
            .map_err(|err| PyValueError::new_err(format!("from_rdkit build failed: {err}")))?;
        let inner = if matches!(sanitize, Some(true)) {
            mol.sanitize()
                .map_err(|err| PyValueError::new_err(err.to_string()))?
        } else {
            mol
        };
        Ok(Self { inner })
    }

    #[classmethod]
    #[pyo3(signature = (path, sanitize=None, coordinate_dim="auto"))]
    #[doc = r#"
Read the first molecule record from an SDF file.

This uses the SDF reader, so SDF data fields after the molfile ``M  END`` line
are parsed as record metadata. Use ``read_mol()`` for RDKit
``MolFromMolBlock``-style molfile-only parsing.

Parameters
----------
path : str
    SDF file path.
sanitize : bool, optional
    Optional molecule preparation flag.
coordinate_dim : {"auto", "2d", "3d"}, optional
    Coordinate interpretation mode. ``"auto"`` preserves the molfile header.
"#]
    fn read_sdf(
        _cls: &Bound<'_, PyType>,
        path: &str,
        sanitize: Option<bool>,
        coordinate_dim: &str,
    ) -> PyResult<Self> {
        reject_unsanitized_mol_reader(sanitize)?;
        let expanded_path = expand_user_path(path)?;
        let coordinate_mode = parse_coordinate_mode(Some(coordinate_dim))?;
        let file = File::open(&expanded_path)
            .map_err(|e| PyValueError::new_err(format!("read_sdf open failed: {e}")))?;
        let mut reader = SdfReader::with_coordinate_mode(BufReader::new(file), coordinate_mode);
        let Some(record) = reader
            .next_record()
            .map_err(|e| PyValueError::new_err(format!("read_sdf parse failed: {e:?}")))?
        else {
            return Err(PyValueError::new_err("read_sdf found no molecule record"));
        };
        Ok(Self {
            inner: record.molecule,
        })
    }

    #[classmethod]
    #[pyo3(signature = (path, sanitize=None, coordinate_dim="auto"))]
    #[doc = r#"
Read one molecule from an MDL molfile.

The parser follows RDKit ``MolFromMolBlock`` record boundaries: it reads the
molfile CTAB through the first ``M  END`` line and ignores unread trailing text,
including SDF data fields and ``$$$$`` record separators. Use ``read_sdf()`` or
``SdfDataset`` when SDF data fields must be parsed.

Parameters
----------
path : str
    Molfile path.
sanitize : bool, optional
    Optional molecule preparation flag.
coordinate_dim : {"auto", "2d", "3d"}, optional
    Coordinate interpretation mode. ``"auto"`` preserves the molfile header.
"#]
    fn read_mol(
        _cls: &Bound<'_, PyType>,
        path: &str,
        sanitize: Option<bool>,
        coordinate_dim: &str,
    ) -> PyResult<Self> {
        reject_unsanitized_mol_reader(sanitize)?;
        let expanded_path = expand_user_path(path)?;
        let coordinate_mode = parse_coordinate_mode(Some(coordinate_dim))?;
        let record = cosmolkit_core::io::molfile::read_mol_file_with_params(
            &expanded_path,
            cosmolkit_core::SdfReadParams {
                coordinate_mode,
                ..Default::default()
            },
        )
        .map_err(|e| PyValueError::new_err(format!("read_mol failed: {e}")))?;
        Ok(Self {
            inner: record.molecule,
        })
    }

    #[classmethod]
    #[pyo3(signature = (mol_text, sanitize=None, coordinate_dim="auto"))]
    #[doc = r#"
Read one molecule from an MDL molfile string.

The parser follows RDKit ``MolFromMolBlock`` record boundaries: it reads the
molfile CTAB through the first ``M  END`` line and ignores unread trailing text,
including SDF data fields and ``$$$$`` record separators. Use
``read_sdf_from_str()`` when SDF data fields must be parsed.
"#]
    fn read_mol_from_str(
        _cls: &Bound<'_, PyType>,
        mol_text: &str,
        sanitize: Option<bool>,
        coordinate_dim: &str,
    ) -> PyResult<Self> {
        reject_unsanitized_mol_reader(sanitize)?;
        let coordinate_mode = parse_coordinate_mode(Some(coordinate_dim))?;
        let record = cosmolkit_core::io::molfile::read_mol_record_from_str_with_params(
            mol_text,
            cosmolkit_core::SdfReadParams {
                coordinate_mode,
                ..Default::default()
            },
        )
        .map_err(|e| PyValueError::new_err(format!("read_mol_from_str failed: {e:?}")))?;
        Ok(Self {
            inner: record.molecule,
        })
    }

    #[classmethod]
    #[pyo3(signature = (sdf_text, sanitize=None, coordinate_dim="auto"))]
    #[doc = r#"
Read one molecule from an SDF record string.

This uses the SDF reader, so data fields after the molfile ``M  END`` line are
parsed as SDF record metadata. Use ``read_mol_from_str()`` for RDKit
``MolFromMolBlock``-style molfile-only parsing that ignores trailing SDF text.
"#]
    fn read_sdf_from_str(
        _cls: &Bound<'_, PyType>,
        sdf_text: &str,
        sanitize: Option<bool>,
        coordinate_dim: &str,
    ) -> PyResult<Self> {
        reject_unsanitized_mol_reader(sanitize)?;
        let coordinate_mode = parse_coordinate_mode(Some(coordinate_dim))?;
        let record = cosmolkit_core::io::sdf::read_sdf_from_str_with_coordinate_mode(
            sdf_text,
            coordinate_mode,
        )
        .map_err(|e| PyValueError::new_err(format!("read_sdf_from_str failed: {e:?}")))?;
        Ok(Self {
            inner: record.molecule,
        })
    }

    #[classmethod]
    #[pyo3(signature = (path, *, sanitize=true, remove_hs=true, variant="corina", cleanup_substructures=true))]
    #[doc = r#"
Read one molecule from a Tripos MOL2 file.

The reader follows the source-ported RDKit ``Mol2FileToMol``/``MolFromMol2File``
profile. The exposed parameters map to RDKit ``Mol2ParserParams``:
``sanitize``, ``removeHs``, ``variant``, and ``cleanupSubstructures``. The only
currently supported variant is ``"corina"``, matching RDKit's public enum.

Parameters
----------
path : str
    MOL2 file path.
sanitize : bool, optional
    Run RDKit-style MOL2 sanitization after parsing.
remove_hs : bool, optional
    Remove explicit hydrogens during MOL2 finalization.
variant : {"corina"}, optional
    MOL2 atom-type definition profile.
cleanup_substructures : bool, optional
    Run RDKit-style cleanup of common MOL2 substructures before charge
    assignment when formal charges are not present.
"#]
    fn read_mol2(
        _cls: &Bound<'_, PyType>,
        path: &str,
        sanitize: bool,
        remove_hs: bool,
        variant: &str,
        cleanup_substructures: bool,
    ) -> PyResult<Self> {
        let expanded_path = expand_user_path(path)?;
        let params = cosmolkit_core::Mol2ReadParams {
            sanitize,
            remove_hs,
            variant: parse_mol2_variant(variant)?,
            cleanup_substructures,
        };
        let Some(record) = cosmolkit_core::read_mol2_file_with_params(&expanded_path, params)
            .map_err(|e| PyValueError::new_err(format!("read_mol2 failed: {e}")))?
        else {
            return Err(PyValueError::new_err("read_mol2 found no molecule record"));
        };
        Ok(Self {
            inner: record.molecule,
        })
    }

    #[classmethod]
    #[pyo3(signature = (mol2_text, *, sanitize=true, remove_hs=true, variant="corina", cleanup_substructures=true))]
    #[doc = r#"
Read one molecule from a Tripos MOL2 string.

The reader follows the source-ported RDKit ``Mol2BlockToMol``/``MolFromMol2Block``
profile. The exposed parameters map to RDKit ``Mol2ParserParams``:
``sanitize``, ``removeHs``, ``variant``, and ``cleanupSubstructures``. The only
currently supported variant is ``"corina"``, matching RDKit's public enum.
"#]
    fn read_mol2_from_str(
        _cls: &Bound<'_, PyType>,
        mol2_text: &str,
        sanitize: bool,
        remove_hs: bool,
        variant: &str,
        cleanup_substructures: bool,
    ) -> PyResult<Self> {
        let params = cosmolkit_core::Mol2ReadParams {
            sanitize,
            remove_hs,
            variant: parse_mol2_variant(variant)?,
            cleanup_substructures,
        };
        let Some(record) = cosmolkit_core::read_mol2_from_str_with_params(mol2_text, params)
            .map_err(|e| PyValueError::new_err(format!("read_mol2_from_str failed: {e}")))?
        else {
            return Err(PyValueError::new_err(
                "read_mol2_from_str found no molecule record",
            ));
        };
        Ok(Self {
            inner: record.molecule,
        })
    }

    #[doc = r#"
Return a new molecule with explicit hydrogens added.

The original ``Molecule`` value is left unchanged.
"#]
    fn with_hydrogens(&self) -> PyResult<Self> {
        let out = self
            .inner
            .with_hydrogens()
            .map_err(|err| PyValueError::new_err(format!("with_hydrogens failed: {err:?}")))?;
        Ok(Self { inner: out })
    }

    #[doc = r#"
Add explicit hydrogens in place.

This is the in-place version of ``with_hydrogens()``.

All public in-place ``Molecule`` methods end with ``_``. If this method returns
an error, the receiver is not guaranteed to equal its pre-call value; use
``with_hydrogens()`` when failure-preserving value semantics are required.
"#]
    fn add_hydrogens_(&mut self) -> PyResult<()> {
        self.inner
            .add_hydrogens_()
            .map_err(|err| PyValueError::new_err(format!("add_hydrogens_ failed: {err:?}")))
    }

    #[pyo3(signature = (sanitize=None))]
    #[doc = r#"
Return a new molecule with explicit hydrogens removed.

The original ``Molecule`` value is left unchanged.
"#]
    fn without_hydrogens(&self, sanitize: Option<bool>) -> PyResult<Self> {
        let out = self
            .inner
            .without_hydrogens_with_sanitize(sanitize.unwrap_or(true))
            .map_err(|err| PyValueError::new_err(format!("without_hydrogens failed: {err:?}")))?;
        Ok(Self { inner: out })
    }

    #[pyo3(signature = (sanitize=None))]
    #[doc = r#"
Remove explicit hydrogens in place.

This is the in-place version of ``without_hydrogens()``.
"#]
    fn remove_hydrogens_(&mut self, sanitize: Option<bool>) -> PyResult<()> {
        self.inner
            .remove_hydrogens_with_sanitize_(sanitize.unwrap_or(true))
            .map_err(|err| PyValueError::new_err(format!("remove_hydrogens_ failed: {err:?}")))
    }

    #[pyo3(signature = (clear_aromatic_flags=None))]
    #[doc = r#"
Return a new molecule with aromatic bonds converted to an explicit Kekule form.

The original ``Molecule`` value is left unchanged.
"#]
    fn with_kekulized_bonds(&self, clear_aromatic_flags: Option<bool>) -> PyResult<Self> {
        let out = self
            .inner
            .with_kekulized_bonds(clear_aromatic_flags.unwrap_or(true))
            .map_err(|err| {
                PyValueError::new_err(format!("with_kekulized_bonds failed: {err:?}"))
            })?;
        Ok(Self { inner: out })
    }

    #[pyo3(signature = (clear_aromatic_flags=None))]
    #[doc = r#"
Convert aromatic bonds to an explicit Kekule form in place.

This is the in-place version of ``with_kekulized_bonds()``.
"#]
    fn kekulize_(&mut self, clear_aromatic_flags: Option<bool>) -> PyResult<()> {
        self.inner
            .kekulize_(clear_aromatic_flags.unwrap_or(true))
            .map_err(|err| PyValueError::new_err(format!("kekulize_ failed: {err:?}")))
    }

    #[doc = r#"
Return the number of atoms.
"#]
    fn num_atoms(&self) -> usize {
        self.inner.atoms().len()
    }

    #[doc = r#"
Return the number of bonds.
"#]
    fn num_bonds(&self) -> usize {
        self.inner.bonds().len()
    }

    #[doc = r#"
Return read-only atom feature records.
"#]
    fn atoms(&self) -> Vec<Atom> {
        let assignment =
            cosmolkit_core::assign_valence(&self.inner, cosmolkit_core::ValenceModel::RdkitLike)
                .ok();
        let mut degrees = vec![0usize; self.inner.atoms().len()];
        for bond in self.inner.bonds() {
            degrees[bond.begin().index()] += 1;
            degrees[bond.end().index()] += 1;
        }
        self.inner
            .atoms()
            .iter()
            .map(|atom| Atom {
                idx: atom.id().index(),
                atomic_num: atom.atomic_number() as usize,
                formal_charge: atom.formal_charge(),
                chiral_tag_name: atom.chiral_tag().rdkit_name().to_string(),
                chiral_tag_code: atom.chiral_tag().python_code(),
                isotope: atom.isotope(),
                atom_map_num: atom.atom_map(),
                is_aromatic: atom.is_aromatic(),
                explicit_hydrogens: atom.explicit_hydrogens() as usize,
                no_implicit: atom.no_implicit(),
                num_radical_electrons: atom.radical_electrons() as usize,
                degree: degrees[atom.id().index()],
                explicit_valence: assignment
                    .as_ref()
                    .map(|v| v.explicit_valence[atom.id().index()] as usize),
                implicit_hydrogens: assignment
                    .as_ref()
                    .map(|v| v.implicit_hydrogens[atom.id().index()] as usize),
                total_num_hs: assignment.as_ref().map(|v| {
                    atom.explicit_hydrogens() as usize
                        + v.implicit_hydrogens[atom.id().index()] as usize
                }),
                total_valence: assignment.as_ref().map(|v| {
                    v.explicit_valence[atom.id().index()] as usize
                        + v.implicit_hydrogens[atom.id().index()] as usize
                }),
            })
            .collect()
    }

    #[doc = r#"
Return read-only bond feature records.
"#]
    fn bonds(&self) -> Vec<Bond> {
        self.inner
            .bonds()
            .iter()
            .map(|bond| Bond {
                idx: bond.id().index(),
                begin_atom_idx: bond.begin().index(),
                end_atom_idx: bond.end().index(),
                bond_type_name: bond.order().rdkit_name().to_string(),
                bond_type_code: bond.order().python_code(),
                bond_dir_name: bond.direction().rdkit_name().to_string(),
                bond_dir_code: bond.direction().python_code(),
                stereo_name: bond.stereo().rdkit_name().to_string(),
                stereo_code: bond.stereo().python_code(),
                stereo_atoms: bond
                    .stereo_atoms()
                    .map(|refs| refs.map(|id| id.index()).to_vec())
                    .unwrap_or_default(),
                is_aromatic: bond.is_aromatic(),
            })
            .collect()
    }

    #[pyo3(signature = (include_unassigned=true))]
    #[doc = r#"
Return chiral center labels.

Parameters
----------
include_unassigned : bool, default True
    Include atoms with unspecified tetrahedral chirality.
"#]
    fn find_chiral_centers(&self, include_unassigned: bool) -> Vec<(usize, String)> {
        self.inner
            .atoms()
            .iter()
            .filter_map(|atom| match atom.chiral_tag() {
                cosmolkit_core::ChiralTag::Unspecified => {
                    if include_unassigned {
                        Some((atom.id().index(), "?".to_string()))
                    } else {
                        None
                    }
                }
                cosmolkit_core::ChiralTag::TetrahedralCw => {
                    Some((atom.id().index(), "CHI_TETRAHEDRAL_CW".to_string()))
                }
                cosmolkit_core::ChiralTag::TetrahedralCcw => {
                    Some((atom.id().index(), "CHI_TETRAHEDRAL_CCW".to_string()))
                }
                cosmolkit_core::ChiralTag::TrigonalBipyramidal => {
                    Some((atom.id().index(), "CHI_TRIGONALBIPYRAMIDAL".to_string()))
                }
                _ => None,
            })
            .collect()
    }

    #[doc = r#"
Return ordered tetrahedral stereo ligand records.

Each record is ``(center_atom_index, ordered_ligands)``. Implicit hydrogen is
represented as ``None``.
"#]
    fn tetrahedral_stereo(&self) -> PyResult<Vec<(usize, Vec<Option<usize>>)>> {
        to_python_tetrahedral_stereo(&self.inner)
    }

    #[doc = r#"
Return a new molecule with 2D coordinates.
"#]
    fn with_2d_coordinates(&self) -> PyResult<Self> {
        let out = self
            .inner
            .with_2d_coordinates()
            .map_err(|err| PyValueError::new_err(format!("with_2d_coordinates failed: {err}")))?;
        Ok(Self { inner: out })
    }

    #[doc = r#"
Compute 2D coordinates in place.

This is the in-place version of ``with_2d_coordinates()``.
"#]
    fn compute_2d_coordinates_(&mut self) -> PyResult<()> {
        self.inner
            .compute_2d_coordinates_()
            .map_err(|err| PyValueError::new_err(format!("compute_2d_coordinates_ failed: {err}")))
    }

    #[pyo3(signature = (params=None))]
    #[doc = r#"
Return a new molecule with one generated 3D conformer.

Parameters
----------
params : EmbedParameters, optional
    Distance-geometry embedding parameters. The default is ``EmbedParameters.etkdg_v3()``.

Returns
-------
Molecule
    A new molecule value containing one additional 3D conformer.
"#]
    fn with_3d_conformer(&self, params: Option<PyRefMut<'_, PyEmbedParameters>>) -> PyResult<Self> {
        Ok(self.with_3d_conformer_result(params)?.molecule)
    }

    #[pyo3(signature = (params=None))]
    #[doc = r#"
Generate one 3D conformer in place.

This is the in-place version of ``with_3d_conformer()``.
"#]
    fn embed_3d_conformer_(
        &mut self,
        params: Option<PyRefMut<'_, PyEmbedParameters>>,
    ) -> PyResult<()> {
        self.inner = self.embed_3d_conformer_result_(params)?.molecule.inner;
        Ok(())
    }

    #[pyo3(signature = (params=None))]
    #[doc = r#"
Return an embedding result object for one generated 3D conformer.

The result keeps the embedded molecule, the returned conformer id, and the
final parameter snapshot so callers can inspect status and failure counters
without relying on side effects on the input ``EmbedParameters`` object.
"#]
    fn with_3d_conformer_result(
        &self,
        mut params: Option<PyRefMut<'_, PyEmbedParameters>>,
    ) -> PyResult<EmbedMoleculeResult> {
        let mut embed_params = params
            .as_ref()
            .map(|value| value.inner.clone())
            .unwrap_or_else(cosmolkit_core::EmbedParameters::etkdg_v3);
        let result = cosmolkit_core::embed_molecule_result(&self.inner, &mut embed_params)
            .map_err(distgeom_pyerr)?;
        if let Some(value) = params.as_mut() {
            value.inner = embed_params.clone();
        }
        Ok(result.into())
    }

    #[pyo3(signature = (params=None))]
    #[doc = r#"
Generate one 3D conformer in place and return the embedding result object.

This is the in-place version of ``with_3d_conformer_result()``.
"#]
    fn embed_3d_conformer_result_(
        &mut self,
        params: Option<PyRefMut<'_, PyEmbedParameters>>,
    ) -> PyResult<EmbedMoleculeResult> {
        let result = self.with_3d_conformer_result(params)?;
        self.inner = result.molecule.inner.clone();
        Ok(result)
    }

    #[pyo3(signature = (num_confs, params=None))]
    #[doc = r#"
Return a new molecule with multiple generated 3D conformers.

Parameters
----------
num_confs : int
    Number of conformers to request.
params : EmbedParameters, optional
    Distance-geometry embedding parameters.

Returns
-------
Molecule
    A new molecule value containing the generated 3D conformers.
"#]
    fn with_3d_conformers(
        &self,
        num_confs: u32,
        params: Option<PyRefMut<'_, PyEmbedParameters>>,
    ) -> PyResult<Self> {
        Ok(self.with_3d_conformers_result(num_confs, params)?.molecule)
    }

    #[pyo3(signature = (num_confs, params=None))]
    #[doc = r#"
Generate multiple 3D conformers in place.

This is the in-place version of ``with_3d_conformers()``.
"#]
    fn embed_3d_conformers_(
        &mut self,
        num_confs: u32,
        params: Option<PyRefMut<'_, PyEmbedParameters>>,
    ) -> PyResult<()> {
        self.inner = self
            .embed_3d_conformers_result_(num_confs, params)?
            .molecule
            .inner;
        Ok(())
    }

    #[pyo3(signature = (num_confs, params=None))]
    #[doc = r#"
Return an embedding result object for multiple generated 3D conformers.

The result keeps the embedded molecule, the kept conformer ids, and the final
parameter snapshot so callers can inspect pruning and tracked failures without
reconstructing that state manually.
"#]
    fn with_3d_conformers_result(
        &self,
        num_confs: u32,
        mut params: Option<PyRefMut<'_, PyEmbedParameters>>,
    ) -> PyResult<EmbedMultipleConfsResult> {
        let mut embed_params = params
            .as_ref()
            .map(|value| value.inner.clone())
            .unwrap_or_else(cosmolkit_core::EmbedParameters::etkdg_v3);
        let result =
            cosmolkit_core::embed_multiple_confs_result(&self.inner, num_confs, &mut embed_params)
                .map_err(distgeom_pyerr)?;
        if let Some(value) = params.as_mut() {
            value.inner = embed_params.clone();
        }
        Ok(result.into())
    }

    #[pyo3(signature = (num_confs, params=None))]
    #[doc = r#"
Generate multiple 3D conformers in place and return the embedding result object.

This is the in-place version of ``with_3d_conformers_result()``.
"#]
    fn embed_3d_conformers_result_(
        &mut self,
        num_confs: u32,
        params: Option<PyRefMut<'_, PyEmbedParameters>>,
    ) -> PyResult<EmbedMultipleConfsResult> {
        let result = self.with_3d_conformers_result(num_confs, params)?;
        self.inner = result.molecule.inner.clone();
        Ok(result)
    }

    #[doc = r#"
Return the number of stored 3D conformers.
"#]
    fn num_conformers(&self) -> usize {
        self.inner.conformers_3d().len()
    }

    #[doc = r#"
Return whether the molecule has 2D coordinates.
"#]
    fn has_2d_coordinates(&self) -> bool {
        self.inner.coordinates_2d().is_some()
    }

    #[gen_stub(override_return_type(type_repr = "numpy.ndarray[typing.Any, typing.Any]", imports = ("numpy", "typing")))]
    #[doc = r#"
Return 2D coordinates as a NumPy array with shape ``(num_atoms, 3)``.

The z column is zero-filled.
"#]
    fn coordinates_2d<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let Some(coords) = self.inner.coordinates_2d() else {
            return Err(PyValueError::new_err(
                "no 2D coordinates present; call with_2d_coordinates() first",
            ));
        };
        let rows: Vec<Vec<f64>> = coords.iter().map(|p| vec![p[0], p[1], 0.0]).collect();
        PyArray2::from_vec2(py, &rows)
            .map(|array| array.into_any())
            .map_err(|err| PyValueError::new_err(format!("Molecule.coordinates_2d failed: {err}")))
    }

    #[pyo3(signature = (conformer_index=0))]
    #[gen_stub(override_return_type(type_repr = "numpy.ndarray[typing.Any, typing.Any]", imports = ("numpy", "typing")))]
    #[doc = r#"
Return 3D coordinates as a NumPy array with shape ``(num_atoms, 3)``.
"#]
    fn coordinates_3d<'py>(
        &self,
        py: Python<'py>,
        conformer_index: usize,
    ) -> PyResult<Bound<'py, PyAny>> {
        let Some(coords) = self.inner.conformers_3d().get(conformer_index) else {
            return Err(PyValueError::new_err(format!(
                "no 3D conformer present at index {conformer_index}"
            )));
        };
        let rows: Vec<Vec<f64>> = coords
            .coordinates()
            .iter()
            .map(|p| vec![p[0], p[1], p[2]])
            .collect();
        PyArray2::from_vec2(py, &rows)
            .map(|array| array.into_any())
            .map_err(|err| PyValueError::new_err(format!("Molecule.coordinates_3d failed: {err}")))
    }

    #[doc = r#"
Return whether UFF parameters are available for every atom in this molecule.
"#]
    fn has_uff_params(&self) -> PyResult<bool> {
        cosmolkit_core::uff_has_all_molecule_params(&self.inner).map_err(forcefield_pyerr)
    }

    #[pyo3(signature = (max_iters=1000, vdw_thresh=10.0, conf_id=-1, ignore_interfrag_interactions=true))]
    #[doc = r#"
Return a UFF optimization result with a new optimized molecule value.

The source molecule is not mutated. The molecule must already contain a 3D
conformer, for example from a 3D SDF, MOL, MOL2, or XYZ input.
"#]
    fn with_uff_optimized(
        &self,
        max_iters: usize,
        vdw_thresh: f64,
        conf_id: isize,
        ignore_interfrag_interactions: bool,
    ) -> PyResult<UffOptimizeMoleculeResult> {
        uff_optimize_molecule(
            self,
            max_iters,
            vdw_thresh,
            conf_id,
            ignore_interfrag_interactions,
        )
    }

    #[pyo3(signature = (num_threads=1, max_iters=1000, vdw_thresh=10.0, ignore_interfrag_interactions=true))]
    #[doc = r#"
Return UFF optimization results for all 3D conformers as a new molecule value.
"#]
    fn with_uff_optimized_confs(
        &self,
        num_threads: i32,
        max_iters: usize,
        vdw_thresh: f64,
        ignore_interfrag_interactions: bool,
    ) -> PyResult<UffOptimizeMoleculeConfsResult> {
        uff_optimize_molecule_confs(
            self,
            num_threads,
            max_iters,
            vdw_thresh,
            ignore_interfrag_interactions,
        )
    }

    #[doc = r#"
Return whether MMFF94 parameters are available for this molecule.
"#]
    fn has_mmff_params(&self) -> PyResult<bool> {
        cosmolkit_core::mmff_has_all_molecule_params(&self.inner).map_err(forcefield_pyerr)
    }

    #[pyo3(signature = (mmff_variant="MMFF94", max_iters=200, non_bonded_thresh=100.0, conf_id=-1, ignore_interfrag_interactions=true))]
    #[doc = r#"
Return an MMFF optimization result with a new optimized molecule value.

The source molecule is not mutated. The molecule must already contain a 3D
conformer. Supported variants follow the Rust core parser, including
``"MMFF94"`` and ``"MMFF94S"``.
"#]
    fn with_mmff_optimized(
        &self,
        mmff_variant: &str,
        max_iters: usize,
        non_bonded_thresh: f64,
        conf_id: isize,
        ignore_interfrag_interactions: bool,
    ) -> PyResult<MmffOptimizeMoleculeResult> {
        mmff_optimize_molecule(
            self,
            mmff_variant,
            max_iters,
            non_bonded_thresh,
            conf_id,
            ignore_interfrag_interactions,
        )
    }

    #[pyo3(signature = (num_threads=1, max_iters=1000, mmff_variant="MMFF94", non_bonded_thresh=10.0, ignore_interfrag_interactions=true))]
    #[doc = r#"
Return MMFF optimization results for all 3D conformers as a new molecule value.
"#]
    fn with_mmff_optimized_confs(
        &self,
        num_threads: i32,
        max_iters: usize,
        mmff_variant: &str,
        non_bonded_thresh: f64,
        ignore_interfrag_interactions: bool,
    ) -> PyResult<MmffOptimizeMoleculeConfsResult> {
        mmff_optimize_molecule_confs(
            self,
            num_threads,
            max_iters,
            mmff_variant,
            non_bonded_thresh,
            ignore_interfrag_interactions,
        )
    }

    #[gen_stub(override_return_type(type_repr = "numpy.ndarray[typing.Any, typing.Any]", imports = ("numpy", "typing")))]
    #[doc = r#"
Return the distance-geometry bounds matrix as a NumPy array.

The returned array uses shape ``(num_atoms, num_atoms)``.
"#]
    fn dg_bounds_matrix<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let rows = self.inner.dg_bounds_matrix().map_err(|err| {
            PyValueError::new_err(format!("Molecule.dg_bounds_matrix failed: {err}"))
        })?;
        PyArray2::from_vec2(py, &rows)
            .map(|array| array.into_any())
            .map_err(|err| {
                PyValueError::new_err(format!("Molecule.dg_bounds_matrix failed: {err}"))
            })
    }

    #[pyo3(signature = (width=300, height=300))]
    #[doc = r#"
Render the molecule to an SVG string.
"#]
    fn to_svg(&self, width: u32, height: u32) -> PyResult<String> {
        self.inner
            .to_svg(width, height)
            .map_err(|err| PyNotImplementedError::new_err(format!("Molecule.to_svg failed: {err}")))
    }

    #[pyo3(signature = (path, width=300, height=300))]
    #[doc = r#"
Write an SVG depiction to a file.
"#]
    fn write_svg(&self, path: &str, width: u32, height: u32) -> PyResult<()> {
        let expanded_path = expand_user_path(path)?;
        let svg = self
            .inner
            .to_svg(width, height)
            .map_err(|err| PyNotImplementedError::new_err(format!("write_svg failed: {err}")))?;
        let mut f = File::create(&expanded_path)
            .map_err(|e| PyValueError::new_err(format!("write_svg create failed: {e}")))?;
        f.write_all(svg.as_bytes())
            .map_err(|e| PyValueError::new_err(format!("write_svg write failed: {e}")))?;
        Ok(())
    }

    #[pyo3(signature = (path, width=300, height=300))]
    #[doc = r#"
Write a PNG depiction to a file.
"#]
    fn write_png(&self, path: &str, width: u32, height: u32) -> PyResult<()> {
        let expanded_path = expand_user_path(path)?;
        let png = self
            .inner
            .to_png(width, height)
            .map_err(|err| PyNotImplementedError::new_err(format!("write_png failed: {err}")))?;
        let mut f = File::create(&expanded_path)
            .map_err(|e| PyValueError::new_err(format!("write_png create failed: {e}")))?;
        f.write_all(&png)
            .map_err(|e| PyValueError::new_err(format!("write_png write failed: {e}")))?;
        Ok(())
    }

    #[pyo3(signature = (width=300, height=300))]
    #[doc = r#"
Render the molecule to PNG bytes.
"#]
    fn to_png<'py>(
        &self,
        py: Python<'py>,
        width: u32,
        height: u32,
    ) -> PyResult<Bound<'py, PyBytes>> {
        let png = self.inner.to_png(width, height).map_err(svg_draw_pyerr)?;
        Ok(PyBytes::new(py, &png))
    }

    #[pyo3(signature = (
        isomeric_smiles=true,
        canonical=true,
        kekule=false,
        clean_stereo=true,
        all_bonds_explicit=false,
        all_hs_explicit=false,
        include_dative_bonds=true,
        ignore_atom_map_numbers=false,
        rooted_at_atom=None
    ))]
    #[doc = r#"
Return a SMILES string.

Parameters
----------
isomeric_smiles : bool, default True
    Include stereochemical and isotopic information when available.
canonical : bool, default True
    Return a canonical SMILES when supported.
kekule : bool, default False
    Write aromatic systems using Kekule bond notation.
clean_stereo : bool, default True
    Normalize stereo annotations before writing.
all_bonds_explicit : bool, default False
    Write explicit bond symbols.
all_hs_explicit : bool, default False
    Write explicit hydrogens.
include_dative_bonds : bool, default True
    Include dative bond notation.
ignore_atom_map_numbers : bool, default False
    Omit atom map numbers from canonical decisions.
rooted_at_atom : int, optional
    Start traversal from a selected atom index.
"#]
    #[allow(clippy::too_many_arguments)]
    fn to_smiles(
        &self,
        isomeric_smiles: bool,
        canonical: bool,
        kekule: bool,
        clean_stereo: bool,
        all_bonds_explicit: bool,
        all_hs_explicit: bool,
        include_dative_bonds: bool,
        ignore_atom_map_numbers: bool,
        rooted_at_atom: Option<usize>,
    ) -> PyResult<String> {
        let params = make_smiles_write_params(
            isomeric_smiles,
            canonical,
            kekule,
            clean_stereo,
            all_bonds_explicit,
            all_hs_explicit,
            include_dative_bonds,
            ignore_atom_map_numbers,
            rooted_at_atom,
        );
        self.inner.to_smiles_with_params(&params).map_err(|err| {
            PyNotImplementedError::new_err(format!(
                "Molecule.to_smiles(...) is not implemented for these parameters yet: {err}"
            ))
        })
    }

    #[pyo3(signature = (path, format=None, include_stereo=true, kekulize=true))]
    #[doc = r#"
Write the molecule as one SDF record.
"#]
    fn write_sdf(
        &self,
        path: &str,
        format: Option<&str>,
        include_stereo: bool,
        kekulize: bool,
    ) -> PyResult<()> {
        let expanded_path = expand_user_path(path)?;
        let fmt = parse_sdf_format(format)?;
        let block = molecule_to_sdf_record_string(&self.inner, fmt, include_stereo, kekulize)
            .map_err(|err| PyValueError::new_err(format!("write_sdf failed: {err}")))?;
        let mut f = File::create(&expanded_path)
            .map_err(|e| PyValueError::new_err(format!("write_sdf create failed: {e}")))?;
        f.write_all(block.as_bytes())
            .map_err(|e| PyValueError::new_err(format!("write_sdf write failed: {e}")))?;
        Ok(())
    }

    #[pyo3(signature = (format=None, include_stereo=true, kekulize=true))]
    #[doc = r#"
Return the molecule as a 2D SDF record string.

If the molecule does not already have 2D coordinates, they are generated for
this export. The original ``Molecule`` value is left unchanged.
"#]
    fn to_2d_sdf_string(
        &self,
        format: Option<&str>,
        include_stereo: bool,
        kekulize: bool,
    ) -> PyResult<String> {
        let fmt = parse_sdf_format(format)?;
        molecule_to_2d_sdf_record_string(&self.inner, fmt, include_stereo, kekulize)
            .map_err(|err| PyValueError::new_err(format!("to_2d_sdf_string failed: {err}")))
    }

    #[pyo3(signature = (format=None, include_stereo=true, kekulize=true))]
    #[doc = r#"
Return the molecule as a 3D SDF record string.

The molecule must already have a 3D conformer, for example from a 3D SDF record.
"#]
    fn to_3d_sdf_string(
        &self,
        format: Option<&str>,
        include_stereo: bool,
        kekulize: bool,
    ) -> PyResult<String> {
        let fmt = parse_sdf_format(format)?;
        molecule_to_3d_sdf_record_string(&self.inner, fmt, include_stereo, kekulize)
            .map_err(|err| PyValueError::new_err(format!("to_3d_sdf_string failed: {err}")))
    }

    #[pyo3(signature = (directory, file_name=None, format=None, include_stereo=true, kekulize=true))]
    #[doc = r#"
Write the molecule as one SDF record inside a directory.

Returns
-------
str
    The output path.
"#]
    fn write_sdf_to_directory(
        &self,
        directory: &str,
        file_name: Option<&str>,
        format: Option<&str>,
        include_stereo: bool,
        kekulize: bool,
    ) -> PyResult<String> {
        let expanded_directory = expand_user_path(directory)?;
        let dir = expanded_directory.as_path();
        if !dir.exists() {
            return Err(PyValueError::new_err(format!(
                "directory does not exist: {directory}"
            )));
        }
        if !dir.is_dir() {
            return Err(PyValueError::new_err(format!(
                "path is not a directory: {directory}"
            )));
        }
        let name = file_name.unwrap_or("molecule.sdf");
        if name.trim().is_empty() {
            return Err(PyValueError::new_err("file_name cannot be empty"));
        }
        let output = dir.join(name);
        let output_str = output
            .to_str()
            .ok_or_else(|| PyValueError::new_err("output path is not valid UTF-8"))?;
        self.write_sdf(output_str, format, include_stereo, kekulize)?;
        Ok(output_str.to_string())
    }

    #[doc = r#"
Create an explicit edit context for this molecule.

The edit context is useful when several changes should be staged and committed
as one new molecule value.
"#]
    fn edit(&self) -> MoleculeEdit {
        MoleculeEdit {
            builder: self.inner.to_builder(),
        }
    }

    #[pyo3(signature = (
        min_path=1,
        max_path=7,
        n_bits=2048,
        n_bits_per_hash=1,
        use_bond_order=true,
        use_hs=false,
        tautomeric_fingerprint=false,
        from_atoms=None
    ))]
    #[doc = r#"
Return an Avalon fingerprint.
"#]
    fn avalon_fingerprint(
        &self,
        min_path: u32,
        max_path: u32,
        n_bits: usize,
        n_bits_per_hash: u32,
        use_bond_order: bool,
        use_hs: bool,
        tautomeric_fingerprint: bool,
        from_atoms: Option<Vec<usize>>,
    ) -> PyResult<Fingerprint> {
        let params = make_avalon_fingerprint_params(
            min_path,
            max_path,
            n_bits,
            n_bits_per_hash,
            use_bond_order,
            use_hs,
            tautomeric_fingerprint,
            from_atoms,
        );
        self.inner
            .avalon_fingerprint(&params)
            .map(|inner| Fingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        radius=2,
        n_bits=2048,
        include_chirality=false,
        use_bond_types=true,
        count_simulation=false,
        count_bounds=None,
        only_nonzero_invariants=false,
        include_redundant_environments=false,
        from_atoms=None,
        ignore_atoms=None,
        custom_atom_invariants=None,
        custom_bond_invariants=None,
        atom_invariants_generator=None,
        atom_invariants_include_ring_membership=true,
        bond_invariants_generator=None,
        bond_invariants_use_bond_types=true,
        bond_invariants_use_chirality=false,
        num_bits_per_feature=1
    ))]
    #[doc = r#"
Return an RDKit-style Morgan fingerprint.

Parameters
----------
radius : int, default 2
    Morgan neighborhood radius.
n_bits : int, default 2048
    Output bit vector size.
include_chirality : bool, default False
    Include atom chirality in invariant updates.
use_bond_types : bool, default True
    Include bond order in invariant updates.
count_simulation : bool, default False
    Apply RDKit count-simulation bit expansion.
count_bounds : list[int], optional
    Count-simulation thresholds. Defaults to ``[1, 2, 4, 8]``.
only_nonzero_invariants : bool, default False
    Skip atoms whose starting invariant is zero.
include_redundant_environments : bool, default False
    Retain duplicate environments instead of applying RDKit redundancy checks.
from_atoms : list[int], optional
    Restrict environments to these root atoms.
ignore_atoms : list[int], optional
    Accepted for RDKit API parity; Morgan currently ignores this input.
custom_atom_invariants : list[int], optional
    Per-atom starting invariants.
custom_bond_invariants : list[int], optional
    Per-bond invariants.
atom_invariants_generator : {"connectivity", "morgan", "feature", "fcfp"}, optional
    Explicit atom invariant generator. ``None`` uses the Morgan connectivity default.
atom_invariants_include_ring_membership : bool, default True
    Include ring membership for the connectivity invariant generator.
bond_invariants_generator : {"morgan", "default", "bond"}, optional
    Explicit Morgan bond invariant generator. ``None`` uses the fingerprint defaults.
bond_invariants_use_bond_types : bool, default True
    Include bond order in the explicit bond invariant generator.
bond_invariants_use_chirality : bool, default False
    Include bond stereo in the explicit bond invariant generator.
num_bits_per_feature : int, default 1
    Number of bits set for each feature.
"#]
    fn fingerprint_morgan(
        &self,
        radius: u32,
        n_bits: usize,
        include_chirality: bool,
        use_bond_types: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        only_nonzero_invariants: bool,
        include_redundant_environments: bool,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        atom_invariants_generator: Option<&str>,
        atom_invariants_include_ring_membership: bool,
        bond_invariants_generator: Option<&str>,
        bond_invariants_use_bond_types: bool,
        bond_invariants_use_chirality: bool,
        num_bits_per_feature: u32,
    ) -> PyResult<Fingerprint> {
        let params = make_morgan_fingerprint_params(
            radius,
            n_bits,
            include_chirality,
            use_bond_types,
            count_simulation,
            count_bounds,
            only_nonzero_invariants,
            include_redundant_environments,
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            atom_invariants_generator,
            atom_invariants_include_ring_membership,
            bond_invariants_generator,
            bond_invariants_use_bond_types,
            bond_invariants_use_chirality,
            num_bits_per_feature,
            false,
        )?;
        self.inner
            .morgan_fingerprint(&params)
            .map(|inner| Fingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        radius=2,
        n_bits=2048,
        include_chirality=false,
        use_bond_types=true,
        count_simulation=false,
        count_bounds=None,
        only_nonzero_invariants=false,
        include_redundant_environments=false,
        from_atoms=None,
        ignore_atoms=None,
        custom_atom_invariants=None,
        custom_bond_invariants=None,
        atom_invariants_generator=None,
        atom_invariants_include_ring_membership=true,
        bond_invariants_generator=None,
        bond_invariants_use_bond_types=true,
        bond_invariants_use_chirality=false,
        num_bits_per_feature=1
    ))]
    #[doc = r#"
Return a Morgan fingerprint together with allocated RDKit-style additional output.
"#]
    fn fingerprint_morgan_with_output(
        &self,
        radius: u32,
        n_bits: usize,
        include_chirality: bool,
        use_bond_types: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        only_nonzero_invariants: bool,
        include_redundant_environments: bool,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        atom_invariants_generator: Option<&str>,
        atom_invariants_include_ring_membership: bool,
        bond_invariants_generator: Option<&str>,
        bond_invariants_use_bond_types: bool,
        bond_invariants_use_chirality: bool,
        num_bits_per_feature: u32,
    ) -> PyResult<MorganFingerprintResult> {
        let params = make_morgan_fingerprint_params(
            radius,
            n_bits,
            include_chirality,
            use_bond_types,
            count_simulation,
            count_bounds,
            only_nonzero_invariants,
            include_redundant_environments,
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            atom_invariants_generator,
            atom_invariants_include_ring_membership,
            bond_invariants_generator,
            bond_invariants_use_bond_types,
            bond_invariants_use_chirality,
            num_bits_per_feature,
            true,
        )?;
        let output = self
            .inner
            .morgan_fingerprint_with_output(&params)
            .map_err(fingerprint_pyerr)?;
        Ok(MorganFingerprintResult {
            fingerprint: Fingerprint {
                inner: output.fingerprint,
            },
            additional_output: output.additional_output.map(Into::into),
        })
    }

    #[pyo3(signature = (
        min_path=1,
        max_path=7,
        n_bits=2048,
        n_bits_per_hash=2,
        use_bond_types=true,
        from_atoms=None,
        ignore_atoms=None
    ))]
    #[doc = r#"
Return a topological fingerprint.
"#]
    fn topological_fingerprint(
        &self,
        min_path: u32,
        max_path: u32,
        n_bits: usize,
        n_bits_per_hash: u32,
        use_bond_types: bool,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
    ) -> Fingerprint {
        let params = make_topological_fingerprint_params(
            min_path,
            max_path,
            n_bits,
            n_bits_per_hash,
            use_bond_types,
            from_atoms,
            ignore_atoms,
        );
        Fingerprint {
            inner: self.inner.topological_fingerprint(&params),
        }
    }

    #[pyo3(signature = (n_bits=166))]
    #[doc = r#"
Return a MACCS fingerprint.
"#]
    fn maccs_fingerprint(&self, n_bits: usize) -> Fingerprint {
        let params = cosmolkit_core::fingerprint::MaccsFingerprintParams { n_bits };
        Fingerprint {
            inner: self.inner.maccs_fingerprint(&params),
        }
    }

    #[doc = r#"
Return a hash of the molecule.
"#]
    fn hash(&self) -> PyResult<u64> {
        self.inner.hash().map_err(hash_pyerr)
    }

    #[doc = r#"
Return a hash of the molecule using the provided atom ranks.
"#]
    fn hash_with_ranks(&self, ranks: Vec<u32>) -> PyResult<u64> {
        self.inner.hash_with_ranks(&ranks).map_err(hash_pyerr)
    }

    #[doc = r#"
Return the connected fragments as separate molecules.
"#]
    fn fragments(&self) -> PyResult<Vec<Molecule>> {
        self.inner
            .fragments()
            .map(|fragments| {
                fragments
                    .into_iter()
                    .map(|inner| Molecule { inner })
                    .collect()
            })
            .map_err(fragment_pyerr)
    }

    #[doc = r#"
Return the largest connected fragment.
"#]
    fn largest_fragment(&self) -> PyResult<Molecule> {
        self.inner
            .largest_fragment()
            .map(|inner| Molecule { inner })
            .map_err(fragment_pyerr)
    }

    #[doc = r#"
Return the Murcko scaffold.
"#]
    fn murcko_scaffold(&self) -> PyResult<Molecule> {
        self.inner
            .murcko_scaffold()
            .map(|inner| Molecule { inner })
            .map_err(hash_pyerr)
    }

    #[doc = r#"
Return the net scaffold.
"#]
    fn net_scaffold(&self) -> PyResult<Molecule> {
        self.inner
            .net_scaffold()
            .map(|inner| Molecule { inner })
            .map_err(hash_pyerr)
    }

    #[pyo3(signature = (conf_id=-1, flavor=0))]
    #[doc = r#"
Return a PDB block string.
"#]
    fn to_pdb_block(&self, conf_id: i32, flavor: u32) -> String {
        self.inner.to_pdb_block(conf_id, flavor)
    }

    #[doc = r#"
Perceive stereochemistry and validate stereo processing for this molecule.
"#]
    fn perceive_stereochemistry(&self) -> PyResult<()> {
        self.inner.perceive_stereochemistry().map_err(stereo_pyerr)
    }

    #[doc = r#"
Serialize the molecule to COSMolKit binary form.
"#]
    fn mol_to_binary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let data = cosmolkit_core::mol_to_binary(&self.inner).map_err(pickle_pyerr)?;
        Ok(PyBytes::new(py, &data))
    }

    #[classmethod]
    #[doc = r#"
Deserialize a molecule from COSMolKit binary data.
"#]
    fn mol_from_binary(_cls: &Bound<'_, PyType>, data: &Bound<'_, PyBytes>) -> PyResult<Self> {
        let inner = cosmolkit_core::mol_from_binary(data.as_bytes()).map_err(pickle_pyerr)?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (strict=None))]
    fn sanitize(&self, strict: Option<bool>) -> PyResult<Self> {
        reject_non_strict_sanitize(strict)?;
        self.inner
            .sanitize()
            .map(|inner| Self { inner })
            .map_err(|err| PyValueError::new_err(err.to_string()))
    }

    #[pyo3(signature = (strict=None))]
    fn sanitize_(&mut self, strict: Option<bool>) -> PyResult<()> {
        reject_non_strict_sanitize(strict)?;
        self.inner
            .sanitize_()
            .map_err(|err| PyValueError::new_err(err.to_string()))
    }

    fn __len__(&self) -> usize {
        self.inner.atoms().len()
    }

    fn __repr__(&self) -> String {
        format!(
            "Molecule(num_atoms={}, num_bonds={}, has_2d_coordinates={})",
            self.inner.atoms().len(),
            self.inner.bonds().len(),
            self.inner.coordinates_2d().is_some()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl Protein {
    #[classmethod]
    #[doc = r#"
Read a PDB file as a protein-focused structural value.

The returned ``Protein`` keeps amino-acid residues and exposes chain, residue,
and atom traversal. Use ``Molecule.from_pdb_block()`` instead when the desired
result is a RDKit-compatible molecule conversion.
"#]
    fn from_pdb(_cls: &Bound<'_, PyType>, path: &str) -> PyResult<Self> {
        let path = expand_user_path(path)?;
        let inner = cosmolkit_core::Protein::from_pdb(
            path.to_str()
                .ok_or_else(|| PyValueError::new_err("path must be valid UTF-8"))?,
        )
        .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self { inner })
    }

    #[classmethod]
    #[doc = r#"
Read PDB text as a protein-focused structural value.

This is the in-memory counterpart to ``Protein.from_pdb()``.
"#]
    fn from_pdb_str(_cls: &Bound<'_, PyType>, text: &str) -> PyResult<Self> {
        let inner = cosmolkit_core::Protein::from_pdb_str(text)
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self { inner })
    }

    #[classmethod]
    #[doc = r#"
Read an mmCIF file as a protein-focused structural value.

The result uses the same protein projection as ``Protein.from_pdb()``.
"#]
    fn from_mmcif(_cls: &Bound<'_, PyType>, path: &str) -> PyResult<Self> {
        let path = expand_user_path(path)?;
        let inner = cosmolkit_core::Protein::from_mmcif(
            path.to_str()
                .ok_or_else(|| PyValueError::new_err("path must be valid UTF-8"))?,
        )
        .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self { inner })
    }

    #[classmethod]
    #[doc = r#"
Read mmCIF text as a protein-focused structural value.

``path`` is used for format context and diagnostic messages.
"#]
    fn from_mmcif_str(_cls: &Bound<'_, PyType>, text: &str, path: Option<&str>) -> PyResult<Self> {
        let inner = cosmolkit_core::Protein::from_mmcif_str(text, path.unwrap_or("input.cif"))
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self { inner })
    }

    #[doc = "Return the number of coordinate models in the protein structure."]
    fn num_models(&self) -> usize {
        self.inner.num_models()
    }

    #[doc = "Return the number of protein chains."]
    fn num_chains(&self) -> usize {
        self.inner.num_chains()
    }

    #[doc = "Return the number of protein residues."]
    fn num_residues(&self) -> usize {
        self.inner.num_residues()
    }

    #[doc = "Return the number of protein atoms."]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }

    #[doc = "Return the protein chains as ``ProteinChain`` views."]
    fn chains(&self) -> Vec<ProteinChain> {
        (0..self.inner.num_chains())
            .map(|index| ProteinChain {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }

    #[doc = "Return all protein residues as ``ProteinResidue`` views."]
    fn residues(&self) -> Vec<ProteinResidue> {
        (0..self.inner.num_residues())
            .map(|index| ProteinResidue {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }

    #[doc = "Return all protein atoms as ``ProteinAtom`` views."]
    fn atoms(&self) -> Vec<ProteinAtom> {
        (0..self.inner.num_atoms())
            .map(|index| ProteinAtom {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }

    fn __getitem__(&self, index: isize) -> PyResult<ProteinChain> {
        let len = self.inner.num_chains() as isize;
        let index = if index < 0 { len + index } else { index };
        if index < 0 || index >= len {
            return Err(PyIndexError::new_err("Protein chain index out of range"));
        }
        Ok(ProteinChain {
            inner: self.inner.clone(),
            index: index as usize,
        })
    }

    fn __len__(&self) -> usize {
        self.inner.num_chains()
    }

    fn __repr__(&self) -> String {
        format!(
            "Protein(chains={}, residues={}, atoms={})",
            self.inner.num_chains(),
            self.inner.num_residues(),
            self.inner.num_atoms()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ProteinChain {
    #[doc = "Return the zero-based chain index."]
    fn index(&self) -> usize {
        self.index
    }

    #[doc = "Return the chain kind, for example ``Protein``."]
    fn kind(&self) -> String {
        format!(
            "{:?}",
            self.inner
                .chain(self.index)
                .expect("valid protein chain")
                .kind()
        )
    }

    #[doc = "Return residues belonging to this chain."]
    fn residues(&self) -> Vec<ProteinResidue> {
        self.inner
            .chain(self.index)
            .expect("valid protein chain")
            .residues()
            .map(|residue| ProteinResidue {
                inner: self.inner.clone(),
                index: residue.id().index() as usize,
            })
            .collect()
    }

    #[doc = "Return atoms belonging to this chain."]
    fn atoms(&self) -> Vec<ProteinAtom> {
        self.inner
            .chain(self.index)
            .expect("valid protein chain")
            .atoms()
            .map(|atom| ProteinAtom {
                inner: self.inner.clone(),
                index: atom.id().index() as usize,
            })
            .collect()
    }

    fn __len__(&self) -> usize {
        self.inner
            .chain(self.index)
            .expect("valid protein chain")
            .residues()
            .count()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl ProteinResidue {
    #[doc = "Return the zero-based residue index."]
    fn index(&self) -> usize {
        self.index
    }

    #[doc = "Return the residue name, for example ``ALA``."]
    fn name(&self) -> String {
        self.inner
            .residues()
            .nth(self.index)
            .expect("valid protein residue")
            .name()
            .to_string()
    }

    #[doc = "Return the residue kind."]
    fn kind(&self) -> String {
        format!(
            "{:?}",
            self.inner
                .residues()
                .nth(self.index)
                .expect("valid protein residue")
                .kind()
        )
    }

    #[doc = "Return the Gemmi tabulated residue code as ``ResidueCode``."]
    #[gen_stub(override_return_type(type_repr = "ResidueCode"))]
    fn code<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let code = self
            .inner
            .residues()
            .nth(self.index)
            .expect("valid protein residue")
            .code();
        residue_code_enum_member(py, code)
    }

    #[doc = "Return the Gemmi-derived tabulated residue information."]
    fn info(&self) -> PyResidueInfo {
        PyResidueInfo {
            inner: self
                .inner
                .residues()
                .nth(self.index)
                .expect("valid protein residue")
                .info(),
        }
    }

    #[doc = "Return Gemmi's one-letter code for this residue."]
    fn one_letter_code(&self) -> String {
        self.inner
            .residues()
            .nth(self.index)
            .expect("valid protein residue")
            .one_letter_code()
            .to_string()
    }

    #[doc = "Return Gemmi's FASTA code for this residue."]
    fn fasta_code(&self) -> String {
        self.inner
            .residues()
            .nth(self.index)
            .expect("valid protein residue")
            .fasta_code()
            .to_string()
    }

    #[doc = "Return whether Gemmi marks this residue as standard."]
    fn is_standard(&self) -> bool {
        self.inner
            .residues()
            .nth(self.index)
            .expect("valid protein residue")
            .is_standard()
    }

    #[doc = "Return atoms belonging to this residue."]
    fn atoms(&self) -> Vec<ProteinAtom> {
        self.inner
            .residues()
            .nth(self.index)
            .expect("valid protein residue")
            .atoms()
            .map(|atom| ProteinAtom {
                inner: self.inner.clone(),
                index: atom.id().index() as usize,
            })
            .collect()
    }

    fn __len__(&self) -> usize {
        self.inner
            .residues()
            .nth(self.index)
            .expect("valid protein residue")
            .atoms()
            .count()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyResidueInfo {
    #[doc = "Return the tabulated residue code as ``ResidueCode``."]
    #[gen_stub(override_return_type(type_repr = "ResidueCode"))]
    fn code<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        residue_code_enum_member(py, self.inner.code)
    }

    #[doc = "Return the tabulated residue name."]
    fn name(&self) -> String {
        self.inner.name.to_string()
    }

    #[doc = "Return the Gemmi residue-info kind as ``ResidueInfoKind``."]
    #[gen_stub(override_return_type(type_repr = "ResidueInfoKind"))]
    fn kind<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        residue_info_kind_enum_member(py, self.inner.kind)
    }

    #[doc = "Return the Gemmi residue-info kind name."]
    fn kind_name(&self) -> String {
        residue_info_kind_name(self.inner.kind).to_string()
    }

    fn linking_type(&self) -> u8 {
        self.inner.linking_type
    }

    fn one_letter_code(&self) -> String {
        self.inner.one_letter_code.to_string()
    }

    fn hydrogen_count(&self) -> u8 {
        self.inner.hydrogen_count
    }

    fn weight(&self) -> f32 {
        self.inner.weight
    }

    fn found(&self) -> bool {
        self.inner.found()
    }

    fn is_water(&self) -> bool {
        self.inner.is_water()
    }

    fn is_dna(&self) -> bool {
        self.inner.is_dna()
    }

    fn is_rna(&self) -> bool {
        self.inner.is_rna()
    }

    fn is_nucleic_acid(&self) -> bool {
        self.inner.is_nucleic_acid()
    }

    fn is_amino_acid(&self) -> bool {
        self.inner.is_amino_acid()
    }

    fn is_buffer_or_water(&self) -> bool {
        self.inner.is_buffer_or_water()
    }

    fn is_standard(&self) -> bool {
        self.inner.is_standard()
    }

    fn fasta_code(&self) -> String {
        self.inner.fasta_code().to_string()
    }

    fn is_peptide_linking(&self) -> bool {
        self.inner.is_peptide_linking()
    }

    fn is_na_linking(&self) -> bool {
        self.inner.is_na_linking()
    }

    fn __repr__(&self) -> String {
        format!(
            "ResidueInfo(name='{}', code='{:?}', kind='{}')",
            self.inner.name,
            self.inner.code,
            residue_info_kind_name(self.inner.kind)
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ProteinAtom {
    #[doc = "Return the zero-based atom index."]
    fn index(&self) -> usize {
        self.index
    }

    #[doc = "Return the atomic number as a string."]
    fn element(&self) -> String {
        self.inner
            .atoms()
            .nth(self.index)
            .expect("valid protein atom")
            .row()
            .element
            .atomic_number()
            .to_string()
    }

    #[doc = "Return the atom name, for example ``CA``."]
    fn name(&self) -> String {
        let name = self
            .inner
            .atoms()
            .nth(self.index)
            .expect("valid protein atom")
            .row()
            .name
            .0;
        String::from_utf8_lossy(&name).trim().to_string()
    }

    #[doc = "Return ``(x, y, z)`` coordinates, or ``None`` when absent."]
    fn position(&self) -> Option<(f32, f32, f32)> {
        self.inner
            .atoms()
            .nth(self.index)
            .expect("valid protein atom")
            .position()
            .map(|[x, y, z]| (x, y, z))
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl Atom {
    fn idx(&self) -> usize {
        self.idx
    }
    fn atomic_num(&self) -> usize {
        self.atomic_num
    }
    fn formal_charge(&self) -> i8 {
        self.formal_charge
    }
    #[gen_stub(override_return_type(type_repr = "ChiralTag"))]
    fn chiral_tag<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "ChiralTag", self.chiral_tag_code)
    }
    fn chiral_tag_code(&self) -> i64 {
        self.chiral_tag_code
    }
    fn chiral_tag_name(&self) -> String {
        self.chiral_tag_name.clone()
    }
    fn isotope(&self) -> Option<u16> {
        self.isotope
    }
    fn atom_map_num(&self) -> Option<u32> {
        self.atom_map_num
    }
    fn is_aromatic(&self) -> bool {
        self.is_aromatic
    }
    fn explicit_hydrogens(&self) -> usize {
        self.explicit_hydrogens
    }
    fn no_implicit(&self) -> bool {
        self.no_implicit
    }
    fn num_radical_electrons(&self) -> usize {
        self.num_radical_electrons
    }
    fn degree(&self) -> usize {
        self.degree
    }
    fn explicit_valence(&self) -> Option<usize> {
        self.explicit_valence
    }
    fn implicit_hydrogens(&self) -> Option<usize> {
        self.implicit_hydrogens
    }
    fn total_num_hs(&self) -> Option<usize> {
        self.total_num_hs
    }
    fn total_valence(&self) -> Option<usize> {
        self.total_valence
    }

    fn __repr__(&self) -> String {
        format!(
            "Atom(idx={}, atomic_num={}, formal_charge={}, chiral_tag='{}', isotope={}, is_aromatic={}, degree={})",
            self.idx,
            self.atomic_num,
            self.formal_charge,
            self.chiral_tag_name,
            self.isotope
                .map(|x| x.to_string())
                .unwrap_or_else(|| "None".to_string()),
            self.is_aromatic,
            self.degree
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl Bond {
    fn idx(&self) -> usize {
        self.idx
    }
    fn begin_atom_idx(&self) -> usize {
        self.begin_atom_idx
    }
    fn end_atom_idx(&self) -> usize {
        self.end_atom_idx
    }
    #[gen_stub(override_return_type(type_repr = "BondOrder"))]
    fn bond_type<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "BondOrder", self.bond_type_code)
    }
    fn bond_type_code(&self) -> i64 {
        self.bond_type_code
    }
    fn bond_type_name(&self) -> String {
        self.bond_type_name.clone()
    }
    #[gen_stub(override_return_type(type_repr = "BondDirection"))]
    fn bond_dir<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "BondDirection", self.bond_dir_code)
    }
    fn bond_dir_code(&self) -> i64 {
        self.bond_dir_code
    }
    fn bond_dir_name(&self) -> String {
        self.bond_dir_name.clone()
    }
    #[gen_stub(override_return_type(type_repr = "BondStereo"))]
    fn stereo<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "BondStereo", self.stereo_code)
    }
    fn stereo_code(&self) -> i64 {
        self.stereo_code
    }
    fn stereo_name(&self) -> String {
        self.stereo_name.clone()
    }
    fn stereo_atoms(&self) -> Vec<usize> {
        self.stereo_atoms.clone()
    }
    fn is_aromatic(&self) -> bool {
        self.is_aromatic
    }

    fn __repr__(&self) -> String {
        format!(
            "Bond(idx={}, begin_atom_idx={}, end_atom_idx={}, bond_type='{}', bond_dir='{}', stereo='{}')",
            self.idx,
            self.begin_atom_idx,
            self.end_atom_idx,
            self.bond_type_name,
            self.bond_dir_name,
            self.stereo_name
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass]
#[doc = r#"
An explicit molecule editing context.

Use ``Molecule.edit()`` to create an editor, apply changes, and call
``commit()`` to receive a new ``Molecule``.

Examples
--------
Create an editor with ``mol.edit()``, apply atom and bond changes, then call
``commit()`` to produce a new ``Molecule``.
"#]
struct MoleculeEdit {
    builder: cosmolkit_core::MoleculeBuilder,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MoleculeEdit {
    #[doc = r#"
Add an atom by element symbol and return its atom index.
"#]
    fn add_atom(&mut self, element: &str) -> PyResult<usize> {
        let Some(atomic_num) = atomic_number_from_element(element) else {
            return Err(PyValueError::new_err(format!(
                "unsupported element symbol '{element}'"
            )));
        };
        let element = cosmolkit_core::Element::from_atomic_number(atomic_num).ok_or_else(|| {
            PyValueError::new_err(format!("unsupported atomic number {atomic_num}"))
        })?;
        Ok(self
            .builder
            .add_atom(cosmolkit_core::AtomSpec::new(element))
            .index())
    }

    #[pyo3(signature = (begin, end, order))]
    #[doc = r#"
Add a bond between two atom indices.

Parameters
----------
begin : int
    Begin atom index.
end : int
    End atom index.
order : {"single", "double", "triple", "aromatic", "dative", "unspecified"}
    Bond order.
"#]
    fn add_bond(&mut self, begin: usize, end: usize, order: &str) -> PyResult<()> {
        let order = match order.to_ascii_lowercase().as_str() {
            "single" => cosmolkit_core::BondOrder::Single,
            "double" => cosmolkit_core::BondOrder::Double,
            "triple" => cosmolkit_core::BondOrder::Triple,
            "aromatic" => cosmolkit_core::BondOrder::Aromatic,
            "dative" => cosmolkit_core::BondOrder::Dative,
            "hydrogen" => cosmolkit_core::BondOrder::Hydrogen,
            "unspecified" | "null" => cosmolkit_core::BondOrder::Null,
            _ => {
                return Err(PyValueError::new_err(format!(
                    "unsupported bond order '{order}'"
                )));
            }
        };
        self.builder
            .add_bond(cosmolkit_core::BondSpec::new(
                cosmolkit_core::AtomId::new(begin),
                cosmolkit_core::AtomId::new(end),
                order,
            ))
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(())
    }

    #[doc = r#"
Set an atom formal charge.
"#]
    fn set_atom_charge(&mut self, atom_index: usize, charge: i32) -> PyResult<()> {
        let charge =
            i8::try_from(charge).map_err(|_| PyValueError::new_err("charge out of i8 range"))?;
        let built = self
            .builder
            .clone()
            .build()
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        if atom_index >= built.num_atoms() {
            return Err(PyValueError::new_err("atom index out of range"));
        }
        self.builder = built.to_builder();
        self.builder
            .set_atom_formal_charge(cosmolkit_core::AtomId::new(atom_index), charge)
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(())
    }

    #[pyo3(signature = (sanitize=None))]
    #[doc = r#"
Commit staged edits and return a new molecule.
"#]
    fn commit(&mut self, sanitize: Option<bool>) -> PyResult<Molecule> {
        let inner = if sanitize.unwrap_or(true) {
            self.builder
                .clone()
                .build()
                .map_err(|err| PyValueError::new_err(err.to_string()))?
                .sanitize()
                .map_err(|err| PyValueError::new_err(err.to_string()))?
        } else {
            self.builder
                .clone()
                .build()
                .map_err(|err| PyValueError::new_err(err.to_string()))?
        };
        Ok(Molecule { inner })
    }

    fn __repr__(&self) -> String {
        format!(
            "MoleculeEdit(num_atoms={}, num_bonds={})",
            self.builder
                .clone()
                .build()
                .map(|m| m.num_atoms())
                .unwrap_or(0),
            self.builder
                .clone()
                .build()
                .map(|m| m.num_bonds())
                .unwrap_or(0)
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct Fingerprint {
    inner: cosmolkit_core::Fingerprint,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Fingerprint {
    #[doc = r#"
Return the fingerprint bit-vector length.
"#]
    fn n_bits(&self) -> usize {
        self.inner.n_bits()
    }

    #[doc = r#"
Return the sorted indexes of all on bits.
"#]
    fn on_bits(&self) -> Vec<usize> {
        self.inner.on_bits()
    }

    #[doc = r#"
Return the Tanimoto similarity to another fingerprint.
"#]
    fn tanimoto(&self, other: &Fingerprint) -> PyResult<f64> {
        self.inner.tanimoto(&other.inner).map_err(fingerprint_pyerr)
    }

    fn __len__(&self) -> usize {
        self.inner.n_bits()
    }

    fn __repr__(&self) -> String {
        format!(
            "Fingerprint(n_bits={}, on_bits={})",
            self.inner.n_bits(),
            self.inner.on_bits().len()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct MorganAdditionalOutput {
    atom_counts: Vec<u32>,
    atom_to_bits: Vec<Vec<usize>>,
    bit_info_map: BTreeMap<usize, Vec<(usize, u32)>>,
    atoms_per_bit: BTreeMap<usize, Vec<Vec<usize>>>,
}

impl From<cosmolkit_core::MorganAdditionalOutput> for MorganAdditionalOutput {
    fn from(value: cosmolkit_core::MorganAdditionalOutput) -> Self {
        Self {
            atom_counts: value.atom_counts,
            atom_to_bits: value.atom_to_bits,
            bit_info_map: value.bit_info_map,
            atoms_per_bit: value.atoms_per_bit,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganAdditionalOutput {
    #[doc = r#"
Return the number of fingerprint environments rooted at each atom.
"#]
    fn atom_counts(&self) -> Vec<u32> {
        self.atom_counts.clone()
    }

    #[doc = r#"
Return the fingerprint bit indexes associated with each atom.
"#]
    fn atom_to_bits(&self) -> Vec<Vec<usize>> {
        self.atom_to_bits.clone()
    }

    #[doc = r#"
Return ``{bit: [(atom_idx, radius), ...]}`` Morgan bit provenance.
"#]
    fn bit_info_map(&self) -> BTreeMap<usize, Vec<(usize, u32)>> {
        self.bit_info_map.clone()
    }

    #[doc = r#"
Return ``{bit: [[atom_idx, ...], ...]}`` atom environments per bit.
"#]
    fn atoms_per_bit(&self) -> BTreeMap<usize, Vec<Vec<usize>>> {
        self.atoms_per_bit.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "MorganAdditionalOutput(atom_counts={}, bit_info_bits={}, atoms_per_bit={})",
            self.atom_counts.len(),
            self.bit_info_map.len(),
            self.atoms_per_bit.len()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct MorganFingerprintResult {
    fingerprint: Fingerprint,
    additional_output: Option<MorganAdditionalOutput>,
}

impl From<cosmolkit_core::MorganFingerprintOutput> for MorganFingerprintResult {
    fn from(value: cosmolkit_core::MorganFingerprintOutput) -> Self {
        Self {
            fingerprint: Fingerprint {
                inner: value.fingerprint,
            },
            additional_output: value.additional_output.map(Into::into),
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganFingerprintResult {
    #[doc = r#"
Return the computed fingerprint.
"#]
    fn fingerprint(&self) -> Fingerprint {
        self.fingerprint.clone()
    }

    #[doc = r#"
Return the Morgan additional output collected by ``fingerprint_morgan_with_output()``.
"#]
    fn additional_output(&self) -> PyResult<MorganAdditionalOutput> {
        self.additional_output.clone().ok_or_else(|| {
            PyValueError::new_err(
                "Morgan additional output was not collected for this fingerprint result",
            )
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "MorganFingerprintResult(n_bits={}, has_additional_output={})",
            self.fingerprint.inner.n_bits(),
            self.additional_output.is_some()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct EmbedMoleculeResult {
    molecule: Molecule,
    conf_id: i32,
    params: PyEmbedParameters,
}

impl From<cosmolkit_core::EmbedMoleculeResult> for EmbedMoleculeResult {
    fn from(value: cosmolkit_core::EmbedMoleculeResult) -> Self {
        Self {
            molecule: Molecule {
                inner: value.molecule,
            },
            conf_id: value.conf_id,
            params: PyEmbedParameters::from_inner(value.params),
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl EmbedMoleculeResult {
    #[doc = r#"
Return the embedded molecule value.
"#]
    fn molecule(&self) -> Molecule {
        self.molecule.clone()
    }

    #[doc = r#"
Return the generated conformer id, or ``-1`` when embedding did not produce a conformer.
"#]
    fn conf_id(&self) -> i32 {
        self.conf_id
    }

    #[doc = r#"
Return whether embedding produced a conformer.
"#]
    fn ok(&self) -> bool {
        self.conf_id >= 0
    }

    #[doc = r#"
Return the final embedding parameters snapshot, including tracked failures.
"#]
    fn params(&self) -> PyEmbedParameters {
        self.params.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "EmbedMoleculeResult(conf_id={}, ok={})",
            self.conf_id,
            self.conf_id >= 0
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct EmbedMultipleConfsResult {
    molecule: Molecule,
    conf_ids: Vec<i32>,
    requested_num_confs: u32,
    params: PyEmbedParameters,
}

impl From<cosmolkit_core::EmbedMultipleConfsResult> for EmbedMultipleConfsResult {
    fn from(value: cosmolkit_core::EmbedMultipleConfsResult) -> Self {
        Self {
            molecule: Molecule {
                inner: value.molecule,
            },
            conf_ids: value.conf_ids,
            requested_num_confs: value.requested_num_confs,
            params: PyEmbedParameters::from_inner(value.params),
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl EmbedMultipleConfsResult {
    #[doc = r#"
Return the embedded molecule value.
"#]
    fn molecule(&self) -> Molecule {
        self.molecule.clone()
    }

    #[doc = r#"
Return the conformer ids that were kept on the returned molecule.
"#]
    fn conf_ids(&self) -> Vec<i32> {
        self.conf_ids.clone()
    }

    #[doc = r#"
Return the requested conformer count.
"#]
    fn requested_num_confs(&self) -> u32 {
        self.requested_num_confs
    }

    #[doc = r#"
Return the number of conformers kept on the returned molecule.
"#]
    fn generated_count(&self) -> usize {
        self.conf_ids.len()
    }

    #[doc = r#"
Return the final embedding parameters snapshot, including tracked failures.
"#]
    fn params(&self) -> PyEmbedParameters {
        self.params.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "EmbedMultipleConfsResult(requested_num_confs={}, generated_count={})",
            self.requested_num_confs,
            self.conf_ids.len()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct UffOptimizeMoleculeResult {
    molecule: Molecule,
    needs_more: i32,
    energy: f64,
}

impl From<cosmolkit_core::UffOptimizeMoleculeResult> for UffOptimizeMoleculeResult {
    fn from(value: cosmolkit_core::UffOptimizeMoleculeResult) -> Self {
        Self {
            molecule: Molecule {
                inner: value.molecule,
            },
            needs_more: value.needs_more,
            energy: value.energy,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl UffOptimizeMoleculeResult {
    #[doc = r#"
Return the optimized molecule value.
"#]
    fn molecule(&self) -> Molecule {
        self.molecule.clone()
    }

    #[doc = r#"
Return whether another minimization pass would still be needed.
"#]
    fn needs_more(&self) -> bool {
        self.needs_more > 0
    }

    #[doc = r#"
Return the raw RDKit-style minimization status code.
"#]
    fn status_code(&self) -> i32 {
        self.needs_more
    }

    #[doc = r#"
Return the final UFF force-field energy.
"#]
    fn energy(&self) -> f64 {
        self.energy
    }

    fn __repr__(&self) -> String {
        format!(
            "UffOptimizeMoleculeResult(needs_more={}, energy={})",
            self.needs_more, self.energy
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct UffOptimizeMoleculeConfResult {
    needs_more: i32,
    energy: f64,
}

impl From<cosmolkit_core::UffOptimizeMoleculeConfResult> for UffOptimizeMoleculeConfResult {
    fn from(value: cosmolkit_core::UffOptimizeMoleculeConfResult) -> Self {
        Self {
            needs_more: value.needs_more,
            energy: value.energy,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl UffOptimizeMoleculeConfResult {
    fn needs_more(&self) -> bool {
        self.needs_more > 0
    }

    fn status_code(&self) -> i32 {
        self.needs_more
    }

    fn energy(&self) -> f64 {
        self.energy
    }

    fn __repr__(&self) -> String {
        format!(
            "UffOptimizeMoleculeConfResult(needs_more={}, energy={})",
            self.needs_more, self.energy
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct UffOptimizeMoleculeConfsResult {
    molecule: Molecule,
    conformer_results: Vec<UffOptimizeMoleculeConfResult>,
}

impl From<cosmolkit_core::UffOptimizeMoleculeConfsResult> for UffOptimizeMoleculeConfsResult {
    fn from(value: cosmolkit_core::UffOptimizeMoleculeConfsResult) -> Self {
        Self {
            molecule: Molecule {
                inner: value.molecule,
            },
            conformer_results: value
                .conformer_results
                .into_iter()
                .map(Into::into)
                .collect(),
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl UffOptimizeMoleculeConfsResult {
    fn molecule(&self) -> Molecule {
        self.molecule.clone()
    }

    fn conformer_results(&self) -> Vec<UffOptimizeMoleculeConfResult> {
        self.conformer_results.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "UffOptimizeMoleculeConfsResult(conformers={})",
            self.conformer_results.len()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct MmffOptimizeMoleculeResult {
    molecule: Molecule,
    needs_more: i32,
}

impl From<cosmolkit_core::MmffOptimizeMoleculeResult> for MmffOptimizeMoleculeResult {
    fn from(value: cosmolkit_core::MmffOptimizeMoleculeResult) -> Self {
        Self {
            molecule: Molecule {
                inner: value.molecule,
            },
            needs_more: value.needs_more,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MmffOptimizeMoleculeResult {
    fn molecule(&self) -> Molecule {
        self.molecule.clone()
    }

    fn needs_more(&self) -> bool {
        self.needs_more > 0
    }

    fn status_code(&self) -> i32 {
        self.needs_more
    }

    fn __repr__(&self) -> String {
        format!("MmffOptimizeMoleculeResult(needs_more={})", self.needs_more)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct MmffOptimizeMoleculeConfResult {
    needs_more: i32,
    energy: f64,
}

impl From<cosmolkit_core::MmffOptimizeMoleculeConfResult> for MmffOptimizeMoleculeConfResult {
    fn from(value: cosmolkit_core::MmffOptimizeMoleculeConfResult) -> Self {
        Self {
            needs_more: value.needs_more,
            energy: value.energy,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MmffOptimizeMoleculeConfResult {
    fn needs_more(&self) -> bool {
        self.needs_more > 0
    }

    fn status_code(&self) -> i32 {
        self.needs_more
    }

    fn energy(&self) -> f64 {
        self.energy
    }

    fn __repr__(&self) -> String {
        format!(
            "MmffOptimizeMoleculeConfResult(needs_more={}, energy={})",
            self.needs_more, self.energy
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct MmffOptimizeMoleculeConfsResult {
    molecule: Molecule,
    conformer_results: Vec<MmffOptimizeMoleculeConfResult>,
}

impl From<cosmolkit_core::MmffOptimizeMoleculeConfsResult> for MmffOptimizeMoleculeConfsResult {
    fn from(value: cosmolkit_core::MmffOptimizeMoleculeConfsResult) -> Self {
        Self {
            molecule: Molecule {
                inner: value.molecule,
            },
            conformer_results: value
                .conformer_results
                .into_iter()
                .map(Into::into)
                .collect(),
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MmffOptimizeMoleculeConfsResult {
    fn molecule(&self) -> Molecule {
        self.molecule.clone()
    }

    fn conformer_results(&self) -> Vec<MmffOptimizeMoleculeConfResult> {
        self.conformer_results.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "MmffOptimizeMoleculeConfsResult(conformers={})",
            self.conformer_results.len()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "SmartsMolecule", skip_from_py_object)]
#[derive(Clone)]
struct PySmartsMolecule {
    inner: cosmolkit_core::smarts_parse::SmartsMolecule,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySmartsMolecule {
    /// Return the number of atom query nodes in the SMARTS pattern.
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }

    /// Return the number of bond query nodes in the SMARTS pattern.
    fn num_bonds(&self) -> usize {
        self.inner.bond_queries.len()
    }

    /// Return ring-closure records as ``(closure_number, atom_index)`` tuples.
    fn ring_closures(&self) -> Vec<(u8, usize)> {
        self.inner.ring_closures.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "SmartsMolecule(num_atoms={}, num_bonds={}, ring_closures={})",
            self.num_atoms(),
            self.num_bonds(),
            self.inner.ring_closures.len()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct SubstructMatchResult {
    atom_mapping: Vec<usize>,
    bond_mapping: Vec<usize>,
}

impl From<cosmolkit_core::SubstructMatchResult> for SubstructMatchResult {
    fn from(value: cosmolkit_core::SubstructMatchResult) -> Self {
        Self {
            atom_mapping: value.atom_mapping,
            bond_mapping: value.bond_mapping,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SubstructMatchResult {
    fn atom_mapping(&self) -> Vec<usize> {
        self.atom_mapping.clone()
    }
    fn bond_mapping(&self) -> Vec<usize> {
        self.bond_mapping.clone()
    }
}

#[pyfunction]
fn version() -> &'static str {
    cosmolkit_core::version()
}

#[pyfunction]
#[doc = r#"
Serialize a molecule to COSMolKit binary bytes.

Use ``mol_to_binary()`` / ``mol_from_binary()`` or the matching ``Molecule``
methods when you need an exact COSMolKit round-trip format instead of text IO.
"#]
fn mol_to_binary<'py>(py: Python<'py>, mol: &Molecule) -> PyResult<Bound<'py, PyBytes>> {
    let data = cosmolkit_core::mol_to_binary(&mol.inner).map_err(pickle_pyerr)?;
    Ok(PyBytes::new(py, &data))
}

#[pyfunction]
#[doc = r#"
Deserialize a molecule from COSMolKit binary bytes.
"#]
fn mol_from_binary(data: &[u8]) -> PyResult<Molecule> {
    let inner = cosmolkit_core::mol_from_binary(data).map_err(pickle_pyerr)?;
    Ok(Molecule { inner })
}

#[pyfunction]
#[doc = r#"
Parse SMARTS text into a ``SmartsMolecule`` query-tree value.

This exposes SMARTS parse metadata in Python. Direct SMARTS query matching is
not yet a Python API.
"#]
fn parse_smarts(smarts: &str) -> PyResult<PySmartsMolecule> {
    let inner = cosmolkit_core::smarts_parse::parse_smarts(smarts).map_err(smarts_parse_pyerr)?;
    Ok(PySmartsMolecule { inner })
}

#[pyfunction]
#[doc = r#"
Return whether UFF parameters are available for every atom in a molecule.
"#]
fn uff_has_all_molecule_params(mol: &Molecule) -> PyResult<bool> {
    cosmolkit_core::uff_has_all_molecule_params(&mol.inner).map_err(forcefield_pyerr)
}

#[pyfunction]
#[pyo3(signature = (mol, max_iters=1000, vdw_thresh=10.0, conf_id=-1, ignore_interfrag_interactions=true))]
#[doc = r#"
Optimize one existing 3D conformer with UFF and return a result object.

The input molecule is not mutated.
"#]
fn uff_optimize_molecule(
    mol: &Molecule,
    max_iters: usize,
    vdw_thresh: f64,
    conf_id: isize,
    ignore_interfrag_interactions: bool,
) -> PyResult<UffOptimizeMoleculeResult> {
    cosmolkit_core::uff_optimize_molecule(
        &mol.inner,
        max_iters,
        vdw_thresh,
        conf_id,
        ignore_interfrag_interactions,
    )
    .map(Into::into)
    .map_err(forcefield_pyerr)
}

#[pyfunction]
#[pyo3(signature = (mol, num_threads=1, max_iters=1000, vdw_thresh=10.0, ignore_interfrag_interactions=true))]
#[doc = r#"
Optimize all existing 3D conformers with UFF and return a result object.

The input molecule is not mutated.
"#]
fn uff_optimize_molecule_confs(
    mol: &Molecule,
    num_threads: i32,
    max_iters: usize,
    vdw_thresh: f64,
    ignore_interfrag_interactions: bool,
) -> PyResult<UffOptimizeMoleculeConfsResult> {
    cosmolkit_core::uff_optimize_molecule_confs(
        &mol.inner,
        num_threads,
        max_iters,
        vdw_thresh,
        ignore_interfrag_interactions,
    )
    .map(Into::into)
    .map_err(forcefield_pyerr)
}

#[pyfunction]
#[doc = r#"
Return whether MMFF94 parameters are available for a molecule.
"#]
fn mmff_has_all_molecule_params(mol: &Molecule) -> PyResult<bool> {
    cosmolkit_core::mmff_has_all_molecule_params(&mol.inner).map_err(forcefield_pyerr)
}

#[pyfunction]
#[pyo3(signature = (mol, mmff_variant="MMFF94", max_iters=200, non_bonded_thresh=100.0, conf_id=-1, ignore_interfrag_interactions=true))]
#[doc = r#"
Optimize one existing 3D conformer with MMFF and return a result object.

The input molecule is not mutated. Supported variants include ``"MMFF94"``
and ``"MMFF94S"``.
"#]
fn mmff_optimize_molecule(
    mol: &Molecule,
    mmff_variant: &str,
    max_iters: usize,
    non_bonded_thresh: f64,
    conf_id: isize,
    ignore_interfrag_interactions: bool,
) -> PyResult<MmffOptimizeMoleculeResult> {
    cosmolkit_core::mmff_optimize_molecule(
        &mol.inner,
        mmff_variant,
        max_iters,
        non_bonded_thresh,
        conf_id,
        ignore_interfrag_interactions,
    )
    .map(Into::into)
    .map_err(forcefield_pyerr)
}

#[pyfunction]
#[pyo3(signature = (mol, num_threads=1, max_iters=1000, mmff_variant="MMFF94", non_bonded_thresh=10.0, ignore_interfrag_interactions=true))]
#[doc = r#"
Optimize all existing 3D conformers with MMFF and return a result object.

The input molecule is not mutated. Supported variants include ``"MMFF94"``
and ``"MMFF94S"``.
"#]
fn mmff_optimize_molecule_confs(
    mol: &Molecule,
    num_threads: i32,
    max_iters: usize,
    mmff_variant: &str,
    non_bonded_thresh: f64,
    ignore_interfrag_interactions: bool,
) -> PyResult<MmffOptimizeMoleculeConfsResult> {
    cosmolkit_core::mmff_optimize_molecule_confs(
        &mol.inner,
        num_threads,
        max_iters,
        mmff_variant,
        non_bonded_thresh,
        ignore_interfrag_interactions,
    )
    .map(Into::into)
    .map_err(forcefield_pyerr)
}

#[pyfunction]
fn has_substruct_match(mol: &Molecule, query: &Molecule) -> bool {
    cosmolkit_core::has_substruct_match(&mol.inner, &query.inner)
}

#[pyfunction]
fn get_substruct_match(mol: &Molecule, query: &Molecule) -> Option<SubstructMatchResult> {
    cosmolkit_core::get_substruct_match(&mol.inner, &query.inner).map(Into::into)
}

#[pyfunction]
fn get_substruct_matches(mol: &Molecule, query: &Molecule) -> Vec<SubstructMatchResult> {
    cosmolkit_core::get_substruct_matches(&mol.inner, &query.inner)
        .into_iter()
        .map(Into::into)
        .collect()
}

#[pyfunction]
#[pyo3(signature = (mol, query, max_matches=1000, uniquify=true))]
fn get_substruct_matches_with_params(
    mol: &Molecule,
    query: &Molecule,
    max_matches: usize,
    uniquify: bool,
) -> Vec<SubstructMatchResult> {
    let params = cosmolkit_core::SubstructMatchParams {
        max_matches,
        uniquify,
    };
    cosmolkit_core::get_substruct_matches_with_params(&mol.inner, &query.inner, &params)
        .into_iter()
        .map(Into::into)
        .collect()
}

#[pyfunction]
#[doc = r#"
Return Gemmi-derived tabulated residue information for a residue name.
"#]
fn find_tabulated_residue(name: &str) -> PyResidueInfo {
    PyResidueInfo {
        inner: cosmolkit_core::find_tabulated_residue(name),
    }
}

#[pyfunction]
#[doc = r#"
Return the Gemmi tabulated residue index for a residue name.
"#]
fn find_tabulated_residue_idx(name: &str) -> usize {
    cosmolkit_core::find_tabulated_residue_idx(name)
}

#[pyfunction]
#[doc = r#"
Return Gemmi-derived tabulated residue information by table index.
"#]
fn get_residue_info(idx: usize) -> PyResult<PyResidueInfo> {
    let Some(inner) = cosmolkit_core::get_residue_info_checked(idx) else {
        return Err(PyIndexError::new_err(format!(
            "residue info index {idx} out of range"
        )));
    };
    Ok(PyResidueInfo { inner })
}

#[pyfunction]
#[doc = r#"
Return the Gemmi tabulated residue code for a residue name.
"#]
fn residue_code_from_name<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
    residue_code_enum_member(py, cosmolkit_core::residue_code_from_name(name))
}

#[pyfunction]
#[doc = r#"
Expand a deprecated Gemmi protein one-letter residue code alias.
"#]
fn expand_protein_one_letter(code: &str) -> PyResult<Option<String>> {
    let mut chars = code.chars();
    let Some(c) = chars.next() else {
        return Err(PyValueError::new_err(
            "code must contain exactly one character",
        ));
    };
    if chars.next().is_some() {
        return Err(PyValueError::new_err(
            "code must contain exactly one character",
        ));
    }
    Ok(cosmolkit_core::expand_protein_one_letter(c).map(str::to_string))
}

#[pyfunction]
#[doc = r#"
Expand a one-letter amino-acid, RNA, or DNA residue sequence using Gemmi's table.
"#]
fn expand_one_letter_sequence(seq: &str, kind: i64) -> PyResult<Vec<String>> {
    let kind = residue_info_kind_from_code(kind)?;
    cosmolkit_core::expand_one_letter_sequence(seq, kind).map_err(residue_info_sequence_pyerr)
}

#[pyfunction]
#[doc = r#"
Expand a deprecated Gemmi protein one-letter residue sequence alias.
"#]
fn expand_protein_one_letter_string(seq: &str) -> PyResult<Vec<String>> {
    cosmolkit_core::expand_protein_one_letter_string(seq).map_err(residue_info_sequence_pyerr)
}

#[pyfunction]
#[doc = r#"
Expand a one-letter amino-acid, RNA, or DNA residue code using Gemmi's table.
"#]
fn expand_one_letter(code: &str, kind: i64) -> PyResult<Option<String>> {
    let mut chars = code.chars();
    let Some(c) = chars.next() else {
        return Err(PyValueError::new_err(
            "code must contain exactly one character",
        ));
    };
    if chars.next().is_some() {
        return Err(PyValueError::new_err(
            "code must contain exactly one character",
        ));
    }
    let kind = residue_info_kind_from_code(kind)?;
    Ok(cosmolkit_core::expand_one_letter(c, kind).map(str::to_string))
}

#[pymodule]
fn cosmolkit(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    add_public_enums(m)?;
    add_batch_validation_error_class(m)?;
    m.add_class::<Molecule>()?;
    m.add_class::<PyEmbedParameters>()?;
    m.add_class::<Protein>()?;
    m.add_class::<ProteinChain>()?;
    m.add_class::<ProteinResidue>()?;
    m.add_class::<PyResidueInfo>()?;
    m.add_class::<ProteinAtom>()?;
    m.add_class::<MoleculeBatch>()?;
    m.add_class::<PySdfDataset>()?;
    m.add_class::<PySdfReader>()?;
    m.add_class::<PySdfRecord>()?;
    m.add_class::<PySdfRecordMetadata>()?;
    m.add_class::<PySdfDatasetIterator>()?;
    m.add_class::<PySdfBatchIterator>()?;
    m.add_class::<PySdfReaderBatchIterator>()?;
    m.add_class::<PyBatchError>()?;
    m.add_class::<PyBatchExportReport>()?;
    m.add_class::<Atom>()?;
    m.add_class::<Bond>()?;
    m.add_class::<MoleculeEdit>()?;
    m.add_class::<Fingerprint>()?;
    m.add_class::<MorganAdditionalOutput>()?;
    m.add_class::<MorganFingerprintResult>()?;
    m.add_class::<EmbedMoleculeResult>()?;
    m.add_class::<EmbedMultipleConfsResult>()?;
    m.add_class::<UffOptimizeMoleculeResult>()?;
    m.add_class::<UffOptimizeMoleculeConfResult>()?;
    m.add_class::<UffOptimizeMoleculeConfsResult>()?;
    m.add_class::<MmffOptimizeMoleculeResult>()?;
    m.add_class::<MmffOptimizeMoleculeConfResult>()?;
    m.add_class::<MmffOptimizeMoleculeConfsResult>()?;
    m.add_class::<PySmartsMolecule>()?;
    m.add_class::<SubstructMatchResult>()?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(mol_to_binary, m)?)?;
    m.add_function(wrap_pyfunction!(mol_from_binary, m)?)?;
    m.add_function(wrap_pyfunction!(parse_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(uff_has_all_molecule_params, m)?)?;
    m.add_function(wrap_pyfunction!(uff_optimize_molecule, m)?)?;
    m.add_function(wrap_pyfunction!(uff_optimize_molecule_confs, m)?)?;
    m.add_function(wrap_pyfunction!(mmff_has_all_molecule_params, m)?)?;
    m.add_function(wrap_pyfunction!(mmff_optimize_molecule, m)?)?;
    m.add_function(wrap_pyfunction!(mmff_optimize_molecule_confs, m)?)?;
    m.add_function(wrap_pyfunction!(has_substruct_match, m)?)?;
    m.add_function(wrap_pyfunction!(get_substruct_match, m)?)?;
    m.add_function(wrap_pyfunction!(get_substruct_matches, m)?)?;
    m.add_function(wrap_pyfunction!(get_substruct_matches_with_params, m)?)?;
    m.add_function(wrap_pyfunction!(find_tabulated_residue, m)?)?;
    m.add_function(wrap_pyfunction!(find_tabulated_residue_idx, m)?)?;
    m.add_function(wrap_pyfunction!(get_residue_info, m)?)?;
    m.add_function(wrap_pyfunction!(residue_code_from_name, m)?)?;
    m.add_function(wrap_pyfunction!(expand_one_letter, m)?)?;
    m.add_function(wrap_pyfunction!(expand_protein_one_letter, m)?)?;
    m.add_function(wrap_pyfunction!(expand_one_letter_sequence, m)?)?;
    m.add_function(wrap_pyfunction!(expand_protein_one_letter_string, m)?)?;
    Ok(())
}

#[cfg(feature = "stubgen")]
define_stub_info_gatherer!(stub_info);
