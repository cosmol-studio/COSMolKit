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
    PyAny, PyBool, PyDict, PyIterator, PyList, PyMapping, PyMappingProxy, PySlice, PySliceMethods,
    PyType,
};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::define_stub_info_gatherer;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
#[cfg(not(feature = "stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;
use rayon::ThreadPoolBuilder;

fn parse_batch_error_mode(errors: Option<&Bound<'_, PyAny>>) -> PyResult<BatchErrorMode> {
    let Some(errors) = errors else {
        return Ok(BatchErrorMode::Raise);
    };
    if let Ok(value) = errors.extract::<String>() {
        return match value.to_ascii_lowercase().as_str() {
            "raise" => Ok(BatchErrorMode::Raise),
            "keep" => Ok(BatchErrorMode::Keep),
            "skip" => Ok(BatchErrorMode::Skip),
            _ => Err(PyValueError::new_err(format!(
                "unsupported errors mode '{value}', expected one of: raise, keep, skip"
            ))),
        };
    }
    match errors.extract::<i64>()? {
        1 => Ok(BatchErrorMode::Raise),
        2 => Ok(BatchErrorMode::Keep),
        3 => Ok(BatchErrorMode::Skip),
        value => Err(PyValueError::new_err(format!(
            "unsupported errors mode code {value}, expected BatchErrorMode.RAISE, KEEP, or SKIP"
        ))),
    }
}

fn batch_error_type_code(error_type: &str) -> i64 {
    match error_type {
        "SmilesParseError" => 1,
        "SdfReadError" => 2,
        "AddHydrogensError" => 3,
        "RemoveHydrogensError" => 4,
        "SanitizeError" => 5,
        "KekulizeError" => 6,
        "CoordinateGenerationError" => 7,
        "SmilesWriteError" => 8,
        "DistanceGeometryError" => 9,
        "FingerprintError" => 10,
        "SvgDrawError" => 11,
        "PreparedDrawError" => 12,
        "SdfWriteError" => 13,
        "ImageExportError" => 14,
        "IoError" => 15,
        "FilenameError" => 16,
        _ => 0,
    }
}

fn batch_error_type_member_name(error_type: &str) -> &'static str {
    match error_type {
        "SmilesParseError" => "SMILES_PARSE",
        "SdfReadError" => "SDF_READ",
        "AddHydrogensError" => "ADD_HYDROGENS",
        "RemoveHydrogensError" => "REMOVE_HYDROGENS",
        "SanitizeError" => "SANITIZE",
        "KekulizeError" => "KEKULIZE",
        "CoordinateGenerationError" => "COORDINATE_GENERATION",
        "SmilesWriteError" => "SMILES_WRITE",
        "DistanceGeometryError" => "DISTANCE_GEOMETRY",
        "FingerprintError" => "FINGERPRINT",
        "SvgDrawError" => "SVG_DRAW",
        "PreparedDrawError" => "PREPARED_DRAW",
        "SdfWriteError" => "SDF_WRITE",
        "ImageExportError" => "IMAGE_EXPORT",
        "IoError" => "IO",
        "FilenameError" => "FILENAME",
        _ => "UNKNOWN",
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
        all_hs_explicit,
        include_dative_bonds,
        ignore_atom_map_numbers,
        rooted_at_atom,
        ..Default::default()
    }
}

fn run_batch_with_n_jobs<T, F>(
    n_jobs: Option<usize>,
    f: F,
) -> PyResult<Result<T, cosmolkit_core::BatchValidationError>>
where
    T: Send,
    F: FnOnce() -> Result<T, cosmolkit_core::BatchValidationError> + Send,
{
    match n_jobs {
        Some(0) => Err(PyValueError::new_err("n_jobs must be >= 1")),
        Some(1) => Ok(f()),
        Some(n) => Ok(ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .map_err(|err| PyValueError::new_err(format!("failed to build rayon pool: {err}")))?
            .install(f)),
        None => Ok(f()),
    }
}

fn validate_n_jobs(n_jobs: Option<usize>) -> PyResult<Option<usize>> {
    if matches!(n_jobs, Some(0)) {
        return Err(PyValueError::new_err("n_jobs must be >= 1"));
    }
    Ok(n_jobs)
}

fn batch_validation_pyerr(error: cosmolkit_core::BatchValidationError) -> PyErr {
    let message = format_batch_errors(&error.errors);
    Python::attach(|py| {
        let py_errors = error
            .errors
            .into_iter()
            .map(|error| Py::new(py, PyBatchError::from(error)))
            .collect::<PyResult<Vec<_>>>();
        match py_errors.and_then(|errors| {
            let errors = PyList::new(py, errors)?;
            let cls = py.import("cosmolkit")?.getattr("BatchValidationError")?;
            cls.call1((message, errors))
        }) {
            Ok(instance) => PyErr::from_value(instance),
            Err(error) => error,
        }
    })
}

fn sanitize_pyerr(error: cosmolkit_core::SanitizeError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn fingerprint_pyerr(error: cosmolkit_core::FingerprintError) -> PyErr {
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

fn format_batch_errors(errors: &[BatchRecordError]) -> String {
    let mut message = format!("batch validation failed with {} error(s)", errors.len());
    for error in errors.iter().take(5) {
        message.push_str(&format!(
            "; index={} stage={} type={} message={}",
            error.index, error.stage, error.error_type, error.message
        ));
    }
    if errors.len() > 5 {
        message.push_str("; ...");
    }
    message
}

fn json_escape(value: &str) -> String {
    let mut out = String::new();
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            _ => out.push(ch),
        }
    }
    out
}

fn write_batch_report(path: &str, report: &cosmolkit_core::BatchExportReport) -> PyResult<()> {
    let expanded_path = expand_user_path(path)?;
    let ext = expanded_path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("json")
        .to_ascii_lowercase();
    let content = if ext == "csv" {
        let mut content = String::from("index,input,stage,error_type,message\n");
        for error in &report.errors {
            content.push_str(&format!(
                "{},\"{}\",\"{}\",\"{}\",\"{}\"\n",
                error.index,
                json_escape(error.input.as_deref().unwrap_or("")),
                json_escape(&error.stage),
                json_escape(&error.error_type),
                json_escape(&error.message)
            ));
        }
        content
    } else {
        let mut content = format!(
            "{{\n  \"total\": {},\n  \"success\": {},\n  \"failed\": {},\n  \"errors\": [",
            report.total, report.success, report.failed
        );
        for (i, error) in report.errors.iter().enumerate() {
            if i > 0 {
                content.push(',');
            }
            content.push_str(&format!(
                "\n    {{\"index\": {}, \"input\": {}, \"stage\": \"{}\", \"error_type\": \"{}\", \"message\": \"{}\"}}",
                error.index,
                error
                    .input
                    .as_ref()
                    .map(|s| format!("\"{}\"", json_escape(s)))
                    .unwrap_or_else(|| "null".to_string()),
                json_escape(&error.stage),
                json_escape(&error.error_type),
                json_escape(&error.message)
            ));
        }
        content.push_str("\n  ]\n}\n");
        content
    };
    fs::write(&expanded_path, content)
        .map_err(|err| PyValueError::new_err(format!("write error report failed: {err}")))
}

fn to_python_tetrahedral_stereo(
    mol: &cosmolkit_core::Molecule,
) -> Vec<(usize, Vec<Option<usize>>)> {
    mol.tetrahedral_stereo()
        .into_iter()
        .map(|stereo| {
            let ligands = stereo
                .ligands
                .into_iter()
                .map(|ligand| match ligand {
                    cosmolkit_core::LigandRef::Atom(index) => Some(index),
                    cosmolkit_core::LigandRef::ImplicitH => None,
                })
                .collect();
            (stereo.center, ligands)
        })
        .collect()
}

fn bond_order_name(order: cosmolkit_core::BondOrder) -> &'static str {
    order.rdkit_name()
}

fn bond_order_code(order: cosmolkit_core::BondOrder) -> i64 {
    order.python_code()
}

fn chiral_tag_name(tag: cosmolkit_core::ChiralTag) -> &'static str {
    tag.rdkit_name()
}

fn chiral_tag_code(tag: cosmolkit_core::ChiralTag) -> i64 {
    tag.python_code()
}

fn bond_direction_name(direction: cosmolkit_core::BondDirection) -> &'static str {
    direction.rdkit_name()
}

fn bond_direction_code(direction: cosmolkit_core::BondDirection) -> i64 {
    direction.python_code()
}

fn bond_stereo_name(stereo: cosmolkit_core::BondStereo) -> &'static str {
    stereo.rdkit_name()
}

fn bond_stereo_code(stereo: cosmolkit_core::BondStereo) -> i64 {
    stereo.python_code()
}

fn enum_member<'py>(py: Python<'py>, enum_name: &str, code: i64) -> PyResult<Bound<'py, PyAny>> {
    let module = py.import("cosmolkit")?;
    module.getattr(enum_name)?.call1((code,))
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
            ("AROMATIC", 5),
            ("DATIVE", 6),
            ("HYDROGEN", 7),
        ],
    )?;
    add_int_enum(
        m,
        "BondDirection",
        "BOND_DIRECTION_MAP",
        &[
            ("NONE", 0),
            ("ENDUPRIGHT", 1),
            ("ENDDOWNRIGHT", 2),
            ("UNKNOWN", 3),
        ],
    )?;
    add_int_enum(
        m,
        "BondStereo",
        "BOND_STEREO_MAP",
        &[
            ("STEREONONE", 0),
            ("STEREOANY", 1),
            ("STEREOCIS", 2),
            ("STEREOTRANS", 3),
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
            ("CHI_TRIGONALBIPYRAMIDAL", 3),
        ],
    )?;
    add_int_enum_with_map_keys(
        m,
        "BatchErrorMode",
        "BATCH_ERROR_MODE_MAP",
        &[
            ("RAISE", 1, "raise"),
            ("KEEP", 2, "keep"),
            ("SKIP", 3, "skip"),
        ],
    )?;
    add_int_enum_with_map_keys(
        m,
        "BatchErrorType",
        "BATCH_ERROR_TYPE_MAP",
        &[
            ("UNKNOWN", 0, "UnknownError"),
            ("SMILES_PARSE", 1, "SmilesParseError"),
            ("SDF_READ", 2, "SdfReadError"),
            ("ADD_HYDROGENS", 3, "AddHydrogensError"),
            ("REMOVE_HYDROGENS", 4, "RemoveHydrogensError"),
            ("SANITIZE", 5, "SanitizeError"),
            ("KEKULIZE", 6, "KekulizeError"),
            ("COORDINATE_GENERATION", 7, "CoordinateGenerationError"),
            ("SMILES_WRITE", 8, "SmilesWriteError"),
            ("DISTANCE_GEOMETRY", 9, "DistanceGeometryError"),
            ("FINGERPRINT", 10, "FingerprintError"),
            ("SVG_DRAW", 11, "SvgDrawError"),
            ("PREPARED_DRAW", 12, "PreparedDrawError"),
            ("SDF_WRITE", 13, "SdfWriteError"),
            ("IMAGE_EXPORT", 14, "ImageExportError"),
            ("IO", 15, "IoError"),
            ("FILENAME", 16, "FilenameError"),
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

    def __init__(self, message, errors=None):
        super().__init__(message)
        self._errors = list(errors or [])

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
        None => Ok(SdfFormat::Auto),
        Some(v) if v == "auto" => Ok(SdfFormat::Auto),
        Some(v) if v == "v2000" || v == "v2k" => Ok(SdfFormat::V2000),
        Some(v) if v == "v3000" || v == "v3k" => Ok(SdfFormat::V3000),
        Some(v) => Err(PyValueError::new_err(format!(
            "unsupported SDF format '{v}', expected one of: auto, v2000, v3000"
        ))),
    }
}

fn parse_sdf_coordinate_mode(
    coordinate_dim: Option<&str>,
) -> PyResult<cosmolkit_core::io::sdf::SdfCoordinateMode> {
    match coordinate_dim.map(|s| s.to_ascii_lowercase()) {
        None => Ok(cosmolkit_core::io::sdf::SdfCoordinateMode::Auto),
        Some(v) if v == "auto" => Ok(cosmolkit_core::io::sdf::SdfCoordinateMode::Auto),
        Some(v) if v == "2d" => Ok(cosmolkit_core::io::sdf::SdfCoordinateMode::Force2D),
        Some(v) if v == "3d" => Ok(cosmolkit_core::io::sdf::SdfCoordinateMode::Force3D),
        Some(v) => Err(PyValueError::new_err(format!(
            "unsupported coordinate_dim '{v}', expected one of: auto, 2d, 3d"
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
        include_stereo,
        kekulize,
    };
    if mol.coords_2d().is_some() {
        molblock::mol_to_2d_sdf_record_with_params(mol, &params)
    } else if mol.coords_3d().is_some() {
        molblock::mol_to_3d_sdf_record_with_params(mol, &params)
    } else {
        Err(
            cosmolkit_core::io::molblock::MolWriteError::UnsupportedSubset(
                "SDF writing requires coordinates; call with_2d_coords() or read a molecule with 2D/3D coordinates before writing SDF",
            ),
        )
    }
}

fn molecule_to_2d_sdf_record_string(
    mol: &cosmolkit_core::Molecule,
    format: SdfFormat,
    include_stereo: bool,
    kekulize: bool,
) -> Result<String, cosmolkit_core::io::molblock::MolWriteError> {
    let params = cosmolkit_core::io::molblock::MolBlockWriteParams {
        format,
        include_stereo,
        kekulize,
    };
    if mol.coords_2d().is_some() {
        molblock::mol_to_2d_sdf_record_with_params(mol, &params)
    } else {
        let with_coords = mol.with_2d_coords()?;
        molblock::mol_to_2d_sdf_record_with_params(&with_coords, &params)
    }
}

fn molecule_to_3d_sdf_record_string(
    mol: &cosmolkit_core::Molecule,
    format: SdfFormat,
    include_stereo: bool,
    kekulize: bool,
) -> Result<String, cosmolkit_core::io::molblock::MolWriteError> {
    let params = cosmolkit_core::io::molblock::MolBlockWriteParams {
        format,
        include_stereo,
        kekulize,
    };
    if mol.coords_3d().is_some() {
        molblock::mol_to_3d_sdf_record_with_params(mol, &params)
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
    cosmolkit_core::ChiralTag::from_rdkit_name(name).ok_or_else(|| {
        PyValueError::new_err(format!("from_rdkit unsupported atom chiral tag '{name}'"))
    })
}

fn rdkit_bond_order_from_name(name: &str) -> PyResult<cosmolkit_core::BondOrder> {
    cosmolkit_core::BondOrder::from_rdkit_name(name)
        .ok_or_else(|| PyValueError::new_err(format!("from_rdkit unsupported bond type '{name}'")))
}

fn rdkit_bond_direction_from_name(name: &str) -> PyResult<cosmolkit_core::BondDirection> {
    cosmolkit_core::BondDirection::from_rdkit_name(name).ok_or_else(|| {
        PyValueError::new_err(format!("from_rdkit unsupported bond direction '{name}'"))
    })
}

fn rdkit_bond_stereo_from_name(name: &str) -> PyResult<cosmolkit_core::BondStereo> {
    cosmolkit_core::BondStereo::from_rdkit_name(name).ok_or_else(|| {
        PyValueError::new_err(format!("from_rdkit unsupported bond stereo '{name}'"))
    })
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
#[doc = r#"
A molecule value.

``Molecule`` stores atoms, bonds, stereochemistry, and optional coordinate data.
Transformation methods such as ``with_hydrogens()``, ``without_hydrogens()``,
``with_kekulized_bonds()``, and ``with_2d_coords()`` return new molecule values.
The original molecule is left unchanged.

Examples
--------
Create molecules with ``Molecule.from_smiles()``, transform them with value
methods such as ``with_2d_coords()``, then export strings, arrays, or depiction
files.
"#]
struct Molecule {
    inner: cosmolkit_core::Molecule,
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

Batch methods can keep invalid records when ``errors="keep"`` is used. In that
case, ``MoleculeBatch.errors()`` returns ``BatchError`` objects describing the
input index, processing stage, error type, and message.
"#]
struct PyBatchError {
    index: usize,
    input: Option<String>,
    stage: String,
    error_type_name: String,
    error_type_code: i64,
    error_type_member_name: String,
    message: String,
}

impl From<BatchRecordError> for PyBatchError {
    fn from(error: BatchRecordError) -> Self {
        let error_type_code = batch_error_type_code(&error.error_type);
        let error_type_member_name = batch_error_type_member_name(&error.error_type).to_string();
        Self {
            index: error.index,
            input: error.input,
            stage: error.stage,
            error_type_name: error.error_type,
            error_type_code,
            error_type_member_name,
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
    total: usize,
    success: usize,
    failed: usize,
    errors: Vec<PyBatchError>,
}

impl From<cosmolkit_core::BatchExportReport> for PyBatchExportReport {
    fn from(report: cosmolkit_core::BatchExportReport) -> Self {
        Self {
            total: report.total,
            success: report.success,
            failed: report.failed,
            errors: report.errors.into_iter().map(PyBatchError::from).collect(),
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
- ``"keep"`` keeps failed records and exposes them through ``errors()``.
- ``"skip"`` omits failed records from the returned batch or export.

Examples
--------
Construct a batch with ``MoleculeBatch.from_smiles_list()``, choose an
``errors`` mode for invalid records, and use ``with_parallel_jobs()`` when the
same worker count should apply to later batch operations.
"#]
struct MoleculeBatch {
    inner: cosmolkit_core::MoleculeBatch,
    parallel_jobs: Option<usize>,
    progress_bar: Option<bool>,
}

impl MoleculeBatch {
    fn effective_n_jobs(&self, n_jobs: Option<usize>) -> Option<usize> {
        n_jobs.or(self.parallel_jobs)
    }

    fn effective_progress_bar(&self, progress_bar: Option<bool>) -> bool {
        progress_bar.or(self.progress_bar).unwrap_or(false)
    }

    fn with_inner(&self, inner: cosmolkit_core::MoleculeBatch) -> Self {
        Self {
            inner,
            parallel_jobs: self.parallel_jobs,
            progress_bar: self.progress_bar,
        }
    }

    fn records_as_molecules(&self) -> Vec<Option<Molecule>> {
        self.inner
            .records
            .iter()
            .map(|record| match record {
                BatchRecord::Valid(molecule) => Some(Molecule {
                    inner: molecule.clone(),
                }),
                BatchRecord::Invalid(_) => None,
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
            .map(|index| self.inner.records[*index].clone())
            .collect();
        Self {
            inner: cosmolkit_core::MoleculeBatch::new(records),
            parallel_jobs: self.parallel_jobs,
            progress_bar: self.progress_bar,
        }
    }

    fn selected_batch_pyobject(&self, py: Python<'_>, indices: &[usize]) -> PyResult<Py<PyAny>> {
        Ok(Py::new(py, self.select_records(indices))?.into_any())
    }

    fn get_record_pyobject(&self, py: Python<'_>, index: usize) -> PyResult<Py<PyAny>> {
        match &self.inner.records[index] {
            BatchRecord::Valid(molecule) => Ok(Py::new(
                py,
                Molecule {
                    inner: molecule.clone(),
                },
            )?
            .into_any()),
            BatchRecord::Invalid(_) => Ok(py.None()),
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
Return the original input value when available.
"#]
    fn input(&self) -> Option<String> {
        self.input.clone()
    }

    #[doc = r#"
Return the processing stage where the error occurred.
"#]
    fn stage(&self) -> String {
        self.stage.clone()
    }

    #[doc = r#"
Return a stable machine-readable error category.
"#]
    #[gen_stub(override_return_type(type_repr = "BatchErrorType"))]
    fn error_type<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "BatchErrorType", self.error_type_code)
    }

    #[doc = r#"
Return the integer code for ``error_type()``.
"#]
    fn error_type_code(&self) -> i64 {
        self.error_type_code
    }

    #[doc = r#"
Return the external string name for this error category.
"#]
    fn error_type_name(&self) -> String {
        self.error_type_name.clone()
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
            ("input".to_string(), self.input.clone().unwrap_or_default()),
            ("stage".to_string(), self.stage.clone()),
            ("error_type".to_string(), self.error_type_name.clone()),
            ("message".to_string(), self.message.clone()),
        ]
    }

    fn __repr__(&self) -> String {
        format!(
            "BatchError(index={}, stage='{}', error_type='{}', message='{}')",
            self.index, self.stage, self.error_type_member_name, self.message
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyBatchExportReport {
    #[doc = r#"
Return the total number of records considered for export.
"#]
    fn total(&self) -> usize {
        self.total
    }

    #[doc = r#"
Return the number of records exported successfully.
"#]
    fn success(&self) -> usize {
        self.success
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
            "BatchExportReport(total={}, success={}, failed={})",
            self.total, self.success, self.failed
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
errors : {"raise", "keep", "skip"}, optional
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
        run_batch_with_n_jobs(n_jobs, move || {
            cosmolkit_core::MoleculeBatch::from_smiles_list_with_sanitize(&smiles, sanitize, mode)
                .map(|inner| Self {
                    inner,
                    parallel_jobs: None,
                    progress_bar: None,
                })
        })?
        .map_err(batch_validation_pyerr)
    }

    #[classmethod]
    #[pyo3(signature = (sdf_text, coordinate_dim=None, errors=None, n_jobs=None))]
    #[doc = r#"
Read all molecule records from an SDF string.

Parameters
----------
sdf_text : str
    SDF text containing one or more records.
coordinate_dim : {"auto", "2d", "3d"}, optional
    How coordinate columns should be interpreted.
errors : {"raise", "keep", "skip"}, optional
    Invalid-record handling mode. The default is ``"raise"``.
n_jobs : int, optional
    Number of worker threads to use.
"#]
    fn read_sdf_records_from_str(
        _cls: &Bound<'_, PyType>,
        sdf_text: &str,
        coordinate_dim: Option<&str>,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
    ) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let coordinate_mode = parse_sdf_coordinate_mode(coordinate_dim)?;
        let sdf_text = sdf_text.to_owned();
        run_batch_with_n_jobs(n_jobs, move || {
            cosmolkit_core::MoleculeBatch::read_sdf_records_from_str(
                &sdf_text,
                coordinate_mode,
                mode,
            )
            .map(|inner| Self {
                inner,
                parallel_jobs: None,
                progress_bar: None,
            })
        })?
        .map_err(batch_validation_pyerr)
    }

    #[pyo3(signature = (n_jobs))]
    #[doc = r#"
Return a new batch configured to use this worker count by default.

Pass ``None`` to clear the batch-level default and let rayon decide. Method-level
``n_jobs`` arguments still override this setting for that one call.
"#]
    fn with_parallel_jobs(&self, n_jobs: Option<usize>) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.clone(),
            parallel_jobs: validate_n_jobs(n_jobs)?,
            progress_bar: self.progress_bar,
        })
    }

    #[doc = r#"
Return the batch-level default worker count, or ``None`` when unset.
"#]
    fn parallel_jobs(&self) -> Option<usize> {
        self.parallel_jobs
    }

    #[pyo3(signature = (progress_bar))]
    #[doc = r#"
Return a new batch configured to show Rust-side progress bars by default.

Pass ``None`` to clear the batch-level default. Method-level ``progress_bar``
arguments still override this setting for that one call.
"#]
    fn with_progress_bar(&self, progress_bar: Option<bool>) -> Self {
        Self {
            inner: self.inner.clone(),
            parallel_jobs: self.parallel_jobs,
            progress_bar,
        }
    }

    #[doc = r#"
Return the batch-level progress-bar default, or ``None`` when unset.
"#]
    fn progress_bar(&self) -> Option<bool> {
        self.progress_bar
    }

    #[pyo3(signature = (errors=None, n_jobs=None, progress_bar=None))]
    #[doc = r#"
Return a new batch with explicit hydrogens added to each valid molecule.
"#]
    fn add_hydrogens(
        &self,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let inner = &self.inner;
        let parallel_jobs = self.parallel_jobs;
        let progress_bar_setting = self.progress_bar;
        let progress_bar = cosmolkit_core::batch_progress_bar(
            self.effective_progress_bar(progress_bar),
            self.inner.len(),
            "add_hydrogens",
        );
        let result = {
            let progress = progress_bar
                .as_ref()
                .map(cosmolkit_core::BatchProgressBar::callback);
            let progress = progress.as_deref();
            run_batch_with_n_jobs(self.effective_n_jobs(n_jobs), move || {
                inner
                    .add_hydrogens_with_progress(mode, progress)
                    .map(|inner| Self {
                        inner,
                        parallel_jobs,
                        progress_bar: progress_bar_setting,
                    })
            })
        };
        if let Some(progress_bar) = progress_bar.as_ref() {
            progress_bar.finish();
        }
        result?.map_err(batch_validation_pyerr)
    }

    #[pyo3(signature = (errors=None, n_jobs=None, progress_bar=None))]
    #[doc = r#"
Return a new batch with explicit hydrogens removed from each valid molecule.
"#]
    fn remove_hydrogens(
        &self,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let inner = &self.inner;
        let parallel_jobs = self.parallel_jobs;
        let progress_bar_setting = self.progress_bar;
        let progress_bar = cosmolkit_core::batch_progress_bar(
            self.effective_progress_bar(progress_bar),
            self.inner.len(),
            "remove_hydrogens",
        );
        let result = {
            let progress = progress_bar
                .as_ref()
                .map(cosmolkit_core::BatchProgressBar::callback);
            let progress = progress.as_deref();
            run_batch_with_n_jobs(self.effective_n_jobs(n_jobs), move || {
                inner
                    .remove_hydrogens_with_progress(mode, progress)
                    .map(|inner| Self {
                        inner,
                        parallel_jobs,
                        progress_bar: progress_bar_setting,
                    })
            })
        };
        if let Some(progress_bar) = progress_bar.as_ref() {
            progress_bar.finish();
        }
        result?.map_err(batch_validation_pyerr)
    }

    #[pyo3(signature = (strict=None, errors=None, n_jobs=None, progress_bar=None))]
    #[doc = r#"
Return a sanitized batch.

Parameters
----------
strict : bool, optional
    Optional strictness flag for available validation steps.
errors : {"raise", "keep", "skip"}, optional
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
        let inner = &self.inner;
        let parallel_jobs = self.parallel_jobs;
        let progress_bar_setting = self.progress_bar;
        let progress_bar = cosmolkit_core::batch_progress_bar(
            self.effective_progress_bar(progress_bar),
            self.inner.len(),
            "sanitize",
        );
        let result = {
            let progress = progress_bar
                .as_ref()
                .map(cosmolkit_core::BatchProgressBar::callback);
            let progress = progress.as_deref();
            run_batch_with_n_jobs(self.effective_n_jobs(n_jobs), move || {
                inner
                    .sanitize_with_progress(mode, progress)
                    .map(|inner| Self {
                        inner,
                        parallel_jobs,
                        progress_bar: progress_bar_setting,
                    })
            })
        };
        if let Some(progress_bar) = progress_bar.as_ref() {
            progress_bar.finish();
        }
        result?.map_err(batch_validation_pyerr)
    }

    #[pyo3(signature = (sanitize=None, errors=None, n_jobs=None, progress_bar=None))]
    #[doc = r#"
Return a new batch with aromatic bonds converted to an explicit Kekule form.
"#]
    fn with_kekulized_bonds(
        &self,
        sanitize: Option<bool>,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let sanitize = sanitize.unwrap_or(true);
        let mode = parse_batch_error_mode(errors)?;
        let inner = &self.inner;
        let parallel_jobs = self.parallel_jobs;
        let progress_bar_setting = self.progress_bar;
        let progress_bar = cosmolkit_core::batch_progress_bar(
            self.effective_progress_bar(progress_bar),
            self.inner.len(),
            "kekulize",
        );
        let result = {
            let progress = progress_bar
                .as_ref()
                .map(cosmolkit_core::BatchProgressBar::callback);
            let progress = progress.as_deref();
            run_batch_with_n_jobs(self.effective_n_jobs(n_jobs), move || {
                inner
                    .kekulize_with_sanitize_and_progress(sanitize, mode, progress)
                    .map(|inner| Self {
                        inner,
                        parallel_jobs,
                        progress_bar: progress_bar_setting,
                    })
            })
        };
        if let Some(progress_bar) = progress_bar.as_ref() {
            progress_bar.finish();
        }
        result?.map_err(batch_validation_pyerr)
    }

    #[pyo3(signature = (errors=None, n_jobs=None, progress_bar=None))]
    #[doc = r#"
Return a new batch with 2D coordinates computed for each valid molecule.
"#]
    fn compute_2d_coords(
        &self,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let inner = &self.inner;
        let parallel_jobs = self.parallel_jobs;
        let progress_bar_setting = self.progress_bar;
        let progress_bar = cosmolkit_core::batch_progress_bar(
            self.effective_progress_bar(progress_bar),
            self.inner.len(),
            "compute_2d_coords",
        );
        let result = {
            let progress = progress_bar
                .as_ref()
                .map(cosmolkit_core::BatchProgressBar::callback);
            let progress = progress.as_deref();
            run_batch_with_n_jobs(self.effective_n_jobs(n_jobs), move || {
                inner
                    .compute_2d_coords_with_progress(mode, progress)
                    .map(|inner| Self {
                        inner,
                        parallel_jobs,
                        progress_bar: progress_bar_setting,
                    })
            })
        };
        if let Some(progress_bar) = progress_bar.as_ref() {
            progress_bar.finish();
        }
        result?.map_err(batch_validation_pyerr)
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
        let inner = &self.inner;
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
        let progress_bar = cosmolkit_core::batch_progress_bar(
            self.effective_progress_bar(progress_bar),
            self.inner.len(),
            "to_smiles_list",
        );
        let result = {
            let progress = progress_bar
                .as_ref()
                .map(cosmolkit_core::BatchProgressBar::callback);
            let progress = progress.as_deref();
            run_batch_with_n_jobs(self.effective_n_jobs(n_jobs), move || {
                inner.to_smiles_list_with_params_and_progress(&params, progress)
            })
        };
        if let Some(progress_bar) = progress_bar.as_ref() {
            progress_bar.finish();
        }
        result?.map_err(batch_validation_pyerr)
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
        let inner = &self.inner;
        let progress_bar = cosmolkit_core::batch_progress_bar(
            self.effective_progress_bar(progress_bar),
            self.inner.len(),
            "dg_bounds_matrix_list",
        );
        let result = {
            let progress = progress_bar
                .as_ref()
                .map(cosmolkit_core::BatchProgressBar::callback);
            let progress = progress.as_deref();
            run_batch_with_n_jobs(self.effective_n_jobs(n_jobs), move || {
                inner.dg_bounds_matrix_list_with_progress(progress)
            })
        };
        if let Some(progress_bar) = progress_bar.as_ref() {
            progress_bar.finish();
        }
        let values = result?.map_err(batch_validation_pyerr)?;
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
        let inner = &self.inner;
        let progress_bar = cosmolkit_core::batch_progress_bar(
            self.effective_progress_bar(progress_bar),
            self.inner.len(),
            "fingerprint_morgan_list",
        );
        let result = {
            let progress = progress_bar
                .as_ref()
                .map(cosmolkit_core::BatchProgressBar::callback);
            let progress = progress.as_deref();
            run_batch_with_n_jobs(self.effective_n_jobs(n_jobs), move || {
                inner
                    .morgan_fingerprint_list_with_progress(&params, progress)
                    .map(|values| {
                        values
                            .into_iter()
                            .map(|value| value.map(|inner| Fingerprint { inner }))
                            .collect()
                    })
            })
        };
        if let Some(progress_bar) = progress_bar.as_ref() {
            progress_bar.finish();
        }
        result?.map_err(batch_validation_pyerr)
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
        let inner = &self.inner;
        let progress_bar = cosmolkit_core::batch_progress_bar(
            self.effective_progress_bar(progress_bar),
            self.inner.len(),
            "fingerprint_morgan_with_output_list",
        );
        let result = {
            let progress = progress_bar
                .as_ref()
                .map(cosmolkit_core::BatchProgressBar::callback);
            let progress = progress.as_deref();
            run_batch_with_n_jobs(self.effective_n_jobs(n_jobs), move || {
                inner
                    .morgan_fingerprint_with_output_list_with_progress(&params, progress)
                    .map(|values| {
                        values
                            .into_iter()
                            .map(|value| value.map(Into::into))
                            .collect()
                    })
            })
        };
        if let Some(progress_bar) = progress_bar.as_ref() {
            progress_bar.finish();
        }
        result?.map_err(batch_validation_pyerr)
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
        let inner = &self.inner;
        let progress_bar = cosmolkit_core::batch_progress_bar(
            self.effective_progress_bar(progress_bar),
            self.inner.len(),
            "to_svg_list",
        );
        let result = {
            let progress = progress_bar
                .as_ref()
                .map(cosmolkit_core::BatchProgressBar::callback);
            let progress = progress.as_deref();
            run_batch_with_n_jobs(self.effective_n_jobs(n_jobs), move || {
                inner.to_svg_list_with_progress(width, height, progress)
            })
        };
        if let Some(progress_bar) = progress_bar.as_ref() {
            progress_bar.finish();
        }
        result?.map_err(batch_validation_pyerr)
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
errors : {"raise", "keep", "skip"}, optional
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
        let inner = &self.inner;
        let progress_bar = cosmolkit_core::batch_progress_bar(
            self.effective_progress_bar(progress_bar),
            self.inner.len(),
            "to_images",
        );
        let result = {
            let progress = progress_bar
                .as_ref()
                .map(cosmolkit_core::BatchProgressBar::callback);
            let progress = progress.as_deref();
            run_batch_with_n_jobs(self.effective_n_jobs(n_jobs), move || {
                inner.write_images_with_progress(
                    out_dir.as_path(),
                    &image_format,
                    width,
                    height,
                    mode,
                    filenames.as_deref(),
                    progress,
                )
            })
        };
        if let Some(progress_bar) = progress_bar.as_ref() {
            progress_bar.finish();
        }
        let report = result?.map_err(batch_validation_pyerr)?;
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
errors : {"raise", "keep", "skip"}, optional
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
        let inner = &self.inner;
        let progress_bar = cosmolkit_core::batch_progress_bar(
            self.effective_progress_bar(progress_bar),
            self.inner.len(),
            "to_sdf",
        );
        let result = {
            let progress = progress_bar
                .as_ref()
                .map(cosmolkit_core::BatchProgressBar::callback);
            let progress = progress.as_deref();
            run_batch_with_n_jobs(self.effective_n_jobs(n_jobs), move || {
                inner.write_sdf_with_progress(path.as_path(), sdf_format, mode, progress)
            })
        };
        if let Some(progress_bar) = progress_bar.as_ref() {
            progress_bar.finish();
        }
        let report = result?.map_err(batch_validation_pyerr)?;
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
errors : {"raise", "keep", "skip"}, optional
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
        let inner = &self.inner;
        let progress_bar = cosmolkit_core::batch_progress_bar(
            self.effective_progress_bar(progress_bar),
            self.inner.len(),
            "to_sdf_files",
        );
        let result = {
            let progress = progress_bar
                .as_ref()
                .map(cosmolkit_core::BatchProgressBar::callback);
            let progress = progress.as_deref();
            run_batch_with_n_jobs(self.effective_n_jobs(n_jobs), move || {
                inner.write_sdf_files_with_progress(
                    out_dir.as_path(),
                    sdf_format,
                    mode,
                    filenames.as_deref(),
                    progress,
                )
            })
        };
        if let Some(progress_bar) = progress_bar.as_ref() {
            progress_bar.finish();
        }
        let report = result?.map_err(batch_validation_pyerr)?;
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
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        let list = PyList::new(py, self.records_as_molecules())?;
        PyIterator::from_object(list.as_any())
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
            self.parallel_jobs,
            self.progress_bar
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
    #[pyo3(signature = (rdmol, sanitize=None))]
    #[doc = r#"
Create a molecule from an RDKit molecule object.

Parameters
----------
rdmol : object
    An object compatible with RDKit's molecule API.
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
        let mut mol = cosmolkit_core::Molecule::new();
        let mut explicit_bond_stereo = Vec::with_capacity(bond_count);

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

            mol.add_atom(cosmolkit_core::Atom {
                index: 0,
                atomic_num,
                is_aromatic: py_method_extract(&atom, "GetIsAromatic")?,
                formal_charge,
                explicit_hydrogens,
                no_implicit: py_method_extract(&atom, "GetNoImplicit")?,
                num_radical_electrons,
                chiral_tag,
                isotope: (isotope_raw != 0).then_some(isotope_raw),
                atom_map_num: (atom_map_raw != 0).then_some(atom_map_raw),
                props: Default::default(),
                query: None,
                rdkit_cip_rank: None,
            });
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
            let molfile_query_bond_code = if py_method_extract::<bool>(&bond, "HasQuery")?
                && py_method_str(&bond, "DescribeQuery")?.trim() == "BondNull"
            {
                Some(8)
            } else {
                None
            };

            explicit_bond_stereo.push(stereo);
            mol.add_bond(cosmolkit_core::Bond {
                index: 0,
                begin_atom,
                end_atom,
                order,
                is_aromatic,
                direction,
                stereo: cosmolkit_core::BondStereo::None,
                stereo_atoms,
                molfile_query_bond_code,
                props: Default::default(),
                query: None,
            });
        }

        cosmolkit_core::assign_double_bond_stereo_from_directions(&mut mol);
        for (bond, stereo) in mol.bonds_mut().iter_mut().zip(explicit_bond_stereo) {
            if !matches!(stereo, cosmolkit_core::BondStereo::None) {
                bond.stereo = stereo;
            }
        }
        mol.rebuild_adjacency();
        let inner = if matches!(sanitize, Some(true)) {
            mol.sanitize().map_err(sanitize_pyerr)?
        } else {
            mol
        };
        Ok(Self { inner })
    }

    #[classmethod]
    #[pyo3(signature = (path, sanitize=None, coordinate_dim=None))]
    #[doc = r#"
Read the first molecule record from an SDF file.

Parameters
----------
path : str
    SDF file path.
sanitize : bool, optional
    Optional molecule preparation flag.
coordinate_dim : {"auto", "2d", "3d"}, optional
    How coordinate columns should be interpreted.
"#]
    fn read_sdf(
        _cls: &Bound<'_, PyType>,
        path: &str,
        sanitize: Option<bool>,
        coordinate_dim: Option<&str>,
    ) -> PyResult<Self> {
        reject_unsanitized_mol_reader(sanitize)?;
        let coordinate_mode = parse_sdf_coordinate_mode(coordinate_dim)?;
        let expanded_path = expand_user_path(path)?;
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
    #[pyo3(signature = (path, sanitize=None, coordinate_dim=None))]
    #[doc = r#"
Read one molecule from an MDL molfile.

Parameters
----------
path : str
    Molfile path.
sanitize : bool, optional
    Optional molecule preparation flag.
coordinate_dim : {"auto", "2d", "3d"}, optional
    How coordinate columns should be interpreted.
"#]
    fn read_mol(
        _cls: &Bound<'_, PyType>,
        path: &str,
        sanitize: Option<bool>,
        coordinate_dim: Option<&str>,
    ) -> PyResult<Self> {
        reject_unsanitized_mol_reader(sanitize)?;
        let coordinate_mode = parse_sdf_coordinate_mode(coordinate_dim)?;
        let expanded_path = expand_user_path(path)?;
        let record = cosmolkit_core::io::molfile::read_mol_file_with_coordinate_mode(
            &expanded_path,
            coordinate_mode,
        )
        .map_err(|e| PyValueError::new_err(format!("read_mol failed: {e:?}")))?;
        Ok(Self {
            inner: record.molecule,
        })
    }

    #[classmethod]
    #[pyo3(signature = (mol_text, sanitize=None, coordinate_dim=None))]
    #[doc = r#"
Read one molecule from an MDL molfile string.
"#]
    fn read_mol_from_str(
        _cls: &Bound<'_, PyType>,
        mol_text: &str,
        sanitize: Option<bool>,
        coordinate_dim: Option<&str>,
    ) -> PyResult<Self> {
        reject_unsanitized_mol_reader(sanitize)?;
        let coordinate_mode = parse_sdf_coordinate_mode(coordinate_dim)?;
        let record = cosmolkit_core::io::molfile::read_mol_record_from_str_with_coordinate_mode(
            mol_text,
            coordinate_mode,
        )
        .map_err(|e| PyValueError::new_err(format!("read_mol_from_str failed: {e:?}")))?;
        Ok(Self {
            inner: record.molecule,
        })
    }

    #[classmethod]
    #[pyo3(signature = (sdf_text, sanitize=None, coordinate_dim=None))]
    #[doc = r#"
Read one molecule from an SDF record string.
"#]
    fn read_sdf_from_str(
        _cls: &Bound<'_, PyType>,
        sdf_text: &str,
        sanitize: Option<bool>,
        coordinate_dim: Option<&str>,
    ) -> PyResult<Self> {
        reject_unsanitized_mol_reader(sanitize)?;
        let coordinate_mode = parse_sdf_coordinate_mode(coordinate_dim)?;
        let record = cosmolkit_core::io::sdf::read_sdf_from_str_with_coordinate_mode(
            sdf_text,
            coordinate_mode,
        )
        .map_err(|e| PyValueError::new_err(format!("read_sdf_from_str failed: {e:?}")))?;
        Ok(Self {
            inner: record.molecule,
        })
    }

    #[doc = r#"
Return a new molecule with explicit hydrogens added.
"#]
    fn with_hydrogens(&self) -> PyResult<Self> {
        let out = self
            .inner
            .with_hydrogens()
            .map_err(|err| PyValueError::new_err(format!("with_hydrogens failed: {err:?}")))?;
        Ok(Self { inner: out })
    }

    #[pyo3(signature = (sanitize=None))]
    #[doc = r#"
Return a new molecule with explicit hydrogens removed.
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
Return a new molecule with aromatic bonds converted to an explicit Kekule form.
"#]
    fn with_kekulized_bonds(&self, sanitize: Option<bool>) -> PyResult<Self> {
        let out = self
            .inner
            .with_kekulized_bonds(sanitize.unwrap_or(true))
            .map_err(|err| {
                PyValueError::new_err(format!("with_kekulized_bonds failed: {err:?}"))
            })?;
        Ok(Self { inner: out })
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
            degrees[bond.begin_atom] += 1;
            degrees[bond.end_atom] += 1;
        }
        self.inner
            .atoms()
            .iter()
            .map(|atom| Atom {
                idx: atom.index,
                atomic_num: atom.atomic_num as usize,
                formal_charge: atom.formal_charge,
                chiral_tag_name: chiral_tag_name(atom.chiral_tag).to_string(),
                chiral_tag_code: chiral_tag_code(atom.chiral_tag),
                isotope: atom.isotope,
                atom_map_num: atom.atom_map_num,
                is_aromatic: atom.is_aromatic,
                explicit_hydrogens: atom.explicit_hydrogens as usize,
                no_implicit: atom.no_implicit,
                num_radical_electrons: atom.num_radical_electrons as usize,
                degree: degrees[atom.index],
                explicit_valence: assignment
                    .as_ref()
                    .map(|v| v.explicit_valence[atom.index] as usize),
                implicit_hydrogens: assignment
                    .as_ref()
                    .map(|v| v.implicit_hydrogens[atom.index] as usize),
                total_num_hs: assignment.as_ref().map(|v| {
                    atom.explicit_hydrogens as usize + v.implicit_hydrogens[atom.index] as usize
                }),
                total_valence: assignment.as_ref().map(|v| {
                    v.explicit_valence[atom.index] as usize
                        + v.implicit_hydrogens[atom.index] as usize
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
                idx: bond.index,
                begin_atom_idx: bond.begin_atom,
                end_atom_idx: bond.end_atom,
                bond_type_name: bond_order_name(bond.order).to_string(),
                bond_type_code: bond_order_code(bond.order),
                bond_dir_name: bond_direction_name(bond.direction).to_string(),
                bond_dir_code: bond_direction_code(bond.direction),
                stereo_name: bond_stereo_name(bond.stereo).to_string(),
                stereo_code: bond_stereo_code(bond.stereo),
                stereo_atoms: bond.stereo_atoms.clone(),
                is_aromatic: bond.is_aromatic,
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
            .filter_map(|atom| match atom.chiral_tag {
                cosmolkit_core::ChiralTag::Unspecified => {
                    if include_unassigned {
                        Some((atom.index, "?".to_string()))
                    } else {
                        None
                    }
                }
                cosmolkit_core::ChiralTag::TetrahedralCw => {
                    Some((atom.index, "CHI_TETRAHEDRAL_CW".to_string()))
                }
                cosmolkit_core::ChiralTag::TetrahedralCcw => {
                    Some((atom.index, "CHI_TETRAHEDRAL_CCW".to_string()))
                }
                cosmolkit_core::ChiralTag::TrigonalBipyramidal => {
                    Some((atom.index, "CHI_TRIGONALBIPYRAMIDAL".to_string()))
                }
            })
            .collect()
    }

    #[doc = r#"
Return ordered tetrahedral stereo ligand records.

Each record is ``(center_atom_index, ordered_ligands)``. Implicit hydrogen is
represented as ``None``.
"#]
    fn tetrahedral_stereo(&self) -> Vec<(usize, Vec<Option<usize>>)> {
        to_python_tetrahedral_stereo(&self.inner)
    }

    #[doc = r#"
Return a new molecule with 2D coordinates.
"#]
    fn with_2d_coords(&self) -> PyResult<Self> {
        let out = self
            .inner
            .with_2d_coords()
            .map_err(|err| PyValueError::new_err(format!("with_2d_coords failed: {err}")))?;
        Ok(Self { inner: out })
    }

    #[doc = r#"
Return the number of stored 3D conformers.
"#]
    fn num_conformers(&self) -> usize {
        self.inner.num_3d_conformers()
    }

    #[doc = r#"
Return whether the molecule has 2D coordinates.
"#]
    fn has_2d_coords(&self) -> bool {
        self.inner.coords_2d().is_some()
    }

    #[gen_stub(override_return_type(type_repr = "numpy.ndarray[typing.Any, typing.Any]", imports = ("numpy", "typing")))]
    #[doc = r#"
Return 2D coordinates as a NumPy array with shape ``(num_atoms, 3)``.

The z column is zero-filled.
"#]
    fn coords_2d<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let Some(coords) = self.inner.coords_2d() else {
            return Err(PyValueError::new_err(
                "no 2D coordinates present; call with_2d_coords() first",
            ));
        };
        let rows: Vec<Vec<f64>> = coords.iter().map(|p| vec![p.x, p.y, 0.0]).collect();
        PyArray2::from_vec2(py, &rows)
            .map_err(|err| PyValueError::new_err(format!("Molecule.coords_2d failed: {err}")))
    }

    #[pyo3(signature = (conformer_index=0))]
    #[gen_stub(override_return_type(type_repr = "numpy.ndarray[typing.Any, typing.Any]", imports = ("numpy", "typing")))]
    #[doc = r#"
Return 3D coordinates as a NumPy array with shape ``(num_atoms, 3)``.
"#]
    fn coords_3d<'py>(
        &self,
        py: Python<'py>,
        conformer_index: usize,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let Some(coords) = self.inner.conformer_3d(conformer_index) else {
            return Err(PyValueError::new_err(format!(
                "no 3D conformer present at index {conformer_index}"
            )));
        };
        let rows: Vec<Vec<f64>> = coords.iter().map(|p| vec![p.x, p.y, p.z]).collect();
        PyArray2::from_vec2(py, &rows)
            .map_err(|err| PyValueError::new_err(format!("Molecule.coords_3d failed: {err}")))
    }

    #[gen_stub(override_return_type(type_repr = "numpy.ndarray[typing.Any, typing.Any]", imports = ("numpy", "typing")))]
    #[doc = r#"
Return the distance-geometry bounds matrix as a NumPy array.
"#]
    fn dg_bounds_matrix<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let rows = self.inner.dg_bounds_matrix().map_err(|err| {
            PyValueError::new_err(format!("Molecule.dg_bounds_matrix failed: {err}"))
        })?;
        PyArray2::from_vec2(py, &rows).map_err(|err| {
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
Write the molecule as an SDF file.
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
Write the molecule to an SDF file inside a directory.

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
            working: self.inner.clone(),
        }
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

    #[pyo3(signature = (strict=None))]
    fn sanitize(&self, strict: Option<bool>) -> PyResult<Self> {
        reject_non_strict_sanitize(strict)?;
        self.inner
            .sanitize()
            .map(|inner| Self { inner })
            .map_err(sanitize_pyerr)
    }

    fn __len__(&self) -> usize {
        self.inner.atoms().len()
    }

    fn __repr__(&self) -> String {
        format!(
            "Molecule(num_atoms={}, num_bonds={}, has_2d_coords={})",
            self.inner.atoms().len(),
            self.inner.bonds().len(),
            self.inner.coords_2d().is_some()
        )
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
    working: cosmolkit_core::Molecule,
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
        let idx = self.working.add_atom(cosmolkit_core::Atom {
            index: 0,
            atomic_num,
            is_aromatic: false,
            formal_charge: 0,
            explicit_hydrogens: 0,
            no_implicit: false,
            num_radical_electrons: 0,
            chiral_tag: cosmolkit_core::ChiralTag::Unspecified,
            isotope: None,
            atom_map_num: None,
            props: Default::default(),
            query: None,
            rdkit_cip_rank: None,
        });
        Ok(idx)
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
        if begin >= self.working.atoms().len() || end >= self.working.atoms().len() {
            return Err(PyValueError::new_err("bond atom index out of range"));
        }
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
        self.working.add_bond(cosmolkit_core::Bond {
            index: 0,
            begin_atom: begin,
            end_atom: end,
            order,
            is_aromatic: matches!(order, cosmolkit_core::BondOrder::Aromatic),
            direction: cosmolkit_core::BondDirection::None,
            stereo: cosmolkit_core::BondStereo::None,
            stereo_atoms: Vec::new(),
            molfile_query_bond_code: None,
            props: Default::default(),
            query: None,
        });
        Ok(())
    }

    #[doc = r#"
Set an atom formal charge.
"#]
    fn set_atom_charge(&mut self, atom_index: usize, charge: i32) -> PyResult<()> {
        let atom = self
            .working
            .atoms_mut()
            .get_mut(atom_index)
            .ok_or_else(|| PyValueError::new_err("atom index out of range"))?;
        let charge =
            i8::try_from(charge).map_err(|_| PyValueError::new_err("charge out of i8 range"))?;
        atom.formal_charge = charge;
        Ok(())
    }

    #[pyo3(signature = (sanitize=None))]
    #[doc = r#"
Commit staged edits and return a new molecule.
"#]
    fn commit(&mut self, sanitize: Option<bool>) -> PyResult<Molecule> {
        let inner = if sanitize.unwrap_or(true) {
            self.working.sanitize().map_err(sanitize_pyerr)?
        } else {
            self.working.clone()
        };
        Ok(Molecule { inner })
    }

    fn __repr__(&self) -> String {
        format!(
            "MoleculeEdit(num_atoms={}, num_bonds={})",
            self.working.atoms().len(),
            self.working.bonds().len()
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

#[pymodule]
fn cosmolkit(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    add_public_enums(m)?;
    add_batch_validation_error_class(m)?;
    m.add_class::<Molecule>()?;
    m.add_class::<MoleculeBatch>()?;
    m.add_class::<PyBatchError>()?;
    m.add_class::<PyBatchExportReport>()?;
    m.add_class::<Atom>()?;
    m.add_class::<Bond>()?;
    m.add_class::<MoleculeEdit>()?;
    m.add_class::<Fingerprint>()?;
    m.add_class::<MorganAdditionalOutput>()?;
    m.add_class::<MorganFingerprintResult>()?;
    Ok(())
}

#[cfg(feature = "stubgen")]
define_stub_info_gatherer!(stub_info);
