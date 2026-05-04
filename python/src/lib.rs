use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufReader, Write};
use std::path::PathBuf;

use cosmolkit_core::io::molblock::{self, SdfFormat};
use cosmolkit_core::io::sdf::SdfReader;
use cosmolkit_core::{BatchErrorMode, BatchRecordError, SmilesWriteParams};
use numpy::PyArray2;
use pyo3::PyErr;
use pyo3::create_exception;
use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyType};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::define_stub_info_gatherer;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
#[cfg(not(feature = "stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;
use rayon::ThreadPoolBuilder;

create_exception!(cosmolkit, BatchValidationError, PyValueError);

fn parse_batch_error_mode(errors: Option<&str>) -> PyResult<BatchErrorMode> {
    match errors.map(|s| s.to_ascii_lowercase()) {
        None => Ok(BatchErrorMode::Raise),
        Some(v) if v == "raise" => Ok(BatchErrorMode::Raise),
        Some(v) if v == "keep" => Ok(BatchErrorMode::Keep),
        Some(v) if v == "skip" => Ok(BatchErrorMode::Skip),
        Some(v) => Err(PyValueError::new_err(format!(
            "unsupported errors mode '{v}', expected one of: raise, keep, skip"
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
        all_hs_explicit,
        include_dative_bonds,
        ignore_atom_map_numbers,
        rooted_at_atom,
        ..Default::default()
    }
}

fn run_with_n_jobs<T, F>(n_jobs: Option<usize>, f: F) -> PyResult<T>
where
    T: Send,
    F: FnOnce() -> PyResult<T> + Send,
{
    match n_jobs {
        Some(0) => Err(PyValueError::new_err("n_jobs must be >= 1")),
        Some(1) => f(),
        Some(n) => ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .map_err(|err| PyValueError::new_err(format!("failed to build rayon pool: {err}")))?
            .install(f),
        None => f(),
    }
}

fn batch_validation_pyerr(error: cosmolkit_core::BatchValidationError) -> PyErr {
    BatchValidationError::new_err(format_batch_errors(&error.errors))
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

fn reject_unsanitized_sdf(sanitize: Option<bool>) -> PyResult<()> {
    if matches!(sanitize, Some(false)) {
        return Err(PyValueError::new_err(
            "sanitize=False is not implemented for SDF readers; MolBlock/SDF parsing currently finalizes chemistry during read",
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

fn unimplemented_api(name: &str) -> PyErr {
    PyNotImplementedError::new_err(format!(
        "{name} is not implemented yet in the cosmolkit Python package"
    ))
}

fn bond_order_name(order: cosmolkit_core::BondOrder) -> &'static str {
    match order {
        cosmolkit_core::BondOrder::Null => "UNSPECIFIED",
        cosmolkit_core::BondOrder::Single => "SINGLE",
        cosmolkit_core::BondOrder::Double => "DOUBLE",
        cosmolkit_core::BondOrder::Triple => "TRIPLE",
        cosmolkit_core::BondOrder::Quadruple => "QUADRUPLE",
        cosmolkit_core::BondOrder::Aromatic => "AROMATIC",
        cosmolkit_core::BondOrder::Dative => "DATIVE",
    }
}

fn chiral_tag_name(tag: cosmolkit_core::ChiralTag) -> &'static str {
    match tag {
        cosmolkit_core::ChiralTag::Unspecified => "CHI_UNSPECIFIED",
        cosmolkit_core::ChiralTag::TetrahedralCw => "CHI_TETRAHEDRAL_CW",
        cosmolkit_core::ChiralTag::TetrahedralCcw => "CHI_TETRAHEDRAL_CCW",
        cosmolkit_core::ChiralTag::TrigonalBipyramidal => "CHI_TRIGONALBIPYRAMIDAL",
    }
}

fn bond_direction_name(direction: cosmolkit_core::BondDirection) -> &'static str {
    match direction {
        cosmolkit_core::BondDirection::None => "NONE",
        cosmolkit_core::BondDirection::EndUpRight => "ENDUPRIGHT",
        cosmolkit_core::BondDirection::EndDownRight => "ENDDOWNRIGHT",
    }
}

fn bond_stereo_name(stereo: cosmolkit_core::BondStereo) -> &'static str {
    match stereo {
        cosmolkit_core::BondStereo::None => "STEREONONE",
        cosmolkit_core::BondStereo::Any => "STEREOANY",
        cosmolkit_core::BondStereo::Cis => "STEREOCIS",
        cosmolkit_core::BondStereo::Trans => "STEREOTRANS",
    }
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
) -> Result<String, cosmolkit_core::io::molblock::MolWriteError> {
    if mol.coords_2d().is_some() {
        molblock::mol_to_2d_sdf_record(mol, format)
    } else {
        molblock::mol_to_3d_sdf_record(mol, format)
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
    match name {
        "CHI_UNSPECIFIED" => Ok(cosmolkit_core::ChiralTag::Unspecified),
        "CHI_TETRAHEDRAL_CW" => Ok(cosmolkit_core::ChiralTag::TetrahedralCw),
        "CHI_TETRAHEDRAL_CCW" => Ok(cosmolkit_core::ChiralTag::TetrahedralCcw),
        other => Err(PyValueError::new_err(format!(
            "from_rdkit unsupported atom chiral tag '{other}'"
        ))),
    }
}

fn rdkit_bond_order_from_name(name: &str) -> PyResult<cosmolkit_core::BondOrder> {
    match name {
        "UNSPECIFIED" | "ZERO" => Ok(cosmolkit_core::BondOrder::Null),
        "SINGLE" => Ok(cosmolkit_core::BondOrder::Single),
        "DOUBLE" => Ok(cosmolkit_core::BondOrder::Double),
        "TRIPLE" => Ok(cosmolkit_core::BondOrder::Triple),
        "QUADRUPLE" => Ok(cosmolkit_core::BondOrder::Quadruple),
        "AROMATIC" => Ok(cosmolkit_core::BondOrder::Aromatic),
        "DATIVE" | "DATIVEL" | "DATIVER" => Ok(cosmolkit_core::BondOrder::Dative),
        other => Err(PyValueError::new_err(format!(
            "from_rdkit unsupported bond type '{other}'"
        ))),
    }
}

fn rdkit_bond_direction_from_name(name: &str) -> PyResult<cosmolkit_core::BondDirection> {
    match name {
        "NONE" => Ok(cosmolkit_core::BondDirection::None),
        "ENDUPRIGHT" => Ok(cosmolkit_core::BondDirection::EndUpRight),
        "ENDDOWNRIGHT" => Ok(cosmolkit_core::BondDirection::EndDownRight),
        other => Err(PyValueError::new_err(format!(
            "from_rdkit unsupported bond direction '{other}'"
        ))),
    }
}

fn rdkit_bond_stereo_from_name(name: &str) -> PyResult<cosmolkit_core::BondStereo> {
    match name {
        "STEREONONE" => Ok(cosmolkit_core::BondStereo::None),
        "STEREOANY" => Ok(cosmolkit_core::BondStereo::Any),
        "STEREOCIS" | "STEREOZ" => Ok(cosmolkit_core::BondStereo::Cis),
        "STEREOTRANS" | "STEREOE" => Ok(cosmolkit_core::BondStereo::Trans),
        other => Err(PyValueError::new_err(format!(
            "from_rdkit unsupported bond stereo '{other}'"
        ))),
    }
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
    chiral_tag: String,
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
    bond_type: String,
    bond_dir: String,
    stereo: String,
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
    error_type: String,
    message: String,
}

impl From<BatchRecordError> for PyBatchError {
    fn from(error: BatchRecordError) -> Self {
        Self {
            index: error.index,
            input: error.input,
            stage: error.stage,
            error_type: error.error_type,
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
``errors`` mode for invalid records, and pass ``n_jobs`` to methods that expose
parallel execution.
"#]
struct MoleculeBatch {
    inner: cosmolkit_core::MoleculeBatch,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
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
Return a short machine-readable error category.
"#]
    fn error_type(&self) -> String {
        self.error_type.clone()
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
            ("error_type".to_string(), self.error_type.clone()),
            ("message".to_string(), self.message.clone()),
        ]
    }

    fn __repr__(&self) -> String {
        format!(
            "BatchError(index={}, stage='{}', error_type='{}', message='{}')",
            self.index, self.stage, self.error_type, self.message
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
        errors: Option<&str>,
        n_jobs: Option<usize>,
    ) -> PyResult<Self> {
        let sanitize = sanitize.unwrap_or(true);
        let mode = parse_batch_error_mode(errors)?;
        run_with_n_jobs(n_jobs, move || {
            cosmolkit_core::MoleculeBatch::from_smiles_list_with_sanitize(&smiles, sanitize, mode)
                .map(|inner| Self { inner })
                .map_err(batch_validation_pyerr)
        })
    }

    #[classmethod]
    #[pyo3(signature = (sdf_records, coordinate_dim=None, errors=None, n_jobs=None))]
    #[doc = r#"
Create a batch from SDF record strings.

Parameters
----------
sdf_records : list[str]
    Individual SDF records.
coordinate_dim : {"auto", "2d", "3d"}, optional
    How coordinate columns should be interpreted.
errors : {"raise", "keep", "skip"}, optional
    Invalid-record handling mode. The default is ``"raise"``.
n_jobs : int, optional
    Number of worker threads to use.
"#]
    fn from_sdf_record_strings(
        _cls: &Bound<'_, PyType>,
        sdf_records: Vec<String>,
        coordinate_dim: Option<&str>,
        errors: Option<&str>,
        n_jobs: Option<usize>,
    ) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let coordinate_mode = parse_sdf_coordinate_mode(coordinate_dim)?;
        run_with_n_jobs(n_jobs, move || {
            cosmolkit_core::MoleculeBatch::from_sdf_record_strings(
                &sdf_records,
                coordinate_mode,
                mode,
            )
            .map(|inner| Self { inner })
            .map_err(batch_validation_pyerr)
        })
    }

    #[pyo3(signature = (errors=None, n_jobs=None))]
    #[doc = r#"
Return a new batch with explicit hydrogens added to each valid molecule.
"#]
    fn add_hydrogens(&self, errors: Option<&str>, n_jobs: Option<usize>) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let inner = &self.inner;
        run_with_n_jobs(n_jobs, move || {
            inner
                .add_hydrogens(mode)
                .map(|inner| Self { inner })
                .map_err(batch_validation_pyerr)
        })
    }

    #[pyo3(signature = (errors=None, n_jobs=None))]
    #[doc = r#"
Return a new batch with explicit hydrogens removed from each valid molecule.
"#]
    fn remove_hydrogens(&self, errors: Option<&str>, n_jobs: Option<usize>) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let inner = &self.inner;
        run_with_n_jobs(n_jobs, move || {
            inner
                .remove_hydrogens(mode)
                .map(|inner| Self { inner })
                .map_err(batch_validation_pyerr)
        })
    }

    #[pyo3(signature = (strict=None, errors=None, n_jobs=None))]
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
        errors: Option<&str>,
        n_jobs: Option<usize>,
    ) -> PyResult<Self> {
        reject_non_strict_sanitize(strict)?;
        let mode = parse_batch_error_mode(errors)?;
        let inner = &self.inner;
        run_with_n_jobs(n_jobs, move || {
            inner
                .sanitize(mode)
                .map(|inner| Self { inner })
                .map_err(batch_validation_pyerr)
        })
    }

    #[pyo3(signature = (sanitize=None, errors=None, n_jobs=None))]
    #[doc = r#"
Return a new batch with aromatic bonds converted to an explicit Kekule form.
"#]
    fn with_kekulized_bonds(
        &self,
        sanitize: Option<bool>,
        errors: Option<&str>,
        n_jobs: Option<usize>,
    ) -> PyResult<Self> {
        let sanitize = sanitize.unwrap_or(true);
        let mode = parse_batch_error_mode(errors)?;
        let inner = &self.inner;
        run_with_n_jobs(n_jobs, move || {
            inner
                .kekulize_with_sanitize(sanitize, mode)
                .map(|inner| Self { inner })
                .map_err(batch_validation_pyerr)
        })
    }

    #[pyo3(signature = (errors=None, n_jobs=None))]
    #[doc = r#"
Return a new batch with 2D coordinates computed for each valid molecule.
"#]
    fn compute_2d_coords(&self, errors: Option<&str>, n_jobs: Option<usize>) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let inner = &self.inner;
        run_with_n_jobs(n_jobs, move || {
            inner
                .compute_2d_coords(mode)
                .map(|inner| Self { inner })
                .map_err(batch_validation_pyerr)
        })
    }

    #[doc = r#"
Return a batch containing only valid molecules.
"#]
    fn filter_valid(&self) -> Self {
        Self {
            inner: self.inner.filter_valid(),
        }
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
        n_jobs=None
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
        run_with_n_jobs(n_jobs, move || {
            inner
                .to_smiles_list_with_params(&params)
                .map_err(batch_validation_pyerr)
        })
    }

    #[pyo3(signature = (n_jobs=None))]
    #[doc = r#"
Return distance-geometry bounds matrices for all valid records.
"#]
    fn dg_bounds_matrix_list(&self, n_jobs: Option<usize>) -> PyResult<Vec<Option<Vec<Vec<f64>>>>> {
        let inner = &self.inner;
        run_with_n_jobs(n_jobs, move || {
            inner
                .dg_bounds_matrix_list()
                .map_err(batch_validation_pyerr)
        })
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
        n_jobs=None
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
        run_with_n_jobs(n_jobs, move || {
            inner
                .morgan_fingerprint_list(&params)
                .map(|values| {
                    values
                        .into_iter()
                        .map(|value| value.map(|inner| Fingerprint { inner }))
                        .collect()
                })
                .map_err(batch_validation_pyerr)
        })
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
        n_jobs=None
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
        run_with_n_jobs(n_jobs, move || {
            inner
                .morgan_fingerprint_with_output_list(&params)
                .map(|values| {
                    values
                        .into_iter()
                        .map(|value| value.map(Into::into))
                        .collect()
                })
                .map_err(batch_validation_pyerr)
        })
    }

    #[pyo3(signature = (width=300, height=300, n_jobs=None))]
    #[doc = r#"
Render each valid molecule to an SVG string.
"#]
    fn to_svg_list(
        &self,
        width: u32,
        height: u32,
        n_jobs: Option<usize>,
    ) -> PyResult<Vec<Option<String>>> {
        let inner = &self.inner;
        run_with_n_jobs(n_jobs, move || {
            inner
                .to_svg_list(width, height)
                .map_err(batch_validation_pyerr)
        })
    }

    #[pyo3(signature = (out_dir, format=None, size=None, n_jobs=None, errors=None, report_path=None))]
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
        errors: Option<&str>,
        report_path: Option<&str>,
    ) -> PyResult<PyBatchExportReport> {
        let mode = parse_batch_error_mode(errors)?;
        let image_format = format.unwrap_or("png").to_string();
        let (width, height) = size.unwrap_or((300, 300));
        let out_dir = expand_user_path(out_dir)?;
        let inner = &self.inner;
        let report = run_with_n_jobs(n_jobs, move || {
            inner
                .write_images(out_dir.as_path(), &image_format, width, height, mode)
                .map_err(batch_validation_pyerr)
        })?;
        if let Some(path) = report_path {
            write_batch_report(path, &report)?;
        }
        Ok(report.into())
    }

    #[pyo3(signature = (path, format=None, errors=None, n_jobs=None, report_path=None))]
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
        errors: Option<&str>,
        n_jobs: Option<usize>,
        report_path: Option<&str>,
    ) -> PyResult<PyBatchExportReport> {
        let mode = parse_batch_error_mode(errors)?;
        let sdf_format = parse_sdf_format(format)?;
        let path = expand_user_path(path)?;
        let inner = &self.inner;
        let report = run_with_n_jobs(n_jobs, move || {
            inner
                .write_sdf(path.as_path(), sdf_format, mode)
                .map_err(batch_validation_pyerr)
        })?;
        if let Some(report_path) = report_path {
            write_batch_report(report_path, &report)?;
        }
        Ok(report.into())
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "MoleculeBatch(n={}, valid={}, invalid={})",
            self.inner.len(),
            self.inner.valid_count(),
            self.inner.invalid_count()
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
        reject_unsanitized_sdf(sanitize)?;
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
    #[pyo3(signature = (sdf_text, sanitize=None, coordinate_dim=None))]
    #[doc = r#"
Read one molecule from an SDF record string.
"#]
    fn read_sdf_record_from_str(
        _cls: &Bound<'_, PyType>,
        sdf_text: &str,
        sanitize: Option<bool>,
        coordinate_dim: Option<&str>,
    ) -> PyResult<Self> {
        reject_unsanitized_sdf(sanitize)?;
        let coordinate_mode = parse_sdf_coordinate_mode(coordinate_dim)?;
        let record = cosmolkit_core::io::sdf::read_sdf_record_from_str_with_coordinate_mode(
            sdf_text,
            coordinate_mode,
        )
        .map_err(|e| PyValueError::new_err(format!("read_sdf_record_from_str failed: {e:?}")))?;
        Ok(Self {
            inner: record.molecule,
        })
    }

    #[classmethod]
    #[pyo3(signature = (sdf_text, sanitize=None, coordinate_dim=None))]
    #[doc = r#"
Read all molecule records from an SDF string.
"#]
    fn read_sdf_records_from_str(
        _cls: &Bound<'_, PyType>,
        sdf_text: &str,
        sanitize: Option<bool>,
        coordinate_dim: Option<&str>,
    ) -> PyResult<Vec<Self>> {
        reject_unsanitized_sdf(sanitize)?;
        let coordinate_mode = parse_sdf_coordinate_mode(coordinate_dim)?;
        let records = cosmolkit_core::io::sdf::read_sdf_records_from_str_with_coordinate_mode(
            sdf_text,
            coordinate_mode,
        )
        .map_err(|e| PyValueError::new_err(format!("read_sdf_records_from_str failed: {e:?}")))?;
        Ok(records
            .into_iter()
            .map(|record| Self {
                inner: record.molecule,
            })
            .collect())
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
                chiral_tag: chiral_tag_name(atom.chiral_tag).to_string(),
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
                bond_type: bond_order_name(bond.order).to_string(),
                bond_dir: bond_direction_name(bond.direction).to_string(),
                stereo: bond_stereo_name(bond.stereo).to_string(),
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

    #[gen_stub(override_return_type(type_repr = "numpy.ndarray", imports = ("numpy")))]
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
    #[gen_stub(override_return_type(type_repr = "numpy.ndarray", imports = ("numpy")))]
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

    #[gen_stub(override_return_type(type_repr = "numpy.ndarray", imports = ("numpy")))]
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

    #[pyo3(signature = (isomeric_smiles=true))]
    #[doc = r#"
Return a SMILES string.

Parameters
----------
isomeric_smiles : bool, default True
    Include stereochemical and isotopic information when available.
"#]
    fn to_smiles(&self, isomeric_smiles: bool) -> PyResult<String> {
        self.inner.to_smiles(isomeric_smiles).map_err(|err| {
            PyNotImplementedError::new_err(format!(
                "Molecule.to_smiles(isomeric_smiles={isomeric_smiles}) is not implemented yet: {err}"
            ))
        })
    }

    #[pyo3(signature = (path, format=None))]
    #[doc = r#"
Write the molecule as an SDF file.
"#]
    fn write_sdf(&self, path: &str, format: Option<&str>) -> PyResult<()> {
        let expanded_path = expand_user_path(path)?;
        let fmt = parse_sdf_format(format)?;
        let block = molecule_to_sdf_record_string(&self.inner, fmt)
            .map_err(|err| PyValueError::new_err(format!("write_sdf failed: {err}")))?;
        let mut f = File::create(&expanded_path)
            .map_err(|e| PyValueError::new_err(format!("write_sdf create failed: {e}")))?;
        f.write_all(block.as_bytes())
            .map_err(|e| PyValueError::new_err(format!("write_sdf write failed: {e}")))?;
        Ok(())
    }

    #[pyo3(signature = (format=None))]
    #[doc = r#"
Return the molecule as an SDF record string.
"#]
    fn to_sdf_string(&self, format: Option<&str>) -> PyResult<String> {
        let fmt = parse_sdf_format(format)?;
        molecule_to_sdf_record_string(&self.inner, fmt)
            .map_err(|err| PyValueError::new_err(format!("to_sdf_string failed: {err}")))
    }

    #[pyo3(signature = (directory, file_name=None, format=None))]
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
        self.write_sdf(output_str, format)?;
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

    fn perceive_rings(&self) -> PyResult<Self> {
        Err(unimplemented_api("Molecule.perceive_rings"))
    }

    fn perceive_aromaticity(&self) -> PyResult<Self> {
        Err(unimplemented_api("Molecule.perceive_aromaticity"))
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
    fn chiral_tag(&self) -> String {
        self.chiral_tag.clone()
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
            self.chiral_tag,
            self.isotope
                .map(|x| x.to_string())
                .unwrap_or_else(|| "None".to_string()),
            self.is_aromatic,
            self.degree
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
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
    fn bond_type(&self) -> String {
        self.bond_type.clone()
    }
    fn bond_dir(&self) -> String {
        self.bond_dir.clone()
    }
    fn stereo(&self) -> String {
        self.stereo.clone()
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
            self.bond_type,
            self.bond_dir,
            self.stereo
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
#[pyclass]
struct QueryMolecule;

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
Return Morgan additional output when it was collected.
"#]
    fn additional_output(&self) -> Option<MorganAdditionalOutput> {
        self.additional_output.clone()
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
#[pyclass]
struct SubstructureMatch;

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass]
struct ValenceReport;

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass]
struct Alignment;

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Alignment {
    #[classmethod]
    #[pyo3(signature = (reference, candidates))]
    fn find_most_similar_fragment(
        _cls: &Bound<'_, PyType>,
        reference: &Molecule,
        candidates: Vec<Py<Molecule>>,
    ) -> PyResult<AlignmentResult> {
        let _ = (reference, candidates);
        Err(unimplemented_api("Alignment.find_most_similar_fragment"))
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass]
struct AlignmentResult;

#[pymodule]
fn cosmolkit(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add(
        "BatchValidationError",
        m.py().get_type::<BatchValidationError>(),
    )?;
    m.add_class::<Molecule>()?;
    m.add_class::<MoleculeBatch>()?;
    m.add_class::<PyBatchError>()?;
    m.add_class::<PyBatchExportReport>()?;
    m.add_class::<Atom>()?;
    m.add_class::<Bond>()?;
    m.add_class::<MoleculeEdit>()?;
    m.add_class::<QueryMolecule>()?;
    m.add_class::<Fingerprint>()?;
    m.add_class::<MorganAdditionalOutput>()?;
    m.add_class::<MorganFingerprintResult>()?;
    m.add_class::<SubstructureMatch>()?;
    m.add_class::<ValenceReport>()?;
    m.add_class::<Alignment>()?;
    m.add_class::<AlignmentResult>()?;
    Ok(())
}

#[cfg(feature = "stubgen")]
define_stub_info_gatherer!(stub_info);
