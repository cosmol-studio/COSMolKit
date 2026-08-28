use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufReader, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex, MutexGuard};

use cosmolkit_core::io::molblock::{self, SdfFormat};
use cosmolkit_core::io::sdf::SdfReader;
use cosmolkit_core::{BatchErrorMode, BatchRecord, BatchRecordError, SmilesWriteParams};
use numpy::{
    AllowTypeChange, PyArray2, PyArrayLike, PyUntypedArray, PyUntypedArrayMethods, ndarray::Ix2,
};
use pyo3::PyErr;
use pyo3::exceptions::{
    PyIndexError, PyNotImplementedError, PyOSError, PyOverflowError, PyRuntimeError, PyTypeError,
    PyValueError,
};
use pyo3::prelude::*;
use pyo3::types::{
    PyAny, PyBool, PyBytes, PyDict, PyIterator, PyList, PyMapping, PyMappingProxy, PySlice,
    PySliceMethods, PyTuple, PyType,
};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::define_stub_info_gatherer;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyfunction, gen_stub_pymethods};
#[cfg(not(feature = "stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;

mod confseq_py;

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

fn make_sdf_read_params(
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    strict_parsing: Option<bool>,
    coordinate_dim: &str,
) -> PyResult<cosmolkit_core::SdfReadParams> {
    Ok(cosmolkit_core::SdfReadParams {
        sanitize: sanitize.unwrap_or(true),
        remove_hs: remove_hs.unwrap_or(true),
        strict_parsing: strict_parsing.unwrap_or(true),
        coordinate_mode: parse_coordinate_mode(Some(coordinate_dim))?,
        ..Default::default()
    })
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

fn inchi_error_kind_name(kind: cosmolkit_core::InchiErrorKind) -> &'static str {
    match kind {
        cosmolkit_core::InchiErrorKind::AllocationFailed => "allocation_failed",
        cosmolkit_core::InchiErrorKind::UnsupportedState => "unsupported_state",
        cosmolkit_core::InchiErrorKind::InvalidInput => "invalid_input",
        cosmolkit_core::InchiErrorKind::InvalidSourceOutput => "invalid_source_output",
        cosmolkit_core::InchiErrorKind::SanitizeFailed => "sanitize_failed",
        cosmolkit_core::InchiErrorKind::Toolkit => "toolkit",
        cosmolkit_core::InchiErrorKind::SourcePort => "source_port",
    }
}

fn inchi_pyerr(error: cosmolkit_core::InchiError) -> PyErr {
    let message = error.to_string();
    Python::attach(|py| {
        let class_name = match error.kind {
            cosmolkit_core::InchiErrorKind::AllocationFailed => "InchiAllocationError",
            cosmolkit_core::InchiErrorKind::UnsupportedState => "InchiUnsupportedStateError",
            _ => "InchiError",
        };
        match (|| -> PyResult<Bound<'_, PyAny>> {
            let cls = py.import("cosmolkit")?.getattr(class_name)?;
            cls.call1((
                message,
                error.operation,
                inchi_error_kind_name(error.kind),
                error.detail,
            ))
        })() {
            Ok(instance) => PyErr::from_value(instance),
            Err(error) => error,
        }
    })
}

fn emit_inchi_diagnostics(diagnostics: &[cosmolkit_core::InchiDiagnostic]) -> PyResult<()> {
    if diagnostics.is_empty() {
        return Ok(());
    }
    Python::attach(|py| {
        let warning_class = py.import("cosmolkit")?.getattr("InchiDiagnosticWarning")?;
        let warn = py.import("warnings")?.getattr("warn")?;
        for diagnostic in diagnostics {
            let level = match diagnostic.level {
                cosmolkit_core::InchiDiagnosticLevel::Warning => "warning",
                cosmolkit_core::InchiDiagnosticLevel::Error => "error",
            };
            let warning = warning_class.call1((level, diagnostic.message.as_str()))?;
            warn.call1((warning,))?;
        }
        Ok(())
    })
}

fn inchi_output_string(
    operation: &'static str,
    field: &'static str,
    bytes: Vec<u8>,
) -> PyResult<String> {
    String::from_utf8(bytes).map_err(|error| {
        inchi_pyerr(cosmolkit_core::InchiError {
            operation,
            kind: cosmolkit_core::InchiErrorKind::InvalidSourceOutput,
            detail: format!("{field} is not valid UTF-8: {error}"),
        })
    })
}

fn fingerprint_pyerr(error: cosmolkit_core::FingerprintError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn normalize_fingerprint_indices<T>(values: Option<Vec<T>>) -> Option<Vec<T>> {
    values.filter(|values| !values.is_empty())
}

fn descriptor_pyerr(error: cosmolkit_core::DescriptorError) -> PyErr {
    PyNotImplementedError::new_err(error.to_string())
}

fn parse_rotatable_bonds_mode(mode: &str) -> PyResult<cosmolkit_core::NumRotatableBondsOptions> {
    match mode {
        "default" => Ok(cosmolkit_core::NumRotatableBondsOptions::Default),
        "non_strict" => Ok(cosmolkit_core::NumRotatableBondsOptions::NonStrict),
        "strict" => Ok(cosmolkit_core::NumRotatableBondsOptions::Strict),
        "strict_linkages" => Ok(cosmolkit_core::NumRotatableBondsOptions::StrictLinkages),
        _ => Err(PyValueError::new_err(format!(
            "unsupported rotatable-bond mode '{mode}', expected one of: default, non_strict, strict, strict_linkages"
        ))),
    }
}

fn svg_draw_pyerr(error: cosmolkit_core::SvgDrawError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn pdb_molecule_pyerr(error: cosmolkit_core::StructureMoleculeConversionError) -> PyErr {
    match error {
        cosmolkit_core::StructureMoleculeConversionError::Unsupported(message) => {
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

fn molecule_pickle_state<'py>(
    py: Python<'py>,
    inner: &cosmolkit_core::Molecule,
) -> PyResult<Bound<'py, PyDict>> {
    let data = cosmolkit_core::mol_to_binary(inner).map_err(pickle_pyerr)?;
    let state = PyDict::new(py);
    state.set_item("kind", "cosmolkit.Molecule")?;
    state.set_item("pickle_schema", 1u16)?;
    state.set_item("cosmolkit_version", env!("CARGO_PKG_VERSION"))?;
    state.set_item("core_format", "cosmolkit-molecule-archive")?;
    state.set_item("payload", PyBytes::new(py, &data))?;
    Ok(state)
}

fn stereo_pyerr(error: cosmolkit_core::StereoError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn enumeration_pyerr(error: cosmolkit_core::EnumerationError) -> PyErr {
    match &error {
        cosmolkit_core::EnumerationError::Operation(
            cosmolkit_core::OperationError::Unsupported { .. }
            | cosmolkit_core::OperationError::UnsupportedFeature { .. },
        ) => PyNotImplementedError::new_err(error.to_string()),
        _ => PyValueError::new_err(error.to_string()),
    }
}

fn potential_stereo_pyerr(error: cosmolkit_core::PotentialStereoError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn chiral_tag_assignment_pyerr(error: cosmolkit_core::OperationError) -> PyErr {
    match &error {
        cosmolkit_core::OperationError::Unsupported { .. }
        | cosmolkit_core::OperationError::UnsupportedFeature { .. } => {
            PyNotImplementedError::new_err(error.to_string())
        }
        cosmolkit_core::OperationError::Stereo {
            source: cosmolkit_core::StereoError::ConformerNotFound { .. },
            ..
        }
        | cosmolkit_core::OperationError::InvalidInput { .. }
        | cosmolkit_core::OperationError::Precondition { .. }
        | cosmolkit_core::OperationError::Valence { .. } => {
            PyValueError::new_err(error.to_string())
        }
        cosmolkit_core::OperationError::Stereo { .. }
        | cosmolkit_core::OperationError::InvariantViolation { .. } => {
            PyRuntimeError::new_err(error.to_string())
        }
        _ => PyValueError::new_err(error.to_string()),
    }
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

#[allow(clippy::too_many_arguments)]
fn make_atom_pair_fingerprint_params(
    n_bits: usize,
    min_distance: u32,
    max_distance: u32,
    use_2d: bool,
    include_chirality: bool,
    count_simulation: bool,
    count_bounds: Option<Vec<u32>>,
    num_bits_per_feature: u32,
    from_atoms: Option<Vec<usize>>,
    ignore_atoms: Option<Vec<usize>>,
    conformer_id: i32,
    custom_atom_invariants: Option<Vec<u32>>,
    collect_additional_output: bool,
) -> cosmolkit_core::AtomPairFingerprintParams {
    cosmolkit_core::AtomPairFingerprintParams {
        n_bits,
        min_distance,
        max_distance,
        use_2d,
        use_chirality: include_chirality,
        count_simulation,
        count_bounds: count_bounds.unwrap_or_else(|| vec![1, 2, 4, 8]),
        num_bits_per_feature,
        from_atoms,
        ignore_atoms,
        conformer_id,
        custom_atom_invariants,
        collect_additional_output,
    }
}

fn atom_pair_call_arguments(
    params: &cosmolkit_core::AtomPairFingerprintParams,
) -> cosmolkit_core::FingerprintFuncArguments {
    cosmolkit_core::FingerprintFuncArguments {
        from_atoms: params.from_atoms.clone(),
        ignore_atoms: params.ignore_atoms.clone(),
        conf_id: params.conformer_id,
        custom_atom_invariants: params.custom_atom_invariants.clone(),
        ..Default::default()
    }
}

fn make_avalon_fingerprint_params(
    n_bits: u32,
    is_query: bool,
    bit_flags: u32,
) -> cosmolkit_core::avalon_fingerprint::AvalonFingerprintParams {
    cosmolkit_core::avalon_fingerprint::AvalonFingerprintParams {
        n_bits,
        is_query,
        bit_flags: cosmolkit_core::AvalonFingerprintFlags::from_bits_retain(bit_flags),
    }
}

#[allow(clippy::too_many_arguments)]
fn make_layered_fingerprint_params(
    layers: u32,
    min_path: u32,
    max_path: u32,
    fp_size: u32,
    atom_counts: Option<Vec<u32>>,
    set_only_bits: Option<&Fingerprint>,
    branched_paths: bool,
    from_atoms: Option<Vec<u32>>,
) -> cosmolkit_core::LayeredFingerprintParams {
    cosmolkit_core::LayeredFingerprintParams {
        layers: cosmolkit_core::LayeredFingerprintLayers::from_bits_retain(layers),
        min_path,
        max_path,
        fp_size,
        atom_counts,
        set_only_bits: set_only_bits.map(|fingerprint| fingerprint.inner.clone()),
        branched_paths,
        from_atoms,
    }
}

#[allow(clippy::too_many_arguments)]
fn make_topological_fingerprint_params(
    min_path: u32,
    max_path: u32,
    fp_size: u32,
    num_bits_per_feature: u32,
    use_hs: bool,
    target_density: f64,
    min_size: u32,
    branched_paths: bool,
    use_bond_order: bool,
    atom_invariants: Option<Vec<u32>>,
    from_atoms: Option<Vec<u32>>,
) -> cosmolkit_core::fingerprint::TopologicalFingerprintParams {
    cosmolkit_core::fingerprint::TopologicalFingerprintParams {
        min_path,
        max_path,
        fp_size,
        num_bits_per_feature,
        use_hs,
        target_density,
        min_size,
        branched_paths,
        use_bond_order,
        atom_invariants,
        from_atoms,
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
    kind.name()
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

fn element_enum_member<'py>(
    py: Python<'py>,
    element: cosmolkit_core::Element,
) -> PyResult<Bound<'py, PyAny>> {
    enum_member(py, "Element", i64::from(element.atomic_number()))
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
    let element_members = cosmolkit_core::ELEMENTS_WITH_DUMMY
        .iter()
        .map(|element| {
            let name = if *element == cosmolkit_core::Element::DUMMY {
                "DUMMY".to_string()
            } else {
                element.symbol().to_ascii_uppercase()
            };
            (name, i64::from(element.atomic_number()))
        })
        .collect::<Vec<_>>();
    let element_member_refs = element_members
        .iter()
        .map(|(name, value)| (name.as_str(), *value))
        .collect::<Vec<_>>();
    let mut element_aliases = cosmolkit_core::ELEMENTS_WITH_DUMMY
        .iter()
        .map(|element| {
            let member_name = if *element == cosmolkit_core::Element::DUMMY {
                "DUMMY".to_string()
            } else {
                element.symbol().to_ascii_uppercase()
            };
            (element.symbol().to_string(), member_name)
        })
        .collect::<Vec<_>>();
    element_aliases.extend([
        ("Uut".to_string(), "NH".to_string()),
        ("Uup".to_string(), "MC".to_string()),
    ]);
    let element_alias_refs = element_aliases
        .iter()
        .map(|(symbol, name)| (symbol.as_str(), name.as_str()))
        .collect::<Vec<_>>();
    add_int_enum_with_map_aliases(
        m,
        "Element",
        "ELEMENT_MAP",
        &element_member_refs,
        &element_alias_refs,
    )?;

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

fn add_inchi_error_classes(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let globals = PyDict::new(py);
    globals.set_item("ValueError", py.get_type::<PyValueError>())?;
    globals.set_item(
        "UserWarning",
        py.import("builtins")?.getattr("UserWarning")?,
    )?;
    let code = r#"
class InchiError(ValueError):
    __module__ = "cosmolkit"

    def __init__(self, message, operation, kind, detail):
        super().__init__(message)
        self.operation = operation
        self.kind = kind
        self.detail = detail

class InchiAllocationError(InchiError):
    __module__ = "cosmolkit"

class InchiUnsupportedStateError(InchiError):
    __module__ = "cosmolkit"

class InchiDiagnosticWarning(UserWarning):
    __module__ = "cosmolkit"

    def __init__(self, level, message):
        super().__init__(message)
        self.level = level
        self.message = message
"#;
    py.import("builtins")?
        .getattr("exec")?
        .call1((code, &globals))?;
    for class_name in [
        "InchiError",
        "InchiAllocationError",
        "InchiUnsupportedStateError",
        "InchiDiagnosticWarning",
    ] {
        let cls = globals
            .get_item(class_name)?
            .ok_or_else(|| PyValueError::new_err(format!("failed to create {class_name} class")))?;
        m.add(class_name, cls)?;
    }
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

fn parse_z_policy(value: &str) -> PyResult<&'static str> {
    match value.to_ascii_lowercase().as_str() {
        "ignore" => Ok("ignore"),
        "require_zero" => Ok("require_zero"),
        "error" => Ok("error"),
        _ => Err(PyValueError::new_err(format!(
            "unsupported z_policy '{value}', expected one of: ignore, require_zero, error"
        ))),
    }
}

fn extract_coordinate_matrix(
    coords: &Bound<'_, PyAny>,
    expected_rows: usize,
    expected_columns: &[usize],
    label: &str,
) -> PyResult<Vec<Vec<f64>>> {
    if let Ok(array) = coords.cast::<PyUntypedArray>()
        && let [actual_rows, _] = array.shape()
        && *actual_rows != expected_rows
    {
        return Err(PyValueError::new_err(format!(
            "{label} row count mismatch: expected {expected_rows}, got {actual_rows}"
        )));
    }
    let array_like = coords
        .extract::<PyArrayLike<'_, f64, Ix2, AllowTypeChange>>()
        .map_err(|err| {
            PyTypeError::new_err(format!("{label} must be a 2D numeric array: {err}"))
        })?;
    let array = array_like.as_array();
    let shape = array.shape();
    if shape[0] != expected_rows {
        return Err(PyValueError::new_err(format!(
            "{label} row count mismatch: expected {expected_rows}, got {}",
            shape[0]
        )));
    }
    if !expected_columns.contains(&shape[1]) {
        let columns = expected_columns
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(" or ");
        return Err(PyValueError::new_err(format!(
            "{label} must have shape (num_atoms, {columns}); got ({}, {})",
            shape[0], shape[1]
        )));
    }
    let mut out = Vec::with_capacity(shape[0]);
    for (row_idx, row) in array.outer_iter().enumerate() {
        let mut values = Vec::with_capacity(shape[1]);
        for (col_idx, value) in row.iter().enumerate() {
            if !value.is_finite() {
                return Err(PyValueError::new_err(format!(
                    "{label} contains a non-finite value at row {row_idx}, column {col_idx}"
                )));
            }
            values.push(*value);
        }
        out.push(values);
    }
    Ok(out)
}

fn extract_2d_coordinates(
    coords: &Bound<'_, PyAny>,
    expected_rows: usize,
    z_policy: &str,
) -> PyResult<Vec<[f64; 2]>> {
    let z_policy = parse_z_policy(z_policy)?;
    let rows = extract_coordinate_matrix(coords, expected_rows, &[2, 3], "2D coordinates")?;
    let mut out = Vec::with_capacity(rows.len());
    for (row_idx, row) in rows.into_iter().enumerate() {
        if row.len() == 3 {
            match z_policy {
                "ignore" => {}
                "require_zero" if row[2].abs() <= 1.0e-12 => {}
                "require_zero" => {
                    return Err(PyValueError::new_err(format!(
                        "2D coordinates require zero z values but row {row_idx} has z={}",
                        row[2]
                    )));
                }
                "error" => {
                    return Err(PyValueError::new_err(
                        "2D coordinates received 3 columns while z_policy='error'",
                    ));
                }
                _ => unreachable!("z_policy is validated above"),
            }
        }
        out.push([row[0], row[1]]);
    }
    Ok(out)
}

fn extract_3d_coordinates(
    coords: &Bound<'_, PyAny>,
    expected_rows: usize,
) -> PyResult<Vec<[f64; 3]>> {
    extract_coordinate_matrix(coords, expected_rows, &[3], "3D coordinates").map(|rows| {
        rows.into_iter()
            .map(|row| [row[0], row[1], row[2]])
            .collect()
    })
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

fn rdkit_hybridization_from_name(name: &str) -> PyResult<cosmolkit_core::Hybridization> {
    match name {
        "UNSPECIFIED" => Ok(cosmolkit_core::Hybridization::Unspecified),
        "S" => Ok(cosmolkit_core::Hybridization::S),
        "SP" => Ok(cosmolkit_core::Hybridization::Sp),
        "SP2" => Ok(cosmolkit_core::Hybridization::Sp2),
        "SP3" => Ok(cosmolkit_core::Hybridization::Sp3),
        "SP2D" => Ok(cosmolkit_core::Hybridization::Sp2d),
        "SP3D" => Ok(cosmolkit_core::Hybridization::Sp3d),
        "SP3D2" => Ok(cosmolkit_core::Hybridization::Sp3d2),
        "OTHER" => Ok(cosmolkit_core::Hybridization::Other),
        _ => Err(PyValueError::new_err(format!(
            "from_rdkit unsupported atom hybridization '{name}'"
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
pub(crate) struct Molecule {
    pub(crate) inner: cosmolkit_core::Molecule,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "AlignmentAtomMap", get_all, set_all, from_py_object)]
#[derive(Clone)]
#[doc = "A probe-to-reference atom pair used by molecular alignment."]
struct PyAlignmentAtomMap {
    probe_atom: usize,
    reference_atom: usize,
}

impl From<cosmolkit_core::AlignmentAtomMap> for PyAlignmentAtomMap {
    fn from(value: cosmolkit_core::AlignmentAtomMap) -> Self {
        Self {
            probe_atom: value.probe_atom,
            reference_atom: value.reference_atom,
        }
    }
}

impl From<&PyAlignmentAtomMap> for cosmolkit_core::AlignmentAtomMap {
    fn from(value: &PyAlignmentAtomMap) -> Self {
        Self {
            probe_atom: value.probe_atom,
            reference_atom: value.reference_atom,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyAlignmentAtomMap {
    #[new]
    fn new(probe_atom: usize, reference_atom: usize) -> Self {
        Self {
            probe_atom,
            reference_atom,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "AlignmentAtomMap(probe_atom={}, reference_atom={})",
            self.probe_atom, self.reference_atom
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "AlignmentParameters", get_all, set_all, from_py_object)]
#[derive(Clone)]
#[doc = "Parameters for an explicit or first-match molecular alignment."]
struct PyAlignmentParameters {
    probe_conformer_id: i32,
    reference_conformer_id: i32,
    atom_map: Option<Vec<PyAlignmentAtomMap>>,
    weights: Option<Vec<f64>>,
    reflect: bool,
    max_iterations: u32,
}

impl PyAlignmentParameters {
    fn core_parameters(&self) -> cosmolkit_core::AlignmentParameters {
        cosmolkit_core::AlignmentParameters {
            probe_conformer_id: self.probe_conformer_id,
            reference_conformer_id: self.reference_conformer_id,
            atom_map: self
                .atom_map
                .as_ref()
                .map(|map| map.iter().map(Into::into).collect()),
            weights: self.weights.clone(),
            reflect: self.reflect,
            max_iterations: self.max_iterations,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyAlignmentParameters {
    #[new]
    #[pyo3(signature = (probe_conformer_id=-1, reference_conformer_id=-1, atom_map=None, weights=None, reflect=false, max_iterations=50))]
    fn new(
        probe_conformer_id: i32,
        reference_conformer_id: i32,
        atom_map: Option<Vec<PyAlignmentAtomMap>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
    ) -> Self {
        Self {
            probe_conformer_id,
            reference_conformer_id,
            atom_map,
            weights,
            reflect,
            max_iterations,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "BestAlignmentParameters", get_all, set_all, from_py_object)]
#[derive(Clone)]
#[doc = "Parameters for source-compatible best molecular alignment and RMSD."]
struct PyBestAlignmentParameters {
    probe_conformer_id: i32,
    reference_conformer_id: i32,
    atom_maps: Vec<Vec<PyAlignmentAtomMap>>,
    weights: Option<Vec<f64>>,
    reflect: bool,
    max_iterations: u32,
    max_matches: i32,
    symmetrize_conjugated_terminal_groups: bool,
    ignore_hydrogens: bool,
    num_threads: i32,
}

impl PyBestAlignmentParameters {
    fn core_parameters(&self) -> cosmolkit_core::BestAlignmentParameters {
        cosmolkit_core::BestAlignmentParameters {
            probe_conformer_id: self.probe_conformer_id,
            reference_conformer_id: self.reference_conformer_id,
            atom_maps: self
                .atom_maps
                .iter()
                .map(|map| map.iter().map(Into::into).collect())
                .collect(),
            weights: self.weights.clone(),
            reflect: self.reflect,
            max_iterations: self.max_iterations,
            max_matches: self.max_matches,
            symmetrize_conjugated_terminal_groups: self.symmetrize_conjugated_terminal_groups,
            ignore_hydrogens: self.ignore_hydrogens,
            num_threads: self.num_threads,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyBestAlignmentParameters {
    #[new]
    #[pyo3(signature = (probe_conformer_id=-1, reference_conformer_id=-1, atom_maps=None, weights=None, reflect=false, max_iterations=50, max_matches=1_000_000, symmetrize_conjugated_terminal_groups=true, ignore_hydrogens=true, num_threads=1))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        probe_conformer_id: i32,
        reference_conformer_id: i32,
        atom_maps: Option<Vec<Vec<PyAlignmentAtomMap>>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
        max_matches: i32,
        symmetrize_conjugated_terminal_groups: bool,
        ignore_hydrogens: bool,
        num_threads: i32,
    ) -> Self {
        Self {
            probe_conformer_id,
            reference_conformer_id,
            atom_maps: atom_maps.unwrap_or_default(),
            weights,
            reflect,
            max_iterations,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            ignore_hydrogens,
            num_threads,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "AllConformerRmsdParameters", get_all, set_all, from_py_object)]
#[derive(Clone)]
#[doc = "Parameters accepted by source-compatible all-conformer best RMSD."]
struct PyAllConformerRmsdParameters {
    atom_maps: Vec<Vec<PyAlignmentAtomMap>>,
    weights: Option<Vec<f64>>,
    max_matches: i32,
    symmetrize_conjugated_terminal_groups: bool,
    ignore_hydrogens: bool,
    num_threads: i32,
}

impl PyAllConformerRmsdParameters {
    fn core_parameters(&self) -> cosmolkit_core::AllConformerRmsdParameters {
        cosmolkit_core::AllConformerRmsdParameters {
            atom_maps: self
                .atom_maps
                .iter()
                .map(|map| map.iter().map(Into::into).collect())
                .collect(),
            weights: self.weights.clone(),
            max_matches: self.max_matches,
            symmetrize_conjugated_terminal_groups: self.symmetrize_conjugated_terminal_groups,
            ignore_hydrogens: self.ignore_hydrogens,
            num_threads: self.num_threads,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyAllConformerRmsdParameters {
    #[new]
    #[pyo3(signature = (atom_maps=None, weights=None, max_matches=1_000_000, symmetrize_conjugated_terminal_groups=true, ignore_hydrogens=true, num_threads=1))]
    fn new(
        atom_maps: Option<Vec<Vec<PyAlignmentAtomMap>>>,
        weights: Option<Vec<f64>>,
        max_matches: i32,
        symmetrize_conjugated_terminal_groups: bool,
        ignore_hydrogens: bool,
        num_threads: i32,
    ) -> Self {
        Self {
            atom_maps: atom_maps.unwrap_or_default(),
            weights,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            ignore_hydrogens,
            num_threads,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "CoordinateRmsdParameters", get_all, set_all, from_py_object)]
#[derive(Clone)]
#[doc = "Parameters for RMSD measurement in the existing coordinate frame."]
struct PyCoordinateRmsdParameters {
    probe_conformer_id: i32,
    reference_conformer_id: i32,
    atom_maps: Vec<Vec<PyAlignmentAtomMap>>,
    weights: Option<Vec<f64>>,
    max_matches: i32,
    symmetrize_conjugated_terminal_groups: bool,
}

impl PyCoordinateRmsdParameters {
    fn core_parameters(&self) -> cosmolkit_core::CoordinateRmsdParameters {
        cosmolkit_core::CoordinateRmsdParameters {
            probe_conformer_id: self.probe_conformer_id,
            reference_conformer_id: self.reference_conformer_id,
            atom_maps: self
                .atom_maps
                .iter()
                .map(|map| map.iter().map(Into::into).collect())
                .collect(),
            weights: self.weights.clone(),
            max_matches: self.max_matches,
            symmetrize_conjugated_terminal_groups: self.symmetrize_conjugated_terminal_groups,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyCoordinateRmsdParameters {
    #[new]
    #[pyo3(signature = (probe_conformer_id=-1, reference_conformer_id=-1, atom_maps=None, weights=None, max_matches=1_000_000, symmetrize_conjugated_terminal_groups=true))]
    fn new(
        probe_conformer_id: i32,
        reference_conformer_id: i32,
        atom_maps: Option<Vec<Vec<PyAlignmentAtomMap>>>,
        weights: Option<Vec<f64>>,
        max_matches: i32,
        symmetrize_conjugated_terminal_groups: bool,
    ) -> Self {
        Self {
            probe_conformer_id,
            reference_conformer_id,
            atom_maps: atom_maps.unwrap_or_default(),
            weights,
            max_matches,
            symmetrize_conjugated_terminal_groups,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ConformerAlignmentParameters",
    get_all,
    set_all,
    from_py_object
)]
#[derive(Clone)]
#[doc = "Parameters for aligning selected or all conformers of one molecule."]
struct PyConformerAlignmentParameters {
    atom_indices: Option<Vec<usize>>,
    conformer_ids: Option<Vec<usize>>,
    weights: Option<Vec<f64>>,
    reflect: bool,
    max_iterations: u32,
}

impl PyConformerAlignmentParameters {
    fn core_parameters(&self) -> cosmolkit_core::ConformerAlignmentParameters {
        cosmolkit_core::ConformerAlignmentParameters {
            atom_indices: self.atom_indices.clone(),
            conformer_ids: self.conformer_ids.clone(),
            weights: self.weights.clone(),
            reflect: self.reflect,
            max_iterations: self.max_iterations,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyConformerAlignmentParameters {
    #[new]
    #[pyo3(signature = (atom_indices=None, conformer_ids=None, weights=None, reflect=false, max_iterations=50))]
    fn new(
        atom_indices: Option<Vec<usize>>,
        conformer_ids: Option<Vec<usize>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
    ) -> Self {
        Self {
            atom_indices,
            conformer_ids,
            weights,
            reflect,
            max_iterations,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "StereoisomerOptions", skip_from_py_object)]
struct PyStereoisomerOptions {
    try_embedding: bool,
    only_unassigned: bool,
    only_stereo_groups: bool,
    max_isomers: usize,
    rand: Option<Py<PyAny>>,
    unique: bool,
}

impl Default for PyStereoisomerOptions {
    fn default() -> Self {
        let options = cosmolkit_core::StereoisomerOptions::default();
        Self {
            try_embedding: options.try_embedding,
            only_unassigned: options.only_unassigned,
            only_stereo_groups: options.only_stereo_groups,
            max_isomers: options.max_isomers,
            rand: None,
            unique: options.unique,
        }
    }
}

impl PyStereoisomerOptions {
    fn core_options(&self) -> cosmolkit_core::StereoisomerOptions {
        cosmolkit_core::StereoisomerOptions {
            try_embedding: self.try_embedding,
            only_unassigned: self.only_unassigned,
            only_stereo_groups: self.only_stereo_groups,
            max_isomers: self.max_isomers,
            random_seed: None,
            unique: self.unique,
        }
    }

    fn random_source(&self, py: Python<'_>) -> PyResult<Option<Py<PyAny>>> {
        let Some(rand) = &self.rand else {
            return Ok(None);
        };
        let random_class = py.import("random")?.getattr("Random")?;
        if rand.bind(py).is_instance(&random_class)? {
            Ok(Some(rand.clone_ref(py)))
        } else {
            Ok(Some(random_class.call1((rand.bind(py),))?.unbind()))
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyStereoisomerOptions {
    #[new]
    #[pyo3(signature = (try_embedding=false, only_unassigned=true, max_isomers=1024, rand=None, unique=true, only_stereo_groups=false))]
    fn new(
        try_embedding: bool,
        only_unassigned: bool,
        max_isomers: usize,
        rand: Option<Py<PyAny>>,
        unique: bool,
        only_stereo_groups: bool,
    ) -> Self {
        Self {
            try_embedding,
            only_unassigned,
            only_stereo_groups,
            max_isomers,
            rand,
            unique,
        }
    }

    #[getter]
    fn try_embedding(&self) -> bool {
        self.try_embedding
    }

    #[setter]
    fn set_try_embedding(&mut self, value: bool) {
        self.try_embedding = value;
    }

    #[getter]
    fn only_unassigned(&self) -> bool {
        self.only_unassigned
    }

    #[setter]
    fn set_only_unassigned(&mut self, value: bool) {
        self.only_unassigned = value;
    }

    #[getter]
    fn only_stereo_groups(&self) -> bool {
        self.only_stereo_groups
    }

    #[setter]
    fn set_only_stereo_groups(&mut self, value: bool) {
        self.only_stereo_groups = value;
    }

    #[getter]
    fn max_isomers(&self) -> usize {
        self.max_isomers
    }

    #[setter]
    fn set_max_isomers(&mut self, value: usize) {
        self.max_isomers = value;
    }

    #[getter]
    fn rand(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.rand.as_ref().map(|rand| rand.clone_ref(py))
    }

    #[setter]
    fn set_rand(&mut self, value: Option<Py<PyAny>>) {
        self.rand = value;
    }

    #[getter]
    fn unique(&self) -> bool {
        self.unique
    }

    #[setter]
    fn set_unique(&mut self, value: bool) {
        self.unique = value;
    }

    fn __repr__(&self) -> String {
        format!(
            "StereoisomerOptions(try_embedding={}, only_unassigned={}, max_isomers={}, rand={}, unique={}, only_stereo_groups={})",
            self.try_embedding,
            self.only_unassigned,
            self.max_isomers,
            if self.rand.is_some() { "..." } else { "None" },
            self.unique,
            self.only_stereo_groups,
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "PotentialStereoInfo", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyPotentialStereoInfo {
    inner: cosmolkit_core::StereoInfo,
}

fn stereo_type_name(value: cosmolkit_core::StereoType) -> &'static str {
    match value {
        cosmolkit_core::StereoType::Unspecified => "unspecified",
        cosmolkit_core::StereoType::AtomTetrahedral => "atom_tetrahedral",
        cosmolkit_core::StereoType::AtomSquarePlanar => "atom_square_planar",
        cosmolkit_core::StereoType::AtomTrigonalBipyramidal => "atom_trigonal_bipyramidal",
        cosmolkit_core::StereoType::AtomOctahedral => "atom_octahedral",
        cosmolkit_core::StereoType::BondDouble => "bond_double",
        cosmolkit_core::StereoType::BondEvenCumulene => "bond_even_cumulene",
        cosmolkit_core::StereoType::BondAtropisomer => "bond_atropisomer",
    }
}

fn stereo_specified_name(value: cosmolkit_core::StereoSpecified) -> &'static str {
    match value {
        cosmolkit_core::StereoSpecified::Unspecified => "unspecified",
        cosmolkit_core::StereoSpecified::Specified => "specified",
        cosmolkit_core::StereoSpecified::Unknown => "unknown",
    }
}

fn stereo_descriptor_name(value: cosmolkit_core::StereoDescriptor) -> &'static str {
    match value {
        cosmolkit_core::StereoDescriptor::None => "none",
        cosmolkit_core::StereoDescriptor::TetrahedralClockwise => "tetrahedral_clockwise",
        cosmolkit_core::StereoDescriptor::TetrahedralCounterclockwise => {
            "tetrahedral_counterclockwise"
        }
        cosmolkit_core::StereoDescriptor::BondCis => "bond_cis",
        cosmolkit_core::StereoDescriptor::BondTrans => "bond_trans",
        cosmolkit_core::StereoDescriptor::BondAtropClockwise => "bond_atrop_clockwise",
        cosmolkit_core::StereoDescriptor::BondAtropCounterclockwise => {
            "bond_atrop_counterclockwise"
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyPotentialStereoInfo {
    #[getter]
    fn stereo_type(&self) -> &'static str {
        stereo_type_name(self.inner.stereo_type())
    }

    #[getter]
    fn specified(&self) -> &'static str {
        stereo_specified_name(self.inner.specified())
    }

    #[getter]
    fn center_kind(&self) -> &'static str {
        match self.inner.center() {
            cosmolkit_core::StereoCenter::Missing => "missing",
            cosmolkit_core::StereoCenter::Atom(_) => "atom",
            cosmolkit_core::StereoCenter::Bond(_) => "bond",
        }
    }

    #[getter]
    fn center_index(&self) -> Option<usize> {
        match self.inner.center() {
            cosmolkit_core::StereoCenter::Missing => None,
            cosmolkit_core::StereoCenter::Atom(atom) => Some(atom.index()),
            cosmolkit_core::StereoCenter::Bond(bond) => Some(bond.index()),
        }
    }

    #[getter]
    fn descriptor(&self) -> &'static str {
        stereo_descriptor_name(self.inner.descriptor())
    }

    #[getter]
    fn permutation(&self) -> u32 {
        self.inner.permutation()
    }

    #[getter]
    fn controlling_atoms(&self) -> Vec<Option<usize>> {
        self.inner
            .controlling_atoms()
            .iter()
            .map(|atom| match atom {
                cosmolkit_core::ControllingAtom::Missing => None,
                cosmolkit_core::ControllingAtom::Atom(atom) => Some(atom.index()),
            })
            .collect()
    }

    fn __repr__(&self) -> String {
        let center_index = self
            .center_index()
            .map_or_else(|| "None".to_owned(), |index| index.to_string());
        format!(
            "PotentialStereoInfo(stereo_type='{}', specified='{}', center_kind='{}', center_index={}, descriptor='{}', permutation={})",
            self.stereo_type(),
            self.specified(),
            self.center_kind(),
            center_index,
            self.descriptor(),
            self.permutation(),
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "PotentialStereoAnalysis", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyPotentialStereoAnalysis {
    molecule: cosmolkit_core::Molecule,
    stereo_info: Vec<cosmolkit_core::StereoInfo>,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyPotentialStereoAnalysis {
    #[getter]
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.molecule.clone(),
        }
    }

    #[getter]
    fn stereo_info(&self) -> Vec<PyPotentialStereoInfo> {
        self.stereo_info
            .iter()
            .cloned()
            .map(|inner| PyPotentialStereoInfo { inner })
            .collect()
    }

    fn __len__(&self) -> usize {
        self.stereo_info.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "PotentialStereoAnalysis(records={})",
            self.stereo_info.len()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "StereoisomerIterator", skip_from_py_object)]
struct PyStereoisomerIterator {
    inner: cosmolkit_core::StereoisomerIterator,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyStereoisomerIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    #[gen_stub(override_return_type(type_repr = "Molecule"))]
    fn __next__(&mut self) -> PyResult<Option<Molecule>> {
        self.inner
            .next()
            .transpose()
            .map(|molecule| molecule.map(|inner| Molecule { inner }))
            .map_err(enumeration_pyerr)
    }

    #[getter]
    fn yielded_count(&self) -> usize {
        self.inner.yielded_count()
    }

    fn __repr__(&self) -> String {
        format!(
            "StereoisomerIterator(yielded_count={})",
            self.inner.yielded_count()
        )
    }
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
    inner: Arc<cosmolkit_core::Protein>,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
struct ProteinChain {
    inner: Arc<cosmolkit_core::Protein>,
    index: usize,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
struct ProteinResidue {
    inner: Arc<cosmolkit_core::Protein>,
    index: usize,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "ElementInfo", skip_from_py_object)]
#[derive(Clone)]
#[doc = r#"
Source-aligned periodic-table metadata used by COSMolKit.

The ``valences()`` values preserve the source periodic-table sentinels and are
not an exhaustive oxidation-state table.
"#]
struct PyElementInfo {
    inner: cosmolkit_core::ElementInfo,
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
    inner: Arc<cosmolkit_core::Protein>,
    index: usize,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
#[doc = r#"
A complete biomolecular structural value.

``BioStructure`` retains all modeled models, chains, residues, atoms, entities,
ligands, waters, nucleic acids, assemblies, and crystallographic metadata.
Use ``Protein`` only when an amino-acid-only projection is intended, and use
``Molecule`` when the desired result is a cheminformatics graph.
"#]
struct BioStructure {
    inner: Arc<cosmolkit_core::BioStructure>,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "MmcifOutputGroups", get_all, set_all, from_py_object)]
#[derive(Clone)]
#[doc = "Category switches for BioStructure mmCIF output; all categories are enabled by default."]
struct PyMmcifOutputGroups {
    atoms: bool,
    block_name: bool,
    entry: bool,
    database_status: bool,
    author: bool,
    cell: bool,
    symmetry: bool,
    entity: bool,
    entity_poly: bool,
    struct_ref: bool,
    chem_comp: bool,
    exptl: bool,
    diffrn: bool,
    reflns: bool,
    refine: bool,
    title_keywords: bool,
    ncs: bool,
    struct_asym: bool,
    origx: bool,
    struct_conf: bool,
    struct_sheet: bool,
    struct_biol: bool,
    assembly: bool,
    conn: bool,
    cis: bool,
    modres: bool,
    scale: bool,
    atom_type: bool,
    entity_poly_seq: bool,
    tls: bool,
    software: bool,
    group_pdb: bool,
    auth_all: bool,
}

impl From<cosmolkit_core::io::bio::MmcifOutputGroups> for PyMmcifOutputGroups {
    fn from(groups: cosmolkit_core::io::bio::MmcifOutputGroups) -> Self {
        Self {
            atoms: groups.atoms,
            block_name: groups.block_name,
            entry: groups.entry,
            database_status: groups.database_status,
            author: groups.author,
            cell: groups.cell,
            symmetry: groups.symmetry,
            entity: groups.entity,
            entity_poly: groups.entity_poly,
            struct_ref: groups.struct_ref,
            chem_comp: groups.chem_comp,
            exptl: groups.exptl,
            diffrn: groups.diffrn,
            reflns: groups.reflns,
            refine: groups.refine,
            title_keywords: groups.title_keywords,
            ncs: groups.ncs,
            struct_asym: groups.struct_asym,
            origx: groups.origx,
            struct_conf: groups.struct_conf,
            struct_sheet: groups.struct_sheet,
            struct_biol: groups.struct_biol,
            assembly: groups.assembly,
            conn: groups.conn,
            cis: groups.cis,
            modres: groups.modres,
            scale: groups.scale,
            atom_type: groups.atom_type,
            entity_poly_seq: groups.entity_poly_seq,
            tls: groups.tls,
            software: groups.software,
            group_pdb: groups.group_pdb,
            auth_all: groups.auth_all,
        }
    }
}

impl From<&PyMmcifOutputGroups> for cosmolkit_core::io::bio::MmcifOutputGroups {
    fn from(groups: &PyMmcifOutputGroups) -> Self {
        Self {
            atoms: groups.atoms,
            block_name: groups.block_name,
            entry: groups.entry,
            database_status: groups.database_status,
            author: groups.author,
            cell: groups.cell,
            symmetry: groups.symmetry,
            entity: groups.entity,
            entity_poly: groups.entity_poly,
            struct_ref: groups.struct_ref,
            chem_comp: groups.chem_comp,
            exptl: groups.exptl,
            diffrn: groups.diffrn,
            reflns: groups.reflns,
            refine: groups.refine,
            title_keywords: groups.title_keywords,
            ncs: groups.ncs,
            struct_asym: groups.struct_asym,
            origx: groups.origx,
            struct_conf: groups.struct_conf,
            struct_sheet: groups.struct_sheet,
            struct_biol: groups.struct_biol,
            assembly: groups.assembly,
            conn: groups.conn,
            cis: groups.cis,
            modres: groups.modres,
            scale: groups.scale,
            atom_type: groups.atom_type,
            entity_poly_seq: groups.entity_poly_seq,
            tls: groups.tls,
            software: groups.software,
            group_pdb: groups.group_pdb,
            auth_all: groups.auth_all,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyMmcifOutputGroups {
    #[new]
    #[pyo3(signature = (all=true))]
    fn new(all: bool) -> Self {
        cosmolkit_core::io::bio::MmcifOutputGroups::all(all).into()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "MmcifWriteOptions", get_all, set_all, from_py_object)]
#[derive(Clone)]
#[doc = "Options for canonical Gemmi-aligned BioStructure mmCIF serialization."]
struct PyMmcifWriteOptions {
    groups: PyMmcifOutputGroups,
    prefer_pairs: bool,
    compact: bool,
    misuse_hash: bool,
    align_pairs: u16,
    align_loops: u16,
}

impl PyMmcifWriteOptions {
    fn core_options(&self) -> cosmolkit_core::io::bio::MmcifWriteOptions {
        cosmolkit_core::io::bio::MmcifWriteOptions {
            groups: (&self.groups).into(),
            prefer_pairs: self.prefer_pairs,
            compact: self.compact,
            misuse_hash: self.misuse_hash,
            align_pairs: self.align_pairs,
            align_loops: self.align_loops,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyMmcifWriteOptions {
    #[new]
    fn new() -> Self {
        let options = cosmolkit_core::io::bio::MmcifWriteOptions::default();
        Self {
            groups: options.groups.into(),
            prefer_pairs: options.prefer_pairs,
            compact: options.compact,
            misuse_hash: options.misuse_hash,
            align_pairs: options.align_pairs,
            align_loops: options.align_loops,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
struct StructureModel {
    inner: Arc<cosmolkit_core::BioStructure>,
    index: usize,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
struct StructureChain {
    inner: Arc<cosmolkit_core::BioStructure>,
    index: usize,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
struct StructureResidue {
    inner: Arc<cosmolkit_core::BioStructure>,
    index: usize,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
struct StructureAtom {
    inner: Arc<cosmolkit_core::BioStructure>,
    index: usize,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
struct StructureEntity {
    inner: Arc<cosmolkit_core::BioStructure>,
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
    cip_code: Option<String>,
    cip_neighbor_order: Option<Vec<usize>>,
    cip_rank: Option<u32>,
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
    cip_code: Option<String>,
    cip_neighbor_order: Option<Vec<usize>>,
}

fn atom_snapshots(molecule: &cosmolkit_core::Molecule) -> Vec<Atom> {
    let assignment = cosmolkit_core::cached_valence_assignment(molecule)
        .cloned()
        .or_else(|| {
            cosmolkit_core::assign_valence(molecule, cosmolkit_core::ValenceModel::RdkitLike).ok()
        });
    let mut degrees = vec![0usize; molecule.atoms().len()];
    for bond in molecule.bonds() {
        degrees[bond.begin().index()] += 1;
        degrees[bond.end().index()] += 1;
    }
    molecule
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
            cip_code: atom.prop("_CIPCode").map(str::to_owned),
            cip_neighbor_order: atom
                .prop("_CIPNeighborOrder")
                .and_then(|value| serde_json::from_str::<Vec<u32>>(value).ok())
                .map(|values| values.into_iter().map(|value| value as usize).collect()),
            cip_rank: atom
                .prop("_CIPRank")
                .and_then(|value| value.parse::<u32>().ok()),
        })
        .collect()
}

fn bond_snapshots(molecule: &cosmolkit_core::Molecule) -> Vec<Bond> {
    molecule
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
            cip_code: bond.prop("_CIPCode").map(str::to_owned),
            cip_neighbor_order: bond
                .prop("_CIPNeighborOrder")
                .and_then(|value| serde_json::from_str::<Vec<u32>>(value).ok())
                .map(|values| values.into_iter().map(|value| value as usize).collect()),
        })
        .collect()
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
    params: cosmolkit_core::SdfReadParams,
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
    params: cosmolkit_core::SdfReadParams,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "SdfDatasetIterator", skip_from_py_object)]
struct PySdfDatasetIterator {
    dataset: cosmolkit_core::SdfDataset,
    params: cosmolkit_core::SdfReadParams,
    position: usize,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "SdfBatchIterator", skip_from_py_object)]
struct PySdfBatchIterator {
    dataset: cosmolkit_core::SdfDataset,
    params: cosmolkit_core::SdfReadParams,
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
    params: cosmolkit_core::SdfReadParams,
    start: usize,
    end: usize,
    errors: BatchErrorMode,
    progress_bar: Option<&cosmolkit_core::BatchProgressBar>,
) -> Result<cosmolkit_core::MoleculeBatch, cosmolkit_core::BatchValidationError> {
    let mut records = Vec::with_capacity(end.saturating_sub(start));
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
    params: cosmolkit_core::SdfReadParams,
    indices: Vec<usize>,
) -> Result<cosmolkit_core::MoleculeBatch, cosmolkit_core::BatchValidationError> {
    let mut records = Vec::with_capacity(indices.len());
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
    #[pyo3(signature = (path, index=None, build=None, coordinate_dim="auto", *, sanitize=None, remove_hs=None, strict_parsing=None))]
    fn open(
        _cls: &Bound<'_, PyType>,
        path: &str,
        index: Option<&Bound<'_, PyAny>>,
        build: Option<&str>,
        coordinate_dim: &str,
        sanitize: Option<bool>,
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
    ) -> PyResult<Self> {
        let expanded_path = expand_user_path(path)?;
        let params = make_sdf_read_params(sanitize, remove_hs, strict_parsing, coordinate_dim)?;
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
        let inner = cosmolkit_core::SdfDataset::open_with_params(expanded_path, params)
            .map_err(|err| PyValueError::new_err(format!("SdfDataset.open failed: {err}")))?;
        Ok(Self { inner, params })
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
                    .record_with_params(index, self.params)
                    .map_err(|err| {
                        PyValueError::new_err(format!("SdfDataset read failed: {err}"))
                    })?;
                Ok(Py::new(py, sdf_record_py(index, record))?.into_any())
            }
            Ok(indices) => {
                let inner = sdf_batch_from_indices(&self.inner, self.params, indices)
                    .map_err(batch_validation_pyerr)?;
                Ok(Py::new(py, MoleculeBatch { inner })?.into_any())
            }
        }
    }

    fn __iter__(&self) -> PySdfDatasetIterator {
        PySdfDatasetIterator {
            dataset: self.inner.clone(),
            params: self.params,
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
            params: self.params,
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
    #[pyo3(signature = (path, coordinate_dim="auto", *, sanitize=None, remove_hs=None, strict_parsing=None))]
    fn open(
        _cls: &Bound<'_, PyType>,
        path: &str,
        coordinate_dim: &str,
        sanitize: Option<bool>,
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
    ) -> PyResult<Self> {
        Ok(Self {
            path: expand_user_path(path)?,
            params: make_sdf_read_params(sanitize, remove_hs, strict_parsing, coordinate_dim)?,
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
            reader: cosmolkit_core::io::sdf::SdfReader::with_params(
                BufReader::new(file),
                self.params,
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
            .record_with_params(index, self.params)
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
            self.params,
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
    #[pyo3(signature = (sdf_text, errors=None, n_jobs=None, coordinate_dim="auto", *, sanitize=None, remove_hs=None, strict_parsing=None))]
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
        sanitize: Option<bool>,
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
    ) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let params = make_sdf_read_params(sanitize, remove_hs, strict_parsing, coordinate_dim)?;
        let inner =
            cosmolkit_core::MoleculeBatch::read_sdf_records_from_str_with_params_and_options(
                sdf_text,
                params,
                mode,
                validate_n_jobs(n_jobs)?,
                None,
            )
            .map_err(batch_validation_pyerr)?;
        Ok(Self { inner })
    }

    #[classmethod]
    #[pyo3(signature = (path, errors=None, n_jobs=None, progress_bar=false, coordinate_dim="auto", *, sanitize=None, remove_hs=None, strict_parsing=None))]
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
        sanitize: Option<bool>,
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
    ) -> PyResult<Self> {
        let mode = parse_batch_error_mode(errors)?;
        let params = make_sdf_read_params(sanitize, remove_hs, strict_parsing, coordinate_dim)?;
        let expanded_path = expand_user_path(path)?;
        if progress_bar {
            let dataset = cosmolkit_core::SdfDataset::open_with_params(&expanded_path, params)
                .map_err(|err| PyValueError::new_err(format!("read_sdf index failed: {err}")))?;
            let inner = cosmolkit_core::MoleculeBatch::read_sdf_dataset_with_params_and_options(
                &dataset,
                params,
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
        let inner =
            cosmolkit_core::MoleculeBatch::read_sdf_records_from_reader_with_params_and_options(
                BufReader::new(file),
                params,
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

    #[pyo3(signature = (n_bits=2048, tautomeric=false, n_jobs=None, progress_bar=None))]
    #[doc = r#"
Return ordered Pattern fingerprints for valid batch records.

Invalid input records remain ``None`` at their original positions. The
fingerprints use the same source-backed core and compile-once Pattern query
table as ``Molecule.pattern_fingerprint``. This ordered batch operation returns
one result per input; it is not RDKit's distinct ``MolBundle`` intersection
overload.
"#]
    fn pattern_fingerprint_list(
        &self,
        n_bits: usize,
        tautomeric: bool,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<Fingerprint>>> {
        let params = cosmolkit_core::PatternFingerprintParams { n_bits, tautomeric };
        self.inner
            .pattern_fingerprint_list_with_options(&params, validate_n_jobs(n_jobs)?, progress_bar)
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
        n_bits=2048,
        min_distance=1,
        max_distance=30,
        use_2d=true,
        include_chirality=false,
        count_simulation=true,
        count_bounds=None,
        num_bits_per_feature=1,
        from_atoms=None,
        ignore_atoms=None,
        conformer_id=-1,
        custom_atom_invariants=None,
        n_jobs=None,
        progress_bar=None
    ))]
    #[doc = r#"
Return ordered explicit-bit AtomPair fingerprints for valid batch records.

Invalid input records remain ``None`` at their original positions.
"#]
    fn fingerprint_atom_pair_list(
        &self,
        n_bits: usize,
        min_distance: u32,
        max_distance: u32,
        use_2d: bool,
        include_chirality: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        num_bits_per_feature: u32,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conformer_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<Fingerprint>>> {
        let params = make_atom_pair_fingerprint_params(
            n_bits,
            min_distance,
            max_distance,
            use_2d,
            include_chirality,
            count_simulation,
            count_bounds,
            num_bits_per_feature,
            from_atoms,
            ignore_atoms,
            conformer_id,
            custom_atom_invariants,
            false,
        );
        self.inner
            .atom_pair_fingerprint_list_with_options(
                &params,
                validate_n_jobs(n_jobs)?,
                progress_bar,
            )
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
        n_bits=2048,
        min_distance=1,
        max_distance=30,
        use_2d=true,
        include_chirality=false,
        count_simulation=true,
        count_bounds=None,
        num_bits_per_feature=1,
        from_atoms=None,
        ignore_atoms=None,
        conformer_id=-1,
        custom_atom_invariants=None,
        n_jobs=None,
        progress_bar=None
    ))]
    fn fingerprint_atom_pair_sparse_count_list(
        &self,
        n_bits: usize,
        min_distance: u32,
        max_distance: u32,
        use_2d: bool,
        include_chirality: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        num_bits_per_feature: u32,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conformer_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<PySparseCountFingerprint>>> {
        let params = make_atom_pair_fingerprint_params(
            n_bits,
            min_distance,
            max_distance,
            use_2d,
            include_chirality,
            count_simulation,
            count_bounds,
            num_bits_per_feature,
            from_atoms,
            ignore_atoms,
            conformer_id,
            custom_atom_invariants,
            false,
        );
        self.inner
            .atom_pair_sparse_count_fingerprint_list_with_options(
                &params,
                validate_n_jobs(n_jobs)?,
                progress_bar,
            )
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| PySparseCountFingerprint { inner }))
                    .collect()
            })
            .map_err(batch_validation_pyerr)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        n_bits=2048,
        min_distance=1,
        max_distance=30,
        use_2d=true,
        include_chirality=false,
        count_simulation=true,
        count_bounds=None,
        num_bits_per_feature=1,
        from_atoms=None,
        ignore_atoms=None,
        conformer_id=-1,
        custom_atom_invariants=None,
        n_jobs=None,
        progress_bar=None
    ))]
    fn fingerprint_atom_pair_count_list(
        &self,
        n_bits: usize,
        min_distance: u32,
        max_distance: u32,
        use_2d: bool,
        include_chirality: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        num_bits_per_feature: u32,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conformer_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<PySparseCountFingerprint>>> {
        let params = make_atom_pair_fingerprint_params(
            n_bits,
            min_distance,
            max_distance,
            use_2d,
            include_chirality,
            count_simulation,
            count_bounds,
            num_bits_per_feature,
            from_atoms,
            ignore_atoms,
            conformer_id,
            custom_atom_invariants,
            false,
        );
        self.inner
            .atom_pair_count_fingerprint_list_with_options(
                &params,
                validate_n_jobs(n_jobs)?,
                progress_bar,
            )
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| PySparseCountFingerprint { inner }))
                    .collect()
            })
            .map_err(batch_validation_pyerr)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        n_bits=2048,
        min_distance=1,
        max_distance=30,
        use_2d=true,
        include_chirality=false,
        count_simulation=true,
        count_bounds=None,
        num_bits_per_feature=1,
        from_atoms=None,
        ignore_atoms=None,
        conformer_id=-1,
        custom_atom_invariants=None,
        n_jobs=None,
        progress_bar=None
    ))]
    fn fingerprint_atom_pair_sparse_bits_list(
        &self,
        n_bits: usize,
        min_distance: u32,
        max_distance: u32,
        use_2d: bool,
        include_chirality: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        num_bits_per_feature: u32,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conformer_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<PySparseBitFingerprint>>> {
        let params = make_atom_pair_fingerprint_params(
            n_bits,
            min_distance,
            max_distance,
            use_2d,
            include_chirality,
            count_simulation,
            count_bounds,
            num_bits_per_feature,
            from_atoms,
            ignore_atoms,
            conformer_id,
            custom_atom_invariants,
            false,
        );
        self.inner
            .atom_pair_sparse_bit_fingerprint_list_with_options(
                &params,
                validate_n_jobs(n_jobs)?,
                progress_bar,
            )
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| PySparseBitFingerprint { inner }))
                    .collect()
            })
            .map_err(batch_validation_pyerr)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        n_bits=2048,
        min_distance=1,
        max_distance=30,
        use_2d=true,
        include_chirality=false,
        count_simulation=true,
        count_bounds=None,
        num_bits_per_feature=1,
        from_atoms=None,
        ignore_atoms=None,
        conformer_id=-1,
        custom_atom_invariants=None,
        n_jobs=None,
        progress_bar=None
    ))]
    fn fingerprint_atom_pair_with_output_list(
        &self,
        n_bits: usize,
        min_distance: u32,
        max_distance: u32,
        use_2d: bool,
        include_chirality: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        num_bits_per_feature: u32,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conformer_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<AtomPairFingerprintResult>>> {
        let params = make_atom_pair_fingerprint_params(
            n_bits,
            min_distance,
            max_distance,
            use_2d,
            include_chirality,
            count_simulation,
            count_bounds,
            num_bits_per_feature,
            from_atoms,
            ignore_atoms,
            conformer_id,
            custom_atom_invariants,
            true,
        );
        self.inner
            .atom_pair_fingerprint_with_output_list_with_options(
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

    #[allow(clippy::too_many_arguments)]
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
        layers=0xffff_ffff,
        min_path=1,
        max_path=7,
        fp_size=2048,
        atom_counts=None,
        set_only_bits=None,
        branched_paths=true,
        from_atoms=None,
        n_jobs=None,
        progress_bar=None
    ))]
    #[doc = r#"
Return ordered Layered fingerprints for valid batch records.

Invalid input records remain ``None`` at their original positions. All
fingerprints delegate to the same Rust scalar core.
"#]
    fn fingerprint_layered_list(
        &self,
        layers: u32,
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        atom_counts: Option<Vec<u32>>,
        set_only_bits: Option<&Fingerprint>,
        branched_paths: bool,
        from_atoms: Option<Vec<u32>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<Fingerprint>>> {
        let params = make_layered_fingerprint_params(
            layers,
            min_path,
            max_path,
            fp_size,
            atom_counts,
            set_only_bits,
            branched_paths,
            from_atoms,
        );
        self.inner
            .layered_fingerprint_list_with_options(&params, validate_n_jobs(n_jobs)?, progress_bar)
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
        layers=0xffff_ffff,
        min_path=1,
        max_path=7,
        fp_size=2048,
        atom_counts=None,
        set_only_bits=None,
        branched_paths=true,
        from_atoms=None,
        n_jobs=None,
        progress_bar=None
    ))]
    #[doc = r#"
Return ordered Layered fingerprints and optional updated atom counts.

Invalid input records remain ``None`` at their original positions.
"#]
    fn fingerprint_layered_with_output_list(
        &self,
        layers: u32,
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        atom_counts: Option<Vec<u32>>,
        set_only_bits: Option<&Fingerprint>,
        branched_paths: bool,
        from_atoms: Option<Vec<u32>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<LayeredFingerprintResult>>> {
        let params = make_layered_fingerprint_params(
            layers,
            min_path,
            max_path,
            fp_size,
            atom_counts,
            set_only_bits,
            branched_paths,
            from_atoms,
        );
        self.inner
            .layered_fingerprint_with_output_list_with_options(
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
    #[pyo3(signature = (inchi, *, sanitize=true, remove_hs=true))]
    #[doc = r#"
Create a molecule from an InChI string.

Returns ``None`` when the source API returns no graph or molecule
sanitization rejects the parsed graph.
"#]
    fn from_inchi(
        _cls: &Bound<'_, PyType>,
        inchi: &str,
        sanitize: bool,
        remove_hs: bool,
    ) -> PyResult<Option<Self>> {
        let output = match cosmolkit_core::mol_from_inchi(inchi.as_bytes(), sanitize, remove_hs) {
            Ok(output) => output,
            Err(error) if error.kind == cosmolkit_core::InchiErrorKind::SanitizeFailed => {
                return Ok(None);
            }
            Err(error) => return Err(inchi_pyerr(error)),
        };
        emit_inchi_diagnostics(&output.diagnostics)?;
        Ok(output.molecule.map(|inner| Self { inner }))
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
        let options = cosmolkit_core::StructureMoleculeOptions {
            sanitize,
            remove_hs,
            flavor,
            proximity_bonding,
        };
        let inner = cosmolkit_core::Molecule::from_pdb_block_with_options(text, options)
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
        let options = cosmolkit_core::StructureMoleculeOptions {
            sanitize,
            remove_hs,
            flavor,
            proximity_bonding,
        };
        let inner = cosmolkit_core::Molecule::from_mmcif_block_with_options(text, options)
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
    By default, preserve the copied RDKit graph and prepare its valence cache.
    Pass ``True`` to run full sanitization or ``False`` to retain an
    unprepared graph without a computed valence cache.

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
            let hybridization =
                rdkit_hybridization_from_name(&py_method_str(&atom, "GetHybridization")?)?;
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
            .with_radical_electrons(num_radical_electrons)
            .with_hybridization(hybridization);
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
        let inner = match sanitize {
            Some(true) => mol
                .sanitize()
                .map_err(|err| PyValueError::new_err(err.to_string()))?,
            Some(false) => mol,
            None => mol
                .with_assigned_valence()
                .map_err(|err| PyValueError::new_err(err.to_string()))?,
        };
        Ok(Self { inner })
    }

    #[classmethod]
    #[pyo3(signature = (path, sanitize=None, coordinate_dim="auto", *, remove_hs=None, strict_parsing=None))]
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
remove_hs : bool, optional
    Optional hydrogen removal flag.
strict_parsing : bool, optional
    Optional strict SDF parsing flag.
coordinate_dim : {"auto", "2d", "3d"}, optional
    Coordinate interpretation mode. ``"auto"`` preserves the molfile header.
"#]
    fn read_sdf(
        _cls: &Bound<'_, PyType>,
        path: &str,
        sanitize: Option<bool>,
        coordinate_dim: &str,
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
    ) -> PyResult<Self> {
        let expanded_path = expand_user_path(path)?;
        let coordinate_mode = parse_coordinate_mode(Some(coordinate_dim))?;
        let file = File::open(&expanded_path)
            .map_err(|e| PyValueError::new_err(format!("read_sdf open failed: {e}")))?;
        let mut reader = SdfReader::with_params(
            BufReader::new(file),
            cosmolkit_core::SdfReadParams {
                sanitize: sanitize.unwrap_or(true),
                remove_hs: remove_hs.unwrap_or(true),
                strict_parsing: strict_parsing.unwrap_or(true),
                coordinate_mode,
                ..Default::default()
            },
        );
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
    #[pyo3(signature = (path, sanitize=None, coordinate_dim="auto", *, remove_hs=None, strict_parsing=None))]
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
remove_hs : bool, optional
    Optional hydrogen removal flag.
strict_parsing : bool, optional
    Optional strict molfile parsing flag.
coordinate_dim : {"auto", "2d", "3d"}, optional
    Coordinate interpretation mode. ``"auto"`` preserves the molfile header.
"#]
    fn read_mol(
        _cls: &Bound<'_, PyType>,
        path: &str,
        sanitize: Option<bool>,
        coordinate_dim: &str,
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
    ) -> PyResult<Self> {
        let expanded_path = expand_user_path(path)?;
        let coordinate_mode = parse_coordinate_mode(Some(coordinate_dim))?;
        let record = cosmolkit_core::io::molfile::read_mol_file_with_params(
            &expanded_path,
            cosmolkit_core::SdfReadParams {
                sanitize: sanitize.unwrap_or(true),
                remove_hs: remove_hs.unwrap_or(true),
                strict_parsing: strict_parsing.unwrap_or(true),
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
    #[pyo3(signature = (mol_text, sanitize=None, coordinate_dim="auto", *, remove_hs=None, strict_parsing=None))]
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
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
    ) -> PyResult<Self> {
        let coordinate_mode = parse_coordinate_mode(Some(coordinate_dim))?;
        let record = cosmolkit_core::io::molfile::read_mol_record_from_str_with_params(
            mol_text,
            cosmolkit_core::SdfReadParams {
                sanitize: sanitize.unwrap_or(true),
                remove_hs: remove_hs.unwrap_or(true),
                strict_parsing: strict_parsing.unwrap_or(true),
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
    #[pyo3(signature = (sdf_text, sanitize=None, coordinate_dim="auto", *, remove_hs=None, strict_parsing=None))]
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
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
    ) -> PyResult<Self> {
        let coordinate_mode = parse_coordinate_mode(Some(coordinate_dim))?;
        let record = cosmolkit_core::io::sdf::read_sdf_from_str_with_params(
            sdf_text,
            cosmolkit_core::SdfReadParams {
                sanitize: sanitize.unwrap_or(true),
                remove_hs: remove_hs.unwrap_or(true),
                strict_parsing: strict_parsing.unwrap_or(true),
                coordinate_mode,
                ..Default::default()
            },
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
an error, it does not roll back and may retain partial changes, but its internal
storage remains complete. Use ``with_hydrogens()`` when failure-preserving value
semantics are required.
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

    #[pyo3(
        signature = (conf_id=-1, replace_existing_tags=true),
        text_signature = "($self, conf_id=-1, replace_existing_tags=True)"
    )]
    #[doc = r#"
Return a new molecule with atom chiral tags assigned from 3D coordinates.

The selected conformer, atom and bond ordering, coordinates, and unrelated
properties are preserved. ``conf_id=-1`` selects the default conformer.
Existing atom chiral tags are replaced unless ``replace_existing_tags`` is
false. The original molecule is left unchanged, including on error.

This stable API has exact full-state parity with RDKit 2026.03.1
``assignChiralTypesFrom3D`` across all 77 fixed oracle records. The covered
surface includes tetrahedral C/S/Se centers, enabled square-planar,
trigonal-bipyramidal, and octahedral centers, property updates, no-op paths,
and defined errors. It does not perform ``assignStereochemistryFrom3D``, 3D
double-bond direction or E/Z assignment, CIP orchestration, or
distinct-substituent validation.
"#]
    fn with_chiral_tags_from_structure(
        &self,
        conf_id: i32,
        replace_existing_tags: bool,
    ) -> PyResult<Self> {
        self.inner
            .with_chiral_tags_from_structure(conf_id, replace_existing_tags)
            .map(|inner| Self { inner })
            .map_err(chiral_tag_assignment_pyerr)
    }

    #[pyo3(
        signature = (conf_id=-1, replace_existing_tags=true),
        text_signature = "($self, conf_id=-1, replace_existing_tags=True)"
    )]
    #[doc = r#"
Assign atom chiral tags from 3D coordinates in place.

This is the in-place form of ``with_chiral_tags_from_structure()``. All public
in-place ``Molecule`` methods end with ``_``. It has the same stable,
pinned-RDKit parity scope as the value-style form. Failures are transactional
and leave the molecule unchanged.
"#]
    fn assign_chiral_tags_from_structure_(
        &mut self,
        conf_id: i32,
        replace_existing_tags: bool,
    ) -> PyResult<()> {
        self.inner
            .assign_chiral_tags_from_structure_(conf_id, replace_existing_tags)
            .map_err(chiral_tag_assignment_pyerr)
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

    #[pyo3(signature = (atoms=None, bonds=None, max_recursive_iterations=0))]
    #[doc = r#"
Return a new molecule with source-backed modern CIP labels assigned.

When both ``atoms`` and ``bonds`` are omitted or empty, the full molecule is
labeled. Once either selection is non-empty, an omitted or empty category
selects no entries, matching the pinned RDKit wrapper dispatch. Assignment is
molecular-context dependent; query descriptors from ``mol.atoms()[i]`` or
``mol.bonds()[i]`` after this call.
"#]
    fn with_cip_labels(
        &self,
        atoms: Option<Vec<usize>>,
        bonds: Option<Vec<usize>>,
        max_recursive_iterations: u32,
    ) -> PyResult<Self> {
        let mut options = cosmolkit_core::CipLabelOptions::default()
            .with_max_recursive_iterations(max_recursive_iterations);
        if let Some(indices) = atoms {
            options = options.with_atoms(indices.into_iter().map(cosmolkit_core::AtomId::new));
        }
        if let Some(indices) = bonds {
            options = options.with_bonds(indices.into_iter().map(cosmolkit_core::BondId::new));
        }
        self.inner
            .with_cip_labels_with_options(options)
            .map(|inner| Self { inner })
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    #[pyo3(signature = (atoms=None, bonds=None, max_recursive_iterations=0))]
    #[doc = r#"
Assign source-backed modern CIP labels in place.

The operation follows COSMolKit's explicit in-place policy. On an error, the
source-backed operation may retain partial source state; internal storage
remains complete and the error is raised to Python.
"#]
    fn assign_cip_labels_(
        &mut self,
        atoms: Option<Vec<usize>>,
        bonds: Option<Vec<usize>>,
        max_recursive_iterations: u32,
    ) -> PyResult<()> {
        let mut options = cosmolkit_core::CipLabelOptions::default()
            .with_max_recursive_iterations(max_recursive_iterations);
        if let Some(indices) = atoms {
            options = options.with_atoms(indices.into_iter().map(cosmolkit_core::AtomId::new));
        }
        if let Some(indices) = bonds {
            options = options.with_bonds(indices.into_iter().map(cosmolkit_core::BondId::new));
        }
        self.inner
            .assign_cip_labels_with_options_(options)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    #[doc = "Return whether the molecule has modern CIP assignment state."]
    fn cip_computed(&self) -> bool {
        self.inner.is_prop_computed("_CIPComputed") && self.inner.prop("_CIPComputed") == Some("1")
    }

    #[doc = r#"
Return read-only atom feature records.
"#]
    fn atoms(&self) -> Vec<Atom> {
        atom_snapshots(&self.inner)
    }

    #[doc = r#"
Return read-only bond feature records.
"#]
    fn bonds(&self) -> Vec<Bond> {
        bond_snapshots(&self.inner)
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

Each record is ``(center_atom_index, ordered_ligands)``. The ligand order is
the stereochemical value, not a plain adjacency listing: opposite tetrahedral
configurations can have the same ligand set but different ligand order.
Equivalent even permutations are canonicalized to one numeric representative;
odd permutations remain distinct because they encode the opposite handedness.
``None`` does not mean that a ligand is absent. It represents a hydrogen ligand
that exists chemically but is implicit in the current molecule graph and
therefore has no atom index.

Specification:
https://github.com/cosmol-studio/COSMolKit/blob/main/dev/tetrahedral_stereo.md
"#]
    fn tetrahedral_stereo(&self) -> PyResult<Vec<(usize, Vec<Option<usize>>)>> {
        to_python_tetrahedral_stereo(&self.inner)
    }

    #[pyo3(signature = (reference, params=None))]
    #[doc = r#"
Compute the transform aligning this molecule to ``reference`` without mutation.

The returned result contains the RMSD, 4x4 transform, and selected atom map.
Use ``with_alignment_to()`` or ``align_to_()`` to apply the transform.
"#]
    fn alignment_transform_to(
        &self,
        reference: &Molecule,
        params: Option<&PyAlignmentParameters>,
    ) -> PyResult<PyAlignmentResult> {
        let params = params
            .map(PyAlignmentParameters::core_parameters)
            .unwrap_or_default();
        self.inner
            .alignment_transform_to(&reference.inner, &params)
            .map(Into::into)
            .map_err(|err| PyValueError::new_err(format!("alignment_transform_to failed: {err}")))
    }

    #[pyo3(signature = (reference, params=None))]
    #[doc = "Return the best source-compatible alignment result without mutating either molecule."]
    fn best_alignment_to(
        &self,
        reference: &Molecule,
        params: Option<&PyBestAlignmentParameters>,
    ) -> PyResult<PyAlignmentResult> {
        let params = params
            .map(PyBestAlignmentParameters::core_parameters)
            .unwrap_or_default();
        self.inner
            .best_alignment_to(&reference.inner, &params)
            .map(Into::into)
            .map_err(|err| PyValueError::new_err(format!("best_alignment_to failed: {err}")))
    }

    #[pyo3(signature = (reference, params=None))]
    #[doc = "Return the best aligned RMSD without changing either molecule's coordinates."]
    fn best_rmsd_to(
        &self,
        reference: &Molecule,
        params: Option<&PyBestAlignmentParameters>,
    ) -> PyResult<f64> {
        let params = params
            .map(PyBestAlignmentParameters::core_parameters)
            .unwrap_or_default();
        self.inner
            .best_rmsd_to(&reference.inner, &params)
            .map_err(|err| PyValueError::new_err(format!("best_rmsd_to failed: {err}")))
    }

    #[pyo3(signature = (reference, params=None))]
    #[doc = r#"
Measure RMSD in the existing coordinate frame without alignment or mutation.

This method corresponds to RDKit ``CalcRMS`` semantics, including map
enumeration and optional terminal-group symmetrization.
"#]
    fn coordinate_rmsd_to(
        &self,
        reference: &Molecule,
        params: Option<&PyCoordinateRmsdParameters>,
    ) -> PyResult<f64> {
        let params = params
            .map(PyCoordinateRmsdParameters::core_parameters)
            .unwrap_or_default();
        self.inner
            .coordinate_rmsd_to(&reference.inner, &params)
            .map_err(|err| PyValueError::new_err(format!("coordinate_rmsd_to failed: {err}")))
    }

    #[pyo3(signature = (params=None))]
    #[doc = "Return best RMSD values for every ordered triangular conformer pair without mutation."]
    fn all_conformer_best_rmsds(
        &self,
        params: Option<&PyAllConformerRmsdParameters>,
    ) -> PyResult<Vec<PyConformerRmsd>> {
        let params = params
            .map(PyAllConformerRmsdParameters::core_parameters)
            .unwrap_or_default();
        self.inner
            .all_conformer_best_rmsds(&params)
            .map(|values| values.into_iter().map(Into::into).collect())
            .map_err(|err| PyValueError::new_err(format!("all_conformer_best_rmsds failed: {err}")))
    }

    #[pyo3(signature = (reference, params=None))]
    #[doc = r#"
Return a new molecule aligned to ``reference`` together with its alignment result.

The source and reference molecules remain unchanged.
"#]
    fn with_alignment_to(
        &self,
        reference: &Molecule,
        params: Option<&PyAlignmentParameters>,
    ) -> PyResult<(Molecule, PyAlignmentResult)> {
        let params = params
            .map(PyAlignmentParameters::core_parameters)
            .unwrap_or_default();
        self.inner
            .with_alignment_to(&reference.inner, &params)
            .map(|(molecule, result)| (Molecule { inner: molecule }, result.into()))
            .map_err(|err| PyValueError::new_err(format!("with_alignment_to failed: {err}")))
    }

    #[pyo3(signature = (reference, params=None))]
    #[doc = "Align this molecule to ``reference`` in place and return the applied result."]
    fn align_to_<'py>(
        mut slf: PyRefMut<'py, Self>,
        #[gen_stub(override_type(type_repr = "Molecule"))] reference: &Bound<'py, PyAny>,
        params: Option<&PyAlignmentParameters>,
    ) -> PyResult<PyAlignmentResult> {
        let params = params
            .map(PyAlignmentParameters::core_parameters)
            .unwrap_or_default();
        let reference_inner = if slf.as_ptr() == reference.as_ptr() {
            slf.inner.clone()
        } else {
            reference.extract::<PyRef<'_, Molecule>>()?.inner.clone()
        };
        slf.inner
            .align_to_(&reference_inner, &params)
            .map(Into::into)
            .map_err(|err| PyValueError::new_err(format!("align_to_ failed: {err}")))
    }

    #[pyo3(signature = (params=None))]
    #[doc = "Return a molecule with aligned conformers and the ordered source RMS report."]
    fn with_aligned_conformers(
        &self,
        params: Option<&PyConformerAlignmentParameters>,
    ) -> PyResult<(Molecule, PyConformerAlignmentReport)> {
        let params = params
            .map(PyConformerAlignmentParameters::core_parameters)
            .unwrap_or_default();
        self.inner
            .with_aligned_conformers_with_params(params)
            .map(|(molecule, report)| (Molecule { inner: molecule }, report.into()))
            .map_err(|err| PyValueError::new_err(format!("with_aligned_conformers failed: {err}")))
    }

    #[pyo3(signature = (params=None))]
    #[doc = "Align selected or all conformers in place and return the ordered source RMS report."]
    fn align_conformers_(
        &mut self,
        params: Option<&PyConformerAlignmentParameters>,
    ) -> PyResult<PyConformerAlignmentReport> {
        let params = params
            .map(PyConformerAlignmentParameters::core_parameters)
            .unwrap_or_default();
        self.inner
            .align_conformers_with_params_(params)
            .map(Into::into)
            .map_err(|err| PyValueError::new_err(format!("align_conformers_ failed: {err}")))
    }

    #[pyo3(signature = (coords=None, *, z_policy="ignore"))]
    #[doc = r#"
Return a new molecule with 2D coordinates.

When ``coords`` is omitted, COSMolKit computes 2D coordinates. When ``coords``
is provided, it must be a numeric array-like object with shape ``(num_atoms, 2)``
or ``(num_atoms, 3)``. Three-column input uses ``z_policy``:

``"ignore"``
    Use x/y columns and ignore z values.
``"require_zero"``
    Require all z values to be zero.
``"error"``
    Reject three-column input.
"#]
    fn with_2d_coordinates(
        &self,
        coords: Option<&Bound<'_, PyAny>>,
        z_policy: &str,
    ) -> PyResult<Self> {
        let out = if let Some(coords) = coords {
            let coords = extract_2d_coordinates(coords, self.inner.atoms().len(), z_policy)?;
            self.inner.with_2d_coordinate_block(coords).map_err(|err| {
                PyValueError::new_err(format!("with_2d_coordinates failed: {err}"))
            })?
        } else {
            self.inner.with_2d_coordinates().map_err(|err| {
                PyValueError::new_err(format!("with_2d_coordinates failed: {err}"))
            })?
        };
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

    #[pyo3(signature = (coords, *, z_policy="ignore"))]
    #[doc = r#"
Set 2D coordinates in place.

``coords`` must be a numeric array-like object with shape ``(num_atoms, 2)`` or
``(num_atoms, 3)``. Three-column input follows the same ``z_policy`` values as
``with_2d_coordinates(coords=...)``.
"#]
    fn set_2d_coordinates_(&mut self, coords: &Bound<'_, PyAny>, z_policy: &str) -> PyResult<()> {
        let coords = extract_2d_coordinates(coords, self.inner.atoms().len(), z_policy)?;
        self.inner
            .set_2d_coordinate_block_(coords)
            .map_err(|err| PyValueError::new_err(format!("set_2d_coordinates_ failed: {err}")))
    }

    #[pyo3(signature = (coords, conformer_index=0))]
    #[doc = r#"
Return a new molecule with an existing 3D conformer's coordinates replaced.

``coords`` must be a numeric array-like object with shape ``(num_atoms, 3)``.
The source molecule must already have a conformer at ``conformer_index``.
"#]
    fn with_3d_coordinates(
        &self,
        coords: &Bound<'_, PyAny>,
        conformer_index: usize,
    ) -> PyResult<Self> {
        let coords = extract_3d_coordinates(coords, self.inner.atoms().len())?;
        let out = self
            .inner
            .with_3d_coordinates(coords, conformer_index)
            .map_err(|err| PyValueError::new_err(format!("with_3d_coordinates failed: {err}")))?;
        Ok(Self { inner: out })
    }

    #[pyo3(signature = (coords, conformer_index=0))]
    #[doc = r#"
Replace an existing 3D conformer's coordinates in place.

``coords`` must be a numeric array-like object with shape ``(num_atoms, 3)``.
"#]
    fn set_3d_coordinates_(
        &mut self,
        coords: &Bound<'_, PyAny>,
        conformer_index: usize,
    ) -> PyResult<()> {
        let coords = extract_3d_coordinates(coords, self.inner.atoms().len())?;
        self.inner
            .set_3d_coordinates_(coords, conformer_index)
            .map_err(|err| PyValueError::new_err(format!("set_3d_coordinates_ failed: {err}")))
    }

    #[doc = r#"
Return a new molecule with all 3D conformers removed.

2D coordinates, topology, and properties are preserved.
"#]
    fn with_cleared_3d_conformers(&self) -> PyResult<Self> {
        let out = self.inner.with_cleared_3d_conformers().map_err(|err| {
            PyValueError::new_err(format!("with_cleared_3d_conformers failed: {err}"))
        })?;
        Ok(Self { inner: out })
    }

    #[doc = r#"
Remove all 3D conformers in place.

This is the in-place version of ``with_cleared_3d_conformers()``.
"#]
    fn clear_3d_conformers_(&mut self) -> PyResult<()> {
        self.inner
            .clear_3d_conformers_()
            .map_err(|err| PyValueError::new_err(format!("clear_3d_conformers_ failed: {err}")))
    }

    #[pyo3(signature = (coords, *, is_3d=true))]
    #[doc = r#"
Return a new molecule with one additional 3D conformer.

``coords`` must be a numeric array-like object with shape ``(num_atoms, 3)``.
"#]
    fn with_added_3d_conformer(&self, coords: &Bound<'_, PyAny>, is_3d: bool) -> PyResult<Self> {
        let coords = extract_3d_coordinates(coords, self.inner.atoms().len())?;
        let out = self
            .inner
            .with_added_3d_conformer(coords, is_3d)
            .map_err(|err| {
                PyValueError::new_err(format!("with_added_3d_conformer failed: {err}"))
            })?;
        Ok(Self { inner: out })
    }

    #[pyo3(signature = (coords, *, is_3d=true))]
    #[doc = r#"
Add one 3D conformer in place and return its conformer id.

``coords`` must be a numeric array-like object with shape ``(num_atoms, 3)``.
"#]
    fn add_3d_conformer_(&mut self, coords: &Bound<'_, PyAny>, is_3d: bool) -> PyResult<usize> {
        let coords = extract_3d_coordinates(coords, self.inner.atoms().len())?;
        let conformer_index = self.inner.conformers_3d().len();
        self.inner
            .add_3d_conformer_(coords, is_3d)
            .map_err(|err| PyValueError::new_err(format!("add_3d_conformer_ failed: {err}")))?;
        Ok(conformer_index)
    }

    #[pyo3(signature = (coords, *, is_3d=true))]
    #[doc = r#"
Return a new molecule with exactly one 3D conformer.

Existing 3D conformers are removed before ``coords`` is stored. This is the
COSMolKit equivalent of RDKit ``RemoveAllConformers(); AddConformer(...)`` for
manual coordinate assignment.
"#]
    fn with_only_3d_conformer(&self, coords: &Bound<'_, PyAny>, is_3d: bool) -> PyResult<Self> {
        let coords = extract_3d_coordinates(coords, self.inner.atoms().len())?;
        let out = self
            .inner
            .with_only_3d_conformer(coords, is_3d)
            .map_err(|err| {
                PyValueError::new_err(format!("with_only_3d_conformer failed: {err}"))
            })?;
        Ok(Self { inner: out })
    }

    #[pyo3(signature = (coords, *, is_3d=true))]
    #[doc = r#"
Replace all 3D conformers in place with exactly one conformer.

``coords`` must be a numeric array-like object with shape ``(num_atoms, 3)``.
"#]
    fn set_only_3d_conformer_(
        &mut self,
        coords: &Bound<'_, PyAny>,
        is_3d: bool,
    ) -> PyResult<usize> {
        let coords = extract_3d_coordinates(coords, self.inner.atoms().len())?;
        self.inner
            .set_only_3d_conformer_(coords, is_3d)
            .map_err(|err| {
                PyValueError::new_err(format!("set_only_3d_conformer_ failed: {err}"))
            })?;
        Ok(0)
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

    #[pyo3(signature = (isomeric_smarts=true, rooted_at_atom=None))]
    #[doc = r#"
Return a SMARTS string for this molecule or compiled query.

``rooted_at_atom`` selects the traversal root when provided.
"#]
    fn to_smarts(&self, isomeric_smarts: bool, rooted_at_atom: Option<usize>) -> PyResult<String> {
        // RDKit✔️✔️:   python::def("MolToSmarts",
        // RDKit✔️✔️:               (std::string(*)(const ROMol &, bool, int))RDKit::MolToSmarts,
        // RDKit✔️✔️:               (python::arg("mol"), python::arg("isomericSmiles") = true,
        // RDKit✔️✔️:                python::arg("rootedAtAtom") = -1),
        // RDKit✔️✔️:               docString.c_str());
        // Complexity review: this method builds one constant-size parameter
        // value and invokes the sole core writer without Python-side traversal.
        let params = SmilesWriteParams {
            do_isomeric_smiles: isomeric_smarts,
            rooted_at_atom,
            ..Default::default()
        };
        cosmolkit_core::mol_to_smarts(&self.inner, &params)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    #[pyo3(signature = (isomeric_smarts=true))]
    #[doc = r#"
Return a CXSMARTS string for this molecule or compiled query.
"#]
    fn to_cx_smarts(&self, isomeric_smarts: bool) -> PyResult<String> {
        // RDKit✔️✔️:   python::def("MolToCXSmarts",
        // RDKit✔️✔️:               (std::string(*)(const ROMol &, bool))RDKit::MolToCXSmarts,
        // RDKit✔️✔️:               (python::arg("mol"), python::arg("isomericSmiles") = true),
        // RDKit✔️✔️:               docString.c_str());
        // Complexity review: this method invokes the canonical SMARTS and CX
        // serializers once each and performs no conversion or reparsing.
        let params = SmilesWriteParams {
            do_isomeric_smiles: isomeric_smarts,
            ..Default::default()
        };
        cosmolkit_core::mol_to_cx_smarts(&self.inner, &params)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    #[pyo3(signature = (options = ""))]
    #[doc = r#"
Return the molecule's InChI without mutating the molecule.
"#]
    fn to_inchi(&self, options: &str) -> PyResult<String> {
        let options = (!options.is_empty()).then_some(options.as_bytes());
        let output = cosmolkit_core::mol_to_inchi(&self.inner, options).map_err(inchi_pyerr)?;
        emit_inchi_diagnostics(&output.diagnostics)?;
        inchi_output_string("mol_to_inchi", "InChI", output.inchi)
    }

    #[pyo3(signature = (options = ""))]
    #[doc = r#"
Return the molecule's InChIKey without mutating the molecule.
"#]
    fn to_inchi_key(&self, options: &str) -> PyResult<String> {
        let options = (!options.is_empty()).then_some(options.as_bytes());
        let output = cosmolkit_core::mol_to_inchi_key(&self.inner, options).map_err(inchi_pyerr)?;
        emit_inchi_diagnostics(&output.diagnostics)?;
        inchi_output_string("mol_to_inchi_key", "InChIKey", output.key)
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

    #[pyo3(signature = (n_bits=512, is_query=false, bit_flags=0xf07fff))]
    #[doc = r#"
Return the source-backed Avalon explicit bit fingerprint.

The parameters follow the pinned RDKit Python adapter: ``n_bits`` controls
the public vector size, ``is_query`` selects query-molecule preprocessing,
and ``bit_flags`` selects the Avalon feature families. The result does not
mutate the source molecule.
"#]
    fn avalon_fingerprint(
        &self,
        n_bits: u32,
        is_query: bool,
        bit_flags: u32,
    ) -> PyResult<Fingerprint> {
        let params = make_avalon_fingerprint_params(n_bits, is_query, bit_flags);
        self.inner
            .avalon_fingerprint(&params)
            .map(|inner| Fingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[pyo3(signature = (n_bits=2048, tautomeric=false))]
    #[doc = r#"
Return the source-backed Pattern fingerprint without mutating the molecule.

``tautomeric=True`` enables the pinned source's tautomer-aware structural
hashing. ``n_bits`` must be greater than zero. Query-bearing molecules follow
the source's Pattern-specific query suppression rules, and all calls reuse one
compile-once table of 13 built-in SMARTS queries.

RDKit identifies Pattern fingerprint version 1.0.0 as experimental. COSMolKit
preserves that upstream metadata while validating this ordinary-molecule
boundary exactly. The source's inert ``atomCounts`` and ``setOnlyBits``
arguments are intentionally omitted, and the distinct ``MolBundle``
intersection overload is not represented by this scalar API.
"#]
    fn pattern_fingerprint(&self, n_bits: usize, tautomeric: bool) -> PyResult<Fingerprint> {
        let params = cosmolkit_core::PatternFingerprintParams { n_bits, tautomeric };
        self.inner
            .pattern_fingerprint(&params)
            .map(|inner| Fingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        layers=0xffff_ffff,
        min_path=1,
        max_path=7,
        fp_size=2048,
        atom_counts=None,
        set_only_bits=None,
        branched_paths=true,
        from_atoms=None
    ))]
    #[doc = r#"
Return the source-backed RDKit Layered fingerprint.

``layers`` retains the source ``unsigned int`` flag value, including inactive
high bits. ``set_only_bits`` masks projected bits with another explicit bit
vector. ``from_atoms=None`` uses the unrooted source branch, while an empty
list is a present empty root selection and therefore yields no paths. The
source molecule is never mutated.
"#]
    fn fingerprint_layered(
        &self,
        layers: u32,
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        atom_counts: Option<Vec<u32>>,
        set_only_bits: Option<&Fingerprint>,
        branched_paths: bool,
        from_atoms: Option<Vec<u32>>,
    ) -> PyResult<Fingerprint> {
        let params = make_layered_fingerprint_params(
            layers,
            min_path,
            max_path,
            fp_size,
            atom_counts,
            set_only_bits,
            branched_paths,
            from_atoms,
        );
        self.inner
            .layered_fingerprint(&params)
            .map(|inner| Fingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        layers=0xffff_ffff,
        min_path=1,
        max_path=7,
        fp_size=2048,
        atom_counts=None,
        set_only_bits=None,
        branched_paths=true,
        from_atoms=None
    ))]
    #[doc = r#"
Return a source-backed Layered fingerprint and the optional updated atom counts.

When ``atom_counts`` is provided its values seed the source count vector and
the returned counts contain the source increments. Omitting it preserves the
source null-pointer branch and returns ``None`` for ``atom_counts()``.
"#]
    fn fingerprint_layered_with_output(
        &self,
        layers: u32,
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        atom_counts: Option<Vec<u32>>,
        set_only_bits: Option<&Fingerprint>,
        branched_paths: bool,
        from_atoms: Option<Vec<u32>>,
    ) -> PyResult<LayeredFingerprintResult> {
        let params = make_layered_fingerprint_params(
            layers,
            min_path,
            max_path,
            fp_size,
            atom_counts,
            set_only_bits,
            branched_paths,
            from_atoms,
        );
        self.inner
            .layered_fingerprint_with_output(&params)
            .map(Into::into)
            .map_err(fingerprint_pyerr)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        n_bits=2048,
        min_distance=1,
        max_distance=30,
        use_2d=true,
        include_chirality=false,
        count_simulation=true,
        count_bounds=None,
        num_bits_per_feature=1,
        from_atoms=None,
        ignore_atoms=None,
        conformer_id=-1,
        custom_atom_invariants=None
    ))]
    #[doc = r#"
Return the explicit-bit AtomPair fingerprint.

The implementation follows the pinned source generator for 2D or conformer
distances, chirality, count simulation, atom filters, and custom invariants.
The molecule is not mutated.
"#]
    fn fingerprint_atom_pair(
        &self,
        n_bits: usize,
        min_distance: u32,
        max_distance: u32,
        use_2d: bool,
        include_chirality: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        num_bits_per_feature: u32,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conformer_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
    ) -> PyResult<Fingerprint> {
        let params = make_atom_pair_fingerprint_params(
            n_bits,
            min_distance,
            max_distance,
            use_2d,
            include_chirality,
            count_simulation,
            count_bounds,
            num_bits_per_feature,
            from_atoms,
            ignore_atoms,
            conformer_id,
            custom_atom_invariants,
            false,
        );
        self.inner
            .atom_pair_fingerprint(&params)
            .map(|inner| Fingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        n_bits=2048,
        min_distance=1,
        max_distance=30,
        use_2d=true,
        include_chirality=false,
        count_simulation=true,
        count_bounds=None,
        num_bits_per_feature=1,
        from_atoms=None,
        ignore_atoms=None,
        conformer_id=-1,
        custom_atom_invariants=None
    ))]
    #[doc = r#"
Return the source-width sparse-count AtomPair fingerprint.
"#]
    fn fingerprint_atom_pair_sparse_count(
        &self,
        n_bits: usize,
        min_distance: u32,
        max_distance: u32,
        use_2d: bool,
        include_chirality: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        num_bits_per_feature: u32,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conformer_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
    ) -> PyResult<PySparseCountFingerprint> {
        let params = make_atom_pair_fingerprint_params(
            n_bits,
            min_distance,
            max_distance,
            use_2d,
            include_chirality,
            count_simulation,
            count_bounds,
            num_bits_per_feature,
            from_atoms,
            ignore_atoms,
            conformer_id,
            custom_atom_invariants,
            false,
        );
        let generator = cosmolkit_core::AtomPairFingerprintGenerator::new(&params)
            .map_err(fingerprint_pyerr)?;
        generator
            .sparse_count_fingerprint(&self.inner, &mut atom_pair_call_arguments(&params))
            .map(|inner| PySparseCountFingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        n_bits=2048,
        min_distance=1,
        max_distance=30,
        use_2d=true,
        include_chirality=false,
        count_simulation=true,
        count_bounds=None,
        num_bits_per_feature=1,
        from_atoms=None,
        ignore_atoms=None,
        conformer_id=-1,
        custom_atom_invariants=None
    ))]
    #[doc = r#"
Return the folded-count AtomPair fingerprint.
"#]
    fn fingerprint_atom_pair_count(
        &self,
        n_bits: usize,
        min_distance: u32,
        max_distance: u32,
        use_2d: bool,
        include_chirality: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        num_bits_per_feature: u32,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conformer_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
    ) -> PyResult<PySparseCountFingerprint> {
        let params = make_atom_pair_fingerprint_params(
            n_bits,
            min_distance,
            max_distance,
            use_2d,
            include_chirality,
            count_simulation,
            count_bounds,
            num_bits_per_feature,
            from_atoms,
            ignore_atoms,
            conformer_id,
            custom_atom_invariants,
            false,
        );
        let generator = cosmolkit_core::AtomPairFingerprintGenerator::new(&params)
            .map_err(fingerprint_pyerr)?;
        generator
            .count_fingerprint(&self.inner, &mut atom_pair_call_arguments(&params))
            .map(|inner| PySparseCountFingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        n_bits=2048,
        min_distance=1,
        max_distance=30,
        use_2d=true,
        include_chirality=false,
        count_simulation=true,
        count_bounds=None,
        num_bits_per_feature=1,
        from_atoms=None,
        ignore_atoms=None,
        conformer_id=-1,
        custom_atom_invariants=None
    ))]
    #[doc = r#"
Return the sparse-bit AtomPair fingerprint.
"#]
    fn fingerprint_atom_pair_sparse_bits(
        &self,
        n_bits: usize,
        min_distance: u32,
        max_distance: u32,
        use_2d: bool,
        include_chirality: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        num_bits_per_feature: u32,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conformer_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
    ) -> PyResult<PySparseBitFingerprint> {
        let params = make_atom_pair_fingerprint_params(
            n_bits,
            min_distance,
            max_distance,
            use_2d,
            include_chirality,
            count_simulation,
            count_bounds,
            num_bits_per_feature,
            from_atoms,
            ignore_atoms,
            conformer_id,
            custom_atom_invariants,
            false,
        );
        let generator = cosmolkit_core::AtomPairFingerprintGenerator::new(&params)
            .map_err(fingerprint_pyerr)?;
        generator
            .sparse_bit_fingerprint(&self.inner, &mut atom_pair_call_arguments(&params))
            .map(|inner| PySparseBitFingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        n_bits=2048,
        min_distance=1,
        max_distance=30,
        use_2d=true,
        include_chirality=false,
        count_simulation=true,
        count_bounds=None,
        num_bits_per_feature=1,
        from_atoms=None,
        ignore_atoms=None,
        conformer_id=-1,
        custom_atom_invariants=None
    ))]
    #[doc = r#"
Return the explicit-bit AtomPair fingerprint with exact provenance output.
"#]
    fn fingerprint_atom_pair_with_output(
        &self,
        n_bits: usize,
        min_distance: u32,
        max_distance: u32,
        use_2d: bool,
        include_chirality: bool,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        num_bits_per_feature: u32,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conformer_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
    ) -> PyResult<AtomPairFingerprintResult> {
        let params = make_atom_pair_fingerprint_params(
            n_bits,
            min_distance,
            max_distance,
            use_2d,
            include_chirality,
            count_simulation,
            count_bounds,
            num_bits_per_feature,
            from_atoms,
            ignore_atoms,
            conformer_id,
            custom_atom_invariants,
            true,
        );
        self.inner
            .atom_pair_fingerprint_with_output(&params)
            .map(Into::into)
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
Return a Morgan fingerprint for the supported RDKit bit-identical branches.

The exposed Morgan branches are checked against RDKit exact-bit parity. A
99.9% match, similarity correlation, or structurally similar hashing is not a
passing state.

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
    Passed through to the RDKit source-backed generator path.
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
Return a Morgan fingerprint with RDKit bit-identical provenance output.

The fingerprint and exposed AdditionalOutput fields are checked against RDKit
exact-bit and exact-field parity for the supported branches.
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
        fp_size=2048,
        num_bits_per_feature=2,
        use_hs=true,
        target_density=0.0,
        min_size=128,
        branched_paths=true,
        use_bond_order=true,
        atom_invariants=None,
        from_atoms=None
    ))]
    #[doc = r#"
Return the source-backed RDKit topological fingerprint.

The parameters and bit ordering follow the pinned ``RDKFingerprintMol``
boundary. Unsupported argument ranges raise ``ValueError``.
"#]
    fn topological_fingerprint(
        &self,
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        num_bits_per_feature: u32,
        use_hs: bool,
        target_density: f64,
        min_size: u32,
        branched_paths: bool,
        use_bond_order: bool,
        atom_invariants: Option<Vec<u32>>,
        from_atoms: Option<Vec<u32>>,
    ) -> PyResult<Fingerprint> {
        let params = make_topological_fingerprint_params(
            min_path,
            max_path,
            fp_size,
            num_bits_per_feature,
            use_hs,
            target_density,
            min_size,
            branched_paths,
            use_bond_order,
            atom_invariants,
            from_atoms,
        );
        self.inner
            .topological_fingerprint(&params)
            .map(|inner| Fingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[pyo3(signature = (
        min_path=1,
        max_path=7,
        fp_size=2048,
        num_bits_per_feature=2,
        use_hs=true,
        target_density=0.0,
        min_size=128,
        branched_paths=true,
        use_bond_order=true,
        atom_invariants=None,
        from_atoms=None,
        atom_bits=false,
        bit_info=false
    ))]
    #[doc = r#"
Return an RDKit topological fingerprint with typed optional provenance.

``atom_bits`` and ``bit_info`` request the corresponding source
``AdditionalOutput`` branches. Provenance bit identifiers retain the source
pre-folding values when density folding is enabled.
"#]
    fn topological_fingerprint_with_output(
        &self,
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        num_bits_per_feature: u32,
        use_hs: bool,
        target_density: f64,
        min_size: u32,
        branched_paths: bool,
        use_bond_order: bool,
        atom_invariants: Option<Vec<u32>>,
        from_atoms: Option<Vec<u32>>,
        atom_bits: bool,
        bit_info: bool,
    ) -> PyResult<TopologicalFingerprintResult> {
        let params = make_topological_fingerprint_params(
            min_path,
            max_path,
            fp_size,
            num_bits_per_feature,
            use_hs,
            target_density,
            min_size,
            branched_paths,
            use_bond_order,
            atom_invariants,
            from_atoms,
        );
        self.inner
            .topological_fingerprint_with_output(
                &params,
                cosmolkit_core::TopologicalFingerprintOutputRequest {
                    atom_bits,
                    bit_info,
                },
            )
            .map(Into::into)
            .map_err(fingerprint_pyerr)
    }

    #[pyo3(signature = (n_bits=166))]
    #[doc = r#"
Return a MACCS fingerprint using RDKit bit-identical key generation.

COSMolKit exposes the public 166-bit projection of RDKit's raw 167-bit MACCS
vector, where RDKit bit 0 is unused and raw bits 1..166 map to public bits
0..165.
"#]
    fn maccs_fingerprint(&self, n_bits: usize) -> PyResult<Fingerprint> {
        let params = cosmolkit_core::fingerprint::MaccsFingerprintParams { n_bits };
        self.inner
            .maccs_fingerprint(&params)
            .map(|inner| Fingerprint { inner })
            .map_err(fingerprint_pyerr)
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

    #[pyo3(signature = (clean=false, flag_possible=true))]
    #[doc = r#"
Analyze potential stereochemistry without mutating this molecule.

The returned analysis contains the isolated molecule state produced by the
source-defined cleanup mode and ordered typed potential-stereo records.
"#]
    fn analyze_potential_stereo(
        &self,
        clean: bool,
        flag_possible: bool,
    ) -> PyResult<PyPotentialStereoAnalysis> {
        self.inner
            .analyze_potential_stereo(cosmolkit_core::PotentialStereoOptions {
                clean,
                flag_possible,
            })
            .map(|analysis| PyPotentialStereoAnalysis {
                molecule: analysis.molecule,
                stereo_info: analysis.stereo_info,
            })
            .map_err(potential_stereo_pyerr)
    }

    #[pyo3(signature = (options=None))]
    #[doc = r#"
Return a lazy iterator over source-ordered stereoisomers.

The source molecule remains unchanged. ``options`` defaults to
``StereoisomerOptions()``. A ``random.Random`` instance or subclass supplied
through ``options.rand`` is consumed lazily through its ``getrandbits()``
method; other seed objects follow Python ``random.Random(seed)`` semantics.
"#]
    fn stereoisomers(
        &self,
        py: Python<'_>,
        options: Option<&PyStereoisomerOptions>,
    ) -> PyResult<PyStereoisomerIterator> {
        let default_options;
        let options = if let Some(options) = options {
            options
        } else {
            default_options = PyStereoisomerOptions::default();
            &default_options
        };
        let core_options = options.core_options();
        let inner = if let Some(random) = options.random_source(py)? {
            cosmolkit_core::enumerate_stereoisomers_with_random_bits(
                &self.inner,
                core_options,
                move |bit_count| {
                    Python::attach(|py| {
                        random
                            .bind(py)
                            .call_method1("getrandbits", (bit_count,))
                            .and_then(|value| value.extract::<num_bigint::BigUint>())
                            .map_err(|error| error.to_string())
                    })
                },
            )
        } else {
            self.inner.stereoisomers(core_options)
        }
        .map_err(enumeration_pyerr)?;
        Ok(PyStereoisomerIterator { inner })
    }

    #[pyo3(signature = (options=None))]
    #[doc = "Return the source-defined upper-bound stereoisomer count."]
    #[gen_stub(override_return_type(type_repr = "builtins.int", imports = ("builtins")))]
    fn stereoisomer_count(
        &self,
        options: Option<&PyStereoisomerOptions>,
    ) -> PyResult<num_bigint::BigUint> {
        let default_options;
        let options = if let Some(options) = options {
            options
        } else {
            default_options = PyStereoisomerOptions::default();
            &default_options
        };
        self.inner
            .stereoisomer_count(&options.core_options())
            .map_err(enumeration_pyerr)
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

    fn __reduce_ex__<'py>(&self, py: Python<'py>, _protocol: u8) -> PyResult<Bound<'py, PyAny>> {
        let module = py.import("cosmolkit.cosmolkit")?;
        let rebuild = module.getattr("_rebuild_molecule_from_pickle")?;
        let state = molecule_pickle_state(py, &self.inner)?;
        let args = PyTuple::new(py, [state])?;
        Ok(PyTuple::new(py, [rebuild.into_any(), args.into_any()])?.into_any())
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

fn bio_write_pyerr(error: cosmolkit_core::io::bio::BioWriteError) -> PyErr {
    match error {
        cosmolkit_core::io::bio::BioWriteError::Io(error) => PyOSError::new_err(error.to_string()),
        error => PyValueError::new_err(error.to_string()),
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl BioStructure {
    #[classmethod]
    #[doc = "Read a PDB file into the complete biomolecular structural model."]
    fn from_pdb(_cls: &Bound<'_, PyType>, path: &str) -> PyResult<Self> {
        let path = expand_user_path(path)?;
        let inner = cosmolkit_core::BioStructure::from_pdb(&path)
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self {
            inner: Arc::new(inner),
        })
    }

    #[classmethod]
    #[doc = "Read PDB text into the complete biomolecular structural model."]
    fn from_pdb_str(_cls: &Bound<'_, PyType>, text: &str) -> PyResult<Self> {
        let inner = cosmolkit_core::BioStructure::from_pdb_str(text)
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self {
            inner: Arc::new(inner),
        })
    }

    #[classmethod]
    #[doc = "Read an mmCIF file into the complete biomolecular structural model."]
    fn from_mmcif(_cls: &Bound<'_, PyType>, path: &str) -> PyResult<Self> {
        let path = expand_user_path(path)?;
        let inner = cosmolkit_core::BioStructure::from_mmcif(&path)
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self {
            inner: Arc::new(inner),
        })
    }

    #[classmethod]
    #[pyo3(signature = (text, path="input.cif"))]
    #[doc = "Read mmCIF text into the complete biomolecular structural model."]
    fn from_mmcif_str(_cls: &Bound<'_, PyType>, text: &str, path: &str) -> PyResult<Self> {
        let inner = cosmolkit_core::BioStructure::from_mmcif_str(text, path)
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self {
            inner: Arc::new(inner),
        })
    }

    #[classmethod]
    #[pyo3(signature = (text, path="input"))]
    #[doc = "Read structural text after detecting PDB, mmCIF, or mmJSON format."]
    fn from_structure_str(_cls: &Bound<'_, PyType>, text: &str, path: &str) -> PyResult<Self> {
        let inner = cosmolkit_core::BioStructure::from_structure_str(text, path)
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self {
            inner: Arc::new(inner),
        })
    }

    #[pyo3(signature = (options=None))]
    #[doc = "Serialize this complete structural model as Gemmi-aligned mmCIF without mutating it."]
    fn to_mmcif(&self, options: Option<&PyMmcifWriteOptions>) -> PyResult<String> {
        let options = options.map_or_else(
            cosmolkit_core::io::bio::MmcifWriteOptions::default,
            PyMmcifWriteOptions::core_options,
        );
        self.inner
            .to_mmcif_with_options(options)
            .map_err(bio_write_pyerr)
    }

    #[pyo3(signature = (path, options=None))]
    #[doc = "Write this complete structural model as Gemmi-aligned mmCIF without mutating it."]
    fn write_mmcif(&self, path: &str, options: Option<&PyMmcifWriteOptions>) -> PyResult<()> {
        let path = expand_user_path(path)?;
        let options = options.map_or_else(
            cosmolkit_core::io::bio::MmcifWriteOptions::default,
            PyMmcifWriteOptions::core_options,
        );
        self.inner
            .write_mmcif_with_options(path, options)
            .map_err(bio_write_pyerr)
    }

    #[doc = "Return the structure name or input data-block name."]
    fn name(&self) -> String {
        self.inner.name().to_string()
    }

    #[doc = "Return the detected input format name."]
    fn input_format(&self) -> String {
        format!("{:?}", self.inner.input_format())
    }

    #[doc = "Return the number of coordinate models."]
    fn num_models(&self) -> usize {
        self.inner.num_models()
    }

    #[doc = "Return the number of chains across all models."]
    fn num_chains(&self) -> usize {
        self.inner.num_chains()
    }

    #[doc = "Return the number of residues of every modeled kind."]
    fn num_residues(&self) -> usize {
        self.inner.num_residues()
    }

    #[doc = "Return the number of atoms of every modeled kind."]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }

    #[doc = "Return the number of structural entities."]
    fn num_entities(&self) -> usize {
        self.inner.num_entities()
    }

    #[doc = "Return all coordinate models as shared read-only views."]
    fn models(&self) -> Vec<StructureModel> {
        (0..self.inner.num_models())
            .map(|index| StructureModel {
                inner: Arc::clone(&self.inner),
                index,
            })
            .collect()
    }

    #[doc = "Return all chains as shared read-only views."]
    fn chains(&self) -> Vec<StructureChain> {
        (0..self.inner.num_chains())
            .map(|index| StructureChain {
                inner: Arc::clone(&self.inner),
                index,
            })
            .collect()
    }

    #[doc = "Return all residues, including ligands, waters, and nucleic acids."]
    fn residues(&self) -> Vec<StructureResidue> {
        (0..self.inner.num_residues())
            .map(|index| StructureResidue {
                inner: Arc::clone(&self.inner),
                index,
            })
            .collect()
    }

    #[doc = "Return all atoms as shared read-only views."]
    fn atoms(&self) -> Vec<StructureAtom> {
        (0..self.inner.num_atoms())
            .map(|index| StructureAtom {
                inner: Arc::clone(&self.inner),
                index,
            })
            .collect()
    }

    #[doc = "Return all structural entities as shared read-only views."]
    fn entities(&self) -> Vec<StructureEntity> {
        (0..self.inner.num_entities())
            .map(|index| StructureEntity {
                inner: Arc::clone(&self.inner),
                index,
            })
            .collect()
    }

    #[doc = r#"
Return an amino-acid-only ``Protein`` projection.

The returned value intentionally excludes ligands, waters, and nucleic acids;
the source ``BioStructure`` remains unchanged.
"#]
    fn protein(&self) -> Protein {
        Protein {
            inner: Arc::new(cosmolkit_core::Protein::project_from_bio_structure(
                &self.inner,
            )),
        }
    }

    #[pyo3(signature = (sanitize=true, remove_hs=true, flavor=0, proximity_bonding=true))]
    #[doc = r#"
Convert the structural rows to a cheminformatics ``Molecule``.

This is an explicit, potentially lossy model conversion. It follows the same
RDKit-compatible graph construction options as ``Molecule.from_pdb_block()``.
"#]
    fn to_molecule(
        &self,
        sanitize: bool,
        remove_hs: bool,
        flavor: u32,
        proximity_bonding: bool,
    ) -> PyResult<Molecule> {
        let options = cosmolkit_core::StructureMoleculeOptions {
            sanitize,
            remove_hs,
            flavor,
            proximity_bonding,
        };
        self.inner
            .to_molecule_with_options(options)
            .map(|inner| Molecule { inner })
            .map_err(pdb_molecule_pyerr)
    }

    fn __getitem__(&self, index: isize) -> PyResult<StructureModel> {
        let index = normalize_python_index(index, self.inner.num_models(), "model")?;
        Ok(StructureModel {
            inner: Arc::clone(&self.inner),
            index,
        })
    }

    fn __len__(&self) -> usize {
        self.inner.num_models()
    }

    fn __repr__(&self) -> String {
        format!(
            "BioStructure(models={}, chains={}, residues={}, atoms={}, entities={})",
            self.inner.num_models(),
            self.inner.num_chains(),
            self.inner.num_residues(),
            self.inner.num_atoms(),
            self.inner.num_entities()
        )
    }
}

fn normalize_python_index(index: isize, len: usize, kind: &str) -> PyResult<usize> {
    let len = len as isize;
    let index = if index < 0 { len + index } else { index };
    if index < 0 || index >= len {
        return Err(PyIndexError::new_err(format!(
            "BioStructure {kind} index out of range"
        )));
    }
    Ok(index as usize)
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl StructureModel {
    fn index(&self) -> usize {
        self.index
    }

    fn source_model_number(&self) -> Option<i32> {
        self.inner.models()[self.index].source_model_number
    }

    fn chains(&self) -> Vec<StructureChain> {
        let span = self.inner.models()[self.index].chain_span;
        (span.start as usize..span.end() as usize)
            .map(|index| StructureChain {
                inner: Arc::clone(&self.inner),
                index,
            })
            .collect()
    }

    fn __len__(&self) -> usize {
        self.inner.models()[self.index].chain_span.len as usize
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl StructureChain {
    fn index(&self) -> usize {
        self.index
    }

    fn model_index(&self) -> usize {
        self.inner.chains()[self.index].model_id.index() as usize
    }

    fn entity_index(&self) -> Option<usize> {
        self.inner.chains()[self.index]
            .entity_id
            .map(|id| id.index() as usize)
    }

    fn kind(&self) -> String {
        format!("{:?}", self.inner.chains()[self.index].kind)
    }

    fn auth_chain_id(&self) -> Option<String> {
        self.inner.chains()[self.index]
            .source
            .auth_chain_id
            .map(|id| id.as_str().to_string())
    }

    fn label_asym_id(&self) -> Option<String> {
        self.inner.chains()[self.index]
            .source
            .label_asym_id
            .map(|id| id.as_str().to_string())
    }

    fn residues(&self) -> Vec<StructureResidue> {
        let span = self.inner.chains()[self.index].residue_span;
        (span.start as usize..span.end() as usize)
            .map(|index| StructureResidue {
                inner: Arc::clone(&self.inner),
                index,
            })
            .collect()
    }

    fn atoms(&self) -> Vec<StructureAtom> {
        let residue_span = self.inner.chains()[self.index].residue_span;
        let residue_range = residue_span.start as usize..residue_span.end() as usize;
        let atom_count = residue_range
            .clone()
            .map(|index| self.inner.residues()[index].atom_span.len as usize)
            .sum();
        let mut atoms = Vec::with_capacity(atom_count);
        for residue_index in residue_range {
            let atom_span = self.inner.residues()[residue_index].atom_span;
            atoms.extend(
                (atom_span.start as usize..atom_span.end() as usize).map(|index| StructureAtom {
                    inner: Arc::clone(&self.inner),
                    index,
                }),
            );
        }
        atoms
    }

    fn __len__(&self) -> usize {
        self.inner.chains()[self.index].residue_span.len as usize
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl StructureResidue {
    fn index(&self) -> usize {
        self.index
    }

    fn chain_index(&self) -> usize {
        self.inner.residues()[self.index].chain_id.index() as usize
    }

    fn name(&self) -> String {
        self.inner.residues()[self.index].name.as_str().to_string()
    }

    fn kind(&self) -> String {
        format!("{:?}", self.inner.residues()[self.index].kind)
    }

    fn entity_kind(&self) -> String {
        format!("{:?}", self.inner.residues()[self.index].entity_kind)
    }

    #[gen_stub(override_return_type(type_repr = "ResidueCode"))]
    fn code<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let code =
            cosmolkit_core::residue_code_from_name(self.inner.residues()[self.index].name.as_str());
        residue_code_enum_member(py, code)
    }

    fn info(&self) -> PyResidueInfo {
        PyResidueInfo {
            inner: cosmolkit_core::find_tabulated_residue(
                self.inner.residues()[self.index].name.as_str(),
            ),
        }
    }

    fn source_sequence_number(&self) -> Option<i32> {
        self.inner.residues()[self.index]
            .source
            .seq_id
            .map(|id| id.seq_num)
    }

    fn insertion_code(&self) -> Option<String> {
        self.inner.residues()[self.index]
            .source
            .seq_id
            .and_then(|id| id.ins_code)
            .map(|code| char::from(code).to_string())
    }

    fn atoms(&self) -> Vec<StructureAtom> {
        let span = self.inner.residues()[self.index].atom_span;
        (span.start as usize..span.end() as usize)
            .map(|index| StructureAtom {
                inner: Arc::clone(&self.inner),
                index,
            })
            .collect()
    }

    fn __len__(&self) -> usize {
        self.inner.residues()[self.index].atom_span.len as usize
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl StructureAtom {
    fn index(&self) -> usize {
        self.index
    }

    fn residue_index(&self) -> usize {
        self.inner.atoms()[self.index].residue_id.index() as usize
    }

    fn name(&self) -> String {
        String::from_utf8_lossy(&self.inner.atoms()[self.index].name.0)
            .trim()
            .to_string()
    }

    #[gen_stub(override_return_type(type_repr = "Element"))]
    fn element<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        element_enum_member(py, self.inner.atoms()[self.index].element)
    }

    fn element_symbol(&self) -> String {
        self.inner.atoms()[self.index].element.symbol().to_string()
    }

    fn position(&self) -> Option<(f64, f64, f64)> {
        self.inner
            .atom_position(cosmolkit_core::BioAtomId::new(self.index as u32))
            .map(|[x, y, z]| (x, y, z))
    }

    fn altloc(&self) -> Option<String> {
        self.inner.atoms()[self.index]
            .altloc
            .map(|value| char::from(value.0).to_string())
    }

    fn occupancy(&self) -> Option<f32> {
        self.inner.atoms()[self.index].occupancy
    }

    fn b_factor(&self) -> Option<f32> {
        self.inner.atoms()[self.index].b_iso
    }

    fn formal_charge(&self) -> Option<i8> {
        self.inner.atoms()[self.index].formal_charge
    }

    fn source_serial(&self) -> Option<i32> {
        self.inner.atoms()[self.index]
            .source
            .serial
            .map(|serial| serial.0)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl StructureEntity {
    fn index(&self) -> usize {
        self.index
    }

    fn source_id(&self) -> String {
        self.inner.entities()[self.index]
            .source
            .source_entity_id
            .clone()
    }

    fn kind(&self) -> String {
        format!("{:?}", self.inner.entities()[self.index].kind)
    }

    fn polymer_kind(&self) -> String {
        format!("{:?}", self.inner.entities()[self.index].polymer_kind)
    }

    fn sequence(&self) -> Vec<String> {
        self.inner.entities()[self.index].sequence.clone()
    }

    fn subchains(&self) -> Vec<String> {
        self.inner.entities()[self.index]
            .subchains
            .iter()
            .map(|id| id.as_str().to_string())
            .collect()
    }

    fn __len__(&self) -> usize {
        self.inner.entities()[self.index].sequence.len()
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
        Ok(Self {
            inner: Arc::new(inner),
        })
    }

    #[classmethod]
    #[doc = r#"
Read PDB text as a protein-focused structural value.

This is the in-memory counterpart to ``Protein.from_pdb()``.
"#]
    fn from_pdb_str(_cls: &Bound<'_, PyType>, text: &str) -> PyResult<Self> {
        let inner = cosmolkit_core::Protein::from_pdb_str(text)
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self {
            inner: Arc::new(inner),
        })
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
        Ok(Self {
            inner: Arc::new(inner),
        })
    }

    #[classmethod]
    #[doc = r#"
Read mmCIF text as a protein-focused structural value.

``path`` is used for format context and diagnostic messages.
"#]
    fn from_mmcif_str(_cls: &Bound<'_, PyType>, text: &str, path: Option<&str>) -> PyResult<Self> {
        let inner = cosmolkit_core::Protein::from_mmcif_str(text, path.unwrap_or("input.cif"))
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self {
            inner: Arc::new(inner),
        })
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

    #[doc = "Return the canonical one-letter amino-acid code, if table-defined."]
    fn canonical_one_letter_code(&self) -> Option<String> {
        self.inner
            .residues()
            .nth(self.index)
            .expect("valid protein residue")
            .info()
            .canonical_one_letter_code()
            .map(|code| code.to_string())
    }

    #[doc = "Return the table-defined standard parent residue code, if available."]
    #[gen_stub(override_return_type(type_repr = "ResidueCode | None"))]
    fn parent_standard_code<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.inner
            .residues()
            .nth(self.index)
            .expect("valid protein residue")
            .info()
            .parent_standard_code()
            .map(|code| residue_code_enum_member(py, code))
            .transpose()
    }

    #[doc = "Return whether this is a non-standard amino acid in the residue table."]
    fn is_modified_amino_acid(&self) -> bool {
        self.inner
            .residues()
            .nth(self.index)
            .expect("valid protein residue")
            .info()
            .is_modified_amino_acid()
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

    fn canonical_one_letter_code(&self) -> Option<String> {
        self.inner
            .canonical_one_letter_code()
            .map(|code| code.to_string())
    }

    #[gen_stub(override_return_type(type_repr = "ResidueCode | None"))]
    fn parent_standard_code<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.inner
            .parent_standard_code()
            .map(|code| residue_code_enum_member(py, code))
            .transpose()
    }

    fn is_modified_amino_acid(&self) -> bool {
        self.inner.is_modified_amino_acid()
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
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl ProteinAtom {
    #[doc = "Return the zero-based atom index."]
    fn index(&self) -> usize {
        self.index
    }

    #[doc = "Return the atom's element as ``Element``."]
    #[gen_stub(override_return_type(type_repr = "Element"))]
    fn element<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let element = self
            .inner
            .atoms()
            .nth(self.index)
            .expect("valid protein atom")
            .row()
            .element;
        element_enum_member(py, element)
    }

    #[doc = "Return the canonical element symbol."]
    fn element_symbol(&self) -> String {
        self.inner
            .atoms()
            .nth(self.index)
            .expect("valid protein atom")
            .row()
            .element
            .symbol()
            .to_string()
    }

    #[doc = "Return the atomic number."]
    fn atomic_num(&self) -> u8 {
        self.inner
            .atoms()
            .nth(self.index)
            .expect("valid protein atom")
            .row()
            .element
            .atomic_number()
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
    fn position(&self) -> Option<(f64, f64, f64)> {
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
    #[gen_stub(override_return_type(type_repr = "Element"))]
    fn element<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let atomic_number = u8::try_from(self.atomic_num)
            .map_err(|_| PyValueError::new_err("atom atomic number is outside u8 range"))?;
        let element =
            cosmolkit_core::Element::from_atomic_number(atomic_number).ok_or_else(|| {
                PyValueError::new_err(format!(
                    "atom atomic number {atomic_number} is outside the Element domain"
                ))
            })?;
        element_enum_member(py, element)
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
    fn cip_descriptor(&self) -> Option<&str> {
        self.cip_code.as_deref()
    }
    fn cip_neighbor_order(&self) -> Option<Vec<usize>> {
        self.cip_neighbor_order.clone()
    }
    fn cip_rank(&self) -> Option<u32> {
        self.cip_rank
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
impl PyElementInfo {
    #[gen_stub(override_return_type(type_repr = "Element"))]
    fn element<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        element_enum_member(py, self.inner.element)
    }

    fn symbol(&self) -> String {
        self.inner.symbol.to_string()
    }

    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number
    }

    fn period(&self) -> u8 {
        self.inner.period
    }

    fn outer_electrons(&self) -> i32 {
        self.inner.outer_electrons
    }

    fn valences(&self) -> Vec<i32> {
        self.inner.valences.to_vec()
    }

    fn atomic_weight(&self) -> f64 {
        self.inner.atomic_weight
    }

    fn __repr__(&self) -> String {
        format!(
            "ElementInfo(symbol='{}', atomic_number={}, period={})",
            self.inner.symbol, self.inner.atomic_number, self.inner.period
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
    fn cip_descriptor(&self) -> Option<&str> {
        self.cip_code.as_deref()
    }
    fn cip_neighbor_order(&self) -> Option<Vec<usize>> {
        self.cip_neighbor_order.clone()
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
#[pyclass(skip_from_py_object)]
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
        let Some(element) = cosmolkit_core::Element::from_symbol(element) else {
            return Err(PyValueError::new_err(format!(
                "unsupported element symbol '{element}'"
            )));
        };
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
#[pyclass(name = "SparseCountFingerprint", skip_from_py_object)]
#[derive(Clone)]
struct PySparseCountFingerprint {
    inner: cosmolkit_core::SparseCountFingerprint,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySparseCountFingerprint {
    fn size(&self) -> u64 {
        self.inner.size()
    }

    fn get_value(&self, bit_id: u64) -> i32 {
        self.inner.get_val(bit_id)
    }

    fn value(&self, bit_id: u64) -> i32 {
        self.inner.get_val(bit_id)
    }

    fn nonzero_elements(&self) -> BTreeMap<u64, i32> {
        self.inner.nonzero_elements().clone()
    }

    fn __len__(&self) -> PyResult<usize> {
        usize::try_from(self.inner.size())
            .map_err(|_| PyOverflowError::new_err("fingerprint size exceeds Python platform size"))
    }

    fn __repr__(&self) -> String {
        format!(
            "SparseCountFingerprint(size={}, nonzero={})",
            self.inner.size(),
            self.inner.nonzero_elements().len()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "SparseBitFingerprint", skip_from_py_object)]
#[derive(Clone)]
struct PySparseBitFingerprint {
    inner: cosmolkit_core::SparseBitFingerprint,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySparseBitFingerprint {
    fn size(&self) -> u64 {
        self.inner.size()
    }

    fn on_bits(&self) -> Vec<u64> {
        self.inner.on_bits().iter().copied().collect()
    }

    fn __len__(&self) -> PyResult<usize> {
        usize::try_from(self.inner.size())
            .map_err(|_| PyOverflowError::new_err("fingerprint size exceeds Python platform size"))
    }

    fn __repr__(&self) -> String {
        format!(
            "SparseBitFingerprint(size={}, on_bits={})",
            self.inner.size(),
            self.inner.on_bits().len()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "AdditionalOutput", skip_from_py_object)]
#[derive(Clone, Default)]
struct PyAdditionalOutput {
    inner: cosmolkit_core::AdditionalOutput,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyAdditionalOutput {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    fn allocate_atom_to_bits(&mut self) {
        self.inner.allocate_atom_to_bits();
    }

    fn collect_atom_to_bits(&mut self) {
        self.inner.allocate_atom_to_bits();
    }

    fn allocate_bit_info_map(&mut self) {
        self.inner.allocate_bit_info_map();
    }

    fn collect_bit_info_map(&mut self) {
        self.inner.allocate_bit_info_map();
    }

    fn allocate_bit_paths(&mut self) {
        self.inner.allocate_bit_paths();
    }

    fn collect_bit_paths(&mut self) {
        self.inner.allocate_bit_paths();
    }

    fn allocate_atom_counts(&mut self) {
        self.inner.allocate_atom_counts();
    }

    fn collect_atom_counts(&mut self) {
        self.inner.allocate_atom_counts();
    }

    fn allocate_atoms_per_bit(&mut self) {
        self.inner.allocate_atoms_per_bit();
    }

    fn collect_atoms_per_bit(&mut self) {
        self.inner.allocate_atoms_per_bit();
    }

    fn atom_to_bits(&self) -> Option<Vec<Vec<u64>>> {
        self.inner.atom_to_bits.clone()
    }

    fn bit_info_map(&self) -> Option<BTreeMap<u64, Vec<(u32, u32)>>> {
        self.inner.bit_info_map.clone()
    }

    fn bit_paths(&self) -> Option<BTreeMap<u64, Vec<Vec<usize>>>> {
        self.inner.bit_paths.clone()
    }

    fn atom_counts(&self) -> Option<Vec<u32>> {
        self.inner.atom_counts.clone()
    }

    fn atoms_per_bit(&self) -> Option<BTreeMap<u64, Vec<Vec<usize>>>> {
        self.inner.atoms_per_bit.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "AdditionalOutput(atom_to_bits={}, bit_info_map={}, bit_paths={}, atom_counts={}, atoms_per_bit={})",
            self.inner.atom_to_bits.is_some(),
            self.inner.bit_info_map.is_some(),
            self.inner.bit_paths.is_some(),
            self.inner.atom_counts.is_some(),
            self.inner.atoms_per_bit.is_some()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "AtomPairsParameters", frozen, skip_from_py_object)]
struct PyAtomPairsParameters;

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyAtomPairsParameters {
    #[classattr]
    const VERSION: &'static str = cosmolkit_core::ATOM_PAIRS_VERSION;

    #[classattr]
    const NUM_TYPE_BITS: u32 = cosmolkit_core::ATOM_PAIR_NUM_TYPE_BITS;

    #[classattr]
    const NUM_PI_BITS: u32 = cosmolkit_core::ATOM_PAIR_NUM_PI_BITS;

    #[classattr]
    const NUM_BRANCH_BITS: u32 = cosmolkit_core::ATOM_PAIR_NUM_BRANCH_BITS;

    #[classattr]
    const NUM_CHIRAL_BITS: u32 = cosmolkit_core::ATOM_PAIR_NUM_CHIRAL_BITS;

    #[classattr]
    const CODE_SIZE: u32 = cosmolkit_core::ATOM_PAIR_CODE_SIZE;

    #[classattr]
    const NUM_PATH_BITS: u32 = cosmolkit_core::ATOM_PAIR_NUM_PATH_BITS;

    #[classattr]
    const MAX_PATH_LENGTH: u32 = cosmolkit_core::ATOM_PAIR_MAX_PATH_LENGTH;

    #[classattr]
    const NUM_ATOM_PAIR_FINGERPRINT_BITS: u32 = cosmolkit_core::ATOM_PAIR_NUM_FINGERPRINT_BITS;

    #[classattr]
    fn atom_types() -> Vec<u32> {
        cosmolkit_core::ATOM_PAIR_ATOM_NUMBER_TYPES.to_vec()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "AtomPairAtomInvariantsGenerator", skip_from_py_object)]
#[derive(Clone)]
struct PyAtomPairAtomInvariantsGenerator {
    inner: cosmolkit_core::fingerprint::AtomPairAtomInvGenerator,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyAtomPairAtomInvariantsGenerator {
    #[new]
    #[pyo3(signature = (include_chirality=false, topological_torsion_correction=false))]
    fn new(include_chirality: bool, topological_torsion_correction: bool) -> Self {
        Self {
            inner: cosmolkit_core::fingerprint::AtomPairAtomInvGenerator::new(
                include_chirality,
                topological_torsion_correction,
            ),
        }
    }

    fn info_string(&self) -> String {
        self.inner.infoString()
    }

    fn to_json(&self) -> String {
        self.inner.toJSON()
    }

    fn __repr__(&self) -> String {
        format!("AtomPairAtomInvariantsGenerator({})", self.info_string())
    }
}

#[derive(Clone)]
struct PyTopologicalTorsionGeneratorState {
    arguments: cosmolkit_core::fingerprint::TopologicalTorsionArguments,
    atom_invariants_generator: cosmolkit_core::fingerprint::AtomPairAtomInvGenerator,
}

impl PyTopologicalTorsionGeneratorState {
    fn build(
        &self,
    ) -> Result<
        cosmolkit_core::TopologicalTorsionFingerprintGenerator,
        cosmolkit_core::FingerprintError,
    > {
        cosmolkit_core::fingerprint::getTopologicalTorsionGenerator(
            &self.arguments,
            Some(self.atom_invariants_generator.clone()),
            true,
        )
    }
}

type SharedTopologicalTorsionGenerator = Arc<Mutex<PyTopologicalTorsionGeneratorState>>;

fn lock_topological_torsion_generator(
    state: &SharedTopologicalTorsionGenerator,
) -> PyResult<MutexGuard<'_, PyTopologicalTorsionGeneratorState>> {
    state.lock().map_err(|_| {
        PyRuntimeError::new_err("Topological Torsion generator state lock was poisoned")
    })
}

fn topological_torsion_call_arguments(
    from_atoms: Option<Vec<usize>>,
    ignore_atoms: Option<Vec<usize>>,
    conf_id: i32,
    custom_atom_invariants: Option<Vec<u32>>,
    custom_bond_invariants: Option<Vec<u32>>,
    additional_output: Option<&PyAdditionalOutput>,
) -> cosmolkit_core::fingerprint::FingerprintFuncArguments {
    // RDKit source: FingerprintGeneratorWrapper.cpp lines 43-82.
    // RDKit✔️✔️: if (!py_fromAtoms.is_none()) {
    // RDKit✔️✔️:   unsigned int len = python::len(py_fromAtoms);
    // RDKit✔️✔️:   if (len) { fromAtoms.reset(new std::vector<std::uint32_t>()); ... }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (!py_ignoreAtoms.is_none()) { ... if (len) { ... } }
    // RDKit✔️✔️: if (!py_atomInvs.is_none()) { ... if (len) { ... } }
    // RDKit✔️✔️: if (!py_bondInvs.is_none()) { ... if (len) { ... } }
    cosmolkit_core::fingerprint::FingerprintFuncArguments {
        from_atoms: normalize_fingerprint_indices(from_atoms),
        ignore_atoms: normalize_fingerprint_indices(ignore_atoms),
        conf_id,
        custom_atom_invariants: normalize_fingerprint_indices(custom_atom_invariants),
        custom_bond_invariants: normalize_fingerprint_indices(custom_bond_invariants),
        additional_output: additional_output.map(|output| output.inner.clone()),
    }
}

fn restore_topological_torsion_additional_output(
    arguments: &mut cosmolkit_core::fingerprint::FingerprintFuncArguments,
    output: &mut Option<PyRefMut<'_, PyAdditionalOutput>>,
) {
    if let (Some(output), Some(updated)) = (output.as_mut(), arguments.additional_output.take()) {
        output.inner = updated;
    }
}

fn require_complete_fingerprint_bulk<T>(values: Vec<Option<T>>) -> PyResult<Vec<T>> {
    values
        .into_iter()
        .enumerate()
        .map(|(index, value)| {
            value.ok_or_else(|| {
                PyValueError::new_err(format!("fingerprint bulk result at index {index} is null"))
            })
        })
        .collect()
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "TopologicalTorsionFingerprintOptions", skip_from_py_object)]
#[derive(Clone)]
struct PyTopologicalTorsionFingerprintOptions {
    state: SharedTopologicalTorsionGenerator,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyTopologicalTorsionFingerprintOptions {
    #[getter]
    fn count_simulation(&self) -> PyResult<bool> {
        Ok(lock_topological_torsion_generator(&self.state)?
            .arguments
            .fingerprint_arguments
            .df_count_simulation)
    }

    #[setter]
    fn set_count_simulation(&self, value: bool) -> PyResult<()> {
        lock_topological_torsion_generator(&self.state)?
            .arguments
            .fingerprint_arguments
            .df_count_simulation = value;
        Ok(())
    }

    #[getter]
    fn include_chirality(&self) -> PyResult<bool> {
        Ok(lock_topological_torsion_generator(&self.state)?
            .arguments
            .fingerprint_arguments
            .df_include_chirality)
    }

    #[setter]
    fn set_include_chirality(&self, value: bool) -> PyResult<()> {
        lock_topological_torsion_generator(&self.state)?
            .arguments
            .fingerprint_arguments
            .df_include_chirality = value;
        Ok(())
    }

    #[getter]
    fn fp_size(&self) -> PyResult<u32> {
        Ok(lock_topological_torsion_generator(&self.state)?
            .arguments
            .fingerprint_arguments
            .d_fp_size)
    }

    #[setter]
    fn set_fp_size(&self, value: u32) -> PyResult<()> {
        lock_topological_torsion_generator(&self.state)?
            .arguments
            .fingerprint_arguments
            .d_fp_size = value;
        Ok(())
    }

    #[getter]
    fn num_bits_per_feature(&self) -> PyResult<u32> {
        Ok(lock_topological_torsion_generator(&self.state)?
            .arguments
            .fingerprint_arguments
            .d_num_bits_per_feature)
    }

    #[setter]
    fn set_num_bits_per_feature(&self, value: u32) -> PyResult<()> {
        lock_topological_torsion_generator(&self.state)?
            .arguments
            .fingerprint_arguments
            .d_num_bits_per_feature = value;
        Ok(())
    }

    #[getter]
    fn torsion_atom_count(&self) -> PyResult<u32> {
        Ok(lock_topological_torsion_generator(&self.state)?
            .arguments
            .d_torsion_atom_count)
    }

    #[setter]
    fn set_torsion_atom_count(&self, value: u32) -> PyResult<()> {
        lock_topological_torsion_generator(&self.state)?
            .arguments
            .d_torsion_atom_count = value;
        Ok(())
    }

    #[getter]
    fn only_shortest_paths(&self) -> PyResult<bool> {
        Ok(lock_topological_torsion_generator(&self.state)?
            .arguments
            .df_only_shortest_paths)
    }

    #[setter]
    fn set_only_shortest_paths(&self, value: bool) -> PyResult<()> {
        lock_topological_torsion_generator(&self.state)?
            .arguments
            .df_only_shortest_paths = value;
        Ok(())
    }

    #[getter]
    fn count_bounds(&self) -> PyResult<Vec<u32>> {
        Ok(lock_topological_torsion_generator(&self.state)?
            .arguments
            .fingerprint_arguments
            .d_count_bounds
            .clone())
    }

    fn set_count_bounds(&self, bounds: Vec<u32>) -> PyResult<()> {
        lock_topological_torsion_generator(&self.state)?
            .arguments
            .fingerprint_arguments
            .d_count_bounds = bounds;
        Ok(())
    }

    fn __repr__(&self) -> PyResult<String> {
        let state = lock_topological_torsion_generator(&self.state)?;
        Ok(format!(
            "TopologicalTorsionFingerprintOptions(torsion_atom_count={}, fp_size={}, count_simulation={}, include_chirality={}, only_shortest_paths={})",
            state.arguments.d_torsion_atom_count,
            state.arguments.fingerprint_arguments.d_fp_size,
            state.arguments.fingerprint_arguments.df_count_simulation,
            state.arguments.fingerprint_arguments.df_include_chirality,
            state.arguments.df_only_shortest_paths
        ))
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "TopologicalTorsionFingerprintGenerator", skip_from_py_object)]
#[derive(Clone)]
struct PyTopologicalTorsionFingerprintGenerator {
    state: SharedTopologicalTorsionGenerator,
}

impl PyTopologicalTorsionFingerprintGenerator {
    fn core_generator(&self) -> PyResult<cosmolkit_core::TopologicalTorsionFingerprintGenerator> {
        lock_topological_torsion_generator(&self.state)?
            .build()
            .map_err(fingerprint_pyerr)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyTopologicalTorsionFingerprintGenerator {
    #[pyo3(signature = (molecule, from_atoms=None, ignore_atoms=None, conf_id=-1, custom_atom_invariants=None, custom_bond_invariants=None, additional_output=None))]
    fn get_sparse_count_fingerprint(
        &self,
        molecule: &Molecule,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conf_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        mut additional_output: Option<PyRefMut<'_, PyAdditionalOutput>>,
    ) -> PyResult<PySparseCountFingerprint> {
        let mut arguments = topological_torsion_call_arguments(
            from_atoms,
            ignore_atoms,
            conf_id,
            custom_atom_invariants,
            custom_bond_invariants,
            additional_output.as_deref(),
        );
        let result = self
            .core_generator()?
            .getSparseCountFingerprint(&molecule.inner, &mut arguments);
        restore_topological_torsion_additional_output(&mut arguments, &mut additional_output);
        result
            .map(|inner| PySparseCountFingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[pyo3(signature = (molecule, from_atoms=None, ignore_atoms=None, conf_id=-1, custom_atom_invariants=None, custom_bond_invariants=None, additional_output=None))]
    fn get_sparse_fingerprint(
        &self,
        molecule: &Molecule,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conf_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        mut additional_output: Option<PyRefMut<'_, PyAdditionalOutput>>,
    ) -> PyResult<PySparseBitFingerprint> {
        let mut arguments = topological_torsion_call_arguments(
            from_atoms,
            ignore_atoms,
            conf_id,
            custom_atom_invariants,
            custom_bond_invariants,
            additional_output.as_deref(),
        );
        let result = self
            .core_generator()?
            .getSparseFingerprint(&molecule.inner, &mut arguments);
        restore_topological_torsion_additional_output(&mut arguments, &mut additional_output);
        result
            .map(|inner| PySparseBitFingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[pyo3(signature = (molecule, from_atoms=None, ignore_atoms=None, conf_id=-1, custom_atom_invariants=None, custom_bond_invariants=None, additional_output=None))]
    fn get_count_fingerprint(
        &self,
        molecule: &Molecule,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conf_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        mut additional_output: Option<PyRefMut<'_, PyAdditionalOutput>>,
    ) -> PyResult<PySparseCountFingerprint> {
        let mut arguments = topological_torsion_call_arguments(
            from_atoms,
            ignore_atoms,
            conf_id,
            custom_atom_invariants,
            custom_bond_invariants,
            additional_output.as_deref(),
        );
        let result = self
            .core_generator()?
            .getCountFingerprint(&molecule.inner, &mut arguments);
        restore_topological_torsion_additional_output(&mut arguments, &mut additional_output);
        result
            .map(|inner| PySparseCountFingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[pyo3(signature = (molecule, from_atoms=None, ignore_atoms=None, conf_id=-1, custom_atom_invariants=None, custom_bond_invariants=None, additional_output=None))]
    fn get_fingerprint(
        &self,
        molecule: &Molecule,
        from_atoms: Option<Vec<usize>>,
        ignore_atoms: Option<Vec<usize>>,
        conf_id: i32,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        mut additional_output: Option<PyRefMut<'_, PyAdditionalOutput>>,
    ) -> PyResult<Fingerprint> {
        let mut arguments = topological_torsion_call_arguments(
            from_atoms,
            ignore_atoms,
            conf_id,
            custom_atom_invariants,
            custom_bond_invariants,
            additional_output.as_deref(),
        );
        let result = self
            .core_generator()?
            .getFingerprint(&molecule.inner, &mut arguments);
        restore_topological_torsion_additional_output(&mut arguments, &mut additional_output);
        result
            .map(|inner| Fingerprint { inner })
            .map_err(fingerprint_pyerr)
    }

    #[pyo3(signature = (molecules, num_threads=1))]
    fn get_sparse_count_fingerprints(
        &self,
        molecules: Vec<Molecule>,
        num_threads: i32,
    ) -> PyResult<Vec<PySparseCountFingerprint>> {
        let refs = molecules
            .iter()
            .map(|molecule| Some(&molecule.inner))
            .collect::<Vec<_>>();
        require_complete_fingerprint_bulk(
            self.core_generator()?
                .getSparseCountFingerprints(&refs, num_threads)
                .map_err(fingerprint_pyerr)?,
        )
        .map(|values| {
            values
                .into_iter()
                .map(|inner| PySparseCountFingerprint { inner })
                .collect()
        })
    }

    #[pyo3(signature = (molecules, num_threads=1))]
    fn get_sparse_fingerprints(
        &self,
        molecules: Vec<Molecule>,
        num_threads: i32,
    ) -> PyResult<Vec<PySparseBitFingerprint>> {
        let refs = molecules
            .iter()
            .map(|molecule| Some(&molecule.inner))
            .collect::<Vec<_>>();
        require_complete_fingerprint_bulk(
            self.core_generator()?
                .getSparseFingerprints(&refs, num_threads)
                .map_err(fingerprint_pyerr)?,
        )
        .map(|values| {
            values
                .into_iter()
                .map(|inner| PySparseBitFingerprint { inner })
                .collect()
        })
    }

    #[pyo3(signature = (molecules, num_threads=1))]
    fn get_count_fingerprints(
        &self,
        molecules: Vec<Molecule>,
        num_threads: i32,
    ) -> PyResult<Vec<PySparseCountFingerprint>> {
        let refs = molecules
            .iter()
            .map(|molecule| Some(&molecule.inner))
            .collect::<Vec<_>>();
        require_complete_fingerprint_bulk(
            self.core_generator()?
                .getCountFingerprints(&refs, num_threads)
                .map_err(fingerprint_pyerr)?,
        )
        .map(|values| {
            values
                .into_iter()
                .map(|inner| PySparseCountFingerprint { inner })
                .collect()
        })
    }

    #[pyo3(signature = (molecules, num_threads=1))]
    fn get_fingerprints(
        &self,
        molecules: Vec<Molecule>,
        num_threads: i32,
    ) -> PyResult<Vec<Fingerprint>> {
        let refs = molecules
            .iter()
            .map(|molecule| Some(&molecule.inner))
            .collect::<Vec<_>>();
        require_complete_fingerprint_bulk(
            self.core_generator()?
                .getFingerprints(&refs, num_threads)
                .map_err(fingerprint_pyerr)?,
        )
        .map(|values| {
            values
                .into_iter()
                .map(|inner| Fingerprint { inner })
                .collect()
        })
    }

    fn get_options(&self) -> PyTopologicalTorsionFingerprintOptions {
        PyTopologicalTorsionFingerprintOptions {
            state: Arc::clone(&self.state),
        }
    }

    fn get_info_string(&self) -> PyResult<String> {
        // RDKit source: FingerprintGenerator.cpp lines 157-173.
        // RDKit✔️✔️: return dp_fingerprintArguments->commonArgumentsString() + separator +
        // RDKit✔️✔️:        dp_fingerprintArguments->infoString() + separator +
        // RDKit✔️✔️:        dp_atomEnvironmentGenerator->infoString() + separator +
        // RDKit✔️✔️:        (dp_atomInvariantsGenerator ? ... : ...) +
        // RDKit✔️✔️:        (dp_bondInvariantsGenerator ? ... : "No bond invariants generator");
        let state = lock_topological_torsion_generator(&self.state)?;
        Ok(format!(
            "{} --- {} --- {} --- {} --- No bond invariants generator",
            state
                .arguments
                .fingerprint_arguments
                .common_arguments_string(),
            state.arguments.infoString(),
            cosmolkit_core::fingerprint::TopologicalTorsionEnvGenerator::new().infoString(),
            state.atom_invariants_generator.infoString(),
        ))
    }

    fn to_json(&self) -> PyResult<String> {
        let generator = self.core_generator()?;
        Ok(cosmolkit_core::fingerprint::generatorToJSON(
            &cosmolkit_core::fingerprint::TypedFingerprintGenerator::TopologicalTorsion(generator),
        ))
    }

    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        topological_torsion_generator_from_json(json)
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "TopologicalTorsionFingerprintGenerator({})",
            self.get_info_string()?
        ))
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
    #[staticmethod]
    #[doc = r#"
Create an explicit bit fingerprint from its width and on-bit indexes.

Every index must be smaller than ``n_bits``.
"#]
    fn from_on_bits(n_bits: usize, on_bits: Vec<usize>) -> PyResult<Self> {
        if let Some(bit) = on_bits.iter().copied().find(|&bit| bit >= n_bits) {
            return Err(PyValueError::new_err(format!(
                "fingerprint bit {bit} is outside n_bits={n_bits}"
            )));
        }
        Ok(Self {
            inner: cosmolkit_core::Fingerprint::from_on_bits(n_bits, on_bits),
        })
    }

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
struct LayeredFingerprintResult {
    fingerprint: Fingerprint,
    atom_counts: Option<Vec<u32>>,
}

impl From<cosmolkit_core::LayeredFingerprintResult> for LayeredFingerprintResult {
    fn from(value: cosmolkit_core::LayeredFingerprintResult) -> Self {
        Self {
            fingerprint: Fingerprint {
                inner: value.fingerprint,
            },
            atom_counts: value.atom_counts,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl LayeredFingerprintResult {
    fn fingerprint(&self) -> Fingerprint {
        self.fingerprint.clone()
    }

    fn atom_counts(&self) -> Option<Vec<u32>> {
        self.atom_counts.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "LayeredFingerprintResult(n_bits={}, has_atom_counts={})",
            self.fingerprint.inner.n_bits(),
            self.atom_counts.is_some()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct AtomPairFingerprintResult {
    fingerprint: Fingerprint,
    additional_output: Option<PyAdditionalOutput>,
}

impl From<cosmolkit_core::AtomPairFingerprintOutput> for AtomPairFingerprintResult {
    fn from(value: cosmolkit_core::AtomPairFingerprintOutput) -> Self {
        Self {
            fingerprint: Fingerprint {
                inner: value.fingerprint,
            },
            additional_output: value
                .additional_output
                .map(|inner| PyAdditionalOutput { inner }),
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomPairFingerprintResult {
    fn fingerprint(&self) -> Fingerprint {
        self.fingerprint.clone()
    }

    fn additional_output(&self) -> PyResult<PyAdditionalOutput> {
        self.additional_output.clone().ok_or_else(|| {
            PyValueError::new_err(
                "AtomPair additional output was not collected for this fingerprint result",
            )
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "AtomPairFingerprintResult(n_bits={}, has_additional_output={})",
            self.fingerprint.inner.n_bits(),
            self.additional_output.is_some()
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
struct TopologicalFingerprintResult {
    fingerprint: Fingerprint,
    atom_bits: Option<Vec<Vec<u32>>>,
    bit_info: Option<BTreeMap<u32, Vec<Vec<i32>>>>,
}

impl From<cosmolkit_core::TopologicalFingerprintResult> for TopologicalFingerprintResult {
    fn from(value: cosmolkit_core::TopologicalFingerprintResult) -> Self {
        Self {
            fingerprint: Fingerprint {
                inner: value.fingerprint,
            },
            atom_bits: value.output.atom_bits,
            bit_info: value.output.bit_info,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalFingerprintResult {
    fn fingerprint(&self) -> Fingerprint {
        self.fingerprint.clone()
    }

    fn atom_bits(&self) -> PyResult<Vec<Vec<u32>>> {
        self.atom_bits.clone().ok_or_else(|| {
            PyValueError::new_err(
                "topological atom_bits output was not requested for this fingerprint result",
            )
        })
    }

    fn bit_info(&self) -> PyResult<BTreeMap<u32, Vec<Vec<i32>>>> {
        self.bit_info.clone().ok_or_else(|| {
            PyValueError::new_err(
                "topological bit_info output was not requested for this fingerprint result",
            )
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "TopologicalFingerprintResult(n_bits={}, has_atom_bits={}, has_bit_info={})",
            self.fingerprint.inner.n_bits(),
            self.atom_bits.is_some(),
            self.bit_info.is_some()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "AlignmentTransform", skip_from_py_object)]
#[derive(Clone)]
struct PyAlignmentTransform {
    matrix: [[f64; 4]; 4],
}

impl From<cosmolkit_core::AlignmentTransform> for PyAlignmentTransform {
    fn from(value: cosmolkit_core::AlignmentTransform) -> Self {
        Self {
            matrix: value.matrix,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyAlignmentTransform {
    #[doc = "Return the source-oriented 4x4 homogeneous transform matrix."]
    fn matrix(&self) -> Vec<Vec<f64>> {
        self.matrix.iter().map(|row| row.to_vec()).collect()
    }

    fn __repr__(&self) -> String {
        "AlignmentTransform(matrix=4x4)".to_string()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "AlignmentResult", skip_from_py_object)]
#[derive(Clone)]
struct PyAlignmentResult {
    inner: cosmolkit_core::AlignmentResult,
}

impl From<cosmolkit_core::AlignmentResult> for PyAlignmentResult {
    fn from(inner: cosmolkit_core::AlignmentResult) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyAlignmentResult {
    #[doc = "Return the aligned root-mean-square deviation."]
    fn rmsd(&self) -> f64 {
        self.inner.rmsd
    }

    #[doc = "Return the transform mapping probe coordinates onto the reference."]
    fn transform(&self) -> PyAlignmentTransform {
        self.inner.transform.into()
    }

    #[doc = "Return the selected probe-to-reference atom map."]
    fn atom_map(&self) -> Vec<PyAlignmentAtomMap> {
        self.inner
            .atom_map
            .iter()
            .copied()
            .map(Into::into)
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "AlignmentResult(rmsd={}, mapped_atoms={})",
            self.inner.rmsd,
            self.inner.atom_map.len()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "ConformerRmsd", skip_from_py_object)]
#[derive(Clone)]
struct PyConformerRmsd {
    inner: cosmolkit_core::ConformerRmsd,
}

impl From<cosmolkit_core::ConformerRmsd> for PyConformerRmsd {
    fn from(inner: cosmolkit_core::ConformerRmsd) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyConformerRmsd {
    fn probe_conformer_id(&self) -> usize {
        self.inner.probe_conformer_id
    }

    fn reference_conformer_id(&self) -> usize {
        self.inner.reference_conformer_id
    }

    fn rmsd(&self) -> f64 {
        self.inner.rmsd
    }

    fn __repr__(&self) -> String {
        format!(
            "ConformerRmsd(probe_conformer_id={}, reference_conformer_id={}, rmsd={})",
            self.inner.probe_conformer_id, self.inner.reference_conformer_id, self.inner.rmsd
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "ConformerAlignmentReport", skip_from_py_object)]
#[derive(Clone)]
struct PyConformerAlignmentReport {
    rmsds: Vec<f64>,
}

impl From<cosmolkit_core::ConformerAlignmentReport> for PyConformerAlignmentReport {
    fn from(value: cosmolkit_core::ConformerAlignmentReport) -> Self {
        Self { rmsds: value.rmsds }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyConformerAlignmentReport {
    #[doc = "Return RMS values in the source conformer-selection order."]
    fn rmsds(&self) -> Vec<f64> {
        self.rmsds.clone()
    }

    fn __repr__(&self) -> String {
        format!("ConformerAlignmentReport(rmsds={})", self.rmsds.len())
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
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct SubstructMatchResult {
    atom_mapping: Vec<usize>,
    bond_mapping: Vec<usize>,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object)]
#[derive(Clone)]
struct SmartsParserParams {
    inner: cosmolkit_core::smarts_parse::SmartsParseParams,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SmartsParserParams {
    // RDKit✔️✔️:   python::class_<RDKit::SmartsParserParams, boost::noncopyable>(
    // RDKit✔️✔️:       "SmartsParserParams", "Parameters controlling SMARTS Parsing")
    // RDKit✔️✔️:       .def_readwrite("debugParse", &RDKit::SmartsParserParams::debugParse,
    // RDKit✔️✔️:                      "controls the amount of debugging information produced")
    // RDKit✔️✔️:       .def_readwrite("parseName", &RDKit::SmartsParserParams::parseName,
    // RDKit✔️✔️:                      "controls whether or not the molecule name is also parsed")
    // RDKit✔️✔️:       .def_readwrite(
    // RDKit✔️✔️:           "allowCXSMILES", &RDKit::SmartsParserParams::allowCXSMILES,
    // RDKit✔️✔️:           "controls whether or not the CXSMILES extensions are parsed")
    // RDKit✔️✔️:       .def_readwrite("strictCXSMILES",
    // RDKit✔️✔️:                      &RDKit::SmartsParserParams::strictCXSMILES,
    // RDKit✔️✔️:                      "controls whether or not problems in CXSMILES parsing "
    // RDKit✔️✔️:                      "causes molecule parsing to fail")
    // RDKit✔️✔️:       .def_readwrite(
    // RDKit✔️✔️:           "mergeHs", &RDKit::SmartsParserParams::mergeHs,
    // RDKit✔️✔️:           "toggles merging H atoms in the SMARTS into neighboring atoms")
    // RDKit✔️✔️:       .def("__setattr__", &safeSetattr);
    // Complexity review: registration creates one Python type with constant
    // field descriptors. Snake-case names follow the COSMolKit Python API and
    // every field maps directly onto the sole core parameter value.
    #[new]
    fn new() -> Self {
        Self {
            inner: cosmolkit_core::smarts_parse::SmartsParseParams::default(),
        }
    }

    #[getter]
    fn debug_parse(&self) -> bool {
        self.inner.debug_parse
    }

    #[setter]
    fn set_debug_parse(&mut self, value: bool) {
        self.inner.debug_parse = value;
    }

    #[getter]
    fn parse_name(&self) -> bool {
        self.inner.parse_name
    }

    #[setter]
    fn set_parse_name(&mut self, value: bool) {
        self.inner.parse_name = value;
    }

    #[getter]
    fn allow_cxsmiles(&self) -> bool {
        self.inner.allow_cxsmiles
    }

    #[setter]
    fn set_allow_cxsmiles(&mut self, value: bool) {
        self.inner.allow_cxsmiles = value;
    }

    #[getter]
    fn strict_cxsmiles(&self) -> bool {
        self.inner.strict_cxsmiles
    }

    #[setter]
    fn set_strict_cxsmiles(&mut self, value: bool) {
        self.inner.strict_cxsmiles = value;
    }

    #[getter]
    fn merge_hs(&self) -> bool {
        self.inner.merge_hs
    }

    #[setter]
    fn set_merge_hs(&mut self, value: bool) {
        self.inner.merge_hs = value;
    }

    #[getter]
    fn skip_cleanup(&self) -> bool {
        self.inner.skip_cleanup
    }

    #[setter]
    fn set_skip_cleanup(&mut self, value: bool) {
        self.inner.skip_cleanup = value;
    }

    #[getter]
    fn replacements(&self) -> BTreeMap<String, String> {
        self.inner.replacements.clone()
    }

    #[setter]
    fn set_replacements(&mut self, value: BTreeMap<String, String>) {
        self.inner.replacements = value;
    }
}

impl From<cosmolkit_core::SubstructMatchResult> for SubstructMatchResult {
    fn from(value: cosmolkit_core::SubstructMatchResult) -> Self {
        Self {
            atom_mapping: convert_matches(value.atom_mapping),
            bond_mapping: value.bond_mapping,
        }
    }
}

fn convert_matches(matches: Vec<usize>) -> Vec<usize> {
    // RDKit✔️🔝: inline PyObject *convertMatches(const MatchVectType &matches) {
    // RDKit✔️🔝:   PyObject *res = PyTuple_New(matches.size());
    // RDKit✔️🔝:   std::for_each(matches.begin(), matches.end(), [res](const auto &pair) {
    // RDKit✔️🔝:     PyTuple_SetItem(res, pair.first, PyInt_FromLong(pair.second));
    // RDKit✔️🔝:   });
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    // COSMolKit's core match result already stores the target atom ids as a
    // dense query-indexed vector, which is exactly the source function's
    // output representation. Moving that vector preserves all ordering and
    // empty-result semantics while avoiding RDKit's O(query atoms) tuple fill.
    matches
}

fn convert_matches_to_tuple_of_pairs(matches: &[usize]) -> Vec<(usize, usize)> {
    // RDKit✔️✔️: inline PyObject *convertMatchesToTupleOfPairs(const MatchVectType &matches) {
    // RDKit✔️✔️:   PyObject *res = PyTuple_New(matches.size());
    // RDKit✔️✔️:   std::for_each(matches.begin(), matches.end(),
    // RDKit✔️✔️:                 [res, &matches](const auto &pair) {
    // RDKit✔️✔️:                   PyObject *pyPair = PyTuple_New(2);
    // RDKit✔️✔️:                   PyTuple_SetItem(pyPair, 0, PyInt_FromLong(pair.first));
    // RDKit✔️✔️:                   PyTuple_SetItem(pyPair, 1, PyInt_FromLong(pair.second));
    // RDKit✔️✔️:                   PyTuple_SetItem(res, &pair - &matches.front(), pyPair);
    // RDKit✔️✔️:                 });
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Complexity review: one O(query atoms) pass allocates one pair per atom,
    // matching the source tuple-of-pairs construction and preserving order.
    matches.iter().copied().enumerate().collect()
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
    fn atom_pairs(&self) -> Vec<(usize, usize)> {
        convert_matches_to_tuple_of_pairs(&self.atom_mapping)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn version() -> &'static str {
    cosmolkit_core::version()
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
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

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pyfunction]
#[doc = r#"
Deserialize a molecule from COSMolKit binary bytes.
"#]
fn mol_from_binary(
    #[gen_stub(override_type(type_repr = "builtins.bytes", imports = ("builtins")))] data: &[u8],
) -> PyResult<Molecule> {
    let inner = cosmolkit_core::mol_from_binary(data).map_err(pickle_pyerr)?;
    Ok(Molecule { inner })
}

#[pyfunction]
#[doc(hidden)]
fn _rebuild_molecule_from_pickle(state: &Bound<'_, PyAny>) -> PyResult<Molecule> {
    let dict = state
        .cast::<PyDict>()
        .map_err(|_| PyValueError::new_err("invalid Molecule pickle state: expected dict"))?;
    let kind: String = dict
        .get_item("kind")?
        .ok_or_else(|| PyValueError::new_err("invalid Molecule pickle state: missing kind"))?
        .extract()?;
    if kind != "cosmolkit.Molecule" {
        return Err(PyValueError::new_err(format!(
            "invalid Molecule pickle state kind: {kind}"
        )));
    }
    let schema: u16 = dict
        .get_item("pickle_schema")?
        .ok_or_else(|| {
            PyValueError::new_err("invalid Molecule pickle state: missing pickle_schema")
        })?
        .extract()?;
    if schema != 1 {
        return Err(PyValueError::new_err(format!(
            "unsupported Molecule pickle schema: {schema}"
        )));
    }
    let core_format: String = dict
        .get_item("core_format")?
        .ok_or_else(|| PyValueError::new_err("invalid Molecule pickle state: missing core_format"))?
        .extract()?;
    if core_format != "cosmolkit-molecule-archive" {
        return Err(PyValueError::new_err(format!(
            "unsupported Molecule pickle core_format: {core_format}"
        )));
    }
    let payload = dict
        .get_item("payload")?
        .ok_or_else(|| PyValueError::new_err("invalid Molecule pickle state: missing payload"))?;
    let payload = payload.cast::<PyBytes>().map_err(|_| {
        PyValueError::new_err("invalid Molecule pickle state: payload must be bytes")
    })?;
    let inner = cosmolkit_core::mol_from_binary(payload.as_bytes()).map_err(pickle_pyerr)?;
    Ok(Molecule { inner })
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[doc = r#"
Compile SMARTS text into a query-bearing ``Molecule``.
"#]
#[pyo3(signature = (smarts, merge_hs=false, replacements=None))]
fn parse_smarts(
    smarts: &str,
    merge_hs: bool,
    replacements: Option<BTreeMap<String, String>>,
) -> PyResult<Molecule> {
    // RDKit✔️✔️: ROMol *MolFromSmarts(python::object ismarts, bool mergeHs,
    // RDKit✔️✔️:                      python::dict replDict) {
    // RDKit✔️✔️:   std::map<std::string, std::string> replacements;
    // RDKit✔️✔️:   const auto items = replDict.items();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < python::len(items); ++i) {
    // RDKit✔️✔️:     const auto item = items[i];
    // RDKit✔️✔️:     replacements[python::extract<std::string>(item[0])] =
    // RDKit✔️✔️:         python::extract<std::string>(item[1]);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::string smarts = pyObjectToString(ismarts);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RWMol *newM;
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     newM = SmartsToMol(smarts, 0, mergeHs, &replacements);
    // RDKit✔️✔️:   } catch (...) {
    // RDKit✔️✔️:     newM = nullptr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return static_cast<ROMol *>(newM);
    // RDKit✔️✔️: }
    // Complexity review: dictionary extraction is linear in replacements and
    // the canonical compiler is linear in SMARTS plus graph size. Rust returns
    // a structured parse error instead of RDKit's Python nullptr convention.
    let params = cosmolkit_core::smarts_parse::SmartsParseParams {
        merge_hs,
        replacements: replacements.unwrap_or_default(),
        ..Default::default()
    };
    let inner = cosmolkit_core::smarts_parse::mol_from_smarts(smarts, &params)
        .map_err(smarts_parse_pyerr)?;
    Ok(Molecule { inner })
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn parse_smarts_with_params(smarts: &str, params: &SmartsParserParams) -> PyResult<Molecule> {
    // RDKit✔️✔️: ROMol *MolFromSmartsHelper(python::object ismiles,
    // RDKit✔️✔️:                            const SmartsParserParams &params) {
    // RDKit✔️✔️:   std::string smiles = pyObjectToString(ismiles);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     return SmartsToMol(smiles, params);
    // RDKit✔️✔️:   } catch (...) {
    // RDKit✔️✔️:     return nullptr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Complexity review: both convert one Python string and invoke the same
    // canonical compiler once. The cloned parameter value is bounded by the
    // replacements map and introduces no query graph copy or alternate parser.
    let inner = cosmolkit_core::smarts_parse::mol_from_smarts(smarts, &params.inner)
        .map_err(smarts_parse_pyerr)?;
    Ok(Molecule { inner })
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[doc = r#"
Return whether UFF parameters are available for every atom in a molecule.
"#]
fn uff_has_all_molecule_params(mol: &Molecule) -> PyResult<bool> {
    cosmolkit_core::uff_has_all_molecule_params(&mol.inner).map_err(forcefield_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
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

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
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

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[doc = r#"
Return whether MMFF94 parameters are available for a molecule.
"#]
fn mmff_has_all_molecule_params(mol: &Molecule) -> PyResult<bool> {
    cosmolkit_core::mmff_has_all_molecule_params(&mol.inner).map_err(forcefield_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
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

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
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

fn py_substruct_helper(
    mol: &cosmolkit_core::Molecule,
    query: &cosmolkit_core::Molecule,
    params: &cosmolkit_core::SubstructMatchParams,
) -> Vec<cosmolkit_core::SubstructMatchResult> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: void pySubstructHelper(T1 &mol, T2 &query,
    // RDKit✔️✔️:                        const SubstructMatchParameters &params,
    // RDKit✔️✔️:                        std::vector<MatchVectType> &matches) {
    // RDKit✔️✔️:   // it's safe to release the GIL here since the functors wrapping the python
    // RDKit✔️✔️:   // callbacks will reacquire it as needed
    // RDKit✔️✔️:   NOGIL gil;
    // RDKit✔️✔️:   matches = SubstructMatch(mol, query, params);
    // RDKit✔️✔️: }
    // Complexity review: GIL release and one core matcher call add constant
    // overhead; result allocation and callback dispatch match the source path.
    Python::attach(|py| {
        py.detach(|| cosmolkit_core::get_substruct_matches_with_params(mol, query, params))
    })
}

fn help_has_substruct_match(
    mol: &cosmolkit_core::Molecule,
    query: &cosmolkit_core::Molecule,
    params: &cosmolkit_core::SubstructMatchParams,
) -> bool {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: bool helpHasSubstructMatch(T1 &mol, T2 &query,
    // RDKit✔️✔️:                            const SubstructMatchParameters &params) {
    // RDKit✔️✔️:   SubstructMatchParameters ps = params;
    // RDKit✔️✔️:   ps.maxMatches = 1;
    // RDKit✔️✔️:   std::vector<MatchVectType> matches;
    // RDKit✔️✔️:   pySubstructHelper(mol, query, params, matches);
    // RDKit✔️✔️:   return matches.size() != 0;
    // RDKit✔️✔️: }
    // Complexity review: this intentionally preserves RDKit's use of the
    // original params after preparing `ps`; search cost therefore follows the
    // caller's max_matches instead of being forcibly bounded to one.
    let mut _ps = params.clone();
    _ps.max_matches = 1;
    !py_substruct_helper(mol, query, params).is_empty()
}

fn help_get_substruct_match(
    mol: &cosmolkit_core::Molecule,
    query: &cosmolkit_core::Molecule,
    params: &cosmolkit_core::SubstructMatchParams,
) -> Option<cosmolkit_core::SubstructMatchResult> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: PyObject *helpGetSubstructMatch(T1 &mol, T2 &query,
    // RDKit✔️✔️:                                 const SubstructMatchParameters &params) {
    // RDKit✔️✔️:   SubstructMatchParameters ps = params;
    // RDKit✔️✔️:   ps.maxMatches = 1;
    // RDKit✔️✔️:   std::vector<MatchVectType> matches;
    // RDKit✔️✔️:   pySubstructHelper(mol, query, params, matches);
    // RDKit✔️✔️:   MatchVectType match;
    // RDKit✔️✔️:   if (matches.size()) {
    // RDKit✔️✔️:     match = matches[0];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return convertMatches(match);
    // RDKit✔️✔️: }
    // Complexity review: preserves the source's original-params call and
    // returns the first already-allocated result without an extra scan.
    let mut _ps = params.clone();
    _ps.max_matches = 1;
    py_substruct_helper(mol, query, params).into_iter().next()
}

fn help_get_substruct_matches(
    mol: &cosmolkit_core::Molecule,
    query: &cosmolkit_core::Molecule,
    params: &cosmolkit_core::SubstructMatchParams,
) -> Vec<cosmolkit_core::SubstructMatchResult> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: PyObject *helpGetSubstructMatches(T1 &mol, T2 &query,
    // RDKit✔️✔️:                                   const SubstructMatchParameters &params) {
    // RDKit✔️✔️:   std::vector<MatchVectType> matches;
    // RDKit✔️✔️:   pySubstructHelper(mol, query, params, matches);
    // RDKit✔️✔️:   PyObject *res = PyTuple_New(matches.size());
    // RDKit✔️✔️:   for (size_t idx = 0; idx < matches.size(); idx++) {
    // RDKit✔️✔️:     PyTuple_SetItem(res, idx, convertMatches(matches[idx]));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Complexity review: delegates one search and moves each result once,
    // matching the source's linear result conversion.
    py_substruct_helper(mol, query, params)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (
    mol,
    query,
    recursion_possible=true,
    use_chirality=false,
    use_query_query_matches=false
))]
#[doc = r#"
Return whether a molecule contains a molecule-query substructure.

The query must be the canonical query-bearing ``Molecule`` returned by
``parse_smarts``. The ordinary-molecule SMARTS and substructure boundary is
covered by the pinned RDKit parity corpus; reaction and database/container
SMARTS remain outside this API.
"#]
fn has_substruct_match(
    mol: &Molecule,
    query: &Molecule,
    recursion_possible: bool,
    use_chirality: bool,
    use_query_query_matches: bool,
) -> bool {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: bool HasSubstructMatch(T1 &mol, T2 &query, bool recursionPossible = true,
    // RDKit✔️✔️:                        bool useChirality = false,
    // RDKit✔️✔️:                        bool useQueryQueryMatches = false) {
    // RDKit✔️✔️:   NOGIL gil;
    // RDKit✔️✔️:   MatchVectType res;
    // RDKit✔️✔️:   return SubstructMatch(mol, query, res, recursionPossible, useChirality,
    // RDKit✔️✔️:                         useQueryQueryMatches);
    // RDKit✔️✔️: }
    // Complexity review: both paths run one bounded substructure search and
    // stop after the first match. Parameter construction is constant-space.
    let params = cosmolkit_core::SubstructMatchParams {
        max_matches: 1,
        recursion_possible,
        use_chirality,
        use_query_query_matches,
        ..Default::default()
    };
    help_has_substruct_match(&mol.inner, &query.inner, &params)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (mol, query, use_chirality=false, use_query_query_matches=false))]
#[doc = r#"
Return the first molecule-query substructure match, if present.

The ordinary-molecule SMARTS and substructure boundary is covered by the
pinned RDKit parity corpus. Reaction and database/container SMARTS remain
outside this API.
"#]
fn get_substruct_match(
    mol: &Molecule,
    query: &Molecule,
    use_chirality: bool,
    use_query_query_matches: bool,
) -> Option<SubstructMatchResult> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: PyObject *GetSubstructMatch(T1 &mol, T2 &query, bool useChirality = false,
    // RDKit✔️✔️:                             bool useQueryQueryMatches = false) {
    // RDKit✔️✔️:   MatchVectType matches;
    // RDKit✔️✔️:   {
    // RDKit✔️✔️:     NOGIL gil;
    // RDKit✔️✔️:     SubstructMatch(mol, query, matches, true, useChirality,
    // RDKit✔️✔️:                    useQueryQueryMatches);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return convertMatches(matches);
    // RDKit✔️✔️: }
    // Complexity review: both paths bound the search to one result and then
    // move its dense mapping into the Python result in O(1) auxiliary work.
    let params = cosmolkit_core::SubstructMatchParams {
        max_matches: 1,
        use_chirality,
        use_query_query_matches,
        ..Default::default()
    };
    help_get_substruct_match(&mol.inner, &query.inner, &params).map(Into::into)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (
    mol,
    query,
    uniquify=true,
    use_chirality=false,
    use_query_query_matches=false,
    max_matches=1000
))]
#[doc = r#"
Return molecule-query substructure matches.

The ordinary-molecule SMARTS and substructure boundary is covered by the
pinned RDKit parity corpus. Reaction and database/container SMARTS remain
outside this API.
"#]
fn get_substruct_matches(
    mol: &Molecule,
    query: &Molecule,
    uniquify: bool,
    use_chirality: bool,
    use_query_query_matches: bool,
    max_matches: usize,
) -> Vec<SubstructMatchResult> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: PyObject *GetSubstructMatches(T1 &mol, T2 &query, bool uniquify = true,
    // RDKit✔️✔️:                               bool useChirality = false,
    // RDKit✔️✔️:                               bool useQueryQueryMatches = false,
    // RDKit✔️✔️:                               unsigned int maxMatches = 1000) {
    // RDKit✔️✔️:   std::vector<MatchVectType> matches;
    // RDKit✔️✔️:   int matched;
    // RDKit✔️✔️:   {
    // RDKit✔️✔️:     NOGIL gil;
    // RDKit✔️✔️:     matched = SubstructMatch(mol, query, matches, uniquify, true, useChirality,
    // RDKit✔️✔️:                              useQueryQueryMatches, maxMatches);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   PyObject *res = PyTuple_New(matched);
    // RDKit✔️✔️:   for (int idx = 0; idx < matched; idx++) {
    // RDKit✔️✔️:     PyTuple_SetItem(res, idx, convertMatches(matches[idx]));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Complexity review: the core search has the same match limit and
    // uniquification shape; converting each result is linear in result count.
    let params = cosmolkit_core::SubstructMatchParams {
        max_matches,
        uniquify,
        use_chirality,
        use_query_query_matches,
        ..Default::default()
    };
    help_get_substruct_matches(&mol.inner, &query.inner, &params)
        .into_iter()
        .map(Into::into)
        .collect()
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (
    mol,
    query,
    max_matches=1000,
    uniquify=true,
    final_match=None,
    atom_match=None,
    bond_match=None
))]
#[doc = r#"
Return molecule-query substructure matches with explicit limits.

The ordinary-molecule SMARTS and substructure boundary is covered by the
pinned RDKit parity corpus. Reaction and database/container SMARTS remain
outside this API.
"#]
fn get_substruct_matches_with_params(
    mol: &Molecule,
    query: &Molecule,
    max_matches: usize,
    uniquify: bool,
    final_match: Option<Py<PyAny>>,
    atom_match: Option<Py<PyAny>>,
    bond_match: Option<Py<PyAny>>,
) -> PyResult<Vec<SubstructMatchResult>> {
    let callback_error = Arc::new(Mutex::new(None::<PyErr>));
    let extra_final_check = final_match.map(|callback| {
        // RDKit✔️✔️: class pyFinalMatchFunctor : public pyFunctor {
        // RDKit✔️✔️:  public:
        // RDKit✔️✔️:   pyFinalMatchFunctor(python::object obj) : dp_obj(std::move(obj)) {}
        // RDKit✔️✔️:   ~pyFinalMatchFunctor() = default;
        // RDKit✔️✔️:   bool operator()(const ROMol &m, std::span<const unsigned int> match) {
        // RDKit✔️✔️:     // grab the GIL
        // RDKit✔️✔️:     PyGILStateHolder h;
        // RDKit✔️✔️:     // boost::python doesn't handle std::span, so we need to convert the span to
        // RDKit✔️✔️:     // a vector before calling into python:
        // RDKit✔️✔️:     std::vector<unsigned int> matchVec(match.begin(), match.end());
        // RDKit✔️✔️:     return python::extract<bool>(dp_obj(boost::ref(m), boost::ref(matchVec)));
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:  private:
        // RDKit✔️✔️:   python::object dp_obj;
        // RDKit✔️✔️: };
        // Complexity review: both copy the O(query atoms) match span once and
        // call Python once per completed candidate. Rust clones Molecule's
        // Arc-backed blocks cheaply and preserves callback exceptions.
        let callback_error = Arc::clone(&callback_error);
        Arc::new(
            move |matched: &cosmolkit_core::Molecule, atom_ids: &[usize]| {
                Python::attach(|py| {
                    let result = Py::new(
                        py,
                        Molecule {
                            inner: matched.clone(),
                        },
                    )
                    .and_then(|py_molecule| {
                        callback
                            .call1(py, (py_molecule, atom_ids.to_vec()))
                            .and_then(|value| value.extract::<bool>(py))
                    });
                    match result {
                        Ok(accepted) => accepted,
                        Err(error) => {
                            *callback_error
                                .lock()
                                .expect("callback error mutex poisoned") = Some(error);
                            false
                        }
                    }
                })
            },
        ) as cosmolkit_core::ExtraFinalCheck
    });
    let extra_atom_check = atom_match.map(|callback| {
        // RDKit✔️❌: template <typename T>
        // RDKit✔️❌: class pyMatchFunctor : public pyFunctor {
        // RDKit✔️❌:  public:
        // RDKit✔️❌:   pyMatchFunctor(python::object obj) : dp_obj(std::move(obj)) {}
        // RDKit✔️❌:   ~pyMatchFunctor() = default;
        // RDKit✔️❌:   bool operator()(const T &a1, const T &a2) {
        // RDKit✔️❌:     // grab the GIL
        // RDKit✔️❌:     PyGILStateHolder h;
        // RDKit✔️❌:     return python::extract<bool>(dp_obj(boost::ref(a1), boost::ref(a2)));
        // RDKit✔️❌:   }
        // RDKit✔️❌:
        // RDKit✔️❌:  private:
        // RDKit✔️❌:   python::object dp_obj;
        // RDKit✔️❌: };
        // Complexity review: callback dispatch remains O(1), but the Python
        // binding owns read-only Atom snapshots instead of zero-copy borrowed
        // wrappers. Snapshots are built once per match invocation and cloned
        // once per callback, adding allocation relative to Boost.Python.
        let callback_error = Arc::clone(&callback_error);
        let query_atoms = Arc::new(atom_snapshots(&query.inner));
        let target_atoms = Arc::new(atom_snapshots(&mol.inner));
        Arc::new(
            move |_query_molecule: &cosmolkit_core::Molecule,
                  query_atom: &cosmolkit_core::Atom,
                  _target_molecule: &cosmolkit_core::Molecule,
                  target_atom: &cosmolkit_core::Atom| {
                Python::attach(|py| {
                    let result = Py::new(py, query_atoms[query_atom.id().index()].clone())
                        .and_then(|py_query_atom| {
                            Py::new(py, target_atoms[target_atom.id().index()].clone()).and_then(
                                |py_target_atom| {
                                    callback
                                        .call1(py, (py_query_atom, py_target_atom))
                                        .and_then(|value| value.extract::<bool>(py))
                                },
                            )
                        });
                    match result {
                        Ok(accepted) => accepted,
                        Err(error) => {
                            let mut stored = callback_error
                                .lock()
                                .expect("callback error mutex poisoned");
                            if stored.is_none() {
                                *stored = Some(error);
                            }
                            false
                        }
                    }
                })
            },
        ) as cosmolkit_core::ExtraAtomCheck
    });
    let extra_bond_check = bond_match.map(|callback| {
        // The verbatim pyMatchFunctor source body is anchored in the atom
        // specialization immediately above; this is its Bond specialization.
        // Complexity review: the same O(1) dispatch and one snapshot clone per
        // argument apply here; this remains allocation-heavier than RDKit's
        // borrowed Bond wrappers.
        let callback_error = Arc::clone(&callback_error);
        let query_bonds = Arc::new(bond_snapshots(&query.inner));
        let target_bonds = Arc::new(bond_snapshots(&mol.inner));
        Arc::new(
            move |query_bond: &cosmolkit_core::Bond, target_bond: &cosmolkit_core::Bond| {
                Python::attach(|py| {
                    let result = Py::new(py, query_bonds[query_bond.id().index()].clone())
                        .and_then(|py_query_bond| {
                            Py::new(py, target_bonds[target_bond.id().index()].clone()).and_then(
                                |py_target_bond| {
                                    callback
                                        .call1(py, (py_query_bond, py_target_bond))
                                        .and_then(|value| value.extract::<bool>(py))
                                },
                            )
                        });
                    match result {
                        Ok(accepted) => accepted,
                        Err(error) => {
                            let mut stored = callback_error
                                .lock()
                                .expect("callback error mutex poisoned");
                            if stored.is_none() {
                                *stored = Some(error);
                            }
                            false
                        }
                    }
                })
            },
        ) as cosmolkit_core::ExtraBondCheck
    });
    let params = cosmolkit_core::SubstructMatchParams {
        max_matches,
        uniquify,
        use_chirality: false,
        specified_stereo_query_matches_unspecified: false,
        extra_atom_check,
        extra_bond_check,
        extra_final_check,
        ..Default::default()
    };
    let matches = help_get_substruct_matches(&mol.inner, &query.inner, &params)
        .into_iter()
        .map(Into::into)
        .collect();
    if let Some(error) = callback_error
        .lock()
        .expect("callback error mutex poisoned")
        .take()
    {
        return Err(error);
    }
    Ok(matches)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, only_heavy=false))]
#[doc = r#"
Return the RDKit-aligned average molecular weight.

Set ``only_heavy=True`` to omit hydrogen atoms and implicit hydrogen mass.
The input molecule is not mutated.
"#]
fn calc_mol_wt(molecule: PyRef<'_, Molecule>, only_heavy: bool) -> PyResult<f64> {
    cosmolkit_core::calc_mol_wt(&molecule.inner, only_heavy).map_err(descriptor_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, only_heavy=false))]
#[doc = r#"
Return the RDKit-aligned exact molecular weight.

Set ``only_heavy=True`` to omit hydrogen atoms and implicit hydrogen mass.
The input molecule is not mutated.
"#]
fn calc_exact_mol_wt(molecule: PyRef<'_, Molecule>, only_heavy: bool) -> PyResult<f64> {
    cosmolkit_core::calc_exact_mol_wt(&molecule.inner, only_heavy).map_err(descriptor_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, separate_isotopes=false, abbreviate_h_isotopes=true))]
#[doc = r#"
Return the RDKit-aligned molecular formula.

``separate_isotopes`` emits isotope-specific terms. When enabled,
``abbreviate_h_isotopes`` writes hydrogen-2 and hydrogen-3 as D and T.
The input molecule is not mutated.
"#]
fn calc_mol_formula(
    molecule: PyRef<'_, Molecule>,
    separate_isotopes: bool,
    abbreviate_h_isotopes: bool,
) -> PyResult<String> {
    cosmolkit_core::calc_mol_formula(&molecule.inner, separate_isotopes, abbreviate_h_isotopes)
        .map_err(descriptor_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[doc = "Return the RDKit-aligned hydrogen-bond donor count without mutating the molecule."]
fn calc_num_hbd(molecule: PyRef<'_, Molecule>) -> PyResult<u32> {
    cosmolkit_core::calc_num_hbd(&molecule.inner).map_err(descriptor_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[doc = "Return the RDKit-aligned hydrogen-bond acceptor count without mutating the molecule."]
fn calc_num_hba(molecule: PyRef<'_, Molecule>) -> PyResult<u32> {
    cosmolkit_core::calc_num_hba(&molecule.inner).map_err(descriptor_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[doc = "Return the RDKit-aligned fraction of carbon atoms that are sp3 without mutating the molecule."]
fn calc_fraction_csp3(molecule: PyRef<'_, Molecule>) -> PyResult<f64> {
    cosmolkit_core::calc_fraction_csp3(&molecule.inner).map_err(descriptor_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pyfunction]
#[gen_stub(override_return_type(type_repr = "tuple[builtins.float, builtins.float]", imports = ("builtins")))]
#[pyo3(signature = (molecule, include_hs=true, force=false))]
#[doc = r#"
Return ``(logp, molar_refractivity)`` using the RDKit Crippen descriptor path.

The chemical graph and user-visible molecule properties are not mutated. The
source-compatible computed Crippen descriptor cache may be populated on the
input molecule, matching RDKit's property-cache behavior.
"#]
fn calc_crippen_descriptors(
    molecule: PyRef<'_, Molecule>,
    include_hs: bool,
    force: bool,
) -> PyResult<(f64, f64)> {
    cosmolkit_core::calc_crippen_descriptors(&molecule.inner, include_hs, force)
        .map(|values| (values.logp, values.molar_refractivity))
        .map_err(descriptor_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, force=false, include_sandp=false))]
#[doc = r#"
Return the RDKit-aligned topological polar surface area.

Set ``include_sandp=True`` to include sulfur and phosphorus contributions.
The input molecule is not mutated.
"#]
fn calc_tpsa(molecule: PyRef<'_, Molecule>, force: bool, include_sandp: bool) -> PyResult<f64> {
    cosmolkit_core::calc_tpsa(&molecule.inner, force, include_sandp).map_err(descriptor_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[doc = "Return the RDKit-aligned aromatic-ring count without mutating the molecule."]
fn calc_num_aromatic_rings(molecule: PyRef<'_, Molecule>) -> PyResult<u32> {
    cosmolkit_core::calc_num_aromatic_rings(&molecule.inner).map_err(descriptor_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pyfunction]
#[pyo3(signature = (molecule, mode="default"))]
#[doc = r#"
Return the RDKit-aligned rotatable-bond count.

``mode`` must be ``"default"``, ``"non_strict"``, ``"strict"``, or
``"strict_linkages"``. The input molecule is not mutated.
"#]
fn calc_num_rotatable_bonds(
    molecule: PyRef<'_, Molecule>,
    #[gen_stub(override_type(type_repr = "typing.Literal['default', 'non_strict', 'strict', 'strict_linkages']", imports = ("typing")))]
    mode: &str,
) -> PyResult<u32> {
    let options = parse_rotatable_bonds_mode(mode)?;
    cosmolkit_core::calc_num_rotatable_bonds(&molecule.inner, options).map_err(descriptor_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[doc = "Return the RDKit-aligned quantitative estimate of drug-likeness without mutating the molecule."]
fn calc_qed(molecule: PyRef<'_, Molecule>) -> PyResult<f64> {
    cosmolkit_core::calc_qed(&molecule.inner).map_err(descriptor_pyerr)
}

macro_rules! python_infallible_float_descriptor {
    ($name:ident, $doc:literal) => {
        #[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
        #[pyfunction]
        #[doc = $doc]
        fn $name(molecule: PyRef<'_, Molecule>) -> f64 {
            cosmolkit_core::$name(&molecule.inner)
        }
    };
}

macro_rules! python_float_descriptor {
    ($name:ident, $doc:literal) => {
        #[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
        #[pyfunction]
        #[doc = $doc]
        fn $name(molecule: PyRef<'_, Molecule>) -> PyResult<f64> {
            cosmolkit_core::$name(&molecule.inner).map_err(descriptor_pyerr)
        }
    };
}

macro_rules! python_count_descriptor {
    ($name:ident, $doc:literal) => {
        #[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
        #[pyfunction]
        #[doc = $doc]
        fn $name(molecule: PyRef<'_, Molecule>) -> PyResult<u32> {
            cosmolkit_core::$name(&molecule.inner).map_err(descriptor_pyerr)
        }
    };
}

python_infallible_float_descriptor!(calc_chi_0, "Return the graph-degree Chi0 descriptor.");
python_infallible_float_descriptor!(calc_chi_1, "Return the graph-degree Chi1 descriptor.");

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, order, force=false))]
#[doc = "Return the order-N valence connectivity descriptor."]
fn calc_chi_nv(molecule: PyRef<'_, Molecule>, order: usize, force: bool) -> PyResult<f64> {
    cosmolkit_core::calc_chi_nv(&molecule.inner, order, force).map_err(descriptor_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, order, force=false))]
#[doc = "Return the order-N principal-quantum connectivity descriptor."]
fn calc_chi_nn(molecule: PyRef<'_, Molecule>, order: usize, force: bool) -> PyResult<f64> {
    cosmolkit_core::calc_chi_nn(&molecule.inner, order, force).map_err(descriptor_pyerr)
}

macro_rules! python_fixed_chi_descriptor {
    ($name:ident) => {
        #[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
        #[pyfunction]
        #[pyo3(signature = (molecule, force=false))]
        fn $name(molecule: PyRef<'_, Molecule>, force: bool) -> PyResult<f64> {
            cosmolkit_core::$name(&molecule.inner, force).map_err(descriptor_pyerr)
        }
    };
}

python_fixed_chi_descriptor!(calc_chi_0v);
python_fixed_chi_descriptor!(calc_chi_1v);
python_fixed_chi_descriptor!(calc_chi_2v);
python_fixed_chi_descriptor!(calc_chi_3v);
python_fixed_chi_descriptor!(calc_chi_4v);
python_fixed_chi_descriptor!(calc_chi_0n);
python_fixed_chi_descriptor!(calc_chi_1n);
python_fixed_chi_descriptor!(calc_chi_2n);
python_fixed_chi_descriptor!(calc_chi_3n);
python_fixed_chi_descriptor!(calc_chi_4n);

python_infallible_float_descriptor!(
    calc_hall_kier_alpha,
    "Return the Hall-Kier alpha descriptor."
);

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[doc = "Return Hall-Kier alpha and atom-index-aligned contributions."]
fn calc_hall_kier_alpha_with_contributions(molecule: PyRef<'_, Molecule>) -> (f64, Vec<f64>) {
    let values = cosmolkit_core::calc_hall_kier_alpha_with_contributions(&molecule.inner);
    (values.alpha, values.atom_contributions)
}

python_infallible_float_descriptor!(calc_kappa_1, "Return the first Kappa shape index.");
python_float_descriptor!(calc_kappa_2, "Return the second Kappa shape index.");
python_float_descriptor!(calc_kappa_3, "Return the third Kappa shape index.");
python_float_descriptor!(calc_phi, "Return the molecular flexibility Phi descriptor.");

python_count_descriptor!(
    calc_lipinski_hba,
    "Return the direct Lipinski nitrogen/oxygen acceptor count."
);
python_count_descriptor!(
    calc_lipinski_hbd,
    "Return the direct Lipinski nitrogen/oxygen donor-site count."
);
python_count_descriptor!(calc_num_heteroatoms, "Return the heteroatom count.");
python_count_descriptor!(calc_num_amide_bonds, "Return the amide-bond count.");
python_count_descriptor!(calc_num_heavy_atoms, "Return the heavy-atom count.");
python_count_descriptor!(
    calc_num_atoms,
    "Return the total atom count including implicit hydrogens."
);
python_count_descriptor!(calc_num_rings, "Return the SSSR ring count.");
python_count_descriptor!(calc_num_heterocycles, "Return the heterocycle count.");
python_count_descriptor!(calc_num_saturated_rings, "Return the saturated-ring count.");
python_count_descriptor!(calc_num_aliphatic_rings, "Return the aliphatic-ring count.");
python_count_descriptor!(
    calc_num_aromatic_heterocycles,
    "Return the aromatic heterocycle count."
);
python_count_descriptor!(
    calc_num_aromatic_carbocycles,
    "Return the aromatic carbocycle count."
);
python_count_descriptor!(
    calc_num_aliphatic_heterocycles,
    "Return the aliphatic heterocycle count."
);
python_count_descriptor!(
    calc_num_aliphatic_carbocycles,
    "Return the aliphatic carbocycle count."
);
python_count_descriptor!(
    calc_num_saturated_heterocycles,
    "Return the saturated heterocycle count."
);
python_count_descriptor!(
    calc_num_saturated_carbocycles,
    "Return the saturated carbocycle count."
);
python_count_descriptor!(calc_num_spiro_atoms, "Return the spiro-atom count.");
python_count_descriptor!(
    calc_num_bridgehead_atoms,
    "Return the bridgehead-atom count."
);
python_count_descriptor!(
    calc_num_atom_stereo_centers,
    "Return the number of possible atom stereocenters."
);
python_count_descriptor!(
    calc_num_unspecified_atom_stereo_centers,
    "Return the number of possible atom stereocenters without a specified chiral tag."
);

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[doc = "Return the fixed-order 42-component molecular quantum number vector."]
fn calc_mqns(molecule: PyRef<'_, Molecule>) -> PyResult<Vec<u32>> {
    cosmolkit_core::calc_mqns(&molecule.inner)
        .map(|values| values.to_vec())
        .map_err(descriptor_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, include_hydrogens=true, force=false))]
#[doc = "Return the Labute approximate surface area."]
fn calc_labute_asa(molecule: PyRef<'_, Molecule>, include_hydrogens: bool, force: bool) -> f64 {
    cosmolkit_core::calc_labute_asa(&molecule.inner, include_hydrogens, force)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, include_hydrogens=true, force=false))]
#[doc = "Return Labute ASA, atom-index-aligned contributions, and the hydrogen contribution."]
fn calc_labute_asa_contributions(
    molecule: PyRef<'_, Molecule>,
    include_hydrogens: bool,
    force: bool,
) -> (f64, Vec<f64>, f64) {
    let values =
        cosmolkit_core::calc_labute_asa_contributions(&molecule.inner, include_hydrogens, force);
    (
        values.asa,
        values.atom_contributions,
        values.hydrogen_contribution,
    )
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, bins=None, force=false))]
#[doc = "Return the SlogP-VSA vector using default or caller-provided bin boundaries."]
fn calc_slogp_vsa(
    molecule: PyRef<'_, Molecule>,
    bins: Option<Vec<f64>>,
    force: bool,
) -> PyResult<Vec<f64>> {
    match bins {
        Some(bins) => cosmolkit_core::calc_slogp_vsa_with_bins(&molecule.inner, &bins, force),
        None => {
            cosmolkit_core::calc_slogp_vsa(&molecule.inner, force).map(|values| values.to_vec())
        }
    }
    .map_err(descriptor_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, bins=None, force=false))]
#[doc = "Return the SMR-VSA vector using default or caller-provided bin boundaries."]
fn calc_smr_vsa(
    molecule: PyRef<'_, Molecule>,
    bins: Option<Vec<f64>>,
    force: bool,
) -> PyResult<Vec<f64>> {
    match bins {
        Some(bins) => cosmolkit_core::calc_smr_vsa_with_bins(&molecule.inner, &bins, force),
        None => cosmolkit_core::calc_smr_vsa(&molecule.inner, force).map(|values| values.to_vec()),
    }
    .map_err(descriptor_pyerr)
}

python_float_descriptor!(calc_slogp_vsa_1, "Return SlogP-VSA bin 1.");
python_float_descriptor!(calc_slogp_vsa_2, "Return SlogP-VSA bin 2.");
python_float_descriptor!(calc_slogp_vsa_3, "Return SlogP-VSA bin 3.");
python_float_descriptor!(calc_slogp_vsa_4, "Return SlogP-VSA bin 4.");
python_float_descriptor!(calc_slogp_vsa_5, "Return SlogP-VSA bin 5.");
python_float_descriptor!(calc_slogp_vsa_6, "Return SlogP-VSA bin 6.");
python_float_descriptor!(calc_slogp_vsa_7, "Return SlogP-VSA bin 7.");
python_float_descriptor!(calc_slogp_vsa_8, "Return SlogP-VSA bin 8.");
python_float_descriptor!(calc_slogp_vsa_9, "Return SlogP-VSA bin 9.");
python_float_descriptor!(calc_slogp_vsa_10, "Return SlogP-VSA bin 10.");
python_float_descriptor!(calc_slogp_vsa_11, "Return SlogP-VSA bin 11.");
python_float_descriptor!(calc_slogp_vsa_12, "Return SlogP-VSA bin 12.");
python_float_descriptor!(calc_smr_vsa_1, "Return SMR-VSA bin 1.");
python_float_descriptor!(calc_smr_vsa_2, "Return SMR-VSA bin 2.");
python_float_descriptor!(calc_smr_vsa_3, "Return SMR-VSA bin 3.");
python_float_descriptor!(calc_smr_vsa_4, "Return SMR-VSA bin 4.");
python_float_descriptor!(calc_smr_vsa_5, "Return SMR-VSA bin 5.");
python_float_descriptor!(calc_smr_vsa_6, "Return SMR-VSA bin 6.");
python_float_descriptor!(calc_smr_vsa_7, "Return SMR-VSA bin 7.");
python_float_descriptor!(calc_smr_vsa_8, "Return SMR-VSA bin 8.");
python_float_descriptor!(calc_smr_vsa_9, "Return SMR-VSA bin 9.");
python_float_descriptor!(calc_smr_vsa_10, "Return SMR-VSA bin 10.");

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pyfunction]
#[gen_stub(override_return_type(type_repr = "Element"))]
#[doc = "Return the element represented by a canonical or source-recognized symbol."]
fn element_from_symbol<'py>(py: Python<'py>, symbol: &str) -> PyResult<Bound<'py, PyAny>> {
    let element = cosmolkit_core::Element::from_symbol(symbol)
        .ok_or_else(|| PyValueError::new_err(format!("unknown element symbol '{symbol}'")))?;
    element_enum_member(py, element)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pyfunction]
#[doc = "Return source-aligned periodic-table metadata for an atomic number."]
fn get_element_info(atomic_number: u8) -> PyResult<PyElementInfo> {
    let element = cosmolkit_core::Element::from_atomic_number(atomic_number).ok_or_else(|| {
        PyValueError::new_err(format!(
            "atomic number {atomic_number} is outside the Element domain 0..=118"
        ))
    })?;
    Ok(PyElementInfo {
        inner: element.info(),
    })
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[doc = r#"
Return Gemmi-derived tabulated residue information for a residue name.
"#]
fn find_tabulated_residue(name: &str) -> PyResidueInfo {
    PyResidueInfo {
        inner: cosmolkit_core::find_tabulated_residue(name),
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[doc = r#"
Return the Gemmi tabulated residue index for a residue name.
"#]
fn find_tabulated_residue_idx(name: &str) -> usize {
    cosmolkit_core::find_tabulated_residue_idx(name)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
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

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pyfunction]
#[gen_stub(override_return_type(type_repr = "ResidueCode"))]
#[doc = r#"
Return the Gemmi tabulated residue code for a residue name.
"#]
fn residue_code_from_name<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
    residue_code_enum_member(py, cosmolkit_core::residue_code_from_name(name))
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
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

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pyfunction]
#[doc = r#"
Expand a one-letter amino-acid, RNA, or DNA residue sequence using Gemmi's table.
"#]
fn expand_one_letter_sequence(
    seq: &str,
    #[gen_stub(override_type(type_repr = "ResidueInfoKind"))] kind: i64,
) -> PyResult<Vec<String>> {
    let kind = residue_info_kind_from_code(kind)?;
    cosmolkit_core::expand_one_letter_sequence(seq, kind).map_err(residue_info_sequence_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[doc = r#"
Expand a deprecated Gemmi protein one-letter residue sequence alias.
"#]
fn expand_protein_one_letter_string(seq: &str) -> PyResult<Vec<String>> {
    cosmolkit_core::expand_protein_one_letter_string(seq).map_err(residue_info_sequence_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pyfunction]
#[doc = r#"
Expand a one-letter amino-acid, RNA, or DNA residue code using Gemmi's table.
"#]
fn expand_one_letter(
    code: &str,
    #[gen_stub(override_type(type_repr = "ResidueInfoKind"))] kind: i64,
) -> PyResult<Option<String>> {
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

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pyfunction]
#[doc = r#"
Generate an InChIKey directly from an InChI string.
"#]
fn inchi_to_key(inchi: &str) -> PyResult<Option<String>> {
    let output = cosmolkit_core::inchi_to_inchi_key(inchi.as_bytes()).map_err(inchi_pyerr)?;
    let failed = output.key.is_empty()
        && output
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.level == cosmolkit_core::InchiDiagnosticLevel::Error);
    emit_inchi_diagnostics(&output.diagnostics)?;
    if failed {
        return Ok(None);
    }
    inchi_output_string("inchi_to_inchi_key", "InChIKey", output.key).map(Some)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (include_chirality=false, torsion_atom_count=4, count_simulation=true, count_bounds=None, fp_size=2048, atom_invariants_generator=None))]
fn get_topological_torsion_generator(
    include_chirality: bool,
    torsion_atom_count: u32,
    count_simulation: bool,
    count_bounds: Option<Vec<u32>>,
    fp_size: u32,
    atom_invariants_generator: Option<&PyAtomPairAtomInvariantsGenerator>,
) -> PyResult<PyTopologicalTorsionFingerprintGenerator> {
    // RDKit source: TopologicalTorsionWrapper.cpp lines 21-45.
    // RDKit✔️✔️: FingerprintGenerator<OutputType> *getTopologicalTorsionFPGenerator(
    // RDKit✔️✔️:     const bool includeChirality, const uint32_t torsionAtomCount,
    // RDKit✔️✔️:     const bool countSimulation, python::object &py_countBounds,
    // RDKit✔️✔️:     const std::uint32_t fpSize, python::object &py_atomInvGen) {
    // RDKit✔️✔️:   AtomInvariantsGenerator *atomInvariantsGenerator = nullptr;
    // RDKit✔️✔️:   ... atomInvariantsGenerator = atomInvariantsGenerator->clone();
    // RDKit✔️✔️:   std::vector<std::uint32_t> countBounds = {1, 2, 4, 8};
    // RDKit✔️✔️:   if (py_countBounds) { ... countBounds = *tmp; }
    // RDKit✔️✔️:   return TopologicalTorsion::getTopologicalTorsionGenerator<OutputType>(
    // RDKit✔️✔️:       includeChirality, torsionAtomCount, atomInvariantsGenerator,
    // RDKit✔️✔️:       countSimulation, fpSize, countBounds, false);
    // RDKit✔️✔️: }
    let arguments = cosmolkit_core::fingerprint::TopologicalTorsionArguments::new(
        include_chirality,
        torsion_atom_count,
        count_simulation,
        count_bounds.unwrap_or_else(|| vec![1, 2, 4, 8]),
        fp_size,
    )
    .map_err(fingerprint_pyerr)?;
    let atom_invariants_generator = atom_invariants_generator.map_or_else(
        || cosmolkit_core::fingerprint::AtomPairAtomInvGenerator::new(include_chirality, true),
        |generator| generator.inner.clone(),
    );
    let state = PyTopologicalTorsionGeneratorState {
        arguments,
        atom_invariants_generator,
    };
    state.build().map_err(fingerprint_pyerr)?;
    Ok(PyTopologicalTorsionFingerprintGenerator {
        state: Arc::new(Mutex::new(state)),
    })
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn topological_torsion_generator_from_json(
    json: &str,
) -> PyResult<PyTopologicalTorsionFingerprintGenerator> {
    let generator =
        cosmolkit_core::fingerprint::generatorFromJSON(json).map_err(fingerprint_pyerr)?;
    let cosmolkit_core::fingerprint::TypedFingerprintGenerator::TopologicalTorsion(generator) =
        generator
    else {
        return Err(PyValueError::new_err(
            "JSON does not describe a Topological Torsion fingerprint generator",
        ));
    };
    Ok(PyTopologicalTorsionFingerprintGenerator {
        state: Arc::new(Mutex::new(PyTopologicalTorsionGeneratorState {
            arguments: generator.fingerprint_arguments,
            atom_invariants_generator: generator.atom_invariants_generator,
        })),
    })
}

fn topological_torsion_legacy_params(
    kind: cosmolkit_core::TopologicalTorsionLegacyKind,
    n_bits: u32,
    torsion_atom_count: u32,
    from_atoms: Option<Vec<usize>>,
    ignore_atoms: Option<Vec<usize>>,
    atom_invariants: Option<Vec<u32>>,
    n_bits_per_entry: u32,
    include_chirality: bool,
) -> cosmolkit_core::TopologicalTorsionLegacyParams {
    cosmolkit_core::TopologicalTorsionLegacyParams {
        kind,
        n_bits,
        torsion_atom_count,
        from_atoms: normalize_fingerprint_indices(from_atoms),
        ignore_atoms: normalize_fingerprint_indices(ignore_atoms),
        atom_invariants: normalize_fingerprint_indices(atom_invariants),
        n_bits_per_entry,
        include_chirality,
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, torsion_atom_count=4, from_atoms=None, ignore_atoms=None, atom_invariants=None, include_chirality=false))]
fn get_topological_torsion_fingerprint(
    molecule: &Molecule,
    torsion_atom_count: u32,
    from_atoms: Option<Vec<usize>>,
    ignore_atoms: Option<Vec<usize>>,
    atom_invariants: Option<Vec<u32>>,
    include_chirality: bool,
) -> PyResult<PySparseCountFingerprint> {
    // RDKit source: rdMolDescriptors.cpp lines 250-274.
    // RDKit✔️✔️: RDKit::SparseIntVect<boost::int64_t> *res;
    // RDKit✔️✔️: res = RDKit::AtomPairs::getTopologicalTorsionFingerprint(
    // RDKit✔️✔️:     mol, targetSize, fvect.get(), ivect.get(), invvect.get(),
    // RDKit✔️✔️:     includeChirality);
    let params = topological_torsion_legacy_params(
        cosmolkit_core::TopologicalTorsionLegacyKind::UnfoldedCount,
        2048,
        torsion_atom_count,
        from_atoms,
        ignore_atoms,
        atom_invariants,
        4,
        include_chirality,
    );
    match cosmolkit_core::topological_torsion_legacy_fingerprint(&molecule.inner, &params)
        .map_err(fingerprint_pyerr)?
    {
        cosmolkit_core::TopologicalTorsionLegacyResult::SparseCount(inner) => {
            Ok(PySparseCountFingerprint { inner })
        }
        cosmolkit_core::TopologicalTorsionLegacyResult::Bit(_) => Err(PyRuntimeError::new_err(
            "legacy unfolded Topological Torsion returned an unexpected bit vector",
        )),
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, n_bits=2048, torsion_atom_count=4, from_atoms=None, ignore_atoms=None, atom_invariants=None, include_chirality=false))]
fn get_hashed_topological_torsion_fingerprint(
    molecule: &Molecule,
    n_bits: u32,
    torsion_atom_count: u32,
    from_atoms: Option<Vec<usize>>,
    ignore_atoms: Option<Vec<usize>>,
    atom_invariants: Option<Vec<u32>>,
    include_chirality: bool,
) -> PyResult<PySparseCountFingerprint> {
    // RDKit source: rdMolDescriptors.cpp lines 276-292.
    // RDKit✔️✔️: res = RDKit::AtomPairs::getHashedTopologicalTorsionFingerprint(
    // RDKit✔️✔️:     mol, nBits, targetSize, fvect.get(), ivect.get(), invvect.get(),
    // RDKit✔️✔️:     includeChirality);
    let params = topological_torsion_legacy_params(
        cosmolkit_core::TopologicalTorsionLegacyKind::HashedCount,
        n_bits,
        torsion_atom_count,
        from_atoms,
        ignore_atoms,
        atom_invariants,
        4,
        include_chirality,
    );
    match cosmolkit_core::topological_torsion_legacy_fingerprint(&molecule.inner, &params)
        .map_err(fingerprint_pyerr)?
    {
        cosmolkit_core::TopologicalTorsionLegacyResult::SparseCount(inner) => {
            Ok(PySparseCountFingerprint { inner })
        }
        cosmolkit_core::TopologicalTorsionLegacyResult::Bit(_) => Err(PyRuntimeError::new_err(
            "legacy hashed-count Topological Torsion returned an unexpected bit vector",
        )),
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, n_bits=2048, torsion_atom_count=4, from_atoms=None, ignore_atoms=None, atom_invariants=None, n_bits_per_entry=4, include_chirality=false))]
fn get_hashed_topological_torsion_fingerprint_as_bit_vect(
    molecule: &Molecule,
    n_bits: u32,
    torsion_atom_count: u32,
    from_atoms: Option<Vec<usize>>,
    ignore_atoms: Option<Vec<usize>>,
    atom_invariants: Option<Vec<u32>>,
    n_bits_per_entry: u32,
    include_chirality: bool,
) -> PyResult<Fingerprint> {
    // RDKit source: rdMolDescriptors.cpp lines 294-309.
    // RDKit✔️✔️: res = RDKit::AtomPairs::getHashedTopologicalTorsionFingerprintAsBitVect(
    // RDKit✔️✔️:     mol, nBits, targetSize, fvect.get(), ivect.get(), invvect.get(),
    // RDKit✔️✔️:     nBitsPerEntry, includeChirality);
    let params = topological_torsion_legacy_params(
        cosmolkit_core::TopologicalTorsionLegacyKind::HashedBit,
        n_bits,
        torsion_atom_count,
        from_atoms,
        ignore_atoms,
        atom_invariants,
        n_bits_per_entry,
        include_chirality,
    );
    match cosmolkit_core::topological_torsion_legacy_fingerprint(&molecule.inner, &params)
        .map_err(fingerprint_pyerr)?
    {
        cosmolkit_core::TopologicalTorsionLegacyResult::Bit(inner) => Ok(Fingerprint { inner }),
        cosmolkit_core::TopologicalTorsionLegacyResult::SparseCount(_) => {
            Err(PyRuntimeError::new_err(
                "legacy hashed-bit Topological Torsion returned an unexpected count vector",
            ))
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, atom_index, branch_subtract=0, include_chirality=false))]
fn get_atom_code(
    molecule: &Molecule,
    atom_index: usize,
    branch_subtract: u32,
    include_chirality: bool,
) -> PyResult<u32> {
    // RDKit source: rdMolDescriptors.cpp lines 978-982.
    // RDKit✔️✔️: python::def("GetAtomPairAtomCode", RDKit::AtomPairs::getAtomCode,
    // RDKit✔️✔️:             (python::arg("atom"), python::arg("branchSubtract") = 0,
    // RDKit✔️✔️:              python::arg("includeChirality") = false),
    // RDKit✔️✔️:             docString.c_str());
    cosmolkit_core::get_atom_code(
        &molecule.inner,
        cosmolkit_core::AtomId::new(atom_index),
        branch_subtract,
        include_chirality,
    )
    .map_err(fingerprint_pyerr)
}

fn atom_pair_symbol(type_index: usize) -> String {
    let atomic_number = cosmolkit_core::ATOM_PAIR_ATOM_NUMBER_TYPES
        .get(type_index)
        .copied();
    atomic_number
        .and_then(|value| u8::try_from(value).ok())
        .and_then(cosmolkit_core::Element::from_atomic_number)
        .map_or_else(|| "X".to_string(), |element| element.symbol().to_string())
}

fn explain_atom_code_fields(
    mut code: u32,
    branch_subtract: u32,
    include_chirality: bool,
) -> (String, u32, u32, Option<&'static str>) {
    // RDKit source: rdkit/Chem/AtomPairs/Utils.py lines 75-99.
    // RDKit✔️✔️: typeMask = (1 << rdMolDescriptors.AtomPairsParameters.numTypeBits) - 1
    // RDKit✔️✔️: branchMask = (1 << rdMolDescriptors.AtomPairsParameters.numBranchBits) - 1
    // RDKit✔️✔️: piMask = (1 << rdMolDescriptors.AtomPairsParameters.numPiBits) - 1
    // RDKit✔️✔️: chiMask = (1 << rdMolDescriptors.AtomPairsParameters.numChiralBits) - 1
    // RDKit✔️✔️: nBranch = int(code & branchMask)
    // RDKit✔️✔️: code = code >> rdMolDescriptors.AtomPairsParameters.numBranchBits
    // RDKit✔️✔️: nPi = int(code & piMask)
    // RDKit✔️✔️: code = code >> rdMolDescriptors.AtomPairsParameters.numPiBits
    // RDKit✔️✔️: typeIdx = int(code & typeMask)
    // RDKit✔️✔️: if not includeChirality: return (atomSymbol, nBranch, nPi)
    // RDKit✔️✔️: code = code >> rdMolDescriptors.AtomPairsParameters.numTypeBits
    // RDKit✔️✔️: chiDict = {0: '', 1: 'R', 2: 'S'}
    let branch_mask = (1 << cosmolkit_core::ATOM_PAIR_NUM_BRANCH_BITS) - 1;
    let pi_mask = (1 << cosmolkit_core::ATOM_PAIR_NUM_PI_BITS) - 1;
    let type_mask = (1 << cosmolkit_core::ATOM_PAIR_NUM_TYPE_BITS) - 1;
    let n_branch = (code & branch_mask).wrapping_add(branch_subtract);
    code >>= cosmolkit_core::ATOM_PAIR_NUM_BRANCH_BITS;
    let n_pi = code & pi_mask;
    code >>= cosmolkit_core::ATOM_PAIR_NUM_PI_BITS;
    let type_index = (code & type_mask) as usize;
    let symbol = atom_pair_symbol(type_index);
    if !include_chirality {
        return (symbol, n_branch, n_pi, None);
    }
    code >>= cosmolkit_core::ATOM_PAIR_NUM_TYPE_BITS;
    let chi_mask = (1 << cosmolkit_core::ATOM_PAIR_NUM_CHIRAL_BITS) - 1;
    let chirality = match code & chi_mask {
        0 => "",
        1 => "R",
        2 => "S",
        _ => "",
    };
    (symbol, n_branch, n_pi, Some(chirality))
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pyfunction]
#[gen_stub(override_return_type(
    type_repr = "typing.Union[tuple[builtins.str, builtins.int, builtins.int], tuple[builtins.str, builtins.int, builtins.int, builtins.str]]",
    imports = ("builtins", "typing")
))]
#[pyo3(signature = (code, branch_subtract=0, include_chirality=false))]
fn explain_atom_code<'py>(
    py: Python<'py>,
    code: u32,
    branch_subtract: u32,
    include_chirality: bool,
) -> PyResult<Bound<'py, PyTuple>> {
    let (symbol, branches, pi_electrons, chirality) =
        explain_atom_code_fields(code, branch_subtract, include_chirality);
    if let Some(chirality) = chirality {
        PyTuple::new(
            py,
            [
                symbol.into_pyobject(py)?.into_any(),
                branches.into_pyobject(py)?.into_any(),
                pi_electrons.into_pyobject(py)?.into_any(),
                chirality.into_pyobject(py)?.into_any(),
            ],
        )
    } else {
        PyTuple::new(
            py,
            [
                symbol.into_pyobject(py)?.into_any(),
                branches.into_pyobject(py)?.into_any(),
                pi_electrons.into_pyobject(py)?.into_any(),
            ],
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, path, size, atom_codes=None))]
fn py_score_path(
    molecule: &Molecule,
    path: Vec<usize>,
    size: usize,
    atom_codes: Option<Vec<u32>>,
) -> PyResult<u64> {
    // RDKit source: rdkit/Chem/AtomPairs/Torsions.py lines 31-75.
    // RDKit✔️✔️: codes = [None] * size
    // RDKit✔️✔️: for i in range(size):
    // RDKit✔️✔️:   if i == 0 or i == size - 1: sub = 1
    // RDKit✔️✔️:   else: sub = 2
    // RDKit✔️✔️:   if not atomCodes:
    // RDKit✔️✔️:     codes[i] = Utils.GetAtomCode(mol.GetAtomWithIdx(path[i]), sub)
    // RDKit✔️✔️:   else: codes[i] = atomCodes[path[i]] - sub
    // RDKit✔️✔️: ... canonicalize ... accum |= code << (codeSize * i)
    if size == 0 {
        return Err(PyValueError::new_err("size must be greater than zero"));
    }
    if path.len() < size {
        return Err(PyValueError::new_err(
            "path must contain at least size atom indices",
        ));
    }
    let atom_codes = normalize_fingerprint_indices(atom_codes);
    if atom_codes
        .as_ref()
        .is_some_and(|codes| codes.len() < molecule.inner.num_atoms())
    {
        return Err(PyValueError::new_err(
            "atom_codes must contain at least one entry per molecule atom",
        ));
    }
    let mut codes = Vec::with_capacity(size);
    for (position, &atom_index) in path.iter().take(size).enumerate() {
        if atom_index >= molecule.inner.num_atoms() {
            return Err(PyIndexError::new_err(format!(
                "atom index {atom_index} is out of range"
            )));
        }
        let subtract = if position == 0 || position + 1 == size {
            1
        } else {
            2
        };
        let code = if let Some(atom_codes) = atom_codes.as_ref() {
            atom_codes[atom_index]
                .checked_sub(subtract)
                .ok_or_else(|| {
                    PyValueError::new_err("atom code is smaller than the path branch subtraction")
                })?
        } else {
            cosmolkit_core::get_atom_code(
                &molecule.inner,
                cosmolkit_core::AtomId::new(atom_index),
                subtract,
                false,
            )
            .map_err(fingerprint_pyerr)?
        };
        codes.push(code);
    }
    cosmolkit_core::get_topological_torsion_code(&codes, false).map_err(fingerprint_pyerr)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pyfunction]
#[gen_stub(override_return_type(
    type_repr = "tuple[tuple[builtins.str, builtins.int, builtins.int], ...]",
    imports = ("builtins")
))]
#[pyo3(signature = (score, size=4))]
fn explain_path_score<'py>(
    py: Python<'py>,
    mut score: u64,
    size: usize,
) -> PyResult<Bound<'py, PyTuple>> {
    // RDKit source: rdkit/Chem/AtomPairs/Torsions.py lines 78-150.
    // RDKit✔️✔️: codeSize = rdMolDescriptors.AtomPairsParameters.codeSize
    // RDKit✔️✔️: codeMask = (1 << codeSize) - 1
    // RDKit✔️✔️: for i in range(size):
    // RDKit✔️✔️:   if i == 0 or i == size - 1: sub = 1
    // RDKit✔️✔️:   else: sub = 2
    // RDKit✔️✔️:   code = score & codeMask; score = score >> codeSize
    // RDKit✔️✔️:   symb, nBranch, nPi = Utils.ExplainAtomCode(code)
    // RDKit✔️✔️:   res[i] = (symb, nBranch + sub, nPi)
    let code_mask = (1_u64 << cosmolkit_core::ATOM_PAIR_CODE_SIZE) - 1;
    let mut entries = Vec::with_capacity(size);
    for position in 0..size {
        let subtract = if position == 0 || position + 1 == size {
            1
        } else {
            2
        };
        let code = (score & code_mask) as u32;
        score >>= cosmolkit_core::ATOM_PAIR_CODE_SIZE;
        let (symbol, branches, pi_electrons, _) = explain_atom_code_fields(code, subtract, false);
        entries.push((symbol, branches, pi_electrons));
    }
    PyTuple::new(py, entries)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (molecule, torsion_atom_count=4))]
fn get_topological_torsion_fingerprint_as_ids(
    molecule: &Molecule,
    torsion_atom_count: u32,
) -> PyResult<Vec<u64>> {
    // RDKit source: rdkit/Chem/AtomPairs/Torsions.py lines 153-159.
    // RDKit✔️✔️: nonZeroElements = GetTopologicalTorsionFingerprint(mol, targetSize).GetNonzeroElements()
    // RDKit✔️✔️: frequencies = sorted(nonZeroElements.items())
    // RDKit✔️✔️: res = []
    // RDKit✔️✔️: for k, v in frequencies: res.extend([k] * v)
    let fingerprint =
        get_topological_torsion_fingerprint(molecule, torsion_atom_count, None, None, None, false)?;
    let mut ids = Vec::new();
    for (&bit_id, &count) in fingerprint.inner.nonzero_elements() {
        let count = usize::try_from(count).map_err(|_| {
            PyValueError::new_err("Topological Torsion count cannot be expanded as ids")
        })?;
        ids.extend(std::iter::repeat_n(bit_id, count));
    }
    Ok(ids)
}

#[pymodule]
fn cosmolkit(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    add_public_enums(m)?;
    add_batch_validation_error_class(m)?;
    add_inchi_error_classes(m)?;
    m.add_class::<Molecule>()?;
    m.add_class::<PyAlignmentAtomMap>()?;
    m.add_class::<PyAlignmentParameters>()?;
    m.add_class::<PyBestAlignmentParameters>()?;
    m.add_class::<PyAllConformerRmsdParameters>()?;
    m.add_class::<PyCoordinateRmsdParameters>()?;
    m.add_class::<PyConformerAlignmentParameters>()?;
    m.add_class::<PyAlignmentTransform>()?;
    m.add_class::<PyAlignmentResult>()?;
    m.add_class::<PyConformerRmsd>()?;
    m.add_class::<PyConformerAlignmentReport>()?;
    m.add_class::<PyStereoisomerOptions>()?;
    m.add_class::<PyStereoisomerIterator>()?;
    m.add_class::<PyPotentialStereoInfo>()?;
    m.add_class::<PyPotentialStereoAnalysis>()?;
    m.add_class::<PyEmbedParameters>()?;
    m.add_class::<BioStructure>()?;
    m.add_class::<PyMmcifOutputGroups>()?;
    m.add_class::<PyMmcifWriteOptions>()?;
    m.add_class::<StructureModel>()?;
    m.add_class::<StructureChain>()?;
    m.add_class::<StructureResidue>()?;
    m.add_class::<StructureAtom>()?;
    m.add_class::<StructureEntity>()?;
    m.add_class::<Protein>()?;
    m.add_class::<ProteinChain>()?;
    m.add_class::<ProteinResidue>()?;
    m.add_class::<PyElementInfo>()?;
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
    m.add_class::<PySparseCountFingerprint>()?;
    m.add_class::<PySparseBitFingerprint>()?;
    m.add_class::<PyAdditionalOutput>()?;
    m.add_class::<PyAtomPairsParameters>()?;
    m.add_class::<PyAtomPairAtomInvariantsGenerator>()?;
    m.add_class::<PyTopologicalTorsionFingerprintOptions>()?;
    m.add_class::<PyTopologicalTorsionFingerprintGenerator>()?;
    m.add_class::<Fingerprint>()?;
    m.add_class::<LayeredFingerprintResult>()?;
    m.add_class::<AtomPairFingerprintResult>()?;
    m.add_class::<MorganAdditionalOutput>()?;
    m.add_class::<MorganFingerprintResult>()?;
    m.add_class::<TopologicalFingerprintResult>()?;
    m.add_class::<EmbedMoleculeResult>()?;
    m.add_class::<EmbedMultipleConfsResult>()?;
    m.add_class::<UffOptimizeMoleculeResult>()?;
    m.add_class::<UffOptimizeMoleculeConfResult>()?;
    m.add_class::<UffOptimizeMoleculeConfsResult>()?;
    m.add_class::<MmffOptimizeMoleculeResult>()?;
    m.add_class::<MmffOptimizeMoleculeConfResult>()?;
    m.add_class::<MmffOptimizeMoleculeConfsResult>()?;
    m.add_class::<SubstructMatchResult>()?;
    // RDKit✔️✔️:   python::class_<RDKit::SmartsParserParams, boost::noncopyable>(
    // RDKit✔️✔️:       "SmartsParserParams", "Parameters controlling SMARTS Parsing")
    // RDKit✔️✔️:       .def_readwrite("debugParse", &RDKit::SmartsParserParams::debugParse,
    // RDKit✔️✔️:                      "controls the amount of debugging information produced")
    // RDKit✔️✔️:       .def_readwrite("parseName", &RDKit::SmartsParserParams::parseName,
    // RDKit✔️✔️:                      "controls whether or not the molecule name is also parsed")
    // RDKit✔️✔️:       .def_readwrite("allowCXSMILES",
    // RDKit✔️✔️:                      &RDKit::SmartsParserParams::allowCXSMILES,
    // RDKit✔️✔️:                      "controls whether or not the CXSMILES extensions are parsed")
    // RDKit✔️✔️:       .def_readwrite("strictCXSMILES",
    // RDKit✔️✔️:                      &RDKit::SmartsParserParams::strictCXSMILES,
    // RDKit✔️✔️:                      "controls whether or not problems in CXSMILES parsing "
    // RDKit✔️✔️:                      "causes molecule parsing to fail")
    // RDKit✔️✔️:       .def_readwrite("mergeHs", &RDKit::SmartsParserParams::mergeHs,
    // RDKit✔️✔️:                      "toggles merging H atoms in the SMARTS into neighboring atoms")
    // RDKit✔️✔️:       .def("__setattr__", &safeSetattr);
    m.add_class::<SmartsParserParams>()?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(mol_to_binary, m)?)?;
    m.add_function(wrap_pyfunction!(mol_from_binary, m)?)?;
    m.setattr(
        "_rebuild_molecule_from_pickle",
        wrap_pyfunction!(_rebuild_molecule_from_pickle, m)?,
    )?;
    // RDKit✔️✔️:   python::def("MolFromSmarts", RDKit::MolFromSmarts,
    // RDKit✔️✔️:               (python::arg("SMARTS"), python::arg("mergeHs") = false,
    // RDKit✔️✔️:                python::arg("replacements") = python::dict()),
    // RDKit✔️✔️:               docString.c_str(),
    // RDKit✔️✔️:               python::return_value_policy<python::manage_new_object>());
    // RDKit✔️✔️:   python::def("MolFromSmarts", MolFromSmartsHelper,
    // RDKit✔️✔️:               (python::arg("SMARTS"), python::arg("params")), docString.c_str(),
    // RDKit✔️✔️:               python::return_value_policy<python::manage_new_object>());
    // Complexity review: PyO3 registers two constant-time call shims under
    // distinct snake-case names because Rust does not overload functions.
    // Both shims invoke the same core compiler and return the same Molecule.
    m.add_function(wrap_pyfunction!(parse_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(parse_smarts_with_params, m)?)?;
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
    m.add_function(wrap_pyfunction!(calc_mol_wt, m)?)?;
    m.add_function(wrap_pyfunction!(calc_exact_mol_wt, m)?)?;
    m.add_function(wrap_pyfunction!(calc_mol_formula, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_hbd, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_hba, m)?)?;
    m.add_function(wrap_pyfunction!(calc_fraction_csp3, m)?)?;
    m.add_function(wrap_pyfunction!(calc_crippen_descriptors, m)?)?;
    m.add_function(wrap_pyfunction!(calc_tpsa, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_aromatic_rings, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_rotatable_bonds, m)?)?;
    m.add_function(wrap_pyfunction!(calc_qed, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_0, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_1, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_nv, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_nn, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_0v, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_1v, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_2v, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_3v, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_4v, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_0n, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_1n, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_2n, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_3n, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_4n, m)?)?;
    m.add_function(wrap_pyfunction!(calc_hall_kier_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(
        calc_hall_kier_alpha_with_contributions,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(calc_kappa_1, m)?)?;
    m.add_function(wrap_pyfunction!(calc_kappa_2, m)?)?;
    m.add_function(wrap_pyfunction!(calc_kappa_3, m)?)?;
    m.add_function(wrap_pyfunction!(calc_phi, m)?)?;
    m.add_function(wrap_pyfunction!(calc_lipinski_hba, m)?)?;
    m.add_function(wrap_pyfunction!(calc_lipinski_hbd, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_heteroatoms, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_amide_bonds, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_heavy_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_rings, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_heterocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_saturated_rings, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_aliphatic_rings, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_aromatic_heterocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_aromatic_carbocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_aliphatic_heterocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_aliphatic_carbocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_saturated_heterocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_saturated_carbocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_spiro_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_bridgehead_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_atom_stereo_centers, m)?)?;
    m.add_function(wrap_pyfunction!(
        calc_num_unspecified_atom_stereo_centers,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(calc_mqns, m)?)?;
    m.add_function(wrap_pyfunction!(calc_labute_asa, m)?)?;
    m.add_function(wrap_pyfunction!(calc_labute_asa_contributions, m)?)?;
    m.add_function(wrap_pyfunction!(calc_slogp_vsa, m)?)?;
    m.add_function(wrap_pyfunction!(calc_smr_vsa, m)?)?;
    m.add_function(wrap_pyfunction!(calc_slogp_vsa_1, m)?)?;
    m.add_function(wrap_pyfunction!(calc_slogp_vsa_2, m)?)?;
    m.add_function(wrap_pyfunction!(calc_slogp_vsa_3, m)?)?;
    m.add_function(wrap_pyfunction!(calc_slogp_vsa_4, m)?)?;
    m.add_function(wrap_pyfunction!(calc_slogp_vsa_5, m)?)?;
    m.add_function(wrap_pyfunction!(calc_slogp_vsa_6, m)?)?;
    m.add_function(wrap_pyfunction!(calc_slogp_vsa_7, m)?)?;
    m.add_function(wrap_pyfunction!(calc_slogp_vsa_8, m)?)?;
    m.add_function(wrap_pyfunction!(calc_slogp_vsa_9, m)?)?;
    m.add_function(wrap_pyfunction!(calc_slogp_vsa_10, m)?)?;
    m.add_function(wrap_pyfunction!(calc_slogp_vsa_11, m)?)?;
    m.add_function(wrap_pyfunction!(calc_slogp_vsa_12, m)?)?;
    m.add_function(wrap_pyfunction!(calc_smr_vsa_1, m)?)?;
    m.add_function(wrap_pyfunction!(calc_smr_vsa_2, m)?)?;
    m.add_function(wrap_pyfunction!(calc_smr_vsa_3, m)?)?;
    m.add_function(wrap_pyfunction!(calc_smr_vsa_4, m)?)?;
    m.add_function(wrap_pyfunction!(calc_smr_vsa_5, m)?)?;
    m.add_function(wrap_pyfunction!(calc_smr_vsa_6, m)?)?;
    m.add_function(wrap_pyfunction!(calc_smr_vsa_7, m)?)?;
    m.add_function(wrap_pyfunction!(calc_smr_vsa_8, m)?)?;
    m.add_function(wrap_pyfunction!(calc_smr_vsa_9, m)?)?;
    m.add_function(wrap_pyfunction!(calc_smr_vsa_10, m)?)?;
    m.add_function(wrap_pyfunction!(element_from_symbol, m)?)?;
    m.add_function(wrap_pyfunction!(get_element_info, m)?)?;
    m.add_function(wrap_pyfunction!(find_tabulated_residue, m)?)?;
    m.add_function(wrap_pyfunction!(find_tabulated_residue_idx, m)?)?;
    m.add_function(wrap_pyfunction!(get_residue_info, m)?)?;
    m.add_function(wrap_pyfunction!(residue_code_from_name, m)?)?;
    m.add_function(wrap_pyfunction!(expand_one_letter, m)?)?;
    m.add_function(wrap_pyfunction!(expand_protein_one_letter, m)?)?;
    m.add_function(wrap_pyfunction!(expand_one_letter_sequence, m)?)?;
    m.add_function(wrap_pyfunction!(expand_protein_one_letter_string, m)?)?;
    m.add_function(wrap_pyfunction!(inchi_to_key, m)?)?;
    m.add_function(wrap_pyfunction!(get_topological_torsion_generator, m)?)?;
    m.add_function(wrap_pyfunction!(
        topological_torsion_generator_from_json,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(get_topological_torsion_fingerprint, m)?)?;
    m.add_function(wrap_pyfunction!(
        get_hashed_topological_torsion_fingerprint,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(
        get_hashed_topological_torsion_fingerprint_as_bit_vect,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(get_atom_code, m)?)?;
    m.add_function(wrap_pyfunction!(explain_atom_code, m)?)?;
    m.add_function(wrap_pyfunction!(py_score_path, m)?)?;
    m.add_function(wrap_pyfunction!(explain_path_score, m)?)?;
    m.add_function(wrap_pyfunction!(
        get_topological_torsion_fingerprint_as_ids,
        m
    )?)?;
    confseq_py::add_confseq_module(m)?;
    Ok(())
}

#[cfg(feature = "stubgen")]
define_stub_info_gatherer!(stub_info);
