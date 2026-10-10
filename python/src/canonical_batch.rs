//! Ordered batch and supplier projections of the canonical public facade.
use crate::canonical_batch_fingerprint_values::BatchFingerprintOutput;
use crate::canonical_layered::LayeredFingerprintResult;
use crate::canonical_sdf::{SdfRecord, sdf_pyerr as sdf_error};
use crate::canonical_values::{
    SparseBitFingerprint, SparseCountFingerprint, SparseCountFingerprint32,
};
use crate::drawing_binding::Molecule;
use crate::fingerprint_numpy::FingerprintBatch;
use crate::text_path::TextPath;
use ::cosmolkit as ck;
use pyo3::exceptions::{PyIndexError, PyNotImplementedError, PyTypeError, PyValueError};
use pyo3::prelude::*;

// Scalar and selection branches share runtime dispatch but have distinct types.
#[cfg(all(test, feature = "stubgen"))]
#[test]
fn index_stubs_distinguish_scalar_and_selection_results() {
    pyo3::Python::initialize();
    let info = crate::stub_info().unwrap();
    let text = info.modules["cosmolkit"].format_with_config(info.config.use_type_statement);
    for (name, scalar) in [
        ("MoleculeBatch", "Molecule | None"),
        ("SdfDataset", "SdfRecord"),
    ] {
        let start = text.find(&format!("class {name}:")).unwrap();
        let class = text[start..].split("\nclass ").next().unwrap();
        assert_eq!(class.matches("def __getitem__").count(), 2, "{class}");
        assert!(class.contains(&format!(") -> {scalar}:")), "{class}");
        assert!(class.contains(") -> MoleculeBatch:"), "{class}");
        assert!(class.contains("typing.Sequence[builtins.int]"), "{class}");
    }
}

#[cfg(feature = "stubgen")]
pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! {
        r#"
        class MoleculeBatch:
            @overload
            def __getitem__(self, key: int) -> Molecule | None: ...
            @overload
            def __getitem__(self, key: slice | typing.Sequence[int] | typing.Sequence[bool]) -> MoleculeBatch: ...
        "#
    }
}
#[cfg(feature = "stubgen")]
pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! {
        r#"
        class SdfDataset:
            @overload
            def __getitem__(self, key: int) -> SdfRecord: ...
            @overload
            def __getitem__(self, key: slice | typing.Sequence[int] | typing.Sequence[bool]) -> MoleculeBatch: ...
        "#
    }
}
use pyo3::types::{PyAny, PyBool, PyDict, PySlice, PySliceMethods, PyType};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
pyo3::create_exception!(cosmolkit, BatchImageError, PyValueError);
fn image_path_error(
    py: Python<'_>,
    source: crate::user_path::ImagePathError<ck::BatchValidationError>,
) -> PyErr {
    match source {
        crate::user_path::ImagePathError::Directory(source)
        | crate::user_path::ImagePathError::Report(source) => {
            PyValueError::new_err(source.to_string())
        }
        crate::user_path::ImagePathError::Export(source) => batch_error(py, source),
        crate::user_path::ImagePathError::ReportWrite(source) => {
            // COSMolKit❗✔️: pinned d892ec3 python/src/lib.rs:552-553.
            //     fs::write(&expanded_path, content)
            //         .map_err(|err| PyValueError::new_err(format!("write error report failed: {err}")))
            // The sole writer retains the actual IO cause in its record error.
            if let Some(cause) = source
                .record_errors
                .first()
                .and_then(std::error::Error::source)
            {
                crate::canonical_values::annotate(
                    py,
                    PyValueError::new_err(format!("write error report failed: {cause}")),
                    "batch",
                    "ReportWrite",
                    &source,
                )
            } else {
                batch_error(py, source)
            }
        }
    }
}

pub(crate) fn batch_image_error(py: Python<'_>, source: &ck::BatchImageError) -> PyErr {
    let kind = match source {
        ck::BatchImageError::InvalidFormat(_) => "InvalidFormat",
        ck::BatchImageError::Write(_) => "Write",
    };
    let error = crate::canonical_values::annotate(
        py,
        BatchImageError::new_err(source.to_string()),
        "batch",
        kind,
        source,
    );
    if let ck::BatchImageError::InvalidFormat(format) = source {
        if let Err(error) = error.value(py).setattr("format", format) {
            return error;
        }
    }
    error
}

pub(crate) fn io_error(py: Python<'_>, source: ck::MolecularIoError) -> PyErr {
    crate::canonical_molecular_io::error_pyerr(py, source)
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct SdfRecordMetadata {
    pub(crate) inner: ck::SdfRecordMetadata,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SdfRecordMetadata {
    fn index(&self) -> usize {
        self.inner.index()
    }
    fn byte_offset(&self) -> u64 {
        self.inner.byte_offset()
    }
    fn byte_len(&self) -> u64 {
        self.inner.byte_len()
    }
    fn byte_range(&self) -> (u64, u64) {
        self.inner.byte_range()
    }
    fn line_range(&self) -> (usize, usize) {
        self.inner.line_range()
    }
    fn title(&self) -> Option<String> {
        self.inner.title().map(str::to_owned)
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

fn add_batch_validation_error_class(module: &Bound<'_, PyModule>) -> PyResult<()> {
    // COSMolKit❗✔️: pinned d892ec3 python/src/lib.rs::add_batch_validation_error_class; canonical projection preserves source branches/defaults.
    // fn add_batch_validation_error_class(m: &Bound<'_, PyModule>) -> PyResult<()> {
    //     let py = m.py();
    //     let globals = PyDict::new(py);
    //     globals.set_item("ValueError", py.get_type::<PyValueError>())?;
    //     let code = r#"
    // class BatchValidationError(ValueError):
    //     __module__ = "cosmolkit"
    //
    //     def __init__(self, message, error_count=0, reason=None, record_errors=None):
    //         super().__init__(message)
    //         self.error_count = int(error_count)
    //         self.reason = reason
    //         self._errors = list(record_errors or [])
    //
    //     def errors(self):
    //         return list(self._errors)
    // "#;
    //     py.import("builtins")?
    //         .getattr("exec")?
    //         .call1((code, &globals))?;
    //     let cls = globals
    //         .get_item("BatchValidationError")?
    //         .ok_or_else(|| PyValueError::new_err("failed to create BatchValidationError class"))?;
    //     m.add("BatchValidationError", cls)?;
    //     Ok(())
    // }

    let py = module.py();
    let globals = pyo3::types::PyDict::new(py);
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
    module.add("BatchValidationError", cls)?;
    Ok(())
}
pub(crate) fn batch_error(
    py: Python<'_>,
    source: impl std::borrow::Borrow<ck::BatchValidationError>,
) -> PyErr {
    let source = source.borrow();
    // COSMolKit❗✔️: pinned d892ec3 python/src/lib.rs::batch_validation_pyerr; canonical projection preserves source branches/defaults.
    // fn batch_validation_pyerr(error: cosmolkit_core::BatchValidationError) -> PyErr {
    //     let message = error.to_string();
    //     Python::attach(|py| {
    //         let error_count = error.errors;
    //         let reason = error.reason.map(|value| value.to_string());
    //         let record_errors: Vec<PyBatchError> =
    //             error.record_errors.into_iter().map(Into::into).collect();
    //         match (|| -> PyResult<Bound<'_, PyAny>> {
    //             let cls = py.import("cosmolkit")?.getattr("BatchValidationError")?;
    //             cls.call1((message, error_count, reason, record_errors))
    //         })() {
    //             Ok(instance) => PyErr::from_value(instance),
    //             Err(error) => error,
    //         }
    //     })
    // }

    let instance = (|| -> PyResult<Bound<'_, PyAny>> {
        let cls = py.import("cosmolkit")?.getattr("BatchValidationError")?;
        let record_errors: Vec<BatchError> = source
            .record_errors
            .iter()
            .cloned()
            .map(|inner| BatchError { inner })
            .collect();
        cls.call1((source.to_string(), source.errors, py.None(), record_errors))
    })();
    match instance {
        Ok(instance) => {
            let error = PyErr::from_value(instance);
            error.set_cause(
                py,
                source
                    .record_errors
                    .first()
                    .map(|cause| crate::canonical_values::source_pyerr(py, cause)),
            );
            error
        }
        Err(error) => error,
    }
}
pub(crate) fn error_mode(value: Option<&Bound<'_, PyAny>>) -> PyResult<ck::BatchErrorMode> {
    let Some(value) = value else {
        return Ok(ck::BatchErrorMode::Strict);
    };
    if let Ok(text) = value.extract::<String>() {
        return match text.to_ascii_lowercase().as_str() {
            "raise" => Ok(ck::BatchErrorMode::Strict),
            "keep" => Ok(ck::BatchErrorMode::KeepErrors),
            _ => Err(PyValueError::new_err(format!(
                "unsupported errors mode '{text}', expected one of: raise, keep"
            ))),
        };
    }
    match value.extract::<i64>()? {
        1 => Ok(ck::BatchErrorMode::Strict),
        2 => Ok(ck::BatchErrorMode::KeepErrors),
        code => Err(PyValueError::new_err(format!(
            "unsupported errors mode code {code}, expected BatchErrorMode.RAISE or KEEP"
        ))),
    }
}
pub(crate) fn n_jobs(value: Option<usize>) -> PyResult<Option<usize>> {
    if value == Some(0) {
        Err(PyValueError::new_err("n_jobs must be >= 1"))
    } else {
        Ok(value)
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct BatchError {
    inner: ck::BatchError,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BatchError {
    fn index(&self) -> usize {
        self.inner.index
    }
    fn operation(&self) -> String {
        self.inner.operation.to_owned()
    }
    fn message(&self) -> String {
        self.inner.message.clone()
    }
    /// Return a fresh dictionary with an integer input index and string details.
    #[gen_stub(override_return_type(type_repr = "dict[str, int | str]"))]
    fn as_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        result.set_item("index", self.inner.index)?;
        result.set_item("operation", self.inner.operation)?;
        result.set_item("message", &self.inner.message)?;
        Ok(result)
    }
    fn __repr__(&self) -> String {
        format!(
            "BatchError(index={}, operation='{}', message='{}')",
            self.inner.index, self.inner.operation, self.inner.message
        )
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct MoleculeBatch {
    pub(crate) inner: ck::MoleculeBatch,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MoleculeBatch {
    #[staticmethod]
    #[pyo3(signature=(text,read,errors=None,n_jobs=None))]
    fn from_sdf_records_with_params(
        py: Python<'_>,
        text: &str,
        read: &crate::canonical_sdf::SdfReadParams,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
    ) -> PyResult<Self> {
        ck::MoleculeBatch::from_sdf_records_with_params(
            text,
            &read.inner,
            error_mode(errors)?,
            n_jobs,
        )
        .map(|inner| Self { inner })
        .map_err(|e| batch_error(py, e))
    }
    #[staticmethod]
    #[pyo3(signature=(path,read,errors=None,n_jobs=None,progress_bar=false))]
    fn read_sdf_with_params(
        py: Python<'_>,
        path: TextPath,
        read: &crate::canonical_sdf::SdfReadParams,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: bool,
    ) -> PyResult<Self> {
        ck::MoleculeBatch::read_sdf_with_params(
            path.as_str(),
            &read.inner,
            error_mode(errors)?,
            n_jobs,
            progress_bar,
        )
        .map(|inner| Self { inner })
        .map_err(|e| batch_error(py, e))
    }
    #[staticmethod]
    fn from_dataset_indices(
        py: Python<'_>,
        dataset: &SdfDataset,
        indices: Vec<usize>,
        mode: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        ck::MoleculeBatch::from_dataset_indices(&dataset.inner, &indices, error_mode(Some(mode))?)
            .map(|inner| Self { inner })
            .map_err(|e| batch_error(py, e))
    }
    #[staticmethod]
    fn from_records(
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Sequence[Molecule | BatchError]",imports=("typing")))]
        records: Vec<Py<PyAny>>,
        mode: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let records = records
            .iter()
            .map(|record| {
                if let Ok(molecule) = record.extract::<PyRef<'_, Molecule>>(py) {
                    Ok(ck::BatchRecord::Molecule(molecule.inner.clone()))
                } else if let Ok(error) = record.extract::<PyRef<'_, BatchError>>(py) {
                    Ok(ck::BatchRecord::Error(error.inner.clone()))
                } else {
                    Err(pyo3::exceptions::PyTypeError::new_err(
                        "records must contain Molecule or BatchError values",
                    ))
                }
            })
            .collect::<PyResult<Vec<_>>>()?;
        ck::MoleculeBatch::from_records(records, error_mode(Some(mode))?)
            .map(|inner| Self { inner })
            .map_err(|e| batch_error(py, e))
    }
    #[gen_stub(override_return_type(type_repr = "Molecule | BatchError | None"))]
    fn get(&self, py: Python<'_>, index: usize) -> PyResult<Option<Py<PyAny>>> {
        self.inner
            .get(index)
            .map(|record| batch_record_object(py, record))
            .transpose()
    }
    #[gen_stub(override_return_type(type_repr = "list[Molecule | BatchError]"))]
    fn records(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
        self.inner
            .records()
            .iter()
            .map(|record| batch_record_object(py, record))
            .collect()
    }
    #[pyo3(signature=(path,params,report_path))]
    fn write_sdf_with_params(
        &self,
        py: Python<'_>,
        path: TextPath,
        params: &crate::canonical_batch_params::BatchExportParams,
        report_path: Option<TextPath>,
    ) -> PyResult<BatchExportReport> {
        self.inner
            .write_sdf_with_params(
                path.as_str(),
                &params.inner,
                report_path.as_ref().map(TextPath::as_str),
            )
            .map(|inner| BatchExportReport { inner })
            .map_err(|e| batch_error(py, e))
    }
    #[pyo3(signature=(directory,params,filenames,report_path))]
    fn write_sdf_files_with_params(
        &self,
        py: Python<'_>,
        directory: TextPath,
        params: &crate::canonical_batch_params::BatchExportParams,
        filenames: Option<Vec<Option<String>>>,
        report_path: Option<TextPath>,
    ) -> PyResult<BatchExportReport> {
        self.inner
            .write_sdf_files_with_params(
                directory.as_str(),
                &params.inner,
                filenames.as_deref(),
                report_path.as_ref().map(TextPath::as_str),
            )
            .map(|inner| BatchExportReport { inner })
            .map_err(|e| batch_error(py, e))
    }
    fn fingerprint_topological_torsion_list(&self, py: Python<'_>) -> PyResult<FingerprintBatch> {
        self.inner
            .fingerprint_topological_torsion_list()
            .map_err(|e| batch_error(py, e))
            .and_then(|values| FingerprintBatch::from_values(py, values))
    }
    fn fingerprint_topological_torsion_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_fingerprint_values::TopologicalTorsionFingerprintParams,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<FingerprintBatch> {
        params
            .execute(py, |execution| {
                self.inner
                    .fingerprint_topological_torsion_list_with_params(&options.inner, execution)
            })
            .and_then(|values| FingerprintBatch::from_values(py, values))
    }
    fn fingerprint_atom_pair_list(&self, py: Python<'_>) -> PyResult<FingerprintBatch> {
        self.inner
            .fingerprint_atom_pair_list()
            .map_err(|source| batch_error(py, source))
            .and_then(|values| FingerprintBatch::from_values(py, values))
    }
    fn fingerprint_atom_pair_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_fingerprint_values::AtomPairFingerprintParams,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<FingerprintBatch> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn fingerprint_atom_pair_list(
        //         &self,
        //         n_bits: usize,
        //         min_distance: u32,
        //         max_distance: u32,
        //         use_2d: bool,
        //         include_chirality: bool,
        //         count_simulation: bool,
        //         count_bounds: Option<Vec<u32>>,
        //         num_bits_per_feature: u32,
        //         from_atoms: Option<Vec<usize>>,
        //         ignore_atoms: Option<Vec<usize>>,
        //         conformer_id: i32,
        //         custom_atom_invariants: Option<Vec<u32>>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Vec<Option<Fingerprint>>> {
        //         let params = make_atom_pair_fingerprint_params(
        //             n_bits,
        //             min_distance,
        //             max_distance,
        //             use_2d,
        //             include_chirality,
        //             count_simulation,
        //             count_bounds,
        //             num_bits_per_feature,
        //             from_atoms,
        //             ignore_atoms,
        //             conformer_id,
        //             custom_atom_invariants,
        //             false,
        //         );
        //         self.inner
        //             .atom_pair_fingerprint_list_with_options(
        //                 &params,
        //                 validate_n_jobs(n_jobs)?,
        //                 progress_bar,
        //             )
        //             .map(|values| {
        //                 values
        //                     .into_iter()
        //                     .map(|value| value.map(|inner| Fingerprint { inner }))
        //                     .collect()
        //             })
        //             .map_err(batch_validation_pyerr)
        //     }
        params
            .execute(py, |execution| {
                self.inner
                    .fingerprint_atom_pair_list_with_params(&options.inner, execution)
            })
            .and_then(|values| FingerprintBatch::from_values(py, values))
    }

    fn fingerprint_atom_pair_sparse_count_list(
        &self,
        py: Python<'_>,
    ) -> PyResult<Vec<Option<SparseCountFingerprint>>> {
        self.inner
            .fingerprint_atom_pair_sparse_count_list()
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| SparseCountFingerprint { inner }))
                    .collect()
            })
            .map_err(|source| batch_error(py, source))
    }
    fn fingerprint_atom_pair_sparse_count_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_fingerprint_values::AtomPairFingerprintParams,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<SparseCountFingerprint>>> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn fingerprint_atom_pair_sparse_count_list(
        //         &self,
        //         n_bits: usize,
        //         min_distance: u32,
        //         max_distance: u32,
        //         use_2d: bool,
        //         include_chirality: bool,
        //         count_simulation: bool,
        //         count_bounds: Option<Vec<u32>>,
        //         num_bits_per_feature: u32,
        //         from_atoms: Option<Vec<usize>>,
        //         ignore_atoms: Option<Vec<usize>>,
        //         conformer_id: i32,
        //         custom_atom_invariants: Option<Vec<u32>>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Vec<Option<PySparseCountFingerprint>>> {
        //         let params = make_atom_pair_fingerprint_params(
        //             n_bits,
        //             min_distance,
        //             max_distance,
        //             use_2d,
        //             include_chirality,
        //             count_simulation,
        //             count_bounds,
        //             num_bits_per_feature,
        //             from_atoms,
        //             ignore_atoms,
        //             conformer_id,
        //             custom_atom_invariants,
        //             false,
        //         );
        //         self.inner
        //             .atom_pair_sparse_count_fingerprint_list_with_options(
        //                 &params,
        //                 validate_n_jobs(n_jobs)?,
        //                 progress_bar,
        //             )
        //             .map(|values| {
        //                 values
        //                     .into_iter()
        //                     .map(|value| value.map(|inner| PySparseCountFingerprint { inner }))
        //                     .collect()
        //             })
        //             .map_err(batch_validation_pyerr)
        //     }
        params
            .execute(py, |execution| {
                self.inner
                    .fingerprint_atom_pair_sparse_count_list_with_params(&options.inner, execution)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| SparseCountFingerprint { inner }))
                    .collect()
            })
    }

    fn fingerprint_atom_pair_count_list(
        &self,
        py: Python<'_>,
    ) -> PyResult<Vec<Option<SparseCountFingerprint32>>> {
        self.inner
            .fingerprint_atom_pair_count_list()
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| SparseCountFingerprint32 { inner }))
                    .collect()
            })
            .map_err(|source| batch_error(py, source))
    }
    fn fingerprint_atom_pair_count_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_fingerprint_values::AtomPairFingerprintParams,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<SparseCountFingerprint32>>> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn fingerprint_atom_pair_count_list(
        //         &self,
        //         n_bits: usize,
        //         min_distance: u32,
        //         max_distance: u32,
        //         use_2d: bool,
        //         include_chirality: bool,
        //         count_simulation: bool,
        //         count_bounds: Option<Vec<u32>>,
        //         num_bits_per_feature: u32,
        //         from_atoms: Option<Vec<usize>>,
        //         ignore_atoms: Option<Vec<usize>>,
        //         conformer_id: i32,
        //         custom_atom_invariants: Option<Vec<u32>>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Vec<Option<PySparseCountFingerprint>>> {
        //         let params = make_atom_pair_fingerprint_params(
        //             n_bits,
        //             min_distance,
        //             max_distance,
        //             use_2d,
        //             include_chirality,
        //             count_simulation,
        //             count_bounds,
        //             num_bits_per_feature,
        //             from_atoms,
        //             ignore_atoms,
        //             conformer_id,
        //             custom_atom_invariants,
        //             false,
        //         );
        //         self.inner
        //             .atom_pair_count_fingerprint_list_with_options(
        //                 &params,
        //                 validate_n_jobs(n_jobs)?,
        //                 progress_bar,
        //             )
        //             .map(|values| {
        //                 values
        //                     .into_iter()
        //                     .map(|value| value.map(|inner| PySparseCountFingerprint { inner }))
        //                     .collect()
        //             })
        //             .map_err(batch_validation_pyerr)
        //     }
        params
            .execute(py, |execution| {
                self.inner
                    .fingerprint_atom_pair_count_list_with_params(&options.inner, execution)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| SparseCountFingerprint32 { inner }))
                    .collect()
            })
    }

    fn fingerprint_atom_pair_sparse_bits_list(
        &self,
        py: Python<'_>,
    ) -> PyResult<Vec<Option<SparseBitFingerprint>>> {
        self.inner
            .fingerprint_atom_pair_sparse_bits_list()
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| SparseBitFingerprint { inner }))
                    .collect()
            })
            .map_err(|source| batch_error(py, source))
    }
    fn fingerprint_atom_pair_sparse_bits_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_fingerprint_values::AtomPairFingerprintParams,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<SparseBitFingerprint>>> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn fingerprint_atom_pair_sparse_bits_list(
        //         &self,
        //         n_bits: usize,
        //         min_distance: u32,
        //         max_distance: u32,
        //         use_2d: bool,
        //         include_chirality: bool,
        //         count_simulation: bool,
        //         count_bounds: Option<Vec<u32>>,
        //         num_bits_per_feature: u32,
        //         from_atoms: Option<Vec<usize>>,
        //         ignore_atoms: Option<Vec<usize>>,
        //         conformer_id: i32,
        //         custom_atom_invariants: Option<Vec<u32>>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Vec<Option<PySparseBitFingerprint>>> {
        //         let params = make_atom_pair_fingerprint_params(
        //             n_bits,
        //             min_distance,
        //             max_distance,
        //             use_2d,
        //             include_chirality,
        //             count_simulation,
        //             count_bounds,
        //             num_bits_per_feature,
        //             from_atoms,
        //             ignore_atoms,
        //             conformer_id,
        //             custom_atom_invariants,
        //             false,
        //         );
        //         self.inner
        //             .atom_pair_sparse_bit_fingerprint_list_with_options(
        //                 &params,
        //                 validate_n_jobs(n_jobs)?,
        //                 progress_bar,
        //             )
        //             .map(|values| {
        //                 values
        //                     .into_iter()
        //                     .map(|value| value.map(|inner| PySparseBitFingerprint { inner }))
        //                     .collect()
        //             })
        //             .map_err(batch_validation_pyerr)
        //     }
        params
            .execute(py, |execution| {
                self.inner
                    .fingerprint_atom_pair_sparse_bits_list_with_params(&options.inner, execution)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| SparseBitFingerprint { inner }))
                    .collect()
            })
    }

    fn fingerprint_atom_pair_with_output_list(
        &self,
        py: Python<'_>,
    ) -> PyResult<Vec<Option<BatchFingerprintOutput>>> {
        self.inner
            .fingerprint_atom_pair_with_output_list()
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| BatchFingerprintOutput { inner }))
                    .collect()
            })
            .map_err(|source| batch_error(py, source))
    }
    #[pyo3(signature=(options, collect_additional_output=true, params=None))]
    fn fingerprint_atom_pair_with_output_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_fingerprint_values::AtomPairFingerprintParams,
        collect_additional_output: bool,
        params: Option<&crate::canonical_batch_params::BatchQueryParams>,
    ) -> PyResult<Vec<Option<BatchFingerprintOutput>>> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn fingerprint_atom_pair_with_output_list(
        //         &self,
        //         n_bits: usize,
        //         min_distance: u32,
        //         max_distance: u32,
        //         use_2d: bool,
        //         include_chirality: bool,
        //         count_simulation: bool,
        //         count_bounds: Option<Vec<u32>>,
        //         num_bits_per_feature: u32,
        //         from_atoms: Option<Vec<usize>>,
        //         ignore_atoms: Option<Vec<usize>>,
        //         conformer_id: i32,
        //         custom_atom_invariants: Option<Vec<u32>>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Vec<Option<AtomPairFingerprintResult>>> {
        //         let params = make_atom_pair_fingerprint_params(
        //             n_bits,
        //             min_distance,
        //             max_distance,
        //             use_2d,
        //             include_chirality,
        //             count_simulation,
        //             count_bounds,
        //             num_bits_per_feature,
        //             from_atoms,
        //             ignore_atoms,
        //             conformer_id,
        //             custom_atom_invariants,
        //             true,
        //         );
        //         self.inner
        //             .atom_pair_fingerprint_with_output_list_with_options(
        //                 &params,
        //                 validate_n_jobs(n_jobs)?,
        //                 progress_bar,
        //             )
        //             .map(|values| {
        //                 values
        //                     .into_iter()
        //                     .map(|value| value.map(Into::into))
        //                     .collect()
        //             })
        //             .map_err(batch_validation_pyerr)
        //     }
        let defaults = crate::canonical_batch_params::BatchQueryParams::default();
        params
            .unwrap_or(&defaults)
            .execute(py, |execution| {
                self.inner
                    .fingerprint_atom_pair_with_output_list_with_params(
                        &options.inner,
                        collect_additional_output,
                        execution,
                    )
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| BatchFingerprintOutput { inner }))
                    .collect()
            })
    }

    fn fingerprint_layered_list(&self, py: Python<'_>) -> PyResult<FingerprintBatch> {
        self.inner
            .fingerprint_layered_list()
            .map_err(|source| batch_error(py, source))
            .and_then(|values| FingerprintBatch::from_values(py, values))
    }
    fn fingerprint_layered_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_layered::LayeredFingerprintParams,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<FingerprintBatch> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn fingerprint_layered_list(
        //         &self,
        //         layers: u32,
        //         min_path: u32,
        //         max_path: u32,
        //         fp_size: u32,
        //         atom_counts: Option<Vec<u32>>,
        //         set_only_bits: Option<&Fingerprint>,
        //         branched_paths: bool,
        //         from_atoms: Option<Vec<u32>>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Vec<Option<Fingerprint>>> {
        //         let params = make_layered_fingerprint_params(
        //             layers,
        //             min_path,
        //             max_path,
        //             fp_size,
        //             atom_counts,
        //             set_only_bits,
        //             branched_paths,
        //             from_atoms,
        //         );
        //         self.inner
        //             .layered_fingerprint_list_with_options(&params, validate_n_jobs(n_jobs)?, progress_bar)
        //             .map(|values| {
        //                 values
        //                     .into_iter()
        //                     .map(|value| value.map(|inner| Fingerprint { inner }))
        //                     .collect()
        //             })
        //             .map_err(batch_validation_pyerr)
        //     }
        params
            .execute(py, |execution| {
                self.inner
                    .fingerprint_layered_list_with_params(&options.inner, execution)
            })
            .and_then(|values| FingerprintBatch::from_values(py, values))
    }

    fn fingerprint_layered_with_output_list(
        &self,
        py: Python<'_>,
    ) -> PyResult<Vec<Option<LayeredFingerprintResult>>> {
        self.inner
            .fingerprint_layered_with_output_list()
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| LayeredFingerprintResult { inner }))
                    .collect()
            })
            .map_err(|source| batch_error(py, source))
    }
    fn fingerprint_layered_with_output_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_layered::LayeredFingerprintParams,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<LayeredFingerprintResult>>> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn fingerprint_layered_with_output_list(
        //         &self,
        //         layers: u32,
        //         min_path: u32,
        //         max_path: u32,
        //         fp_size: u32,
        //         atom_counts: Option<Vec<u32>>,
        //         set_only_bits: Option<&Fingerprint>,
        //         branched_paths: bool,
        //         from_atoms: Option<Vec<u32>>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Vec<Option<LayeredFingerprintResult>>> {
        //         let params = make_layered_fingerprint_params(
        //             layers,
        //             min_path,
        //             max_path,
        //             fp_size,
        //             atom_counts,
        //             set_only_bits,
        //             branched_paths,
        //             from_atoms,
        //         );
        //         self.inner
        //             .layered_fingerprint_with_output_list_with_options(
        //                 &params,
        //                 validate_n_jobs(n_jobs)?,
        //                 progress_bar,
        //             )
        //             .map(|values| {
        //                 values
        //                     .into_iter()
        //                     .map(|value| value.map(Into::into))
        //                     .collect()
        //             })
        //             .map_err(batch_validation_pyerr)
        //     }
        params
            .execute(py, |execution| {
                self.inner
                    .fingerprint_layered_with_output_list_with_params(&options.inner, execution)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| LayeredFingerprintResult { inner }))
                    .collect()
            })
    }

    fn fingerprint_pattern_list(&self, py: Python<'_>) -> PyResult<FingerprintBatch> {
        self.inner
            .fingerprint_pattern_list()
            .map_err(|source| batch_error(py, source))
            .and_then(|values| FingerprintBatch::from_values(py, values))
    }
    fn fingerprint_pattern_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_pattern::PatternFingerprintParams,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<FingerprintBatch> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn pattern_fingerprint_list(
        //         &self,
        //         n_bits: usize,
        //         tautomeric: bool,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Vec<Option<Fingerprint>>> {
        //         let params = cosmolkit_core::PatternFingerprintParams { n_bits, tautomeric };
        //         self.inner
        //             .pattern_fingerprint_list_with_options(&params, validate_n_jobs(n_jobs)?, progress_bar)
        //             .map(|values| {
        //                 values
        //                     .into_iter()
        //                     .map(|value| value.map(|inner| Fingerprint { inner }))
        //                     .collect()
        //             })
        //             .map_err(batch_validation_pyerr)
        //     }
        params
            .execute(py, |execution| {
                self.inner
                    .fingerprint_pattern_list_with_params(&options.inner, execution)
            })
            .and_then(|values| FingerprintBatch::from_values(py, values))
    }

    /// Return Morgan fingerprints in input order (default radius 3, 2048 bits).
    ///
    /// Existing failed records produce None; new calculation errors raise
    /// BatchValidationError. The source batch is unchanged.
    /// Accepts MorganFingerprintParams as options and BatchQueryParams as params,
    /// or their fields as keywords. For example:
    /// batch.fingerprint_morgan_list(generator=ck.MorganParams(radius=2), n_jobs=1).
    /// Omitted n_jobs uses the batch setting, or one worker if unset.
    fn fingerprint_morgan_list(&self, py: Python<'_>) -> PyResult<FingerprintBatch> {
        self.inner
            .fingerprint_morgan_list()
            .map_err(|source| batch_error(py, source))
            .and_then(|values| FingerprintBatch::from_values(py, values))
    }
    /// Return ordered Morgan fingerprints using explicit options and execution params.
    /// Existing failed records produce None; new calculation errors raise
    /// BatchValidationError. Does not modify the source batch.
    fn fingerprint_morgan_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_fingerprint_values::MorganFingerprintParams,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<FingerprintBatch> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn fingerprint_morgan_list(
        //         &self,
        //         radius: u32,
        //         n_bits: usize,
        //         include_chirality: bool,
        //         use_bond_types: bool,
        //         count_simulation: bool,
        //         count_bounds: Option<Vec<u32>>,
        //         only_nonzero_invariants: bool,
        //         include_redundant_environments: bool,
        //         from_atoms: Option<Vec<usize>>,
        //         ignore_atoms: Option<Vec<usize>>,
        //         custom_atom_invariants: Option<Vec<u32>>,
        //         custom_bond_invariants: Option<Vec<u32>>,
        //         atom_invariants_generator: Option<&str>,
        //         atom_invariants_include_ring_membership: bool,
        //         bond_invariants_generator: Option<&str>,
        //         bond_invariants_use_bond_types: bool,
        //         bond_invariants_use_chirality: bool,
        //         num_bits_per_feature: u32,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Vec<Option<Fingerprint>>> {
        //         let params = make_morgan_fingerprint_params(
        //             radius,
        //             n_bits,
        //             include_chirality,
        //             use_bond_types,
        //             count_simulation,
        //             count_bounds,
        //             only_nonzero_invariants,
        //             include_redundant_environments,
        //             from_atoms,
        //             ignore_atoms,
        //             custom_atom_invariants,
        //             custom_bond_invariants,
        //             atom_invariants_generator,
        //             atom_invariants_include_ring_membership,
        //             bond_invariants_generator,
        //             bond_invariants_use_bond_types,
        //             bond_invariants_use_chirality,
        //             num_bits_per_feature,
        //             false,
        //         )?;
        //         self.inner
        //             .morgan_fingerprint_list_with_options(&params, validate_n_jobs(n_jobs)?, progress_bar)
        //             .map(|values| {
        //                 values
        //                     .into_iter()
        //                     .map(|value| value.map(|inner| Fingerprint { inner }))
        //                     .collect()
        //             })
        //             .map_err(batch_validation_pyerr)
        //     }
        params
            .execute(py, |execution| {
                self.inner
                    .fingerprint_morgan_list_with_params(&options.inner, execution)
            })
            .and_then(|values| FingerprintBatch::from_values(py, values))
    }

    fn fingerprint_morgan_with_output_list(
        &self,
        py: Python<'_>,
    ) -> PyResult<Vec<Option<BatchFingerprintOutput>>> {
        self.inner
            .fingerprint_morgan_with_output_list()
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| BatchFingerprintOutput { inner }))
                    .collect()
            })
            .map_err(|source| batch_error(py, source))
    }
    #[pyo3(signature=(options, collect_additional_output=true, params=None))]
    fn fingerprint_morgan_with_output_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_fingerprint_values::MorganFingerprintParams,
        collect_additional_output: bool,
        params: Option<&crate::canonical_batch_params::BatchQueryParams>,
    ) -> PyResult<Vec<Option<BatchFingerprintOutput>>> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn fingerprint_morgan_with_output_list(
        //         &self,
        //         radius: u32,
        //         n_bits: usize,
        //         include_chirality: bool,
        //         use_bond_types: bool,
        //         count_simulation: bool,
        //         count_bounds: Option<Vec<u32>>,
        //         only_nonzero_invariants: bool,
        //         include_redundant_environments: bool,
        //         from_atoms: Option<Vec<usize>>,
        //         ignore_atoms: Option<Vec<usize>>,
        //         custom_atom_invariants: Option<Vec<u32>>,
        //         custom_bond_invariants: Option<Vec<u32>>,
        //         atom_invariants_generator: Option<&str>,
        //         atom_invariants_include_ring_membership: bool,
        //         bond_invariants_generator: Option<&str>,
        //         bond_invariants_use_bond_types: bool,
        //         bond_invariants_use_chirality: bool,
        //         num_bits_per_feature: u32,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Vec<Option<MorganFingerprintResult>>> {
        //         let params = make_morgan_fingerprint_params(
        //             radius,
        //             n_bits,
        //             include_chirality,
        //             use_bond_types,
        //             count_simulation,
        //             count_bounds,
        //             only_nonzero_invariants,
        //             include_redundant_environments,
        //             from_atoms,
        //             ignore_atoms,
        //             custom_atom_invariants,
        //             custom_bond_invariants,
        //             atom_invariants_generator,
        //             atom_invariants_include_ring_membership,
        //             bond_invariants_generator,
        //             bond_invariants_use_bond_types,
        //             bond_invariants_use_chirality,
        //             num_bits_per_feature,
        //             true,
        //         )?;
        //         self.inner
        //             .morgan_fingerprint_with_output_list_with_options(
        //                 &params,
        //                 validate_n_jobs(n_jobs)?,
        //                 progress_bar,
        //             )
        //             .map(|values| {
        //                 values
        //                     .into_iter()
        //                     .map(|value| value.map(Into::into))
        //                     .collect()
        //             })
        //             .map_err(batch_validation_pyerr)
        //     }
        let defaults = crate::canonical_batch_params::BatchQueryParams::default();
        params
            .unwrap_or(&defaults)
            .execute(py, |execution| {
                self.inner.fingerprint_morgan_with_output_list_with_params(
                    &options.inner,
                    collect_additional_output,
                    execution,
                )
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| BatchFingerprintOutput { inner }))
                    .collect()
            })
    }

    /// Parse SMILES strings into a MoleculeBatch, preserving input order.
    ///
    /// Defaults: sanitize=True, errors="raise", one worker. Invalid SMILES raise
    /// BatchValidationError; errors="keep" retains failed slots and errors().
    /// Accepts SmilesParseParams as parse and BatchParams as params,
    /// or their fields as keywords. For example:
    /// ck.MoleculeBatch.from_smiles_list(["CCO", "C1CC"], errors="keep", n_jobs=1).
    #[staticmethod]
    fn from_smiles_list(py: Python<'_>, smiles: Vec<String>) -> PyResult<Self> {
        ck::MoleculeBatch::from_smiles_list(&smiles)
            .map(|inner| Self { inner })
            .map_err(|source| batch_error(py, source))
    }
    /// Parse SMILES using explicit parse and execution params, preserving input order.
    /// BatchParams.errors controls whether failures raise BatchValidationError
    /// or remain as failed records accessible through errors().
    #[staticmethod]
    fn from_smiles_list_with_params(
        py: Python<'_>,
        smiles: Vec<String>,
        parse: &crate::canonical_values::SmilesParseParams,
        params: &crate::canonical_batch_params::BatchParams,
    ) -> PyResult<Self> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn from_smiles_list(
        //         _cls: &Bound<'_, PyType>,
        //         smiles: Vec<String>,
        //         sanitize: Option<bool>,
        //         errors: Option<&Bound<'_, PyAny>>,
        //         n_jobs: Option<usize>,
        //     ) -> PyResult<Self> {
        //         let sanitize = sanitize.unwrap_or(true);
        //         let mode = parse_batch_error_mode(errors)?;
        //         let inner = cosmolkit_core::MoleculeBatch::from_smiles_list_with_sanitize_and_options(
        //             &smiles,
        //             sanitize,
        //             mode,
        //             validate_n_jobs(n_jobs)?,
        //             None,
        //         )
        //         .map_err(batch_validation_pyerr)?;
        //         Ok(Self { inner })
        //     }
        ck::MoleculeBatch::from_smiles_list_with_params(&smiles, &parse.inner, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|source| batch_error(py, source))
    }
    fn invalid_mask(&self) -> Vec<bool> {
        self.inner.invalid_mask()
    }

    fn fingerprint_morgan_list_with_generator_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_fingerprint_values::MorganParams,
        atom_invariants: Option<
            &crate::canonical_fingerprint_values::MorganAtomInvariantsGenerator,
        >,
        bond_invariants: Option<
            &crate::canonical_fingerprint_values::MorganBondInvariantsGenerator,
        >,
        call: &crate::canonical_fingerprint_values::MorganCallParams,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<FingerprintBatch> {
        params
            .execute(py, |execution| {
                self.inner.fingerprint_morgan_list_with_generator_params(
                    &options.inner,
                    atom_invariants.map(|value| &value.inner),
                    bond_invariants.map(|value| &value.inner),
                    &call.inner,
                    execution,
                )
            })
            .and_then(|values| FingerprintBatch::from_values(py, values))
    }
    fn fingerprint_morgan_with_output_list_with_generator_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_fingerprint_values::MorganParams,
        atom_invariants: Option<
            &crate::canonical_fingerprint_values::MorganAtomInvariantsGenerator,
        >,
        bond_invariants: Option<
            &crate::canonical_fingerprint_values::MorganBondInvariantsGenerator,
        >,
        call: &crate::canonical_fingerprint_values::MorganCallParams,
        collect_additional_output: bool,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<BatchFingerprintOutput>>> {
        params
            .execute(py, |execution| {
                self.inner
                    .fingerprint_morgan_with_output_list_with_generator_params(
                        &options.inner,
                        atom_invariants.map(|value| &value.inner),
                        bond_invariants.map(|value| &value.inner),
                        &call.inner,
                        collect_additional_output,
                        execution,
                    )
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| BatchFingerprintOutput { inner }))
                    .collect()
            })
    }
    fn len(&self) -> usize {
        self.inner.len()
    }
    /// Return the error policy inherited by transforms and exports.
    #[gen_stub(override_return_type(type_repr = "BatchErrorMode"))]
    fn error_mode(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let name = match self.inner.error_mode() {
            ck::BatchErrorMode::Strict => "RAISE",
            ck::BatchErrorMode::KeepErrors => "KEEP",
        };
        Ok(py
            .import("cosmolkit")?
            .getattr("BatchErrorMode")?
            .getattr(name)?
            .unbind())
    }
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    fn with_valid_records(&self) -> Self {
        Self {
            inner: self.inner.with_valid_records(),
        }
    }

    /// Return a new batch with explicit hydrogens; leave the source batch unchanged.
    ///
    /// Input order and the batch error policy are preserved. Override errors="raise"
    /// to reject failures or errors="keep" to retain existing and new failed records.
    /// Accepts AddHsParams as options and BatchParams as params,
    /// or their fields as keywords: batch.with_hydrogens(errors="keep", n_jobs=1).
    /// Omitted n_jobs uses the batch setting, or one worker if unset.
    fn with_hydrogens(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_hydrogens()
            .map(|inner| Self { inner })
            .map_err(|source| batch_error(py, source))
    }
    /// Return a new hydrogen-added batch using explicit options and execution params.
    /// Input order is preserved; the source batch is unchanged. BatchParams.errors
    /// controls whether existing and new failures raise or remain as failed records.
    fn with_hydrogens_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_chemistry_values::AddHsParams,
        params: &crate::canonical_batch_params::BatchParams,
    ) -> PyResult<Self> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn with_hydrogens(
        //         &self,
        //         errors: Option<&Bound<'_, PyAny>>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Self> {
        //         let mode = parse_batch_error_mode(errors)?;
        //         let inner = self
        //             .inner
        //             .with_hydrogens_with_options(mode, validate_n_jobs(n_jobs)?, progress_bar)
        //             .map_err(batch_validation_pyerr)?;
        //         Ok(Self { inner })
        //     }
        self.inner
            .with_hydrogens_with_params(&options.inner, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|source| batch_error(py, source))
    }

    fn without_hydrogens(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .without_hydrogens()
            .map(|inner| Self { inner })
            .map_err(|source| batch_error(py, source))
    }
    fn without_hydrogens_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_chemistry_values::RemoveHsParams,
        params: &crate::canonical_batch_params::BatchParams,
    ) -> PyResult<Self> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn without_hydrogens(
        //         &self,
        //         errors: Option<&Bound<'_, PyAny>>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Self> {
        //         let mode = parse_batch_error_mode(errors)?;
        //         let inner = self
        //             .inner
        //             .without_hydrogens_with_options(mode, validate_n_jobs(n_jobs)?, progress_bar)
        //             .map_err(batch_validation_pyerr)?;
        //         Ok(Self { inner })
        //     }
        self.inner
            .without_hydrogens_with_params(&options.inner, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|source| batch_error(py, source))
    }

    fn with_2d_coordinates(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_2d_coordinates()
            .map(|inner| Self { inner })
            .map_err(|source| batch_error(py, source))
    }
    fn with_2d_coordinates_with_params(
        &self,
        py: Python<'_>,
        options: &crate::drawing_binding::Coordinate2DParams,
        params: &crate::canonical_batch_params::BatchParams,
    ) -> PyResult<Self> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn with_2d_coordinates(
        //         &self,
        //         errors: Option<&Bound<'_, PyAny>>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Self> {
        //         let mode = parse_batch_error_mode(errors)?;
        //         let inner = self
        //             .inner
        //             .with_2d_coordinates_with_options(mode, validate_n_jobs(n_jobs)?, progress_bar)
        //             .map_err(batch_validation_pyerr)?;
        //         Ok(Self { inner })
        //     }
        self.inner
            .with_2d_coordinates_with_params(&options.inner, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|source| batch_error(py, source))
    }

    fn sanitize(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .sanitize()
            .map(|inner| Self { inner })
            .map_err(|source| batch_error(py, source))
    }
    fn sanitize_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_chemistry_values::SanitizeParams,
        params: &crate::canonical_batch_params::BatchParams,
    ) -> PyResult<Self> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855 source; typed transport only.
        //     fn sanitize(
        //         &self,
        //         strict: Option<bool>,
        //         errors: Option<&Bound<'_, PyAny>>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Self> {
        //         reject_non_strict_sanitize(strict)?;
        //         let mode = parse_batch_error_mode(errors)?;
        //         let inner = self
        //             .inner
        //             .sanitize_with_options(mode, validate_n_jobs(n_jobs)?, progress_bar)
        //             .map_err(batch_validation_pyerr)?;
        //         Ok(Self { inner })
        //     }
        self.inner
            .sanitize_with_params(&options.inner, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|source| batch_error(py, source))
    }

    fn with_kekulized_bonds(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_kekulized_bonds()
            .map(|inner| Self { inner })
            .map_err(|source| batch_error(py, source))
    }
    fn with_kekulized_bonds_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_chemistry_values::KekulizeParams,
        params: &crate::canonical_batch_params::BatchParams,
    ) -> PyResult<Self> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn with_kekulized_bonds(
        //         &self,
        //         clear_aromatic_flags: Option<bool>,
        //         errors: Option<&Bound<'_, PyAny>>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Self> {
        //         let clear_aromatic_flags = clear_aromatic_flags.unwrap_or(true);
        //         let mode = parse_batch_error_mode(errors)?;
        //         let inner = self
        //             .inner
        //             .with_kekulized_bonds_with_options(
        //                 clear_aromatic_flags,
        //                 mode,
        //                 validate_n_jobs(n_jobs)?,
        //                 progress_bar,
        //             )
        //             .map_err(batch_validation_pyerr)?;
        //         Ok(Self { inner })
        //     }
        self.inner
            .with_kekulized_bonds_with_params(&options.inner, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|source| batch_error(py, source))
    }

    /// Return SMILES strings in input order; defaults to canonical, isomeric SMILES.
    ///
    /// Existing failed records produce None; new writing errors raise
    /// BatchValidationError. The source batch is unchanged.
    /// Accepts SmilesWriteParams as options and BatchQueryParams as params,
    /// or their fields as keywords: batch.to_smiles_list(canonical=False, n_jobs=1).
    /// Omitted n_jobs uses the batch setting, or one worker if unset.
    fn to_smiles_list(&self, py: Python<'_>) -> PyResult<Vec<Option<String>>> {
        self.inner
            .to_smiles_list()
            .map_err(|source| batch_error(py, source))?
            .iter()
            .map(|value| {
                value
                    .as_ref()
                    .map(|text| crate::canonical_sdf::decode_source_text(py, text))
                    .transpose()
            })
            .collect()
    }
    /// Return ordered SMILES strings using explicit options and execution params.
    /// Existing failed records produce None; new writing errors raise
    /// BatchValidationError. Does not modify the source batch.
    fn to_smiles_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::canonical_values::SmilesWriteParams,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<String>>> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn to_smiles_list(
        //         &self,
        //         isomeric_smiles: bool,
        //         canonical: bool,
        //         kekule: bool,
        //         clean_stereo: bool,
        //         all_bonds_explicit: bool,
        //         all_hs_explicit: bool,
        //         include_dative_bonds: bool,
        //         ignore_atom_map_numbers: bool,
        //         rooted_at_atom: Option<usize>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Vec<Option<String>>> {
        //         let params = make_smiles_write_params(
        //             isomeric_smiles,
        //             canonical,
        //             kekule,
        //             clean_stereo,
        //             all_bonds_explicit,
        //             all_hs_explicit,
        //             include_dative_bonds,
        //             ignore_atom_map_numbers,
        //             rooted_at_atom,
        //         );
        //         self.inner
        //             .to_smiles_optional_list_with_params_and_options(
        //                 &params,
        //                 validate_n_jobs(n_jobs)?,
        //                 progress_bar,
        //             )
        //             .map_err(batch_validation_pyerr)
        //     }
        params
            .execute(py, |execution| {
                self.inner
                    .to_smiles_list_with_params(&options.inner, execution)
            })?
            .iter()
            .map(|value| {
                value
                    .as_ref()
                    .map(|text| crate::canonical_sdf::decode_source_text(py, text))
                    .transpose()
            })
            .collect()
    }

    #[gen_stub(override_return_type(type_repr="list[numpy.ndarray | None]",imports=("numpy")))]
    fn dg_bounds_matrix_list<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, pyo3::types::PyList>> {
        bounds_values(
            py,
            self.inner
                .dg_bounds_matrix_list()
                .map_err(|source| batch_error(py, source))?,
        )
    }
    #[gen_stub(override_return_type(type_repr="list[numpy.ndarray | None]",imports=("numpy")))]
    fn dg_bounds_matrix_list_with_params<'py>(
        &self,
        py: Python<'py>,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<Bound<'py, pyo3::types::PyList>> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn dg_bounds_matrix_list<'py>(
        //         &self,
        //         py: Python<'py>,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Bound<'py, PyList>> {
        //         let values = self
        //             .inner
        //             .dg_bounds_matrix_list_with_options(validate_n_jobs(n_jobs)?, progress_bar)
        //             .map_err(batch_validation_pyerr)?;
        //         let out = PyList::empty(py);
        //         for value in values {
        //             if let Some(matrix) = value {
        //                 out.append(PyArray2::from_vec2(py, &matrix).map_err(|err| {
        //                     PyValueError::new_err(format!(
        //                         "MoleculeBatch.dg_bounds_matrix_list failed: {err}"
        //                     ))
        //                 })?)?;
        //             } else {
        //                 out.append(py.None())?;
        //             }
        //         }
        //         Ok(out)
        //     }
        bounds_values(
            py,
            params.execute(py, |execution| {
                self.inner.dg_bounds_matrix_list_with_params(execution)
            })?,
        )
    }

    fn to_svg_list(
        &self,
        py: Python<'_>,
        width: u32,
        height: u32,
    ) -> PyResult<Vec<Option<String>>> {
        self.inner
            .to_svg_list(width, height)
            .map_err(|source| batch_error(py, source))
    }
    fn to_svg_list_with_params(
        &self,
        py: Python<'_>,
        width: u32,
        height: u32,
        params: &crate::canonical_batch_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<String>>> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn to_svg_list(
        //         &self,
        //         width: u32,
        //         height: u32,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<Vec<Option<String>>> {
        //         self.inner
        //             .to_svg_list_with_options(width, height, validate_n_jobs(n_jobs)?, progress_bar)
        //             .map_err(batch_validation_pyerr)
        //     }
        params.execute(py, |execution| {
            self.inner.to_svg_list_with_params(width, height, execution)
        })
    }

    fn write_images(&self, py: Python<'_>, directory: &str) -> PyResult<BatchExportReport> {
        crate::user_path::with_image_user_paths(
            directory,
            None,
            || std::env::var_os("HOME"),
            |directory| self.inner.write_images(directory),
            |path, report| report.write_report(path),
        )
        .map(|inner| BatchExportReport { inner })
        .map_err(|source| image_path_error(py, source))
    }
    fn write_images_with_params(
        &self,
        py: Python<'_>,
        directory: &str,
        options: &crate::canonical_batch_params::BatchImageParams,
    ) -> PyResult<BatchExportReport> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855;
        // original values/defaults/order transported through canonical typed parameters.
        //     fn to_images(
        //         &self,
        //         out_dir: &str,
        //         format: Option<&str>,
        //         size: Option<(u32, u32)>,
        //         n_jobs: Option<usize>,
        //         errors: Option<&Bound<'_, PyAny>>,
        //         report_path: Option<&str>,
        //         filenames: Option<Vec<Option<String>>>,
        //         progress_bar: Option<bool>,
        //     ) -> PyResult<PyBatchExportReport> {
        //         let mode = parse_batch_error_mode(errors)?;
        //         let image_format = format.unwrap_or("png").to_string();
        //         let (width, height) = size.unwrap_or((300, 300));
        //         let out_dir = expand_user_path(out_dir)?;
        //         let filenames = complete_batch_filenames(filenames, self.inner.len(), &image_format)?;
        //         let report = self
        //             .inner
        //             .write_images_with_options(
        //                 out_dir.as_path(),
        //                 &image_format,
        //                 width,
        //                 height,
        //                 mode,
        //                 filenames.as_deref(),
        //                 validate_n_jobs(n_jobs)?,
        //                 progress_bar,
        //             )
        //             .map_err(batch_validation_pyerr)?;
        //         if let Some(path) = report_path {
        //             write_batch_report(path, &report)?;
        //         }
        //         Ok(report.into())
        //     }
        // BatchImageParams::new accepts Option<String> and stores it directly
        // as PathBuf: this conversion is lossless for every Python input.
        // Expanded non-Unicode HOME bytes stay in PathBuf below and are never
        // converted back to text. Borrowed frozen options remain unchanged.
        let report_path = options
            .inner
            .report_path
            .as_ref()
            .map(|path| path.to_string_lossy());
        crate::user_path::with_image_user_paths(
            directory,
            report_path.as_deref(),
            || std::env::var_os("HOME"),
            |directory| {
                let options = ck::BatchImageParams {
                    report_path: None,
                    ..options.inner.clone()
                };
                self.inner.write_images_with_params(directory, &options)
            },
            |path, report| report.write_report(path),
        )
        .map(|inner| BatchExportReport { inner })
        .map_err(|source| image_path_error(py, source))
    }

    #[classmethod]
    #[pyo3(signature=(text,errors=None,n_jobs=None))]
    fn from_sdf_records(
        _cls: &Bound<'_, PyType>,
        py: Python<'_>,
        text: &str,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
    ) -> PyResult<Self> {
        ck::MoleculeBatch::from_sdf_records_with_params(
            text,
            &ck::SdfReadParams::default(),
            error_mode(errors)?,
            self::n_jobs(n_jobs)?,
        )
        .map(|inner| Self { inner })
        .map_err(|e| batch_error(py, e))
    }
    #[classmethod]
    #[pyo3(signature=(path,errors=None,n_jobs=None,progress_bar=false))]
    fn read_sdf(
        _cls: &Bound<'_, PyType>,
        py: Python<'_>,
        path: TextPath,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: bool,
    ) -> PyResult<Self> {
        ck::MoleculeBatch::read_sdf_with_params(
            path.as_str(),
            &ck::SdfReadParams::default(),
            error_mode(errors)?,
            self::n_jobs(n_jobs)?,
            progress_bar,
        )
        .map(|inner| Self { inner })
        .map_err(|e| batch_error(py, e))
    }
    #[pyo3(signature=(path,format=None,errors=None,n_jobs=None,report_path=None,progress_bar=None))]
    fn write_sdf(
        &self,
        py: Python<'_>,
        path: TextPath,
        format: Option<&str>,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        report_path: Option<TextPath>,
        progress_bar: Option<bool>,
    ) -> PyResult<BatchExportReport> {
        let mode = errors.map(|value| error_mode(Some(value))).transpose()?;
        let format = write_params(format, true, true)?.format;
        let params = ck::BatchExportParams {
            format,
            errors: mode,
            n_jobs: self::n_jobs(n_jobs)?,
            progress_bar,
        };
        self.inner
            .write_sdf_with_params(
                path.as_str(),
                &params,
                report_path.as_ref().map(TextPath::as_str),
            )
            .map(|inner| BatchExportReport { inner })
            .map_err(|e| batch_error(py, e))
    }
    #[pyo3(signature=(out_dir,format=None,errors=None,n_jobs=None,report_path=None,filenames=None,progress_bar=None))]
    fn write_sdf_files(
        &self,
        py: Python<'_>,
        out_dir: TextPath,
        format: Option<&str>,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        report_path: Option<TextPath>,
        filenames: Option<Vec<Option<String>>>,
        progress_bar: Option<bool>,
    ) -> PyResult<BatchExportReport> {
        let mode = errors.map(|value| error_mode(Some(value))).transpose()?;
        let format = write_params(format, true, true)?.format;
        let params = ck::BatchExportParams {
            format,
            errors: mode,
            n_jobs: self::n_jobs(n_jobs)?,
            progress_bar,
        };
        self.inner
            .write_sdf_files_with_params(
                out_dir.as_str(),
                &params,
                filenames.as_deref(),
                report_path.as_ref().map(TextPath::as_str),
            )
            .map(|inner| BatchExportReport { inner })
            .map_err(|e| batch_error(py, e))
    }
    fn to_list(&self) -> Vec<Option<Molecule>> {
        self.inner
            .records()
            .iter()
            .map(|r| match r {
                ck::BatchRecord::Molecule(m) => Some(Molecule::from_inner(m.clone())),
                ck::BatchRecord::Error(_) => None,
            })
            .collect()
    }
    #[gen_stub(override_return_type(type_repr="typing.Iterator[Molecule | None]",imports=("typing")))]
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let list = pyo3::types::PyList::new(py, self.to_list())?;
        Ok(pyo3::types::PyIterator::from_object(list.as_any())?.into_any())
    }
    fn __len__(&self) -> usize {
        self.inner.len()
    }
    fn valid_mask(&self) -> Vec<bool> {
        self.inner.valid_mask()
    }
    fn valid_count(&self) -> usize {
        self.inner.valid_count()
    }
    fn invalid_count(&self) -> usize {
        self.inner.invalid_count()
    }
    fn errors(&self) -> Vec<BatchError> {
        self.inner
            .errors()
            .into_iter()
            .map(|inner| BatchError { inner })
            .collect()
    }
    fn parallel_jobs(&self) -> Option<usize> {
        self.inner.parallel_jobs()
    }
    fn progress_bar(&self) -> Option<bool> {
        self.inner.progress_bar()
    }
    #[pyo3(signature=(n_jobs))]
    fn with_parallel_jobs(&self, py: Python<'_>, n_jobs: Option<usize>) -> PyResult<Self> {
        self.inner
            .clone()
            .with_parallel_jobs(self::n_jobs(n_jobs)?)
            .map(|inner| Self { inner })
            .map_err(|e| batch_error(py, e))
    }
    #[pyo3(signature=(progress_bar))]
    fn with_progress_bar(&self, progress_bar: Option<bool>) -> Self {
        Self {
            inner: self.inner.clone().with_progress_bar(progress_bar),
        }
    }
    #[gen_stub(skip)]
    fn __getitem__(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        // COSMolKit❗✔️: pinned d892ec3 python/src/lib.rs::__getitem__; canonical projection preserves source branches/defaults.
        // fn __getitem__(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        //         if let Ok(slice) = key.cast::<PySlice>() {
        //             let indices = self.slice_indices(slice)?;
        //             return self.selected_batch_pyobject(py, &indices);
        //         }
        //         if !key.is_exact_instance_of::<PyBool>() {
        //             if let Ok(index) = key.extract::<isize>() {
        //                 return self.get_record_pyobject(py, self.normalize_index(index)?);
        //             }
        //         }
        //         match self.sequence_indices(key) {
        //             Ok(indices) => return self.selected_batch_pyobject(py, &indices),
        //             Err(error) => {
        //                 if key.extract::<Vec<Py<PyAny>>>().is_ok() {
        //                     return Err(error);
        //                 }
        //             }
        //         }
        //         Err(PyTypeError::new_err(
        //             "MoleculeBatch indices must be integers, slices, integer lists, or boolean masks",
        //         ))
        //     }

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

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct SdfDataset {
    inner: ck::SdfDataset,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl SdfDataset {
    #[staticmethod]
    fn open_with_params(
        py: Python<'_>,
        path: TextPath,
        params: &crate::canonical_sdf::SdfReadParams,
    ) -> PyResult<Self> {
        ck::SdfDataset::open_with_params(path.as_str(), &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    fn len(&self) -> usize {
        self.inner.len()
    }
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    fn record(&self, py: Python<'_>, index: usize) -> PyResult<SdfRecord> {
        self.inner
            .record(index)
            .map(|inner| SdfRecord { inner })
            .map_err(|e| crate::canonical_sdf::sdf_pyerr(py, e))
    }
    fn record_with_params(
        &self,
        py: Python<'_>,
        index: usize,
        params: &crate::canonical_sdf::SdfReadParams,
    ) -> PyResult<SdfRecord> {
        self.inner
            .record_with_params(index, &params.inner)
            .map(|inner| SdfRecord { inner })
            .map_err(|e| crate::canonical_sdf::sdf_pyerr(py, e))
    }
    fn record_text(&self, py: Python<'_>, index: usize) -> PyResult<String> {
        self.inner
            .record_text(index)
            .map_err(|e| crate::canonical_sdf::sdf_pyerr(py, e))
    }
    fn iter(&self) -> SdfDatasetIterator {
        SdfDatasetIterator {
            inner: self.inner.iter(),
        }
    }
    #[classmethod]
    #[pyo3(signature=(path,index=None,build=None))]
    fn open(
        _cls: &Bound<'_, PyType>,
        py: Python<'_>,
        path: TextPath,
        index: Option<&Bound<'_, PyAny>>,
        build: Option<&str>,
    ) -> PyResult<Self> {
        if let Some(build) = build {
            if !matches!(build, "auto" | "always" | "never") {
                return Err(PyValueError::new_err(
                    "unsupported build mode, expected one of: auto, always, never",
                ));
            }
        }
        if let Some(index) = index {
            if let Ok(value) = index.extract::<String>() {
                if !matches!(value.as_str(), "auto" | "memory") {
                    return Err(PyNotImplementedError::new_err(
                        "persistent SDF sidecar indexes are not implemented yet; use index='auto', index='memory', or None",
                    ));
                }
            }
        }
        ck::SdfDataset::open_with_params(path.as_str(), &ck::SdfReadParams::default())
            .map(|inner| Self { inner })
            .map_err(|e| io_error(py, e))
    }
    fn path(&self) -> String {
        self.inner.path().to_string_lossy().into_owned()
    }
    fn __len__(&self) -> usize {
        self.inner.len()
    }
    fn metadata(&self, index: isize) -> PyResult<SdfRecordMetadata> {
        let len = self.inner.len() as isize;
        let index = if index < 0 { len + index } else { index };
        if index < 0 || index >= len {
            return Err(PyIndexError::new_err("SdfDataset index out of range"));
        }
        Ok(SdfRecordMetadata {
            inner: self
                .inner
                .metadata(index as usize)
                .expect("checked metadata index")
                .clone(),
        })
    }
    #[gen_stub(skip)]
    fn __getitem__(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match sdf_indices_from_key(self.inner.len(), key)? {
            Err(index) => self
                .inner
                .record(index)
                .map_err(|e| sdf_error(py, e))
                .and_then(|inner| Ok(Py::new(py, SdfRecord { inner })?.into_any())),
            Ok(indices) => ck::MoleculeBatch::from_dataset_indices(
                &self.inner,
                &indices,
                ck::BatchErrorMode::Strict,
            )
            .map_err(|e| batch_error(py, e))
            .and_then(|inner| Ok(Py::new(py, MoleculeBatch { inner })?.into_any())),
        }
    }
    fn __iter__(&self) -> SdfDatasetIterator {
        SdfDatasetIterator {
            inner: self.inner.iter(),
        }
    }
    #[pyo3(signature=(size=1024,indices=None,errors=None,n_jobs=None,progress_bar=false))]
    fn batches(
        &self,
        py: Python<'_>,
        size: usize,
        indices: Option<Vec<usize>>,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: bool,
    ) -> PyResult<SdfBatchIterator> {
        if size == 0 {
            return Err(PyValueError::new_err("size must be >= 1"));
        }
        if indices.is_some() {
            return Err(PyNotImplementedError::new_err(
                "indices=... for SdfDataset.batches() is not implemented yet",
            ));
        }
        self.inner
            .batches(
                size,
                error_mode(errors)?,
                self::n_jobs(n_jobs)?,
                progress_bar,
            )
            .map(|inner| SdfBatchIterator { inner })
            .map_err(|e| batch_error(py, e))
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct SdfDatasetIterator {
    inner: ck::SdfDatasetIterator,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl SdfDatasetIterator {
    fn __length_hint__(&self) -> usize {
        self.inner.len()
    }
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    #[gen_stub(override_return_type(type_repr = "SdfRecord"))]
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<SdfRecord>> {
        self.inner
            .next()
            .transpose()
            .map(|v| v.map(|inner| SdfRecord { inner }))
            .map_err(|e| sdf_error(py, e))
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct SdfBatchIterator {
    inner: ck::SdfBatchIterator,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl SdfBatchIterator {
    fn next_batch(&mut self, py: Python<'_>) -> PyResult<Option<MoleculeBatch>> {
        self.inner
            .next_batch()
            .map(|v| v.map(|inner| MoleculeBatch { inner }))
            .map_err(|e| batch_error(py, e))
    }
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    #[gen_stub(override_return_type(type_repr = "MoleculeBatch"))]
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<MoleculeBatch>> {
        self.inner
            .next_batch()
            .map(|v| v.map(|inner| MoleculeBatch { inner }))
            .map_err(|e| batch_error(py, e))
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct SdfReader {
    inner: ck::SdfReader,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SdfReader {
    #[staticmethod]
    fn open_with_params(
        py: Python<'_>,
        path: TextPath,
        params: &crate::canonical_sdf::SdfReadParams,
    ) -> PyResult<Self> {
        ck::SdfReader::open_with_params(path.as_str(), &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    fn path(&self) -> std::path::PathBuf {
        self.inner.path().to_path_buf()
    }
    fn params(&self) -> crate::canonical_sdf::SdfReadParams {
        crate::canonical_sdf::SdfReadParams {
            inner: *self.inner.params(),
        }
    }
    #[classmethod]
    #[pyo3(signature=(path))]
    fn open(_cls: &Bound<'_, PyType>, py: Python<'_>, path: TextPath) -> PyResult<Self> {
        // Python's source API stores a path/parameter configuration. Opening the
        // actual file is deferred to batches(), which delegates to Rust.
        ck::SdfReader::open_with_params(path.as_str(), &ck::SdfReadParams::default())
            .map(|inner| Self { inner })
            .map_err(|e| io_error(py, e))
    }
    #[pyo3(signature=(size=1024,errors=None,n_jobs=None,progress_bar=false))]
    fn batches(
        &self,
        py: Python<'_>,
        size: usize,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: bool,
    ) -> PyResult<SdfReaderBatchIterator> {
        if size == 0 {
            return Err(PyValueError::new_err("size must be >= 1"));
        }
        if progress_bar {
            return Err(PyNotImplementedError::new_err(
                "SdfReader.batches(progress_bar=True) cannot show an accurate total for forward-only streams; use SdfDataset.batches() for indexed progress",
            ));
        }
        self.inner
            .batches(size, error_mode(errors)?, self::n_jobs(n_jobs)?)
            .map(|inner| SdfReaderBatchIterator { inner })
            .map_err(|e| batch_error(py, e))
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct SdfReaderBatchIterator {
    pub(crate) inner: ck::SdfReaderBatchIterator<std::io::BufReader<std::fs::File>>,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl SdfReaderBatchIterator {
    fn next_batch(&mut self, py: Python<'_>) -> PyResult<Option<MoleculeBatch>> {
        self.inner
            .next_batch()
            .map(|v| v.map(|inner| MoleculeBatch { inner }))
            .map_err(|e| batch_error(py, e))
    }
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    #[gen_stub(override_return_type(type_repr = "MoleculeBatch"))]
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<MoleculeBatch>> {
        self.inner
            .next_batch()
            .map(|v| v.map(|inner| MoleculeBatch { inner }))
            .map_err(|e| batch_error(py, e))
    }
}

pub(crate) fn write_params(
    format: Option<&str>,
    include_stereo: bool,
    kekulize: bool,
) -> PyResult<ck::MolBlockWriteParams> {
    let format = match format.map(str::to_ascii_lowercase).as_deref() {
        None | Some("v2000" | "v2k") => ck::SdfFormat::V2000,
        Some("v3000" | "v3k") => ck::SdfFormat::V3000,
        Some(value) => {
            return Err(PyValueError::new_err(format!(
                "unsupported SDF format '{value}', expected one of: v2000, v3000"
            )));
        }
    };
    Ok(ck::MolBlockWriteParams {
        format,
        include_stereo,
        kekulize,
        ..Default::default()
    })
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct BatchExportReport {
    inner: ck::BatchExportReport,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BatchExportReport {
    fn write_report(&self, py: Python<'_>, path: TextPath) -> PyResult<()> {
        self.inner
            .write_report(std::path::Path::new(path.as_str()))
            .map_err(|e| batch_error(py, e))
    }
    fn total(&self) -> usize {
        self.inner.total()
    }
    fn success(&self) -> usize {
        self.inner.success()
    }
    fn failed(&self) -> usize {
        self.inner.failed()
    }
    fn errors(&self) -> Vec<BatchError> {
        self.inner
            .errors()
            .iter()
            .cloned()
            .map(|inner| BatchError { inner })
            .collect()
    }
    fn __repr__(&self) -> String {
        format!(
            "BatchExportReport(written={}, failed={})",
            self.inner.written, self.inner.failed
        )
    }
}

impl MoleculeBatch {
    fn normalize_index(&self, index: isize) -> PyResult<usize> {
        // COSMolKit❗✔️: pinned d892ec3 python/src/lib.rs::normalize_index; canonical projection preserves source branches/defaults.
        // fn normalize_index(&self, index: isize) -> PyResult<usize> {
        //         let len = self.inner.len() as isize;
        //         let index = if index < 0 { len + index } else { index };
        //         if index < 0 || index >= len {
        //             return Err(PyIndexError::new_err("MoleculeBatch index out of range"));
        //         }
        //         Ok(index as usize)
        //     }

        let n = self.inner.len() as isize;
        let index = if index < 0 { n + index } else { index };
        if index < 0 || index >= n {
            return Err(PyIndexError::new_err("MoleculeBatch index out of range"));
        }
        Ok(index as usize)
    }
    fn select_records(&self, indices: &[usize]) -> PyResult<Self> {
        // COSMolKit❗✔️: pinned d892ec3 python/src/lib.rs::select_records; canonical projection preserves source branches/defaults.
        // fn select_records(&self, indices: &[usize]) -> Self {
        //         let records = indices
        //             .iter()
        //             .filter_map(|index| self.inner.get(*index).cloned())
        //             .collect();
        //         let inner = cosmolkit_core::MoleculeBatch::new(records)
        //             .with_parallel_jobs(self.inner.parallel_jobs())
        //             .with_progress_bar(self.inner.progress_bar());
        //         Self { inner }
        //     }

        let records = indices
            .iter()
            .map(|index| self.inner.get(*index).expect("normalized index").clone())
            .collect();
        ck::MoleculeBatch::from_records(records, self.inner.error_mode())
            .and_then(|batch| batch.with_parallel_jobs(self.inner.parallel_jobs()))
            .map(|inner| Self {
                inner: inner.with_progress_bar(self.inner.progress_bar()),
            })
            .map_err(|e| Python::attach(|py| batch_error(py, e)))
    }
    fn selected_batch_pyobject(&self, py: Python<'_>, indices: &[usize]) -> PyResult<Py<PyAny>> {
        Ok(Py::new(py, self.select_records(indices)?)?.into_any())
    }
    fn get_record_pyobject(&self, py: Python<'_>, index: usize) -> PyResult<Py<PyAny>> {
        match self.inner.get(index) {
            Some(ck::BatchRecord::Molecule(m)) => {
                Ok(Py::new(py, Molecule::from_inner(m.clone()))?.into_any())
            }
            _ => Ok(py.None()),
        }
    }
    fn slice_indices(&self, slice: &Bound<'_, PySlice>) -> PyResult<Vec<usize>> {
        // COSMolKit❗✔️: pinned d892ec3 python/src/lib.rs::slice_indices; canonical projection preserves source branches/defaults.
        // fn slice_indices(&self, slice: &Bound<'_, PySlice>) -> PyResult<Vec<usize>> {
        //         let indices = slice.indices(self.inner.len() as isize)?;
        //         let mut out = Vec::with_capacity(indices.slicelength);
        //         let mut index = indices.start;
        //         for _ in 0..indices.slicelength {
        //             out.push(index as usize);
        //             index += indices.step;
        //         }
        //         Ok(out)
        //     }

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
        // COSMolKit❗✔️: pinned d892ec3 python/src/lib.rs::sequence_indices; canonical projection preserves source branches/defaults.
        // fn sequence_indices(&self, key: &Bound<'_, PyAny>) -> PyResult<Vec<usize>> {
        //         let items = key.extract::<Vec<Py<PyAny>>>()?;
        //         if items.is_empty() {
        //             return Ok(Vec::new());
        //         }
        //
        //         let py = key.py();
        //         let bool_mask = items
        //             .iter()
        //             .all(|item| item.bind(py).is_exact_instance_of::<PyBool>());
        //         if bool_mask {
        //             if items.len() != self.inner.len() {
        //                 return Err(PyIndexError::new_err(format!(
        //                     "boolean mask length {} does not match MoleculeBatch length {}",
        //                     items.len(),
        //                     self.inner.len()
        //                 )));
        //             }
        //             let mut indices = Vec::new();
        //             for (index, item) in items.iter().enumerate() {
        //                 if item.bind(py).extract::<bool>()? {
        //                     indices.push(index);
        //                 }
        //             }
        //             return Ok(indices);
        //         }
        //
        //         let mut indices = Vec::with_capacity(items.len());
        //         for item in items {
        //             let item = item.bind(py);
        //             if item.is_exact_instance_of::<PyBool>() {
        //                 return Err(PyTypeError::new_err(
        //                     "MoleculeBatch index lists must not mix bool and int values",
        //                 ));
        //             }
        //             indices.push(self.normalize_index(item.extract::<isize>()?)?);
        //         }
        //         Ok(indices)
        //     }

        let items = key.extract::<Vec<Py<PyAny>>>()?;
        if items.is_empty() {
            return Ok(Vec::new());
        }
        let py = key.py();
        if items
            .iter()
            .all(|item| item.bind(py).is_exact_instance_of::<PyBool>())
        {
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
fn bounds_values<'py>(
    py: Python<'py>,
    values: Vec<Option<Vec<Vec<f64>>>>,
) -> PyResult<Bound<'py, pyo3::types::PyList>> {
    let out = pyo3::types::PyList::empty(py);
    for value in values {
        match value {
            Some(matrix) => {
                out.append(numpy::PyArray2::from_vec2(py, &matrix).map_err(|source| {
                    PyValueError::new_err(format!(
                        "MoleculeBatch.dg_bounds_matrix_list failed: {source}"
                    ))
                })?)?
            }
            None => out.append(py.None())?,
        }
    }
    Ok(out)
}
fn batch_record_object(py: Python<'_>, record: &ck::BatchRecord) -> PyResult<Py<PyAny>> {
    match record {
        ck::BatchRecord::Molecule(molecule) => {
            Ok(Py::new(py, Molecule::from_inner(molecule.clone()))?.into_any())
        }
        ck::BatchRecord::Error(error) => Ok(Py::new(
            py,
            BatchError {
                inner: error.clone(),
            },
        )?
        .into_any()),
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("BatchImageError", module.py().get_type::<BatchImageError>())?;
    module.add_class::<SdfRecordMetadata>()?;
    module.add_class::<MoleculeBatch>()?;
    module.add_class::<BatchError>()?;
    module.add_class::<BatchExportReport>()?;
    module.add_class::<SdfDataset>()?;
    module.add_class::<SdfDatasetIterator>()?;
    module.add_class::<SdfBatchIterator>()?;
    module.add_class::<SdfReader>()?;
    module.add_class::<SdfReaderBatchIterator>()?;
    add_batch_validation_error_class(module)?;
    let modes = module
        .py()
        .import("enum")?
        .getattr("IntEnum")?
        .call1(("BatchErrorMode", vec![("RAISE", 1), ("KEEP", 2)]))?;
    modes.setattr("__module__", "cosmolkit")?;
    let aliases = pyo3::types::PyDict::new(module.py());
    aliases.set_item("raise", modes.getattr("RAISE")?)?;
    aliases.set_item("keep", modes.getattr("KEEP")?)?;
    let readonly = module
        .py()
        .import("types")?
        .getattr("MappingProxyType")?
        .call1((aliases,))?;
    module.add("BATCH_ERROR_MODE_MAP", readonly)?;
    module.add("BatchErrorMode", modes)?;
    Ok(())
}
