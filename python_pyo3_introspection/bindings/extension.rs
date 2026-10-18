//! SMILES configuration, module calls and the native molecule batch owner.
use crate::arguments::Argument;
use crate::io_batch_params::BatchErrorMode;
use crate::native::Molecule;
use crate::parameters::parameters;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

parameters!(SmilesParseParams, ck::SmilesParseParams, {
    sanitize: bool,
    allow_cxsmiles: bool,
    strict_cxsmiles: bool,
    parse_name: bool,
    remove_hs: bool,
    skip_cleanup: bool,
    debug_parse: bool,
}, string_maps = [replacements]);

#[pyfunction]
#[pyo3(signature = (text, params=Argument::Missing, *, sanitize=Argument::Missing, allow_cxsmiles=Argument::Missing, strict_cxsmiles=Argument::Missing, parse_name=Argument::Missing, remove_hs=Argument::Missing, skip_cleanup=Argument::Missing, debug_parse=Argument::Missing, replacements=Argument::Missing))]
pub fn mol_from_smiles(
    text: &str,
    params: Argument<Option<PyRef<'_, SmilesParseParams>>>,
    sanitize: Argument<Option<bool>>,
    allow_cxsmiles: Argument<Option<bool>>,
    strict_cxsmiles: Argument<Option<bool>>,
    parse_name: Argument<Option<bool>>,
    remove_hs: Argument<Option<bool>>,
    skip_cleanup: Argument<Option<bool>>,
    debug_parse: Argument<Option<bool>>,
    replacements: Argument<Option<BTreeMap<String, String>>>,
) -> PyResult<Molecule> {
    crate::arguments::configuration_arguments!(params; sanitize, allow_cxsmiles, strict_cxsmiles, parse_name, remove_hs, skip_cleanup, debug_parse, replacements);
    let params = params.as_deref();
    // Only these explicitly bound keywords are supported. The generator does
    // not infer call options from the fields or constructor of SmilesParseParams.
    if let Some(params) = params {
        _mol_from_smiles_params(text, Some(params))
    } else {
        _mol_from_smiles_keywords(
            text,
            sanitize,
            allow_cxsmiles,
            strict_cxsmiles,
            parse_name,
            remove_hs,
            skip_cleanup,
            debug_parse,
            replacements,
        )
    }
}

#[pyfunction]
#[doc = r#"{"role": "overload", "target": "mol_from_smiles"}"#]
#[pyo3(signature = (text, params=None))]
pub fn _mol_from_smiles_params(
    text: &str,
    params: Option<&SmilesParseParams>,
) -> PyResult<Molecule> {
    let defaults;
    let effective = if let Some(params) = params {
        &params.inner
    } else {
        defaults = ck::SmilesParseParams::default();
        &defaults
    };
    ck::Molecule::from_smiles_with_params(text, effective)
        .map(|inner| Molecule { inner })
        .map_err(crate::errors::SmilesError::from_rust)
}

#[pyfunction]
#[doc = r#"{"role": "overload", "target": "mol_from_smiles"}"#]
#[pyo3(signature = (text, *, sanitize=None, allow_cxsmiles=None, strict_cxsmiles=None, parse_name=None, remove_hs=None, skip_cleanup=None, debug_parse=None, replacements=None))]
pub fn _mol_from_smiles_keywords(
    text: &str,
    sanitize: Option<bool>,
    allow_cxsmiles: Option<bool>,
    strict_cxsmiles: Option<bool>,
    parse_name: Option<bool>,
    remove_hs: Option<bool>,
    skip_cleanup: Option<bool>,
    debug_parse: Option<bool>,
    replacements: Option<BTreeMap<String, String>>,
) -> PyResult<Molecule> {
    let mut value = ck::SmilesParseParams::default();
    if let Some(option) = sanitize {
        value.sanitize = option;
    }
    if let Some(option) = allow_cxsmiles {
        value.allow_cxsmiles = option;
    }
    if let Some(option) = strict_cxsmiles {
        value.strict_cxsmiles = option;
    }
    if let Some(option) = parse_name {
        value.parse_name = option;
    }
    if let Some(option) = remove_hs {
        value.remove_hs = option;
    }
    if let Some(option) = skip_cleanup {
        value.skip_cleanup = option;
    }
    if let Some(option) = debug_parse {
        value.debug_parse = option;
    }
    if let Some(option) = replacements {
        value.replacements = option;
    }
    _mol_from_smiles_params(text, Some(&SmilesParseParams { inner: value }))
}

#[pyclass(module = "fixed_projection._native", frozen)]
pub struct MoleculeBatch {
    pub(crate) inner: ck::MoleculeBatch,
}

#[pymethods]
impl MoleculeBatch {
    #[pyo3(signature=(path,params=None,*,format=None,errors=None,n_jobs=None,progress_bar=None,report_path=None))]
    fn write_sdf(
        &self,
        py: Python<'_>,
        path: std::path::PathBuf,
        params: Option<&crate::io_batch_params::BatchExportParams>,
        format: Option<crate::io_batch_params::SdfFormatInput>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
        report_path: Option<std::path::PathBuf>,
    ) -> PyResult<crate::batch_errors::BatchExportReport> {
        if params.is_some()
            && (format.is_some() || errors.is_some() || n_jobs.is_some() || progress_bar.is_some())
        {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Params and export keywords are mutually exclusive",
            ));
        }
        let defaults = ck::BatchExportParams::default();
        let resolved = ck::BatchExportParams {
            format: format
                .map(crate::io_batch_params::SdfFormatInput::into_rust)
                .transpose()?
                .unwrap_or(defaults.format),
            errors: errors
                .map(crate::io_batch_params::BatchErrorModeInput::into_rust)
                .transpose()?,
            n_jobs,
            progress_bar,
        };
        let path = crate::molecular_io_calls::path_text(&path)?;
        let report_path = report_path
            .as_deref()
            .map(crate::molecular_io_calls::path_text)
            .transpose()?;
        py.detach(|| {
            self.inner.write_sdf_with_params(
                path,
                params.map(|v| &v.inner).unwrap_or(&resolved),
                report_path,
            )
        })
        .map(|inner| crate::batch_errors::BatchExportReport { inner })
        .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[pyo3(signature=(path,params,report_path))]
    fn write_sdf_with_params(
        &self,
        py: Python<'_>,
        path: std::path::PathBuf,
        params: &crate::io_batch_params::BatchExportParams,
        report_path: Option<std::path::PathBuf>,
    ) -> PyResult<crate::batch_errors::BatchExportReport> {
        self.write_sdf(py, path, Some(params), None, None, None, None, report_path)
    }
    #[pyo3(signature=(directory,params=None,*,format=None,errors=None,n_jobs=None,progress_bar=None,filenames=None,report_path=None))]
    fn write_sdf_files(
        &self,
        py: Python<'_>,
        directory: std::path::PathBuf,
        params: Option<&crate::io_batch_params::BatchExportParams>,
        format: Option<crate::io_batch_params::SdfFormatInput>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
        filenames: Option<Vec<Option<String>>>,
        report_path: Option<std::path::PathBuf>,
    ) -> PyResult<crate::batch_errors::BatchExportReport> {
        if params.is_some()
            && (format.is_some() || errors.is_some() || n_jobs.is_some() || progress_bar.is_some())
        {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Params and export keywords are mutually exclusive",
            ));
        }
        let defaults = ck::BatchExportParams::default();
        let resolved = ck::BatchExportParams {
            format: format
                .map(crate::io_batch_params::SdfFormatInput::into_rust)
                .transpose()?
                .unwrap_or(defaults.format),
            errors: errors
                .map(crate::io_batch_params::BatchErrorModeInput::into_rust)
                .transpose()?,
            n_jobs,
            progress_bar,
        };
        let path = crate::molecular_io_calls::path_text(&directory)?;
        let report_path = report_path
            .as_deref()
            .map(crate::molecular_io_calls::path_text)
            .transpose()?;
        py.detach(|| {
            self.inner.write_sdf_files_with_params(
                path,
                params.map(|v| &v.inner).unwrap_or(&resolved),
                filenames.as_deref(),
                report_path,
            )
        })
        .map(|inner| crate::batch_errors::BatchExportReport { inner })
        .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[pyo3(signature=(directory,params,filenames,report_path))]
    fn write_sdf_files_with_params(
        &self,
        py: Python<'_>,
        directory: std::path::PathBuf,
        params: &crate::io_batch_params::BatchExportParams,
        filenames: Option<Vec<Option<String>>>,
        report_path: Option<std::path::PathBuf>,
    ) -> PyResult<crate::batch_errors::BatchExportReport> {
        self.write_sdf_files(
            py,
            directory,
            Some(params),
            None,
            None,
            None,
            None,
            filenames,
            report_path,
        )
    }
    #[pyo3(signature=(directory,options=None,*,format=None,width=None,height=None,execution=None,filenames=None,report_path=None))]
    fn write_images(
        &self,
        py: Python<'_>,
        directory: std::path::PathBuf,
        options: Option<&crate::io_batch_params::BatchImageParams>,
        format: Option<String>,
        width: Option<u32>,
        height: Option<u32>,
        execution: Option<&crate::io_batch_params::BatchParams>,
        filenames: Option<Vec<Option<String>>>,
        report_path: Option<std::path::PathBuf>,
    ) -> PyResult<crate::batch_errors::BatchExportReport> {
        if options.is_some()
            && (format.is_some()
                || width.is_some()
                || height.is_some()
                || execution.is_some()
                || filenames.is_some()
                || report_path.is_some())
        {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Params and image keywords are mutually exclusive",
            ));
        }
        let d = ck::BatchImageParams::default();
        let resolved = ck::BatchImageParams {
            format: format.unwrap_or(d.format),
            width: width.unwrap_or(d.width),
            height: height.unwrap_or(d.height),
            execution: execution.map(|v| v.inner.clone()).unwrap_or(d.execution),
            filenames,
            report_path,
        };
        py.detach(|| {
            self.inner.write_images_with_params(
                &directory,
                options.map(|v| &v.inner).unwrap_or(&resolved),
            )
        })
        .map(|inner| crate::batch_errors::BatchExportReport { inner })
        .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    fn write_images_with_params(
        &self,
        py: Python<'_>,
        directory: std::path::PathBuf,
        options: &crate::io_batch_params::BatchImageParams,
    ) -> PyResult<crate::batch_errors::BatchExportReport> {
        self.write_images(
            py,
            directory,
            Some(options),
            None,
            None,
            None,
            None,
            None,
            None,
        )
    }

    #[pyo3(signature=(options,atom_invariants,bond_invariants,call,params))]
    fn fingerprint_morgan_list_with_generator_params(
        &self,
        py: Python<'_>,
        options: &crate::fingerprint_configuration::MorganParams,
        atom_invariants: Option<&crate::morgan_configuration::MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&crate::fingerprint_providers::MorganBondInvariantsGenerator>,
        call: &crate::params::MorganCallParams,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<crate::native::Fingerprint>>> {
        params
            .with_callback(py, |runtime| {
                self.inner.fingerprint_morgan_list_with_generator_params(
                    &options.inner,
                    atom_invariants.map(|v| &v.inner),
                    bond_invariants.map(|v| &v.inner),
                    &call.inner,
                    runtime,
                )
            })?
            .map(|v| {
                v.into_iter()
                    .map(|v| v.map(|inner| crate::native::Fingerprint { inner }))
                    .collect()
            })
            .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }

    #[pyo3(signature=(options,atom_invariants,bond_invariants,call,collect_additional_output,params))]
    fn fingerprint_morgan_with_output_list_with_generator_params(
        &self,
        py: Python<'_>,
        options: &crate::fingerprint_configuration::MorganParams,
        atom_invariants: Option<&crate::morgan_configuration::MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&crate::fingerprint_providers::MorganBondInvariantsGenerator>,
        call: &crate::params::MorganCallParams,
        collect_additional_output: bool,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<crate::batch_fingerprint_values::BatchFingerprintOutput>>> {
        params
            .with_callback(py, |runtime| {
                self.inner
                    .fingerprint_morgan_with_output_list_with_generator_params(
                        &options.inner,
                        atom_invariants.map(|v| &v.inner),
                        bond_invariants.map(|v| &v.inner),
                        &call.inner,
                        collect_additional_output,
                        runtime,
                    )
            })?
            .map(|v| {
                v.into_iter()
                    .map(|v| {
                        v.map(|inner| {
                            crate::batch_fingerprint_values::BatchFingerprintOutput::from_rust(
                                inner,
                            )
                        })
                    })
                    .collect()
            })
            .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }

    #[pyo3(signature=(options=None,params=None,*,n_jobs=None,progress_bar=None))]
    fn to_smiles_list(
        &self,
        py: Python<'_>,
        options: Option<&crate::chemical_params::SmilesWriteParams>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<String>>> {
        let defaults = ck::SmilesWriteParams::default();
        let values = crate::batch_calls::query(py, params, n_jobs, progress_bar, |runtime| {
            self.inner
                .to_smiles_list_with_params(options.map(|v| &v.inner).unwrap_or(&defaults), runtime)
        })?;
        values
            .into_iter()
            .map(|v| {
                v.map(|v| crate::property_values::decode_source_text(py, &v))
                    .transpose()
            })
            .collect()
    }
    fn to_smiles_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::chemical_params::SmilesWriteParams,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<String>>> {
        self.to_smiles_list(py, Some(options), Some(params), None, None)
    }
    #[pyo3(signature=(params=None,*,n_jobs=None,progress_bar=None))]
    fn dg_bounds_matrix_list(
        &self,
        py: Python<'_>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<Vec<Vec<f64>>>>> {
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |p| {
            self.inner.dg_bounds_matrix_list_with_params(p)
        })
    }
    fn dg_bounds_matrix_list_with_params(
        &self,
        py: Python<'_>,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<Vec<Vec<f64>>>>> {
        self.dg_bounds_matrix_list(py, Some(params), None, None)
    }
    #[pyo3(signature=(width,height,params=None,*,n_jobs=None,progress_bar=None))]
    fn to_svg_list(
        &self,
        py: Python<'_>,
        width: u32,
        height: u32,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<String>>> {
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |p| {
            self.inner.to_svg_list_with_params(width, height, p)
        })
    }
    fn to_svg_list_with_params(
        &self,
        py: Python<'_>,
        width: u32,
        height: u32,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<String>>> {
        self.to_svg_list(py, width, height, Some(params), None, None)
    }

    #[pyo3(signature=(options=None,params=None,*,n_jobs=None,progress_bar=None))]
    fn fingerprint_morgan_list(
        &self,
        py: Python<'_>,
        options: Option<&crate::morgan_configuration::MorganFingerprintParams>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<crate::native::Fingerprint>>> {
        let defaults = ck::MorganFingerprintParams::default();
        let selected = options.map(|v| &v.inner).unwrap_or(&defaults);
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |runtime| {
            self.inner
                .fingerprint_morgan_list_with_params(selected, runtime)
        })
        .map(|v| {
            v.into_iter()
                .map(|v| v.map(|inner| crate::native::Fingerprint { inner }))
                .collect()
        })
    }
    fn fingerprint_morgan_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::morgan_configuration::MorganFingerprintParams,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<crate::native::Fingerprint>>> {
        self.fingerprint_morgan_list(py, Some(options), Some(params), None, None)
    }

    #[pyo3(signature=(options=None,params=None,*,n_jobs=None,progress_bar=None))]
    fn fingerprint_atom_pair_list(
        &self,
        py: Python<'_>,
        options: Option<&crate::fingerprint_configuration::AtomPairFingerprintParams>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<crate::native::Fingerprint>>> {
        let defaults = ck::AtomPairFingerprintParams::default();
        let selected = options.map(|v| &v.inner).unwrap_or(&defaults);
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |runtime| {
            self.inner
                .fingerprint_atom_pair_list_with_params(selected, runtime)
        })
        .map(|v| {
            v.into_iter()
                .map(|v| v.map(|inner| crate::native::Fingerprint { inner }))
                .collect()
        })
    }
    fn fingerprint_atom_pair_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::fingerprint_configuration::AtomPairFingerprintParams,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<crate::native::Fingerprint>>> {
        self.fingerprint_atom_pair_list(py, Some(options), Some(params), None, None)
    }

    #[pyo3(signature=(options=None,params=None,*,n_jobs=None,progress_bar=None))]
    fn fingerprint_atom_pair_sparse_count_list(
        &self,
        py: Python<'_>,
        options: Option<&crate::fingerprint_configuration::AtomPairFingerprintParams>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<crate::fingerprint_values::SparseCountFingerprint>>> {
        let defaults = ck::AtomPairFingerprintParams::default();
        let selected = options.map(|v| &v.inner).unwrap_or(&defaults);
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |runtime| {
            self.inner
                .fingerprint_atom_pair_sparse_count_list_with_params(selected, runtime)
        })
        .map(|v| {
            v.into_iter()
                .map(|v| v.map(|inner| crate::fingerprint_values::SparseCountFingerprint { inner }))
                .collect()
        })
    }
    fn fingerprint_atom_pair_sparse_count_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::fingerprint_configuration::AtomPairFingerprintParams,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<crate::fingerprint_values::SparseCountFingerprint>>> {
        self.fingerprint_atom_pair_sparse_count_list(py, Some(options), Some(params), None, None)
    }

    #[pyo3(signature=(options=None,params=None,*,n_jobs=None,progress_bar=None))]
    fn fingerprint_atom_pair_count_list(
        &self,
        py: Python<'_>,
        options: Option<&crate::fingerprint_configuration::AtomPairFingerprintParams>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<crate::fingerprint_values::SparseCountFingerprint32>>> {
        let defaults = ck::AtomPairFingerprintParams::default();
        let selected = options.map(|v| &v.inner).unwrap_or(&defaults);
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |runtime| {
            self.inner
                .fingerprint_atom_pair_count_list_with_params(selected, runtime)
        })
        .map(|v| {
            v.into_iter()
                .map(|v| {
                    v.map(|inner| crate::fingerprint_values::SparseCountFingerprint32 { inner })
                })
                .collect()
        })
    }
    fn fingerprint_atom_pair_count_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::fingerprint_configuration::AtomPairFingerprintParams,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<crate::fingerprint_values::SparseCountFingerprint32>>> {
        self.fingerprint_atom_pair_count_list(py, Some(options), Some(params), None, None)
    }

    #[pyo3(signature=(options=None,params=None,*,n_jobs=None,progress_bar=None))]
    fn fingerprint_atom_pair_sparse_bits_list(
        &self,
        py: Python<'_>,
        options: Option<&crate::fingerprint_configuration::AtomPairFingerprintParams>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<crate::fingerprint_values::SparseBitFingerprint>>> {
        let defaults = ck::AtomPairFingerprintParams::default();
        let selected = options.map(|v| &v.inner).unwrap_or(&defaults);
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |runtime| {
            self.inner
                .fingerprint_atom_pair_sparse_bits_list_with_params(selected, runtime)
        })
        .map(|v| {
            v.into_iter()
                .map(|v| v.map(|inner| crate::fingerprint_values::SparseBitFingerprint { inner }))
                .collect()
        })
    }
    fn fingerprint_atom_pair_sparse_bits_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::fingerprint_configuration::AtomPairFingerprintParams,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<crate::fingerprint_values::SparseBitFingerprint>>> {
        self.fingerprint_atom_pair_sparse_bits_list(py, Some(options), Some(params), None, None)
    }

    #[pyo3(signature=(options=None,params=None,*,n_jobs=None,progress_bar=None))]
    fn fingerprint_topological_torsion_list(
        &self,
        py: Python<'_>,
        options: Option<&crate::fingerprint_configuration::TopologicalTorsionFingerprintParams>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<crate::native::Fingerprint>>> {
        let defaults = ck::TopologicalTorsionFingerprintParams::default();
        let selected = options.map(|v| &v.inner).unwrap_or(&defaults);
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |runtime| {
            self.inner
                .fingerprint_topological_torsion_list_with_params(selected, runtime)
        })
        .map(|v| {
            v.into_iter()
                .map(|v| v.map(|inner| crate::native::Fingerprint { inner }))
                .collect()
        })
    }
    fn fingerprint_topological_torsion_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::fingerprint_configuration::TopologicalTorsionFingerprintParams,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<crate::native::Fingerprint>>> {
        self.fingerprint_topological_torsion_list(py, Some(options), Some(params), None, None)
    }

    #[pyo3(signature=(options=None,params=None,*,n_jobs=None,progress_bar=None))]
    fn fingerprint_layered_list(
        &self,
        py: Python<'_>,
        options: Option<&crate::collection_params::LayeredFingerprintParams>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<crate::native::Fingerprint>>> {
        let defaults = ck::LayeredFingerprintParams::default();
        let selected = options.map(|v| &v.inner).unwrap_or(&defaults);
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |runtime| {
            self.inner
                .fingerprint_layered_list_with_params(selected, runtime)
        })
        .map(|v| {
            v.into_iter()
                .map(|v| v.map(|inner| crate::native::Fingerprint { inner }))
                .collect()
        })
    }
    fn fingerprint_layered_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::collection_params::LayeredFingerprintParams,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<crate::native::Fingerprint>>> {
        self.fingerprint_layered_list(py, Some(options), Some(params), None, None)
    }

    #[pyo3(signature=(options=None,params=None,*,n_jobs=None,progress_bar=None))]
    fn fingerprint_layered_with_output_list(
        &self,
        py: Python<'_>,
        options: Option<&crate::collection_params::LayeredFingerprintParams>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<crate::ordered_results::LayeredFingerprintResult>>> {
        let defaults = ck::LayeredFingerprintParams::default();
        let selected = options.map(|v| &v.inner).unwrap_or(&defaults);
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |runtime| {
            self.inner
                .fingerprint_layered_with_output_list_with_params(selected, runtime)
        })
        .map(|v| {
            v.into_iter()
                .map(|v| v.map(|inner| crate::ordered_results::LayeredFingerprintResult { inner }))
                .collect()
        })
    }
    fn fingerprint_layered_with_output_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::collection_params::LayeredFingerprintParams,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<crate::ordered_results::LayeredFingerprintResult>>> {
        self.fingerprint_layered_with_output_list(py, Some(options), Some(params), None, None)
    }

    #[pyo3(signature=(options=None,params=None,*,n_jobs=None,progress_bar=None))]
    fn fingerprint_pattern_list(
        &self,
        py: Python<'_>,
        options: Option<&crate::params::PatternFingerprintParams>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<crate::native::Fingerprint>>> {
        let defaults = ck::PatternFingerprintParams::default();
        let selected = options.map(|v| &v.inner).unwrap_or(&defaults);
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |runtime| {
            self.inner
                .fingerprint_pattern_list_with_params(selected, runtime)
        })
        .map(|v| {
            v.into_iter()
                .map(|v| v.map(|inner| crate::native::Fingerprint { inner }))
                .collect()
        })
    }
    fn fingerprint_pattern_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::params::PatternFingerprintParams,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<crate::native::Fingerprint>>> {
        self.fingerprint_pattern_list(py, Some(options), Some(params), None, None)
    }

    #[pyo3(signature=(options=None,params=None,*,n_jobs=None,progress_bar=None,collect_additional_output=true))]
    fn fingerprint_morgan_with_output_list(
        &self,
        py: Python<'_>,
        options: Option<&crate::morgan_configuration::MorganFingerprintParams>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
        collect_additional_output: bool,
    ) -> PyResult<Vec<Option<crate::batch_fingerprint_values::BatchFingerprintOutput>>> {
        let defaults = ck::MorganFingerprintParams::default();
        let selected = options.map(|v| &v.inner).unwrap_or(&defaults);
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |runtime| {
            self.inner.fingerprint_morgan_with_output_list_with_params(
                selected,
                collect_additional_output,
                runtime,
            )
        })
        .map(|v| {
            v.into_iter()
                .map(|v| {
                    v.map(|inner| {
                        crate::batch_fingerprint_values::BatchFingerprintOutput::from_rust(inner)
                    })
                })
                .collect()
        })
    }
    fn fingerprint_morgan_with_output_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::morgan_configuration::MorganFingerprintParams,
        collect_additional_output: bool,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<crate::batch_fingerprint_values::BatchFingerprintOutput>>> {
        self.fingerprint_morgan_with_output_list(
            py,
            Some(options),
            Some(params),
            None,
            None,
            collect_additional_output,
        )
    }

    #[pyo3(signature=(options=None,params=None,*,n_jobs=None,progress_bar=None,collect_additional_output=true))]
    fn fingerprint_atom_pair_with_output_list(
        &self,
        py: Python<'_>,
        options: Option<&crate::fingerprint_configuration::AtomPairFingerprintParams>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
        collect_additional_output: bool,
    ) -> PyResult<Vec<Option<crate::batch_fingerprint_values::BatchFingerprintOutput>>> {
        let defaults = ck::AtomPairFingerprintParams::default();
        let selected = options.map(|v| &v.inner).unwrap_or(&defaults);
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |runtime| {
            self.inner
                .fingerprint_atom_pair_with_output_list_with_params(
                    selected,
                    collect_additional_output,
                    runtime,
                )
        })
        .map(|v| {
            v.into_iter()
                .map(|v| {
                    v.map(|inner| {
                        crate::batch_fingerprint_values::BatchFingerprintOutput::from_rust(inner)
                    })
                })
                .collect()
        })
    }
    fn fingerprint_atom_pair_with_output_list_with_params(
        &self,
        py: Python<'_>,
        options: &crate::fingerprint_configuration::AtomPairFingerprintParams,
        collect_additional_output: bool,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<crate::batch_fingerprint_values::BatchFingerprintOutput>>> {
        self.fingerprint_atom_pair_with_output_list(
            py,
            Some(options),
            Some(params),
            None,
            None,
            collect_additional_output,
        )
    }

    #[staticmethod]
    #[pyo3(signature=(smiles,parse=None,params=None,*,errors=None,n_jobs=None,progress_bar=None))]
    pub(crate) fn from_smiles_list(
        py: Python<'_>,
        smiles: Vec<String>,
        parse: Option<&SmilesParseParams>,
        params: Option<&crate::io_batch_params::BatchParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        if params.is_some() && (errors.is_some() || n_jobs.is_some() || progress_bar.is_some()) {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Params and runtime keywords are mutually exclusive",
            ));
        }
        let runtime;
        let parse_default = ck::SmilesParseParams::default();
        let unconfigured = parse.is_none()
            && params.is_none()
            && errors.is_none()
            && n_jobs.is_none()
            && progress_bar.is_none();
        let exec = match params {
            Some(p) => &p.inner,
            None => {
                runtime = crate::batch_calls::execution(errors, n_jobs, progress_bar)?;
                &runtime.inner
            }
        };
        py.detach(|| {
            if unconfigured {
                ck::MoleculeBatch::from_smiles_list(&smiles)
            } else {
                ck::MoleculeBatch::from_smiles_list_with_params(
                    &smiles,
                    parse.map(|v| &v.inner).unwrap_or(&parse_default),
                    exec,
                )
            }
        })
        .map(|inner| Self { inner })
        .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[staticmethod]
    fn from_smiles_list_with_params(
        py: Python<'_>,
        smiles: Vec<String>,
        parse: &SmilesParseParams,
        params: &crate::io_batch_params::BatchParams,
    ) -> PyResult<Self> {
        Self::from_smiles_list(py, smiles, Some(parse), Some(params), None, None, None)
    }
    #[staticmethod]
    #[pyo3(signature=(text,read=None,*,mode=None,n_jobs=None,progress_bar=None))]
    pub(crate) fn from_sdf_records(
        py: Python<'_>,
        text: String,
        read: Option<&crate::selectors::SdfReadParams>,
        mode: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let unconfigured =
            read.is_none() && mode.is_none() && n_jobs.is_none() && progress_bar.is_none();
        let defaults = ck::SdfReadParams::default();
        let selected = mode
            .map(crate::io_batch_params::BatchErrorModeInput::into_rust)
            .transpose()?
            .unwrap_or_default();
        let input = &text;
        py.detach(|| {
            if unconfigured {
                ck::MoleculeBatch::from_sdf_records(input)
            } else {
                ck::MoleculeBatch::from_sdf_records_with_params(
                    input,
                    read.map(|v| &v.inner).unwrap_or(&defaults),
                    selected,
                    n_jobs,
                )
            }
        })
        .map(|inner| Self { inner })
        .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[staticmethod]

    fn from_sdf_records_with_params(
        py: Python<'_>,
        text: String,
        read: &crate::selectors::SdfReadParams,
        mode: crate::io_batch_params::BatchErrorModeInput,
        n_jobs: Option<usize>,
    ) -> PyResult<Self> {
        Self::from_sdf_records(py, text, Some(read), Some(mode), n_jobs, None)
    }
    #[staticmethod]
    #[pyo3(signature=(path,read=None,*,mode=None,n_jobs=None,progress_bar=None))]
    pub(crate) fn read_sdf(
        py: Python<'_>,
        path: std::path::PathBuf,
        read: Option<&crate::selectors::SdfReadParams>,
        mode: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let unconfigured =
            read.is_none() && mode.is_none() && n_jobs.is_none() && progress_bar.is_none();
        let defaults = ck::SdfReadParams::default();
        let selected = mode
            .map(crate::io_batch_params::BatchErrorModeInput::into_rust)
            .transpose()?
            .unwrap_or_default();
        let input = crate::molecular_io_calls::path_text(&path)?;
        py.detach(|| {
            if unconfigured {
                ck::MoleculeBatch::read_sdf(input)
            } else {
                ck::MoleculeBatch::read_sdf_with_params(
                    input,
                    read.map(|v| &v.inner).unwrap_or(&defaults),
                    selected,
                    n_jobs,
                    progress_bar.unwrap_or(false),
                )
            }
        })
        .map(|inner| Self { inner })
        .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[staticmethod]
    #[pyo3(signature=(path,read,mode,n_jobs,progress_bar))]
    fn read_sdf_with_params(
        py: Python<'_>,
        path: std::path::PathBuf,
        read: &crate::selectors::SdfReadParams,
        mode: crate::io_batch_params::BatchErrorModeInput,
        n_jobs: Option<usize>,
        progress_bar: bool,
    ) -> PyResult<Self> {
        Self::read_sdf(py, path, Some(read), Some(mode), n_jobs, Some(progress_bar))
    }

    #[pyo3(signature=(options=None,params=None,*,errors=None,n_jobs=None,progress_bar=None))]
    fn sanitize(
        &self,
        py: Python<'_>,
        options: Option<&crate::sanitize_values::SanitizeParams>,
        params: Option<&crate::io_batch_params::BatchParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        if params.is_some() && (errors.is_some() || n_jobs.is_some() || progress_bar.is_some()) {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Params and runtime keywords are mutually exclusive",
            ));
        }
        if errors.is_some() || n_jobs.is_some() || progress_bar.is_some() {
            self._batch_sanitize_keywords(py, options, errors, n_jobs, progress_bar)
        } else {
            self._batch_sanitize_params(py, options, params)
        }
    }
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.sanitize"}"#]
    #[pyo3(signature=(options=None,params=None))]
    fn _batch_sanitize_params(
        &self,
        py: Python<'_>,
        options: Option<&crate::sanitize_values::SanitizeParams>,
        params: Option<&crate::io_batch_params::BatchParams>,
    ) -> PyResult<Self> {
        let defaults = ck::SanitizeParams::default();
        let execution = ck::BatchParams::default();
        py.detach(|| {
            if options.is_none() && params.is_none() {
                self.inner.sanitize()
            } else {
                self.inner.sanitize_with_params(
                    options.map(|v| &v.inner).unwrap_or(&defaults),
                    params.map(|v| &v.inner).unwrap_or(&execution),
                )
            }
        })
        .map(|inner| Self { inner })
        .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.sanitize"}"#]
    #[pyo3(signature=(options=None,*,errors=None,n_jobs=None,progress_bar=None))]
    fn _batch_sanitize_keywords(
        &self,
        py: Python<'_>,
        options: Option<&crate::sanitize_values::SanitizeParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let execution = crate::batch_calls::execution(errors, n_jobs, progress_bar)?;
        self._batch_sanitize_params(py, options, Some(&execution))
    }
    fn sanitize_with_params(
        &self,
        py: Python<'_>,
        options: &crate::sanitize_values::SanitizeParams,
        params: &crate::io_batch_params::BatchParams,
    ) -> PyResult<Self> {
        self._batch_sanitize_params(py, Some(options), Some(params))
    }
    #[pyo3(signature=(options=None,params=None,*,errors=None,n_jobs=None,progress_bar=None))]
    fn with_hydrogens(
        &self,
        py: Python<'_>,
        options: Option<&crate::chemical_params::AddHsParams>,
        params: Option<&crate::io_batch_params::BatchParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        if params.is_some() && (errors.is_some() || n_jobs.is_some() || progress_bar.is_some()) {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Params and runtime keywords are mutually exclusive",
            ));
        }
        if errors.is_some() || n_jobs.is_some() || progress_bar.is_some() {
            self._batch_with_hydrogens_keywords(py, options, errors, n_jobs, progress_bar)
        } else {
            self._batch_with_hydrogens_params(py, options, params)
        }
    }
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.with_hydrogens"}"#]
    #[pyo3(signature=(options=None,params=None))]
    fn _batch_with_hydrogens_params(
        &self,
        py: Python<'_>,
        options: Option<&crate::chemical_params::AddHsParams>,
        params: Option<&crate::io_batch_params::BatchParams>,
    ) -> PyResult<Self> {
        let defaults = ck::AddHsParams::default();
        let execution = ck::BatchParams::default();
        py.detach(|| {
            if options.is_none() && params.is_none() {
                self.inner.with_hydrogens()
            } else {
                self.inner.with_hydrogens_with_params(
                    options.map(|v| &v.inner).unwrap_or(&defaults),
                    params.map(|v| &v.inner).unwrap_or(&execution),
                )
            }
        })
        .map(|inner| Self { inner })
        .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.with_hydrogens"}"#]
    #[pyo3(signature=(options=None,*,errors=None,n_jobs=None,progress_bar=None))]
    fn _batch_with_hydrogens_keywords(
        &self,
        py: Python<'_>,
        options: Option<&crate::chemical_params::AddHsParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let execution = crate::batch_calls::execution(errors, n_jobs, progress_bar)?;
        self._batch_with_hydrogens_params(py, options, Some(&execution))
    }
    fn with_hydrogens_with_params(
        &self,
        py: Python<'_>,
        options: &crate::chemical_params::AddHsParams,
        params: &crate::io_batch_params::BatchParams,
    ) -> PyResult<Self> {
        self._batch_with_hydrogens_params(py, Some(options), Some(params))
    }
    #[pyo3(signature=(options=None,params=None,*,errors=None,n_jobs=None,progress_bar=None))]
    fn without_hydrogens(
        &self,
        py: Python<'_>,
        options: Option<&crate::params::RemoveHsParams>,
        params: Option<&crate::io_batch_params::BatchParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        if params.is_some() && (errors.is_some() || n_jobs.is_some() || progress_bar.is_some()) {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Params and runtime keywords are mutually exclusive",
            ));
        }
        if errors.is_some() || n_jobs.is_some() || progress_bar.is_some() {
            self._batch_without_hydrogens_keywords(py, options, errors, n_jobs, progress_bar)
        } else {
            self._batch_without_hydrogens_params(py, options, params)
        }
    }
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.without_hydrogens"}"#]
    #[pyo3(signature=(options=None,params=None))]
    fn _batch_without_hydrogens_params(
        &self,
        py: Python<'_>,
        options: Option<&crate::params::RemoveHsParams>,
        params: Option<&crate::io_batch_params::BatchParams>,
    ) -> PyResult<Self> {
        let defaults = ck::RemoveHsParams::default();
        let execution = ck::BatchParams::default();
        py.detach(|| {
            if options.is_none() && params.is_none() {
                self.inner.without_hydrogens()
            } else {
                self.inner.without_hydrogens_with_params(
                    options.map(|v| &v.inner).unwrap_or(&defaults),
                    params.map(|v| &v.inner).unwrap_or(&execution),
                )
            }
        })
        .map(|inner| Self { inner })
        .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.without_hydrogens"}"#]
    #[pyo3(signature=(options=None,*,errors=None,n_jobs=None,progress_bar=None))]
    fn _batch_without_hydrogens_keywords(
        &self,
        py: Python<'_>,
        options: Option<&crate::params::RemoveHsParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let execution = crate::batch_calls::execution(errors, n_jobs, progress_bar)?;
        self._batch_without_hydrogens_params(py, options, Some(&execution))
    }
    fn without_hydrogens_with_params(
        &self,
        py: Python<'_>,
        options: &crate::params::RemoveHsParams,
        params: &crate::io_batch_params::BatchParams,
    ) -> PyResult<Self> {
        self._batch_without_hydrogens_params(py, Some(options), Some(params))
    }
    #[pyo3(signature=(options=None,params=None,*,errors=None,n_jobs=None,progress_bar=None))]
    fn with_kekulized_bonds(
        &self,
        py: Python<'_>,
        options: Option<&crate::params::KekulizeParams>,
        params: Option<&crate::io_batch_params::BatchParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        if params.is_some() && (errors.is_some() || n_jobs.is_some() || progress_bar.is_some()) {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Params and runtime keywords are mutually exclusive",
            ));
        }
        if errors.is_some() || n_jobs.is_some() || progress_bar.is_some() {
            self._batch_with_kekulized_bonds_keywords(py, options, errors, n_jobs, progress_bar)
        } else {
            self._batch_with_kekulized_bonds_params(py, options, params)
        }
    }
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.with_kekulized_bonds"}"#]
    #[pyo3(signature=(options=None,params=None))]
    fn _batch_with_kekulized_bonds_params(
        &self,
        py: Python<'_>,
        options: Option<&crate::params::KekulizeParams>,
        params: Option<&crate::io_batch_params::BatchParams>,
    ) -> PyResult<Self> {
        let defaults = ck::KekulizeParams::default();
        let execution = ck::BatchParams::default();
        py.detach(|| {
            if options.is_none() && params.is_none() {
                self.inner.with_kekulized_bonds()
            } else {
                self.inner.with_kekulized_bonds_with_params(
                    options.map(|v| &v.inner).unwrap_or(&defaults),
                    params.map(|v| &v.inner).unwrap_or(&execution),
                )
            }
        })
        .map(|inner| Self { inner })
        .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.with_kekulized_bonds"}"#]
    #[pyo3(signature=(options=None,*,errors=None,n_jobs=None,progress_bar=None))]
    fn _batch_with_kekulized_bonds_keywords(
        &self,
        py: Python<'_>,
        options: Option<&crate::params::KekulizeParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let execution = crate::batch_calls::execution(errors, n_jobs, progress_bar)?;
        self._batch_with_kekulized_bonds_params(py, options, Some(&execution))
    }
    fn with_kekulized_bonds_with_params(
        &self,
        py: Python<'_>,
        options: &crate::params::KekulizeParams,
        params: &crate::io_batch_params::BatchParams,
    ) -> PyResult<Self> {
        self._batch_with_kekulized_bonds_params(py, Some(options), Some(params))
    }
    #[pyo3(signature=(options=None,params=None,*,errors=None,n_jobs=None,progress_bar=None))]
    fn with_2d_coordinates(
        &self,
        py: Python<'_>,
        options: Option<&crate::optional_configuration::Coordinate2DParams>,
        params: Option<&crate::io_batch_params::BatchParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        if params.is_some() && (errors.is_some() || n_jobs.is_some() || progress_bar.is_some()) {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Params and runtime keywords are mutually exclusive",
            ));
        }
        if errors.is_some() || n_jobs.is_some() || progress_bar.is_some() {
            self._batch_with_2d_coordinates_keywords(py, options, errors, n_jobs, progress_bar)
        } else {
            self._batch_with_2d_coordinates_params(py, options, params)
        }
    }
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.with_2d_coordinates"}"#]
    #[pyo3(signature=(options=None,params=None))]
    fn _batch_with_2d_coordinates_params(
        &self,
        py: Python<'_>,
        options: Option<&crate::optional_configuration::Coordinate2DParams>,
        params: Option<&crate::io_batch_params::BatchParams>,
    ) -> PyResult<Self> {
        let defaults = ck::Coordinate2DParams::default();
        let execution = ck::BatchParams::default();
        py.detach(|| {
            if options.is_none() && params.is_none() {
                self.inner.with_2d_coordinates()
            } else {
                self.inner.with_2d_coordinates_with_params(
                    options.map(|v| &v.inner).unwrap_or(&defaults),
                    params.map(|v| &v.inner).unwrap_or(&execution),
                )
            }
        })
        .map(|inner| Self { inner })
        .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.with_2d_coordinates"}"#]
    #[pyo3(signature=(options=None,*,errors=None,n_jobs=None,progress_bar=None))]
    fn _batch_with_2d_coordinates_keywords(
        &self,
        py: Python<'_>,
        options: Option<&crate::optional_configuration::Coordinate2DParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let execution = crate::batch_calls::execution(errors, n_jobs, progress_bar)?;
        self._batch_with_2d_coordinates_params(py, options, Some(&execution))
    }
    fn with_2d_coordinates_with_params(
        &self,
        py: Python<'_>,
        options: &crate::optional_configuration::Coordinate2DParams,
        params: &crate::io_batch_params::BatchParams,
    ) -> PyResult<Self> {
        self._batch_with_2d_coordinates_params(py, Some(options), Some(params))
    }

    fn records(&self, py: Python<'_>) -> PyResult<Vec<crate::batch_calls::Record>> {
        self.inner
            .records()
            .iter()
            .map(|v| crate::batch_calls::record(py, v))
            .collect()
    }
    fn get(&self, py: Python<'_>, index: usize) -> PyResult<Option<crate::batch_calls::Record>> {
        self.inner
            .get(index)
            .map(|v| crate::batch_calls::record(py, v))
            .transpose()
    }
    fn errors(&self, py: Python<'_>) -> PyResult<Vec<Py<crate::batch_errors::BatchError>>> {
        self.inner
            .errors()
            .into_iter()
            .map(|v| {
                Ok(crate::batch_errors::BatchError::from_rust(v)
                    .value(py)
                    .extract()?)
            })
            .collect()
    }
    #[pyo3(signature=(n_jobs=None))]
    fn with_parallel_jobs(&self, py: Python<'_>, n_jobs: Option<usize>) -> PyResult<Self> {
        py.detach(|| self.inner.clone().with_parallel_jobs(n_jobs))
            .map(|inner| Self { inner })
            .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[pyo3(signature=(progress_bar=None))]
    fn with_progress_bar(&self, progress_bar: Option<bool>) -> Self {
        Self {
            inner: self.inner.clone().with_progress_bar(progress_bar),
        }
    }
    #[staticmethod]
    #[pyo3(signature=(records,mode=None))]
    pub(crate) fn from_records(
        py: Python<'_>,
        records: Vec<crate::batch_calls::Record>,
        mode: Option<crate::io_batch_params::BatchErrorModeInput>,
    ) -> PyResult<Self> {
        let mode = mode
            .map(crate::io_batch_params::BatchErrorModeInput::into_rust)
            .transpose()?
            .unwrap_or_default();
        let records = records
            .into_iter()
            .map(|v| crate::batch_calls::native_record(py, v))
            .collect();
        py.detach(|| ck::MoleculeBatch::from_records(records, mode))
            .map(|inner| Self { inner })
            .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[staticmethod]
    #[pyo3(signature=(dataset,indices,mode=None))]
    fn from_dataset_indices(
        py: Python<'_>,
        dataset: &crate::sdf_sources::SdfDataset,
        indices: Vec<usize>,
        mode: Option<crate::io_batch_params::BatchErrorModeInput>,
    ) -> PyResult<Self> {
        let mode = mode
            .map(crate::io_batch_params::BatchErrorModeInput::into_rust)
            .transpose()?
            .unwrap_or_default();
        py.detach(|| ck::MoleculeBatch::from_dataset_indices(&dataset.inner, &indices, mode))
            .map(|inner| Self { inner })
            .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }

    #[new]
    fn new() -> Self {
        Self {
            inner: ck::MoleculeBatch::default(),
        }
    }
    fn __len__(&self) -> usize {
        self.inner.len()
    }
    fn len(&self) -> usize {
        self.inner.len()
    }
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    fn error_mode(&self) -> BatchErrorMode {
        self.inner.error_mode().into()
    }
    fn valid_mask(&self) -> Vec<bool> {
        self.inner.valid_mask()
    }
    fn invalid_mask(&self) -> Vec<bool> {
        self.inner.invalid_mask()
    }
    fn valid_count(&self) -> usize {
        self.inner.valid_count()
    }
    fn invalid_count(&self) -> usize {
        self.inner.invalid_count()
    }
    fn parallel_jobs(&self) -> Option<usize> {
        self.inner.parallel_jobs()
    }
    fn progress_bar(&self) -> Option<bool> {
        self.inner.progress_bar()
    }
    fn to_list(&self) -> Vec<Option<Molecule>> {
        self.inner
            .to_list()
            .into_iter()
            .map(|value| value.map(|inner| Molecule { inner }))
            .collect()
    }
    fn with_valid_records(&self) -> Self {
        Self {
            inner: self.inner.with_valid_records(),
        }
    }
}

#[cfg(test)]
mod batch_owner_tests {
    use super::*;
    use std::error::Error;

    #[test]
    fn actual_batch_queries_preserve_error_records_native_counts_and_cow() {
        let native = ck::MoleculeBatch::from_smiles_list_with_params(
            &["CCO".into(), "C(".into(), "N".into()],
            &ck::SmilesParseParams::default(),
            &ck::BatchParams {
                errors: Some(ck::BatchErrorMode::KeepErrors),
                n_jobs: Some(2),
                progress_bar: Some(false),
            },
        )
        .unwrap()
        .with_parallel_jobs(Some(2))
        .unwrap()
        .with_progress_bar(Some(false));
        assert_eq!(native.valid_mask(), vec![true, false, true]);
        let ck::BatchRecord::Error(error) = native.get(1).unwrap() else {
            panic!("expected actual native error record")
        };
        assert!(error.source().unwrap().is::<ck::SmilesError>());
        let display = error.to_string();
        let expected = native.to_list();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let batch = Py::new(py, MoleculeBatch { inner: native })?;
            assert_eq!(batch.bind(py).call_method0("len")?.extract::<usize>()?, 3);
            assert_eq!(batch.bind(py).len()?, 3);
            assert!(!batch.bind(py).call_method0("is_empty")?.extract::<bool>()?);
            assert_eq!(
                batch
                    .bind(py)
                    .call_method0("valid_count")?
                    .extract::<usize>()?,
                2
            );
            assert_eq!(
                batch
                    .bind(py)
                    .call_method0("invalid_count")?
                    .extract::<usize>()?,
                1
            );
            assert_eq!(
                batch
                    .bind(py)
                    .call_method0("valid_mask")?
                    .extract::<Vec<bool>>()?,
                vec![true, false, true]
            );
            assert_eq!(
                batch
                    .bind(py)
                    .call_method0("invalid_mask")?
                    .extract::<Vec<bool>>()?,
                vec![false, true, false]
            );
            assert!(
                batch
                    .bind(py)
                    .call_method0("error_mode")?
                    .extract::<BatchErrorMode>()?
                    == BatchErrorMode::KeepErrors
            );
            assert_eq!(
                batch
                    .bind(py)
                    .call_method0("parallel_jobs")?
                    .extract::<Option<usize>>()?,
                Some(2)
            );
            assert_eq!(
                batch
                    .bind(py)
                    .call_method0("progress_bar")?
                    .extract::<Option<bool>>()?,
                Some(false)
            );
            let filtered = batch.bind(py).call_method0("with_valid_records")?;
            assert_eq!(filtered.call_method0("valid_count")?.extract::<usize>()?, 2);
            assert_eq!(
                filtered.call_method0("invalid_count")?.extract::<usize>()?,
                0
            );
            let values = batch
                .bind(py)
                .call_method0("to_list")?
                .extract::<Vec<Option<Py<Molecule>>>>()?;
            for (actual, original) in values.iter().zip(&expected) {
                match (actual, original) {
                    (Some(actual), Some(original)) => {
                        assert_eq!(actual.borrow(py).inner, *original);
                        assert!(std::ptr::eq(
                            actual.borrow(py).inner.topology(),
                            original.topology()
                        ));
                    }
                    (None, None) => {}
                    _ => panic!("native batch row order changed"),
                }
            }
            let owner = batch.borrow(py);
            let ck::BatchRecord::Error(error) = owner.inner.get(1).unwrap() else {
                panic!("error record lost")
            };
            assert_eq!(error.to_string(), display);
            assert!(error.source().unwrap().is::<ck::SmilesError>());
            drop(owner);
            drop(filtered);
            drop(batch);
            assert_eq!(
                values[0]
                    .as_ref()
                    .unwrap()
                    .bind(py)
                    .call_method0("num_atoms")?
                    .extract::<usize>()?,
                3
            );
            Ok(())
        })
        .unwrap();
    }
}
