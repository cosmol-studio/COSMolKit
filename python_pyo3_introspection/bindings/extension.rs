//! SMILES configuration, module calls and the native molecule batch owner.
use crate::arguments::Argument;
use crate::parameters::{ConfigurationField, FieldInput, Scalar};
#[path = "batch_fingerprint_methods.rs"]
mod fingerprint_methods;
use crate::io_batch_params::BatchErrorMode;
use crate::native::Molecule;
use cosmolkit as ck;
use fingerprint_methods::batch_fingerprint_methods;
use pyo3::prelude::*;
use std::collections::BTreeMap;

ck::project_configuration!(ck::SmilesParseParams);

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
        defaults = Python::attach(|py| params.to_rust(py))?;
        &defaults
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
    let value = parse_keywords(
        sanitize,
        allow_cxsmiles,
        strict_cxsmiles,
        parse_name,
        remove_hs,
        skip_cleanup,
        debug_parse,
        replacements,
    );
    ck::Molecule::from_smiles_with_params(text, &value)
        .map(|inner| Molecule { inner })
        .map_err(crate::errors::SmilesError::from_rust)
}

pub(crate) fn parse_keywords(
    sanitize: Option<bool>,
    allow_cxsmiles: Option<bool>,
    strict_cxsmiles: Option<bool>,
    parse_name: Option<bool>,
    remove_hs: Option<bool>,
    skip_cleanup: Option<bool>,
    debug_parse: Option<bool>,
    replacements: Option<BTreeMap<String, String>>,
) -> ck::SmilesParseParams {
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
    value
}

fn batch_from_smiles(
    py: Python<'_>,
    smiles: Vec<String>,
    parse: &ck::SmilesParseParams,
    params: Option<&crate::io_batch_params::BatchParams>,
    errors: Option<crate::io_batch_params::BatchErrorModeInput>,
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
) -> PyResult<MoleculeBatch> {
    if params.is_some() && (errors.is_some() || n_jobs.is_some() || progress_bar.is_some()) {
        return Err(pyo3::exceptions::PyTypeError::new_err(
            "Params and runtime keywords are mutually exclusive",
        ));
    }
    let runtime = match params {
        Some(value) => value.to_rust(py)?,
        None => crate::batch_calls::execution(errors, n_jobs, progress_bar)?.to_rust(py)?,
    };
    py.detach(|| ck::MoleculeBatch::from_smiles_list_with_params(&smiles, parse, &runtime))
        .map(|inner| MoleculeBatch { inner })
        .map_err(crate::batch_errors::BatchValidationError::from_rust)
}

#[pyclass(module = "cosmolkit", frozen)]
pub struct MoleculeBatch {
    pub(crate) inner: ck::MoleculeBatch,
}

impl crate::indexing::Collection for MoleculeBatch {
    type Item = Option<Molecule>;
    type Selected = Self;
    fn length(&self) -> usize {
        self.inner.len()
    }
    fn item(&self, _py: Python<'_>, index: usize) -> PyResult<Self::Item> {
        Ok(match self.inner.get(index) {
            Some(ck::BatchRecord::Molecule(inner)) => Some(Molecule {
                inner: inner.clone(),
            }),
            Some(ck::BatchRecord::Error(_)) => None,
            None => return Err(pyo3::exceptions::PyIndexError::new_err(index)),
        })
    }
    fn select(&self, _py: Python<'_>, indices: &[usize]) -> PyResult<Self> {
        self.inner
            .select(indices)
            .map(|inner| Self { inner })
            .ok_or_else(|| pyo3::exceptions::PyIndexError::new_err("invalid selection"))
    }
}

batch_fingerprint_methods! {
[crate::morgan_configuration::MorganFingerprintParams, py, crate::fingerprint_calls::morgan_keywords(radius,include_chirality,use_bond_types,include_ring_membership,only_nonzero_invariants,include_redundant_environments,fp_size,count_simulation,bits_per_feature,count_bounds,from_atoms,ignore_atoms,custom_atom_invariants,custom_bond_invariants,conformer_id,invariants),
{radius: FieldInput<Scalar<u32>> = FieldInput::<Scalar<u32>>::Default,
include_chirality: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default,
use_bond_types: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default,
include_ring_membership: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default,
only_nonzero_invariants: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default,
include_redundant_environments: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default,
fp_size: FieldInput<Scalar<u32>> = FieldInput::<Scalar<u32>>::Default,
count_simulation: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default,
bits_per_feature: FieldInput<Scalar<u32>> = FieldInput::<Scalar<u32>>::Default,
count_bounds: FieldInput<Vec<u32>> = FieldInput::<Vec<u32>>::Default,
from_atoms: FieldInput<Option<Vec<u32>>> = FieldInput::<Option<Vec<u32>>>::Default,
ignore_atoms: FieldInput<Option<Vec<u32>>> = FieldInput::<Option<Vec<u32>>>::Default,
custom_atom_invariants: FieldInput<Option<Vec<u32>>> = FieldInput::<Option<Vec<u32>>>::Default,
custom_bond_invariants: FieldInput<Option<Vec<u32>>> = FieldInput::<Option<Vec<u32>>>::Default,
conformer_id: FieldInput<Scalar<i32>> = FieldInput::<Scalar<i32>>::Default,
invariants: FieldInput<ConfigurationField<crate::morgan_configuration::MorganInvariants>> = FieldInput::<ConfigurationField<crate::morgan_configuration::MorganInvariants>>::Default}
{
fingerprint_morgan, _fingerprint_morgan_params, _fingerprint_morgan_keywords, fingerprint_morgan_with_params,
crate::sequence::Output<ck::FingerprintBatch>, crate::sequence::Output, {};
fingerprint_morgan_with_output, _fingerprint_morgan_with_output_params, _fingerprint_morgan_with_output_keywords, fingerprint_morgan_with_output_with_params,
Vec<Option<crate::batch_fingerprint_values::BatchFingerprintOutput>>, |rows| rows.into_iter().map(|row| row.map(crate::batch_fingerprint_values::BatchFingerprintOutput::from_rust)).collect(), {collect_additional_output: bool = true};
}

crate::fingerprint_configuration::AtomPairFingerprintParams, py, crate::fingerprint_calls::atom_pair_keywords(min_distance,max_distance,include_chirality,use_2d,count_simulation,fp_size,bits_per_feature,count_bounds,from_atoms,ignore_atoms,custom_atom_invariants,custom_bond_invariants,conformer_id,atom_invariants_generator,use_legacy_stereo_perception),
{min_distance: FieldInput<Scalar<u32>> = FieldInput::<Scalar<u32>>::Default,
max_distance: FieldInput<Scalar<u32>> = FieldInput::<Scalar<u32>>::Default,
include_chirality: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default,
use_2d: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default,
count_simulation: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default,
fp_size: FieldInput<Scalar<u32>> = FieldInput::<Scalar<u32>>::Default,
bits_per_feature: FieldInput<Scalar<u32>> = FieldInput::<Scalar<u32>>::Default,
count_bounds: FieldInput<Vec<u32>> = FieldInput::<Vec<u32>>::Default,
from_atoms: FieldInput<Option<Vec<u32>>> = FieldInput::<Option<Vec<u32>>>::Default,
ignore_atoms: FieldInput<Option<Vec<u32>>> = FieldInput::<Option<Vec<u32>>>::Default,
custom_atom_invariants: FieldInput<Option<Vec<u32>>> = FieldInput::<Option<Vec<u32>>>::Default,
custom_bond_invariants: FieldInput<Option<Vec<u32>>> = FieldInput::<Option<Vec<u32>>>::Default,
conformer_id: FieldInput<Scalar<i32>> = FieldInput::<Scalar<i32>>::Default,
atom_invariants_generator: FieldInput<Option<ConfigurationField<crate::fingerprint_providers::AtomPairAtomInvariantsGenerator>>> = FieldInput::<Option<ConfigurationField<crate::fingerprint_providers::AtomPairAtomInvariantsGenerator>>>::Default,
use_legacy_stereo_perception: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default}
{
fingerprint_atom_pair, _fingerprint_atom_pair_params, _fingerprint_atom_pair_keywords, fingerprint_atom_pair_with_params,
crate::sequence::Output<ck::FingerprintBatch>, crate::sequence::Output, {};
fingerprint_atom_pair_count, _fingerprint_atom_pair_count_params, _fingerprint_atom_pair_count_keywords, fingerprint_atom_pair_count_with_params,
Vec<Option<crate::fingerprint_values::SparseCountFingerprint32>>, |rows| rows.into_iter().map(|row| row.map(|inner| crate::fingerprint_values::SparseCountFingerprint32 { inner })).collect(), {};
fingerprint_atom_pair_sparse, _fingerprint_atom_pair_sparse_params, _fingerprint_atom_pair_sparse_keywords, fingerprint_atom_pair_sparse_with_params,
Vec<Option<crate::fingerprint_values::SparseBitFingerprint>>, |rows| rows.into_iter().map(|row| row.map(|inner| crate::fingerprint_values::SparseBitFingerprint { inner })).collect(), {};
fingerprint_atom_pair_sparse_count, _fingerprint_atom_pair_sparse_count_params, _fingerprint_atom_pair_sparse_count_keywords, fingerprint_atom_pair_sparse_count_with_params,
Vec<Option<crate::fingerprint_values::SparseCountFingerprint>>, |rows| rows.into_iter().map(|row| row.map(|inner| crate::fingerprint_values::SparseCountFingerprint { inner })).collect(), {};
fingerprint_atom_pair_with_output, _fingerprint_atom_pair_with_output_params, _fingerprint_atom_pair_with_output_keywords, fingerprint_atom_pair_with_output_with_params,
Vec<Option<crate::batch_fingerprint_values::BatchFingerprintOutput>>, |rows| rows.into_iter().map(|row| row.map(crate::batch_fingerprint_values::BatchFingerprintOutput::from_rust)).collect(), {collect_additional_output: bool = true};
}

crate::fingerprint_configuration::TopologicalTorsionFingerprintParams, py, crate::fingerprint_calls::topological_torsion_keywords(torsion_atom_count,only_shortest_paths,include_chirality,count_simulation,fp_size,bits_per_feature,count_bounds,from_atoms,ignore_atoms,custom_atom_invariants,custom_bond_invariants,conformer_id,atom_invariants_generator,use_legacy_stereo_perception),
{torsion_atom_count: FieldInput<Scalar<u32>> = FieldInput::<Scalar<u32>>::Default,
only_shortest_paths: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default,
include_chirality: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default,
count_simulation: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default,
fp_size: FieldInput<Scalar<u32>> = FieldInput::<Scalar<u32>>::Default,
bits_per_feature: FieldInput<Scalar<u32>> = FieldInput::<Scalar<u32>>::Default,
count_bounds: FieldInput<Vec<u32>> = FieldInput::<Vec<u32>>::Default,
from_atoms: FieldInput<Option<Vec<u32>>> = FieldInput::<Option<Vec<u32>>>::Default,
ignore_atoms: FieldInput<Option<Vec<u32>>> = FieldInput::<Option<Vec<u32>>>::Default,
custom_atom_invariants: FieldInput<Option<Vec<u32>>> = FieldInput::<Option<Vec<u32>>>::Default,
custom_bond_invariants: FieldInput<Option<Vec<u32>>> = FieldInput::<Option<Vec<u32>>>::Default,
conformer_id: FieldInput<Scalar<i32>> = FieldInput::<Scalar<i32>>::Default,
atom_invariants_generator: FieldInput<Option<ConfigurationField<crate::fingerprint_providers::AtomPairAtomInvariantsGenerator>>> = FieldInput::<Option<ConfigurationField<crate::fingerprint_providers::AtomPairAtomInvariantsGenerator>>>::Default,
use_legacy_stereo_perception: FieldInput<Scalar<bool>> = FieldInput::<Scalar<bool>>::Default}
{
fingerprint_topological_torsion, _fingerprint_topological_torsion_params, _fingerprint_topological_torsion_keywords, fingerprint_topological_torsion_with_params,
crate::sequence::Output<ck::FingerprintBatch>, crate::sequence::Output, {};
}

crate::collection_params::LayeredFingerprintParams, py, crate::layered_producers::keyword_params(layers,min_path,max_path,fp_size,atom_counts,set_only_bits.as_deref(),branched_paths,from_atoms),
{layers: Option<u32> = None,
min_path: Option<u32> = None,
max_path: Option<u32> = None,
fp_size: Option<u32> = None,
atom_counts: Option<Vec<u32>> = None,
set_only_bits: Option<PyRef<'_, crate::native::Fingerprint>> = None,
branched_paths: Option<bool> = None,
from_atoms: Option<Vec<u32>> = None}
{
fingerprint_layered, _fingerprint_layered_params, _fingerprint_layered_keywords, fingerprint_layered_with_params,
crate::sequence::Output<ck::FingerprintBatch>, crate::sequence::Output, {};
fingerprint_layered_with_output, _fingerprint_layered_with_output_params, _fingerprint_layered_with_output_keywords, fingerprint_layered_with_output_with_params,
Vec<Option<crate::ordered_results::LayeredFingerprintResult>>, |rows| rows.into_iter().map(|row| row.map(|inner| crate::ordered_results::LayeredFingerprintResult { inner })).collect(), {};
}

crate::params::PatternFingerprintParams, py, crate::pattern_producers::keyword_params(n_bits,tautomeric).to_rust(py)?,
{n_bits: Option<usize> = None,
tautomeric: Option<bool> = None}
{
fingerprint_pattern, _fingerprint_pattern_params, _fingerprint_pattern_keywords, fingerprint_pattern_with_params,
crate::sequence::Output<ck::FingerprintBatch>, crate::sequence::Output, {};
}]
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
        let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

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
                params.as_ref().unwrap_or(&resolved),
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
        let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

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
                params.as_ref().unwrap_or(&resolved),
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
            execution: execution
                .map(|v| v.to_rust(py))
                .transpose()?
                .unwrap_or(d.execution),
            filenames,
            report_path,
        };
        let options = options.map(|value| value.to_rust(py)).transpose()?;
        py.detach(|| {
            self.inner
                .write_images_with_params(&directory, options.as_ref().unwrap_or(&resolved))
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

    #[pyo3(signature=(options,atom_invariants,bond_invariants,call,*,execution=None))]
    fn fingerprint_morgan_with_generator_params(
        &self,
        py: Python<'_>,
        options: &crate::fingerprint_configuration::MorganGeneratorParams,
        atom_invariants: Option<&crate::morgan_configuration::MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&crate::fingerprint_providers::MorganBondInvariantsGenerator>,
        call: &crate::params::MorganCallParams,
        execution: Option<&crate::batch_query_params::BatchQueryParams>,
    ) -> PyResult<crate::sequence::Output<ck::FingerprintBatch>> {
        let bond_invariants =
            pyo3::Python::attach(|py| bond_invariants.map(|value| value.to_rust(py)).transpose())?;

        let call = call.to_rust(py)?;
        let options = options.to_rust(py)?;
        crate::batch_calls::query(py, execution, None, None, |runtime| {
                self.inner.fingerprint_morgan_with_generator_params(
                    &options,
                    atom_invariants.map(|v| &v.inner),
                    bond_invariants.as_ref(),
                    &call,
                    runtime,
                )
            })
            .map(crate::sequence::Output)
    }

    #[pyo3(signature=(options,atom_invariants,bond_invariants,call,collect_additional_output,*,execution=None))]
    fn fingerprint_morgan_with_output_with_generator_params(
        &self,
        py: Python<'_>,
        options: &crate::fingerprint_configuration::MorganGeneratorParams,
        atom_invariants: Option<&crate::morgan_configuration::MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&crate::fingerprint_providers::MorganBondInvariantsGenerator>,
        call: &crate::params::MorganCallParams,
        collect_additional_output: bool,
        execution: Option<&crate::batch_query_params::BatchQueryParams>,
    ) -> PyResult<Vec<Option<crate::batch_fingerprint_values::BatchFingerprintOutput>>> {
        let bond_invariants =
            pyo3::Python::attach(|py| bond_invariants.map(|value| value.to_rust(py)).transpose())?;

        let call = call.to_rust(py)?;
        let options = options.to_rust(py)?;
        crate::batch_calls::query(py, execution, None, None, |runtime| {
                self.inner
                    .fingerprint_morgan_with_output_with_generator_params(
                        &options,
                        atom_invariants.map(|v| &v.inner),
                        bond_invariants.as_ref(),
                        &call,
                        collect_additional_output,
                        runtime,
                    )
            })
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
    }

    #[pyo3(signature=(options=Argument::Missing,params=None,*,n_jobs=None,progress_bar=None,isomeric_smiles=Argument::Missing,kekule=Argument::Missing,canonical=Argument::Missing,clean_stereo=Argument::Missing,rooted_at_atom=Argument::Missing,all_bonds_explicit=Argument::Missing,all_hydrogens_explicit=Argument::Missing,include_dative_bonds=Argument::Missing,ignore_atom_map_numbers=Argument::Missing))]
    fn to_smiles(
        &self, py: Python<'_>,
        options: Argument<Option<PyRef<'_, crate::chemical_params::SmilesWriteParams>>>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
        isomeric_smiles: Argument<Option<bool>>,
        kekule: Argument<Option<bool>>,
        canonical: Argument<Option<bool>>,
        clean_stereo: Argument<Option<bool>>,
        rooted_at_atom: Argument<Option<usize>>,
        all_bonds_explicit: Argument<Option<bool>>,
        all_hydrogens_explicit: Argument<Option<bool>>,
        include_dative_bonds: Argument<Option<bool>>,
        ignore_atom_map_numbers: Argument<Option<bool>>,
    ) -> PyResult<Vec<Option<String>>> {
        crate::arguments::configuration_arguments!(options; isomeric_smiles, kekule, canonical, clean_stereo, rooted_at_atom, all_bonds_explicit, all_hydrogens_explicit, include_dative_bonds, ignore_atom_map_numbers);
        if let Some(options) = options.as_deref() {
            self._to_smiles_params(py, Some(options), params, n_jobs, progress_bar)
        } else {
            self._to_smiles_keywords(py, params, n_jobs, progress_bar, isomeric_smiles, kekule, canonical, clean_stereo, rooted_at_atom, all_bonds_explicit, all_hydrogens_explicit, include_dative_bonds, ignore_atom_map_numbers)
        }
    }
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.to_smiles"}"#]
    #[pyo3(signature=(options=None,params=None,*,n_jobs=None,progress_bar=None))]
    fn _to_smiles_params(
        &self, py: Python<'_>,
        options: Option<&crate::chemical_params::SmilesWriteParams>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<String>>> {
        let options = options.map(|value| value.to_rust(py)).transpose()?.unwrap_or_default();
        let values = crate::batch_calls::query(py, params, n_jobs, progress_bar, |runtime| {
            self.inner.to_smiles_with_params(&options, runtime)
        })?;
        values.into_iter().map(|value| value.map(|text| crate::property_values::decode_source_text(py, &text)).transpose()).collect()
    }
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.to_smiles"}"#]
    #[pyo3(signature=(*,params=None,n_jobs=None,progress_bar=None,isomeric_smiles=None,kekule=None,canonical=None,clean_stereo=None,rooted_at_atom=None,all_bonds_explicit=None,all_hydrogens_explicit=None,include_dative_bonds=None,ignore_atom_map_numbers=None))]
    fn _to_smiles_keywords(
        &self, py: Python<'_>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
        isomeric_smiles: Option<bool>,
        kekule: Option<bool>,
        canonical: Option<bool>,
        clean_stereo: Option<bool>,
        rooted_at_atom: Option<usize>,
        all_bonds_explicit: Option<bool>,
        all_hydrogens_explicit: Option<bool>,
        include_dative_bonds: Option<bool>,
        ignore_atom_map_numbers: Option<bool>,
    ) -> PyResult<Vec<Option<String>>> {
        let options = crate::smiles_writer::keywords(isomeric_smiles, kekule, canonical, clean_stereo, rooted_at_atom, all_bonds_explicit, all_hydrogens_explicit, include_dative_bonds, ignore_atom_map_numbers);
        self._to_smiles_params(py, Some(&options), params, n_jobs, progress_bar)
    }
    fn to_smiles_with_params(
        &self, py: Python<'_>,
        options: &crate::chemical_params::SmilesWriteParams,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<String>>> {
        self._to_smiles_params(py, Some(options), Some(params), None, None)
    }
    #[pyo3(signature=(params=None,*,n_jobs=None,progress_bar=None))]
    fn dg_bounds_matrix(
        &self,
        py: Python<'_>,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<Vec<Vec<f64>>>>> {
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |p| {
            self.inner.dg_bounds_matrix_with_params(p)
        })
    }
    fn dg_bounds_matrix_with_params(
        &self,
        py: Python<'_>,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<Vec<Vec<f64>>>>> {
        self.dg_bounds_matrix(py, Some(params), None, None)
    }
    #[pyo3(signature=(width,height,params=None,*,n_jobs=None,progress_bar=None))]
    fn to_svg(
        &self,
        py: Python<'_>,
        width: u32,
        height: u32,
        params: Option<&crate::batch_query_params::BatchQueryParams>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Vec<Option<String>>> {
        crate::batch_calls::query(py, params, n_jobs, progress_bar, |p| {
            self.inner.to_svg_with_params(width, height, p)
        })
    }
    fn to_svg_with_params(
        &self,
        py: Python<'_>,
        width: u32,
        height: u32,
        params: &crate::batch_query_params::BatchQueryParams,
    ) -> PyResult<Vec<Option<String>>> {
        self.to_svg(py, width, height, Some(params), None, None)
    }

    #[staticmethod]
    #[pyo3(signature=(smiles,parse=Argument::Missing,params=None,*,errors=None,n_jobs=None,progress_bar=None,sanitize=Argument::Missing,allow_cxsmiles=Argument::Missing,strict_cxsmiles=Argument::Missing,parse_name=Argument::Missing,remove_hs=Argument::Missing,skip_cleanup=Argument::Missing,debug_parse=Argument::Missing,replacements=Argument::Missing))]
    pub(crate) fn from_smiles_list(
        py: Python<'_>, smiles: Vec<String>,
        parse: Argument<Option<PyRef<'_, SmilesParseParams>>>,
        params: Option<&crate::io_batch_params::BatchParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
        sanitize: Argument<Option<bool>>,
        allow_cxsmiles: Argument<Option<bool>>,
        strict_cxsmiles: Argument<Option<bool>>,
        parse_name: Argument<Option<bool>>,
        remove_hs: Argument<Option<bool>>,
        skip_cleanup: Argument<Option<bool>>,
        debug_parse: Argument<Option<bool>>,
        replacements: Argument<Option<BTreeMap<String, String>>>,
    ) -> PyResult<Self> {
        crate::arguments::configuration_arguments!(parse; sanitize, allow_cxsmiles, strict_cxsmiles, parse_name, remove_hs, skip_cleanup, debug_parse, replacements);
        if let Some(parse) = parse.as_deref() {
            Self::_from_smiles_list_params(py, smiles, Some(parse), params, errors, n_jobs, progress_bar)
        } else {
            Self::_from_smiles_list_keywords(py, smiles, params, errors, n_jobs, progress_bar, sanitize, allow_cxsmiles, strict_cxsmiles, parse_name, remove_hs, skip_cleanup, debug_parse, replacements)
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.from_smiles_list"}"#]
    #[pyo3(signature=(smiles,parse=None,params=None,*,errors=None,n_jobs=None,progress_bar=None))]
    pub(crate) fn _from_smiles_list_params(
        py: Python<'_>, smiles: Vec<String>, parse: Option<&SmilesParseParams>,
        params: Option<&crate::io_batch_params::BatchParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let parse = parse.map(|value| value.to_rust(py)).transpose()?.unwrap_or_default();
        batch_from_smiles(py, smiles, &parse, params, errors, n_jobs, progress_bar)
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.from_smiles_list"}"#]
    #[pyo3(signature=(smiles,*,params=None,errors=None,n_jobs=None,progress_bar=None,sanitize=None,allow_cxsmiles=None,strict_cxsmiles=None,parse_name=None,remove_hs=None,skip_cleanup=None,debug_parse=None,replacements=None))]
    pub(crate) fn _from_smiles_list_keywords(
        py: Python<'_>, smiles: Vec<String>,
        params: Option<&crate::io_batch_params::BatchParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
        sanitize: Option<bool>,
        allow_cxsmiles: Option<bool>,
        strict_cxsmiles: Option<bool>,
        parse_name: Option<bool>,
        remove_hs: Option<bool>,
        skip_cleanup: Option<bool>,
        debug_parse: Option<bool>,
        replacements: Option<BTreeMap<String, String>>,
    ) -> PyResult<Self> {
        let parse = parse_keywords(sanitize, allow_cxsmiles, strict_cxsmiles, parse_name, remove_hs, skip_cleanup, debug_parse, replacements);
        batch_from_smiles(py, smiles, &parse, params, errors, n_jobs, progress_bar)
    }
    #[staticmethod]
    fn from_smiles_list_with_params(
        py: Python<'_>, smiles: Vec<String>, parse: &SmilesParseParams,
        params: &crate::io_batch_params::BatchParams,
    ) -> PyResult<Self> {
        Self::_from_smiles_list_params(py, smiles, Some(parse), Some(params), None, None, None)
    }
    #[staticmethod]
    #[pyo3(signature=(text,read=crate::arguments::Argument::Missing,*,errors=None,n_jobs=None,progress_bar=None,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
    pub(crate) fn from_sdf_records(
        py: Python<'_>, text: String,
        read: crate::arguments::Argument<Option<PyRef<'_, crate::selectors::SdfReadParams>>>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
        sanitize: crate::arguments::Argument<Option<bool>>,
        remove_hs: crate::arguments::Argument<Option<bool>>,
        strict_parsing: crate::arguments::Argument<Option<bool>>,
        expand_attachment_points: crate::arguments::Argument<Option<bool>>,
        process_property_lists: crate::arguments::Argument<Option<bool>>,
        coordinate_mode: crate::arguments::Argument<Option<crate::selectors::SdfCoordinateModeInput>>,
    ) -> PyResult<Self> {
        let keyword_form = sanitize.is_supplied() || remove_hs.is_supplied() || strict_parsing.is_supplied() || expand_attachment_points.is_supplied() || process_property_lists.is_supplied() || coordinate_mode.is_supplied();
        crate::arguments::configuration_arguments!(read; sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode);
        if keyword_form {
            Self::_from_sdf_records_keywords(py, text, errors, n_jobs, progress_bar, sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode)
        } else {
            Self::_from_sdf_records_params(py, text, read.as_deref(), errors, n_jobs, progress_bar)
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.from_sdf_records"}"#]
    #[pyo3(signature=(text,read=None,*,errors=None,n_jobs=None,progress_bar=None))]
    pub(crate) fn _from_sdf_records_params(
        py: Python<'_>, text: String, read: Option<&crate::selectors::SdfReadParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let read = read.map(|value| value.to_rust(py)).transpose()?.unwrap_or_default();
        let errors = errors.map(crate::io_batch_params::BatchErrorModeInput::into_rust).transpose()?.unwrap_or_default();
        let input = text.as_str();
        py.detach(|| ck::MoleculeBatch::from_sdf_records_with_params(input, &read, errors, n_jobs, progress_bar.unwrap_or(false)))
            .map(|inner| Self { inner })
            .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.from_sdf_records"}"#]
    #[pyo3(signature=(text,*,errors=None,n_jobs=None,progress_bar=None,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
    pub(crate) fn _from_sdf_records_keywords(
        py: Python<'_>, text: String,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
        sanitize: Option<bool>,
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
        expand_attachment_points: Option<bool>,
        process_property_lists: Option<bool>,
        coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
    ) -> PyResult<Self> {
        let read = crate::molecular_io_calls::sdf_keywords(sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode)?;
        Self::_from_sdf_records_params(py, text, Some(&read), errors, n_jobs, progress_bar)
    }
    #[staticmethod]
    #[pyo3(signature=(text,read,errors=None,n_jobs=None,progress_bar=None))]
    fn from_sdf_records_with_params(
        py: Python<'_>, text: String, read: &crate::selectors::SdfReadParams,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        Self::_from_sdf_records_params(py, text, Some(read), errors, n_jobs, progress_bar)
    }

    #[staticmethod]
    #[pyo3(signature=(path,read=crate::arguments::Argument::Missing,*,errors=None,n_jobs=None,progress_bar=None,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
    pub(crate) fn read_sdf(
        py: Python<'_>, path: std::path::PathBuf,
        read: crate::arguments::Argument<Option<PyRef<'_, crate::selectors::SdfReadParams>>>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
        sanitize: crate::arguments::Argument<Option<bool>>,
        remove_hs: crate::arguments::Argument<Option<bool>>,
        strict_parsing: crate::arguments::Argument<Option<bool>>,
        expand_attachment_points: crate::arguments::Argument<Option<bool>>,
        process_property_lists: crate::arguments::Argument<Option<bool>>,
        coordinate_mode: crate::arguments::Argument<Option<crate::selectors::SdfCoordinateModeInput>>,
    ) -> PyResult<Self> {
        let keyword_form = sanitize.is_supplied() || remove_hs.is_supplied() || strict_parsing.is_supplied() || expand_attachment_points.is_supplied() || process_property_lists.is_supplied() || coordinate_mode.is_supplied();
        crate::arguments::configuration_arguments!(read; sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode);
        if keyword_form {
            Self::_read_sdf_keywords(py, path, errors, n_jobs, progress_bar, sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode)
        } else {
            Self::_read_sdf_params(py, path, read.as_deref(), errors, n_jobs, progress_bar)
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.read_sdf"}"#]
    #[pyo3(signature=(path,read=None,*,errors=None,n_jobs=None,progress_bar=None))]
    pub(crate) fn _read_sdf_params(
        py: Python<'_>, path: std::path::PathBuf, read: Option<&crate::selectors::SdfReadParams>,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        let read = read.map(|value| value.to_rust(py)).transpose()?.unwrap_or_default();
        let errors = errors.map(crate::io_batch_params::BatchErrorModeInput::into_rust).transpose()?.unwrap_or_default();
        let input = crate::molecular_io_calls::path_text(&path)?;
        py.detach(|| ck::MoleculeBatch::read_sdf_with_params(input, &read, errors, n_jobs, progress_bar.unwrap_or(false)))
            .map(|inner| Self { inner })
            .map_err(crate::batch_errors::BatchValidationError::from_rust)
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"MoleculeBatch.read_sdf"}"#]
    #[pyo3(signature=(path,*,errors=None,n_jobs=None,progress_bar=None,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
    pub(crate) fn _read_sdf_keywords(
        py: Python<'_>, path: std::path::PathBuf,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
        sanitize: Option<bool>,
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
        expand_attachment_points: Option<bool>,
        process_property_lists: Option<bool>,
        coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
    ) -> PyResult<Self> {
        let read = crate::molecular_io_calls::sdf_keywords(sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode)?;
        Self::_read_sdf_params(py, path, Some(&read), errors, n_jobs, progress_bar)
    }
    #[staticmethod]
    #[pyo3(signature=(path,read,errors=None,n_jobs=None,progress_bar=None))]
    fn read_sdf_with_params(
        py: Python<'_>, path: std::path::PathBuf, read: &crate::selectors::SdfReadParams,
        errors: Option<crate::io_batch_params::BatchErrorModeInput>,
        n_jobs: Option<usize>, progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        Self::_read_sdf_params(py, path, Some(read), errors, n_jobs, progress_bar)
    }

    #[pyo3(signature=(options=None,params=crate::arguments::Argument::Missing,*,errors=crate::arguments::Argument::Missing,n_jobs=crate::arguments::Argument::Missing,progress_bar=crate::arguments::Argument::Missing))]
    fn sanitize(
        &self,
        py: Python<'_>,
        options: Option<&crate::sanitize_values::SanitizeParams>,
        params: crate::arguments::Argument<Option<PyRef<'_, crate::io_batch_params::BatchParams>>>,
        errors: crate::arguments::Argument<Option<crate::io_batch_params::BatchErrorModeInput>>,
        n_jobs: crate::arguments::Argument<Option<usize>>,
        progress_bar: crate::arguments::Argument<Option<bool>>,
    ) -> PyResult<Self> {
        let keyword_form =
            errors.is_supplied() || n_jobs.is_supplied() || progress_bar.is_supplied();
        crate::arguments::configuration_arguments!(params; errors, n_jobs, progress_bar);
        let params = params.as_deref();
        if keyword_form {
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
        let options =
            pyo3::Python::attach(|py| options.map(|value| value.to_rust(py)).transpose())?;
        let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

        let defaults = ck::SanitizeParams::default();
        let execution = ck::BatchParams::default();
        py.detach(|| {
            if options.is_none() && params.is_none() {
                self.inner.sanitize()
            } else {
                self.inner.sanitize_with_params(
                    options.as_ref().unwrap_or(&defaults),
                    params.as_ref().unwrap_or(&execution),
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
    #[pyo3(signature=(options=None,params=crate::arguments::Argument::Missing,*,errors=crate::arguments::Argument::Missing,n_jobs=crate::arguments::Argument::Missing,progress_bar=crate::arguments::Argument::Missing))]
    fn with_hydrogens(
        &self,
        py: Python<'_>,
        options: Option<&crate::chemical_params::AddHsParams>,
        params: crate::arguments::Argument<Option<PyRef<'_, crate::io_batch_params::BatchParams>>>,
        errors: crate::arguments::Argument<Option<crate::io_batch_params::BatchErrorModeInput>>,
        n_jobs: crate::arguments::Argument<Option<usize>>,
        progress_bar: crate::arguments::Argument<Option<bool>>,
    ) -> PyResult<Self> {
        let keyword_form =
            errors.is_supplied() || n_jobs.is_supplied() || progress_bar.is_supplied();
        crate::arguments::configuration_arguments!(params; errors, n_jobs, progress_bar);
        let params = params.as_deref();
        if keyword_form {
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
        let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

        let defaults = ck::AddHsParams::default();
        let execution = ck::BatchParams::default();
        let options = options.map(|value| value.to_rust(py)).transpose()?;
        py.detach(|| {
            if options.is_none() && params.is_none() {
                self.inner.with_hydrogens()
            } else {
                self.inner.with_hydrogens_with_params(
                    options.as_ref().unwrap_or(&defaults),
                    params.as_ref().unwrap_or(&execution),
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
    #[pyo3(signature=(options=None,params=crate::arguments::Argument::Missing,*,errors=crate::arguments::Argument::Missing,n_jobs=crate::arguments::Argument::Missing,progress_bar=crate::arguments::Argument::Missing))]
    fn without_hydrogens(
        &self,
        py: Python<'_>,
        options: Option<&crate::params::RemoveHsParams>,
        params: crate::arguments::Argument<Option<PyRef<'_, crate::io_batch_params::BatchParams>>>,
        errors: crate::arguments::Argument<Option<crate::io_batch_params::BatchErrorModeInput>>,
        n_jobs: crate::arguments::Argument<Option<usize>>,
        progress_bar: crate::arguments::Argument<Option<bool>>,
    ) -> PyResult<Self> {
        let keyword_form =
            errors.is_supplied() || n_jobs.is_supplied() || progress_bar.is_supplied();
        crate::arguments::configuration_arguments!(params; errors, n_jobs, progress_bar);
        let params = params.as_deref();
        if keyword_form {
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
        let options =
            pyo3::Python::attach(|py| options.map(|value| value.to_rust(py)).transpose())?;
        let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

        let defaults = ck::RemoveHsParams::default();
        let execution = ck::BatchParams::default();
        py.detach(|| {
            if options.is_none() && params.is_none() {
                self.inner.without_hydrogens()
            } else {
                self.inner.without_hydrogens_with_params(
                    options.as_ref().unwrap_or(&defaults),
                    params.as_ref().unwrap_or(&execution),
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
    #[pyo3(signature=(options=None,params=crate::arguments::Argument::Missing,*,errors=crate::arguments::Argument::Missing,n_jobs=crate::arguments::Argument::Missing,progress_bar=crate::arguments::Argument::Missing))]
    fn with_kekulized_bonds(
        &self,
        py: Python<'_>,
        options: Option<&crate::params::KekulizeParams>,
        params: crate::arguments::Argument<Option<PyRef<'_, crate::io_batch_params::BatchParams>>>,
        errors: crate::arguments::Argument<Option<crate::io_batch_params::BatchErrorModeInput>>,
        n_jobs: crate::arguments::Argument<Option<usize>>,
        progress_bar: crate::arguments::Argument<Option<bool>>,
    ) -> PyResult<Self> {
        let keyword_form =
            errors.is_supplied() || n_jobs.is_supplied() || progress_bar.is_supplied();
        crate::arguments::configuration_arguments!(params; errors, n_jobs, progress_bar);
        let params = params.as_deref();
        if keyword_form {
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
        let options =
            pyo3::Python::attach(|py| options.map(|value| value.to_rust(py)).transpose())?;
        let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

        let defaults = ck::KekulizeParams::default();
        let execution = ck::BatchParams::default();
        py.detach(|| {
            if options.is_none() && params.is_none() {
                self.inner.with_kekulized_bonds()
            } else {
                self.inner.with_kekulized_bonds_with_params(
                    options.as_ref().unwrap_or(&defaults),
                    params.as_ref().unwrap_or(&execution),
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
    #[pyo3(signature=(options=None,params=crate::arguments::Argument::Missing,*,errors=crate::arguments::Argument::Missing,n_jobs=crate::arguments::Argument::Missing,progress_bar=crate::arguments::Argument::Missing))]
    fn with_2d_coordinates(
        &self,
        py: Python<'_>,
        options: Option<&crate::optional_configuration::Coordinate2DParams>,
        params: crate::arguments::Argument<Option<PyRef<'_, crate::io_batch_params::BatchParams>>>,
        errors: crate::arguments::Argument<Option<crate::io_batch_params::BatchErrorModeInput>>,
        n_jobs: crate::arguments::Argument<Option<usize>>,
        progress_bar: crate::arguments::Argument<Option<bool>>,
    ) -> PyResult<Self> {
        let keyword_form =
            errors.is_supplied() || n_jobs.is_supplied() || progress_bar.is_supplied();
        crate::arguments::configuration_arguments!(params; errors, n_jobs, progress_bar);
        let params = params.as_deref();
        if keyword_form {
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
        let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

        let defaults = ck::Coordinate2DParams::default();
        let execution = ck::BatchParams::default();
        let options = options.map(|value| value.to_rust(py)).transpose()?;
        py.detach(|| {
            if options.is_none() && params.is_none() {
                self.inner.with_2d_coordinates()
            } else {
                self.inner.with_2d_coordinates_with_params(
                    options.as_ref().unwrap_or(&defaults),
                    params.as_ref().unwrap_or(&execution),
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

    #[getter]
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
    fn select(&self, indices: Vec<usize>) -> Option<Self> {
        self.inner.select(&indices).map(|inner| Self { inner })
    }
    #[doc = r#"{"role":"indexing","protocol":"integer_selection","selection":"Self"}"#]
    fn __getitem__(&self, py: Python<'_>, index: crate::indexing::Index<'_>)
        -> PyResult<crate::indexing::Indexed<Option<Molecule>, Self>> {
        crate::indexing::get(self, py, index)
    }
    #[getter]
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
    #[getter]
    fn len(&self) -> usize {
        self.inner.len()
    }
    #[getter]
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    #[getter]
    fn error_mode(&self) -> BatchErrorMode {
        self.inner.error_mode().into()
    }
    #[getter]
    fn valid_mask(&self) -> Vec<bool> {
        self.inner.valid_mask()
    }
    #[getter]
    fn invalid_mask(&self) -> Vec<bool> {
        self.inner.invalid_mask()
    }
    #[getter]
    fn valid_count(&self) -> usize {
        self.inner.valid_count()
    }
    #[getter]
    fn invalid_count(&self) -> usize {
        self.inner.invalid_count()
    }
    #[getter]
    fn parallel_jobs(&self) -> Option<usize> {
        self.inner.parallel_jobs()
    }
    #[getter]
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
