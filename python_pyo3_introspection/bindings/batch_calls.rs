//! Native batch call conversion helpers; no scheduling or chemistry is reimplemented.
#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::{
        ffi::c_str,
        types::{PyDict, PyModule},
    };
    #[test]
    fn complete_batch_calls_preserve_native_order_errors_cow_and_runtime_callbacks() {
        Python::initialize();
        Python::attach(|py|->PyResult<()>{
let native=ck::MoleculeBatch::from_smiles_list_with_params(&["CO".into(),"C(".into(),"N".into()],&Default::default(),&ck::BatchParams{errors:Some(ck::BatchErrorMode::KeepErrors),n_jobs:Some(1),progress_bar:Some(false)}).unwrap();
let batch=Py::new(py,MoleculeBatch{inner:native.clone()})?;
let values=batch.bind(py).call_method0("records")?.extract::<Vec<Record>>()?;assert_eq!(values.len(),3);
let Source::First(first)=&values[0]else{panic!("valid record lost")};
let ck::BatchRecord::Molecule(original)=native.get(0).unwrap()else{panic!("source")};
assert!(std::ptr::eq(first.borrow(py).inner.topology(),original.topology()));
let actual_blocks=first.borrow(py).inner.to_builder();let original_blocks=original.to_builder();
assert!(std::ptr::eq(actual_blocks.coordinates(),original_blocks.coordinates()));
assert!(std::ptr::eq(first.borrow(py).inner.properties(),original.properties()));
let errors=batch.bind(py).call_method0("errors")?.extract::<Vec<Py<BatchError>>>()?;assert_eq!(errors[0].bind(py).getattr("index")?.extract::<usize>()?,1);
let cause=errors[0].bind(py).getattr("source")?;assert!(cause.is_instance_of::<crate::operation_errors::SmilesError>());assert!(cause.is(&errors[0].bind(py).getattr("__cause__")?));
let renamed=batch.bind(py).call_method1("with_parallel_jobs",(Some(1),))?.extract::<Py<MoleculeBatch>>()?;
assert_eq!(renamed.borrow(py).inner.len(),native.len());let ck::BatchRecord::Molecule(copy)=renamed.borrow(py).inner.get(0).unwrap().clone()else{panic!("valid")};assert!(std::ptr::eq(copy.topology(),original.topology()));assert!(std::ptr::eq(copy.properties(),original.properties()));
let invalid=batch.bind(py).call_method1("with_parallel_jobs",(Some(0),)).unwrap_err();assert!(invalid.is_instance_of::<BatchValidationError>(py));
let transformed=batch.bind(py).call_method0("with_hydrogens")?.extract::<Py<MoleculeBatch>>()?;assert_eq!(transformed.borrow(py).inner.to_list(),native.with_hydrogens().unwrap().to_list());
macro_rules! transform{($name:ident,$alias:ident,$module:ident,$params:ident)=>{{
let options=ck::$params::default();let params=py.get_type::<crate::$module::$params>().call0()?;let runtime=Py::new(py,BatchParams{inner:crate::parameters::NativeField::Owned(ck::BatchParams{errors:Some(ck::BatchErrorMode::KeepErrors),n_jobs:Some(1),progress_bar:Some(false)})})?;
let kw=PyDict::new(py);kw.set_item("n_jobs",1)?;kw.set_item("errors",crate::io_batch_params::BatchErrorMode::KeepErrors)?;kw.set_item("progress_bar",false)?;
let want=native.$alias(&options,&runtime.borrow(py).to_rust(py)?).unwrap();
for actual in [batch.bind(py).call_method0(stringify!($name))?,batch.bind(py).call_method(stringify!($name),(&params,),Some(&kw))?,batch.bind(py).call_method1(stringify!($alias),(&params,&runtime))?]{assert_eq!(actual.extract::<Py<MoleculeBatch>>()?.borrow(py).inner.to_list(),want.to_list());}
assert!(batch.bind(py).call_method(stringify!($name),(&params,&runtime),Some(&kw)).unwrap_err().is_instance_of::<pyo3::exceptions::PyTypeError>(py));
}}}
transform!(sanitize,sanitize_with_params,sanitize_values,SanitizeParams);
transform!(with_hydrogens,with_hydrogens_with_params,chemical_params,AddHsParams);
transform!(without_hydrogens,without_hydrogens_with_params,params,RemoveHsParams);
transform!(with_kekulized_bonds,with_kekulized_bonds_with_params,params,KekulizeParams);
transform!(with_2d_coordinates,with_2d_coordinates_with_params,optional_configuration,Coordinate2DParams);
macro_rules! bits{($name:ident,$alias:ident,$module:ident,$params:ident)=>{{
let defaults=ck::$params::default();let options=py.get_type::<crate::$module::$params>().call0()?.extract::<Py<crate::$module::$params>>()?;
let runtime=py.get_type::<crate::batch_query_params::BatchQueryParams>().call0()?;
let want=native.$alias(&defaults,&Default::default()).unwrap();
for actual in [batch.bind(py).call_method0(stringify!($name))?,batch.bind(py).call_method1(stringify!($alias),(&options,&runtime))?]{let got=actual.extract::<Vec<Option<Py<crate::native::Fingerprint>>>>()?;assert_eq!(got.len(),want.len());for(got,want)in got.iter().zip(&want){assert_eq!(got.as_ref().map(|v|v.borrow(py).inner.on_bits()),want.as_ref().map(ck::Fingerprint::on_bits));}}
}}}
bits!(fingerprint_morgan,fingerprint_morgan_with_params,morgan_configuration,MorganFingerprintParams);
bits!(fingerprint_atom_pair,fingerprint_atom_pair_with_params,fingerprint_configuration,AtomPairFingerprintParams);
bits!(fingerprint_topological_torsion,fingerprint_topological_torsion_with_params,fingerprint_configuration,TopologicalTorsionFingerprintParams);
bits!(fingerprint_layered,fingerprint_layered_with_params,collection_params,LayeredFingerprintParams);
bits!(fingerprint_pattern,fingerprint_pattern_with_params,params,PatternFingerprintParams);
for name in ["fingerprint_atom_pair_sparse_count","fingerprint_atom_pair_count","fingerprint_atom_pair_sparse","fingerprint_layered_with_output","fingerprint_morgan_with_output","fingerprint_atom_pair_with_output","dg_bounds_matrix"]{let result=batch.bind(py).call_method0(name)?;assert_eq!(result.len()?,3);assert!(result.get_item(1)?.is_none());}
let callbacks=PyModule::from_code(py,c_str!("seen=[]\nmarker=LookupError('original batch progress')\ndef notify(): seen.append(1)\ndef fail(): raise marker\n"),c_str!("batch-progress.py"),c_str!("batch_progress"))?;
let kw=PyDict::new(py);kw.set_item("n_jobs",1)?;kw.set_item("progress_callback",callbacks.getattr("notify")?)?;
let runtime=py.get_type::<crate::batch_query_params::BatchQueryParams>().call((),Some(&kw))?;
batch.bind(py).call_method1("to_smiles",(py.None(),&runtime))?;assert_eq!(callbacks.getattr("seen")?.len()?,3);
runtime.setattr("progress_callback",callbacks.getattr("fail")?)?;
let fail=batch.bind(py).call_method1("to_smiles",(py.None(),&runtime)).unwrap_err();assert!(fail.value(py).is(&callbacks.getattr("marker")?));
callbacks.add("BatchQueryParams",py.get_type::<crate::batch_query_params::BatchQueryParams>())?;
py.run(c_str!("import gc, weakref\nclass Cycle:\n    def __init__(self, params): self.params=params\n    def __call__(self): pass\np=BatchQueryParams()\nc=Cycle(p)\nobserver=weakref.ref(c)\np.progress_callback=c\ndel c,p\ngc.collect()\nassert observer() is None\n"),Some(&callbacks.dict()),None)?;
let directory=std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("artifacts").join(format!("batch-complete-{}",std::process::id()));std::fs::create_dir_all(&directory).unwrap();
let report=batch.bind(py).call_method1("write_sdf",(directory.join("out.sdf"),))?.extract::<Py<crate::batch_errors::BatchExportReport>>()?;assert_eq!(report.bind(py).call_method0("success")?.extract::<usize>()?,2);assert_eq!(report.bind(py).call_method0("failed")?.extract::<usize>()?,1);
report.bind(py).call_method1("write_report",(directory.join("report.json"),))?;assert!(std::fs::read_to_string(directory.join("report.json")).unwrap().contains("written"));
let kw=PyDict::new(py);kw.set_item("format","svg")?;kw.set_item("width",120)?;kw.set_item("height",100)?;
let images=batch.bind(py).call_method("write_images",(directory.join("images"),),Some(&kw))?;assert_eq!(images.call_method0("success")?.extract::<usize>()?,2);
let svg=batch.bind(py).call_method1("to_svg",(120,100))?;assert_eq!(svg.len()?,3);assert!(svg.get_item(1)?.is_none());
let keep=batch.bind(py).call_method0("errors")?.get_item(0)?;drop(batch);drop(native);assert_eq!(keep.getattr("index")?.extract::<usize>()?,1);assert!(keep.getattr("source")?.is_instance_of::<crate::operation_errors::SmilesError>());
Ok(())}).unwrap();
    }
}

#[pyfunction]
#[pyo3(signature=(smiles,parse=crate::arguments::Argument::Missing,params=None,*,errors=None,n_jobs=None,progress_bar=None,sanitize=crate::arguments::Argument::Missing,allow_cxsmiles=crate::arguments::Argument::Missing,strict_cxsmiles=crate::arguments::Argument::Missing,parse_name=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,skip_cleanup=crate::arguments::Argument::Missing,debug_parse=crate::arguments::Argument::Missing,replacements=crate::arguments::Argument::Missing))]
pub fn mols_from_smiles_list(
    py: Python<'_>,
    smiles: Vec<String>,
    parse: crate::arguments::Argument<Option<PyRef<'_, crate::extension::SmilesParseParams>>>,
    params: Option<&BatchParams>,
    errors: Option<BatchErrorModeInput>,
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
    sanitize: crate::arguments::Argument<Option<bool>>,
    allow_cxsmiles: crate::arguments::Argument<Option<bool>>,
    strict_cxsmiles: crate::arguments::Argument<Option<bool>>,
    parse_name: crate::arguments::Argument<Option<bool>>,
    remove_hs: crate::arguments::Argument<Option<bool>>,
    skip_cleanup: crate::arguments::Argument<Option<bool>>,
    debug_parse: crate::arguments::Argument<Option<bool>>,
    replacements: crate::arguments::Argument<Option<std::collections::BTreeMap<String, String>>>,
) -> PyResult<MoleculeBatch> {
    MoleculeBatch::from_smiles_list(
        py,
        smiles,
        parse,
        params,
        errors,
        n_jobs,
        progress_bar,
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
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mols_from_smiles_list"}"#]
#[pyo3(signature=(smiles,parse=None,params=None,*,errors=None,n_jobs=None,progress_bar=None))]
pub fn _mols_from_smiles_list_params(
    py: Python<'_>,
    smiles: Vec<String>,
    parse: Option<&crate::extension::SmilesParseParams>,
    params: Option<&BatchParams>,
    errors: Option<BatchErrorModeInput>,
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
) -> PyResult<MoleculeBatch> {
    MoleculeBatch::_from_smiles_list_params(py, smiles, parse, params, errors, n_jobs, progress_bar)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mols_from_smiles_list"}"#]
#[pyo3(signature=(smiles,*,params=None,errors=None,n_jobs=None,progress_bar=None,sanitize=None,allow_cxsmiles=None,strict_cxsmiles=None,parse_name=None,remove_hs=None,skip_cleanup=None,debug_parse=None,replacements=None))]
pub fn _mols_from_smiles_list_keywords(
    py: Python<'_>,
    smiles: Vec<String>,
    params: Option<&BatchParams>,
    errors: Option<BatchErrorModeInput>,
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
    sanitize: Option<bool>,
    allow_cxsmiles: Option<bool>,
    strict_cxsmiles: Option<bool>,
    parse_name: Option<bool>,
    remove_hs: Option<bool>,
    skip_cleanup: Option<bool>,
    debug_parse: Option<bool>,
    replacements: Option<std::collections::BTreeMap<String, String>>,
) -> PyResult<MoleculeBatch> {
    MoleculeBatch::_from_smiles_list_keywords(
        py,
        smiles,
        params,
        errors,
        n_jobs,
        progress_bar,
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
#[pyfunction]
#[pyo3(signature=(records,mode=None))]
pub fn mols_from_records(
    py: Python<'_>,
    records: Vec<Record>,
    mode: Option<BatchErrorModeInput>,
) -> PyResult<MoleculeBatch> {
    MoleculeBatch::from_records(py, records, mode)
}
#[pyfunction]
#[pyo3(signature=(text,read=crate::arguments::Argument::Missing,*,errors=None,n_jobs=None,progress_bar=None,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
pub fn mols_from_sdf_records(
    py: Python<'_>,
    text: String,
    read: crate::arguments::Argument<Option<PyRef<'_, crate::selectors::SdfReadParams>>>,
    errors: Option<BatchErrorModeInput>,
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
    sanitize: crate::arguments::Argument<Option<bool>>,
    remove_hs: crate::arguments::Argument<Option<bool>>,
    strict_parsing: crate::arguments::Argument<Option<bool>>,
    expand_attachment_points: crate::arguments::Argument<Option<bool>>,
    process_property_lists: crate::arguments::Argument<Option<bool>>,
    coordinate_mode: crate::arguments::Argument<Option<crate::selectors::SdfCoordinateModeInput>>,
) -> PyResult<MoleculeBatch> {
    MoleculeBatch::from_sdf_records(
        py,
        text,
        read,
        errors,
        n_jobs,
        progress_bar,
        sanitize,
        remove_hs,
        strict_parsing,
        expand_attachment_points,
        process_property_lists,
        coordinate_mode,
    )
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mols_from_sdf_records"}"#]
#[pyo3(signature=(text,read=None,*,errors=None,n_jobs=None,progress_bar=None))]
pub fn _mols_from_sdf_records_params(
    py: Python<'_>,
    text: String,
    read: Option<&crate::selectors::SdfReadParams>,
    errors: Option<BatchErrorModeInput>,
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
) -> PyResult<MoleculeBatch> {
    MoleculeBatch::_from_sdf_records_params(py, text, read, errors, n_jobs, progress_bar)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mols_from_sdf_records"}"#]
#[pyo3(signature=(text,*,errors=None,n_jobs=None,progress_bar=None,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
pub fn _mols_from_sdf_records_keywords(
    py: Python<'_>,
    text: String,
    errors: Option<BatchErrorModeInput>,
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    strict_parsing: Option<bool>,
    expand_attachment_points: Option<bool>,
    process_property_lists: Option<bool>,
    coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
) -> PyResult<MoleculeBatch> {
    MoleculeBatch::_from_sdf_records_keywords(
        py,
        text,
        errors,
        n_jobs,
        progress_bar,
        sanitize,
        remove_hs,
        strict_parsing,
        expand_attachment_points,
        process_property_lists,
        coordinate_mode,
    )
}

use crate::{
    batch_errors::{BatchError, BatchValidationError},
    error_framework::Source,
    extension::MoleculeBatch,
    io_batch_params::{BatchErrorModeInput, BatchParams},
    native::Molecule,
};
use cosmolkit as ck;
use pyo3::prelude::*;
pub(crate) type Record = Source<Py<Molecule>, Py<BatchError>>;
pub(crate) fn record(py: Python<'_>, v: &ck::BatchRecord) -> PyResult<Record> {
    match v {
        ck::BatchRecord::Molecule(inner) => Ok(Source::First(Py::new(
            py,
            Molecule {
                inner: inner.clone(),
            },
        )?)),
        ck::BatchRecord::Error(e) => Ok(Source::Second(
            BatchError::from_rust(e.clone()).value(py).extract()?,
        )),
    }
}
pub(crate) fn native_record(py: Python<'_>, v: Record) -> ck::BatchRecord {
    match v {
        Source::First(v) => ck::BatchRecord::Molecule(v.borrow(py).inner.clone()),
        Source::Second(v) => ck::BatchRecord::Error(v.borrow(py).inner.as_ref().clone()),
    }
}
pub(crate) fn execution(
    errors: Option<BatchErrorModeInput>,
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
) -> PyResult<BatchParams> {
    Ok(BatchParams {
        inner: crate::parameters::NativeField::Owned(ck::BatchParams {
            errors: errors.map(BatchErrorModeInput::into_rust).transpose()?,
            n_jobs,
            progress_bar,
        }),
    })
}
pub(crate) fn query<T: Send>(
    py: Python<'_>,
    params: Option<&crate::batch_query_params::BatchQueryParams>,
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
    run: impl FnOnce(&ck::BatchQueryParams) -> Result<T, ck::BatchValidationError> + Send,
) -> PyResult<T> {
    if params.is_some() && (n_jobs.is_some() || progress_bar.is_some()) {
        return Err(pyo3::exceptions::PyTypeError::new_err(
            "Params and runtime keywords are mutually exclusive",
        ));
    }
    let result = match params {
        Some(p) => p.with_callback(py, run)?,
        None => {
            let p = ck::BatchQueryParams {
                n_jobs,
                progress_bar,
                progress_callback: None,
            };
            py.detach(|| run(&p))
        }
    };
    result.map_err(BatchValidationError::from_rust)
}
