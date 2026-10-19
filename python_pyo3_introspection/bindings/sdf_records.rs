//! Finalized native SDF records: one graph owner, exact fields and reader policy.
#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::ffi::c_str;
    use pyo3::types::{PyDict, PyModule};

    #[test]
    fn finalized_readers_records_and_causes_preserve_native_values_fields_and_owners() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CO")
                .unwrap()
                .with_2d_coordinates()
                .unwrap();
            let mol = molecule.to_mol().unwrap();
            let sdf = format!("{mol}>  <ID>\nfirst\n\n>  <ID>\nsecond\n\n$$$$\n");
            let params = ck::SdfReadParams {
                sanitize: false,
                remove_hs: false,
                strict_parsing: false,
                expand_attachment_points: true,
                process_property_lists: false,
                coordinate_mode: ck::SdfCoordinateMode::Require3D,
            };
            let native = ck::SdfRecord::from_sdf_with_params(&sdf, &params).unwrap();
            let default = ck::Molecule::from_sdf(&sdf).unwrap().to_binary().unwrap();
            let configured = native.molecule().unwrap().to_binary().unwrap();
            let default_mol = ck::Molecule::from_mol(&mol).unwrap().to_binary().unwrap();
            let configured_mol = ck::Molecule::from_mol_with_params(&mol, &params).unwrap().to_binary().unwrap();
            let query = ck::parse_smarts("C~N").unwrap();
            let properties = ck::MoleculeProperties::default()
                .with_name("query title")
                .with_sdf_data_field("ID", "first")
                .with_sdf_data_field("ID", "second");
            let query_record =
                ck::SdfRecord::from_query_graph(query.clone(), properties.clone()).unwrap();
            let query_sdf = query_record.to_sdf().unwrap();
            let native_query = ck::SdfRecord::from_sdf(&query_sdf).unwrap();
            assert!(native_query.query_graph().is_ok());
            let api = PyModule::new(py, "sdf_record_test")?;
            api.add_class::<SdfRecord>()?;
            api.add_class::<SdfReadParams>()?;
            api.add_class::<SdfError>()?;
            api.add_class::<SdfErrorKind>()?;
            api.add_class::<crate::selectors::SdfCoordinateMode>()?;
            api.add_function(wrap_pyfunction!(mol_from_sdf, &api)?)?;
            api.add_function(wrap_pyfunction!(mol_from_mol, &api)?)?;
            let locals = PyDict::new(py);
            locals.set_item("api", api)?;
            locals.set_item("Molecule", py.get_type::<Molecule>())?;
            locals.set_item("sdf", sdf.clone())?;
            locals.set_item("mol", mol)?;
            locals.set_item("default", default)?;
            locals.set_item("configured", configured)?;
            locals.set_item("default_mol", default_mol)?;
            locals.set_item("configured_mol", configured_mol)?;
            locals.set_item("query_sdf", query_sdf)?;
            locals.set_item("query", Py::new(py, QueryGraph::from_rust(query))?)?;
            locals.set_item(
                "properties",
                Py::new(
                    py,
                    MoleculeProperties {
                        inner: properties.clone(),
                    },
                )?,
            )?;
            py.run(
                c_str!(
                    r#"
params=api.SdfReadParams(sanitize=False,remove_hs=False,strict_parsing=False,
    expand_attachment_points=True,process_property_lists=False,coordinate_mode='require_3d')
kwargs=dict(sanitize=False,remove_hs=False,strict_parsing=False,
    expand_attachment_points=True,process_property_lists=False,coordinate_mode='require_3d')
readers=[(Molecule.from_sdf,sdf,Molecule.from_sdf_with_params,default,configured),
    (Molecule.from_mol,mol,Molecule.from_mol_with_params,default_mol,configured_mol),
    (api.mol_from_sdf,sdf,None,default,configured),(api.mol_from_mol,mol,None,default_mol,configured_mol)]
for reader,text,alias,want_default,want_configured in readers:
    assert bytes(reader(text).to_binary())==bytes(want_default)
    assert bytes(reader(text,params).to_binary())==bytes(want_configured)
    assert bytes(reader(text,**kwargs).to_binary())==bytes(want_configured)
    assert bytes(reader(text,coordinate_mode=api.SdfCoordinateMode.Require3D,
        sanitize=False,remove_hs=False,strict_parsing=False,
        expand_attachment_points=True,process_property_lists=False).to_binary())==bytes(want_configured)
    if alias: assert bytes(alias(text,params).to_binary())==bytes(want_configured)
    try: reader(text,params,sanitize=False)
    except TypeError: pass
    else: raise AssertionError('mixed params accepted')
record=api.SdfRecord.from_sdf(sdf,params)
for other in [api.SdfRecord.from_sdf(sdf,**kwargs),api.SdfRecord.from_sdf_with_params(sdf,params)]:
    assert bytes(other.molecule().to_binary())==bytes(configured)
assert record.index()==0
assert record.data_fields()==[('ID','first'),('ID','second')]
assert record.data_field('ID')=='first' and record.data_field('missing') is None
assert record.properties().sdf_data_fields()==record.data_fields()
assert record.source_coordinate_dim() is not None
assert record.substance_groups()==[]
view=record.graph().value
assert bytes(view.to_binary())==bytes(configured)
view.set_atom_position_(0,[3.0,4.0,5.0])
assert bytes(record.molecule().to_binary())==bytes(configured)
try: record.query_graph()
except api.SdfError as error:
    assert error.variant==api.SdfErrorKind.WrongGraphKind
    assert error.expected=='query_graph' and error.actual=='molecule'
    assert error.__cause__ is None and error.source is None
    try: error.expected='other'
    except AttributeError: pass
    else: raise AssertionError('mutable native error payload')
else: raise AssertionError('concrete record returned query')
qrecord=api.SdfRecord.from_query_graph(query,properties)
assert qrecord.title()=='query title' and qrecord.data_fields()==[('ID','first'),('ID','second')]
assert qrecord.source_coordinate_dim() is None
assert qrecord.properties().name()=='query title'
qview=qrecord.query_graph()
qgraph=qrecord.graph().value
del query,qrecord
assert qview.num_atoms()==2 and qgraph.num_bonds()==1
read_query=api.SdfRecord.from_sdf(query_sdf)
assert read_query.query_graph().num_atoms()==2
try: read_query.molecule()
except api.SdfError as error:
    assert error.expected=='molecule' and error.actual=='query_graph'
else: raise AssertionError('query record lowered')
try: Molecule.from_sdf(query_sdf)
except api.SdfError as error:
    assert error.variant==api.SdfErrorKind.QueryRecord and error.__cause__ is None
else: raise AssertionError('query silently converted')
try: Molecule.from_mol('')
except api.SdfError as error:
    retained=error.source
    assert retained is error.__cause__ and retained.source is None
    assert error.variant==api.SdfErrorKind.Read
else: raise AssertionError('empty MOL accepted')
del record,read_query
assert str(retained)
"#
                ),
                Some(&locals),
                None,
            )?;
            // Borrowed query projections retain this exact immutable graph,
            // not a deep copy or a view whose owner can be released.
            let record = Arc::new(query_record);
            let value = QueryGraph::from_sdf_record(Arc::clone(&record));
            assert!(std::ptr::eq(value.inner(), record.query_graph().unwrap()));
            drop(record);
            assert_eq!(value.inner().num_atoms(), 2);
            for native in [
                ck::SdfError::Construction(ck::OperationError::ReservedAtomPropertyKey {
                    key: "exact\0key".into(),
                }),
                ck::SdfError::Post(ck::MolPostError::MissingCapability("search")),
                ck::SdfError::QueryGraph(ck::QueryGraphError::AdjacencyMismatch),
            ] {
                let display = native.to_string();
                let error = SdfError::from_rust(native);
                assert_eq!(error.value(py).str()?.to_str()?, display);
                assert!(
                    error
                        .value(py)
                        .getattr("source")?
                        .is(&error.value(py).getattr("__cause__")?)
                );
                assert!(!error.value(py).getattr("source")?.is_none());
            }
            let invalid = ck::SdfRecord::from_query_graph(
                ck::parse_smarts("C").unwrap(),
                ck::MoleculeProperties::default().with_name(ck::PropertyText::from_bytes(&[0xff])),
            )
            .unwrap();
            let invalid = Py::new(py, SdfRecord::from_rust(invalid))?;
            let error = invalid.bind(py).call_method0("title").unwrap_err();
            assert!(error.is_instance_of::<pyo3::exceptions::PyUnicodeDecodeError>(py));
            Ok(())
        })
        .unwrap();
    }
}
use crate::embedding_values::CoordinateDimension;
use crate::error_framework::{Source, error_projection};
use crate::group_values::SubstanceGroup;
use crate::metadata_values::MoleculeProperties;
use crate::mol_post_errors::MolPostError;
use crate::native::Molecule;
use crate::operation_errors::OperationError;
use crate::property_values::decode_source_text;
use crate::query_compiler::QueryGraph;
use crate::query_model_errors::QueryGraphError;
use crate::retained_error::RetainedError;
use crate::sdf_read_errors::SdfReadError;
use crate::sdf_values::SdfGraph;
use crate::selectors::{SdfCoordinateModeInput, SdfReadParams};
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;
use std::sync::Arc;

error_projection!(retained SdfError, SdfErrorKind, ck::SdfError, "io", {
    Read(_), Post(_), Construction(_), QueryGraph(_), QueryRecord,
    WrongGraphKind { expected: &'static str, actual: &'static str },
}, cause=|inner: &RetainedError<ck::SdfError>| {
    let source=inner.node().source()?;
    Some(match inner.as_ref() {
        ck::SdfError::Read(value)=>SdfReadError::from_rust(value.clone()),
        ck::SdfError::Post(value)=>MolPostError::from_rust(value.clone()),
        ck::SdfError::Construction(_)=>OperationError::from_ref(source),
        ck::SdfError::QueryGraph(value)=>QueryGraphError::from_rust(value.clone()),
        ck::SdfError::QueryRecord|ck::SdfError::WrongGraphKind{..}=>unreachable!("native SdfError source"),
    })
}, source=Option<Source<Py<SdfReadError>,Source<Py<MolPostError>,Source<Py<OperationError>,Py<QueryGraphError>>>>>);

pub(crate) fn configuration(
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    strict_parsing: Option<bool>,
    expand_attachment_points: Option<bool>,
    process_property_lists: Option<bool>,
    coordinate_mode: Option<SdfCoordinateModeInput>,
) -> PyResult<SdfReadParams> {
    let defaults = ck::SdfReadParams::default();
    Ok(SdfReadParams {
        inner: ck::SdfReadParams {
            sanitize: sanitize.unwrap_or(defaults.sanitize),
            remove_hs: remove_hs.unwrap_or(defaults.remove_hs),
            strict_parsing: strict_parsing.unwrap_or(defaults.strict_parsing),
            expand_attachment_points: expand_attachment_points
                .unwrap_or(defaults.expand_attachment_points),
            process_property_lists: process_property_lists
                .unwrap_or(defaults.process_property_lists),
            coordinate_mode: coordinate_mode
                .map(SdfCoordinateModeInput::into_rust)
                .transpose()?
                .unwrap_or(defaults.coordinate_mode),
        },
    })
}

pub(crate) fn molecule_sdf(
    py: Python<'_>,
    input: &str,
    params: Option<&SdfReadParams>,
) -> PyResult<Molecule> {
    py.detach(|| match params {
        Some(params) => ck::Molecule::from_sdf_with_params(input, &params.inner),
        None => ck::Molecule::from_sdf(input),
    })
    .map(|inner| Molecule { inner })
    .map_err(SdfError::from_rust)
}
pub(crate) fn molecule_mol(
    py: Python<'_>,
    text: &str,
    params: Option<&SdfReadParams>,
) -> PyResult<Molecule> {
    py.detach(|| match params {
        Some(params) => ck::Molecule::from_mol_with_params(text, &params.inner),
        None => ck::Molecule::from_mol(text),
    })
    .map(|inner| Molecule { inner })
    .map_err(SdfError::from_rust)
}

#[pyclass(module = "cosmolkit", frozen)]
pub struct SdfRecord {
    pub(crate) inner: Arc<ck::SdfRecord>,
}
impl SdfRecord {
    pub(crate) fn from_rust(inner: ck::SdfRecord) -> Self {
        Self {
            inner: Arc::new(inner),
        }
    }
}
#[pymethods]
impl SdfRecord {
    #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,format=crate::arguments::Argument::Missing,force_2d=crate::arguments::Argument::Missing,include_stereo=crate::arguments::Argument::Missing,kekulize=crate::arguments::Argument::Missing,precision=crate::arguments::Argument::Missing,coordinate_selection=crate::arguments::Argument::Missing,include_coordinates=crate::arguments::Argument::Missing))]
    pub(crate) fn to_mol(
        &self,
        py: Python<'_>,
        params: crate::arguments::Argument<
            Option<PyRef<'_, crate::io_batch_params::MolBlockWriteParams>>,
        >,
        format: crate::arguments::Argument<Option<crate::io_batch_params::SdfFormatInput>>,
        force_2d: crate::arguments::Argument<Option<bool>>,
        include_stereo: crate::arguments::Argument<Option<bool>>,
        kekulize: crate::arguments::Argument<Option<bool>>,
        precision: crate::arguments::Argument<Option<usize>>,
        coordinate_selection: crate::arguments::Argument<
            Option<crate::io_batch_params::MolCoordinateSelectionInput>,
        >,
        include_coordinates: crate::arguments::Argument<Option<bool>>,
    ) -> PyResult<String> {
        let keyword_form = format.is_supplied()
            || force_2d.is_supplied()
            || include_stereo.is_supplied()
            || kekulize.is_supplied()
            || precision.is_supplied()
            || coordinate_selection.is_supplied()
            || include_coordinates.is_supplied();
        crate::arguments::configuration_arguments!(params; format, force_2d, include_stereo, kekulize, precision, coordinate_selection, include_coordinates);
        let params = params.as_deref();
        if keyword_form {
            self._io_to_mol_keywords(
                py,
                format,
                force_2d,
                include_stereo,
                kekulize,
                precision,
                coordinate_selection,
                include_coordinates,
            )
        } else {
            self._io_to_mol_params(py, params)
        }
    }
    #[doc = r#"{"role":"overload","target":"SdfRecord.to_mol"}"#]
    #[pyo3(signature=(params=None))]
    fn _io_to_mol_params(
        &self,
        py: Python<'_>,
        params: Option<&crate::io_batch_params::MolBlockWriteParams>,
    ) -> PyResult<String> {
        py.detach(|| match params {
            Some(params) => self.inner.to_mol_with_params(&params.inner),
            None => self.inner.to_mol(),
        })
        .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"SdfRecord.to_mol"}"#]
    #[pyo3(signature=(*,format=None,force_2d=None,include_stereo=None,kekulize=None,precision=None,coordinate_selection=None,include_coordinates=None))]
    fn _io_to_mol_keywords(
        &self,
        py: Python<'_>,
        format: Option<crate::io_batch_params::SdfFormatInput>,
        force_2d: Option<bool>,
        include_stereo: Option<bool>,
        kekulize: Option<bool>,
        precision: Option<usize>,
        coordinate_selection: Option<crate::io_batch_params::MolCoordinateSelectionInput>,
        include_coordinates: Option<bool>,
    ) -> PyResult<String> {
        let params = crate::molecular_io_calls::mol_write_keywords(
            format,
            force_2d,
            include_stereo,
            kekulize,
            precision,
            coordinate_selection,
            include_coordinates,
        )?;
        self._io_to_mol_params(py, Some(&params))
    }
    pub(crate) fn to_mol_with_params(
        &self,
        py: Python<'_>,
        params: &crate::io_batch_params::MolBlockWriteParams,
    ) -> PyResult<String> {
        self._io_to_mol_params(py, Some(params))
    }
    #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,format=crate::arguments::Argument::Missing,force_2d=crate::arguments::Argument::Missing,include_stereo=crate::arguments::Argument::Missing,kekulize=crate::arguments::Argument::Missing,precision=crate::arguments::Argument::Missing,coordinate_selection=crate::arguments::Argument::Missing,include_coordinates=crate::arguments::Argument::Missing))]
    pub(crate) fn to_sdf(
        &self,
        py: Python<'_>,
        params: crate::arguments::Argument<
            Option<PyRef<'_, crate::io_batch_params::MolBlockWriteParams>>,
        >,
        format: crate::arguments::Argument<Option<crate::io_batch_params::SdfFormatInput>>,
        force_2d: crate::arguments::Argument<Option<bool>>,
        include_stereo: crate::arguments::Argument<Option<bool>>,
        kekulize: crate::arguments::Argument<Option<bool>>,
        precision: crate::arguments::Argument<Option<usize>>,
        coordinate_selection: crate::arguments::Argument<
            Option<crate::io_batch_params::MolCoordinateSelectionInput>,
        >,
        include_coordinates: crate::arguments::Argument<Option<bool>>,
    ) -> PyResult<String> {
        let keyword_form = format.is_supplied()
            || force_2d.is_supplied()
            || include_stereo.is_supplied()
            || kekulize.is_supplied()
            || precision.is_supplied()
            || coordinate_selection.is_supplied()
            || include_coordinates.is_supplied();
        crate::arguments::configuration_arguments!(params; format, force_2d, include_stereo, kekulize, precision, coordinate_selection, include_coordinates);
        let params = params.as_deref();
        if keyword_form {
            self._io_to_sdf_keywords(
                py,
                format,
                force_2d,
                include_stereo,
                kekulize,
                precision,
                coordinate_selection,
                include_coordinates,
            )
        } else {
            self._io_to_sdf_params(py, params)
        }
    }
    #[doc = r#"{"role":"overload","target":"SdfRecord.to_sdf"}"#]
    #[pyo3(signature=(params=None))]
    fn _io_to_sdf_params(
        &self,
        py: Python<'_>,
        params: Option<&crate::io_batch_params::MolBlockWriteParams>,
    ) -> PyResult<String> {
        py.detach(|| match params {
            Some(params) => self.inner.to_sdf_with_params(&params.inner),
            None => self.inner.to_sdf(),
        })
        .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"SdfRecord.to_sdf"}"#]
    #[pyo3(signature=(*,format=None,force_2d=None,include_stereo=None,kekulize=None,precision=None,coordinate_selection=None,include_coordinates=None))]
    fn _io_to_sdf_keywords(
        &self,
        py: Python<'_>,
        format: Option<crate::io_batch_params::SdfFormatInput>,
        force_2d: Option<bool>,
        include_stereo: Option<bool>,
        kekulize: Option<bool>,
        precision: Option<usize>,
        coordinate_selection: Option<crate::io_batch_params::MolCoordinateSelectionInput>,
        include_coordinates: Option<bool>,
    ) -> PyResult<String> {
        let params = crate::molecular_io_calls::mol_write_keywords(
            format,
            force_2d,
            include_stereo,
            kekulize,
            precision,
            coordinate_selection,
            include_coordinates,
        )?;
        self._io_to_sdf_params(py, Some(&params))
    }
    pub(crate) fn to_sdf_with_params(
        &self,
        py: Python<'_>,
        params: &crate::io_batch_params::MolBlockWriteParams,
    ) -> PyResult<String> {
        self._io_to_sdf_params(py, Some(params))
    }

    #[staticmethod]
    #[pyo3(signature=(input,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
    pub(crate) fn from_sdf(
        py: Python<'_>,
        input: &str,
        params: crate::arguments::Argument<Option<PyRef<'_, SdfReadParams>>>,
        sanitize: crate::arguments::Argument<Option<bool>>,
        remove_hs: crate::arguments::Argument<Option<bool>>,
        strict_parsing: crate::arguments::Argument<Option<bool>>,
        expand_attachment_points: crate::arguments::Argument<Option<bool>>,
        process_property_lists: crate::arguments::Argument<Option<bool>>,
        coordinate_mode: crate::arguments::Argument<
            Option<crate::selectors::SdfCoordinateModeInput>,
        >,
    ) -> PyResult<Self> {
        let keyword_form = sanitize.is_supplied()
            || remove_hs.is_supplied()
            || strict_parsing.is_supplied()
            || expand_attachment_points.is_supplied()
            || process_property_lists.is_supplied()
            || coordinate_mode.is_supplied();
        crate::arguments::configuration_arguments!(params; sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode);
        let params = params.as_deref();
        if keyword_form {
            Self::_from_sdf_keywords(
                py,
                input,
                sanitize,
                remove_hs,
                strict_parsing,
                expand_attachment_points,
                process_property_lists,
                coordinate_mode,
            )
        } else {
            Self::_from_sdf_params(py, input, params)
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"SdfRecord.from_sdf"}"#]
    #[pyo3(signature=(input,params=None))]
    pub(crate) fn _from_sdf_params(
        py: Python<'_>,
        input: &str,
        params: Option<&SdfReadParams>,
    ) -> PyResult<Self> {
        py.detach(|| match params {
            Some(value) => ck::SdfRecord::from_sdf_with_params(input, &value.inner),
            None => ck::SdfRecord::from_sdf(input),
        })
        .map(Self::from_rust)
        .map_err(SdfError::from_rust)
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"SdfRecord.from_sdf"}"#]
    #[pyo3(signature=(input,*,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
    pub(crate) fn _from_sdf_keywords(
        py: Python<'_>,
        input: &str,
        sanitize: Option<bool>,
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
        expand_attachment_points: Option<bool>,
        process_property_lists: Option<bool>,
        coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
    ) -> PyResult<Self> {
        let params = configuration(
            sanitize,
            remove_hs,
            strict_parsing,
            expand_attachment_points,
            process_property_lists,
            coordinate_mode,
        )?;
        Self::_from_sdf_params(py, input, Some(&params))
    }
    #[staticmethod]
    pub(crate) fn from_sdf_with_params(
        py: Python<'_>,
        input: &str,
        params: &SdfReadParams,
    ) -> PyResult<Self> {
        Self::_from_sdf_params(py, input, Some(params))
    }

    #[staticmethod]
    fn from_query_graph(query: &QueryGraph, properties: &MoleculeProperties) -> PyResult<Self> {
        ck::SdfRecord::from_query_graph(query.inner().clone(), properties.inner.clone())
            .map(Self::from_rust)
            .map_err(SdfError::from_rust)
    }
    fn graph(&self) -> SdfGraph {
        match self.inner.graph() {
            ck::SdfGraph::Molecule(value) => SdfGraph::Molecule {
                value: Molecule {
                    inner: value.clone(),
                },
            },
            ck::SdfGraph::Query(_) => SdfGraph::Query {
                value: QueryGraph::from_sdf_record(Arc::clone(&self.inner)),
            },
        }
    }
    fn molecule(&self) -> PyResult<Molecule> {
        self.inner
            .molecule()
            .map(|value| Molecule {
                inner: value.clone(),
            })
            .map_err(SdfError::from_rust)
    }
    fn query_graph(&self) -> PyResult<QueryGraph> {
        self.inner.query_graph().map_err(SdfError::from_rust)?;
        Ok(QueryGraph::from_sdf_record(Arc::clone(&self.inner)))
    }
    fn data_fields(&self, py: Python<'_>) -> PyResult<Vec<(String, String)>> {
        self.inner
            .data_fields()
            .iter()
            .map(|(key, value)| Ok((decode_source_text(py, key)?, decode_source_text(py, value)?)))
            .collect()
    }
    fn properties(&self) -> MoleculeProperties {
        MoleculeProperties {
            inner: self.inner.properties().clone(),
        }
    }
    fn substance_groups(&self) -> Vec<SubstanceGroup> {
        self.inner
            .substance_groups()
            .iter()
            .cloned()
            .map(|inner| SubstanceGroup { inner })
            .collect()
    }
    fn source_coordinate_dim(&self) -> Option<CoordinateDimension> {
        self.inner.source_coordinate_dim().map(Into::into)
    }
    fn index(&self) -> usize {
        self.inner.index()
    }
    fn title(&self, py: Python<'_>) -> PyResult<Option<String>> {
        self.inner
            .title()
            .map(|value| decode_source_text(py, value))
            .transpose()
    }
    fn data_field(&self, py: Python<'_>, name: &str) -> PyResult<Option<String>> {
        self.inner
            .data_field(name)
            .map(|value| decode_source_text(py, value))
            .transpose()
    }
}
#[pyfunction]
#[pyo3(signature=(input,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
pub fn mol_from_sdf(
    py: Python<'_>,
    input: &str,
    params: crate::arguments::Argument<Option<PyRef<'_, SdfReadParams>>>,
    sanitize: crate::arguments::Argument<Option<bool>>,
    remove_hs: crate::arguments::Argument<Option<bool>>,
    strict_parsing: crate::arguments::Argument<Option<bool>>,
    expand_attachment_points: crate::arguments::Argument<Option<bool>>,
    process_property_lists: crate::arguments::Argument<Option<bool>>,
    coordinate_mode: crate::arguments::Argument<Option<crate::selectors::SdfCoordinateModeInput>>,
) -> PyResult<Molecule> {
    let keyword_form = sanitize.is_supplied()
        || remove_hs.is_supplied()
        || strict_parsing.is_supplied()
        || expand_attachment_points.is_supplied()
        || process_property_lists.is_supplied()
        || coordinate_mode.is_supplied();
    crate::arguments::configuration_arguments!(params; sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode);
    let params = params.as_deref();
    if keyword_form {
        _mol_from_sdf_keywords(
            py,
            input,
            sanitize,
            remove_hs,
            strict_parsing,
            expand_attachment_points,
            process_property_lists,
            coordinate_mode,
        )
    } else {
        _mol_from_sdf_params(py, input, params)
    }
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mol_from_sdf"}"#]
#[pyo3(signature=(input,params=None))]
pub fn _mol_from_sdf_params(
    py: Python<'_>,
    input: &str,
    params: Option<&SdfReadParams>,
) -> PyResult<Molecule> {
    molecule_sdf(py, input, params)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mol_from_sdf"}"#]
#[pyo3(signature=(input,*,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
pub fn _mol_from_sdf_keywords(
    py: Python<'_>,
    input: &str,
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    strict_parsing: Option<bool>,
    expand_attachment_points: Option<bool>,
    process_property_lists: Option<bool>,
    coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
) -> PyResult<Molecule> {
    let params = configuration(
        sanitize,
        remove_hs,
        strict_parsing,
        expand_attachment_points,
        process_property_lists,
        coordinate_mode,
    )?;
    _mol_from_sdf_params(py, input, Some(&params))
}

#[pyfunction]
#[pyo3(signature=(text,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
pub fn mol_from_mol(
    py: Python<'_>,
    text: &str,
    params: crate::arguments::Argument<Option<PyRef<'_, SdfReadParams>>>,
    sanitize: crate::arguments::Argument<Option<bool>>,
    remove_hs: crate::arguments::Argument<Option<bool>>,
    strict_parsing: crate::arguments::Argument<Option<bool>>,
    expand_attachment_points: crate::arguments::Argument<Option<bool>>,
    process_property_lists: crate::arguments::Argument<Option<bool>>,
    coordinate_mode: crate::arguments::Argument<Option<crate::selectors::SdfCoordinateModeInput>>,
) -> PyResult<Molecule> {
    let keyword_form = sanitize.is_supplied()
        || remove_hs.is_supplied()
        || strict_parsing.is_supplied()
        || expand_attachment_points.is_supplied()
        || process_property_lists.is_supplied()
        || coordinate_mode.is_supplied();
    crate::arguments::configuration_arguments!(params; sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode);
    let params = params.as_deref();
    if keyword_form {
        _mol_from_mol_keywords(
            py,
            text,
            sanitize,
            remove_hs,
            strict_parsing,
            expand_attachment_points,
            process_property_lists,
            coordinate_mode,
        )
    } else {
        _mol_from_mol_params(py, text, params)
    }
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mol_from_mol"}"#]
#[pyo3(signature=(text,params=None))]
pub fn _mol_from_mol_params(
    py: Python<'_>,
    text: &str,
    params: Option<&SdfReadParams>,
) -> PyResult<Molecule> {
    molecule_mol(py, text, params)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mol_from_mol"}"#]
#[pyo3(signature=(text,*,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
pub fn _mol_from_mol_keywords(
    py: Python<'_>,
    text: &str,
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    strict_parsing: Option<bool>,
    expand_attachment_points: Option<bool>,
    process_property_lists: Option<bool>,
    coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
) -> PyResult<Molecule> {
    let params = configuration(
        sanitize,
        remove_hs,
        strict_parsing,
        expand_attachment_points,
        process_property_lists,
        coordinate_mode,
    )?;
    _mol_from_mol_params(py, text, Some(&params))
}
