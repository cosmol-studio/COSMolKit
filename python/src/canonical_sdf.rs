//! SDF projections use the finalized facade record and its tagged graph.
use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

use crate::canonical_group_values::SubstanceGroup;
use crate::canonical_property_values::MoleculeProperties;
use crate::canonical_search::QueryGraph;
use crate::drawing_binding::Molecule;

pyo3::create_exception!(cosmolkit, SdfError, PyValueError);

pub(crate) fn sdf_pyerr(py: Python<'_>, source: ck::SdfError) -> PyErr {
    let kind = match &source {
        ck::SdfError::Read(_) => "Read",
        ck::SdfError::Post(_) => "Post",
        ck::SdfError::Construction(_) => "Construction",
        ck::SdfError::QueryRecord => "QueryRecord",
        ck::SdfError::WrongGraphKind { .. } => "WrongGraphKind",
    };
    let error = crate::canonical_values::annotate(
        py,
        SdfError::new_err(source.to_string()),
        "io",
        kind,
        &source,
    );
    if let ck::SdfError::WrongGraphKind { expected, actual } = source {
        if let Err(detail) = error
            .value(py)
            .setattr("expected", expected)
            .and_then(|()| error.value(py).setattr("actual", actual))
        {
            return detail;
        }
    }
    error
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum CoordinateDimension {
    TwoD,
    ThreeD,
}

impl From<ck::CoordinateDimension> for CoordinateDimension {
    fn from(value: ck::CoordinateDimension) -> Self {
        match value {
            ck::CoordinateDimension::TwoD => Self::TwoD,
            ck::CoordinateDimension::ThreeD => Self::ThreeD,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum SdfCoordinateMode {
    Preserve,
    Require2D,
    Require3D,
}

impl From<SdfCoordinateMode> for ck::SdfCoordinateMode {
    fn from(value: SdfCoordinateMode) -> Self {
        match value {
            SdfCoordinateMode::Preserve => Self::Preserve,
            SdfCoordinateMode::Require2D => Self::Require2D,
            SdfCoordinateMode::Require3D => Self::Require3D,
        }
    }
}

impl From<ck::SdfCoordinateMode> for SdfCoordinateMode {
    fn from(value: ck::SdfCoordinateMode) -> Self {
        match value {
            ck::SdfCoordinateMode::Preserve => Self::Preserve,
            ck::SdfCoordinateMode::Require2D => Self::Require2D,
            ck::SdfCoordinateMode::Require3D => Self::Require3D,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SdfReadParams {
    pub(crate) inner: ck::SdfReadParams,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SdfReadParams {
    #[new]
    #[pyo3(signature = (*, sanitize=true, remove_hydrogens=true, strict_parsing=true,
        expand_attachment_points=false, process_property_lists=true,
        coordinate_mode=SdfCoordinateMode::Preserve))]
    fn new(
        sanitize: bool,
        remove_hydrogens: bool,
        strict_parsing: bool,
        expand_attachment_points: bool,
        process_property_lists: bool,
        coordinate_mode: SdfCoordinateMode,
    ) -> Self {
        Self {
            inner: ck::SdfReadParams {
                sanitize,
                remove_hydrogens,
                strict_parsing,
                expand_attachment_points,
                process_property_lists,
                coordinate_mode: coordinate_mode.into(),
            },
        }
    }

    #[getter]
    fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    #[getter]
    fn remove_hydrogens(&self) -> bool {
        self.inner.remove_hydrogens
    }
    #[getter]
    fn strict_parsing(&self) -> bool {
        self.inner.strict_parsing
    }
    #[getter]
    fn expand_attachment_points(&self) -> bool {
        self.inner.expand_attachment_points
    }
    #[getter]
    fn process_property_lists(&self) -> bool {
        self.inner.process_property_lists
    }
    #[getter]
    fn coordinate_mode(&self) -> SdfCoordinateMode {
        self.inner.coordinate_mode.into()
    }
}

/// This tag wraps an existing Molecule or QueryGraph; it never lowers a query.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SdfGraph {
    inner: ck::SdfGraph,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SdfGraph {
    #[getter]
    fn kind(&self) -> &'static str {
        match &self.inner {
            ck::SdfGraph::Molecule(_) => "molecule",
            ck::SdfGraph::Query(_) => "query_graph",
        }
    }

    #[getter]
    fn molecule(&self) -> Option<Molecule> {
        match &self.inner {
            ck::SdfGraph::Molecule(value) => Some(Molecule::from_inner(value.clone())),
            ck::SdfGraph::Query(_) => None,
        }
    }

    #[getter]
    fn query_graph(&self) -> Option<QueryGraph> {
        match &self.inner {
            ck::SdfGraph::Query(value) => Some(QueryGraph {
                inner: value.clone(),
            }),
            ck::SdfGraph::Molecule(_) => None,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SdfRecord {
    inner: ck::SdfRecord,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SdfRecord {
    #[staticmethod]
    fn from_sdf(py: Python<'_>, input: &str) -> PyResult<Self> {
        ck::SdfRecord::from_sdf(input)
            .map(|inner| Self { inner })
            .map_err(|error| sdf_pyerr(py, error))
    }

    #[staticmethod]
    fn from_sdf_with_params(py: Python<'_>, input: &str, params: &SdfReadParams) -> PyResult<Self> {
        ck::SdfRecord::from_sdf_with_params(input, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|error| sdf_pyerr(py, error))
    }

    fn graph(&self) -> SdfGraph {
        SdfGraph {
            inner: self.inner.graph().clone(),
        }
    }

    fn molecule(&self, py: Python<'_>) -> PyResult<Molecule> {
        self.inner
            .molecule()
            .map(|value| Molecule::from_inner(value.clone()))
            .map_err(|error| sdf_pyerr(py, error))
    }

    fn query_graph(&self, py: Python<'_>) -> PyResult<QueryGraph> {
        self.inner
            .query_graph()
            .map(|value| QueryGraph {
                inner: value.clone(),
            })
            .map_err(|error| sdf_pyerr(py, error))
    }

    fn data_fields(&self) -> Vec<(String, String)> {
        self.inner.data_fields().to_vec()
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
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<CoordinateDimension>()?;
    module.add_class::<SdfCoordinateMode>()?;
    module.add_class::<SdfReadParams>()?;
    module.add_class::<SdfGraph>()?;
    module.add_class::<SdfRecord>()?;
    module.add("SdfError", module.py().get_type::<SdfError>())?;
    Ok(())
}
