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

pyo3::create_exception!(
    cosmolkit,
    SdfError,
    PyValueError,
    "An SDF record could not be parsed, converted or serialized."
);

pub(crate) fn sdf_pyerr(py: Python<'_>, source: ck::SdfError) -> PyErr {
    let kind = match &source {
        ck::SdfError::Read(_) => "Read",
        ck::SdfError::Post(_) => "Post",
        ck::SdfError::Construction(_) => "Construction",
        ck::SdfError::QueryGraph(_) => "QueryGraph",
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

/// Coordinate dimensionality: TwoD denotes XY rows and ThreeD denotes XYZ rows.
///
/// Declared values: ``TwoD``, ``ThreeD``.
#[cosmolkit_macros::python_enum]
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

/// Explicit interpretation requested for the one conformer read from a
/// MolBlock/SDF record.
///
/// Declared values: ``Preserve``, ``Require2D``, ``Require3D``.
#[cosmolkit_macros::python_enum]
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

/// Writable configuration for MOL/SDF parsing, query handling and coordinate retention.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct SdfReadParams {
    pub(crate) inner: ck::SdfReadParams,
}

#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SdfReadParams {
    /// Configure MOL/SDF parsing, query handling and coordinate retention; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature = (*, sanitize=true, remove_hs=true, strict_parsing=true,
        expand_attachment_points=false, process_property_lists=true,
        coordinate_mode=SdfCoordinateMode::Preserve))]
    fn new(
        sanitize: bool,
        remove_hs: bool,
        strict_parsing: bool,
        expand_attachment_points: bool,
        process_property_lists: bool,
        coordinate_mode: SdfCoordinateMode,
    ) -> Self {
        Self {
            inner: ck::SdfReadParams {
                sanitize,
                remove_hs,
                strict_parsing,
                expand_attachment_points,
                process_property_lists,
                coordinate_mode: coordinate_mode.into(),
            },
        }
    }

    /// Return a new result that will perform the selected chemical sanitization stages. The source molecule is unchanged.
    #[getter]
    fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    /// Whether removable explicit hydrogens are removed during input conversion.
    #[getter]
    fn remove_hs(&self) -> bool {
        self.inner.remove_hs
    }
    /// Whether malformed or inconsistent format fields are rejected.
    #[getter]
    fn strict_parsing(&self) -> bool {
        self.inner.strict_parsing
    }
    /// Whether MOL/SDF attachment-point annotations are expanded during parsing.
    #[getter]
    fn expand_attachment_points(&self) -> bool {
        self.inner.expand_attachment_points
    }
    /// Whether SDF atom/bond property lists are interpreted; strict count mismatches raise an error.
    #[getter]
    fn process_property_lists(&self) -> bool {
        self.inner.process_property_lists
    }
    /// Coordinate retention policy: "preserve", "require_2d", or "require_3d"; accepts SdfCoordinateMode too.
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
    /// Return whether this record stores a concrete molecule or a query graph.
    #[getter]
    fn kind(&self) -> &'static str {
        match &self.inner {
            ck::SdfGraph::Molecule(_) => "molecule",
            ck::SdfGraph::Query(_) => "query_graph",
        }
    }

    /// Return the concrete molecule if this is a molecule record; a query remains a query.
    #[getter]
    fn molecule(&self) -> Option<Molecule> {
        match &self.inner {
            ck::SdfGraph::Molecule(value) => Some(Molecule::from_inner(value.clone())),
            ck::SdfGraph::Query(_) => None,
        }
    }

    /// Return the QueryGraph if this is a query record.
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

/// One SDF record with a concrete Molecule or QueryGraph and ordered data fields.
///
/// Query predicates are retained. Data-field names may repeat; their order and
/// raw values are distinct from interpreted atom/bond property lists.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SdfRecord {
    pub(crate) inner: ck::SdfRecord,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SdfRecord {
    /// Original zero-based input record index; failed records retain their index.
    fn index(&self) -> usize {
        self.inner.index()
    }
    /// Return a MOL text block using the selected writer options; does not write a file or install generated drawing coordinates on this object.
    fn to_mol(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_mol()
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Return a MOL text block using the selected writer options; does not write a file or install generated drawing coordinates on this object. Uses the supplied configuration object.
    fn to_mol_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_molecular_io::MolBlockWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_mol_with_params(&params.inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Return an SDF text record using the selected writer options; does not write a file.
    fn to_sdf(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_sdf()
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Return an SDF text record using the selected writer options; does not write a file. Uses the supplied configuration object.
    fn to_sdf_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_molecular_io::MolBlockWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_sdf_with_params(&params.inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Parse one SDF record, retaining a concrete or query graph plus ordered data fields; strict property-list count mismatches raise an error.
    #[staticmethod]
    fn from_sdf(py: Python<'_>, input: &str) -> PyResult<Self> {
        ck::SdfRecord::from_sdf(input)
            .map(|inner| Self { inner })
            .map_err(|error| sdf_pyerr(py, error))
    }

    /// Parse one SDF record, retaining a concrete or query graph plus ordered data fields; strict property-list count mismatches raise an error. Uses the supplied configuration object.
    #[staticmethod]
    fn from_sdf_with_params(py: Python<'_>, input: &str, params: &SdfReadParams) -> PyResult<Self> {
        ck::SdfRecord::from_sdf_with_params(input, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|error| sdf_pyerr(py, error))
    }

    /// Return the concrete Molecule or QueryGraph stored in this record without discarding query predicates.
    fn graph(&self) -> SdfGraph {
        SdfGraph {
            inner: self.inner.graph().clone(),
        }
    }

    /// Return the concrete molecule when the record contains one; query records are not silently converted to concrete molecules.
    fn molecule(&self, py: Python<'_>) -> PyResult<Molecule> {
        self.inner
            .molecule()
            .map(|value| Molecule::from_inner(value.clone()))
            .map_err(|error| sdf_pyerr(py, error))
    }

    /// Return the query graph when the record contains one; concrete records are not silently converted into queries.
    fn query_graph(&self, py: Python<'_>) -> PyResult<QueryGraph> {
        self.inner
            .query_graph()
            .map(|value| QueryGraph {
                inner: value.clone(),
            })
            .map_err(|error| sdf_pyerr(py, error))
    }

    /// Return ordered (name, value) string pairs from the SDF data section, preserving duplicate names.
    fn data_fields(&self, py: Python<'_>) -> PyResult<Vec<(String, String)>> {
        self.inner
            .data_fields()
            .iter()
            .map(|(key, value)| Ok((decode_source_text(py, key)?, decode_source_text(py, value)?)))
            .collect()
    }

    /// Validate the canonical query graph and retain every ordered source field.
    #[staticmethod]
    fn from_query_graph(
        py: Python<'_>,
        query: &QueryGraph,
        properties: &MoleculeProperties,
    ) -> PyResult<Self> {
        ck::SdfRecord::from_query_graph(query.inner.clone(), properties.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|error| sdf_pyerr(py, error))
    }

    /// Title recorded in the source header.
    fn title(&self, py: Python<'_>) -> PyResult<Option<String>> {
        self.inner
            .title()
            .map(|text| decode_source_text(py, text))
            .transpose()
    }

    /// Return the requested SDF field value if present; ordered duplicate fields remain available through data_fields().
    fn data_field(&self, py: Python<'_>, name: &str) -> PyResult<Option<String>> {
        self.inner
            .data_field(name)
            .map(|text| decode_source_text(py, text))
            .transpose()
    }

    /// Return a detached snapshot of the stored properties.
    fn properties(&self) -> MoleculeProperties {
        MoleculeProperties {
            inner: self.inner.properties().clone(),
        }
    }

    /// Substance-group annotations referencing graph atoms and bonds.
    fn substance_groups(&self) -> Vec<SubstanceGroup> {
        self.inner
            .substance_groups()
            .iter()
            .cloned()
            .map(|inner| SubstanceGroup { inner })
            .collect()
    }

    /// Coordinate dimensionality recorded by the source format, when present.
    fn source_coordinate_dim(&self) -> Option<CoordinateDimension> {
        self.inner.source_coordinate_dim().map(Into::into)
    }
}

// Decode only at the Python string boundary. CPython retains the original
// counted bytes on UnicodeDecodeError; chemistry storage never decodes them.
pub(crate) fn decode_source_text(py: Python<'_>, text: &ck::PropertyText) -> PyResult<String> {
    pyo3::types::PyBytes::new(py, text.as_bytes())
        .call_method1("decode", ("utf-8", "strict"))?
        .extract::<String>()
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
