use crate::chemical_values::{EnumValue, enum_value};
use crate::coordinates;
use cosmolkit as ck;
use numpy::IntoPyArray;
use pyo3::prelude::*;

enum_value!(ChainKind, {
    Protein, Dna, Rna, ProteinDnaComplex, ProteinRnaComplex, LigandOnly,
    WaterOnly, Mixed, Unknown,
});
enum_value!(BioCoordinateFormat, input = BioCoordinateFormatInput, { Unknown, Detect, Pdb, Mmcif, Mmjson, ChemComp });

#[derive(FromPyObject, IntoPyObject)]
pub enum BioCoordinateFormatInput {
    Enum(BioCoordinateFormat),
    Name(String),
    Code(u8),
}
impl BioCoordinateFormatInput {
    pub(crate) fn into_rust(self) -> PyResult<ck::BioCoordinateFormat> {
        use ck::BioCoordinateFormat as F;
        match self {
            Self::Enum(value) => Ok(value.into()),
            Self::Name(name) => match name.as_str() {
                "unknown" => Ok(F::Unknown),
                "detect" => Ok(F::Detect),
                "pdb" => Ok(F::Pdb),
                "mmcif" => Ok(F::Mmcif),
                "mmjson" => Ok(F::Mmjson),
                "chemcomp" => Ok(F::ChemComp),
                _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "unknown BioCoordinateFormat: {name}"
                ))),
            },
            Self::Code(code) => match code {
                0 => Ok(F::Unknown),
                1 => Ok(F::Detect),
                2 => Ok(F::Pdb),
                3 => Ok(F::Mmcif),
                4 => Ok(F::Mmjson),
                5 => Ok(F::ChemComp),
                _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "invalid BioCoordinateFormat code: {code}"
                ))),
            },
        }
    }
}
ck::project_configuration!(ck::BioReadParams);

#[derive(FromPyObject, IntoPyObject)]
pub enum ChainKindInput {
    Enum(ChainKind),
    Name(String),
}
impl ChainKindInput {
    fn into_rust(self) -> PyResult<ck::ChainKind> {
        match self {
            Self::Enum(value) => Ok(value.into()),
            Self::Name(name) => ChainKind::from_name(&name).map(Into::into),
        }
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Author and label chain identifiers retained from structural coordinate files.
    pub struct ChainSourceIds {
        pub(crate) inner: ck::ChainSourceIds,
    }
}
#[pymethods]
impl ChainSourceIds {
    #[new]
    #[pyo3(signature = (auth_chain_id, label_asym_id))]
    fn new(auth_chain_id: Option<&str>, label_asym_id: Option<String>) -> PyResult<Self> {
        let auth_chain_id = auth_chain_id
            .map(|name| {
                ck::PdbChainId::from_ascii(name.as_bytes()).ok_or_else(|| {
                    pyo3::exceptions::PyValueError::new_err(
                        "auth_chain_id is not a valid PdbChainId",
                    )
                })
            })
            .transpose()?;
        Ok(Self {
            inner: ck::ChainSourceIds::new(auth_chain_id, label_asym_id),
        })
    }
    #[getter]
    /// Author-provided chain identifier, or None if absent.
    fn auth_chain_id(&self) -> Option<String> {
        self.inner
            .auth_chain_id()
            .map(|value| value.as_str().to_owned())
    }
    #[getter]
    /// mmCIF label_asym_id identifying the asymmetric unit chain, or None.
    fn label_asym_id(&self) -> Option<&str> {
        self.inner.label_asym_id()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Detached chain row containing parent identifiers and a residue span, without an owning structure.
    pub struct BioChainRow {
        pub(crate) inner: ck::BioChainRow,
    }
}
#[pymethods]
impl BioChainRow {
    #[new]
    #[pyo3(signature = (model_id, entity_id, residue_span, kind, source))]
    fn new(
        model_id: u32,
        entity_id: Option<u32>,
        residue_span: &BioRowSpan,
        kind: ChainKindInput,
        source: &ChainSourceIds,
    ) -> PyResult<Self> {
        let span = residue_span
            .as_rust()
            .map_err(crate::bio_errors::BioStructureError::from_rust)?;
        Ok(Self {
            inner: ck::BioChainRow::new(
                ck::BioModelId::new(model_id),
                entity_id.map(ck::BioEntityId::new),
                span,
                kind.into_rust()?,
                source.inner.clone(),
            ),
        })
    }
    #[getter]
    /// Local index of the parent model row.
    fn model_id(&self) -> u32 {
        self.inner.model_id().value()
    }
    #[getter]
    /// Local entity row identifier, or None if unassigned.
    fn entity_id(&self) -> Option<u32> {
        self.inner.entity_id().map(|id| id.value())
    }
    #[getter]
    /// Half-open range of residue rows belonging to this chain.
    fn residue_span(&self) -> PyResult<BioRowSpan> {
        BioRowSpan::from_rust(self.inner.residue_span())
            .map_err(crate::bio_errors::BioStructureError::from_rust)
    }
    #[getter]
    /// Stored chain classification.
    fn kind(&self) -> ChainKind {
        self.inner.kind().into()
    }
    #[getter]
    /// Source-format identifiers, distinct from local row indices.
    fn source(&self) -> ChainSourceIds {
        ChainSourceIds {
            inner: self.inner.source().clone(),
        }
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Checked half-open range of rows represented by an inclusive start and row count.
    pub struct BioRowSpan {
        pub(crate) inner: ck::BioRowSpan<ck::BioAtomId>,
    }
}

impl BioRowSpan {
    pub(crate) fn as_rust<I>(&self) -> Result<ck::BioRowSpan<I>, ck::BioStructureError> {
        ck::BioRowSpan::new(self.inner.start(), self.inner.len())
    }
    pub(crate) fn from_rust<I>(value: ck::BioRowSpan<I>) -> Result<Self, ck::BioStructureError> {
        ck::BioRowSpan::new(value.start(), value.len()).map(|inner| Self { inner })
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Detached model row containing its chain span and optional source model number.
    pub struct BioModelRow {
        pub(crate) inner: ck::BioModelRow,
    }
}
#[pymethods]
impl BioModelRow {
    #[new]
    #[pyo3(signature = (chain_span, source_model_number))]
    fn new(chain_span: &BioRowSpan, source_model_number: Option<i32>) -> PyResult<Self> {
        chain_span
            .as_rust()
            .map(|span| Self {
                inner: ck::BioModelRow::new(span, source_model_number),
            })
            .map_err(crate::bio_errors::BioStructureError::from_rust)
    }
    #[getter]
    /// Half-open range of chain rows belonging to this model.
    fn chain_span(&self) -> PyResult<BioRowSpan> {
        BioRowSpan::from_rust(self.inner.chain_span())
            .map_err(crate::bio_errors::BioStructureError::from_rust)
    }
    #[getter]
    /// Model number from the input file, or None if absent.
    fn source_model_number(&self) -> Option<i32> {
        self.inner.source_model_number()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Validated atom-name text preserving the native structural-file representation.
    pub struct AtomName {
        pub(crate) inner: ck::AtomName,
    }
}
#[pymethods]
impl AtomName {
    #[staticmethod]
    /// Construct a validated atom name from ASCII bytes; return None if invalid.
    fn from_ascii(bytes: &[u8]) -> Option<Self> {
        ck::AtomName::from_ascii(bytes).map(|inner| Self { inner })
    }
    /// Atom name as text.
    fn as_str(&self) -> &str {
        self.inner.as_str()
    }
    /// Atom name as an independent bytes value.
    fn as_bytes<'py>(&self, py: Python<'py>) -> Bound<'py, pyo3::types::PyBytes> {
        pyo3::types::PyBytes::new(py, self.inner.as_bytes())
    }
}
#[pymethods]
impl BioRowSpan {
    #[new]
    fn new(start: u32, len: u32) -> PyResult<Self> {
        ck::BioRowSpan::new(start, len)
            .map(|inner| Self { inner })
            .map_err(crate::bio_errors::BioStructureError::from_rust)
    }
    #[staticmethod]
    /// Create a checked row span from platform-sized start and length values.
    fn from_usize(start: usize, len: usize) -> PyResult<Self> {
        ck::BioRowSpan::from_usize(start, len)
            .map(|inner| Self { inner })
            .map_err(crate::bio_errors::BioStructureError::from_rust)
    }
    #[getter]
    /// Inclusive starting row index.
    fn start(&self) -> u32 {
        self.inner.start()
    }
    #[getter]
    /// Number of rows in the span.
    fn len(&self) -> u32 {
        self.inner.len()
    }
    #[getter]
    /// Whether the span contains zero rows.
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    #[getter]
    /// Exclusive ending row index.
    fn end(&self) -> u32 {
        self.inner.end()
    }
}

ck::project_configuration!(ck::BioCrystalCell);

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// Three-dimensional affine transform with a 3-by-3 matrix and three-component translation.
    pub struct BioTransform {
        pub(crate) inner: ck::BioTransform,
    }
}
#[pymethods]
impl BioTransform {
    #[new]
    #[pyo3(signature = (matrix=*ck::BioTransform::default().matrix(), translation=*ck::BioTransform::default().translation()))]
    fn new(matrix: [[f64; 3]; 3], translation: [f64; 3]) -> Self {
        Self {
            inner: ck::BioTransform::new(matrix, translation),
        }
    }
    /// Whether matrix and translation components match within the supplied native tolerance.
    fn approx(&self, other: &BioTransform, epsilon: f64) -> bool {
        self.inner.approx(&other.inner, epsilon)
    }
    #[getter]
    #[doc = r#"{"role": "return_type", "imports": ["numpy"], "type": "numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]]"}"#]
    /// Independent float64 NumPy array of shape (3, 3) containing the linear transform.
    fn matrix<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        coordinates::array(py, self.inner.matrix())
    }
    #[getter]
    #[doc = r#"{"role": "return_type", "imports": ["numpy"], "type": "numpy.ndarray[tuple[int], numpy.dtype[numpy.float64]]"}"#]
    /// Independent float64 NumPy array of shape (3,) containing the translation vector.
    fn translation<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<f64>> {
        self.inner.translation().to_vec().into_pyarray(py)
    }
}
