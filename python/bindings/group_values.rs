use crate::chemical_values::spelling;
use crate::property_values::decode_source_text;
use cosmolkit as ck;
use pyo3::prelude::*;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// Local zero-based identifier of a substance group.
    pub struct SubstanceGroupId {
        pub(crate) inner: ck::SubstanceGroupId,
    }
}
#[pymethods]
impl SubstanceGroupId {
    #[new]
    fn construct(index: usize) -> Self {
        Self::new(index)
    }
    #[staticmethod]
    /// Create a local substance-group identifier from a zero-based index.
    fn new(index: usize) -> Self {
        Self {
            inner: ck::SubstanceGroupId::new(index),
        }
    }
    #[getter]
    /// Zero-based substance-group identifier.
    fn index(&self) -> usize {
        self.inner.index()
    }
}

macro_rules! group_kinds {
    ($($variant:ident),+ $(,)?) => {
        crate::configuration::type_projection! {
            #[pyclass(module="cosmolkit", frozen, from_py_object, eq)]
            #[derive(Clone, PartialEq)]
            /// Native substance-group classification, including explicitly named generic groups.
            pub struct SubstanceGroupKind {
                pub(crate) inner: ck::SubstanceGroupKind
            }
        }
        #[pymethods]
        impl SubstanceGroupKind {
            $(#[classattr] #[allow(non_snake_case)]
            fn $variant() -> Self { Self { inner: ck::SubstanceGroupKind::$variant } })+
            #[staticmethod] #[pyo3(name="Generic")]
            /// Create a substance-group kind with an explicit custom name.
            fn generic(value: String) -> Self { Self { inner: ck::SubstanceGroupKind::Generic(value.into()) } }
            #[getter]
            /// Custom group name for Generic groups; None for a standard group kind.
            fn generic_value(&self, py: Python<'_>) -> PyResult<Option<String>> {
                match &self.inner {
                    ck::SubstanceGroupKind::Generic(value) => decode_source_text(py, value).map(Some),
                    $(ck::SubstanceGroupKind::$variant => Ok(None)),+
                }
            }
        }
        impl SubstanceGroupKindInput {
            fn into_rust(self) -> PyResult<ck::SubstanceGroupKind> {
                match self {
                    Self::Kind(value) => Ok(value.inner),
                    Self::Name(value) => {
                        $(if value == spelling(stringify!($variant)) { return Ok(ck::SubstanceGroupKind::$variant); })+
                        Err(pyo3::exceptions::PyValueError::new_err(format!("unknown SubstanceGroupKind: {value}; use Generic(text) for an explicit generic kind")))
                    }
                }
            }
        }
    };
}
#[derive(FromPyObject, IntoPyObject)]
pub enum SubstanceGroupKindInput {
    Kind(SubstanceGroupKind),
    Name(String),
}
group_kinds! { Data, Superatom, MultipleGroup, StructuralRepeatUnit, Monomer,
Copolymer, Crosslink, Graft, Modification, Mer, AnyPolymer, MixtureComponent, Mixture, Formulation }

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// Three three-dimensional points defining an S-group display bracket.
    pub struct SGroupBracket {
        pub(crate) inner: ck::SGroupBracket,
    }
}
#[pymethods]
impl SGroupBracket {
    #[new]
    fn construct(points: [[f64; 3]; 3]) -> Self {
        Self::new(points)
    }
    #[staticmethod]
    /// Create a display bracket from three three-dimensional points.
    fn new(points: [[f64; 3]; 3]) -> Self {
        Self {
            inner: ck::SGroupBracket::new(points),
        }
    }
    #[getter]
    /// Three three-component points defining the bracket geometry.
    fn points(&self) -> [[f64; 3]; 3] {
        *self.inner.points()
    }
}
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// Crossing-bond identifier and associated three-component C-state vector.
    pub struct SGroupCState {
        pub(crate) inner: ck::SGroupCState,
    }
}
#[pymethods]
impl SGroupCState {
    #[new]
    fn construct(bond: usize, vector: [f64; 3]) -> Self {
        Self::new(bond, vector)
    }
    #[staticmethod]
    /// Create a crossing-bond C-state using a zero-based bond index and three-component vector.
    fn new(bond: usize, vector: [f64; 3]) -> Self {
        Self {
            inner: ck::SGroupCState::new(ck::BondId::new(bond), vector),
        }
    }
    #[getter]
    /// Zero-based bond identifier associated with the C-state vector.
    fn bond(&self) -> usize {
        self.inner.bond().index()
    }
    #[getter]
    /// Three-component crossing-bond C-state vector.
    fn vector(&self) -> [f64; 3] {
        *self.inner.vector()
    }
}
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Detached substance-group bracket and label-placement metadata.
    pub struct SGroupDisplay {
        pub(crate) inner: ck::SGroupDisplay,
    }
}
#[pymethods]
impl SGroupDisplay {
    #[new]
    #[pyo3(signature=(*, brackets=None, field_position=None, display_tag=None))]
    fn new(
        brackets: Option<Vec<SGroupBracket>>,
        field_position: Option<[f64; 2]>,
        display_tag: Option<String>,
    ) -> Self {
        Self {
            inner: ck::SGroupDisplay {
                brackets: brackets
                    .map(|values| values.into_iter().map(|value| value.inner).collect())
                    .unwrap_or_else(|| ck::SGroupDisplay::default().brackets),
                field_position,
                display_tag: display_tag.map(Into::into),
            },
        }
    }
    #[getter]
    /// Independent list of substance-group display brackets.
    fn brackets(&self) -> Vec<SGroupBracket> {
        self.inner
            .brackets()
            .iter()
            .copied()
            .map(|inner| SGroupBracket { inner })
            .collect()
    }
    #[getter]
    /// Two-dimensional label position, or None if unspecified.
    fn field_position(&self) -> Option<[f64; 2]> {
        self.inner.field_position
    }
    #[getter]
    /// Source display tag, or None if unspecified.
    fn display_tag(&self, py: Python<'_>) -> PyResult<Option<String>> {
        self.inner
            .display_tag
            .as_ref()
            .map(|value| decode_source_text(py, value))
            .transpose()
    }
}
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Detached substance-group identity, type, crossing-bond and display information.
    pub struct SubstanceGroup {
        pub(crate) inner: ck::SubstanceGroup,
    }
}
#[pymethods]
impl SubstanceGroup {
    #[new]
    fn construct(id: &SubstanceGroupId, kind: SubstanceGroupKindInput) -> PyResult<Self> {
        Self::new(id, kind)
    }
    #[staticmethod]
    /// Create an empty substance group with the supplied identifier and kind.
    fn new(id: &SubstanceGroupId, kind: SubstanceGroupKindInput) -> PyResult<Self> {
        Ok(Self {
            inner: ck::SubstanceGroup::new(id.inner, kind.into_rust()?),
        })
    }
    #[getter]
    /// Identifier of the substance group.
    fn id(&self) -> SubstanceGroupId {
        SubstanceGroupId {
            inner: self.inner.id(),
        }
    }
    #[getter]
    /// Native substance-group classification.
    fn kind(&self) -> SubstanceGroupKind {
        SubstanceGroupKind {
            inner: self.inner.kind().clone(),
        }
    }
    #[getter]
    /// Detached display settings, or None if absent.
    fn display(&self) -> Option<SGroupDisplay> {
        self.inner
            .display()
            .cloned()
            .map(|inner| SGroupDisplay { inner })
    }
    #[getter]
    /// Independent crossing-bond C-state records.
    fn cstates(&self) -> Vec<SGroupCState> {
        self.inner
            .cstates()
            .iter()
            .copied()
            .map(|inner| SGroupCState { inner })
            .collect()
    }
    #[getter]
    /// Zero-based identifiers of head crossing bonds.
    fn head_crossing_bonds(&self) -> Vec<usize> {
        self.inner
            .head_crossing_bonds()
            .iter()
            .copied()
            .map(ck::BondId::index)
            .collect()
    }
    #[getter]
    /// Stored correspondence between crossing bonds.
    fn crossing_bond_correspondence(&self) -> Vec<usize> {
        self.inner
            .crossing_bond_correspondence()
            .iter()
            .copied()
            .map(ck::BondId::index)
            .collect()
    }
}
