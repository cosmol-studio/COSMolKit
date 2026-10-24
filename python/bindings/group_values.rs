use crate::chemical_values::spelling;
use crate::property_values::decode_source_text;
use cosmolkit as ck;
use pyo3::prelude::*;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
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
    fn new(index: usize) -> Self {
        Self {
            inner: ck::SubstanceGroupId::new(index),
        }
    }
    #[getter]
    fn index(&self) -> usize {
        self.inner.index()
    }
}

macro_rules! group_kinds {
    ($($variant:ident),+ $(,)?) => {
        crate::configuration::type_projection! {
            #[pyclass(module="cosmolkit", frozen, from_py_object, eq)]
            #[derive(Clone, PartialEq)]
            pub struct SubstanceGroupKind {
                pub(crate) inner: ck::SubstanceGroupKind
            }
        }
        #[pymethods]
        impl SubstanceGroupKind {
            $(#[classattr] #[allow(non_snake_case)]
            fn $variant() -> Self { Self { inner: ck::SubstanceGroupKind::$variant } })+
            #[staticmethod] #[pyo3(name="Generic")]
            fn generic(value: String) -> Self { Self { inner: ck::SubstanceGroupKind::Generic(value.into()) } }
            #[getter]
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
    fn new(points: [[f64; 3]; 3]) -> Self {
        Self {
            inner: ck::SGroupBracket::new(points),
        }
    }
    #[getter]
    fn points(&self) -> [[f64; 3]; 3] {
        *self.inner.points()
    }
}
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
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
    fn new(bond: usize, vector: [f64; 3]) -> Self {
        Self {
            inner: ck::SGroupCState::new(ck::BondId::new(bond), vector),
        }
    }
    #[getter]
    fn bond(&self) -> usize {
        self.inner.bond().index()
    }
    #[getter]
    fn vector(&self) -> [f64; 3] {
        *self.inner.vector()
    }
}
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
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
    fn brackets(&self) -> Vec<SGroupBracket> {
        self.inner
            .brackets()
            .iter()
            .copied()
            .map(|inner| SGroupBracket { inner })
            .collect()
    }
    #[getter]
    fn field_position(&self) -> Option<[f64; 2]> {
        self.inner.field_position
    }
    #[getter]
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
    fn new(id: &SubstanceGroupId, kind: SubstanceGroupKindInput) -> PyResult<Self> {
        Ok(Self {
            inner: ck::SubstanceGroup::new(id.inner, kind.into_rust()?),
        })
    }
    #[getter]
    fn id(&self) -> SubstanceGroupId {
        SubstanceGroupId {
            inner: self.inner.id(),
        }
    }
    #[getter]
    fn kind(&self) -> SubstanceGroupKind {
        SubstanceGroupKind {
            inner: self.inner.kind().clone(),
        }
    }
    #[getter]
    fn display(&self) -> Option<SGroupDisplay> {
        self.inner
            .display()
            .cloned()
            .map(|inner| SGroupDisplay { inner })
    }
    #[getter]
    fn cstates(&self) -> Vec<SGroupCState> {
        self.inner
            .cstates()
            .iter()
            .copied()
            .map(|inner| SGroupCState { inner })
            .collect()
    }
    #[getter]
    fn head_crossing_bonds(&self) -> Vec<usize> {
        self.inner
            .head_crossing_bonds()
            .iter()
            .copied()
            .map(ck::BondId::index)
            .collect()
    }
    #[getter]
    fn crossing_bond_correspondence(&self) -> Vec<usize> {
        self.inner
            .crossing_bond_correspondence()
            .iter()
            .copied()
            .map(ck::BondId::index)
            .collect()
    }
}
