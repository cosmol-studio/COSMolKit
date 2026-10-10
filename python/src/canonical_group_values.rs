//! Detached group values copied through the canonical facade only.
use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct SubstanceGroupId {
    pub(crate) inner: ck::SubstanceGroupId,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SubstanceGroupId {
    #[staticmethod]
    fn new(index: usize) -> Self {
        Self {
            inner: ck::SubstanceGroupId::new(index),
        }
    }
    fn index(&self) -> usize {
        self.inner.index()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct SubstanceGroupKind {
    pub(crate) inner: ck::SubstanceGroupKind,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SubstanceGroupKind {
    #[classattr]
    #[pyo3(name = "Data")]
    fn variant_data() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::Data,
        }
    }
    #[classattr]
    #[pyo3(name = "Superatom")]
    fn variant_superatom() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::Superatom,
        }
    }
    #[classattr]
    #[pyo3(name = "MultipleGroup")]
    fn variant_multiplegroup() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::MultipleGroup,
        }
    }
    #[classattr]
    #[pyo3(name = "StructuralRepeatUnit")]
    fn variant_structuralrepeatunit() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::StructuralRepeatUnit,
        }
    }
    #[classattr]
    #[pyo3(name = "Monomer")]
    fn variant_monomer() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::Monomer,
        }
    }
    #[classattr]
    #[pyo3(name = "Copolymer")]
    fn variant_copolymer() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::Copolymer,
        }
    }
    #[classattr]
    #[pyo3(name = "Crosslink")]
    fn variant_crosslink() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::Crosslink,
        }
    }
    #[classattr]
    #[pyo3(name = "Graft")]
    fn variant_graft() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::Graft,
        }
    }
    #[classattr]
    #[pyo3(name = "Modification")]
    fn variant_modification() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::Modification,
        }
    }
    #[classattr]
    #[pyo3(name = "Mer")]
    fn variant_mer() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::Mer,
        }
    }
    #[classattr]
    #[pyo3(name = "AnyPolymer")]
    fn variant_anypolymer() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::AnyPolymer,
        }
    }
    #[classattr]
    #[pyo3(name = "MixtureComponent")]
    fn variant_mixturecomponent() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::MixtureComponent,
        }
    }
    #[classattr]
    #[pyo3(name = "Mixture")]
    fn variant_mixture() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::Mixture,
        }
    }
    #[classattr]
    #[pyo3(name = "Formulation")]
    fn variant_formulation() -> SubstanceGroupKind {
        Self {
            inner: ck::SubstanceGroupKind::Formulation,
        }
    }
    #[staticmethod]
    #[pyo3(name = "Generic")]
    fn generic(value: String) -> Self {
        Self {
            inner: ck::SubstanceGroupKind::Generic(value.into()),
        }
    }
    #[getter]
    fn generic_value(&self, py: Python<'_>) -> PyResult<Option<String>> {
        match &self.inner {
            ck::SubstanceGroupKind::Generic(value) => {
                // Decode only at the Python str projection. The native decode
                // error carries original bytes; no lossy chemistry conversion.
                let text = pyo3::types::PyBytes::new(py, value.as_bytes())
                    .call_method1("decode", ("utf-8", "strict"))?
                    .extract::<String>()?;
                Ok(Some(text))
            }
            _ => Ok(None),
        }
    }
    fn __eq__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match other.extract::<PyRef<'_, Self>>() {
            Ok(other) => Ok((self.inner == other.inner)
                .into_pyobject(py)?
                .to_owned()
                .into_any()
                .unbind()),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct SGroupBracket {
    pub(crate) inner: ck::SGroupBracket,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SGroupBracket {
    #[staticmethod]
    fn new(points: [[f64; 3]; 3]) -> Self {
        Self {
            inner: ck::SGroupBracket::new(points),
        }
    }
    fn points(&self) -> [[f64; 3]; 3] {
        *self.inner.points()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct SGroupCState {
    pub(crate) inner: ck::SGroupCState,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SGroupCState {
    #[staticmethod]
    fn new(bond: usize, vector: [f64; 3]) -> Self {
        Self {
            inner: ck::SGroupCState::new(ck::BondId::new(bond), vector),
        }
    }
    fn bond(&self) -> usize {
        self.inner.bond().index()
    }
    fn vector(&self) -> [f64; 3] {
        *self.inner.vector()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct SGroupDisplay {
    pub(crate) inner: ck::SGroupDisplay,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SGroupDisplay {
    fn brackets(&self) -> Vec<SGroupBracket> {
        self.inner
            .brackets()
            .iter()
            .copied()
            .map(|inner| SGroupBracket { inner })
            .collect()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct SubstanceGroup {
    pub(crate) inner: ck::SubstanceGroup,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SubstanceGroup {
    #[staticmethod]
    fn new(id: &SubstanceGroupId, kind: &SubstanceGroupKind) -> Self {
        Self {
            inner: ck::SubstanceGroup::new(id.inner, kind.inner.clone()),
        }
    }
    fn id(&self) -> SubstanceGroupId {
        SubstanceGroupId {
            inner: self.inner.id(),
        }
    }
    fn kind(&self) -> SubstanceGroupKind {
        SubstanceGroupKind {
            inner: self.inner.kind().clone(),
        }
    }
    fn display(&self) -> Option<SGroupDisplay> {
        self.inner
            .display()
            .cloned()
            .map(|inner| SGroupDisplay { inner })
    }
    fn cstates(&self) -> Vec<SGroupCState> {
        self.inner
            .cstates()
            .iter()
            .copied()
            .map(|inner| SGroupCState { inner })
            .collect()
    }
    fn head_crossing_bonds(&self) -> Vec<usize> {
        self.inner
            .head_crossing_bonds()
            .iter()
            .map(|id| id.index())
            .collect()
    }
    fn crossing_bond_correspondence(&self) -> Vec<usize> {
        self.inner
            .crossing_bond_correspondence()
            .iter()
            .map(|id| id.index())
            .collect()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct TemplateAttachment {
    pub(crate) inner: ck::TemplateAttachment,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TemplateAttachment {
    fn target(&self) -> usize {
        self.inner.target().index()
    }
    fn label(&self) -> &str {
        self.inner.label()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct TemplateAttachmentOrder {
    pub(crate) inner: ck::TemplateAttachmentOrder,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TemplateAttachmentOrder {
    fn entries(&self) -> Vec<TemplateAttachment> {
        self.inner
            .entries()
            .iter()
            .cloned()
            .map(|inner| TemplateAttachment { inner })
            .collect()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct StereoGroup {
    pub(crate) inner: ck::StereoGroup,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl StereoGroup {
    fn id(&self) -> Option<u32> {
        self.inner.id()
    }
    fn atoms(&self) -> Vec<usize> {
        self.inner.atoms().iter().map(|id| id.index()).collect()
    }
    fn bonds(&self) -> Vec<usize> {
        self.inner.bonds().iter().map(|id| id.index()).collect()
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<SubstanceGroupId>()?;
    module.add_class::<SubstanceGroupKind>()?;
    module.add_class::<SGroupBracket>()?;
    module.add_class::<SGroupCState>()?;
    module.add_class::<SGroupDisplay>()?;
    module.add_class::<SubstanceGroup>()?;
    module.add_class::<TemplateAttachment>()?;
    module.add_class::<TemplateAttachmentOrder>()?;
    module.add_class::<StereoGroup>()?;
    Ok(())
}
