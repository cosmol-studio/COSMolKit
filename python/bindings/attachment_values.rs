use cosmolkit as ck;
use pyo3::prelude::*;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct TemplateAttachment {
        pub(crate) inner: ck::TemplateAttachment,
    }
}
#[pymethods]
impl TemplateAttachment {
    #[new]
    fn new(target: usize, label: String) -> Self {
        Self {
            inner: ck::TemplateAttachment::new(ck::AtomId::new(target), label),
        }
    }
    #[getter]
    fn target(&self) -> usize {
        self.inner.target().index()
    }
    #[getter]
    fn label(&self) -> &str {
        self.inner.label()
    }
}
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct TemplateAttachmentOrder {
        pub(crate) inner: ck::TemplateAttachmentOrder,
    }
}
#[pymethods]
impl TemplateAttachmentOrder {
    #[new]
    fn new(entries: Vec<TemplateAttachment>) -> PyResult<Self> {
        ck::TemplateAttachmentOrder::new(entries.into_iter().map(|value| value.inner).collect())
            .map(|inner| Self { inner })
            .map_err(TemplateAttachmentOrderError::from_rust)
    }
    #[getter]
    fn entries(&self) -> Vec<TemplateAttachment> {
        self.inner
            .entries()
            .iter()
            .cloned()
            .map(|inner| TemplateAttachment { inner })
            .collect()
    }
    fn validate_for_atom_count(&self, atom_count: usize) -> PyResult<()> {
        self.inner
            .validate_for_atom_count(atom_count)
            .map_err(TemplateAttachmentOrderError::from_rust)
    }
}

crate::enums::binding_enum!(TemplateAttachmentOrderErrorKind, {
    Empty,
    DuplicateTarget,
    DuplicateLabel,
    TargetOutOfRange,
    MappingTargetOutOfRange,
    TargetRemoved,
});
crate::error_framework::error_projection! {record
pub struct TemplateAttachmentOrderError {
    inner: ck::TemplateAttachmentOrderError,
}
}
crate::error_framework::error_projection! {methods impl TemplateAttachmentOrderError {
    #[getter]
    fn domain(&self) -> &'static str {
        "model"
    }
    #[getter]
    fn variant(&self) -> TemplateAttachmentOrderErrorKind {
        use TemplateAttachmentOrderErrorKind as K;
        use ck::TemplateAttachmentOrderError as E;
        match &self.inner {
            E::Empty => K::Empty,
            E::DuplicateTarget {
                first_position: _,
                duplicate_position: _,
                target: _,
            } => K::DuplicateTarget,
            E::DuplicateLabel {
                first_position: _,
                duplicate_position: _,
                label: _,
            } => K::DuplicateLabel,
            E::TargetOutOfRange {
                position: _,
                target: _,
                atom_count: _,
            } => K::TargetOutOfRange,
            E::MappingTargetOutOfRange {
                position: _,
                target: _,
                mapping_len: _,
            } => K::MappingTargetOutOfRange,
            E::TargetRemoved {
                position: _,
                target: _,
            } => K::TargetRemoved,
        }
    }
    #[getter]
    fn first_position(&self) -> Option<usize> {
        match self.inner {
            ck::TemplateAttachmentOrderError::DuplicateTarget { first_position, .. }
            | ck::TemplateAttachmentOrderError::DuplicateLabel { first_position, .. } => {
                Some(first_position)
            }
            _ => None,
        }
    }
    #[getter]
    fn duplicate_position(&self) -> Option<usize> {
        match self.inner {
            ck::TemplateAttachmentOrderError::DuplicateTarget {
                duplicate_position, ..
            }
            | ck::TemplateAttachmentOrderError::DuplicateLabel {
                duplicate_position, ..
            } => Some(duplicate_position),
            _ => None,
        }
    }
    #[getter]
    fn target(&self) -> Option<usize> {
        match self.inner {
            ck::TemplateAttachmentOrderError::DuplicateTarget { target, .. }
            | ck::TemplateAttachmentOrderError::TargetOutOfRange { target, .. }
            | ck::TemplateAttachmentOrderError::MappingTargetOutOfRange { target, .. }
            | ck::TemplateAttachmentOrderError::TargetRemoved { target, .. } => {
                Some(target.index())
            }
            _ => None,
        }
    }
    #[getter]
    fn position(&self) -> Option<usize> {
        match self.inner {
            ck::TemplateAttachmentOrderError::TargetOutOfRange { position, .. }
            | ck::TemplateAttachmentOrderError::MappingTargetOutOfRange { position, .. }
            | ck::TemplateAttachmentOrderError::TargetRemoved { position, .. } => Some(position),
            _ => None,
        }
    }
    #[getter]
    fn label(&self) -> Option<&str> {
        match &self.inner {
            ck::TemplateAttachmentOrderError::DuplicateLabel { label, .. } => Some(label),
            _ => None,
        }
    }
    #[getter]
    fn atom_count(&self) -> Option<usize> {
        match self.inner {
            ck::TemplateAttachmentOrderError::TargetOutOfRange { atom_count, .. } => {
                Some(atom_count)
            }
            _ => None,
        }
    }
    #[getter]
    fn mapping_len(&self) -> Option<usize> {
        match self.inner {
            ck::TemplateAttachmentOrderError::MappingTargetOutOfRange { mapping_len, .. } => {
                Some(mapping_len)
            }
            _ => None,
        }
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
impl TemplateAttachmentOrderError {
    pub(crate) fn from_rust(inner: ck::TemplateAttachmentOrderError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => crate::error_framework::from_payload(value),
            Err(error) => error,
        })
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    pub struct AtomMetadata {
        pub(crate) inner: ck::AtomMetadata,
    }
}
#[pymethods]
impl AtomMetadata {
    #[new]
    fn new(
        degree: usize,
        explicit_valence: i32,
        implicit_hydrogens: i32,
        total_hydrogens: i32,
        total_valence: i32,
    ) -> Self {
        Self {
            inner: ck::AtomMetadata {
                degree,
                explicit_valence,
                implicit_hydrogens,
                total_hydrogens,
                total_valence,
            },
        }
    }
    #[getter]
    fn degree(&self) -> usize {
        self.inner.degree
    }
    #[getter]
    fn explicit_valence(&self) -> i32 {
        self.inner.explicit_valence
    }
    #[getter]
    fn implicit_hydrogens(&self) -> i32 {
        self.inner.implicit_hydrogens
    }
    #[getter]
    fn total_hydrogens(&self) -> i32 {
        self.inner.total_hydrogens
    }
    #[getter]
    fn total_valence(&self) -> i32 {
        self.inner.total_valence
    }
}
