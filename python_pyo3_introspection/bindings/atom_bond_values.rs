use crate::chemical_values::{
    BondDirection, BondOrder, BondStereo, ChiralTag, Hybridization, enum_value,
};
use crate::elements::Element;
use crate::property_values::{PropertyValueError, decode_source_text};
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

enum_value!(CipDescriptor, { R, S, LowerR, LowerS, E, Z, LowerE, LowerZ, M, P, LowerM, LowerP });
#[pymethods]
impl CipDescriptor {
    #[new]
    fn new(value: &str) -> PyResult<Self> {
        value
            .parse::<ck::CipDescriptor>()
            .map(Into::into)
            .map_err(CipDescriptorError::from_rust)
    }
    fn as_str(&self) -> &'static str {
        ck::CipDescriptor::from(*self).as_str()
    }
}

#[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum CipDescriptorErrorKind {
    InvalidNeighborOrder,
    Property,
    InvalidStoredDescriptor,
}
#[pyclass(module="fixed_projection._native", extends=PyValueError)]
pub struct CipDescriptorError {
    inner: ck::CipDescriptorError,
}
#[pymethods]
impl CipDescriptorError {
    #[getter]
    fn domain(&self) -> &'static str {
        "cip"
    }
    #[getter]
    fn variant(&self) -> CipDescriptorErrorKind {
        match &self.inner {
            ck::CipDescriptorError::InvalidNeighborOrder {
                value: _,
                detail: _,
            } => CipDescriptorErrorKind::InvalidNeighborOrder,
            ck::CipDescriptorError::Property(_) => CipDescriptorErrorKind::Property,
            ck::CipDescriptorError::InvalidStoredDescriptor { value: _ } => {
                CipDescriptorErrorKind::InvalidStoredDescriptor
            }
        }
    }
    #[getter]
    fn value(&self, py: Python<'_>) -> PyResult<Option<String>> {
        match &self.inner {
            ck::CipDescriptorError::InvalidNeighborOrder { value, .. }
            | ck::CipDescriptorError::InvalidStoredDescriptor { value } => {
                decode_source_text(py, value).map(Some)
            }
            ck::CipDescriptorError::Property(_) => Ok(None),
        }
    }
    #[getter]
    fn detail(&self) -> Option<&str> {
        match &self.inner {
            ck::CipDescriptorError::InvalidNeighborOrder { detail, .. } => Some(detail),
            _ => None,
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn property_error(&self, py: Python<'_>) -> PyResult<Option<Py<PropertyValueError>>> {
        match &self.inner {
            ck::CipDescriptorError::Property(value) => PropertyValueError::from_rust(value.clone())
                .value(py)
                .extract::<Py<PropertyValueError>>()
                .map(Some)
                .map_err(Into::into),
            _ => Ok(None),
        }
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
impl CipDescriptorError {
    fn from_rust(inner: ck::CipDescriptorError) -> PyErr {
        Python::attach(|py| {
            // Native Property is #[error(transparent)]: its source is the
            // source-free PropertyValueError's source, not that payload itself.
            match Py::new(py, Self { inner }) {
                Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
                Err(error) => error,
            }
        })
    }
}

#[cfg(test)]
mod source_tests {
    use super::*;
    use crate::property_values::PropertyValueKind;
    use std::error::Error;

    #[test]
    fn real_cip_property_failure_keeps_payload_without_inventing_a_cause() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let spec = ck::AtomSpec::new(ck::Element::C)
                .with_prop("_CIPNeighborOrder", ck::PropertyValue::Bool(true))
                .unwrap();
            let atom = ck::Atom::from_spec(ck::AtomId::new(0), spec);
            let native = atom.cip_neighbor_order().unwrap_err();
            assert!(matches!(native, ck::CipDescriptorError::Property(_)));
            assert!(native.source().is_none());
            let expected_display = native.to_string();
            let atom = Py::new(py, Atom { inner: atom })?;
            let error = atom
                .bind(py)
                .call_method0("cip_neighbor_order")
                .unwrap_err();
            assert!(error.is_instance_of::<CipDescriptorError>(py));
            assert!(error.cause(py).is_none());
            assert_eq!(error.value(py).str()?.to_str()?, expected_display);
            let payload = error.value(py).getattr("property_error")?;
            assert!(payload.is_instance_of::<PropertyValueError>());
            assert!(
                payload
                    .call_method0("actual")?
                    .extract::<PropertyValueKind>()?
                    == PropertyValueKind::Bool
            );
            assert!(
                payload
                    .call_method0("expected")?
                    .extract::<PropertyValueKind>()?
                    == PropertyValueKind::IntVector
            );
            assert!(payload.getattr("__cause__")?.is_none());
            assert!(
                error
                    .value(py)
                    .setattr("property_error", py.None())
                    .unwrap_err()
                    .is_instance_of::<pyo3::exceptions::PyAttributeError>(py)
            );
            Ok(())
        })
        .unwrap();
    }
}

#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone)]
pub struct Atom {
    pub(crate) inner: ck::Atom,
}
#[pymethods]
impl Atom {
    fn id(&self) -> usize {
        self.inner.id().index()
    }
    fn element(&self) -> Element {
        Element {
            inner: self.inner.element(),
        }
    }
    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }
    fn formal_charge(&self) -> i8 {
        self.inner.formal_charge()
    }
    fn chiral_tag(&self) -> ChiralTag {
        self.inner.chiral_tag().into()
    }
    fn chiral_tag_code(&self) -> i64 {
        self.inner.chiral_tag_code()
    }
    fn chiral_tag_name(&self) -> &'static str {
        self.inner.chiral_tag_name()
    }
    fn isotope(&self) -> Option<u16> {
        self.inner.isotope()
    }
    fn atom_map(&self) -> Option<u32> {
        self.inner.atom_map()
    }
    fn is_aromatic(&self) -> bool {
        self.inner.is_aromatic()
    }
    fn explicit_hydrogens(&self) -> u8 {
        self.inner.explicit_hydrogens()
    }
    fn no_implicit(&self) -> bool {
        self.inner.no_implicit()
    }
    fn radical_electrons(&self) -> u8 {
        self.inner.radical_electrons()
    }
    fn hybridization(&self) -> Hybridization {
        self.inner.hybridization().into()
    }
    fn template_attachment_order(
        &self,
    ) -> Option<crate::attachment_values::TemplateAttachmentOrder> {
        self.inner
            .template_attachment_order()
            .cloned()
            .map(|inner| crate::attachment_values::TemplateAttachmentOrder { inner })
    }
    fn cip_descriptor(&self) -> PyResult<Option<CipDescriptor>> {
        self.inner
            .cip_descriptor()
            .map(|value| value.map(Into::into))
            .map_err(CipDescriptorError::from_rust)
    }
    fn cip_neighbor_order(&self) -> PyResult<Option<Vec<u32>>> {
        self.inner
            .cip_neighbor_order()
            .map_err(CipDescriptorError::from_rust)
    }
    fn cip_rank(&self) -> PyResult<Option<u32>> {
        self.inner.cip_rank().map_err(PropertyValueError::from_rust)
    }
}

#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone)]
pub struct Bond {
    pub(crate) inner: ck::Bond,
}
#[pymethods]
impl Bond {
    fn id(&self) -> usize {
        self.inner.id().index()
    }
    fn begin(&self) -> usize {
        self.inner.begin().index()
    }
    fn end(&self) -> usize {
        self.inner.end().index()
    }
    fn order(&self) -> BondOrder {
        self.inner.order().into()
    }
    fn order_code(&self) -> i64 {
        self.inner.order_code()
    }
    fn order_name(&self) -> &'static str {
        self.inner.order_name()
    }
    fn direction(&self) -> BondDirection {
        self.inner.direction().into()
    }
    fn direction_code(&self) -> i64 {
        self.inner.direction_code()
    }
    fn direction_name(&self) -> &'static str {
        self.inner.direction_name()
    }
    fn stereo(&self) -> BondStereo {
        self.inner.stereo().into()
    }
    fn stereo_code(&self) -> i64 {
        self.inner.stereo_code()
    }
    fn stereo_name(&self) -> &'static str {
        self.inner.stereo_name()
    }
    fn stereo_atoms(&self) -> Option<[usize; 2]> {
        self.inner
            .stereo_atoms()
            .map(|atoms| atoms.map(ck::AtomId::index))
    }
    fn is_aromatic(&self) -> bool {
        self.inner.is_aromatic()
    }
    fn cip_descriptor(&self) -> PyResult<Option<CipDescriptor>> {
        self.inner
            .cip_descriptor()
            .map(|value| value.map(Into::into))
            .map_err(CipDescriptorError::from_rust)
    }
    fn cip_neighbor_order(&self) -> PyResult<Option<Vec<u32>>> {
        self.inner
            .cip_neighbor_order()
            .map_err(CipDescriptorError::from_rust)
    }
}
