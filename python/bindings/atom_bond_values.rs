use crate::chemical_values::{
    BondDirection, BondOrder, BondStereo, ChiralTag, Hybridization, enum_value,
};
use crate::elements::Element;
use crate::property_values::{PropertyValueError, decode_source_text};
use cosmolkit as ck;
use pyo3::prelude::*;

enum_value!(CipDescriptor, { R, S, LowerR [r], LowerS [s], E, Z, LowerE [e], LowerZ [z], M, P, LowerM [m], LowerP [p] }, constructor = custom, methods = {
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
});

crate::enums::binding_enum!(CipDescriptorErrorKind, {
    InvalidNeighborOrder,
    Property,
    InvalidStoredDescriptor,
});
crate::error_framework::error_projection! {record
pub struct CipDescriptorError {
    inner: ck::CipDescriptorError,
}
}
crate::error_framework::error_projection! {methods impl CipDescriptorError {
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
    fn property_error(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PropertyValueError>>> {
        match &self.inner {
            ck::CipDescriptorError::Property(value) => PropertyValueError::from_rust(value.clone())
                .value(py)
                .extract::<crate::error_framework::ErrorObject<PropertyValueError>>()
                .map(Some)
                .map_err(Into::into),
            _ => Ok(None),
        }
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
impl CipDescriptorError {
    fn from_rust(inner: ck::CipDescriptorError) -> PyErr {
        Python::attach(|py| {
            // Native Property is #[error(transparent)]: its source is the
            // source-free PropertyValueError's source, not that payload itself.
            match Py::new(py, Self { inner }) {
                Ok(value) => crate::error_framework::from_payload(value),
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
            assert!(error.is_instance(py, <CipDescriptorError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert!(error.cause(py).is_none());
            assert_eq!(error.value(py).str()?.to_str()?, expected_display);
            let payload = error.value(py).getattr("property_error")?;
            assert!(payload.is_instance(<PropertyValueError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(
                payload
                    .getattr("actual")?
                    .extract::<PropertyValueKind>()?
                    == PropertyValueKind::Bool
            );
            assert!(
                payload
                    .getattr("expected")?
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

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct Atom {
        pub(crate) inner: ck::Atom,
    }
}
macro_rules! atom_accessors {
    ($type:ident $(, $extra:item)*) => {
#[pymethods]
impl $type {
    #[getter]
    fn id(&self) -> usize {
        self.inner.id().index()
    }
    #[getter]
    fn element(&self) -> Element {
        Element {
            inner: self.inner.element(),
        }
    }
    #[getter]
    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }
    #[getter]
    fn formal_charge(&self) -> i8 {
        self.inner.formal_charge()
    }
    #[getter]
    fn chiral_tag(&self) -> ChiralTag {
        self.inner.chiral_tag().into()
    }
    #[getter]
    fn chiral_tag_code(&self) -> i64 {
        self.inner.chiral_tag_code()
    }
    #[getter]
    fn chiral_tag_name(&self) -> &'static str {
        self.inner.chiral_tag_name()
    }
    #[getter]
    fn isotope(&self) -> Option<u16> {
        self.inner.isotope()
    }
    #[getter]
    fn atom_map(&self) -> Option<u32> {
        self.inner.atom_map()
    }
    #[getter]
    fn is_aromatic(&self) -> bool {
        self.inner.is_aromatic()
    }
    #[getter]
    fn explicit_hydrogens(&self) -> u8 {
        self.inner.explicit_hydrogens()
    }
    #[getter]
    fn no_implicit(&self) -> bool {
        self.inner.no_implicit()
    }
    #[getter]
    fn radical_electrons(&self) -> u8 {
        self.inner.radical_electrons()
    }
    #[getter]
    fn hybridization(&self) -> Hybridization {
        self.inner.hybridization().into()
    }
    #[getter]
    fn template_attachment_order(
        &self,
    ) -> Option<crate::attachment_values::TemplateAttachmentOrder> {
        self.inner
            .template_attachment_order()
            .cloned()
            .map(|inner| crate::attachment_values::TemplateAttachmentOrder { inner })
    }
    #[getter]
    fn cip_descriptor(&self) -> PyResult<Option<CipDescriptor>> {
        self.inner
            .cip_descriptor()
            .map(|value| value.map(Into::into))
            .map_err(CipDescriptorError::from_rust)
    }
    #[getter]
    fn cip_neighbor_order(&self) -> PyResult<Option<Vec<u32>>> {
        self.inner
            .cip_neighbor_order()
            .map_err(CipDescriptorError::from_rust)
    }
    #[getter]
    fn cip_rank(&self) -> PyResult<Option<u32>> {
        self.inner.cip_rank().map_err(PropertyValueError::from_rust)
    }
    $($extra)*
}
    };
}
atom_accessors!(Atom);

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct MoleculeAtom {
        pub(crate) inner: ck::MoleculeAtom,
    }
}
atom_accessors!(
    MoleculeAtom,
    fn value(&self) -> Atom {
        Atom {
            inner: self.inner.value().clone(),
        }
    },
    #[getter]
    fn degree(&self) -> usize {
        self.inner.degree()
    },
    #[getter]
    fn explicit_valence(&self) -> PyResult<i32> {
        self.inner
            .explicit_valence()
            .map_err(crate::valence_errors::ValenceError::from_rust)
    },
    #[getter]
    fn implicit_hydrogens(&self) -> PyResult<i32> {
        self.inner
            .implicit_hydrogens()
            .map_err(crate::valence_errors::ValenceError::from_rust)
    },
    #[getter]
    fn total_hydrogens(&self) -> PyResult<i32> {
        self.inner
            .total_hydrogens()
            .map_err(crate::valence_errors::ValenceError::from_rust)
    },
    #[getter]
    fn total_valence(&self) -> PyResult<i32> {
        self.inner
            .total_valence()
            .map_err(crate::valence_errors::ValenceError::from_rust)
    }
);

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct Bond {
        pub(crate) inner: ck::Bond,
    }
}
#[pymethods]
impl Bond {
    #[getter]
    fn id(&self) -> usize {
        self.inner.id().index()
    }
    #[getter]
    fn begin(&self) -> usize {
        self.inner.begin().index()
    }
    #[getter]
    fn end(&self) -> usize {
        self.inner.end().index()
    }
    #[getter]
    fn order(&self) -> BondOrder {
        self.inner.order().into()
    }
    #[getter]
    fn order_code(&self) -> i64 {
        self.inner.order_code()
    }
    #[getter]
    fn order_name(&self) -> &'static str {
        self.inner.order_name()
    }
    #[getter]
    fn direction(&self) -> BondDirection {
        self.inner.direction().into()
    }
    #[getter]
    fn direction_code(&self) -> i64 {
        self.inner.direction_code()
    }
    #[getter]
    fn direction_name(&self) -> &'static str {
        self.inner.direction_name()
    }
    #[getter]
    fn stereo(&self) -> BondStereo {
        self.inner.stereo().into()
    }
    #[getter]
    fn stereo_code(&self) -> i64 {
        self.inner.stereo_code()
    }
    #[getter]
    fn stereo_name(&self) -> &'static str {
        self.inner.stereo_name()
    }
    #[getter]
    fn stereo_atoms(&self) -> Option<[usize; 2]> {
        self.inner
            .stereo_atoms()
            .map(|atoms| atoms.map(ck::AtomId::index))
    }
    #[getter]
    fn is_aromatic(&self) -> bool {
        self.inner.is_aromatic()
    }
    #[getter]
    fn cip_descriptor(&self) -> PyResult<Option<CipDescriptor>> {
        self.inner
            .cip_descriptor()
            .map(|value| value.map(Into::into))
            .map_err(CipDescriptorError::from_rust)
    }
    #[getter]
    fn cip_neighbor_order(&self) -> PyResult<Option<Vec<u32>>> {
        self.inner
            .cip_neighbor_order()
            .map_err(CipDescriptorError::from_rust)
    }
}

crate::configuration::type_projection! {
    impl NativeScalar for cosmolkit::AtomId {
        type Repr = usize;
        fn from_scalar(index: usize) -> Self {
            Self::new(index)
        }
        fn scalar(self) -> usize {
            self.index()
        }
    }
}

crate::configuration::type_projection! {
    impl NativeScalar for cosmolkit::BondId {
        type Repr = usize;
        fn from_scalar(index: usize) -> Self {
            Self::new(index)
        }
        fn scalar(self) -> usize {
            self.index()
        }
    }
}
