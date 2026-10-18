//! Exact property/canonical failures needed by tautomer scoring and notation.
use crate::attachment_values::TemplateAttachmentOrderError;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::MoleculePropertyError;
use crate::model_errors::{StereoGroupError, TopologyValidationError};
use crate::property_values::{
    MissingPropertyError, PropertyValue, PropertyValueKind, decode_source_text,
};
use crate::ring_errors::RingFindingError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use std::error::Error;

#[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum PropertyStringErrorKind {
    UnsupportedKind,
}
#[pyclass(module="fixed_projection._native", extends=PyValueError)]
pub struct PropertyStringError {
    inner: ck::PropertyStringError,
}
#[pymethods]
impl PropertyStringError {
    #[getter]
    fn domain(&self) -> &'static str {
        "property"
    }
    #[getter]
    fn variant(&self) -> PyResult<PropertyStringErrorKind> {
        match self.inner {
            ck::PropertyStringError::UnsupportedKind { kind: _ } => {
                Ok(PropertyStringErrorKind::UnsupportedKind)
            }
            _ => Err(PyTypeError::new_err(
                "unrecognized native PropertyStringError variant",
            )),
        }
    }
    /// Preserve the actual canonical Rust method; this is not a discriminator.
    fn kind(&self) -> PropertyValueKind {
        self.inner.kind().into()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
impl PropertyStringError {
    pub(crate) fn from_rust(inner: ck::PropertyStringError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}

#[pyfunction]
pub fn property_value_to_text(py: Python<'_>, value: &PropertyValue) -> PyResult<String> {
    let native =
        ck::property_value_to_text(&value.inner).map_err(PropertyStringError::from_rust)?;
    decode_source_text(py, &native)
}

#[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum RequiredPropertyStringErrorKind {
    Missing,
    Conversion,
}
#[pyclass(module="fixed_projection._native", extends=PyValueError)]
pub struct RequiredPropertyStringError {
    inner: ck::RequiredPropertyStringError,
}
#[pymethods]
impl RequiredPropertyStringError {
    #[getter]
    fn domain(&self) -> &'static str {
        "property"
    }
    #[getter]
    fn variant(&self) -> PyResult<RequiredPropertyStringErrorKind> {
        match self.inner {
            ck::RequiredPropertyStringError::Missing(_) => {
                Ok(RequiredPropertyStringErrorKind::Missing)
            }
            ck::RequiredPropertyStringError::Conversion(_) => {
                Ok(RequiredPropertyStringErrorKind::Conversion)
            }
            _ => Err(PyTypeError::new_err(
                "unrecognized native RequiredPropertyStringError variant",
            )),
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn missing(&self, py: Python<'_>) -> PyResult<Option<Py<MissingPropertyError>>> {
        match &self.inner {
            ck::RequiredPropertyStringError::Missing(value) => {
                MissingPropertyError::from_rust(value.clone())
                    .value(py)
                    .extract()
                    .map(Some)
                    .map_err(Into::into)
            }
            _ => Ok(None),
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn conversion(&self, py: Python<'_>) -> PyResult<Option<Py<PropertyStringError>>> {
        match &self.inner {
            ck::RequiredPropertyStringError::Conversion(value) => {
                PropertyStringError::from_rust(*value)
                    .value(py)
                    .extract()
                    .map(Some)
                    .map_err(Into::into)
            }
            _ => Ok(None),
        }
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
impl RequiredPropertyStringError {
    pub(crate) fn from_rust(inner: ck::RequiredPropertyStringError) -> PyErr {
        // Both current variants are transparent, and their values have no
        // native source. Reject future non-exhaustive variants explicitly.
        match &inner {
            ck::RequiredPropertyStringError::Missing(_)
            | ck::RequiredPropertyStringError::Conversion(_) => {}
            _ => {
                return PyTypeError::new_err(
                    "unrecognized native RequiredPropertyStringError variant",
                );
            }
        }
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}

type CanonicalSource = Source<
    Py<StereoGroupError>,
    Source<
        Py<TemplateAttachmentOrderError>,
        Source<Py<MoleculePropertyError>, Py<TopologyValidationError>>,
    >,
>;

error_projection!(CanonicalRankError, CanonicalRankErrorKind, ck::CanonicalRankError, "canonical_rank", {
    HanoiScratchTooSmall, StereoGroup(_),
    UnsignedRankOverflow { atom_index, property, value: u32 },
    InvalidPropertyKind { atom_index, property, kind },
    ComponentAtomOutOfRange { atom_index, atom_count }, DuplicateComponentAtom { atom_index },
    AtomMaskLength { expected, actual }, BondMaskLength { expected, actual },
    AtomSymbolLength { expected, actual }, BondSymbolLength { expected, actual },
    PreparedValenceLength { atom_count, explicit_len: usize, implicit_len: usize },
    PreparedValenceInvalid { atom_index },
    PreparedRingLength { expected_atoms: usize, actual_atoms: usize, expected_bonds: usize, actual_bonds: usize },
    ProtocolDebt { branch: &'static str, reason: &'static str },
    TemplateAttachmentRemap { carrier: usize => |value: &ck::AtomId| value.index(), source },
    RingFinding(_), Valence(_),
}, cause = |inner: &ck::CanonicalRankError| Python::attach(|py| match inner {
    ck::CanonicalRankError::StereoGroup(value) => Some(StereoGroupError::from_rust(*value)),
    ck::CanonicalRankError::TemplateAttachmentRemap { source, .. } => Some(TemplateAttachmentOrderError::from_rust(source.clone())),
    ck::CanonicalRankError::RingFinding(value) => RingFindingError::from_rust(value.clone()).cause(py),
    ck::CanonicalRankError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::CanonicalRankError::HanoiScratchTooSmall | ck::CanonicalRankError::UnsignedRankOverflow { .. }
    | ck::CanonicalRankError::InvalidPropertyKind { .. } | ck::CanonicalRankError::ComponentAtomOutOfRange { .. }
    | ck::CanonicalRankError::DuplicateComponentAtom { .. } | ck::CanonicalRankError::AtomMaskLength { .. }
    | ck::CanonicalRankError::BondMaskLength { .. } | ck::CanonicalRankError::AtomSymbolLength { .. }
    | ck::CanonicalRankError::BondSymbolLength { .. } | ck::CanonicalRankError::PreparedValenceLength { .. }
    | ck::CanonicalRankError::PreparedValenceInvalid { .. } | ck::CanonicalRankError::PreparedRingLength { .. }
    | ck::CanonicalRankError::ProtocolDebt { .. } => None,
}), source = Option<CanonicalSource>, methods = {
    #[getter]
    fn atom_index(&self) -> Option<usize> { match self.inner {
        ck::CanonicalRankError::UnsignedRankOverflow { atom_index, .. }
        | ck::CanonicalRankError::InvalidPropertyKind { atom_index, .. }
        | ck::CanonicalRankError::ComponentAtomOutOfRange { atom_index, .. }
        | ck::CanonicalRankError::DuplicateComponentAtom { atom_index }
        | ck::CanonicalRankError::PreparedValenceInvalid { atom_index } => Some(atom_index), _ => None,
    } }
    #[getter]
    fn property(&self) -> Option<&'static str> { match self.inner {
        ck::CanonicalRankError::UnsignedRankOverflow { property, .. }
        | ck::CanonicalRankError::InvalidPropertyKind { property, .. } => Some(property), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kind(&self) -> Option<PropertyValueKind> { match self.inner {
        ck::CanonicalRankError::InvalidPropertyKind { kind, .. } => Some(kind.into()), _ => None,
    } }
    #[getter]
    fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::CanonicalRankError::ComponentAtomOutOfRange { atom_count, .. }
        | ck::CanonicalRankError::PreparedValenceLength { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    fn expected(&self) -> Option<usize> { match self.inner {
        ck::CanonicalRankError::AtomMaskLength { expected, .. } | ck::CanonicalRankError::BondMaskLength { expected, .. }
        | ck::CanonicalRankError::AtomSymbolLength { expected, .. } | ck::CanonicalRankError::BondSymbolLength { expected, .. } => Some(expected), _ => None,
    } }
    #[getter]
    fn actual(&self) -> Option<usize> { match self.inner {
        ck::CanonicalRankError::AtomMaskLength { actual, .. } | ck::CanonicalRankError::BondMaskLength { actual, .. }
        | ck::CanonicalRankError::AtomSymbolLength { actual, .. } | ck::CanonicalRankError::BondSymbolLength { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn rings(&self, py: Python<'_>) -> PyResult<Option<Py<RingFindingError>>> {
        match &self.inner { ck::CanonicalRankError::RingFinding(value) =>
            RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<Py<ValenceError>>> {
        match &self.inner { ck::CanonicalRankError::Valence(value) =>
            ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
});

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn native_property_conversion_and_canonical_payload_keep_kind_separate_from_variant() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            for native in [
                ck::PropertyValue::Int(-7),
                ck::PropertyValue::UInt(9),
                ck::PropertyValue::Bool(true),
            ] {
                let expected = ck::property_value_to_text(&native).unwrap();
                let value = PropertyValue { inner: native };
                assert_eq!(
                    property_value_to_text(py, &value)?,
                    decode_source_text(py, &expected)?
                );
            }
            let native = ck::PropertyStringError::UnsupportedKind {
                kind: ck::PropertyValueKind::Bool,
            };
            let error = PropertyStringError::from_rust(native);
            assert!(
                error
                    .value(py)
                    .call_method0("kind")?
                    .extract::<PropertyValueKind>()?
                    == PropertyValueKind::Bool
            );
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<PropertyStringErrorKind>()?
                    == PropertyStringErrorKind::UnsupportedKind
            );
            let native = ck::CanonicalRankError::InvalidPropertyKind {
                atom_index: 4,
                property: "_rank",
                kind: ck::PropertyValueKind::Bool,
            };
            let display = native.to_string();
            let error = CanonicalRankError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(
                error.value(py).getattr("atom_index")?.extract::<usize>()?,
                4
            );
            assert!(
                error
                    .value(py)
                    .getattr("kind")?
                    .extract::<PropertyValueKind>()?
                    == PropertyValueKind::Bool
            );
            assert!(error.cause(py).is_none());
            Ok(())
        })
        .unwrap();
    }
}
