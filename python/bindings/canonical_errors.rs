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
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use std::error::Error;

crate::enums::binding_enum!(PropertyStringErrorKind, {
    UnsupportedKind,
});
crate::error_framework::error_projection! {record
pub struct PropertyStringError {
    inner: ck::PropertyStringError,
}
}
crate::error_framework::error_projection! {methods impl PropertyStringError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "property"
    }
    #[getter]
    /// Exact native error category. Variant-specific payload fields are None when that variant is not active.
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
    #[getter]
    fn kind(&self) -> PropertyValueKind {
        self.inner.kind().into()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
impl PropertyStringError {
    pub(crate) fn from_rust(inner: ck::PropertyStringError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => crate::error_framework::from_payload(value),
            Err(error) => error,
        })
    }
}

#[pyfunction]
/// Convert a typed PropertyValue to text using the canonical native conversion; conversion errors propagate.
pub fn property_value_to_text(py: Python<'_>, value: &PropertyValue) -> PyResult<String> {
    let native =
        ck::property_value_to_text(&value.inner).map_err(PropertyStringError::from_rust)?;
    decode_source_text(py, &native)
}

crate::enums::binding_enum!(RequiredPropertyStringErrorKind, {
    Missing,
    Conversion,
});
crate::error_framework::error_projection! {record
pub struct RequiredPropertyStringError {
    inner: ck::RequiredPropertyStringError,
}
}
crate::error_framework::error_projection! {methods impl RequiredPropertyStringError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "property"
    }
    #[getter]
    /// Exact native error category. Variant-specific payload fields are None when that variant is not active.
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
    #[doc = "Structured MissingPropertyError payload of RequiredPropertyStringError.Missing; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn missing(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<MissingPropertyError>>> {
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
    #[doc = "Structured PropertyStringError payload of RequiredPropertyStringError.Conversion; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn conversion(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PropertyStringError>>> {
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
            Ok(value) => crate::error_framework::from_payload(value),
            Err(error) => error,
        })
    }
}

type CanonicalSource = Source<
    crate::error_framework::ErrorObject<StereoGroupError>,
    Source<
        crate::error_framework::ErrorObject<TemplateAttachmentOrderError>,
        Source<
            crate::error_framework::ErrorObject<MoleculePropertyError>,
            crate::error_framework::ErrorObject<TopologyValidationError>,
        >,
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
    #[doc = "Atom index involved in the native failure. Applies to CanonicalRankError.UnsignedRankOverflow, CanonicalRankError.InvalidPropertyKind, CanonicalRankError.ComponentAtomOutOfRange, CanonicalRankError.DuplicateComponentAtom, CanonicalRankError.PreparedValenceInvalid. None when this payload is absent."] fn atom_index(&self) -> Option<usize> { match self.inner {
        ck::CanonicalRankError::UnsignedRankOverflow { atom_index, .. }
        | ck::CanonicalRankError::InvalidPropertyKind { atom_index, .. }
        | ck::CanonicalRankError::ComponentAtomOutOfRange { atom_index, .. }
        | ck::CanonicalRankError::DuplicateComponentAtom { atom_index }
        | ck::CanonicalRankError::PreparedValenceInvalid { atom_index } => Some(atom_index), _ => None,
    } }
    #[getter]
    #[doc = "Property name involved in the failure. Applies to CanonicalRankError.UnsignedRankOverflow, CanonicalRankError.InvalidPropertyKind. None when this payload is absent."] fn property(&self) -> Option<&'static str> { match self.inner {
        ck::CanonicalRankError::UnsignedRankOverflow { property, .. }
        | ck::CanonicalRankError::InvalidPropertyKind { property, .. } => Some(property), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native failure-specific category or value kind. Applies to CanonicalRankError.InvalidPropertyKind. None when this payload is absent."] fn kind(&self) -> Option<PropertyValueKind> { match self.inner {
        ck::CanonicalRankError::InvalidPropertyKind { kind, .. } => Some(kind.into()), _ => None,
    } }
    #[getter]
    #[doc = "Atom count used when validating the reported atom reference. Applies to CanonicalRankError.ComponentAtomOutOfRange, CanonicalRankError.PreparedValenceLength. None when this payload is absent."] fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::CanonicalRankError::ComponentAtomOutOfRange { atom_count, .. }
        | ck::CanonicalRankError::PreparedValenceLength { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    #[doc = "Required value reported by the native check. Applies to CanonicalRankError.AtomMaskLength, CanonicalRankError.BondMaskLength, CanonicalRankError.AtomSymbolLength, CanonicalRankError.BondSymbolLength. None when this payload is absent."] fn expected(&self) -> Option<usize> { match self.inner {
        ck::CanonicalRankError::AtomMaskLength { expected, .. } | ck::CanonicalRankError::BondMaskLength { expected, .. }
        | ck::CanonicalRankError::AtomSymbolLength { expected, .. } | ck::CanonicalRankError::BondSymbolLength { expected, .. } => Some(expected), _ => None,
    } }
    #[getter]
    #[doc = "Observed value reported by the native check. Applies to CanonicalRankError.AtomMaskLength, CanonicalRankError.BondMaskLength, CanonicalRankError.AtomSymbolLength, CanonicalRankError.BondSymbolLength. None when this payload is absent."] fn actual(&self) -> Option<usize> { match self.inner {
        ck::CanonicalRankError::AtomMaskLength { actual, .. } | ck::CanonicalRankError::BondMaskLength { actual, .. }
        | ck::CanonicalRankError::AtomSymbolLength { actual, .. } | ck::CanonicalRankError::BondSymbolLength { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured RingFindingError payload of CanonicalRankError.RingFinding; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn rings(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<RingFindingError>>> {
        match &self.inner { ck::CanonicalRankError::RingFinding(value) =>
            RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured ValenceError payload of CanonicalRankError.Valence; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn valence(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>> {
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
                    .getattr("kind")?
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
