//! Concrete numeric property causes used by the notation writer.
use crate::error_framework::error_projection;
use crate::property_values::{PropertyValueKind, decode_source_text};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(UIntLexicalReadError, UIntLexicalReadErrorKind, ck::UIntLexicalReadError, "property", {
    Empty, Character { position: usize, byte: u8 }, Overflow,
}, cause = |_: &ck::UIntLexicalReadError| None::<PyErr>);

#[derive(IntoPyObject)]
enum PropertyNumberText {
    Number(i64),
    Text(String),
}

error_projection!(PropertyIntReadError, PropertyIntReadErrorKind, ck::PropertyIntReadError, "property", {
    UnsignedOverflow { value }, InvalidKind { kind },
    SignedTextOverflow { value, magnitude: u32, negative: bool },
    Lexical { value, source },
}, cause = |inner: &ck::PropertyIntReadError| match inner {
    ck::PropertyIntReadError::Lexical { source, .. } => Some(UIntLexicalReadError::from_rust(source.clone())),
    ck::PropertyIntReadError::UnsignedOverflow { .. } | ck::PropertyIntReadError::InvalidKind { .. }
    | ck::PropertyIntReadError::SignedTextOverflow { .. } => None,
}, source = Option<crate::error_framework::ErrorObject<UIntLexicalReadError>>, methods = {
    #[getter]
    #[doc = "Value rejected or reported by the native operation. Applies to PropertyIntReadError.UnsignedOverflow, PropertyIntReadError.SignedTextOverflow, PropertyIntReadError.Lexical, PropertyIntReadError.InvalidKind. None when this payload is absent."] fn value(&self, py: Python<'_>) -> PyResult<Option<PropertyNumberText>> {
        match &self.inner {
            ck::PropertyIntReadError::UnsignedOverflow { value } => Ok(Some(PropertyNumberText::Number(i64::from(*value)))),
            ck::PropertyIntReadError::SignedTextOverflow { value, .. }
            | ck::PropertyIntReadError::Lexical { value, .. } => decode_source_text(py, value).map(|value| Some(PropertyNumberText::Text(value))),
            ck::PropertyIntReadError::InvalidKind { .. } => Ok(None),
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native failure-specific category or value kind. Applies to PropertyIntReadError.InvalidKind. None when this payload is absent."] fn kind(&self) -> Option<PropertyValueKind> {
        match self.inner { ck::PropertyIntReadError::InvalidKind { kind } => Some(kind.into()), _ => None }
    }
});

error_projection!(PropertyUIntReadError, PropertyUIntReadErrorKind, ck::PropertyUIntReadError, "property", {
    Negative { value }, InvalidKind { kind }, Lexical { value, source },
}, cause = |inner: &ck::PropertyUIntReadError| match inner {
    ck::PropertyUIntReadError::Lexical { source, .. } => Some(UIntLexicalReadError::from_rust(source.clone())),
    ck::PropertyUIntReadError::Negative { .. } | ck::PropertyUIntReadError::InvalidKind { .. } => None,
}, source = Option<crate::error_framework::ErrorObject<UIntLexicalReadError>>, methods = {
    #[getter]
    #[doc = "Value rejected or reported by the native operation. Applies to PropertyUIntReadError.Negative, PropertyUIntReadError.Lexical, PropertyUIntReadError.InvalidKind. None when this payload is absent."] fn value(&self, py: Python<'_>) -> PyResult<Option<PropertyNumberText>> {
        match &self.inner {
            ck::PropertyUIntReadError::Negative { value } => Ok(Some(PropertyNumberText::Number(i64::from(*value)))),
            ck::PropertyUIntReadError::Lexical { value, .. } => decode_source_text(py, value).map(|value| Some(PropertyNumberText::Text(value))),
            ck::PropertyUIntReadError::InvalidKind { .. } => Ok(None),
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native failure-specific category or value kind. Applies to PropertyUIntReadError.InvalidKind. None when this payload is absent."] fn kind(&self) -> Option<PropertyValueKind> {
        match self.inner { ck::PropertyUIntReadError::InvalidKind { kind } => Some(kind.into()), _ => None }
    }
});

error_projection!(PropertyULongReadError, PropertyULongReadErrorKind, ck::PropertyULongReadError, "property", {
    Negative { value }, InvalidKind { kind }, LexicalOverflow { value }, Lexical { value, source },
}, cause = |inner: &ck::PropertyULongReadError| match inner {
    ck::PropertyULongReadError::Lexical { source, .. } => Some(UIntLexicalReadError::from_rust(source.clone())),
    ck::PropertyULongReadError::Negative { .. } | ck::PropertyULongReadError::InvalidKind { .. }
    | ck::PropertyULongReadError::LexicalOverflow { .. } => None,
}, source = Option<crate::error_framework::ErrorObject<UIntLexicalReadError>>, methods = {
    #[getter]
    #[doc = "Value rejected or reported by the native operation. Applies to PropertyULongReadError.Negative, PropertyULongReadError.Lexical, PropertyULongReadError.LexicalOverflow, PropertyULongReadError.InvalidKind. None when this payload is absent."] fn value(&self, py: Python<'_>) -> PyResult<Option<PropertyNumberText>> {
        match &self.inner {
            ck::PropertyULongReadError::Negative { value } => Ok(Some(PropertyNumberText::Number(i64::from(*value)))),
            ck::PropertyULongReadError::Lexical { value, .. }
            | ck::PropertyULongReadError::LexicalOverflow { value } => decode_source_text(py, value).map(|value| Some(PropertyNumberText::Text(value))),
            ck::PropertyULongReadError::InvalidKind { .. } => Ok(None),
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native failure-specific category or value kind. Applies to PropertyULongReadError.InvalidKind. None when this payload is absent."] fn kind(&self) -> Option<PropertyValueKind> {
        match self.inner { ck::PropertyULongReadError::InvalidKind { kind } => Some(kind.into()), _ => None }
    }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn real_writer_numeric_failure_preserves_kind_payload_and_typed_lexical_source() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let original = ck::Molecule::from_smiles("[Pt@SP1](Cl)(Br)(F)I").unwrap();
            // This private-source property is deliberately prepared in a
            // detached fixture, never through the user-property operation.
            let mut topology = original.topology().clone();
            topology.atoms[0]
                .set_prop("_chiralPermutation", "bad")
                .unwrap();
            let molecule = ck::Molecule::from_parts(
                topology,
                ck::CoordinateBlock::default(),
                original.properties().clone(),
            )
            .unwrap();
            let source = molecule.to_smiles().unwrap_err();
            let mut next: Option<&(dyn Error + 'static)> = Some(&source);
            let native = loop {
                let error = next.expect("writer must retain its numeric cause");
                if let Some(error) = error.downcast_ref::<ck::PropertyIntReadError>() {
                    break error.clone();
                }
                next = error.source();
            };
            let display = native.to_string();
            let error = PropertyIntReadError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(
                error.value(py).getattr("value")?.extract::<String>()?,
                "bad"
            );
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<PropertyIntReadErrorKind>()?
                    == PropertyIntReadErrorKind::Lexical
            );
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            assert!(
                cause
                    .extract::<crate::error_framework::ErrorObject<UIntLexicalReadError>>()?
                    .borrow(py)
                    .variant()
                    == UIntLexicalReadErrorKind::Character
            );
            assert_eq!(cause.getattr("position")?.extract::<usize>()?, 0);
            assert_eq!(cause.getattr("byte")?.extract::<u8>()?, b'b');
            assert!(cause.getattr("__cause__")?.is_none());
            let error = PropertyUIntReadError::from_rust(ck::PropertyUIntReadError::InvalidKind {
                kind: ck::PropertyValueKind::Bool,
            });
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
                    .extract::<PropertyUIntReadErrorKind>()?
                    == PropertyUIntReadErrorKind::InvalidKind
            );
            assert!(
                error
                    .value(py)
                    .setattr("kind", PropertyValueKind::Int)
                    .is_err()
            );
            Ok(())
        })
        .unwrap();
    }
}
