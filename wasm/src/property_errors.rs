//! Property conversion errors preserve their native numeric and text payloads.
use crate::{
    errors::{enum_error, error_projection, exception},
    model_errors::PropertyValueKind,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};

pub(crate) fn text(value: &ck::PropertyText) -> Result<String, JsValue> {
    std::str::from_utf8(value.as_bytes())
        .map(str::to_owned)
        .map_err(|error| js_sys::TypeError::new(&error.to_string()).into())
}

// Native counted text has a UTF-8 JS projection. Reuse the existing readonly
// string-map shape; validate output text rather than replacing invalid bytes.
impl crate::fields::Field for std::collections::BTreeMap<ck::PropertyText, ck::PropertyText> {
    type Input = <std::collections::BTreeMap<String, String> as crate::fields::Field>::Input;
    type Output = <std::collections::BTreeMap<String, String> as crate::fields::Field>::Output;
    fn input(value: Self::Input) -> Result<Self, JsValue> {
        Ok(std::collections::BTreeMap::<String, String>::input(value)?
            .into_iter()
            .map(|(k, v)| (k.into(), v.into()))
            .collect())
    }
    fn decode(value: &JsValue) -> Result<Self, JsValue> {
        Ok(std::collections::BTreeMap::<String, String>::decode(value)?
            .into_iter()
            .map(|(k, v)| (k.into(), v.into()))
            .collect())
    }
    fn get(value: crate::storage::NativeField<Self>) -> Result<Self::Output, JsValue> {
        let strings = value.read(|v| {
            v.iter()
                .map(|(k, v)| Ok((text(k)?, text(v)?)))
                .collect::<Result<std::collections::BTreeMap<String, String>, JsValue>>()
        })??;
        <std::collections::BTreeMap<String, String> as crate::fields::Field>::get(
            crate::storage::NativeField::owned(strings),
        )
    }
}

enum_error!(UIntLexicalReadError, UIntLexicalReadErrorKind, ck::UIntLexicalReadError, {
    Empty => ck::UIntLexicalReadError::Empty,
    Character => ck::UIntLexicalReadError::Character { .. },
    Overflow => ck::UIntLexicalReadError::Overflow,
}, {
    #[wasm_bindgen(getter)] pub fn position(&self) -> Option<usize> { match *self.inner { ck::UIntLexicalReadError::Character { position, .. } => Some(position), _ => None } }
    #[wasm_bindgen(getter)] pub fn byte(&self) -> Option<u8> { match *self.inner { ck::UIntLexicalReadError::Character { byte, .. } => Some(byte), _ => None } }
});

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "number | string")]
    pub type PropertyNumberText;
    #[wasm_bindgen(typescript_type = "UIntLexicalReadError")]
    #[derive(Clone)]
    pub type LexicalCause;
    #[wasm_bindgen(typescript_type = "MissingPropertyError")]
    pub type MissingPayload;
    #[wasm_bindgen(typescript_type = "PropertyStringError")]
    pub type StringPayload;
}

enum_error!(PropertyIntReadError, PropertyIntReadErrorKind, ck::PropertyIntReadError, {
    UnsignedOverflow => ck::PropertyIntReadError::UnsignedOverflow { .. },
    InvalidKind => ck::PropertyIntReadError::InvalidKind { .. },
    SignedTextOverflow => ck::PropertyIntReadError::SignedTextOverflow { .. },
    Lexical => ck::PropertyIntReadError::Lexical { .. },
}, {
    #[wasm_bindgen(getter)] pub fn value(&self) -> Result<Option<PropertyNumberText>, JsValue> { Ok(match self.inner.as_ref() {
        ck::PropertyIntReadError::UnsignedOverflow { value } => Some(JsValue::from(*value).unchecked_into()),
        ck::PropertyIntReadError::SignedTextOverflow { value, .. } | ck::PropertyIntReadError::Lexical { value, .. } => Some(JsValue::from_str(&text(value)?).unchecked_into()),
        ck::PropertyIntReadError::InvalidKind { .. } => None,
    }) }
    #[wasm_bindgen(getter)] pub fn kind(&self) -> Option<PropertyValueKind> { match *self.inner { ck::PropertyIntReadError::InvalidKind { kind } => Some(kind.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn magnitude(&self) -> Option<u32> { match *self.inner { ck::PropertyIntReadError::SignedTextOverflow { magnitude, .. } => Some(magnitude), _ => None } }
    #[wasm_bindgen(getter)] pub fn negative(&self) -> Option<bool> { match *self.inner { ck::PropertyIntReadError::SignedTextOverflow { negative, .. } => Some(negative), _ => None } }
}, source = LexicalCause, cause = |value: &ck::PropertyIntReadError| match value {
    ck::PropertyIntReadError::Lexical { source, .. } => Some(UIntLexicalReadError::from_rust(source.clone()).unchecked_into()),
    ck::PropertyIntReadError::UnsignedOverflow { .. } | ck::PropertyIntReadError::InvalidKind { .. } | ck::PropertyIntReadError::SignedTextOverflow { .. } => None,
});

enum_error!(PropertyUIntReadError, PropertyUIntReadErrorKind, ck::PropertyUIntReadError, {
    Negative => ck::PropertyUIntReadError::Negative { .. }, InvalidKind => ck::PropertyUIntReadError::InvalidKind { .. }, Lexical => ck::PropertyUIntReadError::Lexical { .. },
}, {
    #[wasm_bindgen(getter)] pub fn value(&self) -> Result<Option<PropertyNumberText>, JsValue> { Ok(match self.inner.as_ref() {
        ck::PropertyUIntReadError::Negative { value } => Some(JsValue::from(*value).unchecked_into()),
        ck::PropertyUIntReadError::Lexical { value, .. } => Some(JsValue::from_str(&text(value)?).unchecked_into()),
        ck::PropertyUIntReadError::InvalidKind { .. } => None,
    }) }
    #[wasm_bindgen(getter)] pub fn kind(&self) -> Option<PropertyValueKind> { match *self.inner { ck::PropertyUIntReadError::InvalidKind { kind } => Some(kind.into()), _ => None } }
}, source = LexicalCause, cause = |value: &ck::PropertyUIntReadError| match value {
    ck::PropertyUIntReadError::Lexical { source, .. } => Some(UIntLexicalReadError::from_rust(source.clone()).unchecked_into()),
    ck::PropertyUIntReadError::Negative { .. } | ck::PropertyUIntReadError::InvalidKind { .. } => None,
});

enum_error!(PropertyULongReadError, PropertyULongReadErrorKind, ck::PropertyULongReadError, {
    Negative => ck::PropertyULongReadError::Negative { .. }, InvalidKind => ck::PropertyULongReadError::InvalidKind { .. },
    LexicalOverflow => ck::PropertyULongReadError::LexicalOverflow { .. }, Lexical => ck::PropertyULongReadError::Lexical { .. },
}, {
    #[wasm_bindgen(getter)] pub fn value(&self) -> Result<Option<PropertyNumberText>, JsValue> { Ok(match self.inner.as_ref() {
        ck::PropertyULongReadError::Negative { value } => Some(JsValue::from(*value).unchecked_into()),
        ck::PropertyULongReadError::LexicalOverflow { value } | ck::PropertyULongReadError::Lexical { value, .. } => Some(JsValue::from_str(&text(value)?).unchecked_into()),
        ck::PropertyULongReadError::InvalidKind { .. } => None,
    }) }
    #[wasm_bindgen(getter)] pub fn kind(&self) -> Option<PropertyValueKind> { match *self.inner { ck::PropertyULongReadError::InvalidKind { kind } => Some(kind.into()), _ => None } }
}, source = LexicalCause, cause = |value: &ck::PropertyULongReadError| match value {
    ck::PropertyULongReadError::Lexical { source, .. } => Some(UIntLexicalReadError::from_rust(source.clone()).unchecked_into()),
    ck::PropertyULongReadError::Negative { .. } | ck::PropertyULongReadError::InvalidKind { .. } | ck::PropertyULongReadError::LexicalOverflow { .. } => None,
});

error_projection!(MissingPropertyError, ck::MissingPropertyError, {
    #[wasm_bindgen(getter)]
    pub fn key(&self) -> Result<String, JsValue> {
        text(self.inner.key())
    }
});
error_projection!(PropertyStringError, ck::PropertyStringError, {
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> PropertyValueKind {
        self.inner.kind().into()
    }
});
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum PropertyStringErrorKind {
    UnsupportedKind,
}
#[wasm_bindgen]
impl PropertyStringError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> Result<PropertyStringErrorKind, JsValue> {
        match *self.inner {
            ck::PropertyStringError::UnsupportedKind { .. } => {
                Ok(PropertyStringErrorKind::UnsupportedKind)
            }
            _ => Err(
                js_sys::TypeError::new("unrecognized native PropertyStringError variant").into(),
            ),
        }
    }
}
error_projection!(
    RequiredPropertyStringError,
    ck::RequiredPropertyStringError,
    {
        #[wasm_bindgen(getter)]
        pub fn missing(&self) -> Option<MissingPayload> {
            match self.inner.as_ref() {
                ck::RequiredPropertyStringError::Missing(value) => {
                    Some(MissingPropertyError::from_rust(value.clone()).unchecked_into())
                }
                _ => None,
            }
        }
        #[wasm_bindgen(getter)]
        pub fn conversion(&self) -> Option<StringPayload> {
            match self.inner.as_ref() {
                ck::RequiredPropertyStringError::Conversion(value) => {
                    Some(PropertyStringError::from_rust(*value).unchecked_into())
                }
                _ => None,
            }
        }
    }
);
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum RequiredPropertyStringErrorKind {
    Missing,
    Conversion,
}
#[wasm_bindgen]
impl RequiredPropertyStringError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> Result<RequiredPropertyStringErrorKind, JsValue> {
        match *self.inner {
            ck::RequiredPropertyStringError::Missing(_) => {
                Ok(RequiredPropertyStringErrorKind::Missing)
            }
            ck::RequiredPropertyStringError::Conversion(_) => {
                Ok(RequiredPropertyStringErrorKind::Conversion)
            }
            _ => Err(js_sys::TypeError::new(
                "unrecognized native RequiredPropertyStringError variant",
            )
            .into()),
        }
    }
}
