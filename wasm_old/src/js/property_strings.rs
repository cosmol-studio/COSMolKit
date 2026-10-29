//! Exact-type property projections; no value-kind inference at the JS boundary.
use crate::Molecule;
use crate::alignment_values::set;
use crate::host_values::usize_value;
use crate::property_values::kind;
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct PropertyStringError {
    inner: ck::PropertyStringError,
}
#[wasm_bindgen]
impl PropertyStringError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "property".into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        JsValue::NULL
    }
    #[wasm_bindgen(unchecked_return_type = "PropertyValueKind")]
    pub fn kind(&self) -> u32 {
        kind(self.inner.kind())
    }
}
pub(crate) fn string_error(source: &ck::PropertyStringError) -> Result<JsValue, JsValue> {
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("PropertyStringError");
    let e: JsValue = e.into();
    set(&e, "domain", "property".into())?;
    set(&e, "kind", kind(source.kind()).into())?;
    set(&e, "detail", PropertyStringError { inner: *source }.into())?;
    Ok(e)
}

#[wasm_bindgen]
pub struct AtomPropertyReadError {
    inner: ck::AtomPropertyReadError,
}
#[wasm_bindgen]
impl AtomPropertyReadError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "property".into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> String {
        match &self.inner {
            ck::AtomPropertyReadError::Operation(_) => "Operation",
            ck::AtomPropertyReadError::Type(_) => "Type",
            ck::AtomPropertyReadError::Utf8(_) => "Utf8",
        }
        .into()
    }
    #[wasm_bindgen(getter, js_name=validUpTo)]
    pub fn valid_up_to(&self) -> Option<usize> {
        match &self.inner {
            ck::AtomPropertyReadError::Utf8(e) => Some(e.valid_up_to()),
            _ => None,
        }
    }
    #[wasm_bindgen(getter, js_name=errorLen)]
    pub fn error_len(&self) -> Option<usize> {
        match &self.inner {
            ck::AtomPropertyReadError::Utf8(e) => e.error_len(),
            _ => None,
        }
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        match &self.inner {
            ck::AtomPropertyReadError::Operation(e) => crate::alignment_values::operation_error(e),
            ck::AtomPropertyReadError::Type(e) => crate::property_values::property_error(e),
            ck::AtomPropertyReadError::Utf8(e) => Ok(js_sys::Error::new(&e.to_string()).into()),
        }
    }
}
fn read_error(source: ck::AtomPropertyReadError) -> JsValue {
    let result = (|| -> Result<JsValue, JsValue> {
        let e = js_sys::Error::new(&source.to_string());
        e.set_name("AtomPropertyReadError");
        let e: JsValue = e.into();
        let detail = AtomPropertyReadError { inner: source };
        set(&e, "domain", "property".into())?;
        set(&e, "variant", detail.variant().into())?;
        set(&e, "cause", detail.cause()?)?;
        set(&e, "detail", detail.into())?;
        Ok(e)
    })();
    result.unwrap_or_else(|e| e)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=atomPropertyInt, unchecked_return_type="number | null")]
    pub fn atom_property_int(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
    ) -> Result<JsValue, JsValue> {
        self.inner
            .atom_property_int(usize_value(&atom, "atom")?, key)
            .map_err(read_error)
            .map(|v| v.map(|v| JsValue::from(v)).unwrap_or(JsValue::NULL))
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=withAtomPropertyInt)]
    pub fn with_atom_property_int(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<Molecule, JsValue> {
        let value = crate::host_values::i32_value(&value, "value")?;
        self.inner
            .with_atom_property_int(usize_value(&atom, "atom")?, key, value)
            .map(|inner| Self {
                inner: std::sync::Arc::new(inner),
            })
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=setAtomPropertyInt_)]
    pub fn set_atom_property_int_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        let value = crate::host_values::i32_value(&value, "value")?;
        self.inner
            .set_atom_property_int_(usize_value(&atom, "atom")?, key, value)
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=atomPropertyUint, unchecked_return_type="number | null")]
    pub fn atom_property_uint(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
    ) -> Result<JsValue, JsValue> {
        self.inner
            .atom_property_uint(usize_value(&atom, "atom")?, key)
            .map_err(read_error)
            .map(|v| v.map(|v| JsValue::from(v)).unwrap_or(JsValue::NULL))
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=withAtomPropertyUint)]
    pub fn with_atom_property_uint(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<Molecule, JsValue> {
        let value = crate::host_values::u32_value(&value, "value")?;
        self.inner
            .with_atom_property_uint(usize_value(&atom, "atom")?, key, value)
            .map(|inner| Self {
                inner: std::sync::Arc::new(inner),
            })
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=setAtomPropertyUint_)]
    pub fn set_atom_property_uint_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        let value = crate::host_values::u32_value(&value, "value")?;
        self.inner
            .set_atom_property_uint_(usize_value(&atom, "atom")?, key, value)
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=atomPropertyDouble, unchecked_return_type="number | null")]
    pub fn atom_property_double(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
    ) -> Result<JsValue, JsValue> {
        self.inner
            .atom_property_double(usize_value(&atom, "atom")?, key)
            .map_err(read_error)
            .map(|v| v.map(|v| JsValue::from(v)).unwrap_or(JsValue::NULL))
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=withAtomPropertyDouble)]
    pub fn with_atom_property_double(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<Molecule, JsValue> {
        let value = value
            .as_f64()
            .ok_or_else(|| crate::host_values::type_error("value: number"))?;
        self.inner
            .with_atom_property_double(usize_value(&atom, "atom")?, key, value)
            .map(|inner| Self {
                inner: std::sync::Arc::new(inner),
            })
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=setAtomPropertyDouble_)]
    pub fn set_atom_property_double_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        let value = value
            .as_f64()
            .ok_or_else(|| crate::host_values::type_error("value: number"))?;
        self.inner
            .set_atom_property_double_(usize_value(&atom, "atom")?, key, value)
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=atomPropertyBool, unchecked_return_type="boolean | null")]
    pub fn atom_property_bool(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
    ) -> Result<JsValue, JsValue> {
        self.inner
            .atom_property_bool(usize_value(&atom, "atom")?, key)
            .map_err(read_error)
            .map(|v| v.map(|v| JsValue::from(v)).unwrap_or(JsValue::NULL))
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=withAtomPropertyBool)]
    pub fn with_atom_property_bool(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<Molecule, JsValue> {
        let value = crate::host_values::bool_value(&value, "value")?;
        self.inner
            .with_atom_property_bool(usize_value(&atom, "atom")?, key, value)
            .map(|inner| Self {
                inner: std::sync::Arc::new(inner),
            })
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=setAtomPropertyBool_)]
    pub fn set_atom_property_bool_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        let value = crate::host_values::bool_value(&value, "value")?;
        self.inner
            .set_atom_property_bool_(usize_value(&atom, "atom")?, key, value)
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=atomPropertyString, unchecked_return_type="string | null")]
    pub fn atom_property_string(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
    ) -> Result<JsValue, JsValue> {
        self.inner
            .atom_property_string(usize_value(&atom, "atom")?, key)
            .map_err(read_error)
            .map(|v| v.map(|v| JsValue::from(v)).unwrap_or(JsValue::NULL))
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=withAtomPropertyString)]
    pub fn with_atom_property_string(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "string")] value: JsValue,
    ) -> Result<Molecule, JsValue> {
        let value = value
            .as_string()
            .ok_or_else(|| crate::host_values::type_error("value: string"))?;
        self.inner
            .with_atom_property_string(usize_value(&atom, "atom")?, key, &value)
            .map(|inner| Self {
                inner: std::sync::Arc::new(inner),
            })
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=setAtomPropertyString_)]
    pub fn set_atom_property_string_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "string")] value: JsValue,
    ) -> Result<(), JsValue> {
        let value = value
            .as_string()
            .ok_or_else(|| crate::host_values::type_error("value: string"))?;
        self.inner
            .set_atom_property_string_(usize_value(&atom, "atom")?, key, &value)
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=atomPropertyIntVector, unchecked_return_type="number[] | null")]
    pub fn atom_property_int_vector(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
    ) -> Result<JsValue, JsValue> {
        self.inner
            .atom_property_int_vector(usize_value(&atom, "atom")?, key)
            .map_err(read_error)
            .map(|v| {
                v.map(|v| {
                    v.into_iter()
                        .map(JsValue::from)
                        .collect::<js_sys::Array>()
                        .into()
                })
                .unwrap_or(JsValue::NULL)
            })
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=withAtomPropertyIntVector)]
    pub fn with_atom_property_int_vector(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "number[]")] value: JsValue,
    ) -> Result<Molecule, JsValue> {
        let value = crate::host_values::sequence(&value, "value")?
            .iter()
            .map(|v| crate::host_values::i32_value(&v, "value item"))
            .collect::<Result<Vec<_>, _>>()?;
        self.inner
            .with_atom_property_int_vector(usize_value(&atom, "atom")?, key, &value)
            .map(|inner| Self {
                inner: std::sync::Arc::new(inner),
            })
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=setAtomPropertyIntVector_)]
    pub fn set_atom_property_int_vector_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "number[]")] value: JsValue,
    ) -> Result<(), JsValue> {
        let value = crate::host_values::sequence(&value, "value")?
            .iter()
            .map(|v| crate::host_values::i32_value(&v, "value item"))
            .collect::<Result<Vec<_>, _>>()?;
        self.inner
            .set_atom_property_int_vector_(usize_value(&atom, "atom")?, key, &value)
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=atomPropertyStringVector, unchecked_return_type="string[] | null")]
    pub fn atom_property_string_vector(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
    ) -> Result<JsValue, JsValue> {
        self.inner
            .atom_property_string_vector(usize_value(&atom, "atom")?, key)
            .map_err(read_error)
            .map(|v| {
                v.map(|v| {
                    v.into_iter()
                        .map(JsValue::from)
                        .collect::<js_sys::Array>()
                        .into()
                })
                .unwrap_or(JsValue::NULL)
            })
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=withAtomPropertyStringVector)]
    pub fn with_atom_property_string_vector(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "string[]")] value: JsValue,
    ) -> Result<Molecule, JsValue> {
        let value = crate::host_values::sequence(&value, "value")?
            .iter()
            .map(|v| {
                v.as_string()
                    .ok_or_else(|| crate::host_values::type_error("value item: string"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        self.inner
            .with_atom_property_string_vector(usize_value(&atom, "atom")?, key, &value)
            .map(|inner| Self {
                inner: std::sync::Arc::new(inner),
            })
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=setAtomPropertyStringVector_)]
    pub fn set_atom_property_string_vector_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        key: &str,
        #[wasm_bindgen(unchecked_param_type = "string[]")] value: JsValue,
    ) -> Result<(), JsValue> {
        let value = crate::host_values::sequence(&value, "value")?
            .iter()
            .map(|v| {
                v.as_string()
                    .ok_or_else(|| crate::host_values::type_error("value item: string"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        self.inner
            .set_atom_property_string_vector_(usize_value(&atom, "atom")?, key, &value)
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=bondPropertyString,unchecked_return_type="string | null")]
    pub fn bond_property_string(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] id: JsValue,
        key: &str,
    ) -> Result<JsValue, JsValue> {
        self.inner
            .bond_property_string(usize_value(&id, "id")?, key)
            .map_err(|e| string_error(&e).unwrap_or_else(|e| e))
            .and_then(|v| crate::host_values::optional_text(v.as_ref()))
    }
}
