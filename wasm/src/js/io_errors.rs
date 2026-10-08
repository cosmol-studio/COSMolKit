//! Canonical SDF and molecular IO error projection, including structured sources.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;

pub(crate) fn filesystem_error(source: &std::io::Error) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("IoError");
    let error: JsValue = error.into();
    set(&error, "domain", "io".into())?;
    set(
        &error,
        "kind",
        JsValue::from_str(&format!("{:?}", source.kind())),
    )?;
    set(
        &error,
        "errno",
        source
            .raw_os_error()
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64)),
    )?;
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    Ok(error)
}
#[wasm_bindgen]
pub struct SdfError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: js_sys::Object,
}
#[wasm_bindgen]
impl SdfError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "io".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter,js_name=expected,unchecked_return_type="string | null")]
    pub fn field_0(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"expected".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=actual,unchecked_return_type="string | null")]
    pub fn field_1(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"actual".into()).unwrap_or(JsValue::NULL)
    }
}
pub(crate) fn sdf_error(source: &ck::SdfError) -> Result<JsValue, JsValue> {
    use ck::SdfError as E;
    let fields = js_sys::Object::new();
    set(&fields, "expected", JsValue::NULL)?;
    set(&fields, "actual", JsValue::NULL)?;
    let kind = match source {
        E::Read(_) => "Read",
        E::Post(_) => "Post",
        E::Construction(_) => "Construction",
        E::QueryGraph(_) => "QueryGraph",
        E::QueryRecord => "QueryRecord",
        E::WrongGraphKind { expected, actual } => {
            set(&fields, "expected", (*expected).into())?;
            set(&fields, "actual", (*actual).into())?;
            "WrongGraphKind"
        }
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("SdfError");
    let e: JsValue = error.into();
    set(&e, "domain", "io".into())?;
    set(&e, "kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    for key in js_sys::Object::keys(&fields).iter() {
        let value = js_sys::Reflect::get(&fields, &key)?;
        if !value.is_null() {
            set(&e, &key.as_string().unwrap(), value)?;
        }
    }
    set(
        &e,
        "detail",
        SdfError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
pub struct MolecularIoError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: js_sys::Object,
}
#[wasm_bindgen]
impl MolecularIoError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "io".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter,js_name=filename,unchecked_return_type="string | null")]
    pub fn field_0(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"filename".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=errno,unchecked_return_type="number | null")]
    pub fn field_1(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"errno".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=ioKind,unchecked_return_type="string | null")]
    pub fn field_2(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"ioKind".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=format,unchecked_return_type="string | null")]
    pub fn field_3(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"format".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=parameterName,unchecked_return_type="string | null")]
    pub fn field_4(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"parameterName".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=parameterDetail,unchecked_return_type="string | null")]
    pub fn field_5(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"parameterDetail".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=stage,unchecked_return_type="string | null")]
    pub fn field_6(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"stage".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=coordinateKind,unchecked_return_type="string | null")]
    pub fn field_7(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"coordinateKind".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=twoD,unchecked_return_type="number | null")]
    pub fn field_8(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"twoD".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=threeD,unchecked_return_type="number | null")]
    pub fn field_9(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"threeD".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=dimension,unchecked_return_type="string | null")]
    pub fn field_10(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"dimension".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=coordinateId,unchecked_return_type="number | null")]
    pub fn field_11(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"coordinateId".into()).unwrap_or(JsValue::NULL)
    }
}
pub(crate) fn io_error(source: &ck::MolecularIoError) -> Result<JsValue, JsValue> {
    use ck::MolecularIoError as E;
    if let E::Sdf(e) = source {
        return sdf_error(e);
    }
    let fields = js_sys::Object::new();
    set(&fields, "filename", JsValue::NULL)?;
    set(&fields, "errno", JsValue::NULL)?;
    set(&fields, "ioKind", JsValue::NULL)?;
    set(&fields, "format", JsValue::NULL)?;
    set(&fields, "parameterName", JsValue::NULL)?;
    set(&fields, "parameterDetail", JsValue::NULL)?;
    set(&fields, "stage", JsValue::NULL)?;
    set(&fields, "coordinateKind", JsValue::NULL)?;
    set(&fields, "twoD", JsValue::NULL)?;
    set(&fields, "threeD", JsValue::NULL)?;
    set(&fields, "dimension", JsValue::NULL)?;
    set(&fields, "coordinateId", JsValue::NULL)?;
    let kind = match source {
        E::OutputUtf8 { .. } => "OutputUtf8",
        E::Io { path, source } => {
            set(&fields, "filename", path.to_string_lossy().as_ref().into())?;
            set(
                &fields,
                "errno",
                source.raw_os_error().map_or(JsValue::NULL, JsValue::from),
            )?;
            set(&fields, "ioKind", format!("{:?}", source.kind()).into())?;
            "Io"
        }
        E::NoRecord { format } => {
            set(&fields, "format", (*format).into())?;
            "NoRecord"
        }
        E::Parameter { name, detail } => {
            set(&fields, "parameterName", (*name).into())?;
            set(&fields, "parameterDetail", (*detail).into())?;
            "Parameter"
        }
        E::Mol2Post(e) => {
            set(&fields, "stage", e.stage.into())?;
            "Mol2Post"
        }
        E::MolWrite(e) => {
            match e {
                ck::MolWriteError::AmbiguousCoordinates { two_d, three_d } => {
                    set(&fields, "coordinateKind", "AmbiguousCoordinates".into())?;
                    set(&fields, "twoD", (*two_d as u32).into())?;
                    set(&fields, "threeD", (*three_d as u32).into())?;
                }
                ck::MolWriteError::MissingCoordinate { dimension, id } => {
                    set(&fields, "coordinateKind", "MissingCoordinate".into())?;
                    set(
                        &fields,
                        "dimension",
                        match dimension {
                            ck::CoordinateDimension::TwoD => "2d",
                            ck::CoordinateDimension::ThreeD => "3d",
                        }
                        .into(),
                    )?;
                    set(&fields, "coordinateId", (*id as u32).into())?;
                }
                ck::MolWriteError::CoordinateDimensionMismatch => set(
                    &fields,
                    "coordinateKind",
                    "CoordinateDimensionMismatch".into(),
                )?,
                _ => {}
            };
            "MolWrite"
        }
        E::XyzRead(_) => "XyzRead",
        E::XyzWrite(_) => "XyzWrite",
        E::Mol2Read(_) => "Mol2Read",
        E::Construction(_) => "Construction",
        E::Sdf(_) => unreachable!("Sdf handled above"),
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("MolecularIoError");
    let e: JsValue = error.into();
    set(&e, "domain", "io".into())?;
    set(&e, "kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    for key in js_sys::Object::keys(&fields).iter() {
        let value = js_sys::Reflect::get(&fields, &key)?;
        if !value.is_null() {
            set(&e, &key.as_string().unwrap(), value)?;
        }
    }
    set(
        &e,
        "detail",
        MolecularIoError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(e)
}
