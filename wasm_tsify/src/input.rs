//! Checked value inputs. Serde owns numeric/type validation, not this binding.
//! Only primitive associations belong here; business types declare their own.
use serde::{Deserialize, de::DeserializeOwned};
use tsify::{Ts, Tsify};
use wasm_bindgen::{JsCast, JsValue};

pub(crate) fn fixed_array<T, const N: usize, E: serde::de::Error>(
    value: Vec<T>,
) -> Result<[T; N], E> {
    let length = value.len();
    value
        .try_into()
        .map_err(|_: Vec<T>| E::invalid_length(length, &"an exact-size tuple"))
}

/// The WASM Serde tuple decoder can ignore trailing array elements. Read the
/// whole sequence before Rust's checked fixed-array conversion.
pub(crate) fn exact_array<'de, D, T, const N: usize>(deserializer: D) -> Result<[T; N], D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    fixed_array(Vec::<T>::deserialize(deserializer)?)
}

pub(crate) fn decode<T: Tsify + DeserializeOwned>(value: Ts<T>) -> Result<T, JsValue>
where
    T::JsType: Clone,
{
    value
        .to_rust()
        .map_err(|error| js_sys::TypeError::new(&error.to_string()).into())
}

pub trait CheckedValue: Sized {
    type Wire: Tsify + DeserializeOwned;
    fn unwrap(value: Self::Wire) -> Self;
    fn check_representation(_value: &JsValue) -> Result<(), JsValue> {
        Ok(())
    }
}

/// The borrowed ABI validates the exported class before entering Rust. Only
/// that extraction phase is an input TypeError; operation/callback exceptions
/// after entry retain their original object and classification.
fn instance_error(error: JsValue, entered: bool, expected: &str) -> JsValue {
    if entered {
        error
    } else {
        js_sys::TypeError::new(&format!("expected {expected}")).into()
    }
}

/// Borrow a concrete exported class through its generated wasm-bindgen bridge.
/// The bridge owns class/pointer validation; the caller owns native read/write
/// policy. The borrowed wrapper cannot escape through the independent result R.
pub(crate) fn with_instance<W, R>(
    raw: &JsValue,
    expected: &str,
    bridge: impl FnOnce(&mut dyn FnMut(&W), &JsValue, &js_sys::Array) -> Result<JsValue, JsValue>,
    apply: impl FnOnce(&W) -> Result<R, JsValue>,
) -> Result<R, JsValue> {
    let mut apply = Some(apply);
    let mut result = None;
    let mut entered = false;
    bridge(
        &mut |value| {
            entered = true;
            result = Some(match apply.take() {
                Some(apply) => apply(value),
                None => {
                    Err(js_sys::TypeError::new("native visitor was invoked more than once").into())
                }
            });
        },
        &JsValue::UNDEFINED,
        &js_sys::Array::of1(raw),
    )
    .map_err(|error| instance_error(error, entered, expected))?;
    result.ok_or_else(|| JsValue::from(js_sys::TypeError::new("native visitor was not invoked")))?
}

pub(crate) fn is_plain_options(value: &JsValue) -> bool {
    if !value.is_object() || value.is_null() || js_sys::Array::is_array(value) {
        return false;
    }
    let object = value.unchecked_ref::<js_sys::Object>();
    let prototype = js_sys::Object::get_prototype_of(object);
    prototype.is_null() || prototype == js_sys::Object::get_prototype_of(&js_sys::Object::new())
}

/// serde-wasm-bindgen reads declared struct fields directly, so unknown keys
/// must be rejected before decoding. Names are generated from native metadata.
pub(crate) fn check_options(raw: &JsValue, fields: &[&str]) -> Result<(), JsValue> {
    if !is_plain_options(raw) {
        return Err(js_sys::TypeError::new("expected an options object").into());
    }
    for key in js_sys::Reflect::own_keys(raw)?.iter() {
        if !key
            .as_string()
            .is_some_and(|key| fields.contains(&key.as_str()))
        {
            return Err(js_sys::TypeError::new("unknown configuration field").into());
        }
    }
    Ok(())
}

pub(crate) fn read<T: CheckedValue>(value: Ts<T::Wire>) -> Result<T, JsValue>
where
    <T::Wire as Tsify>::JsType: Clone,
{
    T::check_representation(&value.js_value())?;
    Ok(T::unwrap(decode(value)?))
}

macro_rules! scalar {
    ($name:ident, $native:ty) => {
        #[derive(Deserialize, Tsify)]
        #[serde(transparent)]
        pub struct $name($native);

        impl CheckedValue for $native {
            type Wire = $name;
            fn unwrap(value: $name) -> Self {
                value.0
            }
        }
    };
}

scalar!(U8Input, u8);
scalar!(U16Input, u16);
scalar!(U32Input, u32);
scalar!(I8Input, i8);
scalar!(I16Input, i16);
scalar!(I32Input, i32);
scalar!(F32Input, f32);
scalar!(F64Input, f64);
scalar!(IndexInput, usize);
scalar!(BoolInput, bool);
scalar!(TextInput, String);
scalar!(CharInput, char);

#[derive(Deserialize, Tsify)]
#[serde(transparent)]
pub struct U64Input(#[tsify(type = "bigint")] u64);

#[derive(Deserialize, Tsify)]
#[serde(transparent)]
pub struct I64Input(#[tsify(type = "bigint")] i64);

impl CheckedValue for i64 {
    type Wire = I64Input;
    fn unwrap(value: I64Input) -> Self {
        value.0
    }
    fn check_representation(value: &JsValue) -> Result<(), JsValue> {
        <u64 as CheckedValue>::check_representation(value)
    }
}

impl CheckedValue for u64 {
    type Wire = U64Input;
    fn unwrap(value: U64Input) -> Self {
        value.0
    }
    fn check_representation(value: &JsValue) -> Result<(), JsValue> {
        // Serde also accepts safe Numbers; this API explicitly declares bigint.
        if value.is_bigint() {
            Ok(())
        } else {
            Err(js_sys::TypeError::new("expected bigint").into())
        }
    }
}

#[derive(Deserialize, Tsify)]
#[serde(transparent)]
pub struct BytesInput(#[tsify(type = "Uint8Array")] Vec<u8>);

impl CheckedValue for Vec<u8> {
    type Wire = BytesInput;
    fn unwrap(value: BytesInput) -> Self {
        value.0
    }
    fn check_representation(value: &JsValue) -> Result<(), JsValue> {
        if value.is_instance_of::<js_sys::Uint8Array>() {
            Ok(())
        } else {
            Err(js_sys::TypeError::new("expected Uint8Array").into())
        }
    }
}
