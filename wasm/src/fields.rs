//! Configuration field shapes. Business associations are declared by bindings.
use crate::{input, storage::NativeField};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use tsify::{Ts, Tsify};
use wasm_bindgen::{JsCast, JsValue};

pub(crate) trait Field: Clone + 'static {
    type Input;
    type Output;
    fn input(value: Self::Input) -> Result<Self, JsValue>;
    fn decode(value: &JsValue) -> Result<Self, JsValue>;
    fn get(value: NativeField<Self>) -> Result<Self::Output, JsValue>;
}

/// A binding declares its detached value association beside the class. The
/// borrowed callback ABI reads the existing JS object without consuming it.
macro_rules! value_field {
    ($native:ty, $wrapper:ident, $input:ident, $name:literal) => {
        crate::fields::value_field!($native, $wrapper, $input, $name, $wrapper);
    };
    ($native:ty, $wrapper:ident, $input:ident, $name:literal, $construct:path) => {
        #[wasm_bindgen::prelude::wasm_bindgen]
        extern "C" {
            #[wasm_bindgen(typescript_type = $name)]
            pub type $input;
        }
        const _: () = {
            #[wasm_bindgen::prelude::wasm_bindgen]
            extern "C" {
                #[wasm_bindgen(catch, js_namespace = Reflect, js_name = apply)]
                fn borrow_value(
                    visit: &mut dyn FnMut(&$wrapper),
                    receiver: &wasm_bindgen::JsValue,
                    arguments: &js_sys::Array,
                ) -> Result<wasm_bindgen::JsValue, wasm_bindgen::JsValue>;
            }
            impl $wrapper {
                pub(crate) fn with_ref<R>(
                    raw: &$input,
                    apply: impl FnOnce(&$native) -> Result<R, wasm_bindgen::JsValue>,
                ) -> Result<R, wasm_bindgen::JsValue> {
                    crate::input::with_instance(
                        raw.as_ref(),
                        stringify!($wrapper),
                        |visit, receiver, args| {
                            borrow_value(&mut |value| visit(value), receiver, args)
                        },
                        |value| apply(&value.0),
                    )
                }
                pub(crate) fn with_optional_ref<R>(
                    raw: Option<&$input>,
                    apply: impl FnOnce(Option<&$native>) -> Result<R, wasm_bindgen::JsValue>,
                ) -> Result<R, wasm_bindgen::JsValue> {
                    match raw {
                        Some(raw) => Self::with_ref(raw, |value| apply(Some(value))),
                        None => apply(None),
                    }
                }
            }
            impl crate::fields::Field for $native {
                type Input = $input;
                type Output = $wrapper;
                fn input(value: Self::Input) -> Result<Self, wasm_bindgen::JsValue> {
                    Self::decode(value.as_ref())
                }
                fn decode(raw: &wasm_bindgen::JsValue) -> Result<Self, wasm_bindgen::JsValue> {
                    crate::input::with_instance(
                        raw,
                        stringify!($wrapper),
                        |visit, receiver, args| {
                            borrow_value(&mut |value| visit(value), receiver, args)
                        },
                        |value| Ok(<$native>::clone(&value.0)),
                    )
                }
                fn get(
                    value: crate::storage::NativeField<Self>,
                ) -> Result<Self::Output, wasm_bindgen::JsValue> {
                    value.read(|value| $construct(value.clone()))
                }
            }
        };
    };
}
pub(crate) use value_field;

/// Readonly arrays of native exported values. The element association is
/// declared next to its binding; no business-type dispatch table is needed.
macro_rules! value_array {
    ($native:ty,$array:ident,$typescript:literal) => {
        #[wasm_bindgen::prelude::wasm_bindgen]
        extern "C" {
            #[wasm_bindgen(typescript_type=$typescript)]
            pub type $array;
        }
        impl crate::fields::Field for Vec<$native> {
            type Input = $array;
            type Output = $array;
            fn input(value: Self::Input) -> Result<Self, wasm_bindgen::JsValue> {
                Self::decode(value.as_ref())
            }
            fn decode(value: &wasm_bindgen::JsValue) -> Result<Self, wasm_bindgen::JsValue> {
                if !js_sys::Array::is_array(value) {
                    return Err(js_sys::TypeError::new("expected an array").into());
                }
                use wasm_bindgen::JsCast;
                value
                    .unchecked_ref::<js_sys::Array>()
                    .iter()
                    .map(|v| <$native as crate::fields::Field>::decode(&v))
                    .collect()
            }
            fn get(
                value: crate::storage::NativeField<Self>,
            ) -> Result<Self::Output, wasm_bindgen::JsValue> {
                use wasm_bindgen::JsCast;
                let array = value.read(|values| {
                    let array = js_sys::Array::new();
                    for value in values {
                        let projected = <$native as crate::fields::Field>::get(
                            crate::storage::NativeField::owned(value.clone()),
                        )?;
                        array.push(&wasm_bindgen::JsValue::from(projected));
                    }
                    Ok::<_, wasm_bindgen::JsValue>(array)
                })??;
                Ok(js_sys::Object::freeze(&array).unchecked_into())
            }
        }
    };
}
pub(crate) use value_array;

macro_rules! optional_value_array {
    ($native:ty,$wrapper:ident,$array:ident,$typescript:literal) => {
        #[wasm_bindgen::prelude::wasm_bindgen]
        extern "C" {
            #[wasm_bindgen(typescript_type=$typescript)]
            pub type $array;
        }
        impl $array {
            pub(crate) fn native_values(
                &self,
            ) -> Result<Vec<Option<$native>>, wasm_bindgen::JsValue> {
                use wasm_bindgen::JsCast;
                if !js_sys::Array::is_array(self.as_ref()) {
                    return Err(js_sys::TypeError::new("expected an array").into());
                }
                self.unchecked_ref::<js_sys::Array>()
                    .iter()
                    .map(|value| {
                        if value.is_null() || value.is_undefined() {
                            Ok(None)
                        } else {
                            <$native as crate::fields::Field>::decode(&value).map(Some)
                        }
                    })
                    .collect()
            }
            pub(crate) fn from_native(values: Vec<Option<$native>>) -> Self {
                use wasm_bindgen::JsCast;
                let result = js_sys::Array::new();
                for value in values {
                    result.push(
                        &value
                            .map(|value| wasm_bindgen::JsValue::from($wrapper(value)))
                            .unwrap_or(wasm_bindgen::JsValue::UNDEFINED),
                    );
                }
                js_sys::Object::freeze(&result);
                wasm_bindgen::JsValue::from(result).unchecked_into()
            }
        }
    };
}
pub(crate) use optional_value_array;

macro_rules! scalar {
    ($($ty:ty),+ $(,)?) => {$ (
        impl Field for $ty {
            type Input = Ts<<Self as input::CheckedValue>::Wire>;
            type Output = Self;
            fn input(value: Self::Input) -> Result<Self, JsValue> { input::read(value) }
            fn decode(value: &JsValue) -> Result<Self, JsValue> {
                input::read(Ts::new_unchecked(value.clone()))
            }
            fn get(value: NativeField<Self>) -> Result<Self, JsValue> { value.read(Clone::clone) }
        }
    )+};
}
scalar!(
    bool, u8, u16, u32, i8, i16, i32, usize, f32, f64, String, u64
);

impl<T: Field> Field for Option<T> {
    type Input = Option<T::Input>;
    type Output = Option<T::Output>;
    fn input(value: Self::Input) -> Result<Self, JsValue> {
        value.map(T::input).transpose()
    }
    fn decode(value: &JsValue) -> Result<Self, JsValue> {
        if value.is_null() || value.is_undefined() {
            Ok(None)
        } else {
            T::decode(value).map(Some)
        }
    }
    fn get(value: NativeField<Self>) -> Result<Self::Output, JsValue> {
        if value.read(Option::is_some)? {
            T::get(value.child(Option::as_ref, Option::as_mut)).map(Some)
        } else {
            Ok(None)
        }
    }
}

macro_rules! container {
    ($wire:ident, $native:ty, $ts:literal, $representation:expr) => {
        #[derive(Deserialize, Serialize, Tsify)]
        #[serde(transparent)]
        #[tsify(hashmap_as_object)]
        pub struct $wire(#[tsify(type = $ts)] pub $native);
        impl Field for $native {
            type Input = Ts<$wire>;
            type Output = Ts<$wire>;
            fn input(value: Self::Input) -> Result<Self, JsValue> {
                if !($representation)(&value.js_value()) {
                    return Err(js_sys::TypeError::new(concat!("expected ", $ts)).into());
                }
                Ok(input::decode(value)?.0)
            }
            fn decode(value: &JsValue) -> Result<Self, JsValue> {
                Self::input(Ts::new_unchecked(value.clone()))
            }
            fn get(value: NativeField<Self>) -> Result<Self::Output, JsValue> {
                let result = $wire(value.read(Clone::clone)?)
                    .into_ts()
                    .map_err(wasm_bindgen::JsError::from)?;
                js_sys::Object::freeze(result.js_value().unchecked_ref::<js_sys::Object>());
                Ok(result)
            }
        }
    };
}
container!(
    IndexArray,
    Vec<usize>,
    "ReadonlyArray<number>",
    js_sys::Array::is_array
);
container!(
    BoolArray,
    Vec<bool>,
    "ReadonlyArray<boolean>",
    js_sys::Array::is_array
);
container!(
    U32Array,
    Vec<u32>,
    "ReadonlyArray<number>",
    js_sys::Array::is_array
);
container!(
    I32Array,
    Vec<i32>,
    "ReadonlyArray<number>",
    js_sys::Array::is_array
);
container!(
    FloatArray,
    Vec<f64>,
    "ReadonlyArray<number>",
    js_sys::Array::is_array
);
container!(
    TextArray,
    Vec<String>,
    "ReadonlyArray<string>",
    js_sys::Array::is_array
);
container!(TextMap, BTreeMap<String, String>, "Readonly<Record<string, string>>", input::is_plain_options);

#[derive(Deserialize, Serialize, Tsify)]
#[serde(transparent)]
pub struct IndexPoint2Map(
    #[tsify(type = "ReadonlyMap<number, readonly [number, number]>")]
    #[serde(deserialize_with = "exact_array_map")]
    pub BTreeMap<usize, [f64; 2]>,
);
fn exact_array_map<'de, D, K, T, const N: usize>(
    deserializer: D,
) -> Result<BTreeMap<K, [T; N]>, D::Error>
where
    D: serde::Deserializer<'de>,
    K: Deserialize<'de> + Ord,
    T: Deserialize<'de>,
{
    // serde-wasm-bindgen's tuple reader can ignore trailing elements. Decode
    // the complete sequence, then use Rust's exact fixed-array conversion.
    BTreeMap::<K, Vec<T>>::deserialize(deserializer)?
        .into_iter()
        .map(|(key, value)| input::fixed_array(value).map(|array| (key, array)))
        .collect()
}

#[derive(Deserialize, Serialize, Tsify)]
#[serde(transparent)]
pub struct I32Point3Map(
    #[tsify(type = "ReadonlyMap<number, readonly [number, number, number]>")]
    #[serde(deserialize_with = "exact_array_map")]
    pub BTreeMap<i32, [f64; 3]>,
);
impl Field for BTreeMap<i32, [f64; 3]> {
    type Input = Ts<I32Point3Map>;
    type Output = Ts<I32Point3Map>;
    fn input(value: Self::Input) -> Result<Self, JsValue> {
        if !value.js_value().is_instance_of::<js_sys::Map>() {
            return Err(js_sys::TypeError::new("expected Map").into());
        }
        Ok(input::decode(value)?.0)
    }
    fn decode(value: &JsValue) -> Result<Self, JsValue> {
        Self::input(Ts::new_unchecked(value.clone()))
    }
    fn get(value: NativeField<Self>) -> Result<Self::Output, JsValue> {
        let value = I32Point3Map(value.read(Clone::clone)?)
            .into_ts()
            .map_err(wasm_bindgen::JsError::from)?;
        Ok(Ts::new_unchecked(readonly_array_map(&value.js_value())?))
    }
}

#[derive(Deserialize, Serialize, Tsify)]
#[serde(transparent)]
pub struct PairFloatMap(
    #[tsify(type = "ReadonlyMap<readonly [number, number], number>")]
    #[serde(deserialize_with = "exact_pair_key_map")]
    pub BTreeMap<(u32, u32), f64>,
);
fn exact_pair_key_map<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<BTreeMap<(u32, u32), f64>, D::Error> {
    BTreeMap::<Vec<u32>, f64>::deserialize(deserializer)?
        .into_iter()
        .map(|(key, value)| {
            let [a, b] = input::fixed_array::<u32, 2, D::Error>(key)?;
            Ok(((a, b), value))
        })
        .collect()
}
impl Field for BTreeMap<(u32, u32), f64> {
    type Input = Ts<PairFloatMap>;
    type Output = Ts<PairFloatMap>;
    fn input(value: Self::Input) -> Result<Self, JsValue> {
        if !value.js_value().is_instance_of::<js_sys::Map>() {
            return Err(js_sys::TypeError::new("expected Map").into());
        }
        Ok(input::decode(value)?.0)
    }
    fn decode(value: &JsValue) -> Result<Self, JsValue> {
        Self::input(Ts::new_unchecked(value.clone()))
    }
    fn get(value: NativeField<Self>) -> Result<Self::Output, JsValue> {
        let map = js_sys::Map::new();
        value.read(|values| {
            for (&(a, b), &value) in values {
                let key = js_sys::Array::of2(&a.into(), &b.into());
                js_sys::Object::freeze(&key);
                map.set(&key, &value.into());
            }
        })?;
        Ok(Ts::new_unchecked(readonly_map(&map)?))
    }
}
impl Field for BTreeMap<usize, [f64; 2]> {
    type Input = Ts<IndexPoint2Map>;
    type Output = Ts<IndexPoint2Map>;
    fn input(value: Self::Input) -> Result<Self, JsValue> {
        if !value.js_value().is_instance_of::<js_sys::Map>() {
            return Err(js_sys::TypeError::new("expected Map").into());
        }
        Ok(input::decode(value)?.0)
    }
    fn decode(value: &JsValue) -> Result<Self, JsValue> {
        Self::input(Ts::new_unchecked(value.clone()))
    }
    fn get(value: NativeField<Self>) -> Result<Self::Output, JsValue> {
        let result = IndexPoint2Map(value.read(Clone::clone)?)
            .into_ts()
            .map_err(wasm_bindgen::JsError::from)?;
        let raw = result.js_value();
        Ok(Ts::new_unchecked(readonly_array_map(&raw)?))
    }
}

/// Readonly maps whose values are numeric arrays share the same ownership
/// boundary, independent of the business field using that shape.
fn readonly_array_map(raw: &JsValue) -> Result<JsValue, JsValue> {
    let map = raw.unchecked_ref::<js_sys::Map>();
    for entry in
        js_sys::try_iter(map)?.ok_or_else(|| js_sys::TypeError::new("Map iterator is absent"))?
    {
        let entry = js_sys::Array::from(&entry?);
        js_sys::Object::freeze(entry.get(1).unchecked_ref::<js_sys::Object>());
    }
    readonly_map(map)
}

#[derive(Deserialize, Serialize, Tsify)]
#[serde(transparent)]
pub struct I32ListMap(
    #[tsify(type = "ReadonlyMap<number, ReadonlyArray<number>>")] pub BTreeMap<i32, Vec<i32>>,
);
impl Field for BTreeMap<i32, Vec<i32>> {
    type Input = Ts<I32ListMap>;
    type Output = Ts<I32ListMap>;
    fn input(value: Self::Input) -> Result<Self, JsValue> {
        if !value.js_value().is_instance_of::<js_sys::Map>() {
            return Err(js_sys::TypeError::new("expected Map").into());
        }
        Ok(input::decode(value)?.0)
    }
    fn decode(value: &JsValue) -> Result<Self, JsValue> {
        Self::input(Ts::new_unchecked(value.clone()))
    }
    fn get(value: NativeField<Self>) -> Result<Self::Output, JsValue> {
        let result = I32ListMap(value.read(Clone::clone)?)
            .into_ts()
            .map_err(wasm_bindgen::JsError::from)?;
        Ok(Ts::new_unchecked(readonly_array_map(&result.js_value())?))
    }
}

#[wasm_bindgen::prelude::wasm_bindgen(module = "/src/readonly_map.js")]
extern "C" {
    #[wasm_bindgen(catch,js_name=readonlyMap)]
    pub(crate) fn readonly_map(map: &js_sys::Map) -> Result<JsValue, JsValue>;
}
