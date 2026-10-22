//! Binding protocol probe; chemistry remains in the public Rust facade.
use cosmolkit as ck;
use serde::{Deserialize, Serialize};
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;

mod parameters;

// The field list and field types come from the native declaration, not a JS list.
ck::__ck_configuration_RingSearchParams!(crate::project_parameters, ck::RingSearchParams, {}, {
    RingSearchParams
});

#[wasm_bindgen]
pub struct Molecule(ck::Molecule);

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(constructor)]
    pub fn new(smiles: &str) -> Result<Molecule, JsError> {
        ck::Molecule::from_smiles(smiles)
            .map(Self)
            .map_err(|e| JsError::new(&e.to_string()))
    }

    pub fn copy(&self) -> Molecule {
        Self(self.0.copy())
    }

    pub fn num_atoms(&self) -> usize {
        self.0.num_atoms()
    }

    pub fn num_rings(&self, params: &RingSearchParams) -> Result<usize, JsError> {
        self.0
            .symmetrized_sssr_with_params(&params.inner.borrow())
            .map(|r| r.num_rings())
            .map_err(|e| JsError::new(&e.to_string()))
    }

    /// Protocol probe: deliver an independently owned COW snapshot, not a borrow.
    pub fn inspect(&self, callback: &InspectCallback) -> Result<bool, JsValue> {
        let value = callback.call1(&JsValue::UNDEFINED, &JsValue::from(self.copy()))?;
        value
            .as_bool()
            .ok_or_else(|| js_sys::TypeError::new("callback must return boolean").into())
    }
}

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(extends = js_sys::Function, typescript_type = "(molecule: Molecule) => boolean")]
    pub type InspectCallback;

    #[wasm_bindgen(typescript_type = "unknown")]
    pub type ExternalValue;
}

// These are protocol samples, not new chemistry API or duplicated native Params.
#[derive(Serialize, Deserialize, Tsify)]
#[serde(rename_all = "snake_case")]
pub enum Mode {
    Preserve,
    Replace,
}

#[derive(Serialize, Deserialize, Tsify)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Outcome {
    Completed { count: u32 },
    Failed { message: String },
}

#[derive(Serialize, Deserialize, Tsify)]
#[serde(deny_unknown_fields)]
pub struct ArrayInput {
    pub point: [f64; 3],
    pub indices: Vec<u32>,
    pub mode: Mode,
}

#[wasm_bindgen]
pub fn roundtrip(input: Ts<ArrayInput>) -> Result<Ts<ArrayInput>, JsError> {
    Ok(input.to_rust()?.into_ts()?)
}

#[wasm_bindgen]
pub fn outcome(input: Ts<Outcome>) -> Result<Ts<Outcome>, JsError> {
    Ok(input.to_rust()?.into_ts()?)
}

/// Structured error transport sample; domain-specific mapping belongs to its binding.
#[derive(Serialize, Deserialize, Tsify)]
#[serde(tag = "kind")]
pub enum Failure {
    Syntax { offset: u32, message: String },
    Limit { limit: u32 },
}

#[wasm_bindgen]
pub fn raise(details: Ts<Failure>, cause: &ExternalValue) -> Result<(), JsValue> {
    // Validate before attaching typed data; do not pretend an unchecked JS value is typed.
    let details = details.to_rust().map_err(JsError::from)?;
    let error = js_sys::Error::new("binding protocol failure");
    js_sys::Reflect::set(
        &error,
        &"details".into(),
        &details.into_ts().map_err(JsError::from)?.js_value(),
    )?;
    js_sys::Reflect::set(&error, &"cause".into(), cause.as_ref())?;
    Err(error.into())
}

#[wasm_bindgen]
pub fn failure_details(error: &js_sys::Error) -> Result<Ts<Failure>, JsError> {
    let value = js_sys::Reflect::get(error, &"details".into())
        .map_err(|_| JsError::new("error details are not accessible"))?;
    Ok(Ts::<Failure>::new_unchecked(value).to_rust()?.into_ts()?)
}

/// Owned typed-array output; mutation in JS cannot modify native storage.
#[wasm_bindgen]
pub fn bytes(input: &[u8]) -> Vec<u8> {
    input.to_vec()
}

/// Explicit 64-bit integer boundary: wasm-bindgen generates bigint, not number.
#[wasm_bindgen]
pub fn integer(input: u64) -> u64 {
    input
}

#[derive(Serialize, Deserialize, Tsify, Clone, Default)]
#[serde(transparent)]
pub struct Indices {
    #[tsify(type = "ReadonlyArray<number>")]
    values: Vec<u32>,
}

/// Container protocol: immutable detached reads, explicit whole-field writes.
#[wasm_bindgen]
pub struct Selection {
    indices: Indices,
}

#[wasm_bindgen]
impl Selection {
    #[wasm_bindgen(constructor)]
    pub fn new(indices: Ts<Indices>) -> Result<Self, JsError> {
        Ok(Self {
            indices: indices.to_rust()?,
        })
    }
    #[wasm_bindgen(getter)]
    pub fn indices(&self) -> Result<Ts<Indices>, JsError> {
        let result = self.indices.into_ts()?;
        js_sys::Object::freeze(&result.js_value().unchecked_into::<js_sys::Object>());
        Ok(result)
    }
    #[wasm_bindgen(setter)]
    pub fn set_indices(&mut self, indices: Ts<Indices>) -> Result<(), JsError> {
        self.indices = indices.to_rust()?;
        Ok(())
    }
}

/// A nested configuration protocol sample using the same native field wrapper.
#[wasm_bindgen]
pub struct NestedParams {
    child: RingSearchParams,
}

#[wasm_bindgen]
impl NestedParams {
    #[wasm_bindgen(constructor)]
    pub fn new(child: &RingSearchParams) -> Self {
        Self {
            child: child.copy(),
        }
    }

    #[wasm_bindgen(getter)]
    pub fn child(&self) -> RingSearchParams {
        RingSearchParams {
            inner: self.child.inner.clone(),
        }
    }

    #[wasm_bindgen(setter)]
    pub fn set_child(&mut self, child: &RingSearchParams) {
        // Preserve existing live child handles while assigning an independent value.
        let value = *child.inner.borrow();
        *self.child.inner.borrow_mut() = value;
    }
}
