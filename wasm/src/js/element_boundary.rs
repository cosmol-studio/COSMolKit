//! Exact JavaScript transport for the canonical element projection.
//!
//! Alef includes this module through `custom_rust_modules`. In particular,
//! wasm-bindgen's native `Option<Class>` ABI returns undefined, whereas the
//! registry's language contract requires null. Keep this conversion here.

use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct Element {
    inner: cosmolkit_wasm::Element,
}

#[wasm_bindgen]
impl Element {
    #[wasm_bindgen(js_name = fromAtomicNumber, unchecked_return_type = "Element | null")]
    pub fn from_atomic_number(atomic_number: f64) -> Result<JsValue, JsValue> {
        if !atomic_number.is_finite()
            || atomic_number.fract() != 0.0
            || !(0.0..=f64::from(u8::MAX)).contains(&atomic_number)
        {
            return Err(
                js_sys::RangeError::new("atomicNumber must be an integer in 0..=255").into(),
            );
        }
        Ok(
            cosmolkit_wasm::Element::from_atomic_number(atomic_number as u8)
                .map_or(JsValue::NULL, |inner| Element { inner }.into()),
        )
    }

    #[wasm_bindgen(js_name = fromSymbol, unchecked_return_type = "Element | null")]
    pub fn from_symbol(symbol: &str) -> JsValue {
        cosmolkit_wasm::Element::from_symbol(symbol)
            .map_or(JsValue::NULL, |inner| Element { inner }.into())
    }

    #[wasm_bindgen(js_name = atomicNumber)]
    pub fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }

    pub fn symbol(&self) -> String {
        self.inner.symbol()
    }
}

#[wasm_bindgen]
pub struct ElementInfo {
    inner: cosmolkit_wasm::ElementInfo,
}

#[wasm_bindgen]
impl ElementInfo {
    pub fn element(&self) -> Element {
        Element {
            inner: self.inner.element(),
        }
    }

    pub fn symbol(&self) -> String {
        self.inner.symbol()
    }

    #[wasm_bindgen(js_name = atomicNumber)]
    pub fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }

    pub fn period(&self) -> u8 {
        self.inner.period()
    }

    #[wasm_bindgen(js_name = outerElectrons)]
    pub fn outer_electrons(&self) -> i32 {
        self.inner.outer_electrons()
    }

    pub fn valences(&self) -> Vec<i32> {
        self.inner.valences()
    }

    pub fn rb0(&self) -> f64 {
        self.inner.rb0()
    }

    #[wasm_bindgen(js_name = atomicWeight)]
    pub fn atomic_weight(&self) -> f64 {
        self.inner.atomic_weight()
    }
}

#[wasm_bindgen(js_name = elementInfo)]
pub fn element_info(element: &Element) -> ElementInfo {
    ElementInfo {
        inner: cosmolkit_wasm::element_info(&element.inner),
    }
}
