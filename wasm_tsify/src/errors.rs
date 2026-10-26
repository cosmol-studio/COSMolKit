//! Error transport: retain the typed Rust value and its JavaScript class.
use wasm_bindgen::{JsCast, JsValue, prelude::*};

/// Make a generated domain-error class an Error subclass without replacing its
/// methods, losing its Rust payload, or editing generated JavaScript/declarations.
pub(crate) fn exception(value: JsValue) -> JsValue {
    let object = value.unchecked_ref::<js_sys::Object>();
    let prototype = js_sys::Object::get_prototype_of(object);
    // Name and message are generated native getters on the domain class.
    // This Error is used only to obtain Error.prototype.
    let native_error = js_sys::Error::new("");
    let error_prototype = js_sys::Object::get_prototype_of(&native_error);
    // Built-in Object operations can fail if application code freezes a class
    // prototype. Preserve that actual JavaScript exception rather than panic.
    match js_sys::Object::try_set_prototype_of(&prototype, &error_prototype) {
        Ok(_) => value,
        Err(error) => error,
    }
}

macro_rules! error_projection {
    ($name:ident, $native:ty, {$($accessor:item)*} $(, source = $source:ty, cause = $cause:expr)?) => {
        paste::paste! {
            #[wasm_bindgen]
            pub struct $name {
                inner: crate::storage::ReadonlyValue<$native>,
                $(source: Option<$source>,)?
            }
            #[wasm_bindgen]
            pub struct [<$name Details>] {
                inner: crate::storage::ReadonlyValue<$native>
            }
            #[wasm_bindgen]
            impl [<$name Details>] {
                $($accessor)*
            }

            #[wasm_bindgen]
            impl $name {
                #[wasm_bindgen(getter)]
                pub fn name(&self) -> String {
                    stringify!($name).into()
                }
                #[wasm_bindgen(getter)]
                pub fn message(&self) -> String {
                    self.inner.to_string()
                }
                #[wasm_bindgen(js_name = toString)]
                pub fn display(&self) -> String {
                    self.inner.to_string()
                }
                $(
                    #[wasm_bindgen(getter)]
                    pub fn cause(&self) -> Option<$source> {
                        self.source.clone()
                    }
                )?
                #[wasm_bindgen(getter)]
                pub fn details(&self) -> [<$name Details>] {
                    [<$name Details>] {
                        inner: self.inner.clone()
                    }
                }
            }

            impl $name {
                /// Transparent native error variants delegate to the child's
                /// actual source, rather than inventing an extra cause node.
                #[allow(dead_code, unused_variables)]
                pub(crate) fn native_cause(
                    inner: &crate::storage::ReadonlyValue<$native>
                ) -> Option<wasm_bindgen::JsValue> {
                    let source: Option<wasm_bindgen::JsValue> = None;
                    $(
                        let source = {
                            let typed: Option<$source> = ($cause)(inner);
                            typed.map(wasm_bindgen::JsValue::from)
                        };
                    )?
                    source
                }
                #[allow(dead_code)]
                pub(crate) fn with_native<R>(
                    raw: &JsValue,
                    apply: impl FnOnce(&$native) -> R
                ) -> Result<R, JsValue> {
                    #[wasm_bindgen]
                    extern "C" {
                        #[wasm_bindgen(catch, js_namespace = Reflect, js_name = apply)]
                        fn visit_error(
                            visit: &mut dyn FnMut(&$name),
                            receiver: &JsValue,
                            arguments: &js_sys::Array
                        ) -> Result<JsValue, JsValue>;
                    }
                    crate::input::with_instance(
                        raw,
                        stringify!($name),
                        |visit, receiver, args| visit_error(&mut |value| visit(value), receiver, args),
                        |value| {
                            Ok(apply(value.inner.as_ref()))
                        }
                    )
                }
                pub(crate) fn from_rust(inner: $native) -> JsValue {
                    Self::from_owner(crate::storage::ReadonlyValue::owned(inner))
                }
                pub(crate) fn from_owner(inner: crate::storage::ReadonlyValue<$native>) -> JsValue {
                    $(
                        let source: Option<$source> = ($cause)(&inner);
                    )?
                    crate::errors::exception(
                        Self {
                            inner,
                            $(source: source.map(|value: $source| value),)?
                        }.into()
                    )
                }
            }
        }
    };
}

pub(crate) use error_projection;

macro_rules! enum_error {
    (
        $name:ident, $kind:ident, $native:path,
        {$($(#[$attribute:meta])* $variant:ident => $pattern:pat),+ $(,)?},
        {$($method:item)*}
        $(, source = $source:ty, cause = $cause:expr)?
    ) => {
        #[wasm_bindgen]
        #[derive(Clone, Copy)]
        pub enum $kind {
            // Keep diagnostic codes stable across feature selections. Only
            // the native match arms depend on whether a source variant exists.
            $($variant),+
        }
        crate::errors::error_projection!(
            $name, $native, { $($method)* }
            $(, source = $source, cause = $cause)?
        );
        #[wasm_bindgen]
        impl $name {
            #[wasm_bindgen(getter)]
            pub fn variant(&self) -> $kind {
                match self.inner.as_ref() {
                    $($(#[$attribute])* $pattern => $kind::$variant),+
                }
            }
        }
    };
}
pub(crate) use enum_error;

error_projection!(ElementParseError, cosmolkit::ElementParseError, {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "element".into()
    }
    #[wasm_bindgen(getter)]
    pub fn input(&self) -> String {
        self.inner.input().to_owned()
    }
});
