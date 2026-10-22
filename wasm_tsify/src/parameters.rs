//! Native scalar configuration projection. Unsupported shapes fail compilation.
//! No business-type map or duplicated native field/default declarations.
#[macro_export]
macro_rules! project_parameters {
    (@resume {$name:ident} @native $native:path { $($field:ident: $ty:ty),* } methods {}) => {
        #[wasm_bindgen::prelude::wasm_bindgen]
        pub struct $name {
            inner: std::rc::Rc<std::cell::RefCell<$native>>,
        }

        #[wasm_bindgen::prelude::wasm_bindgen]
        impl $name {
            #[wasm_bindgen(constructor)]
            pub fn new() -> Self {
                Self { inner: std::rc::Rc::new(std::cell::RefCell::new(<$native>::default())) }
            }
            pub fn copy(&self) -> Self {
                Self { inner: std::rc::Rc::new(std::cell::RefCell::new(self.inner.borrow().clone())) }
            }
        }
        paste::paste! {
            #[wasm_bindgen::prelude::wasm_bindgen]
            impl $name {
                $(
                    #[wasm_bindgen(getter)]
                    pub fn $field(&self) -> $ty { self.inner.borrow().$field.clone() }
                    #[wasm_bindgen(setter)]
                    pub fn [<set_ $field>](&self, value: $ty) { self.inner.borrow_mut().$field = value; }
                )*
            }
        }
    };
}
