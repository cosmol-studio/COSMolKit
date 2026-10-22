use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;

macro_rules! flags {
    ($name:ident, $visit:ident, $js_visit:ident, $script:literal, {$($flag:ident),+ $(,)?}) => {
        #[wasm_bindgen]
        pub struct $name { pub(crate) inner: ck::$name }
        #[wasm_bindgen]
        impl $name {
            $(#[wasm_bindgen(getter)] #[allow(non_snake_case)] pub fn $flag() -> Self { Self { inner: ck::$name::$flag } })+
            #[wasm_bindgen(js_name=fromBits)]
            pub fn from_bits(#[wasm_bindgen(unchecked_param_type="number")] bits: JsValue) -> Result<Self, JsValue> {
                Ok(Self { inner: ck::$name::from_bits(crate::host_values::u32_value(&bits, "PDB flags")?) })
            }
            #[wasm_bindgen(getter)]
            pub fn bits(&self) -> u32 { self.inner.bits() }
            pub fn contains(&self, flags: &Self) -> bool { self.inner.contains(flags.inner) }
            #[wasm_bindgen(getter,js_name=isEmpty)]
            pub fn is_empty(&self) -> bool { self.inner.is_empty() }
            pub fn or(&self, flags: &Self) -> Self { Self { inner: self.inner | flags.inner } }
            pub fn and(&self, flags: &Self) -> Self { Self { inner: self.inner & flags.inner } }
        }
        #[wasm_bindgen(inline_js=$script)]
        extern "C" {
            #[wasm_bindgen(catch,js_name=$js_visit)]
            fn $visit(value: &JsValue, visitor: &mut dyn FnMut(&$name)) -> Result<(), JsValue>;
        }
        impl $name {
            pub(crate) fn extract(value: &JsValue) -> Result<ck::$name, JsValue> {
                let mut flags = None;
                $visit(value, &mut |value: &$name| flags = Some(value.inner))
                    .map_err(|_| crate::host_values::type_error(concat!("expected ", stringify!($name))))?;
                flags.ok_or_else(|| crate::host_values::type_error(concat!("expected ", stringify!($name))))
            }
        }
    };
}
flags!(PdbReadFlags, visit_read_flags, visitPdbReadFlags,
    "export function visitPdbReadFlags(v,f){f(v);}",
    {NONE, KEEP_ALTERNATE_AND_PSEUDO_ATOMS, STANDARD_RESIDUE_BOND_ORDERS});
flags!(PdbWriteFlags, visit_write_flags, visitPdbWriteFlags,
    "export function visitPdbWriteFlags(v,f){f(v);}",
    {NONE, OMIT_SINGLE_CONECT, BIDIRECTIONAL_CONECT, NO_BOND_ORDER, MASTER, TER});
