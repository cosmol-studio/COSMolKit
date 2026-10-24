//! Typed property operations call the canonical exact-kind Rust accessors.
use crate::{
    Molecule,
    errors::{enum_error, exception},
    fields::Field,
    input,
    operation_errors::OperationError,
    storage::NativeField,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "OperationError | PropertyValueError | Utf8Error")]
    #[derive(Clone)]
    pub type AtomPropertyReadCause;
}
fn cause(e: &ck::AtomPropertyReadError) -> Option<AtomPropertyReadCause> {
    Some(
        match e {
            ck::AtomPropertyReadError::Operation(v) => OperationError::from_rust(v.clone()),
            ck::AtomPropertyReadError::Type(v) => {
                crate::model_errors::PropertyValueError::from_rust(v.clone())
            }
            ck::AtomPropertyReadError::Utf8(v) => crate::external_errors::Utf8Error::from_rust(*v),
        }
        .unchecked_into(),
    )
}
enum_error!(AtomPropertyReadError,AtomPropertyReadErrorKind,ck::AtomPropertyReadError,{Operation=>ck::AtomPropertyReadError::Operation(_),Type=>ck::AtomPropertyReadError::Type(_),Utf8=>ck::AtomPropertyReadError::Utf8(_)},{
 #[wasm_bindgen(getter)]pub fn source(&self)->AtomPropertyReadCause{cause(&self.inner).expect("all native variants carry a source")}
},source=AtomPropertyReadCause,cause=cause);
macro_rules! scalar_properties {($(($ty:ty,$read:ident,$read_js:ident,$with:ident,$with_js:ident,$set:ident,$set_js:ident)),+)=>{#[wasm_bindgen]impl Molecule{$(
 #[wasm_bindgen(js_name=$read_js)]pub fn $read(&self,atom:Ts<input::IndexInput>,key:Ts<input::TextInput>)->Result<Option<$ty>,JsValue>{self.0.$read(ck::AtomId::new(input::read(atom)?),&input::read::<String>(key)?).map_err(AtomPropertyReadError::from_rust)}
 #[wasm_bindgen(js_name=$with_js)]pub fn $with(&self,atom:Ts<input::IndexInput>,key:Ts<input::TextInput>,value:<$ty as Field>::Input)->Result<Self,JsValue>{self.0.$with(ck::AtomId::new(input::read(atom)?),&input::read::<String>(key)?,<$ty>::input(value)?).map(Self).map_err(OperationError::from_rust)}
 #[wasm_bindgen(js_name=$set_js)]pub fn $set(&mut self,atom:Ts<input::IndexInput>,key:Ts<input::TextInput>,value:<$ty as Field>::Input)->Result<(),JsValue>{self.0.$set(ck::AtomId::new(input::read(atom)?),&input::read::<String>(key)?,<$ty>::input(value)?).map_err(OperationError::from_rust)}
 )+}};}
scalar_properties!(
    (
        i32,
        atom_property_int,
        atomPropertyInt,
        with_atom_property_int,
        withAtomPropertyInt,
        set_atom_property_int_,
        setAtomPropertyInt_
    ),
    (
        u32,
        atom_property_uint,
        atomPropertyUint,
        with_atom_property_uint,
        withAtomPropertyUint,
        set_atom_property_uint_,
        setAtomPropertyUint_
    ),
    (
        f64,
        atom_property_double,
        atomPropertyDouble,
        with_atom_property_double,
        withAtomPropertyDouble,
        set_atom_property_double_,
        setAtomPropertyDouble_
    ),
    (
        bool,
        atom_property_bool,
        atomPropertyBool,
        with_atom_property_bool,
        withAtomPropertyBool,
        set_atom_property_bool_,
        setAtomPropertyBool_
    )
);
macro_rules! borrowed_properties {($(($ty:ty,$read:ident,$read_js:ident,$with:ident,$with_js:ident,$set:ident,$set_js:ident)),+)=>{#[wasm_bindgen]impl Molecule{$(
 #[wasm_bindgen(js_name=$read_js)]pub fn $read(&self,atom:Ts<input::IndexInput>,key:Ts<input::TextInput>)->Result<<Option<$ty> as Field>::Output,JsValue>{let v=self.0.$read(ck::AtomId::new(input::read(atom)?),&input::read::<String>(key)?).map_err(AtomPropertyReadError::from_rust)?;<Option<$ty>>::get(NativeField::owned(v))}
 #[wasm_bindgen(js_name=$with_js)]pub fn $with(&self,atom:Ts<input::IndexInput>,key:Ts<input::TextInput>,value:<$ty as Field>::Input)->Result<Self,JsValue>{self.0.$with(ck::AtomId::new(input::read(atom)?),&input::read::<String>(key)?,&<$ty>::input(value)?).map(Self).map_err(OperationError::from_rust)}
 #[wasm_bindgen(js_name=$set_js)]pub fn $set(&mut self,atom:Ts<input::IndexInput>,key:Ts<input::TextInput>,value:<$ty as Field>::Input)->Result<(),JsValue>{self.0.$set(ck::AtomId::new(input::read(atom)?),&input::read::<String>(key)?,&<$ty>::input(value)?).map_err(OperationError::from_rust)}
 )+}};}
borrowed_properties!(
    (
        String,
        atom_property_string,
        atomPropertyString,
        with_atom_property_string,
        withAtomPropertyString,
        set_atom_property_string_,
        setAtomPropertyString_
    ),
    (
        Vec<i32>,
        atom_property_int_vector,
        atomPropertyIntVector,
        with_atom_property_int_vector,
        withAtomPropertyIntVector,
        set_atom_property_int_vector_,
        setAtomPropertyIntVector_
    ),
    (
        Vec<String>,
        atom_property_string_vector,
        atomPropertyStringVector,
        with_atom_property_string_vector,
        withAtomPropertyStringVector,
        set_atom_property_string_vector_,
        setAtomPropertyStringVector_
    )
);
