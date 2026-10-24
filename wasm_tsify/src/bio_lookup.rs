//! Explicit alternate-location selection delegated to native structure queries.
use crate::{
    atoms::ElementInput,
    bio_hierarchy::{BioAtom, BioAtomArray, BioStructure},
    bio_read_errors::BioStructureError,
    bio_source_values::{AltLocLabelInput, AtomName, AtomNameInput},
    fields::Field,
    input,
    storage::NativeField,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
pub struct AltLocRequest(ck::AltLocRequest);
crate::fields::value_field!(
    ck::AltLocRequest,
    AltLocRequest,
    AltLocRequestInput,
    "AltLocRequest"
);
#[wasm_bindgen]
impl AltLocRequest {
    pub fn any() -> Self {
        Self(ck::AltLocRequest::Any)
    }
    pub fn exact(label: Option<AltLocLabelInput>) -> Result<Self, JsValue> {
        Ok(Self(ck::AltLocRequest::Exact(
            label.map(<ck::AltLocLabel as Field>::input).transpose()?,
        )))
    }
}
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "readonly [number, BioAtom]")]
    pub type BioAtomMatch;
}
fn found((id, atom): (ck::BioAtomId, ck::BioAtom)) -> BioAtomMatch {
    let array = js_sys::Array::new();
    array.push(&JsValue::from(id.value()));
    array.push(&BioAtom(atom).into());
    js_sys::Object::freeze(&array);
    array.unchecked_into()
}
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(js_name=residueAtoms)]
    pub fn residue_atoms(
        &self,
        residue: Ts<input::U32Input>,
    ) -> Result<Option<BioAtomArray>, JsValue> {
        self.0
            .residue_atoms(ck::BioResidueId::new(input::read::<u32>(residue)?))
            .map(|v| <Vec<ck::BioAtom> as Field>::get(NativeField::owned(v)))
            .transpose()
    }
    #[wasm_bindgen(js_name=findAtom)]
    pub fn find_atom(
        &self,
        residue_id: Ts<input::U32Input>,
        name: AtomNameInput,
        request: AltLocRequestInput,
        element: Option<ElementInput>,
    ) -> Result<Option<BioAtomMatch>, JsValue> {
        AtomName::with_ref(&name, |name| {
            AltLocRequest::with_ref(&request, |request| {
                Ok(self
                    .0
                    .find_atom(
                        ck::BioResidueId::new(input::read::<u32>(residue_id)?),
                        *name,
                        *request,
                        element.map(<ck::Element as Field>::input).transpose()?,
                    )
                    .map(found))
            })
        })
    }
    #[wasm_bindgen(js_name=atomByAltloc)]
    pub fn atom_by_altloc(
        &self,
        residue_id: Ts<input::U32Input>,
        name: AtomNameInput,
        altloc: Option<AltLocLabelInput>,
    ) -> Result<BioAtomMatch, JsValue> {
        AtomName::with_ref(&name, |name| {
            self.0
                .atom_by_altloc(
                    ck::BioResidueId::new(input::read::<u32>(residue_id)?),
                    *name,
                    altloc.map(<ck::AltLocLabel as Field>::input).transpose()?,
                )
                .map(found)
                .map_err(BioStructureError::from_rust)
        })
    }
}
