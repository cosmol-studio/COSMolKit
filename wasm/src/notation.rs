//! CX notation options retain native flag and coordinate-selector values.
use crate::fields::Field;
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct CxSmilesFields(pub(crate) ck::CxSmilesFields);
crate::fields::value_field!(
    ck::CxSmilesFields,
    CxSmilesFields,
    CxSmilesFieldsInput,
    "CxSmilesFields"
);
#[wasm_bindgen]
impl CxSmilesFields {
    #[wasm_bindgen(getter,js_name=NONE)]
    pub fn none() -> Self {
        Self(ck::CxSmilesFields::NONE)
    }
    #[wasm_bindgen(getter,js_name=ATOM_LABELS)]
    pub fn atom_labels() -> Self {
        Self(ck::CxSmilesFields::ATOM_LABELS)
    }
    #[wasm_bindgen(getter,js_name=MOLFILE_VALUES)]
    pub fn molfile_values() -> Self {
        Self(ck::CxSmilesFields::MOLFILE_VALUES)
    }
    #[wasm_bindgen(getter,js_name=COORDS)]
    pub fn coords() -> Self {
        Self(ck::CxSmilesFields::COORDS)
    }
    #[wasm_bindgen(getter,js_name=RADICALS)]
    pub fn radicals() -> Self {
        Self(ck::CxSmilesFields::RADICALS)
    }
    #[wasm_bindgen(getter,js_name=ATOM_PROPS)]
    pub fn atom_props() -> Self {
        Self(ck::CxSmilesFields::ATOM_PROPS)
    }
    #[wasm_bindgen(getter,js_name=LINKNODES)]
    pub fn linknodes() -> Self {
        Self(ck::CxSmilesFields::LINKNODES)
    }
    #[wasm_bindgen(getter,js_name=ENHANCED_STEREO)]
    pub fn enhanced_stereo() -> Self {
        Self(ck::CxSmilesFields::ENHANCED_STEREO)
    }
    #[wasm_bindgen(getter,js_name=SGROUPS)]
    pub fn sgroups() -> Self {
        Self(ck::CxSmilesFields::SGROUPS)
    }
    #[wasm_bindgen(getter,js_name=POLYMER)]
    pub fn polymer() -> Self {
        Self(ck::CxSmilesFields::POLYMER)
    }
    #[wasm_bindgen(getter,js_name=BOND_CFG)]
    pub fn bond_cfg() -> Self {
        Self(ck::CxSmilesFields::BOND_CFG)
    }
    #[wasm_bindgen(getter,js_name=BOND_ATROPISOMER)]
    pub fn bond_atropisomer() -> Self {
        Self(ck::CxSmilesFields::BOND_ATROPISOMER)
    }
    #[wasm_bindgen(getter,js_name=COORDINATE_BONDS)]
    pub fn coordinate_bonds() -> Self {
        Self(ck::CxSmilesFields::COORDINATE_BONDS)
    }
    #[wasm_bindgen(getter,js_name=HYDROGEN_BONDS)]
    pub fn hydrogen_bonds() -> Self {
        Self(ck::CxSmilesFields::HYDROGEN_BONDS)
    }
    #[wasm_bindgen(getter,js_name=ZERO_BONDS)]
    pub fn zero_bonds() -> Self {
        Self(ck::CxSmilesFields::ZERO_BONDS)
    }
    #[wasm_bindgen(getter,js_name=ALL)]
    pub fn all() -> Self {
        Self(ck::CxSmilesFields::ALL)
    }
    #[wasm_bindgen(getter,js_name=ALL_BUT_COORDS)]
    pub fn all_but_coords() -> Self {
        Self(ck::CxSmilesFields::ALL_BUT_COORDS)
    }
    #[wasm_bindgen(getter)]
    pub fn bits(&self) -> u32 {
        self.0.bits()
    }
    pub fn contains(&self, other: CxSmilesFieldsInput) -> Result<bool, JsValue> {
        Ok(self.0.contains(ck::CxSmilesFields::input(other)?))
    }
    pub fn union(&self, other: CxSmilesFieldsInput) -> Result<Self, JsValue> {
        Ok(Self(self.0 | ck::CxSmilesFields::input(other)?))
    }
}
ck::__ck_configuration_CxSmilesWriteParams!(
    crate::project_parameters,
    ck::CxSmilesWriteParams,
    {},
    { CxSmilesWriteParams }
);
ck::__ck_configuration_FragmentCxSmilesWriteParams!(
    crate::project_parameters,
    ck::FragmentCxSmilesWriteParams,
    {},
    { FragmentCxSmilesWriteParams }
);
