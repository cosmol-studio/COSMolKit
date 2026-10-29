//! Native IO configuration declarations and their typed value fields.
use crate::{fields::Field, input};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct MolCoordinateSelection(pub(crate) ck::MolCoordinateSelection);
crate::fields::value_field!(
    ck::MolCoordinateSelection,
    MolCoordinateSelection,
    MolCoordinateSelectionInput,
    "MolCoordinateSelection"
);
#[wasm_bindgen]
pub enum MolCoordinateSelectionKind {
    Auto,
    TwoD,
    ThreeD,
}
#[wasm_bindgen]
impl MolCoordinateSelection {
    pub fn auto() -> Self {
        Self(ck::MolCoordinateSelection::Auto)
    }
    #[wasm_bindgen(js_name=twoD)]
    pub fn two_d(id: tsify::Ts<<usize as input::CheckedValue>::Wire>) -> Result<Self, JsValue> {
        Ok(Self(ck::MolCoordinateSelection::TwoD {
            id: input::read::<usize>(id)?,
        }))
    }
    #[wasm_bindgen(js_name=threeD)]
    pub fn three_d(id: tsify::Ts<<usize as input::CheckedValue>::Wire>) -> Result<Self, JsValue> {
        Ok(Self(ck::MolCoordinateSelection::ThreeD {
            id: input::read::<usize>(id)?,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> MolCoordinateSelectionKind {
        match self.0 {
            ck::MolCoordinateSelection::Auto => MolCoordinateSelectionKind::Auto,
            ck::MolCoordinateSelection::TwoD { .. } => MolCoordinateSelectionKind::TwoD,
            ck::MolCoordinateSelection::ThreeD { .. } => MolCoordinateSelectionKind::ThreeD,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> Option<usize> {
        match self.0 {
            ck::MolCoordinateSelection::Auto => None,
            ck::MolCoordinateSelection::TwoD { id } | ck::MolCoordinateSelection::ThreeD { id } => {
                Some(id)
            }
        }
    }
}
#[wasm_bindgen]
pub struct PdbWriteFlags(pub(crate) ck::PdbWriteFlags);
crate::fields::value_field!(
    ck::PdbWriteFlags,
    PdbWriteFlags,
    PdbWriteFlagsInput,
    "PdbWriteFlags"
);
#[wasm_bindgen]
impl PdbWriteFlags {
    #[wasm_bindgen(getter,js_name=NONE)]
    pub fn none() -> Self {
        Self(ck::PdbWriteFlags::NONE)
    }
    #[wasm_bindgen(getter,js_name=OMIT_SINGLE_CONECT)]
    pub fn omit_single_conect() -> Self {
        Self(ck::PdbWriteFlags::OMIT_SINGLE_CONECT)
    }
    #[wasm_bindgen(getter,js_name=BIDIRECTIONAL_CONECT)]
    pub fn bidirectional_conect() -> Self {
        Self(ck::PdbWriteFlags::BIDIRECTIONAL_CONECT)
    }
    #[wasm_bindgen(getter,js_name=NO_BOND_ORDER)]
    pub fn no_bond_order() -> Self {
        Self(ck::PdbWriteFlags::NO_BOND_ORDER)
    }
    #[wasm_bindgen(getter,js_name=MASTER)]
    pub fn master() -> Self {
        Self(ck::PdbWriteFlags::MASTER)
    }
    #[wasm_bindgen(getter,js_name=TER)]
    pub fn ter() -> Self {
        Self(ck::PdbWriteFlags::TER)
    }
    #[wasm_bindgen(js_name=fromBits)]
    pub fn from_bits(bits: tsify::Ts<input::U32Input>) -> Result<Self, JsValue> {
        Ok(Self(ck::PdbWriteFlags::from_bits(input::read::<u32>(
            bits,
        )?)))
    }
    #[wasm_bindgen(getter)]
    pub fn bits(&self) -> u32 {
        self.0.bits()
    }
    #[wasm_bindgen(getter,js_name=isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    pub fn contains(&self, other: PdbWriteFlagsInput) -> Result<bool, JsValue> {
        Ok(self.0.contains(<ck::PdbWriteFlags as Field>::input(other)?))
    }
    pub fn union(&self, other: PdbWriteFlagsInput) -> Result<Self, JsValue> {
        Ok(Self(self.0 | <ck::PdbWriteFlags as Field>::input(other)?))
    }
}
ck::__ck_configuration_MolBlockWriteParams!(
    crate::project_parameters,
    ck::MolBlockWriteParams,
    {},
    { MolBlockWriteParams }
);
ck::__ck_configuration_PdbWriteParams!(crate::project_parameters, ck::PdbWriteParams, {}, {
    PdbWriteParams
});
ck::__ck_configuration_Mol2ReadParams!(crate::project_parameters, ck::Mol2ReadParams, {}, {
    Mol2ReadParams
});
ck::__ck_configuration_XyzWriteParams!(crate::project_parameters, ck::XyzWriteParams, {}, {
    XyzWriteParams
});
