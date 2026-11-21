//! Native reaction template values and metadata-generated options.
pub(crate) use crate::reaction_validation_values::*;
use crate::{fields::Field, input};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub enum ReactionCoordinateSelectionKind {
    Auto,
    TwoD,
    ThreeD,
}
#[wasm_bindgen]
pub struct ReactionCoordinateSelection(pub(crate) ck::ReactionCoordinateSelection);
crate::fields::value_field!(
    ck::ReactionCoordinateSelection,
    ReactionCoordinateSelection,
    ReactionCoordinateSelectionInput,
    "ReactionCoordinateSelection"
);
crate::fields::value_array!(
    ck::ReactionCoordinateSelection,
    ReactionCoordinateSelectionArray,
    "ReadonlyArray<ReactionCoordinateSelection>"
);
#[wasm_bindgen]
impl ReactionCoordinateSelection {
    #[wasm_bindgen(js_name=auto)]
    pub fn auto() -> Self {
        Self(ck::ReactionCoordinateSelection::Auto)
    }
    #[wasm_bindgen(js_name=twoD)]
    pub fn two_d(id: Ts<input::IndexInput>) -> Result<Self, JsValue> {
        Ok(Self(ck::ReactionCoordinateSelection::TwoD {
            id: input::read::<usize>(id)?,
        }))
    }
    #[wasm_bindgen(js_name=threeD)]
    pub fn three_d(id: Ts<input::IndexInput>) -> Result<Self, JsValue> {
        Ok(Self(ck::ReactionCoordinateSelection::ThreeD {
            id: input::read::<usize>(id)?,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> Option<usize> {
        match self.0 {
            ck::ReactionCoordinateSelection::Auto => None,
            ck::ReactionCoordinateSelection::TwoD { id }
            | ck::ReactionCoordinateSelection::ThreeD { id } => Some(id),
        }
    }
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> ReactionCoordinateSelectionKind {
        match self.0 {
            ck::ReactionCoordinateSelection::Auto => ReactionCoordinateSelectionKind::Auto,
            ck::ReactionCoordinateSelection::TwoD { .. } => ReactionCoordinateSelectionKind::TwoD,
            ck::ReactionCoordinateSelection::ThreeD { .. } => {
                ReactionCoordinateSelectionKind::ThreeD
            }
        }
    }
}
ck::__ck_configuration_ReactionParseParams!(
    crate::project_parameters,
    ck::ReactionParseParams,
    {},
    { ReactionParseParams }
);
ck::__ck_configuration_ReactionWriteParams!(
    crate::project_parameters,
    ck::ReactionWriteParams,
    {},
    { ReactionWriteParams }
);
ck::__ck_configuration_ReactionValidationParams!(
    crate::project_parameters,
    ck::ReactionValidationParams,
    {},
    { ReactionValidationParams }
);
ck::__ck_configuration_ReactionTemplateRemovalParams!(
    crate::project_parameters,
    ck::ReactionTemplateRemovalParams,
    {},
    { ReactionTemplateRemovalParams }
);
ck::__ck_configuration_ReactionRunParams!(crate::project_parameters, ck::ReactionRunParams, {}, {
    ReactionRunParams
});
ck::__ck_configuration_ReactionSingleRunParams!(
    crate::project_parameters,
    ck::ReactionSingleRunParams,
    {},
    { ReactionSingleRunParams }
);
ck::__ck_configuration_ReactionApplyParams!(
    crate::project_parameters,
    ck::ReactionApplyParams,
    {},
    { ReactionApplyParams }
);
