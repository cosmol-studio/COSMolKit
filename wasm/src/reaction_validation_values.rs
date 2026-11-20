use crate::{fields::Field, input};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;
crate::enums::enum_value!(ReactionRole,"ReactionRole"{Reactant="reactant",Product="product",Agent="agent"});
crate::enums::enum_value!(ReactionValidationSeverity,"ReactionValidationSeverity"{Warning="warning",Error="error"});
crate::enums::enum_value!(ReactionValidationIssueKind,"ReactionValidationIssueKind"{MissingReactants="missing_reactants",MissingProducts="missing_products",DuplicateReactantMap="duplicate_reactant_map",UnmappedReactant="unmapped_reactant",DuplicateProductMap="duplicate_product_map",MissingProductMapReactant="missing_product_map_reactant",UnmappedProduct="unmapped_product",UnmappedReactantMaps="unmapped_reactant_maps",MultipleCharge="multiple_charge",MultipleHydrogenCount="multiple_hydrogen_count",MultipleMass="multiple_mass",MultipleIsotope="multiple_isotope"});
#[wasm_bindgen]
pub struct ReactionValidationIssue(pub(crate) ck::ReactionValidationIssue);
crate::fields::value_field!(
    ck::ReactionValidationIssue,
    ReactionValidationIssue,
    ReactionValidationIssueInput,
    "ReactionValidationIssue"
);
crate::fields::value_array!(
    ck::ReactionValidationIssue,
    ReactionValidationIssues,
    "ReadonlyArray<ReactionValidationIssue>"
);
#[wasm_bindgen]
impl ReactionValidationIssue {
    #[wasm_bindgen(constructor)]
    pub fn new(
        kind: ReactionValidationIssueKindInput,
        severity: ReactionValidationSeverityInput,
        role: Option<ReactionRoleInput>,
        template: Option<Ts<input::IndexInput>>,
        atom: Option<Ts<input::IndexInput>>,
        map: Option<Ts<input::I32Input>>,
        maps: Ts<crate::fields::I32Array>,
        detail: Ts<input::TextInput>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::ReactionValidationIssue {
            kind: ck::ReactionValidationIssueKind::input(kind)?,
            severity: ck::ReactionValidationSeverity::input(severity)?,
            role: role.map(ck::ReactionRole::input).transpose()?,
            template: template.map(input::read::<usize>).transpose()?,
            atom: atom
                .map(input::read::<usize>)
                .transpose()?
                .map(ck::AtomId::new),
            map: map.map(input::read::<i32>).transpose()?,
            maps: Vec::<i32>::input(maps)?,
            detail: input::read::<String>(detail)?,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> ReactionValidationIssueKind {
        self.0.kind.into()
    }
    #[wasm_bindgen(getter)]
    pub fn severity(&self) -> ReactionValidationSeverity {
        self.0.severity.into()
    }
    #[wasm_bindgen(getter)]
    pub fn role(&self) -> Option<ReactionRole> {
        self.0.role.map(Into::into)
    }
    #[wasm_bindgen(getter)]
    pub fn template(&self) -> Option<usize> {
        self.0.template
    }
    #[wasm_bindgen(getter)]
    pub fn atom(&self) -> Option<usize> {
        self.0.atom.map(|id| id.index())
    }
    #[wasm_bindgen(getter)]
    pub fn map(&self) -> Option<i32> {
        self.0.map
    }
    #[wasm_bindgen(getter)]
    pub fn maps(&self) -> Vec<i32> {
        self.0.maps.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn detail(&self) -> String {
        self.0.detail.clone()
    }
}
#[wasm_bindgen]
pub struct ReactionValidationReport(pub(crate) ck::ReactionValidationReport);
#[wasm_bindgen]
impl ReactionValidationReport {
    #[wasm_bindgen(constructor)]
    pub fn new(
        warnings: ReactionValidationIssues,
        errors: ReactionValidationIssues,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::ReactionValidationReport {
            warnings: Vec::<ck::ReactionValidationIssue>::input(warnings)?,
            errors: Vec::<ck::ReactionValidationIssue>::input(errors)?,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn warnings(&self) -> Vec<ReactionValidationIssue> {
        self.0
            .warnings
            .iter()
            .cloned()
            .map(ReactionValidationIssue)
            .collect()
    }
    #[wasm_bindgen(getter)]
    pub fn errors(&self) -> Vec<ReactionValidationIssue> {
        self.0
            .errors
            .iter()
            .cloned()
            .map(ReactionValidationIssue)
            .collect()
    }
    #[wasm_bindgen(getter,js_name=numWarnings)]
    pub fn num_warnings(&self) -> usize {
        self.0.num_warnings()
    }
    #[wasm_bindgen(getter,js_name=numErrors)]
    pub fn num_errors(&self) -> usize {
        self.0.num_errors()
    }
    #[wasm_bindgen(getter,js_name=isValid)]
    pub fn is_valid(&self) -> bool {
        self.0.is_valid()
    }
}
