//! Canonical reaction values and operations, without reaction algorithms.
use crate::alignment_values::operation_error;
use crate::query_values::QueryGraph;
use crate::reaction_errors::*;
use crate::reaction_parameters::*;
use crate::search::SubstructMatchParams;
use crate::{Molecule, host_values::*};
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use std::{cell::RefCell, sync::Arc};
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum ReactionRole {
    Reactant = 0,
    Product = 1,
    Agent = 2,
}
impl From<ck::ReactionRole> for ReactionRole {
    fn from(v: ck::ReactionRole) -> Self {
        match v {
            ck::ReactionRole::Reactant => Self::Reactant,
            ck::ReactionRole::Product => Self::Product,
            ck::ReactionRole::Agent => Self::Agent,
        }
    }
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum ReactionValidationSeverity {
    Warning = 0,
    Error = 1,
}
impl From<ck::ReactionValidationSeverity> for ReactionValidationSeverity {
    fn from(v: ck::ReactionValidationSeverity) -> Self {
        match v {
            ck::ReactionValidationSeverity::Warning => Self::Warning,
            ck::ReactionValidationSeverity::Error => Self::Error,
        }
    }
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum ReactionValidationIssueKind {
    MissingReactants = 0,
    MissingProducts = 1,
    DuplicateReactantMap = 2,
    UnmappedReactant = 3,
    DuplicateProductMap = 4,
    MissingProductMapReactant = 5,
    UnmappedProduct = 6,
    UnmappedReactantMaps = 7,
    MultipleCharge = 8,
    MultipleHydrogenCount = 9,
    MultipleMass = 10,
    MultipleIsotope = 11,
}
impl From<ck::ReactionValidationIssueKind> for ReactionValidationIssueKind {
    fn from(v: ck::ReactionValidationIssueKind) -> Self {
        match v {
            ck::ReactionValidationIssueKind::MissingReactants => Self::MissingReactants,
            ck::ReactionValidationIssueKind::MissingProducts => Self::MissingProducts,
            ck::ReactionValidationIssueKind::DuplicateReactantMap => Self::DuplicateReactantMap,
            ck::ReactionValidationIssueKind::UnmappedReactant => Self::UnmappedReactant,
            ck::ReactionValidationIssueKind::DuplicateProductMap => Self::DuplicateProductMap,
            ck::ReactionValidationIssueKind::MissingProductMapReactant => {
                Self::MissingProductMapReactant
            }
            ck::ReactionValidationIssueKind::UnmappedProduct => Self::UnmappedProduct,
            ck::ReactionValidationIssueKind::UnmappedReactantMaps => Self::UnmappedReactantMaps,
            ck::ReactionValidationIssueKind::MultipleCharge => Self::MultipleCharge,
            ck::ReactionValidationIssueKind::MultipleHydrogenCount => Self::MultipleHydrogenCount,
            ck::ReactionValidationIssueKind::MultipleMass => Self::MultipleMass,
            ck::ReactionValidationIssueKind::MultipleIsotope => Self::MultipleIsotope,
        }
    }
}
#[wasm_bindgen]
pub struct Reaction {
    inner: RefCell<ck::Reaction>,
}
#[wasm_bindgen]
impl Reaction {
    #[wasm_bindgen(unchecked_return_type = "Molecule[][]")]
    pub fn run(
        &self,
        #[wasm_bindgen(unchecked_param_type = "Molecule[]")] reactants: JsValue,
        params: &ReactionRunParams,
    ) -> Result<Array, JsValue> {
        let inputs = molecules(&reactants)?;
        let inputs: Vec<_> = inputs.iter().map(|m| &**m).collect();
        cosmolkit_wasm::reaction_run(&mut self.inner.borrow_mut(), &inputs, &params.inner)
            .map(sets).map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            inner: ck::Reaction::new().into(),
        }
    }
    #[wasm_bindgen(js_name=fromSmirks)]
    pub fn from_smirks(text: &str) -> Result<Self, JsValue> {
        ck::Reaction::from_smirks(text)
            .map(|inner| Self {
                inner: inner.into(),
            })
            .map_err(|e| parse_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromSmirksWithParams)]
    pub fn from_smirks_with_params(
        text: &str,
        params: &ReactionParseParams,
    ) -> Result<Self, JsValue> {
        ck::Reaction::from_smirks_with_params(text, &params.inner)
            .map(|inner| Self {
                inner: inner.into(),
            })
            .map_err(|e| parse_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromTemplates)]
    pub fn from_templates(
        #[wasm_bindgen(unchecked_param_type = "QueryGraph[]")] reactants: JsValue,
        #[wasm_bindgen(unchecked_param_type = "QueryGraph[]")] products: JsValue,
        #[wasm_bindgen(unchecked_param_type = "QueryGraph[]")] agents: JsValue,
    ) -> Result<Self, JsValue> {
        ck::Reaction::from_templates(queries(&reactants)?, queries(&products)?, queries(&agents)?)
            .map(|inner| Self {
                inner: inner.into(),
            })
            .map_err(|e| model_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numReactantTemplates)]
    pub fn num_reactant_templates(&self) -> usize {
        self.inner.borrow().num_reactant_templates()
    }
    #[wasm_bindgen(js_name=reactantTemplates,unchecked_return_type="QueryGraph[]")]
    pub fn reactant_templates(&self) -> Array {
        self.inner
            .borrow()
            .reactant_templates()
            .iter()
            .cloned()
            .map(|inner| JsValue::from(QueryGraph { inner }))
            .collect()
    }
    #[wasm_bindgen(js_name=reactantTemplate)]
    pub fn reactant_template(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
    ) -> Result<QueryGraph, JsValue> {
        self.inner
            .borrow()
            .reactant_template(usize_value(&index, "index")?)
            .map(|inner| QueryGraph {
                inner: inner.clone(),
            })
            .map_err(|e| model_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withReactantTemplate)]
    pub fn with_reactant_template(&self, template: &QueryGraph) -> Result<Self, JsValue> {
        self.inner
            .borrow()
            .with_reactant_template(template.inner.clone())
            .map(|inner| Self {
                inner: inner.into(),
            })
            .map_err(|e| model_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numProductTemplates)]
    pub fn num_product_templates(&self) -> usize {
        self.inner.borrow().num_product_templates()
    }
    #[wasm_bindgen(js_name=productTemplates,unchecked_return_type="QueryGraph[]")]
    pub fn product_templates(&self) -> Array {
        self.inner
            .borrow()
            .product_templates()
            .iter()
            .cloned()
            .map(|inner| JsValue::from(QueryGraph { inner }))
            .collect()
    }
    #[wasm_bindgen(js_name=productTemplate)]
    pub fn product_template(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
    ) -> Result<QueryGraph, JsValue> {
        self.inner
            .borrow()
            .product_template(usize_value(&index, "index")?)
            .map(|inner| QueryGraph {
                inner: inner.clone(),
            })
            .map_err(|e| model_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withProductTemplate)]
    pub fn with_product_template(&self, template: &QueryGraph) -> Result<Self, JsValue> {
        self.inner
            .borrow()
            .with_product_template(template.inner.clone())
            .map(|inner| Self {
                inner: inner.into(),
            })
            .map_err(|e| model_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numAgentTemplates)]
    pub fn num_agent_templates(&self) -> usize {
        self.inner.borrow().num_agent_templates()
    }
    #[wasm_bindgen(js_name=agentTemplates,unchecked_return_type="QueryGraph[]")]
    pub fn agent_templates(&self) -> Array {
        self.inner
            .borrow()
            .agent_templates()
            .iter()
            .cloned()
            .map(|inner| JsValue::from(QueryGraph { inner }))
            .collect()
    }
    #[wasm_bindgen(js_name=agentTemplate)]
    pub fn agent_template(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
    ) -> Result<QueryGraph, JsValue> {
        self.inner
            .borrow()
            .agent_template(usize_value(&index, "index")?)
            .map(|inner| QueryGraph {
                inner: inner.clone(),
            })
            .map_err(|e| model_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withAgentTemplate)]
    pub fn with_agent_template(&self, template: &QueryGraph) -> Result<Self, JsValue> {
        self.inner
            .borrow()
            .with_agent_template(template.inner.clone())
            .map(|inner| Self {
                inner: inner.into(),
            })
            .map_err(|e| model_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=isInitialized)]
    pub fn is_initialized(&self) -> bool {
        self.inner.borrow().is_initialized()
    }
    #[wasm_bindgen(js_name=implicitProperties)]
    pub fn implicit_properties(&self) -> bool {
        self.inner.borrow().implicit_properties()
    }
    #[wasm_bindgen(js_name=withImplicitProperties)]
    pub fn with_implicit_properties(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] enabled: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: self
                .inner
                .borrow()
                .with_implicit_properties(bool_value(&enabled, "enabled")?)
                .into(),
        })
    }
    #[wasm_bindgen(js_name=matchParams)]
    pub fn match_params(&self) -> SubstructMatchParams {
        SubstructMatchParams {
            inner: self.inner.borrow().match_params().clone(),
        }
    }
    #[wasm_bindgen(js_name=withMatchParams)]
    pub fn with_match_params(&self, params: &SubstructMatchParams) -> Result<Self, JsValue> {
        self.inner
            .borrow()
            .with_match_params(&params.inner)
            .map(|inner| Self {
                inner: inner.into(),
            })
            .map_err(|e| model_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toSmirks)]
    pub fn to_smirks(&self) -> Result<String, JsValue> {
        let value = self
            .inner
            .borrow()
            .to_smirks()
            .map_err(|e| write_error(&e).unwrap_or_else(|e| e))?;
        text(&value)
    }
    #[wasm_bindgen(js_name=toSmirksWithParams)]
    pub fn to_smirks_with_params(&self, params: &ReactionWriteParams) -> Result<String, JsValue> {
        let value = self
            .inner
            .borrow()
            .to_smirks_with_params(&params.inner)
            .map_err(|e| write_error(&e).unwrap_or_else(|e| e))?;
        text(&value)
    }
    #[wasm_bindgen(js_name=toCxSmirks)]
    pub fn to_cx_smirks(&self) -> Result<String, JsValue> {
        let value = self
            .inner
            .borrow()
            .to_cx_smirks()
            .map_err(|e| write_error(&e).unwrap_or_else(|e| e))?;
        text(&value)
    }
    #[wasm_bindgen(js_name=toCxSmirksWithParams)]
    pub fn to_cx_smirks_with_params(
        &self,
        params: &ReactionWriteParams,
    ) -> Result<String, JsValue> {
        let value = self
            .inner
            .borrow()
            .to_cx_smirks_with_params(&params.inner)
            .map_err(|e| write_error(&e).unwrap_or_else(|e| e))?;
        text(&value)
    }
    #[wasm_bindgen(js_name=validate)]
    pub fn validate(&self) -> Result<ReactionValidationReport, JsValue> {
        self.inner
            .borrow()
            .validate()
            .map(|inner| ReactionValidationReport { inner })
            .map_err(|e| validation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=validateWithParams)]
    pub fn validate_with_params(
        &self,
        params: &ReactionValidationParams,
    ) -> Result<ReactionValidationReport, JsValue> {
        self.inner
            .borrow()
            .validate_with_params(&params.inner)
            .map(|inner| ReactionValidationReport { inner })
            .map_err(|e| validation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withInitialized)]
    pub fn with_initialized(&self) -> Result<Reaction, JsValue> {
        self.inner
            .borrow()
            .with_initialized()
            .map(|inner| Reaction {
                inner: inner.into(),
            })
            .map_err(|e| initialization_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withInitializedWithParams)]
    pub fn with_initialized_with_params(
        &self,
        params: &ReactionValidationParams,
    ) -> Result<Reaction, JsValue> {
        self.inner
            .borrow()
            .with_initialized_with_params(&params.inner)
            .map(|inner| Reaction {
                inner: inner.into(),
            })
            .map_err(|e| initialization_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withoutAgents)]
    pub fn without_agents(&self) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: self.inner.borrow().without_agents(),
        }
    }
    #[wasm_bindgen(js_name=withoutUnmappedReactants)]
    pub fn without_unmapped_reactants(&self) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: self.inner.borrow().without_unmapped_reactants(),
        }
    }
    #[wasm_bindgen(js_name=withoutUnmappedReactantsWithParams)]
    pub fn without_unmapped_reactants_with_params(
        &self,
        params: &ReactionTemplateRemovalParams,
    ) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: self
                .inner
                .borrow()
                .without_unmapped_reactants_with_params(&params.inner),
        }
    }
    #[wasm_bindgen(js_name=withoutUnmappedProducts)]
    pub fn without_unmapped_products(&self) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: self.inner.borrow().without_unmapped_products(),
        }
    }
    #[wasm_bindgen(js_name=withoutUnmappedProductsWithParams)]
    pub fn without_unmapped_products_with_params(
        &self,
        params: &ReactionTemplateRemovalParams,
    ) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: self
                .inner
                .borrow()
                .without_unmapped_products_with_params(&params.inner),
        }
    }
}
#[wasm_bindgen(
    inline_js = "export function visitReactionQuery(v,f){try{f(v);}catch(cause){throw new TypeError('invalid QueryGraph',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitReactionQuery)]
    fn visitReactionQuery(v: &JsValue, f: &mut dyn FnMut(&QueryGraph)) -> Result<(), JsValue>;
}
#[wasm_bindgen(
    inline_js = "export function visitReactionMolecule(v,f){try{f(v);}catch(cause){throw new TypeError('invalid Molecule',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitReactionMolecule)]
    fn visitReactionMolecule(v: &JsValue, f: &mut dyn FnMut(&Molecule)) -> Result<(), JsValue>;
}
fn queries(value: &JsValue) -> Result<Vec<ck::QueryGraph>, JsValue> {
    sequence(value, "templates")?
        .iter()
        .map(|v| {
            let mut out = None;
            visitReactionQuery(&v, &mut |q: &QueryGraph| out = Some(q.inner.clone()))?;
            out.ok_or_else(|| type_error("template"))
        })
        .collect()
}
fn molecules(value: &JsValue) -> Result<Vec<Arc<cosmolkit_wasm::Molecule>>, JsValue> {
    sequence(value, "reactants")?
        .iter()
        .map(|v| {
            let mut out = None;
            visitReactionMolecule(&v, &mut |m: &Molecule| out = Some(m.inner.clone()))?;
            out.ok_or_else(|| type_error("reactant"))
        })
        .collect()
}
fn sets(values: Vec<Vec<cosmolkit_wasm::Molecule>>) -> Array {
    values
        .into_iter()
        .map(|set| {
            JsValue::from(
                set.into_iter()
                    .map(|inner| {
                        JsValue::from(Molecule {
                            inner: Arc::new(inner),
                        })
                    })
                    .collect::<Array>(),
            )
        })
        .collect()
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=reactionProducts,unchecked_return_type="Molecule[][]")]
    pub fn reaction_products(
        &self,
        reaction: &Reaction,
        #[wasm_bindgen(unchecked_param_type = "number")] reactant_template: JsValue,
    ) -> Result<Array, JsValue> {
        self.inner
            .reaction_products(
                &mut reaction.inner.borrow_mut(),
                usize_value(&reactant_template, "reactantTemplate")?,
            )
            .map(sets)
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=reactionProductsWithParams,unchecked_return_type="Molecule[][]")]
    pub fn reaction_products_with_params(
        &self,
        reaction: &Reaction,
        #[wasm_bindgen(unchecked_param_type = "number")] reactant_template: JsValue,
        params: &ReactionSingleRunParams,
    ) -> Result<Array, JsValue> {
        self.inner
            .reaction_products_with_params(
                &mut reaction.inner.borrow_mut(),
                usize_value(&reactant_template, "reactantTemplate")?,
                &params.inner,
            )
            .map(sets)
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=reactionProductsFromInputs,unchecked_return_type="Molecule[][]")]
    pub fn reaction_products_from_inputs(
        &self,
        reaction: &Reaction,
        #[wasm_bindgen(unchecked_param_type = "Molecule[]")] reactants: JsValue,
        params: &ReactionRunParams,
    ) -> Result<Array, JsValue> {
        let inputs = molecules(&reactants)?;
        let inputs: Vec<_> = inputs.iter().map(|m| &**m).collect();
        self.inner
            .reaction_products_from_inputs(&mut reaction.inner.borrow_mut(), &inputs, &params.inner)
            .map(sets)
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=applyReaction)]
    pub fn apply_reaction(&self, reaction: &Reaction) -> Result<ReactionApplyResult, JsValue> {
        self.inner
            .apply_reaction(&mut reaction.inner.borrow_mut())
            .map(|inner| ReactionApplyResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=applyReactionWithParams)]
    pub fn apply_reaction_with_params(
        &self,
        reaction: &Reaction,
        params: &ReactionApplyParams,
    ) -> Result<ReactionApplyResult, JsValue> {
        self.inner
            .apply_reaction_with_params(&mut reaction.inner.borrow_mut(), &params.inner)
            .map(|inner| ReactionApplyResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=applyReaction_)]
    pub fn apply_reaction_(&self, reaction: &Reaction) -> Result<bool, JsValue> {
        self.inner
            .apply_reaction_(&mut reaction.inner.borrow_mut())
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=applyReactionWithParams_)]
    pub fn apply_reaction_with_params_(
        &self,
        reaction: &Reaction,
        params: &ReactionApplyParams,
    ) -> Result<bool, JsValue> {
        self.inner
            .apply_reaction_with_params_(&mut reaction.inner.borrow_mut(), &params.inner)
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen]
pub struct ReactionTemplateRemoval {
    pub(crate) inner: ck::ReactionTemplateRemoval,
}
#[wasm_bindgen]
impl ReactionTemplateRemoval {
    #[wasm_bindgen(getter,js_name=reaction)]
    pub fn reaction(&self) -> Reaction {
        Reaction {
            inner: self.inner.reaction().clone().into(),
        }
    }
    #[wasm_bindgen(getter,js_name=removedTemplates,unchecked_return_type="QueryGraph[]")]
    pub fn removed_templates(&self) -> Array {
        self.inner
            .removed_templates()
            .iter()
            .cloned()
            .map(|inner| JsValue::from(QueryGraph { inner }))
            .collect()
    }
}
#[wasm_bindgen]
pub struct ReactionApplyResult {
    pub(crate) inner: cosmolkit_wasm::ReactionApplyResult,
}
#[wasm_bindgen]
impl ReactionApplyResult {
    #[wasm_bindgen(getter,js_name=molecule)]
    pub fn molecule(&self) -> Molecule {
        Molecule {
            inner: Arc::new(self.inner.molecule()),
        }
    }
    #[wasm_bindgen(getter,js_name=changed)]
    pub fn changed(&self) -> bool {
        self.inner.changed()
    }
}
#[wasm_bindgen]
pub struct ReactionValidationReport {
    pub(crate) inner: ck::ReactionValidationReport,
}
#[wasm_bindgen]
impl ReactionValidationReport {
    #[wasm_bindgen(getter,js_name=warnings,unchecked_return_type="ReactionValidationIssue[]")]
    pub fn warnings(&self) -> Array {
        self.inner
            .warnings()
            .iter()
            .cloned()
            .map(|inner| JsValue::from(ReactionValidationIssue { inner }))
            .collect()
    }
    #[wasm_bindgen(getter,js_name=errors,unchecked_return_type="ReactionValidationIssue[]")]
    pub fn errors(&self) -> Array {
        self.inner
            .errors()
            .iter()
            .cloned()
            .map(|inner| JsValue::from(ReactionValidationIssue { inner }))
            .collect()
    }
    #[wasm_bindgen(getter,js_name=isValid)]
    pub fn is_valid(&self) -> bool {
        self.inner.is_valid()
    }
    #[wasm_bindgen(getter,js_name=numWarnings)]
    pub fn num_warnings(&self) -> usize {
        self.inner.num_warnings()
    }
    #[wasm_bindgen(getter,js_name=numErrors)]
    pub fn num_errors(&self) -> usize {
        self.inner.num_errors()
    }
}
#[wasm_bindgen]
pub struct ReactionValidationIssue {
    pub(crate) inner: ck::ReactionValidationIssue,
}
#[wasm_bindgen]
impl ReactionValidationIssue {
    #[wasm_bindgen(getter,js_name=kind)]
    pub fn kind(&self) -> ReactionValidationIssueKind {
        self.inner.kind().into()
    }
    #[wasm_bindgen(getter,js_name=severity)]
    pub fn severity(&self) -> ReactionValidationSeverity {
        self.inner.severity().into()
    }
    #[wasm_bindgen(getter,js_name=role,unchecked_return_type="ReactionRole | null")]
    pub fn role(&self) -> JsValue {
        self.inner.role().map_or(JsValue::NULL, |v| {
            JsValue::from(ReactionRole::from(v) as u32)
        })
    }
    #[wasm_bindgen(getter,js_name=template,unchecked_return_type="number | null")]
    pub fn template(&self) -> JsValue {
        self.inner
            .template()
            .map_or(JsValue::NULL, |v| JsValue::from(v as f64))
    }
    #[wasm_bindgen(getter,js_name=atom,unchecked_return_type="number | null")]
    pub fn atom(&self) -> JsValue {
        self.inner
            .atom()
            .map_or(JsValue::NULL, |v| JsValue::from(v.index() as f64))
    }
    #[wasm_bindgen(getter,js_name=map,unchecked_return_type="number | null")]
    pub fn map(&self) -> JsValue {
        self.inner.map().map_or(JsValue::NULL, JsValue::from)
    }
    #[wasm_bindgen(getter,js_name=maps)]
    pub fn maps(&self) -> Vec<i32> {
        self.inner.maps().to_vec()
    }
    #[wasm_bindgen(getter,js_name=detail)]
    pub fn detail(&self) -> String {
        self.inner.detail().to_owned()
    }
}
#[wasm_bindgen(js_name=parseSmirks)]
pub fn parse_smirks(text: &str) -> Result<Reaction, JsValue> {
    Reaction::from_smirks(text)
}
#[wasm_bindgen(js_name=parseSmirksWithParams)]
pub fn parse_smirks_with_params(
    text: &str,
    params: &ReactionParseParams,
) -> Result<Reaction, JsValue> {
    Reaction::from_smirks_with_params(text, params)
}
