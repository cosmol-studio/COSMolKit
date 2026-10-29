//! Reaction execution retains the actual native mutable reaction state.
use crate::fragments::Molecules;
use crate::{
    Molecule,
    fields::Field,
    input,
    operation_errors::OperationError,
    reaction::{Reaction, ReactionInput},
    reaction_values::*,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "ReadonlyArray<ReadonlyArray<Molecule>>")]
    pub type ProductSets;
}
fn sets(values: Vec<Vec<ck::Molecule>>) -> ProductSets {
    let result = js_sys::Array::new();
    for value in values {
        let row = value
            .into_iter()
            .map(|m| JsValue::from(Molecule(m)))
            .collect::<js_sys::Array>();
        js_sys::Object::freeze(&row);
        result.push(&row);
    }
    js_sys::Object::freeze(&result).unchecked_into()
}
#[wasm_bindgen]
pub struct ReactionApplyResult(ck::ReactionApplyResult);
#[wasm_bindgen]
impl ReactionApplyResult {
    #[wasm_bindgen(constructor)]
    pub fn new(
        molecule: crate::MoleculeInput,
        changed: Ts<input::BoolInput>,
    ) -> Result<Self, JsValue> {
        Ok(Self(
            (
                ck::Molecule::input(molecule)?,
                input::read::<bool>(changed)?,
            )
                .into(),
        ))
    }
    #[wasm_bindgen(getter)]
    pub fn molecule(&self) -> Molecule {
        Molecule(self.0.molecule().clone())
    }
    #[wasm_bindgen(getter)]
    pub fn changed(&self) -> bool {
        self.0.changed()
    }
}
#[wasm_bindgen]
impl Reaction {
    #[wasm_bindgen]
    pub fn run(
        &mut self,
        reactants: Molecules,
        params: Option<ReactionRunParamsInput>,
    ) -> Result<ProductSets, JsValue> {
        let reactants = Vec::<ck::Molecule>::input(reactants)?;
        let refs = reactants.iter().collect::<Vec<_>>();
        ReactionRunParams::with_input(params.as_ref(), |p| {
            self.0.write(|r| {
                let targets = r
                    .match_params()
                    .extra_final_check
                    .as_ref()
                    .map(|_| reactants.clone());
                crate::callback::run_with_targets(targets, || {
                    r.run(&refs, p).map(sets).map_err(OperationError::from_rust)
                })
            })?
        })
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=reactionProducts)]
    pub fn reaction_products(
        &self,
        reaction: ReactionInput,
        reactant_template: Ts<input::IndexInput>,
        params: Option<ReactionSingleRunParamsInput>,
    ) -> Result<ProductSets, JsValue> {
        let index = input::read::<usize>(reactant_template)?;
        Reaction::with_mut(&reaction, |r| {
            ReactionSingleRunParams::with_input(params.as_ref(), |p| {
                let targets = r
                    .match_params()
                    .extra_final_check
                    .as_ref()
                    .map(|_| vec![self.0.clone()]);
                crate::callback::run_with_targets(targets, || {
                    self.0
                        .reaction_products_with_params(r, index, p)
                        .map(sets)
                        .map_err(OperationError::from_rust)
                })
            })
        })
    }
    #[wasm_bindgen(js_name=reactionProductsFromInputs)]
    pub fn reaction_products_from_inputs(
        &self,
        reaction: ReactionInput,
        reactants: Molecules,
        params: Option<ReactionRunParamsInput>,
    ) -> Result<ProductSets, JsValue> {
        let reactants = Vec::<ck::Molecule>::input(reactants)?;
        let refs = reactants.iter().collect::<Vec<_>>();
        Reaction::with_mut(&reaction, |r| {
            ReactionRunParams::with_input(params.as_ref(), |p| {
                let targets = r
                    .match_params()
                    .extra_final_check
                    .as_ref()
                    .map(|_| reactants.clone());
                crate::callback::run_with_targets(targets, || {
                    self.0
                        .reaction_products_from_inputs(r, &refs, p)
                        .map(sets)
                        .map_err(OperationError::from_rust)
                })
            })
        })
    }
    #[wasm_bindgen(js_name=applyReaction)]
    pub fn apply_reaction(
        &self,
        reaction: ReactionInput,
        params: Option<ReactionApplyParamsInput>,
    ) -> Result<ReactionApplyResult, JsValue> {
        Reaction::with_mut(&reaction, |r| {
            ReactionApplyParams::with_input(params.as_ref(), |p| {
                let targets = r
                    .match_params()
                    .extra_final_check
                    .as_ref()
                    .map(|_| vec![self.0.clone()]);
                crate::callback::run_with_targets(targets, || {
                    self.0
                        .apply_reaction_with_params(r, p)
                        .map(ReactionApplyResult)
                        .map_err(OperationError::from_rust)
                })
            })
        })
    }
    #[wasm_bindgen(js_name=applyReaction_)]
    pub fn apply_reaction_(
        &mut self,
        reaction: ReactionInput,
        params: Option<ReactionApplyParamsInput>,
    ) -> Result<bool, JsValue> {
        Reaction::with_mut(&reaction, |r| {
            ReactionApplyParams::with_input(params.as_ref(), |p| {
                let targets = r
                    .match_params()
                    .extra_final_check
                    .as_ref()
                    .map(|_| vec![self.0.clone()]);
                crate::callback::run_with_targets(targets, || {
                    self.0
                        .apply_reaction_with_params_(r, p)
                        .map_err(OperationError::from_rust)
                })
            })
        })
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = reactionProductsWithParams)]
    pub fn reaction_products_with_params(
        &self,
        reaction: ReactionInput,
        reactant_template: Ts<input::IndexInput>,
        params: ReactionSingleRunParamsInput,
    ) -> Result<ProductSets, JsValue> {
        self.reaction_products(reaction, reactant_template, Some(params))
    }
    #[wasm_bindgen(js_name = applyReactionWithParams)]
    pub fn apply_reaction_with_params(
        &self,
        reaction: ReactionInput,
        params: ReactionApplyParamsInput,
    ) -> Result<ReactionApplyResult, JsValue> {
        self.apply_reaction(reaction, Some(params))
    }
    #[wasm_bindgen(js_name = applyReactionWithParams_)]
    pub fn apply_reaction_with_params_(
        &mut self,
        reaction: ReactionInput,
        params: ReactionApplyParamsInput,
    ) -> Result<bool, JsValue> {
        self.apply_reaction_(reaction, Some(params))
    }
}
