//! Reaction template calls forward the canonical facade.
use crate::{
    fields::Field,
    input,
    query_values::{QueryGraph, QueryGraphArray, QueryGraphInput},
    reaction_errors::*,
    reaction_values::*,
    storage::NativeField,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct Reaction(pub(crate) NativeField<ck::Reaction>);
crate::storage::mutable_owner!(Reaction, ck::Reaction, ReactionInput, "Reaction");
#[wasm_bindgen(js_name = reactionFromSmirks)]
pub fn reaction_from_smirks(
    text: Ts<input::TextInput>,
    params: Option<ReactionParseParamsInput>,
) -> Result<Reaction, JsValue> {
    Reaction::from_smirks(text, params)
}
impl Reaction {
    pub(crate) fn from_native(value: ck::Reaction) -> Self {
        Self(NativeField::owned(value))
    }
}
#[wasm_bindgen]
impl Reaction {
    #[wasm_bindgen(js_name=fromSmirksWithParams)]
    pub fn from_smirks_with_params(
        text: Ts<input::TextInput>,
        params: ReactionParseParamsInput,
    ) -> Result<Self, JsValue> {
        Self::from_smirks(text, Some(params))
    }
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self::from_native(ck::Reaction::new())
    }
    #[wasm_bindgen]
    pub fn copy(&self) -> Result<Self, JsValue> {
        self.0.read(|r| Self::from_native(r.clone()))
    }
    #[wasm_bindgen(js_name=fromSmirks)]
    pub fn from_smirks(
        text: Ts<input::TextInput>,
        params: Option<ReactionParseParamsInput>,
    ) -> Result<Self, JsValue> {
        let text = input::read::<String>(text)?;
        ReactionParseParams::with_input(params.as_ref(), |p| {
            ck::Reaction::from_smirks_with_params(&text, p)
                .map(Self::from_native)
                .map_err(ReactionParseError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=fromTemplates)]
    pub fn from_templates(
        reactants: QueryGraphArray,
        products: QueryGraphArray,
        agents: QueryGraphArray,
    ) -> Result<Self, JsValue> {
        ck::Reaction::from_templates(
            Vec::<ck::QueryGraph>::input(reactants)?,
            Vec::<ck::QueryGraph>::input(products)?,
            Vec::<ck::QueryGraph>::input(agents)?,
        )
        .map(Self::from_native)
        .map_err(ReactionModelError::from_rust)
    }
    #[wasm_bindgen(getter,js_name=isInitialized)]
    pub fn is_initialized(&self) -> Result<bool, JsValue> {
        self.0.read(|r| r.is_initialized())
    }
    #[wasm_bindgen(getter,js_name=implicitProperties)]
    pub fn implicit_properties(&self) -> Result<bool, JsValue> {
        self.0.read(|r| r.implicit_properties())
    }
    #[wasm_bindgen(js_name=withImplicitProperties)]
    pub fn with_implicit_properties(&self, enabled: Ts<input::BoolInput>) -> Result<Self, JsValue> {
        let enabled = input::read::<bool>(enabled)?;
        self.0
            .read(|r| Self::from_native(r.with_implicit_properties(enabled)))
    }
    #[wasm_bindgen]
    pub fn validate(
        &self,
        params: Option<ReactionValidationParamsInput>,
    ) -> Result<ReactionValidationReport, JsValue> {
        ReactionValidationParams::with_input(params.as_ref(), |p| {
            self.0.read(|r| {
                r.validate_with_params(p)
                    .map(ReactionValidationReport)
                    .map_err(ReactionValidationError::from_rust)
            })?
        })
    }
    #[wasm_bindgen(js_name=withInitialized)]
    pub fn with_initialized(
        &self,
        params: Option<ReactionValidationParamsInput>,
    ) -> Result<Self, JsValue> {
        ReactionValidationParams::with_input(params.as_ref(), |p| {
            self.0.read(|r| {
                r.with_initialized_with_params(p)
                    .map(Self::from_native)
                    .map_err(ReactionInitializationError::from_rust)
            })?
        })
    }
    #[wasm_bindgen(js_name=withoutAgents)]
    pub fn without_agents(&self) -> Result<ReactionTemplateRemoval, JsValue> {
        self.0.read(|r| ReactionTemplateRemoval(r.without_agents()))
    }
    #[wasm_bindgen(getter,js_name=reactantTemplates)]
    pub fn reactant_templates(&self) -> Result<QueryGraphArray, JsValue> {
        self.0.read(|r| {
            Vec::<ck::QueryGraph>::get(NativeField::owned(r.reactant_templates().to_vec()))
        })?
    }
    #[wasm_bindgen(getter,js_name=numReactantTemplates)]
    pub fn num_reactant_templates(&self) -> Result<usize, JsValue> {
        self.0.read(|r| r.num_reactant_templates())
    }
    #[wasm_bindgen(js_name=reactantTemplate)]
    pub fn reactant_template(&self, index: Ts<input::IndexInput>) -> Result<QueryGraph, JsValue> {
        let index = input::read::<usize>(index)?;
        self.0.read(|r| {
            r.reactant_template(index)
                .cloned()
                .map(QueryGraph::from_native)
                .map_err(ReactionModelError::from_rust)
        })?
    }
    #[wasm_bindgen(js_name=withReactantTemplate)]
    pub fn with_reactant_template(&self, template: QueryGraphInput) -> Result<Self, JsValue> {
        let template = ck::QueryGraph::input(template)?;
        self.0.read(|r| {
            r.with_reactant_template(template)
                .map(Self::from_native)
                .map_err(ReactionModelError::from_rust)
        })?
    }
    #[wasm_bindgen(getter,js_name=productTemplates)]
    pub fn product_templates(&self) -> Result<QueryGraphArray, JsValue> {
        self.0.read(|r| {
            Vec::<ck::QueryGraph>::get(NativeField::owned(r.product_templates().to_vec()))
        })?
    }
    #[wasm_bindgen(getter,js_name=numProductTemplates)]
    pub fn num_product_templates(&self) -> Result<usize, JsValue> {
        self.0.read(|r| r.num_product_templates())
    }
    #[wasm_bindgen(js_name=productTemplate)]
    pub fn product_template(&self, index: Ts<input::IndexInput>) -> Result<QueryGraph, JsValue> {
        let index = input::read::<usize>(index)?;
        self.0.read(|r| {
            r.product_template(index)
                .cloned()
                .map(QueryGraph::from_native)
                .map_err(ReactionModelError::from_rust)
        })?
    }
    #[wasm_bindgen(js_name=withProductTemplate)]
    pub fn with_product_template(&self, template: QueryGraphInput) -> Result<Self, JsValue> {
        let template = ck::QueryGraph::input(template)?;
        self.0.read(|r| {
            r.with_product_template(template)
                .map(Self::from_native)
                .map_err(ReactionModelError::from_rust)
        })?
    }
    #[wasm_bindgen(getter,js_name=agentTemplates)]
    pub fn agent_templates(&self) -> Result<QueryGraphArray, JsValue> {
        self.0.read(|r| {
            Vec::<ck::QueryGraph>::get(NativeField::owned(r.agent_templates().to_vec()))
        })?
    }
    #[wasm_bindgen(getter,js_name=numAgentTemplates)]
    pub fn num_agent_templates(&self) -> Result<usize, JsValue> {
        self.0.read(|r| r.num_agent_templates())
    }
    #[wasm_bindgen(js_name=agentTemplate)]
    pub fn agent_template(&self, index: Ts<input::IndexInput>) -> Result<QueryGraph, JsValue> {
        let index = input::read::<usize>(index)?;
        self.0.read(|r| {
            r.agent_template(index)
                .cloned()
                .map(QueryGraph::from_native)
                .map_err(ReactionModelError::from_rust)
        })?
    }
    #[wasm_bindgen(js_name=withAgentTemplate)]
    pub fn with_agent_template(&self, template: QueryGraphInput) -> Result<Self, JsValue> {
        let template = ck::QueryGraph::input(template)?;
        self.0.read(|r| {
            r.with_agent_template(template)
                .map(Self::from_native)
                .map_err(ReactionModelError::from_rust)
        })?
    }
    #[wasm_bindgen(js_name=toSmirks)]
    pub fn to_smirks(&self, params: Option<ReactionWriteParamsInput>) -> Result<String, JsValue> {
        ReactionWriteParams::with_input(params.as_ref(), |p| {
            let value = self.0.read(|r| {
                r.to_smirks_with_params(p)
                    .map_err(ReactionWriteError::from_rust)
            })??;
            crate::property_errors::text(&value)
        })
    }
    #[wasm_bindgen(js_name=toCxSmirks)]
    pub fn to_cx_smirks(
        &self,
        params: Option<ReactionWriteParamsInput>,
    ) -> Result<String, JsValue> {
        ReactionWriteParams::with_input(params.as_ref(), |p| {
            let value = self.0.read(|r| {
                r.to_cx_smirks_with_params(p)
                    .map_err(ReactionWriteError::from_rust)
            })??;
            crate::property_errors::text(&value)
        })
    }
    #[wasm_bindgen(js_name=withoutUnmappedReactants)]
    pub fn without_unmapped_reactants(
        &self,
        params: Option<ReactionTemplateRemovalParamsInput>,
    ) -> Result<ReactionTemplateRemoval, JsValue> {
        ReactionTemplateRemovalParams::with_input(params.as_ref(), |p| {
            self.0
                .read(|r| ReactionTemplateRemoval(r.without_unmapped_reactants_with_params(p)))
        })
    }
    #[wasm_bindgen(js_name=withoutUnmappedProducts)]
    pub fn without_unmapped_products(
        &self,
        params: Option<ReactionTemplateRemovalParamsInput>,
    ) -> Result<ReactionTemplateRemoval, JsValue> {
        ReactionTemplateRemovalParams::with_input(params.as_ref(), |p| {
            self.0
                .read(|r| ReactionTemplateRemoval(r.without_unmapped_products_with_params(p)))
        })
    }
}
#[wasm_bindgen]
pub struct ReactionTemplateRemoval(ck::ReactionTemplateRemoval);
#[wasm_bindgen]
impl ReactionTemplateRemoval {
    #[wasm_bindgen(getter)]
    pub fn reaction(&self) -> Reaction {
        Reaction::from_native(self.0.reaction().clone())
    }
    #[wasm_bindgen(getter,js_name=removedTemplates)]
    pub fn removed_templates(&self) -> Result<QueryGraphArray, JsValue> {
        Vec::<ck::QueryGraph>::get(NativeField::owned(self.0.removed_templates().to_vec()))
    }
}
#[wasm_bindgen(js_name=parseSmirks)]
pub fn parse_smirks(
    text: Ts<input::TextInput>,
    params: Option<ReactionParseParamsInput>,
) -> Result<Reaction, JsValue> {
    Reaction::from_smirks(text, params)
}

#[wasm_bindgen(js_name=parseSmirksWithParams)]
pub fn parse_smirks_with_params(
    text: Ts<input::TextInput>,
    params: ReactionParseParamsInput,
) -> Result<Reaction, JsValue> {
    Reaction::from_smirks(text, Some(params))
}
#[wasm_bindgen]
impl Reaction {
    #[wasm_bindgen(js_name=validateWithParams)]
    pub fn validate_with_params(
        &self,
        params: ReactionValidationParamsInput,
    ) -> Result<ReactionValidationReport, JsValue> {
        self.validate(Some(params))
    }
    #[wasm_bindgen(js_name=withInitializedWithParams)]
    pub fn with_initialized_with_params(
        &self,
        params: ReactionValidationParamsInput,
    ) -> Result<Self, JsValue> {
        self.with_initialized(Some(params))
    }
    #[wasm_bindgen(js_name=toSmirksWithParams)]
    pub fn to_smirks_with_params(
        &self,
        params: ReactionWriteParamsInput,
    ) -> Result<String, JsValue> {
        self.to_smirks(Some(params))
    }
    #[wasm_bindgen(js_name=toCxSmirksWithParams)]
    pub fn to_cx_smirks_with_params(
        &self,
        params: ReactionWriteParamsInput,
    ) -> Result<String, JsValue> {
        self.to_cx_smirks(Some(params))
    }
    #[wasm_bindgen(js_name=withoutUnmappedReactantsWithParams)]
    pub fn without_unmapped_reactants_with_params(
        &self,
        params: ReactionTemplateRemovalParamsInput,
    ) -> Result<ReactionTemplateRemoval, JsValue> {
        self.without_unmapped_reactants(Some(params))
    }
    #[wasm_bindgen(js_name=withoutUnmappedProductsWithParams)]
    pub fn without_unmapped_products_with_params(
        &self,
        params: ReactionTemplateRemovalParamsInput,
    ) -> Result<ReactionTemplateRemoval, JsValue> {
        self.without_unmapped_products(Some(params))
    }
}
