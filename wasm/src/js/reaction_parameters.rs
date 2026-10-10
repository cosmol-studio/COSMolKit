//! Typed, writable reaction option values.
use crate::host_values::*;
use crate::notation_writers::CxSmilesFields;
use cosmolkit_wasm::rust as ck;
use js_sys::{Array, Map, Object, Reflect};
use std::collections::BTreeMap;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct ReactionCoordinateSelection {
    pub(crate) inner: ck::ReactionCoordinateSelection,
}
#[wasm_bindgen]
impl ReactionCoordinateSelection {
    pub fn auto() -> Self {
        Self {
            inner: ck::ReactionCoordinateSelection::auto(),
        }
    }
    #[wasm_bindgen(js_name=twoD)]
    pub fn two_d(
        #[wasm_bindgen(unchecked_param_type = "number")] id: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::ReactionCoordinateSelection::two_d(usize_value(&id, "id")?),
        })
    }
    #[wasm_bindgen(js_name=threeD)]
    pub fn three_d(
        #[wasm_bindgen(unchecked_param_type = "number")] id: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::ReactionCoordinateSelection::three_d(usize_value(&id, "id")?),
        })
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn id(&self) -> JsValue {
        self.inner
            .id()
            .map_or(JsValue::NULL, |v| JsValue::from(v as f64))
    }
    #[wasm_bindgen(getter,js_name=isAuto)]
    pub fn is_auto(&self) -> bool {
        self.inner.is_auto()
    }
    #[wasm_bindgen(getter,js_name=is2d)]
    pub fn is_2d(&self) -> bool {
        self.inner.is_2d()
    }
    #[wasm_bindgen(getter,js_name=is3d)]
    pub fn is_3d(&self) -> bool {
        self.inner.is_3d()
    }
}
#[wasm_bindgen(
    inline_js = "export function visitReactionCoordinateSelection(v,f){try{f(v);}catch(cause){throw new TypeError('invalid ReactionCoordinateSelection',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitReactionCoordinateSelection)]
    fn visitReactionCoordinateSelection(
        value: &JsValue,
        f: &mut dyn FnMut(&ReactionCoordinateSelection),
    ) -> Result<(), JsValue>;
}
#[wasm_bindgen(
    inline_js = "export function visitReactionCxSmilesFields(v,f){try{f(v);}catch(cause){throw new TypeError('invalid CxSmilesFields',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitReactionCxSmilesFields)]
    fn visitReactionCxSmilesFields(
        value: &JsValue,
        f: &mut dyn FnMut(&CxSmilesFields),
    ) -> Result<(), JsValue>;
}
fn selection(value: &JsValue) -> Result<ck::ReactionCoordinateSelection, JsValue> {
    let mut out = None;
    visitReactionCoordinateSelection(value, &mut |v: &ReactionCoordinateSelection| {
        out = Some(v.inner)
    })?;
    out.ok_or_else(|| type_error("coordinate selection"))
}
fn selections(value: &JsValue) -> Result<Vec<ck::ReactionCoordinateSelection>, JsValue> {
    sequence(value, "coordinate selections")?
        .iter()
        .map(|v| selection(&v))
        .collect()
}
fn read_replacements(value: &JsValue) -> Result<BTreeMap<String, String>, JsValue> {
    let mut out = BTreeMap::new();
    if let Some(map) = value.dyn_ref::<Map>() {
        for pair in js_sys::try_iter(map)?.ok_or_else(|| type_error("replacements"))? {
            let pair = Array::from(&pair?);
            out.insert(
                pair.get(0)
                    .as_string()
                    .ok_or_else(|| type_error("replacement key"))?,
                pair.get(1)
                    .as_string()
                    .ok_or_else(|| type_error("replacement value"))?,
            );
        }
    } else {
        if !value.is_object() || Array::is_array(value) {
            return Err(type_error("replacements"));
        }
        for key in Object::keys(&Object::from(value.clone())).iter() {
            let v = Reflect::get(value, &key)?;
            out.insert(
                key.as_string()
                    .ok_or_else(|| type_error("replacement key"))?,
                v.as_string()
                    .ok_or_else(|| type_error("replacement value"))?,
            );
        }
    }
    Ok(out)
}
#[wasm_bindgen]
pub struct ReactionParseParams {
    pub(crate) inner: ck::ReactionParseParams,
}
#[cosmolkit_wasm::javascript_options("ReactionParseParamsOptions")]
#[wasm_bindgen]
impl ReactionParseParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_smiles: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] sanitize: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "Map<string,string> | Record<string,string>"
        )]
        replacements: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] allow_cxsmiles: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] strict_cxsmiles: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::ReactionParseParams::default();
        if !use_smiles.is_undefined() {
            inner.use_smiles = bool_value(&use_smiles, "useSmiles")?;
        }
        if !sanitize.is_undefined() {
            inner.sanitize = bool_value(&sanitize, "sanitize")?;
        }
        if !replacements.is_undefined() {
            inner.replacements = read_replacements(&replacements)?;
        }
        if !allow_cxsmiles.is_undefined() {
            inner.allow_cxsmiles = bool_value(&allow_cxsmiles, "allowCxsmiles")?;
        }
        if !strict_cxsmiles.is_undefined() {
            inner.strict_cxsmiles = bool_value(&strict_cxsmiles, "strictCxsmiles")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=useSmiles)]
    pub fn use_smiles(&self) -> bool {
        self.inner.use_smiles()
    }
    #[wasm_bindgen(getter,js_name=sanitize)]
    pub fn sanitize(&self) -> bool {
        self.inner.sanitize()
    }
    #[wasm_bindgen(getter,js_name=replacements,unchecked_return_type="Map<string,string>")]
    pub fn replacements(&self) -> Map {
        {
            let map = Map::new();
            for (k, v) in self.inner.replacements() {
                map.set(&k.clone().into(), &v.clone().into());
            }
            map
        }
    }
    #[wasm_bindgen(getter,js_name=allowCxsmiles)]
    pub fn allow_cxsmiles(&self) -> bool {
        self.inner.allow_cxsmiles()
    }
    #[wasm_bindgen(getter,js_name=strictCxsmiles)]
    pub fn strict_cxsmiles(&self) -> bool {
        self.inner.strict_cxsmiles()
    }
}
#[wasm_bindgen]
pub struct ReactionValidationParams {
    pub(crate) inner: ck::ReactionValidationParams,
}
#[cosmolkit_wasm::javascript_options("ReactionValidationParamsOptions")]
#[wasm_bindgen]
impl ReactionValidationParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] silent: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::ReactionValidationParams::default();
        if !silent.is_undefined() {
            inner.silent = bool_value(&silent, "silent")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=silent)]
    pub fn silent(&self) -> bool {
        self.inner.silent()
    }
}
#[wasm_bindgen]
pub struct ReactionSingleRunParams {
    pub(crate) inner: ck::ReactionSingleRunParams,
}
#[cosmolkit_wasm::javascript_options("ReactionSingleRunParamsOptions")]
#[wasm_bindgen]
impl ReactionSingleRunParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "ReactionCoordinateSelection")]
        coordinate_selection: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::ReactionSingleRunParams::default();
        if !coordinate_selection.is_undefined() {
            inner.coordinate_selection = selection(&coordinate_selection)?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=coordinateSelection)]
    pub fn coordinate_selection(&self) -> ReactionCoordinateSelection {
        ReactionCoordinateSelection {
            inner: self.inner.coordinate_selection(),
        }
    }
}
#[wasm_bindgen]
pub struct ReactionRunParams {
    pub(crate) inner: ck::ReactionRunParams,
}
#[cosmolkit_wasm::javascript_options("ReactionRunParamsOptions")]
#[wasm_bindgen]
impl ReactionRunParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_products: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "ReactionCoordinateSelection[]")]
        coordinate_selections: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] copy_atom_properties: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::ReactionRunParams::default();
        if !max_products.is_undefined() {
            inner.max_products = u32_value(&max_products, "maxProducts")?;
        }
        if !coordinate_selections.is_undefined() {
            inner.coordinate_selections = selections(&coordinate_selections)?;
        }
        if !copy_atom_properties.is_undefined() {
            inner.copy_atom_properties = bool_value(&copy_atom_properties, "copyAtomProperties")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=maxProducts)]
    pub fn max_products(&self) -> u32 {
        self.inner.max_products()
    }
    #[wasm_bindgen(getter,js_name=copyAtomProperties)]
    pub fn copy_atom_properties(&self) -> bool {
        self.inner.copy_atom_properties()
    }
    #[wasm_bindgen(setter,js_name=copyAtomProperties)]
    pub fn set_copy_atom_properties(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.copy_atom_properties = bool_value(&value, "copyAtomProperties")?;
        Ok(())
    }
    #[wasm_bindgen(setter,js_name=maxProducts)]
    pub fn set_max_products(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.max_products = u32_value(&value, "maxProducts")?;
        Ok(())
    }
    #[wasm_bindgen(setter,js_name=coordinateSelections)]
    pub fn set_coordinate_selections(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "ReactionCoordinateSelection[]")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.coordinate_selections = selections(&value)?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=coordinateSelections,unchecked_return_type="ReactionCoordinateSelection[]")]
    pub fn coordinate_selections(&self) -> Array {
        self.inner
            .coordinate_selections()
            .iter()
            .map(|v| JsValue::from(ReactionCoordinateSelection { inner: *v }))
            .collect()
    }
}
#[wasm_bindgen]
pub struct ReactionApplyParams {
    pub(crate) inner: ck::ReactionApplyParams,
}
#[cosmolkit_wasm::javascript_options("ReactionApplyParamsOptions")]
#[wasm_bindgen]
impl ReactionApplyParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_unmatched_atoms: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::ReactionApplyParams::default();
        if !remove_unmatched_atoms.is_undefined() {
            inner.remove_unmatched_atoms =
                bool_value(&remove_unmatched_atoms, "removeUnmatchedAtoms")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=removeUnmatchedAtoms)]
    pub fn remove_unmatched_atoms(&self) -> bool {
        self.inner.remove_unmatched_atoms()
    }
}
#[wasm_bindgen]
pub struct ReactionTemplateRemovalParams {
    pub(crate) inner: ck::ReactionTemplateRemovalParams,
}
#[cosmolkit_wasm::javascript_options("ReactionTemplateRemovalParamsOptions")]
#[wasm_bindgen]
impl ReactionTemplateRemovalParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] threshold_unmapped_atoms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] move_to_agent_templates: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::ReactionTemplateRemovalParams::default();
        if !threshold_unmapped_atoms.is_undefined() {
            inner.threshold_unmapped_atoms = threshold_unmapped_atoms
                .as_f64()
                .ok_or_else(|| type_error("thresholdUnmappedAtoms"))?;
        }
        if !move_to_agent_templates.is_undefined() {
            inner.move_to_agent_templates =
                bool_value(&move_to_agent_templates, "moveToAgentTemplates")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=thresholdUnmappedAtoms)]
    pub fn threshold_unmapped_atoms(&self) -> f64 {
        self.inner.threshold_unmapped_atoms()
    }
    #[wasm_bindgen(getter,js_name=moveToAgentTemplates)]
    pub fn move_to_agent_templates(&self) -> bool {
        self.inner.move_to_agent_templates()
    }
}
#[wasm_bindgen]
pub struct ReactionWriteParams {
    pub(crate) inner: ck::ReactionWriteParams,
}
#[cosmolkit_wasm::javascript_options("ReactionWriteParamsOptions")]
#[wasm_bindgen]
impl ReactionWriteParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] canonical: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] isomeric_smiles: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] rooted_at_atom: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_dative_bonds: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_cx: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "CxSmilesFields")] cx_fields: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "ReactionCoordinateSelection[]")]
        coordinate_selections: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::ReactionWriteParams::default();
        if !canonical.is_undefined() {
            inner.canonical = bool_value(&canonical, "canonical")?;
        }
        if !isomeric_smiles.is_undefined() {
            inner.isomeric_smiles = bool_value(&isomeric_smiles, "isomericSmiles")?;
        }
        if !rooted_at_atom.is_undefined() {
            inner.rooted_at_atom = if rooted_at_atom.is_null() {
                None
            } else {
                Some(usize_value(&rooted_at_atom, "rootedAtAtom")?)
            };
        }
        if !include_dative_bonds.is_undefined() {
            inner.include_dative_bonds = bool_value(&include_dative_bonds, "includeDativeBonds")?;
        }
        if !include_cx.is_undefined() {
            inner.include_cx = bool_value(&include_cx, "includeCx")?;
        }
        if !cx_fields.is_undefined() {
            inner.cx_fields = {
                let mut out = None;
                visitReactionCxSmilesFields(&cx_fields, &mut |v: &CxSmilesFields| {
                    out = Some(v.inner)
                })?;
                out.ok_or_else(|| type_error("cxFields"))?
            };
        }
        if !coordinate_selections.is_undefined() {
            inner.coordinate_selections = selections(&coordinate_selections)?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=canonical)]
    pub fn canonical(&self) -> bool {
        self.inner.canonical()
    }
    #[wasm_bindgen(getter,js_name=isomericSmiles)]
    pub fn isomeric_smiles(&self) -> bool {
        self.inner.isomeric_smiles()
    }
    #[wasm_bindgen(getter,js_name=rootedAtAtom,unchecked_return_type="number | null")]
    pub fn rooted_at_atom(&self) -> JsValue {
        self.inner
            .rooted_at_atom()
            .map_or(JsValue::NULL, |v| JsValue::from(v as f64))
    }
    #[wasm_bindgen(getter,js_name=includeDativeBonds)]
    pub fn include_dative_bonds(&self) -> bool {
        self.inner.include_dative_bonds()
    }
    #[wasm_bindgen(getter,js_name=includeCx)]
    pub fn include_cx(&self) -> bool {
        self.inner.include_cx()
    }
    #[wasm_bindgen(getter,js_name=cxFields)]
    pub fn cx_fields(&self) -> CxSmilesFields {
        CxSmilesFields {
            inner: self.inner.cx_fields(),
        }
    }
    #[wasm_bindgen(getter,js_name=coordinateSelections,unchecked_return_type="ReactionCoordinateSelection[]")]
    pub fn coordinate_selections(&self) -> Array {
        self.inner
            .coordinate_selections()
            .iter()
            .map(|v| JsValue::from(ReactionCoordinateSelection { inner: *v }))
            .collect()
    }
}
