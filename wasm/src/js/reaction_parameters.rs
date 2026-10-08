//! Typed, immutable reaction option values.
use crate::host_values::*;
use cosmolkit_wasm::rust as ck;
use js_sys::{Array, Map, Object, Reflect};
use std::collections::BTreeMap;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct CxSmilesFields {
    pub(crate) inner: ck::CxSmilesFields,
}
#[wasm_bindgen]
impl CxSmilesFields {
    pub fn bits(&self) -> u32 {
        self.inner.bits()
    }
    pub fn contains(&self, other: &CxSmilesFields) -> bool {
        self.inner.contains(other.inner)
    }
    pub fn or(&self, other: &CxSmilesFields) -> Self {
        Self {
            inner: self.inner | other.inner,
        }
    }
    #[wasm_bindgen(getter,js_name=NONE)]
    pub fn flag_none() -> Self {
        Self {
            inner: ck::CxSmilesFields::NONE,
        }
    }
    #[wasm_bindgen(getter,js_name=ATOM_LABELS)]
    pub fn flag_atom_labels() -> Self {
        Self {
            inner: ck::CxSmilesFields::ATOM_LABELS,
        }
    }
    #[wasm_bindgen(getter,js_name=MOLFILE_VALUES)]
    pub fn flag_molfile_values() -> Self {
        Self {
            inner: ck::CxSmilesFields::MOLFILE_VALUES,
        }
    }
    #[wasm_bindgen(getter,js_name=COORDS)]
    pub fn flag_coords() -> Self {
        Self {
            inner: ck::CxSmilesFields::COORDS,
        }
    }
    #[wasm_bindgen(getter,js_name=RADICALS)]
    pub fn flag_radicals() -> Self {
        Self {
            inner: ck::CxSmilesFields::RADICALS,
        }
    }
    #[wasm_bindgen(getter,js_name=ATOM_PROPS)]
    pub fn flag_atom_props() -> Self {
        Self {
            inner: ck::CxSmilesFields::ATOM_PROPS,
        }
    }
    #[wasm_bindgen(getter,js_name=LINKNODES)]
    pub fn flag_linknodes() -> Self {
        Self {
            inner: ck::CxSmilesFields::LINKNODES,
        }
    }
    #[wasm_bindgen(getter,js_name=ENHANCED_STEREO)]
    pub fn flag_enhanced_stereo() -> Self {
        Self {
            inner: ck::CxSmilesFields::ENHANCED_STEREO,
        }
    }
    #[wasm_bindgen(getter,js_name=SGROUPS)]
    pub fn flag_sgroups() -> Self {
        Self {
            inner: ck::CxSmilesFields::SGROUPS,
        }
    }
    #[wasm_bindgen(getter,js_name=POLYMER)]
    pub fn flag_polymer() -> Self {
        Self {
            inner: ck::CxSmilesFields::POLYMER,
        }
    }
    #[wasm_bindgen(getter,js_name=BOND_CFG)]
    pub fn flag_bond_cfg() -> Self {
        Self {
            inner: ck::CxSmilesFields::BOND_CFG,
        }
    }
    #[wasm_bindgen(getter,js_name=BOND_ATROPISOMER)]
    pub fn flag_bond_atropisomer() -> Self {
        Self {
            inner: ck::CxSmilesFields::BOND_ATROPISOMER,
        }
    }
    #[wasm_bindgen(getter,js_name=COORDINATE_BONDS)]
    pub fn flag_coordinate_bonds() -> Self {
        Self {
            inner: ck::CxSmilesFields::COORDINATE_BONDS,
        }
    }
    #[wasm_bindgen(getter,js_name=HYDROGEN_BONDS)]
    pub fn flag_hydrogen_bonds() -> Self {
        Self {
            inner: ck::CxSmilesFields::HYDROGEN_BONDS,
        }
    }
    #[wasm_bindgen(getter,js_name=ZERO_BONDS)]
    pub fn flag_zero_bonds() -> Self {
        Self {
            inner: ck::CxSmilesFields::ZERO_BONDS,
        }
    }
    #[wasm_bindgen(getter,js_name=ALL)]
    pub fn flag_all() -> Self {
        Self {
            inner: ck::CxSmilesFields::ALL,
        }
    }
    #[wasm_bindgen(getter,js_name=ALL_BUT_COORDS)]
    pub fn flag_all_but_coords() -> Self {
        Self {
            inner: ck::CxSmilesFields::ALL_BUT_COORDS,
        }
    }
}
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
    pub fn read_replacements(&self) -> Map {
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
#[wasm_bindgen]
impl ReactionRunParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_products: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "ReactionCoordinateSelection[]")]
        coordinate_selections: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::ReactionRunParams::default();
        if !max_products.is_undefined() {
            inner.max_products = u32_value(&max_products, "maxProducts")?;
        }
        if !coordinate_selections.is_undefined() {
            inner.coordinate_selections = selections(&coordinate_selections)?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=maxProducts)]
    pub fn max_products(&self) -> u32 {
        self.inner.max_products()
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
#[wasm_bindgen]
impl ReactionWriteParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] canonical: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] do_isomeric_smiles: JsValue,
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
        if !do_isomeric_smiles.is_undefined() {
            inner.do_isomeric_smiles = bool_value(&do_isomeric_smiles, "doIsomericSmiles")?;
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
    #[wasm_bindgen(getter,js_name=doIsomericSmiles)]
    pub fn do_isomeric_smiles(&self) -> bool {
        self.inner.do_isomeric_smiles()
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
