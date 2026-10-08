//! Frozen parser parameters, copied into the canonical facade without parser logic.
use crate::host_values::{bool_value, type_error};
use cosmolkit_wasm::rust as ck;
use js_sys::{Map, Object, Reflect};
use std::collections::BTreeMap;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct SmilesParseParams {
    pub(crate) inner: ck::SmilesParseParams,
}

#[wasm_bindgen]
impl SmilesParseParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] sanitize: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] allow_cxsmiles: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] strict_cxsmiles: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] parse_name: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_hydrogens: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] skip_cleanup: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] debug_parse: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "Map<string, string> | Record<string, string> | null"
        )]
        replacements: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::SmilesParseParams::default();
        if !sanitize.is_undefined() {
            inner.sanitize = bool_value(&sanitize, "sanitize")?;
        }
        if !allow_cxsmiles.is_undefined() {
            inner.allow_cxsmiles = bool_value(&allow_cxsmiles, "allowCxsmiles")?;
        }
        if !strict_cxsmiles.is_undefined() {
            inner.strict_cxsmiles = bool_value(&strict_cxsmiles, "strictCxsmiles")?;
        }
        if !parse_name.is_undefined() {
            inner.parse_name = bool_value(&parse_name, "parseName")?;
        }
        if !remove_hydrogens.is_undefined() {
            inner.remove_hydrogens = bool_value(&remove_hydrogens, "removeHydrogens")?;
        }
        if !skip_cleanup.is_undefined() {
            inner.skip_cleanup = bool_value(&skip_cleanup, "skipCleanup")?;
        }
        if !debug_parse.is_undefined() {
            inner.debug_parse = bool_value(&debug_parse, "debugParse")?;
        }
        if !replacements.is_null() && !replacements.is_undefined() {
            let mut values = BTreeMap::new();
            if let Some(map) = replacements.dyn_ref::<Map>() {
                for entry in js_sys::try_iter(map)?.ok_or_else(|| type_error("replacements"))? {
                    let pair = js_sys::Array::from(&entry?);
                    let key = pair
                        .get(0)
                        .as_string()
                        .ok_or_else(|| type_error("replacement key"))?;
                    let value = pair
                        .get(1)
                        .as_string()
                        .ok_or_else(|| type_error("replacement value"))?;
                    values.insert(key, value);
                }
            } else {
                if !replacements.is_object() || js_sys::Array::is_array(&replacements) {
                    return Err(type_error("replacements"));
                }
                let object: &Object = replacements.unchecked_ref();
                for key in Object::keys(object).iter() {
                    let value = Reflect::get(object, &key)?
                        .as_string()
                        .ok_or_else(|| type_error("replacement value"))?;
                    values.insert(
                        key.as_string()
                            .ok_or_else(|| type_error("replacement key"))?,
                        value,
                    );
                }
            }
            inner.replacements = values;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter, js_name = sanitize)]
    pub fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    #[wasm_bindgen(getter, js_name = allowCxsmiles)]
    pub fn allow_cxsmiles(&self) -> bool {
        self.inner.allow_cxsmiles
    }
    #[wasm_bindgen(getter, js_name = strictCxsmiles)]
    pub fn strict_cxsmiles(&self) -> bool {
        self.inner.strict_cxsmiles
    }
    #[wasm_bindgen(getter, js_name = parseName)]
    pub fn parse_name(&self) -> bool {
        self.inner.parse_name
    }
    #[wasm_bindgen(getter, js_name = removeHydrogens)]
    pub fn remove_hydrogens(&self) -> bool {
        self.inner.remove_hydrogens
    }
    #[wasm_bindgen(getter, js_name = skipCleanup)]
    pub fn skip_cleanup(&self) -> bool {
        self.inner.skip_cleanup
    }
    #[wasm_bindgen(getter, js_name = debugParse)]
    pub fn debug_parse(&self) -> bool {
        self.inner.debug_parse
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Map<string, string>")]
    pub fn replacements(&self) -> Map {
        let map = Map::new();
        for (key, value) in &self.inner.replacements {
            map.set(&key.into(), &value.into());
        }
        map
    }
}

#[wasm_bindgen]
pub struct SmilesWriteParams {
    pub(crate) inner: ck::SmilesWriteParams,
}
#[wasm_bindgen]
impl SmilesWriteParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] do_isomeric_smiles: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] do_kekule: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] canonical: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] clean_stereo: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] rooted_at_atom: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] all_bonds_explicit: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] all_hydrogens_explicit: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_dative_bonds: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] ignore_atom_map_numbers: JsValue,
    ) -> Result<Self, JsValue> {
        let default = ck::SmilesWriteParams::default();
        Ok(Self {
            inner: ck::SmilesWriteParams {
                do_isomeric_smiles: if do_isomeric_smiles.is_undefined() {
                    default.do_isomeric_smiles
                } else {
                    crate::host_values::bool_value(&do_isomeric_smiles, "doIsomericSmiles")?
                },
                do_kekule: if do_kekule.is_undefined() {
                    default.do_kekule
                } else {
                    crate::host_values::bool_value(&do_kekule, "doKekule")?
                },
                canonical: if canonical.is_undefined() {
                    default.canonical
                } else {
                    crate::host_values::bool_value(&canonical, "canonical")?
                },
                clean_stereo: if clean_stereo.is_undefined() {
                    default.clean_stereo
                } else {
                    crate::host_values::bool_value(&clean_stereo, "cleanStereo")?
                },
                rooted_at_atom: if rooted_at_atom.is_null() || rooted_at_atom.is_undefined() {
                    None
                } else {
                    Some(ck::AtomId::new(crate::host_values::usize_value(
                        &rooted_at_atom,
                        "rootedAtAtom",
                    )?))
                },
                all_bonds_explicit: if all_bonds_explicit.is_undefined() {
                    default.all_bonds_explicit
                } else {
                    crate::host_values::bool_value(&all_bonds_explicit, "allBondsExplicit")?
                },
                all_hydrogens_explicit: if all_hydrogens_explicit.is_undefined() {
                    default.all_hydrogens_explicit
                } else {
                    crate::host_values::bool_value(&all_hydrogens_explicit, "allHydrogensExplicit")?
                },
                include_dative_bonds: if include_dative_bonds.is_undefined() {
                    default.include_dative_bonds
                } else {
                    crate::host_values::bool_value(&include_dative_bonds, "includeDativeBonds")?
                },
                ignore_atom_map_numbers: if ignore_atom_map_numbers.is_undefined() {
                    default.ignore_atom_map_numbers
                } else {
                    crate::host_values::bool_value(
                        &ignore_atom_map_numbers,
                        "ignoreAtomMapNumbers",
                    )?
                },
            },
        })
    }
    #[wasm_bindgen(getter,js_name = doIsomericSmiles)]
    pub fn do_isomeric_smiles(&self) -> bool {
        self.inner.do_isomeric_smiles
    }
    #[wasm_bindgen(getter,js_name = doKekule)]
    pub fn do_kekule(&self) -> bool {
        self.inner.do_kekule
    }
    #[wasm_bindgen(getter,js_name = canonical)]
    pub fn canonical(&self) -> bool {
        self.inner.canonical
    }
    #[wasm_bindgen(getter,js_name = cleanStereo)]
    pub fn clean_stereo(&self) -> bool {
        self.inner.clean_stereo
    }
    #[wasm_bindgen(getter,js_name = rootedAtAtom,unchecked_return_type = "number | null")]
    pub fn rooted_at_atom(&self) -> JsValue {
        self.inner
            .rooted_at_atom
            .map_or(JsValue::NULL, |id| JsValue::from_f64(id.index() as f64))
    }
    #[wasm_bindgen(getter,js_name = allBondsExplicit)]
    pub fn all_bonds_explicit(&self) -> bool {
        self.inner.all_bonds_explicit
    }
    #[wasm_bindgen(getter,js_name = allHydrogensExplicit)]
    pub fn all_hydrogens_explicit(&self) -> bool {
        self.inner.all_hydrogens_explicit
    }
    #[wasm_bindgen(getter,js_name = includeDativeBonds)]
    pub fn include_dative_bonds(&self) -> bool {
        self.inner.include_dative_bonds
    }
    #[wasm_bindgen(getter,js_name = ignoreAtomMapNumbers)]
    pub fn ignore_atom_map_numbers(&self) -> bool {
        self.inner.ignore_atom_map_numbers
    }
}
