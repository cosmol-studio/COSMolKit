//! Frozen parser parameters, copied into the canonical facade without parser logic.
use crate::host_values::{bool_value, type_error};
use cosmolkit_wasm::rust as ck;
use js_sys::{Map, Object, Reflect};
use std::collections::BTreeMap;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct SmartsParseParams {
    pub(crate) inner: ck::SmartsParseParams,
}

#[wasm_bindgen]
impl SmartsParseParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] allow_cxsmiles: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] strict_cxsmiles: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] parse_name: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] merge_hs: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] skip_cleanup: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] debug_parse: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "Map<string, string> | Record<string, string> | null"
        )]
        replacements: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::SmartsParseParams::default();

        if !allow_cxsmiles.is_undefined() {
            inner.allow_cxsmiles = bool_value(&allow_cxsmiles, "allowCxsmiles")?;
        }
        if !strict_cxsmiles.is_undefined() {
            inner.strict_cxsmiles = bool_value(&strict_cxsmiles, "strictCxsmiles")?;
        }
        if !parse_name.is_undefined() {
            inner.parse_name = bool_value(&parse_name, "parseName")?;
        }
        if !merge_hs.is_undefined() {
            inner.merge_hs = bool_value(&merge_hs, "mergeHs")?;
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
                    values.insert(key.into(), value.into());
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
                            .ok_or_else(|| type_error("replacement key"))?
                            .into(),
                        value.into(),
                    );
                }
            }
            inner.replacements = values;
        }
        Ok(Self { inner })
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
    #[wasm_bindgen(getter, js_name = mergeHs)]
    pub fn merge_hs(&self) -> bool {
        self.inner.merge_hs
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
    pub fn replacements(&self) -> Result<Map, JsValue> {
        let map = Map::new();
        for (key, value) in &self.inner.replacements {
            map.set(
                &crate::host_values::text(key)?.into(),
                &crate::host_values::text(value)?.into(),
            );
        }
        Ok(map)
    }
}

#[wasm_bindgen]
pub struct QueryGraph {
    pub(crate) inner: ck::QueryGraph,
}
#[wasm_bindgen]
impl QueryGraph {
    #[wasm_bindgen(js_name=fromSmarts)]
    pub fn from_smarts(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: canonical_search.rs::QueryGraph::from_smarts:
        //         ck::search::from_smarts(text)
        // No query parsing or chemical model duplication; one owned result.
        ck::search::from_smarts(text)
            .map(|inner| Self { inner })
            .map_err(|e| parse_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromSmartsWithParams)]
    pub fn from_smarts_with_params(
        text: &str,
        params: &SmartsParseParams,
    ) -> Result<Self, JsValue> {
        ck::search::from_smarts_with_params(text, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| parse_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numAtoms)]
    pub fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    #[wasm_bindgen(js_name=numBonds)]
    pub fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }
    #[wasm_bindgen(unchecked_return_type = "string | null")]
    pub fn name(&self) -> Result<JsValue, JsValue> {
        self.inner
            .prop("_Name")
            .map(|value| {
                value
                    .as_string()
                    .map_err(|e| crate::property_values::property_error(&e).unwrap_or_else(|e| e))
                    .and_then(crate::host_values::text)
            })
            .transpose()
            .map(|value| value.map_or(JsValue::NULL, JsValue::from))
    }
}
#[wasm_bindgen(js_name=parseSmarts)]
pub fn parse_smarts(text: &str) -> Result<QueryGraph, JsValue> {
    ck::parse_smarts(text)
        .map(|inner| QueryGraph { inner })
        .map_err(|e| parse_error(&e).unwrap_or_else(|e| e))
}
#[wasm_bindgen(js_name=parseSmartsWithParams)]
pub fn parse_smarts_with_params(
    text: &str,
    params: &SmartsParseParams,
) -> Result<QueryGraph, JsValue> {
    ck::parse_smarts_with_params(text, &params.inner)
        .map(|inner| QueryGraph { inner })
        .map_err(|e| parse_error(&e).unwrap_or_else(|e| e))
}
#[wasm_bindgen]
pub struct SmartsParseError {
    fields: js_sys::Object,
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl SmartsParseError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "search".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter,js_name=position,unchecked_return_type="number | null")]
    pub fn payload_0(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"position".into())
    }
    #[wasm_bindgen(getter,js_name=character,unchecked_return_type="string | null")]
    pub fn payload_1(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"character".into())
    }
    #[wasm_bindgen(getter,js_name=context,unchecked_return_type="string | null")]
    pub fn payload_2(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"context".into())
    }
    #[wasm_bindgen(getter,js_name=primitiveDetail,unchecked_return_type="string | null")]
    pub fn payload_3(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"primitiveDetail".into())
    }
    #[wasm_bindgen(getter,js_name=ring,unchecked_return_type="number | null")]
    pub fn payload_4(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"ring".into())
    }
    #[wasm_bindgen(getter,js_name=atom,unchecked_return_type="number | null")]
    pub fn payload_5(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"atom".into())
    }
    #[wasm_bindgen(getter,js_name=beginAtom,unchecked_return_type="number | null")]
    pub fn payload_6(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"beginAtom".into())
    }
    #[wasm_bindgen(getter,js_name=endAtom,unchecked_return_type="number | null")]
    pub fn payload_7(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"endAtom".into())
    }
    #[wasm_bindgen(getter,js_name=feature,unchecked_return_type="string | null")]
    pub fn payload_8(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"feature".into())
    }
    #[wasm_bindgen(getter,js_name=carrier,unchecked_return_type="number | null")]
    pub fn payload_9(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"carrier".into())
    }
}
pub(crate) fn parse_error(source: &ck::SmartsParseError) -> Result<JsValue, JsValue> {
    use crate::alignment_values::{set, source_error};
    use ck::SmartsParseError as E;
    use std::error::Error as _;
    let kind = match source {
        E::MissingRecursiveQueryGraph => "MissingRecursiveQueryGraph",
        E::CxLowering(..) => "CxLowering",
        E::QueryGraph(..) => "QueryGraph",
        E::ParserCarrier(..) => "ParserCarrier",
        E::AtomProperty(..) => "AtomProperty",
        E::BondProperty(..) => "BondProperty",
        E::MoleculeProperty(..) => "MoleculeProperty",
        E::UnclosedBracket(_) => "UnclosedBracket",
        E::UnexpectedCharacter { .. } => "UnexpectedCharacter",
        E::UnexpectedEnd(_) => "UnexpectedEnd",
        E::InvalidAtomPrimitive { .. } => "InvalidAtomPrimitive",
        E::UnclosedParenthesis(_) => "UnclosedParenthesis",
        E::UnbalancedRingClosure(_) => "UnbalancedRingClosure",
        E::SelfRingClosure { .. } => "SelfRingClosure",
        E::DuplicateRingBond { .. } => "DuplicateRingBond",
        E::CxSmiles(_) => "CxSmiles",
        E::Parse(_) => "Parse",
        E::UnsupportedFeature(_) => "UnsupportedFeature",
        E::TemplateAttachmentRemap { .. } => "TemplateAttachmentRemap",
    };
    let fields = js_sys::Object::new();
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("SmartsParseError");
    let error: JsValue = error.into();
    set(&error, "domain", "search".into())?;
    set(&error, "kind", kind.into())?;
    set(
        &error,
        "detail",
        SmartsParseError {
            fields: fields.clone(),
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    let number = |key: &str, n: usize| set(&error, key, JsValue::from_f64(n as f64));
    match source {
        E::UnclosedBracket(p) | E::UnclosedParenthesis(p) => number("position", *p)?,
        E::UnexpectedCharacter {
            position,
            character,
            context,
        } => {
            number("position", *position)?;
            set(&error, "character", character.to_string().into())?;
            set(&error, "context", context.clone().into())?;
        }
        E::InvalidAtomPrimitive { position, detail } => {
            number("position", *position)?;
            set(&error, "primitiveDetail", detail.clone().into())?;
        }
        E::UnbalancedRingClosure(r) => set(&error, "ring", (*r).into())?,
        E::SelfRingClosure { ring, atom } => {
            set(&error, "ring", (*ring).into())?;
            number("atom", *atom)?;
        }
        E::DuplicateRingBond {
            ring,
            begin_atom,
            end_atom,
        } => {
            set(&error, "ring", (*ring).into())?;
            number("beginAtom", *begin_atom)?;
            number("endAtom", *end_atom)?;
        }
        E::UnsupportedFeature(feature) => set(&error, "feature", (*feature).into())?,
        E::TemplateAttachmentRemap { carrier, .. } => number("carrier", *carrier)?,
        _ => (),
    }
    let value = js_sys::Reflect::get(&error, &"position".into())?;
    set(
        &fields,
        "position",
        if value.is_undefined() {
            JsValue::NULL
        } else {
            value
        },
    )?;
    let value = js_sys::Reflect::get(&error, &"character".into())?;
    set(
        &fields,
        "character",
        if value.is_undefined() {
            JsValue::NULL
        } else {
            value
        },
    )?;
    let value = js_sys::Reflect::get(&error, &"context".into())?;
    set(
        &fields,
        "context",
        if value.is_undefined() {
            JsValue::NULL
        } else {
            value
        },
    )?;
    let value = js_sys::Reflect::get(&error, &"primitiveDetail".into())?;
    set(
        &fields,
        "primitiveDetail",
        if value.is_undefined() {
            JsValue::NULL
        } else {
            value
        },
    )?;
    let value = js_sys::Reflect::get(&error, &"ring".into())?;
    set(
        &fields,
        "ring",
        if value.is_undefined() {
            JsValue::NULL
        } else {
            value
        },
    )?;
    let value = js_sys::Reflect::get(&error, &"atom".into())?;
    set(
        &fields,
        "atom",
        if value.is_undefined() {
            JsValue::NULL
        } else {
            value
        },
    )?;
    let value = js_sys::Reflect::get(&error, &"beginAtom".into())?;
    set(
        &fields,
        "beginAtom",
        if value.is_undefined() {
            JsValue::NULL
        } else {
            value
        },
    )?;
    let value = js_sys::Reflect::get(&error, &"endAtom".into())?;
    set(
        &fields,
        "endAtom",
        if value.is_undefined() {
            JsValue::NULL
        } else {
            value
        },
    )?;
    let value = js_sys::Reflect::get(&error, &"feature".into())?;
    set(
        &fields,
        "feature",
        if value.is_undefined() {
            JsValue::NULL
        } else {
            value
        },
    )?;
    let value = js_sys::Reflect::get(&error, &"carrier".into())?;
    set(
        &fields,
        "carrier",
        if value.is_undefined() {
            JsValue::NULL
        } else {
            value
        },
    )?;
    Ok(error)
}
