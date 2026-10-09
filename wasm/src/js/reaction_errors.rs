//! Structured reaction errors preserve discriminants, payloads and native causes.
use crate::alignment_values::source_error;
use crate::reaction::{ReactionRole, ReactionValidationReport};
use cosmolkit_wasm::rust as ck;
use js_sys::{Error, Object, Reflect};
use std::error::Error as _;
use wasm_bindgen::prelude::*;
fn set_field(target: &JsValue, key: &str, value: JsValue) -> Result<(), JsValue> {
    Reflect::set(target, &key.into(), &value)?;
    Ok(())
}
#[wasm_bindgen]
pub struct ReactionModelError {
    fields: Object,
    message: String,
    kind: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl ReactionModelError {
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter,js_name=role,unchecked_return_type="ReactionRole | null")]
    pub fn role(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"role".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=template,unchecked_return_type="number | null")]
    pub fn template(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"template".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=index,unchecked_return_type="number | null")]
    pub fn index(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"index".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=count,unchecked_return_type="number | null")]
    pub fn count(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"count".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
}
pub(crate) fn model_error(source: &ck::ReactionModelError) -> Result<JsValue, JsValue> {
    use ck::ReactionModelError as E;
    let kind = match source {
        E::Template { .. } => "Template",
        E::TemplateIndex { .. } => "TemplateIndex",
    };
    let error = Error::new(&source.to_string());
    error.set_name("ReactionModelError");
    let error: JsValue = error.into();
    set_field(&error, "domain", "reaction".into())?;
    set_field(&error, "kind", kind.into())?;
    let fields = Object::new();
    match source {
        E::Template { role, template, .. } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
        }
        E::TemplateIndex {
            role, index, count, ..
        } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "index", JsValue::from(*index as f64))?;
            set_field(&fields, "count", JsValue::from(*count as f64))?;
        }
        _ => (),
    }
    let cause = match source {
        _ => source
            .source()
            .map(source_error)
            .transpose()?
            .unwrap_or(JsValue::NULL),
    };
    if !cause.is_null() {
        set_field(&error, "cause", cause.clone())?;
    }
    for key in Object::keys(&fields).iter() {
        let value = Reflect::get(&fields, &key)?;
        Reflect::set(&error, &key, &value)?;
    }
    set_field(
        &error,
        "detail",
        ReactionModelError {
            fields,
            message: source.to_string(),
            kind: kind.into(),
            cause,
        }
        .into(),
    )?;
    Ok(error)
}
#[wasm_bindgen]
pub struct ReactionParseError {
    fields: Object,
    message: String,
    kind: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl ReactionParseError {
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter,js_name=role,unchecked_return_type="ReactionRole | null")]
    pub fn role(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"role".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=template,unchecked_return_type="number | null")]
    pub fn template(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"template".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=atom,unchecked_return_type="number | null")]
    pub fn atom(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"atom".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=atomCount,unchecked_return_type="number | null")]
    pub fn atom_count(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"atomCount".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=reactant,unchecked_return_type="number | null")]
    pub fn reactant(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"reactant".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=product,unchecked_return_type="number | null")]
    pub fn product(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"product".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=flag,unchecked_return_type="number | null")]
    pub fn flag(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"flag".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=count,unchecked_return_type="number | null")]
    pub fn count(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"count".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=start,unchecked_return_type="number | null")]
    pub fn start(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"start".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=end,unchecked_return_type="number | null")]
    pub fn end(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"end".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=text,unchecked_return_type="string | null")]
    pub fn text(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"text".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=startAtom,unchecked_return_type="number | null")]
    pub fn start_atom(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"startAtom".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=startBond,unchecked_return_type="number | null")]
    pub fn start_bond(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"startBond".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
}
pub(crate) fn parse_error(source: &ck::ReactionParseError) -> Result<JsValue, JsValue> {
    use ck::ReactionParseError as E;
    let kind = match source {
        E::ParserCleanup { .. } => "ParserCleanup",
        E::TemplateAtomBounds { .. } => "TemplateAtomBounds",
        E::StereoDegreeArithmetic { .. } => "StereoDegreeArithmetic",
        E::StereoInversionFlag { .. } => "StereoInversionFlag",
        E::TemplateProperty { .. } => "TemplateProperty",
        E::Separators { .. } => "Separators",
        E::MultiStep { .. } => "MultiStep",
        E::ComponentBounds { .. } => "ComponentBounds",
        E::OffsetOverflow => "OffsetOverflow",
        E::Smarts { .. } => "Smarts",
        E::Smiles { .. } => "Smiles",
        E::ComponentEncoding { .. } => "ComponentEncoding",
        E::Model { .. } => "Model",
        E::AgentFragments { .. } => "AgentFragments",
        E::CxParse { .. } => "CxParse",
        E::CxLowering { .. } => "CxLowering",
        E::StereoOrder { .. } => "StereoOrder",
    };
    let error = Error::new(&source.to_string());
    error.set_name("ReactionParseError");
    let error: JsValue = error.into();
    set_field(&error, "domain", "reaction".into())?;
    set_field(&error, "kind", kind.into())?;
    let fields = Object::new();
    match source {
        E::ParserCleanup { role, template, .. } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
        }
        E::TemplateAtomBounds {
            role,
            template,
            atom,
            atom_count,
            ..
        } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
            set_field(&fields, "atomCount", JsValue::from(*atom_count as f64))?;
        }
        E::StereoDegreeArithmetic {
            reactant, product, ..
        } => {
            set_field(&fields, "reactant", JsValue::from(*reactant as f64))?;
            set_field(&fields, "product", JsValue::from(*product as f64))?;
        }
        E::StereoInversionFlag {
            product,
            atom,
            flag,
            ..
        } => {
            set_field(&fields, "product", JsValue::from(*product as f64))?;
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
            set_field(
                &fields,
                "flag",
                flag.map_or(JsValue::NULL, |v| JsValue::from(v as f64)),
            )?;
        }
        E::Separators { count, .. } => {
            set_field(&fields, "count", JsValue::from(*count as f64))?;
        }
        E::MultiStep { count, .. } => {
            set_field(&fields, "count", JsValue::from(*count as f64))?;
        }
        E::ComponentBounds { start, end, .. } => {
            set_field(&fields, "start", JsValue::from(*start as f64))?;
            set_field(&fields, "end", JsValue::from(*end as f64))?;
        }
        E::Smarts {
            role,
            template,
            text,
            ..
        } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
            set_field(&fields, "text", crate::host_values::text(text)?.into())?;
        }
        E::Smiles {
            role,
            template,
            text,
            ..
        } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
            set_field(&fields, "text", crate::host_values::text(text)?.into())?;
        }
        E::ComponentEncoding {
            role,
            template,
            text,
            ..
        } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
            set_field(&fields, "text", crate::host_values::text(text)?.into())?;
        }
        E::Model { role, template, .. } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
        }
        E::CxParse {
            role,
            template,
            start_atom,
            start_bond,
            ..
        } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
            set_field(&fields, "startAtom", JsValue::from(*start_atom as f64))?;
            set_field(&fields, "startBond", JsValue::from(*start_bond as f64))?;
        }
        E::CxLowering {
            role,
            template,
            start_atom,
            start_bond,
            ..
        } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
            set_field(&fields, "startAtom", JsValue::from(*start_atom as f64))?;
            set_field(&fields, "startBond", JsValue::from(*start_bond as f64))?;
        }
        E::StereoOrder { product, atom, .. } => {
            set_field(&fields, "product", JsValue::from(*product as f64))?;
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
        }
        _ => (),
    }
    let cause = match source {
        _ => source
            .source()
            .map(source_error)
            .transpose()?
            .unwrap_or(JsValue::NULL),
    };
    if !cause.is_null() {
        set_field(&error, "cause", cause.clone())?;
    }
    for key in Object::keys(&fields).iter() {
        let value = Reflect::get(&fields, &key)?;
        Reflect::set(&error, &key, &value)?;
    }
    set_field(
        &error,
        "detail",
        ReactionParseError {
            fields,
            message: source.to_string(),
            kind: kind.into(),
            cause,
        }
        .into(),
    )?;
    Ok(error)
}
#[wasm_bindgen]
pub struct ReactionRunError {
    fields: Object,
    message: String,
    kind: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl ReactionRunError {
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter,js_name=set,unchecked_return_type="number | null")]
    pub fn set_field(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"set".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=template,unchecked_return_type="number | null")]
    pub fn template(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"template".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=expected,unchecked_return_type="number | null")]
    pub fn expected(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"expected".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=actual,unchecked_return_type="number | null")]
    pub fn actual(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"actual".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=index,unchecked_return_type="number | null")]
    pub fn index(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"index".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=count,unchecked_return_type="number | null")]
    pub fn count(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"count".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=reactant,unchecked_return_type="number | null")]
    pub fn reactant(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"reactant".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=atom,unchecked_return_type="number | null")]
    pub fn atom(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"atom".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=atomCount,unchecked_return_type="number | null")]
    pub fn atom_count(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"atomCount".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=level,unchecked_return_type="number | null")]
    pub fn level(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"level".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
}
pub(crate) fn run_error(source: &ck::ReactionRunError) -> Result<JsValue, JsValue> {
    use ck::ReactionRunError as E;
    let kind = match source {
        E::Product { .. } => "Product",
        E::CoordinateSelectionArity { .. } => "CoordinateSelectionArity",
        E::Initialization(..) => "Initialization",
        E::NeedsInitialization => "NeedsInitialization",
        E::ReactantArity { .. } => "ReactantArity",
        E::ReactantTemplateIndex { .. } => "ReactantTemplateIndex",
        E::Matching { .. } => "Matching",
        E::MatchAtomIndex { .. } => "MatchAtomIndex",
        E::CombinationLevel { .. } => "CombinationLevel",
        E::CombinationSize { .. } => "CombinationSize",
        E::EmptyCombinationLevels => "EmptyCombinationLevels",
    };
    let error = Error::new(&source.to_string());
    error.set_name("ReactionRunError");
    let error: JsValue = error.into();
    set_field(&error, "domain", "reaction".into())?;
    set_field(&error, "kind", kind.into())?;
    let fields = Object::new();
    match source {
        E::Product { set, template, .. } => {
            set_field(&fields, "set", JsValue::from(*set as f64))?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
        }
        E::CoordinateSelectionArity {
            expected, actual, ..
        } => {
            set_field(&fields, "expected", JsValue::from(*expected as f64))?;
            set_field(&fields, "actual", JsValue::from(*actual as f64))?;
        }
        E::ReactantArity {
            expected, actual, ..
        } => {
            set_field(&fields, "expected", JsValue::from(*expected as f64))?;
            set_field(&fields, "actual", JsValue::from(*actual as f64))?;
        }
        E::ReactantTemplateIndex { index, count, .. } => {
            set_field(&fields, "index", JsValue::from(*index as f64))?;
            set_field(&fields, "count", JsValue::from(*count as f64))?;
        }
        E::Matching {
            reactant, template, ..
        } => {
            set_field(&fields, "reactant", JsValue::from(*reactant as f64))?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
        }
        E::MatchAtomIndex {
            reactant,
            template,
            atom,
            atom_count,
            ..
        } => {
            set_field(&fields, "reactant", JsValue::from(*reactant as f64))?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
            set_field(&fields, "atom", JsValue::from(*atom as f64))?;
            set_field(&fields, "atomCount", JsValue::from(*atom_count as f64))?;
        }
        E::CombinationLevel { level, count, .. } => {
            set_field(&fields, "level", JsValue::from(*level as f64))?;
            set_field(&fields, "count", JsValue::from(*count as f64))?;
        }
        E::CombinationSize {
            expected, actual, ..
        } => {
            set_field(&fields, "expected", JsValue::from(*expected as f64))?;
            set_field(&fields, "actual", JsValue::from(*actual as f64))?;
        }
        _ => (),
    }
    let cause = match source {
        E::Initialization(cause) => source_error(cause)?,
        _ => source
            .source()
            .map(source_error)
            .transpose()?
            .unwrap_or(JsValue::NULL),
    };
    if !cause.is_null() {
        set_field(&error, "cause", cause.clone())?;
    }
    for key in Object::keys(&fields).iter() {
        let value = Reflect::get(&fields, &key)?;
        Reflect::set(&error, &key, &value)?;
    }
    set_field(
        &error,
        "detail",
        ReactionRunError {
            fields,
            message: source.to_string(),
            kind: kind.into(),
            cause,
        }
        .into(),
    )?;
    Ok(error)
}
#[wasm_bindgen]
pub struct ReactionApplyError {
    fields: Object,
    message: String,
    kind: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl ReactionApplyError {
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter,js_name=reactants,unchecked_return_type="number | null")]
    pub fn reactants(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"reactants".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=products,unchecked_return_type="number | null")]
    pub fn products(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"products".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=atom,unchecked_return_type="number | null")]
    pub fn atom(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"atom".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
}
pub(crate) fn apply_error(source: &ck::ReactionApplyError) -> Result<JsValue, JsValue> {
    use ck::ReactionApplyError as E;
    let kind = match source {
        E::SourceText(..) => "SourceText",
        E::Initialization(..) => "Initialization",
        E::ApplicabilityArity { .. } => "ApplicabilityArity",
        E::AddsProductAtom { .. } => "AddsProductAtom",
        E::Matching(..) => "Matching",
        E::Product(..) => "Product",
        E::Edit(..) => "Edit",
    };
    let error = Error::new(&source.to_string());
    error.set_name("ReactionApplyError");
    let error: JsValue = error.into();
    set_field(&error, "domain", "reaction".into())?;
    set_field(&error, "kind", kind.into())?;
    let fields = Object::new();
    match source {
        E::ApplicabilityArity {
            reactants,
            products,
            ..
        } => {
            set_field(&fields, "reactants", JsValue::from(*reactants as f64))?;
            set_field(&fields, "products", JsValue::from(*products as f64))?;
        }
        E::AddsProductAtom { atom, .. } => {
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
        }
        _ => (),
    }
    let cause = match source {
        E::SourceText(cause) => source_error(cause)?,
        E::Initialization(cause) => source_error(cause)?,
        E::Matching(cause) => source_error(cause)?,
        E::Product(cause) => source_error(cause)?,
        E::Edit(cause) => source_error(cause)?,
        _ => source
            .source()
            .map(source_error)
            .transpose()?
            .unwrap_or(JsValue::NULL),
    };
    if !cause.is_null() {
        set_field(&error, "cause", cause.clone())?;
    }
    for key in Object::keys(&fields).iter() {
        let value = Reflect::get(&fields, &key)?;
        Reflect::set(&error, &key, &value)?;
    }
    set_field(
        &error,
        "detail",
        ReactionApplyError {
            fields,
            message: source.to_string(),
            kind: kind.into(),
            cause,
        }
        .into(),
    )?;
    Ok(error)
}
#[wasm_bindgen]
pub struct ReactionProductError {
    fields: Object,
    message: String,
    kind: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl ReactionProductError {
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter,js_name=atom,unchecked_return_type="number | null")]
    pub fn atom(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"atom".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=key,unchecked_return_type="string | null")]
    pub fn key(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"key".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=stage,unchecked_return_type="string | null")]
    pub fn stage(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"stage".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=detail,unchecked_return_type="string | null")]
    pub fn detail(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"detail".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=reactantAtom,unchecked_return_type="number | null")]
    pub fn reactant_atom(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"reactantAtom".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=productAtom,unchecked_return_type="number | null")]
    pub fn product_atom(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"productAtom".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=bond,unchecked_return_type="number | null")]
    pub fn bond(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"bond".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=rowKind,unchecked_return_type="string | null")]
    pub fn row_kind(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"rowKind".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=index,unchecked_return_type="number | null")]
    pub fn index(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"index".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
}
pub(crate) fn product_error(source: &ck::ReactionProductError) -> Result<JsValue, JsValue> {
    use ck::ReactionProductError as E;
    let kind = match source {
        E::StereoGroup(..) => "StereoGroup",
        E::Coordinate(..) => "Coordinate",
        E::TopologyEdit(..) => "TopologyEdit",
        E::Adjacency(..) => "Adjacency",
        E::RingFinding(..) => "RingFinding",
        E::SourceUInt(..) => "SourceUInt",
        E::Valence(..) => "Valence",
        E::StereoOrder(..) => "StereoOrder",
        E::DoubleBondStereo(..) => "DoubleBondStereo",
        E::CoordinateSelection(..) => "CoordinateSelection",
        E::CarrierIdentity(..) => "CarrierIdentity",
        E::MoleculeProperty(..) => "MoleculeProperty",
        E::AtomProperty(..) => "AtomProperty",
        E::BondValue(..) => "BondValue",
        E::TemplateProperty(..) => "TemplateProperty",
        E::StereoGetter { .. } => "StereoGetter",
        E::PropertyInt { .. } => "PropertyInt",
        E::PropertyUInt { .. } => "PropertyUInt",
        E::MissingProperty { .. } => "MissingProperty",
        E::Invariant { .. } => "Invariant",
        E::RowOverflow { .. } => "RowOverflow",
        E::QueryReactantBond { .. } => "QueryReactantBond",
    };
    let error = Error::new(&source.to_string());
    error.set_name("ReactionProductError");
    let error: JsValue = error.into();
    set_field(&error, "domain", "reaction".into())?;
    set_field(&error, "kind", kind.into())?;
    let fields = Object::new();
    match source {
        E::StereoGetter { atom, .. } => {
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
        }
        E::PropertyInt { atom, key, .. } => {
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
            set_field(&fields, "key", (*key).into())?;
        }
        E::PropertyUInt { atom, key, .. } => {
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
            set_field(&fields, "key", (*key).into())?;
        }
        E::MissingProperty { atom, key, .. } => {
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
            set_field(&fields, "key", (*key).into())?;
        }
        E::Invariant {
            stage,
            detail,
            reactant_atom,
            product_atom,
            bond,
            ..
        } => {
            set_field(&fields, "stage", (*stage).into())?;
            set_field(&fields, "detail", (*detail).into())?;
            set_field(
                &fields,
                "reactantAtom",
                reactant_atom.map_or(JsValue::NULL, |v| JsValue::from(v as f64)),
            )?;
            set_field(
                &fields,
                "productAtom",
                product_atom.map_or(JsValue::NULL, |v| JsValue::from(v as f64)),
            )?;
            set_field(
                &fields,
                "bond",
                bond.map_or(JsValue::NULL, |v| JsValue::from(v.index() as f64)),
            )?;
        }
        E::RowOverflow { kind, index, .. } => {
            set_field(&fields, "rowKind", (*kind).into())?;
            set_field(&fields, "index", JsValue::from(*index as f64))?;
        }
        E::QueryReactantBond { bond, .. } => {
            set_field(&fields, "bond", JsValue::from(bond.index() as f64))?;
        }
        _ => (),
    }
    let cause = match source {
        E::Coordinate(cause) => source_error(cause)?,
        E::TopologyEdit(cause) => source_error(cause)?,
        E::Adjacency(cause) => source_error(cause)?,
        E::RingFinding(cause) => source_error(cause)?,
        E::SourceUInt(cause) => source_error(cause)?,
        E::Valence(cause) => source_error(cause)?,
        E::StereoOrder(cause) => source_error(cause)?,
        E::DoubleBondStereo(cause) => source_error(cause)?,
        E::CoordinateSelection(cause) => source_error(cause)?,
        E::CarrierIdentity(cause) => source_error(cause)?,
        E::MoleculeProperty(cause) => source_error(cause)?,
        E::AtomProperty(cause) => source_error(cause)?,
        E::BondValue(cause) => source_error(cause)?,
        E::TemplateProperty(cause) => source_error(cause)?,
        _ => source
            .source()
            .map(source_error)
            .transpose()?
            .unwrap_or(JsValue::NULL),
    };
    if !cause.is_null() {
        set_field(&error, "cause", cause.clone())?;
    }
    for key in Object::keys(&fields).iter() {
        let value = Reflect::get(&fields, &key)?;
        Reflect::set(&error, &key, &value)?;
    }
    set_field(
        &error,
        "detail",
        ReactionProductError {
            fields,
            message: source.to_string(),
            kind: kind.into(),
            cause,
        }
        .into(),
    )?;
    Ok(error)
}
#[wasm_bindgen]
pub struct ReactionWriteError {
    fields: Object,
    message: String,
    kind: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl ReactionWriteError {
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter,js_name=role,unchecked_return_type="ReactionRole | null")]
    pub fn role(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"role".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=template,unchecked_return_type="number | null")]
    pub fn template(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"template".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=expected,unchecked_return_type="number | null")]
    pub fn expected(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"expected".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=actual,unchecked_return_type="number | null")]
    pub fn actual(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"actual".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
}
pub(crate) fn write_error(source: &ck::ReactionWriteError) -> Result<JsValue, JsValue> {
    use ck::ReactionWriteError as E;
    let kind = match source {
        E::Cx(..) => "Cx",
        E::Template { .. } => "Template",
        E::Connectivity { .. } => "Connectivity",
        E::CoordinateSelectionArity { .. } => "CoordinateSelectionArity",
    };
    let error = Error::new(&source.to_string());
    error.set_name("ReactionWriteError");
    let error: JsValue = error.into();
    set_field(&error, "domain", "reaction".into())?;
    set_field(&error, "kind", kind.into())?;
    let fields = Object::new();
    match source {
        E::Template { role, template, .. } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
        }
        E::Connectivity { role, template, .. } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
        }
        E::CoordinateSelectionArity {
            expected, actual, ..
        } => {
            set_field(&fields, "expected", JsValue::from(*expected as f64))?;
            set_field(&fields, "actual", JsValue::from(*actual as f64))?;
        }
        _ => (),
    }
    let cause = match source {
        E::Cx(cause) => source_error(cause)?,
        _ => source
            .source()
            .map(source_error)
            .transpose()?
            .unwrap_or(JsValue::NULL),
    };
    if !cause.is_null() {
        set_field(&error, "cause", cause.clone())?;
    }
    for key in Object::keys(&fields).iter() {
        let value = Reflect::get(&fields, &key)?;
        Reflect::set(&error, &key, &value)?;
    }
    set_field(
        &error,
        "detail",
        ReactionWriteError {
            fields,
            message: source.to_string(),
            kind: kind.into(),
            cause,
        }
        .into(),
    )?;
    Ok(error)
}
#[wasm_bindgen]
pub struct ReactionValidationError {
    fields: Object,
    message: String,
    kind: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl ReactionValidationError {
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter,js_name=role,unchecked_return_type="ReactionRole | null")]
    pub fn role(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"role".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=template,unchecked_return_type="number | null")]
    pub fn template(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"template".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=atom,unchecked_return_type="number | null")]
    pub fn atom(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"atom".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=property,unchecked_return_type="string | null")]
    pub fn property(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"property".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=value,unchecked_return_type="number | null")]
    pub fn value(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"value".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=atomCount,unchecked_return_type="number | null")]
    pub fn atom_count(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"atomCount".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
    #[wasm_bindgen(getter,js_name=map,unchecked_return_type="number | null")]
    pub fn map(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"map".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
}
pub(crate) fn validation_error(source: &ck::ReactionValidationError) -> Result<JsValue, JsValue> {
    use ck::ReactionValidationError as E;
    let kind = match source {
        E::Property { .. } => "Property",
        E::MapOverflow { .. } => "MapOverflow",
        E::AtomBounds { .. } => "AtomBounds",
        E::MissingReactingAtom { .. } => "MissingReactingAtom",
        E::QueryArithmetic { .. } => "QueryArithmetic",
        E::Annotation { .. } => "Annotation",
    };
    let error = Error::new(&source.to_string());
    error.set_name("ReactionValidationError");
    let error: JsValue = error.into();
    set_field(&error, "domain", "reaction".into())?;
    set_field(&error, "kind", kind.into())?;
    let fields = Object::new();
    match source {
        E::Property {
            role,
            template,
            atom,
            property,
            ..
        } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
            set_field(&fields, "property", (*property).into())?;
        }
        E::MapOverflow {
            role,
            template,
            atom,
            property,
            value,
            ..
        } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
            set_field(&fields, "property", (*property).into())?;
            set_field(&fields, "value", JsValue::from(*value as f64))?;
        }
        E::AtomBounds {
            role,
            template,
            atom,
            atom_count,
            ..
        } => {
            set_field(
                &fields,
                "role",
                JsValue::from(ReactionRole::from(*role) as u32),
            )?;
            set_field(&fields, "template", JsValue::from(*template as f64))?;
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
            set_field(&fields, "atomCount", JsValue::from(*atom_count as f64))?;
        }
        E::MissingReactingAtom {
            template,
            atom,
            map,
            ..
        } => {
            set_field(&fields, "template", JsValue::from(*template as f64))?;
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
            set_field(&fields, "map", JsValue::from(*map as f64))?;
        }
        E::QueryArithmetic { template, atom, .. } => {
            set_field(&fields, "template", JsValue::from(*template as f64))?;
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
        }
        E::Annotation {
            template,
            atom,
            property,
            ..
        } => {
            set_field(&fields, "template", JsValue::from(*template as f64))?;
            set_field(&fields, "atom", JsValue::from(atom.index() as f64))?;
            set_field(&fields, "property", (*property).into())?;
        }
        _ => (),
    }
    let cause = match source {
        _ => source
            .source()
            .map(source_error)
            .transpose()?
            .unwrap_or(JsValue::NULL),
    };
    if !cause.is_null() {
        set_field(&error, "cause", cause.clone())?;
    }
    for key in Object::keys(&fields).iter() {
        let value = Reflect::get(&fields, &key)?;
        Reflect::set(&error, &key, &value)?;
    }
    set_field(
        &error,
        "detail",
        ReactionValidationError {
            fields,
            message: source.to_string(),
            kind: kind.into(),
            cause,
        }
        .into(),
    )?;
    Ok(error)
}
#[wasm_bindgen]
pub struct ReactionInitializationError {
    fields: Object,
    message: String,
    kind: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl ReactionInitializationError {
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter,js_name=report,unchecked_return_type="ReactionValidationReport | null")]
    pub fn report(&self) -> Result<JsValue, JsValue> {
        let v = Reflect::get(&self.fields, &"report".into())?;
        Ok(if v.is_undefined() { JsValue::NULL } else { v })
    }
}
pub(crate) fn initialization_error(
    source: &ck::ReactionInitializationError,
) -> Result<JsValue, JsValue> {
    use ck::ReactionInitializationError as E;
    let kind = match source {
        E::Validation { .. } => "Validation",
        E::Invalid { .. } => "Invalid",
    };
    let error = Error::new(&source.to_string());
    error.set_name("ReactionInitializationError");
    let error: JsValue = error.into();
    set_field(&error, "domain", "reaction".into())?;
    set_field(&error, "kind", kind.into())?;
    let fields = Object::new();
    match source {
        E::Invalid { report, .. } => {
            set_field(
                &fields,
                "report",
                JsValue::from(ReactionValidationReport {
                    inner: report.clone(),
                }),
            )?;
        }
        _ => (),
    }
    let cause = match source {
        _ => source
            .source()
            .map(source_error)
            .transpose()?
            .unwrap_or(JsValue::NULL),
    };
    if !cause.is_null() {
        set_field(&error, "cause", cause.clone())?;
    }
    for key in Object::keys(&fields).iter() {
        let value = Reflect::get(&fields, &key)?;
        Reflect::set(&error, &key, &value)?;
    }
    set_field(
        &error,
        "detail",
        ReactionInitializationError {
            fields,
            message: source.to_string(),
            kind: kind.into(),
            cause,
        }
        .into(),
    )?;
    Ok(error)
}
