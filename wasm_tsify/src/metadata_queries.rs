//! Native periodic, atom-metadata and fingerprint explanation values.
use crate::{
    Molecule,
    atoms::{Element, ElementInput},
    fields::Field,
    input,
};
use cosmolkit as ck;
use tsify::{Ts, Tsify};
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
pub struct ElementInfo(ck::ElementInfo);
#[wasm_bindgen]
impl ElementInfo {
    #[wasm_bindgen(getter)]
    pub fn element(&self) -> Element {
        Element(self.0.element)
    }
    #[wasm_bindgen(getter)]
    pub fn symbol(&self) -> String {
        self.0.symbol.into()
    }
    #[wasm_bindgen(getter,js_name=atomicNumber)]
    pub fn atomic_number(&self) -> u8 {
        self.0.atomic_number
    }
    #[wasm_bindgen(getter)]
    pub fn period(&self) -> u8 {
        self.0.period
    }
    #[wasm_bindgen(getter,js_name=outerElectrons)]
    pub fn outer_electrons(&self) -> i32 {
        self.0.outer_electrons
    }
    #[wasm_bindgen(getter)]
    pub fn valences(&self) -> Vec<i32> {
        self.0.valences.to_vec()
    }
    #[wasm_bindgen(getter)]
    pub fn rb0(&self) -> f64 {
        self.0.rb0
    }
    #[wasm_bindgen(getter,js_name=atomicWeight)]
    pub fn atomic_weight(&self) -> f64 {
        self.0.atomic_weight
    }
}
#[wasm_bindgen(js_name=elementInfo)]
pub fn element_info(element: ElementInput) -> Result<ElementInfo, JsValue> {
    Ok(ElementInfo(ck::element_info(ck::Element::input(element)?)))
}
#[wasm_bindgen]
pub struct AtomMetadata(pub(crate) ck::AtomMetadata);
#[wasm_bindgen]
impl AtomMetadata {
    #[wasm_bindgen(constructor)]
    pub fn new(
        degree: Ts<input::IndexInput>,
        explicit_valence: Ts<input::I32Input>,
        implicit_hydrogens: Ts<input::I32Input>,
        total_hydrogens: Ts<input::I32Input>,
        total_valence: Ts<input::I32Input>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::AtomMetadata {
            degree: input::read::<usize>(degree)?,
            explicit_valence: input::read::<i32>(explicit_valence)?,
            implicit_hydrogens: input::read::<i32>(implicit_hydrogens)?,
            total_hydrogens: input::read::<i32>(total_hydrogens)?,
            total_valence: input::read::<i32>(total_valence)?,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn degree(&self) -> usize {
        self.0.degree
    }
    #[wasm_bindgen(getter,js_name=explicitValence)]
    pub fn explicit_valence(&self) -> i32 {
        self.0.explicit_valence
    }
    #[wasm_bindgen(getter,js_name=implicitHydrogens)]
    pub fn implicit_hydrogens(&self) -> i32 {
        self.0.implicit_hydrogens
    }
    #[wasm_bindgen(getter,js_name=totalHydrogens)]
    pub fn total_hydrogens(&self) -> i32 {
        self.0.total_hydrogens
    }
    #[wasm_bindgen(getter,js_name=totalValence)]
    pub fn total_valence(&self) -> i32 {
        self.0.total_valence
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=hasValenceViolation)]
    pub fn has_valence_violation(&self, atom: Ts<input::IndexInput>) -> Result<bool, JsValue> {
        self.0
            .has_valence_violation(ck::AtomId::new(input::read::<usize>(atom)?))
            .map_err(crate::valence_errors::ValenceError::from_rust)
    }
    #[wasm_bindgen(js_name=bondPropertyString)]
    pub fn bond_property_string(
        &self,
        bond: Ts<input::IndexInput>,
        key: Ts<input::TextInput>,
    ) -> Result<Option<String>, JsValue> {
        self.0
            .bond_property_string(
                ck::BondId::new(input::read::<usize>(bond)?),
                &input::read::<String>(key)?,
            )
            .map_err(crate::property_errors::PropertyStringError::from_rust)?
            .as_ref()
            .map(crate::property_errors::text)
            .transpose()
    }
    #[wasm_bindgen(js_name=atomMetadata)]
    pub fn atom_metadata(
        &self,
        recalculate: Option<Ts<input::BoolInput>>,
    ) -> Result<Vec<AtomMetadata>, JsValue> {
        self.0
            .atom_metadata(
                recalculate
                    .map(input::read::<bool>)
                    .transpose()?
                    .unwrap_or(true),
            )
            .map(|v| v.into_iter().map(AtomMetadata).collect())
            .map_err(crate::valence_errors::ValenceError::from_rust)
    }
}
crate::errors::enum_error!(AtomCodeExplanationError,AtomCodeExplanationErrorKind,ck::AtomCodeExplanationError,{UnknownChirality=>ck::AtomCodeExplanationError::UnknownChirality{..}},{#[wasm_bindgen(getter)]pub fn code(&self)->u8{match self.inner.as_ref(){ck::AtomCodeExplanationError::UnknownChirality{code}=>*code}}});
#[wasm_bindgen]
pub struct AtomCodeExplanation(ck::AtomCodeExplanation);
#[wasm_bindgen]
impl AtomCodeExplanation {
    #[wasm_bindgen(js_name=fromCode)]
    pub fn from_code(
        code: Ts<input::U64Input>,
        branch_subtract: Option<Ts<input::I64Input>>,
        include_chirality: Option<Ts<input::BoolInput>>,
    ) -> Result<Self, JsValue> {
        ck::AtomCodeExplanation::from_code(
            input::read::<u64>(code)?,
            branch_subtract
                .map(input::read::<i64>)
                .transpose()?
                .unwrap_or(0),
            include_chirality
                .map(input::read::<bool>)
                .transpose()?
                .unwrap_or(false),
        )
        .map(Self)
        .map_err(AtomCodeExplanationError::from_rust)
    }
    #[wasm_bindgen(getter)]
    pub fn symbol(&self) -> String {
        self.0.symbol().into()
    }
    #[wasm_bindgen(getter,js_name=branchCount)]
    pub fn branch_count(&self) -> u32 {
        self.0.branch_count()
    }
    #[wasm_bindgen(getter,js_name=piElectrons)]
    pub fn pi_electrons(&self) -> u32 {
        self.0.pi_electrons()
    }
    #[wasm_bindgen(getter)]
    pub fn chirality(&self) -> Option<String> {
        self.0.chirality().map(Into::into)
    }
}
#[derive(serde::Serialize, Tsify)]
#[serde(transparent)]
pub struct ExplainedPath(
    #[tsify(type = "ReadonlyArray<readonly [string, number, number]>")]
    Vec<(&'static str, u32, u32)>,
);
#[wasm_bindgen(js_name=explainPathScore)]
pub fn explain_path_score(
    score: Ts<input::U64Input>,
    size: Option<Ts<input::IndexInput>>,
) -> Result<Ts<ExplainedPath>, JsValue> {
    let output = ExplainedPath(ck::explain_path_score(
        input::read::<u64>(score)?,
        size.map(input::read::<usize>).transpose()?.unwrap_or(4),
    ))
    .into_ts()
    .map_err(|e| JsValue::from(js_sys::Error::new(&e.to_string())))?;
    let raw = output.js_value();
    let rows: &js_sys::Array = raw.unchecked_ref();
    for row in rows.iter() {
        js_sys::Object::freeze(row.unchecked_ref::<js_sys::Object>());
    }
    js_sys::Object::freeze(rows);
    Ok(output)
}
