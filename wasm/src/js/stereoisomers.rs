//! Typed options, arbitrary-width random sources and owned lazy enumeration.
use crate::{
    Molecule,
    alignment_values::{set, source_error},
    host_values::*,
};
use cosmolkit_wasm::rust as ck;
use js_sys::{BigInt, Function};
use std::{
    error::Error as RustError,
    sync::{Arc, Mutex},
};
use wasm_bindgen::{JsCast, prelude::*};
type Failure = Arc<Mutex<Option<JsValue>>>;
fn operation<T>(v: Result<T, ck::OperationError>) -> Result<T, JsValue> {
    v.map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
}
fn completed<T>(
    failure: Option<&Failure>,
    result: Result<T, ck::OperationError>,
) -> Result<T, JsValue> {
    if let Some(failure) = failure {
        if let Some(error) = failure.lock().expect("random callback error mutex").take() {
            return Err(error);
        }
    }
    operation(result)
}
fn decimal(value: &JsValue) -> Result<String, JsValue> {
    if !value.is_bigint() {
        return Err(type_error("expected bigint"));
    }
    Ok(value
        .clone()
        .unchecked_into::<BigInt>()
        .to_string(10)?
        .into())
}
fn bits(function: &Function, slot: &Failure, width: usize) -> Result<String, String> {
    function
        .call1(&JsValue::UNDEFINED, &JsValue::from(width as u32))
        .and_then(|value| {
            let text = decimal(&value)?;
            if text.starts_with('-') {
                Err(js_sys::RangeError::new("random bits must be nonnegative").into())
            } else {
                Ok(text)
            }
        })
        .map_err(|error| {
            let mut saved = slot.lock().expect("random callback error mutex");
            if saved.is_none() {
                *saved = Some(error);
            }
            "JavaScript random bit provider failed".into()
        })
}
#[wasm_bindgen]
#[derive(Clone)]
pub struct StereoisomerRandomSource {
    inner: ck::StereoisomerRandomSource,
    failure: Option<Failure>,
}
#[wasm_bindgen]
impl StereoisomerRandomSource {
    #[wasm_bindgen(js_name=fromIntegerSeed)]
    pub fn from_integer_seed(
        #[wasm_bindgen(unchecked_param_type = "bigint")] value: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::StereoisomerRandomSource::from_integer_seed(
                decimal(&value)?
                    .parse()
                    .map_err(|_| type_error("integer seed"))?,
            ),
            failure: None,
        })
    }
    #[wasm_bindgen(js_name=fromRandomBits)]
    pub fn from_random_bits(
        #[wasm_bindgen(unchecked_param_type = "(bitCount: number) => bigint")] function: Function,
    ) -> Self {
        let failure = Arc::new(Mutex::new(None));
        let slot = Arc::clone(&failure);
        Self {
            inner: ck::StereoisomerRandomSource::from_random_bits(Box::new(move |width| {
                bits(&function, &slot, width)?
                    .parse()
                    .map_err(|e| format!("invalid random bits: {e}"))
            })),
            failure: Some(failure),
        }
    }
}
#[wasm_bindgen(
    inline_js = "export function visitStereoisomerRandomSource(value, visit) { visit(value); }"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitStereoisomerRandomSource)]
    fn visit_source(
        value: &JsValue,
        visit: &mut dyn FnMut(&StereoisomerRandomSource),
    ) -> Result<(), JsValue>;
}
fn random_source(value: &JsValue) -> Result<Option<StereoisomerRandomSource>, JsValue> {
    if value.is_undefined() || value.is_null() {
        return Ok(None);
    }
    let mut result = None;
    visit_source(value, &mut |source| result = Some(source.clone()))?;
    result.map(Some).ok_or_else(|| type_error("randomSource"))
}
#[wasm_bindgen]
#[derive(Clone)]
pub struct StereoisomerOptions {
    inner: ck::StereoisomerOptions,
    source: Option<StereoisomerRandomSource>,
}
#[cosmolkit_wasm::javascript_options("StereoisomerOptionsConfiguration", clone_snapshot)]
#[wasm_bindgen]
impl StereoisomerOptions {
    #[wasm_bindgen(getter,js_name=tryEmbedding)]
    pub fn try_embedding(&self) -> bool {
        self.inner.try_embedding()
    }
    #[wasm_bindgen(setter,js_name=tryEmbedding)]
    pub fn set_try_embedding(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner
            .set_try_embedding(bool_value(&value, "tryEmbedding")?);
        Ok(())
    }
    #[wasm_bindgen(js_name=setTryEmbedding)]
    pub fn set_try_embedding_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_try_embedding(value)
    }
    #[wasm_bindgen(getter,js_name=onlyUnassigned)]
    pub fn only_unassigned(&self) -> bool {
        self.inner.only_unassigned()
    }
    #[wasm_bindgen(setter,js_name=onlyUnassigned)]
    pub fn set_only_unassigned(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner
            .set_only_unassigned(bool_value(&value, "onlyUnassigned")?);
        Ok(())
    }
    #[wasm_bindgen(js_name=setOnlyUnassigned)]
    pub fn set_only_unassigned_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_only_unassigned(value)
    }
    #[wasm_bindgen(getter,js_name=maxIsomers)]
    pub fn max_isomers(&self) -> usize {
        self.inner.max_isomers()
    }
    #[wasm_bindgen(setter,js_name=maxIsomers)]
    pub fn set_max_isomers(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner
            .set_max_isomers(usize_value(&value, "maxIsomers")?);
        Ok(())
    }
    #[wasm_bindgen(js_name=setMaxIsomers)]
    pub fn set_max_isomers_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_max_isomers(value)
    }
    #[wasm_bindgen(getter,js_name=unique)]
    pub fn unique(&self) -> bool {
        self.inner.unique()
    }
    #[wasm_bindgen(setter,js_name=unique)]
    pub fn set_unique(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.set_unique(bool_value(&value, "unique")?);
        Ok(())
    }
    #[wasm_bindgen(js_name=setUnique)]
    pub fn set_unique_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_unique(value)
    }
    #[wasm_bindgen(getter,js_name=onlyStereoGroups)]
    pub fn only_stereo_groups(&self) -> bool {
        self.inner.only_stereo_groups()
    }
    #[wasm_bindgen(setter,js_name=onlyStereoGroups)]
    pub fn set_only_stereo_groups(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner
            .set_only_stereo_groups(bool_value(&value, "onlyStereoGroups")?);
        Ok(())
    }
    #[wasm_bindgen(js_name=setOnlyStereoGroups)]
    pub fn set_only_stereo_groups_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_only_stereo_groups(value)
    }
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] try_embedding: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] only_unassigned: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_isomers: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "StereoisomerRandomSource | null")]
        random_source: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] unique: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] only_stereo_groups: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::StereoisomerOptions::default();
        if !try_embedding.is_undefined() {
            inner.set_try_embedding(bool_value(&try_embedding, "tryEmbedding")?);
        }
        if !only_unassigned.is_undefined() {
            inner.set_only_unassigned(bool_value(&only_unassigned, "onlyUnassigned")?);
        }
        if !max_isomers.is_undefined() {
            inner.set_max_isomers(usize_value(&max_isomers, "maxIsomers")?);
        }
        if !unique.is_undefined() {
            inner.set_unique(bool_value(&unique, "unique")?);
        }
        if !only_stereo_groups.is_undefined() {
            inner.set_only_stereo_groups(bool_value(&only_stereo_groups, "onlyStereoGroups")?);
        }
        let source = self::random_source(&random_source)?;
        inner.set_random_source(source.as_ref().map(|v| v.inner.clone()));
        Ok(Self { inner, source })
    }
    #[wasm_bindgen(js_name=default)]
    pub fn defaults() -> Self {
        Self {
            inner: ck::StereoisomerOptions::default(),
            source: None,
        }
    }
    #[wasm_bindgen(getter,js_name=randomSource,unchecked_return_type="StereoisomerRandomSource | null")]
    pub fn random_source(&self) -> JsValue {
        self.source.clone().map_or(JsValue::NULL, Into::into)
    }
    #[wasm_bindgen(setter,js_name=randomSource)]
    pub fn set_random_source(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "StereoisomerRandomSource | null")] value: JsValue,
    ) -> Result<(), JsValue> {
        let source = random_source(&value)?;
        self.inner
            .set_random_source(source.as_ref().map(|v| v.inner.clone()));
        self.source = source;
        Ok(())
    }
    #[wasm_bindgen(js_name=setRandomSource)]
    pub fn set_random_source_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "StereoisomerRandomSource | null")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_random_source(value)
    }
}
#[wasm_bindgen]
pub struct StereoisomerIterator {
    inner: ck::StereoisomerIterator,
    failure: Option<Failure>,
}
#[wasm_bindgen]
impl StereoisomerIterator {
    #[wasm_bindgen(unchecked_return_type = "Molecule | null")]
    pub fn next(&mut self) -> Result<JsValue, JsValue> {
        let result = self.inner.next().transpose();
        if let Some(slot) = &self.failure {
            if let Some(error) = slot.lock().expect("random callback error mutex").take() {
                return Err(error);
            }
        }
        operation(result).map(|v| {
            v.map_or(JsValue::NULL, |inner| {
                Molecule {
                    inner: cosmolkit_wasm::Molecule::from(inner).into(),
                }
                .into()
            })
        })
    }
    #[wasm_bindgen(getter,js_name=yieldedCount)]
    pub fn yielded_count(&self) -> usize {
        self.inner.yielded_count()
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=stereoisomerCount,unchecked_return_type="bigint")]
    pub fn stereoisomer_count(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "StereoisomerOptions | StereoisomerOptionsConfiguration"
        )]
        options: JsValue,
    ) -> Result<JsValue, JsValue> {
        self.stereoisomer_count_with_options(&StereoisomerOptions::from_js_configuration(&options)?)
    }
    #[wasm_bindgen(js_name=stereoisomerCountWithOptions,unchecked_return_type="bigint")]
    pub fn stereoisomer_count_with_options(
        &self,
        options: &StereoisomerOptions,
    ) -> Result<JsValue, JsValue> {
        let value = self
            .inner
            .stereoisomer_count_with_options(&options.inner)
            .map_err(|e| enumeration_source(&e).unwrap_or_else(|e| e))?;
        BigInt::new(&value.into())
            .map(Into::into)
            .map_err(Into::into)
    }
    #[wasm_bindgen(js_name=enumerateStereoisomers)]
    pub fn enumerate_stereoisomers(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "StereoisomerOptions | StereoisomerOptionsConfiguration"
        )]
        options: JsValue,
    ) -> Result<StereoisomerIterator, JsValue> {
        self.enumerate_stereoisomers_with_options(&StereoisomerOptions::from_js_configuration(
            &options,
        )?)
    }
    #[wasm_bindgen(js_name=enumerateStereoisomersWithOptions)]
    pub fn enumerate_stereoisomers_with_options(
        &self,
        options: &StereoisomerOptions,
    ) -> Result<StereoisomerIterator, JsValue> {
        completed(
            options.source.as_ref().and_then(|v| v.failure.as_ref()),
            self.inner
                .enumerate_stereoisomers_with_options(&options.inner),
        )
        .map(|inner| StereoisomerIterator {
            inner,
            failure: options.source.as_ref().and_then(|v| v.failure.clone()),
        })
    }
    #[wasm_bindgen(js_name=enumerateStereoisomersWithRandomBits)]
    pub fn enumerate_stereoisomers_with_random_bits(
        &self,
        options: &StereoisomerOptions,
        #[wasm_bindgen(unchecked_param_type = "(bitCount: number) => bigint")] callback: Function,
    ) -> Result<StereoisomerIterator, JsValue> {
        let failure = Arc::new(Mutex::new(None));
        let slot = Arc::clone(&failure);
        completed(
            Some(&failure),
            self.inner.enumerate_stereoisomers_with_random_bits(
                &options.inner,
                Box::new(move |width| bits(&callback, &slot, width)),
            ),
        )
        .map(|inner| StereoisomerIterator {
            inner,
            failure: Some(failure),
        })
    }
}
#[wasm_bindgen]
pub struct EnumerationError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
pub struct EnumerationRunError {
    inner: ck::EnumerationRunError,
}
#[wasm_bindgen]
impl EnumerationRunError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "stereoisomers".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        "Enumeration".into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
}
pub(crate) fn run_source(source: &ck::EnumerationRunError) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("EnumerationRunError");
    let error: JsValue = error.into();
    set(&error, "domain", "stereoisomers".into())?;
    set(&error, "kind", "Enumeration".into())?;
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    set(
        &error,
        "detail",
        EnumerationRunError {
            inner: source.clone(),
        }
        .into(),
    )?;
    Ok(error)
}
#[wasm_bindgen]
impl EnumerationError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "stereoisomers".into()
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
}
pub(crate) fn enumeration_source(source: &ck::EnumerationError) -> Result<JsValue, JsValue> {
    use ck::EnumerationError as E;
    let kind = match source {
        E::Valence(..) => "Valence",
        E::Rings(..) => "Rings",
        E::Stereo(..) => "Stereo",
        E::DoubleBond(..) => "DoubleBond",
        E::SmilesWrite(..) => "SmilesWrite",
        E::PotentialStereo(..) => "PotentialStereo",
        E::Hydrogens(..) => "Hydrogens",
        E::Embedding(..) => "Embedding",
        E::AtomProperty(..) => "AtomProperty",
        E::BondProperty(..) => "BondProperty",
        E::MoleculeProperty(..) => "MoleculeProperty",
        E::Topology(..) => "Topology",
        E::Coordinates(..) => "Coordinates",
        E::InvalidFlipperAtom { .. } => "InvalidFlipperAtom",
        E::InvalidFlipperBond { .. } => "InvalidFlipperBond",
        E::RandomBitsSource { .. } => "RandomBitsSource",
        E::RandomSourcePoisoned(..) => "RandomSourcePoisoned",
        E::EmbeddingCoordinates => "EmbeddingCoordinates",
    };
    let message = source.to_string();
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&message);
    error.set_name("EnumerationError");
    let error: JsValue = error.into();
    set(&error, "domain", "stereoisomers".into())?;
    set(&error, "kind", kind.into())?;
    if !cause.is_null() {
        set(&error, "cause", cause.clone())?;
    }
    match source {
        E::InvalidFlipperAtom { atom } => set(&error, "atom", (atom.index() as f64).into())?,
        E::InvalidFlipperBond { bond } => set(&error, "bond", (bond.index() as f64).into())?,
        E::RandomBitsSource { bit_count, message } => {
            set(&error, "bitCount", (*bit_count as f64).into())?;
            set(&error, "message", message.clone().into())?;
        }
        _ => {}
    }
    set(
        &error,
        "detail",
        EnumerationError {
            kind: kind.into(),
            message,
            cause,
        }
        .into(),
    )?;
    Ok(error)
}
