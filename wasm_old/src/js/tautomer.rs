//! Owned, restricted tautomer snapshots and canonical configuration forwarding.
use crate::{
    Molecule,
    alignment_values::{set, source_error},
    host_values::*,
};
use cosmolkit_wasm::rust as ck;
use js_sys::{Array, Function};
use std::{
    cell::RefCell,
    error::Error as RustError,
    rc::Rc,
    sync::{Arc, Mutex},
};
use wasm_bindgen::{JsCast, prelude::*};

fn operation<T>(result: Result<T, ck::OperationError>) -> Result<T, JsValue> {
    result.map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
}
fn molecule(inner: ck::Molecule) -> Molecule {
    Molecule {
        inner: cosmolkit_wasm::Molecule::from(inner).into(),
    }
}
fn array_ids(ids: impl Iterator<Item = usize>) -> JsValue {
    ids.map(|id| JsValue::from(id as u32))
        .collect::<Array>()
        .into()
}
#[wasm_bindgen(
    inline_js = "export function visitTautomerMolecule(v,f){try{f(v);}catch(cause){throw new TypeError('expected Molecule',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitTautomerMolecule)]
    fn visit_molecule(v: &JsValue, f: &mut dyn FnMut(&Molecule)) -> Result<(), JsValue>;
}
#[wasm_bindgen(js_name=canonicalTautomerFromMolecules)]
pub fn canonical_tautomer_from_molecules(
    #[wasm_bindgen(unchecked_param_type = "Molecule[]")] molecules: JsValue,
    #[wasm_bindgen(unchecked_optional_param_type = "TautomerParams | TautomerParamsOptions")]
    params: JsValue,
) -> Result<Molecule, JsValue> {
    canonical_tautomer_from_molecules_with_params(
        molecules,
        &TautomerParams::from_js_configuration(&params)?,
    )
}
#[wasm_bindgen(js_name=canonicalTautomerFromMoleculesWithParams)]
pub fn canonical_tautomer_from_molecules_with_params(
    #[wasm_bindgen(unchecked_param_type = "Molecule[]")] molecules: JsValue,
    params: &TautomerParams,
) -> Result<Molecule, JsValue> {
    let values = sequence(&molecules, "molecules")?
        .iter()
        .map(|v| {
            let mut out = None;
            visit_molecule(&v, &mut |m: &Molecule| out = Some(Arc::clone(&m.inner)))?;
            out.ok_or_else(|| type_error("molecule"))
        })
        .collect::<Result<Vec<_>, JsValue>>()?;
    params
        .run(|p| cosmolkit_wasm::canonical_tautomer_from_molecules_with_params(&values, p))
        .map(|inner| Molecule {
            inner: inner.into(),
        })
}
fn checked_callback(value: JsValue, name: &str) -> Result<Option<Function>, JsValue> {
    if value.is_null() || value.is_undefined() {
        Ok(None)
    } else {
        value.dyn_into().map(Some).map_err(|_| type_error(name))
    }
}

#[wasm_bindgen]
pub struct TautomerRunError {
    inner: ck::TautomerRunError,
}
#[wasm_bindgen]
impl TautomerRunError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "tautomer".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        run_kind(&self.inner).into()
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
fn run_kind(e: &ck::TautomerRunError) -> &'static str {
    use ck::TautomerRunError as E;
    match e {
        E::Smiles(..) => "Smiles",
        E::Valence(..) => "Valence",
        E::Rings(..) => "Rings",
        E::Kekulize(..) => "Kekulize",
        E::Sanitize(..) => "Sanitize",
        E::LegacyStereo(..) => "LegacyStereo",
        E::QueryCompile(..) => "QueryCompile",
        E::Match(..) => "Match",
        E::MatchContext(..) => "MatchContext",
        E::Coordinates(..) => "Coordinates",
        E::Topology(..) => "Topology",
        E::AtomProperty(..) => "AtomProperty",
        E::MoleculeProperty(..) => "MoleculeProperty",
        E::BondValue(..) => "BondValue",
        E::PropertyString(..) => "PropertyString",
        E::AtomMappingCount { .. } => "AtomMappingCount",
        E::BondMappingCount { .. } => "BondMappingCount",
        E::AtomCountMismatch { .. } => "AtomCountMismatch",
        E::BondCountMismatch { .. } => "BondCountMismatch",
        E::AtomOutOfRange { .. } => "AtomOutOfRange",
        E::BondOutOfRange { .. } => "BondOutOfRange",
        E::EmptyTransformMatch => "EmptyTransformMatch",
        E::EditCount { .. } => "EditCount",
        E::HydrogenCountOutOfRange { .. } => "HydrogenCountOutOfRange",
        E::FormalChargeOutOfRange { .. } => "FormalChargeOutOfRange",
        E::MissingMappedBond { .. } => "MissingMappedBond",
        E::MissingValence { .. } => "MissingValence",
        E::NoCanonicalTautomer => "NoCanonicalTautomer",
        E::SourceCacheCommit => "SourceCacheCommit",
        E::ScoreCacheTargetMismatch => "ScoreCacheTargetMismatch",
        E::Control(..) => "Control",
        E::Callback(..) => "Callback",
    }
}
pub(crate) fn run_source(e: &ck::TautomerRunError) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&e.to_string());
    error.set_name("TautomerRunError");
    let error: JsValue = error.into();
    set(&error, "domain", "tautomer".into())?;
    set(&error, "kind", run_kind(e).into())?;
    use ck::TautomerRunError as E;
    match e {
        E::AtomMappingCount { expected, actual }
        | E::BondMappingCount { expected, actual }
        | E::AtomCountMismatch { expected, actual }
        | E::BondCountMismatch { expected, actual } => {
            set(&error, "expected", (*expected as f64).into())?;
            set(&error, "actual", (*actual as f64).into())?;
        }
        E::AtomOutOfRange { atom, atom_count } => {
            set(&error, "atom", (*atom as f64).into())?;
            set(&error, "atomCount", (*atom_count as f64).into())?;
        }
        E::BondOutOfRange { bond, bond_count } => {
            set(&error, "bond", (*bond as f64).into())?;
            set(&error, "bondCount", (*bond_count as f64).into())?;
        }
        E::EditCount {
            field,
            expected,
            actual,
        } => {
            set(&error, "field", (*field).into())?;
            set(&error, "expected", (*expected as f64).into())?;
            set(&error, "actual", (*actual as f64).into())?;
        }
        E::HydrogenCountOutOfRange { atom, count } => {
            set(&error, "atom", (atom.index() as f64).into())?;
            set(&error, "count", (*count).into())?;
        }
        E::FormalChargeOutOfRange { atom, charge } => {
            set(&error, "atom", (atom.index() as f64).into())?;
            set(&error, "charge", (*charge).into())?;
        }
        E::MissingMappedBond { bond } => set(&error, "bond", (*bond as f64).into())?,
        E::MissingValence { atom } => set(&error, "atom", (atom.index() as f64).into())?,
        E::Control(reason) | E::Callback(reason) => set(&error, "reason", reason.clone().into())?,
        _ => {}
    }
    if let Some(cause) = e.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    set(
        &error,
        "detail",
        TautomerRunError { inner: e.clone() }.into(),
    )?;
    Ok(error)
}
fn run_error(e: ck::TautomerRunError) -> JsValue {
    run_source(&e).unwrap_or_else(|e| e)
}
#[wasm_bindgen]
pub struct TautomerCatalogError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl TautomerCatalogError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "tautomer_catalog".into()
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
pub(crate) fn catalog_source(e: &ck::TautomerCatalogError) -> Result<JsValue, JsValue> {
    use ck::TautomerCatalogError as E;
    let kind = match e {
        E::BadInputFile { .. } => "BadInputFile",
        E::BadStreamContents(..) => "BadStreamContents",
        E::InvalidUtf8(..) => "InvalidUtf8",
        E::Transform(..) => "Transform",
        E::TransformIndexOutOfRange { .. } => "TransformIndexOutOfRange",
        E::DeserializationUnderConstruction => "DeserializationUnderConstruction",
    };
    let message = e.to_string();
    let cause = e.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&message);
    error.set_name("TautomerCatalogError");
    let error: JsValue = error.into();
    set(&error, "domain", "tautomer_catalog".into())?;
    set(&error, "kind", kind.into())?;
    if !cause.is_null() {
        set(&error, "cause", cause.clone())?;
    }
    match e {
        E::BadInputFile { path, .. } => {
            set(&error, "path", path.to_string_lossy().to_string().into())?
        }
        E::TransformIndexOutOfRange { index, len } => {
            set(&error, "index", (*index as f64).into())?;
            set(&error, "len", (*len as f64).into())?;
        }
        _ => {}
    }
    set(
        &error,
        "detail",
        TautomerCatalogError {
            kind: kind.into(),
            message,
            cause,
        }
        .into(),
    )?;
    Ok(error)
}
fn catalog_error(e: ck::TautomerCatalogError) -> JsValue {
    catalog_source(&e).unwrap_or_else(|e| e)
}

#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum TautomerEnumerationStatus {
    Completed,
    MaxTautomersReached,
    MaxTransformsReached,
    Canceled,
}
fn status(value: ck::TautomerEnumerationStatus) -> TautomerEnumerationStatus {
    match value {
        ck::TautomerEnumerationStatus::Completed => TautomerEnumerationStatus::Completed,
        ck::TautomerEnumerationStatus::MaxTautomersReached => {
            TautomerEnumerationStatus::MaxTautomersReached
        }
        ck::TautomerEnumerationStatus::MaxTransformsReached => {
            TautomerEnumerationStatus::MaxTransformsReached
        }
        ck::TautomerEnumerationStatus::Canceled => TautomerEnumerationStatus::Canceled,
    }
}

#[wasm_bindgen]
#[derive(Clone)]
pub struct TautomerScoreTerm {
    inner: ck::TautomerScoreTerm,
}
#[wasm_bindgen]
impl TautomerScoreTerm {
    #[wasm_bindgen(constructor)]
    pub fn new(name: String, smarts: String, score: JsValue) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::TautomerScoreTerm::new(name, smarts, i32_value(&score, "score")?),
        })
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.inner.name().into()
    }
    #[wasm_bindgen(getter)]
    pub fn smarts(&self) -> String {
        self.inner.smarts().into()
    }
    #[wasm_bindgen(getter)]
    pub fn score(&self) -> i32 {
        self.inner.score()
    }
}
#[wasm_bindgen(
    inline_js = "export function visitTautomerScoreTerm(value, visit) { visit(value); }"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitTautomerScoreTerm)]
    fn visit_term(
        value: &JsValue,
        visit: &mut dyn FnMut(&TautomerScoreTerm),
    ) -> Result<(), JsValue>;
}
fn term(value: &JsValue) -> Result<ck::TautomerScoreTerm, JsValue> {
    let mut out = None;
    visit_term(value, &mut |term| out = Some(term.inner.clone()))?;
    out.ok_or_else(|| type_error("TautomerScoreTerm"))
}
#[wasm_bindgen(js_name=defaultTautomerScoreTerms,unchecked_return_type="TautomerScoreTerm[]")]
pub fn default_terms() -> JsValue {
    ck::default_tautomer_score_terms()
        .iter()
        .cloned()
        .map(|inner| JsValue::from(TautomerScoreTerm { inner }))
        .collect::<Array>()
        .into()
}

#[wasm_bindgen]
#[derive(Clone)]
pub struct TautomerScoreParams {
    inner: ck::TautomerScoreParams,
}
#[cosmolkit_wasm::javascript_options("TautomerScoreParamsOptions", clone_snapshot)]
#[wasm_bindgen]
impl TautomerScoreParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "TautomerScoreTerm[] | null")]
        terms: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::TautomerScoreParams {
                terms: if terms.is_null() || terms.is_undefined() {
                    None
                } else {
                    Some(
                        sequence(&terms, "terms")?
                            .iter()
                            .map(|v| term(&v))
                            .collect::<Result<_, _>>()?,
                    )
                },
            },
        })
    }
    #[wasm_bindgen(getter, unchecked_return_type = "TautomerScoreTerm[] | null")]
    pub fn terms(&self) -> JsValue {
        self.inner.terms.as_ref().map_or(JsValue::NULL, |v| {
            v.iter()
                .cloned()
                .map(|inner| JsValue::from(TautomerScoreTerm { inner }))
                .collect::<Array>()
                .into()
        })
    }
}
#[wasm_bindgen]
pub struct TautomerScore {
    inner: ck::TautomerScore,
}
#[wasm_bindgen]
impl TautomerScore {
    #[wasm_bindgen(getter)]
    pub fn ring(&self) -> i32 {
        self.inner.ring()
    }
    #[wasm_bindgen(getter)]
    pub fn substructure(&self) -> i32 {
        self.inner.substructure()
    }
    #[wasm_bindgen(getter,js_name=heteroHydrogen)]
    pub fn hetero_hydrogen(&self) -> i32 {
        self.inner.hetero_hydrogen()
    }
    #[wasm_bindgen(getter)]
    pub fn total(&self) -> i32 {
        self.inner.total()
    }
}
type View = Rc<RefCell<ck::TautomerMoleculeView<'static>>>;
#[wasm_bindgen]
pub struct TautomerMoleculeView {
    inner: View,
}
#[wasm_bindgen]
impl TautomerMoleculeView {
    #[wasm_bindgen(js_name=toOwned)]
    pub fn to_owned(&self) -> Self {
        Self {
            inner: Rc::new(RefCell::new(self.inner.borrow().to_owned())),
        }
    }
    #[wasm_bindgen(js_name=numAtoms)]
    pub fn num_atoms(&self) -> usize {
        self.inner.borrow().num_atoms()
    }
    #[wasm_bindgen(js_name=numBonds)]
    pub fn num_bonds(&self) -> usize {
        self.inner.borrow().num_bonds()
    }
    #[wasm_bindgen(js_name=toSmiles)]
    pub fn to_smiles(&self) -> Result<String, JsValue> {
        text(&self.inner.borrow().to_smiles().map_err(run_error)?)
    }
    #[wasm_bindgen(js_name=atomDegree)]
    pub fn atom_degree(&self, id: JsValue) -> Result<Option<usize>, JsValue> {
        Ok(self
            .inner
            .borrow()
            .atom_degree(ck::AtomId::new(usize_value(&id, "atom")?)))
    }
    #[wasm_bindgen(unchecked_return_type = "Atom[]")]
    pub fn atoms(&self) -> JsValue {
        self.inner
            .borrow()
            .atoms()
            .iter()
            .cloned()
            .map(|inner| JsValue::from(crate::atom_bond::Atom { inner }))
            .collect::<Array>()
            .into()
    }
    #[wasm_bindgen(unchecked_return_type = "Bond[]")]
    pub fn bonds(&self) -> JsValue {
        self.inner
            .borrow()
            .bonds()
            .iter()
            .cloned()
            .map(|inner| JsValue::from(crate::atom_bond::Bond { inner }))
            .collect::<Array>()
            .into()
    }
    #[wasm_bindgen(unchecked_return_type = "Atom | null")]
    pub fn atom(&self, id: JsValue) -> Result<JsValue, JsValue> {
        Ok(self
            .inner
            .borrow()
            .atom(ck::AtomId::new(usize_value(&id, "atom")?))
            .cloned()
            .map_or(JsValue::NULL, |inner| {
                crate::atom_bond::Atom { inner }.into()
            }))
    }
    #[wasm_bindgen(unchecked_return_type = "Bond | null")]
    pub fn bond(&self, id: JsValue) -> Result<JsValue, JsValue> {
        Ok(self
            .inner
            .borrow()
            .bond(ck::BondId::new(usize_value(&id, "bond")?))
            .cloned()
            .map_or(JsValue::NULL, |inner| {
                crate::atom_bond::Bond { inner }.into()
            }))
    }
    #[wasm_bindgen(js_name=atomMetadata,unchecked_return_type="AtomMetadata[]")]
    pub fn atom_metadata(&self) -> Result<JsValue, JsValue> {
        self.inner
            .borrow()
            .atom_metadata()
            .map(|v| {
                v.into_iter()
                    .map(|inner| JsValue::from(crate::atom_bond::AtomMetadata { inner }))
                    .collect::<Array>()
                    .into()
            })
            .map_err(|e| crate::valence_errors::valence_error(&e).unwrap_or_else(|e| e))
    }
    pub fn properties(&self) -> crate::property_values::MoleculeProperties {
        crate::property_values::MoleculeProperties {
            inner: self.inner.borrow().properties().clone(),
        }
    }
    #[wasm_bindgen(js_name=tautomerScore)]
    pub fn tautomer_score(&self) -> Result<TautomerScore, JsValue> {
        self.inner
            .borrow_mut()
            .tautomer_score()
            .map(|inner| TautomerScore { inner })
            .map_err(run_error)
    }
}
#[wasm_bindgen]
pub struct TautomerProgress {
    status: ck::TautomerEnumerationStatus,
    num_transforms: u32,
    modified_atoms: Vec<usize>,
    modified_bonds: Vec<usize>,
    entries: Vec<(String, View)>,
}
impl TautomerProgress {
    fn snapshot(progress: &mut ck::TautomerProgress<'_>) -> Result<Self, JsValue> {
        Ok(Self {
            status: progress.status(),
            num_transforms: progress.num_transforms(),
            modified_atoms: progress
                .modified_atoms()
                .iter()
                .map(|id| id.index())
                .collect(),
            modified_bonds: progress
                .modified_bonds()
                .iter()
                .map(|id| id.index())
                .collect(),
            entries: progress
                .entries()
                .map(|(key, view)| Ok((text(key)?, Rc::new(RefCell::new(view.to_owned())))))
                .collect::<Result<_, JsValue>>()?,
        })
    }
}
#[wasm_bindgen]
impl TautomerProgress {
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    #[wasm_bindgen(js_name=toOwned)]
    pub fn to_owned(&self) -> Self {
        Self {
            status: self.status,
            num_transforms: self.num_transforms,
            modified_atoms: self.modified_atoms.clone(),
            modified_bonds: self.modified_bonds.clone(),
            entries: self
                .entries
                .iter()
                .map(|(key, view)| (key.clone(), Rc::new(RefCell::new(view.borrow().to_owned()))))
                .collect(),
        }
    }
    #[wasm_bindgen(js_name=isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    #[wasm_bindgen(getter)]
    pub fn status(&self) -> TautomerEnumerationStatus {
        status(self.status)
    }
    #[wasm_bindgen(getter,js_name=numTransforms)]
    pub fn num_transforms(&self) -> u32 {
        self.num_transforms
    }
    #[wasm_bindgen(getter,js_name=modifiedAtoms,unchecked_return_type="number[]")]
    pub fn modified_atoms(&self) -> JsValue {
        array_ids(self.modified_atoms.iter().copied())
    }
    #[wasm_bindgen(getter,js_name=modifiedBonds,unchecked_return_type="number[]")]
    pub fn modified_bonds(&self) -> JsValue {
        array_ids(self.modified_bonds.iter().copied())
    }
    #[wasm_bindgen(unchecked_return_type = "Array<[string, TautomerMoleculeView]>")]
    pub fn entries(&self) -> JsValue {
        self.entries
            .iter()
            .map(|(key, view)| {
                let row = Array::new();
                row.push(&key.clone().into());
                row.push(
                    &TautomerMoleculeView {
                        inner: Rc::clone(view),
                    }
                    .into(),
                );
                JsValue::from(row)
            })
            .collect::<Array>()
            .into()
    }
}
#[derive(Default)]
struct CallbackFailure(Mutex<Option<JsValue>>);
impl CallbackFailure {
    fn retain(&self, error: JsValue) -> ck::TautomerRunError {
        let mut slot = self.0.lock().expect("callback error mutex");
        if slot.is_none() {
            *slot = Some(error);
        }
        ck::TautomerRunError::Callback("JavaScript callback failed".into())
    }
}
struct JsCallback {
    function: Function,
    failure: Arc<CallbackFailure>,
}
impl ck::TautomerEnumerationCallback for JsCallback {
    fn should_continue(
        &self,
        mut source: ck::TautomerMoleculeView<'_>,
        mut progress: ck::TautomerProgress<'_>,
    ) -> Result<bool, ck::TautomerRunError> {
        let source_snapshot = Rc::new(RefCell::new(source.to_owned()));
        let saved =
            TautomerProgress::snapshot(&mut progress).map_err(|e| self.failure.retain(e))?;
        let entries = saved.entries.clone();
        let result = self
            .function
            .call2(
                &JsValue::UNDEFINED,
                &TautomerMoleculeView {
                    inner: Rc::clone(&source_snapshot),
                }
                .into(),
                &saved.into(),
            )
            .and_then(|v| bool_value(&v, "callback result"))
            .map_err(|e| self.failure.retain(e));
        // Match the Python restricted cache loan: publish only completed score
        // cache transitions into these exact invocation views, never live input.
        source.retain_score_cache_from(&source_snapshot.borrow())?;
        let actual = progress.entries();
        if actual.len() != entries.len() {
            return Err(ck::TautomerRunError::ScoreCacheTargetMismatch);
        }
        for ((key, mut target), (saved_key, snapshot)) in actual.zip(entries) {
            if key.as_bytes() != saved_key.as_bytes() {
                return Err(ck::TautomerRunError::ScoreCacheTargetMismatch);
            }
            target.retain_score_cache_from(&snapshot.borrow())?;
        }
        result
    }
}
struct JsScorer {
    function: Function,
    failure: Arc<CallbackFailure>,
}
impl ck::TautomerScorer for JsScorer {
    fn score(&self, mut source: ck::TautomerMoleculeView<'_>) -> Result<i32, ck::TautomerRunError> {
        let snapshot = Rc::new(RefCell::new(source.to_owned()));
        let result = self
            .function
            .call1(
                &JsValue::UNDEFINED,
                &TautomerMoleculeView {
                    inner: Rc::clone(&snapshot),
                }
                .into(),
            )
            .and_then(|v| i32_value(&v, "scorer result"))
            .map_err(|e| self.failure.retain(e));
        source.retain_score_cache_from(&snapshot.borrow())?;
        result
    }
}
#[wasm_bindgen]
#[derive(Clone)]
pub struct TautomerParams {
    inner: ck::TautomerParams,
    callback: Option<Function>,
    scorer: Option<Function>,
}
#[cosmolkit_wasm::javascript_options("TautomerParamsOptions", clone_snapshot)]
#[wasm_bindgen]
impl TautomerParams {
    #[wasm_bindgen(getter,js_name=maxTautomers)]
    pub fn max_tautomers(&self) -> u32 {
        self.inner.max_tautomers()
    }
    #[wasm_bindgen(setter,js_name=maxTautomers)]
    pub fn set_max_tautomers(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner
            .set_max_tautomers(u32_value(&value, "maxTautomers")?);
        Ok(())
    }
    #[wasm_bindgen(js_name=setMaxTautomers)]
    pub fn set_max_tautomers_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_max_tautomers(value)
    }
    #[wasm_bindgen(js_name=withMaxTautomers)]
    pub fn with_max_tautomers(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<Self, JsValue> {
        let mut out = self.clone();
        out.set_max_tautomers(value)?;
        Ok(out)
    }
    #[wasm_bindgen(getter,js_name=maxTransforms)]
    pub fn max_transforms(&self) -> u32 {
        self.inner.max_transforms()
    }
    #[wasm_bindgen(setter,js_name=maxTransforms)]
    pub fn set_max_transforms(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner
            .set_max_transforms(u32_value(&value, "maxTransforms")?);
        Ok(())
    }
    #[wasm_bindgen(js_name=setMaxTransforms)]
    pub fn set_max_transforms_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_max_transforms(value)
    }
    #[wasm_bindgen(js_name=withMaxTransforms)]
    pub fn with_max_transforms(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<Self, JsValue> {
        let mut out = self.clone();
        out.set_max_transforms(value)?;
        Ok(out)
    }
    #[wasm_bindgen(getter,js_name=removeSp3Stereo)]
    pub fn remove_sp3_stereo(&self) -> bool {
        self.inner.remove_sp3_stereo()
    }
    #[wasm_bindgen(setter,js_name=removeSp3Stereo)]
    pub fn set_remove_sp3_stereo(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner
            .set_remove_sp3_stereo(bool_value(&value, "removeSp3Stereo")?);
        Ok(())
    }
    #[wasm_bindgen(js_name=setRemoveSp3Stereo)]
    pub fn set_remove_sp3_stereo_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_remove_sp3_stereo(value)
    }
    #[wasm_bindgen(js_name=withRemoveSp3Stereo)]
    pub fn with_remove_sp3_stereo(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<Self, JsValue> {
        let mut out = self.clone();
        out.set_remove_sp3_stereo(value)?;
        Ok(out)
    }
    #[wasm_bindgen(getter,js_name=removeBondStereo)]
    pub fn remove_bond_stereo(&self) -> bool {
        self.inner.remove_bond_stereo()
    }
    #[wasm_bindgen(setter,js_name=removeBondStereo)]
    pub fn set_remove_bond_stereo(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner
            .set_remove_bond_stereo(bool_value(&value, "removeBondStereo")?);
        Ok(())
    }
    #[wasm_bindgen(js_name=setRemoveBondStereo)]
    pub fn set_remove_bond_stereo_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_remove_bond_stereo(value)
    }
    #[wasm_bindgen(js_name=withRemoveBondStereo)]
    pub fn with_remove_bond_stereo(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<Self, JsValue> {
        let mut out = self.clone();
        out.set_remove_bond_stereo(value)?;
        Ok(out)
    }
    #[wasm_bindgen(getter,js_name=removeIsotopicHydrogens)]
    pub fn remove_isotopic_hydrogens(&self) -> bool {
        self.inner.remove_isotopic_hydrogens()
    }
    #[wasm_bindgen(setter,js_name=removeIsotopicHydrogens)]
    pub fn set_remove_isotopic_hydrogens(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner
            .set_remove_isotopic_hydrogens(bool_value(&value, "removeIsotopicHydrogens")?);
        Ok(())
    }
    #[wasm_bindgen(js_name=setRemoveIsotopicHydrogens)]
    pub fn set_remove_isotopic_hydrogens_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_remove_isotopic_hydrogens(value)
    }
    #[wasm_bindgen(js_name=withRemoveIsotopicHydrogens)]
    pub fn with_remove_isotopic_hydrogens(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<Self, JsValue> {
        let mut out = self.clone();
        out.set_remove_isotopic_hydrogens(value)?;
        Ok(out)
    }
    #[wasm_bindgen(getter,js_name=reassignStereo)]
    pub fn reassign_stereo(&self) -> bool {
        self.inner.reassign_stereo()
    }
    #[wasm_bindgen(setter,js_name=reassignStereo)]
    pub fn set_reassign_stereo(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner
            .set_reassign_stereo(bool_value(&value, "reassignStereo")?);
        Ok(())
    }
    #[wasm_bindgen(js_name=setReassignStereo)]
    pub fn set_reassign_stereo_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_reassign_stereo(value)
    }
    #[wasm_bindgen(js_name=withReassignStereo)]
    pub fn with_reassign_stereo(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<Self, JsValue> {
        let mut out = self.clone();
        out.set_reassign_stereo(value)?;
        Ok(out)
    }
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_tautomers: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_transforms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_sp3_stereo: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_bond_stereo: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        remove_isotopic_hydrogens: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] reassign_stereo: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "((source: TautomerMoleculeView, progress: TautomerProgress) => boolean) | null"
        )]
        callback: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "((molecule: TautomerMoleculeView) => number) | null"
        )]
        scorer: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "TautomerScoreParams | TautomerScoreParamsOptions | null"
        )]
        score_params: JsValue,
    ) -> Result<Self, JsValue> {
        let mut out = Self {
            inner: ck::TautomerParams::default(),
            callback: checked_callback(callback, "callback")?,
            scorer: checked_callback(scorer, "scorer")?,
        };
        if !max_tautomers.is_undefined() {
            out.inner
                .set_max_tautomers(u32_value(&max_tautomers, "maxTautomers")?);
        }
        if !max_transforms.is_undefined() {
            out.inner
                .set_max_transforms(u32_value(&max_transforms, "maxTransforms")?);
        }
        if !remove_sp3_stereo.is_undefined() {
            out.inner
                .set_remove_sp3_stereo(bool_value(&remove_sp3_stereo, "removeSp3Stereo")?);
        }
        if !remove_bond_stereo.is_undefined() {
            out.inner
                .set_remove_bond_stereo(bool_value(&remove_bond_stereo, "removeBondStereo")?);
        }
        if !remove_isotopic_hydrogens.is_undefined() {
            out.inner.set_remove_isotopic_hydrogens(bool_value(
                &remove_isotopic_hydrogens,
                "removeIsotopicHydrogens",
            )?);
        }
        if !reassign_stereo.is_undefined() {
            out.inner
                .set_reassign_stereo(bool_value(&reassign_stereo, "reassignStereo")?);
        }
        if !score_params.is_null() {
            out.inner.score_params =
                TautomerScoreParams::from_js_configuration(&score_params)?.inner;
        }
        Ok(out)
    }
    pub fn v1() -> Result<Self, JsValue> {
        ck::TautomerParams::v1()
            .map(|inner| Self {
                inner,
                callback: None,
                scorer: None,
            })
            .map_err(catalog_error)
    }
    #[wasm_bindgen(js_name=fromTransformData)]
    pub fn from_transform_data(data: JsValue) -> Result<Self, JsValue> {
        let rows = sequence(&data, "transform data")?
            .iter()
            .map(|v| {
                let row = sequence(&v, "transform row")?;
                if row.length() != 4 {
                    return Err(type_error("transform row must have four strings"));
                }
                row.iter()
                    .map(|v| v.as_string().ok_or_else(|| type_error("transform field")))
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        let borrowed = rows
            .iter()
            .map(|r| (r[0].as_str(), r[1].as_str(), r[2].as_str(), r[3].as_str()))
            .collect::<Vec<_>>();
        ck::TautomerParams::from_transform_data(&borrowed)
            .map(|inner| Self {
                inner,
                callback: None,
                scorer: None,
            })
            .map_err(catalog_error)
    }
    #[wasm_bindgen(js_name=fromTransformFile)]
    pub fn from_transform_file(path: String) -> Result<Self, JsValue> {
        ck::TautomerParams::from_transform_file(path)
            .map(|inner| Self {
                inner,
                callback: None,
                scorer: None,
            })
            .map_err(catalog_error)
    }
    #[wasm_bindgen(getter,js_name=transformCount)]
    pub fn transform_count(&self) -> usize {
        self.inner.transform_count()
    }
    #[wasm_bindgen(getter,js_name=callback,unchecked_return_type="((source: TautomerMoleculeView, progress: TautomerProgress) => boolean) | null")]
    pub fn callback(&self) -> JsValue {
        self.callback.clone().map_or(JsValue::NULL, Into::into)
    }
    #[wasm_bindgen(setter,js_name=callback)]
    pub fn set_callback_property(
        &mut self,
        #[wasm_bindgen(
            unchecked_param_type = "((source: TautomerMoleculeView, progress: TautomerProgress) => boolean) | null"
        )]
        value: JsValue,
    ) -> Result<(), JsValue> {
        self.callback = checked_callback(value, "callback")?;
        Ok(())
    }
    #[wasm_bindgen(js_name=setCallback)]
    pub fn set_callback(
        &mut self,
        #[wasm_bindgen(
            unchecked_param_type = "((source: TautomerMoleculeView, progress: TautomerProgress) => boolean) | null"
        )]
        value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_callback_property(value)
    }
    #[wasm_bindgen(getter,js_name=scorer,unchecked_return_type="((molecule: TautomerMoleculeView) => number) | null")]
    pub fn scorer(&self) -> JsValue {
        self.scorer.clone().map_or(JsValue::NULL, Into::into)
    }
    #[wasm_bindgen(setter,js_name=scorer)]
    pub fn set_scorer_property(
        &mut self,
        #[wasm_bindgen(
            unchecked_param_type = "((molecule: TautomerMoleculeView) => number) | null"
        )]
        value: JsValue,
    ) -> Result<(), JsValue> {
        self.scorer = checked_callback(value, "scorer")?;
        Ok(())
    }
    #[wasm_bindgen(js_name=setScorer)]
    pub fn set_scorer(
        &mut self,
        #[wasm_bindgen(
            unchecked_param_type = "((molecule: TautomerMoleculeView) => number) | null"
        )]
        value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_scorer_property(value)
    }
    #[wasm_bindgen(getter,js_name=scoreParams)]
    pub fn score_params(&self) -> TautomerScoreParams {
        TautomerScoreParams {
            inner: self.inner.score_params.clone(),
        }
    }
    #[wasm_bindgen(setter,js_name=scoreParams)]
    pub fn set_score_params(
        &mut self,
        #[wasm_bindgen(
            unchecked_param_type = "TautomerScoreParams | TautomerScoreParamsOptions | null"
        )]
        value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.score_params = if value.is_null() {
            Default::default()
        } else {
            TautomerScoreParams::from_js_configuration(&value)?.inner
        };
        Ok(())
    }
}
impl TautomerParams {
    fn run<T>(
        &self,
        run: impl FnOnce(&ck::TautomerParams) -> Result<T, ck::OperationError>,
    ) -> Result<T, JsValue> {
        let mut params = self.inner.clone();
        let failure = Arc::new(CallbackFailure::default());
        params.set_callback(self.callback.as_ref().map(|function| {
            Arc::new(JsCallback {
                function: function.clone(),
                failure: Arc::clone(&failure),
            }) as Arc<dyn ck::TautomerEnumerationCallback>
        }));
        params.set_scorer(self.scorer.as_ref().map(|function| {
            Arc::new(JsScorer {
                function: function.clone(),
                failure: Arc::clone(&failure),
            }) as Arc<dyn ck::TautomerScorer>
        }));
        let result = run(&params);
        if let Some(error) = failure.0.lock().expect("callback error mutex").take() {
            return Err(error);
        }
        operation(result)
    }
}
#[wasm_bindgen]
pub struct TautomerEnumeration {
    inner: ck::TautomerEnumeration,
}
#[wasm_bindgen]
impl TautomerEnumeration {
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    #[wasm_bindgen(js_name=isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    #[wasm_bindgen(getter)]
    pub fn status(&self) -> TautomerEnumerationStatus {
        status(self.inner.status())
    }
    #[wasm_bindgen(getter,js_name=modifiedAtoms,unchecked_return_type="number[]")]
    pub fn modified_atoms(&self) -> JsValue {
        array_ids(self.inner.modified_atoms().iter().map(|id| id.index()))
    }
    #[wasm_bindgen(getter,js_name=modifiedBonds,unchecked_return_type="number[]")]
    pub fn modified_bonds(&self) -> JsValue {
        array_ids(self.inner.modified_bonds().iter().map(|id| id.index()))
    }
    #[wasm_bindgen(unchecked_return_type = "Molecule | null")]
    pub fn get(&self, index: JsValue) -> Result<JsValue, JsValue> {
        Ok(self
            .inner
            .get(usize_value(&index, "index")?)
            .cloned()
            .map_or(JsValue::NULL, |v| molecule(v).into()))
    }
    #[wasm_bindgen(js_name=canonicalSmiles,unchecked_return_type="string[]")]
    pub fn canonical_smiles(&self) -> Result<JsValue, JsValue> {
        self.inner
            .canonical_smiles()
            .into_iter()
            .map(|v| text(v).map(JsValue::from))
            .collect::<Result<Array, _>>()
            .map(Into::into)
    }
    #[wasm_bindgen(unchecked_return_type = "Molecule[]")]
    pub fn iter(&self) -> JsValue {
        self.inner
            .iter()
            .cloned()
            .map(|v| JsValue::from(molecule(v)))
            .collect::<Array>()
            .into()
    }
    #[wasm_bindgen(unchecked_return_type = "Array<[string, Molecule]>")]
    pub fn entries(&self) -> Result<JsValue, JsValue> {
        self.inner
            .entries()
            .map(|(key, mol)| {
                let row = Array::new();
                row.push(&text(key)?.into());
                row.push(&molecule(mol.clone()).into());
                Ok(JsValue::from(row))
            })
            .collect::<Result<Array, JsValue>>()
            .map(Into::into)
    }
    #[wasm_bindgen(js_name=canonicalTautomer)]
    pub fn canonical_tautomer(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "TautomerParams | TautomerParamsOptions")]
        params: JsValue,
    ) -> Result<Molecule, JsValue> {
        self.canonical_tautomer_with_params(&TautomerParams::from_js_configuration(&params)?)
    }
    #[wasm_bindgen(js_name=canonicalTautomerWithParams)]
    pub fn canonical_tautomer_with_params(
        &self,
        params: &TautomerParams,
    ) -> Result<Molecule, JsValue> {
        params
            .run(|p| self.inner.canonical_tautomer_with_params(p))
            .map(molecule)
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=enumerateTautomers)]
    pub fn enumerate_tautomers(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "TautomerParams | TautomerParamsOptions")]
        params: JsValue,
    ) -> Result<TautomerEnumeration, JsValue> {
        self.enumerate_tautomers_with_params(&TautomerParams::from_js_configuration(&params)?)
    }
    #[wasm_bindgen(js_name=enumerateTautomersWithParams)]
    pub fn enumerate_tautomers_with_params(
        &self,
        params: &TautomerParams,
    ) -> Result<TautomerEnumeration, JsValue> {
        params
            .run(|p| self.inner.enumerate_tautomers_with_params(p))
            .map(|inner| TautomerEnumeration { inner })
    }
    #[wasm_bindgen(js_name=canonicalTautomer)]
    pub fn canonical_tautomer(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "TautomerParams | TautomerParamsOptions")]
        params: JsValue,
    ) -> Result<Molecule, JsValue> {
        self.canonical_tautomer_with_params(&TautomerParams::from_js_configuration(&params)?)
    }
    #[wasm_bindgen(js_name=canonicalTautomerWithParams)]
    pub fn canonical_tautomer_with_params(
        &self,
        params: &TautomerParams,
    ) -> Result<Molecule, JsValue> {
        params
            .run(|p| self.inner.canonical_tautomer_with_params(p))
            .map(|inner| Molecule {
                inner: inner.into(),
            })
    }
    #[wasm_bindgen(js_name=tautomerScore)]
    pub fn tautomer_score(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "TautomerScoreParams | TautomerScoreParamsOptions"
        )]
        params: JsValue,
    ) -> Result<TautomerScore, JsValue> {
        self.tautomer_score_with_params(&TautomerScoreParams::from_js_configuration(&params)?)
    }
    #[wasm_bindgen(js_name=tautomerScoreWithParams)]
    pub fn tautomer_score_with_params(
        &self,
        params: &TautomerScoreParams,
    ) -> Result<TautomerScore, JsValue> {
        operation(self.inner.tautomer_score_with_params(&params.inner))
            .map(|inner| TautomerScore { inner })
    }
}
