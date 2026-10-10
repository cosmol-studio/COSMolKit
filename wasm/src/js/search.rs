//! Full search policies and detached query/result projection.
use crate::Molecule;
use crate::host_values::*;
use crate::query_values::QueryGraph;
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use std::sync::{Arc, Mutex};
use wasm_bindgen::prelude::*;
/// Read-only query carrier, including wildcard atomic identities.
#[wasm_bindgen]
pub struct QueryAtom {
    inner: ck::QueryAtom,
}
#[wasm_bindgen]
impl QueryAtom {
    #[wasm_bindgen(js_name = id)]
    pub fn id(&self) -> usize {
        self.inner.id().index()
    }
    #[wasm_bindgen(js_name = atomicNumber)]
    pub fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }
    #[wasm_bindgen(js_name = formalCharge)]
    pub fn formal_charge(&self) -> i8 {
        self.inner.formal_charge()
    }
    #[wasm_bindgen(js_name = explicitHydrogens)]
    pub fn explicit_hydrogens(&self) -> u8 {
        self.inner.explicit_hydrogens()
    }
    #[wasm_bindgen(js_name = isotope, unchecked_return_type = "number | null")]
    pub fn isotope(&self) -> JsValue {
        self.inner
            .isotope()
            .map_or(JsValue::NULL, |v| JsValue::from(v as u32))
    }
    #[wasm_bindgen(js_name = atomMap, unchecked_return_type = "number | null")]
    pub fn atom_map(&self) -> JsValue {
        self.inner.atom_map().map_or(JsValue::NULL, JsValue::from)
    }
    #[wasm_bindgen(js_name = isAromatic)]
    pub fn is_aromatic(&self) -> bool {
        self.inner.is_aromatic()
    }
    #[wasm_bindgen(js_name = noImplicit)]
    pub fn no_implicit(&self) -> bool {
        self.inner.no_implicit()
    }
    #[wasm_bindgen(js_name = radicalElectrons)]
    pub fn radical_electrons(&self) -> u8 {
        self.inner.radical_electrons()
    }
}
#[wasm_bindgen]
#[derive(Clone)]
pub struct SubstructMatchParams {
    pub(crate) inner: ck::SubstructMatchParams,
    final_match: Option<js_sys::Function>,
    atom_match: Option<js_sys::Function>,
    bond_match: Option<js_sys::Function>,
}
#[cosmolkit_wasm::javascript_options("SubstructMatchOptions", clone_snapshot)]
#[wasm_bindgen]
impl SubstructMatchParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_matches: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] uniquify: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_chirality: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_enhanced_stereo: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        specified_stereo_query_matches_unspecified: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_query_query_matches: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] recursion_possible: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_recursive_matches: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_threads: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        aromatic_matches_conjugated: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        aromatic_matches_single_or_double: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "string[] | null")] atom_properties: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "string[] | null")] bond_properties: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        extra_atom_check_overrides_default_check: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        extra_bond_check_overrides_default_check: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_generic_matchers: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "((target: Molecule, atoms: number[]) => boolean) | null"
        )]
        final_match: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "((query: QueryAtom, target: Atom) => boolean) | null"
        )]
        atom_match: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "((query: Bond, target: Bond) => boolean) | null"
        )]
        bond_match: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::SubstructMatchParams::default();
        if !max_matches.is_undefined() {
            inner.max_matches = usize_value(&max_matches, "maxMatches")?;
        }
        if !uniquify.is_undefined() {
            inner.uniquify = bool_value(&uniquify, "uniquify")?;
        }
        if !use_chirality.is_undefined() {
            inner.use_chirality = bool_value(&use_chirality, "useChirality")?;
        }
        if !use_enhanced_stereo.is_undefined() {
            inner.use_enhanced_stereo = bool_value(&use_enhanced_stereo, "useEnhancedStereo")?;
        }
        if !specified_stereo_query_matches_unspecified.is_undefined() {
            inner.specified_stereo_query_matches_unspecified = bool_value(
                &specified_stereo_query_matches_unspecified,
                "specifiedStereoQueryMatchesUnspecified",
            )?;
        }
        if !use_query_query_matches.is_undefined() {
            inner.use_query_query_matches =
                bool_value(&use_query_query_matches, "useQueryQueryMatches")?;
        }
        if !recursion_possible.is_undefined() {
            inner.recursion_possible = bool_value(&recursion_possible, "recursionPossible")?;
        }
        if !max_recursive_matches.is_undefined() {
            inner.max_recursive_matches =
                usize_value(&max_recursive_matches, "maxRecursiveMatches")?;
        }
        if !num_threads.is_undefined() {
            inner.num_threads = i32_value(&num_threads, "numThreads")?;
        }
        if !aromatic_matches_conjugated.is_undefined() {
            inner.aromatic_matches_conjugated =
                bool_value(&aromatic_matches_conjugated, "aromaticMatchesConjugated")?;
        }
        if !aromatic_matches_single_or_double.is_undefined() {
            inner.aromatic_matches_single_or_double = bool_value(
                &aromatic_matches_single_or_double,
                "aromaticMatchesSingleOrDouble",
            )?;
        }
        if !atom_properties.is_undefined() {
            inner.atom_properties = strings(&atom_properties, "atomProperties")?;
        }
        if !bond_properties.is_undefined() {
            inner.bond_properties = strings(&bond_properties, "bondProperties")?;
        }
        if !extra_atom_check_overrides_default_check.is_undefined() {
            inner.extra_atom_check_overrides_default_check = bool_value(
                &extra_atom_check_overrides_default_check,
                "extraAtomCheckOverridesDefaultCheck",
            )?;
        }
        if !extra_bond_check_overrides_default_check.is_undefined() {
            inner.extra_bond_check_overrides_default_check = bool_value(
                &extra_bond_check_overrides_default_check,
                "extraBondCheckOverridesDefaultCheck",
            )?;
        }
        if !use_generic_matchers.is_undefined() {
            inner.use_generic_matchers = bool_value(&use_generic_matchers, "useGenericMatchers")?;
        }
        Ok(Self {
            inner,
            final_match: callback(&final_match, "finalMatch")?,
            atom_match: callback(&atom_match, "atomMatch")?,
            bond_match: callback(&bond_match, "bondMatch")?,
        })
    }
    #[wasm_bindgen(getter,js_name=maxMatches)]
    pub fn max_matches(&self) -> usize {
        self.inner.max_matches
    }
    #[wasm_bindgen(getter,js_name=uniquify)]
    pub fn uniquify(&self) -> bool {
        self.inner.uniquify
    }
    #[wasm_bindgen(getter,js_name=useChirality)]
    pub fn use_chirality(&self) -> bool {
        self.inner.use_chirality
    }
    #[wasm_bindgen(getter,js_name=useEnhancedStereo)]
    pub fn use_enhanced_stereo(&self) -> bool {
        self.inner.use_enhanced_stereo
    }
    #[wasm_bindgen(getter,js_name=specifiedStereoQueryMatchesUnspecified)]
    pub fn specified_stereo_query_matches_unspecified(&self) -> bool {
        self.inner.specified_stereo_query_matches_unspecified
    }
    #[wasm_bindgen(getter,js_name=useQueryQueryMatches)]
    pub fn use_query_query_matches(&self) -> bool {
        self.inner.use_query_query_matches
    }
    #[wasm_bindgen(getter,js_name=recursionPossible)]
    pub fn recursion_possible(&self) -> bool {
        self.inner.recursion_possible
    }
    #[wasm_bindgen(getter,js_name=maxRecursiveMatches)]
    pub fn max_recursive_matches(&self) -> usize {
        self.inner.max_recursive_matches
    }
    #[wasm_bindgen(getter,js_name=numThreads)]
    pub fn num_threads(&self) -> i32 {
        self.inner.num_threads
    }
    #[wasm_bindgen(getter,js_name=aromaticMatchesConjugated)]
    pub fn aromatic_matches_conjugated(&self) -> bool {
        self.inner.aromatic_matches_conjugated
    }
    #[wasm_bindgen(getter,js_name=aromaticMatchesSingleOrDouble)]
    pub fn aromatic_matches_single_or_double(&self) -> bool {
        self.inner.aromatic_matches_single_or_double
    }
    #[wasm_bindgen(getter,js_name=atomProperties,unchecked_return_type="string[]")]
    pub fn atom_properties(&self) -> Array {
        self.inner
            .atom_properties
            .iter()
            .map(JsValue::from)
            .collect()
    }
    #[wasm_bindgen(getter,js_name=bondProperties,unchecked_return_type="string[]")]
    pub fn bond_properties(&self) -> Array {
        self.inner
            .bond_properties
            .iter()
            .map(JsValue::from)
            .collect()
    }
    #[wasm_bindgen(getter,js_name=extraAtomCheckOverridesDefaultCheck)]
    pub fn extra_atom_check_overrides_default_check(&self) -> bool {
        self.inner.extra_atom_check_overrides_default_check
    }
    #[wasm_bindgen(getter,js_name=extraBondCheckOverridesDefaultCheck)]
    pub fn extra_bond_check_overrides_default_check(&self) -> bool {
        self.inner.extra_bond_check_overrides_default_check
    }
    #[wasm_bindgen(getter,js_name=useGenericMatchers)]
    pub fn use_generic_matchers(&self) -> bool {
        self.inner.use_generic_matchers
    }
    #[wasm_bindgen(getter, js_name = finalMatch, unchecked_return_type = "((target: Molecule, atoms: number[]) => boolean) | null")]
    pub fn final_match(&self) -> JsValue {
        self.final_match
            .as_ref()
            .map_or(JsValue::NULL, |v| v.clone().into())
    }
    #[wasm_bindgen(getter, js_name = atomMatch, unchecked_return_type = "((query: QueryAtom, target: Atom) => boolean) | null")]
    pub fn atom_match(&self) -> JsValue {
        self.atom_match
            .as_ref()
            .map_or(JsValue::NULL, |v| v.clone().into())
    }
    #[wasm_bindgen(getter, js_name = bondMatch, unchecked_return_type = "((query: Bond, target: Bond) => boolean) | null")]
    pub fn bond_match(&self) -> JsValue {
        self.bond_match
            .as_ref()
            .map_or(JsValue::NULL, |v| v.clone().into())
    }
}
fn callback(value: &JsValue, name: &str) -> Result<Option<js_sys::Function>, JsValue> {
    if value.is_null() || value.is_undefined() {
        return Ok(None);
    }
    value
        .clone()
        .dyn_into()
        .map(Some)
        .map_err(|_| type_error(name))
}

#[derive(Default)]
struct CallbackErrors(Mutex<Option<JsValue>>);
impl CallbackErrors {
    fn call(&self, call: impl FnOnce() -> Result<JsValue, JsValue>) -> bool {
        // Do not call more user code after its first exception. Return false to
        // the canonical matcher, then rethrow that same exception at the boundary.
        if self.0.lock().expect("callback error mutex").is_some() {
            return false;
        }
        let result = call().and_then(|v| {
            v.as_bool()
                .ok_or_else(|| type_error("callback result must be boolean"))
        });
        match result {
            Ok(value) => value,
            Err(error) => {
                *self.0.lock().expect("callback error mutex") = Some(error);
                false
            }
        }
    }
}
fn strings(v: &JsValue, name: &str) -> Result<Vec<String>, JsValue> {
    if v.is_null() {
        return Ok(Vec::new());
    }
    sequence(v, name)?
        .iter()
        .map(|v| v.as_string().ok_or_else(|| type_error(name)))
        .collect()
}
#[wasm_bindgen]
pub struct SmartsWriteParams {
    pub(crate) inner: ck::SmartsWriteParams,
}
#[cosmolkit_wasm::javascript_options("SmartsWriteParamsOptions")]
#[wasm_bindgen]
impl SmartsWriteParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_atom_maps: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] isomeric_smiles: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_dative_bonds: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] rooted_at_atom: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::SmartsWriteParams::default();
        if !include_atom_maps.is_undefined() {
            inner.include_atom_maps = bool_value(&include_atom_maps, "includeAtomMaps")?;
        }
        if !isomeric_smiles.is_undefined() {
            inner.isomeric_smiles = bool_value(&isomeric_smiles, "isomericSmiles")?;
        }
        if !include_dative_bonds.is_undefined() {
            inner.include_dative_bonds = bool_value(&include_dative_bonds, "includeDativeBonds")?;
        }
        if !rooted_at_atom.is_null() && !rooted_at_atom.is_undefined() {
            inner.rooted_at_atom = Some(usize_value(&rooted_at_atom, "rootedAtAtom")?);
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=includeAtomMaps)]
    pub fn include_atom_maps(&self) -> bool {
        self.inner.include_atom_maps
    }
    #[wasm_bindgen(getter,js_name=isomericSmiles)]
    pub fn isomeric_smiles(&self) -> bool {
        self.inner.isomeric_smiles
    }
    #[wasm_bindgen(getter,js_name=includeDativeBonds)]
    pub fn include_dative_bonds(&self) -> bool {
        self.inner.include_dative_bonds
    }
    #[wasm_bindgen(getter,js_name=rootedAtAtom,unchecked_return_type="number | null")]
    pub fn rooted_at_atom(&self) -> JsValue {
        self.inner
            .rooted_at_atom
            .map_or(JsValue::NULL, |v| JsValue::from(v as u32))
    }
}
#[wasm_bindgen]
pub struct MatchResult {
    inner: ck::MatchResult,
}
#[wasm_bindgen]
impl MatchResult {
    #[wasm_bindgen(js_name=atomMapping,unchecked_return_type="number[]")]
    pub fn atom_mapping(&self) -> Array {
        self.inner
            .atom_mapping
            .iter()
            .map(|v| JsValue::from(*v as u32))
            .collect()
    }
    #[wasm_bindgen(js_name=bondMapping,unchecked_return_type="number[]")]
    pub fn bond_mapping(&self) -> Array {
        self.inner
            .bond_mapping
            .iter()
            .map(|v| JsValue::from(*v as u32))
            .collect()
    }
    #[wasm_bindgen(js_name=atomPairs,unchecked_return_type="[number, number][]")]
    pub fn atom_pairs(&self) -> Array {
        self.inner
            .atom_mapping
            .iter()
            .enumerate()
            .map(|(i, v)| {
                let a = Array::new();
                a.push(&(i as u32).into());
                a.push(&(*v as u32).into());
                JsValue::from(a)
            })
            .collect()
    }
}
#[wasm_bindgen]
pub struct CompiledQuery {
    inner: ck::CompiledQuery,
}
#[wasm_bindgen]
impl CompiledQuery {
    #[wasm_bindgen(js_name=numAtoms)]
    pub fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    #[wasm_bindgen(js_name=numBonds)]
    pub fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }
    #[wasm_bindgen(js_name=atomOrder,unchecked_return_type="number[]")]
    pub fn atom_order(&self) -> Array {
        self.inner
            .atom_order()
            .iter()
            .map(|v| JsValue::from(*v as u32))
            .collect()
    }
    pub fn query(&self) -> QueryGraph {
        QueryGraph {
            inner: self.inner.query().clone(),
        }
    }
}
#[wasm_bindgen(js_name=compileQuery)]
pub fn compile_query(
    #[wasm_bindgen(unchecked_param_type = "QueryGraph")] query: JsValue,
) -> Result<CompiledQuery, JsValue> {
    // COSMolKit❗✔️: ck::compile_query(&query.inner)
    let mut result = None;
    visit_querygraph(&query, &mut |query: &QueryGraph| {
        result = Some(
            ck::compile_query(&query.inner)
                .map(|inner| CompiledQuery { inner })
                .map_err(|e| crate::search_errors::compile_error(&e).unwrap_or_else(|e| e)),
        );
    })?;
    result.ok_or_else(|| type_error("query"))?
}
#[wasm_bindgen(js_name=writeSmarts)]
pub fn write_smarts(
    #[wasm_bindgen(unchecked_param_type = "QueryGraph")] query: JsValue,
    #[wasm_bindgen(unchecked_param_type = "SmartsWriteParams")] params: JsValue,
) -> Result<String, JsValue> {
    // COSMolKit❗✔️: ck::write_smarts(&query.inner,&params.inner)
    let mut result = None;
    visit_querygraph(&query, &mut |query: &QueryGraph| {
        let mut inner_result = None;
        let visited = visit_write_params(&params, &mut |params: &SmartsWriteParams| {
            inner_result = Some(
                ck::write_smarts(&query.inner, &params.inner)
                    .map_err(|e| crate::search_errors::write_error(&e).unwrap_or_else(|e| e))
                    .and_then(|v| crate::host_values::text(&v)),
            );
        });
        result = Some(visited.and_then(|()| inner_result.ok_or_else(|| type_error("params"))?));
    })?;
    result.ok_or_else(|| type_error("query"))?
}
#[wasm_bindgen(js_name=writeCxSmarts)]
pub fn write_cx_smarts(
    #[wasm_bindgen(unchecked_param_type = "QueryGraph")] query: JsValue,
    #[wasm_bindgen(unchecked_param_type = "SmartsWriteParams")] params: JsValue,
) -> Result<String, JsValue> {
    // COSMolKit❗✔️: ck::write_cx_smarts(&query.inner,&params.inner)
    let mut result = None;
    visit_querygraph(&query, &mut |query: &QueryGraph| {
        let mut inner_result = None;
        let visited = visit_write_params(&params, &mut |params: &SmartsWriteParams| {
            inner_result = Some(
                ck::write_cx_smarts(&query.inner, &params.inner)
                    .map_err(|e| crate::search_errors::write_error(&e).unwrap_or_else(|e| e))
                    .and_then(|v| crate::host_values::text(&v)),
            );
        });
        result = Some(visited.and_then(|()| inner_result.ok_or_else(|| type_error("params"))?));
    })?;
    result.ok_or_else(|| type_error("query"))?
}
#[wasm_bindgen(
    inline_js = "export function visitSearchSmartsWriteParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid SmartsWriteParams',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitSearchSmartsWriteParams)]
    fn visit_write_params(
        v: &JsValue,
        f: &mut dyn FnMut(&SmartsWriteParams),
    ) -> Result<(), JsValue>;
}
#[wasm_bindgen(
    inline_js = "export function visitSearchSubstructMatchParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid SubstructMatchParams',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitSearchSubstructMatchParams)]
    fn visit_substructmatchparams(
        v: &JsValue,
        f: &mut dyn FnMut(&SubstructMatchParams),
    ) -> Result<(), JsValue>;
}
impl SubstructMatchParams {
    pub(crate) fn from_inner(inner: ck::SubstructMatchParams) -> Self {
        Self {
            inner,
            final_match: None,
            atom_match: None,
            bond_match: None,
        }
    }
    pub(crate) fn has_callbacks(&self) -> bool {
        self.final_match.is_some() || self.atom_match.is_some() || self.bond_match.is_some()
    }
    fn with_callbacks<T>(
        &self,
        molecule: &cosmolkit_wasm::Molecule,
        run: impl FnOnce(&ck::SubstructMatchParams) -> Result<T, ck::SubstructMatchError>,
    ) -> Result<T, JsValue> {
        let mut params = self.inner.clone();
        let errors = Arc::new(CallbackErrors::default());
        if let Some(callback) = &self.final_match {
            let callback = callback.clone();
            let errors = Arc::clone(&errors);
            // Match Python's owned COW callback snapshot; never give a JS
            // callback a mutable reference to the receiver's live runtime.
            let snapshot = Mutex::new(molecule.clone());
            params.extra_final_check = Some(Arc::new(move |_, ids| {
                errors.call(|| {
                    let target: JsValue = Molecule {
                        inner: Arc::new(snapshot.lock().expect("callback snapshot mutex").clone()),
                    }
                    .into();
                    let ids: Array = ids.iter().map(|id| JsValue::from(*id as u32)).collect();
                    callback.call2(&JsValue::UNDEFINED, &target, &ids)
                })
            }));
        }
        if let Some(callback) = &self.atom_match {
            let callback = callback.clone();
            let errors = Arc::clone(&errors);
            params.extra_atom_check = Some(Arc::new(move |_, query, _, target| {
                errors.call(|| {
                    let query: JsValue = QueryAtom {
                        inner: query.clone(),
                    }
                    .into();
                    let target: JsValue = crate::atom_bond::Atom {
                        inner: target.clone(),
                    }
                    .into();
                    callback.call2(&JsValue::UNDEFINED, &query, &target)
                })
            }));
        }
        if let Some(callback) = &self.bond_match {
            let callback = callback.clone();
            let errors = Arc::clone(&errors);
            params.extra_bond_check = Some(Arc::new(move |query, target| {
                errors.call(|| {
                    let query: JsValue = crate::atom_bond::Bond {
                        inner: query.clone(),
                    }
                    .into();
                    let target: JsValue = crate::atom_bond::Bond {
                        inner: target.clone(),
                    }
                    .into();
                    callback.call2(&JsValue::UNDEFINED, &query, &target)
                })
            }));
        }
        let result = run(&params);
        let error = errors.0.lock().expect("callback error mutex").take();
        match error {
            Some(error) => Err(error),
            None => {
                result.map_err(|e| crate::search_errors::substruct_error(&e).unwrap_or_else(|e| e))
            }
        }
    }
}
#[wasm_bindgen(
    inline_js = "export function visitSearchQueryGraph(v,f){try{f(v);}catch(cause){throw new TypeError('invalid QueryGraph',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitSearchQueryGraph)]
    fn visit_querygraph(v: &JsValue, f: &mut dyn FnMut(&QueryGraph)) -> Result<(), JsValue>;
}
#[wasm_bindgen(
    inline_js = "export function visitSearchCompiledQuery(v,f){try{f(v);}catch(cause){throw new TypeError('invalid CompiledQuery',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitSearchCompiledQuery)]
    fn visit_compiledquery(v: &JsValue, f: &mut dyn FnMut(&CompiledQuery)) -> Result<(), JsValue>;
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=toSmarts)]
    pub fn to_smarts(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "SmartsWriteParams | SmartsWriteParamsOptions"
        )]
        params: JsValue,
    ) -> Result<String, JsValue> {
        self.to_smarts_with_params(&SmartsWriteParams::from_js_configuration(&params)?)
    }
    #[wasm_bindgen(js_name=toSmartsWithParams)]
    pub fn to_smarts_with_params(&self, params: &SmartsWriteParams) -> Result<String, JsValue> {
        self.inner
            .to_smarts_with_params(&params.inner)
            .map_err(|e| crate::search_errors::write_error(&e).unwrap_or_else(|e| e))
            .and_then(|v| text(&v))
    }
    #[wasm_bindgen(js_name=toCxSmarts)]
    pub fn to_cx_smarts(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "SmartsWriteParams | SmartsWriteParamsOptions"
        )]
        params: JsValue,
    ) -> Result<String, JsValue> {
        self.to_cx_smarts_with_params(&SmartsWriteParams::from_js_configuration(&params)?)
    }
    #[wasm_bindgen(js_name=toCxSmartsWithParams)]
    pub fn to_cx_smarts_with_params(&self, params: &SmartsWriteParams) -> Result<String, JsValue> {
        self.inner
            .to_cx_smarts_with_params(&params.inner)
            .map_err(|e| crate::search_errors::write_error(&e).unwrap_or_else(|e| e))
            .and_then(|v| text(&v))
    }
    #[wasm_bindgen(js_name=substructMatch,unchecked_return_type="MatchResult | null")]
    pub fn substruct_match(
        &self,
        #[wasm_bindgen(unchecked_param_type = "QueryGraph")] query: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "SubstructMatchParams | SubstructMatchOptions"
        )]
        params: JsValue,
    ) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.inner.substruct_match(&q.inner)
        let mut result = None;
        let configured = SubstructMatchParams::from_js_configuration(&params)?;
        visit_querygraph(&query, &mut |q: &QueryGraph| {
            result = Some(
                configured
                    .with_callbacks(&self.inner, |params| {
                        self.inner.substruct_match_with_params(&q.inner, params)
                    })
                    .map(|v| v.map_or(JsValue::NULL, |inner| MatchResult { inner }.into())),
            );
        })?;
        result.ok_or_else(|| type_error("query"))?
    }
    #[wasm_bindgen(js_name=substructMatches)]
    pub fn substruct_matches(
        &self,
        #[wasm_bindgen(unchecked_param_type = "QueryGraph")] query: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "SubstructMatchParams | SubstructMatchOptions"
        )]
        params: JsValue,
    ) -> Result<Vec<MatchResult>, JsValue> {
        // COSMolKit❗✔️: self.inner.substruct_matches(&q.inner)
        let mut result = None;
        let configured = SubstructMatchParams::from_js_configuration(&params)?;
        visit_querygraph(&query, &mut |q: &QueryGraph| {
            result = Some(
                configured
                    .with_callbacks(&self.inner, |params| {
                        self.inner.substruct_matches_with_params(&q.inner, params)
                    })
                    .map(|v| v.into_iter().map(|inner| MatchResult { inner }).collect()),
            );
        })?;
        result.ok_or_else(|| type_error("query"))?
    }
    #[wasm_bindgen(js_name=hasSubstructMatch)]
    pub fn has_substruct_match(
        &self,
        #[wasm_bindgen(unchecked_param_type = "QueryGraph")] query: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "SubstructMatchParams | SubstructMatchOptions"
        )]
        params: JsValue,
    ) -> Result<bool, JsValue> {
        // COSMolKit❗✔️: self.inner.has_substruct_match(&q.inner)
        let mut result = None;
        let configured = SubstructMatchParams::from_js_configuration(&params)?;
        visit_querygraph(&query, &mut |q: &QueryGraph| {
            result = Some(configured.with_callbacks(&self.inner, |params| {
                self.inner.has_substruct_match_with_params(&q.inner, params)
            }));
        })?;
        result.ok_or_else(|| type_error("query"))?
    }
    #[wasm_bindgen(js_name=substructMatchesWithParams)]
    pub fn substruct_matches_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "QueryGraph")] query: JsValue,
        #[wasm_bindgen(unchecked_param_type = "SubstructMatchParams")] params: JsValue,
    ) -> Result<Vec<MatchResult>, JsValue> {
        // COSMolKit❗✔️: self.inner.substruct_matches_with_params(&q.inner,&p.inner)
        let mut result = None;
        visit_querygraph(&query, &mut |q: &QueryGraph| {
            let mut inner_result = None;
            let visited = visit_substructmatchparams(&params, &mut |p: &SubstructMatchParams| {
                inner_result = Some(
                    p.with_callbacks(&self.inner, |params| {
                        self.inner.substruct_matches_with_params(&q.inner, params)
                    })
                    .map(|v| v.into_iter().map(|inner| MatchResult { inner }).collect()),
                );
            });
            result = Some(visited.and_then(|()| inner_result.ok_or_else(|| type_error("params"))?));
        })?;
        result.ok_or_else(|| type_error("query"))?
    }
    #[wasm_bindgen(js_name=substructMatchWithParams,unchecked_return_type="MatchResult | null")]
    pub fn substruct_match_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "QueryGraph")] query: JsValue,
        #[wasm_bindgen(unchecked_param_type = "SubstructMatchParams")] params: JsValue,
    ) -> Result<JsValue, JsValue> {
        let mut result = None;
        visit_querygraph(&query, &mut |q: &QueryGraph| {
            let mut inner_result = None;
            let visited = visit_substructmatchparams(&params, &mut |p: &SubstructMatchParams| {
                inner_result = Some(
                    p.with_callbacks(&self.inner, |params| {
                        self.inner.substruct_match_with_params(&q.inner, params)
                    })
                    .map(|v| v.map_or(JsValue::NULL, |inner| MatchResult { inner }.into())),
                );
            });
            result = Some(visited.and_then(|()| inner_result.ok_or_else(|| type_error("params"))?));
        })?;
        result.ok_or_else(|| type_error("query"))?
    }
    #[wasm_bindgen(js_name=hasSubstructMatchWithParams)]
    pub fn has_substruct_match_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "QueryGraph")] query: JsValue,
        #[wasm_bindgen(unchecked_param_type = "SubstructMatchParams")] params: JsValue,
    ) -> Result<bool, JsValue> {
        let mut result = None;
        visit_querygraph(&query, &mut |q: &QueryGraph| {
            let mut inner_result = None;
            let visited = visit_substructmatchparams(&params, &mut |p: &SubstructMatchParams| {
                inner_result = Some(p.with_callbacks(&self.inner, |params| {
                    self.inner.has_substruct_match_with_params(&q.inner, params)
                }));
            });
            result = Some(visited.and_then(|()| inner_result.ok_or_else(|| type_error("params"))?));
        })?;
        result.ok_or_else(|| type_error("query"))?
    }
    #[wasm_bindgen(js_name=substructMatchesCompiled)]
    pub fn substruct_matches_compiled(
        &self,
        #[wasm_bindgen(unchecked_param_type = "CompiledQuery")] query: JsValue,
    ) -> Result<Vec<MatchResult>, JsValue> {
        // COSMolKit❗✔️: self.inner.substruct_matches_compiled(&q.inner)
        let mut result = None;
        visit_compiledquery(&query, &mut |q: &CompiledQuery| {
            result = Some(
                self.inner
                    .substruct_matches_compiled(&q.inner)
                    .map(|v| v.into_iter().map(|inner| MatchResult { inner }).collect())
                    .map_err(|e| crate::search_errors::match_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("query"))?
    }
}
