//! Full search policies and detached query/result projection.
use crate::Molecule;
use crate::host_values::*;
use crate::query_values::QueryGraph;
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct SubstructMatchParams {
    pub(crate) inner: ck::SubstructMatchParams,
}
#[cosmolkit_wasm::javascript_options("SubstructMatchOptions")]
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
        Ok(Self { inner })
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
    inner: ck::SmartsWriteParams,
}
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
    fn from_configuration(value: &JsValue) -> Result<Self, JsValue> {
        let mut inner = None;
        // The generated visitor checks the actual native parameter class.
        // Its callback performs only a detached parameter copy, no operation.
        if visit_substructmatchparams(value, &mut |params: &SubstructMatchParams| {
            inner = Some(params.inner.clone());
        }).is_ok() {
            return inner.map(|inner| Self { inner }).ok_or_else(|| type_error("params"));
        }
        Self::from_js_options(value)
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
    #[wasm_bindgen(js_name=substructMatch,unchecked_return_type="MatchResult | null")]
    pub fn substruct_match(
        &self,
        #[wasm_bindgen(unchecked_param_type = "QueryGraph")] query: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "SubstructMatchParams | SubstructMatchOptions")] params: JsValue,
    ) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.inner.substruct_match(&q.inner)
        let mut result = None;
        let configured = (!params.is_undefined()).then(|| SubstructMatchParams::from_configuration(&params)).transpose()?;
        visit_querygraph(&query, &mut |q: &QueryGraph| {
            result = Some(
                configured.as_ref().map_or_else(
                    || self.inner.substruct_match(&q.inner),
                    |params| self.inner.substruct_match_with_params(&q.inner, &params.inner),
                )
                    .map(|v| v.map_or(JsValue::NULL, |inner| MatchResult { inner }.into()))
                    .map_err(|e| crate::search_errors::substruct_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("query"))?
    }
    #[wasm_bindgen(js_name=substructMatches)]
    pub fn substruct_matches(
        &self,
        #[wasm_bindgen(unchecked_param_type = "QueryGraph")] query: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "SubstructMatchParams | SubstructMatchOptions")] params: JsValue,
    ) -> Result<Vec<MatchResult>, JsValue> {
        // COSMolKit❗✔️: self.inner.substruct_matches(&q.inner)
        let mut result = None;
        let configured = (!params.is_undefined()).then(|| SubstructMatchParams::from_configuration(&params)).transpose()?;
        visit_querygraph(&query, &mut |q: &QueryGraph| {
            result = Some(
                configured.as_ref().map_or_else(
                    || self.inner.substruct_matches(&q.inner),
                    |params| self.inner.substruct_matches_with_params(&q.inner, &params.inner),
                )
                    .map(|v| v.into_iter().map(|inner| MatchResult { inner }).collect())
                    .map_err(|e| crate::search_errors::substruct_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("query"))?
    }
    #[wasm_bindgen(js_name=hasSubstructMatch)]
    pub fn has_substruct_match(
        &self,
        #[wasm_bindgen(unchecked_param_type = "QueryGraph")] query: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "SubstructMatchParams | SubstructMatchOptions")] params: JsValue,
    ) -> Result<bool, JsValue> {
        // COSMolKit❗✔️: self.inner.has_substruct_match(&q.inner)
        let mut result = None;
        let configured = (!params.is_undefined()).then(|| SubstructMatchParams::from_configuration(&params)).transpose()?;
        visit_querygraph(&query, &mut |q: &QueryGraph| {
            result = Some(
                configured.as_ref().map_or_else(
                    || self.inner.has_substruct_match(&q.inner),
                    |params| self.inner.has_substruct_match_with_params(&q.inner, &params.inner),
                )
                    .map_err(|e| crate::search_errors::substruct_error(&e).unwrap_or_else(|e| e)),
            );
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
                    self.inner
                        .substruct_matches_with_params(&q.inner, &p.inner)
                        .map(|v| v.into_iter().map(|inner| MatchResult { inner }).collect())
                        .map_err(|e| {
                            crate::search_errors::substruct_error(&e).unwrap_or_else(|e| e)
                        }),
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
                    self.inner
                        .substruct_match_with_params(&q.inner, &p.inner)
                        .map(|v| v.map_or(JsValue::NULL, |inner| MatchResult { inner }.into()))
                        .map_err(|e| crate::search_errors::substruct_error(&e).unwrap_or_else(|e| e)),
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
                inner_result = Some(
                    self.inner
                        .has_substruct_match_with_params(&q.inner, &p.inner)
                        .map_err(|e| crate::search_errors::substruct_error(&e).unwrap_or_else(|e| e)),
                );
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
