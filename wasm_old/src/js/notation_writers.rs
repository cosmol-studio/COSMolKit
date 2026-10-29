//! Typed notation options and canonical read-only writer calls.
use crate::Molecule;
use crate::host_values::*;
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
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
pub struct CxCoordinateSelection {
    inner: ck::CxCoordinateSelection,
}
#[wasm_bindgen]
impl CxCoordinateSelection {
    pub fn auto() -> Self {
        Self {
            inner: ck::CxCoordinateSelection::Auto,
        }
    }
    #[wasm_bindgen(js_name=twoD)]
    pub fn two_d(
        #[wasm_bindgen(unchecked_param_type = "number")] id: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::CxCoordinateSelection::TwoD {
                id: usize_value(&id, "id")?,
            },
        })
    }
    #[wasm_bindgen(js_name=threeD)]
    pub fn three_d(
        #[wasm_bindgen(unchecked_param_type = "number")] id: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::CxCoordinateSelection::ThreeD {
                id: usize_value(&id, "id")?,
            },
        })
    }
}

#[wasm_bindgen(
    inline_js = "export function visitNotationFields(v,f){try{f(v);}catch(cause){throw new TypeError('invalid CxSmilesFields',{cause});}} export function visitNotationSelection(v,f){try{f(v);}catch(cause){throw new TypeError('invalid CxCoordinateSelection',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitNotationFields)]
    fn visit_fields(value: &JsValue, f: &mut dyn FnMut(&CxSmilesFields)) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitNotationSelection)]
    fn visit_selection(
        value: &JsValue,
        f: &mut dyn FnMut(&CxCoordinateSelection),
    ) -> Result<(), JsValue>;
}

#[wasm_bindgen]
pub struct CxSmilesWriteParams {
    inner: ck::CxSmilesWriteParams,
}
#[cosmolkit_wasm::javascript_options("CxSmilesWriteParamsOptions")]
#[wasm_bindgen]
impl CxSmilesWriteParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(
            unchecked_optional_param_type = "SmilesWriteParams | SmilesWriteParamsOptions"
        )]
        smiles: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "CxSmilesFields")] fields: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "CxCoordinateSelection")]
        coordinate_selection: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::CxSmilesWriteParams::default();
        if !smiles.is_undefined() {
            inner.smiles =
                crate::smiles_parameters::SmilesWriteParams::from_js_configuration(&smiles)?.inner;
        }
        if !fields.is_undefined() {
            visit_fields(&fields, &mut |v| inner.fields = v.inner)?;
        }
        if !coordinate_selection.is_undefined() {
            visit_selection(&coordinate_selection, &mut |v| {
                inner.coordinate_selection = v.inner
            })?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter)]
    pub fn smiles(&self) -> crate::smiles_parameters::SmilesWriteParams {
        crate::smiles_parameters::SmilesWriteParams {
            inner: self.inner.smiles,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn fields(&self) -> CxSmilesFields {
        CxSmilesFields {
            inner: self.inner.fields,
        }
    }
    #[wasm_bindgen(getter,js_name=coordinateSelection)]
    pub fn coordinate_selection(&self) -> CxCoordinateSelection {
        CxCoordinateSelection {
            inner: self.inner.coordinate_selection,
        }
    }
}

#[wasm_bindgen]
pub struct RandomSmilesWriteParams {
    inner: ck::RandomSmilesWriteParams,
}
#[cosmolkit_wasm::javascript_options("RandomSmilesWriteParamsOptions")]
#[wasm_bindgen]
impl RandomSmilesWriteParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] isomeric_smiles: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] kekule: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] all_bonds_explicit: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] all_hydrogens_explicit: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::RandomSmilesWriteParams::default();
        if !isomeric_smiles.is_undefined() {
            inner.isomeric_smiles = bool_value(&isomeric_smiles, "isomericSmiles")?;
        }
        if !kekule.is_undefined() {
            inner.kekule = bool_value(&kekule, "kekule")?;
        }
        if !all_bonds_explicit.is_undefined() {
            inner.all_bonds_explicit = bool_value(&all_bonds_explicit, "allBondsExplicit")?;
        }
        if !all_hydrogens_explicit.is_undefined() {
            inner.all_hydrogens_explicit =
                bool_value(&all_hydrogens_explicit, "allHydrogensExplicit")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=isomericSmiles)]
    pub fn isomeric_smiles(&self) -> bool {
        self.inner.isomeric_smiles
    }
    #[wasm_bindgen(getter)]
    pub fn kekule(&self) -> bool {
        self.inner.kekule
    }
    #[wasm_bindgen(getter,js_name=allBondsExplicit)]
    pub fn all_bonds_explicit(&self) -> bool {
        self.inner.all_bonds_explicit
    }
    #[wasm_bindgen(getter,js_name=allHydrogensExplicit)]
    pub fn all_hydrogens_explicit(&self) -> bool {
        self.inner.all_hydrogens_explicit
    }
}

fn indices(value: &JsValue, name: &str) -> Result<Vec<usize>, JsValue> {
    sequence(value, name)?
        .iter()
        .map(|v| usize_value(&v, name))
        .collect()
}
fn optional_indices(value: &JsValue, name: &str) -> Result<Option<Vec<usize>>, JsValue> {
    if value.is_null() || value.is_undefined() {
        Ok(None)
    } else {
        indices(value, name).map(Some)
    }
}
fn optional_symbols(value: &JsValue, name: &str) -> Result<Option<Vec<String>>, JsValue> {
    if value.is_null() || value.is_undefined() {
        return Ok(None);
    }
    sequence(value, name)?
        .iter()
        .map(|v| v.as_string().ok_or_else(|| type_error(name)))
        .collect::<Result<_, _>>()
        .map(Some)
}
fn indices_js(value: impl Iterator<Item = usize>) -> JsValue {
    value
        .map(|v| JsValue::from(v as u32))
        .collect::<Array>()
        .into()
}
fn bonds_js(value: &Option<Vec<ck::BondId>>) -> JsValue {
    value
        .as_ref()
        .map_or(JsValue::NULL, |v| indices_js(v.iter().map(|id| id.index())))
}
fn symbols_js(value: &Option<Vec<String>>) -> JsValue {
    value.as_ref().map_or(JsValue::NULL, |v| {
        v.iter().map(JsValue::from).collect::<Array>().into()
    })
}

#[wasm_bindgen]
pub struct FragmentSmilesWriteParams {
    inner: ck::FragmentSmilesWriteParams,
}
#[cosmolkit_wasm::javascript_options("FragmentSmilesWriteParamsOptions")]
#[wasm_bindgen]
impl FragmentSmilesWriteParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number[]")] atoms: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "SmilesWriteParams | SmilesWriteParamsOptions"
        )]
        smiles: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | null")] bonds: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "string[] | null")] atom_symbols: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "string[] | null")] bond_symbols: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::FragmentSmilesWriteParams {
                atoms: if atoms.is_undefined() {
                    Vec::new()
                } else {
                    indices(&atoms, "atoms")?
                        .into_iter()
                        .map(ck::AtomId::new)
                        .collect()
                },
                smiles: crate::smiles_parameters::SmilesWriteParams::from_js_configuration(
                    &smiles,
                )?
                .inner,
                bonds: optional_indices(&bonds, "bonds")?
                    .map(|v| v.into_iter().map(ck::BondId::new).collect()),
                atom_symbols: optional_symbols(&atom_symbols, "atomSymbols")?,
                bond_symbols: optional_symbols(&bond_symbols, "bondSymbols")?,
            },
        })
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number[]")]
    pub fn atoms(&self) -> JsValue {
        indices_js(self.inner.atoms.iter().map(|id| id.index()))
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number[] | null")]
    pub fn bonds(&self) -> JsValue {
        bonds_js(&self.inner.bonds)
    }
    #[wasm_bindgen(getter,js_name=atomSymbols,unchecked_return_type="string[] | null")]
    pub fn atom_symbols(&self) -> JsValue {
        symbols_js(&self.inner.atom_symbols)
    }
    #[wasm_bindgen(getter,js_name=bondSymbols,unchecked_return_type="string[] | null")]
    pub fn bond_symbols(&self) -> JsValue {
        symbols_js(&self.inner.bond_symbols)
    }
    #[wasm_bindgen(getter)]
    pub fn smiles(&self) -> crate::smiles_parameters::SmilesWriteParams {
        crate::smiles_parameters::SmilesWriteParams {
            inner: self.inner.smiles,
        }
    }
}

#[wasm_bindgen]
pub struct FragmentCxSmilesWriteParams {
    inner: ck::FragmentCxSmilesWriteParams,
}
#[cosmolkit_wasm::javascript_options("FragmentCxSmilesWriteParamsOptions")]
#[wasm_bindgen]
impl FragmentCxSmilesWriteParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number[]")] atoms: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "CxSmilesWriteParams | CxSmilesWriteParamsOptions"
        )]
        cx: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | null")] bonds: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "string[] | null")] atom_symbols: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "string[] | null")] bond_symbols: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::FragmentCxSmilesWriteParams {
                atoms: if atoms.is_undefined() {
                    Vec::new()
                } else {
                    indices(&atoms, "atoms")?
                        .into_iter()
                        .map(ck::AtomId::new)
                        .collect()
                },
                cx: CxSmilesWriteParams::from_js_configuration(&cx)?.inner,
                bonds: optional_indices(&bonds, "bonds")?
                    .map(|v| v.into_iter().map(ck::BondId::new).collect()),
                atom_symbols: optional_symbols(&atom_symbols, "atomSymbols")?,
                bond_symbols: optional_symbols(&bond_symbols, "bondSymbols")?,
            },
        })
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number[]")]
    pub fn atoms(&self) -> JsValue {
        indices_js(self.inner.atoms.iter().map(|id| id.index()))
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number[] | null")]
    pub fn bonds(&self) -> JsValue {
        bonds_js(&self.inner.bonds)
    }
    #[wasm_bindgen(getter,js_name=atomSymbols,unchecked_return_type="string[] | null")]
    pub fn atom_symbols(&self) -> JsValue {
        symbols_js(&self.inner.atom_symbols)
    }
    #[wasm_bindgen(getter,js_name=bondSymbols,unchecked_return_type="string[] | null")]
    pub fn bond_symbols(&self) -> JsValue {
        symbols_js(&self.inner.bond_symbols)
    }
    #[wasm_bindgen(getter)]
    pub fn cx(&self) -> CxSmilesWriteParams {
        CxSmilesWriteParams {
            inner: self.inner.cx,
        }
    }
}

fn writer_result(value: Result<ck::PropertyText, ck::SmilesWriteError>) -> Result<String, JsValue> {
    value
        .map_err(|e| crate::smiles_boundary::smiles_write_error(&e).unwrap_or_else(|e| e))
        .and_then(|v| text(&v))
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=toCxSmiles)]
    pub fn to_cx_smiles(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "CxSmilesWriteParams | CxSmilesWriteParamsOptions"
        )]
        params: JsValue,
    ) -> Result<String, JsValue> {
        self.to_cx_smiles_with_params(&CxSmilesWriteParams::from_js_configuration(&params)?)
    }
    #[wasm_bindgen(js_name=toCxSmilesWithParams)]
    pub fn to_cx_smiles_with_params(
        &self,
        params: &CxSmilesWriteParams,
    ) -> Result<String, JsValue> {
        writer_result(self.inner.to_cx_smiles_with_params(&params.inner))
    }
    #[wasm_bindgen(js_name=toFragmentSmiles)]
    pub fn to_fragment_smiles(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[]")] atoms: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "FragmentSmilesWriteParams | FragmentSmilesWriteParamsOptions"
        )]
        params: JsValue,
    ) -> Result<String, JsValue> {
        let mut params = FragmentSmilesWriteParams::from_js_configuration(&params)?;
        params.inner.atoms = indices(&atoms, "atoms")?
            .into_iter()
            .map(ck::AtomId::new)
            .collect();
        self.to_fragment_smiles_with_params(&params)
    }
    #[wasm_bindgen(js_name=toFragmentSmilesWithParams)]
    pub fn to_fragment_smiles_with_params(
        &self,
        params: &FragmentSmilesWriteParams,
    ) -> Result<String, JsValue> {
        writer_result(self.inner.to_fragment_smiles_with_params(&params.inner))
    }
    #[wasm_bindgen(js_name=toFragmentCxSmiles)]
    pub fn to_fragment_cx_smiles(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[]")] atoms: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "FragmentCxSmilesWriteParams | FragmentCxSmilesWriteParamsOptions"
        )]
        params: JsValue,
    ) -> Result<String, JsValue> {
        let mut params = FragmentCxSmilesWriteParams::from_js_configuration(&params)?;
        params.inner.atoms = indices(&atoms, "atoms")?
            .into_iter()
            .map(ck::AtomId::new)
            .collect();
        self.to_fragment_cx_smiles_with_params(&params)
    }
    #[wasm_bindgen(js_name=toFragmentCxSmilesWithParams)]
    pub fn to_fragment_cx_smiles_with_params(
        &self,
        params: &FragmentCxSmilesWriteParams,
    ) -> Result<String, JsValue> {
        writer_result(self.inner.to_fragment_cx_smiles_with_params(&params.inner))
    }
    #[wasm_bindgen(js_name=toRandomSmiles,unchecked_return_type="string[]")]
    pub fn to_random_smiles(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] count: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] seed: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "RandomSmilesWriteParams | RandomSmilesWriteParamsOptions"
        )]
        params: JsValue,
    ) -> Result<JsValue, JsValue> {
        self.to_random_smiles_with_params(
            count,
            seed,
            &RandomSmilesWriteParams::from_js_configuration(&params)?,
        )
    }
    #[wasm_bindgen(js_name=toRandomSmilesWithParams,unchecked_return_type="string[]")]
    pub fn to_random_smiles_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] count: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] seed: JsValue,
        params: &RandomSmilesWriteParams,
    ) -> Result<JsValue, JsValue> {
        let values = self
            .inner
            .to_random_smiles_with_params(
                u32_value(&count, "count")?,
                u32_value(&seed, "seed")?,
                &params.inner,
            )
            .map_err(|e| crate::smiles_boundary::smiles_write_error(&e).unwrap_or_else(|e| e))?;
        values
            .iter()
            .map(text)
            .map(|v| v.map(JsValue::from))
            .collect::<Result<Array, _>>()
            .map(Into::into)
    }
}
