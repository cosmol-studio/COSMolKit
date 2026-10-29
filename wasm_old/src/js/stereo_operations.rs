//! Typed stereo options and snapshots; all perception remains in cosmolkit.
use crate::{Molecule, host_values::*};
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use wasm_bindgen::prelude::*;
fn operation<T>(result: Result<T, ck::OperationError>) -> Result<T, JsValue> {
    result.map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
}
fn indices(value: &JsValue, name: &str) -> Result<Option<Vec<usize>>, JsValue> {
    if value.is_null() || value.is_undefined() {
        Ok(None)
    } else {
        sequence(value, name)?
            .iter()
            .map(|v| usize_value(&v, name))
            .collect::<Result<_, _>>()
            .map(Some)
    }
}
fn ids(value: impl Iterator<Item = usize>) -> JsValue {
    value
        .map(|v| JsValue::from(v as u32))
        .collect::<Array>()
        .into()
}

#[wasm_bindgen]
pub struct StructureTagParams {
    inner: ck::StructureTagParams,
}
#[cosmolkit_wasm::javascript_options("StructureTagParamsOptions")]
#[wasm_bindgen]
impl StructureTagParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] replace_existing_tags: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::StructureTagParams::default();
        if !conformer_id.is_undefined() {
            inner.conformer_id = i32_value(&conformer_id, "conformerId")?;
        }
        if !replace_existing_tags.is_undefined() {
            inner.replace_existing_tags =
                bool_value(&replace_existing_tags, "replaceExistingTags")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=conformerId)]
    pub fn conformer_id(&self) -> i32 {
        self.inner.conformer_id
    }
    #[wasm_bindgen(getter,js_name=replaceExistingTags)]
    pub fn replace_existing_tags(&self) -> bool {
        self.inner.replace_existing_tags
    }
}
#[wasm_bindgen]
pub struct PotentialStereoParams {
    inner: ck::PotentialStereoParams,
}
#[cosmolkit_wasm::javascript_options("PotentialStereoParamsOptions")]
#[wasm_bindgen]
impl PotentialStereoParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] clean: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] flag_possible: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] allow_nontetrahedral: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::PotentialStereoParams::default();
        if !clean.is_undefined() {
            inner.clean = bool_value(&clean, "clean")?;
        }
        if !flag_possible.is_undefined() {
            inner.flag_possible = bool_value(&flag_possible, "flagPossible")?;
        }
        if !allow_nontetrahedral.is_undefined() {
            inner.allow_nontetrahedral = bool_value(&allow_nontetrahedral, "allowNontetrahedral")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter)]
    pub fn clean(&self) -> bool {
        self.inner.clean
    }
    #[wasm_bindgen(getter,js_name=flagPossible)]
    pub fn flag_possible(&self) -> bool {
        self.inner.flag_possible
    }
    #[wasm_bindgen(getter,js_name=allowNontetrahedral)]
    pub fn allow_nontetrahedral(&self) -> bool {
        self.inner.allow_nontetrahedral
    }
}
#[wasm_bindgen]
pub struct CipLabelOptions {
    inner: ck::CipLabelOptions,
}
#[cosmolkit_wasm::javascript_options("CipLabelOptionsOptions")]
#[wasm_bindgen]
impl CipLabelOptions {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | null")] atoms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | null")] bonds: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_recursive_iterations: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::CipLabelOptions::default();
        if let Some(atoms) = indices(&atoms, "atoms")? {
            inner = inner.with_atoms(atoms.into_iter().map(ck::AtomId::new));
        }
        if let Some(bonds) = indices(&bonds, "bonds")? {
            inner = inner.with_bonds(bonds.into_iter().map(ck::BondId::new));
        }
        if !max_recursive_iterations.is_undefined() {
            inner = inner.with_max_recursive_iterations(u32_value(
                &max_recursive_iterations,
                "maxRecursiveIterations",
            )?);
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number[] | null")]
    pub fn atoms(&self) -> JsValue {
        self.inner
            .atoms()
            .map_or(JsValue::NULL, |rows| ids(rows.iter().map(|id| id.index())))
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number[] | null")]
    pub fn bonds(&self) -> JsValue {
        self.inner
            .bonds()
            .map_or(JsValue::NULL, |rows| ids(rows.iter().map(|id| id.index())))
    }
    #[wasm_bindgen(getter,js_name=maxRecursiveIterations)]
    pub fn max_recursive_iterations(&self) -> u32 {
        self.inner.max_recursive_iterations()
    }
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum PotentialStereoType {
    AtomTetrahedral,
    AtomSquarePlanar,
    AtomTrigonalBipyramidal,
    AtomOctahedral,
    BondDouble,
    BondCumuleneEven,
    BondAtropisomer,
}
fn stereo_type(v: ck::PotentialStereoType) -> PotentialStereoType {
    match v {
        ck::PotentialStereoType::AtomTetrahedral => PotentialStereoType::AtomTetrahedral,
        ck::PotentialStereoType::AtomSquarePlanar => PotentialStereoType::AtomSquarePlanar,
        ck::PotentialStereoType::AtomTrigonalBipyramidal => {
            PotentialStereoType::AtomTrigonalBipyramidal
        }
        ck::PotentialStereoType::AtomOctahedral => PotentialStereoType::AtomOctahedral,
        ck::PotentialStereoType::BondDouble => PotentialStereoType::BondDouble,
        ck::PotentialStereoType::BondCumuleneEven => PotentialStereoType::BondCumuleneEven,
        ck::PotentialStereoType::BondAtropisomer => PotentialStereoType::BondAtropisomer,
    }
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum PotentialStereoSpecified {
    Unspecified,
    Specified,
    Unknown,
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum PotentialStereoDescriptor {
    None,
    TetrahedralClockwise,
    TetrahedralCounterclockwise,
    BondCis,
    BondTrans,
    BondAtropCw,
    BondAtropCcw,
}
fn descriptor(v: ck::PotentialStereoDescriptor) -> PotentialStereoDescriptor {
    match v {
        ck::PotentialStereoDescriptor::None => PotentialStereoDescriptor::None,
        ck::PotentialStereoDescriptor::TetrahedralClockwise => {
            PotentialStereoDescriptor::TetrahedralClockwise
        }
        ck::PotentialStereoDescriptor::TetrahedralCounterclockwise => {
            PotentialStereoDescriptor::TetrahedralCounterclockwise
        }
        ck::PotentialStereoDescriptor::BondCis => PotentialStereoDescriptor::BondCis,
        ck::PotentialStereoDescriptor::BondTrans => PotentialStereoDescriptor::BondTrans,
        ck::PotentialStereoDescriptor::BondAtropCw => PotentialStereoDescriptor::BondAtropCw,
        ck::PotentialStereoDescriptor::BondAtropCcw => PotentialStereoDescriptor::BondAtropCcw,
    }
}
#[wasm_bindgen]
pub struct PotentialStereoCenter {
    inner: ck::PotentialStereoCenter,
}
#[wasm_bindgen]
impl PotentialStereoCenter {
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        match self.inner {
            ck::PotentialStereoCenter::Atom(..) => "atom",
            ck::PotentialStereoCenter::Bond(..) => "bond",
        }
        .into()
    }
    #[wasm_bindgen(getter)]
    pub fn index(&self) -> usize {
        match self.inner {
            ck::PotentialStereoCenter::Atom(id) => id.index(),
            ck::PotentialStereoCenter::Bond(id) => id.index(),
        }
    }
}
#[wasm_bindgen]
pub struct PotentialStereoInfo {
    inner: ck::PotentialStereoInfo,
}
#[wasm_bindgen]
impl PotentialStereoInfo {
    #[wasm_bindgen(getter,js_name=stereoType)]
    pub fn stereo_type(&self) -> PotentialStereoType {
        stereo_type(self.inner.stereo_type)
    }
    #[wasm_bindgen(getter)]
    pub fn specified(&self) -> PotentialStereoSpecified {
        match self.inner.specified {
            ck::PotentialStereoSpecified::Unspecified => PotentialStereoSpecified::Unspecified,
            ck::PotentialStereoSpecified::Specified => PotentialStereoSpecified::Specified,
            ck::PotentialStereoSpecified::Unknown => PotentialStereoSpecified::Unknown,
        }
    }
    #[wasm_bindgen(getter,js_name=centeredOn)]
    pub fn centered_on(&self) -> PotentialStereoCenter {
        PotentialStereoCenter {
            inner: self.inner.centered_on,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn descriptor(&self) -> PotentialStereoDescriptor {
        descriptor(self.inner.descriptor)
    }
    #[wasm_bindgen(getter)]
    pub fn permutation(&self) -> u32 {
        self.inner.permutation
    }
    #[wasm_bindgen(getter,js_name=controllingAtoms,unchecked_return_type="Array<number | null>")]
    pub fn controlling_atoms(&self) -> JsValue {
        self.inner
            .controlling_atoms
            .iter()
            .map(|id| id.map_or(JsValue::NULL, |id| JsValue::from(id.index() as u32)))
            .collect::<Array>()
            .into()
    }
}
#[wasm_bindgen]
pub struct RingStereoRelation {
    inner: ck::RingStereoRelation,
}
#[wasm_bindgen]
impl RingStereoRelation {
    #[wasm_bindgen(getter)]
    pub fn atom(&self) -> usize {
        self.inner.atom.index()
    }
    #[wasm_bindgen(getter)]
    pub fn other(&self) -> usize {
        self.inner.other.index()
    }
    #[wasm_bindgen(getter,js_name=sameOrientation)]
    pub fn same_orientation(&self) -> bool {
        self.inner.same_orientation
    }
}
#[wasm_bindgen]
pub struct PotentialStereoResult {
    inner: ck::PotentialStereoResult,
}
#[wasm_bindgen]
impl PotentialStereoResult {
    #[wasm_bindgen(getter, unchecked_return_type = "PotentialStereoInfo[]")]
    pub fn stereo(&self) -> JsValue {
        self.inner
            .stereo
            .iter()
            .cloned()
            .map(|inner| JsValue::from(PotentialStereoInfo { inner }))
            .collect::<Array>()
            .into()
    }
    #[wasm_bindgen(getter,js_name=atomRanks,unchecked_return_type="number[]")]
    pub fn atom_ranks(&self) -> JsValue {
        self.inner
            .atom_ranks
            .iter()
            .map(|v| JsValue::from(*v))
            .collect::<Array>()
            .into()
    }
    #[wasm_bindgen(getter,js_name=ringRelations,unchecked_return_type="RingStereoRelation[]")]
    pub fn ring_relations(&self) -> JsValue {
        self.inner
            .ring_relations
            .iter()
            .copied()
            .map(|inner| JsValue::from(RingStereoRelation { inner }))
            .collect::<Array>()
            .into()
    }
    #[wasm_bindgen(getter,js_name=cleanedMolecule,unchecked_return_type="Molecule | null")]
    pub fn cleaned_molecule(&self) -> JsValue {
        self.inner
            .cleaned_molecule
            .clone()
            .map_or(JsValue::NULL, |inner| {
                Molecule {
                    inner: cosmolkit_wasm::Molecule::from(inner).into(),
                }
                .into()
            })
    }
}
#[wasm_bindgen]
pub struct LigandRef {
    inner: ck::LigandRef,
}
#[wasm_bindgen]
impl LigandRef {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] atom: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: if atom.is_null() || atom.is_undefined() {
                ck::LigandRef::ImplicitHydrogen
            } else {
                ck::LigandRef::Atom(ck::AtomId::new(usize_value(&atom, "atom")?))
            },
        })
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn atom(&self) -> JsValue {
        match self.inner {
            ck::LigandRef::Atom(id) => JsValue::from(id.index() as u32),
            ck::LigandRef::ImplicitHydrogen => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=isImplicitHydrogen)]
    pub fn is_implicit_hydrogen(&self) -> bool {
        matches!(self.inner, ck::LigandRef::ImplicitHydrogen)
    }
}
#[wasm_bindgen]
pub struct TetrahedralStereo {
    inner: ck::TetrahedralStereo,
}
#[wasm_bindgen]
impl TetrahedralStereo {
    #[wasm_bindgen(getter)]
    pub fn center(&self) -> usize {
        self.inner.center.index()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "LigandRef[]")]
    pub fn ligands(&self) -> JsValue {
        self.inner
            .ligands
            .iter()
            .copied()
            .map(|inner| JsValue::from(LigandRef { inner }))
            .collect::<Array>()
            .into()
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=potentialStereo)]
    pub fn potential_stereo(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "PotentialStereoParams | PotentialStereoParamsOptions"
        )]
        params: JsValue,
    ) -> Result<PotentialStereoResult, JsValue> {
        self.potential_stereo_with_params(&PotentialStereoParams::from_js_configuration(&params)?)
    }
    #[wasm_bindgen(js_name=potentialStereoWithParams)]
    pub fn potential_stereo_with_params(
        &self,
        params: &PotentialStereoParams,
    ) -> Result<PotentialStereoResult, JsValue> {
        operation(self.inner.potential_stereo_with_params(&params.inner))
            .map(|inner| PotentialStereoResult { inner })
    }
    #[wasm_bindgen(js_name=withChiralTagsFromStructure)]
    pub fn with_chiral_tags_from_structure(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "StructureTagParams | StructureTagParamsOptions"
        )]
        params: JsValue,
    ) -> Result<Molecule, JsValue> {
        self.with_chiral_tags_from_structure_with_params(
            &StructureTagParams::from_js_configuration(&params)?,
        )
    }
    #[wasm_bindgen(js_name=withChiralTagsFromStructureWithParams)]
    pub fn with_chiral_tags_from_structure_with_params(
        &self,
        params: &StructureTagParams,
    ) -> Result<Molecule, JsValue> {
        operation(
            self.inner
                .with_chiral_tags_from_structure_with_params(&params.inner),
        )
        .map(|inner| Molecule {
            inner: inner.into(),
        })
    }
    #[wasm_bindgen(js_name=assignChiralTagsFromStructure_)]
    pub fn assign_chiral_tags_from_structure_(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "StructureTagParams | StructureTagParamsOptions"
        )]
        params: JsValue,
    ) -> Result<(), JsValue> {
        self.assign_chiral_tags_from_structure_with_params_(
            &StructureTagParams::from_js_configuration(&params)?,
        )
    }
    #[wasm_bindgen(js_name=assignChiralTagsFromStructureWithParams_)]
    pub fn assign_chiral_tags_from_structure_with_params_(
        &self,
        params: &StructureTagParams,
    ) -> Result<(), JsValue> {
        operation(
            self.inner
                .assign_chiral_tags_from_structure_with_params_(&params.inner),
        )
    }
    #[wasm_bindgen(js_name=withCipLabels)]
    pub fn with_cip_labels(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "CipLabelOptions | CipLabelOptionsOptions")]
        options: JsValue,
    ) -> Result<Molecule, JsValue> {
        self.with_cip_labels_with_options(&CipLabelOptions::from_js_configuration(&options)?)
    }
    #[wasm_bindgen(js_name=withCipLabelsWithOptions)]
    pub fn with_cip_labels_with_options(
        &self,
        options: &CipLabelOptions,
    ) -> Result<Molecule, JsValue> {
        operation(self.inner.with_cip_labels_with_options(&options.inner)).map(|inner| Molecule {
            inner: inner.into(),
        })
    }
    #[wasm_bindgen(js_name=assignCipLabels_)]
    pub fn assign_cip_labels_(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "CipLabelOptions | CipLabelOptionsOptions")]
        options: JsValue,
    ) -> Result<(), JsValue> {
        self.assign_cip_labels_with_options_(&CipLabelOptions::from_js_configuration(&options)?)
    }
    #[wasm_bindgen(js_name=assignCipLabelsWithOptions_)]
    pub fn assign_cip_labels_with_options_(
        &self,
        options: &CipLabelOptions,
    ) -> Result<(), JsValue> {
        operation(self.inner.assign_cip_labels_with_options_(&options.inner))
    }
    #[wasm_bindgen(js_name=cipComputed)]
    pub fn cip_computed(&self) -> Result<bool, JsValue> {
        self.inner
            .cip_computed()
            .map_err(|e| crate::property_values::property_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=tetrahedralStereo,unchecked_return_type="TetrahedralStereo[]")]
    pub fn tetrahedral_stereo(&self) -> Result<JsValue, JsValue> {
        self.inner
            .tetrahedral_stereo()
            .map(|v| {
                v.into_iter()
                    .map(|inner| JsValue::from(TetrahedralStereo { inner }))
                    .collect::<Array>()
                    .into()
            })
            .map_err(|e| crate::stereo_queries::read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=perceiveStereochemistry)]
    pub fn perceive_stereochemistry(&self) -> Result<(), JsValue> {
        self.inner
            .perceive_stereochemistry()
            .map_err(|e| crate::stereo_queries::read_error(&e).unwrap_or_else(|e| e))
    }
}
