//! Configuration declarations and value/in-place native assignment operations.
use crate::{Molecule, fields::Field, operation_errors::OperationError, storage::NativeField};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(getter,js_name=cipComputed)]
    pub fn cip_computed(&self) -> Result<bool, JsValue> {
        self.0
            .cip_computed()
            .map_err(crate::model_errors::PropertyValueError::from_rust)
    }
}
ck::__ck_configuration_ValenceParams!(crate::project_parameters, ck::ValenceParams, {}, {
    ValenceParams
});
ck::__ck_configuration_AromaticityParams!(crate::project_parameters, ck::AromaticityParams, {}, {
    AromaticityParams
});
ck::__ck_configuration_CipLabelOptions!(crate::project_parameters, ck::CipLabelOptions, {}, {
    CipLabelOptions
});
ck::__ck_configuration_PotentialStereoParams!(
    crate::project_parameters,
    ck::PotentialStereoParams,
    {},
    { PotentialStereoParams }
);
macro_rules! assignments{($(($params:ident,$input:ident,$with:ident,$js:ident,$native:ident,$set:ident,$setjs:ident,$nativeset:ident)),+)=>{#[wasm_bindgen]impl Molecule{$(
 #[wasm_bindgen(js_name=$js)]pub fn $with(&self,params:Option<$input>)->Result<Self,JsValue>{$params::with_input(params.as_ref(),|p|self.0.$native(p).map(Self).map_err(OperationError::from_rust))}
 #[wasm_bindgen(js_name=$setjs)]pub fn $set(&mut self,params:Option<$input>)->Result<(),JsValue>{$params::with_input(params.as_ref(),|p|self.0.$nativeset(p).map_err(OperationError::from_rust))}
)+}};}
assignments!(
    (
        ValenceParams,
        ValenceParamsInput,
        with_assigned_valence,
        withAssignedValence,
        with_assigned_valence_with_params,
        assign_valence_,
        assignValence_,
        assign_valence_with_params_
    ),
    (
        AromaticityParams,
        AromaticityParamsInput,
        with_assigned_aromaticity,
        withAssignedAromaticity,
        with_assigned_aromaticity_with_params,
        assign_aromaticity_,
        assignAromaticity_,
        assign_aromaticity_with_params_
    ),
    (
        CipLabelOptions,
        CipLabelOptionsInput,
        with_cip_labels,
        withCipLabels,
        with_cip_labels_with_options,
        assign_cip_labels_,
        assignCipLabels_,
        assign_cip_labels_with_options_
    )
);
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=withAssignedRadicals)]
    pub fn with_assigned_radicals(&self) -> Result<Self, JsValue> {
        self.0
            .with_assigned_radicals()
            .map(Self)
            .map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=assignRadicals_)]
    pub fn assign_radicals(&mut self) -> Result<(), JsValue> {
        self.0.assign_radicals_().map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=potentialStereo)]
    pub fn potential_stereo(
        &self,
        params: Option<PotentialStereoParamsInput>,
    ) -> Result<PotentialStereoResult, JsValue> {
        PotentialStereoParams::with_input(params.as_ref(), |p| {
            self.0
                .potential_stereo_with_params(p)
                .map(PotentialStereoResult)
                .map_err(OperationError::from_rust)
        })
    }
}
#[derive(serde::Serialize, serde::Deserialize, tsify::Tsify)]
#[serde(tag = "kind", content = "index", rename_all = "camelCase")]
pub enum PotentialStereoCenter {
    Atom(usize),
    Bond(usize),
}
#[wasm_bindgen]
pub struct PotentialStereoInfo(pub(crate) ck::PotentialStereoInfo);
crate::fields::value_field!(
    ck::PotentialStereoInfo,
    PotentialStereoInfo,
    PotentialStereoInfoInput,
    "PotentialStereoInfo"
);
crate::fields::value_array!(
    ck::PotentialStereoInfo,
    PotentialStereoInfos,
    "ReadonlyArray<PotentialStereoInfo>"
);
#[wasm_bindgen]
impl PotentialStereoInfo {
    #[wasm_bindgen(constructor)]
    pub fn new(
        stereo_type: crate::enums::PotentialStereoTypeInput,
        specified: crate::enums::PotentialStereoSpecifiedInput,
        center: tsify::Ts<PotentialStereoCenter>,
        descriptor: crate::enums::PotentialStereoDescriptorInput,
        permutation: tsify::Ts<crate::input::U32Input>,
        controlling_atoms: tsify::Ts<ControllingAtoms>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::PotentialStereoInfo {
            stereo_type: ck::PotentialStereoType::input(stereo_type)?,
            specified: ck::PotentialStereoSpecified::input(specified)?,
            centered_on: match crate::input::decode(center)? {
                PotentialStereoCenter::Atom(i) => {
                    ck::PotentialStereoCenter::Atom(ck::AtomId::new(i))
                }
                PotentialStereoCenter::Bond(i) => {
                    ck::PotentialStereoCenter::Bond(ck::BondId::new(i))
                }
            },
            descriptor: ck::PotentialStereoDescriptor::input(descriptor)?,
            permutation: crate::input::read(permutation)?,
            controlling_atoms: crate::input::decode(controlling_atoms)?
                .0
                .into_iter()
                .map(|v| v.map(ck::AtomId::new))
                .collect(),
        }))
    }
    #[wasm_bindgen(getter,js_name=stereoType)]
    pub fn stereo_type(&self) -> crate::enums::PotentialStereoType {
        self.0.stereo_type.into()
    }
    #[wasm_bindgen(getter)]
    pub fn specified(&self) -> crate::enums::PotentialStereoSpecified {
        self.0.specified.into()
    }
    #[wasm_bindgen(getter)]
    pub fn descriptor(&self) -> crate::enums::PotentialStereoDescriptor {
        self.0.descriptor.into()
    }
    #[wasm_bindgen(getter)]
    pub fn permutation(&self) -> u32 {
        self.0.permutation
    }
    #[wasm_bindgen(getter,js_name=centeredOn)]
    pub fn centered_on(&self) -> Result<tsify::Ts<PotentialStereoCenter>, JsValue> {
        use tsify::Tsify;
        match self.0.centered_on {
            ck::PotentialStereoCenter::Atom(v) => PotentialStereoCenter::Atom(v.index()),
            ck::PotentialStereoCenter::Bond(v) => PotentialStereoCenter::Bond(v.index()),
        }
        .into_ts()
        .map_err(wasm_bindgen::JsError::from)
        .map_err(Into::into)
    }
    #[wasm_bindgen(getter,js_name=controllingAtoms)]
    pub fn controlling_atoms(&self) -> Result<tsify::Ts<ControllingAtoms>, JsValue> {
        use tsify::Tsify;
        ControllingAtoms(
            self.0
                .controlling_atoms
                .iter()
                .map(|v| v.map(ck::AtomId::index))
                .collect(),
        )
        .into_ts()
        .map_err(wasm_bindgen::JsError::from)
        .map_err(Into::into)
    }
}
#[derive(serde::Serialize, serde::Deserialize, tsify::Tsify)]
#[serde(transparent)]
pub struct ControllingAtoms(Vec<Option<usize>>);
#[wasm_bindgen]
pub struct RingStereoRelation(pub(crate) ck::RingStereoRelation);
crate::fields::value_field!(
    ck::RingStereoRelation,
    RingStereoRelation,
    RingStereoRelationInput,
    "RingStereoRelation"
);
crate::fields::value_array!(
    ck::RingStereoRelation,
    RingStereoRelations,
    "ReadonlyArray<RingStereoRelation>"
);
#[wasm_bindgen]
impl RingStereoRelation {
    #[wasm_bindgen(constructor)]
    pub fn new(
        atom: tsify::Ts<crate::input::IndexInput>,
        other: tsify::Ts<crate::input::IndexInput>,
        same_orientation: tsify::Ts<crate::input::BoolInput>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::RingStereoRelation {
            atom: ck::AtomId::new(crate::input::read(atom)?),
            other: ck::AtomId::new(crate::input::read(other)?),
            same_orientation: crate::input::read(same_orientation)?,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn atom(&self) -> usize {
        self.0.atom.index()
    }
    #[wasm_bindgen(getter)]
    pub fn other(&self) -> usize {
        self.0.other.index()
    }
    #[wasm_bindgen(getter,js_name=sameOrientation)]
    pub fn same_orientation(&self) -> bool {
        self.0.same_orientation
    }
}
#[wasm_bindgen]
pub struct PotentialStereoResult(ck::PotentialStereoResult);
#[wasm_bindgen]
impl PotentialStereoResult {
    #[wasm_bindgen(constructor)]
    pub fn new(
        stereo: PotentialStereoInfos,
        atom_ranks: <Vec<u32> as Field>::Input,
        ring_relations: RingStereoRelations,
        cleaned_molecule: Option<crate::MoleculeInput>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::PotentialStereoResult {
            stereo: Vec::<ck::PotentialStereoInfo>::input(stereo)?,
            atom_ranks: Vec::<u32>::input(atom_ranks)?,
            ring_relations: Vec::<ck::RingStereoRelation>::input(ring_relations)?,
            cleaned_molecule: cleaned_molecule.map(ck::Molecule::input).transpose()?,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn stereo(&self) -> Result<PotentialStereoInfos, JsValue> {
        Vec::<ck::PotentialStereoInfo>::get(NativeField::owned(self.0.stereo.clone()))
    }
    #[wasm_bindgen(getter,js_name=atomRanks)]
    pub fn atom_ranks(&self) -> Vec<u32> {
        self.0.atom_ranks.clone()
    }
    #[wasm_bindgen(getter,js_name=ringRelations)]
    pub fn ring_relations(&self) -> Result<RingStereoRelations, JsValue> {
        Vec::<ck::RingStereoRelation>::get(NativeField::owned(self.0.ring_relations.clone()))
    }
    #[wasm_bindgen(getter,js_name=cleanedMolecule)]
    pub fn cleaned_molecule(&self) -> Option<Molecule> {
        self.0.cleaned_molecule.as_ref().map(|m| Molecule(m.copy()))
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = potentialStereoWithParams)]
    pub fn potential_stereo_with_params(
        &self,
        params: PotentialStereoParamsInput,
    ) -> Result<PotentialStereoResult, JsValue> {
        self.potential_stereo(Some(params))
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=withAssignedValenceWithParams)]
    pub fn with_assigned_valence_with_params(
        &self,
        params: ValenceParamsInput,
    ) -> Result<Self, JsValue> {
        self.with_assigned_valence(Some(params))
    }
    #[wasm_bindgen(js_name=assignValenceWithParams_)]
    pub fn assign_valence_with_params_(
        &mut self,
        params: ValenceParamsInput,
    ) -> Result<(), JsValue> {
        self.assign_valence_(Some(params))
    }
    #[wasm_bindgen(js_name=withAssignedAromaticityWithParams)]
    pub fn with_assigned_aromaticity_with_params(
        &self,
        params: AromaticityParamsInput,
    ) -> Result<Self, JsValue> {
        self.with_assigned_aromaticity(Some(params))
    }
    #[wasm_bindgen(js_name=assignAromaticityWithParams_)]
    pub fn assign_aromaticity_with_params_(
        &mut self,
        params: AromaticityParamsInput,
    ) -> Result<(), JsValue> {
        self.assign_aromaticity_(Some(params))
    }
    #[wasm_bindgen(js_name=withCipLabelsWithOptions)]
    pub fn with_cip_labels_with_options(
        &self,
        params: CipLabelOptionsInput,
    ) -> Result<Self, JsValue> {
        self.with_cip_labels(Some(params))
    }
    #[wasm_bindgen(js_name=assignCipLabelsWithOptions_)]
    pub fn assign_cip_labels_with_options_(
        &mut self,
        params: CipLabelOptionsInput,
    ) -> Result<(), JsValue> {
        self.assign_cip_labels_(Some(params))
    }
}
