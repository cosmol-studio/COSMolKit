//! Native alignment inputs and reports; matrix dimensions are checked on input.
use crate::{fields::Field, input, storage::NativeField};
use cosmolkit as ck;
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct AlignmentAtomMap(pub(crate) ck::AlignmentAtomMap);
crate::fields::value_field!(
    ck::AlignmentAtomMap,
    AlignmentAtomMap,
    AlignmentAtomMapInput,
    "AlignmentAtomMap"
);
crate::fields::value_array!(
    ck::AlignmentAtomMap,
    AlignmentAtomMaps,
    "ReadonlyArray<AlignmentAtomMap>"
);
crate::fields::value_array!(
    Vec<ck::AlignmentAtomMap>,
    AlignmentMaps,
    "ReadonlyArray<ReadonlyArray<AlignmentAtomMap>>"
);
#[wasm_bindgen]
impl AlignmentAtomMap {
    #[wasm_bindgen(constructor)]
    pub fn new(
        probe_atom: Ts<input::IndexInput>,
        reference_atom: Ts<input::IndexInput>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::AlignmentAtomMap::new(
            input::read(probe_atom)?,
            input::read(reference_atom)?,
        )))
    }
    #[wasm_bindgen(getter,js_name=probeAtom)]
    pub fn probe_atom(&self) -> usize {
        self.0.probe_atom
    }
    #[wasm_bindgen(getter,js_name=referenceAtom)]
    pub fn reference_atom(&self) -> usize {
        self.0.reference_atom
    }
}
ck::__ck_configuration_AlignmentParameters!(
    crate::project_parameters,
    ck::AlignmentParameters,
    {},
    { AlignmentParameters }
);
ck::__ck_configuration_BestAlignmentParameters!(
    crate::project_parameters,
    ck::BestAlignmentParameters,
    {},
    { BestAlignmentParameters }
);
ck::__ck_configuration_CoordinateRmsdParameters!(
    crate::project_parameters,
    ck::CoordinateRmsdParameters,
    {},
    { CoordinateRmsdParameters }
);
ck::__ck_configuration_AllConformerRmsdParameters!(
    crate::project_parameters,
    ck::AllConformerRmsdParameters,
    {},
    { AllConformerRmsdParameters }
);
ck::__ck_configuration_ConformerAlignmentParameters!(
    crate::project_parameters,
    ck::ConformerAlignmentParameters,
    {},
    { ConformerAlignmentParameters }
);
#[derive(serde::Serialize, serde::Deserialize, Tsify)]
#[serde(transparent)]
pub struct TransformRow(#[serde(deserialize_with = "input::exact_array")] [f64; 4]);
#[derive(serde::Serialize, serde::Deserialize, Tsify)]
#[serde(transparent)]
pub struct TransformMatrix(#[serde(deserialize_with = "input::exact_array")] [TransformRow; 4]);
#[wasm_bindgen]
pub struct AlignmentTransform(pub(crate) ck::AlignmentTransform);
crate::fields::value_field!(
    ck::AlignmentTransform,
    AlignmentTransform,
    AlignmentTransformInput,
    "AlignmentTransform"
);
#[wasm_bindgen]
impl AlignmentTransform {
    #[wasm_bindgen(constructor)]
    pub fn new(matrix: Ts<TransformMatrix>) -> Result<Self, JsValue> {
        Ok(Self(ck::AlignmentTransform {
            matrix: input::decode(matrix)?.0.map(|row| row.0),
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn matrix(&self) -> Result<Ts<TransformMatrix>, JsValue> {
        TransformMatrix(self.0.matrix.map(TransformRow))
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
}
#[wasm_bindgen]
pub struct AlignmentResult(pub(crate) ck::AlignmentResult);
crate::fields::value_field!(
    ck::AlignmentResult,
    AlignmentResult,
    AlignmentResultInput,
    "AlignmentResult"
);
#[wasm_bindgen]
impl AlignmentResult {
    #[wasm_bindgen(constructor)]
    pub fn new(
        rmsd: Ts<input::F64Input>,
        transform: AlignmentTransformInput,
        atom_map: AlignmentAtomMaps,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::AlignmentResult {
            rmsd: input::read(rmsd)?,
            transform: ck::AlignmentTransform::input(transform)?,
            atom_map: Vec::<ck::AlignmentAtomMap>::input(atom_map)?,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn rmsd(&self) -> f64 {
        self.0.rmsd()
    }
    #[wasm_bindgen(getter)]
    pub fn transform(&self) -> AlignmentTransform {
        AlignmentTransform(*self.0.transform())
    }
    #[wasm_bindgen(getter,js_name=atomMap)]
    pub fn atom_map(&self) -> Result<AlignmentAtomMaps, JsValue> {
        Vec::<ck::AlignmentAtomMap>::get(NativeField::owned(self.0.atom_map().to_vec()))
    }
}
#[wasm_bindgen]
pub struct ConformerRmsd(pub(crate) ck::ConformerRmsd);
crate::fields::value_field!(
    ck::ConformerRmsd,
    ConformerRmsd,
    ConformerRmsdInput,
    "ConformerRmsd"
);
crate::fields::value_array!(
    ck::ConformerRmsd,
    ConformerRmsds,
    "ReadonlyArray<ConformerRmsd>"
);
#[wasm_bindgen]
impl ConformerRmsd {
    #[wasm_bindgen(constructor)]
    pub fn new(
        probe_conformer_id: Ts<input::IndexInput>,
        reference_conformer_id: Ts<input::IndexInput>,
        rmsd: Ts<input::F64Input>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::ConformerRmsd {
            probe_conformer_id: input::read(probe_conformer_id)?,
            reference_conformer_id: input::read(reference_conformer_id)?,
            rmsd: input::read(rmsd)?,
        }))
    }
    #[wasm_bindgen(getter,js_name=probeConformerId)]
    pub fn probe_conformer_id(&self) -> usize {
        self.0.probe_conformer_id()
    }
    #[wasm_bindgen(getter,js_name=referenceConformerId)]
    pub fn reference_conformer_id(&self) -> usize {
        self.0.reference_conformer_id()
    }
    #[wasm_bindgen(getter)]
    pub fn rmsd(&self) -> f64 {
        self.0.rmsd()
    }
}
#[wasm_bindgen]
pub struct ConformerAlignmentReport(pub(crate) ck::ConformerAlignmentReport);
#[wasm_bindgen]
impl ConformerAlignmentReport {
    #[wasm_bindgen(constructor)]
    pub fn new(rmsds: <Vec<f64> as Field>::Input) -> Result<Self, JsValue> {
        Ok(Self(ck::ConformerAlignmentReport {
            rmsds: Vec::<f64>::input(rmsds)?,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn rmsds(&self) -> Result<<Vec<f64> as Field>::Output, JsValue> {
        Vec::<f64>::get(NativeField::owned(self.0.rmsds().to_vec()))
    }
}
