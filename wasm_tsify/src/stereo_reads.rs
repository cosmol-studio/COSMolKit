//! Native read-only stereo values and CIP labels.
use crate::{
    Molecule,
    errors::{enum_error, exception},
    fields::Field,
    input,
    storage::NativeField,
};
use cosmolkit as ck;
use tsify::{Ts, Tsify};
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | TemplateAttachmentOrderError | ValenceError | RingFindingError | PotentialStereoError | CipLabelerError | PropertyStringError | FromUtf8Error"
    )]
    #[derive(Clone)]
    pub type StereoReadCause;
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | ValenceError | RingFindingError | PotentialStereoError | CipLabelerError | PropertyStringError | FromUtf8Error"
    )]
    pub type StereoReadPayload;
    #[wasm_bindgen(typescript_type = "readonly [LigandRef, LigandRef, LigandRef, LigandRef]")]
    pub type LigandTuple;
}
fn cause(e: &ck::StereoReadError) -> Option<StereoReadCause> {
    use ck::StereoReadError as E;
    Some(
        match e {
            E::InvalidTopology(v) => {
                return crate::model_errors::topology_cause(v).map(JsCast::unchecked_into);
            }
            E::Valence(v) => crate::valence_errors::ValenceError::from_rust(v.clone()),
            E::Rings(v) => crate::ring_errors::RingFindingError::from_rust(v.clone()),
            E::PotentialStereo(v) => {
                crate::potential_stereo_errors::PotentialStereoError::from_rust(v.clone())
            }
            E::CipLabeler(v) => crate::cip_labeler_errors::CipLabelerError::from_rust(v.clone()),
            E::PropertyString(v) => crate::property_errors::PropertyStringError::from_rust(*v),
            E::CipLabelEncoding { source, .. } => {
                crate::external_errors::FromUtf8Error::from_rust(source.clone())
            }
        }
        .unchecked_into(),
    )
}
enum_error!(StereoReadError,StereoReadErrorKind,ck::StereoReadError,{
InvalidTopology=>ck::StereoReadError::InvalidTopology(_),Valence=>ck::StereoReadError::Valence(_),Rings=>ck::StereoReadError::Rings(_),PotentialStereo=>ck::StereoReadError::PotentialStereo(_),CipLabeler=>ck::StereoReadError::CipLabeler(_),PropertyString=>ck::StereoReadError::PropertyString(_),CipLabelEncoding=>ck::StereoReadError::CipLabelEncoding{..}
},{
 #[wasm_bindgen(getter)]pub fn atom(&self)->Option<usize>{match *self.inner{ck::StereoReadError::CipLabelEncoding{atom,..}=>Some(atom.index()),_=>None}}
 #[wasm_bindgen(getter)]pub fn source(&self)->StereoReadPayload{use ck::StereoReadError as E;match self.inner.as_ref(){E::InvalidTopology(v)=>crate::model_errors::TopologyValidationError::from_rust(v.clone()),E::Valence(v)=>crate::valence_errors::ValenceError::from_rust(v.clone()),E::Rings(v)=>crate::ring_errors::RingFindingError::from_rust(v.clone()),E::PotentialStereo(v)=>crate::potential_stereo_errors::PotentialStereoError::from_rust(v.clone()),E::CipLabeler(v)=>crate::cip_labeler_errors::CipLabelerError::from_rust(v.clone()),E::PropertyString(v)=>crate::property_errors::PropertyStringError::from_rust(*v),E::CipLabelEncoding{source,..}=>crate::external_errors::FromUtf8Error::from_rust(source.clone())}.unchecked_into()}
},source=StereoReadCause,cause=cause);
#[wasm_bindgen]
pub struct LigandRef(pub(crate) ck::LigandRef);
crate::fields::value_field!(ck::LigandRef, LigandRef, LigandRefInput, "LigandRef");
crate::fields::value_array!(ck::LigandRef, LigandRefs, "ReadonlyArray<LigandRef>");
#[wasm_bindgen]
impl LigandRef {
    #[wasm_bindgen(constructor)]
    pub fn new(atom: Option<Ts<input::IndexInput>>) -> Result<Self, JsValue> {
        Ok(Self(match atom {
            Some(v) => ck::LigandRef::Atom(ck::AtomId::new(input::read(v)?)),
            None => ck::LigandRef::ImplicitHydrogen,
        }))
    }
    #[wasm_bindgen(js_name=Atom)]
    pub fn atom(atom: Ts<input::IndexInput>) -> Result<Self, JsValue> {
        Self::new(Some(atom))
    }
    #[wasm_bindgen(getter,js_name=ImplicitHydrogen)]
    pub fn implicit_hydrogen() -> Self {
        Self(ck::LigandRef::ImplicitHydrogen)
    }
    #[wasm_bindgen(getter,js_name=atomId)]
    pub fn atom_id(&self) -> Option<usize> {
        match self.0 {
            ck::LigandRef::Atom(v) => Some(v.index()),
            ck::LigandRef::ImplicitHydrogen => None,
        }
    }
}
#[wasm_bindgen]
pub struct TetrahedralStereo(pub(crate) ck::TetrahedralStereo);
crate::fields::value_field!(
    ck::TetrahedralStereo,
    TetrahedralStereo,
    TetrahedralStereoInput,
    "TetrahedralStereo"
);
crate::fields::value_array!(
    ck::TetrahedralStereo,
    TetrahedralStereoArray,
    "ReadonlyArray<TetrahedralStereo>"
);
#[wasm_bindgen]
impl TetrahedralStereo {
    #[wasm_bindgen(constructor)]
    pub fn new(center: Ts<input::IndexInput>, ligands: LigandTuple) -> Result<Self, JsValue> {
        let ligands = Vec::<ck::LigandRef>::decode(ligands.as_ref())?
            .try_into()
            .map_err(|_: Vec<ck::LigandRef>| js_sys::TypeError::new("expected four ligands"))?;
        Ok(Self(ck::TetrahedralStereo {
            center: ck::AtomId::new(input::read(center)?),
            ligands,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn center(&self) -> usize {
        self.0.center.index()
    }
    #[wasm_bindgen(getter)]
    pub fn ligands(&self) -> Result<LigandTuple, JsValue> {
        Vec::<ck::LigandRef>::get(NativeField::owned(self.0.ligands.to_vec()))
            .map(JsCast::unchecked_into)
    }
}
#[derive(serde::Serialize, Tsify)]
#[serde(transparent)]
pub struct ChiralCenters(Vec<(usize, String)>);
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=tetrahedralStereo)]
    pub fn tetrahedral_stereo(&self) -> Result<TetrahedralStereoArray, JsValue> {
        Vec::<ck::TetrahedralStereo>::get(NativeField::owned(
            self.0
                .tetrahedral_stereo()
                .map_err(StereoReadError::from_rust)?,
        ))
    }
    #[wasm_bindgen(js_name=perceiveStereochemistry)]
    pub fn perceive_stereochemistry(&self) -> Result<(), JsValue> {
        self.0
            .perceive_stereochemistry()
            .map_err(StereoReadError::from_rust)
    }
    #[wasm_bindgen(js_name=findChiralCenters)]
    pub fn find_chiral_centers(
        &self,
        include_unassigned: Option<Ts<input::BoolInput>>,
    ) -> Result<Ts<ChiralCenters>, JsValue> {
        ChiralCenters(
            self.0
                .find_chiral_centers(
                    include_unassigned
                        .map(input::read::<bool>)
                        .transpose()?
                        .unwrap_or(false),
                )
                .map_err(StereoReadError::from_rust)?,
        )
        .into_ts()
        .map_err(wasm_bindgen::JsError::from)
        .map_err(Into::into)
    }
}
