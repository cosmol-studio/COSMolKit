//! BIO-to-molecule conversion preserves transparent native source semantics.
use crate::errors::{enum_error, exception};
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "AtomPropertyError | TopologyValidationError | CoordinateValidationError | ValenceError | TopologyEditError | PropertyUIntReadError"
    )]
    #[derive(Clone)]
    pub type PdbPostCause;
}
fn post_source(owner: &ReadonlyValue<ck::PdbPostprocessError>) -> Option<PdbPostCause> {
    Some(
        match owner.as_ref() {
            ck::PdbPostprocessError::AtomProperty(_) => {
                crate::model_errors::AtomPropertyError::from_owner(owner.child(|e| match e {
                    ck::PdbPostprocessError::AtomProperty(v) => v,
                    _ => unreachable!("immutable PDB postprocess variant"),
                }))
            }
            ck::PdbPostprocessError::Topology(_) => {
                crate::model_errors::TopologyValidationError::from_owner(owner.child(|e| match e {
                    ck::PdbPostprocessError::Topology(v) => v,
                    _ => unreachable!("immutable PDB postprocess variant"),
                }))
            }
            ck::PdbPostprocessError::Coordinates(_) => {
                crate::model_errors::CoordinateValidationError::from_owner(owner.child(
                    |e| match e {
                        ck::PdbPostprocessError::Coordinates(v) => v,
                        _ => unreachable!("immutable PDB postprocess variant"),
                    },
                ))
            }
            ck::PdbPostprocessError::Valence(_) => {
                crate::valence_errors::ValenceError::from_owner(owner.child(|e| match e {
                    ck::PdbPostprocessError::Valence(v) => v,
                    _ => unreachable!("immutable PDB postprocess variant"),
                }))
            }
            ck::PdbPostprocessError::TopologyEdit(_) => {
                crate::topology_errors::TopologyEditError::from_owner(owner.child(|e| match e {
                    ck::PdbPostprocessError::TopologyEdit(v) => v,
                    _ => unreachable!("immutable PDB postprocess variant"),
                }))
            }
            ck::PdbPostprocessError::PropertyUInt(_) => {
                crate::property_errors::PropertyUIntReadError::from_owner(owner.child(
                    |e| match e {
                        ck::PdbPostprocessError::PropertyUInt(v) => v,
                        _ => unreachable!("immutable PDB postprocess variant"),
                    },
                ))
            }
            ck::PdbPostprocessError::Unsupported(_) | ck::PdbPostprocessError::InvalidState(_) => {
                return None;
            }
        }
        .unchecked_into(),
    )
}
enum_error!(PdbPostprocessError,PdbPostprocessErrorKind,ck::PdbPostprocessError,{AtomProperty=>ck::PdbPostprocessError::AtomProperty(_),Topology=>ck::PdbPostprocessError::Topology(_),Coordinates=>ck::PdbPostprocessError::Coordinates(_),Valence=>ck::PdbPostprocessError::Valence(_),TopologyEdit=>ck::PdbPostprocessError::TopologyEdit(_),PropertyUInt=>ck::PdbPostprocessError::PropertyUInt(_),Unsupported=>ck::PdbPostprocessError::Unsupported(_),InvalidState=>ck::PdbPostprocessError::InvalidState(_)},{
#[wasm_bindgen(getter)]pub fn unsupported(&self)->Option<String>{match self.inner.as_ref(){ck::PdbPostprocessError::Unsupported(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=invalidState)]pub fn invalid_state(&self)->Option<String>{match self.inner.as_ref(){ck::PdbPostprocessError::InvalidState(v)=>Some(v.clone()),_=>None}}
},source=PdbPostCause,cause=post_source);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "AtomPropertyError | TopologyValidationError | CoordinateValidationError | ValenceError | TopologyEditError | PropertyUIntReadError | MoleculePropertyError | BondValueError | QueryStateError | CleanupError | PropertyCacheError | RingFindingError | KekulizeError | RadicalError | AromaticityError | ConjugationError | HybridizationError | AtropisomerError | StereoError | AdjustHsError | StereoOrderError | MappingValidationError | PropertyIntReadError | TemplateAttachmentOrderError"
    )]
    #[derive(Clone)]
    pub type BioConversionCause;
    #[wasm_bindgen(typescript_type = "BioStructureError")]
    pub type BioStructurePayload;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type BioTopologyPayload;
    #[wasm_bindgen(typescript_type = "PdbPostprocessError")]
    pub type BioChemistryPayload;
    #[wasm_bindgen(typescript_type = "SanitizeError")]
    pub type BioSanitizePayload;
    #[wasm_bindgen(typescript_type = "HydrogenError")]
    pub type BioHydrogensPayload;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type BioValencePayload;
    #[wasm_bindgen(typescript_type = "StereoError")]
    pub type BioStereoPayload;
}
fn conversion_source(
    owner: &ReadonlyValue<ck::BioMoleculeConversionError>,
) -> Option<BioConversionCause> {
    match owner.as_ref() {
        ck::BioMoleculeConversionError::Structure(_)
        | ck::BioMoleculeConversionError::Topology(_) => None,
        ck::BioMoleculeConversionError::Chemistry(_) => post_source(&owner.child(|e| match e {
            ck::BioMoleculeConversionError::Chemistry(v) => v,
            _ => unreachable!("immutable chemistry source"),
        }))
        .map(JsCast::unchecked_into),
        ck::BioMoleculeConversionError::Sanitize(v) => {
            crate::sanitize_errors::cause(v).map(JsCast::unchecked_into)
        }
        ck::BioMoleculeConversionError::Hydrogens(v) => {
            crate::hydrogen_errors::cause(v).map(JsCast::unchecked_into)
        }
        ck::BioMoleculeConversionError::Valence(v) => {
            crate::valence_errors::cause(v).map(JsCast::unchecked_into)
        }
        ck::BioMoleculeConversionError::Stereo(v) => {
            crate::stereo_errors::cause(v).map(JsCast::unchecked_into)
        }
    }
}
enum_error!(BioMoleculeConversionError,BioMoleculeConversionErrorKind,ck::BioMoleculeConversionError,{Structure=>ck::BioMoleculeConversionError::Structure(_),Topology=>ck::BioMoleculeConversionError::Topology(_),Chemistry=>ck::BioMoleculeConversionError::Chemistry(_),Sanitize=>ck::BioMoleculeConversionError::Sanitize(_),Hydrogens=>ck::BioMoleculeConversionError::Hydrogens(_),Valence=>ck::BioMoleculeConversionError::Valence(_),Stereo=>ck::BioMoleculeConversionError::Stereo(_)},{
#[wasm_bindgen(getter,js_name=structure)]pub fn structure(&self)->Option<BioStructurePayload>{match self.inner.as_ref(){ck::BioMoleculeConversionError::Structure(_)=>Some(crate::bio_read_errors::BioStructureError::from_owner(self.inner.child(|e|match e{ck::BioMoleculeConversionError::Structure(v)=>v,_=>unreachable!("immutable BIO conversion variant")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=topology)]pub fn topology(&self)->Option<BioTopologyPayload>{match self.inner.as_ref(){ck::BioMoleculeConversionError::Topology(_)=>Some(crate::model_errors::TopologyValidationError::from_owner(self.inner.child(|e|match e{ck::BioMoleculeConversionError::Topology(v)=>v,_=>unreachable!("immutable BIO conversion variant")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=chemistry)]pub fn chemistry(&self)->Option<BioChemistryPayload>{match self.inner.as_ref(){ck::BioMoleculeConversionError::Chemistry(_)=>Some(crate::bio_conversion_errors::PdbPostprocessError::from_owner(self.inner.child(|e|match e{ck::BioMoleculeConversionError::Chemistry(v)=>v,_=>unreachable!("immutable BIO conversion variant")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=sanitize)]pub fn sanitize(&self)->Option<BioSanitizePayload>{match self.inner.as_ref(){ck::BioMoleculeConversionError::Sanitize(_)=>Some(crate::sanitize_errors::SanitizeError::from_owner(self.inner.child(|e|match e{ck::BioMoleculeConversionError::Sanitize(v)=>v,_=>unreachable!("immutable BIO conversion variant")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=hydrogens)]pub fn hydrogens(&self)->Option<BioHydrogensPayload>{match self.inner.as_ref(){ck::BioMoleculeConversionError::Hydrogens(_)=>Some(crate::hydrogen_errors::HydrogenError::from_owner(self.inner.child(|e|match e{ck::BioMoleculeConversionError::Hydrogens(v)=>v,_=>unreachable!("immutable BIO conversion variant")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=valence)]pub fn valence(&self)->Option<BioValencePayload>{match self.inner.as_ref(){ck::BioMoleculeConversionError::Valence(_)=>Some(crate::valence_errors::ValenceError::from_owner(self.inner.child(|e|match e{ck::BioMoleculeConversionError::Valence(v)=>v,_=>unreachable!("immutable BIO conversion variant")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=stereo)]pub fn stereo(&self)->Option<BioStereoPayload>{match self.inner.as_ref(){ck::BioMoleculeConversionError::Stereo(_)=>Some(crate::stereo_errors::StereoError::from_owner(self.inner.child(|e|match e{ck::BioMoleculeConversionError::Stereo(v)=>v,_=>unreachable!("immutable BIO conversion variant")})).unchecked_into()),_=>None}}
},source=BioConversionCause,cause=conversion_source);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "BioMoleculeConversionError | OperationError")]
    #[derive(Clone)]
    pub type BioMoleculeCause;
}
enum_error!(BioMoleculeError,BioMoleculeErrorKind,ck::BioMoleculeError,{Conversion=>ck::BioMoleculeError::Conversion(_),Construction=>ck::BioMoleculeError::Construction(_)},{},source=BioMoleculeCause,cause=|owner:&ReadonlyValue<ck::BioMoleculeError>|Some(match owner.as_ref(){ck::BioMoleculeError::Conversion(_)=>BioMoleculeConversionError::from_owner(owner.child(|e|match e{ck::BioMoleculeError::Conversion(v)=>v,_=>unreachable!("immutable conversion cause")})),ck::BioMoleculeError::Construction(_)=>crate::operation_errors::OperationError::from_owner(owner.child(|e|match e{ck::BioMoleculeError::Construction(v)=>v,_=>unreachable!("immutable construction cause")}))}.unchecked_into()));
