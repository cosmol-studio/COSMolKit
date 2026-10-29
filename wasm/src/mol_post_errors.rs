//! Source-ordered MOL finalization error values, without interpretation.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "ValenceError | PropertyIntReadError | PropertyStringError | RingFindingError | StereoError | DoubleBondStereoError | BondDirectionStereoError | AtropisomerError | SanitizeError | HydrogenError | LegacyStereoError | AttachmentExpansionError | QueryAtomConversionError | QueryGraphError | QueryStateError | TopologyValidationError | TopologyEditError | PropertyValueError | AtomPropertyError | BondValueError | CoordinateValidationError | MoleculePropertyError"
    )]
    #[derive(Clone)]
    pub type MolProcessingCause;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type ProcessingValence;
    #[wasm_bindgen(typescript_type = "PropertyIntReadError")]
    pub type ProcessingPropertyInt;
    #[wasm_bindgen(typescript_type = "PropertyStringError")]
    pub type ProcessingPropertyString;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    pub type ProcessingRings;
    #[wasm_bindgen(typescript_type = "StereoError")]
    pub type ProcessingStereo;
    #[wasm_bindgen(typescript_type = "DoubleBondStereoError")]
    pub type ProcessingDoubleBondStereo;
    #[wasm_bindgen(typescript_type = "BondDirectionStereoError")]
    pub type ProcessingBondDirectionStereo;
    #[wasm_bindgen(typescript_type = "AtropisomerError")]
    pub type ProcessingAtropisomer;
    #[wasm_bindgen(typescript_type = "SanitizeError")]
    pub type ProcessingSanitize;
    #[wasm_bindgen(typescript_type = "HydrogenError")]
    pub type ProcessingHydrogen;
    #[wasm_bindgen(typescript_type = "LegacyStereoError")]
    pub type ProcessingLegacyStereo;
    #[wasm_bindgen(typescript_type = "AttachmentExpansionError")]
    pub type ProcessingAttachmentExpansion;
    #[wasm_bindgen(typescript_type = "QueryAtomConversionError")]
    pub type ProcessingQueryAtom;
    #[wasm_bindgen(typescript_type = "QueryGraphError")]
    pub type ProcessingQueryGraph;
    #[wasm_bindgen(typescript_type = "QueryStateError")]
    pub type ProcessingQueryState;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type ProcessingTopology;
    #[wasm_bindgen(typescript_type = "TopologyEditError")]
    pub type ProcessingTopologyEdit;
    #[wasm_bindgen(typescript_type = "PropertyValueError")]
    pub type ProcessingPropertyValue;
    #[wasm_bindgen(typescript_type = "AtomPropertyError")]
    pub type ProcessingAtomProperty;
    #[wasm_bindgen(typescript_type = "BondValueError")]
    pub type ProcessingBondValue;
    #[wasm_bindgen(typescript_type = "CoordinateValidationError")]
    pub type ProcessingCoordinates;
    #[wasm_bindgen(typescript_type = "MoleculePropertyError")]
    pub type ProcessingMoleculeProperty;
    #[wasm_bindgen(typescript_type = "MolProcessingError")]
    #[derive(Clone)]
    pub type MolPostCause;
    #[wasm_bindgen(typescript_type = "MolProcessingError")]
    pub type PostProcessing;
    #[wasm_bindgen(typescript_type = "QueryAtomConversionError")]
    pub type PostQueryAtom;
    #[wasm_bindgen(typescript_type = "number | string")]
    pub type PostValue;
}
fn processing_cause(e: &ck::MolProcessingError) -> Option<MolProcessingCause> {
    Some(match e {
        ck::MolProcessingError::Valence(v) => {
            crate::valence_errors::ValenceError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::PropertyInt(v) => {
            crate::property_errors::PropertyIntReadError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::PropertyString(v) => {
            crate::property_errors::PropertyStringError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::Rings(v) => {
            crate::ring_errors::RingFindingError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::Stereo(v) => {
            crate::stereo_errors::StereoError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::DoubleBondStereo(v) => {
            crate::double_stereo_errors::DoubleBondStereoError::from_rust(v.clone())
                .unchecked_into()
        }
        ck::MolProcessingError::BondDirectionStereo(v) => {
            crate::mol_post_support_errors::BondDirectionStereoError::from_rust(v.clone())
                .unchecked_into()
        }
        ck::MolProcessingError::Atropisomer(v) => {
            crate::atropisomer_errors::AtropisomerError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::Sanitize(v) => {
            crate::sanitize_errors::SanitizeError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::Hydrogen(v) => {
            crate::hydrogen_errors::HydrogenError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::LegacyStereo(v) => {
            crate::legacy_stereo_errors::LegacyStereoError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::AttachmentExpansion(v) => {
            crate::mol_post_support_errors::AttachmentExpansionError::from_rust(v.clone())
                .unchecked_into()
        }
        ck::MolProcessingError::QueryAtom(v) => {
            crate::query_model_errors::QueryAtomConversionError::from_rust(v.clone())
                .unchecked_into()
        }
        ck::MolProcessingError::QueryGraph(v) => {
            crate::query_model_errors::QueryGraphError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::QueryState(v) => {
            crate::query_errors::QueryStateError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::Topology(v) => {
            crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::TopologyEdit(v) => {
            crate::topology_errors::TopologyEditError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::PropertyValue(v) => {
            crate::model_errors::PropertyValueError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::AtomProperty(v) => {
            crate::model_errors::AtomPropertyError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::BondValue(v) => {
            crate::model_errors::BondValueError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::Coordinates(v) => {
            crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into()
        }
        ck::MolProcessingError::MoleculeProperty(v) => {
            crate::model_errors::MoleculePropertyError::from_rust(v.clone()).unchecked_into()
        }
    })
}
enum_error!(MolProcessingError,MolProcessingErrorKind,ck::MolProcessingError,{
Valence=>ck::MolProcessingError::Valence(_),
PropertyInt=>ck::MolProcessingError::PropertyInt(_),
PropertyString=>ck::MolProcessingError::PropertyString(_),
Rings=>ck::MolProcessingError::Rings(_),
Stereo=>ck::MolProcessingError::Stereo(_),
DoubleBondStereo=>ck::MolProcessingError::DoubleBondStereo(_),
BondDirectionStereo=>ck::MolProcessingError::BondDirectionStereo(_),
Atropisomer=>ck::MolProcessingError::Atropisomer(_),
Sanitize=>ck::MolProcessingError::Sanitize(_),
Hydrogen=>ck::MolProcessingError::Hydrogen(_),
LegacyStereo=>ck::MolProcessingError::LegacyStereo(_),
AttachmentExpansion=>ck::MolProcessingError::AttachmentExpansion(_),
QueryAtom=>ck::MolProcessingError::QueryAtom(_),
QueryGraph=>ck::MolProcessingError::QueryGraph(_),
QueryState=>ck::MolProcessingError::QueryState(_),
Topology=>ck::MolProcessingError::Topology(_),
TopologyEdit=>ck::MolProcessingError::TopologyEdit(_),
PropertyValue=>ck::MolProcessingError::PropertyValue(_),
AtomProperty=>ck::MolProcessingError::AtomProperty(_),
BondValue=>ck::MolProcessingError::BondValue(_),
Coordinates=>ck::MolProcessingError::Coordinates(_),
MoleculeProperty=>ck::MolProcessingError::MoleculeProperty(_),
},{
#[wasm_bindgen(getter,js_name=valence)] pub fn valence(&self)->Option<ProcessingValence>{match self.inner.as_ref(){ck::MolProcessingError::Valence(v)=>Some(crate::valence_errors::ValenceError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=propertyInt)] pub fn property_int(&self)->Option<ProcessingPropertyInt>{match self.inner.as_ref(){ck::MolProcessingError::PropertyInt(v)=>Some(crate::property_errors::PropertyIntReadError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=propertyString)] pub fn property_string(&self)->Option<ProcessingPropertyString>{match self.inner.as_ref(){ck::MolProcessingError::PropertyString(v)=>Some(crate::property_errors::PropertyStringError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=rings)] pub fn rings(&self)->Option<ProcessingRings>{match self.inner.as_ref(){ck::MolProcessingError::Rings(v)=>Some(crate::ring_errors::RingFindingError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=stereo)] pub fn stereo(&self)->Option<ProcessingStereo>{match self.inner.as_ref(){ck::MolProcessingError::Stereo(v)=>Some(crate::stereo_errors::StereoError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=doubleBondStereo)] pub fn double_bond_stereo(&self)->Option<ProcessingDoubleBondStereo>{match self.inner.as_ref(){ck::MolProcessingError::DoubleBondStereo(v)=>Some(crate::double_stereo_errors::DoubleBondStereoError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=bondDirectionStereo)] pub fn bond_direction_stereo(&self)->Option<ProcessingBondDirectionStereo>{match self.inner.as_ref(){ck::MolProcessingError::BondDirectionStereo(v)=>Some(crate::mol_post_support_errors::BondDirectionStereoError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=atropisomer)] pub fn atropisomer(&self)->Option<ProcessingAtropisomer>{match self.inner.as_ref(){ck::MolProcessingError::Atropisomer(v)=>Some(crate::atropisomer_errors::AtropisomerError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=sanitize)] pub fn sanitize(&self)->Option<ProcessingSanitize>{match self.inner.as_ref(){ck::MolProcessingError::Sanitize(v)=>Some(crate::sanitize_errors::SanitizeError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=hydrogen)] pub fn hydrogen(&self)->Option<ProcessingHydrogen>{match self.inner.as_ref(){ck::MolProcessingError::Hydrogen(v)=>Some(crate::hydrogen_errors::HydrogenError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=legacyStereo)] pub fn legacy_stereo(&self)->Option<ProcessingLegacyStereo>{match self.inner.as_ref(){ck::MolProcessingError::LegacyStereo(v)=>Some(crate::legacy_stereo_errors::LegacyStereoError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=attachmentExpansion)] pub fn attachment_expansion(&self)->Option<ProcessingAttachmentExpansion>{match self.inner.as_ref(){ck::MolProcessingError::AttachmentExpansion(v)=>Some(crate::mol_post_support_errors::AttachmentExpansionError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=queryAtom)] pub fn query_atom(&self)->Option<ProcessingQueryAtom>{match self.inner.as_ref(){ck::MolProcessingError::QueryAtom(v)=>Some(crate::query_model_errors::QueryAtomConversionError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=queryGraph)] pub fn query_graph(&self)->Option<ProcessingQueryGraph>{match self.inner.as_ref(){ck::MolProcessingError::QueryGraph(v)=>Some(crate::query_model_errors::QueryGraphError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=queryState)] pub fn query_state(&self)->Option<ProcessingQueryState>{match self.inner.as_ref(){ck::MolProcessingError::QueryState(v)=>Some(crate::query_errors::QueryStateError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=topology)] pub fn topology(&self)->Option<ProcessingTopology>{match self.inner.as_ref(){ck::MolProcessingError::Topology(v)=>Some(crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=topologyEdit)] pub fn topology_edit(&self)->Option<ProcessingTopologyEdit>{match self.inner.as_ref(){ck::MolProcessingError::TopologyEdit(v)=>Some(crate::topology_errors::TopologyEditError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=propertyValue)] pub fn property_value(&self)->Option<ProcessingPropertyValue>{match self.inner.as_ref(){ck::MolProcessingError::PropertyValue(v)=>Some(crate::model_errors::PropertyValueError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=atomProperty)] pub fn atom_property(&self)->Option<ProcessingAtomProperty>{match self.inner.as_ref(){ck::MolProcessingError::AtomProperty(v)=>Some(crate::model_errors::AtomPropertyError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=bondValue)] pub fn bond_value(&self)->Option<ProcessingBondValue>{match self.inner.as_ref(){ck::MolProcessingError::BondValue(v)=>Some(crate::model_errors::BondValueError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=coordinates)] pub fn coordinates(&self)->Option<ProcessingCoordinates>{match self.inner.as_ref(){ck::MolProcessingError::Coordinates(v)=>Some(crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=moleculeProperty)] pub fn molecule_property(&self)->Option<ProcessingMoleculeProperty>{match self.inner.as_ref(){ck::MolProcessingError::MoleculeProperty(v)=>Some(crate::model_errors::MoleculePropertyError::from_rust(v.clone()).unchecked_into()),_=>None}}
},source=MolProcessingCause,cause=processing_cause);
fn post_cause(e: &ck::MolPostError) -> Option<MolPostCause> {
    match e {
        ck::MolPostError::AttachmentExpansion(v) | ck::MolPostError::Processing(v) => {
            Some(MolProcessingError::from_rust(v.clone()).unchecked_into())
        }
        _ => None,
    }
}
enum_error!(MolPostError,MolPostErrorKind,ck::MolPostError,{
MissingCapability=>ck::MolPostError::MissingCapability(_),
Representation=>ck::MolPostError::Representation(_),
AttachmentExpansion=>ck::MolPostError::AttachmentExpansion(_),
QueryAtomConversion=>ck::MolPostError::QueryAtomConversion(_),
Processing=>ck::MolPostError::Processing(_),
UnsignedPropertyOverflow=>ck::MolPostError::UnsignedPropertyOverflow{..},
InvalidPropertyKind=>ck::MolPostError::InvalidPropertyKind{..},
AttachmentValue=>ck::MolPostError::AttachmentValue{..},
SGroupPropertyKind=>ck::MolPostError::SGroupPropertyKind{..},
DataFieldNumber=>ck::MolPostError::DataFieldNumber{..},
},{
#[wasm_bindgen(getter)] pub fn atom(&self)->Option<usize>{match self.inner.as_ref(){ck::MolPostError::UnsignedPropertyOverflow{atom,..}|ck::MolPostError::InvalidPropertyKind{atom,..}|ck::MolPostError::AttachmentValue{atom,..}=>Some(atom.index()),_=>None}}
#[wasm_bindgen(getter)] pub fn property(&self)->Option<String>{match self.inner.as_ref(){ck::MolPostError::UnsignedPropertyOverflow{property,..}|ck::MolPostError::InvalidPropertyKind{property,..}|ck::MolPostError::SGroupPropertyKind{property,..}=>Some(property.to_string()),_=>None}}
#[wasm_bindgen(getter)] pub fn kind(&self)->Option<crate::model_errors::PropertyValueKind>{match self.inner.as_ref(){ck::MolPostError::InvalidPropertyKind{kind,..}|ck::MolPostError::SGroupPropertyKind{kind,..}=>Some((*kind).into()),_=>None}}
#[wasm_bindgen(getter)] pub fn value(&self)->Result<Option<PostValue>,JsValue>{match self.inner.as_ref(){ck::MolPostError::UnsignedPropertyOverflow{value,..}=>Ok(Some(JsValue::from(*value).unchecked_into())),ck::MolPostError::AttachmentValue{value,..}|ck::MolPostError::DataFieldNumber{value,..}=>crate::property_errors::text(value).map(|v|Some(JsValue::from(v).unchecked_into())),_=>Ok(None)}}
#[wasm_bindgen(getter)] pub fn group(&self)->Option<usize>{match self.inner.as_ref(){ck::MolPostError::SGroupPropertyKind{group,..}=>Some(group.index()),_=>None}}
#[wasm_bindgen(getter)] pub fn field(&self)->Option<String>{match self.inner.as_ref(){ck::MolPostError::DataFieldNumber{field,..}=>Some(field.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=missingCapability)] pub fn missing_capability(&self)->Option<String>{match self.inner.as_ref(){ck::MolPostError::MissingCapability(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=representation)] pub fn representation(&self)->Option<String>{match self.inner.as_ref(){ck::MolPostError::Representation(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=attachmentExpansion)] pub fn attachment_expansion(&self)->Option<PostProcessing>{match self.inner.as_ref(){ck::MolPostError::AttachmentExpansion(v)=>Some(MolProcessingError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=processing)] pub fn processing(&self)->Option<PostProcessing>{match self.inner.as_ref(){ck::MolPostError::Processing(v)=>Some(MolProcessingError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=queryAtomConversion)] pub fn query_atom_conversion(&self)->Option<PostQueryAtom>{match self.inner.as_ref(){ck::MolPostError::QueryAtomConversion(v)=>Some(crate::query_model_errors::QueryAtomConversionError::from_rust(*v).unchecked_into()),_=>None}}
},source=MolPostCause,cause=post_cause);
