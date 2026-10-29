//! MMFF table/property errors preserve native transparent causes and payloads.
use crate::{
    errors::child,
    errors::{enum_error, error_projection, exception},
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
crate::enums::enum_value!(MmffParamTable,"MmffParamTable"{Def="def",Prop="prop",Pbci="pbci",Chg="chg",Bond="bond",Bndk="bndk",Stbn="stbn",Dfsb="dfsb",Angle="angle",Oop="oop",Tor="tor",Vdw="vdw",HerschbachLaurie="herschbach_laurie",CovRadPauEle="cov_rad_pau_ele"});
#[wasm_bindgen]
pub enum MmffParamParseCauseKind {
    InvalidUnsigned,
    InvalidFloat,
    MissingToken,
    EmptyProcessedLine,
}
#[wasm_bindgen]
pub struct MmffParamParseCause(ck::MmffParamParseCause);
#[wasm_bindgen]
impl MmffParamParseCause {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> MmffParamParseCauseKind {
        match &self.0 {
            ck::MmffParamParseCause::InvalidUnsigned { .. } => {
                MmffParamParseCauseKind::InvalidUnsigned
            }
            ck::MmffParamParseCause::InvalidFloat { .. } => MmffParamParseCauseKind::InvalidFloat,
            ck::MmffParamParseCause::MissingToken => MmffParamParseCauseKind::MissingToken,
            ck::MmffParamParseCause::EmptyProcessedLine => {
                MmffParamParseCauseKind::EmptyProcessedLine
            }
        }
    }
    #[wasm_bindgen(getter)]
    pub fn cell(&self) -> Option<String> {
        match &self.0 {
            ck::MmffParamParseCause::InvalidUnsigned { cell }
            | ck::MmffParamParseCause::InvalidFloat { cell } => Some(cell.clone()),
            _ => None,
        }
    }
}
error_projection!(MmffParamParseError, ck::MmffParamParseError, {
    #[wasm_bindgen(getter)]
    pub fn table(&self) -> MmffParamTable {
        self.inner.table().into()
    }
    #[wasm_bindgen(getter)]
    pub fn line(&self) -> usize {
        self.inner.line()
    }
    #[wasm_bindgen(getter)]
    pub fn column(&self) -> usize {
        self.inner.column()
    }
    #[wasm_bindgen(getter)]
    pub fn cause(&self) -> MmffParamParseCause {
        MmffParamParseCause(self.inner.cause().clone())
    }
});
enum_error!(MmffAngleLookupError,MmffAngleLookupErrorKind,ck::MmffAngleLookupError,{MissingDefinition=>ck::MmffAngleLookupError::MissingDefinition{..}},{#[wasm_bindgen(getter,js_name=atomType)]pub fn atom_type(&self)->u32{let ck::MmffAngleLookupError::MissingDefinition{atom_type}=*self.inner;atom_type}});
enum_error!(MmffOopLookupError,MmffOopLookupErrorKind,ck::MmffOopLookupError,{MissingDefinition=>ck::MmffOopLookupError::MissingDefinition{..}},{#[wasm_bindgen(getter,js_name=atomType)]pub fn atom_type(&self)->u32{let ck::MmffOopLookupError::MissingDefinition{atom_type}=*self.inner;atom_type}});
enum_error!(MmffTorLookupError,MmffTorLookupErrorKind,ck::MmffTorLookupError,{MissingDefinition=>ck::MmffTorLookupError::MissingDefinition{..},InvalidEquivalentLevel=>ck::MmffTorLookupError::InvalidEquivalentLevel{..}},{
#[wasm_bindgen(getter,js_name=atomType)]pub fn atom_type(&self)->Option<u32>{match *self.inner{ck::MmffTorLookupError::MissingDefinition{atom_type}=>Some(atom_type),_=>None}}
#[wasm_bindgen(getter)]pub fn level(&self)->Option<usize>{match *self.inner{ck::MmffTorLookupError::InvalidEquivalentLevel{level}=>Some(level),_=>None}}
});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "MmffParamParseError")]
    #[derive(Clone)]
    pub type MmffParseCause;
    #[wasm_bindgen(typescript_type = "MmffVdwParseError")]
    pub type VdwPayload;
    #[wasm_bindgen(typescript_type = "MmffAngleLookupError")]
    pub type AnglePayload;
    #[wasm_bindgen(typescript_type = "MmffTorLookupError")]
    pub type TorsionPayload;
    #[wasm_bindgen(typescript_type = "MmffOopLookupError")]
    pub type OopPayload;
    #[wasm_bindgen(typescript_type = "MmffParamError")]
    pub type ParamsPayload;
    #[wasm_bindgen(typescript_type = "KekulizeError")]
    pub type KekulizePayload;
    #[wasm_bindgen(typescript_type = "AromaticityError")]
    pub type AromaticityPayload;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    pub type RingsPayload;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type ValencePayload;
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | MappingValidationError | MoleculePropertyError | StereoGroupError | TemplateAttachmentOrderError | MmffParamParseError"
    )]
    #[derive(Clone)]
    pub type PropertiesCause;
}
enum_error!(MmffVdwParseError,MmffVdwParseErrorKind,ck::MmffVdwParseError,{Table=>ck::MmffVdwParseError::Table(_),MissingHeader=>ck::MmffVdwParseError::MissingHeader},{},source=MmffParseCause,cause=|owner:&ReadonlyValue<ck::MmffVdwParseError>|match owner.as_ref(){ck::MmffVdwParseError::Table(_)=>Some(child!(owner,ck::MmffVdwParseError::Table,MmffParamParseError)),_=>None});
pub(crate) fn param_cause(owner: &ReadonlyValue<ck::MmffParamError>) -> Option<MmffParseCause> {
    match owner.as_ref() {
        ck::MmffParamError::Vdw(ck::MmffVdwParseError::Table(_)) => Some(
            MmffParamParseError::from_owner(owner.child(|e| match e {
                ck::MmffParamError::Vdw(ck::MmffVdwParseError::Table(v)) => v,
                _ => unreachable!("immutable MMFF table cause"),
            }))
            .unchecked_into(),
        ),
        _ => None,
    }
}
enum_error!(MmffParamError,MmffParamErrorKind,ck::MmffParamError,{Table=>ck::MmffParamError::Table(_),Vdw=>ck::MmffParamError::Vdw(_),Angle=>ck::MmffParamError::Angle(_),Torsion=>ck::MmffParamError::Torsion(_),Oop=>ck::MmffParamError::Oop(_),MissingDefaultData=>ck::MmffParamError::MissingDefaultData{..}},{
#[wasm_bindgen(getter)]pub fn symbol(&self)->Option<String>{match self.inner.as_ref(){ck::MmffParamError::MissingDefaultData{symbol}=>Some((*symbol).into()),_=>None}}
#[wasm_bindgen(getter)]pub fn table(&self)->Option<MmffParseCause>{match self.inner.as_ref(){ck::MmffParamError::Table(_)=>Some(child!(&self.inner,ck::MmffParamError::Table,MmffParamParseError)),_=>None}}
#[wasm_bindgen(getter)]pub fn vdw(&self)->Option<VdwPayload>{match self.inner.as_ref(){ck::MmffParamError::Vdw(_)=>Some(child!(&self.inner,ck::MmffParamError::Vdw,MmffVdwParseError)),_=>None}}
#[wasm_bindgen(getter)]pub fn angle(&self)->Option<AnglePayload>{match self.inner.as_ref(){ck::MmffParamError::Angle(_)=>Some(child!(&self.inner,ck::MmffParamError::Angle,MmffAngleLookupError)),_=>None}}
#[wasm_bindgen(getter)]pub fn torsion(&self)->Option<TorsionPayload>{match self.inner.as_ref(){ck::MmffParamError::Torsion(_)=>Some(child!(&self.inner,ck::MmffParamError::Torsion,MmffTorLookupError)),_=>None}}
#[wasm_bindgen(getter)]pub fn oop(&self)->Option<OopPayload>{match self.inner.as_ref(){ck::MmffParamError::Oop(_)=>Some(child!(&self.inner,ck::MmffParamError::Oop,MmffOopLookupError)),_=>None}}
},source=MmffParseCause,cause=param_cause);
pub(crate) fn properties_cause(
    owner: &ReadonlyValue<ck::MmffMolPropertiesError>,
) -> Option<PropertiesCause> {
    match owner.as_ref() {
        ck::MmffMolPropertiesError::Params(_) => param_cause(&owner.child(|e| match e {
            ck::MmffMolPropertiesError::Params(v) => v,
            _ => unreachable!("immutable MMFF parameter error"),
        }))
        .map(JsCast::unchecked_into),
        ck::MmffMolPropertiesError::Kekulize(v) => {
            crate::kekulize_errors::cause(v).map(JsCast::unchecked_into)
        }
        ck::MmffMolPropertiesError::Aromaticity(v) => {
            crate::aromaticity_errors::cause(v).map(JsCast::unchecked_into)
        }
        ck::MmffMolPropertiesError::RingFinding(v) => {
            crate::ring_errors::cause(v).map(JsCast::unchecked_into)
        }
        ck::MmffMolPropertiesError::Valence(v) => {
            crate::valence_errors::cause(v).map(JsCast::unchecked_into)
        }
        _ => None,
    }
}
enum_error!(MmffMolPropertiesError,MmffMolPropertiesErrorKind,ck::MmffMolPropertiesError,{Params=>ck::MmffMolPropertiesError::Params(_),Kekulize=>ck::MmffMolPropertiesError::Kekulize(_),Aromaticity=>ck::MmffMolPropertiesError::Aromaticity(_),RingFinding=>ck::MmffMolPropertiesError::RingFinding(_),Valence=>ck::MmffMolPropertiesError::Valence(_),AtomIndexOutOfRange=>ck::MmffMolPropertiesError::AtomIndexOutOfRange{..},AtomTypePropertiesMissing=>ck::MmffMolPropertiesError::AtomTypePropertiesMissing{..},AtomTypePbciMissing=>ck::MmffMolPropertiesError::AtomTypePbciMissing{..}},{
#[wasm_bindgen(getter,js_name=atomIndex)]pub fn atom_index(&self)->Option<usize>{match self.inner.as_ref(){ck::MmffMolPropertiesError::AtomIndexOutOfRange{atom_index,..}=>Some(*atom_index),_=>None}}
#[wasm_bindgen(getter)]pub fn atoms(&self)->Option<usize>{match self.inner.as_ref(){ck::MmffMolPropertiesError::AtomIndexOutOfRange{atoms,..}=>Some(*atoms),_=>None}}
#[wasm_bindgen(getter,js_name=atomType)]pub fn atom_type(&self)->Option<u8>{match self.inner.as_ref(){ck::MmffMolPropertiesError::AtomTypePropertiesMissing{atom_type}|ck::MmffMolPropertiesError::AtomTypePbciMissing{atom_type}=>Some(*atom_type),_=>None}}
#[wasm_bindgen(getter)]pub fn params(&self)->Option<ParamsPayload>{match self.inner.as_ref(){ck::MmffMolPropertiesError::Params(_)=>Some(child!(&self.inner,ck::MmffMolPropertiesError::Params,MmffParamError)),_=>None}}
#[wasm_bindgen(getter)]pub fn kekulize(&self)->Option<KekulizePayload>{match self.inner.as_ref(){ck::MmffMolPropertiesError::Kekulize(_)=>Some(child!(&self.inner,ck::MmffMolPropertiesError::Kekulize,crate::kekulize_errors::KekulizeError)),_=>None}}
#[wasm_bindgen(getter)]pub fn aromaticity(&self)->Option<AromaticityPayload>{match self.inner.as_ref(){ck::MmffMolPropertiesError::Aromaticity(_)=>Some(child!(&self.inner,ck::MmffMolPropertiesError::Aromaticity,crate::aromaticity_errors::AromaticityError)),_=>None}}
#[wasm_bindgen(getter,js_name=ringFinding)]pub fn ring_finding(&self)->Option<RingsPayload>{match self.inner.as_ref(){ck::MmffMolPropertiesError::RingFinding(_)=>Some(child!(&self.inner,ck::MmffMolPropertiesError::RingFinding,crate::ring_errors::RingFindingError)),_=>None}}
#[wasm_bindgen(getter)]pub fn valence(&self)->Option<ValencePayload>{match self.inner.as_ref(){ck::MmffMolPropertiesError::Valence(_)=>Some(child!(&self.inner,ck::MmffMolPropertiesError::Valence,crate::valence_errors::ValenceError)),_=>None}}
},source=PropertiesCause,cause=properties_cause);
