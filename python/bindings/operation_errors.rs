//! Complete operation errors and readonly context, without runtime authority.
use crate::alignment_errors::AlignmentError;
use crate::cip_atom_code_errors::AtomCodeError;
use crate::cip_atom_code_errors::CipLabelerError;
use crate::conformer_execution::ConformerRunError;
use crate::conformer_execution::EnumerationRunError;
use crate::coordinate_2d_errors::Coordinate2DError;
use crate::coordinate_input_errors::CoordinateInputError;
use crate::error_framework::{Source, error_projection};
use crate::forcefield_execution_errors::MmffOptimizationError;
use crate::forcefield_execution_errors::UffOptimizationError;
use crate::fragment_errors::MoleculeFragmentsError;
use crate::hydrogen_errors::HydrogenError;
use crate::kekulize_errors::KekulizeError;
use crate::model_edit_errors::AtomPropertyError;
use crate::model_edit_errors::BondValueError;
use crate::model_edit_errors::MoleculePropertyError;
use crate::model_errors::CoordinateValidationError;
use crate::model_errors::TopologyValidationError;
use crate::operation_metadata::MoleculeOpSpec;
use crate::potential_stereo_errors::PotentialStereoError;
use crate::reaction_execution_errors::ReactionApplyError;
use crate::reaction_execution_errors::ReactionRunError;
use crate::retained_error::RetainedError;
use crate::ring_errors::RingFindingError;
use crate::sanitize_errors::SanitizeError;
use crate::sanitize_leaf_errors::RadicalError;
use crate::sanitize_support_errors::AromaticityError;
use crate::scaffold_errors::ScaffoldError;
use crate::smiles_parse_errors::SmilesParseError;
use crate::smiles_stereo_errors::SmilesStereoError;
use crate::stereo_errors::StereoError;
use crate::tautomer_errors::TautomerRunError;
use crate::topology_edit_errors::TopologyEditError;
use crate::transform_errors::TransformError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;
type OperationSource=Source<crate::error_framework::ErrorObject<MoleculeFragmentsError>,Source<crate::error_framework::ErrorObject<ScaffoldError>,Source<crate::error_framework::ErrorObject<ReactionRunError>,Source<crate::error_framework::ErrorObject<ReactionApplyError>,Source<crate::error_framework::ErrorObject<EnumerationRunError>,Source<crate::error_framework::ErrorObject<AlignmentError>,Source<crate::error_framework::ErrorObject<ConformerRunError>,Source<crate::error_framework::ErrorObject<TautomerRunError>,Source<crate::error_framework::ErrorObject<TopologyValidationError>,Source<crate::error_framework::ErrorObject<TopologyEditError>,Source<crate::error_framework::ErrorObject<CoordinateValidationError>,Source<crate::error_framework::ErrorObject<MoleculePropertyError>,Source<crate::error_framework::ErrorObject<AtomPropertyError>,Source<crate::error_framework::ErrorObject<BondValueError>,Source<crate::error_framework::ErrorObject<ValenceError>,Source<crate::error_framework::ErrorObject<RadicalError>,Source<crate::error_framework::ErrorObject<RingFindingError>,Source<crate::error_framework::ErrorObject<PotentialStereoError>,Source<crate::error_framework::ErrorObject<StereoError>,Source<crate::error_framework::ErrorObject<CipLabelerError>,Source<crate::error_framework::ErrorObject<AtomCodeError>,Source<crate::error_framework::ErrorObject<TransformError>,Source<crate::error_framework::ErrorObject<CoordinateInputError>,Source<crate::error_framework::ErrorObject<Coordinate2DError>,Source<crate::error_framework::ErrorObject<UffOptimizationError>,Source<crate::error_framework::ErrorObject<MmffOptimizationError>,Source<crate::error_framework::ErrorObject<KekulizeError>,Source<crate::error_framework::ErrorObject<AromaticityError>,Source<crate::error_framework::ErrorObject<SanitizeError>,crate::error_framework::ErrorObject<HydrogenError>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>;
fn project_operation_source(inner: &RetainedError<ck::OperationError>) -> Option<PyErr> {
    let node = inner.node().source()?;
    Some(match inner.as_ref() {
        ck::OperationError::Fragments(value) => MoleculeFragmentsError::from_rust(value.clone()),
        ck::OperationError::Scaffold(value) => ScaffoldError::from_rust(value.clone()),
        ck::OperationError::ReactionRun(value) => ReactionRunError::from_rust(value.clone()),
        ck::OperationError::ReactionApply(value) => ReactionApplyError::from_rust(value.clone()),
        ck::OperationError::Enumeration(value) => EnumerationRunError::from_ref(node),
        ck::OperationError::Alignment(value) => AlignmentError::from_rust(value.clone()),
        ck::OperationError::Conformer(value) => ConformerRunError::from_ref(node),
        ck::OperationError::Tautomer(value) => TautomerRunError::from_rust(value.clone()),
        ck::OperationError::InvalidTopology(value) => {
            TopologyValidationError::from_rust(value.clone())
        }
        ck::OperationError::InvalidTopologyEdit(value) => {
            TopologyEditError::from_rust(value.clone())
        }
        ck::OperationError::InvalidCoordinates(value) => {
            CoordinateValidationError::from_rust(value.clone())
        }
        ck::OperationError::InvalidProperty(value) => {
            MoleculePropertyError::from_rust(value.clone())
        }
        ck::OperationError::AtomProperty(value) => AtomPropertyError::from_rust(value.clone()),
        ck::OperationError::BondProperty(value) => BondValueError::from_rust(value.clone()),
        ck::OperationError::Valence(value) => ValenceError::from_rust(value.clone()),
        ck::OperationError::Radical(value) => RadicalError::from_rust(value.clone()),
        ck::OperationError::Rings(value) => RingFindingError::from_rust(value.clone()),
        ck::OperationError::PotentialStereo(value) => {
            PotentialStereoError::from_rust(value.clone())
        }
        ck::OperationError::Stereo(value) => StereoError::from_rust(value.clone()),
        ck::OperationError::CipLabeler(value) => CipLabelerError::from_rust(value.clone()),
        ck::OperationError::AtomCode(value) => AtomCodeError::from_rust(value.clone()),
        ck::OperationError::Transform(value) => TransformError::from_rust(value.clone()),
        ck::OperationError::CoordinateInput(value) => {
            CoordinateInputError::from_rust(value.clone())
        }
        ck::OperationError::Coordinate2D(value) => Coordinate2DError::from_ref(node),
        ck::OperationError::UffOptimization(value) => UffOptimizationError::from_ref(node),
        ck::OperationError::MmffOptimization(value) => MmffOptimizationError::from_ref(node),
        ck::OperationError::Kekulize(value) => KekulizeError::from_rust(value.clone()),
        ck::OperationError::Aromaticity(value) => AromaticityError::from_rust(value.clone()),
        ck::OperationError::Sanitize(value) => SanitizeError::from_rust(value.clone()),
        ck::OperationError::Hydrogen(value) => HydrogenError::from_rust(value.clone()),
        _ => unreachable!("native OperationError::source and projection disagree"),
    })
}
use crate::chemical_values::enum_value;
enum_value!(MoleculeOpOutput, { Single, Multiple, LazyMultiple });
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct UnsupportedFeatureError {
    owner: RetainedError<ck::OperationError>,
}
#[pymethods]
impl UnsupportedFeatureError {
    #[getter]
    fn feature(&self) -> &'static str {
        let ck::OperationError::UnsupportedFeature { source, .. } = self.owner.as_ref() else {
            unreachable!()
        };
        source.feature
    }
    #[getter]
    fn reason(&self) -> &'static str {
        let ck::OperationError::UnsupportedFeature { source, .. } = self.owner.as_ref() else {
            unreachable!()
        };
        source.reason
    }
}

#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct OperationDerivedState {
    owner: RetainedError<ck::OperationError>,
}
#[pymethods]
impl OperationDerivedState {
    #[getter]
    fn bits(&self) -> u16 {
        match self.owner.as_ref() {
            ck::OperationError::DerivedEffectContract { states, .. }
            | ck::OperationError::DerivedPreservationMismatch { states, .. } => {
                states.bits() as u16
            }
            _ => unreachable!(),
        }
    }
    fn __repr__(&self) -> String {
        match self.owner.as_ref() {
            ck::OperationError::DerivedEffectContract { states, .. }
            | ck::OperationError::DerivedPreservationMismatch { states, .. } => {
                format!("{:?}", states)
            }
            _ => unreachable!(),
        }
    }
    fn __eq__(&self, other: &Self) -> bool {
        match (self.owner.as_ref(), other.owner.as_ref()) {
            (
                ck::OperationError::DerivedEffectContract { states: a, .. },
                ck::OperationError::DerivedEffectContract { states: b, .. },
            ) => a == b,
            (
                ck::OperationError::DerivedEffectContract { states: a, .. },
                ck::OperationError::DerivedPreservationMismatch { states: b, .. },
            ) => a == b,
            (
                ck::OperationError::DerivedPreservationMismatch { states: a, .. },
                ck::OperationError::DerivedEffectContract { states: b, .. },
            ) => a == b,
            (
                ck::OperationError::DerivedPreservationMismatch { states: a, .. },
                ck::OperationError::DerivedPreservationMismatch { states: b, .. },
            ) => a == b,
            _ => unreachable!(),
        }
    }
}

#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct OperationSemanticPreconditions {
    owner: RetainedError<ck::OperationError>,
}
#[pymethods]
impl OperationSemanticPreconditions {
    #[getter]
    fn bits(&self) -> u16 {
        match self.owner.as_ref() {
            ck::OperationError::SemanticPreconditionContract { missing, .. } => {
                missing.bits() as u16
            }
            _ => unreachable!(),
        }
    }
    fn __repr__(&self) -> String {
        match self.owner.as_ref() {
            ck::OperationError::SemanticPreconditionContract { missing, .. } => {
                format!("{:?}", missing)
            }
            _ => unreachable!(),
        }
    }
    fn __eq__(&self, other: &Self) -> bool {
        match (self.owner.as_ref(), other.owner.as_ref()) {
            (
                ck::OperationError::SemanticPreconditionContract { missing: a, .. },
                ck::OperationError::SemanticPreconditionContract { missing: b, .. },
            ) => a == b,
            _ => unreachable!(),
        }
    }
}

#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct OperationCipStatePolicy {
    owner: RetainedError<ck::OperationError>,
}
#[pymethods]
impl OperationCipStatePolicy {
    #[getter]
    fn discriminant(&self) -> isize {
        match self.owner.as_ref() {
            ck::OperationError::CipStateContract { policy, .. } => *policy as isize,
            _ => unreachable!(),
        }
    }
    fn __repr__(&self) -> String {
        match self.owner.as_ref() {
            ck::OperationError::CipStateContract { policy, .. } => format!("{:?}", policy),
            _ => unreachable!(),
        }
    }
    fn __eq__(&self, other: &Self) -> bool {
        match (self.owner.as_ref(), other.owner.as_ref()) {
            (
                ck::OperationError::CipStateContract { policy: a, .. },
                ck::OperationError::CipStateContract { policy: b, .. },
            ) => a == b,
            _ => unreachable!(),
        }
    }
}

#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct OperationMappingRequirement {
    owner: RetainedError<ck::OperationError>,
}
#[pymethods]
impl OperationMappingRequirement {
    #[getter]
    fn discriminant(&self) -> isize {
        match self.owner.as_ref() {
            ck::OperationError::MappingContract { requirement, .. } => *requirement as isize,
            _ => unreachable!(),
        }
    }
    fn __repr__(&self) -> String {
        match self.owner.as_ref() {
            ck::OperationError::MappingContract { requirement, .. } => format!("{:?}", requirement),
            _ => unreachable!(),
        }
    }
    fn __eq__(&self, other: &Self) -> bool {
        match (self.owner.as_ref(), other.owner.as_ref()) {
            (
                ck::OperationError::MappingContract { requirement: a, .. },
                ck::OperationError::MappingContract { requirement: b, .. },
            ) => a == b,
            _ => unreachable!(),
        }
    }
}

error_projection!(retained OperationError,OperationErrorKind,ck::OperationError,"operation",{
 AtomPropertyIndex{atom:usize=>|value:&ck::AtomId|value.index(),atom_count:usize},
 ReservedAtomPropertyKey{key:String},
 Fragments(_),Scaffold(_),ReactionRun(_),ReactionApply(_),Enumeration(_),Alignment(_),Conformer(_),Tautomer(_),InvalidTopology(_),InvalidTopologyEdit(_),InvalidCoordinates(_),InvalidProperty(_),AtomProperty(_),BondProperty(_),Valence(_),Radical(_),Rings(_),PotentialStereo(_),Stereo(_),CipLabeler(_),AtomCode(_),Transform(_),CoordinateInput(_),Coordinate2D(_),UffOptimization(_),MmffOptimization(_),Kekulize(_),Aromaticity(_),Sanitize(_),Hydrogen(_),
 EmptyFragments, UnsupportedFeature{operation,source},Unsupported{operation},
 OutputMismatch{operation,expected,actual},AccessDenied{operation,block},
 BlockCheckedOut{operation,block},BlockNotCheckedOut{operation,block},
 IncompleteCommit{operation,block},MappingContract{operation,issue,requirement},
 InvalidTopologyMapping{operation,source},AutoRemapContract{operation,block,issue},
 OperationContract{operation,field,issue,expected,actual},
 SemanticPreconditionContract{operation,missing,issue},CoordinateAppendRequiresValues{operation},
 DerivedEffectContract{operation,action:&'static str,states,issue},
 DerivedPreservationMismatch{operation,states,field,expected,actual},
 CipStateContract{operation,policy,issue},
 InvalidPropertyList{target:&'static str,name,values:usize,expected},
 InvalidDerivedCache{state:&'static str,field,actual,expected},
 InvalidAlgorithmResult{operation,field,actual,expected},
 InvalidReconstructionOrigin{operation,entity:&'static str,destination:usize,input:usize,row:usize,input_count:usize,row_count:Option<usize>},
 Algorithm{operation,detail:String},
},cause=project_operation_source,methods={
 #[staticmethod] #[doc=r#"{"role":"error_source"}"#]
 fn cause<'py>(slf:Bound<'py,pyo3::exceptions::PyBaseException>)->PyResult<Option<OperationSource>>{slf.getattr("__cause__")?.extract().map_err(Into::into)}
 #[getter] #[doc=r#"{"role":"error_payload"}"#]
 fn source(&self,py:Python<'_>)->PyResult<Option<Source<UnsupportedFeatureError,crate::error_framework::ErrorObject<crate::model_edit_errors::MappingValidationError>>>>{
  match self.inner.as_ref(){
   ck::OperationError::UnsupportedFeature{..}=>Ok(Some(Source::First(UnsupportedFeatureError{owner:self.inner.clone()}))),
   ck::OperationError::InvalidTopologyMapping{source,..}=>crate::model_edit_errors::MappingValidationError::from_rust(source.clone()).value(py).extract().map(|v|Some(Source::Second(v))).map_err(Into::into),
   _=>Ok(None),
  }
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#]
 fn operation(&self)->Option<Source<&'static str,MoleculeOpSpec>>{match self.inner.as_ref(){
 ck::OperationError::UnsupportedFeature{operation,..}=>Some(Source::Second(MoleculeOpSpec::from_native(operation))),
 ck::OperationError::Unsupported{operation,..}|ck::OperationError::OutputMismatch{operation,..}|ck::OperationError::AccessDenied{operation,..}|ck::OperationError::BlockCheckedOut{operation,..}|ck::OperationError::BlockNotCheckedOut{operation,..}|ck::OperationError::IncompleteCommit{operation,..}|ck::OperationError::MappingContract{operation,..}|ck::OperationError::InvalidTopologyMapping{operation,..}|ck::OperationError::AutoRemapContract{operation,..}|ck::OperationError::OperationContract{operation,..}|ck::OperationError::SemanticPreconditionContract{operation,..}|ck::OperationError::CoordinateAppendRequiresValues{operation,..}|ck::OperationError::DerivedEffectContract{operation,..}|ck::OperationError::DerivedPreservationMismatch{operation,..}|ck::OperationError::CipStateContract{operation,..}|ck::OperationError::InvalidAlgorithmResult{operation,..}|ck::OperationError::InvalidReconstructionOrigin{operation,..}|ck::OperationError::Algorithm{operation,..}=>Some(Source::First(*operation)),_=>None}}
#[getter] fn block(&self)->Option<&'static str>{match self.inner.as_ref(){ck::OperationError::AccessDenied{block,..}|ck::OperationError::BlockCheckedOut{block,..}|ck::OperationError::BlockNotCheckedOut{block,..}|ck::OperationError::IncompleteCommit{block,..}|ck::OperationError::AutoRemapContract{block,..}=>Some(block),_=>None}}
#[getter] fn issue(&self)->Option<&'static str>{match self.inner.as_ref(){ck::OperationError::MappingContract{issue,..}|ck::OperationError::AutoRemapContract{issue,..}|ck::OperationError::OperationContract{issue,..}|ck::OperationError::SemanticPreconditionContract{issue,..}|ck::OperationError::DerivedEffectContract{issue,..}|ck::OperationError::CipStateContract{issue,..}=>Some(issue),_=>None}}
#[getter] fn field(&self)->Option<&str>{match self.inner.as_ref(){
 ck::OperationError::DerivedPreservationMismatch{field,..}=>Some(field.as_str()),
 ck::OperationError::OperationContract{field,..}|ck::OperationError::InvalidDerivedCache{field,..}|ck::OperationError::InvalidAlgorithmResult{field,..}=>Some(field),_=>None}}
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn expected(&self)->Option<Source<usize,Source<String,MoleculeOpOutput>>>{match self.inner.as_ref(){
 ck::OperationError::InvalidDerivedCache{expected,..}|ck::OperationError::InvalidAlgorithmResult{expected,..}|ck::OperationError::InvalidPropertyList{expected,..}=>Some(Source::First(*expected)),
 ck::OperationError::OperationContract{expected,..}=>Some(Source::First(*expected as usize)),
 ck::OperationError::DerivedPreservationMismatch{expected,..}=>Some(Source::Second(Source::First(expected.clone()))),
 ck::OperationError::OutputMismatch{expected,..}=>Some(Source::Second(Source::Second((*expected).into()))),_=>None}}
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn actual(&self)->Option<Source<usize,Source<String,MoleculeOpOutput>>>{match self.inner.as_ref(){
 ck::OperationError::InvalidDerivedCache{actual,..}|ck::OperationError::InvalidAlgorithmResult{actual,..}=>Some(Source::First(*actual)),
 ck::OperationError::OperationContract{actual,..}=>Some(Source::First(*actual as usize)),
 ck::OperationError::DerivedPreservationMismatch{actual,..}=>Some(Source::Second(Source::First(actual.clone()))),
 ck::OperationError::OutputMismatch{actual,..}=>Some(Source::Second(Source::Second((*actual).into()))),_=>None}}
#[getter] fn name(&self,py:Python<'_>)->PyResult<Option<String>>{match self.inner.as_ref(){ck::OperationError::InvalidPropertyList{name,..}=>crate::property_values::decode_source_text(py,name).map(Some),_=>Ok(None)}}
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn states(&self)->Option<OperationDerivedState>{match self.inner.as_ref(){ck::OperationError::DerivedEffectContract{..}|ck::OperationError::DerivedPreservationMismatch{..}=>Some(OperationDerivedState{owner:self.inner.clone()}),_=>None}}
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn missing(&self)->Option<OperationSemanticPreconditions>{match self.inner.as_ref(){ck::OperationError::SemanticPreconditionContract{..}=>Some(OperationSemanticPreconditions{owner:self.inner.clone()}),_=>None}}
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn policy(&self)->Option<OperationCipStatePolicy>{match self.inner.as_ref(){ck::OperationError::CipStateContract{..}=>Some(OperationCipStatePolicy{owner:self.inner.clone()}),_=>None}}
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn requirement(&self)->Option<OperationMappingRequirement>{match self.inner.as_ref(){ck::OperationError::MappingContract{..}=>Some(OperationMappingRequirement{owner:self.inner.clone()}),_=>None}}
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn fragments(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<MoleculeFragmentsError>>>{
 match self.inner.as_ref(){ck::OperationError::Fragments(value)=>{
 let projected=MoleculeFragmentsError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn scaffold(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<ScaffoldError>>>{
 match self.inner.as_ref(){ck::OperationError::Scaffold(value)=>{
 let projected=ScaffoldError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn reaction_run(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<ReactionRunError>>>{
 match self.inner.as_ref(){ck::OperationError::ReactionRun(value)=>{
 let projected=ReactionRunError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn reaction_apply(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<ReactionApplyError>>>{
 match self.inner.as_ref(){ck::OperationError::ReactionApply(value)=>{
 let projected=ReactionApplyError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn enumeration(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<EnumerationRunError>>>{
 match self.inner.as_ref(){ck::OperationError::Enumeration(value)=>{
 let projected=EnumerationRunError::from_ref(self.inner.node().project(|error|{let value=error.downcast_ref::<ck::OperationError>().expect("operation payload owner");let ck::OperationError::Enumeration(child)=value else{unreachable!()};child}));
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn alignment(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<AlignmentError>>>{
 match self.inner.as_ref(){ck::OperationError::Alignment(value)=>{
 let projected=AlignmentError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn conformer(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<ConformerRunError>>>{
 match self.inner.as_ref(){ck::OperationError::Conformer(value)=>{
 let projected=ConformerRunError::from_ref(self.inner.node().project(|error|{let value=error.downcast_ref::<ck::OperationError>().expect("operation payload owner");let ck::OperationError::Conformer(child)=value else{unreachable!()};child}));
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn tautomer(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<TautomerRunError>>>{
 match self.inner.as_ref(){ck::OperationError::Tautomer(value)=>{
 let projected=TautomerRunError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn invalid_topology(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<TopologyValidationError>>>{
 match self.inner.as_ref(){ck::OperationError::InvalidTopology(value)=>{
 let projected=TopologyValidationError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn invalid_topology_edit(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<TopologyEditError>>>{
 match self.inner.as_ref(){ck::OperationError::InvalidTopologyEdit(value)=>{
 let projected=TopologyEditError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn invalid_coordinates(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<CoordinateValidationError>>>{
 match self.inner.as_ref(){ck::OperationError::InvalidCoordinates(value)=>{
 let projected=CoordinateValidationError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn invalid_property(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<MoleculePropertyError>>>{
 match self.inner.as_ref(){ck::OperationError::InvalidProperty(value)=>{
 let projected=MoleculePropertyError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn atom_property(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<AtomPropertyError>>>{
 match self.inner.as_ref(){ck::OperationError::AtomProperty(value)=>{
 let projected=AtomPropertyError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn bond_property(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<BondValueError>>>{
 match self.inner.as_ref(){ck::OperationError::BondProperty(value)=>{
 let projected=BondValueError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn valence(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>>{
 match self.inner.as_ref(){ck::OperationError::Valence(value)=>{
 let projected=ValenceError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn radical(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<RadicalError>>>{
 match self.inner.as_ref(){ck::OperationError::Radical(value)=>{
 let projected=RadicalError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn rings(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<RingFindingError>>>{
 match self.inner.as_ref(){ck::OperationError::Rings(value)=>{
 let projected=RingFindingError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn potential_stereo(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<PotentialStereoError>>>{
 match self.inner.as_ref(){ck::OperationError::PotentialStereo(value)=>{
 let projected=PotentialStereoError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn stereo(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<StereoError>>>{
 match self.inner.as_ref(){ck::OperationError::Stereo(value)=>{
 let projected=StereoError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn cip_labeler(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<CipLabelerError>>>{
 match self.inner.as_ref(){ck::OperationError::CipLabeler(value)=>{
 let projected=CipLabelerError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn atom_code(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<AtomCodeError>>>{
 match self.inner.as_ref(){ck::OperationError::AtomCode(value)=>{
 let projected=AtomCodeError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn transform(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<TransformError>>>{
 match self.inner.as_ref(){ck::OperationError::Transform(value)=>{
 let projected=TransformError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn coordinate_input(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<CoordinateInputError>>>{
 match self.inner.as_ref(){ck::OperationError::CoordinateInput(value)=>{
 let projected=CoordinateInputError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn coordinate2_d(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<Coordinate2DError>>>{
 match self.inner.as_ref(){ck::OperationError::Coordinate2D(value)=>{
 let projected=Coordinate2DError::from_ref(self.inner.node().project(|error|{let value=error.downcast_ref::<ck::OperationError>().expect("operation payload owner");let ck::OperationError::Coordinate2D(child)=value else{unreachable!()};child}));
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn uff_optimization(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<UffOptimizationError>>>{
 match self.inner.as_ref(){ck::OperationError::UffOptimization(value)=>{
 let projected=UffOptimizationError::from_ref(self.inner.node().project(|error|{let value=error.downcast_ref::<ck::OperationError>().expect("operation payload owner");let ck::OperationError::UffOptimization(child)=value else{unreachable!()};child}));
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn mmff_optimization(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<MmffOptimizationError>>>{
 match self.inner.as_ref(){ck::OperationError::MmffOptimization(value)=>{
 let projected=MmffOptimizationError::from_ref(self.inner.node().project(|error|{let value=error.downcast_ref::<ck::OperationError>().expect("operation payload owner");let ck::OperationError::MmffOptimization(child)=value else{unreachable!()};child}));
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn kekulize(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<KekulizeError>>>{
 match self.inner.as_ref(){ck::OperationError::Kekulize(value)=>{
 let projected=KekulizeError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn aromaticity(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<AromaticityError>>>{
 match self.inner.as_ref(){ck::OperationError::Aromaticity(value)=>{
 let projected=AromaticityError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn sanitize(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<SanitizeError>>>{
 match self.inner.as_ref(){ck::OperationError::Sanitize(value)=>{
 let projected=SanitizeError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
#[getter] #[doc=r#"{"role":"error_payload"}"#] fn hydrogen(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<HydrogenError>>>{
 match self.inner.as_ref(){ck::OperationError::Hydrogen(value)=>{
 let projected=HydrogenError::from_rust(value.clone());
 projected.value(py).extract().map(Some).map_err(Into::into)},_=>Ok(None)}
 }
});
type SmilesErrorSource = Source<
    crate::error_framework::ErrorObject<SmilesParseError>,
    Source<
        crate::error_framework::ErrorObject<HydrogenError>,
        Source<
            crate::error_framework::ErrorObject<SanitizeError>,
            Source<
                crate::error_framework::ErrorObject<SmilesStereoError>,
                crate::error_framework::ErrorObject<OperationError>,
            >,
        >,
    >,
>;
error_projection!(retained SmilesError,SmilesErrorKind,ck::SmilesError,"smiles",{
 Parse(_),Hydrogen(_),Sanitize(_),Stereo(_),Construction(_),
},cause=|inner:&RetainedError<ck::SmilesError>|{
 let node=inner.node().source().expect("all native SMILES variants have an actual source");
 Some(match inner.as_ref(){
 ck::SmilesError::Parse(value)=>SmilesParseError::from_rust(value.clone()),
 ck::SmilesError::Hydrogen(value)=>HydrogenError::from_rust(value.clone()),
 ck::SmilesError::Sanitize(value)=>SanitizeError::from_rust(value.clone()),
 ck::SmilesError::Stereo(value)=>SmilesStereoError::from_rust(value.clone()),
 ck::SmilesError::Construction(_)=>OperationError::from_ref(node),
 })
},source=Option<SmilesErrorSource>);
