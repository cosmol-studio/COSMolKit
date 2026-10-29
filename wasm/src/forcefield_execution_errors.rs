//! Native force-field execution failures retain actual sources and payloads.
use crate::errors::{enum_error, error_projection, exception};
use crate::forcefield_kernel_errors::ForceFieldKernelError;
use crate::fragment_errors::MoleculeFragmentsError;
use crate::match_error::MatchError;
use crate::mmff_errors::{MmffMolPropertiesError, MmffParamParseError};
use crate::model_errors::{
    AtomPropertyError, BondValueError, MappingValidationError, MoleculePropertyError,
    StereoGroupError, TemplateAttachmentOrderError, TopologyValidationError,
};
use crate::path_errors::PathError;
use crate::query_compile_error::QueryCompileError;
use crate::ring_errors::RingFindingError;
use crate::smarts_parse_errors::SmartsParseError;
use crate::storage::ReadonlyValue;
use crate::uff_errors::{
    UffBuilderError, UffDefaultTorsionQueryError, UffParamError, UffParameterError,
    UffTorsionBondQueryError, UffTypingError,
};
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use std::error::Error;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "FragmentCoordinateViewError | NonbondedFragmentMappingError | NonbondedAssemblyError | ForceFieldConstructionError | AutomaticForceFieldConstructionError | OptimizationStageError | SerialConformerOptimizationError | SingleConformerOptimizationError | SerialUffOptimizationError | UffThreadCountError | UffPreparedOptimizationError | MmffBuilderError | MmffOptimizationFailure | SerialResultCapacityError | SerialCoordinateCountError | UffSingleError | UffConformerError | MmffOptimizationSourceError | MmffOptimizationError | UffEvaluationError | UffOptimizationError | MoleculeFragmentsError | UffBuilderError | UffParamError | UffTypingError | UffParameterError | UffTorsionBondQueryError | UffDefaultTorsionQueryError | ForceFieldKernelError | MmffMolPropertiesError | MmffParamParseError | SmartsParseError | QueryCompileError | MatchError | PathError | TopologyValidationError | MappingValidationError | MoleculePropertyError | StereoGroupError | TemplateAttachmentOrderError | AtomPropertyError | BondValueError | RingFindingError | ValenceError"
    )]
    #[derive(Clone)]
    pub type ExecutionCause;
}

// Only actual native sources from this execution graph enter this projector.
fn project(owner: ReadonlyValue<dyn Error>) -> ExecutionCause {
    if owner.is::<ck::FragmentCoordinateViewError>() {
        return FragmentCoordinateViewError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::FragmentCoordinateViewError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::NonbondedFragmentMappingError>() {
        return NonbondedFragmentMappingError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::NonbondedFragmentMappingError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::NonbondedAssemblyError>() {
        return NonbondedAssemblyError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::NonbondedAssemblyError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::ForceFieldConstructionError>() {
        return ForceFieldConstructionError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::ForceFieldConstructionError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::AutomaticForceFieldConstructionError>() {
        return AutomaticForceFieldConstructionError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::AutomaticForceFieldConstructionError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::OptimizationStageError>() {
        return OptimizationStageError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::OptimizationStageError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::SerialConformerOptimizationError>() {
        return SerialConformerOptimizationError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::SerialConformerOptimizationError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::SingleConformerOptimizationError>() {
        return SingleConformerOptimizationError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::SingleConformerOptimizationError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::SerialUffOptimizationError>() {
        return SerialUffOptimizationError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::SerialUffOptimizationError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::UffThreadCountError>() {
        return UffThreadCountError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::UffThreadCountError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::UffPreparedOptimizationError>() {
        return UffPreparedOptimizationError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::UffPreparedOptimizationError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::MmffBuilderError>() {
        return MmffBuilderError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::MmffBuilderError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::MmffOptimizationFailure>() {
        return MmffOptimizationFailure::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::MmffOptimizationFailure>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::SerialResultCapacityError>() {
        return SerialResultCapacityError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::SerialResultCapacityError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::SerialCoordinateCountError>() {
        return SerialCoordinateCountError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::SerialCoordinateCountError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::UffSingleError>() {
        return UffSingleError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::UffSingleError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::UffConformerError>() {
        return UffConformerError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::UffConformerError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::MmffOptimizationSourceError>() {
        return MmffOptimizationSourceError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::MmffOptimizationSourceError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::MmffOptimizationError>() {
        return MmffOptimizationError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::MmffOptimizationError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::UffEvaluationError>() {
        return UffEvaluationError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::UffEvaluationError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::UffOptimizationError>() {
        return UffOptimizationError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::UffOptimizationError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::MoleculeFragmentsError>() {
        return MoleculeFragmentsError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::MoleculeFragmentsError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::UffBuilderError>() {
        return UffBuilderError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::UffBuilderError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::UffParamError>() {
        return UffParamError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::UffParamError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::UffTypingError>() {
        return UffTypingError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::UffTypingError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::UffParameterError>() {
        return UffParameterError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::UffParameterError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::UffTorsionBondQueryError>() {
        return UffTorsionBondQueryError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::UffTorsionBondQueryError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::UffDefaultTorsionQueryError>() {
        return UffDefaultTorsionQueryError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::UffDefaultTorsionQueryError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::ForceFieldKernelError>() {
        return ForceFieldKernelError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::ForceFieldKernelError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::MmffMolPropertiesError>() {
        return MmffMolPropertiesError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::MmffMolPropertiesError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::MmffParamParseError>() {
        return MmffParamParseError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::MmffParamParseError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::SmartsParseError>() {
        return SmartsParseError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::SmartsParseError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::QueryCompileError>() {
        return QueryCompileError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::QueryCompileError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::MatchError>() {
        return MatchError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::MatchError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::PathError>() {
        return PathError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::PathError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::TopologyValidationError>() {
        return TopologyValidationError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::TopologyValidationError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::MappingValidationError>() {
        return MappingValidationError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::MappingValidationError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::MoleculePropertyError>() {
        return MoleculePropertyError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::MoleculePropertyError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::StereoGroupError>() {
        return StereoGroupError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::StereoGroupError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::TemplateAttachmentOrderError>() {
        return TemplateAttachmentOrderError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::TemplateAttachmentOrderError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::AtomPropertyError>() {
        return AtomPropertyError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::AtomPropertyError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::BondValueError>() {
        return BondValueError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::BondValueError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::RingFindingError>() {
        return RingFindingError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::RingFindingError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    if owner.is::<ck::ValenceError>() {
        return ValenceError::from_owner(owner.child(|e| {
            e.downcast_ref::<ck::ValenceError>()
                .expect("selected native source")
        }))
        .unchecked_into();
    }
    unreachable!("unmodeled concrete force-field source")
}
fn cause<T: Error + 'static>(owner: &ReadonlyValue<T>) -> Option<ExecutionCause> {
    owner
        .source()
        .map(|_| project(owner.child(|e| e.source().expect("immutable native source"))))
}
enum_error!(FragmentCoordinateViewError,FragmentCoordinateViewErrorKind,ck::FragmentCoordinateViewError,{SelectedConformerMissing=>ck::FragmentCoordinateViewError::SelectedConformerMissing{..},SelectedConformerDuplicated=>ck::FragmentCoordinateViewError::SelectedConformerDuplicated{..},KernelCoordinateRowCount=>ck::FragmentCoordinateViewError::KernelCoordinateRowCount{..},KernelCoordinateRowWidth=>ck::FragmentCoordinateViewError::KernelCoordinateRowWidth{..}},{
#[wasm_bindgen(getter,js_name=id)]pub fn id(&self)->Option<usize>{match self.inner.as_ref(){ck::FragmentCoordinateViewError::SelectedConformerMissing{id,..}=>Some(*id),ck::FragmentCoordinateViewError::SelectedConformerDuplicated{id,..}=>Some(*id),ck::FragmentCoordinateViewError::KernelCoordinateRowCount{id,..}=>Some(*id),ck::FragmentCoordinateViewError::KernelCoordinateRowWidth{id,..}=>Some(*id),_=>None}}
#[wasm_bindgen(getter,js_name=expected)]pub fn expected(&self)->Option<usize>{match self.inner.as_ref(){ck::FragmentCoordinateViewError::KernelCoordinateRowCount{expected,..}=>Some(*expected),_=>None}}
#[wasm_bindgen(getter,js_name=actual)]pub fn actual(&self)->Option<usize>{match self.inner.as_ref(){ck::FragmentCoordinateViewError::KernelCoordinateRowCount{actual,..}=>Some(*actual),ck::FragmentCoordinateViewError::KernelCoordinateRowWidth{actual,..}=>Some(*actual),_=>None}}
#[wasm_bindgen(getter,js_name=row)]pub fn row(&self)->Option<usize>{match self.inner.as_ref(){ck::FragmentCoordinateViewError::KernelCoordinateRowWidth{row,..}=>Some(*row),_=>None}}
},source=ExecutionCause,cause=cause::<ck::FragmentCoordinateViewError>);
enum_error!(NonbondedFragmentMappingError,NonbondedFragmentMappingErrorKind,ck::NonbondedFragmentMappingError,{FragmentCopy=>ck::NonbondedFragmentMappingError::FragmentCopy(_),CoordinateView=>ck::NonbondedFragmentMappingError::CoordinateView(_),RequestedMappingMissing=>ck::NonbondedFragmentMappingError::RequestedMappingMissing,AtomOutOfRange=>ck::NonbondedFragmentMappingError::AtomOutOfRange{..},DuplicateAtom=>ck::NonbondedFragmentMappingError::DuplicateAtom{..},MissingAtom=>ck::NonbondedFragmentMappingError::MissingAtom{..}},{
#[wasm_bindgen(getter,js_name=componentIndex)]pub fn component_index(&self)->Option<usize>{match self.inner.as_ref(){ck::NonbondedFragmentMappingError::AtomOutOfRange{component_index,..}=>Some(*component_index),_=>None}}
#[wasm_bindgen(getter,js_name=atomIndex)]pub fn atom_index(&self)->Option<usize>{match self.inner.as_ref(){ck::NonbondedFragmentMappingError::AtomOutOfRange{atom_index,..}=>Some(*atom_index),ck::NonbondedFragmentMappingError::DuplicateAtom{atom_index,..}=>Some(*atom_index),ck::NonbondedFragmentMappingError::MissingAtom{atom_index,..}=>Some(*atom_index),_=>None}}
#[wasm_bindgen(getter,js_name=atomCount)]pub fn atom_count(&self)->Option<usize>{match self.inner.as_ref(){ck::NonbondedFragmentMappingError::AtomOutOfRange{atom_count,..}=>Some(*atom_count),_=>None}}
#[wasm_bindgen(getter,js_name=firstComponent)]pub fn first_component(&self)->Option<usize>{match self.inner.as_ref(){ck::NonbondedFragmentMappingError::DuplicateAtom{first_component,..}=>Some(*first_component),_=>None}}
#[wasm_bindgen(getter,js_name=secondComponent)]pub fn second_component(&self)->Option<usize>{match self.inner.as_ref(){ck::NonbondedFragmentMappingError::DuplicateAtom{second_component,..}=>Some(*second_component),_=>None}}
#[wasm_bindgen(getter,js_name=fragmentCopy)]pub fn fragment_copy(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::NonbondedFragmentMappingError::FragmentCopy(_)=>Some(MoleculeFragmentsError::from_owner(self.inner.child(|e|match e{ck::NonbondedFragmentMappingError::FragmentCopy(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=coordinateView)]pub fn coordinate_view(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::NonbondedFragmentMappingError::CoordinateView(_)=>Some(FragmentCoordinateViewError::from_owner(self.inner.child(|e|match e{ck::NonbondedFragmentMappingError::CoordinateView(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
},source=ExecutionCause,cause=cause::<ck::NonbondedFragmentMappingError>);
enum_error!(NonbondedAssemblyError,NonbondedAssemblyErrorKind,ck::NonbondedAssemblyError,{FragmentMapping=>ck::NonbondedAssemblyError::FragmentMapping(_),PairBuilder=>ck::NonbondedAssemblyError::PairBuilder(_)},{
#[wasm_bindgen(getter,js_name=fragmentMapping)]pub fn fragment_mapping(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::NonbondedAssemblyError::FragmentMapping(_)=>Some(NonbondedFragmentMappingError::from_owner(self.inner.child(|e|match e{ck::NonbondedAssemblyError::FragmentMapping(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=pairBuilder)]pub fn pair_builder(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::NonbondedAssemblyError::PairBuilder(_)=>Some(UffBuilderError::from_owner(self.inner.child(|e|match e{ck::NonbondedAssemblyError::PairBuilder(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
},source=ExecutionCause,cause=cause::<ck::NonbondedAssemblyError>);
enum_error!(ForceFieldConstructionError,ForceFieldConstructionErrorKind,ck::ForceFieldConstructionError,{Builder=>ck::ForceFieldConstructionError::Builder(_),Nonbonded=>ck::ForceFieldConstructionError::Nonbonded(_)},{
#[wasm_bindgen(getter,js_name=builder)]pub fn builder(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::ForceFieldConstructionError::Builder(_)=>Some(UffBuilderError::from_owner(self.inner.child(|e|match e{ck::ForceFieldConstructionError::Builder(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=nonbonded)]pub fn nonbonded(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::ForceFieldConstructionError::Nonbonded(_)=>Some(NonbondedAssemblyError::from_owner(self.inner.child(|e|match e{ck::ForceFieldConstructionError::Nonbonded(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
},source=ExecutionCause,cause=cause::<ck::ForceFieldConstructionError>);
enum_error!(AutomaticForceFieldConstructionError,AutomaticForceFieldConstructionErrorKind,ck::AutomaticForceFieldConstructionError,{ParameterTable=>ck::AutomaticForceFieldConstructionError::ParameterTable(_),Typing=>ck::AutomaticForceFieldConstructionError::Typing(_),Construction=>ck::AutomaticForceFieldConstructionError::Construction(_)},{
#[wasm_bindgen(getter,js_name=parameterTable)]pub fn parameter_table(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::AutomaticForceFieldConstructionError::ParameterTable(_)=>Some(UffParamError::from_owner(self.inner.child(|e|match e{ck::AutomaticForceFieldConstructionError::ParameterTable(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=typing)]pub fn typing(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::AutomaticForceFieldConstructionError::Typing(_)=>Some(UffTypingError::from_owner(self.inner.child(|e|match e{ck::AutomaticForceFieldConstructionError::Typing(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=construction)]pub fn construction(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::AutomaticForceFieldConstructionError::Construction(_)=>Some(ForceFieldConstructionError::from_owner(self.inner.child(|e|match e{ck::AutomaticForceFieldConstructionError::Construction(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
},source=ExecutionCause,cause=cause::<ck::AutomaticForceFieldConstructionError>);
enum_error!(OptimizationStageError,OptimizationStageErrorKind,ck::OptimizationStageError,{Initialize=>ck::OptimizationStageError::Initialize(_),Minimize=>ck::OptimizationStageError::Minimize(_),FinalEnergy=>ck::OptimizationStageError::FinalEnergy(_)},{
#[wasm_bindgen(getter,js_name=initialize)]pub fn initialize(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::OptimizationStageError::Initialize(_)=>Some(ForceFieldKernelError::from_owner(self.inner.child(|e|match e{ck::OptimizationStageError::Initialize(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=minimize)]pub fn minimize(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::OptimizationStageError::Minimize(_)=>Some(ForceFieldKernelError::from_owner(self.inner.child(|e|match e{ck::OptimizationStageError::Minimize(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=finalEnergy)]pub fn final_energy(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::OptimizationStageError::FinalEnergy(_)=>Some(ForceFieldKernelError::from_owner(self.inner.child(|e|match e{ck::OptimizationStageError::FinalEnergy(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
},source=ExecutionCause,cause=cause::<ck::OptimizationStageError>);
enum_error!(SerialConformerOptimizationError,SerialConformerOptimizationErrorKind,ck::SerialConformerOptimizationError,{ResultCapacity=>ck::SerialConformerOptimizationError::ResultCapacity(_),CoordinateCount=>ck::SerialConformerOptimizationError::CoordinateCount{..},Optimization=>ck::SerialConformerOptimizationError::Optimization{..}},{
#[wasm_bindgen(getter,js_name=inputIndex)]pub fn input_index(&self)->Option<usize>{match self.inner.as_ref(){ck::SerialConformerOptimizationError::CoordinateCount{input_index,..}=>Some(*input_index),ck::SerialConformerOptimizationError::Optimization{input_index,..}=>Some(*input_index),_=>None}}
#[wasm_bindgen(getter,js_name=conformerId)]pub fn conformer_id(&self)->Option<usize>{match self.inner.as_ref(){ck::SerialConformerOptimizationError::Optimization{conformer_id,..}=>Some(*conformer_id),_=>None}}
#[wasm_bindgen(getter,js_name=resultCapacity)]pub fn result_capacity(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::SerialConformerOptimizationError::ResultCapacity(_)=>Some(SerialResultCapacityError::from_owner(self.inner.child(|e|match e{ck::SerialConformerOptimizationError::ResultCapacity(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter)]pub fn source(&self)->Option<ExecutionCause>{cause(&self.inner)}
},source=ExecutionCause,cause=cause::<ck::SerialConformerOptimizationError>);
enum_error!(SingleConformerOptimizationError,SingleConformerOptimizationErrorKind,ck::SingleConformerOptimizationError,{Construction=>ck::SingleConformerOptimizationError::Construction(_),Optimization=>ck::SingleConformerOptimizationError::Optimization(_)},{
#[wasm_bindgen(getter,js_name=construction)]pub fn construction(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::SingleConformerOptimizationError::Construction(_)=>Some(AutomaticForceFieldConstructionError::from_owner(self.inner.child(|e|match e{ck::SingleConformerOptimizationError::Construction(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=optimization)]pub fn optimization(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::SingleConformerOptimizationError::Optimization(_)=>Some(OptimizationStageError::from_owner(self.inner.child(|e|match e{ck::SingleConformerOptimizationError::Optimization(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
},source=ExecutionCause,cause=cause::<ck::SingleConformerOptimizationError>);
enum_error!(SerialUffOptimizationError,SerialUffOptimizationErrorKind,ck::SerialUffOptimizationError,{Construction=>ck::SerialUffOptimizationError::Construction(_),Optimization=>ck::SerialUffOptimizationError::Optimization(_)},{
#[wasm_bindgen(getter,js_name=construction)]pub fn construction(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::SerialUffOptimizationError::Construction(_)=>Some(AutomaticForceFieldConstructionError::from_owner(self.inner.child(|e|match e{ck::SerialUffOptimizationError::Construction(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=optimization)]pub fn optimization(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::SerialUffOptimizationError::Optimization(_)=>Some(SerialConformerOptimizationError::from_owner(self.inner.child(|e|match e{ck::SerialUffOptimizationError::Optimization(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
},source=ExecutionCause,cause=cause::<ck::SerialUffOptimizationError>);
enum_error!(UffThreadCountError,UffThreadCountErrorKind,ck::UffThreadCountError,{UndefinedSignedNegation=>ck::UffThreadCountError::UndefinedSignedNegation},{
},source=ExecutionCause,cause=cause::<ck::UffThreadCountError>);
enum_error!(UffPreparedOptimizationError,UffPreparedOptimizationErrorKind,ck::UffPreparedOptimizationError,{Preparation=>ck::UffPreparedOptimizationError::Preparation(_),Single=>ck::UffPreparedOptimizationError::Single(_),Serial=>ck::UffPreparedOptimizationError::Serial(_),ConformerStage=>ck::UffPreparedOptimizationError::ConformerStage(_)},{
#[wasm_bindgen(getter,js_name=preparation)]pub fn preparation(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::UffPreparedOptimizationError::Preparation(_)=>Some(UffParameterError::from_owner(self.inner.child(|e|match e{ck::UffPreparedOptimizationError::Preparation(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=single)]pub fn single(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::UffPreparedOptimizationError::Single(_)=>Some(SingleConformerOptimizationError::from_owner(self.inner.child(|e|match e{ck::UffPreparedOptimizationError::Single(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=serial)]pub fn serial(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::UffPreparedOptimizationError::Serial(_)=>Some(SerialUffOptimizationError::from_owner(self.inner.child(|e|match e{ck::UffPreparedOptimizationError::Serial(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=conformerStage)]pub fn conformer_stage(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::UffPreparedOptimizationError::ConformerStage(_)=>Some(SerialConformerOptimizationError::from_owner(self.inner.child(|e|match e{ck::UffPreparedOptimizationError::ConformerStage(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
},source=ExecutionCause,cause=cause::<ck::UffPreparedOptimizationError>);
enum_error!(MmffBuilderError,MmffBuilderErrorKind,ck::MmffBuilderError,{MolProperties=>ck::MmffBuilderError::MolProperties(_),ForceFieldKernel=>ck::MmffBuilderError::ForceFieldKernel(_),TorsionBondQuery=>ck::MmffBuilderError::TorsionBondQuery(_),Fragment=>ck::MmffBuilderError::Fragment(_),Components=>ck::MmffBuilderError::Components(_),Valence=>ck::MmffBuilderError::Valence(_),Rings=>ck::MmffBuilderError::Rings(_),Missing3dConformer=>ck::MmffBuilderError::Missing3dConformer{..},InvalidConformerId=>ck::MmffBuilderError::InvalidConformerId{..}},{
#[wasm_bindgen(getter,js_name=confId)]pub fn conf_id(&self)->Option<isize>{match self.inner.as_ref(){ck::MmffBuilderError::Missing3dConformer{conf_id,..}=>Some(*conf_id),ck::MmffBuilderError::InvalidConformerId{conf_id,..}=>Some(*conf_id),_=>None}}
#[wasm_bindgen(getter,js_name=molProperties)]pub fn mol_properties(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::MmffBuilderError::MolProperties(_)=>Some(MmffMolPropertiesError::from_owner(self.inner.child(|e|match e{ck::MmffBuilderError::MolProperties(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=forceFieldKernel)]pub fn force_field_kernel(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::MmffBuilderError::ForceFieldKernel(_)=>Some(ForceFieldKernelError::from_owner(self.inner.child(|e|match e{ck::MmffBuilderError::ForceFieldKernel(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=torsionBondQuery)]pub fn torsion_bond_query(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::MmffBuilderError::TorsionBondQuery(_)=>Some(UffTorsionBondQueryError::from_owner(self.inner.child(|e|match e{ck::MmffBuilderError::TorsionBondQuery(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=fragment)]pub fn fragment(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::MmffBuilderError::Fragment(_)=>Some(NonbondedFragmentMappingError::from_owner(self.inner.child(|e|match e{ck::MmffBuilderError::Fragment(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=components)]pub fn components(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::MmffBuilderError::Components(_)=>Some(PathError::from_owner(self.inner.child(|e|match e{ck::MmffBuilderError::Components(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=valence)]pub fn valence(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::MmffBuilderError::Valence(_)=>Some(ValenceError::from_owner(self.inner.child(|e|match e{ck::MmffBuilderError::Valence(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=rings)]pub fn rings(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::MmffBuilderError::Rings(_)=>Some(RingFindingError::from_owner(self.inner.child(|e|match e{ck::MmffBuilderError::Rings(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
},source=ExecutionCause,cause=cause::<ck::MmffBuilderError>);
enum_error!(MmffOptimizationFailure,MmffOptimizationFailureKind,ck::MmffOptimizationFailure,{Properties=>ck::MmffOptimizationFailure::Properties(_),Builder=>ck::MmffOptimizationFailure::Builder(_),Kernel=>ck::MmffOptimizationFailure::Kernel(_),Stage=>ck::MmffOptimizationFailure::Stage(_),Serial=>ck::MmffOptimizationFailure::Serial(_)},{
#[wasm_bindgen(getter,js_name=properties)]pub fn properties(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::MmffOptimizationFailure::Properties(_)=>Some(MmffMolPropertiesError::from_owner(self.inner.child(|e|match e{ck::MmffOptimizationFailure::Properties(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=builder)]pub fn builder(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::MmffOptimizationFailure::Builder(_)=>Some(MmffBuilderError::from_owner(self.inner.child(|e|match e{ck::MmffOptimizationFailure::Builder(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=kernel)]pub fn kernel(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::MmffOptimizationFailure::Kernel(_)=>Some(ForceFieldKernelError::from_owner(self.inner.child(|e|match e{ck::MmffOptimizationFailure::Kernel(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=stage)]pub fn stage(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::MmffOptimizationFailure::Stage(_)=>Some(OptimizationStageError::from_owner(self.inner.child(|e|match e{ck::MmffOptimizationFailure::Stage(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=serial)]pub fn serial(&self)->Option<ExecutionCause>{match self.inner.as_ref(){ck::MmffOptimizationFailure::Serial(_)=>Some(SerialConformerOptimizationError::from_owner(self.inner.child(|e|match e{ck::MmffOptimizationFailure::Serial(v)=>v,_=>unreachable!("immutable payload")})).unchecked_into()),_=>None}}
},source=ExecutionCause,cause=cause::<ck::MmffOptimizationFailure>);
error_projection!(SerialResultCapacityError, ck::SerialResultCapacityError, {
    #[wasm_bindgen(getter,js_name=conformers)]
    pub fn conformers(&self) -> usize {
        self.inner.conformers()
    }
    #[wasm_bindgen(getter,js_name=resultSlots)]
    pub fn result_slots(&self) -> usize {
        self.inner.result_slots()
    }
});
error_projection!(
    SerialCoordinateCountError,
    ck::SerialCoordinateCountError,
    {
        #[wasm_bindgen(getter,js_name=conformerId)]
        pub fn conformer_id(&self) -> usize {
            self.inner.conformer_id()
        }
        #[wasm_bindgen(getter,js_name=atoms)]
        pub fn atoms(&self) -> usize {
            self.inner.atoms()
        }
        #[wasm_bindgen(getter,js_name=coordinates)]
        pub fn coordinates(&self) -> usize {
            self.inner.coordinates()
        }
    }
);
error_projection!(
    UffSingleError,
    ck::UffSingleError,
    {},
    source = ExecutionCause,
    cause = cause::<ck::UffSingleError>
);
error_projection!(
    UffConformerError,
    ck::UffConformerError,
    {},
    source = ExecutionCause,
    cause = cause::<ck::UffConformerError>
);
error_projection!(
    MmffOptimizationSourceError,
    ck::MmffOptimizationSourceError,
    {},
    source = ExecutionCause,
    cause = cause::<ck::MmffOptimizationSourceError>
);
error_projection!(
    MmffOptimizationError,
    ck::MmffOptimizationError,
    {},
    source = ExecutionCause,
    cause = cause::<ck::MmffOptimizationError>
);

#[wasm_bindgen]
pub enum UffEvaluationFailureKind {
    MissingConformer,
    Preparation,
    Construction,
    Kernel,
}
#[wasm_bindgen]
pub struct UffEvaluationFailure(ReadonlyValue<ck::UffEvaluationFailure>);
#[wasm_bindgen]
impl UffEvaluationFailure {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> UffEvaluationFailureKind {
        match self.0.as_ref() {
            ck::UffEvaluationFailure::MissingConformer { .. } => {
                UffEvaluationFailureKind::MissingConformer
            }
            ck::UffEvaluationFailure::Preparation(_) => UffEvaluationFailureKind::Preparation,
            ck::UffEvaluationFailure::Construction(_) => UffEvaluationFailureKind::Construction,
            ck::UffEvaluationFailure::Kernel(_) => UffEvaluationFailureKind::Kernel,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn requested(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::UffEvaluationFailure::MissingConformer { requested } => *requested,
            _ => None,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn preparation(&self) -> Option<ExecutionCause> {
        match self.0.as_ref() {
            ck::UffEvaluationFailure::Preparation(_) => Some(
                UffParameterError::from_owner(self.0.child(|e| match e {
                    ck::UffEvaluationFailure::Preparation(v) => v,
                    _ => unreachable!("immutable failure"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn construction(&self) -> Option<ExecutionCause> {
        match self.0.as_ref() {
            ck::UffEvaluationFailure::Construction(_) => Some(
                AutomaticForceFieldConstructionError::from_owner(self.0.child(|e| match e {
                    ck::UffEvaluationFailure::Construction(v) => v,
                    _ => unreachable!("immutable failure"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn kernel(&self) -> Option<ExecutionCause> {
        match self.0.as_ref() {
            ck::UffEvaluationFailure::Kernel(_) => Some(
                ForceFieldKernelError::from_owner(self.0.child(|e| match e {
                    ck::UffEvaluationFailure::Kernel(v) => v,
                    _ => unreachable!("immutable failure"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
}
error_projection!(
    UffEvaluationError,
    ck::UffEvaluationError,
    {
        #[wasm_bindgen(getter)]
        pub fn failure(&self) -> UffEvaluationFailure {
            UffEvaluationFailure(self.inner.child(|e| e.failure()))
        }
    },
    source = ExecutionCause,
    cause = cause::<ck::UffEvaluationError>
);
#[wasm_bindgen]
pub enum UffOptimizationErrorVariant {
    MissingConformer,
    Rings,
    Optimization,
    ConformerOptimization,
    Evaluation,
}
#[wasm_bindgen]
pub struct UffOptimizationErrorKind(ck::UffOptimizationErrorKind);
#[wasm_bindgen]
impl UffOptimizationErrorKind {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> UffOptimizationErrorVariant {
        match self.0 {
            ck::UffOptimizationErrorKind::MissingConformer { .. } => {
                UffOptimizationErrorVariant::MissingConformer
            }
            ck::UffOptimizationErrorKind::Rings => UffOptimizationErrorVariant::Rings,
            ck::UffOptimizationErrorKind::Optimization => UffOptimizationErrorVariant::Optimization,
            ck::UffOptimizationErrorKind::ConformerOptimization => {
                UffOptimizationErrorVariant::ConformerOptimization
            }
            ck::UffOptimizationErrorKind::Evaluation => UffOptimizationErrorVariant::Evaluation,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn requested(&self) -> Option<usize> {
        match self.0 {
            ck::UffOptimizationErrorKind::MissingConformer { requested } => requested,
            _ => None,
        }
    }
}
error_projection!(
    UffOptimizationError,
    ck::UffOptimizationError,
    {
        #[wasm_bindgen(getter)]
        pub fn kind(&self) -> UffOptimizationErrorKind {
            UffOptimizationErrorKind(self.inner.kind())
        }
    },
    source = ExecutionCause,
    cause = cause::<ck::UffOptimizationError>
);
