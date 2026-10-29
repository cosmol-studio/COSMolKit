//! Exact native UFF preparation, query and parameter failures.
use crate::{
    errors::child,
    errors::{enum_error, error_projection, exception},
    forcefield_kernel_errors::*,
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use std::error::Error;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "ValenceError | PropertyStringError | UffBondMathError")]
    #[derive(Clone)]
    pub type TypingCause;
    #[wasm_bindgen(
        typescript_type = "UffTorsionBondQueryError | UffInversionContributionError | ValenceError | ForceFieldKernelError | TopologyValidationError"
    )]
    #[derive(Clone)]
    pub type BuilderCause;
    #[wasm_bindgen(typescript_type = "number | Hybridization")]
    pub type BuilderActual;
    #[wasm_bindgen(typescript_type = "SmartsParseError | QueryCompileError")]
    #[derive(Clone)]
    pub type DefaultQueryCause;
    #[wasm_bindgen(
        typescript_type = "UffDefaultTorsionQueryError | SmartsParseError | QueryCompileError | MatchError"
    )]
    #[derive(Clone)]
    pub type TorsionQueryCause;
    #[wasm_bindgen(typescript_type = "UffBuilderError | UffParamError | UffTypingError")]
    #[derive(Clone)]
    pub type ParameterCause;
}
crate::enums::enum_value!(UffParameterErrorKind,"UffParameterErrorKind"{Preparation="preparation",ParameterTable="parameter_table",Typing="typing"});
enum_error!(UffTypingError,UffTypingErrorKind,ck::UffTypingError,{CoreValence=>ck::UffTypingError::CoreValence(_),CorePropertyString=>ck::UffTypingError::CorePropertyString(_),BondMath=>ck::UffTypingError::BondMath(_),PreparedStateLength=>ck::UffTypingError::PreparedStateLength{..},AtomIdOutOfBounds=>ck::UffTypingError::AtomIdOutOfBounds{..},AtomIndexOutOfBounds=>ck::UffTypingError::AtomIndexOutOfBounds{..},BondIdOutOfBounds=>ck::UffTypingError::BondIdOutOfBounds{..},CenterBondMissingAfterSourceSuccess=>ck::UffTypingError::CenterBondMissingAfterSourceSuccess{..},CentralParameterSlotMissingAfterSourceSuccess=>ck::UffTypingError::CentralParameterSlotMissingAfterSourceSuccess{..}},{
#[wasm_bindgen(getter,js_name=input)]pub fn input(&self)->Option<UffTypingInput>{match self.inner.as_ref(){ck::UffTypingError::PreparedStateLength{input,..}=>Some((*input).into()),_=>None}}
#[wasm_bindgen(getter,js_name=expected)]pub fn expected(&self)->Option<usize>{match self.inner.as_ref(){ck::UffTypingError::PreparedStateLength{expected,..}=>Some(*expected),_=>None}}
#[wasm_bindgen(getter,js_name=actual)]pub fn actual(&self)->Option<usize>{match self.inner.as_ref(){ck::UffTypingError::PreparedStateLength{actual,..}=>Some(*actual),_=>None}}
#[wasm_bindgen(getter,js_name=atomId)]pub fn atom_id(&self)->Option<usize>{match self.inner.as_ref(){ck::UffTypingError::AtomIdOutOfBounds{atom_id,..}=>Some(atom_id.index()),_=>None}}
#[wasm_bindgen(getter,js_name=atomCount)]pub fn atom_count(&self)->Option<usize>{match self.inner.as_ref(){ck::UffTypingError::AtomIdOutOfBounds{atom_count,..}=>Some(*atom_count),ck::UffTypingError::AtomIndexOutOfBounds{atom_count,..}=>Some(*atom_count),_=>None}}
#[wasm_bindgen(getter,js_name=index)]pub fn index(&self)->Option<usize>{match self.inner.as_ref(){ck::UffTypingError::AtomIndexOutOfBounds{index,..}=>Some(*index),_=>None}}
#[wasm_bindgen(getter,js_name=bondId)]pub fn bond_id(&self)->Option<usize>{match self.inner.as_ref(){ck::UffTypingError::BondIdOutOfBounds{bond_id,..}=>Some(bond_id.index()),_=>None}}
#[wasm_bindgen(getter,js_name=bondCount)]pub fn bond_count(&self)->Option<usize>{match self.inner.as_ref(){ck::UffTypingError::BondIdOutOfBounds{bond_count,..}=>Some(*bond_count),_=>None}}
#[wasm_bindgen(getter,js_name=idx2)]pub fn idx2(&self)->Option<usize>{match self.inner.as_ref(){ck::UffTypingError::CenterBondMissingAfterSourceSuccess{idx2,..}=>Some(*idx2),_=>None}}
#[wasm_bindgen(getter,js_name=idx3)]pub fn idx3(&self)->Option<usize>{match self.inner.as_ref(){ck::UffTypingError::CenterBondMissingAfterSourceSuccess{idx3,..}=>Some(*idx3),_=>None}}
#[wasm_bindgen(getter,js_name=atomIndex)]pub fn atom_index(&self)->Option<usize>{match self.inner.as_ref(){ck::UffTypingError::CentralParameterSlotMissingAfterSourceSuccess{atom_index,..}=>Some(*atom_index),_=>None}}
},source=TypingCause,cause=|owner:&ReadonlyValue<ck::UffTypingError>|match owner.as_ref(){ck::UffTypingError::CoreValence(_)=>Some(child!(owner,ck::UffTypingError::CoreValence,crate::valence_errors::ValenceError)),ck::UffTypingError::CorePropertyString(_)=>Some(child!(owner,ck::UffTypingError::CorePropertyString,crate::property_errors::PropertyStringError)),ck::UffTypingError::BondMath(_)=>Some(child!(owner,ck::UffTypingError::BondMath,UffBondMathError)),_=>None});
enum_error!(UffBuilderError,UffBuilderErrorKind,ck::UffBuilderError,{NeighborMatrixIndexOutOfRange=>ck::UffBuilderError::NeighborMatrixIndexOutOfRange{..},NeighborMatrixIndexOverflow=>ck::UffBuilderError::NeighborMatrixIndexOverflow{..},NeighborMatrixAllocationFailed=>ck::UffBuilderError::NeighborMatrixAllocationFailed{..},NeighborMatrixCellPositionOverflow=>ck::UffBuilderError::NeighborMatrixCellPositionOverflow{..},NeighborMatrixStorageOutOfRange=>ck::UffBuilderError::NeighborMatrixStorageOutOfRange{..},ValenceAssignmentLengthMismatch=>ck::UffBuilderError::ValenceAssignmentLengthMismatch{..},SourceValencePrecondition=>ck::UffBuilderError::SourceValencePrecondition{..},SourceValenceOutOfRange=>ck::UffBuilderError::SourceValenceOutOfRange{..},ParamsLengthMismatch=>ck::UffBuilderError::ParamsLengthMismatch{..},SelectedThreeDimensionalConformerNotFound=>ck::UffBuilderError::SelectedThreeDimensionalConformerNotFound{..},SelectedConformerCoordinateCountMismatch=>ck::UffBuilderError::SelectedConformerCoordinateCountMismatch{..},SourceBondIndexOverflow=>ck::UffBuilderError::SourceBondIndexOverflow{..},SourceTorsionIndexOverflow=>ck::UffBuilderError::SourceTorsionIndexOverflow{..},SourceAngleIndexOverflow=>ck::UffBuilderError::SourceAngleIndexOverflow{..},SourceInversionIndexOverflow=>ck::UffBuilderError::SourceInversionIndexOverflow{..},SourceDirectionVectorBelowTolerance=>ck::UffBuilderError::SourceDirectionVectorBelowTolerance{..},SourceTbpAtomPrecondition=>ck::UffBuilderError::SourceTbpAtomPrecondition{..},SourceTbpHybridizationPrecondition=>ck::UffBuilderError::SourceTbpHybridizationPrecondition{..},SourceTbpDegreePrecondition=>ck::UffBuilderError::SourceTbpDegreePrecondition{..},SourceTbpAxialBondNotFound=>ck::UffBuilderError::SourceTbpAxialBondNotFound{..},SourceTbpEquatorialBondNotFound=>ck::UffBuilderError::SourceTbpEquatorialBondNotFound{..},SourceTbpCenterParamsMissing=>ck::UffBuilderError::SourceTbpCenterParamsMissing{..},TorsionBondQuery=>ck::UffBuilderError::TorsionBondQuery(_),InversionContribution=>ck::UffBuilderError::InversionContribution(_),Valence=>ck::UffBuilderError::Valence(_),ForceFieldKernel=>ck::UffBuilderError::ForceFieldKernel(_),TopologyValidation=>ck::UffBuilderError::TopologyValidation(_)},{
#[wasm_bindgen(getter,js_name=nAtoms)]pub fn n_atoms(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::NeighborMatrixIndexOutOfRange{n_atoms,..}=>Some(*n_atoms),ck::UffBuilderError::NeighborMatrixIndexOverflow{n_atoms,..}=>Some(*n_atoms),_=>None}}
#[wasm_bindgen(getter,js_name=i)]pub fn i(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::NeighborMatrixIndexOutOfRange{i,..}=>Some(*i),_=>None}}
#[wasm_bindgen(getter,js_name=j)]pub fn j(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::NeighborMatrixIndexOutOfRange{j,..}=>Some(*j),_=>None}}
#[wasm_bindgen(getter,js_name=byteCount)]pub fn byte_count(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::NeighborMatrixAllocationFailed{byte_count,..}=>Some(*byte_count),_=>None}}
#[wasm_bindgen(getter,js_name=position)]pub fn position(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::NeighborMatrixCellPositionOverflow{position,..}=>Some(*position),ck::UffBuilderError::NeighborMatrixStorageOutOfRange{position,..}=>Some(*position),_=>None}}
#[wasm_bindgen(getter,js_name=byteIndex)]pub fn byte_index(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::NeighborMatrixStorageOutOfRange{byte_index,..}=>Some(*byte_index),_=>None}}
#[wasm_bindgen(getter,js_name=storageLen)]pub fn storage_len(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::NeighborMatrixStorageOutOfRange{storage_len,..}=>Some(*storage_len),_=>None}}
#[wasm_bindgen(getter,js_name=field)]pub fn field(&self)->Option<UffPreparedValenceField>{match self.inner.as_ref(){ck::UffBuilderError::ValenceAssignmentLengthMismatch{field,..}=>Some((*field).into()),ck::UffBuilderError::SourceValencePrecondition{field,..}=>Some((*field).into()),ck::UffBuilderError::SourceValenceOutOfRange{field,..}=>Some((*field).into()),_=>None}}
#[wasm_bindgen(getter,js_name=expected)]pub fn expected(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::ValenceAssignmentLengthMismatch{expected,..}=>Some(*expected),_=>None}}
#[wasm_bindgen(getter,js_name=actual)]pub fn actual(&self)->Option<BuilderActual>{match self.inner.as_ref(){ck::UffBuilderError::ValenceAssignmentLengthMismatch{actual,..}=>Some(JsValue::from(*actual as u32).unchecked_into()),ck::UffBuilderError::SourceTbpHybridizationPrecondition{actual,..}=>Some(JsValue::from(crate::enums::Hybridization::from(*actual)).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=atomId)]pub fn atom_id(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::SourceValencePrecondition{atom_id,..}=>Some(atom_id.index()),ck::UffBuilderError::SourceValenceOutOfRange{atom_id,..}=>Some(atom_id.index()),_=>None}}
#[wasm_bindgen(getter,js_name=value)]pub fn value(&self)->Option<i32>{match self.inner.as_ref(){ck::UffBuilderError::SourceValencePrecondition{value,..}=>Some(*value),ck::UffBuilderError::SourceValenceOutOfRange{value,..}=>Some(*value),_=>None}}
#[wasm_bindgen(getter,js_name=atoms)]pub fn atoms(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::ParamsLengthMismatch{atoms,..}=>Some(*atoms),ck::UffBuilderError::SelectedConformerCoordinateCountMismatch{atoms,..}=>Some(*atoms),_=>None}}
#[wasm_bindgen(getter,js_name=params)]pub fn params(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::ParamsLengthMismatch{params,..}=>Some(*params),_=>None}}
#[wasm_bindgen(getter,js_name=conformerId)]pub fn conformer_id(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::SelectedThreeDimensionalConformerNotFound{conformer_id,..}=>Some(*conformer_id),ck::UffBuilderError::SelectedConformerCoordinateCountMismatch{conformer_id,..}=>Some(*conformer_id),_=>None}}
#[wasm_bindgen(getter,js_name=coordinates)]pub fn coordinates(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::SelectedConformerCoordinateCountMismatch{coordinates,..}=>Some(*coordinates),_=>None}}
#[wasm_bindgen(getter,js_name=atomIndex)]pub fn atom_index(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::SourceBondIndexOverflow{atom_index,..}=>Some(*atom_index),ck::UffBuilderError::SourceTorsionIndexOverflow{atom_index,..}=>Some(*atom_index),ck::UffBuilderError::SourceAngleIndexOverflow{atom_index,..}=>Some(*atom_index),ck::UffBuilderError::SourceInversionIndexOverflow{atom_index,..}=>Some(*atom_index),_=>None}}
#[wasm_bindgen(getter,js_name=centerAtomIndex)]pub fn center_atom_index(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::SourceDirectionVectorBelowTolerance{center_atom_index,..}=>Some(*center_atom_index),ck::UffBuilderError::SourceTbpAtomPrecondition{center_atom_index,..}=>Some(*center_atom_index),ck::UffBuilderError::SourceTbpHybridizationPrecondition{center_atom_index,..}=>Some(*center_atom_index),ck::UffBuilderError::SourceTbpDegreePrecondition{center_atom_index,..}=>Some(*center_atom_index),ck::UffBuilderError::SourceTbpAxialBondNotFound{center_atom_index,..}=>Some(*center_atom_index),ck::UffBuilderError::SourceTbpEquatorialBondNotFound{center_atom_index,..}=>Some(*center_atom_index),ck::UffBuilderError::SourceTbpCenterParamsMissing{center_atom_index,..}=>Some(*center_atom_index),_=>None}}
#[wasm_bindgen(getter,js_name=neighborAtomIndex)]pub fn neighbor_atom_index(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::SourceDirectionVectorBelowTolerance{neighbor_atom_index,..}=>Some(*neighbor_atom_index),_=>None}}
#[wasm_bindgen(getter,js_name=actualDegree)]pub fn actual_degree(&self)->Option<usize>{match self.inner.as_ref(){ck::UffBuilderError::SourceTbpDegreePrecondition{actual_degree,..}=>Some(*actual_degree),_=>None}}
#[wasm_bindgen(getter,js_name=role)]pub fn role(&self)->Option<u8>{match self.inner.as_ref(){ck::UffBuilderError::SourceTbpEquatorialBondNotFound{role,..}=>Some(*role),_=>None}}
},source=BuilderCause,cause=|owner:&ReadonlyValue<ck::UffBuilderError>|match owner.as_ref(){ck::UffBuilderError::TorsionBondQuery(_)=>Some(child!(owner,ck::UffBuilderError::TorsionBondQuery,UffTorsionBondQueryError)),ck::UffBuilderError::InversionContribution(_)=>Some(child!(owner,ck::UffBuilderError::InversionContribution,UffInversionContributionError)),ck::UffBuilderError::Valence(_)=>Some(child!(owner,ck::UffBuilderError::Valence,crate::valence_errors::ValenceError)),ck::UffBuilderError::ForceFieldKernel(_)=>Some(child!(owner,ck::UffBuilderError::ForceFieldKernel,ForceFieldKernelError)),ck::UffBuilderError::TopologyValidation(_)=>Some(child!(owner,ck::UffBuilderError::TopologyValidation,crate::model_errors::TopologyValidationError)),_=>None});
enum_error!(UffParamError,UffParamErrorKind,ck::UffParamError,{EmptyLine=>ck::UffParamError::EmptyLine{..},MalformedLine=>ck::UffParamError::MalformedLine{..},ParseFloat=>ck::UffParamError::ParseFloat{..}},{
#[wasm_bindgen(getter,js_name=lineNumber)]pub fn line_number(&self)->Option<usize>{match self.inner.as_ref(){ck::UffParamError::EmptyLine{line_number,..}=>Some(*line_number),ck::UffParamError::MalformedLine{line_number,..}=>Some(*line_number),ck::UffParamError::ParseFloat{line_number,..}=>Some(*line_number),_=>None}}
#[wasm_bindgen(getter,js_name=columnCount)]pub fn column_count(&self)->Option<usize>{match self.inner.as_ref(){ck::UffParamError::MalformedLine{column_count,..}=>Some(*column_count),_=>None}}
#[wasm_bindgen(getter,js_name=columnName)]pub fn column_name(&self)->Option<String>{match self.inner.as_ref(){ck::UffParamError::ParseFloat{column_name,..}=>Some((*column_name).into()),_=>None}}
#[wasm_bindgen(getter,js_name=value)]pub fn value(&self)->Option<String>{match self.inner.as_ref(){ck::UffParamError::ParseFloat{value,..}=>Some(value.clone()),_=>None}}
});
enum_error!(UffBondMathError,UffBondMathErrorKind,ck::UffBondMathError,{InvalidBondOrder=>ck::UffBondMathError::InvalidBondOrder{..}},{#[wasm_bindgen(getter,js_name=bondOrder)]pub fn bond_order(&self)->f64{let ck::UffBondMathError::InvalidBondOrder{bond_order}=*self.inner;bond_order}});
enum_error!(UffInversionContributionError,UffInversionContributionErrorKind,ck::UffInversionContributionError,{IndexOutOfRange=>ck::UffInversionContributionError::IndexOutOfRange{..}},{#[wasm_bindgen(getter)]pub fn argument(&self)->UffInversionIndexArgument{let ck::UffInversionContributionError::IndexOutOfRange{argument,..}=*self.inner;argument.into()}#[wasm_bindgen(getter)]pub fn index(&self)->u32{let ck::UffInversionContributionError::IndexOutOfRange{index,..}=*self.inner;index}#[wasm_bindgen(getter,js_name=upperBound)]pub fn upper_bound(&self)->usize{let ck::UffInversionContributionError::IndexOutOfRange{upper_bound,..}=*self.inner;upper_bound}});
enum_error!(UffDefaultTorsionQueryError,UffDefaultTorsionQueryErrorKind,ck::UffDefaultTorsionQueryError,{Parse=>ck::UffDefaultTorsionQueryError::Parse(_),Compile=>ck::UffDefaultTorsionQueryError::Compile(_)},{},source=DefaultQueryCause,cause=|owner:&ReadonlyValue<ck::UffDefaultTorsionQueryError>|Some(match owner.as_ref(){ck::UffDefaultTorsionQueryError::Parse(_)=>child!(owner,ck::UffDefaultTorsionQueryError::Parse,crate::smarts_parse_errors::SmartsParseError),ck::UffDefaultTorsionQueryError::Compile(_)=>child!(owner,ck::UffDefaultTorsionQueryError::Compile,crate::query_compile_error::QueryCompileError)}));
enum_error!(UffTorsionBondQueryError,UffTorsionBondQueryErrorKind,ck::UffTorsionBondQueryError,{DefaultQuery=>ck::UffTorsionBondQueryError::DefaultQuery(_),Parse=>ck::UffTorsionBondQueryError::Parse(_),Compile=>ck::UffTorsionBondQueryError::Compile(_),Match=>ck::UffTorsionBondQueryError::Match(_),MatchArity=>ck::UffTorsionBondQueryError::MatchArity{..},MissingMatchedBond=>ck::UffTorsionBondQueryError::MissingMatchedBond{..}},{
#[wasm_bindgen(getter)]pub fn smarts(&self)->Option<String>{match self.inner.as_ref(){ck::UffTorsionBondQueryError::MatchArity{smarts,..}=>Some(smarts.clone()),_=>None}}
#[wasm_bindgen(getter)]pub fn actual(&self)->Option<usize>{match self.inner.as_ref(){ck::UffTorsionBondQueryError::MatchArity{actual,..}=>Some(*actual),_=>None}}
#[wasm_bindgen(getter,js_name=beginAtomIndex)]pub fn begin_atom_index(&self)->Option<usize>{match self.inner.as_ref(){ck::UffTorsionBondQueryError::MissingMatchedBond{begin_atom_index,..}=>Some(*begin_atom_index),_=>None}}
#[wasm_bindgen(getter,js_name=endAtomIndex)]pub fn end_atom_index(&self)->Option<usize>{match self.inner.as_ref(){ck::UffTorsionBondQueryError::MissingMatchedBond{end_atom_index,..}=>Some(*end_atom_index),_=>None}}
},source=TorsionQueryCause,cause=|owner:&ReadonlyValue<ck::UffTorsionBondQueryError>|match owner.as_ref(){ck::UffTorsionBondQueryError::DefaultQuery(_)=>Some(child!(owner,ck::UffTorsionBondQueryError::DefaultQuery,UffDefaultTorsionQueryError)),ck::UffTorsionBondQueryError::Parse(_)=>Some(child!(owner,ck::UffTorsionBondQueryError::Parse,crate::smarts_parse_errors::SmartsParseError)),ck::UffTorsionBondQueryError::Compile(_)=>Some(child!(owner,ck::UffTorsionBondQueryError::Compile,crate::query_compile_error::QueryCompileError)),ck::UffTorsionBondQueryError::Match(_)=>Some(child!(owner,ck::UffTorsionBondQueryError::Match,crate::match_error::MatchError)),_=>None});
error_projection!(
    UffParameterError,
    ck::UffParameterError,
    {
        #[wasm_bindgen(getter)]
        pub fn kind(&self) -> UffParameterErrorKind {
            self.inner.kind().into()
        }
    },
    source = ParameterCause,
    cause = |owner: &ReadonlyValue<ck::UffParameterError>| Some(
        match owner.kind() {
            ck::UffParameterErrorKind::Preparation =>
                UffBuilderError::from_owner(owner.child(|e| {
                    e.source()
                        .unwrap()
                        .downcast_ref::<ck::UffBuilderError>()
                        .expect("native preparation cause")
                })),
            ck::UffParameterErrorKind::ParameterTable =>
                UffParamError::from_owner(owner.child(|e| {
                    e.source()
                        .unwrap()
                        .downcast_ref::<ck::UffParamError>()
                        .expect("native table cause")
                })),
            ck::UffParameterErrorKind::Typing => UffTypingError::from_owner(owner.child(|e| {
                e.source()
                    .unwrap()
                    .downcast_ref::<ck::UffTypingError>()
                    .expect("native typing cause")
            })),
        }
        .unchecked_into()
    )
);
