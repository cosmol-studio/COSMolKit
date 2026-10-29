//! Fragment failures preserve the native transparent-source graph.
use crate::{
    errors::{enum_error, error_projection, exception},
    model_errors::*,
    property_errors::*,
    topology_errors::TopologyEditError,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | TemplateAttachmentOrderError | UIntLexicalReadError | BondValueError | MoleculePropertyError | AtomPropertyError | NonNullable<SanitizeError['cause']> | AtomPathSubsetCopyError | FullCopyComponentError | SanitizeError"
    )]
    #[derive(Clone)]
    pub type FragmentCause;
    #[wasm_bindgen(
        typescript_type = "BondValueError | PropertyUIntReadError | PropertyStringError | TopologyEditError | MappingValidationError | CoordinateValidationError | StereoGroupError | MoleculePropertyError | AtomPropertyError | SelectedBondCopyError | SelectedSubstanceGroupCopyError | TopologyValidationError | SanitizeError | PathError | AtomPathSubsetCopyError | FullCopyComponentError | OrderedFragmentBuildError"
    )]
    pub type FragmentPayload;
    #[wasm_bindgen(typescript_type = "MoleculeFragmentsFailure")]
    #[derive(Clone)]
    pub type FragmentsFailure;
}
fn bond_cause(e: &ck::SelectedBondCopyError) -> Option<FragmentCause> {
    match e {
        ck::SelectedBondCopyError::Stereo(v) => {
            Some(BondValueError::from_rust(v.clone()).unchecked_into())
        }
        _ => None,
    }
}
enum_error!(SelectedBondCopyError,SelectedBondCopyErrorKind,ck::SelectedBondCopyError,{Stereo=>ck::SelectedBondCopyError::Stereo(_),MissingEndpointMapping=>ck::SelectedBondCopyError::MissingEndpointMapping{..}},{
    #[wasm_bindgen(getter)]pub fn bond(&self)->Option<usize>{match *self.inner{ck::SelectedBondCopyError::MissingEndpointMapping{bond,..}=>Some(bond.index()),_=>None}}
    #[wasm_bindgen(getter)]pub fn atom(&self)->Option<usize>{match *self.inner{ck::SelectedBondCopyError::MissingEndpointMapping{atom,..}=>Some(atom.index()),_=>None}}
    #[wasm_bindgen(getter)]pub fn stereo(&self)->Option<FragmentPayload>{match self.inner.as_ref(){ck::SelectedBondCopyError::Stereo(v)=>Some(BondValueError::from_rust(v.clone()).unchecked_into()),_=>None}}
},source=FragmentCause,cause=bond_cause);
enum_error!(SelectedSubstanceGroupCopyError,SelectedSubstanceGroupCopyErrorKind,ck::SelectedSubstanceGroupCopyError,{
 MissingAtomMapping=>ck::SelectedSubstanceGroupCopyError::MissingAtomMapping{..},AtomOutOfRange=>ck::SelectedSubstanceGroupCopyError::AtomOutOfRange{..},MissingBondMapping=>ck::SelectedSubstanceGroupCopyError::MissingBondMapping{..},BondOutOfRange=>ck::SelectedSubstanceGroupCopyError::BondOutOfRange{..}
},{
 #[wasm_bindgen(getter)]pub fn group(&self)->usize{use ck::SelectedSubstanceGroupCopyError as E;match *self.inner{E::MissingAtomMapping{group,..}|E::AtomOutOfRange{group,..}|E::MissingBondMapping{group,..}|E::BondOutOfRange{group,..}=>group.index()}}
 #[wasm_bindgen(getter)]pub fn list(&self)->Option<crate::enums::SubstanceGroupAtomList>{use ck::SelectedSubstanceGroupCopyError as E;match *self.inner{E::MissingAtomMapping{list,..}|E::AtomOutOfRange{list,..}=>Some(list.into()),_=>None}}
 #[wasm_bindgen(getter)]pub fn atom(&self)->Option<usize>{match *self.inner{ck::SelectedSubstanceGroupCopyError::MissingAtomMapping{atom,..}=>Some(atom.index()),_=>None}}
 #[wasm_bindgen(getter)]pub fn bond(&self)->Option<usize>{match *self.inner{ck::SelectedSubstanceGroupCopyError::MissingBondMapping{bond,..}=>Some(bond.index()),_=>None}}
 #[wasm_bindgen(getter,js_name=sourceAtom)]pub fn source_atom(&self)->Option<usize>{match *self.inner{ck::SelectedSubstanceGroupCopyError::AtomOutOfRange{source_atom,..}=>Some(source_atom.index()),_=>None}}
 #[wasm_bindgen(getter,js_name=mappedAtom)]pub fn mapped_atom(&self)->Option<usize>{match *self.inner{ck::SelectedSubstanceGroupCopyError::AtomOutOfRange{mapped_atom,..}=>Some(mapped_atom.index()),_=>None}}
 #[wasm_bindgen(getter,js_name=atomCount)]pub fn atom_count(&self)->Option<usize>{match *self.inner{ck::SelectedSubstanceGroupCopyError::AtomOutOfRange{atom_count,..}=>Some(atom_count),_=>None}}
 #[wasm_bindgen(getter,js_name=sourceBond)]pub fn source_bond(&self)->Option<usize>{match *self.inner{ck::SelectedSubstanceGroupCopyError::BondOutOfRange{source_bond,..}=>Some(source_bond.index()),_=>None}}
 #[wasm_bindgen(getter,js_name=mappedBond)]pub fn mapped_bond(&self)->Option<usize>{match *self.inner{ck::SelectedSubstanceGroupCopyError::BondOutOfRange{mapped_bond,..}=>Some(mapped_bond.index()),_=>None}}
 #[wasm_bindgen(getter,js_name=bondCount)]pub fn bond_count(&self)->Option<usize>{match *self.inner{ck::SelectedSubstanceGroupCopyError::BondOutOfRange{bond_count,..}=>Some(bond_count),_=>None}}
});
fn full_cause(e: &ck::FullCopyComponentError) -> Option<FragmentCause> {
    use ck::FullCopyComponentError as E;
    match e {
        E::SourceUInt(ck::PropertyUIntReadError::Lexical { source, .. }) => {
            Some(UIntLexicalReadError::from_rust(source.clone()).unchecked_into())
        }
        E::TopologyEdit(v) => crate::topology_errors::cause(v).map(JsCast::unchecked_into),
        _ => None,
    }
}
enum_error!(FullCopyComponentError,FullCopyComponentErrorKind,ck::FullCopyComponentError,{
SourceUInt=>ck::FullCopyComponentError::SourceUInt(_),SourceText=>ck::FullCopyComponentError::SourceText(_),TopologyEdit=>ck::FullCopyComponentError::TopologyEdit(_),MappingValidation=>ck::FullCopyComponentError::MappingValidation(_),CoordinateValidation=>ck::FullCopyComponentError::CoordinateValidation(_),AtomMaskSize=>ck::FullCopyComponentError::AtomMaskSize{..},AddedAtomRow=>ck::FullCopyComponentError::AddedAtomRow{..},MissingCoordinateRow=>ck::FullCopyComponentError::MissingCoordinateRow{..}
},{
 #[wasm_bindgen(getter)]pub fn source(&self)->Option<FragmentPayload>{use ck::FullCopyComponentError as E;Some(match self.inner.as_ref(){E::SourceUInt(v)=>PropertyUIntReadError::from_rust(v.clone()),E::SourceText(v)=>PropertyStringError::from_rust(*v),E::TopologyEdit(v)=>TopologyEditError::from_rust(v.clone()),E::MappingValidation(v)=>MappingValidationError::from_rust(v.clone()),E::CoordinateValidation(v)=>CoordinateValidationError::from_rust(v.clone()),_=>return None}.unchecked_into())}
 #[wasm_bindgen(getter)]pub fn expected(&self)->Option<usize>{match *self.inner{ck::FullCopyComponentError::AtomMaskSize{expected,..}=>Some(expected),_=>None}}
 #[wasm_bindgen(getter)]pub fn actual(&self)->Option<usize>{match *self.inner{ck::FullCopyComponentError::AtomMaskSize{actual,..}=>Some(actual),_=>None}}
 #[wasm_bindgen(getter)]pub fn row(&self)->Option<usize>{match *self.inner{ck::FullCopyComponentError::AddedAtomRow{row}=>Some(row),_=>None}}
 #[wasm_bindgen(getter,js_name=sourceRow)]pub fn source_row(&self)->Option<usize>{match *self.inner{ck::FullCopyComponentError::MissingCoordinateRow{source_row,..}=>Some(source_row),_=>None}}
 #[wasm_bindgen(getter,js_name=conformerId)]pub fn conformer_id(&self)->Option<usize>{match *self.inner{ck::FullCopyComponentError::MissingCoordinateRow{conformer_id,..}=>Some(conformer_id),_=>None}}
},source=FragmentCause,cause=full_cause);
fn subset_cause(e: &ck::AtomPathSubsetCopyError) -> Option<FragmentCause> {
    use ck::AtomPathSubsetCopyError as E;
    match e {
        E::StereoGroup(v) => Some(StereoGroupError::from_rust(*v).unchecked_into()),
        E::ComputedProperties(v) => crate::topology_errors::cause(v).map(JsCast::unchecked_into),
        E::MoleculeProperty(v) => {
            Some(MoleculePropertyError::from_rust(v.clone()).unchecked_into())
        }
        E::BondProperty(v) => Some(BondValueError::from_rust(v.clone()).unchecked_into()),
        E::AtomProperty(v) => Some(AtomPropertyError::from_rust(v.clone()).unchecked_into()),
        E::BondCopy(v) => bond_cause(v),
        E::TopologyValidation(v) => {
            crate::model_errors::topology_cause(v).map(JsCast::unchecked_into)
        }
        E::Sanitize(v) => crate::sanitize_errors::cause(v).map(JsCast::unchecked_into),
        E::SubstanceGroupCopy(_) | E::MappingValidation(_) => None,
    }
}
enum_error!(AtomPathSubsetCopyError,AtomPathSubsetCopyErrorKind,ck::AtomPathSubsetCopyError,{
StereoGroup=>ck::AtomPathSubsetCopyError::StereoGroup(_),ComputedProperties=>ck::AtomPathSubsetCopyError::ComputedProperties(_),MoleculeProperty=>ck::AtomPathSubsetCopyError::MoleculeProperty(_),BondProperty=>ck::AtomPathSubsetCopyError::BondProperty(_),AtomProperty=>ck::AtomPathSubsetCopyError::AtomProperty(_),BondCopy=>ck::AtomPathSubsetCopyError::BondCopy(_),SubstanceGroupCopy=>ck::AtomPathSubsetCopyError::SubstanceGroupCopy(_),TopologyValidation=>ck::AtomPathSubsetCopyError::TopologyValidation(_),MappingValidation=>ck::AtomPathSubsetCopyError::MappingValidation(_),Sanitize=>ck::AtomPathSubsetCopyError::Sanitize(_)
},{
#[wasm_bindgen(getter)]pub fn source(&self)->FragmentPayload{use ck::AtomPathSubsetCopyError as E;match self.inner.as_ref(){E::StereoGroup(v)=>StereoGroupError::from_rust(*v),E::ComputedProperties(v)=>TopologyEditError::from_rust(v.clone()),E::MoleculeProperty(v)=>MoleculePropertyError::from_rust(v.clone()),E::BondProperty(v)=>BondValueError::from_rust(v.clone()),E::AtomProperty(v)=>AtomPropertyError::from_rust(v.clone()),E::BondCopy(v)=>SelectedBondCopyError::from_rust(v.clone()),E::SubstanceGroupCopy(v)=>SelectedSubstanceGroupCopyError::from_rust(v.clone()),E::TopologyValidation(v)=>TopologyValidationError::from_rust(v.clone()),E::MappingValidation(v)=>MappingValidationError::from_rust(v.clone()),E::Sanitize(v)=>crate::sanitize_errors::SanitizeError::from_rust(v.clone())}.unchecked_into()}
},source=FragmentCause,cause=subset_cause);
fn build_cause(e: &ck::OrderedFragmentBuildError) -> Option<FragmentCause> {
    use ck::OrderedFragmentBuildError as E;
    match e {
        E::ComponentLabels(v) => crate::path_errors::source(v).map(JsCast::unchecked_into),
        E::FastSubset { source, .. } => {
            Some(AtomPathSubsetCopyError::from_rust(source.clone()).unchecked_into())
        }
        E::SlowFullCopy { source, .. } => {
            Some(FullCopyComponentError::from_rust(source.clone()).unchecked_into())
        }
    }
}
enum_error!(OrderedFragmentBuildError,OrderedFragmentBuildErrorKind,ck::OrderedFragmentBuildError,{ComponentLabels=>ck::OrderedFragmentBuildError::ComponentLabels(_),FastSubset=>ck::OrderedFragmentBuildError::FastSubset{..},SlowFullCopy=>ck::OrderedFragmentBuildError::SlowFullCopy{..}},{
#[wasm_bindgen(getter,js_name=componentIndex)]pub fn component_index(&self)->Option<usize>{match *self.inner{ck::OrderedFragmentBuildError::FastSubset{component_index,..}|ck::OrderedFragmentBuildError::SlowFullCopy{component_index,..}=>Some(component_index),_=>None}}
#[wasm_bindgen(getter)]pub fn source(&self)->FragmentPayload{match self.inner.as_ref(){ck::OrderedFragmentBuildError::ComponentLabels(v)=>crate::path_errors::PathError::from_rust(v.clone()),ck::OrderedFragmentBuildError::FastSubset{source,..}=>AtomPathSubsetCopyError::from_rust(source.clone()),ck::OrderedFragmentBuildError::SlowFullCopy{source,..}=>FullCopyComponentError::from_rust(source.clone())}.unchecked_into()}
},source=FragmentCause,cause=build_cause);
fn failure_cause(e: &ck::MoleculeFragmentsFailure) -> Option<FragmentCause> {
    match e {
        ck::MoleculeFragmentsFailure::Build(v) => build_cause(v),
        ck::MoleculeFragmentsFailure::FinalComputedProperties { source, .. } => {
            Some(MoleculePropertyError::from_rust(source.clone()).unchecked_into())
        }
        ck::MoleculeFragmentsFailure::FinalSanitize { source, .. } => {
            Some(crate::sanitize_errors::SanitizeError::from_rust(source.clone()).unchecked_into())
        }
    }
}
enum_error!(MoleculeFragmentsFailure,MoleculeFragmentsFailureKind,ck::MoleculeFragmentsFailure,{Build=>ck::MoleculeFragmentsFailure::Build(_),FinalComputedProperties=>ck::MoleculeFragmentsFailure::FinalComputedProperties{..},FinalSanitize=>ck::MoleculeFragmentsFailure::FinalSanitize{..}},{
#[wasm_bindgen(getter,js_name=componentIndex)]pub fn component_index(&self)->Option<usize>{match *self.inner{ck::MoleculeFragmentsFailure::FinalComputedProperties{component_index,..}|ck::MoleculeFragmentsFailure::FinalSanitize{component_index,..}=>Some(component_index),_=>None}}
#[wasm_bindgen(getter)]pub fn source(&self)->FragmentPayload{match self.inner.as_ref(){ck::MoleculeFragmentsFailure::Build(v)=>OrderedFragmentBuildError::from_rust(v.clone()),ck::MoleculeFragmentsFailure::FinalComputedProperties{source,..}=>MoleculePropertyError::from_rust(source.clone()),ck::MoleculeFragmentsFailure::FinalSanitize{source,..}=>crate::sanitize_errors::SanitizeError::from_rust(source.clone())}.unchecked_into()}
},source=FragmentCause,cause=failure_cause);
error_projection!(
    MoleculeFragmentsError,
    ck::MoleculeFragmentsError,
    {
        #[wasm_bindgen(getter,js_name=componentIndex)]
        pub fn component_index(&self) -> Option<usize> {
            self.inner.component_index()
        }
        #[wasm_bindgen(getter)]
        pub fn failure(&self) -> FragmentsFailure {
            MoleculeFragmentsFailure::from_rust(self.inner.failure().clone()).unchecked_into()
        }
    },
    source = FragmentsFailure,
    cause = |e: &ck::MoleculeFragmentsError| Some(
        MoleculeFragmentsFailure::from_rust(e.failure().clone()).unchecked_into()
    )
);
