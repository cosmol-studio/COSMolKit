//! Exact native MCS search failures and transparent sources.
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
crate::errors::enum_error!(SearchMcsError,SearchMcsErrorKind,ck::SearchMcsError,{
TooFewInputs=>ck::SearchMcsError::TooFewInputs{..},
ThresholdAboveOne=>ck::SearchMcsError::ThresholdAboveOne,
AtomOutOfRange=>ck::SearchMcsError::AtomOutOfRange{..},
MissingConformer=>ck::SearchMcsError::MissingConformer{..},
CoordinateOutOfRange=>ck::SearchMcsError::CoordinateOutOfRange{..},
MissingRingInfo=>ck::SearchMcsError::MissingRingInfo{..},
MissingValence=>ck::SearchMcsError::MissingValence{..},
ValenceOutOfRange=>ck::SearchMcsError::ValenceOutOfRange{..},
BondOutOfRange=>ck::SearchMcsError::BondOutOfRange{..},
SeedExcludedBondOutOfRange=>ck::SearchMcsError::SeedExcludedBondOutOfRange{..},
SeedBondAlreadyExcluded=>ck::SearchMcsError::SeedBondAlreadyExcluded{..},
SeedAtomMappingMissing=>ck::SearchMcsError::SeedAtomMappingMissing{..},
TargetAtomMappingMissing=>ck::SearchMcsError::TargetAtomMappingMissing{..},
SeedAtomNotBondEndpoint=>ck::SearchMcsError::SeedAtomNotBondEndpoint{..},
SeedNewAtomMissing=>ck::SearchMcsError::SeedNewAtomMissing{..},
SeedRemainingBondUnderflow=>ck::SearchMcsError::SeedRemainingBondUnderflow{..},
SeedRemainingAtomUnderflow=>ck::SearchMcsError::SeedRemainingAtomUnderflow{..}
},{
#[wasm_bindgen(getter)] pub fn count(&self)->Option<usize>{match self.inner.as_ref(){ck::SearchMcsError::TooFewInputs{count,..}|ck::SearchMcsError::SeedExcludedBondOutOfRange{count,..}=>Some(*count),_=>None}}
#[wasm_bindgen(getter)] pub fn side(&self)->Option<String>{match self.inner.as_ref(){ck::SearchMcsError::AtomOutOfRange{side,..}|ck::SearchMcsError::MissingConformer{side,..}|ck::SearchMcsError::CoordinateOutOfRange{side,..}|ck::SearchMcsError::MissingRingInfo{side,..}|ck::SearchMcsError::MissingValence{side,..}|ck::SearchMcsError::ValenceOutOfRange{side,..}|ck::SearchMcsError::BondOutOfRange{side,..}=>Some(side.to_string()),_=>None}}
#[wasm_bindgen(getter)] pub fn atom(&self)->Option<usize>{match self.inner.as_ref(){ck::SearchMcsError::AtomOutOfRange{atom,..}|ck::SearchMcsError::CoordinateOutOfRange{atom,..}|ck::SearchMcsError::ValenceOutOfRange{atom,..}|ck::SearchMcsError::SeedAtomMappingMissing{atom,..}|ck::SearchMcsError::TargetAtomMappingMissing{atom,..}|ck::SearchMcsError::SeedAtomNotBondEndpoint{atom,..}|ck::SearchMcsError::SeedNewAtomMissing{atom,..}=>Some(*atom),_=>None}}
#[wasm_bindgen(getter)] pub fn bond(&self)->Option<usize>{match self.inner.as_ref(){ck::SearchMcsError::BondOutOfRange{bond,..}|ck::SearchMcsError::SeedExcludedBondOutOfRange{bond,..}|ck::SearchMcsError::SeedBondAlreadyExcluded{bond,..}|ck::SearchMcsError::SeedAtomNotBondEndpoint{bond,..}=>Some(*bond),_=>None}}
#[wasm_bindgen(getter)] pub fn remaining(&self)->Option<usize>{match self.inner.as_ref(){ck::SearchMcsError::SeedRemainingBondUnderflow{remaining,..}|ck::SearchMcsError::SeedRemainingAtomUnderflow{remaining,..}=>Some(*remaining),_=>None}}
#[wasm_bindgen(getter)] pub fn frontier(&self)->Option<usize>{match self.inner.as_ref(){ck::SearchMcsError::SeedRemainingBondUnderflow{frontier,..}=>Some(*frontier),_=>None}}
#[wasm_bindgen(getter)] pub fn added(&self)->Option<usize>{match self.inner.as_ref(){ck::SearchMcsError::SeedRemainingAtomUnderflow{added,..}=>Some(*added),_=>None}}
});
crate::errors::enum_error!(McsProgressError,McsProgressErrorKind,ck::McsProgressError,{Callback=>ck::McsProgressError::Callback{..}},{#[wasm_bindgen(getter)]pub fn message(&self)->String{match self.inner.as_ref(){ck::McsProgressError::Callback{message}=>message.clone()}}});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "PropertyUIntReadError | AtomPropertyError")]
    #[derive(Clone)]
    pub type McsCause;
}
fn cause(owner: &ReadonlyValue<ck::McsError>) -> Option<McsCause> {
    match owner.as_ref() {
        ck::McsError::StereoOrder(_) => {
            crate::stereo_order_errors::StereoOrderError::native_cause(&owner.child(|v| match v {
                ck::McsError::StereoOrder(v) => v,
                _ => unreachable!(),
            }))
            .map(JsCast::unchecked_into)
        }
        _ => None,
    }
}
crate::errors::enum_error!(McsError,McsErrorKind,ck::McsError,{
TargetTableCount=>ck::McsError::TargetTableCount{..},
ThresholdCountOutOfRange=>ck::McsError::ThresholdCountOutOfRange{..},
MatchTableOutOfRange=>ck::McsError::MatchTableOutOfRange{..},
RingMembershipMissing=>ck::McsError::RingMembershipMissing{..},
RingOutOfRange=>ck::McsError::RingOutOfRange{..},
MappedBondMissing=>ck::McsError::MappedBondMissing{..},
StereoOrder=>ck::McsError::StereoOrder(_),
TooManyNewBonds=>ck::McsError::TooManyNewBonds{..},
InitialSeedParse=>ck::McsError::InitialSeedParse{..},
InitialSeedMatch=>ck::McsError::InitialSeedMatch{..},
InitialSeedBondMissing=>ck::McsError::InitialSeedBondMissing{..},
ResultValueOutOfRange=>ck::McsError::ResultValueOutOfRange{..},
ResultQueryGraph=>ck::McsError::ResultQueryGraph{..},
ResultSmarts=>ck::McsError::ResultSmarts{..},
SeedReconstructionBondMissing=>ck::McsError::SeedReconstructionBondMissing{..},
ResultContextMissing=>ck::McsError::ResultContextMissing,
Progress=>ck::McsError::Progress(_),
State=>ck::McsError::State(_)
},{
#[wasm_bindgen(getter)] pub fn targets(&self)->Option<usize>{match self.inner.as_ref(){ck::McsError::TargetTableCount{targets,..}|ck::McsError::ThresholdCountOutOfRange{targets,..}=>Some(*targets),_=>None}}
#[wasm_bindgen(getter)] pub fn tables(&self)->Option<usize>{match self.inner.as_ref(){ck::McsError::TargetTableCount{tables,..}=>Some(*tables),_=>None}}
#[wasm_bindgen(getter)] pub fn threshold(&self)->Option<usize>{match self.inner.as_ref(){ck::McsError::ThresholdCountOutOfRange{threshold,..}=>Some(*threshold),_=>None}}
#[wasm_bindgen(getter)] pub fn kind(&self)->Option<String>{match self.inner.as_ref(){ck::McsError::MatchTableOutOfRange{kind,..}|ck::McsError::ResultValueOutOfRange{kind,..}=>Some(kind.to_string()),_=>None}}
#[wasm_bindgen(getter)] pub fn row(&self)->Option<usize>{match self.inner.as_ref(){ck::McsError::MatchTableOutOfRange{row,..}=>Some(*row),_=>None}}
#[wasm_bindgen(getter)] pub fn column(&self)->Option<usize>{match self.inner.as_ref(){ck::McsError::MatchTableOutOfRange{column,..}=>Some(*column),_=>None}}
#[wasm_bindgen(getter)] pub fn bond(&self)->Option<usize>{match self.inner.as_ref(){ck::McsError::RingMembershipMissing{bond,..}|ck::McsError::InitialSeedBondMissing{bond,..}=>Some(*bond),_=>None}}
#[wasm_bindgen(getter)] pub fn ring(&self)->Option<usize>{match self.inner.as_ref(){ck::McsError::RingOutOfRange{ring,..}=>Some(*ring),_=>None}}
#[wasm_bindgen(getter)] pub fn side(&self)->Option<String>{match self.inner.as_ref(){ck::McsError::MappedBondMissing{side,..}=>Some(side.to_string()),_=>None}}
#[wasm_bindgen(getter)] pub fn begin(&self)->Option<usize>{match self.inner.as_ref(){ck::McsError::MappedBondMissing{begin,..}|ck::McsError::SeedReconstructionBondMissing{begin,..}=>Some(*begin),_=>None}}
#[wasm_bindgen(getter)] pub fn end(&self)->Option<usize>{match self.inner.as_ref(){ck::McsError::MappedBondMissing{end,..}|ck::McsError::SeedReconstructionBondMissing{end,..}=>Some(*end),_=>None}}
#[wasm_bindgen(getter)] pub fn count(&self)->Option<usize>{match self.inner.as_ref(){ck::McsError::TooManyNewBonds{count,..}=>Some(*count),_=>None}}
#[wasm_bindgen(getter)] pub fn message(&self)->Option<String>{match self.inner.as_ref(){ck::McsError::InitialSeedParse{message,..}|ck::McsError::InitialSeedMatch{message,..}|ck::McsError::ResultQueryGraph{message,..}|ck::McsError::ResultSmarts{message,..}=>Some(message.to_string()),_=>None}}
#[wasm_bindgen(getter)] pub fn value(&self)->Option<usize>{match self.inner.as_ref(){ck::McsError::ResultValueOutOfRange{value,..}=>Some(*value),_=>None}}
#[wasm_bindgen(getter,js_name=stereoOrder)]pub fn stereoorder(&self)->Option<StereoOrderErrorPayload>{match self.inner.as_ref(){ck::McsError::StereoOrder(_)=>Some(crate::stereo_order_errors::StereoOrderError::from_owner(self.inner.child(|e|match e{ck::McsError::StereoOrder(v)=>v,_=>unreachable!()})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=progress)]pub fn progress(&self)->Option<McsProgressErrorPayload>{match self.inner.as_ref(){ck::McsError::Progress(_)=>Some(McsProgressError::from_owner(self.inner.child(|e|match e{ck::McsError::Progress(v)=>v,_=>unreachable!()})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=state)]pub fn state(&self)->Option<SearchMcsErrorPayload>{match self.inner.as_ref(){ck::McsError::State(_)=>Some(SearchMcsError::from_owner(self.inner.child(|e|match e{ck::McsError::State(v)=>v,_=>unreachable!()})).unchecked_into()),_=>None}}
}, source=McsCause,cause=cause);

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "StereoOrderError")]
    pub type StereoOrderErrorPayload;
}

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "McsProgressError")]
    pub type McsProgressErrorPayload;
}

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "SearchMcsError")]
    pub type SearchMcsErrorPayload;
}
