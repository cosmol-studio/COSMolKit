//! Native transparent reaction execution errors with all context retained.
use crate::errors::{enum_error, exception};
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "ReactionProductError | ReactionRunError | ReactionApplyError | StereoGroupError | CoordinateValidationError | TopologyEditError | AdjacencyError | RingFindingError | PropertyUIntReadError | ValenceError | StereoOrderError | DoubleBondStereoError | SmilesParseError | QueryAtomConversionError | MoleculePropertyError | AtomPropertyError | BondValueError | ReactionValidationError | SubstructMatchError | PropertyIntReadError | ReactionInitializationError | PropertyStringError | PropertyULongReadError | LegacyStereoError | PotentialStereoError | ParserStereoOrderError | PropertyValueError | ParserCarrierError | KekulizeError | SanitizeError | HydrogenError | SmilesStereoError | AtropisomerError | RequiredPropertyStringError | TopologyValidationError | TemplateAttachmentOrderError | RingDecomposerError | UIntLexicalReadError"
    )]
    #[derive(Clone)]
    pub type ReactionExecutionCause;
}
fn reaction_product_error_cause(
    owner: &ReadonlyValue<ck::ReactionProductError>,
) -> Option<ReactionExecutionCause> {
    match owner.as_ref() {
        ck::ReactionProductError::StereoGroup(_) => Some(
            crate::model_errors::StereoGroupError::from_owner(owner.child(|e| match e {
                ck::ReactionProductError::StereoGroup(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionProductError::Coordinate(_) => {
            crate::model_errors::CoordinateValidationError::native_cause(&owner.child(
                |e| match e {
                    ck::ReactionProductError::Coordinate(source) => source,
                    _ => unreachable!("immutable native reaction cause"),
                },
            ))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::TopologyEdit(_) => {
            crate::topology_errors::TopologyEditError::native_cause(&owner.child(|e| match e {
                ck::ReactionProductError::TopologyEdit(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::Adjacency(_) => {
            crate::model_errors::AdjacencyError::native_cause(&owner.child(|e| match e {
                ck::ReactionProductError::Adjacency(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::RingFinding(_) => {
            crate::ring_errors::RingFindingError::native_cause(&owner.child(|e| match e {
                ck::ReactionProductError::RingFinding(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::SourceUInt(_) => {
            crate::property_errors::PropertyUIntReadError::native_cause(&owner.child(|e| match e {
                ck::ReactionProductError::SourceUInt(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::Valence(_) => {
            crate::valence_errors::ValenceError::native_cause(&owner.child(|e| match e {
                ck::ReactionProductError::Valence(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::StereoOrder(_) => {
            crate::stereo_order_errors::StereoOrderError::native_cause(&owner.child(|e| match e {
                ck::ReactionProductError::StereoOrder(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::DoubleBondStereo(_) => {
            crate::double_stereo_errors::DoubleBondStereoError::native_cause(&owner.child(|e| {
                match e {
                    ck::ReactionProductError::DoubleBondStereo(source) => source,
                    _ => unreachable!("immutable native reaction cause"),
                }
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::CoordinateSelection(_) => {
            crate::smiles_parse_errors::SmilesParseError::native_cause(&owner.child(|e| match e {
                ck::ReactionProductError::CoordinateSelection(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::CarrierIdentity(_) => {
            crate::query_model_errors::QueryAtomConversionError::native_cause(&owner.child(|e| {
                match e {
                    ck::ReactionProductError::CarrierIdentity(source) => source,
                    _ => unreachable!("immutable native reaction cause"),
                }
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::MoleculeProperty(_) => {
            crate::model_errors::MoleculePropertyError::native_cause(&owner.child(|e| match e {
                ck::ReactionProductError::MoleculeProperty(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::AtomProperty(_) => {
            crate::model_errors::AtomPropertyError::native_cause(&owner.child(|e| match e {
                ck::ReactionProductError::AtomProperty(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::BondValue(_) => {
            crate::model_errors::BondValueError::native_cause(&owner.child(|e| match e {
                ck::ReactionProductError::BondValue(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::TemplateProperty(_) => {
            crate::reaction_errors::ReactionValidationError::native_cause(&owner.child(
                |e| match e {
                    ck::ReactionProductError::TemplateProperty(source) => source,
                    _ => unreachable!("immutable native reaction cause"),
                },
            ))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionProductError::StereoGetter { .. } => Some(
            crate::match_errors::SubstructMatchError::from_owner(owner.child(|e| match e {
                ck::ReactionProductError::StereoGetter { source, .. } => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionProductError::PropertyInt { .. } => Some(
            crate::property_errors::PropertyIntReadError::from_owner(owner.child(|e| match e {
                ck::ReactionProductError::PropertyInt { source, .. } => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionProductError::PropertyUInt { .. } => Some(
            crate::property_errors::PropertyUIntReadError::from_owner(owner.child(|e| match e {
                ck::ReactionProductError::PropertyUInt { source, .. } => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .unchecked_into(),
        ),
        _ => None,
    }
}
enum_error!(ReactionProductError,ReactionProductErrorKind,ck::ReactionProductError,{StereoGroup=>ck::ReactionProductError::StereoGroup(_),Coordinate=>ck::ReactionProductError::Coordinate(_),TopologyEdit=>ck::ReactionProductError::TopologyEdit(_),Adjacency=>ck::ReactionProductError::Adjacency(_),RingFinding=>ck::ReactionProductError::RingFinding(_),SourceUInt=>ck::ReactionProductError::SourceUInt(_),Valence=>ck::ReactionProductError::Valence(_),StereoOrder=>ck::ReactionProductError::StereoOrder(_),DoubleBondStereo=>ck::ReactionProductError::DoubleBondStereo(_),CoordinateSelection=>ck::ReactionProductError::CoordinateSelection(_),CarrierIdentity=>ck::ReactionProductError::CarrierIdentity(_),MoleculeProperty=>ck::ReactionProductError::MoleculeProperty(_),AtomProperty=>ck::ReactionProductError::AtomProperty(_),BondValue=>ck::ReactionProductError::BondValue(_),TemplateProperty=>ck::ReactionProductError::TemplateProperty(_),StereoGetter=>ck::ReactionProductError::StereoGetter{..},PropertyInt=>ck::ReactionProductError::PropertyInt{..},PropertyUInt=>ck::ReactionProductError::PropertyUInt{..},MissingProperty=>ck::ReactionProductError::MissingProperty{..},Invariant=>ck::ReactionProductError::Invariant{..},RowOverflow=>ck::ReactionProductError::RowOverflow{..},QueryReactantBond=>ck::ReactionProductError::QueryReactantBond{..}},{
#[wasm_bindgen(getter,js_name=atom)]pub fn atom(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionProductError::StereoGetter{atom,..}=>Some(atom.index()),ck::ReactionProductError::PropertyInt{atom,..}=>Some(atom.index()),ck::ReactionProductError::PropertyUInt{atom,..}=>Some(atom.index()),ck::ReactionProductError::MissingProperty{atom,..}=>Some(atom.index()),_=>None}}
#[wasm_bindgen(getter,js_name=key)]pub fn key(&self)->Option<String>{match self.inner.as_ref(){ck::ReactionProductError::PropertyInt{key,..}=>Some(key.to_string()),ck::ReactionProductError::PropertyUInt{key,..}=>Some(key.to_string()),ck::ReactionProductError::MissingProperty{key,..}=>Some(key.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=stage)]pub fn stage(&self)->Option<String>{match self.inner.as_ref(){ck::ReactionProductError::Invariant{stage,..}=>Some(stage.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=detail)]pub fn detail(&self)->Option<String>{match self.inner.as_ref(){ck::ReactionProductError::Invariant{detail,..}=>Some(detail.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=reactantAtom)]pub fn reactant_atom(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionProductError::Invariant{reactant_atom,..}=>*reactant_atom,_=>None}}
#[wasm_bindgen(getter,js_name=productAtom)]pub fn product_atom(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionProductError::Invariant{product_atom,..}=>*product_atom,_=>None}}
#[wasm_bindgen(getter,js_name=bond)]pub fn bond(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionProductError::Invariant{bond,..}=>bond.map(|v|v.index()),ck::ReactionProductError::QueryReactantBond{bond,..}=>Some(bond.index()),_=>None}}
#[wasm_bindgen(getter,js_name=kind)]pub fn kind(&self)->Option<String>{match self.inner.as_ref(){ck::ReactionProductError::RowOverflow{kind,..}=>Some(kind.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=index)]pub fn index(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionProductError::RowOverflow{index,..}=>Some(*index as usize),_=>None}}
#[wasm_bindgen(getter,js_name=stereoGroup)]pub fn stereo_group(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::StereoGroup(_)=>Some(crate::model_errors::StereoGroupError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::StereoGroup(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=coordinate)]pub fn coordinate(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::Coordinate(_)=>Some(crate::model_errors::CoordinateValidationError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::Coordinate(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=topologyEdit)]pub fn topology_edit(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::TopologyEdit(_)=>Some(crate::topology_errors::TopologyEditError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::TopologyEdit(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=adjacency)]pub fn adjacency(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::Adjacency(_)=>Some(crate::model_errors::AdjacencyError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::Adjacency(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=ringFinding)]pub fn ring_finding(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::RingFinding(_)=>Some(crate::ring_errors::RingFindingError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::RingFinding(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=sourceUInt)]pub fn source_u_int(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::SourceUInt(_)=>Some(crate::property_errors::PropertyUIntReadError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::SourceUInt(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=valence)]pub fn valence(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::Valence(_)=>Some(crate::valence_errors::ValenceError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::Valence(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=stereoOrder)]pub fn stereo_order(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::StereoOrder(_)=>Some(crate::stereo_order_errors::StereoOrderError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::StereoOrder(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=doubleBondStereo)]pub fn double_bond_stereo(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::DoubleBondStereo(_)=>Some(crate::double_stereo_errors::DoubleBondStereoError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::DoubleBondStereo(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=coordinateSelection)]pub fn coordinate_selection(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::CoordinateSelection(_)=>Some(crate::smiles_parse_errors::SmilesParseError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::CoordinateSelection(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=carrierIdentity)]pub fn carrier_identity(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::CarrierIdentity(_)=>Some(crate::query_model_errors::QueryAtomConversionError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::CarrierIdentity(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=moleculeProperty)]pub fn molecule_property(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::MoleculeProperty(_)=>Some(crate::model_errors::MoleculePropertyError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::MoleculeProperty(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=atomProperty)]pub fn atom_property(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::AtomProperty(_)=>Some(crate::model_errors::AtomPropertyError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::AtomProperty(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=bondValue)]pub fn bond_value(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::BondValue(_)=>Some(crate::model_errors::BondValueError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::BondValue(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=templateProperty)]pub fn template_property(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionProductError::TemplateProperty(_)=>Some(crate::reaction_errors::ReactionValidationError::from_owner(self.inner.child(|e|match e{ck::ReactionProductError::TemplateProperty(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
},source=ReactionExecutionCause,cause=reaction_product_error_cause);
fn reaction_run_error_cause(
    owner: &ReadonlyValue<ck::ReactionRunError>,
) -> Option<ReactionExecutionCause> {
    match owner.as_ref() {
        ck::ReactionRunError::Product { .. } => Some(
            ReactionProductError::from_owner(owner.child(|e| match e {
                ck::ReactionRunError::Product { source, .. } => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionRunError::Initialization(_) => {
            crate::reaction_errors::ReactionInitializationError::native_cause(&owner.child(|e| {
                match e {
                    ck::ReactionRunError::Initialization(source) => source,
                    _ => unreachable!("immutable native reaction cause"),
                }
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionRunError::Matching { .. } => Some(
            crate::match_errors::SubstructMatchError::from_owner(owner.child(|e| match e {
                ck::ReactionRunError::Matching { source, .. } => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .unchecked_into(),
        ),
        _ => None,
    }
}
enum_error!(ReactionRunError,ReactionRunErrorKind,ck::ReactionRunError,{Product=>ck::ReactionRunError::Product{..},CoordinateSelectionArity=>ck::ReactionRunError::CoordinateSelectionArity{..},Initialization=>ck::ReactionRunError::Initialization(_),NeedsInitialization=>ck::ReactionRunError::NeedsInitialization,ReactantArity=>ck::ReactionRunError::ReactantArity{..},ReactantTemplateIndex=>ck::ReactionRunError::ReactantTemplateIndex{..},Matching=>ck::ReactionRunError::Matching{..},MatchAtomIndex=>ck::ReactionRunError::MatchAtomIndex{..},CombinationLevel=>ck::ReactionRunError::CombinationLevel{..},CombinationSize=>ck::ReactionRunError::CombinationSize{..},EmptyCombinationLevels=>ck::ReactionRunError::EmptyCombinationLevels},{
#[wasm_bindgen(getter,js_name=set)]pub fn set(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionRunError::Product{set,..}=>Some(*set as usize),_=>None}}
#[wasm_bindgen(getter,js_name=template)]pub fn template(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionRunError::Product{template,..}=>Some(*template as usize),ck::ReactionRunError::Matching{template,..}=>Some(*template as usize),ck::ReactionRunError::MatchAtomIndex{template,..}=>Some(*template as usize),_=>None}}
#[wasm_bindgen(getter,js_name=expected)]pub fn expected(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionRunError::CoordinateSelectionArity{expected,..}=>Some(*expected as usize),ck::ReactionRunError::ReactantArity{expected,..}=>Some(*expected as usize),ck::ReactionRunError::CombinationSize{expected,..}=>Some(*expected as usize),_=>None}}
#[wasm_bindgen(getter,js_name=actual)]pub fn actual(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionRunError::CoordinateSelectionArity{actual,..}=>Some(*actual as usize),ck::ReactionRunError::ReactantArity{actual,..}=>Some(*actual as usize),ck::ReactionRunError::CombinationSize{actual,..}=>Some(*actual as usize),_=>None}}
#[wasm_bindgen(getter,js_name=index)]pub fn index(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionRunError::ReactantTemplateIndex{index,..}=>Some(*index as usize),_=>None}}
#[wasm_bindgen(getter,js_name=count)]pub fn count(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionRunError::ReactantTemplateIndex{count,..}=>Some(*count as usize),ck::ReactionRunError::CombinationLevel{count,..}=>Some(*count as usize),_=>None}}
#[wasm_bindgen(getter,js_name=reactant)]pub fn reactant(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionRunError::Matching{reactant,..}=>Some(*reactant as usize),ck::ReactionRunError::MatchAtomIndex{reactant,..}=>Some(*reactant as usize),_=>None}}
#[wasm_bindgen(getter,js_name=atom)]pub fn atom(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionRunError::MatchAtomIndex{atom,..}=>Some(*atom as usize),_=>None}}
#[wasm_bindgen(getter,js_name=atomCount)]pub fn atom_count(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionRunError::MatchAtomIndex{atom_count,..}=>Some(*atom_count as usize),_=>None}}
#[wasm_bindgen(getter,js_name=level)]pub fn level(&self)->Option<u32>{match self.inner.as_ref(){ck::ReactionRunError::CombinationLevel{level,..}=>Some(*level),_=>None}}
#[wasm_bindgen(getter,js_name=initialization)]pub fn initialization(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionRunError::Initialization(_)=>Some(crate::reaction_errors::ReactionInitializationError::from_owner(self.inner.child(|e|match e{ck::ReactionRunError::Initialization(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
},source=ReactionExecutionCause,cause=reaction_run_error_cause);
fn reaction_apply_error_cause(
    owner: &ReadonlyValue<ck::ReactionApplyError>,
) -> Option<ReactionExecutionCause> {
    match owner.as_ref() {
        ck::ReactionApplyError::SourceText(_) => {
            crate::property_errors::PropertyStringError::native_cause(&owner.child(|e| match e {
                ck::ReactionApplyError::SourceText(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionApplyError::Initialization(_) => {
            crate::reaction_errors::ReactionInitializationError::native_cause(&owner.child(|e| {
                match e {
                    ck::ReactionApplyError::Initialization(source) => source,
                    _ => unreachable!("immutable native reaction cause"),
                }
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionApplyError::Matching(_) => {
            ReactionRunError::native_cause(&owner.child(|e| match e {
                ck::ReactionApplyError::Matching(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionApplyError::Product(_) => {
            ReactionProductError::native_cause(&owner.child(|e| match e {
                ck::ReactionApplyError::Product(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        ck::ReactionApplyError::Edit(_) => {
            crate::topology_errors::TopologyEditError::native_cause(&owner.child(|e| match e {
                ck::ReactionApplyError::Edit(source) => source,
                _ => unreachable!("immutable native reaction cause"),
            }))
            .map(JsCast::unchecked_into)
        }
        _ => None,
    }
}
enum_error!(ReactionApplyError,ReactionApplyErrorKind,ck::ReactionApplyError,{SourceText=>ck::ReactionApplyError::SourceText(_),Initialization=>ck::ReactionApplyError::Initialization(_),ApplicabilityArity=>ck::ReactionApplyError::ApplicabilityArity{..},AddsProductAtom=>ck::ReactionApplyError::AddsProductAtom{..},Matching=>ck::ReactionApplyError::Matching(_),Product=>ck::ReactionApplyError::Product(_),Edit=>ck::ReactionApplyError::Edit(_)},{
#[wasm_bindgen(getter,js_name=reactants)]pub fn reactants(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionApplyError::ApplicabilityArity{reactants,..}=>Some(*reactants as usize),_=>None}}
#[wasm_bindgen(getter,js_name=products)]pub fn products(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionApplyError::ApplicabilityArity{products,..}=>Some(*products as usize),_=>None}}
#[wasm_bindgen(getter,js_name=atom)]pub fn atom(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionApplyError::AddsProductAtom{atom,..}=>Some(atom.index()),_=>None}}
#[wasm_bindgen(getter,js_name=sourceText)]pub fn source_text(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionApplyError::SourceText(_)=>Some(crate::property_errors::PropertyStringError::from_owner(self.inner.child(|e|match e{ck::ReactionApplyError::SourceText(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=initialization)]pub fn initialization(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionApplyError::Initialization(_)=>Some(crate::reaction_errors::ReactionInitializationError::from_owner(self.inner.child(|e|match e{ck::ReactionApplyError::Initialization(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=matching)]pub fn matching(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionApplyError::Matching(_)=>Some(ReactionRunError::from_owner(self.inner.child(|e|match e{ck::ReactionApplyError::Matching(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=product)]pub fn product(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionApplyError::Product(_)=>Some(ReactionProductError::from_owner(self.inner.child(|e|match e{ck::ReactionApplyError::Product(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=edit)]pub fn edit(&self)->Option<ReactionExecutionCause>{match self.inner.as_ref(){ck::ReactionApplyError::Edit(_)=>Some(crate::topology_errors::TopologyEditError::from_owner(self.inner.child(|e|match e{ck::ReactionApplyError::Edit(v)=>v,_=>unreachable!("immutable native payload")})).unchecked_into()),_=>None}}
},source=ReactionExecutionCause,cause=reaction_apply_error_cause);
