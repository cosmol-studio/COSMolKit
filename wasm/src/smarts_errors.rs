//! SMARTS failures retain their native root, fields, and source graph.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "ValenceError | LegacyStereoError | SmilesParseError")]
    #[derive(Clone)]
    pub type SmartsTraversalCause;
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | SmartsTraversalError | CoordinateValidationError | RingFindingError | WedgeError | SmilesParseError | PropertyUIntReadError | MoleculePropertyError | AtropisomerError | AtomPropertyError | PropertyIntReadError | PropertyValueError | RequiredPropertyStringError | QueryGraphError | PropertyStringError | ValenceError"
    )]
    #[derive(Clone)]
    pub type SmartsWriteCause;
    #[wasm_bindgen(typescript_type = "number | PropertyValueKind")]
    pub type SmartsActual;
    #[wasm_bindgen(typescript_type = "number | string")]
    pub type SmartsExpected;
    #[wasm_bindgen(typescript_type = "string | PropertyValueKind")]
    pub type SmartsKind;
    #[wasm_bindgen(typescript_type = "AtomQueryPredicate | BondQueryPredicate")]
    pub type SmartsPredicate;
}
fn traversal_source(e: &ck::SmartsTraversalError) -> Option<SmartsTraversalCause> {
    Some(
        match e {
            ck::SmartsTraversalError::Valence(v) => {
                crate::valence_errors::ValenceError::from_rust(v.clone())
            }
            ck::SmartsTraversalError::Stereo(v) => {
                crate::legacy_stereo_errors::LegacyStereoError::from_rust(v.as_ref().clone())
            }
            ck::SmartsTraversalError::Traversal(v) => {
                crate::smiles_parse_errors::SmilesParseError::from_rust(v.clone())
            }
        }
        .unchecked_into(),
    )
}
enum_error!(SmartsTraversalError,SmartsTraversalErrorKind,ck::SmartsTraversalError,{
    Valence=>ck::SmartsTraversalError::Valence(_),
    Stereo=>ck::SmartsTraversalError::Stereo(_),
    Traversal=>ck::SmartsTraversalError::Traversal(_),
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<SmartsTraversalCause>{traversal_source(&self.inner)}
},source=SmartsTraversalCause,cause=traversal_source);

pub(crate) fn write_source(e: &ck::SmartsWriteError) -> Option<SmartsWriteCause> {
    use ck::SmartsWriteError as E;
    let value = match e {
        E::StereoGroup(v) => crate::model_errors::StereoGroupError::from_rust(v.clone()),
        E::Traversal(v) => SmartsTraversalError::from_rust(v.clone()),
        E::CxCoordinates(v)
        | E::CxCoordinateSource(v)
        | E::CxCoordinateStorage { source: v, .. } => {
            crate::model_errors::CoordinateValidationError::from_rust(v.clone())
        }
        E::CxRingInfo(v) => crate::ring_errors::RingFindingError::from_rust(v.clone()),
        E::CxWedge(v) => crate::wedge_errors::WedgeError::from_rust(v.clone()),
        E::CxAtomPropertyOutput(v)
        | E::CxCoordinateOutput(v)
        | E::CanonicalTraversal(v)
        | E::CxCoordinateSelection { source: v, .. } => {
            crate::smiles_parse_errors::SmilesParseError::from_rust(v.clone())
        }
        E::CxMoleculePropertyUInt { source: v, .. }
        | E::CxAtomPropertyUInt { source: v, .. }
        | E::CxBondPropertyUInt { source: v, .. }
        | E::CxSgroupPropertyUInt { source: v, .. } => {
            crate::property_errors::PropertyUIntReadError::from_rust(v.clone())
        }
        E::CxSgroupPropertyWrite { source: v, .. } | E::MoleculePropertyWrite(v) => {
            crate::model_errors::MoleculePropertyError::from_rust(v.clone())
        }
        E::CxStereoGroup(v) => crate::atropisomer_errors::AtropisomerError::from_rust(v.clone()),
        E::AtomPropertyWrite { source: v, .. }
        | E::CxPropertyList { source: v, .. }
        | E::CxAtomPropertyWrite { source: v, .. } => {
            crate::model_errors::AtomPropertyError::from_rust(v.clone())
        }
        E::AtomMapInt { source: v, .. } | E::CxAtomPropertyInt { source: v, .. } => {
            crate::property_errors::PropertyIntReadError::from_rust(v.clone())
        }
        E::PropertyValue(v) => crate::model_errors::PropertyValueError::from_rust(v.clone()),
        E::CxRequiredProperty { source: v, .. } => {
            crate::property_errors::RequiredPropertyStringError::from_rust(v.clone())
        }
        E::CxComposition(v) => crate::query_model_errors::QueryGraphError::from_rust(v.clone()),
        E::Property(v) => crate::property_errors::PropertyStringError::from_rust(*v),
        E::Valence(v) => crate::valence_errors::ValenceError::from_rust(v.clone()),
        _ => return None,
    };
    Some(value.unchecked_into())
}
enum_error!(SmartsWriteError,SmartsWriteErrorKind,ck::SmartsWriteError,{
    StereoGroup=>ck::SmartsWriteError::StereoGroup(_),Traversal=>ck::SmartsWriteError::Traversal(_),
    CxCoordinates=>ck::SmartsWriteError::CxCoordinates(_),CxRingInfo=>ck::SmartsWriteError::CxRingInfo(_),
    CxRingAtomOrderIndex=>ck::SmartsWriteError::CxRingAtomOrderIndex{..},CxRingStereoReferenceMissing=>ck::SmartsWriteError::CxRingStereoReferenceMissing{..},
    CxWedge=>ck::SmartsWriteError::CxWedge(_),CxBondConfigAtropMissingCarriers=>ck::SmartsWriteError::CxBondConfigAtropMissingCarriers{..},
    CxAtomPropertyOutput=>ck::SmartsWriteError::CxAtomPropertyOutput(_),CxMissingConformer=>ck::SmartsWriteError::CxMissingConformer,
    CxCoordinateSource=>ck::SmartsWriteError::CxCoordinateSource(_),CxCoordinateOutput=>ck::SmartsWriteError::CxCoordinateOutput(_),
    CxSourceBondOutOfRange=>ck::SmartsWriteError::CxSourceBondOutOfRange{..},CxMoleculePropertyUInt=>ck::SmartsWriteError::CxMoleculePropertyUInt{..},
    CxSgroupVectorCast=>ck::SmartsWriteError::CxSgroupVectorCast{..},CxSgroupPropertyWrite=>ck::SmartsWriteError::CxSgroupPropertyWrite{..},
    CxStereoGroup=>ck::SmartsWriteError::CxStereoGroup(_),CxSourceAtomOutOfRange=>ck::SmartsWriteError::CxSourceAtomOutOfRange{..},
    MoleculePropertyWrite=>ck::SmartsWriteError::MoleculePropertyWrite(_),SourceAtomCount=>ck::SmartsWriteError::SourceAtomCount{..},
    AtomPropertyWrite=>ck::SmartsWriteError::AtomPropertyWrite{..},CanonicalTraversal=>ck::SmartsWriteError::CanonicalTraversal(_),
    UnwritableBondQuery=>ck::SmartsWriteError::UnwritableBondQuery{..},SourceAtomToLeftIndex=>ck::SmartsWriteError::SourceAtomToLeftIndex{..},
    SourceBondBeginIndex=>ck::SmartsWriteError::SourceBondBeginIndex{..},AtomMapInt=>ck::SmartsWriteError::AtomMapInt{..},
    AtomTypeAtomicNumber=>ck::SmartsWriteError::AtomTypeAtomicNumber{..},ChargeMagnitudeOverflow=>ck::SmartsWriteError::ChargeMagnitudeOverflow{..},
    PropertyValue=>ck::SmartsWriteError::PropertyValue(_),CxRequiredProperty=>ck::SmartsWriteError::CxRequiredProperty{..},
    CxPropertyList=>ck::SmartsWriteError::CxPropertyList{..},CxCoordinateSelectionArity=>ck::SmartsWriteError::CxCoordinateSelectionArity{..},
    CxMissingOutputOrder=>ck::SmartsWriteError::CxMissingOutputOrder{..},CxOutputOrderPropertyType=>ck::SmartsWriteError::CxOutputOrderPropertyType{..},
    CxAtomPropertyKind=>ck::SmartsWriteError::CxAtomPropertyKind{..},CxAtomPropertyUInt=>ck::SmartsWriteError::CxAtomPropertyUInt{..},
    CxAtomPropertyWrite=>ck::SmartsWriteError::CxAtomPropertyWrite{..},CxCoordinateStorage=>ck::SmartsWriteError::CxCoordinateStorage{..},
    CxCoordinateSelection=>ck::SmartsWriteError::CxCoordinateSelection{..},CxComposition=>ck::SmartsWriteError::CxComposition(_),
    CxOutputOrder=>ck::SmartsWriteError::CxOutputOrder{..},CxRowCount=>ck::SmartsWriteError::CxRowCount{..},
    CxBondPropertyUInt=>ck::SmartsWriteError::CxBondPropertyUInt{..},InvalidPropertyKind=>ck::SmartsWriteError::InvalidPropertyKind{..},
    CxAtomPropertyInt=>ck::SmartsWriteError::CxAtomPropertyInt{..},CxSgroupPropertyUInt=>ck::SmartsWriteError::CxSgroupPropertyUInt{..},
    Property=>ck::SmartsWriteError::Property(_),Valence=>ck::SmartsWriteError::Valence(_),InvalidGraph=>ck::SmartsWriteError::InvalidGraph(_),
    QueryGraphTraversalUnsupported=>ck::SmartsWriteError::QueryGraphTraversalUnsupported{..},OrAboveAndBelowAnd=>ck::SmartsWriteError::OrAboveAndBelowAnd,
    UnknownCombination=>ck::SmartsWriteError::UnknownCombination{..},MissingRecursiveQueryMolecule=>ck::SmartsWriteError::MissingRecursiveQueryMolecule,
    SourceBondDirection=>ck::SmartsWriteError::SourceBondDirection{..},UnsupportedBondQuery=>ck::SmartsWriteError::UnsupportedBondQuery{..},
    UnsupportedAtomQuery=>ck::SmartsWriteError::UnsupportedAtomQuery{..},CompositeChildCount=>ck::SmartsWriteError::CompositeChildCount{..},XorComposite=>ck::SmartsWriteError::XorComposite,
    QueryGraphCxExtensionsUnsupported=>ck::SmartsWriteError::QueryGraphCxExtensionsUnsupported{..},RootedAtomOutOfRange=>ck::SmartsWriteError::RootedAtomOutOfRange{..},
    EmptyAtomSelection=>ck::SmartsWriteError::EmptyAtomSelection,EmptyBondSelection=>ck::SmartsWriteError::EmptyBondSelection,
    FragmentAtomOutOfRange=>ck::SmartsWriteError::FragmentAtomOutOfRange{..},FragmentBondOutOfRange=>ck::SmartsWriteError::FragmentBondOutOfRange{..},BondAtomNotEndpoint=>ck::SmartsWriteError::BondAtomNotEndpoint{..},
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<SmartsWriteCause>{write_source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn index(&self)->Option<usize>{use ck::SmartsWriteError as E;match self.inner.as_ref(){E::CxRingAtomOrderIndex{index,..}|E::SourceAtomToLeftIndex{index}|E::SourceBondBeginIndex{index}=>Some(*index),_=>None}}
    #[wasm_bindgen(getter)] pub fn count(&self)->Option<usize>{use ck::SmartsWriteError as E;match self.inner.as_ref(){E::CxRingAtomOrderIndex{count,..}|E::SourceAtomCount{count}|E::CxOutputOrder{count,..}|E::CxRowCount{count,..}=>Some(*count),_=>None}}
    #[wasm_bindgen(getter)] pub fn side(&self)->Option<usize>{match self.inner.as_ref(){ck::SmartsWriteError::CxRingStereoReferenceMissing{side,..}=>Some(*side),_=>None}}
    #[wasm_bindgen(getter)] pub fn bond(&self)->Option<usize>{use ck::SmartsWriteError as E;match self.inner.as_ref(){E::CxRingStereoReferenceMissing{bond,..}|E::CxBondConfigAtropMissingCarriers{bond}|E::CxSourceBondOutOfRange{bond,..}|E::CxBondPropertyUInt{bond,..}|E::InvalidPropertyKind{bond,..}=>Some(bond.index()),E::FragmentBondOutOfRange{bond}|E::BondAtomNotEndpoint{bond,..}=>Some(*bond),_=>None}}
    #[wasm_bindgen(getter)] pub fn atom(&self)->Option<usize>{use ck::SmartsWriteError as E;match self.inner.as_ref(){E::CxSourceAtomOutOfRange{atom,..}|E::AtomPropertyWrite{atom,..}|E::AtomMapInt{atom,..}|E::CxRequiredProperty{atom,..}|E::CxPropertyList{atom,..}|E::CxAtomPropertyKind{atom,..}|E::CxAtomPropertyUInt{atom,..}|E::CxAtomPropertyWrite{atom,..}|E::CxAtomPropertyInt{atom,..}=>Some(atom.index()),E::RootedAtomOutOfRange{atom}|E::FragmentAtomOutOfRange{atom}|E::BondAtomNotEndpoint{atom,..}=>Some(*atom),_=>None}}
    #[wasm_bindgen(getter,js_name=bondCount)] pub fn bond_count(&self)->Option<usize>{match self.inner.as_ref(){ck::SmartsWriteError::CxSourceBondOutOfRange{bond_count,..}=>Some(*bond_count),_=>None}}
    #[wasm_bindgen(getter,js_name=atomCount)] pub fn atom_count(&self)->Option<usize>{match self.inner.as_ref(){ck::SmartsWriteError::CxSourceAtomOutOfRange{atom_count,..}=>Some(*atom_count),_=>None}}
    #[wasm_bindgen(getter)] pub fn group(&self)->Option<usize>{use ck::SmartsWriteError as E;match self.inner.as_ref(){E::CxSgroupVectorCast{group,..}|E::CxSgroupPropertyWrite{group,..}|E::CxSgroupPropertyUInt{group,..}=>Some(group.index()),_=>None}}
    #[wasm_bindgen(getter)] pub fn property(&self)->Option<String>{use ck::SmartsWriteError as E;match self.inner.as_ref(){E::CxMoleculePropertyUInt{property,..}|E::CxSgroupVectorCast{property,..}|E::CxSgroupPropertyWrite{property,..}|E::CxOutputOrderPropertyType{property,..}|E::CxAtomPropertyKind{property,..}|E::CxAtomPropertyUInt{property,..}|E::CxAtomPropertyWrite{property,..}|E::CxBondPropertyUInt{property,..}|E::InvalidPropertyKind{property,..}|E::CxAtomPropertyInt{property,..}|E::CxSgroupPropertyUInt{property,..}=>Some((*property).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn actual(&self)->Option<SmartsActual>{use ck::SmartsWriteError as E;match self.inner.as_ref(){E::CxSgroupVectorCast{actual,..}=>Some(JsValue::from(crate::model_errors::PropertyValueKind::from(*actual) as u32).unchecked_into()),E::CxCoordinateSelectionArity{actual,..}=>Some(JsValue::from(*actual as u32).unchecked_into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn expected(&self)->Option<SmartsExpected>{use ck::SmartsWriteError as E;match self.inner.as_ref(){E::CxSgroupVectorCast{expected,..}=>Some(JsValue::from(*expected).unchecked_into()),E::CxCoordinateSelectionArity{expected,..}=>Some(JsValue::from(*expected as u32).unchecked_into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn kind(&self)->Option<SmartsKind>{use ck::SmartsWriteError as E;match self.inner.as_ref(){E::CxAtomPropertyKind{kind,..}|E::InvalidPropertyKind{kind,..}=>Some(JsValue::from(crate::model_errors::PropertyValueKind::from(*kind) as u32).unchecked_into()),E::CxOutputOrder{kind,..}|E::CxRowCount{kind,..}|E::CompositeChildCount{kind}=>Some(JsValue::from(*kind).unchecked_into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn template(&self)->Option<usize>{use ck::SmartsWriteError as E;match self.inner.as_ref(){E::CxMissingOutputOrder{template}|E::CxOutputOrderPropertyType{template,..}|E::CxAtomPropertyKind{template,..}|E::CxAtomPropertyUInt{template,..}|E::CxAtomPropertyWrite{template,..}|E::CxCoordinateStorage{template,..}|E::CxCoordinateSelection{template,..}|E::CxOutputOrder{template,..}=>Some(*template),_=>None}}
    #[wasm_bindgen(getter)] pub fn row(&self)->Option<usize>{match self.inner.as_ref(){ck::SmartsWriteError::CxOutputOrder{row,..}=>Some(*row),_=>None}}
    #[wasm_bindgen(getter)] pub fn description(&self)->Option<String>{use ck::SmartsWriteError as E;match self.inner.as_ref(){E::UnwritableBondQuery{description}=>Some((*description).into()),E::UnknownCombination{description}=>Some(description.clone()),_=>None}}
    #[wasm_bindgen(getter)] pub fn detail(&self)->Option<String>{use ck::SmartsWriteError as E;match self.inner.as_ref(){E::QueryGraphTraversalUnsupported{detail}|E::QueryGraphCxExtensionsUnsupported{detail}=>Some((*detail).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=atomicNumber)] pub fn atomic_number(&self)->Option<i32>{match self.inner.as_ref(){ck::SmartsWriteError::AtomTypeAtomicNumber{atomic_number}=>Some(*atomic_number),_=>None}}
    #[wasm_bindgen(getter)] pub fn value(&self)->Option<i32>{match self.inner.as_ref(){ck::SmartsWriteError::ChargeMagnitudeOverflow{value}=>Some(*value),_=>None}}
    #[wasm_bindgen(getter,js_name=invalidGraph)] pub fn invalid_graph(&self)->Option<String>{match self.inner.as_ref(){ck::SmartsWriteError::InvalidGraph(v)=>Some(v.clone()),_=>None}}
    #[wasm_bindgen(getter)] pub fn direction(&self)->Option<crate::enums::BondDirection>{match self.inner.as_ref(){ck::SmartsWriteError::SourceBondDirection{direction}=>Some((*direction).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn predicate(&self)->Option<SmartsPredicate>{match self.inner.as_ref(){ck::SmartsWriteError::UnsupportedAtomQuery{..}=>Some(JsValue::from(crate::query_values::AtomQueryPredicate(self.inner.child(|e|match e{ck::SmartsWriteError::UnsupportedAtomQuery{predicate}=>predicate,_=>unreachable!("immutable selected error variant")}))).unchecked_into()),ck::SmartsWriteError::UnsupportedBondQuery{..}=>Some(JsValue::from(crate::query_values::BondQueryPredicate(self.inner.child(|e|match e{ck::SmartsWriteError::UnsupportedBondQuery{predicate}=>predicate,_=>unreachable!("immutable selected error variant")}))).unchecked_into()),_=>None}}
},source=SmartsWriteCause,cause=write_source);
