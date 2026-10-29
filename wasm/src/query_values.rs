//! Readonly query payloads retain native owners; no query interpretation.
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
pub struct QueryAtom(pub(crate) ReadonlyValue<ck::QueryAtom>);
#[wasm_bindgen]
impl QueryAtom {
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> usize {
        self.0.id().index()
    }
    #[wasm_bindgen(getter,js_name=atomicNumber)]
    pub fn atomic_number(&self) -> u8 {
        self.0.atomic_number()
    }
    #[wasm_bindgen(getter,js_name=formalCharge)]
    pub fn formal_charge(&self) -> i8 {
        self.0.formal_charge()
    }
    #[wasm_bindgen(getter,js_name=explicitHydrogens)]
    pub fn explicit_hydrogens(&self) -> u8 {
        self.0.explicit_hydrogens()
    }
    #[wasm_bindgen(getter)]
    pub fn isotope(&self) -> Option<u16> {
        self.0.isotope()
    }
    #[wasm_bindgen(getter,js_name=atomMap)]
    pub fn atom_map(&self) -> Option<u32> {
        self.0.atom_map()
    }
    #[wasm_bindgen(getter,js_name=isAromatic)]
    pub fn is_aromatic(&self) -> bool {
        self.0.is_aromatic()
    }
    #[wasm_bindgen(getter,js_name=noImplicit)]
    pub fn no_implicit(&self) -> bool {
        self.0.no_implicit()
    }
    #[wasm_bindgen(getter,js_name=radicalElectrons)]
    pub fn radical_electrons(&self) -> u8 {
        self.0.radical_electrons()
    }
}
#[wasm_bindgen]
pub struct QueryGraph(pub(crate) ReadonlyValue<ck::QueryGraph>);
impl QueryGraph {
    pub(crate) fn from_native(value: ck::QueryGraph) -> Self {
        Self(ReadonlyValue::owned(value))
    }
}
crate::fields::value_field!(
    ck::QueryGraph,
    QueryGraph,
    QueryGraphInput,
    "QueryGraph",
    QueryGraph::from_native
);
crate::fields::value_array!(ck::QueryGraph, QueryGraphArray, "ReadonlyArray<QueryGraph>");
#[wasm_bindgen]
impl QueryGraph {
    #[wasm_bindgen(getter,js_name=numAtoms)]
    pub fn num_atoms(&self) -> usize {
        self.0.num_atoms()
    }
    #[wasm_bindgen(getter,js_name=numBonds)]
    pub fn num_bonds(&self) -> usize {
        self.0.num_bonds()
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> Result<Option<String>, JsValue> {
        self.0
            .name()
            .map_err(crate::model_errors::PropertyValueError::from_rust)?
            .map(crate::property_errors::text)
            .transpose()
    }
}
#[wasm_bindgen]
pub struct RecursiveStructureQuery(ReadonlyValue<ck::RecursiveStructureQuery>);
#[wasm_bindgen]
impl RecursiveStructureQuery {
    #[wasm_bindgen(getter,js_name=queryGraph)]
    pub fn query_graph(&self) -> Option<QueryGraph> {
        self.0.query_graph().map(|_| {
            QueryGraph(
                self.0
                    .child(|v| v.query_graph().expect("existing immutable recursive graph")),
            )
        })
    }
    #[wasm_bindgen(getter,js_name=sourceSmarts)]
    pub fn source_smarts(&self) -> Result<Option<String>, JsValue> {
        self.0
            .source_smarts()
            .map(crate::property_errors::text)
            .transpose()
    }
    #[wasm_bindgen(getter,js_name=serialNumber)]
    pub fn serial_number(&self) -> u32 {
        self.0.serial_number()
    }
    #[wasm_bindgen(getter)]
    pub fn initialized(&self) -> bool {
        self.0.is_initialized()
    }
}
#[wasm_bindgen]
pub struct AtomRangeQuery(ck::AtomRangeQuery);
#[wasm_bindgen]
impl AtomRangeQuery {
    #[wasm_bindgen(getter)]
    pub fn bounds(&self) -> AtomRangeBounds {
        AtomRangeBounds(self.0.bounds())
    }
    #[wasm_bindgen(getter,js_name=dataFunction)]
    pub fn data_function(&self) -> AtomRangeDataFunction {
        AtomRangeDataFunction(self.0.data_function())
    }
}
#[wasm_bindgen]
pub struct AtomRangeBounds(ck::AtomRangeBounds);
#[wasm_bindgen]
pub enum AtomRangeBoundsKind {
    LessEqual,
    GreaterEqual,
    Inclusive,
}
#[wasm_bindgen]
impl AtomRangeBounds {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> AtomRangeBoundsKind {
        match &self.0 {
            ck::AtomRangeBounds::LessEqual(_) => AtomRangeBoundsKind::LessEqual,
            ck::AtomRangeBounds::GreaterEqual(_) => AtomRangeBoundsKind::GreaterEqual,
            ck::AtomRangeBounds::Inclusive { .. } => AtomRangeBoundsKind::Inclusive,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn value(&self) -> Option<i32> {
        match self.0 {
            ck::AtomRangeBounds::LessEqual(value) => Some(value),
            ck::AtomRangeBounds::GreaterEqual(value) => Some(value),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=lower)]
    pub fn lower(&self) -> Option<i32> {
        match self.0 {
            ck::AtomRangeBounds::Inclusive { lower, .. } => Some(lower),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=upper)]
    pub fn upper(&self) -> Option<i32> {
        match self.0 {
            ck::AtomRangeBounds::Inclusive { upper, .. } => Some(upper),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=lowerOpen)]
    pub fn lower_open(&self) -> Option<bool> {
        match self.0 {
            ck::AtomRangeBounds::Inclusive { lower_open, .. } => Some(lower_open),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=upperOpen)]
    pub fn upper_open(&self) -> Option<bool> {
        match self.0 {
            ck::AtomRangeBounds::Inclusive { upper_open, .. } => Some(upper_open),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct AtomRangeDataFunction(ck::AtomRangeDataFunction);
#[wasm_bindgen]
pub enum AtomRangeDataFunctionKind {
    AtomRingSize,
    ExplicitDegree,
    NonHydrogenDegree,
    TotalDegree,
    TotalValence,
    NumAtomRings,
    NumHeteroatomNeighbors,
    NumAliphaticHeteroatomNeighbors,
    MinRingSize,
    RingBondCount,
    ImplicitHydrogenCount,
    FormalCharge,
    NegativeFormalCharge,
}
#[wasm_bindgen]
impl AtomRangeDataFunction {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> AtomRangeDataFunctionKind {
        match &self.0 {
            ck::AtomRangeDataFunction::AtomRingSize { .. } => {
                AtomRangeDataFunctionKind::AtomRingSize
            }
            ck::AtomRangeDataFunction::ExplicitDegree => AtomRangeDataFunctionKind::ExplicitDegree,
            ck::AtomRangeDataFunction::NonHydrogenDegree => {
                AtomRangeDataFunctionKind::NonHydrogenDegree
            }
            ck::AtomRangeDataFunction::TotalDegree => AtomRangeDataFunctionKind::TotalDegree,
            ck::AtomRangeDataFunction::TotalValence => AtomRangeDataFunctionKind::TotalValence,
            ck::AtomRangeDataFunction::NumAtomRings => AtomRangeDataFunctionKind::NumAtomRings,
            ck::AtomRangeDataFunction::NumHeteroatomNeighbors => {
                AtomRangeDataFunctionKind::NumHeteroatomNeighbors
            }
            ck::AtomRangeDataFunction::NumAliphaticHeteroatomNeighbors => {
                AtomRangeDataFunctionKind::NumAliphaticHeteroatomNeighbors
            }
            ck::AtomRangeDataFunction::MinRingSize => AtomRangeDataFunctionKind::MinRingSize,
            ck::AtomRangeDataFunction::RingBondCount => AtomRangeDataFunctionKind::RingBondCount,
            ck::AtomRangeDataFunction::ImplicitHydrogenCount => {
                AtomRangeDataFunctionKind::ImplicitHydrogenCount
            }
            ck::AtomRangeDataFunction::FormalCharge => AtomRangeDataFunctionKind::FormalCharge,
            ck::AtomRangeDataFunction::NegativeFormalCharge => {
                AtomRangeDataFunctionKind::NegativeFormalCharge
            }
        }
    }
    #[wasm_bindgen(getter,js_name=lower)]
    pub fn lower(&self) -> Option<i32> {
        match self.0 {
            ck::AtomRangeDataFunction::AtomRingSize { lower, .. } => Some(lower),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=upper)]
    pub fn upper(&self) -> Option<i32> {
        match self.0 {
            ck::AtomRangeDataFunction::AtomRingSize { upper, .. } => Some(upper),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=lowerOpen)]
    pub fn lower_open(&self) -> Option<bool> {
        match self.0 {
            ck::AtomRangeDataFunction::AtomRingSize { lower_open, .. } => Some(lower_open),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=upperOpen)]
    pub fn upper_open(&self) -> Option<bool> {
        match self.0 {
            ck::AtomRangeDataFunction::AtomRingSize { upper_open, .. } => Some(upper_open),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct AtomQueryPredicate(pub(crate) ReadonlyValue<ck::AtomQueryPredicate>);
#[wasm_bindgen]
pub enum AtomQueryPredicateKind {
    AtomicNumber,
    FormalCharge,
    NegativeFormalCharge,
    NumRadicalElectrons,
    Isotope,
    HydrogenCount,
    ImplicitHydrogenCount,
    ImplicitHydrogenCountLessEqual,
    ImplicitValence,
    ExplicitValence,
    ExplicitDegree,
    ExplicitDegreeLessEqual,
    NonHydrogenDegree,
    NonHydrogenDegreeLessEqual,
    NonHydrogenDegreeGreaterEqual,
    HeavyAtomDegree,
    NumHeteroatomNeighbors,
    NumAliphaticHeteroatomNeighbors,
    RingBondCount,
    RingBondCountLessEqual,
    RGroupLabel,
    TotalDegree,
    TotalDegreeLessEqual,
    TotalDegreeGreaterEqual,
    TotalValence,
    TotalValenceLessEqual,
    TotalValenceGreaterEqual,
    NumAtomRings,
    InRingOfSize,
    InRingOfSizeLessEqual,
    InRingOfSizeGreaterEqual,
    SmallestRingSize,
    SmallestRingSizeLessEqual,
    SmallestRingSizeGreaterEqual,
    Mass,
    ChiralPermutationMatch,
    DegreeLessEqual,
    DegreeGreaterEqual,
    Any,
    HasChiralTag,
    MissingChiralTag,
    HasImplicitHydrogen,
    HasHeteroatomNeighbors,
    HasAliphaticHeteroatomNeighbors,
    HasRingBond,
    IsBridgehead,
    IsUnsaturated,
    InRing,
    AtomType,
    PropertyValue,
    AtomicNumberIn,
    AtomicNumberNotIn,
    IsAromatic,
    RecursiveSmarts,
    HasProperty,
    MolFileAlias,
    HybridizationMatch,
    ChiralTagMatch,
    Range,
    UnsupportedFeature,
}
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "number | boolean | string | ReadonlyArray<number> | RecursiveStructureQuery | AtomRangeQuery"
    )]
    pub type AtomQueryPredicateValue;
}
#[wasm_bindgen]
impl AtomQueryPredicate {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> AtomQueryPredicateKind {
        match self.0.as_ref() {
            ck::AtomQueryPredicate::AtomicNumber(_) => AtomQueryPredicateKind::AtomicNumber,
            ck::AtomQueryPredicate::FormalCharge(_) => AtomQueryPredicateKind::FormalCharge,
            ck::AtomQueryPredicate::NegativeFormalCharge(_) => {
                AtomQueryPredicateKind::NegativeFormalCharge
            }
            ck::AtomQueryPredicate::NumRadicalElectrons(_) => {
                AtomQueryPredicateKind::NumRadicalElectrons
            }
            ck::AtomQueryPredicate::Isotope(_) => AtomQueryPredicateKind::Isotope,
            ck::AtomQueryPredicate::HydrogenCount(_) => AtomQueryPredicateKind::HydrogenCount,
            ck::AtomQueryPredicate::ImplicitHydrogenCount(_) => {
                AtomQueryPredicateKind::ImplicitHydrogenCount
            }
            ck::AtomQueryPredicate::ImplicitHydrogenCountLessEqual(_) => {
                AtomQueryPredicateKind::ImplicitHydrogenCountLessEqual
            }
            ck::AtomQueryPredicate::ImplicitValence(_) => AtomQueryPredicateKind::ImplicitValence,
            ck::AtomQueryPredicate::ExplicitValence(_) => AtomQueryPredicateKind::ExplicitValence,
            ck::AtomQueryPredicate::ExplicitDegree(_) => AtomQueryPredicateKind::ExplicitDegree,
            ck::AtomQueryPredicate::ExplicitDegreeLessEqual(_) => {
                AtomQueryPredicateKind::ExplicitDegreeLessEqual
            }
            ck::AtomQueryPredicate::NonHydrogenDegree(_) => {
                AtomQueryPredicateKind::NonHydrogenDegree
            }
            ck::AtomQueryPredicate::NonHydrogenDegreeLessEqual(_) => {
                AtomQueryPredicateKind::NonHydrogenDegreeLessEqual
            }
            ck::AtomQueryPredicate::NonHydrogenDegreeGreaterEqual(_) => {
                AtomQueryPredicateKind::NonHydrogenDegreeGreaterEqual
            }
            ck::AtomQueryPredicate::HeavyAtomDegree(_) => AtomQueryPredicateKind::HeavyAtomDegree,
            ck::AtomQueryPredicate::NumHeteroatomNeighbors(_) => {
                AtomQueryPredicateKind::NumHeteroatomNeighbors
            }
            ck::AtomQueryPredicate::NumAliphaticHeteroatomNeighbors(_) => {
                AtomQueryPredicateKind::NumAliphaticHeteroatomNeighbors
            }
            ck::AtomQueryPredicate::RingBondCount(_) => AtomQueryPredicateKind::RingBondCount,
            ck::AtomQueryPredicate::RingBondCountLessEqual(_) => {
                AtomQueryPredicateKind::RingBondCountLessEqual
            }
            ck::AtomQueryPredicate::RGroupLabel(_) => AtomQueryPredicateKind::RGroupLabel,
            ck::AtomQueryPredicate::TotalDegree(_) => AtomQueryPredicateKind::TotalDegree,
            ck::AtomQueryPredicate::TotalDegreeLessEqual(_) => {
                AtomQueryPredicateKind::TotalDegreeLessEqual
            }
            ck::AtomQueryPredicate::TotalDegreeGreaterEqual(_) => {
                AtomQueryPredicateKind::TotalDegreeGreaterEqual
            }
            ck::AtomQueryPredicate::TotalValence(_) => AtomQueryPredicateKind::TotalValence,
            ck::AtomQueryPredicate::TotalValenceLessEqual(_) => {
                AtomQueryPredicateKind::TotalValenceLessEqual
            }
            ck::AtomQueryPredicate::TotalValenceGreaterEqual(_) => {
                AtomQueryPredicateKind::TotalValenceGreaterEqual
            }
            ck::AtomQueryPredicate::NumAtomRings(_) => AtomQueryPredicateKind::NumAtomRings,
            ck::AtomQueryPredicate::InRingOfSize(_) => AtomQueryPredicateKind::InRingOfSize,
            ck::AtomQueryPredicate::InRingOfSizeLessEqual(_) => {
                AtomQueryPredicateKind::InRingOfSizeLessEqual
            }
            ck::AtomQueryPredicate::InRingOfSizeGreaterEqual(_) => {
                AtomQueryPredicateKind::InRingOfSizeGreaterEqual
            }
            ck::AtomQueryPredicate::SmallestRingSize(_) => AtomQueryPredicateKind::SmallestRingSize,
            ck::AtomQueryPredicate::SmallestRingSizeLessEqual(_) => {
                AtomQueryPredicateKind::SmallestRingSizeLessEqual
            }
            ck::AtomQueryPredicate::SmallestRingSizeGreaterEqual(_) => {
                AtomQueryPredicateKind::SmallestRingSizeGreaterEqual
            }
            ck::AtomQueryPredicate::Mass(_) => AtomQueryPredicateKind::Mass,
            ck::AtomQueryPredicate::ChiralPermutationMatch(_) => {
                AtomQueryPredicateKind::ChiralPermutationMatch
            }
            ck::AtomQueryPredicate::DegreeLessEqual(_) => AtomQueryPredicateKind::DegreeLessEqual,
            ck::AtomQueryPredicate::DegreeGreaterEqual(_) => {
                AtomQueryPredicateKind::DegreeGreaterEqual
            }
            ck::AtomQueryPredicate::AtomicNumberIn(_) => AtomQueryPredicateKind::AtomicNumberIn,
            ck::AtomQueryPredicate::AtomicNumberNotIn(_) => {
                AtomQueryPredicateKind::AtomicNumberNotIn
            }
            ck::AtomQueryPredicate::IsAromatic(_) => AtomQueryPredicateKind::IsAromatic,
            ck::AtomQueryPredicate::RecursiveSmarts(_) => AtomQueryPredicateKind::RecursiveSmarts,
            ck::AtomQueryPredicate::HasProperty(_) => AtomQueryPredicateKind::HasProperty,
            ck::AtomQueryPredicate::MolFileAlias(_) => AtomQueryPredicateKind::MolFileAlias,
            ck::AtomQueryPredicate::HybridizationMatch(_) => {
                AtomQueryPredicateKind::HybridizationMatch
            }
            ck::AtomQueryPredicate::ChiralTagMatch(_) => AtomQueryPredicateKind::ChiralTagMatch,
            ck::AtomQueryPredicate::Range(_) => AtomQueryPredicateKind::Range,
            ck::AtomQueryPredicate::UnsupportedFeature(_) => {
                AtomQueryPredicateKind::UnsupportedFeature
            }
            ck::AtomQueryPredicate::Any => AtomQueryPredicateKind::Any,
            ck::AtomQueryPredicate::HasChiralTag => AtomQueryPredicateKind::HasChiralTag,
            ck::AtomQueryPredicate::MissingChiralTag => AtomQueryPredicateKind::MissingChiralTag,
            ck::AtomQueryPredicate::HasImplicitHydrogen => {
                AtomQueryPredicateKind::HasImplicitHydrogen
            }
            ck::AtomQueryPredicate::HasHeteroatomNeighbors => {
                AtomQueryPredicateKind::HasHeteroatomNeighbors
            }
            ck::AtomQueryPredicate::HasAliphaticHeteroatomNeighbors => {
                AtomQueryPredicateKind::HasAliphaticHeteroatomNeighbors
            }
            ck::AtomQueryPredicate::HasRingBond => AtomQueryPredicateKind::HasRingBond,
            ck::AtomQueryPredicate::IsBridgehead => AtomQueryPredicateKind::IsBridgehead,
            ck::AtomQueryPredicate::IsUnsaturated => AtomQueryPredicateKind::IsUnsaturated,
            ck::AtomQueryPredicate::InRing => AtomQueryPredicateKind::InRing,
            ck::AtomQueryPredicate::AtomType { .. } => AtomQueryPredicateKind::AtomType,
            ck::AtomQueryPredicate::PropertyValue { .. } => AtomQueryPredicateKind::PropertyValue,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> Option<String> {
        match self.0.as_ref() {
            ck::AtomQueryPredicate::PropertyValue { name, .. } => Some(name.clone()),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=atomicNumber)]
    pub fn atomic_number(&self) -> Option<u8> {
        match self.0.as_ref() {
            ck::AtomQueryPredicate::AtomType { atomic_number, .. } => Some(*atomic_number),
            _ => None,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn aromatic(&self) -> Option<bool> {
        match self.0.as_ref() {
            ck::AtomQueryPredicate::AtomType { aromatic, .. } => Some(*aromatic),
            _ => None,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn value(&self) -> Option<AtomQueryPredicateValue> {
        match self.0.as_ref() {
            ck::AtomQueryPredicate::AtomicNumber(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::FormalCharge(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::NegativeFormalCharge(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::NumRadicalElectrons(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::Isotope(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::HydrogenCount(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::ImplicitHydrogenCount(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::ImplicitHydrogenCountLessEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::ImplicitValence(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::ExplicitValence(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::ExplicitDegree(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::ExplicitDegreeLessEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::NonHydrogenDegree(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::NonHydrogenDegreeLessEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::NonHydrogenDegreeGreaterEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::HeavyAtomDegree(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::NumHeteroatomNeighbors(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::NumAliphaticHeteroatomNeighbors(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::RingBondCount(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::RingBondCountLessEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::RGroupLabel(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::TotalDegree(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::TotalDegreeLessEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::TotalDegreeGreaterEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::TotalValence(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::TotalValenceLessEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::TotalValenceGreaterEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::NumAtomRings(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::InRingOfSize(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::InRingOfSizeLessEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::InRingOfSizeGreaterEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::SmallestRingSize(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::SmallestRingSizeLessEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::SmallestRingSizeGreaterEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::Mass(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::ChiralPermutationMatch(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::DegreeLessEqual(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::AtomQueryPredicate::DegreeGreaterEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::AtomQueryPredicate::AtomicNumberIn(v) => Some(
                ({
                    let values: js_sys::Array =
                        v.iter().map(|value| JsValue::from(*value as u32)).collect();
                    js_sys::Object::freeze(&values);
                    JsValue::from(values)
                })
                .unchecked_into(),
            ),
            ck::AtomQueryPredicate::AtomicNumberNotIn(v) => Some(
                ({
                    let values: js_sys::Array =
                        v.iter().map(|value| JsValue::from(*value as u32)).collect();
                    js_sys::Object::freeze(&values);
                    JsValue::from(values)
                })
                .unchecked_into(),
            ),
            ck::AtomQueryPredicate::IsAromatic(v) => Some((JsValue::from(*v)).unchecked_into()),
            ck::AtomQueryPredicate::RecursiveSmarts(_) => Some(
                (JsValue::from(RecursiveStructureQuery(self.0.child(|p| match p {
                    ck::AtomQueryPredicate::RecursiveSmarts(v) => v,
                    _ => unreachable!("immutable recursive predicate"),
                }))))
                .unchecked_into(),
            ),
            ck::AtomQueryPredicate::HasProperty(v) => {
                Some((JsValue::from(v.to_string())).unchecked_into())
            }
            ck::AtomQueryPredicate::MolFileAlias(v) => {
                Some((JsValue::from(v.to_string())).unchecked_into())
            }
            ck::AtomQueryPredicate::HybridizationMatch(v) => {
                Some((JsValue::from(*v as u32)).unchecked_into())
            }
            ck::AtomQueryPredicate::ChiralTagMatch(v) => {
                Some((JsValue::from(*v as u32)).unchecked_into())
            }
            ck::AtomQueryPredicate::Range(v) => {
                Some((JsValue::from(AtomRangeQuery(v.clone()))).unchecked_into())
            }
            ck::AtomQueryPredicate::UnsupportedFeature(v) => {
                Some((JsValue::from(v.to_string())).unchecked_into())
            }
            ck::AtomQueryPredicate::PropertyValue { value, .. } => {
                Some(JsValue::from(value.clone()).unchecked_into())
            }
            ck::AtomQueryPredicate::Any => None,
            ck::AtomQueryPredicate::HasChiralTag => None,
            ck::AtomQueryPredicate::MissingChiralTag => None,
            ck::AtomQueryPredicate::HasImplicitHydrogen => None,
            ck::AtomQueryPredicate::HasHeteroatomNeighbors => None,
            ck::AtomQueryPredicate::HasAliphaticHeteroatomNeighbors => None,
            ck::AtomQueryPredicate::HasRingBond => None,
            ck::AtomQueryPredicate::IsBridgehead => None,
            ck::AtomQueryPredicate::IsUnsaturated => None,
            ck::AtomQueryPredicate::InRing => None,
            ck::AtomQueryPredicate::AtomType { .. } => None,
        }
    }
}
#[wasm_bindgen]
pub struct BondQueryPredicate(pub(crate) ReadonlyValue<ck::BondQueryPredicate>);
#[wasm_bindgen]
pub enum BondQueryPredicateKind {
    NumRingBonds,
    InRingOfSize,
    MinRingSize,
    NumRingBondsGreaterEqual,
    NumRingBondsLessEqual,
    MolFileQueryCode,
    Any,
    HasStereo,
    IsConjugated,
    PropertyValue,
    Order,
    OrderIn,
    IsAromatic,
    IsInRing,
    Direction,
    Stereo,
    HasProperty,
    UnsupportedFeature,
}
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "number | boolean | string | ReadonlyArray<number>")]
    pub type BondQueryPredicateValue;
}
#[wasm_bindgen]
impl BondQueryPredicate {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> BondQueryPredicateKind {
        match self.0.as_ref() {
            ck::BondQueryPredicate::NumRingBonds(_) => BondQueryPredicateKind::NumRingBonds,
            ck::BondQueryPredicate::InRingOfSize(_) => BondQueryPredicateKind::InRingOfSize,
            ck::BondQueryPredicate::MinRingSize(_) => BondQueryPredicateKind::MinRingSize,
            ck::BondQueryPredicate::NumRingBondsGreaterEqual(_) => {
                BondQueryPredicateKind::NumRingBondsGreaterEqual
            }
            ck::BondQueryPredicate::NumRingBondsLessEqual(_) => {
                BondQueryPredicateKind::NumRingBondsLessEqual
            }
            ck::BondQueryPredicate::MolFileQueryCode(_) => BondQueryPredicateKind::MolFileQueryCode,
            ck::BondQueryPredicate::Order(_) => BondQueryPredicateKind::Order,
            ck::BondQueryPredicate::OrderIn(_) => BondQueryPredicateKind::OrderIn,
            ck::BondQueryPredicate::IsAromatic(_) => BondQueryPredicateKind::IsAromatic,
            ck::BondQueryPredicate::IsInRing(_) => BondQueryPredicateKind::IsInRing,
            ck::BondQueryPredicate::Direction(_) => BondQueryPredicateKind::Direction,
            ck::BondQueryPredicate::Stereo(_) => BondQueryPredicateKind::Stereo,
            ck::BondQueryPredicate::HasProperty(_) => BondQueryPredicateKind::HasProperty,
            ck::BondQueryPredicate::UnsupportedFeature(_) => {
                BondQueryPredicateKind::UnsupportedFeature
            }
            ck::BondQueryPredicate::Any => BondQueryPredicateKind::Any,
            ck::BondQueryPredicate::HasStereo => BondQueryPredicateKind::HasStereo,
            ck::BondQueryPredicate::IsConjugated => BondQueryPredicateKind::IsConjugated,
            ck::BondQueryPredicate::PropertyValue { .. } => BondQueryPredicateKind::PropertyValue,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> Option<String> {
        match self.0.as_ref() {
            ck::BondQueryPredicate::PropertyValue { name, .. } => Some(name.clone()),
            _ => None,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn value(&self) -> Option<BondQueryPredicateValue> {
        match self.0.as_ref() {
            ck::BondQueryPredicate::NumRingBonds(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::BondQueryPredicate::InRingOfSize(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::BondQueryPredicate::MinRingSize(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::BondQueryPredicate::NumRingBondsGreaterEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::BondQueryPredicate::NumRingBondsLessEqual(v) => {
                Some(JsValue::from(*v).unchecked_into())
            }
            ck::BondQueryPredicate::MolFileQueryCode(v) => Some(JsValue::from(*v).unchecked_into()),
            ck::BondQueryPredicate::Order(v) => Some((JsValue::from(*v as u32)).unchecked_into()),
            ck::BondQueryPredicate::OrderIn(v) => Some(
                ({
                    let values: js_sys::Array =
                        v.iter().map(|value| JsValue::from(*value as u32)).collect();
                    js_sys::Object::freeze(&values);
                    JsValue::from(values)
                })
                .unchecked_into(),
            ),
            ck::BondQueryPredicate::IsAromatic(v) => Some((JsValue::from(*v)).unchecked_into()),
            ck::BondQueryPredicate::IsInRing(v) => Some((JsValue::from(*v)).unchecked_into()),
            ck::BondQueryPredicate::Direction(v) => {
                Some((JsValue::from(*v as u32)).unchecked_into())
            }
            ck::BondQueryPredicate::Stereo(v) => Some((JsValue::from(*v as u32)).unchecked_into()),
            ck::BondQueryPredicate::HasProperty(v) => {
                Some((JsValue::from(v.to_string())).unchecked_into())
            }
            ck::BondQueryPredicate::UnsupportedFeature(v) => {
                Some((JsValue::from(v.to_string())).unchecked_into())
            }
            ck::BondQueryPredicate::PropertyValue { value, .. } => {
                Some(JsValue::from(value.clone()).unchecked_into())
            }
            ck::BondQueryPredicate::Any => None,
            ck::BondQueryPredicate::HasStereo => None,
            ck::BondQueryPredicate::IsConjugated => None,
        }
    }
}
