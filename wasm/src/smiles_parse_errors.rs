//! Native SMILES errors, including typed writer and parser causes.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
pub struct CxCoordinateSelection(pub(crate) ck::CxCoordinateSelection);
#[wasm_bindgen]
pub enum CxCoordinateSelectionKind {
    Auto,
    TwoD,
    ThreeD,
}
#[wasm_bindgen]
impl CxCoordinateSelection {
    pub fn auto() -> Self {
        Self(ck::CxCoordinateSelection::Auto)
    }
    #[wasm_bindgen(js_name = twoD)]
    pub fn two_d(id: tsify::Ts<crate::input::IndexInput>) -> Result<Self, JsValue> {
        Ok(Self(ck::CxCoordinateSelection::TwoD {
            id: crate::input::read::<usize>(id)?,
        }))
    }
    #[wasm_bindgen(js_name = threeD)]
    pub fn three_d(id: tsify::Ts<crate::input::IndexInput>) -> Result<Self, JsValue> {
        Ok(Self(ck::CxCoordinateSelection::ThreeD {
            id: crate::input::read::<usize>(id)?,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> CxCoordinateSelectionKind {
        match self.0 {
            ck::CxCoordinateSelection::Auto => CxCoordinateSelectionKind::Auto,
            ck::CxCoordinateSelection::TwoD { .. } => CxCoordinateSelectionKind::TwoD,
            ck::CxCoordinateSelection::ThreeD { .. } => CxCoordinateSelectionKind::ThreeD,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> Option<usize> {
        match self.0 {
            ck::CxCoordinateSelection::Auto => None,
            ck::CxCoordinateSelection::TwoD { id } | ck::CxCoordinateSelection::ThreeD { id } => {
                Some(id)
            }
        }
    }
}
crate::fields::value_field!(
    ck::CxCoordinateSelection,
    CxCoordinateSelection,
    CxCoordinateSelectionInput,
    "CxCoordinateSelection"
);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | TopologyEditError | PropertyIntReadError | PropertyULongReadError | RingFindingError | LegacyStereoError | PotentialStereoError | ParserStereoOrderError | StereoOrderError | PropertyValueError | CoordinateValidationError | PropertyUIntReadError | AtomPropertyError | BondValueError | MoleculePropertyError | ParserCarrierError | KekulizeError | PropertyStringError | SanitizeError | HydrogenError | SmilesStereoError | AtropisomerError | RequiredPropertyStringError"
    )]
    #[derive(Clone)]
    pub type SmilesParseCause;
    #[wasm_bindgen(typescript_type = "StereoGroupError")]
    pub type StereoGroupErrorPayload;
    #[wasm_bindgen(typescript_type = "TopologyEditError")]
    pub type TopologyEditErrorPayload;
    #[wasm_bindgen(typescript_type = "PropertyIntReadError")]
    pub type PropertyIntReadErrorPayload;
    #[wasm_bindgen(typescript_type = "PropertyULongReadError")]
    pub type PropertyULongReadErrorPayload;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    pub type RingFindingErrorPayload;
    #[wasm_bindgen(typescript_type = "LegacyStereoError")]
    pub type LegacyStereoErrorPayload;
    #[wasm_bindgen(typescript_type = "PotentialStereoError")]
    pub type PotentialStereoErrorPayload;
    #[wasm_bindgen(typescript_type = "ParserStereoOrderError")]
    pub type ParserStereoOrderErrorPayload;
    #[wasm_bindgen(typescript_type = "StereoOrderError")]
    pub type StereoOrderErrorPayload;
    #[wasm_bindgen(typescript_type = "PropertyValueError")]
    pub type PropertyValueErrorPayload;
    #[wasm_bindgen(typescript_type = "CoordinateValidationError")]
    pub type CoordinateValidationErrorPayload;
    #[wasm_bindgen(typescript_type = "PropertyUIntReadError")]
    pub type PropertyUIntReadErrorPayload;
    #[wasm_bindgen(typescript_type = "AtomPropertyError")]
    pub type AtomPropertyErrorPayload;
    #[wasm_bindgen(typescript_type = "BondValueError")]
    pub type BondValueErrorPayload;
    #[wasm_bindgen(typescript_type = "MoleculePropertyError")]
    pub type MoleculePropertyErrorPayload;
    #[wasm_bindgen(typescript_type = "ParserCarrierError")]
    pub type ParserCarrierErrorPayload;
    #[wasm_bindgen(typescript_type = "KekulizeError")]
    pub type KekulizeErrorPayload;
    #[wasm_bindgen(typescript_type = "PropertyStringError")]
    pub type PropertyStringErrorPayload;
    #[wasm_bindgen(typescript_type = "SanitizeError")]
    pub type SanitizeErrorPayload;
    #[wasm_bindgen(typescript_type = "HydrogenError")]
    pub type HydrogenErrorPayload;
    #[wasm_bindgen(typescript_type = "SmilesStereoError")]
    pub type SmilesStereoErrorPayload;
    #[wasm_bindgen(typescript_type = "AtropisomerError")]
    pub type AtropisomerErrorPayload;
    #[wasm_bindgen(typescript_type = "RequiredPropertyStringError")]
    pub type RequiredPropertyStringErrorPayload;
}
fn cause(error: &ck::SmilesParseError) -> Option<SmilesParseCause> {
    match error {
        ck::SmilesParseError::StereoGroup(v) => {
            Some(crate::model_errors::StereoGroupError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::TopologyEdit(v) => {
            Some(crate::topology_errors::TopologyEditError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::WriterInt(v) => Some(
            crate::property_errors::PropertyIntReadError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::SmilesParseError::WriterULong(v) => Some(
            crate::property_errors::PropertyULongReadError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::SmilesParseError::WriterRings(v) => {
            Some(crate::ring_errors::RingFindingError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::WriterLegacyStereo(v) => Some(
            crate::legacy_stereo_errors::LegacyStereoError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::SmilesParseError::WriterPotentialStereo(v) => Some(
            crate::potential_stereo_errors::PotentialStereoError::from_rust(v.clone())
                .unchecked_into(),
        ),
        ck::SmilesParseError::WriterParserStereoOrder(v) => Some(
            crate::parser_errors::ParserStereoOrderError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::SmilesParseError::WriterStereoOrder(v) => Some(
            crate::stereo_order_errors::StereoOrderError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::SmilesParseError::WriterPropertyKind(v) => {
            Some(crate::model_errors::PropertyValueError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::Coordinates(v) => Some(
            crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::SmilesParseError::WriterNumeric(v) => Some(
            crate::property_errors::PropertyUIntReadError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::SmilesParseError::AtomProperty(v) => {
            Some(crate::model_errors::AtomPropertyError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::BondProperty(v) => {
            Some(crate::model_errors::BondValueError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::MoleculeProperty(v) => {
            Some(crate::model_errors::MoleculePropertyError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::PropertyKind(v) => {
            Some(crate::model_errors::PropertyValueError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::ParserCarrier(v) => {
            Some(crate::parser_errors::ParserCarrierError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::WriterKekulize(v) => {
            Some(crate::kekulize_errors::KekulizeError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::WriterProperty(v) => {
            Some(crate::property_errors::PropertyStringError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::ParserSanitize(v) => {
            Some(crate::sanitize_errors::SanitizeError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::ParserRemoveHydrogens(v) => {
            Some(crate::hydrogen_errors::HydrogenError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::ParserStereo(v) => Some(
            crate::smiles_stereo_errors::SmilesStereoError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::SmilesParseError::ParserAtropisomer(v) => {
            Some(crate::atropisomer_errors::AtropisomerError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::WriterStereoBond(v) => {
            Some(crate::model_errors::BondValueError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesParseError::WriterRequiredProperty { source, .. } => Some(
            crate::property_errors::RequiredPropertyStringError::from_rust(source.clone())
                .unchecked_into(),
        ),
        ck::SmilesParseError::WriterPropertyList { source, .. } => {
            Some(crate::model_errors::AtomPropertyError::from_rust(source.clone()).unchecked_into())
        }
        _ => None,
    }
}
enum_error!(SmilesParseError,SmilesParseErrorKind,ck::SmilesParseError,{
StereoGroup=>ck::SmilesParseError::StereoGroup(_),
TopologyEdit=>ck::SmilesParseError::TopologyEdit(_),
WriterInt=>ck::SmilesParseError::WriterInt(_),
WriterULong=>ck::SmilesParseError::WriterULong(_),
WriterRings=>ck::SmilesParseError::WriterRings(_),
WriterLegacyStereo=>ck::SmilesParseError::WriterLegacyStereo(_),
WriterPotentialStereo=>ck::SmilesParseError::WriterPotentialStereo(_),
WriterParserStereoOrder=>ck::SmilesParseError::WriterParserStereoOrder(_),
WriterStereoOrder=>ck::SmilesParseError::WriterStereoOrder(_),
WriterPropertyKind=>ck::SmilesParseError::WriterPropertyKind(_),
Coordinates=>ck::SmilesParseError::Coordinates(_),
WriterNumeric=>ck::SmilesParseError::WriterNumeric(_),
AtomProperty=>ck::SmilesParseError::AtomProperty(_),
BondProperty=>ck::SmilesParseError::BondProperty(_),
MoleculeProperty=>ck::SmilesParseError::MoleculeProperty(_),
PropertyKind=>ck::SmilesParseError::PropertyKind(_),
ParserCarrier=>ck::SmilesParseError::ParserCarrier(_),
WriterKekulize=>ck::SmilesParseError::WriterKekulize(_),
WriterProperty=>ck::SmilesParseError::WriterProperty(_),
ParserSanitize=>ck::SmilesParseError::ParserSanitize(_),
ParserRemoveHydrogens=>ck::SmilesParseError::ParserRemoveHydrogens(_),
ParserStereo=>ck::SmilesParseError::ParserStereo(_),
ParserAtropisomer=>ck::SmilesParseError::ParserAtropisomer(_),
WriterStereoBond=>ck::SmilesParseError::WriterStereoBond(_),
CxLinkAtomOutOfRange=>ck::SmilesParseError::CxLinkAtomOutOfRange{..},
CxLinkOrderOutOfRange=>ck::SmilesParseError::CxLinkOrderOutOfRange{..},
CxZeroBondOutOfRange=>ck::SmilesParseError::CxZeroBondOutOfRange{..},
CxTypedBondOutOfRange=>ck::SmilesParseError::CxTypedBondOutOfRange{..},
CxAtomPropertyAtomOutOfRange=>ck::SmilesParseError::CxAtomPropertyAtomOutOfRange{..},
CxCoordinateAtomOutOfRange=>ck::SmilesParseError::CxCoordinateAtomOutOfRange{..},
TraversalStateIndex=>ck::SmilesParseError::TraversalStateIndex{..},
WriterRequiredProperty=>ck::SmilesParseError::WriterRequiredProperty{..},
WriterPropertyList=>ck::SmilesParseError::WriterPropertyList{..},
Unsupported=>ck::SmilesParseError::Unsupported{..},
Syntax=>ck::SmilesParseError::Syntax{..},
Atom=>ck::SmilesParseError::Atom{..},
UnclosedRing=>ck::SmilesParseError::UnclosedRing{..},
AmbiguousCoordinateSelection=>ck::SmilesParseError::AmbiguousCoordinateSelection{..},
MissingCoordinateSelection=>ck::SmilesParseError::MissingCoordinateSelection{..},
WriterRootAtomOutOfRange=>ck::SmilesParseError::WriterRootAtomOutOfRange{..},
ReplacementCycle=>ck::SmilesParseError::ReplacementCycle{..},
WriterCanonicalInvariant=>ck::SmilesParseError::WriterCanonicalInvariant(_),
Cx=>ck::SmilesParseError::Cx(_),
UnsupportedCx=>ck::SmilesParseError::UnsupportedCx(_),
UnsupportedWriter=>ck::SmilesParseError::UnsupportedWriter(_),
CanonicalRank=>ck::SmilesParseError::CanonicalRank(_),
WriterValence=>ck::SmilesParseError::WriterValence(_),
WriterStereo=>ck::SmilesParseError::WriterStereo(_),
Model=>ck::SmilesParseError::Model(_),
TraversalTooManyOpenRings=>ck::SmilesParseError::TraversalTooManyOpenRings,
EmptyReplacementKey=>ck::SmilesParseError::EmptyReplacementKey,
},{
#[wasm_bindgen(getter,js_name=stereoGroup)] pub fn stereo_group(&self)->Option<StereoGroupErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::StereoGroup(v)=>Some(crate::model_errors::StereoGroupError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=topologyEdit)] pub fn topology_edit(&self)->Option<TopologyEditErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::TopologyEdit(v)=>Some(crate::topology_errors::TopologyEditError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writerInt)] pub fn writer_int(&self)->Option<PropertyIntReadErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::WriterInt(v)=>Some(crate::property_errors::PropertyIntReadError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writerULong)] pub fn writer_u_long(&self)->Option<PropertyULongReadErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::WriterULong(v)=>Some(crate::property_errors::PropertyULongReadError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writerRings)] pub fn writer_rings(&self)->Option<RingFindingErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::WriterRings(v)=>Some(crate::ring_errors::RingFindingError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writerLegacyStereo)] pub fn writer_legacy_stereo(&self)->Option<LegacyStereoErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::WriterLegacyStereo(v)=>Some(crate::legacy_stereo_errors::LegacyStereoError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writerPotentialStereo)] pub fn writer_potential_stereo(&self)->Option<PotentialStereoErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::WriterPotentialStereo(v)=>Some(crate::potential_stereo_errors::PotentialStereoError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writerParserStereoOrder)] pub fn writer_parser_stereo_order(&self)->Option<ParserStereoOrderErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::WriterParserStereoOrder(v)=>Some(crate::parser_errors::ParserStereoOrderError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writerStereoOrder)] pub fn writer_stereo_order(&self)->Option<StereoOrderErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::WriterStereoOrder(v)=>Some(crate::stereo_order_errors::StereoOrderError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writerPropertyKind)] pub fn writer_property_kind(&self)->Option<PropertyValueErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::WriterPropertyKind(v)=>Some(crate::model_errors::PropertyValueError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=coordinates)] pub fn coordinates(&self)->Option<CoordinateValidationErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::Coordinates(v)=>Some(crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writerNumeric)] pub fn writer_numeric(&self)->Option<PropertyUIntReadErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::WriterNumeric(v)=>Some(crate::property_errors::PropertyUIntReadError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=atomProperty)] pub fn atom_property(&self)->Option<AtomPropertyErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::AtomProperty(v)=>Some(crate::model_errors::AtomPropertyError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=bondProperty)] pub fn bond_property(&self)->Option<BondValueErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::BondProperty(v)=>Some(crate::model_errors::BondValueError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=moleculeProperty)] pub fn molecule_property(&self)->Option<MoleculePropertyErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::MoleculeProperty(v)=>Some(crate::model_errors::MoleculePropertyError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=propertyKind)] pub fn property_kind(&self)->Option<PropertyValueErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::PropertyKind(v)=>Some(crate::model_errors::PropertyValueError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=parserCarrier)] pub fn parser_carrier(&self)->Option<ParserCarrierErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::ParserCarrier(v)=>Some(crate::parser_errors::ParserCarrierError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writerKekulize)] pub fn writer_kekulize(&self)->Option<KekulizeErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::WriterKekulize(v)=>Some(crate::kekulize_errors::KekulizeError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writerProperty)] pub fn writer_property(&self)->Option<PropertyStringErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::WriterProperty(v)=>Some(crate::property_errors::PropertyStringError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=parserSanitize)] pub fn parser_sanitize(&self)->Option<SanitizeErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::ParserSanitize(v)=>Some(crate::sanitize_errors::SanitizeError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=parserRemoveHydrogens)] pub fn parser_remove_hydrogens(&self)->Option<HydrogenErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::ParserRemoveHydrogens(v)=>Some(crate::hydrogen_errors::HydrogenError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=parserStereo)] pub fn parser_stereo(&self)->Option<SmilesStereoErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::ParserStereo(v)=>Some(crate::smiles_stereo_errors::SmilesStereoError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=parserAtropisomer)] pub fn parser_atropisomer(&self)->Option<AtropisomerErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::ParserAtropisomer(v)=>Some(crate::atropisomer_errors::AtropisomerError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writerStereoBond)] pub fn writer_stereo_bond(&self)->Option<BondValueErrorPayload> {match self.inner.as_ref() {ck::SmilesParseError::WriterStereoBond(v)=>Some(crate::model_errors::BondValueError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writerCanonicalInvariant)] pub fn writer_canonical_invariant(&self)->Option<String> {match self.inner.as_ref() {ck::SmilesParseError::WriterCanonicalInvariant(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=cx)] pub fn cx(&self)->Option<String> {match self.inner.as_ref() {ck::SmilesParseError::Cx(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=unsupportedCx)] pub fn unsupported_cx(&self)->Option<String> {match self.inner.as_ref() {ck::SmilesParseError::UnsupportedCx(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=unsupportedWriter)] pub fn unsupported_writer(&self)->Option<String> {match self.inner.as_ref() {ck::SmilesParseError::UnsupportedWriter(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=canonicalRank)] pub fn canonical_rank(&self)->Option<String> {match self.inner.as_ref() {ck::SmilesParseError::CanonicalRank(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=writerValence)] pub fn writer_valence(&self)->Option<String> {match self.inner.as_ref() {ck::SmilesParseError::WriterValence(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=writerStereo)] pub fn writer_stereo(&self)->Option<String> {match self.inner.as_ref() {ck::SmilesParseError::WriterStereo(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=model)] pub fn model(&self)->Option<String> {match self.inner.as_ref() {ck::SmilesParseError::Model(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=atom)] pub fn atom(&self)->Option<usize>{match self.inner.as_ref(){
ck::SmilesParseError::CxLinkAtomOutOfRange{atom,..}=>Some(*atom),
ck::SmilesParseError::CxAtomPropertyAtomOutOfRange{atom,..}=>Some(atom.index()),
ck::SmilesParseError::CxCoordinateAtomOutOfRange{atom,..}=>Some(atom.index()),
ck::SmilesParseError::WriterRequiredProperty{atom,..}=>Some(atom.index()),
ck::SmilesParseError::WriterPropertyList{atom,..}=>Some(atom.index()),
_=>None}}
#[wasm_bindgen(getter,js_name=atomCount)] pub fn atom_count(&self)->Option<usize>{match self.inner.as_ref(){
ck::SmilesParseError::CxLinkAtomOutOfRange{atom_count,..}=>Some(*atom_count),
ck::SmilesParseError::CxAtomPropertyAtomOutOfRange{atom_count,..}=>Some(*atom_count),
ck::SmilesParseError::CxCoordinateAtomOutOfRange{atom_count,..}=>Some(*atom_count),
ck::SmilesParseError::WriterRootAtomOutOfRange{atom_count,..}=>Some(*atom_count),
_=>None}}
#[wasm_bindgen(getter,js_name=index)] pub fn index(&self)->Option<usize>{match self.inner.as_ref(){
ck::SmilesParseError::CxLinkOrderOutOfRange{index,..}=>Some(*index),
ck::SmilesParseError::TraversalStateIndex{index,..}=>Some(*index),
ck::SmilesParseError::UnclosedRing{index,..}=>Some(*index as usize),
_=>None}}
#[wasm_bindgen(getter,js_name=count)] pub fn count(&self)->Option<usize>{match self.inner.as_ref(){
ck::SmilesParseError::CxLinkOrderOutOfRange{count,..}=>Some(*count),
ck::SmilesParseError::TraversalStateIndex{count,..}=>Some(*count),
_=>None}}
#[wasm_bindgen(getter,js_name=bond)] pub fn bond(&self)->Option<usize>{match self.inner.as_ref(){
ck::SmilesParseError::CxZeroBondOutOfRange{bond,..}=>Some(bond.index()),
ck::SmilesParseError::CxTypedBondOutOfRange{bond,..}=>Some(bond.index()),
_=>None}}
#[wasm_bindgen(getter,js_name=bondCount)] pub fn bond_count(&self)->Option<usize>{match self.inner.as_ref(){
ck::SmilesParseError::CxZeroBondOutOfRange{bond_count,..}=>Some(*bond_count),
ck::SmilesParseError::CxTypedBondOutOfRange{bond_count,..}=>Some(*bond_count),
_=>None}}
#[wasm_bindgen(getter,js_name=state)] pub fn state(&self)->Option<String>{match self.inner.as_ref(){
ck::SmilesParseError::TraversalStateIndex{state,..}=>Some(state.to_string()),
_=>None}}
#[wasm_bindgen(getter,js_name=token)] pub fn token(&self)->Option<String>{match self.inner.as_ref(){
ck::SmilesParseError::Unsupported{token,..}=>Some(token.to_string()),
_=>None}}
#[wasm_bindgen(getter,js_name=offset)] pub fn offset(&self)->Option<usize>{match self.inner.as_ref(){
ck::SmilesParseError::Unsupported{offset,..}=>Some(*offset),
ck::SmilesParseError::Syntax{offset,..}=>Some(*offset),
ck::SmilesParseError::Atom{offset,..}=>Some(*offset),
_=>None}}
#[wasm_bindgen(getter,js_name=message)] pub fn message(&self)->Option<String>{match self.inner.as_ref(){
ck::SmilesParseError::Syntax{message,..}=>Some(message.to_string()),
_=>None}}
#[wasm_bindgen(getter,js_name=text)] pub fn text(&self)->Option<String>{match self.inner.as_ref(){
ck::SmilesParseError::Atom{text,..}=>Some(text.to_string()),
_=>None}}
#[wasm_bindgen(getter,js_name=twoDCount)] pub fn two_d_count(&self)->Option<usize>{match self.inner.as_ref(){
ck::SmilesParseError::AmbiguousCoordinateSelection{two_d_count,..}=>Some(*two_d_count),
_=>None}}
#[wasm_bindgen(getter,js_name=threeDCount)] pub fn three_d_count(&self)->Option<usize>{match self.inner.as_ref(){
ck::SmilesParseError::AmbiguousCoordinateSelection{three_d_count,..}=>Some(*three_d_count),
_=>None}}
#[wasm_bindgen(getter,js_name=selection)] pub fn selection(&self)->Option<CxCoordinateSelection>{match self.inner.as_ref(){
ck::SmilesParseError::MissingCoordinateSelection{selection,..}=>Some(CxCoordinateSelection(*selection)),
_=>None}}
#[wasm_bindgen(getter,js_name=atomIndex)] pub fn atom_index(&self)->Option<usize>{match self.inner.as_ref(){
ck::SmilesParseError::WriterRootAtomOutOfRange{atom_index,..}=>Some(*atom_index),
_=>None}}
#[wasm_bindgen(getter,js_name=key)] pub fn key(&self)->Option<String>{match self.inner.as_ref(){
ck::SmilesParseError::ReplacementCycle{key,..}=>Some(key.to_string()),
_=>None}}
#[wasm_bindgen(getter)] pub fn source(&self)->Option<SmilesParseCause> {cause(self.inner.as_ref())}
},source=SmilesParseCause,cause=cause);
