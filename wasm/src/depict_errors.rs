//! Readonly depiction errors and their existing native source graph.
use crate::errors::{enum_error, exception};
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "SmartsParseError | TopologyValidationError | QueryAtomConversionError | RingFindingError | PathError | IoError"
    )]
    #[derive(Clone)]
    pub type Coordinate2DTemplateCause;
    #[wasm_bindgen(
        typescript_type = "PropertyUIntReadError | SubstructMatchError | PathError | MatrixError"
    )]
    #[derive(Clone)]
    pub type Coordinate2DLayoutCause;
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | ValenceError | RingFindingError | LegacyStereoError | Coordinate2DTemplateError | Coordinate2DLayoutError | CoordinateValidationError"
    )]
    #[derive(Clone)]
    pub type Coordinate2DCause;
}
fn template_source(
    e: &ReadonlyValue<ck::Coordinate2DTemplateError>,
) -> Option<Coordinate2DTemplateCause> {
    use ck::Coordinate2DTemplateError as E;
    Some(
        match e.as_ref() {
            E::InvalidDefaultRow { source, .. } | E::ExternalInvalidSmarts { source, .. } => {
                crate::smarts_parse_errors::SmartsParseError::from_rust(source.clone())
            }
            E::InvalidTopology { source, .. } => {
                crate::model_errors::TopologyValidationError::from_rust(source.clone())
            }
            E::NonElementIdentity { source, .. } => {
                crate::query_model_errors::QueryAtomConversionError::from_rust(source.clone())
            }
            E::RingInitialization { source, .. } => {
                crate::ring_errors::RingFindingError::from_rust(source.clone())
            }
            E::ConnectedComponents { source, .. } => {
                crate::path_errors::PathError::from_rust(source.clone())
            }
            E::ExternalOpen { .. } | E::ExternalRead { .. } => {
                crate::external_errors::IoError::from_owner(e.child(|e| match e {
                    E::ExternalOpen { source, .. } | E::ExternalRead { source, .. } => {
                        source.as_ref()
                    }
                    _ => unreachable!("immutable selected external source"),
                }))
            }
            _ => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(Coordinate2DTemplateError,Coordinate2DTemplateErrorKind,ck::Coordinate2DTemplateError,{
    InvalidDefaultRow=>ck::Coordinate2DTemplateError::InvalidDefaultRow{..},InvalidTopology=>ck::Coordinate2DTemplateError::InvalidTopology{..},NonElementIdentity=>ck::Coordinate2DTemplateError::NonElementIdentity{..},RingInitialization=>ck::Coordinate2DTemplateError::RingInitialization{..},ConnectedComponents=>ck::Coordinate2DTemplateError::ConnectedComponents{..},ExternalOpen=>ck::Coordinate2DTemplateError::ExternalOpen{..},ExternalRead=>ck::Coordinate2DTemplateError::ExternalRead{..},ExternalInvalidSmarts=>ck::Coordinate2DTemplateError::ExternalInvalidSmarts{..},MissingCoordinates=>ck::Coordinate2DTemplateError::MissingCoordinates{..},ThreeDimensionalCoordinates=>ck::Coordinate2DTemplateError::ThreeDimensionalCoordinates{..},MultipleFragments=>ck::Coordinate2DTemplateError::MultipleFragments{..},NotRingSystem=>ck::Coordinate2DTemplateError::NotRingSystem{..},
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<Coordinate2DTemplateCause>{template_source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn index(&self)->Option<usize>{use ck::Coordinate2DTemplateError as E;match self.inner.as_ref(){E::InvalidDefaultRow{index,..}|E::InvalidTopology{index,..}|E::NonElementIdentity{index,..}|E::RingInitialization{index,..}|E::ConnectedComponents{index,..}=>Some(*index),_=>None}}
    #[wasm_bindgen(getter)] pub fn path(&self)->Option<String>{use ck::Coordinate2DTemplateError as E;match self.inner.as_ref(){E::ExternalOpen{path,..}|E::ExternalRead{path,..}|E::ExternalInvalidSmarts{path,..}=>Some(path.clone()),_=>None}}
    #[wasm_bindgen(getter)] pub fn line(&self)->Option<usize>{match self.inner.as_ref(){ck::Coordinate2DTemplateError::ExternalInvalidSmarts{line,..}=>Some(*line),_=>None}}
    #[wasm_bindgen(getter)] pub fn row(&self)->Option<String>{use ck::Coordinate2DTemplateError as E;match self.inner.as_ref(){E::MissingCoordinates{row}|E::ThreeDimensionalCoordinates{row}|E::MultipleFragments{row}|E::NotRingSystem{row}=>Some(row.clone()),_=>None}}
},source=Coordinate2DTemplateCause,cause=template_source);
fn layout_source(e: &ck::Coordinate2DLayoutError) -> Option<Coordinate2DLayoutCause> {
    use ck::Coordinate2DLayoutError as E;
    Some(
        match e {
            E::InvalidRankProperty { source, .. } => {
                crate::property_errors::PropertyUIntReadError::from_rust(source.clone())
            }
            E::TemplateMatch(v) => crate::match_errors::SubstructMatchError::from_rust(v.clone()),
            E::GraphPath(v) => crate::path_errors::PathError::from_rust(v.clone()),
            E::GraphDistance(v) => crate::matrix_errors::MatrixError::from_rust(v.clone()),
            _ => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(Coordinate2DLayoutError,Coordinate2DLayoutErrorKind,ck::Coordinate2DLayoutError,{
    AtomIndexOutOfRange=>ck::Coordinate2DLayoutError::AtomIndexOutOfRange{..},AtomAlreadyEmbedded=>ck::Coordinate2DLayoutError::AtomAlreadyEmbedded{..},AtomNotEmbedded=>ck::Coordinate2DLayoutError::AtomNotEmbedded{..},NotEnoughEmbeddedNeighbors=>ck::Coordinate2DLayoutError::NotEnoughEmbeddedNeighbors{..},CoincidentPoints=>ck::Coordinate2DLayoutError::CoincidentPoints,InvalidAngle=>ck::Coordinate2DLayoutError::InvalidAngle,NoCommonAtoms=>ck::Coordinate2DLayoutError::NoCommonAtoms,MismatchedTopology=>ck::Coordinate2DLayoutError::MismatchedTopology,EmptyAttachment=>ck::Coordinate2DLayoutError::EmptyAttachment{..},NonTetrahedralNoLigand=>ck::Coordinate2DLayoutError::NonTetrahedralNoLigand{..},NonTetrahedralLigandOverflow=>ck::Coordinate2DLayoutError::NonTetrahedralLigandOverflow{..},CisTransBondInvalid=>ck::Coordinate2DLayoutError::CisTransBondInvalid{..},CollisionBondInvalid=>ck::Coordinate2DLayoutError::CollisionBondInvalid{..},UndefinedSamplingDistance=>ck::Coordinate2DLayoutError::UndefinedSamplingDistance{..},AtomCountTooLarge=>ck::Coordinate2DLayoutError::AtomCountTooLarge{..},InvalidRankProperty=>ck::Coordinate2DLayoutError::InvalidRankProperty{..},GeometryNotEnoughNeighbors=>ck::Coordinate2DLayoutError::GeometryNotEnoughNeighbors{..},TemplateMatch=>ck::Coordinate2DLayoutError::TemplateMatch(_),GraphPath=>ck::Coordinate2DLayoutError::GraphPath(_),GraphDistance=>ck::Coordinate2DLayoutError::GraphDistance(_),
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<Coordinate2DLayoutCause>{layout_source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn atom(&self)->Option<usize>{use ck::Coordinate2DLayoutError as E;match self.inner.as_ref(){E::AtomIndexOutOfRange{atom,..}|E::AtomAlreadyEmbedded{atom}|E::AtomNotEmbedded{atom}|E::NotEnoughEmbeddedNeighbors{atom,..}|E::EmptyAttachment{atom}|E::InvalidRankProperty{atom,..}|E::GeometryNotEnoughNeighbors{atom,..}=>Some(*atom),_=>None}}
    #[wasm_bindgen(getter,js_name=atomCount)] pub fn atom_count(&self)->Option<usize>{use ck::Coordinate2DLayoutError as E;match self.inner.as_ref(){E::AtomIndexOutOfRange{atom_count,..}|E::AtomCountTooLarge{atom_count}=>Some(*atom_count),_=>None}}
    #[wasm_bindgen(getter)] pub fn count(&self)->Option<usize>{use ck::Coordinate2DLayoutError as E;match self.inner.as_ref(){E::NotEnoughEmbeddedNeighbors{count,..}|E::GeometryNotEnoughNeighbors{count,..}=>Some(*count),_=>None}}
    #[wasm_bindgen(getter)] pub fn centre(&self)->Option<usize>{use ck::Coordinate2DLayoutError as E;match self.inner.as_ref(){E::NonTetrahedralNoLigand{centre}|E::NonTetrahedralLigandOverflow{centre}=>Some(*centre),_=>None}}
    #[wasm_bindgen(getter)] pub fn bond(&self)->Option<usize>{use ck::Coordinate2DLayoutError as E;match self.inner.as_ref(){E::CisTransBondInvalid{bond}|E::CollisionBondInvalid{bond}=>Some(*bond),_=>None}}
    #[wasm_bindgen(getter)] pub fn first(&self)->Option<usize>{match self.inner.as_ref(){ck::Coordinate2DLayoutError::UndefinedSamplingDistance{first,..}=>Some(*first),_=>None}}
    #[wasm_bindgen(getter)] pub fn second(&self)->Option<usize>{match self.inner.as_ref(){ck::Coordinate2DLayoutError::UndefinedSamplingDistance{second,..}=>Some(*second),_=>None}}
    #[wasm_bindgen(getter)] pub fn key(&self)->Option<String>{match self.inner.as_ref(){ck::Coordinate2DLayoutError::InvalidRankProperty{key,..}=>Some((*key).into()),_=>None}}
},source=Coordinate2DLayoutCause,cause=layout_source);
fn depict_source(e: &ck::Coordinate2DError) -> Option<Coordinate2DCause> {
    use ck::Coordinate2DError as E;
    Some(
        match e {
            E::InvalidTopology(v) => {
                crate::model_errors::TopologyValidationError::from_rust(v.clone())
            }
            E::PropertyCache(v) => crate::valence_errors::ValenceError::from_rust(v.clone()),
            E::RingFinding(v) => crate::ring_errors::RingFindingError::from_rust(v.clone()),
            E::StereoAssignment(v) => {
                crate::legacy_stereo_errors::LegacyStereoError::from_rust(v.clone())
            }
            E::TemplateLoading(v) => Coordinate2DTemplateError::from_rust(v.clone()),
            E::Fragment(v) => Coordinate2DLayoutError::from_rust(v.clone()),
            E::CoordinateValidation(v) => {
                crate::model_errors::CoordinateValidationError::from_rust(v.clone())
            }
            _ => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(Coordinate2DError,Coordinate2DErrorKind,ck::Coordinate2DError,{
    InvalidTopology=>ck::Coordinate2DError::InvalidTopology(_),PropertyCache=>ck::Coordinate2DError::PropertyCache(_),RingFinding=>ck::Coordinate2DError::RingFinding(_),StereoAssignment=>ck::Coordinate2DError::StereoAssignment(_),TemplateLoading=>ck::Coordinate2DError::TemplateLoading(_),Fragment=>ck::Coordinate2DError::Fragment(_),CoordinateValidation=>ck::Coordinate2DError::CoordinateValidation(_),ConformerIdOverflow=>ck::Coordinate2DError::ConformerIdOverflow{..},CoordGenUnavailable=>ck::Coordinate2DError::CoordGenUnavailable,
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<Coordinate2DCause>{depict_source(&self.inner)}
    #[wasm_bindgen(getter,js_name=maxId)] pub fn max_id(&self)->Option<usize>{match self.inner.as_ref(){ck::Coordinate2DError::ConformerIdOverflow{max_id}=>Some(*max_id),_=>None}}
},source=Coordinate2DCause,cause=depict_source);
pub(crate) fn cause(e: &ck::Coordinate2DError) -> Option<Coordinate2DCause> {
    depict_source(e)
}
