//! Detached 2D layout and depiction boundaries.

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;
use std::sync::Arc;

use cosmolkit_core::{
    LegacyStereoError, RingFindingError, RingSearchParams, ValenceError, ValenceModel,
    ValenceParams, assign_legacy_stereochemistry_for_depiction, assign_valence, symmetrized_sssr,
};
use cosmolkit_model::{
    Conformer2D, CoordinateBlock, CoordinateValidationError, MoleculeProperties, TopologyBlock,
    TopologyValidationError,
};

use crate::embedded_frag::{
    EmbeddedFrag, FragmentError, embed_cis_trans_systems, embed_fused_systems,
    orient_and_shift_fragments, seed_coordinate_constraints,
    translate_single_coordinate_constraint,
};
use crate::geometry::{GeometryError, PointMap, atom_depict_rank};
use crate::nontetrahedral::embed_nontetrahedral_stereo;
use crate::templates::{CoordinateTemplates, TemplateError};

mod draw;
mod draw_prepare;
mod embedded_frag;
mod geometry;
mod nontetrahedral;
mod raster;
mod templates;

#[cfg(test)]
mod disconnected_stereo_tests;

#[derive(Debug, Clone, PartialEq)]
pub enum DepictError {
    InvalidTopology(TopologyValidationError),
    PropertyCache(ValenceError),
    RingFinding(RingFindingError),
    StereoAssignment(LegacyStereoError),
    TemplateLoading(Coordinate2DTemplateError),
    Fragment(Coordinate2DLayoutError),
    CoordinateValidation(CoordinateValidationError),
    ConformerIdOverflow { max_id: usize },
    CoordGenUnavailable,
}

impl fmt::Display for DepictError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidTopology(source) => write!(formatter, "invalid topology: {source}"),
            Self::PropertyCache(detail) => {
                write!(formatter, "property-cache assignment failed: {detail}")
            }
            Self::RingFinding(detail) => write!(formatter, "ring finding failed: {detail}"),
            Self::StereoAssignment(detail) => {
                write!(formatter, "stereochemistry assignment failed: {detail}")
            }
            Self::TemplateLoading(detail) => write!(formatter, "template loading failed: {detail}"),
            Self::Fragment(detail) => write!(formatter, "fragment layout failed: {detail}"),
            Self::CoordinateValidation(detail) => {
                write!(formatter, "coordinate validation failed: {detail}")
            }
            Self::ConformerIdOverflow { max_id } => write!(
                formatter,
                "cannot append a 2D conformer after the maximum identifier {max_id}"
            ),
            Self::CoordGenUnavailable => formatter.write_str(
                "CoordGen is not available; request the RDKit depiction route explicitly",
            ),
        }
    }
}

impl Error for DepictError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::InvalidTopology(source) => Some(source),
            Self::PropertyCache(source) => Some(source),
            Self::RingFinding(source) => Some(source),
            Self::StereoAssignment(source) => Some(source),
            Self::TemplateLoading(source) => Some(source),
            Self::Fragment(source) => Some(source),
            Self::CoordinateValidation(source) => Some(source),
            Self::ConformerIdOverflow { .. } | Self::CoordGenUnavailable => None,
        }
    }
}

/// Public, typed projection of failures while loading coordinate templates.
///
/// The depiction implementation keeps its template loader private while
/// preserving its stable categories, row/path context, and typed causes.
#[derive(Debug, Clone)]
pub enum Coordinate2DTemplateError {
    InvalidDefaultRow {
        index: usize,
        source: cosmolkit_search::SmartsParseError,
    },
    InvalidTopology {
        index: usize,
        source: cosmolkit_model::TopologyValidationError,
    },
    NonElementIdentity {
        index: usize,
        source: cosmolkit_model::QueryAtomConversionError,
    },
    RingInitialization {
        index: usize,
        source: cosmolkit_core::RingFindingError,
    },
    ConnectedComponents {
        index: usize,
        source: cosmolkit_core::PathError,
    },
    ExternalOpen {
        path: String,
        source: Arc<std::io::Error>,
    },
    ExternalRead {
        path: String,
        source: Arc<std::io::Error>,
    },
    ExternalInvalidSmarts {
        path: String,
        line: usize,
        source: cosmolkit_search::SmartsParseError,
    },
    MissingCoordinates {
        row: String,
    },
    ThreeDimensionalCoordinates {
        row: String,
    },
    MultipleFragments {
        row: String,
    },
    NotRingSystem {
        row: String,
    },
}

impl PartialEq for Coordinate2DTemplateError {
    fn eq(&self, other: &Self) -> bool {
        use Coordinate2DTemplateError as Error;
        match (self, other) {
            (
                Error::InvalidDefaultRow {
                    index: left_index,
                    source: left_source,
                },
                Error::InvalidDefaultRow {
                    index: right_index,
                    source: right_source,
                },
            ) => left_index == right_index && left_source == right_source,
            (
                Error::InvalidTopology {
                    index: left_index,
                    source: left_source,
                },
                Error::InvalidTopology {
                    index: right_index,
                    source: right_source,
                },
            ) => left_index == right_index && left_source == right_source,
            (
                Error::NonElementIdentity {
                    index: left_index,
                    source: left_source,
                },
                Error::NonElementIdentity {
                    index: right_index,
                    source: right_source,
                },
            ) => left_index == right_index && left_source == right_source,
            (
                Error::RingInitialization {
                    index: left_index,
                    source: left_source,
                },
                Error::RingInitialization {
                    index: right_index,
                    source: right_source,
                },
            ) => left_index == right_index && left_source == right_source,
            (
                Error::ConnectedComponents {
                    index: left_index,
                    source: left_source,
                },
                Error::ConnectedComponents {
                    index: right_index,
                    source: right_source,
                },
            ) => left_index == right_index && left_source == right_source,
            (
                Error::ExternalOpen {
                    path: left_path,
                    source: left_source,
                },
                Error::ExternalOpen {
                    path: right_path,
                    source: right_source,
                },
            )
            | (
                Error::ExternalRead {
                    path: left_path,
                    source: left_source,
                },
                Error::ExternalRead {
                    path: right_path,
                    source: right_source,
                },
            ) => {
                left_path == right_path
                    && left_source.kind() == right_source.kind()
                    && left_source.raw_os_error() == right_source.raw_os_error()
                    && left_source.to_string() == right_source.to_string()
            }
            (
                Error::ExternalInvalidSmarts {
                    path: left_path,
                    line: left_line,
                    source: left_source,
                },
                Error::ExternalInvalidSmarts {
                    path: right_path,
                    line: right_line,
                    source: right_source,
                },
            ) => left_path == right_path && left_line == right_line && left_source == right_source,
            (Error::MissingCoordinates { row: left }, Error::MissingCoordinates { row: right })
            | (
                Error::ThreeDimensionalCoordinates { row: left },
                Error::ThreeDimensionalCoordinates { row: right },
            )
            | (Error::MultipleFragments { row: left }, Error::MultipleFragments { row: right })
            | (Error::NotRingSystem { row: left }, Error::NotRingSystem { row: right }) => {
                left == right
            }
            _ => false,
        }
    }
}

impl fmt::Display for Coordinate2DTemplateError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDefaultRow { index, source } => {
                write!(formatter, "invalid default template row {index}: {source}")
            }
            Self::InvalidTopology { index, source } => {
                write!(
                    formatter,
                    "invalid topology for template row {index}: {source}"
                )
            }
            Self::NonElementIdentity { index, source } => write!(
                formatter,
                "non-Element atom identity in template row {index}: {source}"
            ),
            Self::RingInitialization { index, source } => {
                write!(
                    formatter,
                    "ring initialization failed for template row {index}: {source}"
                )
            }
            Self::ConnectedComponents { index, source } => write!(
                formatter,
                "connected-component analysis failed for template row {index}: {source}"
            ),
            Self::ExternalOpen { path, source } => {
                write!(formatter, "could not open template file {path}: {source}")
            }
            Self::ExternalRead { path, source } => {
                write!(formatter, "could not read template file {path}: {source}")
            }
            Self::ExternalInvalidSmarts { path, line, source } => {
                write!(
                    formatter,
                    "invalid SMARTS in {path} at line {line}: {source}"
                )
            }
            Self::MissingCoordinates { row } => {
                write!(formatter, "template has no coordinates: {row}")
            }
            Self::ThreeDimensionalCoordinates { row } => {
                write!(formatter, "template coordinates are 3D: {row}")
            }
            Self::MultipleFragments { row } => {
                write!(formatter, "template has multiple fragments: {row}")
            }
            Self::NotRingSystem { row } => {
                write!(formatter, "template is not a ring system: {row}")
            }
        }
    }
}

impl Error for Coordinate2DTemplateError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::InvalidDefaultRow { source, .. } | Self::ExternalInvalidSmarts { source, .. } => {
                Some(source)
            }
            Self::InvalidTopology { source, .. } => Some(source),
            Self::NonElementIdentity { source, .. } => Some(source),
            Self::RingInitialization { source, .. } => Some(source),
            Self::ConnectedComponents { source, .. } => Some(source),
            Self::ExternalOpen { source, .. } | Self::ExternalRead { source, .. } => {
                Some(source.as_ref())
            }
            Self::MissingCoordinates { .. }
            | Self::ThreeDimensionalCoordinates { .. }
            | Self::MultipleFragments { .. }
            | Self::NotRingSystem { .. } => None,
        }
    }
}

impl From<TemplateError> for Coordinate2DTemplateError {
    fn from(error: TemplateError) -> Self {
        match error {
            TemplateError::InvalidDefaultRow { index, source } => {
                Self::InvalidDefaultRow { index, source }
            }
            TemplateError::InvalidTopology { index, source } => {
                Self::InvalidTopology { index, source }
            }
            TemplateError::NonElementIdentity { index, source } => {
                Self::NonElementIdentity { index, source }
            }
            TemplateError::RingInitialization { index, source } => {
                Self::RingInitialization { index, source }
            }
            TemplateError::ConnectedComponents { index, source } => {
                Self::ConnectedComponents { index, source }
            }
            TemplateError::ExternalOpen { path, source } => Self::ExternalOpen {
                path,
                source: Arc::new(source),
            },
            TemplateError::ExternalRead { path, source } => Self::ExternalRead {
                path,
                source: Arc::new(source),
            },
            TemplateError::ExternalInvalidSmarts { path, line, source } => {
                Self::ExternalInvalidSmarts { path, line, source }
            }
            TemplateError::MissingCoordinates { row } => Self::MissingCoordinates { row },
            TemplateError::ThreeDimensionalCoordinates { row } => {
                Self::ThreeDimensionalCoordinates { row }
            }
            TemplateError::MultipleFragments { row } => Self::MultipleFragments { row },
            TemplateError::NotRingSystem { row } => Self::NotRingSystem { row },
        }
    }
}

/// Public, typed projection of failures while arranging detached fragments.
///
/// Private depiction geometry details are represented by their stable fields;
/// already-public search, path, and matrix causes remain source-chain leaves.
#[derive(Debug, Clone, PartialEq)]
pub enum Coordinate2DLayoutError {
    AtomIndexOutOfRange {
        atom: usize,
        atom_count: usize,
    },
    AtomAlreadyEmbedded {
        atom: usize,
    },
    AtomNotEmbedded {
        atom: usize,
    },
    NotEnoughEmbeddedNeighbors {
        atom: usize,
        count: usize,
    },
    CoincidentPoints,
    InvalidAngle,
    NoCommonAtoms,
    MismatchedTopology,
    EmptyAttachment {
        atom: usize,
    },
    NonTetrahedralNoLigand {
        centre: usize,
    },
    NonTetrahedralLigandOverflow {
        centre: usize,
    },
    CisTransBondInvalid {
        bond: usize,
    },
    CollisionBondInvalid {
        bond: usize,
    },
    UndefinedSamplingDistance {
        first: usize,
        second: usize,
    },
    AtomCountTooLarge {
        atom_count: usize,
    },
    InvalidRankProperty {
        atom: usize,
        key: &'static str,
        source: cosmolkit_core::PropertyUIntReadError,
    },
    GeometryNotEnoughNeighbors {
        atom: usize,
        count: usize,
    },
    TemplateMatch(cosmolkit_search::SubstructMatchError),
    GraphPath(cosmolkit_core::PathError),
    GraphDistance(cosmolkit_core::MatrixError),
}

impl fmt::Display for Coordinate2DLayoutError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::AtomIndexOutOfRange { atom, atom_count } => {
                write!(
                    formatter,
                    "atom {atom} is out of range for {atom_count} atoms"
                )
            }
            Self::AtomAlreadyEmbedded { atom } => {
                write!(formatter, "atom {atom} is already embedded")
            }
            Self::AtomNotEmbedded { atom } => write!(formatter, "atom {atom} is not embedded"),
            Self::NotEnoughEmbeddedNeighbors { atom, count } => {
                write!(formatter, "atom {atom} has only {count} embedded neighbors")
            }
            Self::CoincidentPoints => formatter.write_str("reference points are coincident"),
            Self::InvalidAngle => formatter.write_str("invalid fragment angle"),
            Self::NoCommonAtoms => formatter.write_str("fragments have no common atoms"),
            Self::MismatchedTopology => formatter.write_str("fragment topology does not match"),
            Self::EmptyAttachment { atom } => {
                write!(formatter, "atom {atom} has an empty attachment")
            }
            Self::NonTetrahedralNoLigand { centre } => {
                write!(formatter, "non-tetrahedral center {centre} has no ligand")
            }
            Self::NonTetrahedralLigandOverflow { centre } => {
                write!(
                    formatter,
                    "non-tetrahedral center {centre} has too many ligands"
                )
            }
            Self::CisTransBondInvalid { bond } => {
                write!(formatter, "invalid cis/trans bond {bond}")
            }
            Self::CollisionBondInvalid { bond } => {
                write!(formatter, "invalid collision bond {bond}")
            }
            Self::UndefinedSamplingDistance { first, second } => write!(
                formatter,
                "sampling distance for atom pair {first}-{second} is undefined"
            ),
            Self::AtomCountTooLarge { atom_count } => {
                write!(
                    formatter,
                    "atom count {atom_count} exceeds depiction limits"
                )
            }
            Self::InvalidRankProperty { atom, key, source } => {
                write!(
                    formatter,
                    "atom {atom} has invalid rank property {key}: {source}"
                )
            }
            Self::GeometryNotEnoughNeighbors { atom, count } => {
                write!(formatter, "atom {atom} has only {count} geometry neighbors")
            }
            Self::TemplateMatch(source) => write!(formatter, "template matching failed: {source}"),
            Self::GraphPath(source) => {
                write!(formatter, "fragment path calculation failed: {source}")
            }
            Self::GraphDistance(source) => {
                write!(formatter, "fragment distance calculation failed: {source}")
            }
        }
    }
}

impl Error for Coordinate2DLayoutError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::InvalidRankProperty { source, .. } => Some(source),
            Self::TemplateMatch(source) => Some(source),
            Self::GraphPath(source) => Some(source),
            Self::GraphDistance(source) => Some(source),
            Self::AtomIndexOutOfRange { .. }
            | Self::AtomAlreadyEmbedded { .. }
            | Self::AtomNotEmbedded { .. }
            | Self::NotEnoughEmbeddedNeighbors { .. }
            | Self::CoincidentPoints
            | Self::InvalidAngle
            | Self::NoCommonAtoms
            | Self::MismatchedTopology
            | Self::EmptyAttachment { .. }
            | Self::NonTetrahedralNoLigand { .. }
            | Self::NonTetrahedralLigandOverflow { .. }
            | Self::CisTransBondInvalid { .. }
            | Self::CollisionBondInvalid { .. }
            | Self::UndefinedSamplingDistance { .. }
            | Self::AtomCountTooLarge { .. }
            | Self::GeometryNotEnoughNeighbors { .. } => None,
        }
    }
}

impl From<GeometryError> for Coordinate2DLayoutError {
    fn from(error: GeometryError) -> Self {
        match error {
            GeometryError::AtomIndexOutOfRange { atom, atom_count } => {
                Self::AtomIndexOutOfRange { atom, atom_count }
            }
            GeometryError::AtomCountTooLarge { atom_count } => {
                Self::AtomCountTooLarge { atom_count }
            }
            GeometryError::InvalidRankProperty { atom, key, source } => {
                Self::InvalidRankProperty { atom, key, source }
            }
            GeometryError::NotEnoughNeighbors { atom, count } => {
                Self::GeometryNotEnoughNeighbors { atom, count }
            }
        }
    }
}

impl From<FragmentError> for Coordinate2DLayoutError {
    fn from(error: FragmentError) -> Self {
        match error {
            FragmentError::AtomIndexOutOfRange { atom, atom_count } => {
                Self::AtomIndexOutOfRange { atom, atom_count }
            }
            FragmentError::AtomAlreadyEmbedded { atom } => Self::AtomAlreadyEmbedded { atom },
            FragmentError::AtomNotEmbedded { atom } => Self::AtomNotEmbedded { atom },
            FragmentError::NotEnoughEmbeddedNeighbors { atom, count } => {
                Self::NotEnoughEmbeddedNeighbors { atom, count }
            }
            FragmentError::CoincidentPoints => Self::CoincidentPoints,
            FragmentError::InvalidAngle => Self::InvalidAngle,
            FragmentError::NoCommonAtoms => Self::NoCommonAtoms,
            FragmentError::MismatchedTopology => Self::MismatchedTopology,
            FragmentError::EmptyAttachment { atom } => Self::EmptyAttachment { atom },
            FragmentError::NonTetrahedralNoLigand { centre } => {
                Self::NonTetrahedralNoLigand { centre }
            }
            FragmentError::NonTetrahedralLigandOverflow { centre } => {
                Self::NonTetrahedralLigandOverflow { centre }
            }
            FragmentError::CisTransBondInvalid { bond } => Self::CisTransBondInvalid { bond },
            FragmentError::CollisionBondInvalid { bond } => Self::CollisionBondInvalid { bond },
            FragmentError::UndefinedSamplingDistance { first, second } => {
                Self::UndefinedSamplingDistance { first, second }
            }
            FragmentError::Geometry(error) => error.into(),
            FragmentError::TemplateMatch(source) => Self::TemplateMatch(source),
            FragmentError::GraphPath(source) => Self::GraphPath(source),
            FragmentError::GraphDistance(source) => Self::GraphDistance(source),
        }
    }
}

impl From<FragmentError> for DepictError {
    fn from(value: FragmentError) -> Self {
        Self::Fragment(value.into())
    }
}

impl From<TemplateError> for DepictError {
    fn from(value: TemplateError) -> Self {
        Self::TemplateLoading(value.into())
    }
}

/// Frozen detached counterpart of RDKit's `Compute2DCoordParameters`.
///
/// `clear_existing_2d` is consumed by the live runtime when this detached
/// result is installed. It deliberately does not affect independent 3D rows.
#[derive(Debug, Clone, PartialEq)]
pub struct Compute2DCoordinatesParams {
    pub coordinate_map: BTreeMap<usize, [f64; 2]>,
    pub canonical_orientation: bool,
    pub clear_existing_2d: bool,
    pub flips_per_sample: u32,
    pub samples: u32,
    pub sample_seed: i32,
    pub permute_degree_four: bool,
    pub force_rdkit: bool,
    pub use_ring_templates: bool,
}

impl Default for Compute2DCoordinatesParams {
    fn default() -> Self {
        Self {
            coordinate_map: BTreeMap::new(),
            canonical_orientation: false,
            clear_existing_2d: true,
            flips_per_sample: 0,
            samples: 0,
            sample_seed: 0,
            permute_degree_four: false,
            force_rdkit: false,
            use_ring_templates: false,
        }
    }
}

fn largest_unfinished_fragment(fragments: &[EmbeddedFrag<'_>]) -> Option<usize> {
    // RDKit❗✔️: std::list<EmbeddedFrag>::iterator _findLargestFrag(
    // RDKit❗✔️:     std::list<EmbeddedFrag> &efrags) {
    // RDKit❗✔️:   std::list<EmbeddedFrag>::iterator mfri;
    // RDKit❗✔️:   int msiz = 0;
    // RDKit❗✔️:   for (auto efri = efrags.begin(); efri != efrags.end(); ++efri) {
    // RDKit❗✔️:     if ((!efri->isDone()) && (efri->Size() > msiz)) {
    // RDKit❗✔️:       msiz = efri->Size();
    // RDKit❗✔️:       mfri = efri;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (msiz == 0) { mfri = efrags.end(); }
    // RDKit❗✔️:   return mfri;
    // RDKit❗✔️: }
    // Behavior review: strict-greater replacement retains the first largest
    // unfinished fragment in source list order.
    // Complexity review: one linear scan without allocation, matching source.
    let mut selected = None;
    let mut size = 0;
    for (index, fragment) in fragments.iter().enumerate() {
        if !fragment.done && fragment.atoms.len() > size {
            selected = Some(index);
            size = fragment.atoms.len();
        }
    }
    selected
}

fn compute_initial_coordinates<'a>(
    topology: &'a TopologyBlock,
    rings: &'a cosmolkit_core::RingInfo,
    coordinate_map: Option<&PointMap>,
    use_ring_templates: bool,
) -> Result<Vec<EmbeddedFrag<'a>>, DepictError> {
    // BEGIN RECOVERY DEP-04 SOURCE compute_initial_coordinates
    // RDKit❗❌: void computeInitialCoords(RDKit::ROMol &mol,
    // RDKit❗❌:                           const RDGeom::INT_POINT2D_MAP *coordMap,
    // RDKit❗❌:                           std::list<EmbeddedFrag> &efrags,
    // RDKit❗❌:                           bool useRingTemplates) {
    // RDKit❗❌:   std::vector<int> atomRanks;
    // RDKit❗❌:   atomRanks.resize(mol.getNumAtoms());
    // RDKit❗❌:   for (auto i = 0u; i < mol.getNumAtoms(); ++i) {
    // RDKit❗❌:     atomRanks[i] = getAtomDepictRank(mol.getAtomWithIdx(i));
    // RDKit❗❌:   }
    // RDKit❗❌:   RDKit::VECT_INT_VECT arings;
    // RDKit❗❌:
    // RDKit❗❌:   // first find all the rings
    // RDKit❗❌:   bool includeDativeBonds = true;
    // RDKit❗❌:   RDKit::MolOps::symmetrizeSSSR(mol, arings, includeDativeBonds);
    // RDKit❗❌:
    // RDKit❗❌:   // do stereochemistry
    // RDKit❗❌:   RDKit::MolOps::assignStereochemistry(mol, false);
    // RDKit❗❌:
    // RDKit❗❌:   efrags.clear();
    // RDKit❗❌:
    // RDKit❗❌:   // user-specified coordinates exist
    // RDKit❗❌:   bool preSpec = false;
    // RDKit❗❌:   // first embed any atoms for which the coordinates have been specified.
    // RDKit❗❌:   if ((coordMap) && (coordMap->size() > 1)) {
    // RDKit❗❌:     EmbeddedFrag efrag(&mol, *coordMap);
    // RDKit❗❌:     // add this to the list of embedded fragments
    // RDKit❗❌:     efrags.push_back(efrag);
    // RDKit❗❌:     preSpec = true;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (arings.size() > 0) {
    // RDKit❗❌:     // first deal with the fused rings
    // RDKit❗❌:     DepictorLocal::embedFusedSystems(mol, arings, efrags, coordMap,
    // RDKit❗❌:                                      useRingTemplates);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // do non-tetrahedral stereo
    // RDKit❗❌:   DepictorLocal::embedNontetrahedralStereo(mol, efrags, atomRanks);
    // RDKit❗❌:
    // RDKit❗❌:   // deal with any cis/trans systems
    // RDKit❗❌:   DepictorLocal::embedCisTransSystems(mol, efrags);
    // RDKit❗❌:   // now get the atoms that are not yet embedded in either a cis/trans system
    // RDKit❗❌:   // or a ring system (or simply the first atom)
    // RDKit❗❌:   auto nratms = DepictorLocal::getNonEmbeddedAtoms(mol, efrags);
    // RDKit❗❌:   std::list<EmbeddedFrag>::iterator mri;
    // RDKit❗❌:   if (preSpec) {
    // RDKit❗❌:     // if the user specified coordinates on some of the atoms use that as
    // RDKit❗❌:     // as the starting fragment and it should be at the beginning of the vector
    // RDKit❗❌:     mri = efrags.begin();
    // RDKit❗❌:   } else {
    // RDKit❗❌:     // otherwise - find the largest fragment that was embedded
    // RDKit❗❌:     mri = DepictorLocal::_findLargestFrag(efrags);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   while ((mri != efrags.end()) || (nratms.size() > 0)) {
    // RDKit❗❌:     if (mri == efrags.end()) {
    // RDKit❗❌:       // we are out of embedded fragments, if there are any
    // RDKit❗❌:       // non embedded atoms use them to start a fragment
    // RDKit❗❌:       auto mrank = RDKit::MAX_INT;
    // RDKit❗❌:       auto mnri = nratms.end();
    // RDKit❗❌:       for (auto nri = nratms.begin(); nri != nratms.end(); ++nri) {
    // RDKit❗❌:         auto rank = atomRanks.at(*nri);
    // RDKit❗❌:         rank *= mol.getNumAtoms();
    // RDKit❗❌:         // use the atom index as well so that we at least
    // RDKit❗❌:         // get reproducible depictions in cases where things
    // RDKit❗❌:         // have identical ranks.
    // RDKit❗❌:         rank += *nri;
    // RDKit❗❌:         if (rank < mrank) {
    // RDKit❗❌:           mrank = rank;
    // RDKit❗❌:           mnri = nri;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       EmbeddedFrag efrag((*mnri), &mol);
    // RDKit❗❌:       nratms.erase(mnri);
    // RDKit❗❌:       efrags.push_back(efrag);
    // RDKit❗❌:       mri = efrags.end();
    // RDKit❗❌:       --mri;
    // RDKit❗❌:     }
    // RDKit❗❌:     mri->markDone();
    // RDKit❗❌:     mri->expandEfrag(nratms, efrags);
    // RDKit❗❌:     mri = DepictorLocal::_findLargestFrag(efrags);
    // RDKit❗❌:   }
    // RDKit❗❌:   // at this point any remaining efrags should belong individual fragments in
    // RDKit❗❌:   // the molecule
    // RDKit❗❌: }
    // END RECOVERY DEP-04 SOURCE compute_initial_coordinates

    // Behavior review: the caller supplies source-shaped cache, SymmSSSR and
    // cleanIt=false stereo preparation. This preserves helper append order,
    // first-largest ties, wrapped signed ranks and disconnected expansion.
    // Complexity review: graph/rank scans and fragment expansion match source;
    // Vec removal shifts handles but adds no graph-scale nested traversal.
    let atom_ranks = (0..topology.atoms.len())
        .map(|atom| atom_depict_rank(topology, atom))
        .collect::<Result<Vec<_>, _>>()
        .map_err(FragmentError::from)?;

    let mut fragments = Vec::new();
    let prespecified = seed_coordinate_constraints(topology, rings, coordinate_map)?;
    let has_prespecified = prespecified.is_some();
    if let Some(fragment) = prespecified {
        fragments.push(fragment);
    }
    let mut templates = if use_ring_templates {
        CoordinateTemplates::new_default()?
    } else {
        CoordinateTemplates::default()
    };
    fragments.extend(embed_fused_systems(
        topology,
        rings,
        coordinate_map,
        use_ring_templates,
        &mut templates,
    )?);

    fragments.extend(embed_nontetrahedral_stereo(topology, rings)?);

    fragments.extend(embed_cis_trans_systems(topology, rings)?);

    let mut embedded = vec![false; topology.atoms.len()];
    for fragment in &fragments {
        for &atom in fragment.atoms.keys() {
            embedded[atom] = true;
        }
    }
    let mut nonembedded: Vec<_> = embedded
        .iter()
        .enumerate()
        .filter_map(|(atom, &done)| (!done).then_some(atom))
        .collect();

    let mut selected = has_prespecified
        .then_some(0)
        .or_else(|| largest_unfinished_fragment(&fragments));
    while selected.is_some() || !nonembedded.is_empty() {
        let index = if let Some(index) = selected {
            index
        } else {
            let atom_count = topology.atoms.len() as i32;
            let (position, &atom) = nonembedded
                .iter()
                .enumerate()
                .min_by_key(|(_, atom)| {
                    atom_ranks[**atom]
                        .wrapping_mul(atom_count)
                        .wrapping_add(**atom as i32)
                })
                .expect("loop requires a nonempty atom list");
            let seed_stage = format!("initial_seed_atom_{atom}");

            let seed = EmbeddedFrag::from_single(atom, topology, rings)?;

            nonembedded.remove(position);

            fragments.push(seed);

            fragments.len() - 1
        };
        let mut fragment = fragments.remove(index);

        fragment.done = true;

        // RDKit✔️❌: mri->markDone();
        // RDKit✔️❌: mri->expandEfrag(nratms, efrags);
        // A list iterator retains its position when preceding fragments are
        // erased. Track that position while the selected Vec element is moved
        // out for exclusive access; reinserting at its old index changes the
        // disconnected-component packing order. Vec erasure still shifts rows,
        // unlike the source list's constant-time erasure.
        let mut insertion_index = index;
        fragment.expand_fragment(&mut nonembedded, &mut fragments, &mut insertion_index)?;

        fragments.insert(insertion_index, fragment);
        selected = largest_unfinished_fragment(&fragments);
    }

    Ok(fragments)
}

/// Compute one detached atom-ordered 2D conformer.
pub fn compute_2d_coordinates(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    params: &Compute2DCoordinatesParams,
) -> Result<Conformer2D, DepictError> {
    // BEGIN RECOVERY DEP-04 SOURCE compute_2d_coordinates
    // RDKit❗❌: unsigned int compute2DCoords(RDKit::ROMol &mol,
    // RDKit❗❌:                              const Compute2DCoordParameters &params) {
    // RDKit❗❌:   if (mol.needsUpdatePropertyCache()) {
    // RDKit❗❌:     mol.updatePropertyCache(false);
    // RDKit❗❌:   }
    // RDKit❗❌: #ifdef RDK_BUILD_COORDGEN_SUPPORT
    // RDKit❗❌:   // default to use CoordGen if we have it installed
    // RDKit❗❌:   if (!params.forceRDKit && preferCoordGen) {
    // RDKit❗❌:     RDKit::CoordGen::CoordGenParams coordgen_params;
    // RDKit❗❌:     if (params.coordMap) {
    // RDKit❗❌:       coordgen_params.coordMap = *params.coordMap;
    // RDKit❗❌:     }
    // RDKit❗❌:     auto cid = RDKit::CoordGen::addCoords(mol, &coordgen_params);
    // RDKit❗❌:     return cid;
    // RDKit❗❌:   };
    // RDKit❗❌: #endif
    // RDKit❗❌:
    // RDKit❗❌:   RDKit::ROMol cp(mol);
    // RDKit❗❌:   // storage for pieces of a molecule/s that are embedded in 2D
    // RDKit❗❌:   std::list<EmbeddedFrag> efrags;
    // RDKit❗❌:   computeInitialCoords(cp, params.coordMap, efrags, params.useRingTemplates);
    // RDKit❗❌:
    // RDKit❗❌: #if 1
    // RDKit❗❌:   // perform random sampling here to improve the density
    // RDKit❗❌:   for (auto &eri : efrags) {
    // RDKit❗❌:     // either sample the 2D space by randomly flipping rotatable
    // RDKit❗❌:     // bonds in the structure or flip only bonds along the shortest
    // RDKit❗❌:     // path between colliding atoms - don't do both
    // RDKit❗❌:     if ((params.nSamples > 0) && (params.nFlipsPerSample > 0)) {
    // RDKit❗❌:       eri.randomSampleFlipsAndPermutations(
    // RDKit❗❌:           params.nFlipsPerSample, params.nSamples, params.sampleSeed, nullptr,
    // RDKit❗❌:           0.0, params.permuteDeg4Nodes);
    // RDKit❗❌:     } else {
    // RDKit❗❌:       eri.removeCollisionsBondAndSpiroFlip();
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   for (auto &eri : efrags) {
    // RDKit❗❌:     // if there are any remaining collisions
    // RDKit❗❌:     eri.removeCollisionsOpenAngles();
    // RDKit❗❌:     eri.removeCollisionsShortenBonds();
    // RDKit❗❌:   }
    // RDKit❗❌:   if (!params.coordMap || !params.coordMap->size()) {
    // RDKit❗❌:     if (params.canonOrient && efrags.size()) {
    // RDKit❗❌:       // if we do not have any prespecified coordinates - canonicalize
    // RDKit❗❌:       // the orientation of the fragment so that the longest axes fall
    // RDKit❗❌:       // along the x-axis etc.
    // RDKit❗❌:       for (auto &eri : efrags) {
    // RDKit❗❌:         eri.canonicalizeOrientation();
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   DepictorLocal::_shiftCoords(efrags);
    // RDKit❗❌: #endif
    // RDKit❗❌:   // create a conformation on the molecule and copy the coordinates
    // RDKit❗❌:   auto cid = copyCoordinate(mol, efrags, params.clearConfs);
    // RDKit❗❌:
    // RDKit❗❌:   // special case for a single-atom coordMap template
    // RDKit❗❌:   if ((params.coordMap) && (params.coordMap->size() == 1)) {
    // RDKit❗❌:     auto &conf = mol.getConformer(cid);
    // RDKit❗❌:     auto cRef = params.coordMap->begin();
    // RDKit❗❌:     const auto &confPos = conf.getAtomPos(cRef->first);
    // RDKit❗❌:     auto refPos = cRef->second;
    // RDKit❗❌:     refPos.x -= confPos.x;
    // RDKit❗❌:     refPos.y -= confPos.y;
    // RDKit❗❌:     for (auto i = 0u; i < conf.getNumAtoms(); ++i) {
    // RDKit❗❌:       auto confPos = conf.getAtomPos(i);
    // RDKit❗❌:       confPos.x += refPos.x;
    // RDKit❗❌:       confPos.y += refPos.y;
    // RDKit❗❌:       conf.setAtomPos(i, confPos);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return cid;
    // RDKit❗❌: }
    // END RECOVERY DEP-04 SOURCE compute_2d_coordinates

    // Behavior review: CK fixes `preferCoordGen=false`, the pinned upstream
    // default; no silent alternate algorithm is selected. `force_rdkit` is
    // consequently accepted but does not change this fixed route. Runtime
    // alone applies dimension-separated clear/append semantics.
    // Complexity review: one temporary topology clone and final coordinate
    // allocation match source clone/copy shape; helper costs are unchanged.
    topology.validate().map_err(DepictError::InvalidTopology)?;
    let valence = assign_valence(
        topology,
        &ValenceParams {
            model: ValenceModel::RdkitLike,
            strict: false,
        },
    )
    .map_err(DepictError::PropertyCache)?;

    let rings = symmetrized_sssr(
        topology,
        &RingSearchParams {
            include_dative_bonds: true,
            include_hydrogen_bonds: false,
        },
    )
    .map_err(DepictError::RingFinding)?;
    let copied_topology = topology.clone();

    // Pinned Chirality.cpp:2889-2905; MolOps.h defaults force=false and
    // flagPossibleStereoCenters=false. The fixed source profile is legacy.
    // Behavior: presence alone returns the existing copied topology. Values,
    // types and computed membership do not control this branch. The existing
    // core owner handles the absent branch; its copy-local done publication
    // is not returned or published to the borrowed live properties. Existing
    // valence/ring preparation ordering and absent cache validity remain
    // qualified separately; this is not complete ROMol prepared-state parity.
    // Cost: borrowed BTreeMap lookup is O(log P) versus source Dict's O(P)
    // presence scan, with no allocation; the same single clone is moved into
    // either branch. No chemistry or extra preparation is duplicated.
    let working = if properties.prop("_StereochemDone").is_some() {
        copied_topology
    } else {
        #[cfg(test)]
        d2_probe1_tests::prepared_properties::observe_dispatch();
        assign_legacy_stereochemistry_for_depiction(copied_topology, &valence, &rings)
            .map_err(DepictError::StereoAssignment)?
    };
    #[cfg(test)]
    d2_probe1_tests::prepared_properties::observe_working(&working);
    let coordinate_map = Some(&params.coordinate_map);
    let mut fragments =
        compute_initial_coordinates(&working, &rings, coordinate_map, params.use_ring_templates)?;

    for fragment in &mut fragments {
        if params.samples > 0 && params.flips_per_sample > 0 {
            fragment.random_sample_flips_and_permutations(
                params.flips_per_sample,
                params.samples,
                params.sample_seed,
                None,
                0.0,
                params.permute_degree_four,
            )?;
        } else {
            fragment.remove_collisions_bond_and_spiro_flip()?;
        }
    }

    for fragment_index in 0..fragments.len() {
        fragments[fragment_index].remove_collisions_open_angles()?;

        fragments[fragment_index].remove_collisions_shorten_bonds()?;
    }
    orient_and_shift_fragments(
        &mut fragments,
        params.canonical_orientation,
        Some(params.coordinate_map.len()),
    );

    translate_single_coordinate_constraint(&working, &mut fragments, coordinate_map)?;

    // Behavior review: rows start at source default origin and are overwritten
    // in fragment/map order; detached ID zero is provisional for the runtime.
    // Complexity review: one atom-sized allocation and ordered fragment copy.
    let mut coordinates = vec![[0.0, 0.0]; working.atoms.len()];
    for fragment in &fragments {
        for (&atom, embedded) in &fragment.atoms {
            coordinates[atom] = embedded.loc;
        }
    }
    let conformer = Conformer2D::new(0, coordinates);
    conformer
        .validate_for_atom_count(working.atoms.len())
        .map_err(DepictError::CoordinateValidation)?;

    Ok(conformer)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DepictOptions {
    pub width: u32,
    pub height: u32,
}

#[deprecated(note = "use compute_2d_coordinates")]
pub fn layout_2d(
    topology: &TopologyBlock,
    _options: &DepictOptions,
) -> Result<CoordinateBlock, DepictError> {
    // This deprecated topology-only helper explicitly has absent MOL properties.
    let absent_properties = MoleculeProperties::default();
    let conformer = compute_2d_coordinates(
        topology,
        &absent_properties,
        &Compute2DCoordinatesParams::default(),
    )?;
    Ok(CoordinateBlock {
        conformers_2d: vec![conformer],
        ..Default::default()
    })
}

pub use draw::DrawingError;

/// Borrowed input to detached drawing preparation. Stored 2D layouts take
/// priority; missing 2D is generated without replacing stored 3D conformers.
pub struct DrawingInput<'a> {
    pub topology: &'a TopologyBlock,
    pub coordinates: &'a CoordinateBlock,
    pub properties: &'a cosmolkit_model::MoleculeProperties,
    pub valence: Option<&'a cosmolkit_core::ValenceAssignment>,
    pub rings: Option<&'a cosmolkit_core::RingInfo>,
}

pub fn render_svg(
    input: DrawingInput<'_>,
    options: &DepictOptions,
) -> Result<Vec<u8>, DrawingError> {
    render_svg_with_identity(input, options, draw::SvgIdentity::PinnedSource)
}

/// Render the public COSMolKit SVG product with its declared website namespace.
/// Scene preparation, geometry and glyph emission use the same detached owner.
pub fn render_cosmolkit_svg(
    input: DrawingInput<'_>,
    options: &DepictOptions,
) -> Result<Vec<u8>, DrawingError> {
    render_svg_with_identity(input, options, draw::SvgIdentity::Cosmolkit)
}

fn render_svg_with_identity(
    input: DrawingInput<'_>,
    options: &DepictOptions,
    identity: draw::SvgIdentity,
) -> Result<Vec<u8>, DrawingError> {
    if options.width == 0 || options.height == 0 {
        return Err(DrawingError::InvalidDimensions {
            width: options.width,
            height: options.height,
        });
    }
    let prepared = draw_prepare::prepare(input)?;
    draw::render_prepared_svg_with_identity(
        &prepared.borrow(),
        options.width,
        options.height,
        identity,
    )
}

pub fn render_png(
    input: DrawingInput<'_>,
    options: &DepictOptions,
) -> Result<Vec<u8>, DrawingError> {
    let svg = render_svg(input, options)?;
    raster::svg_to_png(&svg)
}

#[cfg(test)]
mod error_projection_tests {
    use super::*;
    use cosmolkit_model::QueryAtomConversionError;

    #[test]
    fn non_element_template_identity_projects_typed_row_error() {
        let source = QueryAtomConversionError::NonElementAtomicNumber {
            atom: cosmolkit_model::AtomId::new(0),
            atomic_number: 119,
        };
        let projected = crate::Coordinate2DTemplateError::from(TemplateError::NonElementIdentity {
            index: 7,
            source,
        });
        assert_eq!(
            projected,
            crate::Coordinate2DTemplateError::NonElementIdentity { index: 7, source }
        );

        let same = crate::Coordinate2DTemplateError::from(TemplateError::NonElementIdentity {
            index: 7,
            source,
        });
        let different_row =
            crate::Coordinate2DTemplateError::from(TemplateError::NonElementIdentity {
                index: 8,
                source,
            });
        let different_source =
            crate::Coordinate2DTemplateError::from(TemplateError::NonElementIdentity {
                index: 7,
                source: QueryAtomConversionError::NonElementAtomicNumber {
                    atom: cosmolkit_model::AtomId::new(0),
                    atomic_number: 118,
                },
            });
        assert_eq!(projected, same);
        assert_ne!(projected, different_row);
        assert_ne!(projected, different_source);
        assert!(projected.to_string().contains("template row 7"));
        assert!(projected.to_string().contains("119"));

        let public_source = std::error::Error::source(&projected)
            .expect("typed conversion cause is retained as the public source");
        assert_eq!(
            public_source.downcast_ref::<QueryAtomConversionError>(),
            Some(&source)
        );
    }
}

#[cfg(test)]
mod d2_probe1_tests {
    use std::str::FromStr;

    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondDirection, BondId, BondOrder, BondSpec, BondStereo,
        ChiralTag, Element, Hybridization, MoleculeProperties, PropertyValue, TopologyBlock,
    };

    use super::{Compute2DCoordinatesParams, compute_2d_coordinates};

    pub(super) mod prepared_properties {
        use super::*;
        use std::cell::{Cell, RefCell};
        use std::collections::BTreeMap;

        thread_local! {
            static DISPATCHES: Cell<usize> = const { Cell::new(0) };
            static CAPTURE: RefCell<(bool, Option<TopologyBlock>)> = const { RefCell::new((false, None)) };
        }

        pub(crate) fn observe_dispatch() {
            DISPATCHES.with(|count| count.set(count.get() + 1));
        }

        pub(crate) fn observe_working(topology: &TopologyBlock) {
            CAPTURE.with(|capture| {
                let mut capture = capture.borrow_mut();
                if capture.0 {
                    capture.1 = Some(topology.clone());
                }
            });
        }

        fn captured_topology(block: &str) -> TopologyBlock {
            let trace = format!(
                "D2\tconstructor_final\tMETA\tfixed\n{}",
                block
                    .lines()
                    .map(|line| format!("D2\tconstructor_final\t{line}\n"))
                    .collect::<String>()
            );
            super::last_boundary::transport(&trace).0
        }

        fn captured_properties(block: &str) -> MoleculeProperties {
            let mut properties = MoleculeProperties::default();
            for line in block.lines() {
                let row = line.split('\t').collect::<Vec<_>>();
                assert_eq!(row.len(), 8);
                assert_eq!(&row[..2], &["PROP", "MOL"]);
                let key = super::last_boundary::unhex(row[4]);
                if key == "__computedProps" {
                    assert_eq!(row[5], "12");
                    // Vector representation remains raw source metadata. Scalar
                    // computed membership is taken from each actual capture.
                    continue;
                }
                assert!(matches!(row[5], "1" | "3"));
                let value = super::last_boundary::unhex(row[7]);
                if row[6] == "1" {
                    properties.set_computed_prop(key, value).unwrap();
                } else {
                    properties.set_prop(key, value).unwrap();
                }
            }
            properties
        }

        #[test]
        fn d2_prepared_property_domain_all_48_native_cells_twice() {
            let fixed = include_str!(
                "../../../testdata/depict_2d/expected/rdkit/prepared_property_presence.tsv"
            );
            let mut blocks = BTreeMap::new();
            for section in fixed.split("BLOCK\t").skip(1) {
                let (key, content) = section.split_once('\n').unwrap();
                let (content, _) = content.split_once("END\n").unwrap();
                assert!(blocks.insert(key, content).is_none());
            }
            let cells = fixed
                .lines()
                .filter(|line| line.starts_with("CELL\t"))
                .collect::<Vec<_>>();
            assert_eq!(cells.len(), 48);
            assert_eq!(blocks.len(), 39);
            let mut discrepancies = Vec::new();
            let mut calls = 0;
            let mut errors = 0;
            let mut preserved = 0;
            let mut dispatch_matches = 0;
            let mut working_matches = 0;
            let mut coordinate_matches = 0;
            let mut equal_bits = 0;
            let mut scalar_bits = 0;
            CAPTURE.with(|capture| capture.borrow_mut().0 = true);
            for cell in cells {
                let row = cell.split('\t').collect::<Vec<_>>();
                assert_eq!(row.len(), 7);
                let label = row[1];
                let expected_working = captured_topology(blocks[row[4]]);
                let expected_xy = blocks[row[6]]
                    .lines()
                    .map(|line| {
                        let fields = line.split('\t').collect::<Vec<_>>();
                        assert_eq!(fields.len(), 4);
                        assert_eq!(fields[0], "XY");
                        (
                            parse::<usize>(fields[1]),
                            [parse::<u64>(fields[2]), parse::<u64>(fields[3])],
                        )
                    })
                    .collect::<Vec<_>>();
                for repeat in 0..2 {
                    let topology = captured_topology(blocks[row[2]]);
                    let properties = captured_properties(blocks[row[3]]);
                    let params = Compute2DCoordinatesParams {
                        canonical_orientation: label.ends_with(":O1"),
                        ..Default::default()
                    };
                    let before_topology = topology.clone();
                    let before_properties = properties.clone();
                    let before_params = params.clone();
                    let before_values = format!("{topology:?}|{properties:?}|{params:?}");
                    let before_addresses = (
                        topology.atoms.as_ptr(),
                        topology.bonds.as_ptr(),
                        properties.props() as *const _,
                        &params as *const _,
                    );
                    let before_map_bits = params
                        .coordinate_map
                        .iter()
                        .map(|(id, xy)| (*id, xy.map(f64::to_bits)))
                        .collect::<Vec<_>>();
                    let baseline = DISPATCHES.with(Cell::get);
                    let result = compute_2d_coordinates(&topology, &properties, &params);
                    calls += 1;
                    // Inspect input preservation BEFORE inspecting any Result.
                    let preservation = topology == before_topology
                        && properties == before_properties
                        && params == before_params
                        && before_values == format!("{topology:?}|{properties:?}|{params:?}")
                        && before_addresses
                            == (
                                topology.atoms.as_ptr(),
                                topology.bonds.as_ptr(),
                                properties.props() as *const _,
                                &params as *const _,
                            )
                        && before_map_bits
                            == params
                                .coordinate_map
                                .iter()
                                .map(|(id, xy)| (*id, xy.map(f64::to_bits)))
                                .collect::<Vec<_>>();
                    if preservation {
                        preserved += 1;
                    } else {
                        discrepancies.push(format!("{label}/{repeat}: input preservation"));
                    }
                    let delta = DISPATCHES.with(Cell::get) - baseline;
                    let expected_delta = usize::from(properties.prop("_StereochemDone").is_none());
                    if delta == expected_delta {
                        dispatch_matches += 1;
                    } else {
                        discrepancies.push(format!("{label}/{repeat}: actual dispatch delta {delta} expected {expected_delta}"));
                    }
                    let working = CAPTURE.with(|capture| capture.borrow_mut().1.take());
                    if working.as_ref() == Some(&expected_working) {
                        working_matches += 1;
                    } else {
                        discrepancies.push(format!("{label}/{repeat}: complete working topology differs\nactual={working:?}\nexpected={expected_working:?}"));
                    }
                    match result {
                        Err(error) => {
                            errors += 1;
                            discrepancies.push(format!("{label}/{repeat}: {error:?}"));
                        }
                        Ok(conformer) => {
                            let mut matches = conformer.coordinates().len() == expected_xy.len();
                            for &(atom, expected) in &expected_xy {
                                if let Some(actual) = conformer.coordinates().get(atom) {
                                    for axis in 0..2 {
                                        scalar_bits += 1;
                                        equal_bits +=
                                            usize::from(actual[axis].to_bits() == expected[axis]);
                                        if (actual[axis] - f64::from_bits(expected[axis])).abs()
                                            > 1e-8
                                        {
                                            matches = false;
                                            discrepancies.push(format!("{label}/{repeat}: XY {atom}/{axis} actual={} expected={} bits={}/{}", actual[axis], f64::from_bits(expected[axis]), actual[axis].to_bits(), expected[axis]));
                                        }
                                    }
                                } else {
                                    matches = false;
                                }
                            }
                            if matches {
                                coordinate_matches += 1;
                            } else {
                                discrepancies.push(format!(
                                    "{label}/{repeat}: ordered coordinate rows differ"
                                ));
                            }
                        }
                    }
                }
            }
            CAPTURE.with(|capture| capture.borrow_mut().0 = false);
            println!(
                "D2-PREPARED domain_calls={calls} errors={errors} preservation={preserved}/96 dispatch={dispatch_matches}/96 working={working_matches}/96 coordinates={coordinate_matches}/96 equal_bits={equal_bits}/{scalar_bits}"
            );
            assert_eq!(calls, 96);
            assert!(discrepancies.is_empty(), "{}", discrepancies.join("\n"));
        }
    }
    #[derive(Default)]
    struct SourceInput {
        atoms: Vec<Vec<String>>,
        bonds: Vec<Vec<String>>,
    }

    fn parse<T>(field: &str) -> T
    where
        T: FromStr,
        T::Err: std::fmt::Debug,
    {
        field.parse().expect("valid pinned source trace field")
    }

    fn source_inputs(text: &str) -> Vec<SourceInput> {
        let mut inputs = Vec::new();
        let mut current = None;
        for line in text.lines() {
            let fields = line.split('\t').collect::<Vec<_>>();
            if fields.len() < 3 || fields[0] != "D2" {
                continue;
            }
            if fields[1] == "constructor_final" && fields[2] == "META" {
                if let Some(input) = current.take() {
                    inputs.push(input);
                }
                current = Some(SourceInput::default());
                continue;
            }
            if fields[1] != "constructor_final" {
                continue;
            }
            let Some(input) = current.as_mut() else {
                continue;
            };
            match fields[2] {
                "ATOM" => input.atoms.push(
                    fields[3..]
                        .iter()
                        .map(|field| (*field).to_owned())
                        .collect(),
                ),
                "BOND" => input.bonds.push(
                    fields[3..]
                        .iter()
                        .map(|field| (*field).to_owned())
                        .collect(),
                ),
                _ => {}
            }
        }
        if let Some(input) = current {
            inputs.push(input);
        }
        inputs
    }

    fn topology(input: &SourceInput) -> TopologyBlock {
        let mut atoms = Vec::with_capacity(input.atoms.len());
        for row in &input.atoms {
            assert_eq!(row.len(), 11);
            let index: usize = parse(&row[0]);
            assert_eq!(index, atoms.len(), "source atom table index is ordered");
            let element = Element::from_atomic_number(parse(&row[1]))
                .expect("pinned source atomic number is modeled");
            let hybridization = Hybridization::from_rdkit_code(parse(&row[8]))
                .expect("pinned source hybridization is modeled");
            let chiral_tag = ChiralTag::from_rdkit_code(parse(&row[9]))
                .expect("pinned source chiral tag is modeled");
            let mut spec = AtomSpec::new(element)
                .with_isotope(parse(&row[2]))
                .with_formal_charge(parse(&row[3]))
                .with_explicit_hydrogens(parse(&row[4]))
                .with_no_implicit(parse::<u8>(&row[5]) != 0)
                .with_radical_electrons(parse(&row[6]))
                .with_aromatic(parse::<u8>(&row[7]) != 0)
                .with_hybridization(hybridization)
                .with_chiral_tag(chiral_tag);
            let atom_map: u32 = parse(&row[10]);
            if atom_map != 0 {
                spec = spec.with_atom_map(atom_map);
            }
            atoms.push(Atom::from_spec(AtomId::new(index), spec));
        }

        let mut bonds = Vec::with_capacity(input.bonds.len());
        for row in &input.bonds {
            assert_eq!(row.len(), 9);
            let index: usize = parse(&row[0]);
            assert_eq!(index, bonds.len(), "source bond table index is ordered");
            let begin = AtomId::new(parse(&row[1]));
            let end = AtomId::new(parse(&row[2]));
            let order = BondOrder::from_rdkit_code(parse(&row[3]))
                .expect("pinned source bond order is modeled");
            let direction = BondDirection::from_rdkit_code(parse(&row[6]))
                .expect("pinned source bond direction is modeled");
            let stereo = BondStereo::from_rdkit_code(parse(&row[7]))
                .expect("pinned source bond stereo is modeled");
            let mut spec = BondSpec::new(begin, end, order)
                .with_aromatic(parse::<u8>(&row[4]) != 0)
                .with_conjugated(parse::<u8>(&row[5]) != 0)
                .with_direction(direction)
                .with_stereo(stereo);
            if !row[8].is_empty() {
                let stereo_atoms = row[8]
                    .split(',')
                    .map(|field| AtomId::new(parse(field)))
                    .collect::<Vec<_>>();
                assert_eq!(stereo_atoms.len(), 2, "source stereo references are paired");
                spec = spec.with_stereo_atoms(stereo_atoms[0], stereo_atoms[1]);
            }
            bonds.push(Bond::from_spec(BondId::new(index), spec));
        }
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("pinned constructor-final topology is modeled")
    }

    mod last_boundary {
        use super::*;
        pub(super) fn unhex(text: &str) -> String {
            assert_eq!(text.len() % 2, 0);
            String::from_utf8(
                text.as_bytes()
                    .chunks_exact(2)
                    .map(|pair| u8::from_str_radix(std::str::from_utf8(pair).unwrap(), 16).unwrap())
                    .collect(),
            )
            .expect("native property bytes are UTF-8 in these fixed cases")
        }
        pub(super) fn transport(case_trace: &str) -> (TopologyBlock, MoleculeProperties) {
            let inputs = source_inputs(case_trace);
            assert_eq!(inputs.len(), 1);
            let mut topology = topology(&inputs[0]);
            for line in case_trace.lines() {
                let fields = line.split('\t').collect::<Vec<_>>();
                if fields.len() < 3 || fields[..2] != ["D2", "constructor_final"] {
                    continue;
                }
                if fields[2] == "QUERY" {
                    assert_eq!(fields[5], "0", "selected constructor has no query");
                }
                if fields[2] != "PROP" || fields[3] == "MOL" {
                    continue;
                }
                assert_eq!(fields.len(), 10);
                let key = unhex(fields[6]);
                if key == "__computedProps" {
                    assert_eq!(fields[7], "12");
                    // Source vector metadata is retained in the raw input;
                    // actual scalar computed membership is transported below.
                    continue;
                }
                let value = unhex(fields[9]);
                let value = match fields[7] {
                    "1" => PropertyValue::Int(parse(&value)),
                    "6" => PropertyValue::UInt(parse::<u32>(&value)),
                    "2" => PropertyValue::Double(f64::from_bits(parse(&value))),
                    "3" => PropertyValue::String(value.into()),
                    "5" => PropertyValue::Bool(match value.as_str() {
                        "0" => false,
                        "1" => true,
                        _ => panic!("native boolean is 0/1"),
                    }),
                    other => panic!("unmodeled actual scalar property tag {other}"),
                };
                let index: usize = parse(fields[4]);
                let computed = fields[8] == "1";
                match fields[3] {
                    "ATOM" => {
                        if computed {
                            topology.atoms[index].set_computed_prop(key, value).unwrap();
                        } else {
                            topology.atoms[index].set_prop(key, value).unwrap();
                        }
                    }
                    "BOND" => {
                        if computed {
                            topology.bonds[index].set_computed_prop(key, value).unwrap();
                        } else {
                            topology.bonds[index].set_prop(key, value).unwrap();
                        }
                    }
                    other => panic!("unknown property owner {other}"),
                }
            }
            let mut properties = MoleculeProperties::default();
            for row in case_trace
                .lines()
                .filter(|row| row.starts_with("D2\tconstructor_final\tPROP\tMOL\t"))
            {
                let fields = row.split('\t').collect::<Vec<_>>();
                assert_eq!(fields.len(), 10);
                let key = unhex(fields[6]);
                if key == "__computedProps" {
                    assert_eq!(fields[7], "12");
                    // Native vector remains raw metadata, scalar membership
                    // is transported through the existing model markers.
                    continue;
                }
                assert!(matches!(fields[7], "1" | "3"));
                let value = unhex(fields[9]);
                if fields[8] == "1" {
                    properties.set_computed_prop(key, value).unwrap();
                } else {
                    properties.set_prop(key, value).unwrap();
                }
            }
            (topology, properties)
        }
    }
}
