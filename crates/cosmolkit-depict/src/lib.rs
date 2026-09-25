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
    Conformer2D, CoordinateBlock, CoordinateValidationError, TopologyBlock, TopologyValidationError,
};

use crate::embedded_frag::{
    EmbeddedFrag, FragmentError, embed_cis_trans_systems, embed_fused_systems,
    orient_and_shift_fragments, seed_coordinate_constraints,
    translate_single_coordinate_constraint,
};
use crate::geometry::{GeometryError, PointMap, atom_depict_rank};
use crate::nontetrahedral::embed_nontetrahedral_stereo;
use crate::templates::{CoordinateTemplates, TemplateError};

mod embedded_frag;
mod geometry;
mod nontetrahedral;
mod templates;

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
    AtomIndexOutOfRange { atom: usize, atom_count: usize },
    AtomAlreadyEmbedded { atom: usize },
    AtomNotEmbedded { atom: usize },
    NotEnoughEmbeddedNeighbors { atom: usize, count: usize },
    CoincidentPoints,
    InvalidAngle,
    NoCommonAtoms,
    MismatchedTopology,
    EmptyAttachment { atom: usize },
    NonTetrahedralNoLigand { centre: usize },
    NonTetrahedralLigandOverflow { centre: usize },
    CisTransBondInvalid { bond: usize },
    CollisionBondInvalid { bond: usize },
    UndefinedSamplingDistance { first: usize, second: usize },
    AtomCountTooLarge { atom_count: usize },
    InvalidRankProperty { atom: usize, key: &'static str },
    GeometryNotEnoughNeighbors { atom: usize, count: usize },
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
            Self::InvalidRankProperty { atom, key } => {
                write!(formatter, "atom {atom} has invalid rank property {key}")
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
            | Self::InvalidRankProperty { .. }
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
            GeometryError::InvalidRankProperty { atom, key } => {
                Self::InvalidRankProperty { atom, key }
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
    // RDKit❗✔️: void computeInitialCoords(RDKit::ROMol &mol,
    // RDKit❗✔️:                           const RDGeom::INT_POINT2D_MAP *coordMap,
    // RDKit❗✔️:                           std::list<EmbeddedFrag> &efrags,
    // RDKit❗✔️:                           bool useRingTemplates) {
    // RDKit❗✔️:   std::vector<int> atomRanks;
    // RDKit❗✔️:   atomRanks.resize(mol.getNumAtoms());
    // RDKit❗✔️:   for (auto i = 0u; i < mol.getNumAtoms(); ++i) {
    // RDKit❗✔️:     atomRanks[i] = getAtomDepictRank(mol.getAtomWithIdx(i));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   RDKit::VECT_INT_VECT arings;
    // RDKit❗✔️:   bool includeDativeBonds = true;
    // RDKit❗✔️:   RDKit::MolOps::symmetrizeSSSR(mol, arings, includeDativeBonds);
    // RDKit❗✔️:   RDKit::MolOps::assignStereochemistry(mol, false);
    // RDKit❗✔️:   efrags.clear();
    // RDKit❗✔️:   bool preSpec = false;
    // RDKit❗✔️:   if ((coordMap) && (coordMap->size() > 1)) {
    // RDKit❗✔️:     EmbeddedFrag efrag(&mol, *coordMap);
    // RDKit❗✔️:     efrags.push_back(efrag);
    // RDKit❗✔️:     preSpec = true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (arings.size() > 0) {
    // RDKit❗✔️:     DepictorLocal::embedFusedSystems(mol, arings, efrags, coordMap,
    // RDKit❗✔️:                                      useRingTemplates);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   DepictorLocal::embedNontetrahedralStereo(mol, efrags, atomRanks);
    // RDKit❗✔️:   DepictorLocal::embedCisTransSystems(mol, efrags);
    // RDKit❗✔️:   auto nratms = DepictorLocal::getNonEmbeddedAtoms(mol, efrags);
    // RDKit❗✔️:   std::list<EmbeddedFrag>::iterator mri;
    // RDKit❗✔️:   if (preSpec) { mri = efrags.begin(); }
    // RDKit❗✔️:   else { mri = DepictorLocal::_findLargestFrag(efrags); }
    // RDKit❗✔️:   while ((mri != efrags.end()) || (nratms.size() > 0)) {
    // RDKit❗✔️:     if (mri == efrags.end()) {
    // RDKit❗✔️:       auto mrank = static_cast<int>(RDKit::MAX_INT);
    // RDKit❗✔️:       RDKit::INT_LIST_I mnri;
    // RDKit❗✔️:       for (auto nri = nratms.begin(); nri != nratms.end(); ++nri) {
    // RDKit❗✔️:         auto rank = atomRanks.at(*nri);
    // RDKit❗✔️:         rank *= mol.getNumAtoms();
    // RDKit❗✔️:         rank += *nri;
    // RDKit❗✔️:         if (rank < mrank) { mrank = rank; mnri = nri; }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       EmbeddedFrag efrag((*mnri), &mol);
    // RDKit❗✔️:       nratms.erase(mnri);
    // RDKit❗✔️:       efrags.push_back(efrag);
    // RDKit❗✔️:       mri = efrags.end(); --mri;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     mri->markDone();
    // RDKit❗✔️:     mri->expandEfrag(nratms, efrags);
    // RDKit❗✔️:     mri = DepictorLocal::_findLargestFrag(efrags);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
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
            nonembedded.remove(position);
            fragments.push(EmbeddedFrag::from_single(atom, topology, rings)?);
            fragments.len() - 1
        };
        let mut fragment = fragments.remove(index);
        fragment.done = true;
        fragment.expand_fragment(&mut nonembedded, &mut fragments)?;
        fragments.insert(index.min(fragments.len()), fragment);
        selected = largest_unfinished_fragment(&fragments);
    }
    Ok(fragments)
}

/// Compute one detached atom-ordered 2D conformer.
pub fn compute_2d_coordinates(
    topology: &TopologyBlock,
    params: &Compute2DCoordinatesParams,
) -> Result<Conformer2D, DepictError> {
    // RDKit❗✔️: unsigned int compute2DCoords(RDKit::ROMol &mol,
    // RDKit❗✔️:                              const Compute2DCoordParameters &params) {
    // RDKit❗✔️:   if (mol.needsUpdatePropertyCache()) {
    // RDKit❗✔️:     mol.updatePropertyCache(false);
    // RDKit❗✔️:   }
    // RDKit❗✔️: #ifdef RDK_BUILD_COORDGEN_SUPPORT
    // RDKit❗✔️:   if (!params.forceRDKit && preferCoordGen) {
    // RDKit❗✔️:     RDKit::CoordGen::CoordGenParams coordgen_params;
    // RDKit❗✔️:     if (params.coordMap) { coordgen_params.coordMap = *params.coordMap; }
    // RDKit❗✔️:     auto cid = RDKit::CoordGen::addCoords(mol, &coordgen_params);
    // RDKit❗✔️:     return cid;
    // RDKit❗✔️:   };
    // RDKit❗✔️: #endif
    // RDKit❗✔️:   RDKit::ROMol cp(mol);
    // RDKit❗✔️:   std::list<EmbeddedFrag> efrags;
    // RDKit❗✔️:   computeInitialCoords(cp, params.coordMap, efrags, params.useRingTemplates);
    // RDKit❗✔️:   for (auto &eri : efrags) {
    // RDKit❗✔️:     if ((params.nSamples > 0) && (params.nFlipsPerSample > 0)) {
    // RDKit❗✔️:       eri.randomSampleFlipsAndPermutations(params.nFlipsPerSample,
    // RDKit❗✔️:           params.nSamples, params.sampleSeed, nullptr, 0.0,
    // RDKit❗✔️:           params.permuteDeg4Nodes);
    // RDKit❗✔️:     } else { eri.removeCollisionsBondFlip(); }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (auto &eri : efrags) {
    // RDKit❗✔️:     eri.removeCollisionsOpenAngles();
    // RDKit❗✔️:     eri.removeCollisionsShortenBonds();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!params.coordMap || !params.coordMap->size()) {
    // RDKit❗✔️:     if (params.canonOrient && efrags.size()) {
    // RDKit❗✔️:       for (auto &eri : efrags) { eri.canonicalizeOrientation(); }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   DepictorLocal::_shiftCoords(efrags);
    // RDKit❗✔️:   auto cid = copyCoordinate(mol, efrags, params.clearConfs);
    // RDKit❗✔️:   if ((params.coordMap) && (params.coordMap->size() == 1)) {
    // RDKit❗✔️:     // translate every copied conformer row to the singleton anchor
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return cid;
    // RDKit❗✔️: }
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
    let working = assign_legacy_stereochemistry_for_depiction(topology.clone(), &valence, &rings)
        .map_err(DepictError::StereoAssignment)?;
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
            fragment.remove_collisions_bond_flip()?;
        }
    }
    for fragment in &mut fragments {
        fragment.remove_collisions_open_angles()?;
        fragment.remove_collisions_shorten_bonds()?;
    }
    orient_and_shift_fragments(
        &mut fragments,
        params.canonical_orientation,
        Some(params.coordinate_map.len()),
    );
    translate_single_coordinate_constraint(&working, &mut fragments, coordinate_map)?;

    // RDKit❗✔️: auto *conf = new RDKit::Conformer(mol.getNumAtoms());
    // RDKit❗✔️: conf->set3D(false);
    // RDKit❗✔️: for (const auto &efrag : efrags) {
    // RDKit❗✔️:   for (const auto &eai : efrag.GetEmbeddedAtoms()) {
    // RDKit❗✔️:     const auto &cr = eai.second.loc;
    // RDKit❗✔️:     RDGeom::Point3D fcr(cr.x, cr.y, 0.0);
    // RDKit❗✔️:     conf->setAtomPos(eai.first, fcr);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
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
    let conformer = compute_2d_coordinates(topology, &Compute2DCoordinatesParams::default())?;
    Ok(CoordinateBlock {
        conformers_2d: vec![conformer],
        ..Default::default()
    })
}

pub fn render_svg(
    _topology: &TopologyBlock,
    _coordinates: &CoordinateBlock,
    _options: &DepictOptions,
) -> Result<String, DepictError> {
    Err(DepictError::CoordGenUnavailable)
}
