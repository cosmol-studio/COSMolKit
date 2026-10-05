//! Source-backed molecular alignment over explicit detached values.

mod support;

use std::borrow::Cow;
use thiserror::Error;

use cosmolkit_core::{RingInfo, ValenceAssignment, align_points, alignment_transform_point};
use cosmolkit_model::{Conformer3D, CoordinateBlock, CoordinateValidationError, TopologyBlock};
use cosmolkit_search::{SearchTarget, SubstructMatchParams};

type Transform3D = [[f64; 4]; 4];

/// Explicit detached input; no runtime capability or live molecule crosses this boundary.
#[derive(Clone, Copy)]
pub struct AlignmentInput<'a> {
    pub topology: &'a TopologyBlock,
    pub coordinates: &'a CoordinateBlock,
    pub rings: Option<&'a RingInfo>,
    pub valence: Option<&'a ValenceAssignment>,
}
impl AlignmentInput<'_> {
    fn conformers_3d(&self) -> &[Conformer3D] {
        &self.coordinates.conformers_3d
    }
    fn atoms(&self) -> &[cosmolkit_model::Atom] {
        &self.topology.atoms
    }
    fn search_target(&self) -> SearchTarget<'_> {
        SearchTarget::new(
            self.topology,
            self.coordinates,
            &self.topology.stereo_groups,
            self.rings,
            self.valence,
        )
    }
}

const DEFAULT_MAX_ITERATIONS: u32 = 50;
const DEFAULT_MAX_MATCHES: i32 = 1_000_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentAtomMap {
    pub probe_atom: usize,
    pub reference_atom: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlignmentParameters {
    pub probe_conformer_id: i32,
    pub reference_conformer_id: i32,
    pub atom_map: Option<Vec<AlignmentAtomMap>>,
    pub weights: Option<Vec<f64>>,
    pub reflect: bool,
    pub max_iterations: u32,
}

impl Default for AlignmentParameters {
    fn default() -> Self {
        Self {
            probe_conformer_id: -1,
            reference_conformer_id: -1,
            atom_map: None,
            weights: None,
            reflect: false,
            max_iterations: DEFAULT_MAX_ITERATIONS,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BestAlignmentParameters {
    pub probe_conformer_id: i32,
    pub reference_conformer_id: i32,
    pub atom_maps: Vec<Vec<AlignmentAtomMap>>,
    pub weights: Option<Vec<f64>>,
    pub reflect: bool,
    pub max_iterations: u32,
    pub max_matches: i32,
    pub symmetrize_conjugated_terminal_groups: bool,
    pub ignore_hydrogens: bool,
    pub num_threads: i32,
}

impl Default for BestAlignmentParameters {
    fn default() -> Self {
        Self {
            probe_conformer_id: -1,
            reference_conformer_id: -1,
            atom_maps: Vec::new(),
            weights: None,
            reflect: false,
            max_iterations: DEFAULT_MAX_ITERATIONS,
            max_matches: DEFAULT_MAX_MATCHES,
            symmetrize_conjugated_terminal_groups: true,
            ignore_hydrogens: true,
            num_threads: 1,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AllConformerRmsdParameters {
    pub atom_maps: Vec<Vec<AlignmentAtomMap>>,
    pub weights: Option<Vec<f64>>,
    pub max_matches: i32,
    pub symmetrize_conjugated_terminal_groups: bool,
    pub ignore_hydrogens: bool,
    pub num_threads: i32,
}

impl Default for AllConformerRmsdParameters {
    fn default() -> Self {
        Self {
            atom_maps: Vec::new(),
            weights: None,
            max_matches: DEFAULT_MAX_MATCHES,
            symmetrize_conjugated_terminal_groups: true,
            ignore_hydrogens: true,
            num_threads: 1,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CoordinateRmsdParameters {
    pub probe_conformer_id: i32,
    pub reference_conformer_id: i32,
    pub atom_maps: Vec<Vec<AlignmentAtomMap>>,
    pub weights: Option<Vec<f64>>,
    pub max_matches: i32,
    pub symmetrize_conjugated_terminal_groups: bool,
}

impl Default for CoordinateRmsdParameters {
    fn default() -> Self {
        Self {
            probe_conformer_id: -1,
            reference_conformer_id: -1,
            atom_maps: Vec::new(),
            weights: None,
            max_matches: DEFAULT_MAX_MATCHES,
            symmetrize_conjugated_terminal_groups: true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlignmentTransform {
    pub matrix: Transform3D,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlignmentResult {
    pub rmsd: f64,
    pub transform: AlignmentTransform,
    pub atom_map: Vec<AlignmentAtomMap>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConformerRmsd {
    pub probe_conformer_id: usize,
    pub reference_conformer_id: usize,
    pub rmsd: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ConformerAlignmentReport {
    pub rmsds: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq, Error)]
pub enum AlignmentError {
    #[error(transparent)]
    InvalidCoordinates(#[from] CoordinateValidationError),
    #[error(transparent)]
    Matching(#[from] cosmolkit_search::SubstructMatchError),
    #[error(transparent)]
    QueryGraph(#[from] cosmolkit_model::QueryGraphError),
    #[error(transparent)]
    QueryParse(#[from] cosmolkit_search::SmartsParseError),
    #[error(transparent)]
    ThreadSelection(#[from] ThreadSelectionError),
    #[error("molecule has no 3D conformers")]
    NoConformers,
    #[error("conformer id {id} was not found")]
    ConformerNotFound { id: i32 },
    #[error("alignment atom map is empty")]
    EmptyAtomMap,
    #[error("probe atom index {index} is out of range for {atom_count} atoms")]
    ProbeAtomOutOfRange { index: usize, atom_count: usize },
    #[error("reference atom index {index} is out of range for {atom_count} atoms")]
    ReferenceAtomOutOfRange { index: usize, atom_count: usize },
    #[error("alignment has {map_len} mapped atoms but {weight_len} weights")]
    WeightCountMismatch { map_len: usize, weight_len: usize },
    #[error("alignment weight at index {index} must be positive")]
    NonPositiveWeight { index: usize },
    #[error("no substructure match found between the reference and probe molecules")]
    NoSubstructureMatch,
    #[error("terminal-group symmetrization failed: {message}")]
    TerminalGroupSymmetrization { message: &'static str },
    #[error("alignment numerical precondition failed: {message}")]
    NumericalPrecondition { message: &'static str },
    #[error("alignment worker terminated unexpectedly")]
    WorkerTerminated,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ConformerAlignmentParameters {
    pub atom_indices: Option<Vec<usize>>,
    pub conformer_ids: Option<Vec<usize>>,
    pub weights: Option<Vec<f64>>,
    pub reflect: bool,
    pub max_iterations: u32,
}

impl Default for ConformerAlignmentParameters {
    fn default() -> Self {
        Self {
            atom_indices: None,
            conformer_ids: None,
            weights: None,
            reflect: false,
            max_iterations: DEFAULT_MAX_ITERATIONS,
        }
    }
}

#[derive(Clone, Debug)]
pub struct ThreadSelectionError(std::sync::Arc<cosmolkit_core::ThreadCountError>);
impl PartialEq for ThreadSelectionError {
    fn eq(&self, other: &Self) -> bool {
        std::sync::Arc::ptr_eq(&self.0, &other.0)
    }
}
impl std::fmt::Display for ThreadSelectionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}
impl std::error::Error for ThreadSelectionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(self.0.as_ref())
    }
}

fn conformer_by_id<'a>(
    molecule: &'a AlignmentInput<'_>,
    id: i32,
) -> Result<&'a Conformer3D, AlignmentError> {
    // Verbatim pinned RDKit Code/GraphMol/ROMol.cpp, ROMol::getConformer.
    // RDKit✔️❌: const Conformer &ROMol::getConformer(int id) const {
    // RDKit✔️❌:   // make sure we have more than one conformation
    // RDKit✔️❌:   if (d_confs.size() == 0) {
    // RDKit✔️❌:     throw ConformerException("No conformations available on the molecule");
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   if (id < 0) {
    // RDKit✔️❌:     return *(d_confs.front());
    // RDKit✔️❌:   }
    // RDKit✔️❌:   auto cid = (unsigned int)id;
    // RDKit✔️❌:   for (auto conf : d_confs) {
    // RDKit✔️❌:     if (conf->getId() == cid) {
    // RDKit✔️❌:       return *conf;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   // we did not find a conformation with the specified ID
    // RDKit✔️❌:   std::string mesg = "Can't find conformation with ID: ";
    // RDKit✔️❌:   mesg += id;
    // RDKit✔️❌:   throw ConformerException(mesg);
    // RDKit✔️❌: }
    // Source selection is O(C). The detached public boundary additionally
    // validates the selected f64 rows in O(N); this cost is recorded, not hidden.

    let conformers = molecule.conformers_3d();
    let index = if id < 0 {
        (!conformers.is_empty()).then_some(0)
    } else {
        conformers.iter().position(|c| c.id() == id as usize)
    }
    .ok_or_else(|| {
        if conformers.is_empty() {
            AlignmentError::NoConformers
        } else {
            AlignmentError::ConformerNotFound { id }
        }
    })?;
    conformers[index].validate_for_atom_count(molecule.topology.atoms.len())?;
    Ok(&conformers[index])
}

fn first_alignment_map(
    probe: &AlignmentInput<'_>,
    reference: &AlignmentInput<'_>,
    params: &AlignmentParameters,
) -> Result<Vec<AlignmentAtomMap>, AlignmentError> {
    // Performance loss: source borrows the ROMol for ordinary matching; this
    // detached query adapter clones O(V+E) atoms/bonds/stereo groups.

    if let Some(map) = &params.atom_map {
        return Ok(map.clone());
    }
    // RDKit✔️❌: double getAlignmentTransform(const ROMol &prbMol, const ROMol &refMol,
    // RDKit✔️❌:                              RDGeom::Transform3D &trans, int prbCid, int refCid,
    // RDKit✔️❌:                              const MatchVectType *atomMap,
    // RDKit✔️❌:                              const RDNumeric::DoubleVector *weights,
    // RDKit✔️❌:                              bool reflect, unsigned int maxIterations) {
    // RDKit✔️❌:   const Conformer &prbCnf = prbMol.getConformer(prbCid);
    // RDKit✔️❌:   const Conformer &refCnf = refMol.getConformer(refCid);
    // RDKit✔️❌:   MatchVectType match;
    // RDKit✔️❌:   if (!atomMap) {
    // RDKit✔️❌:     // we have to figure out the mapping between the two molecule
    // RDKit✔️❌:     const bool recursionPossible = true;
    // RDKit✔️❌:     const bool useChirality = false;
    // RDKit✔️❌:     const bool useQueryQueryMatches = true;
    // RDKit✔️❌:     if (SubstructMatch(refMol, prbMol, match, recursionPossible, useChirality,
    // RDKit✔️❌:                        useQueryQueryMatches)) {
    // RDKit✔️❌:       atomMap = &match;
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       throw MolAlignException(
    // RDKit✔️❌:           "No sub-structure match found between the probe and query mol");
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   double msd = alignConfsOnAtomMap(prbCnf, refCnf, *atomMap, trans, weights,
    // RDKit✔️❌:                                    reflect, maxIterations);
    // RDKit✔️❌:   return sqrt(msd);
    // RDKit✔️❌: }
    let matcher_params = SubstructMatchParams {
        max_matches: 1,
        use_chirality: false,
        use_query_query_matches: true,
        ..SubstructMatchParams::default()
    };
    let matched = cosmolkit_search::try_get_substruct_matches_with_params(
        &reference.search_target(),
        &support::query_for(probe)?,
        &matcher_params,
    )?
    .into_iter()
    .next()
    .ok_or(AlignmentError::NoSubstructureMatch)?;
    Ok(matched
        .atom_mapping
        .into_iter()
        .enumerate()
        .map(|(probe_atom, reference_atom)| AlignmentAtomMap {
            probe_atom,
            reference_atom,
        })
        .collect())
}

fn validate_map(
    map: &[AlignmentAtomMap],
    probe: &Conformer3D,
    reference: &Conformer3D,
    weights: Option<&[f64]>,
    positive_weights: bool,
) -> Result<(), AlignmentError> {
    if map.is_empty() {
        return Err(AlignmentError::EmptyAtomMap);
    }
    for entry in map {
        if entry.probe_atom >= probe.coordinates().len() {
            return Err(AlignmentError::ProbeAtomOutOfRange {
                index: entry.probe_atom,
                atom_count: probe.coordinates().len(),
            });
        }
        if entry.reference_atom >= reference.coordinates().len() {
            return Err(AlignmentError::ReferenceAtomOutOfRange {
                index: entry.reference_atom,
                atom_count: reference.coordinates().len(),
            });
        }
    }
    if let Some(weights) = weights {
        if weights.len() != map.len() {
            return Err(AlignmentError::WeightCountMismatch {
                map_len: map.len(),
                weight_len: weights.len(),
            });
        }
        if positive_weights && let Some(index) = weights.iter().position(|weight| !(*weight > 0.0))
        {
            return Err(AlignmentError::NonPositiveWeight { index });
        }
    }
    Ok(())
}

fn automatic_maps(
    probe: &AlignmentInput<'_>,
    reference: &AlignmentInput<'_>,
    max_matches: i32,
    symmetrize_conjugated_terminal_groups: bool,
    ignore_hydrogens: bool,
) -> Result<Vec<Vec<AlignmentAtomMap>>, AlignmentError> {
    // Performance loss: source borrows the ROMol for ordinary matching; this
    // detached query adapter clones O(V+E) atoms/bonds/stereo groups.

    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8: Code/GraphMol/MolAlign/AlignMolecules.cpp
    // BEGIN VERBATIM CPP getAllMatchesPrbRef
    // RDKit✔️❌: void getAllMatchesPrbRef(const ROMol &prbMol, const ROMol &refMol,
    // RDKit✔️❌:                          std::vector<MatchVectType> &matches, int maxMatches,
    // RDKit✔️❌:                          bool symmetrizeConjugatedTerminalGroups,
    // RDKit✔️❌:                          bool ignoreHs = false) {
    // RDKit✔️❌:   bool uniquify = false;
    // RDKit✔️❌:   bool recursionPossible = true;
    // RDKit✔️❌:   bool useChirality = false;
    // RDKit✔️❌:   bool useQueryQueryMatches = false;
    // RDKit✔️❌:
    // RDKit✔️❌:   std::unique_ptr<RWMol> prbMolSymm;
    // RDKit✔️❌:   if (symmetrizeConjugatedTerminalGroups) {
    // RDKit✔️❌:     prbMolSymm.reset(new RWMol(prbMol));
    // RDKit✔️❌:     details::symmetrizeTerminalAtoms(*prbMolSymm);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   const auto &prbMolForMatch = prbMolSymm ? *prbMolSymm : prbMol;
    // RDKit✔️❌:   SubstructMatch(refMol, prbMolForMatch, matches, uniquify, recursionPossible,
    // RDKit✔️❌:                  useChirality, useQueryQueryMatches, maxMatches);
    // RDKit✔️❌:
    // RDKit✔️❌:   if (matches.empty()) {
    // RDKit✔️❌:     throw MolAlignException(
    // RDKit✔️❌:         "No sub-structure match found between the reference and probe mol");
    // RDKit✔️❌:   }
    // Independent rdWarningLog stream/configuration is not represented in
    // this detached API; enumeration, errors and filtering remain supported.
    // RDKit❌❌:   if (matches.size() > 1e6) {
    // RDKit❌❌:     std::string name;
    // RDKit❌❌:     prbMol.getPropIfPresent(common_properties::_Name, name);
    // RDKit❌❌:     BOOST_LOG(rdWarningLog)
    // RDKit❌❌:         << "Warning in " << __FUNCTION__ << ": " << matches.size()
    // RDKit❌❌:         << " matches detected for molecule " << name << ", this may "
    // RDKit❌❌:         << "lead to a performance slowdown." << std::endl;
    // RDKit❌❌:   }
    // RDKit✔️❌:   if (ignoreHs) {
    // RDKit✔️❌:     // filter Hs out of the matches
    // RDKit✔️❌:     for (auto &match : matches) {
    // RDKit✔️❌:       std::erase_if(match, [&prbMolForMatch](const auto &mi) {
    // RDKit✔️❌:         return prbMolForMatch.getAtomWithIdx(mi.first)->getAtomicNum() == 1;
    // RDKit✔️❌:       });
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END VERBATIM CPP getAllMatchesPrbRef

    let symmetrized_probe;
    let probe_for_match = if symmetrize_conjugated_terminal_groups {
        symmetrized_probe = support::symmetrize_terminal_atoms(probe)?;
        &symmetrized_probe
    } else {
        &support::query_for(probe)?
    };
    let matcher_params = SubstructMatchParams {
        uniquify: false,
        // RDKit's legacy overload accepts unsigned int; BestAlignmentParams
        // stores int and relies on the same modulo-2^32 conversion here.
        max_matches: max_matches as u32 as usize,
        use_chirality: false,
        ..Default::default()
    };
    let maps: Vec<_> = cosmolkit_search::try_get_substruct_matches_with_params(
        &reference.search_target(),
        probe_for_match,
        &matcher_params,
    )?
    .into_iter()
    .map(|matched| {
        matched
            .atom_mapping
            .into_iter()
            .enumerate()
            .filter(|(probe_atom, _)| {
                !ignore_hydrogens || probe.atoms()[*probe_atom].atomic_number() != 1
            })
            .map(|(probe_atom, reference_atom)| AlignmentAtomMap {
                probe_atom,
                reference_atom,
            })
            .collect::<Vec<_>>()
    })
    .collect();
    if maps.is_empty() {
        Err(AlignmentError::NoSubstructureMatch)
    } else {
        Ok(maps)
    }
}

fn maps_for<'a>(
    probe: &AlignmentInput<'_>,
    reference: &AlignmentInput<'_>,
    atom_maps: &'a [Vec<AlignmentAtomMap>],
    max_matches: i32,
    symmetrize_conjugated_terminal_groups: bool,
    ignore_hydrogens: bool,
) -> Result<Cow<'a, [Vec<AlignmentAtomMap>]>, AlignmentError> {
    // Performance: automatic ordinary matching allocates an O(V+E) query clone;
    // explicit maps remain borrowed without that adapter allocation.
    if atom_maps.is_empty() {
        automatic_maps(
            probe,
            reference,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            ignore_hydrogens,
        )
        .map(Cow::Owned)
    } else {
        Ok(Cow::Borrowed(atom_maps))
    }
}

fn aligned_result_with_msd(
    probe: &Conformer3D,
    reference: &Conformer3D,
    map: &[AlignmentAtomMap],
    weights: Option<&[f64]>,
    reflect: bool,
    max_iterations: u32,
) -> Result<(f64, AlignmentResult), AlignmentError> {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8: Code/GraphMol/MolAlign/AlignMolecules.cpp
    // BEGIN VERBATIM CPP alignConfsOnAtomMap
    // RDKit✔️✔️: double alignConfsOnAtomMap(const Conformer &prbCnf, const Conformer &refCnf,
    // RDKit✔️✔️:                            const MatchVectType &atomMap,
    // RDKit✔️✔️:                            RDGeom::Transform3D &trans,
    // RDKit✔️✔️:                            const RDNumeric::DoubleVector *weights, bool reflect,
    // RDKit✔️✔️:                            unsigned int maxIterations) {
    // RDKit✔️✔️:   RDGeom::Point3DConstPtrVect refPoints, prbPoints;
    // RDKit✔️✔️:   for (const auto &mi : atomMap) {
    // RDKit✔️✔️:     prbPoints.push_back(&prbCnf.getAtomPos(mi.first));
    // RDKit✔️✔️:     refPoints.push_back(&refCnf.getAtomPos(mi.second));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double ssr = RDNumeric::Alignments::AlignPoints(
    // RDKit✔️✔️:       refPoints, prbPoints, trans, weights, reflect, maxIterations);
    // RDKit✔️✔️:   return ssr / static_cast<double>(prbPoints.size());
    // RDKit✔️✔️: }
    // END VERBATIM CPP alignConfsOnAtomMap

    validate_map(map, probe, reference, weights, true)?;
    let probe_points: Vec<_> = map
        .iter()
        .map(|entry| probe.coordinates()[entry.probe_atom])
        .collect();
    let reference_points: Vec<_> = map
        .iter()
        .map(|entry| reference.coordinates()[entry.reference_atom])
        .collect();
    let (ssr, matrix) = align_points(
        &reference_points,
        &probe_points,
        weights,
        reflect,
        max_iterations as usize,
    )
    .map_err(|message| AlignmentError::NumericalPrecondition { message })?;
    Ok((
        ssr / map.len() as f64,
        AlignmentResult {
            rmsd: (ssr / map.len() as f64).sqrt(),
            transform: AlignmentTransform { matrix },
            atom_map: Vec::new(),
        },
    ))
}

fn aligned_result(
    probe: &Conformer3D,
    reference: &Conformer3D,
    map: &[AlignmentAtomMap],
    weights: Option<&[f64]>,
    reflect: bool,
    max_iterations: u32,
) -> Result<AlignmentResult, AlignmentError> {
    let mut result =
        aligned_result_with_msd(probe, reference, map, weights, reflect, max_iterations)?.1;
    result.atom_map = map.to_vec();
    Ok(result)
}

fn num_threads_to_use(target: i32) -> Result<usize, AlignmentError> {
    cosmolkit_core::rdkit_thread_count(target)
        .map(|count| count.get() as usize)
        .map_err(|cause| ThreadSelectionError(std::sync::Arc::new(cause)).into())
}

fn best_aligned_result(
    probe: &Conformer3D,
    reference: &Conformer3D,
    maps: &[Vec<AlignmentAtomMap>],
    weights: Option<&[f64]>,
    reflect: bool,
    max_iterations: u32,
    num_threads: i32,
) -> Result<AlignmentResult, AlignmentError> {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8: Code/GraphMol/MolAlign/AlignMolecules.cpp
    // BEGIN VERBATIM CPP getBestRMSInternal
    // RDKit✔️✔️: double getBestRMSInternal(const ROMol &prbMol, const ROMol &refMol, int prbCid,
    // RDKit✔️✔️:                           int refCid, const std::vector<MatchVectType> &matches,
    // RDKit✔️✔️:                           RDGeom::Transform3D *trans, MatchVectType *bestMatch,
    // RDKit✔️✔️:                           const RDNumeric::DoubleVector *weights, bool reflect,
    // RDKit✔️✔️:                           unsigned int maxIters, unsigned int numThreads) {
    // RDKit✔️✔️:   PRECONDITION(!matches.empty(), "matches must not be empty");
    // RDKit✔️✔️: #ifndef RDK_BUILD_THREADSAFE_SSS
    // RDKit✔️✔️:   numThreads = 1;
    // RDKit✔️✔️: #endif
    // RDKit✔️✔️:   double msdBest = std::numeric_limits<double>::max();
    // RDKit✔️✔️:   const Conformer &prbCnf = prbMol.getConformer(prbCid);
    // RDKit✔️✔️:   const Conformer &refCnf = refMol.getConformer(refCid);
    // RDKit✔️✔️:   const MatchVectType *bestMatchPtr = &matches[0];
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (numThreads == 1) {
    // RDKit✔️✔️:     for (const auto &matche : matches) {
    // RDKit✔️✔️:       RDGeom::Transform3D tmpTrans;
    // RDKit✔️✔️:       double msd = trans ? alignConfsOnAtomMap(prbCnf, refCnf, matche, tmpTrans,
    // RDKit✔️✔️:                                                weights, reflect, maxIters)
    // RDKit✔️✔️:                          : calcMSDInternal(prbCnf, refCnf, matche, weights);
    // RDKit✔️✔️:       if (msd < msdBest) {
    // RDKit✔️✔️:         msdBest = msd;
    // RDKit✔️✔️:         bestMatchPtr = &matche;
    // RDKit✔️✔️:         if (trans) {
    // RDKit✔️✔️:           trans->assign(tmpTrans);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit✔️✔️:   else {
    // RDKit✔️✔️:     std::vector<std::thread> tg;
    // RDKit✔️✔️:     std::vector<
    // RDKit✔️✔️:         std::vector<std::tuple<double, unsigned int, RDGeom::Transform3D>>>
    // RDKit✔️✔️:         rmsds(numThreads);
    // RDKit✔️✔️:     for (auto ti = 0u; ti < numThreads; ++ti) {
    // RDKit✔️✔️:       auto func = [&](unsigned int tidx) {
    // RDKit✔️✔️:         for (auto midx = tidx; midx < matches.size(); midx += numThreads) {
    // RDKit✔️✔️:           auto matche = matches[midx];
    // RDKit✔️✔️:           RDGeom::Transform3D tmpTrans;
    // RDKit✔️✔️:           auto msd = trans
    // RDKit✔️✔️:                          ? alignConfsOnAtomMap(prbCnf, refCnf, matche, tmpTrans,
    // RDKit✔️✔️:                                                weights, reflect, maxIters)
    // RDKit✔️✔️:                          : calcMSDInternal(prbCnf, refCnf, matche, weights);
    // RDKit✔️✔️:           rmsds[tidx].emplace_back(msd, midx, tmpTrans);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       };
    // RDKit✔️✔️:       tg.emplace_back(std::thread(func, ti));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (auto &thread : tg) {
    // RDKit✔️✔️:       if (thread.joinable()) {
    // RDKit✔️✔️:         thread.join();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (const auto &rv : rmsds) {
    // RDKit✔️✔️:       for (const auto &res : rv) {
    // RDKit✔️✔️:         const auto &[msd, midx, tf] = res;
    // RDKit✔️✔️:         if (msd < msdBest) {
    // RDKit✔️✔️:           msdBest = msd;
    // RDKit✔️✔️:           bestMatchPtr = &matches[midx];
    // RDKit✔️✔️:           if (trans) {
    // RDKit✔️✔️:             trans->assign(tf);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: #endif
    // RDKit✔️✔️:   if (bestMatch) {
    // RDKit✔️✔️:     *bestMatch = *bestMatchPtr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return sqrt(msdBest);
    // RDKit✔️✔️: }
    // END VERBATIM CPP getBestRMSInternal

    if maps.is_empty() {
        return Err(AlignmentError::NoSubstructureMatch);
    }
    for map in maps {
        validate_map(map, probe, reference, weights, true)?;
    }
    let thread_count = num_threads_to_use(num_threads)?;
    if thread_count == 1 {
        let mut best_msd = f64::MAX;
        let mut best = AlignmentResult {
            rmsd: best_msd.sqrt(),
            transform: AlignmentTransform {
                matrix: std::array::from_fn(|r| {
                    std::array::from_fn(|c| {
                        cosmolkit_core::Transform3D::identity().values()[4 * r + c]
                    })
                }),
            },
            atom_map: maps[0].clone(),
        };
        for map in maps {
            let (msd, mut candidate) =
                aligned_result_with_msd(probe, reference, map, weights, reflect, max_iterations)?;
            if msd < best_msd {
                best_msd = msd;
                candidate.atom_map = map.to_vec();
                best = candidate;
            }
        }
        return Ok(best);
    }

    let buckets = std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(thread_count);
        for thread_index in 0..thread_count {
            handles.push(scope.spawn(move || {
                let mut bucket = Vec::new();
                let mut map_index = thread_index;
                while map_index < maps.len() {
                    bucket.push((
                        map_index,
                        aligned_result_with_msd(
                            probe,
                            reference,
                            &maps[map_index],
                            weights,
                            reflect,
                            max_iterations,
                        ),
                    ));
                    map_index += thread_count;
                }
                bucket
            }));
        }
        handles
            .into_iter()
            .map(|handle| handle.join())
            .collect::<Vec<_>>()
    });
    let mut best_msd = f64::MAX;
    let mut best = AlignmentResult {
        rmsd: best_msd.sqrt(),
        transform: AlignmentTransform {
            matrix: std::array::from_fn(|r| {
                std::array::from_fn(|c| cosmolkit_core::Transform3D::identity().values()[4 * r + c])
            }),
        },
        atom_map: maps[0].clone(),
    };
    for bucket in buckets {
        let bucket = bucket.map_err(|_| AlignmentError::WorkerTerminated)?;
        for (map_index, candidate) in bucket {
            let (msd, mut candidate) = candidate?;
            if msd < best_msd {
                best_msd = msd;
                candidate.atom_map = maps[map_index].clone();
                best = candidate;
            }
        }
    }
    Ok(best)
}

fn coordinate_rmsd(
    probe: &Conformer3D,
    reference: &Conformer3D,
    map: &[AlignmentAtomMap],
    weights: Option<&[f64]>,
) -> Result<f64, AlignmentError> {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8: Code/GraphMol/MolAlign/AlignMolecules.cpp
    // BEGIN VERBATIM CPP calcMSDInternal
    // RDKit✔️✔️: double calcMSDInternal(const Conformer &prbCnf, const Conformer &refCnf,
    // RDKit✔️✔️:                        const MatchVectType &atomMap,
    // RDKit✔️✔️:                        const RDNumeric::DoubleVector *weights) {
    // RDKit✔️✔️:   unsigned int npt = atomMap.size();
    // RDKit✔️✔️:   std::unique_ptr<RDNumeric::DoubleVector> unitWeights;
    // RDKit✔️✔️:   if (!weights) {
    // RDKit✔️✔️:     unitWeights.reset(new RDNumeric::DoubleVector(npt, 1.0));
    // RDKit✔️✔️:     weights = unitWeights.get();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     PRECONDITION(npt == weights->size(), "Mismatch in number of weights");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   RDGeom::Point3DConstPtrVect refPoints, prbPoints;
    // RDKit✔️✔️:   for (const auto &mi : atomMap) {
    // RDKit✔️✔️:     prbPoints.push_back(&prbCnf.getAtomPos(mi.first));
    // RDKit✔️✔️:     refPoints.push_back(&refCnf.getAtomPos(mi.second));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double ssr = 0.;
    // RDKit✔️✔️:   const RDGeom::Point3D *rpt;
    // RDKit✔️✔️:   const RDGeom::Point3D *ppt;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < npt; ++i) {
    // RDKit✔️✔️:     rpt = refPoints[i];
    // RDKit✔️✔️:     ppt = prbPoints[i];
    // RDKit✔️✔️:     ssr += (*weights)[i] * (*ppt - *rpt).lengthSq();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return ssr / static_cast<double>(npt);
    // RDKit✔️✔️: }
    // END VERBATIM CPP calcMSDInternal

    validate_map(map, probe, reference, weights, false)?;
    let mut sum = 0.0;
    for (index, entry) in map.iter().enumerate() {
        let p = probe.coordinates()[entry.probe_atom];
        let r = reference.coordinates()[entry.reference_atom];
        let weight = weights.map_or(1.0, |weights| weights[index]);
        sum += weight * ((p[0] - r[0]).powi(2) + (p[1] - r[1]).powi(2) + (p[2] - r[2]).powi(2));
    }
    Ok(sum / map.len() as f64)
}

pub fn align_conformers(
    coordinates: &mut CoordinateBlock,
    atom_count: usize,
    params: &ConformerAlignmentParameters,
) -> Result<Vec<f64>, AlignmentError> {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8: Code/GraphMol/MolAlign/AlignMolecules.cpp
    // BEGIN VERBATIM CPP alignMolConformers
    // RDKit✔️✔️: void alignMolConformers(ROMol &mol, const std::vector<unsigned int> *atomIds,
    // RDKit✔️✔️:                         const std::vector<unsigned int> *confIds,
    // RDKit✔️✔️:                         const RDNumeric::DoubleVector *weights, bool reflect,
    // RDKit✔️✔️:                         unsigned int maxIters, std::vector<double> *RMSlist) {
    // RDKit✔️✔️:   if (mol.getNumConformers() == 0) {
    // RDKit✔️✔️:     // nothing to be done ;
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3DConstPtrVect refPoints, prbPoints;
    // RDKit✔️✔️:   int cid = -1;
    // RDKit✔️✔️:   if ((confIds != nullptr) && (confIds->size() > 0)) {
    // RDKit✔️✔️:     cid = confIds->front();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const Conformer &refCnf = mol.getConformer(cid);
    // RDKit✔️✔️:   _fillAtomPositions(refPoints, refCnf, atomIds);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // now loop throught the remaininf conformations and transform them
    // RDKit✔️✔️:   RDGeom::Transform3D trans;
    // RDKit✔️✔️:   double ssd;
    // RDKit✔️✔️:   if (confIds == nullptr) {
    // RDKit✔️✔️:     unsigned int i = 0;
    // RDKit✔️✔️:     ROMol::ConformerIterator cnfi;
    // RDKit✔️✔️:     // Conformer *conf;
    // RDKit✔️✔️:     for (cnfi = mol.beginConformers(); cnfi != mol.endConformers(); cnfi++) {
    // RDKit✔️✔️:       // conf = (*cnfi);
    // RDKit✔️✔️:       i += 1;
    // RDKit✔️✔️:       if (i == 1) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       _fillAtomPositions(prbPoints, *(*cnfi), atomIds);
    // RDKit✔️✔️:       ssd = RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans,
    // RDKit✔️✔️:                                                weights, reflect, maxIters);
    // RDKit✔️✔️:       if (RMSlist) {
    // RDKit✔️✔️:         ssd /= (prbPoints.size());
    // RDKit✔️✔️:         RMSlist->push_back(sqrt(ssd));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       MolTransforms::transformConformer(*(*cnfi), trans);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     std::vector<unsigned int>::const_iterator cai;
    // RDKit✔️✔️:     unsigned int i = 0;
    // RDKit✔️✔️:     for (cai = confIds->begin(); cai != confIds->end(); cai++) {
    // RDKit✔️✔️:       i += 1;
    // RDKit✔️✔️:       if (i == 1) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       Conformer &conf = mol.getConformer(*cai);
    // RDKit✔️✔️:       _fillAtomPositions(prbPoints, conf, atomIds);
    // RDKit✔️✔️:       ssd = RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans,
    // RDKit✔️✔️:                                                weights, reflect, maxIters);
    // RDKit✔️✔️:       if (RMSlist) {
    // RDKit✔️✔️:         ssd /= (prbPoints.size());
    // RDKit✔️✔️:         RMSlist->push_back(sqrt(ssd));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       MolTransforms::transformConformer(conf, trans);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END VERBATIM CPP alignMolConformers

    coordinates.validate_for_atom_count(atom_count)?;
    if coordinates.conformers_3d.is_empty() {
        return Ok(Vec::new());
    }
    // Verbatim Python-wrapper normalization and detached position selection.
    // RDKit✔️✔️: void alignMolConfs(ROMol &mol, python::object atomIds, python::object confIds,
    // RDKit✔️✔️:                    python::object weights, bool reflect, unsigned int maxIters,
    // RDKit✔️✔️:                    python::object RMSlist) {
    // RDKit✔️✔️:   std::unique_ptr<RDNumeric::DoubleVector> wtsVec(translateDoubleSeq(weights));
    // RDKit✔️✔️:   std::unique_ptr<std::vector<unsigned int>> aIds(translateIntSeq(atomIds));
    // RDKit✔️✔️:   std::unique_ptr<std::vector<unsigned int>> cIds(translateIntSeq(confIds));
    // RDKit✔️✔️:   std::unique_ptr<std::vector<double>> RMSvector;
    // RDKit✔️✔️:   if (RMSlist != python::object()) {
    // RDKit✔️✔️:     RMSvector.reset(new std::vector<double>());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   {
    // RDKit✔️✔️:     NOGIL gil;
    // RDKit✔️✔️:     MolAlign::alignMolConformers(mol, aIds.get(), cIds.get(), wtsVec.get(),
    // RDKit✔️✔️:                                  reflect, maxIters, RMSvector.get());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (RMSvector) {
    // RDKit✔️✔️:     auto &pyl = static_cast<python::list &>(RMSlist);
    // RDKit✔️✔️:     for (double i : *RMSvector) {
    // RDKit✔️✔️:       pyl.append(i);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: std::vector<unsigned int> *translateIntSeq(const python::object &intSeq) {
    // RDKit✔️✔️:   PySequenceHolder<unsigned int> ints(intSeq);
    // RDKit✔️✔️:   std::vector<unsigned int> *intVec = nullptr;
    // RDKit✔️✔️:   if (ints.size() > 0) {
    // RDKit✔️✔️:     intVec = new std::vector<unsigned int>;
    // RDKit✔️✔️:     for (unsigned int i = 0; i < ints.size(); ++i) {
    // RDKit✔️✔️:       intVec->push_back(ints[i]);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return intVec;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: void _fillAtomPositions(RDGeom::Point3DConstPtrVect &pts, const Conformer &conf,
    // RDKit✔️✔️:                         const std::vector<unsigned int> *atomIds = nullptr) {
    // RDKit✔️✔️:   unsigned int na = conf.getNumAtoms();
    // RDKit✔️✔️:   pts.clear();
    // RDKit✔️✔️:   if (atomIds == nullptr) {
    // RDKit✔️✔️:     unsigned int ai;
    // RDKit✔️✔️:     pts.reserve(na);
    // RDKit✔️✔️:     for (ai = 0; ai < na; ++ai) {
    // RDKit✔️✔️:       pts.push_back(&conf.getAtomPos(ai));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     pts.reserve(atomIds->size());
    // RDKit✔️✔️:     std::vector<unsigned int>::const_iterator cai;
    // RDKit✔️✔️:     for (cai = atomIds->begin(); cai != atomIds->end(); cai++) {
    // RDKit✔️✔️:       pts.push_back(&conf.getAtomPos(*cai));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // translateIntSeq([]) returns nullptr: both None and [] select all atoms.
    let atoms: Vec<usize> = match &params.atom_indices {
        Some(indices) if !indices.is_empty() => indices.clone(),
        _ => (0..atom_count).collect(),
    };
    if atoms.is_empty() {
        return Err(AlignmentError::EmptyAtomMap);
    }
    if let Some(&index) = atoms.iter().find(|&&index| index >= atom_count) {
        return Err(AlignmentError::ProbeAtomOutOfRange { index, atom_count });
    }
    let selected: Vec<usize> = match &params.conformer_ids {
        None => (0..coordinates.conformers_3d.len()).collect(),
        Some(ids) if ids.is_empty() => (0..coordinates.conformers_3d.len()).collect(),
        Some(ids) => ids
            .iter()
            .map(|&id| {
                coordinates
                    .conformers_3d
                    .iter()
                    .position(|conformer| conformer.id() == id)
                    .ok_or(AlignmentError::ConformerNotFound { id: id as i32 })
            })
            .collect::<Result<_, _>>()?,
    };
    if selected.is_empty() {
        return Ok(Vec::new());
    }
    // AlignPoints checks weights only when there is a probe conformer.
    if selected.len() > 1 {
        if let Some(weights) = &params.weights {
            if weights.len() != atoms.len() {
                return Err(AlignmentError::WeightCountMismatch {
                    map_len: atoms.len(),
                    weight_len: weights.len(),
                });
            }
            if let Some(index) = weights.iter().position(|weight| !(*weight > 0.0)) {
                return Err(AlignmentError::NonPositiveWeight { index });
            }
        }
    }
    let reference_points: Vec<_> = atoms
        .iter()
        .map(|&atom| coordinates.conformers_3d[selected[0]].coordinates()[atom])
        .collect();
    let mut rmsds = Vec::with_capacity(selected.len().saturating_sub(1));
    for &conformer_index in selected.iter().skip(1) {
        let probe_points: Vec<_> = atoms
            .iter()
            .map(|&atom| coordinates.conformers_3d[conformer_index].coordinates()[atom])
            .collect();
        let (ssr, transform) = align_points(
            &reference_points,
            &probe_points,
            params.weights.as_deref(),
            params.reflect,
            params.max_iterations as usize,
        )
        .map_err(|message| AlignmentError::NumericalPrecondition { message })?;
        rmsds.push((ssr / atoms.len() as f64).sqrt());
        for point in coordinates.conformers_3d[conformer_index].coordinates_mut() {
            *point = alignment_transform_point(&transform, *point);
        }
    }
    Ok(rmsds)
}

pub fn apply_alignment(
    coordinates: &mut CoordinateBlock,
    probe_conformer_id: i32,
    result: &AlignmentResult,
) -> Result<(), AlignmentError> {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8: Code/GraphMol/MolAlign/AlignMolecules.cpp
    // BEGIN VERBATIM CPP alignMol
    // RDKit✔️✔️: double alignMol(ROMol &prbMol, const ROMol &refMol, int prbCid, int refCid,
    // RDKit✔️✔️:                 const MatchVectType *atomMap,
    // RDKit✔️✔️:                 const RDNumeric::DoubleVector *weights, bool reflect,
    // RDKit✔️✔️:                 unsigned int maxIterations) {
    // RDKit✔️✔️:   RDGeom::Transform3D trans;
    // RDKit✔️✔️:   double res = getAlignmentTransform(prbMol, refMol, trans, prbCid, refCid,
    // RDKit✔️✔️:                                      atomMap, weights, reflect, maxIterations);
    // RDKit✔️✔️:   // now transform the relevant conformation on prbMol
    // RDKit✔️✔️:   Conformer &conf = prbMol.getConformer(prbCid);
    // RDKit✔️✔️:   MolTransforms::transformConformer(conf, trans);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END VERBATIM CPP alignMol

    let conformer_index = if probe_conformer_id < 0 {
        (!coordinates.conformers_3d.is_empty()).then_some(0)
    } else {
        coordinates
            .conformers_3d
            .iter()
            .position(|c| c.id() == probe_conformer_id as usize)
    }
    .ok_or_else(|| {
        if coordinates.conformers_3d.is_empty() {
            AlignmentError::NoConformers
        } else {
            AlignmentError::ConformerNotFound {
                id: probe_conformer_id,
            }
        }
    })?;
    for point in coordinates.conformers_3d[conformer_index].coordinates_mut() {
        *point = alignment_transform_point(&result.transform.matrix, *point);
    }
    Ok(())
}

impl AlignmentInput<'_> {
    pub fn alignment_transform_to(
        &self,
        reference: &Self,
        params: &AlignmentParameters,
    ) -> Result<AlignmentResult, AlignmentError> {
        // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8: Code/GraphMol/MolAlign/AlignMolecules.cpp
        // BEGIN VERBATIM CPP getAlignmentTransform
        // RDKit✔️✔️: double getAlignmentTransform(const ROMol &prbMol, const ROMol &refMol,
        // RDKit✔️✔️:                              RDGeom::Transform3D &trans, int prbCid, int refCid,
        // RDKit✔️✔️:                              const MatchVectType *atomMap,
        // RDKit✔️✔️:                              const RDNumeric::DoubleVector *weights,
        // RDKit✔️✔️:                              bool reflect, unsigned int maxIterations) {
        // RDKit✔️✔️:   const Conformer &prbCnf = prbMol.getConformer(prbCid);
        // RDKit✔️✔️:   const Conformer &refCnf = refMol.getConformer(refCid);
        // RDKit✔️✔️:   MatchVectType match;
        // RDKit✔️✔️:   if (!atomMap) {
        // RDKit✔️✔️:     // we have to figure out the mapping between the two molecule
        // RDKit✔️✔️:     const bool recursionPossible = true;
        // RDKit✔️✔️:     const bool useChirality = false;
        // RDKit✔️✔️:     const bool useQueryQueryMatches = true;
        // RDKit✔️✔️:     if (SubstructMatch(refMol, prbMol, match, recursionPossible, useChirality,
        // RDKit✔️✔️:                        useQueryQueryMatches)) {
        // RDKit✔️✔️:       atomMap = &match;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       throw MolAlignException(
        // RDKit✔️✔️:           "No sub-structure match found between the probe and query mol");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   double msd = alignConfsOnAtomMap(prbCnf, refCnf, *atomMap, trans, weights,
        // RDKit✔️✔️:                                    reflect, maxIterations);
        // RDKit✔️✔️:   return sqrt(msd);
        // RDKit✔️✔️: }
        // END VERBATIM CPP getAlignmentTransform

        let probe = conformer_by_id(self, params.probe_conformer_id)?;
        let reference_conformer = conformer_by_id(reference, params.reference_conformer_id)?;
        let map = first_alignment_map(self, reference, params)?;
        aligned_result(
            probe,
            reference_conformer,
            &map,
            params.weights.as_deref(),
            params.reflect,
            params.max_iterations,
        )
    }

    pub fn best_alignment_to(
        &self,
        reference: &Self,
        params: &BestAlignmentParameters,
    ) -> Result<AlignmentResult, AlignmentError> {
        // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8: Code/GraphMol/MolAlign/AlignMolecules.cpp
        // BEGIN VERBATIM CPP getBestAlignmentTransform
        // RDKit✔️✔️: double getBestAlignmentTransform(const ROMol &prbMol, const ROMol &refMol,
        // RDKit✔️✔️:                                  RDGeom::Transform3D &bestTrans,
        // RDKit✔️✔️:                                  MatchVectType &bestMatch,
        // RDKit✔️✔️:                                  const BestAlignmentParams &params, int prbCid,
        // RDKit✔️✔️:                                  int refCid, bool reflect,
        // RDKit✔️✔️:                                  unsigned int maxIters) {
        // RDKit✔️✔️:   std::vector<MatchVectType> allMatches;
        // RDKit✔️✔️:   if (params.map.empty()) {
        // RDKit✔️✔️:     getAllMatchesPrbRef(prbMol, refMol, allMatches, params.maxMatches,
        // RDKit✔️✔️:                         params.symmetrizeConjugatedTerminalGroups,
        // RDKit✔️✔️:                         params.ignoreHs);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   const auto &matches = params.map.empty() ? allMatches : params.map;
        // RDKit✔️✔️:   auto bestRMS = getBestRMSInternal(
        // RDKit✔️✔️:       prbMol, refMol, prbCid, refCid, matches, &bestTrans, &bestMatch,
        // RDKit✔️✔️:       params.weights, reflect, maxIters, getNumThreadsToUse(params.numThreads));
        // RDKit✔️✔️:   return bestRMS;
        // RDKit✔️✔️: }
        // END VERBATIM CPP getBestAlignmentTransform

        let maps = maps_for(
            self,
            reference,
            &params.atom_maps,
            params.max_matches,
            params.symmetrize_conjugated_terminal_groups,
            params.ignore_hydrogens,
        )?;
        let probe = conformer_by_id(self, params.probe_conformer_id)?;
        let reference_conformer = conformer_by_id(reference, params.reference_conformer_id)?;
        best_aligned_result(
            probe,
            reference_conformer,
            &maps,
            params.weights.as_deref(),
            params.reflect,
            params.max_iterations,
            params.num_threads,
        )
    }

    pub fn best_rmsd_to(
        &self,
        reference: &Self,
        params: &BestAlignmentParameters,
    ) -> Result<f64, AlignmentError> {
        Ok(self.best_alignment_to(reference, params)?.rmsd)
    }

    pub fn coordinate_rmsd_to(
        &self,
        reference: &Self,
        params: &CoordinateRmsdParameters,
    ) -> Result<f64, AlignmentError> {
        // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8: Code/GraphMol/MolAlign/AlignMolecules.cpp
        // BEGIN VERBATIM CPP CalcRMS
        // RDKit✔️✔️: double CalcRMS(ROMol &prbMol, const ROMol &refMol, int prbCid, int refCid,
        // RDKit✔️✔️:                const std::vector<MatchVectType> &map, int maxMatches,
        // RDKit✔️✔️:                bool symmetrizeConjugatedTerminalGroups,
        // RDKit✔️✔️:                const RDNumeric::DoubleVector *weights) {
        // RDKit✔️✔️:   std::vector<MatchVectType> allMatches;
        // RDKit✔️✔️:   if (map.empty()) {
        // RDKit✔️✔️:     getAllMatchesPrbRef(prbMol, refMol, allMatches, maxMatches,
        // RDKit✔️✔️:                         symmetrizeConjugatedTerminalGroups);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   const auto &matches = map.empty() ? allMatches : map;
        // RDKit✔️✔️:   bool reflect = false;
        // RDKit✔️✔️:   unsigned int maxIters = 50;
        // RDKit✔️✔️:   unsigned int numThreads = 1;
        // RDKit✔️✔️:   return getBestRMSInternal(prbMol, refMol, prbCid, refCid, matches, nullptr,
        // RDKit✔️✔️:                             nullptr, weights, reflect, maxIters, numThreads);
        // RDKit✔️✔️: }
        // END VERBATIM CPP CalcRMS

        let maps = maps_for(
            self,
            reference,
            &params.atom_maps,
            params.max_matches,
            params.symmetrize_conjugated_terminal_groups,
            false,
        )?;
        let probe = conformer_by_id(self, params.probe_conformer_id)?;
        let reference_conformer = conformer_by_id(reference, params.reference_conformer_id)?;
        let mut best = f64::MAX;
        for map in maps.iter() {
            let value =
                coordinate_rmsd(probe, reference_conformer, &map, params.weights.as_deref())?;
            if value < best {
                best = value;
            }
        }
        Ok(best.sqrt())
    }

    pub fn all_conformer_best_rmsds(
        &self,
        params: &AllConformerRmsdParameters,
    ) -> Result<Vec<ConformerRmsd>, AlignmentError> {
        // Performance loss: source serial mode streams the nested pair loops;
        // this implementation additionally allocates O(C^2) pair indices.

        // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8: Code/GraphMol/MolAlign/AlignMolecules.cpp
        // BEGIN VERBATIM CPP getAllConformerBestRMS
        // RDKit✔️❌: std::vector<double> getAllConformerBestRMS(const ROMol &mol,
        // RDKit✔️❌:                                            const BestAlignmentParams &params) {
        // RDKit✔️❌:   auto numThreads = getNumThreadsToUse(params.numThreads);
        // RDKit✔️❌:   std::vector<MatchVectType> allMatches;
        // RDKit✔️❌:   if (params.map.empty()) {
        // RDKit✔️❌:     getAllMatchesPrbRef(mol, mol, allMatches, params.maxMatches,
        // RDKit✔️❌:                         params.symmetrizeConjugatedTerminalGroups,
        // RDKit✔️❌:                         params.ignoreHs);
        // RDKit✔️❌:   }
        // RDKit✔️❌:   const auto &matches = params.map.empty() ? allMatches : params.map;
        // RDKit✔️❌:
        // RDKit✔️❌:   std::vector<double> res;
        // RDKit✔️❌:   RDGeom::Transform3D trans;
        // RDKit✔️❌:   bool reflect = false;
        // RDKit✔️❌:   unsigned int maxIters = 50;
        // RDKit✔️❌:   std::vector<int> cids;
        // RDKit✔️❌:   for (auto cit = mol.beginConformers(); cit != mol.endConformers(); ++cit) {
        // RDKit✔️❌:     cids.push_back((*cit)->getId());
        // RDKit✔️❌:   }
        // RDKit✔️❌:   if (numThreads == 1) {
        // RDKit✔️❌:     for (auto ci = 0u; ci < mol.getNumConformers(); ++ci) {
        // RDKit✔️❌:       for (auto cj = 0u; cj < ci; ++cj) {
        // RDKit✔️❌:         res.push_back(getBestRMSInternal(mol, mol, cids[ci], cids[cj], matches,
        // RDKit✔️❌:                                          &trans, nullptr, params.weights,
        // RDKit✔️❌:                                          reflect, maxIters, 1));
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // RDKit✔️❌: #ifdef RDK_BUILD_THREADSAFE_SSS
        // RDKit✔️❌:   else {
        // RDKit✔️❌:     std::vector<std::pair<unsigned int, unsigned int>> pairs;
        // RDKit✔️❌:     for (auto ci = 0u; ci < mol.getNumConformers(); ++ci) {
        // RDKit✔️❌:       for (auto cj = 0u; cj < ci; ++cj) {
        // RDKit✔️❌:         pairs.emplace_back(cids[ci], cids[cj]);
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:     std::vector<std::vector<std::pair<unsigned int, double>>> rmsds(numThreads);
        // RDKit✔️❌:     auto func = [&](unsigned int tidx) {
        // RDKit✔️❌:       RDGeom::Transform3D trans;
        // RDKit✔️❌:       bool reflect = false;
        // RDKit✔️❌:       unsigned int maxIters = 50;
        // RDKit✔️❌:       for (auto i = tidx; i < pairs.size(); i += numThreads) {
        // RDKit✔️❌:         auto rms = getBestRMSInternal(mol, mol, pairs[i].first, pairs[i].second,
        // RDKit✔️❌:                                       matches, &trans, nullptr, params.weights,
        // RDKit✔️❌:                                       reflect, maxIters, 1);
        // RDKit✔️❌:         rmsds[tidx].emplace_back(i, rms);
        // RDKit✔️❌:       }
        // RDKit✔️❌:     };
        // RDKit✔️❌:     std::vector<std::thread> tg;
        // RDKit✔️❌:     for (auto ti = 0u; ti < numThreads; ++ti) {
        // RDKit✔️❌:       tg.emplace_back(std::thread(func, ti));
        // RDKit✔️❌:     }
        // RDKit✔️❌:     for (auto &thread : tg) {
        // RDKit✔️❌:       if (thread.joinable()) {
        // RDKit✔️❌:         thread.join();
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:     res.resize(pairs.size());
        // RDKit✔️❌:     for (const auto &tres : rmsds) {
        // RDKit✔️❌:       for (const auto &v : tres) {
        // RDKit✔️❌:         res[v.first] = v.second;
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // RDKit✔️❌: #endif
        // RDKit✔️❌:   return res;
        // RDKit✔️❌: }
        // END VERBATIM CPP getAllConformerBestRMS

        let conformers = self.conformers_3d();
        let maps = maps_for(
            self,
            self,
            &params.atom_maps,
            params.max_matches,
            params.symmetrize_conjugated_terminal_groups,
            params.ignore_hydrogens,
        )?;
        let mut pairs =
            Vec::with_capacity(conformers.len() * conformers.len().saturating_sub(1) / 2);
        for probe_index in 0..conformers.len() {
            for reference_index in 0..probe_index {
                pairs.push((probe_index, reference_index));
            }
        }
        let evaluate = |pair_index: usize| -> Result<ConformerRmsd, AlignmentError> {
            let (probe_index, reference_index) = pairs[pair_index];
            let aligned = best_aligned_result(
                &conformers[probe_index],
                &conformers[reference_index],
                &maps,
                params.weights.as_deref(),
                false,
                DEFAULT_MAX_ITERATIONS,
                1,
            )?;
            Ok(ConformerRmsd {
                probe_conformer_id: conformers[probe_index].id(),
                reference_conformer_id: conformers[reference_index].id(),
                rmsd: aligned.rmsd,
            })
        };
        let thread_count = num_threads_to_use(params.num_threads)?;
        if thread_count == 1 {
            return (0..pairs.len()).map(evaluate).collect();
        }
        let pair_count = pairs.len();
        let buckets = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(thread_count);
            for thread_index in 0..thread_count {
                let evaluate = &evaluate;
                handles.push(scope.spawn(move || {
                    let mut bucket = Vec::new();
                    let mut pair_index = thread_index;
                    while pair_index < pair_count {
                        bucket.push((pair_index, evaluate(pair_index)));
                        pair_index += thread_count;
                    }
                    bucket
                }));
            }
            handles
                .into_iter()
                .map(|handle| handle.join())
                .collect::<Vec<_>>()
        });
        let mut ordered = vec![None; pair_count];
        for bucket in buckets {
            for (pair_index, result) in bucket.map_err(|_| AlignmentError::WorkerTerminated)? {
                ordered[pair_index] = Some(result);
            }
        }
        ordered
            .into_iter()
            .map(|entry| entry.ok_or(AlignmentError::WorkerTerminated)?)
            .collect()
    }
}

impl AlignmentResult {
    pub const fn rmsd(&self) -> f64 {
        self.rmsd
    }
    pub const fn transform(&self) -> &AlignmentTransform {
        &self.transform
    }
    pub fn atom_map(&self) -> &[AlignmentAtomMap] {
        &self.atom_map
    }
}
impl AlignmentTransform {
    pub const fn matrix(&self) -> &[[f64; 4]; 4] {
        &self.matrix
    }
}
impl ConformerRmsd {
    pub const fn probe_conformer_id(&self) -> usize {
        self.probe_conformer_id
    }
    pub const fn reference_conformer_id(&self) -> usize {
        self.reference_conformer_id
    }
    pub const fn rmsd(&self) -> f64 {
        self.rmsd
    }
}
impl ConformerAlignmentReport {
    pub fn rmsds(&self) -> &[f64] {
        &self.rmsds
    }
}

impl AlignmentAtomMap {
    pub const fn new(probe_atom: usize, reference_atom: usize) -> Self {
        Self {
            probe_atom,
            reference_atom,
        }
    }
}
impl AlignmentParameters {
    pub fn new(
        probe_conformer_id: i32,
        reference_conformer_id: i32,
        atom_map: Option<Vec<AlignmentAtomMap>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
    ) -> Self {
        Self {
            probe_conformer_id,
            reference_conformer_id,
            atom_map,
            weights,
            reflect,
            max_iterations,
        }
    }
}
impl BestAlignmentParameters {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        probe_conformer_id: i32,
        reference_conformer_id: i32,
        atom_maps: Vec<Vec<AlignmentAtomMap>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
        max_matches: i32,
        symmetrize_conjugated_terminal_groups: bool,
        ignore_hydrogens: bool,
        num_threads: i32,
    ) -> Self {
        Self {
            probe_conformer_id,
            reference_conformer_id,
            atom_maps,
            weights,
            reflect,
            max_iterations,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            ignore_hydrogens,
            num_threads,
        }
    }
}
impl CoordinateRmsdParameters {
    pub fn new(
        probe_conformer_id: i32,
        reference_conformer_id: i32,
        atom_maps: Vec<Vec<AlignmentAtomMap>>,
        weights: Option<Vec<f64>>,
        max_matches: i32,
        symmetrize_conjugated_terminal_groups: bool,
    ) -> Self {
        Self {
            probe_conformer_id,
            reference_conformer_id,
            atom_maps,
            weights,
            max_matches,
            symmetrize_conjugated_terminal_groups,
        }
    }
}
impl AllConformerRmsdParameters {
    pub fn new(
        atom_maps: Vec<Vec<AlignmentAtomMap>>,
        weights: Option<Vec<f64>>,
        max_matches: i32,
        symmetrize_conjugated_terminal_groups: bool,
        ignore_hydrogens: bool,
        num_threads: i32,
    ) -> Self {
        Self {
            atom_maps,
            weights,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            ignore_hydrogens,
            num_threads,
        }
    }
}
impl ConformerAlignmentParameters {
    pub fn new(
        atom_indices: Option<Vec<usize>>,
        conformer_ids: Option<Vec<usize>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
    ) -> Self {
        Self {
            atom_indices,
            conformer_ids,
            weights,
            reflect,
            max_iterations,
        }
    }
}
