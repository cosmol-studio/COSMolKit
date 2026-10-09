//! Source-backed Layered fingerprint over detached concrete topology.
use crate::{Fingerprint, FingerprintError, hash::hash_range};
use cosmolkit_core::{
    GraphPath, PathError, PathRepresentation, PathSearchParams, SubgraphSearchParams,
    all_paths_in_range, all_subgraphs_in_range, query_atom_paths_in_range,
    query_bond_paths_in_range, query_subgraphs_in_range,
};
use cosmolkit_core::{
    RingFindingError, RingInfo, RingSearchParams, find_sssr, find_sssr_from_parts,
};
use cosmolkit_model::BondOrder;
use cosmolkit_model::{
    AdjacencyList, AtomId, Bond, CoordinateBlock, QueryGraph, QueryGraphError, TopologyBlock,
    TopologyValidationError,
};
use cosmolkit_search::{
    SearchTarget, is_atom_aromatic, is_complex_atom_query, is_complex_bond_query,
    is_complex_concrete_bond_query, is_query_atom_aromatic,
};
use std::borrow::Cow;
use std::collections::BTreeMap;
use std::fmt;
use std::ops::{BitOr, BitOrAssign};

#[derive(Debug)]
pub enum LayeredFingerprintError {
    InvalidArguments { reason: &'static str },
    Topology(TopologyValidationError),
    Query(QueryGraphError),
    Rings(RingFindingError),
    Paths(PathError),
    Value(FingerprintError),
}
impl fmt::Display for LayeredFingerprintError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidArguments { reason } => f.write_str(reason),
            Self::Topology(e) => e.fmt(f),
            Self::Query(e) => e.fmt(f),
            Self::Rings(e) => e.fmt(f),
            Self::Paths(e) => e.fmt(f),
            Self::Value(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for LayeredFingerprintError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Topology(e) => Some(e),
            Self::Query(e) => Some(e),
            Self::Rings(e) => Some(e),
            Self::Paths(e) => Some(e),
            Self::Value(e) => Some(e),
            Self::InvalidArguments { .. } => None,
        }
    }
}
impl From<TopologyValidationError> for LayeredFingerprintError {
    fn from(e: TopologyValidationError) -> Self {
        Self::Topology(e)
    }
}
impl From<RingFindingError> for LayeredFingerprintError {
    fn from(e: RingFindingError) -> Self {
        Self::Rings(e)
    }
}
impl From<QueryGraphError> for LayeredFingerprintError {
    fn from(e: QueryGraphError) -> Self {
        Self::Query(e)
    }
}
impl From<PathError> for LayeredFingerprintError {
    fn from(e: PathError) -> Self {
        Self::Paths(e)
    }
}
impl From<FingerprintError> for LayeredFingerprintError {
    fn from(e: FingerprintError) -> Self {
        Self::Value(e)
    }
}

// BEGIN RDKIT CPP CONSTANTS LayeredFingerprintMol metadata
// RDKit✔️✔️: const unsigned int maxFingerprintLayers = 10;
pub const LAYERED_FINGERPRINT_MAX_LAYERS: usize = 10;
// RDKit✔️✔️: const std::string LayeredFingerprintMolVersion = "0.7.0";
pub const LAYERED_FINGERPRINT_VERSION: &str = "0.7.0";
// RDKit✔️✔️: const unsigned int substructLayers = 0x07;
pub const LAYERED_FINGERPRINT_SUBSTRUCTURE_LAYERS: u32 = 0x07;
// END RDKIT CPP CONSTANTS LayeredFingerprintMol metadata

/// Source layer flags for RDKit's experimental Layered fingerprint algorithm.
///
/// Unknown/high source bits are retained because the C++ API accepts the full
/// `unsigned int` value and silently produces no components for unimplemented
/// layer slots. Only the six named layers currently emit components:
/// topology (`0x01`), bond order (`0x02`), atom type (`0x04`), ring presence
/// (`0x08`), minimum ring size (`0x10`), and aromaticity (`0x20`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LayeredFingerprintLayers(u32);

impl LayeredFingerprintLayers {
    pub const TOPOLOGY: Self = Self(0x01);
    pub const BOND_ORDER: Self = Self(0x02);
    pub const ATOM_TYPE: Self = Self(0x04);
    pub const RING_PRESENCE: Self = Self(0x08);
    pub const RING_SIZE: Self = Self(0x10);
    pub const AROMATICITY: Self = Self(0x20);
    pub const ACTIVE: Self = Self(0x3f);
    pub const SUBSTRUCTURE: Self = Self(LAYERED_FINGERPRINT_SUBSTRUCTURE_LAYERS);
    pub const ALL_SOURCE_BITS: Self = Self(u32::MAX);

    #[must_use]
    pub const fn bits(self) -> u32 {
        self.0
    }

    #[must_use]
    pub const fn from_bits_retain(bits: u32) -> Self {
        Self(bits)
    }

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        self.0 & other.0 == other.0
    }
}

impl BitOr for LayeredFingerprintLayers {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl BitOrAssign for LayeredFingerprintLayers {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

/// Parameters for the source-backed Layered fingerprint API.
///
/// The defaults reproduce the source wrapper: all source flag bits, bond-path
/// lengths 1 through 7, 2,048 output bits, branched path enumeration, no atom
/// counts, no output-bit mask, and no root selection. `from_atoms: None`
/// selects the whole graph, while `Some(Vec::new())` is a present but empty
/// root selection and therefore enumerates no paths.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LayeredFingerprintParams {
    /// Layer flags. Unknown high bits are retained and emit no components.
    pub layers: LayeredFingerprintLayers,
    /// Minimum path length. The unrooted linear source branch counts atoms;
    /// branched and rooted paths count bonds. Zero is rejected as `minPath==0`.
    pub min_path: u32,
    /// Maximum path length, in the same units as `min_path`.
    /// Values below `min_path` are rejected.
    pub max_path: u32,
    /// Explicit output width. Zero is rejected.
    pub fp_size: u32,
    /// Optional seeded source `atomCounts` vector. Values are incremented and
    /// returned without clearing the caller-provided seed.
    pub atom_counts: Option<Vec<u32>>,
    /// Optional projection mask, which must have exactly `fp_size` bits.
    pub set_only_bits: Option<Fingerprint>,
    /// Enumerate branched subgraphs when true and linear paths when false.
    /// Unrooted linear paths preserve the pinned source's atom-index selection
    /// and subsequent bond-index interpretation; invalid accesses return an error.
    pub branched_paths: bool,
    /// `None` is an absent source pointer; `Some(Vec::new())` is a present
    /// empty selection and therefore enumerates no paths.
    pub from_atoms: Option<Vec<u32>>,
}

impl Default for LayeredFingerprintParams {
    fn default() -> Self {
        // RDKit✔️✔️:     const ROMol &mol, unsigned int layerFlags = 0xFFFFFFFF,
        // RDKit✔️✔️:     unsigned int minPath = 1, unsigned int maxPath = 7,
        // RDKit✔️✔️:     unsigned int fpSize = 2048, std::vector<unsigned int> *atomCounts = nullptr,
        // RDKit✔️✔️:     ExplicitBitVect *setOnlyBits = nullptr, bool branchedPaths = true,
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *fromAtoms = nullptr);
        Self {
            layers: LayeredFingerprintLayers::ALL_SOURCE_BITS,
            min_path: 1,
            max_path: 7,
            fp_size: 2048,
            atom_counts: None,
            set_only_bits: None,
            branched_paths: true,
            from_atoms: None,
        }
    }
}

impl LayeredFingerprintParams {
    pub fn validate(&self) -> Result<(), LayeredFingerprintError> {
        // RDKit✔️✔️:   PRECONDITION(minPath != 0, "minPath==0");
        // RDKit✔️✔️:   PRECONDITION(maxPath >= minPath, "maxPath<minPath");
        // RDKit✔️✔️:   PRECONDITION(fpSize != 0, "fpSize==0");
        if self.min_path == 0 {
            return Err(LayeredFingerprintError::InvalidArguments {
                reason: "minPath==0",
            });
        }
        if self.max_path < self.min_path {
            return Err(LayeredFingerprintError::InvalidArguments {
                reason: "maxPath<minPath",
            });
        }
        if self.fp_size == 0 {
            return Err(LayeredFingerprintError::InvalidArguments {
                reason: "fpSize==0",
            });
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LayeredFingerprintResult {
    /// The fixed-width Layered bit vector.
    pub fingerprint: Fingerprint,
    /// Updated seeded counts, or `None` when counts were not requested.
    ///
    /// Every atom in an accepted path is incremented once for that path, even
    /// when several active layers set bits or several projections collide.
    pub atom_counts: Option<Vec<u32>>,
}

impl LayeredFingerprintResult {
    pub fn fingerprint(&self) -> &Fingerprint {
        &self.fingerprint
    }
    pub fn atom_counts(&self) -> Option<&[u32]> {
        self.atom_counts.as_deref()
    }
}

#[derive(Clone, Copy)]
pub(super) enum LayeredGraphInput<'a> {
    Concrete(&'a TopologyBlock),
    Query(&'a QueryGraph),
}
impl<'a> LayeredGraphInput<'a> {
    fn atom_count(self) -> usize {
        match self {
            Self::Concrete(t) => t.atoms.len(),
            Self::Query(q) => q.num_atoms(),
        }
    }
    fn bond_count(self) -> usize {
        match self {
            Self::Concrete(t) => t.bonds.len(),
            Self::Query(q) => q.num_bonds(),
        }
    }
    fn bond(self, index: usize) -> &'a Bond {
        match self {
            Self::Concrete(t) => &t.bonds[index],
            Self::Query(q) => q.bonds()[index].bond(),
        }
    }
}

fn enumerate_fingerprint_paths_for_root(
    graph: LayeredGraphInput<'_>,
    lower: usize,
    upper: usize,
    use_hs: bool,
    branched_paths: bool,
    root: Option<u32>,
    ignore_atoms: Option<&[bool]>,
) -> Result<BTreeMap<usize, Vec<Vec<usize>>>, LayeredFingerprintError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::RDKitFPUtils::enumerateAllPaths (Release_2026_03_6)
    // RDKit❗✔️: void enumerateAllPaths(const ROMol &mol, INT_PATH_LIST_MAP &allPaths,
    // RDKit❗✔️:                        const std::vector<std::uint32_t> *fromAtoms,
    // RDKit❗✔️:                        bool branchedPaths, bool useHs, unsigned int minPath,
    // RDKit❗✔️:                        unsigned int maxPath,
    // RDKit❗✔️:                        boost::dynamic_bitset<> *ignoreAtoms) {
    // RDKit❗✔️:   PRECONDITION(!ignoreAtoms || ignoreAtoms->size() == mol.getNumAtoms(),
    // RDKit❗✔️:                "bad ignoreAtoms size");
    // RDKit❗✔️:   if (!fromAtoms) {
    // RDKit❗✔️:     if (branchedPaths) {
    // RDKit❗✔️:       int rootedAtAtom = -1;
    // RDKit❗✔️:       allPaths = findAllSubgraphsOfLengthsMtoN(mol, minPath, maxPath, useHs,
    // RDKit❗✔️:                                                rootedAtAtom, ignoreAtoms);
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       bool useBonds = true;
    // RDKit❗✔️:       int rootedAtAtom = -1;
    // RDKit❗✔️:       bool onlyShortestPaths = false;
    // RDKit❗✔️:       allPaths = findAllPathsOfLengthsMtoN(mol, minPath, maxPath, useBonds,
    // RDKit❗✔️:                                            useHs, rootedAtAtom,
    // RDKit❗✔️:                                            onlyShortestPaths, ignoreAtoms);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     for (auto aidx : *fromAtoms) {
    // RDKit❗✔️:       INT_PATH_LIST_MAP tPaths;
    // RDKit❗✔️:       if (branchedPaths) {
    // RDKit❗✔️:         tPaths = findAllSubgraphsOfLengthsMtoN(mol, minPath, maxPath, useHs,
    // RDKit❗✔️:                                                aidx, ignoreAtoms);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         bool useBonds = true;
    // RDKit❗✔️:         bool onlyShortestPaths = false;
    // RDKit❗✔️:         tPaths =
    // RDKit❗✔️:             findAllPathsOfLengthsMtoN(mol, minPath, maxPath, useBonds, useHs,
    // RDKit❗✔️:                                       aidx, onlyShortestPaths, ignoreAtoms);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       for (INT_PATH_LIST_MAP::const_iterator tpit = tPaths.begin();
    // RDKit❗✔️:            tpit != tPaths.end(); ++tpit) {
    // RDKit❗✔️: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗✔️:         std::cerr << "paths from " << aidx << " size: " << tpit->first
    // RDKit❗✔️:                   << std::endl;
    // RDKit❗✔️:         for (auto path : tpit->second) {
    // RDKit❗✔️:           std::cerr << " path: ";
    // RDKit❗✔️:           std::copy(path.begin(), path.end(),
    // RDKit❗✔️:                     std::ostream_iterator<int>(std::cerr, ", "));
    // RDKit❗✔️:           std::cerr << std::endl;
    // RDKit❗✔️:         }
    // RDKit❗✔️: #endif
    // RDKit❗✔️:
    // RDKit❗✔️:         allPaths[tpit->first].insert(allPaths[tpit->first].begin(),
    // RDKit❗✔️:                                      tpit->second.begin(), tpit->second.end());
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::RDKitFPUtils::enumerateAllPaths
    // Source enumeration belongs to CORE paths; IDs are moved, never cloned.
    let rooted_at_atom = root.map(|index| AtomId::new(index as usize));
    if branched_paths {
        let params = SubgraphSearchParams {
            use_hydrogens: use_hs,
            rooted_at_atom,
            ignore_atoms,
        };
        let paths = match graph {
            LayeredGraphInput::Concrete(topology) => {
                all_subgraphs_in_range(topology, lower, upper, &params)?
            }
            LayeredGraphInput::Query(query) => {
                query_subgraphs_in_range(query, lower, upper, &params)?
            }
        };
        return Ok(paths
            .into_iter()
            .map(|(size, rows)| {
                (
                    size,
                    rows.into_iter()
                        .map(|row| row.into_iter().map(|id| id.index()).collect())
                        .collect(),
                )
            })
            .collect());
    }
    let paths = match graph {
        LayeredGraphInput::Concrete(topology) => all_paths_in_range(
            topology,
            lower,
            upper,
            &PathSearchParams {
                use_hydrogens: use_hs,
                rooted_at_atom,
                ignore_atoms,
                ..PathSearchParams::default()
            },
        )?,
        LayeredGraphInput::Query(query) => query_bond_paths_in_range(
            query,
            lower,
            upper,
            &SubgraphSearchParams {
                use_hydrogens: use_hs,
                rooted_at_atom,
                ignore_atoms,
            },
        )?,
    };
    paths
        .into_iter()
        .map(|(size, rows)| {
            let rows = rows
                .into_iter()
                .map(|row| match row {
                    GraphPath::Bonds(ids) => Ok(ids.into_iter().map(|id| id.index()).collect()),
                    GraphPath::Atoms(_) => Err(LayeredFingerprintError::InvalidArguments {
                        reason: "bond path enumeration returned atom path",
                    }),
                })
                .collect::<Result<Vec<_>, _>>()?;
            Ok((size, rows))
        })
        .collect()
}

pub(super) fn enumerate_fingerprint_paths(
    graph: LayeredGraphInput<'_>,
    min_path: u32,
    max_path: u32,
    use_hs: bool,
    branched_paths: bool,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[bool]>,
) -> Result<BTreeMap<usize, Vec<Vec<usize>>>, LayeredFingerprintError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::RDKitFPUtils::enumerateAllPaths (Release_2026_03_6)
    // RDKit❗✔️: void enumerateAllPaths(const ROMol &mol, INT_PATH_LIST_MAP &allPaths,
    // RDKit❗✔️:                        const std::vector<std::uint32_t> *fromAtoms,
    // RDKit❗✔️:                        bool branchedPaths, bool useHs, unsigned int minPath,
    // RDKit❗✔️:                        unsigned int maxPath,
    // RDKit❗✔️:                        boost::dynamic_bitset<> *ignoreAtoms) {
    // RDKit❗✔️:   PRECONDITION(!ignoreAtoms || ignoreAtoms->size() == mol.getNumAtoms(),
    // RDKit❗✔️:                "bad ignoreAtoms size");
    // RDKit❗✔️:   if (!fromAtoms) {
    // RDKit❗✔️:     if (branchedPaths) {
    // RDKit❗✔️:       int rootedAtAtom = -1;
    // RDKit❗✔️:       allPaths = findAllSubgraphsOfLengthsMtoN(mol, minPath, maxPath, useHs,
    // RDKit❗✔️:                                                rootedAtAtom, ignoreAtoms);
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       bool useBonds = true;
    // RDKit❗✔️:       int rootedAtAtom = -1;
    // RDKit❗✔️:       bool onlyShortestPaths = false;
    // RDKit❗✔️:       allPaths = findAllPathsOfLengthsMtoN(mol, minPath, maxPath, useBonds,
    // RDKit❗✔️:                                            useHs, rootedAtAtom,
    // RDKit❗✔️:                                            onlyShortestPaths, ignoreAtoms);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     for (auto aidx : *fromAtoms) {
    // RDKit❗✔️:       INT_PATH_LIST_MAP tPaths;
    // RDKit❗✔️:       if (branchedPaths) {
    // RDKit❗✔️:         tPaths = findAllSubgraphsOfLengthsMtoN(mol, minPath, maxPath, useHs,
    // RDKit❗✔️:                                                aidx, ignoreAtoms);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         bool useBonds = true;
    // RDKit❗✔️:         bool onlyShortestPaths = false;
    // RDKit❗✔️:         tPaths =
    // RDKit❗✔️:             findAllPathsOfLengthsMtoN(mol, minPath, maxPath, useBonds, useHs,
    // RDKit❗✔️:                                       aidx, onlyShortestPaths, ignoreAtoms);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       for (INT_PATH_LIST_MAP::const_iterator tpit = tPaths.begin();
    // RDKit❗✔️:            tpit != tPaths.end(); ++tpit) {
    // RDKit❗✔️: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗✔️:         std::cerr << "paths from " << aidx << " size: " << tpit->first
    // RDKit❗✔️:                   << std::endl;
    // RDKit❗✔️:         for (auto path : tpit->second) {
    // RDKit❗✔️:           std::cerr << " path: ";
    // RDKit❗✔️:           std::copy(path.begin(), path.end(),
    // RDKit❗✔️:                     std::ostream_iterator<int>(std::cerr, ", "));
    // RDKit❗✔️:           std::cerr << std::endl;
    // RDKit❗✔️:         }
    // RDKit❗✔️: #endif
    // RDKit❗✔️:
    // RDKit❗✔️:         allPaths[tpit->first].insert(allPaths[tpit->first].begin(),
    // RDKit❗✔️:                                      tpit->second.begin(), tpit->second.end());
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::RDKitFPUtils::enumerateAllPaths
    // Local complexity review: each source call and Rust helper enumerates the
    // same path/subgraph state once per requested root. Both prepend each root
    // group to the per-length vector; no molecule or completed path map is
    // cloned beyond the source-equivalent recursive path state.
    if min_path == 0 || max_path < min_path {
        return Err(LayeredFingerprintError::InvalidArguments {
            reason: "invalid path lengths",
        });
    }
    if ignore_atoms.is_some_and(|mask| mask.len() != graph.atom_count()) {
        return Err(LayeredFingerprintError::InvalidArguments {
            reason: "bad ignoreAtoms size",
        });
    }
    let lower = min_path as usize;
    let upper = max_path as usize;
    let Some(roots) = from_atoms else {
        return enumerate_fingerprint_paths_for_root(
            graph,
            lower,
            upper,
            use_hs,
            branched_paths,
            None,
            ignore_atoms,
        );
    };

    let mut result: BTreeMap<usize, Vec<Vec<usize>>> = BTreeMap::new();
    for &root in roots {
        let rooted_paths = enumerate_fingerprint_paths_for_root(
            graph,
            lower,
            upper,
            use_hs,
            branched_paths,
            Some(root),
            ignore_atoms,
        )?;
        for (length, paths) in rooted_paths {
            result.entry(length).or_default().splice(0..0, paths);
        }
    }
    Ok(result)
}

#[derive(Debug)]
struct LayeredFingerprintPreparation<'a> {
    ring_info: Cow<'a, RingInfo>,
    bond_cache: Vec<&'a Bond>,
    query_masks: Vec<u8>,
    aromatic_atoms: Vec<bool>,
    atomic_numbers: Vec<u32>,
}

fn prepare_layered_fingerprint<'a>(
    graph: LayeredGraphInput<'a>,
    cached_rings: Option<&'a RingInfo>,
    min_path: u32,
    max_path: u32,
    fp_size: usize,
    atom_counts: Option<&[u32]>,
    set_only_bits: Option<&Fingerprint>,
) -> Result<LayeredFingerprintPreparation<'a>, LayeredFingerprintError> {
    // BEGIN RDKIT CPP FUNCTION LayeredFingerprintMol preparation
    // RDKit✔️✔️:   PRECONDITION(minPath != 0, "minPath==0");
    // RDKit✔️✔️:   PRECONDITION(maxPath >= minPath, "maxPath<minPath");
    // RDKit✔️✔️:   PRECONDITION(fpSize != 0, "fpSize==0");
    // RDKit✔️✔️:   PRECONDITION(!atomCounts || atomCounts->size() >= mol.getNumAtoms(),
    // RDKit✔️✔️:                "bad atomCounts size");
    // RDKit✔️✔️:   PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
    // RDKit✔️✔️:                "bad setOnlyBits size");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isInitialized()) {
    // RDKit✔️✔️:     MolOps::findSSSR(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<const Bond *> bondCache;
    // RDKit✔️✔️:   bondCache.resize(mol.getNumBonds());
    // RDKit✔️✔️:   std::vector<short> isQueryBond(mol.getNumBonds(), 0);
    // RDKit✔️✔️:   ROMol::EDGE_ITER firstB, lastB;
    // RDKit✔️✔️:   boost::tie(firstB, lastB) = mol.getEdges();
    // RDKit✔️✔️:   while (firstB != lastB) {
    // RDKit✔️✔️:     const Bond *bond = mol[*firstB];
    // RDKit✔️✔️:     isQueryBond[bond->getIdx()] = 0x0;
    // RDKit✔️✔️:     bondCache[bond->getIdx()] = bond;
    // RDKit✔️✔️:     if (isComplexQuery(bond)) {
    // RDKit✔️✔️:       isQueryBond[bond->getIdx()] = 0x1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (isComplexQuery(bond->getBeginAtom())) {
    // RDKit✔️✔️:       isQueryBond[bond->getIdx()] |= 0x2;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (isComplexQuery(bond->getEndAtom())) {
    // RDKit✔️✔️:       isQueryBond[bond->getIdx()] |= 0x4;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++firstB;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<bool> aromaticAtoms(mol.getNumAtoms(), false);
    // RDKit✔️✔️:   std::vector<int> anums(mol.getNumAtoms(), 0);
    // RDKit✔️✔️:   ROMol::VERTEX_ITER firstA, lastA;
    // RDKit✔️✔️:   boost::tie(firstA, lastA) = mol.getVertices();
    // RDKit✔️✔️:   while (firstA != lastA) {
    // RDKit✔️✔️:     const Atom *atom = mol[*firstA];
    // RDKit✔️✔️:     if (isAtomAromatic(atom)) {
    // RDKit✔️✔️:       aromaticAtoms[atom->getIdx()] = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     anums[atom->getIdx()] = atom->getAtomicNum();
    // RDKit✔️✔️:     ++firstA;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION LayeredFingerprintMol preparation
    // Local complexity review: the five parameter preconditions are O(1).
    // Structural topology/query validation below is O(A+B). An initialized
    // ring cache is borrowed; cold preparation runs CORE exact SSSR. Bond/atom
    // caches require O(B)/O(A) fills; cold query preparation has the additional
    // carrier and adjacency allocation described at its branch below.
    if min_path == 0 {
        return Err(LayeredFingerprintError::InvalidArguments {
            reason: "minPath==0",
        });
    }
    if max_path < min_path {
        return Err(LayeredFingerprintError::InvalidArguments {
            reason: "maxPath<minPath",
        });
    }
    if fp_size == 0 {
        return Err(LayeredFingerprintError::InvalidArguments {
            reason: "fpSize==0",
        });
    }
    if atom_counts.is_some_and(|counts| counts.len() < graph.atom_count()) {
        return Err(LayeredFingerprintError::InvalidArguments {
            reason: "bad atomCounts size",
        });
    }
    if set_only_bits.is_some_and(|bits| bits.n_bits() as usize != fp_size) {
        return Err(LayeredFingerprintError::InvalidArguments {
            reason: "bad setOnlyBits size",
        });
    }

    match graph {
        LayeredGraphInput::Concrete(topology) => topology.validate()?,
        LayeredGraphInput::Query(query) => query.validate()?,
    }
    let ring_info = match cached_rings.filter(|rings| rings.is_initialized()) {
        Some(rings) => Cow::Borrowed(rings),
        None => Cow::Owned(match graph {
            LayeredGraphInput::Concrete(topology) => {
                find_sssr(topology, &RingSearchParams::default())?
            }
            LayeredGraphInput::Query(query) => {
                // RDKit✔️❌:     MolOps::findSSSR(mol);
                // Existing ring owner requires a contiguous Bond slice and
                // canonical adjacency. The cold QueryGraph branch copies B
                // bond carriers and allocates indexed adjacency O(A+B); it
                // preserves all real identities/predicates in the original
                // graph and never creates a concrete Atom or Molecule.
                // This extra cold scratch is a source-cost difference.
                let bonds: Vec<_> = query.bonds().iter().map(|row| row.bond().clone()).collect();
                let adjacency = AdjacencyList::try_from_topology(query.num_atoms(), &bonds)
                    .map_err(|e| LayeredFingerprintError::Rings(RingFindingError::Adjacency(e)))?;
                find_sssr_from_parts(query.num_atoms(), &bonds, &adjacency)?
            }
        }),
    };
    if ring_info.atom_row_count() != graph.atom_count()
        || ring_info.bond_row_count() != graph.bond_count()
    {
        return Err(LayeredFingerprintError::InvalidArguments {
            reason: "ring rows do not match topology",
        });
    }
    let bond_cache = (0..graph.bond_count())
        .map(|index| graph.bond(index))
        .collect();
    let (query_masks, aromatic_atoms, atomic_numbers) = match graph {
        LayeredGraphInput::Concrete(topology) => {
            // Concrete carriers retain explicit parser bond-query identity;
            // SEARCH alone evaluates source complexity from that borrowed root.
            let coordinates = CoordinateBlock::default();
            let target = SearchTarget::new(
                topology,
                &coordinates,
                &topology.stereo_groups,
                Some(&ring_info),
                None,
            );
            (
                topology
                    .bonds
                    .iter()
                    .map(|bond| u8::from(is_complex_concrete_bond_query(bond)))
                    .collect(),
                topology
                    .atoms
                    .iter()
                    .map(|atom| is_atom_aromatic(atom, &target))
                    .collect(),
                topology
                    .atoms
                    .iter()
                    .map(|atom| u32::from(atom.atomic_number()))
                    .collect(),
            )
        }
        LayeredGraphInput::Query(query) => {
            let masks = query
                .bonds()
                .iter()
                .map(|bond| {
                    u8::from(is_complex_bond_query(bond))
                        | (u8::from(is_complex_atom_query(&query.atoms()[bond.begin().index()]))
                            << 1)
                        | (u8::from(is_complex_atom_query(&query.atoms()[bond.end().index()])) << 2)
                })
                .collect();
            (
                masks,
                query
                    .atoms()
                    .iter()
                    .map(|atom| is_query_atom_aromatic(atom, query))
                    .collect(),
                query
                    .atoms()
                    .iter()
                    .map(|atom| u32::from(atom.atomic_number()))
                    .collect(),
            )
        }
    };

    Ok(LayeredFingerprintPreparation {
        ring_info,
        bond_cache,
        query_masks,
        aromatic_atoms,
        atomic_numbers,
    })
}

#[inline]
fn layered_topology_hash(
    bond_neighbor_count: u32,
    begin_atom_degree: u32,
    end_atom_degree: u32,
) -> u32 {
    // BEGIN RDKIT CPP BLOCK LayeredFingerprintMol layer 1
    // RDKit✔️✔️:         if (layerFlags & 0x1) {
    // RDKit✔️✔️:           // layer 1: straight topology
    // RDKit✔️✔️:           unsigned int a1Deg, a2Deg;
    // RDKit✔️✔️:           a1Deg = atomDegrees[bi->getBeginAtomIdx()];
    // RDKit✔️✔️:           a2Deg = atomDegrees[bi->getEndAtomIdx()];
    // RDKit✔️✔️:           if (a1Deg < a2Deg) {
    // RDKit✔️✔️:             std::swap(a1Deg, a2Deg);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           ourHash = bondNbrs[i] % 8;  // 3 bits here
    // RDKit✔️✔️:           ourHash |= (a1Deg % 8) << 3;
    // RDKit✔️✔️:           ourHash |= (a2Deg % 8) << 6;
    // RDKit✔️✔️:           hashLayers[0].push_back(ourHash);
    // RDKit✔️✔️:         }
    // END RDKIT CPP BLOCK LayeredFingerprintMol layer 1
    // Local complexity review: both forms perform three modulo operations,
    // one conditional swap, two shifts, and two ORs in O(1), with no lookup,
    // allocation, clone, or branch beyond the source branch.
    let (larger_degree, smaller_degree) = if begin_atom_degree < end_atom_degree {
        (end_atom_degree, begin_atom_degree)
    } else {
        (begin_atom_degree, end_atom_degree)
    };
    (bond_neighbor_count % 8) | ((larger_degree % 8) << 3) | ((smaller_degree % 8) << 6)
}

#[inline]
fn layered_bond_order_hash(
    bond: &Bond,
    bond_neighbor_count: u32,
    begin_atom_degree: u32,
    end_atom_degree: u32,
    path_queries: u8,
) -> Option<u32> {
    // BEGIN RDKIT CPP BLOCK LayeredFingerprintMol layer 2
    // RDKit✔️✔️:         if (layerFlags & 0x2 && !(pathQueries & 0x1)) {
    // RDKit✔️✔️:           // layer 2: include bond orders:
    // RDKit✔️✔️:           unsigned int bondHash;
    // RDKit✔️✔️:           // makes sure aromatic bonds and single bonds  always hash the same:
    // RDKit✔️✔️:           if (!bi->getIsAromatic() && bi->getBondType() != Bond::SINGLE &&
    // RDKit✔️✔️:               bi->getBondType() != Bond::AROMATIC) {
    // RDKit✔️✔️:             bondHash = bi->getBondType();
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             bondHash = Bond::SINGLE;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           unsigned int a1Deg, a2Deg;
    // RDKit✔️✔️:           a1Deg = atomDegrees[bi->getBeginAtomIdx()];
    // RDKit✔️✔️:           a2Deg = atomDegrees[bi->getEndAtomIdx()];
    // RDKit✔️✔️:           if (a1Deg < a2Deg) {
    // RDKit✔️✔️:             std::swap(a1Deg, a2Deg);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           ourHash = bondHash % 8;
    // RDKit✔️✔️:           ourHash |= (bondNbrs[i] % 8) << 3;
    // RDKit✔️✔️:           ourHash |= (a1Deg % 8) << 6;
    // RDKit✔️✔️:           ourHash |= (a2Deg % 8) << 9;
    // RDKit✔️✔️:
    // RDKit✔️✔️:           hashLayers[1].push_back(ourHash);
    // RDKit✔️✔️:         }
    // END RDKIT CPP BLOCK LayeredFingerprintMol layer 2
    // Local complexity review: source and Rust each use one constant-time
    // query-mask branch, one bond-state branch, one degree canonicalization,
    // four modulo/packing fields, and no allocation or graph traversal.
    if path_queries & 0x1 != 0 {
        return None;
    }
    let bond_hash = if !bond.is_aromatic()
        && bond.order() != BondOrder::Single
        && bond.order() != BondOrder::Aromatic
    {
        bond.order().rdkit_code() as u32
    } else {
        BondOrder::Single.rdkit_code() as u32
    };
    let (larger_degree, smaller_degree) = if begin_atom_degree < end_atom_degree {
        (end_atom_degree, begin_atom_degree)
    } else {
        (begin_atom_degree, end_atom_degree)
    };
    Some(
        (bond_hash % 8)
            | ((bond_neighbor_count % 8) << 3)
            | ((larger_degree % 8) << 6)
            | ((smaller_degree % 8) << 9),
    )
}

#[inline]
fn layered_atom_type_hash(
    begin_atomic_number: u32,
    end_atomic_number: u32,
    begin_atom_degree: u32,
    end_atom_degree: u32,
    bond_neighbor_count: u32,
    path_queries: u8,
) -> Option<u32> {
    // BEGIN RDKIT CPP BLOCK LayeredFingerprintMol layer 3
    // RDKit✔️✔️:         if (layerFlags & 0x4 && !(pathQueries & 0x6)) {
    // RDKit✔️✔️:           // std::cerr<<" consider: "<<bi->getBeginAtomIdx()<<" - "
    // RDKit✔️✔️:           // <<bi->getEndAtomIdx()<<std::endl;
    // RDKit✔️✔️:           // layer 3: include atom types:
    // RDKit✔️✔️:           unsigned int a1Hash, a2Hash;
    // RDKit✔️✔️:           a1Hash = (anums[bi->getBeginAtomIdx()] % 128);
    // RDKit✔️✔️:           a2Hash = (anums[bi->getEndAtomIdx()] % 128);
    // RDKit✔️✔️:           unsigned int a1Deg, a2Deg;
    // RDKit✔️✔️:           a1Deg = atomDegrees[bi->getBeginAtomIdx()];
    // RDKit✔️✔️:           a2Deg = atomDegrees[bi->getEndAtomIdx()];
    // RDKit✔️✔️:           if (a1Hash < a2Hash) {
    // RDKit✔️✔️:             std::swap(a1Hash, a2Hash);
    // RDKit✔️✔️:             std::swap(a1Deg, a2Deg);
    // RDKit✔️✔️:           } else if (a1Hash == a2Hash && a1Deg < a2Deg) {
    // RDKit✔️✔️:             std::swap(a1Deg, a2Deg);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           ourHash = a1Hash;
    // RDKit✔️✔️:           ourHash |= a2Hash << 7;
    // RDKit✔️✔️:           ourHash |= (a1Deg % 8) << 14;
    // RDKit✔️✔️:           ourHash |= (a2Deg % 8) << 17;
    // RDKit✔️✔️:           ourHash |= (bondNbrs[i] % 8) << 20;
    // RDKit✔️✔️:           hashLayers[2].push_back(ourHash);
    // RDKit✔️✔️:         }
    // END RDKIT CPP BLOCK LayeredFingerprintMol layer 3
    // Local complexity review: source and Rust each use one suppression test,
    // two modulo-normalized atom keys, lexicographic endpoint ordering, and
    // five fixed-width packed fields in O(1), without allocation or lookup.
    if path_queries & 0x6 != 0 {
        return None;
    }
    let mut first_hash = begin_atomic_number % 128;
    let mut second_hash = end_atomic_number % 128;
    let mut first_degree = begin_atom_degree;
    let mut second_degree = end_atom_degree;
    if first_hash < second_hash {
        std::mem::swap(&mut first_hash, &mut second_hash);
        std::mem::swap(&mut first_degree, &mut second_degree);
    } else if first_hash == second_hash && first_degree < second_degree {
        std::mem::swap(&mut first_degree, &mut second_degree);
    }
    Some(
        first_hash
            | (second_hash << 7)
            | ((first_degree % 8) << 14)
            | ((second_degree % 8) << 17)
            | ((bond_neighbor_count % 8) << 20),
    )
}

#[inline]
fn layered_aromaticity_hash(
    begin_aromatic: bool,
    end_aromatic: bool,
    bond_neighbor_count: u32,
    path_queries: u8,
) -> Option<u32> {
    // BEGIN RDKIT CPP BLOCK LayeredFingerprintMol layer 6
    // RDKit✔️✔️:         if (layerFlags & 0x20 && !(pathQueries & 0x6)) {
    // RDKit✔️✔️:           // std::cerr<<" consider: "<<bi->getBeginAtomIdx()<<" - "
    // RDKit✔️✔️:           // <<bi->getEndAtomIdx()<<std::endl;
    // RDKit✔️✔️:           // layer 6: aromaticity:
    // RDKit✔️✔️:           bool a1Hash = aromaticAtoms[bi->getBeginAtomIdx()];
    // RDKit✔️✔️:           bool a2Hash = aromaticAtoms[bi->getEndAtomIdx()];
    // RDKit✔️✔️:
    // RDKit✔️✔️:           if ((!a1Hash) && a2Hash) {
    // RDKit✔️✔️:             std::swap(a1Hash, a2Hash);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           ourHash = a1Hash;
    // RDKit✔️✔️:           ourHash |= a2Hash << 1;
    // RDKit✔️✔️:           ourHash |= (bondNbrs[i] % 8) << 5;
    // RDKit✔️✔️:           hashLayers[5].push_back(ourHash);
    // RDKit✔️✔️:         }
    // END RDKIT CPP BLOCK LayeredFingerprintMol layer 6
    // Local complexity review: both implementations perform one suppression
    // branch, one endpoint canonicalization, one modulo, and three fixed-width
    // packs in O(1), with no allocation, clone, lookup, or traversal.
    if path_queries & 0x6 != 0 {
        return None;
    }
    let (first_aromatic, second_aromatic) = if !begin_aromatic && end_aromatic {
        (end_aromatic, begin_aromatic)
    } else {
        (begin_aromatic, end_aromatic)
    };
    Some(
        u32::from(first_aromatic)
            | (u32::from(second_aromatic) << 1)
            | ((bond_neighbor_count % 8) << 5),
    )
}

#[inline]
fn layered_ring_presence_hash(bond: &Bond, ring_info: &RingInfo, path_queries: u8) -> Option<u32> {
    // BEGIN RDKIT CPP BLOCK LayeredFingerprintMol layer 4
    // RDKit✔️✔️:         if (layerFlags & 0x8 && !(pathQueries & 0x6)) {
    // RDKit✔️✔️:           // layer 4: include ring information
    // RDKit✔️✔️:           if (queryIsBondInRing(bi)) {
    // RDKit✔️✔️:             hashLayers[3].push_back(1);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // END RDKIT CPP BLOCK LayeredFingerprintMol layer 4
    // Local complexity review: both forms perform the same mask branch and
    // O(1) indexed ring-membership lookup, allocate nothing, and deliberately
    // omit rather than encode non-ring bonds.
    if path_queries & 0x6 != 0 {
        return None;
    }
    (cosmolkit_search::query_is_bond_in_ring(bond, ring_info) != 0).then_some(1)
}

#[inline]
fn layered_min_ring_size_hash(bond: &Bond, ring_info: &RingInfo, path_queries: u8) -> Option<u32> {
    // BEGIN RDKIT CPP BLOCK LayeredFingerprintMol layer 5
    // RDKit✔️✔️:         if (layerFlags & 0x10 && !(pathQueries & 0x6)) {
    // RDKit✔️✔️:           // layer 5: include ring size information
    // RDKit✔️✔️:           ourHash = (queryBondMinRingSize(bi) % 8);
    // RDKit✔️✔️:           hashLayers[4].push_back(ourHash);
    // RDKit✔️✔️:         }
    // END RDKIT CPP BLOCK LayeredFingerprintMol layer 5
    // Local complexity review: source and Rust each scan only this bond's
    // ring-membership list to select the minimum and apply one modulo. Both
    // are O(R_bond), use O(1) auxiliary space, and allocate or clone nothing.
    if path_queries & 0x6 != 0 {
        return None;
    }
    Some((cosmolkit_search::query_bond_min_ring_size(bond, ring_info) % 8) as u32)
}

fn project_layered_path(
    hash_layers: &mut [Vec<u32>],
    atoms_in_path: &[bool],
    fp_size: usize,
    set_only_bits: Option<&Fingerprint>,
    result: &mut Fingerprint,
    mut atom_counts: Option<&mut [u32]>,
) -> Result<(), LayeredFingerprintError> {
    // BEGIN RDKIT CPP BLOCK LayeredFingerprintMol projection
    // RDKit✔️❌:       unsigned int l = 0;
    // RDKit✔️❌:       bool flaggedPath = false;
    // RDKit✔️❌:       for (auto layerIt = hashLayers.begin(); layerIt != hashLayers.end();
    // RDKit✔️❌:            ++layerIt, ++l) {
    // RDKit✔️❌:         if (!layerIt->size()) {
    // RDKit✔️❌:           continue;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         // ----
    // RDKit✔️❌:         std::sort(layerIt->begin(), layerIt->end());
    // RDKit✔️❌:
    // RDKit✔️❌:         // finally, we will add the number of distinct atoms in the path at the
    // RDKit✔️❌:         // end
    // RDKit✔️❌:         // of the vect. This allows us to distinguish C1CC1 from CC(C)C
    // RDKit✔️❌:         layerIt->push_back(static_cast<unsigned int>(atomsInPath.count()));
    // RDKit✔️❌:
    // RDKit✔️❌:         layerIt->push_back(l + 1);
    // RDKit✔️❌:
    // RDKit✔️❌:         // hash the path to generate a seed:
    // RDKit✔️❌:         unsigned long seed =
    // RDKit✔️❌:             gboost::hash_range(layerIt->begin(), layerIt->end());
    // RDKit✔️❌:
    // RDKit✔️❌:         unsigned int bitId = seed % fpSize;
    // RDKit✔️❌:         if (!setOnlyBits || (*setOnlyBits)[bitId]) {
    // RDKit✔️❌:           res->setBit(bitId);
    // RDKit✔️❌:           if (atomCounts && !flaggedPath) {
    // RDKit✔️❌:             for (unsigned int aIdx = 0; aIdx < atomsInPath.size(); ++aIdx) {
    // RDKit✔️❌:               if (atomsInPath[aIdx]) {
    // RDKit✔️❌:                 (*atomCounts)[aIdx] += 1;
    // RDKit✔️❌:               }
    // RDKit✔️❌:             }
    // RDKit✔️❌:             flaggedPath = true;
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // END RDKIT CPP BLOCK LayeredFingerprintMol projection
    // Local complexity review: each nonempty layer is sorted, extended by two
    // scalars, hashed once and projected into an indexed packed fingerprint.
    // The path atom mask is a byte-per-atom Vec<bool>, unlike the source packed
    // dynamic_bitset: counting scans O(A) bytes instead of packed words. Count
    // updates likewise inspect the byte mask. This is a material memory and
    // scan cost; the output fingerprint itself still uses packed words.
    debug_assert!(fp_size > 0);
    debug_assert_eq!(result.n_bits() as usize, fp_size);
    debug_assert!(set_only_bits.is_none_or(|bits| bits.n_bits() as usize == fp_size));
    debug_assert!(
        atom_counts
            .as_ref()
            .is_none_or(|counts| counts.len() >= atoms_in_path.len())
    );

    let distinct_atom_count = atoms_in_path.iter().filter(|&&present| present).count() as u32;
    let mut flagged_path = false;
    for (layer_index, layer) in hash_layers.iter_mut().enumerate() {
        if layer.is_empty() {
            continue;
        }
        layer.sort_unstable();
        layer.push(distinct_atom_count);
        layer.push(layer_index as u32 + 1);

        let bit_id = hash_range(layer) as usize % fp_size;
        let accepted = match set_only_bits {
            Some(mask) => mask.get_bit(bit_id as u32)?,
            None => true,
        };
        if !accepted {
            continue;
        }
        result.set_bit(bit_id as u32)?;
        if !flagged_path {
            if let Some(counts) = atom_counts.as_deref_mut() {
                for atom_index in 0..atoms_in_path.len() {
                    if atoms_in_path[atom_index] {
                        counts[atom_index] = counts[atom_index].wrapping_add(1);
                    }
                }
                flagged_path = true;
            }
        }
    }
    Ok(())
}

pub fn layered_fingerprint(
    topology: &TopologyBlock,
    cached_rings: Option<&RingInfo>,
    params: &LayeredFingerprintParams,
) -> Result<Fingerprint, LayeredFingerprintError> {
    Ok(layered_fingerprint_with_output(topology, cached_rings, params)?.fingerprint)
}

/// Compute a source-backed Layered fingerprint and optional atom counts.
///
/// This read-only operation neither mutates the molecule nor stores ring or
/// fingerprint intermediates in it. Calls can be repeated, interleaved with
/// other fingerprint families, or run concurrently on shared molecules.
/// Invalid path bounds, widths, count lengths, mask widths, and roots return a
/// structured [`LayeredFingerprintError`].
pub fn layered_fingerprint_with_output(
    topology: &TopologyBlock,
    cached_rings: Option<&RingInfo>,
    params: &LayeredFingerprintParams,
) -> Result<LayeredFingerprintResult, LayeredFingerprintError> {
    layered_fingerprint_impl(LayeredGraphInput::Concrete(topology), cached_rings, params)
}

pub fn layered_query_fingerprint_with_output(
    query: &QueryGraph,
    cached_rings: Option<&RingInfo>,
    params: &LayeredFingerprintParams,
) -> Result<LayeredFingerprintResult, LayeredFingerprintError> {
    layered_fingerprint_impl(LayeredGraphInput::Query(query), cached_rings, params)
}

pub fn layered_query_fingerprint(
    query: &QueryGraph,
    cached_rings: Option<&RingInfo>,
    params: &LayeredFingerprintParams,
) -> Result<Fingerprint, LayeredFingerprintError> {
    Ok(layered_query_fingerprint_with_output(query, cached_rings, params)?.fingerprint)
}

fn layered_fingerprint_impl(
    graph: LayeredGraphInput<'_>,
    cached_rings: Option<&RingInfo>,
    params: &LayeredFingerprintParams,
) -> Result<LayeredFingerprintResult, LayeredFingerprintError> {
    params.validate()?;
    if params.from_atoms.as_ref().is_some_and(|roots| {
        roots
            .iter()
            .any(|&root| root as usize >= graph.atom_count())
    }) {
        return Err(LayeredFingerprintError::InvalidArguments {
            reason: "fromAtoms contains atom index out of range",
        });
    }

    // Preserve the caller seed and return detached updated counts: cloning
    // allocates and copies O(seed length), including any extra seed rows.
    // The source updates its caller-owned vector without this extra copy.
    let mut atom_counts = params.atom_counts.clone();
    let prepared = prepare_layered_fingerprint(
        graph,
        cached_rings,
        params.min_path,
        params.max_path,
        params.fp_size as usize,
        atom_counts.as_deref(),
        params.set_only_bits.as_ref(),
    )?;

    // BEGIN RDKIT CPP BLOCK LayeredFingerprintMol path selection
    // RDKit✔️✔️:   auto *res = new ExplicitBitVect(fpSize);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   INT_PATH_LIST_MAP allPaths;
    // RDKit✔️✔️:   if (!fromAtoms) {
    // RDKit✔️✔️:     if (branchedPaths) {
    // RDKit✔️✔️:       allPaths = findAllSubgraphsOfLengthsMtoN(mol, minPath, maxPath, false);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       allPaths = findAllPathsOfLengthsMtoN(mol, minPath, maxPath, false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     for (auto aidx : *fromAtoms) {
    // RDKit✔️✔️:       INT_PATH_LIST_MAP tPaths;
    // RDKit✔️✔️:       if (branchedPaths) {
    // RDKit✔️✔️:         tPaths =
    // RDKit✔️✔️:             findAllSubgraphsOfLengthsMtoN(mol, minPath, maxPath, false, aidx);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         tPaths =
    // RDKit✔️✔️:             findAllPathsOfLengthsMtoN(mol, minPath, maxPath, true, false, aidx);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       for (INT_PATH_LIST_MAP::const_iterator tpit = tPaths.begin();
    // RDKit✔️✔️:            tpit != tPaths.end(); ++tpit) {
    // RDKit✔️✔️:         allPaths[tpit->first].insert(allPaths[tpit->first].begin(),
    // RDKit✔️✔️:                                      tpit->second.begin(), tpit->second.end());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP BLOCK LayeredFingerprintMol path selection
    // Preserve the pinned source's atom IDs, atom-length limits and ordering
    // for this branch, including their later interpretation as bond IDs.
    // Do not replace it with the rooted branch's corrected bond paths.
    // Defined native accesses reproduce literally; the existing checked bond
    // lookup below rejects paths which would access C++ storage out of bounds.
    // Local cost: reuse CORE's same enumerator and move each ID vector into
    // the existing per-length map; no graph clone or second enumeration.
    let all_paths = if !params.branched_paths && params.from_atoms.is_none() {
        let paths = match graph {
            LayeredGraphInput::Concrete(topology) => all_paths_in_range(
                topology,
                params.min_path as usize,
                params.max_path as usize,
                &PathSearchParams {
                    representation: PathRepresentation::Atoms,
                    ..Default::default()
                },
            )?,
            LayeredGraphInput::Query(query) => query_atom_paths_in_range(
                query,
                params.min_path as usize,
                params.max_path as usize,
                &SubgraphSearchParams::default(),
            )?,
        };
        paths
            .into_iter()
            .map(|(size, rows)| {
                let rows = rows
                    .into_iter()
                    .map(|row| match row {
                        GraphPath::Atoms(ids) => Ok(ids.into_iter().map(|id| id.index()).collect()),
                        GraphPath::Bonds(_) => Err(LayeredFingerprintError::InvalidArguments {
                            reason: "atom path enumeration returned bond path",
                        }),
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                Ok((size, rows))
            })
            .collect::<Result<BTreeMap<_, _>, LayeredFingerprintError>>()?
    } else {
        enumerate_fingerprint_paths(
            graph,
            params.min_path,
            params.max_path,
            false,
            params.branched_paths,
            params.from_atoms.as_deref(),
            None,
        )?
    };
    let fp_size = params.fp_size as usize;
    let mut fingerprint = Fingerprint::new(params.fp_size);

    for paths in all_paths.values() {
        for path in paths {
            // BEGIN RDKIT CPP BLOCK LayeredFingerprintMol path intermediates
            // RDKit✔️❌:       std::vector<std::vector<unsigned int>> hashLayers(maxFingerprintLayers);
            // RDKit✔️❌:       for (unsigned int i = 0; i < maxFingerprintLayers; ++i) {
            // RDKit✔️❌:         if (layerFlags & (0x1 << i)) {
            // RDKit✔️❌:           hashLayers[i].reserve(maxPath);
            // RDKit✔️❌:         }
            // RDKit✔️❌:       }
            // RDKit✔️❌:
            // RDKit✔️❌:       // details about what kinds of query features appear on the path:
            // RDKit✔️❌:       unsigned int pathQueries = 0;
            // RDKit✔️❌:       for (int pIt : path) {
            // RDKit✔️❌:         pathQueries |= isQueryBond[pIt];
            // RDKit✔️❌:       }
            // RDKit✔️❌:
            // RDKit✔️❌:       // calculate the number of neighbors each bond has in the path:
            // RDKit✔️❌:       std::vector<unsigned int> bondNbrs(path.size(), 0);
            // RDKit✔️❌:       atomsInPath.reset();
            // RDKit✔️❌:
            // RDKit✔️❌:       std::vector<unsigned int> atomDegrees(mol.getNumAtoms(), 0);
            // RDKit✔️❌:       for (int i : path) {
            // RDKit✔️❌:         const Bond *bi = bondCache[i];
            // RDKit✔️❌:         atomDegrees[bi->getBeginAtomIdx()]++;
            // RDKit✔️❌:         atomDegrees[bi->getEndAtomIdx()]++;
            // RDKit✔️❌:         atomsInPath.set(bi->getBeginAtomIdx());
            // RDKit✔️❌:         atomsInPath.set(bi->getEndAtomIdx());
            // RDKit✔️❌:       }
            // RDKit✔️❌:
            // RDKit✔️❌:       for (unsigned int i = 0; i < path.size(); ++i) {
            // RDKit✔️❌:         const Bond *bi = bondCache[path[i]];
            // RDKit✔️❌:         for (unsigned int j = i + 1; j < path.size(); ++j) {
            // RDKit✔️❌:           const Bond *bj = bondCache[path[j]];
            // RDKit✔️❌:           if (bi->getBeginAtomIdx() == bj->getBeginAtomIdx() ||
            // RDKit✔️❌:               bi->getBeginAtomIdx() == bj->getEndAtomIdx() ||
            // RDKit✔️❌:               bi->getEndAtomIdx() == bj->getBeginAtomIdx() ||
            // RDKit✔️❌:               bi->getEndAtomIdx() == bj->getEndAtomIdx()) {
            // RDKit✔️❌:             ++bondNbrs[i];
            // RDKit✔️❌:             ++bondNbrs[j];
            // RDKit✔️❌:           }
            // RDKit✔️❌:         }
            // END RDKIT CPP BLOCK LayeredFingerprintMol path intermediates
            // Local complexity review: source and Rust both allocate ten
            // layer vectors, scan P bonds for query/degree state, perform the
            // same O(P^2) pairwise adjacency count, and encode each bond once.
            // Rust additionally allocates a new O(A) byte mask for each path;
            // the source allocates one packed dynamic_bitset before all paths
            // and resets/reuses it. Per-path allocation and byte-mask memory
            // costs prevent performance equivalence for this block.
            let mut hash_layers = vec![Vec::new(); LAYERED_FINGERPRINT_MAX_LAYERS];
            for (layer_index, layer) in hash_layers.iter_mut().enumerate() {
                if params.layers.bits() & (1u32 << layer_index) != 0 {
                    layer.reserve(params.max_path as usize);
                }
            }

            let mut path_queries = 0u8;
            let mut atom_degrees = vec![0u32; graph.atom_count()];
            let mut atoms_in_path = vec![false; graph.atom_count()];
            for &bond_index in path {
                let bond = *prepared.bond_cache.get(bond_index).ok_or(
                    LayeredFingerprintError::InvalidArguments {
                        reason: "enumerated path contains invalid bond index",
                    },
                )?;
                path_queries |= prepared.query_masks[bond_index];
                let begin = bond.begin().index();
                let end = bond.end().index();
                atom_degrees[begin] = atom_degrees[begin].wrapping_add(1);
                atom_degrees[end] = atom_degrees[end].wrapping_add(1);
                atoms_in_path[begin] = true;
                atoms_in_path[end] = true;
            }

            let mut bond_neighbors = vec![0u32; path.len()];
            for first_position in 0..path.len() {
                let first = prepared.bond_cache[path[first_position]];
                for second_position in (first_position + 1)..path.len() {
                    let second = prepared.bond_cache[path[second_position]];
                    if first.begin() == second.begin()
                        || first.begin() == second.end()
                        || first.end() == second.begin()
                        || first.end() == second.end()
                    {
                        bond_neighbors[first_position] =
                            bond_neighbors[first_position].wrapping_add(1);
                        bond_neighbors[second_position] =
                            bond_neighbors[second_position].wrapping_add(1);
                    }
                }
            }

            for (path_position, &bond_index) in path.iter().enumerate() {
                let bond = prepared.bond_cache[bond_index];
                let begin = bond.begin().index();
                let end = bond.end().index();
                let neighbor_count = bond_neighbors[path_position];
                let begin_degree = atom_degrees[begin];
                let end_degree = atom_degrees[end];

                if params.layers.contains(LayeredFingerprintLayers::TOPOLOGY) {
                    hash_layers[0].push(layered_topology_hash(
                        neighbor_count,
                        begin_degree,
                        end_degree,
                    ));
                }
                if params.layers.contains(LayeredFingerprintLayers::BOND_ORDER) {
                    if let Some(hash) = layered_bond_order_hash(
                        bond,
                        neighbor_count,
                        begin_degree,
                        end_degree,
                        path_queries,
                    ) {
                        hash_layers[1].push(hash);
                    }
                }
                if params.layers.contains(LayeredFingerprintLayers::ATOM_TYPE) {
                    if let Some(hash) = layered_atom_type_hash(
                        prepared.atomic_numbers[begin],
                        prepared.atomic_numbers[end],
                        begin_degree,
                        end_degree,
                        neighbor_count,
                        path_queries,
                    ) {
                        hash_layers[2].push(hash);
                    }
                }
                if params
                    .layers
                    .contains(LayeredFingerprintLayers::RING_PRESENCE)
                {
                    if let Some(hash) =
                        layered_ring_presence_hash(bond, prepared.ring_info.as_ref(), path_queries)
                    {
                        hash_layers[3].push(hash);
                    }
                }
                if params.layers.contains(LayeredFingerprintLayers::RING_SIZE) {
                    if let Some(hash) =
                        layered_min_ring_size_hash(bond, prepared.ring_info.as_ref(), path_queries)
                    {
                        hash_layers[4].push(hash);
                    }
                }
                if params
                    .layers
                    .contains(LayeredFingerprintLayers::AROMATICITY)
                {
                    if let Some(hash) = layered_aromaticity_hash(
                        prepared.aromatic_atoms[begin],
                        prepared.aromatic_atoms[end],
                        neighbor_count,
                        path_queries,
                    ) {
                        hash_layers[5].push(hash);
                    }
                }
            }

            project_layered_path(
                &mut hash_layers,
                &atoms_in_path,
                fp_size,
                params.set_only_bits.as_ref(),
                &mut fingerprint,
                atom_counts.as_deref_mut(),
            )?;
        }
    }

    Ok(LayeredFingerprintResult {
        fingerprint,
        atom_counts,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cco() -> TopologyBlock {
        cosmolkit_smiles::parse_smiles("CCO", &Default::default())
            .unwrap()
            .topology
    }
    #[test]
    fn original_default_all_bits_and_seeded_counts_are_preserved() {
        let topology = cco();
        let original = topology.clone();
        let params = LayeredFingerprintParams {
            atom_counts: Some(vec![0; 3]),
            ..Default::default()
        };
        let result = layered_fingerprint_with_output(&topology, None, &params).unwrap();
        assert_eq!(result.fingerprint.n_bits(), 2048);
        assert_eq!(
            result.fingerprint.on_bits(),
            [92, 360, 596, 610, 611, 674, 867, 1044, 1111, 1783, 1784]
        );
        assert_eq!(result.atom_counts, Some(vec![2, 3, 2]));
        assert_eq!(topology, original);
        assert_eq!(params.atom_counts, Some(vec![0; 3]));
        assert_eq!(
            layered_fingerprint(
                &topology,
                None,
                &LayeredFingerprintParams {
                    layers: LayeredFingerprintLayers::ACTIVE,
                    ..params
                }
            )
            .unwrap(),
            result.fingerprint
        );
    }
    #[test]
    fn original_masks_empty_roots_and_rooted_linear_results_are_preserved() {
        let topology = cco();
        let masked = LayeredFingerprintParams {
            atom_counts: Some(vec![10, 20, 30]),
            set_only_bits: Some(Fingerprint::from_on_bits(2048, [674]).unwrap()),
            ..Default::default()
        };
        let result = layered_fingerprint_with_output(&topology, None, &masked).unwrap();
        assert_eq!(result.fingerprint.on_bits(), [674]);
        assert_eq!(result.atom_counts, Some(vec![11, 22, 31]));
        let empty = LayeredFingerprintParams {
            from_atoms: Some(vec![]),
            atom_counts: Some(vec![0; 3]),
            ..Default::default()
        };
        let result = layered_fingerprint_with_output(&topology, None, &empty).unwrap();
        assert!(result.fingerprint.on_bits().is_empty());
        assert_eq!(result.atom_counts, Some(vec![0; 3]));
        let linear = LayeredFingerprintParams {
            branched_paths: false,
            from_atoms: Some(vec![0]),
            atom_counts: Some(vec![0; 3]),
            ..Default::default()
        };
        let result = layered_fingerprint_with_output(&topology, None, &linear).unwrap();
        assert_eq!(
            result.fingerprint.on_bits(),
            [360, 596, 610, 611, 674, 867, 1044, 1111, 1783, 1784]
        );
        assert_eq!(result.atom_counts, Some(vec![2, 2, 1]));
    }

    #[test]
    fn unrooted_linear_uses_pinned_atom_paths_for_concrete_and_query_graphs() {
        // Pinned RDKit 2026.03.6 LayeredFingerprint, not CK-generated values:
        // phenol, layerFlags=0xffffffff, minPath=2, maxPath=4, fpSize=4096,
        // branchedPaths=False, zero atomCounts, even setOnlyBits.
        let expected_bits = [
            138, 170, 354, 366, 470, 978, 1590, 1672, 1750, 1784, 2016, 2096, 2152, 2366, 2610,
            3758, 3870, 3926,
        ];
        let expected_counts = vec![5, 14, 14, 11, 10, 10, 11];
        let topology = cosmolkit_smiles::parse_smiles("Oc1ccccc1", &Default::default())
            .unwrap()
            .topology;
        let query = cosmolkit_search::parse_smarts("Oc1ccccc1", &Default::default()).unwrap();
        let params = LayeredFingerprintParams {
            min_path: 2,
            max_path: 4,
            fp_size: 4096,
            branched_paths: false,
            layers: LayeredFingerprintLayers::ALL_SOURCE_BITS,
            atom_counts: Some(vec![0; 7]),
            set_only_bits: Some(Fingerprint::from_on_bits(4096, (0..4096).step_by(2)).unwrap()),
            ..Default::default()
        };
        let original = topology.clone();
        let query_before = format!("{query:?}");
        for graph in [
            LayeredGraphInput::Concrete(&topology),
            LayeredGraphInput::Query(&query),
        ] {
            let result = layered_fingerprint_impl(graph, None, &params).unwrap();
            assert_eq!(result.fingerprint.on_bits(), expected_bits);
            assert_eq!(result.atom_counts.as_ref(), Some(&expected_counts));
        }
        assert_eq!(topology, original);
        assert_eq!(format!("{query:?}"), query_before);
        assert_eq!(params.atom_counts, Some(vec![0; 7]));
    }

    #[test]
    fn unrooted_linear_native_invalid_bond_access_is_a_checked_error() {
        let topology = cco();
        let original = topology.clone();
        let params = LayeredFingerprintParams {
            branched_paths: false,
            ..Default::default()
        };
        assert!(matches!(
            layered_fingerprint(&topology, None, &params),
            Err(LayeredFingerprintError::InvalidArguments {
                reason: "enumerated path contains invalid bond index"
            })
        ));
        assert_eq!(topology, original);
        let query = cosmolkit_search::parse_smarts("CCO", &Default::default()).unwrap();
        assert!(matches!(
            layered_query_fingerprint(&query, None, &params),
            Err(LayeredFingerprintError::InvalidArguments {
                reason: "enumerated path contains invalid bond index"
            })
        ));
    }
    #[test]
    fn counts_wrap_once_per_path_and_keep_extra_seed_entries() {
        let topology = cco();
        let params = LayeredFingerprintParams {
            fp_size: 1,
            atom_counts: Some(vec![u32::MAX, u32::MAX, u32::MAX, 55]),
            ..Default::default()
        };
        let result = layered_fingerprint_with_output(&topology, None, &params).unwrap();
        assert_eq!(result.fingerprint.on_bits(), [0]);
        assert_eq!(result.atom_counts, Some(vec![1, 2, 1, 55]));
        let duplicated = LayeredFingerprintParams {
            from_atoms: Some(vec![0, 0]),
            atom_counts: Some(vec![0; 3]),
            ..Default::default()
        };
        assert_eq!(
            layered_fingerprint_with_output(&topology, None, &duplicated)
                .unwrap()
                .atom_counts,
            Some(vec![4, 4, 2])
        );
    }
    #[test]
    fn original_high_flags_and_all_preconditions_are_preserved() {
        let topology = cco();
        let high = LayeredFingerprintParams {
            layers: LayeredFingerprintLayers::from_bits_retain(0xffff_ffc0),
            atom_counts: Some(vec![5, 6, 7]),
            ..Default::default()
        };
        let result = layered_fingerprint_with_output(&topology, None, &high).unwrap();
        assert!(result.fingerprint.on_bits().is_empty());
        assert_eq!(result.atom_counts, Some(vec![5, 6, 7]));
        for (params, reason) in [
            (
                LayeredFingerprintParams {
                    min_path: 0,
                    ..Default::default()
                },
                "minPath==0",
            ),
            (
                LayeredFingerprintParams {
                    min_path: 3,
                    max_path: 2,
                    ..Default::default()
                },
                "maxPath<minPath",
            ),
            (
                LayeredFingerprintParams {
                    fp_size: 0,
                    ..Default::default()
                },
                "fpSize==0",
            ),
            (
                LayeredFingerprintParams {
                    atom_counts: Some(vec![0; 2]),
                    ..Default::default()
                },
                "bad atomCounts size",
            ),
            (
                LayeredFingerprintParams {
                    set_only_bits: Some(Fingerprint::new(64)),
                    ..Default::default()
                },
                "bad setOnlyBits size",
            ),
            (
                LayeredFingerprintParams {
                    from_atoms: Some(vec![3]),
                    ..Default::default()
                },
                "fromAtoms contains atom index out of range",
            ),
        ] {
            assert!(matches!(layered_fingerprint(&topology, None, &params),
                Err(LayeredFingerprintError::InvalidArguments { reason: actual }) if actual == reason));
        }
    }
    #[test]
    fn original_query22_fixture_preserves_masks_aromaticity_and_complete_bits() {
        let fixture: serde_json::Value = serde_json::from_str(include_str!(
            "../../../testdata/fingerprint/fixtures/rdkit/layered_fingerprint_query_cases.json"
        ))
        .unwrap();
        let params = LayeredFingerprintParams {
            layers: LayeredFingerprintLayers::ACTIVE,
            fp_size: 512,
            ..Default::default()
        };
        let mut count = 0;
        for case in fixture["complexity_masks"]
            .as_array()
            .unwrap()
            .iter()
            .chain(fixture["aromaticity_branches"].as_array().unwrap())
        {
            let id = case["case_id"].as_str().unwrap();
            let input = case["input"].as_str().unwrap();
            let masks: Vec<u8> = serde_json::from_value(case["query_masks"].clone()).unwrap();
            let aromatic: Vec<bool> =
                serde_json::from_value(case["aromatic_atoms"].clone()).unwrap();
            let bits: Vec<u32> = serde_json::from_value(case["on_bits"].clone()).unwrap();
            if case["notation"] == "smarts" {
                let query = cosmolkit_search::parse_smarts(input, &Default::default()).unwrap();
                // RecursiveStructureQuery::copy quick-copies nested molecules,
                // so Clone is not a lossless storage snapshot (ROMol.cpp).
                let original = format!("{query:?}");
                let prepared = prepare_layered_fingerprint(
                    LayeredGraphInput::Query(&query),
                    None,
                    1,
                    7,
                    512,
                    None,
                    None,
                )
                .unwrap();
                assert_eq!(prepared.query_masks, masks, "{id}");
                assert_eq!(prepared.aromatic_atoms, aromatic, "{id}");
                assert_eq!(
                    layered_query_fingerprint(&query, None, &params)
                        .unwrap()
                        .on_bits(),
                    bits,
                    "{id}"
                );
                assert_eq!(format!("{query:?}"), original);
            } else {
                let topology = cosmolkit_smiles::parse_smiles(input, &Default::default())
                    .unwrap()
                    .topology;
                let prepared = prepare_layered_fingerprint(
                    LayeredGraphInput::Concrete(&topology),
                    None,
                    1,
                    7,
                    512,
                    None,
                    None,
                )
                .unwrap();
                assert_eq!(prepared.query_masks, masks, "{id}");
                assert_eq!(prepared.aromatic_atoms, aromatic, "{id}");
                assert_eq!(
                    layered_fingerprint(&topology, None, &params)
                        .unwrap()
                        .on_bits(),
                    bits,
                    "{id}"
                );
            }
            count += 1;
        }
        assert_eq!(count, 21);
    }
    #[test]
    fn original_additional_typed_aromaticity_assertions_are_preserved() {
        use cosmolkit_model::{
            Atom, AtomQueryPredicate as A, AtomSpec, Element, QueryAtom, QueryNode as Q,
        };
        let number = || Q::predicate(A::AtomicNumber(6));
        let aromatic = || Q::predicate(A::IsAromatic(true));
        let charge = || Q::predicate(A::FormalCharge(0));
        for (predicate, expected) in [
            (number(), true),
            (Q::or(vec![aromatic()]), false),
            (Q::xor(vec![aromatic(), number()]), false),
            (Q::and(vec![aromatic(), number()]), false),
            (Q::and(vec![number(), charge()]), false),
            (Q::not(Q::and(vec![number(), aromatic()])), false),
            (Q::not(Q::not(aromatic())), false),
            (charge(), false),
        ] {
            let atom = QueryAtom::from_parts(
                Atom::from_spec(
                    AtomId::new(0),
                    AtomSpec::new(Element::C).with_aromatic(true),
                ),
                predicate,
            );
            let graph = QueryGraph::from_parts(
                vec![atom],
                vec![],
                Vec::<(
                    cosmolkit_model::PropertyText,
                    cosmolkit_model::PropertyValue,
                )>::new(),
                vec![],
                vec![],
                vec![],
            )
            .unwrap();
            assert_eq!(is_query_atom_aromatic(&graph.atoms()[0], &graph), expected);
        }
    }
    #[test]
    fn concrete_source_null_query_suppresses_only_bond_order_layer() {
        use cosmolkit_model::{
            AtomSpec, Bond, BondId, BondOrder, BondQueryPredicate, BondSpec, Element, QueryNode,
        };
        let atoms = (0..2)
            .map(|i| cosmolkit_model::Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect::<Vec<_>>();
        let spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Unspecified);
        let ordinary = TopologyBlock::try_from_parts(
            atoms.clone(),
            vec![Bond::from_spec(BondId::new(0), spec.clone())],
            vec![],
            vec![],
        )
        .unwrap();
        let explicit = TopologyBlock::try_from_parts(
            atoms,
            vec![Bond::from_spec(
                BondId::new(0),
                spec.with_query(QueryNode::predicate(BondQueryPredicate::Any)),
            )],
            vec![],
            vec![],
        )
        .unwrap();
        let params = LayeredFingerprintParams {
            layers: LayeredFingerprintLayers::BOND_ORDER,
            ..Default::default()
        };
        assert_eq!(
            layered_fingerprint(&ordinary, None, &params)
                .unwrap()
                .on_bits(),
            vec![507]
        );
        assert!(
            layered_fingerprint(&explicit, None, &params)
                .unwrap()
                .on_bits()
                .is_empty()
        );
        let params = LayeredFingerprintParams {
            layers: LayeredFingerprintLayers::ACTIVE,
            ..Default::default()
        };
        assert_eq!(
            layered_fingerprint(&explicit, None, &params)
                .unwrap()
                .on_bits(),
            vec![610, 611, 674, 1044]
        );
    }
}
