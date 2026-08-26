use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashSet, VecDeque};
use std::ops::{BitOr, BitOrAssign};
use std::sync::OnceLock;

use crate::chemistry::ciplabeler::assign_cip_labels;
use crate::chemistry::valence::rdkit_atomic_mass;
use crate::search::smarts_parse::{SmartsParseParams, mol_from_smarts};
use crate::{AdjacencyList, AtomId, BondOrder, ChiralTag, Molecule};
use serde_json::Value;

mod atom_pair;
pub(crate) mod generator;

pub(crate) use atom_pair::atom_pair_function_arguments;
pub use atom_pair::{
    AtomPairAtomInvariantsGenerator as AtomPairAtomInvGenerator, AtomPairFingerprintGenerator,
    AtomPairFingerprintOutput, AtomPairFingerprintParams, atom_pair_fingerprint,
    atom_pair_fingerprint_with_output,
};

// RDKit marker convention defined in dev/source_reproduction_protocol.md.
// Copied source lines appear as:  // RDKit<beh><perf>: ...

// RDKit source file: FingerprintUtil.cpp

// BEGIN RDKIT CPP CONSTANTS AtomPairs atom-code layout (FingerprintUtil.h)
// RDKit✔️✔️: const unsigned int numTypeBits = 4;
pub const ATOM_PAIR_NUM_TYPE_BITS: u32 = atom_pair::NUM_TYPE_BITS;
// RDKit✔️✔️: const unsigned int atomNumberTypes[1 << numTypeBits] = {
// RDKit✔️✔️:     5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53};
pub const ATOM_PAIR_ATOM_NUMBER_TYPES: [u32; 1 << ATOM_PAIR_NUM_TYPE_BITS] =
    atom_pair::ATOM_NUMBER_TYPES;
// RDKit✔️✔️: const unsigned int numPiBits = 2;
pub const ATOM_PAIR_NUM_PI_BITS: u32 = atom_pair::NUM_PI_BITS;
// RDKit✔️✔️: const unsigned int maxNumPi = (1 << numPiBits) - 1;
pub const ATOM_PAIR_MAX_NUM_PI: u32 = atom_pair::MAX_NUM_PI;
// RDKit✔️✔️: const unsigned int numBranchBits = 3;
pub const ATOM_PAIR_NUM_BRANCH_BITS: u32 = atom_pair::NUM_BRANCH_BITS;
// RDKit✔️✔️: const unsigned int maxNumBranches = (1 << numBranchBits) - 1;
pub const ATOM_PAIR_MAX_NUM_BRANCHES: u32 = atom_pair::MAX_NUM_BRANCHES;
// RDKit✔️✔️: const unsigned int numChiralBits = 2;
pub const ATOM_PAIR_NUM_CHIRAL_BITS: u32 = atom_pair::NUM_CHIRAL_BITS;
// RDKit✔️✔️: const unsigned int codeSize = numTypeBits + numPiBits + numBranchBits;
pub const ATOM_PAIR_CODE_SIZE: u32 = atom_pair::CODE_SIZE;
// RDKit✔️✔️: const unsigned int numPathBits = 5;
pub const ATOM_PAIR_NUM_PATH_BITS: u32 = atom_pair::NUM_PATH_BITS;
// RDKit✔️✔️: const unsigned int maxPathLen = (1 << numPathBits) - 1;
pub const ATOM_PAIR_MAX_PATH_LENGTH: u32 = atom_pair::MAX_PATH_LEN;
// RDKit✔️✔️: const unsigned int numAtomPairFingerprintBits =
// RDKit✔️✔️:     numPathBits + 2 * codeSize;
pub const ATOM_PAIR_NUM_FINGERPRINT_BITS: u32 = atom_pair::NUM_ATOM_PAIR_FINGERPRINT_BITS;
// RDKit source: AtomPairs.h line 44.
// RDKit✔️✔️: const std::string atomPairsVersion = "1.1.0";
pub const ATOM_PAIRS_VERSION: &str = "1.1.0";
// END RDKIT CPP CONSTANTS AtomPairs atom-code layout (FingerprintUtil.h)

pub fn get_atom_code(
    molecule: &Molecule,
    atom_id: AtomId,
    branch_subtract: u32,
    include_chirality: bool,
) -> Result<u32, FingerprintError> {
    atom_pair::get_atom_code(molecule, atom_id, branch_subtract, include_chirality)
        .map_err(FingerprintError::from)
}

#[must_use]
const fn topological_torsion_correct_atom_invariant(invariant: u32) -> u32 {
    invariant.wrapping_sub(2)
}

#[must_use]
fn topological_torsion_reverse(path_codes: &[u32]) -> bool {
    let mut i = 0usize;
    let mut j = path_codes.len() - 1;
    while i < j {
        if path_codes[i] > path_codes[j] {
            return true;
        } else if path_codes[i] < path_codes[j] {
            return false;
        }
        i += 1;
        j -= 1;
    }
    false
}

pub fn get_topological_torsion_code(
    path_codes: &[u32],
    include_chirality: bool,
) -> Result<u64, FingerprintError> {
    // BEGIN RDKIT CPP FUNCTION AtomPairs::getTopologicalTorsionCode
    // RDKit✔️✔️: std::uint64_t getTopologicalTorsionCode(
    // RDKit✔️✔️:     const std::vector<std::uint32_t> &pathCodes, bool includeChirality) {
    // RDKit✔️✔️:   bool reverseIt = false;
    // RDKit✔️✔️:   unsigned int i = 0;
    // RDKit✔️✔️:   unsigned int j = pathCodes.size() - 1;
    // RDKit✔️✔️:   while (i < j) {
    // RDKit✔️✔️:     if (pathCodes[i] > pathCodes[j]) {
    // RDKit✔️✔️:       reverseIt = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     } else if (pathCodes[i] < pathCodes[j]) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++i;
    // RDKit✔️✔️:     --j;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int shiftSize = codeSize + (includeChirality ? numChiralBits : 0);
    // RDKit✔️✔️:   std::uint64_t res = 0;
    // RDKit✔️✔️:   if (reverseIt) {
    // RDKit✔️✔️:     for (unsigned int i = 0; i < pathCodes.size(); ++i) {
    // RDKit✔️✔️:       res |= static_cast<std::uint64_t>(pathCodes[pathCodes.size() - i - 1])
    // RDKit✔️✔️:              << (shiftSize * i);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     for (unsigned int i = 0; i < pathCodes.size(); ++i) {
    // RDKit✔️✔️:       res |= static_cast<std::uint64_t>(pathCodes[i]) << (shiftSize * i);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AtomPairs::getTopologicalTorsionCode
    // The source underflows `size() - 1` for an empty vector and has undefined
    // behavior when a left-shift count reaches 64. COSMolKit fails closed for
    // those inputs while preserving the source's exact u64 packing otherwise.
    if path_codes.is_empty() {
        return Err(FingerprintError::InvalidArguments {
            reason: "topological torsion code path must not be empty",
        });
    }

    let shift_size = (ATOM_PAIR_CODE_SIZE
        + if include_chirality {
            ATOM_PAIR_NUM_CHIRAL_BITS
        } else {
            0
        }) as usize;
    let last_shift =
        shift_size
            .checked_mul(path_codes.len() - 1)
            .ok_or(FingerprintError::InvalidArguments {
                reason: "topological torsion code shift exceeds platform limits",
            })?;
    if last_shift >= u64::BITS as usize {
        return Err(FingerprintError::InvalidArguments {
            reason: "topological torsion code shift must be less than 64 bits",
        });
    }

    let reverse_it = topological_torsion_reverse(path_codes);
    let mut result = 0u64;
    if reverse_it {
        for (i, path_code) in path_codes.iter().rev().enumerate() {
            result |= u64::from(*path_code) << (shift_size * i);
        }
    } else {
        for (i, path_code) in path_codes.iter().enumerate() {
            result |= u64::from(*path_code) << (shift_size * i);
        }
    }
    Ok(result)
}

pub fn get_topological_torsion_hash(path_codes: &[u32]) -> Result<u32, FingerprintError> {
    // BEGIN RDKIT CPP FUNCTION AtomPairs::getTopologicalTorsionHash
    // RDKit✔️✔️: std::uint32_t getTopologicalTorsionHash(
    // RDKit✔️✔️:     const std::vector<std::uint32_t> &pathCodes) {
    // RDKit✔️✔️:   bool reverseIt = false;
    // RDKit✔️✔️:   unsigned int i = 0;
    // RDKit✔️✔️:   unsigned int j = pathCodes.size() - 1;
    // RDKit✔️✔️:   while (i < j) {
    // RDKit✔️✔️:     if (pathCodes[i] > pathCodes[j]) {
    // RDKit✔️✔️:       reverseIt = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     } else if (pathCodes[i] < pathCodes[j]) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++i;
    // RDKit✔️✔️:     --j;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::uint32_t res = 0;
    // RDKit✔️✔️:   if (reverseIt) {
    // RDKit✔️✔️:     for (unsigned int i = 0; i < pathCodes.size(); ++i) {
    // RDKit✔️✔️:       gboost::hash_combine(res, pathCodes[pathCodes.size() - i - 1]);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     for (unsigned int pathCode : pathCodes) {
    // RDKit✔️✔️:       gboost::hash_combine(res, pathCode);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AtomPairs::getTopologicalTorsionHash
    // As in the code-packing function, COSMolKit rejects the source's empty
    // path underflow instead of attempting to reproduce undefined behavior.
    if path_codes.is_empty() {
        return Err(FingerprintError::InvalidArguments {
            reason: "topological torsion hash path must not be empty",
        });
    }

    let reverse_it = topological_torsion_reverse(path_codes);
    let mut result = 0u32;
    if reverse_it {
        for path_code in path_codes.iter().rev() {
            hash_combine(&mut result, *path_code);
        }
    } else {
        for path_code in path_codes {
            hash_combine(&mut result, *path_code);
        }
    }
    Ok(result)
}

// RDKit source: FingerprintUtil.cpp lines 98-112 (smartsPatterns)
// RDKit✔️✔️: const char *smartsPatterns[6] = {
// RDKit✔️✔️:     "[$([N;!H0;v3,v4&+1]),\
// RDKit✔️✔️: $([O,S;H1;+0]),\
// RDKit✔️✔️: n&H1&+0]",                                                  // Donor
// RDKit✔️✔️:     "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),\
// RDKit✔️✔️: $([O,S;H0;v2]),\
// RDKit✔️✔️: $([O,S;-]),\
// RDKit✔️✔️: $([N;v3;!$(N-*=[O,N,P,S])]),\
// RDKit✔️✔️: n&H0&+0,\
// RDKit✔️✔️: $([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]",                    // Acceptor
// RDKit✔️✔️:     "[a]",                                                  // Aromatic
// RDKit✔️✔️:     "[F,Cl,Br,I]",                                          // Halogen
// RDKit✔️✔️:     "[#7;+,\
// RDKit✔️✔️: $([N;H2&+0][$([C,a]);!$([C,a](=O))]),\
// RDKit✔️✔️: $([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),\
// RDKit✔️✔️: $([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",  // Basic
// RDKit✔️✔️:     "[$([C,S](=[O,S,P])-[O;H1,-1])]"                        // Acidic
// RDKit✔️✔️: };
// RDKit source: FingerprintUtil.cpp lines 132-133
// RDKit✔️✔️: std::vector<std::string> defaultFeatureSmarts(smartsPatterns,
// RDKit✔️✔️:                                               smartsPatterns + 6);
const DEFAULT_FEATURE_SMARTS: [&str; 6] = [
    "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]",
    "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]",
    "[a]",
    "[F,Cl,Br,I]",
    "[#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",
    "[$([C,S](=[O,S,P])-[O;H1,-1])]",
];

// RDKit source: FingerprintUtil.cpp lines 115-129 (ss_matcher)
// RDKit✔️✔️: const RDKit::ROMol *ss_matcher::getMatcher() const { return m_matcher.get(); }
// RDKit✔️✔️: ss_matcher::ss_matcher() {};
// RDKit✔️✔️: ss_matcher::ss_matcher(const std::string &pattern) {
// RDKit✔️✔️:   RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
// RDKit✔️✔️:   TEST_ASSERT(p);
// RDKit✔️✔️:   m_matcher.reset(p);
// RDKit✔️✔️: };
//
// COSMolKit owns the equivalent query-bearing `Molecule` directly instead of
// caching an `ROMOL_SPTR`; matcher construction and access retain the same
// observable behavior.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SsMatcher {
    matcher: Molecule,
}

impl SsMatcher {
    #[must_use]
    pub fn try_new(pattern: &str) -> Result<Self, FingerprintError> {
        // RDKit✔️✔️:   RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
        // RDKit✔️✔️:   TEST_ASSERT(p);
        // RDKit✔️✔️:   m_matcher.reset(p);
        // Local complexity review: one canonical linear SMARTS compilation
        // replaces the transitional conversion and stores its query graph.
        let matcher = mol_from_smarts(pattern, &SmartsParseParams::default()).map_err(|error| {
            FingerprintError::InvalidSmartsPattern {
                pattern: pattern.to_string(),
                reason: error.to_string(),
            }
        })?;
        Ok(Self { matcher })
    }

    #[must_use]
    pub fn from_pattern(pattern: &str) -> Result<Self, FingerprintError> {
        Self::try_new(pattern)
    }

    #[must_use]
    #[allow(non_snake_case)]
    pub fn getMatcher(&self) -> &Molecule {
        &self.matcher
    }
}

#[must_use]
pub(crate) fn default_feature_smarts() -> &'static [&'static str; 6] {
    &DEFAULT_FEATURE_SMARTS
}

#[must_use]
pub(crate) fn default_feature_matchers() -> Result<Vec<SsMatcher>, FingerprintError> {
    DEFAULT_FEATURE_SMARTS
        .iter()
        .map(|pattern| SsMatcher::from_pattern(pattern))
        .collect()
}

fn cached_default_feature_matchers() -> Result<&'static [SsMatcher], FingerprintError> {
    static CACHE: OnceLock<Result<Vec<SsMatcher>, FingerprintError>> = OnceLock::new();
    match CACHE.get_or_init(default_feature_matchers) {
        Ok(matchers) => Ok(matchers.as_slice()),
        Err(error) => Err(error.clone()),
    }
}

fn cached_rdkit_maccs_pattern_matchers() -> Result<&'static [(usize, SsMatcher)], FingerprintError>
{
    static CACHE: OnceLock<Result<Vec<(usize, SsMatcher)>, FingerprintError>> = OnceLock::new();
    match CACHE.get_or_init(|| {
        RDKIT_MACCS_PATTERNS
            .iter()
            .map(|pattern| {
                SsMatcher::from_pattern(pattern.smarts).map(|matcher| (pattern.bit, matcher))
            })
            .collect()
    }) {
        Ok(matchers) => Ok(matchers.as_slice()),
        Err(error) => Err(error.clone()),
    }
}

/// RDKit-style scratch container for optional fingerprint provenance outputs.
///
/// RDKit source: `FingerprintGenerator.h` lines 27-75 (`AdditionalOutput`).
/// RDKit✔️✔️: atomToBits / bitInfoMap / bitPaths / atomCounts / atomsPerBit
/// RDKit✔️✔️: are optional allocations that start null and are populated on
/// RDKit✔️✔️: demand through explicit allocate* methods.
/// RDKit source: `FingerprintGenerator.cpp` lines 134-147 (`reinitAdditionalOutput`).
/// RDKit✔️✔️: allocated additional-output containers are reinitialized for the
/// RDKit✔️✔️: current atom count before fingerprint generation.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AdditionalOutput {
    pub atom_to_bits: Option<Vec<Vec<u64>>>,
    pub bit_info_map: Option<BTreeMap<u64, Vec<(u32, u32)>>>,
    pub bit_paths: Option<BTreeMap<u64, Vec<Vec<usize>>>>,
    pub atom_counts: Option<Vec<u32>>,
    pub atoms_per_bit: Option<BTreeMap<u64, Vec<Vec<usize>>>>,
}

impl AdditionalOutput {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// RDKit✔️✔️: void allocateAtomToBits() { atomToBitsHolder.reset(new ...); }
    pub fn allocate_atom_to_bits(&mut self) {
        self.atom_to_bits = Some(Vec::new());
    }

    /// RDKit✔️✔️: void allocateBitInfoMap() { bitInfoMapHolder.reset(new ...); }
    pub fn allocate_bit_info_map(&mut self) {
        self.bit_info_map = Some(BTreeMap::new());
    }

    /// RDKit✔️✔️: void allocateBitPaths() { bitPathsHolder.reset(new ...); }
    pub fn allocate_bit_paths(&mut self) {
        self.bit_paths = Some(BTreeMap::new());
    }

    /// RDKit✔️✔️: void allocateAtomCounts() { atomCountsHolder.reset(new ...); }
    pub fn allocate_atom_counts(&mut self) {
        self.atom_counts = Some(Vec::new());
    }

    /// RDKit✔️✔️: void allocateAtomsPerBit() { atomsPerBitHolder.reset(new ...); }
    pub fn allocate_atoms_per_bit(&mut self) {
        self.atoms_per_bit = Some(BTreeMap::new());
    }

    /// RDKit source: FingerprintGenerator.cpp lines 134-147
    /// RDKit✔️✔️: void reinitAdditionalOutput(AdditionalOutput &ao, size_t numAtoms) {
    /// RDKit✔️✔️:   if (ao.atomCounts) {
    /// RDKit✔️✔️:     ao.atomCounts->resize(numAtoms);
    /// RDKit✔️✔️:     std::fill(ao.atomCounts->begin(), ao.atomCounts->end(), 0);
    /// RDKit✔️✔️:   }
    /// RDKit✔️✔️:   if (ao.atomToBits) {
    /// RDKit✔️✔️:     ao.atomToBits->resize(numAtoms);
    /// RDKit✔️✔️:     std::fill(ao.atomToBits->begin(), ao.atomToBits->end(),
    /// RDKit✔️✔️:               std::vector<std::uint64_t>());
    /// RDKit✔️✔️:   }
    /// RDKit✔️✔️:   if (ao.bitInfoMap) ao.bitInfoMap->clear();
    /// RDKit✔️✔️:   if (ao.bitPaths) ao.bitPaths->clear();
    /// RDKit✔️✔️: }
    pub fn reset_for_atom_count(&mut self, num_atoms: usize) {
        if let Some(atom_counts) = self.atom_counts.as_mut() {
            atom_counts.resize(num_atoms, 0);
            atom_counts.fill(0);
        }
        if let Some(atom_to_bits) = self.atom_to_bits.as_mut() {
            atom_to_bits.resize_with(num_atoms, Vec::new);
            for entry in atom_to_bits.iter_mut() {
                entry.clear();
            }
        }
        if let Some(bit_info_map) = self.bit_info_map.as_mut() {
            bit_info_map.clear();
        }
        if let Some(bit_paths) = self.bit_paths.as_mut() {
            bit_paths.clear();
        }
    }
}

/// RDKit-style common fingerprint arguments shared by Morgan/MACCS wrappers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FingerprintArguments {
    pub df_count_simulation: bool,
    pub df_include_chirality: bool,
    pub d_count_bounds: Vec<u32>,
    pub d_fp_size: u32,
    pub d_num_bits_per_feature: u32,
}

impl Default for FingerprintArguments {
    fn default() -> Self {
        Self {
            df_count_simulation: false,
            df_include_chirality: false,
            d_count_bounds: Vec::new(),
            d_fp_size: 2048,
            d_num_bits_per_feature: 1,
        }
    }
}

impl FingerprintArguments {
    /// RDKit✔️✔️: FingerprintArguments(bool countSimulation, const std::vector<std::uint32_t> countBounds,
    /// RDKit✔️✔️:                       std::uint32_t fpSize, std::uint32_t numBitsPerFeature = 1,
    /// RDKit✔️✔️:                       bool includeChirality = false);
    #[must_use]
    pub fn new(
        count_simulation: bool,
        count_bounds: Vec<u32>,
        fp_size: u32,
        num_bits_per_feature: u32,
        include_chirality: bool,
    ) -> Result<Self, FingerprintError> {
        // RDKit source: FingerprintGenerator.cpp lines 36-48
        // RDKit✔️✔️: FingerprintArguments::FingerprintArguments(
        // RDKit✔️✔️:     const bool countSimulation, const std::vector<std::uint32_t> countBounds,
        // RDKit✔️✔️:     std::uint32_t fpSize, std::uint32_t numBitsPerFeature,
        // RDKit✔️✔️:     bool includeChirality)
        // RDKit✔️✔️:     : df_countSimulation(countSimulation),
        // RDKit✔️✔️:       df_includeChirality(includeChirality),
        // RDKit✔️✔️:       d_countBounds(countBounds),
        // RDKit✔️✔️:       d_fpSize(fpSize),
        // RDKit✔️✔️:       d_numBitsPerFeature(numBitsPerFeature) {
        // RDKit✔️✔️:   PRECONDITION(!countSimulation || !countBounds.empty(),
        // RDKit✔️✔️:                "bad count bounds provided");
        // RDKit✔️✔️:   PRECONDITION(d_numBitsPerFeature > 0, "numBitsPerFeature must be >0");
        // RDKit✔️✔️: }
        if count_simulation && count_bounds.is_empty() {
            return Err(FingerprintError::InvalidArguments {
                reason: "countSimulation requires non-empty countBounds",
            });
        }
        if num_bits_per_feature == 0 {
            return Err(FingerprintError::InvalidArguments {
                reason: "numBitsPerFeature must be > 0",
            });
        }
        Ok(Self {
            df_count_simulation: count_simulation,
            df_include_chirality: include_chirality,
            d_count_bounds: count_bounds,
            d_fp_size: fp_size,
            d_num_bits_per_feature: num_bits_per_feature,
        })
    }

    /// RDKit✔️✔️: return "Common arguments : countSimulation=... fpSize=... bitsPerFeature=... includeChirality=...";
    #[must_use]
    pub fn common_arguments_string(&self) -> String {
        // RDKit source: FingerprintGenerator.cpp lines 50-54
        // RDKit✔️✔️: std::string FingerprintArguments::commonArgumentsString() const {
        // RDKit✔️✔️:   return "Common arguments : countSimulation=" +
        // RDKit✔️✔️:          std::to_string(df_countSimulation) +
        // RDKit✔️✔️:          " fpSize=" + std::to_string(d_fpSize) +
        // RDKit✔️✔️:          " bitsPerFeature=" + std::to_string(d_numBitsPerFeature) +
        // RDKit✔️✔️:          " includeChirality=" + std::to_string(df_includeChirality);
        // RDKit✔️✔️: }
        format!(
            "Common arguments : countSimulation={} fpSize={} bitsPerFeature={} includeChirality={}",
            self.df_count_simulation as u8,
            self.d_fp_size,
            self.d_num_bits_per_feature,
            self.df_include_chirality as u8
        )
    }

    #[allow(non_snake_case)]
    #[must_use]
    pub fn commonArgumentsString(&self) -> String {
        self.common_arguments_string()
    }

    /// RDKit✔️✔️: pt.put("countSimulation", df_countSimulation); ... pt.add_child("countBounds", countBoundsNode);
    #[must_use]
    pub fn to_json(&self) -> String {
        // RDKit source: FingerprintGenerator.cpp lines 58-71
        // RDKit✔️✔️: void FingerprintArguments::toJSON(boost::property_tree::ptree &pt) const {
        // RDKit✔️✔️:   pt.put("countSimulation", df_countSimulation);
        // RDKit✔️✔️:   pt.put("fpSize", d_fpSize);
        // RDKit✔️✔️:   pt.put("numBitsPerFeature", d_numBitsPerFeature);
        // RDKit✔️✔️:   pt.put("includeChirality", df_includeChirality);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   boost::property_tree::ptree countBoundsNode;
        // RDKit✔️✔️:   for (const auto &bound : d_countBounds) {
        // RDKit✔️✔️:     boost::property_tree::ptree boundNode;
        // RDKit✔️✔️:     boundNode.put("", bound);
        // RDKit✔️✔️:     countBoundsNode.push_back(std::make_pair("", boundNode));
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   pt.add_child("countBounds", countBoundsNode);
        // RDKit✔️✔️: }
        let count_bounds = self
            .d_count_bounds
            .iter()
            .map(|bound| format!("\"{bound}\""))
            .collect::<Vec<_>>()
            .join(",");
        // Boost property_tree stores every leaf as text, so write_json quotes
        // numeric and boolean values even though fromJSON reads typed values.
        format!(
            "{{\"countSimulation\":\"{}\",\"fpSize\":\"{}\",\"numBitsPerFeature\":\"{}\",\"includeChirality\":\"{}\",\"countBounds\":[{}]}}",
            self.df_count_simulation,
            self.d_fp_size,
            self.d_num_bits_per_feature,
            self.df_include_chirality,
            count_bounds
        )
    }

    #[allow(non_snake_case)]
    #[must_use]
    pub fn toJSON(&self) -> String {
        self.to_json()
    }

    // RDKit✔️✔️: df_countSimulation = pt.get<bool>("countSimulation", df_countSimulation);
    /// RDKit✔️✔️: d_countBounds.clear();
    /// RDKit✔️✔️: if (countBoundsNode) { ... push_back(bound); }
    pub fn from_json(&mut self, json: &str) -> Result<(), FingerprintError> {
        // RDKit source: FingerprintGenerator.cpp lines 73-90
        // RDKit✔️✔️: void FingerprintArguments::fromJSON(const boost::property_tree::ptree &pt) {
        // RDKit✔️✔️:   df_countSimulation = pt.get<bool>("countSimulation", df_countSimulation);
        // RDKit✔️✔️:   d_fpSize = pt.get<std::uint32_t>("fpSize", d_fpSize);
        // RDKit✔️✔️:   d_numBitsPerFeature =
        // RDKit✔️✔️:       pt.get<std::uint32_t>("numBitsPerFeature", d_numBitsPerFeature);
        // RDKit✔️✔️:   df_includeChirality = pt.get<bool>("includeChirality", df_includeChirality);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   d_countBounds.clear();
        // RDKit✔️✔️:   auto countBoundsNode = pt.get_child_optional("countBounds");
        // RDKit✔️✔️:   if (countBoundsNode) {
        // RDKit✔️✔️:     for (const auto &boundNode : *countBoundsNode) {
        // RDKit✔️✔️:       d_countBounds.push_back(boundNode.second.get_value<std::uint32_t>());
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        if json.trim().is_empty() {
            return Ok(());
        }
        let value: Value = serde_json::from_str(json)
            .map_err(|error| FingerprintError::InvalidArgumentsJson(error.to_string()))?;
        self.from_json_value(&value)
    }

    #[allow(non_snake_case)]
    pub fn fromJSON(&mut self, json: &str) -> Result<(), FingerprintError> {
        self.from_json(json)
    }

    fn from_json_value(&mut self, value: &Value) -> Result<(), FingerprintError> {
        let object = value.as_object().ok_or_else(|| {
            FingerprintError::InvalidArgumentsJson("expected JSON object".to_string())
        })?;

        if let Some(field) = object.get("countSimulation") {
            self.df_count_simulation = json_value_as_bool("countSimulation", field)?;
        }
        if let Some(field) = object.get("fpSize") {
            self.d_fp_size = json_value_as_u32("fpSize", field)?;
        }
        if let Some(field) = object.get("numBitsPerFeature") {
            self.d_num_bits_per_feature = json_value_as_u32("numBitsPerFeature", field)?;
        }
        if let Some(field) = object.get("includeChirality") {
            self.df_include_chirality = json_value_as_bool("includeChirality", field)?;
        }

        self.d_count_bounds.clear();
        if let Some(field) = object.get("countBounds") {
            let bounds = field.as_array().ok_or_else(|| {
                FingerprintError::InvalidArgumentsJson("countBounds must be an array".to_string())
            })?;
            for bound in bounds {
                self.d_count_bounds
                    .push(json_value_as_u32("countBounds entry", bound)?);
            }
        }
        Ok(())
    }
}

/// RDKit Topological Torsion fingerprint arguments layered over the shared
/// fingerprint-generator arguments.
// RDKit source: TopologicalTorsionGenerator.h lines 21-50
// RDKit✔️✔️: class RDKIT_FINGERPRINTS_EXPORT TopologicalTorsionArguments
// RDKit✔️✔️:     : public FingerprintArguments {
// RDKit✔️✔️:  public:
// RDKit✔️✔️:   uint32_t d_torsionAtomCount = 4;
// RDKit✔️✔️:   bool df_onlyShortestPaths = false;
// RDKit✔️✔️:
// RDKit✔️✔️:   std::string infoString() const override;
// RDKit✔️✔️:   void toJSON(boost::property_tree::ptree &pt) const override;
// RDKit✔️✔️:   void fromJSON(const boost::property_tree::ptree &pt) override;
// RDKit✔️✔️:
// RDKit✔️✔️:   TopologicalTorsionArguments(
// RDKit✔️✔️:       const bool includeChirality = false, const uint32_t torsionAtomCount = 4,
// RDKit✔️✔️:       const bool countSimulation = true,
// RDKit✔️✔️:       const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
// RDKit✔️✔️:       const std::uint32_t fpSize = 2048);
// RDKit✔️✔️: };
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologicalTorsionArguments {
    pub fingerprint_arguments: FingerprintArguments,
    pub d_torsion_atom_count: u32,
    pub df_only_shortest_paths: bool,
}

impl Default for TopologicalTorsionArguments {
    fn default() -> Self {
        Self {
            fingerprint_arguments: FingerprintArguments {
                df_count_simulation: true,
                df_include_chirality: false,
                d_count_bounds: vec![1, 2, 4, 8],
                d_fp_size: 2048,
                d_num_bits_per_feature: 1,
            },
            d_torsion_atom_count: 4,
            df_only_shortest_paths: false,
        }
    }
}

impl TopologicalTorsionArguments {
    #[must_use]
    pub fn new(
        include_chirality: bool,
        torsion_atom_count: u32,
        count_simulation: bool,
        count_bounds: Vec<u32>,
        fp_size: u32,
    ) -> Result<Self, FingerprintError> {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 25-32
        // RDKit✔️✔️: TopologicalTorsionArguments::TopologicalTorsionArguments(
        // RDKit✔️✔️:     const bool includeChirality, const uint32_t torsionAtomCount,
        // RDKit✔️✔️:     const bool countSimulation, const std::vector<std::uint32_t> countBounds,
        // RDKit✔️✔️:     const std::uint32_t fpSize)
        // RDKit✔️✔️:     : FingerprintArguments(countSimulation, countBounds, fpSize, 1,
        // RDKit✔️✔️:                            includeChirality),
        // RDKit✔️✔️:       d_torsionAtomCount(torsionAtomCount) {};
        let fingerprint_arguments = FingerprintArguments::new(
            count_simulation,
            count_bounds,
            fp_size,
            1,
            include_chirality,
        )?;
        Ok(Self {
            fingerprint_arguments,
            d_torsion_atom_count: torsion_atom_count,
            df_only_shortest_paths: false,
        })
    }

    #[must_use]
    pub fn info_string(&self) -> String {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 47-51
        // RDKit✔️✔️: std::string TopologicalTorsionArguments::infoString() const {
        // RDKit✔️✔️:   return "TopologicalTorsionArguments torsionAtomCount=" +
        // RDKit✔️✔️:          std::to_string(d_torsionAtomCount) +
        // RDKit✔️✔️:          " onlyShortestPaths=" + std::to_string(df_onlyShortestPaths);
        // RDKit✔️✔️: };
        format!(
            "TopologicalTorsionArguments torsionAtomCount={} onlyShortestPaths={}",
            self.d_torsion_atom_count, self.df_only_shortest_paths as u8
        )
    }

    #[allow(non_snake_case)]
    #[must_use]
    pub fn infoString(&self) -> String {
        self.info_string()
    }

    #[must_use]
    pub fn to_json(&self) -> String {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 52-58
        // RDKit✔️✔️: void TopologicalTorsionArguments::toJSON(
        // RDKit✔️✔️:     boost::property_tree::ptree &pt) const {
        // RDKit✔️✔️:   pt.put("type", "TopologicalTorsionArguments");
        // RDKit✔️✔️:   pt.put("torsionAtomCount", d_torsionAtomCount);
        // RDKit✔️✔️:   pt.put("onlyShortestPaths", df_onlyShortestPaths);
        // RDKit✔️✔️:   FingerprintArguments::toJSON(pt);
        // RDKit✔️✔️: }
        let common = self.fingerprint_arguments.to_json();
        let common_body = &common[1..common.len().saturating_sub(1)];
        let mut json = format!(
            "{{\"type\":\"TopologicalTorsionArguments\",\"torsionAtomCount\":\"{}\",\"onlyShortestPaths\":\"{}\"",
            self.d_torsion_atom_count, self.df_only_shortest_paths
        );
        if !common_body.is_empty() {
            json.push(',');
            json.push_str(common_body);
        }
        json.push('}');
        json
    }

    #[allow(non_snake_case)]
    #[must_use]
    pub fn toJSON(&self) -> String {
        self.to_json()
    }

    pub fn from_json(&mut self, json: &str) -> Result<(), FingerprintError> {
        if json.trim().is_empty() {
            return Ok(());
        }
        let value: Value = serde_json::from_str(json)
            .map_err(|error| FingerprintError::InvalidArgumentsJson(error.to_string()))?;
        self.from_json_value(&value)
    }

    #[allow(non_snake_case)]
    pub fn fromJSON(&mut self, json: &str) -> Result<(), FingerprintError> {
        self.from_json(json)
    }

    fn from_json_value(&mut self, value: &Value) -> Result<(), FingerprintError> {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 59-65
        // RDKit✔️✔️: void TopologicalTorsionArguments::fromJSON(
        // RDKit✔️✔️:     const boost::property_tree::ptree &pt) {
        // RDKit✔️✔️:   d_torsionAtomCount = pt.get<uint32_t>("torsionAtomCount", d_torsionAtomCount);
        // RDKit✔️✔️:   df_onlyShortestPaths =
        // RDKit✔️✔️:       pt.get<bool>("onlyShortestPaths", df_onlyShortestPaths);
        // RDKit✔️✔️:   FingerprintArguments::fromJSON(pt);
        // RDKit✔️✔️: }
        let object = value.as_object().ok_or_else(|| {
            FingerprintError::InvalidArgumentsJson("expected JSON object".to_string())
        })?;
        if let Some(field) = object.get("torsionAtomCount") {
            self.d_torsion_atom_count = json_value_as_u32("torsionAtomCount", field)?;
        }
        if let Some(field) = object.get("onlyShortestPaths") {
            self.df_only_shortest_paths = json_value_as_bool("onlyShortestPaths", field)?;
        }
        self.fingerprint_arguments.from_json_value(value)
    }
}

/// One source `TopologicalTorsionAtomEnv` produced by the Topological Torsion
/// environment generator.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologicalTorsionAtomEnv {
    bit_id: u64,
    atom_path: Vec<usize>,
}

impl TopologicalTorsionAtomEnv {
    #[must_use]
    pub fn new(bit_id: u64, atom_path: Vec<usize>) -> Self {
        // RDKit source: TopologicalTorsionGenerator.h lines 69-70
        // RDKit✔️✔️: TopologicalTorsionAtomEnv(OutputType bitId, INT_VECT atomPath)
        // RDKit✔️✔️:     : d_bitId(bitId), d_atomPath(std::move(atomPath)) {}
        Self { bit_id, atom_path }
    }

    #[must_use]
    #[allow(non_snake_case)]
    pub fn getBitId(&self) -> u64 {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 89-101
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: OutputType TopologicalTorsionAtomEnv<OutputType>::getBitId(
        // RDKit✔️✔️:     FingerprintArguments *,              // arguments
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *,  // atomInvariants
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *,  // bondInvariants
        // RDKit✔️✔️:     AdditionalOutput *,                  // additionalOutput,
        // RDKit✔️✔️:     const bool,                          // hashResults
        // RDKit✔️✔️:     const std::uint64_t                  // fpSize
        // RDKit✔️✔️: ) const {
        // RDKit✔️✔️:   return d_bitId;
        // RDKit✔️✔️: };
        self.bit_id
    }

    #[allow(non_snake_case)]
    pub fn updateAdditionalOutput(&self, additional_output: &mut AdditionalOutput, bit_id: u64) {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 68-87
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: void TopologicalTorsionAtomEnv<OutputType>::updateAdditionalOutput(
        // RDKit✔️✔️:     AdditionalOutput *additionalOutput, size_t bitId) const {
        // RDKit✔️✔️:   PRECONDITION(additionalOutput, "bad output pointer");
        // RDKit✔️✔️:   if (additionalOutput->atomToBits || additionalOutput->atomCounts) {
        // RDKit✔️✔️:     for (auto aid : d_atomPath) {
        // RDKit✔️✔️:       if (additionalOutput->atomToBits) {
        // RDKit✔️✔️:         additionalOutput->atomToBits->at(aid).push_back(bitId);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (additionalOutput->atomCounts) {
        // RDKit✔️✔️:         additionalOutput->atomCounts->at(aid)++;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (additionalOutput->bitPaths) {
        // RDKit✔️✔️:     (*additionalOutput->bitPaths)[bitId].push_back(d_atomPath);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (additionalOutput->atomsPerBit) {
        // RDKit✔️✔️:     (*additionalOutput->atomsPerBit)[bitId].push_back(d_atomPath);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        if additional_output.atom_to_bits.is_some() || additional_output.atom_counts.is_some() {
            for &atom_id in &self.atom_path {
                if let Some(atom_to_bits) = additional_output.atom_to_bits.as_mut() {
                    atom_to_bits[atom_id].push(bit_id);
                }
                if let Some(atom_counts) = additional_output.atom_counts.as_mut() {
                    atom_counts[atom_id] += 1;
                }
            }
        }
        if let Some(bit_paths) = additional_output.bit_paths.as_mut() {
            bit_paths
                .entry(bit_id)
                .or_default()
                .push(self.atom_path.clone());
        }
        if let Some(atoms_per_bit) = additional_output.atoms_per_bit.as_mut() {
            atoms_per_bit
                .entry(bit_id)
                .or_default()
                .push(self.atom_path.clone());
        }
    }
}

/// RDKit Topological Torsion atom-environment generator. The source stores a
/// borrowed arguments pointer in the generator base; Rust passes that borrow
/// explicitly so there is only one owned arguments value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
enum TopologicalTorsionAtomCodeMode {
    #[default]
    Modern,
    LegacyUnfolded,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct TopologicalTorsionEnvGenerator {
    atom_code_mode: TopologicalTorsionAtomCodeMode,
}

impl TopologicalTorsionEnvGenerator {
    #[must_use]
    pub const fn new() -> Self {
        Self {
            atom_code_mode: TopologicalTorsionAtomCodeMode::Modern,
        }
    }

    fn use_legacy_unfolded_atom_codes(&mut self) {
        self.atom_code_mode = TopologicalTorsionAtomCodeMode::LegacyUnfolded;
    }

    #[must_use]
    #[allow(non_snake_case)]
    pub fn infoString(&self) -> String {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 196-200
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: std::string TopologicalTorsionEnvGenerator<OutputType>::infoString() const {
        // RDKit✔️✔️:   return "TopologicalTorsionEnvGenerator";
        // RDKit✔️✔️: };
        "TopologicalTorsionEnvGenerator".to_string()
    }

    #[must_use]
    #[allow(non_snake_case)]
    pub fn toJSON(&self) -> String {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 201-206
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: void TopologicalTorsionEnvGenerator<OutputType>::toJSON(
        // RDKit✔️✔️:     boost::property_tree::ptree &pt) const {
        // RDKit✔️✔️:   pt.put("type", "TopologicalTorsionEnvGenerator");
        // RDKit✔️✔️:   AtomEnvironmentGenerator<OutputType>::toJSON(pt);
        // RDKit✔️✔️: };
        r#"{"type":"TopologicalTorsionEnvGenerator"}"#.to_string()
    }

    #[allow(non_snake_case)]
    pub fn fromJSON(&mut self, json: &str) -> Result<(), FingerprintError> {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 207-211
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: void TopologicalTorsionEnvGenerator<OutputType>::fromJSON(
        // RDKit✔️✔️:     const boost::property_tree::ptree &pt) {
        // RDKit✔️✔️:   AtomEnvironmentGenerator<OutputType>::fromJSON(pt);
        // RDKit✔️✔️: };
        validate_stateless_environment_generator_json(json)
    }

    #[allow(non_snake_case)]
    pub fn getEnvironments(
        &self,
        molecule: &Molecule,
        arguments: &TopologicalTorsionArguments,
        from_atoms: Option<&[usize]>,
        ignore_atoms: Option<&[usize]>,
        atom_invariants: &[u32],
        hash_results: bool,
    ) -> Result<Vec<TopologicalTorsionAtomEnv>, FingerprintError> {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 101-194
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: std::vector<AtomEnvironment<OutputType> *>
        // RDKit✔️✔️: TopologicalTorsionEnvGenerator<OutputType>::getEnvironments(
        // RDKit✔️✔️:     const ROMol &mol, FingerprintArguments *arguments,
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *fromAtoms,
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *ignoreAtoms,
        // RDKit✔️✔️:     const int,                 // confId
        // RDKit✔️✔️:     const AdditionalOutput *,  // additionalOutput
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *atomInvariants,
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *,  // bondInvariants
        // RDKit✔️✔️:     const bool hashResults) const {
        // RDKit✔️✔️:   auto *topologicalTorsionArguments =
        // RDKit✔️✔️:       dynamic_cast<TopologicalTorsionArguments *>(arguments);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   std::vector<AtomEnvironment<OutputType> *> result;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   boost::dynamic_bitset<> *fromAtomsBV = nullptr;
        // RDKit✔️✔️:   if (fromAtoms) {
        // RDKit✔️✔️:     fromAtomsBV = new boost::dynamic_bitset<>(mol.getNumAtoms());
        // RDKit✔️✔️:     for (auto fAt : *fromAtoms) {
        // RDKit✔️✔️:       fromAtomsBV->set(fAt);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   boost::dynamic_bitset<> *ignoreAtomsBV = nullptr;
        // RDKit✔️✔️:   if (ignoreAtoms) {
        // RDKit✔️✔️:     ignoreAtomsBV = new boost::dynamic_bitset<>(mol.getNumAtoms());
        // RDKit✔️✔️:     for (auto fAt : *ignoreAtoms) {
        // RDKit✔️✔️:       ignoreAtomsBV->set(fAt);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   boost::dynamic_bitset<> pAtoms(mol.getNumAtoms());
        // RDKit✔️✔️:   bool useBonds = false;
        // RDKit✔️✔️:   bool useHs = false;
        // RDKit✔️✔️:   int rootedAtAtom = -1;
        // RDKit✔️✔️:   PATH_LIST paths = findAllPathsOfLengthN(
        // RDKit✔️✔️:       mol, topologicalTorsionArguments->d_torsionAtomCount, useBonds, useHs,
        // RDKit✔️✔️:       rootedAtAtom, topologicalTorsionArguments->df_onlyShortestPaths);
        // RDKit✔️✔️:   for (PATH_LIST::const_iterator pathIt = paths.begin(); pathIt != paths.end();
        // RDKit✔️✔️:        ++pathIt) {
        // RDKit✔️✔️:     bool keepIt = true;
        // RDKit✔️✔️:     if (fromAtomsBV) {
        // RDKit✔️✔️:       keepIt = false;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     std::vector<std::uint32_t> pathCodes;
        // RDKit✔️✔️:     const PATH_TYPE &path = *pathIt;
        // RDKit✔️✔️:     if (fromAtomsBV) {
        // RDKit✔️✔️:       if (fromAtomsBV->test(static_cast<std::uint32_t>(path.front())) ||
        // RDKit✔️✔️:           fromAtomsBV->test(static_cast<std::uint32_t>(path.back()))) {
        // RDKit✔️✔️:         keepIt = true;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (keepIt && ignoreAtomsBV) {
        // RDKit✔️✔️:       for (int pElem : path) {
        // RDKit✔️✔️:         if (ignoreAtomsBV->test(pElem)) {
        // RDKit✔️✔️:           keepIt = false;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (keepIt) {
        // RDKit✔️✔️:       pAtoms.reset();
        // RDKit✔️✔️:       for (auto pIt = path.begin(); pIt < path.end(); ++pIt) {
        // RDKit✔️✔️:         // look for a cycle that doesn't start at the first atom
        // RDKit✔️✔️:         // we can't effectively canonicalize these at the moment
        // RDKit✔️✔️:         // (was github #811)
        // RDKit✔️✔️:         if (pIt != path.begin() && *pIt != *(path.begin()) && pAtoms[*pIt]) {
        // RDKit✔️✔️:           pathCodes.clear();
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:         pAtoms.set(*pIt);
        // RDKit✔️✔️:         unsigned int code = (*atomInvariants)[*pIt] % ((1 << codeSize) - 1) + 1;
        // RDKit✔️✔️:         // subtract off the branching number:
        // RDKit✔️✔️:         if (pIt != path.begin() && pIt + 1 != path.end()) {
        // RDKit✔️✔️:           --code;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:         pathCodes.push_back(code);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (pathCodes.size()) {
        // RDKit✔️✔️:         OutputType code;
        // RDKit✔️✔️:         if (hashResults) {
        // RDKit✔️✔️:           code = getTopologicalTorsionHash(pathCodes);
        // RDKit✔️✔️:         } else {
        // RDKit✔️✔️:           code = getTopologicalTorsionCode(
        // RDKit✔️✔️:               pathCodes, topologicalTorsionArguments->df_includeChirality);
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:         result.push_back(new TopologicalTorsionAtomEnv<OutputType>(code, path));
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   delete fromAtomsBV;
        // RDKit✔️✔️:   delete ignoreAtomsBV;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return result;
        // RDKit✔️✔️: };
        let atom_count = molecule.num_atoms();
        if atom_invariants.len() < atom_count {
            return Err(FingerprintError::InvalidArguments {
                reason: "bad atom invariants size",
            });
        }
        if from_atoms
            .into_iter()
            .chain(ignore_atoms)
            .flatten()
            .any(|&atom_id| atom_id >= atom_count)
        {
            return Err(FingerprintError::InvalidArguments {
                reason: "atom selection contains an atom index outside the molecule",
            });
        }

        let mut from_atoms_bits = from_atoms.map(|_| vec![false; atom_count]);
        if let (Some(bits), Some(atom_ids)) = (from_atoms_bits.as_mut(), from_atoms) {
            for &atom_id in atom_ids {
                bits[atom_id] = true;
            }
        }
        let mut ignore_atoms_bits = ignore_atoms.map(|_| vec![false; atom_count]);
        if let (Some(bits), Some(atom_ids)) = (ignore_atoms_bits.as_mut(), ignore_atoms) {
            for &atom_id in atom_ids {
                bits[atom_id] = true;
            }
        }

        let paths = find_all_paths_of_length_n(
            molecule,
            arguments.d_torsion_atom_count as usize,
            false,
            false,
            -1,
            arguments.df_only_shortest_paths,
        )?;
        let mut result = Vec::new();
        let mut path_atoms = vec![false; atom_count];
        let code_modulus = (1_u32 << ATOM_PAIR_CODE_SIZE) - 1;
        for path in paths {
            let mut keep = from_atoms_bits.is_none();
            if let Some(bits) = from_atoms_bits.as_ref() {
                let (Some(&first), Some(&last)) = (path.first(), path.last()) else {
                    return Err(FingerprintError::InvalidArguments {
                        reason: "enumerated topological torsion path is empty",
                    });
                };
                keep = bits[first] || bits[last];
            }
            if keep
                && ignore_atoms_bits
                    .as_ref()
                    .is_some_and(|bits| path.iter().any(|&atom_id| bits[atom_id]))
            {
                keep = false;
            }
            if !keep {
                continue;
            }

            path_atoms.fill(false);
            let mut path_codes = Vec::new();
            for (position, &atom_id) in path.iter().enumerate() {
                if position != 0 && atom_id != path[0] && path_atoms[atom_id] {
                    path_codes.clear();
                    break;
                }
                path_atoms[atom_id] = true;
                // The modern generator and deprecated unfolded API share this
                // traversal but intentionally use different source formulas.
                // RDKit✔️✔️: unsigned int code = (*atomInvariants)[*pIt] % ((1 << codeSize) - 1) + 1;
                // RDKit✔️✔️: unsigned int code = atomCodes[*pIt] - 1;
                let mut code = match self.atom_code_mode {
                    TopologicalTorsionAtomCodeMode::Modern => {
                        atom_invariants[atom_id] % code_modulus + 1
                    }
                    TopologicalTorsionAtomCodeMode::LegacyUnfolded => {
                        atom_invariants[atom_id].wrapping_add(1)
                    }
                };
                if position != 0 && position + 1 != path.len() {
                    code -= 1;
                }
                path_codes.push(code);
            }
            if path_codes.is_empty() {
                continue;
            }
            let bit_id = if hash_results {
                u64::from(get_topological_torsion_hash(&path_codes)?)
            } else {
                get_topological_torsion_code(
                    &path_codes,
                    arguments.fingerprint_arguments.df_include_chirality,
                )?
            };
            result.push(TopologicalTorsionAtomEnv::new(bit_id, path));
        }
        Ok(result)
    }

    pub fn get_result_size(
        &self,
        arguments: &TopologicalTorsionArguments,
    ) -> Result<u64, FingerprintError> {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 33-45
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: OutputType TopologicalTorsionEnvGenerator<OutputType>::getResultSize() const {
        // RDKit✔️✔️:   OutputType result = 1;
        // RDKit✔️✔️:   return (result << ((
        // RDKit✔️✔️:               dynamic_cast<const TopologicalTorsionArguments *>(
        // RDKit✔️✔️:                   this->dp_fingerprintArguments)
        // RDKit✔️✔️:                   ->d_torsionAtomCount *
        // RDKit✔️✔️:               (codeSize + (dynamic_cast<const TopologicalTorsionArguments *>(
        // RDKit✔️✔️:                                this->dp_fingerprintArguments)
        // RDKit✔️✔️:                                    ->df_includeChirality
        // RDKit✔️✔️:                                ? numChiralBits
        // RDKit✔️✔️:                                : 0)))));
        // RDKit✔️✔️: };
        let bits_per_atom = ATOM_PAIR_CODE_SIZE
            + if arguments.fingerprint_arguments.df_include_chirality {
                ATOM_PAIR_NUM_CHIRAL_BITS
            } else {
                0
            };
        let shift = arguments
            .d_torsion_atom_count
            .checked_mul(bits_per_atom)
            .ok_or(FingerprintError::InvalidArguments {
                reason: "topological torsion result-size width overflow",
            })?;
        if shift >= u64::BITS {
            return Err(FingerprintError::InvalidArguments {
                reason: "topological torsion result-size shift must be less than 64 bits",
            });
        }
        Ok(1_u64 << shift)
    }

    #[allow(non_snake_case)]
    pub fn getResultSize(
        &self,
        arguments: &TopologicalTorsionArguments,
    ) -> Result<u64, FingerprintError> {
        self.get_result_size(arguments)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganFingerprintParams {
    pub radius: u32,
    pub n_bits: usize,
    pub use_chirality: bool,
    pub use_bond_types: bool,
    pub count_simulation: bool,
    pub count_bounds: Vec<u32>,
    pub only_nonzero_invariants: bool,
    pub include_ring_membership: bool,
    pub include_redundant_environments: bool,
    pub from_atoms: Option<Vec<usize>>,
    pub ignore_atoms: Option<Vec<usize>>,
    pub custom_atom_invariants: Option<Vec<u32>>,
    pub custom_bond_invariants: Option<Vec<u32>>,
    pub atom_invariants_generator: MorganAtomInvariantsGenerator,
    pub bond_invariants_generator: Option<MorganBondInvariantsGenerator>,
    pub num_bits_per_feature: u32,
    pub collect_additional_output: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct FingerprintFuncArguments {
    pub from_atoms: Option<Vec<usize>>,
    pub ignore_atoms: Option<Vec<usize>>,
    pub conf_id: i32,
    pub custom_atom_invariants: Option<Vec<u32>>,
    pub custom_bond_invariants: Option<Vec<u32>>,
    pub additional_output: Option<AdditionalOutput>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseCountFingerprint {
    size: u64,
    nonzero_elements: BTreeMap<u64, i32>,
}

impl SparseCountFingerprint {
    #[must_use]
    pub fn new(size: u64) -> Self {
        Self {
            size,
            nonzero_elements: BTreeMap::new(),
        }
    }

    #[must_use]
    pub fn size(&self) -> u64 {
        self.size
    }

    #[must_use]
    pub fn get_val(&self, bit_id: u64) -> i32 {
        self.nonzero_elements.get(&bit_id).copied().unwrap_or(0)
    }

    pub fn set_val(&mut self, bit_id: u64, value: i32) {
        if value == 0 {
            self.nonzero_elements.remove(&bit_id);
        } else {
            self.nonzero_elements.insert(bit_id, value);
        }
    }

    #[must_use]
    pub fn nonzero_elements(&self) -> &BTreeMap<u64, i32> {
        &self.nonzero_elements
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseBitFingerprint {
    size: u64,
    on_bits: BTreeSet<u64>,
}

impl SparseBitFingerprint {
    #[must_use]
    pub fn new(size: u64) -> Self {
        Self {
            size,
            on_bits: BTreeSet::new(),
        }
    }

    #[must_use]
    pub fn size(&self) -> u64 {
        self.size
    }

    pub fn set_bit(&mut self, bit_id: u64) {
        self.on_bits.insert(bit_id);
    }

    #[must_use]
    pub fn on_bits(&self) -> &BTreeSet<u64> {
        &self.on_bits
    }
}

struct RdkitFingerprintMtRng {
    state: [u32; 4],
    index: usize,
}

impl RdkitFingerprintMtRng {
    const W: usize = 32;
    const N: usize = 4;
    const M: usize = 2;
    const R: usize = 31;
    const A: u32 = 0x9908_b0df;
    const U: usize = 11;
    const D: u32 = u32::MAX;
    const S: usize = 7;
    const B: u32 = 0x9d2c_5680;
    const T: usize = 15;
    const C: u32 = 0xefc6_0000;
    const L: usize = 18;
    const F: u32 = 1_812_433_253;

    fn new(seed: u32) -> Self {
        let mut rng = Self {
            state: [0; Self::N],
            index: Self::N,
        };
        rng.seed(seed);
        rng
    }

    fn seed(&mut self, seed: u32) {
        // Boost✔️✔️: const UIntType mask = (max)();
        // Boost adapter note: RDKit's deprecated boost::random::mersenne_twister
        // template forwards to mersenne_twister_engine with f=1812433253.
        // Boost✔️✔️: x[0] = value & mask;
        // Boost✔️✔️: for (i = 1; i < n; i++) {
        // Boost✔️✔️:   x[i] = (f * (x[i-1] ^ (x[i-1] >> (w-2))) + i) & mask;
        // Boost✔️✔️: }
        // Boost✔️✔️: normalize_state();
        self.state[0] = seed;
        for idx in 1..Self::N {
            self.state[idx] = Self::F
                .wrapping_mul(self.state[idx - 1] ^ (self.state[idx - 1] >> (Self::W - 2)))
                .wrapping_add(idx as u32);
        }
        self.index = Self::N;
        self.normalize_state();
    }

    fn normalize_state(&mut self) {
        // Boost✔️✔️: const UIntType upper_mask = (~static_cast<UIntType>(0)) << r;
        // Boost✔️✔️: const UIntType lower_mask = ~upper_mask;
        // Boost✔️✔️: UIntType y0 = x[m-1] ^ x[n-1];
        // Boost✔️✔️: if(y0 & (static_cast<UIntType>(1) << (w-1))) {
        // Boost✔️✔️:   y0 = ((y0 ^ a) << 1) | 1;
        // Boost✔️✔️: } else {
        // Boost✔️✔️:   y0 = y0 << 1;
        // Boost✔️✔️: }
        // Boost✔️✔️: x[0] = (x[0] & upper_mask) | (y0 & lower_mask);
        // Boost✔️✔️: for(std::size_t j = 0; j < n; ++j) {
        // Boost✔️✔️:   if(x[j] != 0) return;
        // Boost✔️✔️: }
        // Boost✔️✔️: x[0] = static_cast<UIntType>(1) << (w-1);
        let upper_mask = u32::MAX << Self::R;
        let lower_mask = !upper_mask;
        let mut y0 = self.state[Self::M - 1] ^ self.state[Self::N - 1];
        if y0 & (1_u32 << (Self::W - 1)) != 0 {
            y0 = ((y0 ^ Self::A) << 1) | 1;
        } else {
            y0 <<= 1;
        }
        self.state[0] = (self.state[0] & upper_mask) | (y0 & lower_mask);
        if self.state.iter().all(|&value| value == 0) {
            self.state[0] = 1_u32 << (Self::W - 1);
        }
    }

    fn twist(&mut self) {
        // Boost✔️✔️: const UIntType upper_mask = (~static_cast<UIntType>(0)) << r;
        // Boost✔️✔️: const UIntType lower_mask = ~upper_mask;
        // Boost✔️✔️: UIntType y = (x[j] & upper_mask) | (x[j+1] & lower_mask);
        // Boost✔️✔️: x[j] = x[...] ^ (y >> 1) ^ ((x[...]&1) * a);
        let upper_mask = u32::MAX << Self::R;
        let lower_mask = !upper_mask;
        for idx in 0..(Self::N - Self::M) {
            let y = (self.state[idx] & upper_mask) | (self.state[idx + 1] & lower_mask);
            self.state[idx] = self.state[idx + Self::M]
                ^ (y >> 1)
                ^ ((self.state[idx + 1] & 1).wrapping_mul(Self::A));
        }
        for idx in (Self::N - Self::M)..(Self::N - 1) {
            let y = (self.state[idx] & upper_mask) | (self.state[idx + 1] & lower_mask);
            self.state[idx] = self.state[idx - (Self::N - Self::M)]
                ^ (y >> 1)
                ^ ((self.state[idx + 1] & 1).wrapping_mul(Self::A));
        }
        let y = (self.state[Self::N - 1] & upper_mask) | (self.state[0] & lower_mask);
        self.state[Self::N - 1] =
            self.state[Self::M - 1] ^ (y >> 1) ^ ((self.state[0] & 1).wrapping_mul(Self::A));
        self.index = 0;
    }

    fn next_u32(&mut self) -> u32 {
        // Boost✔️✔️: if(i == n) twist();
        // Boost✔️✔️: UIntType z = x[i]; ++i;
        // Boost✔️✔️: z ^= ((z >> u) & d);
        // Boost✔️✔️: z ^= ((z << s) & b);
        // Boost✔️✔️: z ^= ((z << t) & c);
        // Boost✔️✔️: z ^= (z >> l);
        // Boost✔️✔️: return z;
        if self.index == Self::N {
            self.twist();
        }
        let mut z = self.state[self.index];
        self.index += 1;
        z ^= (z >> Self::U) & Self::D;
        z ^= (z << Self::S) & Self::B;
        z ^= (z << Self::T) & Self::C;
        z ^= z >> Self::L;
        z
    }

    fn uniform_int_0_to_i32_max(&mut self) -> u32 {
        // Boost✔️✔️: dist.reset(new distrib_type(0, INT_MAX));
        // Boost✔️✔️: mixed_range_type bucket_size = brange / (range + 1);
        // Boost✔️✔️: result = eng(); result /= bucket_size;
        // Boost✔️✔️: if(result <= range) return result;
        loop {
            let result = self.next_u32() / 2;
            if result <= i32::MAX as u32 {
                return result;
            }
        }
    }
}

impl Default for MorganFingerprintParams {
    fn default() -> Self {
        Self {
            radius: 2,
            n_bits: 2048,
            use_chirality: false,
            use_bond_types: true,
            count_simulation: false,
            count_bounds: vec![1, 2, 4, 8],
            only_nonzero_invariants: false,
            include_ring_membership: true,
            include_redundant_environments: false,
            from_atoms: None,
            ignore_atoms: None,
            custom_atom_invariants: None,
            custom_bond_invariants: None,
            atom_invariants_generator: MorganAtomInvariantsGenerator::Connectivity {
                include_ring_membership: true,
            },
            bond_invariants_generator: None,
            num_bits_per_feature: 1,
            collect_additional_output: false,
        }
    }
}

/// RDKit-style Morgan fingerprint argument bundle.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganArguments {
    pub fingerprint_arguments: FingerprintArguments,
    pub df_only_nonzero_invariants: bool,
    pub d_radius: u32,
    pub df_include_redundant_environments: bool,
    pub df_use_bond_types: bool,
}

impl Default for MorganArguments {
    fn default() -> Self {
        Self {
            fingerprint_arguments: FingerprintArguments {
                df_count_simulation: false,
                df_include_chirality: false,
                d_count_bounds: vec![1, 2, 4, 8],
                d_fp_size: 2048,
                d_num_bits_per_feature: 1,
            },
            df_only_nonzero_invariants: false,
            d_radius: 3,
            df_include_redundant_environments: false,
            df_use_bond_types: true,
        }
    }
}

impl MorganArguments {
    /// RDKit✔️✔️: MorganArguments(unsigned int radius = 3, bool countSimulation = false,
    /// RDKit✔️✔️:                   bool includeChirality = false, bool onlyNonzeroInvariants = false,
    /// RDKit✔️✔️:                   std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
    /// RDKit✔️✔️:                   std::uint32_t fpSize = 2048,
    /// RDKit✔️✔️:                   bool includeRedundantEnvironments = false,
    /// RDKit✔️✔️:                   bool useBondTypes = true)
    #[must_use]
    pub fn new(
        radius: u32,
        count_simulation: bool,
        include_chirality: bool,
        only_nonzero_invariants: bool,
        count_bounds: Vec<u32>,
        fp_size: u32,
        include_redundant_environments: bool,
        use_bond_types: bool,
    ) -> Result<Self, FingerprintError> {
        let fingerprint_arguments = FingerprintArguments::new(
            count_simulation,
            count_bounds,
            fp_size,
            1,
            include_chirality,
        )?;
        Ok(Self {
            fingerprint_arguments,
            df_only_nonzero_invariants: only_nonzero_invariants,
            d_radius: radius,
            df_include_redundant_environments: include_redundant_environments,
            df_use_bond_types: use_bond_types,
        })
    }

    /// RDKit✔️✔️: return "MorganArguments onlyNonzeroInvariants=" + ... + " radius=" + ...;
    #[must_use]
    pub fn infoString(&self) -> String {
        // RDKit source: MorganGenerator.cpp lines 222-227
        // RDKit✔️✔️: std::string MorganArguments::infoString() const {
        // RDKit✔️✔️:   return "MorganArguments onlyNonzeroInvariants=" +
        // RDKit✔️✔️:          std::to_string(df_onlyNonzeroInvariants) +
        // RDKit✔️✔️:          " radius=" + std::to_string(d_radius);
        // RDKit✔️✔️: }
        format!(
            "MorganArguments onlyNonzeroInvariants={} radius={}",
            self.df_only_nonzero_invariants as u8, self.d_radius
        )
    }

    /// RDKit✔️✔️: pt.put("type", "MorganArguments"); pt.put("onlyNonzeroInvariants", ...); pt.put("radius", ...); FingerprintArguments::toJSON(pt);
    #[must_use]
    pub fn toJSON(&self) -> String {
        // RDKit source: MorganGenerator.cpp lines 222-227
        // RDKit✔️✔️: void MorganArguments::toJSON(boost::property_tree::ptree &pt) const {
        // RDKit✔️✔️:   pt.put("type", "MorganArguments");
        // RDKit✔️✔️:   pt.put("onlyNonzeroInvariants", df_onlyNonzeroInvariants);
        // RDKit✔️✔️:   pt.put("radius", d_radius);
        // RDKit✔️✔️:   FingerprintArguments::toJSON(pt);
        // RDKit✔️✔️: }
        let fp = self.fingerprint_arguments.to_json();
        let fp_body = &fp[1..fp.len().saturating_sub(1)];
        let mut fields = vec![format!(
            "\"type\":\"MorganArguments\",\"onlyNonzeroInvariants\":{},\"radius\":{}",
            self.df_only_nonzero_invariants, self.d_radius
        )];
        if !fp_body.is_empty() {
            fields.push(fp_body.to_string());
        }
        format!("{{{}}}", fields.join(","))
    }

    /// RDKit✔️✔️: d_radius = pt.get<std::uint32_t>("radius", d_radius); df_onlyNonzeroInvariants = ...; FingerprintArguments::fromJSON(pt);
    pub fn fromJSON(&mut self, json: &str) -> Result<(), FingerprintError> {
        // RDKit source: MorganGenerator.cpp lines 228-234
        // RDKit✔️✔️: void MorganArguments::fromJSON(const boost::property_tree::ptree &pt) {
        // RDKit✔️✔️:   d_radius = pt.get<std::uint32_t>("radius", d_radius);
        // RDKit✔️✔️:   df_onlyNonzeroInvariants =
        // RDKit✔️✔️:       pt.get<bool>("onlyNonzeroInvariants", df_onlyNonzeroInvariants);
        // RDKit✔️✔️:   FingerprintArguments::fromJSON(pt);
        // RDKit✔️✔️: }
        if json.trim().is_empty() {
            return Ok(());
        }
        let value: Value = serde_json::from_str(json)
            .map_err(|error| FingerprintError::InvalidArgumentsJson(error.to_string()))?;
        self.from_json_value(&value)
    }

    fn from_json_value(&mut self, value: &Value) -> Result<(), FingerprintError> {
        let object = value.as_object().ok_or_else(|| {
            FingerprintError::InvalidArgumentsJson("expected JSON object".to_string())
        })?;
        if let Some(field) = object.get("radius") {
            self.d_radius = json_value_as_u32("radius", field)?;
        }
        if let Some(field) = object.get("onlyNonzeroInvariants") {
            self.df_only_nonzero_invariants = json_value_as_bool("onlyNonzeroInvariants", field)?;
        }
        self.fingerprint_arguments.from_json_value(value)?;
        Ok(())
    }
}

/// RDKit-style Morgan atom-invariant generator.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganAtomInvGenerator {
    pub include_ring_membership: bool,
}

impl MorganAtomInvGenerator {
    /// RDKit✔️✔️: MorganAtomInvGenerator(const bool includeRingMembership = true);
    #[must_use]
    pub fn new(include_ring_membership: bool) -> Self {
        Self {
            include_ring_membership,
        }
    }

    /// RDKit✔️✔️: std::vector<std::uint32_t> *MorganAtomInvGenerator::getAtomInvariants(const ROMol &mol) const
    #[must_use]
    pub fn getAtomInvariants(&self, molecule: &Molecule) -> Vec<u32> {
        // RDKit source: MorganGenerator.cpp lines 39-48
        // RDKit✔️✔️: MorganAtomInvGenerator::MorganAtomInvGenerator(const bool includeRingMembership)
        // RDKit✔️✔️:     : df_includeRingMembership(includeRingMembership) {}
        // RDKit✔️✔️: std::vector<std::uint32_t> *MorganAtomInvGenerator::getAtomInvariants(
        // RDKit✔️✔️:     const ROMol &mol) const {
        // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
        // RDKit✔️✔️:   std::unique_ptr<std::vector<std::uint32_t>> atomInvariants(
        // RDKit✔️✔️:       new std::vector<std::uint32_t>(nAtoms));
        // RDKit✔️✔️:   getConnectivityInvariants(mol, *atomInvariants, df_includeRingMembership);
        // RDKit✔️✔️:   return atomInvariants.release();
        // RDKit✔️✔️: }
        let adjacency = molecule.topology_block().adjacency.clone();
        compute_connectivity_invariants(
            molecule,
            &adjacency,
            &MorganFingerprintParams {
                include_ring_membership: self.include_ring_membership,
                ..Default::default()
            },
        )
    }

    /// RDKit✔️✔️: std::string MorganAtomInvGenerator::infoString() const
    #[must_use]
    pub fn infoString(&self) -> String {
        format!(
            "MorganInvariantGenerator includeRingMembership={}",
            self.include_ring_membership as u8
        )
    }

    /// RDKit✔️✔️: void MorganAtomInvGenerator::toJSON(boost::property_tree::ptree &pt) const
    #[must_use]
    pub fn toJSON(&self) -> String {
        format!(
            "{{\"type\":\"MorganAtomInvGenerator\",\"includeRingMembership\":{}}}",
            self.include_ring_membership
        )
    }

    /// RDKit✔️✔️: void MorganAtomInvGenerator::fromJSON(const boost::property_tree::ptree &pt)
    pub fn fromJSON(&mut self, json: &str) -> Result<(), FingerprintError> {
        if json.trim().is_empty() {
            return Ok(());
        }
        let value: Value = serde_json::from_str(json)
            .map_err(|error| FingerprintError::InvalidArgumentsJson(error.to_string()))?;
        let object = value.as_object().ok_or_else(|| {
            FingerprintError::InvalidArgumentsJson("expected JSON object".to_string())
        })?;
        if let Some(field) = object.get("includeRingMembership") {
            self.include_ring_membership = json_value_as_bool("includeRingMembership", field)?;
        }
        Ok(())
    }

    /// RDKit✔️✔️: MorganAtomInvGenerator *clone() const
    #[must_use]
    pub fn clone(&self) -> Self {
        Self {
            include_ring_membership: self.include_ring_membership,
        }
    }
}

/// RDKit-style Morgan feature atom-invariant generator.
///
/// `patterns == None` preserves RDKit's `dp_patterns == nullptr` branch and
/// uses the default Gobbi/Poppinger feature SMARTS. `Some(patterns)` preserves
/// RDKit's supplied-pattern branch and uses only those queries.
#[derive(Debug, PartialEq)]
pub struct MorganFeatureAtomInvGenerator {
    patterns: Option<Vec<SsMatcher>>,
    pattern_smarts: Option<Vec<String>>,
}

impl MorganFeatureAtomInvGenerator {
    /// RDKit✔️✔️: MorganFeatureAtomInvGenerator(const std::vector<const ROMol *> *patterns = nullptr);
    #[must_use]
    pub fn new() -> Self {
        // RDKit source: MorganGenerator.cpp lines 66-75
        // RDKit✔️✔️: MorganFeatureAtomInvGenerator::MorganFeatureAtomInvGenerator(
        // RDKit✔️✔️:     const std::vector<const ROMol *> *patterns) {
        // RDKit✔️✔️:   if (patterns) {
        // RDKit✔️✔️:     dp_patterns = new std::vector<const ROMol *>;
        // RDKit✔️✔️:     dp_patterns->reserve(patterns->size());
        // RDKit✔️✔️:     for (auto pattern : *patterns) {
        // RDKit✔️✔️:       dp_patterns->push_back(new ROMol(*pattern));
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        Self {
            patterns: None,
            pattern_smarts: None,
        }
    }

    /// Constructing supplied feature SMARTS from strings is not source-ported.
    pub fn from_smarts_patterns(patterns: &[&str]) -> Result<Self, FingerprintError> {
        if patterns.is_empty() {
            return Ok(Self::new());
        }
        Err(FingerprintError::UnsupportedOption {
            option: "MorganFeatureAtomInvGenerator.patternSMARTS",
            reason: "RDKit SmartsToMol supplied-pattern branch is not source-ported",
        })
    }

    /// RDKit✔️✔️: void MorganFeatureAtomInvGenerator::cleanUpPatterns()
    fn clean_up_patterns(&mut self) {
        // RDKit source: MorganGenerator.cpp lines 77-85
        // RDKit✔️✔️: void MorganFeatureAtomInvGenerator::cleanUpPatterns() {
        // RDKit✔️✔️:   if (dp_patterns) {
        // RDKit✔️✔️:     for (auto mol : *dp_patterns) {
        // RDKit✔️✔️:       delete mol;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     delete dp_patterns;
        // RDKit✔️✔️:     dp_patterns = nullptr;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        self.patterns = None;
        self.pattern_smarts = None;
    }

    /// RDKit✔️✔️: std::vector<std::uint32_t> *MorganFeatureAtomInvGenerator::getAtomInvariants(const ROMol &mol) const
    #[must_use]
    #[allow(non_snake_case)]
    pub fn getAtomInvariants(&self, molecule: &Molecule) -> Result<Vec<u32>, FingerprintError> {
        // RDKit source: MorganGenerator.cpp lines 132-140
        // RDKit✔️✔️: std::vector<std::uint32_t> *MorganFeatureAtomInvGenerator::getAtomInvariants(
        // RDKit✔️✔️:     const ROMol &mol) const {
        // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
        // RDKit✔️✔️:   std::vector<std::uint32_t> *result = new std::vector<std::uint32_t>(nAtoms);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   getFeatureInvariants(mol, *result, dp_patterns);
        // RDKit✔️✔️:   return result;
        // RDKit✔️✔️: }
        compute_feature_invariants_with_patterns(molecule, self.patterns.as_deref())
    }

    /// RDKit✔️✔️: std::string MorganFeatureAtomInvGenerator::infoString() const
    #[must_use]
    #[allow(non_snake_case)]
    pub fn infoString(&self) -> String {
        // RDKit source: MorganGenerator.cpp lines 91-93
        // RDKit✔️✔️: std::string MorganFeatureAtomInvGenerator::infoString() const {
        // RDKit✔️✔️:   return "MorganFeatureInvariantGenerator";
        // RDKit✔️✔️: }
        "MorganFeatureInvariantGenerator".to_string()
    }

    /// RDKit✔️✔️: void MorganFeatureAtomInvGenerator::toJSON(boost::property_tree::ptree &pt) const
    #[must_use]
    #[allow(non_snake_case)]
    pub fn toJSON(&self) -> String {
        // RDKit source: MorganGenerator.cpp lines 94-108
        // RDKit✔️✔️: void MorganFeatureAtomInvGenerator::toJSON(
        // RDKit✔️✔️:     boost::property_tree::ptree &pt) const {
        // RDKit✔️✔️:   pt.put("type", "MorganFeatureAtomInvGenerator");
        // RDKit✔️✔️:   if (dp_patterns) {
        // RDKit✔️✔️:     boost::property_tree::ptree patternsNode;
        // RDKit✔️✔️:     for (const auto &pattern : *dp_patterns) {
        // RDKit✔️✔️:       boost::property_tree::ptree patternNode;
        // RDKit✔️✔️:       std::string smarts = MolToSmarts(*pattern);
        // RDKit✔️✔️:       patternNode.put("", smarts);
        // RDKit✔️✔️:       patternsNode.push_back(std::make_pair("", patternNode));
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     pt.add_child("patternSMARTS", patternsNode);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   AtomInvariantsGenerator::toJSON(pt);
        // RDKit✔️✔️: }
        let mut fields = vec![r#""type":"MorganFeatureAtomInvGenerator""#.to_string()];
        if let Some(pattern_smarts) = &self.pattern_smarts {
            let encoded = pattern_smarts
                .iter()
                .map(|pattern| Value::String(pattern.clone()).to_string())
                .collect::<Vec<_>>()
                .join(",");
            fields.push(format!(r#""patternSMARTS":[{encoded}]"#));
        }
        format!("{{{}}}", fields.join(","))
    }

    /// RDKit✔️✔️: void MorganFeatureAtomInvGenerator::fromJSON(const boost::property_tree::ptree &pt)
    #[allow(non_snake_case)]
    pub fn fromJSON(&mut self, json: &str) -> Result<(), FingerprintError> {
        // RDKit source: MorganGenerator.cpp lines 109-124
        // RDKit✔️✔️: void MorganFeatureAtomInvGenerator::fromJSON(
        // RDKit✔️✔️:     const boost::property_tree::ptree &pt) {
        // RDKit✔️✔️:   if (pt.get_child_optional("patternSMARTS")) {
        // RDKit✔️✔️:     const auto &patternsNode = pt.get_child("patternSMARTS");
        // RDKit✔️✔️:     cleanUpPatterns();
        // RDKit✔️✔️:     dp_patterns = new std::vector<const ROMol *>();
        // RDKit✔️✔️:     for (const auto &patternNode : patternsNode) {
        // RDKit✔️✔️:       std::string smarts = patternNode.second.get_value<std::string>();
        // RDKit❌❌:       ROMol *patternMol = SmartsToMol(smarts);
        // RDKit❌❌:       if (patternMol) {
        // RDKit✔️✔️:         dp_patterns->push_back(patternMol);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   AtomInvariantsGenerator::fromJSON(pt);
        // RDKit✔️✔️: }
        if json.trim().is_empty() {
            return Ok(());
        }
        let value: Value = serde_json::from_str(json)
            .map_err(|error| FingerprintError::InvalidArgumentsJson(error.to_string()))?;
        let object = value.as_object().ok_or_else(|| {
            FingerprintError::InvalidArgumentsJson("expected JSON object".to_string())
        })?;
        let Some(field) = object.get("patternSMARTS") else {
            return Ok(());
        };
        let patterns = field.as_array().ok_or_else(|| {
            FingerprintError::InvalidArgumentsJson("patternSMARTS must be an array".to_string())
        })?;
        let mut pattern_refs = Vec::with_capacity(patterns.len());
        for pattern in patterns {
            let pattern = pattern.as_str().ok_or_else(|| {
                FingerprintError::InvalidArgumentsJson(
                    "patternSMARTS entries must be strings".to_string(),
                )
            })?;
            pattern_refs.push(pattern);
        }
        Self::from_smarts_patterns(&pattern_refs).map(|updated| {
            *self = updated;
        })
    }

    /// RDKit✔️✔️: MorganFeatureAtomInvGenerator *clone() const
    #[must_use]
    pub fn clone(&self) -> Self {
        // RDKit source: MorganGenerator.cpp lines 126-129
        // RDKit✔️✔️: MorganFeatureAtomInvGenerator *MorganFeatureAtomInvGenerator::clone() const {
        // RDKit✔️✔️:   return new MorganFeatureAtomInvGenerator(dp_patterns);
        // RDKit✔️✔️: }
        Clone::clone(self)
    }
}

impl Clone for MorganFeatureAtomInvGenerator {
    fn clone(&self) -> Self {
        Self {
            patterns: self.patterns.clone(),
            pattern_smarts: self.pattern_smarts.clone(),
        }
    }
}

impl Default for MorganFeatureAtomInvGenerator {
    fn default() -> Self {
        Self::new()
    }
}

impl Drop for MorganFeatureAtomInvGenerator {
    fn drop(&mut self) {
        // RDKit source: MorganGenerator.cpp lines 87-89
        // RDKit✔️✔️: MorganFeatureAtomInvGenerator::~MorganFeatureAtomInvGenerator() {
        // RDKit✔️✔️:   cleanUpPatterns();
        // RDKit✔️✔️: }
        self.clean_up_patterns();
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MorganAtomInvariantsGenerator {
    Connectivity { include_ring_membership: bool },
    Feature,
}

impl MorganAtomInvariantsGenerator {
    fn get_atom_invariants(&self, molecule: &Molecule) -> Result<Vec<u32>, FingerprintError> {
        match self {
            Self::Connectivity {
                include_ring_membership,
            } => {
                Ok(MorganAtomInvGenerator::new(*include_ring_membership)
                    .getAtomInvariants(molecule))
            }
            Self::Feature => MorganFeatureAtomInvGenerator::new().getAtomInvariants(molecule),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganBondInvariantsGenerator {
    pub use_bond_types: bool,
    pub use_chirality: bool,
}

pub type MorganBondInvGenerator = MorganBondInvariantsGenerator;

fn rdkit_use_legacy_stereo_perception() -> bool {
    // RDKit✔️✔️: constexpr bool useLegacyStereoDefaultVal =
    // RDKit✔️✔️:     true;  //!< whether or not the legacy stereo perception code is used by
    // RDKit✔️✔️:            //!< default
    // RDKit✔️✔️: bool getUseLegacyStereoPerception() {
    // RDKit✔️✔️:   return getValFromEnvironment(useLegacyStereoEnvVar,
    // RDKit✔️✔️:                                useLegacyStereoDefaultVal);
    // RDKit✔️✔️: }
    true
}

impl MorganBondInvariantsGenerator {
    /// RDKit✔️✔️: MorganBondInvGenerator(const bool useBondTypes = true, const bool useChirality = false);
    #[must_use]
    pub fn new(use_bond_types: bool, use_chirality: bool) -> Self {
        // RDKit source: MorganGenerator.cpp lines 142-144
        // RDKit✔️✔️: MorganBondInvGenerator::MorganBondInvGenerator(const bool useBondTypes,
        // RDKit✔️✔️:                                                const bool useChirality)
        // RDKit✔️✔️:     : df_useBondTypes(useBondTypes), df_useChirality(useChirality) {}
        Self {
            use_bond_types,
            use_chirality,
        }
    }

    /// RDKit✔️✔️: std::vector<std::uint32_t> *MorganBondInvGenerator::getBondInvariants(const ROMol &mol) const
    #[must_use]
    #[allow(non_snake_case)]
    pub fn getBondInvariants(&self, molecule: &Molecule) -> Result<Vec<u32>, FingerprintError> {
        self.try_get_bond_invariants(molecule)
    }

    fn try_get_bond_invariants(&self, molecule: &Molecule) -> Result<Vec<u32>, FingerprintError> {
        // RDKit source: MorganGenerator.cpp lines 146-183
        // RDKit✔️✔️: std::vector<std::uint32_t> *MorganBondInvGenerator::getBondInvariants(
        // RDKit✔️✔️:     const ROMol &mol) const {
        // RDKit✔️✔️:   std::vector<std::uint32_t> *result =
        // RDKit✔️✔️:       new std::vector<std::uint32_t>(mol.getNumBonds());
        // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
        // RDKit✔️✔️:     Bond const *bond = mol.getBondWithIdx(i);
        // RDKit✔️✔️:     int32_t bondInvariant = 1;
        // RDKit✔️✔️:     if (df_useBondTypes) {
        // RDKit✔️✔️:       if (!df_useChirality || bond->getBondType() != Bond::DOUBLE ||
        // RDKit✔️✔️:           bond->getStereo() == Bond::STEREONONE) {
        // RDKit✔️✔️:         bondInvariant = static_cast<int32_t>(bond->getBondType());
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         auto bondStereo = static_cast<int32_t>(bond->getStereo());
        // RDKit✔️✔️:         if (!Chirality::getUseLegacyStereoPerception()) {
        // RDKit✔️✔️:           if (!mol.hasProp(common_properties::_CIPComputed)) {
        // RDKit✔️✔️:             CIPLabeler::assignCIPLabels(const_cast<ROMol &>(mol));
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:           std::string cipCode;
        // RDKit✔️✔️:           if (bond->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
        // RDKit✔️✔️:             if (cipCode == "E") {
        // RDKit✔️✔️:               bondStereo = static_cast<int32_t>(Bond::STEREOE);
        // RDKit✔️✔️:             } else if (cipCode == "Z") {
        // RDKit✔️✔️:               bondStereo = static_cast<int32_t>(Bond::STEREOZ);
        // RDKit✔️✔️:             }
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:         const int32_t stereoOffset = 100;
        // RDKit✔️✔️:         const int32_t bondTypeOffset = 10;
        // RDKit✔️✔️:         bondInvariant =
        // RDKit✔️✔️:             stereoOffset +
        // RDKit✔️✔️:             bondTypeOffset * static_cast<int32_t>(bond->getBondType()) +
        // RDKit✔️✔️:             bondStereo;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     (*result)[bond->getIdx()] = static_cast<int32_t>(bondInvariant);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return result;
        // RDKit✔️✔️: }
        let needs_cip_labels = self.use_bond_types
            && self.use_chirality
            && molecule.bonds().iter().any(|bond| {
                bond.order() == crate::BondOrder::Double && bond.stereo() != crate::BondStereo::None
            });
        let labeled;
        let molecule = if needs_cip_labels
            && !rdkit_use_legacy_stereo_perception()
            && molecule.prop("_CIPComputed").is_none()
        {
            labeled = assign_cip_labels(molecule, 0)?;
            &labeled
        } else {
            molecule
        };
        Ok(molecule
            .bonds()
            .iter()
            .map(|bond| self.bond_invariant(bond))
            .collect())
    }

    fn bond_invariant(&self, bond: &crate::Bond) -> u32 {
        let mut bond_invariant = 1u32;
        if self.use_bond_types {
            let bond_type = rdkit_bond_type_code(bond.order());
            if !self.use_chirality
                || bond.order() != crate::BondOrder::Double
                || bond.stereo() == crate::BondStereo::None
            {
                bond_invariant = bond_type;
            } else {
                let bond_stereo = match bond.prop("_CIPCode") {
                    Some("E") => crate::BondStereo::E,
                    Some("Z") => crate::BondStereo::Z,
                    _ => bond.stereo(),
                };
                let stereo_offset = 100u32;
                let bond_type_offset = 10u32;
                bond_invariant = stereo_offset
                    + bond_type_offset * bond_type
                    + rdkit_bond_stereo_code(bond_stereo);
            }
        }
        bond_invariant
    }

    /// RDKit✔️✔️: std::string MorganBondInvGenerator::infoString() const
    #[must_use]
    #[allow(non_snake_case)]
    pub fn infoString(&self) -> String {
        // RDKit source: MorganGenerator.cpp lines 185-189
        // RDKit✔️✔️: std::string MorganBondInvGenerator::infoString() const {
        // RDKit✔️✔️:   return "MorganInvariantGenerator useBondTypes=" +
        // RDKit✔️✔️:          std::to_string(df_useBondTypes) +
        // RDKit✔️✔️:          " useChirality=" + std::to_string(df_useChirality);
        // RDKit✔️✔️: }
        format!(
            "MorganInvariantGenerator useBondTypes={} useChirality={}",
            self.use_bond_types as u8, self.use_chirality as u8
        )
    }

    /// RDKit✔️✔️: void MorganBondInvGenerator::toJSON(boost::property_tree::ptree &pt) const
    #[must_use]
    #[allow(non_snake_case)]
    pub fn toJSON(&self) -> String {
        // RDKit source: MorganGenerator.cpp lines 190-195
        // RDKit✔️✔️: void MorganBondInvGenerator::toJSON(boost::property_tree::ptree &pt) const {
        // RDKit✔️✔️:   pt.put("type", "MorganBondInvGenerator");
        // RDKit✔️✔️:   pt.put("useBondTypes", df_useBondTypes);
        // RDKit✔️✔️:   pt.put("useChirality", df_useChirality);
        // RDKit✔️✔️:   BondInvariantsGenerator::toJSON(pt);
        // RDKit✔️✔️: }
        format!(
            "{{\"type\":\"MorganBondInvGenerator\",\"useBondTypes\":{},\"useChirality\":{}}}",
            self.use_bond_types, self.use_chirality
        )
    }

    /// RDKit✔️✔️: void MorganBondInvGenerator::fromJSON(const boost::property_tree::ptree &pt)
    #[allow(non_snake_case)]
    pub fn fromJSON(&mut self, json: &str) -> Result<(), FingerprintError> {
        // RDKit source: MorganGenerator.cpp lines 196-200
        // RDKit✔️✔️: void MorganBondInvGenerator::fromJSON(const boost::property_tree::ptree &pt) {
        // RDKit✔️✔️:   df_useBondTypes = pt.get<bool>("useBondTypes", df_useBondTypes);
        // RDKit✔️✔️:   df_useChirality = pt.get<bool>("useChirality", df_useChirality);
        // RDKit✔️✔️:   BondInvariantsGenerator::fromJSON(pt);
        // RDKit✔️✔️: }
        if json.trim().is_empty() {
            return Ok(());
        }
        let value: Value = serde_json::from_str(json)
            .map_err(|error| FingerprintError::InvalidArgumentsJson(error.to_string()))?;
        let object = value.as_object().ok_or_else(|| {
            FingerprintError::InvalidArgumentsJson("expected JSON object".to_string())
        })?;
        if let Some(field) = object.get("useBondTypes") {
            self.use_bond_types = json_value_as_bool("useBondTypes", field)?;
        }
        if let Some(field) = object.get("useChirality") {
            self.use_chirality = json_value_as_bool("useChirality", field)?;
        }
        Ok(())
    }

    /// RDKit✔️✔️: MorganBondInvGenerator *clone() const
    #[must_use]
    pub fn clone(&self) -> Self {
        // RDKit source: MorganGenerator.cpp lines 202-204
        // RDKit✔️✔️: MorganBondInvGenerator *MorganBondInvGenerator::clone() const {
        // RDKit✔️✔️:   return new MorganBondInvGenerator(df_useBondTypes, df_useChirality);
        // RDKit✔️✔️: }
        Self::new(self.use_bond_types, self.use_chirality)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct MorganAdditionalOutput {
    pub atom_counts: Vec<u32>,
    pub atom_to_bits: Vec<Vec<usize>>,
    pub bit_info_map: BTreeMap<usize, Vec<(usize, u32)>>,
    pub atoms_per_bit: BTreeMap<usize, Vec<Vec<usize>>>,
}

/// RDKit-style Morgan atom environment.
pub struct MorganAtomEnv {
    code: u64,
    atom_id: usize,
    layer: u32,
    atom_count: usize,
    adjacency: AdjacencyList,
}

impl MorganAtomEnv {
    /// RDKit✔️✔️: MorganAtomEnv(const std::uint32_t code, const unsigned int atomId,
    /// RDKit✔️✔️:               const unsigned int layer, const ROMol *mol)
    #[must_use]
    pub fn new(code: u64, atom_id: usize, layer: u32, molecule: &Molecule) -> Self {
        Self {
            code,
            atom_id,
            layer,
            atom_count: molecule.num_atoms(),
            adjacency: molecule.topology_block().adjacency.clone(),
        }
    }

    // RDKit✔️✔️: OutputType MorganAtomEnv<OutputType>::getBitId(...) const
    #[must_use]
    #[allow(non_snake_case)]
    pub fn getBitId(&self) -> u64 {
        // RDKit source: MorganGenerator.cpp lines 267-277
        // RDKit✔️✔️: OutputType MorganAtomEnv<OutputType>::getBitId(
        // RDKit✔️✔️:     FingerprintArguments *,              // arguments
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *,  // atomInvariants
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *,  // bondInvariants
        // RDKit✔️✔️:     AdditionalOutput *,                  // additional Output
        // RDKit✔️✔️:     const bool,                          // hashResults
        // RDKit✔️✔️:     const std::uint64_t                  // fpSize
        // RDKit✔️✔️: ) const {
        // RDKit✔️✔️:   return d_code;
        // RDKit✔️✔️: }
        self.code
    }

    // RDKit✔️❌: void MorganAtomEnv<OutputType>::updateAdditionalOutput(AdditionalOutput *additionalOutput, size_t bitId) const
    #[allow(non_snake_case)]
    pub fn updateAdditionalOutput(&self, additional_output: &mut AdditionalOutput, bit_id: u64) {
        // RDKit source: MorganGenerator.cpp lines 236-264
        // RDKit✔️❌: void MorganAtomEnv<OutputType>::updateAdditionalOutput(
        // RDKit✔️❌:     AdditionalOutput *additionalOutput, size_t bitId) const {
        // RDKit✔️❌:   PRECONDITION(additionalOutput, "bad output pointer");
        // RDKit✔️❌:   PRECONDITION(d_mol, "bad mol pointer");
        // RDKit✔️❌:   if (additionalOutput->bitInfoMap) {
        // RDKit✔️❌:     (*additionalOutput->bitInfoMap)[bitId].emplace_back(d_atomId, d_layer);
        // RDKit✔️❌:   }
        // RDKit✔️❌:   if (additionalOutput->atomCounts) {
        // RDKit✔️❌:     (*additionalOutput->atomCounts)[d_atomId]++;
        // RDKit✔️❌:   }
        // RDKit✔️❌:   if (additionalOutput->atomToBits) {
        // RDKit✔️❌:     (*additionalOutput->atomToBits)[d_atomId].push_back(bitId);
        // RDKit✔️❌:   }
        // RDKit✔️❌:   if (additionalOutput->atomsPerBit) {
        // RDKit✔️❌:     std::vector<int> atomsInvolved;
        // RDKit✔️❌:     atomsInvolved.push_back(d_atomId);
        // RDKit✔️❌:     if (d_layer > 0) {
        // RDKit✔️❌:       const auto dm = MolOps::getDistanceMat(*d_mol);
        // RDKit✔️❌:       for (unsigned int i = 0; i < d_mol->getNumAtoms(); ++i) {
        // RDKit✔️❌:         if (static_cast<unsigned int>(dm[d_atomId * d_mol->getNumAtoms() + i] +
        // RDKit✔️❌:                                       .1) <= d_layer &&
        // RDKit✔️❌:             i != d_atomId) {
        // RDKit✔️❌:           atomsInvolved.push_back(i);
        // RDKit✔️❌:         }
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:     (*additionalOutput->atomsPerBit)[bitId].push_back(std::move(atomsInvolved));
        // RDKit✔️❌:   }
        // RDKit✔️❌: }
        //
        // COSMolKit uses bounded BFS from the center atom instead of building
        // RDKit's full distance matrix for each atomsPerBit update.
        if let Some(bit_info_map) = additional_output.bit_info_map.as_mut() {
            bit_info_map
                .entry(bit_id)
                .or_default()
                .push((self.atom_id as u32, self.layer));
        }
        if let Some(atom_counts) = additional_output.atom_counts.as_mut() {
            atom_counts[self.atom_id] += 1;
        }
        if let Some(atom_to_bits) = additional_output.atom_to_bits.as_mut() {
            atom_to_bits[self.atom_id].push(bit_id);
        }
        if let Some(atoms_per_bit) = additional_output.atoms_per_bit.as_mut() {
            atoms_per_bit
                .entry(bit_id)
                .or_default()
                .push(self.atoms_involved());
        }
    }

    fn atoms_involved(&self) -> Vec<usize> {
        let mut atoms = vec![self.atom_id];
        if self.layer == 0 {
            return atoms;
        }

        let mut distances = vec![u32::MAX; self.atom_count];
        let mut queue = VecDeque::new();
        distances[self.atom_id] = 0;
        queue.push_back(self.atom_id);
        while let Some(current) = queue.pop_front() {
            let next_distance = distances[current] + 1;
            if next_distance > self.layer {
                continue;
            }
            for neighbor in self.adjacency.neighbors_of(current) {
                if distances[neighbor.atom_index] != u32::MAX {
                    continue;
                }
                distances[neighbor.atom_index] = next_distance;
                queue.push_back(neighbor.atom_index);
            }
        }

        for (atom_idx, distance) in distances.into_iter().enumerate() {
            if atom_idx != self.atom_id && distance <= self.layer {
                atoms.push(atom_idx);
            }
        }
        atoms
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct MorganBondEnvironment {
    bit_count: usize,
    blocks: Vec<usize>,
}

impl MorganBondEnvironment {
    fn new(bit_count: usize) -> Self {
        let block_count = bit_count.div_ceil(usize::BITS as usize);
        Self {
            bit_count,
            blocks: vec![0; block_count],
        }
    }

    fn set(&mut self, bit: usize) {
        if bit >= self.bit_count {
            return;
        }
        let block = bit / usize::BITS as usize;
        let offset = bit % usize::BITS as usize;
        self.blocks[block] |= 1usize << offset;
    }

    fn union_with(&mut self, other: &Self) {
        for (block, other_block) in self.blocks.iter_mut().zip(&other.blocks) {
            *block |= *other_block;
        }
    }
}

impl Ord for MorganBondEnvironment {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Boost dynamic_bitset source:
        //   return (a.m_num_bits == b.m_num_bits) && (a.m_bits == b.m_bits);
        //   for (size_type ii = a.num_blocks(); ii > 0; --ii) {
        //     size_type i = ii-1;
        //     if (a.m_bits[i] < b.m_bits[i]) return true;
        //     else if (a.m_bits[i] > b.m_bits[i]) return false;
        //   }
        // RDKit✔️✔️: Morgan AccumTuple ordering compares boost::dynamic_bitset
        // RDKit✔️✔️: environments before the generated code and atom index.
        match self.bit_count.cmp(&other.bit_count) {
            std::cmp::Ordering::Equal => {
                for (left, right) in self.blocks.iter().zip(&other.blocks).rev() {
                    match left.cmp(right) {
                        std::cmp::Ordering::Equal => continue,
                        ordering => return ordering,
                    }
                }
                std::cmp::Ordering::Equal
            }
            ordering => ordering,
        }
    }
}

impl PartialOrd for MorganBondEnvironment {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// RDKit-style Morgan atom-environment generator.
#[derive(Debug, Clone, Copy, Default)]
pub struct MorganEnvGenerator;

impl MorganEnvGenerator {
    #[must_use]
    pub fn new() -> Self {
        Self
    }

    // RDKit✔️✔️: std::string MorganEnvGenerator<OutputType>::infoString() const
    #[must_use]
    #[allow(non_snake_case)]
    pub fn infoString(&self) -> String {
        // RDKit source: MorganGenerator.cpp lines 461-464
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: std::string MorganEnvGenerator<OutputType>::infoString() const {
        // RDKit✔️✔️:   return "MorganEnvironmentGenerator";
        // RDKit✔️✔️: }
        "MorganEnvironmentGenerator".to_string()
    }

    // RDKit✔️✔️: void MorganEnvGenerator<OutputType>::toJSON(boost::property_tree::ptree &pt) const
    #[must_use]
    #[allow(non_snake_case)]
    pub fn toJSON(&self) -> String {
        // RDKit source: MorganGenerator.cpp lines 465-471
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: void MorganEnvGenerator<OutputType>::toJSON(
        // RDKit✔️✔️:     boost::property_tree::ptree &pt) const {
        // RDKit✔️✔️:   pt.put("type", "MorganEnvGenerator");
        // RDKit✔️✔️:   AtomEnvironmentGenerator<OutputType>::toJSON(pt);
        // RDKit✔️✔️: }
        r#"{"type":"MorganEnvGenerator"}"#.to_string()
    }

    // RDKit✔️✔️: void MorganEnvGenerator<OutputType>::fromJSON(const boost::property_tree::ptree &pt)
    #[allow(non_snake_case)]
    pub fn fromJSON(&mut self, json: &str) -> Result<(), FingerprintError> {
        // RDKit source: MorganGenerator.cpp lines 473-477
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: void MorganEnvGenerator<OutputType>::fromJSON(
        // RDKit✔️✔️:     const boost::property_tree::ptree &pt) {
        // RDKit✔️✔️:   AtomEnvironmentGenerator<OutputType>::fromJSON(pt);
        // RDKit✔️✔️: }
        validate_stateless_environment_generator_json(json)
    }

    // RDKit✔️✔️: OutputType MorganEnvGenerator<OutputType>::getResultSize() const
    #[must_use]
    #[allow(non_snake_case)]
    pub fn getResultSize(&self) -> u64 {
        // RDKit source: MorganGenerator.cpp lines 212-215
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: OutputType MorganEnvGenerator<OutputType>::getResultSize() const {
        // RDKit✔️✔️:   return std::numeric_limits<OutputType>::max();
        // RDKit✔️✔️: }
        u64::from(u32::MAX)
    }

    // RDKit✔️❌: std::vector<AtomEnvironment<OutputType> *> MorganEnvGenerator<OutputType>::getEnvironments(...)
    #[allow(non_snake_case)]
    pub fn getEnvironments(
        &self,
        molecule: &Molecule,
        arguments: &MorganArguments,
        from_atoms: Option<&[usize]>,
        _ignore_atoms: Option<&[usize]>,
        atom_invariants: &[u32],
        bond_invariants: &[u32],
    ) -> Result<Vec<MorganAtomEnv>, FingerprintError> {
        // RDKit source: MorganGenerator.cpp lines 281-459
        // RDKit✔️✔️: std::vector<AtomEnvironment<OutputType> *>
        // RDKit✔️✔️: MorganEnvGenerator<OutputType>::getEnvironments(
        // RDKit✔️✔️:     const ROMol &mol, FingerprintArguments *arguments,
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *fromAtoms,
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *,  // ignoreAtoms
        // RDKit✔️✔️:     const int,                           // confId
        // RDKit✔️✔️:     const AdditionalOutput *,            // additionalOutput
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *atomInvariants,
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *bondInvariants,
        // RDKit✔️✔️:     const bool  // hashResults
        // RDKit✔️✔️: ) const {
        // RDKit✔️✔️:   PRECONDITION(atomInvariants && (atomInvariants->size() >= mol.getNumAtoms()),
        // RDKit✔️✔️:                "bad atom invariants size");
        // RDKit✔️✔️:   PRECONDITION(bondInvariants && (bondInvariants->size() >= mol.getNumBonds()),
        // RDKit✔️✔️:                "bad bond invariants size");
        if atom_invariants.len() < molecule.num_atoms() {
            return Err(FingerprintError::InvalidArguments {
                reason: "bad atom invariants size",
            });
        }
        if bond_invariants.len() < molecule.num_bonds() {
            return Err(FingerprintError::InvalidArguments {
                reason: "bad bond invariants size",
            });
        }
        if from_atoms
            .into_iter()
            .flatten()
            .any(|&atom_idx| atom_idx >= molecule.num_atoms())
        {
            return Err(FingerprintError::InvalidArguments {
                reason: "fromAtoms contains an atom index outside the molecule",
            });
        }

        // RDKit✔️✔️:   auto *morganArguments = dynamic_cast<MorganArguments *>(arguments);
        // RDKit✔️✔️:   PRECONDITION(morganArguments, "bad arguments type");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
        // RDKit✔️✔️:   const unsigned int maxNumResults = (morganArguments->d_radius + 1) * nAtoms;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   std::vector<AtomEnvironment<OutputType> *> result =
        // RDKit✔️✔️:       std::vector<AtomEnvironment<OutputType> *>();
        // RDKit✔️✔️:   result.reserve(maxNumResults);
        let n_atoms = molecule.num_atoms();
        let max_num_results = (arguments.d_radius as usize + 1) * n_atoms;
        let mut result = Vec::with_capacity(max_num_results);

        // RDKit✔️✔️:   // if we are using chirality, we need to make sure the atoms have R/S labels
        // RDKit✔️✔️:   if (morganArguments->df_includeChirality &&
        // RDKit✔️✔️:       !Chirality::getUseLegacyStereoPerception() &&
        // RDKit✔️✔️:       !mol.hasProp(common_properties::_CIPComputed)) {
        // RDKit✔️✔️:     CIPLabeler::assignCIPLabels(const_cast<ROMol &>(mol));
        // RDKit✔️✔️:   }
        let labeled;
        let molecule = if arguments.fingerprint_arguments.df_include_chirality
            && !rdkit_use_legacy_stereo_perception()
            && molecule.prop("_CIPComputed").is_none()
        {
            labeled = assign_cip_labels(molecule, 0)?;
            &labeled
        } else {
            molecule
        };

        // RDKit✔️✔️:   std::vector<OutputType> currentInvariants(atomInvariants->size());
        // RDKit✔️✔️:   std::copy(atomInvariants->begin(), atomInvariants->end(),
        // RDKit✔️✔️:             currentInvariants.begin());
        // RDKit✔️✔️:   std::vector<OutputType> nextLayerInvariants(nAtoms);
        // RDKit✔️✔️:   std::vector<std::pair<int32_t, uint32_t>> neighborhoodInvariants;
        // RDKit✔️✔️:   neighborhoodInvariants.reserve(8);
        let mut current_invariants = atom_invariants.to_vec();
        let mut next_layer_invariants = vec![0u32; n_atoms];
        let mut neighborhood_invariants: Vec<(i32, u32)> = Vec::with_capacity(8);

        // RDKit✔️✔️:   boost::dynamic_bitset<> includeAtoms(nAtoms);
        // RDKit✔️✔️:   if (fromAtoms) {
        // RDKit✔️✔️:     for (auto idx : *fromAtoms) {
        // RDKit✔️✔️:       includeAtoms.set(idx, 1);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     includeAtoms.set();
        // RDKit✔️✔️:   }
        let mut include_atoms = vec![from_atoms.is_none(); n_atoms];
        if let Some(from_atoms) = from_atoms {
            for &atom_idx in from_atoms {
                include_atoms[atom_idx] = true;
            }
        }

        // RDKit✔️❌:   boost::dynamic_bitset<> chiralAtoms(nAtoms);
        // RDKit✔️❌:   std::unordered_set<boost::dynamic_bitset<>> neighborhoods;
        // RDKit✔️❌:   neighborhoods.reserve(maxNumResults);
        // RDKit✔️❌:   std::vector<boost::dynamic_bitset<>> atomNeighborhoods(
        // RDKit✔️❌:       nAtoms, boost::dynamic_bitset<>(mol.getNumBonds()));
        // RDKit✔️❌:   std::vector<boost::dynamic_bitset<>> roundAtomNeighborhoods =
        // RDKit✔️❌:       atomNeighborhoods;
        // RDKit✔️❌:   boost::dynamic_bitset<> deadAtoms(nAtoms);
        //
        // `MorganBondEnvironment` is a sparse set of bond indices with an
        // Ord implementation matching boost::dynamic_bitset's high-bit-first
        // comparison for equal-size bitsets. It avoids allocating a dense
        // bitset per atom while preserving the source ordering and equality.
        let mut chiral_atoms = vec![false; n_atoms];
        let mut neighborhoods: HashSet<MorganBondEnvironment> =
            HashSet::with_capacity(max_num_results);
        let empty_environment = MorganBondEnvironment::new(molecule.num_bonds());
        let mut atom_neighborhoods = vec![empty_environment.clone(); n_atoms];
        let mut round_atom_neighborhoods = atom_neighborhoods.clone();
        let mut dead_atoms = vec![false; n_atoms];

        // RDKit✔️✔️:   std::vector<unsigned int> atomOrder(nAtoms);
        // RDKit✔️✔️:   if (morganArguments->df_onlyNonzeroInvariants) {
        // RDKit✔️✔️:     std::vector<std::pair<int32_t, uint32_t>> ordering;
        // RDKit✔️✔️:     for (unsigned int i = 0; i < nAtoms; ++i) {
        // RDKit✔️✔️:       if (!currentInvariants[i]) {
        // RDKit✔️✔️:         ordering.emplace_back(1, i);
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         ordering.emplace_back(0, i);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     std::sort(ordering.begin(), ordering.end());
        // RDKit✔️✔️:     for (unsigned int i = 0; i < nAtoms; ++i) {
        // RDKit✔️✔️:       atomOrder[i] = ordering[i].second;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     for (unsigned int i = 0; i < nAtoms; ++i) {
        // RDKit✔️✔️:       atomOrder[i] = i;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        let mut atom_order: Vec<usize> = (0..n_atoms).collect();
        if arguments.df_only_nonzero_invariants {
            atom_order.sort_by_key(|&atom_idx| {
                (
                    if current_invariants[atom_idx] == 0 {
                        1
                    } else {
                        0
                    },
                    atom_idx,
                )
            });
        }

        // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
        // RDKit✔️✔️:     if (includeAtoms[i]) {
        // RDKit✔️✔️:       if (!morganArguments->df_onlyNonzeroInvariants || currentInvariants[i]) {
        // RDKit✔️✔️:         result.push_back(
        // RDKit✔️✔️:             new MorganAtomEnv<OutputType>(currentInvariants[i], i, 0, &mol));
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        for atom_idx in 0..n_atoms {
            if include_atoms[atom_idx]
                && (!arguments.df_only_nonzero_invariants || current_invariants[atom_idx] != 0)
            {
                result.push(MorganAtomEnv::new(
                    u64::from(current_invariants[atom_idx]),
                    atom_idx,
                    0,
                    molecule,
                ));
            }
        }

        // RDKit✔️✔️:   for (unsigned int layer = 0; layer < morganArguments->d_radius; ++layer) {
        // RDKit✔️✔️:     std::vector<AccumTuple> allNeighborhoodsThisRound;
        for layer in 0..arguments.d_radius {
            let mut all_neighborhoods_this_round: Vec<(MorganBondEnvironment, u32, usize)> =
                Vec::new();

            // RDKit✔️✔️:     for (auto atomIdx : atomOrder) {
            // RDKit✔️✔️:       if (!deadAtoms[atomIdx]) {
            for &atom_idx in &atom_order {
                if dead_atoms[atom_idx] {
                    continue;
                }

                // RDKit✔️✔️:         const Atom *tAtom = mol.getAtomWithIdx(atomIdx);
                // RDKit✔️✔️:         if (!tAtom->getDegree()) {
                // RDKit✔️✔️:           deadAtoms.set(atomIdx, 1);
                // RDKit✔️✔️:           continue;
                // RDKit✔️✔️:         }
                let neighbors = molecule.topology_block().adjacency.neighbors_of(atom_idx);
                if neighbors.is_empty() {
                    dead_atoms[atom_idx] = true;
                    continue;
                }

                // RDKit✔️✔️:         ROMol::OEDGE_ITER beg, end;
                // RDKit✔️✔️:         boost::tie(beg, end) = mol.getAtomBonds(tAtom);
                // RDKit✔️✔️:         neighborhoodInvariants.clear();
                // RDKit✔️✔️:         while (beg != end) {
                // RDKit✔️✔️:           const Bond *bond = mol[*beg];
                // RDKit✔️❌:           roundAtomNeighborhoods[atomIdx][bond->getIdx()] = 1;
                // RDKit✔️✔️:           unsigned int oIdx = bond->getOtherAtomIdx(atomIdx);
                // RDKit✔️❌:           roundAtomNeighborhoods[atomIdx] |= atomNeighborhoods[oIdx];
                // RDKit✔️✔️:           auto bt = static_cast<int32_t>((*bondInvariants)[bond->getIdx()]);
                // RDKit✔️✔️:           neighborhoodInvariants.push_back(
                // RDKit✔️✔️:               std::make_pair(bt, currentInvariants[oIdx]));
                // RDKit✔️✔️:           ++beg;
                // RDKit✔️✔️:         }
                neighborhood_invariants.clear();
                for neighbor in neighbors {
                    let bond_idx = neighbor.bond.index();
                    round_atom_neighborhoods[atom_idx].set(bond_idx);
                    round_atom_neighborhoods[atom_idx]
                        .union_with(&atom_neighborhoods[neighbor.atom_index]);
                    neighborhood_invariants.push((
                        bond_invariants[bond_idx] as i32,
                        current_invariants[neighbor.atom_index],
                    ));
                }

                // RDKit✔️✔️:         std::sort(neighborhoodInvariants.begin(), neighborhoodInvariants.end());
                // RDKit✔️✔️:         std::uint32_t invar = layer;
                // RDKit✔️✔️:         gboost::hash_combine(invar, currentInvariants[atomIdx]);
                // RDKit✔️✔️:         bool looksChiral = (tAtom->getChiralTag() != Atom::CHI_UNSPECIFIED);
                neighborhood_invariants.sort_unstable();
                let mut invar = layer;
                hash_combine(&mut invar, current_invariants[atom_idx]);
                let atom = &molecule.atoms()[atom_idx];
                let mut looks_chiral = atom.chiral_tag() != ChiralTag::Unspecified;

                // RDKit✔️✔️:         for (std::vector<std::pair<int32_t, uint32_t>>::const_iterator it =
                // RDKit✔️✔️:                  neighborhoodInvariants.begin();
                // RDKit✔️✔️:              it != neighborhoodInvariants.end(); ++it) {
                // RDKit✔️✔️:           gboost::hash_combine(invar, *it);
                // RDKit✔️✔️:           if (morganArguments->df_includeChirality && looksChiral &&
                // RDKit✔️✔️:               !chiralAtoms[atomIdx]) {
                // RDKit✔️✔️:             if (it->first != static_cast<int32_t>(Bond::SINGLE)) {
                // RDKit✔️✔️:               looksChiral = false;
                // RDKit✔️✔️:             } else if (it != neighborhoodInvariants.begin() &&
                // RDKit✔️✔️:                        it->second == (it - 1)->second) {
                // RDKit✔️✔️:               looksChiral = false;
                // RDKit✔️✔️:             }
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:         }
                for (idx, &(bond_inv, neighbor_inv)) in neighborhood_invariants.iter().enumerate() {
                    let mut pair_hash = 0u32;
                    hash_combine(&mut pair_hash, bond_inv as u32);
                    hash_combine(&mut pair_hash, neighbor_inv);
                    hash_combine(&mut invar, pair_hash);

                    if arguments.fingerprint_arguments.df_include_chirality
                        && looks_chiral
                        && !chiral_atoms[atom_idx]
                    {
                        if bond_inv != rdkit_bond_type_code(BondOrder::Single) as i32 {
                            looks_chiral = false;
                        } else if idx > 0 && neighbor_inv == neighborhood_invariants[idx - 1].1 {
                            looks_chiral = false;
                        }
                    }
                }

                // RDKit✔️✔️:         if (morganArguments->df_includeChirality && looksChiral) {
                // RDKit✔️✔️:           chiralAtoms[atomIdx] = 1;
                // RDKit✔️✔️:           std::string cip = "";
                // RDKit✔️✔️:           tAtom->getPropIfPresent(common_properties::_CIPCode, cip);
                // RDKit✔️✔️:           if (cip == "R") {
                // RDKit✔️✔️:             gboost::hash_combine(invar, 3);
                // RDKit✔️✔️:           } else if (cip == "S") {
                // RDKit✔️✔️:             gboost::hash_combine(invar, 2);
                // RDKit✔️✔️:           } else {
                // RDKit✔️✔️:             gboost::hash_combine(invar, 1);
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:         }
                if arguments.fingerprint_arguments.df_include_chirality && looks_chiral {
                    chiral_atoms[atom_idx] = true;
                    match atom.prop("_CIPCode").unwrap_or("") {
                        "R" => hash_combine(&mut invar, 3),
                        "S" => hash_combine(&mut invar, 2),
                        _ => hash_combine(&mut invar, 1),
                    }
                }

                // RDKit✔️✔️:         nextLayerInvariants[atomIdx] = static_cast<OutputType>(invar);
                // RDKit✔️✔️:         allNeighborhoodsThisRound.push_back(
                // RDKit✔️❌:             std::make_tuple(roundAtomNeighborhoods[atomIdx],
                // RDKit✔️✔️:                             static_cast<OutputType>(invar), atomIdx));
                next_layer_invariants[atom_idx] = invar;
                all_neighborhoods_this_round.push((
                    round_atom_neighborhoods[atom_idx].clone(),
                    invar,
                    atom_idx,
                ));
            }

            // RDKit✔️❌:     std::sort(allNeighborhoodsThisRound.begin(),
            // RDKit✔️❌:               allNeighborhoodsThisRound.end());
            // RDKit✔️✔️:     for (std::vector<AccumTuple>::const_iterator iter =
            // RDKit✔️✔️:              allNeighborhoodsThisRound.begin();
            // RDKit✔️✔️:          iter != allNeighborhoodsThisRound.end(); ++iter) {
            all_neighborhoods_this_round.sort_unstable();
            for (environment, code, atom_idx) in all_neighborhoods_this_round {
                // RDKit✔️✔️:       if (morganArguments->df_includeRedundantEnvironments ||
                // RDKit✔️❌:           neighborhoods.count(std::get<0>(*iter)) == 0) {
                // RDKit✔️✔️:         if (!morganArguments->df_onlyNonzeroInvariants ||
                // RDKit✔️✔️:             (*atomInvariants)[std::get<2>(*iter)]) {
                // RDKit✔️✔️:           if (includeAtoms[std::get<2>(*iter)]) {
                // RDKit✔️✔️:             result.push_back(new MorganAtomEnv<OutputType>(
                // RDKit✔️✔️:                 std::get<1>(*iter), std::get<2>(*iter), layer + 1, &mol));
                // RDKit✔️❌:             neighborhoods.insert(std::get<0>(*iter));
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:         }
                // RDKit✔️✔️:       } else {
                // RDKit✔️✔️:         deadAtoms[std::get<2>(*iter)] = 1;
                // RDKit✔️✔️:       }
                if arguments.df_include_redundant_environments
                    || !neighborhoods.contains(&environment)
                {
                    if (!arguments.df_only_nonzero_invariants || atom_invariants[atom_idx] != 0)
                        && include_atoms[atom_idx]
                    {
                        result.push(MorganAtomEnv::new(
                            u64::from(code),
                            atom_idx,
                            layer + 1,
                            molecule,
                        ));
                        neighborhoods.insert(environment);
                    }
                } else {
                    dead_atoms[atom_idx] = true;
                }
            }

            // RDKit✔️✔️:     currentInvariants.swap(nextLayerInvariants);
            // RDKit✔️✔️:     std::fill(nextLayerInvariants.begin(), nextLayerInvariants.end(), 0);
            // RDKit✔️❌:     atomNeighborhoods = roundAtomNeighborhoods;
            // RDKit✔️✔️:   }
            std::mem::swap(&mut current_invariants, &mut next_layer_invariants);
            next_layer_invariants.fill(0);
            atom_neighborhoods.clone_from(&round_atom_neighborhoods);
        }

        // RDKit✔️✔️:   return result;
        // RDKit✔️✔️: }
        Ok(result)
    }
}

#[derive(Debug, Clone)]
pub struct MorganFingerprintGenerator {
    pub atom_environment_generator: MorganEnvGenerator,
    pub fingerprint_arguments: MorganArguments,
    pub atom_invariants_generator: MorganAtomInvariantsGenerator,
    pub bond_invariants_generator: MorganBondInvariantsGenerator,
    pub owns_atom_invariants_generator: bool,
    pub owns_bond_invariants_generator: bool,
}

impl generator::FingerprintEnvironment for MorganAtomEnv {
    fn bit_id(
        &self,
        _arguments: &FingerprintArguments,
        _atom_invariants: &[u32],
        _bond_invariants: &[u32],
        _hash_results: bool,
        _fp_size: u64,
    ) -> Result<u64, FingerprintError> {
        Ok(self.getBitId())
    }

    fn update_additional_output(&self, output: &mut AdditionalOutput, bit_id: u64) {
        self.updateAdditionalOutput(output, bit_id);
    }
}

impl generator::FingerprintFamily for MorganFingerprintGenerator {
    type Environment = MorganAtomEnv;

    fn common_arguments(&self) -> &FingerprintArguments {
        &self.fingerprint_arguments.fingerprint_arguments
    }

    fn result_size(&self) -> Result<u64, FingerprintError> {
        Ok(self.atom_environment_generator.getResultSize())
    }

    fn arguments_info_string(&self) -> String {
        self.fingerprint_arguments.infoString()
    }

    fn environment_info_string(&self) -> String {
        self.atom_environment_generator.infoString()
    }

    fn atom_invariants_info_string(&self) -> Option<String> {
        Some(match &self.atom_invariants_generator {
            MorganAtomInvariantsGenerator::Connectivity {
                include_ring_membership,
            } => MorganAtomInvGenerator::new(*include_ring_membership).infoString(),
            MorganAtomInvariantsGenerator::Feature => {
                MorganFeatureAtomInvGenerator::new().infoString()
            }
        })
    }

    fn bond_invariants_info_string(&self) -> Option<String> {
        Some(self.bond_invariants_generator.infoString())
    }

    fn arguments_json(&self) -> String {
        self.fingerprint_arguments.toJSON()
    }

    fn environment_json(&self) -> String {
        self.atom_environment_generator.toJSON()
    }

    fn atom_invariants_json(&self) -> Option<String> {
        Some(match &self.atom_invariants_generator {
            MorganAtomInvariantsGenerator::Connectivity {
                include_ring_membership,
            } => MorganAtomInvGenerator::new(*include_ring_membership).toJSON(),
            MorganAtomInvariantsGenerator::Feature => MorganFeatureAtomInvGenerator::new().toJSON(),
        })
    }

    fn bond_invariants_json(&self) -> Option<String> {
        Some(self.bond_invariants_generator.toJSON())
    }

    fn atom_invariants(&self, molecule: &Molecule) -> Result<Vec<u32>, FingerprintError> {
        self.atom_invariants_generator.get_atom_invariants(molecule)
    }

    fn bond_invariants(&self, molecule: &Molecule) -> Result<Vec<u32>, FingerprintError> {
        self.bond_invariants_generator
            .try_get_bond_invariants(molecule)
    }

    fn environments(
        &self,
        molecule: &Molecule,
        from_atoms: Option<&[usize]>,
        ignore_atoms: Option<&[usize]>,
        _conf_id: i32,
        atom_invariants: &[u32],
        bond_invariants: &[u32],
        _hash_results: bool,
    ) -> Result<Vec<Self::Environment>, FingerprintError> {
        self.atom_environment_generator.getEnvironments(
            molecule,
            &self.fingerprint_arguments,
            from_atoms,
            ignore_atoms,
            atom_invariants,
            bond_invariants,
        )
    }
}

/// RDKit's Topological Torsion generator is instantiated only with a 64-bit
/// output type; Rust fixes the bit-id and sparse-vector domain to `u64`.
// RDKit✔️✔️: // Topological torsion fingerprint does not support 32 bit output yet
#[derive(Debug, Clone)]
pub struct TopologicalTorsionFingerprintGenerator {
    pub atom_environment_generator: TopologicalTorsionEnvGenerator,
    pub fingerprint_arguments: TopologicalTorsionArguments,
    pub atom_invariants_generator: AtomPairAtomInvGenerator,
    pub owns_atom_invariants_generator: bool,
}

impl generator::FingerprintEnvironment for TopologicalTorsionAtomEnv {
    fn bit_id(
        &self,
        _arguments: &FingerprintArguments,
        _atom_invariants: &[u32],
        _bond_invariants: &[u32],
        _hash_results: bool,
        _fp_size: u64,
    ) -> Result<u64, FingerprintError> {
        Ok(self.getBitId())
    }

    fn update_additional_output(&self, output: &mut AdditionalOutput, bit_id: u64) {
        self.updateAdditionalOutput(output, bit_id);
    }
}

impl generator::FingerprintFamily for TopologicalTorsionFingerprintGenerator {
    type Environment = TopologicalTorsionAtomEnv;

    fn common_arguments(&self) -> &FingerprintArguments {
        &self.fingerprint_arguments.fingerprint_arguments
    }

    fn result_size(&self) -> Result<u64, FingerprintError> {
        self.atom_environment_generator
            .get_result_size(&self.fingerprint_arguments)
    }

    fn arguments_info_string(&self) -> String {
        self.fingerprint_arguments.infoString()
    }

    fn environment_info_string(&self) -> String {
        self.atom_environment_generator.infoString()
    }

    fn atom_invariants_info_string(&self) -> Option<String> {
        Some(self.atom_invariants_generator.infoString())
    }

    fn bond_invariants_info_string(&self) -> Option<String> {
        None
    }

    fn arguments_json(&self) -> String {
        self.fingerprint_arguments.toJSON()
    }

    fn environment_json(&self) -> String {
        self.atom_environment_generator.toJSON()
    }

    fn atom_invariants_json(&self) -> Option<String> {
        Some(self.atom_invariants_generator.toJSON())
    }

    fn bond_invariants_json(&self) -> Option<String> {
        None
    }

    fn atom_invariants(&self, molecule: &Molecule) -> Result<Vec<u32>, FingerprintError> {
        self.atom_invariants_generator.getAtomInvariants(molecule)
    }

    fn bond_invariants(&self, _molecule: &Molecule) -> Result<Vec<u32>, FingerprintError> {
        Ok(Vec::new())
    }

    fn environments(
        &self,
        molecule: &Molecule,
        from_atoms: Option<&[usize]>,
        ignore_atoms: Option<&[usize]>,
        _conf_id: i32,
        atom_invariants: &[u32],
        _bond_invariants: &[u32],
        hash_results: bool,
    ) -> Result<Vec<Self::Environment>, FingerprintError> {
        self.atom_environment_generator.getEnvironments(
            molecule,
            &self.fingerprint_arguments,
            from_atoms,
            ignore_atoms,
            atom_invariants,
            hash_results,
        )
    }
}

macro_rules! fingerprint_generator_api {
    ($generator:ty) => {
        impl $generator {
            #[must_use]
            pub fn info_string(&self) -> String {
                generator::FingerprintGenerator::new(self).info_string()
            }

            #[must_use]
            pub fn to_json(&self) -> String {
                generator::FingerprintGenerator::new(self).to_json()
            }

            pub fn sparse_count_fingerprint(
                &self,
                molecule: &Molecule,
                args: &mut FingerprintFuncArguments,
            ) -> Result<SparseCountFingerprint, FingerprintError> {
                generator::FingerprintGenerator::new(self)
                    .get_sparse_count_fingerprint(molecule, args)
            }

            pub fn sparse_fingerprint(
                &self,
                molecule: &Molecule,
                args: &mut FingerprintFuncArguments,
            ) -> Result<SparseBitFingerprint, FingerprintError> {
                generator::FingerprintGenerator::new(self).get_sparse_fingerprint(molecule, args)
            }

            pub fn count_fingerprint(
                &self,
                molecule: &Molecule,
                args: &mut FingerprintFuncArguments,
            ) -> Result<SparseCountFingerprint, FingerprintError> {
                generator::FingerprintGenerator::new(self).get_count_fingerprint(molecule, args)
            }

            pub fn fingerprint(
                &self,
                molecule: &Molecule,
                args: &mut FingerprintFuncArguments,
            ) -> Result<Fingerprint, FingerprintError> {
                generator::FingerprintGenerator::new(self).get_fingerprint(molecule, args)
            }

            #[allow(non_snake_case)]
            pub fn getSparseCountFingerprint(
                &self,
                molecule: &Molecule,
                args: &mut FingerprintFuncArguments,
            ) -> Result<SparseCountFingerprint, FingerprintError> {
                self.sparse_count_fingerprint(molecule, args)
            }

            #[allow(non_snake_case)]
            pub fn getSparseFingerprint(
                &self,
                molecule: &Molecule,
                args: &mut FingerprintFuncArguments,
            ) -> Result<SparseBitFingerprint, FingerprintError> {
                self.sparse_fingerprint(molecule, args)
            }

            #[allow(non_snake_case)]
            pub fn getCountFingerprint(
                &self,
                molecule: &Molecule,
                args: &mut FingerprintFuncArguments,
            ) -> Result<SparseCountFingerprint, FingerprintError> {
                self.count_fingerprint(molecule, args)
            }

            #[allow(non_snake_case)]
            pub fn getFingerprint(
                &self,
                molecule: &Molecule,
                args: &mut FingerprintFuncArguments,
            ) -> Result<Fingerprint, FingerprintError> {
                self.fingerprint(molecule, args)
            }

            #[allow(non_snake_case)]
            pub fn getFingerprintHelper(
                &self,
                molecule: &Molecule,
                args: &mut FingerprintFuncArguments,
                fp_size: u64,
            ) -> Result<SparseCountFingerprint, FingerprintError> {
                generator::FingerprintGenerator::new(self)
                    .get_fingerprint_helper(molecule, args, fp_size)
            }

            #[allow(non_snake_case)]
            pub fn getFingerprints(
                &self,
                molecules: &[Option<&Molecule>],
                num_threads: i32,
            ) -> Result<Vec<Option<Fingerprint>>, FingerprintError> {
                generator::FingerprintGenerator::new(self).get_fingerprints(molecules, num_threads)
            }

            #[allow(non_snake_case)]
            pub fn getSparseFingerprints(
                &self,
                molecules: &[Option<&Molecule>],
                num_threads: i32,
            ) -> Result<Vec<Option<SparseBitFingerprint>>, FingerprintError> {
                generator::FingerprintGenerator::new(self)
                    .get_sparse_fingerprints(molecules, num_threads)
            }

            #[allow(non_snake_case)]
            pub fn getCountFingerprints(
                &self,
                molecules: &[Option<&Molecule>],
                num_threads: i32,
            ) -> Result<Vec<Option<SparseCountFingerprint>>, FingerprintError> {
                generator::FingerprintGenerator::new(self)
                    .get_count_fingerprints(molecules, num_threads)
            }

            #[allow(non_snake_case)]
            pub fn getSparseCountFingerprints(
                &self,
                molecules: &[Option<&Molecule>],
                num_threads: i32,
            ) -> Result<Vec<Option<SparseCountFingerprint>>, FingerprintError> {
                generator::FingerprintGenerator::new(self)
                    .get_sparse_count_fingerprints(molecules, num_threads)
            }
        }
    };
}

fingerprint_generator_api!(MorganFingerprintGenerator);
fingerprint_generator_api!(TopologicalTorsionFingerprintGenerator);

/// Typed source generator restored by the shared JSON and public dispatch path.
pub type TypedFingerprintGenerator = generator::RestoredFingerprintGenerator;

/// RDKit's default fingerprint-family selector.
// RDKit source: FingerprintGenerator.h lines 485-490
// RDKit✔️✔️:   AtomPairFP,
// RDKit✔️✔️:   MorganFP,
// RDKit❌❌:   RDKitFP,
// RDKit✔️✔️:   TopologicalTorsionFP
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FPType {
    AtomPairFP,
    MorganFP,
    RDKitFP,
    TopologicalTorsionFP,
}

fn default_generator_for_fp_type(
    fingerprint_type: FPType,
) -> Result<TypedFingerprintGenerator, FingerprintError> {
    // RDKit source: FingerprintGenerator.cpp lines 850-874; the same switch is
    // repeated for all four source bulk output functions.
    // RDKit✔️✔️:   switch (fPType) {
    // RDKit✔️✔️:     case FPType::AtomPairFP: {
    // RDKit✔️✔️:       generator = AtomPair::getAtomPairGenerator<std::uint64_t>();
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     case FPType::MorganFP: {
    // RDKit✔️✔️:       generator = MorganFingerprint::getMorganGenerator<std::uint64_t>(2);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit❌❌:     case FPType::RDKitFP: {
    // RDKit❌❌:       generator = RDKitFP::getRDKitFPGenerator<std::uint64_t>();
    // RDKit❌❌:       break;
    // RDKit❌❌:     }
    // RDKit✔️✔️:     case FPType::TopologicalTorsionFP: {
    // RDKit✔️✔️:       generator =
    // RDKit✔️✔️:           TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     default: {
    // RDKit✔️✔️:       throw UnimplementedFPException(
    // RDKit✔️✔️:           "Fingerprint type not implemented for getSparseCountFP");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    match fingerprint_type {
        FPType::MorganFP => Ok(TypedFingerprintGenerator::Morgan(getMorganGenerator(
            &MorganArguments::default(),
            None,
            None,
            true,
            true,
        ))),
        FPType::TopologicalTorsionFP => Ok(TypedFingerprintGenerator::TopologicalTorsion(
            getTopologicalTorsionGenerator(&TopologicalTorsionArguments::default(), None, true)?,
        )),
        FPType::AtomPairFP => Ok(TypedFingerprintGenerator::AtomPair(
            atom_pair::atom_pair_generator(&atom_pair::AtomPairArguments::default(), None, true),
        )),
        FPType::RDKitFP => Err(FingerprintError::UnsupportedOption {
            option: "FPType::RDKitFP",
            reason: "the modern RDKitFP generator is outside the modeled shared generator core",
        }),
    }
}

fn scalar_dispatch_result<T>(mut results: Vec<Option<T>>) -> Result<T, FingerprintError> {
    results
        .pop()
        .flatten()
        .ok_or(FingerprintError::InvalidArguments {
            reason: "scalar fingerprint dispatch produced no result",
        })
}

#[allow(non_snake_case)]
pub fn getSparseCountFP(
    molecule: &Molecule,
    fingerprint_type: FPType,
) -> Result<SparseCountFingerprint, FingerprintError> {
    // RDKit source: FingerprintGenerator.cpp lines 826-830
    // RDKit✔️✔️: SparseIntVect<std::uint64_t> *getSparseCountFP(const ROMol &mol,
    // RDKit✔️✔️:                                                FPType fPType) {
    // RDKit✔️✔️:   std::vector<const ROMol *> tempVect(1, &mol);
    // RDKit✔️✔️:   return (*getSparseCountFPBulk(tempVect, fPType))[0];
    // RDKit✔️✔️: }
    scalar_dispatch_result(getSparseCountFPBulk(&[Some(molecule)], fingerprint_type)?)
}

#[allow(non_snake_case)]
pub fn getSparseFP(
    molecule: &Molecule,
    fingerprint_type: FPType,
) -> Result<SparseBitFingerprint, FingerprintError> {
    // RDKit source: FingerprintGenerator.cpp lines 832-835
    // RDKit✔️✔️: SparseBitVect *getSparseFP(const ROMol &mol, FPType fPType) {
    // RDKit✔️✔️:   std::vector<const ROMol *> tempVect(1, &mol);
    // RDKit✔️✔️:   return (*getSparseFPBulk(tempVect, fPType))[0];
    // RDKit✔️✔️: }
    scalar_dispatch_result(getSparseFPBulk(&[Some(molecule)], fingerprint_type)?)
}

#[allow(non_snake_case)]
pub fn getCountFP(
    molecule: &Molecule,
    fingerprint_type: FPType,
) -> Result<SparseCountFingerprint, FingerprintError> {
    // RDKit source: FingerprintGenerator.cpp lines 837-840
    // RDKit✔️✔️: SparseIntVect<std::uint32_t> *getCountFP(const ROMol &mol, FPType fPType) {
    // RDKit✔️✔️:   std::vector<const ROMol *> tempVect(1, &mol);
    // RDKit✔️✔️:   return (*getCountFPBulk(tempVect, fPType))[0];
    // RDKit✔️✔️: }
    scalar_dispatch_result(getCountFPBulk(&[Some(molecule)], fingerprint_type)?)
}

#[allow(non_snake_case)]
pub fn getFP(
    molecule: &Molecule,
    fingerprint_type: FPType,
) -> Result<Fingerprint, FingerprintError> {
    // RDKit source: FingerprintGenerator.cpp lines 842-845
    // RDKit✔️✔️: ExplicitBitVect *getFP(const ROMol &mol, FPType fPType) {
    // RDKit✔️✔️:   std::vector<const ROMol *> tempVect(1, &mol);
    // RDKit✔️✔️:   return (*getFPBulk(tempVect, fPType))[0];
    // RDKit✔️✔️: }
    scalar_dispatch_result(getFPBulk(&[Some(molecule)], fingerprint_type)?)
}

#[allow(non_snake_case)]
pub fn getSparseCountFPBulk(
    molecules: &[Option<&Molecule>],
    fingerprint_type: FPType,
) -> Result<Vec<Option<SparseCountFingerprint>>, FingerprintError> {
    // RDKit source: FingerprintGenerator.cpp lines 847-881
    // RDKit✔️✔️: std::vector<SparseIntVect<std::uint64_t> *> *getSparseCountFPBulk(
    // RDKit✔️✔️:     const std::vector<const ROMol *> molVector, FPType fPType) {
    // RDKit✔️✔️:   FingerprintGenerator<std::uint64_t> *generator = nullptr;
    // RDKit✔️✔️:   auto *res = new std::vector<SparseIntVect<std::uint64_t> *>();
    // RDKit✔️✔️:   for (const auto *mol : molVector) {
    // RDKit✔️✔️:     res->push_back(generator->getSparseCountFingerprint(*mol));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   delete generator;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    default_generator_for_fp_type(fingerprint_type)?.getSparseCountFingerprints(molecules, 1)
}

#[allow(non_snake_case)]
pub fn getSparseFPBulk(
    molecules: &[Option<&Molecule>],
    fingerprint_type: FPType,
) -> Result<Vec<Option<SparseBitFingerprint>>, FingerprintError> {
    // RDKit source: FingerprintGenerator.cpp lines 883-917
    // RDKit✔️✔️: std::vector<SparseBitVect *> *getSparseFPBulk(
    // RDKit✔️✔️:     const std::vector<const ROMol *> molVector, FPType fPType) {
    // RDKit✔️✔️:   FingerprintGenerator<std::uint64_t> *generator = nullptr;
    // RDKit✔️✔️:   auto *res = new std::vector<SparseBitVect *>();
    // RDKit✔️✔️:   for (const auto *mol : molVector) {
    // RDKit✔️✔️:     res->push_back(generator->getSparseFingerprint(*mol));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   delete generator;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    default_generator_for_fp_type(fingerprint_type)?.getSparseFingerprints(molecules, 1)
}

#[allow(non_snake_case)]
pub fn getCountFPBulk(
    molecules: &[Option<&Molecule>],
    fingerprint_type: FPType,
) -> Result<Vec<Option<SparseCountFingerprint>>, FingerprintError> {
    // RDKit source: FingerprintGenerator.cpp lines 919-953
    // RDKit✔️✔️: std::vector<SparseIntVect<std::uint32_t> *> *getCountFPBulk(
    // RDKit✔️✔️:     const std::vector<const ROMol *> molVector, FPType fPType) {
    // RDKit✔️✔️:   FingerprintGenerator<std::uint64_t> *generator = nullptr;
    // RDKit✔️✔️:   auto *res = new std::vector<SparseIntVect<std::uint32_t> *>();
    // RDKit✔️✔️:   for (const auto *mol : molVector) {
    // RDKit✔️✔️:     res->push_back(generator->getCountFingerprint(*mol));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   delete generator;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    default_generator_for_fp_type(fingerprint_type)?.getCountFingerprints(molecules, 1)
}

#[allow(non_snake_case)]
pub fn getFPBulk(
    molecules: &[Option<&Molecule>],
    fingerprint_type: FPType,
) -> Result<Vec<Option<Fingerprint>>, FingerprintError> {
    // RDKit source: FingerprintGenerator.cpp lines 955-989
    // RDKit✔️✔️: std::vector<ExplicitBitVect *> *getFPBulk(
    // RDKit✔️✔️:     const std::vector<const ROMol *> molVector, FPType fPType) {
    // RDKit✔️✔️:   FingerprintGenerator<std::uint64_t> *generator = nullptr;
    // RDKit✔️✔️:   auto *res = new std::vector<ExplicitBitVect *>();
    // RDKit✔️✔️:   for (const auto *mol : molVector) {
    // RDKit✔️✔️:     res->push_back(generator->getFingerprint(*mol));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   delete generator;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    default_generator_for_fp_type(fingerprint_type)?.getFingerprints(molecules, 1)
}

#[allow(non_snake_case)]
pub fn getMorganGenerator(
    args: &MorganArguments,
    atom_invariants_generator: Option<MorganAtomInvariantsGenerator>,
    bond_invariants_generator: Option<MorganBondInvariantsGenerator>,
    owns_atom_inv_gen: bool,
    owns_bond_inv_gen: bool,
) -> MorganFingerprintGenerator {
    // RDKit source: MorganGenerator.cpp lines 479-502
    // RDKit✔️✔️: template <typename OutputType>
    // RDKit✔️✔️: FingerprintGenerator<OutputType> *getMorganGenerator(
    // RDKit✔️✔️:     const MorganArguments &args,
    // RDKit✔️✔️:     AtomInvariantsGenerator *atomInvariantsGenerator,
    // RDKit✔️✔️:     BondInvariantsGenerator *bondInvariantsGenerator, bool ownsAtomInvGen,
    // RDKit✔️✔️:     bool ownsBondInvGen) {
    // RDKit✔️✔️:   AtomEnvironmentGenerator<OutputType> *morganEnvGenerator =
    // RDKit✔️✔️:       new MorganEnvGenerator<OutputType>();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bool ownsAtomInvGenerator = ownsAtomInvGen;
    let atom_environment_generator = MorganEnvGenerator::new();
    let mut owns_atom_invariants_generator = owns_atom_inv_gen;

    // RDKit✔️✔️:   if (!atomInvariantsGenerator) {
    // RDKit✔️✔️:     atomInvariantsGenerator = new MorganAtomInvGenerator();
    // RDKit✔️✔️:     ownsAtomInvGenerator = true;
    // RDKit✔️✔️:   }
    let atom_invariants_generator = match atom_invariants_generator {
        Some(generator) => generator,
        None => {
            owns_atom_invariants_generator = true;
            MorganAtomInvariantsGenerator::Connectivity {
                include_ring_membership: true,
            }
        }
    };

    // RDKit✔️✔️:   bool ownsBondInvGenerator = ownsBondInvGen;
    let mut owns_bond_invariants_generator = owns_bond_inv_gen;

    // RDKit✔️✔️:   if (!bondInvariantsGenerator) {
    // RDKit✔️✔️:     bondInvariantsGenerator = new MorganBondInvGenerator(
    // RDKit✔️✔️:         args.df_useBondTypes, args.df_includeChirality);
    // RDKit✔️✔️:     ownsBondInvGenerator = true;
    // RDKit✔️✔️:   }
    let bond_invariants_generator = match bond_invariants_generator {
        Some(generator) => generator,
        None => {
            owns_bond_invariants_generator = true;
            MorganBondInvariantsGenerator::new(
                args.df_use_bond_types,
                args.fingerprint_arguments.df_include_chirality,
            )
        }
    };

    // RDKit✔️✔️:   return new FingerprintGenerator<OutputType>(
    // RDKit✔️✔️:       morganEnvGenerator, new MorganArguments(args), atomInvariantsGenerator,
    // RDKit✔️✔️:       bondInvariantsGenerator, ownsAtomInvGenerator, ownsBondInvGenerator);
    // RDKit✔️✔️: }
    MorganFingerprintGenerator {
        atom_environment_generator,
        fingerprint_arguments: args.clone(),
        atom_invariants_generator,
        bond_invariants_generator,
        owns_atom_invariants_generator,
        owns_bond_invariants_generator,
    }
}

#[allow(non_snake_case)]
pub fn getTopologicalTorsionGenerator(
    args: &TopologicalTorsionArguments,
    atom_invariants_generator: Option<AtomPairAtomInvGenerator>,
    owns_atom_inv_gen: bool,
) -> Result<TopologicalTorsionFingerprintGenerator, FingerprintError> {
    // RDKit source: TopologicalTorsionGenerator.cpp lines 212-229
    // RDKit✔️✔️: template <typename OutputType>
    // RDKit✔️✔️: FingerprintGenerator<OutputType> *getTopologicalTorsionGenerator(
    // RDKit✔️✔️:     const TopologicalTorsionArguments &args,
    // RDKit✔️✔️:     AtomInvariantsGenerator *atomInvariantsGenerator,
    // RDKit✔️✔️:     const bool ownsAtomInvGen) {
    // RDKit✔️✔️:   auto *envGenerator = new TopologicalTorsionEnvGenerator<OutputType>();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bool ownsAtomInvGenerator = ownsAtomInvGen;
    // RDKit✔️✔️:   if (!atomInvariantsGenerator) {
    // RDKit✔️✔️:     atomInvariantsGenerator =
    // RDKit✔️✔️:         new AtomPair::AtomPairAtomInvGenerator(args.df_includeChirality, true);
    // RDKit✔️✔️:     ownsAtomInvGenerator = true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return new FingerprintGenerator<OutputType>(
    // RDKit✔️✔️:       envGenerator, new TopologicalTorsionArguments(args),
    // RDKit✔️✔️:       atomInvariantsGenerator, nullptr, ownsAtomInvGenerator, false);
    // RDKit✔️✔️: };
    let atom_environment_generator = TopologicalTorsionEnvGenerator::new();
    atom_environment_generator.get_result_size(args)?;
    let mut owns_atom_invariants_generator = owns_atom_inv_gen;
    let atom_invariants_generator = match atom_invariants_generator {
        Some(generator) => generator,
        None => {
            owns_atom_invariants_generator = true;
            AtomPairAtomInvGenerator::new(args.fingerprint_arguments.df_include_chirality, true)
        }
    };

    Ok(TopologicalTorsionFingerprintGenerator {
        atom_environment_generator,
        fingerprint_arguments: args.clone(),
        atom_invariants_generator,
        owns_atom_invariants_generator,
    })
}

#[allow(clippy::too_many_arguments, non_snake_case)]
pub fn getTopologicalTorsionGeneratorWithParams(
    include_chirality: bool,
    torsion_atom_count: u32,
    atom_invariants_generator: Option<AtomPairAtomInvGenerator>,
    count_simulation: bool,
    fp_size: u32,
    count_bounds: Vec<u32>,
    owns_atom_inv_gen: bool,
) -> Result<TopologicalTorsionFingerprintGenerator, FingerprintError> {
    // RDKit source: TopologicalTorsionGenerator.cpp lines 230-240
    // RDKit✔️✔️: template <typename OutputType>
    // RDKit✔️✔️: FingerprintGenerator<OutputType> *getTopologicalTorsionGenerator(
    // RDKit✔️✔️:     bool includeChirality, uint32_t torsionAtomCount,
    // RDKit✔️✔️:     AtomInvariantsGenerator *atomInvariantsGenerator, bool countSimulation,
    // RDKit✔️✔️:     std::uint32_t fpSize, std::vector<std::uint32_t> countBounds,
    // RDKit✔️✔️:     bool ownsAtomInvGen) {
    // RDKit✔️✔️:   TopologicalTorsionArguments arguments(includeChirality, torsionAtomCount,
    // RDKit✔️✔️:                                         countSimulation, countBounds, fpSize);
    // RDKit✔️✔️:   return getTopologicalTorsionGenerator<OutputType>(
    // RDKit✔️✔️:       arguments, atomInvariantsGenerator, ownsAtomInvGen);
    // RDKit✔️✔️: };
    let arguments = TopologicalTorsionArguments::new(
        include_chirality,
        torsion_atom_count,
        count_simulation,
        count_bounds,
        fp_size,
    )?;
    getTopologicalTorsionGenerator(&arguments, atom_invariants_generator, owns_atom_inv_gen)
}

/// Legacy un-hashed Topological Torsion fingerprint compatibility entry.
#[deprecated(note = "please use TopologicalTorsionFingerprintGenerator")]
#[allow(clippy::too_many_arguments, non_snake_case)]
pub fn getTopologicalTorsionFingerprint(
    molecule: &Molecule,
    target_size: u32,
    from_atoms: Option<&[usize]>,
    ignore_atoms: Option<&[usize]>,
    atom_invariants: Option<&[u32]>,
    include_chirality: bool,
) -> Result<SparseCountFingerprint, FingerprintError> {
    // RDKit source: AtomPairs.cpp lines 159-265
    // A Rust deprecation attribute provides the source warning at API use
    // time; it is compile-time instead of RDKit's runtime log emission.
    // RDKit❗✔️:   RDLog::deprecationWarning("please use TopologicalTorsionGenerator");
    // RDKit✔️✔️:   PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit✔️✔️:                "bad atomInvariants size");
    // RDKit✔️✔️:   const ROMol *lmol = &mol;
    // RDKit✔️✔️:   std::unique_ptr<ROMol> tmol;
    // RDKit✔️✔️:   if (includeChirality && !mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit✔️✔️:     tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    // RDKit✔️✔️:     MolOps::assignStereochemistry(*tmol);
    // RDKit✔️✔️:     lmol = tmol.get();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   boost::uint64_t sz = 1;
    // RDKit✔️✔️:   sz = (sz << (targetSize *
    // RDKit✔️✔️:                (codeSize + (includeChirality ? numChiralBits : 0))));
    // RDKit✔️✔️:   // NOTE: this -1 is incorrect but it's needed for backwards compatibility.
    // RDKit✔️✔️:   //  hopefully we'll never have a case with a torsion that hits this.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //  mmm, bug compatible.
    // RDKit✔️✔️:   sz -= 1;
    // RDKit✔️✔️:   auto *res = new SparseIntVect<boost::int64_t>(sz);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<std::uint32_t> atomCodes;
    // RDKit✔️✔️:   atomCodes.reserve(lmol->getNumAtoms());
    // RDKit✔️✔️:   for (ROMol::ConstAtomIterator atomItI = lmol->beginAtoms();
    // RDKit✔️✔️:        atomItI != lmol->endAtoms(); ++atomItI) {
    // RDKit✔️✔️:     if (!atomInvariants) {
    // RDKit✔️✔️:       atomCodes.push_back(getAtomCode(*atomItI, 0, includeChirality));
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       // need to add to the atomCode here because we subtract off up to 2 below
    // RDKit✔️✔️:       // as part of the branch correction
    // RDKit✔️✔️:       atomCodes.push_back(
    // RDKit✔️✔️:           (*atomInvariants)[(*atomItI)->getIdx()] % ((1 << codeSize) - 1) + 2);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   boost::dynamic_bitset<> *fromAtomsBV = nullptr;
    // RDKit✔️✔️:   if (fromAtoms) {
    // RDKit✔️✔️:     fromAtomsBV = new boost::dynamic_bitset<>(lmol->getNumAtoms());
    // RDKit✔️✔️:     for (auto fAt : *fromAtoms) {
    // RDKit✔️✔️:       fromAtomsBV->set(fAt);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   boost::dynamic_bitset<> *ignoreAtomsBV = nullptr;
    // RDKit✔️✔️:   if (ignoreAtoms) {
    // RDKit✔️✔️:     ignoreAtomsBV = new boost::dynamic_bitset<>(mol.getNumAtoms());
    // RDKit✔️✔️:     for (auto fAt : *ignoreAtoms) {
    // RDKit✔️✔️:       ignoreAtomsBV->set(fAt);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   boost::dynamic_bitset<> pAtoms(lmol->getNumAtoms());
    // RDKit✔️✔️:   PATH_LIST paths = findAllPathsOfLengthN(*lmol, targetSize, false);
    // RDKit✔️✔️:   for (PATH_LIST::const_iterator pathIt = paths.begin(); pathIt != paths.end();
    // RDKit✔️✔️:        ++pathIt) {
    // RDKit✔️✔️:     bool keepIt = true;
    // RDKit✔️✔️:     if (fromAtomsBV) {
    // RDKit✔️✔️:       keepIt = false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::vector<std::uint32_t> pathCodes;
    // RDKit✔️✔️:     const PATH_TYPE &path = *pathIt;
    // RDKit✔️✔️:     if (fromAtomsBV) {
    // RDKit✔️✔️:       if (fromAtomsBV->test(static_cast<std::uint32_t>(path.front())) ||
    // RDKit✔️✔️:           fromAtomsBV->test(static_cast<std::uint32_t>(path.back()))) {
    // RDKit✔️✔️:         keepIt = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (keepIt && ignoreAtomsBV) {
    // RDKit✔️✔️:       for (auto pElem : path) {
    // RDKit✔️✔️:         if (ignoreAtomsBV->test(pElem)) {
    // RDKit✔️✔️:           keepIt = false;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (keepIt) {
    // RDKit✔️✔️:       pAtoms.reset();
    // RDKit✔️✔️:       for (auto pIt = path.begin(); pIt < path.end(); ++pIt) {
    // RDKit✔️✔️:         // look for a cycle that doesn't start at the first atom
    // RDKit✔️✔️:         // we can't effectively canonicalize these at the moment
    // RDKit✔️✔️:         // (was github #811)
    // RDKit✔️✔️:         if (pIt != path.begin() && *pIt != *(path.begin()) && pAtoms[*pIt]) {
    // RDKit✔️✔️:           pathCodes.clear();
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         pAtoms.set(*pIt);
    // RDKit✔️✔️:         unsigned int code = atomCodes[*pIt] - 1;
    // RDKit✔️✔️:         // subtract off the branching number:
    // RDKit✔️✔️:         if (pIt != path.begin() && pIt + 1 != path.end()) {
    // RDKit✔️✔️:           --code;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         pathCodes.push_back(code);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (pathCodes.size()) {
    // RDKit✔️✔️:         boost::int64_t code =
    // RDKit✔️✔️:             getTopologicalTorsionCode(pathCodes, includeChirality);
    // RDKit✔️✔️:         updateElement(*res, code);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   delete fromAtomsBV;
    // RDKit✔️✔️:   delete ignoreAtomsBV;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    if atom_invariants.is_some_and(|invariants| invariants.len() < molecule.num_atoms()) {
        return Err(FingerprintError::InvalidArguments {
            reason: "bad atom invariants size",
        });
    }

    let arguments = TopologicalTorsionArguments::new(
        include_chirality,
        target_size,
        false,
        vec![1, 2, 4, 8],
        2048,
    )?;
    let mut generator = getTopologicalTorsionGenerator(&arguments, None, true)?;
    if atom_invariants.is_none() {
        // AtomPairAtomInvGenerator has already applied its source `- 2`
        // correction. The legacy environment mode adds one at endpoints and
        // leaves internal atoms unchanged, reproducing legacy `atomCode - 1`
        // followed by the internal branch correction without duplicating the
        // path-enumeration core.
        generator
            .atom_environment_generator
            .use_legacy_unfolded_atom_codes();
    }
    let prepared_molecule;
    let fingerprint_molecule = if include_chirality && molecule.prop("_StereochemDone").is_none() {
        // RDKit source: AtomPairs.cpp lines 169-174
        // RDKit❗✔️: if (includeChirality && !mol.hasProp(common_properties::_StereochemDone)) {
        // RDKit❗✔️:   tmol = std::unique_ptr<ROMol>(new ROMol(mol));
        // RDKit❗✔️:   MolOps::assignStereochemistry(*tmol);
        // RDKit❗✔️:   lmol = tmol.get();
        // RDKit❗✔️: }
        prepared_molecule = {
            let mut copy = molecule.clone();
            crate::smiles::assign_stereochemistry_cleanup_subset(&mut copy, false)?;
            copy
        };
        &prepared_molecule
    } else {
        molecule
    };
    let mut call_arguments = FingerprintFuncArguments {
        from_atoms: from_atoms.map(<[usize]>::to_vec),
        ignore_atoms: ignore_atoms.map(<[usize]>::to_vec),
        custom_atom_invariants: atom_invariants.map(<[u32]>::to_vec),
        ..FingerprintFuncArguments::default()
    };
    let mut result =
        generator.getSparseCountFingerprint(fingerprint_molecule, &mut call_arguments)?;
    let compatibility_size = generator
        .atom_environment_generator
        .get_result_size(&generator.fingerprint_arguments)?
        .checked_sub(1)
        .ok_or(FingerprintError::InvalidArguments {
            reason: "legacy topological torsion compatibility size underflow",
        })?;
    if result
        .nonzero_elements()
        .keys()
        .any(|&bit_id| bit_id >= compatibility_size)
    {
        return Err(FingerprintError::InvalidArguments {
            reason: "legacy topological torsion id reaches the compatibility vector boundary",
        });
    }
    result.size = compatibility_size;
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn torsion_fp_calc(
    molecule: &Molecule,
    n_bits: u32,
    target_size: u32,
    from_atoms: Option<&[usize]>,
    ignore_atoms: Option<&[usize]>,
    atom_invariants: Option<&[u32]>,
    include_chirality: bool,
) -> Result<SparseCountFingerprint, FingerprintError> {
    // RDKit source: AtomPairs.cpp lines 267-297
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: void TorsionFpCalc(T *res, const ROMol &mol, unsigned int nBits,
    // RDKit✔️✔️:                    unsigned int targetSize,
    // RDKit✔️✔️:                    const std::vector<std::uint32_t> *fromAtoms,
    // RDKit✔️✔️:                    const std::vector<std::uint32_t> *ignoreAtoms,
    // RDKit✔️✔️:                    const std::vector<std::uint32_t> *atomInvariants,
    // RDKit✔️✔️:                    bool includeChirality) {
    // RDKit✔️✔️:   PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit✔️✔️:                "bad atomInvariants size");
    // RDKit✔️✔️:   const ROMol *lmol = &mol;
    // RDKit✔️✔️:   std::unique_ptr<ROMol> tmol;
    // RDKit✔️✔️:   if (includeChirality && !mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit✔️✔️:     tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    // RDKit✔️✔️:     MolOps::assignStereochemistry(*tmol);
    // RDKit✔️✔️:     lmol = tmol.get();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpgen{
    // RDKit✔️✔️:       RDKit::TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>(
    // RDKit✔️✔️:           includeChirality, targetSize, nullptr, true, nBits)};
    // RDKit✔️✔️:   FingerprintFuncArguments args;
    // RDKit✔️✔️:   args.fromAtoms = fromAtoms;
    // RDKit✔️✔️:   args.ignoreAtoms = ignoreAtoms;
    // RDKit✔️✔️:   args.customAtomInvariants = atomInvariants;
    // RDKit✔️✔️:   auto siv = fpgen->getCountFingerprint(*lmol, args);
    // RDKit✔️🔝:   for (auto v : siv->getNonzeroElements()) {
    // RDKit✔️🔝:     res->setVal(v.first, v.second);
    // RDKit✔️🔝:   }
    // RDKit✔️✔️: }
    // Returning the shared generator's owned count vector preserves every
    // element and its declared size while eliminating the source adapter's
    // second sparse-container allocation and O(k) element copy.
    if n_bits == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "legacy hashed topological torsion nBits must be greater than zero",
        });
    }
    if atom_invariants.is_some_and(|invariants| invariants.len() < molecule.num_atoms()) {
        return Err(FingerprintError::InvalidArguments {
            reason: "bad atom invariants size",
        });
    }

    let generator = getTopologicalTorsionGeneratorWithParams(
        include_chirality,
        target_size,
        None,
        true,
        n_bits,
        vec![1, 2, 4, 8],
        true,
    )?;
    let mut args = FingerprintFuncArguments {
        from_atoms: from_atoms.map(<[usize]>::to_vec),
        ignore_atoms: ignore_atoms.map(<[usize]>::to_vec),
        custom_atom_invariants: atom_invariants.map(<[u32]>::to_vec),
        ..FingerprintFuncArguments::default()
    };
    generator.getCountFingerprint(molecule, &mut args)
}

#[deprecated(note = "please use TopologicalTorsionGenerator")]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub fn getHashedTopologicalTorsionFingerprint(
    molecule: &Molecule,
    n_bits: u32,
    target_size: u32,
    from_atoms: Option<&[usize]>,
    ignore_atoms: Option<&[usize]>,
    atom_invariants: Option<&[u32]>,
    include_chirality: bool,
) -> Result<SparseCountFingerprint, FingerprintError> {
    // RDKit source: AtomPairs.cpp lines 299-310
    // RDKit✔️✔️: SparseIntVect<boost::int64_t> *getHashedTopologicalTorsionFingerprint(
    // RDKit✔️✔️:     const ROMol &mol, unsigned int nBits, unsigned int targetSize,
    // RDKit✔️✔️:     const std::vector<std::uint32_t> *fromAtoms,
    // RDKit✔️✔️:     const std::vector<std::uint32_t> *ignoreAtoms,
    // RDKit✔️✔️:     const std::vector<std::uint32_t> *atomInvariants, bool includeChirality) {
    // RDKit❗✔️:   RDLog::deprecationWarning("please use TopologicalTorsionGenerator");
    // RDKit✔️✔️:   PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit✔️✔️:                "bad atomInvariants size");
    // RDKit✔️🔝:   auto *res = new SparseIntVect<boost::int64_t>(nBits);
    // RDKit✔️🔝:   TorsionFpCalc(res, mol, nBits, targetSize, fromAtoms, ignoreAtoms,
    // RDKit✔️🔝:                 atomInvariants, includeChirality);
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    // Rust's deprecation attribute reports use at compile time instead of at
    // runtime. `torsion_fp_calc` already returns the shared generator's owned
    // sparse vector, so this adapter preserves the source result while avoiding
    // the source wrapper's additional allocation and element-copy pass.
    torsion_fp_calc(
        molecule,
        n_bits,
        target_size,
        from_atoms,
        ignore_atoms,
        atom_invariants,
        include_chirality,
    )
}

#[deprecated(note = "please use TopologicalTorsionGenerator")]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub fn getHashedTopologicalTorsionFingerprintAsBitVect(
    molecule: &Molecule,
    n_bits: u32,
    target_size: u32,
    from_atoms: Option<&[usize]>,
    ignore_atoms: Option<&[usize]>,
    atom_invariants: Option<&[u32]>,
    n_bits_per_entry: u32,
    include_chirality: bool,
) -> Result<Fingerprint, FingerprintError> {
    // RDKit source: AtomPairs.cpp lines 312-347
    // RDKit✔️✔️: ExplicitBitVect *getHashedTopologicalTorsionFingerprintAsBitVect(
    // RDKit✔️✔️:     const ROMol &mol, unsigned int nBits, unsigned int targetSize,
    // RDKit✔️✔️:     const std::vector<std::uint32_t> *fromAtoms,
    // RDKit✔️✔️:     const std::vector<std::uint32_t> *ignoreAtoms,
    // RDKit✔️✔️:     const std::vector<std::uint32_t> *atomInvariants,
    // RDKit✔️✔️:     unsigned int nBitsPerEntry, bool includeChirality) {
    // RDKit❗✔️:   RDLog::deprecationWarning("please use TopologicalTorsionGenerator");
    // RDKit✔️✔️:   PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit✔️✔️:                "bad atomInvariants size");
    // RDKit✔️✔️:   static int bounds[4] = {1, 2, 4, 8};
    // RDKit✔️✔️:   unsigned int blockLength = nBits / nBitsPerEntry;
    // RDKit✔️✔️:   auto *sres = new SparseIntVect<boost::int64_t>(blockLength);
    // RDKit✔️✔️:   TorsionFpCalc(sres, mol, blockLength, targetSize, fromAtoms, ignoreAtoms,
    // RDKit✔️✔️:                 atomInvariants, includeChirality);
    // RDKit✔️✔️:   auto *res = new ExplicitBitVect(nBits);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (nBitsPerEntry != 4) {
    // RDKit✔️✔️:     for (auto val : sres->getNonzeroElements()) {
    // RDKit✔️✔️:       for (unsigned int i = 0; i < nBitsPerEntry; ++i) {
    // RDKit✔️✔️:         if (val.second > static_cast<int>(i)) {
    // RDKit✔️✔️:           res->setBit(val.first * nBitsPerEntry + i);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     for (auto val : sres->getNonzeroElements()) {
    // RDKit✔️✔️:       for (unsigned int i = 0; i < nBitsPerEntry; ++i) {
    // RDKit✔️✔️:         if (val.second >= bounds[i]) {
    // RDKit✔️✔️:           res->setBit(val.first * nBitsPerEntry + i);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   delete sres;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Rust's deprecation attribute reports use at compile time. Invalid zero
    // divisors and zero-sized hash blocks are returned as structured errors
    // instead of entering source-undefined arithmetic/generator states.
    if n_bits_per_entry == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "legacy topological torsion nBitsPerEntry must be greater than zero",
        });
    }
    let block_length = n_bits / n_bits_per_entry;
    if block_length == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "legacy topological torsion bit vector requires at least one complete entry block",
        });
    }

    let sparse_counts = torsion_fp_calc(
        molecule,
        block_length,
        target_size,
        from_atoms,
        ignore_atoms,
        atom_invariants,
        include_chirality,
    )?;
    let mut on_bits = Vec::new();
    const BOUNDS: [i32; 4] = [1, 2, 4, 8];
    for (&bit_id, &count) in sparse_counts.nonzero_elements() {
        for i in 0..n_bits_per_entry {
            let set = if n_bits_per_entry == 4 {
                count >= BOUNDS[i as usize]
            } else {
                i64::from(count) > i64::from(i)
            };
            if set {
                let output_bit = bit_id
                    .checked_mul(u64::from(n_bits_per_entry))
                    .and_then(|base| base.checked_add(u64::from(i)))
                    .ok_or(FingerprintError::InvalidArguments {
                        reason: "legacy topological torsion bit index overflow",
                    })?;
                on_bits.push(usize::try_from(output_bit).map_err(|_| {
                    FingerprintError::InvalidArguments {
                        reason: "legacy topological torsion bit index exceeds platform size",
                    }
                })?);
            }
        }
    }
    let output_size = usize::try_from(n_bits).map_err(|_| FingerprintError::InvalidArguments {
        reason: "legacy topological torsion bit-vector size exceeds platform size",
    })?;
    Ok(Fingerprint::from_on_bits(output_size, on_bits))
}

/// Rust-native parameters for the RDKit Topological Torsion generator.
///
/// This family is intentionally distinct from [`TopologicalFingerprintParams`],
/// which configures RDKit's path/subgraph `RDKFingerprintMol` algorithm.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologicalTorsionFingerprintParams {
    pub include_chirality: bool,
    pub torsion_atom_count: u32,
    pub count_simulation: bool,
    pub count_bounds: Vec<u32>,
    pub fp_size: u32,
    pub num_bits_per_feature: u32,
    pub only_shortest_paths: bool,
    pub from_atoms: Option<Vec<usize>>,
    pub ignore_atoms: Option<Vec<usize>>,
    pub custom_atom_invariants: Option<Vec<u32>>,
    pub atom_invariants_generator: Option<AtomPairAtomInvGenerator>,
}

impl Default for TopologicalTorsionFingerprintParams {
    fn default() -> Self {
        Self {
            include_chirality: false,
            torsion_atom_count: 4,
            count_simulation: true,
            count_bounds: vec![1, 2, 4, 8],
            fp_size: 2048,
            num_bits_per_feature: 1,
            only_shortest_paths: false,
            from_atoms: None,
            ignore_atoms: None,
            custom_atom_invariants: None,
            atom_invariants_generator: None,
        }
    }
}

impl TopologicalTorsionFingerprintParams {
    fn generator_arguments(&self) -> Result<TopologicalTorsionArguments, FingerprintError> {
        if self.num_bits_per_feature == 0 {
            return Err(FingerprintError::InvalidArguments {
                reason: "Topological Torsion num_bits_per_feature must be greater than zero",
            });
        }
        let mut arguments = TopologicalTorsionArguments::new(
            self.include_chirality,
            self.torsion_atom_count,
            self.count_simulation,
            self.count_bounds.clone(),
            self.fp_size,
        )?;
        arguments.fingerprint_arguments.d_num_bits_per_feature = self.num_bits_per_feature;
        arguments.df_only_shortest_paths = self.only_shortest_paths;
        Ok(arguments)
    }

    fn call_arguments(&self) -> FingerprintFuncArguments {
        FingerprintFuncArguments {
            from_atoms: self.from_atoms.clone(),
            ignore_atoms: self.ignore_atoms.clone(),
            custom_atom_invariants: self.custom_atom_invariants.clone(),
            ..FingerprintFuncArguments::default()
        }
    }
}

/// Vector form requested from [`topological_torsion_fingerprint_with_output`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TopologicalTorsionFingerprintVector {
    SparseCount,
    SparseBit,
    Count,
    #[default]
    Bit,
}

/// Typed request for a Topological Torsion vector and shared provenance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct TopologicalTorsionFingerprintOutputRequest {
    pub vector: TopologicalTorsionFingerprintVector,
    pub atom_to_bits: bool,
    pub atom_counts: bool,
    pub bit_paths: bool,
    pub atoms_per_bit: bool,
}

/// One of the four source-supported Topological Torsion vector forms.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TopologicalTorsionFingerprintValue {
    SparseCount(SparseCountFingerprint),
    SparseBit(SparseBitFingerprint),
    Count(SparseCountFingerprint),
    Bit(Fingerprint),
}

/// Rust-native Topological Torsion result with optional shared provenance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologicalTorsionFingerprintResult {
    pub fingerprint: TopologicalTorsionFingerprintValue,
    pub additional_output: Option<AdditionalOutput>,
}

/// Construct the public Rust Topological Torsion generator over the sole
/// shared fingerprint-generator core.
pub fn topological_torsion_generator(
    params: &TopologicalTorsionFingerprintParams,
) -> Result<TopologicalTorsionFingerprintGenerator, FingerprintError> {
    let arguments = params.generator_arguments()?;
    getTopologicalTorsionGenerator(&arguments, params.atom_invariants_generator.clone(), true)
}

fn allocate_topological_torsion_output(
    request: TopologicalTorsionFingerprintOutputRequest,
) -> Option<AdditionalOutput> {
    if !(request.atom_to_bits || request.atom_counts || request.bit_paths || request.atoms_per_bit)
    {
        return None;
    }
    let mut output = AdditionalOutput::new();
    if request.atom_to_bits {
        output.allocate_atom_to_bits();
    }
    if request.atom_counts {
        output.allocate_atom_counts();
    }
    if request.bit_paths {
        output.allocate_bit_paths();
    }
    if request.atoms_per_bit {
        output.allocate_atoms_per_bit();
    }
    Some(output)
}

/// Generate any source-supported Topological Torsion vector form and the
/// requested shared provenance without mutating `molecule`.
pub fn topological_torsion_fingerprint_with_output(
    molecule: &Molecule,
    params: &TopologicalTorsionFingerprintParams,
    request: TopologicalTorsionFingerprintOutputRequest,
) -> Result<TopologicalTorsionFingerprintResult, FingerprintError> {
    if params.fp_size == 0
        && matches!(
            request.vector,
            TopologicalTorsionFingerprintVector::SparseBit
                | TopologicalTorsionFingerprintVector::Bit
        )
    {
        return Err(FingerprintError::InvalidArguments {
            reason: "Topological Torsion bit-vector convenience output requires fp_size > 0",
        });
    }
    let generator = topological_torsion_generator(params)?;
    let mut arguments = params.call_arguments();
    arguments.additional_output = allocate_topological_torsion_output(request);
    let fingerprint = match request.vector {
        TopologicalTorsionFingerprintVector::SparseCount => {
            TopologicalTorsionFingerprintValue::SparseCount(
                generator.getSparseCountFingerprint(molecule, &mut arguments)?,
            )
        }
        TopologicalTorsionFingerprintVector::SparseBit => {
            TopologicalTorsionFingerprintValue::SparseBit(
                generator.getSparseFingerprint(molecule, &mut arguments)?,
            )
        }
        TopologicalTorsionFingerprintVector::Count => TopologicalTorsionFingerprintValue::Count(
            generator.getCountFingerprint(molecule, &mut arguments)?,
        ),
        TopologicalTorsionFingerprintVector::Bit => TopologicalTorsionFingerprintValue::Bit(
            generator.getFingerprint(molecule, &mut arguments)?,
        ),
    };
    Ok(TopologicalTorsionFingerprintResult {
        fingerprint,
        additional_output: arguments.additional_output,
    })
}

pub fn topological_torsion_sparse_count_fingerprint(
    molecule: &Molecule,
    params: &TopologicalTorsionFingerprintParams,
) -> Result<SparseCountFingerprint, FingerprintError> {
    let generator = topological_torsion_generator(params)?;
    generator.getSparseCountFingerprint(molecule, &mut params.call_arguments())
}

pub fn topological_torsion_sparse_fingerprint(
    molecule: &Molecule,
    params: &TopologicalTorsionFingerprintParams,
) -> Result<SparseBitFingerprint, FingerprintError> {
    if params.fp_size == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "Topological Torsion sparse-bit convenience output requires fp_size > 0",
        });
    }
    let generator = topological_torsion_generator(params)?;
    generator.getSparseFingerprint(molecule, &mut params.call_arguments())
}

pub fn topological_torsion_count_fingerprint(
    molecule: &Molecule,
    params: &TopologicalTorsionFingerprintParams,
) -> Result<SparseCountFingerprint, FingerprintError> {
    let generator = topological_torsion_generator(params)?;
    generator.getCountFingerprint(molecule, &mut params.call_arguments())
}

pub fn topological_torsion_fingerprint(
    molecule: &Molecule,
    params: &TopologicalTorsionFingerprintParams,
) -> Result<Fingerprint, FingerprintError> {
    if params.fp_size == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "Topological Torsion bit-vector convenience output requires fp_size > 0",
        });
    }
    let generator = topological_torsion_generator(params)?;
    generator.getFingerprint(molecule, &mut params.call_arguments())
}

/// Legacy compatibility form selected by [`TopologicalTorsionLegacyParams`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TopologicalTorsionLegacyKind {
    #[default]
    UnfoldedCount,
    HashedCount,
    HashedBit,
}

/// Rust-native parameters for the three deprecated RDKit compatibility paths.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologicalTorsionLegacyParams {
    pub kind: TopologicalTorsionLegacyKind,
    pub n_bits: u32,
    pub torsion_atom_count: u32,
    pub from_atoms: Option<Vec<usize>>,
    pub ignore_atoms: Option<Vec<usize>>,
    pub atom_invariants: Option<Vec<u32>>,
    pub n_bits_per_entry: u32,
    pub include_chirality: bool,
}

impl Default for TopologicalTorsionLegacyParams {
    fn default() -> Self {
        Self {
            kind: TopologicalTorsionLegacyKind::UnfoldedCount,
            n_bits: 2048,
            torsion_atom_count: 4,
            from_atoms: None,
            ignore_atoms: None,
            atom_invariants: None,
            n_bits_per_entry: 4,
            include_chirality: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TopologicalTorsionLegacyResult {
    SparseCount(SparseCountFingerprint),
    Bit(Fingerprint),
}

/// Call one of the three source-compatible legacy adapters. This function is
/// a naming/typing layer only and does not contain a second chemistry path.
#[allow(deprecated)]
pub fn topological_torsion_legacy_fingerprint(
    molecule: &Molecule,
    params: &TopologicalTorsionLegacyParams,
) -> Result<TopologicalTorsionLegacyResult, FingerprintError> {
    let result = match params.kind {
        TopologicalTorsionLegacyKind::UnfoldedCount => {
            TopologicalTorsionLegacyResult::SparseCount(getTopologicalTorsionFingerprint(
                molecule,
                params.torsion_atom_count,
                params.from_atoms.as_deref(),
                params.ignore_atoms.as_deref(),
                params.atom_invariants.as_deref(),
                params.include_chirality,
            )?)
        }
        TopologicalTorsionLegacyKind::HashedCount => {
            TopologicalTorsionLegacyResult::SparseCount(getHashedTopologicalTorsionFingerprint(
                molecule,
                params.n_bits,
                params.torsion_atom_count,
                params.from_atoms.as_deref(),
                params.ignore_atoms.as_deref(),
                params.atom_invariants.as_deref(),
                params.include_chirality,
            )?)
        }
        TopologicalTorsionLegacyKind::HashedBit => {
            TopologicalTorsionLegacyResult::Bit(getHashedTopologicalTorsionFingerprintAsBitVect(
                molecule,
                params.n_bits,
                params.torsion_atom_count,
                params.from_atoms.as_deref(),
                params.ignore_atoms.as_deref(),
                params.atom_invariants.as_deref(),
                params.n_bits_per_entry,
                params.include_chirality,
            )?)
        }
    };
    Ok(result)
}

#[allow(non_snake_case)]
pub fn generatorToJSON(generator: &TypedFingerprintGenerator) -> String {
    generator.to_json()
}

#[allow(non_snake_case)]
pub fn generatorFromJSON(json: &str) -> Result<TypedFingerprintGenerator, FingerprintError> {
    generator::generator_from_json(json)
}

#[allow(clippy::too_many_arguments, non_snake_case)]
pub fn getMorganGeneratorWithParams(
    radius: u32,
    count_simulation: bool,
    include_chirality: bool,
    use_bond_types: bool,
    only_nonzero_invariants: bool,
    include_redundant_environments: bool,
    atom_invariants_generator: Option<MorganAtomInvariantsGenerator>,
    bond_invariants_generator: Option<MorganBondInvariantsGenerator>,
    fp_size: u32,
    count_bounds: Vec<u32>,
    owns_atom_inv_gen: bool,
    owns_bond_inv_gen: bool,
) -> Result<MorganFingerprintGenerator, FingerprintError> {
    // RDKit source: MorganGenerator.cpp lines 504-520
    // RDKit✔️✔️: template <typename OutputType>
    // RDKit✔️✔️: FingerprintGenerator<OutputType> *getMorganGenerator(
    // RDKit✔️✔️:     unsigned int radius, bool countSimulation, bool includeChirality,
    // RDKit✔️✔️:     bool useBondTypes, bool onlyNonzeroInvariants,
    // RDKit✔️✔️:     bool includeRedundantEnvironments,
    // RDKit✔️✔️:     AtomInvariantsGenerator *atomInvariantsGenerator,
    // RDKit✔️✔️:     BondInvariantsGenerator *bondInvariantsGenerator, std::uint32_t fpSize,
    // RDKit✔️✔️:     std::vector<std::uint32_t> countBounds, bool ownsAtomInvGen,
    // RDKit✔️✔️:     bool ownsBondInvGen) {
    // RDKit✔️✔️:   MorganArguments arguments(radius, countSimulation, includeChirality,
    // RDKit✔️✔️:                             onlyNonzeroInvariants, countBounds, fpSize,
    // RDKit✔️✔️:                             includeRedundantEnvironments, useBondTypes);
    let arguments = MorganArguments::new(
        radius,
        count_simulation,
        include_chirality,
        only_nonzero_invariants,
        count_bounds,
        fp_size,
        include_redundant_environments,
        use_bond_types,
    )?;

    // RDKit✔️✔️:   return getMorganGenerator<OutputType>(arguments, atomInvariantsGenerator,
    // RDKit✔️✔️:                                         bondInvariantsGenerator, ownsAtomInvGen,
    // RDKit✔️✔️:                                         ownsBondInvGen);
    // RDKit✔️✔️: }
    Ok(getMorganGenerator(
        &arguments,
        atom_invariants_generator,
        bond_invariants_generator,
        owns_atom_inv_gen,
        owns_bond_inv_gen,
    ))
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganFingerprintOutput {
    pub fingerprint: Fingerprint,
    pub additional_output: Option<MorganAdditionalOutput>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganSparseFingerprintOutput {
    pub fingerprint: SparseCountFingerprint,
    pub atoms_setting_bits: Option<BTreeMap<usize, Vec<(usize, u32)>>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganBitFingerprintOutput {
    pub fingerprint: Fingerprint,
    pub atoms_setting_bits: Option<BTreeMap<usize, Vec<(usize, u32)>>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fingerprint {
    bits: Vec<u64>,
    n_bits: usize,
}

impl Fingerprint {
    pub(crate) fn from_lsb_bytes(n_bits: usize, bytes: &[u8]) -> Self {
        let mut bits = vec![0; n_bits.div_ceil(64)];
        for (byte_index, &byte) in bytes.iter().enumerate() {
            for bit_index in 0..8 {
                let bit = byte_index * 8 + bit_index;
                if bit < n_bits && byte & (1 << bit_index) != 0 {
                    bits[bit / 64] |= 1u64 << (bit % 64);
                }
            }
        }
        Self { bits, n_bits }
    }

    #[must_use]
    pub fn from_on_bits(n_bits: usize, on_bits: impl IntoIterator<Item = usize>) -> Self {
        let mut bits = vec![0; n_bits.div_ceil(64)];
        for bit in on_bits {
            assert!(
                bit < n_bits,
                "fingerprint bit {bit} is outside n_bits={n_bits}"
            );
            bits[bit / 64] |= 1u64 << (bit % 64);
        }
        Self { bits, n_bits }
    }

    #[must_use]
    pub fn n_bits(&self) -> usize {
        self.n_bits
    }

    #[must_use]
    pub fn on_bits(&self) -> Vec<usize> {
        let mut out = Vec::new();
        for (word_idx, word) in self.bits.iter().copied().enumerate() {
            let mut remaining = word;
            while remaining != 0 {
                let offset = remaining.trailing_zeros() as usize;
                let bit = word_idx * 64 + offset;
                if bit < self.n_bits {
                    out.push(bit);
                }
                remaining &= remaining - 1;
            }
        }
        out
    }

    pub fn tanimoto(&self, other: &Self) -> Result<f64, FingerprintError> {
        if self.n_bits != other.n_bits {
            return Err(FingerprintError::BitLengthMismatch {
                left: self.n_bits,
                right: other.n_bits,
            });
        }
        let mut intersection = 0u32;
        let mut union = 0u32;
        for (left, right) in self.bits.iter().zip(&other.bits) {
            intersection += (left & right).count_ones();
            union += (left | right).count_ones();
        }
        Ok(if union == 0 {
            0.0
        } else {
            f64::from(intersection) / f64::from(union)
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum FingerprintError {
    #[error("fingerprint requires n_bits > 0")]
    EmptyFingerprint,
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
    #[error("invalid fingerprint arguments: {reason}")]
    InvalidArguments { reason: &'static str },
    #[error("invalid fingerprint arguments JSON: {0}")]
    InvalidArgumentsJson(String),
    #[error("CIPLabeler failed while preparing Morgan fingerprint chirality: {reason}")]
    CipLabeler { reason: String },
    #[error("stereochemistry preparation failed while generating a fingerprint: {reason}")]
    StereoPreparation { reason: String },
    #[error("ring preparation failed while generating a fingerprint: {reason}")]
    RingPreparation { reason: String },
    #[error("AtomPair fingerprint generation failed: {reason}")]
    AtomPair { reason: String },
    #[error("sparse fingerprint index {index} is outside vector length {size}")]
    SparseIndexOutOfRange { index: u64, size: u64 },
    #[error(transparent)]
    Valence(#[from] crate::ValenceError),
    #[error("invalid SMARTS pattern '{pattern}': {reason}")]
    InvalidSmartsPattern { pattern: String, reason: String },
    #[error("unsupported fingerprint option {option}: {reason}")]
    UnsupportedOption {
        option: &'static str,
        reason: &'static str,
    },
    #[error("fingerprint parallel execution failed: {0}")]
    ParallelExecution(String),
    #[error("Avalon REACCS conversion failed: {reason}")]
    AvalonConversion { reason: String },
    #[error("fingerprint bit length mismatch: {left} != {right}")]
    BitLengthMismatch { left: usize, right: usize },
}

impl From<crate::RingFindingError> for FingerprintError {
    fn from(error: crate::RingFindingError) -> Self {
        Self::RingPreparation {
            reason: error.to_string(),
        }
    }
}

impl From<crate::chemistry::ciplabeler::CipLabelerError> for FingerprintError {
    fn from(error: crate::chemistry::ciplabeler::CipLabelerError) -> Self {
        Self::CipLabeler {
            reason: error.to_string(),
        }
    }
}

impl From<crate::StereoError> for FingerprintError {
    fn from(error: crate::StereoError) -> Self {
        Self::StereoPreparation {
            reason: error.to_string(),
        }
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
    /// Minimum bond-path length. Zero is rejected as `minPath==0`.
    pub min_path: u32,
    /// Maximum bond-path length. Values below `min_path` are rejected.
    pub max_path: u32,
    /// Explicit output width. Zero is rejected.
    pub fp_size: u32,
    /// Optional seeded source `atomCounts` vector. Values are incremented and
    /// returned without clearing the caller-provided seed.
    pub atom_counts: Option<Vec<u32>>,
    /// Optional projection mask, which must have exactly `fp_size` bits.
    pub set_only_bits: Option<Fingerprint>,
    /// Enumerate branched subgraphs when true and linear bond paths when false.
    pub branched_paths: bool,
    /// `None` is an absent source pointer; `Some(Vec::new())` is a present
    /// empty selection and therefore enumerates no paths.
    pub from_atoms: Option<Vec<u32>>,
}

impl Default for LayeredFingerprintParams {
    fn default() -> Self {
        // RDKit✔️✔️: const ROMol &mol, unsigned int layerFlags = 0xFFFFFFFF,
        // RDKit✔️✔️: unsigned int minPath = 1, unsigned int maxPath = 7,
        // RDKit✔️✔️: unsigned int fpSize = 2048, std::vector<unsigned int> *atomCounts = nullptr,
        // RDKit✔️✔️: ExplicitBitVect *setOnlyBits = nullptr, bool branchedPaths = true,
        // RDKit✔️✔️: const std::vector<std::uint32_t> *fromAtoms = nullptr);
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
    pub fn validate(&self) -> Result<(), FingerprintError> {
        if self.min_path == 0 {
            return Err(FingerprintError::InvalidArguments {
                reason: "minPath==0",
            });
        }
        if self.max_path < self.min_path {
            return Err(FingerprintError::InvalidArguments {
                reason: "maxPath<minPath",
            });
        }
        if self.fp_size == 0 {
            return Err(FingerprintError::InvalidArguments {
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

fn copied_atoms_setting_bits(
    args: &FingerprintFuncArguments,
) -> Option<BTreeMap<usize, Vec<(usize, u32)>>> {
    args.additional_output
        .as_ref()
        .and_then(|output| output.bit_info_map.as_ref())
        .map(|bit_info_map| {
            bit_info_map
                .iter()
                .map(|(&bit, entries)| {
                    (
                        bit as usize,
                        entries
                            .iter()
                            .map(|&(atom_id, layer)| (atom_id as usize, layer))
                            .collect(),
                    )
                })
                .collect()
        })
}

// ---------------------------------------------------------------------------
// Public API — true read-only descriptor computation, follows the same
// &Molecule convention as mol_to_smiles and symmetrize_sssr.
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
pub fn morgan_get_fingerprint(
    molecule: &Molecule,
    radius: u32,
    invariants: Option<Vec<u32>>,
    from_atoms: Option<Vec<usize>>,
    use_chirality: bool,
    use_bond_types: bool,
    use_counts: bool,
    only_nonzero_invariants: bool,
    collect_atoms_setting_bits: bool,
    include_redundant_environments: bool,
) -> Result<MorganSparseFingerprintOutput, FingerprintError> {
    // RDKit source: MorganFingerprints.cpp lines 45-83
    // RDKit✔️✔️: SparseIntVect<uint32_t> *getFingerprint(
    // RDKit✔️✔️:     const ROMol &mol, unsigned int radius, std::vector<uint32_t> *invariants,
    // RDKit✔️✔️:     const std::vector<uint32_t> *fromAtoms, bool useChirality,
    // RDKit✔️✔️:     bool useBondTypes, bool useCounts, bool onlyNonzeroInvariants,
    // RDKit✔️✔️:     BitInfoMap *atomsSettingBits, bool includeRedundantEnvironments) {
    // RDKit✔️✔️:   bool countSimulation = false;
    // RDKit✔️✔️:   std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
    // RDKit✔️✔️:       MorganFingerprint::getMorganGenerator<std::uint32_t>(
    // RDKit✔️✔️:           radius, countSimulation, useChirality, useBondTypes,
    // RDKit✔️✔️:           onlyNonzeroInvariants, includeRedundantEnvironments));
    let generator = getMorganGeneratorWithParams(
        radius,
        false,
        use_chirality,
        use_bond_types,
        only_nonzero_invariants,
        include_redundant_environments,
        None,
        None,
        2048,
        vec![1, 2, 4, 8],
        false,
        false,
    )?;

    // RDKit✔️✔️:   RDKit::FingerprintFuncArguments args;
    // RDKit✔️✔️:   args.fromAtoms = fromAtoms;
    // RDKit✔️✔️:   args.customAtomInvariants = invariants;
    // RDKit✔️✔️:   AdditionalOutput ao;
    // RDKit✔️✔️:   if (atomsSettingBits) {
    // RDKit✔️✔️:     args.additionalOutput = &ao;
    // RDKit✔️✔️:     ao.allocateBitInfoMap();
    // RDKit✔️✔️:   }
    let mut args = FingerprintFuncArguments {
        from_atoms,
        custom_atom_invariants: invariants,
        ..Default::default()
    };
    if collect_atoms_setting_bits {
        let mut additional_output = AdditionalOutput::new();
        additional_output.allocate_bit_info_map();
        args.additional_output = Some(additional_output);
    }

    // RDKit✔️✔️:   SparseIntVect<uint32_t> *res;
    // RDKit✔️✔️:   if (!useCounts) {
    // RDKit✔️✔️:     auto tmp = fpgen->getSparseFingerprint(mol, args);
    // RDKit✔️✔️:     res = new SparseIntVect<uint32_t>(std::numeric_limits<uint32_t>::max());
    // RDKit✔️✔️:     for (auto idx : *(tmp->dp_bits)) {
    // RDKit✔️✔️:       res->setVal(idx, 1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = fpgen->getSparseCountFingerprint(mol, args).release();
    // RDKit✔️✔️:   }
    let fingerprint = if use_counts {
        generator.getSparseCountFingerprint(molecule, &mut args)?
    } else {
        let sparse_bits = generator.getSparseFingerprint(molecule, &mut args)?;
        let mut result = SparseCountFingerprint::new(u64::from(u32::MAX));
        for &bit_id in sparse_bits.on_bits() {
            result.set_val(bit_id, 1);
        }
        result
    };

    // RDKit✔️✔️:   if (atomsSettingBits) {
    // RDKit✔️✔️:     atomsSettingBits->clear();
    // RDKit✔️✔️:     for (const auto &pr : *(ao.bitInfoMap)) {
    // RDKit✔️✔️:       (*atomsSettingBits)[pr.first] = pr.second;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let atoms_setting_bits = if collect_atoms_setting_bits {
        copied_atoms_setting_bits(&args)
    } else {
        None
    };

    Ok(MorganSparseFingerprintOutput {
        fingerprint,
        atoms_setting_bits,
    })
}

#[allow(clippy::too_many_arguments)]
pub fn morgan_get_hashed_fingerprint(
    molecule: &Molecule,
    radius: u32,
    n_bits: u32,
    invariants: Option<Vec<u32>>,
    from_atoms: Option<Vec<usize>>,
    use_chirality: bool,
    use_bond_types: bool,
    only_nonzero_invariants: bool,
    collect_atoms_setting_bits: bool,
    include_redundant_environments: bool,
) -> Result<MorganSparseFingerprintOutput, FingerprintError> {
    // RDKit source: MorganFingerprints.cpp lines 85-119
    // RDKit✔️✔️: SparseIntVect<uint32_t> *getHashedFingerprint(
    // RDKit✔️✔️:     const ROMol &mol, unsigned int radius, unsigned int nBits,
    // RDKit✔️✔️:     std::vector<uint32_t> *invariants, const std::vector<uint32_t> *fromAtoms,
    // RDKit✔️✔️:     bool useChirality, bool useBondTypes, bool onlyNonzeroInvariants,
    // RDKit✔️✔️:     BitInfoMap *atomsSettingBits, bool includeRedundantEnvironments) {
    // RDKit✔️✔️:   if (nBits == 0) {
    // RDKit✔️✔️:     throw ValueErrorException("nBits can not be zero");
    // RDKit✔️✔️:   }
    if n_bits == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "nBits can not be zero",
        });
    }

    // RDKit✔️✔️:   bool countSimulation = false;
    // RDKit✔️✔️:   std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
    // RDKit✔️✔️:       MorganFingerprint::getMorganGenerator<std::uint32_t>(
    // RDKit✔️✔️:           radius, countSimulation, useChirality, useBondTypes,
    // RDKit✔️✔️:           onlyNonzeroInvariants, includeRedundantEnvironments, nullptr, nullptr,
    // RDKit✔️✔️:           nBits));
    let generator = getMorganGeneratorWithParams(
        radius,
        false,
        use_chirality,
        use_bond_types,
        only_nonzero_invariants,
        include_redundant_environments,
        None,
        None,
        n_bits,
        vec![1, 2, 4, 8],
        false,
        false,
    )?;

    // RDKit✔️✔️:   RDKit::FingerprintFuncArguments args;
    // RDKit✔️✔️:   args.fromAtoms = fromAtoms;
    // RDKit✔️✔️:   args.customAtomInvariants = invariants;
    // RDKit✔️✔️:   AdditionalOutput ao;
    // RDKit✔️✔️:   if (atomsSettingBits) {
    // RDKit✔️✔️:     args.additionalOutput = &ao;
    // RDKit✔️✔️:     ao.allocateBitInfoMap();
    // RDKit✔️✔️:   }
    let mut args = FingerprintFuncArguments {
        from_atoms,
        custom_atom_invariants: invariants,
        ..Default::default()
    };
    if collect_atoms_setting_bits {
        let mut additional_output = AdditionalOutput::new();
        additional_output.allocate_bit_info_map();
        args.additional_output = Some(additional_output);
    }

    // RDKit✔️✔️:   auto res = fpgen->getCountFingerprint(mol, args).release();
    let fingerprint = generator.getCountFingerprint(molecule, &mut args)?;

    // RDKit✔️✔️:   if (atomsSettingBits) {
    // RDKit✔️✔️:     atomsSettingBits->clear();
    // RDKit✔️✔️:     for (const auto &pr : *(ao.bitInfoMap)) {
    // RDKit✔️✔️:       (*atomsSettingBits)[pr.first] = pr.second;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let atoms_setting_bits = if collect_atoms_setting_bits {
        copied_atoms_setting_bits(&args)
    } else {
        None
    };

    Ok(MorganSparseFingerprintOutput {
        fingerprint,
        atoms_setting_bits,
    })
}

#[allow(clippy::too_many_arguments)]
pub fn morgan_get_fingerprint_as_bit_vect(
    molecule: &Molecule,
    radius: u32,
    n_bits: u32,
    invariants: Option<Vec<u32>>,
    from_atoms: Option<Vec<usize>>,
    use_chirality: bool,
    use_bond_types: bool,
    only_nonzero_invariants: bool,
    collect_atoms_setting_bits: bool,
    include_redundant_environments: bool,
) -> Result<MorganBitFingerprintOutput, FingerprintError> {
    // RDKit source: MorganFingerprints.cpp lines 121-155
    // RDKit✔️✔️: ExplicitBitVect *getFingerprintAsBitVect(
    // RDKit✔️✔️:     const ROMol &mol, unsigned int radius, unsigned int nBits,
    // RDKit✔️✔️:     std::vector<uint32_t> *invariants, const std::vector<uint32_t> *fromAtoms,
    // RDKit✔️✔️:     bool useChirality, bool useBondTypes, bool onlyNonzeroInvariants,
    // RDKit✔️✔️:     BitInfoMap *atomsSettingBits, bool includeRedundantEnvironments) {
    // RDKit✔️✔️:   if (nBits == 0) {
    // RDKit✔️✔️:     throw ValueErrorException("nBits can not be zero");
    // RDKit✔️✔️:   }
    if n_bits == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "nBits can not be zero",
        });
    }

    // RDKit✔️✔️:   bool countSimulation = false;
    // RDKit✔️✔️:   std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
    // RDKit✔️✔️:       MorganFingerprint::getMorganGenerator<std::uint32_t>(
    // RDKit✔️✔️:           radius, countSimulation, useChirality, useBondTypes,
    // RDKit✔️✔️:           onlyNonzeroInvariants, includeRedundantEnvironments, nullptr, nullptr,
    // RDKit✔️✔️:           nBits));
    let generator = getMorganGeneratorWithParams(
        radius,
        false,
        use_chirality,
        use_bond_types,
        only_nonzero_invariants,
        include_redundant_environments,
        None,
        None,
        n_bits,
        vec![1, 2, 4, 8],
        false,
        false,
    )?;

    // RDKit✔️✔️:   RDKit::FingerprintFuncArguments args;
    // RDKit✔️✔️:   args.fromAtoms = fromAtoms;
    // RDKit✔️✔️:   args.customAtomInvariants = invariants;
    // RDKit✔️✔️:   AdditionalOutput ao;
    // RDKit✔️✔️:   if (atomsSettingBits) {
    // RDKit✔️✔️:     args.additionalOutput = &ao;
    // RDKit✔️✔️:     ao.allocateBitInfoMap();
    // RDKit✔️✔️:   }
    let mut args = FingerprintFuncArguments {
        from_atoms,
        custom_atom_invariants: invariants,
        ..Default::default()
    };
    if collect_atoms_setting_bits {
        let mut additional_output = AdditionalOutput::new();
        additional_output.allocate_bit_info_map();
        args.additional_output = Some(additional_output);
    }

    // RDKit✔️✔️:   auto res = fpgen->getFingerprint(mol, args).release();
    let fingerprint = generator.getFingerprint(molecule, &mut args)?;

    // RDKit✔️✔️:   if (atomsSettingBits) {
    // RDKit✔️✔️:     atomsSettingBits->clear();
    // RDKit✔️✔️:     for (const auto &pr : *(ao.bitInfoMap)) {
    // RDKit✔️✔️:       (*atomsSettingBits)[pr.first] = pr.second;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let atoms_setting_bits = if collect_atoms_setting_bits {
        copied_atoms_setting_bits(&args)
    } else {
        None
    };

    Ok(MorganBitFingerprintOutput {
        fingerprint,
        atoms_setting_bits,
    })
}

pub fn morgan_fingerprint(
    molecule: &Molecule,
    params: &MorganFingerprintParams,
) -> Result<Fingerprint, FingerprintError> {
    let output = morgan_fingerprint_with_output(molecule, params)?;
    Ok(output.fingerprint)
}

pub fn morgan_fingerprint_with_output(
    molecule: &Molecule,
    params: &MorganFingerprintParams,
) -> Result<MorganFingerprintOutput, FingerprintError> {
    validate_morgan_params(params)?;
    if params.n_bits == 0 {
        return Err(FingerprintError::EmptyFingerprint);
    }
    let atom_invariants_generator = match params.atom_invariants_generator.clone() {
        MorganAtomInvariantsGenerator::Connectivity {
            include_ring_membership,
        } => Some(MorganAtomInvariantsGenerator::Connectivity {
            include_ring_membership,
        }),
        MorganAtomInvariantsGenerator::Feature => Some(MorganAtomInvariantsGenerator::Feature),
    };
    let bond_invariants_generator =
        Some(params.bond_invariants_generator.clone().unwrap_or_else(|| {
            MorganBondInvariantsGenerator {
                use_bond_types: params.use_bond_types,
                use_chirality: params.use_chirality,
            }
        }));
    let mut generator = getMorganGeneratorWithParams(
        params.radius,
        params.count_simulation,
        params.use_chirality,
        params.use_bond_types,
        params.only_nonzero_invariants,
        params.include_redundant_environments,
        atom_invariants_generator,
        bond_invariants_generator,
        params.n_bits as u32,
        params.count_bounds.clone(),
        true,
        true,
    )?;
    // RDKit✔️✔️: fpgen->getOptions()->d_numBitsPerFeature = nBitsPerHash;
    generator
        .fingerprint_arguments
        .fingerprint_arguments
        .d_num_bits_per_feature = params.num_bits_per_feature;

    let mut args = FingerprintFuncArguments {
        from_atoms: params.from_atoms.clone(),
        ignore_atoms: params.ignore_atoms.clone(),
        custom_atom_invariants: params.custom_atom_invariants.clone(),
        custom_bond_invariants: params.custom_bond_invariants.clone(),
        ..Default::default()
    };
    if params.collect_additional_output {
        let mut additional_output = AdditionalOutput::new();
        additional_output.allocate_atom_counts();
        additional_output.allocate_atom_to_bits();
        additional_output.allocate_bit_info_map();
        additional_output.allocate_atoms_per_bit();
        args.additional_output = Some(additional_output);
    }

    let fingerprint = generator.getFingerprint(molecule, &mut args)?;
    let additional_output = args
        .additional_output
        .map(morgan_additional_output_from_rdkit_output);

    Ok(MorganFingerprintOutput {
        fingerprint,
        additional_output,
    })
}

// ---------------------------------------------------------------------------
// Initial invariants
// ---------------------------------------------------------------------------

// RDKit✔️✔️: void getConnectivityInvariants(const ROMol &mol,
// RDKit✔️✔️:     std::vector<uint32_t> &invars,
// RDKit✔️✔️:     bool includeRingMembership) {
// RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
// RDKit✔️✔️:   PRECONDITION(invars.size() >= nAtoms, "vector too small");
// RDKit✔️✔️:   gboost::hash<std::vector<uint32_t>> vectHasher;
// RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
// RDKit✔️✔️:     Atom const *atom = mol.getAtomWithIdx(i);
// RDKit✔️✔️:     std::vector<uint32_t> components;
// RDKit✔️✔️:     components.push_back(atom->getAtomicNum());
// RDKit✔️✔️:     components.push_back(atom->getTotalDegree());
// RDKit✔️✔️:     components.push_back(atom->getTotalNumHs(true));
// RDKit✔️✔️:     components.push_back(atom->getFormalCharge());
// RDKit✔️✔️:     int deltaMass = static_cast<int>(
// RDKit✔️✔️:         atom->getMass() -
// RDKit✔️✔️:         PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
// RDKit✔️✔️:     components.push_back(deltaMass);
// RDKit✔️✔️:     if (includeRingMembership &&
// RDKit✔️✔️:         atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx())) {
// RDKit✔️✔️:       components.push_back(1);
// RDKit✔️✔️:     }
// RDKit✔️✔️:     invars[i] = vectHasher(components);
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT FUNCTION getConnectivityInvariants
//
// COSMolKit uses the same component vector and hashing algorithm as RDKit:
// atomic number, total degree, total Hs, formal charge, delta mass,
// optional ring membership flag, all hashed via boost::hash_combine.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn compute_connectivity_invariants(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    params: &MorganFingerprintParams,
) -> Vec<u32> {
    let num_atoms = molecule.num_atoms();

    let ring_info = if params.include_ring_membership {
        molecule.derived_cache().rings.as_ref()
    } else {
        None
    };

    let valence = molecule.derived_cache().valence.as_ref();

    let mut invariants = Vec::with_capacity(num_atoms);
    for i in 0..num_atoms {
        let atom = &molecule.atoms()[i];
        let degree = adjacency.neighbors_of(i).len() as u32;

        let implicit_hs = valence
            .and_then(|v| v.implicit_hydrogens.get(i).copied())
            .unwrap_or(0) as u32;

        // RDKit source: Atom.cpp lines 277-289
        // RDKit✔️✔️: unsigned int Atom::getTotalDegree() const {
        // RDKit✔️✔️:   unsigned int res = this->getTotalNumHs(false) + this->getDegree();
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
        // RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
        // RDKit✔️✔️:   if (includeNeighbors && dp_mol) {
        // RDKit✔️✔️:     auto nbrs = dp_mol->atomNeighbors(this);
        // RDKit✔️✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
        // RDKit✔️✔️:       return (nbr->getAtomicNum() == 1);
        // RDKit✔️✔️:     });
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        let attached_hs = atom.explicit_hydrogens() as u32 + implicit_hs;
        let neighbor_hs = adjacency
            .neighbors_of(i)
            .iter()
            .filter(|neighbor| molecule.atoms()[neighbor.atom_index].atomic_number() == 1)
            .count() as u32;
        let total_degree = degree + attached_hs;
        let total_hs = attached_hs + neighbor_hs;

        // RDKit's deltaMass = static_cast<int>(atom->getMass() - atomicWeight).
        let mass = rdkit_atomic_mass(atom.atomic_number(), atom.isotope());
        let atomic_weight = rdkit_atomic_mass(atom.atomic_number(), None);
        let delta_mass = (mass - atomic_weight).trunc() as i32;

        // Build the same component vector as RDKit, then hash via hash_combine.
        let mut components: Vec<u32> = Vec::with_capacity(6);
        components.push(atom.atomic_number() as u32);
        components.push(total_degree);
        components.push(total_hs);
        components.push(atom.formal_charge() as u32);
        components.push(delta_mass as u32);

        if let Some(rings) = ring_info {
            if rings.num_atom_rings(AtomId::new(i)) > 0 {
                components.push(1);
            }
        }

        // gboost::hash<std::vector<uint32_t>> uses successive hash_combine calls.
        let mut inv = 0u32;
        for &c in &components {
            hash_combine(&mut inv, c);
        }
        invariants.push(inv);
    }
    invariants
}

// RDKit source: FingerprintUtil.cpp lines 135-160 (getFeatureInvariants)
// RDKit✔️✔️: void getFeatureInvariants(const ROMol &mol,
// RDKit✔️✔️:     std::vector<uint32_t> &invars,
// RDKit✔️✔️:     const std::vector<const ROMol *> *patterns) {
// RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
// RDKit✔️✔️:   PRECONDITION(invars.size() >= nAtoms, "vector too small");
// RDKit✔️✔️:
// RDKit✔️✔️:   auto useLocalPatterns = patterns == nullptr;
// RDKit✔️✔️:   std::vector<const ROMol *> featureMatchers;
// RDKit✔️✔️:   if (useLocalPatterns) {
// RDKit✔️✔️:     featureMatchers.reserve(defaultFeatureSmarts.size());
// RDKit✔️✔️:     for (const auto &smaIt : defaultFeatureSmarts) {
// RDKit✔️✔️:       const ROMol *matcher = pattern_flyweight(smaIt).get().getMatcher();
// RDKit✔️✔️:       CHECK_INVARIANT(matcher, "bad smarts");
// RDKit✔️✔️:       featureMatchers.push_back(matcher);
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   std::fill(invars.begin(), invars.end(), 0);
// RDKit✔️✔️:   auto &queries = (useLocalPatterns ? featureMatchers : *patterns);
// RDKit✔️✔️:   for (unsigned int i = 0; i < queries.size(); ++i) {
// RDKit✔️✔️:     unsigned int mask = 1 << i;
// RDKit✔️✔️:     std::vector<MatchVectType> matchVect;
// RDKit✔️✔️:     // to maintain thread safety, we have to copy the pattern
// RDKit✔️✔️:     // molecules:
// RDKit✔️✔️:     SubstructMatch(mol, ROMol(*queries[i], true), matchVect);
// RDKit✔️✔️:     for (const auto &mvIt : matchVect) {
// RDKit✔️✔️:       for (const auto &mIt : mvIt) {
// RDKit✔️✔️:         invars[mIt.second] |= mask;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }  // end of getFeatureInvariants()
//
// COSMolKit uses the same 6 SMARTS patterns and source-backed substructure
// matching instead of heuristic element/property classification.
fn compute_feature_invariants_with_patterns(
    molecule: &Molecule,
    supplied_patterns: Option<&[SsMatcher]>,
) -> Result<Vec<u32>, FingerprintError> {
    let num_atoms = molecule.num_atoms();
    let mut invariants = vec![0u32; num_atoms];
    let feature_matchers = match supplied_patterns {
        Some(patterns) => patterns,
        None => cached_default_feature_matchers()?,
    };
    for (feature_idx, matcher) in feature_matchers.iter().enumerate() {
        let query = matcher.getMatcher();
        let mask = 1u32 << feature_idx;
        for matched in crate::get_substruct_matches(molecule, query) {
            for atom_idx in matched.atom_mapping {
                invariants[atom_idx] |= mask;
            }
        }
    }
    Ok(invariants)
}

fn compute_feature_invariants(molecule: &Molecule) -> Result<Vec<u32>, FingerprintError> {
    MorganFeatureAtomInvGenerator::new().getAtomInvariants(molecule)
}

// ---------------------------------------------------------------------------
// Invariant propagation — hash_combine matches boost::hash_combine pattern
// ---------------------------------------------------------------------------

fn validate_morgan_params(params: &MorganFingerprintParams) -> Result<(), FingerprintError> {
    // Morgan options exposed through MorganFingerprintParams are checked by the
    // RDKit branch-matrix parity target. Source branches not represented by
    // this parameter object, such as JSON-supplied feature SMARTS patterns,
    // fail closed at their helper boundary instead of returning local guesses.
    let _ = params;
    Ok(())
}

fn rdkit_bond_type_code(order: crate::BondOrder) -> u32 {
    match order {
        crate::BondOrder::Null | crate::BondOrder::Unspecified => 0,
        crate::BondOrder::Single => 1,
        crate::BondOrder::Double => 2,
        crate::BondOrder::Triple => 3,
        crate::BondOrder::Quadruple => 4,
        crate::BondOrder::Quintuple => 5,
        crate::BondOrder::Hextuple => 6,
        crate::BondOrder::OneAndHalf => 7,
        crate::BondOrder::TwoAndHalf => 8,
        crate::BondOrder::ThreeAndHalf => 9,
        crate::BondOrder::FourAndHalf => 10,
        crate::BondOrder::FiveAndHalf => 11,
        crate::BondOrder::Aromatic => 12,
        crate::BondOrder::Ionic => 13,
        crate::BondOrder::Hydrogen => 14,
        crate::BondOrder::ThreeCenter => 15,
        crate::BondOrder::DativeOne => 16,
        crate::BondOrder::Dative => 17,
        crate::BondOrder::DativeLeft => 18,
        crate::BondOrder::DativeRight => 19,
        crate::BondOrder::Other => 20,
        crate::BondOrder::Zero => 21,
    }
}

fn rdkit_bond_stereo_code(stereo: crate::BondStereo) -> u32 {
    match stereo {
        crate::BondStereo::None => 0,
        crate::BondStereo::Any => 1,
        crate::BondStereo::Z => 2,
        crate::BondStereo::E => 3,
        crate::BondStereo::Cis => 4,
        crate::BondStereo::Trans => 5,
        crate::BondStereo::AtropCw => 6,
        crate::BondStereo::AtropCcw => 7,
    }
}

pub(crate) fn hash_combine(seed: &mut u32, value: u32) {
    // BEGIN RDKIT CPP FUNCTION hash_combine
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: inline void hash_combine(std::hash_result_t& seed, T const& v)
    // RDKit✔️✔️: {
    // RDKit✔️✔️:   gboost::hash<T> hasher;
    // RDKit✔️✔️:   seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION hash_combine
    // `value` is the identity hash for the unsigned-int inputs used by the
    // fingerprint paths. Local complexity review: both implementations make
    // the same fixed number of 32-bit scalar operations in O(1), allocate
    // nothing, and explicitly wrap modulo 2^32.
    *seed ^= value
        .wrapping_add(0x9e3779b9u32)
        .wrapping_add(seed.wrapping_shl(6))
        .wrapping_add(seed.wrapping_shr(2));
}

pub(crate) fn hash_range(values: &[u32]) -> u32 {
    // BEGIN RDKIT CPP TYPE hash_result_t
    // RDKit✔️✔️: namespace std {
    // RDKit✔️✔️: typedef std::uint32_t hash_result_t;
    // RDKit✔️✔️: }
    // END RDKIT CPP TYPE hash_result_t
    // BEGIN RDKIT CPP FUNCTION hash_range
    // RDKit✔️✔️: template <class It>
    // RDKit✔️✔️: inline std::hash_result_t hash_range(It first, It last) {
    // RDKit✔️✔️:   std::hash_result_t seed = 0;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (; first != last; ++first) {
    // RDKit✔️✔️:     hash_combine(seed, *first);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return seed;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION hash_range
    // Local complexity review: both implementations make one ordered O(n)
    // pass, retain one 32-bit seed, and allocate or clone nothing. A slice
    // iterator has the same contiguous access pattern as the source vector
    // iterators; the explicit u32 result prevents native pointer width from
    // changing wrapping behavior.
    let mut seed = 0u32;
    for &value in values {
        hash_combine(&mut seed, value);
    }
    seed
}

fn morgan_additional_output_from_rdkit_output(output: AdditionalOutput) -> MorganAdditionalOutput {
    let atom_counts = output.atom_counts.unwrap_or_default();
    let atom_to_bits = output
        .atom_to_bits
        .unwrap_or_default()
        .into_iter()
        .map(|bits| bits.into_iter().map(|bit| bit as usize).collect())
        .collect();
    let bit_info_map = output
        .bit_info_map
        .unwrap_or_default()
        .into_iter()
        .map(|(bit, entries)| {
            (
                bit as usize,
                entries
                    .into_iter()
                    .map(|(atom_id, layer)| (atom_id as usize, layer))
                    .collect(),
            )
        })
        .collect();
    let atoms_per_bit = output
        .atoms_per_bit
        .unwrap_or_default()
        .into_iter()
        .map(|(bit, atoms)| (bit as usize, atoms))
        .collect();
    MorganAdditionalOutput {
        atom_counts,
        atom_to_bits,
        bit_info_map,
        atoms_per_bit,
    }
}

fn additional_output_allocation_mismatch(
    left: bool,
    right: bool,
    reason: &'static str,
) -> Result<(), FingerprintError> {
    if left ^ right {
        Err(FingerprintError::InvalidArguments { reason })
    } else {
        Ok(())
    }
}

fn duplicate_additional_output_bit(
    old_output: &AdditionalOutput,
    new_output: &mut AdditionalOutput,
    orig_bit_id: u64,
    new_bit_id: u64,
) -> Result<(), FingerprintError> {
    // RDKit source: FingerprintGenerator.cpp lines 438-479
    // RDKit✔️✔️: template <typename OutputType>
    // RDKit✔️✔️: void duplicateAdditionalOutputBit(AdditionalOutput &oldAO,
    // RDKit✔️✔️:                                   AdditionalOutput &newAO, OutputType origBitId,
    // RDKit✔️✔️:                                   OutputType newBitId) {
    // RDKit✔️✔️:   PRECONDITION(!((oldAO.bitInfoMap != nullptr) ^ (newAO.bitInfoMap != nullptr)),
    // RDKit✔️✔️:                "bitInfoMap not allocated");
    // RDKit✔️✔️:   PRECONDITION(!((oldAO.atomToBits != nullptr) ^ (newAO.atomToBits != nullptr)),
    // RDKit✔️✔️:                "atomToBits not allocated");
    // RDKit✔️✔️:   PRECONDITION(!((oldAO.bitPaths != nullptr) ^ (newAO.bitPaths != nullptr)),
    // RDKit✔️✔️:                "bitPaths not allocated");
    additional_output_allocation_mismatch(
        old_output.bit_info_map.is_some(),
        new_output.bit_info_map.is_some(),
        "bitInfoMap not allocated",
    )?;
    additional_output_allocation_mismatch(
        old_output.atom_to_bits.is_some(),
        new_output.atom_to_bits.is_some(),
        "atomToBits not allocated",
    )?;
    additional_output_allocation_mismatch(
        old_output.bit_paths.is_some(),
        new_output.bit_paths.is_some(),
        "bitPaths not allocated",
    )?;

    // RDKit✔️✔️:   // we don't need to do anything with atomCounts

    // RDKit✔️✔️:   if (oldAO.atomToBits) {
    // RDKit✔️✔️:     if (newAO.atomToBits->empty()) {
    // RDKit✔️✔️:       newAO.atomToBits->resize(oldAO.atomToBits->size());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned int i = 0; i < oldAO.atomToBits->size(); ++i) {
    // RDKit✔️✔️:       const auto &nv = oldAO.atomToBits->at(i);
    // RDKit✔️✔️:       if (std::find(nv.begin(), nv.end(), origBitId) != nv.end()) {
    // RDKit✔️✔️:         newAO.atomToBits->at(i).push_back(newBitId);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if let (Some(old_atom_to_bits), Some(new_atom_to_bits)) = (
        old_output.atom_to_bits.as_ref(),
        new_output.atom_to_bits.as_mut(),
    ) {
        if new_atom_to_bits.is_empty() {
            new_atom_to_bits.resize_with(old_atom_to_bits.len(), Vec::new);
        }
        for (idx, old_bits) in old_atom_to_bits.iter().enumerate() {
            if old_bits.contains(&orig_bit_id) {
                new_atom_to_bits[idx].push(new_bit_id);
            }
        }
    }

    // RDKit✔️✔️:   if (oldAO.bitInfoMap) {
    // RDKit✔️✔️:     const auto v = oldAO.bitInfoMap->find(origBitId);
    // RDKit✔️✔️:     if (v != oldAO.bitInfoMap->end()) {
    // RDKit✔️✔️:       (*newAO.bitInfoMap)[newBitId] = v->second;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if let (Some(old_bit_info_map), Some(new_bit_info_map)) = (
        old_output.bit_info_map.as_ref(),
        new_output.bit_info_map.as_mut(),
    ) {
        if let Some(value) = old_bit_info_map.get(&orig_bit_id) {
            new_bit_info_map.insert(new_bit_id, value.clone());
        }
    }

    // RDKit✔️✔️:   if (oldAO.bitPaths) {
    // RDKit✔️✔️:     const auto v = oldAO.bitPaths->find(origBitId);
    // RDKit✔️✔️:     if (v != oldAO.bitPaths->end()) {
    // RDKit✔️✔️:       (*newAO.bitPaths)[newBitId] = v->second;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if let (Some(old_bit_paths), Some(new_bit_paths)) =
        (old_output.bit_paths.as_ref(), new_output.bit_paths.as_mut())
    {
        if let Some(value) = old_bit_paths.get(&orig_bit_id) {
            new_bit_paths.insert(new_bit_id, value.clone());
        }
    }

    // RDKit✔️✔️:   if (oldAO.atomsPerBit) {
    // RDKit✔️✔️:     const auto v = oldAO.atomsPerBit->find(origBitId);
    // RDKit✔️✔️:     if (v != oldAO.atomsPerBit->end()) {
    // RDKit✔️✔️:       (*newAO.atomsPerBit)[newBitId] = v->second;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if let (Some(old_atoms_per_bit), Some(new_atoms_per_bit)) = (
        old_output.atoms_per_bit.as_ref(),
        new_output.atoms_per_bit.as_mut(),
    ) {
        if let Some(value) = old_atoms_per_bit.get(&orig_bit_id) {
            new_atoms_per_bit.insert(new_bit_id, value.clone());
        }
    }

    Ok(())
}

fn setup_temp_additional_output(
    args: &mut FingerprintFuncArguments,
    count_simulation_output: &mut AdditionalOutput,
    num_atoms: usize,
) {
    // RDKit source: FingerprintGenerator.cpp lines 481-499
    // RDKit✔️✔️: void setupTempAdditionalOutput(RDKit::FingerprintFuncArguments &args,
    // RDKit✔️✔️:                                AdditionalOutput &countSimulationOutput,
    // RDKit✔️✔️:                                size_t numAtoms) {
    let Some(additional_output) = args.additional_output.as_mut() else {
        return;
    };

    // RDKit✔️✔️:   if (args.additionalOutput->atomToBits) {
    // RDKit✔️✔️:     countSimulationOutput.allocateAtomToBits();
    // RDKit✔️✔️:   }
    if additional_output.atom_to_bits.is_some() {
        count_simulation_output.allocate_atom_to_bits();
    }
    // RDKit✔️✔️:   if (args.additionalOutput->atomCounts) {
    // RDKit✔️✔️:     countSimulationOutput.allocateAtomCounts();
    // RDKit✔️✔️:   }
    if additional_output.atom_counts.is_some() {
        count_simulation_output.allocate_atom_counts();
    }
    // RDKit✔️✔️:   if (args.additionalOutput->bitInfoMap) {
    // RDKit✔️✔️:     countSimulationOutput.allocateBitInfoMap();
    // RDKit✔️✔️:   }
    if additional_output.bit_info_map.is_some() {
        count_simulation_output.allocate_bit_info_map();
    }
    // RDKit✔️✔️:   if (args.additionalOutput->bitPaths) {
    // RDKit✔️✔️:     countSimulationOutput.allocateBitPaths();
    // RDKit✔️✔️:   }
    if additional_output.bit_paths.is_some() {
        count_simulation_output.allocate_bit_paths();
    }
    // RDKit✔️✔️:   if (args.additionalOutput->atomsPerBit) {
    // RDKit✔️✔️:     countSimulationOutput.allocateAtomsPerBit();
    // RDKit✔️✔️:   }
    if additional_output.atoms_per_bit.is_some() {
        count_simulation_output.allocate_atoms_per_bit();
    }
    // RDKit✔️✔️:   reinitAdditionalOutput(*args.additionalOutput, numAtoms);
    // RDKit✔️✔️: }
    additional_output.reset_for_atom_count(num_atoms);
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn validate_stateless_environment_generator_json(json: &str) -> Result<(), FingerprintError> {
    if json.trim().is_empty() {
        return Ok(());
    }
    let value: Value = serde_json::from_str(json)
        .map_err(|error| FingerprintError::InvalidArgumentsJson(error.to_string()))?;
    value.as_object().ok_or_else(|| {
        FingerprintError::InvalidArgumentsJson("expected JSON object".to_string())
    })?;
    Ok(())
}

fn json_value_as_bool(name: &str, value: &Value) -> Result<bool, FingerprintError> {
    if let Some(flag) = value.as_bool() {
        return Ok(flag);
    }
    if let Some(number) = value.as_u64() {
        return match number {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(FingerprintError::InvalidArgumentsJson(format!(
                "{name} must be a boolean"
            ))),
        };
    }
    if let Some(text) = value.as_str() {
        if let Ok(flag) = text.parse::<bool>() {
            return Ok(flag);
        }
        if let Ok(number) = text.parse::<u64>() {
            return match number {
                0 => Ok(false),
                1 => Ok(true),
                _ => Err(FingerprintError::InvalidArgumentsJson(format!(
                    "{name} must be a boolean"
                ))),
            };
        }
    }
    Err(FingerprintError::InvalidArgumentsJson(format!(
        "{name} must be a boolean"
    )))
}

fn json_value_as_u32(name: &str, value: &Value) -> Result<u32, FingerprintError> {
    if let Some(number) = value.as_u64() {
        return u32::try_from(number).map_err(|_| {
            FingerprintError::InvalidArgumentsJson(format!("{name} must be a 32-bit integer"))
        });
    }
    if let Some(text) = value.as_str() {
        return text.parse::<u32>().map_err(|_| {
            FingerprintError::InvalidArgumentsJson(format!("{name} must be a 32-bit integer"))
        });
    }
    Err(FingerprintError::InvalidArgumentsJson(format!(
        "{name} must be a 32-bit integer"
    )))
}

// ---------------------------------------------------------------------------
// Topological (Path-Based) Fingerprint
// RDKit source: GraphMol/Fingerprints/Fingerprints.h
// RDKit source: GraphMol/Fingerprints/FingerprintUtil.cpp
// ---------------------------------------------------------------------------

/// RDKitFP's default atom invariant and the bond-hash input preparation are
/// kept as explicit helpers so the later path/environment port can consume the
/// exact source intermediate state.
#[must_use]
pub(crate) fn rdkit_fp_atom_invariants(molecule: &Molecule) -> Vec<u32> {
    // RDKit source: FingerprintUtil.cpp lines 271-280
    // RDKit✔️✔️: void buildDefaultRDKitFingerprintAtomInvariants(
    // RDKit✔️✔️:     const ROMol &mol, std::vector<std::uint32_t> &lAtomInvariants) {
    // RDKit✔️✔️:   lAtomInvariants.clear();
    // RDKit✔️✔️:   lAtomInvariants.reserve(mol.getNumAtoms());
    // RDKit✔️✔️:   for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
    // RDKit✔️✔️:        atomIt != mol.endAtoms(); ++atomIt) {
    // RDKit✔️✔️:     unsigned int aHash = ((*atomIt)->getAtomicNum() % 128) << 1 |
    // RDKit✔️✔️:                          static_cast<unsigned int>((*atomIt)->getIsAromatic());
    // RDKit✔️✔️:     lAtomInvariants.push_back(aHash);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    molecule
        .atoms()
        .iter()
        .map(|atom| ((u32::from(atom.atomic_number()) % 128) << 1) | u32::from(atom.is_aromatic()))
        .collect()
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct RdkitFpBondHashInputs {
    pub atoms_in_path: Vec<bool>,
    pub bond_hashes: Vec<u32>,
}

pub(crate) fn rdkit_fp_generate_bond_hash_inputs(
    molecule: &Molecule,
    path: &[usize],
    use_bond_order: bool,
    atom_invariants: &[u32],
) -> Result<RdkitFpBondHashInputs, FingerprintError> {
    // RDKit source: FingerprintUtil.cpp lines 345-370
    // RDKit✔️✔️: std::vector<unsigned int> generateBondHashes(
    // RDKit✔️✔️:     const ROMol &mol, boost::dynamic_bitset<> &atomsInPath,
    // RDKit✔️✔️:     const std::vector<const Bond *> &bondCache,
    // RDKit✔️✔️:     const std::vector<short> &isQueryBond, const PATH_TYPE &path,
    // RDKit✔️✔️:     bool useBondOrder, const std::vector<std::uint32_t> *atomInvariants) {
    // RDKit✔️✔️:   PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit✔️✔️:                "bad atomInvariants size");
    // RDKit✔️✔️:   std::vector<unsigned int> bondHashes;
    // RDKit✔️✔️:   atomsInPath.reset();
    // RDKit✔️✔️:   bool queryInPath = false;
    // RDKit✔️✔️:   std::vector<unsigned int> atomDegrees(mol.getNumAtoms(), 0);
    // RDKit✔️✔️:   for (unsigned int i = 0; i < path.size() && !queryInPath; ++i) {
    // RDKit✔️✔️:     const Bond *bi = bondCache[path[i]];
    // RDKit✔️✔️:     CHECK_INVARIANT(bi, "bond not in cache");
    // RDKit✔️✔️:     atomDegrees[bi->getBeginAtomIdx()]++;
    // RDKit✔️✔️:     atomDegrees[bi->getEndAtomIdx()]++;
    // RDKit✔️✔️:     atomsInPath.set(bi->getBeginAtomIdx());
    // RDKit✔️✔️:     atomsInPath.set(bi->getEndAtomIdx());
    // RDKit✔️✔️:     if (isQueryBond[path[i]]) {
    // RDKit✔️✔️:       queryInPath = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (queryInPath) {
    // RDKit✔️✔️:     return bondHashes;
    // RDKit✔️✔️:   }
    if atom_invariants.len() < molecule.num_atoms() {
        return Err(FingerprintError::InvalidArguments {
            reason: "bad atomInvariants size",
        });
    }

    let mut atoms_in_path = vec![false; molecule.num_atoms()];
    let mut atom_degrees = vec![0_u32; molecule.num_atoms()];
    let mut query_in_path = false;
    for &bond_index in path {
        let bond = molecule
            .bonds()
            .get(bond_index)
            .ok_or(FingerprintError::InvalidArguments {
                reason: "bond not in cache",
            })?;
        let begin = bond.begin().index();
        let end = bond.end().index();
        atom_degrees[begin] += 1;
        atom_degrees[end] += 1;
        atoms_in_path[begin] = true;
        atoms_in_path[end] = true;
        if bond.query().is_some()
            || molecule.atoms()[begin].query().is_some()
            || molecule.atoms()[end].query().is_some()
        {
            query_in_path = true;
        }
    }
    if query_in_path {
        return Ok(RdkitFpBondHashInputs {
            atoms_in_path,
            bond_hashes: Vec::new(),
        });
    }

    // RDKit source: FingerprintUtil.cpp lines 372-434
    // RDKit✔️✔️:   std::vector<unsigned int> bondNbrs(path.size(), 0);
    // RDKit✔️✔️:   bondHashes.reserve(path.size() + 1);
    // RDKit✔️✔️:   for (unsigned int i = 0; i < path.size(); ++i) {
    // RDKit✔️✔️:     const Bond *bi = bondCache[path[i]];
    // RDKit✔️✔️:     for (unsigned int j = i + 1; j < path.size(); ++j) {
    // RDKit✔️✔️:       const Bond *bj = bondCache[path[j]];
    // RDKit✔️✔️:       if (bi->getBeginAtomIdx() == bj->getBeginAtomIdx() ||
    // RDKit✔️✔️:           bi->getBeginAtomIdx() == bj->getEndAtomIdx() ||
    // RDKit✔️✔️:           bi->getEndAtomIdx() == bj->getBeginAtomIdx() ||
    // RDKit✔️✔️:           bi->getEndAtomIdx() == bj->getEndAtomIdx()) {
    // RDKit✔️✔️:         ++bondNbrs[i];
    // RDKit✔️✔️:         ++bondNbrs[j];
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     unsigned int a1Hash = (*atomInvariants)[bi->getBeginAtomIdx()];
    // RDKit✔️✔️:     unsigned int a2Hash = (*atomInvariants)[bi->getEndAtomIdx()];
    // RDKit✔️✔️:     unsigned int deg1 = atomDegrees[bi->getBeginAtomIdx()];
    // RDKit✔️✔️:     unsigned int deg2 = atomDegrees[bi->getEndAtomIdx()];
    // RDKit✔️✔️:     if (a1Hash < a2Hash) {
    // RDKit✔️✔️:       std::swap(a1Hash, a2Hash);
    // RDKit✔️✔️:       std::swap(deg1, deg2);
    // RDKit✔️✔️:     } else if (a1Hash == a2Hash && deg1 < deg2) {
    // RDKit✔️✔️:       std::swap(deg1, deg2);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     unsigned int bondHash = 1;
    // RDKit✔️✔️:     if (useBondOrder) {
    // RDKit✔️✔️:       if (bi->getIsAromatic() || bi->getBondType() == Bond::AROMATIC) {
    // RDKit✔️✔️:         bondHash = Bond::AROMATIC;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         bondHash = bi->getBondType();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::uint32_t ourHash = bondNbrs[i];
    // RDKit✔️✔️:     gboost::hash_combine(ourHash, bondHash);
    // RDKit✔️✔️:     gboost::hash_combine(ourHash, a1Hash);
    // RDKit✔️✔️:     gboost::hash_combine(ourHash, deg1);
    // RDKit✔️✔️:     gboost::hash_combine(ourHash, a2Hash);
    // RDKit✔️✔️:     gboost::hash_combine(ourHash, deg2);
    // RDKit✔️✔️:     bondHashes.push_back(ourHash);
    // RDKit✔️✔️:   }
    let mut bond_neighbors = vec![0_u32; path.len()];
    for i in 0..path.len() {
        let first = &molecule.bonds()[path[i]];
        for j in (i + 1)..path.len() {
            let second = &molecule.bonds()[path[j]];
            if first.begin() == second.begin()
                || first.begin() == second.end()
                || first.end() == second.begin()
                || first.end() == second.end()
            {
                bond_neighbors[i] += 1;
                bond_neighbors[j] += 1;
            }
        }
    }
    let mut bond_hashes = Vec::with_capacity(path.len() + 1);
    for (i, &bond_index) in path.iter().enumerate() {
        let bond = &molecule.bonds()[bond_index];
        let begin = bond.begin().index();
        let end = bond.end().index();
        let mut first_hash = atom_invariants[begin];
        let mut second_hash = atom_invariants[end];
        let mut first_degree = atom_degrees[begin];
        let mut second_degree = atom_degrees[end];
        if first_hash < second_hash {
            std::mem::swap(&mut first_hash, &mut second_hash);
            std::mem::swap(&mut first_degree, &mut second_degree);
        } else if first_hash == second_hash && first_degree < second_degree {
            std::mem::swap(&mut first_degree, &mut second_degree);
        }
        let bond_hash = if use_bond_order {
            if bond.is_aromatic() || bond.order() == crate::BondOrder::Aromatic {
                rdkit_bond_type_code(crate::BondOrder::Aromatic)
            } else {
                rdkit_bond_type_code(bond.order())
            }
        } else {
            1
        };
        let mut our_hash = bond_neighbors[i];
        hash_combine(&mut our_hash, bond_hash);
        hash_combine(&mut our_hash, first_hash);
        hash_combine(&mut our_hash, first_degree);
        hash_combine(&mut our_hash, second_hash);
        hash_combine(&mut our_hash, second_degree);
        bond_hashes.push(our_hash);
    }
    Ok(RdkitFpBondHashInputs {
        atoms_in_path,
        bond_hashes,
    })
}

fn rdkit_fp_bond_neighbors(molecule: &Molecule, use_hs: bool) -> Vec<Vec<usize>> {
    let mut neighbors = vec![Vec::new(); molecule.num_bonds()];
    // RDKit source: Subgraphs.cpp lines 23-68 (`getNbrsList`).
    // RDKit✔️✔️: for (unsigned int atomIdx = 0; atomIdx < mol.getNumAtoms(); ++atomIdx) {
    // RDKit✔️✔️:   const Atom *atom = mol.getAtomWithIdx(atomIdx);
    // RDKit✔️✔️:   ROMol::OEDGE_ITER bIt1, end;
    // RDKit✔️✔️:   boost::tie(bIt1, end) = mol.getAtomBonds(atom);
    // RDKit✔️✔️:   while (bIt1 != end) {
    // RDKit✔️✔️:     const Bond *bond1 = mol[*bIt1];
    // RDKit✔️✔️:     if (useHs || bond1->getOtherAtom(atom)->getAtomicNum() != 1) {
    // RDKit✔️✔️:       int bid1 = bond1->getIdx();
    // RDKit✔️✔️:       ROMol::OEDGE_ITER bIt2 = mol.getAtomBonds(atom).first;
    // RDKit✔️✔️:       while (bIt2 != end) {
    // RDKit✔️✔️:         const Bond *bond2 = mol[*bIt2];
    // RDKit✔️✔️:         int bid2 = bond2->getIdx();
    // RDKit✔️✔️:         if (bid1 != bid2 &&
    // RDKit✔️✔️:             (useHs || bond2->getOtherAtom(atom)->getAtomicNum() != 1)) {
    // RDKit✔️✔️:           nbrs[bid1].push_back(bid2);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         ++bIt2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++bIt1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let adjacency = &molecule.topology_block().adjacency;
    for atom_index in 0..molecule.num_atoms() {
        let incident: Vec<_> = adjacency
            .neighbors_of(atom_index)
            .iter()
            .copied()
            .filter(|neighbor| {
                use_hs
                    || (molecule.atoms()[neighbor.atom_index].atomic_number() != 1
                        && molecule.atoms()[atom_index].atomic_number() != 1)
            })
            .collect();
        for first in &incident {
            for second in &incident {
                if first.bond != second.bond {
                    neighbors[first.bond.index()].push(second.bond.index());
                }
            }
        }
    }
    neighbors
}

fn rdkit_fp_recurse_subgraphs(
    neighbors: &[Vec<usize>],
    path: Vec<usize>,
    mut candidates: Vec<usize>,
    lower_len: usize,
    upper_len: usize,
    mut forbidden: Vec<bool>,
    result: &mut BTreeMap<usize, Vec<Vec<usize>>>,
) {
    // RDKit source: Subgraphs.cpp lines 120-176 (`recurseWalkRange`).
    // RDKit✔️✔️: unsigned int nsize = spath.size();
    // RDKit✔️✔️: if ((nsize >= lowerLen) && (nsize <= upperLen)) {
    // RDKit✔️✔️:   res[nsize].push_back(spath);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (nsize == upperLen || nsize > upperLen) return;
    // RDKit✔️✔️: while (cands.size() != 0) {
    // RDKit✔️✔️:   int next = cands.back(); cands.pop_back();
    // RDKit✔️✔️:   if (!forbidden[next]) {
    // RDKit✔️✔️:     forbidden[next] = 1;
    // RDKit✔️✔️:     INT_VECT tstack = cands;
    // RDKit✔️✔️:     for (int &bid : nbrs[next]) if (!forbidden[bid]) tstack.push_back(bid);
    // RDKit✔️✔️:     PATH_TYPE tpath = spath; tpath.push_back(next);
    // RDKit✔️✔️:     recurseWalkRange(nbrs, tpath, tstack, lowerLen, upperLen,
    // RDKit✔️✔️:                       forbidden, res);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let size = path.len();
    if size >= lower_len && size <= upper_len {
        result.entry(size).or_default().push(path.clone());
    }
    if size >= upper_len {
        return;
    }
    while let Some(next) = candidates.pop() {
        if forbidden.get(next).copied().unwrap_or(true) {
            continue;
        }
        forbidden[next] = true;
        let mut next_candidates = candidates.clone();
        for &candidate in &neighbors[next] {
            if !forbidden[candidate] {
                next_candidates.push(candidate);
            }
        }
        let mut next_path = path.clone();
        next_path.push(next);
        rdkit_fp_recurse_subgraphs(
            neighbors,
            next_path,
            next_candidates,
            lower_len,
            upper_len,
            forbidden.clone(),
            result,
        );
    }
}

fn extend_paths(
    adjacency: &[u8],
    dim: usize,
    paths: &[Vec<usize>],
    allow_ring_closures: i64,
    distance_matrix: Option<&[f64]>,
) -> Result<Vec<Vec<usize>>, FingerprintError> {
    // BEGIN RDKIT CPP FUNCTION Subgraphs::extendPaths
    // RDKit✔️✔️: PATH_LIST
    // RDKit✔️✔️: extendPaths(int *adjMat, unsigned int dim, const PATH_LIST &paths,
    // RDKit✔️✔️:             int allowRingClosures = -1, double *distMat = nullptr) {
    // RDKit✔️✔️:   PRECONDITION(adjMat, "no matrix");
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //  extend each of the currently active paths by adding
    // RDKit✔️✔️:   //   a single adjacent index to the end of each
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   PATH_LIST res;
    // RDKit✔️✔️:   PATH_LIST::const_iterator path;
    // RDKit✔️✔️:   for (path = paths.begin(); path != paths.end(); ++path) {
    // RDKit✔️✔️:     unsigned int endIdx = (*path)[path->size() - 1];
    // RDKit✔️✔️:     unsigned int iTab = endIdx * dim;
    // RDKit✔️✔️:     for (unsigned int otherIdx = 0; otherIdx < dim; otherIdx++) {
    // RDKit✔️✔️:       if (adjMat[iTab + otherIdx] == 1) {
    // RDKit✔️✔️:         if (distMat &&
    // RDKit✔️✔️:             distMat[path->front() * dim + otherIdx] - path->size() < -0.001) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // test 1: make sure the new atom is not already
    // RDKit✔️✔️:         //   in the path
    // RDKit✔️✔️:         auto loc =
    // RDKit✔️✔️:             std::find(path->begin(), path->end(), static_cast<int>(otherIdx));
    // RDKit✔️✔️:         // The two conditions for adding the atom are:
    // RDKit✔️✔️:         //   1) it's not there already
    // RDKit✔️✔️:         //   2) it's there, but ring closures are allowed and this
    // RDKit✔️✔️:         //      will be the last addition to the path.
    // RDKit✔️✔️:         if (loc == path->end()) {
    // RDKit✔️✔️:           // the easy case
    // RDKit✔️✔️:           // PATH_TYPE newPath=*path;
    // RDKit✔️✔️:           // newPath.push_back(otherIdx);
    // RDKit✔️✔️:           // res.push_back(newPath);
    // RDKit✔️✔️:           res.push_back(*path);
    // RDKit✔️✔️:           res.rbegin()->push_back(otherIdx);
    // RDKit✔️✔️:         } else if (allowRingClosures > 2 &&
    // RDKit✔️✔️:                    static_cast<int>(path->size()) == allowRingClosures - 1) {
    // RDKit✔️✔️:           // We *might* be adding the atom, but we need to make sure
    // RDKit✔️✔️:           // that we're not just duplicating the second to last
    // RDKit✔️✔️:           // element of the path:
    // RDKit✔️✔️:           auto rIt = path->rbegin();
    // RDKit✔️✔️:           rIt++;
    // RDKit✔️✔️:           if (*rIt != static_cast<int>(otherIdx)) {
    // RDKit✔️✔️:             // PATH_TYPE newPath=*path;
    // RDKit✔️✔️:             // newPath.push_back(otherIdx);
    // RDKit✔️✔️:             // res.push_back(newPath);
    // RDKit✔️✔️:             res.push_back(*path);
    // RDKit✔️✔️:             res.rbegin()->push_back(otherIdx);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Subgraphs::extendPaths
    let matrix_len = dim
        .checked_mul(dim)
        .ok_or(FingerprintError::InvalidArguments {
            reason: "path adjacency matrix dimensions overflow",
        })?;
    if adjacency.len() != matrix_len {
        return Err(FingerprintError::InvalidArguments {
            reason: "path adjacency matrix has invalid dimensions",
        });
    }
    if distance_matrix.is_some_and(|matrix| matrix.len() != matrix_len) {
        return Err(FingerprintError::InvalidArguments {
            reason: "path distance matrix has invalid dimensions",
        });
    }

    let mut result = Vec::new();
    for path in paths {
        let end = *path.last().ok_or(FingerprintError::InvalidArguments {
            reason: "path to extend must not be empty",
        })?;
        let start = path[0];
        if end >= dim || start >= dim {
            return Err(FingerprintError::InvalidArguments {
                reason: "path atom index is out of range",
            });
        }
        let row_offset = end * dim;
        // The source scans every column of the dense adjacency matrix in
        // atom-index order. Keep that order rather than the molecule's bond
        // insertion order.
        for other_idx in 0..dim {
            if adjacency[row_offset + other_idx] != 1 {
                continue;
            }
            if distance_matrix.is_some_and(|matrix| {
                matrix[start * dim + other_idx] - (path.len() as f64) < -0.001
            }) {
                continue;
            }
            if !path.contains(&other_idx) {
                let mut next = path.clone();
                next.push(other_idx);
                result.push(next);
            } else if allow_ring_closures > 2
                && i64::try_from(path.len()).ok() == Some(allow_ring_closures - 1)
                && path[path.len() - 2] != other_idx
            {
                let mut next = path.clone();
                next.push(other_idx);
                result.push(next);
            }
        }
    }
    Ok(result)
}

fn path_finder_helper(
    adjacency: &[u8],
    dim: usize,
    min_len: usize,
    max_len: usize,
    rooted_at_atom: i64,
    distance_matrix: Option<&[f64]>,
) -> Result<BTreeMap<usize, Vec<Vec<usize>>>, FingerprintError> {
    // BEGIN RDKIT CPP FUNCTION Subgraphs::pathFinderHelper
    // RDKit✔️✔️: INT_PATH_LIST_MAP
    // RDKit✔️✔️: pathFinderHelper(int *adjMat, unsigned int dim, unsigned int minLen,
    // RDKit✔️✔️:                  unsigned int maxLen, int rootedAtAtom, double *distMat) {
    // RDKit✔️✔️:   PRECONDITION(adjMat, "no matrix");
    // RDKit✔️✔️:   PRECONDITION(minLen <= maxLen, "bad lengths provided");
    // RDKit✔️✔️:   // finds all paths of length N using an adjacency matrix,
    // RDKit✔️✔️:   //  which is constructed elsewhere
    // RDKit✔️✔️:   INT_PATH_LIST_MAP res;
    // RDKit✔️✔️:   PATH_LIST paths;
    // RDKit✔️✔️:   paths.clear();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (rootedAtAtom < 0) {
    // RDKit✔️✔️:     // start a path at each possible index
    // RDKit✔️✔️:     for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:       PATH_TYPE tPath;
    // RDKit✔️✔️:       tPath.push_back(i);
    // RDKit✔️✔️:       paths.push_back(tPath);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (rootedAtAtom < static_cast<int>(dim)) {
    // RDKit✔️✔️:     // only start a path at the atom of interest:
    // RDKit✔️✔️:     PATH_TYPE tPath;
    // RDKit✔️✔️:     tPath.push_back(rootedAtAtom);
    // RDKit✔️✔️:     paths.push_back(tPath);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // and build them up one index at a time:
    // RDKit✔️✔️:   for (unsigned int length = 1; length < maxLen; length++) {
    // RDKit✔️✔️:     // extend each path:
    // RDKit✔️✔️:     if (length >= minLen) {
    // RDKit✔️✔️:       res[length] = paths;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     paths = extendPaths(adjMat, dim, paths, maxLen, distMat);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res[maxLen] = paths;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Subgraphs::pathFinderHelper
    if min_len > max_len {
        return Err(FingerprintError::InvalidArguments {
            reason: "minimum path length exceeds maximum path length",
        });
    }
    let allow_ring_closures =
        i64::try_from(max_len).map_err(|_| FingerprintError::InvalidArguments {
            reason: "path length exceeds supported range",
        })?;

    let mut result = BTreeMap::new();
    let mut paths = Vec::new();
    if rooted_at_atom < 0 {
        for atom_index in 0..dim {
            paths.push(vec![atom_index]);
        }
    } else if usize::try_from(rooted_at_atom).is_ok_and(|root| root < dim) {
        paths.push(vec![rooted_at_atom as usize]);
    } else {
        return Ok(result);
    }

    for length in 1..max_len {
        if length >= min_len {
            result.insert(length, paths.clone());
        }
        paths = extend_paths(adjacency, dim, &paths, allow_ring_closures, distance_matrix)?;
    }
    result.insert(max_len, paths);
    Ok(result)
}

fn rdkit_fp_bond_between_atoms(molecule: &Molecule, begin: usize, end: usize) -> Option<usize> {
    // RDKit source: GraphMol/ROMol.cpp lines 338-350.
    // RDKit✔️✔️: const Bond *ROMol::getBondBetweenAtoms(unsigned int idx1,
    // RDKit✔️✔️:                                        unsigned int idx2) const {
    // RDKit✔️✔️:   URANGE_CHECK(idx1, getNumAtoms());
    // RDKit✔️✔️:   URANGE_CHECK(idx2, getNumAtoms());
    // RDKit✔️✔️:   const Bond *res = nullptr;
    // RDKit✔️✔️:   auto [edge, found] = boost::edge(boost::vertex(idx1, d_graph),
    // RDKit✔️✔️:                                    boost::vertex(idx2, d_graph), d_graph);
    // RDKit✔️✔️:   if (found) {
    // RDKit✔️✔️:     res = d_graph[edge];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    if begin >= molecule.num_atoms() || end >= molecule.num_atoms() {
        return None;
    }
    molecule
        .topology_block()
        .adjacency
        .neighbors_of(begin)
        .iter()
        .find(|neighbor| neighbor.atom_index == end)
        .map(|neighbor| neighbor.bond.index())
}

fn find_all_paths_of_lengths_m_to_n(
    molecule: &Molecule,
    mut lower_len: usize,
    mut upper_len: usize,
    use_bonds: bool,
    use_hs: bool,
    rooted_at_atom: i64,
    only_shortest_paths: bool,
) -> Result<BTreeMap<usize, Vec<Vec<usize>>>, FingerprintError> {
    // BEGIN RDKIT CPP FUNCTION findAllPathsOfLengthsMtoN
    // RDKit✔️✔️: INT_PATH_LIST_MAP
    // RDKit✔️✔️: findAllPathsOfLengthsMtoN(const ROMol &mol, unsigned int lowerLen,
    // RDKit✔️✔️:                           unsigned int upperLen, bool useBonds, bool useHs,
    // RDKit✔️✔️:                           int rootedAtAtom, bool onlyShortestPaths) {
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //  We can't be clever here and just use the bond adjacency matrix
    // RDKit✔️✔️:   //  to solve this problem when useBonds is true.  This is because
    // RDKit✔️✔️:   //  the bond adjacency matrices for the molecules C1CC1 and CC(C)C
    // RDKit✔️✔️:   //  are indistinguishable.  In the second case, t-butane (and
    // RDKit✔️✔️:   //  anything else with a T junction), we'll get some subgraphs mixed
    // RDKit✔️✔️:   //  in with the paths.  So we have to construct paths of atoms and
    // RDKit✔️✔️:   //  then convert them into bond paths.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   PRECONDITION(lowerLen <= upperLen, "");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // the molecule owns the distance matrix pointer (if we need to get it)
    // RDKit✔️✔️:   double *distMat = onlyShortestPaths ? MolOps::getDistanceMat(mol) : nullptr;
    // RDKit✔️✔️:   int *adjMat, dim;
    // RDKit✔️✔️:   dim = mol.getNumAtoms();
    // RDKit✔️✔️:   adjMat = new int[dim * dim];
    // RDKit✔️✔️:   memset((void *)adjMat, 0, dim * dim * sizeof(int));
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!distMat) {
    // RDKit✔️✔️:     // generate the adjacency matrix by hand by looping over the bonds
    // RDKit✔️✔️:     ROMol::ConstBondIterator bondIt;
    // RDKit✔️✔️:     for (bondIt = mol.beginBonds(); bondIt != mol.endBonds(); bondIt++) {
    // RDKit✔️✔️:       Atom *beg = (*bondIt)->getBeginAtom();
    // RDKit✔️✔️:       Atom *end = (*bondIt)->getEndAtom();
    // RDKit✔️✔️:       // check for H, which we might be skipping
    // RDKit✔️✔️:       if (useHs || (beg->getAtomicNum() != 1 && end->getAtomicNum() != 1)) {
    // RDKit✔️✔️:         adjMat[beg->getIdx() * dim + end->getIdx()] = 1;
    // RDKit✔️✔️:         adjMat[end->getIdx() * dim + beg->getIdx()] = 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     // if we have the distance matrix, we can just loop over that:
    // RDKit✔️✔️:     for (auto i = 0; i < dim; ++i) {
    // RDKit✔️✔️:       for (auto j = i + 1; j < dim; ++j) {
    // RDKit✔️✔️:         if (fabs(distMat[i * dim + j] - 1) < 1e-4) {
    // RDKit✔️✔️:           adjMat[i * dim + j] = 1;
    // RDKit✔️✔️:           adjMat[j * dim + i] = 1;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // if we're using bonds, we'll need to find paths of length N+1,
    // RDKit✔️✔️:   // then convert them
    // RDKit✔️✔️:   if (useBonds) {
    // RDKit✔️✔️:     ++lowerLen;
    // RDKit✔️✔️:     ++upperLen;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // find the paths themselves
    // RDKit✔️✔️:   INT_PATH_LIST_MAP atomPaths = Subgraphs::pathFinderHelper(
    // RDKit✔️✔️:       adjMat, dim, lowerLen, upperLen, rootedAtAtom, distMat);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // clean up the adjacency matrix
    // RDKit✔️✔️:   delete[] adjMat;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   INT_PATH_LIST_MAP res;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //--------------------------------------------------------
    // RDKit✔️✔️:   // loop through all the paths we have and make sure that there are
    // RDKit✔️✔️:   // no duplicates (duplicate = contains identical bond indices)
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //  We need to use the bond paths for this duplicate finding
    // RDKit✔️✔️:   //  because, in rings, there can be many paths which share atom
    // RDKit✔️✔️:   //  indices but which have different bond compositions. For example,
    // RDKit✔️✔️:   //  there is only one "atom unique" path of length 5 bonds (6 atoms)
    // RDKit✔️✔️:   //  through a 6-ring, but there are six bond paths.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   if (!useBonds && lowerLen >= 1) {
    // RDKit✔️✔️:     res[1] = atomPaths[1];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (useBonds || upperLen > 1) {
    // RDKit✔️✔️:     for (unsigned int i = lowerLen; i <= upperLen; ++i) {
    // RDKit✔️✔️:       if (i <= 1) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       std::vector<boost::dynamic_bitset<>> invars;
    // RDKit✔️✔️:
    // RDKit✔️✔️:       for (PATH_LIST::const_iterator vivI = atomPaths[i].begin();
    // RDKit✔️✔️:            vivI != atomPaths[i].end(); ++vivI) {
    // RDKit✔️✔️:         boost::dynamic_bitset<> invar(mol.getNumBonds());
    // RDKit✔️✔️:         const PATH_TYPE &resi = *vivI;
    // RDKit✔️✔️:         PATH_TYPE locV;
    // RDKit✔️✔️:         locV.reserve(i);
    // RDKit✔️✔️:         for (unsigned int j = 0; j < i - 1; j++) {
    // RDKit✔️✔️:           const Bond *bond = mol.getBondBetweenAtoms(resi[j], resi[j + 1]);
    // RDKit✔️✔️:           locV.push_back(bond->getIdx());
    // RDKit✔️✔️:           invar.set(bond->getIdx());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (std::find(invars.begin(), invars.end(), invar) == invars.end()) {
    // RDKit✔️✔️:           invars.push_back(invar);
    // RDKit✔️✔️:           if (useBonds) {
    // RDKit✔️✔️:             res[i - 1].push_back(locV);
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             res[i].push_back(resi);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION findAllPathsOfLengthsMtoN
    if lower_len > upper_len {
        return Err(FingerprintError::InvalidArguments {
            reason: "minimum path length exceeds maximum path length",
        });
    }

    let distance_matrix = only_shortest_paths
        .then(|| crate::chemistry::matrices::topological_distance_matrix(molecule));
    let dim = molecule.num_atoms();
    let matrix_len = dim
        .checked_mul(dim)
        .ok_or(FingerprintError::InvalidArguments {
            reason: "path adjacency matrix dimensions overflow",
        })?;
    let mut adjacency = vec![0u8; matrix_len];
    if let Some(matrix) = distance_matrix.as_deref() {
        for first in 0..dim {
            for second in (first + 1)..dim {
                if (matrix[first * dim + second] - 1.0).abs() < 1.0e-4 {
                    adjacency[first * dim + second] = 1;
                    adjacency[second * dim + first] = 1;
                }
            }
        }
    } else {
        for bond in molecule.bonds() {
            let begin = bond.begin().index();
            let end = bond.end().index();
            if use_hs
                || (molecule.atoms()[begin].atomic_number() != 1
                    && molecule.atoms()[end].atomic_number() != 1)
            {
                adjacency[begin * dim + end] = 1;
                adjacency[end * dim + begin] = 1;
            }
        }
    }

    if use_bonds {
        lower_len = lower_len
            .checked_add(1)
            .ok_or(FingerprintError::InvalidArguments {
                reason: "minimum bond path length exceeds supported range",
            })?;
        upper_len = upper_len
            .checked_add(1)
            .ok_or(FingerprintError::InvalidArguments {
                reason: "maximum bond path length exceeds supported range",
            })?;
    }
    let atom_paths = path_finder_helper(
        &adjacency,
        dim,
        lower_len,
        upper_len,
        rooted_at_atom,
        distance_matrix.as_deref(),
    )?;

    let mut result = BTreeMap::new();
    if !use_bonds && lower_len >= 1 {
        result.insert(1, atom_paths.get(&1).cloned().unwrap_or_default());
    }
    if use_bonds || upper_len > 1 {
        for path_length in lower_len..=upper_len {
            if path_length <= 1 {
                continue;
            }
            let mut invariants: Vec<Vec<bool>> = Vec::new();
            if let Some(paths) = atom_paths.get(&path_length) {
                for atom_path in paths {
                    if atom_path.len() < path_length {
                        return Err(FingerprintError::InvalidArguments {
                            reason: "enumerated atom path is shorter than its length key",
                        });
                    }
                    let mut invariant = vec![false; molecule.num_bonds()];
                    let mut bond_path = Vec::with_capacity(path_length);
                    for atom_pair in atom_path[..path_length].windows(2) {
                        let bond_index =
                            rdkit_fp_bond_between_atoms(molecule, atom_pair[0], atom_pair[1])
                                .ok_or(FingerprintError::InvalidArguments {
                                    reason: "path contains no connecting bond",
                                })?;
                        bond_path.push(bond_index);
                        invariant[bond_index] = true;
                    }
                    if !invariants.contains(&invariant) {
                        invariants.push(invariant);
                        if use_bonds {
                            result.entry(path_length - 1).or_default().push(bond_path);
                        } else {
                            result
                                .entry(path_length)
                                .or_default()
                                .push(atom_path.clone());
                        }
                    }
                }
            }
        }
    }
    Ok(result)
}

fn find_all_paths_of_length_n(
    molecule: &Molecule,
    target_len: usize,
    use_bonds: bool,
    use_hs: bool,
    rooted_at_atom: i64,
    only_shortest_paths: bool,
) -> Result<Vec<Vec<usize>>, FingerprintError> {
    // BEGIN RDKIT CPP FUNCTION findAllPathsOfLengthN
    // RDKit✔️✔️: PATH_LIST
    // RDKit✔️✔️: findAllPathsOfLengthN(const ROMol &mol, unsigned int targetLen, bool useBonds,
    // RDKit✔️✔️:                       bool useHs, int rootedAtAtom, bool onlyShortestPaths) {
    // RDKit✔️✔️:   return findAllPathsOfLengthsMtoN(mol, targetLen, targetLen, useBonds, useHs,
    // RDKit✔️✔️:                                    rootedAtAtom, onlyShortestPaths)[targetLen];
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION findAllPathsOfLengthN
    let mut paths = find_all_paths_of_lengths_m_to_n(
        molecule,
        target_len,
        target_len,
        use_bonds,
        use_hs,
        rooted_at_atom,
        only_shortest_paths,
    )?;
    Ok(paths.remove(&target_len).unwrap_or_default())
}

fn enumerate_fingerprint_paths_for_root(
    molecule: &Molecule,
    lower: usize,
    upper: usize,
    use_hs: bool,
    branched_paths: bool,
    root: Option<u32>,
) -> Result<BTreeMap<usize, Vec<Vec<usize>>>, FingerprintError> {
    if branched_paths {
        let neighbors = rdkit_fp_bond_neighbors(molecule, use_hs);
        let mut result = (lower..=upper)
            .map(|length| (length, Vec::new()))
            .collect::<BTreeMap<_, _>>();
        let mut forbidden = vec![false; molecule.num_bonds()];
        for bond_index in 0..molecule.num_bonds() {
            let bond = &molecule.bonds()[bond_index];
            if !neighbors[bond_index].is_empty()
                || (use_hs
                    || (molecule.atoms()[bond.begin().index()].atomic_number() != 1
                        && molecule.atoms()[bond.end().index()].atomic_number() != 1))
            {
                if let Some(root) = root.map(|value| value as usize) {
                    if root >= molecule.num_atoms()
                        || (root != bond.begin().index() && root != bond.end().index())
                    {
                        continue;
                    }
                }
                if forbidden[bond_index] {
                    continue;
                }
                forbidden[bond_index] = true;
                rdkit_fp_recurse_subgraphs(
                    &neighbors,
                    vec![bond_index],
                    neighbors[bond_index].clone(),
                    lower,
                    upper,
                    forbidden.clone(),
                    &mut result,
                );
            }
        }
        return Ok(result);
    }

    find_all_paths_of_lengths_m_to_n(
        molecule,
        lower,
        upper,
        true,
        use_hs,
        root.map_or(-1, i64::from),
        false,
    )
}

pub(crate) fn enumerate_fingerprint_paths(
    molecule: &Molecule,
    min_path: u32,
    max_path: u32,
    use_hs: bool,
    branched_paths: bool,
    from_atoms: Option<&[u32]>,
) -> Result<BTreeMap<usize, Vec<Vec<usize>>>, FingerprintError> {
    // BEGIN RDKIT CPP FUNCTION RDKitFPUtils::enumerateAllPaths
    // RDKit✔️✔️: void enumerateAllPaths(const ROMol &mol, INT_PATH_LIST_MAP &allPaths,
    // RDKit✔️✔️:                        const std::vector<std::uint32_t> *fromAtoms,
    // RDKit✔️✔️:                        bool branchedPaths, bool useHs, unsigned int minPath,
    // RDKit✔️✔️:                        unsigned int maxPath) {
    // RDKit✔️✔️:   if (!fromAtoms) {
    // RDKit✔️✔️:     if (branchedPaths) {
    // RDKit✔️✔️:       allPaths = findAllSubgraphsOfLengthsMtoN(mol, minPath, maxPath, useHs);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       allPaths = findAllPathsOfLengthsMtoN(mol, minPath, maxPath, true, useHs);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     for (auto aidx : *fromAtoms) {
    // RDKit✔️✔️:       INT_PATH_LIST_MAP tPaths;
    // RDKit✔️✔️:       if (branchedPaths) {
    // RDKit✔️✔️:         tPaths =
    // RDKit✔️✔️:             findAllSubgraphsOfLengthsMtoN(mol, minPath, maxPath, useHs, aidx);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         tPaths =
    // RDKit✔️✔️:             findAllPathsOfLengthsMtoN(mol, minPath, maxPath, true, useHs, aidx);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       for (INT_PATH_LIST_MAP::const_iterator tpit = tPaths.begin();
    // RDKit✔️✔️:            tpit != tPaths.end(); ++tpit) {
    // RDKit✔️✔️:         allPaths[tpit->first].insert(allPaths[tpit->first].begin(),
    // RDKit✔️✔️:                                          tpit->second.begin(), tpit->second.end());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKitFPUtils::enumerateAllPaths
    // Local complexity review: each source call and Rust helper enumerates the
    // same path/subgraph state once per requested root. Both prepend each root
    // group to the per-length vector; no molecule or completed path map is
    // cloned beyond the source-equivalent recursive path state.
    if min_path == 0 || max_path < min_path {
        return Err(FingerprintError::InvalidArguments {
            reason: "invalid path lengths",
        });
    }
    let lower = min_path as usize;
    let upper = max_path as usize;
    let Some(roots) = from_atoms else {
        return enumerate_fingerprint_paths_for_root(
            molecule,
            lower,
            upper,
            use_hs,
            branched_paths,
            None,
        );
    };

    let mut result: BTreeMap<usize, Vec<Vec<usize>>> = BTreeMap::new();
    for &root in roots {
        let rooted_paths = enumerate_fingerprint_paths_for_root(
            molecule,
            lower,
            upper,
            use_hs,
            branched_paths,
            Some(root),
        )?;
        for (length, paths) in rooted_paths {
            result.entry(length).or_default().splice(0..0, paths);
        }
    }
    Ok(result)
}

#[derive(Debug)]
struct LayeredFingerprintPreparation<'a> {
    ring_info: Cow<'a, crate::RingInfo>,
    bond_cache: Vec<&'a crate::Bond>,
    query_masks: Vec<u8>,
    aromatic_atoms: Vec<bool>,
    atomic_numbers: Vec<u32>,
}

fn prepare_layered_fingerprint<'a>(
    molecule: &'a Molecule,
    min_path: u32,
    max_path: u32,
    fp_size: usize,
    atom_counts: Option<&[u32]>,
    set_only_bits: Option<&Fingerprint>,
) -> Result<LayeredFingerprintPreparation<'a>, FingerprintError> {
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
    // Local complexity review: validation is O(1), ring preparation reuses an
    // initialized cache or performs the same exact-SSSR calculation, and the
    // bond/atom caches are each filled by one ordered O(B)/O(A) pass. The Rust
    // representation performs the same allocations and no molecule clone.
    if min_path == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "minPath==0",
        });
    }
    if max_path < min_path {
        return Err(FingerprintError::InvalidArguments {
            reason: "maxPath<minPath",
        });
    }
    if fp_size == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "fpSize==0",
        });
    }
    if atom_counts.is_some_and(|counts| counts.len() < molecule.num_atoms()) {
        return Err(FingerprintError::InvalidArguments {
            reason: "bad atomCounts size",
        });
    }
    if set_only_bits.is_some_and(|bits| bits.n_bits() != fp_size) {
        return Err(FingerprintError::InvalidArguments {
            reason: "bad setOnlyBits size",
        });
    }

    let ring_info = match molecule
        .derived_cache()
        .rings
        .as_ref()
        .filter(|rings| rings.is_initialized())
    {
        Some(rings) => Cow::Borrowed(rings),
        None => Cow::Owned(crate::rings::find_sssr(molecule)?),
    };

    let mut bond_cache = Vec::with_capacity(molecule.num_bonds());
    let mut query_masks = Vec::with_capacity(molecule.num_bonds());
    for bond in molecule.bonds() {
        let mut query_mask = 0u8;
        if crate::search::query::is_complex_bond_query(bond) {
            query_mask = 0x1;
        }
        if crate::search::query::is_complex_atom_query(&molecule.atoms()[bond.begin().index()]) {
            query_mask |= 0x2;
        }
        if crate::search::query::is_complex_atom_query(&molecule.atoms()[bond.end().index()]) {
            query_mask |= 0x4;
        }
        bond_cache.push(bond);
        query_masks.push(query_mask);
    }

    let mut aromatic_atoms = Vec::with_capacity(molecule.num_atoms());
    let mut atomic_numbers = Vec::with_capacity(molecule.num_atoms());
    for atom in molecule.atoms() {
        aromatic_atoms.push(crate::search::query::is_atom_aromatic(atom, molecule));
        atomic_numbers.push(u32::from(atom.atomic_number()));
    }

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
    bond: &crate::Bond,
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
        rdkit_bond_type_code(bond.order())
    } else {
        rdkit_bond_type_code(BondOrder::Single)
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
fn layered_ring_presence_hash(
    bond: &crate::Bond,
    ring_info: &crate::RingInfo,
    path_queries: u8,
) -> Option<u32> {
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
    (crate::search::query::query_is_bond_in_ring(bond, ring_info) != 0).then_some(1)
}

#[inline]
fn layered_min_ring_size_hash(
    bond: &crate::Bond,
    ring_info: &crate::RingInfo,
    path_queries: u8,
) -> Option<u32> {
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
    Some((crate::search::query::query_bond_min_ring_size(bond, ring_info) % 8) as u32)
}

fn project_layered_path(
    hash_layers: &mut [Vec<u32>],
    atoms_in_path: &[bool],
    fp_size: usize,
    set_only_bits: Option<&Fingerprint>,
    result_words: &mut [u64],
    mut atom_counts: Option<&mut [u32]>,
) {
    // BEGIN RDKIT CPP BLOCK LayeredFingerprintMol projection
    // RDKit✔️✔️:       unsigned int l = 0;
    // RDKit✔️✔️:       bool flaggedPath = false;
    // RDKit✔️✔️:       for (auto layerIt = hashLayers.begin(); layerIt != hashLayers.end();
    // RDKit✔️✔️:            ++layerIt, ++l) {
    // RDKit✔️✔️:         if (!layerIt->size()) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // ----
    // RDKit✔️✔️:         std::sort(layerIt->begin(), layerIt->end());
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // finally, we will add the number of distinct atoms in the path at the
    // RDKit✔️✔️:         // end
    // RDKit✔️✔️:         // of the vect. This allows us to distinguish C1CC1 from CC(C)C
    // RDKit✔️✔️:         layerIt->push_back(static_cast<unsigned int>(atomsInPath.count()));
    // RDKit✔️✔️:
    // RDKit✔️✔️:         layerIt->push_back(l + 1);
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // hash the path to generate a seed:
    // RDKit✔️✔️:         unsigned long seed =
    // RDKit✔️✔️:             gboost::hash_range(layerIt->begin(), layerIt->end());
    // RDKit✔️✔️:
    // RDKit✔️✔️:         unsigned int bitId = seed % fpSize;
    // RDKit✔️✔️:         if (!setOnlyBits || (*setOnlyBits)[bitId]) {
    // RDKit✔️✔️:           res->setBit(bitId);
    // RDKit✔️✔️:           if (atomCounts && !flaggedPath) {
    // RDKit✔️✔️:             for (unsigned int aIdx = 0; aIdx < atomsInPath.size(); ++aIdx) {
    // RDKit✔️✔️:               if (atomsInPath[aIdx]) {
    // RDKit✔️✔️:                 (*atomCounts)[aIdx] += 1;
    // RDKit✔️✔️:               }
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             flaggedPath = true;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // END RDKIT CPP BLOCK LayeredFingerprintMol projection
    // Local complexity review: both implementations sort each nonempty layer,
    // append two scalar suffixes, hash it once, set one indexed bit, and scan
    // the atom mask at most once per accepted path. Rust uses the same dense
    // word representation and performs no additional completed-path clone.
    debug_assert!(fp_size > 0);
    debug_assert_eq!(result_words.len(), fp_size.div_ceil(64));
    debug_assert!(set_only_bits.is_none_or(|bits| bits.n_bits() == fp_size));
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
        let accepted =
            set_only_bits.is_none_or(|mask| mask.bits[bit_id / 64] & (1u64 << (bit_id % 64)) != 0);
        if !accepted {
            continue;
        }
        result_words[bit_id / 64] |= 1u64 << (bit_id % 64);
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
}

pub fn layered_fingerprint(
    molecule: &Molecule,
    params: &LayeredFingerprintParams,
) -> Result<Fingerprint, FingerprintError> {
    Ok(layered_fingerprint_with_output(molecule, params)?.fingerprint)
}

/// Compute a source-backed Layered fingerprint and optional atom counts.
///
/// This read-only operation neither mutates the molecule nor stores ring or
/// fingerprint intermediates in it. Calls can be repeated, interleaved with
/// other fingerprint families, or run concurrently on shared molecules.
/// Invalid path bounds, widths, count lengths, mask widths, and roots return a
/// structured [`FingerprintError`].
pub fn layered_fingerprint_with_output(
    molecule: &Molecule,
    params: &LayeredFingerprintParams,
) -> Result<LayeredFingerprintResult, FingerprintError> {
    params.validate()?;
    if params.from_atoms.as_ref().is_some_and(|roots| {
        roots
            .iter()
            .any(|&root| root as usize >= molecule.num_atoms())
    }) {
        return Err(FingerprintError::InvalidArguments {
            reason: "fromAtoms contains atom index out of range",
        });
    }

    let mut atom_counts = params.atom_counts.clone();
    let prepared = prepare_layered_fingerprint(
        molecule,
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
    // RDKit❗✔️:       allPaths = findAllPathsOfLengthsMtoN(mol, minPath, maxPath, false);
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
    // RDKit✔️✔️:                                          tpit->second.begin(), tpit->second.end());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP BLOCK LayeredFingerprintMol path selection
    // The pinned unrooted linear call passes `false` in the `useBonds` slot
    // and can index atom paths as bonds. Isolated pinned-RDKit calls terminate
    // with SIGSEGV for common acyclic molecules. The maintained Layered
    // contract uses the header-documented bond-path semantics shared by the
    // rooted call; COSMolKit does not reproduce an upstream process crash.
    let all_paths = enumerate_fingerprint_paths(
        molecule,
        params.min_path,
        params.max_path,
        false,
        params.branched_paths,
        params.from_atoms.as_deref(),
    )?;
    let fp_size = params.fp_size as usize;
    let mut result_words = vec![0u64; fp_size.div_ceil(64)];

    for paths in all_paths.values() {
        for path in paths {
            // BEGIN RDKIT CPP BLOCK LayeredFingerprintMol path intermediates
            // RDKit✔️✔️:       std::vector<std::vector<unsigned int>> hashLayers(maxFingerprintLayers);
            // RDKit✔️✔️:       for (unsigned int i = 0; i < maxFingerprintLayers; ++i) {
            // RDKit✔️✔️:         if (layerFlags & (0x1 << i)) {
            // RDKit✔️✔️:           hashLayers[i].reserve(maxPath);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:
            // RDKit✔️✔️:       // details about what kinds of query features appear on the path:
            // RDKit✔️✔️:       unsigned int pathQueries = 0;
            // RDKit✔️✔️:       for (int pIt : path) {
            // RDKit✔️✔️:         pathQueries |= isQueryBond[pIt];
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:
            // RDKit✔️✔️:       // calculate the number of neighbors each bond has in the path:
            // RDKit✔️✔️:       std::vector<unsigned int> bondNbrs(path.size(), 0);
            // RDKit✔️✔️:       atomsInPath.reset();
            // RDKit✔️✔️:
            // RDKit✔️✔️:       std::vector<unsigned int> atomDegrees(mol.getNumAtoms(), 0);
            // RDKit✔️✔️:       for (int i : path) {
            // RDKit✔️✔️:         const Bond *bi = bondCache[i];
            // RDKit✔️✔️:         atomDegrees[bi->getBeginAtomIdx()]++;
            // RDKit✔️✔️:         atomDegrees[bi->getEndAtomIdx()]++;
            // RDKit✔️✔️:         atomsInPath.set(bi->getBeginAtomIdx());
            // RDKit✔️✔️:         atomsInPath.set(bi->getEndAtomIdx());
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:
            // RDKit✔️✔️:       for (unsigned int i = 0; i < path.size(); ++i) {
            // RDKit✔️✔️:         const Bond *bi = bondCache[path[i]];
            // RDKit✔️✔️:         for (unsigned int j = i + 1; j < path.size(); ++j) {
            // RDKit✔️✔️:           const Bond *bj = bondCache[path[j]];
            // RDKit✔️✔️:           if (bi->getBeginAtomIdx() == bj->getBeginAtomIdx() ||
            // RDKit✔️✔️:               bi->getBeginAtomIdx() == bj->getEndAtomIdx() ||
            // RDKit✔️✔️:               bi->getEndAtomIdx() == bj->getBeginAtomIdx() ||
            // RDKit✔️✔️:               bi->getEndAtomIdx() == bj->getEndAtomIdx()) {
            // RDKit✔️✔️:             ++bondNbrs[i];
            // RDKit✔️✔️:             ++bondNbrs[j];
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:         }
            // END RDKIT CPP BLOCK LayeredFingerprintMol path intermediates
            // Local complexity review: source and Rust both allocate ten
            // layer vectors, scan P bonds for query/degree state, perform the
            // same O(P^2) pairwise adjacency count, and encode each bond once.
            let mut hash_layers = vec![Vec::new(); LAYERED_FINGERPRINT_MAX_LAYERS];
            for (layer_index, layer) in hash_layers.iter_mut().enumerate() {
                if params.layers.bits() & (1u32 << layer_index) != 0 {
                    layer.reserve(params.max_path as usize);
                }
            }

            let mut path_queries = 0u8;
            let mut atom_degrees = vec![0u32; molecule.num_atoms()];
            let mut atoms_in_path = vec![false; molecule.num_atoms()];
            for &bond_index in path {
                let bond = *prepared.bond_cache.get(bond_index).ok_or(
                    FingerprintError::InvalidArguments {
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
                &mut result_words,
                atom_counts.as_deref_mut(),
            );
        }
    }

    Ok(LayeredFingerprintResult {
        fingerprint: Fingerprint {
            bits: result_words,
            n_bits: fp_size,
        },
        atom_counts,
    })
}

/// One source `RDKitFPAtomEnv` produced by the RDKitFP environment generator.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct RdkitFpEnvironment {
    bit_id: u64,
    atoms_in_path: Vec<bool>,
    bond_path: Vec<usize>,
}

impl RdkitFpEnvironment {
    #[must_use]
    fn bit_id(&self) -> u64 {
        self.bit_id
    }

    /// RDKit source: RDKitFPGenerator.cpp lines 121-149
    // RDKit✔️✔️: template <typename OutputType>
    // RDKit✔️✔️: void RDKitFPAtomEnv<OutputType>::updateAdditionalOutput(
    // RDKit✔️✔️:     AdditionalOutput *additionalOutput, size_t bitId) const {
    // RDKit✔️✔️:   PRECONDITION(additionalOutput, "bad output pointer");
    // RDKit✔️✔️:   if (additionalOutput->bitPaths) {
    // RDKit✔️✔️:     (*additionalOutput->bitPaths)[bitId].push_back(d_bondPath);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (additionalOutput->atomToBits || additionalOutput->atomCounts ||
    // RDKit✔️✔️:       additionalOutput->atomsPerBit) {
    // RDKit✔️✔️:     if (additionalOutput->atomsPerBit) {
    // RDKit✔️✔️:       (*additionalOutput->atomsPerBit)[bitId].emplace_back();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (size_t i = 0; i < d_atomsInPath.size(); ++i) {
    // RDKit✔️✔️:       if (d_atomsInPath[i]) {
    // RDKit✔️✔️:         if (additionalOutput->atomsPerBit) {
    // RDKit✔️✔️:           (*additionalOutput->atomsPerBit)[bitId].back().push_back(i);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (additionalOutput->atomToBits) {
    // RDKit✔️✔️:           auto &alist = additionalOutput->atomToBits->at(i);
    // RDKit✔️✔️:           if (std::find(alist.begin(), alist.end(), bitId) == alist.end()) {
    // RDKit✔️✔️:             alist.push_back(bitId);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (additionalOutput->atomCounts) {
    // RDKit✔️✔️:           additionalOutput->atomCounts->at(i)++;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    fn update_additional_output(&self, additional_output: &mut AdditionalOutput, bit_id: u64) {
        if let Some(bit_paths) = additional_output.bit_paths.as_mut() {
            bit_paths
                .entry(bit_id)
                .or_default()
                .push(self.bond_path.clone());
        }
        if additional_output.atom_to_bits.is_some()
            || additional_output.atom_counts.is_some()
            || additional_output.atoms_per_bit.is_some()
        {
            if let Some(atoms_per_bit) = additional_output.atoms_per_bit.as_mut() {
                atoms_per_bit.entry(bit_id).or_default().push(Vec::new());
            }
            for (atom_index, &in_path) in self.atoms_in_path.iter().enumerate() {
                if !in_path {
                    continue;
                }
                if let Some(atoms_per_bit) = additional_output.atoms_per_bit.as_mut() {
                    atoms_per_bit
                        .get_mut(&bit_id)
                        .expect("atomsPerBit entry allocated above")
                        .last_mut()
                        .expect("atomsPerBit path allocated above")
                        .push(atom_index);
                }
                if let Some(atom_to_bits) = additional_output.atom_to_bits.as_mut() {
                    let bits = &mut atom_to_bits[atom_index];
                    if !bits.contains(&bit_id) {
                        bits.push(bit_id);
                    }
                }
                if let Some(atom_counts) = additional_output.atom_counts.as_mut() {
                    atom_counts[atom_index] += 1;
                }
            }
        }
    }
}

pub(crate) fn generate_rdkit_fp_environments(
    molecule: &Molecule,
    params: &TopologicalFingerprintParams,
    atom_invariants: &[u32],
) -> Result<Vec<RdkitFpEnvironment>, FingerprintError> {
    if atom_invariants.len() < molecule.num_atoms() {
        return Err(FingerprintError::InvalidArguments {
            reason: "bad atomInvariants size",
        });
    }

    // RDKit source: RDKitFPGenerator.cpp lines 182-236 (`getEnvironments`).
    // RDKit✔️✔️: PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit✔️✔️:              "bad atomInvariants size");
    // RDKit✔️✔️: auto *fpArguments = dynamic_cast<RDKitFPArguments *>(arguments);
    // RDKit✔️✔️: std::vector<AtomEnvironment<OutputType> *> result;
    // RDKit✔️✔️: INT_PATH_LIST_MAP allPaths;
    // RDKit✔️✔️: RDKitFPUtils::enumerateAllPaths(
    // RDKit✔️✔️:     mol, allPaths, fromAtoms, fpArguments->df_branchedPaths,
    // RDKit✔️✔️:     fpArguments->df_useHs, fpArguments->d_minPath, fpArguments->d_maxPath);
    let all_paths = enumerate_fingerprint_paths(
        molecule,
        params.min_path,
        params.max_path,
        params.use_hs,
        params.branched_paths,
        params.from_atoms.as_deref(),
    )?;

    // RDKit source: FingerprintUtil.cpp lines 291-315 (`identifyQueryBonds`).
    // RDKit✔️✔️: std::vector<short> isQueryBond(mol.getNumBonds(), 0);
    // RDKit✔️✔️: std::vector<const Bond *> bondCache;
    // RDKit✔️✔️: RDKitFPUtils::identifyQueryBonds(mol, bondCache, isQueryBond);
    // The typed molecule query fields are inspected directly by
    // `rdkit_fp_generate_bond_hash_inputs`, preserving the same path-local
    // query short circuit without a second cache representation.
    let mut result = Vec::new();
    for paths in all_paths.values() {
        for path in paths {
            // RDKit✔️✔️: std::vector<std::uint32_t> bondHashes = RDKitFPUtils::generateBondHashes(
            // RDKit✔️✔️:     mol, atomsInPath, bondCache, isQueryBond, path,
            // RDKit✔️✔️:     fpArguments->df_useBondOrder, atomInvariants);
            let hash_inputs = rdkit_fp_generate_bond_hash_inputs(
                molecule,
                path,
                params.use_bond_order,
                atom_invariants,
            )?;
            if hash_inputs.bond_hashes.is_empty() {
                continue;
            }

            // RDKit✔️✔️: unsigned long seed;
            // RDKit✔️✔️: if (path.size() > 1) {
            // RDKit✔️✔️:   std::sort(bondHashes.begin(), bondHashes.end());
            // RDKit✔️✔️:   bondHashes.push_back(static_cast<std::uint32_t>(atomsInPath.count()));
            // RDKit✔️✔️:   seed = gboost::hash_range(bondHashes.begin(), bondHashes.end());
            // RDKit✔️✔️: } else {
            // RDKit✔️✔️:   seed = bondHashes[0];
            // RDKit✔️✔️: }
            let seed = if path.len() > 1 {
                let mut bond_hashes = hash_inputs.bond_hashes;
                bond_hashes.sort_unstable();
                bond_hashes.push(
                    hash_inputs
                        .atoms_in_path
                        .iter()
                        .filter(|&&in_path| in_path)
                        .count() as u32,
                );
                hash_range(&bond_hashes)
            } else {
                hash_inputs.bond_hashes[0]
            };
            result.push(RdkitFpEnvironment {
                bit_id: u64::from(seed),
                atoms_in_path: hash_inputs.atoms_in_path,
                bond_path: path.clone(),
            });
        }
    }
    Ok(result)
}

/// Parameters for the source-backed `RDKFingerprintMol` boundary.
///
/// Parameters for the source-backed RDKit generator and legacy
/// `RDKFingerprintMol` boundary.
///
/// # Parameters
/// - `min_path`: minimum path length in bonds (default 1).
/// - `max_path`: maximum path length in bonds (default 7).
/// - `fp_size`: size of the output fingerprint (default 2048).
/// - `num_bits_per_feature`: number of bit positions to set per path hash
///   (default 2).
/// - `use_hs`: include explicit hydrogens in paths (default true).
/// - `target_density`: density threshold for source-compatible factor-of-two
///   folding (default 0.0).
/// - `min_size`: minimum folded size (default 128).
/// - `branched_paths`: include branched subgraphs as well as linear paths
///   (default true).
/// - `use_bond_order`: include bond order in path hashes (default true).
/// - `atom_invariants`: optional source-sized custom atom invariants.
/// - `from_atoms`: if `Some`, only enumerate paths starting from these atoms.
#[derive(Debug, Clone)]
pub struct TopologicalFingerprintParams {
    pub min_path: u32,
    pub max_path: u32,
    pub fp_size: u32,
    pub num_bits_per_feature: u32,
    pub use_hs: bool,
    pub target_density: f64,
    pub min_size: u32,
    pub branched_paths: bool,
    pub use_bond_order: bool,
    pub atom_invariants: Option<Vec<u32>>,
    pub from_atoms: Option<Vec<u32>>,
}

impl Default for TopologicalFingerprintParams {
    fn default() -> Self {
        Self {
            min_path: 1,
            max_path: 7,
            fp_size: 2048,
            num_bits_per_feature: 2,
            use_hs: true,
            target_density: 0.0,
            min_size: 128,
            branched_paths: true,
            use_bond_order: true,
            atom_invariants: None,
            from_atoms: None,
        }
    }
}

impl TopologicalFingerprintParams {
    pub fn validate(&self) -> Result<(), FingerprintError> {
        if self.min_path == 0 {
            return Err(FingerprintError::InvalidArguments {
                reason: "minPath==0",
            });
        }
        if self.max_path < self.min_path {
            return Err(FingerprintError::InvalidArguments {
                reason: "maxPath<minPath",
            });
        }
        if self.fp_size == 0 {
            return Err(FingerprintError::InvalidArguments {
                reason: "fpSize==0",
            });
        }
        if self.num_bits_per_feature == 0 {
            return Err(FingerprintError::InvalidArguments {
                reason: "nBitsPerHash==0",
            });
        }
        if !self.target_density.is_finite() || self.target_density < 0.0 {
            return Err(FingerprintError::InvalidArguments {
                reason: "tgtDensity must be finite and non-negative",
            });
        }
        Ok(())
    }
}

/// Typed request for the two legacy `RDKFingerprintMol` provenance outputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct TopologicalFingerprintOutputRequest {
    pub atom_bits: bool,
    pub bit_info: bool,
}

/// Typed provenance returned by the source `atomBits` and `bitInfo` outputs.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct TopologicalFingerprintOutput {
    pub atom_bits: Option<Vec<Vec<u32>>>,
    pub bit_info: Option<BTreeMap<u32, Vec<Vec<i32>>>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologicalFingerprintResult {
    pub fingerprint: Fingerprint,
    pub output: TopologicalFingerprintOutput,
}

// RDKit source: GraphMol/Fingerprints/Fingerprints.h::RDKFingerprintMol
// RDKit✔️✔️: RDKIT_FINGERPRINTS_EXPORT ExplicitBitVect *RDKFingerprintMol(
// RDKit✔️✔️:     const ROMol &mol, unsigned int minPath = 1, unsigned int maxPath = 7,
// RDKit✔️✔️:     unsigned int fpSize = 2048, unsigned int nBitsPerHash = 2,
// RDKit✔️✔️:     bool useHs = true, double tgtDensity = 0.0, unsigned int minSize = 128,
// RDKit✔️✔️:     bool branchedPaths = true, bool useBondOrder = true,
// RDKit✔️✔️:     std::vector<std::uint32_t> *atomInvariants = nullptr,
// RDKit✔️✔️:     const std::vector<std::uint32_t> *fromAtoms = nullptr,
// RDKit✔️✔️:     std::vector<std::vector<std::uint32_t>> *atomBits = nullptr,
// RDKit✔️✔️:     std::map<std::uint32_t, std::vector<std::vector<int>>> *bitInfo = nullptr);
fn rdkit_fp_fold_bits(on_bits: &mut BTreeSet<usize>, current_size: &mut usize) {
    // RDKit source: DataStructs/BitOps.cpp `FoldFingerprint`.
    // RDKit✔️✔️: for (unsigned int i = 0; i < nBits / factor; ++i) {
    // RDKit✔️✔️:   if (bv.getBit(i) || bv.getBit(i + nBits / factor)) {
    // RDKit✔️✔️:     res->setBit(i);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let folded_size = *current_size / 2;
    let folded = on_bits.iter().map(|&bit| bit % folded_size).collect();
    *on_bits = folded;
    *current_size = folded_size;
}

fn rdkit_fp_output_from_additional(
    additional_output: AdditionalOutput,
    request: TopologicalFingerprintOutputRequest,
) -> TopologicalFingerprintOutput {
    let atom_bits = request.atom_bits.then(|| {
        additional_output
            .atom_to_bits
            .unwrap_or_default()
            .into_iter()
            .map(|bits| bits.into_iter().map(|bit| bit as u32).collect())
            .collect()
    });
    let bit_info = request.bit_info.then(|| {
        additional_output
            .bit_paths
            .unwrap_or_default()
            .into_iter()
            .map(|(bit, paths)| {
                (
                    bit as u32,
                    paths
                        .into_iter()
                        .map(|path| path.into_iter().map(|bond| bond as i32).collect())
                        .collect(),
                )
            })
            .collect()
    });
    TopologicalFingerprintOutput {
        atom_bits,
        bit_info,
    }
}

pub fn topological_fingerprint(
    molecule: &Molecule,
    params: &TopologicalFingerprintParams,
) -> Result<Fingerprint, FingerprintError> {
    Ok(topological_fingerprint_with_output(molecule, params, Default::default())?.fingerprint)
}

pub fn topological_fingerprint_with_output(
    molecule: &Molecule,
    params: &TopologicalFingerprintParams,
    request: TopologicalFingerprintOutputRequest,
) -> Result<TopologicalFingerprintResult, FingerprintError> {
    params.validate()?;

    // RDKit source: Fingerprints.cpp lines 190-231 (`RDKFingerprintMol`).
    // RDKit✔️✔️: std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
    // RDKit✔️✔️:     RDKit::RDKitFP::getRDKitFPGenerator<std::uint32_t>(
    // RDKit✔️✔️:         minPath, maxPath, useHs, branchedPaths, useBondOrder));
    // RDKit✔️✔️: fpgen->getOptions()->d_fpSize = fpSize;
    // RDKit✔️✔️: fpgen->getOptions()->d_numBitsPerFeature = nBitsPerHash;
    // RDKit✔️✔️: FingerprintFuncArguments args;
    // RDKit✔️✔️: args.customAtomInvariants = atomInvariants;
    // RDKit✔️✔️: args.fromAtoms = fromAtoms;
    // RDKit✔️✔️: AdditionalOutput ao;
    // RDKit✔️✔️: if (atomBits) { args.additionalOutput = &ao; ao.allocateAtomToBits(); }
    // RDKit✔️✔️: if (bitInfo) { args.additionalOutput = &ao; ao.allocateBitPaths(); }
    let atom_invariants = params
        .atom_invariants
        .clone()
        .unwrap_or_else(|| rdkit_fp_atom_invariants(molecule));
    let mut additional_output = AdditionalOutput::new();
    if request.atom_bits {
        additional_output.allocate_atom_to_bits();
    }
    if request.bit_info {
        additional_output.allocate_bit_paths();
    }
    if request.atom_bits || request.bit_info {
        additional_output.reset_for_atom_count(molecule.num_atoms());
    }

    let environments = generate_rdkit_fp_environments(molecule, params, &atom_invariants)?;
    let mut on_bits = BTreeSet::new();
    let mut random_source =
        (params.num_bits_per_feature > 1).then(|| RdkitFingerprintMtRng::new(42));
    for environment in environments {
        let seed = environment.bit_id();
        let mut bit_id = (seed % u64::from(params.fp_size)) as usize;
        on_bits.insert(bit_id);
        if request.atom_bits || request.bit_info {
            environment.update_additional_output(&mut additional_output, bit_id as u64);
        }
        if let Some(random_source) = random_source.as_mut() {
            random_source.seed(seed as u32);
            for _ in 1..params.num_bits_per_feature {
                bit_id = (u64::from(random_source.uniform_int_0_to_i32_max())
                    % u64::from(params.fp_size)) as usize;
                on_bits.insert(bit_id);
                if request.atom_bits || request.bit_info {
                    environment.update_additional_output(&mut additional_output, bit_id as u64);
                }
            }
        }
    }

    // RDKit source: Fingerprints.cpp lines 231-242 (density folding).
    // RDKit✔️✔️: if (tgtDensity > 0.0) {
    // RDKit✔️✔️:   while (static_cast<double>(res->getNumOnBits()) / res->getNumBits() <
    // RDKit✔️✔️:              tgtDensity && res->getNumBits() >= 2 * minSize) {
    // RDKit✔️✔️:     ExplicitBitVect *tmpV = FoldFingerprint(*res, 2);
    // RDKit✔️✔️:     delete res;
    // RDKit✔️✔️:     res = tmpV;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let mut fingerprint_size = params.fp_size as usize;
    while params.target_density > 0.0
        && (on_bits.len() as f64) / (fingerprint_size as f64) < params.target_density
        && fingerprint_size >= (params.min_size as usize).saturating_mul(2)
    {
        rdkit_fp_fold_bits(&mut on_bits, &mut fingerprint_size);
    }

    let fingerprint = Fingerprint::from_on_bits(fingerprint_size, on_bits);
    Ok(TopologicalFingerprintResult {
        fingerprint,
        output: rdkit_fp_output_from_additional(additional_output, request),
    })
}

// ---------------------------------------------------------------------------
// MACCS Keys (166-bit structural keys)
// RDKit source: GraphMol/Fingerprints/MACCS.cpp
// ---------------------------------------------------------------------------

/// Parameters for MACCS key generation.
///
/// The MACCS (Molecular ACCess System) keys are 166 predefined structural
/// features encoded as a bit vector. Each bit corresponds to a specific
/// substructure or element property.
///
/// This surface is source-backed against RDKit `MACCSkeys.GenMACCSKeys` for
/// the exposed raw 167-bit vector and COSMolKit's public 166-bit projection.
///
/// # Reference
/// RDKit MACCS implementation: Code/GraphMol/Fingerprints/MACCS.cpp
/// - 166 SMARTS patterns mapped to bits 1-166 (bit 0 unused).
/// - Heuristic key detection is not accepted; the legacy local MACCS heuristic
///   implementation has been removed.
/// - COSMolKit exposes only the RDKit public 166-bit projection. Other output
///   lengths are not RDKit MACCS behavior and fail explicitly.
#[derive(Debug, Clone)]
pub struct MaccsFingerprintParams {
    pub n_bits: usize,
}

impl Default for MaccsFingerprintParams {
    fn default() -> Self {
        Self {
            n_bits: COSMOLKIT_MACCS_PUBLIC_BITS,
        }
    }
}

/// RDKit MACCS fingerprints allocate 167 raw bits and leave bit 0 unused.
pub(crate) const RDKIT_MACCS_RAW_BITS: usize = 167;
pub(crate) const COSMOLKIT_MACCS_PUBLIC_BITS: usize = RDKIT_MACCS_RAW_BITS - 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct MaccsPatternSpec {
    pub(crate) bit: usize,
    pub(crate) smarts: &'static str,
}

// RDKit source: MACCS.cpp::Patterns
// RDKit✔️✔️: struct Patterns {
// RDKit✔️✔️:   std::unique_ptr<RDKit::ROMol> bit_8 =
// RDKit✔️✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6!#1]1~*~*~*~1"));
// RDKit✔️✔️:   ...
// RDKit✔️✔️:   std::unique_ptr<RDKit::ROMol> bit_165 =
// RDKit✔️✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[R]"));
// RDKit✔️✔️: };
#[allow(dead_code)]
pub(crate) const RDKIT_MACCS_PATTERNS: &[MaccsPatternSpec] = &[
    MaccsPatternSpec {
        bit: 8,
        smarts: "[!#6!#1]1~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 11,
        smarts: "*1~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 13,
        smarts: "[#8]~[#7](~[#6])~[#6]",
    },
    MaccsPatternSpec {
        bit: 14,
        smarts: "[#16]-[#16]",
    },
    MaccsPatternSpec {
        bit: 15,
        smarts: "[#8]~[#6](~[#8])~[#8]",
    },
    MaccsPatternSpec {
        bit: 16,
        smarts: "[!#6!#1]1~*~*~1",
    },
    MaccsPatternSpec {
        bit: 17,
        smarts: "[#6]#[#6]",
    },
    MaccsPatternSpec {
        bit: 19,
        smarts: "*1~*~*~*~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 20,
        smarts: "[#14]",
    },
    MaccsPatternSpec {
        bit: 21,
        smarts: "[#6]=[#6](~[!#6!#1])~[!#6!#1]",
    },
    MaccsPatternSpec {
        bit: 22,
        smarts: "*1~*~*~1",
    },
    MaccsPatternSpec {
        bit: 23,
        smarts: "[#7]~[#6](~[#8])~[#8]",
    },
    MaccsPatternSpec {
        bit: 24,
        smarts: "[#7]-[#8]",
    },
    MaccsPatternSpec {
        bit: 25,
        smarts: "[#7]~[#6](~[#7])~[#7]",
    },
    MaccsPatternSpec {
        bit: 26,
        smarts: "[#6]=@[#6](@*)@*",
    },
    MaccsPatternSpec {
        bit: 28,
        smarts: "[!#6!#1]~[CH2]~[!#6!#1]",
    },
    MaccsPatternSpec {
        bit: 30,
        smarts: "[#6]~[!#6!#1](~[#6])(~[#6])~*",
    },
    MaccsPatternSpec {
        bit: 31,
        smarts: "[!#6!#1]~[F,Cl,Br,I]",
    },
    MaccsPatternSpec {
        bit: 32,
        smarts: "[#6]~[#16]~[#7]",
    },
    MaccsPatternSpec {
        bit: 33,
        smarts: "[#7]~[#16]",
    },
    MaccsPatternSpec {
        bit: 34,
        smarts: "[CH2]=*",
    },
    MaccsPatternSpec {
        bit: 36,
        smarts: "[#16R]",
    },
    MaccsPatternSpec {
        bit: 37,
        smarts: "[#7]~[#6](~[#8])~[#7]",
    },
    MaccsPatternSpec {
        bit: 38,
        smarts: "[#7]~[#6](~[#6])~[#7]",
    },
    MaccsPatternSpec {
        bit: 39,
        smarts: "[#8]~[#16](~[#8])~[#8]",
    },
    MaccsPatternSpec {
        bit: 40,
        smarts: "[#16]-[#8]",
    },
    MaccsPatternSpec {
        bit: 41,
        smarts: "[#6]#[#7]",
    },
    MaccsPatternSpec {
        bit: 43,
        smarts: "[!#6!#1!H0]~*~[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 44,
        smarts: "[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]",
    },
    MaccsPatternSpec {
        bit: 45,
        smarts: "[#6]=[#6]~[#7]",
    },
    MaccsPatternSpec {
        bit: 47,
        smarts: "[#16]~*~[#7]",
    },
    MaccsPatternSpec {
        bit: 48,
        smarts: "[#8]~[!#6!#1](~[#8])~[#8]",
    },
    MaccsPatternSpec {
        bit: 49,
        smarts: "[!+0]",
    },
    MaccsPatternSpec {
        bit: 50,
        smarts: "[#6]=[#6](~[#6])~[#6]",
    },
    MaccsPatternSpec {
        bit: 51,
        smarts: "[#6]~[#16]~[#8]",
    },
    MaccsPatternSpec {
        bit: 52,
        smarts: "[#7]~[#7]",
    },
    MaccsPatternSpec {
        bit: 53,
        smarts: "[!#6!#1!H0]~*~*~*~[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 54,
        smarts: "[!#6!#1!H0]~*~*~[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 55,
        smarts: "[#8]~[#16]~[#8]",
    },
    MaccsPatternSpec {
        bit: 56,
        smarts: "[#8]~[#7](~[#8])~[#6]",
    },
    MaccsPatternSpec {
        bit: 57,
        smarts: "[#8R]",
    },
    MaccsPatternSpec {
        bit: 58,
        smarts: "[!#6!#1]~[#16]~[!#6!#1]",
    },
    MaccsPatternSpec {
        bit: 59,
        smarts: "[#16]!:*:*",
    },
    MaccsPatternSpec {
        bit: 60,
        smarts: "[#16]=[#8]",
    },
    MaccsPatternSpec {
        bit: 61,
        smarts: "*~[#16](~*)~*",
    },
    MaccsPatternSpec {
        bit: 62,
        smarts: "*@*!@*@*",
    },
    MaccsPatternSpec {
        bit: 63,
        smarts: "[#7]=[#8]",
    },
    MaccsPatternSpec {
        bit: 64,
        smarts: "*@*!@[#16]",
    },
    MaccsPatternSpec {
        bit: 65,
        smarts: "c:n",
    },
    MaccsPatternSpec {
        bit: 66,
        smarts: "[#6]~[#6](~[#6])(~[#6])~*",
    },
    MaccsPatternSpec {
        bit: 67,
        smarts: "[!#6!#1]~[#16]",
    },
    MaccsPatternSpec {
        bit: 68,
        smarts: "[!#6!#1!H0]~[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 69,
        smarts: "[!#6!#1]~[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 70,
        smarts: "[!#6!#1]~[#7]~[!#6!#1]",
    },
    MaccsPatternSpec {
        bit: 71,
        smarts: "[#7]~[#8]",
    },
    MaccsPatternSpec {
        bit: 72,
        smarts: "[#8]~*~*~[#8]",
    },
    MaccsPatternSpec {
        bit: 73,
        smarts: "[#16]=*",
    },
    MaccsPatternSpec {
        bit: 74,
        smarts: "[CH3]~*~[CH3]",
    },
    MaccsPatternSpec {
        bit: 75,
        smarts: "*!@[#7]@*",
    },
    MaccsPatternSpec {
        bit: 76,
        smarts: "[#6]=[#6](~*)~*",
    },
    MaccsPatternSpec {
        bit: 77,
        smarts: "[#7]~*~[#7]",
    },
    MaccsPatternSpec {
        bit: 78,
        smarts: "[#6]=[#7]",
    },
    MaccsPatternSpec {
        bit: 79,
        smarts: "[#7]~*~*~[#7]",
    },
    MaccsPatternSpec {
        bit: 80,
        smarts: "[#7]~*~*~*~[#7]",
    },
    MaccsPatternSpec {
        bit: 81,
        smarts: "[#16]~*(~*)~*",
    },
    MaccsPatternSpec {
        bit: 82,
        smarts: "*~[CH2]~[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 83,
        smarts: "[!#6!#1]1~*~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 84,
        smarts: "[NH2]",
    },
    MaccsPatternSpec {
        bit: 85,
        smarts: "[#6]~[#7](~[#6])~[#6]",
    },
    MaccsPatternSpec {
        bit: 86,
        smarts: "[C;H2,H3][!#6!#1][C;H2,H3]",
    },
    MaccsPatternSpec {
        bit: 87,
        smarts: "[F,Cl,Br,I]!@*@*",
    },
    MaccsPatternSpec {
        bit: 89,
        smarts: "[#8]~*~*~*~[#8]",
    },
    MaccsPatternSpec {
        bit: 90,
        smarts: "[$([!#6!#1!H0]~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[CH2R]1)]",
    },
    MaccsPatternSpec {
        bit: 91,
        smarts: "[$([!#6!#1!H0]~*~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~*~[R]1@[R]@[CH2R]1)]",
    },
    MaccsPatternSpec {
        bit: 92,
        smarts: "[#8]~[#6](~[#7])~[#6]",
    },
    MaccsPatternSpec {
        bit: 93,
        smarts: "[!#6!#1]~[CH3]",
    },
    MaccsPatternSpec {
        bit: 94,
        smarts: "[!#6!#1]~[#7]",
    },
    MaccsPatternSpec {
        bit: 95,
        smarts: "[#7]~*~*~[#8]",
    },
    MaccsPatternSpec {
        bit: 96,
        smarts: "*1~*~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 97,
        smarts: "[#7]~*~*~*~[#8]",
    },
    MaccsPatternSpec {
        bit: 98,
        smarts: "[!#6!#1]1~*~*~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 99,
        smarts: "[#6]=[#6]",
    },
    MaccsPatternSpec {
        bit: 100,
        smarts: "*~[CH2]~[#7]",
    },
    MaccsPatternSpec {
        bit: 101,
        smarts: "[$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1)]",
    },
    MaccsPatternSpec {
        bit: 102,
        smarts: "[!#6!#1]~[#8]",
    },
    MaccsPatternSpec {
        bit: 104,
        smarts: "[!#6!#1!H0]~*~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 105,
        smarts: "*@*(@*)@*",
    },
    MaccsPatternSpec {
        bit: 106,
        smarts: "[!#6!#1]~*(~[!#6!#1])~[!#6!#1]",
    },
    MaccsPatternSpec {
        bit: 107,
        smarts: "[F,Cl,Br,I]~*(~*)~*",
    },
    MaccsPatternSpec {
        bit: 108,
        smarts: "[CH3]~*~*~*~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 109,
        smarts: "*~[CH2]~[#8]",
    },
    MaccsPatternSpec {
        bit: 110,
        smarts: "[#7]~[#6]~[#8]",
    },
    MaccsPatternSpec {
        bit: 111,
        smarts: "[#7]~*~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 112,
        smarts: "*~*(~*)(~*)~*",
    },
    MaccsPatternSpec {
        bit: 113,
        smarts: "[#8]!:*:*",
    },
    MaccsPatternSpec {
        bit: 114,
        smarts: "[CH3]~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 115,
        smarts: "[CH3]~*~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 116,
        smarts: "[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]",
    },
    MaccsPatternSpec {
        bit: 117,
        smarts: "[#7]~*~[#8]",
    },
    MaccsPatternSpec {
        bit: 118,
        smarts: "[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]",
    },
    MaccsPatternSpec {
        bit: 119,
        smarts: "[#7]=*",
    },
    MaccsPatternSpec {
        bit: 120,
        smarts: "[!#6R]",
    },
    MaccsPatternSpec {
        bit: 121,
        smarts: "[#7R]",
    },
    MaccsPatternSpec {
        bit: 122,
        smarts: "*~[#7](~*)~*",
    },
    MaccsPatternSpec {
        bit: 123,
        smarts: "[#8]~[#6]~[#8]",
    },
    MaccsPatternSpec {
        bit: 124,
        smarts: "[!#6!#1]~[!#6!#1]",
    },
    MaccsPatternSpec {
        bit: 126,
        smarts: "*!@[#8]!@*",
    },
    MaccsPatternSpec {
        bit: 127,
        smarts: "*@*!@[#8]",
    },
    MaccsPatternSpec {
        bit: 128,
        smarts: "[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2R]@[R]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2R]1),$(*~[CH2]~*~[R]1@[R]@[CH2R]1)]",
    },
    MaccsPatternSpec {
        bit: 129,
        smarts: "[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[CH2R]1)]",
    },
    MaccsPatternSpec {
        bit: 131,
        smarts: "[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 132,
        smarts: "[#8]~*~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 133,
        smarts: "*@*!@[#7]",
    },
    MaccsPatternSpec {
        bit: 135,
        smarts: "[#7]!:*:*",
    },
    MaccsPatternSpec {
        bit: 136,
        smarts: "[#8]=*",
    },
    MaccsPatternSpec {
        bit: 137,
        smarts: "[!C!cR]",
    },
    MaccsPatternSpec {
        bit: 138,
        smarts: "[!#6!#1]~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 139,
        smarts: "[O!H0]",
    },
    MaccsPatternSpec {
        bit: 140,
        smarts: "[#8]",
    },
    MaccsPatternSpec {
        bit: 141,
        smarts: "[CH3]",
    },
    MaccsPatternSpec {
        bit: 142,
        smarts: "[#7]",
    },
    MaccsPatternSpec {
        bit: 144,
        smarts: "*!:*:*!:*",
    },
    MaccsPatternSpec {
        bit: 145,
        smarts: "*1~*~*~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 147,
        smarts: "[$(*~[CH2]~[CH2]~*),$([R]1@[CH2R]@[CH2R]1)]",
    },
    MaccsPatternSpec {
        bit: 148,
        smarts: "*~[!#6!#1](~*)~*",
    },
    MaccsPatternSpec {
        bit: 149,
        smarts: "[C;H3,H4]",
    },
    MaccsPatternSpec {
        bit: 150,
        smarts: "*!@*@*!@*",
    },
    MaccsPatternSpec {
        bit: 151,
        smarts: "[#7!H0]",
    },
    MaccsPatternSpec {
        bit: 152,
        smarts: "[#8]~[#6](~[#6])~[#6]",
    },
    MaccsPatternSpec {
        bit: 154,
        smarts: "[#6]=[#8]",
    },
    MaccsPatternSpec {
        bit: 155,
        smarts: "*!@[CH2]!@*",
    },
    MaccsPatternSpec {
        bit: 156,
        smarts: "[#7]~*(~*)~*",
    },
    MaccsPatternSpec {
        bit: 157,
        smarts: "[#6]-[#8]",
    },
    MaccsPatternSpec {
        bit: 158,
        smarts: "[#6]-[#7]",
    },
    MaccsPatternSpec {
        bit: 162,
        smarts: "a",
    },
    MaccsPatternSpec {
        bit: 165,
        smarts: "[R]",
    },
];

#[allow(dead_code)]
pub(crate) fn rdkit_maccs_pattern(bit: usize) -> Option<&'static MaccsPatternSpec> {
    RDKIT_MACCS_PATTERNS
        .iter()
        .find(|pattern| pattern.bit == bit)
}

#[allow(dead_code)]
pub(crate) fn rdkit_maccs_public_index(bit: usize) -> Option<usize> {
    (1..RDKIT_MACCS_RAW_BITS).contains(&bit).then(|| bit - 1)
}

fn rdkit_maccs_generate_fp_raw_bits_001_166(
    molecule: &Molecule,
) -> Result<BTreeSet<usize>, FingerprintError> {
    let mut bits = BTreeSet::new();
    let pattern_matchers = cached_rdkit_maccs_pattern_matchers()?;

    // RDKit source: MACCS.cpp::GenerateFP direct element-key loop.
    // RDKit✔️✔️:   for (atom = mol.beginAtoms(); atom != mol.endAtoms(); ++atom) {
    // RDKit✔️✔️:     switch ((*atom)->getAtomicNum()) {
    for atom in molecule.atoms() {
        match atom.atomic_number() {
            // RDKit✔️✔️:       case 3:
            // RDKit✔️✔️:       case 11:
            // RDKit✔️✔️:       case 19:
            // RDKit✔️✔️:       case 37:
            // RDKit✔️✔️:       case 55:
            // RDKit✔️✔️:       case 87:
            // RDKit✔️✔️:         fp.setBit(35);
            3 | 11 | 19 | 37 | 55 | 87 => {
                bits.insert(35);
            }
            // RDKit✔️✔️:       case 4:
            // RDKit✔️✔️:       case 12:
            // RDKit✔️✔️:       case 20:
            // RDKit✔️✔️:       case 38:
            // RDKit✔️✔️:       case 56:
            // RDKit✔️✔️:       case 88:
            // RDKit✔️✔️:         fp.setBit(10);
            4 | 12 | 20 | 38 | 56 | 88 => {
                bits.insert(10);
            }
            // RDKit✔️✔️:       case 5:
            // RDKit✔️✔️:       case 13:
            // RDKit✔️✔️:       case 31:
            // RDKit✔️✔️:       case 49:
            // RDKit✔️✔️:       case 81:
            // RDKit✔️✔️:         fp.setBit(18);
            5 | 13 | 31 | 49 | 81 => {
                bits.insert(18);
            }
            // RDKit✔️✔️:       case 9:
            // RDKit✔️✔️:         fp.setBit(42);
            // RDKit✔️✔️:         fp.setBit(134);
            9 => {
                bits.insert(42);
                bits.insert(134);
            }
            // RDKit✔️✔️:       case 15:
            // RDKit✔️✔️:         fp.setBit(29);
            15 => {
                bits.insert(29);
            }
            // RDKit✔️✔️:       case 16:
            // RDKit✔️✔️:         fp.setBit(88);
            16 => {
                bits.insert(88);
            }
            // RDKit✔️✔️:       case 17:
            // RDKit✔️✔️:         fp.setBit(103);
            // RDKit✔️✔️:         fp.setBit(134);
            17 => {
                bits.insert(103);
                bits.insert(134);
            }
            // RDKit✔️✔️:       case 21:
            // RDKit✔️✔️:       case 22:
            // RDKit✔️✔️:       case 39:
            // RDKit✔️✔️:       case 40:
            // RDKit✔️✔️:       case 72:
            // RDKit✔️✔️:         fp.setBit(5);
            21 | 22 | 39 | 40 | 72 => {
                bits.insert(5);
            }
            // RDKit✔️✔️:       case 23:
            // RDKit✔️✔️:       case 24:
            // RDKit✔️✔️:       case 25:
            // RDKit✔️✔️:       case 41:
            // RDKit✔️✔️:       case 42:
            // RDKit✔️✔️:       case 43:
            // RDKit✔️✔️:       case 73:
            // RDKit✔️✔️:       case 74:
            // RDKit✔️✔️:       case 75:
            // RDKit✔️✔️:         fp.setBit(7);
            23 | 24 | 25 | 41 | 42 | 43 | 73 | 74 | 75 => {
                bits.insert(7);
            }
            // RDKit✔️✔️:       case 26:
            // RDKit✔️✔️:       case 27:
            // RDKit✔️✔️:       case 28:
            // RDKit✔️✔️:       case 44:
            // RDKit✔️✔️:       case 45:
            // RDKit✔️✔️:       case 46:
            // RDKit✔️✔️:       case 76:
            // RDKit✔️✔️:       case 77:
            // RDKit✔️✔️:       case 78:
            // RDKit✔️✔️:         fp.setBit(9);
            26 | 27 | 28 | 44 | 45 | 46 | 76 | 77 | 78 => {
                bits.insert(9);
            }
            // RDKit✔️✔️:       case 29:
            // RDKit✔️✔️:       case 30:
            // RDKit✔️✔️:       case 47:
            // RDKit✔️✔️:       case 48:
            // RDKit✔️✔️:       case 79:
            // RDKit✔️✔️:       case 80:
            // RDKit✔️✔️:         fp.setBit(12);
            29 | 30 | 47 | 48 | 79 | 80 => {
                bits.insert(12);
            }
            // RDKit✔️✔️:       case 32:
            // RDKit✔️✔️:       case 33:
            // RDKit✔️✔️:       case 34:
            // RDKit✔️✔️:       case 50:
            // RDKit✔️✔️:       case 51:
            // RDKit✔️✔️:       case 52:
            // RDKit✔️✔️:       case 82:
            // RDKit✔️✔️:       case 83:
            // RDKit✔️✔️:       case 84:
            // RDKit✔️✔️:         fp.setBit(3);
            32 | 33 | 34 | 50 | 51 | 52 | 82 | 83 | 84 => {
                bits.insert(3);
            }
            // RDKit✔️✔️:       case 35:
            // RDKit✔️✔️:         fp.setBit(46);
            // RDKit✔️✔️:         fp.setBit(134);
            35 => {
                bits.insert(46);
                bits.insert(134);
            }
            // RDKit✔️✔️:       case 53:
            // RDKit✔️✔️:         fp.setBit(27);
            // RDKit✔️✔️:         fp.setBit(134);
            53 => {
                bits.insert(27);
                bits.insert(134);
            }
            // RDKit✔️✔️:       case 57:
            // RDKit✔️✔️:       case 58:
            // RDKit✔️✔️:       case 59:
            // RDKit✔️✔️:       case 60:
            // RDKit✔️✔️:       case 61:
            // RDKit✔️✔️:       case 62:
            // RDKit✔️✔️:       case 63:
            // RDKit✔️✔️:       case 64:
            // RDKit✔️✔️:       case 65:
            // RDKit✔️✔️:       case 66:
            // RDKit✔️✔️:       case 67:
            // RDKit✔️✔️:       case 68:
            // RDKit✔️✔️:       case 69:
            // RDKit✔️✔️:       case 70:
            // RDKit✔️✔️:       case 71:
            // RDKit✔️✔️:         fp.setBit(6);
            57..=71 => {
                bits.insert(6);
            }
            // RDKit✔️✔️:       case 89:
            // RDKit✔️✔️:       case 90:
            // RDKit✔️✔️:       case 91:
            // RDKit✔️✔️:       case 92:
            // RDKit✔️✔️:       case 93:
            // RDKit✔️✔️:       case 94:
            // RDKit✔️✔️:       case 95:
            // RDKit✔️✔️:       case 96:
            // RDKit✔️✔️:       case 97:
            // RDKit✔️✔️:       case 98:
            // RDKit✔️✔️:       case 99:
            // RDKit✔️✔️:       case 100:
            // RDKit✔️✔️:       case 101:
            // RDKit✔️✔️:       case 102:
            // RDKit✔️✔️:       case 103:
            // RDKit✔️✔️:         fp.setBit(4);
            89..=103 => {
                bits.insert(4);
            }
            // RDKit✔️✔️:       case 104:
            // RDKit✔️✔️:         fp.setBit(2);
            104 => {
                bits.insert(2);
            }
            _ => {}
        }
    }

    // RDKit source: MACCS.cpp::GenerateFP pattern-key checks through bit_120.
    // RDKit✔️✔️:   if (RDKit::SubstructMatch(mol, *pats.bit_8, match, true)) {
    // RDKit✔️✔️:     fp.setBit(8);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ...
    // RDKit✔️✔️:   if (RDKit::SubstructMatch(mol, *pats.bit_117, match, true)) {
    // RDKit✔️✔️:     fp.setBit(117);
    // RDKit✔️✔️:   }
    for &(bit, ref matcher) in pattern_matchers
        .iter()
        .filter(|&&(bit, _)| bit <= 120 && bit != 118 && bit != 120)
    {
        if crate::has_substruct_match(molecule, matcher.getMatcher()) {
            bits.insert(bit);
        }
    }

    // RDKit✔️✔️:   if (RDKit::SubstructMatch(mol, *pats.bit_118, matches, true, true) > 1) {
    // RDKit✔️✔️:     fp.setBit(118);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (RDKit::SubstructMatch(mol, *pats.bit_120, matches, true, true) > 1) {
    // RDKit✔️✔️:     fp.setBit(120);
    // RDKit✔️✔️:   }
    for bit in [118usize, 120] {
        if let Some(pattern) = pattern_matchers
            .iter()
            .find(|&&(pattern_bit, _)| pattern_bit == bit)
            && crate::get_substruct_matches(molecule, pattern.1.getMatcher()).len() > 1
        {
            bits.insert(bit);
        }
    }

    let maccs_match_count = |bit: usize| -> usize {
        pattern_matchers
            .iter()
            .find(|&&(pattern_bit, _)| pattern_bit == bit)
            .map_or(0, |pattern| {
                crate::get_substruct_matches(molecule, pattern.1.getMatcher()).len()
            })
    };

    let has_maccs_match = |bit: usize| -> bool {
        pattern_matchers
            .iter()
            .find(|&&(pattern_bit, _)| pattern_bit == bit)
            .is_some_and(|pattern| crate::has_substruct_match(molecule, pattern.1.getMatcher()))
    };

    // RDKit source: MACCS.cpp::GenerateFP pattern-key checks from bit_121 through bit_165.
    // RDKit✔️✔️:   if (RDKit::SubstructMatch(mol, *pats.bit_121, match, true)) {
    // RDKit✔️✔️:     fp.setBit(121);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ...
    // RDKit✔️✔️:   if (RDKit::SubstructMatch(mol, *pats.bit_165, match, true)) {
    // RDKit✔️✔️:     fp.setBit(165);
    // RDKit✔️✔️:   }
    for bit in [
        121usize, 122, 123, 126, 128, 129, 132, 133, 135, 137, 139, 144, 147, 148, 150, 151, 152,
        154, 155, 156, 157, 158, 162, 165,
    ] {
        if has_maccs_match(bit) {
            bits.insert(bit);
        }
    }

    // RDKit✔️✔️:   count = RDKit::SubstructMatch(mol, *pats.bit_124, matches, true, true);
    // RDKit✔️✔️:   if (count > 0) {
    // RDKit✔️✔️:     fp.setBit(124);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (count > 1) {
    // RDKit✔️✔️:     fp.setBit(130);
    // RDKit✔️✔️:   }
    let count = maccs_match_count(124);
    if count > 0 {
        bits.insert(124);
    }
    if count > 1 {
        bits.insert(130);
    }

    // RDKit✔️✔️:   count = RDKit::SubstructMatch(mol, *pats.bit_127, matches, true, true);
    // RDKit✔️✔️:   if (count > 1) {
    // RDKit✔️✔️:     fp.setBit(127);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (count > 0) {
    // RDKit✔️✔️:     fp.setBit(143);
    // RDKit✔️✔️:   }
    let count = maccs_match_count(127);
    if count > 1 {
        bits.insert(127);
    }
    if count > 0 {
        bits.insert(143);
    }

    // RDKit✔️✔️:   if (RDKit::SubstructMatch(mol, *pats.bit_131, matches, true, true) > 1) {
    // RDKit✔️✔️:     fp.setBit(131);
    // RDKit✔️✔️:   }
    if maccs_match_count(131) > 1 {
        bits.insert(131);
    }

    // RDKit✔️✔️:   if (RDKit::SubstructMatch(mol, *pats.bit_136, matches, true, true) > 1) {
    // RDKit✔️✔️:     fp.setBit(136);
    // RDKit✔️✔️:   }
    if maccs_match_count(136) > 1 {
        bits.insert(136);
    }

    // RDKit✔️✔️:   count = RDKit::SubstructMatch(mol, *pats.bit_138, matches, true, true);
    // RDKit✔️✔️:   if (count > 1) {
    // RDKit✔️✔️:     fp.setBit(138);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (count > 0) {
    // RDKit✔️✔️:     fp.setBit(153);
    // RDKit✔️✔️:   }
    let count = maccs_match_count(138);
    if count > 1 {
        bits.insert(138);
    }
    if count > 0 {
        bits.insert(153);
    }

    // RDKit✔️✔️:   count = RDKit::SubstructMatch(mol, *pats.bit_140, matches, true, true);
    // RDKit✔️✔️:   if (count > 3) {
    // RDKit✔️✔️:     fp.setBit(140);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (count > 2) {
    // RDKit✔️✔️:     fp.setBit(146);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (count > 1) {
    // RDKit✔️✔️:     fp.setBit(159);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (count > 0) {
    // RDKit✔️✔️:     fp.setBit(164);
    // RDKit✔️✔️:   }
    let count = maccs_match_count(140);
    if count > 3 {
        bits.insert(140);
    }
    if count > 2 {
        bits.insert(146);
    }
    if count > 1 {
        bits.insert(159);
    }
    if count > 0 {
        bits.insert(164);
    }

    // RDKit✔️✔️:   if (RDKit::SubstructMatch(mol, *pats.bit_141, matches, true, true) > 2) {
    // RDKit✔️✔️:     fp.setBit(141);
    // RDKit✔️✔️:   }
    if maccs_match_count(141) > 2 {
        bits.insert(141);
    }

    // RDKit✔️✔️:   count = RDKit::SubstructMatch(mol, *pats.bit_142, matches, true, true);
    // RDKit✔️✔️:   if (count > 1) {
    // RDKit✔️✔️:     fp.setBit(142);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (count > 0) {
    // RDKit✔️✔️:     fp.setBit(161);
    // RDKit✔️✔️:   }
    let count = maccs_match_count(142);
    if count > 1 {
        bits.insert(142);
    }
    if count > 0 {
        bits.insert(161);
    }

    // RDKit✔️✔️:   count = RDKit::SubstructMatch(mol, *pats.bit_145, matches, true, true);
    // RDKit✔️✔️:   if (count > 1) {
    // RDKit✔️✔️:     fp.setBit(145);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (count > 0) {
    // RDKit✔️✔️:     fp.setBit(163);
    // RDKit✔️✔️:   }
    let count = maccs_match_count(145);
    if count > 1 {
        bits.insert(145);
    }
    if count > 0 {
        bits.insert(163);
    }

    // RDKit✔️✔️:   count = RDKit::SubstructMatch(mol, *pats.bit_149, matches, true, true);
    // RDKit✔️✔️:   if (count > 1) {
    // RDKit✔️✔️:     fp.setBit(149);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (count > 0) {
    // RDKit✔️✔️:     fp.setBit(160);
    // RDKit✔️✔️:   }
    let count = maccs_match_count(149);
    if count > 1 {
        bits.insert(149);
    }
    if count > 0 {
        bits.insert(160);
    }

    // RDKit✔️✔️:   /* BIT 125 */
    // RDKit✔️❌:   RDKit::RingInfo *info = mol.getRingInfo();
    // RDKit✔️❌:   unsigned int ringcount = info->numRings();
    // RDKit✔️❌:   unsigned int nArom = 0;
    // RDKit✔️❌:   for (unsigned int i = 0; i < ringcount; i++) {
    // RDKit✔️❌:     bool isArom = true;
    // RDKit✔️❌:     const std::vector<int> *ring = &info->bondRings()[i];
    // RDKit✔️❌:     std::vector<int>::const_iterator iter;
    // RDKit✔️❌:     for (iter = ring->begin(); iter != ring->end(); ++iter) {
    // RDKit✔️❌:       if (!mol.getBondWithIdx(*iter)->getIsAromatic()) {
    // RDKit✔️❌:         isArom = false;
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (isArom) {
    // RDKit✔️❌:       if (nArom) {
    // RDKit✔️❌:         fp.setBit(125);
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         nArom++;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    //
    // COSMolKit may compute SSSR when no cached ring info is present; this preserves
    // the observed bit semantics but can add one ring-finding pass relative to RDKit's
    // initialized RingInfo lookup.
    let owned_rings = molecule
        .derived_cache()
        .rings
        .is_none()
        .then(|| crate::symmetrize_sssr(molecule).ok())
        .flatten();
    let ring_info = molecule
        .derived_cache()
        .rings
        .as_ref()
        .or(owned_rings.as_ref());
    if let Some(ring_info) = ring_info {
        let mut aromatic_ring_count = 0usize;
        for bond_ring in ring_info.bond_rings() {
            if bond_ring.iter().all(|bond| {
                molecule
                    .bonds()
                    .get(bond.index())
                    .is_some_and(crate::Bond::is_aromatic)
            }) {
                if aromatic_ring_count > 0 {
                    bits.insert(125);
                    break;
                }
                aromatic_ring_count += 1;
            }
        }
    }

    // RDKit✔️✔️:   /* BIT 166 */
    // RDKit✔️✔️:   std::vector<int> mapping;
    // RDKit✔️✔️:   if (RDKit::MolOps::getMolFrags(mol, mapping) > 1) {
    // RDKit✔️✔️:     fp.setBit(166);
    // RDKit✔️✔️:   }
    let mut fragment_count = 0usize;
    let mut visited = vec![false; molecule.num_atoms()];
    let adjacency = &molecule.topology_block().adjacency;
    for atom_idx in 0..molecule.num_atoms() {
        if visited[atom_idx] {
            continue;
        }
        fragment_count += 1;
        if fragment_count > 1 {
            bits.insert(166);
            break;
        }
        visited[atom_idx] = true;
        let mut queue = VecDeque::from([atom_idx]);
        while let Some(current) = queue.pop_front() {
            for neighbor in adjacency.neighbors_of(current) {
                if !visited[neighbor.atom_index] {
                    visited[neighbor.atom_index] = true;
                    queue.push_back(neighbor.atom_index);
                }
            }
        }
    }

    Ok(bits)
}

// RDKit✔️❌: ExplicitBitVect *getFingerprintAsBitVect(const ROMol &mol)
// RDKit✔️❌:   std::unique_ptr<ExplicitBitVect> fp(new ExplicitBitVect(167));
// RDKit✔️❌:   GenerateFP(mol, *fp);
// RDKit✔️❌:   return fp.release();
// END RDKIT FUNCTION getFingerprintAsBitVect
//
// COSMolKit computes the RDKit raw 167-bit vector, then projects RDKit raw bits
// 1..166 onto public indices 0..165. Bit 0 remains unused in the raw vector.
// The raw implementation can compute SSSR ring info when the molecule has no
// cached ring cache, so performance is marked worse than RDKit's direct
// initialized RingInfo access for that branch.
pub fn maccs_get_fingerprint_as_bit_vect(
    molecule: &Molecule,
) -> Result<Fingerprint, FingerprintError> {
    Ok(Fingerprint::from_on_bits(
        RDKIT_MACCS_RAW_BITS,
        rdkit_maccs_generate_fp_raw_bits_001_166(molecule)?,
    ))
}

pub fn maccs_fingerprint(
    molecule: &Molecule,
    params: &MaccsFingerprintParams,
) -> Result<Fingerprint, FingerprintError> {
    let n_bits = params.n_bits;
    if n_bits != COSMOLKIT_MACCS_PUBLIC_BITS {
        return Err(FingerprintError::UnsupportedOption {
            option: "MaccsFingerprintParams.n_bits",
            reason: "RDKit MACCS exposes a fixed 167-bit raw vector with bit 0 unused; COSMolKit only exposes the exact 166-bit public projection",
        });
    }
    if molecule.num_atoms() == 0 {
        return Ok(Fingerprint::from_on_bits(n_bits, []));
    }

    let projected_on_bits = maccs_get_fingerprint_as_bit_vect(molecule)?
        .on_bits()
        .into_iter()
        .filter_map(rdkit_maccs_public_index)
        .filter(|&public_bit| public_bit < COSMOLKIT_MACCS_PUBLIC_BITS);
    Ok(Fingerprint::from_on_bits(n_bits, projected_on_bits))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod topological_torsion_code_tests;

#[cfg(test)]
mod topological_torsion_path_tests;

#[cfg(test)]
mod topological_torsion_arguments_tests;

#[cfg(test)]
mod topological_torsion_environment_tests;

#[cfg(test)]
mod topological_torsion_generator_tests;

#[cfg(test)]
mod topological_torsion_legacy_tests;

#[cfg(test)]
mod topological_torsion_public_api_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        AtomQueryPredicate, AtomSpec, BondQueryPredicate, BondSpec, BondStereo, Element, Molecule,
        QueryNode,
    };
    use serde::Deserialize;

    fn default_morgan_params(radius: u32, n_bits: usize) -> MorganFingerprintParams {
        MorganFingerprintParams {
            radius,
            n_bits,
            ..Default::default()
        }
    }

    #[test]
    fn topological_fingerprint_empty_molecule_matches_source() {
        let fingerprint =
            topological_fingerprint(&Molecule::new(), &TopologicalFingerprintParams::default())
                .expect("empty molecule fingerprint");
        assert_eq!(fingerprint.n_bits(), 2048);
        assert!(fingerprint.on_bits().is_empty());
    }

    #[test]
    fn topological_params_match_rdkit_boundary_defaults() {
        let params = TopologicalFingerprintParams::default();
        assert_eq!(params.min_path, 1);
        assert_eq!(params.max_path, 7);
        assert_eq!(params.fp_size, 2048);
        assert_eq!(params.num_bits_per_feature, 2);
        assert!(params.use_hs);
        assert_eq!(params.target_density, 0.0);
        assert_eq!(params.min_size, 128);
        assert!(params.branched_paths);
        assert!(params.use_bond_order);
        assert!(params.atom_invariants.is_none());
        assert!(params.from_atoms.is_none());
    }

    #[test]
    fn topological_params_reject_source_precondition_ranges() {
        let mut params = TopologicalFingerprintParams::default();
        params.min_path = 0;
        assert!(matches!(
            params.validate(),
            Err(FingerprintError::InvalidArguments {
                reason: "minPath==0"
            })
        ));
        params = TopologicalFingerprintParams::default();
        params.max_path = 0;
        assert!(matches!(
            params.validate(),
            Err(FingerprintError::InvalidArguments {
                reason: "maxPath<minPath"
            })
        ));
        params = TopologicalFingerprintParams::default();
        params.fp_size = 0;
        assert!(matches!(
            params.validate(),
            Err(FingerprintError::InvalidArguments {
                reason: "fpSize==0"
            })
        ));
        params = TopologicalFingerprintParams::default();
        params.num_bits_per_feature = 0;
        assert!(matches!(
            params.validate(),
            Err(FingerprintError::InvalidArguments {
                reason: "nBitsPerHash==0"
            })
        ));
    }

    #[test]
    fn topological_typed_provenance_request_allocates_empty_source_outputs() {
        let request = TopologicalFingerprintOutputRequest {
            atom_bits: true,
            bit_info: true,
        };
        let result = topological_fingerprint_with_output(
            &Molecule::new(),
            &TopologicalFingerprintParams::default(),
            request,
        )
        .expect("source provenance output");
        assert_eq!(result.output.atom_bits, Some(Vec::new()));
        assert_eq!(result.output.bit_info, Some(BTreeMap::new()));
    }

    #[test]
    fn rdkit_fp_atom_invariants_use_atomic_number_and_aromatic_flag() {
        let aliphatic = Molecule::from_smiles("CCO").expect("aliphatic fixture");
        assert_eq!(rdkit_fp_atom_invariants(&aliphatic), vec![12, 12, 16]);
        let aromatic = Molecule::from_smiles("c1ccccc1").expect("aromatic fixture");
        assert_eq!(rdkit_fp_atom_invariants(&aromatic), vec![13; 6]);
    }

    #[test]
    fn layered_boost_hash_range_matches_empty_order_and_overflow() {
        assert_eq!(hash_range(&[]), 0x0000_0000);
        assert_eq!(hash_range(&[0]), 0x9e37_79b9);
        assert_eq!(hash_range(&[1]), 0x9e37_79ba);
        assert_eq!(hash_range(&[1, 2]), 0xcd94_bf13);
        assert_eq!(hash_range(&[2, 1]), 0xcd94_bf53);
        assert_eq!(hash_range(&[u32::MAX; 4]), 0x4841_0b19);
        assert_eq!(
            hash_range(&[u32::MAX, 0, 0x8000_0000, 0x7fff_ffff]),
            0x6841_6bc8
        );
    }

    #[test]
    fn layered_boost_hash_range_is_pinned_to_uint32_not_native_width() {
        let source_width_hash = hash_range(&[0x1234_5678, 0x9abc_def0, u32::MAX]);
        assert_eq!(source_width_hash, 0xf41b_442d);

        // Applying the same expression to a native-width 64-bit seed retains
        // high bits and therefore does not reproduce RDKit's `hash_result_t`.
        let native_width_hash = 0x0000_0b74_341b_442du64;
        assert_ne!(u64::from(source_width_hash), native_width_hash);
    }

    #[test]
    fn rdkit_fp_bond_hash_inputs_match_boost_hash_combine_order() {
        let molecule = Molecule::from_smiles("CCO").expect("fixture");
        let invariants = rdkit_fp_atom_invariants(&molecule);
        let first = rdkit_fp_generate_bond_hash_inputs(&molecule, &[0], true, &invariants)
            .expect("first bond hash");
        assert_eq!(first.atoms_in_path, vec![true, true, false]);
        assert_eq!(first.bond_hashes, vec![4_275_705_116]);
        let second = rdkit_fp_generate_bond_hash_inputs(&molecule, &[1], true, &invariants)
            .expect("second bond hash");
        assert_eq!(second.bond_hashes, vec![4_274_652_475]);
        let without_order = rdkit_fp_generate_bond_hash_inputs(&molecule, &[0], false, &invariants)
            .expect("order-free bond hash");
        assert_eq!(first.bond_hashes, without_order.bond_hashes);
        let double_bond = Molecule::from_smiles("C=C").expect("double-bond fixture");
        let double_invariants = rdkit_fp_atom_invariants(&double_bond);
        let with_order =
            rdkit_fp_generate_bond_hash_inputs(&double_bond, &[0], true, &double_invariants)
                .expect("double bond hash");
        let without_order =
            rdkit_fp_generate_bond_hash_inputs(&double_bond, &[0], false, &double_invariants)
                .expect("double bond order-free hash");
        assert_ne!(with_order.bond_hashes, without_order.bond_hashes);
    }

    #[test]
    fn rdkit_fp_query_bond_path_is_rejected_at_hash_input_boundary() {
        let query =
            crate::search::smarts_parse::compile_query_fixture("[#6]~[#6]").expect("query fixture");
        let invariants = rdkit_fp_atom_invariants(&query);
        let inputs = rdkit_fp_generate_bond_hash_inputs(&query, &[0], true, &invariants)
            .expect("query path result");
        assert!(inputs.bond_hashes.is_empty());
        assert_eq!(inputs.atoms_in_path, vec![true, true]);
    }

    fn expected_paths(
        molecule: &Molecule,
        min_path: u32,
        max_path: u32,
        use_hs: bool,
        branched_paths: bool,
        from_atoms: Option<&[u32]>,
    ) -> BTreeMap<usize, Vec<Vec<usize>>> {
        enumerate_fingerprint_paths(
            molecule,
            min_path,
            max_path,
            use_hs,
            branched_paths,
            from_atoms,
        )
        .expect("valid source path arguments")
    }

    #[test]
    fn rdkit_fp_linear_and_branched_paths_match_source_order() {
        let linear = Molecule::from_smiles("CCCC").expect("linear fixture");
        assert_eq!(
            expected_paths(&linear, 1, 3, true, false, None),
            BTreeMap::from([
                (1, vec![vec![0], vec![1], vec![2]]),
                (2, vec![vec![0, 1], vec![1, 2]]),
                (3, vec![vec![0, 1, 2]]),
            ])
        );

        let branched = Molecule::from_smiles("CC(C)C").expect("branched fixture");
        assert_eq!(
            expected_paths(&branched, 1, 3, true, false, None),
            BTreeMap::from([
                (1, vec![vec![0], vec![1], vec![2]]),
                (2, vec![vec![0, 1], vec![0, 2], vec![1, 2]]),
            ])
        );
        assert_eq!(
            expected_paths(&branched, 1, 3, true, true, None),
            BTreeMap::from([
                (1, vec![vec![0], vec![1], vec![2]]),
                (2, vec![vec![0, 2], vec![0, 1], vec![1, 2]]),
                (3, vec![vec![0, 2, 1]]),
            ])
        );
    }

    #[test]
    fn layered_shared_path_enumeration_distinguishes_absent_and_empty_roots() {
        let molecule = Molecule::from_smiles("CC(C)C").expect("root-state fixture");
        for branched_paths in [false, true] {
            let absent = expected_paths(&molecule, 1, 3, true, branched_paths, None);
            let empty = expected_paths(&molecule, 1, 3, true, branched_paths, Some(&[]));
            assert!(absent.values().any(|paths| !paths.is_empty()));
            assert_eq!(empty, BTreeMap::new());
        }
    }

    #[test]
    fn layered_shared_path_enumeration_preserves_root_duplicates_and_prepend_order() {
        fn prepend_groups(
            groups: &[&BTreeMap<usize, Vec<Vec<usize>>>],
        ) -> BTreeMap<usize, Vec<Vec<usize>>> {
            let mut result = BTreeMap::new();
            for group in groups {
                for (&length, paths) in *group {
                    result
                        .entry(length)
                        .or_insert_with(Vec::new)
                        .splice(0..0, paths.iter().cloned());
                }
            }
            result
        }

        let molecule = Molecule::from_smiles("CC(C)C").expect("root-order fixture");
        for branched_paths in [false, true] {
            let root_zero = expected_paths(&molecule, 1, 3, true, branched_paths, Some(&[0]));
            let root_two = expected_paths(&molecule, 1, 3, true, branched_paths, Some(&[2]));

            assert_eq!(
                expected_paths(&molecule, 1, 3, true, branched_paths, Some(&[0, 0])),
                prepend_groups(&[&root_zero, &root_zero])
            );
            assert_eq!(
                expected_paths(&molecule, 1, 3, true, branched_paths, Some(&[0, 2])),
                prepend_groups(&[&root_zero, &root_two])
            );
            assert_eq!(
                expected_paths(&molecule, 1, 3, true, branched_paths, Some(&[2, 0])),
                prepend_groups(&[&root_two, &root_zero])
            );
        }
    }

    #[test]
    fn layered_preparation_rejects_each_source_precondition() {
        let molecule = Molecule::from_smiles("CCO").expect("precondition fixture");
        let error =
            |result: Result<LayeredFingerprintPreparation<'_>, FingerprintError>| match result
                .expect_err("source precondition must fail")
            {
                FingerprintError::InvalidArguments { reason } => reason,
                other => panic!("unexpected error: {other}"),
            };

        assert_eq!(
            error(prepare_layered_fingerprint(
                &molecule, 0, 7, 2048, None, None
            )),
            "minPath==0"
        );
        assert_eq!(
            error(prepare_layered_fingerprint(
                &molecule, 2, 1, 2048, None, None
            )),
            "maxPath<minPath"
        );
        assert_eq!(
            error(prepare_layered_fingerprint(&molecule, 1, 7, 0, None, None)),
            "fpSize==0"
        );
        assert_eq!(
            error(prepare_layered_fingerprint(
                &molecule,
                1,
                7,
                2048,
                Some(&[0, 0]),
                None,
            )),
            "bad atomCounts size"
        );
        let wrong_width = Fingerprint::from_on_bits(1024, []);
        assert_eq!(
            error(prepare_layered_fingerprint(
                &molecule,
                1,
                7,
                2048,
                None,
                Some(&wrong_width),
            )),
            "bad setOnlyBits size"
        );
    }

    #[test]
    fn layered_preparation_reuses_initialized_rings_or_computes_exact_sssr() {
        let uncached = Molecule::new();
        let prepared = prepare_layered_fingerprint(&uncached, 1, 7, 2048, None, None)
            .expect("uncached exact SSSR preparation");
        assert!(matches!(prepared.ring_info, Cow::Owned(_)));
        assert!(prepared.ring_info.is_initialized());

        let cached = Molecule::from_smiles_with_sanitize("C1CC1", false)
            .expect("cached-ring fixture")
            .with_assigned_rings()
            .expect("materialize initialized ring cache");
        let prepared = prepare_layered_fingerprint(&cached, 1, 7, 2048, None, None)
            .expect("cached ring preparation");
        assert!(matches!(prepared.ring_info, Cow::Borrowed(_)));
        assert_eq!(prepared.ring_info.num_rings(), 1);
    }

    #[test]
    fn layered_preparation_builds_source_ordered_atom_bond_and_query_caches() {
        let molecule = Molecule::from_smiles("CCO").expect("cache fixture");
        let prepared = prepare_layered_fingerprint(&molecule, 1, 7, 2048, None, None)
            .expect("cache preparation");
        assert_eq!(
            prepared
                .bond_cache
                .iter()
                .map(|bond| bond.id().index())
                .collect::<Vec<_>>(),
            vec![0, 1]
        );
        assert_eq!(prepared.query_masks, vec![0, 0]);
        assert_eq!(prepared.atomic_numbers, vec![6, 6, 8]);
        assert_eq!(prepared.aromatic_atoms, vec![false, false, false]);

        let aromatic = Molecule::from_smiles("c1ccccc1").expect("aromatic cache fixture");
        let prepared = prepare_layered_fingerprint(&aromatic, 1, 7, 2048, None, None)
            .expect("aromatic cache preparation");
        assert_eq!(prepared.atomic_numbers, vec![6; 6]);
        assert_eq!(prepared.aromatic_atoms, vec![true; 6]);

        let mut builder = Molecule::builder();
        let begin = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_query(QueryNode::predicate(AtomQueryPredicate::FormalCharge(0))),
        );
        let end = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_query(QueryNode::predicate(AtomQueryPredicate::FormalCharge(0))),
        );
        builder
            .add_bond(
                BondSpec::new(begin, end, BondOrder::Single)
                    .with_query(QueryNode::predicate(BondQueryPredicate::Any)),
            )
            .expect("query-mask bond");
        let query = builder.build().expect("query-mask molecule");
        let prepared = prepare_layered_fingerprint(&query, 1, 7, 2048, None, None)
            .expect("query-mask preparation");
        assert_eq!(prepared.query_masks, vec![0x7]);
    }

    #[derive(Debug, Deserialize)]
    struct LayeredQueryFixture {
        schema_version: u32,
        reference: LayeredQueryReference,
        parameters: LayeredQueryParameters,
        complexity_masks: Vec<LayeredQueryCase>,
        aromaticity_branches: Vec<LayeredQueryCase>,
    }

    #[derive(Debug, Deserialize)]
    struct LayeredQueryReference {
        name: String,
        version: String,
        source_revision: String,
        source_paths: Vec<String>,
    }

    #[derive(Debug, Deserialize)]
    struct LayeredQueryParameters {
        layer_flags: u32,
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        branched_paths: bool,
    }

    #[derive(Debug, Deserialize)]
    struct LayeredQueryCase {
        case_id: String,
        notation: String,
        input: String,
        query_masks: Vec<u8>,
        aromatic_atoms: Vec<bool>,
        on_bits: Vec<usize>,
    }

    #[test]
    fn layered_query_parity_matches_source_masks_and_aromaticity() {
        let fixture_path = cosmolkit_test_support::repo_root()
            .join("testdata/fingerprint/fixtures/rdkit/layered_fingerprint_query_cases.json");
        let fixture_text = std::fs::read_to_string(&fixture_path)
            .unwrap_or_else(|error| panic!("failed to read {}: {error}", fixture_path.display()));
        let fixture: LayeredQueryFixture = serde_json::from_str(&fixture_text)
            .unwrap_or_else(|error| panic!("failed to parse {}: {error}", fixture_path.display()));

        assert_eq!(fixture.schema_version, 1);
        assert_eq!(fixture.reference.name, "RDKit");
        assert_eq!(fixture.reference.version, "2026.3.1");
        assert_eq!(
            fixture.reference.source_revision,
            "351f8f378f8ad6bbd517980c38896e66bf907af8c"
        );
        assert!(
            fixture
                .reference
                .source_paths
                .iter()
                .any(|path| path.contains("Fingerprints/test1.cpp"))
        );
        assert!(
            fixture
                .reference
                .source_paths
                .iter()
                .any(|path| path.contains("QueryOps.cpp"))
        );

        let params = LayeredFingerprintParams {
            layers: LayeredFingerprintLayers::from_bits_retain(fixture.parameters.layer_flags),
            min_path: fixture.parameters.min_path,
            max_path: fixture.parameters.max_path,
            fp_size: fixture.parameters.fp_size,
            branched_paths: fixture.parameters.branched_paths,
            ..LayeredFingerprintParams::default()
        };

        let assert_case = |case: &LayeredQueryCase| {
            let molecule = match case.notation.as_str() {
                "smiles" => Molecule::from_smiles(&case.input).unwrap_or_else(|error| {
                    panic!("{} ({}) failed to parse: {error}", case.case_id, case.input)
                }),
                "smarts" => mol_from_smarts(&case.input, &SmartsParseParams::default())
                    .unwrap_or_else(|error| {
                        panic!("{} ({}) failed to parse: {error}", case.case_id, case.input)
                    }),
                notation => panic!("{}: unsupported notation {notation}", case.case_id),
            };
            let prepared = prepare_layered_fingerprint(
                &molecule,
                params.min_path,
                params.max_path,
                params.fp_size as usize,
                None,
                None,
            )
            .unwrap_or_else(|error| panic!("{} preparation failed: {error}", case.case_id));
            assert_eq!(
                prepared.query_masks, case.query_masks,
                "{}: exact three-bit bond/endpoint query masks",
                case.case_id
            );
            assert_eq!(
                prepared.aromatic_atoms, case.aromatic_atoms,
                "{}: exact query-aware aromaticity cache",
                case.case_id
            );
            let fingerprint = layered_fingerprint(&molecule, &params)
                .unwrap_or_else(|error| panic!("{} fingerprint failed: {error}", case.case_id));
            assert_eq!(
                fingerprint.on_bits(),
                case.on_bits,
                "{}: exact complete Layered fingerprint",
                case.case_id
            );
        };

        for case in &fixture.complexity_masks {
            assert_case(case);
        }
        assert_eq!(
            fixture
                .complexity_masks
                .iter()
                .map(|case| case.query_masks[0])
                .collect::<Vec<_>>(),
            (0u8..=7).collect::<Vec<_>>(),
            "fixture must cover every source three-bit mask exactly once"
        );
        for case in &fixture.aromaticity_branches {
            assert_case(case);
        }

        // The committed SMARTS fixture covers every source-representable
        // root. These typed cases close the remaining QueryOps dispatch
        // branches without widening the private query API just for tests.
        let query_aromaticity = |query, aromatic| {
            let mut builder = Molecule::builder();
            let atom_id = builder.add_atom(
                AtomSpec::new(Element::C)
                    .with_aromatic(aromatic)
                    .with_query(query),
            );
            let molecule = builder.build().expect("typed aromaticity query fixture");
            crate::search::query::is_atom_aromatic(&molecule.atoms()[atom_id.index()], &molecule)
        };
        let number = || QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));
        let aromatic = || QueryNode::predicate(AtomQueryPredicate::IsAromatic(true));
        let formal_charge = || QueryNode::predicate(AtomQueryPredicate::FormalCharge(0));

        assert!(query_aromaticity(number(), true));
        assert!(!query_aromaticity(QueryNode::or(vec![aromatic()]), true));
        assert!(!query_aromaticity(
            QueryNode::xor(vec![aromatic(), number()]),
            true,
        ));
        assert!(!query_aromaticity(
            QueryNode::and(vec![aromatic(), number()]),
            true,
        ));
        assert!(!query_aromaticity(
            QueryNode::and(vec![number(), formal_charge()]),
            true,
        ));
        assert!(!query_aromaticity(
            QueryNode::not(QueryNode::and(vec![number(), aromatic()])),
            true,
        ));
        assert!(!query_aromaticity(
            QueryNode::not(QueryNode::not(aromatic())),
            true,
        ));
        assert!(!query_aromaticity(formal_charge(), true));
    }

    #[test]
    fn layered_topology_bond_order_encoders_match_source_packing_and_modulo() {
        fn one_bond(order: BondOrder, aromatic: bool) -> Molecule {
            let mut builder = Molecule::builder();
            let begin = builder.add_atom(AtomSpec::new(Element::C));
            let end = builder.add_atom(AtomSpec::new(Element::C));
            builder
                .add_bond(BondSpec::new(begin, end, order).with_aromatic(aromatic))
                .expect("encoder fixture bond");
            builder.build().expect("encoder fixture")
        }

        assert_eq!(layered_topology_hash(9, 2, 5), 169);
        assert_eq!(layered_topology_hash(9, 5, 2), 169);
        assert_eq!(layered_topology_hash(u32::MAX, u32::MAX, u32::MAX), 511);

        let double = one_bond(BondOrder::Double, false);
        assert_eq!(
            layered_bond_order_hash(&double.bonds()[0], 9, 2, 5, 0),
            Some(1_354)
        );
        assert_eq!(
            layered_bond_order_hash(&double.bonds()[0], 9, 5, 2, 0),
            Some(1_354)
        );
        assert_eq!(
            layered_bond_order_hash(&double.bonds()[0], u32::MAX, u32::MAX, u32::MAX, 0,),
            Some(4_090)
        );

        let quadruple = one_bond(BondOrder::Quadruple, false);
        assert_eq!(
            layered_bond_order_hash(&quadruple.bonds()[0], 9, 2, 5, 0),
            Some(1_356)
        );
        let aromatic_flag = one_bond(BondOrder::Double, true);
        assert_eq!(
            layered_bond_order_hash(&aromatic_flag.bonds()[0], 9, 2, 5, 0),
            Some(1_353)
        );
    }

    #[test]
    fn layered_topology_bond_order_encoder_suppresses_only_complex_bond_queries() {
        let molecule = Molecule::from_smiles("C=C").expect("suppression fixture");
        let bond = &molecule.bonds()[0];
        assert_eq!(layered_bond_order_hash(bond, 0, 1, 1, 0x1), None);
        assert_eq!(layered_bond_order_hash(bond, 0, 1, 1, 0x7), None);
        assert_eq!(layered_bond_order_hash(bond, 0, 1, 1, 0x6), Some(578));
        assert_eq!(layered_topology_hash(0, 1, 1), 72);
    }

    #[test]
    fn layered_atom_type_aromaticity_encoders_match_endpoint_order_and_packing() {
        assert_eq!(layered_atom_type_hash(8, 6, 2, 5, 9, 0), Some(1_737_480));
        assert_eq!(layered_atom_type_hash(6, 8, 5, 2, 9, 0), Some(1_737_480));
        assert_eq!(layered_atom_type_hash(6, 6, 2, 5, 9, 0), Some(1_393_414));
        assert_eq!(
            layered_atom_type_hash(u32::MAX, u32::MAX, u32::MAX, u32::MAX, u32::MAX, 0,),
            Some(8_388_607)
        );

        assert_eq!(layered_aromaticity_hash(true, false, 9, 0), Some(33));
        assert_eq!(layered_aromaticity_hash(false, true, 9, 0), Some(33));
        assert_eq!(layered_aromaticity_hash(true, true, 9, 0), Some(35));
        assert_eq!(layered_aromaticity_hash(false, false, 9, 0), Some(32));
        assert_eq!(layered_aromaticity_hash(true, true, u32::MAX, 0), Some(227));
    }

    #[test]
    fn layered_atom_type_aromaticity_encoders_use_query_aware_cache_and_suppression() {
        assert_eq!(layered_atom_type_hash(6, 8, 1, 1, 0, 0x2), None);
        assert_eq!(layered_atom_type_hash(6, 8, 1, 1, 0, 0x4), None);
        assert_eq!(layered_atom_type_hash(6, 8, 1, 1, 0, 0x1), Some(148_232));
        assert_eq!(layered_aromaticity_hash(true, false, 0, 0x6), None);
        assert_eq!(layered_aromaticity_hash(true, false, 0, 0x1), Some(1));

        let mut builder = Molecule::builder();
        let aromatic = builder.add_atom(AtomSpec::new(Element::C).with_query(
            QueryNode::predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: true,
            }),
        ));
        let aliphatic = builder.add_atom(AtomSpec::new(Element::C).with_query(
            QueryNode::predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            }),
        ));
        builder
            .add_bond(BondSpec::new(aromatic, aliphatic, BondOrder::Single))
            .expect("query-aware aromaticity bond");
        let molecule = builder.build().expect("query-aware aromaticity fixture");
        let prepared = prepare_layered_fingerprint(&molecule, 1, 7, 2048, None, None)
            .expect("query-aware aromaticity preparation");
        assert_eq!(prepared.query_masks, vec![0]);
        assert_eq!(prepared.aromatic_atoms, vec![true, false]);
        assert_eq!(
            layered_aromaticity_hash(
                prepared.aromatic_atoms[0],
                prepared.aromatic_atoms[1],
                0,
                prepared.query_masks[0],
            ),
            Some(1)
        );
    }

    #[test]
    fn layered_ring_encoders_distinguish_omission_zero_and_source_sssr_size() {
        let chain = Molecule::from_smiles("CC").expect("acyclic ring fixture");
        let chain_rings = crate::rings::find_sssr(&chain).expect("acyclic exact SSSR");
        let chain_bond = &chain.bonds()[0];
        assert_eq!(
            layered_ring_presence_hash(chain_bond, &chain_rings, 0),
            None
        );
        assert_eq!(
            layered_min_ring_size_hash(chain_bond, &chain_rings, 0),
            Some(0)
        );

        let triangle = Molecule::from_smiles("C1CC1").expect("ring fixture");
        let triangle_rings = crate::rings::find_sssr(&triangle).expect("ring exact SSSR");
        assert_eq!(
            layered_ring_presence_hash(&triangle.bonds()[0], &triangle_rings, 0),
            Some(1)
        );
        assert_eq!(
            layered_min_ring_size_hash(&triangle.bonds()[0], &triangle_rings, 0),
            Some(3)
        );

        let nine_ring = Molecule::from_smiles("C1CCCCCCCC1").expect("modulo ring fixture");
        let nine_ring_info = crate::rings::find_sssr(&nine_ring).expect("modulo exact SSSR");
        assert_eq!(
            layered_min_ring_size_hash(&nine_ring.bonds()[0], &nine_ring_info, 0),
            Some(1)
        );
    }

    #[test]
    fn layered_ring_encoders_match_fused_membership_and_endpoint_query_suppression() {
        let fused = Molecule::from_smiles("c1ccc2ccccc2c1").expect("fused-ring fixture");
        let ring_info = crate::rings::find_sssr(&fused).expect("fused exact SSSR");
        let shared_bond = fused
            .bonds()
            .iter()
            .find(|bond| ring_info.num_bond_rings(bond.id()) > 1)
            .expect("fused exact SSSR shared bond");
        assert_eq!(ring_info.min_bond_ring_size(shared_bond.id()), 6);
        assert_eq!(
            layered_ring_presence_hash(shared_bond, &ring_info, 0),
            Some(1)
        );
        assert_eq!(
            layered_min_ring_size_hash(shared_bond, &ring_info, 0),
            Some(6)
        );
        assert_eq!(
            layered_ring_presence_hash(shared_bond, &ring_info, 0x1),
            Some(1)
        );
        assert_eq!(
            layered_ring_presence_hash(shared_bond, &ring_info, 0x2),
            None
        );
        assert_eq!(
            layered_min_ring_size_hash(shared_bond, &ring_info, 0x4),
            None
        );
    }

    #[test]
    fn layered_projection_sorts_suffixes_hashes_and_applies_sparse_mask() {
        let mut layers = vec![Vec::new(); 10];
        layers[0] = vec![3, 1, 2];
        layers[1] = vec![7];
        let atoms_in_path = vec![true, true, true];
        let mask = Fingerprint::from_on_bits(1000, [146]);
        let mut words = vec![0u64; 1000usize.div_ceil(64)];
        let mut counts = vec![4, 10, u32::MAX];

        project_layered_path(
            &mut layers,
            &atoms_in_path,
            1000,
            Some(&mask),
            &mut words,
            Some(&mut counts),
        );

        assert_eq!(layers[0], vec![1, 2, 3, 3, 1]);
        assert_eq!(hash_range(&layers[0]), 0xfee9_dcd9);
        assert_eq!(hash_range(&layers[0]) % 1000, 289);
        assert_eq!(layers[1], vec![7, 3, 2]);
        assert_eq!(hash_range(&layers[1]), 0xfb5d_90da);
        assert_eq!(hash_range(&layers[1]) % 1000, 146);
        assert_eq!(
            Fingerprint {
                bits: words,
                n_bits: 1000,
            }
            .on_bits(),
            vec![146]
        );
        assert_eq!(counts, vec![5, 11, 0]);
    }

    #[test]
    fn layered_projection_counts_once_per_accepted_path_across_collisions_and_repeats() {
        let atoms_in_path = vec![true, false, true];
        let allow_all = Fingerprint::from_on_bits(1, [0]);
        let deny_all = Fingerprint::from_on_bits(1, []);
        let mut words = vec![0u64; 1];
        let mut counts = vec![0, 5, 0];

        let mut colliding_layers = vec![Vec::new(); 10];
        colliding_layers[0] = vec![1];
        colliding_layers[5] = vec![1];
        project_layered_path(
            &mut colliding_layers,
            &atoms_in_path,
            1,
            Some(&allow_all),
            &mut words,
            Some(&mut counts),
        );
        assert_eq!(words, vec![1]);
        assert_eq!(counts, vec![1, 5, 1]);

        let mut repeated_path_layers = vec![Vec::new(); 10];
        repeated_path_layers[0] = vec![1];
        project_layered_path(
            &mut repeated_path_layers,
            &atoms_in_path,
            1,
            None,
            &mut words,
            Some(&mut counts),
        );
        assert_eq!(counts, vec![2, 5, 2]);

        let mut rejected_layers = vec![Vec::new(); 10];
        rejected_layers[0] = vec![1];
        project_layered_path(
            &mut rejected_layers,
            &atoms_in_path,
            1,
            Some(&deny_all),
            &mut words,
            Some(&mut counts),
        );
        assert_eq!(counts, vec![2, 5, 2]);
    }

    #[test]
    fn layered_end_to_end_defaults_and_active_layers_match_pinned_source() {
        let molecule = Molecule::from_smiles("CCO").expect("end-to-end fixture");
        let original = molecule.clone();
        let params = LayeredFingerprintParams {
            atom_counts: Some(vec![0; 3]),
            ..LayeredFingerprintParams::default()
        };
        let result = layered_fingerprint_with_output(&molecule, &params)
            .expect("default Layered fingerprint");
        assert_eq!(result.fingerprint.n_bits(), 2048);
        assert_eq!(
            result.fingerprint.on_bits(),
            vec![92, 360, 596, 610, 611, 674, 867, 1044, 1111, 1783, 1784]
        );
        assert_eq!(result.atom_counts, Some(vec![2, 3, 2]));
        assert_eq!(molecule, original);

        let active = LayeredFingerprintParams {
            layers: LayeredFingerprintLayers::ACTIVE,
            ..params.clone()
        };
        assert_eq!(
            layered_fingerprint(&molecule, &active).expect("active Layered fingerprint"),
            result.fingerprint
        );
        assert_eq!(
            molecule
                .layered_fingerprint(&params)
                .expect("Molecule Layered method"),
            result.fingerprint
        );
    }

    #[test]
    fn layered_end_to_end_layers_masks_counts_and_roots_match_pinned_source() {
        let molecule = Molecule::from_smiles("CCO").expect("option fixture");
        let topology = LayeredFingerprintParams {
            layers: LayeredFingerprintLayers::TOPOLOGY,
            atom_counts: Some(vec![0; 3]),
            ..LayeredFingerprintParams::default()
        };
        let result = layered_fingerprint_with_output(&molecule, &topology)
            .expect("topology-only Layered fingerprint");
        assert_eq!(result.fingerprint.on_bits(), vec![674, 867]);
        assert_eq!(result.atom_counts, Some(vec![2, 3, 2]));

        let masked = LayeredFingerprintParams {
            atom_counts: Some(vec![10, 20, 30]),
            set_only_bits: Some(Fingerprint::from_on_bits(2048, [674])),
            ..LayeredFingerprintParams::default()
        };
        let result = layered_fingerprint_with_output(&molecule, &masked)
            .expect("masked Layered fingerprint");
        assert_eq!(result.fingerprint.on_bits(), vec![674]);
        assert_eq!(result.atom_counts, Some(vec![11, 22, 31]));

        let empty_roots = LayeredFingerprintParams {
            atom_counts: Some(vec![0; 3]),
            from_atoms: Some(Vec::new()),
            ..LayeredFingerprintParams::default()
        };
        let result = layered_fingerprint_with_output(&molecule, &empty_roots)
            .expect("present-empty root selection");
        assert!(result.fingerprint.on_bits().is_empty());
        assert_eq!(result.atom_counts, Some(vec![0, 0, 0]));

        let rooted_linear = LayeredFingerprintParams {
            atom_counts: Some(vec![0; 3]),
            branched_paths: false,
            from_atoms: Some(vec![0]),
            ..LayeredFingerprintParams::default()
        };
        let result = layered_fingerprint_with_output(&molecule, &rooted_linear)
            .expect("rooted linear Layered fingerprint");
        assert_eq!(
            result.fingerprint.on_bits(),
            vec![360, 596, 610, 611, 674, 867, 1044, 1111, 1783, 1784]
        );
        assert_eq!(result.atom_counts, Some(vec![2, 2, 1]));
    }

    #[test]
    fn layered_public_api_preserves_source_metadata_high_bits_and_errors() {
        let defaults = LayeredFingerprintParams::default();
        assert_eq!(defaults.layers.bits(), u32::MAX);
        assert_eq!(defaults.min_path, 1);
        assert_eq!(defaults.max_path, 7);
        assert_eq!(defaults.fp_size, 2048);
        assert!(defaults.branched_paths);
        assert_eq!(LAYERED_FINGERPRINT_MAX_LAYERS, 10);
        assert_eq!(LAYERED_FINGERPRINT_VERSION, "0.7.0");
        assert_eq!(LayeredFingerprintLayers::SUBSTRUCTURE.bits(), 0x07);

        let molecule = Molecule::from_smiles("CCO").expect("metadata fixture");
        let high_only = LayeredFingerprintParams {
            layers: LayeredFingerprintLayers::from_bits_retain(0xffff_ffc0),
            atom_counts: Some(vec![5, 6, 7]),
            ..LayeredFingerprintParams::default()
        };
        let result = layered_fingerprint_with_output(&molecule, &high_only)
            .expect("source-accepted high flags");
        assert!(result.fingerprint.on_bits().is_empty());
        assert_eq!(result.atom_counts, Some(vec![5, 6, 7]));

        let invalid_root = LayeredFingerprintParams {
            from_atoms: Some(vec![3]),
            ..LayeredFingerprintParams::default()
        };
        assert!(matches!(
            layered_fingerprint(&molecule, &invalid_root),
            Err(FingerprintError::InvalidArguments {
                reason: "fromAtoms contains atom index out of range"
            })
        ));
        assert_eq!(
            crate::LAYERED_FINGERPRINT_FEATURE.status,
            crate::SupportStatus::Experimental
        );
        assert!(
            crate::PUBLIC_FEATURES
                .iter()
                .any(|feature| **feature == crate::LAYERED_FINGERPRINT_FEATURE)
        );
    }

    #[test]
    fn rdkit_fp_ring_and_fused_ring_paths_preserve_bond_order() {
        let ring = Molecule::from_smiles("C1CCCCC1").expect("ring fixture");
        assert_eq!(
            expected_paths(&ring, 1, 3, true, false, None),
            BTreeMap::from([
                (
                    1,
                    vec![vec![0], vec![5], vec![1], vec![2], vec![3], vec![4]]
                ),
                (
                    2,
                    vec![
                        vec![0, 1],
                        vec![5, 4],
                        vec![0, 5],
                        vec![1, 2],
                        vec![2, 3],
                        vec![3, 4]
                    ]
                ),
                (
                    3,
                    vec![
                        vec![0, 1, 2],
                        vec![5, 4, 3],
                        vec![0, 5, 4],
                        vec![1, 2, 3],
                        vec![1, 0, 5],
                        vec![2, 3, 4]
                    ]
                ),
            ])
        );
        let fused = Molecule::from_smiles("c1ccc2ccccc2c1").expect("fused-ring fixture");
        assert_eq!(
            expected_paths(&fused, 2, 2, true, false, None)[&2],
            vec![
                vec![0, 1],
                vec![9, 8],
                vec![0, 9],
                vec![1, 2],
                vec![2, 3],
                vec![2, 10],
                vec![3, 4],
                vec![10, 7],
                vec![10, 8],
                vec![3, 10],
                vec![4, 5],
                vec![5, 6],
                vec![6, 7],
                vec![7, 8],
            ]
        );
    }

    #[test]
    fn rdkit_fp_paths_handle_disconnected_and_restricted_roots() {
        let disconnected = Molecule::from_smiles("CC.CC").expect("disconnected fixture");
        assert_eq!(
            expected_paths(&disconnected, 1, 3, true, false, None),
            BTreeMap::from([(1, vec![vec![0], vec![1]])])
        );

        let chain = Molecule::from_smiles("CCCC").expect("restricted fixture");
        assert_eq!(
            expected_paths(&chain, 1, 3, true, false, Some(&[2])),
            BTreeMap::from([(1, vec![vec![1], vec![2]]), (2, vec![vec![1, 0]]),])
        );
        assert_eq!(
            expected_paths(&chain, 1, 3, true, true, Some(&[2])),
            BTreeMap::from([
                (1, vec![vec![1], vec![2]]),
                (2, vec![vec![1, 2], vec![1, 0]]),
                (3, vec![vec![1, 2, 0]]),
            ])
        );
        let invalid_root = expected_paths(&chain, 1, 2, true, false, Some(&[99]));
        assert_eq!(invalid_root, BTreeMap::new());
    }

    #[test]
    fn rdkit_fp_explicit_hydrogen_filter_matches_source() {
        let molecule = Molecule::from_smiles("CC")
            .expect("explicit-H fixture")
            .with_hydrogens()
            .expect("materialize explicit hydrogens");
        assert_eq!(molecule.num_atoms(), 8);
        assert_eq!(
            expected_paths(&molecule, 1, 1, false, false, None)[&1],
            vec![vec![0]]
        );
        assert_eq!(
            expected_paths(&molecule, 1, 1, true, false, None)[&1],
            vec![
                vec![0],
                vec![1],
                vec![2],
                vec![3],
                vec![4],
                vec![5],
                vec![6]
            ]
        );
    }

    #[test]
    fn rdkit_fp_environment_hashes_and_additional_output_match_source() {
        let molecule = Molecule::from_smiles("CCO").expect("environment fixture");
        let params = TopologicalFingerprintParams {
            min_path: 1,
            max_path: 2,
            branched_paths: false,
            ..TopologicalFingerprintParams::default()
        };
        let invariants = rdkit_fp_atom_invariants(&molecule);
        let environments = generate_rdkit_fp_environments(&molecule, &params, &invariants)
            .expect("environment generation");
        assert_eq!(
            environments
                .iter()
                .map(RdkitFpEnvironment::bit_id)
                .collect::<Vec<_>>(),
            vec![4_275_705_116, 4_274_652_475, 1_524_090_560]
        );
        assert_eq!(
            environments
                .iter()
                .map(|environment| environment.bond_path.clone())
                .collect::<Vec<_>>(),
            vec![vec![0], vec![1], vec![0, 1]]
        );

        let mut output = AdditionalOutput::new();
        output.allocate_atom_to_bits();
        output.allocate_bit_paths();
        output.allocate_atom_counts();
        output.allocate_atoms_per_bit();
        output.reset_for_atom_count(molecule.num_atoms());
        for environment in &environments {
            environment.update_additional_output(&mut output, environment.bit_id());
        }
        environments[0].update_additional_output(&mut output, environments[0].bit_id());
        assert_eq!(output.atom_counts, Some(vec![3, 4, 2]));
        assert_eq!(
            output.atom_to_bits,
            Some(vec![
                vec![4_275_705_116, 1_524_090_560],
                vec![4_275_705_116, 4_274_652_475, 1_524_090_560],
                vec![4_274_652_475, 1_524_090_560],
            ])
        );
        assert_eq!(
            output.bit_paths,
            Some(BTreeMap::from([
                (1_524_090_560, vec![vec![0, 1]]),
                (4_274_652_475, vec![vec![1]]),
                (4_275_705_116, vec![vec![0], vec![0]]),
            ]))
        );
        assert_eq!(
            output.atoms_per_bit,
            Some(BTreeMap::from([
                (1_524_090_560, vec![vec![0, 1, 2]]),
                (4_274_652_475, vec![vec![1, 2]]),
                (4_275_705_116, vec![vec![0, 1], vec![0, 1]]),
            ]))
        );
    }

    #[test]
    fn fingerprint_smarts_parse_errors_are_returned() {
        let err = SsMatcher::from_pattern("[#6")
            .expect_err("invalid SMARTS must return an error instead of panicking");
        assert!(matches!(err, FingerprintError::InvalidSmartsPattern { .. }));
    }

    fn methane() -> Molecule {
        Molecule::from_smiles_with_sanitize("C", false).unwrap()
    }

    fn ethane() -> Molecule {
        Molecule::from_smiles_with_sanitize("CC", false).unwrap()
    }

    fn benzene() -> Molecule {
        Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap()
    }

    fn allocated_morgan_additional_output(num_atoms: usize) -> AdditionalOutput {
        let mut output = AdditionalOutput::new();
        output.allocate_atom_counts();
        output.allocate_atom_to_bits();
        output.allocate_bit_info_map();
        output.allocate_atoms_per_bit();
        output.reset_for_atom_count(num_atoms);
        output
    }

    fn morgan_arguments_for_env(
        radius: u32,
        include_redundant_environments: bool,
        include_chirality: bool,
        use_bond_types: bool,
    ) -> MorganArguments {
        MorganArguments::new(
            radius,
            false,
            include_chirality,
            false,
            vec![1, 2, 4, 8],
            0,
            include_redundant_environments,
            use_bond_types,
        )
        .unwrap()
    }

    fn rdkit_maccs_patterns_oracle() -> &'static [(usize, &'static str)] {
        &[
            (8usize, "[!#6!#1]1~*~*~*~1"),
            (11usize, "*1~*~*~*~1"),
            (13usize, "[#8]~[#7](~[#6])~[#6]"),
            (14usize, "[#16]-[#16]"),
            (15usize, "[#8]~[#6](~[#8])~[#8]"),
            (16usize, "[!#6!#1]1~*~*~1"),
            (17usize, "[#6]#[#6]"),
            (19usize, "*1~*~*~*~*~*~*~1"),
            (20usize, "[#14]"),
            (21usize, "[#6]=[#6](~[!#6!#1])~[!#6!#1]"),
            (22usize, "*1~*~*~1"),
            (23usize, "[#7]~[#6](~[#8])~[#8]"),
            (24usize, "[#7]-[#8]"),
            (25usize, "[#7]~[#6](~[#7])~[#7]"),
            (26usize, "[#6]=@[#6](@*)@*"),
            (28usize, "[!#6!#1]~[CH2]~[!#6!#1]"),
            (30usize, "[#6]~[!#6!#1](~[#6])(~[#6])~*"),
            (31usize, "[!#6!#1]~[F,Cl,Br,I]"),
            (32usize, "[#6]~[#16]~[#7]"),
            (33usize, "[#7]~[#16]"),
            (34usize, "[CH2]=*"),
            (36usize, "[#16R]"),
            (37usize, "[#7]~[#6](~[#8])~[#7]"),
            (38usize, "[#7]~[#6](~[#6])~[#7]"),
            (39usize, "[#8]~[#16](~[#8])~[#8]"),
            (40usize, "[#16]-[#8]"),
            (41usize, "[#6]#[#7]"),
            (43usize, "[!#6!#1!H0]~*~[!#6!#1!H0]"),
            (
                44usize,
                "[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]",
            ),
            (45usize, "[#6]=[#6]~[#7]"),
            (47usize, "[#16]~*~[#7]"),
            (48usize, "[#8]~[!#6!#1](~[#8])~[#8]"),
            (49usize, "[!+0]"),
            (50usize, "[#6]=[#6](~[#6])~[#6]"),
            (51usize, "[#6]~[#16]~[#8]"),
            (52usize, "[#7]~[#7]"),
            (53usize, "[!#6!#1!H0]~*~*~*~[!#6!#1!H0]"),
            (54usize, "[!#6!#1!H0]~*~*~[!#6!#1!H0]"),
            (55usize, "[#8]~[#16]~[#8]"),
            (56usize, "[#8]~[#7](~[#8])~[#6]"),
            (57usize, "[#8R]"),
            (58usize, "[!#6!#1]~[#16]~[!#6!#1]"),
            (59usize, "[#16]!:*:*"),
            (60usize, "[#16]=[#8]"),
            (61usize, "*~[#16](~*)~*"),
            (62usize, "*@*!@*@*"),
            (63usize, "[#7]=[#8]"),
            (64usize, "*@*!@[#16]"),
            (65usize, "c:n"),
            (66usize, "[#6]~[#6](~[#6])(~[#6])~*"),
            (67usize, "[!#6!#1]~[#16]"),
            (68usize, "[!#6!#1!H0]~[!#6!#1!H0]"),
            (69usize, "[!#6!#1]~[!#6!#1!H0]"),
            (70usize, "[!#6!#1]~[#7]~[!#6!#1]"),
            (71usize, "[#7]~[#8]"),
            (72usize, "[#8]~*~*~[#8]"),
            (73usize, "[#16]=*"),
            (74usize, "[CH3]~*~[CH3]"),
            (75usize, "*!@[#7]@*"),
            (76usize, "[#6]=[#6](~*)~*"),
            (77usize, "[#7]~*~[#7]"),
            (78usize, "[#6]=[#7]"),
            (79usize, "[#7]~*~*~[#7]"),
            (80usize, "[#7]~*~*~*~[#7]"),
            (81usize, "[#16]~*(~*)~*"),
            (82usize, "*~[CH2]~[!#6!#1!H0]"),
            (83usize, "[!#6!#1]1~*~*~*~*~1"),
            (84usize, "[NH2]"),
            (85usize, "[#6]~[#7](~[#6])~[#6]"),
            (86usize, "[C;H2,H3][!#6!#1][C;H2,H3]"),
            (87usize, "[F,Cl,Br,I]!@*@*"),
            (89usize, "[#8]~*~*~*~[#8]"),
            (
                90usize,
                "[$([!#6!#1!H0]~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[CH2R]1)]",
            ),
            (
                91usize,
                "[$([!#6!#1!H0]~*~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~*~[R]1@[R]@[CH2R]1)]",
            ),
            (92usize, "[#8]~[#6](~[#7])~[#6]"),
            (93usize, "[!#6!#1]~[CH3]"),
            (94usize, "[!#6!#1]~[#7]"),
            (95usize, "[#7]~*~*~[#8]"),
            (96usize, "*1~*~*~*~*~1"),
            (97usize, "[#7]~*~*~*~[#8]"),
            (98usize, "[!#6!#1]1~*~*~*~*~*~1"),
            (99usize, "[#6]=[#6]"),
            (100usize, "*~[CH2]~[#7]"),
            (
                101usize,
                "[$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1)]",
            ),
            (102usize, "[!#6!#1]~[#8]"),
            (104usize, "[!#6!#1!H0]~*~[CH2]~*"),
            (105usize, "*@*(@*)@*"),
            (106usize, "[!#6!#1]~*(~[!#6!#1])~[!#6!#1]"),
            (107usize, "[F,Cl,Br,I]~*(~*)~*"),
            (108usize, "[CH3]~*~*~*~[CH2]~*"),
            (109usize, "*~[CH2]~[#8]"),
            (110usize, "[#7]~[#6]~[#8]"),
            (111usize, "[#7]~*~[CH2]~*"),
            (112usize, "*~*(~*)(~*)~*"),
            (113usize, "[#8]!:*:*"),
            (114usize, "[CH3]~[CH2]~*"),
            (115usize, "[CH3]~*~[CH2]~*"),
            (116usize, "[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]"),
            (117usize, "[#7]~*~[#8]"),
            (118usize, "[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]"),
            (119usize, "[#7]=*"),
            (120usize, "[!#6R]"),
            (121usize, "[#7R]"),
            (122usize, "*~[#7](~*)~*"),
            (123usize, "[#8]~[#6]~[#8]"),
            (124usize, "[!#6!#1]~[!#6!#1]"),
            (126usize, "*!@[#8]!@*"),
            (127usize, "*@*!@[#8]"),
            (
                128usize,
                "[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2R]@[R]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2R]1),$(*~[CH2]~*~[R]1@[R]@[CH2R]1)]",
            ),
            (
                129usize,
                "[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[CH2R]1)]",
            ),
            (131usize, "[!#6!#1!H0]"),
            (132usize, "[#8]~*~[CH2]~*"),
            (133usize, "*@*!@[#7]"),
            (135usize, "[#7]!:*:*"),
            (136usize, "[#8]=*"),
            (137usize, "[!C!cR]"),
            (138usize, "[!#6!#1]~[CH2]~*"),
            (139usize, "[O!H0]"),
            (140usize, "[#8]"),
            (141usize, "[CH3]"),
            (142usize, "[#7]"),
            (144usize, "*!:*:*!:*"),
            (145usize, "*1~*~*~*~*~*~1"),
            (147usize, "[$(*~[CH2]~[CH2]~*),$([R]1@[CH2R]@[CH2R]1)]"),
            (148usize, "*~[!#6!#1](~*)~*"),
            (149usize, "[C;H3,H4]"),
            (150usize, "*!@*@*!@*"),
            (151usize, "[#7!H0]"),
            (152usize, "[#8]~[#6](~[#6])~[#6]"),
            (154usize, "[#6]=[#8]"),
            (155usize, "*!@[CH2]!@*"),
            (156usize, "[#7]~*(~*)~*"),
            (157usize, "[#6]-[#8]"),
            (158usize, "[#6]-[#7]"),
            (162usize, "a"),
            (165usize, "[R]"),
        ]
    }

    #[test]
    fn maccs_pattern_table_matches_rdkit_source_patterns() {
        assert_eq!(RDKIT_MACCS_RAW_BITS, 167);
        assert_eq!(COSMOLKIT_MACCS_PUBLIC_BITS, 166);
        assert_eq!(rdkit_maccs_public_index(0), None);
        assert_eq!(rdkit_maccs_public_index(1), Some(0));
        assert_eq!(rdkit_maccs_public_index(166), Some(165));
        assert_eq!(rdkit_maccs_public_index(167), None);

        let expected = rdkit_maccs_patterns_oracle();
        assert_eq!(RDKIT_MACCS_PATTERNS.len(), expected.len());
        assert_eq!(expected.len(), 136);

        for (actual, &(expected_bit, expected_smarts)) in
            RDKIT_MACCS_PATTERNS.iter().zip(expected.iter())
        {
            assert_eq!(actual.bit, expected_bit);
            assert_eq!(actual.smarts, expected_smarts);

            let looked_up =
                rdkit_maccs_pattern(expected_bit).expect("RDKit MACCS pattern bit is present");
            assert_eq!(looked_up.bit, expected_bit);
            assert_eq!(looked_up.smarts, expected_smarts);
        }

        for bit in 0..=RDKIT_MACCS_RAW_BITS {
            let expected_entry = expected
                .iter()
                .find(|&&(expected_bit, _)| expected_bit == bit);
            assert_eq!(
                rdkit_maccs_pattern(bit).map(|pattern| (pattern.bit, pattern.smarts)),
                expected_entry.copied(),
                "pattern lookup mismatch for RDKit MACCS bit {bit}"
            );
        }
    }

    fn single_atom_molecule(atomic_number: u8) -> Molecule {
        let mut builder = Molecule::builder();
        let element = Element::from_atomic_number(atomic_number).expect("valid element");
        builder.add_atom(AtomSpec::new(element));
        builder.build().unwrap()
    }

    fn assert_maccs_keys_001_040(smiles: &str, expected_public_bits: &[usize]) {
        let mol = Molecule::from_smiles(smiles).unwrap_or_else(|err| {
            panic!("MACCS keys 001-040 fixture {smiles} should parse: {err}")
        });
        assert_maccs_keys_001_040_for_mol(smiles, &mol, expected_public_bits);
    }

    fn assert_maccs_keys_001_040_for_mol(
        label: &str,
        mol: &Molecule,
        expected_public_bits: &[usize],
    ) {
        let params = MaccsFingerprintParams::default();
        let actual: Vec<usize> = maccs_fingerprint(mol, &params)
            .expect("MACCS fingerprint")
            .on_bits()
            .into_iter()
            .filter(|&bit| bit < 40)
            .collect();
        assert_eq!(
            actual, expected_public_bits,
            "RDKit MACCS keys 001-040 public projection mismatch for {label}"
        );
    }

    fn assert_maccs_keys_041_080(smiles: &str, expected_public_bits: &[usize]) {
        let mol = Molecule::from_smiles(smiles).unwrap_or_else(|err| {
            panic!("MACCS keys 041-080 fixture {smiles} should parse: {err}")
        });
        let params = MaccsFingerprintParams::default();
        let actual: Vec<usize> = maccs_fingerprint(&mol, &params)
            .expect("MACCS fingerprint")
            .on_bits()
            .into_iter()
            .filter(|&bit| (40..80).contains(&bit))
            .collect();
        assert_eq!(
            actual, expected_public_bits,
            "RDKit MACCS keys 041-080 public projection mismatch for {smiles}"
        );
    }

    fn assert_maccs_keys_081_120(smiles: &str, expected_public_bits: &[usize]) {
        let mol = Molecule::from_smiles(smiles).unwrap_or_else(|err| {
            panic!("MACCS keys 081-120 fixture {smiles} should parse: {err}")
        });
        let params = MaccsFingerprintParams::default();
        let actual: Vec<usize> = maccs_fingerprint(&mol, &params)
            .expect("MACCS fingerprint")
            .on_bits()
            .into_iter()
            .filter(|&bit| (80..120).contains(&bit))
            .collect();
        assert_eq!(
            actual, expected_public_bits,
            "RDKit MACCS keys 081-120 public projection mismatch for {smiles}"
        );
    }

    fn assert_maccs_keys_121_166(smiles: &str, expected_public_bits: &[usize]) {
        let mol = Molecule::from_smiles(smiles).unwrap_or_else(|err| {
            panic!("MACCS keys 121-166 fixture {smiles} should parse: {err}")
        });
        let params = MaccsFingerprintParams::default();
        let actual: Vec<usize> = maccs_fingerprint(&mol, &params)
            .expect("MACCS fingerprint")
            .on_bits()
            .into_iter()
            .filter(|&bit| (120..166).contains(&bit))
            .collect();
        assert_eq!(
            actual, expected_public_bits,
            "RDKit MACCS keys 121-166 public projection mismatch for {smiles}"
        );
    }

    fn assert_maccs_full_vector_for_mol(
        label: &str,
        mol: &Molecule,
        expected_raw_bits: &[usize],
        expected_public_bits: &[usize],
    ) {
        let raw = maccs_get_fingerprint_as_bit_vect(mol).expect("raw MACCS fingerprint");
        assert_eq!(raw.n_bits(), RDKIT_MACCS_RAW_BITS);
        assert!(
            !raw.on_bits().contains(&0),
            "RDKit MACCS raw bit 0 must stay unused for {label}"
        );
        assert_eq!(
            raw.on_bits(),
            expected_raw_bits,
            "RDKit MACCS raw 167-bit vector mismatch for {label}"
        );

        let params = MaccsFingerprintParams::default();
        let public = maccs_fingerprint(mol, &params).expect("MACCS fingerprint");
        assert_eq!(public.n_bits(), COSMOLKIT_MACCS_PUBLIC_BITS);
        assert_eq!(
            public.on_bits(),
            expected_public_bits,
            "COSMolKit MACCS public 166-bit projection mismatch for {label}"
        );

        let projected_from_raw: Vec<usize> = raw
            .on_bits()
            .into_iter()
            .filter_map(rdkit_maccs_public_index)
            .collect();
        assert_eq!(
            projected_from_raw, expected_public_bits,
            "RDKit raw-to-public MACCS projection mismatch for {label}"
        );
    }

    fn assert_maccs_full_vector(
        label: &str,
        smiles: &str,
        expected_raw_bits: &[usize],
        expected_public_bits: &[usize],
    ) {
        let mol = Molecule::from_smiles(smiles)
            .unwrap_or_else(|err| panic!("MACCS full-vector fixture {smiles} should parse: {err}"));
        assert_maccs_full_vector_for_mol(label, &mol, expected_raw_bits, expected_public_bits);
    }

    #[test]
    fn maccs_keys_001_040_direct_element_keys_match_rdkit() {
        let fixtures: &[(&str, u8, &[usize])] = &[
            ("key_002_rf", 104, &[1]),
            ("key_003_ge", 32, &[2]),
            ("key_004_ac", 89, &[3]),
            ("key_005_sc", 21, &[4]),
            ("key_006_la", 57, &[5]),
            ("key_007_v", 23, &[6]),
            ("key_009_fe", 26, &[8]),
            ("key_010_be", 4, &[9]),
            ("key_012_cu", 29, &[11]),
            ("key_018_boron", 5, &[17]),
            ("key_020_silicon", 14, &[19]),
            ("key_027_iodine", 53, &[26]),
            ("key_029_phosphorus", 15, &[28]),
            ("key_035_lithium", 3, &[34]),
        ];
        for &(label, atomic_number, expected_public_bits) in fixtures {
            let mol = single_atom_molecule(atomic_number);
            assert_maccs_keys_001_040_for_mol(label, &mol, expected_public_bits);
        }

        assert_maccs_keys_001_040("C", &[]);
    }

    #[test]
    fn maccs_fingerprint_full_bit_vectors_match_rdkit_raw_and_public_projection() {
        let empty = Molecule::new();
        assert_maccs_full_vector_for_mol("empty", &empty, &[], &[]);

        let fixtures: &[(&str, &str, &[usize], &[usize])] = &[
            ("methane", "C", &[160], &[159]),
            ("fluorine_atom", "F", &[42, 134], &[41, 133]),
            ("sulfur_atom", "S", &[88], &[87]),
            ("chlorine_atom", "Cl", &[103, 134], &[102, 133]),
            (
                "salt_fragments",
                "CCO.Cl",
                &[
                    82, 103, 109, 114, 131, 134, 139, 153, 155, 157, 160, 164, 166,
                ],
                &[
                    81, 102, 108, 113, 130, 133, 138, 152, 154, 156, 159, 163, 165,
                ],
            ),
            ("benzene", "c1ccccc1", &[162, 163, 165], &[161, 162, 164]),
            (
                "biphenyl",
                "c1ccccc1c2ccccc2",
                &[62, 125, 145, 162, 163, 165],
                &[61, 124, 144, 161, 162, 164],
            ),
            (
                "pyridine",
                "c1ncccc1",
                &[65, 98, 121, 137, 161, 162, 163, 165],
                &[64, 97, 120, 136, 160, 161, 162, 164],
            ),
            (
                "morpholine",
                "O1CCNCC1",
                &[
                    57, 82, 86, 91, 95, 98, 100, 104, 109, 111, 118, 120, 121, 128, 129, 132, 137,
                    138, 147, 151, 153, 157, 158, 161, 163, 164, 165,
                ],
                &[
                    56, 81, 85, 90, 94, 97, 99, 103, 108, 110, 117, 119, 120, 127, 128, 131, 136,
                    137, 146, 150, 152, 156, 157, 160, 162, 163, 164,
                ],
            ),
            ("ammonium", "[NH4+]", &[49, 151, 161], &[48, 150, 160]),
            (
                "acetate",
                "CC(=O)[O-]",
                &[49, 123, 154, 157, 159, 160, 164],
                &[48, 122, 153, 156, 158, 159, 163],
            ),
            ("isotopic_methane", "[13CH4]", &[160], &[159]),
            (
                "nitro",
                "N=O",
                &[63, 69, 71, 94, 102, 119, 124, 151, 161, 164],
                &[62, 68, 70, 93, 101, 118, 123, 150, 160, 163],
            ),
            (
                "cyclopropanol",
                "C1CC1O",
                &[22, 90, 104, 127, 132, 139, 143, 147, 152, 157, 164, 165],
                &[21, 89, 103, 126, 131, 138, 142, 146, 151, 156, 163, 164],
            ),
            (
                "fragment_methanes",
                "C.C",
                &[149, 160, 166],
                &[148, 159, 165],
            ),
            (
                "all_key_low_mix",
                "NCCO",
                &[
                    54, 82, 84, 95, 100, 104, 109, 111, 118, 131, 132, 138, 139, 147, 151, 153,
                    155, 157, 158, 161, 164,
                ],
                &[
                    53, 81, 83, 94, 99, 103, 108, 110, 117, 130, 131, 137, 138, 146, 150, 152, 154,
                    156, 157, 160, 163,
                ],
            ),
            (
                "all_key_high_mix",
                "OCOCOCO",
                &[
                    28, 82, 86, 89, 90, 109, 123, 126, 128, 131, 138, 139, 140, 146, 153, 155, 157,
                    159, 164,
                ],
                &[
                    27, 81, 85, 88, 89, 108, 122, 125, 127, 130, 137, 138, 139, 145, 152, 154, 156,
                    158, 163,
                ],
            ),
        ];

        for &(label, smiles, expected_raw_bits, expected_public_bits) in fixtures {
            assert_maccs_full_vector(label, smiles, expected_raw_bits, expected_public_bits);
        }

        let err = maccs_fingerprint(
            &Molecule::from_smiles("NCCO").unwrap(),
            &MaccsFingerprintParams { n_bits: 64 },
        )
        .unwrap_err();
        assert!(matches!(
            err,
            FingerprintError::UnsupportedOption {
                option: "MaccsFingerprintParams.n_bits",
                ..
            }
        ));
    }

    #[test]
    fn maccs_keys_001_040_pattern_keys_match_rdkit() {
        let fixtures: &[(&str, &str, &[usize])] = &[
            ("key_008_hetero_four_ring", "O1CCC1", &[7, 10]),
            ("key_011_four_ring", "C1CCC1", &[10]),
            ("key_013_o_n_c_c", "ON(C)C", &[12, 23]),
            ("key_014_disulfide", "CSSC", &[13]),
            ("key_015_o_c_o_o", "O=C(O)O", &[14]),
            ("key_016_hetero_three_ring", "O1CC1", &[15, 21]),
            ("key_017_alkyne", "C#C", &[16]),
            ("key_019_seven_ring", "C1CCCCCC1", &[18]),
            ("key_021_alkene_dihetero", "C=C(O)O", &[20, 33]),
            ("key_022_three_ring", "C1CC1", &[21]),
            ("key_023_n_c_o_o", "NC(=O)O", &[22]),
            ("key_024_n_o", "ON(C)C", &[12, 23]),
            ("key_025_n_c_n_n", "NC(N)N", &[24]),
            ("key_026_cyclic_alkene", "C1=C2CCCC2C1", &[10, 18, 25]),
            ("key_028_hetero_ch2_hetero", "OCO", &[27]),
            ("key_030_c_hetero_c_c_any", "C[S](C)(C)C", &[29]),
            ("key_031_hetero_halogen", "N[Pt](Cl)(Cl)N", &[8, 30]),
            ("key_032_c_s_n", "CSN", &[31, 32]),
            ("key_033_n_s", "CSN", &[31, 32]),
            ("key_034_ch2_double", "C=C", &[33]),
            ("key_036_s_ring", "S1CC1", &[15, 21, 35]),
            ("key_037_n_c_o_n", "NC(=O)N", &[36]),
            ("key_038_n_c_c_n", "NC(C)N", &[37]),
            ("key_039_o_s_o_o", "COS(=O)(=O)O", &[38, 39]),
            ("key_040_s_o", "CSO", &[39]),
        ];
        for &(label, smiles, expected_public_bits) in fixtures {
            assert_maccs_keys_001_040(smiles, expected_public_bits);
            assert!(
                expected_public_bits.iter().any(|&bit| bit + 1 <= 40),
                "{label} should exercise at least one RDKit raw key in 1..=40"
            );
        }
    }

    #[test]
    fn maccs_keys_041_080_match_rdkit() {
        let fixtures: &[(&str, &str, &[usize])] = &[
            ("key_041_c_n_triple", "C#N", &[40]),
            ("key_042_fluorine", "F", &[41]),
            ("key_043_hetero_bridge_h", "OCO", &[42]),
            ("key_044_exotic_element", "[SeH2]", &[43]),
            ("key_045_c_c_n", "C=CN", &[44]),
            ("key_046_bromine", "Br", &[45]),
            ("key_047_s_x_n", "SCN", &[42, 46]),
            (
                "key_048_o_hetero_o_o",
                "COS(=O)(=O)O",
                &[47, 54, 57, 59, 60, 66, 68, 72],
            ),
            ("key_049_charged", "[NH4+]", &[48]),
            ("key_050_substituted_alkene", "CC(C)=C", &[49, 73, 75]),
            ("key_051_c_s_o", "CSO", &[50, 66, 68]),
            ("key_052_n_n", "NNO", &[42, 51, 67, 68, 69, 70]),
            ("key_053_hetero_bridge_3", "NCCCO", &[52]),
            ("key_054_hetero_bridge_2", "NCCN", &[53, 78]),
            (
                "key_055_o_s_o",
                "CS(=O)(=O)C",
                &[50, 54, 57, 59, 60, 66, 72, 73],
            ),
            ("key_056_o_n_o_c", "ON(O)C", &[42, 55, 68, 69, 70]),
            ("key_057_o_ring", "O1CC1", &[56]),
            (
                "key_058_hetero_s_hetero",
                "CS(=O)(=O)C",
                &[50, 54, 57, 59, 60, 66, 72, 73],
            ),
            ("key_059_s_aromatic_chain", "Sc1ccccc1", &[58, 63]),
            ("key_060_s_o_double", "CS(=O)C", &[50, 59, 60, 66, 72, 73]),
            ("key_061_s_three_neighbors", "C[S](C)(C)C", &[60, 73]),
            ("key_062_ring_nonring_ring", "C1CC1C1CC1", &[61]),
            ("key_063_n_o_double", "N=O", &[62, 68, 70]),
            ("key_064_ring_to_s", "Sc1ccccc1", &[58, 63]),
            ("key_065_aromatic_c_n", "c1ncccc1", &[64]),
            ("key_066_quaternary_carbon", "CC(C)(C)C", &[65, 73]),
            ("key_067_hetero_s", "CSSC", &[66]),
            (
                "key_068_hetero_h_hetero_h",
                "NNO",
                &[42, 51, 67, 68, 69, 70],
            ),
            ("key_069_hetero_hetero_h", "CSN", &[66, 68]),
            ("key_070_hetero_n_hetero", "NNO", &[42, 51, 67, 68, 69, 70]),
            ("key_071_n_o", "N=O", &[62, 68, 70]),
            ("key_072_o_x_x_o", "OCCO", &[53, 71]),
            ("key_073_s_double", "CS(=O)C", &[50, 59, 60, 66, 72, 73]),
            ("key_074_methyl_bridge", "CCC", &[73]),
            ("key_075_exocyclic_n_ring", "CN1CC1", &[74]),
            ("key_076_substituted_alkene", "CC(C)=C", &[49, 73, 75]),
            ("key_077_n_x_n", "NC(=O)N", &[42, 76]),
            ("key_078_c_n_double", "C=N", &[77]),
            ("key_079_n_x_x_n", "NCCN", &[53, 78]),
            ("key_080_n_x_x_x_n", "NCCCN", &[52, 79]),
        ];
        for &(label, smiles, expected_public_bits) in fixtures {
            assert_maccs_keys_041_080(smiles, expected_public_bits);
            assert!(
                expected_public_bits
                    .iter()
                    .any(|&bit| (40..80).contains(&bit)),
                "{label} should exercise at least one RDKit raw key in 41..=80"
            );
        }
    }

    #[test]
    fn maccs_keys_081_120_match_rdkit() {
        let fixtures: &[(&str, &str, &[usize])] = &[
            (
                "key_081_s_three_neighbors",
                "C[S](C)(C)C",
                &[85, 87, 92, 111],
            ),
            ("key_082_ch2_hetero_h", "CCN", &[81, 83, 99, 113]),
            (
                "key_083_hetero_five_ring",
                "O1CCCC1",
                &[82, 85, 95, 108, 117],
            ),
            ("key_084_nh2", "N", &[]),
            ("key_085_c_n_c_c", "CN(C)C", &[84, 85, 92]),
            ("key_086_c_h2_h3_hetero_c", "CN(C)C", &[84, 85, 92]),
            ("key_087_halogen_ring_chain", "FC1CC1", &[86, 106]),
            ("key_088_sulfur", "S", &[87]),
            ("key_089_o_bridge_3", "OCCCO", &[81, 88, 89, 103, 108, 117]),
            (
                "key_090_hetero_ch2_bridge",
                "NCCO",
                &[81, 83, 94, 99, 103, 108, 110, 117],
            ),
            (
                "key_091_hetero_ch2_bridge_4",
                "NCCCO",
                &[81, 83, 89, 96, 99, 103, 108, 110, 117],
            ),
            ("key_092_o_c_n_c", "OC(N)C", &[83, 91, 109, 116]),
            ("key_093_hetero_methyl", "CN", &[83, 92]),
            ("key_094_hetero_n", "CN", &[83, 92]),
            (
                "key_095_n_bridge_o",
                "NCCO",
                &[81, 83, 94, 99, 103, 108, 110, 117],
            ),
            ("key_096_five_ring", "C1CCCC1", &[95, 117]),
            (
                "key_097_n_bridge3_o",
                "NCCCO",
                &[81, 83, 89, 96, 99, 103, 108, 110, 117],
            ),
            ("key_098_hetero_six_ring", "O1CCCCC1", &[85, 97, 108, 117]),
            ("key_099_alkene", "C=C", &[98]),
            ("key_100_ch2_n", "CCN", &[81, 83, 99, 113]),
            ("key_101_large_ring", "C1CCCCCCC1", &[100, 117]),
            ("key_102_hetero_o", "CO", &[92]),
            ("key_103_chlorine", "Cl", &[102]),
            ("key_104_hetero_ch2_chain", "NCC", &[81, 83, 99, 113]),
            ("key_105_ring_branch_ring", "C1CC(C1)C1CC1", &[117]),
            (
                "key_106_hetero_three_neighbors",
                "N(O)(O)O",
                &[93, 101, 105],
            ),
            (
                "key_107_halogen_three_neighbors",
                "FC(F)(F)F",
                &[105, 106, 111],
            ),
            ("key_108_methyl_chain_ch2", "CCCC", &[113, 114, 117]),
            ("key_109_ch2_o", "CCO", &[81, 108, 113]),
            ("key_110_n_c_o", "NCO", &[81, 83, 99, 108, 109, 116]),
            ("key_111_n_ch2_chain", "NCC", &[81, 83, 99, 113]),
            ("key_112_quaternary_any", "CC(C)(C)C", &[111]),
            ("key_113_o_aromatic_chain", "Oc1ccccc1", &[112]),
            ("key_114_methyl_ch2_any", "CCC", &[113]),
            ("key_115_methyl_any_ch2_any", "CC(C)C", &[]),
            ("key_116_methyl_ch2_bridge", "CCCC", &[113, 114, 117]),
            ("key_117_n_x_o", "NCO", &[81, 83, 99, 108, 109, 116]),
            ("key_118_two_ch2_paths", "CCCC", &[113, 114, 117]),
            ("key_119_n_double_any", "N=O", &[93, 101, 118]),
            (
                "key_120_two_noncarbon_ring_atoms",
                "N1CCO1",
                &[81, 89, 93, 94, 99, 101, 103, 108, 110, 117, 119],
            ),
        ];
        for &(label, smiles, expected_public_bits) in fixtures {
            assert_maccs_keys_081_120(smiles, expected_public_bits);
            if !expected_public_bits.is_empty() {
                assert!(
                    expected_public_bits
                        .iter()
                        .any(|&bit| (80..120).contains(&bit)),
                    "{label} should exercise at least one RDKit raw key in 81..=120"
                );
            }
        }
    }

    #[test]
    fn maccs_keys_121_166_match_rdkit() {
        let fixtures: &[(&str, &str, &[usize])] = &[
            (
                "key_121_n_ring",
                "N1CCCC1",
                &[120, 128, 136, 137, 146, 150, 152, 157, 160, 164],
            ),
            (
                "key_122_n_three_neighbors",
                "CN(C)C",
                &[121, 140, 147, 148, 157, 159, 160],
            ),
            ("key_123_o_c_o", "COC", &[125, 148, 156, 159, 163]),
            ("key_124_two_hetero", "NO", &[123, 130, 138, 150, 160, 163]),
            (
                "key_125_two_aromatic_rings",
                "c1ccccc1c2ccccc2",
                &[124, 144, 161, 162, 164],
            ),
            ("key_126_o_nonring", "COC", &[125, 148, 156, 159, 163]),
            (
                "key_127_ring_nonring_o",
                "C1CC1O",
                &[126, 131, 138, 142, 146, 151, 156, 163, 164],
            ),
            (
                "key_128_ch2_bridge_5",
                "CCCCCCC",
                &[127, 128, 146, 148, 154, 159],
            ),
            ("key_129_ch2_bridge_4", "CCCCCC", &[128, 146, 148, 154, 159]),
            (
                "key_130_two_hetero_pairs",
                "NOON",
                &[123, 125, 129, 130, 141, 150, 158, 160, 163],
            ),
            (
                "key_131_two_hetero_h",
                "NNO",
                &[123, 129, 130, 138, 141, 150, 160, 163],
            ),
            (
                "key_132_o_ch2_chain",
                "OCCC",
                &[131, 138, 146, 152, 154, 156, 159, 163],
            ),
            (
                "key_133_ring_nonring_n",
                "C1CC1N",
                &[132, 146, 150, 155, 157, 160, 164],
            ),
            ("key_134_halogen", "F", &[133]),
            (
                "key_135_n_aromatic_chain",
                "Nc1ccccc1",
                &[132, 134, 150, 155, 157, 160, 161, 162, 164],
            ),
            ("key_136_two_o_double", "O=CC=O", &[135, 153, 158, 163]),
            (
                "key_137_noncarbon_ring",
                "O1CC1",
                &[136, 146, 152, 156, 163, 164],
            ),
            (
                "key_138_two_hetero_ch2",
                "NCCO",
                &[130, 131, 137, 138, 146, 150, 152, 154, 156, 157, 160, 163],
            ),
            ("key_139_o_no_h", "COC", &[125, 148, 156, 159, 163]),
            (
                "key_140_four_oxygens",
                "OCOCOCO",
                &[
                    122, 125, 127, 130, 137, 138, 139, 145, 152, 154, 156, 158, 163,
                ],
            ),
            ("key_141_three_methyl", "CC(C)(C)C", &[140, 148, 159]),
            ("key_142_two_nitrogens", "NN", &[123, 130, 141, 150, 160]),
            (
                "key_143_ring_nonring_o_once",
                "C1CC1O",
                &[126, 131, 138, 142, 146, 151, 156, 163, 164],
            ),
            ("key_144_aromatic_chain_three", "c1ccccc1", &[161, 162, 164]),
            (
                "key_145_two_six_rings",
                "C1CCCCC1C1CCCCC1",
                &[127, 128, 144, 146, 162, 164],
            ),
            (
                "key_146_three_oxygens",
                "OCOCO",
                &[122, 125, 130, 137, 138, 145, 152, 154, 156, 158, 163],
            ),
            ("key_147_two_ch2", "CCCC", &[146, 148, 154, 159]),
            (
                "key_148_hetero_three_neighbors",
                "N(C)(C)C",
                &[121, 140, 147, 148, 157, 159, 160],
            ),
            ("key_149_two_methyl", "CCC", &[148, 154, 159]),
            (
                "key_150_nonring_ring_path",
                "C1CC1CC1CC1",
                &[127, 128, 146, 154, 164],
            ),
            (
                "key_151_n_no_h",
                "CN(C)C",
                &[121, 140, 147, 148, 157, 159, 160],
            ),
            ("key_152_o_c_c_c", "OC(C)C", &[138, 148, 151, 156, 159, 163]),
            (
                "key_153_hetero_ch2_once",
                "CCN",
                &[150, 152, 154, 157, 159, 160],
            ),
            ("key_154_carbonyl", "C=O", &[153, 163]),
            ("key_155_ch2_nonring", "CCC", &[148, 154, 159]),
            (
                "key_156_n_three_neighbors",
                "N(C)(C)C",
                &[121, 140, 147, 148, 157, 159, 160],
            ),
            ("key_157_c_o_single", "CO", &[138, 156, 159, 163]),
            ("key_158_c_n_single", "CN", &[150, 157, 159, 160]),
            (
                "key_159_two_oxygens",
                "OCO",
                &[122, 130, 138, 152, 154, 156, 158, 163],
            ),
            ("key_160_methyl_once", "CC", &[148, 159]),
            ("key_161_n_once", "CN", &[150, 157, 159, 160]),
            ("key_162_aromatic", "c1ccccc1", &[161, 162, 164]),
            (
                "key_163_six_ring_once",
                "C1CCCCC1",
                &[127, 128, 146, 162, 164],
            ),
            ("key_164_o_once", "CO", &[138, 156, 159, 163]),
            ("key_165_ring", "C1CC1", &[146, 164]),
            ("key_166_fragments", "C.C", &[148, 159, 165]),
        ];
        for &(label, smiles, expected_public_bits) in fixtures {
            assert_maccs_keys_121_166(smiles, expected_public_bits);
            assert!(
                expected_public_bits
                    .iter()
                    .any(|&bit| (120..166).contains(&bit)),
                "{label} should exercise at least one RDKit raw key in 121..=166"
            );
        }
    }

    fn rdkit_morgan_env_output(
        mol: &Molecule,
        arguments: &MorganArguments,
        from_atoms: Option<&[usize]>,
        ignore_atoms: Option<&[usize]>,
    ) -> Result<MorganAdditionalOutput, FingerprintError> {
        let atom_invariants = MorganAtomInvGenerator::new(true).getAtomInvariants(mol);
        let bond_invariants = MorganBondInvGenerator::new(arguments.df_use_bond_types, false)
            .getBondInvariants(mol)?;
        let envs = MorganEnvGenerator::new().getEnvironments(
            mol,
            arguments,
            from_atoms,
            ignore_atoms,
            &atom_invariants,
            &bond_invariants,
        )?;
        let mut output = allocated_morgan_additional_output(mol.num_atoms());
        for env in envs {
            let bit_id = env.getBitId();
            env.updateAdditionalOutput(&mut output, bit_id);
        }
        Ok(morgan_additional_output_from_rdkit_output(output))
    }

    fn bimap(entries: &[(usize, &[(usize, u32)])]) -> BTreeMap<usize, Vec<(usize, u32)>> {
        entries
            .iter()
            .map(|(bit, values)| (*bit, values.to_vec()))
            .collect()
    }

    fn atoms_per_bit(entries: &[(usize, &[&[usize]])]) -> BTreeMap<usize, Vec<Vec<usize>>> {
        entries
            .iter()
            .map(|(bit, values)| (*bit, values.iter().map(|atoms| atoms.to_vec()).collect()))
            .collect()
    }

    fn rdkit_morgan_oracle_mol(smiles: &str) -> Molecule {
        Molecule::from_smiles(smiles).unwrap()
    }

    fn two_atom_molecule_with_bond(order: BondOrder) -> Molecule {
        two_atom_molecule_with_bond_spec(BondSpec::new(AtomId::new(0), AtomId::new(1), order))
    }

    fn two_atom_molecule_with_bond_spec(spec: BondSpec) -> Molecule {
        let mut builder = Molecule::builder();
        let left = builder.add_atom(AtomSpec::new(Element::C));
        let right = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(left, right, spec.order()).with_stereo(spec.stereo()))
            .unwrap();
        builder.build().unwrap()
    }

    fn double_bond_stereo_molecule(stereo: BondStereo) -> Molecule {
        let mut builder = Molecule::builder();
        let left = builder.add_atom(AtomSpec::new(Element::C));
        let right = builder.add_atom(AtomSpec::new(Element::C));
        let left_ref = builder.add_atom(AtomSpec::new(Element::C));
        let right_ref = builder.add_atom(AtomSpec::new(Element::C));
        let mut double_bond = BondSpec::new(left, right, BondOrder::Double).with_stereo(stereo);
        if matches!(stereo, BondStereo::Cis | BondStereo::Trans) {
            double_bond = double_bond.with_stereo_atoms(left_ref, right_ref);
        }
        builder.add_bond(double_bond).unwrap();
        builder
            .add_bond(BondSpec::new(left, left_ref, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(right, right_ref, BondOrder::Single))
            .unwrap();
        builder.build().unwrap()
    }

    fn stereo_done_double_bond_without_cip_computed(stereo: BondStereo) -> Molecule {
        double_bond_stereo_molecule(stereo).with_prop("_StereochemDone", "1")
    }

    #[test]
    fn morgan_fingerprint_empty_molecule_returns_empty_fingerprint() {
        let mol = Molecule::from_smiles_with_sanitize("", false).unwrap();
        let fp = morgan_fingerprint(&mol, &default_morgan_params(2, 2048)).unwrap();
        assert_eq!(fp.on_bits(), Vec::<usize>::new());
    }

    #[test]
    fn morgan_fingerprint_empty_params_n_bits_zero_returns_error() {
        let mol = methane();
        let params = MorganFingerprintParams {
            n_bits: 0,
            ..Default::default()
        };
        assert!(matches!(
            morgan_fingerprint(&mol, &params),
            Err(FingerprintError::EmptyFingerprint)
        ));
    }

    #[test]
    fn morgan_fingerprint_methane_radius0_produces_deterministic_fingerprint() {
        let mol = methane();
        let fp_a = morgan_fingerprint(&mol, &default_morgan_params(0, 2048)).unwrap();
        let fp_b = morgan_fingerprint(&mol, &default_morgan_params(0, 2048)).unwrap();
        assert_eq!(fp_a, fp_b);
        assert!(!fp_a.on_bits().is_empty(), "expected at least one on-bit");
    }

    #[test]
    fn morgan_fingerprint_tanimoto_self_is_one() {
        let mol = benzene();
        let fp = morgan_fingerprint(&mol, &default_morgan_params(2, 2048)).unwrap();
        let similarity = fp.tanimoto(&fp).unwrap();
        assert!(
            (similarity - 1.0).abs() < 1e-9,
            "tanimoto of fingerprint with itself should be 1.0, got {similarity}"
        );
    }

    #[test]
    fn morgan_fingerprint_n_bits_matches_param() {
        let mol = benzene();
        for n_bits in [64, 256, 1024, 2048] {
            let fp = morgan_fingerprint(&mol, &default_morgan_params(2, n_bits)).unwrap();
            assert_eq!(fp.n_bits(), n_bits);
        }
    }

    #[test]
    fn morgan_fingerprint_radius_determinism() {
        let mol = benzene();
        for radius in 0..=3 {
            let fp_a = morgan_fingerprint(&mol, &default_morgan_params(radius, 2048)).unwrap();
            let fp_b = morgan_fingerprint(&mol, &default_morgan_params(radius, 2048)).unwrap();
            assert_eq!(fp_a, fp_b, "radius={radius} should be deterministic");
        }
    }

    #[test]
    fn morgan_fingerprint_ethane_and_methane_differ() {
        let m = methane();
        let e = ethane();
        let fp_m = morgan_fingerprint(&m, &default_morgan_params(0, 2048)).unwrap();
        let fp_e = morgan_fingerprint(&e, &default_morgan_params(0, 2048)).unwrap();
        assert_ne!(
            fp_m, fp_e,
            "methane and ethane should have different fingerprints"
        );
    }

    #[test]
    fn morgan_fingerprint_benzene_and_ethane_differ() {
        let b = benzene();
        let e = ethane();
        let fp_b = morgan_fingerprint(&b, &default_morgan_params(2, 2048)).unwrap();
        let fp_e = morgan_fingerprint(&e, &default_morgan_params(2, 2048)).unwrap();
        assert_ne!(
            fp_b, fp_e,
            "benzene and ethane should have different fingerprints"
        );
    }

    #[test]
    fn morgan_fingerprint_radius_increases_on_bits() {
        let mol = ethane();
        let fp_r0 = morgan_fingerprint(&mol, &default_morgan_params(0, 2048)).unwrap();
        let fp_r2 = morgan_fingerprint(&mol, &default_morgan_params(2, 2048)).unwrap();
        assert!(
            fp_r2.on_bits().len() >= fp_r0.on_bits().len(),
            "larger radius should produce at least as many on-bits"
        );
    }

    #[test]
    fn morgan_fingerprint_with_output_produces_additional_data() {
        let mol = ethane();
        let params = MorganFingerprintParams {
            radius: 1,
            n_bits: 2048,
            collect_additional_output: true,
            ..Default::default()
        };
        let output = morgan_fingerprint_with_output(&mol, &params).unwrap();
        assert!(output.additional_output.is_some());
        let extra = output.additional_output.unwrap();
        assert_eq!(extra.atom_counts.len(), mol.num_atoms());
        assert!(!extra.bit_info_map.is_empty());
    }

    #[test]
    fn morgan_fingerprint_from_atoms_filters_by_allowed_indices() {
        let mol = ethane();
        let params = MorganFingerprintParams {
            radius: 0,
            n_bits: 2048,
            from_atoms: Some(vec![0]),
            ..Default::default()
        };
        let fp = morgan_fingerprint(&mol, &params).unwrap();
        assert!(!fp.on_bits().is_empty());

        let params_empty = MorganFingerprintParams {
            radius: 0,
            n_bits: 2048,
            from_atoms: Some(vec![]),
            ..Default::default()
        };
        let fp_empty = morgan_fingerprint(&mol, &params_empty).unwrap();
        assert!(
            fp_empty.on_bits().is_empty(),
            "no from_atoms → empty fingerprint"
        );
    }

    #[test]
    fn morgan_fingerprint_ignore_atoms_excludes_indices() {
        let mol = ethane();
        let params_full = default_morgan_params(0, 2048);
        let params_exclude = MorganFingerprintParams {
            radius: 0,
            n_bits: 2048,
            ignore_atoms: Some(vec![1]),
            ..Default::default()
        };
        let fp_full = morgan_fingerprint(&mol, &params_full).unwrap();
        let fp_excluded = morgan_fingerprint(&mol, &params_exclude).unwrap();
        assert_ne!(fp_full.on_bits().len(), 0);
        assert!(
            fp_excluded.on_bits().len() <= fp_full.on_bits().len(),
            "excluding an atom should not increase on-bits"
        );
    }

    #[test]
    fn morgan_fingerprint_feature_generator_produces_deterministic_fingerprint() {
        let mol = benzene();
        let params = MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            atom_invariants_generator: MorganAtomInvariantsGenerator::Feature,
            ..Default::default()
        };
        // Feature invariants now use element/property classification.
        let fp_a = morgan_fingerprint(&mol, &params).unwrap();
        let fp_b = morgan_fingerprint(&mol, &params).unwrap();
        assert_eq!(fp_a, fp_b, "feature invariants should be deterministic");
        assert!(
            !fp_a.on_bits().is_empty(),
            "expected on-bits from feature invariants"
        );
    }

    #[test]
    fn morgan_fingerprint_custom_invariants_override_default() {
        let mol = ethane();
        let custom = vec![42u32, 99u32];
        let params = MorganFingerprintParams {
            radius: 0,
            n_bits: 2048,
            custom_atom_invariants: Some(custom),
            ..Default::default()
        };
        let fp_a = morgan_fingerprint(&mol, &params).unwrap();
        let fp_b = morgan_fingerprint(&mol, &params).unwrap();
        assert_eq!(fp_a, fp_b);
    }

    #[test]
    fn morgan_fingerprint_zero_bonds_molecule_does_not_panic() {
        let mol = Molecule::from_smiles_with_sanitize("[H][H]", false).unwrap();
        let fp = morgan_fingerprint(&mol, &default_morgan_params(2, 2048));
        assert!(fp.is_ok());
    }

    #[test]
    fn morgan_fingerprint_count_simulation_runs() {
        let mol = benzene();
        let params = MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            count_simulation: true,
            count_bounds: vec![1, 2, 4, 8],
            ..Default::default()
        };
        let fp = morgan_fingerprint(&mol, &params).unwrap();
        assert!(!fp.on_bits().is_empty());
        let std_fp = morgan_fingerprint(&mol, &default_morgan_params(2, 2048)).unwrap();
        assert!(
            fp.on_bits().len() >= std_fp.on_bits().len(),
            "count-simulation should set at least as many bits as standard mode"
        );
    }

    #[test]
    fn morgan_fingerprint_uses_topology_adjacency_without_derived_cache() {
        let mol = benzene();
        let fp = morgan_fingerprint(&mol, &default_morgan_params(2, 2048)).unwrap();
        assert!(!fp.on_bits().is_empty());
    }

    #[test]
    fn generator_chirality_preparation_morgan_clones_unprepared_input_without_regression() {
        let prepared = Molecule::from_smiles("C[C@@H](O)CC").unwrap();
        assert_eq!(prepared.prop("_StereochemDone"), Some("1"));
        let mut unprepared = prepared.clone();
        unprepared.properties_mut().clear_prop("_StereochemDone");
        let source_snapshot = unprepared.clone();
        let params = MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            use_chirality: true,
            ..Default::default()
        };

        let prepared_fingerprint = morgan_fingerprint(&prepared, &params).unwrap();
        let unprepared_fingerprint = morgan_fingerprint(&unprepared, &params).unwrap();

        assert_eq!(unprepared_fingerprint, prepared_fingerprint);
        assert_eq!(unprepared, source_snapshot);
        assert_eq!(unprepared.prop("_StereochemDone"), None);
    }

    #[test]
    fn generator_chirality_preparation_atom_pair_preserves_r_s_and_no_label_states() {
        use atom_pair::{AtomPairArguments, atom_pair_generator};

        let r_molecule = Molecule::from_smiles("C[C@@H](O)CC").unwrap();
        let s_molecule = Molecule::from_smiles("C[C@H](O)CC").unwrap();
        let r_with_cip_computed = assign_cip_labels(&r_molecule, 0).unwrap();
        assert_eq!(r_with_cip_computed.prop("_CIPComputed"), Some("1"));
        let mut no_label = r_molecule.clone();
        no_label.properties_mut().clear_prop("_CIPComputed");
        for atom in &mut no_label.topology_block_mut().atoms {
            atom.clear_prop("_CIPCode");
            atom.clear_prop("_CIPRank");
        }
        let mut unprepared_no_label = no_label.clone();
        unprepared_no_label
            .properties_mut()
            .clear_prop("_StereochemDone");
        let source_snapshot = unprepared_no_label.clone();

        let arguments = AtomPairArguments::new(false, true, true, 1, 30, Vec::new(), 2048)
            .expect("valid chiral AtomPair arguments");
        let generator = atom_pair_generator(&arguments, None, false);
        let sparse_count = |molecule: &Molecule| {
            generator
                .sparse_count_fingerprint(molecule, &mut FingerprintFuncArguments::default())
                .unwrap()
        };

        let r_fingerprint = sparse_count(&r_molecule);
        let r_with_cip_computed_fingerprint = sparse_count(&r_with_cip_computed);
        let s_fingerprint = sparse_count(&s_molecule);
        let no_label_fingerprint = sparse_count(&no_label);
        let unprepared_no_label_fingerprint = sparse_count(&unprepared_no_label);

        assert_ne!(r_fingerprint, s_fingerprint);
        assert_eq!(r_fingerprint, r_with_cip_computed_fingerprint);
        assert_ne!(r_fingerprint, no_label_fingerprint);
        assert_eq!(unprepared_no_label_fingerprint, no_label_fingerprint);
        assert_eq!(unprepared_no_label, source_snapshot);
        assert_eq!(unprepared_no_label.prop("_StereochemDone"), None);
        assert_eq!(unprepared_no_label.prop("_CIPComputed"), None);
    }

    #[test]
    fn generator_chirality_preparation_custom_atom_invariants_remain_authoritative() {
        use atom_pair::{AtomPairArguments, atom_pair_generator};

        let r_molecule = Molecule::from_smiles("C[C@@H](O)CC").unwrap();
        let mut s_molecule = Molecule::from_smiles("C[C@H](O)CC").unwrap();
        s_molecule.properties_mut().clear_prop("_StereochemDone");
        let source_snapshot = s_molecule.clone();
        let arguments = AtomPairArguments::new(false, true, true, 1, 30, Vec::new(), 2048)
            .expect("valid chiral AtomPair arguments");
        let generator = atom_pair_generator(&arguments, None, false);
        let custom_atom_invariants = vec![17; r_molecule.num_atoms()];
        let fingerprint = |molecule: &Molecule| {
            generator
                .sparse_count_fingerprint(
                    molecule,
                    &mut FingerprintFuncArguments {
                        custom_atom_invariants: Some(custom_atom_invariants.clone()),
                        ..Default::default()
                    },
                )
                .unwrap()
        };

        assert_eq!(fingerprint(&r_molecule), fingerprint(&s_molecule));
        assert_eq!(s_molecule, source_snapshot);
    }

    #[test]
    fn morgan_fingerprint_chirality_disabled_produces_identical_fingerprints() {
        // Without chirality, R and S enantiomers produce the same fingerprint.
        let r_mol = Molecule::from_smiles_with_sanitize("C[C@@H](O)CC", false).unwrap();
        let s_mol = Molecule::from_smiles_with_sanitize("C[C@H](O)CC", false).unwrap();
        let params = MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            use_chirality: false,
            ..Default::default()
        };
        let fp_r = morgan_fingerprint(&r_mol, &params).unwrap();
        let fp_s = morgan_fingerprint(&s_mol, &params).unwrap();
        assert_eq!(
            fp_r, fp_s,
            "R and S should have same fingerprint without chirality"
        );
    }

    #[test]
    fn morgan_fingerprint_custom_bond_invariants_override_default() {
        let mol = ethane();
        let params = MorganFingerprintParams {
            radius: 0,
            n_bits: 2048,
            custom_bond_invariants: Some(vec![5u32]),
            ..Default::default()
        };
        let fp = morgan_fingerprint(&mol, &params).unwrap();
        assert!(!fp.on_bits().is_empty());
    }

    #[test]
    fn additional_output_allocates_empty_owned_containers() {
        let mut ao = AdditionalOutput::new();

        assert_eq!(ao, AdditionalOutput::default());

        ao.allocate_atom_to_bits();
        ao.allocate_bit_info_map();
        ao.allocate_bit_paths();
        ao.allocate_atom_counts();
        ao.allocate_atoms_per_bit();

        assert_eq!(ao.atom_to_bits, Some(Vec::new()));
        assert_eq!(ao.bit_info_map, Some(BTreeMap::new()));
        assert_eq!(ao.bit_paths, Some(BTreeMap::new()));
        assert_eq!(ao.atom_counts, Some(Vec::new()));
        assert_eq!(ao.atoms_per_bit, Some(BTreeMap::new()));
    }

    #[test]
    fn additional_output_reset_reinitializes_allocated_containers() {
        let mut ao = AdditionalOutput::new();
        ao.allocate_atom_to_bits();
        ao.allocate_bit_info_map();
        ao.allocate_bit_paths();
        ao.allocate_atom_counts();
        ao.allocate_atoms_per_bit();

        if let Some(atom_to_bits) = ao.atom_to_bits.as_mut() {
            atom_to_bits.push(vec![1, 2, 3]);
        }
        if let Some(bit_info_map) = ao.bit_info_map.as_mut() {
            bit_info_map.insert(7, vec![(4, 2)]);
        }
        if let Some(bit_paths) = ao.bit_paths.as_mut() {
            bit_paths.insert(11, vec![vec![0, 1]]);
        }
        if let Some(atom_counts) = ao.atom_counts.as_mut() {
            atom_counts.extend([3, 4]);
        }
        if let Some(atoms_per_bit) = ao.atoms_per_bit.as_mut() {
            atoms_per_bit.insert(13, vec![vec![2, 5]]);
        }

        ao.reset_for_atom_count(4);

        assert_eq!(
            ao.atom_to_bits,
            Some(vec![Vec::new(), Vec::new(), Vec::new(), Vec::new()])
        );
        assert_eq!(ao.bit_info_map, Some(BTreeMap::new()));
        assert_eq!(ao.bit_paths, Some(BTreeMap::new()));
        assert_eq!(ao.atom_counts, Some(vec![0, 0, 0, 0]));
        assert_eq!(
            ao.atoms_per_bit,
            Some(BTreeMap::from([(13, vec![vec![2, 5]])]))
        );
    }

    #[test]
    fn additional_output_reset_leaves_unallocated_fields_unset() {
        let mut ao = AdditionalOutput::new();
        ao.reset_for_atom_count(3);
        assert_eq!(ao, AdditionalOutput::default());
    }

    #[test]
    fn morgan_atom_env_get_bit_id_returns_code_identity() {
        let mol = ethane();
        let env = MorganAtomEnv::new(0x1234_5678_9abc_def0, 1, 0, &mol);

        assert_eq!(env.getBitId(), 0x1234_5678_9abc_def0);
    }

    #[test]
    fn morgan_atom_env_update_additional_output_updates_allocated_rdkit_fields() {
        let mol = ethane();
        let mut output = allocated_morgan_additional_output(mol.num_atoms());
        let env = MorganAtomEnv::new(0x55, 1, 0, &mol);

        env.updateAdditionalOutput(&mut output, 77);
        env.updateAdditionalOutput(&mut output, 91);

        assert_eq!(output.atom_counts, Some(vec![0, 2]));
        assert_eq!(output.atom_to_bits, Some(vec![Vec::new(), vec![77, 91]]));
        assert_eq!(
            output.bit_info_map,
            Some(BTreeMap::from([(77, vec![(1, 0)]), (91, vec![(1, 0)])]))
        );
        assert_eq!(
            output.atoms_per_bit,
            Some(BTreeMap::from([(77, vec![vec![1]]), (91, vec![vec![1]])]))
        );
        assert_eq!(output.bit_paths, None);
    }

    #[test]
    fn morgan_atom_env_atoms_per_bit_uses_center_then_atom_index_distance_order() {
        let mol = Molecule::from_smiles_with_sanitize("CC(C)O", false).unwrap();
        let mut output = allocated_morgan_additional_output(mol.num_atoms());
        let env = MorganAtomEnv::new(0x99, 1, 1, &mol);

        env.updateAdditionalOutput(&mut output, 123);

        assert_eq!(
            output.atoms_per_bit,
            Some(BTreeMap::from([(123, vec![vec![1, 0, 2, 3]])]))
        );
        assert_eq!(
            output.bit_info_map,
            Some(BTreeMap::from([(123, vec![(1, 1)])]))
        );
        assert_eq!(output.atom_counts, Some(vec![0, 1, 0, 0]));
        assert_eq!(
            output.atom_to_bits,
            Some(vec![Vec::new(), vec![123], Vec::new(), Vec::new()])
        );
    }

    #[test]
    fn morgan_env_generator_metadata_json_and_result_size_match_rdkit_shape() {
        let mut generator = MorganEnvGenerator::new();

        assert_eq!(generator.infoString(), "MorganEnvironmentGenerator");
        assert_eq!(generator.toJSON(), r#"{"type":"MorganEnvGenerator"}"#);
        assert_eq!(generator.getResultSize(), u64::from(u32::MAX));

        generator.fromJSON("").unwrap();
        generator
            .fromJSON(r#"{"type":"MorganEnvGenerator"}"#)
            .unwrap();
        generator.fromJSON(r#"{"type":"OtherGenerator"}"#).unwrap();
        assert_eq!(generator.infoString(), "MorganEnvironmentGenerator");
        assert_eq!(generator.toJSON(), r#"{"type":"MorganEnvGenerator"}"#);
        assert!(matches!(
            generator.fromJSON("[1,2,3]"),
            Err(FingerprintError::InvalidArgumentsJson(_))
        ));
        assert!(matches!(
            generator.fromJSON("{"),
            Err(FingerprintError::InvalidArgumentsJson(_))
        ));
    }

    #[test]
    fn get_morgan_generator_defaults_allocate_owned_invariant_generators() {
        let args =
            MorganArguments::new(2, false, true, true, vec![1, 2, 4], 1024, false, false).unwrap();

        let generator = getMorganGenerator(&args, None, None, false, false);

        assert_eq!(generator.fingerprint_arguments, args);
        assert_eq!(
            generator.atom_invariants_generator,
            MorganAtomInvariantsGenerator::Connectivity {
                include_ring_membership: true
            }
        );
        assert_eq!(
            generator.bond_invariants_generator,
            MorganBondInvariantsGenerator::new(false, true)
        );
        assert!(generator.owns_atom_invariants_generator);
        assert!(generator.owns_bond_invariants_generator);
        assert_eq!(
            generator.atom_environment_generator.infoString(),
            "MorganEnvironmentGenerator"
        );
    }

    #[test]
    fn get_morgan_generator_preserves_custom_generators_and_ownership_flags() {
        let args = MorganArguments::default();
        let atom_generator = MorganAtomInvariantsGenerator::Feature;
        let bond_generator = MorganBondInvariantsGenerator::new(false, true);

        let generator = getMorganGenerator(
            &args,
            Some(atom_generator.clone()),
            Some(bond_generator.clone()),
            false,
            true,
        );

        assert_eq!(generator.fingerprint_arguments, args);
        assert_eq!(generator.atom_invariants_generator, atom_generator);
        assert_eq!(generator.bond_invariants_generator, bond_generator);
        assert!(!generator.owns_atom_invariants_generator);
        assert!(generator.owns_bond_invariants_generator);
    }

    #[test]
    fn get_morgan_generator_with_params_wires_all_morgan_arguments() {
        let generator = getMorganGeneratorWithParams(
            5,
            true,
            true,
            false,
            true,
            true,
            Some(MorganAtomInvariantsGenerator::Feature),
            None,
            4096,
            vec![1, 3, 7],
            false,
            false,
        )
        .unwrap();

        let args = &generator.fingerprint_arguments;
        assert_eq!(args.d_radius, 5);
        assert!(args.fingerprint_arguments.df_count_simulation);
        assert!(args.fingerprint_arguments.df_include_chirality);
        assert_eq!(args.fingerprint_arguments.d_count_bounds, vec![1, 3, 7]);
        assert_eq!(args.fingerprint_arguments.d_fp_size, 4096);
        assert_eq!(args.fingerprint_arguments.d_num_bits_per_feature, 1);
        assert!(args.df_only_nonzero_invariants);
        assert!(args.df_include_redundant_environments);
        assert!(!args.df_use_bond_types);
        assert_eq!(
            generator.atom_invariants_generator,
            MorganAtomInvariantsGenerator::Feature
        );
        assert_eq!(
            generator.bond_invariants_generator,
            MorganBondInvariantsGenerator::new(false, true)
        );
        assert!(!generator.owns_atom_invariants_generator);
        assert!(generator.owns_bond_invariants_generator);
    }

    #[test]
    fn get_morgan_generator_with_params_rejects_empty_count_bounds_when_counting() {
        let error = getMorganGeneratorWithParams(
            2,
            true,
            false,
            true,
            false,
            false,
            None,
            None,
            2048,
            Vec::new(),
            false,
            false,
        )
        .unwrap_err();

        assert_eq!(
            error,
            FingerprintError::InvalidArguments {
                reason: "countSimulation requires non-empty countBounds"
            }
        );
    }

    #[test]
    fn fingerprint_generator_helper_morgan_default_generator_matches_rdkit_env_golden() {
        let mol = rdkit_morgan_oracle_mol("CC");
        let generator = getMorganGeneratorWithParams(
            1,
            false,
            false,
            true,
            false,
            false,
            None,
            None,
            2048,
            vec![1, 2, 4, 8],
            false,
            false,
        )
        .unwrap();
        let mut args = FingerprintFuncArguments::default();

        let result = generator
            .getFingerprintHelper(&mol, &mut args, 2048)
            .unwrap();

        assert_eq!(result.size(), 2048);
        assert_eq!(
            result.nonzero_elements(),
            &BTreeMap::from([(1057, 2), (1275, 1)])
        );
    }

    #[test]
    fn fingerprint_generator_helper_morgan_custom_invariants_drive_environment_bits() {
        let mol = rdkit_morgan_oracle_mol("CC");
        let generator = getMorganGeneratorWithParams(
            1,
            false,
            false,
            true,
            false,
            true,
            None,
            None,
            2048,
            vec![1, 2, 4, 8],
            false,
            false,
        )
        .unwrap();
        let mut args = FingerprintFuncArguments {
            custom_atom_invariants: Some(vec![10, 20]),
            custom_bond_invariants: Some(vec![7]),
            ..Default::default()
        };

        let result = generator
            .getFingerprintHelper(&mol, &mut args, 2048)
            .unwrap();

        assert_eq!(
            result.nonzero_elements(),
            &BTreeMap::from([(10, 1), (20, 1), (805, 1), (1170, 1)])
        );
    }

    #[test]
    fn fingerprint_generator_helper_morgan_from_atoms_filters_output_centers() {
        let mol = rdkit_morgan_oracle_mol("CC");
        let generator = getMorganGeneratorWithParams(
            1,
            false,
            false,
            true,
            false,
            true,
            None,
            None,
            2048,
            vec![1, 2, 4, 8],
            false,
            false,
        )
        .unwrap();
        let mut args = FingerprintFuncArguments {
            from_atoms: Some(vec![1]),
            custom_atom_invariants: Some(vec![10, 20]),
            custom_bond_invariants: Some(vec![7]),
            ..Default::default()
        };

        let result = generator
            .getFingerprintHelper(&mol, &mut args, 2048)
            .unwrap();

        assert_eq!(
            result.nonzero_elements(),
            &BTreeMap::from([(20, 1), (1170, 1)])
        );
    }

    #[test]
    fn fingerprint_generator_helper_morgan_ignore_atoms_matches_rdkit_unused_argument() {
        let mol = rdkit_morgan_oracle_mol("CC");
        let generator = getMorganGeneratorWithParams(
            1,
            false,
            false,
            true,
            false,
            true,
            None,
            None,
            2048,
            vec![1, 2, 4, 8],
            false,
            false,
        )
        .unwrap();
        let mut args = FingerprintFuncArguments {
            ignore_atoms: Some(vec![0]),
            custom_atom_invariants: Some(vec![10, 20]),
            custom_bond_invariants: Some(vec![7]),
            ..Default::default()
        };

        let result = generator
            .getFingerprintHelper(&mol, &mut args, 2048)
            .unwrap();

        assert_eq!(
            result.nonzero_elements(),
            &BTreeMap::from([(10, 1), (20, 1), (805, 1), (1170, 1)])
        );
    }

    #[test]
    fn fingerprint_generator_helper_morgan_updates_additional_output_in_bit_order() {
        let mol = rdkit_morgan_oracle_mol("CC");
        let generator = getMorganGeneratorWithParams(
            1,
            false,
            false,
            true,
            false,
            true,
            None,
            None,
            2048,
            vec![1, 2, 4, 8],
            false,
            false,
        )
        .unwrap();
        let mut additional_output = allocated_morgan_additional_output(mol.num_atoms());
        if let Some(atom_counts) = additional_output.atom_counts.as_mut() {
            atom_counts.push(99);
        }
        let mut args = FingerprintFuncArguments {
            custom_atom_invariants: Some(vec![10, 20]),
            custom_bond_invariants: Some(vec![7]),
            additional_output: Some(additional_output),
            ..Default::default()
        };

        let result = generator
            .getFingerprintHelper(&mol, &mut args, 2048)
            .unwrap();

        assert_eq!(
            result.nonzero_elements(),
            &BTreeMap::from([(10, 1), (20, 1), (805, 1), (1170, 1)])
        );
        let output = args.additional_output.unwrap();
        let output = morgan_additional_output_from_rdkit_output(output);
        assert_eq!(output.atom_counts, vec![2, 2]);
        assert_eq!(output.atom_to_bits, vec![vec![10, 805], vec![20, 1170]]);
        assert_eq!(
            output.bit_info_map,
            bimap(&[
                (10, &[(0, 0)]),
                (20, &[(1, 0)]),
                (805, &[(0, 1)]),
                (1170, &[(1, 1)])
            ])
        );
        assert_eq!(
            output.atoms_per_bit,
            atoms_per_bit(&[
                (10, &[&[0]]),
                (20, &[&[1]]),
                (805, &[&[0, 1]]),
                (1170, &[&[1, 0]])
            ])
        );
    }

    #[test]
    fn fingerprint_count_simulation_additional_output_duplicates_rdkit_fields() {
        let old_output = AdditionalOutput {
            atom_to_bits: Some(vec![vec![5, 8], vec![8], Vec::new()]),
            bit_info_map: Some(BTreeMap::from([(8, vec![(0, 1), (1, 1)])])),
            bit_paths: Some(BTreeMap::from([(8, vec![vec![0, 1]])])),
            atom_counts: Some(vec![9, 7, 3]),
            atoms_per_bit: Some(BTreeMap::from([(8, vec![vec![0, 1], vec![1, 0]])])),
        };
        let mut new_output = AdditionalOutput {
            atom_to_bits: Some(Vec::new()),
            bit_info_map: Some(BTreeMap::new()),
            bit_paths: Some(BTreeMap::new()),
            atom_counts: Some(vec![99]),
            atoms_per_bit: Some(BTreeMap::new()),
        };

        duplicate_additional_output_bit(&old_output, &mut new_output, 8, 23).unwrap();

        assert_eq!(
            new_output.atom_to_bits,
            Some(vec![vec![23], vec![23], Vec::new()])
        );
        assert_eq!(
            new_output.bit_info_map,
            Some(BTreeMap::from([(23, vec![(0, 1), (1, 1)])]))
        );
        assert_eq!(
            new_output.bit_paths,
            Some(BTreeMap::from([(23, vec![vec![0, 1]])]))
        );
        assert_eq!(new_output.atom_counts, Some(vec![99]));
        assert_eq!(
            new_output.atoms_per_bit,
            Some(BTreeMap::from([(23, vec![vec![0, 1], vec![1, 0]])]))
        );
    }

    #[test]
    fn fingerprint_count_simulation_additional_output_rejects_mismatched_allocations() {
        let old_output = AdditionalOutput {
            bit_info_map: Some(BTreeMap::new()),
            ..Default::default()
        };
        let mut new_output = AdditionalOutput::default();

        let error = duplicate_additional_output_bit(&old_output, &mut new_output, 1, 2)
            .expect_err("RDKit PRECONDITION mismatch should become a structured error");

        assert_eq!(
            error,
            FingerprintError::InvalidArguments {
                reason: "bitInfoMap not allocated"
            }
        );
    }

    #[test]
    fn fingerprint_count_simulation_additional_output_setup_allocates_matching_temp_fields() {
        let mut original = AdditionalOutput {
            atom_to_bits: Some(vec![vec![7], vec![8, 9], vec![10]]),
            bit_info_map: Some(BTreeMap::from([(7, vec![(0, 0)])])),
            bit_paths: None,
            atom_counts: Some(vec![3, 4, 5]),
            atoms_per_bit: Some(BTreeMap::from([(7, vec![vec![0]])])),
        };
        let mut args = FingerprintFuncArguments {
            additional_output: Some(original.clone()),
            ..Default::default()
        };
        let mut count_simulation_output = AdditionalOutput::default();

        setup_temp_additional_output(&mut args, &mut count_simulation_output, 2);

        assert_eq!(
            count_simulation_output,
            AdditionalOutput {
                atom_to_bits: Some(Vec::new()),
                bit_info_map: Some(BTreeMap::new()),
                bit_paths: None,
                atom_counts: Some(Vec::new()),
                atoms_per_bit: Some(BTreeMap::new()),
            }
        );
        original.reset_for_atom_count(2);
        assert_eq!(args.additional_output, Some(original));
    }

    #[test]
    fn morgan_fingerprint_generator_outputs_sparse_count_sparse_bit_count_and_explicit_ethane() {
        let mol = rdkit_morgan_oracle_mol("CC");
        let generator = getMorganGeneratorWithParams(
            1,
            false,
            false,
            true,
            false,
            false,
            None,
            None,
            2048,
            vec![1, 2, 4, 8],
            false,
            false,
        )
        .unwrap();

        let sparse_count = generator
            .getSparseCountFingerprint(&mol, &mut FingerprintFuncArguments::default())
            .unwrap();
        let sparse_bit = generator
            .getSparseFingerprint(&mol, &mut FingerprintFuncArguments::default())
            .unwrap();
        let hashed_count = generator
            .getCountFingerprint(&mol, &mut FingerprintFuncArguments::default())
            .unwrap();
        let explicit_bit = generator
            .getFingerprint(&mol, &mut FingerprintFuncArguments::default())
            .unwrap();

        assert_eq!(sparse_count.size(), u64::from(u32::MAX));
        assert_eq!(
            sparse_count.nonzero_elements(),
            &BTreeMap::from([(2246728737, 2), (3545175291, 1)])
        );
        assert_eq!(sparse_bit.size(), u64::from(u32::MAX));
        assert_eq!(
            sparse_bit.on_bits(),
            &BTreeSet::from([2246728737, 3545175291])
        );
        assert_eq!(hashed_count.size(), 2048);
        assert_eq!(
            hashed_count.nonzero_elements(),
            &BTreeMap::from([(1057, 2), (1275, 1)])
        );
        assert_eq!(explicit_bit.n_bits(), 2048);
        assert_eq!(explicit_bit.on_bits(), vec![1057, 1275]);
    }

    #[test]
    fn fingerprint_generator_architecture_has_one_family_neutral_projector() {
        let family_source = include_str!("fingerprint.rs");
        let production_family_source = family_source
            .split("// Tests")
            .next()
            .expect("fingerprint production source");
        let generator_source = include_str!("fingerprint/generator.rs");

        assert_eq!(
            generator_source
                .matches("for environment in environments")
                .count(),
            1,
            "the common core must own the only environment projection loop"
        );
        assert_eq!(
            generator_source
                .matches("pub(crate) fn get_sparse_count_fingerprint(")
                .count(),
            1
        );
        assert_eq!(
            generator_source
                .matches("pub(crate) fn get_sparse_fingerprint(")
                .count(),
            1
        );
        assert_eq!(
            generator_source
                .matches("pub(crate) fn get_count_fingerprint(")
                .count(),
            1
        );
        assert_eq!(
            generator_source
                .matches("pub(crate) fn get_fingerprint(")
                .count(),
            1
        );
        assert!(
            family_source.contains(
                "generator::FingerprintGenerator::new(self).get_sparse_count_fingerprint"
            )
        );
        assert!(
            family_source
                .contains("generator::FingerprintGenerator::new(self).get_sparse_fingerprint")
        );
        assert!(
            family_source
                .contains("generator::FingerprintGenerator::new(self).get_count_fingerprint")
        );
        assert!(
            family_source.contains("generator::FingerprintGenerator::new(self).get_fingerprint")
        );
        assert!(!production_family_source.contains("fn build_fingerprint("));
        assert!(!production_family_source.contains("fn compute_initial_invariants("));
        assert!(!production_family_source.contains("fn fold_invariant("));
        assert!(!production_family_source.contains("fn atom_is_excluded("));
    }

    #[test]
    fn shared_generator_preserves_morgan_json_outputs_and_all_provenance_fields() {
        let molecule = rdkit_morgan_oracle_mol("CCO");
        let mut original_arguments =
            MorganArguments::new(2, true, false, false, vec![1, 2, 4, 8], 256, false, true)
                .unwrap();
        original_arguments
            .fingerprint_arguments
            .d_num_bits_per_feature = 2;
        let json = original_arguments.toJSON();
        let mut restored_arguments = MorganArguments::default();
        restored_arguments.fromJSON(&json).unwrap();
        assert_eq!(restored_arguments, original_arguments);

        let original = getMorganGenerator(&original_arguments, None, None, false, false);
        let restored = getMorganGenerator(&restored_arguments, None, None, false, false);

        fn all_outputs(atom_count: usize) -> AdditionalOutput {
            let mut output = allocated_morgan_additional_output(atom_count);
            output.allocate_bit_paths();
            output
        }

        let mut original_args = FingerprintFuncArguments {
            additional_output: Some(all_outputs(molecule.num_atoms())),
            ..Default::default()
        };
        let mut restored_args = original_args.clone();
        assert_eq!(
            original
                .getSparseCountFingerprint(&molecule, &mut original_args)
                .unwrap(),
            restored
                .getSparseCountFingerprint(&molecule, &mut restored_args)
                .unwrap()
        );
        assert_eq!(
            original_args.additional_output,
            restored_args.additional_output
        );

        let mut original_args = FingerprintFuncArguments {
            additional_output: Some(all_outputs(molecule.num_atoms())),
            ..Default::default()
        };
        let mut restored_args = original_args.clone();
        assert_eq!(
            original
                .getSparseFingerprint(&molecule, &mut original_args)
                .unwrap(),
            restored
                .getSparseFingerprint(&molecule, &mut restored_args)
                .unwrap()
        );
        assert_eq!(
            original_args.additional_output,
            restored_args.additional_output
        );

        let mut original_args = FingerprintFuncArguments {
            additional_output: Some(all_outputs(molecule.num_atoms())),
            ..Default::default()
        };
        let mut restored_args = original_args.clone();
        assert_eq!(
            original
                .getCountFingerprint(&molecule, &mut original_args)
                .unwrap(),
            restored
                .getCountFingerprint(&molecule, &mut restored_args)
                .unwrap()
        );
        assert_eq!(
            original_args.additional_output,
            restored_args.additional_output
        );

        let mut original_args = FingerprintFuncArguments {
            additional_output: Some(all_outputs(molecule.num_atoms())),
            ..Default::default()
        };
        let mut restored_args = original_args.clone();
        assert_eq!(
            original
                .getFingerprint(&molecule, &mut original_args)
                .unwrap(),
            restored
                .getFingerprint(&molecule, &mut restored_args)
                .unwrap()
        );
        assert_eq!(
            original_args.additional_output,
            restored_args.additional_output
        );
    }

    #[test]
    fn morgan_fingerprint_generator_outputs_radius_and_fp_size_match_rdkit_golden() {
        let mol = rdkit_morgan_oracle_mol("CCC");
        let generator = getMorganGeneratorWithParams(
            2,
            false,
            true,
            true,
            false,
            false,
            None,
            None,
            128,
            vec![1, 2, 4, 8],
            false,
            false,
        )
        .unwrap();

        let sparse_count = generator
            .getSparseCountFingerprint(&mol, &mut FingerprintFuncArguments::default())
            .unwrap();
        let sparse_bit = generator
            .getSparseFingerprint(&mol, &mut FingerprintFuncArguments::default())
            .unwrap();
        let hashed_count = generator
            .getCountFingerprint(&mol, &mut FingerprintFuncArguments::default())
            .unwrap();
        let explicit_bit = generator
            .getFingerprint(&mol, &mut FingerprintFuncArguments::default())
            .unwrap();

        assert_eq!(
            sparse_count.nonzero_elements(),
            &BTreeMap::from([
                (2068133184, 1),
                (2245384272, 1),
                (2246728737, 2),
                (3542456614, 2),
            ])
        );
        assert_eq!(
            sparse_bit.on_bits(),
            &BTreeSet::from([2068133184, 2245384272, 2246728737, 3542456614])
        );
        assert_eq!(hashed_count.size(), 128);
        assert_eq!(
            hashed_count.nonzero_elements(),
            &BTreeMap::from([(33, 2), (38, 2), (64, 1), (80, 1)])
        );
        assert_eq!(explicit_bit.n_bits(), 128);
        assert_eq!(explicit_bit.on_bits(), vec![33, 38, 64, 80]);
    }

    #[test]
    fn morgan_fingerprint_generator_outputs_count_simulation_and_additional_output() {
        let mol = rdkit_morgan_oracle_mol("CC");
        let generator = getMorganGeneratorWithParams(
            1,
            true,
            false,
            true,
            false,
            false,
            None,
            None,
            2048,
            vec![1, 2, 4, 8],
            false,
            false,
        )
        .unwrap();
        let mut args = FingerprintFuncArguments {
            additional_output: Some(allocated_morgan_additional_output(mol.num_atoms())),
            ..Default::default()
        };

        let explicit_bit = generator.getFingerprint(&mol, &mut args).unwrap();

        assert_eq!(explicit_bit.on_bits(), vec![132, 133, 1004]);
        let output = morgan_additional_output_from_rdkit_output(args.additional_output.unwrap());
        assert_eq!(output.atom_counts, vec![2, 1]);
        assert_eq!(
            output.atom_to_bits,
            vec![vec![132, 133, 1004], vec![132, 133]]
        );
        assert_eq!(
            output.bit_info_map,
            bimap(&[
                (132, &[(0, 0), (1, 0)]),
                (133, &[(0, 0), (1, 0)]),
                (1004, &[(0, 1)])
            ])
        );
        assert_eq!(
            output.atoms_per_bit,
            atoms_per_bit(&[
                (132, &[&[0], &[1]]),
                (133, &[&[0], &[1]]),
                (1004, &[&[0, 1]])
            ])
        );
    }

    #[test]
    fn morgan_fingerprint_generator_outputs_sparse_count_simulation_additional_output() {
        let mol = rdkit_morgan_oracle_mol("CC");
        let generator = getMorganGeneratorWithParams(
            1,
            true,
            false,
            true,
            false,
            false,
            None,
            None,
            2048,
            vec![1, 2, 4, 8],
            false,
            false,
        )
        .unwrap();
        let mut args = FingerprintFuncArguments {
            additional_output: Some(allocated_morgan_additional_output(mol.num_atoms())),
            ..Default::default()
        };

        let sparse_bit = generator.getSparseFingerprint(&mol, &mut args).unwrap();

        assert_eq!(
            sparse_bit.on_bits(),
            &BTreeSet::from([396980364, 396980365, 1295799288])
        );
        let output = morgan_additional_output_from_rdkit_output(args.additional_output.unwrap());
        assert_eq!(output.atom_counts, vec![2, 1]);
        assert_eq!(
            output.atom_to_bits,
            vec![
                vec![396980364, 396980365, 1295799288],
                vec![396980364, 396980365]
            ]
        );
        assert_eq!(
            output.bit_info_map,
            bimap(&[
                (396980364, &[(0, 0), (1, 0)]),
                (396980365, &[(0, 0), (1, 0)]),
                (1295799288, &[(0, 1)])
            ])
        );
        assert_eq!(
            output.atoms_per_bit,
            atoms_per_bit(&[
                (396980364, &[&[0], &[1]]),
                (396980365, &[&[0], &[1]]),
                (1295799288, &[&[0, 1]])
            ])
        );
    }

    #[test]
    fn morgan_fingerprint_generator_outputs_count_simulation_redundant_environments() {
        let mol = rdkit_morgan_oracle_mol("CC");
        let generator = getMorganGeneratorWithParams(
            1,
            true,
            false,
            true,
            false,
            true,
            None,
            None,
            2048,
            vec![1, 2, 4, 8],
            false,
            false,
        )
        .unwrap();
        let mut args = FingerprintFuncArguments {
            additional_output: Some(allocated_morgan_additional_output(mol.num_atoms())),
            ..Default::default()
        };

        let explicit_bit = generator.getFingerprint(&mol, &mut args).unwrap();

        assert_eq!(explicit_bit.on_bits(), vec![132, 133, 1004, 1005]);
        let output = morgan_additional_output_from_rdkit_output(args.additional_output.unwrap());
        assert_eq!(output.atom_counts, vec![2, 2]);
        assert_eq!(
            output.atom_to_bits,
            vec![vec![132, 133, 1004, 1005], vec![132, 133, 1004, 1005]]
        );
        assert_eq!(
            output.bit_info_map,
            bimap(&[
                (132, &[(0, 0), (1, 0)]),
                (133, &[(0, 0), (1, 0)]),
                (1004, &[(0, 1), (1, 1)]),
                (1005, &[(0, 1), (1, 1)])
            ])
        );
        assert_eq!(
            output.atoms_per_bit,
            atoms_per_bit(&[
                (132, &[&[0], &[1]]),
                (133, &[&[0], &[1]]),
                (1004, &[&[0, 1], &[1, 0]]),
                (1005, &[&[0, 1], &[1, 0]])
            ])
        );
    }

    #[test]
    fn morgan_fingerprint_include_chirality_uses_legacy_stereo_labels_by_default() {
        let missing_cip = stereo_done_double_bond_without_cip_computed(BondStereo::Trans);
        assert_eq!(missing_cip.prop("_StereochemDone"), Some("1"));
        assert_eq!(missing_cip.prop("_CIPComputed"), None);
        assert_eq!(missing_cip.bonds()[0].prop("_CIPCode"), None);

        let precomputed_cip = assign_cip_labels(&missing_cip, 0).unwrap();
        assert_eq!(precomputed_cip.prop("_CIPComputed"), Some("1"));
        assert_eq!(precomputed_cip.bonds()[0].prop("_CIPCode"), Some("E"));

        let generator = getMorganGeneratorWithParams(
            2,
            false,
            true,
            true,
            false,
            false,
            None,
            None,
            2048,
            vec![1, 2, 4, 8],
            false,
            false,
        )
        .unwrap();

        let sparse_count_missing = generator
            .getSparseCountFingerprint(&missing_cip, &mut FingerprintFuncArguments::default())
            .unwrap();
        let sparse_count_precomputed = generator
            .getSparseCountFingerprint(&precomputed_cip, &mut FingerprintFuncArguments::default())
            .unwrap();
        assert_ne!(sparse_count_missing, sparse_count_precomputed);

        let sparse_bit_missing = generator
            .getSparseFingerprint(&missing_cip, &mut FingerprintFuncArguments::default())
            .unwrap();
        let sparse_bit_precomputed = generator
            .getSparseFingerprint(&precomputed_cip, &mut FingerprintFuncArguments::default())
            .unwrap();
        assert_ne!(sparse_bit_missing, sparse_bit_precomputed);

        let hashed_count_missing = generator
            .getCountFingerprint(&missing_cip, &mut FingerprintFuncArguments::default())
            .unwrap();
        let hashed_count_precomputed = generator
            .getCountFingerprint(&precomputed_cip, &mut FingerprintFuncArguments::default())
            .unwrap();
        assert_ne!(hashed_count_missing, hashed_count_precomputed);

        let mut explicit_missing_args = FingerprintFuncArguments {
            additional_output: Some(allocated_morgan_additional_output(missing_cip.num_atoms())),
            ..Default::default()
        };
        let explicit_bit_missing = generator
            .getFingerprint(&missing_cip, &mut explicit_missing_args)
            .unwrap();
        let explicit_missing_output = morgan_additional_output_from_rdkit_output(
            explicit_missing_args.additional_output.unwrap(),
        );

        let mut explicit_precomputed_args = FingerprintFuncArguments {
            additional_output: Some(allocated_morgan_additional_output(
                precomputed_cip.num_atoms(),
            )),
            ..Default::default()
        };
        let explicit_bit_precomputed = generator
            .getFingerprint(&precomputed_cip, &mut explicit_precomputed_args)
            .unwrap();
        let explicit_precomputed_output = morgan_additional_output_from_rdkit_output(
            explicit_precomputed_args.additional_output.unwrap(),
        );

        assert_ne!(explicit_bit_missing, explicit_bit_precomputed);
        assert_ne!(explicit_missing_output, explicit_precomputed_output);
        assert!(!sparse_count_missing.nonzero_elements().is_empty());
        assert!(!sparse_bit_missing.on_bits().is_empty());
        assert!(!hashed_count_missing.nonzero_elements().is_empty());
        assert!(!explicit_bit_missing.on_bits().is_empty());
        assert!(!explicit_missing_output.bit_info_map.is_empty());
    }

    #[test]
    fn morgan_fingerprints_get_fingerprint_sparse_counts_and_atoms_setting_bits_match_rdkit() {
        let mol = rdkit_morgan_oracle_mol("CC");

        let output =
            morgan_get_fingerprint(&mol, 1, None, None, false, true, true, false, true, false)
                .unwrap();

        assert_eq!(output.fingerprint.size(), u64::from(u32::MAX));
        assert_eq!(
            output.fingerprint.nonzero_elements(),
            &BTreeMap::from([(2246728737, 2), (3545175291, 1)])
        );
        assert_eq!(
            output.atoms_setting_bits,
            Some(bimap(&[
                (2246728737, &[(0, 0), (1, 0)]),
                (3545175291, &[(0, 1)])
            ]))
        );
    }

    #[test]
    fn morgan_fingerprints_get_fingerprint_sparse_bits_are_returned_as_sparse_int_values() {
        let mol = rdkit_morgan_oracle_mol("CC");

        let output =
            morgan_get_fingerprint(&mol, 1, None, None, false, true, false, false, true, false)
                .unwrap();

        assert_eq!(output.fingerprint.size(), u64::from(u32::MAX));
        assert_eq!(
            output.fingerprint.nonzero_elements(),
            &BTreeMap::from([(2246728737, 1), (3545175291, 1)])
        );
        assert_eq!(
            output.atoms_setting_bits,
            Some(bimap(&[
                (2246728737, &[(0, 0), (1, 0)]),
                (3545175291, &[(0, 1)])
            ]))
        );
    }

    #[test]
    fn morgan_fingerprints_get_fingerprint_custom_invariants_and_from_atoms_match_rdkit() {
        let mol = rdkit_morgan_oracle_mol("CC");

        let output = morgan_get_fingerprint(
            &mol,
            1,
            Some(vec![10, 20]),
            Some(vec![1]),
            false,
            true,
            true,
            false,
            true,
            true,
        )
        .unwrap();

        assert_eq!(
            output.fingerprint.nonzero_elements(),
            &BTreeMap::from([(20, 1), (3205493690, 1)])
        );
        assert_eq!(
            output.atoms_setting_bits,
            Some(bimap(&[(20, &[(1, 0)]), (3205493690, &[(1, 1)])]))
        );
    }

    #[test]
    fn morgan_fingerprints_get_hashed_fingerprint_matches_rdkit_counts_and_atoms_setting_bits() {
        let mol = rdkit_morgan_oracle_mol("CC");

        let output = morgan_get_hashed_fingerprint(
            &mol, 1, 2048, None, None, false, true, false, true, false,
        )
        .unwrap();

        assert_eq!(output.fingerprint.size(), 2048);
        assert_eq!(
            output.fingerprint.nonzero_elements(),
            &BTreeMap::from([(1057, 2), (1275, 1)])
        );
        assert_eq!(
            output.atoms_setting_bits,
            Some(bimap(&[(1057, &[(0, 0), (1, 0)]), (1275, &[(0, 1)])]))
        );
    }

    #[test]
    fn morgan_fingerprints_get_hashed_fingerprint_custom_invariants_and_from_atoms_match_rdkit() {
        let mol = rdkit_morgan_oracle_mol("CC");

        let output = morgan_get_hashed_fingerprint(
            &mol,
            1,
            2048,
            Some(vec![10, 20]),
            Some(vec![1]),
            false,
            true,
            false,
            true,
            true,
        )
        .unwrap();

        assert_eq!(
            output.fingerprint.nonzero_elements(),
            &BTreeMap::from([(20, 1), (954, 1)])
        );
        assert_eq!(
            output.atoms_setting_bits,
            Some(bimap(&[(20, &[(1, 0)]), (954, &[(1, 1)])]))
        );
    }

    #[test]
    fn morgan_fingerprints_get_hashed_fingerprint_rejects_zero_n_bits() {
        let mol = rdkit_morgan_oracle_mol("CC");

        let error =
            morgan_get_hashed_fingerprint(&mol, 1, 0, None, None, false, true, false, false, false)
                .unwrap_err();

        assert_eq!(
            error,
            FingerprintError::InvalidArguments {
                reason: "nBits can not be zero"
            }
        );
    }

    #[test]
    fn morgan_fingerprints_get_fingerprint_as_bit_vect_matches_rdkit_bits_and_atoms_setting_bits() {
        let mol = rdkit_morgan_oracle_mol("CC");

        let output = morgan_get_fingerprint_as_bit_vect(
            &mol, 1, 2048, None, None, false, true, false, true, false,
        )
        .unwrap();

        assert_eq!(output.fingerprint.n_bits(), 2048);
        assert_eq!(output.fingerprint.on_bits(), vec![1057, 1275]);
        assert_eq!(
            output.atoms_setting_bits,
            Some(bimap(&[(1057, &[(0, 0), (1, 0)]), (1275, &[(0, 1)])]))
        );
    }

    #[test]
    fn morgan_fingerprints_get_fingerprint_as_bit_vect_custom_invariants_and_from_atoms_match_rdkit()
     {
        let mol = rdkit_morgan_oracle_mol("CC");

        let output = morgan_get_fingerprint_as_bit_vect(
            &mol,
            1,
            2048,
            Some(vec![10, 20]),
            Some(vec![1]),
            false,
            true,
            false,
            true,
            true,
        )
        .unwrap();

        assert_eq!(output.fingerprint.on_bits(), vec![20, 954]);
        assert_eq!(
            output.atoms_setting_bits,
            Some(bimap(&[(20, &[(1, 0)]), (954, &[(1, 1)])]))
        );
    }

    #[test]
    fn morgan_fingerprints_get_fingerprint_as_bit_vect_rejects_zero_n_bits() {
        let mol = rdkit_morgan_oracle_mol("CC");

        let error = morgan_get_fingerprint_as_bit_vect(
            &mol, 1, 0, None, None, false, true, false, false, false,
        )
        .unwrap_err();

        assert_eq!(
            error,
            FingerprintError::InvalidArguments {
                reason: "nBits can not be zero"
            }
        );
    }

    #[test]
    fn morgan_env_generator_radius_zero_matches_rdkit_oracle_for_ethane() {
        let mol = rdkit_morgan_oracle_mol("CC");
        let args = morgan_arguments_for_env(0, false, false, true);

        let output = rdkit_morgan_env_output(&mol, &args, None, None).unwrap();

        assert_eq!(output.atom_counts, vec![1, 1]);
        assert_eq!(
            output.atom_to_bits,
            vec![vec![2246728737], vec![2246728737]]
        );
        assert_eq!(
            output.bit_info_map,
            bimap(&[(2246728737, &[(0, 0), (1, 0)])])
        );
        assert_eq!(
            output.atoms_per_bit,
            atoms_per_bit(&[(2246728737, &[&[0], &[1]])])
        );
    }

    #[test]
    fn morgan_env_generator_duplicate_environments_follow_redundant_flag() {
        let mol = rdkit_morgan_oracle_mol("CC");

        let nonredundant = rdkit_morgan_env_output(
            &mol,
            &morgan_arguments_for_env(1, false, false, true),
            None,
            None,
        )
        .unwrap();
        assert_eq!(nonredundant.atom_counts, vec![2, 1]);
        assert_eq!(
            nonredundant.bit_info_map,
            bimap(&[(2246728737, &[(0, 0), (1, 0)]), (3545175291, &[(0, 1)]),])
        );
        assert_eq!(
            nonredundant.atoms_per_bit,
            atoms_per_bit(&[(2246728737, &[&[0], &[1]]), (3545175291, &[&[0, 1]]),])
        );

        let redundant = rdkit_morgan_env_output(
            &mol,
            &morgan_arguments_for_env(1, true, false, true),
            None,
            None,
        )
        .unwrap();
        assert_eq!(redundant.atom_counts, vec![2, 2]);
        assert_eq!(
            redundant.bit_info_map,
            bimap(&[
                (2246728737, &[(0, 0), (1, 0)]),
                (3545175291, &[(0, 1), (1, 1)]),
            ])
        );
        assert_eq!(
            redundant.atoms_per_bit,
            atoms_per_bit(&[
                (2246728737, &[&[0], &[1]]),
                (3545175291, &[&[0, 1], &[1, 0]]),
            ])
        );
    }

    #[test]
    fn morgan_env_generator_radius_two_dead_atoms_and_layer_semantics_match_rdkit() {
        let mol = rdkit_morgan_oracle_mol("CCC");
        let args = morgan_arguments_for_env(2, false, false, true);

        let output = rdkit_morgan_env_output(&mol, &args, None, None).unwrap();

        assert_eq!(output.atom_counts, vec![2, 2, 2]);
        assert_eq!(
            output.bit_info_map,
            bimap(&[
                (2068133184, &[(1, 1)]),
                (2245384272, &[(1, 0)]),
                (2246728737, &[(0, 0), (2, 0)]),
                (3542456614, &[(0, 1), (2, 1)]),
            ])
        );
        assert_eq!(
            output.atoms_per_bit,
            atoms_per_bit(&[
                (2068133184, &[&[1, 0, 2]]),
                (2245384272, &[&[1]]),
                (2246728737, &[&[0], &[2]]),
                (3542456614, &[&[0, 1], &[2, 1]]),
            ])
        );
    }

    #[test]
    fn morgan_env_generator_handles_isolated_atoms_like_rdkit() {
        let mol = rdkit_morgan_oracle_mol("[He]");
        let args = morgan_arguments_for_env(2, false, false, true);

        let output = rdkit_morgan_env_output(&mol, &args, None, None).unwrap();

        assert_eq!(output.atom_counts, vec![1]);
        assert_eq!(output.atom_to_bits, vec![vec![2312954353]]);
        assert_eq!(output.bit_info_map, bimap(&[(2312954353, &[(0, 0)])]));
        assert_eq!(
            output.atoms_per_bit,
            atoms_per_bit(&[(2312954353, &[&[0]])])
        );
    }

    #[test]
    fn morgan_env_generator_from_atoms_filters_output_centers_only() {
        let mol = rdkit_morgan_oracle_mol("CCO");
        let args = morgan_arguments_for_env(1, false, false, true);

        let output = rdkit_morgan_env_output(&mol, &args, Some(&[0]), None).unwrap();

        assert_eq!(output.atom_counts, vec![2, 0, 0]);
        assert_eq!(
            output.atom_to_bits,
            vec![vec![2246728737, 3542456614], Vec::new(), Vec::new()]
        );
        assert_eq!(
            output.bit_info_map,
            bimap(&[(2246728737, &[(0, 0)]), (3542456614, &[(0, 1)])])
        );
        assert_eq!(
            output.atoms_per_bit,
            atoms_per_bit(&[(2246728737, &[&[0]]), (3542456614, &[&[0, 1]])])
        );
    }

    #[test]
    fn morgan_env_generator_ignore_atoms_argument_is_unused_by_rdkit_source() {
        let mol = rdkit_morgan_oracle_mol("CCO");
        let args = morgan_arguments_for_env(1, false, false, true);

        let output = rdkit_morgan_env_output(&mol, &args, None, Some(&[0])).unwrap();

        assert_eq!(output.atom_counts, vec![2, 2, 2]);
        assert_eq!(
            output.bit_info_map,
            bimap(&[
                (864662311, &[(2, 0)]),
                (1535166686, &[(2, 1)]),
                (2245384272, &[(1, 0)]),
                (2246728737, &[(0, 0)]),
                (3542456614, &[(0, 1)]),
                (4018048386, &[(1, 1)]),
            ])
        );
        assert_eq!(
            output.atoms_per_bit,
            atoms_per_bit(&[
                (864662311, &[&[2]]),
                (1535166686, &[&[2, 1]]),
                (2245384272, &[&[1]]),
                (2246728737, &[&[0]]),
                (3542456614, &[&[0, 1]]),
                (4018048386, &[&[1, 0, 2]]),
            ])
        );
    }

    #[test]
    fn morgan_env_generator_bond_type_branch_changes_double_bond_environment() {
        let mol = rdkit_morgan_oracle_mol("C=C");

        let with_bond_types = rdkit_morgan_env_output(
            &mol,
            &morgan_arguments_for_env(1, true, false, true),
            None,
            None,
        )
        .unwrap();
        let without_bond_types = rdkit_morgan_env_output(
            &mol,
            &morgan_arguments_for_env(1, true, false, false),
            None,
            None,
        )
        .unwrap();

        assert_eq!(
            with_bond_types.bit_info_map,
            bimap(&[
                (2246997334, &[(0, 0), (1, 0)]),
                (3695448525, &[(0, 1), (1, 1)]),
            ])
        );
        assert_eq!(
            without_bond_types.bit_info_map,
            bimap(&[
                (2246997334, &[(0, 0), (1, 0)]),
                (3695449228, &[(0, 1), (1, 1)]),
            ])
        );
    }

    #[test]
    fn morgan_env_generator_chirality_preserves_legacy_double_bond_inputs_by_default() {
        let missing_cip = stereo_done_double_bond_without_cip_computed(BondStereo::Trans);
        let args = morgan_arguments_for_env(1, false, true, true);

        let legacy_output = rdkit_morgan_env_output(&missing_cip, &args, None, None).unwrap();

        let with_cip = assign_cip_labels(&missing_cip, 0).unwrap();
        let output = rdkit_morgan_env_output(&with_cip, &args, None, None).unwrap();
        assert_eq!(legacy_output, output);
        assert_eq!(with_cip.prop("_CIPComputed"), Some("1"));
        assert!(
            with_cip.bonds()[0].prop("_CIPCode").is_some(),
            "ported CIPLabeler remains available for non-legacy stereo mode"
        );
    }

    #[test]
    fn fingerprint_arguments_defaults_match_rdkit_shape() {
        let args = FingerprintArguments::default();
        assert!(!args.df_count_simulation);
        assert!(!args.df_include_chirality);
        assert_eq!(args.d_count_bounds, Vec::<u32>::new());
        assert_eq!(args.d_fp_size, 2048);
        assert_eq!(args.d_num_bits_per_feature, 1);
        assert_eq!(
            args.commonArgumentsString(),
            "Common arguments : countSimulation=0 fpSize=2048 bitsPerFeature=1 includeChirality=0"
        );
    }

    #[test]
    fn fingerprint_arguments_rejects_invalid_count_bounds_and_feature_width() {
        assert!(matches!(
            FingerprintArguments::new(true, Vec::new(), 2048, 1, false),
            Err(FingerprintError::InvalidArguments { .. })
        ));
        assert!(matches!(
            FingerprintArguments::new(false, Vec::new(), 2048, 0, false),
            Err(FingerprintError::InvalidArguments { .. })
        ));
    }

    #[test]
    fn fingerprint_arguments_to_json_matches_rdkit_field_shape() {
        let args = FingerprintArguments::new(true, vec![1, 2, 4, 8], 4096, 3, true).unwrap();
        assert_eq!(
            args.toJSON(),
            r#"{"countSimulation":"true","fpSize":"4096","numBitsPerFeature":"3","includeChirality":"true","countBounds":["1","2","4","8"]}"#
        );
    }

    #[test]
    fn fingerprint_arguments_from_json_updates_present_fields_and_clears_bounds() {
        let mut args = FingerprintArguments::new(true, vec![1, 2, 4, 8], 4096, 3, true).unwrap();
        args.fromJSON(r#"{"countSimulation":false,"fpSize":1024,"numBitsPerFeature":2,"includeChirality":false,"countBounds":[2,5,9]}"#)
            .unwrap();
        assert!(!args.df_count_simulation);
        assert!(!args.df_include_chirality);
        assert_eq!(args.d_fp_size, 1024);
        assert_eq!(args.d_num_bits_per_feature, 2);
        assert_eq!(args.d_count_bounds, vec![2, 5, 9]);
    }

    #[test]
    fn fingerprint_arguments_from_json_preserves_unmentioned_fields_and_rejects_invalid_json() {
        let mut args = FingerprintArguments::new(false, vec![1, 2], 2048, 1, true).unwrap();
        args.fromJSON(r#"{"fpSize":512}"#).unwrap();
        assert!(!args.df_count_simulation);
        assert!(args.df_include_chirality);
        assert_eq!(args.d_fp_size, 512);
        assert_eq!(args.d_count_bounds, Vec::<u32>::new());
        assert!(args.fromJSON("{").is_err());
    }

    #[test]
    fn morgan_arguments_defaults_match_rdkit_shape() {
        let args = MorganArguments::default();
        assert_eq!(args.d_radius, 3);
        assert!(!args.df_only_nonzero_invariants);
        assert!(!args.df_include_redundant_environments);
        assert!(args.df_use_bond_types);
        assert!(!args.fingerprint_arguments.df_count_simulation);
        assert!(!args.fingerprint_arguments.df_include_chirality);
        assert_eq!(args.fingerprint_arguments.d_count_bounds, vec![1, 2, 4, 8]);
        assert_eq!(args.fingerprint_arguments.d_fp_size, 2048);
        assert_eq!(args.fingerprint_arguments.d_num_bits_per_feature, 1);
        assert_eq!(
            args.infoString(),
            "MorganArguments onlyNonzeroInvariants=0 radius=3"
        );
    }

    #[test]
    fn morgan_arguments_to_json_matches_rdkit_field_shape() {
        let args =
            MorganArguments::new(2, true, true, true, vec![1, 2, 4, 8], 1024, true, false).unwrap();
        assert_eq!(
            args.toJSON(),
            r#"{"type":"MorganArguments","onlyNonzeroInvariants":true,"radius":2,"countSimulation":"true","fpSize":"1024","numBitsPerFeature":"1","includeChirality":"true","countBounds":["1","2","4","8"]}"#
        );
    }

    #[test]
    fn morgan_arguments_constructor_sets_redundant_environments_and_bond_types() {
        let args =
            MorganArguments::new(4, false, false, false, vec![1, 2, 4, 8], 2048, true, false)
                .unwrap();
        assert_eq!(args.d_radius, 4);
        assert!(args.df_include_redundant_environments);
        assert!(!args.df_use_bond_types);
    }

    #[test]
    fn morgan_arguments_from_json_updates_present_fields_and_round_trips() {
        let mut args = MorganArguments::default();
        args.fromJSON(
            r#"{"radius":5,"onlyNonzeroInvariants":true,"countSimulation":true,"fpSize":512,"numBitsPerFeature":2,"includeChirality":true,"countBounds":[3,9]}"#,
        )
        .unwrap();
        assert_eq!(args.d_radius, 5);
        assert!(args.df_only_nonzero_invariants);
        assert!(args.fingerprint_arguments.df_count_simulation);
        assert_eq!(args.fingerprint_arguments.d_fp_size, 512);
        assert_eq!(args.fingerprint_arguments.d_num_bits_per_feature, 2);
        assert!(args.fingerprint_arguments.df_include_chirality);
        assert_eq!(args.fingerprint_arguments.d_count_bounds, vec![3, 9]);
        let json = args.toJSON();
        let mut round_trip = MorganArguments::default();
        round_trip.fromJSON(&json).unwrap();
        assert_eq!(round_trip, args);
    }

    #[test]
    fn morgan_atom_inv_generator_defaults_and_round_trip_match_rdkit_shape() {
        let mut generator = MorganAtomInvGenerator::new(true);
        assert_eq!(
            generator.infoString(),
            "MorganInvariantGenerator includeRingMembership=1"
        );
        assert_eq!(
            generator.toJSON(),
            r#"{"type":"MorganAtomInvGenerator","includeRingMembership":true}"#
        );
        let clone = generator.clone();
        assert_eq!(clone, generator);

        generator
            .fromJSON(r#"{"includeRingMembership":false}"#)
            .unwrap();
        assert_eq!(
            generator.toJSON(),
            r#"{"type":"MorganAtomInvGenerator","includeRingMembership":false}"#
        );
        let mut round_trip = MorganAtomInvGenerator::new(true);
        round_trip.fromJSON(&generator.toJSON()).unwrap();
        assert_eq!(round_trip, generator);
    }

    #[test]
    fn morgan_atom_inv_generator_matches_connectivity_branch_shape() {
        let mol = Molecule::from_smiles("C1CC1C").unwrap();
        let ring_generator = MorganAtomInvGenerator::new(true);
        let plain_generator = MorganAtomInvGenerator::new(false);
        let ring_invariants = ring_generator.getAtomInvariants(&mol);
        let plain_invariants = plain_generator.getAtomInvariants(&mol);
        assert_eq!(ring_invariants.len(), mol.num_atoms());
        assert_eq!(plain_invariants.len(), mol.num_atoms());
        assert!(
            ring_invariants
                .iter()
                .zip(&plain_invariants)
                .any(|(left, right)| left != right)
        );
    }

    #[test]
    fn morgan_atom_inv_generator_matches_rdkit_golden_connectivity_invariants() {
        let cases = [
            ("C", true, vec![2246733040]),
            ("C", false, vec![2246733040]),
            ("[NH4+]", true, vec![847680145]),
            ("[2H]", true, vec![4269929704]),
            ("[2H]O[2H]", true, vec![4277593707, 864666390, 4277593707]),
            ("CCO", true, vec![2246728737, 2245384272, 864662311]),
            (
                "c1ccccc1",
                true,
                vec![
                    3218693969, 3218693969, 3218693969, 3218693969, 3218693969, 3218693969,
                ],
            ),
            (
                "c1ccccc1",
                false,
                vec![
                    2246703798, 2246703798, 2246703798, 2246703798, 2246703798, 2246703798,
                ],
            ),
        ];

        for (smiles, include_ring_membership, expected) in cases {
            let mol = Molecule::from_smiles(smiles).unwrap();
            let invariants =
                MorganAtomInvGenerator::new(include_ring_membership).getAtomInvariants(&mol);
            assert_eq!(
                invariants, expected,
                "{smiles} include_ring={include_ring_membership}"
            );
        }
    }

    #[test]
    fn smarts_consumer_fingerprint_patterns() {
        assert_eq!(
            default_feature_smarts(),
            &[
                "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]",
                "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]",
                "[a]",
                "[F,Cl,Br,I]",
                "[#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",
                "[$([C,S](=[O,S,P])-[O;H1,-1])]",
            ]
        );

        let matchers = default_feature_matchers().expect("default feature SMARTS must parse");
        assert_eq!(matchers.len(), 6);
        assert_eq!(matchers.iter().map(SsMatcher::getMatcher).count(), 6);
    }

    #[test]
    fn default_feature_smarts_match_expected_atoms_on_rdkit_source_fixtures() {
        let fixtures = [
            ("CCO", 0, vec![2]),
            ("CC(=O)C", 1, vec![2]),
            ("c1ccccc1", 2, vec![0, 1, 2, 3, 4, 5]),
            ("CCl", 3, vec![1]),
            ("[NH4+]", 4, vec![0]),
            ("CC(=O)O", 5, vec![1]),
        ];

        let matchers = default_feature_matchers().expect("default feature SMARTS must parse");
        for (smiles, feature_idx, expected_atoms) in fixtures {
            let mol = Molecule::from_smiles_with_sanitize(smiles, false).unwrap();
            let query = matchers[feature_idx].getMatcher();
            let matches = crate::get_substruct_matches(&mol, query);
            let mut atom_indices: Vec<usize> = matches
                .iter()
                .flat_map(|matched| matched.atom_mapping.iter().copied())
                .collect();
            atom_indices.sort_unstable();
            atom_indices.dedup();
            let mut expected = expected_atoms.clone();
            expected.sort_unstable();
            expected.dedup();
            assert_eq!(
                atom_indices, expected,
                "feature {feature_idx} mismatch for {smiles}"
            );
        }
    }

    #[test]
    fn morgan_feature_invariants_generator_default_patterns_match_rdkit_golden_vectors() {
        let cases = [
            ("CCO", vec![0, 0, 0x03]),
            ("CC(=O)C", vec![0, 0, 0x02, 0]),
            ("c1ccccc1", vec![0x04, 0x04, 0x04, 0x04, 0x04, 0x04]),
            ("CCl", vec![0, 0x08]),
            ("[NH4+]", vec![0x11]),
            ("CC(=O)O", vec![0, 0x20, 0x02, 0x01]),
        ];

        let generator = MorganFeatureAtomInvGenerator::new();
        for (smiles, expected) in cases {
            let mol = Molecule::from_smiles_with_sanitize(smiles, false).unwrap();
            assert_eq!(
                generator
                    .getAtomInvariants(&mol)
                    .expect("feature invariants"),
                expected,
                "default feature invariant mismatch for {smiles}"
            );
        }
    }

    #[test]
    fn morgan_feature_invariants_generator_metadata_json_and_clone_match_rdkit_shape() {
        let generator = MorganFeatureAtomInvGenerator::new();
        assert_eq!(generator.infoString(), "MorganFeatureInvariantGenerator");
        assert_eq!(
            generator.toJSON(),
            r#"{"type":"MorganFeatureAtomInvGenerator"}"#
        );
        assert_eq!(generator.clone(), generator);

        let mut round_trip = MorganFeatureAtomInvGenerator::new();
        round_trip.fromJSON(&generator.toJSON()).unwrap();
        assert_eq!(round_trip, generator);
        assert_eq!(round_trip.clone(), generator);
    }

    #[test]
    fn morgan_feature_invariants_generator_supplied_patterns_fail_closed_until_source_ported() {
        let err = MorganFeatureAtomInvGenerator::from_smarts_patterns(&["[#6]"]).unwrap_err();
        assert!(matches!(
            err,
            FingerprintError::UnsupportedOption {
                option: "MorganFeatureAtomInvGenerator.patternSMARTS",
                ..
            }
        ));

        let mut generator = MorganFeatureAtomInvGenerator::new();
        let err = generator
            .fromJSON(r#"{"patternSMARTS":["[#6]"]}"#)
            .unwrap_err();
        assert!(matches!(
            err,
            FingerprintError::UnsupportedOption {
                option: "MorganFeatureAtomInvGenerator.patternSMARTS",
                ..
            }
        ));
    }

    #[test]
    fn morgan_bond_invariants_generator_metadata_json_and_clone_match_rdkit_shape() {
        let mut generator = MorganBondInvGenerator::new(true, false);
        assert_eq!(
            generator.infoString(),
            "MorganInvariantGenerator useBondTypes=1 useChirality=0"
        );
        assert_eq!(
            generator.toJSON(),
            r#"{"type":"MorganBondInvGenerator","useBondTypes":true,"useChirality":false}"#
        );
        assert_eq!(generator.clone(), generator);

        generator
            .fromJSON(r#"{"useBondTypes":false,"useChirality":true}"#)
            .unwrap();
        assert_eq!(generator, MorganBondInvGenerator::new(false, true));
        let mut round_trip = MorganBondInvGenerator::new(true, false);
        round_trip.fromJSON(&generator.toJSON()).unwrap();
        assert_eq!(round_trip, generator);
    }

    #[test]
    fn morgan_bond_invariants_generator_matches_rdkit_bond_type_codes() {
        let cases = [
            (BondOrder::Unspecified, 0),
            (BondOrder::Single, 1),
            (BondOrder::Double, 2),
            (BondOrder::Triple, 3),
            (BondOrder::Quadruple, 4),
            (BondOrder::Quintuple, 5),
            (BondOrder::Hextuple, 6),
            (BondOrder::OneAndHalf, 7),
            (BondOrder::TwoAndHalf, 8),
            (BondOrder::ThreeAndHalf, 9),
            (BondOrder::FourAndHalf, 10),
            (BondOrder::FiveAndHalf, 11),
            (BondOrder::Aromatic, 12),
            (BondOrder::Ionic, 13),
            (BondOrder::Hydrogen, 14),
            (BondOrder::ThreeCenter, 15),
            (BondOrder::DativeOne, 16),
            (BondOrder::Dative, 17),
            (BondOrder::DativeLeft, 18),
            (BondOrder::DativeRight, 19),
            (BondOrder::Other, 20),
            (BondOrder::Zero, 21),
            (BondOrder::Null, 0),
        ];
        let generator = MorganBondInvGenerator::new(true, false);

        for (order, expected) in cases {
            let mol = two_atom_molecule_with_bond(order);
            assert_eq!(
                generator.getBondInvariants(&mol).expect("bond invariants"),
                vec![expected],
                "bond order {order:?}"
            );
        }
    }

    #[test]
    fn morgan_bond_invariants_generator_no_bond_type_mode_uses_rdkit_constant_one() {
        let generator = MorganBondInvGenerator::new(false, false);

        for order in [
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Aromatic,
            BondOrder::Dative,
            BondOrder::Zero,
        ] {
            let mol = two_atom_molecule_with_bond(order);
            assert_eq!(
                generator.getBondInvariants(&mol).expect("bond invariants"),
                vec![1],
                "{order:?}"
            );
        }
    }

    #[test]
    fn morgan_bond_invariants_generator_double_bond_chirality_matches_rdkit_formula() {
        let cases = [
            (BondStereo::None, 2),
            (BondStereo::Any, 121),
            (BondStereo::Cis, 124),
            (BondStereo::Trans, 125),
        ];
        let generator = MorganBondInvGenerator::new(true, true);

        for (stereo, expected) in cases {
            let mol = if matches!(stereo, BondStereo::Cis | BondStereo::Trans) {
                stereo_done_double_bond_without_cip_computed(stereo)
            } else {
                double_bond_stereo_molecule(stereo)
            };
            let invariants = generator.getBondInvariants(&mol).expect("bond invariants");
            assert_eq!(invariants[0], expected, "double bond stereo {stereo:?}");
            assert_eq!(&invariants[1..], &[1, 1], "side bonds for {stereo:?}");
        }
    }

    #[test]
    fn morgan_bond_invariants_generator_chirality_ignored_for_non_double_bonds() {
        let generator = MorganBondInvGenerator::new(true, true);
        let mol = two_atom_molecule_with_bond_spec(
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_stereo(BondStereo::E),
        );

        assert_eq!(
            generator.getBondInvariants(&mol).expect("bond invariants"),
            vec![1]
        );
    }

    #[test]
    fn morgan_feature_atom_invariants_use_source_feature_matchers() {
        let mol = Molecule::from_smiles_with_sanitize("OCCN", false).unwrap();
        let params = MorganFingerprintParams {
            radius: 0,
            n_bits: 2048,
            atom_invariants_generator: MorganAtomInvariantsGenerator::Feature,
            ..Default::default()
        };
        let output = morgan_fingerprint_with_output(&mol, &params).unwrap();
        assert!(!output.fingerprint.on_bits().is_empty());
        let invariants = compute_feature_invariants(&mol).expect("feature invariants");
        assert_eq!(invariants.len(), mol.num_atoms());
        assert!(invariants.iter().any(|inv| *inv != 0));
    }

    #[test]
    fn topological_torsion_atom_invariants_custom_correction_wraps_u32() {
        assert_eq!(topological_torsion_correct_atom_invariant(0), u32::MAX - 1);
        assert_eq!(topological_torsion_correct_atom_invariant(1), u32::MAX);
        assert_eq!(topological_torsion_correct_atom_invariant(2), 0);
        assert_eq!(
            topological_torsion_correct_atom_invariant(u32::MAX),
            u32::MAX - 2
        );
    }
}
