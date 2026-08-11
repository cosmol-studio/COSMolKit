use std::collections::{BTreeMap, BTreeSet, HashSet, VecDeque};
use std::sync::OnceLock;

use crate::chemistry::ciplabeler::assign_cip_labels;
use crate::chemistry::valence::rdkit_atomic_mass;
use crate::search::smarts_parse::build_query_molecule;
use crate::{AdjacencyList, AtomId, BondOrder, ChiralTag, Molecule};
use serde_json::Value;

// RDKit marker convention defined in dev/source_reproduction_protocol.md.
// Copied source lines appear as:  // RDKit<beh><perf>: ...

// RDKit source file: FingerprintUtil.cpp

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
        let matcher = build_query_molecule(pattern).map_err(|reason| {
            FingerprintError::InvalidSmartsPattern {
                pattern: pattern.to_string(),
                reason,
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

    /// RDKit✔️✔️: reinitAdditionalOutput(AdditionalOutput &ao, size_t numAtoms)
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
        if let Some(atoms_per_bit) = self.atoms_per_bit.as_mut() {
            atoms_per_bit.clear();
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
            .map(u32::to_string)
            .collect::<Vec<_>>()
            .join(",");
        format!(
            "{{\"countSimulation\":{},\"fpSize\":{},\"numBitsPerFeature\":{},\"includeChirality\":{},\"countBounds\":[{}]}}",
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

    /// RDKit✔️✔️: df_countSimulation = pt.get<bool>("countSimulation", df_countSimulation);
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

    /// RDKit✔️✔️: OutputType MorganAtomEnv<OutputType>::getBitId(...) const
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

    /// RDKit✔️❌: void MorganAtomEnv<OutputType>::updateAdditionalOutput(AdditionalOutput *additionalOutput, size_t bitId) const
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

    /// RDKit✔️✔️: std::string MorganEnvGenerator<OutputType>::infoString() const
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

    /// RDKit✔️✔️: void MorganEnvGenerator<OutputType>::toJSON(boost::property_tree::ptree &pt) const
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

    /// RDKit✔️✔️: void MorganEnvGenerator<OutputType>::fromJSON(const boost::property_tree::ptree &pt)
    #[allow(non_snake_case)]
    pub fn fromJSON(&mut self, json: &str) -> Result<(), FingerprintError> {
        // RDKit source: MorganGenerator.cpp lines 473-477
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: void MorganEnvGenerator<OutputType>::fromJSON(
        // RDKit✔️✔️:     const boost::property_tree::ptree &pt) {
        // RDKit✔️✔️:   AtomEnvironmentGenerator<OutputType>::fromJSON(pt);
        // RDKit✔️✔️: }
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

    /// RDKit✔️✔️: OutputType MorganEnvGenerator<OutputType>::getResultSize() const
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

    /// RDKit✔️❌: std::vector<AtomEnvironment<OutputType> *> MorganEnvGenerator<OutputType>::getEnvironments(...)
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

impl MorganFingerprintGenerator {
    #[allow(non_snake_case)]
    pub fn getSparseCountFingerprint(
        &self,
        molecule: &Molecule,
        args: &mut FingerprintFuncArguments,
    ) -> Result<SparseCountFingerprint, FingerprintError> {
        // RDKit source: FingerprintGenerator.cpp lines 505-508
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: std::unique_ptr<SparseIntVect<OutputType>>
        // RDKit✔️✔️: FingerprintGenerator<OutputType>::getSparseCountFingerprint(
        // RDKit✔️✔️:     const ROMol &mol, FingerprintFuncArguments &args) const {
        // RDKit✔️✔️:   return getFingerprintHelper(mol, args);
        // RDKit✔️✔️: }
        self.getFingerprintHelper(molecule, args, 0)
    }

    #[allow(non_snake_case)]
    pub fn getSparseFingerprint(
        &self,
        molecule: &Molecule,
        args: &mut FingerprintFuncArguments,
    ) -> Result<SparseBitFingerprint, FingerprintError> {
        // RDKit source: FingerprintGenerator.cpp lines 510-575
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: std::unique_ptr<SparseBitVect>
        // RDKit✔️✔️: FingerprintGenerator<OutputType>::getSparseFingerprint(
        // RDKit✔️✔️:     const ROMol &mol, FingerprintFuncArguments &args) const {
        // RDKit✔️✔️:   std::uint32_t resultSize =
        // RDKit✔️✔️:       std::min((std::uint64_t)std::numeric_limits<std::uint32_t>::max(),
        // RDKit✔️✔️:                (std::uint64_t)dp_atomEnvironmentGenerator->getResultSize());
        let result_size = self
            .atom_environment_generator
            .getResultSize()
            .min(u64::from(u32::MAX)) as u32;

        // RDKit✔️✔️:   std::uint32_t effectiveSize = resultSize;
        // RDKit✔️✔️:   if (dp_fingerprintArguments->df_countSimulation) {
        // RDKit✔️✔️:     effectiveSize /= dp_fingerprintArguments->d_countBounds.size();
        // RDKit✔️✔️:   }
        let count_simulation = self
            .fingerprint_arguments
            .fingerprint_arguments
            .df_count_simulation;
        let count_bounds = &self
            .fingerprint_arguments
            .fingerprint_arguments
            .d_count_bounds;
        let effective_size = if count_simulation {
            result_size / count_bounds.len() as u32
        } else {
            result_size
        };

        // RDKit✔️✔️:   AdditionalOutput countSimulationOutput;
        // RDKit✔️✔️:   AdditionalOutput *origAO = nullptr;
        // RDKit✔️✔️:   if (dp_fingerprintArguments->df_countSimulation && args.additionalOutput) {
        // RDKit✔️✔️:     setupTempAdditionalOutput(args, countSimulationOutput, mol.getNumAtoms());
        // RDKit✔️✔️:     origAO = args.additionalOutput;
        // RDKit✔️✔️:     args.additionalOutput = &countSimulationOutput;
        // RDKit✔️✔️:   }
        let mut original_additional_output = None;
        if count_simulation && args.additional_output.is_some() {
            let mut count_simulation_output = AdditionalOutput::default();
            setup_temp_additional_output(args, &mut count_simulation_output, molecule.num_atoms());
            original_additional_output = args.additional_output.take();
            args.additional_output = Some(count_simulation_output);
        }

        // RDKit✔️✔️:   auto tempResult = getFingerprintHelper(mol, args, effectiveSize);
        // RDKit✔️✔️:   auto result = std::make_unique<SparseBitVect>(resultSize);
        let temp_result = self.getFingerprintHelper(molecule, args, u64::from(effective_size))?;
        let mut result = SparseBitFingerprint::new(u64::from(result_size));

        // RDKit✔️✔️:   for (auto val : tempResult->getNonzeroElements()) {
        for (&bit_id, &count) in temp_result.nonzero_elements() {
            // RDKit✔️✔️:     if (dp_fingerprintArguments->df_countSimulation) {
            if count_simulation {
                for (idx, &bound) in count_bounds.iter().enumerate() {
                    // RDKit✔️✔️:       if (val.second >= static_cast<int>(bounds_count[i])) {
                    // RDKit✔️✔️:         OutputType nBitId = val.first * bounds_count.size() + i;
                    // RDKit✔️✔️:         result->setBit(nBitId);
                    // RDKit✔️✔️:         if (args.additionalOutput) {
                    // RDKit✔️✔️:           duplicateAdditionalOutputBit(*args.additionalOutput, *origAO,
                    // RDKit✔️✔️:                                          static_cast<OutputType>(val.first),
                    // RDKit✔️✔️:                                          nBitId);
                    // RDKit✔️✔️:         }
                    // RDKit✔️✔️:       }
                    if count >= bound as i32 {
                        let new_bit_id = bit_id * count_bounds.len() as u64 + idx as u64;
                        result.set_bit(new_bit_id);
                        if let (Some(temp_output), Some(orig_output)) = (
                            args.additional_output.as_ref(),
                            original_additional_output.as_mut(),
                        ) {
                            duplicate_additional_output_bit(
                                temp_output,
                                orig_output,
                                bit_id,
                                new_bit_id,
                            )?;
                        }
                    }
                }
            } else {
                // RDKit✔️✔️:     } else {
                // RDKit✔️✔️:       result->setBit(val.first);
                // RDKit✔️✔️:     }
                result.set_bit(bit_id);
            }
        }

        // RDKit✔️✔️:   if (origAO) {
        // RDKit✔️✔️:     if (origAO->atomCounts) {
        // RDKit✔️✔️:       *origAO->atomCounts = *countSimulationOutput.atomCounts;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     args.additionalOutput = origAO;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return result;
        // RDKit✔️✔️: }
        if let Some(mut orig_output) = original_additional_output {
            if orig_output.atom_counts.is_some() {
                orig_output.atom_counts = args
                    .additional_output
                    .as_ref()
                    .and_then(|output| output.atom_counts.clone());
            }
            args.additional_output = Some(orig_output);
        }
        Ok(result)
    }

    #[allow(non_snake_case)]
    pub fn getCountFingerprint(
        &self,
        molecule: &Molecule,
        args: &mut FingerprintFuncArguments,
    ) -> Result<SparseCountFingerprint, FingerprintError> {
        // RDKit source: FingerprintGenerator.cpp lines 577-590
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: std::unique_ptr<SparseIntVect<std::uint32_t>>
        // RDKit✔️✔️: FingerprintGenerator<OutputType>::getCountFingerprint(
        // RDKit✔️✔️:     const ROMol &mol, FingerprintFuncArguments &args) const {
        // RDKit✔️✔️:   auto tempResult =
        // RDKit✔️✔️:       getFingerprintHelper(mol, args, dp_fingerprintArguments->d_fpSize);
        let fp_size = self.fingerprint_arguments.fingerprint_arguments.d_fp_size;
        let temp_result = self.getFingerprintHelper(molecule, args, u64::from(fp_size))?;

        // RDKit✔️✔️:   auto result = std::make_unique<SparseIntVect<std::uint32_t>>(
        // RDKit✔️✔️:       dp_fingerprintArguments->d_fpSize);
        // RDKit✔️✔️:   for (auto val : tempResult->getNonzeroElements()) {
        // RDKit✔️✔️:     result->setVal(val.first, val.second);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return result;
        // RDKit✔️✔️: }
        let mut result = SparseCountFingerprint::new(u64::from(fp_size));
        for (&bit_id, &count) in temp_result.nonzero_elements() {
            result.set_val(bit_id, count);
        }
        Ok(result)
    }

    #[allow(non_snake_case)]
    pub fn getFingerprint(
        &self,
        molecule: &Molecule,
        args: &mut FingerprintFuncArguments,
    ) -> Result<Fingerprint, FingerprintError> {
        // RDKit source: FingerprintGenerator.cpp lines 592-650
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: std::unique_ptr<ExplicitBitVect>
        // RDKit✔️✔️: FingerprintGenerator<OutputType>::getFingerprint(
        // RDKit✔️✔️:     const ROMol &mol, FingerprintFuncArguments &args) const {
        // RDKit✔️✔️:   std::uint32_t effectiveSize = dp_fingerprintArguments->d_fpSize;
        let fp_args = &self.fingerprint_arguments.fingerprint_arguments;
        let mut effective_size = fp_args.d_fp_size;

        // RDKit✔️✔️:   if (dp_fingerprintArguments->df_countSimulation) {
        // RDKit✔️✔️:     if (dp_fingerprintArguments->d_countBounds.empty()) {
        // RDKit✔️✔️:       throw ValueErrorException("Count bounds are empty");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (dp_fingerprintArguments->d_countBounds.size() >= effectiveSize) {
        // RDKit✔️✔️:       throw ValueErrorException("Count bounds size is >= fingerprint size");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     effectiveSize /= dp_fingerprintArguments->d_countBounds.size();
        // RDKit✔️✔️:   }
        if fp_args.df_count_simulation {
            if fp_args.d_count_bounds.is_empty() {
                return Err(FingerprintError::InvalidArguments {
                    reason: "Count bounds are empty",
                });
            }
            if fp_args.d_count_bounds.len() >= effective_size as usize {
                return Err(FingerprintError::InvalidArguments {
                    reason: "Count bounds size is >= fingerprint size",
                });
            }
            effective_size /= fp_args.d_count_bounds.len() as u32;
        }

        // RDKit✔️✔️:   AdditionalOutput countSimulationOutput;
        // RDKit✔️✔️:   AdditionalOutput *origAO = nullptr;
        // RDKit✔️✔️:   if (dp_fingerprintArguments->df_countSimulation && args.additionalOutput) {
        // RDKit✔️✔️:     setupTempAdditionalOutput(args, countSimulationOutput, mol.getNumAtoms());
        // RDKit✔️✔️:     origAO = args.additionalOutput;
        // RDKit✔️✔️:     args.additionalOutput = &countSimulationOutput;
        // RDKit✔️✔️:   }
        let mut original_additional_output = None;
        if fp_args.df_count_simulation && args.additional_output.is_some() {
            let mut count_simulation_output = AdditionalOutput::default();
            setup_temp_additional_output(args, &mut count_simulation_output, molecule.num_atoms());
            original_additional_output = args.additional_output.take();
            args.additional_output = Some(count_simulation_output);
        }

        // RDKit✔️✔️:   auto tempResult = getFingerprintHelper(mol, args, effectiveSize);
        // RDKit✔️✔️:   auto result =
        // RDKit✔️✔️:       std::make_unique<ExplicitBitVect>(dp_fingerprintArguments->d_fpSize);
        let temp_result = self.getFingerprintHelper(molecule, args, u64::from(effective_size))?;
        let mut on_bits = Vec::new();

        // RDKit✔️✔️:   for (auto val : tempResult->getNonzeroElements()) {
        for (&bit_id, &count) in temp_result.nonzero_elements() {
            // RDKit✔️✔️:     if (dp_fingerprintArguments->df_countSimulation) {
            if fp_args.df_count_simulation {
                for (idx, &bound) in fp_args.d_count_bounds.iter().enumerate() {
                    // RDKit✔️✔️:       if (val.second >= static_cast<int>(bounds_count[i])) {
                    // RDKit✔️✔️:         OutputType nBitId = val.first * bounds_count.size() + i;
                    // RDKit✔️✔️:         result->setBit(nBitId);
                    // RDKit✔️✔️:         if (args.additionalOutput) {
                    // RDKit✔️✔️:           duplicateAdditionalOutputBit(*args.additionalOutput, *origAO,
                    // RDKit✔️✔️:                                          static_cast<OutputType>(val.first),
                    // RDKit✔️✔️:                                          nBitId);
                    // RDKit✔️✔️:         }
                    // RDKit✔️✔️:       }
                    if count >= bound as i32 {
                        let new_bit_id = bit_id * fp_args.d_count_bounds.len() as u64 + idx as u64;
                        on_bits.push(new_bit_id as usize);
                        if let (Some(temp_output), Some(orig_output)) = (
                            args.additional_output.as_ref(),
                            original_additional_output.as_mut(),
                        ) {
                            duplicate_additional_output_bit(
                                temp_output,
                                orig_output,
                                bit_id,
                                new_bit_id,
                            )?;
                        }
                    }
                }
            } else {
                // RDKit✔️✔️:     } else {
                // RDKit✔️✔️:       result->setBit(val.first);
                // RDKit✔️✔️:     }
                on_bits.push(bit_id as usize);
            }
        }

        // RDKit✔️✔️:   if (origAO) {
        // RDKit✔️✔️:     if (origAO->atomCounts) {
        // RDKit✔️✔️:       *origAO->atomCounts = *countSimulationOutput.atomCounts;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     args.additionalOutput = origAO;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return result;
        // RDKit✔️✔️: }
        if let Some(mut orig_output) = original_additional_output {
            if orig_output.atom_counts.is_some() {
                orig_output.atom_counts = args
                    .additional_output
                    .as_ref()
                    .and_then(|output| output.atom_counts.clone());
            }
            args.additional_output = Some(orig_output);
        }

        Ok(Fingerprint::from_on_bits(
            fp_args.d_fp_size as usize,
            on_bits,
        ))
    }

    #[allow(non_snake_case)]
    pub fn getFingerprintHelper(
        &self,
        molecule: &Molecule,
        args: &mut FingerprintFuncArguments,
        fp_size: u64,
    ) -> Result<SparseCountFingerprint, FingerprintError> {
        // RDKit source: FingerprintGenerator.cpp lines 325-435
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: std::unique_ptr<SparseIntVect<OutputType>>
        // RDKit✔️✔️: FingerprintGenerator<OutputType>::getFingerprintHelper(
        // RDKit✔️✔️:     const ROMol &mol, FingerprintFuncArguments &args,
        // RDKit✔️✔️:     const std::uint64_t fpSize) const {
        // RDKit❌❌:   const ROMol *lmol = &mol;
        // RDKit❌❌:   std::unique_ptr<ROMol> tmol;
        // RDKit❌❌:   if (dp_fingerprintArguments->df_includeChirality &&
        // RDKit❌❌:       !mol.hasProp(common_properties::_StereochemDone)) {
        // RDKit❌❌:     tmol = std::unique_ptr<ROMol>(new ROMol(mol));
        // RDKit❌❌:     MolOps::assignStereochemistry(*tmol);
        // RDKit❌❌:     lmol = tmol.get();
        // RDKit❌❌:   }
        if self
            .fingerprint_arguments
            .fingerprint_arguments
            .df_include_chirality
            && molecule.prop("_StereochemDone").is_none()
        {
            return Err(FingerprintError::UnsupportedOption {
                option: "includeChirality",
                reason: "FingerprintGenerator requires RDKit assignStereochemistry parity before recomputing missing stereochemistry labels",
            });
        }

        let labeled;
        let molecule = if self
            .fingerprint_arguments
            .fingerprint_arguments
            .df_include_chirality
            && !rdkit_use_legacy_stereo_perception()
            && molecule.prop("_CIPComputed").is_none()
        {
            labeled = assign_cip_labels(molecule, 0)?;
            &labeled
        } else {
            molecule
        };

        // RDKit✔️✔️:   if (args.additionalOutput) {
        // RDKit✔️✔️:     reinitAdditionalOutput(*args.additionalOutput, mol.getNumAtoms());
        // RDKit✔️✔️:   }
        if let Some(additional_output) = args.additional_output.as_mut() {
            additional_output.reset_for_atom_count(molecule.num_atoms());
        }

        // RDKit✔️✔️:   bool hashResults = false;
        // RDKit✔️✔️:   if (fpSize != 0) {
        // RDKit✔️✔️:     hashResults = true;
        // RDKit✔️✔️:   }
        let hash_results = fp_size != 0;

        // RDKit✔️✔️:   std::unique_ptr<std::vector<std::uint32_t>> atomInvariants = nullptr;
        // RDKit✔️✔️:   if (args.customAtomInvariants) {
        // RDKit✔️✔️:     atomInvariants.reset(
        // RDKit✔️✔️:         new std::vector<std::uint32_t>(*args.customAtomInvariants));
        // RDKit✔️✔️:   } else if (dp_atomInvariantsGenerator) {
        // RDKit✔️✔️:     atomInvariants.reset(dp_atomInvariantsGenerator->getAtomInvariants(mol));
        // RDKit✔️✔️:   }
        let atom_invariants = if let Some(custom) = &args.custom_atom_invariants {
            custom.clone()
        } else {
            self.atom_invariants_generator
                .get_atom_invariants(molecule)?
        };

        // RDKit✔️✔️:   std::unique_ptr<std::vector<std::uint32_t>> bondInvariants = nullptr;
        // RDKit✔️✔️:   if (args.customBondInvariants) {
        // RDKit✔️✔️:     bondInvariants.reset(
        // RDKit✔️✔️:         new std::vector<std::uint32_t>(*args.customBondInvariants));
        // RDKit✔️✔️:   } else if (dp_bondInvariantsGenerator) {
        // RDKit✔️✔️:     bondInvariants.reset(dp_bondInvariantsGenerator->getBondInvariants(mol));
        // RDKit✔️✔️:   }
        let bond_invariants = if let Some(custom) = &args.custom_bond_invariants {
            custom.clone()
        } else {
            self.bond_invariants_generator
                .try_get_bond_invariants(molecule)?
        };

        // RDKit✔️✔️:   auto atomEnvironments = dp_atomEnvironmentGenerator->getEnvironments(
        // RDKit✔️✔️:       *lmol, dp_fingerprintArguments, args.fromAtoms, args.ignoreAtoms,
        // RDKit✔️✔️:       args.confId, args.additionalOutput, atomInvariants.get(),
        // RDKit✔️✔️:       bondInvariants.get(), hashResults);
        let atom_environments = self.atom_environment_generator.getEnvironments(
            molecule,
            &self.fingerprint_arguments,
            args.from_atoms.as_deref(),
            args.ignore_atoms.as_deref(),
            &atom_invariants,
            &bond_invariants,
        )?;

        // RDKit✔️✔️:   auto res = std::make_unique<SparseIntVect<OutputType>>(
        // RDKit✔️✔️:       fpSize ? fpSize : dp_atomEnvironmentGenerator->getResultSize());
        let result_size = if fp_size == 0 {
            self.atom_environment_generator.getResultSize()
        } else {
            fp_size
        };
        let mut result = SparseCountFingerprint::new(result_size);

        // RDKit✔️✔️:   typedef boost::random::mersenne_twister<std::uint32_t, 32, 4, 2, 31,
        // RDKit✔️✔️:                                           0x9908b0df, 11, 7, 0x9d2c5680, 15,
        // RDKit✔️✔️:                                           0xefc60000, 18, 3346425566U>
        // RDKit✔️✔️:       rng_type;
        // RDKit✔️✔️:   typedef boost::uniform_int<> distrib_type;
        // RDKit✔️✔️:   typedef boost::variate_generator<rng_type &, distrib_type> source_type;
        // RDKit✔️✔️:   std::unique_ptr<rng_type> generator;
        // RDKit✔️✔️:   std::unique_ptr<distrib_type> dist;
        // RDKit✔️✔️:   std::unique_ptr<source_type> randomSource;
        // RDKit✔️✔️:   if (dp_fingerprintArguments->d_numBitsPerFeature > 1) {
        // RDKit✔️✔️:     generator.reset(new rng_type(42u));
        // RDKit✔️✔️:     dist.reset(new distrib_type(0, INT_MAX));
        // RDKit✔️✔️:     randomSource.reset(new source_type(*generator, *dist));
        // RDKit✔️✔️:   }
        let mut random_source = (self
            .fingerprint_arguments
            .fingerprint_arguments
            .d_num_bits_per_feature
            > 1)
        .then(|| RdkitFingerprintMtRng::new(42));

        // RDKit✔️✔️:   for (const auto env : atomEnvironments) {
        // RDKit✔️✔️:     OutputType seed = env->getBitId(dp_fingerprintArguments,
        // RDKit✔️✔️:                                     atomInvariants.get(), bondInvariants.get(),
        // RDKit✔️✔️:                                     args.additionalOutput, hashResults, fpSize);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     auto bitId = seed;
        // RDKit✔️✔️:     if (fpSize != 0) {
        // RDKit✔️✔️:       bitId %= fpSize;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     res->setVal(bitId, res->getVal(bitId) + 1);
        // RDKit✔️✔️:     if (args.additionalOutput) {
        // RDKit✔️✔️:       env->updateAdditionalOutput(args.additionalOutput, bitId);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (dp_fingerprintArguments->d_numBitsPerFeature > 1) {
        // RDKit✔️✔️:       generator->seed(static_cast<rng_type::result_type>(seed));
        // RDKit✔️✔️:
        // RDKit✔️✔️:       for (boost::uint32_t bitN = 1;
        // RDKit✔️✔️:            bitN < dp_fingerprintArguments->d_numBitsPerFeature; ++bitN) {
        // RDKit✔️✔️:         bitId = (*randomSource)();
        // RDKit✔️✔️:         if (fpSize != 0) {
        // RDKit✔️✔️:           bitId %= fpSize;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:         res->setVal(bitId, res->getVal(bitId) + 1);
        // RDKit✔️✔️:         if (args.additionalOutput) {
        // RDKit✔️✔️:           env->updateAdditionalOutput(args.additionalOutput, bitId);
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     delete env;
        // RDKit✔️✔️:   }
        for env in atom_environments {
            let seed = env.getBitId();
            let bit_id = if hash_results { seed % fp_size } else { seed };
            result.set_val(bit_id, result.get_val(bit_id) + 1);
            if let Some(additional_output) = args.additional_output.as_mut() {
                env.updateAdditionalOutput(additional_output, bit_id);
            }
            if let Some(random_source) = random_source.as_mut() {
                random_source.seed(seed as u32);
                for _ in 1..self
                    .fingerprint_arguments
                    .fingerprint_arguments
                    .d_num_bits_per_feature
                {
                    let random_bit_id = u64::from(random_source.uniform_int_0_to_i32_max());
                    let random_bit_id = if hash_results {
                        random_bit_id % fp_size
                    } else {
                        random_bit_id
                    };
                    result.set_val(random_bit_id, result.get_val(random_bit_id) + 1);
                    if let Some(additional_output) = args.additional_output.as_mut() {
                        env.updateAdditionalOutput(additional_output, random_bit_id);
                    }
                }
            }
        }

        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        Ok(result)
    }
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
    #[error("Morgan fingerprint requires n_bits > 0")]
    EmptyFingerprint,
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
    #[error("invalid fingerprint arguments: {reason}")]
    InvalidArguments { reason: &'static str },
    #[error("invalid fingerprint arguments JSON: {0}")]
    InvalidArgumentsJson(String),
    #[error("CIPLabeler failed while preparing Morgan fingerprint chirality: {reason}")]
    CipLabeler { reason: String },
    #[error("invalid SMARTS pattern '{pattern}': {reason}")]
    InvalidSmartsPattern { pattern: String, reason: String },
    #[error("unsupported fingerprint option {option}: {reason}")]
    UnsupportedOption {
        option: &'static str,
        reason: &'static str,
    },
    #[error("fingerprint bit length mismatch: {left} != {right}")]
    BitLengthMismatch { left: usize, right: usize },
}

impl From<crate::chemistry::ciplabeler::CipLabelerError> for FingerprintError {
    fn from(error: crate::chemistry::ciplabeler::CipLabelerError) -> Self {
        Self::CipLabeler {
            reason: error.to_string(),
        }
    }
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
    let generator = getMorganGeneratorWithParams(
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

fn compute_initial_invariants(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    params: &MorganFingerprintParams,
) -> Result<Vec<u32>, FingerprintError> {
    let invariants = match &params.atom_invariants_generator {
        MorganAtomInvariantsGenerator::Connectivity { .. } => {
            compute_connectivity_invariants(molecule, adjacency, params)
        }
        MorganAtomInvariantsGenerator::Feature => compute_feature_invariants(molecule)?,
    };

    if let Some(custom) = &params.custom_atom_invariants {
        let mut overridden = invariants;
        for (i, inv) in overridden.iter_mut().enumerate() {
            if let Some(c) = custom.get(i) {
                *inv = *c;
            }
        }
        Ok(overridden)
    } else {
        Ok(invariants)
    }
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

fn morgan_bond_invariant(
    bond_idx: usize,
    bond: &crate::Bond,
    params: &MorganFingerprintParams,
) -> u32 {
    // custom_bond_invariants override: if provided and this bond index has a
    // value, use it directly instead of computing from bond order.
    if let Some(custom) = &params.custom_bond_invariants {
        if let Some(&inv) = custom.get(bond_idx) {
            return inv;
        }
    }
    let use_bond_types = params
        .bond_invariants_generator
        .as_ref()
        .map_or(params.use_bond_types, |generator| generator.use_bond_types);
    let use_chirality = params
        .bond_invariants_generator
        .as_ref()
        .map_or(params.use_chirality, |generator| generator.use_chirality);
    MorganBondInvGenerator::new(use_bond_types, use_chirality).bond_invariant(bond)
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

// RDKit uses:  gboost::hash_combine(seed, value)
// which expands to:  seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
pub(crate) fn hash_combine(seed: &mut u32, value: u32) {
    *seed ^= value
        .wrapping_add(0x9e3779b9u32)
        .wrapping_add(seed.wrapping_shl(6))
        .wrapping_add(seed.wrapping_shr(2));
}

// ---------------------------------------------------------------------------
// Fingerprint construction
// ---------------------------------------------------------------------------

fn fold_invariant(invariant: u32, n_bits: usize) -> usize {
    invariant as usize % n_bits
}

fn build_fingerprint(
    molecule: &Molecule,
    all_rounds: &[Vec<u32>],
    params: &MorganFingerprintParams,
) -> Result<MorganFingerprintOutput, FingerprintError> {
    let n_bits = params.n_bits;
    let collect = params.collect_additional_output;
    let mut additional_output = if collect {
        let mut output = AdditionalOutput::new();
        output.allocate_atom_counts();
        output.allocate_atom_to_bits();
        output.allocate_bit_info_map();
        output.allocate_atoms_per_bit();
        output.reset_for_atom_count(molecule.num_atoms());
        Some(output)
    } else {
        None
    };

    let mut on_bits = Vec::new();

    // Track seen invariants for include_redundant_environments=false.
    // When disabled, each unique invariant value is only allowed to
    // contribute bits once across all atoms and all rounds, matching
    // RDKit's behavior of deduplicating redundant environments.
    let mut seen_invariants: std::collections::HashSet<u32> = std::collections::HashSet::new();

    for (round_idx, round_invs) in all_rounds.iter().enumerate() {
        let round = round_idx as u32;
        for atom_idx in 0..molecule.num_atoms() {
            if atom_is_excluded(atom_idx, params) {
                continue;
            }
            let inv = round_invs[atom_idx];

            // When include_redundant_environments is false, skip invariants
            // that have already contributed bits in a previous round or atom.
            if !params.include_redundant_environments && !seen_invariants.insert(inv) {
                continue;
            }

            if params.only_nonzero_invariants && inv == 0 {
                continue;
            }

            if params.count_simulation {
                let env = MorganAtomEnv::new(u64::from(inv), atom_idx, round, molecule);
                let bit = fold_invariant(env.getBitId() as u32, n_bits);
                on_bits.push(bit);

                if let Some(output) = additional_output.as_mut() {
                    env.updateAdditionalOutput(output, bit as u64);
                }
            } else {
                for chunk in 0..params.num_bits_per_feature {
                    let code = inv.wrapping_add(chunk.wrapping_mul(0x517cc1b7));
                    let env = MorganAtomEnv::new(u64::from(code), atom_idx, round, molecule);
                    let bit = fold_invariant(env.getBitId() as u32, n_bits);
                    on_bits.push(bit);

                    if let Some(output) = additional_output.as_mut() {
                        env.updateAdditionalOutput(output, bit as u64);
                    }
                }
            }
        }
    }

    // Count-simulation folding: each unique bit's count is compared against
    // count_bounds to set additional offset bits.
    if params.count_simulation && !params.count_bounds.is_empty() {
        let mut counts_per_bit: BTreeMap<usize, u32> = BTreeMap::new();
        for &bit in &on_bits {
            *counts_per_bit.entry(bit).or_insert(0) += 1;
        }
        on_bits.clear();
        for (&bit, &count) in &counts_per_bit {
            on_bits.push(bit);
            for (bound_idx, &bound) in params.count_bounds.iter().enumerate().skip(1) {
                if count >= bound {
                    let offset_bit = (bit + bound_idx * n_bits) % n_bits;
                    on_bits.push(offset_bit);
                }
            }
        }
    }

    let fingerprint = Fingerprint::from_on_bits(n_bits, on_bits.iter().copied());

    let additional_output = additional_output.map(morgan_additional_output_from_rdkit_output);

    Ok(MorganFingerprintOutput {
        fingerprint,
        additional_output,
    })
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

fn atom_is_excluded(index: usize, params: &MorganFingerprintParams) -> bool {
    if let Some(from) = &params.from_atoms {
        return !from.contains(&index);
    }
    if let Some(ignore) = &params.ignore_atoms {
        return ignore.contains(&index);
    }
    false
}

// ---------------------------------------------------------------------------
// Topological (Path-Based) Fingerprint
// RDKit source: GraphMol/Fingerprints/Fingerprints.h
// RDKit source: GraphMol/Fingerprints/FingerprintUtil.cpp
// ---------------------------------------------------------------------------

/// Parameters reserved for the future source-backed `RDKFingerprintMol` port.
///
/// The current public shape is retained for API planning only. The operation
/// returns `FingerprintError::UnsupportedOption` until the complete RDKit
/// generator and option surface is ported.
///
/// # Parameters
/// - `min_path`: minimum path length in bonds (default 1).
/// - `max_path`: maximum path length in bonds (default 7).
/// - `n_bits`: size of the output fingerprint (default 2048).
/// - `n_bits_per_hash`: number of bit positions to set per path hash (default 2).
/// - `use_bond_types`: maps to RDKit's `useBondOrder` option.
/// - `from_atoms`: if `Some`, only enumerate paths starting from these atoms.
/// - `ignore_atoms`: reserved compatibility parameter; RDKit's exposed
///   `RDKFingerprintMol` signature has no corresponding option.
#[derive(Debug, Clone)]
pub struct TopologicalFingerprintParams {
    pub min_path: u32,
    pub max_path: u32,
    pub n_bits: usize,
    pub n_bits_per_hash: u32,
    pub use_bond_types: bool,
    pub from_atoms: Option<Vec<usize>>,
    pub ignore_atoms: Option<Vec<usize>>,
}

impl Default for TopologicalFingerprintParams {
    fn default() -> Self {
        Self {
            min_path: 1,
            max_path: 7,
            n_bits: 2048,
            n_bits_per_hash: 2,
            use_bond_types: true,
            from_atoms: None,
            ignore_atoms: None,
        }
    }
}

// RDKit source: GraphMol/Fingerprints/Fingerprints.h::RDKFingerprintMol
// RDKit❌❌: RDKIT_FINGERPRINTS_EXPORT ExplicitBitVect *RDKFingerprintMol(
// RDKit❌❌:     const ROMol &mol, unsigned int minPath = 1, unsigned int maxPath = 7,
// RDKit❌❌:     unsigned int fpSize = 2048, unsigned int nBitsPerHash = 2,
// RDKit❌❌:     bool useHs = true, double tgtDensity = 0.0, unsigned int minSize = 128,
// RDKit❌❌:     bool branchedPaths = true, bool useBondOrder = true,
// RDKit❌❌:     std::vector<std::uint32_t> *atomInvariants = nullptr,
// RDKit❌❌:     const std::vector<std::uint32_t> *fromAtoms = nullptr,
// RDKit❌❌:     std::vector<std::vector<std::uint32_t>> *atomBits = nullptr,
// RDKit❌❌:     std::map<std::uint32_t, std::vector<std::vector<int>>> *bitInfo = nullptr);
// RDKit source: GraphMol/Fingerprints/Fingerprints.cpp::RDKFingerprintMol
// RDKit❌❌: std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
// RDKit❌❌:     RDKit::RDKitFP::getRDKitFPGenerator<std::uint32_t>(
// RDKit❌❌:         minPath, maxPath, useHs, branchedPaths, useBondOrder));
// RDKit❌❌: fpgen->getOptions()->d_fpSize = fpSize;
// RDKit❌❌: fpgen->getOptions()->d_numBitsPerFeature = nBitsPerHash;
//
// No approximate bit vector is returned. The previous local DFS/hash
// implementation was removed because it was not RDKit RDKFingerprint behavior.
pub fn topological_fingerprint(
    _molecule: &Molecule,
    _params: &TopologicalFingerprintParams,
) -> Result<Fingerprint, FingerprintError> {
    Err(FingerprintError::UnsupportedOption {
        option: "topological_fingerprint",
        reason: "RDKit RDKFingerprint exact-bit source port is not implemented",
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
mod tests {
    use super::*;
    use crate::{AtomSpec, BondSpec, BondStereo, Element, Molecule};

    fn default_morgan_params(radius: u32, n_bits: usize) -> MorganFingerprintParams {
        MorganFingerprintParams {
            radius,
            n_bits,
            ..Default::default()
        }
    }

    #[test]
    fn unported_topological_fingerprint_fails_closed() {
        let err =
            topological_fingerprint(&Molecule::new(), &TopologicalFingerprintParams::default())
                .expect_err("unfinished RDKFingerprint must not return approximate bits");
        assert!(matches!(
            err,
            FingerprintError::UnsupportedOption {
                option: "topological_fingerprint",
                ..
            }
        ));
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
    fn morgan_fingerprint_include_chirality_without_stereochem_done_fails_closed() {
        // RDKit FingerprintGenerator first calls MolOps::assignStereochemistry()
        // when includeChirality=true and _StereochemDone is absent. That
        // branch remains explicitly unsupported until the RDKit source path is
        // ported; do not silently substitute legacy/partial stereo labels.
        let r_mol = Molecule::from_smiles_with_sanitize("C[C@@H](O)CC", false).unwrap();
        let params = MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            use_chirality: true,
            ..Default::default()
        };
        let err = morgan_fingerprint(&r_mol, &params).unwrap_err();
        assert!(matches!(
            err,
            FingerprintError::UnsupportedOption {
                option: "includeChirality",
                ..
            }
        ));
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
        assert_eq!(ao.atoms_per_bit, Some(BTreeMap::new()));
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
    fn morgan_fingerprint_generator_outputs_radius_and_fp_size_match_rdkit_golden() {
        let mol = rdkit_morgan_oracle_mol("CCC");
        let generator = getMorganGeneratorWithParams(
            2,
            false,
            false,
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
            r#"{"countSimulation":true,"fpSize":4096,"numBitsPerFeature":3,"includeChirality":true,"countBounds":[1,2,4,8]}"#
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
            r#"{"type":"MorganArguments","onlyNonzeroInvariants":true,"radius":2,"countSimulation":true,"fpSize":1024,"numBitsPerFeature":1,"includeChirality":true,"countBounds":[1,2,4,8]}"#
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
    fn default_feature_smarts_matchers_parse_source_patterns() {
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
}
