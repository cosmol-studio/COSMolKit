//! RDKFingerprint over the sole detached graph values. Pinned RDKit351f8f.
use crate::generator::{
    FingerprintArguments, FingerprintEnvironment, accumulate_sparse_counts_into,
    project_fingerprint,
};
use crate::hash::{hash_combine, hash_range};
use crate::layered::{LayeredGraphInput as GraphInput, enumerate_fingerprint_paths};
use crate::{Fingerprint, FingerprintAdditionalOutput, FingerprintError, SparseCountFingerprint};
use cosmolkit_model::{Bond, QueryGraph, TopologyBlock};
use cosmolkit_search::{
    is_complex_atom_query, is_complex_bond_query, is_complex_concrete_bond_query,
};
use std::{borrow::Cow, collections::BTreeMap};
#[derive(Debug)]
pub enum TopologicalFingerprintError {
    InvalidArguments { reason: &'static str },
    OutputNotRequested { field: &'static str },
    Topology(cosmolkit_model::TopologyValidationError),
    Query(cosmolkit_model::QueryGraphError),
    Paths(crate::LayeredFingerprintError),
    Value(FingerprintError),
}
impl std::fmt::Display for TopologicalFingerprintError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::OutputNotRequested { field } => write!(
                f,
                "topological {field} output was not requested for this fingerprint result"
            ),
            Self::InvalidArguments { reason } => f.write_str(reason),
            Self::Topology(e) => e.fmt(f),
            Self::Query(e) => e.fmt(f),
            Self::Paths(e) => e.fmt(f),
            Self::Value(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for TopologicalFingerprintError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidArguments { .. } | Self::OutputNotRequested { .. } => None,
            Self::Topology(e) => Some(e),
            Self::Query(e) => Some(e),
            Self::Paths(e) => Some(e),
            Self::Value(e) => Some(e),
        }
    }
}
impl From<FingerprintError> for TopologicalFingerprintError {
    fn from(e: FingerprintError) -> Self {
        Self::Value(e)
    }
}
impl From<crate::LayeredFingerprintError> for TopologicalFingerprintError {
    fn from(e: crate::LayeredFingerprintError) -> Self {
        Self::Paths(e)
    }
}
impl From<cosmolkit_model::TopologyValidationError> for TopologicalFingerprintError {
    fn from(e: cosmolkit_model::TopologyValidationError) -> Self {
        Self::Topology(e)
    }
}
impl From<cosmolkit_model::QueryGraphError> for TopologicalFingerprintError {
    fn from(e: cosmolkit_model::QueryGraphError) -> Self {
        Self::Query(e)
    }
}
#[derive(Debug, Clone, PartialEq)]
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
    pub ignore_atoms: Option<Vec<u32>>,
}

impl Default for TopologicalFingerprintParams {
    fn default() -> Self {
        // RDKit❗✔️: RDKIT_FINGERPRINTS_EXPORT ExplicitBitVect *RDKFingerprintMol(
        // RDKit❗✔️:     const ROMol &mol, unsigned int minPath = 1, unsigned int maxPath = 7,
        // RDKit❗✔️:     unsigned int fpSize = 2048, unsigned int nBitsPerHash = 2,
        // RDKit❗✔️:     bool useHs = true, double tgtDensity = 0.0, unsigned int minSize = 128,
        // RDKit❗✔️:     bool branchedPaths = true, bool useBondOrder = true,
        // RDKit❗✔️:     std::vector<std::uint32_t> *atomInvariants = nullptr,
        // RDKit❗✔️:     const std::vector<std::uint32_t> *fromAtoms = nullptr,
        // RDKit❗✔️:     std::vector<std::vector<std::uint32_t>> *atomBits = nullptr,
        // RDKit❗✔️:     std::map<std::uint32_t, std::vector<std::vector<int>>> *bitInfo = nullptr);
        // Complexity: constant defaults, no runtime allocation or chemistry work.
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
            ignore_atoms: None,
        }
    }
}

impl TopologicalFingerprintParams {
    pub fn validate(&self) -> Result<(), TopologicalFingerprintError> {
        // RDKit❗✔️:   PRECONDITION(minPath != 0, "minPath==0");
        // RDKit❗✔️:   PRECONDITION(maxPath >= minPath, "maxPath<minPath");
        // RDKit❗✔️:   PRECONDITION(fpSize != 0, "fpSize==0");
        // RDKit❗✔️:   PRECONDITION(nBitsPerHash != 0, "nBitsPerHash==0");
        // Complexity: constant source argument checks; retained original Rust
        // finite/nonnegative-density validation is a separate public boundary.
        if self.min_path == 0 {
            return Err(TopologicalFingerprintError::InvalidArguments {
                reason: "minPath==0",
            });
        }
        if self.max_path < self.min_path {
            return Err(TopologicalFingerprintError::InvalidArguments {
                reason: "maxPath<minPath",
            });
        }
        if self.fp_size == 0 {
            return Err(TopologicalFingerprintError::InvalidArguments {
                reason: "fpSize==0",
            });
        }
        if self.num_bits_per_feature == 0 {
            return Err(TopologicalFingerprintError::InvalidArguments {
                reason: "nBitsPerHash==0",
            });
        }
        if !self.target_density.is_finite() || self.target_density < 0.0 {
            return Err(TopologicalFingerprintError::InvalidArguments {
                reason: "tgtDensity must be finite and non-negative",
            });
        }
        Ok(())
    }
}

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

impl TopologicalFingerprintResult {
    /// Borrow the calculated value without regenerating the fingerprint.
    pub fn fingerprint(&self) -> &Fingerprint {
        &self.fingerprint
    }
    /// Source atomBits ordering and original pre-fold bit identifiers.
    pub fn atom_bits(&self) -> Result<&[Vec<u32>], TopologicalFingerprintError> {
        self.output
            .atom_bits
            .as_deref()
            .ok_or(TopologicalFingerprintError::OutputNotRequested { field: "atom_bits" })
    }
    /// Source bitInfo ordering and original pre-fold bit identifiers.
    pub fn bit_info(&self) -> Result<&BTreeMap<u32, Vec<Vec<i32>>>, TopologicalFingerprintError> {
        self.output
            .bit_info
            .as_ref()
            .ok_or(TopologicalFingerprintError::OutputNotRequested { field: "bit_info" })
    }
}
impl std::fmt::Display for TopologicalFingerprintResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "TopologicalFingerprintResult(n_bits={}, has_atom_bits={}, has_bit_info={})",
            self.fingerprint.n_bits(),
            self.output.atom_bits.is_some(),
            self.output.bit_info.is_some()
        )
    }
}

fn atom_count(graph: GraphInput<'_>) -> usize {
    match graph {
        GraphInput::Concrete(t) => t.atoms.len(),
        GraphInput::Query(q) => q.num_atoms(),
    }
}
fn bonds(graph: GraphInput<'_>) -> usize {
    match graph {
        GraphInput::Concrete(t) => t.bonds.len(),
        GraphInput::Query(q) => q.num_bonds(),
    }
}
fn bond(graph: GraphInput<'_>, index: usize) -> &Bond {
    match graph {
        GraphInput::Concrete(t) => &t.bonds[index],
        GraphInput::Query(q) => q.bonds()[index].bond(),
    }
}
fn validate_graph(graph: GraphInput<'_>) -> Result<(), TopologicalFingerprintError> {
    match graph {
        GraphInput::Concrete(t) => t.validate().map_err(Into::into),
        GraphInput::Query(q) => q.validate().map_err(Into::into),
    }
}
fn rdkit_fp_atom_invariants(graph: GraphInput<'_>) -> Vec<u32> {
    // RDKit❗✔️: std::vector<std::uint32_t> *RDKitFPAtomInvGenerator::getAtomInvariants(
    // RDKit❗✔️:     const ROMol &mol) const {
    // RDKit❗✔️:   auto *result = new std::vector<std::uint32_t>();
    // RDKit❗✔️:   result->reserve(mol.getNumAtoms());
    // RDKit❗✔️:   for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
    // RDKit❗✔️:        atomIt != mol.endAtoms(); ++atomIt) {
    // RDKit❗✔️:     unsigned int aHash = ((*atomIt)->getAtomicNum() % 128) << 1 |
    // RDKit❗✔️:                          static_cast<unsigned int>((*atomIt)->getIsAromatic());
    // RDKit❗✔️:     result->push_back(aHash);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return result;
    // RDKit❗✔️: }

    // Complexity: one reserved vector and one O(A) pass, borrowing source identity.
    match graph {
        GraphInput::Concrete(t) => t
            .atoms
            .iter()
            .map(|a| ((u32::from(a.atomic_number()) % 128) << 1) | u32::from(a.is_aromatic()))
            .collect(),
        GraphInput::Query(q) => q
            .atoms()
            .iter()
            .map(|a| ((u32::from(a.atomic_number()) % 128) << 1) | u32::from(a.is_aromatic()))
            .collect(),
    }
}
fn identify_query_bonds(graph: GraphInput<'_>) -> Vec<u8> {
    // RDKit❗✔️: void identifyQueryBonds(const ROMol &mol, std::vector<const Bond *> &bondCache,
    // RDKit❗✔️:                         std::vector<short> &isQueryBond) {
    // RDKit❗✔️:   bondCache.resize(mol.getNumBonds());
    // RDKit❗✔️:   ROMol::EDGE_ITER firstB, lastB;
    // RDKit❗✔️:   boost::tie(firstB, lastB) = mol.getEdges();
    // RDKit❗✔️:   while (firstB != lastB) {
    // RDKit❗✔️:     const Bond *bond = mol[*firstB];
    // RDKit❗✔️:     isQueryBond[bond->getIdx()] = 0x0;
    // RDKit❗✔️:     bondCache[bond->getIdx()] = bond;
    // RDKit❗✔️:     if (isComplexQuery(bond)) {
    // RDKit❗✔️:       isQueryBond[bond->getIdx()] = 0x1;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (isComplexQuery(bond->getBeginAtom())) {
    // RDKit❗✔️:       isQueryBond[bond->getIdx()] |= 0x2;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (isComplexQuery(bond->getEndAtom())) {
    // RDKit❗✔️:       isQueryBond[bond->getIdx()] |= 0x4;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     ++firstB;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }

    // Complexity: O(B), one byte per bond (source short), no query AST clone.
    (0..bonds(graph))
        .map(|i| match graph {
            GraphInput::Concrete(t) => u8::from(is_complex_concrete_bond_query(&t.bonds[i])),
            GraphInput::Query(q) => {
                let b = &q.bonds()[i];
                u8::from(is_complex_bond_query(b))
                    | (u8::from(is_complex_atom_query(&q.atoms()[b.begin().index()])) << 1)
                    | (u8::from(is_complex_atom_query(&q.atoms()[b.end().index()])) << 2)
            }
        })
        .collect()
}
#[derive(Debug, Clone, PartialEq, Eq)]
struct RdkitFpBondHashInputs {
    atoms_in_path: Vec<bool>,
    bond_hashes: Vec<u32>,
}
fn rdkit_fp_generate_bond_hash_inputs(
    graph: GraphInput<'_>,
    path: &[usize],
    use_bond_order: bool,
    atom_invariants: &[u32],
    query_bonds: &[u8],
) -> Result<RdkitFpBondHashInputs, TopologicalFingerprintError> {
    // RDKit❗❌: std::vector<unsigned int> generateBondHashes(
    // RDKit❗❌:     const ROMol &mol, boost::dynamic_bitset<> &atomsInPath,
    // RDKit❗❌:     const std::vector<const Bond *> &bondCache,
    // RDKit❗❌:     const std::vector<short> &isQueryBond, const PATH_TYPE &path,
    // RDKit❗❌:     bool useBondOrder, const std::vector<std::uint32_t> *atomInvariants) {
    // RDKit❗❌:   PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit❗❌:                "bad atomInvariants size");
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<unsigned int> bondHashes;
    // RDKit❗❌:   atomsInPath.reset();
    // RDKit❗❌:   bool queryInPath = false;
    // RDKit❗❌:   std::vector<unsigned int> atomDegrees(mol.getNumAtoms(), 0);
    // RDKit❗❌:   for (unsigned int i = 0; i < path.size() && !queryInPath; ++i) {
    // RDKit❗❌:     const Bond *bi = bondCache[path[i]];
    // RDKit❗❌:     CHECK_INVARIANT(bi, "bond not in cache");
    // RDKit❗❌:     atomDegrees[bi->getBeginAtomIdx()]++;
    // RDKit❗❌:     atomDegrees[bi->getEndAtomIdx()]++;
    // RDKit❗❌:     atomsInPath.set(bi->getBeginAtomIdx());
    // RDKit❗❌:     atomsInPath.set(bi->getEndAtomIdx());
    // RDKit❗❌:     if (isQueryBond[path[i]]) {
    // RDKit❗❌:       queryInPath = true;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (queryInPath) {
    // RDKit❗❌:     return bondHashes;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // -----------------
    // RDKit❗❌:   // calculate the bond hashes:
    // RDKit❗❌:   std::vector<unsigned int> bondNbrs(path.size(), 0);
    // RDKit❗❌:   bondHashes.reserve(path.size() + 1);
    // RDKit❗❌:
    // RDKit❗❌:   for (unsigned int i = 0; i < path.size(); ++i) {
    // RDKit❗❌:     const Bond *bi = bondCache[path[i]];
    // RDKit❗❌: #ifdef REPORT_FP_STATS
    // RDKit❗❌:     if (std::find(atomsToUse.begin(), atomsToUse.end(),
    // RDKit❗❌:                   bi->getBeginAtomIdx()) == atomsToUse.end()) {
    // RDKit❗❌:       atomsToUse.push_back(bi->getBeginAtomIdx());
    // RDKit❗❌:     }
    // RDKit❗❌:     if (std::find(atomsToUse.begin(), atomsToUse.end(), bi->getEndAtomIdx()) ==
    // RDKit❗❌:         atomsToUse.end()) {
    // RDKit❗❌:       atomsToUse.push_back(bi->getEndAtomIdx());
    // RDKit❗❌:     }
    // RDKit❗❌: #endif
    // RDKit❗❌:     for (unsigned int j = i + 1; j < path.size(); ++j) {
    // RDKit❗❌:       const Bond *bj = bondCache[path[j]];
    // RDKit❗❌:       if (bi->getBeginAtomIdx() == bj->getBeginAtomIdx() ||
    // RDKit❗❌:           bi->getBeginAtomIdx() == bj->getEndAtomIdx() ||
    // RDKit❗❌:           bi->getEndAtomIdx() == bj->getBeginAtomIdx() ||
    // RDKit❗❌:           bi->getEndAtomIdx() == bj->getEndAtomIdx()) {
    // RDKit❗❌:         ++bondNbrs[i];
    // RDKit❗❌:         ++bondNbrs[j];
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗❌:     std::cerr << "   bond(" << i << "):" << bondNbrs[i] << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:     // we have the count of neighbors for bond bi, compute its hash:
    // RDKit❗❌:     unsigned int a1Hash = (*atomInvariants)[bi->getBeginAtomIdx()];
    // RDKit❗❌:     unsigned int a2Hash = (*atomInvariants)[bi->getEndAtomIdx()];
    // RDKit❗❌:     unsigned int deg1 = atomDegrees[bi->getBeginAtomIdx()];
    // RDKit❗❌:     unsigned int deg2 = atomDegrees[bi->getEndAtomIdx()];
    // RDKit❗❌:     if (a1Hash < a2Hash) {
    // RDKit❗❌:       std::swap(a1Hash, a2Hash);
    // RDKit❗❌:       std::swap(deg1, deg2);
    // RDKit❗❌:     } else if (a1Hash == a2Hash && deg1 < deg2) {
    // RDKit❗❌:       std::swap(deg1, deg2);
    // RDKit❗❌:     }
    // RDKit❗❌:     unsigned int bondHash = 1;
    // RDKit❗❌:     if (useBondOrder) {
    // RDKit❗❌:       if (bi->getIsAromatic() || bi->getBondType() == Bond::AROMATIC) {
    // RDKit❗❌:         // makes sure aromatic bonds always hash as aromatic
    // RDKit❗❌:         bondHash = Bond::AROMATIC;
    // RDKit❗❌:       } else {
    // RDKit❗❌:         bondHash = bi->getBondType();
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     std::uint32_t ourHash = bondNbrs[i];
    // RDKit❗❌:     gboost::hash_combine(ourHash, bondHash);
    // RDKit❗❌:     gboost::hash_combine(ourHash, a1Hash);
    // RDKit❗❌:     gboost::hash_combine(ourHash, deg1);
    // RDKit❗❌:     gboost::hash_combine(ourHash, a2Hash);
    // RDKit❗❌:     gboost::hash_combine(ourHash, deg2);
    // RDKit❗❌:     bondHashes.push_back(ourHash);
    // RDKit❗❌:     // std::cerr<<"    "<<bi->getIdx()<<"
    // RDKit❗❌:     // "<<a1Hash<<"("<<deg1<<")"<<"-"<<a2Hash<<"("<<deg2<<")"<<" "<<bondHash<<"
    // RDKit❗❌:     // -> "<<ourHash<<std::endl;
    // RDKit❗❌:   }
    // RDKit❗❌:   return bondHashes;
    // RDKit❗❌: }

    // Complexity: source O(A+P²) path degrees/neighbors retained. New byte mask
    // per path replaces reused packed source bits; this allocation/cost is worse.
    if atom_invariants.len() < atom_count(graph) {
        return Err(TopologicalFingerprintError::InvalidArguments {
            reason: "bad atomInvariants size",
        });
    }
    let mut atoms_in_path = vec![false; atom_count(graph)];
    let mut degrees = vec![0u32; atom_count(graph)];
    for &i in path {
        if i >= bonds(graph) {
            return Err(TopologicalFingerprintError::InvalidArguments {
                reason: "bond not in cache",
            });
        }
        let b = bond(graph, i);
        let a = b.begin().index();
        let z = b.end().index();
        degrees[a] = degrees[a].wrapping_add(1);
        degrees[z] = degrees[z].wrapping_add(1);
        atoms_in_path[a] = true;
        atoms_in_path[z] = true;
        if query_bonds[i] != 0 {
            return Ok(RdkitFpBondHashInputs {
                atoms_in_path,
                bond_hashes: Vec::new(),
            });
        }
    }
    let mut nbrs = vec![0u32; path.len()];
    for i in 0..path.len() {
        for j in i + 1..path.len() {
            let a = bond(graph, path[i]);
            let b = bond(graph, path[j]);
            if a.begin() == b.begin()
                || a.begin() == b.end()
                || a.end() == b.begin()
                || a.end() == b.end()
            {
                nbrs[i] = nbrs[i].wrapping_add(1);
                nbrs[j] = nbrs[j].wrapping_add(1);
            }
        }
    }
    let mut hashes = Vec::with_capacity(path.len() + 1);
    for (i, &id) in path.iter().enumerate() {
        let b = bond(graph, id);
        let a = b.begin().index();
        let z = b.end().index();
        let (mut ah, mut zh, mut ad, mut zd) = (
            atom_invariants[a],
            atom_invariants[z],
            degrees[a],
            degrees[z],
        );
        if ah < zh {
            std::mem::swap(&mut ah, &mut zh);
            std::mem::swap(&mut ad, &mut zd);
        } else if ah == zh && ad < zd {
            std::mem::swap(&mut ad, &mut zd);
        }
        let order = if use_bond_order {
            if b.is_aromatic() {
                cosmolkit_model::BondOrder::Aromatic.rdkit_code() as u32
            } else {
                b.order().rdkit_code() as u32
            }
        } else {
            1
        };
        let mut h = nbrs[i];
        for v in [order, ah, ad, zh, zd] {
            hash_combine(&mut h, v);
        }
        hashes.push(h);
    }
    Ok(RdkitFpBondHashInputs {
        atoms_in_path,
        bond_hashes: hashes,
    })
}
#[derive(Debug, Clone, PartialEq, Eq)]
struct RdkitFpEnvironment {
    bit_id: u32,
    atoms_in_path: Vec<bool>,
    bond_path: Vec<usize>,
}
impl RdkitFpEnvironment {
    fn bit_id(&self) -> u32 {
        self.bit_id
    }
    fn update_additional_output(&self, output: &mut FingerprintAdditionalOutput, bit_id: u64) {
        // RDKit❗✔️: void RDKitFPAtomEnv<OutputType>::updateAdditionalOutput(
        // RDKit❗✔️:     AdditionalOutput *additionalOutput, size_t bitId) const {
        // RDKit❗✔️:   PRECONDITION(additionalOutput, "bad output pointer");
        // RDKit❗✔️:   if (additionalOutput->bitPaths) {
        // RDKit❗✔️:     (*additionalOutput->bitPaths)[bitId].push_back(d_bondPath);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (additionalOutput->atomToBits || additionalOutput->atomCounts ||
        // RDKit❗✔️:       additionalOutput->atomsPerBit) {
        // RDKit❗✔️:     if (additionalOutput->atomsPerBit) {
        // RDKit❗✔️:       (*additionalOutput->atomsPerBit)[bitId].emplace_back();
        // RDKit❗✔️:     }
        // RDKit❗✔️:     for (size_t i = 0; i < d_atomsInPath.size(); ++i) {
        // RDKit❗✔️:       if (d_atomsInPath[i]) {
        // RDKit❗✔️:         if (additionalOutput->atomsPerBit) {
        // RDKit❗✔️:           (*additionalOutput->atomsPerBit)[bitId].back().push_back(i);
        // RDKit❗✔️:         }
        // RDKit❗✔️:         if (additionalOutput->atomToBits) {
        // RDKit❗✔️:           auto &alist = additionalOutput->atomToBits->at(i);
        // RDKit❗✔️:           if (std::find(alist.begin(), alist.end(), bitId) == alist.end()) {
        // RDKit❗✔️:             alist.push_back(bitId);
        // RDKit❗✔️:           }
        // RDKit❗✔️:         }
        // RDKit❗✔️:         if (additionalOutput->atomCounts) {
        // RDKit❗✔️:           additionalOutput->atomCounts->at(i)++;
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }

        // Complexity: source map insert O(logK)+path copy, atom mask scan O(A),
        // atom-to-bit linear dedup preserved; source unsigned counts wrap explicitly.
        if let Some(paths) = output.bit_paths.as_mut() {
            paths
                .entry(bit_id)
                .or_default()
                .push(self.bond_path.iter().map(|&i| i as i32).collect());
        }
        if output.atom_to_bits.is_some()
            || output.atom_counts.is_some()
            || output.atoms_per_bit.is_some()
        {
            if let Some(rows) = output.atoms_per_bit.as_mut() {
                rows.entry(bit_id).or_default().push(Vec::new());
            }
            for (atom, &yes) in self.atoms_in_path.iter().enumerate() {
                if !yes {
                    continue;
                }
                if let Some(rows) = output.atoms_per_bit.as_mut() {
                    rows.get_mut(&bit_id)
                        .expect("inserted source row")
                        .last_mut()
                        .expect("inserted occurrence")
                        .push(atom as i32);
                }
                if let Some(rows) = output.atom_to_bits.as_mut() {
                    if !rows[atom].contains(&bit_id) {
                        rows[atom].push(bit_id);
                    }
                }
                if let Some(rows) = output.atom_counts.as_mut() {
                    rows[atom] = rows[atom].wrapping_add(1);
                }
            }
        }
    }
}
impl FingerprintEnvironment<TopologicalFingerprintError> for RdkitFpEnvironment {
    type Output = u32;
    type State = ();
    fn bit_id(
        &self,
        _args: &FingerprintArguments,
        _atoms: &[u32],
        _bonds: &[u32],
        _output: Option<&mut FingerprintAdditionalOutput>,
        _hash: bool,
        _size: u64,
    ) -> Result<u32, TopologicalFingerprintError> {
        // RDKit❗✔️: OutputType RDKitFPAtomEnv<OutputType>::getBitId(
        // RDKit❗✔️:     FingerprintArguments *,              // arguments
        // RDKit❗✔️:     const std::vector<std::uint32_t> *,  // atomInvariants
        // RDKit❗✔️:     const std::vector<std::uint32_t> *,  // bondInvariants
        // RDKit❗✔️:     AdditionalOutput *,                  // additional Output
        // RDKit❗✔️:     const bool,                          // hashResults
        // RDKit❗✔️:     const std::uint64_t                  // fpSize
        // RDKit❗✔️: ) const {
        // RDKit❗✔️:   return d_bitId;
        // RDKit❗✔️: }

        Ok(self.bit_id)
    }
    fn update_output(
        &self,
        o: &mut FingerprintAdditionalOutput,
        b: u64,
        _s: &mut (),
    ) -> Result<(), TopologicalFingerprintError> {
        self.update_additional_output(o, b);
        Ok(())
    }
}
fn generate_rdkit_fp_environments(
    graph: GraphInput<'_>,
    params: &TopologicalFingerprintParams,
    atoms: &[u32],
) -> Result<Vec<RdkitFpEnvironment>, TopologicalFingerprintError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::RDKitFP::RDKitFPEnvGenerator<OutputType>::getEnvironments (Release_2026_03_6)
    // RDKit❗❌: std::vector<AtomEnvironment<OutputType> *>
    // RDKit❗❌: RDKitFPEnvGenerator<OutputType>::getEnvironments(
    // RDKit❗❌:     const ROMol &mol, FingerprintArguments *arguments,
    // RDKit❗❌:     const std::vector<std::uint32_t> *fromAtoms,
    // RDKit❗❌:     const std::vector<std::uint32_t> *ignoreAtoms,
    // RDKit❗❌:     const int,                 // confId
    // RDKit❗❌:     const AdditionalOutput *,  // additionalOutput
    // RDKit❗❌:     const std::vector<std::uint32_t> *atomInvariants,
    // RDKit❗❌:     const std::vector<std::uint32_t> *,  // bondInvariants
    // RDKit❗❌:     const bool                           // hashResults
    // RDKit❗❌: ) const {
    // RDKit❗❌:   PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit❗❌:                "bad atomInvariants size");
    // RDKit❗❌:
    // RDKit❗❌:   auto *fpArguments = dynamic_cast<RDKitFPArguments *>(arguments);
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<AtomEnvironment<OutputType> *> result;
    // RDKit❗❌:
    // RDKit❗❌:   // get all paths
    // RDKit❗❌:   INT_PATH_LIST_MAP allPaths;
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> ignoreAtomsBitset;
    // RDKit❗❌:   if (ignoreAtoms) {
    // RDKit❗❌:     ignoreAtomsBitset.resize(mol.getNumAtoms());
    // RDKit❗❌:     std::ranges::for_each(*ignoreAtoms, [&](const auto atomIdx) {
    // RDKit❗❌:       ignoreAtomsBitset.set(atomIdx);
    // RDKit❗❌:     });
    // RDKit❗❌:   }
    // RDKit❗❌:   RDKitFPUtils::enumerateAllPaths(
    // RDKit❗❌:       mol, allPaths, fromAtoms, fpArguments->df_branchedPaths,
    // RDKit❗❌:       fpArguments->df_useHs, fpArguments->d_minPath, fpArguments->d_maxPath,
    // RDKit❗❌:       ignoreAtoms ? &ignoreAtomsBitset : nullptr);
    // RDKit❗❌:
    // RDKit❗❌:   // identify query bonds
    // RDKit❗❌:   std::vector<short> isQueryBond(mol.getNumBonds(), 0);
    // RDKit❗❌:   std::vector<const Bond *> bondCache;
    // RDKit❗❌:   RDKitFPUtils::identifyQueryBonds(mol, bondCache, isQueryBond);
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> atomsInPath(mol.getNumAtoms());
    // RDKit❗❌:   for (INT_PATH_LIST_MAP_CI paths = allPaths.begin(); paths != allPaths.end();
    // RDKit❗❌:        paths++) {
    // RDKit❗❌:     for (const auto &path : paths->second) {
    // RDKit❗❌:       // the bond hashes of the path
    // RDKit❗❌:       std::vector<std::uint32_t> bondHashes = RDKitFPUtils::generateBondHashes(
    // RDKit❗❌:           mol, atomsInPath, bondCache, isQueryBond, path,
    // RDKit❗❌:           fpArguments->df_useBondOrder, atomInvariants);
    // RDKit❗❌:       if (!bondHashes.size()) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       // hash the path to generate a seed:
    // RDKit❗❌:       unsigned long seed;
    // RDKit❗❌:       if (path.size() > 1) {
    // RDKit❗❌:         std::sort(bondHashes.begin(), bondHashes.end());
    // RDKit❗❌:
    // RDKit❗❌:         // finally, we will add the number of distinct atoms in the path at the
    // RDKit❗❌:         // end
    // RDKit❗❌:         // of the vect. This allows us to distinguish C1CC1 from CC(C)C
    // RDKit❗❌:         bondHashes.push_back(static_cast<std::uint32_t>(atomsInPath.count()));
    // RDKit❗❌:         seed = gboost::hash_range(bondHashes.begin(), bondHashes.end());
    // RDKit❗❌:       } else {
    // RDKit❗❌:         seed = bondHashes[0];
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       result.push_back(new RDKitFPAtomEnv<OutputType>(
    // RDKit❗❌:           static_cast<OutputType>(seed), atomsInPath, path));
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return result;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RDKit::RDKitFP::RDKitFPEnvGenerator<OutputType>::getEnvironments

    // Complexity: CORE shared path enumeration, owned path map then environment
    // vector retains source traversal/order. Byte masks add memory vs packed bits.
    if atoms.len() < atom_count(graph) {
        return Err(TopologicalFingerprintError::InvalidArguments {
            reason: "bad atomInvariants size",
        });
    }
    let mut ignore_mask = params
        .ignore_atoms
        .as_ref()
        .map(|_| vec![false; atom_count(graph)]);
    if let (Some(indices), Some(mask)) = (params.ignore_atoms.as_ref(), ignore_mask.as_mut()) {
        for &index in indices {
            let slot = mask.get_mut(index as usize).ok_or(
                TopologicalFingerprintError::InvalidArguments {
                    reason: "ignoreAtoms atom out of range",
                },
            )?;
            *slot = true;
        }
    }
    let paths = enumerate_fingerprint_paths(
        graph,
        params.min_path,
        params.max_path,
        params.use_hs,
        params.branched_paths,
        params.from_atoms.as_deref(),
        ignore_mask.as_deref(),
    )?;
    let query_bonds = identify_query_bonds(graph);
    let mut result = Vec::new();
    for rows in paths.into_values() {
        for path in rows {
            let inputs = rdkit_fp_generate_bond_hash_inputs(
                graph,
                &path,
                params.use_bond_order,
                atoms,
                &query_bonds,
            )?;
            if inputs.bond_hashes.is_empty() {
                continue;
            }
            let seed = if path.len() > 1 {
                let mut h = inputs.bond_hashes;
                h.sort_unstable();
                h.push(inputs.atoms_in_path.iter().filter(|&&a| a).count() as u32);
                hash_range(&h)
            } else {
                inputs.bond_hashes[0]
            };
            result.push(RdkitFpEnvironment {
                bit_id: seed,
                atoms_in_path: inputs.atoms_in_path,
                bond_path: path,
            });
        }
    }
    Ok(result)
}
pub fn topological_fingerprint(
    topology: &TopologyBlock,
    params: &TopologicalFingerprintParams,
) -> Result<Fingerprint, TopologicalFingerprintError> {
    Ok(
        fingerprint_for_graph(GraphInput::Concrete(topology), params, Default::default())?
            .fingerprint,
    )
}
pub fn topological_fingerprint_with_output(
    topology: &TopologyBlock,
    params: &TopologicalFingerprintParams,
    request: TopologicalFingerprintOutputRequest,
) -> Result<TopologicalFingerprintResult, TopologicalFingerprintError> {
    fingerprint_for_graph(GraphInput::Concrete(topology), params, request)
}
pub fn topological_query_fingerprint(
    query: &QueryGraph,
    params: &TopologicalFingerprintParams,
) -> Result<Fingerprint, TopologicalFingerprintError> {
    Ok(fingerprint_for_graph(GraphInput::Query(query), params, Default::default())?.fingerprint)
}
pub fn topological_query_fingerprint_with_output(
    query: &QueryGraph,
    params: &TopologicalFingerprintParams,
    request: TopologicalFingerprintOutputRequest,
) -> Result<TopologicalFingerprintResult, TopologicalFingerprintError> {
    fingerprint_for_graph(GraphInput::Query(query), params, request)
}
fn fingerprint_for_graph(
    graph: GraphInput<'_>,
    params: &TopologicalFingerprintParams,
    request: TopologicalFingerprintOutputRequest,
) -> Result<TopologicalFingerprintResult, TopologicalFingerprintError> {
    // RDKit❗❌: ExplicitBitVect *RDKFingerprintMol(
    // RDKit❗❌:     const ROMol &mol, unsigned int minPath, unsigned int maxPath,
    // RDKit❗❌:     unsigned int fpSize, unsigned int nBitsPerHash, bool useHs,
    // RDKit❗❌:     double tgtDensity, unsigned int minSize, bool branchedPaths,
    // RDKit❗❌:     bool useBondOrder, std::vector<std::uint32_t> *atomInvariants,
    // RDKit❗❌:     const std::vector<std::uint32_t> *fromAtoms,
    // RDKit❗❌:     std::vector<std::vector<std::uint32_t>> *atomBits,
    // RDKit❗❌:     std::map<std::uint32_t, std::vector<std::vector<int>>> *bitInfo) {
    // RDKit❗❌:   PRECONDITION(minPath != 0, "minPath==0");
    // RDKit❗❌:   PRECONDITION(maxPath >= minPath, "maxPath<minPath");
    // RDKit❗❌:   PRECONDITION(fpSize != 0, "fpSize==0");
    // RDKit❗❌:   PRECONDITION(nBitsPerHash != 0, "nBitsPerHash==0");
    // RDKit❗❌:   PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit❗❌:                "bad atomInvariants size");
    // RDKit❗❌:   PRECONDITION(!atomBits || atomBits->size() >= mol.getNumAtoms(),
    // RDKit❗❌:                "bad atomBits size");
    // RDKit❗❌:
    // RDKit❗❌:   std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
    // RDKit❗❌:       RDKit::RDKitFP::getRDKitFPGenerator<std::uint32_t>(
    // RDKit❗❌:           minPath, maxPath, useHs, branchedPaths, useBondOrder));
    // RDKit❗❌:   fpgen->getOptions()->d_fpSize = fpSize;
    // RDKit❗❌:   fpgen->getOptions()->d_numBitsPerFeature = nBitsPerHash;
    // RDKit❗❌:
    // RDKit❗❌:   FingerprintFuncArguments args;
    // RDKit❗❌:   args.customAtomInvariants = atomInvariants;
    // RDKit❗❌:   args.fromAtoms = fromAtoms;
    // RDKit❗❌:
    // RDKit❗❌:   AdditionalOutput ao;
    // RDKit❗❌:   if (atomBits) {
    // RDKit❗❌:     args.additionalOutput = &ao;
    // RDKit❗❌:     ao.allocateAtomToBits();
    // RDKit❗❌:   }
    // RDKit❗❌:   if (bitInfo) {
    // RDKit❗❌:     args.additionalOutput = &ao;
    // RDKit❗❌:     ao.allocateBitPaths();
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   auto res = fpgen->getFingerprint(mol, args).release();
    // RDKit❗❌:
    // RDKit❗❌:   if (atomBits) {
    // RDKit❗❌:     atomBits->clear();
    // RDKit❗❌:     for (const auto &abl : *ao.atomToBits) {
    // RDKit❗❌:       std::vector<std::uint32_t> uv;
    // RDKit❗❌:       uv.reserve(abl.size());
    // RDKit❗❌:       for (auto l : abl) {
    // RDKit❗❌:         uv.push_back(static_cast<std::uint32_t>(l));
    // RDKit❗❌:       }
    // RDKit❗❌:       atomBits->emplace_back(std::move(uv));
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (bitInfo) {
    // RDKit❗❌:     bitInfo->clear();
    // RDKit❗❌:     for (const auto &abl : *ao.bitPaths) {
    // RDKit❗❌:       (*bitInfo)[static_cast<std::uint32_t>(abl.first)] = abl.second;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // EFF: this could be faster by folding by more than a factor
    // RDKit❗❌:   // of 2 each time, but we're not going to be spending much
    // RDKit❗❌:   // time here anyway
    // RDKit❗❌:   if (tgtDensity > 0.0) {
    // RDKit❗❌:     while (static_cast<double>(res->getNumOnBits()) / res->getNumBits() <
    // RDKit❗❌:                tgtDensity &&
    // RDKit❗❌:            res->getNumBits() >= 2 * minSize) {
    // RDKit❗❌:       ExplicitBitVect *tmpV = FoldFingerprint(*res, 2);
    // RDKit❗❌:       delete res;
    // RDKit❗❌:       res = tmpV;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }

    // Complexity: borrow custom invariants; shared source generator accumulation,
    // RNG, dense projection and folding owners reused without duplicate loops.
    // Graph validation and byte environment masks add O(A+B)/per-path memory.
    params.validate()?;
    if params
        .atom_invariants
        .as_ref()
        .is_some_and(|v| v.len() < atom_count(graph))
    {
        return Err(TopologicalFingerprintError::InvalidArguments {
            reason: "bad atomInvariants size",
        });
    }
    validate_graph(graph)?;
    let atoms = match params.atom_invariants.as_deref() {
        Some(a) => Cow::Borrowed(a),
        None => Cow::Owned(rdkit_fp_atom_invariants(graph)),
    };
    let mut ao = FingerprintAdditionalOutput::new();
    if request.atom_bits {
        ao.allocate_atom_to_bits();
    }
    if request.bit_info {
        ao.allocate_bit_paths();
    }
    let enabled = request.atom_bits || request.bit_info;
    if enabled {
        ao.reinitialize(atom_count(graph));
    }
    let common = FingerprintArguments {
        fp_size: params.fp_size,
        bits_per_feature: params.num_bits_per_feature,
        ..Default::default()
    };
    let envs = generate_rdkit_fp_environments(graph, params, &atoms)?;
    let mut fp = project_fingerprint(
        &common,
        atom_count(graph),
        enabled.then_some(&mut ao),
        |size, output| {
            let mut counts = SparseCountFingerprint::new(size);
            accumulate_sparse_counts_into(envs, &common, &atoms, &[], size, output, &mut counts)?;
            Ok::<_, TopologicalFingerprintError>(counts)
        },
    )?;
    // C++2*minSize is unsigned32 arithmetic; preserve wrap, delegate invalid
    // fold-factor source errors, and never substitute a truncated/empty result.
    while params.target_density > 0.0
        && f64::from(fp.num_on_bits()) / f64::from(fp.n_bits()) < params.target_density
        && fp.n_bits() >= params.min_size.wrapping_mul(2)
    {
        fp = crate::folding::fold_fingerprint(&fp, 2)?;
    }
    Ok(TopologicalFingerprintResult {
        fingerprint: fp,
        output: TopologicalFingerprintOutput {
            atom_bits: ao.atom_to_bits.map(|rows| {
                rows.into_iter()
                    .map(|r| r.into_iter().map(|b| b as u32).collect())
                    .collect()
            }),
            bit_info: ao
                .bit_paths
                .map(|m| m.into_iter().map(|(b, p)| (b as u32, p)).collect()),
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::TestMolecule;
    enum Molecule {
        Concrete(TestMolecule),
        Query(QueryGraph),
    }
    impl Molecule {
        fn new() -> Self {
            Self::Concrete(TestMolecule::new())
        }
        fn from_smiles(text: &str) -> Result<Self, cosmolkit_smiles::SmilesParseError> {
            TestMolecule::from_smiles(text).map(Self::Concrete)
        }
        fn graph(&self) -> GraphInput<'_> {
            match self {
                Self::Concrete(m) => GraphInput::Concrete(&m.topology),
                Self::Query(q) => GraphInput::Query(q),
            }
        }
        fn num_atoms(&self) -> usize {
            atom_count(self.graph())
        }
        fn with_hydrogens(self) -> Result<Self, cosmolkit_core::HydrogenError> {
            let Self::Concrete(m) = self else {
                panic!("ordinary source hydrogen fixture")
            };
            let r = cosmolkit_core::add_hydrogens_impl(m.topology, m.coordinates, m.properties)?;
            Ok(Self::Concrete(TestMolecule::prepared(
                r.topology,
                r.coordinates,
            )))
        }
    }
    fn query_fixture(s: &str) -> Result<Molecule, cosmolkit_search::SmartsParseError> {
        cosmolkit_search::parse_smarts(s, &Default::default()).map(Molecule::Query)
    }
    fn topological_fingerprint(
        m: &Molecule,
        p: &TopologicalFingerprintParams,
    ) -> Result<Fingerprint, TopologicalFingerprintError> {
        Ok(fingerprint_for_graph(m.graph(), p, Default::default())?.fingerprint)
    }
    fn topological_fingerprint_with_output(
        m: &Molecule,
        p: &TopologicalFingerprintParams,
        r: TopologicalFingerprintOutputRequest,
    ) -> Result<TopologicalFingerprintResult, TopologicalFingerprintError> {
        fingerprint_for_graph(m.graph(), p, r)
    }
    fn expected_paths(
        m: &Molecule,
        a: u32,
        z: u32,
        h: bool,
        b: bool,
        roots: Option<&[u32]>,
    ) -> BTreeMap<usize, Vec<Vec<usize>>> {
        enumerate_fingerprint_paths(m.graph(), a, z, h, b, roots, None)
            .expect("valid source path arguments")
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
            Err(TopologicalFingerprintError::InvalidArguments {
                reason: "minPath==0"
            })
        ));
        params = TopologicalFingerprintParams::default();
        params.max_path = 0;
        assert!(matches!(
            params.validate(),
            Err(TopologicalFingerprintError::InvalidArguments {
                reason: "maxPath<minPath"
            })
        ));
        params = TopologicalFingerprintParams::default();
        params.fp_size = 0;
        assert!(matches!(
            params.validate(),
            Err(TopologicalFingerprintError::InvalidArguments {
                reason: "fpSize==0"
            })
        ));
        params = TopologicalFingerprintParams::default();
        params.num_bits_per_feature = 0;
        assert!(matches!(
            params.validate(),
            Err(TopologicalFingerprintError::InvalidArguments {
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
        assert_eq!(
            rdkit_fp_atom_invariants(aliphatic.graph()),
            vec![12, 12, 16]
        );
        let aromatic = Molecule::from_smiles("c1ccccc1").expect("aromatic fixture");
        assert_eq!(rdkit_fp_atom_invariants(aromatic.graph()), vec![13; 6]);
    }

    #[test]
    fn rdkit_fp_bond_hash_inputs_match_boost_hash_combine_order() {
        let molecule = Molecule::from_smiles("CCO").expect("fixture");
        let invariants = rdkit_fp_atom_invariants(molecule.graph());
        let first = rdkit_fp_generate_bond_hash_inputs(
            molecule.graph(),
            &[0],
            true,
            &invariants,
            &identify_query_bonds(molecule.graph()),
        )
        .expect("first bond hash");
        assert_eq!(first.atoms_in_path, vec![true, true, false]);
        assert_eq!(first.bond_hashes, vec![4_275_705_116]);
        let second = rdkit_fp_generate_bond_hash_inputs(
            molecule.graph(),
            &[1],
            true,
            &invariants,
            &identify_query_bonds(molecule.graph()),
        )
        .expect("second bond hash");
        assert_eq!(second.bond_hashes, vec![4_274_652_475]);
        let without_order = rdkit_fp_generate_bond_hash_inputs(
            molecule.graph(),
            &[0],
            false,
            &invariants,
            &identify_query_bonds(molecule.graph()),
        )
        .expect("order-free bond hash");
        assert_eq!(first.bond_hashes, without_order.bond_hashes);
        let double_bond = Molecule::from_smiles("C=C").expect("double-bond fixture");
        let double_invariants = rdkit_fp_atom_invariants(double_bond.graph());
        let with_order = rdkit_fp_generate_bond_hash_inputs(
            double_bond.graph(),
            &[0],
            true,
            &double_invariants,
            &identify_query_bonds(double_bond.graph()),
        )
        .expect("double bond hash");
        let without_order = rdkit_fp_generate_bond_hash_inputs(
            double_bond.graph(),
            &[0],
            false,
            &double_invariants,
            &identify_query_bonds(double_bond.graph()),
        )
        .expect("double bond order-free hash");
        assert_ne!(with_order.bond_hashes, without_order.bond_hashes);
    }

    #[test]
    fn rdkit_fp_query_bond_path_is_rejected_at_hash_input_boundary() {
        let query = query_fixture("[#6]~[#6]").expect("query fixture");
        let invariants = rdkit_fp_atom_invariants(query.graph());
        let inputs = rdkit_fp_generate_bond_hash_inputs(
            query.graph(),
            &[0],
            true,
            &invariants,
            &identify_query_bonds(query.graph()),
        )
        .expect("query path result");
        assert!(inputs.bond_hashes.is_empty());
        assert_eq!(inputs.atoms_in_path, vec![true, true]);
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
        let invariants = rdkit_fp_atom_invariants(molecule.graph());
        let environments = generate_rdkit_fp_environments(molecule.graph(), &params, &invariants)
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

        let mut output = FingerprintAdditionalOutput::new();
        output.allocate_atom_to_bits();
        output.allocate_bit_paths();
        output.allocate_atom_counts();
        output.allocate_atoms_per_bit();
        output.reinitialize(molecule.num_atoms());
        for environment in &environments {
            environment.update_additional_output(&mut output, u64::from(environment.bit_id()));
        }
        environments[0].update_additional_output(&mut output, u64::from(environments[0].bit_id()));
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
    fn search01_fingerprint_ignored_ids_cover_concrete_query_empty_and_duplicates() {
        for m in [
            Molecule::from_smiles("CC(C)C").unwrap(),
            query_fixture("CC(C)C").unwrap(),
        ] {
            for branched in [false, true] {
                let mut p = TopologicalFingerprintParams {
                    branched_paths: branched,
                    ..Default::default()
                };
                let baseline = topological_fingerprint(&m, &p).unwrap();
                assert!(!baseline.on_bits().is_empty());
                p.ignore_atoms = Some(vec![]);
                assert_eq!(topological_fingerprint(&m, &p).unwrap(), baseline);
                p.ignore_atoms = Some(vec![1, 1]);
                assert!(
                    topological_fingerprint(&m, &p)
                        .unwrap()
                        .on_bits()
                        .is_empty()
                );
                p.from_atoms = Some(vec![1, 1]);
                assert!(
                    topological_fingerprint(&m, &p)
                        .unwrap()
                        .on_bits()
                        .is_empty()
                );
            }
        }
    }
    #[test]
    fn search01_fingerprint_root_prepend_order_and_ignored_endpoint_mask() {
        for m in [
            Molecule::from_smiles("CC(C)C").unwrap(),
            query_fixture("CC(C)C").unwrap(),
        ] {
            let mask = [true, false, false, false];
            for branched in [false, true] {
                let paths = enumerate_fingerprint_paths(
                    m.graph(),
                    1,
                    1,
                    true,
                    branched,
                    Some(&[2, 3, 2]),
                    Some(&mask),
                )
                .unwrap();
                assert_eq!(paths[&1], vec![vec![1], vec![2], vec![1]]);
                let no_roots =
                    enumerate_fingerprint_paths(m.graph(), 1, 1, true, branched, None, Some(&mask))
                        .unwrap();
                assert_eq!(no_roots[&1], vec![vec![1], vec![2]]);
                assert!(
                    enumerate_fingerprint_paths(
                        m.graph(),
                        1,
                        1,
                        true,
                        branched,
                        Some(&[0]),
                        Some(&mask)
                    )
                    .unwrap()
                    .values()
                    .all(Vec::is_empty)
                );
            }
        }
    }
    #[test]
    fn search01_generator_checks_atom_invariants_before_ignored_index_and_empty_roots() {
        let m = Molecule::from_smiles("CC").unwrap();
        let p = TopologicalFingerprintParams {
            ignore_atoms: Some(vec![2]),
            ..Default::default()
        };
        assert!(matches!(
            generate_rdkit_fp_environments(m.graph(), &p, &[]),
            Err(TopologicalFingerprintError::InvalidArguments {
                reason: "bad atomInvariants size"
            })
        ));
        assert!(matches!(
            generate_rdkit_fp_environments(m.graph(), &p, &[12, 12]),
            Err(TopologicalFingerprintError::InvalidArguments {
                reason: "ignoreAtoms atom out of range"
            })
        ));
        let p = TopologicalFingerprintParams {
            from_atoms: Some(vec![]),
            ignore_atoms: Some(vec![0]),
            ..Default::default()
        };
        assert!(
            generate_rdkit_fp_environments(m.graph(), &p, &[12, 12])
                .unwrap()
                .is_empty()
        );
    }
}
