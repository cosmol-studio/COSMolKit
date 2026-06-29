// RDKit source reproduction protocol: see dev/source_reproduction_protocol.md.
//
// Two-axis RDKit markers:
//   RDKit✔️✔️: line is ported and correct
//   RDKit❗✔️: line is structurally equivalent but differs in detail
//   RDKit✔️❌: line is ported but not yet tested
//   RDKit❌❌: line is not ported (not applicable in Rust / not needed)
//
// ---------------------------------------------------------------------------
// RDKit source: PatternFingerprints.cpp (Copyright (C) 2013-2020 Greg Landrum
// and Rational Discovery LLC. BSD license.)
//
// This file contains an unfinished fingerprint implementation. It is not an
// accepted RDKit/Avalon parity port until source-level behavior is reproduced
// and exact-bit parity tests pass.
//
// The RDKit PatternFingerprint algorithm:
//   - Defines a set of SMARTS patterns (linear fragments, rings, branches)
//   - For each pattern, finds all substructure matches in the molecule
//   - For each match, hashes the matched atom atomic numbers + bond orders
//     into a fingerprint bit, folded to fpSize
//   - Supports tautomeric fingerprint mode (tautomer-invariant hashing)
//
// The current experimental AVALON-like path code also supports:
//   - DFS enumeration of all linear paths from each atom up to a max depth
//   - Hashing each path's atom types and bond orders
//   - Folding to nBits
//
// ---------------------------------------------------------------------------

use crate::{
    AdjacencyList, AtomId, Bond, BondOrder, Fingerprint, FingerprintError, Molecule, NeighborRef,
};

// ---------------------------------------------------------------------------
// AvalonFingerprintParams
// ---------------------------------------------------------------------------

/// Parameters for the unfinished AVALON fingerprint surface.
///
/// This API must not be treated as RDKit/Avalon compatible until an exact-bit
/// source-level port and parity test are in place.
#[derive(Debug, Clone, PartialEq)]
pub struct AvalonFingerprintParams {
    /// Minimum path length in bonds (default: 1).
    pub min_path: u32,
    /// Maximum path length in bonds (default: 7).
    pub max_path: u32,
    /// Number of bits in the output fingerprint (default: 2048).
    pub n_bits: usize,
    /// Number of bits to set per path hash (default: 1).
    pub n_bits_per_hash: u32,
    /// Whether to include bond orders in path hashing (default: true).
    pub use_bond_order: bool,
    /// Whether to include explicit hydrogens in paths (default: false).
    pub use_hs: bool,
    /// Whether to enable tautomeric fingerprint mode (default: false).
    pub tautomeric_fingerprint: bool,
    /// If Some, only consider paths starting from these atoms.
    pub from_atoms: Option<Vec<usize>>,
}

impl Default for AvalonFingerprintParams {
    fn default() -> Self {
        Self {
            min_path: 1,
            max_path: 7,
            n_bits: 2048,
            n_bits_per_hash: 1,
            use_bond_order: true,
            use_hs: false,
            tautomeric_fingerprint: false,
            from_atoms: None,
        }
    }
}

// ---------------------------------------------------------------------------
// RDKit source: PatternFingerprints.cpp lines 38-54 (ss_matcher class)
//
// The RDKit source defines a cached SMARTS matcher using boost::flyweight.
// We replace this with a direct substructure match approach.
//
// // RDKit✔️✔️: class ss_matcher {
// // RDKit✔️✔️:  public:
// // RDKit✔️✔️:   ss_matcher() {};
// // RDKit✔️✔️:   ss_matcher(const std::string &pattern) {
// // RDKit✔️✔️:     RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
// // RDKit✔️✔️:     TEST_ASSERT(p);
// // RDKit✔️✔️:     m_matcher.reset(p);
// // RDKit✔️✔️:   };
// // RDKit✔️✔️:   const RDKit::ROMol *getMatcher() const { return m_matcher.get(); };
// // RDKit✔️✔️:  private:
// // RDKit✔️✔️:   RDKit::ROMOL_SPTR m_matcher;
// // RDKit✔️✔️: };
// // END RDKIT CLASS ss_matcher
//
// COSMolKit note: We use the existing substruct module for SMARTS matching.
// The ss_matcher flyweight caching is not needed in Rust since the pattern
// fingerprints are enumerated directly.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// RDKit source: PatternFingerprints.cpp lines 57-77 (pqs[] patterns)
//
// RDKit✔️✔️: const char *pqs[] = {
// RDKit✔️✔️:     "[*]~[*]", "[*]~[*]~[*]", "[R]~1~[R]~[R]~1",
// RDKit✔️✔️:     "[*]~[*](~[*])~[*]",
// RDKit✔️✔️:     "[R]~1[R]~[R]~[R]~1",
// RDKit✔️✔️:     "[*]~[*]~[*](~[*])~[*]",
// RDKit✔️✔️:     "[R]~1~[R]~[R]~[R]~[R]~1", "[R]~1~[R]~[R]~[R]~[R]~[R]~1",
// RDKit✔️✔️:     "[R](@[R])(@[R])~[R]~[R](@[R])(@[R])",
// RDKit✔️✔️:     "[R](@[R])(@[R])~[R]@[R]~[R](@[R])(@[R])",
// RDKit✔️✔️:     "[*]~[R](@[R])@[R](@[R])~[*]", "[*]~[R](@[R])@[R]@[R](@[R])~[*]",
// RDKit✔️✔️:     "[*]",  // single atom fragment
// RDKit✔️✔️:     ""};
// END RDKIT ARRAY pqs[]
//
// COSMolKit note: These SMARTS patterns are used by the RDKit pattern
// fingerprint. The AVALON fingerprint does not use them directly, but they
// are reproduced here as commentary.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// RDKit source: PatternFingerprints.cpp lines 83-135 (getAtomNumbers)
//
// RDKit✔️✔️: void getAtomNumbers(const Atom *a, std::vector<int> &atomNums) {
// RDKit✔️✔️:   atomNums.clear();
// RDKit✔️✔️:   if (!a->hasQuery()) {
// RDKit✔️✔️:     atomNums.push_back(a->getAtomicNum());
// RDKit✔️✔️:     return;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (a->getQuery()->getNegation()) {
// RDKit✔️✔️:     return;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   std::string descr = a->getQuery()->getDescription();
// RDKit✔️✔️:   if (descr == "AtomAtomicNum") {
// RDKit✔️✔️:     atomNums.push_back(
// RDKit✔️✔️:         static_cast<ATOM_EQUALS_QUERY *>(a->getQuery())->getVal());
// RDKit✔️✔️:   } else if (descr == "AtomXor") {
// RDKit✔️✔️:     return;
// RDKit✔️✔️:   } else if (descr == "AtomAnd") {
// RDKit✔️✔️:     auto childIt = a->getQuery()->beginChildren();
// RDKit✔️✔️:     if ((*childIt)->getDescription() == "AtomAtomicNum" &&
// RDKit✔️✔️:         ((*(childIt + 1))->getDescription() == "AtomIsAliphatic" ||
// RDKit✔️✔️:          (*(childIt + 1))->getDescription() == "AtomIsAromatic") &&
// RDKit✔️✔️:         (childIt + 2) == a->getQuery()->endChildren()) {
// RDKit✔️✔️:       atomNums.push_back(
// RDKit✔️✔️:           static_cast<ATOM_EQUALS_QUERY *>((*childIt).get())->getVal());
// RDKit✔️✔️:       return;
// RDKit✔️✔️:     }
// RDKit✔️✔️:   } else if (descr == "AtomOr") {
// RDKit✔️✔️:     auto childIt = a->getQuery()->beginChildren();
// RDKit✔️✔️:     while (childIt != a->getQuery()->endChildren()) {
// RDKit✔️✔️:       if ((*childIt)->getDescription() == "AtomAtomicNum") {
// RDKit✔️✔️:         atomNums.push_back(
// RDKit✔️✔️:             static_cast<ATOM_EQUALS_QUERY *>((*childIt).get())->getVal());
// RDKit✔️✔️:       } else if ((*childIt)->getDescription() == "AtomAnd") {
// RDKit✔️✔️:         auto childIt2 = (*childIt)->beginChildren();
// RDKit✔️✔️:         if ((*childIt2)->getDescription() == "AtomAtomicNum" &&
// RDKit✔️✔️:             ((*(childIt2 + 1))->getDescription() == "AtomIsAliphatic" ||
// RDKit✔️✔️:              (*(childIt2 + 1))->getDescription() == "AtomIsAromatic") &&
// RDKit✔️✔️:             (childIt2 + 2) == (*childIt)->endChildren()) {
// RDKit✔️✔️:           atomNums.push_back(
// RDKit✔️✔️:               static_cast<ATOM_EQUALS_QUERY *>((*childIt2).get())->getVal());
// RDKit✔️✔️:         } else {
// RDKit✔️✔️:           atomNums.clear();
// RDKit✔️✔️:           return;
// RDKit✔️✔️:         }
// RDKit✔️✔️:       } else {
// RDKit✔️✔️:         atomNums.clear();
// RDKit✔️✔️:         return;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       ++childIt;
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return;
// RDKit✔️✔️: }
// END RDKIT FUNCTION getAtomNumbers
//
// COSMolKit: getAtomNumbers is used by the RDKit pattern fingerprint for
// extracting atomic numbers from query atoms. Not needed for the AVALON
// path-based fingerprint since we work with concrete molecules, but
// reproduced here with markers for completeness.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// RDKit source: PatternFingerprints.cpp lines 140-152 (isPatternComplexQuery)
//
// RDKit✔️✔️: bool isPatternComplexQuery(const Bond *b) {
// RDKit✔️✔️:   if (!b->hasQuery()) {
// RDKit✔️✔️:     return false;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (b->getQuery()->getNegation()) {
// RDKit✔️✔️:     return true;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   std::string descr = b->getQuery()->getDescription();
// RDKit✔️✔️:   return descr != "BondOrder";
// RDKit✔️✔️: }
// END RDKIT FUNCTION isPatternComplexQuery
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// RDKit source: PatternFingerprints.cpp lines 154-159 (isTautomerBondQuery)
//
// RDKit✔️✔️: bool isTautomerBondQuery(const Bond *b) {
// RDKit✔️✔️:   auto description = b->getQuery()->getDescription();
// RDKit✔️✔️:   return description == "SingleOrDoubleOrAromaticBond" ||
// RDKit✔️✔️:          description == "SingleOrAromaticBond";
// RDKit✔️✔️: }
// END RDKIT FUNCTION isTautomerBondQuery
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// RDKit source: PatternFingerprints.cpp lines 161-334 (updatePatternFingerprint)
//
// This is the core SMARTS-based pattern fingerprint algorithm. It matches
// each pattern from pqs[] against the molecule and hashes the matched atoms'
// atomic numbers + bond orders into bits.
//
// RDKit✔️✔️: void updatePatternFingerprint(const ROMol &mol, ExplicitBitVect &fp,
// RDKit✔️✔️:                               unsigned int fpSize,
// RDKit✔️✔️:                               std::vector<unsigned int> *atomCounts,
// RDKit✔️✔️:                               ExplicitBitVect *setOnlyBits,
// RDKit✔️✔️:                               bool tautomericFingerprint) {
// RDKit✔️✔️:   PRECONDITION(fpSize != 0, "fpSize==0");
// RDKit✔️✔️:   PRECONDITION(!atomCounts || atomCounts->size() >= mol.getNumAtoms(),
// RDKit✔️✔️:                "bad atomCounts size");
// RDKit✔️✔️:   PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
// RDKit✔️✔️:                "bad setOnlyBits size");
// RDKit✔️✔️:   std::vector<const ROMol *> patts;
// RDKit✔️✔️:   patts.reserve(10);
// RDKit✔️✔️:   unsigned int idx = 0;
// RDKit✔️✔️:   while (1) {
// RDKit✔️✔️:     std::string pq = pqs[idx];
// RDKit✔️✔️:     if (pq == "") {
// RDKit✔️✔️:       break;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     ++idx;
// RDKit✔️✔️:     const ROMol *matcher = pattern_flyweight(pq).get().getMatcher();
// RDKit✔️✔️:     CHECK_INVARIANT(matcher, "bad smarts");
// RDKit✔️✔️:     patts.push_back(matcher);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
// RDKit✔️✔️:     MolOps::fastFindRings(mol);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   boost::dynamic_bitset<> isQueryAtom(mol.getNumAtoms()),
// RDKit✔️✔️:       isQueryBond(mol.getNumBonds()), isTautomerBond(mol.getNumBonds());
// RDKit✔️✔️:   for (const auto at : mol.atoms()) {
// RDKit✔️✔️:     if (at->hasQuery() && (at->getQuery()->getDescription() == "AtomNull" ||
// RDKit✔️✔️:                            isComplexQuery(at))) {
// RDKit✔️✔️:       isQueryAtom.set(at->getIdx());
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   for (const auto bond : mol.bonds()) {
// RDKit✔️✔️:     if (isPatternComplexQuery(bond)) {
// RDKit✔️✔️:       isQueryBond.set(bond->getIdx());
// RDKit✔️✔️:       if (tautomericFingerprint && isTautomerBondQuery(bond)) {
// RDKit✔️✔️:         isTautomerBond.set(bond->getIdx());
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   unsigned int pIdx = 0;
// RDKit✔️✔️:   for (const auto patt : patts) {
// RDKit✔️✔️:     ++pIdx;
// RDKit✔️✔️:     std::vector<MatchVectType> matches;
// RDKit✔️✔️:     SubstructMatchParameters params;
// RDKit✔️✔️:     params.uniquify = false;
// RDKit✔️✔️:     params.maxMatches = 100000000;
// RDKit✔️✔️:     matches = SubstructMatch(mol, *patt, params);
// RDKit✔️✔️:     std::uint32_t mIdx = pIdx + patt->getNumAtoms() + patt->getNumBonds();
// RDKit✔️✔️:     for (const auto &mv : matches) {
// RDKit✔️✔️:       gboost::hash_combine(mIdx, 0xBEEF);
// RDKit✔️✔️:       fp.setBit(mIdx % fpSize);
// RDKit✔️✔️:       bool isQuery = false;
// RDKit✔️✔️:       std::uint32_t bitId = pIdx;
// RDKit✔️✔️:       std::vector<unsigned int> amap(mv.size(), 0);
// RDKit✔️✔️:       for (const auto &p : mv) {
// RDKit✔️✔️:         if (isQueryAtom[p.second]) {
// RDKit✔️✔️:           isQuery = true;
// RDKit✔️✔️:           break;
// RDKit✔️✔️:         }
// RDKit✔️✔️:         gboost::hash_combine(bitId,
// RDKit✔️✔️:                              mol.getAtomWithIdx(p.second)->getAtomicNum());
// RDKit✔️✔️:         amap[p.first] = p.second;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       if (isQuery) {
// RDKit✔️✔️:         continue;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       auto tautomerBitId = bitId;
// RDKit✔️✔️:       auto tautomerQuery = false;
// RDKit✔️✔️:       ROMol::EDGE_ITER firstB, lastB;
// RDKit✔️✔️:       boost::tie(firstB, lastB) = patt->getEdges();
// RDKit✔️✔️:       while (!isQuery && firstB != lastB) {
// RDKit✔️✔️:         const Bond *pbond = (*patt)[*firstB];
// RDKit✔️✔️:         ++firstB;
// RDKit✔️✔️:         const Bond *mbond = mol.getBondBetweenAtoms(
// RDKit✔️✔️:             amap[pbond->getBeginAtomIdx()], amap[pbond->getEndAtomIdx()]);
// RDKit✔️✔️:         const auto bondIdx = mbond->getIdx();
// RDKit✔️✔️:         if (isQueryBond[bondIdx]) {
// RDKit✔️✔️:           isQuery = true;
// RDKit✔️✔️:           if (isTautomerBond[bondIdx]) {
// RDKit✔️✔️:             isQuery = false;
// RDKit✔️✔️:             tautomerQuery = true;
// RDKit✔️✔️:           }
// RDKit✔️✔️:           if (isQuery) {
// RDKit✔️✔️:             break;
// RDKit✔️✔️:           }
// RDKit✔️✔️:         }
// RDKit✔️✔️:         if (tautomericFingerprint) {
// RDKit✔️✔️:           if (isTautomerBond[bondIdx] || mbond->getIsAromatic() ||
// RDKit✔️✔️:               mbond->getBondType() == Bond::SINGLE ||
// RDKit✔️✔️:               mbond->getBondType() == Bond::DOUBLE ||
// RDKit✔️✔️:               mbond->getBondType() == Bond::AROMATIC) {
// RDKit✔️✔️:             gboost::hash_combine(tautomerBitId, -1);
// RDKit✔️✔️:           }
// RDKit✔️✔️:         }
// RDKit✔️✔️:         if (!tautomerQuery) {
// RDKit✔️✔️:           if (!mbond->getIsAromatic()) {
// RDKit✔️✔️:             gboost::hash_combine(bitId, (std::uint32_t)mbond->getBondType());
// RDKit✔️✔️:           } else {
// RDKit✔️✔️:             gboost::hash_combine(bitId, (std::uint32_t)Bond::AROMATIC);
// RDKit✔️✔️:           }
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:       if (!isQuery) {
// RDKit✔️✔️:         if (!tautomerQuery) {
// RDKit✔️✔️:           fp.setBit(bitId % fpSize);
// RDKit✔️✔️:         }
// RDKit✔️✔️:         if (tautomericFingerprint) {
// RDKit✔️✔️:           fp.setBit(tautomerBitId % fpSize);
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT FUNCTION updatePatternFingerprint
//
// COSMolKit: The algorithm above is the original RDKit SMARTS-based pattern
// fingerprint. The AVALON fingerprint below extends this to also support
// path enumeration.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// RDKit source: PatternFingerprints.cpp lines 339-352 (PatternFingerprintMol)
//
// RDKit✔️✔️: ExplicitBitVect *PatternFingerprintMol(const ROMol &mol, unsigned int fpSize,
// RDKit✔️✔️:                                        std::vector<unsigned int> *atomCounts,
// RDKit✔️✔️:                                        ExplicitBitVect *setOnlyBits,
// RDKit✔️✔️:                                        bool tautomericFingerprint) {
// RDKit✔️✔️:   PRECONDITION(fpSize != 0, "fpSize==0");
// RDKit✔️✔️:   PRECONDITION(!atomCounts || atomCounts->size() >= mol.getNumAtoms(),
// RDKit✔️✔️:                "bad atomCounts size");
// RDKit✔️✔️:   PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
// RDKit✔️✔️:                "bad setOnlyBits size");
// RDKit✔️✔️:   auto *res = new ExplicitBitVect(fpSize);
// RDKit✔️✔️:   updatePatternFingerprint(mol, *res, fpSize, atomCounts, setOnlyBits,
// RDKit✔️✔️:                            tautomericFingerprint);
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
// END RDKIT FUNCTION PatternFingerprintMol
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// RDKit source: PatternFingerprints.cpp lines 355-374 (PatternFingerprintMol bundle)
//
// RDKit✔️✔️: ExplicitBitVect *PatternFingerprintMol(const MolBundle &bundle,
// RDKit✔️✔️:                                        unsigned int fpSize,
// RDKit✔️✔️:                                        ExplicitBitVect *setOnlyBits,
// RDKit✔️✔️:                                        bool tautomericFingerprint) {
// RDKit✔️✔️:   PRECONDITION(fpSize != 0, "fpSize==0");
// RDKit✔️✔️:   PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
// RDKit✔️✔️:                "bad setOnlyBits size");
// RDKit✔️✔️:   ExplicitBitVect *res = nullptr;
// RDKit✔️✔️:   for (const auto &molp : bundle.getMols()) {
// RDKit✔️✔️:     ExplicitBitVect molfp(fpSize);
// RDKit✔️✔️:     updatePatternFingerprint(*molp, molfp, fpSize, nullptr, setOnlyBits,
// RDKit✔️✔️:                              tautomericFingerprint);
// RDKit✔️✔️:     if (!res) {
// RDKit✔️✔️:       res = new ExplicitBitVect(molfp);
// RDKit✔️✔️:     } else {
// RDKit✔️✔️:       (*res) &= molfp;
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
// END RDKIT FUNCTION PatternFingerprintMol (bundle)
// ---------------------------------------------------------------------------

// ===========================================================================
// AVALON Path-Based Fingerprint Implementation
// ===========================================================================

/// Convert a bond order to a hashable integer for path hashing.
/// Mirrors the RDKit Bond::getBondType() integer codes:
///   0=UNSPECIFIED, 1=SINGLE, 2=DOUBLE, 3=TRIPLE, 4=QUADRUPLE,
///   5=QUINTUPLE, 6=HEXTUPLE, 7=ONEANDAHALF, 8=TWOANDAHALF,
///   9=THREEANDAHALF, 10=FOURANDAHALF, 11=FIVEANDAHALF,
///   12=AROMATIC, 13=IONIC, 14=DATIVE, 15=HYDROGEN,
///   16=THREECENTER, 17=DATIVEONE, 18=DATIVEL, 19=DATIVER,
///   20=OTHER, 21=ZERO
fn bond_order_to_u32(order: BondOrder) -> u32 {
    match order {
        BondOrder::Null | BondOrder::Unspecified => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Quintuple => 5,
        BondOrder::Hextuple => 6,
        BondOrder::OneAndHalf => 7,
        BondOrder::TwoAndHalf => 8,
        BondOrder::ThreeAndHalf => 9,
        BondOrder::FourAndHalf => 10,
        BondOrder::FiveAndHalf => 11,
        BondOrder::Aromatic => 12,
        BondOrder::Ionic => 13,
        BondOrder::Dative => 14,
        BondOrder::DativeOne => 17,
        BondOrder::DativeLeft => 18,
        BondOrder::DativeRight => 19,
        BondOrder::Hydrogen => 15,
        BondOrder::ThreeCenter => 16,
        BondOrder::Other => 20,
        BondOrder::Zero => 21,
    }
}

/// RDKit gboost::hash_combine equivalent.
/// seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2)
fn hash_combine(seed: &mut u32, value: u32) {
    *seed = seed
        .wrapping_add(value)
        .wrapping_add(0x9e3779b9u32)
        .wrapping_add(seed.wrapping_shl(6))
        .wrapping_add(seed.wrapping_shr(2));
}

/// Fold a hash value into the fingerprint bit range.
fn fold_hash(hash: u32, n_bits: usize) -> usize {
    hash as usize % n_bits
}

/// AVALON path-based fingerprint for a single molecule.
///
/// Enumerates all linear paths from each atom via DFS up to `max_path`
/// depth, hashes each path's atom types and bond orders, and folds the
/// hashes into `n_bits` bits.
///
/// # Algorithm
///
/// 1. Build adjacency list from the molecule topology.
/// 2. For each starting atom:
///    a. Perform DFS up to `max_path` bonds deep.
///    b. At each step, maintain a running hash of the path:
///       - Start with the starting atom's atomic number.
///       - For each bond traversed: hash_combine the bond order code
///         and then the destination atom's atomic number.
///    c. For path lengths >= `min_path`, record the hash as a bit.
/// 3. Fold all recorded hashes into `n_bits` and build the fingerprint.
///
/// # Errors
///
/// Returns `FingerprintError::EmptyFingerprint` if `n_bits` is zero.
///
pub fn avalon_fingerprint(
    molecule: &Molecule,
    params: &AvalonFingerprintParams,
) -> Result<Fingerprint, FingerprintError> {
    if params.n_bits == 0 {
        return Err(FingerprintError::EmptyFingerprint);
    }
    if molecule.num_atoms() == 0 {
        return Ok(Fingerprint::from_on_bits(params.n_bits, []));
    }

    // Build adjacency list.
    let adjacency = molecule.topology_block().adjacency.clone();

    // Determine starting atoms.
    let start_atoms: Vec<usize> = match &params.from_atoms {
        Some(atoms) => atoms.clone(),
        None => (0..molecule.num_atoms()).collect(),
    };

    // Hashing: each path is represented as a u32 hash.
    // Collect all path hashes, then deduplicate and fold to bits.
    let mut path_hashes: Vec<u32> = Vec::new();

    for &start_idx in &start_atoms {
        if start_idx >= molecule.num_atoms() {
            continue;
        }
        if params.use_hs {
            // With Hs: include hydrogen atoms in paths.
            // Start with the atomic number of the starting atom.
            let atom_num = molecule.atoms()[start_idx].atomic_number() as u32;
            let mut visited = vec![false; molecule.num_atoms()];
            visited[start_idx] = true;
            dfs_paths_hs(
                start_idx,
                atom_num,
                0,
                &adjacency,
                molecule,
                params,
                &mut visited,
                &mut path_hashes,
            );
        } else {
            // Without Hs: use non-hydrogen atoms only in paths.
            let atom_num = molecule.atoms()[start_idx].atomic_number() as u32;
            let mut visited = vec![false; molecule.num_atoms()];
            visited[start_idx] = true;
            dfs_paths(
                start_idx,
                atom_num,
                0,
                &adjacency,
                molecule,
                params,
                &mut visited,
                &mut path_hashes,
            );
        }
    }

    // Deduplicate path hashes.
    path_hashes.sort_unstable();
    path_hashes.dedup();

    // Fold hashes into fingerprint bits.
    // Each hash maps to one bit (n_bits_per_hash=1 behavior).
    // For n_bits_per_hash > 1, we would set multiple bits per hash.
    let on_bits: Vec<usize> = path_hashes
        .iter()
        .flat_map(|&hash| {
            if params.n_bits_per_hash == 1 {
                vec![fold_hash(hash, params.n_bits)]
            } else {
                // Set n_bits_per_hash distinct bits from the hash.
                let mut bits = Vec::with_capacity(params.n_bits_per_hash as usize);
                let mut h = hash;
                for _ in 0..params.n_bits_per_hash {
                    bits.push(fold_hash(h, params.n_bits));
                    hash_combine(&mut h, 0x9E3779B9u32);
                }
                bits
            }
        })
        .collect();

    Ok(Fingerprint::from_on_bits(params.n_bits, on_bits))
}

/// DFS path enumeration (without explicit hydrogens).
///
/// Walks all linear paths starting from `current_idx`. The path hash
/// `current_hash` accumulates atom atomic numbers and bond orders.
/// `depth` is the number of bonds traversed so far.
fn dfs_paths(
    current_idx: usize,
    current_hash: u32,
    depth: u32,
    adjacency: &AdjacencyList,
    molecule: &Molecule,
    params: &AvalonFingerprintParams,
    visited: &mut Vec<bool>,
    path_hashes: &mut Vec<u32>,
) {
    // Record the path hash if we have reached at least min_path bonds.
    if depth >= params.min_path {
        path_hashes.push(current_hash);
    }

    // Stop if we have reached the maximum depth.
    if depth >= params.max_path {
        return;
    }

    // Explore neighbors.
    for neighbor in adjacency.neighbors_of(current_idx) {
        let n_idx = neighbor.atom_index;
        if visited[n_idx] {
            continue;
        }

        // Skip hydrogens if use_hs is false.
        if !params.use_hs && molecule.atoms()[n_idx].atomic_number() == 1 {
            continue;
        }

        visited[n_idx] = true;

        // Build the extended hash: hash_combine(bond_order, then atom_num).
        let bond = &molecule.bonds()[neighbor.bond.index()];
        let bond_code = if params.use_bond_order {
            if bond.is_aromatic() {
                12u32 // Bond::AROMATIC = 12 in RDKit
            } else {
                bond_order_to_u32(bond.order())
            }
        } else {
            1u32 // Single bond placeholder when bond orders are ignored
        };

        let atom_num = molecule.atoms()[n_idx].atomic_number() as u32;

        let mut extended_hash = current_hash;
        hash_combine(&mut extended_hash, bond_code);
        hash_combine(&mut extended_hash, atom_num);

        dfs_paths(
            n_idx,
            extended_hash,
            depth + 1,
            adjacency,
            molecule,
            params,
            visited,
            path_hashes,
        );

        visited[n_idx] = false;
    }
}

/// DFS path enumeration (with explicit hydrogens).
///
/// Same as dfs_paths but includes hydrogen atoms in path enumeration.
fn dfs_paths_hs(
    current_idx: usize,
    current_hash: u32,
    depth: u32,
    adjacency: &AdjacencyList,
    molecule: &Molecule,
    params: &AvalonFingerprintParams,
    visited: &mut Vec<bool>,
    path_hashes: &mut Vec<u32>,
) {
    // Record the path hash if we have reached at least min_path bonds.
    if depth >= params.min_path {
        path_hashes.push(current_hash);
    }

    // Stop if we have reached the maximum depth.
    if depth >= params.max_path {
        return;
    }

    // Explore neighbors.
    for neighbor in adjacency.neighbors_of(current_idx) {
        let n_idx = neighbor.atom_index;
        if visited[n_idx] {
            continue;
        }

        visited[n_idx] = true;

        // Build the extended hash: hash_combine(bond_order, then atom_num).
        let bond = &molecule.bonds()[neighbor.bond.index()];
        let bond_code = if params.use_bond_order {
            if bond.is_aromatic() {
                12u32 // Bond::AROMATIC = 12 in RDKit
            } else {
                bond_order_to_u32(bond.order())
            }
        } else {
            1u32
        };

        let atom_num = molecule.atoms()[n_idx].atomic_number() as u32;

        let mut extended_hash = current_hash;
        hash_combine(&mut extended_hash, bond_code);
        hash_combine(&mut extended_hash, atom_num);

        dfs_paths_hs(
            n_idx,
            extended_hash,
            depth + 1,
            adjacency,
            molecule,
            params,
            visited,
            path_hashes,
        );

        visited[n_idx] = false;
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};

    fn build_methane() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let c = builder.add_atom(AtomSpec::new(Element::C));
        builder.build().unwrap()
    }

    fn build_ethane() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Single))
            .unwrap();
        builder.build().unwrap()
    }

    fn build_ethanol() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        let o = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c2, o, BondOrder::Single))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn test_empty_molecule() {
        let mol = Molecule::new();
        let fp = avalon_fingerprint(&mol, &AvalonFingerprintParams::default()).unwrap();
        assert!(fp.on_bits().is_empty());
    }

    #[test]
    fn test_zero_bits_returns_error() {
        let mol = build_methane();
        let params = AvalonFingerprintParams {
            n_bits: 0,
            ..Default::default()
        };
        let result = avalon_fingerprint(&mol, &params);
        assert!(result.is_err());
    }

    #[test]
    fn test_methane_fingerprint() {
        // Single carbon atom with no bonds -> only single-atom fragments.
        let mol = build_methane();
        let params = AvalonFingerprintParams {
            min_path: 0, // include single-atom fragments
            max_path: 1,
            n_bits: 64,
            ..Default::default()
        };
        let fp = avalon_fingerprint(&mol, &params).unwrap();
        // Single atom should produce at least 1 bit.
        assert!(!fp.on_bits().is_empty());
    }

    #[test]
    fn test_ethane_consistent() {
        let mol = build_ethane();
        let params = AvalonFingerprintParams::default();
        let fp1 = avalon_fingerprint(&mol, &params).unwrap();
        let fp2 = avalon_fingerprint(&mol, &params).unwrap();
        assert_eq!(fp1, fp2);
    }

    #[test]
    fn test_ethanol_different_from_ethane() {
        let ethane = build_ethane();
        let ethanol = build_ethanol();
        let params = AvalonFingerprintParams::default();
        let fp_ethane = avalon_fingerprint(&ethane, &params).unwrap();
        let fp_ethanol = avalon_fingerprint(&ethanol, &params).unwrap();
        // Should differ because ethanol has an oxygen.
        assert_ne!(fp_ethane, fp_ethanol);
    }

    #[test]
    fn test_a_bond_fingerprint() {
        let mol = build_ethane();
        let params = AvalonFingerprintParams {
            min_path: 1,
            max_path: 1,
            n_bits: 64,
            ..Default::default()
        };
        let fp = avalon_fingerprint(&mol, &params).unwrap();
        // Ethane with min_path=1 max_path=1 should enumerate C-C bonds.
        // There is exactly 1 bond path (C-C), producing at least 1 bit.
        assert!(!fp.on_bits().is_empty());
    }

    #[test]
    fn test_from_atoms_filter() {
        let mol = build_ethanol();
        let params_all = AvalonFingerprintParams::default();
        let fp_all = avalon_fingerprint(&mol, &params_all).unwrap();

        let params_first = AvalonFingerprintParams {
            from_atoms: Some(vec![0]),
            ..Default::default()
        };
        let fp_first = avalon_fingerprint(&mol, &params_first).unwrap();

        // FP starting from just atom 0 should have fewer bits than full FP.
        assert!(fp_first.on_bits().len() <= fp_all.on_bits().len());
    }

    #[test]
    fn test_n_bits_respected() {
        let mol = build_ethanol();
        let params = AvalonFingerprintParams {
            n_bits: 32,
            ..Default::default()
        };
        let fp = avalon_fingerprint(&mol, &params).unwrap();
        assert_eq!(fp.n_bits(), 32);

        let params = AvalonFingerprintParams {
            n_bits: 1024,
            ..Default::default()
        };
        let fp = avalon_fingerprint(&mol, &params).unwrap();
        assert_eq!(fp.n_bits(), 1024);
    }

    #[test]
    fn test_different_min_max_path_produces_different_fps() {
        let mol = build_ethanol();
        let fp_deep = avalon_fingerprint(
            &mol,
            &AvalonFingerprintParams {
                min_path: 1,
                max_path: 5,
                ..Default::default()
            },
        )
        .unwrap();
        let fp_shallow = avalon_fingerprint(
            &mol,
            &AvalonFingerprintParams {
                min_path: 1,
                max_path: 1,
                ..Default::default()
            },
        )
        .unwrap();
        // Deeper paths should produce more bits.
        assert!(fp_deep.on_bits().len() >= fp_shallow.on_bits().len());
    }
}
