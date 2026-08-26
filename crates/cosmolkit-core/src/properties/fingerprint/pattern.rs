use std::sync::OnceLock;

use super::{
    Fingerprint, FingerprintError, SsMatcher, hash_combine, rdkit_bond_type_code,
    rdkit_fp_bond_between_atoms,
};
use crate::search::query::{
    AtomQueryPredicate, BondQueryPredicate, QueryNode, build_query_match_context,
    is_complex_atom_query,
};
use crate::search::substruct::{
    SubstructMatchParams, try_get_substruct_matches_with_params_and_context,
};
use crate::{Bond, BondOrder, Molecule};

// RDKit✔️✔️: const std::string PatternFingerprintMolVersion = "1.0.0";
pub const PATTERN_FINGERPRINT_VERSION: &str = "1.0.0";

/// Parameters for the source-compatible Pattern fingerprint.
///
/// The default is a 2,048-bit ordinary Pattern fingerprint. Set `tautomeric`
/// to reproduce the source's tautomer-aware structural hashing. RDKit labels
/// Pattern fingerprint version 1.0.0 experimental; COSMolKit preserves that
/// upstream metadata while validating the modeled ordinary-molecule boundary
/// exactly against the pinned source revision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PatternFingerprintParams {
    /// Number of bits in the explicit result. Zero is rejected.
    pub n_bits: usize,
    /// Whether single, double, and aromatic bonds use tautomer-aware hashing.
    pub tautomeric: bool,
}

impl Default for PatternFingerprintParams {
    fn default() -> Self {
        // RDKit✔️✔️: const ROMol &mol, unsigned int fpSize = 2048,
        // RDKit✔️✔️: std::vector<unsigned int> *atomCounts = nullptr,
        // RDKit✔️✔️: ExplicitBitVect *setOnlyBits = nullptr, bool tautomericFingerprint = false);
        Self {
            n_bits: 2048,
            tautomeric: false,
        }
    }
}

/// Computes an explicit Pattern fingerprint without mutating `molecule`.
///
/// Query-bearing molecules follow RDKit's Pattern-specific atom and bond
/// suppression rules. The implementation compiles the fixed 13-query table
/// once and shares it across scalar, batch, and concurrent calls.
///
/// RDKit's ordinary overload also accepts `atomCounts` and `setOnlyBits`, but
/// the pinned implementation only validates their sizes and otherwise leaves
/// both arguments inert. They are intentionally absent here instead of being
/// exposed with behavior that the source does not implement. The separate
/// RDKit `MolBundle` intersection overload is outside this API because a
/// bundle is not equivalent to an ordered molecule batch.
///
/// # Errors
///
/// Returns [`FingerprintError::EmptyFingerprint`] when `n_bits` is zero and
/// preserves structured SMARTS or substructure failures from the shared
/// Pattern core.
///
/// # Example
///
/// ```
/// use cosmolkit_core::{Molecule, PatternFingerprintParams};
///
/// let molecule = Molecule::from_smiles("c1ccccc1O")?;
/// let fingerprint = molecule.pattern_fingerprint(&PatternFingerprintParams {
///     n_bits: 2048,
///     tautomeric: false,
/// })?;
/// assert_eq!(fingerprint.n_bits(), 2048);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn pattern_fingerprint(
    molecule: &Molecule,
    params: &PatternFingerprintParams,
) -> Result<Fingerprint, FingerprintError> {
    // RDKit✔️❌: ExplicitBitVect *PatternFingerprintMol(const ROMol &mol, unsigned int fpSize,
    // RDKit✔️❌:                                        std::vector<unsigned int> *atomCounts,
    // RDKit✔️❌:                                        ExplicitBitVect *setOnlyBits,
    // RDKit✔️❌:                                        bool tautomericFingerprint) {
    // RDKit✔️❌:   PRECONDITION(fpSize != 0, "fpSize==0");
    // RDKit✔️❌:   PRECONDITION(!atomCounts || atomCounts->size() >= mol.getNumAtoms(),
    // RDKit✔️❌:                "bad atomCounts size");
    // RDKit✔️❌:   PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
    // RDKit✔️❌:                "bad setOnlyBits size");
    // RDKit✔️❌:   auto *res = new ExplicitBitVect(fpSize);
    // RDKit✔️❌:   updatePatternFingerprint(mol, *res, fpSize, atomCounts, setOnlyBits,
    // RDKit✔️❌:                            tautomericFingerprint);
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    //
    // The project-native facade omits the two source arguments proven inert
    // for this overload. Allocation and dispatch remain one zeroed explicit
    // vector plus the sole Pattern core. Its performance marker follows that
    // core until the canonical matcher allocation gap is resolved.
    if params.n_bits == 0 {
        return Err(FingerprintError::EmptyFingerprint);
    }
    let mut fingerprint = Fingerprint::zeroed(params.n_bits);
    update_pattern_fingerprint(
        molecule,
        &mut fingerprint,
        params.n_bits,
        None,
        None,
        params.tautomeric,
    )?;
    Ok(fingerprint)
}

// RDKit✔️✔️: const char *pqs[] = {
// RDKit✔️✔️:     "[*]~[*]", "[*]~[*]~[*]", "[R]~1~[R]~[R]~1",
// RDKit✔️✔️:     //"[*]~[*]~[*]~[*]",
// RDKit✔️✔️:     "[*]~[*](~[*])~[*]",
// RDKit✔️✔️:     //"[*]~[R]~1[R]~[R]~1",
// RDKit✔️✔️:     "[R]~1[R]~[R]~[R]~1",
// RDKit✔️✔️:     //"[*]~[*]~[*]~[*]~[*]",
// RDKit✔️✔️:     "[*]~[*]~[*](~[*])~[*]",
// RDKit✔️✔️:     //"[*]~[R]~1[R]~[R]~1~[*]",
// RDKit✔️✔️:     "[R]~1~[R]~[R]~[R]~[R]~1", "[R]~1~[R]~[R]~[R]~[R]~[R]~1",
// RDKit✔️✔️:     //"[R2]~[R1]~[R2]", Github #151: can't have ring counts in an SSS pattern
// RDKit✔️✔️:     //"[R2]~[R1]~[R1]~[R2]",  Github #151: can't have ring counts in an SSS
// RDKit✔️✔️:     // pattern
// RDKit✔️✔️:     "[R](@[R])(@[R])~[R]~[R](@[R])(@[R])",
// RDKit✔️✔️:     "[R](@[R])(@[R])~[R]@[R]~[R](@[R])(@[R])",
// RDKit✔️✔️:
// RDKit✔️✔️:     //"[*]!@[R]~[R]!@[*]",  Github #151: can't have !@ in an SSS pattern
// RDKit✔️✔️:     //"[*]!@[R]~[R]~[R]!@[*]", Github #151: can't have !@ in an SSS pattern
// RDKit✔️✔️:     "[*]~[R](@[R])@[R](@[R])~[*]", "[*]~[R](@[R])@[R]@[R](@[R])~[*]",
// RDKit✔️✔️:     "[*]",  // single atom fragment
// RDKit✔️✔️:     ""};
//
// The empty C-string is a loop sentinel rather than an active matcher. A
// fixed Rust slice preserves the 13 active entries and their one-based source
// order without retaining a runtime sentinel.
pub(super) const PATTERN_FINGERPRINT_SMARTS: [&str; 13] = [
    "[*]~[*]",
    "[*]~[*]~[*]",
    "[R]~1~[R]~[R]~1",
    "[*]~[*](~[*])~[*]",
    "[R]~1[R]~[R]~[R]~1",
    "[*]~[*]~[*](~[*])~[*]",
    "[R]~1~[R]~[R]~[R]~[R]~1",
    "[R]~1~[R]~[R]~[R]~[R]~[R]~1",
    "[R](@[R])(@[R])~[R]~[R](@[R])(@[R])",
    "[R](@[R])(@[R])~[R]@[R]~[R](@[R])(@[R])",
    "[*]~[R](@[R])@[R](@[R])~[*]",
    "[*]~[R](@[R])@[R]@[R](@[R])~[*]",
    "[*]",
];

pub(super) fn compiled_pattern_fingerprint_queries()
-> Result<&'static [SsMatcher], FingerprintError> {
    // RDKit✔️✔️: typedef boost::flyweight<boost::flyweights::key_value<std::string, ss_matcher>,
    // RDKit✔️✔️:                          boost::flyweights::no_tracking>
    // RDKit✔️✔️:     pattern_flyweight;
    //
    // `OnceLock` provides process-lifetime, thread-safe, no-eviction ownership
    // for this fixed table. Compilation is linear once; later access is O(1)
    // and performs no parsing, allocation, molecule clone, or query clone.
    static CACHE: OnceLock<Result<Vec<SsMatcher>, FingerprintError>> = OnceLock::new();
    match CACHE.get_or_init(|| {
        PATTERN_FINGERPRINT_SMARTS
            .iter()
            .map(|pattern| SsMatcher::try_new(pattern))
            .collect()
    }) {
        Ok(matchers) => Ok(matchers.as_slice()),
        Err(error) => Err(error.clone()),
    }
}

#[inline]
pub(super) fn is_pattern_complex_query(bond: &Bond) -> bool {
    // RDKit✔️✔️: bool isPatternComplexQuery(const Bond *b) {
    // RDKit✔️✔️:   if (!b->hasQuery()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // negated things are always complex:
    // RDKit✔️✔️:   if (b->getQuery()->getNegation()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::string descr = b->getQuery()->getDescription();
    // RDKit✔️✔️:   // std::cerr<<"   !!!!!! "<<b->getIdx()<<"
    // RDKit✔️✔️:   // "<<b->getBeginAtomIdx()<<"-"<<b->getEndAtomIdx()<<" "<<descr<<std::endl;
    // RDKit✔️✔️:   return descr != "BondOrder";
    // RDKit✔️✔️: }
    //
    // `Order` is the sole typed identity for RDKit's `BondOrder`
    // description. A root `Not` is the typed equivalent of the source query's
    // independent negation flag. Local complexity review: both versions make
    // one optional/root query inspection in O(1), without traversal,
    // allocation, cloning, keyed lookup, or molecule access.
    match bond.query() {
        None => false,
        Some(QueryNode::Not(_)) => true,
        Some(QueryNode::Predicate(BondQueryPredicate::Order(_))) => false,
        Some(_) => true,
    }
}

#[inline]
pub(super) fn is_tautomer_bond_query(bond: &Bond) -> bool {
    // RDKit✔️✔️: bool isTautomerBondQuery(const Bond *b) {
    // RDKit✔️✔️:   // assumes we have already tested true for isPatternComplexQuery
    // RDKit✔️✔️:   auto description = b->getQuery()->getDescription();
    // RDKit✔️✔️:   return description == "SingleOrDoubleOrAromaticBond" ||
    // RDKit✔️✔️:          description == "SingleOrAromaticBond";
    // RDKit✔️✔️: }
    //
    // RDKit stores negation beside the query description, so a root `Not`
    // does not change the description inspected here. `OrderIn` is the sole
    // typed fixed-order-set representation, with source-order vectors
    // preserving the two named descriptions. Local complexity review: both
    // implementations perform O(1) root inspection and at most three O(1)
    // enum comparisons without traversal, allocation, cloning, or lookup.
    let query = match bond.query() {
        Some(QueryNode::Not(child)) => child.as_ref(),
        Some(query) => query,
        None => return false,
    };
    matches!(
        query,
        QueryNode::Predicate(BondQueryPredicate::OrderIn(orders))
            if orders.as_slice() == [BondOrder::Single, BondOrder::Aromatic]
                || orders.as_slice()
                    == [BondOrder::Single, BondOrder::Double, BondOrder::Aromatic]
    )
}

fn update_pattern_fingerprint(
    molecule: &Molecule,
    fingerprint: &mut Fingerprint,
    fingerprint_size: usize,
    atom_counts: Option<&mut [u32]>,
    set_only_bits: Option<&Fingerprint>,
    tautomeric_fingerprint: bool,
) -> Result<(), FingerprintError> {
    update_pattern_fingerprint_impl(
        molecule,
        fingerprint,
        fingerprint_size,
        atom_counts,
        set_only_bits,
        tautomeric_fingerprint,
        |_| {},
    )
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum PatternTraceEvent {
    PatternMatches {
        pattern_index: u32,
        count: usize,
    },
    CountBit {
        pattern_index: u32,
        occurrence: usize,
        seed: u32,
        bit: usize,
    },
    AtomHash {
        pattern_index: u32,
        atomic_number: u32,
        seed: u32,
    },
    BondHash {
        pattern_index: u32,
        bond_code: u32,
        seed: u32,
    },
    QueryAtomSuppressed {
        pattern_index: u32,
        atom_index: usize,
    },
    QueryBondSuppressed {
        pattern_index: u32,
        bond_index: usize,
    },
    TautomerQueryBond {
        pattern_index: u32,
        bond_index: usize,
    },
    StructureBit {
        pattern_index: u32,
        seed: u32,
        bit: usize,
    },
    TautomerBit {
        pattern_index: u32,
        seed: u32,
        bit: usize,
    },
}

#[allow(clippy::too_many_arguments)]
fn update_pattern_fingerprint_impl(
    molecule: &Molecule,
    fingerprint: &mut Fingerprint,
    fingerprint_size: usize,
    atom_counts: Option<&mut [u32]>,
    set_only_bits: Option<&Fingerprint>,
    tautomeric_fingerprint: bool,
    mut trace: impl FnMut(PatternTraceEvent),
) -> Result<(), FingerprintError> {
    // RDKit✔️❌: void updatePatternFingerprint(const ROMol &mol, ExplicitBitVect &fp,
    // RDKit✔️❌:                               unsigned int fpSize,
    // RDKit✔️❌:                               std::vector<unsigned int> *atomCounts,
    // RDKit✔️❌:                               ExplicitBitVect *setOnlyBits,
    // RDKit✔️❌:                               bool tautomericFingerprint) {
    // RDKit✔️❌:   PRECONDITION(fpSize != 0, "fpSize==0");
    // RDKit✔️❌:   PRECONDITION(!atomCounts || atomCounts->size() >= mol.getNumAtoms(),
    // RDKit✔️❌:                "bad atomCounts size");
    // RDKit✔️❌:   PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
    // RDKit✔️❌:                "bad setOnlyBits size");
    // RDKit✔️❌:
    // RDKit✔️❌:   std::vector<const ROMol *> patts;
    // RDKit✔️❌:   patts.reserve(10);
    // RDKit✔️❌:   unsigned int idx = 0;
    // RDKit✔️❌:   while (1) {
    // RDKit✔️❌:     std::string pq = pqs[idx];
    // RDKit✔️❌:     if (pq == "") {
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ++idx;
    // RDKit✔️❌:     const ROMol *matcher = pattern_flyweight(pq).get().getMatcher();
    // RDKit✔️❌:     CHECK_INVARIANT(matcher, "bad smarts");
    // RDKit✔️❌:     patts.push_back(matcher);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️❌:     MolOps::fastFindRings(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   boost::dynamic_bitset<> isQueryAtom(mol.getNumAtoms()),
    // RDKit✔️❌:       isQueryBond(mol.getNumBonds()), isTautomerBond(mol.getNumBonds());
    // RDKit✔️❌:   for (const auto at : mol.atoms()) {
    // RDKit✔️❌:     // isComplexQuery() no longer considers "AtomNull" to be complex, but for
    // RDKit✔️❌:     // the purposes of the pattern FP, it definitely needs to be treated as a
    // RDKit✔️❌:     // query feature.
    // RDKit✔️❌:     if (at->hasQuery() && (at->getQuery()->getDescription() == "AtomNull" ||
    // RDKit✔️❌:                            isComplexQuery(at))) {
    // RDKit✔️❌:       isQueryAtom.set(at->getIdx());
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (const auto bond : mol.bonds()) {
    // RDKit✔️❌:     if (isPatternComplexQuery(bond)) {
    // RDKit✔️❌:       isQueryBond.set(bond->getIdx());
    // RDKit✔️❌:       if (tautomericFingerprint && isTautomerBondQuery(bond)) {
    // RDKit✔️❌:         isTautomerBond.set(bond->getIdx());
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   unsigned int pIdx = 0;
    // RDKit✔️❌:   for (const auto patt : patts) {
    // RDKit✔️❌:     ++pIdx;
    // RDKit✔️❌:     std::vector<MatchVectType> matches;
    // RDKit✔️❌:     // uniquify matches?
    // RDKit✔️❌:     //   time for 10K molecules w/ uniquify: 5.24s
    // RDKit✔️❌:     //   time for 10K molecules w/o uniquify: 4.87s
    // RDKit✔️❌:
    // RDKit✔️❌:     SubstructMatchParameters params;
    // RDKit✔️❌:     params.uniquify = false;
    // RDKit✔️❌:     // raise maxMatches really high. This was the cause for github #2614.
    // RDKit✔️❌:     // if we end up with more matches than this, we're completely hosed: :-)
    // RDKit✔️❌:     params.maxMatches = 100000000;
    // RDKit✔️❌:     matches = SubstructMatch(mol, *patt, params);
    // RDKit✔️❌:
    // RDKit✔️❌:     std::uint32_t mIdx = pIdx + patt->getNumAtoms() + patt->getNumBonds();
    // RDKit✔️❌:     for (const auto &mv : matches) {
    // RDKit✔️❌:       // collect bits counting the number of occurrences of the pattern:
    // RDKit✔️❌:       gboost::hash_combine(mIdx, 0xBEEF);
    // RDKit✔️❌:       fp.setBit(mIdx % fpSize);
    // RDKit✔️❌:
    // RDKit✔️❌:       bool isQuery = false;
    // RDKit✔️❌:       std::uint32_t bitId = pIdx;
    // RDKit✔️❌:       std::vector<unsigned int> amap(mv.size(), 0);
    // RDKit✔️❌:       for (const auto &p : mv) {
    // RDKit✔️❌:         if (isQueryAtom[p.second]) {
    // RDKit✔️❌:           isQuery = true;
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         gboost::hash_combine(bitId,
    // RDKit✔️❌:                              mol.getAtomWithIdx(p.second)->getAtomicNum());
    // RDKit✔️❌:         amap[p.first] = p.second;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (isQuery) {
    // RDKit✔️❌:         continue;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       auto tautomerBitId = bitId;
    // RDKit✔️❌:       auto tautomerQuery = false;
    // RDKit✔️❌:       ROMol::EDGE_ITER firstB, lastB;
    // RDKit✔️❌:       boost::tie(firstB, lastB) = patt->getEdges();
    // RDKit✔️❌:       while (!isQuery && firstB != lastB) {
    // RDKit✔️❌:         const Bond *pbond = (*patt)[*firstB];
    // RDKit✔️❌:         ++firstB;
    // RDKit✔️❌:         const Bond *mbond = mol.getBondBetweenAtoms(
    // RDKit✔️❌:             amap[pbond->getBeginAtomIdx()], amap[pbond->getEndAtomIdx()]);
    // RDKit✔️❌:         const auto bondIdx = mbond->getIdx();
    // RDKit✔️❌:
    // RDKit✔️❌:         if (isQueryBond[bondIdx]) {
    // RDKit✔️❌:           isQuery = true;
    // RDKit✔️❌:           if (isTautomerBond[bondIdx]) {
    // RDKit✔️❌:             isQuery = false;
    // RDKit✔️❌:             tautomerQuery = true;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           if (isQuery) {
    // RDKit✔️❌:             break;
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:
    // RDKit✔️❌:         if (tautomericFingerprint) {
    // RDKit✔️❌:           if (isTautomerBond[bondIdx] || mbond->getIsAromatic() ||
    // RDKit✔️❌:               mbond->getBondType() == Bond::SINGLE ||
    // RDKit✔️❌:               mbond->getBondType() == Bond::DOUBLE ||
    // RDKit✔️❌:               mbond->getBondType() == Bond::AROMATIC) {
    // RDKit✔️❌:             gboost::hash_combine(tautomerBitId, -1);
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:
    // RDKit✔️❌:         if (!tautomerQuery) {
    // RDKit✔️❌:           if (!mbond->getIsAromatic()) {
    // RDKit✔️❌:             gboost::hash_combine(bitId, (std::uint32_t)mbond->getBondType());
    // RDKit✔️❌:           } else {
    // RDKit✔️❌:             gboost::hash_combine(bitId, (std::uint32_t)Bond::AROMATIC);
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:
    // RDKit✔️❌:       if (!isQuery) {
    // RDKit✔️❌:         if (!tautomerQuery) {
    // RDKit✔️❌:           fp.setBit(bitId % fpSize);
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (tautomericFingerprint) {
    // RDKit✔️❌:           fp.setBit(tautomerBitId % fpSize);
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    //
    // Behavioral mapping follows the source state machine line-for-line. The
    // canonical matcher receives one target query context shared across all
    // 13 immutable patterns, matching RDKit's single prepared ROMol state.
    // Both paths perform one target-state preparation, 13 VF2 calls, and one
    // pass over every non-unique result; neither reparses or clones a query or
    // target molecule. Hashing and bit insertion are O(1) per source event.
    if fingerprint_size == 0 {
        return Err(FingerprintError::EmptyFingerprint);
    }
    if fingerprint.n_bits() != fingerprint_size {
        return Err(FingerprintError::BitLengthMismatch {
            left: fingerprint.n_bits(),
            right: fingerprint_size,
        });
    }
    if atom_counts
        .as_ref()
        .is_some_and(|counts| counts.len() < molecule.num_atoms())
    {
        return Err(FingerprintError::InvalidArguments {
            reason: "Pattern atom_counts length is smaller than molecule atom count",
        });
    }
    if set_only_bits.is_some_and(|bits| bits.n_bits() != fingerprint_size) {
        return Err(FingerprintError::InvalidArguments {
            reason: "Pattern set_only_bits length differs from fingerprint size",
        });
    }

    let patterns = compiled_pattern_fingerprint_queries()?;
    let is_query_atom: Vec<bool> = molecule
        .atoms()
        .iter()
        .map(|atom| {
            atom.query().is_some()
                && (matches!(
                    atom.query(),
                    Some(QueryNode::Predicate(AtomQueryPredicate::Any))
                ) || is_complex_atom_query(atom))
        })
        .collect();
    let mut is_query_bond = vec![false; molecule.num_bonds()];
    let mut is_tautomer_bond = vec![false; molecule.num_bonds()];
    for bond in molecule.bonds() {
        if is_pattern_complex_query(bond) {
            is_query_bond[bond.id().index()] = true;
            if tautomeric_fingerprint && is_tautomer_bond_query(bond) {
                is_tautomer_bond[bond.id().index()] = true;
            }
        }
    }

    let mut params = SubstructMatchParams::default();
    params.uniquify = false;
    params.max_matches = 100_000_000;
    let query_context = build_query_match_context(molecule);
    for (pattern_offset, matcher) in patterns.iter().enumerate() {
        let pattern_index = pattern_offset as u32 + 1;
        let pattern = matcher.getMatcher();
        let matches = try_get_substruct_matches_with_params_and_context(
            molecule,
            pattern,
            &params,
            &query_context,
        )
        .map_err(|error| FingerprintError::Pattern {
            reason: error.to_string(),
        })?;
        trace(PatternTraceEvent::PatternMatches {
            pattern_index,
            count: matches.len(),
        });

        let mut match_index = pattern_index
            .wrapping_add(pattern.num_atoms() as u32)
            .wrapping_add(pattern.num_bonds() as u32);
        for (occurrence, matched) in matches.into_iter().enumerate() {
            hash_combine(&mut match_index, 0xBEEF);
            let count_bit = match_index as usize % fingerprint_size;
            fingerprint.set_bit(count_bit);
            trace(PatternTraceEvent::CountBit {
                pattern_index,
                occurrence,
                seed: match_index,
                bit: count_bit,
            });

            let mut is_query = false;
            let mut bit_id = pattern_index;
            let mut atom_map = vec![0usize; matched.atom_mapping.len()];
            for (query_atom_index, &molecule_atom_index) in matched.atom_mapping.iter().enumerate()
            {
                if is_query_atom[molecule_atom_index] {
                    is_query = true;
                    trace(PatternTraceEvent::QueryAtomSuppressed {
                        pattern_index,
                        atom_index: molecule_atom_index,
                    });
                    break;
                }
                let atomic_number =
                    u32::from(molecule.atoms()[molecule_atom_index].atomic_number());
                hash_combine(&mut bit_id, atomic_number);
                trace(PatternTraceEvent::AtomHash {
                    pattern_index,
                    atomic_number,
                    seed: bit_id,
                });
                atom_map[query_atom_index] = molecule_atom_index;
            }
            if is_query {
                continue;
            }

            let mut tautomer_bit_id = bit_id;
            let mut tautomer_query = false;
            for pattern_bond in pattern.bonds() {
                let molecule_bond_index = rdkit_fp_bond_between_atoms(
                    molecule,
                    atom_map[pattern_bond.begin().index()],
                    atom_map[pattern_bond.end().index()],
                )
                .expect("substructure atom mapping must preserve query bonds");
                let molecule_bond = &molecule.bonds()[molecule_bond_index];

                if is_query_bond[molecule_bond_index] {
                    is_query = true;
                    if is_tautomer_bond[molecule_bond_index] {
                        is_query = false;
                        tautomer_query = true;
                        trace(PatternTraceEvent::TautomerQueryBond {
                            pattern_index,
                            bond_index: molecule_bond_index,
                        });
                    }
                    if is_query {
                        trace(PatternTraceEvent::QueryBondSuppressed {
                            pattern_index,
                            bond_index: molecule_bond_index,
                        });
                        break;
                    }
                }

                if tautomeric_fingerprint
                    && (is_tautomer_bond[molecule_bond_index]
                        || molecule_bond.is_aromatic()
                        || matches!(
                            molecule_bond.order(),
                            BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
                        ))
                {
                    hash_combine(&mut tautomer_bit_id, u32::MAX);
                }

                if !tautomer_query {
                    let bond_code = if molecule_bond.is_aromatic() {
                        rdkit_bond_type_code(BondOrder::Aromatic)
                    } else {
                        rdkit_bond_type_code(molecule_bond.order())
                    };
                    hash_combine(&mut bit_id, bond_code);
                    trace(PatternTraceEvent::BondHash {
                        pattern_index,
                        bond_code,
                        seed: bit_id,
                    });
                }
            }

            if !is_query {
                if !tautomer_query {
                    let bit = bit_id as usize % fingerprint_size;
                    fingerprint.set_bit(bit);
                    trace(PatternTraceEvent::StructureBit {
                        pattern_index,
                        seed: bit_id,
                        bit,
                    });
                }
                if tautomeric_fingerprint {
                    let bit = tautomer_bit_id as usize % fingerprint_size;
                    fingerprint.set_bit(bit);
                    trace(PatternTraceEvent::TautomerBit {
                        pattern_index,
                        seed: tautomer_bit_id,
                        bit,
                    });
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::sync::{Arc, Barrier};

    use super::*;
    use crate::{AtomId, BondId, BondSpec};

    fn test_bond(query: Option<QueryNode<BondQueryPredicate>>) -> Bond {
        let spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Unspecified);
        let spec = match query {
            Some(query) => spec.with_query(query),
            None => spec,
        };
        Bond::from_spec(BondId::new(0), spec)
    }

    #[test]
    fn pattern_compiled_queries_preserve_source_order_and_share_one_cache() {
        let first = compiled_pattern_fingerprint_queries().expect("compile Pattern SMARTS");
        let second = compiled_pattern_fingerprint_queries().expect("reuse Pattern SMARTS");

        assert_eq!(first.len(), PATTERN_FINGERPRINT_SMARTS.len());
        assert!(std::ptr::eq(first.as_ptr(), second.as_ptr()));
        for (matcher, expected) in first.iter().zip(PATTERN_FINGERPRINT_SMARTS) {
            let reparsed = crate::search::smarts_parse::mol_from_smarts(
                expected,
                &crate::search::smarts_parse::SmartsParseParams::default(),
            )
            .expect("reparse source Pattern SMARTS");
            assert_eq!(matcher.getMatcher().num_atoms(), reparsed.num_atoms());
            assert_eq!(matcher.getMatcher().num_bonds(), reparsed.num_bonds());
        }
    }

    #[test]
    fn pattern_compiled_queries_are_shared_during_concurrent_access() {
        const THREADS: usize = 16;
        let barrier = Arc::new(Barrier::new(THREADS));
        let handles: Vec<_> = (0..THREADS)
            .map(|_| {
                let barrier = Arc::clone(&barrier);
                std::thread::spawn(move || {
                    barrier.wait();
                    let matchers =
                        compiled_pattern_fingerprint_queries().expect("concurrent Pattern cache");
                    (matchers.as_ptr() as usize, matchers.len())
                })
            })
            .collect();
        let results: Vec<_> = handles
            .into_iter()
            .map(|handle| handle.join().expect("Pattern cache thread"))
            .collect();

        assert!(results.iter().all(|result| *result == results[0]));
        assert_eq!(results[0].1, PATTERN_FINGERPRINT_SMARTS.len());
    }

    #[test]
    fn pattern_query_classifiers_reproduce_description_and_negation_semantics() {
        let order = |order| QueryNode::predicate(BondQueryPredicate::Order(order));
        let order_in = |orders| QueryNode::predicate(BondQueryPredicate::OrderIn(orders));

        assert!(!is_pattern_complex_query(&test_bond(None)));
        assert!(!is_pattern_complex_query(&test_bond(Some(order(
            BondOrder::Single,
        )))));
        assert!(is_pattern_complex_query(&test_bond(Some(QueryNode::not(
            order(BondOrder::Single),
        )))));
        assert!(is_pattern_complex_query(&test_bond(Some(order_in(vec![
            BondOrder::Single,
            BondOrder::Aromatic,
        ])))));
        assert!(is_pattern_complex_query(&test_bond(Some(QueryNode::or(
            vec![order(BondOrder::Single), order(BondOrder::Aromatic)],
        )))));

        for query in [
            order_in(vec![BondOrder::Single, BondOrder::Aromatic]),
            order_in(vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Aromatic,
            ]),
            QueryNode::not(order_in(vec![BondOrder::Single, BondOrder::Aromatic])),
        ] {
            assert!(is_tautomer_bond_query(&test_bond(Some(query))));
        }

        for query in [
            order(BondOrder::Single),
            order_in(vec![BondOrder::Aromatic, BondOrder::Single]),
            order_in(vec![BondOrder::Double, BondOrder::Aromatic]),
            order_in(vec![BondOrder::Single, BondOrder::Double]),
            QueryNode::or(vec![order(BondOrder::Single), order(BondOrder::Aromatic)]),
        ] {
            assert!(!is_tautomer_bond_query(&test_bond(Some(query))));
        }
        assert!(!is_tautomer_bond_query(&test_bond(None)));
    }

    fn traced_pattern_fingerprint(
        molecule: &Molecule,
        fingerprint_size: usize,
        tautomeric: bool,
    ) -> (Fingerprint, Vec<PatternTraceEvent>) {
        let mut fingerprint = Fingerprint::zeroed(fingerprint_size);
        let mut events = Vec::new();
        update_pattern_fingerprint_impl(
            molecule,
            &mut fingerprint,
            fingerprint_size,
            None,
            None,
            tautomeric,
            |event| events.push(event),
        )
        .expect("Pattern fingerprint");
        (fingerprint, events)
    }

    #[test]
    fn pattern_core_matches_exact_ethane_hash_evolution_and_tautomer_bits() {
        let molecule = Molecule::from_smiles("CC").expect("ethane");
        let (ordinary, ordinary_events) = traced_pattern_fingerprint(&molecule, 2048, false);
        assert_eq!(ordinary.on_bits(), vec![429, 778, 1022, 1061, 1236, 1295]);
        assert!(ordinary_events.contains(&PatternTraceEvent::CountBit {
            pattern_index: 1,
            occurrence: 0,
            seed: 2_654_484_909,
            bit: 429,
        }));
        assert!(ordinary_events.contains(&PatternTraceEvent::CountBit {
            pattern_index: 1,
            occurrence: 1,
            seed: 3_454_831_614,
            bit: 1022,
        }));
        assert!(ordinary_events.contains(&PatternTraceEvent::AtomHash {
            pattern_index: 1,
            atomic_number: 6,
            seed: 2_654_435_838,
        }));
        assert!(ordinary_events.contains(&PatternTraceEvent::BondHash {
            pattern_index: 1,
            bond_code: 1,
            seed: 4_217_150_218,
        }));
        assert!(ordinary_events.contains(&PatternTraceEvent::StructureBit {
            pattern_index: 1,
            seed: 4_217_150_218,
            bit: 778,
        }));

        let (tautomeric, tautomeric_events) = traced_pattern_fingerprint(&molecule, 2048, true);
        assert_eq!(
            tautomeric.on_bits(),
            vec![429, 776, 778, 1022, 1061, 1236, 1295]
        );
        assert!(tautomeric_events.contains(&PatternTraceEvent::TautomerBit {
            pattern_index: 1,
            seed: 4_217_150_216,
            bit: 776,
        }));
    }

    #[test]
    fn pattern_core_exercises_every_source_pattern_with_non_unique_matches() {
        let fixtures = [
            "CC",
            "CCC",
            "C1CC1",
            "CC(C)C",
            "C1CCC1",
            "CC(C)CC",
            "C1CCCC1",
            "C1CCCCC1",
            "C1CC2CCC1C2",
            "C12C3C4C1C5C2C3C45",
            "c1ccc2ccccc2c1",
            "C1C2CC3CC1CC(C2)C3",
            "[Na+]",
        ];
        let mut maximum_counts = [0usize; 13];
        for smiles in fixtures {
            let molecule = Molecule::from_smiles(smiles).expect("Pattern fixture");
            let (_, events) = traced_pattern_fingerprint(&molecule, 2048, false);
            for event in events {
                if let PatternTraceEvent::PatternMatches {
                    pattern_index,
                    count,
                } = event
                {
                    maximum_counts[pattern_index as usize - 1] =
                        maximum_counts[pattern_index as usize - 1].max(count);
                }
            }
        }
        assert!(
            maximum_counts.iter().all(|&count| count > 0),
            "every source pattern must be exercised: {maximum_counts:?}"
        );
        assert!(
            maximum_counts.iter().all(|&count| count != 1),
            "non-unique matching must preserve symmetry multiplicity: {maximum_counts:?}"
        );
    }

    #[test]
    fn pattern_core_normalizes_aromatic_bonds_and_keeps_single_atom_pattern() {
        let benzene = Molecule::from_smiles("c1ccccc1").expect("benzene");
        let (_, events) = traced_pattern_fingerprint(&benzene, 2048, false);
        assert!(
            events
                .iter()
                .any(|event| matches!(event, PatternTraceEvent::BondHash { bond_code: 12, .. }))
        );

        let sodium = Molecule::from_smiles("[Na+]").expect("sodium");
        let (fingerprint, events) = traced_pattern_fingerprint(&sodium, 2048, false);
        assert!(events.contains(&PatternTraceEvent::PatternMatches {
            pattern_index: 13,
            count: 1,
        }));
        assert!(events.iter().any(|event| matches!(
            event,
            PatternTraceEvent::StructureBit {
                pattern_index: 13,
                ..
            }
        )));
        assert!(!fingerprint.on_bits().is_empty());
    }

    #[test]
    fn pattern_core_suppresses_query_atoms_and_non_tautomer_query_bonds() {
        let query_atom = crate::search::smarts_parse::mol_from_smarts(
            "[*]",
            &crate::search::smarts_parse::SmartsParseParams::default(),
        )
        .expect("query atom");
        let (_, atom_events) = traced_pattern_fingerprint(&query_atom, 2048, false);
        assert!(
            atom_events
                .iter()
                .any(|event| matches!(event, PatternTraceEvent::QueryAtomSuppressed { .. }))
        );
        assert!(!atom_events.iter().any(|event| matches!(
            event,
            PatternTraceEvent::StructureBit {
                pattern_index: 13,
                ..
            }
        )));

        let query_bond = crate::search::smarts_parse::mol_from_smarts(
            "C~C",
            &crate::search::smarts_parse::SmartsParseParams::default(),
        )
        .expect("query bond");
        let (_, bond_events) = traced_pattern_fingerprint(&query_bond, 2048, false);
        assert!(bond_events.iter().any(|event| matches!(
            event,
            PatternTraceEvent::QueryBondSuppressed {
                pattern_index: 1,
                ..
            }
        )));
        assert_eq!(
            traced_pattern_fingerprint(&query_bond, 257, false)
                .0
                .on_bits(),
            vec![14, 44, 132, 136, 146]
        );

        let any_query = crate::search::smarts_parse::mol_from_smarts(
            "C~N",
            &crate::search::smarts_parse::SmartsParseParams::default(),
        )
        .expect("any query bond");
        assert_eq!(
            traced_pattern_fingerprint(&any_query, 257, false)
                .0
                .on_bits(),
            vec![14, 43, 44, 132, 136, 146]
        );

        let order_query = crate::search::smarts_parse::mol_from_smarts(
            "C-,=N",
            &crate::search::smarts_parse::SmartsParseParams::default(),
        )
        .expect("single-or-double query bond");
        assert_eq!(
            traced_pattern_fingerprint(&order_query, 257, false)
                .0
                .on_bits(),
            vec![14, 43, 44, 132, 136, 146]
        );
    }

    #[test]
    fn pattern_core_tautomer_query_uses_u32_max_hash_and_suppresses_structure_bit() {
        let query = crate::search::smarts_parse::mol_from_smarts(
            "CC",
            &crate::search::smarts_parse::SmartsParseParams::default(),
        )
        .expect("single-or-aromatic query bond");
        let (_, events) = traced_pattern_fingerprint(&query, 2048, true);
        assert!(events.iter().any(|event| matches!(
            event,
            PatternTraceEvent::TautomerQueryBond {
                pattern_index: 1,
                ..
            }
        )));
        assert!(!events.iter().any(|event| matches!(
            event,
            PatternTraceEvent::StructureBit {
                pattern_index: 1,
                ..
            }
        )));
        assert!(events.iter().any(|event| matches!(
            event,
            PatternTraceEvent::TautomerBit {
                pattern_index: 1,
                ..
            }
        )));
    }

    #[test]
    fn pattern_core_width_collision_and_inert_arguments_match_source() {
        let molecule = Molecule::from_smiles("CCC").expect("propane");
        let (width_one, _) = traced_pattern_fingerprint(&molecule, 1, true);
        assert_eq!(width_one.on_bits(), vec![0]);

        let mut baseline = Fingerprint::zeroed(127);
        update_pattern_fingerprint(&molecule, &mut baseline, 127, None, None, true)
            .expect("baseline");
        let mut atom_counts = vec![17; molecule.num_atoms()];
        let set_only_bits = Fingerprint::from_on_bits(127, [0, 5, 126]);
        let mut with_inert_arguments = Fingerprint::zeroed(127);
        update_pattern_fingerprint(
            &molecule,
            &mut with_inert_arguments,
            127,
            Some(&mut atom_counts),
            Some(&set_only_bits),
            true,
        )
        .expect("inert arguments");
        assert_eq!(with_inert_arguments, baseline);
        assert_eq!(atom_counts, vec![17; molecule.num_atoms()]);

        let mut invalid = Fingerprint::zeroed(127);
        assert!(matches!(
            update_pattern_fingerprint(
                &molecule,
                &mut invalid,
                127,
                Some(&mut [0; 2]),
                None,
                false,
            ),
            Err(FingerprintError::InvalidArguments { .. })
        ));
        assert!(matches!(
            update_pattern_fingerprint(
                &molecule,
                &mut invalid,
                127,
                None,
                Some(&Fingerprint::zeroed(128)),
                false,
            ),
            Err(FingerprintError::InvalidArguments { .. })
        ));
        let mut empty = Fingerprint::zeroed(0);
        assert_eq!(
            update_pattern_fingerprint(&molecule, &mut empty, 0, None, None, false),
            Err(FingerprintError::EmptyFingerprint)
        );
    }
}
