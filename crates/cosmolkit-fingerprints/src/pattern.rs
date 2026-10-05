//! Pinned Pattern fingerprint; SEARCH owns matching and query classification.
use crate::{Fingerprint, FingerprintError, hash::hash_combine};
use cosmolkit_core::{RingFindingError, RingInfo};
use cosmolkit_model::{
    AtomQueryPredicate, Bond, BondOrder, BondQueryPredicate, CoordinateBlock, Element,
    QueryAtomConversionError, QueryAtomIdentity, QueryGraph, QueryGraphError, QueryNode,
    TopologyBlock, TopologyValidationError,
};
use cosmolkit_search::{
    CompiledQuery, QueryCompileError, QueryMatchContextError, SearchTarget, SmartsParseError,
    SubstructMatchError, SubstructMatchParams, build_ring_only_query_match_context,
    is_complex_atom_query, is_pattern_complex_query, is_tautomer_bond_query, parse_smarts,
    try_get_substruct_atom_matches_with_compiled_query_and_context,
};
use std::{borrow::Cow, error::Error, fmt, sync::OnceLock};

#[derive(Debug, Clone, PartialEq)]
pub enum PatternFingerprintError {
    EmptyFingerprint,
    InvalidArguments { reason: &'static str },
    BitLengthMismatch { left: usize, right: usize },
    Topology(TopologyValidationError),
    Query(QueryGraphError),
    QueryCarrier(QueryAtomConversionError),
    Rings(RingFindingError),
    Smarts(SmartsParseError),
    QueryCompile(QueryCompileError),
    QueryContext(QueryMatchContextError),
    Match(SubstructMatchError),
    Value(FingerprintError),
}
impl fmt::Display for PatternFingerprintError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyFingerprint => f.write_str("fingerprint requires n_bits > 0"),
            Self::InvalidArguments { reason } => {
                write!(f, "invalid fingerprint arguments: {reason}")
            }
            Self::BitLengthMismatch { left, right } => {
                write!(f, "fingerprint bit lengths differ: {left} and {right}")
            }
            Self::Topology(e) => e.fmt(f),
            Self::Query(e) => e.fmt(f),
            Self::QueryCarrier(e) => e.fmt(f),
            Self::Rings(e) => e.fmt(f),
            Self::Smarts(e) => e.fmt(f),
            Self::QueryCompile(e) => e.fmt(f),
            Self::QueryContext(e) => e.fmt(f),
            Self::Match(e) => e.fmt(f),
            Self::Value(e) => e.fmt(f),
        }
    }
}
impl Error for PatternFingerprintError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Topology(e) => Some(e),
            Self::Query(e) => Some(e),
            Self::QueryCarrier(e) => Some(e),
            Self::Rings(e) => Some(e),
            Self::Smarts(e) => Some(e),
            Self::QueryCompile(e) => Some(e),
            Self::QueryContext(e) => Some(e),
            Self::Match(e) => Some(e),
            Self::Value(e) => Some(e),
            Self::EmptyFingerprint
            | Self::InvalidArguments { .. }
            | Self::BitLengthMismatch { .. } => None,
        }
    }
}
macro_rules! pattern_error_from {
    ($ty:ty, $variant:ident) => {
        impl From<$ty> for PatternFingerprintError {
            fn from(e: $ty) -> Self {
                Self::$variant(e)
            }
        }
    };
}
pattern_error_from!(TopologyValidationError, Topology);
pattern_error_from!(QueryGraphError, Query);
pattern_error_from!(QueryAtomConversionError, QueryCarrier);
pattern_error_from!(RingFindingError, Rings);
pattern_error_from!(SmartsParseError, Smarts);
pattern_error_from!(QueryCompileError, QueryCompile);
pattern_error_from!(QueryMatchContextError, QueryContext);
pattern_error_from!(SubstructMatchError, Match);
pattern_error_from!(FingerprintError, Value);

#[derive(Clone, Copy)]
enum PatternGraphInput<'a> {
    Concrete(&'a TopologyBlock),
    Query(&'a QueryGraph),
}
impl<'a> PatternGraphInput<'a> {
    fn num_atoms(self) -> usize {
        match self {
            Self::Concrete(t) => t.atoms.len(),
            Self::Query(q) => q.num_atoms(),
        }
    }
    fn num_bonds(self) -> usize {
        match self {
            Self::Concrete(t) => t.bonds.len(),
            Self::Query(q) => q.num_bonds(),
        }
    }
    fn atomic_number(self, i: usize) -> u8 {
        match self {
            Self::Concrete(t) => t.atoms[i].atomic_number(),
            Self::Query(q) => q.atoms()[i].atomic_number(),
        }
    }
    fn is_query_atom(self, i: usize) -> bool {
        match self {
            Self::Concrete(_) => false,
            Self::Query(q) => {
                let atom = &q.atoms()[i];
                !atom.predicate_is_carrier_derived()
                    && (matches!(
                        atom.predicate(),
                        QueryNode::Predicate(AtomQueryPredicate::Any)
                    ) || is_complex_atom_query(atom))
            }
        }
    }
    fn bond_query(self, i: usize) -> Option<&'a QueryNode<BondQueryPredicate>> {
        match self {
            Self::Concrete(t) => t.bonds[i].query(),
            Self::Query(q) => {
                let b = &q.bonds()[i];
                (!b.predicate_is_carrier_derived()).then(|| b.predicate())
            }
        }
    }
    fn matching_topology(self) -> Result<Cow<'a, TopologyBlock>, PatternFingerprintError> {
        match self {
            Self::Concrete(t) => {
                t.validate()?;
                Ok(Cow::Borrowed(t))
            }
            Self::Query(q) => {
                q.validate()?;
                // Fixed queries use only Any, ring membership and bond order,
                // with source useQueryQueryMatches=false. Their target access
                // needs carrier attributes, never target ASTs or atomic number.
                // A non-Element numeric query identity is projected explicitly
                // only for this Element-only matcher carrier. The source u8
                // identity remains authoritative in atomic_number() hashing;
                // it is never overwritten in the canonical query value.
                let atoms = q
                    .atoms()
                    .iter()
                    .map(|a| match a.identity() {
                        QueryAtomIdentity::Element(_) => a.try_to_atom(),
                        QueryAtomIdentity::AtomicNumber(_) => a
                            .clone()
                            .with_identity(QueryAtomIdentity::Element(Element::DUMMY))
                            .try_to_atom(),
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let bonds = q.bonds().iter().map(|b| b.bond().clone()).collect();
                Ok(Cow::Owned(TopologyBlock::try_from_parts(
                    atoms,
                    bonds,
                    vec![],
                    vec![],
                )?))
            }
        }
    }
}

// RDKit❗✔️: const std::string PatternFingerprintMolVersion = "1.0.0";
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
        // RDKit❗✔️:     const ROMol &mol, unsigned int fpSize = 2048,
        // RDKit❗✔️:     std::vector<unsigned int> *atomCounts = nullptr,
        // RDKit❗✔️:     ExplicitBitVect *setOnlyBits = nullptr, bool tautomericFingerprint = false);
        Self {
            n_bits: 2048,
            tautomeric: false,
        }
    }
}

/// Compute the pinned experimental Pattern fingerprint from detached topology.
/// The source ordinary overload's atomCounts and setOnlyBits are size-checked
/// but inert; they are private validation inputs, not public configuration.
pub fn pattern_fingerprint(
    topology: &TopologyBlock,
    rings: Option<&RingInfo>,
    params: &PatternFingerprintParams,
) -> Result<Fingerprint, PatternFingerprintError> {
    fingerprint_for_graph(PatternGraphInput::Concrete(topology), rings, params)
}

/// Compute the same source fingerprint for the canonical detached query value.
pub fn pattern_query_fingerprint(
    query: &QueryGraph,
    rings: Option<&RingInfo>,
    params: &PatternFingerprintParams,
) -> Result<Fingerprint, PatternFingerprintError> {
    fingerprint_for_graph(PatternGraphInput::Query(query), rings, params)
}

fn fingerprint_for_graph(
    graph: PatternGraphInput<'_>,
    rings: Option<&RingInfo>,
    params: &PatternFingerprintParams,
) -> Result<Fingerprint, PatternFingerprintError> {
    // RDKit❗❌: ExplicitBitVect *PatternFingerprintMol(const ROMol &mol, unsigned int fpSize,
    // RDKit❗❌:                                        std::vector<unsigned int> *atomCounts,
    // RDKit❗❌:                                        ExplicitBitVect *setOnlyBits,
    // RDKit❗❌:                                        bool tautomericFingerprint) {
    // RDKit❗❌:   PRECONDITION(fpSize != 0, "fpSize==0");
    // RDKit❗❌:   PRECONDITION(!atomCounts || atomCounts->size() >= mol.getNumAtoms(),
    // RDKit❗❌:                "bad atomCounts size");
    // RDKit❗❌:   PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
    // RDKit❗❌:                "bad setOnlyBits size");
    // RDKit❗❌:   auto *res = new ExplicitBitVect(fpSize);
    // RDKit❗❌:   updatePatternFingerprint(mol, *res, fpSize, atomCounts, setOnlyBits,
    // RDKit❗❌:                            tautomericFingerprint);
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    if params.n_bits == 0 {
        return Err(PatternFingerprintError::EmptyFingerprint);
    }
    let width =
        u32::try_from(params.n_bits).map_err(|_| PatternFingerprintError::InvalidArguments {
            reason: "Pattern n_bits exceeds source unsigned int width",
        })?;
    let mut fingerprint = Fingerprint::new(width);
    update_pattern_fingerprint_impl(
        graph,
        rings,
        &mut fingerprint,
        params.n_bits,
        None,
        None,
        params.tautomeric,
        |_| {},
    )?;
    Ok(fingerprint)
}

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
-> Result<&'static [CompiledQuery], PatternFingerprintError> {
    // RDKit❗✔️: const char *pqs[] = {
    // RDKit❗✔️:     "[*]~[*]", "[*]~[*]~[*]", "[R]~1~[R]~[R]~1",
    // RDKit❗✔️:     //"[*]~[*]~[*]~[*]",
    // RDKit❗✔️:     "[*]~[*](~[*])~[*]",
    // RDKit❗✔️:     //"[*]~[R]~1[R]~[R]~1",
    // RDKit❗✔️:     "[R]~1[R]~[R]~[R]~1",
    // RDKit❗✔️:     //"[*]~[*]~[*]~[*]~[*]",
    // RDKit❗✔️:     "[*]~[*]~[*](~[*])~[*]",
    // RDKit❗✔️:     //"[*]~[R]~1[R]~[R]~1~[*]",
    // RDKit❗✔️:     "[R]~1~[R]~[R]~[R]~[R]~1", "[R]~1~[R]~[R]~[R]~[R]~[R]~1",
    // RDKit❗✔️:     //"[R2]~[R1]~[R2]", Github #151: can't have ring counts in an SSS pattern
    // RDKit❗✔️:     //"[R2]~[R1]~[R1]~[R2]",  Github #151: can't have ring counts in an SSS
    // RDKit❗✔️:     // pattern
    // RDKit❗✔️:     "[R](@[R])(@[R])~[R]~[R](@[R])(@[R])",
    // RDKit❗✔️:     "[R](@[R])(@[R])~[R]@[R]~[R](@[R])(@[R])",
    // RDKit❗✔️:
    // RDKit❗✔️:     //"[*]!@[R]~[R]!@[*]",  Github #151: can't have !@ in an SSS pattern
    // RDKit❗✔️:     //"[*]!@[R]~[R]~[R]!@[*]", Github #151: can't have !@ in an SSS pattern
    // RDKit❗✔️:     "[*]~[R](@[R])@[R](@[R])~[*]", "[*]~[R](@[R])@[R]@[R](@[R])~[*]",
    // RDKit❗✔️:     "[*]",  // single atom fragment
    // RDKit❗✔️:     ""};
    //
    // The empty C-string is a loop sentinel rather than an active matcher. A
    // fixed Rust slice preserves the 13 active entries and their one-based source
    // order without retaining a runtime sentinel.
    // RDKit❗✔️: typedef boost::flyweight<boost::flyweights::key_value<std::string, ss_matcher>,
    // RDKit❗✔️:                          boost::flyweights::no_tracking>
    // RDKit❗✔️:     pattern_flyweight;
    //
    // `OnceLock` provides process-lifetime, thread-safe, no-eviction ownership
    // for this fixed table. Compilation is linear once; later access is O(1)
    // and performs no parsing, allocation, molecule clone, or query clone.
    static CACHE: OnceLock<Result<Vec<CompiledQuery>, PatternFingerprintError>> = OnceLock::new();
    match CACHE.get_or_init(|| {
        PATTERN_FINGERPRINT_SMARTS
            .iter()
            .map(|pattern| compile_pattern_query(pattern))
            .collect()
    }) {
        Ok(matchers) => Ok(matchers.as_slice()),
        Err(error) => Err(error.clone()),
    }
}

fn compile_pattern_query(pattern: &str) -> Result<CompiledQuery, PatternFingerprintError> {
    // BEGIN RDKIT CPP FUNCTION ss_matcher::ss_matcher
    // RDKit❗✔️:   ss_matcher(const std::string &pattern) {
    // RDKit❗✔️:     RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
    // RDKit❗✔️:     TEST_ASSERT(p);
    // RDKit❗✔️:     m_matcher.reset(p);
    // RDKit❗✔️:   };
    // END RDKIT CPP FUNCTION ss_matcher::ss_matcher
    // SEARCH owns parsing and compilation once. Its typed failures are retained;
    // no parse/compile, query clone or allocation occurs on warmed cache reads.
    Ok(CompiledQuery::compile(parse_smarts(
        pattern,
        &Default::default(),
    )?)?)
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
    molecule: PatternGraphInput<'_>,
    rings: Option<&RingInfo>,
    fingerprint: &mut Fingerprint,
    fingerprint_size: usize,
    atom_counts: Option<&mut [u32]>,
    set_only_bits: Option<&Fingerprint>,
    tautomeric_fingerprint: bool,
    mut trace: impl FnMut(PatternTraceEvent),
) -> Result<(), PatternFingerprintError> {
    // RDKit❗❌: void updatePatternFingerprint(const ROMol &mol, ExplicitBitVect &fp,
    // RDKit❗❌:                               unsigned int fpSize,
    // RDKit❗❌:                               std::vector<unsigned int> *atomCounts,
    // RDKit❗❌:                               ExplicitBitVect *setOnlyBits,
    // RDKit❗❌:                               bool tautomericFingerprint) {
    // RDKit❗❌:   PRECONDITION(fpSize != 0, "fpSize==0");
    // RDKit❗❌:   PRECONDITION(!atomCounts || atomCounts->size() >= mol.getNumAtoms(),
    // RDKit❗❌:                "bad atomCounts size");
    // RDKit❗❌:   PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
    // RDKit❗❌:                "bad setOnlyBits size");
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<const ROMol *> patts;
    // RDKit❗❌:   patts.reserve(10);
    // RDKit❗❌:   unsigned int idx = 0;
    // RDKit❗❌:   while (1) {
    // RDKit❗❌:     std::string pq = pqs[idx];
    // RDKit❗❌:     if (pq == "") {
    // RDKit❗❌:       break;
    // RDKit❗❌:     }
    // RDKit❗❌:     ++idx;
    // RDKit❗❌:     const ROMol *matcher = pattern_flyweight(pq).get().getMatcher();
    // RDKit❗❌:     CHECK_INVARIANT(matcher, "bad smarts");
    // RDKit❗❌:     patts.push_back(matcher);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit❗❌:     MolOps::fastFindRings(mol);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> isQueryAtom(mol.getNumAtoms()),
    // RDKit❗❌:       isQueryBond(mol.getNumBonds()), isTautomerBond(mol.getNumBonds());
    // RDKit❗❌:   for (const auto at : mol.atoms()) {
    // RDKit❗❌:     // isComplexQuery() no longer considers "AtomNull" to be complex, but for
    // RDKit❗❌:     // the purposes of the pattern FP, it definitely needs to be treated as a
    // RDKit❗❌:     // query feature.
    // RDKit❗❌:     if (at->hasQuery() && (at->getQuery()->getDescription() == "AtomNull" ||
    // RDKit❗❌:                            isComplexQuery(at))) {
    // RDKit❗❌:       isQueryAtom.set(at->getIdx());
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (const auto bond : mol.bonds()) {
    // RDKit❗❌:     if (isPatternComplexQuery(bond)) {
    // RDKit❗❌:       isQueryBond.set(bond->getIdx());
    // RDKit❗❌:       if (tautomericFingerprint && isTautomerBondQuery(bond)) {
    // RDKit❗❌:         isTautomerBond.set(bond->getIdx());
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   unsigned int pIdx = 0;
    // RDKit❗❌:   for (const auto patt : patts) {
    // RDKit❗❌:     ++pIdx;
    // RDKit❗❌:     std::vector<MatchVectType> matches;
    // RDKit❗❌:     // uniquify matches?
    // RDKit❗❌:     //   time for 10K molecules w/ uniquify: 5.24s
    // RDKit❗❌:     //   time for 10K molecules w/o uniquify: 4.87s
    // RDKit❗❌:
    // RDKit❗❌:     SubstructMatchParameters params;
    // RDKit❗❌:     params.uniquify = false;
    // RDKit❗❌:     // raise maxMatches really high. This was the cause for github #2614.
    // RDKit❗❌:     // if we end up with more matches than this, we're completely hosed: :-)
    // RDKit❗❌:     params.maxMatches = 100000000;
    // RDKit❗❌:     matches = SubstructMatch(mol, *patt, params);
    // RDKit❗❌:
    // RDKit❗❌:     std::uint32_t mIdx = pIdx + patt->getNumAtoms() + patt->getNumBonds();
    // RDKit❗❌:     for (const auto &mv : matches) {
    // RDKit❗❌: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗❌:       std::cerr << "\nPatt: " << pIdx << " | ";
    // RDKit❗❌: #endif
    // RDKit❗❌:       // collect bits counting the number of occurrences of the pattern:
    // RDKit❗❌:       gboost::hash_combine(mIdx, 0xBEEF);
    // RDKit❗❌:       fp.setBit(mIdx % fpSize);
    // RDKit❗❌: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗❌:       std::cerr << "count: " << mIdx % fpSize << " | ";
    // RDKit❗❌: #endif
    // RDKit❗❌:
    // RDKit❗❌:       bool isQuery = false;
    // RDKit❗❌:       std::uint32_t bitId = pIdx;
    // RDKit❗❌:       std::vector<unsigned int> amap(mv.size(), 0);
    // RDKit❗❌:       for (const auto &p : mv) {
    // RDKit❗❌: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗❌:         std::cerr << p.second << " ";
    // RDKit❗❌: #endif
    // RDKit❗❌:         if (isQueryAtom[p.second]) {
    // RDKit❗❌:           isQuery = true;
    // RDKit❗❌: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗❌:           std::cerr << "atom query.";
    // RDKit❗❌: #endif
    // RDKit❗❌:           break;
    // RDKit❗❌:         }
    // RDKit❗❌:         gboost::hash_combine(bitId,
    // RDKit❗❌:                              mol.getAtomWithIdx(p.second)->getAtomicNum());
    // RDKit❗❌:         amap[p.first] = p.second;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (isQuery) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:       auto tautomerBitId = bitId;
    // RDKit❗❌:       auto tautomerQuery = false;
    // RDKit❗❌:       ROMol::EDGE_ITER firstB, lastB;
    // RDKit❗❌:       boost::tie(firstB, lastB) = patt->getEdges();
    // RDKit❗❌: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗❌:       std::cerr << " bs:|| ";
    // RDKit❗❌: #endif
    // RDKit❗❌:       while (!isQuery && firstB != lastB) {
    // RDKit❗❌:         const Bond *pbond = (*patt)[*firstB];
    // RDKit❗❌:         ++firstB;
    // RDKit❗❌:         const Bond *mbond = mol.getBondBetweenAtoms(
    // RDKit❗❌:             amap[pbond->getBeginAtomIdx()], amap[pbond->getEndAtomIdx()]);
    // RDKit❗❌:         const auto bondIdx = mbond->getIdx();
    // RDKit❗❌:
    // RDKit❗❌:         if (isQueryBond[bondIdx]) {
    // RDKit❗❌:           isQuery = true;
    // RDKit❗❌:           if (isTautomerBond[bondIdx]) {
    // RDKit❗❌:             isQuery = false;
    // RDKit❗❌:             tautomerQuery = true;
    // RDKit❗❌: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗❌:             std::cerr << "tautomer query: " << mbond->getIdx();
    // RDKit❗❌: #endif
    // RDKit❗❌:           }
    // RDKit❗❌:           if (isQuery) {
    // RDKit❗❌: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗❌:             std::cerr << "bond query: " << mbond->getIdx();
    // RDKit❗❌: #endif
    // RDKit❗❌:             break;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:         if (tautomericFingerprint) {
    // RDKit❗❌:           if (isTautomerBond[bondIdx] || mbond->getIsAromatic() ||
    // RDKit❗❌:               mbond->getBondType() == Bond::SINGLE ||
    // RDKit❗❌:               mbond->getBondType() == Bond::DOUBLE ||
    // RDKit❗❌:               mbond->getBondType() == Bond::AROMATIC) {
    // RDKit❗❌:             gboost::hash_combine(tautomerBitId, -1);
    // RDKit❗❌: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗❌:             std::cerr << "T ";
    // RDKit❗❌: #endif
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:         if (!tautomerQuery) {
    // RDKit❗❌:           if (!mbond->getIsAromatic()) {
    // RDKit❗❌:             gboost::hash_combine(bitId, (std::uint32_t)mbond->getBondType());
    // RDKit❗❌: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗❌:             std::cerr << mbond->getBondType() << " ";
    // RDKit❗❌: #endif
    // RDKit❗❌:           } else {
    // RDKit❗❌:             gboost::hash_combine(bitId, (std::uint32_t)Bond::AROMATIC);
    // RDKit❗❌: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗❌:             std::cerr << Bond::AROMATIC << " ";
    // RDKit❗❌: #endif
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       if (!isQuery) {
    // RDKit❗❌:         if (!tautomerQuery) {
    // RDKit❗❌: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗❌:           std::cerr << " set: " << bitId << " " << bitId % fpSize;
    // RDKit❗❌: #endif
    // RDKit❗❌:           fp.setBit(bitId % fpSize);
    // RDKit❗❌:         }
    // RDKit❗❌:         if (tautomericFingerprint) {
    // RDKit❗❌: #ifdef VERBOSE_FINGERPRINTING
    // RDKit❗❌:           std::cerr << " tset: " << tautomerBitId << " "
    // RDKit❗❌:                     << tautomerBitId % fpSize;
    // RDKit❗❌: #endif
    // RDKit❗❌:           fp.setBit(tautomerBitId % fpSize);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // RDKit❗❌:
    // RDKit❗❌: }  // namespace
    // RDKit❗❌:
    // RDKit❗❌: // caller owns the result, it must be deleted
    // Complexity: SEARCH owns the matcher and all query classification. One
    // borrowed ring-only context is shared across the 13 compiled queries.
    // Vec<bool> masks use byte storage rather than source bitset packing;
    // structural validation and query-carrier projection add O(A+B) work.
    // Query inputs clone common carrier properties once, not per match. The
    // shared SEARCH matcher allocation costs are not claimed source-equivalent.
    if fingerprint_size == 0 {
        return Err(PatternFingerprintError::EmptyFingerprint);
    }
    if fingerprint.n_bits() as usize != fingerprint_size {
        return Err(PatternFingerprintError::BitLengthMismatch {
            left: fingerprint.n_bits() as usize,
            right: fingerprint_size,
        });
    }
    if atom_counts
        .as_ref()
        .is_some_and(|counts| counts.len() < molecule.num_atoms())
    {
        return Err(PatternFingerprintError::InvalidArguments {
            reason: "Pattern atom_counts length is smaller than molecule atom count",
        });
    }
    if set_only_bits.is_some_and(|bits| bits.n_bits() as usize != fingerprint_size) {
        return Err(PatternFingerprintError::InvalidArguments {
            reason: "Pattern set_only_bits length differs from fingerprint size",
        });
    }

    let patterns = compiled_pattern_fingerprint_queries()?;
    let topology = molecule.matching_topology()?;
    let owned_rings;
    let rings = match rings.filter(|r| r.is_find_fast_or_better()) {
        Some(r) => r,
        None => {
            owned_rings = cosmolkit_core::fast_find_rings(&topology)?;
            &owned_rings
        }
    };
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(&topology, &coordinates, &[], Some(rings), None);
    let query_context = build_ring_only_query_match_context(&topology, rings)?;
    let is_query_atom: Vec<bool> = (0..molecule.num_atoms())
        .map(|i| molecule.is_query_atom(i))
        .collect();
    let mut is_query_bond = vec![false; molecule.num_bonds()];
    let mut is_tautomer_bond = vec![false; molecule.num_bonds()];
    for i in 0..molecule.num_bonds() {
        let root = molecule.bond_query(i);
        if is_pattern_complex_query(root) {
            is_query_bond[i] = true;
            is_tautomer_bond[i] = tautomeric_fingerprint && is_tautomer_bond_query(root);
        }
    }

    let mut params = SubstructMatchParams::default();
    params.uniquify = false;
    params.max_matches = 100_000_000;
    for (pattern_offset, matcher) in patterns.iter().enumerate() {
        let pattern_index = pattern_offset as u32 + 1;
        let pattern = matcher.query();
        let matches = try_get_substruct_atom_matches_with_compiled_query_and_context(
            &target,
            matcher,
            &params,
            &query_context,
        )?;
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
            fingerprint.set_bit(count_bit as u32)?;
            trace(PatternTraceEvent::CountBit {
                pattern_index,
                occurrence,
                seed: match_index,
                bit: count_bit,
            });

            let mut is_query = false;
            let mut bit_id = pattern_index;
            let mut atom_map = vec![0usize; matched.len()];
            for (query_atom_index, &molecule_atom_index) in matched.iter().enumerate() {
                if is_query_atom[molecule_atom_index] {
                    is_query = true;
                    trace(PatternTraceEvent::QueryAtomSuppressed {
                        pattern_index,
                        atom_index: molecule_atom_index,
                    });
                    break;
                }
                let atomic_number = u32::from(molecule.atomic_number(molecule_atom_index));
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
                let molecule_bond_index = topology
                    .adjacency
                    .neighbors_of(atom_map[pattern_bond.begin().index()])
                    .iter()
                    .find(|n| n.atom_index == atom_map[pattern_bond.end().index()])
                    .ok_or(PatternFingerprintError::InvalidArguments {
                        reason: "Pattern substructure mapping does not preserve query bond",
                    })?
                    .bond
                    .index();
                let molecule_bond = &topology.bonds[molecule_bond_index];

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
                        BondOrder::Aromatic.rdkit_code() as u32
                    } else {
                        molecule_bond.order().rdkit_code() as u32
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
                    fingerprint.set_bit(bit as u32)?;
                    trace(PatternTraceEvent::StructureBit {
                        pattern_index,
                        seed: bit_id,
                        bit,
                    });
                }
                if tautomeric_fingerprint {
                    let bit = tautomer_bit_id as usize % fingerprint_size;
                    fingerprint.set_bit(bit as u32)?;
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
    use cosmolkit_model::{AtomId, BondId, BondSpec};

    enum Fixture {
        Concrete(TopologyBlock),
        Query(QueryGraph),
    }
    impl Fixture {
        fn from_smiles(s: &str) -> Result<Self, cosmolkit_smiles::SmilesParseError> {
            Ok(Self::Concrete(
                cosmolkit_smiles::parse_smiles(s, &Default::default())?.topology,
            ))
        }
        fn input(&self) -> PatternGraphInput<'_> {
            match self {
                Self::Concrete(t) => PatternGraphInput::Concrete(t),
                Self::Query(q) => PatternGraphInput::Query(q),
            }
        }
        fn num_atoms(&self) -> usize {
            self.input().num_atoms()
        }
    }
    fn fixture_from_smarts(
        s: &str,
        p: &cosmolkit_search::SmartsParseParams,
    ) -> Result<Fixture, SmartsParseError> {
        Ok(Fixture::Query(parse_smarts(s, p)?))
    }
    fn update_pattern_fingerprint(
        m: &Fixture,
        f: &mut Fingerprint,
        n: usize,
        a: Option<&mut [u32]>,
        b: Option<&Fingerprint>,
        t: bool,
    ) -> Result<(), PatternFingerprintError> {
        super::update_pattern_fingerprint_impl(m.input(), None, f, n, a, b, t, |_| {})
    }
    fn is_pattern_complex_query(b: &Bond) -> bool {
        cosmolkit_search::is_pattern_complex_query(b.query())
    }
    fn is_tautomer_bond_query(b: &Bond) -> bool {
        cosmolkit_search::is_tautomer_bond_query(b.query())
    }

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
            let reparsed =
                fixture_from_smarts(expected, &cosmolkit_search::SmartsParseParams::default())
                    .expect("reparse source Pattern SMARTS");
            assert_eq!(matcher.query().num_atoms(), reparsed.input().num_atoms());
            assert_eq!(matcher.query().num_bonds(), reparsed.input().num_bonds());
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
        molecule: &Fixture,
        fingerprint_size: usize,
        tautomeric: bool,
    ) -> (Fingerprint, Vec<PatternTraceEvent>) {
        let mut fingerprint = Fingerprint::new(fingerprint_size as u32);
        let mut events = Vec::new();
        update_pattern_fingerprint_impl(
            molecule.input(),
            None,
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
        let molecule = Fixture::from_smiles("CC").expect("ethane");
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
            let molecule = Fixture::from_smiles(smiles).expect("Pattern fixture");
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
        let benzene = Fixture::from_smiles("c1ccccc1").expect("benzene");
        let (_, events) = traced_pattern_fingerprint(&benzene, 2048, false);
        assert!(
            events
                .iter()
                .any(|event| matches!(event, PatternTraceEvent::BondHash { bond_code: 12, .. }))
        );

        let sodium = Fixture::from_smiles("[Na+]").expect("sodium");
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
        let query_atom =
            fixture_from_smarts("[*]", &cosmolkit_search::SmartsParseParams::default())
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

        let query_bond =
            fixture_from_smarts("C~C", &cosmolkit_search::SmartsParseParams::default())
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

        let any_query = fixture_from_smarts("C~N", &cosmolkit_search::SmartsParseParams::default())
            .expect("any query bond");
        assert_eq!(
            traced_pattern_fingerprint(&any_query, 257, false)
                .0
                .on_bits(),
            vec![14, 43, 44, 132, 136, 146]
        );

        let order_query =
            fixture_from_smarts("C-,=N", &cosmolkit_search::SmartsParseParams::default())
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
        let query = fixture_from_smarts("CC", &cosmolkit_search::SmartsParseParams::default())
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
        let molecule = Fixture::from_smiles("CCC").expect("propane");
        let (width_one, _) = traced_pattern_fingerprint(&molecule, 1, true);
        assert_eq!(width_one.on_bits(), vec![0]);

        let mut baseline = Fingerprint::new(127);
        update_pattern_fingerprint(&molecule, &mut baseline, 127, None, None, true)
            .expect("baseline");
        let mut atom_counts = vec![17; molecule.num_atoms()];
        let set_only_bits = Fingerprint::from_on_bits(127, [0, 5, 126]).unwrap();
        let mut with_inert_arguments = Fingerprint::new(127);
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

        let mut invalid = Fingerprint::new(127);
        assert!(matches!(
            update_pattern_fingerprint(
                &molecule,
                &mut invalid,
                127,
                Some(&mut [0; 2]),
                None,
                false,
            ),
            Err(PatternFingerprintError::InvalidArguments { .. })
        ));
        assert!(matches!(
            update_pattern_fingerprint(
                &molecule,
                &mut invalid,
                127,
                None,
                Some(&Fingerprint::new(128)),
                false,
            ),
            Err(PatternFingerprintError::InvalidArguments { .. })
        ));
        let mut empty = Fingerprint::new(0);
        assert_eq!(
            update_pattern_fingerprint(&molecule, &mut empty, 0, None, None, false),
            Err(PatternFingerprintError::EmptyFingerprint)
        );
    }
}

#[cfg(test)]
#[path = "pattern_screening_tests.rs"]
mod screening_tests;
