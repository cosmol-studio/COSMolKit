//! Direct Lipinski descriptor owners (RDKit `Lipinski.cpp`).

use crate::{DescriptorError, DescriptorInput, DescriptorResult};
use cosmolkit_core::{ValenceAssignment, total_hydrogen_count_from_validated};
use cosmolkit_model::TopologyBlock;

/// Shared kernel for the direct Lipinski HBA count.
///
/// Behavior review: reproduces the source's single atom loop and its exact
/// `getAtomicNum() == 7 || getAtomicNum() == 8` predicate; charge, aromatic
/// flags, hydrogens, and degrees play no role, and the function has no error
/// path (the source returns an unsigned count). The preceding
/// `validate_topology` is a COSMolKit representational-safety precheck the
/// source performs implicitly by operating on a valid `ROMol`.
/// Complexity review: one O(n) pass over atom rows with an iterator
/// filter+count and zero allocation, matching the source's ConstAtomIterator
/// loop; no per-atom lookup or buffering is introduced.
pub(crate) fn lipinski_hba_kernel(topology: &TopologyBlock) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcLipinskiHBA
    // RDKit✔️✔️: unsigned int calcLipinskiHBA(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (ROMol::ConstAtomIterator iter = mol.beginAtoms(); iter != mol.endAtoms();
    // RDKit✔️✔️:        ++iter) {
    // RDKit✔️✔️:     if ((*iter)->getAtomicNum() == 7 || (*iter)->getAtomicNum() == 8) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcLipinskiHBA
    crate::validate_topology(topology, "lipinski_hba")?;
    u32::try_from(
        topology
            .atoms
            .iter()
            .filter(|atom| matches!(atom.atomic_number(), 7 | 8))
            .count(),
    )
    .map_err(|_| DescriptorError::Unsupported {
        function: "lipinski_hba",
        detail: "acceptor count exceeds the RDKit unsigned result model".to_owned(),
    })
}

/// Direct Lipinski hydrogen-bond acceptor count over prepared FINAL input.
///
/// This is the DIRECT N/O atom count (`Lipinski.cpp:77-86`), not the general
/// SMARTS-based `NumHBA` (Q03): quaternary `[N+]` still counts here even
/// though the recursive acceptor pattern excludes it.
pub fn lipinski_hba_prepared(input: &DescriptorInput) -> DescriptorResult<u32> {
    lipinski_hba_kernel(input.topology())
}

/// Shared kernel for the direct Lipinski HBD hydrogen sum.
///
/// Behavior review: reproduces the source loop — the N/O predicate then
/// `getTotalNumHs(true)`, i.e. explicit H count + implicit H rows + hydrogen
/// ATOM neighbors (the neighbor predicate is `getAtomicNum() == 1`, so
/// isotopic H neighbors such as `[2H]` count too). The borrowed assignment
/// is shape-checked, never recomputed; hydrogen rows flow through the
/// canonical core owner preserving typed errors (negative implicit rows stay
/// errors, never zero).
/// Complexity review: one O(n) atom pass with O(degree) neighbor scans only
/// on N/O rows through the pinned core owner; zero allocation beyond the
/// borrow. The source reads cached atom state; this kernel borrows the
/// caller's prepared assignment instead of reassigning valence.
pub(crate) fn lipinski_hbd_kernel(
    topology: &TopologyBlock,
    assignment: &ValenceAssignment,
) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcLipinskiHBD
    // RDKit✔️✔️: unsigned int calcLipinskiHBD(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (ROMol::ConstAtomIterator iter = mol.beginAtoms(); iter != mol.endAtoms();
    // RDKit✔️✔️:        ++iter) {
    // RDKit✔️✔️:     if (((*iter)->getAtomicNum() == 7 || (*iter)->getAtomicNum() == 8)) {
    // RDKit✔️✔️:       res += (*iter)->getTotalNumHs(true);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcLipinskiHBD
    crate::validate_topology(topology, "lipinski_hbd")?;
    let assignment = crate::prepared_valence(topology, Some(assignment), "lipinski_hbd")?;
    let mut result = 0_u32;
    for atom in &topology.atoms {
        if matches!(atom.atomic_number(), 7 | 8) {
            let hydrogens =
                total_hydrogen_count_from_validated(topology, assignment.as_ref(), atom.id(), true)
                    .map_err(|source| DescriptorError::Valence {
                        function: "lipinski_hbd",
                        source,
                    })?;
            result = result
                .checked_add(hydrogens)
                .ok_or_else(|| DescriptorError::Unsupported {
                    function: "lipinski_hbd",
                    detail: "donor hydrogen count exceeds the RDKit unsigned result model"
                        .to_owned(),
                })?;
        }
    }
    Ok(result)
}

/// Direct Lipinski hydrogen-bond donor hydrogen sum over prepared FINAL
/// input.
///
/// Borrows the prepared valence rows; never reassigns valence.
pub fn lipinski_hbd_prepared(input: &DescriptorInput) -> DescriptorResult<u32> {
    lipinski_hbd_kernel(input.topology(), input.valence())
}

/// Fraction of carbon rows whose SOURCE total degree (explicit bond
/// neighbors + implicit/explicit-property H, `includeNeighbors=false`) is
/// exactly four.
///
/// Behavior owner for `RDKit::Descriptors::calcFractionCSP3`
/// (Lipinski.cpp:212-232) including the source-defined zero-carbon fallback.
/// The criterion is total degree, never a hybridization flag. The borrowed
/// assignment is shape-checked, never recomputed.
///
/// Complexity review: one O(n) row pass; per carbon an O(1) neighbor-length
/// read plus O(1) hydrogen-row lookup through the pinned core owner; zero
/// allocation beyond the result f64.
pub(crate) fn fraction_csp3_kernel(
    topology: &TopologyBlock,
    assignment: &ValenceAssignment,
) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcFractionCSP3
    // RDKit✔️✔️: double calcFractionCSP3(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int nCSP3 = 0;
    // RDKit✔️✔️:   unsigned int nC = 0;
    // RDKit✔️✔️:   ROMol::VERTEX_ITER atBegin, atEnd;
    // RDKit✔️✔️:   boost::tie(atBegin, atEnd) = mol.getVertices();
    // RDKit✔️✔️:   while (atBegin != atEnd) {
    // RDKit✔️✔️:     const Atom *at = mol[*atBegin];
    // RDKit✔️✔️:     if (at->getAtomicNum() == 6) {
    // RDKit✔️✔️:       ++nC;
    // RDKit✔️✔️:       if (at->getTotalDegree() == 4) {
    // RDKit✔️✔️:         ++nCSP3;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++atBegin;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!nC) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return static_cast<double>(nCSP3) / nC;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcFractionCSP3
    // BEGIN RDKIT CPP FUNCTION: RDKit::Atom::getTotalDegree
    // RDKit✔️✔️: unsigned int Atom::getTotalDegree() const {
    // RDKit✔️✔️:   unsigned int res = this->getTotalNumHs(false) + this->getDegree();
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Atom::getTotalDegree
    crate::validate_topology(topology, "fraction_csp3")?;
    let assignment = crate::prepared_valence(topology, Some(assignment), "fraction_csp3")?;
    let mut carbon_count = 0_u32;
    let mut sp3_count = 0_u32;
    for atom in &topology.atoms {
        if atom.atomic_number() != 6 {
            continue;
        }
        carbon_count += 1;
        let hydrogens =
            total_hydrogen_count_from_validated(topology, assignment.as_ref(), atom.id(), false)
                .map_err(|source| DescriptorError::Valence {
                    function: "fraction_csp3",
                    source,
                })?;
        let degree = u32::try_from(topology.adjacency.neighbors_of(atom.id().index()).len())
            .map_err(|_| DescriptorError::Unsupported {
                function: "fraction_csp3",
                detail: "atom degree exceeds the RDKit unsigned result model".to_owned(),
            })?;
        if degree.saturating_add(hydrogens) == 4 {
            sp3_count += 1;
        }
    }
    if carbon_count == 0 {
        return Ok(0.0);
    }
    Ok(f64::from(sp3_count) / f64::from(carbon_count))
}

/// Prepared form of [`crate::fraction_csp3`].
///
/// Borrows the prepared valence rows; never reassigns valence.
pub fn fraction_csp3_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<f64> {
    fraction_csp3_kernel(input.topology(), input.valence())
}

/// Fixed SMARTS for the general hydrogen-bond donor count
/// (`Lipinski.cpp:197`). Single-atom query, comma-OR of five alternatives;
/// non-recursive (no `$`), so the source takes the shared-query branch.
const NUM_HBD_PATTERN: &str = "[N&!H0&v3,N&!H0&+1&v4,O&H1&+0,S&H1&+0,n&H1&+0]";

/// Version string the source exports next to `calcNumHBD`.
pub const NUM_HBD_VERSION: &str = "2.0.1";

/// General SMARTS-based hydrogen-bond donor count over prepared FINAL
/// input.
///
/// Behavior review: reproduces the `SMARTSCOUNTFUNC(NumHBD, ...)` expansion
/// — one retained pattern (`Q01`), default `SubstructMatchParameters`, and
/// `matches.size()`. The bracket is ONE single-atom query whose inner commas
/// are SMARTS OR: (aliphatic N, total H>0, v3) | (aliphatic N, H>0, +1, v4)
/// | (aliphatic O, H1, +0) | (aliphatic S, H1, +0) | (aromatic n, H1, +0);
/// `H` folds implicit and explicit-neighbor hydrogens through the prepared
/// valence rows, `v` is total valence, `+` formal charge. The count is the
/// number of matching atoms. This is NOT the direct Lipinski donor sum
/// ([`lipinski_hbd_prepared`]): sulfur donors match here, `[NH4+]` counts
/// one atom here but four attached hydrogens there.
///
/// Complexity review: one retained-query acquisition (map hit after the
/// first process call) plus one matcher pass; identical cost class to the
/// source flyweight + `SubstructMatch`, with the per-call parse the
/// historical port performed removed by construction.
pub(crate) fn num_hbd_with_valence(
    topology: &cosmolkit_model::TopologyBlock,
    valence: &cosmolkit_core::ValenceAssignment,
) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP MACRO EXPANSION: SMARTSCOUNTFUNC -> calcNumHBD
    // RDKit✔️✔️: #define SMARTSCOUNTFUNC(nm, pattern, vers)         \
    // RDKit✔️✔️:   const std::string nm##Version = vers;            \
    // RDKit✔️✔️:   unsigned int calc##nm(const RDKit::ROMol &mol) { \
    // RDKit✔️✔️:     pattern_flyweight m(pattern);                  \
    // RDKit✔️✔️:     return m.get().countMatches(mol);              \
    // RDKit✔️✔️:   }                                                \
    // RDKit✔️✔️:   extern int no_such_variable
    // RDKit✔️✔️: SMARTSCOUNTFUNC(NumHBD, "[N&!H0&v3,N&!H0&+1&v4,O&H1&+0,S&H1&+0,n&H1&+0]",
    // RDKit✔️✔️:                 "2.0.1");
    // END RDKIT CPP MACRO EXPANSION: SMARTSCOUNTFUNC -> calcNumHBD
    //
    // HBD-PUBLIC narrow owner: the fixed pattern reads hydrogen-count and
    // valence predicates but NO ring predicate (ROOT-resolved from the
    // pinned source), so this owner runs through the ONE narrow
    // borrowed-valence context and matcher path. The supplied prepared
    // valence rows are BORROWED — never recomputed, extended, installed or
    // cloned; no ring state is read, gated on or fabricated. Borrow review:
    // borrowing the rows does not erase the real validation/context/
    // adjacency/matcher cost; no allocation-free/O(1) claim is made here.
    crate::patterns::count_pattern_matches_with_valence(
        topology,
        valence,
        "num_hbd",
        NUM_HBD_PATTERN,
    )
}

/// General SMARTS-based hydrogen-bond donor count over prepared FINAL
/// input (thin routing form).
///
/// Routes through the ONE narrow HBD owner
/// [`num_hbd_with_valence`] using the prepared input's supplied FINAL
/// topology and valence rows; the pattern reads no ring predicate, so the
/// prepared ring rows are neither required nor consulted on this path.
/// Signature and prepared-input contract are unchanged.
pub fn num_hbd_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_hbd_with_valence(input.topology(), input.valence())
}

/// Fixed SMARTS for the general hydrogen-bond acceptor count
/// (`Lipinski.cpp:199-202`). Single-atom query whose predicate is the
/// comma-OR of five RECURSIVE `$()` subqueries, so the source takes the
/// `m_needCopies` (per-call query copy) branch; see the 🔝 note below.
const NUM_HBA_PATTERN: &str = "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0X2,o,s;+0])]";

/// Version string the source exports next to `calcNumHBA`.
pub const NUM_HBA_VERSION: &str = "2.0.2";

/// General SMARTS-based hydrogen-bond acceptor count over prepared FINAL
/// input.
///
/// Behavior review: reproduces the `SMARTSCOUNTFUNC(NumHBA, ...)`
/// expansion — one retained pattern (`Q01`), default
/// `SubstructMatchParameters`, `matches.size()`. The bracket is ONE
/// single-atom query; an atom matches if any recursive branch holds:
/// (1) O/S with exactly one total H and total valence 2 bonded to a
/// neighbor that is NOT `=[O,N,P,S]` (hydroxyl/thiol O–H / S–H on plain
/// carbon; excludes acid/ester/amide-adjacent OH), (2) O/S with zero H
/// and total valence 2 (ethers AND carbonyl oxygens — a `C=O` oxygen has
/// order-2 valence, confirmed by RDKit test.cpp:2227-2234), (3) anionic
/// O/S (carboxylate/sulfonate), (4) trivalent N whose single-or-aromatic
/// bond does not reach an `=!@`-bonded O/N/P/S (excludes amide-like N;
/// nitrile N still matches), (5) neutral aromatic `n` (no H, ring
/// connectivity 2), `o`, `s`. This is NOT the direct Lipinski N+O sum
/// ([`lipinski_hba_prepared`]): thiophene S counts here, acid OH does
/// not, and a carbonyl O counts here only because of branch 2.
///
/// Complexity review: one retained-query acquisition (map hit after the
/// first process call) plus one matcher pass; per-atom evaluation of the
/// five recursive subqueries is the same cost class as the source
/// flyweight + `SubstructMatch`, with the per-call parse the historical
/// port performed removed by construction.
///
/// 🔝 note: the source's `m_needCopies` deep-copies the query on EVERY
/// call for recursive-query thread safety (Lipinski.cpp:43-47). Rust
/// shares one immutable `Arc<QueryGraph>` instead; matching never mutates
/// the query, so semantics are unchanged while the per-call O(query)
/// copy is removed (same recorded deviation as the `Q01` rotatable
/// pattern).
pub fn num_hba_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP MACRO EXPANSION: SMARTSCOUNTFUNC -> calcNumHBA
    // RDKit✔️✔️: #define SMARTSCOUNTFUNC(nm, pattern, vers)         \
    // RDKit✔️✔️:   const std::string nm##Version = vers;            \
    // RDKit✔️✔️:   unsigned int calc##nm(const RDKit::ROMol &mol) { \
    // RDKit✔️✔️:     pattern_flyweight m(pattern);                  \
    // RDKit✔️✔️:     return m.get().countMatches(mol);              \
    // RDKit✔️✔️:   }                                                \
    // RDKit✔️✔️:   extern int no_such_variable
    // RDKit✔️✔️: SMARTSCOUNTFUNC(NumHBA,
    // RDKit✔️✔️:                 "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$("
    // RDKit✔️✔️:                 "[N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0X2,o,s;+0])]",
    // RDKit✔️✔️:                 "2.0.2");
    // (C++ adjacent-string-literal concatenation resolves to the single
    // NUM_HBA_PATTERN string above.)
    // END RDKIT CPP MACRO EXPANSION: SMARTSCOUNTFUNC -> calcNumHBA
    crate::patterns::count_pattern_matches(input, "num_hba", NUM_HBA_PATTERN)
}

/// Fixed SMARTS for the general heteroatom count (`Lipinski.cpp:203`).
/// Single-atom query, `;`-AND of two negated atomic-number primitives;
/// non-recursive (no `$`), so the source takes the shared-query branch.
/// Public because the pattern itself is frozen source data asserted by
/// the regressions.
pub const NUM_HETEROATOMS_PATTERN: &str = "[!#6;!#1]";

/// Version string the source exports next to `calcNumHeteroatoms`.
pub const NUM_HETEROATOMS_VERSION: &str = "1.0.1";

/// General SMARTS-based heteroatom count over prepared FINAL input.
///
/// Behavior review: reproduces the
/// `SMARTSCOUNTFUNC(NumHeteroatoms, "[!#6;!#1]", "1.0.1")` expansion —
/// one retained pattern (`Q01`), default `SubstructMatchParameters`,
/// `matches.size()`. The bracket is ONE single-atom query: atomic number
/// NOT 6 (carbon) AND NOT 1 (hydrogen). Only graph rows are examined:
/// implicit hydrogens never count; explicit `[2H]` rows are `#1` and are
/// excluded; dummy atoms (atomic number 0) satisfy both negations and DO
/// count — a pure pattern consequence, not an element whitelist.
///
/// Complexity review: one retained-query acquisition (map hit after the
/// first process call) plus one matcher pass over a non-recursive
/// single-atom query; identical cost class to the source flyweight +
/// `SubstructMatch`, with the per-call parse the historical port
/// performed removed by construction.
pub fn num_heteroatoms_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP MACRO EXPANSION: SMARTSCOUNTFUNC -> calcNumHeteroatoms
    // RDKit✔️✔️: #define SMARTSCOUNTFUNC(nm, pattern, vers)         \
    // RDKit✔️✔️:   const std::string nm##Version = vers;            \
    // RDKit✔️✔️:   unsigned int calc##nm(const RDKit::ROMol &mol) { \
    // RDKit✔️✔️:     pattern_flyweight m(pattern);                  \
    // RDKit✔️✔️:     return m.get().countMatches(mol);              \
    // RDKit✔️✔️:   }                                                \
    // RDKit✔️✔️:   extern int no_such_variable
    // RDKit✔️✔️: SMARTSCOUNTFUNC(NumHeteroatoms, "[!#6;!#1]", "1.0.1");
    // END RDKIT CPP MACRO EXPANSION: SMARTSCOUNTFUNC -> calcNumHeteroatoms
    crate::patterns::count_pattern_matches(input, "num_heteroatoms", NUM_HETEROATOMS_PATTERN)
}

/// General SMARTS-based heteroatom count over a bare TOPOLOGY (no
/// chemistry rows).
///
/// D-A resolution (HETERO-PUBLIC): the `[!#6;!#1]` pattern reads ONLY
/// atomic-number rows, so this entry routes through the ONE common
/// retained-pattern match owner with the narrow topology-only query
/// context (validation + borrowed adjacency; ring_info None, valence
/// None). NO valence assignment, ring find, fabricated chemistry rows or
/// direct element counting; default `SubstructMatchParameters`
/// (uniquify=true, maxMatches=1000) exactly as the prepared path.
///
/// Complexity review: one O(V+E) validation + one retained-query
/// acquisition + one matcher pass — the same cost class as the source
/// flyweight + `SubstructMatch` on a mol with precomputed state.
pub fn num_heteroatoms_topology(topology: &TopologyBlock) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP MACRO EXPANSION: SMARTSCOUNTFUNC -> calcNumHeteroatoms
    // RDKit✔️✔️: #define SMARTSCOUNTFUNC(nm, pattern, vers)         \
    // RDKit✔️✔️:   const std::string nm##Version = vers;            \
    // RDKit✔️✔️:   unsigned int calc##nm(const RDKit::ROMol &mol) { \
    // RDKit✔️✔️:     pattern_flyweight m(pattern);                  \
    // RDKit✔️✔️:     return m.get().countMatches(mol);              \
    // RDKit✔️✔️:   }                                                \
    // RDKit✔️✔️:   extern int no_such_variable
    // RDKit✔️✔️: SMARTSCOUNTFUNC(NumHeteroatoms, "[!#6;!#1]", "1.0.1");
    // END RDKIT CPP MACRO EXPANSION: SMARTSCOUNTFUNC -> calcNumHeteroatoms
    crate::patterns::count_pattern_matches_topology_only(
        topology,
        "num_heteroatoms",
        NUM_HETEROATOMS_PATTERN,
    )
}

/// Fixed SMARTS for the general amide-bond count (`Lipinski.cpp:204`).
/// THREE-atom connected query; non-recursive (no `$`), so the source takes
/// the shared-query branch.
const NUM_AMIDE_BONDS_PATTERN: &str = "C(=[O;!R])N";

/// Version string the source exports next to `calcNumAmideBonds`.
pub const NUM_AMIDE_BONDS_VERSION: &str = "1.0.0";

/// General SMARTS-based amide-bond count over prepared FINAL input.
///
/// Behavior review: reproduces the
/// `SMARTSCOUNTFUNC(NumAmideBonds, "C(=[O;!R])N", "1.0.0")` expansion —
/// one retained pattern (`Q01`), default `SubstructMatchParameters`
/// (uniquify on), `matches.size()`. The query is aliphatic C with a
/// double-bonded O that is NOT in a ring (`!R` — the carbonyl O of a
/// lactam is acyclic, so cyclic amides still count) and a default
/// single-or-aromatic bond to an aliphatic N with no H/charge
/// constraint. The count is the number of UNIQUE full match vectors:
/// urea `O=C(N)N` yields two matches (the two N mappings are distinct
/// vectors); sulfonamides do not match because the carbonyl atom must
/// be carbon.
///
/// Complexity review: one retained-query acquisition (map hit after the
/// first process call) plus one matcher pass over a non-recursive
/// three-atom query; identical cost class to the source flyweight +
/// `SubstructMatch`, with the per-call parse the historical port
/// performed removed by construction.
pub fn num_amide_bonds_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP MACRO EXPANSION: SMARTSCOUNTFUNC -> calcNumAmideBonds
    // RDKit✔️✔️: #define SMARTSCOUNTFUNC(nm, pattern, vers)         \
    // RDKit✔️✔️:   const std::string nm##Version = vers;            \
    // RDKit✔️✔️:   unsigned int calc##nm(const RDKit::ROMol &mol) { \
    // RDKit✔️✔️:     pattern_flyweight m(pattern);                  \
    // RDKit✔️✔️:     return m.get().countMatches(mol);              \
    // RDKit✔️✔️:   }                                                \
    // RDKit✔️✔️:   extern int no_such_variable
    // RDKit✔️✔️: SMARTSCOUNTFUNC(NumAmideBonds, "C(=[O;!R])N", "1.0.0");
    // END RDKIT CPP MACRO EXPANSION: SMARTSCOUNTFUNC -> calcNumAmideBonds
    crate::patterns::count_pattern_matches(input, "num_amide_bonds", NUM_AMIDE_BONDS_PATTERN)
}
