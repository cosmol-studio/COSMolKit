//! Subgraph isomorphism matching (VF2) for molecule pattern matching.
//!
//! ## RDKit provenance (protocol: dev/source_reproduction_protocol.md)
//!
//! This module reproduces RDKit's substructure matching from:
//! - `third_party/rdkit/Code/GraphMol/Substruct/vf2.hpp` (~682 lines C++)
//! - `third_party/rdkit/Code/GraphMol/Substruct/SubstructMatch.cpp` (~735 lines C++)
//!
//! The VF2 algorithm implementation is adapted from vflib-2.0 by P. Foggia,
//! extensively modified by Greg Landrum, ported to Rust with depth-based
//! term_1/term_2 tracking (BackTrack decrements counters instead of
//! recomputing from scratch).
//!
//! ## Marker convention
//!
//! Each copied C++ block below uses the two-axis status marker:
//! - RDKit✔️✔️: fully reproduced behavior and performance
//! - RDKit✔️❌: functionally correct, but with a known performance gap
//! - RDKit❗✔️: unfinished behavior that must not be presented as parity
//! - RDKit❌❌: not yet ported

use crate::chemistry::forcefield::crystalff::build_crystalff_query_molecule;
use crate::search::query::{
    QueryMatchContext, atom_predicate_matches_with_context, bond_predicate_matches_with_context,
    build_query_match_context,
};
use crate::{
    Atom, AtomQueryPredicate, Bond, BondOrder, BondQueryPredicate, BondStereo, ChiralTag, Molecule,
    StereoGroupKind,
};
use std::collections::BTreeMap;
use std::env;
use std::time::Instant;

fn row61_torsion349_substruct_trace_enabled(
    mol: &Molecule,
    query: &Molecule,
    params: &SubstructMatchParams,
) -> bool {
    env::var("RDKIT_ROW61_TRACE").ok().as_deref() == Some("1")
        && mol.num_atoms() == 26
        && query.num_atoms() == 5
        && query.num_bonds() == 4
        && !params.uniquify
}

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

/// Result of a single substructure match.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SubstructMatchResult {
    /// Mapping from query atom index to molecule atom index.
    pub atom_mapping: Vec<usize>,
    /// Mapping from query bond index to molecule bond index.
    pub bond_mapping: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SubstructMatchError {
    #[error(
        "RDKit substructure matching branch {branch} is unsupported until {rdkit_function} is source-ported"
    )]
    Unsupported {
        branch: &'static str,
        rdkit_function: &'static str,
    },
}

type SubstructMatchResultList = Result<Vec<SubstructMatchResult>, SubstructMatchError>;

/// Parameters controlling substructure matching behaviour.
#[derive(Debug, Clone)]
pub struct SubstructMatchParams {
    /// Maximum number of matches to return (default: 1000).
    pub max_matches: usize,
    /// Whether to uniquify results (default: true).
    pub uniquify: bool,
    /// Whether atom/bond stereochemistry participates in matching.
    pub use_chirality: bool,
    /// Whether specified query stereo may match unspecified molecule stereo.
    pub specified_stereo_query_matches_unspecified: bool,
}

type RecursiveQueryMatchCache = BTreeMap<String, Vec<bool>>;

impl Default for SubstructMatchParams {
    fn default() -> Self {
        // RDKit✔️✔️: bool useChirality = false;  //!< Use chirality in determining whether or not
        // RDKit✔️✔️:                             //!< atoms/bonds match
        // RDKit✔️✔️: bool uniquify = true;            //!< uniquify (by atom index) match results
        // RDKit✔️✔️: unsigned int maxMatches = 1000;  //!< maximum number of matches to return
        // RDKit✔️✔️: bool specifiedStereoQueryMatchesUnspecified =
        // RDKit✔️✔️:     false;  //!< If set, query atoms and bonds with specified stereochemistry
        // RDKit✔️✔️:             //!< will match atoms and bonds with unspecified stereochemistry
        Self {
            max_matches: 1000,
            uniquify: true,
            use_chirality: false,
            specified_stereo_query_matches_unspecified: false,
        }
    }
}

// ---------------------------------------------------------------------------
// Internal minimum-degree graph representation
// ---------------------------------------------------------------------------

/// Minimal adjacency info needed for VF2.
///
/// `nbrs[i]` is a slice into `edges`.
#[derive(Debug, Clone)]
struct Vf2Graph {
    n_atoms: usize,
    n_bonds: usize,
    /// For each atom index, the neighbor indices and bond ids.
    adjacency: Vec<Vec<(usize, usize)>>, // (neighbor_atom_index, bond_index)
}

/// Build a VF2-compatible adjacency view from a molecule.
///
/// The C++ code iterates `out_edges` via Boost graph. We pre-build adjacency
/// once and use raw index lookups.
fn build_vf2_graph(mol: &Molecule) -> Vf2Graph {
    // RDKit source (implicit in vf2.hpp usage of out_edges):
    //   The VF2 state stores Graph *g1, *g2 and calls:
    //     boost::out_edges(node, *g)
    //     boost::out_degree(node, *g)
    //     boost::adjacent_vertices(node, *g)
    //   These are all O(1) in Boost adjacency_list.
    //
    // RDKit✔️❌: We build a flat adjacency Vec<(usize, usize)> per atom.
    //   This adds a one-time O(V+E) allocation vs the Boost inline storage,
    //   but lookups are O(degree) which matches the original hot-path cost.
    let n_atoms = mol.num_atoms();
    let mut adjacency: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n_atoms];
    for (bond_idx, bond) in mol.bonds().iter().enumerate() {
        let b = bond.begin().index();
        let e = bond.end().index();
        adjacency[b].push((e, bond_idx));
        adjacency[e].push((b, bond_idx));
    }
    Vf2Graph {
        n_atoms,
        n_bonds: mol.num_bonds(),
        adjacency,
    }
}

impl Vf2Graph {
    fn out_degree(&self, node: usize) -> usize {
        self.adjacency[node].len()
    }

    fn out_edges(&self, node: usize) -> &[(usize, usize)] {
        &self.adjacency[node]
    }
}

// ---------------------------------------------------------------------------
// Atom and bond matching functors
// ---------------------------------------------------------------------------

// RDKit source (SubstructMatch.cpp):
//   class AtomLabelFunctor {
//    public:
//     AtomLabelFunctor(const ROMol &query, const ROMol &mol,
//                      const SubstructMatchParameters &ps)
//         : d_query(query), d_mol(mol), d_params(ps) {};
//     bool operator()(unsigned int i, unsigned int j) const {
//       bool res = false;
//       if (d_params.useChirality) {
//         const Atom *qAt = d_query.getAtomWithIdx(i);
//         if (qAt->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
//             qAt->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
//           const Atom *mAt = d_mol.getAtomWithIdx(j);
//           if (!d_params.specifiedStereoQueryMatchesUnspecified &&
//               mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW &&
//               mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW) {
//             return false;
//           }
//         }
//       }
//       res = atomCompat(d_query[i], d_mol[j], d_params);
//       return res;
//     }
//    private:
//     const ROMol &d_query;
//     const ROMol &d_mol;
//     const SubstructMatchParameters &d_params;
//   };
//
// RDKit❗✔️: AtomLabelFunctor is ported as plain functions. The
//   useChirality specified/unspecified precheck is wired below; the final
//   tetrahedral parity check remains in MolMatchFinalCheckFunctor.

fn has_rdkit_tetrahedral_chiral_label(atom: &Atom) -> bool {
    matches!(
        atom.chiral_tag(),
        ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
    )
}

fn atom_label_chirality_matches(
    query_atom: &Atom,
    mol_atom: &Atom,
    params: &SubstructMatchParams,
) -> bool {
    if !params.use_chirality {
        return true;
    }
    // RDKit✔️✔️:     if (d_params.useChirality) {
    // RDKit✔️✔️:       const Atom *qAt = d_query.getAtomWithIdx(i);
    // RDKit✔️✔️:       if (qAt->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:           qAt->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:         const Atom *mAt = d_mol.getAtomWithIdx(j);
    // RDKit✔️✔️:         if (!d_params.specifiedStereoQueryMatchesUnspecified &&
    // RDKit✔️✔️:             mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW &&
    // RDKit✔️✔️:             mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    if has_rdkit_tetrahedral_chiral_label(query_atom)
        && !params.specified_stereo_query_matches_unspecified
        && !has_rdkit_tetrahedral_chiral_label(mol_atom)
    {
        return false;
    }
    true
}

/// Unfinished atom compatibility check for RDKit's atomCompat.
///
/// RDKit `atomCompat` checks atomic number, aromaticity, isotope,
/// formal charge, and query atom queries. This remains unfinished until the
/// full source behavior and parity surface are closed.
fn atom_matches(query_atom: &Atom, query_mol: &Molecule, mol_atom: &Atom, mol: &Molecule) -> bool {
    // RDKit source (atomCompat, unfinished):
    //   if (query->hasQuery()) {
    //     return query->Match(mol);
    //   }
    //   if (query->getAtomicNum() != mol->getAtomicNum() &&
    //       query->getAtomicNum() > 0) return false;
    //   if (query->getIsAromatic() != mol->getIsAromatic() &&
    //       query->getIsAromatic()) return false;
    //   if (query->getIsotope() && query->getIsotope() != mol->getIsotope())
    //     return false;
    //   if (query->getFormalCharge() != mol->getFormalCharge() &&
    //       query->getFormalCharge() != 0) return false;
    //   return true;

    // RDKit❗✔️: If the query atom has a query predicate tree, evaluate it
    // through the non-cached path. The recursive-SMARTS cached path is wired
    // separately through `atom_matches_with_recursive_cache()`.
    if let Some(query_node) = query_atom.query() {
        let query_ctx = build_query_match_context(mol);
        return match query_node {
            crate::QueryNode::Predicate(pred) => atom_query_predicate_matches_for_substruct(
                mol_atom,
                pred,
                mol,
                &SubstructMatchParams::default(),
                None,
                &query_ctx,
            ),
            crate::QueryNode::And(children) => children.iter().all(|child| match child {
                crate::QueryNode::Predicate(pred) => atom_query_predicate_matches_for_substruct(
                    mol_atom,
                    pred,
                    mol,
                    &SubstructMatchParams::default(),
                    None,
                    &query_ctx,
                ),
                _ => evaluate_atom_query(
                    child,
                    mol_atom,
                    mol,
                    &SubstructMatchParams::default(),
                    None,
                    &query_ctx,
                ),
            }),
            crate::QueryNode::Or(children) => children.iter().any(|child| match child {
                crate::QueryNode::Predicate(pred) => atom_query_predicate_matches_for_substruct(
                    mol_atom,
                    pred,
                    mol,
                    &SubstructMatchParams::default(),
                    None,
                    &query_ctx,
                ),
                _ => evaluate_atom_query(
                    child,
                    mol_atom,
                    mol,
                    &SubstructMatchParams::default(),
                    None,
                    &query_ctx,
                ),
            }),
            crate::QueryNode::Not(child) => !evaluate_atom_query(
                child,
                mol_atom,
                mol,
                &SubstructMatchParams::default(),
                None,
                &query_ctx,
            ),
        };
    }

    let q_an = query_atom.atomic_number();
    let m_an = mol_atom.atomic_number();

    // If query specifies a (non-dummy) atomic number, it must match.
    if q_an != 0 && q_an != m_an {
        return false;
    }

    // Aromaticity check: if query atom is aromatic, mol atom must also be aromatic.
    if query_atom.is_aromatic() && !mol_atom.is_aromatic() {
        return false;
    }

    // Isotope check: if query specifies an isotope, it must match.
    if let Some(q_iso) = query_atom.isotope() {
        if mol_atom.isotope() != Some(q_iso) {
            return false;
        }
    }

    // Formal charge check: if query has a non-zero charge, it must match.
    let q_charge = query_atom.formal_charge();
    if q_charge != 0 && q_charge != mol_atom.formal_charge() {
        return false;
    }

    let _ = query_mol;
    true
}

fn recursive_smarts_root_matches(
    atom: &Atom,
    recursive_smarts: &str,
    mol: &Molecule,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.h :: RecursiveStructureQuery
    // RDKit✔️✔️: class RDKIT_GRAPHMOL_EXPORT RecursiveStructureQuery
    // RDKit✔️✔️:     : public Queries::SetQuery<int, Atom const *, true> {
    // RDKit✔️✔️:   RecursiveStructureQuery(ROMol const *query, unsigned int serialNumber = 0)
    // RDKit✔️✔️:       : Queries::SetQuery<int, Atom const *, true>(),
    // RDKit✔️✔️:         d_serialNumber(serialNumber) {
    // RDKit✔️✔️:     setQueryMol(query);
    // RDKit✔️✔️:     setDataFunc(getAtIdx);
    // RDKit✔️✔️:     setDescription("RecursiveStructure");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   static inline int getAtIdx(Atom const *at) {
    // RDKit✔️✔️:     PRECONDITION(at, "bad atom argument");
    // RDKit✔️✔️:     return at->getIdx();
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    //
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructMatch.cpp :: detail::RecursiveMatcher
    // RDKit✔️✔️:   if (!query.hasProp(common_properties::_queryRootAtom)) {
    // RDKit✔️✔️:     matches.push_back(pairs.begin()->second);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     int rootIdx;
    // RDKit✔️✔️:     query.getProp(common_properties::_queryRootAtom, rootIdx);
    // RDKit✔️✔️:     bool found = false;
    // RDKit✔️✔️:     for (const auto &pairIter : pairs) {
    // RDKit✔️✔️:       if (pairIter.first == static_cast<unsigned int>(rootIdx)) {
    // RDKit✔️✔️:         matches.push_back(pairIter.second);
    // RDKit✔️✔️:         found = true;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    //
    // COSMolKit currently parses the recursive SMARTS used by Lipinski NumHBA
    // without `_queryRootAtom`; matching therefore uses RDKit's first mapped
    // query atom as the recursive root and tests membership in the cached
    // RecursiveStructureQuery atom-index set.
    if let Some(cache) = recursive_cache
        && let Some(match_starts) = cache.get(recursive_smarts)
    {
        return match_starts
            .get(atom.id().index())
            .copied()
            .unwrap_or(false);
    }

    let Some(inner) = recursive_smarts
        .strip_prefix("$(")
        .and_then(|value| value.strip_suffix(')'))
    else {
        return false;
    };
    let Ok(query) = build_crystalff_query_molecule(inner) else {
        return false;
    };
    substruct_match_impl(
        mol,
        &query,
        &SubstructMatchParams {
            max_matches: 1000,
            uniquify: false,
            use_chirality: false,
            specified_stereo_query_matches_unspecified: false,
        },
    )
    .unwrap_or_default()
    .into_iter()
    .any(|matched| matched.atom_mapping.first().copied() == Some(atom.id().index()))
}

fn atom_query_predicate_matches_for_substruct(
    atom: &Atom,
    pred: &AtomQueryPredicate,
    mol: &Molecule,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
    query_ctx: &QueryMatchContext,
) -> bool {
    match pred {
        // RDKit✔️✔️: Chiral SMARTS labels are not ordinary atom-compatibility
        // constraints when `useChirality` is false. AtomLabelFunctor and
        // MolMatchFinalCheckFunctor handle stereochemistry explicitly.
        AtomQueryPredicate::ChiralTagMatch(_) if !params.use_chirality => true,
        AtomQueryPredicate::RecursiveSmarts(recursive_smarts) => {
            recursive_smarts_root_matches(atom, recursive_smarts, mol, recursive_cache)
        }
        _ => atom_predicate_matches_with_context(atom, pred, mol, query_ctx),
    }
}

/// RDKit❗✔️: Evaluation of an atom query node for the SMARTS subset currently
/// modeled by COSMolKit.
///
/// Recursive SMARTS are evaluated through the recursive match cache used by
/// SubstructMatch; unsupported predicate leaves still evaluate false.
fn evaluate_atom_query(
    query: &crate::QueryNode<AtomQueryPredicate>,
    atom: &Atom,
    mol: &Molecule,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
    query_ctx: &QueryMatchContext,
) -> bool {
    match query {
        crate::QueryNode::Predicate(pred) => atom_query_predicate_matches_for_substruct(
            atom,
            pred,
            mol,
            params,
            recursive_cache,
            query_ctx,
        ),
        crate::QueryNode::And(children) => children
            .iter()
            .all(|c| evaluate_atom_query(c, atom, mol, params, recursive_cache, query_ctx)),
        crate::QueryNode::Or(children) => children
            .iter()
            .any(|c| evaluate_atom_query(c, atom, mol, params, recursive_cache, query_ctx)),
        crate::QueryNode::Not(child) => {
            !evaluate_atom_query(child, atom, mol, params, recursive_cache, query_ctx)
        }
    }
}

// RDKit source (SubstructMatch.cpp):
//   class BondLabelFunctor {
//    public:
//     BondLabelFunctor(const ROMol &query, const ROMol &mol,
//                      const SubstructMatchParameters &ps)
//         : d_query(query), d_mol(mol), d_params(ps) {};
//     bool operator()(MolGraph::edge_descriptor i,
//                     MolGraph::edge_descriptor j) const {
//       if (d_params.useChirality) {
//         const Bond *qBnd = d_query[i];
//         if (qBnd->getBondType() == Bond::DOUBLE &&
//             qBnd->getStereo() > Bond::STEREOANY) {
//           const Bond *mBnd = d_mol[j];
//           if (mBnd->getBondType() == Bond::DOUBLE &&
//               !d_params.specifiedStereoQueryMatchesUnspecified &&
//               mBnd->getStereo() <= Bond::STEREOANY) {
//             return false;
//           }
//         }
//       }
//       bool res = bondCompat(d_query[i], d_mol[j], d_params);
//       return res;
//     }
//    private:
//     const ROMol &d_query;
//     const ROMol &d_mol;
//     const SubstructMatchParameters &d_params;
//   };

/// RDKit❗✔️: bond compatibility check for RDKit's bondCompat.
///
/// Handles query-bond predicate trees, bond order matching, and aromatic
/// fallback rules required by the default Morgan feature SMARTS patterns.
fn bond_matches(query_bond: &Bond, query_mol: &Molecule, mol_bond: &Bond, mol: &Molecule) -> bool {
    // RDKit source (bondCompat, unfinished):
    //   if (qBnd->hasQuery()) return qBnd->Match(mBnd);
    //   if (qBnd->getIsAromatic() && mBnd->getIsAromatic()) return true;
    //   return qBnd->getBondType() == mBnd->getBondType();

    // RDKit❗✔️: If query bond has a query predicate tree, try to evaluate.
    if let Some(query_node) = query_bond.query() {
        let query_ctx = build_query_match_context(mol);
        return evaluate_bond_query(query_node, mol_bond, mol, &query_ctx);
    }

    let q_aromatic = query_bond.is_aromatic();
    let m_aromatic = mol_bond.is_aromatic();

    // Both aromatic: match regardless of bond order.
    if q_aromatic && m_aromatic {
        return true;
    }

    // Either could be aromatic, use order comparison.
    // RDKit rule: aromatic query bond matches aromatic mol bond (handled above)
    // or single bond in mol (for Kekule forms).
    // Also: query single bond matches aromatic mol bond.
    let q_order = query_bond.order();
    let m_order = mol_bond.order();

    if q_order == m_order {
        return true;
    }

    // Aromatic/Single interchange.
    if q_aromatic && m_order == BondOrder::Single {
        return true;
    }
    if m_aromatic && q_order == BondOrder::Single {
        return true;
    }

    let _ = query_mol;
    false
}

fn rdkit_bond_stereo_is_above_any(stereo: BondStereo) -> bool {
    !matches!(stereo, BondStereo::None | BondStereo::Any)
}

fn bond_label_chirality_matches(
    query_bond: &Bond,
    mol_bond: &Bond,
    params: &SubstructMatchParams,
) -> bool {
    if !params.use_chirality {
        return true;
    }
    // RDKit✔️✔️:     if (d_params.useChirality) {
    // RDKit✔️✔️:       const Bond *qBnd = d_query[i];
    // RDKit✔️✔️:       if (qBnd->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:           qBnd->getStereo() > Bond::STEREOANY) {
    // RDKit✔️✔️:         const Bond *mBnd = d_mol[j];
    // RDKit✔️✔️:         if (mBnd->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:             !d_params.specifiedStereoQueryMatchesUnspecified &&
    // RDKit✔️✔️:             mBnd->getStereo() <= Bond::STEREOANY) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    if query_bond.order() == BondOrder::Double
        && rdkit_bond_stereo_is_above_any(query_bond.stereo())
        && mol_bond.order() == BondOrder::Double
        && !params.specified_stereo_query_matches_unspecified
        && !rdkit_bond_stereo_is_above_any(mol_bond.stereo())
    {
        return false;
    }
    true
}

/// RDKit❗✔️: Evaluation of a bond query node for the currently modeled SMARTS
/// bond predicate subset.
fn evaluate_bond_query(
    query: &crate::QueryNode<BondQueryPredicate>,
    bond: &Bond,
    mol: &Molecule,
    query_ctx: &QueryMatchContext,
) -> bool {
    match query {
        crate::QueryNode::Predicate(pred) => {
            bond_predicate_matches_with_context(bond, pred, mol, query_ctx)
        }
        crate::QueryNode::And(children) => children
            .iter()
            .all(|c| evaluate_bond_query(c, bond, mol, query_ctx)),
        crate::QueryNode::Or(children) => children
            .iter()
            .any(|c| evaluate_bond_query(c, bond, mol, query_ctx)),
        crate::QueryNode::Not(child) => !evaluate_bond_query(child, bond, mol, query_ctx),
    }
}

// ---------------------------------------------------------------------------
// VF2 State Machine
// ---------------------------------------------------------------------------
//
// ## RDKit source reproduction: vf2.hpp
//
// The following section reproduces the VF2SubState class from vf2.hpp.
// The C++ code is shown as verbatim comments with RDKit markers.
//
// ### Key design differences from RDKit:
//
// 1. `core_1`/`core_2`: Same role — mapping from query atom idx → mol atom idx
//    and vice versa. Uses `Option<usize>` instead of NULL_NODE sentinel.
//
// 2. `term_1`/`term_2`: Stores the core_len *depth* at which each atom was
//    added to the terminal set, exactly as in vf2.hpp. BackTrack decrements
//    counters keyed by depth, not recomputes from scratch.
//
// 3. No shared_ptr copy semantics: VF2SubState in RDKit uses COW with
//    `share_count`. Rust's Clone+Vf2State avoids raw pointer sharing.
//    This means each VF2 recursive branch owns its state, which is
//    semantically correct but allocates O(depth * n) instead of
//    O(n) shared storage. For typical molecule sizes (<1000 atoms) this
//    is negligible; for very large searches the COW approach could be
//    reinstated with Arc<Vec<NodeId>>.
//
// 4. No boost graph: hand-rolled Vf2Graph adjacency.

// RDKit source (vf2.hpp):
//   typedef std::uint32_t node_id;
//   const node_id NULL_NODE = 0xFFFFFFFF;

type NodeId = usize;
const NULL_NODE: NodeId = usize::MAX;

// RDKit source (vf2.hpp):
//   template <class Graph>
//   struct Pair {
//     node_id n1, n2;
//     bool hasiter{false};
//     RDK_ADJ_ITER nbrbeg, nbrend;
//     Pair() : n1(NULL_NODE), n2(NULL_NODE) {}
//   };

#[derive(Debug, Clone)]
struct Vf2Pair {
    n1: NodeId,
    n2: NodeId,
    hasiter: bool,
    /// VF2+ source atom in the mol graph (g2) whose adjacency drives the
    /// neighbor iterator.
    nbr_node: NodeId,
    /// VF2+ neighbor iterator over mol graph (g2) neighbors.
    nbr_cursor: usize,
    nbr_end: usize,
}

impl Vf2Pair {
    fn new() -> Self {
        Self {
            n1: NULL_NODE,
            n2: NULL_NODE,
            hasiter: false,
            nbr_node: NULL_NODE,
            nbr_cursor: 0,
            nbr_end: 0,
        }
    }
}

// RDKit source (vf2.hpp):
//   /**
//    * The ordering by in/out degree
//    */
//   static bool nodeInfoComp1(const NodeInfo &a, const NodeInfo &b) {
//     if (a.out < b.out) { return true; }
//     if (a.out > b.out) { return false; }
//     if (a.in < b.in) { return true; }
//     if (a.in > b.in) { return false; }
//     return false;
//   }

#[derive(Debug, Clone, Copy)]
struct NodeInfo {
    id: usize,
    in_deg: usize,
    out_deg: usize,
}

// RDKit✔️✔️: nodeInfoComp1 — sort ascending by out-degree then in-degree.
fn node_info_cmp1(a: &NodeInfo, b: &NodeInfo) -> std::cmp::Ordering {
    // RDKit✔️✔️: if (a.out < b.out) { return true; }
    // RDKit✔️✔️: if (a.out > b.out) { return false; }
    // RDKit✔️✔️: if (a.in < b.in) { return true; }
    // RDKit✔️✔️: if (a.in > b.in) { return false; }
    // RDKit✔️✔️: return false;
    a.out_deg
        .cmp(&b.out_deg)
        .then_with(|| a.in_deg.cmp(&b.in_deg))
}

// RDKit source (vf2.hpp):
//   static int nodeInfoComp2(const NodeInfo &a, const NodeInfo &b) {
//     if (!a.in && b.in) return 1;
//     if (a.in && !b.in) return -1;
//     if (a.out < b.out) return -1;
//     if (a.out > b.out) return 1;
//     if (a.in < b.in) return -1;
//     if (a.in > b.in) return 1;
//     return 0;
//   }

// RDKit✔️✔️: nodeInfoComp2 — sort ascending by frequency (out=run count) then
//   valence (in=valence sum), except zero-valence nodes sort after non-zero.
fn node_info_cmp2(a: &NodeInfo, b: &NodeInfo) -> std::cmp::Ordering {
    // RDKit✔️✔️: if (!a.in && b.in) return 1;
    // RDKit✔️✔️: if (a.in && !b.in) return -1;
    if a.in_deg == 0 && b.in_deg != 0 {
        return std::cmp::Ordering::Greater;
    }
    if a.in_deg != 0 && b.in_deg == 0 {
        return std::cmp::Ordering::Less;
    }
    // RDKit✔️✔️: if (a.out < b.out) return -1;
    // RDKit✔️✔️: if (a.out > b.out) return 1;
    // RDKit✔️✔️: if (a.in < b.in) return -1;
    // RDKit✔️✔️: if (a.in > b.in) return 1;
    // RDKit✔️✔️: return 0;
    a.out_deg
        .cmp(&b.out_deg)
        .then_with(|| a.in_deg.cmp(&b.in_deg))
}

// RDKit source (vf2.hpp), SortNodesByFrequency:
//   Sorts the nodes of a graphs, returning a heap-allocated vector
//   with the node ids in the proper orders.
//   The sorting criterion takes into account:
//     1 - The number of nodes with the same in/out degree.
//     2 - The valence of the nodes.
//   The nodes at the beginning of the vector are the most singular,
//   from which the matching should start.

// RDKit✔️✔️: SortNodesByFrequency — returns query node ordering for VF2.
fn sort_nodes_by_frequency(g: &Vf2Graph) -> Vec<NodeId> {
    // RDKit✔️✔️: template <class Graph>
    // RDKit✔️✔️: node_id *SortNodesByFrequency(const Graph *g) {
    // RDKit✔️✔️:   std::vector<NodeInfo> vect;
    // RDKit✔️✔️:   vect.reserve(boost::num_vertices(*g));
    // RDKit✔️✔️:   typename Graph::vertex_iterator bNode, eNode;
    // RDKit✔️✔️:   boost::tie(bNode, eNode) = boost::vertices(*g);
    // RDKit✔️✔️:   while (bNode != eNode) {
    // RDKit✔️✔️:     NodeInfo t;
    // RDKit✔️✔️:     t.id = vect.size();
    // RDKit✔️✔️:     t.in = boost::out_degree(*bNode, *g);
    // RDKit✔️✔️:     t.out = boost::out_degree(*bNode, *g);
    // RDKit✔️✔️:     vect.push_back(t);
    // RDKit✔️✔️:     ++bNode;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::sort(vect.begin(), vect.end(), nodeInfoComp1);
    let mut vect: Vec<NodeInfo> = (0..g.n_atoms)
        .map(|i| {
            let deg = g.out_degree(i);
            NodeInfo {
                id: i,
                in_deg: deg,
                out_deg: deg,
            }
        })
        .collect();
    vect.sort_unstable_by(node_info_cmp1);

    // RDKit✔️✔️:   unsigned int run = 1;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < vect.size(); i += run) {
    // RDKit✔️✔️:     for (run = 1; i+run < vect.size() && vect[i+run].in == vect[i].in
    // RDKit✔️✔️:             && vect[i+run].out == vect[i].out; ++run) { ; }
    // RDKit✔️✔️:     for (unsigned int j = 0; j < run; ++j) {
    // RDKit✔️✔️:       vect[i+j].in += vect[i+j].out;  // in=valence sum
    // RDKit✔️✔️:       vect[i+j].out = run;            // out=frequency
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut i = 0;
    while i < vect.len() {
        let mut run = 1;
        while i + run < vect.len()
            && vect[i + run].in_deg == vect[i].in_deg
            && vect[i + run].out_deg == vect[i].out_deg
        {
            run += 1;
        }
        for j in 0..run {
            vect[i + j].in_deg += vect[i + j].out_deg; // valence sum
            vect[i + j].out_deg = run; // frequency
        }
        i += run;
    }

    // RDKit✔️✔️:   std::sort(vect.begin(), vect.end(), nodeInfoComp2);
    vect.sort_unstable_by(node_info_cmp2);

    // RDKit✔️✔️:   node_id *nodes = new node_id[vect.size()];
    // RDKit✔️✔️:   for (unsigned int i = 0; i < vect.size(); ++i) {
    // RDKit✔️✔️:     nodes[i] = vect[i].id;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nodes;
    vect.iter().map(|ni| ni.id).collect()
}

// RDKit source (vf2.hpp), VF2SubState class:
//   template <class Graph, class VertexCompatible, class EdgeCompatible,
//             class MatchChecking>
//   class VF2SubState {
//    private:
//     Graph *g1, *g2;
//     VertexCompatible &vc;
//     EdgeCompatible &ec;
//     MatchChecking &mc;
//     unsigned int n1, n2;
//     unsigned int core_len;
//     unsigned int t1_len;
//     unsigned int t2_len;  // Core nodes are also counted by these...
//     node_id *core_1;
//     node_id *core_2;
//     node_id *term_1;
//     node_id *term_2;
//     node_id *order;
//     long *share_count;
//     int *vs_compared;

/// RDKit❗✔️: VF2 subgraph isomorphism state.
///
/// g1 = query graph, g2 = molecule graph.
/// core_1[i] = mapping from query atom i -> mol atom j (or None).
/// core_2[j] = mapping from mol atom j -> query atom i (or None).
/// term_1[i] = depth (core_len) when atom i entered terminal set (0 = not terminal).
/// term_2[j] = same for mol atoms.
struct Vf2SubState<'a> {
    g1: &'a Vf2Graph,
    g2: &'a Vf2Graph,
    n1: usize,
    n2: usize,
    core_len: usize,
    t1_len: usize,
    t2_len: usize,
    core_1: Vec<NodeId>,
    core_2: Vec<NodeId>,
    term_1: Vec<usize>,
    term_2: Vec<usize>,
    order: Option<Vec<NodeId>>,
}

impl<'a> Vf2SubState<'a> {
    // RDKit source (vf2.hpp):
    //   VF2SubState(Graph *ag1, Graph *ag2, VertexCompatible &avc,
    //               EdgeCompatible &aec, MatchChecking &amc, bool sortNodes=false)
    //       : g1(ag1), g2(ag2), vc(avc), ec(aec), mc(amc),
    //         n1(num_vertices(*ag1)), n2(num_vertices(*ag2)) {
    //     if (sortNodes) { order = SortNodesByFrequency(ag1); }
    //     else { order = nullptr; }
    //     core_len = 0; t1_len = 0; t2_len = 0;
    //     core_1 = new node_id[n1]; core_2 = new node_id[n2];
    //     term_1 = new node_id[n1]; term_2 = new node_id[n2];
    //     for (unsigned int i = 0; i < n1; i++) {
    //       core_1[i] = NULL_NODE; term_1[i] = 0;
    //     }
    //     for (unsigned int i = 0; i < n2; i++) {
    //       core_2[i] = NULL_NODE; term_2[i] = 0;
    //     }
    //   }

    /// RDKit✔️❌: Construct VF2 state. `sort_nodes` enables frequency-based
    ///   query node ordering. Performance gap: Rust Vec allocation vs C++
    ///   raw new[] — negligible for typical molecule sizes.
    fn new(g1: &'a Vf2Graph, g2: &'a Vf2Graph, sort_nodes: bool) -> Self {
        let n1 = g1.n_atoms;
        let n2 = g2.n_atoms;
        let order = if sort_nodes {
            Some(sort_nodes_by_frequency(g1))
        } else {
            None
        };

        // RDKit✔️✔️: core_len = 0; t1_len = 0; t2_len = 0;
        // RDKit✔️✔️: core_1[i] = NULL_NODE; term_1[i] = 0;
        // RDKit✔️✔️: core_2[j] = NULL_NODE; term_2[j] = 0;
        Self {
            g1,
            g2,
            n1,
            n2,
            core_len: 0,
            t1_len: 0,
            t2_len: 0,
            core_1: vec![NULL_NODE; n1],
            core_2: vec![NULL_NODE; n2],
            term_1: vec![0usize; n1],
            term_2: vec![0usize; n2],
            order,
        }
    }

    fn debug_order(&self) -> Option<&[NodeId]> {
        self.order.as_deref()
    }

    // RDKit source (vf2.hpp):
    //   bool IsGoal() { return core_len == n1; }
    //   bool IsDead() { return n1 > n2 || t1_len > t2_len; }

    /// RDKit✔️✔️: IsGoal — all query atoms matched.
    fn is_goal(&self) -> bool {
        self.core_len == self.n1
    }

    /// RDKit✔️✔️: IsDead — more query atoms than mol atoms, or more
    ///   terminal query atoms than terminal mol atoms.
    fn is_dead(&self) -> bool {
        self.n1 > self.n2 || self.t1_len > self.t2_len
    }

    // RDKit source (vf2.hpp):
    //   bool NextPair(Pair<Graph> &pair) {
    //     if (pair.n1 == NULL_NODE) { pair.n1 = 0; }
    //     if (pair.n2 == NULL_NODE) { pair.n2 = 0; }
    //     else { pair.n2++; }
    //     ...
    //     if (t1_len > core_len && t2_len > core_len) {
    //       while (pair.n1 < n1 &&
    //              (core_1[pair.n1] != NULL_NODE || term_1[pair.n1] == 0)) {
    //         pair.n1++; pair.n2 = 0;
    //       }
    //       ...
    //     } else if (pair.n1 == 0 && order != nullptr) {
    //       // Optimisation: ...
    //       unsigned int i = 0;
    //       while (i < n1 && core_1[pair.n1 = order[i]] != NULL_NODE) { i++; }
    //       ...
    //     } else {
    //       while (pair.n1 < n1 && core_1[pair.n1] != NULL_NODE) {
    //         pair.n1++; pair.n2 = 0;
    //       }
    //     }
    //     // VF2 Plus iterator ...
    //     if (pair.hasiter) { ... }
    //     else if (t1_len > core_len && t2_len > core_len) {
    //       while (pair.n2 < n2 &&
    //              (core_2[pair.n2] != NULL_NODE || term_2[pair.n2] == 0)) {
    //         pair.n2++;
    //       }
    //     } else {
    //       while (pair.n2 < n2 && core_2[pair.n2] != NULL_NODE) { pair.n2++; }
    //     }
    //     return pair.n1 < n1 && pair.n2 < n2;
    //   }

    /// RDKit✔️❌: NextPair — find the next candidate pair (n1 from query,
    ///   n2 from mol) to try matching.
    ///
    /// Uses terminal-set-based iteration from vf2.hpp, including the VF2+
    /// neighbor iterator that restricts mol-side candidates to neighbors of
    /// the already-mapped terminal predecessor.
    fn next_pair(&self, pair: &mut Vf2Pair) -> bool {
        // RDKit✔️✔️: if (pair.n1 == NULL_NODE) pair.n1 = 0;
        // RDKit✔️✔️: if (pair.n2 == NULL_NODE) pair.n2 = 0;
        // RDKit✔️✔️: else pair.n2++;
        if pair.n1 == NULL_NODE {
            pair.n1 = 0;
        }
        if pair.n2 == NULL_NODE {
            pair.n2 = 0;
        } else {
            pair.n2 += 1;
        }

        // --- Select query node (n1) ---
        // RDKit✔️✔️: if (t1_len > core_len && t2_len > core_len) {
        if self.t1_len > self.core_len && self.t2_len > self.core_len {
            // RDKit✔️✔️: while (pair.n1 < n1 &&
            // RDKit✔️✔️:   (core_1[pair.n1] != NULL_NODE || term_1[pair.n1] == 0)) {
            // RDKit✔️✔️:   pair.n1++; pair.n2 = 0;
            // RDKit✔️✔️: }
            while pair.n1 < self.n1
                && (self.core_1[pair.n1] != NULL_NODE || self.term_1[pair.n1] == 0)
            {
                pair.n1 += 1;
                pair.n2 = 0;
            }
            // RDKit✔️✔️: /* Initialize VF2 Plus neighbor iterator.
            // RDKit✔️✔️:  * The next query node (pair.n1) has been selected from the terminal
            // RDKit✔️✔️:  * set and is therefore adjacent to an already mapped atom (in
            // RDKit✔️✔️:  * core_1). Rather than select pair.n2 from all atoms (0...n2) we can
            // RDKit✔️✔️:  * select it from the neighbors of this mapped atom (0...deg(nbor))
            // RDKit✔️✔️:  * since it must also be adajcent to this mapped atom!
            // RDKit✔️✔️:  */
            // RDKit✔️✔️: if (!pair.hasiter) {
            // RDKit✔️✔️:   boost::tie(n1iter_beg, n1iter_end) =
            // RDKit✔️✔️:       boost::adjacent_vertices(pair.n1, *g1);
            // RDKit✔️✔️:   while (n1iter_beg != n1iter_end && core_1[*n1iter_beg] == NULL_NODE) {
            // RDKit✔️✔️:     ++n1iter_beg;
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   assert(n1iter_beg != n1iter_end);
            // RDKit✔️✔️:   boost::tie(pair.nbrbeg, pair.nbrend) =
            // RDKit✔️✔️:       boost::adjacent_vertices(core_1[*n1iter_beg], *g2);
            // RDKit✔️✔️:   pair.hasiter = true;
            // RDKit✔️✔️: }
            if !pair.hasiter {
                let mut mapped_terminal_neighbor = NULL_NODE;
                for &(query_neighbor, _) in self.g1.out_edges(pair.n1) {
                    if self.core_1[query_neighbor] != NULL_NODE {
                        mapped_terminal_neighbor = self.core_1[query_neighbor];
                        break;
                    }
                }
                debug_assert_ne!(mapped_terminal_neighbor, NULL_NODE);
                if mapped_terminal_neighbor != NULL_NODE {
                    pair.nbr_node = mapped_terminal_neighbor;
                    pair.nbr_cursor = 0;
                    pair.nbr_end = self.g2.out_edges(mapped_terminal_neighbor).len();
                    pair.hasiter = true;
                }
            }
        } else if pair.n1 == 0 {
            // RDKit✔️✔️: } else if (pair.n1 == 0 && order != nullptr) {
            if let Some(order) = &self.order {
                // RDKit✔️✔️:   unsigned int i = 0;
                // RDKit✔️✔️:   while (i < n1 && core_1[pair.n1 = order[i]] != NULL_NODE) { i++; }
                // RDKit✔️✔️:   if (i == n1) pair.n1 = n1;
                let mut i = 0;
                while i < self.n1 {
                    let candidate = order[i];
                    if self.core_1[candidate] == NULL_NODE {
                        pair.n1 = candidate;
                        break;
                    }
                    i += 1;
                }
                if i == self.n1 {
                    pair.n1 = self.n1;
                }
            } else {
                // RDKit✔️✔️: } else {
                // RDKit✔️✔️:   while (pair.n1 < n1 && core_1[pair.n1] != NULL_NODE) {
                // RDKit✔️✔️:     pair.n1++; pair.n2 = 0;
                // RDKit✔️✔️:   }
                while pair.n1 < self.n1 && self.core_1[pair.n1] != NULL_NODE {
                    pair.n1 += 1;
                    pair.n2 = 0;
                }
            }
        } else {
            // RDKit✔️✔️: } else {
            // RDKit✔️✔️:   while (pair.n1 < n1 && core_1[pair.n1] != NULL_NODE) {
            // RDKit✔️✔️:     pair.n1++; pair.n2 = 0;
            // RDKit✔️✔️:   }
            while pair.n1 < self.n1 && self.core_1[pair.n1] != NULL_NODE {
                pair.n1 += 1;
                pair.n2 = 0;
            }
        }

        // --- Select mol node (n2) ---
        // RDKit✔️✔️: if (pair.hasiter) { ... }
        if pair.hasiter {
            // RDKit✔️✔️: while (pair.nbrbeg < pair.nbrend && core_2[*pair.nbrbeg] != NULL_NODE) {
            // RDKit✔️✔️:   ++pair.nbrbeg;
            // RDKit✔️✔️: }
            let neighbors = self.g2.out_edges(pair.nbr_node);
            while pair.nbr_cursor < pair.nbr_end
                && self.core_2[neighbors[pair.nbr_cursor].0] != NULL_NODE
            {
                pair.nbr_cursor += 1;
            }
            // RDKit✔️✔️: if (pair.nbrbeg < pair.nbrend) {
            // RDKit✔️✔️:   pair.n2 = *pair.nbrbeg;
            // RDKit✔️✔️:   ++pair.nbrbeg;
            // RDKit✔️✔️: } else {
            // RDKit✔️✔️:   pair.n2 = n2;
            // RDKit✔️✔️: }
            if pair.nbr_cursor < pair.nbr_end {
                pair.n2 = neighbors[pair.nbr_cursor].0;
                pair.nbr_cursor += 1;
            } else {
                pair.n2 = self.n2;
            }
        } else if self.t1_len > self.core_len && self.t2_len > self.core_len {
            // RDKit✔️✔️: } else if (t1_len > core_len && t2_len > core_len) {
            // RDKit✔️✔️:   while (pair.n2 < n2 &&
            // RDKit✔️✔️:     (core_2[pair.n2] != NULL_NODE || term_2[pair.n2] == 0)) {
            // RDKit✔️✔️:     pair.n2++;
            // RDKit✔️✔️:   }
            while pair.n2 < self.n2
                && (self.core_2[pair.n2] != NULL_NODE || self.term_2[pair.n2] == 0)
            {
                pair.n2 += 1;
            }
        } else {
            // RDKit✔️✔️: } else {
            // RDKit✔️✔️:   while (pair.n2 < n2 && core_2[pair.n2] != NULL_NODE) { pair.n2++; }
            // RDKit✔️✔️: }
            while pair.n2 < self.n2 && self.core_2[pair.n2] != NULL_NODE {
                pair.n2 += 1;
            }
        }

        // RDKit✔️✔️: return pair.n1 < n1 && pair.n2 < n2;
        let ok = pair.n1 < self.n1 && pair.n2 < self.n2;
        if ok
            && env::var("RDKIT_ROW61_TRACE").ok().as_deref() == Some("1")
            && self.n1 == 5
            && self.n2 == 26
        {
            println!(
                "row61_substruct_next_pair core_len={} t1_len={} t2_len={} hasiter={} n1={} n2={}",
                self.core_len, self.t1_len, self.t2_len, pair.hasiter, pair.n1, pair.n2
            );
        }
        ok
    }

    // RDKit source (vf2.hpp), IsFeasiblePair:
    //   bool IsFeasiblePair(node_id node1, node_id node2) {
    //     assert(node1 < n1); assert(node2 < n2);
    //     assert(core_1[node1] == NULL_NODE); assert(core_2[node2] == NULL_NODE);
    //
    //     // O(1) check for adjacency list
    //     if (boost::out_degree(node1, *g1) > boost::out_degree(node2, *g2)) {
    //       return false;
    //     }
    //     if (!vc(node1, node2)) { return false; }
    //
    //     unsigned int other1, other2;
    //     // Check the out edges of node1
    //     typename Graph::out_edge_iterator bNbrs, eNbrs;
    //     boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
    //     while (bNbrs != eNbrs) {
    //       other1 = getOtherIdx(*g1, *bNbrs, node1);
    //       if (core_1[other1] != NULL_NODE) {
    //         other2 = core_1[other1];
    //         typename Graph::edge_descriptor oEdge;
    //         bool found;
    //         boost::tie(oEdge, found) = boost::edge(node2, other2, *g2);
    //         if (!found || !ec(*bNbrs, oEdge)) { return false; }
    //       }
    //       ++bNbrs;
    //     }
    //     return true;
    //   }

    /// RDKit✔️❌: IsFeasiblePair — check if (node1, node2) can be added.
    ///
    /// Performs degree check, vertex compatibility, and edge compatibility
    /// for already-matched neighbors. RDK_VF2_PRUNING (terminal count
    /// pre-check) is not enabled — the C++ code also has it behind an
    /// ifdef that is not defined at the top of vf2.hpp.
    fn is_feasible_pair(
        &self,
        node1: NodeId,
        node2: NodeId,
        atom_fn: &impl Fn(usize, usize) -> bool,
        bond_fn: &impl Fn(usize, usize) -> bool,
    ) -> bool {
        // RDKit✔️✔️: assert(node1 < n1); assert(node2 < n2);
        // RDKit✔️✔️: assert(core_1[node1] == NULL_NODE);
        // RDKit✔️✔️: assert(core_2[node2] == NULL_NODE);
        debug_assert!(node1 < self.n1);
        debug_assert!(node2 < self.n2);
        debug_assert_eq!(self.core_1[node1], NULL_NODE);
        debug_assert_eq!(self.core_2[node2], NULL_NODE);
        if self.core_1[node1] != NULL_NODE || self.core_2[node2] != NULL_NODE {
            return false;
        }

        // RDKit✔️✔️: if (boost::out_degree(node1, *g1) > boost::out_degree(node2, *g2)) {
        // RDKit✔️✔️:   return false;
        // RDKit✔️✔️: }
        if self.g1.out_degree(node1) > self.g2.out_degree(node2) {
            return false;
        }

        // RDKit✔️✔️: if (!vc(node1, node2)) { return false; }
        if !atom_fn(node1, node2) {
            return false;
        }

        // RDKit✔️✔️: // Check the out edges of node1
        // RDKit✔️✔️: boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
        // RDKit✔️✔️: while (bNbrs != eNbrs) {
        // RDKit✔️✔️:   other1 = getOtherIdx(*g1, *bNbrs, node1);
        // RDKit✔️✔️:   if (core_1[other1] != NULL_NODE) {
        // RDKit✔️✔️:     other2 = core_1[other1];
        // RDKit✔️✔️:     if (!found || !ec(*bNbrs, oEdge)) { return false; }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   ++bNbrs;
        // RDKit✔️✔️: }
        for &(other1, edge_idx1) in self.g1.out_edges(node1) {
            if other1 == node1 {
                continue;
            }
            if self.core_1[other1] != NULL_NODE {
                let other2 = self.core_1[other1];
                // Check that (node2, other2) has a matching bond.
                let bond_found = self.find_bond(node2, other2);
                match bond_found {
                    Some(edge_idx2) => {
                        if !bond_fn(edge_idx1, edge_idx2) {
                            return false;
                        }
                    }
                    None => return false,
                }
            }
        }

        true
    }

    /// Find a bond between atom `a` and `b` in the molecule graph (g2).
    fn find_bond(&self, a: NodeId, b: NodeId) -> Option<usize> {
        for &(nbr, bond_idx) in self.g2.out_edges(a) {
            if nbr == b {
                return Some(bond_idx);
            }
        }
        None
    }

    // RDKit source (vf2.hpp), AddPair:
    //   void AddPair(node_id node1, node_id node2) {
    //     assert(node1 < n1); assert(node2 < n2);
    //     assert(core_len < n1); assert(core_len < n2);
    //     ++core_len;
    //     if (!term_1[node1]) { term_1[node1] = core_len; ++t1_len; }
    //     if (!term_2[node2]) { term_2[node2] = core_len; ++t2_len; }
    //     core_1[node1] = node2; core_2[node2] = node1;
    //     // FIX: this is explicitly ignoring directionality
    //     // out_edges of node1...
    //     while (bNbrs != eNbrs) {
    //       unsigned int other = getOtherIdx(*g1, *bNbrs, node1);
    //       if (!term_1[other]) { term_1[other] = core_len; ++t1_len; }
    //       ++bNbrs;
    //     }
    //     // out_edges of node2...
    //     while (bNbrs != eNbrs) {
    //       unsigned int other = getOtherIdx(*g2, *bNbrs, node2);
    //       if (!term_2[other]) { term_2[other] = core_len; ++t2_len; }
    //       ++bNbrs;
    //     }
    //   }

    /// RDKit✔️✔️: AddPair — add (node1, node2) to the mapping.
    ///   Updates terminal sets with depth (core_len) tracking.
    fn add_pair(&mut self, node1: NodeId, node2: NodeId) {
        // RDKit✔️✔️: ++core_len;
        self.core_len += 1;
        let depth = self.core_len;

        // RDKit✔️✔️: if (!term_1[node1]) { term_1[node1] = core_len; ++t1_len; }
        if self.term_1[node1] == 0 {
            self.term_1[node1] = depth;
            self.t1_len += 1;
        }

        // RDKit✔️✔️: if (!term_2[node2]) { term_2[node2] = core_len; ++t2_len; }
        if self.term_2[node2] == 0 {
            self.term_2[node2] = depth;
            self.t2_len += 1;
        }

        // RDKit✔️✔️: core_1[node1] = node2; core_2[node2] = node1;
        self.core_1[node1] = node2;
        self.core_2[node2] = node1;

        // RDKit✔️✔️: // FIX: explicitly ignoring directionality
        // RDKit✔️✔️: boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
        // RDKit✔️✔️: while (bNbrs != eNbrs) {
        // RDKit✔️✔️:   unsigned int other = getOtherIdx(*g1, *bNbrs, node1);
        // RDKit✔️✔️:   if (!term_1[other]) { term_1[other] = core_len; ++t1_len; }
        // RDKit✔️✔️:   ++bNbrs;
        // RDKit✔️✔️: }
        for &(other, _) in self.g1.out_edges(node1) {
            if other == node1 {
                continue;
            }
            if self.term_1[other] == 0 {
                self.term_1[other] = depth;
                self.t1_len += 1;
            }
        }

        // RDKit✔️✔️: boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
        // RDKit✔️✔️: while (bNbrs != eNbrs) {
        // RDKit✔️✔️:   unsigned int other = getOtherIdx(*g2, *bNbrs, node2);
        // RDKit✔️✔️:   if (!term_2[other]) { term_2[other] = core_len; ++t2_len; }
        // RDKit✔️✔️:   ++bNbrs;
        // RDKit✔️✔️: }
        for &(other, _) in self.g2.out_edges(node2) {
            if other == node2 {
                continue;
            }
            if self.term_2[other] == 0 {
                self.term_2[other] = depth;
                self.t2_len += 1;
            }
        }
    }

    // RDKit source (vf2.hpp), BackTrack:
    //   void BackTrack(node_id node1, node_id node2) {
    //     if (term_1[node1] == core_len) { term_1[node1] = 0; --t1_len; }
    //     boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
    //     while (bNbrs != eNbrs) {
    //       unsigned int other = getOtherIdx(*g1, *bNbrs, node1);
    //       if (term_1[other] == core_len) { term_1[other] = 0; --t1_len; }
    //       ++bNbrs;
    //     }
    //     if (term_2[node2] == core_len) { term_2[node2] = 0; --t2_len; }
    //     boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
    //     while (bNbrs != eNbrs) {
    //       unsigned int other = getOtherIdx(*g2, *bNbrs, node2);
    //       if (term_2[other] == core_len) { term_2[other] = 0; --t2_len; }
    //       ++bNbrs;
    //     }
    //     core_1[node1] = NULL_NODE;
    //     core_2[node2] = NULL_NODE;
    //     --core_len;
    //   }

    /// RDKit✔️✔️: BackTrack — remove (node1, node2) from mapping.
    ///   Uses depth-based term_X tracking: only clears entries tagged
    ///   with the current core_len. This avoids recomputing terminal
    ///   sets from scratch.
    fn back_track(&mut self, node1: NodeId, node2: NodeId) {
        let depth = self.core_len;

        // RDKit✔️✔️: if (term_1[node1] == core_len) { term_1[node1] = 0; --t1_len; }
        if self.term_1[node1] == depth {
            self.term_1[node1] = 0;
            self.t1_len -= 1;
        }

        // RDKit✔️✔️: boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
        // RDKit✔️✔️: while (bNbrs != eNbrs) {
        // RDKit✔️✔️:   unsigned int other = getOtherIdx(*g1, *bNbrs, node1);
        // RDKit✔️✔️:   if (term_1[other] == core_len) { term_1[other] = 0; --t1_len; }
        // RDKit✔️✔️:   ++bNbrs;
        // RDKit✔️✔️: }
        for &(other, _) in self.g1.out_edges(node1) {
            if other == node1 {
                continue;
            }
            if self.term_1[other] == depth {
                self.term_1[other] = 0;
                self.t1_len -= 1;
            }
        }

        // RDKit✔️✔️: if (term_2[node2] == core_len) { term_2[node2] = 0; --t2_len; }
        if self.term_2[node2] == depth {
            self.term_2[node2] = 0;
            self.t2_len -= 1;
        }

        // RDKit✔️✔️: boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
        // RDKit✔️✔️: while (bNbrs != eNbrs) {
        // RDKit✔️✔️:   unsigned int other = getOtherIdx(*g2, *bNbrs, node2);
        // RDKit✔️✔️:   if (term_2[other] == core_len) { term_2[other] = 0; --t2_len; }
        // RDKit✔️✔️:   ++bNbrs;
        // RDKit✔️✔️: }
        for &(other, _) in self.g2.out_edges(node2) {
            if other == node2 {
                continue;
            }
            if self.term_2[other] == depth {
                self.term_2[other] = 0;
                self.t2_len -= 1;
            }
        }

        // RDKit✔️✔️: core_1[node1] = NULL_NODE;
        // RDKit✔️✔️: core_2[node2] = NULL_NODE;
        // RDKit✔️✔️: --core_len;
        self.core_1[node1] = NULL_NODE;
        self.core_2[node2] = NULL_NODE;
        self.core_len -= 1;
    }

    // RDKit source (vf2.hpp), GetCoreSet:
    //   void GetCoreSet(node_id c1[], node_id c2[]) {
    //     unsigned int i, j;
    //     for (i = 0, j = 0; i < n1; ++i) {
    //       if (core_1[i] != NULL_NODE) { c1[j] = i; c2[j] = core_1[i]; ++j; }
    //     }
    //   }

    /// RDKit✔️✔️: GetCoreSet — extract current mapping arrays.
    fn get_core_set(&self) -> (Vec<NodeId>, Vec<NodeId>) {
        let mut c1 = Vec::with_capacity(self.core_len);
        let mut c2 = Vec::with_capacity(self.core_len);
        // RDKit✔️✔️: for (i = 0, j = 0; i < n1; ++i) {
        // RDKit✔️✔️:   if (core_1[i] != NULL_NODE) { c1[j] = i; c2[j] = core_1[i]; ++j; }
        // RDKit✔️✔️: }
        for i in 0..self.n1 {
            if self.core_1[i] != NULL_NODE {
                c1.push(i);
                c2.push(self.core_1[i]);
            }
        }
        (c1, c2)
    }
}

// ---------------------------------------------------------------------------
// VF2 recursive matching
// ---------------------------------------------------------------------------
//
// RDKit source (vf2.hpp):
//   bool Match(node_id c1[], node_id c2[]) {
//     if (IsGoal()) { GetCoreSet(c1, c2); if (MatchChecks(c1, c2)) return true; }
//     if (IsDead()) return false;
//     Pair<Graph> pair;
//     while (NextPair(pair)) {
//       if (IsFeasiblePair(pair.n1, pair.n2)) {
//         AddPair(pair.n1, pair.n2);
//         if (Match(c1, c2)) return true;  // recurse
//         BackTrack(pair.n1, pair.n2);
//       }
//     }
//     return false;
//   }

/// RDKit✔️❌: Match — find first match via VF2 recursion.
///
/// Matches RDKit's `Match(c1, c2)` entry point. `match_check` allows
/// final verification (like MolMatchFinalCheckFunctor). If None, all
/// completed matches are accepted.
fn vf2_match(
    state: &mut Vf2SubState,
    atom_fn: &impl Fn(usize, usize) -> bool,
    bond_fn: &impl Fn(usize, usize) -> bool,
    mut match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
) -> Option<(Vec<NodeId>, Vec<NodeId>)> {
    // RDKit✔️✔️: if (IsGoal()) { GetCoreSet(c1, c2); if (MatchChecks(c1, c2)) return true; }
    if state.is_goal() {
        let (c1, c2) = state.get_core_set();
        if match_check.as_mut().map_or(true, |check| check(&c1, &c2)) {
            return Some((c1, c2));
        }
    }

    // RDKit✔️✔️: if (IsDead()) return false;
    if state.is_dead() {
        return None;
    }

    // RDKit✔️✔️: Pair<Graph> pair;
    // RDKit✔️✔️: while (NextPair(pair)) {
    // RDKit✔️✔️:   if (IsFeasiblePair(pair.n1, pair.n2)) {
    // RDKit✔️✔️:     AddPair(pair.n1, pair.n2);
    // RDKit✔️✔️:     if (Match(c1, c2)) return true;  // recurse
    // RDKit✔️✔️:     BackTrack(pair.n1, pair.n2);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return false;
    let mut pair = Vf2Pair::new();
    while state.next_pair(&mut pair) {
        if state.is_feasible_pair(pair.n1, pair.n2, atom_fn, bond_fn) {
            state.add_pair(pair.n1, pair.n2);
            if let Some(result) = vf2_match(state, atom_fn, bond_fn, match_check.as_deref_mut()) {
                return Some(result);
            }
            state.back_track(pair.n1, pair.n2);
        }
    }
    None
}

// RDKit source (vf2.hpp), MatchAll:
//   template <class DoubleBackInsertionSequence>
//   bool MatchAll(node_id c1[], node_id c2[], DoubleBackInsertionSequence &res,
//                 unsigned int lim = 0) {
//     if (IsGoal()) {
//       GetCoreSet(c1, c2);
//       if (MatchChecks(c1, c2)) {
//         typename DoubleBackInsertionSequence::value_type newSeq;
//         newSeq.reserve(core_len);
//         for (unsigned int i = 0; i < core_len; ++i) {
//           newSeq.emplace_back(c1[i], c2[i]);
//         }
//         res.push_back(newSeq);
//         return lim && res.size() >= lim;
//       }
//     }
//     if (IsDead()) return false;
//     Pair<Graph> pair;
//     while (NextPair(pair)) {
//       if (IsFeasiblePair(pair.n1, pair.n2)) {
//         AddPair(pair.n1, pair.n2);
//         if (MatchAll(c1, c2, res, lim)) return true;  // recurse
//         BackTrack(pair.n1, pair.n2);
//       }
//     }
//     return false;
//   }

/// RDKit✔️❌: MatchAll — find all matches up to `max_matches`.
///
/// Collects matches into `results` as (c1, c2) pairs.
/// Returns true when the limit has been reached, signaling the caller
/// to stop.
fn vf2_match_all(
    state: &mut Vf2SubState,
    atom_fn: &impl Fn(usize, usize) -> bool,
    bond_fn: &impl Fn(usize, usize) -> bool,
    mut match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
    results: &mut Vec<(Vec<NodeId>, Vec<NodeId>)>,
    max_matches: usize,
) -> bool {
    // RDKit✔️✔️: if (IsGoal()) { ... }
    if state.is_goal() {
        let (c1, c2) = state.get_core_set();
        if match_check.as_mut().map_or(true, |check| check(&c1, &c2)) {
            // RDKit✔️✔️: newSeq.reserve(core_len);
            // RDKit✔️✔️: for (unsigned int i = 0; i < core_len; ++i) {
            // RDKit✔️✔️:   newSeq.emplace_back(c1[i], c2[i]);
            // RDKit✔️✔️: }
            // RDKit✔️✔️: res.push_back(newSeq);
            results.push((c1, c2));
            // RDKit✔️✔️: return lim && res.size() >= lim;
            return max_matches > 0 && results.len() >= max_matches;
        }
    }

    // RDKit✔️✔️: if (IsDead()) return false;
    if state.is_dead() {
        return false;
    }

    // RDKit✔️✔️: Pair<Graph> pair;
    // RDKit✔️✔️: while (NextPair(pair)) { ... }
    let mut pair = Vf2Pair::new();
    while state.next_pair(&mut pair) {
        if state.is_feasible_pair(pair.n1, pair.n2, atom_fn, bond_fn) {
            state.add_pair(pair.n1, pair.n2);
            // RDKit✔️✔️: if (MatchAll(c1, c2, res, lim)) return true;
            if vf2_match_all(
                state,
                atom_fn,
                bond_fn,
                match_check.as_deref_mut(),
                results,
                max_matches,
            ) {
                return true;
            }
            state.back_track(pair.n1, pair.n2);
        }
    }
    false
}

// ---------------------------------------------------------------------------
// Final match check (simplified MolMatchFinalCheckFunctor)
// ---------------------------------------------------------------------------
//
// RDKit source (SubstructMatch.cpp):
//   bool MolMatchFinalCheckFunctor::operator()(const std::uint32_t q_c[],
//                                              const std::uint32_t m_c[]) {
//     if (d_params.extraFinalCheck || d_params.useGenericMatchers) { ... }
//     HashedStorageType match;
//     if (d_params.uniquify) {
//       match.resize(d_mol.getNumAtoms());
//       std::fill(match.begin(), match.end(), 0);
//       for (unsigned int i = 0; i < d_query.getNumAtoms(); ++i) {
//         match[m_c[i]] = 1;
//       }
//       if (matchesSeen.find(match) != matchesSeen.end()) { return false; }
//     }
//     if (!d_params.useChirality) {
//       if (d_params.uniquify) { matchesSeen.insert(match); }
//       return true;
//     }
//     // ... chirality checks ...
//   }

/// RDKit✔️✔️: Final match atom-set mask used for uniquification.
fn match_mask(atom_mapping: &[usize], mol_num_atoms: usize) -> Vec<bool> {
    let mut mask = vec![false; mol_num_atoms];
    for &ma in atom_mapping {
        if ma < mol_num_atoms {
            mask[ma] = true;
        }
    }
    mask
}

fn count_swaps_to_interconvert_i32(reference: &[i32], probe: &[i32]) -> Option<u32> {
    // BEGIN RDKIT CPP FUNCTION countSwapsToInterconvert
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: unsigned int countSwapsToInterconvert(const T &ref, T probe) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   PRECONDITION(ref.size() == probe.size(), "bad vector sizes");
    // RDKit✔️✔️:   for (unsigned int i = 0; i < ref.size(); ++i) {
    // RDKit✔️✔️:     if (ref[i] != probe[i]) {
    // RDKit✔️✔️:       unsigned int j = i + 1;
    // RDKit✔️✔️:       while (probe[j] != ref[i]) {
    // RDKit✔️✔️:         ++j;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       std::swap(probe[i], probe[j]);
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    if reference.len() != probe.len() {
        return None;
    }
    let mut work = probe.to_vec();
    let mut swaps = 0_u32;
    for (idx, expected) in reference.iter().copied().enumerate() {
        if work[idx] == expected {
            continue;
        }
        let found_idx = work[idx + 1..]
            .iter()
            .position(|value| *value == expected)
            .map(|offset| idx + 1 + offset)?;
        work.swap(idx, found_idx);
        swaps = swaps.saturating_add(1);
    }
    Some(swaps)
}

fn rdkit_atom_perturbation_order_from_bond_indices(
    mol: &Molecule,
    atom_idx: usize,
    probe: &[i32],
) -> Result<u32, SubstructMatchError> {
    // BEGIN RDKIT CPP FUNCTION Atom::getPerturbationOrder
    // RDKit✔️✔️: int Atom::getPerturbationOrder(const INT_LIST &probe) const {
    // RDKit✔️✔️:   INT_LIST ref;
    // RDKit✔️✔️:   for (const auto bond : getOwningMol().atomBonds(this)) {
    // RDKit✔️✔️:     ref.push_back(bond->getIdx());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return static_cast<int>(countSwapsToInterconvert(probe, ref));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    let reference: Vec<i32> = mol
        .topology_block()
        .adjacency
        .neighbors_of(atom_idx)
        .iter()
        .map(|neighbor| i32::try_from(neighbor.bond.index()))
        .collect::<Result<_, _>>()
        .map_err(|_| SubstructMatchError::Unsupported {
            branch: "MolMatchFinalCheckFunctor/Atom::getPerturbationOrder/bond-index-overflow",
            rdkit_function: "Atom::getPerturbationOrder",
        })?;
    count_swaps_to_interconvert_i32(probe, &reference).ok_or(SubstructMatchError::Unsupported {
        branch: "MolMatchFinalCheckFunctor/Atom::getPerturbationOrder/unmodeled-bond-ordering",
        rdkit_function: "Atom::getPerturbationOrder/countSwapsToInterconvert",
    })
}

fn rdkit_translate_ez_label_to_cis_trans(stereo: BondStereo) -> BondStereo {
    match stereo {
        BondStereo::E => BondStereo::Trans,
        BondStereo::Z => BondStereo::Cis,
        other => other,
    }
}

fn find_bond_between(mol: &Molecule, begin: usize, end: usize) -> Option<&Bond> {
    mol.bonds().iter().find(|bond| {
        let b = bond.begin().index();
        let e = bond.end().index();
        (b == begin && e == end) || (b == end && e == begin)
    })
}

fn rdkit_match_final_check(
    mol: &Molecule,
    query: &Molecule,
    params: &SubstructMatchParams,
    c1: &[NodeId],
    c2: &[NodeId],
    matches_seen: &mut Vec<Vec<bool>>,
) -> Result<bool, SubstructMatchError> {
    // BEGIN RDKIT CPP FUNCTION MolMatchFinalCheckFunctor::operator()
    // RDKit✔️✔️: bool MolMatchFinalCheckFunctor::operator()(const std::uint32_t q_c[],
    // RDKit✔️✔️:                                            const std::uint32_t m_c[]) {
    // RDKit✔️✔️:   HashedStorageType match;
    // RDKit✔️✔️:   if (d_params.uniquify) {
    // RDKit✔️✔️:     match.resize(d_mol.getNumAtoms());
    // RDKit✔️✔️:     std::fill(match.begin(), match.end(), 0);
    // RDKit✔️✔️:     for (unsigned int i = 0; i < d_query.getNumAtoms(); ++i) {
    // RDKit✔️✔️:       match[m_c[i]] = 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (matchesSeen.find(match) != matchesSeen.end()) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut q_to_mol = vec![NULL_NODE; query.num_atoms()];
    for (&qa, &ma) in c1.iter().zip(c2.iter()) {
        if qa < q_to_mol.len() {
            q_to_mol[qa] = ma;
        }
    }
    let match_key = if params.uniquify {
        let mask = match_mask(&q_to_mol, mol.num_atoms());
        if matches_seen.iter().any(|existing| *existing == mask) {
            return Ok(false);
        }
        Some(mask)
    } else {
        None
    };

    // RDKit✔️✔️:   if (!d_params.useChirality) {
    // RDKit✔️✔️:     if (d_params.uniquify) {
    // RDKit✔️✔️:       matchesSeen.insert(match);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    if !params.use_chirality {
        if let Some(mask) = match_key {
            matches_seen.push(mask);
        }
        return Ok(true);
    }

    // RDKit✔️✔️:   std::unordered_map<unsigned int, bool> matches;
    // RDKit❌❌:   if (d_params.useEnhancedStereo) {
    // RDKit❌❌:     if (!detail::enhancedStereoIsOK(...)) { return false; }
    // RDKit❌❌:   }
    //
    // COSMolKit stores enhanced stereo groups, but `useEnhancedStereo` and
    // `enhancedStereoIsOK()` are not yet modeled in substructure matching.
    // Absolute groups are ignored by the RDKit constructor for this path; OR
    // and AND groups would alter final acceptance and must fail closed.
    if mol
        .stereo_groups()
        .iter()
        .any(|group| !matches!(group.kind(), StereoGroupKind::Absolute))
    {
        return Err(SubstructMatchError::Unsupported {
            branch: "MolMatchFinalCheckFunctor/enhancedStereoIsOK",
            rdkit_function: "detail::enhancedStereoIsOK",
        });
    }

    // RDKit✔️✔️:   // check chiral atoms:
    // RDKit✔️✔️:   for (unsigned int i = 0; i < d_query.getNumAtoms(); ++i) {
    // RDKit✔️✔️:     const Atom *qAt = d_query.getAtomWithIdx(q_c[i]);
    // RDKit✔️✔️:     if (qAt->getDegree() < 3 || !detail::hasChiralLabel(qAt)) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    for qi in 0..query.num_atoms() {
        let q_at = &query.atoms()[qi];
        if query.topology_block().adjacency.neighbors_of(qi).len() < 3
            || !has_rdkit_tetrahedral_chiral_label(q_at)
        {
            continue;
        }
        let mi = q_to_mol[qi];
        let m_at = &mol.atoms()[mi];
        // RDKit✔️✔️:     if (!detail::hasChiralLabel(mAt)) {
        // RDKit✔️✔️:       if (d_params.specifiedStereoQueryMatchesUnspecified) {
        // RDKit✔️✔️:         continue;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        if !has_rdkit_tetrahedral_chiral_label(m_at) {
            if params.specified_stereo_query_matches_unspecified {
                continue;
            }
            return Ok(false);
        }
        // RDKit✔️✔️:     if (qAt->getDegree() > mAt->getDegree()) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        if query.topology_block().adjacency.neighbors_of(qi).len()
            > mol.topology_block().adjacency.neighbors_of(mi).len()
        {
            return Ok(false);
        }

        // RDKit✔️✔️:     INT_LIST qOrder;
        // RDKit✔️✔️:     INT_LIST mOrder;
        // RDKit✔️✔️:     for (unsigned int j = 0; j < d_query.getNumAtoms(); ++j) {
        // RDKit✔️✔️:       const Bond *qB = d_query.getBondBetweenAtoms(q_c[i], q_c[j]);
        // RDKit✔️✔️:       const Bond *mB = d_mol.getBondBetweenAtoms(m_c[i], m_c[j]);
        // RDKit✔️✔️:       if (qB && mB) {
        // RDKit✔️✔️:         mOrder.push_back(mB->getIdx());
        // RDKit✔️✔️:         qOrder.push_back(qB->getIdx());
        // RDKit✔️✔️:         if (mOrder.size() == qAt->getDegree()) {
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let mut q_order: Vec<i32> = Vec::new();
        let mut m_order: Vec<i32> = Vec::new();
        for qj in 0..query.num_atoms() {
            let Some(q_bond) = find_bond_between(query, qi, qj) else {
                continue;
            };
            let mj = q_to_mol[qj];
            let Some(m_bond) = find_bond_between(mol, mi, mj) else {
                continue;
            };
            q_order.push(i32::try_from(q_bond.id().index()).map_err(|_| {
                SubstructMatchError::Unsupported {
                    branch: "MolMatchFinalCheckFunctor/qOrder/bond-index-overflow",
                    rdkit_function: "Atom::getPerturbationOrder",
                }
            })?);
            m_order.push(i32::try_from(m_bond.id().index()).map_err(|_| {
                SubstructMatchError::Unsupported {
                    branch: "MolMatchFinalCheckFunctor/mOrder/bond-index-overflow",
                    rdkit_function: "countSwapsToInterconvert",
                }
            })?);
            if m_order.len() == query.topology_block().adjacency.neighbors_of(qi).len() {
                break;
            }
        }
        if q_order.len() != query.topology_block().adjacency.neighbors_of(qi).len()
            || q_order.len() != m_order.len()
        {
            return Err(SubstructMatchError::Unsupported {
                branch: "MolMatchFinalCheckFunctor/chiral-atom-missing-matched-neighbors",
                rdkit_function: "MolMatchFinalCheckFunctor::operator()",
            });
        }
        // RDKit✔️✔️:     int qPermCount = qAt->getPerturbationOrder(qOrder);
        let q_perm_count = rdkit_atom_perturbation_order_from_bond_indices(query, qi, &q_order)?;

        // RDKit✔️✔️:     unsigned unmatchedNeighbors = mAt->getDegree() - mOrder.size();
        // RDKit✔️✔️:     mOrder.insert(mOrder.end(), unmatchedNeighbors, -1);
        let unmatched_neighbors = mol
            .topology_block()
            .adjacency
            .neighbors_of(mi)
            .len()
            .saturating_sub(m_order.len());
        m_order.extend(std::iter::repeat_n(-1, unmatched_neighbors));

        // RDKit✔️✔️:     INT_LIST moOrder;
        // RDKit✔️✔️:     for (const auto &bond : d_mol.atomBonds(mAt)) {
        // RDKit✔️✔️:       const int dbidx = bond->getIdx();
        // RDKit✔️✔️:       if (std::find(mOrder.begin(), mOrder.end(), dbidx) != mOrder.end()) {
        // RDKit✔️✔️:         moOrder.push_back(dbidx);
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         moOrder.push_back(-1);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let mo_order: Vec<i32> = mol
            .topology_block()
            .adjacency
            .neighbors_of(mi)
            .iter()
            .map(|neighbor| {
                i32::try_from(neighbor.bond.index()).map(|bond_idx| {
                    if m_order.contains(&bond_idx) {
                        bond_idx
                    } else {
                        -1
                    }
                })
            })
            .collect::<Result<_, _>>()
            .map_err(|_| SubstructMatchError::Unsupported {
                branch: "MolMatchFinalCheckFunctor/moOrder/bond-index-overflow",
                rdkit_function: "countSwapsToInterconvert",
            })?;
        // RDKit✔️✔️:     const int mPermCount =
        // RDKit✔️✔️:         static_cast<int>(countSwapsToInterconvert(moOrder, mOrder));
        let m_perm_count = count_swaps_to_interconvert_i32(&mo_order, &m_order).ok_or(
            SubstructMatchError::Unsupported {
                branch: "MolMatchFinalCheckFunctor/mPermCount/unmodeled-bond-ordering",
                rdkit_function: "countSwapsToInterconvert",
            },
        )?;

        // RDKit✔️✔️:     const bool requireMatch = qPermCount % 2 == mPermCount % 2;
        // RDKit✔️✔️:     const bool labelsMatch = qAt->getChiralTag() == mAt->getChiralTag();
        // RDKit✔️✔️:     const bool matchOK = requireMatch == labelsMatch;
        // RDKit✔️✔️:     if (!matchOK) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        let require_match = q_perm_count % 2 == m_perm_count % 2;
        let labels_match = q_at.chiral_tag() == m_at.chiral_tag();
        if require_match != labels_match {
            return Ok(false);
        }
    }

    // RDKit✔️✔️:   // now check double bonds
    // RDKit✔️✔️:   for (const auto &qBnd : d_query.bonds()) {
    // RDKit✔️✔️:     if (qBnd->getBondType() != Bond::DOUBLE ||
    // RDKit✔️✔️:         qBnd->getStereo() <= Bond::STEREOANY) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    for q_bnd in query.bonds() {
        if q_bnd.order() != BondOrder::Double || !rdkit_bond_stereo_is_above_any(q_bnd.stereo()) {
            continue;
        }
        // RDKit✔️✔️:     if (qBnd->getStereoAtoms().size() != 2) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        let Some(q_stereo_atoms) = q_bnd.stereo_atoms() else {
            continue;
        };
        // RDKit✔️✔️:     const Bond *mBnd = d_mol.getBondBetweenAtoms(
        // RDKit✔️✔️:         q_to_mol[qBnd->getBeginAtomIdx()], q_to_mol[qBnd->getEndAtomIdx()]);
        let q_begin_mol = q_to_mol[q_bnd.begin().index()];
        let q_end_mol = q_to_mol[q_bnd.end().index()];
        let Some(m_bnd) = find_bond_between(mol, q_begin_mol, q_end_mol) else {
            return Err(SubstructMatchError::Unsupported {
                branch: "MolMatchFinalCheckFunctor/double-bond-matching-bond-missing",
                rdkit_function: "MolMatchFinalCheckFunctor::operator()",
            });
        };
        // RDKit✔️✔️:     if (mBnd->getBondType() != Bond::DOUBLE) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if m_bnd.order() != BondOrder::Double {
            continue;
        }
        // RDKit✔️✔️:     if (!d_params.specifiedStereoQueryMatchesUnspecified &&
        // RDKit✔️✔️:         mBnd->getStereo() <= Bond::STEREOANY) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        if !params.specified_stereo_query_matches_unspecified
            && !rdkit_bond_stereo_is_above_any(m_bnd.stereo())
        {
            return Ok(false);
        }
        // RDKit✔️✔️:     if (mBnd->getStereoAtoms().size() != 2) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        let Some(m_stereo_atoms) = m_bnd.stereo_atoms() else {
            continue;
        };

        // RDKit✔️✔️:     unsigned int end1Matches = 0;
        // RDKit✔️✔️:     unsigned int end2Matches = 0;
        // RDKit✔️✔️:     if (q_to_mol[qBnd->getBeginAtomIdx()] == mBnd->getBeginAtomIdx()) {
        // RDKit✔️✔️:       if (q_to_mol[qBnd->getStereoAtoms()[0]] ==
        // RDKit✔️✔️:           static_cast<unsigned>(mBnd->getStereoAtoms()[0])) {
        // RDKit✔️✔️:         end1Matches = 1;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (q_to_mol[qBnd->getStereoAtoms()[1]] ==
        // RDKit✔️✔️:           static_cast<unsigned>(mBnd->getStereoAtoms()[1])) {
        // RDKit✔️✔️:         end2Matches = 1;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       if (q_to_mol[qBnd->getStereoAtoms()[0]] ==
        // RDKit✔️✔️:           static_cast<unsigned>(mBnd->getStereoAtoms()[1])) {
        // RDKit✔️✔️:         end1Matches = 1;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (q_to_mol[qBnd->getStereoAtoms()[1]] ==
        // RDKit✔️✔️:           static_cast<unsigned>(mBnd->getStereoAtoms()[0])) {
        // RDKit✔️✔️:         end2Matches = 1;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let mut end1_matches = 0_u32;
        let mut end2_matches = 0_u32;
        if q_begin_mol == m_bnd.begin().index() {
            if q_to_mol[q_stereo_atoms[0].index()] == m_stereo_atoms[0].index() {
                end1_matches = 1;
            }
            if q_to_mol[q_stereo_atoms[1].index()] == m_stereo_atoms[1].index() {
                end2_matches = 1;
            }
        } else {
            if q_to_mol[q_stereo_atoms[0].index()] == m_stereo_atoms[1].index() {
                end1_matches = 1;
            }
            if q_to_mol[q_stereo_atoms[1].index()] == m_stereo_atoms[0].index() {
                end2_matches = 1;
            }
        }

        // RDKit✔️✔️:     const unsigned totalMatches = end1Matches + end2Matches;
        // RDKit✔️✔️:     const auto mStereo =
        // RDKit✔️✔️:         Chirality::translateEZLabelToCisTrans(mBnd->getStereo());
        // RDKit✔️✔️:     const auto qStereo =
        // RDKit✔️✔️:         Chirality::translateEZLabelToCisTrans(qBnd->getStereo());
        // RDKit✔️✔️:     if (mStereo == qStereo && totalMatches == 1) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (mStereo != qStereo && totalMatches != 1) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        let total_matches = end1_matches + end2_matches;
        let m_stereo = rdkit_translate_ez_label_to_cis_trans(m_bnd.stereo());
        let q_stereo = rdkit_translate_ez_label_to_cis_trans(q_bnd.stereo());
        if m_stereo == q_stereo && total_matches == 1 {
            return Ok(false);
        }
        if m_stereo != q_stereo && total_matches != 1 {
            return Ok(false);
        }
    }

    // RDKit✔️✔️:   if (d_params.uniquify) {
    // RDKit✔️✔️:     matchesSeen.insert(match);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    if let Some(mask) = match_key {
        matches_seen.push(mask);
    }
    Ok(true)
}

// ---------------------------------------------------------------------------
// Bond mapping builder
// ---------------------------------------------------------------------------

/// Build the bond mapping for a match result.
///
/// For each query bond (by index), find the corresponding molecular bond
/// that connects the matched query endpoints.
#[allow(dead_code)]
fn build_bond_mapping(
    query_atom_to_mol: &[Option<usize>],
    query: &Vf2Graph,
    mol: &Vf2Graph,
) -> Vec<usize> {
    let mut bond_mapping = Vec::with_capacity(query.n_bonds);
    for bond_idx in 0..query.n_bonds {
        // Find the query atoms connected by this bond.
        let mut q_begin = NULL_NODE;
        let mut q_end = NULL_NODE;
        for qa in 0..query.n_atoms {
            for &(nbr, eidx) in &query.adjacency[qa] {
                if eidx == bond_idx {
                    q_begin = qa;
                    q_end = nbr;
                    break;
                }
            }
            if q_begin != NULL_NODE {
                break;
            }
        }

        if q_begin != NULL_NODE {
            let m_begin = query_atom_to_mol[q_begin];
            let m_end = query_atom_to_mol[q_end];
            if let (Some(mb), Some(me)) = (m_begin, m_end) {
                // Find bond between mb and me in mol.
                let mut mol_bond_idx = NULL_NODE;
                for &(nbr, eidx) in &mol.adjacency[mb] {
                    if nbr == me {
                        mol_bond_idx = eidx;
                        break;
                    }
                }
                bond_mapping.push(mol_bond_idx);
            } else {
                bond_mapping.push(NULL_NODE);
            }
        } else {
            bond_mapping.push(NULL_NODE);
        }
    }
    bond_mapping
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// RDKit❌❌: `SubstructMatch` entry point.
///
/// Wires up atom/bond functors and runs VF2. This is a simplified port of
/// RDKit's `SubstructMatch()` from SubstructMatch.cpp (lines 481-525).
///
/// Returns all matches (up to params.max_matches).
fn collect_recursive_smarts_from_atom_query(
    query: &crate::QueryNode<AtomQueryPredicate>,
    recursive_smarts: &mut Vec<String>,
) {
    match query {
        crate::QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(smarts)) => {
            recursive_smarts.push(smarts.clone());
        }
        crate::QueryNode::Predicate(_) => {}
        crate::QueryNode::And(children) | crate::QueryNode::Or(children) => {
            for child in children {
                collect_recursive_smarts_from_atom_query(child, recursive_smarts);
            }
        }
        crate::QueryNode::Not(child) => {
            collect_recursive_smarts_from_atom_query(child, recursive_smarts);
        }
    }
}

fn populate_recursive_query_match_cache(
    mol: &Molecule,
    query: &Molecule,
    params: &SubstructMatchParams,
    recursive_cache: &mut RecursiveQueryMatchCache,
) -> Result<(), SubstructMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructMatch.cpp :: SubstructMatch
    // RDKit✔️✔️:   if (params.recursionPossible) {
    // RDKit✔️✔️:     detail::SUBQUERY_MAP subqueryMap;
    // RDKit✔️✔️:     ROMol::ConstAtomIterator atIt;
    // RDKit✔️✔️:     for (const auto atom : query.atoms()) {
    // RDKit✔️✔️:       if (atom->hasQuery()) {
    // RDKit✔️✔️:         detail::MatchSubqueries(mol, atom->getQuery(), params, subqueryMap,
    // RDKit✔️✔️:                                 locker.locked);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    //
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructMatch.cpp :: detail::MatchSubqueries
    // RDKit✔️✔️:   if (query->getDescription() == "RecursiveStructure") {
    // RDKit✔️✔️:     auto *rsq = (RecursiveStructureQuery *)query;
    // RDKit✔️✔️:     rsq->clear();
    // RDKit✔️✔️:     bool matchDone = false;
    // RDKit✔️✔️:     if (rsq->getSerialNumber() &&
    // RDKit✔️✔️:         subqueryMap.find(rsq->getSerialNumber()) != subqueryMap.end()) {
    // RDKit✔️✔️:       matchDone = true;
    // RDKit✔️✔️:       auto orsq =
    // RDKit✔️✔️:           (const RecursiveStructureQuery *)subqueryMap[rsq->getSerialNumber()];
    // RDKit✔️✔️:       for (auto setIter = orsq->beginSet(); setIter != orsq->endSet();
    // RDKit✔️✔️:            ++setIter) {
    // RDKit✔️✔️:         rsq->insert(*setIter);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (!matchDone) {
    // RDKit✔️✔️:       ROMol const *queryMol = rsq->getQueryMol();
    // RDKit✔️✔️:       if (queryMol) {
    // RDKit✔️✔️:         std::vector<int> matchStarts;
    // RDKit✔️✔️:         unsigned int res = RecursiveMatcher(mol, *queryMol, matchStarts,
    // RDKit✔️✔️:                                             subqueryMap, params, locked);
    // RDKit✔️✔️:         if (res) {
    // RDKit✔️✔️:           for (int &matchStart : matchStarts) {
    // RDKit✔️✔️:             rsq->insert(matchStart);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (rsq->getSerialNumber()) {
    // RDKit✔️✔️:         subqueryMap[rsq->getSerialNumber()] = query;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto childIt = query->beginChildren(); childIt != query->endChildren();
    // RDKit✔️✔️:        ++childIt) {
    // RDKit✔️✔️:     MatchSubqueries(mol, childIt->get(), params, subqueryMap, locked);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    //
    // RDKit serial-number sharing is represented by string-deduplicated cached
    // recursive SMARTS. Each cache value is the RecursiveStructureQuery set of
    // molecule atom indices that can serve as the inner query root.
    let mut recursive_smarts = Vec::new();
    for atom in query.atoms() {
        if let Some(query_node) = atom.query() {
            collect_recursive_smarts_from_atom_query(query_node, &mut recursive_smarts);
        }
    }
    recursive_smarts.sort();
    recursive_smarts.dedup();

    for recursive_smarts in recursive_smarts {
        if recursive_cache.contains_key(&recursive_smarts) {
            continue;
        }
        let Some(inner) = recursive_smarts
            .strip_prefix("$(")
            .and_then(|value| value.strip_suffix(')'))
        else {
            recursive_cache.insert(recursive_smarts, vec![false; mol.num_atoms()]);
            continue;
        };
        let Ok(inner_query) = build_crystalff_query_molecule(inner) else {
            recursive_cache.insert(recursive_smarts, vec![false; mol.num_atoms()]);
            continue;
        };
        populate_recursive_query_match_cache(mol, &inner_query, params, recursive_cache)?;
        let matches = substruct_match_impl_with_recursive_cache(
            mol,
            &inner_query,
            &SubstructMatchParams {
                max_matches: params.max_matches.max(1000),
                uniquify: false,
                use_chirality: params.use_chirality,
                specified_stereo_query_matches_unspecified: params
                    .specified_stereo_query_matches_unspecified,
            },
            Some(recursive_cache),
        )?;
        let mut match_starts = vec![false; mol.num_atoms()];
        for matched in matches {
            if let Some(&root_atom_idx) = matched.atom_mapping.first()
                && root_atom_idx != NULL_NODE
                && root_atom_idx < match_starts.len()
            {
                match_starts[root_atom_idx] = true;
            }
        }
        recursive_cache.insert(recursive_smarts, match_starts);
    }
    Ok(())
}

fn substruct_match_impl_with_recursive_cache(
    mol: &Molecule,
    query: &Molecule,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
) -> SubstructMatchResultList {
    let m_num_atoms = mol.num_atoms();
    let q_num_atoms = query.num_atoms();
    let trace_row64_substruct =
        env::var("RDKIT_ROW64_SUBSTRUCT_TRACE").ok().as_deref() == Some("1") && m_num_atoms == 106;
    let trace_row61_torsion349 = row61_torsion349_substruct_trace_enabled(mol, query, params);
    let query_ctx = build_query_match_context(mol);

    // RDKit source (SubstructMatch.cpp):
    //   if (!mNumAtoms || !qNumAtoms || qNumAtoms > mNumAtoms) {
    //     return matches;
    //   }
    if m_num_atoms == 0 || q_num_atoms == 0 || q_num_atoms > m_num_atoms {
        return Ok(Vec::new());
    }

    // Build VF2 graphs.
    let graph_start = trace_row64_substruct.then(Instant::now);
    let q_graph = build_vf2_graph(query);
    let m_graph = build_vf2_graph(mol);
    let graph_elapsed = graph_start.map(|start| start.elapsed().as_secs_f64());

    // Build atom matching closure.
    // RDKit source:
    //   detail::AtomLabelFunctor atomLabeler(query, mol, params);
    //   detail::BondLabelFunctor bondLabeler(query, mol, params);
    //   MolMatchFinalCheckFunctor matchChecker(query, mol, params);
    let atom_fn = |qi: usize, mj: usize| -> bool {
        let qa = &query.atoms()[qi];
        let ma = &mol.atoms()[mj];
        atom_label_chirality_matches(qa, ma, params)
            && atom_matches_with_recursive_cache(
                qa,
                query,
                ma,
                mol,
                params,
                recursive_cache,
                &query_ctx,
            )
    };

    let bond_fn = |qei: usize, mei: usize| -> bool {
        let qb = &query.bonds()[qei];
        let mb = &mol.bonds()[mei];
        bond_label_chirality_matches(qb, mb, params)
            && if let Some(query_node) = qb.query() {
                evaluate_bond_query(query_node, mb, mol, &query_ctx)
            } else {
                bond_matches(qb, query, mb, mol)
            }
    };

    // RDKit source:
    //   bool found = boost::vf2_all(query.getTopology(), mol.getTopology(),
    //                               atomLabeler, bondLabeler, matchChecker,
    //                               pms, params.maxMatches);
    let mut state = Vf2SubState::new(&q_graph, &m_graph, /*sort_nodes=*/ false);
    if trace_row61_torsion349 {
        println!(
            "row61_substruct_begin q_num_atoms={} q_num_bonds={} m_num_atoms={} order={:?}",
            q_num_atoms,
            query.num_bonds(),
            m_num_atoms,
            state.debug_order()
        );
    }
    let mut raw_matches: Vec<(Vec<NodeId>, Vec<NodeId>)> = Vec::new();
    let mut matches_seen: Vec<Vec<bool>> = Vec::new();
    let mut final_check_error: Option<SubstructMatchError> = None;
    let mut check_fn = |c1: &[NodeId], c2: &[NodeId]| -> bool {
        match rdkit_match_final_check(mol, query, params, c1, c2, &mut matches_seen) {
            Ok(accepted) => accepted,
            Err(err) => {
                final_check_error = Some(err);
                false
            }
        }
    };

    let vf2_start = trace_row64_substruct.then(Instant::now);
    vf2_match_all(
        &mut state,
        &atom_fn,
        &bond_fn,
        Some(&mut check_fn),
        &mut raw_matches,
        params.max_matches,
    );
    if let Some(err) = final_check_error {
        return Err(err);
    }
    if let Some(vf2_start) = vf2_start {
        println!(
            "row64_substruct_core q_atoms={} q_bonds={} graph_build={:.6} vf2={:.6} raw_matches={}",
            q_num_atoms,
            query.num_bonds(),
            graph_elapsed.unwrap_or(0.0),
            vf2_start.elapsed().as_secs_f64(),
            raw_matches.len()
        );
    }

    // RDKit source (SubstructMatch.cpp):
    //   if (found) {
    //     const unsigned int nQueryAtoms = query.getNumAtoms();
    //     matches.reserve(pms.size());
    //     MatchVectType matchVect(nQueryAtoms);
    //     for (const auto &pairs : pms) {
    //       for (const auto &pair : pairs) {
    //         matchVect[pair.first] = pair;
    //       }
    //       matches.push_back(matchVect);
    //     }
    //   }
    let mut results: Vec<SubstructMatchResult> = Vec::new();

    for (c1, c2) in &raw_matches {
        if trace_row61_torsion349 {
            println!("row61_substruct_raw_core c1={c1:?} c2={c2:?}");
        }
        // Build atom_mapping: query_atom_index -> mol_atom_index.
        // RDKit uses MatchVectType (vector<pair<int,int>>) where
        // pair.second is the mol atom index and pair.first is query atom index.
        let mut atom_to_mol: Vec<Option<usize>> = vec![None; q_num_atoms];
        for (&qa, &ma) in c1.iter().zip(c2.iter()) {
            if qa < q_num_atoms {
                atom_to_mol[qa] = Some(ma);
            }
        }

        // Build bond mapping by looking up bonds between matched atoms.
        let mut bond_mapping = Vec::with_capacity(query.num_bonds());
        for qbond in query.bonds() {
            let q_begin = qbond.begin().index();
            let q_end = qbond.end().index();
            let m_begin = atom_to_mol[q_begin];
            let m_end = atom_to_mol[q_end];
            match (m_begin, m_end) {
                (Some(mb), Some(me)) => {
                    // Find bond between mb and me in mol.
                    let found = m_graph.adjacency[mb]
                        .iter()
                        .find(|&&(nbr, _)| nbr == me)
                        .map(|&(_, eidx)| eidx);
                    bond_mapping.push(found.unwrap_or(NULL_NODE));
                }
                _ => {
                    bond_mapping.push(NULL_NODE);
                }
            }
        }

        results.push(SubstructMatchResult {
            atom_mapping: atom_to_mol
                .into_iter()
                .map(|x| x.unwrap_or(NULL_NODE))
                .collect(),
            bond_mapping,
        });
    }

    Ok(results)
}

fn atom_matches_with_recursive_cache(
    query_atom: &Atom,
    query_mol: &Molecule,
    mol_atom: &Atom,
    mol: &Molecule,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
    query_ctx: &QueryMatchContext,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructUtils.cpp :: atomCompat
    // RDKit✔️✔️:   if (ps.useQueryQueryMatches && a1->hasQuery() && a2->hasQuery()) {
    // RDKit✔️✔️:     res = static_cast<const QueryAtom *>(a1)->QueryMatch(
    // RDKit✔️✔️:         static_cast<const QueryAtom *>(a2));
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = a1->Match(a2);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!res) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    if let Some(query_node) = query_atom.query() {
        return match query_node {
            crate::QueryNode::Predicate(pred) => atom_query_predicate_matches_for_substruct(
                mol_atom,
                pred,
                mol,
                params,
                recursive_cache,
                query_ctx,
            ),
            crate::QueryNode::And(children) => children.iter().all(|child| match child {
                crate::QueryNode::Predicate(pred) => atom_query_predicate_matches_for_substruct(
                    mol_atom,
                    pred,
                    mol,
                    params,
                    recursive_cache,
                    query_ctx,
                ),
                _ => evaluate_atom_query(child, mol_atom, mol, params, recursive_cache, query_ctx),
            }),
            crate::QueryNode::Or(children) => children.iter().any(|child| match child {
                crate::QueryNode::Predicate(pred) => atom_query_predicate_matches_for_substruct(
                    mol_atom,
                    pred,
                    mol,
                    params,
                    recursive_cache,
                    query_ctx,
                ),
                _ => evaluate_atom_query(child, mol_atom, mol, params, recursive_cache, query_ctx),
            }),
            crate::QueryNode::Not(child) => {
                !evaluate_atom_query(child, mol_atom, mol, params, recursive_cache, query_ctx)
            }
        };
    }
    atom_matches(query_atom, query_mol, mol_atom, mol)
}

fn substruct_match_impl(
    mol: &Molecule,
    query: &Molecule,
    params: &SubstructMatchParams,
) -> SubstructMatchResultList {
    let trace_row64_substruct = env::var("RDKIT_ROW64_SUBSTRUCT_TRACE").ok().as_deref()
        == Some("1")
        && mol.num_atoms() == 106;
    let recursive_cache_start = trace_row64_substruct.then(Instant::now);
    let mut recursive_cache = RecursiveQueryMatchCache::new();
    populate_recursive_query_match_cache(mol, query, params, &mut recursive_cache)?;
    let recursive_cache_elapsed = recursive_cache_start.map(|start| start.elapsed().as_secs_f64());
    let match_start = trace_row64_substruct.then(Instant::now);
    let results =
        substruct_match_impl_with_recursive_cache(mol, query, params, Some(&recursive_cache))?;
    if let Some(match_start) = match_start {
        println!(
            "row64_substruct_timing q_atoms={} q_bonds={} recursive_cache={:.6} match_core={:.6} matches={}",
            query.num_atoms(),
            query.num_bonds(),
            recursive_cache_elapsed.unwrap_or(0.0),
            match_start.elapsed().as_secs_f64(),
            results.len()
        );
    }
    Ok(results)
}

/// Check if a molecule contains a substructure match for the given query.
///
/// This is the public API for `has_substruct_match`.
/// RDKit✔️❌: VF2-based substructure matching ported from vf2.hpp + SubstructMatch.cpp.
pub fn has_substruct_match(mol: &Molecule, query: &Molecule) -> bool {
    let params = SubstructMatchParams::default();
    let mut params = params;
    params.max_matches = 1;
    substruct_match_impl(mol, query, &params)
        .map(|matches| !matches.is_empty())
        .unwrap_or(false)
}

/// Get the first substructure match, if any.
///
/// This is the public API for `get_substruct_match`.
/// RDKit✔️❌: VF2-based substructure matching ported from vf2.hpp + SubstructMatch.cpp.
pub fn get_substruct_match(mol: &Molecule, query: &Molecule) -> Option<SubstructMatchResult> {
    let params = SubstructMatchParams::default();
    let mut params = params;
    params.max_matches = 1;
    substruct_match_impl(mol, query, &params)
        .ok()
        .and_then(|matches| matches.into_iter().next())
}

/// Get all substructure matches with default parameters.
///
/// This is the public API for `get_substruct_matches`.
/// RDKit✔️❌: VF2-based substructure matching ported from vf2.hpp + SubstructMatch.cpp.
pub fn get_substruct_matches(mol: &Molecule, query: &Molecule) -> Vec<SubstructMatchResult> {
    let params = SubstructMatchParams::default();
    substruct_match_impl(mol, query, &params).unwrap_or_default()
}

/// Get all substructure matches with custom parameters.
///
/// This is the public API for `get_substruct_matches_with_params`.
/// RDKit✔️❌: VF2-based substructure matching ported from vf2.hpp + SubstructMatch.cpp.
pub fn get_substruct_matches_with_params(
    mol: &Molecule,
    query: &Molecule,
    params: &SubstructMatchParams,
) -> Vec<SubstructMatchResult> {
    substruct_match_impl(mol, query, params).unwrap_or_default()
}

/// Get all substructure matches with custom parameters and structured
/// unsupported-feature errors for source-porting callers.
pub fn try_get_substruct_matches_with_params(
    mol: &Molecule,
    query: &Molecule,
    params: &SubstructMatchParams,
) -> SubstructMatchResultList {
    substruct_match_impl(mol, query, params)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MoleculeBuilder;

    fn make_mol_c() -> Molecule {
        // Methane: C
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder.build().expect("build methane")
    }

    fn make_mol_cc() -> Molecule {
        // Ethane: CC
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c0, c1, BondOrder::Single))
            .expect("add bond");
        builder.build().expect("build ethane")
    }

    fn make_mol_cco() -> Molecule {
        // Ethanol: CCO
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let o = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(c0, c1, BondOrder::Single))
            .expect("add C-C bond");
        builder
            .add_bond(crate::BondSpec::new(c1, o, BondOrder::Single))
            .expect("add C-O bond");
        builder.build().expect("build ethanol")
    }

    fn make_mol_coc() -> Molecule {
        // Dimethyl ether: COC
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let o = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c0, o, BondOrder::Single))
            .expect("add C-O bond 1");
        builder
            .add_bond(crate::BondSpec::new(o, c1, BondOrder::Single))
            .expect("add O-C bond 2");
        builder.build().expect("build dimethyl ether")
    }

    #[test]
    fn test_has_substruct_match_self() {
        let c = make_mol_c();
        assert!(
            has_substruct_match(&c, &c),
            "a molecule should match itself"
        );
    }

    #[test]
    fn test_has_substruct_match_cc_in_cco() {
        let cc = make_mol_cc();
        let cco = make_mol_cco();
        assert!(
            has_substruct_match(&cco, &cc),
            "CCO should contain CC as substructure"
        );
    }

    #[test]
    fn test_has_substruct_match_no_match() {
        let c = make_mol_c();
        let cco = make_mol_cco();
        assert!(
            !has_substruct_match(&c, &cco),
            "a single carbon should not contain CCO"
        );
    }

    #[test]
    fn test_get_substruct_match_self() {
        let cco = make_mol_cco();
        let result = get_substruct_match(&cco, &cco);
        assert!(result.is_some(), "self-match should return Some");
        let result = result.unwrap();
        assert_eq!(result.atom_mapping.len(), 3);
        // Identity mapping: 0->0, 1->1, 2->2
        for (qa, ma) in result.atom_mapping.iter().enumerate() {
            assert_eq!(*ma, qa, "self-match should have identity mapping");
        }
    }

    #[test]
    fn test_get_substruct_match_cc_in_cco() {
        let cc = make_mol_cc();
        let cco = make_mol_cco();
        let result = get_substruct_match(&cco, &cc);
        assert!(result.is_some(), "CC should match in CCO");
    }

    #[test]
    fn test_get_substruct_match_no_match() {
        let c = make_mol_c();
        let cco = make_mol_cco();
        let result = get_substruct_match(&c, &cco);
        assert!(
            result.is_none(),
            "C should not match CCO (query larger than mol)"
        );
    }

    #[test]
    fn test_get_substruct_matches_cco_in_cco() {
        let cco = make_mol_cco();
        let matches = get_substruct_matches(&cco, &cco);
        assert!(!matches.is_empty(), "should find at least self-match");
    }

    #[test]
    fn test_substruct_coc_matches_cco() {
        // COC (dimethyl ether) should not match CCO (ethanol) — different topology.
        let coc = make_mol_coc();
        let cco = make_mol_cco();
        assert!(
            !has_substruct_match(&cco, &coc),
            "CCO should not match COC topology"
        );
        // But CO should match CCO (CO is a substructure of CCO).
        let mut builder = MoleculeBuilder::new();
        let c = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let o = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(c, o, BondOrder::Single))
            .expect("add CO bond");
        let co = builder.build().expect("build CO");
        assert!(has_substruct_match(&cco, &co), "CCO should match CO");
    }

    #[test]
    fn test_has_substruct_match_empty_mol() {
        let empty = Molecule::new();
        let c = make_mol_c();
        assert!(
            !has_substruct_match(&empty, &c),
            "empty molecule should not match anything"
        );
        assert!(
            !has_substruct_match(&c, &empty),
            "molecule should not match empty query"
        );
    }

    #[test]
    fn test_substruct_match_params_max_matches() {
        let c = make_mol_c();
        let params = SubstructMatchParams {
            max_matches: 1,
            uniquify: true,
            use_chirality: false,
            specified_stereo_query_matches_unspecified: false,
        };
        let matches = get_substruct_matches_with_params(&c, &c, &params);
        assert_eq!(matches.len(), 1, "max_matches=1 should return one match");
    }

    #[test]
    fn feature_smarts_substruct_matches_required_query_semantics() {
        let cases = [
            (
                "Donor",
                "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]",
                "CCO",
                vec![2],
                vec![2],
            ),
            (
                "Acceptor",
                "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]",
                "CC(=O)C",
                vec![2],
                vec![2],
            ),
            (
                "Aromatic",
                "[a]",
                "c1ccccc1",
                vec![0],
                vec![0, 1, 2, 3, 4, 5],
            ),
            ("Halogen", "[F,Cl,Br,I]", "CCl", vec![1], vec![1]),
            (
                "Basic",
                "[#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",
                "[NH4+]",
                vec![0],
                vec![0],
            ),
            (
                "Acidic",
                "[$([C,S](=[O,S,P])-[O;H1,-1])]",
                "CC(=O)O",
                vec![1],
                vec![1],
            ),
        ];

        for (name, smarts, smiles, expected_first, expected_atoms) in cases {
            let mol = Molecule::from_smiles_with_sanitize(smiles, false)
                .unwrap_or_else(|_| panic!("{name} molecule should parse"));
            let query = build_crystalff_query_molecule(smarts)
                .unwrap_or_else(|_| panic!("{name} SMARTS should build query molecule"));
            let matches = get_substruct_matches(&mol, &query);
            assert!(
                !matches.is_empty(),
                "{name} should produce at least one match"
            );
            assert_eq!(
                matches[0].atom_mapping, expected_first,
                "{name} first match atom mapping"
            );
            let mut atom_indices: Vec<usize> = matches
                .iter()
                .flat_map(|matched| matched.atom_mapping.iter().copied())
                .filter(|idx| *idx != NULL_NODE)
                .collect();
            atom_indices.sort_unstable();
            atom_indices.dedup();
            assert_eq!(atom_indices, expected_atoms, "{name} feature SMARTS");
        }
    }

    #[test]
    fn lipinski_hba_recursive_smarts_matches_rdkit_root_semantics() {
        const HBA: &str = "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0X2,o,s;+0])]";
        let cases = [
            (
                "alcohol oxygen recursive branch",
                "[O,S;H1;v2]-[!$(*=[O,N,P,S])]",
                "CCO",
                vec![vec![2, 1]],
            ),
            (
                "carboxylic acid oxygen rejected by negated recursive neighbor",
                "[O,S;H1;v2]-[!$(*=[O,N,P,S])]",
                "CC(=O)O",
                Vec::<Vec<usize>>::new(),
            ),
            (
                "amine nitrogen recursive branch",
                "[N;v3;!$(N-*=!@[O,N,P,S])]",
                "CCN",
                vec![vec![2]],
            ),
            (
                "amide nitrogen rejected by mixed bond recursive query",
                "[N;v3;!$(N-*=!@[O,N,P,S])]",
                "CC(=O)N",
                Vec::<Vec<usize>>::new(),
            ),
            (
                "mixed single-double non-ring bond query",
                "N-*=!@[O,N,P,S]",
                "CC(=O)N",
                vec![vec![3, 1, 2]],
            ),
            ("full HBA ethanol", HBA, "CCO", vec![vec![2]]),
            ("full HBA carboxylic acid", HBA, "CC(=O)O", vec![vec![2]]),
            ("full HBA amide", HBA, "CC(=O)N", vec![vec![2]]),
            ("full HBA pyridine", HBA, "c1ccncc1", vec![vec![3]]),
            ("full HBA furan", HBA, "c1ccoc1", vec![vec![3]]),
        ];

        for (name, smarts, smiles, expected) in cases {
            let mol = Molecule::from_smiles_with_sanitize(smiles, true)
                .unwrap_or_else(|_| panic!("{name} molecule should parse"));
            let query = build_crystalff_query_molecule(smarts)
                .unwrap_or_else(|_| panic!("{name} SMARTS should build"));
            let matches = get_substruct_matches(&mol, &query);
            let atom_mappings = matches
                .iter()
                .map(|matched| matched.atom_mapping.clone())
                .collect::<Vec<_>>();
            assert_eq!(atom_mappings, expected, "{name}");
        }
    }

    #[test]
    fn maccs_patterns_substruct_matches_required_topology_semantics() {
        let cases = [
            (
                "four-membered ring",
                "*1~*~*~*~1",
                "C1CCC1",
                true,
                vec![0, 1, 2, 3],
            ),
            (
                "four-membered ring rejects chain",
                "*1~*~*~*~1",
                "CCCC",
                false,
                Vec::new(),
            ),
            ("ring bond", "*@*(@*)@*", "C12CC1C2", true, vec![0, 2, 1, 3]),
            (
                "non-ring oxygen bridge",
                "*!@[#8]!@*",
                "COC",
                true,
                vec![0, 1, 2],
            ),
            (
                "branch degree",
                "*~*(~*)(~*)~*",
                "CC(C)(C)C",
                true,
                vec![0, 1, 2, 3, 4],
            ),
            (
                "recursive ring closure",
                "[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]",
                "CC1CC1",
                true,
                vec![0],
            ),
        ];

        for (name, smarts, smiles, expected_match, expected_first) in cases {
            let mol = Molecule::from_smiles_with_sanitize(smiles, true)
                .unwrap_or_else(|_| panic!("{name} molecule should parse"));
            let query = build_crystalff_query_molecule(smarts)
                .unwrap_or_else(|_| panic!("{name} MACCS SMARTS should build query"));
            let matches = get_substruct_matches(&mol, &query);
            assert_eq!(
                !matches.is_empty(),
                expected_match,
                "{name} match truth value"
            );
            if expected_match {
                assert_eq!(
                    matches[0].atom_mapping, expected_first,
                    "{name} first match"
                );
            }
        }
    }

    struct MaccsPatternGolden {
        bit: u16,
        smarts: &'static str,
        smiles: &'static str,
        first_match: &'static [usize],
    }

    fn rdkit_maccs_pattern_positive_goldens() -> &'static [MaccsPatternGolden] {
        // RDKit source: MACCS.cpp::Patterns initializes these SMARTS strings
        // with `RDKit::SmartsToMol(...)`. `first_match` values were generated
        // from pinned RDKit 2026.03.1 using `Mol.GetSubstructMatch()`.
        &[
            MaccsPatternGolden {
                bit: 8,
                smarts: "[!#6!#1]1~*~*~*~1",
                smiles: "O1CCC1",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 11,
                smarts: "*1~*~*~*~1",
                smiles: "C1CCC1",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 13,
                smarts: "[#8]~[#7](~[#6])~[#6]",
                smiles: "ON(C)C",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 14,
                smarts: "[#16]-[#16]",
                smiles: "CSSC",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 15,
                smarts: "[#8]~[#6](~[#8])~[#8]",
                smiles: "O=C(O)O",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 16,
                smarts: "[!#6!#1]1~*~*~1",
                smiles: "O1CC1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 17,
                smarts: "[#6]#[#6]",
                smiles: "C#C",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 19,
                smarts: "*1~*~*~*~*~*~*~1",
                smiles: "C1CCCCCC1",
                first_match: &[0, 1, 2, 3, 4, 5, 6],
            },
            MaccsPatternGolden {
                bit: 20,
                smarts: "[#14]",
                smiles: "[SiH4]",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 21,
                smarts: "[#6]=[#6](~[!#6!#1])~[!#6!#1]",
                smiles: "C=C(O)O",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 22,
                smarts: "*1~*~*~1",
                smiles: "C1CC1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 23,
                smarts: "[#7]~[#6](~[#8])~[#8]",
                smiles: "NC(=O)O",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 24,
                smarts: "[#7]-[#8]",
                smiles: "ON(C)C",
                first_match: &[1, 0],
            },
            MaccsPatternGolden {
                bit: 25,
                smarts: "[#7]~[#6](~[#7])~[#7]",
                smiles: "NC(N)N",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 26,
                smarts: "[#6]=@[#6](@*)@*",
                smiles: "C1=C2CCCC2C1",
                first_match: &[0, 1, 2, 5],
            },
            MaccsPatternGolden {
                bit: 28,
                smarts: "[!#6!#1]~[CH2]~[!#6!#1]",
                smiles: "OCO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 30,
                smarts: "[#6]~[!#6!#1](~[#6])(~[#6])~*",
                smiles: "C[S](C)(C)C",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 31,
                smarts: "[!#6!#1]~[F,Cl,Br,I]",
                smiles: "N[Pt](Cl)(Cl)N",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 32,
                smarts: "[#6]~[#16]~[#7]",
                smiles: "CSN",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 33,
                smarts: "[#7]~[#16]",
                smiles: "CSN",
                first_match: &[2, 1],
            },
            MaccsPatternGolden {
                bit: 34,
                smarts: "[CH2]=*",
                smiles: "C=C",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 36,
                smarts: "[#16R]",
                smiles: "S1CC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 37,
                smarts: "[#7]~[#6](~[#8])~[#7]",
                smiles: "NC(=O)N",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 38,
                smarts: "[#7]~[#6](~[#6])~[#7]",
                smiles: "NC(C)N",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 39,
                smarts: "[#8]~[#16](~[#8])~[#8]",
                smiles: "COS(=O)(=O)O",
                first_match: &[1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 40,
                smarts: "[#16]-[#8]",
                smiles: "CSO",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 41,
                smarts: "[#6]#[#7]",
                smiles: "C#N",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 43,
                smarts: "[!#6!#1!H0]~*~[!#6!#1!H0]",
                smiles: "OCO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 44,
                smarts: "[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]",
                smiles: "[SeH2]",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 45,
                smarts: "[#6]=[#6]~[#7]",
                smiles: "C=CN",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 47,
                smarts: "[#16]~*~[#7]",
                smiles: "SCN",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 48,
                smarts: "[#8]~[!#6!#1](~[#8])~[#8]",
                smiles: "COS(=O)(=O)O",
                first_match: &[1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 49,
                smarts: "[!+0]",
                smiles: "O=N(=O)O",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 50,
                smarts: "[#6]=[#6](~[#6])~[#6]",
                smiles: "CC(C)=C",
                first_match: &[3, 1, 0, 2],
            },
            MaccsPatternGolden {
                bit: 51,
                smarts: "[#6]~[#16]~[#8]",
                smiles: "CSO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 52,
                smarts: "[#7]~[#7]",
                smiles: "NNO",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 53,
                smarts: "[!#6!#1!H0]~*~*~*~[!#6!#1!H0]",
                smiles: "NCCCO",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 54,
                smarts: "[!#6!#1!H0]~*~*~[!#6!#1!H0]",
                smiles: "NCCN",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 55,
                smarts: "[#8]~[#16]~[#8]",
                smiles: "CS(=O)(=O)C",
                first_match: &[2, 1, 3],
            },
            MaccsPatternGolden {
                bit: 56,
                smarts: "[#8]~[#7](~[#8])~[#6]",
                smiles: "ON(O)C",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 57,
                smarts: "[#8R]",
                smiles: "O1CC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 58,
                smarts: "[!#6!#1]~[#16]~[!#6!#1]",
                smiles: "CS(=O)(=O)C",
                first_match: &[2, 1, 3],
            },
            MaccsPatternGolden {
                bit: 59,
                smarts: "[#16]!:*:*",
                smiles: "Sc1ccccc1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 60,
                smarts: "[#16]=[#8]",
                smiles: "CS(=O)C",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 61,
                smarts: "*~[#16](~*)~*",
                smiles: "C[S](C)(C)C",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 62,
                smarts: "*@*!@*@*",
                smiles: "C1CC1C1CC1",
                first_match: &[0, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 63,
                smarts: "[#7]=[#8]",
                smiles: "N=O",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 64,
                smarts: "*@*!@[#16]",
                smiles: "Sc1ccccc1",
                first_match: &[2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 65,
                smarts: "c:n",
                smiles: "c1ncccc1",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 66,
                smarts: "[#6]~[#6](~[#6])(~[#6])~*",
                smiles: "CC(C)(C)C",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 67,
                smarts: "[!#6!#1]~[#16]",
                smiles: "CSSC",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 68,
                smarts: "[!#6!#1!H0]~[!#6!#1!H0]",
                smiles: "NNO",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 69,
                smarts: "[!#6!#1]~[!#6!#1!H0]",
                smiles: "CSN",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 70,
                smarts: "[!#6!#1]~[#7]~[!#6!#1]",
                smiles: "NNO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 71,
                smarts: "[#7]~[#8]",
                smiles: "N=O",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 72,
                smarts: "[#8]~*~*~[#8]",
                smiles: "OCCO",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 73,
                smarts: "[#16]=*",
                smiles: "CS(=O)C",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 74,
                smarts: "[CH3]~*~[CH3]",
                smiles: "CCC",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 75,
                smarts: "*!@[#7]@*",
                smiles: "CN1CC1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 76,
                smarts: "[#6]=[#6](~*)~*",
                smiles: "CC(C)=C",
                first_match: &[3, 1, 0, 2],
            },
            MaccsPatternGolden {
                bit: 77,
                smarts: "[#7]~*~[#7]",
                smiles: "NC(=O)N",
                first_match: &[0, 1, 3],
            },
            MaccsPatternGolden {
                bit: 78,
                smarts: "[#6]=[#7]",
                smiles: "C=N",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 79,
                smarts: "[#7]~*~*~[#7]",
                smiles: "NCCN",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 80,
                smarts: "[#7]~*~*~*~[#7]",
                smiles: "NCCCN",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 81,
                smarts: "[#16]~*(~*)~*",
                smiles: "Sc1ccccc1",
                first_match: &[0, 1, 2, 6],
            },
            MaccsPatternGolden {
                bit: 82,
                smarts: "*~[CH2]~[!#6!#1!H0]",
                smiles: "CCO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 83,
                smarts: "[!#6!#1]1~*~*~*~*~1",
                smiles: "O1CCCC1",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 84,
                smarts: "[NH2]",
                smiles: "CCN",
                first_match: &[2],
            },
            MaccsPatternGolden {
                bit: 85,
                smarts: "[#6]~[#7](~[#6])~[#6]",
                smiles: "CN(C)C",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 86,
                smarts: "[C;H2,H3][!#6!#1][C;H2,H3]",
                smiles: "COC",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 87,
                smarts: "[F,Cl,Br,I]!@*@*",
                smiles: "Clc1ccccc1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 89,
                smarts: "[#8]~*~*~*~[#8]",
                smiles: "OCCCO",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 90,
                smarts: "[$([!#6!#1!H0]~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[CH2R]1)]",
                smiles: "N1CCC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 91,
                smarts: "[$([!#6!#1!H0]~*~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~*~[R]1@[R]@[CH2R]1)]",
                smiles: "NCCCCCN",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 92,
                smarts: "[#8]~[#6](~[#7])~[#6]",
                smiles: "CC(=O)N",
                first_match: &[2, 1, 3, 0],
            },
            MaccsPatternGolden {
                bit: 93,
                smarts: "[!#6!#1]~[CH3]",
                smiles: "COC",
                first_match: &[1, 0],
            },
            MaccsPatternGolden {
                bit: 94,
                smarts: "[!#6!#1]~[#7]",
                smiles: "CSN",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 95,
                smarts: "[#7]~*~*~[#8]",
                smiles: "NCCO",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 96,
                smarts: "*1~*~*~*~*~1",
                smiles: "C1CCCC1",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 97,
                smarts: "[#7]~*~*~*~[#8]",
                smiles: "NCCCO",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 98,
                smarts: "[!#6!#1]1~*~*~*~*~*~1",
                smiles: "O1CCCCC1",
                first_match: &[0, 1, 2, 3, 4, 5],
            },
            MaccsPatternGolden {
                bit: 99,
                smarts: "[#6]=[#6]",
                smiles: "C=C",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 100,
                smarts: "*~[CH2]~[#7]",
                smiles: "CCN",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 101,
                smarts: "[$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1)]",
                smiles: "C1CCCCCCC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 102,
                smarts: "[!#6!#1]~[#8]",
                smiles: "CSO",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 104,
                smarts: "[!#6!#1!H0]~*~[CH2]~*",
                smiles: "CCCO",
                first_match: &[3, 2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 105,
                smarts: "*@*(@*)@*",
                smiles: "C12CC1C2",
                first_match: &[0, 2, 1, 3],
            },
            MaccsPatternGolden {
                bit: 106,
                smarts: "[!#6!#1]~*(~[!#6!#1])~[!#6!#1]",
                smiles: "COS(=O)(=O)O",
                first_match: &[1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 107,
                smarts: "[F,Cl,Br,I]~*(~*)~*",
                smiles: "CC(C)(C)Cl",
                first_match: &[4, 1, 0, 2],
            },
            MaccsPatternGolden {
                bit: 108,
                smarts: "[CH3]~*~*~*~[CH2]~*",
                smiles: "CCCCCC",
                first_match: &[0, 1, 2, 3, 4, 5],
            },
            MaccsPatternGolden {
                bit: 109,
                smarts: "*~[CH2]~[#8]",
                smiles: "CCO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 110,
                smarts: "[#7]~[#6]~[#8]",
                smiles: "CC(=O)N",
                first_match: &[3, 1, 2],
            },
            MaccsPatternGolden {
                bit: 111,
                smarts: "[#7]~*~[CH2]~*",
                smiles: "CCCN",
                first_match: &[3, 2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 112,
                smarts: "*~*(~*)(~*)~*",
                smiles: "CC(C)(C)C",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 113,
                smarts: "[#8]!:*:*",
                smiles: "Oc1ccccc1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 114,
                smarts: "[CH3]~[CH2]~*",
                smiles: "CCC",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 115,
                smarts: "[CH3]~*~[CH2]~*",
                smiles: "CCCC",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 116,
                smarts: "[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]",
                smiles: "CCCCC",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 117,
                smarts: "[#7]~*~[#8]",
                smiles: "CC(=O)N",
                first_match: &[3, 1, 2],
            },
            MaccsPatternGolden {
                bit: 118,
                smarts: "[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]",
                smiles: "CCCC",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 119,
                smarts: "[#7]=*",
                smiles: "N=O",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 120,
                smarts: "[!#6R]",
                smiles: "O1CC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 121,
                smarts: "[#7R]",
                smiles: "N1CC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 122,
                smarts: "*~[#7](~*)~*",
                smiles: "ON(C)C",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 123,
                smarts: "[#8]~[#6]~[#8]",
                smiles: "OCO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 124,
                smarts: "[!#6!#1]~[!#6!#1]",
                smiles: "CSSC",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 126,
                smarts: "*!@[#8]!@*",
                smiles: "COC",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 127,
                smarts: "*@*!@[#8]",
                smiles: "Oc1ccccc1",
                first_match: &[2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 128,
                smarts: "[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2R]@[R]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2R]1),$(*~[CH2]~*~[R]1@[R]@[CH2R]1)]",
                smiles: "CCCCCCC",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 129,
                smarts: "[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[CH2R]1)]",
                smiles: "CCCCCC",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 131,
                smarts: "[!#6!#1!H0]",
                smiles: "CCO",
                first_match: &[2],
            },
            MaccsPatternGolden {
                bit: 132,
                smarts: "[#8]~*~[CH2]~*",
                smiles: "CCCO",
                first_match: &[3, 2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 133,
                smarts: "*@*!@[#7]",
                smiles: "Nc1ccccc1",
                first_match: &[2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 135,
                smarts: "[#7]!:*:*",
                smiles: "Nc1ccccc1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 136,
                smarts: "[#8]=*",
                smiles: "CS(=O)C",
                first_match: &[2, 1],
            },
            MaccsPatternGolden {
                bit: 137,
                smarts: "[!C!cR]",
                smiles: "O1CC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 138,
                smarts: "[!#6!#1]~[CH2]~*",
                smiles: "CCO",
                first_match: &[2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 139,
                smarts: "[O!H0]",
                smiles: "CCO",
                first_match: &[2],
            },
            MaccsPatternGolden {
                bit: 140,
                smarts: "[#8]",
                smiles: "CCO",
                first_match: &[2],
            },
            MaccsPatternGolden {
                bit: 141,
                smarts: "[CH3]",
                smiles: "CC",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 142,
                smarts: "[#7]",
                smiles: "C#N",
                first_match: &[1],
            },
            MaccsPatternGolden {
                bit: 144,
                smarts: "*!:*:*!:*",
                smiles: "Cc1ccccc1C",
                first_match: &[0, 1, 6, 7],
            },
            MaccsPatternGolden {
                bit: 145,
                smarts: "*1~*~*~*~*~*~1",
                smiles: "C1CCCCC1",
                first_match: &[0, 1, 2, 3, 4, 5],
            },
            MaccsPatternGolden {
                bit: 147,
                smarts: "[$(*~[CH2]~[CH2]~*),$([R]1@[CH2R]@[CH2R]1)]",
                smiles: "CCCC",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 148,
                smarts: "*~[!#6!#1](~*)~*",
                smiles: "C[S](C)(C)C",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 149,
                smarts: "[C;H3,H4]",
                smiles: "C",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 150,
                smarts: "*!@*@*!@*",
                smiles: "Cc1ccccc1C",
                first_match: &[0, 1, 6, 7],
            },
            MaccsPatternGolden {
                bit: 151,
                smarts: "[#7!H0]",
                smiles: "CCN",
                first_match: &[2],
            },
            MaccsPatternGolden {
                bit: 152,
                smarts: "[#8]~[#6](~[#6])~[#6]",
                smiles: "CC(C)(C)O",
                first_match: &[4, 1, 0, 2],
            },
            MaccsPatternGolden {
                bit: 154,
                smarts: "[#6]=[#8]",
                smiles: "CC(=O)O",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 155,
                smarts: "*!@[CH2]!@*",
                smiles: "CCC",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 156,
                smarts: "[#7]~*(~*)~*",
                smiles: "CC(C)N",
                first_match: &[3, 1, 0, 2],
            },
            MaccsPatternGolden {
                bit: 157,
                smarts: "[#6]-[#8]",
                smiles: "CCO",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 158,
                smarts: "[#6]-[#7]",
                smiles: "CCN",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 162,
                smarts: "a",
                smiles: "c1ccccc1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 165,
                smarts: "[R]",
                smiles: "C1CC1",
                first_match: &[0],
            },
        ]
    }

    #[test]
    fn maccs_patterns_match_rdkit_positive_truth_and_first_atom_maps() {
        let goldens = rdkit_maccs_pattern_positive_goldens();
        assert_eq!(goldens.len(), 136);

        for golden in goldens {
            let mol =
                Molecule::from_smiles_with_sanitize(golden.smiles, true).unwrap_or_else(|error| {
                    panic!("MACCS bit {} target SMILES failed: {error}", golden.bit)
                });
            let query = build_crystalff_query_molecule(golden.smarts)
                .unwrap_or_else(|error| panic!("MACCS bit {} SMARTS failed: {error}", golden.bit));
            let matches = get_substruct_matches(&mol, &query);
            assert!(
                !matches.is_empty(),
                "MACCS bit {} should match RDKit positive target {} with SMARTS {}",
                golden.bit,
                golden.smiles,
                golden.smarts
            );
            assert_eq!(
                matches[0].atom_mapping, golden.first_match,
                "MACCS bit {} first RDKit atom map for {}",
                golden.bit, golden.smiles
            );
        }
    }

    #[test]
    fn maccs_bit_030_requires_four_explicit_neighbors_like_rdkit() {
        let mol = Molecule::from_smiles("ON(C)C").expect("fixture should parse");
        let query = build_crystalff_query_molecule("[#6]~[!#6!#1](~[#6])(~[#6])~*")
            .expect("MACCS bit 30 SMARTS should build");

        assert_eq!(query.num_atoms(), 5);
        assert_eq!(query.num_bonds(), 4);
        assert!(
            !has_substruct_match(&mol, &query),
            "RDKit does not match MACCS bit 30 against ON(C)C through the max-matches=1 path"
        );
        assert!(
            get_substruct_matches(&mol, &query).is_empty(),
            "RDKit does not match MACCS bit 30 against ON(C)C"
        );
    }

    #[test]
    fn test_atom_matches_basic() {
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c0, c1, BondOrder::Single))
            .expect("add bond");
        let mol = builder.build().expect("build");
        let q_atom = &mol.atoms()[0];
        let m_atom = &mol.atoms()[1];
        assert!(
            atom_matches(q_atom, &mol, m_atom, &mol),
            "two carbons should match"
        );
    }

    #[test]
    fn test_atom_matches_different_elements() {
        let mut builder = MoleculeBuilder::new();
        let c = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let o = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(c, o, BondOrder::Single))
            .expect("add bond");
        let mol = builder.build().expect("build");
        // Query atom = C, target atom = O — should not match (C != O, and
        // atomic number check should reject since 6 != 8).
        // But our atom_matches uses query atomic number: query=6, mol=8.
        // If query atomic number != 0, it must match. So C does not match O.
        let c_atom = &mol.atoms()[0]; // C
        let o_atom = &mol.atoms()[1]; // O
        assert!(
            !atom_matches(c_atom, &mol, o_atom, &mol),
            "C should not match O via basic atomic number check"
        );
    }

    #[test]
    fn test_bond_matches_single() {
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c0, c1, BondOrder::Single))
            .expect("add bond");
        let mol = builder.build().expect("build");
        assert!(bond_matches(&mol.bonds()[0], &mol, &mol.bonds()[0], &mol));
    }

    #[test]
    fn test_vf2_graph_building() {
        let cc = make_mol_cc();
        let g = build_vf2_graph(&cc);
        assert_eq!(g.n_atoms, 2);
        assert_eq!(g.n_bonds, 1);
        assert_eq!(g.out_degree(0), 1);
        assert_eq!(g.out_degree(1), 1);
        assert_eq!(g.out_edges(0)[0].0, 1);
        assert_eq!(g.out_edges(1)[0].0, 0);
    }

    #[test]
    fn test_sort_nodes_by_frequency_small() {
        let cc = make_mol_cc();
        let g = build_vf2_graph(&cc);
        let order = sort_nodes_by_frequency(&g);
        // Both nodes have degree 1, so order depends on sort stability.
        assert_eq!(order.len(), 2);
        // Both should be present.
        assert!(order.contains(&0));
        assert!(order.contains(&1));
    }

    #[test]
    fn test_vf2_state_initial() {
        let cc = make_mol_cc();
        let g = build_vf2_graph(&cc);
        let state = Vf2SubState::new(&g, &g, true);
        assert!(!state.is_goal());
        assert!(!state.is_dead());
        assert_eq!(state.core_len, 0);
    }

    #[test]
    fn test_vf2_next_pair_initial() {
        let cc = make_mol_cc();
        let g = build_vf2_graph(&cc);
        let state = Vf2SubState::new(&g, &g, true);
        let mut pair = Vf2Pair::new();
        let has_next = state.next_pair(&mut pair);
        assert!(has_next, "should find a next pair");
    }
}
