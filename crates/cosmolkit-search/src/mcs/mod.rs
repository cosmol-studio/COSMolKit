//! Maximum common substructure search over detached model values.

use crate::query_behavior::{
    CompositeQueryType, make_atom_in_ring_query, make_atom_isotope_query,
    make_atom_min_ring_size_query, make_atom_num_query, make_bond_is_in_ring_query,
    make_bond_order_equals_query, query_atom_expand_query, query_bond_expand_query,
};
use crate::{SearchTarget, SearchTargetAccess};
use cosmolkit_core::{RingFindType, RingInfo};
use cosmolkit_model::{
    AtomId, AtomQueryPredicate, Bond, BondId, BondSpec, QueryAtom, QueryAtomIdentity, QueryBond,
    QueryGraph, QueryNode,
};
use cosmolkit_types::{BondOrder, BondStereo, ChiralTag};
use std::collections::{BTreeMap, BTreeSet};

/// Built-in atom comparison modes used by FMCS.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AtomComparator {
    AtomCompareAny,
    #[default]
    AtomCompareElements,
    AtomCompareIsotopes,
    AtomCompareAnyHeavyAtom,
}

/// Built-in bond comparison modes used by FMCS.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BondComparator {
    BondCompareAny,
    #[default]
    BondCompareOrder,
    BondCompareOrderExact,
}

/// Ring-fusion comparison modes accepted by the source convenience entrypoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RingComparator {
    #[default]
    IgnoreRingFusion,
    PermissiveRingFusion,
    StrictRingFusion,
}

/// Atom-side FMCS comparison options.
#[derive(Debug, Clone, PartialEq)]
pub struct McsAtomCompareParameters {
    pub match_valences: bool,
    pub match_chiral_tag: bool,
    pub match_formal_charge: bool,
    pub ring_matches_ring_only: bool,
    pub complete_rings_only: bool,
    pub match_isotope: bool,
    pub max_distance: f64,
}

impl Default for McsAtomCompareParameters {
    fn default() -> Self {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.h :: MCSAtomCompareParameters member initializers
        // RDKit✔️✔️:   bool MatchValences = false;
        // RDKit✔️✔️:   bool MatchChiralTag = false;
        // RDKit✔️✔️:   bool MatchFormalCharge = false;
        // RDKit✔️✔️:   bool RingMatchesRingOnly = false;
        // RDKit✔️✔️:   bool CompleteRingsOnly = false;
        // RDKit✔️✔️:   bool MatchIsotope = false;
        // RDKit✔️✔️:   double MaxDistance = -1.0;
        // END RDKIT CPP FUNCTION
        // Local complexity review: fixed-size scalar initialization is O(1)
        // and allocates no heap state, matching the source member initializers.
        Self {
            match_valences: false,
            match_chiral_tag: false,
            match_formal_charge: false,
            ring_matches_ring_only: false,
            complete_rings_only: false,
            match_isotope: false,
            max_distance: -1.0,
        }
    }
}

/// Bond-side FMCS comparison options.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct McsBondCompareParameters {
    pub ring_matches_ring_only: bool,
    pub complete_rings_only: bool,
    pub match_fused_rings: bool,
    pub match_fused_rings_strict: bool,
    pub match_stereo: bool,
}

impl Default for McsBondCompareParameters {
    fn default() -> Self {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.h :: MCSBondCompareParameters member initializers
        // RDKit✔️✔️:   bool RingMatchesRingOnly = false;
        // RDKit✔️✔️:   bool CompleteRingsOnly = false;
        // RDKit✔️✔️:   bool MatchFusedRings = false;
        // RDKit✔️✔️:   bool MatchFusedRingsStrict = false;
        // RDKit✔️✔️:   bool MatchStereo = false;
        // END RDKIT CPP FUNCTION
        // Local complexity review: fixed-size scalar initialization is O(1)
        // and allocates no heap state, matching the source member initializers.
        Self {
            ring_matches_ring_only: false,
            complete_rings_only: false,
            match_fused_rings: false,
            match_fused_rings_strict: false,
            match_stereo: false,
        }
    }
}

/// Detached FMCS parameter value.
#[derive(Debug, Clone, PartialEq)]
pub struct McsParameters {
    pub store_all: bool,
    pub maximize_bonds: bool,
    pub threshold: f64,
    pub timeout: u32,
    pub verbose: bool,
    pub atom_compare_parameters: McsAtomCompareParameters,
    pub bond_compare_parameters: McsBondCompareParameters,
    pub atom_comparator: AtomComparator,
    pub bond_comparator: BondComparator,
    pub initial_seed: String,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum McsError {
    #[error("MCS requires at least two inputs, received {count}")]
    TooFewInputs { count: usize },
    #[error("MCS threshold must be at most 1.0")]
    ThresholdAboveOne,
    #[error("{side} MCS atom index {atom} is out of range")]
    AtomOutOfRange { side: &'static str, atom: usize },
    #[error("{side} MCS target has no 3D conformer")]
    MissingConformer { side: &'static str },
    #[error("{side} MCS conformer has no coordinate for atom {atom}")]
    CoordinateOutOfRange { side: &'static str, atom: usize },
    #[error("{side} MCS target has no ring information")]
    MissingRingInfo { side: &'static str },
    #[error("{side} MCS target has no valence assignment")]
    MissingValence { side: &'static str },
    #[error("{side} MCS valence assignment has no row for atom {atom}")]
    ValenceOutOfRange { side: &'static str, atom: usize },
    #[error("{side} MCS bond index {bond} is out of range")]
    BondOutOfRange { side: &'static str, bond: usize },
    #[error("MCS seed excluded-bond mask has {count} rows, so source bond {bond} is out of range")]
    SeedExcludedBondOutOfRange { bond: usize, count: usize },
    #[error("MCS seed source bond {bond} is already excluded")]
    SeedBondAlreadyExcluded { bond: usize },
    #[error("MCS seed has no mapped row for source atom {atom}")]
    SeedAtomMappingMissing { atom: usize },
    #[error("MCS target match has no target atom mapped for query atom {atom}")]
    TargetAtomMappingMissing { atom: usize },
    #[error("source atom {atom} is not an endpoint of MCS source bond {bond}")]
    SeedAtomNotBondEndpoint { atom: usize, bond: usize },
    #[error("MCS frontier source atom {atom} has no new-atom payload")]
    SeedNewAtomMissing { atom: usize },
    #[error("MCS remaining bond bound {remaining} is smaller than frontier size {frontier}")]
    SeedRemainingBondUnderflow { remaining: usize, frontier: usize },
    #[error("MCS remaining atom bound {remaining} is smaller than new-atom count {added}")]
    SeedRemainingAtomUnderflow { remaining: usize, added: usize },
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct McsInputOrder {
    molecule_indices: Vec<usize>,
    threshold_count: usize,
    start_index: usize,
    end_index: usize,
}

fn prepare_mcs_input_order(
    molecules: &[SearchTarget<'_>],
    threshold: f64,
) -> Result<McsInputOrder, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: MaximumCommonSubgraph::find input ordering
    // RDKit✔️✔️:   if (src_mols.size() < 2) {
    // RDKit✔️✔️:     throw std::runtime_error(
    // RDKit✔️✔️:         "FMCS. Invalid argument. mols.size() must be at least 2");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (Parameters.Threshold > 1.0) {
    // RDKit✔️✔️:     throw std::runtime_error(
    // RDKit✔️✔️:         "FMCS. Invalid argument. Parameter Threshold must be 1.0 or "
    // RDKit✔️✔️:         "less.");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // minimal required number of matched targets:
    // RDKit✔️✔️:   // at least one target, max all targets
    // RDKit✔️✔️:   ThresholdCount = static_cast<unsigned int>(std::min(
    // RDKit✔️✔️:       static_cast<int>(src_mols.size()) - 1,
    // RDKit✔️✔️:       std::max(1, static_cast<int>(ceil(static_cast<double>(src_mols.size()) *
    // RDKit✔️✔️:                                         Parameters.Threshold)) -
    // RDKit✔️✔️:                       1)));
    // RDKit✔️✔️:   // sort source set of molecules by their 'size' and assume the smallest
    // RDKit✔️✔️:   // molecule as a query
    // RDKit✔️✔️:   std::stable_sort(Molecules.begin(), Molecules.end(), molPtr_NumBondLess);
    // RDKit✔️✔️:   size_t startIdx = 0;
    // RDKit✔️✔️:   size_t endIdx = Molecules.size() - ThresholdCount;
    // RDKit✔️✔️:   while (startIdx < endIdx && !Molecules.at(startIdx)->getNumAtoms()) {
    // RDKit✔️✔️:     ++startIdx;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    // Local complexity review: both implementations retain one O(n) pointer/
    // index vector, use stable O(n log n) bond-count sorting, and scan only the
    // leading atom-empty prefix. No molecule or topology is cloned.
    if molecules.len() < 2 {
        return Err(McsError::TooFewInputs {
            count: molecules.len(),
        });
    }
    if threshold > 1.0 {
        return Err(McsError::ThresholdAboveOne);
    }
    let threshold_count = (molecules.len() - 1)
        .min(1_isize.max((molecules.len() as f64 * threshold).ceil() as isize - 1) as usize);
    let mut molecule_indices = (0..molecules.len()).collect::<Vec<_>>();
    molecule_indices.sort_by_key(|&index| molecules[index].bonds().len());
    let mut start_index = 0;
    let end_index = molecules.len() - threshold_count;
    while start_index < end_index && molecules[molecule_indices[start_index]].atoms().is_empty() {
        start_index += 1;
    }
    Ok(McsInputOrder {
        molecule_indices,
        threshold_count,
        start_index,
        end_index,
    })
}

fn mcs_atom<'a>(
    target: &'a SearchTarget<'_>,
    atom: usize,
    side: &'static str,
) -> Result<&'a cosmolkit_model::Atom, McsError> {
    target
        .atoms()
        .get(atom)
        .ok_or(McsError::AtomOutOfRange { side, atom })
}

fn mcs_bond<'a>(
    target: &'a SearchTarget<'_>,
    bond: usize,
    side: &'static str,
) -> Result<&'a cosmolkit_model::Bond, McsError> {
    target
        .bonds()
        .get(bond)
        .ok_or(McsError::BondOutOfRange { side, bond })
}

fn check_atom_chirality(
    left: &SearchTarget<'_>,
    left_atom: usize,
    right: &SearchTarget<'_>,
    right_atom: usize,
) -> Result<bool, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: checkAtomChirality
    // RDKit✔️✔️:   const auto a1 = mol1.getAtomWithIdx(atom1);
    // RDKit✔️✔️:   const auto a2 = mol2.getAtomWithIdx(atom2);
    // RDKit✔️✔️:   const auto ac1 = a1->getChiralTag();
    // RDKit✔️✔️:   const auto ac2 = a2->getChiralTag();
    // RDKit✔️✔️:   if (ac1 == Atom::CHI_TETRAHEDRAL_CW || ac1 == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:     return (ac2 == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:             ac2 == Atom::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // END RDKIT CPP FUNCTION
    // Local complexity review: two indexed atom reads and constant-time tag
    // comparisons reproduce the source O(1) helper without allocation.
    let left_tag = mcs_atom(left, left_atom, "left")?.chiral_tag();
    let right_tag = mcs_atom(right, right_atom, "right")?.chiral_tag();
    if matches!(
        left_tag,
        ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
    ) {
        Ok(matches!(
            right_tag,
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ))
    } else {
        Ok(true)
    }
}

fn mcs_coordinate(
    target: &SearchTarget<'_>,
    atom: usize,
    side: &'static str,
) -> Result<[f64; 3], McsError> {
    let conformer = target
        .conformers_3d()
        .first()
        .ok_or(McsError::MissingConformer { side })?;
    conformer
        .coordinates()
        .get(atom)
        .copied()
        .ok_or(McsError::CoordinateOutOfRange { side, atom })
}

fn check_atom_distance(
    params: &McsAtomCompareParameters,
    left: &SearchTarget<'_>,
    left_atom: usize,
    right: &SearchTarget<'_>,
    right_atom: usize,
) -> Result<bool, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: checkAtomDistance
    // RDKit✔️✔️:   const auto &ci1 = mol1.getConformer();
    // RDKit✔️✔️:   const auto &ci2 = mol2.getConformer();
    // RDKit✔️✔️:   const auto &pos1 = ci1.getAtomPos(atom1);
    // RDKit✔️✔️:   const auto &pos2 = ci2.getAtomPos(atom2);
    // RDKit✔️✔️:   bool withinRange = (pos1 - pos2).length() <= p.MaxDistance;
    // RDKit✔️✔️:   return withinRange;
    // END RDKIT CPP FUNCTION
    // Local complexity review: two direct coordinate reads and three scalar
    // differences reproduce the source O(1) distance calculation.
    let left_position = mcs_coordinate(left, left_atom, "left")?;
    let right_position = mcs_coordinate(right, right_atom, "right")?;
    let distance = left_position
        .iter()
        .zip(right_position)
        .map(|(left, right)| {
            let delta = left - right;
            delta * delta
        })
        .sum::<f64>()
        .sqrt();
    Ok(distance <= params.max_distance)
}

fn check_atom_stereo_and_distance(
    params: &McsAtomCompareParameters,
    left: &SearchTarget<'_>,
    left_atom: usize,
    right: &SearchTarget<'_>,
    right_atom: usize,
) -> Result<bool, McsError> {
    // RDKit✔️✔️:   if (p.MatchChiralTag && !checkAtomChirality(p, mol1, atom1, mol2, atom2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p.MaxDistance > 0 && !checkAtomDistance(p, mol1, atom1, mol2, atom2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // Complexity review: the option gates preserve the source constant-time
    // short circuit and avoid coordinate access when distance is disabled.
    if params.match_chiral_tag && !check_atom_chirality(left, left_atom, right, right_atom)? {
        return Ok(false);
    }
    if params.max_distance > 0.0
        && !check_atom_distance(params, left, left_atom, right, right_atom)?
    {
        return Ok(false);
    }
    Ok(true)
}

fn check_atom_charge(
    left: &SearchTarget<'_>,
    left_atom: usize,
    right: &SearchTarget<'_>,
    right_atom: usize,
) -> Result<bool, McsError> {
    // RDKit✔️✔️:   const auto a1 = mol1.getAtomWithIdx(atom1);
    // RDKit✔️✔️:   const auto a2 = mol2.getAtomWithIdx(atom2);
    // RDKit✔️✔️:   return a1->getFormalCharge() == a2->getFormalCharge();
    // Complexity review: two direct carrier reads reproduce the source O(1)
    // comparison without consulting query predicates.
    Ok(mcs_atom(left, left_atom, "left")?.formal_charge()
        == mcs_atom(right, right_atom, "right")?.formal_charge())
}

fn check_atom_ring_match(
    params: &McsAtomCompareParameters,
    left: &SearchTarget<'_>,
    left_atom: usize,
    right: &SearchTarget<'_>,
    right_atom: usize,
) -> Result<bool, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: checkAtomRingMatch
    // RDKit✔️✔️:   if (p.RingMatchesRingOnly) {
    // RDKit✔️✔️:     const auto ri1 = mol1.getRingInfo();
    // RDKit✔️✔️:     const auto ri2 = mol2.getRingInfo();
    // RDKit✔️✔️:     bool atom1inRing = (ri1->numAtomRings(atom1) > 0);
    // RDKit✔️✔️:     bool atom2inRing = (ri2->numAtomRings(atom2) > 0);
    // RDKit✔️✔️:     return atom1inRing == atom2inRing;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    // Local complexity review: ring assignments provide the same O(1) row
    // membership counts as source RingInfo and require no graph traversal.
    if !params.ring_matches_ring_only {
        return Ok(true);
    }
    let left_atom = mcs_atom(left, left_atom, "left")?;
    let right_atom = mcs_atom(right, right_atom, "right")?;
    let left_rings = left
        .ring_info()
        .ok_or(McsError::MissingRingInfo { side: "left" })?;
    let right_rings = right
        .ring_info()
        .ok_or(McsError::MissingRingInfo { side: "right" })?;
    Ok((left_rings.num_atom_rings(left_atom.id()) > 0)
        == (right_rings.num_atom_rings(right_atom.id()) > 0))
}

fn mcs_atom_compare_any(
    params: &McsAtomCompareParameters,
    left: &SearchTarget<'_>,
    left_atom: usize,
    right: &SearchTarget<'_>,
    right_atom: usize,
) -> Result<bool, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: MCSAtomCompareAny
    // RDKit✔️✔️:   if (p.MatchChiralTag && !checkAtomChirality(p, mol1, atom1, mol2, atom2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p.MatchFormalCharge && !checkAtomCharge(p, mol1, atom1, mol2, atom2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p.MaxDistance > 0 && !checkAtomDistance(p, mol1, atom1, mol2, atom2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p.RingMatchesRingOnly) {
    // RDKit✔️✔️:     return checkAtomRingMatch(p, mol1, atom1, mol2, atom2);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // END RDKIT CPP FUNCTION
    // Local complexity review: source short-circuit order is retained and
    // every enabled helper is O(1), so the comparator remains O(1).
    if params.match_chiral_tag && !check_atom_chirality(left, left_atom, right, right_atom)? {
        return Ok(false);
    }
    if params.match_formal_charge && !check_atom_charge(left, left_atom, right, right_atom)? {
        return Ok(false);
    }
    if params.max_distance > 0.0
        && !check_atom_distance(params, left, left_atom, right, right_atom)?
    {
        return Ok(false);
    }
    check_atom_ring_match(params, left, left_atom, right, right_atom)
}

fn mcs_atom_compare_any_heavy(
    params: &McsAtomCompareParameters,
    left: &SearchTarget<'_>,
    left_atom: usize,
    right: &SearchTarget<'_>,
    right_atom: usize,
) -> Result<bool, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: MCSAtomCompareAnyHeavyAtom
    // RDKit✔️✔️:   const auto a1 = mol1.getAtomWithIdx(atom1);
    // RDKit✔️✔️:   const auto a2 = mol2.getAtomWithIdx(atom2);
    // RDKit✔️✔️:   // Any atom, including H, matches another atom of the same type,  according to
    // RDKit✔️✔️:   // the other flags
    // RDKit✔️✔️:   if (a1->getAtomicNum() == a2->getAtomicNum() ||
    // RDKit✔️✔️:       (a1->getAtomicNum() > 1 && a2->getAtomicNum() > 1)) {
    // RDKit✔️✔️:     return MCSAtomCompareAny(p, mol1, atom1, mol2, atom2, nullptr);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // END RDKIT CPP FUNCTION
    // Local complexity review: two carrier atomic-number reads precede the
    // O(1) any-atom comparator, matching source branching and allocation cost.
    let left_number = mcs_atom(left, left_atom, "left")?.atomic_number();
    let right_number = mcs_atom(right, right_atom, "right")?.atomic_number();
    if left_number == right_number || (left_number > 1 && right_number > 1) {
        mcs_atom_compare_any(params, left, left_atom, right, right_atom)
    } else {
        Ok(false)
    }
}

fn mcs_total_valence(
    target: &SearchTarget<'_>,
    atom: usize,
    side: &'static str,
) -> Result<i32, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atom.cpp :: Atom::getTotalValence
    // RDKit✔️✔️: unsigned int Atom::getTotalValence() const {
    // RDKit✔️✔️:   return getValence(ValenceType::EXPLICIT) + getValence(ValenceType::IMPLICIT);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    // Local complexity review: both implementations read two indexed cached
    // scalar values and add them in O(1), without traversal or allocation.
    mcs_atom(target, atom, side)?;
    let valence = target.valence().ok_or(McsError::MissingValence { side })?;
    valence
        .explicit_valence
        .get(atom)
        .zip(valence.implicit_hydrogens.get(atom))
        .map(|(explicit, implicit)| explicit + implicit)
        .ok_or(McsError::ValenceOutOfRange { side, atom })
}

fn mcs_atom_compare_elements(
    params: &McsAtomCompareParameters,
    left: &SearchTarget<'_>,
    left_atom: usize,
    right: &SearchTarget<'_>,
    right_atom: usize,
) -> Result<bool, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: MCSAtomCompareElements
    // RDKit✔️✔️:   const auto a1 = mol1.getAtomWithIdx(atom1);
    // RDKit✔️✔️:   const auto a2 = mol2.getAtomWithIdx(atom2);
    // RDKit✔️✔️:   if (a1->getAtomicNum() != a2->getAtomicNum()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p.MatchValences && a1->getTotalValence() != a2->getTotalValence()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p.MatchChiralTag && !checkAtomChirality(p, mol1, atom1, mol2, atom2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p.MatchFormalCharge && !checkAtomCharge(p, mol1, atom1, mol2, atom2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p.MaxDistance > 0 && !checkAtomDistance(p, mol1, atom1, mol2, atom2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p.RingMatchesRingOnly) {
    // RDKit✔️✔️:     return checkAtomRingMatch(p, mol1, atom1, mol2, atom2);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // END RDKIT CPP FUNCTION
    // Local complexity review: the comparator performs only direct indexed
    // reads and source-ordered O(1) helper calls, with the same short circuits
    // and no allocation.
    let left_number = mcs_atom(left, left_atom, "left")?.atomic_number();
    let right_number = mcs_atom(right, right_atom, "right")?.atomic_number();
    if left_number != right_number {
        return Ok(false);
    }
    if params.match_valences
        && mcs_total_valence(left, left_atom, "left")?
            != mcs_total_valence(right, right_atom, "right")?
    {
        return Ok(false);
    }
    if params.match_chiral_tag && !check_atom_chirality(left, left_atom, right, right_atom)? {
        return Ok(false);
    }
    if params.match_formal_charge && !check_atom_charge(left, left_atom, right, right_atom)? {
        return Ok(false);
    }
    if params.max_distance > 0.0
        && !check_atom_distance(params, left, left_atom, right, right_atom)?
    {
        return Ok(false);
    }
    check_atom_ring_match(params, left, left_atom, right, right_atom)
}

fn mcs_atom_compare_isotopes(
    params: &McsAtomCompareParameters,
    left: &SearchTarget<'_>,
    left_atom: usize,
    right: &SearchTarget<'_>,
    right_atom: usize,
) -> Result<bool, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: MCSAtomCompareIsotopes
    // RDKit✔️✔️:   // ignore everything except isotope information:
    // RDKit✔️✔️:   const auto a1 = mol1.getAtomWithIdx(atom1);
    // RDKit✔️✔️:   const auto a2 = mol2.getAtomWithIdx(atom2);
    // RDKit✔️✔️:   if (a1->getIsotope() != a2->getIsotope()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p.MatchChiralTag && !checkAtomChirality(p, mol1, atom1, mol2, atom2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p.MatchFormalCharge && !checkAtomCharge(p, mol1, atom1, mol2, atom2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p.MaxDistance > 0 && !checkAtomDistance(p, mol1, atom1, mol2, atom2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p.RingMatchesRingOnly) {
    // RDKit✔️✔️:     return checkAtomRingMatch(p, mol1, atom1, mol2, atom2);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // END RDKIT CPP FUNCTION
    // Local complexity review: direct isotope and optional carrier/helper
    // reads preserve source O(1) short-circuit behavior without allocation.
    // RDKit's unset isotope carrier reads as zero, so `None` and explicit
    // zero share the same source identity here.
    let left_isotope = mcs_atom(left, left_atom, "left")?.isotope().unwrap_or(0);
    let right_isotope = mcs_atom(right, right_atom, "right")?.isotope().unwrap_or(0);
    if left_isotope != right_isotope {
        return Ok(false);
    }
    if params.match_chiral_tag && !check_atom_chirality(left, left_atom, right, right_atom)? {
        return Ok(false);
    }
    if params.match_formal_charge && !check_atom_charge(left, left_atom, right, right_atom)? {
        return Ok(false);
    }
    if params.max_distance > 0.0
        && !check_atom_distance(params, left, left_atom, right, right_atom)?
    {
        return Ok(false);
    }
    check_atom_ring_match(params, left, left_atom, right, right_atom)
}

fn check_bond_stereo(
    _params: &McsBondCompareParameters,
    left: &SearchTarget<'_>,
    left_bond: usize,
    right: &SearchTarget<'_>,
    right_bond: usize,
) -> Result<bool, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: checkBondStereo
    // RDKit✔️✔️:   const auto b1 = mol1.getBondWithIdx(bond1);
    // RDKit✔️✔️:   const auto b2 = mol2.getBondWithIdx(bond2);
    // RDKit✔️✔️:   auto bs1 = b1->getStereo();
    // RDKit✔️✔️:   auto bs2 = b2->getStereo();
    // RDKit✔️✔️:   if (b1->getBondType() == Bond::DOUBLE && b2->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:     if (bs1 > Bond::STEREOANY && !(bs2 > Bond::STEREOANY)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return bs1 == bs2;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // END RDKIT CPP FUNCTION
    // Local complexity review: two indexed bond reads and fixed enum/order
    // comparisons reproduce the source O(1) helper without allocation.
    // The source parameter is unused; MatchStereo is applied by each caller.
    let left = mcs_bond(left, left_bond, "left")?;
    let right = mcs_bond(right, right_bond, "right")?;
    let left_stereo = left.stereo();
    let right_stereo = right.stereo();
    if left.order() == BondOrder::Double && right.order() == BondOrder::Double {
        if left_stereo.rdkit_code() > BondStereo::Any.rdkit_code()
            && right_stereo.rdkit_code() <= BondStereo::Any.rdkit_code()
        {
            Ok(false)
        } else {
            Ok(left_stereo == right_stereo)
        }
    } else {
        Ok(true)
    }
}

fn mcs_ring_is_fused(ring_info: &cosmolkit_core::RingInfo, ring: usize) -> bool {
    ring_info.bond_rings()[ring]
        .iter()
        .any(|bond| ring_info.num_bond_rings(*bond) > 1)
}

fn have_pair_of_compatible_rings(
    _params: &McsBondCompareParameters,
    left: &SearchTarget<'_>,
    left_bond: usize,
    right: &SearchTarget<'_>,
    right_bond: usize,
) -> Result<bool, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: havePairOfCompatibleRings
    // RDKit✔️❌:   const auto ri1 = mol1.getRingInfo();
    // RDKit✔️❌:   const auto ri2 = mol2.getRingInfo();
    // RDKit✔️❌:   const auto &bondRings1 = ri1->bondRings();
    // RDKit✔️❌:   const auto &bondRings2 = ri2->bondRings();
    // RDKit✔️❌:   for (unsigned int ringIdx1 : ri1->bondMembers(bond1)) {
    // RDKit✔️❌:     const auto &ring1 = bondRings1.at(ringIdx1);
    // RDKit✔️❌:     bool isRing1Fused = ri1->isRingFused(ringIdx1);
    // RDKit✔️❌:     for (unsigned int ringIdx2 : ri2->bondMembers(bond2)) {
    // RDKit✔️❌:       const auto &ring2 = bondRings2.at(ringIdx2);
    // RDKit✔️❌:       if (ring1.size() == ring2.size()) {
    // RDKit✔️❌:         return true;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (isRing1Fused && ring2.size() > ring1.size()) {
    // RDKit✔️❌:         return true;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       bool isRing2Fused = ri2->isRingFused(ringIdx2);
    // RDKit✔️❌:       if (isRing2Fused && ring1.size() > ring2.size()) {
    // RDKit✔️❌:         return true;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return false;
    // END RDKIT CPP FUNCTION
    // Local complexity review: immutable detached RingInfo exposes exact ring
    // memberships but not its lazy fused-ring cache. Fusion is therefore
    // recomputed from shared bond membership for each inspected ring, which
    // preserves behavior but can add a ring-size scan relative to cached source
    // calls.
    let left_bond = mcs_bond(left, left_bond, "left")?;
    let right_bond = mcs_bond(right, right_bond, "right")?;
    let left_ring_info = left
        .ring_info()
        .ok_or(McsError::MissingRingInfo { side: "left" })?;
    let right_ring_info = right
        .ring_info()
        .ok_or(McsError::MissingRingInfo { side: "right" })?;
    let left_bond_rings = left_ring_info.bond_members(left_bond.id());
    let right_bond_rings = right_ring_info.bond_members(right_bond.id());
    for &left_ring in left_bond_rings {
        let left_size = left_ring_info.bond_rings()[left_ring].len();
        let left_is_fused = mcs_ring_is_fused(left_ring_info, left_ring);
        for &right_ring in right_bond_rings {
            let right_size = right_ring_info.bond_rings()[right_ring].len();
            if left_size == right_size {
                return Ok(true);
            }
            if left_is_fused && right_size > left_size {
                return Ok(true);
            }
            let right_is_fused = mcs_ring_is_fused(right_ring_info, right_ring);
            if right_is_fused && left_size > right_size {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn check_bond_ring_match(
    params: &McsBondCompareParameters,
    left: &SearchTarget<'_>,
    left_bond: usize,
    right: &SearchTarget<'_>,
    right_bond: usize,
) -> Result<bool, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: checkBondRingMatch
    // RDKit✔️❌:   const auto ri1 = mol1.getRingInfo();
    // RDKit✔️❌:   const auto ri2 = mol2.getRingInfo();
    // RDKit✔️❌:   // indices of rings in the query molecule
    // RDKit✔️❌:   const auto &ringIndices1 = ri1->bondMembers(bond1);
    // RDKit✔️❌:   // indices of rings in the target molecule
    // RDKit✔️❌:   const auto &ringIndices2 = ri2->bondMembers(bond2);
    // RDKit✔️❌:   bool bond1inRing = !ringIndices1.empty();
    // RDKit✔️❌:   bool bond2inRing = !ringIndices2.empty();
    // RDKit✔️❌:   bool res = (bond1inRing == bond2inRing);
    // RDKit✔️❌:   // if rings should be complete, we need to check upfront that there
    // RDKit✔️❌:   // is at least one pair of compatible rings; if there isn't, there
    // RDKit✔️❌:   // will never be a chance of complete match, so we should fail early
    // RDKit✔️❌:   if (p.CompleteRingsOnly && bond1inRing && bond2inRing) {
    // RDKit✔️❌:     res = havePairOfCompatibleRings(p, mol1, bond1, mol2, bond2);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   // bond are both either in a ring or not
    // RDKit✔️❌:   return res;
    // END RDKIT CPP FUNCTION
    // Local complexity review: direct membership checks remain O(1). The
    // CompleteRingsOnly branch delegates to the behavior-exact Q111 helper,
    // whose immutable fusion lookup has the documented extra ring-size scan.
    let left_bond_ref = mcs_bond(left, left_bond, "left")?;
    let right_bond_ref = mcs_bond(right, right_bond, "right")?;
    let left_ring_info = left
        .ring_info()
        .ok_or(McsError::MissingRingInfo { side: "left" })?;
    let right_ring_info = right
        .ring_info()
        .ok_or(McsError::MissingRingInfo { side: "right" })?;
    let left_in_ring = !left_ring_info.bond_members(left_bond_ref.id()).is_empty();
    let right_in_ring = !right_ring_info.bond_members(right_bond_ref.id()).is_empty();
    if params.complete_rings_only && left_in_ring && right_in_ring {
        have_pair_of_compatible_rings(params, left, left_bond, right, right_bond)
    } else {
        Ok(left_in_ring == right_in_ring)
    }
}

fn mcs_bond_compare_any(
    params: &McsBondCompareParameters,
    left: &SearchTarget<'_>,
    left_bond: usize,
    right: &SearchTarget<'_>,
    right_bond: usize,
) -> Result<bool, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: MCSBondCompareAny
    // RDKit✔️✔️:   if (p.MatchStereo && !checkBondStereo(p, mol1, bond1, mol2, bond2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️❌:   if (p.RingMatchesRingOnly) {
    // RDKit✔️❌:     return checkBondRingMatch(p, mol1, bond1, mol2, bond2);
    // RDKit✔️❌:   }
    // RDKit✔️✔️:   return true;
    // END RDKIT CPP FUNCTION
    // Local complexity review: stereo comparison and option short circuiting
    // remain O(1). Ring comparison inherits Q111's documented extra ring scan
    // only when complete-ring compatibility reaches fused-ring inspection.
    if params.match_stereo && !check_bond_stereo(params, left, left_bond, right, right_bond)? {
        return Ok(false);
    }
    if params.ring_matches_ring_only {
        return check_bond_ring_match(params, left, left_bond, right, right_bond);
    }
    Ok(true)
}

fn bond_orders_match(left: BondOrder, right: BondOrder, ignore_aromatization: bool) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: BondMatchOrderMatrix::BondMatchOrderMatrix/isEqual
    // RDKit✔️✔️:     memset(MatchMatrix, 0, sizeof(MatchMatrix));
    // RDKit✔️✔️:     // fill cells of the same and unspecified type
    // RDKit✔️✔️:     for (size_t i = 0; i <= Bond::ZERO; ++i) {
    // RDKit✔️✔️:       MatchMatrix[i][i] = true;
    // RDKit✔️✔️:       MatchMatrix[Bond::UNSPECIFIED][i] = MatchMatrix[i][Bond::UNSPECIFIED] =
    // RDKit✔️✔️:           true;
    // RDKit✔️✔️:       MatchMatrix[Bond::ZERO][i] = MatchMatrix[i][Bond::ZERO] = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (ignoreAromatization) {
    // RDKit✔️✔️:       MatchMatrix[Bond::SINGLE][Bond::AROMATIC] =
    // RDKit✔️✔️:           MatchMatrix[Bond::AROMATIC][Bond::SINGLE] = true;
    // RDKit✔️✔️:       MatchMatrix[Bond::SINGLE][Bond::ONEANDAHALF] =
    // RDKit✔️✔️:           MatchMatrix[Bond::ONEANDAHALF][Bond::SINGLE] = true;
    // RDKit✔️✔️:       MatchMatrix[Bond::DOUBLE][Bond::TWOANDAHALF] =
    // RDKit✔️✔️:           MatchMatrix[Bond::TWOANDAHALF][Bond::DOUBLE] = true;
    // RDKit✔️✔️:       MatchMatrix[Bond::TRIPLE][Bond::THREEANDAHALF] =
    // RDKit✔️✔️:           MatchMatrix[Bond::THREEANDAHALF][Bond::TRIPLE] = true;
    // RDKit✔️✔️:       MatchMatrix[Bond::QUADRUPLE][Bond::FOURANDAHALF] =
    // RDKit✔️✔️:           MatchMatrix[Bond::FOURANDAHALF][Bond::QUADRUPLE] = true;
    // RDKit✔️✔️:       MatchMatrix[Bond::QUINTUPLE][Bond::FIVEANDAHALF] =
    // RDKit✔️✔️:           MatchMatrix[Bond::FIVEANDAHALF][Bond::QUINTUPLE] = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   inline bool isEqual(unsigned int i, unsigned int j) const {
    // RDKit✔️✔️:     return MatchMatrix[i][j];
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    // Local complexity review: the model preserves the full closed RDKit bond
    // enum. A fixed set of scalar comparisons reproduces the static matrix's
    // O(1) lookup without allocation or data-dependent traversal.
    if left == right
        || matches!(left, BondOrder::Unspecified | BondOrder::Zero)
        || matches!(right, BondOrder::Unspecified | BondOrder::Zero)
    {
        return true;
    }
    ignore_aromatization
        && matches!(
            (left, right),
            (BondOrder::Single, BondOrder::Aromatic)
                | (BondOrder::Aromatic, BondOrder::Single)
                | (BondOrder::Single, BondOrder::OneAndHalf)
                | (BondOrder::OneAndHalf, BondOrder::Single)
                | (BondOrder::Double, BondOrder::TwoAndHalf)
                | (BondOrder::TwoAndHalf, BondOrder::Double)
                | (BondOrder::Triple, BondOrder::ThreeAndHalf)
                | (BondOrder::ThreeAndHalf, BondOrder::Triple)
                | (BondOrder::Quadruple, BondOrder::FourAndHalf)
                | (BondOrder::FourAndHalf, BondOrder::Quadruple)
                | (BondOrder::Quintuple, BondOrder::FiveAndHalf)
                | (BondOrder::FiveAndHalf, BondOrder::Quintuple)
        )
}

fn mcs_bond_compare_order(
    params: &McsBondCompareParameters,
    left: &SearchTarget<'_>,
    left_bond: usize,
    right: &SearchTarget<'_>,
    right_bond: usize,
) -> Result<bool, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: MCSBondCompareOrder
    // RDKit✔️✔️:   static const BondMatchOrderMatrix match(true);  // ignore Aromatization
    // RDKit✔️✔️:   const auto b1 = mol1.getBondWithIdx(bond1);
    // RDKit✔️✔️:   const auto b2 = mol2.getBondWithIdx(bond2);
    // RDKit✔️✔️:   auto t1 = b1->getBondType();
    // RDKit✔️✔️:   auto t2 = b2->getBondType();
    // RDKit✔️❌:   if (match.isEqual(t1, t2)) {
    // RDKit✔️❌:     if (p.MatchStereo && !checkBondStereo(p, mol1, bond1, mol2, bond2)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (p.RingMatchesRingOnly) {
    // RDKit✔️❌:       return checkBondRingMatch(p, mol1, bond1, mol2, bond2);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     return true;
    // RDKit✔️❌:   }
    // RDKit✔️✔️:   return false;
    // END RDKIT CPP FUNCTION
    // Local complexity review: order and stereo comparisons are O(1); the
    // ring branch inherits Q111's documented fused-ring scan only when enabled.
    let left_order = mcs_bond(left, left_bond, "left")?.order();
    let right_order = mcs_bond(right, right_bond, "right")?.order();
    if bond_orders_match(left_order, right_order, true) {
        if params.match_stereo && !check_bond_stereo(params, left, left_bond, right, right_bond)? {
            return Ok(false);
        }
        if params.ring_matches_ring_only {
            return check_bond_ring_match(params, left, left_bond, right, right_bond);
        }
        return Ok(true);
    }
    Ok(false)
}

fn mcs_bond_compare_order_exact(
    params: &McsBondCompareParameters,
    left: &SearchTarget<'_>,
    left_bond: usize,
    right: &SearchTarget<'_>,
    right_bond: usize,
) -> Result<bool, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: MCSBondCompareOrderExact
    // RDKit✔️✔️:   static const BondMatchOrderMatrix match(false);  // AROMATIC != SINGLE
    // RDKit✔️✔️:   const auto b1 = mol1.getBondWithIdx(bond1);
    // RDKit✔️✔️:   const auto b2 = mol2.getBondWithIdx(bond2);
    // RDKit✔️✔️:   auto t1 = b1->getBondType();
    // RDKit✔️✔️:   auto t2 = b2->getBondType();
    // RDKit✔️❌:   if (match.isEqual(t1, t2)) {
    // RDKit✔️❌:     if (p.MatchStereo && !checkBondStereo(p, mol1, bond1, mol2, bond2)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (p.RingMatchesRingOnly) {
    // RDKit✔️❌:       return checkBondRingMatch(p, mol1, bond1, mol2, bond2);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     return true;
    // RDKit✔️❌:   }
    // RDKit✔️✔️:   return false;
    // END RDKIT CPP FUNCTION
    // Local complexity review: order and stereo comparisons are O(1); the
    // ring branch inherits Q111's documented fused-ring scan only when enabled.
    let left_order = mcs_bond(left, left_bond, "left")?.order();
    let right_order = mcs_bond(right, right_bond, "right")?.order();
    if bond_orders_match(left_order, right_order, false) {
        if params.match_stereo && !check_bond_stereo(params, left, left_bond, right, right_bond)? {
            return Ok(false);
        }
        if params.ring_matches_ring_only {
            return check_bond_ring_match(params, left, left_bond, right, right_bond);
        }
        return Ok(true);
    }
    Ok(false)
}

fn mcs_final_candidate_accept<E>(
    params: &McsParameters,
    ring_fusion_check: &mut dyn FnMut() -> Result<bool, E>,
    chirality_check: &mut dyn FnMut() -> Result<bool, E>,
    mut user_final_check: Option<&mut dyn FnMut() -> Result<bool, E>>,
) -> Result<bool, E> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: FinalMatchCheckFunction
    // RDKit✔️✔️:   PRECONDITION(p, "p must not be NULL");
    // RDKit✔️✔️:   if ((p->BondCompareParameters.MatchFusedRings ||
    // RDKit✔️✔️:        p->BondCompareParameters.MatchFusedRingsStrict) &&
    // RDKit✔️✔️:       !ringFusionCheck(c1, c2, mol1, query, mol2, target, *p)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (p->AtomCompareParameters.MatchChiralTag &&
    // RDKit✔️✔️:       !FinalChiralityCheckFunction(c1, c2, mol1, query, mol2, target, p)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const auto ip = dynamic_cast<const detail::MCSParametersInternal *>(p);
    // RDKit✔️✔️:   if (ip && ip->UserFinalMatchChecker) {
    // RDKit✔️✔️:     return ip->UserFinalMatchChecker(c1, c2, mol1, query, mol2, target, p);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // END RDKIT CPP FUNCTION
    // Local complexity review: each enabled check is invoked at most once in
    // source order, with constant dispatcher overhead and no allocation.
    if (params.bond_compare_parameters.match_fused_rings
        || params.bond_compare_parameters.match_fused_rings_strict)
        && !ring_fusion_check()?
    {
        return Ok(false);
    }
    if params.atom_compare_parameters.match_chiral_tag && !chirality_check()? {
        return Ok(false);
    }
    if let Some(user_final_check) = user_final_check.as_mut() {
        return user_final_check();
    }
    Ok(true)
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct McsProgressData {
    num_atoms: usize,
    num_bonds: usize,
    seed_processed: u32,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum McsProgressError {
    #[error("MCS progress callback failed: {message}")]
    Callback { message: String },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct McsProgressOutcome<'a> {
    progress: McsProgressData,
    canceled: bool,
    best_atoms: &'a [usize],
    best_bonds: &'a [usize],
}

fn mcs_progress_callback_timeout(params: &McsParameters, start: u64, now: u64) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: MCSProgressCallbackTimeout
    // RDKit✔️✔️:   PRECONDITION(userData, "userData must not be NULL");
    // RDKit✔️✔️:   auto t0 = static_cast<unsigned long long *>(userData);
    // RDKit✔️✔️:   unsigned long long t = nanoClock();
    // RDKit✔️✔️:   return !params.Timeout || (t - *t0 <= params.Timeout * 1000000ULL);
    // END RDKIT CPP FUNCTION
    // Local complexity review: fixed-width arithmetic and one comparison are
    // constant time and allocate no state, matching the source callback.
    params.timeout == 0 || now.wrapping_sub(start) <= u64::from(params.timeout) * 1_000_000
}

fn mcs_progress_after_seed<'a>(
    params: &McsParameters,
    start: u64,
    now: u64,
    best: &'a McsMoleculeFragment,
    progress: &mut McsProgressData,
    callback: Option<
        &mut dyn FnMut(&McsProgressData, &McsParameters) -> Result<bool, McsProgressError>,
    >,
) -> Result<McsProgressOutcome<'a>, McsProgressError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: MaximumCommonSubgraph::growSeeds progress branch
    // RDKit✔️✔️:     if (Parameters.ProgressCallback) {
    // RDKit✔️✔️:       Stat.NumAtoms = getMaxNumberAtoms();
    // RDKit✔️✔️:       Stat.NumBonds = getMaxNumberBonds();
    // RDKit✔️✔️:       if (!Parameters.ProgressCallback(Stat, Parameters,
    // RDKit✔️✔️:                                        Parameters.ProgressCallbackUserData)) {
    // RDKit✔️✔️:         canceled = true;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   if (mcsFound) {  // postponed copy of current set of molecules for
    // RDKit✔️✔️:                    // threshold < 1.
    // RDKit✔️✔️:     McsIdx.QueryMolecule = QueryMolecule;
    // RDKit✔️✔️:     McsIdx.Targets = Targets;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return !canceled;
    // END RDKIT CPP FUNCTION
    // RDKit✔️✔️:   bool Canceled{false};  // interrupted by timeout or user defined progress
    // RDKit✔️✔️:                            // callback. Contains valid current MCS !
    // Local complexity review: progress updates and dispatch are O(1). The
    // outcome borrows the already stored best rows, avoiding a per-step clone.
    progress.num_atoms = best.atoms.len();
    progress.num_bonds = best.bonds.len();
    let keep_running = if let Some(callback) = callback {
        callback(progress, params)?
    } else {
        mcs_progress_callback_timeout(params, start, now)
    };
    Ok(McsProgressOutcome {
        progress: *progress,
        canceled: !keep_running,
        best_atoms: &best.atoms,
        best_bonds: &best.bonds,
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct McsMatchTable {
    rows: usize,
    columns: usize,
    values: Vec<bool>,
}

impl McsMatchTable {
    fn new(rows: usize, columns: usize) -> Self {
        Self {
            rows,
            columns,
            values: vec![false; rows * columns],
        }
    }

    fn set(&mut self, row: usize, column: usize, value: bool) {
        self.values[row * self.columns + column] = value;
    }

    fn get(&self, row: usize, column: usize) -> Option<bool> {
        (row < self.rows && column < self.columns).then(|| self.values[row * self.columns + column])
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct McsMatchTables {
    atoms: McsMatchTable,
    bonds: McsMatchTable,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct McsMoleculeFragment {
    atoms: Vec<usize>,
    bonds: Vec<usize>,
    seed_atom_index_map: BTreeMap<usize, usize>,
}

impl Default for McsMoleculeFragment {
    fn default() -> Self {
        // BEGIN RDKIT CPP TYPE: third_party/rdkit/Code/GraphMol/FMCS/Seed.h :: MolFragment
        // RDKit✔️✔️:   std::vector<const Atom *> Atoms;
        // RDKit✔️✔️:   std::vector<const Bond *> Bonds;
        // RDKit✔️✔️:   std::map<unsigned int, unsigned int> SeedAtomIdxMap;
        // END RDKIT CPP TYPE
        // Local complexity review: source and Rust construct three empty owned
        // containers in constant time without allocating element storage.
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            seed_atom_index_map: BTreeMap::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct McsSeedTopologyBond {
    source_bond: usize,
    begin_seed_atom: usize,
    end_seed_atom: usize,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct McsSeedTopology {
    source_atoms: Vec<usize>,
    bonds: Vec<McsSeedTopologyBond>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct McsNewBond {
    bond_index: usize,
    new_atom_index: usize,
    end_atom_index: Option<usize>,
    new_atom: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct McsTargetMatch {
    empty: bool,
    matched_atom_size: usize,
    matched_bond_size: usize,
    target_atom_indices: Vec<usize>,
    target_bond_indices: Vec<usize>,
    visited_target_bonds: Vec<bool>,
    visited_target_atoms: Vec<bool>,
}

impl Default for McsTargetMatch {
    fn default() -> Self {
        // BEGIN RDKIT CPP TYPE: third_party/rdkit/Code/GraphMol/FMCS/TargetMatch.h :: TargetMatch fields/default constructor
        // RDKit✔️✔️:   bool Empty{true};
        // RDKit✔️✔️:   size_t MatchedAtomSize{0};
        // RDKit✔️✔️:   size_t MatchedBondSize{0};
        // RDKit✔️✔️:   std::vector<unsigned int> TargetAtomIdx;
        // RDKit✔️✔️:   std::vector<unsigned int> TargetBondIdx;
        // RDKit✔️✔️:   boost::dynamic_bitset<> VisitedTargetBonds;
        // RDKit✔️✔️:   boost::dynamic_bitset<> VisitedTargetAtoms;  // for checking rings
        // RDKit✔️✔️:   TargetMatch() {}
        // END RDKIT CPP TYPE
        // Local complexity review: scalar initialization and empty vector
        // construction are constant time without element allocation.
        Self {
            empty: true,
            matched_atom_size: 0,
            matched_bond_size: 0,
            target_atom_indices: Vec::new(),
            target_bond_indices: Vec::new(),
            visited_target_bonds: Vec::new(),
            visited_target_atoms: Vec::new(),
        }
    }
}

impl McsTargetMatch {
    fn clear(&mut self) {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/TargetMatch.h :: TargetMatch::clear
        // RDKit✔️✔️:   void clear() {
        // RDKit✔️✔️:     Empty = true;
        // RDKit✔️✔️:
        // RDKit✔️✔️:     TargetAtomIdx.clear();
        // RDKit✔️✔️:     TargetBondIdx.clear();
        // RDKit✔️✔️:     VisitedTargetBonds.clear();
        // RDKit✔️✔️:     VisitedTargetAtoms.clear();
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION
        // Local complexity review: the four owned buffers are released once,
        // matching the source clear operations without scanning their values.
        self.empty = true;
        self.target_atom_indices.clear();
        self.target_bond_indices.clear();
        self.visited_target_bonds.clear();
        self.visited_target_atoms.clear();
    }

    fn init(
        &mut self,
        seed: &McsSeed,
        matches: &[(usize, usize)],
        query: &SearchTarget<'_>,
        target: &SearchTarget<'_>,
    ) -> Result<(), McsError> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/TargetMatch.cpp :: TargetMatch::init
        // RDKit✔️✔️:   TargetAtomIdx.clear();
        // RDKit✔️✔️:   TargetAtomIdx.resize(query.getNumAtoms(), NotSet);
        // RDKit✔️✔️:   TargetBondIdx.clear();
        // RDKit✔️✔️:   TargetBondIdx.resize(query.getNumBonds(), NotSet);
        // RDKit✔️✔️:   VisitedTargetBonds.resize(target.Molecule->getNumBonds());
        // RDKit✔️✔️:   VisitedTargetAtoms.resize(target.Molecule->getNumAtoms());
        // RDKit✔️✔️:   VisitedTargetBonds.reset();
        // RDKit✔️✔️:   VisitedTargetAtoms.reset();
        // RDKit✔️✔️:
        // RDKit✔️✔️:   MatchedAtomSize = match.size();
        // RDKit✔️✔️:   for (const auto &m : match) {
        // RDKit✔️✔️:     TargetAtomIdx[seed.MoleculeFragment.Atoms.at(m.first)->getIdx()] = m.second;
        // RDKit✔️✔️:     VisitedTargetAtoms.set(m.second);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   MatchedBondSize = 0;
        // RDKit✔️✔️:   for (const auto bond : seed.MoleculeFragment.Bonds) {
        // RDKit✔️✔️:     unsigned int i = bond->getBeginAtomIdx();
        // RDKit✔️✔️:     unsigned int j = bond->getEndAtomIdx();
        // RDKit✔️✔️:     unsigned int ti = TargetAtomIdx.at(i);
        // RDKit✔️✔️:     unsigned int tj = TargetAtomIdx.at(j);
        // RDKit✔️✔️:     const auto tb = target.Molecule->getBondBetweenAtoms(ti, tj);
        // RDKit✔️✔️:     if (tb) {
        // RDKit✔️✔️:       ++MatchedBondSize;
        // RDKit✔️✔️:       TargetBondIdx[bond->getIdx()] = tb->getIdx();  // add
        // RDKit✔️✔️:       VisitedTargetBonds.set(tb->getIdx());
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   Empty = false;
        // END RDKIT CPP FUNCTION
        // Local complexity review: the four arrays allocate once at their
        // source sizes. Atom rows and seed bonds are each visited once, while
        // target bond lookup scans only one endpoint's adjacency, matching the
        // source getBondBetweenAtoms degree-bounded lookup.
        self.target_atom_indices.clear();
        self.target_atom_indices
            .resize(query.num_atoms(), usize::MAX);
        self.target_bond_indices.clear();
        self.target_bond_indices
            .resize(query.num_bonds(), usize::MAX);
        self.visited_target_bonds.clear();
        self.visited_target_bonds.resize(target.num_bonds(), false);
        self.visited_target_atoms.clear();
        self.visited_target_atoms.resize(target.num_atoms(), false);

        self.matched_atom_size = matches.len();
        for &(seed_atom, target_atom) in matches {
            let query_atom = seed.molecule_fragment.atoms.get(seed_atom).copied().ok_or(
                McsError::AtomOutOfRange {
                    side: "seed match",
                    atom: seed_atom,
                },
            )?;
            let mapped =
                self.target_atom_indices
                    .get_mut(query_atom)
                    .ok_or(McsError::AtomOutOfRange {
                        side: "query",
                        atom: query_atom,
                    })?;
            *mapped = target_atom;
            let visited =
                self.visited_target_atoms
                    .get_mut(target_atom)
                    .ok_or(McsError::AtomOutOfRange {
                        side: "target match",
                        atom: target_atom,
                    })?;
            *visited = true;
        }

        self.matched_bond_size = 0;
        for &query_bond in &seed.molecule_fragment.bonds {
            let bond = query
                .bonds()
                .get(query_bond)
                .ok_or(McsError::BondOutOfRange {
                    side: "query",
                    bond: query_bond,
                })?;
            let begin = bond.begin().index();
            let end = bond.end().index();
            let target_begin = self
                .target_atom_indices
                .get(begin)
                .copied()
                .filter(|atom| *atom != usize::MAX)
                .ok_or(McsError::TargetAtomMappingMissing { atom: begin })?;
            let target_end = self
                .target_atom_indices
                .get(end)
                .copied()
                .filter(|atom| *atom != usize::MAX)
                .ok_or(McsError::TargetAtomMappingMissing { atom: end })?;
            if target_begin >= target.num_atoms() {
                return Err(McsError::AtomOutOfRange {
                    side: "target match",
                    atom: target_begin,
                });
            }
            if target_end >= target.num_atoms() {
                return Err(McsError::AtomOutOfRange {
                    side: "target match",
                    atom: target_end,
                });
            }
            let target_bond = target
                .adjacency()
                .neighbors_of(target_begin)
                .iter()
                .find(|neighbor| neighbor.atom_index == target_end)
                .map(|neighbor| neighbor.bond.index());
            if let Some(target_bond) = target_bond {
                self.matched_bond_size += 1;
                self.target_bond_indices[query_bond] = target_bond;
                self.visited_target_bonds[target_bond] = true;
            }
        }
        self.empty = false;
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum McsCandidateMatchError {
    #[error("MCS has {targets} targets but {tables} target match-table pairs")]
    TargetTableCount { targets: usize, tables: usize },
    #[error("MCS threshold target count {threshold} exceeds target count {targets}")]
    ThresholdCountOutOfRange { threshold: usize, targets: usize },
    #[error("MCS {kind} match table has no row {row}, column {column}")]
    MatchTableOutOfRange {
        kind: &'static str,
        row: usize,
        column: usize,
    },
    #[error("MCS source bond {bond} has no ring membership")]
    RingMembershipMissing { bond: usize },
    #[error("MCS ring membership references missing ring {ring}")]
    RingOutOfRange { ring: usize },
    #[error("MCS mapped {side} seed edge has no bond between atoms {begin} and {end}")]
    MappedBondMissing {
        side: &'static str,
        begin: usize,
        end: usize,
    },
    #[error(transparent)]
    StereoOrder(#[from] cosmolkit_core::StereoOrderError),
    #[error("MCS seed has {count} outgoing bonds; the source limit is 64")]
    TooManyNewBonds { count: usize },
    #[error("invalid MCS initial SMARTS: {message}")]
    InitialSeedParse { message: String },
    #[error("MCS initial SMARTS matching failed: {message}")]
    InitialSeedMatch { message: String },
    #[error("MCS initial SMARTS bond {bond} has no mapped query bond")]
    InitialSeedBondMissing { bond: usize },
    #[error("MCS result {kind} value {value} exceeds the modeled integer range")]
    ResultValueOutOfRange { kind: &'static str, value: usize },
    #[error("MCS result query graph is invalid: {message}")]
    ResultQueryGraph { message: String },
    #[error("MCS result SMARTS serialization failed: {message}")]
    ResultSmarts { message: String },
    #[error("MCS seed reconstruction found no target bond between atoms {begin} and {end}")]
    SeedReconstructionBondMissing { begin: usize, end: usize },
    #[error("MCS retained result has no source query and target context")]
    ResultContextMissing,
    #[error(transparent)]
    Progress(#[from] McsProgressError),
    #[error(transparent)]
    State(#[from] McsError),
}

fn mcs_target_bond_between(target: &SearchTarget<'_>, begin: usize, end: usize) -> Option<usize> {
    target
        .adjacency()
        .neighbors_of(begin)
        .iter()
        .find(|neighbor| neighbor.atom_index == end)
        .map(|neighbor| neighbor.bond.index())
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct McsPackedBits {
    words: Vec<u64>,
    len: usize,
}

impl McsPackedBits {
    fn new(len: usize) -> Self {
        Self {
            words: vec![0; len.div_ceil(64)],
            len,
        }
    }

    fn set(&mut self, index: usize, value: bool) {
        debug_assert!(index < self.len);
        let mask = 1_u64 << (index % 64);
        if value {
            self.words[index / 64] |= mask;
        } else {
            self.words[index / 64] &= !mask;
        }
    }

    fn test(&self, index: usize) -> bool {
        debug_assert!(index < self.len);
        self.words[index / 64] & (1_u64 << (index % 64)) != 0
    }

    fn count(&self) -> usize {
        self.words
            .iter()
            .map(|word| word.count_ones() as usize)
            .sum()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct McsRingBondCount {
    all: McsPackedBits,
    nonfused: McsPackedBits,
    fused: McsPackedBits,
    nonfused_count_pass1: usize,
    fused_count_pass1: usize,
}

struct McsRingBondCountVect<'a, 'b> {
    molecule: &'a SearchTarget<'b>,
    ring_info: &'a RingInfo,
    side: &'static str,
    rings: Vec<McsRingBondCount>,
    mcs_bonds: McsPackedBits,
}

impl<'a, 'b> McsRingBondCountVect<'a, 'b> {
    fn new(
        molecule: &'a SearchTarget<'b>,
        side: &'static str,
    ) -> Result<Self, McsCandidateMatchError> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: RingBondCountVect::RingBondCountVect
        // RDKit✔️✔️:   RingBondCountVect(const ROMol &mol)
        // RDKit✔️✔️:       : d_mol(mol), d_ringInfo(mol.getRingInfo()) {
        // RDKit✔️✔️:     d_ringBondCountVect.resize(d_ringInfo->numRings());
        // RDKit✔️✔️:     d_isMCSBond.resize(mol.getNumBonds());
        // RDKit✔️✔️:     for (auto &ringBondCount : d_ringBondCountVect) {
        // RDKit✔️✔️:       ringBondCount.isMCSRingBond.resize(mol.getNumBonds());
        // RDKit✔️✔️:       ringBondCount.isMCSRingBondNonFused.resize(mol.getNumBonds());
        // RDKit✔️✔️:       ringBondCount.isMCSRingBondFused.resize(mol.getNumBonds());
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION
        // Packed u64 masks retain the source bitset's O(bonds/word) storage;
        // construction initializes three per-ring masks and one global mask.
        let ring_info = molecule
            .ring_info()
            .filter(|info| info.is_initialized())
            .ok_or(McsError::MissingRingInfo { side })?;
        let bond_count = molecule.num_bonds();
        let rings = (0..ring_info.num_rings())
            .map(|_| McsRingBondCount {
                all: McsPackedBits::new(bond_count),
                nonfused: McsPackedBits::new(bond_count),
                fused: McsPackedBits::new(bond_count),
                nonfused_count_pass1: 0,
                fused_count_pass1: 0,
            })
            .collect();
        Ok(Self {
            molecule,
            ring_info,
            side,
            rings,
            mcs_bonds: McsPackedBits::new(bond_count),
        })
    }

    fn set_mcs_bond_bits_pass1(
        &mut self,
        begin_seed_atom: usize,
        end_seed_atom: usize,
        mapping: &[usize],
        graph_vertices: &[usize],
    ) -> Result<(), McsCandidateMatchError> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: RingBondCountVect::setMCSBondBitsPass1
        // RDKit✔️✔️:   void setMCSBondBitsPass1(unsigned int beginAtomIdx, unsigned int endAtomIdx,
        // RDKit✔️✔️:                            const std::uint32_t c[], const FMCS::Graph &graph) {
        // RDKit✔️✔️:     const auto bond =
        // RDKit✔️✔️:         d_mol.getBondBetweenAtoms(graph[c[beginAtomIdx]], graph[c[endAtomIdx]]);
        // RDKit✔️✔️:     CHECK_INVARIANT(bond, "");
        // RDKit✔️✔️:     const auto bi = bond->getIdx();
        // RDKit✔️✔️:     d_isMCSBond.set(bi);
        // RDKit✔️✔️:     if (!d_ringInfo->numBondRings(bi)) {
        // RDKit✔️✔️:       return;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (d_ringInfo->numBondRings(bi) == 1) {
        // RDKit✔️✔️:       const auto ringIdx = d_ringInfo->bondMembers(bi).front();
        // RDKit✔️✔️:       d_ringBondCountVect[ringIdx].isMCSRingBond.set(bi);
        // RDKit✔️✔️:       d_ringBondCountVect[ringIdx].isMCSRingBondNonFused.set(bi);
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       for (const auto &ringIdx : d_ringInfo->bondMembers(bi)) {
        // RDKit✔️✔️:         d_ringBondCountVect[ringIdx].isMCSRingBond.set(bi);
        // RDKit✔️✔️:         d_ringBondCountVect[ringIdx].isMCSRingBondFused.set(bi);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION
        // Each mapped edge uses indexed rows and one original adjacency lookup;
        // ring marking visits only this bond's membership list.
        let mapped_atom = |seed_atom| -> Result<usize, McsCandidateMatchError> {
            let vertex = *mapping.get(seed_atom).ok_or(McsError::AtomOutOfRange {
                side: "final mapping",
                atom: seed_atom,
            })?;
            let original = *graph_vertices.get(vertex).ok_or(McsError::AtomOutOfRange {
                side: "final graph vertex",
                atom: vertex,
            })?;
            if original >= self.molecule.num_atoms() {
                return Err(McsError::AtomOutOfRange {
                    side: self.side,
                    atom: original,
                }
                .into());
            }
            Ok(original)
        };
        let begin = mapped_atom(begin_seed_atom)?;
        let end = mapped_atom(end_seed_atom)?;
        let bond = mcs_target_bond_between(self.molecule, begin, end).ok_or(
            McsCandidateMatchError::MappedBondMissing {
                side: self.side,
                begin,
                end,
            },
        )?;
        if bond >= self.molecule.num_bonds() {
            return Err(McsError::BondOutOfRange {
                side: self.side,
                bond,
            }
            .into());
        }
        self.mcs_bonds.set(bond, true);
        let memberships = self.ring_info.bond_members(BondId::new(bond));
        match memberships {
            [] => {}
            [ring] => {
                let row = self
                    .rings
                    .get_mut(*ring)
                    .ok_or(McsCandidateMatchError::RingOutOfRange { ring: *ring })?;
                row.all.set(bond, true);
                row.nonfused.set(bond, true);
            }
            rings => {
                for &ring in rings {
                    let row = self
                        .rings
                        .get_mut(ring)
                        .ok_or(McsCandidateMatchError::RingOutOfRange { ring })?;
                    row.all.set(bond, true);
                    row.fused.set(bond, true);
                }
            }
        }
        Ok(())
    }

    fn set_mcs_bond_bits_pass2(&mut self) -> Result<(), McsCandidateMatchError> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: RingBondCountVect::setMCSBondBitsPass2
        // RDKit✔️✔️:   void setMCSBondBitsPass2() {
        // RDKit✔️✔️:     for (auto &ringBondCount : d_ringBondCountVect) {
        // RDKit✔️✔️:       ringBondCount.nonFusedCountPass1 =
        // RDKit✔️✔️:           ringBondCount.isMCSRingBondNonFused.count();
        // RDKit✔️✔️:       ringBondCount.fusedCountPass1 = ringBondCount.isMCSRingBondFused.count();
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     for (unsigned int bi = 0; bi < d_mol.getNumBonds(); ++bi) {
        // RDKit✔️✔️:       if (!d_isMCSBond.test(bi)) {
        // RDKit✔️✔️:         continue;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       int fusedBondRingIdx = -1;
        // RDKit✔️✔️:       unsigned int fusedBondCount = 0;
        // RDKit✔️✔️:       for (auto &ringBondCount : d_ringBondCountVect) {
        // RDKit✔️✔️:         if (!ringBondCount.isMCSRingBondFused.test(bi)) {
        // RDKit✔️✔️:           continue;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:         auto ringIdx = &ringBondCount - &d_ringBondCountVect.front();
        // RDKit✔️✔️:         if (ringBondCount.nonFusedCountPass1 == 0 &&
        // RDKit✔️✔️:             ringBondCount.fusedCountPass1 <
        // RDKit✔️✔️:                 d_ringInfo->bondRings().at(ringIdx).size()) {
        // RDKit✔️✔️:           ringBondCount.isMCSRingBondFused.set(bi, false);
        // RDKit✔️✔️:         } else {
        // RDKit✔️✔️:           ++fusedBondCount;
        // RDKit✔️✔️:           fusedBondRingIdx = ringIdx;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (fusedBondCount == 1) {
        // RDKit✔️✔️:         d_ringBondCountVect[fusedBondRingIdx].isMCSRingBondNonFused.set(bi);
        // RDKit✔️✔️:         d_ringBondCountVect[fusedBondRingIdx].isMCSRingBondFused.set(bi, false);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION
        // One snapshot scan, then source bond-major/ring-minor traversal over
        // packed masks; no new graph search or mask clone is introduced.
        for row in &mut self.rings {
            row.nonfused_count_pass1 = row.nonfused.count();
            row.fused_count_pass1 = row.fused.count();
        }
        for bond in 0..self.molecule.num_bonds() {
            if !self.mcs_bonds.test(bond) {
                continue;
            }
            let mut fused_ring = None;
            let mut fused_count = 0;
            for (ring, row) in self.rings.iter_mut().enumerate() {
                if !row.fused.test(bond) {
                    continue;
                }
                let ring_size = self
                    .ring_info
                    .bond_rings()
                    .get(ring)
                    .ok_or(McsCandidateMatchError::RingOutOfRange { ring })?
                    .len();
                if row.nonfused_count_pass1 == 0 && row.fused_count_pass1 < ring_size {
                    row.fused.set(bond, false);
                } else {
                    fused_count += 1;
                    fused_ring = Some(ring);
                }
            }
            if fused_count == 1 {
                let ring = fused_ring.expect("one surviving fused row has an index");
                self.rings[ring].nonfused.set(bond, true);
                self.rings[ring].fused.set(bond, false);
            }
        }
        Ok(())
    }

    fn is_ring_fusion_honored(&self) -> Result<bool, McsCandidateMatchError> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: RingBondCountVect::isRingFusionHonored
        // RDKit✔️✔️:   bool isRingFusionHonored() {
        // RDKit✔️✔️:     for (const auto &ringBondCount : d_ringBondCountVect) {
        // RDKit✔️✔️:       unsigned int ringIdx = &ringBondCount - &d_ringBondCountVect.front();
        // RDKit✔️✔️:       const auto &bondRings = d_ringInfo->bondRings().at(ringIdx);
        // RDKit✔️✔️:       const auto numRingBondsInMCS = ringBondCount.isMCSRingBond.count();
        // RDKit✔️✔️:       if (!numRingBondsInMCS || numRingBondsInMCS == bondRings.size()) {
        // RDKit✔️✔️:         continue;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       const auto numNonFusedRingBondsInMCS =
        // RDKit✔️✔️:           ringBondCount.isMCSRingBondNonFused.count();
        // RDKit✔️✔️:       const auto numFusedRingBondsInMCS =
        // RDKit✔️✔️:           ringBondCount.isMCSRingBondFused.count();
        // RDKit✔️✔️:       const auto numMissingFusedBonds = std::count_if(
        // RDKit✔️✔️:           bondRings.begin(), bondRings.end(),
        // RDKit✔️✔️:           [this, &ringBondCount](const auto bi) {
        // RDKit✔️✔️:             return (d_ringInfo->numBondRings(bi) > 1 &&
        // RDKit✔️✔️:                     !ringBondCount.isMCSRingBondNonFused.test(bi) &&
        // RDKit✔️✔️:                     !ringBondCount.isMCSRingBondFused.test(bi));
        // RDKit✔️✔️:           });
        // RDKit✔️✔️:       if (numMissingFusedBonds + numFusedRingBondsInMCS +
        // RDKit✔️✔️:               numNonFusedRingBondsInMCS ==
        // RDKit✔️✔️:           bondRings.size()) {
        // RDKit✔️✔️:         return false;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     return true;
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION
        // Each ring scans only its stored bond row; bit probes stay O(1).
        for (ring, row) in self.rings.iter().enumerate() {
            let ring_bonds = self
                .ring_info
                .bond_rings()
                .get(ring)
                .ok_or(McsCandidateMatchError::RingOutOfRange { ring })?;
            let selected = row.all.count();
            if selected == 0 || selected == ring_bonds.len() {
                continue;
            }
            let mut missing_fused = 0;
            for bond_id in ring_bonds {
                let bond = bond_id.index();
                if bond >= self.molecule.num_bonds() {
                    return Err(McsError::BondOutOfRange {
                        side: self.side,
                        bond,
                    }
                    .into());
                }
                if self.ring_info.num_bond_rings(*bond_id) > 1
                    && !row.nonfused.test(bond)
                    && !row.fused.test(bond)
                {
                    missing_fused += 1;
                }
            }
            if missing_fused + row.fused.count() + row.nonfused.count() == ring_bonds.len() {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

#[allow(clippy::too_many_arguments)]
fn mcs_ring_fusion_check(
    c1: &[usize],
    c2: &[usize],
    query_molecule: &SearchTarget<'_>,
    query_graph: &McsSeedTopology,
    target_molecule: &SearchTarget<'_>,
    target_graph_vertices: &[usize],
    params: &McsParameters,
) -> Result<bool, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: ringFusionCheck
    // RDKit✔️✔️: inline bool ringFusionCheck(const std::uint32_t c1[], const std::uint32_t c2[],
    // RDKit✔️✔️:                             const ROMol &mol1, const FMCS::Graph &query,
    // RDKit✔️✔️:                             const ROMol &mol2, const FMCS::Graph &target,
    // RDKit✔️✔️:                             const MCSParameters &p) {
    // RDKit✔️✔️:   bool res = true;
    // RDKit✔️✔️:   if (boost::num_edges(target) < boost::num_edges(query)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   RingBondCountVect mol1RingBondCountVect(mol1);
    // RDKit✔️✔️:   RingBondCountVect mol2RingBondCountVect(mol2);
    // RDKit✔️✔️:   auto queryEdges = boost::edges(query);
    // RDKit✔️✔️:   std::for_each(queryEdges.first, queryEdges.second,
    // RDKit✔️✔️:                 [&c1, &c2, &query, &target, &mol1RingBondCountVect,
    // RDKit✔️✔️:                  &mol2RingBondCountVect](const auto &edge) {
    // RDKit✔️✔️:                   const auto beginAtomIdx = boost::source(edge, query);
    // RDKit✔️✔️:                   const auto endAtomIdx = boost::target(edge, query);
    // RDKit✔️✔️:                   mol1RingBondCountVect.setMCSBondBitsPass1(
    // RDKit✔️✔️:                       beginAtomIdx, endAtomIdx, c1, query);
    // RDKit✔️✔️:                   mol2RingBondCountVect.setMCSBondBitsPass1(
    // RDKit✔️✔️:                       beginAtomIdx, endAtomIdx, c2, target);
    // RDKit✔️✔️:                 });
    // RDKit✔️✔️:   mol1RingBondCountVect.setMCSBondBitsPass2();
    // RDKit✔️✔️:   mol2RingBondCountVect.setMCSBondBitsPass2();
    // RDKit✔️✔️:   bool mol1Honored = mol1RingBondCountVect.isRingFusionHonored();
    // RDKit✔️✔️:   bool mol2Honored = mol2RingBondCountVect.isRingFusionHonored();
    // RDKit✔️✔️:   if (p.BondCompareParameters.MatchFusedRingsStrict) {
    // RDKit✔️✔️:     res = mol1Honored && mol2Honored;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = mol1Honored || mol2Honored;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    // Both packed counters are built once, and every seed edge is visited once;
    // the later pass2 scans the source bond/ring rows in the source order.
    if target_molecule.num_bonds() < query_graph.bonds.len() {
        return Ok(true);
    }
    let mut query_counts = McsRingBondCountVect::new(query_molecule, "query")?;
    let mut target_counts = McsRingBondCountVect::new(target_molecule, "target")?;
    for edge in &query_graph.bonds {
        query_counts.set_mcs_bond_bits_pass1(
            edge.begin_seed_atom,
            edge.end_seed_atom,
            c1,
            &query_graph.source_atoms,
        )?;
        target_counts.set_mcs_bond_bits_pass1(
            edge.begin_seed_atom,
            edge.end_seed_atom,
            c2,
            target_graph_vertices,
        )?;
    }
    query_counts.set_mcs_bond_bits_pass2()?;
    target_counts.set_mcs_bond_bits_pass2()?;
    let query_honored = query_counts.is_ring_fusion_honored()?;
    let target_honored = target_counts.is_ring_fusion_honored()?;
    Ok(if params.bond_compare_parameters.match_fused_rings_strict {
        query_honored && target_honored
    } else {
        query_honored || target_honored
    })
}

#[allow(clippy::too_many_arguments)]
fn mcs_final_tetrahedral_check(
    c1: &[usize],
    c2: &[usize],
    query_molecule: &SearchTarget<'_>,
    query_graph: &McsSeedTopology,
    target_molecule: &SearchTarget<'_>,
    target_graph_vertices: &[usize],
) -> Result<bool, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: FinalChiralityCheckFunction (active uint32 tetrahedral portion)
    // RDKit❗✔️: bool FinalChiralityCheckFunction(const std::uint32_t c1[],
    // RDKit❗✔️:                                  const std::uint32_t c2[], const ROMol &mol1,
    // RDKit❗✔️:                                  const FMCS::Graph &query, const ROMol &mol2,
    // RDKit❗✔️:                                  const FMCS::Graph &target,
    // RDKit❗✔️:                                  const MCSParameters * /*unused*/) {
    // RDKit❗✔️:   const unsigned int qna = boost::num_vertices(query);  // getNumAtoms()
    // RDKit❗✔️:   // check chiral atoms only:
    // RDKit❗✔️:   for (unsigned int i = 0; i < qna; ++i) {
    // RDKit❗✔️:     const auto a1 = mol1.getAtomWithIdx(query[c1[i]]);
    // RDKit❗✔️:     const auto ac1 = a1->getChiralTag();
    // RDKit❗✔️:
    // RDKit❗✔️:     const auto a2 = mol2.getAtomWithIdx(target[c2[i]]);
    // RDKit❗✔️:     const auto ac2 = a2->getChiralTag();
    // RDKit❗✔️:
    // RDKit❗✔️:     ///*------------------ OLD Code :
    // RDKit❗✔️:     // ???: non chiral query atoms ARE ALLOWED TO MATCH to Chiral target atoms
    // RDKit❗✔️:     // (see test for issue 481)
    // RDKit❗✔️:     if (a1->getDegree() <
    // RDKit❗✔️:             3 ||  // #688: doesn't deal with "explicit" Hs properly
    // RDKit❗✔️:         !(ac1 == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit❗✔️:           ac1 == Atom::CHI_TETRAHEDRAL_CCW)) {
    // RDKit❗✔️:       continue;  // skip non chiral center QUERY atoms
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (!(ac2 == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit❗✔️:           ac2 == Atom::CHI_TETRAHEDRAL_CCW)) {
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     //--------------------
    // RDKit❗✔️:     /* More accurate check:
    // RDKit❗✔️:
    // RDKit❗✔️:             if( !(ac1 == Atom::CHI_TETRAHEDRAL_CW || ac1 ==
    // RDKit❗✔️:        Atom::CHI_TETRAHEDRAL_CCW)
    // RDKit❗✔️:              && !(ac2 == Atom::CHI_TETRAHEDRAL_CW || ac2 ==
    // RDKit❗✔️:        Atom::CHI_TETRAHEDRAL_CCW))
    // RDKit❗✔️:                 continue; // skip check if both atoms are non chiral center
    // RDKit❗✔️:
    // RDKit❗✔️:             if(!(   (ac1 == Atom::CHI_TETRAHEDRAL_CW || ac1 ==
    // RDKit❗✔️:        Atom::CHI_TETRAHEDRAL_CCW)
    // RDKit❗✔️:                  && (ac2 == Atom::CHI_TETRAHEDRAL_CW || ac2 ==
    // RDKit❗✔️:        Atom::CHI_TETRAHEDRAL_CCW)))//ac2 != ac1)
    // RDKit❗✔️:                  return false; // both atoms must be chiral or not without a
    // RDKit❗✔️:        query priority
    // RDKit❗✔️:     */
    // RDKit❗✔️:     const unsigned int a1Degree =
    // RDKit❗✔️:         boost::out_degree(c1[i], query);  // a1.getDegree();
    // RDKit❗✔️:     // number of all connected atoms in a seed
    // RDKit❗✔️:     if (a1Degree > a2->getDegree()) {  // #688 was != . // FIX issue 631
    // RDKit❗✔️:       // printf("atoms Degree (%u, %u) %u [%u], %u\n", query[c1[i]],
    // RDKit❗✔️:       // target[c2[i]], a1Degree, a1.getDegree(), a2.getDegree());
    // RDKit❗✔️:       if (1 == a1Degree && a1->getDegree() == a2->getDegree()) {
    // RDKit❗✔️:         continue;  // continue to grow the seed
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         return false;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     INT_LIST qOrder;
    // RDKit❗✔️:     for (unsigned int j = 0; j < qna && qOrder.size() != a1Degree; ++j) {
    // RDKit❗✔️:       const auto qB = mol1.getBondBetweenAtoms(query[c1[i]], query[c1[j]]);
    // RDKit❗✔️:       if (qB) {
    // RDKit❗✔️:         qOrder.push_back(qB->getIdx());
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // #688
    // RDKit❗✔️:     INT_LIST qmoOrder;
    // RDKit❗✔️:     {
    // RDKit❗✔️:       for (const auto &nbri :
    // RDKit❗✔️:            boost::make_iterator_range(mol1.getAtomBonds(a1))) {
    // RDKit❗✔️:         int dbidx = mol1[nbri]->getIdx();
    // RDKit❗✔️:         if (std::find(qOrder.begin(), qOrder.end(), dbidx) != qOrder.end()) {
    // RDKit❗✔️:           qmoOrder.push_back(dbidx);
    // RDKit❗✔️:         }
    // RDKit❗✔️:         //            else
    // RDKit❗✔️:         //                qmoOrder.push_back(-1);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     int qPermCount =  // was: a1.getPerturbationOrder(qOrder);
    // RDKit❗✔️:         static_cast<int>(countSwapsToInterconvert(qmoOrder, qOrder));
    // RDKit❗✔️:
    // RDKit❗✔️:     INT_LIST mOrder;
    // RDKit❗✔️:     for (unsigned int j = 0; j < qna && mOrder.size() != a2->getDegree(); ++j) {
    // RDKit❗✔️:       const auto mB = mol2.getBondBetweenAtoms(target[c2[i]], target[c2[j]]);
    // RDKit❗✔️:       if (mB) {
    // RDKit❗✔️:         mOrder.push_back(mB->getIdx());
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // #688
    // RDKit❗✔️:     while (mOrder.size() < a2->getDegree()) {
    // RDKit❗✔️:       mOrder.push_back(-1);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     INT_LIST moOrder;
    // RDKit❗✔️:     for (const auto &nbri : boost::make_iterator_range(mol2.getAtomBonds(a2))) {
    // RDKit❗✔️:       int dbidx = mol2[nbri]->getIdx();
    // RDKit❗✔️:       if (std::find(mOrder.begin(), mOrder.end(), dbidx) != mOrder.end()) {
    // RDKit❗✔️:         moOrder.push_back(dbidx);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         moOrder.push_back(-1);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     int mPermCount =  // was: a2.getPerturbationOrder(mOrder);
    // RDKit❗✔️:         static_cast<int>(countSwapsToInterconvert(moOrder, mOrder));
    // RDKit❗✔️:     //----
    // RDKit❗✔️:
    // RDKit❗✔️:     if ((qPermCount % 2 == mPermCount % 2 &&
    // RDKit❗✔️:          a1->getChiralTag() != a2->getChiralTag()) ||
    // RDKit❗✔️:         (qPermCount % 2 != mPermCount % 2 &&
    // RDKit❗✔️:          a1->getChiralTag() == a2->getChiralTag())) {
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION
    // Each mapped neighbor lookup scans source-order adjacency, like the
    // source getBondBetweenAtoms call. The core swap helper retains first-match
    // semantics even when multiple absent ligands use the same sentinel.
    let qna = query_graph.source_atoms.len();
    for i in 0..qna {
        let query_vertex = *c1.get(i).ok_or(McsError::AtomOutOfRange {
            side: "query seed mapping",
            atom: i,
        })?;
        let target_vertex = *c2.get(i).ok_or(McsError::AtomOutOfRange {
            side: "target seed mapping",
            atom: i,
        })?;
        let query_atom_index =
            *query_graph
                .source_atoms
                .get(query_vertex)
                .ok_or(McsError::AtomOutOfRange {
                    side: "query graph vertex",
                    atom: query_vertex,
                })?;
        let target_atom_index =
            *target_graph_vertices
                .get(target_vertex)
                .ok_or(McsError::AtomOutOfRange {
                    side: "target graph vertex",
                    atom: target_vertex,
                })?;
        let query_atom = mcs_atom(query_molecule, query_atom_index, "query molecule")?;
        let target_atom = mcs_atom(target_molecule, target_atom_index, "target molecule")?;
        let query_tag = query_atom.chiral_tag();
        let target_tag = target_atom.chiral_tag();
        let query_full_degree = query_molecule
            .adjacency()
            .neighbors_of(query_atom_index)
            .len();
        if query_full_degree < 3
            || !matches!(
                query_tag,
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            )
        {
            continue;
        }
        if !matches!(
            target_tag,
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            return Ok(false);
        }
        let seed_degree = query_graph
            .bonds
            .iter()
            .filter(|edge| {
                edge.begin_seed_atom == query_vertex || edge.end_seed_atom == query_vertex
            })
            .count();
        let target_full_degree = target_molecule
            .adjacency()
            .neighbors_of(target_atom_index)
            .len();
        if seed_degree > target_full_degree {
            if seed_degree == 1 && query_full_degree == target_full_degree {
                continue;
            }
            return Ok(false);
        }
        let mut query_order = Vec::with_capacity(seed_degree);
        for j in 0..qna {
            if query_order.len() == seed_degree {
                break;
            }
            let neighbor_vertex = *c1.get(j).ok_or(McsError::AtomOutOfRange {
                side: "query seed mapping",
                atom: j,
            })?;
            let neighbor_atom =
                *query_graph
                    .source_atoms
                    .get(neighbor_vertex)
                    .ok_or(McsError::AtomOutOfRange {
                        side: "query graph vertex",
                        atom: neighbor_vertex,
                    })?;
            if let Some(bond) =
                mcs_target_bond_between(query_molecule, query_atom_index, neighbor_atom)
            {
                query_order.push(bond);
            }
        }
        let query_molecule_order: Vec<usize> = query_molecule
            .adjacency()
            .neighbors_of(query_atom_index)
            .iter()
            .map(|neighbor| neighbor.bond.index())
            .filter(|bond| query_order.contains(bond))
            .collect();
        let query_swaps =
            cosmolkit_core::count_swaps_to_interconvert(&query_molecule_order, &query_order)?;

        let mut target_order = Vec::with_capacity(target_full_degree);
        for j in 0..qna {
            if target_order.len() == target_full_degree {
                break;
            }
            let neighbor_vertex = *c2.get(j).ok_or(McsError::AtomOutOfRange {
                side: "target seed mapping",
                atom: j,
            })?;
            let neighbor_atom =
                *target_graph_vertices
                    .get(neighbor_vertex)
                    .ok_or(McsError::AtomOutOfRange {
                        side: "target graph vertex",
                        atom: neighbor_vertex,
                    })?;
            if let Some(bond) =
                mcs_target_bond_between(target_molecule, target_atom_index, neighbor_atom)
            {
                target_order.push(Some(bond));
            }
        }
        target_order.resize(target_full_degree, None);
        let target_molecule_order: Vec<Option<usize>> = target_molecule
            .adjacency()
            .neighbors_of(target_atom_index)
            .iter()
            .map(|neighbor| {
                let bond = neighbor.bond.index();
                target_order.contains(&Some(bond)).then_some(bond)
            })
            .collect();
        let target_swaps =
            cosmolkit_core::count_swaps_to_interconvert(&target_molecule_order, &target_order)?;
        if (query_swaps % 2 == target_swaps % 2 && query_tag != target_tag)
            || (query_swaps % 2 != target_swaps % 2 && query_tag == target_tag)
        {
            return Ok(false);
        }
    }
    Ok(true)
}

#[allow(clippy::too_many_arguments)]
fn mcs_final_chirality_check(
    c1: &[usize],
    c2: &[usize],
    query_molecule: &SearchTarget<'_>,
    query_graph: &McsSeedTopology,
    target_molecule: &SearchTarget<'_>,
    target_graph_vertices: &[usize],
) -> Result<bool, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: FinalChiralityCheckFunction (active uint32 double-bond portion)
    // RDKit❗✔️:   // check double bonds ONLY (why ???)
    // RDKit❗✔️:   std::map<unsigned int, unsigned int> qMap;
    // RDKit❗✔️:   for (unsigned int j = 0; j < qna; ++j) {
    // RDKit❗✔️:     qMap[query[c1[j]]] = j;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (const auto &bondIdx : boost::make_iterator_range(boost::edges(query))) {
    // RDKit❗✔️:     const auto qBnd = mol1.getBondWithIdx(query[bondIdx]);
    // RDKit❗✔️:     if (qBnd->getBondType() != Bond::DOUBLE ||
    // RDKit❗✔️:         qBnd->getStereo() <= Bond::STEREOANY) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // don't think this can actually happen, but check to be sure:
    // RDKit❗✔️:     if (qBnd->getStereoAtoms().size() != 2) {  // MUST check it in the seed, not
    // RDKit❗✔️:                                                // in full query molecule, but
    // RDKit❗✔️:                                                // never happens !!!
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     const auto mBnd =
    // RDKit❗✔️:         mol2.getBondBetweenAtoms(target[c2[qMap[qBnd->getBeginAtomIdx()]]],
    // RDKit❗✔️:                                  target[c2[qMap[qBnd->getEndAtomIdx()]]]);
    // RDKit❗✔️:     CHECK_INVARIANT(mBnd, "Matching bond not found");
    // RDKit❗✔️:     if (mBnd->getBondType() != Bond::DOUBLE ||
    // RDKit❗✔️:         mBnd->getStereo() <= Bond::STEREOANY) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // don't think this can actually happen, but check to be sure:
    // RDKit❗✔️:     if (mBnd->getStereoAtoms().size() != 2) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     unsigned int end1Matches = 0;
    // RDKit❗✔️:     unsigned int end2Matches = 0;
    // RDKit❗✔️:     if (target[c2[qMap[qBnd->getBeginAtomIdx()]]] ==
    // RDKit❗✔️:         rdcast<unsigned int>(mBnd->getBeginAtomIdx())) {
    // RDKit❗✔️:       // query Begin == mol Begin
    // RDKit❗✔️:       if (target[c2[qMap[qBnd->getStereoAtoms()[0]]]] ==
    // RDKit❗✔️:           rdcast<unsigned int>(mBnd->getStereoAtoms()[0])) {
    // RDKit❗✔️:         end1Matches = 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (target[c2[qMap[qBnd->getStereoAtoms()[1]]]] ==
    // RDKit❗✔️:           rdcast<unsigned int>(mBnd->getStereoAtoms()[1])) {
    // RDKit❗✔️:         end2Matches = 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       // query End == mol Begin
    // RDKit❗✔️:       if (target[c2[qMap[qBnd->getStereoAtoms()[0]]]] ==
    // RDKit❗✔️:           rdcast<unsigned int>(mBnd->getStereoAtoms()[1])) {
    // RDKit❗✔️:         end1Matches = 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (target[c2[qMap[qBnd->getStereoAtoms()[1]]]] ==
    // RDKit❗✔️:           rdcast<unsigned int>(mBnd->getStereoAtoms()[0])) {
    // RDKit❗✔️:         end2Matches = 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // std::cerr<<"  bnd: "<<qBnd->getIdx()<<":"<<qBnd->getStereo()<<" -
    // RDKit❗✔️:     // "<<mBnd->getIdx()<<":"<<mBnd->getStereo()<<"  --  "<<end1Matches<<"
    // RDKit❗✔️:     // "<<end2Matches<<std::endl;
    // RDKit❗✔️:     if (mBnd->getStereo() == qBnd->getStereo() &&
    // RDKit❗✔️:         (end1Matches + end2Matches) == 1) {
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (mBnd->getStereo() != qBnd->getStereo() &&
    // RDKit❗✔️:         (end1Matches + end2Matches) != 1) {
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return true;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION
    // The map remains ordered as in std::map; entry insertion reproduces
    // operator[]'s zero default for stereo atoms outside the seed. No graph or
    // molecule state is cloned, and edge traversal stays source ordered.
    if !mcs_final_tetrahedral_check(
        c1,
        c2,
        query_molecule,
        query_graph,
        target_molecule,
        target_graph_vertices,
    )? {
        return Ok(false);
    }
    let mut query_row_by_atom = BTreeMap::new();
    for j in 0..query_graph.source_atoms.len() {
        let vertex = *c1.get(j).ok_or(McsError::AtomOutOfRange {
            side: "query seed mapping",
            atom: j,
        })?;
        let atom = *query_graph
            .source_atoms
            .get(vertex)
            .ok_or(McsError::AtomOutOfRange {
                side: "query graph vertex",
                atom: vertex,
            })?;
        mcs_atom(query_molecule, atom, "query molecule")?;
        query_row_by_atom.insert(atom, j);
    }
    let mut mapped_target_atom =
        |original_query_atom: usize| -> Result<usize, McsCandidateMatchError> {
            let row = *query_row_by_atom.entry(original_query_atom).or_insert(0);
            let vertex = *c2.get(row).ok_or(McsError::AtomOutOfRange {
                side: "target seed mapping",
                atom: row,
            })?;
            let atom = *target_graph_vertices
                .get(vertex)
                .ok_or(McsError::AtomOutOfRange {
                    side: "target graph vertex",
                    atom: vertex,
                })?;
            mcs_atom(target_molecule, atom, "target molecule")?;
            Ok(atom)
        };
    for edge in &query_graph.bonds {
        let query_bond = mcs_bond(query_molecule, edge.source_bond, "query molecule")?;
        if query_bond.order() != BondOrder::Double
            || query_bond.stereo().rdkit_code() <= BondStereo::Any.rdkit_code()
        {
            continue;
        }
        let Some(query_stereo_atoms) = query_bond.stereo_atoms() else {
            continue;
        };
        let begin = mapped_target_atom(query_bond.begin().index())?;
        let end = mapped_target_atom(query_bond.end().index())?;
        let target_bond_index = mcs_target_bond_between(target_molecule, begin, end).ok_or(
            McsCandidateMatchError::MappedBondMissing {
                side: "target",
                begin,
                end,
            },
        )?;
        let target_bond = mcs_bond(target_molecule, target_bond_index, "target molecule")?;
        if target_bond.order() != BondOrder::Double
            || target_bond.stereo().rdkit_code() <= BondStereo::Any.rdkit_code()
        {
            continue;
        }
        let Some(target_stereo_atoms) = target_bond.stereo_atoms() else {
            continue;
        };
        let first = mapped_target_atom(query_stereo_atoms[0].index())?;
        let second = mapped_target_atom(query_stereo_atoms[1].index())?;
        let matched = if begin == target_bond.begin().index() {
            usize::from(first == target_stereo_atoms[0].index())
                + usize::from(second == target_stereo_atoms[1].index())
        } else {
            usize::from(first == target_stereo_atoms[1].index())
                + usize::from(second == target_stereo_atoms[0].index())
        };
        if (query_bond.stereo() == target_bond.stereo() && matched == 1)
            || (query_bond.stereo() != target_bond.stereo() && matched != 1)
        {
            return Ok(false);
        }
    }
    Ok(true)
}

fn mcs_final_mapping_accept(
    topology: &McsSeedTopology,
    query: &SearchTarget<'_>,
    target: &SearchTarget<'_>,
    params: &McsParameters,
    mapping: &[(usize, usize)],
    user_final_check: Option<&mut dyn FnMut() -> Result<bool, McsCandidateMatchError>>,
) -> Result<bool, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/SubstructMatchCustom.cpp :: MolMatchFinalCheckFunctor::operator()
    // RDKit❗❌:   bool operator()(const boost::detail::node_id c1[],
    // RDKit❗❌:                   const boost::detail::node_id c2[]) const {
    // RDKit❗❌:     if (static_cast<unsigned int>(c1[0]) >=
    // RDKit❗❌:         boost::num_vertices(QueryTopology)) {
    // RDKit❗❌:       return false;
    // RDKit❗❌:     }
    // RDKit❗❌:     auto compare = Parameters ? Parameters->FinalMatchChecker : nullptr;
    // RDKit❗❌:     return compare ? compare(c1, c2, d_query, QueryTopology, d_mol,
    // RDKit❗❌:                              TargetTopology, Parameters)
    // RDKit❗❌:                    : true;
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION
    // The mapped pairs are VF2 seed vertex rows and target original atom IDs.
    // The target topology graph has one vertex per original target atom in
    // index order; constructing its identity vector costs O(target atoms).
    let query_vertices = topology.source_atoms.len();
    let c1 = (0..query_vertices).collect::<Vec<_>>();
    let target_vertices = (0..target.num_atoms()).collect::<Vec<_>>();
    let mut c2 = vec![usize::MAX; query_vertices];
    for &(query_vertex, target_atom) in mapping {
        if query_vertex >= query_vertices {
            return Err(McsError::AtomOutOfRange {
                side: "query seed mapping",
                atom: query_vertex,
            }
            .into());
        }
        mcs_atom(target, target_atom, "target molecule")?;
        c2[query_vertex] = target_atom;
    }
    for (query_vertex, target_atom) in c2.iter().copied().enumerate() {
        if target_atom == usize::MAX {
            return Err(McsError::TargetAtomMappingMissing {
                atom: topology.source_atoms[query_vertex],
            }
            .into());
        }
    }
    let mut ring =
        || mcs_ring_fusion_check(&c1, &c2, query, topology, target, &target_vertices, params);
    let mut chirality =
        || mcs_final_chirality_check(&c1, &c2, query, topology, target, &target_vertices);
    mcs_final_candidate_accept(params, &mut ring, &mut chirality, user_final_check)
}

fn mcs_find_full_mapping(
    seed: &McsSeed,
    target: &SearchTarget<'_>,
    tables: &McsMatchTables,
    final_check: &mut dyn FnMut(&[(usize, usize)]) -> Result<bool, McsCandidateMatchError>,
) -> Result<Option<Vec<(usize, usize)>>, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/SubstructMatchCustom.cpp :: SubstructMatchCustomTable
    // RDKit✔️❌:   if (query.m_vertices.size() > target.m_vertices.size()  // query > target
    // RDKit✔️❌:       || query.m_edges.size() > target.m_edges.size()) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   MolMatchFinalCheckFunctor mc(query, target, querySrc, mol, p);
    // RDKit✔️❌:
    // RDKit✔️❌:   AtomTableCompareFunctor ac(query, target, atomMatchTable);
    // RDKit✔️❌:   BondTableCompareFunctor bc(query, target, bondMatchTable);
    // RDKit✔️❌:
    // RDKit✔️❌:   match_V_t dummy_match;
    // RDKit✔️❌:   if (!match) {
    // RDKit✔️❌:     match = &dummy_match;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return boost::vf2(query, target, ac, bc, mc, *match);
    // END RDKIT CPP FUNCTION
    // Local complexity review: this preserves exact subgraph feasibility and
    // source row ordering, but lacks Boost VF2's frontier pruning and can
    // explore more candidates. The performance marker therefore remains ❌.
    let atom_count = seed.topology.source_atoms.len();
    if atom_count > target.num_atoms() || seed.topology.bonds.len() > target.num_bonds() {
        return Ok(None);
    }

    let mut seed_adjacency = vec![Vec::<(usize, usize)>::new(); atom_count];
    for edge in &seed.topology.bonds {
        if edge.begin_seed_atom >= atom_count {
            return Err(McsError::AtomOutOfRange {
                side: "seed topology",
                atom: edge.begin_seed_atom,
            }
            .into());
        }
        if edge.end_seed_atom >= atom_count {
            return Err(McsError::AtomOutOfRange {
                side: "seed topology",
                atom: edge.end_seed_atom,
            }
            .into());
        }
        seed_adjacency[edge.begin_seed_atom].push((edge.end_seed_atom, edge.source_bond));
        seed_adjacency[edge.end_seed_atom].push((edge.begin_seed_atom, edge.source_bond));
    }

    fn visit(
        depth: usize,
        seed: &McsSeed,
        target: &SearchTarget<'_>,
        tables: &McsMatchTables,
        seed_adjacency: &[Vec<(usize, usize)>],
        mapping: &mut [usize],
        used_target_atoms: &mut [bool],
        final_check: &mut dyn FnMut(&[(usize, usize)]) -> Result<bool, McsCandidateMatchError>,
    ) -> Result<bool, McsCandidateMatchError> {
        if depth == mapping.len() {
            let pairs = mapping.iter().copied().enumerate().collect::<Vec<_>>();
            return final_check(&pairs);
        }

        let query_atom = seed.topology.source_atoms[depth];
        for target_atom in 0..target.num_atoms() {
            if used_target_atoms[target_atom] {
                continue;
            }
            let atom_matches = tables.atoms.get(query_atom, target_atom).ok_or(
                McsCandidateMatchError::MatchTableOutOfRange {
                    kind: "atom",
                    row: query_atom,
                    column: target_atom,
                },
            )?;
            if !atom_matches {
                continue;
            }

            let mut feasible = true;
            for &(other_seed_atom, query_bond) in &seed_adjacency[depth] {
                let other_target_atom = mapping[other_seed_atom];
                if other_target_atom == usize::MAX {
                    continue;
                }
                let Some(target_bond) =
                    mcs_target_bond_between(target, target_atom, other_target_atom)
                else {
                    feasible = false;
                    break;
                };
                let bond_matches = tables.bonds.get(query_bond, target_bond).ok_or(
                    McsCandidateMatchError::MatchTableOutOfRange {
                        kind: "bond",
                        row: query_bond,
                        column: target_bond,
                    },
                )?;
                if !bond_matches {
                    feasible = false;
                    break;
                }
            }
            if !feasible {
                continue;
            }

            mapping[depth] = target_atom;
            used_target_atoms[target_atom] = true;
            if visit(
                depth + 1,
                seed,
                target,
                tables,
                seed_adjacency,
                mapping,
                used_target_atoms,
                final_check,
            )? {
                return Ok(true);
            }
            used_target_atoms[target_atom] = false;
            mapping[depth] = usize::MAX;
        }
        Ok(false)
    }

    let mut mapping = vec![usize::MAX; atom_count];
    let mut used_target_atoms = vec![false; target.num_atoms()];
    if visit(
        0,
        seed,
        target,
        tables,
        &seed_adjacency,
        &mut mapping,
        &mut used_target_atoms,
        final_check,
    )? {
        Ok(Some(mapping.into_iter().enumerate().collect()))
    } else {
        Ok(None)
    }
}

fn mcs_match_incremental_fast(
    seed: &McsSeed,
    target_match: &mut McsTargetMatch,
    target: &SearchTarget<'_>,
    tables: &McsMatchTables,
    final_check: &mut dyn FnMut(&[(usize, usize)]) -> Result<bool, McsCandidateMatchError>,
) -> Result<bool, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: MaximumCommonSubgraph::matchIncrementalFast
    // RDKit✔️✔️:   if (match.empty()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool matched = false;
    // RDKit✔️✔️:   for (unsigned int newBondSeedIdx = match.MatchedBondSize;
    // RDKit✔️✔️:        newBondSeedIdx < seed.getNumBonds(); newBondSeedIdx++) {
    // RDKit✔️✔️:     matched = false;
    // RDKit✔️✔️:     bool atomAdded = false;
    // RDKit✔️✔️:     const auto newBond = seed.MoleculeFragment.Bonds.at(newBondSeedIdx);
    // RDKit✔️✔️:     unsigned int newBondQueryIdx = newBond->getIdx();
    // RDKit✔️✔️:     unsigned int newBondSourceAtomSeedIdx;
    // RDKit✔️✔️:     unsigned int newBondOtherAtomSeedIdx;
    // RDKit✔️✔️:     unsigned int i =
    // RDKit✔️✔️:         seed.MoleculeFragment.SeedAtomIdxMap.at(newBond->getBeginAtomIdx());
    // RDKit✔️✔️:     unsigned int j =
    // RDKit✔️✔️:         seed.MoleculeFragment.SeedAtomIdxMap.at(newBond->getEndAtomIdx());
    // RDKit✔️✔️:     if (i >= match.MatchedAtomSize) {
    // RDKit✔️✔️:       newBondSourceAtomSeedIdx = j;
    // RDKit✔️✔️:       newBondOtherAtomSeedIdx = i;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       newBondSourceAtomSeedIdx = i;
    // RDKit✔️✔️:       newBondOtherAtomSeedIdx = j;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     unsigned int newBondOtherAtomQueryIdx =
    // RDKit✔️✔️:         seed.MoleculeFragment.Atoms.at(newBondOtherAtomSeedIdx)->getIdx();
    // RDKit✔️✔️:     unsigned int newBondSourceAtomQueryIdx =
    // RDKit✔️✔️:         seed.MoleculeFragment.Atoms.at(newBondSourceAtomSeedIdx)->getIdx();
    // RDKit✔️✔️:     unsigned int newBondSourceAtomTargetIdx =
    // RDKit✔️✔️:         match.TargetAtomIdx.at(newBondSourceAtomQueryIdx);
    // RDKit✔️✔️:     const Bond *tb = nullptr;
    // RDKit✔️✔️:     unsigned int newBondOtherAtomTargetIdx = NotSet;
    // RDKit✔️✔️:     if (newBondOtherAtomSeedIdx < match.MatchedAtomSize) {
    // RDKit✔️✔️:       newBondOtherAtomTargetIdx =
    // RDKit✔️✔️:           match.TargetAtomIdx.at(newBondOtherAtomQueryIdx);
    // RDKit✔️✔️:       tb = target.Molecule->getBondBetweenAtoms(newBondSourceAtomTargetIdx,
    // RDKit✔️✔️:                                                 newBondOtherAtomTargetIdx);
    // RDKit✔️✔️:       if (tb) {
    // RDKit✔️✔️:         unsigned int tbi = tb->getIdx();
    // RDKit✔️✔️:         unsigned int qbi =
    // RDKit✔️✔️:             seed.MoleculeFragment.Bonds.at(newBondSeedIdx)->getIdx();
    // RDKit✔️✔️:         if (!match.VisitedTargetBonds.test(tbi)) {
    // RDKit✔️✔️:           matched = target.BondMatchTable.at(qbi, tbi);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       for (const auto &nbri :
    // RDKit✔️✔️:            boost::make_iterator_range(target.Molecule->getAtomBonds(atom))) {
    // RDKit✔️✔️:         tb = (*target.Molecule)[nbri];
    // RDKit✔️✔️:         if (match.VisitedTargetBonds.test(tb->getIdx())) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         newBondOtherAtomTargetIdx = tb->getBeginAtomIdx();
    // RDKit✔️✔️:         if (newBondSourceAtomTargetIdx == newBondOtherAtomTargetIdx) {
    // RDKit✔️✔️:           newBondOtherAtomTargetIdx = tb->getEndAtomIdx();
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (match.VisitedTargetAtoms.test(newBondOtherAtomTargetIdx)) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         matched = target.AtomMatchTable.at(newBondOtherAtomQueryIdx,
    // RDKit✔️✔️:                                            newBondOtherAtomTargetIdx) &&
    // RDKit✔️✔️:                   target.BondMatchTable.at(newBondQueryIdx, tb->getIdx());
    // RDKit✔️✔️:         if (matched) {
    // RDKit✔️✔️:           atomAdded = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (matched) {
    // RDKit✔️✔️:       if (atomAdded) {
    // RDKit✔️✔️:         match.MatchedAtomSize++;
    // RDKit✔️✔️:         match.TargetAtomIdx[newBondOtherAtomQueryIdx] =
    // RDKit✔️✔️:             newBondOtherAtomTargetIdx;
    // RDKit✔️✔️:         match.VisitedTargetAtoms.set(newBondOtherAtomTargetIdx);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       match.MatchedBondSize++;
    // RDKit✔️✔️:       match.TargetBondIdx[newBondQueryIdx] = tb->getIdx();
    // RDKit✔️✔️:       match.VisitedTargetBonds.set(tb->getIdx());
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       match.clear();
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (match.MatchedAtomSize != seed.getNumAtoms() ||
    // RDKit✔️✔️:       match.MatchedBondSize != seed.getNumBonds()) {
    // RDKit✔️✔️:     match.clear();
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (matched && Parameters.FinalMatchChecker) {
    // RDKit✔️✔️:     matched = Parameters.FinalMatchChecker(c1.data(), c2.data(), *QueryMolecule,
    // RDKit✔️✔️:                                            seed.Topology, *target.Molecule,
    // RDKit✔️✔️:                                            target.Topology, &Parameters);
    // RDKit✔️✔️:     if (!matched) {
    // RDKit✔️✔️:       match.clear();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return matched;
    // END RDKIT CPP FUNCTION
    // Local complexity review: each new seed bond is visited once; a new atom
    // scans only one target adjacency list, while all cache/table accesses are
    // indexed. Final mapping allocation matches the source c1/c2 vectors.
    if target_match.empty {
        return Ok(false);
    }
    let mut matched = false;
    for new_bond_seed_index in target_match.matched_bond_size..seed.molecule_fragment.bonds.len() {
        matched = false;
        let mut atom_added = false;
        let query_bond = seed.molecule_fragment.bonds[new_bond_seed_index];
        let edge =
            seed.topology
                .bonds
                .get(new_bond_seed_index)
                .ok_or(McsError::BondOutOfRange {
                    side: "seed topology",
                    bond: new_bond_seed_index,
                })?;
        let (source_seed_atom, other_seed_atom) =
            if edge.begin_seed_atom >= target_match.matched_atom_size {
                (edge.end_seed_atom, edge.begin_seed_atom)
            } else {
                (edge.begin_seed_atom, edge.end_seed_atom)
            };
        let other_query_atom = seed
            .molecule_fragment
            .atoms
            .get(other_seed_atom)
            .copied()
            .ok_or(McsError::AtomOutOfRange {
                side: "seed",
                atom: other_seed_atom,
            })?;
        let source_query_atom = seed
            .molecule_fragment
            .atoms
            .get(source_seed_atom)
            .copied()
            .ok_or(McsError::AtomOutOfRange {
                side: "seed",
                atom: source_seed_atom,
            })?;
        let source_target_atom = target_match
            .target_atom_indices
            .get(source_query_atom)
            .copied()
            .filter(|atom| *atom != usize::MAX)
            .ok_or(McsError::TargetAtomMappingMissing {
                atom: source_query_atom,
            })?;
        let mut other_target_atom = usize::MAX;
        let mut target_bond = None;

        if other_seed_atom < target_match.matched_atom_size {
            other_target_atom = target_match
                .target_atom_indices
                .get(other_query_atom)
                .copied()
                .filter(|atom| *atom != usize::MAX)
                .ok_or(McsError::TargetAtomMappingMissing {
                    atom: other_query_atom,
                })?;
            target_bond = mcs_target_bond_between(target, source_target_atom, other_target_atom);
            if let Some(candidate_bond) = target_bond {
                let visited = target_match
                    .visited_target_bonds
                    .get(candidate_bond)
                    .copied()
                    .ok_or(McsError::BondOutOfRange {
                        side: "target match cache",
                        bond: candidate_bond,
                    })?;
                if !visited {
                    matched = tables.bonds.get(query_bond, candidate_bond).ok_or(
                        McsCandidateMatchError::MatchTableOutOfRange {
                            kind: "bond",
                            row: query_bond,
                            column: candidate_bond,
                        },
                    )?;
                }
            }
        } else {
            for neighbor in target.adjacency().neighbors_of(source_target_atom) {
                let candidate_bond = neighbor.bond.index();
                if target_match
                    .visited_target_bonds
                    .get(candidate_bond)
                    .copied()
                    .ok_or(McsError::BondOutOfRange {
                        side: "target match cache",
                        bond: candidate_bond,
                    })?
                {
                    continue;
                }
                let candidate_atom = neighbor.atom_index;
                if target_match
                    .visited_target_atoms
                    .get(candidate_atom)
                    .copied()
                    .ok_or(McsError::AtomOutOfRange {
                        side: "target match cache",
                        atom: candidate_atom,
                    })?
                {
                    continue;
                }
                let atom_matches = tables.atoms.get(other_query_atom, candidate_atom).ok_or(
                    McsCandidateMatchError::MatchTableOutOfRange {
                        kind: "atom",
                        row: other_query_atom,
                        column: candidate_atom,
                    },
                )?;
                let bond_matches = tables.bonds.get(query_bond, candidate_bond).ok_or(
                    McsCandidateMatchError::MatchTableOutOfRange {
                        kind: "bond",
                        row: query_bond,
                        column: candidate_bond,
                    },
                )?;
                matched = atom_matches && bond_matches;
                if matched {
                    atom_added = true;
                    other_target_atom = candidate_atom;
                    target_bond = Some(candidate_bond);
                    break;
                }
            }
        }

        if matched {
            if atom_added {
                target_match.matched_atom_size += 1;
                target_match.target_atom_indices[other_query_atom] = other_target_atom;
                target_match.visited_target_atoms[other_target_atom] = true;
            }
            let target_bond = target_bond.expect("matched incremental bond has target identity");
            target_match.matched_bond_size += 1;
            target_match.target_bond_indices[query_bond] = target_bond;
            target_match.visited_target_bonds[target_bond] = true;
        } else {
            target_match.clear();
            return Ok(false);
        }
    }

    if target_match.matched_atom_size != seed.molecule_fragment.atoms.len()
        || target_match.matched_bond_size != seed.molecule_fragment.bonds.len()
    {
        target_match.clear();
        return Ok(false);
    }
    if matched {
        let mapping = seed
            .topology
            .source_atoms
            .iter()
            .copied()
            .enumerate()
            .map(|(seed_atom, query_atom)| {
                target_match
                    .target_atom_indices
                    .get(query_atom)
                    .copied()
                    .filter(|target_atom| *target_atom != usize::MAX)
                    .map(|target_atom| (seed_atom, target_atom))
                    .ok_or(McsError::TargetAtomMappingMissing { atom: query_atom })
            })
            .collect::<Result<Vec<_>, _>>()?;
        if !final_check(&mapping)? {
            target_match.clear();
            matched = false;
        }
    }
    Ok(matched)
}

fn mcs_match_full_candidate(
    seed: &mut McsSeed,
    query: &SearchTarget<'_>,
    targets: &[SearchTarget<'_>],
    target_tables: &[McsMatchTables],
    threshold_count: usize,
    params: Option<&McsParameters>,
    mut final_check: Option<&mut dyn FnMut(usize, &[(usize, usize)]) -> Result<bool, McsError>>,
) -> Result<bool, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: MaximumCommonSubgraph::match
    // RDKit✔️❌:   unsigned int max_miss = Targets.size() - ThresholdCount;
    // RDKit✔️❌:   unsigned int missing = 0;
    // RDKit✔️❌:   unsigned int passed = 0;
    // RDKit✔️❌:
    // RDKit✔️❌:   for (const auto &tag : Targets) {
    // RDKit✔️❌:     unsigned int itarget = &tag - &Targets.front();
    // RDKit✔️❌:     bool target_matched = false;
    // RDKit✔️❌:     if (!seed.MatchResult.empty() && !seed.MatchResult.at(itarget).empty()) {
    // RDKit✔️❌:       target_matched = matchIncrementalFast(seed, itarget);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (!target_matched) {  // slow full match
    // RDKit✔️❌:       match_V_t match;      // THERE IS NO Bonds match INFO !!!!
    // RDKit✔️❌:       target_matched = SubstructMatchCustomTable(
    // RDKit✔️❌:           tag.Topology, *tag.Molecule, seed.Topology, *QueryMolecule,
    // RDKit✔️❌:           tag.AtomMatchTable, tag.BondMatchTable, &Parameters, &match);
    // RDKit✔️❌:       // save current match info
    // RDKit✔️❌:       if (target_matched) {
    // RDKit✔️❌:         if (seed.MatchResult.empty()) {
    // RDKit✔️❌:           seed.MatchResult.resize(Targets.size());
    // RDKit✔️❌:         }
    // RDKit✔️❌:         seed.MatchResult[itarget].init(seed, match, *QueryMolecule, tag);
    // RDKit✔️❌:       } else if (!seed.MatchResult.empty()) {
    // RDKit✔️❌:         seed.MatchResult[itarget].clear();  //.Empty = true; // == fast clear();
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (target_matched) {
    // RDKit✔️❌:       if (++passed >= ThresholdCount) {  // it's enough
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else {  // mismatched
    // RDKit✔️❌:       if (++missing > max_miss) {
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (missing <= max_miss) {
    // RDKit✔️❌:     return true;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return false;
    // END RDKIT CPP FUNCTION
    // The cache-first branch falls through to the full match only when the
    // source incremental check returns false. Target order and early
    // success/failure thresholds remain source-exact.
    // Its complexity inherits mcs_find_full_mapping's documented VF2 gap.
    if targets.len() != target_tables.len() {
        return Err(McsCandidateMatchError::TargetTableCount {
            targets: targets.len(),
            tables: target_tables.len(),
        });
    }
    if threshold_count > targets.len() {
        return Err(McsCandidateMatchError::ThresholdCountOutOfRange {
            threshold: threshold_count,
            targets: targets.len(),
        });
    }
    let max_miss = targets.len() - threshold_count;
    let mut missing = 0;
    let mut passed = 0;

    for (target_index, (target, tables)) in targets.iter().zip(target_tables).enumerate() {
        let mut accept_mapping = |mapping: &[(usize, usize)]| {
            if let Some(params) = params {
                let has_user_hook = final_check.is_some();
                let mut user_hook = || {
                    final_check.as_deref_mut().expect("present user final hook")(
                        target_index,
                        mapping,
                    )
                    .map_err(McsCandidateMatchError::from)
                };
                let user = if has_user_hook {
                    Some(&mut user_hook as &mut dyn FnMut() -> Result<bool, McsCandidateMatchError>)
                } else {
                    None
                };
                mcs_final_mapping_accept(&seed.topology, query, target, params, mapping, user)
            } else {
                // The source temporarily clears FinalMatchChecker while
                // seeding; this suppresses built-ins and the optional hook.
                Ok(true)
            }
        };
        let mut target_matched = false;
        if !seed.match_result.is_empty() && !seed.match_result[target_index].empty {
            let mut target_match = std::mem::take(&mut seed.match_result[target_index]);
            let incremental = mcs_match_incremental_fast(
                seed,
                &mut target_match,
                target,
                tables,
                &mut accept_mapping,
            );
            seed.match_result[target_index] = target_match;
            target_matched = incremental?;
        }
        let full_mapping = if target_matched {
            None
        } else {
            mcs_find_full_mapping(seed, target, tables, &mut accept_mapping)?
        };
        if target_matched {
            passed += 1;
            if passed >= threshold_count {
                break;
            }
        } else if let Some(mapping) = full_mapping {
            if seed.match_result.is_empty() {
                seed.match_result
                    .resize_with(targets.len(), McsTargetMatch::default);
            }
            let mut target_match = McsTargetMatch::default();
            target_match.init(seed, &mapping, query, target)?;
            seed.match_result[target_index] = target_match;
            passed += 1;
            if passed >= threshold_count {
                break;
            }
        } else {
            if !seed.match_result.is_empty() {
                seed.match_result[target_index].clear();
            }
            missing += 1;
            if missing > max_miss {
                break;
            }
        }
    }
    Ok(missing <= max_miss)
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct McsSeedQueue {
    seeds: Vec<McsSeed>,
}

impl McsSeedQueue {
    fn add(&mut self, seed: &McsSeed) -> usize {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/SeedSet.h :: SeedSet::add
        // RDKit✔️✔️:   Value &add(const Value &seed) {
        // RDKit✔️✔️:     iterator where;
        // RDKit✔️✔️:     for (where = Seeds.begin(); where != Seeds.end();
        // RDKit✔️✔️:          where++) {  // find position in sorted list
        // RDKit✔️✔️:       if (where->getNumBonds() < seed.getNumBonds()) {
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     iterator it = Seeds.insert(where, EmptySeed);
        // RDKit✔️✔️:     Value &val = *it;
        // RDKit✔️✔️:     val.setMoleculeFragment(seed);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     return val;
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Seed.h :: Seed::operator=
        // RDKit✔️✔️:     NewBonds = src.NewBonds;
        // RDKit✔️✔️:     GrowingStage = src.GrowingStage;
        // RDKit✔️✔️:     MoleculeFragment = src.MoleculeFragment;
        // RDKit✔️✔️:     Topology = src.Topology;
        // RDKit✔️✔️:     ExcludedBonds = src.ExcludedBonds;
        // RDKit✔️✔️:     LastAddedAtomsBeginIdx = src.LastAddedAtomsBeginIdx;
        // RDKit✔️✔️:     LastAddedBondsBeginIdx = src.LastAddedBondsBeginIdx;
        // RDKit✔️✔️:     RemainingBonds = src.RemainingBonds;
        // RDKit✔️✔️:     RemainingAtoms = src.RemainingAtoms;
        // RDKit✔️✔️:     StoreAllDegenerateMCS = src.StoreAllDegenerateMCS;
        // RDKit✔️✔️:     MatchResult = src.MatchResult;
        // RDKit✔️✔️:     CopyComplete = true;  // LAST
        // END RDKIT CPP FUNCTION
        // Local complexity review: the insertion search is linear over the
        // source queue and the accepted seed is cloned once, matching the
        // list scan plus source assignment. Equal bond counts remain stable.
        let bond_count = seed.molecule_fragment.bonds.len();
        let position = self
            .seeds
            .iter()
            .position(|queued| queued.molecule_fragment.bonds.len() < bond_count)
            .unwrap_or(self.seeds.len());
        let mut inserted = seed.clone();
        inserted.copy_complete = true;
        self.seeds.insert(position, inserted);
        position
    }
}

fn mcs_check_if_match_and_append(
    seed: &mut McsSeed,
    queue: &mut McsSeedQueue,
    query: &SearchTarget<'_>,
    targets: &[SearchTarget<'_>],
    target_tables: &[McsMatchTables],
    threshold_count: usize,
    params: Option<&McsParameters>,
    final_check: Option<&mut dyn FnMut(usize, &[(usize, usize)]) -> Result<bool, McsError>>,
) -> Result<bool, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: MaximumCommonSubgraph::checkIfMatchAndAppend
    // RDKit✔️❌:   bool found = foundInCache;
    // RDKit✔️❌:
    // RDKit✔️❌:   if (!found) {
    // RDKit✔️❌:     found = match(seed);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   Seed *newSeed = nullptr;
    // RDKit✔️❌:
    // RDKit✔️❌:   {
    // RDKit✔️❌:     if (found) {  // Store new generated seed, if found in cache or in
    // RDKit✔️❌:                       // all(- threshold) targets
    // RDKit✔️❌:       newSeed = &Seeds.add(seed);
    // RDKit✔️❌:       newSeed->CopyComplete = false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (newSeed) {
    // RDKit✔️❌:     *newSeed = seed;  // non-blocking copy for MULTI_THREAD and best CPU
    // RDKit✔️❌:                         // utilization
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   return found;  // new matched seed has been actually added
    // END RDKIT CPP FUNCTION
    // Compile-time duplicate/hash caches are outside the modeled state. The
    // source no-cache path calls the completed matcher once and clones only an
    // accepted seed. Overall performance inherits the full match VF2 gap.
    let found = mcs_match_full_candidate(
        seed,
        query,
        targets,
        target_tables,
        threshold_count,
        params,
        final_check,
    )?;
    if found {
        queue.add(seed);
    }
    Ok(found)
}

fn mcs_can_add_all_non_fused_ring_bonds_connected_to_bond(
    parent: &McsSeed,
    query: &SearchTarget<'_>,
    source_atom: usize,
    source_bond: usize,
    targets: &[SearchTarget<'_>],
    target_tables: &[McsMatchTables],
    threshold_count: usize,
    params: &McsParameters,
    mut final_check: Option<
        &mut dyn FnMut(usize, &[(usize, usize)], &McsParameters) -> Result<bool, McsError>,
    >,
) -> Result<bool, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Seed.cpp :: Seed::canAddAllNonFusedRingBondsConnectedToBond
    // RDKit✔️❌:   const auto &mol = bond.getOwningMol();
    // RDKit✔️❌:   const auto ri = mol.getRingInfo();
    // RDKit✔️❌:   int bondIdx = bond.getIdx();
    // RDKit✔️❌:   const auto &bondRings = ri->bondRings().at(ri->bondMembers(bondIdx).front());
    // RDKit✔️❌:   std::set<unsigned int> nonFusedRingBondIndices;
    // RDKit✔️❌:   boost::dynamic_bitset<> connectedAtomIndices(mol.getNumAtoms());
    // RDKit✔️❌:   Seed seed;
    // RDKit✔️❌:   seed.createFromParent(this);
    // RDKit✔️❌:   for (const auto &bi : bondRings) {
    // RDKit✔️❌:     if (bi != bondIdx && ri->numBondRings(bi) == 1) {
    // RDKit✔️❌:       nonFusedRingBondIndices.insert(bi);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   auto currAtom = &srcAtom;
    // RDKit✔️❌:   auto currBond = &bond;
    // RDKit✔️❌:   auto excludedBonds = seed.ExcludedBonds;
    // RDKit✔️❌:   while (currBond) {
    // RDKit✔️❌:     if (!seed.ExcludedBonds.test(currBond->getIdx())) {
    // RDKit✔️❌:       connectedAtomIndices.set(currBond->getBeginAtomIdx());
    // RDKit✔️❌:       connectedAtomIndices.set(currBond->getEndAtomIdx());
    // RDKit✔️❌:       seed.addNewBondFromAtom(*currAtom, *currBond);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     currBond = nullptr;
    // RDKit✔️❌:     for (const auto &candBondIdx : nonFusedRingBondIndices) {
    // RDKit✔️❌:       const auto candBond = mol.getBondWithIdx(candBondIdx);
    // RDKit✔️❌:       if (connectedAtomIndices.test(candBond->getBeginAtomIdx())) {
    // RDKit✔️❌:         currAtom = candBond->getBeginAtom();
    // RDKit✔️❌:       } else if (connectedAtomIndices.test(candBond->getEndAtomIdx())) {
    // RDKit✔️❌:         currAtom = candBond->getEndAtom();
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         continue;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       nonFusedRingBondIndices.erase(candBondIdx);
    // RDKit✔️❌:       currBond = candBond;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (seed.NewBonds.empty()) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   seed.addNewBondsToSeed(mol, seed);
    // RDKit✔️❌:   seed.MatchResult = MatchResult;
    // RDKit✔️❌:   seed.ExcludedBonds = excludedBonds;
    // RDKit✔️❌:   MCSParameters &p = mcs.parameters();
    // RDKit✔️❌:   bool origMatchFusedRings = p.BondCompareParameters.MatchFusedRings;
    // RDKit✔️❌:   bool origMatchFusedRingsStrict =
    // RDKit✔️❌:       p.BondCompareParameters.MatchFusedRingsStrict;
    // RDKit✔️❌:   p.BondCompareParameters.MatchFusedRings = false;
    // RDKit✔️❌:   p.BondCompareParameters.MatchFusedRingsStrict = false;
    // RDKit✔️❌:   bool res = mcs.match(seed);
    // RDKit✔️❌:   p.BondCompareParameters.MatchFusedRings = origMatchFusedRings;
    // RDKit✔️❌:   p.BondCompareParameters.MatchFusedRingsStrict = origMatchFusedRingsStrict;
    // RDKit✔️❌:   return res;
    // END RDKIT CPP FUNCTION
    // Local complexity review: ring candidates retain std::set ordering and
    // O(log n) removal through BTreeSet. The detached immutable parameter
    // boundary requires cloning the parameter value for the temporary fusion
    // override, and the completed full matcher retains its documented VF2
    // frontier-pruning gap, so the performance marker remains ❌.
    let bond = query
        .bonds()
        .get(source_bond)
        .ok_or(McsError::BondOutOfRange {
            side: "query",
            bond: source_bond,
        })?;
    let ring_info = query
        .ring_info()
        .ok_or(McsError::MissingRingInfo { side: "query" })?;
    let ring = ring_info
        .bond_members(bond.id())
        .first()
        .copied()
        .ok_or(McsCandidateMatchError::RingMembershipMissing { bond: source_bond })?;
    let ring_bonds = ring_info
        .bond_rings()
        .get(ring)
        .ok_or(McsCandidateMatchError::RingOutOfRange { ring })?;
    let mut non_fused_ring_bonds = BTreeSet::new();
    for ring_bond in ring_bonds {
        let ring_bond = ring_bond.index();
        if ring_bond != source_bond && ring_info.num_bond_rings(query.bonds()[ring_bond].id()) == 1
        {
            non_fused_ring_bonds.insert(ring_bond);
        }
    }

    let mut connected_atoms = vec![false; query.num_atoms()];
    let mut candidate = McsSeed::default();
    candidate.create_from_parent(parent);
    let saved_excluded_bonds = candidate.excluded_bonds.clone();
    let mut current = Some((source_atom, source_bond));
    while let Some((current_atom, current_bond)) = current.take() {
        let excluded = candidate.excluded_bonds.get(current_bond).copied().ok_or(
            McsError::SeedExcludedBondOutOfRange {
                bond: current_bond,
                count: candidate.excluded_bonds.len(),
            },
        )?;
        if !excluded {
            let current_bond_ref =
                query
                    .bonds()
                    .get(current_bond)
                    .ok_or(McsError::BondOutOfRange {
                        side: "query",
                        bond: current_bond,
                    })?;
            for endpoint in [
                current_bond_ref.begin().index(),
                current_bond_ref.end().index(),
            ] {
                let connected =
                    connected_atoms
                        .get_mut(endpoint)
                        .ok_or(McsError::AtomOutOfRange {
                            side: "query",
                            atom: endpoint,
                        })?;
                *connected = true;
            }
            candidate.add_new_bond_from_atom(query, current_atom, current_bond)?;
        }

        let mut next = None;
        for &candidate_bond in &non_fused_ring_bonds {
            let candidate_bond_ref =
                query
                    .bonds()
                    .get(candidate_bond)
                    .ok_or(McsError::BondOutOfRange {
                        side: "query",
                        bond: candidate_bond,
                    })?;
            let begin = candidate_bond_ref.begin().index();
            let end = candidate_bond_ref.end().index();
            if connected_atoms.get(begin).copied().unwrap_or(false) {
                next = Some((begin, candidate_bond));
                break;
            }
            if connected_atoms.get(end).copied().unwrap_or(false) {
                next = Some((end, candidate_bond));
                break;
            }
        }
        if let Some((next_atom, next_bond)) = next {
            non_fused_ring_bonds.remove(&next_bond);
            current = Some((next_atom, next_bond));
        }
    }
    if candidate.new_bonds.is_empty() {
        return Ok(false);
    }

    let frontier = McsSeed {
        new_bonds: std::mem::take(&mut candidate.new_bonds),
        remaining_bonds: candidate.remaining_bonds,
        remaining_atoms: candidate.remaining_atoms,
        ..McsSeed::default()
    };
    frontier.add_new_bonds_to_seed(query, &mut candidate)?;
    candidate.new_bonds = frontier.new_bonds;
    candidate.match_result = parent.match_result.clone();
    candidate.excluded_bonds = saved_excluded_bonds;

    let mut relaxed_params = params.clone();
    relaxed_params.bond_compare_parameters.match_fused_rings = false;
    relaxed_params
        .bond_compare_parameters
        .match_fused_rings_strict = false;
    if let Some(check) = final_check.as_deref_mut() {
        let mut relaxed_check =
            |target: usize, mapping: &[(usize, usize)]| check(target, mapping, &relaxed_params);
        mcs_match_full_candidate(
            &mut candidate,
            query,
            targets,
            target_tables,
            threshold_count,
            Some(&relaxed_params),
            Some(&mut relaxed_check),
        )
    } else {
        mcs_match_full_candidate(
            &mut candidate,
            query,
            targets,
            target_tables,
            threshold_count,
            Some(&relaxed_params),
            None,
        )
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct McsInitialSeeds {
    queue: McsSeedQueue,
    query_matched_bonds: usize,
    query_matched_atoms: usize,
    query_single_matched_atom: Option<usize>,
}

fn mcs_make_initial_seeds(
    query: &SearchTarget<'_>,
    targets: &[SearchTarget<'_>],
    target_tables: &[McsMatchTables],
    threshold_count: usize,
    params: &McsParameters,
    run_final_checks: bool,
    mut final_check: Option<&mut dyn FnMut(usize, &[(usize, usize)]) -> Result<bool, McsError>>,
    mut should_accept_single: Option<
        &mut dyn FnMut(usize, usize, usize, &McsParameters) -> Result<bool, McsError>,
    >,
) -> Result<McsInitialSeeds, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: MaximumCommonSubgraph::makeInitialSeeds
    // RDKit✔️❌:   // build a set of initial seeds as "all" single bonds from query
    // RDKit✔️❌:   // molecule
    // RDKit✔️❌:   boost::dynamic_bitset<> excludedBonds(QueryMolecule->getNumBonds());
    // RDKit✔️❌:
    // RDKit✔️❌:   Seeds.clear();
    // RDKit✔️❌:   QueryMoleculeMatchedBonds = 0;
    // RDKit✔️❌:   QueryMoleculeMatchedAtoms = 0;
    // RDKit✔️❌:   QueryMoleculeSingleMatchedAtom = nullptr;
    // RDKit✔️❌:   if (!Parameters.InitialSeed.empty()) {  // make user defined seed
    // RDKit✔️❌:     std::unique_ptr<const ROMol> initialSeedMolecule(
    // RDKit✔️❌:         static_cast<const ROMol *>(SmartsToMol(Parameters.InitialSeed)));
    // RDKit✔️❌:     // make a set of of seed as indices and pointers to current query
    // RDKit✔️❌:     // molecule items based on matching results
    // RDKit✔️❌:     std::vector<MatchVectType> matching_substructs;
    // RDKit✔️❌:     SubstructMatch(*QueryMolecule, *initialSeedMolecule, matching_substructs);
    // RDKit✔️❌:     // loop throw all fragments of Query matched to initial seed
    // RDKit✔️❌:     for (const auto &ms : matching_substructs) {
    // RDKit✔️❌:       Seed seed;
    // RDKit✔️❌:       seed.setStoreAllDegenerateMCS(Parameters.StoreAll);
    // RDKit✔️❌:       seed.ExcludedBonds = excludedBonds;
    // RDKit✔️❌:       seed.MatchResult.resize(Targets.size());
    // RDKit✔️❌:       // add all matched atoms of the matched query fragment
    // RDKit✔️❌:       std::map<unsigned int, unsigned int> initialSeedToQueryAtom;
    // RDKit✔️❌:       for (const auto &msb : ms) {
    // RDKit✔️❌:         unsigned int qai = msb.second;
    // RDKit✔️❌:         unsigned int sai = msb.first;
    // RDKit✔️❌:         seed.addAtom(QueryMolecule->getAtomWithIdx(qai));
    // RDKit✔️❌:         initialSeedToQueryAtom[sai] = qai;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       // add all bonds (existed in initial seed !!!) between all matched
    // RDKit✔️❌:       // atoms in query
    // RDKit✔️❌:       for (const auto &msb : ms) {
    // RDKit✔️❌:         const auto atom = initialSeedMolecule->getAtomWithIdx(msb.first);
    // RDKit✔️❌:         for (const auto &nbri : boost::make_iterator_range(
    // RDKit✔️❌:                  initialSeedMolecule->getAtomBonds(atom))) {
    // RDKit✔️❌:           const auto initialBond = (*initialSeedMolecule)[nbri];
    // RDKit✔️❌:           unsigned int qai1 =
    // RDKit✔️❌:               initialSeedToQueryAtom.at(initialBond->getBeginAtomIdx());
    // RDKit✔️❌:           unsigned int qai2 =
    // RDKit✔️❌:               initialSeedToQueryAtom.at(initialBond->getEndAtomIdx());
    // RDKit✔️❌:
    // RDKit✔️❌:           const auto b = QueryMolecule->getBondBetweenAtoms(qai1, qai2);
    // RDKit✔️❌:           CHECK_INVARIANT(b, "bond must not be NULL");
    // RDKit✔️❌:           if (!seed.ExcludedBonds.test(b->getIdx())) {
    // RDKit✔️❌:             seed.addBond(b);
    // RDKit✔️❌:             seed.ExcludedBonds.set(b->getIdx());
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       seed.computeRemainingSize(*QueryMolecule);
    // RDKit✔️❌:
    // RDKit✔️❌:       if (checkIfMatchAndAppend(seed)) {
    // RDKit✔️❌:         QueryMoleculeMatchedBonds = seed.getNumBonds();
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (Seeds.empty()) {  // create a set of seeds from each query bond
    // RDKit✔️❌:     std::vector<WeightedBond> wbVec;
    // RDKit✔️❌:     wbVec.reserve(QueryMolecule->getNumBonds());
    // RDKit✔️❌:     for (const auto &bond : QueryMolecule->bonds()) {
    // RDKit✔️❌:       wbVec.emplace_back(bond);
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     for (const auto &wb : wbVec) {
    // RDKit✔️❌:       Seed seed;
    // RDKit✔️❌:       seed.setStoreAllDegenerateMCS(Parameters.StoreAll);
    // RDKit✔️❌:       seed.MatchResult.resize(Targets.size());
    // RDKit✔️❌:       seed.addAtom(wb.BondPtr->getBeginAtom());
    // RDKit✔️❌:       seed.addAtom(wb.BondPtr->getEndAtom());
    // RDKit✔️❌:       seed.ExcludedBonds = excludedBonds;  // all bonds from first to current
    // RDKit✔️❌:       seed.addBond(wb.BondPtr);
    // RDKit✔️❌:       excludedBonds.set(wb.BondPtr->getIdx());
    // RDKit✔️❌:
    // RDKit✔️❌:       seed.computeRemainingSize(*QueryMolecule);
    // RDKit✔️❌:
    // RDKit✔️❌:       if (checkIfMatchAndAppend(seed)) {
    // RDKit✔️❌:         ++QueryMoleculeMatchedBonds;
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         // disable (mark as already processed) mismatched bond in all
    // RDKit✔️❌:         // seeds
    // RDKit✔️❌:         for (auto &Seed : Seeds) {
    // RDKit✔️❌:           Seed.ExcludedBonds.set(wb.BondPtr->getIdx());
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   auto nq = QueryMolecule->getNumAtoms();
    // RDKit✔️❌:   MatchVectType singleAtomPairMatch(1);
    // RDKit✔️❌:   MatchVectType emptyBondPairMatch;
    // RDKit✔️❌:   for (size_t i = 0; i < nq; i++) {  // all query's atoms
    // RDKit✔️❌:     const auto queryMolAtom = QueryMolecule->getAtomWithIdx(i);
    // RDKit✔️❌:     bool isQueryMolAtomInRing = queryIsAtomInRing(queryMolAtom);
    // RDKit✔️❌:     unsigned int matched = 0;
    // RDKit✔️❌:     const Atom *candQueryMoleculeSingleMatchedAtom = nullptr;
    // RDKit✔️❌:     for (const auto &tag : Targets) {
    // RDKit✔️❌:       auto nt = tag.Molecule->getNumAtoms();
    // RDKit✔️❌:       for (size_t aj = 0; aj < nt; aj++) {
    // RDKit✔️❌:         if (tag.AtomMatchTable.at(i, aj)) {
    // RDKit✔️❌:           const auto targetMolAtom = tag.Molecule->getAtomWithIdx(aj);
    // RDKit✔️❌:           bool isTargetMolAtomInRing = queryIsAtomInRing(targetMolAtom);
    // RDKit✔️❌:           ++matched;
    // RDKit✔️❌:           if (!(Parameters.BondCompareParameters.CompleteRingsOnly &&
    // RDKit✔️❌:                 (isQueryMolAtomInRing || isTargetMolAtomInRing))) {
    // RDKit✔️❌:             bool shouldAccept = !Parameters.ShouldAcceptMCS;
    // RDKit✔️❌:             if (!shouldAccept) {
    // RDKit✔️❌:               shouldAccept = Parameters.ShouldAcceptMCS(
    // RDKit✔️❌:                   *QueryMolecule, *tag.Molecule, singleAtomPairMatch,
    // RDKit✔️❌:                   emptyBondPairMatch, &Parameters);
    // RDKit✔️❌:             }
    // RDKit✔️❌:             if (shouldAccept) {
    // RDKit✔️❌:               candQueryMoleculeSingleMatchedAtom = queryMolAtom;
    // RDKit✔️❌:             }
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (matched && matched >= ThresholdCount) {
    // RDKit✔️❌:       ++QueryMoleculeMatchedAtoms;
    // RDKit✔️❌:       if (candQueryMoleculeSingleMatchedAtom) {
    // RDKit✔️❌:         if (!QueryMoleculeSingleMatchedAtom) {
    // RDKit✔️❌:           QueryMoleculeSingleMatchedAtom = candQueryMoleculeSingleMatchedAtom;
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           QueryMoleculeSingleMatchedAtom =
    // RDKit✔️❌:               (std::max)(candQueryMoleculeSingleMatchedAtom,
    // RDKit✔️❌:                          QueryMoleculeSingleMatchedAtom,
    // RDKit✔️❌:                          [](const Atom *a, const Atom *b) {
    // RDKit✔️❌:                            if (a->getDegree() != b->getDegree()) {
    // RDKit✔️❌:                              return (a->getDegree() < b->getDegree());
    // RDKit✔️❌:                            } else if (a->getFormalCharge() !=
    // RDKit✔️❌:                                       b->getFormalCharge()) {
    // RDKit✔️❌:                              return (a->getFormalCharge() <
    // RDKit✔️❌:                                      b->getFormalCharge());
    // RDKit✔️❌:                            } else if (a->getAtomicNum() != b->getAtomicNum()) {
    // RDKit✔️❌:                              return (a->getAtomicNum() < b->getAtomicNum());
    // RDKit✔️❌:                            }
    // RDKit✔️❌:                            return (a->getIdx() < b->getIdx());
    // RDKit✔️❌:                          });
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION
    // Local complexity review: parsing and substructure matching reuse the
    // canonical query parser/VF2 path, whose result allocation gap is already
    // documented. Seed construction, source bond traversal, exclusion updates
    // and single-atom scans retain the source asymptotic costs and ordering.
    let mut result = McsInitialSeeds::default();
    let mut excluded_bonds = vec![false; query.num_bonds()];

    if !params.initial_seed.is_empty() {
        let initial =
            crate::parse_smarts(&params.initial_seed, &crate::SmartsParseParams::default())
                .map_err(|error| McsCandidateMatchError::InitialSeedParse {
                    message: error.to_string(),
                })?;
        let matches = crate::try_get_substruct_matches_with_params(
            query,
            &initial,
            &crate::SubstructMatchParams::default(),
        )
        .map_err(|error| McsCandidateMatchError::InitialSeedMatch {
            message: error.to_string(),
        })?;
        for matched in matches {
            let mut seed = McsSeed {
                store_all_degenerate_mcs: params.store_all,
                excluded_bonds: excluded_bonds.clone(),
                match_result: vec![McsTargetMatch::default(); targets.len()],
                ..McsSeed::default()
            };
            for &query_atom in &matched.atom_mapping {
                if query_atom >= query.num_atoms() {
                    return Err(McsError::AtomOutOfRange {
                        side: "initial SMARTS match",
                        atom: query_atom,
                    }
                    .into());
                }
                seed.add_atom(query_atom);
            }
            for initial_atom in 0..initial.num_atoms() {
                for &(_, initial_bond) in &initial.adjacency()[initial_atom] {
                    let bond = &initial.bonds()[initial_bond];
                    let query_begin = matched.atom_mapping[bond.begin().index()];
                    let query_end = matched.atom_mapping[bond.end().index()];
                    let query_bond = mcs_target_bond_between(query, query_begin, query_end).ok_or(
                        McsCandidateMatchError::InitialSeedBondMissing { bond: initial_bond },
                    )?;
                    if !seed.excluded_bonds[query_bond] {
                        seed.add_bond(query, query_bond)?;
                    }
                }
            }
            seed.compute_remaining_size(query)?;
            let accepted = if let Some(check) = final_check.as_deref_mut() {
                mcs_check_if_match_and_append(
                    &mut seed,
                    &mut result.queue,
                    query,
                    targets,
                    target_tables,
                    threshold_count,
                    run_final_checks.then_some(params),
                    Some(check),
                )?
            } else {
                mcs_check_if_match_and_append(
                    &mut seed,
                    &mut result.queue,
                    query,
                    targets,
                    target_tables,
                    threshold_count,
                    run_final_checks.then_some(params),
                    None,
                )?
            };
            if accepted {
                result.query_matched_bonds = seed.molecule_fragment.bonds.len();
            }
        }
    }

    if result.queue.seeds.is_empty() {
        for source_bond in 0..query.num_bonds() {
            let bond = &query.bonds()[source_bond];
            let mut seed = McsSeed {
                store_all_degenerate_mcs: params.store_all,
                match_result: vec![McsTargetMatch::default(); targets.len()],
                ..McsSeed::default()
            };
            seed.add_atom(bond.begin().index());
            seed.add_atom(bond.end().index());
            seed.excluded_bonds = excluded_bonds.clone();
            seed.add_bond(query, source_bond)?;
            excluded_bonds[source_bond] = true;
            seed.compute_remaining_size(query)?;
            let accepted = if let Some(check) = final_check.as_deref_mut() {
                mcs_check_if_match_and_append(
                    &mut seed,
                    &mut result.queue,
                    query,
                    targets,
                    target_tables,
                    threshold_count,
                    run_final_checks.then_some(params),
                    Some(check),
                )?
            } else {
                mcs_check_if_match_and_append(
                    &mut seed,
                    &mut result.queue,
                    query,
                    targets,
                    target_tables,
                    threshold_count,
                    run_final_checks.then_some(params),
                    None,
                )?
            };
            if accepted {
                result.query_matched_bonds += 1;
            } else {
                for queued in &mut result.queue.seeds {
                    queued.excluded_bonds[source_bond] = true;
                }
            }
        }
    }

    for query_atom in 0..query.num_atoms() {
        let query_in_ring = query
            .ring_info()
            .is_some_and(|rings| rings.num_atom_rings(query.atoms()[query_atom].id()) > 0);
        let mut matched_targets = 0;
        let mut candidate = None;
        for (target_index, (target, tables)) in targets.iter().zip(target_tables).enumerate() {
            for target_atom in 0..target.num_atoms() {
                if tables.atoms.get(query_atom, target_atom).ok_or(
                    McsCandidateMatchError::MatchTableOutOfRange {
                        kind: "atom",
                        row: query_atom,
                        column: target_atom,
                    },
                )? {
                    let target_in_ring = target.ring_info().is_some_and(|rings| {
                        rings.num_atom_rings(target.atoms()[target_atom].id()) > 0
                    });
                    matched_targets += 1;
                    if !(params.bond_compare_parameters.complete_rings_only
                        && (query_in_ring || target_in_ring))
                    {
                        let accepted = if let Some(check) = should_accept_single.as_deref_mut() {
                            check(query_atom, target_index, target_atom, params)?
                        } else {
                            true
                        };
                        if accepted {
                            candidate = Some(query_atom);
                        }
                    }
                    break;
                }
            }
        }
        if matched_targets > 0 && matched_targets >= threshold_count {
            result.query_matched_atoms += 1;
            if let Some(candidate) = candidate {
                let better = |left: usize, right: usize| {
                    let left_atom = &query.atoms()[left];
                    let right_atom = &query.atoms()[right];
                    (
                        query.adjacency().neighbors_of(left).len(),
                        left_atom.formal_charge(),
                        left_atom.atomic_number(),
                        left,
                    ) > (
                        query.adjacency().neighbors_of(right).len(),
                        right_atom.formal_charge(),
                        right_atom.atomic_number(),
                        right,
                    )
                };
                if result
                    .query_single_matched_atom
                    .is_none_or(|current| better(candidate, current))
                {
                    result.query_single_matched_atom = Some(candidate);
                }
            }
        }
    }
    Ok(result)
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct McsDuplicateSeedKey {
    atom_indices: Vec<usize>,
    bond_indices: Vec<usize>,
}

impl McsDuplicateSeedKey {
    fn add_atom(&mut self, atom: usize) {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/DuplicatedSeedCache.h :: DuplicatedSeedCache::TKey::addAtom
        // RDKit✔️✔️:       auto it = std::lower_bound(AtomIdx.begin(), AtomIdx.end(), i);
        // RDKit✔️✔️:       AtomIdx.insert(it, i);
        // END RDKIT CPP FUNCTION
        // Local complexity review: binary search plus vector insertion has the
        // same O(log n) comparison and O(n) movement cost as the source.
        let position = self.atom_indices.partition_point(|current| *current < atom);
        self.atom_indices.insert(position, atom);
    }

    fn add_bond(&mut self, bond: usize) {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/DuplicatedSeedCache.h :: DuplicatedSeedCache::TKey::addBond
        // RDKit✔️✔️:       auto it = std::lower_bound(BondIdx.begin(), BondIdx.end(), i);
        // RDKit✔️✔️:       BondIdx.insert(it, i);
        // END RDKIT CPP FUNCTION
        // Local complexity review: binary search plus vector insertion has the
        // same O(log n) comparison and O(n) movement cost as the source.
        let position = self.bond_indices.partition_point(|current| *current < bond);
        self.bond_indices.insert(position, bond);
    }

    fn source_equal(&self, right: &Self) -> bool {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/DuplicatedSeedCache.h :: DuplicatedSeedCache::TKey::operator==
        // RDKit✔️✔️:       return AtomIdx.size() == right.AtomIdx.size() &&
        // RDKit✔️✔️:              BondIdx.size() == right.BondIdx.size() &&
        // RDKit✔️✔️:              0 == std::memcmp(&AtomIdx[0], &right.AtomIdx[0],
        // RDKit✔️✔️:                               AtomIdx.size() * sizeof(unsigned int)) &&
        // RDKit✔️✔️:              0 == std::memcmp(&BondIdx[0], &right.BondIdx[0],
        // RDKit✔️✔️:                               BondIdx.size() * sizeof(unsigned int));
        // END RDKIT CPP FUNCTION
        // Local complexity review: slice equality performs the same length
        // checks and linear contiguous-value comparison without allocation.
        self.atom_indices == right.atom_indices && self.bond_indices == right.bond_indices
    }
}

type McsBitSet = u64;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct McsComposition2N {
    bits: McsBitSet,
    inverse_bits: McsBitSet,
    max_value: McsBitSet,
    value_mask: McsBitSet,
}

impl McsComposition2N {
    fn new(max_value: McsBitSet, value_mask: McsBitSet) -> Self {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Composition2N.h :: Composition2N::Composition2N
        // RDKit✔️✔️: typedef unsigned long long BitSet;
        // RDKit✔️✔️: class Composition2N {  // generator of 2^N-1 possible bit combinations
        // RDKit✔️✔️:   BitSet Bits, InverseBits;
        // RDKit✔️✔️:   BitSet MaxValue, ValueMask;  // need for inverse bitset must be 2^N-1
        // RDKit✔️✔️:  public:
        // RDKit✔️✔️:   Composition2N(BitSet maxValue, BitSet valueMask)
        // RDKit✔️✔️:       : Bits(0), InverseBits(0), MaxValue(maxValue), ValueMask(valueMask) {}
        // END RDKIT CPP FUNCTION
        // Local complexity review: both representations contain four 64-bit
        // scalars and construction performs the same four constant-time stores.
        Self {
            bits: 0,
            inverse_bits: 0,
            max_value,
            value_mask,
        }
    }

    fn compute_2n(power: u32) -> McsBitSet {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Composition2N.h :: Composition2N::compute2N
        // RDKit✔️✔️:   static void compute2N(unsigned int power, BitSet &value) {
        // RDKit✔️✔️:     value = 1uLL << power;
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION
        // Local complexity review: the Rust shift is the same constant-time
        // 64-bit operation for the source-supported power range below 64.
        1_u64 << power
    }

    const fn bit_set(&self) -> McsBitSet {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Composition2N.h :: Composition2N::getBitSet
        // RDKit✔️✔️:   BitSet getBitSet() const {
        // RDKit✔️✔️:     return InverseBits;  // inverse to generate biggest seed first and then
        // RDKit✔️✔️:                          // decrease number of external bonds
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION
        self.inverse_bits
    }

    fn generate_next(&mut self) -> bool {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Composition2N.h :: Composition2N::generateNext
        // RDKit✔️✔️:   bool generateNext() {
        // RDKit✔️✔️:     if ((++Bits) <= MaxValue) {
        // RDKit✔️✔️:       InverseBits = (~Bits + 1) & ValueMask;
        // RDKit✔️✔️:       return true;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION
        // Local complexity review: both implementations use one 64-bit
        // increment, comparison, two's-complement negation and mask in O(1).
        self.bits = self.bits.wrapping_add(1);
        if self.bits <= self.max_value {
            self.inverse_bits = (!self.bits).wrapping_add(1) & self.value_mask;
            true
        } else {
            false
        }
    }

    fn is_2_power(&self) -> bool {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Composition2N.h :: Composition2N::is2Power
        // RDKit✔️✔️:   bool is2Power() const {  // one bit is set only
        // RDKit✔️✔️:     BitSet bits = getBitSet();
        // RDKit✔️✔️:     unsigned int n = 0;
        // RDKit✔️✔️:     while (0 == (bits & 1uLL) &&
        // RDKit✔️✔️:            ++n < sizeof(bits) * 8) {  // find lowest bitwise 1
        // RDKit✔️✔️:       bits >>= 1u;                    // shift all zero lower bits
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (0 != (bits & 1uLL)) {
        // RDKit✔️✔️:       bits >>= 1u;  // shift first set bit too
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     return 0 == bits;  // remained bits except lowest 1
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION
        // Local complexity review: the loop performs at most 64 single-bit
        // shifts exactly as the source and allocates no temporary container.
        let mut bits = self.bit_set();
        let mut n = 0_u32;
        while bits & 1 == 0 && {
            n += 1;
            n < McsBitSet::BITS
        } {
            bits >>= 1;
        }
        if bits & 1 != 0 {
            bits >>= 1;
        }
        bits == 0
    }

    fn is_set(&self, bit: u32) -> bool {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Composition2N.h :: Composition2N::isSet
        // RDKit✔️✔️:   bool isSet(unsigned int bit) const {
        // RDKit✔️✔️:     return 0 != (getBitSet() & (1uLL << bit));
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION
        // Local complexity review: both implementations perform one 64-bit
        // shift, mask and comparison in constant time.
        self.bit_set() & (1_u64 << bit) != 0
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct McsSeed {
    new_bonds: Vec<McsNewBond>,
    store_all_degenerate_mcs: bool,
    copy_complete: bool,
    growing_stage: u32,
    molecule_fragment: McsMoleculeFragment,
    topology: McsSeedTopology,
    excluded_bonds: Vec<bool>,
    last_added_atoms_begin_index: usize,
    last_added_bonds_begin_index: usize,
    remaining_bonds: usize,
    remaining_atoms: usize,
    duplicate_key: McsDuplicateSeedKey,
    match_result: Vec<McsTargetMatch>,
}

impl Default for McsSeed {
    fn default() -> Self {
        // BEGIN RDKIT CPP TYPE: third_party/rdkit/Code/GraphMol/FMCS/Seed.h :: Seed fields/default constructor
        // RDKit✔️✔️:   mutable std::vector<NewBond> NewBonds;
        // RDKit✔️✔️:   bool StoreAllDegenerateMCS = false;
        // RDKit✔️✔️:   bool CopyComplete{false};
        // RDKit✔️✔️:   mutable unsigned int GrowingStage{0};
        // RDKit✔️✔️:   MolFragment MoleculeFragment;
        // RDKit✔️✔️:   Graph Topology;
        // RDKit✔️✔️:   boost::dynamic_bitset<> ExcludedBonds;
        // RDKit✔️✔️:   unsigned int LastAddedAtomsBeginIdx{0};
        // RDKit✔️✔️:   unsigned int LastAddedBondsBeginIdx{0};
        // RDKit✔️✔️:   unsigned int RemainingBonds{0};
        // RDKit✔️✔️:   unsigned int RemainingAtoms{0};
        // RDKit✔️✔️:   DuplicatedSeedCache::TKey DupCacheKey;
        // RDKit✔️✔️:   std::vector<TargetMatch> MatchResult;
        // RDKit✔️✔️:   Seed()
        // RDKit✔️✔️:   {}
        // END RDKIT CPP TYPE
        // Local complexity review: all scalar defaults and empty owned
        // containers match source initialization in O(1), without allocation.
        Self {
            new_bonds: Vec::new(),
            store_all_degenerate_mcs: false,
            copy_complete: false,
            growing_stage: 0,
            molecule_fragment: McsMoleculeFragment::default(),
            topology: McsSeedTopology::default(),
            excluded_bonds: Vec::new(),
            last_added_atoms_begin_index: 0,
            last_added_bonds_begin_index: 0,
            remaining_bonds: 0,
            remaining_atoms: 0,
            duplicate_key: McsDuplicateSeedKey::default(),
            match_result: Vec::new(),
        }
    }
}

impl McsSeed {
    fn create_from_parent(&mut self, parent: &Self) {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Seed.h :: Seed::createFromParent
        // RDKit✔️✔️:     MoleculeFragment = parent->MoleculeFragment;
        // RDKit✔️✔️:     Topology = parent->Topology;
        // RDKit✔️✔️:     ExcludedBonds = parent->ExcludedBonds;
        // RDKit✔️✔️:     RemainingBonds = parent->RemainingBonds;
        // RDKit✔️✔️:     RemainingAtoms = parent->RemainingAtoms;
        // RDKit✔️✔️:     StoreAllDegenerateMCS = parent->StoreAllDegenerateMCS;
        // RDKit✔️✔️:     DupCacheKey = parent->DupCacheKey;
        // RDKit✔️✔️:     LastAddedAtomsBeginIdx = getNumAtoms();  // previous size
        // RDKit✔️✔️:     LastAddedBondsBeginIdx = getNumBonds();  // previous size
        // RDKit✔️✔️:     GrowingStage = 0;
        // END RDKIT CPP FUNCTION
        // Local complexity review: each source-owned vector/map/bitset copy is
        // cloned once, with the same linear element cost and no extra scans.
        self.molecule_fragment = parent.molecule_fragment.clone();
        self.topology = parent.topology.clone();
        self.excluded_bonds = parent.excluded_bonds.clone();
        self.remaining_bonds = parent.remaining_bonds;
        self.remaining_atoms = parent.remaining_atoms;
        self.store_all_degenerate_mcs = parent.store_all_degenerate_mcs;
        self.duplicate_key = parent.duplicate_key.clone();
        self.last_added_atoms_begin_index = self.molecule_fragment.atoms.len();
        self.last_added_bonds_begin_index = self.molecule_fragment.bonds.len();
        self.growing_stage = 0;
    }

    fn add_atom(&mut self, source_atom: usize) -> usize {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Seed.cpp :: Seed::addAtom
        // RDKit✔️✔️:   unsigned int i = MoleculeFragment.Atoms.size();
        // RDKit✔️✔️:   unsigned int aqi = atom->getIdx();
        // RDKit✔️✔️:   MoleculeFragment.Atoms.push_back(atom);
        // RDKit✔️✔️:   MoleculeFragment.SeedAtomIdxMap[aqi] = i;
        // RDKit✔️✔️:   Topology.addAtom(aqi);
        // RDKit✔️✔️:   return i;
        // END RDKIT CPP FUNCTION
        // Local complexity review: vector appends remain amortized O(1), and
        // BTreeMap insertion matches the source std::map O(log n) update.
        let seed_atom = self.molecule_fragment.atoms.len();
        self.molecule_fragment.atoms.push(source_atom);
        self.molecule_fragment
            .seed_atom_index_map
            .insert(source_atom, seed_atom);
        self.topology.source_atoms.push(source_atom);
        // RDKit✔️✔️:   DupCacheKey.addAtom(aqi);
        self.duplicate_key.add_atom(source_atom);
        seed_atom
    }

    fn add_bond(
        &mut self,
        query: &SearchTarget<'_>,
        source_bond: usize,
    ) -> Result<usize, McsError> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Seed.cpp :: Seed::addBond
        // RDKit✔️✔️:   unsigned int bi = bond->getIdx();
        // RDKit✔️✔️:   CHECK_INVARIANT(!ExcludedBonds.test(bi), "");
        // RDKit✔️✔️:   ExcludedBonds.set(bi);
        // RDKit✔️✔️:   MoleculeFragment.Bonds.push_back(bond);
        // RDKit✔️✔️:   // remap idx to seed's indices:
        // RDKit✔️✔️:   unsigned int i = MoleculeFragment.SeedAtomIdxMap.at(bond->getBeginAtomIdx());
        // RDKit✔️✔️:   unsigned int j = MoleculeFragment.SeedAtomIdxMap.at(bond->getEndAtomIdx());
        // RDKit✔️✔️:   Topology.addBond(bi, i, j);
        // RDKit✔️✔️:   return getNumBonds();
        // END RDKIT CPP FUNCTION
        // Local complexity review: source bond and exclusion rows use direct
        // indexed access, fragment/topology appends remain amortized O(1), and
        // BTreeMap endpoint lookups match the source std::map O(log n) cost.
        let bond = query
            .bonds()
            .get(source_bond)
            .ok_or(McsError::BondOutOfRange {
                side: "query",
                bond: source_bond,
            })?;
        let excluded_bond_count = self.excluded_bonds.len();
        let excluded = self.excluded_bonds.get_mut(source_bond).ok_or(
            McsError::SeedExcludedBondOutOfRange {
                bond: source_bond,
                count: excluded_bond_count,
            },
        )?;
        if *excluded {
            return Err(McsError::SeedBondAlreadyExcluded { bond: source_bond });
        }
        *excluded = true;
        self.molecule_fragment.bonds.push(source_bond);

        let begin_source_atom = bond.begin().index();
        let begin_seed_atom = self
            .molecule_fragment
            .seed_atom_index_map
            .get(&begin_source_atom)
            .copied()
            .ok_or(McsError::SeedAtomMappingMissing {
                atom: begin_source_atom,
            })?;
        let end_source_atom = bond.end().index();
        let end_seed_atom = self
            .molecule_fragment
            .seed_atom_index_map
            .get(&end_source_atom)
            .copied()
            .ok_or(McsError::SeedAtomMappingMissing {
                atom: end_source_atom,
            })?;
        self.topology.bonds.push(McsSeedTopologyBond {
            source_bond,
            begin_seed_atom,
            end_seed_atom,
        });
        // RDKit✔️✔️:   DupCacheKey.addBond(bi);
        self.duplicate_key.add_bond(source_bond);
        Ok(self.molecule_fragment.bonds.len())
    }

    fn add_new_bond_from_atom(
        &mut self,
        query: &SearchTarget<'_>,
        source_atom: usize,
        source_bond: usize,
    ) -> Result<(), McsError> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Seed.cpp :: Seed::addNewBondFromAtom
        // RDKit✔️✔️:   const auto end_atom = bond.getOtherAtom(&srcAtom);
        // RDKit✔️✔️:   unsigned int end_atom_idx = NotSet;
        // RDKit✔️✔️:   for (unsigned int i = 0; i < getNumAtoms(); ++i) {
        // RDKit✔️✔️:     // already exists in this seed
        // RDKit✔️✔️:     if (end_atom == MoleculeFragment.Atoms.at(i)) {
        // RDKit✔️✔️:       end_atom_idx = i;
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   NewBonds.emplace_back(bond.getIdx(), end_atom->getIdx(), end_atom_idx,
        // RDKit✔️✔️:                         NotSet == end_atom_idx ? end_atom : nullptr);
        // END RDKIT CPP FUNCTION
        // Local complexity review: the endpoint validation is constant time,
        // the first matching seed atom is found by the same single linear scan,
        // and the frontier append remains amortized O(1).
        if source_atom >= query.num_atoms() {
            return Err(McsError::AtomOutOfRange {
                side: "query",
                atom: source_atom,
            });
        }
        let bond = query
            .bonds()
            .get(source_bond)
            .ok_or(McsError::BondOutOfRange {
                side: "query",
                bond: source_bond,
            })?;
        let begin = bond.begin().index();
        let end = bond.end().index();
        let new_atom_index = if source_atom == begin {
            end
        } else if source_atom == end {
            begin
        } else {
            return Err(McsError::SeedAtomNotBondEndpoint {
                atom: source_atom,
                bond: source_bond,
            });
        };
        let end_atom_index = self
            .molecule_fragment
            .atoms
            .iter()
            .position(|atom| *atom == new_atom_index);
        self.new_bonds.push(McsNewBond {
            bond_index: source_bond,
            new_atom_index,
            end_atom_index,
            new_atom: end_atom_index.is_none().then_some(new_atom_index),
        });
        Ok(())
    }

    fn add_new_bonds_to_seed(
        &self,
        query: &SearchTarget<'_>,
        seed: &mut Self,
    ) -> Result<Vec<bool>, McsError> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Seed.cpp :: Seed::addNewBondsToSeed
        // RDKit✔️✔️:   boost::dynamic_bitset<> newAtomsSet(qmol.getNumAtoms());
        // RDKit✔️✔️:   for (const auto &newBond : NewBonds) {
        // RDKit✔️✔️:     unsigned int aIdx = newBond.EndAtomIdx;
        // RDKit✔️✔️:     if (NotSet == aIdx) {  // new atom
        // RDKit✔️✔️:       // check if new bonds simultaneously close a ring
        // RDKit✔️✔️:       if (!newAtomsSet.test(newBond.NewAtomIdx)) {
        // RDKit✔️✔️:         const auto end_atom = newBond.NewAtom;
        // RDKit✔️✔️:         aIdx = seed.addAtom(end_atom);
        // RDKit✔️✔️:         newAtomsSet.set(newBond.NewAtomIdx);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     const auto src_bond = qmol.getBondWithIdx(newBond.BondIdx);
        // RDKit✔️✔️:     seed.addBond(src_bond);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   seed.RemainingBonds = RemainingBonds - NewBonds.size();  // Added ALL !!!
        // RDKit✔️✔️:   seed.RemainingAtoms =
        // RDKit✔️✔️:       RemainingAtoms - newAtomsSet.count();  // new atoms added to seed
        // RDKit✔️✔️:   return newAtomsSet;
        // END RDKIT CPP FUNCTION
        // Local complexity review: one query-sized bit vector is allocated;
        // each frontier row is visited once and delegates to the source-shaped
        // amortized vector/map operations without an additional global scan.
        let mut new_atoms_set = vec![false; query.num_atoms()];
        for new_bond in &self.new_bonds {
            if new_bond.end_atom_index.is_none() {
                let was_added = new_atoms_set.get_mut(new_bond.new_atom_index).ok_or(
                    McsError::AtomOutOfRange {
                        side: "frontier",
                        atom: new_bond.new_atom_index,
                    },
                )?;
                if !*was_added {
                    let end_atom = new_bond.new_atom.ok_or(McsError::SeedNewAtomMissing {
                        atom: new_bond.new_atom_index,
                    })?;
                    seed.add_atom(end_atom);
                    *was_added = true;
                }
            }
            seed.add_bond(query, new_bond.bond_index)?;
        }
        seed.remaining_bonds = self
            .remaining_bonds
            .checked_sub(self.new_bonds.len())
            .ok_or(McsError::SeedRemainingBondUnderflow {
                remaining: self.remaining_bonds,
                frontier: self.new_bonds.len(),
            })?;
        let added_atoms = new_atoms_set.iter().filter(|added| **added).count();
        seed.remaining_atoms = self.remaining_atoms.checked_sub(added_atoms).ok_or(
            McsError::SeedRemainingAtomUnderflow {
                remaining: self.remaining_atoms,
                added: added_atoms,
            },
        )?;
        Ok(new_atoms_set)
    }

    fn can_grow_bigger_than(&self, max_bonds: usize, max_atoms: usize) -> bool {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Seed.h :: Seed::canGrowBiggerThan
        // RDKit✔️✔️:   bool canGrowBiggerThan(unsigned int maxBonds, unsigned int maxAtoms) const {
        // RDKit✔️✔️:     return RemainingBonds + getNumBonds() > maxBonds ||
        // RDKit✔️✔️:            (RemainingBonds + getNumBonds() == maxBonds &&
        // RDKit✔️✔️:             (RemainingAtoms + getNumAtoms() > maxAtoms ||
        // RDKit✔️✔️:              (StoreAllDegenerateMCS &&
        // RDKit✔️✔️:               RemainingAtoms + getNumAtoms() == maxAtoms)));
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION
        // Local complexity review: both implementations perform a constant
        // number of scalar additions and comparisons without allocation.
        let possible_bonds = self
            .remaining_bonds
            .wrapping_add(self.molecule_fragment.bonds.len());
        let possible_atoms = self
            .remaining_atoms
            .wrapping_add(self.molecule_fragment.atoms.len());
        possible_bonds > max_bonds
            || (possible_bonds == max_bonds
                && (possible_atoms > max_atoms
                    || (self.store_all_degenerate_mcs && possible_atoms == max_atoms)))
    }

    #[allow(clippy::too_many_arguments)]
    fn grow_extensions(
        &mut self,
        queue: &mut McsSeedQueue,
        query: &SearchTarget<'_>,
        targets: &[SearchTarget<'_>],
        target_tables: &[McsMatchTables],
        threshold_count: usize,
        max_bonds: usize,
        max_atoms: usize,
        params: &McsParameters,
        mut final_check: Option<&mut dyn FnMut(usize, &[(usize, usize)]) -> Result<bool, McsError>>,
    ) -> Result<(), McsCandidateMatchError> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Seed.cpp :: Seed::grow extension generation
        // RDKit✔️❌:   if (0 == GrowingStage) {
        // RDKit✔️❌:     // 1. Check and add the biggest child seed with all outgoing bonds added:
        // RDKit✔️❌:     // Add all bonds at first (build the biggest child seed). All new atoms are
        // RDKit✔️❌:     // already in the seed
        // RDKit✔️❌:     Seed seed;
        // RDKit✔️❌:     seed.createFromParent(this);
        // RDKit✔️❌:     newAtomsSet = addNewBondsToSeed(qmol, seed);
        // RDKit✔️❌:     if (!seed.canGrowBiggerThan(mcs.getMaxNumberBonds(),
        // RDKit✔️❌:                                 mcs.getMaxNumberAtoms())) {
        // RDKit✔️❌:       GrowingStage = NotSet;
        // RDKit✔️❌:       return;
        // RDKit✔️❌:     }
        // RDKit✔️❌:     seed.MatchResult = MatchResult;
        // RDKit✔️❌:     bool allMatched = mcs.checkIfMatchAndAppend(seed);
        // RDKit✔️❌:     GrowingStage = 1;
        // RDKit✔️❌:     if (allMatched && NewBonds.size() > 1) {
        // RDKit✔️❌:       return;  // grow deep first. postpone next growing steps
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // RDKit✔️❌:   if (1 == NewBonds.size()) {
        // RDKit✔️❌:     GrowingStage = NotSet;
        // RDKit✔️❌:     return;  // everything has been done
        // RDKit✔️❌:   }
        // END RDKIT CPP FUNCTION
        // Matching delegates to the completed exact-behavior candidate helper,
        // whose documented VF2 frontier-pruning gap keeps this performance
        // marker at ❌. Seed creation and source branch order otherwise retain
        // the same allocations and traversal shape.
        if self.growing_stage == 0 {
            let mut seed = McsSeed::default();
            seed.create_from_parent(self);
            self.add_new_bonds_to_seed(query, &mut seed)?;
            if !seed.can_grow_bigger_than(max_bonds, max_atoms) {
                self.growing_stage = u32::MAX;
                return Ok(());
            }
            seed.match_result.clone_from(&self.match_result);
            let all_matched = if let Some(check) = final_check.as_deref_mut() {
                mcs_check_if_match_and_append(
                    &mut seed,
                    queue,
                    query,
                    targets,
                    target_tables,
                    threshold_count,
                    Some(params),
                    Some(check),
                )?
            } else {
                mcs_check_if_match_and_append(
                    &mut seed,
                    queue,
                    query,
                    targets,
                    target_tables,
                    threshold_count,
                    Some(params),
                    None,
                )?
            };
            self.growing_stage = 1;
            if all_matched && self.new_bonds.len() > 1 {
                return Ok(());
            }
        }
        if self.new_bonds.len() == 1 {
            self.growing_stage = u32::MAX;
            return Ok(());
        }

        // RDKit✔️❌:   unsigned int numErasedNewBonds = 0;
        // RDKit✔️❌:   for (auto &newBond : NewBonds) {
        // RDKit✔️❌:     Seed seed;
        // RDKit✔️❌:     seed.createFromParent(this);
        // RDKit✔️❌:     unsigned int aIdx = newBond.EndAtomIdx;
        // RDKit✔️❌:     if (NotSet == aIdx) {  // new atom
        // RDKit✔️❌:       const auto end_atom = newBond.NewAtom;
        // RDKit✔️❌:       aIdx = seed.addAtom(end_atom);
        // RDKit✔️❌:     }
        // RDKit✔️❌:     const auto src_bond = qmol.getBondWithIdx(newBond.BondIdx);
        // RDKit✔️❌:     seed.addBond(src_bond);
        // RDKit✔️❌:     seed.computeRemainingSize(qmol);
        // RDKit✔️❌:     if (seed.canGrowBiggerThan(mcs.getMaxNumberBonds(),
        // RDKit✔️❌:                                mcs.getMaxNumberAtoms())) {
        // RDKit✔️❌:       if (!MatchResult.empty()) {
        // RDKit✔️❌:         seed.MatchResult = MatchResult;
        // RDKit✔️❌:       }
        // RDKit✔️❌:       if (!mcs.checkIfMatchAndAppend(seed)) {
        // RDKit✔️❌:         newBond.BondIdx = NotSet;
        // RDKit✔️❌:         ++numErasedNewBonds;
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        let mut erased = vec![false; self.new_bonds.len()];
        let mut num_erased_new_bonds = 0_usize;
        for (new_bond_index, new_bond) in self.new_bonds.iter().cloned().enumerate() {
            let mut seed = McsSeed::default();
            seed.create_from_parent(self);
            if new_bond.end_atom_index.is_none() {
                seed.add_atom(new_bond.new_atom.ok_or(McsError::SeedNewAtomMissing {
                    atom: new_bond.new_atom_index,
                })?);
            }
            seed.add_bond(query, new_bond.bond_index)?;
            seed.compute_remaining_size(query)?;
            if seed.can_grow_bigger_than(max_bonds, max_atoms) {
                if !self.match_result.is_empty() {
                    seed.match_result.clone_from(&self.match_result);
                }
                let matched = if let Some(check) = final_check.as_deref_mut() {
                    mcs_check_if_match_and_append(
                        &mut seed,
                        queue,
                        query,
                        targets,
                        target_tables,
                        threshold_count,
                        Some(params),
                        Some(check),
                    )?
                } else {
                    mcs_check_if_match_and_append(
                        &mut seed,
                        queue,
                        query,
                        targets,
                        target_tables,
                        threshold_count,
                        Some(params),
                        None,
                    )?
                };
                if !matched {
                    erased[new_bond_index] = true;
                    num_erased_new_bonds += 1;
                }
            }
        }

        // RDKit✔️✔️:   if (numErasedNewBonds > 0) {
        // RDKit✔️✔️:     std::vector<NewBond> dirtyNewBonds;
        // RDKit✔️✔️:     dirtyNewBonds.reserve(NewBonds.size());
        // RDKit✔️✔️:     dirtyNewBonds.swap(NewBonds);
        // RDKit✔️✔️:     for (const auto &dirtyNewBond : dirtyNewBonds) {
        // RDKit✔️✔️:       if (NotSet != dirtyNewBond.BondIdx) {
        // RDKit✔️✔️:         NewBonds.push_back(dirtyNewBond);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        if num_erased_new_bonds > 0 {
            self.new_bonds = self
                .new_bonds
                .drain(..)
                .zip(erased)
                .filter_map(|(bond, erased)| (!erased).then_some(bond))
                .collect();
        }

        // RDKit✔️❌:   if (NewBonds.size() > 1) {
        // RDKit✔️❌:     if (sizeof(unsigned long long) * 8 < NewBonds.size()) {
        // RDKit✔️❌:       throw std::runtime_error(
        // RDKit✔️❌:           "Max number of new external bonds of a seed >64");
        // RDKit✔️❌:     }
        // RDKit✔️❌:     BitSet maxCompositionValue;
        // RDKit✔️❌:     Composition2N::compute2N(NewBonds.size(), maxCompositionValue);
        // RDKit✔️❌:     --maxCompositionValue;  // 2^N-1
        // RDKit✔️❌:     Composition2N composition(maxCompositionValue, maxCompositionValue);
        // RDKit✔️❌:     while (composition.generateNext()) {
        // RDKit✔️❌:       if (composition.is2Power()) {
        // RDKit✔️❌:         continue;
        // RDKit✔️❌:       }
        // RDKit✔️❌:       if (0 == numErasedNewBonds &&
        // RDKit✔️❌:           composition.getBitSet() == maxCompositionValue) {
        // RDKit✔️❌:         continue;
        // RDKit✔️❌:       }
        // RDKit✔️❌:       Seed seed;
        // RDKit✔️❌:       seed.createFromParent(this);
        // RDKit✔️❌:       newAtomsSet.reset();
        // RDKit✔️❌:       for (const auto &newBond : NewBonds) {
        // RDKit✔️❌:         const auto i = &newBond - &NewBonds.front();
        // RDKit✔️❌:         if (composition.isSet(i)) {
        // RDKit✔️❌:           unsigned int aIdx = newBond.EndAtomIdx;
        // RDKit✔️❌:           if (NotSet == aIdx) {  // new atom
        // RDKit✔️❌:             if (!newAtomsSet.test(newBond.NewAtomIdx)) {
        // RDKit✔️❌:               const auto end_atom = newBond.NewAtom;
        // RDKit✔️❌:               aIdx = seed.addAtom(end_atom);
        // RDKit✔️❌:               newAtomsSet.set(newBond.NewAtomIdx);
        // RDKit✔️❌:             }
        // RDKit✔️❌:           }
        // RDKit✔️❌:           const auto src_bond = qmol.getBondWithIdx(newBond.BondIdx);
        // RDKit✔️❌:           seed.addBond(src_bond);
        // RDKit✔️❌:         }
        // RDKit✔️❌:       }
        // RDKit✔️❌:       seed.computeRemainingSize(qmol);
        // RDKit✔️❌:       if (seed.canGrowBiggerThan(mcs.getMaxNumberBonds(),
        // RDKit✔️❌:                                   mcs.getMaxNumberAtoms())) {
        // RDKit✔️❌:         seed.MatchResult = MatchResult;
        // RDKit✔️❌:         bool found = mcs.checkIfMatchAndAppend(seed);
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        if self.new_bonds.len() > 1 {
            if self.new_bonds.len() > McsBitSet::BITS as usize {
                return Err(McsCandidateMatchError::TooManyNewBonds {
                    count: self.new_bonds.len(),
                });
            }
            let max_composition_value = if self.new_bonds.len() == McsBitSet::BITS as usize {
                McsBitSet::MAX
            } else {
                McsComposition2N::compute_2n(self.new_bonds.len() as u32) - 1
            };
            let mut composition =
                McsComposition2N::new(max_composition_value, max_composition_value);
            let mut new_atoms_set = vec![false; query.num_atoms()];
            while composition.generate_next() {
                if composition.is_2_power() {
                    continue;
                }
                if num_erased_new_bonds == 0 && composition.bit_set() == max_composition_value {
                    continue;
                }
                let mut seed = McsSeed::default();
                seed.create_from_parent(self);
                new_atoms_set.fill(false);
                for (new_bond_index, new_bond) in self.new_bonds.iter().enumerate() {
                    if !composition.is_set(new_bond_index as u32) {
                        continue;
                    }
                    if new_bond.end_atom_index.is_none() {
                        let added = new_atoms_set.get_mut(new_bond.new_atom_index).ok_or(
                            McsError::AtomOutOfRange {
                                side: "frontier",
                                atom: new_bond.new_atom_index,
                            },
                        )?;
                        if !*added {
                            seed.add_atom(new_bond.new_atom.ok_or(
                                McsError::SeedNewAtomMissing {
                                    atom: new_bond.new_atom_index,
                                },
                            )?);
                            *added = true;
                        }
                    }
                    seed.add_bond(query, new_bond.bond_index)?;
                }
                seed.compute_remaining_size(query)?;
                if seed.can_grow_bigger_than(max_bonds, max_atoms) {
                    seed.match_result.clone_from(&self.match_result);
                    if let Some(check) = final_check.as_deref_mut() {
                        mcs_check_if_match_and_append(
                            &mut seed,
                            queue,
                            query,
                            targets,
                            target_tables,
                            threshold_count,
                            Some(params),
                            Some(check),
                        )?;
                    } else {
                        mcs_check_if_match_and_append(
                            &mut seed,
                            queue,
                            query,
                            targets,
                            target_tables,
                            threshold_count,
                            Some(params),
                            None,
                        )?;
                    }
                }
            }
        }
        // RDKit✔️✔️:   GrowingStage = NotSet;  // finished
        self.growing_stage = u32::MAX;
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn grow(
        &mut self,
        queue: &mut McsSeedQueue,
        query: &SearchTarget<'_>,
        targets: &[SearchTarget<'_>],
        target_tables: &[McsMatchTables],
        threshold_count: usize,
        max_bonds: usize,
        max_atoms: usize,
        params: &McsParameters,
        mut final_check: Option<
            &mut dyn FnMut(usize, &[(usize, usize)], &McsParameters) -> Result<bool, McsError>,
        >,
    ) -> Result<(), McsCandidateMatchError> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Seed.cpp :: Seed::grow ring/completion entry branches
        // RDKit✔️❌: void Seed::grow(MaximumCommonSubgraph &mcs) const {
        // RDKit✔️❌:   if (!canGrowBiggerThan(mcs.getMaxNumberBonds(), mcs.getMaxNumberAtoms())) {
        // RDKit✔️❌:     GrowingStage = NotSet;  // finished
        // RDKit✔️❌:     return;
        // RDKit✔️❌:   }
        // RDKit✔️❌:
        // RDKit✔️❌:   const auto &qmol = mcs.getQueryMolecule();
        // RDKit✔️❌:   boost::dynamic_bitset<> newAtomsSet(
        // RDKit✔️❌:       qmol.getNumAtoms());  // keep track of newly added atoms
        // RDKit✔️❌:
        // RDKit✔️❌:   if (0 == GrowingStage) {
        // RDKit✔️❌:     // 0. Fill out list of all directly connected outgoing bonds
        // RDKit✔️❌:     // non const method, multistage growing optimisation
        // RDKit✔️❌:     fillNewBonds(qmol, &mcs);
        // RDKit✔️❌:     if (NewBonds.empty()) {
        // RDKit✔️❌:       GrowingStage = NotSet;  // finished
        // RDKit✔️❌:       return;
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // END RDKIT CPP FUNCTION
        // Local complexity review: the constant-time remaining-bound gate and
        // one stage-zero frontier scan match the source. Complete-ring frontier
        // checks reuse Q128a, including its documented VF2 gap; immutable
        // parameters preserve the source's temporary fused-option rollback.
        if !self.can_grow_bigger_than(max_bonds, max_atoms) {
            self.growing_stage = u32::MAX;
            return Ok(());
        }

        if self.growing_stage == 0 {
            if let Some(check) = final_check.as_deref_mut() {
                self.fill_new_bonds(
                    query,
                    Some(params),
                    targets,
                    target_tables,
                    threshold_count,
                    Some(check),
                )?;
            } else {
                self.fill_new_bonds(
                    query,
                    Some(params),
                    targets,
                    target_tables,
                    threshold_count,
                    None,
                )?;
            }
            if self.new_bonds.is_empty() {
                self.growing_stage = u32::MAX;
                return Ok(());
            }
        }

        if final_check.is_some() {
            let mut check_without_params = |target_index: usize, mapping: &[(usize, usize)]| {
                final_check
                    .as_deref_mut()
                    .expect("checked final MCS callback presence")(
                    target_index, mapping, params
                )
            };
            self.grow_extensions(
                queue,
                query,
                targets,
                target_tables,
                threshold_count,
                max_bonds,
                max_atoms,
                params,
                Some(&mut check_without_params),
            )
        } else {
            self.grow_extensions(
                queue,
                query,
                targets,
                target_tables,
                threshold_count,
                max_bonds,
                max_atoms,
                params,
                None,
            )
        }
    }

    fn fill_new_bonds(
        &mut self,
        query: &SearchTarget<'_>,
        params: Option<&McsParameters>,
        targets: &[SearchTarget<'_>],
        target_tables: &[McsMatchTables],
        threshold_count: usize,
        mut final_check: Option<
            &mut dyn FnMut(usize, &[(usize, usize)], &McsParameters) -> Result<bool, McsError>,
        >,
    ) -> Result<(), McsCandidateMatchError> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Seed.cpp :: Seed::fillNewBonds
        // RDKit✔️❌:   auto excludedBonds = ExcludedBonds;
        // RDKit✔️❌:   const auto ri = qmol.getRingInfo();
        // RDKit✔️❌:   // all atoms added on previous growing only
        // RDKit✔️❌:   for (unsigned int srcAtomIdx = LastAddedAtomsBeginIdx;
        // RDKit✔️❌:        srcAtomIdx < getNumAtoms(); ++srcAtomIdx) {
        // RDKit✔️❌:     const auto atom = MoleculeFragment.Atoms.at(srcAtomIdx);
        // RDKit✔️❌:     for (const auto &nbri :
        // RDKit✔️❌:          boost::make_iterator_range(qmol.getAtomBonds(atom))) {
        // RDKit✔️❌:       const auto bond = qmol[nbri];
        // RDKit✔️❌:       const auto bi = bond->getIdx();
        // RDKit✔️❌:       // already in the seed or NewBonds list from another atom in a RING
        // RDKit✔️❌:       if (excludedBonds.test(bi)) {
        // RDKit✔️❌:         continue;
        // RDKit✔️❌:       }
        // RDKit✔️❌:       excludedBonds.set(bi);
        // RDKit✔️❌:       if (mcs && mcs->parameters().BondCompareParameters.CompleteRingsOnly &&
        // RDKit✔️❌:           ri->numBondRings(bi) == 1 &&
        // RDKit✔️❌:           !canAddAllNonFusedRingBondsConnectedToBond(*atom, *bond, *mcs)) {
        // RDKit✔️❌:         continue;
        // RDKit✔️❌:       }
        // RDKit✔️❌:       addNewBondFromAtom(*atom, *bond);
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // END RDKIT CPP FUNCTION
        // Local complexity review: the copied exclusion mask, last-frontier
        // scan and adjacency traversal match the source costs. Complete-ring
        // candidates delegate to Q128a and inherit its documented full-match
        // VF2 frontier-pruning gap, so the performance marker remains ❌.
        let mut excluded_bonds = self.excluded_bonds.clone();
        for source_seed_atom in
            self.last_added_atoms_begin_index..self.molecule_fragment.atoms.len()
        {
            let source_atom = self.molecule_fragment.atoms[source_seed_atom];
            if source_atom >= query.num_atoms() {
                return Err(McsError::AtomOutOfRange {
                    side: "query",
                    atom: source_atom,
                }
                .into());
            }
            for neighbor in query.adjacency().neighbors_of(source_atom) {
                let source_bond = neighbor.bond.index();
                let excluded_count = excluded_bonds.len();
                let excluded = excluded_bonds.get_mut(source_bond).ok_or(
                    McsError::SeedExcludedBondOutOfRange {
                        bond: source_bond,
                        count: excluded_count,
                    },
                )?;
                if *excluded {
                    continue;
                }
                *excluded = true;

                let mut accepted = true;
                if let Some(params) = params
                    && params.bond_compare_parameters.complete_rings_only
                {
                    let ring_info = query
                        .ring_info()
                        .ok_or(McsError::MissingRingInfo { side: "query" })?;
                    if ring_info.num_bond_rings(neighbor.bond) == 1 {
                        accepted = if let Some(check) = final_check.as_deref_mut() {
                            mcs_can_add_all_non_fused_ring_bonds_connected_to_bond(
                                self,
                                query,
                                source_atom,
                                source_bond,
                                targets,
                                target_tables,
                                threshold_count,
                                params,
                                Some(check),
                            )?
                        } else {
                            mcs_can_add_all_non_fused_ring_bonds_connected_to_bond(
                                self,
                                query,
                                source_atom,
                                source_bond,
                                targets,
                                target_tables,
                                threshold_count,
                                params,
                                None,
                            )?
                        };
                    }
                }
                if !accepted {
                    continue;
                }
                self.add_new_bond_from_atom(query, source_atom, source_bond)?;
            }
        }
        Ok(())
    }

    fn compute_remaining_size(&mut self, query: &SearchTarget<'_>) -> Result<(), McsError> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/Seed.cpp :: Seed::computeRemainingSize
        // RDKit✔️✔️:   RemainingBonds = RemainingAtoms = 0;
        // RDKit✔️✔️:   std::vector<unsigned int> end_atom_stack;
        // RDKit✔️✔️:   auto visitedBonds = ExcludedBonds;
        // RDKit✔️✔️:   boost::dynamic_bitset<> visitedAtoms(qmol.getNumAtoms());
        // RDKit✔️✔️:   std::for_each(
        // RDKit✔️✔️:       MoleculeFragment.Atoms.begin(), MoleculeFragment.Atoms.end(),
        // RDKit✔️✔️:       [&visitedAtoms](const auto &atom) { visitedAtoms.set(atom->getIdx()); });
        // RDKit✔️✔️:   // SDF all paths
        // RDKit✔️✔️:   // 1. direct neighbours
        // RDKit✔️✔️:   for (unsigned int seedAtomIdx = LastAddedAtomsBeginIdx;
        // RDKit✔️✔️:        seedAtomIdx < getNumAtoms(); ++seedAtomIdx) {
        // RDKit✔️✔️:     const auto atom = MoleculeFragment.Atoms.at(seedAtomIdx);
        // RDKit✔️✔️:     for (const auto &nbri :
        // RDKit✔️✔️:          boost::make_iterator_range(qmol.getAtomBonds(atom))) {
        // RDKit✔️✔️:       const auto bond = qmol[nbri];
        // RDKit✔️✔️:       if (!visitedBonds.test(bond->getIdx())) {
        // RDKit✔️✔️:         ++RemainingBonds;
        // RDKit✔️✔️:         visitedBonds.set(bond->getIdx());
        // RDKit✔️✔️:         unsigned int end_atom_idx = (atom == bond->getBeginAtom())
        // RDKit✔️✔️:                                         ? bond->getEndAtomIdx()
        // RDKit✔️✔️:                                         : bond->getBeginAtomIdx();
        // RDKit✔️✔️:         if (!visitedAtoms.test(end_atom_idx)) {  // check RING/CYCLE
        // RDKit✔️✔️:           ++RemainingAtoms;
        // RDKit✔️✔️:           visitedAtoms.set(end_atom_idx);
        // RDKit✔️✔️:           end_atom_stack.push_back(end_atom_idx);
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   // 2. go deep
        // RDKit✔️✔️:   while (!end_atom_stack.empty()) {
        // RDKit✔️✔️:     unsigned int ai = end_atom_stack.back();
        // RDKit✔️✔️:     end_atom_stack.pop_back();
        // RDKit✔️✔️:     const auto atom = qmol.getAtomWithIdx(ai);
        // RDKit✔️✔️:     for (const auto &nbri :
        // RDKit✔️✔️:          boost::make_iterator_range(qmol.getAtomBonds(atom))) {
        // RDKit✔️✔️:       const auto bond = qmol[nbri];
        // RDKit✔️✔️:       if (!visitedBonds.test(bond->getIdx())) {
        // RDKit✔️✔️:         ++RemainingBonds;
        // RDKit✔️✔️:         visitedBonds.set(bond->getIdx());
        // RDKit✔️✔️:         unsigned int end_atom_idx = (ai == bond->getBeginAtomIdx())
        // RDKit✔️✔️:                                         ? bond->getEndAtomIdx()
        // RDKit✔️✔️:                                         : bond->getBeginAtomIdx();
        // RDKit✔️✔️:         if (!visitedAtoms.test(end_atom_idx)) {  // check RING/CYCLE
        // RDKit✔️✔️:           ++RemainingAtoms;
        // RDKit✔️✔️:           visitedAtoms.set(end_atom_idx);
        // RDKit✔️✔️:           end_atom_stack.push_back(end_atom_idx);
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION
        // Local complexity review: exclusion/visited bit vectors and the DFS
        // stack are linear in query size; every reachable atom and bond row is
        // processed at most once, matching the source O(V+E) traversal.
        self.remaining_bonds = 0;
        self.remaining_atoms = 0;
        let mut end_atom_stack = Vec::new();
        let mut visited_bonds = self.excluded_bonds.clone();
        let mut visited_atoms = vec![false; query.num_atoms()];

        for &atom in &self.molecule_fragment.atoms {
            let visited = visited_atoms
                .get_mut(atom)
                .ok_or(McsError::AtomOutOfRange { side: "seed", atom })?;
            *visited = true;
        }

        for seed_atom_index in self.last_added_atoms_begin_index..self.molecule_fragment.atoms.len()
        {
            let atom = self.molecule_fragment.atoms[seed_atom_index];
            for neighbor in query.adjacency().neighbors_of(atom) {
                let bond = neighbor.bond.index();
                let visited =
                    visited_bonds
                        .get_mut(bond)
                        .ok_or(McsError::SeedExcludedBondOutOfRange {
                            bond,
                            count: self.excluded_bonds.len(),
                        })?;
                if !*visited {
                    self.remaining_bonds += 1;
                    *visited = true;
                    if !visited_atoms[neighbor.atom_index] {
                        self.remaining_atoms += 1;
                        visited_atoms[neighbor.atom_index] = true;
                        end_atom_stack.push(neighbor.atom_index);
                    }
                }
            }
        }

        while let Some(atom) = end_atom_stack.pop() {
            for neighbor in query.adjacency().neighbors_of(atom) {
                let bond = neighbor.bond.index();
                let visited =
                    visited_bonds
                        .get_mut(bond)
                        .ok_or(McsError::SeedExcludedBondOutOfRange {
                            bond,
                            count: self.excluded_bonds.len(),
                        })?;
                if !*visited {
                    self.remaining_bonds += 1;
                    *visited = true;
                    if !visited_atoms[neighbor.atom_index] {
                        self.remaining_atoms += 1;
                        visited_atoms[neighbor.atom_index] = true;
                        end_atom_stack.push(neighbor.atom_index);
                    }
                }
            }
        }
        Ok(())
    }
}

fn mcs_check_if_rings_are_closed(
    seed: &McsSeed,
    query: &SearchTarget<'_>,
    no_lone_ring_atoms: bool,
) -> Result<bool, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: checkIfRingsAreClosed
    // RDKit✔️✔️: bool checkIfRingsAreClosed(const Seed &fs, bool noLoneRingAtoms) {
    // RDKit✔️✔️:   if (fs.MoleculeFragment.Bonds.empty() && fs.MoleculeFragment.Atoms.empty()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const auto &om = fs.MoleculeFragment.Atoms.front()->getOwningMol();
    // RDKit✔️✔️:   const auto ri = om.getRingInfo();
    // RDKit✔️✔️:   if (!ri->numRings()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   boost::dynamic_bitset<> mcsBonds(om.getNumBonds());
    // RDKit✔️✔️:   boost::dynamic_bitset<> mcsNonFusedRings(ri->numRings());
    // RDKit✔️✔️:   boost::dynamic_bitset<> mcsFusedRings(ri->numRings());
    // RDKit✔️✔️:   for (const auto &bond : fs.MoleculeFragment.Bonds) {
    // RDKit✔️✔️:     auto bi = bond->getIdx();
    // RDKit✔️✔️:     mcsBonds.set(bi);
    // RDKit✔️✔️:     if (ri->numBondRings(bi) == 1) {
    // RDKit✔️✔️:       mcsNonFusedRings.set(ri->bondMembers(bi).front());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (unsigned int ringIdx = 0; ringIdx < mcsNonFusedRings.size(); ++ringIdx) {
    // RDKit✔️✔️:     if (!mcsNonFusedRings.test(ringIdx)) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (const auto &bi : ri->bondRings().at(ringIdx)) {
    // RDKit✔️✔️:       bool keepBond = false;
    // RDKit✔️✔️:       for (unsigned int memberOf : ri->bondMembers(bi)) {
    // RDKit✔️✔️:         if (memberOf == ringIdx) {
    // RDKit✔️✔️:           keepBond = true;
    // RDKit✔️✔️:         } else if (mcsNonFusedRings.test(memberOf)) {
    // RDKit✔️✔️:           keepBond = false;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (keepBond && !mcsBonds.test(bi)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (noLoneRingAtoms) {
    // RDKit✔️✔️:     for (const auto &atom : fs.MoleculeFragment.Atoms) {
    // RDKit✔️✔️:       auto ai = atom->getIdx();
    // RDKit✔️✔️:       const auto &ringIndices = ri->atomMembers(ai);
    // RDKit✔️✔️:       if (!ringIndices.empty() &&
    // RDKit✔️✔️:           !std::any_of(ringIndices.begin(), ringIndices.end(),
    // RDKit✔️✔️:                        [&mcsNonFusedRings](const auto &ringIdx) {
    // RDKit✔️✔️:                          return mcsNonFusedRings.test(ringIdx);
    // RDKit✔️✔️:                        })) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (mcsNonFusedRings.none()) {
    // RDKit✔️✔️:     for (const auto &bond : fs.MoleculeFragment.Bonds) {
    // RDKit✔️✔️:       auto bi = bond->getIdx();
    // RDKit✔️✔️:       if (ri->numBondRings(bi) > 1) {
    // RDKit✔️✔️:         for (auto ringIdx : ri->bondMembers(bi)) {
    // RDKit✔️✔️:           mcsFusedRings.set(ringIdx);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (mcsFusedRings.any()) {
    // RDKit✔️✔️:     for (unsigned int ringIdx = 0; ringIdx < mcsFusedRings.size(); ++ringIdx) {
    // RDKit✔️✔️:       if (!mcsFusedRings.test(ringIdx)) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       const auto &ringBondIndices = ri->bondRings().at(ringIdx);
    // RDKit✔️✔️:       if (std::all_of(
    // RDKit✔️✔️:               ringBondIndices.begin(), ringBondIndices.end(),
    // RDKit✔️✔️:               [&mcsBonds](const auto &bi) { return mcsBonds.test(bi); })) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    // Local complexity review: the three source-sized bit sets are represented
    // by three boolean vectors. Selected bonds, marked rings and their member
    // rows are traversed in the same order and with the same asymptotic cost.
    if seed.molecule_fragment.bonds.is_empty() && seed.molecule_fragment.atoms.is_empty() {
        return Ok(true);
    }
    let ring_info = query
        .ring_info()
        .ok_or(McsError::MissingRingInfo { side: "query" })?;
    if ring_info.num_rings() == 0 {
        return Ok(true);
    }

    let mut mcs_bonds = vec![false; query.num_bonds()];
    let mut mcs_non_fused_rings = vec![false; ring_info.num_rings()];
    let mut mcs_fused_rings = vec![false; ring_info.num_rings()];
    for &source_bond in &seed.molecule_fragment.bonds {
        let bond = query
            .bonds()
            .get(source_bond)
            .ok_or(McsError::BondOutOfRange {
                side: "query",
                bond: source_bond,
            })?;
        mcs_bonds[source_bond] = true;
        if ring_info.num_bond_rings(bond.id()) == 1 {
            let ring = ring_info
                .bond_members(bond.id())
                .first()
                .copied()
                .ok_or(McsCandidateMatchError::RingMembershipMissing { bond: source_bond })?;
            *mcs_non_fused_rings
                .get_mut(ring)
                .ok_or(McsCandidateMatchError::RingOutOfRange { ring })? = true;
        }
    }

    for ring in 0..mcs_non_fused_rings.len() {
        if !mcs_non_fused_rings[ring] {
            continue;
        }
        let ring_bonds = ring_info
            .bond_rings()
            .get(ring)
            .ok_or(McsCandidateMatchError::RingOutOfRange { ring })?;
        for ring_bond in ring_bonds {
            let source_bond = ring_bond.index();
            let bond = query
                .bonds()
                .get(source_bond)
                .ok_or(McsError::BondOutOfRange {
                    side: "query",
                    bond: source_bond,
                })?;
            let mut keep_bond = false;
            for &member_of in ring_info.bond_members(bond.id()) {
                if member_of == ring {
                    keep_bond = true;
                } else if mcs_non_fused_rings
                    .get(member_of)
                    .copied()
                    .ok_or(McsCandidateMatchError::RingOutOfRange { ring: member_of })?
                {
                    keep_bond = false;
                    break;
                }
            }
            if keep_bond && !mcs_bonds[source_bond] {
                return Ok(false);
            }
        }
    }

    if no_lone_ring_atoms {
        for &source_atom in &seed.molecule_fragment.atoms {
            let atom = query
                .atoms()
                .get(source_atom)
                .ok_or(McsError::AtomOutOfRange {
                    side: "query",
                    atom: source_atom,
                })?;
            let memberships = ring_info.atom_members(atom.id());
            if !memberships.is_empty()
                && !memberships
                    .iter()
                    .copied()
                    .any(|ring| mcs_non_fused_rings.get(ring).copied().unwrap_or(false))
            {
                return Ok(false);
            }
        }
    }

    if !mcs_non_fused_rings.iter().any(|selected| *selected) {
        for &source_bond in &seed.molecule_fragment.bonds {
            let bond = &query.bonds()[source_bond];
            if ring_info.num_bond_rings(bond.id()) > 1 {
                for &ring in ring_info.bond_members(bond.id()) {
                    *mcs_fused_rings
                        .get_mut(ring)
                        .ok_or(McsCandidateMatchError::RingOutOfRange { ring })? = true;
                }
            }
        }
    }
    if mcs_fused_rings.iter().any(|selected| *selected) {
        for (ring, selected) in mcs_fused_rings.iter().copied().enumerate() {
            if !selected {
                continue;
            }
            let ring_bonds = ring_info
                .bond_rings()
                .get(ring)
                .ok_or(McsCandidateMatchError::RingOutOfRange { ring })?;
            if ring_bonds
                .iter()
                .all(|bond| mcs_bonds.get(bond.index()).copied().unwrap_or(false))
            {
                return Ok(true);
            }
        }
        return Ok(false);
    }
    Ok(true)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct McsResultContext {
    query_input: usize,
    target_inputs: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct McsRetainedResult {
    fragment: McsMoleculeFragment,
    context: McsResultContext,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct McsSearchState {
    // BEGIN RDKIT CPP TYPE: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.h :: MCS and search-global members
    // RDKit✔️✔️:     std::vector<const Atom *> Atoms;
    // RDKit✔️✔️:     std::vector<const Bond *> Bonds;
    // RDKit✔️✔️:     const ROMol *QueryMolecule;
    // RDKit✔️✔️:     std::vector<Target> Targets;
    // RDKit✔️✔️:   MCS McsIdx;
    // RDKit✔️✔️:   std::map<std::vector<unsigned int>, MCS> DegenerateMcsMap;
    // END RDKIT CPP TYPE
    // Stable input indices retain the borrowed query and ordered target
    // identity without cloning a molecule or any match table.
    best: McsMoleculeFragment,
    best_context: Option<McsResultContext>,
    degenerate: BTreeMap<Vec<usize>, McsRetainedResult>,
    progress: McsProgressData,
}

impl McsSearchState {
    fn clear_best_for_zero_bond_fallback(&mut self) {
        // RDKit✔️✔️:       McsIdx = MCS();      // clear
        // The source does not clear DegenerateMcsMap at this point.
        self.best = McsMoleculeFragment::default();
        self.best_context = None;
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct McsGrowSeedsOutcome {
    mcs_found: bool,
    canceled: bool,
}

#[allow(clippy::too_many_arguments)]
fn mcs_grow_seeds(
    queue: &mut McsSeedQueue,
    state: &mut McsSearchState,
    context: &McsResultContext,
    query: &SearchTarget<'_>,
    targets: &[SearchTarget<'_>],
    target_tables: &[McsMatchTables],
    threshold_count: usize,
    query_matched_bonds: usize,
    params: &McsParameters,
    start: u64,
    now: &mut dyn FnMut() -> u64,
    mut final_check: Option<
        &mut dyn FnMut(usize, &[(usize, usize)], &McsParameters) -> Result<bool, McsError>,
    >,
    mut should_accept_mcs: Option<
        &mut dyn FnMut(&McsMoleculeFragment, &McsParameters) -> Result<bool, McsError>,
    >,
    mut progress_callback: Option<
        &mut dyn FnMut(&McsProgressData, &McsParameters) -> Result<bool, McsProgressError>,
    >,
) -> Result<McsGrowSeedsOutcome, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: MaximumCommonSubgraph::growSeeds
    // RDKit✔️❌: bool MaximumCommonSubgraph::growSeeds() {
    // RDKit✔️❌:   bool mcsFound = false;
    // RDKit✔️❌:   bool canceled = false;
    // RDKit✔️❌:   // Find MCS -- SDF Seed growing OPTIMISATION (it works in 3 times
    // RDKit✔️❌:   // faster)
    // RDKit✔️❌:   while (!Seeds.empty()) {
    // RDKit✔️❌:     if (getMaxNumberBonds() == QueryMoleculeMatchedBonds) {  // MCS == Query
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     auto si = Seeds.begin();
    // RDKit✔️❌:     si->grow(*this);
    // RDKit✔️❌:     {
    // RDKit✔️❌:       const Seed &fs = Seeds.front();
    // RDKit✔️❌:       // bigger substructure found
    // RDKit✔️❌:       if (fs.CopyComplete) {
    // RDKit✔️❌:         bool possibleMCS = false;
    // RDKit✔️❌:         if (!Parameters.MaximizeBonds) {
    // RDKit✔️❌:           possibleMCS = (fs.getNumAtoms() > getMaxNumberAtoms() ||
    // RDKit✔️❌:                          (fs.getNumAtoms() == getMaxNumberAtoms() &&
    // RDKit✔️❌:                           fs.getNumBonds() > getMaxNumberBonds()));
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           possibleMCS = (fs.getNumBonds() > getMaxNumberBonds() ||
    // RDKit✔️❌:                          (fs.getNumBonds() == getMaxNumberBonds() &&
    // RDKit✔️❌:                           fs.getNumAtoms() > getMaxNumberAtoms()));
    // RDKit✔️❌:         }
    // RDKit✔️❌:         bool isDegenerateMCS = (fs.getNumBonds() == getMaxNumberBonds() &&
    // RDKit✔️❌:                                 fs.getNumAtoms() == getMaxNumberAtoms());
    // RDKit✔️❌:         if (!possibleMCS && Parameters.StoreAll) {
    // RDKit✔️❌:           possibleMCS = isDegenerateMCS;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         // #945: test here to see if the MCS actually has all rings closed
    // RDKit✔️❌:         if (possibleMCS && Parameters.BondCompareParameters.CompleteRingsOnly) {
    // RDKit✔️❌:           possibleMCS = checkIfRingsAreClosed(
    // RDKit✔️❌:               fs, Parameters.AtomCompareParameters.CompleteRingsOnly);
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (possibleMCS) {
    // RDKit✔️❌:           possibleMCS = checkIfShouldAcceptMCS(
    // RDKit✔️❌:               fs.MoleculeFragment, *QueryMolecule, Targets, Parameters);
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (possibleMCS) {
    // RDKit✔️❌:           mcsFound = true;
    // RDKit✔️❌:           McsIdx.Atoms = fs.MoleculeFragment.Atoms;
    // RDKit✔️❌:           McsIdx.Bonds = fs.MoleculeFragment.Bonds;
    // RDKit✔️❌:           if (Parameters.StoreAll) {
    // RDKit✔️❌:             if (!isDegenerateMCS) {
    // RDKit✔️❌:               DegenerateMcsMap.clear();
    // RDKit✔️❌:             }
    // RDKit✔️❌:             std::vector<unsigned int> key(McsIdx.Bonds.size());
    // RDKit✔️❌:             std::transform(McsIdx.Bonds.begin(), McsIdx.Bonds.end(),
    // RDKit✔️❌:                            key.begin(),
    // RDKit✔️❌:                            [](const auto bond) { return bond->getIdx(); });
    // RDKit✔️❌:             std::sort(key.begin(), key.end());
    // RDKit✔️❌:             MCS value(McsIdx);
    // RDKit✔️❌:             value.QueryMolecule = QueryMolecule;
    // RDKit✔️❌:             value.Targets = Targets;
    // RDKit✔️❌:             DegenerateMcsMap[key] = value;
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (NotSet == si->GrowingStage) {  // finished
    // RDKit✔️❌:       Seeds.erase(si);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (Parameters.ProgressCallback) {
    // RDKit✔️❌:       Stat.NumAtoms = getMaxNumberAtoms();
    // RDKit✔️❌:       Stat.NumBonds = getMaxNumberBonds();
    // RDKit✔️❌:       if (!Parameters.ProgressCallback(Stat, Parameters,
    // RDKit✔️❌:                                        Parameters.ProgressCallbackUserData)) {
    // RDKit✔️❌:         canceled = true;
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return !canceled;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION
    // Local complexity review: queue insertion and active-seed reinsertion
    // retain the source linear ordered-list scans. Matching inherits the
    // completed full-candidate helper's documented VF2 frontier-pruning gap.
    let mut outcome = McsGrowSeedsOutcome::default();
    let work = (|| -> Result<(), McsCandidateMatchError> {
        while !queue.seeds.is_empty() {
            if state.best.bonds.len() == query_matched_bonds {
                break;
            }

            let mut active_seed = queue.seeds.remove(0);
            let grow_result = if let Some(check) = final_check.as_deref_mut() {
                active_seed.grow(
                    queue,
                    query,
                    targets,
                    target_tables,
                    threshold_count,
                    state.best.bonds.len(),
                    state.best.atoms.len(),
                    params,
                    Some(check),
                )
            } else {
                active_seed.grow(
                    queue,
                    query,
                    targets,
                    target_tables,
                    threshold_count,
                    state.best.bonds.len(),
                    state.best.atoms.len(),
                    params,
                    None,
                )
            };

            let active_bonds = active_seed.molecule_fragment.bonds.len();
            let active_position = queue
                .seeds
                .iter()
                .position(|seed| seed.molecule_fragment.bonds.len() <= active_bonds)
                .unwrap_or(queue.seeds.len());
            let active_finished = active_seed.growing_stage == u32::MAX;
            queue.seeds.insert(active_position, active_seed);
            grow_result?;

            let front = &queue.seeds[0];
            if front.copy_complete {
                let atom_count = front.molecule_fragment.atoms.len();
                let bond_count = front.molecule_fragment.bonds.len();
                let mut possible_mcs = if params.maximize_bonds {
                    bond_count > state.best.bonds.len()
                        || (bond_count == state.best.bonds.len()
                            && atom_count > state.best.atoms.len())
                } else {
                    atom_count > state.best.atoms.len()
                        || (atom_count == state.best.atoms.len()
                            && bond_count > state.best.bonds.len())
                };
                let is_degenerate_mcs =
                    bond_count == state.best.bonds.len() && atom_count == state.best.atoms.len();
                if !possible_mcs && params.store_all {
                    possible_mcs = is_degenerate_mcs;
                }
                if possible_mcs && params.bond_compare_parameters.complete_rings_only {
                    possible_mcs = mcs_check_if_rings_are_closed(
                        front,
                        query,
                        params.atom_compare_parameters.complete_rings_only,
                    )?;
                }
                if possible_mcs {
                    possible_mcs = if let Some(check) = should_accept_mcs.as_deref_mut() {
                        check(&front.molecule_fragment, params)?
                    } else {
                        true
                    };
                }
                if possible_mcs {
                    outcome.mcs_found = true;
                    state.best = front.molecule_fragment.clone();
                    if params.store_all {
                        if !is_degenerate_mcs {
                            state.degenerate.clear();
                        }
                        let mut key = state.best.bonds.clone();
                        key.sort_unstable();
                        state.degenerate.insert(
                            key,
                            McsRetainedResult {
                                fragment: state.best.clone(),
                                context: context.clone(),
                            },
                        );
                    }
                }
            }

            if active_finished {
                queue.seeds.remove(active_position);
            }
            let observed_now = if progress_callback.is_some() {
                start
            } else {
                now()
            };
            let canceled = if let Some(callback) = progress_callback.as_deref_mut() {
                mcs_progress_after_seed(
                    params,
                    start,
                    observed_now,
                    &state.best,
                    &mut state.progress,
                    Some(callback),
                )?
                .canceled
            } else {
                mcs_progress_after_seed(
                    params,
                    start,
                    observed_now,
                    &state.best,
                    &mut state.progress,
                    None,
                )?
                .canceled
            };
            if canceled {
                outcome.canceled = true;
                break;
            }
        }
        Ok(())
    })();
    if outcome.mcs_found {
        // RDKit✔️✔️:   if (mcsFound) {  // postponed copy of current set of molecules for
        // RDKit✔️✔️:   // threshold < 1.
        // RDKit✔️✔️:     McsIdx.QueryMolecule = QueryMolecule;
        // RDKit✔️✔️:     McsIdx.Targets = Targets;
        // RDKit✔️✔️:   }
        state.best_context = Some(context.clone());
    }
    work?;
    Ok(outcome)
}

#[allow(clippy::too_many_arguments)]
fn mcs_build_result_query_graph(
    fragment: &McsMoleculeFragment,
    query: &SearchTarget<'_>,
    targets: &[SearchTarget<'_>],
    target_tables: &[McsMatchTables],
    params: &McsParameters,
    mut final_check: Option<
        &mut dyn FnMut(usize, &[(usize, usize)], &McsParameters) -> Result<bool, McsError>,
    >,
) -> Result<QueryGraph, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: MaximumCommonSubgraph::generateResultSMARTSAndQueryMol query construction
    // RDKit✔️❌:   // match the result MCS with all targets to check if it is exact match
    // RDKit✔️❌:   // or template
    // RDKit✔️❌:   Seed seed;  // result MCS
    // RDKit✔️❌:   seed.setStoreAllDegenerateMCS(Parameters.StoreAll);
    // RDKit✔️❌:   seed.ExcludedBonds.resize(mcsIdx.QueryMolecule->getNumBonds(), false);
    // RDKit✔️❌:   std::vector<AtomMatchSet> atomMatchResult(mcsIdx.Targets.size());
    // RDKit✔️❌:   std::vector<unsigned int> atomIdxMap(mcsIdx.QueryMolecule->getNumAtoms());
    // RDKit✔️❌:   std::vector<std::map<unsigned int, const Bond *>> bondMatchSet(
    // RDKit✔️❌:       mcsIdx.Bonds.size());  // key is unique BondType
    // RDKit✔️❌:   std::vector<std::map<unsigned int, const Atom *>> atomMatchSet(
    // RDKit✔️❌:       mcsIdx.Atoms.size());  // key is unique atomic number
    // RDKit✔️❌:
    // RDKit✔️❌:   for (const auto &atom : mcsIdx.Atoms) {
    // RDKit✔️❌:     atomIdxMap[atom->getIdx()] = seed.getNumAtoms();
    // RDKit✔️❌:     seed.addAtom(atom);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (const auto &bond : mcsIdx.Bonds) {
    // RDKit✔️❌:     seed.addBond(bond);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   std::vector<unsigned int> matchedTargetIndices;
    // RDKit✔️❌:   if (!mcsIdx.Bonds.empty()) {
    // RDKit✔️❌:     for (const auto &tag : mcsIdx.Targets) {
    // RDKit✔️❌:       match_V_t match;  // THERE IS NO Bonds match INFO !!!!
    // RDKit✔️❌:       bool target_matched = SubstructMatchCustomTable(
    // RDKit✔️❌:           tag.Topology, *tag.Molecule, seed.Topology, *QueryMolecule,
    // RDKit✔️❌:           tag.AtomMatchTable, tag.BondMatchTable, &Parameters, &match);
    // RDKit✔️❌:       if (!target_matched) {
    // RDKit✔️❌:         continue;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       unsigned int itarget = &tag - &mcsIdx.Targets.front();
    // RDKit✔️❌:       matchedTargetIndices.push_back(itarget);
    // RDKit✔️❌:       atomMatchResult.at(itarget).resize(seed.getNumAtoms());
    // RDKit✔️❌:       for (const auto &m : match) {
    // RDKit✔️❌:         const auto ai = m.first;  // SeedAtomIdx
    // RDKit✔️❌:         atomMatchResult.at(itarget).at(ai).QueryAtomIdx =
    // RDKit✔️❌:             seed.Topology[m.first];
    // RDKit✔️❌:         atomMatchResult.at(itarget).at(ai).TargetAtomIdx =
    // RDKit✔️❌:             tag.Topology[m.second];
    // RDKit✔️❌:         const auto ta = tag.Molecule->getAtomWithIdx(tag.Topology[m.second]);
    // RDKit✔️❌:         if (ta->getAtomicNum() !=
    // RDKit✔️❌:             seed.MoleculeFragment.Atoms.at(ai)->getAtomicNum()) {
    // RDKit✔️❌:           atomMatchSet[ai][ta->getAtomicNum()] = ta;  // add
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       // AND BUILD BOND MATCH INFO
    // RDKit✔️❌:       for (const auto &bond : mcsIdx.Bonds) {
    // RDKit✔️❌:         const auto bi = &bond - &mcsIdx.Bonds.front();
    // RDKit✔️❌:         const auto i = atomIdxMap.at(bond->getBeginAtomIdx());
    // RDKit✔️❌:         const auto j = atomIdxMap.at(bond->getEndAtomIdx());
    // RDKit✔️❌:         const auto ti = atomMatchResult.at(itarget).at(i).TargetAtomIdx;
    // RDKit✔️❌:         const auto tj = atomMatchResult.at(itarget).at(j).TargetAtomIdx;
    // RDKit✔️❌:         const auto tb = tag.Molecule->getBondBetweenAtoms(ti, tj);
    // RDKit✔️❌:         if (tb && bond->getBondType() != tb->getBondType()) {
    // RDKit✔️❌:           bondMatchSet[bi][tb->getBondType()] = tb;  // add
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // create molecule from MCS for MolToSmarts()
    // RDKit✔️❌:   auto mol = new RWMol();
    // RDKit✔️❌:   ROMOL_SPTR molSptr(mol);
    // RDKit✔️❌:   const auto ri = mcsIdx.QueryMolecule->getRingInfo();
    // RDKit✔️❌:   boost::dynamic_bitset<> mcsRingIsComplete;
    // RDKit✔️❌:   bool needAtomRingQueries =
    // RDKit✔️❌:       (Parameters.AtomCompareParameters.RingMatchesRingOnly ||
    // RDKit✔️❌:        Parameters.BondCompareParameters.MatchFusedRingsStrict);
    // RDKit✔️❌:   if (needAtomRingQueries) {
    // RDKit✔️❌:     mcsRingIsComplete.resize(ri->numRings(), true);
    // RDKit✔️❌:     boost::dynamic_bitset<> queryBondInMcs(mcsIdx.QueryMolecule->getNumBonds());
    // RDKit✔️❌:     for (const auto &bond : mcsIdx.Bonds) {
    // RDKit✔️❌:       queryBondInMcs.set(bond->getIdx());
    // RDKit✔️❌:     }
    // RDKit✔️❌:     const auto &bondRings = ri->bondRings();
    // RDKit✔️❌:     for (const auto &bondRing : bondRings) {
    // RDKit✔️❌:       auto ringIdx = &bondRing - &bondRings.front();
    // RDKit✔️❌:       for (const auto &bondIdx : bondRing) {
    // RDKit✔️❌:         if (!queryBondInMcs.test(bondIdx)) {
    // RDKit✔️❌:           mcsRingIsComplete.reset(ringIdx);
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (const auto &atom : mcsIdx.Atoms) {
    // RDKit✔️❌:     auto queryAtomIdx = atom->getIdx();
    // RDKit✔️❌:     auto numAtomRings = ri->numAtomRings(queryAtomIdx);
    // RDKit✔️❌:     QueryAtom a;
    // RDKit✔️❌:     const auto ai = &atom - &mcsIdx.Atoms.front();
    // RDKit✔️❌:     if (Parameters.AtomTyper == MCSAtomCompareIsotopes ||
    // RDKit✔️❌:         Parameters.AtomCompareParameters
    // RDKit✔️❌:             .MatchIsotope) {  // do '[0*]-[0*]-[13*]' for CC[13NH2]
    // RDKit✔️❌:       a.setQuery(makeAtomIsotopeQuery(static_cast<int>(atom->getIsotope())));
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       // generate [#6] instead of C or c !
    // RDKit✔️❌:       a.setQuery(makeAtomNumQuery(atom->getAtomicNum()));
    // RDKit✔️❌:       // for all atomMatchSet[ai] items add atom query to template like
    // RDKit✔️❌:       // [#6,#17,#9, ... ]
    // RDKit✔️❌:       for (const auto &am : atomMatchSet.at(ai)) {
    // RDKit✔️❌:         a.expandQuery(makeAtomNumQuery(am.second->getAtomicNum()),
    // RDKit✔️❌:                       Queries::COMPOSITE_OR);
    // RDKit✔️❌:         if (Parameters.AtomCompareParameters.MatchChiralTag &&
    // RDKit✔️❌:             (am.second->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️❌:              am.second->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW)) {
    // RDKit✔️❌:           a.setChiralTag(am.second->getChiralTag());
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (needAtomRingQueries) {
    // RDKit✔️❌:       const auto &ringIndicesAtomIsMemberOf = ri->atomMembers(queryAtomIdx);
    // RDKit✔️❌:       auto numCompleteRings = std::count_if(
    // RDKit✔️❌:           ringIndicesAtomIsMemberOf.begin(), ringIndicesAtomIsMemberOf.end(),
    // RDKit✔️❌:           [&mcsRingIsComplete](const auto &ringIdx) {
    // RDKit✔️❌:             return mcsRingIsComplete.test(ringIdx);
    // RDKit✔️❌:           });
    // RDKit✔️❌:       if (Parameters.AtomCompareParameters.RingMatchesRingOnly &&
    // RDKit✔️❌:           !numCompleteRings) {
    // RDKit✔️❌:         auto q = makeAtomInRingQuery();
    // RDKit✔️❌:         q->setNegation(!numAtomRings);
    // RDKit✔️❌:         a.expandQuery(q, Queries::COMPOSITE_AND, true);
    // RDKit✔️❌:       } else if (Parameters.BondCompareParameters.MatchFusedRingsStrict &&
    // RDKit✔️❌:                  numAtomRings == 1 && numCompleteRings == 1) {
    // RDKit✔️❌:         auto ringSize =
    // RDKit✔️❌:             ri->atomRings().at(ringIndicesAtomIsMemberOf.front()).size();
    // RDKit✔️❌:         auto q = new ATOM_OR_QUERY;
    // RDKit✔️❌:         q->setDescription("AtomOr");
    // RDKit✔️❌:         q->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(
    // RDKit✔️❌:             makeAtomMinRingSizeQuery(ringSize)));
    // RDKit✔️❌:         auto q2 = makeAtomInNRingsQuery(1);
    // RDKit✔️❌:         q2->setNegation(true);
    // RDKit✔️❌:         q->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(q2));
    // RDKit✔️❌:         a.expandQuery(q, Queries::COMPOSITE_AND, true);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     mol->addAtom(&a, true, false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (const auto &bond : mcsIdx.Bonds) {
    // RDKit✔️❌:     QueryBond b;
    // RDKit✔️❌:     const auto bi = &bond - &mcsIdx.Bonds.front();
    // RDKit✔️❌:     const auto beginAtomIdx = atomIdxMap.at(bond->getBeginAtomIdx());
    // RDKit✔️❌:     const auto endAtomIdx = atomIdxMap.at(bond->getEndAtomIdx());
    // RDKit✔️❌:     b.setBeginAtomIdx(beginAtomIdx);
    // RDKit✔️❌:     b.setEndAtomIdx(endAtomIdx);
    // RDKit✔️❌:     b.setQuery(makeBondOrderEqualsQuery(bond->getBondType()));
    // RDKit✔️❌:     // add OR template if need
    // RDKit✔️❌:     for (const auto &bm : bondMatchSet.at(bi)) {
    // RDKit✔️❌:       b.expandQuery(makeBondOrderEqualsQuery(bm.second->getBondType()),
    // RDKit✔️❌:                     Queries::COMPOSITE_OR);
    // RDKit✔️❌:       if (Parameters.BondCompareParameters.MatchStereo &&
    // RDKit✔️❌:           bm.second->getStereo() > Bond::STEREOANY) {
    // RDKit✔️❌:         b.setStereo(bm.second->getStereo());
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (Parameters.BondCompareParameters.RingMatchesRingOnly ||
    // RDKit✔️❌:         Parameters.BondCompareParameters.MatchFusedRingsStrict) {
    // RDKit✔️❌:       const auto numBondRings = ri->numBondRings(bond->getIdx());
    // RDKit✔️❌:       auto q = makeBondIsInRingQuery();
    // RDKit✔️❌:       q->setNegation(!numBondRings);
    // RDKit✔️❌:       b.expandQuery(q, Queries::COMPOSITE_AND, true);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     mol->addBond(&b, false);
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION
    // Local complexity review: result rematching uses the completed full
    // matcher and inherits its documented VF2 pruning gap. Source-row maps,
    // ordered alternative maps, ring scans and graph construction otherwise
    // retain the source asymptotic costs and allocation shape.
    let mut seed = McsSeed {
        store_all_degenerate_mcs: params.store_all,
        excluded_bonds: vec![false; query.num_bonds()],
        ..McsSeed::default()
    };
    let mut atom_index_map = vec![usize::MAX; query.num_atoms()];
    for &source_atom in &fragment.atoms {
        if source_atom >= query.num_atoms() {
            return Err(McsError::AtomOutOfRange {
                side: "MCS result query",
                atom: source_atom,
            }
            .into());
        }
        atom_index_map[source_atom] = seed.molecule_fragment.atoms.len();
        seed.add_atom(source_atom);
    }
    for &source_bond in &fragment.bonds {
        seed.add_bond(query, source_bond)?;
    }

    let mut atom_match_sets = vec![BTreeMap::<u8, (usize, usize)>::new(); fragment.atoms.len()];
    let mut bond_match_sets =
        vec![BTreeMap::<i64, (BondOrder, usize, usize)>::new(); fragment.bonds.len()];
    if !fragment.bonds.is_empty() {
        if targets.len() != target_tables.len() {
            return Err(McsCandidateMatchError::TargetTableCount {
                targets: targets.len(),
                tables: target_tables.len(),
            });
        }
        for (target_index, (target, tables)) in targets.iter().zip(target_tables).enumerate() {
            let mut accept_mapping = |mapping: &[(usize, usize)]| {
                let has_user_hook = final_check.is_some();
                let mut user_hook = || {
                    final_check
                        .as_deref_mut()
                        .expect("present result final hook")(
                        target_index, mapping, params
                    )
                    .map_err(McsCandidateMatchError::from)
                };
                let user = if has_user_hook {
                    Some(&mut user_hook as &mut dyn FnMut() -> Result<bool, McsCandidateMatchError>)
                } else {
                    None
                };
                mcs_final_mapping_accept(&seed.topology, query, target, params, mapping, user)
            };
            let mapping = mcs_find_full_mapping(&seed, target, tables, &mut accept_mapping)?;
            let Some(mapping) = mapping else {
                continue;
            };
            let mut target_atoms = vec![usize::MAX; fragment.atoms.len()];
            for (seed_atom, target_atom) in mapping {
                target_atoms[seed_atom] = target_atom;
                let source_atom = seed.topology.source_atoms[seed_atom];
                let source_number = query.atoms()[source_atom].atomic_number();
                let target_number = target.atoms()[target_atom].atomic_number();
                if target_number != source_number {
                    atom_match_sets[seed_atom].insert(target_number, (target_index, target_atom));
                }
            }
            for (result_bond, &source_bond) in fragment.bonds.iter().enumerate() {
                let bond = query
                    .bonds()
                    .get(source_bond)
                    .ok_or(McsError::BondOutOfRange {
                        side: "MCS result query",
                        bond: source_bond,
                    })?;
                let begin = atom_index_map[bond.begin().index()];
                let end = atom_index_map[bond.end().index()];
                let target_begin = target_atoms[begin];
                let target_end = target_atoms[end];
                if let Some(target_bond) = mcs_target_bond_between(target, target_begin, target_end)
                {
                    let target_order = target.bonds()[target_bond].order();
                    if target_order != bond.order() {
                        bond_match_sets[result_bond].insert(
                            target_order.rdkit_code(),
                            (target_order, target_index, target_bond),
                        );
                    }
                }
            }
        }
    }

    let need_atom_ring_queries = params.atom_compare_parameters.ring_matches_ring_only
        || params.bond_compare_parameters.match_fused_rings_strict;
    let need_bond_ring_queries = params.bond_compare_parameters.ring_matches_ring_only
        || params.bond_compare_parameters.match_fused_rings_strict;
    let ring_info = if need_atom_ring_queries || need_bond_ring_queries {
        Some(
            query
                .ring_info()
                .ok_or(McsError::MissingRingInfo { side: "query" })?,
        )
    } else {
        None
    };
    let mut mcs_ring_is_complete = Vec::new();
    if need_atom_ring_queries {
        let rings = ring_info.expect("required ring information is present");
        let mut query_bond_in_mcs = vec![false; query.num_bonds()];
        for &source_bond in &fragment.bonds {
            query_bond_in_mcs[source_bond] = true;
        }
        mcs_ring_is_complete = rings
            .bond_rings()
            .iter()
            .map(|ring| {
                ring.iter().all(|bond| {
                    query_bond_in_mcs
                        .get(bond.index())
                        .copied()
                        .unwrap_or(false)
                })
            })
            .collect();
    }

    let mut atoms = Vec::with_capacity(fragment.atoms.len());
    for (result_atom, &source_atom) in fragment.atoms.iter().enumerate() {
        let source = &query.atoms()[source_atom];
        let mut chiral_tag = ChiralTag::Unspecified;
        let mut predicate = if params.atom_comparator == AtomComparator::AtomCompareIsotopes
            || params.atom_compare_parameters.match_isotope
        {
            make_atom_isotope_query(i32::from(source.isotope().unwrap_or(0)))
        } else {
            let mut predicate = make_atom_num_query(source.atomic_number());
            for &(target_index, target_atom) in atom_match_sets[result_atom].values() {
                let target = &targets[target_index].atoms()[target_atom];
                query_atom_expand_query(
                    &mut predicate,
                    make_atom_num_query(target.atomic_number()),
                    CompositeQueryType::Or,
                    false,
                );
                if params.atom_compare_parameters.match_chiral_tag
                    && matches!(
                        target.chiral_tag(),
                        ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
                    )
                {
                    chiral_tag = target.chiral_tag();
                }
            }
            predicate
        };

        if need_atom_ring_queries {
            let rings = ring_info.expect("required ring information is present");
            let memberships = rings.atom_members(source.id());
            let num_atom_rings = memberships.len();
            let num_complete_rings = memberships
                .iter()
                .filter(|ring| mcs_ring_is_complete.get(**ring).copied().unwrap_or(false))
                .count();
            if params.atom_compare_parameters.ring_matches_ring_only && num_complete_rings == 0 {
                let mut ring_query = make_atom_in_ring_query();
                ring_query.set_negation(num_atom_rings == 0);
                query_atom_expand_query(&mut predicate, ring_query, CompositeQueryType::And, true);
            } else if params.bond_compare_parameters.match_fused_rings_strict
                && num_atom_rings == 1
                && num_complete_rings == 1
            {
                let ring = memberships[0];
                let ring_size = rings
                    .atom_rings()
                    .get(ring)
                    .ok_or(McsCandidateMatchError::RingOutOfRange { ring })?
                    .len();
                let ring_size = i32::try_from(ring_size).map_err(|_| {
                    McsCandidateMatchError::ResultValueOutOfRange {
                        kind: "ring size",
                        value: ring_size,
                    }
                })?;
                let strict_ring_query = QueryNode::or(vec![
                    make_atom_min_ring_size_query(ring_size),
                    QueryNode::not(QueryNode::predicate(AtomQueryPredicate::NumAtomRings(1))),
                ]);
                query_atom_expand_query(
                    &mut predicate,
                    strict_ring_query,
                    CompositeQueryType::And,
                    true,
                );
            }
        }

        let mut atom = QueryAtom::from_identity_parts(
            AtomId::new(result_atom),
            QueryAtomIdentity::from_atomic_number(0),
            predicate,
        );
        atom.set_chiral_tag(chiral_tag);
        atoms.push(atom);
    }

    let mut bonds = Vec::with_capacity(fragment.bonds.len());
    for (result_bond, &source_bond) in fragment.bonds.iter().enumerate() {
        let source = &query.bonds()[source_bond];
        let begin = atom_index_map[source.begin().index()];
        let end = atom_index_map[source.end().index()];
        let mut predicate = make_bond_order_equals_query(source.order());
        let mut stereo = BondStereo::None;
        for &(target_order, target_index, target_bond) in bond_match_sets[result_bond].values() {
            query_bond_expand_query(
                &mut predicate,
                make_bond_order_equals_query(target_order),
                CompositeQueryType::Or,
                false,
            );
            let candidate_stereo = targets[target_index].bonds()[target_bond].stereo();
            if params.bond_compare_parameters.match_stereo
                && candidate_stereo.rdkit_code() > BondStereo::Any.rdkit_code()
            {
                stereo = candidate_stereo;
            }
        }
        if need_bond_ring_queries {
            let rings = ring_info.expect("required ring information is present");
            let mut ring_query = make_bond_is_in_ring_query();
            ring_query.set_negation(rings.num_bond_rings(source.id()) == 0);
            query_bond_expand_query(&mut predicate, ring_query, CompositeQueryType::And, true);
        }
        let carrier = Bond::from_spec(
            BondId::new(result_bond),
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Unspecified)
                .with_stereo(stereo),
        );
        bonds.push(QueryBond::from_parts(carrier, predicate));
    }

    QueryGraph::from_parts(
        atoms,
        bonds,
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .map_err(|error| McsCandidateMatchError::ResultQueryGraph {
        message: error.to_string(),
    })
}

#[allow(clippy::too_many_arguments)]
fn mcs_generate_result_smarts_and_query_graph(
    fragment: &McsMoleculeFragment,
    query: &SearchTarget<'_>,
    targets: &[SearchTarget<'_>],
    target_tables: &[McsMatchTables],
    params: &McsParameters,
    final_check: Option<
        &mut dyn FnMut(usize, &[(usize, usize)], &McsParameters) -> Result<bool, McsError>,
    >,
) -> Result<(String, QueryGraph), McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: MaximumCommonSubgraph::generateResultSMARTSAndQueryMol serialization
    // RDKit✔️✔️:   return std::make_pair(MolToSmarts(*mol, true), molSptr);
    // END RDKIT CPP FUNCTION
    // Local complexity review: both paths serialize the completed result query
    // exactly once. The detached writer owns its traversal buffers and returns
    // one SMARTS string while the already built QueryGraph is moved unchanged.
    let graph =
        mcs_build_result_query_graph(fragment, query, targets, target_tables, params, final_check)?;
    let smarts = crate::query_graph_to_smarts(&graph, &crate::SmartsWriteParams::default())
        .map_err(|error| McsCandidateMatchError::ResultSmarts {
            message: error.to_string(),
        })?;
    Ok((smarts, graph))
}

/// Find the source-shaped MCS over borrowed detached inputs.
pub fn find_mcs(
    molecules: &[SearchTarget<'_>],
    params: &McsParameters,
) -> Result<crate::McsResult, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: MaximumCommonSubgraph::find
    // RDKit❗❌: MCSResult MaximumCommonSubgraph::find(const std::vector<ROMOL_SPTR> &src_mols) {
    // RDKit❗❌:   clear();
    // RDKit❗❌:   MCSResult res;
    // RDKit❗❌:
    // RDKit❗❌:   if (src_mols.size() < 2) {
    // RDKit❗❌:     throw std::runtime_error(
    // RDKit❗❌:         "FMCS. Invalid argument. mols.size() must be at least 2");
    // RDKit❗❌:   }
    // RDKit❗❌:   if (Parameters.Threshold > 1.0) {
    // RDKit❗❌:     throw std::runtime_error(
    // RDKit❗❌:         "FMCS. Invalid argument. Parameter Threshold must be 1.0 or "
    // RDKit❗❌:         "less.");
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // minimal required number of matched targets:
    // RDKit❗❌:   // at least one target, max all targets
    // RDKit❗❌:   ThresholdCount = static_cast<unsigned int>(std::min(
    // RDKit❗❌:       static_cast<int>(src_mols.size()) - 1,
    // RDKit❗❌:       std::max(1, static_cast<int>(ceil(static_cast<double>(src_mols.size()) *
    // RDKit❗❌:                                         Parameters.Threshold)) -
    // RDKit❗❌:                       1)));
    // RDKit❗❌:
    // RDKit❗❌:   // AtomCompareParameters.CompleteRingsOnly implies
    // RDKit❗❌:   // BondCompareParameters.CompleteRingsOnly
    // RDKit❗❌:   if (Parameters.AtomCompareParameters.CompleteRingsOnly) {
    // RDKit❗❌:     Parameters.BondCompareParameters.CompleteRingsOnly = true;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // Selecting CompleteRingsOnly option also enables
    // RDKit❗❌:   // --ring-matches-ring-only. ring--ring and chain bonds only match chain
    // RDKit❗❌:   // bonds.
    // RDKit❗❌:   if (Parameters.BondCompareParameters.CompleteRingsOnly) {
    // RDKit❗❌:     Parameters.BondCompareParameters.RingMatchesRingOnly = true;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (Parameters.AtomCompareParameters.CompleteRingsOnly) {
    // RDKit❗❌:     Parameters.AtomCompareParameters.RingMatchesRingOnly = true;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   unsigned int i = 0;
    // RDKit❗❌:   boost::dynamic_bitset<> faked_ring_info(src_mols.size());
    // RDKit❗❌:   for (const auto &src_mol : src_mols) {
    // RDKit❗❌:     Molecules.push_back(src_mol.get());
    // RDKit❗❌:     if (!Molecules.back()->getRingInfo()->isInitialized()) {
    // RDKit❗❌:       Molecules.back()->getRingInfo()->initialize();  // but do not fill out !!!
    // RDKit❗❌:       faked_ring_info.set(i);
    // RDKit❗❌:     }
    // RDKit❗❌:     ++i;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // sort source set of molecules by their 'size' and assume the smallest
    // RDKit❗❌:   // molecule as a query
    // RDKit❗❌:   std::stable_sort(Molecules.begin(), Molecules.end(), molPtr_NumBondLess);
    // RDKit❗❌:   size_t startIdx = 0;
    // RDKit❗❌:   size_t endIdx = Molecules.size() - ThresholdCount;
    // RDKit❗❌:   while (startIdx < endIdx && !Molecules.at(startIdx)->getNumAtoms()) {
    // RDKit❗❌:     ++startIdx;
    // RDKit❗❌:   }
    // RDKit❗❌:   bool areSeedsEmpty = false;
    // RDKit❗❌:   for (size_t i = startIdx; i < endIdx && !areSeedsEmpty && !res.Canceled;
    // RDKit❗❌:        ++i) {
    // RDKit❗❌:     init(startIdx);
    // RDKit❗❌:     if (Targets.empty()) {
    // RDKit❗❌:       break;
    // RDKit❗❌:     }
    // RDKit❗❌:     MCSFinalMatchCheckFunction tff = Parameters.FinalMatchChecker;
    // RDKit❗❌:     // skip final match check for initial seed to allow future growing
    // RDKit❗❌:     Parameters.FinalMatchChecker = nullptr;
    // RDKit❗❌:     makeInitialSeeds();
    // RDKit❗❌:     Parameters.FinalMatchChecker = tff;  // restore final functor
    // RDKit❗❌:
    // RDKit❗❌:     if (Parameters.Verbose) {
    // RDKit❗❌:       std::cout << "Query " << MolToSmiles(*QueryMolecule) << " "
    // RDKit❗❌:                 << QueryMolecule->getNumAtoms() << "("
    // RDKit❗❌:                 << QueryMoleculeMatchedAtoms << ") atoms, "
    // RDKit❗❌:                 << QueryMolecule->getNumBonds() << "("
    // RDKit❗❌:                 << QueryMoleculeMatchedBonds << ") bonds\n";
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     areSeedsEmpty = Seeds.empty();
    // RDKit❗❌:     res.Canceled = !(areSeedsEmpty || growSeeds());
    // RDKit❗❌:     // verify what MCS is equal to one of initial seed for chirality match
    // RDKit❗❌:     if (getMaxNumberBonds() == 0) {
    // RDKit❗❌:       McsIdx = MCS();      // clear
    // RDKit❗❌:       makeInitialSeeds();  // check all possible initial seeds
    // RDKit❗❌:       if (!areSeedsEmpty) {
    // RDKit❗❌:         const Seed &fs = Seeds.front();
    // RDKit❗❌:         if ((1 == getMaxNumberBonds() ||
    // RDKit❗❌:              !(Parameters.BondCompareParameters.CompleteRingsOnly &&
    // RDKit❗❌:                fs.MoleculeFragment.Bonds.size() == 1 &&
    // RDKit❗❌:                queryIsBondInRing(fs.MoleculeFragment.Bonds.front()))) &&
    // RDKit❗❌:             checkIfShouldAcceptMCS(fs.MoleculeFragment, *QueryMolecule, Targets,
    // RDKit❗❌:                                    Parameters)) {
    // RDKit❗❌:           McsIdx.QueryMolecule = QueryMolecule;
    // RDKit❗❌:           McsIdx.Targets = Targets;
    // RDKit❗❌:           McsIdx.Atoms = fs.MoleculeFragment.Atoms;
    // RDKit❗❌:           McsIdx.Bonds = fs.MoleculeFragment.Bonds;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       if (!McsIdx.QueryMolecule && QueryMoleculeSingleMatchedAtom) {
    // RDKit❗❌:         McsIdx.QueryMolecule = QueryMolecule;
    // RDKit❗❌:         McsIdx.Targets = Targets;
    // RDKit❗❌:         McsIdx.Atoms =
    // RDKit❗❌:             std::vector<const Atom *>{QueryMoleculeSingleMatchedAtom};
    // RDKit❗❌:         McsIdx.Bonds = std::vector<const Bond *>();
    // RDKit❗❌:       }
    // RDKit❗❌:     } else if (i + 1 < endIdx) {
    // RDKit❗❌:       Seed seed;
    // RDKit❗❌:       if (createSeedFromMCS(i, seed)) {  // MCS is matched with new query
    // RDKit❗❌:         Seeds.push_back(seed);
    // RDKit❗❌:       }
    // RDKit❗❌:       std::swap(
    // RDKit❗❌:           Molecules.at(startIdx),
    // RDKit❗❌:           Molecules.at(i + 1));  // change query molecule for threshold < 1.
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   res.NumAtoms = getMaxNumberAtoms();
    // RDKit❗❌:   if (!res.NumAtoms && QueryMoleculeSingleMatchedAtom) {
    // RDKit❗❌:     res.NumAtoms = 1;
    // RDKit❗❌:   }
    // RDKit❗❌:   res.NumBonds = getMaxNumberBonds();
    // RDKit❗❌:
    // RDKit❗❌:   if (res.NumBonds > 0 || QueryMoleculeSingleMatchedAtom) {
    // RDKit❗❌:     if (!Parameters.StoreAll) {
    // RDKit❗❌:       auto smartsQueryMolPair = generateResultSMARTSAndQueryMol(McsIdx);
    // RDKit❗❌:       res.SmartsString = std::move(smartsQueryMolPair.first);
    // RDKit❗❌:       res.QueryMol = std::move(smartsQueryMolPair.second);
    // RDKit❗❌:     } else {
    // RDKit❗❌:       std::transform(DegenerateMcsMap.begin(), DegenerateMcsMap.end(),
    // RDKit❗❌:                      std::inserter(res.DegenerateSmartsQueryMolDict,
    // RDKit❗❌:                                    res.DegenerateSmartsQueryMolDict.end()),
    // RDKit❗❌:                      [this](const auto &pair) {
    // RDKit❗❌:                        return generateResultSMARTSAndQueryMol(pair.second);
    // RDKit❗❌:                      });
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌: #ifdef VERBOSE_STATISTICS_ON
    // RDKit❗❌:   if (Parameters.Verbose && res.NumAtoms > 0) {
    // RDKit❗❌:     for (const auto &tag : Targets) {
    // RDKit❗❌:       unsigned int itarget = &tag - &Targets.front();
    // RDKit❗❌:       MatchVectType match;
    // RDKit❗❌:
    // RDKit❗❌:       bool target_matched =
    // RDKit❗❌:           res.QueryMol && SubstructMatch(*tag.Molecule, *res.QueryMol, match)
    // RDKit❗❌:               ? true
    // RDKit❗❌:               : false;
    // RDKit❗❌:       if (!target_matched) {
    // RDKit❗❌:         std::cout << "Target " << itarget + 1
    // RDKit❗❌:                   << (target_matched ? " matched " : " MISMATCHED ")
    // RDKit❗❌:                   << MolToSmiles(*tag.Molecule) << "\n";
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     std::cout << "STATISTICS:\n";
    // RDKit❗❌:     std::cout << "Total Growing Steps  = " << VerboseStatistics.TotalSteps
    // RDKit❗❌:               << ", MCS found on " << VerboseStatistics.MCSFoundStep << " step";
    // RDKit❗❌:     if (VerboseStatistics.MCSFoundTime - To > 0) {
    // RDKit❗❌:       printf(", for %.4lf seconds\n",
    // RDKit❗❌:              double(VerboseStatistics.MCSFoundTime - To) / 1000000.);
    // RDKit❗❌:     } else {
    // RDKit❗❌:       std::cout << ", for less than 1 second\n";
    // RDKit❗❌:     }
    // RDKit❗❌:     std::cout << "Initial   Seeds      = " << VerboseStatistics.InitialSeed
    // RDKit❗❌:               << ",  Mismatched " << VerboseStatistics.MismatchedInitialSeed
    // RDKit❗❌:               << "\n";
    // RDKit❗❌:     std::cout << "Inspected Seeds      = " << VerboseStatistics.Seed << "\n";
    // RDKit❗❌:     std::cout << "Rejected by BestSize = "
    // RDKit❗❌:               << VerboseStatistics.RemainingSizeRejected << "\n";
    // RDKit❗❌:     std::cout << "IndividualBondExcluded   = "
    // RDKit❗❌:               << VerboseStatistics.IndividualBondExcluded << "\n";
    // RDKit❗❌: #ifdef EXCLUDE_WRONG_COMPOSITION
    // RDKit❗❌:     std::cout << "Rejected by WrongComposition = "
    // RDKit❗❌:               << VerboseStatistics.WrongCompositionRejected << " [ "
    // RDKit❗❌:               << VerboseStatistics.WrongCompositionDetected << " Detected ]\n";
    // RDKit❗❌: #endif
    // RDKit❗❌:     std::cout << "MatchCheck Seeds     = " << VerboseStatistics.SeedCheck
    // RDKit❗❌:               << "\n";
    // RDKit❗❌:     std::cout  //<< "\n"
    // RDKit❗❌:         << "     MatchCalls = " << VerboseStatistics.MatchCall << "\n"
    // RDKit❗❌:         << "     MatchFound = " << VerboseStatistics.MatchCallTrue << "\n";
    // RDKit❗❌:     std::cout << " fastMatchCalls = " << VerboseStatistics.FastMatchCall << "\n"
    // RDKit❗❌:               << " fastMatchFound = " << VerboseStatistics.FastMatchCallTrue
    // RDKit❗❌:               << "\n";
    // RDKit❗❌:     std::cout << " slowMatchCalls = "
    // RDKit❗❌:               << VerboseStatistics.MatchCall -
    // RDKit❗❌:                      VerboseStatistics.FastMatchCallTrue
    // RDKit❗❌:               << "\n"
    // RDKit❗❌:               << " slowMatchFound = " << VerboseStatistics.SlowMatchCallTrue
    // RDKit❗❌:               << "\n";
    // RDKit❗❌:
    // RDKit❗❌: #ifdef VERBOSE_STATISTICS_FASTCALLS_ON
    // RDKit❗❌:     std::cout << "AtomFunctorCalls = " << VerboseStatistics.AtomFunctorCalls
    // RDKit❗❌:               << "\n";
    // RDKit❗❌:     std::cout << "BondCompareCalls = " << VerboseStatistics.BondCompareCalls
    // RDKit❗❌:               << "\n";
    // RDKit❗❌: #endif
    // RDKit❗❌:     std::cout << "  DupCacheFound = " << VerboseStatistics.DupCacheFound
    // RDKit❗❌:               << "   " << VerboseStatistics.DupCacheFoundMatch << " matched, "
    // RDKit❗❌:               << VerboseStatistics.DupCacheFound -
    // RDKit❗❌:                      VerboseStatistics.DupCacheFoundMatch
    // RDKit❗❌:               << " mismatched\n";
    // RDKit❗❌: #ifdef FAST_SUBSTRUCT_CACHE
    // RDKit❗❌:     std::cout << "HashCache size  = " << HashCache.keyssize() << " keys\n";
    // RDKit❗❌:     std::cout << "HashCache size  = " << HashCache.fullsize() << " entries\n";
    // RDKit❗❌:     std::cout << "FindHashInCache = " << VerboseStatistics.FindHashInCache
    // RDKit❗❌:               << "\n";
    // RDKit❗❌:     std::cout << "HashFoundInCache= " << VerboseStatistics.HashKeyFoundInCache
    // RDKit❗❌:               << "\n";
    // RDKit❗❌:     std::cout << "ExactMatchCalls = " << VerboseStatistics.ExactMatchCall
    // RDKit❗❌:               << "\n"
    // RDKit❗❌:               << "ExactMatchFound = " << VerboseStatistics.ExactMatchCallTrue
    // RDKit❗❌:               << "\n";
    // RDKit❗❌: #endif
    // RDKit❗❌:   }
    // RDKit❗❌: #endif
    // RDKit❗❌:
    // RDKit❗❌:   auto pos = faked_ring_info.find_first();
    // RDKit❗❌:   while (pos != boost::dynamic_bitset<>::npos) {
    // RDKit❗❌:     src_mols[pos]->getRingInfo()->reset();
    // RDKit❗❌:     pos = faked_ring_info.find_next(pos);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   clear();
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION
    // Local complexity review: stable index sorting and result-context indices
    // avoid molecule clones. Every target table is rebuilt for each source
    // init, and full matching retains its documented VF2-pruning gap.
    let order = prepare_mcs_input_order(molecules, params.threshold)?;
    let started = std::time::Instant::now();
    let mut effective_params = params.clone();
    if effective_params.atom_compare_parameters.match_chiral_tag {
        effective_params.bond_compare_parameters.match_stereo = true;
    }
    if effective_params.atom_compare_parameters.complete_rings_only {
        effective_params.bond_compare_parameters.complete_rings_only = true;
    }
    if effective_params.bond_compare_parameters.complete_rings_only {
        effective_params
            .bond_compare_parameters
            .ring_matches_ring_only = true;
    }
    if effective_params.atom_compare_parameters.complete_rings_only {
        effective_params
            .atom_compare_parameters
            .ring_matches_ring_only = true;
    }

    // The source initializes an absent RingInfo without finding rings. Keep
    // that empty state local and preserve every other borrowed target field.
    let empty_rings = molecules
        .iter()
        .map(|molecule| {
            molecule
                .ring_info()
                .filter(|info| info.is_initialized())
                .is_none()
                .then(|| {
                    RingInfo::new(
                        RingFindType::OtherOrUnknown,
                        molecule.num_atoms(),
                        molecule.num_bonds(),
                    )
                })
        })
        .collect::<Vec<_>>();
    let inputs = molecules
        .iter()
        .enumerate()
        .map(|(index, molecule)| {
            empty_rings[index]
                .as_ref()
                .map_or(*molecule, |rings| molecule.with_ring_info(rings))
        })
        .collect::<Vec<_>>();
    let mut sorted = order.molecule_indices;
    let mut state = McsSearchState::default();
    let mut canceled = false;
    let mut are_seeds_empty = false;
    let mut query_single_matched_atom = None;

    for index in order.start_index..order.end_index {
        if are_seeds_empty || canceled {
            break;
        }
        let query_input = sorted[order.start_index];
        let query = &inputs[query_input];
        // Source init(startIdx) always uses position startIdx as query, and
        // Targets always begin at sorted position one, even after swaps.
        let target_inputs = sorted[1..].to_vec();
        let targets = target_inputs
            .iter()
            .map(|&target_input| inputs[target_input])
            .collect::<Vec<_>>();
        if targets.is_empty() {
            break;
        }
        let target_tables = targets
            .iter()
            .map(|target| build_query_target_match_tables(&effective_params, query, target))
            .collect::<Result<Vec<_>, _>>()?;
        let context = McsResultContext {
            query_input,
            target_inputs,
        };
        let mut initial = mcs_make_initial_seeds(
            query,
            &targets,
            &target_tables,
            order.threshold_count,
            &effective_params,
            false,
            None,
            None,
        )?;
        query_single_matched_atom = initial.query_single_matched_atom;
        are_seeds_empty = initial.queue.seeds.is_empty();
        if !are_seeds_empty {
            let mut now = || started.elapsed().as_nanos() as u64;
            canceled = mcs_grow_seeds(
                &mut initial.queue,
                &mut state,
                &context,
                query,
                &targets,
                &target_tables,
                order.threshold_count,
                initial.query_matched_bonds,
                &effective_params,
                0,
                &mut now,
                None,
                None,
                None,
            )?
            .canceled;
        }
        if state.best.bonds.is_empty() {
            state.clear_best_for_zero_bond_fallback();
            let second = mcs_make_initial_seeds(
                query,
                &targets,
                &target_tables,
                order.threshold_count,
                &effective_params,
                true,
                None,
                None,
            )?;
            query_single_matched_atom = second.query_single_matched_atom;
            if !are_seeds_empty {
                if let Some(first) = second.queue.seeds.first() {
                    let reject_single_ring =
                        effective_params.bond_compare_parameters.complete_rings_only
                            && first.molecule_fragment.bonds.len() == 1
                            && query.ring_info().is_some_and(|rings| {
                                rings.num_bond_rings(
                                    query.bonds()[first.molecule_fragment.bonds[0]].id(),
                                ) > 0
                            });
                    if !reject_single_ring {
                        state.best = first.molecule_fragment.clone();
                        state.best_context = Some(context.clone());
                    }
                }
            }
            if state.best_context.is_none() {
                if let Some(atom) = query_single_matched_atom {
                    state.best.atoms = vec![atom];
                    state.best.bonds.clear();
                    state.best_context = Some(context.clone());
                }
            }
        } else if index + 1 < order.end_index {
            let retained = state
                .best_context
                .as_ref()
                .ok_or(McsCandidateMatchError::ResultContextMissing)?;
            let &new_query_input = retained.target_inputs.get(index).ok_or(
                McsCandidateMatchError::TargetTableCount {
                    targets: retained.target_inputs.len(),
                    tables: index + 1,
                },
            )?;
            let old_query = &inputs[retained.query_input];
            let new_query = &inputs[new_query_input];
            let new_query_table =
                build_query_target_match_tables(&effective_params, old_query, new_query)?;
            if let Some(seed) = mcs_create_seed_from_mcs(
                &state.best,
                old_query,
                new_query,
                &new_query_table,
                &effective_params,
                None,
            )? {
                initial.queue.seeds.push(seed);
            }
            sorted.swap(order.start_index, index + 1);
        }
    }

    let atom_count = state.best.atoms.len().max(usize::from(
        state.best.atoms.is_empty() && query_single_matched_atom.is_some(),
    ));
    let bond_count = state.best.bonds.len();
    let mut result = crate::McsResult {
        query: None,
        atom_count,
        bond_count,
        completed: !canceled,
        smarts: String::new(),
        degenerate: BTreeMap::new(),
    };
    if bond_count > 0 || query_single_matched_atom.is_some() {
        if !effective_params.store_all {
            let context = state
                .best_context
                .as_ref()
                .ok_or(McsCandidateMatchError::ResultContextMissing)?;
            let query = &inputs[context.query_input];
            let targets = context
                .target_inputs
                .iter()
                .map(|&input| inputs[input])
                .collect::<Vec<_>>();
            let tables = targets
                .iter()
                .map(|target| build_query_target_match_tables(&effective_params, query, target))
                .collect::<Result<Vec<_>, _>>()?;
            let (smarts, graph) = mcs_generate_result_smarts_and_query_graph(
                &state.best,
                query,
                &targets,
                &tables,
                &effective_params,
                None,
            )?;
            result.smarts = smarts;
            result.query = Some(graph);
        } else {
            for retained in state.degenerate.values() {
                let query = &inputs[retained.context.query_input];
                let targets = retained
                    .context
                    .target_inputs
                    .iter()
                    .map(|&input| inputs[input])
                    .collect::<Vec<_>>();
                let tables = targets
                    .iter()
                    .map(|target| build_query_target_match_tables(&effective_params, query, target))
                    .collect::<Result<Vec<_>, _>>()?;
                let (smarts, graph) = mcs_generate_result_smarts_and_query_graph(
                    &retained.fragment,
                    query,
                    &targets,
                    &tables,
                    &effective_params,
                    None,
                )?;
                result.degenerate.entry(smarts).or_insert(graph);
            }
        }
    }
    Ok(result)
}

fn mcs_create_seed_from_mcs(
    fragment: &McsMoleculeFragment,
    query: &SearchTarget<'_>,
    new_query: &SearchTarget<'_>,
    new_query_tables: &McsMatchTables,
    params: &McsParameters,
    mut final_check: Option<
        &mut dyn FnMut(&[(usize, usize)], &McsParameters) -> Result<bool, McsError>,
    >,
) -> Result<Option<McsSeed>, McsCandidateMatchError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: MaximumCommonSubgraph::createSeedFromMCS
    // RDKit✔️❌: bool MaximumCommonSubgraph::createSeedFromMCS(size_t newQueryTarget,
    // RDKit✔️❌:                                               Seed &newSeed) {
    // RDKit✔️❌:   Seed mcs;
    // RDKit✔️❌:   mcs.setStoreAllDegenerateMCS(Parameters.StoreAll);
    // RDKit✔️❌:   mcs.ExcludedBonds.resize(McsIdx.QueryMolecule->getNumBonds(), false);
    // RDKit✔️❌:   std::vector<unsigned int> mcsAtomIdxMap(McsIdx.QueryMolecule->getNumAtoms());
    // RDKit✔️❌:
    // RDKit✔️❌:   for (const auto &atom : McsIdx.Atoms) {
    // RDKit✔️❌:     mcsAtomIdxMap[atom->getIdx()] = mcs.addAtom(atom);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (const auto &bond : McsIdx.Bonds) {
    // RDKit✔️❌:     mcs.addBond(bond);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   const Target &newQuery = McsIdx.Targets.at(newQueryTarget);
    // RDKit✔️❌:
    // RDKit✔️❌:   match_V_t match;
    // RDKit✔️❌:   bool target_matched = SubstructMatchCustomTable(
    // RDKit✔️❌:       newQuery.Topology, *newQuery.Molecule, mcs.Topology,
    // RDKit✔️❌:       *McsIdx.QueryMolecule, newQuery.AtomMatchTable, newQuery.BondMatchTable,
    // RDKit✔️❌:       &Parameters, &match);
    // RDKit✔️❌:   if (!target_matched) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   AtomMatchSet atomMatchResult(mcs.getNumAtoms());
    // RDKit✔️❌:
    // RDKit✔️❌:   newSeed.ExcludedBonds.resize(newQuery.Molecule->getNumBonds(), false);
    // RDKit✔️❌:
    // RDKit✔️❌:   for (const auto &m : match) {
    // RDKit✔️❌:     unsigned int ai = m.first;  // SeedAtomIdx in mcs seed
    // RDKit✔️❌:     atomMatchResult[ai].QueryAtomIdx = mcs.Topology[m.first];
    // RDKit✔️❌:     atomMatchResult[ai].TargetAtomIdx = newQuery.Topology[m.second];
    // RDKit✔️❌:     const auto ta =
    // RDKit✔️❌:         newQuery.Molecule->getAtomWithIdx(newQuery.Topology[m.second]);
    // RDKit✔️❌:     newSeed.addAtom(ta);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (const auto &bond : McsIdx.Bonds) {
    // RDKit✔️❌:     unsigned int i = mcsAtomIdxMap.at(bond->getBeginAtomIdx());
    // RDKit✔️❌:     unsigned int j = mcsAtomIdxMap.at(bond->getEndAtomIdx());
    // RDKit✔️❌:     unsigned int ti = atomMatchResult.at(i).TargetAtomIdx;
    // RDKit✔️❌:     unsigned int tj = atomMatchResult.at(j).TargetAtomIdx;
    // RDKit✔️❌:     const auto tb = newQuery.Molecule->getBondBetweenAtoms(ti, tj);
    // RDKit✔️❌:     CHECK_INVARIANT(tb, "tb most not be NULL");
    // RDKit✔️❌:     newSeed.addBond(tb);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   newSeed.computeRemainingSize(*newQuery.Molecule);
    // RDKit✔️❌:   return true;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION
    // Local complexity review: seed construction, source-to-seed indexing,
    // target reconstruction and remaining-size traversal match the source.
    // Full rematching inherits the documented backtracking matcher gap versus
    // Boost VF2, so the performance marker remains ❌.
    let mut mcs = McsSeed {
        store_all_degenerate_mcs: params.store_all,
        excluded_bonds: vec![false; query.num_bonds()],
        ..McsSeed::default()
    };
    let mut mcs_atom_index_map = vec![usize::MAX; query.num_atoms()];
    for &source_atom in &fragment.atoms {
        if source_atom >= query.num_atoms() {
            return Err(McsError::AtomOutOfRange {
                side: "MCS seed reconstruction query",
                atom: source_atom,
            }
            .into());
        }
        mcs_atom_index_map[source_atom] = mcs.add_atom(source_atom);
    }
    for &source_bond in &fragment.bonds {
        mcs.add_bond(query, source_bond)?;
    }

    let mut accept_mapping = |mapping: &[(usize, usize)]| {
        let has_user_hook = final_check.is_some();
        let mut user_hook = || {
            final_check
                .as_deref_mut()
                .expect("present reconstruction final hook")(mapping, params)
            .map_err(McsCandidateMatchError::from)
        };
        let user = if has_user_hook {
            Some(&mut user_hook as &mut dyn FnMut() -> Result<bool, McsCandidateMatchError>)
        } else {
            None
        };
        mcs_final_mapping_accept(&mcs.topology, query, new_query, params, mapping, user)
    };
    let mapping = mcs_find_full_mapping(&mcs, new_query, new_query_tables, &mut accept_mapping)?;
    let Some(mapping) = mapping else {
        return Ok(None);
    };

    let mut target_atoms = vec![usize::MAX; mcs.topology.source_atoms.len()];
    for (seed_atom, target_atom) in mapping {
        let mapped = target_atoms
            .get_mut(seed_atom)
            .ok_or(McsError::AtomOutOfRange {
                side: "MCS seed reconstruction mapping",
                atom: seed_atom,
            })?;
        if target_atom >= new_query.num_atoms() {
            return Err(McsError::AtomOutOfRange {
                side: "MCS seed reconstruction target",
                atom: target_atom,
            }
            .into());
        }
        *mapped = target_atom;
    }

    let mut new_seed = McsSeed {
        store_all_degenerate_mcs: params.store_all,
        excluded_bonds: vec![false; new_query.num_bonds()],
        ..McsSeed::default()
    };
    for &target_atom in &target_atoms {
        if target_atom == usize::MAX {
            return Err(McsError::TargetAtomMappingMissing {
                atom: new_seed.molecule_fragment.atoms.len(),
            }
            .into());
        }
        new_seed.add_atom(target_atom);
    }

    for &source_bond in &fragment.bonds {
        let bond = query
            .bonds()
            .get(source_bond)
            .ok_or(McsError::BondOutOfRange {
                side: "MCS seed reconstruction query",
                bond: source_bond,
            })?;
        let begin_seed_atom = mcs_atom_index_map[bond.begin().index()];
        let end_seed_atom = mcs_atom_index_map[bond.end().index()];
        let target_begin = target_atoms[begin_seed_atom];
        let target_end = target_atoms[end_seed_atom];
        let target_bond = mcs_target_bond_between(new_query, target_begin, target_end).ok_or(
            McsCandidateMatchError::SeedReconstructionBondMissing {
                begin: target_begin,
                end: target_end,
            },
        )?;
        new_seed.add_bond(new_query, target_bond)?;
    }
    new_seed.compute_remaining_size(new_query)?;
    Ok(Some(new_seed))
}

fn build_initial_match_tables(
    params: &McsParameters,
    query: &SearchTarget<'_>,
) -> Result<McsMatchTables, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: MaximumCommonSubgraph::init match tables
    // RDKit✔️❌:   // fill out match tables
    // RDKit✔️❌:   nq = QueryMolecule->getNumAtoms();
    // RDKit✔️❌:   QueryAtomMatchTable.resize(nq, nq);
    // RDKit✔️❌:   for (size_t aj = 0; aj < nq; aj++) {
    // RDKit✔️❌:     for (size_t ai = 0; ai < nq; ai++) {
    // RDKit✔️❌:       QueryAtomMatchTable.set(
    // RDKit✔️❌:           ai, aj,
    // RDKit✔️❌:           Parameters.AtomTyper(Parameters.AtomCompareParameters, *QueryMolecule,
    // RDKit✔️❌:                                ai, *QueryMolecule, aj,
    // RDKit✔️❌:                                Parameters.CompareFunctionsUserData));
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   nq = QueryMolecule->getNumBonds();
    // RDKit✔️❌:   QueryBondMatchTable.resize(nq, nq);
    // RDKit✔️❌:   for (size_t aj = 0; aj < nq; aj++) {
    // RDKit✔️❌:     for (size_t ai = 0; ai < nq; ai++) {
    // RDKit✔️❌:       QueryBondMatchTable.set(
    // RDKit✔️❌:           ai, aj,
    // RDKit✔️❌:           Parameters.BondTyper(Parameters.BondCompareParameters, *QueryMolecule,
    // RDKit✔️❌:                                ai, *QueryMolecule, aj,
    // RDKit✔️❌:                                Parameters.CompareFunctionsUserData));
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION
    // Source query-self tables are the rectangular constructor with both
    // operands borrowing the same query. No table or comparator is duplicated.
    build_query_target_match_tables(params, query, query)
}

fn build_query_target_match_tables(
    params: &McsParameters,
    query: &SearchTarget<'_>,
    target: &SearchTarget<'_>,
) -> Result<McsMatchTables, McsError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.cpp :: MaximumCommonSubgraph::init target match tables
    // RDKit✔️❌:     size_t nq = QueryMolecule->getNumAtoms();
    // RDKit✔️❌:     size_t nt = Targets.at(i).Molecule->getNumAtoms();
    // RDKit✔️❌:     Targets[i].AtomMatchTable.resize(nq, nt);
    // RDKit✔️❌:     for (size_t aj = 0; aj < nt; aj++) {
    // RDKit✔️❌:       for (size_t ai = 0; ai < nq; ai++) {
    // RDKit✔️❌:         Targets[i].AtomMatchTable.set(
    // RDKit✔️❌:             ai, aj,
    // RDKit✔️❌:             Parameters.AtomTyper(Parameters.AtomCompareParameters,
    // RDKit✔️❌:                                  *QueryMolecule, ai, *Targets.at(i).Molecule,
    // RDKit✔️❌:                                  aj, Parameters.CompareFunctionsUserData));
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     nq = QueryMolecule->getNumBonds();
    // RDKit✔️❌:     nt = Targets.at(i).Molecule->getNumBonds();
    // RDKit✔️❌:     Targets[i].BondMatchTable.resize(nq, nt);
    // RDKit✔️❌:     for (size_t aj = 0; aj < nt; aj++) {
    // RDKit✔️❌:       for (size_t ai = 0; ai < nq; ai++) {
    // RDKit✔️❌:         Targets[i].BondMatchTable.set(
    // RDKit✔️❌:             ai, aj,
    // RDKit✔️❌:             Parameters.BondTyper(Parameters.BondCompareParameters,
    // RDKit✔️❌:                                  *QueryMolecule, ai, *Targets.at(i).Molecule,
    // RDKit✔️❌:                                  aj, Parameters.CompareFunctionsUserData));
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // END RDKIT CPP FUNCTION
    // Local complexity review: each table is one row-major Vec<bool>
    // allocation. Construction is O(query atoms * target atoms + query bonds
    // * target bonds), plus the costs of the existing comparators. The source
    // TArray2D<bool> uses vector<bool> and packs bits; Rust Vec<bool> stores
    // one byte per cell, so the storage cost is materially higher.
    let mut atoms = McsMatchTable::new(query.num_atoms(), target.num_atoms());
    for right_atom in 0..target.num_atoms() {
        for left_atom in 0..query.num_atoms() {
            let matches = match params.atom_comparator {
                AtomComparator::AtomCompareAny => mcs_atom_compare_any(
                    &params.atom_compare_parameters,
                    query,
                    left_atom,
                    target,
                    right_atom,
                )?,
                AtomComparator::AtomCompareElements => mcs_atom_compare_elements(
                    &params.atom_compare_parameters,
                    query,
                    left_atom,
                    target,
                    right_atom,
                )?,
                AtomComparator::AtomCompareIsotopes => mcs_atom_compare_isotopes(
                    &params.atom_compare_parameters,
                    query,
                    left_atom,
                    target,
                    right_atom,
                )?,
                AtomComparator::AtomCompareAnyHeavyAtom => mcs_atom_compare_any_heavy(
                    &params.atom_compare_parameters,
                    query,
                    left_atom,
                    target,
                    right_atom,
                )?,
            };
            atoms.set(left_atom, right_atom, matches);
        }
    }

    let mut bonds = McsMatchTable::new(query.num_bonds(), target.num_bonds());
    for right_bond in 0..target.num_bonds() {
        for left_bond in 0..query.num_bonds() {
            let matches = match params.bond_comparator {
                BondComparator::BondCompareAny => mcs_bond_compare_any(
                    &params.bond_compare_parameters,
                    query,
                    left_bond,
                    target,
                    right_bond,
                )?,
                BondComparator::BondCompareOrder => mcs_bond_compare_order(
                    &params.bond_compare_parameters,
                    query,
                    left_bond,
                    target,
                    right_bond,
                )?,
                BondComparator::BondCompareOrderExact => mcs_bond_compare_order_exact(
                    &params.bond_compare_parameters,
                    query,
                    left_bond,
                    target,
                    right_bond,
                )?,
            };
            bonds.set(left_bond, right_bond, matches);
        }
    }
    Ok(McsMatchTables { atoms, bonds })
}

impl Default for McsParameters {
    fn default() -> Self {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.h :: MCSParameters member initializers
        // RDKit✔️✔️:   bool StoreAll = false;
        // RDKit✔️✔️:   bool MaximizeBonds = true;
        // RDKit✔️✔️:   double Threshold = 1.0;    // match all molecules
        // RDKit✔️✔️:   unsigned int Timeout = 0;  // in seconds
        // RDKit✔️✔️:   bool Verbose = false;
        // RDKit✔️✔️:   MCSAtomCompareFunction AtomTyper = MCSAtomCompareElements;
        // RDKit✔️✔️:   MCSBondCompareFunction BondTyper = MCSBondCompareOrder;
        // RDKit✔️✔️:   std::string InitialSeed = "";  // user defined or empty string (default)
        // END RDKIT CPP FUNCTION
        // Local complexity review: the source and Rust defaults initialize a
        // fixed set of scalars plus one empty string, all in O(1).
        Self {
            store_all: false,
            maximize_bonds: true,
            threshold: 1.0,
            timeout: 0,
            verbose: false,
            atom_compare_parameters: McsAtomCompareParameters::default(),
            bond_compare_parameters: McsBondCompareParameters::default(),
            atom_comparator: AtomComparator::default(),
            bond_comparator: BondComparator::default(),
            initial_seed: String::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum McsParametersJsonError {
    #[error("invalid MCS parameter JSON: {0}")]
    InvalidJson(String),
    #[error("invalid MCS parameter field `{field}`")]
    InvalidField { field: &'static str },
}

impl From<serde_json::Error> for McsParametersJsonError {
    fn from(error: serde_json::Error) -> Self {
        Self::InvalidJson(error.to_string())
    }
}

fn json_bool(
    object: &serde_json::Map<String, serde_json::Value>,
    field: &'static str,
) -> Result<Option<bool>, McsParametersJsonError> {
    let Some(value) = object.get(field) else {
        return Ok(None);
    };
    match value {
        serde_json::Value::Bool(value) => Ok(Some(*value)),
        serde_json::Value::Number(value) if value.as_u64() == Some(1) => Ok(Some(true)),
        serde_json::Value::Number(value) if value.as_u64() == Some(0) => Ok(Some(false)),
        serde_json::Value::String(value) if value == "true" || value == "1" => Ok(Some(true)),
        serde_json::Value::String(value) if value == "false" || value == "0" => Ok(Some(false)),
        _ => Err(McsParametersJsonError::InvalidField { field }),
    }
}

fn json_f64(
    object: &serde_json::Map<String, serde_json::Value>,
    field: &'static str,
) -> Result<Option<f64>, McsParametersJsonError> {
    let Some(value) = object.get(field) else {
        return Ok(None);
    };
    let parsed = match value {
        serde_json::Value::Number(value) => value.as_f64(),
        serde_json::Value::String(value) => value.parse().ok(),
        _ => None,
    };
    parsed
        .map(Some)
        .ok_or(McsParametersJsonError::InvalidField { field })
}

fn json_u32(
    object: &serde_json::Map<String, serde_json::Value>,
    field: &'static str,
) -> Result<Option<u32>, McsParametersJsonError> {
    let Some(value) = object.get(field) else {
        return Ok(None);
    };
    let parsed = match value {
        serde_json::Value::Number(value) => value.as_u64().and_then(|value| value.try_into().ok()),
        serde_json::Value::String(value) => value.parse().ok(),
        _ => None,
    };
    parsed
        .map(Some)
        .ok_or(McsParametersJsonError::InvalidField { field })
}

fn json_string<'a>(
    object: &'a serde_json::Map<String, serde_json::Value>,
    field: &'static str,
) -> Result<Option<&'a str>, McsParametersJsonError> {
    let Some(value) = object.get(field) else {
        return Ok(None);
    };
    value
        .as_str()
        .map(Some)
        .ok_or(McsParametersJsonError::InvalidField { field })
}

/// Apply the source JSON parameter projection to an existing parameter value.
pub fn update_mcs_parameters_from_json(
    params: &mut McsParameters,
    json: &str,
) -> Result<(), McsParametersJsonError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FMCS/FMCS.cpp :: parseMCSParametersJSON
    // RDKit✔️✔️:   if (!params || !json || !strlen(json)) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   p.MaximizeBonds = pt.get<bool>("MaximizeBonds", p.MaximizeBonds);
    // RDKit✔️✔️:   p.Threshold = pt.get<double>("Threshold", p.Threshold);
    // RDKit✔️✔️:   p.Timeout = pt.get<unsigned int>("Timeout", p.Timeout);
    // RDKit✔️✔️:   p.AtomCompareParameters.MatchValences =
    // RDKit✔️✔️:       pt.get<bool>("MatchValences", p.AtomCompareParameters.MatchValences);
    // RDKit✔️✔️:   p.AtomCompareParameters.MatchChiralTag =
    // RDKit✔️✔️:       pt.get<bool>("MatchChiralTag", p.AtomCompareParameters.MatchChiralTag);
    // RDKit✔️✔️:   p.AtomCompareParameters.MatchFormalCharge = pt.get<bool>(
    // RDKit✔️✔️:       "MatchFormalCharge", p.AtomCompareParameters.MatchFormalCharge);
    // RDKit✔️✔️:   p.AtomCompareParameters.RingMatchesRingOnly = pt.get<bool>(
    // RDKit✔️✔️:       "RingMatchesRingOnly", p.AtomCompareParameters.RingMatchesRingOnly);
    // RDKit✔️✔️:   p.AtomCompareParameters.MaxDistance =
    // RDKit✔️✔️:       pt.get<double>("MaxDistance", p.AtomCompareParameters.MaxDistance);
    // RDKit✔️✔️:   p.BondCompareParameters.RingMatchesRingOnly = pt.get<bool>(
    // RDKit✔️✔️:       "RingMatchesRingOnly", p.BondCompareParameters.RingMatchesRingOnly);
    // RDKit✔️✔️:   p.AtomCompareParameters.RingMatchesRingOnly = pt.get<bool>(
    // RDKit✔️✔️:       "AtomRingMatchesRingOnly", p.AtomCompareParameters.RingMatchesRingOnly);
    // RDKit✔️✔️:   p.BondCompareParameters.RingMatchesRingOnly = pt.get<bool>(
    // RDKit✔️✔️:       "BondRingMatchesRingOnly", p.BondCompareParameters.RingMatchesRingOnly);
    // RDKit✔️✔️:   p.BondCompareParameters.CompleteRingsOnly = pt.get<bool>(
    // RDKit✔️✔️:       "CompleteRingsOnly", p.BondCompareParameters.CompleteRingsOnly);
    // RDKit✔️✔️:   p.AtomCompareParameters.CompleteRingsOnly = pt.get<bool>(
    // RDKit✔️✔️:       "AtomCompleteRingsOnly", p.AtomCompareParameters.CompleteRingsOnly);
    // RDKit✔️✔️:   p.BondCompareParameters.CompleteRingsOnly = pt.get<bool>(
    // RDKit✔️✔️:       "BondCompleteRingsOnly", p.BondCompareParameters.CompleteRingsOnly);
    // RDKit✔️✔️:   p.BondCompareParameters.MatchFusedRings =
    // RDKit✔️✔️:       pt.get<bool>("MatchFusedRings", p.BondCompareParameters.MatchFusedRings);
    // RDKit✔️✔️:   p.BondCompareParameters.MatchFusedRingsStrict = pt.get<bool>(
    // RDKit✔️✔️:       "MatchFusedRingsStrict", p.BondCompareParameters.MatchFusedRingsStrict);
    // RDKit✔️✔️:   p.BondCompareParameters.MatchStereo =
    // RDKit✔️✔️:       pt.get<bool>("MatchStereo", p.BondCompareParameters.MatchStereo);
    // RDKit✔️✔️:   p.StoreAll = pt.get<bool>("StoreAll", p.StoreAll);
    // RDKit✔️✔️:   p.setMCSAtomTyperFromConstChar(
    // RDKit✔️✔️:       pt.get<std::string>("AtomCompare", "def").c_str());
    // RDKit✔️✔️:   p.setMCSBondTyperFromConstChar(
    // RDKit✔️✔️:       pt.get<std::string>("BondCompare", "def").c_str());
    // RDKit✔️✔️:   p.InitialSeed = pt.get<std::string>("InitialSeed", "");
    // END RDKIT CPP FUNCTION
    // Local complexity review: parsing is O(input length), then performs a
    // fixed number of expected O(1) object lookups. Assignments occur in the
    // source order, so an invalid later field retains preceding updates.
    if json.is_empty() {
        return Ok(());
    }
    let value: serde_json::Value = serde_json::from_str(json)?;
    let object = value
        .as_object()
        .ok_or(McsParametersJsonError::InvalidField { field: "root" })?;

    if let Some(value) = json_bool(object, "MaximizeBonds")? {
        params.maximize_bonds = value;
    }
    if let Some(value) = json_f64(object, "Threshold")? {
        params.threshold = value;
    }
    if let Some(value) = json_u32(object, "Timeout")? {
        params.timeout = value;
    }
    if let Some(value) = json_bool(object, "MatchValences")? {
        params.atom_compare_parameters.match_valences = value;
    }
    if let Some(value) = json_bool(object, "MatchChiralTag")? {
        params.atom_compare_parameters.match_chiral_tag = value;
    }
    if let Some(value) = json_bool(object, "MatchFormalCharge")? {
        params.atom_compare_parameters.match_formal_charge = value;
    }
    if let Some(value) = json_bool(object, "RingMatchesRingOnly")? {
        params.atom_compare_parameters.ring_matches_ring_only = value;
    }
    if let Some(value) = json_f64(object, "MaxDistance")? {
        params.atom_compare_parameters.max_distance = value;
    }
    if let Some(value) = json_bool(object, "RingMatchesRingOnly")? {
        params.bond_compare_parameters.ring_matches_ring_only = value;
    }
    if let Some(value) = json_bool(object, "AtomRingMatchesRingOnly")? {
        params.atom_compare_parameters.ring_matches_ring_only = value;
    }
    if let Some(value) = json_bool(object, "BondRingMatchesRingOnly")? {
        params.bond_compare_parameters.ring_matches_ring_only = value;
    }
    if let Some(value) = json_bool(object, "CompleteRingsOnly")? {
        params.bond_compare_parameters.complete_rings_only = value;
    }
    if let Some(value) = json_bool(object, "AtomCompleteRingsOnly")? {
        params.atom_compare_parameters.complete_rings_only = value;
    }
    if let Some(value) = json_bool(object, "BondCompleteRingsOnly")? {
        params.bond_compare_parameters.complete_rings_only = value;
    }
    if let Some(value) = json_bool(object, "MatchFusedRings")? {
        params.bond_compare_parameters.match_fused_rings = value;
    }
    if let Some(value) = json_bool(object, "MatchFusedRingsStrict")? {
        params.bond_compare_parameters.match_fused_rings_strict = value;
    }
    if let Some(value) = json_bool(object, "MatchStereo")? {
        params.bond_compare_parameters.match_stereo = value;
    }
    if let Some(value) = json_bool(object, "StoreAll")? {
        params.store_all = value;
    }
    if let Some(value) = json_string(object, "AtomCompare")? {
        params.atom_comparator = match value {
            "Any" => AtomComparator::AtomCompareAny,
            "Elements" => AtomComparator::AtomCompareElements,
            "Isotopes" => AtomComparator::AtomCompareIsotopes,
            "AnyHeavy" => AtomComparator::AtomCompareAnyHeavyAtom,
            _ => params.atom_comparator,
        };
    }
    if let Some(value) = json_string(object, "BondCompare")? {
        params.bond_comparator = match value {
            "Any" => BondComparator::BondCompareAny,
            "Order" => BondComparator::BondCompareOrder,
            "OrderExact" => BondComparator::BondCompareOrderExact,
            _ => params.bond_comparator,
        };
    }
    params.initial_seed = json_string(object, "InitialSeed")?
        .unwrap_or_default()
        .to_owned();
    Ok(())
}

#[cfg(test)]
mod q134_state_ {
    #[test]
    fn best_and_provenance_survive_later_calls() {
        super::tests::q134_state_best_and_provenance_survive_later_worse_empty_and_canceled_calls();
    }

    #[test]
    fn store_all_keeps_tie_context_and_source_key_overwrite() {
        super::tests::q134_state_store_all_keeps_tie_context_and_source_key_overwrite();
    }

    #[test]
    fn objective_and_zero_bond_reset_follow_source() {
        super::tests::q134_state_objective_and_zero_bond_reset_follow_source();
    }
}

#[cfg(test)]
mod final_ring_pass1_ {
    #[test]
    fn zero_one_and_multiple_memberships_use_original_bond_ids() {
        super::tests::final_ring_pass1_zero_one_and_multiple_memberships_use_original_bond_ids();
    }

    #[test]
    fn partial_and_nonidentity_mapping_stays_on_mapped_edge() {
        super::tests::final_ring_pass1_partial_and_nonidentity_mapping_stays_on_mapped_edge();
    }

    #[test]
    fn reports_missing_state_and_invalid_mapped_edges() {
        super::tests::final_ring_pass1_reports_missing_state_and_invalid_mapped_edges();
    }
}

#[cfg(test)]
mod final_ring_refine_ {
    #[test]
    fn pass2_removes_incomplete_fusion_and_preserves_snapshot() {
        super::tests::final_ring_refine_pass2_removes_incomplete_fusion_and_preserves_snapshot();
    }

    #[test]
    fn pass2_reclassifies_one_survivor_and_retains_two() {
        super::tests::final_ring_refine_pass2_reclassifies_one_survivor_and_retains_two();
    }

    #[test]
    fn honor_check_distinguishes_none_all_partial_and_missing_fused() {
        super::tests::final_ring_refine_honor_check_distinguishes_none_all_partial_and_missing_fused();
    }
}

#[cfg(test)]
mod final_ring_fusion_ {
    #[test]
    fn neither_one_and_both_molecule_checks_follow_mode() {
        super::tests::final_ring_fusion_neither_one_and_both_molecule_checks_follow_mode();
    }

    #[test]
    fn target_smaller_shortcut_precedes_required_ring_state() {
        super::tests::final_ring_fusion_target_smaller_shortcut_precedes_required_ring_state();
    }

    #[test]
    fn distinct_query_and_target_vertex_identities_are_preserved() {
        super::tests::final_ring_fusion_distinct_query_and_target_vertex_identities_are_preserved();
    }
}

#[cfg(test)]
mod final_chiral_atoms_ {
    #[test]
    fn query_skip_and_target_tag_and_degree_branches() {
        super::tests::final_chiral_atoms_query_skip_and_target_tag_and_degree_branches();
    }

    #[test]
    fn mapped_neighbor_permutations_use_tag_parity() {
        super::tests::final_chiral_atoms_mapped_neighbor_permutations_use_tag_parity();
    }

    #[test]
    fn repeated_missing_ligands_and_nonidentity_vertices() {
        super::tests::final_chiral_atoms_repeated_missing_ligands_and_nonidentity_vertices();
    }
}

#[cfg(test)]
mod final_chiral_bonds_ {
    #[test]
    fn source_skip_branches_and_target_orientation() {
        super::tests::final_chiral_bonds_source_skip_branches_and_target_orientation();
    }

    #[test]
    fn stereo_label_and_neighbor_match_counts() {
        super::tests::final_chiral_bonds_stereo_label_and_neighbor_match_counts();
    }

    #[test]
    fn missing_query_map_key_and_typed_invariants() {
        super::tests::final_chiral_bonds_missing_query_map_key_and_typed_invariants();
    }
}

#[cfg(test)]
mod final_mapping_ {
    #[test]
    fn option_gates_and_ring_error_precede_user_hook() {
        super::tests::final_mapping_option_gates_and_ring_error_precede_user_hook();
    }

    #[test]
    fn mapped_chirality_rejects_before_user_hook() {
        super::tests::final_mapping_mapped_chirality_rejects_before_user_hook();
    }

    #[test]
    fn initial_seed_suppression_and_growth_restore_checks() {
        super::tests::final_mapping_initial_seed_suppression_and_growth_restore_checks();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_core::ValenceAssignment;
    use cosmolkit_model::{
        Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondQueryPredicate, BondSpec,
        Conformer3D, CoordinateBlock, QueryAtom, QueryNode, QueryStateRef, TopologyBlock,
    };
    use cosmolkit_types::Element;

    fn q107_topology(chiral_tags: &[ChiralTag]) -> TopologyBlock {
        TopologyBlock::try_from_parts(
            chiral_tags
                .iter()
                .copied()
                .enumerate()
                .map(|(index, chiral_tag)| {
                    Atom::from_spec(
                        AtomId::new(index),
                        AtomSpec::new(Element::C).with_chiral_tag(chiral_tag),
                    )
                })
                .collect(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed atom-only topology is valid")
    }

    fn q107_target<'a>(
        topology: &'a TopologyBlock,
        coordinates: &'a CoordinateBlock,
    ) -> SearchTarget<'a> {
        SearchTarget::new(topology, coordinates, &topology.stereo_groups, None, None)
    }

    fn q104_topology(specs: impl IntoIterator<Item = AtomSpec>) -> TopologyBlock {
        TopologyBlock::try_from_parts(
            specs
                .into_iter()
                .enumerate()
                .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
                .collect(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed Q104 atom-only topology is valid")
    }

    fn q110_topology(order: BondOrder, stereo: BondStereo) -> TopologyBlock {
        TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            ],
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), order).with_stereo(stereo),
            )],
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed Q110 one-bond topology is valid")
    }

    fn q111_cycle(size: usize) -> TopologyBlock {
        let atoms = (0..size)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let bonds = (0..size)
            .map(|index| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(
                        AtomId::new(index),
                        AtomId::new((index + 1) % size),
                        BondOrder::Single,
                    ),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed Q111 cycle topology is valid")
    }

    fn q111_fused_triangles() -> TopologyBlock {
        let atoms = (0..4)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let endpoints = [(0, 1), (1, 2), (2, 0), (1, 3), (3, 0)];
        let bonds = endpoints
            .into_iter()
            .enumerate()
            .map(|(index, (begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed Q111 fused topology is valid")
    }

    fn q113_chain(atom_count: usize) -> TopologyBlock {
        let atoms = (0..atom_count)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let bonds = (1..atom_count)
            .map(|index| {
                Bond::from_spec(
                    BondId::new(index - 1),
                    BondSpec::new(
                        AtomId::new(index - 1),
                        AtomId::new(index),
                        BondOrder::Single,
                    ),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed Q113 chain topology is valid")
    }

    fn q114_topology() -> TopologyBlock {
        TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C).with_isotope(13)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::O).with_isotope(13)),
                Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::H).with_isotope(1)),
            ],
            vec![
                Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
                ),
                Bond::from_spec(
                    BondId::new(1),
                    BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Aromatic),
                ),
            ],
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed Q114 topology is valid")
    }

    #[test]
    fn q103_parameter_defaults_and_enum_names_match_source() {
        assert_eq!(
            [
                AtomComparator::AtomCompareAny,
                AtomComparator::AtomCompareElements,
                AtomComparator::AtomCompareIsotopes,
                AtomComparator::AtomCompareAnyHeavyAtom,
            ],
            [
                AtomComparator::AtomCompareAny,
                AtomComparator::default(),
                AtomComparator::AtomCompareIsotopes,
                AtomComparator::AtomCompareAnyHeavyAtom,
            ]
        );
        assert_eq!(
            [
                BondComparator::BondCompareAny,
                BondComparator::BondCompareOrder,
                BondComparator::BondCompareOrderExact,
            ],
            [
                BondComparator::BondCompareAny,
                BondComparator::default(),
                BondComparator::BondCompareOrderExact,
            ]
        );
        assert_eq!(
            [
                RingComparator::IgnoreRingFusion,
                RingComparator::PermissiveRingFusion,
                RingComparator::StrictRingFusion,
            ],
            [
                RingComparator::default(),
                RingComparator::PermissiveRingFusion,
                RingComparator::StrictRingFusion,
            ]
        );

        assert_eq!(
            McsParameters::default(),
            McsParameters {
                store_all: false,
                maximize_bonds: true,
                threshold: 1.0,
                timeout: 0,
                verbose: false,
                atom_compare_parameters: McsAtomCompareParameters {
                    match_valences: false,
                    match_chiral_tag: false,
                    match_formal_charge: false,
                    ring_matches_ring_only: false,
                    complete_rings_only: false,
                    match_isotope: false,
                    max_distance: -1.0,
                },
                bond_compare_parameters: McsBondCompareParameters {
                    ring_matches_ring_only: false,
                    complete_rings_only: false,
                    match_fused_rings: false,
                    match_fused_rings_strict: false,
                    match_stereo: false,
                },
                atom_comparator: AtomComparator::AtomCompareElements,
                bond_comparator: BondComparator::BondCompareOrder,
                initial_seed: String::new(),
            }
        );
    }

    #[test]
    fn q103_json_updates_exact_source_fields_and_precedence() {
        let mut params = McsParameters {
            verbose: true,
            initial_seed: "old".to_owned(),
            atom_compare_parameters: McsAtomCompareParameters {
                match_isotope: true,
                ..McsAtomCompareParameters::default()
            },
            ..McsParameters::default()
        };
        update_mcs_parameters_from_json(
            &mut params,
            r#"{
                "MaximizeBonds": false,
                "Threshold": 0.75,
                "Timeout": 17,
                "MatchValences": true,
                "MatchChiralTag": true,
                "MatchFormalCharge": true,
                "RingMatchesRingOnly": true,
                "AtomRingMatchesRingOnly": false,
                "BondRingMatchesRingOnly": false,
                "MaxDistance": 2.5,
                "CompleteRingsOnly": true,
                "AtomCompleteRingsOnly": true,
                "BondCompleteRingsOnly": false,
                "MatchFusedRings": true,
                "MatchFusedRingsStrict": true,
                "MatchStereo": true,
                "StoreAll": true,
                "AtomCompare": "AnyHeavy",
                "BondCompare": "OrderExact",
                "InitialSeed": "[#6]-[#8]",
                "Unknown": 12
            }"#,
        )
        .expect("fixed source-shaped JSON is valid");

        assert!(!params.maximize_bonds);
        assert_eq!(params.threshold, 0.75);
        assert_eq!(params.timeout, 17);
        assert!(params.store_all);
        assert!(params.verbose, "Verbose is a field but not a JSON key");
        assert!(params.atom_compare_parameters.match_valences);
        assert!(params.atom_compare_parameters.match_chiral_tag);
        assert!(params.atom_compare_parameters.match_formal_charge);
        assert!(!params.atom_compare_parameters.ring_matches_ring_only);
        assert!(params.atom_compare_parameters.complete_rings_only);
        assert!(params.atom_compare_parameters.match_isotope);
        assert_eq!(params.atom_compare_parameters.max_distance, 2.5);
        assert!(!params.bond_compare_parameters.ring_matches_ring_only);
        assert!(!params.bond_compare_parameters.complete_rings_only);
        assert!(params.bond_compare_parameters.match_fused_rings);
        assert!(params.bond_compare_parameters.match_fused_rings_strict);
        assert!(params.bond_compare_parameters.match_stereo);
        assert_eq!(
            params.atom_comparator,
            AtomComparator::AtomCompareAnyHeavyAtom
        );
        assert_eq!(
            params.bond_comparator,
            BondComparator::BondCompareOrderExact
        );
        assert_eq!(params.initial_seed, "[#6]-[#8]");

        update_mcs_parameters_from_json(
            &mut params,
            r#"{"AtomCompare":"unknown","BondCompare":"def"}"#,
        )
        .expect("source ignores unknown comparator names");
        assert_eq!(
            params.atom_comparator,
            AtomComparator::AtomCompareAnyHeavyAtom
        );
        assert_eq!(
            params.bond_comparator,
            BondComparator::BondCompareOrderExact
        );
        assert!(params.initial_seed.is_empty());
    }

    #[test]
    fn q103_json_reports_typed_errors_and_retains_prior_source_updates() {
        let mut params = McsParameters::default();
        assert!(matches!(
            update_mcs_parameters_from_json(&mut params, "{"),
            Err(McsParametersJsonError::InvalidJson(_))
        ));
        assert_eq!(
            update_mcs_parameters_from_json(&mut params, "[]"),
            Err(McsParametersJsonError::InvalidField { field: "root" })
        );
        assert_eq!(
            update_mcs_parameters_from_json(
                &mut params,
                r#"{"MaximizeBonds":false,"Threshold":"not-a-number"}"#,
            ),
            Err(McsParametersJsonError::InvalidField { field: "Threshold" })
        );
        assert!(!params.maximize_bonds);
        assert_eq!(params.threshold, 1.0);
        assert_eq!(
            update_mcs_parameters_from_json(&mut params, r#"{"Timeout":-1}"#),
            Err(McsParametersJsonError::InvalidField { field: "Timeout" })
        );
        assert_eq!(params.timeout, 0);
        update_mcs_parameters_from_json(&mut params, "").expect("empty JSON is a no-op");
    }

    #[test]
    fn q107_chirality_uses_directional_source_mapping_and_option_gate() {
        let left_topology = q107_topology(&[
            ChiralTag::Unspecified,
            ChiralTag::TetrahedralCw,
            ChiralTag::Other,
        ]);
        let right_topology = q107_topology(&[
            ChiralTag::Unspecified,
            ChiralTag::TetrahedralCcw,
            ChiralTag::SquarePlanar,
        ]);
        let coordinates = CoordinateBlock::default();
        let left = q107_target(&left_topology, &coordinates);
        let right = q107_target(&right_topology, &coordinates);

        assert!(check_atom_chirality(&left, 1, &right, 1).unwrap());
        assert!(!check_atom_chirality(&left, 1, &right, 0).unwrap());
        assert!(check_atom_chirality(&left, 0, &right, 1).unwrap());
        assert!(check_atom_chirality(&left, 2, &right, 0).unwrap());
        assert_eq!(
            check_atom_chirality(&left, 3, &right, 0),
            Err(McsError::AtomOutOfRange {
                side: "left",
                atom: 3,
            })
        );

        let disabled = McsAtomCompareParameters::default();
        assert!(check_atom_stereo_and_distance(&disabled, &left, 1, &right, 0).unwrap());
        let enabled = McsAtomCompareParameters {
            match_chiral_tag: true,
            ..McsAtomCompareParameters::default()
        };
        assert!(!check_atom_stereo_and_distance(&enabled, &left, 1, &right, 0).unwrap());
    }

    #[test]
    fn q107_distance_uses_first_conformer_mapping_and_inclusive_threshold() {
        let topology = q107_topology(&[ChiralTag::Unspecified, ChiralTag::Unspecified]);
        let left_coordinates = CoordinateBlock {
            conformers_3d: vec![
                Conformer3D::new(7, vec![[9.0, 9.0, 9.0], [0.0, 0.0, 0.0]], true),
                Conformer3D::new(9, vec![[100.0, 0.0, 0.0], [100.0, 0.0, 0.0]], true),
            ],
            ..CoordinateBlock::default()
        };
        let right_coordinates = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(
                3,
                vec![[3.0, 4.0, 0.0], [50.0, 50.0, 50.0]],
                true,
            )],
            ..CoordinateBlock::default()
        };
        let left = q107_target(&topology, &left_coordinates);
        let right = q107_target(&topology, &right_coordinates);
        let exact = McsAtomCompareParameters {
            max_distance: 5.0,
            ..McsAtomCompareParameters::default()
        };
        assert!(check_atom_distance(&exact, &left, 1, &right, 0).unwrap());
        assert!(check_atom_stereo_and_distance(&exact, &left, 1, &right, 0).unwrap());
        let short = McsAtomCompareParameters {
            max_distance: 4.999,
            ..McsAtomCompareParameters::default()
        };
        assert!(!check_atom_distance(&short, &left, 1, &right, 0).unwrap());
    }

    #[test]
    fn q107_distance_gate_preserves_missing_coordinate_errors() {
        let topology = q107_topology(&[ChiralTag::Unspecified, ChiralTag::Unspecified]);
        let coordinates = CoordinateBlock::default();
        let target = q107_target(&topology, &coordinates);
        let disabled = McsAtomCompareParameters::default();
        assert!(check_atom_stereo_and_distance(&disabled, &target, 0, &target, 0).unwrap());

        let enabled = McsAtomCompareParameters {
            max_distance: 1.0,
            ..McsAtomCompareParameters::default()
        };
        assert_eq!(
            check_atom_stereo_and_distance(&enabled, &target, 0, &target, 0),
            Err(McsError::MissingConformer { side: "left" })
        );

        let short_coordinates = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(0, vec![[0.0, 0.0, 0.0]], true)],
            ..CoordinateBlock::default()
        };
        let short = q107_target(&topology, &short_coordinates);
        assert_eq!(
            check_atom_distance(&enabled, &short, 1, &short, 0),
            Err(McsError::CoordinateOutOfRange {
                side: "left",
                atom: 1,
            })
        );
    }

    #[test]
    fn q104_any_and_any_heavy_preserve_hydrogen_and_isotope_rules() {
        let left_topology = q104_topology([
            AtomSpec::new(Element::H).with_isotope(1),
            AtomSpec::new(Element::C)
                .with_isotope(13)
                .with_formal_charge(1),
        ]);
        let right_topology = q104_topology([
            AtomSpec::new(Element::H).with_isotope(2),
            AtomSpec::new(Element::O).with_isotope(18),
        ]);
        let coordinates = CoordinateBlock::default();
        let left = q107_target(&left_topology, &coordinates);
        let right = q107_target(&right_topology, &coordinates);
        let params = McsAtomCompareParameters::default();

        assert!(mcs_atom_compare_any(&params, &left, 0, &right, 1).unwrap());
        assert!(mcs_atom_compare_any(&params, &left, 1, &right, 0).unwrap());
        assert!(mcs_atom_compare_any_heavy(&params, &left, 0, &right, 0).unwrap());
        assert!(!mcs_atom_compare_any_heavy(&params, &left, 0, &right, 1).unwrap());
        assert!(mcs_atom_compare_any_heavy(&params, &left, 1, &right, 1).unwrap());
    }

    #[test]
    fn q104_any_applies_charge_flag_and_ignores_query_predicate_identity() {
        let left_topology = q104_topology([
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::C).with_formal_charge(1),
        ]);
        let right_topology = q104_topology([AtomSpec::new(Element::H), AtomSpec::new(Element::O)]);
        let left_rows = vec![
            QueryAtom::from_parts(
                left_topology.atoms[0].clone(),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
            ),
            QueryAtom::from_parts(
                left_topology.atoms[1].clone(),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
            ),
        ];
        let right_rows = vec![
            QueryAtom::from_parts(
                right_topology.atoms[0].clone(),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            ),
            QueryAtom::from_parts(
                right_topology.atoms[1].clone(),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
            ),
        ];
        let left_state = QueryStateRef::try_for_topology(&left_rows, &[], &left_topology).unwrap();
        let right_state =
            QueryStateRef::try_for_topology(&right_rows, &[], &right_topology).unwrap();
        let coordinates = CoordinateBlock::default();
        let left = q107_target(&left_topology, &coordinates)
            .try_with_query_state(left_state)
            .unwrap();
        let right = q107_target(&right_topology, &coordinates)
            .try_with_query_state(right_state)
            .unwrap();

        let default = McsAtomCompareParameters::default();
        assert!(mcs_atom_compare_any(&default, &left, 1, &right, 1).unwrap());
        assert!(mcs_atom_compare_any_heavy(&default, &left, 1, &right, 1).unwrap());
        assert!(!mcs_atom_compare_any_heavy(&default, &left, 0, &right, 1).unwrap());

        let charge = McsAtomCompareParameters {
            match_formal_charge: true,
            ..McsAtomCompareParameters::default()
        };
        assert!(!mcs_atom_compare_any(&charge, &left, 1, &right, 1).unwrap());
        assert!(mcs_atom_compare_any(&charge, &left, 0, &right, 0).unwrap());
    }

    #[test]
    fn q104_ring_option_requires_and_uses_detached_ring_assignments() {
        let topology = q104_topology([AtomSpec::new(Element::C)]);
        let coordinates = CoordinateBlock::default();
        let plain = q107_target(&topology, &coordinates);
        let params = McsAtomCompareParameters {
            ring_matches_ring_only: true,
            ..McsAtomCompareParameters::default()
        };
        assert_eq!(
            mcs_atom_compare_any(&params, &plain, 0, &plain, 0),
            Err(McsError::MissingRingInfo { side: "left" })
        );

        let rings = cosmolkit_core::fast_find_rings(&topology).unwrap();
        let with_rings = SearchTarget::new(
            &topology,
            &coordinates,
            &topology.stereo_groups,
            Some(&rings),
            None,
        );
        assert!(mcs_atom_compare_any(&params, &with_rings, 0, &with_rings, 0).unwrap());
    }

    #[test]
    fn q105_elements_require_carrier_atomic_identity_before_optional_state() {
        let carbon = q104_topology([AtomSpec::new(Element::C)]);
        let oxygen = q104_topology([AtomSpec::new(Element::O)]);
        let coordinates = CoordinateBlock::default();
        let carbon_target = q107_target(&carbon, &coordinates);
        let oxygen_target = q107_target(&oxygen, &coordinates);

        let valence_enabled = McsAtomCompareParameters {
            match_valences: true,
            ..McsAtomCompareParameters::default()
        };
        assert!(
            !mcs_atom_compare_elements(&valence_enabled, &carbon_target, 0, &oxygen_target, 0,)
                .unwrap()
        );
        assert!(
            mcs_atom_compare_elements(
                &McsAtomCompareParameters::default(),
                &carbon_target,
                0,
                &carbon_target,
                0,
            )
            .unwrap()
        );
    }

    #[test]
    fn q105_elements_apply_total_valence_and_charge_only_when_enabled() {
        let positive = q104_topology([AtomSpec::new(Element::C).with_formal_charge(1)]);
        let neutral = q104_topology([AtomSpec::new(Element::C)]);
        let coordinates = CoordinateBlock::default();
        let left_valence = ValenceAssignment {
            explicit_valence: vec![3],
            implicit_hydrogens: vec![1],
        };
        let right_equal_total = ValenceAssignment {
            explicit_valence: vec![2],
            implicit_hydrogens: vec![2],
        };
        let right_different_total = ValenceAssignment {
            explicit_valence: vec![2],
            implicit_hydrogens: vec![1],
        };
        let left = SearchTarget::new(
            &positive,
            &coordinates,
            &positive.stereo_groups,
            None,
            Some(&left_valence),
        );
        let equal = SearchTarget::new(
            &neutral,
            &coordinates,
            &neutral.stereo_groups,
            None,
            Some(&right_equal_total),
        );
        let different = SearchTarget::new(
            &neutral,
            &coordinates,
            &neutral.stereo_groups,
            None,
            Some(&right_different_total),
        );

        let default = McsAtomCompareParameters::default();
        assert!(mcs_atom_compare_elements(&default, &left, 0, &different, 0).unwrap());

        let valence = McsAtomCompareParameters {
            match_valences: true,
            ..McsAtomCompareParameters::default()
        };
        assert!(mcs_atom_compare_elements(&valence, &left, 0, &equal, 0).unwrap());
        assert!(!mcs_atom_compare_elements(&valence, &left, 0, &different, 0).unwrap());

        let charge = McsAtomCompareParameters {
            match_formal_charge: true,
            ..McsAtomCompareParameters::default()
        };
        assert!(!mcs_atom_compare_elements(&charge, &left, 0, &equal, 0).unwrap());
    }

    #[test]
    fn q105_valence_option_preserves_missing_and_short_assignment_errors() {
        let topology = q104_topology([AtomSpec::new(Element::C)]);
        let coordinates = CoordinateBlock::default();
        let missing = q107_target(&topology, &coordinates);
        let params = McsAtomCompareParameters {
            match_valences: true,
            ..McsAtomCompareParameters::default()
        };
        assert_eq!(
            mcs_atom_compare_elements(&params, &missing, 0, &missing, 0),
            Err(McsError::MissingValence { side: "left" })
        );

        let short_assignment = ValenceAssignment {
            explicit_valence: Vec::new(),
            implicit_hydrogens: Vec::new(),
        };
        let short = SearchTarget::new(
            &topology,
            &coordinates,
            &topology.stereo_groups,
            None,
            Some(&short_assignment),
        );
        assert_eq!(
            mcs_atom_compare_elements(&params, &short, 0, &short, 0),
            Err(McsError::ValenceOutOfRange {
                side: "left",
                atom: 0,
            })
        );
    }

    #[test]
    fn q106_isotopes_treat_zero_as_unspecified_and_ignore_element_identity() {
        let unspecified = q104_topology([AtomSpec::new(Element::C)]);
        let explicit_zero = q104_topology([AtomSpec::new(Element::O).with_isotope(0)]);
        let isotope_13_n = q104_topology([AtomSpec::new(Element::N).with_isotope(13)]);
        let isotope_13_cl = q104_topology([AtomSpec::new(Element::CL).with_isotope(13)]);
        let isotope_14 = q104_topology([AtomSpec::new(Element::N).with_isotope(14)]);
        let coordinates = CoordinateBlock::default();
        let params = McsAtomCompareParameters::default();

        assert!(
            mcs_atom_compare_isotopes(
                &params,
                &q107_target(&unspecified, &coordinates),
                0,
                &q107_target(&explicit_zero, &coordinates),
                0,
            )
            .unwrap()
        );
        assert!(
            mcs_atom_compare_isotopes(
                &params,
                &q107_target(&isotope_13_n, &coordinates),
                0,
                &q107_target(&isotope_13_cl, &coordinates),
                0,
            )
            .unwrap()
        );
        assert!(
            !mcs_atom_compare_isotopes(
                &params,
                &q107_target(&isotope_13_n, &coordinates),
                0,
                &q107_target(&isotope_14, &coordinates),
                0,
            )
            .unwrap()
        );
    }

    #[test]
    fn q106_isotopes_apply_charge_chirality_distance_and_ring_checks() {
        let left_topology = q104_topology([AtomSpec::new(Element::C)
            .with_isotope(13)
            .with_formal_charge(1)
            .with_chiral_tag(ChiralTag::TetrahedralCw)]);
        let right_topology = q104_topology([AtomSpec::new(Element::O).with_isotope(13)]);
        let left_coordinates = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(0, vec![[0.0, 0.0, 0.0]], true)],
            ..CoordinateBlock::default()
        };
        let right_coordinates = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(0, vec![[2.0, 0.0, 0.0]], true)],
            ..CoordinateBlock::default()
        };
        let left = q107_target(&left_topology, &left_coordinates);
        let right = q107_target(&right_topology, &right_coordinates);
        assert!(
            mcs_atom_compare_isotopes(&McsAtomCompareParameters::default(), &left, 0, &right, 0,)
                .unwrap()
        );
        assert!(
            !mcs_atom_compare_isotopes(
                &McsAtomCompareParameters {
                    match_formal_charge: true,
                    ..McsAtomCompareParameters::default()
                },
                &left,
                0,
                &right,
                0,
            )
            .unwrap()
        );
        assert!(
            !mcs_atom_compare_isotopes(
                &McsAtomCompareParameters {
                    match_chiral_tag: true,
                    ..McsAtomCompareParameters::default()
                },
                &left,
                0,
                &right,
                0,
            )
            .unwrap()
        );
        assert!(
            !mcs_atom_compare_isotopes(
                &McsAtomCompareParameters {
                    max_distance: 1.0,
                    ..McsAtomCompareParameters::default()
                },
                &left,
                0,
                &right,
                0,
            )
            .unwrap()
        );

        let left_rings = cosmolkit_core::fast_find_rings(&left_topology).unwrap();
        let right_rings = cosmolkit_core::fast_find_rings(&right_topology).unwrap();
        let left_with_rings = SearchTarget::new(
            &left_topology,
            &left_coordinates,
            &left_topology.stereo_groups,
            Some(&left_rings),
            None,
        );
        let right_with_rings = SearchTarget::new(
            &right_topology,
            &right_coordinates,
            &right_topology.stereo_groups,
            Some(&right_rings),
            None,
        );
        assert!(
            mcs_atom_compare_isotopes(
                &McsAtomCompareParameters {
                    ring_matches_ring_only: true,
                    ..McsAtomCompareParameters::default()
                },
                &left_with_rings,
                0,
                &right_with_rings,
                0,
            )
            .unwrap()
        );
    }

    #[test]
    fn q110_double_bond_stereo_preserves_source_enum_mapping() {
        let z = q110_topology(BondOrder::Double, BondStereo::Z);
        let e = q110_topology(BondOrder::Double, BondStereo::E);
        let any = q110_topology(BondOrder::Double, BondStereo::Any);
        let none = q110_topology(BondOrder::Double, BondStereo::None);
        let coordinates = CoordinateBlock::default();
        let z = q107_target(&z, &coordinates);
        let e = q107_target(&e, &coordinates);
        let any = q107_target(&any, &coordinates);
        let none = q107_target(&none, &coordinates);
        let params = McsBondCompareParameters::default();

        assert!(check_bond_stereo(&params, &z, 0, &z, 0).unwrap());
        assert!(!check_bond_stereo(&params, &z, 0, &e, 0).unwrap());
        assert!(!check_bond_stereo(&params, &z, 0, &any, 0).unwrap());
        assert!(!check_bond_stereo(&params, &none, 0, &any, 0).unwrap());
        assert!(check_bond_stereo(&params, &any, 0, &any, 0).unwrap());
    }

    #[test]
    fn q110_helper_ignores_option_and_non_double_stereo_as_source_does() {
        let double_z = q110_topology(BondOrder::Double, BondStereo::Z);
        let single_e = q110_topology(BondOrder::Single, BondStereo::E);
        let coordinates = CoordinateBlock::default();
        let double_z = q107_target(&double_z, &coordinates);
        let single_e = q107_target(&single_e, &coordinates);
        let disabled = McsBondCompareParameters::default();
        let enabled = McsBondCompareParameters {
            match_stereo: true,
            ..McsBondCompareParameters::default()
        };

        assert!(check_bond_stereo(&disabled, &double_z, 0, &single_e, 0).unwrap());
        assert!(check_bond_stereo(&enabled, &double_z, 0, &single_e, 0).unwrap());
        assert_eq!(
            check_bond_stereo(&enabled, &double_z, 1, &single_e, 0),
            Err(McsError::BondOutOfRange {
                side: "left",
                bond: 1,
            })
        );
    }

    #[test]
    fn q111_ring_pair_accepts_equal_sizes_and_rejects_plain_unequal_sizes() {
        let triangle = q111_cycle(3);
        let square = q111_cycle(4);
        let triangle_rings = cosmolkit_core::fast_find_rings(&triangle).unwrap();
        let square_rings = cosmolkit_core::fast_find_rings(&square).unwrap();
        let coordinates = CoordinateBlock::default();
        let triangle = SearchTarget::new(
            &triangle,
            &coordinates,
            &triangle.stereo_groups,
            Some(&triangle_rings),
            None,
        );
        let square = SearchTarget::new(
            &square,
            &coordinates,
            &square.stereo_groups,
            Some(&square_rings),
            None,
        );
        let params = McsBondCompareParameters::default();

        assert!(have_pair_of_compatible_rings(&params, &triangle, 0, &triangle, 0).unwrap());
        assert!(!have_pair_of_compatible_rings(&params, &triangle, 0, &square, 0).unwrap());
    }

    pub(super) fn final_ring_pass1_zero_one_and_multiple_memberships_use_original_bond_ids() {
        let coordinates = CoordinateBlock::default();
        let chain = q113_chain(3);
        let chain_rings = cosmolkit_core::fast_find_rings(&chain).unwrap();
        let chain_target = SearchTarget::new(
            &chain,
            &coordinates,
            &chain.stereo_groups,
            Some(&chain_rings),
            None,
        );
        let mut chain_counts = McsRingBondCountVect::new(&chain_target, "query").unwrap();
        chain_counts
            .set_mcs_bond_bits_pass1(0, 1, &[0, 1], &[2, 1, 0])
            .unwrap();
        assert_eq!(chain_counts.rings.len(), 0);
        assert!(chain_counts.mcs_bonds.test(1));
        assert!(!chain_counts.mcs_bonds.test(0));

        let triangle = q111_cycle(3);
        let triangle_rings = cosmolkit_core::fast_find_rings(&triangle).unwrap();
        let triangle_target = SearchTarget::new(
            &triangle,
            &coordinates,
            &triangle.stereo_groups,
            Some(&triangle_rings),
            None,
        );
        let mut triangle_counts = McsRingBondCountVect::new(&triangle_target, "query").unwrap();
        triangle_counts
            .set_mcs_bond_bits_pass1(0, 1, &[0, 1], &[0, 1, 2])
            .unwrap();
        assert!(triangle_counts.rings[0].all.test(0));
        assert!(triangle_counts.rings[0].nonfused.test(0));
        assert!(!triangle_counts.rings[0].fused.test(0));

        let fused = q111_fused_triangles();
        let fused_rings = cosmolkit_core::fast_find_rings(&fused).unwrap();
        let fused_target = SearchTarget::new(
            &fused,
            &coordinates,
            &fused.stereo_groups,
            Some(&fused_rings),
            None,
        );
        let mut fused_counts = McsRingBondCountVect::new(&fused_target, "query").unwrap();
        fused_counts
            .set_mcs_bond_bits_pass1(0, 1, &[0, 1], &[0, 1, 2, 3])
            .unwrap();
        assert_eq!(fused_counts.rings.len(), 2);
        assert!(fused_counts.mcs_bonds.test(0));
        for ring in &fused_counts.rings {
            assert!(ring.all.test(0));
            assert!(ring.fused.test(0));
            assert!(!ring.nonfused.test(0));
        }
    }

    pub(super) fn final_ring_pass1_partial_and_nonidentity_mapping_stays_on_mapped_edge() {
        let triangle = q111_cycle(3);
        let rings = cosmolkit_core::fast_find_rings(&triangle).unwrap();
        let coordinates = CoordinateBlock::default();
        let target = SearchTarget::new(
            &triangle,
            &coordinates,
            &triangle.stereo_groups,
            Some(&rings),
            None,
        );
        let mut counts = McsRingBondCountVect::new(&target, "target").unwrap();
        counts
            .set_mcs_bond_bits_pass1(0, 1, &[1, 0], &[0, 1, 2])
            .unwrap();
        assert_eq!(counts.mcs_bonds.count(), 1);
        assert!(counts.mcs_bonds.test(0));
        assert_eq!(counts.rings[0].all.count(), 1);
    }

    pub(super) fn final_ring_pass1_reports_missing_state_and_invalid_mapped_edges() {
        let chain = q113_chain(3);
        let rings = cosmolkit_core::fast_find_rings(&chain).unwrap();
        let coordinates = CoordinateBlock::default();
        let without_rings =
            SearchTarget::new(&chain, &coordinates, &chain.stereo_groups, None, None);
        assert!(matches!(
            McsRingBondCountVect::new(&without_rings, "query"),
            Err(McsCandidateMatchError::State(McsError::MissingRingInfo {
                side: "query"
            }))
        ));
        let target = SearchTarget::new(
            &chain,
            &coordinates,
            &chain.stereo_groups,
            Some(&rings),
            None,
        );
        let mut counts = McsRingBondCountVect::new(&target, "target").unwrap();
        assert!(matches!(
            counts.set_mcs_bond_bits_pass1(0, 1, &[0, 1], &[0, 2, 1]),
            Err(McsCandidateMatchError::MappedBondMissing {
                side: "target",
                begin: 0,
                end: 2,
            })
        ));
        assert!(matches!(
            counts.set_mcs_bond_bits_pass1(0, 1, &[0, 3], &[0, 1, 2]),
            Err(McsCandidateMatchError::State(McsError::AtomOutOfRange {
                side: "final graph vertex",
                atom: 3,
            }))
        ));
        assert_eq!(counts.mcs_bonds.count(), 0);
    }

    fn final_ring_refine_select(
        counts: &mut McsRingBondCountVect<'_, '_>,
        topology: &TopologyBlock,
        bonds: &[usize],
    ) {
        let identity = (0..topology.atoms.len()).collect::<Vec<_>>();
        for &index in bonds {
            let bond = &topology.bonds[index];
            counts
                .set_mcs_bond_bits_pass1(
                    bond.begin().index(),
                    bond.end().index(),
                    &identity,
                    &identity,
                )
                .unwrap();
        }
    }

    pub(super) fn final_ring_refine_pass2_removes_incomplete_fusion_and_preserves_snapshot() {
        let fused = q111_fused_triangles();
        let rings = cosmolkit_core::fast_find_rings(&fused).unwrap();
        let coordinates = CoordinateBlock::default();
        let target = SearchTarget::new(
            &fused,
            &coordinates,
            &fused.stereo_groups,
            Some(&rings),
            None,
        );
        let shared = (0..fused.bonds.len())
            .find(|&bond| rings.num_bond_rings(BondId::new(bond)) == 2)
            .unwrap();
        let mut counts = McsRingBondCountVect::new(&target, "query").unwrap();
        final_ring_refine_select(&mut counts, &fused, &[shared]);
        counts.set_mcs_bond_bits_pass2().unwrap();
        for row in &counts.rings {
            assert_eq!(row.fused_count_pass1, 1);
            assert_eq!(row.nonfused_count_pass1, 0);
            assert!(!row.fused.test(shared));
            assert!(!row.nonfused.test(shared));
            assert!(row.all.test(shared));
        }
        assert!(counts.is_ring_fusion_honored().unwrap());
    }

    pub(super) fn final_ring_refine_pass2_reclassifies_one_survivor_and_retains_two() {
        let fused = q111_fused_triangles();
        let rings = cosmolkit_core::fast_find_rings(&fused).unwrap();
        let coordinates = CoordinateBlock::default();
        let target = SearchTarget::new(
            &fused,
            &coordinates,
            &fused.stereo_groups,
            Some(&rings),
            None,
        );
        let shared = (0..fused.bonds.len())
            .find(|&bond| rings.num_bond_rings(BondId::new(bond)) == 2)
            .unwrap();
        let unique = rings
            .bond_rings()
            .iter()
            .map(|row| {
                row.iter()
                    .find(|bond| bond.index() != shared)
                    .unwrap()
                    .index()
            })
            .collect::<Vec<_>>();
        let mut one = McsRingBondCountVect::new(&target, "query").unwrap();
        final_ring_refine_select(&mut one, &fused, &[shared, unique[0]]);
        one.set_mcs_bond_bits_pass2().unwrap();
        assert_eq!(one.rings[0].nonfused_count_pass1, 1);
        assert_eq!(one.rings[0].fused_count_pass1, 1);
        assert_eq!(one.rings[1].nonfused_count_pass1, 0);
        assert_eq!(one.rings[1].fused_count_pass1, 1);
        assert!(one.rings[0].nonfused.test(shared));
        assert!(!one.rings[0].fused.test(shared));
        assert!(!one.rings[1].fused.test(shared));

        let mut both = McsRingBondCountVect::new(&target, "query").unwrap();
        final_ring_refine_select(&mut both, &fused, &[shared, unique[0], unique[1]]);
        both.set_mcs_bond_bits_pass2().unwrap();
        for row in &both.rings {
            assert_eq!(row.nonfused_count_pass1, 1);
            assert_eq!(row.fused_count_pass1, 1);
            assert!(row.fused.test(shared));
            assert!(!row.nonfused.test(shared));
        }
    }

    pub(super) fn final_ring_refine_honor_check_distinguishes_none_all_partial_and_missing_fused() {
        let fused = q111_fused_triangles();
        let rings = cosmolkit_core::fast_find_rings(&fused).unwrap();
        let coordinates = CoordinateBlock::default();
        let target = SearchTarget::new(
            &fused,
            &coordinates,
            &fused.stereo_groups,
            Some(&rings),
            None,
        );
        let unique = rings.bond_rings()[0]
            .iter()
            .filter(|bond| rings.num_bond_rings(**bond) == 1)
            .map(|bond| bond.index())
            .collect::<Vec<_>>();
        assert_eq!(unique.len(), 2);

        let mut none = McsRingBondCountVect::new(&target, "query").unwrap();
        none.set_mcs_bond_bits_pass2().unwrap();
        assert!(none.is_ring_fusion_honored().unwrap());
        let mut all = McsRingBondCountVect::new(&target, "query").unwrap();
        final_ring_refine_select(
            &mut all,
            &fused,
            &(0..fused.bonds.len()).collect::<Vec<_>>(),
        );
        all.set_mcs_bond_bits_pass2().unwrap();
        assert!(all.is_ring_fusion_honored().unwrap());
        let mut partial = McsRingBondCountVect::new(&target, "query").unwrap();
        final_ring_refine_select(&mut partial, &fused, &unique[..1]);
        partial.set_mcs_bond_bits_pass2().unwrap();
        assert!(partial.is_ring_fusion_honored().unwrap());
        let mut missing_fused = McsRingBondCountVect::new(&target, "query").unwrap();
        final_ring_refine_select(&mut missing_fused, &fused, &unique);
        missing_fused.set_mcs_bond_bits_pass2().unwrap();
        assert_eq!(missing_fused.rings[0].all.count(), 2);
        assert!(!missing_fused.is_ring_fusion_honored().unwrap());
    }

    fn final_ring_fusion_path_graph() -> McsSeedTopology {
        McsSeedTopology {
            source_atoms: vec![1, 2, 0],
            bonds: vec![
                McsSeedTopologyBond {
                    source_bond: 1,
                    begin_seed_atom: 0,
                    end_seed_atom: 1,
                },
                McsSeedTopologyBond {
                    source_bond: 2,
                    begin_seed_atom: 1,
                    end_seed_atom: 2,
                },
            ],
        }
    }

    fn final_ring_fusion_params(strict: bool) -> McsParameters {
        let mut params = McsParameters::default();
        params.bond_compare_parameters.match_fused_rings = true;
        params.bond_compare_parameters.match_fused_rings_strict = strict;
        params
    }

    pub(super) fn final_ring_fusion_neither_one_and_both_molecule_checks_follow_mode() {
        let fused = q111_fused_triangles();
        let chain = q113_chain(4);
        let fused_rings = cosmolkit_core::fast_find_rings(&fused).unwrap();
        let chain_rings = cosmolkit_core::fast_find_rings(&chain).unwrap();
        let coordinates = CoordinateBlock::default();
        let fused_target = SearchTarget::new(
            &fused,
            &coordinates,
            &fused.stereo_groups,
            Some(&fused_rings),
            None,
        );
        let chain_target = SearchTarget::new(
            &chain,
            &coordinates,
            &chain.stereo_groups,
            Some(&chain_rings),
            None,
        );
        let graph = final_ring_fusion_path_graph();
        let strict = final_ring_fusion_params(true);
        let permissive = final_ring_fusion_params(false);
        let fused_vertices = [0, 1, 2, 3];
        let chain_vertices = [0, 1, 2, 3];
        assert_eq!(
            mcs_ring_fusion_check(
                &[0, 1, 2],
                &[1, 2, 0],
                &fused_target,
                &graph,
                &fused_target,
                &fused_vertices,
                &permissive,
            ),
            Ok(false)
        );
        assert_eq!(
            mcs_ring_fusion_check(
                &[0, 1, 2],
                &[0, 1, 2],
                &fused_target,
                &graph,
                &chain_target,
                &chain_vertices,
                &permissive,
            ),
            Ok(true)
        );
        assert_eq!(
            mcs_ring_fusion_check(
                &[0, 1, 2],
                &[0, 1, 2],
                &fused_target,
                &graph,
                &chain_target,
                &chain_vertices,
                &strict,
            ),
            Ok(false)
        );

        let one_edge = McsSeedTopology {
            source_atoms: graph.source_atoms.clone(),
            bonds: graph.bonds[..1].to_vec(),
        };
        assert_eq!(
            mcs_ring_fusion_check(
                &[0, 1, 2],
                &[1, 2, 0],
                &fused_target,
                &one_edge,
                &fused_target,
                &fused_vertices,
                &strict,
            ),
            Ok(true)
        );
    }

    pub(super) fn final_ring_fusion_target_smaller_shortcut_precedes_required_ring_state() {
        let fused = q111_fused_triangles();
        let chain = q113_chain(2);
        let coordinates = CoordinateBlock::default();
        let query = SearchTarget::new(&fused, &coordinates, &fused.stereo_groups, None, None);
        let target = SearchTarget::new(&chain, &coordinates, &chain.stereo_groups, None, None);
        assert_eq!(
            mcs_ring_fusion_check(
                &[0, 1, 2],
                &[0, 1, 0],
                &query,
                &final_ring_fusion_path_graph(),
                &target,
                &[0, 1],
                &final_ring_fusion_params(true),
            ),
            Ok(true)
        );
    }

    pub(super) fn final_ring_fusion_distinct_query_and_target_vertex_identities_are_preserved() {
        let fused = q111_fused_triangles();
        let chain = q113_chain(4);
        let fused_rings = cosmolkit_core::fast_find_rings(&fused).unwrap();
        let chain_rings = cosmolkit_core::fast_find_rings(&chain).unwrap();
        let coordinates = CoordinateBlock::default();
        let query = SearchTarget::new(
            &chain,
            &coordinates,
            &chain.stereo_groups,
            Some(&chain_rings),
            None,
        );
        let target = SearchTarget::new(
            &fused,
            &coordinates,
            &fused.stereo_groups,
            Some(&fused_rings),
            None,
        );
        let graph = McsSeedTopology {
            source_atoms: vec![0, 1, 2],
            bonds: vec![
                McsSeedTopologyBond {
                    source_bond: 0,
                    begin_seed_atom: 0,
                    end_seed_atom: 1,
                },
                McsSeedTopologyBond {
                    source_bond: 1,
                    begin_seed_atom: 1,
                    end_seed_atom: 2,
                },
            ],
        };
        let target_vertices = [2, 1, 0, 3];
        assert_eq!(
            mcs_ring_fusion_check(
                &[0, 1, 2],
                &[1, 0, 2],
                &query,
                &graph,
                &target,
                &target_vertices,
                &final_ring_fusion_params(false),
            ),
            Ok(true)
        );
        assert_eq!(
            mcs_ring_fusion_check(
                &[0, 1, 2],
                &[1, 0, 2],
                &query,
                &graph,
                &target,
                &target_vertices,
                &final_ring_fusion_params(true),
            ),
            Ok(false)
        );
    }

    fn final_chiral_star(atom_count: usize, center_tag: ChiralTag) -> TopologyBlock {
        let atoms = (0..atom_count)
            .map(|index| {
                let tag = if index == 0 {
                    center_tag
                } else {
                    ChiralTag::Unspecified
                };
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(Element::C).with_chiral_tag(tag),
                )
            })
            .collect();
        let bonds = (1..atom_count)
            .map(|end| {
                Bond::from_spec(
                    BondId::new(end - 1),
                    BondSpec::new(AtomId::new(0), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed chiral star topology is valid")
    }

    fn final_chiral_star_graph(atom_count: usize) -> McsSeedTopology {
        McsSeedTopology {
            source_atoms: (0..atom_count).collect(),
            bonds: (1..atom_count)
                .map(|end| McsSeedTopologyBond {
                    source_bond: end - 1,
                    begin_seed_atom: 0,
                    end_seed_atom: end,
                })
                .collect(),
        }
    }

    pub(super) fn final_chiral_atoms_query_skip_and_target_tag_and_degree_branches() {
        let coordinates = CoordinateBlock::default();
        let plain = final_chiral_star(4, ChiralTag::Unspecified);
        let tagged = final_chiral_star(4, ChiralTag::TetrahedralCw);
        let plain_target = q107_target(&plain, &coordinates);
        let tagged_target = q107_target(&tagged, &coordinates);
        let graph = final_chiral_star_graph(4);
        let identity = [0, 1, 2, 3];
        assert_eq!(
            mcs_final_tetrahedral_check(
                &identity,
                &identity,
                &plain_target,
                &graph,
                &tagged_target,
                &identity,
            ),
            Ok(true)
        );
        assert_eq!(
            mcs_final_tetrahedral_check(
                &identity,
                &identity,
                &tagged_target,
                &graph,
                &plain_target,
                &identity,
            ),
            Ok(false)
        );
        let low = final_chiral_star(3, ChiralTag::TetrahedralCw);
        let low_target = q107_target(&low, &coordinates);
        assert_eq!(
            mcs_final_tetrahedral_check(
                &[0, 1, 2],
                &[0, 1, 2],
                &low_target,
                &final_chiral_star_graph(3),
                &plain_target,
                &identity,
            ),
            Ok(true)
        );
        let four = final_chiral_star(5, ChiralTag::TetrahedralCw);
        let four_target = q107_target(&four, &coordinates);
        assert_eq!(
            mcs_final_tetrahedral_check(
                &[0, 1, 2, 3, 4],
                &[0, 1, 2, 3, 0],
                &four_target,
                &final_chiral_star_graph(5),
                &tagged_target,
                &identity,
            ),
            Ok(false)
        );
    }

    pub(super) fn final_chiral_atoms_mapped_neighbor_permutations_use_tag_parity() {
        let coordinates = CoordinateBlock::default();
        let cw = final_chiral_star(4, ChiralTag::TetrahedralCw);
        let ccw = final_chiral_star(4, ChiralTag::TetrahedralCcw);
        let cw_target = q107_target(&cw, &coordinates);
        let ccw_target = q107_target(&ccw, &coordinates);
        let graph = final_chiral_star_graph(4);
        let identity = [0, 1, 2, 3];
        let odd = [0, 2, 1, 3];
        assert_eq!(
            mcs_final_tetrahedral_check(
                &identity, &identity, &cw_target, &graph, &cw_target, &identity,
            ),
            Ok(true)
        );
        assert_eq!(
            mcs_final_tetrahedral_check(
                &identity,
                &identity,
                &cw_target,
                &graph,
                &ccw_target,
                &identity,
            ),
            Ok(false)
        );
        assert_eq!(
            mcs_final_tetrahedral_check(&identity, &odd, &cw_target, &graph, &cw_target, &identity,),
            Ok(false)
        );
        assert_eq!(
            mcs_final_tetrahedral_check(
                &identity,
                &odd,
                &cw_target,
                &graph,
                &ccw_target,
                &identity,
            ),
            Ok(true)
        );
    }

    pub(super) fn final_chiral_atoms_repeated_missing_ligands_and_nonidentity_vertices() {
        let coordinates = CoordinateBlock::default();
        let topology = final_chiral_star(5, ChiralTag::TetrahedralCw);
        let target = q107_target(&topology, &coordinates);
        let partial = McsSeedTopology {
            source_atoms: vec![0, 1, 2],
            bonds: vec![
                McsSeedTopologyBond {
                    source_bond: 0,
                    begin_seed_atom: 0,
                    end_seed_atom: 1,
                },
                McsSeedTopologyBond {
                    source_bond: 1,
                    begin_seed_atom: 0,
                    end_seed_atom: 2,
                },
            ],
        };
        assert_eq!(
            mcs_final_tetrahedral_check(
                &[0, 1, 2],
                &[0, 1, 2],
                &target,
                &partial,
                &target,
                &[0, 1, 2, 3, 4],
            ),
            Ok(true)
        );
        let graph = McsSeedTopology {
            source_atoms: vec![2, 0, 1, 3, 4],
            bonds: vec![
                McsSeedTopologyBond {
                    source_bond: 0,
                    begin_seed_atom: 1,
                    end_seed_atom: 2,
                },
                McsSeedTopologyBond {
                    source_bond: 1,
                    begin_seed_atom: 1,
                    end_seed_atom: 0,
                },
                McsSeedTopologyBond {
                    source_bond: 2,
                    begin_seed_atom: 1,
                    end_seed_atom: 3,
                },
                McsSeedTopologyBond {
                    source_bond: 3,
                    begin_seed_atom: 1,
                    end_seed_atom: 4,
                },
            ],
        };
        let mapping = [1, 2, 0, 3, 4];
        assert_eq!(
            mcs_final_tetrahedral_check(
                &mapping,
                &mapping,
                &target,
                &graph,
                &target,
                &graph.source_atoms,
            ),
            Ok(true)
        );
        assert!(matches!(
            mcs_final_tetrahedral_check(
                &[1, 2, 0, 3, 9],
                &mapping,
                &target,
                &graph,
                &target,
                &graph.source_atoms,
            ),
            Err(McsCandidateMatchError::State(McsError::AtomOutOfRange {
                side: "query graph vertex",
                atom: 9
            }))
        ));
    }

    fn final_chiral_double(
        order: BondOrder,
        stereo: BondStereo,
        reverse: bool,
        stereo_atoms: Option<[usize; 2]>,
    ) -> TopologyBlock {
        let atoms = (0..6)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let (begin, end) = if reverse { (2, 1) } else { (1, 2) };
        let mut center =
            BondSpec::new(AtomId::new(begin), AtomId::new(end), order).with_stereo(stereo);
        if let Some([first, second]) = stereo_atoms {
            center = center.with_stereo_atoms(AtomId::new(first), AtomId::new(second));
        }
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(BondId::new(1), center),
            Bond::from_spec(
                BondId::new(2),
                BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
            ),
        ];
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed double-bond topology is valid")
    }

    fn final_chiral_double_graph() -> McsSeedTopology {
        McsSeedTopology {
            source_atoms: vec![0, 1, 2, 3],
            bonds: vec![McsSeedTopologyBond {
                source_bond: 1,
                begin_seed_atom: 1,
                end_seed_atom: 2,
            }],
        }
    }

    pub(super) fn final_chiral_bonds_source_skip_branches_and_target_orientation() {
        let coordinates = CoordinateBlock::default();
        let z = final_chiral_double(BondOrder::Double, BondStereo::Z, false, Some([0, 3]));
        let reversed = final_chiral_double(BondOrder::Double, BondStereo::Z, true, Some([3, 0]));
        let single = final_chiral_double(BondOrder::Single, BondStereo::Z, false, Some([0, 3]));
        let any = final_chiral_double(BondOrder::Double, BondStereo::Any, false, Some([0, 3]));
        let absent = final_chiral_double(BondOrder::Double, BondStereo::Z, false, None);
        let z_target = q107_target(&z, &coordinates);
        let reverse_target = q107_target(&reversed, &coordinates);
        let single_target = q107_target(&single, &coordinates);
        let any_target = q107_target(&any, &coordinates);
        let absent_target = q107_target(&absent, &coordinates);
        let graph = final_chiral_double_graph();
        let map = [0, 1, 2, 3];
        let target_vertices = [0, 1, 2, 3, 4, 5];
        for query in [&single_target, &any_target, &absent_target] {
            assert_eq!(
                mcs_final_chirality_check(&map, &map, query, &graph, &z_target, &target_vertices),
                Ok(true)
            );
        }
        for target in [&single_target, &any_target, &absent_target] {
            assert_eq!(
                mcs_final_chirality_check(&map, &map, &z_target, &graph, target, &target_vertices),
                Ok(true)
            );
        }
        assert_eq!(
            mcs_final_chirality_check(
                &map,
                &map,
                &z_target,
                &graph,
                &reverse_target,
                &target_vertices,
            ),
            Ok(true),
            "reversed target bond compares opposite stereo neighbor positions"
        );
    }

    pub(super) fn final_chiral_bonds_stereo_label_and_neighbor_match_counts() {
        let coordinates = CoordinateBlock::default();
        let z = final_chiral_double(BondOrder::Double, BondStereo::Z, false, Some([0, 3]));
        let e = final_chiral_double(BondOrder::Double, BondStereo::E, false, Some([0, 3]));
        let z_target = q107_target(&z, &coordinates);
        let e_target = q107_target(&e, &coordinates);
        let graph = final_chiral_double_graph();
        let query_map = [0, 1, 2, 3];
        let target_vertices = [0, 1, 2, 3, 4, 5];
        for (target_map, matches) in [([4, 1, 2, 5], 0), ([4, 1, 2, 3], 1), ([0, 1, 2, 3], 2)] {
            assert_eq!(
                mcs_final_chirality_check(
                    &query_map,
                    &target_map,
                    &z_target,
                    &graph,
                    &z_target,
                    &target_vertices,
                ),
                Ok(matches != 1)
            );
            assert_eq!(
                mcs_final_chirality_check(
                    &query_map,
                    &target_map,
                    &z_target,
                    &graph,
                    &e_target,
                    &target_vertices,
                ),
                Ok(matches == 1)
            );
        }
    }

    pub(super) fn final_chiral_bonds_missing_query_map_key_and_typed_invariants() {
        let coordinates = CoordinateBlock::default();
        let z = final_chiral_double(BondOrder::Double, BondStereo::Z, false, Some([0, 3]));
        let e = final_chiral_double(BondOrder::Double, BondStereo::E, false, Some([0, 3]));
        let z_target = q107_target(&z, &coordinates);
        let e_target = q107_target(&e, &coordinates);
        let partial = McsSeedTopology {
            source_atoms: vec![1, 2],
            bonds: vec![McsSeedTopologyBond {
                source_bond: 1,
                begin_seed_atom: 0,
                end_seed_atom: 1,
            }],
        };
        let target_vertices = [1, 2];
        assert_eq!(
            mcs_final_chirality_check(
                &[0, 1],
                &[0, 1],
                &z_target,
                &partial,
                &z_target,
                &target_vertices,
            ),
            Ok(true),
            "missing stereo neighbor qMap keys insert row zero"
        );
        assert_eq!(
            mcs_final_chirality_check(
                &[0, 1],
                &[0, 1],
                &z_target,
                &partial,
                &e_target,
                &target_vertices,
            ),
            Ok(false)
        );
        let graph = final_chiral_double_graph();
        assert!(matches!(
            mcs_final_chirality_check(
                &[0, 1, 2, 3],
                &[0, 1, 4, 3],
                &z_target,
                &graph,
                &z_target,
                &[0, 1, 2, 3, 4, 5],
            ),
            Err(McsCandidateMatchError::MappedBondMissing {
                side: "target",
                begin: 1,
                end: 4
            })
        ));
        let invalid = McsSeedTopology {
            source_atoms: vec![0, 1, 2, 3],
            bonds: vec![McsSeedTopologyBond {
                source_bond: 99,
                begin_seed_atom: 1,
                end_seed_atom: 2,
            }],
        };
        assert!(matches!(
            mcs_final_chirality_check(
                &[0, 1, 2, 3],
                &[0, 1, 2, 3],
                &z_target,
                &invalid,
                &z_target,
                &[0, 1, 2, 3, 4, 5],
            ),
            Err(McsCandidateMatchError::State(McsError::BondOutOfRange {
                side: "query molecule",
                bond: 99
            }))
        ));
    }

    pub(super) fn final_mapping_option_gates_and_ring_error_precede_user_hook() {
        let topology = q111_cycle(3);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let target = q107_target(&topology, &coordinates);
        let graph = McsSeedTopology {
            source_atoms: vec![0, 1],
            bonds: vec![McsSeedTopologyBond {
                source_bond: 0,
                begin_seed_atom: 0,
                end_seed_atom: 1,
            }],
        };
        let called = std::cell::Cell::new(0);
        let mut user = || -> Result<bool, McsCandidateMatchError> {
            called.set(called.get() + 1);
            Ok(true)
        };
        assert_eq!(
            mcs_final_mapping_accept(
                &graph,
                &query,
                &target,
                &McsParameters::default(),
                &[(0, 0), (1, 1)],
                Some(&mut user),
            ),
            Ok(true)
        );
        assert_eq!(
            called.get(),
            1,
            "disabled built-ins do not require ring state"
        );
        let params = McsParameters {
            bond_compare_parameters: McsBondCompareParameters {
                match_fused_rings: true,
                ..McsBondCompareParameters::default()
            },
            ..McsParameters::default()
        };
        assert_eq!(
            mcs_final_mapping_accept(
                &graph,
                &query,
                &target,
                &params,
                &[(0, 0), (1, 1)],
                Some(&mut user),
            ),
            Err(McsCandidateMatchError::State(McsError::MissingRingInfo {
                side: "query"
            }))
        );
        assert_eq!(called.get(), 1, "ring error must precede the user hook");
    }

    pub(super) fn final_mapping_mapped_chirality_rejects_before_user_hook() {
        let coordinates = CoordinateBlock::default();
        let query_topology = final_chiral_star(4, ChiralTag::TetrahedralCw);
        let target_topology = final_chiral_star(4, ChiralTag::Unspecified);
        let query_rings = cosmolkit_core::fast_find_rings(&query_topology).unwrap();
        let target_rings = cosmolkit_core::fast_find_rings(&target_topology).unwrap();
        let query = SearchTarget::new(
            &query_topology,
            &coordinates,
            &query_topology.stereo_groups,
            Some(&query_rings),
            None,
        );
        let target = SearchTarget::new(
            &target_topology,
            &coordinates,
            &target_topology.stereo_groups,
            Some(&target_rings),
            None,
        );
        let params = McsParameters {
            atom_compare_parameters: McsAtomCompareParameters {
                match_chiral_tag: true,
                ..McsAtomCompareParameters::default()
            },
            bond_compare_parameters: McsBondCompareParameters {
                match_fused_rings_strict: true,
                ..McsBondCompareParameters::default()
            },
            ..McsParameters::default()
        };
        let called = std::cell::Cell::new(0);
        let mut user = || -> Result<bool, McsCandidateMatchError> {
            called.set(called.get() + 1);
            Ok(true)
        };
        assert_eq!(
            mcs_final_mapping_accept(
                &final_chiral_star_graph(4),
                &query,
                &target,
                &params,
                &[(0, 0), (1, 1), (2, 2), (3, 3)],
                Some(&mut user),
            ),
            Ok(false)
        );
        assert_eq!(
            called.get(),
            0,
            "chirality rejection must precede the user hook"
        );
        let option_off = McsParameters::default();
        assert_eq!(
            mcs_final_mapping_accept(
                &final_chiral_star_graph(4),
                &query,
                &target,
                &option_off,
                &[(0, 0), (1, 1), (2, 2), (3, 3)],
                Some(&mut user),
            ),
            Ok(true)
        );
        assert_eq!(called.get(), 1);
    }

    pub(super) fn final_mapping_initial_seed_suppression_and_growth_restore_checks() {
        let coordinates = CoordinateBlock::default();
        let query_topology = final_chiral_star(4, ChiralTag::TetrahedralCw);
        let target_topology = final_chiral_star(4, ChiralTag::Unspecified);
        let query = q107_target(&query_topology, &coordinates);
        let target = q107_target(&target_topology, &coordinates);
        let tables = q126_tables(&query, &target, true, true);
        let params = McsParameters {
            atom_compare_parameters: McsAtomCompareParameters {
                match_chiral_tag: true,
                ..McsAtomCompareParameters::default()
            },
            ..McsParameters::default()
        };
        let initial = mcs_make_initial_seeds(
            &query,
            &[target],
            &[tables.clone()],
            1,
            &params,
            false,
            None,
            None,
        )
        .unwrap();
        assert_eq!(initial.queue.seeds.len(), 3);
        let restored = mcs_make_initial_seeds(
            &query,
            &[target],
            &[tables.clone()],
            1,
            &params,
            true,
            None,
            None,
        )
        .unwrap();
        assert!(restored.queue.seeds.is_empty());
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        seed.add_atom(0);
        seed.add_atom(1);
        seed.add_bond(&query, 0).unwrap();
        assert_eq!(
            mcs_match_full_candidate(
                &mut seed,
                &query,
                &[target],
                &[tables],
                1,
                Some(&params),
                None,
            ),
            Ok(false),
            "growth runs the built-in check without a user callback"
        );
    }

    #[test]
    fn q111_ring_pair_allows_size_difference_on_either_fused_side() {
        let fused = q111_fused_triangles();
        let square = q111_cycle(4);
        let fused_rings = cosmolkit_core::fast_find_rings(&fused).unwrap();
        let square_rings = cosmolkit_core::fast_find_rings(&square).unwrap();
        assert_eq!(fused_rings.num_bond_rings(BondId::new(0)), 2);
        let coordinates = CoordinateBlock::default();
        let fused = SearchTarget::new(
            &fused,
            &coordinates,
            &fused.stereo_groups,
            Some(&fused_rings),
            None,
        );
        let square = SearchTarget::new(
            &square,
            &coordinates,
            &square.stereo_groups,
            Some(&square_rings),
            None,
        );
        let params = McsBondCompareParameters {
            complete_rings_only: true,
            ..McsBondCompareParameters::default()
        };

        assert!(have_pair_of_compatible_rings(&params, &fused, 0, &square, 0).unwrap());
        assert!(have_pair_of_compatible_rings(&params, &square, 0, &fused, 0).unwrap());
    }

    #[test]
    fn q111_ring_pair_requires_detached_ring_information() {
        let triangle = q111_cycle(3);
        let coordinates = CoordinateBlock::default();
        let plain = q107_target(&triangle, &coordinates);
        assert_eq!(
            have_pair_of_compatible_rings(
                &McsBondCompareParameters::default(),
                &plain,
                0,
                &plain,
                0,
            ),
            Err(McsError::MissingRingInfo { side: "left" })
        );
    }

    #[test]
    fn q112_ring_match_compares_membership_for_both_ring_flag_values() {
        let triangle = q111_cycle(3);
        let chain = q110_topology(BondOrder::Single, BondStereo::None);
        let triangle_rings = cosmolkit_core::fast_find_rings(&triangle).unwrap();
        let chain_rings = cosmolkit_core::fast_find_rings(&chain).unwrap();
        let coordinates = CoordinateBlock::default();
        let triangle = SearchTarget::new(
            &triangle,
            &coordinates,
            &triangle.stereo_groups,
            Some(&triangle_rings),
            None,
        );
        let chain = SearchTarget::new(
            &chain,
            &coordinates,
            &chain.stereo_groups,
            Some(&chain_rings),
            None,
        );

        for ring_matches_ring_only in [false, true] {
            let params = McsBondCompareParameters {
                ring_matches_ring_only,
                ..McsBondCompareParameters::default()
            };
            assert!(check_bond_ring_match(&params, &triangle, 0, &triangle, 0).unwrap());
            assert!(check_bond_ring_match(&params, &chain, 0, &chain, 0).unwrap());
            assert!(!check_bond_ring_match(&params, &triangle, 0, &chain, 0).unwrap());
        }
    }

    #[test]
    fn q112_complete_ring_flag_requires_a_compatible_ring_pair() {
        let triangle = q111_cycle(3);
        let square = q111_cycle(4);
        let fused = q111_fused_triangles();
        let triangle_rings = cosmolkit_core::fast_find_rings(&triangle).unwrap();
        let square_rings = cosmolkit_core::fast_find_rings(&square).unwrap();
        let fused_rings = cosmolkit_core::fast_find_rings(&fused).unwrap();
        let coordinates = CoordinateBlock::default();
        let triangle = SearchTarget::new(
            &triangle,
            &coordinates,
            &triangle.stereo_groups,
            Some(&triangle_rings),
            None,
        );
        let square = SearchTarget::new(
            &square,
            &coordinates,
            &square.stereo_groups,
            Some(&square_rings),
            None,
        );
        let fused = SearchTarget::new(
            &fused,
            &coordinates,
            &fused.stereo_groups,
            Some(&fused_rings),
            None,
        );

        assert!(
            check_bond_ring_match(
                &McsBondCompareParameters::default(),
                &triangle,
                0,
                &square,
                0,
            )
            .unwrap()
        );
        let complete = McsBondCompareParameters {
            complete_rings_only: true,
            ..McsBondCompareParameters::default()
        };
        assert!(!check_bond_ring_match(&complete, &triangle, 0, &square, 0).unwrap());
        assert!(check_bond_ring_match(&complete, &fused, 0, &square, 0).unwrap());
    }

    #[test]
    fn q108_any_bond_skips_disabled_prechecks_and_complete_only_flag() {
        let topology = q110_topology(BondOrder::Single, BondStereo::None);
        let coordinates = CoordinateBlock::default();
        let target = q107_target(&topology, &coordinates);

        assert!(
            mcs_bond_compare_any(&McsBondCompareParameters::default(), &target, 7, &target, 9,)
                .unwrap()
        );
        assert!(
            mcs_bond_compare_any(
                &McsBondCompareParameters {
                    complete_rings_only: true,
                    ..McsBondCompareParameters::default()
                },
                &target,
                7,
                &target,
                9,
            )
            .unwrap()
        );
    }

    #[test]
    fn q108_any_bond_applies_ring_membership_only_under_source_gate() {
        let triangle = q111_cycle(3);
        let chain = q110_topology(BondOrder::Single, BondStereo::None);
        let triangle_rings = cosmolkit_core::fast_find_rings(&triangle).unwrap();
        let chain_rings = cosmolkit_core::fast_find_rings(&chain).unwrap();
        let coordinates = CoordinateBlock::default();
        let triangle = SearchTarget::new(
            &triangle,
            &coordinates,
            &triangle.stereo_groups,
            Some(&triangle_rings),
            None,
        );
        let chain = SearchTarget::new(
            &chain,
            &coordinates,
            &chain.stereo_groups,
            Some(&chain_rings),
            None,
        );
        let params = McsBondCompareParameters {
            ring_matches_ring_only: true,
            ..McsBondCompareParameters::default()
        };

        assert!(!mcs_bond_compare_any(&params, &triangle, 0, &chain, 0).unwrap());
        assert!(mcs_bond_compare_any(&params, &triangle, 0, &triangle, 0).unwrap());
    }

    #[test]
    fn q108_any_bond_checks_stereo_before_ring_state() {
        let z = q110_topology(BondOrder::Double, BondStereo::Z);
        let e = q110_topology(BondOrder::Double, BondStereo::E);
        let coordinates = CoordinateBlock::default();
        let z = q107_target(&z, &coordinates);
        let e = q107_target(&e, &coordinates);
        let params = McsBondCompareParameters {
            ring_matches_ring_only: true,
            match_stereo: true,
            ..McsBondCompareParameters::default()
        };

        assert!(!mcs_bond_compare_any(&params, &z, 0, &e, 0).unwrap());
        assert_eq!(
            mcs_bond_compare_any(&params, &z, 0, &z, 0),
            Err(McsError::MissingRingInfo { side: "left" })
        );
        assert!(mcs_bond_compare_any(&McsBondCompareParameters::default(), &z, 0, &e, 0,).unwrap());
    }

    #[test]
    fn q109_order_matrix_preserves_all_source_equivalences() {
        let all_orders = [
            BondOrder::Unspecified,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Triple,
            BondOrder::Quadruple,
            BondOrder::Quintuple,
            BondOrder::Hextuple,
            BondOrder::OneAndHalf,
            BondOrder::TwoAndHalf,
            BondOrder::ThreeAndHalf,
            BondOrder::FourAndHalf,
            BondOrder::FiveAndHalf,
            BondOrder::Aromatic,
            BondOrder::Ionic,
            BondOrder::Hydrogen,
            BondOrder::ThreeCenter,
            BondOrder::DativeOne,
            BondOrder::Dative,
            BondOrder::DativeLeft,
            BondOrder::DativeRight,
            BondOrder::Other,
            BondOrder::Zero,
        ];
        for order in all_orders {
            assert!(bond_orders_match(order, order, false));
            assert!(bond_orders_match(BondOrder::Unspecified, order, false));
            assert!(bond_orders_match(order, BondOrder::Unspecified, false));
            assert!(bond_orders_match(BondOrder::Zero, order, false));
            assert!(bond_orders_match(order, BondOrder::Zero, false));
        }

        let relaxed_pairs = [
            (BondOrder::Single, BondOrder::Aromatic),
            (BondOrder::Single, BondOrder::OneAndHalf),
            (BondOrder::Double, BondOrder::TwoAndHalf),
            (BondOrder::Triple, BondOrder::ThreeAndHalf),
            (BondOrder::Quadruple, BondOrder::FourAndHalf),
            (BondOrder::Quintuple, BondOrder::FiveAndHalf),
        ];
        for (left, right) in relaxed_pairs {
            assert!(bond_orders_match(left, right, true));
            assert!(bond_orders_match(right, left, true));
            assert!(!bond_orders_match(left, right, false));
            assert!(!bond_orders_match(right, left, false));
        }
        assert!(!bond_orders_match(
            BondOrder::Double,
            BondOrder::Aromatic,
            true,
        ));
    }

    #[test]
    fn q109_order_and_exact_comparators_distinguish_aromatic_relaxation() {
        let single = q110_topology(BondOrder::Single, BondStereo::None);
        let aromatic = q110_topology(BondOrder::Aromatic, BondStereo::None);
        let double = q110_topology(BondOrder::Double, BondStereo::None);
        let coordinates = CoordinateBlock::default();
        let single = q107_target(&single, &coordinates);
        let aromatic = q107_target(&aromatic, &coordinates);
        let double = q107_target(&double, &coordinates);
        let params = McsBondCompareParameters::default();

        assert!(mcs_bond_compare_order(&params, &single, 0, &aromatic, 0).unwrap());
        assert!(!mcs_bond_compare_order_exact(&params, &single, 0, &aromatic, 0).unwrap());
        assert!(!mcs_bond_compare_order(&params, &double, 0, &aromatic, 0).unwrap());
        assert!(mcs_bond_compare_order_exact(&params, &double, 0, &double, 0).unwrap());
    }

    #[test]
    fn q109_order_comparators_apply_prechecks_after_order_match() {
        let single = q110_topology(BondOrder::Single, BondStereo::None);
        let double_z = q110_topology(BondOrder::Double, BondStereo::Z);
        let double_e = q110_topology(BondOrder::Double, BondStereo::E);
        let coordinates = CoordinateBlock::default();
        let single = q107_target(&single, &coordinates);
        let double_z = q107_target(&double_z, &coordinates);
        let double_e = q107_target(&double_e, &coordinates);
        let ring_params = McsBondCompareParameters {
            ring_matches_ring_only: true,
            ..McsBondCompareParameters::default()
        };

        assert!(!mcs_bond_compare_order(&ring_params, &single, 0, &double_z, 0).unwrap());
        assert_eq!(
            mcs_bond_compare_order(&ring_params, &single, 0, &single, 0),
            Err(McsError::MissingRingInfo { side: "left" })
        );
        let stereo_params = McsBondCompareParameters {
            match_stereo: true,
            ..McsBondCompareParameters::default()
        };
        assert!(!mcs_bond_compare_order(&stereo_params, &double_z, 0, &double_e, 0,).unwrap());
        assert!(
            !mcs_bond_compare_order_exact(&stereo_params, &double_z, 0, &double_e, 0,).unwrap()
        );
    }

    #[test]
    fn q113_input_order_rejects_empty_single_and_above_one_threshold() {
        assert_eq!(
            prepare_mcs_input_order(&[], 1.0),
            Err(McsError::TooFewInputs { count: 0 })
        );

        let topology = q113_chain(1);
        let coordinates = CoordinateBlock::default();
        let target = q107_target(&topology, &coordinates);
        assert_eq!(
            prepare_mcs_input_order(&[target], 1.0),
            Err(McsError::TooFewInputs { count: 1 })
        );

        let topologies = [q113_chain(1), q113_chain(2)];
        let targets = topologies
            .iter()
            .map(|topology| q107_target(topology, &coordinates))
            .collect::<Vec<_>>();
        assert_eq!(
            prepare_mcs_input_order(&targets, 1.000_001),
            Err(McsError::ThresholdAboveOne)
        );
    }

    #[test]
    fn q113_threshold_count_preserves_source_clamping_and_rounding() {
        let topologies = [q113_chain(1), q113_chain(2), q113_chain(3), q113_chain(4)];
        let coordinates = CoordinateBlock::default();
        let targets = topologies
            .iter()
            .map(|topology| q107_target(topology, &coordinates))
            .collect::<Vec<_>>();

        assert_eq!(
            prepare_mcs_input_order(&targets, 1.0)
                .unwrap()
                .threshold_count,
            3
        );
        assert_eq!(
            prepare_mcs_input_order(&targets, 0.5)
                .unwrap()
                .threshold_count,
            1
        );
        assert_eq!(
            prepare_mcs_input_order(&targets, 0.51)
                .unwrap()
                .threshold_count,
            2
        );
        assert_eq!(
            prepare_mcs_input_order(&targets, 0.0)
                .unwrap()
                .threshold_count,
            1
        );
        assert_eq!(
            prepare_mcs_input_order(&targets, -4.0)
                .unwrap()
                .threshold_count,
            1
        );
    }

    #[test]
    fn q113_input_order_is_stable_and_skips_only_leading_empty_queries() {
        let topologies = [
            q113_chain(3),
            q113_chain(0),
            q113_chain(2),
            q113_chain(1),
            q113_chain(2),
        ];
        let coordinates = CoordinateBlock::default();
        let targets = topologies
            .iter()
            .map(|topology| q107_target(topology, &coordinates))
            .collect::<Vec<_>>();
        let order = prepare_mcs_input_order(&targets, 0.4).unwrap();

        assert_eq!(order.molecule_indices, vec![1, 3, 2, 4, 0]);
        assert_eq!(order.threshold_count, 1);
        assert_eq!(order.end_index, 4);
        assert_eq!(order.start_index, 1);
        assert_eq!(order.molecule_indices[order.start_index], 3);
    }

    #[test]
    fn q114_match_tables_preserve_source_dimensions_and_row_column_lookup() {
        let topology = q114_topology();
        let coordinates = CoordinateBlock::default();
        let target = q107_target(&topology, &coordinates);
        let tables = build_initial_match_tables(&McsParameters::default(), &target).unwrap();

        assert_eq!((tables.atoms.rows, tables.atoms.columns), (3, 3));
        assert_eq!((tables.bonds.rows, tables.bonds.columns), (2, 2));
        assert_eq!(tables.atoms.get(0, 0), Some(true));
        assert_eq!(tables.atoms.get(0, 1), Some(false));
        assert_eq!(tables.atoms.get(3, 0), None);
        assert_eq!(tables.atoms.get(0, 3), None);
        assert_eq!(tables.bonds.get(0, 1), Some(true));
        assert_eq!(tables.bonds.get(1, 0), Some(true));
    }

    #[test]
    fn q114_match_tables_dispatch_every_completed_comparator() {
        let topology = q114_topology();
        let coordinates = CoordinateBlock::default();
        let target = q107_target(&topology, &coordinates);
        let mut params = McsParameters::default();

        params.atom_comparator = AtomComparator::AtomCompareAny;
        params.bond_comparator = BondComparator::BondCompareAny;
        let any = build_initial_match_tables(&params, &target).unwrap();
        assert_eq!(any.atoms.get(0, 1), Some(true));
        assert_eq!(any.bonds.get(0, 1), Some(true));

        params.atom_comparator = AtomComparator::AtomCompareAnyHeavyAtom;
        params.bond_comparator = BondComparator::BondCompareOrderExact;
        let heavy_exact = build_initial_match_tables(&params, &target).unwrap();
        assert_eq!(heavy_exact.atoms.get(0, 1), Some(true));
        assert_eq!(heavy_exact.atoms.get(0, 2), Some(false));
        assert_eq!(heavy_exact.bonds.get(0, 1), Some(false));

        params.atom_comparator = AtomComparator::AtomCompareIsotopes;
        params.bond_comparator = BondComparator::BondCompareOrder;
        let isotope_order = build_initial_match_tables(&params, &target).unwrap();
        assert_eq!(isotope_order.atoms.get(0, 1), Some(true));
        assert_eq!(isotope_order.atoms.get(1, 2), Some(false));
        assert_eq!(isotope_order.bonds.get(0, 1), Some(true));
    }

    #[test]
    fn q114_match_tables_propagate_comparator_state_errors() {
        let topology = q114_topology();
        let coordinates = CoordinateBlock::default();
        let target = q107_target(&topology, &coordinates);
        let params = McsParameters {
            atom_compare_parameters: McsAtomCompareParameters {
                match_valences: true,
                ..McsAtomCompareParameters::default()
            },
            ..McsParameters::default()
        };

        assert_eq!(
            build_initial_match_tables(&params, &target),
            Err(McsError::MissingValence { side: "left" })
        );
    }

    fn q134_table_topology(
        specs: &[AtomSpec],
        edges: &[(usize, usize, BondOrder)],
    ) -> TopologyBlock {
        let atoms = specs
            .iter()
            .cloned()
            .enumerate()
            .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
            .collect();
        let bonds = edges
            .iter()
            .copied()
            .enumerate()
            .map(|(index, (begin, end, order))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }

    #[test]
    fn q134_rectangular_tables_preserve_query_rows_target_columns_and_comparators() {
        let query_topology = q134_table_topology(
            &[
                AtomSpec::new(Element::C).with_isotope(13),
                AtomSpec::new(Element::O).with_isotope(13),
                AtomSpec::new(Element::H).with_isotope(1),
            ],
            &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Double)],
        );
        let target_topology = q134_table_topology(
            &[
                AtomSpec::new(Element::O).with_isotope(13),
                AtomSpec::new(Element::C).with_isotope(13),
                AtomSpec::new(Element::H).with_isotope(1),
                AtomSpec::new(Element::C).with_isotope(12),
            ],
            &[
                (0, 1, BondOrder::Double),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Aromatic),
            ],
        );
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&query_topology, &coordinates);
        let target = q107_target(&target_topology, &coordinates);
        let mut params = McsParameters::default();
        let query_snapshot = query_topology.clone();
        let target_snapshot = target_topology.clone();

        let forward = build_query_target_match_tables(&params, &query, &target).unwrap();
        assert_eq!((forward.atoms.rows, forward.atoms.columns), (3, 4));
        assert_eq!((forward.bonds.rows, forward.bonds.columns), (2, 3));
        assert_eq!(
            forward.atoms.values,
            vec![
                false, true, false, true, true, false, false, false, false, false, true, false
            ]
        );
        assert_eq!(
            forward.bonds.values,
            vec![false, true, true, true, false, false]
        );
        assert_eq!(forward.atoms.get(0, 4), None);
        assert_eq!(forward.bonds.get(2, 0), None);

        let reverse = build_query_target_match_tables(&params, &target, &query).unwrap();
        assert_eq!((reverse.atoms.rows, reverse.atoms.columns), (4, 3));
        assert_eq!((reverse.bonds.rows, reverse.bonds.columns), (3, 2));
        assert_eq!(
            reverse.atoms.values,
            vec![
                false, true, false, true, false, false, false, false, true, true, false, false
            ]
        );
        assert_eq!(
            reverse.bonds.values,
            vec![false, true, true, false, true, false]
        );

        params.atom_comparator = AtomComparator::AtomCompareAny;
        params.bond_comparator = BondComparator::BondCompareAny;
        let any = build_query_target_match_tables(&params, &query, &target).unwrap();
        assert!(any.atoms.values.iter().all(|cell| *cell));
        assert!(any.bonds.values.iter().all(|cell| *cell));

        params.atom_comparator = AtomComparator::AtomCompareAnyHeavyAtom;
        params.bond_comparator = BondComparator::BondCompareOrderExact;
        let heavy_exact = build_query_target_match_tables(&params, &query, &target).unwrap();
        assert_eq!(
            heavy_exact.atoms.values,
            vec![
                true, true, false, true, true, true, false, true, false, false, true, false
            ]
        );
        assert_eq!(
            heavy_exact.bonds.values,
            vec![false, true, false, true, false, false]
        );

        params.atom_comparator = AtomComparator::AtomCompareIsotopes;
        let isotopes = build_query_target_match_tables(&params, &query, &target).unwrap();
        assert_eq!(
            isotopes.atoms.values,
            vec![
                true, true, false, false, true, true, false, false, false, false, true, false
            ]
        );
        assert_eq!(
            build_initial_match_tables(&params, &query).unwrap(),
            build_query_target_match_tables(&params, &query, &query).unwrap()
        );
        assert_eq!(query_topology, query_snapshot);
        assert_eq!(target_topology, target_snapshot);
    }

    #[test]
    fn q134_rectangular_tables_keep_zero_dimensions_and_first_typed_error() {
        let empty = q134_table_topology(&[], &[]);
        let carbon = q134_table_topology(&[AtomSpec::new(Element::C)], &[]);
        let coordinates = CoordinateBlock::default();
        let empty_target = q107_target(&empty, &coordinates);
        let carbon_target = q107_target(&carbon, &coordinates);
        let params = McsParameters::default();
        let empty_to_carbon =
            build_query_target_match_tables(&params, &empty_target, &carbon_target).unwrap();
        assert_eq!(
            (empty_to_carbon.atoms.rows, empty_to_carbon.atoms.columns),
            (0, 1)
        );
        assert_eq!(
            (empty_to_carbon.bonds.rows, empty_to_carbon.bonds.columns),
            (0, 0)
        );
        assert!(empty_to_carbon.atoms.values.is_empty());
        let carbon_to_empty =
            build_query_target_match_tables(&params, &carbon_target, &empty_target).unwrap();
        assert_eq!(
            (carbon_to_empty.atoms.rows, carbon_to_empty.atoms.columns),
            (1, 0)
        );
        assert_eq!(
            (carbon_to_empty.bonds.rows, carbon_to_empty.bonds.columns),
            (0, 0)
        );
        assert!(carbon_to_empty.atoms.values.is_empty());

        let valence = ValenceAssignment {
            explicit_valence: vec![4],
            implicit_hydrogens: vec![0],
        };
        let with_valence = SearchTarget::new(
            &carbon,
            &coordinates,
            &carbon.stereo_groups,
            None,
            Some(&valence),
        );
        let mut params = McsParameters::default();
        params.atom_compare_parameters.match_valences = true;
        assert_eq!(
            build_query_target_match_tables(&params, &carbon_target, &with_valence),
            Err(McsError::MissingValence { side: "left" })
        );
        assert_eq!(
            build_query_target_match_tables(&params, &with_valence, &carbon_target),
            Err(McsError::MissingValence { side: "right" })
        );
    }

    #[test]
    fn q134_table_comparison_retains_source_directional_chirality_gate() {
        let chiral = q107_topology(&[ChiralTag::TetrahedralCw]);
        let plain = q107_topology(&[ChiralTag::Unspecified]);
        let coordinates = CoordinateBlock::default();
        let chiral_target = q107_target(&chiral, &coordinates);
        let plain_target = q107_target(&plain, &coordinates);
        let mut params = McsParameters::default();
        params.atom_compare_parameters.match_chiral_tag = true;
        assert_eq!(
            build_query_target_match_tables(&params, &chiral_target, &plain_target)
                .unwrap()
                .atoms
                .get(0, 0),
            Some(false)
        );
        assert_eq!(
            build_query_target_match_tables(&params, &plain_target, &chiral_target)
                .unwrap()
                .atoms
                .get(0, 0),
            Some(true)
        );
    }

    #[test]
    fn q114a_seed_and_fragment_defaults_match_source_empty_state() {
        let fragment = McsMoleculeFragment::default();
        assert!(fragment.atoms.is_empty());
        assert!(fragment.bonds.is_empty());
        assert!(fragment.seed_atom_index_map.is_empty());

        let seed = McsSeed::default();
        assert!(seed.new_bonds.is_empty());
        assert!(!seed.store_all_degenerate_mcs);
        assert!(!seed.copy_complete);
        assert_eq!(seed.growing_stage, 0);
        assert_eq!(seed.molecule_fragment, fragment);
        assert!(seed.topology.source_atoms.is_empty());
        assert!(seed.topology.bonds.is_empty());
        assert!(seed.excluded_bonds.is_empty());
        assert_eq!(seed.last_added_atoms_begin_index, 0);
        assert_eq!(seed.last_added_bonds_begin_index, 0);
        assert_eq!(seed.remaining_bonds, 0);
        assert_eq!(seed.remaining_atoms, 0);
        assert!(seed.match_result.is_empty());
    }

    #[test]
    fn q114a_seed_clone_owns_state_without_fabricating_completion() {
        let mut seed = McsSeed::default();
        seed.molecule_fragment.atoms.push(4);
        seed.molecule_fragment.seed_atom_index_map.insert(4, 0);
        seed.topology.source_atoms.push(4);
        seed.excluded_bonds = vec![false, true];
        seed.match_result.push(McsTargetMatch {
            empty: false,
            matched_atom_size: 1,
            matched_bond_size: 0,
            target_atom_indices: vec![7],
            target_bond_indices: Vec::new(),
            visited_target_bonds: vec![false],
            visited_target_atoms: vec![false, true],
        });

        let mut copied = seed.clone();
        copied.molecule_fragment.atoms.push(8);
        copied.molecule_fragment.seed_atom_index_map.insert(8, 1);
        copied.topology.source_atoms.push(8);
        copied.excluded_bonds[0] = true;
        copied.match_result[0].target_atom_indices[0] = 9;

        assert_eq!(seed.molecule_fragment.atoms, vec![4]);
        assert_eq!(seed.molecule_fragment.seed_atom_index_map.get(&8), None);
        assert_eq!(seed.topology.source_atoms, vec![4]);
        assert_eq!(seed.excluded_bonds, vec![false, true]);
        assert_eq!(seed.match_result[0].target_atom_indices, vec![7]);
        assert!(!copied.copy_complete);
        assert_eq!(copied.growing_stage, 0);
    }

    #[test]
    fn q115_add_atom_preserves_append_order_and_source_duplicate_mapping() {
        let mut seed = McsSeed::default();

        assert_eq!(seed.add_atom(4), 0);
        assert_eq!(seed.add_atom(2), 1);
        assert_eq!(seed.add_atom(4), 2);

        assert_eq!(seed.molecule_fragment.atoms, vec![4, 2, 4]);
        assert_eq!(seed.topology.source_atoms, vec![4, 2, 4]);
        assert_eq!(seed.molecule_fragment.seed_atom_index_map.get(&2), Some(&1));
        assert_eq!(seed.molecule_fragment.seed_atom_index_map.get(&4), Some(&2));
        assert!(seed.molecule_fragment.bonds.is_empty());
        assert!(seed.topology.bonds.is_empty());
    }

    #[test]
    fn q116_add_bond_maps_source_endpoints_and_rejects_duplicates() {
        let topology = q114_topology();
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        assert_eq!(seed.add_atom(1), 0);
        assert_eq!(seed.add_atom(0), 1);

        assert_eq!(seed.add_bond(&query, 0), Ok(1));
        assert_eq!(seed.molecule_fragment.bonds, vec![0]);
        assert_eq!(seed.excluded_bonds, vec![true, false]);
        assert_eq!(
            seed.topology.bonds,
            vec![McsSeedTopologyBond {
                source_bond: 0,
                begin_seed_atom: 1,
                end_seed_atom: 0,
            }]
        );

        let before_duplicate = seed.clone();
        assert_eq!(
            seed.add_bond(&query, 0),
            Err(McsError::SeedBondAlreadyExcluded { bond: 0 })
        );
        assert_eq!(seed, before_duplicate);
    }

    #[test]
    fn q116_add_bond_validates_source_and_exclusion_rows_before_mutation() {
        let topology = q114_topology();
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };

        let before_missing_source = seed.clone();
        assert_eq!(
            seed.add_bond(&query, 2),
            Err(McsError::BondOutOfRange {
                side: "query",
                bond: 2,
            })
        );
        assert_eq!(seed, before_missing_source);

        seed.excluded_bonds.truncate(1);
        let before_short_mask = seed.clone();
        assert_eq!(
            seed.add_bond(&query, 1),
            Err(McsError::SeedExcludedBondOutOfRange { bond: 1, count: 1 })
        );
        assert_eq!(seed, before_short_mask);
    }

    #[test]
    fn q116_add_bond_preserves_source_mutation_order_on_missing_endpoints() {
        let topology = q114_topology();
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);

        let mut missing_begin = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        missing_begin.add_atom(1);
        assert_eq!(
            missing_begin.add_bond(&query, 0),
            Err(McsError::SeedAtomMappingMissing { atom: 0 })
        );
        assert_eq!(missing_begin.excluded_bonds, vec![true, false]);
        assert_eq!(missing_begin.molecule_fragment.bonds, vec![0]);
        assert!(missing_begin.topology.bonds.is_empty());

        let mut missing_end = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        missing_end.add_atom(0);
        assert_eq!(
            missing_end.add_bond(&query, 0),
            Err(McsError::SeedAtomMappingMissing { atom: 1 })
        );
        assert_eq!(missing_end.excluded_bonds, vec![true, false]);
        assert_eq!(missing_end.molecule_fragment.bonds, vec![0]);
        assert!(missing_end.topology.bonds.is_empty());
    }

    #[test]
    fn q116a_create_from_parent_copies_source_fields_and_resets_growth_boundary() {
        let topology = q114_topology();
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let mut parent = McsSeed {
            store_all_degenerate_mcs: true,
            growing_stage: 9,
            excluded_bonds: vec![false; query.num_bonds()],
            remaining_bonds: 7,
            remaining_atoms: 5,
            ..McsSeed::default()
        };
        parent.add_atom(0);
        parent.add_atom(1);
        parent.add_bond(&query, 0).unwrap();

        let retained_frontier = McsNewBond {
            bond_index: 1,
            new_atom_index: 2,
            end_atom_index: None,
            new_atom: Some(2),
        };
        let mut child = McsSeed {
            new_bonds: vec![retained_frontier.clone()],
            copy_complete: true,
            growing_stage: 4,
            ..McsSeed::default()
        };
        child.create_from_parent(&parent);

        assert_eq!(child.molecule_fragment, parent.molecule_fragment);
        assert_eq!(child.topology, parent.topology);
        assert_eq!(child.excluded_bonds, parent.excluded_bonds);
        assert_eq!(child.remaining_bonds, 7);
        assert_eq!(child.remaining_atoms, 5);
        assert!(child.store_all_degenerate_mcs);
        assert_eq!(child.last_added_atoms_begin_index, 2);
        assert_eq!(child.last_added_bonds_begin_index, 1);
        assert_eq!(child.growing_stage, 0);
        assert_eq!(child.new_bonds, vec![retained_frontier]);
        assert!(child.copy_complete);
    }

    #[test]
    fn q116a_create_from_parent_owns_topology_maps_and_exclusions() {
        let topology = q114_topology();
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let mut parent = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        parent.add_atom(0);
        parent.add_atom(1);
        parent.add_bond(&query, 0).unwrap();

        let mut child = McsSeed::default();
        child.create_from_parent(&parent);
        child.molecule_fragment.atoms[0] = 2;
        child.molecule_fragment.seed_atom_index_map.insert(2, 0);
        child.topology.source_atoms[0] = 2;
        child.excluded_bonds[1] = true;

        assert_eq!(parent.molecule_fragment.atoms, vec![0, 1]);
        assert_eq!(parent.molecule_fragment.seed_atom_index_map.get(&2), None);
        assert_eq!(parent.topology.source_atoms, vec![0, 1]);
        assert_eq!(parent.excluded_bonds, vec![true, false]);
    }

    #[test]
    fn q116b_frontier_distinguishes_new_and_existing_endpoints_in_source_order() {
        let topology = q114_topology();
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let mut seed = McsSeed::default();
        seed.add_atom(0);

        seed.add_new_bond_from_atom(&query, 0, 0).unwrap();
        assert_eq!(
            seed.new_bonds,
            vec![McsNewBond {
                bond_index: 0,
                new_atom_index: 1,
                end_atom_index: None,
                new_atom: Some(1),
            }]
        );

        seed.add_atom(1);
        seed.add_atom(1);
        seed.add_new_bond_from_atom(&query, 0, 0).unwrap();
        seed.add_new_bond_from_atom(&query, 1, 0).unwrap();
        assert_eq!(
            &seed.new_bonds[1..],
            [
                McsNewBond {
                    bond_index: 0,
                    new_atom_index: 1,
                    end_atom_index: Some(1),
                    new_atom: None,
                },
                McsNewBond {
                    bond_index: 0,
                    new_atom_index: 0,
                    end_atom_index: Some(0),
                    new_atom: None,
                },
            ]
        );
    }

    #[test]
    fn q116b_frontier_rejects_invalid_source_identity_without_append() {
        let topology = q114_topology();
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let mut seed = McsSeed::default();

        assert_eq!(
            seed.add_new_bond_from_atom(&query, 2, 0),
            Err(McsError::SeedAtomNotBondEndpoint { atom: 2, bond: 0 })
        );
        assert_eq!(
            seed.add_new_bond_from_atom(&query, 3, 0),
            Err(McsError::AtomOutOfRange {
                side: "query",
                atom: 3,
            })
        );
        assert_eq!(
            seed.add_new_bond_from_atom(&query, 0, 2),
            Err(McsError::BondOutOfRange {
                side: "query",
                bond: 2,
            })
        );
        assert!(seed.new_bonds.is_empty());
    }

    #[test]
    fn q117_accumulation_adds_shared_new_atom_once_and_bonds_in_source_order() {
        let topology = q111_cycle(3);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let mut parent = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            remaining_bonds: 2,
            remaining_atoms: 1,
            ..McsSeed::default()
        };
        parent.add_atom(0);
        parent.add_atom(2);
        parent.add_bond(&query, 2).unwrap();
        parent.add_new_bond_from_atom(&query, 0, 0).unwrap();
        parent.add_new_bond_from_atom(&query, 2, 1).unwrap();

        let mut child = McsSeed::default();
        child.create_from_parent(&parent);
        let added = parent.add_new_bonds_to_seed(&query, &mut child).unwrap();

        assert_eq!(added, vec![false, true, false]);
        assert_eq!(child.molecule_fragment.atoms, vec![0, 2, 1]);
        assert_eq!(child.molecule_fragment.bonds, vec![2, 0, 1]);
        assert_eq!(
            child.topology.bonds,
            vec![
                McsSeedTopologyBond {
                    source_bond: 2,
                    begin_seed_atom: 1,
                    end_seed_atom: 0,
                },
                McsSeedTopologyBond {
                    source_bond: 0,
                    begin_seed_atom: 0,
                    end_seed_atom: 2,
                },
                McsSeedTopologyBond {
                    source_bond: 1,
                    begin_seed_atom: 2,
                    end_seed_atom: 1,
                },
            ]
        );
        assert_eq!(child.excluded_bonds, vec![true, true, true]);
        assert_eq!(child.remaining_bonds, 0);
        assert_eq!(child.remaining_atoms, 0);
    }

    #[test]
    fn q117_accumulation_preserves_existing_endpoint_and_missing_payload_constraints() {
        let topology = q114_topology();
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let mut parent = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            remaining_bonds: 1,
            ..McsSeed::default()
        };
        parent.add_atom(0);
        parent.add_atom(1);
        parent.add_new_bond_from_atom(&query, 0, 0).unwrap();
        let mut child = McsSeed::default();
        child.create_from_parent(&parent);

        assert_eq!(
            parent.add_new_bonds_to_seed(&query, &mut child).unwrap(),
            vec![false, false, false]
        );
        assert_eq!(child.molecule_fragment.atoms, vec![0, 1]);
        assert_eq!(child.molecule_fragment.bonds, vec![0]);

        let invalid = McsSeed {
            new_bonds: vec![McsNewBond {
                bond_index: 0,
                new_atom_index: 1,
                end_atom_index: None,
                new_atom: None,
            }],
            remaining_bonds: 1,
            remaining_atoms: 1,
            ..McsSeed::default()
        };
        let mut invalid_child = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        invalid_child.add_atom(0);
        assert_eq!(
            invalid.add_new_bonds_to_seed(&query, &mut invalid_child),
            Err(McsError::SeedNewAtomMissing { atom: 1 })
        );
        assert_eq!(invalid_child.molecule_fragment.atoms, vec![0]);
        assert!(invalid_child.molecule_fragment.bonds.is_empty());
    }

    #[test]
    fn q119_duplicate_key_sorts_source_ids_and_retains_multiplicity() {
        let mut left = McsDuplicateSeedKey::default();
        for atom in [5, 1, 5] {
            left.add_atom(atom);
        }
        for bond in [7, 2] {
            left.add_bond(bond);
        }
        assert_eq!(left.atom_indices, vec![1, 5, 5]);
        assert_eq!(left.bond_indices, vec![2, 7]);

        let mut reordered = McsDuplicateSeedKey::default();
        for atom in [5, 5, 1] {
            reordered.add_atom(atom);
        }
        for bond in [2, 7] {
            reordered.add_bond(bond);
        }
        assert!(left.source_equal(&reordered));

        let mut missing_duplicate = reordered.clone();
        missing_duplicate.atom_indices.remove(2);
        assert!(!left.source_equal(&missing_duplicate));
        let mut different_bond = reordered;
        different_bond.bond_indices[1] = 8;
        assert!(!left.source_equal(&different_bond));
    }

    #[test]
    fn q119_seed_key_uses_only_sorted_atom_and_bond_source_state() {
        let topology = q111_cycle(3);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let mut left = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        for atom in [2, 0, 1] {
            left.add_atom(atom);
        }
        left.add_bond(&query, 2).unwrap();
        left.add_bond(&query, 0).unwrap();

        let mut right = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            remaining_atoms: 99,
            ..McsSeed::default()
        };
        for atom in [1, 0, 2] {
            right.add_atom(atom);
        }
        right.add_bond(&query, 0).unwrap();
        right.add_bond(&query, 2).unwrap();

        assert_ne!(left.molecule_fragment.atoms, right.molecule_fragment.atoms);
        assert_ne!(left.molecule_fragment.bonds, right.molecule_fragment.bonds);
        assert!(left.duplicate_key.source_equal(&right.duplicate_key));

        let mut child = McsSeed::default();
        child.create_from_parent(&left);
        assert!(child.duplicate_key.source_equal(&left.duplicate_key));
        child.add_atom(2);
        assert!(!child.duplicate_key.source_equal(&left.duplicate_key));
    }

    #[test]
    fn q121_remaining_bound_counts_cycle_bonds_and_atoms_once() {
        let topology = q111_cycle(3);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        seed.add_atom(0);
        seed.add_atom(1);
        seed.add_bond(&query, 0).unwrap();
        seed.last_added_atoms_begin_index = 0;

        seed.compute_remaining_size(&query).unwrap();
        assert_eq!(seed.remaining_bonds, 2);
        assert_eq!(seed.remaining_atoms, 1);
    }

    #[test]
    fn q121_remaining_bound_starts_only_at_last_added_frontier() {
        let topology = TopologyBlock::try_from_parts(
            (0..5)
                .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
                .collect(),
            [(0, 1), (1, 2), (3, 4)]
                .into_iter()
                .enumerate()
                .map(|(index, (begin, end))| {
                    Bond::from_spec(
                        BondId::new(index),
                        BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                    )
                })
                .collect(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed Q121 disconnected topology is valid");
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        seed.add_atom(0);

        seed.compute_remaining_size(&query).unwrap();
        assert_eq!(seed.remaining_bonds, 2);
        assert_eq!(seed.remaining_atoms, 2);

        seed.last_added_atoms_begin_index = seed.molecule_fragment.atoms.len();
        seed.compute_remaining_size(&query).unwrap();
        assert_eq!(seed.remaining_bonds, 0);
        assert_eq!(seed.remaining_atoms, 0);
    }

    #[test]
    fn q121a_empty_and_single_masks_preserve_source_termination() {
        let mut empty = McsComposition2N::new(0, 0);
        assert_eq!(empty.bit_set(), 0);
        assert!(empty.is_2_power());
        assert!(!empty.generate_next());
        assert_eq!(empty.bit_set(), 0);

        let single_max = McsComposition2N::compute_2n(1) - 1;
        let mut single = McsComposition2N::new(single_max, single_max);
        assert!(single.generate_next());
        assert_eq!(single.bit_set(), 1);
        assert!(single.is_2_power());
        assert!(single.is_set(0));
        assert!(!single.generate_next());
        assert_eq!(single.bit_set(), 1);
        assert!(!single.generate_next());
        assert_eq!(single.bit_set(), 1);
    }

    #[test]
    fn q121a_multiple_masks_descend_and_retain_unsigned_long_long_width() {
        assert_eq!(McsBitSet::BITS, 64);
        assert_eq!(McsComposition2N::compute_2n(63), 1_u64 << 63);

        let max = McsComposition2N::compute_2n(3) - 1;
        let mut composition = McsComposition2N::new(max, max);
        let mut values = Vec::new();
        let mut single_bit_values = Vec::new();
        while composition.generate_next() {
            values.push(composition.bit_set());
            if composition.is_2_power() {
                single_bit_values.push(composition.bit_set());
            }
        }
        assert_eq!(values, vec![7, 6, 5, 4, 3, 2, 1]);
        assert_eq!(single_bit_values, vec![4, 2, 1]);
        assert!(composition.is_set(0));
        assert!(!composition.is_set(1));

        let wide_max = McsComposition2N::compute_2n(40) - 1;
        let mut wide = McsComposition2N::new(wide_max, wide_max);
        assert!(wide.generate_next());
        assert_eq!(wide.bit_set(), wide_max);
        assert!(wide.is_set(39));
    }

    #[test]
    fn q122_biggest_child_defers_then_resumes_source_ordered_combinations() {
        let topology = TopologyBlock::try_from_parts(
            (0..4)
                .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
                .collect(),
            [(0, 1), (0, 2), (0, 3)]
                .into_iter()
                .enumerate()
                .map(|(index, (begin, end))| {
                    Bond::from_spec(
                        BondId::new(index),
                        BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                    )
                })
                .collect(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed Q122 star topology is valid");
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let target = q107_target(&topology, &coordinates);
        let tables = q126_tables(&query, &target, true, true);
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            remaining_bonds: query.num_bonds(),
            remaining_atoms: query.num_atoms() - 1,
            ..McsSeed::default()
        };
        seed.add_atom(0);
        seed.fill_new_bonds(&query, None, &[], &[], 0, None)
            .unwrap();
        let mut queue = McsSeedQueue::default();

        seed.grow_extensions(
            &mut queue,
            &query,
            &[target],
            &[tables.clone()],
            1,
            0,
            0,
            &McsParameters::default(),
            None,
        )
        .unwrap();
        assert_eq!(seed.growing_stage, 1);
        assert_eq!(queue.seeds.len(), 1);
        assert_eq!(queue.seeds[0].molecule_fragment.bonds, vec![0, 1, 2]);

        seed.grow_extensions(
            &mut queue,
            &query,
            &[target],
            &[tables],
            1,
            0,
            0,
            &McsParameters::default(),
            None,
        )
        .unwrap();
        assert_eq!(seed.growing_stage, u32::MAX);
        assert_eq!(
            queue
                .seeds
                .iter()
                .map(|child| child.molecule_fragment.bonds.clone())
                .collect::<Vec<_>>(),
            vec![
                vec![0, 1, 2],
                vec![1, 2],
                vec![0, 2],
                vec![0, 1],
                vec![0],
                vec![1],
                vec![2],
            ]
        );
    }

    #[test]
    fn q122_single_extension_finishes_and_store_all_preserves_tie_branch() {
        let topology = q113_chain(2);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let target = q107_target(&topology, &coordinates);
        let tables = q126_tables(&query, &target, true, true);
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            remaining_bonds: 1,
            remaining_atoms: 1,
            ..McsSeed::default()
        };
        seed.add_atom(0);
        seed.fill_new_bonds(&query, None, &[], &[], 0, None)
            .unwrap();
        let mut queue = McsSeedQueue::default();
        seed.grow_extensions(
            &mut queue,
            &query,
            &[target],
            &[tables],
            1,
            0,
            0,
            &McsParameters::default(),
            None,
        )
        .unwrap();
        assert_eq!(seed.growing_stage, u32::MAX);
        assert_eq!(queue.seeds.len(), 1);
        assert_eq!(queue.seeds[0].molecule_fragment.bonds, vec![0]);

        let tied = McsSeed {
            store_all_degenerate_mcs: true,
            remaining_bonds: 0,
            remaining_atoms: 0,
            molecule_fragment: queue.seeds[0].molecule_fragment.clone(),
            ..McsSeed::default()
        };
        assert!(tied.can_grow_bigger_than(1, 2));
        let non_store_all = McsSeed {
            store_all_degenerate_mcs: false,
            ..tied
        };
        assert!(!non_store_all.can_grow_bigger_than(1, 2));
    }

    #[test]
    fn q123_complete_ring_growth_restores_fused_options_before_child_match() {
        let topology = q111_cycle(3);
        let ring_info = cosmolkit_core::fast_find_rings(&topology).unwrap();
        let coordinates = CoordinateBlock::default();
        let query = SearchTarget::new(
            &topology,
            &coordinates,
            &topology.stereo_groups,
            Some(&ring_info),
            None,
        );
        let target = SearchTarget::new(
            &topology,
            &coordinates,
            &topology.stereo_groups,
            Some(&ring_info),
            None,
        );
        let tables = q126_tables(&query, &target, true, true);
        let params = McsParameters {
            bond_compare_parameters: McsBondCompareParameters {
                complete_rings_only: true,
                match_fused_rings: true,
                match_fused_rings_strict: true,
                ..McsBondCompareParameters::default()
            },
            ..McsParameters::default()
        };
        let original_params = params.clone();
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            remaining_bonds: query.num_bonds(),
            remaining_atoms: query.num_atoms() - 1,
            ..McsSeed::default()
        };
        seed.add_atom(0);
        let mut queue = McsSeedQueue::default();
        let mut option_states = Vec::new();
        let mut accept = |_: usize, _: &[(usize, usize)], observed: &McsParameters| {
            option_states.push((
                observed.bond_compare_parameters.match_fused_rings,
                observed.bond_compare_parameters.match_fused_rings_strict,
            ));
            Ok(true)
        };

        seed.grow(
            &mut queue,
            &query,
            &[target],
            &[tables],
            1,
            0,
            0,
            &params,
            Some(&mut accept),
        )
        .unwrap();

        assert_eq!(params, original_params);
        assert_eq!(seed.growing_stage, 1);
        assert_eq!(
            seed.new_bonds
                .iter()
                .map(|bond| bond.bond_index)
                .collect::<Vec<_>>(),
            vec![0, 2]
        );
        assert_eq!(queue.seeds.len(), 1);
        assert_eq!(queue.seeds[0].molecule_fragment.bonds, vec![0, 2]);
        assert!(option_states.iter().any(|state| *state == (false, false)));
        assert_eq!(option_states.last(), Some(&(true, true)));
    }

    #[test]
    fn q123_rejected_complete_ring_frontier_finishes_without_partial_mutation() {
        let query_topology = q111_cycle(3);
        let query_rings = cosmolkit_core::fast_find_rings(&query_topology).unwrap();
        let target_topology = q113_chain(3);
        let coordinates = CoordinateBlock::default();
        let query = SearchTarget::new(
            &query_topology,
            &coordinates,
            &query_topology.stereo_groups,
            Some(&query_rings),
            None,
        );
        let target = q107_target(&target_topology, &coordinates);
        let tables = q126_tables(&query, &target, true, true);
        let params = McsParameters {
            bond_compare_parameters: McsBondCompareParameters {
                complete_rings_only: true,
                match_fused_rings: true,
                match_fused_rings_strict: true,
                ..McsBondCompareParameters::default()
            },
            ..McsParameters::default()
        };
        let original_params = params.clone();
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            remaining_bonds: query.num_bonds(),
            remaining_atoms: query.num_atoms() - 1,
            ..McsSeed::default()
        };
        seed.add_atom(0);
        let original_fragment = seed.molecule_fragment.clone();
        let original_excluded = seed.excluded_bonds.clone();
        let mut queue = McsSeedQueue::default();

        seed.grow(
            &mut queue,
            &query,
            &[target],
            &[tables],
            1,
            0,
            0,
            &params,
            None,
        )
        .unwrap();

        assert_eq!(params, original_params);
        assert_eq!(seed.growing_stage, u32::MAX);
        assert!(seed.new_bonds.is_empty());
        assert_eq!(seed.molecule_fragment, original_fragment);
        assert_eq!(seed.excluded_bonds, original_excluded);
        assert!(queue.seeds.is_empty());
    }

    #[test]
    fn q125_target_match_maps_query_rows_to_exact_target_identities() {
        let query_topology = q114_topology();
        let target_topology = TopologyBlock::try_from_parts(
            (0..4)
                .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
                .collect(),
            [(3, 2), (0, 1), (1, 3)]
                .into_iter()
                .enumerate()
                .map(|(index, (begin, end))| {
                    Bond::from_spec(
                        BondId::new(index),
                        BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                    )
                })
                .collect(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed Q125 target topology is valid");
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&query_topology, &coordinates);
        let target = q107_target(&target_topology, &coordinates);
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        for atom in [1, 0, 2] {
            seed.add_atom(atom);
        }
        seed.add_bond(&query, 0).unwrap();
        seed.add_bond(&query, 1).unwrap();

        let mut target_match = McsTargetMatch::default();
        target_match
            .init(&seed, &[(0, 3), (1, 1), (2, 2)], &query, &target)
            .unwrap();

        assert!(!target_match.empty);
        assert_eq!(target_match.matched_atom_size, 3);
        assert_eq!(target_match.matched_bond_size, 2);
        assert_eq!(target_match.target_atom_indices, vec![1, 3, 2]);
        assert_eq!(target_match.target_bond_indices, vec![2, 0]);
        assert_eq!(
            target_match.visited_target_atoms,
            vec![false, true, true, true]
        );
        assert_eq!(target_match.visited_target_bonds, vec![true, false, true]);
    }

    #[test]
    fn q125_target_match_retains_not_set_for_unmatched_query_rows_and_bonds() {
        let query_topology = q110_topology(BondOrder::Single, BondStereo::None);
        let target_topology = q113_chain(3);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&query_topology, &coordinates);
        let target = q107_target(&target_topology, &coordinates);
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        seed.add_atom(0);
        seed.add_atom(1);
        seed.add_bond(&query, 0).unwrap();

        let mut target_match = McsTargetMatch::default();
        target_match
            .init(&seed, &[(0, 0), (1, 2)], &query, &target)
            .unwrap();

        assert!(!target_match.empty);
        assert_eq!(target_match.matched_atom_size, 2);
        assert_eq!(target_match.matched_bond_size, 0);
        assert_eq!(target_match.target_atom_indices, vec![0, 2]);
        assert_eq!(target_match.target_bond_indices, vec![usize::MAX]);
        assert_eq!(target_match.visited_target_atoms, vec![true, false, true]);
        assert_eq!(target_match.visited_target_bonds, vec![false, false]);
    }

    fn q126_tables(
        query: &SearchTarget<'_>,
        target: &SearchTarget<'_>,
        atom_value: bool,
        bond_value: bool,
    ) -> McsMatchTables {
        McsMatchTables {
            atoms: McsMatchTable {
                rows: query.num_atoms(),
                columns: target.num_atoms(),
                values: vec![atom_value; query.num_atoms() * target.num_atoms()],
            },
            bonds: McsMatchTable {
                rows: query.num_bonds(),
                columns: target.num_bonds(),
                values: vec![bond_value; query.num_bonds() * target.num_bonds()],
            },
        }
    }

    #[test]
    fn q126_full_match_uses_source_target_order_and_threshold_stops() {
        let query_topology = q113_chain(2);
        let miss_topology = q107_topology(&[ChiralTag::Unspecified, ChiralTag::Unspecified]);
        let match_topology = q113_chain(2);
        let later_topology = q113_chain(2);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&query_topology, &coordinates);
        let targets = vec![
            q107_target(&miss_topology, &coordinates),
            q107_target(&match_topology, &coordinates),
            q107_target(&later_topology, &coordinates),
        ];
        let tables = targets
            .iter()
            .map(|target| q126_tables(&query, target, true, true))
            .collect::<Vec<_>>();
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        seed.add_atom(0);
        seed.add_atom(1);
        seed.add_bond(&query, 0).unwrap();
        let mut checked_targets = Vec::new();
        let mut final_check = |target: usize, _: &[(usize, usize)]| {
            checked_targets.push(target);
            Ok(true)
        };

        assert_eq!(
            mcs_match_full_candidate(
                &mut seed,
                &query,
                &targets,
                &tables,
                1,
                Some(&McsParameters::default()),
                Some(&mut final_check),
            ),
            Ok(true)
        );
        assert_eq!(checked_targets, vec![1]);
        assert_eq!(seed.match_result.len(), 3);
        assert!(seed.match_result[0].empty);
        assert!(!seed.match_result[1].empty);
        assert!(seed.match_result[2].empty);

        let two_misses = vec![targets[0], targets[0], targets[2]];
        let two_miss_tables = two_misses
            .iter()
            .map(|target| q126_tables(&query, target, true, true))
            .collect::<Vec<_>>();
        let mut failed_seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        failed_seed.add_atom(0);
        failed_seed.add_atom(1);
        failed_seed.add_bond(&query, 0).unwrap();
        assert_eq!(
            mcs_match_full_candidate(
                &mut failed_seed,
                &query,
                &two_misses,
                &two_miss_tables,
                2,
                Some(&McsParameters::default()),
                None,
            ),
            Ok(false)
        );
        assert!(failed_seed.match_result.is_empty());
    }

    #[test]
    fn q126_full_match_requires_exact_topology_and_match_table_feasibility() {
        let query_topology = q113_chain(2);
        let target_topology = q113_chain(3);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&query_topology, &coordinates);
        let target = q107_target(&target_topology, &coordinates);
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        seed.add_atom(0);
        seed.add_atom(1);
        seed.add_bond(&query, 0).unwrap();

        let mut tables = q126_tables(&query, &target, true, false);
        tables.bonds.set(0, 1, true);
        let mut accept = |_: &[(usize, usize)]| Ok(true);
        let mapping = mcs_find_full_mapping(&seed, &target, &tables, &mut accept)
            .unwrap()
            .expect("the second target bond is the sole compatible embedding");
        let mut target_atoms = mapping.iter().map(|&(_, atom)| atom).collect::<Vec<_>>();
        target_atoms.sort_unstable();
        assert_eq!(target_atoms, vec![1, 2]);
        assert_eq!(
            mcs_target_bond_between(&target, mapping[0].1, mapping[1].1),
            Some(1)
        );

        tables.bonds.set(0, 1, false);
        assert_eq!(
            mcs_find_full_mapping(&seed, &target, &tables, &mut accept),
            Ok(None)
        );
    }

    #[test]
    fn q127_incremental_reuse_extends_new_atom_and_old_atom_closure_state() {
        let chain_topology = q113_chain(3);
        let coordinates = CoordinateBlock::default();
        let chain = q107_target(&chain_topology, &coordinates);
        let chain_tables = q126_tables(&chain, &chain, true, true);
        let mut chain_seed = McsSeed {
            excluded_bonds: vec![false; chain.num_bonds()],
            ..McsSeed::default()
        };
        chain_seed.add_atom(0);
        chain_seed.add_atom(1);
        chain_seed.add_bond(&chain, 0).unwrap();
        let mut chain_match = McsTargetMatch::default();
        chain_match
            .init(&chain_seed, &[(0, 0), (1, 1)], &chain, &chain)
            .unwrap();
        chain_seed.add_atom(2);
        chain_seed.add_bond(&chain, 1).unwrap();
        let mut accept = |_: &[(usize, usize)]| Ok(true);

        assert_eq!(
            mcs_match_incremental_fast(
                &chain_seed,
                &mut chain_match,
                &chain,
                &chain_tables,
                &mut accept,
            ),
            Ok(true)
        );
        assert_eq!(chain_match.matched_atom_size, 3);
        assert_eq!(chain_match.matched_bond_size, 2);
        assert_eq!(chain_match.target_atom_indices, vec![0, 1, 2]);
        assert_eq!(chain_match.target_bond_indices, vec![0, 1]);

        let ring_topology = q111_cycle(3);
        let ring = q107_target(&ring_topology, &coordinates);
        let ring_tables = q126_tables(&ring, &ring, true, true);
        let mut ring_seed = McsSeed {
            excluded_bonds: vec![false; ring.num_bonds()],
            ..McsSeed::default()
        };
        for atom in 0..3 {
            ring_seed.add_atom(atom);
        }
        ring_seed.add_bond(&ring, 0).unwrap();
        ring_seed.add_bond(&ring, 1).unwrap();
        let mut ring_match = McsTargetMatch::default();
        ring_match
            .init(&ring_seed, &[(0, 0), (1, 1), (2, 2)], &ring, &ring)
            .unwrap();
        ring_seed.add_bond(&ring, 2).unwrap();

        assert_eq!(
            mcs_match_incremental_fast(
                &ring_seed,
                &mut ring_match,
                &ring,
                &ring_tables,
                &mut accept,
            ),
            Ok(true)
        );
        assert_eq!(ring_match.matched_atom_size, 3);
        assert_eq!(ring_match.matched_bond_size, 3);
        assert_eq!(ring_match.target_bond_indices, vec![0, 1, 2]);
    }

    #[test]
    fn q127_incremental_rejection_clears_cache_then_uses_source_full_fallback() {
        let topology = q113_chain(3);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let target = q107_target(&topology, &coordinates);
        let tables = q126_tables(&query, &target, true, true);
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            ..McsSeed::default()
        };
        seed.add_atom(0);
        seed.add_atom(1);
        seed.add_bond(&query, 0).unwrap();
        let mut cached = McsTargetMatch::default();
        cached
            .init(&seed, &[(0, 0), (1, 1)], &query, &target)
            .unwrap();
        seed.match_result.push(cached);
        seed.add_atom(2);
        seed.add_bond(&query, 1).unwrap();
        let mut checks = 0;
        let mut reject_incremental = |target_index: usize, _: &[(usize, usize)]| {
            assert_eq!(target_index, 0);
            checks += 1;
            Ok(checks == 2)
        };

        assert_eq!(
            mcs_match_full_candidate(
                &mut seed,
                &query,
                &[target],
                &[tables],
                1,
                Some(&McsParameters::default()),
                Some(&mut reject_incremental),
            ),
            Ok(true)
        );
        assert_eq!(checks, 2);
        assert!(!seed.match_result[0].empty);
        assert_eq!(seed.match_result[0].matched_atom_size, 3);
        assert_eq!(seed.match_result[0].matched_bond_size, 2);

        let mut impossible_tables = q126_tables(&query, &target, true, true);
        impossible_tables.bonds.set(1, 1, false);
        let mut invalidated = seed.match_result[0].clone();
        invalidated.matched_atom_size = 2;
        invalidated.matched_bond_size = 1;
        invalidated.target_atom_indices[2] = usize::MAX;
        invalidated.target_bond_indices[1] = usize::MAX;
        invalidated.visited_target_atoms[2] = false;
        invalidated.visited_target_bonds[1] = false;
        let mut accept = |_: &[(usize, usize)]| Ok(true);
        assert_eq!(
            mcs_match_incremental_fast(
                &seed,
                &mut invalidated,
                &target,
                &impossible_tables,
                &mut accept,
            ),
            Ok(false)
        );
        assert!(invalidated.empty);
        assert!(invalidated.target_atom_indices.is_empty());
        assert!(invalidated.target_bond_indices.is_empty());
        assert!(invalidated.visited_target_atoms.is_empty());
        assert!(invalidated.visited_target_bonds.is_empty());
    }

    #[test]
    fn q128_acceptance_applies_threshold_and_final_check_before_queue_insert() {
        let query_topology = q113_chain(2);
        let match_topology = q113_chain(2);
        let miss_topology = q107_topology(&[ChiralTag::Unspecified, ChiralTag::Unspecified]);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&query_topology, &coordinates);
        let targets = vec![
            q107_target(&match_topology, &coordinates),
            q107_target(&miss_topology, &coordinates),
        ];
        let tables = targets
            .iter()
            .map(|target| q126_tables(&query, target, true, true))
            .collect::<Vec<_>>();
        let make_seed = || {
            let mut seed = McsSeed {
                excluded_bonds: vec![false; query.num_bonds()],
                ..McsSeed::default()
            };
            seed.add_atom(0);
            seed.add_atom(1);
            seed.add_bond(&query, 0).unwrap();
            seed
        };
        let mut queue = McsSeedQueue::default();

        let mut threshold_seed = make_seed();
        assert_eq!(
            mcs_check_if_match_and_append(
                &mut threshold_seed,
                &mut queue,
                &query,
                &targets,
                &tables,
                2,
                Some(&McsParameters::default()),
                None,
            ),
            Ok(false)
        );
        assert!(queue.seeds.is_empty());

        let mut rejected_seed = make_seed();
        let mut reject = |_: usize, _: &[(usize, usize)]| Ok(false);
        assert_eq!(
            mcs_check_if_match_and_append(
                &mut rejected_seed,
                &mut queue,
                &query,
                &targets,
                &tables,
                1,
                Some(&McsParameters::default()),
                Some(&mut reject),
            ),
            Ok(false)
        );
        assert!(queue.seeds.is_empty());

        let mut accepted_seed = make_seed();
        let mut accept = |target: usize, mapping: &[(usize, usize)]| {
            assert_eq!(target, 0);
            assert_eq!(mapping.len(), 2);
            Ok(true)
        };
        assert_eq!(
            mcs_check_if_match_and_append(
                &mut accepted_seed,
                &mut queue,
                &query,
                &targets,
                &tables,
                1,
                Some(&McsParameters::default()),
                Some(&mut accept),
            ),
            Ok(true)
        );
        assert_eq!(queue.seeds.len(), 1);
        assert!(queue.seeds[0].copy_complete);
        assert_eq!(queue.seeds[0].match_result, accepted_seed.match_result);
    }

    #[test]
    fn q128_seed_queue_orders_descending_bonds_and_retains_equal_insertion_order() {
        let queue_seed = |identity: usize, bond_count: usize| McsSeed {
            copy_complete: false,
            molecule_fragment: McsMoleculeFragment {
                atoms: vec![identity],
                bonds: (0..bond_count).collect(),
                seed_atom_index_map: BTreeMap::new(),
            },
            ..McsSeed::default()
        };
        let mut queue = McsSeedQueue::default();
        queue.add(&queue_seed(1, 1));
        queue.add(&queue_seed(3, 3));
        queue.add(&queue_seed(2, 1));

        assert_eq!(
            queue
                .seeds
                .iter()
                .map(|seed| seed.molecule_fragment.atoms[0])
                .collect::<Vec<_>>(),
            vec![3, 1, 2]
        );
        assert_eq!(
            queue
                .seeds
                .iter()
                .map(|seed| seed.molecule_fragment.bonds.len())
                .collect::<Vec<_>>(),
            vec![3, 1, 1]
        );
        assert!(queue.seeds.iter().all(|seed| seed.copy_complete));
    }

    #[test]
    fn q128a_uses_first_ring_membership_connected_order_and_relaxed_fusion_flags() {
        let query_topology = q111_fused_triangles();
        let query_rings = cosmolkit_core::fast_find_rings(&query_topology).unwrap();
        let source_bond = 0;
        let memberships = query_rings.bond_members(BondId::new(source_bond));
        assert_eq!(memberships.len(), 2);
        let first_ring = &query_rings.bond_rings()[memberships[0]];
        assert_eq!(first_ring.len(), 3);

        let mut first_ring_atoms = BTreeSet::new();
        for bond in first_ring {
            let bond = &query_topology.bonds[bond.index()];
            first_ring_atoms.insert(bond.begin().index());
            first_ring_atoms.insert(bond.end().index());
        }
        assert_eq!(first_ring_atoms.len(), 3);
        let target_topology = q111_cycle(first_ring_atoms.len());
        let coordinates = CoordinateBlock::default();
        let query = SearchTarget::new(
            &query_topology,
            &coordinates,
            &query_topology.stereo_groups,
            Some(&query_rings),
            None,
        );
        let target = q107_target(&target_topology, &coordinates);
        let targets = vec![target];
        let mut tables = q126_tables(&query, &targets[0], true, true);
        for query_atom in 0..query.num_atoms() {
            if !first_ring_atoms.contains(&query_atom) {
                for target_atom in 0..targets[0].num_atoms() {
                    tables.atoms.set(query_atom, target_atom, false);
                }
            }
        }

        let source_atom = query_topology.bonds[source_bond].begin().index();
        let mut parent = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            remaining_bonds: query.num_bonds(),
            remaining_atoms: query.num_atoms() - 1,
            ..McsSeed::default()
        };
        parent.add_atom(source_atom);
        let params = McsParameters {
            bond_compare_parameters: McsBondCompareParameters {
                match_fused_rings: true,
                match_fused_rings_strict: true,
                ..McsBondCompareParameters::default()
            },
            ..McsParameters::default()
        };
        let mut checked = 0;
        let mut final_check =
            |target_index: usize, mapping: &[(usize, usize)], seen: &McsParameters| {
                checked += 1;
                assert_eq!(target_index, 0);
                assert_eq!(mapping.len(), first_ring_atoms.len());
                assert!(!seen.bond_compare_parameters.match_fused_rings);
                assert!(!seen.bond_compare_parameters.match_fused_rings_strict);
                Ok(true)
            };

        assert_eq!(
            mcs_can_add_all_non_fused_ring_bonds_connected_to_bond(
                &parent,
                &query,
                source_atom,
                source_bond,
                &targets,
                &[tables],
                1,
                &params,
                Some(&mut final_check),
            ),
            Ok(true)
        );
        assert_eq!(checked, 1);
        assert!(params.bond_compare_parameters.match_fused_rings);
        assert!(params.bond_compare_parameters.match_fused_rings_strict);
    }

    #[test]
    fn q128a_preserves_missing_ring_non_ring_endpoint_and_empty_frontier_failures() {
        let chain_topology = q113_chain(3);
        let chain_rings = cosmolkit_core::fast_find_rings(&chain_topology).unwrap();
        let coordinates = CoordinateBlock::default();
        let plain_chain = q107_target(&chain_topology, &coordinates);
        let ringed_chain = SearchTarget::new(
            &chain_topology,
            &coordinates,
            &chain_topology.stereo_groups,
            Some(&chain_rings),
            None,
        );
        let parent = McsSeed {
            excluded_bonds: vec![false; chain_topology.bonds.len()],
            remaining_bonds: chain_topology.bonds.len(),
            remaining_atoms: chain_topology.atoms.len() - 1,
            ..McsSeed::default()
        };
        assert_eq!(
            mcs_can_add_all_non_fused_ring_bonds_connected_to_bond(
                &parent,
                &plain_chain,
                0,
                0,
                &[],
                &[],
                0,
                &McsParameters::default(),
                None,
            ),
            Err(McsCandidateMatchError::State(McsError::MissingRingInfo {
                side: "query",
            }))
        );
        assert_eq!(
            mcs_can_add_all_non_fused_ring_bonds_connected_to_bond(
                &parent,
                &ringed_chain,
                0,
                0,
                &[],
                &[],
                0,
                &McsParameters::default(),
                None,
            ),
            Err(McsCandidateMatchError::RingMembershipMissing { bond: 0 })
        );

        let ring_topology = q111_cycle(3);
        let ring_info = cosmolkit_core::fast_find_rings(&ring_topology).unwrap();
        let ring = SearchTarget::new(
            &ring_topology,
            &coordinates,
            &ring_topology.stereo_groups,
            Some(&ring_info),
            None,
        );
        let mut endpoint_parent = McsSeed {
            excluded_bonds: vec![false; ring.num_bonds()],
            remaining_bonds: ring.num_bonds(),
            remaining_atoms: ring.num_atoms() - 1,
            ..McsSeed::default()
        };
        endpoint_parent.add_atom(2);
        assert_eq!(
            mcs_can_add_all_non_fused_ring_bonds_connected_to_bond(
                &endpoint_parent,
                &ring,
                2,
                0,
                &[],
                &[],
                0,
                &McsParameters::default(),
                None,
            ),
            Err(McsCandidateMatchError::State(
                McsError::SeedAtomNotBondEndpoint { atom: 2, bond: 0 }
            ))
        );

        let mut excluded_parent = endpoint_parent;
        excluded_parent.excluded_bonds[0] = true;
        assert_eq!(
            mcs_can_add_all_non_fused_ring_bonds_connected_to_bond(
                &excluded_parent,
                &ring,
                0,
                0,
                &[],
                &[],
                0,
                &McsParameters::default(),
                None,
            ),
            Ok(false)
        );
    }

    #[test]
    fn q120_frontier_uses_last_added_atom_and_source_adjacency_order_once() {
        let topology = q113_chain(4);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let mut seed = McsSeed {
            excluded_bonds: vec![false; query.num_bonds()],
            last_added_atoms_begin_index: 0,
            remaining_bonds: query.num_bonds(),
            remaining_atoms: query.num_atoms() - 2,
            ..McsSeed::default()
        };
        seed.add_atom(1);
        seed.add_atom(2);

        seed.fill_new_bonds(&query, None, &[], &[], 0, None)
            .unwrap();
        assert_eq!(
            seed.new_bonds
                .iter()
                .map(|bond| bond.bond_index)
                .collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
        assert_eq!(seed.new_bonds[0].new_atom_index, 0);
        assert_eq!(seed.new_bonds[1].new_atom_index, 2);
        assert_eq!(seed.new_bonds[1].end_atom_index, Some(1));
        assert_eq!(seed.new_bonds[2].new_atom_index, 3);

        let mut last_only = McsSeed {
            excluded_bonds: vec![false, true, false],
            last_added_atoms_begin_index: 1,
            ..McsSeed::default()
        };
        last_only.add_atom(1);
        last_only.add_atom(2);
        last_only
            .fill_new_bonds(&query, None, &[], &[], 0, None)
            .unwrap();
        assert_eq!(
            last_only
                .new_bonds
                .iter()
                .map(|bond| bond.bond_index)
                .collect::<Vec<_>>(),
            vec![2]
        );
    }

    #[test]
    fn q120_complete_ring_gate_eliminates_only_source_rejected_candidates() {
        let query_topology = q111_cycle(3);
        let query_rings = cosmolkit_core::fast_find_rings(&query_topology).unwrap();
        let miss_topology = q113_chain(3);
        let hit_topology = q111_cycle(3);
        let coordinates = CoordinateBlock::default();
        let query = SearchTarget::new(
            &query_topology,
            &coordinates,
            &query_topology.stereo_groups,
            Some(&query_rings),
            None,
        );
        let miss = q107_target(&miss_topology, &coordinates);
        let hit = q107_target(&hit_topology, &coordinates);
        let params = McsParameters {
            bond_compare_parameters: McsBondCompareParameters {
                complete_rings_only: true,
                ..McsBondCompareParameters::default()
            },
            ..McsParameters::default()
        };
        let make_seed = || {
            let mut seed = McsSeed {
                excluded_bonds: vec![false; query.num_bonds()],
                remaining_bonds: query.num_bonds(),
                remaining_atoms: query.num_atoms() - 1,
                ..McsSeed::default()
            };
            seed.add_atom(0);
            seed
        };

        let mut rejected = make_seed();
        rejected
            .fill_new_bonds(
                &query,
                Some(&params),
                &[miss],
                &[q126_tables(&query, &miss, true, true)],
                1,
                None,
            )
            .unwrap();
        assert!(rejected.new_bonds.is_empty());

        let mut accepted = make_seed();
        accepted
            .fill_new_bonds(
                &query,
                Some(&params),
                &[hit],
                &[q126_tables(&query, &hit, true, true)],
                1,
                None,
            )
            .unwrap();
        assert_eq!(
            accepted
                .new_bonds
                .iter()
                .map(|bond| bond.bond_index)
                .collect::<Vec<_>>(),
            vec![0, 2]
        );
    }

    #[test]
    fn q118_initial_smarts_seeds_precede_fallback_and_retain_match_order() {
        let query_topology = q113_chain(4);
        let target_topology = q113_chain(4);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&query_topology, &coordinates);
        let target = q107_target(&target_topology, &coordinates);
        let tables = q126_tables(&query, &target, true, true);
        let params = McsParameters {
            initial_seed: "C-C-C".to_owned(),
            ..McsParameters::default()
        };

        let initial =
            mcs_make_initial_seeds(&query, &[target], &[tables], 1, &params, false, None, None)
                .unwrap();
        assert_eq!(initial.query_matched_bonds, 2);
        assert_eq!(initial.queue.seeds.len(), 2);
        assert_eq!(
            initial
                .queue
                .seeds
                .iter()
                .map(|seed| seed.molecule_fragment.atoms.clone())
                .collect::<Vec<_>>(),
            vec![vec![0, 1, 2], vec![1, 2, 3]]
        );
        assert_eq!(
            initial
                .queue
                .seeds
                .iter()
                .map(|seed| seed.molecule_fragment.bonds.clone())
                .collect::<Vec<_>>(),
            vec![vec![0, 1], vec![1, 2]]
        );
        assert_eq!(initial.query_matched_atoms, 4);
        assert_eq!(initial.query_single_matched_atom, Some(2));
    }

    #[test]
    fn q118_fallback_keeps_source_bond_order_propagates_rejection_and_reports_parse_errors() {
        let query_topology = q113_chain(4);
        let target_topology = q113_chain(4);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&query_topology, &coordinates);
        let target = q107_target(&target_topology, &coordinates);
        let mut tables = q126_tables(&query, &target, true, true);
        for target_bond in 0..target.num_bonds() {
            tables.bonds.set(1, target_bond, false);
        }
        let fallback_params = McsParameters {
            initial_seed: "N".to_owned(),
            ..McsParameters::default()
        };
        let fallback = mcs_make_initial_seeds(
            &query,
            &[target],
            &[tables],
            1,
            &fallback_params,
            false,
            None,
            None,
        )
        .unwrap();
        assert_eq!(fallback.query_matched_bonds, 2);
        assert_eq!(
            fallback
                .queue
                .seeds
                .iter()
                .map(|seed| seed.molecule_fragment.bonds.clone())
                .collect::<Vec<_>>(),
            vec![vec![0], vec![2]]
        );
        assert_eq!(
            fallback.queue.seeds[0].excluded_bonds,
            vec![true, true, false]
        );
        assert_eq!(
            fallback.queue.seeds[1].excluded_bonds,
            vec![true, true, true]
        );

        let invalid = McsParameters {
            initial_seed: "[".to_owned(),
            ..McsParameters::default()
        };
        assert!(matches!(
            mcs_make_initial_seeds(
                &query,
                &[target],
                &[q126_tables(&query, &target, true, true)],
                1,
                &invalid,
                false,
                None,
                None,
            ),
            Err(McsCandidateMatchError::InitialSeedParse { .. })
        ));
    }

    #[test]
    fn q129_final_accept_uses_source_option_gates_and_check_order() {
        let calls = std::cell::RefCell::new(Vec::new());
        let mut ring = || -> Result<bool, McsError> {
            calls.borrow_mut().push("ring");
            Ok(true)
        };
        let mut chirality = || {
            calls.borrow_mut().push("chirality");
            Ok(true)
        };
        let mut user = || {
            calls.borrow_mut().push("user");
            Ok(true)
        };
        assert_eq!(
            mcs_final_candidate_accept(
                &McsParameters::default(),
                &mut ring,
                &mut chirality,
                Some(&mut user),
            ),
            Ok(true)
        );
        assert_eq!(*calls.borrow(), vec!["user"]);

        calls.borrow_mut().clear();
        let strict = McsParameters {
            atom_compare_parameters: McsAtomCompareParameters {
                match_chiral_tag: true,
                ..McsAtomCompareParameters::default()
            },
            bond_compare_parameters: McsBondCompareParameters {
                match_fused_rings_strict: true,
                ..McsBondCompareParameters::default()
            },
            ..McsParameters::default()
        };
        assert_eq!(
            mcs_final_candidate_accept(&strict, &mut ring, &mut chirality, Some(&mut user),),
            Ok(true)
        );
        assert_eq!(*calls.borrow(), vec!["ring", "chirality", "user"]);
    }

    #[test]
    fn q129_final_accept_preserves_rejection_and_error_short_circuits() {
        let calls = std::cell::RefCell::new(Vec::new());
        let params = McsParameters {
            atom_compare_parameters: McsAtomCompareParameters {
                match_chiral_tag: true,
                ..McsAtomCompareParameters::default()
            },
            bond_compare_parameters: McsBondCompareParameters {
                match_fused_rings: true,
                ..McsBondCompareParameters::default()
            },
            ..McsParameters::default()
        };
        let mut ring_reject = || {
            calls.borrow_mut().push("ring");
            Ok(false)
        };
        let mut chirality = || {
            calls.borrow_mut().push("chirality");
            Ok(true)
        };
        let mut user = || {
            calls.borrow_mut().push("user");
            Ok(true)
        };
        assert_eq!(
            mcs_final_candidate_accept(&params, &mut ring_reject, &mut chirality, Some(&mut user),),
            Ok(false)
        );
        assert_eq!(*calls.borrow(), vec!["ring"]);

        calls.borrow_mut().clear();
        let mut ring_error = || {
            calls.borrow_mut().push("ring");
            Err(McsError::MissingRingInfo { side: "left" })
        };
        assert_eq!(
            mcs_final_candidate_accept(&params, &mut ring_error, &mut chirality, Some(&mut user),),
            Err(McsError::MissingRingInfo { side: "left" })
        );
        assert_eq!(*calls.borrow(), vec!["ring"]);

        calls.borrow_mut().clear();
        let no_ring = McsParameters {
            atom_compare_parameters: McsAtomCompareParameters {
                match_chiral_tag: true,
                ..McsAtomCompareParameters::default()
            },
            ..McsParameters::default()
        };
        let mut disabled_ring = || {
            calls.borrow_mut().push("ring");
            Ok(false)
        };
        let mut chirality_reject = || {
            calls.borrow_mut().push("chirality");
            Ok(false)
        };
        assert_eq!(
            mcs_final_candidate_accept(
                &no_ring,
                &mut disabled_ring,
                &mut chirality_reject,
                Some(&mut user),
            ),
            Ok(false)
        );
        assert_eq!(*calls.borrow(), vec!["chirality"]);
    }

    #[test]
    fn q130_timeout_is_disabled_at_zero_and_inclusive_at_deadline() {
        assert!(mcs_progress_callback_timeout(
            &McsParameters::default(),
            10,
            u64::MAX,
        ));

        let params = McsParameters {
            timeout: 2,
            ..McsParameters::default()
        };
        assert!(mcs_progress_callback_timeout(&params, 10, 2_000_010));
        assert!(!mcs_progress_callback_timeout(&params, 10, 2_000_011));
    }

    #[test]
    fn q130_progress_updates_best_before_callback_and_preserves_partial_result() {
        let params = McsParameters {
            timeout: 1,
            ..McsParameters::default()
        };
        let best = McsMoleculeFragment {
            atoms: vec![4, 2],
            bonds: vec![7],
            seed_atom_index_map: BTreeMap::new(),
        };
        let mut progress = McsProgressData {
            seed_processed: 9,
            ..McsProgressData::default()
        };
        let mut observed = None;
        let mut cancel = |data: &McsProgressData, seen_params: &McsParameters| {
            observed = Some(*data);
            assert_eq!(seen_params.timeout, 1);
            Ok(false)
        };

        let outcome = mcs_progress_after_seed(
            &params,
            0,
            1_000_001,
            &best,
            &mut progress,
            Some(&mut cancel),
        )
        .unwrap();

        assert_eq!(
            observed,
            Some(McsProgressData {
                num_atoms: 2,
                num_bonds: 1,
                seed_processed: 9,
            })
        );
        assert!(outcome.canceled);
        assert_eq!(outcome.best_atoms, [4, 2]);
        assert_eq!(outcome.best_bonds, [7]);
        assert_eq!(outcome.progress, progress);
    }

    #[test]
    fn q130_progress_callback_error_follows_count_update() {
        let best = McsMoleculeFragment {
            atoms: vec![3],
            bonds: vec![8, 5],
            seed_atom_index_map: BTreeMap::new(),
        };
        let mut progress = McsProgressData::default();
        let mut callback = |data: &McsProgressData, _: &McsParameters| {
            assert_eq!(data.num_atoms, 1);
            assert_eq!(data.num_bonds, 2);
            Err(McsProgressError::Callback {
                message: "boom".to_owned(),
            })
        };

        assert_eq!(
            mcs_progress_after_seed(
                &McsParameters::default(),
                0,
                0,
                &best,
                &mut progress,
                Some(&mut callback),
            ),
            Err(McsProgressError::Callback {
                message: "boom".to_owned(),
            })
        );
        assert_eq!(progress.num_atoms, 1);
        assert_eq!(progress.num_bonds, 2);
        assert_eq!(progress.seed_processed, 0);
        assert_eq!(best.atoms, vec![3]);
        assert_eq!(best.bonds, vec![8, 5]);
    }

    fn q134_state_seed(atoms: &[usize], bonds: &[usize], source_bonds: usize) -> McsSeed {
        McsSeed {
            copy_complete: true,
            growing_stage: u32::MAX,
            molecule_fragment: McsMoleculeFragment {
                atoms: atoms.to_vec(),
                bonds: bonds.to_vec(),
                seed_atom_index_map: BTreeMap::new(),
            },
            excluded_bonds: vec![true; source_bonds],
            ..McsSeed::default()
        }
    }

    pub(super) fn q134_state_best_and_provenance_survive_later_worse_empty_and_canceled_calls() {
        let topology = q113_chain(4);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let first_context = McsResultContext {
            query_input: 5,
            target_inputs: vec![1, 2],
        };
        let later_context = McsResultContext {
            query_input: 7,
            target_inputs: vec![3, 4],
        };
        let mut state = McsSearchState::default();
        let mut now = || 0;
        let mut first = McsSeedQueue {
            seeds: vec![q134_state_seed(&[0, 1], &[0, 1], 3)],
        };
        let found = mcs_grow_seeds(
            &mut first,
            &mut state,
            &first_context,
            &query,
            &[],
            &[],
            0,
            3,
            &McsParameters::default(),
            0,
            &mut now,
            None,
            None,
            None,
        )
        .unwrap();
        assert!(found.mcs_found);
        assert!(!found.canceled);
        assert_eq!(state.best.bonds, [0, 1]);
        assert_eq!(state.best_context, Some(first_context.clone()));

        let mut worse = McsSeedQueue {
            seeds: vec![q134_state_seed(&[2, 3], &[2], 3)],
        };
        let mut seen_progress = Vec::new();
        let mut cancel = |progress: &McsProgressData, _: &McsParameters| {
            seen_progress.push((progress.num_atoms, progress.num_bonds));
            Ok(false)
        };
        let second = mcs_grow_seeds(
            &mut worse,
            &mut state,
            &later_context,
            &query,
            &[],
            &[],
            0,
            3,
            &McsParameters::default(),
            0,
            &mut now,
            None,
            None,
            Some(&mut cancel),
        )
        .unwrap();
        assert!(!second.mcs_found);
        assert!(second.canceled);
        assert_eq!(seen_progress, [(2, 2)]);
        assert_eq!(state.best.bonds, [0, 1]);
        assert_eq!(state.best_context, Some(first_context.clone()));
        assert_eq!((state.progress.num_atoms, state.progress.num_bonds), (2, 2));

        let mut empty = McsSeedQueue::default();
        let third = mcs_grow_seeds(
            &mut empty,
            &mut state,
            &later_context,
            &query,
            &[],
            &[],
            0,
            3,
            &McsParameters::default(),
            0,
            &mut now,
            None,
            None,
            None,
        )
        .unwrap();
        assert!(!third.mcs_found);
        assert!(!third.canceled);
        assert_eq!(state.best_context, Some(first_context));

        let mut improved = McsSeedQueue {
            seeds: vec![q134_state_seed(&[0, 1, 2], &[0, 1, 2], 3)],
        };
        let mut cancel_after_accept = |progress: &McsProgressData, _: &McsParameters| {
            assert_eq!((progress.num_atoms, progress.num_bonds), (3, 3));
            Ok(false)
        };
        let fourth = mcs_grow_seeds(
            &mut improved,
            &mut state,
            &later_context,
            &query,
            &[],
            &[],
            0,
            3,
            &McsParameters::default(),
            0,
            &mut now,
            None,
            None,
            Some(&mut cancel_after_accept),
        )
        .unwrap();
        assert!(fourth.mcs_found);
        assert!(fourth.canceled);
        assert_eq!(state.best.bonds, [0, 1, 2]);
        assert_eq!(state.best_context, Some(later_context));
    }

    pub(super) fn q134_state_store_all_keeps_tie_context_and_source_key_overwrite() {
        let topology = q113_chain(4);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let params = McsParameters {
            store_all: true,
            ..McsParameters::default()
        };
        let mut state = McsSearchState::default();
        let mut now = || 0;
        for (query_input, targets, atoms, bonds) in [
            (0, vec![1, 2], vec![0, 1], vec![0]),
            (2, vec![0, 1], vec![1, 2], vec![1]),
            (1, vec![2, 0], vec![2, 3], vec![0]),
        ] {
            let context = McsResultContext {
                query_input,
                target_inputs: targets,
            };
            let mut queue = McsSeedQueue {
                seeds: vec![q134_state_seed(&atoms, &bonds, 3)],
            };
            let call = mcs_grow_seeds(
                &mut queue,
                &mut state,
                &context,
                &query,
                &[],
                &[],
                0,
                3,
                &params,
                0,
                &mut now,
                None,
                None,
                None,
            )
            .unwrap();
            assert!(call.mcs_found);
            assert_eq!(state.best_context, Some(context));
        }
        assert_eq!(state.degenerate.len(), 2);
        assert_eq!(state.degenerate[&vec![0]].fragment.atoms, [2, 3]);
        assert_eq!(state.degenerate[&vec![0]].context.query_input, 1);
        assert_eq!(state.degenerate[&vec![0]].context.target_inputs, [2, 0]);
        assert_eq!(state.degenerate[&vec![1]].context.query_input, 2);

        let before_rejection = state.clone();
        let mut rejected = McsSeedQueue {
            seeds: vec![q134_state_seed(&[0, 1, 2], &[0, 1], 3)],
        };
        let mut reject = |_: &McsMoleculeFragment, _: &McsParameters| Ok(false);
        let rejected_call = mcs_grow_seeds(
            &mut rejected,
            &mut state,
            &McsResultContext {
                query_input: 6,
                target_inputs: vec![0, 1],
            },
            &query,
            &[],
            &[],
            0,
            3,
            &params,
            0,
            &mut now,
            None,
            Some(&mut reject),
            None,
        )
        .unwrap();
        assert!(!rejected_call.mcs_found);
        assert_eq!(state.best, before_rejection.best);
        assert_eq!(state.best_context, before_rejection.best_context);
        assert_eq!(state.degenerate, before_rejection.degenerate);

        let improved_context = McsResultContext {
            query_input: 3,
            target_inputs: vec![0, 2],
        };
        let mut improved = McsSeedQueue {
            seeds: vec![q134_state_seed(&[0, 1, 2], &[0, 1], 3)],
        };
        mcs_grow_seeds(
            &mut improved,
            &mut state,
            &improved_context,
            &query,
            &[],
            &[],
            0,
            3,
            &params,
            0,
            &mut now,
            None,
            None,
            None,
        )
        .unwrap();
        assert_eq!(state.degenerate.len(), 1);
        assert_eq!(state.degenerate[&vec![0, 1]].context, improved_context);
    }

    pub(super) fn q134_state_objective_and_zero_bond_reset_follow_source() {
        let topology = q113_chain(4);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let first_context = McsResultContext {
            query_input: 0,
            target_inputs: vec![1],
        };
        let second_context = McsResultContext {
            query_input: 1,
            target_inputs: vec![0],
        };
        for maximize_bonds in [false, true] {
            let params = McsParameters {
                maximize_bonds,
                ..McsParameters::default()
            };
            let mut state = McsSearchState::default();
            let mut now = || 0;
            for (context, atoms, bonds) in [
                (&first_context, vec![0, 1, 2], vec![0]),
                (&second_context, vec![0, 1], vec![0, 1]),
            ] {
                let mut queue = McsSeedQueue {
                    seeds: vec![q134_state_seed(&atoms, &bonds, 3)],
                };
                mcs_grow_seeds(
                    &mut queue,
                    &mut state,
                    context,
                    &query,
                    &[],
                    &[],
                    0,
                    3,
                    &params,
                    0,
                    &mut now,
                    None,
                    None,
                    None,
                )
                .unwrap();
            }
            if maximize_bonds {
                assert_eq!(state.best.bonds, [0, 1]);
                assert_eq!(state.best_context, Some(second_context.clone()));
            } else {
                assert_eq!(state.best.bonds, [0]);
                assert_eq!(state.best_context, Some(first_context.clone()));
            }
        }

        let mut state = McsSearchState {
            best: McsMoleculeFragment {
                atoms: vec![2],
                ..McsMoleculeFragment::default()
            },
            best_context: Some(first_context),
            degenerate: BTreeMap::from([(
                vec![],
                McsRetainedResult {
                    fragment: McsMoleculeFragment::default(),
                    context: second_context,
                },
            )]),
            ..McsSearchState::default()
        };
        state.clear_best_for_zero_bond_fallback();
        assert!(state.best.atoms.is_empty());
        assert!(state.best_context.is_none());
        assert_eq!(state.degenerate.len(), 1);
    }

    #[test]
    fn q124_seed_queue_respects_atom_and_bond_objectives_after_pruning() {
        let topology = q113_chain(4);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let make_seed = |atoms: Vec<usize>, bonds: Vec<usize>| McsSeed {
            copy_complete: true,
            growing_stage: u32::MAX,
            molecule_fragment: McsMoleculeFragment {
                atoms,
                bonds,
                seed_atom_index_map: BTreeMap::new(),
            },
            excluded_bonds: vec![true; query.num_bonds()],
            ..McsSeed::default()
        };
        let make_queue = || McsSeedQueue {
            seeds: vec![
                make_seed(vec![0, 1], vec![0, 1]),
                make_seed(vec![0, 1, 2], vec![2]),
            ],
        };

        let mut atom_queue = make_queue();
        let mut atom_state = McsSearchState::default();
        let context = McsResultContext {
            query_input: 0,
            target_inputs: Vec::new(),
        };
        let mut accepted = Vec::new();
        let mut accept = |fragment: &McsMoleculeFragment, _: &McsParameters| {
            accepted.push((fragment.atoms.len(), fragment.bonds.len()));
            Ok(true)
        };
        let mut progress = Vec::new();
        let mut progress_callback = |data: &McsProgressData, _: &McsParameters| {
            progress.push((data.num_atoms, data.num_bonds));
            Ok(true)
        };
        let mut now = || 0;
        let atom_result = mcs_grow_seeds(
            &mut atom_queue,
            &mut atom_state,
            &context,
            &query,
            &[],
            &[],
            0,
            query.num_bonds(),
            &McsParameters {
                maximize_bonds: false,
                ..McsParameters::default()
            },
            0,
            &mut now,
            None,
            Some(&mut accept),
            Some(&mut progress_callback),
        )
        .unwrap();
        assert_eq!(accepted, vec![(2, 2), (3, 1)]);
        assert_eq!(progress, vec![(2, 2), (3, 1)]);
        assert_eq!(atom_state.best.atoms, vec![0, 1, 2]);
        assert_eq!(atom_state.best.bonds, vec![2]);
        assert!(atom_result.mcs_found);
        assert!(!atom_result.canceled);
        assert!(atom_queue.seeds.is_empty());

        let mut bond_queue = make_queue();
        let mut bond_state = McsSearchState::default();
        let mut now = || 0;
        let bond_result = mcs_grow_seeds(
            &mut bond_queue,
            &mut bond_state,
            &context,
            &query,
            &[],
            &[],
            0,
            query.num_bonds(),
            &McsParameters::default(),
            0,
            &mut now,
            None,
            None,
            None,
        )
        .unwrap();
        assert_eq!(bond_state.best.atoms, vec![0, 1]);
        assert_eq!(bond_state.best.bonds, vec![0, 1]);
        assert!(bond_result.mcs_found);
        assert!(bond_queue.seeds.is_empty());
    }

    #[test]
    fn q124_store_all_retains_ties_and_query_bound_stops_after_progress() {
        let topology = q113_chain(3);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let make_seed = |atoms: Vec<usize>, bonds: Vec<usize>| McsSeed {
            copy_complete: true,
            growing_stage: u32::MAX,
            molecule_fragment: McsMoleculeFragment {
                atoms,
                bonds,
                seed_atom_index_map: BTreeMap::new(),
            },
            excluded_bonds: vec![true; query.num_bonds()],
            store_all_degenerate_mcs: true,
            ..McsSeed::default()
        };
        let mut tie_queue = McsSeedQueue {
            seeds: vec![
                make_seed(vec![0, 1], vec![0]),
                make_seed(vec![1, 2], vec![1]),
            ],
        };
        let mut tie_state = McsSearchState::default();
        let context = McsResultContext {
            query_input: 0,
            target_inputs: Vec::new(),
        };
        let mut accepted = 0;
        let mut accept = |_: &McsMoleculeFragment, _: &McsParameters| {
            accepted += 1;
            Ok(true)
        };
        let mut now = || 0;
        let tie_result = mcs_grow_seeds(
            &mut tie_queue,
            &mut tie_state,
            &context,
            &query,
            &[],
            &[],
            0,
            query.num_bonds(),
            &McsParameters {
                store_all: true,
                ..McsParameters::default()
            },
            0,
            &mut now,
            None,
            Some(&mut accept),
            None,
        )
        .unwrap();
        assert_eq!(accepted, 2);
        assert_eq!(tie_state.best.bonds, vec![1]);
        assert!(tie_result.mcs_found);
        assert_eq!(
            tie_state.degenerate.keys().cloned().collect::<Vec<_>>(),
            vec![vec![0], vec![1]]
        );

        let mut bounded_queue = McsSeedQueue {
            seeds: vec![
                make_seed(vec![0, 1, 2], vec![0, 1]),
                make_seed(vec![0, 1], vec![0]),
            ],
        };
        let mut bounded_state = McsSearchState::default();
        let mut progress_calls = 0;
        let mut progress_callback = |_: &McsProgressData, _: &McsParameters| {
            progress_calls += 1;
            Ok(true)
        };
        let mut now = || 0;
        let bounded_result = mcs_grow_seeds(
            &mut bounded_queue,
            &mut bounded_state,
            &context,
            &query,
            &[],
            &[],
            0,
            2,
            &McsParameters::default(),
            0,
            &mut now,
            None,
            None,
            Some(&mut progress_callback),
        )
        .unwrap();
        assert_eq!(bounded_state.best.bonds, vec![0, 1]);
        assert!(bounded_result.mcs_found);
        assert_eq!(progress_calls, 1);
        assert_eq!(bounded_queue.seeds.len(), 1);
        assert_eq!(bounded_queue.seeds[0].molecule_fragment.bonds, vec![0]);
    }

    #[test]
    fn q124_complete_ring_helper_precedes_acceptance_and_cancellation_retains_best() {
        let topology = q111_cycle(3);
        let rings = cosmolkit_core::fast_find_rings(&topology).unwrap();
        let coordinates = CoordinateBlock::default();
        let query = SearchTarget::new(
            &topology,
            &coordinates,
            &topology.stereo_groups,
            Some(&rings),
            None,
        );
        let make_seed = |bonds: Vec<usize>| McsSeed {
            copy_complete: true,
            growing_stage: u32::MAX,
            molecule_fragment: McsMoleculeFragment {
                atoms: vec![0, 1, 2],
                bonds,
                seed_atom_index_map: BTreeMap::new(),
            },
            excluded_bonds: vec![true; query.num_bonds()],
            ..McsSeed::default()
        };
        let params = McsParameters {
            bond_compare_parameters: McsBondCompareParameters {
                complete_rings_only: true,
                ..McsBondCompareParameters::default()
            },
            ..McsParameters::default()
        };

        let mut partial_queue = McsSeedQueue {
            seeds: vec![make_seed(vec![0, 1])],
        };
        let mut partial_state = McsSearchState::default();
        let context = McsResultContext {
            query_input: 0,
            target_inputs: Vec::new(),
        };
        let mut accept_calls = 0;
        let mut accept = |_: &McsMoleculeFragment, _: &McsParameters| {
            accept_calls += 1;
            Ok(true)
        };
        let mut now = || 0;
        let partial = mcs_grow_seeds(
            &mut partial_queue,
            &mut partial_state,
            &context,
            &query,
            &[],
            &[],
            0,
            query.num_bonds(),
            &params,
            0,
            &mut now,
            None,
            Some(&mut accept),
            None,
        )
        .unwrap();
        assert!(!partial.mcs_found);
        assert!(partial_state.best.bonds.is_empty());
        assert_eq!(accept_calls, 0);

        let mut full_queue = McsSeedQueue {
            seeds: vec![make_seed(vec![0, 1, 2])],
        };
        let mut full_state = McsSearchState::default();
        let mut canceled_progress = |data: &McsProgressData, _: &McsParameters| {
            assert_eq!((data.num_atoms, data.num_bonds), (3, 3));
            Ok(false)
        };
        let mut now = || 0;
        let full = mcs_grow_seeds(
            &mut full_queue,
            &mut full_state,
            &context,
            &query,
            &[],
            &[],
            0,
            query.num_bonds(),
            &params,
            0,
            &mut now,
            None,
            None,
            Some(&mut canceled_progress),
        )
        .unwrap();
        assert_eq!(full_state.best.bonds, vec![0, 1, 2]);
        assert!(full.mcs_found);
        assert!(full.canceled);
    }

    #[test]
    fn q131_result_query_renumbers_fragment_rows_and_collects_target_alternatives() {
        let query_topology = TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::N)),
                Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::O)),
            ],
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(2), AtomId::new(0), BondOrder::Double),
            )],
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed Q131 source topology is valid");
        let target_topology = TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(
                    AtomId::new(0),
                    AtomSpec::new(Element::F).with_chiral_tag(ChiralTag::TetrahedralCw),
                ),
                Atom::from_spec(
                    AtomId::new(1),
                    AtomSpec::new(Element::CL).with_chiral_tag(ChiralTag::TetrahedralCcw),
                ),
            ],
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                    .with_stereo(BondStereo::E),
            )],
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed Q131 target topology is valid");
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&query_topology, &coordinates);
        let target = q107_target(&target_topology, &coordinates);
        let tables = q126_tables(&query, &target, true, true);
        let result = mcs_build_result_query_graph(
            &McsMoleculeFragment {
                atoms: vec![2, 0],
                bonds: vec![0],
                seed_atom_index_map: BTreeMap::new(),
            },
            &query,
            &[target],
            &[tables],
            &McsParameters {
                atom_compare_parameters: McsAtomCompareParameters {
                    match_chiral_tag: true,
                    ..McsAtomCompareParameters::default()
                },
                bond_compare_parameters: McsBondCompareParameters {
                    match_stereo: true,
                    ..McsBondCompareParameters::default()
                },
                ..McsParameters::default()
            },
            None,
        )
        .unwrap();

        assert_eq!(result.num_atoms(), 2);
        assert_eq!(result.num_bonds(), 1);
        assert_eq!(result.atoms()[0].id(), AtomId::new(0));
        assert_eq!(result.atoms()[1].id(), AtomId::new(1));
        assert_eq!(
            result.atoms()[0].identity(),
            QueryAtomIdentity::from_atomic_number(0)
        );
        assert_eq!(
            result.atoms()[1].identity(),
            QueryAtomIdentity::from_atomic_number(0)
        );
        assert_eq!(
            result.atoms()[0].predicate(),
            &QueryNode::or(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(9)),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
            ])
        );
        assert_eq!(
            result.atoms()[1].predicate(),
            &QueryNode::or(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(17)),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            ])
        );
        assert_eq!(result.atoms()[0].chiral_tag(), ChiralTag::TetrahedralCw);
        assert_eq!(result.atoms()[1].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(result.bonds()[0].id(), BondId::new(0));
        assert_eq!(result.bonds()[0].begin(), AtomId::new(0));
        assert_eq!(result.bonds()[0].end(), AtomId::new(1));
        assert_eq!(result.bonds()[0].bond().order(), BondOrder::Unspecified);
        assert_eq!(result.bonds()[0].bond().stereo(), BondStereo::E);
        assert_eq!(
            result.bonds()[0].predicate(),
            &QueryNode::or(vec![
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
            ])
        );
    }

    #[test]
    fn q131_result_query_preserves_isotope_and_complete_ring_predicates() {
        let topology = TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C).with_isotope(13)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::C)),
            ],
            [(0, 1), (1, 2), (2, 0)]
                .into_iter()
                .enumerate()
                .map(|(index, (begin, end))| {
                    Bond::from_spec(
                        BondId::new(index),
                        BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                    )
                })
                .collect(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed Q131 isotope ring topology is valid");
        let rings = cosmolkit_core::fast_find_rings(&topology).unwrap();
        let coordinates = CoordinateBlock::default();
        let query = SearchTarget::new(
            &topology,
            &coordinates,
            &topology.stereo_groups,
            Some(&rings),
            None,
        );
        let result = mcs_build_result_query_graph(
            &McsMoleculeFragment {
                atoms: vec![0, 1, 2],
                bonds: vec![0, 1, 2],
                seed_atom_index_map: BTreeMap::new(),
            },
            &query,
            &[],
            &[],
            &McsParameters {
                atom_comparator: AtomComparator::AtomCompareIsotopes,
                bond_compare_parameters: McsBondCompareParameters {
                    match_fused_rings_strict: true,
                    ..McsBondCompareParameters::default()
                },
                ..McsParameters::default()
            },
            None,
        )
        .unwrap();

        assert_eq!(
            result.atoms()[0].predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
                QueryNode::or(vec![
                    QueryNode::predicate(AtomQueryPredicate::SmallestRingSize(3)),
                    QueryNode::not(QueryNode::predicate(AtomQueryPredicate::NumAtomRings(1))),
                ]),
            ])
        );
        assert_eq!(
            result.atoms()[1].predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::Isotope(0)),
                QueryNode::or(vec![
                    QueryNode::predicate(AtomQueryPredicate::SmallestRingSize(3)),
                    QueryNode::not(QueryNode::predicate(AtomQueryPredicate::NumAtomRings(1))),
                ]),
            ])
        );
        for bond in result.bonds() {
            assert_eq!(bond.bond().order(), BondOrder::Unspecified);
            assert_eq!(
                bond.predicate(),
                &QueryNode::and(vec![
                    QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
                    QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
                ])
            );
        }
    }

    #[test]
    fn q131_result_query_distinguishes_partial_ring_and_nonring_rows() {
        let ring_topology = q111_cycle(3);
        let ring_info = cosmolkit_core::fast_find_rings(&ring_topology).unwrap();
        let coordinates = CoordinateBlock::default();
        let ring_query = SearchTarget::new(
            &ring_topology,
            &coordinates,
            &ring_topology.stereo_groups,
            Some(&ring_info),
            None,
        );
        let params = McsParameters {
            atom_compare_parameters: McsAtomCompareParameters {
                ring_matches_ring_only: true,
                ..McsAtomCompareParameters::default()
            },
            bond_compare_parameters: McsBondCompareParameters {
                ring_matches_ring_only: true,
                ..McsBondCompareParameters::default()
            },
            ..McsParameters::default()
        };
        let partial = mcs_build_result_query_graph(
            &McsMoleculeFragment {
                atoms: vec![0, 1],
                bonds: vec![0],
                seed_atom_index_map: BTreeMap::new(),
            },
            &ring_query,
            &[],
            &[],
            &params,
            None,
        )
        .unwrap();
        assert_eq!(
            partial.atoms()[0].predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::predicate(AtomQueryPredicate::InRing),
            ])
        );
        assert_eq!(
            partial.bonds()[0].predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
                QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
            ])
        );

        let chain_topology = q113_chain(2);
        let chain_rings = cosmolkit_core::fast_find_rings(&chain_topology).unwrap();
        let chain_query = SearchTarget::new(
            &chain_topology,
            &coordinates,
            &chain_topology.stereo_groups,
            Some(&chain_rings),
            None,
        );
        let nonring = mcs_build_result_query_graph(
            &McsMoleculeFragment {
                atoms: vec![0, 1],
                bonds: vec![0],
                seed_atom_index_map: BTreeMap::new(),
            },
            &chain_query,
            &[],
            &[],
            &params,
            None,
        )
        .unwrap();
        assert_eq!(
            nonring.atoms()[0].predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::not(QueryNode::predicate(AtomQueryPredicate::InRing)),
            ])
        );
        assert_eq!(
            nonring.bonds()[0].predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
                QueryNode::not(QueryNode::predicate(BondQueryPredicate::IsInRing(true))),
            ])
        );
    }

    #[test]
    fn q132_result_smarts_uses_remapped_query_graph_and_canonical_writer() {
        let topology = TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::N)),
                Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::O)),
            ],
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(2), AtomId::new(0), BondOrder::Double),
            )],
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed Q132 source topology is valid");
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&topology, &coordinates);
        let (smarts, graph) = mcs_generate_result_smarts_and_query_graph(
            &McsMoleculeFragment {
                atoms: vec![2, 0],
                bonds: vec![0],
                seed_atom_index_map: BTreeMap::new(),
            },
            &query,
            &[],
            &[],
            &McsParameters::default(),
            None,
        )
        .unwrap();

        assert_eq!(smarts, "[#8]=[#6]");
        assert_eq!(
            smarts,
            crate::query_graph_to_smarts(&graph, &crate::SmartsWriteParams::default()).unwrap()
        );
        assert_eq!(
            graph.atoms()[0].identity(),
            QueryAtomIdentity::from_atomic_number(0)
        );
        assert_eq!(
            graph.atoms()[1].identity(),
            QueryAtomIdentity::from_atomic_number(0)
        );
        assert_eq!(graph.bonds()[0].begin(), AtomId::new(0));
        assert_eq!(graph.bonds()[0].end(), AtomId::new(1));
    }

    #[test]
    fn q133_seed_reconstruction_switches_to_new_query_row_identities() {
        let query_topology = TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::N)),
                Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::O)),
            ],
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(2), AtomId::new(0), BondOrder::Double),
            )],
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed Q133 source topology is valid");
        let target_topology = TopologyBlock::try_from_parts(
            (0..4)
                .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
                .collect(),
            [
                (0, 2, BondOrder::Single),
                (3, 1, BondOrder::Double),
                (1, 2, BondOrder::Single),
            ]
            .into_iter()
            .enumerate()
            .map(|(index, (begin, end, order))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
                )
            })
            .collect(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed Q133 target topology is valid");
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&query_topology, &coordinates);
        let target = q107_target(&target_topology, &coordinates);
        let mut atom_table = McsMatchTable {
            rows: query.num_atoms(),
            columns: target.num_atoms(),
            values: vec![false; query.num_atoms() * target.num_atoms()],
        };
        atom_table.set(2, 3, true);
        atom_table.set(0, 1, true);
        let mut bond_table = McsMatchTable {
            rows: query.num_bonds(),
            columns: target.num_bonds(),
            values: vec![false; query.num_bonds() * target.num_bonds()],
        };
        bond_table.set(0, 1, true);
        let tables = McsMatchTables {
            atoms: atom_table,
            bonds: bond_table,
        };
        let params = McsParameters {
            store_all: true,
            ..McsParameters::default()
        };
        let mut checked_mapping = Vec::new();
        let mut final_check = |mapping: &[(usize, usize)], observed: &McsParameters| {
            checked_mapping = mapping.to_vec();
            assert!(observed.store_all);
            Ok(true)
        };
        let seed = mcs_create_seed_from_mcs(
            &McsMoleculeFragment {
                atoms: vec![2, 0],
                bonds: vec![0],
                seed_atom_index_map: BTreeMap::new(),
            },
            &query,
            &target,
            &tables,
            &params,
            Some(&mut final_check),
        )
        .unwrap()
        .expect("fixed Q133 fragment matches the new query");

        assert_eq!(checked_mapping, vec![(0, 3), (1, 1)]);
        assert!(seed.store_all_degenerate_mcs);
        assert_eq!(seed.molecule_fragment.atoms, vec![3, 1]);
        assert_eq!(seed.molecule_fragment.bonds, vec![1]);
        assert_eq!(
            seed.molecule_fragment.seed_atom_index_map,
            BTreeMap::from([(1, 1), (3, 0)])
        );
        assert_eq!(seed.topology.source_atoms, vec![3, 1]);
        assert_eq!(
            seed.topology.bonds,
            vec![McsSeedTopologyBond {
                source_bond: 1,
                begin_seed_atom: 0,
                end_seed_atom: 1,
            }]
        );
        assert_eq!(seed.excluded_bonds, vec![false, true, false]);
        assert_eq!(seed.remaining_bonds, 2);
        assert_eq!(seed.remaining_atoms, 2);
    }

    #[test]
    fn q133_seed_reconstruction_returns_none_for_source_no_match() {
        let query_topology = q113_chain(2);
        let target_topology = q113_chain(2);
        let coordinates = CoordinateBlock::default();
        let query = q107_target(&query_topology, &coordinates);
        let target = q107_target(&target_topology, &coordinates);
        let tables = q126_tables(&query, &target, false, false);
        let mut final_check_calls = 0;
        let mut final_check = |_: &[(usize, usize)], _: &McsParameters| {
            final_check_calls += 1;
            Ok(true)
        };

        assert_eq!(
            mcs_create_seed_from_mcs(
                &McsMoleculeFragment {
                    atoms: vec![0, 1],
                    bonds: vec![0],
                    seed_atom_index_map: BTreeMap::new(),
                },
                &query,
                &target,
                &tables,
                &McsParameters::default(),
                Some(&mut final_check),
            )
            .unwrap(),
            None
        );
        assert_eq!(final_check_calls, 0);
    }
}
