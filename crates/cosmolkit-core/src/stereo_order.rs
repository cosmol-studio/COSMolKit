//! Detached tetrahedral ligand ordering and mapping primitives.
//!
//! This module owns parity-preserving order conversion. It does not perceive
//! stereo, mutate topology, or accept a live molecule/runtime capability.

use cosmolkit_model::{
    Atom, AtomId, Bond, BondId, BondValueError, MappingValidationError, TopologyBlock,
    TopologyMapping, TopologyValidationError,
};
use cosmolkit_types::{BondOrder, BondStereo, ChiralTag};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TetrahedralLigand {
    Bond(BondId),
    Implicit,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TetrahedralRemap {
    pub center: AtomId,
    pub ligands: Vec<TetrahedralLigand>,
    pub tag: ChiralTag,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum StereoOrderError {
    #[error("source chiral permutation read failed: {0}")]
    ChiralPermutationRead(#[from] crate::PropertyUIntReadError),
    #[error("source chiral permutation property write failed: {0}")]
    ChiralPermutationWrite(#[from] cosmolkit_model::AtomPropertyError),
    #[error(
        "atom {atom} has conflicting typed ({typed}) and raw ({raw}) chiral permutation projections"
    )]
    ConflictingChiralPermutation { atom: AtomId, typed: u32, raw: u32 },

    #[error("reference order has {reference} entries but probe order has {probe}")]
    PermutationLength { reference: usize, probe: usize },
    #[error(
        "probe order does not contain the value required at reference position {reference_position}"
    )]
    MissingProbeValue { reference_position: usize },
    #[error("bond index {bond_index} exceeds the source unsigned32 index range")]
    BondIndexSourceWidth { bond_index: usize },
    #[error("invalid topology: {0}")]
    InvalidTopology(TopologyValidationError),
    #[error("invalid topology mapping: {0}")]
    InvalidMapping(MappingValidationError),
    #[error("stereo center {center} is out of range for {atom_count} atoms")]
    CenterOutOfRange { center: AtomId, atom_count: usize },
    #[error("bond {bond} is out of range for {bond_count} bonds")]
    BondOutOfRange { bond: BondId, bond_count: usize },
    #[error("bond {bond} is not incident to stereo center {center}")]
    BondNotIncident { bond: BondId, center: AtomId },
    #[error("tetrahedral order has {actual} ligands; maximum is {maximum}")]
    TooManyLigands { actual: usize, maximum: usize },
    #[error("bond {bond} is duplicated at tetrahedral ligand position {position}")]
    DuplicateBondLigand { position: usize, bond: BondId },
    #[error("implicit ligand is duplicated at positions {first} and {second}")]
    MultipleImplicitLigands { first: usize, second: usize },
    #[error("chiral tag {tag:?} is not a supported explicit tetrahedral CW/CCW tag")]
    UnsupportedTetrahedralTag { tag: ChiralTag },
    #[error("stereo center {center} was removed by the topology mapping")]
    RemovedCenter { center: AtomId },
    #[error("tetrahedral ligand bond {bond} was removed by the topology mapping")]
    RemovedLigand { bond: BondId },
    #[error("implicit replacement bond {bond} is out of range for {old_bond_count} source bonds")]
    ImplicitReplacementOutOfRange { bond: BondId, old_bond_count: usize },
    #[error("implicit replacement bond {bond} was retained as {mapped}")]
    ImplicitReplacementRetained { bond: BondId, mapped: BondId },
    #[error("implicit replacement bond {bond} is not present in the source ligand order")]
    ImplicitReplacementNotInSource { bond: BondId },
}

/// Count the first-match swaps used by RDKit to convert `probe` to `reference`.
pub fn count_swaps_to_interconvert<T: Copy + Eq>(
    reference: &[T],
    probe: &[T],
) -> Result<usize, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION RDGeneral::countSwapsToInterconvert
    // RDKit❗✔️: template <class T>
    // RDKit❗✔️: unsigned int countSwapsToInterconvert(const T &ref, T probe) {
    // RDKit❗✔️:   PRECONDITION(ref.size() == probe.size(), "size mismatch");
    // RDKit❗✔️:   typename T::const_iterator refIt = ref.begin();
    // RDKit❗✔️:   typename T::iterator probeIt = probe.begin();
    // RDKit❗✔️:   typename T::iterator probeIt2;
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int nSwaps = 0;
    // RDKit❗✔️:   while (refIt != ref.end()) {
    // RDKit❗✔️:     if ((*probeIt) != (*refIt)) {
    // RDKit❗✔️:       bool foundIt = false;
    // RDKit❗✔️:       probeIt2 = probeIt;
    // RDKit❗✔️:       while ((*probeIt2) != (*refIt) && probeIt2 != probe.end()) {
    // RDKit❗✔️:         ++probeIt2;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (probeIt2 != probe.end()) {
    // RDKit❗✔️:         foundIt = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       CHECK_INVARIANT(foundIt, "could not find probe element");
    // RDKit❗✔️:
    // RDKit❗✔️:       std::swap(*probeIt, *probeIt2);
    // RDKit❗✔️:       nSwaps++;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     ++probeIt;
    // RDKit❗✔️:     ++refIt;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return nSwaps;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDGeneral::countSwapsToInterconvert
    // Source nSwaps is unsigned32: preserve defined wrap at every increment,
    // then widen its bits to the existing public usize return carrier.
    // The first matching value in the remaining probe suffix is exchanged;
    // duplicate values retain native encounter order. One detached probe copy
    // matches native pass-by-value storage, with O(N^2) worst-case comparisons
    // and O(N) scratch; no graph clone, sorting or alternate parity algorithm.
    // Missing values become the source CHECK_INVARIANT's structural error.
    // Native tests the dereference before its end guard, which is undefined
    // for a missing value. Rust never dereferences past the slice; behavioral
    // status stays ❗ for that source-language edge, not an invented fallback.
    if reference.len() != probe.len() {
        return Err(StereoOrderError::PermutationLength {
            reference: reference.len(),
            probe: probe.len(),
        });
    }

    let mut probe = probe.to_vec();
    let mut swaps = 0_u32;
    for (reference_position, expected) in reference.iter().enumerate() {
        if probe[reference_position] == *expected {
            continue;
        }
        let Some(offset) = probe[reference_position..]
            .iter()
            .position(|candidate| candidate == expected)
        else {
            return Err(StereoOrderError::MissingProbeValue { reference_position });
        };
        probe.swap(reference_position, reference_position + offset);
        swaps = swaps.wrapping_add(1);
    }
    Ok(swaps as usize)
}

/// Count source atom perturbations from its actual incident-bond encounter order.
pub fn atom_perturbation_order(
    probe: &[i32],
    incident_bond_indices: impl IntoIterator<Item = usize>,
) -> Result<i32, StereoOrderError> {
    // BEGIN COMPLETE PINNED SF304 Atom::getPerturbationOrder
    // RDKit❗🔝: int Atom::getPerturbationOrder(const INT_LIST &probe) const {
    // RDKit❗🔝:   INT_LIST ref;
    // RDKit❗🔝:   for (const auto bnd : getOwningMol().atomBonds(this)) {
    // RDKit❗🔝:     ref.push_back(bnd->getIdx());
    // RDKit❗🔝:   }
    // RDKit❗🔝:   return static_cast<int>(countSwapsToInterconvert(probe, ref));
    // RDKit❗🔝: }
    // END COMPLETE PINNED SF304 Atom::getPerturbationOrder
    // The detached iterator supplies the actual owning graph's incident order;
    // consume it once, without sorting, filtering zero/dative bonds or imposing
    // tetrahedral degree limits. Native unsigned32 indices become INT_LIST's
    // signed32 elements, including indices with the high bit set. Wider project
    // indices fail structurally at this transport boundary, never as unsupported.
    // Reuse the sole swap algorithm with native argument orientation. Its
    // unsigned32 result is converted to the source signed32 return bits.
    // Counter behavior retains ❗ for the missing-value native undefined read.
    // Vec storage removes linked-list node allocation and pointer chasing while
    // preserving encounter order: O(D) construction and O(D^2) swap scanning,
    // O(D) scratch, fewer allocations than two source INT_LIST instances.
    let reference = incident_bond_indices
        .into_iter()
        .map(|bond_index| {
            u32::try_from(bond_index)
                .map(|index| index as i32)
                .map_err(|_| StereoOrderError::BondIndexSourceWidth { bond_index })
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(count_swaps_to_interconvert(probe, &reference)? as u32 as i32)
}

/// Invert a source-supported atom stereochemical tag and its permutation state.
pub fn invert_atom_chirality(atom: &mut Atom) -> Result<bool, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION Atom::invertChirality tables and body
    // RDKit❗❌: static const unsigned char octahedral_invert[31] = {
    // RDKit❗❌:     0,   //  0 -> 0
    // RDKit❗❌:     2,   //  1 -> 2
    // RDKit❗❌:     1,   //  2 -> 1
    // RDKit❗❌:     16,  //  3 -> 16
    // RDKit❗❌:     14,  //  4 -> 14
    // RDKit❗❌:     15,  //  5 -> 15
    // RDKit❗❌:     18,  //  6 -> 18
    // RDKit❗❌:     17,  //  7 -> 17
    // RDKit❗❌:     10,  //  8 -> 10
    // RDKit❗❌:     11,  //  9 -> 11
    // RDKit❗❌:     8,   // 10 -> 8
    // RDKit❗❌:     9,   // 11 -> 9
    // RDKit❗❌:     13,  // 12 -> 13
    // RDKit❗❌:     12,  // 13 -> 12
    // RDKit❗❌:     4,   // 14 -> 4
    // RDKit❗❌:     5,   // 15 -> 5
    // RDKit❗❌:     3,   // 16 -> 3
    // RDKit❗❌:     7,   // 17 -> 7
    // RDKit❗❌:     6,   // 18 -> 6
    // RDKit❗❌:     24,  // 19 -> 24
    // RDKit❗❌:     23,  // 20 -> 23
    // RDKit❗❌:     22,  // 21 -> 22
    // RDKit❗❌:     21,  // 22 -> 21
    // RDKit❗❌:     20,  // 23 -> 20
    // RDKit❗❌:     19,  // 24 -> 19
    // RDKit❗❌:     30,  // 25 -> 30
    // RDKit❗❌:     29,  // 26 -> 29
    // RDKit❗❌:     28,  // 27 -> 28
    // RDKit❗❌:     27,  // 28 -> 27
    // RDKit❗❌:     26,  // 29 -> 26
    // RDKit❗❌:     25   // 30 -> 25
    // RDKit❗❌: };
    // RDKit❗❌:
    // RDKit❗❌: static const unsigned char trigonalbipyramidal_invert[21] = {
    // RDKit❗❌:     0,   //  0 -> 0
    // RDKit❗❌:     2,   //  1 -> 2
    // RDKit❗❌:     1,   //  2 -> 1
    // RDKit❗❌:     4,   //  3 -> 4
    // RDKit❗❌:     3,   //  4 -> 3
    // RDKit❗❌:     6,   //  5 -> 6
    // RDKit❗❌:     5,   //  6 -> 5
    // RDKit❗❌:     8,   //  7 -> 8
    // RDKit❗❌:     7,   //  8 -> 7
    // RDKit❗❌:     11,  //  9 -> 11
    // RDKit❗❌:     12,  // 10 -> 12
    // RDKit❗❌:     9,   // 11 -> 9
    // RDKit❗❌:     10,  // 12 -> 10
    // RDKit❗❌:     14,  // 13 -> 14
    // RDKit❗❌:     13,  // 14 -> 13
    // RDKit❗❌:     20,  // 15 -> 20
    // RDKit❗❌:     19,  // 16 -> 19
    // RDKit❗❌:     18,  // 17 -> 28
    // RDKit❗❌:     17,  // 18 -> 17
    // RDKit❗❌:     16,  // 19 -> 16
    // RDKit❗❌:     15   // 20 -> 15
    // RDKit❗❌: };
    // RDKit❗❌: bool Atom::invertChirality() {
    // RDKit❗❌:   unsigned int perm;
    // RDKit❗❌:   switch (getChiralTag()) {
    // RDKit❗❌:     case CHI_TETRAHEDRAL_CW:
    // RDKit❗❌:       setChiralTag(CHI_TETRAHEDRAL_CCW);
    // RDKit❗❌:       return true;
    // RDKit❗❌:     case CHI_TETRAHEDRAL_CCW:
    // RDKit❗❌:       setChiralTag(CHI_TETRAHEDRAL_CW);
    // RDKit❗❌:       return true;
    // RDKit❗❌:     case CHI_TETRAHEDRAL:
    // RDKit❗❌:       if (getPropIfPresent(common_properties::_chiralPermutation, perm)) {
    // RDKit❗❌:         if (perm == 1) {
    // RDKit❗❌:           perm = 2;
    // RDKit❗❌:         } else if (perm == 2) {
    // RDKit❗❌:           perm = 1;
    // RDKit❗❌:         } else {
    // RDKit❗❌:           perm = 0;
    // RDKit❗❌:         }
    // RDKit❗❌:         setProp(common_properties::_chiralPermutation, perm);
    // RDKit❗❌:         return perm != 0;
    // RDKit❗❌:       }
    // RDKit❗❌:       break;
    // RDKit❗❌:     case CHI_TRIGONALBIPYRAMIDAL:
    // RDKit❗❌:       if (getPropIfPresent(common_properties::_chiralPermutation, perm)) {
    // RDKit❗❌:         perm = (perm <= 20) ? trigonalbipyramidal_invert[perm] : 0;
    // RDKit❗❌:         setProp(common_properties::_chiralPermutation, perm);
    // RDKit❗❌:         return perm != 0;
    // RDKit❗❌:       }
    // RDKit❗❌:       break;
    // RDKit❗❌:     case CHI_OCTAHEDRAL:
    // RDKit❗❌:       if (getPropIfPresent(common_properties::_chiralPermutation, perm)) {
    // RDKit❗❌:         perm = (perm <= 30) ? octahedral_invert[perm] : 0;
    // RDKit❗❌:         setProp(common_properties::_chiralPermutation, perm);
    // RDKit❗❌:         return perm != 0;
    // RDKit❗❌:       }
    // RDKit❗❌:       break;
    // RDKit❗❌:     default:
    // RDKit❗❌:       break;
    // RDKit❗❌:   }
    // RDKit❗❌:   return false;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION Atom::invertChirality tables and body
    const OCTAHEDRAL_INVERT: [u32; 31] = [
        0, 2, 1, 16, 14, 15, 18, 17, 10, 11, 8, 9, 13, 12, 4, 5, 3, 7, 6, 24, 23, 22, 21, 20, 19,
        30, 29, 28, 27, 26, 25,
    ];
    const TRIGONAL_BIPYRAMIDAL_INVERT: [u32; 21] = [
        0, 2, 1, 4, 3, 6, 5, 8, 7, 11, 12, 9, 10, 14, 13, 20, 19, 18, 17, 16, 15,
    ];
    // Complete native branch/table behavior. Canonical source UInt reads are
    // performed only for tags whose source switch reads the property. Model's
    // explicit typed source fact and raw dictionary are two existing transport
    // projections; contradictory duplicate facts are a structural input error,
    // never a guessed precedence or absent-value default. Native ownership,
    // dictionary insertion identity for typed-only facts and dual projection
    // synchronization remain model boundary differences, with ❗ behavior.
    // Constant dispatch/table access; same canonical property conversion cost,
    // no graph clone, neighbor scan, sort, or speculative property reads.
    // Known cost gap: PropertyText owns Vec bytes, so raw-key overwrite makes
    // an extra key allocation compared with native borrowed string_view and
    // existing Dict Pair overwrite, despite faster tree lookup for many keys.
    match atom.chiral_tag() {
        ChiralTag::TetrahedralCw => {
            atom.set_chiral_tag(ChiralTag::TetrahedralCcw);
            Ok(true)
        }
        ChiralTag::TetrahedralCcw => {
            atom.set_chiral_tag(ChiralTag::TetrahedralCw);
            Ok(true)
        }
        ChiralTag::Tetrahedral => invert_chiral_permutation(atom, &[0, 2, 1]),
        ChiralTag::TrigonalBipyramidal => {
            invert_chiral_permutation(atom, &TRIGONAL_BIPYRAMIDAL_INVERT)
        }
        ChiralTag::Octahedral => invert_chiral_permutation(atom, &OCTAHEDRAL_INVERT),
        _ => Ok(false),
    }
}

fn invert_chiral_permutation(atom: &mut Atom, table: &[u32]) -> Result<bool, StereoOrderError> {
    // RDKit❗❌: if (getPropIfPresent(common_properties::_chiralPermutation, perm)) {
    // RDKit❗❌:   perm = (perm <= 20) ? trigonalbipyramidal_invert[perm] : 0;
    // RDKit❗❌:   setProp(common_properties::_chiralPermutation, perm);
    // RDKit❗❌:   return perm != 0;
    // RDKit❗❌: }
    // RDKit❗❌: template <typename T>
    // RDKit❗❌:   bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit❗❌:     return d_props.getValIfPresent(key, res);
    // RDKit❗❌:   }
    // Canonical reader implements native UInt numeric/tag/right-trimmed lexical
    // conversion and propagates its actual errors before any atom mutation.
    let raw = atom.prop("_chiralPermutation");
    let has_raw = raw.is_some();
    let permutation = match raw {
        Some(value) => {
            let value = crate::property_value_to_uint(value)?;
            if let Some(typed) = atom.chiral_permutation()
                && typed != value
            {
                return Err(StereoOrderError::ConflictingChiralPermutation {
                    atom: atom.id(),
                    typed,
                    raw: value,
                });
            }
            Some(value)
        }
        None => atom.chiral_permutation(),
    };
    let Some(permutation) = permutation else {
        return Ok(false);
    };
    let inverted = if (permutation as usize) < table.len() {
        table[permutation as usize]
    } else {
        0
    };
    if has_raw {
        atom.set_prop(
            "_chiralPermutation",
            cosmolkit_model::PropertyValue::UInt(inverted),
        )?;
    }
    // Do not manufacture a typed fact for raw-only source state: otherwise a
    // later computed-property clear would leave a ghost permutation behind.
    if !has_raw || atom.chiral_permutation().is_some() {
        atom.set_chiral_permutation(Some(inverted));
    }
    Ok(inverted != 0)
}

/// Invert atrop stereochemistry on a bond while preserving all other bond state.
pub fn invert_bond_chirality(bond: &mut Bond) -> Result<bool, BondValueError> {
    // BEGIN RDKIT CPP FUNCTION Bond::invertChirality
    // RDKit✔️✔️: bool Bond::invertChirality() {
    // RDKit✔️✔️:   switch (getStereo()) {
    // RDKit✔️✔️:     case STEREOATROPCW:
    // RDKit✔️✔️:       setStereo(STEREOATROPCCW);
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     case STEREOATROPCCW:
    // RDKit✔️✔️:       setStereo(STEREOATROPCW);
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Bond::invertChirality
    // Behavior review: only the two source-supported atrop states change; every
    // other stereo state and all unrelated bond fields remain untouched.
    // Complexity review: one enum match and one typed setter call, with no
    // allocation or topology scan, preserves the source constant-time path.
    match bond.stereo() {
        BondStereo::AtropCw => {
            bond.set_stereo(BondStereo::AtropCcw)?;
            Ok(true)
        }
        BondStereo::AtropCcw => {
            bond.set_stereo(BondStereo::AtropCw)?;
            Ok(true)
        }
        _ => Ok(false),
    }
}

pub fn invert_tetrahedral_tag(tag: ChiralTag) -> Result<ChiralTag, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION Atom::invertChirality tetrahedral CW/CCW branches
    // RDKit✔️✔️:     case CHI_TETRAHEDRAL_CW:
    // RDKit✔️✔️:       setChiralTag(CHI_TETRAHEDRAL_CCW);
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     case CHI_TETRAHEDRAL_CCW:
    // RDKit✔️✔️:       setChiralTag(CHI_TETRAHEDRAL_CW);
    // RDKit✔️✔️:       return true;
    // END RDKIT CPP FUNCTION Atom::invertChirality tetrahedral CW/CCW branches
    match tag {
        ChiralTag::TetrahedralCw => Ok(ChiralTag::TetrahedralCcw),
        ChiralTag::TetrahedralCcw => Ok(ChiralTag::TetrahedralCw),
        tag => Err(StereoOrderError::UnsupportedTetrahedralTag { tag }),
    }
}

pub fn bond_affects_atom_chirality(bond: &Bond, center: AtomId) -> Result<bool, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::detail::bondAffectsAtomChirality
    // RDKit✔️✔️: bool bondAffectsAtomChirality(const Bond *bond, const Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond pointer");
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom pointer");
    // RDKit✔️✔️:   if (bond->getBondType() == Bond::BondType::UNSPECIFIED ||
    // RDKit✔️✔️:       bond->getBondType() == Bond::BondType::ZERO ||
    // RDKit✔️✔️:       (bond->getBondType() == Bond::BondType::DATIVE &&
    // RDKit✔️✔️:        bond->getBeginAtomIdx() == atom->getIdx())) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::detail::bondAffectsAtomChirality
    if bond.begin() != center && bond.end() != center {
        return Err(StereoOrderError::BondNotIncident {
            bond: bond.id(),
            center,
        });
    }
    Ok(
        !matches!(bond.order(), BondOrder::Unspecified | BondOrder::Zero)
            && !(bond.order() == BondOrder::Dative && bond.begin() == center),
    )
}

pub fn atom_nonzero_degree(
    topology: &TopologyBlock,
    center: AtomId,
) -> Result<usize, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::detail::getAtomNonzeroDegree
    // RDKit✔️❌: unsigned int getAtomNonzeroDegree(const Atom *atom) {
    // RDKit✔️❌:   PRECONDITION(atom, "bad pointer");
    // RDKit✔️❌:   PRECONDITION(atom->hasOwningMol(), "no owning molecule");
    // RDKit✔️❌:   unsigned int res = 0;
    // RDKit✔️❌:   for (auto bond : atom->getOwningMol().atomBonds(atom)) {
    // RDKit✔️❌:     if (!bondAffectsAtomChirality(bond, atom)) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ++res;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION Chirality::detail::getAtomNonzeroDegree
    // The detached boundary validates the full topology before the source
    // degree scan, changing O(degree) to O(V+E); behavior is preserved but the
    // second marker is intentionally negative.
    topology
        .validate()
        .map_err(StereoOrderError::InvalidTopology)?;
    if center.index() >= topology.atoms.len() {
        return Err(StereoOrderError::CenterOutOfRange {
            center,
            atom_count: topology.atoms.len(),
        });
    }
    let mut degree = 0;
    for neighbor in topology.adjacency.neighbors_of(center.index()) {
        let bond =
            topology
                .bonds
                .get(neighbor.bond.index())
                .ok_or(StereoOrderError::BondOutOfRange {
                    bond: neighbor.bond,
                    bond_count: topology.bonds.len(),
                })?;
        if bond_affects_atom_chirality(bond, center)? {
            degree += 1;
        }
    }
    Ok(degree)
}

pub fn incident_tetrahedral_bond_order(
    topology: &TopologyBlock,
    center: AtomId,
) -> Result<Vec<TetrahedralLigand>, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION Atom::getPerturbationOrder incident storage order
    // RDKit✔️❌:   INT_LIST ref;
    // RDKit✔️❌:   for (const auto bnd : getOwningMol().atomBonds(this)) {
    // RDKit✔️❌:     ref.push_back(bnd->getIdx());
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return static_cast<int>(countSwapsToInterconvert(probe, ref));
    // END RDKIT CPP FUNCTION Atom::getPerturbationOrder incident storage order
    // The detached topology validation adds O(V+E) work before retaining the
    // source incident-bond order, so the second marker is intentionally
    // negative. Non-chirality-affecting bonds are then removed by the separate
    // source-backed predicate above.
    topology
        .validate()
        .map_err(StereoOrderError::InvalidTopology)?;
    if center.index() >= topology.atoms.len() {
        return Err(StereoOrderError::CenterOutOfRange {
            center,
            atom_count: topology.atoms.len(),
        });
    }

    let mut ligands = Vec::new();
    for neighbor in topology.adjacency.neighbors_of(center.index()) {
        let bond =
            topology
                .bonds
                .get(neighbor.bond.index())
                .ok_or(StereoOrderError::BondOutOfRange {
                    bond: neighbor.bond,
                    bond_count: topology.bonds.len(),
                })?;
        if bond_affects_atom_chirality(bond, center)? {
            ligands.push(TetrahedralLigand::Bond(bond.id()));
        }
    }
    validate_ligand_order(&ligands)?;
    Ok(ligands)
}

pub fn tetrahedral_tag_after_order_change(
    tag: ChiralTag,
    reference: &[TetrahedralLigand],
    probe: &[TetrahedralLigand],
) -> Result<ChiralTag, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION FindStereo tetrahedral controlling-order normalization
    // RDKit✔️✔️:       unsigned nSwaps =
    // RDKit✔️✔️:           countSwapsToInterconvert(origNbrOrder, sinfo.controllingAtoms);
    // RDKit✔️✔️:       if (nSwaps % 2) {
    // RDKit✔️✔️:         stereo = (stereo == Atom::ChiralType::CHI_TETRAHEDRAL_CCW
    // RDKit✔️✔️:                       ? Atom::ChiralType::CHI_TETRAHEDRAL_CW
    // RDKit✔️✔️:                       : Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️✔️:       }
    // END RDKIT CPP FUNCTION FindStereo tetrahedral controlling-order normalization
    validate_tetrahedral_tag(tag)?;
    validate_ligand_order(reference)?;
    validate_ligand_order(probe)?;
    let swaps = count_swaps_to_interconvert(reference, probe)?;
    if swaps % 2 == 1 {
        invert_tetrahedral_tag(tag)
    } else {
        Ok(tag)
    }
}

#[allow(clippy::too_many_arguments)]
pub fn remap_tetrahedral_center(
    source_center: AtomId,
    source_ligands: &[TetrahedralLigand],
    source_tag: ChiralTag,
    mapping: &TopologyMapping,
    old_atom_count: usize,
    new_atom_count: usize,
    old_bond_count: usize,
    new_bond_count: usize,
    removed_bond_as_implicit: Option<BondId>,
    target_ligands: &[TetrahedralLigand],
) -> Result<TetrahedralRemap, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION AddHs remove-H tetrahedral order adjustment
    // RDKit✔️❌:       INT_LIST neighborIndices;
    // RDKit✔️❌:       for (const auto &nbnd : mol.atomBonds(heavyAtom)) {
    // RDKit✔️❌:         if (nbnd->getIdx() != bond->getIdx()) {
    // RDKit✔️❌:           neighborIndices.push_back(nbnd->getIdx());
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       neighborIndices.push_back(bond->getIdx());
    // RDKit✔️❌:
    // RDKit✔️❌:       int nSwaps = heavyAtom->getPerturbationOrder(neighborIndices);
    // RDKit✔️❌:       if (nSwaps % 2) {
    // RDKit✔️❌:         heavyAtom->invertChirality();
    // RDKit✔️❌:       }
    // END RDKIT CPP FUNCTION AddHs remove-H tetrahedral order adjustment
    // The complete detached mapping validation is O(V+E), while the source
    // adjustment only scans one center's incident bonds. This intentional
    // fail-closed boundary accounts for the negative performance marker.
    mapping
        .validate_for_counts(
            old_atom_count,
            new_atom_count,
            old_bond_count,
            new_bond_count,
        )
        .map_err(StereoOrderError::InvalidMapping)?;
    validate_tetrahedral_tag(source_tag)?;
    validate_ligand_order(source_ligands)?;
    validate_ligand_order(target_ligands)?;
    validate_bond_ranges(source_ligands, old_bond_count)?;
    validate_bond_ranges(target_ligands, new_bond_count)?;

    if source_center.index() >= old_atom_count {
        return Err(StereoOrderError::CenterOutOfRange {
            center: source_center,
            atom_count: old_atom_count,
        });
    }
    let target_center = mapping.atoms().old_to_new()[source_center.index()].ok_or(
        StereoOrderError::RemovedCenter {
            center: source_center,
        },
    )?;

    if let Some(replacement) = removed_bond_as_implicit {
        if replacement.index() >= old_bond_count {
            return Err(StereoOrderError::ImplicitReplacementOutOfRange {
                bond: replacement,
                old_bond_count,
            });
        }
        if !source_ligands.contains(&TetrahedralLigand::Bond(replacement)) {
            return Err(StereoOrderError::ImplicitReplacementNotInSource { bond: replacement });
        }
        if let Some(mapped) = mapping.bonds().old_to_new()[replacement.index()] {
            return Err(StereoOrderError::ImplicitReplacementRetained {
                bond: replacement,
                mapped,
            });
        }
    }

    let mut mapped_ligands = Vec::with_capacity(source_ligands.len());
    for ligand in source_ligands {
        match *ligand {
            TetrahedralLigand::Implicit => mapped_ligands.push(TetrahedralLigand::Implicit),
            TetrahedralLigand::Bond(bond) => match mapping.bonds().old_to_new()[bond.index()] {
                Some(mapped) => mapped_ligands.push(TetrahedralLigand::Bond(mapped)),
                None if removed_bond_as_implicit == Some(bond) => {
                    mapped_ligands.push(TetrahedralLigand::Implicit);
                }
                None => return Err(StereoOrderError::RemovedLigand { bond }),
            },
        }
    }
    validate_ligand_order(&mapped_ligands)?;

    let tag = tetrahedral_tag_after_order_change(source_tag, &mapped_ligands, target_ligands)?;
    Ok(TetrahedralRemap {
        center: target_center,
        ligands: target_ligands.to_vec(),
        tag,
    })
}

fn validate_tetrahedral_tag(tag: ChiralTag) -> Result<(), StereoOrderError> {
    if matches!(tag, ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw) {
        Ok(())
    } else {
        Err(StereoOrderError::UnsupportedTetrahedralTag { tag })
    }
}

fn validate_ligand_order(ligands: &[TetrahedralLigand]) -> Result<(), StereoOrderError> {
    const MAXIMUM: usize = 4;
    if ligands.len() > MAXIMUM {
        return Err(StereoOrderError::TooManyLigands {
            actual: ligands.len(),
            maximum: MAXIMUM,
        });
    }
    let mut implicit_position = None;
    for (position, ligand) in ligands.iter().copied().enumerate() {
        match ligand {
            TetrahedralLigand::Implicit => {
                if let Some(first) = implicit_position {
                    return Err(StereoOrderError::MultipleImplicitLigands {
                        first,
                        second: position,
                    });
                }
                implicit_position = Some(position);
            }
            TetrahedralLigand::Bond(bond) => {
                if ligands[..position].contains(&TetrahedralLigand::Bond(bond)) {
                    return Err(StereoOrderError::DuplicateBondLigand { position, bond });
                }
            }
        }
    }
    Ok(())
}

fn validate_bond_ranges(
    ligands: &[TetrahedralLigand],
    bond_count: usize,
) -> Result<(), StereoOrderError> {
    for ligand in ligands {
        if let TetrahedralLigand::Bond(bond) = ligand
            && bond.index() >= bond_count
        {
            return Err(StereoOrderError::BondOutOfRange {
                bond: *bond,
                bond_count,
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod stereo_inversion_tests {
    use super::{invert_atom_chirality, invert_bond_chirality};
    use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec};
    use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Element};

    fn atom(tag: ChiralTag, permutation: Option<u32>) -> Atom {
        let spec = AtomSpec::new(Element::C)
            .with_formal_charge(-1)
            .with_explicit_hydrogens(2)
            .with_unknown_stereo(true)
            .with_chiral_tag(tag);
        let mut atom = Atom::from_spec(AtomId::new(0), spec);
        atom.set_chiral_permutation(permutation);
        atom
    }

    fn assert_transition(
        tag: ChiralTag,
        permutation: Option<u32>,
        expected_tag: ChiralTag,
        expected_permutation: Option<u32>,
        expected_changed: bool,
    ) {
        let mut actual = atom(tag, permutation);
        let mut expected = actual.clone();
        expected.set_chiral_tag(expected_tag);
        expected.set_chiral_permutation(expected_permutation);

        assert_eq!(
            invert_atom_chirality(&mut actual).unwrap(),
            expected_changed
        );
        assert_eq!(actual, expected);
    }

    fn bond(stereo: BondStereo) -> Bond {
        let spec = BondSpec::new(AtomId::new(3), AtomId::new(7), BondOrder::Double)
            .with_aromatic(true)
            .with_conjugated(true)
            .with_direction(BondDirection::EndDownRight)
            .with_stereo(stereo)
            .with_stereo_atoms(AtomId::new(2), AtomId::new(8))
            .with_unknown_stereo(true)
            .with_prop("source_prop", "kept".to_owned())
            .expect("valid ordinary bond property")
            .with_computed_prop("computed_prop", "kept".to_owned())
            .expect("valid computed bond property");
        Bond::from_spec(BondId::new(5), spec)
    }

    #[test]
    fn stereo_inversion_atom_tetrahedral_tags_and_permutation_edges() {
        assert_transition(
            ChiralTag::TetrahedralCw,
            None,
            ChiralTag::TetrahedralCcw,
            None,
            true,
        );
        assert_transition(
            ChiralTag::TetrahedralCcw,
            Some(17),
            ChiralTag::TetrahedralCw,
            Some(17),
            true,
        );
        assert_transition(
            ChiralTag::Tetrahedral,
            None,
            ChiralTag::Tetrahedral,
            None,
            false,
        );
        assert_transition(
            ChiralTag::Tetrahedral,
            Some(0),
            ChiralTag::Tetrahedral,
            Some(0),
            false,
        );
        assert_transition(
            ChiralTag::Tetrahedral,
            Some(1),
            ChiralTag::Tetrahedral,
            Some(2),
            true,
        );
        assert_transition(
            ChiralTag::Tetrahedral,
            Some(2),
            ChiralTag::Tetrahedral,
            Some(1),
            true,
        );
        assert_transition(
            ChiralTag::Tetrahedral,
            Some(3),
            ChiralTag::Tetrahedral,
            Some(0),
            false,
        );
        assert_transition(
            ChiralTag::Tetrahedral,
            Some(u32::MAX),
            ChiralTag::Tetrahedral,
            Some(0),
            false,
        );
    }

    #[test]
    fn stereo_inversion_atom_trigonal_bipyramidal_table_rows() {
        const EXPECTED: [u32; 21] = [
            0, 2, 1, 4, 3, 6, 5, 8, 7, 11, 12, 9, 10, 14, 13, 20, 19, 18, 17, 16, 15,
        ];

        for (permutation, expected) in EXPECTED.into_iter().enumerate() {
            assert_transition(
                ChiralTag::TrigonalBipyramidal,
                Some(permutation as u32),
                ChiralTag::TrigonalBipyramidal,
                Some(expected),
                expected != 0,
            );
        }
        for out_of_range in [21, u32::MAX] {
            assert_transition(
                ChiralTag::TrigonalBipyramidal,
                Some(out_of_range),
                ChiralTag::TrigonalBipyramidal,
                Some(0),
                false,
            );
        }
    }

    #[test]
    fn stereo_inversion_atom_octahedral_table_rows() {
        const EXPECTED: [u32; 31] = [
            0, 2, 1, 16, 14, 15, 18, 17, 10, 11, 8, 9, 13, 12, 4, 5, 3, 7, 6, 24, 23, 22, 21, 20,
            19, 30, 29, 28, 27, 26, 25,
        ];

        for (permutation, expected) in EXPECTED.into_iter().enumerate() {
            assert_transition(
                ChiralTag::Octahedral,
                Some(permutation as u32),
                ChiralTag::Octahedral,
                Some(expected),
                expected != 0,
            );
        }
        for out_of_range in [31, u32::MAX] {
            assert_transition(
                ChiralTag::Octahedral,
                Some(out_of_range),
                ChiralTag::Octahedral,
                Some(0),
                false,
            );
        }
    }

    #[test]
    fn stereo_inversion_atom_unsupported_tags_preserve_all_state() {
        for tag in [
            ChiralTag::Unspecified,
            ChiralTag::Other,
            ChiralTag::Allene,
            ChiralTag::SquarePlanar,
        ] {
            assert_transition(tag, Some(19), tag, Some(19), false);
        }
    }

    #[test]
    fn stereo_inversion_bond_atrop_states_and_repeated_inversion() {
        for (initial, inverted) in [
            (BondStereo::AtropCw, BondStereo::AtropCcw),
            (BondStereo::AtropCcw, BondStereo::AtropCw),
        ] {
            let mut actual = bond(initial);
            let original = actual.clone();
            let mut expected = original.clone();
            expected
                .set_stereo(inverted)
                .expect("atrop state is valid without reference atoms");

            assert_eq!(invert_bond_chirality(&mut actual), Ok(true));
            assert_eq!(actual, expected);
            assert_eq!(invert_bond_chirality(&mut actual), Ok(true));
            assert_eq!(actual, original);
        }
    }

    #[test]
    fn stereo_inversion_bond_other_stereo_states_preserve_all_state() {
        for stereo in [
            BondStereo::None,
            BondStereo::Any,
            BondStereo::Z,
            BondStereo::E,
            BondStereo::Cis,
            BondStereo::Trans,
        ] {
            let mut actual = bond(stereo);
            let original = actual.clone();

            assert_eq!(invert_bond_chirality(&mut actual), Ok(false));
            assert_eq!(actual, original);
        }
    }
}

#[cfg(test)]
mod source_atom_invert_chirality_property_complete_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue};
    use cosmolkit_types::Element;
    fn atom(tag: ChiralTag, value: PropertyValue) -> Atom {
        let mut atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C).with_chiral_tag(tag),
        );
        atom.set_prop("_chiralPermutation", value).unwrap();
        atom
    }
    #[test]
    fn source_atom_inversion_reads_native_unsigned_property_coercions_and_writes_uint() {
        for (tag, value, expected) in [
            (ChiralTag::Tetrahedral, PropertyValue::Int(1), 2),
            (
                ChiralTag::Tetrahedral,
                PropertyValue::String("+2  ".into()),
                1,
            ),
            (ChiralTag::Tetrahedral, PropertyValue::UInt(u32::MAX), 0),
            (
                ChiralTag::TrigonalBipyramidal,
                PropertyValue::String("3".into()),
                4,
            ),
            (ChiralTag::Octahedral, PropertyValue::Int(3), 16),
        ] {
            let mut atom = atom(tag, value);
            assert_eq!(invert_atom_chirality(&mut atom).unwrap(), expected != 0);
            assert_eq!(
                atom.prop("_chiralPermutation"),
                Some(&PropertyValue::UInt(expected))
            );
            assert_eq!(atom.chiral_permutation(), None); // raw-only state stays raw-only
        }
    }
    #[test]
    fn source_atom_inversion_wrong_property_kind_or_lexical_value_errors_without_mutation() {
        for value in [
            PropertyValue::Bool(true),
            PropertyValue::String("1x".into()),
            PropertyValue::Int(-1),
            PropertyValue::String("".into()),
        ] {
            let mut atom = atom(ChiralTag::Octahedral, value);
            let before = atom.clone();
            assert!(matches!(
                invert_atom_chirality(&mut atom),
                Err(StereoOrderError::ChiralPermutationRead(_))
            ));
            assert_eq!(atom, before);
        }
    }
    #[test]
    fn source_atom_inversion_cw_and_unhandled_tags_do_not_read_bad_permutation_property() {
        for (tag, expected_tag, changed) in [
            (ChiralTag::TetrahedralCw, ChiralTag::TetrahedralCcw, true),
            (ChiralTag::TetrahedralCcw, ChiralTag::TetrahedralCw, true),
            (ChiralTag::Unspecified, ChiralTag::Unspecified, false),
            (ChiralTag::Other, ChiralTag::Other, false),
        ] {
            let mut atom = atom(tag, PropertyValue::Bool(false));
            atom.set_chiral_permutation(Some(1));
            let mut expected = atom.clone();
            expected.set_chiral_tag(expected_tag);
            assert_eq!(invert_atom_chirality(&mut atom).unwrap(), changed);
            assert_eq!(atom, expected);
        }
    }
    #[test]
    fn source_atom_inversion_noncomputed_write_does_not_read_malformed_computed_list() {
        let mut atom = atom(ChiralTag::TrigonalBipyramidal, PropertyValue::UInt(20));
        atom.set_prop("__computedProps", PropertyValue::Int(7))
            .unwrap();
        assert!(invert_atom_chirality(&mut atom).unwrap());
        assert_eq!(
            atom.prop("_chiralPermutation"),
            Some(&PropertyValue::UInt(15))
        );
        assert_eq!(atom.prop("__computedProps"), Some(&PropertyValue::Int(7)));
    }
    #[test]
    fn source_atom_inversion_conflicting_existing_transport_projections_are_structural_error() {
        let mut atom = atom(ChiralTag::Tetrahedral, PropertyValue::UInt(1));
        atom.set_chiral_permutation(Some(2));
        let before = atom.clone();
        assert_eq!(
            invert_atom_chirality(&mut atom),
            Err(StereoOrderError::ConflictingChiralPermutation {
                atom: AtomId::new(0),
                typed: 2,
                raw: 1
            })
        );
        assert_eq!(atom, before);
    }
    #[test]
    fn source_atom_inversion_raw_computed_permutation_clear_restores_absent_property_branch() {
        let mut atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TrigonalBipyramidal),
        );
        atom.set_computed_prop("_chiralPermutation", PropertyValue::UInt(20))
            .unwrap();
        assert!(invert_atom_chirality(&mut atom).unwrap());
        assert_eq!(
            atom.prop("_chiralPermutation"),
            Some(&PropertyValue::UInt(15))
        );
        assert_eq!(atom.chiral_permutation(), None);
        atom.clear_computed_props().unwrap();
        let before = atom.clone();
        assert!(!invert_atom_chirality(&mut atom).unwrap());
        assert_eq!(atom, before);
    }
}
