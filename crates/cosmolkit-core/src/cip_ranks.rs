use cosmolkit_model::{
    AtomId, BondId, QueryStateError, QueryStateRef, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::BondOrder;

use crate::{ValenceAssignment, most_common_isotope};

const MAX_CIP_BONDS: usize = 16;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum CipRankError {
    #[error("invalid topology: {0}")]
    InvalidTopology(TopologyValidationError),
    #[error(transparent)]
    InvalidQueryState(#[from] QueryStateError),
    #[error("valence row {field} has length {actual}, expected {atom_count}")]
    ValenceRowCount {
        field: &'static str,
        actual: usize,
        atom_count: usize,
    },
    #[error("implicit hydrogen value for atom {atom} must be nonnegative, got {value}")]
    NegativeImplicitHydrogen { atom: AtomId, value: i32 },
    #[error("atom map number {map_number} on atom {atom} exceeds the source signed-int range")]
    AtomMapOutOfRange { atom: AtomId, map_number: u32 },
    #[error("CIP invariant row has length {actual}, expected {atom_count}")]
    InvariantCount { actual: usize, atom_count: usize },
    #[error("CIP invariant {value} for atom {atom} exceeds the source signed-int range")]
    InvariantOutOfRange { atom: AtomId, value: i64 },
    #[error(
        "atom {atom} has {degree} neighbors, but the source CIP feature stride supports at most {maximum_supported}"
    )]
    TooManyNeighbors {
        atom: AtomId,
        degree: usize,
        maximum_supported: usize,
    },
    #[error("bond {bond} has unsupported CIP-ranking order {order:?}")]
    UnsupportedBondOrder { bond: BondId, order: BondOrder },
}

/// Assigns RDKit legacy atom CIP symmetry ranks over detached topology values.
pub fn assign_atom_cip_ranks(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<Vec<u32>, CipRankError> {
    assign_atom_cip_ranks_with_query_state(topology, valence, None)
}

#[doc(hidden)]
pub fn assign_atom_cip_ranks_with_query_state(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<Vec<u32>, CipRankError> {
    // RDKit❗✔️: void assignAtomCIPRanks(const ROMol &mol, UINT_VECT &ranks) {
    // RDKit❗✔️:   PRECONDITION((!ranks.size() || ranks.size() >= mol.getNumAtoms()),
    // RDKit❗✔️:                "bad ranks size");
    // RDKit❗✔️:   if (!ranks.size()) {
    // RDKit❗✔️:     ranks.resize(mol.getNumAtoms());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   unsigned int numAtoms = mol.getNumAtoms();
    // RDKit❗✔️: #ifndef USE_NEW_STEREOCHEMISTRY
    // RDKit❗✔️:   // get the initial invariants:
    // RDKit❗✔️:   DOUBLE_VECT invars(numAtoms, 0);
    // RDKit❗✔️:   buildCIPInvariants(mol, invars);
    // RDKit❗✔️:   iterateCIPRanks(mol, invars, ranks, false);
    // RDKit❌❌: #else
    // RDKit❌❌:   Canon::chiralRankMolAtoms(mol, ranks);
    // RDKit❗✔️: #endif
    // RDKit❌❌:   // copy the ranks onto the atoms:
    // RDKit❌❌:   for (unsigned int i = 0; i < numAtoms; ++i) {
    // RDKit❌❌:     mol[i]->setProp(common_properties::_CIPRank, ranks[i], 1);
    // RDKit❌❌:   }
    // RDKit❗✔️: }
    // The detached API deliberately returns an exactly sized row instead of
    // mutating caller storage or atom properties. The selected legacy branch
    // is otherwise implemented by the source-shaped helpers below.
    validate_inputs(topology, valence)?;
    if let Some(state) = query_state {
        QueryStateRef::try_for_topology(state.atoms(), state.bonds(), topology)?;
    }
    let invariants = build_cip_invariants(topology)?;
    iterate_cip_ranks(topology, valence, &invariants, false, query_state)
}

/// Refines caller-built integer invariants using RDKit's legacy CIP kernel.
///
/// High-level stereo code owns construction of chirality-enriched invariants;
/// this primitive only reproduces `iterateCIPRanks(..., seedWithInvars=true)`.
pub fn refine_atom_cip_ranks_from_invariants(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    invariants: &[i64],
) -> Result<Vec<u32>, CipRankError> {
    refine_atom_cip_ranks_from_invariants_with_query_state(topology, valence, invariants, None)
}

#[doc(hidden)]
pub fn refine_atom_cip_ranks_from_invariants_with_query_state(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    invariants: &[i64],
    query_state: Option<QueryStateRef<'_>>,
) -> Result<Vec<u32>, CipRankError> {
    validate_inputs(topology, valence)?;
    if let Some(state) = query_state {
        QueryStateRef::try_for_topology(state.atoms(), state.bonds(), topology)?;
    }
    if invariants.len() != topology.atoms.len() {
        return Err(CipRankError::InvariantCount {
            actual: invariants.len(),
            atom_count: topology.atoms.len(),
        });
    }
    iterate_cip_ranks(topology, valence, invariants, true, query_state)
}

fn validate_inputs(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<(), CipRankError> {
    topology.validate().map_err(CipRankError::InvalidTopology)?;
    let atom_count = topology.atoms.len();
    if valence.explicit_valence.len() != atom_count {
        return Err(CipRankError::ValenceRowCount {
            field: "explicit_valence",
            actual: valence.explicit_valence.len(),
            atom_count,
        });
    }
    if valence.implicit_hydrogens.len() != atom_count {
        return Err(CipRankError::ValenceRowCount {
            field: "implicit_hydrogens",
            actual: valence.implicit_hydrogens.len(),
            atom_count,
        });
    }
    for (index, &value) in valence.implicit_hydrogens.iter().enumerate() {
        if value < 0 {
            return Err(CipRankError::NegativeImplicitHydrogen {
                atom: AtomId::new(index),
                value,
            });
        }
    }
    for atom in &topology.atoms {
        if let Some(map_number) = atom.atom_map()
            && map_number > i32::MAX as u32
        {
            return Err(CipRankError::AtomMapOutOfRange {
                atom: atom.id(),
                map_number,
            });
        }
        let degree = topology.adjacency.neighbors_of(atom.id().index()).len();
        if degree >= MAX_CIP_BONDS {
            return Err(CipRankError::TooManyNeighbors {
                atom: atom.id(),
                degree,
                maximum_supported: MAX_CIP_BONDS - 1,
            });
        }
    }
    Ok(())
}

fn build_cip_invariants(topology: &TopologyBlock) -> Result<Vec<i64>, CipRankError> {
    // RDKit❗✔️: void buildCIPInvariants(const ROMol &mol, DOUBLE_VECT &res) {
    // RDKit❗✔️:   PRECONDITION(res.size() >= mol.getNumAtoms(), "res vect too small");
    // RDKit❗✔️:   int atsSoFar = 0;
    // RDKit❗✔️:   //
    // RDKit❗✔️:   // NOTE:
    // RDKit❗✔️:   // If you make modifications to this, keep in mind that it is
    // RDKit❗✔️:   // essential that the initial comparison of ranks behave properly.
    // RDKit❗✔️:   // So, though it seems like it would makes sense to include
    // RDKit❗✔️:   // information about the number of Hs (or charge, etc) in the CIP
    // RDKit❗✔️:   // invariants, this will result in bad rankings.  For example, in
    // RDKit❗✔️:   // this molecule: OC[C@H](C)O, including the number of Hs would
    // RDKit❗✔️:   // cause the methyl group (atom 3) to be ranked higher than the CH2
    // RDKit❗✔️:   // connected to O (atom 1).  This is totally wrong.
    // RDKit❗✔️:   //
    // RDKit❗✔️:   // We also don't include any pre-existing stereochemistry information.
    // RDKit❗✔️:   // Though R and S assignments do factor in to the priorities of atoms,
    // RDKit❗✔️:   // we're starting here from scratch and we'll let the R and S stuff
    // RDKit❗✔️:   // be taken into account during the iterations.
    // RDKit❗✔️:   //
    // RDKit❗✔️:   for (const auto atom : mol.atoms()) {
    // RDKit❗✔️:     const unsigned short nMassBits = 10;
    // RDKit❗✔️:     const unsigned short maxMass = 1 << nMassBits;
    // RDKit❗✔️:     unsigned long invariant = 0;
    // RDKit❗✔️:     int num = atom->getAtomicNum() % 128;
    // RDKit❗✔️:     // get an int with the deviation in the mass from the default:
    // RDKit❗✔️:     int mass = 0;
    // RDKit❗✔️:     if (atom->getIsotope()) {
    // RDKit❗✔️:       mass =
    // RDKit❗✔️:           atom->getIsotope() -
    // RDKit❗✔️:           PeriodicTable::getTable()->getMostCommonIsotope(atom->getAtomicNum());
    // RDKit❗✔️:       if (mass >= 0) {
    // RDKit❗✔️:         mass += 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     mass += maxMass / 2;
    // RDKit❗✔️:     if (mass < 0) {
    // RDKit❗✔️:       mass = 0;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       mass = mass % maxMass;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     invariant = num;  // 7 bits here
    // RDKit❗✔️:     invariant = (invariant << nMassBits) | mass;
    // RDKit❗✔️:
    // RDKit❗✔️:     int mapnum = -1;
    // RDKit❗✔️:     atom->getPropIfPresent(common_properties::molAtomMapNumber, mapnum);
    // RDKit❗✔️:     mapnum = (mapnum + 1) % 1024;  // increment to allow map numbers of zero
    // RDKit❗✔️:                                      // (though that would be stupid)
    // RDKit❗✔️:     invariant = (invariant << 10) | mapnum;
    // RDKit❗✔️:
    // RDKit❗✔️:     res[atsSoFar++] = invariant;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let mut result = Vec::with_capacity(topology.atoms.len());
    for atom in &topology.atoms {
        let number = i64::from(atom.atomic_number() % 128);
        let mut mass = 0_i64;
        if let Some(isotope) = atom.isotope() {
            mass = i64::from(isotope) - i64::from(most_common_isotope(atom.element()));
            if mass >= 0 {
                mass += 1;
            }
        }
        mass += 1_i64 << 9;
        if mass < 0 {
            mass = 0;
        } else {
            mass %= 1_i64 << 10;
        }
        let map_number = atom
            .atom_map()
            .map_or(0, |value| (i64::from(value) + 1) % 1024);
        result.push((((number << 10) | mass) << 10) | map_number);
    }
    Ok(result)
}

#[derive(Debug)]
struct BondFeatures {
    counts_and_neighbor_indices: Vec<(u8, usize)>,
    num_neighbors: Vec<usize>,
}

fn compute_bond_features(topology: &TopologyBlock) -> Result<BondFeatures, CipRankError> {
    // RDKit❗✔️: PrecomputedBondFeatures computeBondFeatures(const ROMol &mol) {
    // RDKit❗✔️:   PrecomputedBondFeatures features;
    // RDKit❗✔️:   const unsigned int numAtoms = mol.getNumAtoms();
    // RDKit❗✔️:   features.countsAndNeighborIndices.resize(numAtoms * kMaxBonds);
    // RDKit❗✔️:   features.numNeighbors.resize(numAtoms, 0);
    // RDKit❗✔️:
    // RDKit❗✔️:   for (size_t atomIdx = 0; atomIdx < numAtoms; atomIdx++) {
    // RDKit❗✔️:     int indexOffset = atomIdx * kMaxBonds;
    // RDKit❗✔️:     for (const auto bond : mol.atomBonds(mol[atomIdx])) {
    // RDKit❗✔️:       const unsigned int nbrIdx = bond->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:       features.numNeighbors[nbrIdx]++;
    // RDKit❗✔️:       auto &[count, neighborIndex] =
    // RDKit❗✔️:           features.countsAndNeighborIndices.at(indexOffset);
    // RDKit❗✔️:       neighborIndex = nbrIdx;
    // RDKit❗✔️:
    // RDKit❗✔️:       // put the neighbor in 2N times where N is the bond order as a double.
    // RDKit❗✔️:       // this is to treat aromatic linkages on fair footing. i.e. at least in
    // RDKit❗✔️:       // the first iteration --c(:c):c and --C(=C)-C should look the same.
    // RDKit❗✔️:       // this was part of issue 3009911
    // RDKit❗✔️:
    // RDKit❗✔️:       // a special case for chiral phosphorus compounds
    // RDKit❗✔️:       // (this was leading to incorrect assignment of R/S labels ):
    // RDKit❗✔️:       bool isChiralPhosphorusSpecialCase = false;
    // RDKit❗✔️:       if (bond->getBondType() == Bond::DOUBLE) {
    // RDKit❗✔️:         const Atom *nbr = mol[nbrIdx];
    // RDKit❗✔️:         if (nbr->getAtomicNum() == 15) {
    // RDKit❗✔️:           unsigned int nbrDeg = nbr->getDegree();
    // RDKit❗✔️:           isChiralPhosphorusSpecialCase = nbrDeg == 3 || nbrDeg == 4;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       };
    // RDKit❗✔️:
    // RDKit❗✔️:       // general justification of this is:
    // RDKit❗✔️:       // Paragraph 2.2. in the 1966 article is "Valence-Bond Conventions:
    // RDKit❗✔️:       // Multiple-Bond Unsaturation and Aromaticity". It contains several
    // RDKit❗✔️:       // conventions of which convention (b) is the one applying here:
    // RDKit❗✔️:       // "(b) Contributions by d orbitals to bonds of quadriligant atoms are
    // RDKit❗✔️:       // neglected."
    // RDKit❗✔️:       // FIX: this applies to more than just P
    // RDKit❗✔️:       if (isChiralPhosphorusSpecialCase) {
    // RDKit❗✔️:         count += 1;
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         count += getTwiceBondType(*bond);
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       ++indexOffset;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return features;
    // RDKit❗✔️: }
    let mut result = BondFeatures {
        counts_and_neighbor_indices: vec![(0, 0); topology.atoms.len() * MAX_CIP_BONDS],
        num_neighbors: vec![0; topology.atoms.len()],
    };
    for atom_index in 0..topology.atoms.len() {
        let mut feature_index = atom_index * MAX_CIP_BONDS;
        for neighbor in topology.adjacency.neighbors_of(atom_index) {
            result.num_neighbors[neighbor.atom_index] += 1;
            let bond = &topology.bonds[neighbor.bond.index()];
            let phosphorus_special_case = bond.order() == BondOrder::Double
                && topology.atoms[neighbor.atom_index].atomic_number() == 15
                && matches!(
                    topology.adjacency.neighbors_of(neighbor.atom_index).len(),
                    3 | 4
                );
            result.counts_and_neighbor_indices[feature_index] = (
                if phosphorus_special_case {
                    1
                } else {
                    twice_bond_type(bond.id(), bond.order())?
                },
                neighbor.atom_index,
            );
            feature_index += 1;
        }
    }
    Ok(result)
}

fn twice_bond_type(bond: BondId, order: BondOrder) -> Result<u8, CipRankError> {
    // RDKit✔️✔️: uint8_t getTwiceBondType(const Bond &b) {
    // RDKit✔️✔️:   switch (b.getBondType()) {
    // RDKit✔️✔️:     case Bond::UNSPECIFIED:
    // RDKit✔️✔️:     case Bond::IONIC:
    // RDKit✔️✔️:     case Bond::ZERO:
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::SINGLE:
    // RDKit✔️✔️:       return 2;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::DOUBLE:
    // RDKit✔️✔️:       return 4;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::TRIPLE:
    // RDKit✔️✔️:       return 6;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::QUADRUPLE:
    // RDKit✔️✔️:       return 8;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::QUINTUPLE:
    // RDKit✔️✔️:       return 10;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::HEXTUPLE:
    // RDKit✔️✔️:       return 12;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::ONEANDAHALF:
    // RDKit✔️✔️:       return 3;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::TWOANDAHALF:
    // RDKit✔️✔️:       return 5;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::THREEANDAHALF:
    // RDKit✔️✔️:       return 7;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::FOURANDAHALF:
    // RDKit✔️✔️:       return 9;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::FIVEANDAHALF:
    // RDKit✔️✔️:       return 11;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::AROMATIC:
    // RDKit✔️✔️:       return 3;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::DATIVEONE:
    // RDKit✔️✔️:       return 2;
    // RDKit✔️✔️:       break;  // FIX: this should probably be different
    // RDKit✔️✔️:     case Bond::DATIVE:
    // RDKit✔️✔️:       return 2;
    // RDKit✔️✔️:       break;  // FIX: again probably wrong
    // RDKit✔️✔️:     case Bond::HYDROGEN:
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:       break;
    // RDKit❌❌:     default:
    // RDKit❌❌:       UNDER_CONSTRUCTION("Bad bond type");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    match order {
        BondOrder::Unspecified | BondOrder::Ionic | BondOrder::Zero => Ok(0),
        BondOrder::Single | BondOrder::DativeOne | BondOrder::Dative => Ok(2),
        BondOrder::Double => Ok(4),
        BondOrder::Triple => Ok(6),
        BondOrder::Quadruple => Ok(8),
        BondOrder::Quintuple => Ok(10),
        BondOrder::Hextuple => Ok(12),
        BondOrder::OneAndHalf | BondOrder::Aromatic => Ok(3),
        BondOrder::TwoAndHalf => Ok(5),
        BondOrder::ThreeAndHalf => Ok(7),
        BondOrder::FourAndHalf => Ok(9),
        BondOrder::FiveAndHalf => Ok(11),
        BondOrder::Hydrogen => Ok(0),
        BondOrder::DativeLeft
        | BondOrder::DativeRight
        | BondOrder::ThreeCenter
        | BondOrder::Other => Err(CipRankError::UnsupportedBondOrder { bond, order }),
    }
}

fn iterate_cip_ranks(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    invariants: &[i64],
    seed_with_invariants: bool,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<Vec<u32>, CipRankError> {
    // RDKit❗✔️: void iterateCIPRanks(const ROMol &mol, const DOUBLE_VECT &invars,
    // RDKit❗✔️:                      UINT_VECT &ranks, bool seedWithInvars) {
    // RDKit❗✔️:   PRECONDITION(invars.size() == mol.getNumAtoms(), "bad invars size");
    // RDKit❗✔️:   PRECONDITION(ranks.size() >= mol.getNumAtoms(), "bad ranks size");
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int numAtoms = mol.getNumAtoms();
    // RDKit❗✔️:   CIP_ENTRY_VECT cipEntries(numAtoms);
    // RDKit❗✔️:   for (auto &vec : cipEntries) {
    // RDKit❗✔️:     vec.reserve(16);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<SortableCIPReference> sortableEntries;
    // RDKit❗✔️:   sortableEntries.reserve(numAtoms);
    // RDKit❗✔️:   for (size_t i = 0; i < cipEntries.size(); i++) {
    // RDKit❗✔️:     sortableEntries.emplace_back(&cipEntries[i], i);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (unsigned int i = 0; i < numAtoms; i++) {
    // RDKit❗✔️:     cipEntries[i].push_back(static_cast<int>(invars[i]));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   unsigned int numRanks;
    // RDKit❗✔️:   std::sort(sortableEntries.begin(), sortableEntries.end());
    // RDKit❗✔️:   std::vector<std::pair<int, int>> needsSorting;
    // RDKit❗✔️:   findSegmentsToResort(sortableEntries, needsSorting, numRanks);
    // RDKit❗✔️:   recomputeRanks(sortableEntries, ranks);
    // RDKit❗✔️:
    // RDKit❗✔️:   // Start each atom's rank vector with its atomic number:
    // RDKit❗✔️:   //  Note: in general one should avoid the temptation to
    // RDKit❗✔️:   //  use invariants here, those lead to incorrect answers
    // RDKit❗✔️:   for (unsigned int i = 0; i < numAtoms; i++) {
    // RDKit❗✔️:     if (seedWithInvars) {
    // RDKit❗✔️:       cipEntries[i][0] = static_cast<int>(invars[i]);
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       cipEntries[i][0] = mol[i]->getAtomicNum();
    // RDKit❗✔️:       cipEntries[i].push_back(static_cast<int>(ranks[i]));
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // Based on above seeding, the rank will be set at index 1 or 2.
    // RDKit❗✔️:   const int cipRankIndex = seedWithInvars ? 1 : 2;
    // RDKit❗✔️:
    // RDKit❗✔️:   // Loop until either:
    // RDKit❗✔️:   //   1) all classes are uniquified
    // RDKit❗✔️:   //   2) the number of ranks doesn't change from one iteration to
    // RDKit❗✔️:   //      the next
    // RDKit❗✔️:   //   3) we've gone through maxIts times
    // RDKit❗✔️:   //      maxIts is calculated by dividing the number of atoms
    // RDKit❗✔️:   //      by 2. That's a pessimal version of the
    // RDKit❗✔️:   //      maximum number of steps required for two atoms to
    // RDKit❗✔️:   //      "feel" each other (each influences one additional
    // RDKit❗✔️:   //      neighbor shell per iteration).
    // RDKit❗✔️:   unsigned int maxIts = numAtoms / 2 + 1;
    // RDKit❗✔️:   unsigned int numIts = 0;
    // RDKit❗✔️:   int lastNumRanks = -1;
    // RDKit❗✔️:
    // RDKit❗✔️:   PrecomputedBondFeatures bondFeatures = computeBondFeatures(mol);
    // RDKit❗✔️:
    // RDKit❗✔️:   while (!needsSorting.empty() && numIts < maxIts &&
    // RDKit❗✔️:          (lastNumRanks < 0 ||
    // RDKit❗✔️:           static_cast<unsigned int>(lastNumRanks) < numRanks)) {
    // RDKit❗✔️:     // ----------------------------------------------------
    // RDKit❗✔️:     //
    // RDKit❗✔️:     // for each atom, get a sorted list of its neighbors' ranks:
    // RDKit❗✔️:     //
    // RDKit❗✔️:     for (unsigned int index = 0; index < numAtoms; ++index) {
    // RDKit❗✔️:       const unsigned int indexOffset = kMaxBonds * index;
    // RDKit❗✔️:       const int numNeighbors = bondFeatures.numNeighbors[index];
    // RDKit❗✔️:
    // RDKit❗✔️:       auto *sortBegin = &bondFeatures.countsAndNeighborIndices[indexOffset];
    // RDKit❗✔️:       auto *sortEnd = sortBegin + numNeighbors + 1;
    // RDKit❗✔️:
    // RDKit❗✔️:       // For each of our neighbors' ranks weighted by bond type, copy it N times
    // RDKit❗✔️:       // to our cipEntry in reverse rank order, where N is the weight.
    // RDKit❗✔️:       if (numNeighbors > 1) {  // compare vs 1 for performance.
    // RDKit❗✔️:         std::sort(sortBegin, sortEnd,
    // RDKit❗✔️:                   [&ranks](const std::pair<std::uint8_t, int> &countAndIdx1,
    // RDKit❗✔️:                            const std::pair<std::uint8_t, int> &countAndIdx2) {
    // RDKit❗✔️:                     return ranks[countAndIdx1.second] >
    // RDKit❗✔️:                            ranks[countAndIdx2.second];
    // RDKit❗✔️:                   });
    // RDKit❗✔️:       }
    // RDKit❗✔️:       auto &cipEntry = cipEntries[index];
    // RDKit❗✔️:       for (auto *iter = sortBegin; iter != sortEnd; ++iter) {
    // RDKit❗✔️:         const auto &[count, idx] = *iter;
    // RDKit❗✔️:         cipEntry.insert(cipEntry.end(), count, ranks[idx] + 1);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       // add a zero for each coordinated H as long as we're not a query atom
    // RDKit❗✔️:       if (!mol[index]->hasQuery()) {
    // RDKit❗✔️:         cipEntry.insert(cipEntry.end(), mol[index]->getTotalNumHs(), 0);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // ----------------------------------------------------
    // RDKit❗✔️:     //
    // RDKit❗✔️:     // sort the new ranks and update the list of active indices:
    // RDKit❗✔️:     //
    // RDKit❗✔️:     lastNumRanks = numRanks;
    // RDKit❗✔️:
    // RDKit❗✔️:     // Loop through previously tied atom sections and re-sort.
    // RDKit❗✔️:     for (const auto &[firstIdx, lastIdx] : needsSorting) {
    // RDKit❗✔️:       std::sort(sortableEntries.begin() + firstIdx,
    // RDKit❗✔️:                 sortableEntries.begin() + lastIdx + 1);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     findSegmentsToResort(sortableEntries, needsSorting, numRanks);
    // RDKit❗✔️:     // Map out of order rankings back to the absolute rankings vector.
    // RDKit❗✔️:     recomputeRanks(sortableEntries, ranks);
    // RDKit❗✔️:
    // RDKit❗✔️:     // now truncate each vector and stick the rank at the end
    // RDKit❗✔️:     if (static_cast<unsigned int>(lastNumRanks) != numRanks) {
    // RDKit❗✔️:       for (unsigned int i = 0; i < numAtoms; ++i) {
    // RDKit❗✔️:         cipEntries[i].resize(cipRankIndex + 1);
    // RDKit❗✔️:         cipEntries[i][cipRankIndex] = ranks[i];
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     ++numIts;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let atom_count = topology.atoms.len();
    if invariants.len() != atom_count {
        return Err(CipRankError::InvariantCount {
            actual: invariants.len(),
            atom_count,
        });
    }
    let mut entries = Vec::with_capacity(atom_count);
    for (index, &value) in invariants.iter().enumerate() {
        let value = i32::try_from(value).map_err(|_| CipRankError::InvariantOutOfRange {
            atom: AtomId::new(index),
            value,
        })?;
        entries.push(vec![value]);
    }
    if atom_count == 0 {
        return Ok(Vec::new());
    }
    let mut ranks = vec![0; atom_count];
    let mut order = (0..atom_count).collect::<Vec<_>>();
    order.sort_unstable_by(|left, right| entries[*left].cmp(&entries[*right]));
    let (mut tied_ranges, mut rank_count) = assign_cip_ranks(&entries, &order, &mut ranks);
    for atom_index in 0..atom_count {
        if !seed_with_invariants {
            entries[atom_index][0] = i32::from(topology.atoms[atom_index].atomic_number());
            entries[atom_index].push(i32::try_from(ranks[atom_index]).unwrap_or(i32::MAX));
        }
    }
    let rank_index = usize::from(!seed_with_invariants) + 1;
    let mut features = compute_bond_features(topology)?;
    let max_iterations = atom_count / 2 + 1;
    let mut iteration = 0;
    let mut last_rank_count = None;
    while !tied_ranges.is_empty()
        && iteration < max_iterations
        && last_rank_count.is_none_or(|last| last < rank_count)
    {
        for atom_index in 0..atom_count {
            let start = atom_index * MAX_CIP_BONDS;
            let end = start + features.num_neighbors[atom_index] + 1;
            let neighbors = &mut features.counts_and_neighbor_indices[start..end];
            if features.num_neighbors[atom_index] > 1 {
                neighbors.sort_unstable_by(|left, right| ranks[right.1].cmp(&ranks[left.1]));
            }
            for &(count, neighbor) in neighbors.iter() {
                entries[atom_index].extend(std::iter::repeat_n(
                    i32::try_from(ranks[neighbor]).unwrap_or(i32::MAX) + 1,
                    usize::from(count),
                ));
            }
            let total_hydrogens = usize::from(topology.atoms[atom_index].explicit_hydrogens())
                + usize::try_from(valence.implicit_hydrogens[atom_index])
                    .expect("negative implicit hydrogen rejected before iteration");
            // RDKit✔️✔️:       // add a zero for each coordinated H as long as we're not a query atom
            // RDKit✔️✔️:       if (!mol[index]->hasQuery()) {
            // RDKit✔️✔️:         cipEntry.insert(cipEntry.end(), mol[index]->getTotalNumHs(), 0);
            // RDKit✔️✔️:       }
            // Behavior review: Explicit typed rows reproduce QueryAtom identity;
            // CarrierDerived and absent state reproduce ordinary Atom identity.
            // Complexity review: one O(1) origin lookup replaces the source
            // virtual `hasQuery()` test and does not change entry allocation.
            if !query_state.is_some_and(|state| state.atom_has_query(AtomId::new(atom_index))) {
                entries[atom_index].extend(std::iter::repeat_n(0, total_hydrogens));
            }
        }
        last_rank_count = Some(rank_count);
        for &(first, last) in &tied_ranges {
            order[first..=last]
                .sort_unstable_by(|left, right| entries[*left].cmp(&entries[*right]));
        }
        (tied_ranges, rank_count) = assign_cip_ranks(&entries, &order, &mut ranks);
        if last_rank_count != Some(rank_count) {
            for atom_index in 0..atom_count {
                entries[atom_index].resize(rank_index + 1, 0);
                entries[atom_index][rank_index] =
                    i32::try_from(ranks[atom_index]).unwrap_or(i32::MAX);
            }
        }
        iteration += 1;
    }
    Ok(ranks)
}

fn assign_cip_ranks(
    entries: &[Vec<i32>],
    order: &[usize],
    ranks: &mut [u32],
) -> (Vec<(usize, usize)>, usize) {
    // RDKit❗✔️: void findSegmentsToResort(std::vector<SortableCIPReference> &sortedEntries,
    // RDKit❗✔️:                           std::vector<std::pair<int, int>> &res,
    // RDKit❗✔️:                           unsigned int &numIndependentEntries) {
    // RDKit❗✔️:   res.clear();
    // RDKit❗✔️:   numIndependentEntries = rdcast<unsigned int>(sortedEntries.size());
    // RDKit❗✔️:   SortableCIPReference *current = &sortedEntries.front();
    // RDKit❗✔️:   int runningRank = 0;
    // RDKit❗✔️:   current->currRank = runningRank;
    // RDKit❗✔️:   bool inEqualSection = false;
    // RDKit❗✔️:
    // RDKit❗✔️:   for (size_t i = 1; i < sortedEntries.size(); i++) {
    // RDKit❗✔️:     SortableCIPReference &entry = sortedEntries[i];
    // RDKit❗✔️:     if (*current == entry) {
    // RDKit❗✔️:       entry.currRank = runningRank;
    // RDKit❗✔️:       numIndependentEntries--;
    // RDKit❗✔️:       // Case where we need to open a section
    // RDKit❗✔️:       if (!inEqualSection) {
    // RDKit❗✔️:         inEqualSection = true;
    // RDKit❗✔️:         auto &[firstIndex, _] = res.emplace_back();
    // RDKit❗✔️:         // Go back to the first in this section, we only catch at first + 1
    // RDKit❗✔️:         firstIndex = i - 1;
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         // Case where we are already in a section, nullop
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       // Case where we're closing an open section.
    // RDKit❗✔️:       runningRank++;
    // RDKit❗✔️:       entry.currRank = runningRank;
    // RDKit❗✔️:       current = &entry;
    // RDKit❗✔️:
    // RDKit❗✔️:       if (inEqualSection) {
    // RDKit❗✔️:         auto &[_, finalIndex] = res.back();
    // RDKit❗✔️:         finalIndex = i;
    // RDKit❗✔️:         inEqualSection = false;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // Handle currently open.
    // RDKit❗✔️:   if (inEqualSection) {
    // RDKit❗✔️:     auto &[_, finalIndex] = res.back();
    // RDKit❗✔️:     finalIndex = sortedEntries.size() - 1;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: void recomputeRanks(const std::vector<SortableCIPReference> &sortedEntries,
    // RDKit❗✔️:                     std::vector<unsigned int> &ranks) {
    // RDKit❗✔️:   for (size_t rank = 0; rank < ranks.size(); ++rank) {
    // RDKit❗✔️:     const auto &cipEntry = sortedEntries[rank];
    // RDKit❗✔️:     ranks[cipEntry.atomIdx] = cipEntry.currRank;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let mut tied = Vec::new();
    let mut independent = order.len();
    let mut running_rank = 0_u32;
    let mut current = 0;
    let mut open_tie = None;
    ranks[order[0]] = running_rank;
    for position in 1..order.len() {
        if entries[order[current]] == entries[order[position]] {
            ranks[order[position]] = running_rank;
            independent -= 1;
            open_tie.get_or_insert(position - 1);
        } else {
            running_rank += 1;
            ranks[order[position]] = running_rank;
            current = position;
            if let Some(first) = open_tie.take() {
                tied.push((first, position));
            }
        }
    }
    if let Some(first) = open_tie {
        tied.push((first, order.len() - 1));
    }
    (tied, independent)
}

#[cfg(test)]
mod q05_tests {
    use super::*;
    use cosmolkit_model::{
        Atom, AtomQueryPredicate, AtomSpec, Bond, BondQueryPredicate, BondSpec, QueryAtom,
        QueryBond, QueryNode, QueryStateRef,
    };
    use cosmolkit_types::Element;

    #[test]
    fn q05_core_sanitize_query_consumers_legacy_cip_omits_query_atom_hydrogen_zeros() {
        let atoms = [Element::C, Element::F, Element::C, Element::C]
            .into_iter()
            .enumerate()
            .map(|(index, element)| Atom::from_spec(AtomId::new(index), AtomSpec::new(element)))
            .collect::<Vec<_>>();
        let bonds = [(0, 1), (0, 2), (0, 3)]
            .into_iter()
            .enumerate()
            .map(|(index, (begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect::<Vec<_>>();
        let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        let valence = ValenceAssignment {
            explicit_valence: vec![3, 1, 1, 1],
            implicit_hydrogens: vec![1, 0, 3, 3],
        };
        let carrier_atoms = topology
            .atoms
            .iter()
            .map(|carrier| {
                QueryAtom::from_carrier_parts(
                    carrier.clone(),
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(carrier.atomic_number())),
                )
            })
            .collect::<Vec<_>>();
        let query_bonds = topology
            .bonds
            .iter()
            .map(|carrier| {
                QueryBond::from_carrier_parts(
                    carrier.clone(),
                    QueryNode::predicate(BondQueryPredicate::Order(carrier.order())),
                )
            })
            .collect::<Vec<_>>();
        let carrier_state =
            QueryStateRef::try_for_topology(&carrier_atoms, &query_bonds, &topology).unwrap();
        assert!(!carrier_state.atom_has_query(AtomId::new(2)));
        assert_eq!(
            assign_atom_cip_ranks(&topology, &valence).unwrap(),
            vec![1, 2, 0, 0]
        );

        let mut explicit_atoms = carrier_atoms;
        explicit_atoms[2] = QueryAtom::from_parts(
            topology.atoms[2].clone(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        );
        let explicit_state =
            QueryStateRef::try_for_topology(&explicit_atoms, &query_bonds, &topology).unwrap();
        assert!(explicit_state.atom_has_query(AtomId::new(2)));
        assert_eq!(
            assign_atom_cip_ranks_with_query_state(&topology, &valence, Some(explicit_state))
                .unwrap(),
            vec![2, 3, 0, 1]
        );
    }
}
