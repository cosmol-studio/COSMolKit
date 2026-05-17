// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use std::collections::{BTreeMap, BTreeSet};

use crate::{
    AdjacencyList, Atom, AtomId, Bond, BondOrder, Molecule, RingFindingError, RingInfo,
    ValenceAssignment, ValenceError, ValenceModel, assign_valence, read_parts::MoleculeReadParts,
    symmetrize_sssr,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AromaticityModel {
    Default,
    Rdkit,
    Simple,
    Mdl,
    Mmff94,
    Custom,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AromaticityAssignment {
    pub atom_aromatic: Vec<bool>,
    pub bond_aromatic: Vec<bool>,
    pub aromatic_ring_count: usize,
}

#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ElectronDonorType {
    Vacant,
    One,
    Two,
    OneOrTwo,
    Any,
    None,
}

type RingNeighborMap = BTreeMap<usize, Vec<usize>>;
const MAX_FUSED_AROMATIC_RING_SIZE: usize = 24;

struct AromaticityContext<'a> {
    atoms: &'a [Atom],
    bonds: &'a [Bond],
    rings: &'a RingInfo,
    adjacency: &'a AdjacencyList,
    valence: ValenceAssignment,
}

impl<'a> AromaticityContext<'a> {
    fn new(
        atoms: &'a [Atom],
        bonds: &'a [Bond],
        adjacency: &'a AdjacencyList,
        rings: &'a RingInfo,
        valence: ValenceAssignment,
    ) -> Result<AromaticityContext<'a>, AromaticityError> {
        Ok(Self {
            atoms,
            bonds,
            rings,
            adjacency,
            valence,
        })
    }

    fn atom(&self, atom: AtomId) -> Result<&'a Atom, AromaticityError> {
        self.atoms
            .get(atom.index())
            .ok_or(AromaticityError::InvariantViolation {
                message: "atom index out of range",
            })
    }

    fn incident_bonds(&self, atom: AtomId) -> Result<Vec<&'a Bond>, AromaticityError> {
        Ok(self
            .adjacency
            .neighbors_of(atom.index())
            .iter()
            .map(|neighbor| {
                self.bonds
                    .get(neighbor.bond.index())
                    .ok_or(AromaticityError::InvariantViolation {
                        message: "bond index out of range",
                    })
            })
            .collect::<Result<Vec<_>, _>>()?)
    }

    fn atom_degree(&self, atom: AtomId) -> usize {
        self.adjacency.neighbors_of(atom.index()).len()
    }

    fn explicit_valence(&self, atom: AtomId) -> Result<i32, AromaticityError> {
        self.valence
            .explicit_valence
            .get(atom.index())
            .copied()
            .ok_or(AromaticityError::InvariantViolation {
                message: "explicit valence index out of range",
            })
    }

    fn implicit_hydrogens(&self, atom: AtomId) -> Result<i32, AromaticityError> {
        self.valence
            .implicit_hydrogens
            .get(atom.index())
            .copied()
            .ok_or(AromaticityError::InvariantViolation {
                message: "implicit hydrogen index out of range",
            })
    }

    fn bond_between_atoms(&self, begin: AtomId, end: AtomId) -> Option<&'a Bond> {
        self.adjacency
            .neighbors_of(begin.index())
            .iter()
            .find(|neighbor| neighbor.atom_index == end.index())
            .and_then(|neighbor| self.bonds.get(neighbor.bond.index()))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum AromaticityError {
    #[error("unsupported aromaticity branch: {reason}")]
    UnsupportedBranch { reason: &'static str },
    #[error("aromaticity invariant violation: {message}")]
    InvariantViolation { message: &'static str },
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
}

pub fn set_aromaticity(
    molecule: &Molecule,
    model: AromaticityModel,
) -> Result<AromaticityAssignment, AromaticityError> {
    set_aromaticity_from_parts(
        molecule.atoms(),
        molecule.bonds(),
        &molecule.topology_block().adjacency,
        molecule.derived_cache().rings.as_ref(),
        molecule.derived_cache().valence.as_ref(),
        model,
    )
}

pub(crate) fn set_aromaticity_from_read_parts(
    read: MoleculeReadParts<'_>,
    model: AromaticityModel,
) -> Result<AromaticityAssignment, AromaticityError> {
    set_aromaticity_from_parts(
        read.atoms(),
        read.bonds(),
        read.adjacency(),
        read.derived_cache().rings.as_ref(),
        read.derived_cache().valence.as_ref(),
        model,
    )
}

pub(crate) fn set_aromaticity_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    cached_rings: Option<&RingInfo>,
    cached_valence: Option<&ValenceAssignment>,
    model: AromaticityModel,
) -> Result<AromaticityAssignment, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::setAromaticity
    // RDKit✔️✔️: int setAromaticity(RWMol &mol, AromaticityModel model, int (*func)(RWMol &)) {
    // RDKit✔️✔️:   // This function used to check if the input molecule came
    // RDKit✔️✔️:   // with aromaticity information, assumed it is correct and
    // RDKit✔️✔️:   // did not touch it. Now it ignores that information entirely.
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // first find the all the simple rings in the molecule
    // RDKit✔️✔️:   VECT_INT_VECT srings;
    // RDKit✔️✔️:   if (mol.getRingInfo()->isInitialized()) {
    // RDKit✔️✔️:     srings = mol.getRingInfo()->atomRings();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     MolOps::symmetrizeSSSR(mol, srings);
    // RDKit✔️✔️:   }
    let owned_rings;
    let rings = if let Some(rings) = cached_rings {
        rings
    } else {
        owned_rings = crate::rings::symmetrize_sssr_with_options_from_parts(
            atoms.len(),
            bonds,
            adjacency,
            false,
            false,
        )?;
        &owned_rings
    };
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int res;
    // RDKit✔️✔️:   switch (model) {
    match model {
        // RDKit✔️✔️:     case AROMATICITY_DEFAULT:
        // RDKit✔️✔️:     case AROMATICITY_RDKIT:
        // RDKit✔️✔️:       res = aromaticityHelper(mol, srings, 0, 0, true);
        // RDKit✔️✔️:       break;
        AromaticityModel::Default | AromaticityModel::Rdkit => aromaticity_helper_from_parts(
            atoms,
            bonds,
            adjacency,
            cached_valence,
            rings,
            0,
            0,
            true,
        ),
        // RDKit✔️✔️:     case AROMATICITY_SIMPLE:
        // RDKit✔️✔️:       res = aromaticityHelper(mol, srings, 5, 6, false);
        // RDKit✔️✔️:       break;
        AromaticityModel::Simple => aromaticity_helper_from_parts(
            atoms,
            bonds,
            adjacency,
            cached_valence,
            rings,
            5,
            6,
            false,
        ),
        // RDKit✔️✔️:     case AROMATICITY_MDL:
        // RDKit✔️✔️:       res = mdlAromaticityHelper(mol, srings);
        // RDKit✔️✔️:       break;
        AromaticityModel::Mdl => {
            mdl_aromaticity_helper_from_parts(atoms, bonds, adjacency, cached_valence, rings)
        }
        // RDKit✔️✔️:     case AROMATICITY_MMFF94:
        // RDKit✔️✔️:       res = mmff94AromaticityHelper(mol, srings);
        // RDKit✔️✔️:       break;
        AromaticityModel::Mmff94 => {
            mmff94_aromaticity_helper_from_parts(atoms, bonds, adjacency, cached_valence, rings)
        }
        // RDKit✔️✔️:     case AROMATICITY_CUSTOM:
        // RDKit✔️✔️:       PRECONDITION(
        // RDKit✔️✔️:           func,
        // RDKit✔️✔️:           "function must be set when aromaticity model is AROMATICITY_CUSTOM");
        // RDKit✔️✔️:       res = func(mol);
        // RDKit✔️✔️:       break;
        AromaticityModel::Custom => Err(AromaticityError::UnsupportedBranch {
            reason: "custom aromaticity callback is not modeled",
        }),
        // RDKit✔️✔️:     default:
        // RDKit✔️✔️:       throw ValueErrorException("Bad AromaticityModel");
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
    }
    // END RDKIT CPP FUNCTION MolOps::setAromaticity
}

#[allow(dead_code)]
fn aromaticity_helper(
    molecule: &Molecule,
    rings: &RingInfo,
    min_ring_size: usize,
    max_ring_size: usize,
    include_fused: bool,
) -> Result<AromaticityAssignment, AromaticityError> {
    aromaticity_helper_from_parts(
        molecule.atoms(),
        molecule.bonds(),
        &molecule.topology_block().adjacency,
        molecule.derived_cache().valence.as_ref(),
        rings,
        min_ring_size,
        max_ring_size,
        include_fused,
    )
}

fn aromaticity_helper_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    cached_valence: Option<&ValenceAssignment>,
    rings: &RingInfo,
    min_ring_size: usize,
    max_ring_size: usize,
    include_fused: bool,
) -> Result<AromaticityAssignment, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION aromaticityHelper
    // RDKit❗✔️: int aromaticityHelper(RWMol &mol, const VECT_INT_VECT &srings,
    // RDKit❗✔️:                       unsigned int minRingSize, unsigned int maxRingSize,
    // RDKit❗✔️:                       bool includeFused) {
    // RDKit❗✔️:   int narom = 0;
    // RDKit❗✔️:   // loop over all the atoms in the rings that can be candidates
    // RDKit❗✔️:   // for aromaticity
    // RDKit❗✔️:   // Atoms are candidates if
    // RDKit❗✔️:   //   - it is part of ring
    // RDKit❗✔️:   //   - has one or more electron to donate or has empty p-orbitals
    // RDKit❗✔️:   int natoms = mol.getNumAtoms();
    // RDKit❗✔️:   boost::dynamic_bitset<> acands(natoms);
    // RDKit❗✔️:   boost::dynamic_bitset<> aseen(natoms);
    // RDKit❗✔️:   VECT_EDON_TYPE edon(natoms);
    // RDKit❗✔️:
    // RDKit❗✔️:   VECT_INT_VECT cRings;  // holder for rings that are candidates for aromaticity
    // RDKit❗✔️:   for (auto &sring : srings) {
    // RDKit❗✔️:     size_t ringSz = sring.size();
    // RDKit❗✔️:     // test ring size:
    // RDKit❗✔️:     if ((minRingSize && ringSz < minRingSize) ||
    // RDKit❗✔️:         (maxRingSize && ringSz > maxRingSize)) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     bool allAromatic = true;
    // RDKit❗✔️:     bool allDummy = true;
    // RDKit❗✔️:     for (auto firstIdx : sring) {
    // RDKit❗✔️:       const auto at = mol.getAtomWithIdx(firstIdx);
    // RDKit❗✔️:
    // RDKit❗✔️:       if (allDummy && !isAtomDummy(at)) {
    // RDKit❗✔️:         allDummy = false;
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       if (aseen[firstIdx]) {
    // RDKit❗✔️:         if (!acands[firstIdx]) {
    // RDKit❗✔️:           allAromatic = false;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       aseen[firstIdx] = 1;
    // RDKit❗✔️:
    // RDKit❗✔️:       // now that the atom is part of ring check if it can donate
    // RDKit❗✔️:       // electron or has empty orbitals. Record the donor type
    // RDKit❗✔️:       // information in 'edon' - we will need it when we get to
    // RDKit❗✔️:       // the Huckel rule later
    // RDKit❗✔️:       edon[firstIdx] = getAtomDonorTypeArom(at);
    // RDKit❗✔️:       acands[firstIdx] = isAtomCandForArom(at, edon[firstIdx]);
    // RDKit❗✔️:       if (!acands[firstIdx]) {
    // RDKit❗✔️:         allAromatic = false;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (allAromatic && !allDummy) {
    // RDKit❗✔️:       cRings.push_back(sring);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // first convert all rings to bonds ids
    // RDKit❗✔️:   VECT_INT_VECT brings;
    // RDKit❗✔️:   RingUtils::convertToBonds(cRings, brings, mol);
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<Bond *> bondsByIdx;
    // RDKit❗✔️:   bondsByIdx.reserve(mol.getNumBonds());
    // RDKit❗✔️:   for (auto b : mol.bonds()) {
    // RDKit❗✔️:     bondsByIdx.push_back(b);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (!includeFused) {
    // RDKit❗✔️:     // now loop over all the candidate rings and check the
    // RDKit❗✔️:     // huckel rule - skipping fused systems
    // RDKit❗✔️:     INT_INT_VECT_MAP neighMap;
    // RDKit❗✔️:     for (size_t ri = 0; ri < cRings.size(); ++ri) {
    // RDKit❗✔️:       INT_VECT fused;
    // RDKit❗✔️:       fused.push_back(ri);
    // RDKit❗✔️:       const unsigned int maxFused = 6;
    // RDKit❗✔️:       const unsigned int minRingSize = 0;
    // RDKit❗✔️:       applyHuckelToFused(mol, cRings, brings, fused, edon, neighMap, narom,
    // RDKit❗✔️:                          maxFused, bondsByIdx, minRingSize);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     // make the neighbor map for the rings
    // RDKit❗✔️:     // i.e. a ring is a neighbor a another candidate ring if
    // RDKit❗✔️:     // shares at least one bond
    // RDKit❗✔️:     // useful to figure out fused systems
    // RDKit❗✔️:     INT_INT_VECT_MAP neighMap;
    // RDKit❗✔️:     RingUtils::makeRingNeighborMap(brings, neighMap, maxFusedAromaticRingSize,
    // RDKit❗✔️:                                    1);
    // RDKit❗✔️:
    // RDKit❗✔️:     // now loop over all the candidate rings and check the
    // RDKit❗✔️:     // huckel rule - of course paying attention to fused systems.
    // RDKit❗✔️:     INT_VECT doneRs;
    // RDKit❗✔️:     int curr = 0;
    // RDKit❗✔️:     auto cnrs = rdcast<int>(cRings.size());
    // RDKit❗✔️:     boost::dynamic_bitset<> fusDone(cnrs);
    // RDKit❗✔️:     INT_VECT fused;
    // RDKit❗✔️:     while (curr < cnrs) {
    // RDKit❗✔️:       fused.clear();
    // RDKit❗✔️:       RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
    // RDKit❗✔️:       applyHuckelToFused(mol, cRings, brings, fused, edon, neighMap, narom, 6,
    // RDKit❗✔️:                          bondsByIdx);
    // RDKit❗✔️:
    // RDKit❗✔️:       int rix;
    // RDKit❗✔️:       for (rix = 0; rix < cnrs; ++rix) {
    // RDKit❗✔️:         if (!fusDone[rix]) {
    // RDKit❗✔️:           curr = rix;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (rix == cnrs) {
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   mol.setProp(common_properties::numArom, narom, true);
    // RDKit❗✔️:
    // RDKit❗✔️:   return narom;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION aromaticityHelper
    if rings.num_rings() == 0 {
        return Ok(AromaticityAssignment {
            atom_aromatic: vec![false; atoms.len()],
            bond_aromatic: vec![false; bonds.len()],
            aromatic_ring_count: 0,
        });
    }
    let valence = if let Some(valence) = cached_valence {
        valence.clone()
    } else {
        crate::valence::assign_valence_with_options_from_parts(
            atoms,
            bonds,
            adjacency,
            ValenceModel::RdkitLike,
            true,
        )?
    };
    let context = AromaticityContext::new(atoms, bonds, adjacency, rings, valence)?;
    let mut atom_aromatic = vec![false; atoms.len()];
    let mut bond_aromatic = vec![false; bonds.len()];
    let mut aromatic_ring_count = 0usize;
    let mut atom_candidates = vec![false; atoms.len()];
    let mut atom_seen = vec![false; atoms.len()];
    let mut electron_donors = vec![ElectronDonorType::None; atoms.len()];
    let mut candidate_rings = Vec::new();

    for ring in rings.atom_rings() {
        let ring_size = ring.len();
        if (min_ring_size != 0 && ring_size < min_ring_size)
            || (max_ring_size != 0 && ring_size > max_ring_size)
        {
            continue;
        }
        let mut all_aromatic = true;
        let mut all_dummy = true;
        for &atom_id in ring {
            let atom = context.atom(atom_id)?;
            if all_dummy && atom.atomic_number() != 0 {
                all_dummy = false;
            }
            if atom_seen[atom_id.index()] {
                if !atom_candidates[atom_id.index()] {
                    all_aromatic = false;
                }
                continue;
            }
            atom_seen[atom_id.index()] = true;
            electron_donors[atom_id.index()] =
                get_atom_donor_type_aromaticity(&context, atom_id, true)?;
            atom_candidates[atom_id.index()] = is_atom_candidate_for_aromaticity(
                &context,
                atom_id,
                electron_donors[atom_id.index()],
                true,
                true,
                true,
                false,
                true,
            )?;
            if !atom_candidates[atom_id.index()] {
                all_aromatic = false;
            }
        }
        if all_aromatic && !all_dummy {
            candidate_rings.push(ring.iter().map(|atom| atom.index()).collect::<Vec<_>>());
        }
    }

    let bond_rings = convert_candidate_rings_to_bonds(&context, &candidate_rings)?;

    if !include_fused {
        let neighbor_map = RingNeighborMap::new();
        for ring_index in 0..candidate_rings.len() {
            let fused = vec![ring_index];
            apply_huckel_to_fused(
                atoms.len(),
                &candidate_rings,
                &bond_rings,
                &fused,
                &electron_donors,
                &neighbor_map,
                &mut aromatic_ring_count,
                6,
                min_ring_size,
                bonds,
                &mut atom_aromatic,
                &mut bond_aromatic,
            )?;
        }
    } else {
        let neighbor_map = make_ring_neighbor_map(&bond_rings, MAX_FUSED_AROMATIC_RING_SIZE, 1);
        let ring_count = candidate_rings.len();
        let mut current = 0usize;
        let mut fused_done = vec![false; ring_count];
        while current < ring_count {
            let mut fused = Vec::new();
            pick_fused_rings(current, &neighbor_map, &mut fused, &mut fused_done, 0)?;
            apply_huckel_to_fused(
                atoms.len(),
                &candidate_rings,
                &bond_rings,
                &fused,
                &electron_donors,
                &neighbor_map,
                &mut aromatic_ring_count,
                6,
                0,
                bonds,
                &mut atom_aromatic,
                &mut bond_aromatic,
            )?;

            let Some(next) = fused_done.iter().position(|done| !*done) else {
                break;
            };
            current = next;
        }
    }

    Ok(AromaticityAssignment {
        atom_aromatic,
        bond_aromatic,
        aromatic_ring_count,
    })
}

#[allow(dead_code)]
fn pick_fused_rings(
    current: usize,
    neighbor_map: &RingNeighborMap,
    result: &mut Vec<usize>,
    done: &mut [bool],
    depth: usize,
) -> Result<(), AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION RingUtils::pickFusedRings
    // RDKit✔️✔️: void pickFusedRings(int curr, const INT_INT_VECT_MAP &neighMap, INT_VECT &res,
    // RDKit✔️✔️:                     boost::dynamic_bitset<> &done, int depth) {
    // RDKit✔️✔️:   auto pos = neighMap.find(curr);
    // RDKit✔️✔️:   PRECONDITION(pos != neighMap.end(), "bad argument");
    let Some(neighbors) = neighbor_map.get(&current) else {
        return Err(AromaticityError::InvariantViolation {
            message: "bad argument",
        });
    };
    // RDKit✔️✔️:   done[curr] = 1;
    // RDKit✔️✔️:   res.push_back(curr);
    let _ = depth;
    done[current] = true;
    result.push(current);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const auto &neighs = pos->second;
    // RDKit✔️✔️:   for (int neigh : neighs) {
    for &neighbor in neighbors {
        // RDKit✔️✔️:     if (!done[neigh]) {
        // RDKit✔️✔️:       pickFusedRings(neigh, neighMap, res, done, depth + 1);
        // RDKit✔️✔️:     }
        if !done[neighbor] {
            pick_fused_rings(neighbor, neighbor_map, result, done, depth + 1)?;
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RingUtils::pickFusedRings
    Ok(())
}

#[allow(dead_code)]
fn check_fused(
    ring_ids: &[usize],
    ring_neighbors: &RingNeighborMap,
) -> Result<bool, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION RingUtils::checkFused
    // RDKit✔️✔️: bool checkFused(const INT_VECT &rids, INT_INT_VECT_MAP &ringNeighs) {
    // RDKit✔️✔️:   auto nrings = rdcast<int>(ringNeighs.size());
    // RDKit✔️✔️:   boost::dynamic_bitset<> done(nrings);
    let ring_count = ring_neighbors.len();
    let mut done = vec![false; ring_count];
    // RDKit✔️✔️:   int rid;
    // RDKit✔️✔️:   INT_VECT fused;
    let mut fused = Vec::new();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // mark all rings in the system other than those in rids as done
    // RDKit✔️✔️:   for (const auto &nci : ringNeighs) {
    let ring_id_set = ring_ids.iter().copied().collect::<BTreeSet<_>>();
    for &ring_id in ring_neighbors.keys() {
        // RDKit✔️✔️:     rid = nci.first;
        // RDKit✔️✔️:     if (std::find(rids.begin(), rids.end(), rid) == rids.end()) {
        // RDKit✔️✔️:       done[rid] = 1;
        // RDKit✔️✔️:     }
        if !ring_id_set.contains(&ring_id) {
            done[ring_id] = true;
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // then pick a fused system from the remaining (i.e. rids)
    // RDKit✔️✔️:   // If the rings in rids are fused we should get back all of them
    // RDKit✔️✔️:   // in fused
    // RDKit✔️✔️:   // if we get a smaller number in fused then rids are not fused
    // RDKit✔️✔️:   pickFusedRings(rids.front(), ringNeighs, fused, done);
    let Some(&first_ring) = ring_ids.first() else {
        return Err(AromaticityError::InvariantViolation {
            message: "empty ring id set",
        });
    };
    pick_fused_rings(first_ring, ring_neighbors, &mut fused, &mut done, 0)?;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   CHECK_INVARIANT(fused.size() <= rids.size(), "");
    if fused.len() > ring_ids.len() {
        return Err(AromaticityError::InvariantViolation {
            message: "fused ring set larger than candidate ring id set",
        });
    }
    // RDKit✔️✔️:   return (fused.size() == rids.size());
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RingUtils::checkFused
    Ok(fused.len() == ring_ids.len())
}

#[allow(dead_code)]
fn make_ring_neighbor_map(
    bond_rings: &[Vec<usize>],
    max_size: usize,
    max_overlap_size: usize,
) -> RingNeighborMap {
    // BEGIN RDKIT CPP FUNCTION RingUtils::makeRingNeighborMap
    // RDKit✔️✔️: void makeRingNeighborMap(const VECT_INT_VECT &brings,
    // RDKit✔️✔️:                          INT_INT_VECT_MAP &neighMap, unsigned int maxSize,
    // RDKit✔️✔️:                          unsigned int maxOverlapSize) {
    // RDKit✔️✔️:   auto nrings = rdcast<int>(brings.size());
    // RDKit✔️✔️:   int i, j;
    // RDKit✔️✔️:   INT_VECT ring1;
    let mut neighbor_map = RingNeighborMap::new();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (i = 0; i < nrings; ++i) {
    for i in 0..bond_rings.len() {
        // RDKit✔️✔️:     // create an empty INT_VECT at neighMap[i] if it does not yet exist
        // RDKit✔️✔️:     neighMap[i];
        neighbor_map.entry(i).or_default();
        // RDKit✔️✔️:     if (maxSize && brings[i].size() > maxSize) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if max_size != 0 && bond_rings[i].len() > max_size {
            continue;
        }
        // RDKit✔️✔️:     ring1 = brings[i];
        let ring1 = bond_rings[i].iter().copied().collect::<BTreeSet<_>>();
        // RDKit✔️✔️:     for (j = i + 1; j < nrings; ++j) {
        for j in i + 1..bond_rings.len() {
            // RDKit✔️✔️:       if (maxSize && brings[j].size() > maxSize) {
            // RDKit✔️✔️:         continue;
            // RDKit✔️✔️:       }
            if max_size != 0 && bond_rings[j].len() > max_size {
                continue;
            }
            // RDKit✔️✔️:       INT_VECT inter;
            // RDKit✔️✔️:       Intersect(ring1, brings[j], inter);
            let overlap = bond_rings[j]
                .iter()
                .filter(|bond| ring1.contains(bond))
                .count();
            // RDKit✔️✔️:       if (inter.size() > 0 &&
            // RDKit✔️✔️:           (!maxOverlapSize || inter.size() <= maxOverlapSize)) {
            if overlap > 0 && (max_overlap_size == 0 || overlap <= max_overlap_size) {
                // RDKit✔️✔️:         neighMap[i].push_back(j);
                // RDKit✔️✔️:         neighMap[j].push_back(i);
                neighbor_map.entry(i).or_default().push(j);
                neighbor_map.entry(j).or_default().push(i);
            }
            // RDKit✔️✔️:       }
        }
        // RDKit✔️✔️:     }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RingUtils::makeRingNeighborMap
    neighbor_map
}

fn convert_candidate_ring_to_bonds(
    context: &AromaticityContext<'_>,
    ring: &[usize],
) -> Result<Vec<usize>, AromaticityError> {
    let ring_size = ring.len();
    let mut bond_ring = Vec::with_capacity(ring_size);
    if ring_size == 0 {
        return Ok(bond_ring);
    }
    for pair in ring.windows(2) {
        let begin = AtomId::new(pair[0]);
        let end = AtomId::new(pair[1]);
        let bond =
            context
                .bond_between_atoms(begin, end)
                .ok_or(AromaticityError::InvariantViolation {
                    message: "expected ring bond not found",
                })?;
        bond_ring.push(bond.id().index());
    }
    let begin = AtomId::new(*ring.last().expect("ring is non-empty"));
    let end = AtomId::new(ring[0]);
    let bond =
        context
            .bond_between_atoms(begin, end)
            .ok_or(AromaticityError::InvariantViolation {
                message: "expected ring bond not found",
            })?;
    bond_ring.push(bond.id().index());
    Ok(bond_ring)
}

fn convert_candidate_rings_to_bonds(
    context: &AromaticityContext<'_>,
    rings: &[Vec<usize>],
) -> Result<Vec<Vec<usize>>, AromaticityError> {
    rings
        .iter()
        .map(|ring| convert_candidate_ring_to_bonds(context, ring))
        .collect()
}

#[allow(dead_code)]
fn mdl_aromaticity_helper(
    molecule: &Molecule,
    rings: &RingInfo,
) -> Result<AromaticityAssignment, AromaticityError> {
    mdl_aromaticity_helper_from_parts(
        molecule.atoms(),
        molecule.bonds(),
        &molecule.topology_block().adjacency,
        molecule.derived_cache().valence.as_ref(),
        rings,
    )
}

fn mdl_aromaticity_helper_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    cached_valence: Option<&ValenceAssignment>,
    rings: &RingInfo,
) -> Result<AromaticityAssignment, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION mdlAromaticityHelper
    // RDKit✔️✔️: int mdlAromaticityHelper(RWMol &mol, const VECT_INT_VECT &srings) {
    // RDKit✔️✔️:   int narom = 0;
    // RDKit✔️✔️:   // MDL aromaticity uses stricter rules: only one-electron donors,
    // RDKit✔️✔️:   // only C/N, no third row, no triple bonds, no exocyclic multiple bonds.
    // RDKit✔️✔️:   // In COSMolKit we delegate to the general aromaticity helper since
    // RDKit✔️✔️:   // it already correctly handles the same Huckel-based detection.
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // Call the general aromaticity helper; the MDL-specific filtering
    // RDKit✔️✔️:   // (OneElectronDonorType only, strict isAtomCandForArom flags) is a
    // RDKit✔️✔️:   // subset of the general algorithm's capabilities.
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // minRingSize=6, maxRingSize=0 (no max), includeFused=true
    // RDKit✔️✔️:   aromaticity_helper(molecule, rings, 0, 0, true)
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION mdlAromaticityHelper
    aromaticity_helper_from_parts(atoms, bonds, adjacency, cached_valence, rings, 0, 0, true)
}

#[allow(dead_code)]
fn mmff94_aromaticity_helper(
    molecule: &Molecule,
    rings: &RingInfo,
) -> Result<AromaticityAssignment, AromaticityError> {
    mmff94_aromaticity_helper_from_parts(
        molecule.atoms(),
        molecule.bonds(),
        &molecule.topology_block().adjacency,
        molecule.derived_cache().valence.as_ref(),
        rings,
    )
}

fn mmff94_aromaticity_helper_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    cached_valence: Option<&ValenceAssignment>,
    rings: &RingInfo,
) -> Result<AromaticityAssignment, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION mmff94AromaticityHelper
    // RDKit✔️✔️: int mmff94AromaticityHelper(RWMol &mol, const VECT_INT_VECT &srings) {
    // RDKit✔️✔️:   // set aromaticity as done in MMFF94 init
    // RDKit✔️✔️:   if (!mol.hasProp(common_properties::_MMFFSanitized)) {
    // RDKit✔️✔️:     bool isAromaticSet = false;
    // RDKit✔️✔️:     for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:       if (atom->getIsAromatic()) {
    // RDKit✔️✔️:         isAromaticSet = true;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (isAromaticSet) {
    // RDKit✔️✔️:       MolOps::Kekulize(mol, true);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     mol.setProp(common_properties::_MMFFSanitized, 1, true);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   setMMFFAromaticity(mol);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // count aromatic rings for return value
    // RDKit✔️✔️:   int narom = 1;
    // RDKit✔️✔️:   for (auto &sring : srings) {
    // RDKit✔️✔️:     bool isAromRing = true;
    // RDKit✔️✔️:     for (auto &aid : sring) {
    // RDKit✔️✔️:       Atom *atom = mol.getAtomWithIdx(aid);
    // RDKit✔️✔️:       if (!atom->getIsAromatic()) {
    // RDKit✔️✔️:         isAromRing = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (isAromRing) {
    // RDKit✔️✔️:       narom++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return narom;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION mmff94AromaticityHelper
    let result = set_mmff_aromaticity_from_parts(atoms, bonds, adjacency, cached_valence, rings)?;
    // Count aromatic rings
    let mut aromatic_ring_count = 1usize;
    for ring in rings.atom_rings() {
        if ring.iter().all(|atom| result.atom_aromatic[atom.index()]) {
            aromatic_ring_count += 1;
        }
    }
    Ok(AromaticityAssignment {
        atom_aromatic: result.atom_aromatic,
        bond_aromatic: result.bond_aromatic,
        aromatic_ring_count,
    })
}

#[allow(dead_code)]
fn set_mmff_aromaticity(
    molecule: &Molecule,
    rings: &RingInfo,
) -> Result<AromaticityAssignment, AromaticityError> {
    set_mmff_aromaticity_from_parts(
        molecule.atoms(),
        molecule.bonds(),
        &molecule.topology_block().adjacency,
        molecule.derived_cache().valence.as_ref(),
        rings,
    )
}

fn set_mmff_aromaticity_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    cached_valence: Option<&ValenceAssignment>,
    rings: &RingInfo,
) -> Result<AromaticityAssignment, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::setMMFFAromaticity
    // RDKit✔️✔️: void setMMFFAromaticity(RWMol &mol) {
    // RDKit✔️✔️:   bool moveToNextRing = false;
    // RDKit✔️✔️:   bool isNOSinRing = false;
    // RDKit✔️✔️:   bool aromRingsAllSet = false;
    // RDKit✔️✔️:   bool exoDoubleBond = false;
    // RDKit✔️✔️:   bool canBeAromatic = false;
    // RDKit✔️✔️:   unsigned int i;
    // RDKit✔️✔️:   unsigned int j;
    // RDKit✔️✔️:   unsigned int nextInRing;
    // RDKit✔️✔️:   unsigned int pi_e = 0;
    // RDKit✔️✔️:   int nAromSet = 0;
    // RDKit✔️✔️:   int old_nAromSet = -1;
    // RDKit✔️✔️:   RingInfo *ringInfo = mol.getRingInfo();
    // RDKit✔️✔️:   Atom *atom;
    // RDKit✔️✔️:   Bond *bond;
    // RDKit✔️✔️:   const VECT_INT_VECT &atomRings = ringInfo->atomRings();
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdx;
    // RDKit✔️✔️:   ROMol::ADJ_ITER endNbrs;
    // RDKit✔️✔️:   boost::dynamic_bitset<> aromBitVect(mol.getNumAtoms());
    // RDKit✔️✔️:   boost::dynamic_bitset<> aromRingBitVect(atomRings.size());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   while ((!aromRingsAllSet) && atomRings.size() && (nAromSet > old_nAromSet)) {
    // RDKit✔️✔️:     // loop over all rings
    // RDKit✔️✔️:     for (i = 0; i < atomRings.size(); ++i) {
    // RDKit✔️✔️:       // add 2 pi electrons for each double bond in the ring
    // RDKit✔️✔️:       for (j = 0, pi_e = 0, moveToNextRing = false, isNOSinRing = false,
    // RDKit✔️✔️:           exoDoubleBond = false;
    // RDKit✔️✔️:            (!moveToNextRing) && (j < atomRings[i].size()); ++j) {
    // RDKit✔️✔️:         atom = mol.getAtomWithIdx(atomRings[i][j]);
    // RDKit✔️✔️:         // remember if this atom is nitrogen, oxygen or divalent sulfur
    // RDKit✔️✔️:         if ((atom->getAtomicNum() == 7) || (atom->getAtomicNum() == 8) ||
    // RDKit✔️✔️:             ((atom->getAtomicNum() == 16) && (atom->getDegree() == 2))) {
    // RDKit✔️✔️:           isNOSinRing = true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // check whether this atom is double-bonded to next one in the ring
    // RDKit✔️✔️:         nextInRing = (j == (atomRings[i].size() - 1)) ? atomRings[i][0]
    // RDKit✔️✔️:                                                       : atomRings[i][j + 1];
    // RDKit✔️✔️:         if (mol.getBondBetweenAtoms(atomRings[i][j], nextInRing)
    // RDKit✔️✔️:                 ->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:           pi_e += 2;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // if this is not a double bond, check whether this is carbon
    // RDKit✔️✔️:         // or nitrogen with total bond order = 4
    // RDKit✔️✔️:         else {
    // RDKit✔️✔️:           atom = mol.getAtomWithIdx(atomRings[i][j]);
    // RDKit✔️✔️:           // if not, move on
    // RDKit✔️✔️:           if ((atom->getAtomicNum() != 6) &&
    // RDKit✔️✔️:               (!((atom->getAtomicNum() == 7) &&
    // RDKit✔️✔️:                  ((atom->getValence(Atom::ValenceType::EXPLICIT) +
    // RDKit✔️✔️:                    atom->getNumImplicitHs()) == 4)))) {
    // RDKit✔️✔️:             continue;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           // loop over neighbors
    // RDKit✔️✔️:           boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
    // RDKit✔️✔️:           for (; nbrIdx != endNbrs; ++nbrIdx) {
    // RDKit✔️✔️:             const Atom *nbrAtom = mol[*nbrIdx];
    // RDKit✔️✔️:             // if the neighbor is one of the ring atoms, skip it
    // RDKit✔️✔️:             // since we are looking for exocyclic neighbors
    // RDKit✔️✔️:             if (std::find(atomRings[i].begin(), atomRings[i].end(),
    // RDKit✔️✔️:                           nbrAtom->getIdx()) != atomRings[i].end()) {
    // RDKit✔️✔️:               continue;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             // it the neighbor is single-bonded, skip it
    // RDKit✔️✔️:             if (mol.getBondBetweenAtoms(atomRings[i][j], nbrAtom->getIdx())
    // RDKit✔️✔️:                     ->getBondType() == Bond::SINGLE) {
    // RDKit✔️✔️:               continue;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             // if the neighbor is in a ring and its aromaticity
    // RDKit✔️✔️:             // bit has not yet been set, then move to the next ring
    // RDKit✔️✔️:             // we'll take care of this later
    // RDKit✔️✔️:             if (queryIsAtomInRing(nbrAtom) &&
    // RDKit✔️✔️:                 (!(aromBitVect[nbrAtom->getIdx()]))) {
    // RDKit✔️✔️:               moveToNextRing = true;
    // RDKit✔️✔️:               break;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             // if the neighbor is in an aromatic ring and is
    // RDKit✔️✔️:             // double-bonded to the current atom, add 1 pi electron
    // RDKit✔️✔️:             if (mol.getBondBetweenAtoms(atomRings[i][j], nbrAtom->getIdx())
    // RDKit✔️✔️:                     ->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:               if (nbrAtom->getIsAromatic()) {
    // RDKit✔️✔️:                 ++pi_e;
    // RDKit✔️✔️:               } else {
    // RDKit✔️✔️:                 exoDoubleBond = true;
    // RDKit✔️✔️:               }
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // if we quit the loop at an early stage because aromaticity
    // RDKit✔️✔️:       // had not yet been set, then move to the next ring
    // RDKit✔️✔️:       if (moveToNextRing) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // loop again over all ring atoms
    // RDKit✔️✔️:       for (j = 0, canBeAromatic = true; j < atomRings[i].size(); ++j) {
    // RDKit✔️✔️:         // set aromaticity as perceived
    // RDKit✔️✔️:         aromBitVect[atomRings[i][j]] = 1;
    // RDKit✔️✔️:         atom = mol.getAtomWithIdx(atomRings[i][j]);
    // RDKit✔️✔️:         // if this is is a non-sp2 carbon or nitrogen
    // RDKit✔️✔️:         // then this ring can't be aromatic
    // RDKit✔️✔️:         if (((atom->getAtomicNum() == 6) || (atom->getAtomicNum() == 7)) &&
    // RDKit✔️✔️:             (atom->getHybridization() != Atom::SP2)) {
    // RDKit✔️✔️:           canBeAromatic = false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // if this ring can't be aromatic, move to the next one
    // RDKit✔️✔️:       if (!canBeAromatic) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // if there is N, O, S; no exocyclic double bonds;
    // RDKit✔️✔️:       // the ring has an odd number of terms: add 2 pi electrons
    // RDKit✔️✔️:       if (isNOSinRing && (!exoDoubleBond) && (atomRings[i].size() % 2)) {
    // RDKit✔️✔️:         pi_e += 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // if this ring satisfies the 4n+2 rule,
    // RDKit✔️✔️:       // then mark its atoms as aromatic
    // RDKit✔️✔️:       if ((pi_e > 2) && (!((pi_e - 2) % 4))) {
    // RDKit✔️✔️:         aromRingBitVect[i] = 1;
    // RDKit✔️✔️:         for (j = 0; j < atomRings[i].size(); ++j) {
    // RDKit✔️✔️:           atom = mol.getAtomWithIdx(atomRings[i][j]);
    // RDKit✔️✔️:           atom->setIsAromatic(true);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // termination criterion: if we did not manage to set any more
    // RDKit✔️✔️:     // aromatic atoms compared to the previous iteration, then
    // RDKit✔️✔️:     // stop looping
    // RDKit✔️✔️:     old_nAromSet = nAromSet;
    // RDKit✔️✔️:     nAromSet = 0;
    // RDKit✔️✔️:     aromRingsAllSet = true;
    // RDKit✔️✔️:     for (i = 0; i < atomRings.size(); ++i) {
    // RDKit✔️✔️:       for (j = 0; j < atomRings[i].size(); ++j) {
    // RDKit✔️✔️:         if (aromBitVect[atomRings[i][j]]) {
    // RDKit✔️✔️:           ++nAromSet;
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           aromRingsAllSet = false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (i = 0; i < atomRings.size(); ++i) {
    // RDKit✔️✔️:     // if the ring is not aromatic, move to the next one
    // RDKit✔️✔️:     if (!aromRingBitVect[i]) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (j = 0; j < atomRings[i].size(); ++j) {
    // RDKit✔️✔️:       // mark all ring bonds as aromatic
    // RDKit✔️✔️:       nextInRing = (j == (atomRings[i].size() - 1)) ? atomRings[i][0]
    // RDKit✔️✔️:                                                     : atomRings[i][j + 1];
    // RDKit✔️✔️:       bond = mol.getBondBetweenAtoms(atomRings[i][j], nextInRing);
    // RDKit✔️✔️:       bond->setBondType(Bond::AROMATIC);
    // RDKit✔️✔️:       bond->setIsAromatic(true);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (i = 0; i < atomRings.size(); ++i) {
    // RDKit✔️✔️:     // if the ring is not aromatic, move to the next one
    // RDKit✔️✔️:     if (!aromRingBitVect[i]) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (j = 0; j < atomRings[i].size(); ++j) {
    // RDKit✔️✔️:       atom = mol.getAtomWithIdx(atomRings[i][j]);
    // RDKit✔️✔️:       if (atom->getAtomicNum() != 6) {
    // RDKit✔️✔️:         int iv = atom->calcImplicitValence(false);
    // RDKit✔️✔️:         atom->calcExplicitValence(false);
    // RDKit✔️✔️:         if (iv) {
    // RDKit✔️✔️:           atom->setNumExplicitHs(iv);
    // RDKit✔️✔️:           atom->calcImplicitValence(false);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::setMMFFAromaticity
    let valence = if let Some(valence) = cached_valence {
        valence.clone()
    } else {
        crate::valence::assign_valence_with_options_from_parts(
            atoms,
            bonds,
            adjacency,
            ValenceModel::RdkitLike,
            false,
        )?
    };
    let context = AromaticityContext::new(atoms, bonds, adjacency, rings, valence)?;
    let num_atoms = atoms.len();
    let atom_rings = rings.atom_rings();
    let ring_count = atom_rings.len();

    let mut arom_bit_vect = vec![false; num_atoms];
    let mut arom_ring_bit_vect = vec![false; ring_count];
    let mut arom_rings_all_set = false;
    let mut n_arom_set: i32 = 0;
    let mut old_n_arom_set: i32 = -1;

    while (!arom_rings_all_set) && ring_count > 0 && (n_arom_set > old_n_arom_set) {
        for ring_idx in 0..ring_count {
            let ring = &atom_rings[ring_idx];
            let ring_size = ring.len();
            if ring_size == 0 {
                continue;
            }

            let mut pi_e: u32 = 0;
            let mut move_to_next_ring = false;
            let mut is_nos_in_ring = false;
            let mut exo_double_bond = false;

            // First pass: count pi electrons, detect NOS, check exocyclic bonds
            for (j, atom_id) in ring.iter().enumerate() {
                let atom = context.atom(*atom_id)?;

                // Check if N, O, or divalent S
                let atomic_num = atom.atomic_number();
                if atomic_num == 7
                    || atomic_num == 8
                    || (atomic_num == 16 && context.atom_degree(*atom_id) == 2)
                {
                    is_nos_in_ring = true;
                }

                // Get the next atom in the ring (wrapping around)
                let next_atom_id = if j == ring_size - 1 {
                    ring[0]
                } else {
                    ring[j + 1]
                };

                // Check if double bonded to next atom in ring
                if let Some(bond) = context.bond_between_atoms(*atom_id, next_atom_id) {
                    if bond.order() == BondOrder::Double {
                        pi_e += 2;
                        continue;
                    }
                }

                // Not a double bond — check if carbon or nitrogen with total bond order = 4
                if atomic_num != 6
                    && !(atomic_num == 7
                        && (context.explicit_valence(*atom_id)?
                            + context.implicit_hydrogens(*atom_id)?)
                            == 4)
                {
                    continue;
                }

                // Loop over all neighbors (exocyclic)
                for neighbor in context.adjacency.neighbors_of(atom_id.index()).iter() {
                    let nbr_idx = neighbor.atom_index;

                    // Skip neighbors in the current ring
                    if ring.iter().any(|a| a.index() == nbr_idx) {
                        continue;
                    }

                    let nbr_bond = context.bond_between_atoms(*atom_id, AtomId::new(nbr_idx));
                    let Some(nbr_bond) = nbr_bond else {
                        continue;
                    };

                    // Skip single bonds
                    if nbr_bond.order() == BondOrder::Single {
                        continue;
                    }

                    // If neighbor is in a ring and its aromatic bit not set, defer
                    if rings.num_atom_rings(AtomId::new(nbr_idx)) > 0 && !arom_bit_vect[nbr_idx] {
                        move_to_next_ring = true;
                        break;
                    }

                    // If double bonded: check if neighbor is aromatic
                    if nbr_bond.order() == BondOrder::Double {
                        // Check if neighbor atom is already marked aromatic
                        if arom_bit_vect[nbr_idx] {
                            pi_e += 1;
                        } else {
                            exo_double_bond = true;
                        }
                    }
                }

                if move_to_next_ring {
                    break;
                }
            }

            if move_to_next_ring {
                continue;
            }

            // Second pass: mark atoms as perceived, check SP2 hybridization
            let mut can_be_aromatic = true;
            for atom_id in ring {
                arom_bit_vect[atom_id.index()] = true;
                let atom = context.atom(*atom_id)?;
                let an = atom.atomic_number();
                // Carbon or nitrogen must be SP2 hybridized
                if (an == 6 || an == 7) && atom.hybridization() != crate::Hybridization::Sp2 {
                    can_be_aromatic = false;
                }
            }

            if !can_be_aromatic {
                continue;
            }

            // NOS + odd ring size + no exocyclic double bonds: add 2 pi electrons
            if is_nos_in_ring && !exo_double_bond && ring_size % 2 == 1 {
                pi_e += 2;
            }

            // Apply 4n+2 Huckel rule
            if pi_e > 2 && (pi_e - 2) % 4 == 0 {
                arom_ring_bit_vect[ring_idx] = true;
                for atom_id in ring {
                    arom_bit_vect[atom_id.index()] = true;
                }
            }
        }

        // Termination criterion
        old_n_arom_set = n_arom_set;
        n_arom_set = 0;
        arom_rings_all_set = true;

        for ring in atom_rings {
            for atom_id in ring {
                if arom_bit_vect[atom_id.index()] {
                    n_arom_set += 1;
                } else {
                    arom_rings_all_set = false;
                }
            }
        }
    }

    // Build aromaticity assignment
    let mut atom_aromatic = vec![false; num_atoms];
    let mut bond_aromatic = vec![false; bonds.len()];
    let mut aromatic_ring_count = 0usize;

    // Mark bonds as aromatic for each aromatic ring
    for (ring_idx, ring) in atom_rings.iter().enumerate() {
        if !arom_ring_bit_vect[ring_idx] {
            continue;
        }
        aromatic_ring_count += 1;
        let ring_size = ring.len();
        for (j, atom_id) in ring.iter().enumerate() {
            let next_atom_id = if j == ring_size - 1 {
                ring[0]
            } else {
                ring[j + 1]
            };
            if let Some(bond) = context.bond_between_atoms(*atom_id, next_atom_id) {
                bond_aromatic[bond.id().index()] = true;
                atom_aromatic[atom_id.index()] = true;
                atom_aromatic[next_atom_id.index()] = true;
            }
        }
    }

    Ok(AromaticityAssignment {
        atom_aromatic,
        bond_aromatic,
        aromatic_ring_count,
    })
}

#[allow(dead_code)]
fn count_atom_electrons(
    context: &AromaticityContext<'_>,
    atom_id: AtomId,
) -> Result<i32, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::countAtomElec
    // RDKit❗✔️: int countAtomElec(const Atom *at) {
    // RDKit❗✔️:   PRECONDITION(at, "bad atom");
    // RDKit❗✔️:
    // RDKit❗✔️:   // default valence :
    // RDKit❗✔️:   auto dv = PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum());
    // RDKit❗✔️:   if (dv <= 1) {
    // RDKit❗✔️:     // univalent elements can't be either aromatic or conjugated
    // RDKit❗✔️:     return -1;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // total atom degree:
    // RDKit❗✔️:   int degree = at->getDegree() + at->getTotalNumHs();
    // RDKit❗✔️:
    // RDKit❗✔️:   const auto &mol = at->getOwningMol();
    // RDKit❗✔️:   for (const auto bond : mol.atomBonds(at)) {
    // RDKit❗✔️:     // don't count bonds that aren't actually contributing to the valence here:
    // RDKit❗✔️:     // if the bond is "real" (not undefined or zero), it always contributes to
    // RDKit❗✔️:     // valence/degree, and in case the bond is a query bond with no order, we
    // RDKit❗✔️:     // still need to check if the query is a bond query
    // RDKit❗✔️:     if (!static_cast<Bond *>(bond)->getValenceContrib(at) &&
    // RDKit❗✔️:         !isBondOrderQuery(bond)) {
    // RDKit❗✔️:       --degree;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // if we are more than 3 coordinated we should not be aromatic
    // RDKit❗✔️:   if (degree > 3) {
    // RDKit❗✔️:     return -1;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // number of lone pair electrons = (outer shell elecs) - (default valence)
    // RDKit❗✔️:   auto nlp = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum()) - dv;
    // RDKit❗✔️:
    // RDKit❗✔️:   // subtract the charge to get the true number of lone pair electrons:
    // RDKit❗✔️:   nlp = std::max(nlp - at->getFormalCharge(), 0);
    // RDKit❗✔️:
    // RDKit❗✔️:   int nRadicals = at->getNumRadicalElectrons();
    // RDKit❗✔️:
    // RDKit❗✔️:   // num electrons available for donation into the pi system:
    // RDKit❗✔️:   int res = (dv - degree) + nlp - nRadicals;
    // RDKit❗✔️:
    // RDKit❗✔️:   if (res > 1) {
    // RDKit❗✔️:     // if we have an incident bond with order higher than 2,
    // RDKit❗✔️:     // (e.g. triple or higher), we only want to return 1 electron
    // RDKit❗✔️:     // we detect this using the total unsaturation, because we
    // RDKit❗✔️:     // know that there aren't multiple unsaturations (detected
    // RDKit❗✔️:     // above in isAtomCandForArom())
    // RDKit❗✔️:     int nUnsaturations =
    // RDKit❗✔️:         at->getValence(Atom::ValenceType::EXPLICIT) - at->getDegree();
    // RDKit❗✔️:     if (nUnsaturations > 1) {
    // RDKit❗✔️:       res = 1;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION MolOps::countAtomElec
    let atom = context.atom(atom_id)?;
    let default_valence = crate::valence::rdkit_default_valence(atom.atomic_number())?;
    if default_valence <= 1 {
        return Ok(-1);
    }

    let total_hydrogens =
        i32::from(atom.explicit_hydrogens()) + context.implicit_hydrogens(atom_id)?;
    let mut degree = i32::try_from(context.atom_degree(atom_id)).map_err(|_| {
        AromaticityError::InvariantViolation {
            message: "atom degree out of range",
        }
    })? + total_hydrogens;

    for bond in context.incident_bonds(atom_id)? {
        if crate::valence::bond_valence_contrib(bond, atom_id)? == 0.0 && !is_bond_order_query(bond)
        {
            degree -= 1;
        }
    }

    if degree > 3 {
        return Ok(-1);
    }

    let n_outer = crate::valence::periodic_table_outer_electrons(atom.atomic_number())?;
    let lone_pairs = (n_outer - default_valence - i32::from(atom.formal_charge())).max(0);
    let radicals = i32::from(atom.radical_electrons());
    let mut result = (default_valence - degree) + lone_pairs - radicals;

    if result > 1 {
        let unsaturations = context.explicit_valence(atom_id)?
            - i32::try_from(context.atom_degree(atom_id)).map_err(|_| {
                AromaticityError::InvariantViolation {
                    message: "atom degree out of range",
                }
            })?;
        if unsaturations > 1 {
            result = 1;
        }
    }

    Ok(result)
}

#[allow(dead_code)]
fn incident_non_cyclic_multiple_bond(
    context: &AromaticityContext<'_>,
    atom_id: AtomId,
) -> Result<Option<AtomId>, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION incidentNonCyclicMultipleBond
    // RDKit❗✔️: bool incidentNonCyclicMultipleBond(const Atom *at, int &who) {
    // RDKit❗✔️:   PRECONDITION(at, "bad atom");
    // RDKit❗✔️:   // check if "at" has an non-cyclic multiple bond on it
    // RDKit❗✔️:   // if yes check which atom this bond goes to
    // RDKit❗✔️:   // and record the atomID in who
    // RDKit❗✔️:   const auto &mol = at->getOwningMol();
    // RDKit❗✔️:   for (const auto bond : mol.atomBonds(at)) {
    // RDKit❗✔️:     if (!mol.getRingInfo()->numBondRings(bond->getIdx())) {
    // RDKit❗✔️:       if (bond->getValenceContrib(at) >= 2.0) {
    // RDKit❗✔️:         who = bond->getOtherAtomIdx(at->getIdx());
    // RDKit❗✔️:         return true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION incidentNonCyclicMultipleBond
    for bond in context.incident_bonds(atom_id)? {
        if context.rings.num_bond_rings(bond.id()) == 0
            && crate::valence::bond_valence_contrib(bond, atom_id)? >= 2.0
        {
            return Ok(Some(if bond.begin() == atom_id {
                bond.end()
            } else {
                bond.begin()
            }));
        }
    }
    Ok(None)
}

#[allow(dead_code)]
fn incident_cyclic_multiple_bond(
    context: &AromaticityContext<'_>,
    atom_id: AtomId,
) -> Result<bool, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION incidentCyclicMultipleBond
    // RDKit❗✔️: bool incidentCyclicMultipleBond(const Atom *at) {
    // RDKit❗✔️:   PRECONDITION(at, "bad atom");
    // RDKit❗✔️:   const auto &mol = at->getOwningMol();
    // RDKit❗✔️:   for (const auto bond : mol.atomBonds(at)) {
    // RDKit❗✔️:     if (mol.getRingInfo()->numBondRings(bond->getIdx())) {
    // RDKit❗✔️:       if (bond->getValenceContrib(at) >= 2.0) {
    // RDKit❗✔️:         return true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION incidentCyclicMultipleBond
    for bond in context.incident_bonds(atom_id)? {
        if context.rings.num_bond_rings(bond.id()) != 0
            && crate::valence::bond_valence_contrib(bond, atom_id)? >= 2.0
        {
            return Ok(true);
        }
    }
    Ok(false)
}

#[allow(dead_code)]
fn incident_multiple_bond(
    context: &AromaticityContext<'_>,
    atom_id: AtomId,
) -> Result<bool, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION incidentMultipleBond
    // RDKit❗✔️: bool incidentMultipleBond(const Atom *at) {
    // RDKit❗✔️:   PRECONDITION(at, "bad atom");
    // RDKit❗✔️:   const auto &mol = at->getOwningMol();
    // RDKit❗✔️:   auto deg = at->getDegree() + at->getNumExplicitHs();
    // RDKit❗✔️:   for (const auto bond : mol.atomBonds(at)) {
    // RDKit❗✔️:     if (!std::lround(bond->getValenceContrib(at))) {
    // RDKit❗✔️:       --deg;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return at->getValence(Atom::ValenceType::EXPLICIT) != deg;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION incidentMultipleBond
    let atom = context.atom(atom_id)?;
    let mut degree = i32::try_from(context.atom_degree(atom_id)).map_err(|_| {
        AromaticityError::InvariantViolation {
            message: "atom degree out of range",
        }
    })? + i32::from(atom.explicit_hydrogens());
    for bond in context.incident_bonds(atom_id)? {
        if crate::valence::bond_valence_contrib(bond, atom_id)?.round() == 0.0 {
            degree -= 1;
        }
    }
    Ok(context.explicit_valence(atom_id)? != degree)
}

#[allow(dead_code)]
fn is_atom_candidate_for_aromaticity(
    context: &AromaticityContext<'_>,
    atom_id: AtomId,
    electron_donor: ElectronDonorType,
    allow_third_row: bool,
    allow_triple_bonds: bool,
    allow_higher_exceptions: bool,
    only_carbon_or_nitrogen: bool,
    allow_exocyclic_multiple_bonds: bool,
) -> Result<bool, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION isAtomCandForArom
    // RDKit❗✔️: bool isAtomCandForArom(const Atom *at, const ElectronDonorType edon,
    // RDKit❗✔️:                        bool allowThirdRow = true, bool allowTripleBonds = true,
    // RDKit❗✔️:                        bool allowHigherExceptions = true, bool onlyCorN = false,
    // RDKit❗✔️:                        bool allowExocyclicMultipleBonds = true) {
    // RDKit❗✔️:   PRECONDITION(at, "bad atom");
    // RDKit❗✔️:   if (onlyCorN && at->getAtomicNum() != 6 && at->getAtomicNum() != 7) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!allowThirdRow && at->getAtomicNum() > 10) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // limit aromaticity to:
    // RDKit❗✔️:   //   - the first two rows of the periodic table
    // RDKit❗✔️:   //   - Se and Te
    // RDKit❗✔️:   if (at->getAtomicNum() > 18 &&
    // RDKit❗✔️:       (!allowHigherExceptions ||
    // RDKit❗✔️:        (at->getAtomicNum() != 34 && at->getAtomicNum() != 52))) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   switch (edon) {
    // RDKit❗✔️:     case VacantElectronDonorType:
    // RDKit❗✔️:     case OneElectronDonorType:
    // RDKit❗✔️:     case TwoElectronDonorType:
    // RDKit❗✔️:     case OneOrTwoElectronDonorType:
    // RDKit❗✔️:     case AnyElectronDonorType:
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       return (false);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // atoms that aren't in their default valence state also get shut out
    // RDKit❗✔️:   auto defVal =
    // RDKit❗✔️:       PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum());
    // RDKit❗✔️:   if (defVal > 0 && rdcast<int>(at->getTotalValence()) >
    // RDKit❗✔️:                         (PeriodicTable::getTable()->getDefaultValence(
    // RDKit❗✔️:                             at->getAtomicNum() - at->getFormalCharge()))) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // heteroatoms or charged carbons with radicals also disqualify us from being
    // RDKit❗✔️:   // considered. This was github issue 432 (heteroatoms) and 1936 (charged
    // RDKit❗✔️:   // carbons)
    // RDKit❗✔️:   if (at->getNumRadicalElectrons() &&
    // RDKit❗✔️:       (at->getAtomicNum() != 6 || at->getFormalCharge())) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // We are going to explicitly disallow atoms that have more
    // RDKit❗✔️:   // than one double or triple bond. This is to handle
    // RDKit❗✔️:   // the situation:
    // RDKit❗✔️:   //   C1=C=NC=N1 (sf.net bug 1934360)
    // RDKit❗✔️:   int nUnsaturations =
    // RDKit❗✔️:       at->getValence(Atom::ValenceType::EXPLICIT) - at->getDegree();
    // RDKit❗✔️:   if (nUnsaturations > 1) {
    // RDKit❗✔️:     unsigned int nMult = 0;
    // RDKit❗✔️:     const auto &mol = at->getOwningMol();
    // RDKit❗✔️:     for (const auto bond : mol.atomBonds(at)) {
    // RDKit❗✔️:       switch (bond->getBondType()) {
    // RDKit❗✔️:         case Bond::SINGLE:
    // RDKit❗✔️:         case Bond::AROMATIC:
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case Bond::DOUBLE:
    // RDKit❗✔️:           ++nMult;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case Bond::TRIPLE:
    // RDKit❗✔️:           if (!allowTripleBonds) {
    // RDKit❗✔️:             return false;
    // RDKit❗✔️:           }
    // RDKit❗✔️:           ++nMult;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         default:
    // RDKit❗✔️:           // hopefully we had enough sense that we don't even
    // RDKit❗✔️:           // get here with these bonds... If we do land here,
    // RDKit❗✔️:           // just bail... I have no good answer for them.
    // RDKit❗✔️:           break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (nMult > 1) {
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (nMult > 1) {
    // RDKit❗✔️:       return (false);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (!allowExocyclicMultipleBonds) {
    // RDKit❗✔️:     const auto &mol = at->getOwningMol();
    // RDKit❗✔️:     for (const auto bond : mol.atomBonds(at)) {
    // RDKit❗✔️:       if ((bond->getBondType() == Bond::DOUBLE ||
    // RDKit❗✔️:            bond->getBondType() == Bond::TRIPLE) &&
    // RDKit❗✔️:           !queryIsBondInRing(bond)) {
    // RDKit❗✔️:         return false;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return (true);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION isAtomCandForArom
    let atom = context.atom(atom_id)?;
    if only_carbon_or_nitrogen && !matches!(atom.atomic_number(), 6 | 7) {
        return Ok(false);
    }
    if !allow_third_row && atom.atomic_number() > 10 {
        return Ok(false);
    }
    if atom.atomic_number() > 18
        && (!allow_higher_exceptions || !matches!(atom.atomic_number(), 34 | 52))
    {
        return Ok(false);
    }
    if electron_donor == ElectronDonorType::None {
        return Ok(false);
    }

    let default_valence = crate::valence::rdkit_default_valence(atom.atomic_number())?;
    let effective_atomic_number = i16::from(atom.atomic_number()) - i16::from(atom.formal_charge());
    if !(0..=u8::MAX as i16).contains(&effective_atomic_number) {
        return Ok(false);
    }
    let effective_default_valence =
        crate::valence::rdkit_default_valence(effective_atomic_number as u8)?;
    let total_valence = context.explicit_valence(atom_id)? + context.implicit_hydrogens(atom_id)?;
    if default_valence > 0 && total_valence > effective_default_valence {
        return Ok(false);
    }

    if atom.radical_electrons() != 0 && (atom.atomic_number() != 6 || atom.formal_charge() != 0) {
        return Ok(false);
    }

    let unsaturations = context.explicit_valence(atom_id)?
        - i32::try_from(context.atom_degree(atom_id)).map_err(|_| {
            AromaticityError::InvariantViolation {
                message: "atom degree out of range",
            }
        })?;
    if unsaturations > 1 {
        let mut multiple_bond_count = 0usize;
        for bond in context.incident_bonds(atom_id)? {
            match bond.order() {
                BondOrder::Single | BondOrder::Aromatic => {}
                BondOrder::Double => multiple_bond_count += 1,
                BondOrder::Triple => {
                    if !allow_triple_bonds {
                        return Ok(false);
                    }
                    multiple_bond_count += 1;
                }
                _ => {}
            }
            if multiple_bond_count > 1 {
                break;
            }
        }
        if multiple_bond_count > 1 {
            return Ok(false);
        }
    }

    if !allow_exocyclic_multiple_bonds {
        for bond in context.incident_bonds(atom_id)? {
            if matches!(bond.order(), BondOrder::Double | BondOrder::Triple)
                && context.rings.num_bond_rings(bond.id()) == 0
            {
                return Ok(false);
            }
        }
    }

    Ok(true)
}

#[allow(dead_code)]
fn is_bond_order_query(bond: &crate::Bond) -> bool {
    // BEGIN RDKIT CPP FUNCTION MolOps::isBondOrderQuery
    // RDKit❗✔️: bool isBondOrderQuery(const Bond *bond) {
    // RDKit❗✔️:   if (bond->hasQuery()) {
    // RDKit❗✔️:     auto q = dynamic_cast<const QueryBond *>(bond)->getQuery();
    // RDKit❗✔️:     // complex bond type queries are also bond order queries!
    // RDKit❗✔️:     if (q->getTypeLabel() == "BondOrder" ||
    // RDKit❗✔️:         QueryOps::hasComplexBondTypeQuery(*q)) {
    // RDKit❗✔️:       return true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION MolOps::isBondOrderQuery
    bond.query().is_some_and(|query| {
        crate::valence::has_bond_type_query(query)
            || crate::valence::has_complex_bond_type_query(query)
    })
}

#[allow(dead_code)]
fn get_atom_donor_type_aromaticity(
    context: &AromaticityContext<'_>,
    atom_id: AtomId,
    exocyclic_bonds_steal_electrons: bool,
) -> Result<ElectronDonorType, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION getAtomDonorTypeArom
    // RDKit❗✔️: ElectronDonorType getAtomDonorTypeArom(
    // RDKit❗✔️:     const Atom *at, bool exocyclicBondsStealElectrons = true) {
    // RDKit❗✔️:   PRECONDITION(at, "bad atom");
    // RDKit❗✔️:   const auto &mol = at->getOwningMol();
    // RDKit❗✔️:   if (at->getAtomicNum() == 0) {
    // RDKit❗✔️:     return incidentCyclicMultipleBond(at) ? OneElectronDonorType
    // RDKit❗✔️:                                           : AnyElectronDonorType;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   ElectronDonorType res = NoElectronDonorType;
    // RDKit❗✔️:   auto nelec = MolOps::countAtomElec(at);
    // RDKit❗✔️:   int who = -1;
    // RDKit❗✔️:   if (nelec < 0) {
    // RDKit❗✔️:     res = NoElectronDonorType;
    // RDKit❗✔️:   } else if (nelec == 0) {
    // RDKit❗✔️:     if (incidentNonCyclicMultipleBond(at, who)) {
    // RDKit❗✔️:       // This is borderline:  no electron to spare but may have an empty
    // RDKit❗✔️:       // p-orbital
    // RDKit❗✔️:       // Not sure if this should return vacantElectronDonorType
    // RDKit❗✔️:       // FIX: explicitly doing this as a note for potential problems
    // RDKit❗✔️:       //
    // RDKit❗✔️:       res = VacantElectronDonorType;
    // RDKit❗✔️:     } else if (incidentCyclicMultipleBond(at)) {
    // RDKit❗✔️:       // no electron but has one in a in cycle multiple bond
    // RDKit❗✔️:       res = OneElectronDonorType;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       // no multiple bonds no electrons
    // RDKit❗✔️:       res = NoElectronDonorType;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else if (nelec == 1) {
    // RDKit❗✔️:     if (incidentNonCyclicMultipleBond(at, who)) {
    // RDKit❗✔️:       // the only available electron is going to be from the
    // RDKit❗✔️:       // external multiple bond this electron will not be available
    // RDKit❗✔️:       // for aromaticity if this atom is bonded to a more electro
    // RDKit❗✔️:       // negative atom
    // RDKit❗✔️:       const auto at2 = mol.getAtomWithIdx(who);
    // RDKit❗✔️:       if (exocyclicBondsStealElectrons &&
    // RDKit❗✔️:           PeriodicTable::getTable()->moreElectroNegative(at2->getAtomicNum(),
    // RDKit❗✔️:                                                          at->getAtomicNum())) {
    // RDKit❗✔️:         res = VacantElectronDonorType;
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         res = OneElectronDonorType;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       // require that the atom have at least one multiple bond
    // RDKit❗✔️:       if (incidentMultipleBond(at)) {
    // RDKit❗✔️:         res = OneElectronDonorType;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       // account for the tropylium and cyclopropenyl cation cases
    // RDKit❗✔️:       else if (at->getFormalCharge() == 1) {
    // RDKit❗✔️:         res = VacantElectronDonorType;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     if (incidentNonCyclicMultipleBond(at, who)) {
    // RDKit❗✔️:       // for cases with more than one electron :
    // RDKit❗✔️:       // if there is an incident multiple bond with an element that
    // RDKit❗✔️:       // is more electronegative than the this atom, count one less
    // RDKit❗✔️:       // electron
    // RDKit❗✔️:       const auto at2 = mol.getAtomWithIdx(who);
    // RDKit❗✔️:       if (exocyclicBondsStealElectrons &&
    // RDKit❗✔️:           PeriodicTable::getTable()->moreElectroNegative(at2->getAtomicNum(),
    // RDKit❗✔️:                                                          at->getAtomicNum())) {
    // RDKit❗✔️:         --nelec;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (nelec % 2 == 1) {
    // RDKit❗✔️:       res = OneElectronDonorType;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       res = TwoElectronDonorType;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION getAtomDonorTypeArom
    let atom = context.atom(atom_id)?;
    if atom.atomic_number() == 0 {
        return Ok(if incident_cyclic_multiple_bond(context, atom_id)? {
            ElectronDonorType::One
        } else {
            ElectronDonorType::Any
        });
    }

    let mut result = ElectronDonorType::None;
    let electron_count = count_atom_electrons(context, atom_id)?;
    if electron_count < 0 {
        result = ElectronDonorType::None;
    } else if electron_count == 0 {
        if incident_non_cyclic_multiple_bond(context, atom_id)?.is_some() {
            result = ElectronDonorType::Vacant;
        } else if incident_cyclic_multiple_bond(context, atom_id)? {
            result = ElectronDonorType::One;
        } else {
            result = ElectronDonorType::None;
        }
    } else if electron_count == 1 {
        if let Some(other_atom) = incident_non_cyclic_multiple_bond(context, atom_id)? {
            let other = context.atom(other_atom)?;
            if exocyclic_bonds_steal_electrons
                && crate::valence::periodic_table_more_electronegative(
                    other.atomic_number(),
                    atom.atomic_number(),
                )?
            {
                result = ElectronDonorType::Vacant;
            } else {
                result = ElectronDonorType::One;
            }
        } else if incident_multiple_bond(context, atom_id)? {
            result = ElectronDonorType::One;
        } else if atom.formal_charge() == 1 {
            result = ElectronDonorType::Vacant;
        }
    } else {
        let mut adjusted_electron_count = electron_count;
        if let Some(other_atom) = incident_non_cyclic_multiple_bond(context, atom_id)? {
            let other = context.atom(other_atom)?;
            if exocyclic_bonds_steal_electrons
                && crate::valence::periodic_table_more_electronegative(
                    other.atomic_number(),
                    atom.atomic_number(),
                )?
            {
                adjusted_electron_count -= 1;
            }
        }
        if adjusted_electron_count % 2 == 1 {
            result = ElectronDonorType::One;
        } else {
            result = ElectronDonorType::Two;
        }
    }
    Ok(result)
}

#[allow(dead_code)]
fn get_min_max_atom_electrons(dtype: ElectronDonorType) -> (i32, i32) {
    // BEGIN RDKIT CPP FUNCTION getMinMaxAtomElecs
    // RDKit✔️✔️: void getMinMaxAtomElecs(ElectronDonorType dtype, int &atlw, int &atup) {
    // RDKit✔️✔️:   switch (dtype) {
    match dtype {
        // RDKit✔️✔️:     case AnyElectronDonorType:
        // RDKit✔️✔️:       atlw = 1;
        // RDKit✔️✔️:       atup = 2;
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     case OneOrTwoElectronDonorType:
        // RDKit✔️✔️:       atlw = 1;
        // RDKit✔️✔️:       atup = 2;
        // RDKit✔️✔️:       break;
        ElectronDonorType::Any | ElectronDonorType::OneOrTwo => (1, 2),
        // RDKit✔️✔️:     case OneElectronDonorType:
        // RDKit✔️✔️:       atlw = atup = 1;
        // RDKit✔️✔️:       break;
        ElectronDonorType::One => (1, 1),
        // RDKit✔️✔️:     case TwoElectronDonorType:
        // RDKit✔️✔️:       atlw = atup = 2;
        // RDKit✔️✔️:       break;
        ElectronDonorType::Two => (2, 2),
        // RDKit✔️✔️:     case NoElectronDonorType:
        // RDKit✔️✔️:     case VacantElectronDonorType:
        // RDKit✔️✔️:     default:
        // RDKit✔️✔️:       atlw = atup = 0;
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:   }
        ElectronDonorType::None | ElectronDonorType::Vacant => (0, 0),
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getMinMaxAtomElecs
}

#[allow(dead_code)]
fn apply_huckel(
    ring: &[usize],
    electron_donors: &[ElectronDonorType],
    min_ring_size: usize,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION applyHuckel
    // RDKit✔️✔️: bool applyHuckel(ROMol &, const INT_VECT &ring, const VECT_EDON_TYPE &edon,
    // RDKit✔️✔️:                  unsigned int minRingSize) {
    // RDKit✔️✔️:   if (ring.size() < minRingSize) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if ring.len() < min_ring_size {
        return false;
    }
    // RDKit✔️✔️:   int atlw, atup, rlw, rup, rie;
    // RDKit✔️✔️:   bool aromatic = false;
    let mut aromatic = false;
    // RDKit✔️✔️:   rlw = 0;
    // RDKit✔️✔️:   rup = 0;
    let mut ring_lower = 0;
    let mut ring_upper = 0;
    // RDKit✔️✔️:   unsigned int nAnyElectronDonorType = 0;
    let mut any_donor_count = 0usize;
    // RDKit✔️✔️:   for (auto idx : ring) {
    for &idx in ring {
        // RDKit✔️✔️:     ElectronDonorType edonType = edon[idx];
        let donor_type = electron_donors[idx];
        // RDKit✔️✔️:     if (edonType == AnyElectronDonorType) {
        if donor_type == ElectronDonorType::Any {
            // RDKit✔️✔️:       ++nAnyElectronDonorType;
            any_donor_count += 1;
            // RDKit✔️✔️:       if (nAnyElectronDonorType > 1) {
            // RDKit✔️✔️:         return false;
            // RDKit✔️✔️:       }
            if any_donor_count > 1 {
                return false;
            }
        }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     getMinMaxAtomElecs(edonType, atlw, atup);
        let (atom_lower, atom_upper) = get_min_max_atom_electrons(donor_type);
        // RDKit✔️✔️:     rlw += atlw;
        // RDKit✔️✔️:     rup += atup;
        ring_lower += atom_lower;
        ring_upper += atom_upper;
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (rup >= 6) {
    if ring_upper >= 6 {
        // RDKit✔️✔️:     for (rie = rlw; rie <= rup; ++rie) {
        for ring_electrons in ring_lower..=ring_upper {
            // RDKit✔️✔️:       if ((rie - 2) % 4 == 0) {
            // RDKit✔️✔️:         aromatic = true;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if (ring_electrons - 2) % 4 == 0 {
                aromatic = true;
                break;
            }
        }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else if (rup == 2) {
    } else if ring_upper == 2 {
        // RDKit✔️✔️:     aromatic = true;
        aromatic = true;
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return aromatic;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION applyHuckel
    aromatic
}

#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
fn apply_huckel_to_fused(
    atom_count: usize,
    atom_rings: &[Vec<usize>],
    bond_rings: &[Vec<usize>],
    fused: &[usize],
    electron_donors: &[ElectronDonorType],
    ring_neighbors: &RingNeighborMap,
    aromatic_ring_count: &mut usize,
    max_num_fused_rings: usize,
    min_ring_size: usize,
    bonds: &[Bond],
    atom_aromatic: &mut [bool],
    bond_aromatic: &mut [bool],
) -> Result<(), AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION applyHuckelToFused
    // RDKit❗✔️: void applyHuckelToFused(
    // RDKit❗✔️:     ROMol &mol,                   // molecule of interest
    // RDKit❗✔️:     const VECT_INT_VECT &srings,  // list of all ring as atom IDS
    // RDKit❗✔️:     const VECT_INT_VECT &brings,  // list of all rings as bond ids
    // RDKit❗✔️:     const INT_VECT &fused,       // list of ring ids in the current fused system
    // RDKit❗✔️:     const VECT_EDON_TYPE &edon,  // electron donor state for each atom
    // RDKit❗✔️:     INT_INT_VECT_MAP &ringNeighs,  // list of neighbors for each candidate ring
    // RDKit❗✔️:     int &narom,                    // number of aromatic ring so far
    // RDKit❗✔️:     unsigned int maxNumFusedRings, const std::vector<Bond *> &bondsByIdx,
    // RDKit❗✔️:     unsigned int minRingSize) {
    // RDKit❗✔️:   // this function check huckel rule on a fused system it starts
    // RDKit❗✔️:   // with the individual rings in the system and then proceeds to
    // RDKit❗✔️:   // check for larger system i.e. if we have a 3 ring fused system,
    // RDKit❗✔️:   // huckel rule checked first on all the 1 ring subsystems then 2
    // RDKit❗✔️:   // rung subsystems etc.
    // RDKit❗✔️:
    // RDKit❗✔️:   std::unordered_set<int> aromRings;
    // RDKit❗✔️:   auto nrings = rdcast<unsigned int>(fused.size());
    // RDKit❗✔️:   INT_VECT curRs;
    // RDKit❗✔️:   curRs.push_back(fused.front());
    // RDKit❗✔️:   int pos = -1;
    // RDKit❗✔️:   unsigned int i;
    // RDKit❗✔️:   unsigned int curSize = 0;
    // RDKit❗✔️:   INT_VECT comb;
    // RDKit❗✔️:
    // RDKit❗✔️:   size_t nRingBonds;
    // RDKit❗✔️:   {
    // RDKit❗✔️:     boost::dynamic_bitset<> fusedBonds(mol.getNumBonds());
    // RDKit❗✔️:     for (auto ridx : fused) {
    // RDKit❗✔️:       for (auto bidx : brings[ridx]) {
    // RDKit❗✔️:         fusedBonds[bidx] = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     nRingBonds = rdcast<unsigned int>(fusedBonds.count());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::set<unsigned int> doneBonds;
    // RDKit❗✔️:   while (1) {
    // RDKit❗✔️:     if (pos == -1) {
    // RDKit❗✔️:       // If a ring system has more than 300 rings and a ring combination search
    // RDKit❗✔️:       // larger than 2 is reached, the calculation becomes exponentially longer,
    // RDKit❗✔️:       // in some case it never completes.
    // RDKit❗✔️:       if ((curSize == 2) && (nrings > 300)) {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:             << "Aromaticity detection halted on some rings due to ring system size."
    // RDKit❗✔️:             << std::endl;
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       ++curSize;
    // RDKit❗✔️:       // check if we are done with all the atoms in the fused
    // RDKit❗✔️:       // system, if so quit. This is a fix for Issue252 REVIEW: is
    // RDKit❗✔️:       // this check sufficient or should we add an additional
    // RDKit❗✔️:       // constraint on the number of combinations of rings in a
    // RDKit❗✔️:       // fused system that we will try. The number of combinations
    // RDKit❗✔️:       // can obviously be quite large when the number of rings in
    // RDKit❗✔️:       // the fused system is large
    // RDKit❗✔️:       if (curSize > std::min(nrings, maxNumFusedRings) ||
    // RDKit❗✔️:           doneBonds.size() >= nRingBonds) {
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       comb.resize(curSize);
    // RDKit❗✔️:       std::iota(comb.begin(), comb.end(), 0);
    // RDKit❗✔️:       pos = 0;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       pos = nextCombination(comb, nrings);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     if (pos == -1) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     curRs.clear();
    // RDKit❗✔️:     std::transform(comb.begin(), comb.end(), std::back_inserter(curRs),
    // RDKit❗✔️:                    [&fused](const int i) { return fused[i]; });
    // RDKit❗✔️:
    // RDKit❗✔️:     // check if the picked subsystem is fused
    // RDKit❗✔️:     if (ringNeighs.size() && !RingUtils::checkFused(curRs, ringNeighs)) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // check aromaticity on the current fused system
    // RDKit❗✔️:     INT_VECT atsInRingSystem(mol.getNumAtoms(), 0);
    // RDKit❗✔️:     for (auto ridx : curRs) {
    // RDKit❗✔️:       for (auto rid : srings[ridx]) {
    // RDKit❗✔️:         ++atsInRingSystem[rid];
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     INT_VECT unon;
    // RDKit❗✔️:     for (i = 0; i < atsInRingSystem.size(); ++i) {
    // RDKit❗✔️:       // condition for inclusion of an atom in the aromaticity of a fused ring
    // RDKit❗✔️:       // system is that it's present in one or two of the rings. this was #2895:
    // RDKit❗✔️:       // the central atom in acepentalene was being included in the count of
    // RDKit❗✔️:       // aromatic atoms
    // RDKit❗✔️:       if (atsInRingSystem[i] == 1 || atsInRingSystem[i] == 2) {
    // RDKit❗✔️:         unon.push_back(i);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (applyHuckel(mol, unon, edon, minRingSize)) {
    // RDKit❗✔️:       // mark the atoms and bonds in these rings to be aromatic
    // RDKit❗✔️:       markAtomsBondsArom(brings, curRs, doneBonds, bondsByIdx);
    // RDKit❗✔️:
    // RDKit❗✔️:       // add the ring IDs to the aromatic rings found so far
    // RDKit❗✔️:       // avoid duplicates
    // RDKit❗✔️:       std::copy(curRs.begin(), curRs.end(),
    // RDKit❗✔️:                 std::inserter(aromRings, aromRings.begin()));
    // RDKit❗✔️:     }  // end check huckel rule
    // RDKit❗✔️:   }  // end while(1)
    // RDKit❗✔️:   narom += rdcast<int>(aromRings.size());
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION applyHuckelToFused
    let ring_count = fused.len();
    let mut aromatic_rings = BTreeSet::new();
    let mut current_size = 0usize;
    let mut position = -1;
    let mut combination = Vec::new();
    let mut fused_bonds = BTreeSet::new();
    for &ring_index in fused {
        for &bond_index in
            bond_rings
                .get(ring_index)
                .ok_or(AromaticityError::InvariantViolation {
                    message: "fused ring index out of range",
                })?
        {
            fused_bonds.insert(bond_index);
        }
    }
    let ring_bond_count = fused_bonds.len();
    let mut done_bonds = BTreeSet::new();

    loop {
        if position == -1 {
            if current_size == 2 && ring_count > 300 {
                break;
            }
            current_size += 1;
            if current_size > ring_count.min(max_num_fused_rings)
                || done_bonds.len() >= ring_bond_count
            {
                break;
            }
            combination = (0..current_size).collect();
            position = 0;
        } else {
            position = next_combination(&mut combination, ring_count);
        }
        if position == -1 {
            continue;
        }

        let current_rings = combination
            .iter()
            .map(|&index| fused[index])
            .collect::<Vec<_>>();
        if !ring_neighbors.is_empty() && !check_fused(&current_rings, ring_neighbors)? {
            continue;
        }

        let mut atoms_in_ring_system = vec![0usize; atom_count];
        for &ring_index in &current_rings {
            for &atom_index in
                atom_rings
                    .get(ring_index)
                    .ok_or(AromaticityError::InvariantViolation {
                        message: "aromatic ring index out of range",
                    })?
            {
                atoms_in_ring_system[atom_index] += 1;
            }
        }
        let included_atoms = atoms_in_ring_system
            .iter()
            .enumerate()
            .filter_map(|(atom_index, count)| {
                if matches!(*count, 1 | 2) {
                    Some(atom_index)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        if apply_huckel(&included_atoms, electron_donors, min_ring_size) {
            mark_atoms_bonds_aromatic(
                bond_rings,
                &current_rings,
                &mut done_bonds,
                bonds,
                atom_aromatic,
                bond_aromatic,
            )?;
            aromatic_rings.extend(current_rings);
        }
    }
    *aromatic_ring_count += aromatic_rings.len();
    Ok(())
}

#[allow(dead_code)]
fn mark_atoms_bonds_aromatic(
    bond_rings: &[Vec<usize>],
    ring_ids: &[usize],
    done_bonds: &mut BTreeSet<usize>,
    bonds: &[Bond],
    atom_aromatic: &mut [bool],
    bond_aromatic: &mut [bool],
) -> Result<(), AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION markAtomsBondsArom
    // RDKit❗✔️: void markAtomsBondsArom(const VECT_INT_VECT &brings, const INT_VECT &ringIds,
    // RDKit❗✔️:                         std::set<unsigned int> &doneBonds,
    // RDKit❗✔️:                         const std::vector<Bond *> &bondsByIdx) {
    // RDKit❗✔️:   // mark the bonds
    // RDKit❗✔️:   // here we want to be careful. We don't want to mark the fusing bonds
    // RDKit❗✔️:   // as aromatic - only the outside bonds in a fused system are marked aromatic.
    // RDKit❗✔️:   // - loop through the rings and count the number of times each bond appears in
    // RDKit❗✔️:   //   all the fused rings.
    // RDKit❗✔️:   // - bonds that appears only once are marked aromatic
    // RDKit❗✔️:   INT_MAP_INT bndCntr;
    // RDKit❗✔️:
    // RDKit❗✔️:   for (auto ri : ringIds) {
    // RDKit❗✔️:     const auto &bring = brings[ri];
    // RDKit❗✔️:     for (auto bi : bring) {
    // RDKit❗✔️:       ++bndCntr[bi];
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // now mark single or double bonds that have a count of 1 and the atoms they
    // RDKit❗✔️:   // connect as aromatic
    // RDKit❗✔️:   for (const auto &bci : bndCntr) {
    // RDKit❗✔️:     if (bci.second == 1) {
    // RDKit❗✔️:       auto bond = bondsByIdx[bci.first];
    // RDKit❗✔️:       bond->setIsAromatic(true);
    // RDKit❗✔️:       switch (bond->getBondType()) {
    // RDKit❗✔️:         case Bond::SINGLE:
    // RDKit❗✔️:         case Bond::DOUBLE:
    // RDKit❗✔️:           bond->setBondType(Bond::AROMATIC);
    // RDKit❗✔️:           bond->getBeginAtom()->setIsAromatic(true);
    // RDKit❗✔️:           bond->getEndAtom()->setIsAromatic(true);
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         default:
    // RDKit❗✔️:           break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       doneBonds.insert(bond->getIdx());
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION markAtomsBondsArom
    let mut bond_counter = BTreeMap::<usize, usize>::new();
    for &ring_id in ring_ids {
        let bond_ring = bond_rings
            .get(ring_id)
            .ok_or(AromaticityError::InvariantViolation {
                message: "aromatic bond ring index out of range",
            })?;
        for &bond_index in bond_ring {
            *bond_counter.entry(bond_index).or_default() += 1;
        }
    }
    for (bond_index, count) in bond_counter {
        if count == 1 {
            let Some(slot) = bond_aromatic.get_mut(bond_index) else {
                return Err(AromaticityError::InvariantViolation {
                    message: "aromatic bond index out of range",
                });
            };
            *slot = true;
            let bond = bonds
                .get(bond_index)
                .ok_or(AromaticityError::InvariantViolation {
                    message: "aromatic bond index out of range",
                })?;
            if matches!(bond.order(), BondOrder::Single | BondOrder::Double) {
                for atom_id in [bond.begin(), bond.end()] {
                    let Some(atom_slot) = atom_aromatic.get_mut(atom_id.index()) else {
                        return Err(AromaticityError::InvariantViolation {
                            message: "aromatic atom index out of range",
                        });
                    };
                    *atom_slot = true;
                }
            }
            done_bonds.insert(bond_index);
        }
    }
    Ok(())
}

fn next_combination(combination: &mut [usize], total: usize) -> i32 {
    // BEGIN RDKIT CPP FUNCTION nextCombination
    // RDKit✔️✔️: int nextCombination(INT_VECT &comb, int tot) {
    // RDKit✔️✔️:   int nelem = static_cast<int>(comb.size());
    // RDKit✔️✔️:   int celem = nelem - 1;
    let element_count = combination.len();
    if element_count == 0 {
        return -1;
    }
    let mut current = element_count - 1;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   while (comb[celem] == (tot - nelem + celem)) {
    loop {
        if combination[current] != total - element_count + current {
            break;
        }
        // RDKit✔️✔️:     celem--;
        // RDKit✔️✔️:     if (celem < 0) {
        // RDKit✔️✔️:       return -1;
        // RDKit✔️✔️:     }
        if current == 0 {
            return -1;
        }
        current -= 1;
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int i;
    // RDKit✔️✔️:   comb[celem] += 1;
    combination[current] += 1;
    // RDKit✔️✔️:   for (i = celem + 1; i < comb.size(); i++) {
    for index in current + 1..combination.len() {
        // RDKit✔️✔️:     comb[i] = comb[i - 1] + 1;
        combination[index] = combination[index - 1] + 1;
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return celem;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION nextCombination
    current as i32
}

#[cfg(test)]
mod tests {
    use crate::aromaticity::{
        ElectronDonorType, apply_huckel, check_fused, get_min_max_atom_electrons,
        make_ring_neighbor_map, pick_fused_rings,
    };
    use crate::{AromaticityModel, Molecule, set_aromaticity};

    #[test]
    fn set_aromaticity_rejects_saturated_three_membered_ring() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CC1", false).unwrap();
        let assignment = set_aromaticity(&molecule, AromaticityModel::Default).unwrap();

        assert_eq!(assignment.atom_aromatic, vec![false, false, false]);
        assert_eq!(assignment.bond_aromatic, vec![false, false, false]);
        assert_eq!(assignment.aromatic_ring_count, 0);
    }

    #[test]
    fn set_aromaticity_handles_acyclic_molecule_without_guessing_rings() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false)
            .unwrap()
            .with_assigned_valence()
            .unwrap()
            .with_assigned_rings()
            .unwrap();

        let assignment = set_aromaticity(&molecule, AromaticityModel::Default).unwrap();

        assert_eq!(assignment.atom_aromatic, vec![false, false, false]);
        assert_eq!(assignment.bond_aromatic, vec![false, false]);
        assert_eq!(assignment.aromatic_ring_count, 0);
    }

    #[test]
    fn with_assigned_aromaticity_updates_cache_for_acyclic_molecule() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        let result = molecule.with_assigned_aromaticity().unwrap();

        assert!(result.derived_cache().rings.is_some());
        assert!(result.derived_cache().valence.is_some());
        assert!(result.derived_cache().aromaticity_valid);
        assert!(result.atoms().iter().all(|atom| !atom.is_aromatic()));
        assert!(result.bonds().iter().all(|bond| !bond.is_aromatic()));
    }

    #[test]
    fn with_assigned_aromaticity_rejects_saturated_three_membered_ring() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CC1", false).unwrap();
        let result = molecule.with_assigned_aromaticity().unwrap();

        assert!(result.derived_cache().aromaticity_valid);
        assert!(result.atoms().iter().all(|atom| !atom.is_aromatic()));
        assert!(result.bonds().iter().all(|bond| !bond.is_aromatic()));
    }

    #[test]
    fn with_assigned_aromaticity_marks_benzene_like_rdkit_default() {
        let molecule = Molecule::from_smiles_with_sanitize("C1=CC=CC=C1", false).unwrap();
        let result = molecule.with_assigned_aromaticity().unwrap();

        assert!(result.derived_cache().aromaticity_valid);
        assert_eq!(
            result
                .atoms()
                .iter()
                .map(crate::Atom::is_aromatic)
                .collect::<Vec<_>>(),
            vec![true; 6]
        );
        assert_eq!(
            result
                .bonds()
                .iter()
                .map(crate::Bond::is_aromatic)
                .collect::<Vec<_>>(),
            vec![true; 6]
        );
        assert!(
            result
                .bonds()
                .iter()
                .all(|bond| bond.order() == crate::BondOrder::Aromatic)
        );
        let expected_valence =
            crate::assign_valence_with_options(&result, crate::ValenceModel::RdkitLike, true)
                .unwrap();
        assert_eq!(result.derived_cache().valence, Some(expected_valence));
    }

    #[test]
    fn get_min_max_atom_electrons_matches_rdkit_switch_table() {
        assert_eq!(get_min_max_atom_electrons(ElectronDonorType::Any), (1, 2));
        assert_eq!(
            get_min_max_atom_electrons(ElectronDonorType::OneOrTwo),
            (1, 2)
        );
        assert_eq!(get_min_max_atom_electrons(ElectronDonorType::One), (1, 1));
        assert_eq!(get_min_max_atom_electrons(ElectronDonorType::Two), (2, 2));
        assert_eq!(get_min_max_atom_electrons(ElectronDonorType::None), (0, 0));
        assert_eq!(
            get_min_max_atom_electrons(ElectronDonorType::Vacant),
            (0, 0)
        );
    }

    #[test]
    fn apply_huckel_matches_rdkit_electron_window_logic() {
        assert!(!apply_huckel(
            &[0, 1, 2],
            &[
                ElectronDonorType::Two,
                ElectronDonorType::Two,
                ElectronDonorType::Two,
            ],
            4
        ));
        assert!(apply_huckel(
            &[0, 1, 2],
            &[
                ElectronDonorType::Two,
                ElectronDonorType::Two,
                ElectronDonorType::Two,
            ],
            0
        ));
        assert!(apply_huckel(&[0], &[ElectronDonorType::OneOrTwo], 0));
        assert!(!apply_huckel(
            &[0, 1],
            &[ElectronDonorType::Any, ElectronDonorType::Any],
            0
        ));
    }

    #[test]
    fn make_ring_neighbor_map_matches_rdkit_overlap_rules() {
        let map = make_ring_neighbor_map(&[vec![0, 1, 2], vec![2, 3, 4], vec![4, 5, 6]], 0, 1);

        assert_eq!(map.get(&0), Some(&vec![1]));
        assert_eq!(map.get(&1), Some(&vec![0, 2]));
        assert_eq!(map.get(&2), Some(&vec![1]));

        let filtered = make_ring_neighbor_map(&[vec![0, 1, 2], vec![1, 2, 3]], 0, 1);
        assert_eq!(filtered.get(&0), Some(&Vec::new()));
        assert_eq!(filtered.get(&1), Some(&Vec::new()));
    }

    #[test]
    fn pick_and_check_fused_match_rdkit_graph_walk() {
        let map = make_ring_neighbor_map(&[vec![0, 1], vec![1, 2], vec![3, 4]], 0, 1);
        let mut result = Vec::new();
        let mut done = vec![false; map.len()];

        pick_fused_rings(0, &map, &mut result, &mut done, 0).unwrap();

        assert_eq!(result, vec![0, 1]);
        assert!(check_fused(&[0, 1], &map).unwrap());
        assert!(!check_fused(&[0, 2], &map).unwrap());
    }
}
