//! RDKit-aligned early sanitize cleanup over detached topology values.

use std::cmp::Ordering;

use cosmolkit_model::{AtomId, TopologyBlock, TopologyValidationError};
use cosmolkit_types::BondOrder;

use crate::{
    CanonicalRankError, ValenceError, calculate_explicit_valence_for_topology,
    calculate_implicit_valence_for_topology, rdkit_valence_list,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CleanupParams {
    pub charge_normalization: bool,
    pub organometallics: bool,
}

impl Default for CleanupParams {
    fn default() -> Self {
        Self {
            charge_normalization: true,
            organometallics: true,
        }
    }
}

#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum CleanupError {
    #[error("invalid topology: {0}")]
    InvalidTopology(#[from] TopologyValidationError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
    #[error(transparent)]
    CanonicalRank(#[from] CanonicalRankError),
}

/// Run the independently selectable source cleanup stages in sanitize order.
pub fn cleanup(
    topology: &TopologyBlock,
    params: &CleanupParams,
) -> Result<TopologyBlock, CleanupError> {
    topology.validate()?;
    let mut result = topology.clone();
    if params.charge_normalization {
        cleanup_charges(&mut result)?;
    }
    if params.organometallics {
        cleanup_organometallics(&mut result)?;
    }
    result.validate()?;
    Ok(result)
}

fn cleanup_nitrogens(topology: &mut TopologyBlock) -> Result<(), CleanupError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: nitrogensCleanup
    // RDKit✔️✔️: void nitrogensCleanup(RWMol &mol) {
    // RDKit✔️✔️:   // conversions here:
    // RDKit✔️✔️:   // - neutral 5 coordinate Ns with double bonds to Os to the
    // RDKit✔️✔️:   //   zwitterionic form.  e.g.:
    // RDKit✔️✔️:   //   CN(=O)=O -> C[N+](=O)[O-]
    // RDKit✔️✔️:   //   and:
    // RDKit✔️✔️:   //   C1=CC=CN(=O)=C1 -> C1=CC=C[N+]([O-])=C1
    // RDKit✔️✔️:   // - neutral 5 coordinate Ns with triple bonds to Ns to the
    // RDKit✔️✔️:   //   zwitterionic form.  e.g.:
    // RDKit✔️✔️:   //   C-N=N#N -> C-N=[N+]=[N-]
    // RDKit✔️✔️:
    // RDKit✔️✔️:   boost::dynamic_bitset<> nitrogensToConsider(mol.getNumAtoms());
    let mut nitrogens_to_consider = vec![false; topology.atoms.len()];
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    for atom_index in 0..topology.atoms.len() {
        let atom_id = AtomId::new(atom_index);
        // RDKit✔️✔️:     if (atom->getAtomicNum() != 7) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if topology.atoms[atom_index].atomic_number() != 7 {
            continue;
        }
        // RDKit✔️✔️:     // we only want to do neutrals so that things like this don't get
        // RDKit✔️✔️:     // munged:
        // RDKit✔️✔️:     //  O=[n+]1occcc1
        // RDKit✔️✔️:     // this was sf.net issue 1811276
        // RDKit✔️✔️:     if (atom->getFormalCharge()) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if topology.atoms[atom_index].formal_charge() != 0 {
            continue;
        }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     // NOTE that we are calling calcExplicitValence() here, we do
        // RDKit✔️✔️:     // this because we cannot be sure that it has already been
        // RDKit✔️✔️:     // called on the atom (cleanUp() gets called pretty early in
        // RDKit✔️✔️:     // the sanitization process):
        // RDKit✔️✔️:     if (atom->calcExplicitValence(false) != 5) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if calculate_explicit_valence_for_topology(topology, atom_id, false, false)? != 5 {
            continue;
        }
        // RDKit✔️✔️:     nitrogensToConsider.set(atom->getIdx());
        nitrogens_to_consider[atom_index] = true;
        // RDKit✔️✔️:     // we need to play this little aromaticity game because the
        // RDKit✔️✔️:     // explicit valence code modifies its results for aromatic
        // RDKit✔️✔️:     // atoms.
        // RDKit✔️✔️:     auto aromHolder = atom->getIsAromatic();
        // RDKit✔️✔️:     atom->setIsAromatic(0);
        // RDKit✔️✔️:     unsigned int aid = atom->getIdx();
        // RDKit✔️✔️:     bool updateNeeded = false;
        let neighbors = topology.adjacency.neighbors_of(atom_index).to_vec();
        let mut update = None;
        // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
        for neighbor in neighbors {
            let neighbor_atom = &topology.atoms[neighbor.atom_index];
            let bond = &topology.bonds[neighbor.bond.index()];
            // RDKit✔️✔️:       if ((nbr->getAtomicNum() == 8) && (nbr->getFormalCharge() == 0) &&
            // RDKit✔️✔️:           (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
            // RDKit✔️✔️:            Bond::DOUBLE)) {
            if neighbor_atom.atomic_number() == 8
                && neighbor_atom.formal_charge() == 0
                && bond.order() == BondOrder::Double
            {
                // RDKit✔️✔️:         // here's the double bonded oxygen
                // RDKit✔️✔️:         auto b = mol.getBondBetweenAtoms(aid, nbr->getIdx());
                // RDKit✔️✔️:         b->setBondType(Bond::SINGLE);
                // RDKit✔️✔️:         atom->setFormalCharge(1);
                // RDKit✔️✔️:         nbr->setFormalCharge(-1);
                // RDKit✔️✔️:         updateNeeded = true;
                // RDKit✔️✔️:         break;
                update = Some((neighbor.bond, neighbor.atom_index));
                break;
                // RDKit✔️✔️:       }
            }
        }
        if let Some((bond_id, oxygen_index)) = update {
            topology.bonds[bond_id.index()].set_order(BondOrder::Single);
            topology.atoms[atom_index].set_formal_charge(1);
            topology.atoms[oxygen_index].set_formal_charge(-1);
        }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     // force a recalculation of the explicit valence if we changed anything
        // RDKit✔️✔️:     atom->setIsAromatic(aromHolder);
        // RDKit✔️✔️:     if (updateNeeded) {
        // RDKit✔️✔️:       atom->calcExplicitValence(false);
        // RDKit✔️✔️:     }
        // Detached values have no property cache to recalculate and preserve
        // the source aromatic flag throughout.
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // now repeat for the weird N#N case:
    // RDKit✔️✔️:   for (auto aid = nitrogensToConsider.find_first();
    // RDKit✔️✔️:        aid != boost::dynamic_bitset<>::npos;
    // RDKit✔️✔️:        aid = nitrogensToConsider.find_next(aid)) {
    for (atom_index, selected) in nitrogens_to_consider.into_iter().enumerate() {
        if !selected {
            continue;
        }
        // RDKit✔️✔️:     Atom *atom = mol.getAtomWithIdx(aid);
        // RDKit✔️✔️:     auto aromHolder = atom->getIsAromatic();
        // RDKit✔️✔️:     atom->setIsAromatic(0);
        // RDKit✔️✔️:     bool updateNeeded = false;
        let neighbors = topology.adjacency.neighbors_of(atom_index).to_vec();
        let mut update = None;
        // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
        for neighbor in neighbors {
            let neighbor_atom = &topology.atoms[neighbor.atom_index];
            let bond = &topology.bonds[neighbor.bond.index()];
            // RDKit✔️✔️:       if ((nbr->getAtomicNum() == 7) && (nbr->getFormalCharge() == 0) &&
            // RDKit✔️✔️:           (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
            // RDKit✔️✔️:            Bond::TRIPLE)) {
            if neighbor_atom.atomic_number() == 7
                && neighbor_atom.formal_charge() == 0
                && bond.order() == BondOrder::Triple
            {
                // RDKit✔️✔️:         // here's the triple bonded nitrogen
                // RDKit✔️✔️:         auto b = mol.getBondBetweenAtoms(aid, nbr->getIdx());
                // RDKit✔️✔️:         b->setBondType(Bond::DOUBLE);
                // RDKit✔️✔️:         atom->setFormalCharge(1);
                // RDKit✔️✔️:         nbr->setFormalCharge(-1);
                // RDKit✔️✔️:         updateNeeded = true;
                // RDKit✔️✔️:         break;
                update = Some((neighbor.bond, neighbor.atom_index));
                break;
                // RDKit✔️✔️:       }
            }
        }
        if let Some((bond_id, nitrogen_index)) = update {
            topology.bonds[bond_id.index()].set_order(BondOrder::Double);
            topology.atoms[atom_index].set_formal_charge(1);
            topology.atoms[nitrogen_index].set_formal_charge(-1);
        }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     // force a recalculation of the explicit valence here
        // RDKit✔️✔️:     atom->setIsAromatic(aromHolder);
        // RDKit✔️✔️:     if (updateNeeded) {
        // RDKit✔️✔️:       atom->calcExplicitValence(false);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: nitrogensCleanup
    Ok(())
}

fn cleanup_phosphorus(topology: &mut TopologyBlock, atom_index: usize) -> Result<(), CleanupError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: phosphorusCleanup
    // RDKit✔️✔️: void phosphorusCleanup(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   // conversions here:
    // RDKit✔️✔️:   // - neutral 5 coordinate Ps with one double bonds to an Os
    // RDKit✔️✔️:   //   and one to a C or N to the zwitterionic form.  e.g.:
    // RDKit✔️✔️:   //   C=P(=O)X -> C=[P+]([O-])X
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // we only want to do neutrals
    // RDKit✔️✔️:   if (atom->getFormalCharge()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    if topology.atoms[atom_index].formal_charge() != 0 {
        return Ok(());
    }
    let atom_id = AtomId::new(atom_index);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // NOTE that we are calling calcExplicitValence() here, we do
    // RDKit✔️✔️:   // this because we cannot be sure that it has already been
    // RDKit✔️✔️:   // called on the atom (cleanUp() gets called pretty early in
    // RDKit✔️✔️:   // the sanitization process):
    // RDKit✔️✔️:   if (atom->calcExplicitValence(false) == 5 && atom->getDegree() == 3) {
    if calculate_explicit_valence_for_topology(topology, atom_id, false, false)? == 5
        && topology.adjacency.neighbors_of(atom_index).len() == 3
    {
        // RDKit✔️✔️:     unsigned int aid = atom->getIdx();
        // RDKit✔️✔️:     Bond *dbl_to_O = nullptr;
        // RDKit✔️✔️:     Atom *O_atom = nullptr;
        // RDKit✔️✔️:     bool hasDoubleToCorN = false;
        let mut double_to_oxygen = None;
        let mut has_double_to_carbon_or_nitrogen = false;
        // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
        for neighbor in topology.adjacency.neighbors_of(atom_index) {
            let neighbor_atom = &topology.atoms[neighbor.atom_index];
            let bond = &topology.bonds[neighbor.bond.index()];
            // RDKit✔️✔️:       if ((nbr->getAtomicNum() == 8) && (nbr->getFormalCharge() == 0) &&
            // RDKit✔️✔️:           (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
            // RDKit✔️✔️:            Bond::DOUBLE)) {
            // RDKit✔️✔️:         // here's the double bonded oxygen
            // RDKit✔️✔️:         dbl_to_O = mol.getBondBetweenAtoms(aid, nbr->getIdx());
            // RDKit✔️✔️:         O_atom = nbr;
            if neighbor_atom.atomic_number() == 8
                && neighbor_atom.formal_charge() == 0
                && bond.order() == BondOrder::Double
            {
                double_to_oxygen = Some((neighbor.bond, neighbor.atom_index));
            // RDKit✔️✔️:       } else if ((nbr->getAtomicNum() == 6 || nbr->getAtomicNum() == 7) &&
            // RDKit✔️✔️:                  (nbr->getDegree() >= 2) &&
            // RDKit✔️✔️:                  (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
            // RDKit✔️✔️:                   Bond::DOUBLE)) {
            // RDKit✔️✔️:         hasDoubleToCorN = true;
            } else if matches!(neighbor_atom.atomic_number(), 6 | 7)
                && topology.adjacency.neighbors_of(neighbor.atom_index).len() >= 2
                && bond.order() == BondOrder::Double
            {
                has_double_to_carbon_or_nitrogen = true;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }  // end of loop over the first neigh
        }
        // RDKit✔️✔️:     if (hasDoubleToCorN && dbl_to_O != nullptr) {
        if has_double_to_carbon_or_nitrogen && let Some((bond_id, oxygen_index)) = double_to_oxygen
        {
            // RDKit✔️✔️:       TEST_ASSERT(O_atom != nullptr);
            // RDKit✔️✔️:       O_atom->setFormalCharge(-1);
            // RDKit✔️✔️:       dbl_to_O->setBondType(Bond::SINGLE);
            // RDKit✔️✔️:       atom->setFormalCharge(1);
            topology.atoms[oxygen_index].set_formal_charge(-1);
            topology.bonds[bond_id.index()].set_order(BondOrder::Single);
            topology.atoms[atom_index].set_formal_charge(1);
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   // force a recalculation of the explicit valence here
    // RDKit✔️✔️:   atom->calcExplicitValence(false);
    // Detached values do not persist a property cache.
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: phosphorusCleanup
    Ok(())
}

fn cleanup_halogen(topology: &mut TopologyBlock, atom_index: usize) -> Result<(), CleanupError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: halogenCleanup
    // RDKit✔️✔️: void halogenCleanup(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:   // Conversions done:
    // RDKit✔️✔️:   //    X(=O)(=O)(=O)O -> [X+3]([O-])([O-])([O-])O
    // RDKit✔️✔️:   //    X(=O)(=O)O -> [X+2]([O-])([O-])O
    // RDKit✔️✔️:   //    X(=O)O -> [X+]([O-])O
    // RDKit✔️✔️:   int ev = atom->calcExplicitValence(false);
    let atom_id = AtomId::new(atom_index);
    let explicit_valence =
        calculate_explicit_valence_for_topology(topology, atom_id, false, false)?;
    // RDKit✔️✔️:   if (atom->getFormalCharge() == 0 && (ev == 7 || ev == 5 || ev == 3)) {
    if topology.atoms[atom_index].formal_charge() == 0 && matches!(explicit_valence, 7 | 5 | 3) {
        // RDKit✔️✔️:     bool neighborsAllO = true;
        // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
        // RDKit✔️✔️:       if (nbr->getAtomicNum() != 8) {
        // RDKit✔️✔️:         neighborsAllO = false;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let neighbors = topology.adjacency.neighbors_of(atom_index).to_vec();
        let neighbors_all_oxygen = neighbors
            .iter()
            .all(|neighbor| topology.atoms[neighbor.atom_index].atomic_number() == 8);
        // RDKit✔️✔️:     if (neighborsAllO) {
        if neighbors_all_oxygen {
            // RDKit✔️✔️:       int formalCharge = 0;
            let mut formal_charge = 0i8;
            // RDKit✔️✔️:       for (auto bond : mol.atomBonds(atom)) {
            for neighbor in neighbors {
                // RDKit✔️✔️:         if (bond->getBondType() == Bond::DOUBLE) {
                if topology.bonds[neighbor.bond.index()].order() == BondOrder::Double {
                    // RDKit✔️✔️:           bond->setBondType(Bond::SINGLE);
                    // RDKit✔️✔️:           auto otherAtom = bond->getOtherAtom(atom);
                    // RDKit✔️✔️:           formalCharge++;
                    // RDKit✔️✔️:           otherAtom->setFormalCharge(-1);
                    // RDKit✔️✔️:           otherAtom->calcExplicitValence(false);
                    topology.bonds[neighbor.bond.index()].set_order(BondOrder::Single);
                    formal_charge += 1;
                    topology.atoms[neighbor.atom_index].set_formal_charge(-1);
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:       atom->setFormalCharge(formalCharge);
            // RDKit✔️✔️:       atom->calcExplicitValence(false);
            topology.atoms[atom_index].set_formal_charge(formal_charge);
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: halogenCleanup
    Ok(())
}

fn cleanup_charges(topology: &mut TopologyBlock) -> Result<(), CleanupError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: cleanUp
    // RDKit✔️✔️: void cleanUp(RWMol &mol) {
    // RDKit✔️✔️:   nitrogensCleanup(mol);
    cleanup_nitrogens(topology)?;
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    for atom_index in 0..topology.atoms.len() {
        // RDKit✔️✔️:     switch (atom->getAtomicNum()) {
        match topology.atoms[atom_index].atomic_number() {
            // RDKit✔️✔️:       case 15:
            // RDKit✔️✔️:         phosphorusCleanup(mol, atom);
            // RDKit✔️✔️:         break;
            15 => cleanup_phosphorus(topology, atom_index)?,
            // RDKit✔️✔️:       case 17:
            // RDKit✔️✔️:       case 35:
            // RDKit✔️✔️:       case 53:
            // RDKit✔️✔️:         halogenCleanup(mol, atom);
            // RDKit✔️✔️:         break;
            17 | 35 | 53 => cleanup_halogen(topology, atom_index)?,
            _ => {} // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: cleanUp
    Ok(())
}

fn is_metal(atomic_number: u8) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.cpp :: makeMAtomQuery/isMetal
    // RDKit✔️✔️: ATOM_OR_QUERY *makeMAtomQuery() {
    // RDKit✔️✔️:   // using the definition from Marvin Sketch, which produces the following
    // RDKit✔️✔️:   // SMARTS:
    // RDKit✔️✔️:   // !#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
    // RDKit✔️✔️:   // We expanded this with !#0 as part of #6106
    // RDKit✔️✔️:   // it's easier to define what isn't a metal than what is. :-)
    // RDKit✔️✔️:   ATOM_OR_QUERY *res = makeMHAtomQuery();
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(1)));
    // RDKit✔️✔️:   res->setTypeLabel("M");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: ATOM_OR_QUERY *makeMHAtomQuery() {
    // RDKit✔️✔️:   // using the definition from Marvin Sketch, which produces the following
    // RDKit✔️✔️:   // SMARTS:
    // RDKit✔️✔️:   // !#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
    // RDKit✔️✔️:   // We expanded this with !#0 as part of #6106
    // RDKit✔️✔️:   // it's easier to define what isn't a metal than what is. :-)
    // RDKit✔️✔️:   auto *res = new ATOM_OR_QUERY;
    // RDKit✔️✔️:   res->setDescription("AtomOr");
    // RDKit✔️✔️:   res->setNegation(true);
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(0)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(2)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(5)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(6)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(7)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(8)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(9)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(10)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(14)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(15)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(16)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(17)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(18)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(33)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(34)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(35)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(36)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(52)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(53)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(54)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(85)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(86)));
    // RDKit✔️✔️:   res->setTypeLabel("MH");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: bool isMetal(const Atom &atom) {
    // RDKit✔️✔️:   static const std::unique_ptr<ATOM_OR_QUERY> q(makeMAtomQuery());
    // RDKit✔️✔️:   return q->Match(&atom);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.cpp :: makeMAtomQuery/isMetal
    !matches!(
        atomic_number,
        0 | 1
            | 2
            | 5
            | 6
            | 7
            | 8
            | 9
            | 10
            | 14
            | 15
            | 16
            | 17
            | 18
            | 33
            | 34
            | 35
            | 36
            | 52
            | 53
            | 54
            | 85
            | 86
    )
}

fn is_hypervalent_non_metal(
    topology: &TopologyBlock,
    atom_index: usize,
) -> Result<bool, CleanupError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: isHypervalentNonMetal
    // RDKit✔️✔️: bool isHypervalentNonMetal(Atom *atom) {
    // RDKit✔️✔️:   if (QueryOps::isMetal(*atom)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    let atom = &topology.atoms[atom_index];
    if is_metal(atom.atomic_number()) {
        return Ok(false);
    }
    // RDKit✔️✔️:   atom->updatePropertyCache(false);
    // RDKit✔️✔️:   int ev = atom->getValence(Atom::ValenceType::EXPLICIT);
    let atom_id = AtomId::new(atom_index);
    let explicit_valence =
        calculate_explicit_valence_for_topology(topology, atom_id, false, false)?;
    // RDKit✔️✔️:   // Check the explicit valence of the non-metal against the allowed
    // RDKit✔️✔️:   // valences of the atom, adjusted by its formal charge.  This means that
    // RDKit✔️✔️:   // N+ is treated the same as C, O+ the same as N.  This allows for,
    // RDKit✔️✔️:   // for example, c1cccc[n+]1-[Fe] to be acceptable and not turned into
    // RDKit✔️✔️:   // c1cccc[n+]1->[Fe].  After all, c1cccc[n+]1-C is ok.  Although this is
    // RDKit✔️✔️:   // a poor example because c1ccccn1->[Fe] appears to be the normal
    // RDKit✔️✔️:   // way that pyridine complexes with transition metals.  Heme b in
    // RDKit✔️✔️:   // CHEBI:26355 is an example of when this is required.
    // RDKit✔️✔️:   int effAtomicNum = atom->getAtomicNum() - atom->getFormalCharge();
    let effective_atomic_number = i32::from(atom.atomic_number()) - i32::from(atom.formal_charge());
    // RDKit✔️✔️:   if (effAtomicNum <= 0) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if effective_atomic_number <= 0 {
        return Ok(false);
    }
    let effective_atomic_number =
        u8::try_from(effective_atomic_number).map_err(|_| ValenceError::PeriodicTableLookup {
            atomic_number: atom.atomic_number(),
            field: "effective atomic number valences",
        })?;
    // RDKit✔️✔️:   // atom is a non-metal. If its explicit valence is greater than the
    // RDKit✔️✔️:   // maximum allowed valence then it is hypervalent.
    // RDKit✔️✔️:   //  We have a special case in here for aromatic atoms where the explicit
    // RDKit✔️✔️:   //  valence matches the max allowed and the degree is 4. This is there for
    // RDKit✔️✔️:   //  cases like cyclopentadienyl - metal systems. We need this special case
    // RDKit✔️✔️:   //  because the explicit valence on the C atoms there ends up being 4
    // RDKit✔️✔️:   const auto &otherValens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(effAtomicNum);
    let valences =
        rdkit_valence_list(effective_atomic_number)?.ok_or(ValenceError::PeriodicTableLookup {
            atomic_number: effective_atomic_number,
            field: "valences",
        })?;
    // RDKit✔️✔️:   auto maxV = otherValens.back();
    let maximum_valence = *valences.last().ok_or(ValenceError::PeriodicTableLookup {
        atomic_number: effective_atomic_number,
        field: "valences",
    })?;
    let implicit_valence =
        calculate_implicit_valence_for_topology(topology, atom_id, explicit_valence, false, false)?;
    let total_degree_is_four = implicit_valence >= 0
        && topology
            .adjacency
            .neighbors_of(atom_index)
            .len()
            .checked_add(usize::from(atom.explicit_hydrogens()))
            .and_then(|degree| degree.checked_add(implicit_valence as usize))
            == Some(4);
    // RDKit✔️✔️:   if (maxV > 0 && (ev > maxV || (ev == maxV && atom->getIsAromatic() &&
    // RDKit✔️✔️:                                  atom->getTotalDegree() == 4))) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    if maximum_valence > 0
        && (explicit_valence > maximum_valence
            || (explicit_valence == maximum_valence && atom.is_aromatic() && total_degree_is_four))
    {
        return Ok(true);
    }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: isHypervalentNonMetal
    Ok(false)
}

fn number_of_dative_bonds(topology: &TopologyBlock, atom_index: usize) -> usize {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: numDativeBonds
    // RDKit✔️✔️: int numDativeBonds(const Atom *atom) {
    // RDKit✔️✔️:   int numDatives = 0;
    // RDKit✔️✔️:   auto &mol = atom->getOwningMol();
    // RDKit✔️✔️:   for (auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::BondType::DATIVE ||
    // RDKit✔️✔️:         bond->getBondType() == Bond::BondType::DATIVEONE ||
    // RDKit✔️✔️:         bond->getBondType() == Bond::BondType::DATIVEL ||
    // RDKit✔️✔️:         bond->getBondType() == Bond::BondType::DATIVER) {
    // RDKit✔️✔️:       ++numDatives;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return numDatives;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: numDativeBonds
    topology
        .adjacency
        .neighbors_of(atom_index)
        .iter()
        .filter(|neighbor| {
            matches!(
                topology.bonds[neighbor.bond.index()].order(),
                BondOrder::Dative
                    | BondOrder::DativeOne
                    | BondOrder::DativeLeft
                    | BondOrder::DativeRight
            )
        })
        .count()
}

fn no_dative(atomic_number: u8) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: noDative
    // RDKit✔️✔️: bool noDative(const Atom *a) {
    // RDKit✔️✔️:   static const std::set<int> noD{1, 2, 9, 10};
    // RDKit✔️✔️:   return (noD.find(a->getAtomicNum()) != noD.end());
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: noDative
    matches!(atomic_number, 1 | 2 | 9 | 10)
}

fn cleanup_metal_bond(
    topology: &mut TopologyBlock,
    atom_index: usize,
    ranks: &[usize],
) -> Result<(), CleanupError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: metalBondCleanup
    // RDKit✔️✔️: void metalBondCleanup(RWMol &mol, Atom *atom,
    // RDKit✔️✔️:                       const std::vector<unsigned int> &ranks) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom in metalBondCleanup");
    // RDKit✔️✔️:   // The IUPAC recommendation for ligand->metal coordination bonds is that
    // RDKit✔️✔️:   // they be single.  This upsets the RDKit valence model, as seen in
    // RDKit✔️✔️:   // CHEBI:26355, heme b.  If the valence of a non-metal atom is above the
    // RDKit✔️✔️:   // maximum in the RDKit model, and there are single bonds from it to metal
    // RDKit✔️✔️:   // change those bonds to atom->metal dative.
    // RDKit✔️✔️:   // If the atom is bonded to more than 1 metal atom, choose the one
    // RDKit✔️✔️:   // with the fewer dative bonds incident on it, with the canonical
    // RDKit✔️✔️:   // rank of the atoms as a tie-breaker.
    // RDKit✔️✔️:   if (isHypervalentNonMetal(atom) && !noDative(atom)) {
    if is_hypervalent_non_metal(topology, atom_index)?
        && !no_dative(topology.atoms[atom_index].atomic_number())
    {
        // RDKit✔️✔️:     std::vector<Atom *> metals;
        let mut metals = Vec::new();
        // RDKit✔️✔️:     // see if there are any metals bonded to it by a single bond
        // RDKit✔️✔️:     for (auto bond : mol.atomBonds(atom)) {
        // RDKit✔️✔️:       if (bond->getBondType() == Bond::BondType::SINGLE &&
        // RDKit✔️✔️:           QueryOps::isMetal(*bond->getOtherAtom(atom))) {
        // RDKit✔️✔️:         metals.push_back(bond->getOtherAtom(atom));
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        for neighbor in topology.adjacency.neighbors_of(atom_index) {
            if topology.bonds[neighbor.bond.index()].order() == BondOrder::Single
                && is_metal(topology.atoms[neighbor.atom_index].atomic_number())
            {
                metals.push(neighbor.atom_index);
            }
        }
        // RDKit✔️✔️:     if (!metals.empty()) {
        if !metals.is_empty() {
            // RDKit✔️✔️:       std::sort(metals.begin(), metals.end(),
            // RDKit✔️✔️:                 [&](const Atom *a1, const Atom *a2) -> bool {
            // RDKit✔️✔️:                   int nda1 = numDativeBonds(a1);
            // RDKit✔️✔️:                   int nda2 = numDativeBonds(a2);
            // RDKit✔️✔️:                   if (nda1 == nda2) {
            // RDKit✔️✔️:                     return ranks[a1->getIdx()] > ranks[a2->getIdx()];
            // RDKit✔️✔️:                   } else {
            // RDKit✔️✔️:                     return nda1 < nda2;
            // RDKit✔️✔️:                   }
            // RDKit✔️✔️:                 });
            metals.sort_by(|left, right| {
                let left_count = number_of_dative_bonds(topology, *left);
                let right_count = number_of_dative_bonds(topology, *right);
                match left_count.cmp(&right_count) {
                    Ordering::Equal => ranks[*right].cmp(&ranks[*left]),
                    ordering => ordering,
                }
            });
            // RDKit✔️✔️:       auto bond =
            // RDKit✔️✔️:           mol.getBondBetweenAtoms(atom->getIdx(), metals.front()->getIdx());
            let selected_metal = metals[0];
            let selected_bond = topology
                .adjacency
                .neighbors_of(atom_index)
                .iter()
                .find(|neighbor| neighbor.atom_index == selected_metal)
                .map(|neighbor| neighbor.bond);
            // RDKit✔️✔️:       if (bond) {
            if let Some(bond_id) = selected_bond {
                // RDKit✔️✔️:         bond->setBondType(RDKit::Bond::BondType::DATIVE);
                // RDKit✔️✔️:         bond->setBeginAtom(atom);
                // RDKit✔️✔️:         bond->setEndAtom(metals.front());
                let bond = &mut topology.bonds[bond_id.index()];
                bond.set_order(BondOrder::Dative);
                bond.set_endpoints(AtomId::new(atom_index), AtomId::new(selected_metal));
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: metalBondCleanup
    Ok(())
}

fn cleanup_organometallics(topology: &mut TopologyBlock) -> Result<(), CleanupError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: cleanUpOrganometallics
    // RDKit✔️✔️: void cleanUpOrganometallics(RWMol &mol) {
    // RDKit✔️✔️:   // At present all this does is look for single bonds between
    // RDKit✔️✔️:   // non-metals and metals where the non-metal exceeds one of
    // RDKit✔️✔️:   // its normal valence states, and replaces that bond with
    // RDKit✔️✔️:   // a dative one from the non-metal to the metal.
    // RDKit✔️✔️:   bool needsFixing = false;
    let mut needs_fixing = false;
    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    for atom_index in 0..topology.atoms.len() {
        // RDKit✔️✔️:     if (isHypervalentNonMetal(atom) && !noDative(atom)) {
        if is_hypervalent_non_metal(topology, atom_index)?
            && !no_dative(topology.atoms[atom_index].atomic_number())
        {
            // RDKit✔️✔️:       // see if there are any metals bonded to it by a single bond
            // RDKit✔️✔️:       for (auto bond : mol.atomBonds(atom)) {
            // RDKit✔️✔️:         if (bond->getBondType() == Bond::BondType::SINGLE &&
            // RDKit✔️✔️:             QueryOps::isMetal(*bond->getOtherAtom(atom))) {
            // RDKit✔️✔️:           needsFixing = true;
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            needs_fixing = topology
                .adjacency
                .neighbors_of(atom_index)
                .iter()
                .any(|neighbor| {
                    topology.bonds[neighbor.bond.index()].order() == BondOrder::Single
                        && is_metal(topology.atoms[neighbor.atom_index].atomic_number())
                });
        }
        // RDKit✔️✔️:     if (needsFixing) {
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     }
        if needs_fixing {
            break;
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   if (!needsFixing) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    if !needs_fixing {
        return Ok(());
    }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   mol.updatePropertyCache(false);
    // Canonical ranking below computes the same detached valence state without
    // installing a runtime cache.
    // RDKit✔️✔️:   // First see if anything needs doing
    // RDKit✔️✔️:   std::vector<unsigned int> ranks(mol.getNumAtoms());
    // RDKit✔️✔️:   RDKit::Canon::rankMolAtoms(mol, ranks);
    let ranks = crate::kekulize::rank_mol_atoms(topology)?;
    // RDKit✔️✔️:   std::vector<std::pair<int, int>> atom_ranks;
    // RDKit✔️✔️:   for (size_t i = 0; i < ranks.size(); ++i) {
    // RDKit✔️✔️:     atom_ranks.push_back(std::make_pair(i, ranks[i]));
    // RDKit✔️✔️:   }
    let mut atom_ranks = ranks.iter().copied().enumerate().collect::<Vec<_>>();
    // RDKit✔️✔️:   std::sort(atom_ranks.begin(), atom_ranks.end(),
    // RDKit✔️✔️:             [](const std::pair<int, int> &p1, std::pair<int, int> &p2) -> bool {
    // RDKit✔️✔️:               return p1.second < p2.second;
    // RDKit✔️✔️:             });
    atom_ranks.sort_by_key(|entry| entry.1);
    // RDKit✔️✔️:   for (auto ar : atom_ranks) {
    // RDKit✔️✔️:     auto atom = mol.getAtomWithIdx(ar.first);
    // RDKit✔️✔️:     metalBondCleanup(mol, atom, ranks);
    // RDKit✔️✔️:   }
    for (atom_index, _) in atom_ranks {
        cleanup_metal_bond(topology, atom_index, &ranks)?;
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: cleanUpOrganometallics
    Ok(())
}

#[cfg(test)]
mod tests {
    use cosmolkit_model::{AdjacencyList, Atom, AtomSpec, Bond, BondId, BondSpec};
    use cosmolkit_types::Element;

    use super::*;

    #[test]
    fn dative_variant_count_includes_all_four_source_orders() {
        for order in [
            BondOrder::Dative,
            BondOrder::DativeOne,
            BondOrder::DativeLeft,
            BondOrder::DativeRight,
        ] {
            let atoms = [26, 30]
                .into_iter()
                .enumerate()
                .map(|(index, atomic_number)| {
                    Atom::from_spec(
                        AtomId::new(index),
                        AtomSpec::new(
                            Element::from_atomic_number(atomic_number).expect("test element"),
                        ),
                    )
                })
                .collect::<Vec<_>>();
            let bonds = vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), order),
            )];
            let topology = TopologyBlock {
                adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
                atoms,
                bonds,
                ..TopologyBlock::default()
            };
            topology.validate().unwrap();

            assert_eq!(number_of_dative_bonds(&topology, 0), 1, "{order:?}");
            assert_eq!(number_of_dative_bonds(&topology, 1), 1, "{order:?}");
        }
    }
}
