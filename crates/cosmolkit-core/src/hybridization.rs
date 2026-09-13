//! RDKit-aligned hybridization assignment over detached topology values.

use cosmolkit_model::{Atom, AtomId, TopologyBlock, TopologyValidationError};
use cosmolkit_types::{BondOrder, ChiralTag, Hybridization};

use crate::{ValenceAssignment, ValenceError, periodic_table_outer_electrons};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HybridizationAssignment {
    pub values: Vec<Hybridization>,
}

#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum HybridizationError {
    #[error("invalid topology: {0}")]
    InvalidTopology(#[from] TopologyValidationError),
    #[error("valence assignment field {field} has {actual} rows; expected {expected}")]
    ValenceAssignmentLength {
        field: &'static str,
        actual: usize,
        expected: usize,
    },
    #[error("invalid valence row for atom {atom}: {field}={value}")]
    InvalidValenceRow {
        atom: AtomId,
        field: &'static str,
        value: i32,
    },
    #[error("integer overflow for atom {atom} while computing {field}")]
    IntegerOverflow { atom: AtomId, field: &'static str },
    #[error(transparent)]
    Valence(#[from] ValenceError),
}

fn validate_valence(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<(), HybridizationError> {
    let expected = topology.atoms.len();
    for (field, rows) in [
        ("explicit_valence", &valence.explicit_valence),
        ("implicit_hydrogens", &valence.implicit_hydrogens),
    ] {
        if rows.len() != expected {
            return Err(HybridizationError::ValenceAssignmentLength {
                field,
                actual: rows.len(),
                expected,
            });
        }
        if let Some((row, &value)) = rows.iter().enumerate().find(|(_, value)| **value < 0) {
            return Err(HybridizationError::InvalidValenceRow {
                atom: AtomId::new(row),
                field,
                value,
            });
        }
    }
    Ok(())
}

fn atom(topology: &TopologyBlock, atom_id: AtomId) -> &Atom {
    &topology.atoms[atom_id.index()]
}

fn overflow(atom: AtomId, field: &'static str) -> HybridizationError {
    HybridizationError::IntegerOverflow { atom, field }
}

fn total_degree(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom_id: AtomId,
) -> Result<i32, HybridizationError> {
    let adjacency_degree = i32::try_from(topology.adjacency.neighbors_of(atom_id.index()).len())
        .map_err(|_| overflow(atom_id, "total degree"))?;
    let hydrogens = i32::try_from(crate::hcount::total_hydrogen_count_from_validated(
        topology, valence, atom_id, false,
    )?)
    .map_err(|_| overflow(atom_id, "total degree"))?;
    adjacency_degree
        .checked_add(hydrogens)
        .ok_or_else(|| overflow(atom_id, "total degree"))
}

fn is_dative(order: BondOrder) -> bool {
    // BEGIN RDKIT CPP FUNCTION isDative
    // RDKit✔️✔️: inline bool isDative(const Bond::BondType bt) {
    // RDKit✔️✔️:   return bt == Bond::BondType::DATIVE || bt == Bond::BondType::DATIVEL ||
    // RDKit✔️✔️:          bt == Bond::BondType::DATIVER || bt == Bond::BondType::DATIVEONE;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isDative
    matches!(
        order,
        BondOrder::Dative | BondOrder::DativeLeft | BondOrder::DativeRight | BondOrder::DativeOne
    )
}

fn num_bonds_plus_lone_pairs(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom_id: AtomId,
) -> Result<i32, HybridizationError> {
    // BEGIN RDKIT CPP FUNCTION numBondsPlusLonePairs
    // RDKit✔️✔️: int numBondsPlusLonePairs(Atom *at) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad atom");
    // RDKit✔️✔️:   int deg = at->getTotalDegree();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto &mol = at->getOwningMol();
    // RDKit✔️✔️:   for (const auto bond : mol.atomBonds(at)) {
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::ZERO ||
    // RDKit✔️✔️:         (isDative(*bond) && at->getIdx() != bond->getEndAtomIdx())) {
    // RDKit✔️✔️:       --deg;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (at->getAtomicNum() <= 1) {
    // RDKit✔️✔️:     return deg;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   int nouter = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
    // RDKit✔️✔️:   int totalValence = at->getTotalValence();
    // RDKit✔️✔️:   int chg = at->getFormalCharge();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int numFreeElectrons = nouter - (totalValence + chg);
    // RDKit✔️✔️:   if (totalValence + nouter - chg < 8) {
    // RDKit✔️✔️:     // we're below an octet, so we need to think
    // RDKit✔️✔️:     // about radicals:
    // RDKit✔️✔️:     int numRadicals = at->getNumRadicalElectrons();
    // RDKit✔️✔️:     int numLonePairs = (numFreeElectrons - numRadicals) / 2;
    // RDKit✔️✔️:     return deg + numLonePairs + numRadicals;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     int numLonePairs = numFreeElectrons / 2;
    // RDKit✔️✔️:     return deg + numLonePairs;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION numBondsPlusLonePairs
    let mut degree = total_degree(topology, valence, atom_id)?;
    for neighbor in topology.adjacency.neighbors_of(atom_id.index()) {
        let bond = &topology.bonds[neighbor.bond.index()];
        if bond.order() == BondOrder::Zero || (is_dative(bond.order()) && atom_id != bond.end()) {
            degree = degree
                .checked_sub(1)
                .ok_or_else(|| overflow(atom_id, "dative-adjusted degree"))?;
        }
    }

    let current = atom(topology, atom_id);
    if current.atomic_number() <= 1 {
        return Ok(degree);
    }
    let outer = periodic_table_outer_electrons(current.atomic_number())?;
    let total_valence = valence.explicit_valence[atom_id.index()]
        .checked_add(valence.implicit_hydrogens[atom_id.index()])
        .ok_or_else(|| overflow(atom_id, "total valence"))?;
    let charge = i32::from(current.formal_charge());
    let free_electrons = outer
        .checked_sub(
            total_valence
                .checked_add(charge)
                .ok_or_else(|| overflow(atom_id, "free electrons"))?,
        )
        .ok_or_else(|| overflow(atom_id, "free electrons"))?;
    let octet_sum = total_valence
        .checked_add(outer)
        .and_then(|value| value.checked_sub(charge))
        .ok_or_else(|| overflow(atom_id, "octet test"))?;
    if octet_sum < 8 {
        let radicals = i32::from(current.radical_electrons());
        let lone_pairs = free_electrons
            .checked_sub(radicals)
            .ok_or_else(|| overflow(atom_id, "below-octet lone pairs"))?
            / 2;
        degree
            .checked_add(lone_pairs)
            .and_then(|value| value.checked_add(radicals))
            .ok_or_else(|| overflow(atom_id, "below-octet orbitals"))
    } else {
        degree
            .checked_add(free_electrons / 2)
            .ok_or_else(|| overflow(atom_id, "octet orbitals"))
    }
}

pub fn assign_hybridization(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<HybridizationAssignment, HybridizationError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::setHybridization
    // RDKit✔️❌: void setHybridization(ROMol &mol) {
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     if (atom->getAtomicNum() == 0) {
    // RDKit✔️❌:       atom->setHybridization(Atom::UNSPECIFIED);
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       // if the stereo spec matches the coordination number, this is easy
    // RDKit✔️❌:       switch (atom->getChiralTag()) {
    // RDKit✔️❌:         case Atom::ChiralType::CHI_TETRAHEDRAL:
    // RDKit✔️❌:         case Atom::ChiralType::CHI_TETRAHEDRAL_CW:
    // RDKit✔️❌:         case Atom::ChiralType::CHI_TETRAHEDRAL_CCW:
    // RDKit✔️❌:           if (atom->getTotalDegree() == 4) {
    // RDKit✔️❌:             atom->setHybridization(Atom::HybridizationType::SP3);
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case Atom::ChiralType::CHI_SQUAREPLANAR:
    // RDKit✔️❌:           if (atom->getTotalDegree() <= 4 && atom->getTotalDegree() >= 2) {
    // RDKit✔️❌:             atom->setHybridization(Atom::HybridizationType::SP2D);
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
    // RDKit✔️❌:           if (atom->getTotalDegree() <= 5 && atom->getTotalDegree() >= 2) {
    // RDKit✔️❌:             atom->setHybridization(Atom::HybridizationType::SP3D);
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case Atom::ChiralType::CHI_OCTAHEDRAL:
    // RDKit✔️❌:           if (atom->getTotalDegree() <= 6 && atom->getTotalDegree() >= 2) {
    // RDKit✔️❌:             atom->setHybridization(Atom::HybridizationType::SP3D2);
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         default:
    // RDKit✔️❌:           break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       // otherwise we have to do some work
    // RDKit✔️❌:       int norbs;
    // RDKit✔️❌:       // try to be smart for early elements, but for later
    // RDKit✔️❌:       // ones just use the degree
    // RDKit✔️❌:       // FIX: we should probably also be using the degree for metals
    // RDKit✔️❌:       if (atom->getAtomicNum() < 89) {
    // RDKit✔️❌:         norbs = numBondsPlusLonePairs(atom);
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         norbs = atom->getTotalDegree();
    // RDKit✔️❌:       }
    // RDKit✔️❌:       switch (norbs) {
    // RDKit✔️❌:         case 0:
    // RDKit✔️❌:           // This occurs for things like Na+
    // RDKit✔️❌:           atom->setHybridization(Atom::S);
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 1:
    // RDKit✔️❌:           atom->setHybridization(Atom::S);
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 2:
    // RDKit✔️❌:           atom->setHybridization(Atom::SP);
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 3:
    // RDKit✔️❌:           atom->setHybridization(Atom::SP2);
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 4:
    // RDKit✔️❌:           // potentially SP3, but we'll set it down to SP2
    // RDKit✔️❌:           // if we have a conjugated bond (like the second O
    // RDKit✔️❌:           // in O=CO)
    // RDKit✔️❌:           // we'll also avoid setting the hybridization down to
    // RDKit✔️❌:           // SP2 in the case of an atom with degree higher than 3
    // RDKit✔️❌:           // (e.g. things like CP1(C)=CC=CN=C1C, where the P
    // RDKit✔️❌:           //   has norbs = 4, and a conjugated bond, but clearly should
    // RDKit✔️❌:           //   not be SP2)
    // RDKit✔️❌:           // This is Issue276
    // RDKit✔️❌:           if (atom->getTotalDegree() > 3 ||
    // RDKit✔️❌:               !MolOps::atomHasConjugatedBond(atom)) {
    // RDKit✔️❌:             atom->setHybridization(Atom::SP3);
    // RDKit✔️❌:           } else {
    // RDKit✔️❌:             atom->setHybridization(Atom::SP2);
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 5:
    // RDKit✔️❌:           atom->setHybridization(Atom::SP3D);
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 6:
    // RDKit✔️❌:           atom->setHybridization(Atom::SP3D2);
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         default:
    // RDKit✔️❌:           atom->setHybridization(Atom::UNSPECIFIED);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::setHybridization
    // Returning a detached assignment keeps source mutation outside this
    // algorithm owner but allocates O(V) output storage that RDKit does not.
    topology.validate()?;
    validate_valence(topology, valence)?;
    let mut values = Vec::with_capacity(topology.atoms.len());
    for atom_index in 0..topology.atoms.len() {
        let atom_id = AtomId::new(atom_index);
        let current = atom(topology, atom_id);
        if current.atomic_number() == 0 {
            values.push(Hybridization::Unspecified);
            continue;
        }

        let degree = total_degree(topology, valence, atom_id)?;
        let tagged = match current.chiral_tag() {
            ChiralTag::Tetrahedral | ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
                if degree == 4 =>
            {
                Some(Hybridization::Sp3)
            }
            ChiralTag::SquarePlanar if (2..=4).contains(&degree) => Some(Hybridization::Sp2d),
            ChiralTag::TrigonalBipyramidal if (2..=5).contains(&degree) => {
                Some(Hybridization::Sp3d)
            }
            ChiralTag::Octahedral if (2..=6).contains(&degree) => Some(Hybridization::Sp3d2),
            _ => None,
        };
        if let Some(value) = tagged {
            values.push(value);
            continue;
        }

        let orbitals = if current.atomic_number() < 89 {
            num_bonds_plus_lone_pairs(topology, valence, atom_id)?
        } else {
            degree
        };
        let value = match orbitals {
            0 | 1 => Hybridization::S,
            2 => Hybridization::Sp,
            3 => Hybridization::Sp2,
            4 => {
                if degree > 3
                    || !crate::conjugation::atom_has_conjugated_bond_from_validated(
                        topology, atom_id,
                    )
                {
                    Hybridization::Sp3
                } else {
                    Hybridization::Sp2
                }
            }
            5 => Hybridization::Sp3d,
            6 => Hybridization::Sp3d2,
            _ => Hybridization::Unspecified,
        };
        values.push(value);
    }
    Ok(HybridizationAssignment { values })
}
