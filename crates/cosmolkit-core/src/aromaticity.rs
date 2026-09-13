//! RDKit-aligned aromaticity assignment over detached topology values.
//!
//! This module owns chemistry only. It never accepts a live `Molecule`, an
//! operation context, or runtime cache authority.

use std::collections::BTreeSet;

use cosmolkit_model::{Atom, AtomId, Bond, BondId, TopologyBlock, TopologyValidationError};
use cosmolkit_types::{BondOrder, Hybridization};

use crate::{
    KekulizeError, RingFindingError, RingInfo, ValenceAssignment, ValenceError,
    bond_valence_contrib, get_effective_atomic_num, periodic_table_more_electronegative,
    periodic_table_outer_electrons, rdkit_default_valence,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AromaticityModel {
    Rdkit,
    Simple,
    Mdl,
    Mmff94,
    Custom,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AromaticityParams {
    pub model: AromaticityModel,
}

impl Default for AromaticityParams {
    fn default() -> Self {
        Self {
            model: AromaticityModel::Rdkit,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AromaticityAssignment {
    pub topology: TopologyBlock,
    pub aromatic_ring_count: usize,
}

#[derive(Clone, Debug, Eq, PartialEq, thiserror::Error)]
pub enum AromaticityError {
    #[error("invalid topology: {0}")]
    InvalidTopology(#[from] TopologyValidationError),
    #[error("ring information is not initialized")]
    RingInfoNotInitialized,
    #[error(
        "ring information dimensions ({ring_atom_count} atoms, {ring_bond_count} bonds) do not match topology ({topology_atom_count} atoms, {topology_bond_count} bonds)"
    )]
    RingInfoDimensionMismatch {
        ring_atom_count: usize,
        ring_bond_count: usize,
        topology_atom_count: usize,
        topology_bond_count: usize,
    },
    #[error(
        "ring atom/bond table length mismatch: {atom_rows} atom rows and {bond_rows} bond rows"
    )]
    RingTableLengthMismatch { atom_rows: usize, bond_rows: usize },
    #[error(
        "ring {ring_index} atom/bond length mismatch: {atom_count} atoms and {bond_count} bonds"
    )]
    RingRowLengthMismatch {
        ring_index: usize,
        atom_count: usize,
        bond_count: usize,
    },
    #[error("ring {ring_index} references atom {atom}, outside {atom_count} topology atoms")]
    RingAtomOutOfRange {
        ring_index: usize,
        atom: AtomId,
        atom_count: usize,
    },
    #[error("ring {ring_index} references bond {bond}, outside {bond_count} topology bonds")]
    RingBondOutOfRange {
        ring_index: usize,
        bond: BondId,
        bond_count: usize,
    },
    #[error("atom {atom} is outside {atom_count} topology atoms")]
    AtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("bond {bond} is outside {bond_count} topology bonds")]
    BondOutOfRange { bond: BondId, bond_count: usize },
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
    #[error("expected ring bond not found between atom {begin} and atom {end}")]
    ExpectedRingBondNotFound { begin: AtomId, end: AtomId },
    #[error("unsupported aromaticity model {model:?}: {detail}")]
    UnsupportedModel {
        model: AromaticityModel,
        detail: &'static str,
    },
    #[error("unsupported aromaticity state at {field}: {detail}")]
    UnsupportedState {
        field: &'static str,
        detail: &'static str,
    },
    #[error("integer overflow while computing {field}")]
    IntegerOverflow { field: &'static str },
    #[error(
        "aromaticity changed topology shape from {input_atoms} atoms/{input_bonds} bonds to {output_atoms} atoms/{output_bonds} bonds"
    )]
    UnexpectedTopologyShape {
        input_atoms: usize,
        input_bonds: usize,
        output_atoms: usize,
        output_bonds: usize,
    },
    #[error("aromaticity changed atom identity at row {row} from {expected} to {actual}")]
    AtomIdentityChanged {
        row: usize,
        expected: AtomId,
        actual: AtomId,
    },
    #[error("aromaticity changed bond identity at row {row} from {expected} to {actual}")]
    BondIdentityChanged {
        row: usize,
        expected: BondId,
        actual: BondId,
    },
    #[error(
        "aromaticity model {model:?} returned {actual} aromatic rings; maximum source count is {maximum}"
    )]
    AromaticRingCountOutOfRange {
        model: AromaticityModel,
        actual: usize,
        maximum: usize,
    },
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
    #[error(transparent)]
    Kekulize(#[from] KekulizeError),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ElectronDonorType {
    Vacant,
    One,
    Two,
    OneOrTwo,
    Any,
    None,
}

#[derive(Debug, Clone, Copy)]
struct CandidateOptions {
    allow_third_row: bool,
    allow_triple_bonds: bool,
    allow_higher_exceptions: bool,
    only_c_or_n: bool,
    allow_exocyclic_multiple_bonds: bool,
}

impl Default for CandidateOptions {
    fn default() -> Self {
        Self {
            allow_third_row: true,
            allow_triple_bonds: true,
            allow_higher_exceptions: true,
            only_c_or_n: false,
            allow_exocyclic_multiple_bonds: true,
        }
    }
}

fn validate_inputs(topology: &TopologyBlock, rings: &RingInfo) -> Result<(), AromaticityError> {
    topology.validate()?;
    if !rings.is_initialized() {
        return Err(AromaticityError::RingInfoNotInitialized);
    }
    if rings.atom_row_count() != topology.atoms.len()
        || rings.bond_row_count() != topology.bonds.len()
    {
        return Err(AromaticityError::RingInfoDimensionMismatch {
            ring_atom_count: rings.atom_row_count(),
            ring_bond_count: rings.bond_row_count(),
            topology_atom_count: topology.atoms.len(),
            topology_bond_count: topology.bonds.len(),
        });
    }
    if rings.atom_rings().len() != rings.bond_rings().len() {
        return Err(AromaticityError::RingTableLengthMismatch {
            atom_rows: rings.atom_rings().len(),
            bond_rows: rings.bond_rings().len(),
        });
    }
    for (ring_index, (atoms, bonds)) in rings
        .atom_rings()
        .iter()
        .zip(rings.bond_rings())
        .enumerate()
    {
        if atoms.len() != bonds.len() {
            return Err(AromaticityError::RingRowLengthMismatch {
                ring_index,
                atom_count: atoms.len(),
                bond_count: bonds.len(),
            });
        }
        for &atom in atoms {
            if atom.index() >= topology.atoms.len() {
                return Err(AromaticityError::RingAtomOutOfRange {
                    ring_index,
                    atom,
                    atom_count: topology.atoms.len(),
                });
            }
        }
        for &bond in bonds {
            if bond.index() >= topology.bonds.len() {
                return Err(AromaticityError::RingBondOutOfRange {
                    ring_index,
                    bond,
                    bond_count: topology.bonds.len(),
                });
            }
        }
    }
    Ok(())
}

fn validate_valence(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<(), AromaticityError> {
    for (field, actual) in [
        ("explicit_valence", valence.explicit_valence.len()),
        ("implicit_hydrogens", valence.implicit_hydrogens.len()),
    ] {
        if actual != topology.atoms.len() {
            return Err(AromaticityError::ValenceAssignmentLength {
                field,
                actual,
                expected: topology.atoms.len(),
            });
        }
    }
    for (index, (&explicit, &implicit)) in valence
        .explicit_valence
        .iter()
        .zip(&valence.implicit_hydrogens)
        .enumerate()
    {
        if explicit < 0 {
            return Err(AromaticityError::InvalidValenceRow {
                atom: AtomId::new(index),
                field: "explicit_valence",
                value: explicit,
            });
        }
        if implicit < 0 {
            return Err(AromaticityError::InvalidValenceRow {
                atom: AtomId::new(index),
                field: "implicit_hydrogens",
                value: implicit,
            });
        }
    }
    Ok(())
}

fn atom<'a>(topology: &'a TopologyBlock, id: AtomId) -> Result<&'a Atom, AromaticityError> {
    topology
        .atoms
        .get(id.index())
        .ok_or(AromaticityError::AtomOutOfRange {
            atom: id,
            atom_count: topology.atoms.len(),
        })
}

fn bond<'a>(topology: &'a TopologyBlock, id: BondId) -> Result<&'a Bond, AromaticityError> {
    topology
        .bonds
        .get(id.index())
        .ok_or(AromaticityError::BondOutOfRange {
            bond: id,
            bond_count: topology.bonds.len(),
        })
}

fn valence_row(valence: &ValenceAssignment, atom: AtomId) -> Result<(i32, i32), AromaticityError> {
    let explicit = *valence.explicit_valence.get(atom.index()).ok_or(
        AromaticityError::ValenceAssignmentLength {
            field: "explicit_valence",
            actual: valence.explicit_valence.len(),
            expected: atom.index().saturating_add(1),
        },
    )?;
    let implicit = *valence.implicit_hydrogens.get(atom.index()).ok_or(
        AromaticityError::ValenceAssignmentLength {
            field: "implicit_hydrogens",
            actual: valence.implicit_hydrogens.len(),
            expected: atom.index().saturating_add(1),
        },
    )?;
    Ok((explicit, implicit))
}

fn incident_noncyclic_multiple_bond(
    topology: &TopologyBlock,
    rings: &RingInfo,
    atom_id: AtomId,
) -> Result<Option<AtomId>, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION incidentNonCyclicMultipleBond
    // RDKit✔️✔️: bool incidentNonCyclicMultipleBond(const Atom *at, int &who) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad atom");
    // RDKit✔️✔️:   // check if "at" has an non-cyclic multiple bond on it
    // RDKit✔️✔️:   // if yes check which atom this bond goes to
    // RDKit✔️✔️:   // and record the atomID in who
    // RDKit✔️✔️:   const auto &mol = at->getOwningMol();
    // RDKit✔️✔️:   for (const auto bond : mol.atomBonds(at)) {
    // RDKit✔️✔️:     if (!mol.getRingInfo()->numBondRings(bond->getIdx())) {
    // RDKit✔️✔️:       if (bond->getValenceContrib(at) >= 2.0) {
    // RDKit✔️✔️:         who = bond->getOtherAtomIdx(at->getIdx());
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION incidentNonCyclicMultipleBond
    atom(topology, atom_id)?;
    for neighbor in topology.adjacency.neighbors_of(atom_id.index()) {
        let candidate = bond(topology, neighbor.bond)?;
        if rings.num_bond_rings(neighbor.bond) == 0
            && bond_valence_contrib(candidate, atom_id)? >= 2.0
        {
            return Ok(Some(AtomId::new(neighbor.atom_index)));
        }
    }
    Ok(None)
}

fn incident_cyclic_multiple_bond(
    topology: &TopologyBlock,
    rings: &RingInfo,
    atom_id: AtomId,
) -> Result<bool, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION incidentCyclicMultipleBond
    // RDKit✔️✔️: bool incidentCyclicMultipleBond(const Atom *at) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad atom");
    // RDKit✔️✔️:   const auto &mol = at->getOwningMol();
    // RDKit✔️✔️:   for (const auto bond : mol.atomBonds(at)) {
    // RDKit✔️✔️:     if (mol.getRingInfo()->numBondRings(bond->getIdx())) {
    // RDKit✔️✔️:       if (bond->getValenceContrib(at) >= 2.0) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION incidentCyclicMultipleBond
    atom(topology, atom_id)?;
    for neighbor in topology.adjacency.neighbors_of(atom_id.index()) {
        let candidate = bond(topology, neighbor.bond)?;
        if rings.num_bond_rings(neighbor.bond) != 0
            && bond_valence_contrib(candidate, atom_id)? >= 2.0
        {
            return Ok(true);
        }
    }
    Ok(false)
}

fn incident_multiple_bond(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom_id: AtomId,
) -> Result<bool, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION incidentMultipleBond
    // RDKit✔️✔️: bool incidentMultipleBond(const Atom *at) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad atom");
    // RDKit✔️✔️:   const auto &mol = at->getOwningMol();
    // RDKit✔️✔️:   auto deg = at->getDegree() + at->getNumExplicitHs();
    // RDKit✔️✔️:   for (const auto bond : mol.atomBonds(at)) {
    // RDKit✔️✔️:     if (!std::lround(bond->getValenceContrib(at))) {
    // RDKit✔️✔️:       --deg;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return at->getValence(Atom::ValenceType::EXPLICIT) != deg;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION incidentMultipleBond
    let atom = atom(topology, atom_id)?;
    let mut degree = i32::try_from(topology.adjacency.neighbors_of(atom_id.index()).len())
        .map_err(|_| AromaticityError::IntegerOverflow {
            field: "atom degree",
        })?
        .checked_add(i32::from(atom.explicit_hydrogens()))
        .ok_or(AromaticityError::IntegerOverflow {
            field: "atom degree",
        })?;
    for neighbor in topology.adjacency.neighbors_of(atom_id.index()) {
        if bond_valence_contrib(bond(topology, neighbor.bond)?, atom_id)?.round() == 0.0 {
            degree = degree
                .checked_sub(1)
                .ok_or(AromaticityError::IntegerOverflow {
                    field: "atom degree",
                })?;
        }
    }
    Ok(valence_row(valence, atom_id)?.0 != degree)
}

fn is_bond_order_query(_bond: &Bond) -> bool {
    // BEGIN RDKIT CPP FUNCTION isBondOrderQuery
    // RDKit✔️✔️: bool isBondOrderQuery(const Bond *bond) {
    // RDKit✔️✔️:   if (bond->hasQuery()) {
    // RDKit✔️✔️:     auto q = dynamic_cast<const QueryBond *>(bond)->getQuery();
    // RDKit✔️✔️:     // complex bond type queries are also bond order queries!
    // RDKit✔️✔️:     if (q->getTypeLabel() == "BondOrder" ||
    // RDKit✔️✔️:         QueryOps::hasComplexBondTypeQuery(*q)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isBondOrderQuery
    // `TopologyBlock` contains concrete `Bond` values. Query bonds belong to
    // `QueryGraph`, so the source `hasQuery()` branch is absent from this
    // function's modeled state space rather than silently approximated.
    false
}

pub(crate) fn count_atom_electrons(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom_id: AtomId,
) -> Result<i32, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::countAtomElec
    // RDKit✔️✔️: int countAtomElec(const Atom *at) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad atom");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // default valence :
    // RDKit✔️✔️:   auto dv = PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum());
    // RDKit✔️✔️:   if (dv <= 1) {
    // RDKit✔️✔️:     // univalent elements can't be either aromatic or conjugated
    // RDKit✔️✔️:     return -1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // total atom degree:
    // RDKit✔️✔️:   int degree = at->getDegree() + at->getTotalNumHs();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const auto &mol = at->getOwningMol();
    // RDKit✔️✔️:   for (const auto bond : mol.atomBonds(at)) {
    // RDKit✔️✔️:     // don't count bonds that aren't actually contributing to the valence here:
    // RDKit✔️✔️:     // if the bond is "real" (not undefined or zero), it always contributes to
    // RDKit✔️✔️:     // valence/degree, and in case the bond is a query bond with no order, we
    // RDKit✔️✔️:     // still need to check if the query is a bond query
    // RDKit✔️✔️:     if (!static_cast<Bond *>(bond)->getValenceContrib(at) &&
    // RDKit✔️✔️:         !isBondOrderQuery(bond)) {
    // RDKit✔️✔️:       --degree;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // if we are more than 3 coordinated we should not be aromatic
    // RDKit✔️✔️:   if (degree > 3) {
    // RDKit✔️✔️:     return -1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // number of lone pair electrons = (outer shell elecs) - (default valence)
    // RDKit✔️✔️:   auto nlp = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum()) - dv;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // subtract the charge to get the true number of lone pair electrons:
    // RDKit✔️✔️:   nlp = std::max(nlp - at->getFormalCharge(), 0);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int nRadicals = at->getNumRadicalElectrons();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // num electrons available for donation into the pi system:
    // RDKit✔️✔️:   int res = (dv - degree) + nlp - nRadicals;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (res > 1) {
    // RDKit✔️✔️:     // if we have an incident bond with order higher than 2,
    // RDKit✔️✔️:     // (e.g. triple or higher), we only want to return 1 electron
    // RDKit✔️✔️:     // we detect this using the total unsaturation, because we
    // RDKit✔️✔️:     // know that there aren't multiple unsaturations (detected
    // RDKit✔️✔️:     // above in isAtomCandForArom())
    // RDKit✔️✔️:     int nUnsaturations =
    // RDKit✔️✔️:         at->getValence(Atom::ValenceType::EXPLICIT) - at->getDegree();
    // RDKit✔️✔️:     if (nUnsaturations > 1) {
    // RDKit✔️✔️:       res = 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::countAtomElec
    let atom = atom(topology, atom_id)?;
    let default_valence = rdkit_default_valence(atom.atomic_number())?;
    if default_valence <= 1 {
        return Ok(-1);
    }
    let (explicit_valence, implicit_hydrogens) = valence_row(valence, atom_id)?;
    let mut degree = i32::try_from(topology.adjacency.neighbors_of(atom_id.index()).len())
        .map_err(|_| AromaticityError::IntegerOverflow {
            field: "atom degree",
        })?
        .checked_add(i32::from(atom.explicit_hydrogens()))
        .and_then(|value| value.checked_add(implicit_hydrogens))
        .ok_or(AromaticityError::IntegerOverflow {
            field: "atom degree",
        })?;
    for neighbor in topology.adjacency.neighbors_of(atom_id.index()) {
        let incident = bond(topology, neighbor.bond)?;
        if bond_valence_contrib(incident, atom_id)? == 0.0 && !is_bond_order_query(incident) {
            degree = degree
                .checked_sub(1)
                .ok_or(AromaticityError::IntegerOverflow {
                    field: "atom degree",
                })?;
        }
    }
    if degree > 3 {
        return Ok(-1);
    }
    let lone_pairs = (periodic_table_outer_electrons(atom.atomic_number())?
        - default_valence
        - i32::from(atom.formal_charge()))
    .max(0);
    let mut result = default_valence
        .checked_sub(degree)
        .and_then(|value| value.checked_add(lone_pairs))
        .and_then(|value| value.checked_sub(i32::from(atom.radical_electrons())))
        .ok_or(AromaticityError::IntegerOverflow {
            field: "available pi electrons",
        })?;
    if result > 1 {
        let graph_degree = i32::try_from(topology.adjacency.neighbors_of(atom_id.index()).len())
            .map_err(|_| AromaticityError::IntegerOverflow {
                field: "atom degree",
            })?;
        if explicit_valence - graph_degree > 1 {
            result = 1;
        }
    }
    Ok(result)
}

fn is_atom_candidate(
    topology: &TopologyBlock,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    atom_id: AtomId,
    donor: ElectronDonorType,
    options: CandidateOptions,
) -> Result<bool, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION isAtomCandForArom
    // RDKit✔️✔️: bool isAtomCandForArom(const Atom *at, const ElectronDonorType edon,
    // RDKit✔️✔️:                        bool allowThirdRow = true, bool allowTripleBonds = true,
    // RDKit✔️✔️:                        bool allowHigherExceptions = true, bool onlyCorN = false,
    // RDKit✔️✔️:                        bool allowExocyclicMultipleBonds = true) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad atom");
    // RDKit✔️✔️:   if (onlyCorN && at->getAtomicNum() != 6 && at->getAtomicNum() != 7) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!allowThirdRow && at->getAtomicNum() > 10) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // limit aromaticity to:
    // RDKit✔️✔️:   //   - the first two rows of the periodic table
    // RDKit✔️✔️:   //   - Se and Te
    // RDKit✔️✔️:   if (at->getAtomicNum() > 18 &&
    // RDKit✔️✔️:       (!allowHigherExceptions ||
    // RDKit✔️✔️:        (at->getAtomicNum() != 34 && at->getAtomicNum() != 52))) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   switch (edon) {
    // RDKit✔️✔️:     case VacantElectronDonorType:
    // RDKit✔️✔️:     case OneElectronDonorType:
    // RDKit✔️✔️:     case TwoElectronDonorType:
    // RDKit✔️✔️:     case OneOrTwoElectronDonorType:
    // RDKit✔️✔️:     case AnyElectronDonorType:
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       return (false);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // atoms that aren't in their default valence state also get shut out
    // RDKit✔️✔️:   auto defVal =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum());
    // RDKit✔️✔️:   if (defVal > 0 && rdcast<int>(at->getTotalValence()) >
    // RDKit✔️✔️:                         (PeriodicTable::getTable()->getDefaultValence(
    // RDKit✔️✔️:                             at->getAtomicNum() - at->getFormalCharge()))) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // heteroatoms or charged carbons with radicals also disqualify us from being
    // RDKit✔️✔️:   // considered. This was github issue 432 (heteroatoms) and 1936 (charged
    // RDKit✔️✔️:   // carbons)
    // RDKit✔️✔️:   if (at->getNumRadicalElectrons() &&
    // RDKit✔️✔️:       (at->getAtomicNum() != 6 || at->getFormalCharge())) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // We are going to explicitly disallow atoms that have more
    // RDKit✔️✔️:   // than one double or triple bond. This is to handle
    // RDKit✔️✔️:   // the situation:
    // RDKit✔️✔️:   //   C1=C=NC=N1 (sf.net bug 1934360)
    // RDKit✔️✔️:   int nUnsaturations =
    // RDKit✔️✔️:       at->getValence(Atom::ValenceType::EXPLICIT) - at->getDegree();
    // RDKit✔️✔️:   if (nUnsaturations > 1) {
    // RDKit✔️✔️:     unsigned int nMult = 0;
    // RDKit✔️✔️:     const auto &mol = at->getOwningMol();
    // RDKit✔️✔️:     for (const auto bond : mol.atomBonds(at)) {
    // RDKit✔️✔️:       switch (bond->getBondType()) {
    // RDKit✔️✔️:         case Bond::SINGLE:
    // RDKit✔️✔️:         case Bond::AROMATIC:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Bond::DOUBLE:
    // RDKit✔️✔️:           ++nMult;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Bond::TRIPLE:
    // RDKit✔️✔️:           if (!allowTripleBonds) {
    // RDKit✔️✔️:             return false;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           ++nMult;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           // hopefully we had enough sense that we don't even
    // RDKit✔️✔️:           // get here with these bonds... If we do land here,
    // RDKit✔️✔️:           // just bail... I have no good answer for them.
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (nMult > 1) {
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (nMult > 1) {
    // RDKit✔️✔️:       return (false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!allowExocyclicMultipleBonds) {
    // RDKit✔️✔️:     const auto &mol = at->getOwningMol();
    // RDKit✔️✔️:     for (const auto bond : mol.atomBonds(at)) {
    // RDKit✔️✔️:       if ((bond->getBondType() == Bond::DOUBLE ||
    // RDKit✔️✔️:            bond->getBondType() == Bond::TRIPLE) &&
    // RDKit✔️✔️:           !queryIsBondInRing(bond)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return (true);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isAtomCandForArom
    let atom = atom(topology, atom_id)?;
    if options.only_c_or_n && !matches!(atom.atomic_number(), 6 | 7) {
        return Ok(false);
    }
    if !options.allow_third_row && atom.atomic_number() > 10 {
        return Ok(false);
    }
    if atom.atomic_number() > 18
        && (!options.allow_higher_exceptions || !matches!(atom.atomic_number(), 34 | 52))
    {
        return Ok(false);
    }
    if donor == ElectronDonorType::None {
        return Ok(false);
    }
    let default_valence = rdkit_default_valence(atom.atomic_number())?;
    let (explicit, implicit) = valence_row(valence, atom_id)?;
    let total_valence =
        explicit
            .checked_add(implicit)
            .ok_or(AromaticityError::IntegerOverflow {
                field: "total valence",
            })?;
    let effective_atomic_number = get_effective_atomic_num(atom, true)?;
    if default_valence > 0 && total_valence > rdkit_default_valence(effective_atomic_number)? {
        return Ok(false);
    }
    if atom.radical_electrons() != 0 && (atom.atomic_number() != 6 || atom.formal_charge() != 0) {
        return Ok(false);
    }
    let graph_degree = i32::try_from(topology.adjacency.neighbors_of(atom_id.index()).len())
        .map_err(|_| AromaticityError::IntegerOverflow {
            field: "atom degree",
        })?;
    if explicit - graph_degree > 1 {
        let mut multiple_count = 0usize;
        for neighbor in topology.adjacency.neighbors_of(atom_id.index()) {
            match bond(topology, neighbor.bond)?.order() {
                BondOrder::Single | BondOrder::Aromatic => {}
                BondOrder::Double => multiple_count += 1,
                BondOrder::Triple => {
                    if !options.allow_triple_bonds {
                        return Ok(false);
                    }
                    multiple_count += 1;
                }
                _ => {}
            }
            if multiple_count > 1 {
                break;
            }
        }
        if multiple_count > 1 {
            return Ok(false);
        }
    }
    if !options.allow_exocyclic_multiple_bonds {
        for neighbor in topology.adjacency.neighbors_of(atom_id.index()) {
            let incident = bond(topology, neighbor.bond)?;
            if matches!(incident.order(), BondOrder::Double | BondOrder::Triple)
                && rings.num_bond_rings(incident.id()) == 0
            {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn atom_donor_type(
    topology: &TopologyBlock,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    atom_id: AtomId,
    exocyclic_bonds_steal_electrons: bool,
) -> Result<ElectronDonorType, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION getAtomDonorTypeArom
    // RDKit✔️✔️: ElectronDonorType getAtomDonorTypeArom(
    // RDKit✔️✔️:     const Atom *at, bool exocyclicBondsStealElectrons = true) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad atom");
    // RDKit✔️✔️:   const auto &mol = at->getOwningMol();
    // RDKit✔️✔️:   if (at->getAtomicNum() == 0) {
    // RDKit✔️✔️:     return incidentCyclicMultipleBond(at) ? OneElectronDonorType
    // RDKit✔️✔️:                                           : AnyElectronDonorType;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   ElectronDonorType res = NoElectronDonorType;
    // RDKit✔️✔️:   auto nelec = MolOps::countAtomElec(at);
    // RDKit✔️✔️:   int who = -1;
    // RDKit✔️✔️:   if (nelec < 0) {
    // RDKit✔️✔️:     res = NoElectronDonorType;
    // RDKit✔️✔️:   } else if (nelec == 0) {
    // RDKit✔️✔️:     if (incidentNonCyclicMultipleBond(at, who)) {
    // RDKit✔️✔️:       // This is borderline:  no electron to spare but may have an empty
    // RDKit✔️✔️:       // p-orbital
    // RDKit✔️✔️:       // Not sure if this should return vacantElectronDonorType
    // RDKit✔️✔️:       // FIX: explicitly doing this as a note for potential problems
    // RDKit✔️✔️:       //
    // RDKit✔️✔️:       res = VacantElectronDonorType;
    // RDKit✔️✔️:     } else if (incidentCyclicMultipleBond(at)) {
    // RDKit✔️✔️:       // no electron but has one in a in cycle multiple bond
    // RDKit✔️✔️:       res = OneElectronDonorType;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       // no multiple bonds no electrons
    // RDKit✔️✔️:       res = NoElectronDonorType;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (nelec == 1) {
    // RDKit✔️✔️:     if (incidentNonCyclicMultipleBond(at, who)) {
    // RDKit✔️✔️:       // the only available electron is going to be from the
    // RDKit✔️✔️:       // external multiple bond this electron will not be available
    // RDKit✔️✔️:       // for aromaticity if this atom is bonded to a more electro
    // RDKit✔️✔️:       // negative atom
    // RDKit✔️✔️:       const auto at2 = mol.getAtomWithIdx(who);
    // RDKit✔️✔️:       if (exocyclicBondsStealElectrons &&
    // RDKit✔️✔️:           PeriodicTable::getTable()->moreElectroNegative(at2->getAtomicNum(),
    // RDKit✔️✔️:                                                            at->getAtomicNum())) {
    // RDKit✔️✔️:         res = VacantElectronDonorType;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         res = OneElectronDonorType;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       // require that the atom have at least one multiple bond
    // RDKit✔️✔️:       if (incidentMultipleBond(at)) {
    // RDKit✔️✔️:         res = OneElectronDonorType;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // account for the tropylium and cyclopropenyl cation cases
    // RDKit✔️✔️:       else if (at->getFormalCharge() == 1) {
    // RDKit✔️✔️:         res = VacantElectronDonorType;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     if (incidentNonCyclicMultipleBond(at, who)) {
    // RDKit✔️✔️:       // for cases with more than one electron :
    // RDKit✔️✔️:       // if there is an incident multiple bond with an element that
    // RDKit✔️✔️:       // is more electronegative than the this atom, count one less
    // RDKit✔️✔️:       // electron
    // RDKit✔️✔️:       const auto at2 = mol.getAtomWithIdx(who);
    // RDKit✔️✔️:       if (exocyclicBondsStealElectrons &&
    // RDKit✔️✔️:           PeriodicTable::getTable()->moreElectroNegative(at2->getAtomicNum(),
    // RDKit✔️✔️:                                                            at->getAtomicNum())) {
    // RDKit✔️✔️:         --nelec;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (nelec % 2 == 1) {
    // RDKit✔️✔️:       res = OneElectronDonorType;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       res = TwoElectronDonorType;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getAtomDonorTypeArom
    let current_atom = atom(topology, atom_id)?;
    if current_atom.atomic_number() == 0 {
        return Ok(
            if incident_cyclic_multiple_bond(topology, rings, atom_id)? {
                ElectronDonorType::One
            } else {
                ElectronDonorType::Any
            },
        );
    }
    let mut electrons = count_atom_electrons(topology, valence, atom_id)?;
    let external = incident_noncyclic_multiple_bond(topology, rings, atom_id)?;
    let donor = if electrons < 0 {
        ElectronDonorType::None
    } else if electrons == 0 {
        if external.is_some() {
            ElectronDonorType::Vacant
        } else if incident_cyclic_multiple_bond(topology, rings, atom_id)? {
            ElectronDonorType::One
        } else {
            ElectronDonorType::None
        }
    } else if electrons == 1 {
        if let Some(other) = external {
            if exocyclic_bonds_steal_electrons
                && periodic_table_more_electronegative(
                    atom(topology, other)?.atomic_number(),
                    current_atom.atomic_number(),
                )?
            {
                ElectronDonorType::Vacant
            } else {
                ElectronDonorType::One
            }
        } else if incident_multiple_bond(topology, valence, atom_id)? {
            ElectronDonorType::One
        } else if current_atom.formal_charge() == 1 {
            ElectronDonorType::Vacant
        } else {
            ElectronDonorType::None
        }
    } else {
        if let Some(other) = external
            && exocyclic_bonds_steal_electrons
            && periodic_table_more_electronegative(
                atom(topology, other)?.atomic_number(),
                current_atom.atomic_number(),
            )?
        {
            electrons -= 1;
        }
        if electrons % 2 == 1 {
            ElectronDonorType::One
        } else {
            ElectronDonorType::Two
        }
    };
    Ok(donor)
}

type RingNeighborMap = Vec<Vec<usize>>;

fn pick_fused_rings(
    current: usize,
    neighbor_map: &[Vec<usize>],
    done: &mut [bool],
) -> Result<Vec<usize>, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION RingUtils::pickFusedRings
    // RDKit✔️✔️: void pickFusedRings(int curr, const INT_INT_VECT_MAP &neighMap, INT_VECT &res,
    // RDKit✔️✔️:                     boost::dynamic_bitset<> &done, int depth) {
    // RDKit✔️✔️:   auto pos = neighMap.find(curr);
    // RDKit✔️✔️:   PRECONDITION(pos != neighMap.end(), "bad argument");
    // RDKit✔️✔️:   done[curr] = 1;
    // RDKit✔️✔️:   res.push_back(curr);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const auto &neighs = pos->second;
    // RDKit✔️✔️:   for (int neigh : neighs) {
    // RDKit✔️✔️:     if (!done[neigh]) {
    // RDKit✔️✔️:       pickFusedRings(neigh, neighMap, res, done, depth + 1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RingUtils::pickFusedRings
    if current >= neighbor_map.len() || current >= done.len() {
        return Err(AromaticityError::UnsupportedState {
            field: "fused ring index",
            detail: "ring neighbor map index is out of range",
        });
    }
    done[current] = true;
    let mut result = vec![current];
    let mut frames = vec![(current, 0usize)];
    while let Some((ring, next_neighbor)) = frames.last_mut() {
        let Some(&neighbor) = neighbor_map[*ring].get(*next_neighbor) else {
            frames.pop();
            continue;
        };
        *next_neighbor += 1;
        if neighbor >= done.len() {
            return Err(AromaticityError::UnsupportedState {
                field: "fused ring neighbor",
                detail: "ring neighbor index is out of range",
            });
        }
        if !done[neighbor] {
            done[neighbor] = true;
            result.push(neighbor);
            frames.push((neighbor, 0));
        }
    }
    Ok(result)
}

fn check_fused(
    ring_ids: &[usize],
    ring_neighbors: &[Vec<usize>],
) -> Result<bool, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION RingUtils::checkFused
    // RDKit✔️✔️: bool checkFused(const INT_VECT &rids, INT_INT_VECT_MAP &ringNeighs) {
    // RDKit✔️✔️:   auto nrings = rdcast<int>(ringNeighs.size());
    // RDKit✔️✔️:   boost::dynamic_bitset<> done(nrings);
    // RDKit✔️✔️:   int rid;
    // RDKit✔️✔️:   INT_VECT fused;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // mark all rings in the system other than those in rids as done
    // RDKit✔️✔️:   for (const auto &nci : ringNeighs) {
    // RDKit✔️✔️:     rid = nci.first;
    // RDKit✔️✔️:     if (std::find(rids.begin(), rids.end(), rid) == rids.end()) {
    // RDKit✔️✔️:       done[rid] = 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // then pick a fused system from the remaining (i.e. rids)
    // RDKit✔️✔️:   // If the rings in rids are fused we should get back all of them
    // RDKit✔️✔️:   // in fused
    // RDKit✔️✔️:   // if we get a smaller number in fused then rids are not fused
    // RDKit✔️✔️:   pickFusedRings(rids.front(), ringNeighs, fused, done);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   CHECK_INVARIANT(fused.size() <= rids.size(), "");
    // RDKit✔️✔️:   return (fused.size() == rids.size());
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RingUtils::checkFused
    let Some(&first) = ring_ids.first() else {
        return Ok(false);
    };
    let mut selected = vec![false; ring_neighbors.len()];
    for &ring in ring_ids {
        let Some(slot) = selected.get_mut(ring) else {
            return Err(AromaticityError::UnsupportedState {
                field: "fused ring subset",
                detail: "selected ring index is out of range",
            });
        };
        *slot = true;
    }
    let mut done: Vec<bool> = selected.iter().map(|is_selected| !is_selected).collect();
    let fused = pick_fused_rings(first, ring_neighbors, &mut done)?;
    Ok(fused.len() == ring_ids.len())
}

fn make_ring_neighbor_map(
    bond_rings: &[Vec<BondId>],
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
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (i = 0; i < nrings; ++i) {
    // RDKit✔️✔️:     // create an empty INT_VECT at neighMap[i] if it does not yet exist
    // RDKit✔️✔️:     neighMap[i];
    // RDKit✔️✔️:     if (maxSize && brings[i].size() > maxSize) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ring1 = brings[i];
    // RDKit✔️✔️:     for (j = i + 1; j < nrings; ++j) {
    // RDKit✔️✔️:       if (maxSize && brings[j].size() > maxSize) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       INT_VECT inter;
    // RDKit✔️✔️:       Intersect(ring1, brings[j], inter);
    // RDKit✔️✔️:       if (inter.size() > 0 &&
    // RDKit✔️✔️:           (!maxOverlapSize || inter.size() <= maxOverlapSize)) {
    // RDKit✔️✔️:         neighMap[i].push_back(j);
    // RDKit✔️✔️:         neighMap[j].push_back(i);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RingUtils::makeRingNeighborMap
    let mut neighbors = vec![Vec::new(); bond_rings.len()];
    for first in 0..bond_rings.len() {
        if max_size != 0 && bond_rings[first].len() > max_size {
            continue;
        }
        for second in first + 1..bond_rings.len() {
            if max_size != 0 && bond_rings[second].len() > max_size {
                continue;
            }
            let overlap = bond_rings[first]
                .iter()
                .filter(|bond| bond_rings[second].contains(bond))
                .count();
            if overlap != 0 && (max_overlap_size == 0 || overlap <= max_overlap_size) {
                neighbors[first].push(second);
                neighbors[second].push(first);
            }
        }
    }
    neighbors
}

fn convert_to_bonds(
    topology: &TopologyBlock,
    atom_rings: &[Vec<AtomId>],
) -> Result<Vec<Vec<BondId>>, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION RingUtils::convertToBonds (single and vector overloads)
    // RDKit✔️✔️: void convertToBonds(const INT_VECT &ring, INT_VECT &bondRing,
    // RDKit✔️✔️:                     const ROMol &mol) {
    // RDKit✔️✔️:   const auto rsiz = rdcast<unsigned int>(ring.size());
    // RDKit✔️✔️:   bondRing.resize(rsiz);
    // RDKit✔️✔️:   for (unsigned int i = 0; i < (rsiz - 1); i++) {
    // RDKit✔️✔️:     const Bond *bnd = mol.getBondBetweenAtoms(ring[i], ring[i + 1]);
    // RDKit✔️✔️:     if (!bnd) {
    // RDKit✔️✔️:       throw ValueErrorException("expected bond not found");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     bondRing[i] = bnd->getIdx();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // bond from last to first atom
    // RDKit✔️✔️:   const Bond *bnd = mol.getBondBetweenAtoms(ring[rsiz - 1], ring[0]);
    // RDKit✔️✔️:   if (!bnd) {
    // RDKit✔️✔️:     throw ValueErrorException("expected bond not found");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bondRing[rsiz - 1] = bnd->getIdx();
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: void convertToBonds(const VECT_INT_VECT &res, VECT_INT_VECT &brings,
    // RDKit✔️✔️:                     const ROMol &mol) {
    // RDKit✔️✔️:   brings.reserve(res.size());
    // RDKit✔️✔️:   for (const auto &ring : res) {
    // RDKit✔️✔️:     INT_VECT bring;
    // RDKit✔️✔️:     convertToBonds(ring, bring, mol);
    // RDKit✔️✔️:     brings.push_back(bring);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RingUtils::convertToBonds (single and vector overloads)
    let mut bond_rings = Vec::with_capacity(atom_rings.len());
    for ring in atom_rings {
        if ring.is_empty() {
            bond_rings.push(Vec::new());
            continue;
        }
        let mut bond_ring = Vec::with_capacity(ring.len());
        for index in 0..ring.len() {
            let begin = ring[index];
            let end = ring[(index + 1) % ring.len()];
            let Some(neighbor) = topology
                .adjacency
                .neighbors_of(begin.index())
                .iter()
                .find(|neighbor| neighbor.atom_index == end.index())
            else {
                return Err(AromaticityError::ExpectedRingBondNotFound { begin, end });
            };
            bond_ring.push(neighbor.bond);
        }
        bond_rings.push(bond_ring);
    }
    Ok(bond_rings)
}

fn next_combination(combination: &mut [usize], total: usize) -> Option<usize> {
    // BEGIN RDKIT CPP FUNCTION nextCombination
    // RDKit✔️✔️: int nextCombination(INT_VECT &comb, int tot) {
    // RDKit✔️✔️:   int nelem = static_cast<int>(comb.size());
    // RDKit✔️✔️:   int celem = nelem - 1;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   while (comb[celem] == (tot - nelem + celem)) {
    // RDKit✔️✔️:     celem--;
    // RDKit✔️✔️:     if (celem < 0) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int i;
    // RDKit✔️✔️:   comb[celem] += 1;
    // RDKit✔️✔️:   for (i = celem + 1; i < comb.size(); i++) {
    // RDKit✔️✔️:     comb[i] = comb[i - 1] + 1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return celem;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION nextCombination
    let mut current = combination.len().checked_sub(1)?;
    loop {
        let terminal = total.checked_sub(combination.len())?.checked_add(current)?;
        if combination[current] != terminal {
            break;
        }
        current = current.checked_sub(1)?;
    }
    combination[current] = combination[current].checked_add(1)?;
    for index in current + 1..combination.len() {
        combination[index] = combination[index - 1].checked_add(1)?;
    }
    Some(current)
}

fn get_min_max_atom_electrons(donor: ElectronDonorType) -> (i32, i32) {
    // BEGIN RDKIT CPP FUNCTION getMinMaxAtomElecs
    // RDKit✔️✔️: void getMinMaxAtomElecs(ElectronDonorType dtype, int &atlw, int &atup) {
    // RDKit✔️✔️:   switch (dtype) {
    // RDKit✔️✔️:     case AnyElectronDonorType:
    // RDKit✔️✔️:       atlw = 1;
    // RDKit✔️✔️:       atup = 2;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case OneOrTwoElectronDonorType:
    // RDKit✔️✔️:       atlw = 1;
    // RDKit✔️✔️:       atup = 2;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case OneElectronDonorType:
    // RDKit✔️✔️:       atlw = atup = 1;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case TwoElectronDonorType:
    // RDKit✔️✔️:       atlw = atup = 2;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case NoElectronDonorType:
    // RDKit✔️✔️:     case VacantElectronDonorType:
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       atlw = atup = 0;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getMinMaxAtomElecs
    match donor {
        ElectronDonorType::Any | ElectronDonorType::OneOrTwo => (1, 2),
        ElectronDonorType::One => (1, 1),
        ElectronDonorType::Two => (2, 2),
        ElectronDonorType::None | ElectronDonorType::Vacant => (0, 0),
    }
}

fn apply_huckel(
    ring: &[AtomId],
    donors: &[ElectronDonorType],
    min_ring_size: usize,
) -> Result<bool, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION applyHuckel
    // RDKit✔️✔️: bool applyHuckel(ROMol &, const INT_VECT &ring, const VECT_EDON_TYPE &edon,
    // RDKit✔️✔️:                  unsigned int minRingSize) {
    // RDKit✔️✔️:   if (ring.size() < minRingSize) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   int atlw, atup, rlw, rup, rie;
    // RDKit✔️✔️:   bool aromatic = false;
    // RDKit✔️✔️:   rlw = 0;
    // RDKit✔️✔️:   rup = 0;
    // RDKit✔️✔️:   unsigned int nAnyElectronDonorType = 0;
    // RDKit✔️✔️:   for (auto idx : ring) {
    // RDKit✔️✔️:     ElectronDonorType edonType = edon[idx];
    // RDKit✔️✔️:     if (edonType == AnyElectronDonorType) {
    // RDKit✔️✔️:       ++nAnyElectronDonorType;
    // RDKit✔️✔️:       if (nAnyElectronDonorType > 1) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     getMinMaxAtomElecs(edonType, atlw, atup);
    // RDKit✔️✔️:     rlw += atlw;
    // RDKit✔️✔️:     rup += atup;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (rup >= 6) {
    // RDKit✔️✔️:     for (rie = rlw; rie <= rup; ++rie) {
    // RDKit✔️✔️:       if ((rie - 2) % 4 == 0) {
    // RDKit✔️✔️:         aromatic = true;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (rup == 2) {
    // RDKit✔️✔️:     aromatic = true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return aromatic;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION applyHuckel
    if ring.len() < min_ring_size {
        return Ok(false);
    }
    let mut lower = 0i32;
    let mut upper = 0i32;
    let mut any_count = 0usize;
    for &atom_id in ring {
        let donor = *donors
            .get(atom_id.index())
            .ok_or(AromaticityError::AtomOutOfRange {
                atom: atom_id,
                atom_count: donors.len(),
            })?;
        if donor == ElectronDonorType::Any {
            any_count += 1;
            if any_count > 1 {
                return Ok(false);
            }
        }
        let (atom_lower, atom_upper) = get_min_max_atom_electrons(donor);
        lower = lower
            .checked_add(atom_lower)
            .ok_or(AromaticityError::IntegerOverflow {
                field: "minimum aromatic electrons",
            })?;
        upper = upper
            .checked_add(atom_upper)
            .ok_or(AromaticityError::IntegerOverflow {
                field: "maximum aromatic electrons",
            })?;
    }
    if upper >= 6 {
        Ok((lower..=upper).any(|electrons| (electrons - 2) % 4 == 0))
    } else {
        Ok(upper == 2)
    }
}

fn mark_atoms_bonds_aromatic(
    topology: &mut TopologyBlock,
    bond_rings: &[Vec<BondId>],
    ring_ids: &[usize],
    done_bonds: &mut BTreeSet<usize>,
) -> Result<(), AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION markAtomsBondsArom
    // RDKit✔️✔️: void markAtomsBondsArom(const VECT_INT_VECT &brings, const INT_VECT &ringIds,
    // RDKit✔️✔️:                         std::set<unsigned int> &doneBonds,
    // RDKit✔️✔️:                         const std::vector<Bond *> &bondsByIdx) {
    // RDKit✔️✔️:   // mark the bonds
    // RDKit✔️✔️:   // here we want to be careful. We don't want to mark the fusing bonds
    // RDKit✔️✔️:   // as aromatic - only the outside bonds in a fused system are marked aromatic.
    // RDKit✔️✔️:   // - loop through the rings and count the number of times each bond appears in
    // RDKit✔️✔️:   //   all the fused rings.
    // RDKit✔️✔️:   // - bonds that appears only once are marked aromatic
    // RDKit✔️✔️:   INT_MAP_INT bndCntr;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto ri : ringIds) {
    // RDKit✔️✔️:     const auto &bring = brings[ri];
    // RDKit✔️✔️:     for (auto bi : bring) {
    // RDKit✔️✔️:       ++bndCntr[bi];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // now mark single or double bonds that have a count of 1 and the atoms they
    // RDKit✔️✔️:   // connect as aromatic
    // RDKit✔️✔️:   for (const auto &bci : bndCntr) {
    // RDKit✔️✔️:     if (bci.second == 1) {
    // RDKit✔️✔️:       auto bond = bondsByIdx[bci.first];
    // RDKit✔️✔️:       bond->setIsAromatic(true);
    // RDKit✔️✔️:       switch (bond->getBondType()) {
    // RDKit✔️✔️:         case Bond::SINGLE:
    // RDKit✔️✔️:         case Bond::DOUBLE:
    // RDKit✔️✔️:           bond->setBondType(Bond::AROMATIC);
    // RDKit✔️✔️:           bond->getBeginAtom()->setIsAromatic(true);
    // RDKit✔️✔️:           bond->getEndAtom()->setIsAromatic(true);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       doneBonds.insert(bond->getIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION markAtomsBondsArom
    let mut counts = vec![0usize; topology.bonds.len()];
    for &ring_id in ring_ids {
        let ring = bond_rings
            .get(ring_id)
            .ok_or(AromaticityError::UnsupportedState {
                field: "aromatic ring id",
                detail: "ring index is out of range",
            })?;
        for &bond_id in ring {
            let count =
                counts
                    .get_mut(bond_id.index())
                    .ok_or(AromaticityError::BondOutOfRange {
                        bond: bond_id,
                        bond_count: topology.bonds.len(),
                    })?;
            *count = count
                .checked_add(1)
                .ok_or(AromaticityError::IntegerOverflow {
                    field: "fused bond occurrence count",
                })?;
        }
    }
    for (bond_index, count) in counts.into_iter().enumerate() {
        if count != 1 {
            continue;
        }
        let bond_id = BondId::new(bond_index);
        let (begin, end, order) = {
            let current = bond(topology, bond_id)?;
            (current.begin(), current.end(), current.order())
        };
        let current =
            topology
                .bonds
                .get_mut(bond_index)
                .ok_or(AromaticityError::BondOutOfRange {
                    bond: bond_id,
                    bond_count: bond_index,
                })?;
        current.set_aromatic(true);
        if matches!(order, BondOrder::Single | BondOrder::Double) {
            current.set_order(BondOrder::Aromatic);
            let atom_count = topology.atoms.len();
            topology
                .atoms
                .get_mut(begin.index())
                .ok_or(AromaticityError::AtomOutOfRange {
                    atom: begin,
                    atom_count,
                })?
                .set_aromatic(true);
            topology
                .atoms
                .get_mut(end.index())
                .ok_or(AromaticityError::AtomOutOfRange {
                    atom: end,
                    atom_count,
                })?
                .set_aromatic(true);
        }
        done_bonds.insert(bond_index);
    }
    Ok(())
}

fn apply_huckel_to_fused(
    topology: &mut TopologyBlock,
    atom_rings: &[Vec<AtomId>],
    bond_rings: &[Vec<BondId>],
    fused: &[usize],
    donors: &[ElectronDonorType],
    ring_neighbors: &[Vec<usize>],
    max_fused_rings: usize,
    min_ring_size: usize,
) -> Result<usize, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION applyHuckelToFused
    // RDKit✔️✔️: void applyHuckelToFused(
    // RDKit✔️✔️:     ROMol &mol,
    // RDKit✔️✔️:     const VECT_INT_VECT &srings,
    // RDKit✔️✔️:     const VECT_INT_VECT &brings,
    // RDKit✔️✔️:     const INT_VECT &fused,
    // RDKit✔️✔️:     const VECT_EDON_TYPE &edon,
    // RDKit✔️✔️:     INT_INT_VECT_MAP &ringNeighs,
    // RDKit✔️✔️:     int &narom,
    // RDKit✔️✔️:     unsigned int maxNumFusedRings, const std::vector<Bond *> &bondsByIdx,
    // RDKit✔️✔️:     unsigned int minRingSize) {
    // RDKit✔️✔️:   std::unordered_set<int> aromRings;
    // RDKit✔️✔️:   auto nrings = rdcast<unsigned int>(fused.size());
    // RDKit✔️✔️:   INT_VECT curRs;
    // RDKit✔️✔️:   curRs.push_back(fused.front());
    // RDKit✔️✔️:   int pos = -1;
    // RDKit✔️✔️:   unsigned int i;
    // RDKit✔️✔️:   unsigned int curSize = 0;
    // RDKit✔️✔️:   INT_VECT comb;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   size_t nRingBonds;
    // RDKit✔️✔️:   {
    // RDKit✔️✔️:     boost::dynamic_bitset<> fusedBonds(mol.getNumBonds());
    // RDKit✔️✔️:     for (auto ridx : fused) {
    // RDKit✔️✔️:       for (auto bidx : brings[ridx]) {
    // RDKit✔️✔️:         fusedBonds[bidx] = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     nRingBonds = rdcast<unsigned int>(fusedBonds.count());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::set<unsigned int> doneBonds;
    // RDKit✔️✔️:   while (1) {
    // RDKit✔️✔️:     if (pos == -1) {
    // RDKit✔️✔️:       if ((curSize == 2) && (nrings > 300)) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "Aromaticity detection halted on some rings due to ring system size."
    // RDKit✔️✔️:             << std::endl;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       ++curSize;
    // RDKit✔️✔️:       if (curSize > std::min(nrings, maxNumFusedRings) ||
    // RDKit✔️✔️:           doneBonds.size() >= nRingBonds) {
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       comb.resize(curSize);
    // RDKit✔️✔️:       std::iota(comb.begin(), comb.end(), 0);
    // RDKit✔️✔️:       pos = 0;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       pos = nextCombination(comb, nrings);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (pos == -1) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     curRs.clear();
    // RDKit✔️✔️:     std::transform(comb.begin(), comb.end(), std::back_inserter(curRs),
    // RDKit✔️✔️:                    [&fused](const int i) { return fused[i]; });
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (ringNeighs.size() && !RingUtils::checkFused(curRs, ringNeighs)) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     INT_VECT atsInRingSystem(mol.getNumAtoms(), 0);
    // RDKit✔️✔️:     for (auto ridx : curRs) {
    // RDKit✔️✔️:       for (auto rid : srings[ridx]) {
    // RDKit✔️✔️:         ++atsInRingSystem[rid];
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     INT_VECT unon;
    // RDKit✔️✔️:     for (i = 0; i < atsInRingSystem.size(); ++i) {
    // RDKit✔️✔️:       if (atsInRingSystem[i] == 1 || atsInRingSystem[i] == 2) {
    // RDKit✔️✔️:         unon.push_back(i);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (applyHuckel(mol, unon, edon, minRingSize)) {
    // RDKit✔️✔️:       markAtomsBondsArom(brings, curRs, doneBonds, bondsByIdx);
    // RDKit✔️✔️:       std::copy(curRs.begin(), curRs.end(),
    // RDKit✔️✔️:                 std::inserter(aromRings, aromRings.begin()));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   narom += rdcast<int>(aromRings.size());
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION applyHuckelToFused
    if fused.is_empty() {
        return Ok(0);
    }
    let mut fused_bonds = BTreeSet::new();
    for &ring_id in fused {
        let ring = bond_rings
            .get(ring_id)
            .ok_or(AromaticityError::UnsupportedState {
                field: "fused ring id",
                detail: "ring index is out of range",
            })?;
        fused_bonds.extend(ring.iter().map(|bond| bond.index()));
    }
    let mut aromatic_rings = BTreeSet::new();
    let mut done_bonds = BTreeSet::new();
    let maximum = fused.len().min(max_fused_rings);
    for current_size in 1..=maximum {
        if current_size == 3 && fused.len() > 300 {
            break;
        }
        if done_bonds.len() >= fused_bonds.len() {
            break;
        }
        let mut combination: Vec<usize> = (0..current_size).collect();
        loop {
            let current_rings: Vec<usize> = combination.iter().map(|&index| fused[index]).collect();
            if ring_neighbors.is_empty() || check_fused(&current_rings, ring_neighbors)? {
                let mut occurrences = vec![0u8; topology.atoms.len()];
                for &ring_id in &current_rings {
                    let ring =
                        atom_rings
                            .get(ring_id)
                            .ok_or(AromaticityError::UnsupportedState {
                                field: "fused atom ring id",
                                detail: "ring index is out of range",
                            })?;
                    for &atom_id in ring {
                        let count = occurrences.get_mut(atom_id.index()).ok_or(
                            AromaticityError::AtomOutOfRange {
                                atom: atom_id,
                                atom_count: topology.atoms.len(),
                            },
                        )?;
                        *count = count
                            .checked_add(1)
                            .ok_or(AromaticityError::IntegerOverflow {
                                field: "fused atom occurrence count",
                            })?;
                    }
                }
                let perimeter: Vec<AtomId> = occurrences
                    .into_iter()
                    .enumerate()
                    .filter_map(|(index, count)| {
                        matches!(count, 1 | 2).then_some(AtomId::new(index))
                    })
                    .collect();
                if apply_huckel(&perimeter, donors, min_ring_size)? {
                    mark_atoms_bonds_aromatic(
                        topology,
                        bond_rings,
                        &current_rings,
                        &mut done_bonds,
                    )?;
                    aromatic_rings.extend(current_rings);
                }
            }
            if next_combination(&mut combination, fused.len()).is_none() {
                break;
            }
        }
    }
    Ok(aromatic_rings.len())
}

fn aromaticity_helper(
    topology: &TopologyBlock,
    rings: &RingInfo,
    min_ring_size: usize,
    max_ring_size: usize,
    include_fused: bool,
) -> Result<AromaticityAssignment, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION aromaticityHelper
    // RDKit✔️✔️: int aromaticityHelper(RWMol &mol, const VECT_INT_VECT &srings,
    // RDKit✔️✔️:                       unsigned int minRingSize, unsigned int maxRingSize,
    // RDKit✔️✔️:                       bool includeFused) {
    // RDKit✔️✔️:   int narom = 0;
    // RDKit✔️✔️:   int natoms = mol.getNumAtoms();
    // RDKit✔️✔️:   boost::dynamic_bitset<> acands(natoms);
    // RDKit✔️✔️:   boost::dynamic_bitset<> aseen(natoms);
    // RDKit✔️✔️:   VECT_EDON_TYPE edon(natoms);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   VECT_INT_VECT cRings;
    // RDKit✔️✔️:   for (auto &sring : srings) {
    // RDKit✔️✔️:     size_t ringSz = sring.size();
    // RDKit✔️✔️:     if ((minRingSize && ringSz < minRingSize) ||
    // RDKit✔️✔️:         (maxRingSize && ringSz > maxRingSize)) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     bool allAromatic = true;
    // RDKit✔️✔️:     bool allDummy = true;
    // RDKit✔️✔️:     for (auto firstIdx : sring) {
    // RDKit✔️✔️:       const auto at = mol.getAtomWithIdx(firstIdx);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (allDummy && !isAtomDummy(at)) {
    // RDKit✔️✔️:         allDummy = false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (aseen[firstIdx]) {
    // RDKit✔️✔️:         if (!acands[firstIdx]) {
    // RDKit✔️✔️:           allAromatic = false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       aseen[firstIdx] = 1;
    // RDKit✔️✔️:
    // RDKit✔️✔️:       edon[firstIdx] = getAtomDonorTypeArom(at);
    // RDKit✔️✔️:       acands[firstIdx] = isAtomCandForArom(at, edon[firstIdx]);
    // RDKit✔️✔️:       if (!acands[firstIdx]) {
    // RDKit✔️✔️:         allAromatic = false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (allAromatic && !allDummy) {
    // RDKit✔️✔️:       cRings.push_back(sring);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   VECT_INT_VECT brings;
    // RDKit✔️✔️:   RingUtils::convertToBonds(cRings, brings, mol);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<Bond *> bondsByIdx;
    // RDKit✔️✔️:   bondsByIdx.reserve(mol.getNumBonds());
    // RDKit✔️✔️:   for (auto b : mol.bonds()) {
    // RDKit✔️✔️:     bondsByIdx.push_back(b);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!includeFused) {
    // RDKit✔️✔️:     INT_INT_VECT_MAP neighMap;
    // RDKit✔️✔️:     for (size_t ri = 0; ri < cRings.size(); ++ri) {
    // RDKit✔️✔️:       INT_VECT fused;
    // RDKit✔️✔️:       fused.push_back(ri);
    // RDKit✔️✔️:       const unsigned int maxFused = 6;
    // RDKit✔️✔️:       const unsigned int minRingSize = 0;
    // RDKit✔️✔️:       applyHuckelToFused(mol, cRings, brings, fused, edon, neighMap, narom,
    // RDKit✔️✔️:                          maxFused, bondsByIdx, minRingSize);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     INT_INT_VECT_MAP neighMap;
    // RDKit✔️✔️:     RingUtils::makeRingNeighborMap(brings, neighMap, maxFusedAromaticRingSize,
    // RDKit✔️✔️:                                    1);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     INT_VECT doneRs;
    // RDKit✔️✔️:     int curr = 0;
    // RDKit✔️✔️:     auto cnrs = rdcast<int>(cRings.size());
    // RDKit✔️✔️:     boost::dynamic_bitset<> fusDone(cnrs);
    // RDKit✔️✔️:     INT_VECT fused;
    // RDKit✔️✔️:     while (curr < cnrs) {
    // RDKit✔️✔️:       fused.clear();
    // RDKit✔️✔️:       RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
    // RDKit✔️✔️:       applyHuckelToFused(mol, cRings, brings, fused, edon, neighMap, narom, 6,
    // RDKit✔️✔️:                          bondsByIdx);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       int rix;
    // RDKit✔️✔️:       for (rix = 0; rix < cnrs; ++rix) {
    // RDKit✔️✔️:         if (!fusDone[rix]) {
    // RDKit✔️✔️:           curr = rix;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (rix == cnrs) {
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   mol.setProp(common_properties::numArom, narom, true);
    // RDKit✔️✔️:   return narom;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION aromaticityHelper
    validate_inputs(topology, rings)?;
    let valence = crate::assign_valence_for_topology(topology, crate::ValenceModel::RdkitLike)?;
    validate_valence(topology, &valence)?;
    let mut seen = vec![false; topology.atoms.len()];
    let mut candidates = vec![false; topology.atoms.len()];
    let mut donors = vec![ElectronDonorType::None; topology.atoms.len()];
    let mut candidate_rings = Vec::new();
    for source_ring in rings.atom_rings() {
        let ring_size = source_ring.len();
        if (min_ring_size != 0 && ring_size < min_ring_size)
            || (max_ring_size != 0 && ring_size > max_ring_size)
        {
            continue;
        }
        let mut all_aromatic = true;
        let mut all_dummy = true;
        for &atom_id in source_ring {
            let current_atom = atom(topology, atom_id)?;
            if all_dummy && current_atom.atomic_number() != 0 {
                all_dummy = false;
            }
            if seen[atom_id.index()] {
                if !candidates[atom_id.index()] {
                    all_aromatic = false;
                }
                continue;
            }
            seen[atom_id.index()] = true;
            donors[atom_id.index()] = atom_donor_type(topology, rings, &valence, atom_id, true)?;
            candidates[atom_id.index()] = is_atom_candidate(
                topology,
                rings,
                &valence,
                atom_id,
                donors[atom_id.index()],
                CandidateOptions::default(),
            )?;
            if !candidates[atom_id.index()] {
                all_aromatic = false;
            }
        }
        if all_aromatic && !all_dummy {
            candidate_rings.push(source_ring.clone());
        }
    }
    let bond_rings = convert_to_bonds(topology, &candidate_rings)?;
    let mut working = topology.clone();
    let mut aromatic_ring_count = 0usize;
    if include_fused {
        let neighbor_map = make_ring_neighbor_map(&bond_rings, 24, 1);
        let mut done = vec![false; candidate_rings.len()];
        while let Some(current) = done.iter().position(|is_done| !*is_done) {
            let fused = pick_fused_rings(current, &neighbor_map, &mut done)?;
            aromatic_ring_count = aromatic_ring_count
                .checked_add(apply_huckel_to_fused(
                    &mut working,
                    &candidate_rings,
                    &bond_rings,
                    &fused,
                    &donors,
                    &neighbor_map,
                    6,
                    0,
                )?)
                .ok_or(AromaticityError::IntegerOverflow {
                    field: "aromatic ring count",
                })?;
        }
    } else {
        for ring_index in 0..candidate_rings.len() {
            aromatic_ring_count = aromatic_ring_count
                .checked_add(apply_huckel_to_fused(
                    &mut working,
                    &candidate_rings,
                    &bond_rings,
                    &[ring_index],
                    &donors,
                    &[],
                    6,
                    0,
                )?)
                .ok_or(AromaticityError::IntegerOverflow {
                    field: "aromatic ring count",
                })?;
        }
    }
    working.validate()?;
    Ok(AromaticityAssignment {
        topology: working,
        aromatic_ring_count,
    })
}

fn mdl_aromaticity_helper(
    topology: &TopologyBlock,
    rings: &RingInfo,
) -> Result<AromaticityAssignment, AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION mdlAromaticityHelper
    // RDKit✔️✔️: int mdlAromaticityHelper(RWMol &mol, const VECT_INT_VECT &srings) {
    // RDKit✔️✔️:   int narom = 0;
    // RDKit✔️✔️:   int natoms = mol.getNumAtoms();
    // RDKit✔️✔️:   boost::dynamic_bitset<> acands(natoms);
    // RDKit✔️✔️:   boost::dynamic_bitset<> aseen(natoms);
    // RDKit✔️✔️:   VECT_EDON_TYPE edon(natoms);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   VECT_INT_VECT cRings;
    // RDKit✔️✔️:   for (auto &sring : srings) {
    // RDKit✔️✔️:     bool allAromatic = true;
    // RDKit✔️✔️:     bool allDummy = true;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     for (auto firstIdx : sring) {
    // RDKit✔️✔️:       const auto at = mol.getAtomWithIdx(firstIdx);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (allDummy && at->getAtomicNum() != 0) {
    // RDKit✔️✔️:         allDummy = false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (aseen[firstIdx]) {
    // RDKit✔️✔️:         if (!acands[firstIdx]) {
    // RDKit✔️✔️:           allAromatic = false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       aseen[firstIdx] = 1;
    // RDKit✔️✔️:
    // RDKit✔️✔️:       edon[firstIdx] = getAtomDonorTypeArom(at, false);
    // RDKit✔️✔️:       if (edon[firstIdx] != OneElectronDonorType) {
    // RDKit✔️✔️:         allAromatic = false;
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       const bool allowThirdRow = false;
    // RDKit✔️✔️:       const bool allowTripleBonds = false;
    // RDKit✔️✔️:       const bool allowHigherExceptions = false;
    // RDKit✔️✔️:       const bool onlyCorN = true;
    // RDKit✔️✔️:       const bool allowExocyclicMultipleBonds = false;
    // RDKit✔️✔️:       acands[firstIdx] = isAtomCandForArom(
    // RDKit✔️✔️:           at, edon[firstIdx], allowThirdRow, allowTripleBonds,
    // RDKit✔️✔️:           allowHigherExceptions, onlyCorN, allowExocyclicMultipleBonds);
    // RDKit✔️✔️:       if (!acands[firstIdx]) {
    // RDKit✔️✔️:         allAromatic = false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (allAromatic && !allDummy) {
    // RDKit✔️✔️:       cRings.push_back(sring);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   VECT_INT_VECT brings;
    // RDKit✔️✔️:   RingUtils::convertToBonds(cRings, brings, mol);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   INT_INT_VECT_MAP neighMap;
    // RDKit✔️✔️:   RingUtils::makeRingNeighborMap(brings, neighMap, maxFusedAromaticRingSize, 1);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   INT_VECT doneRs;
    // RDKit✔️✔️:   int curr = 0;
    // RDKit✔️✔️:   auto cnrs = rdcast<int>(cRings.size());
    // RDKit✔️✔️:   boost::dynamic_bitset<> fusDone(cnrs);
    // RDKit✔️✔️:   INT_VECT fused;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<Bond *> bondsByIdx;
    // RDKit✔️✔️:   bondsByIdx.reserve(mol.getNumBonds());
    // RDKit✔️✔️:   for (const auto b : mol.bonds()) {
    // RDKit✔️✔️:     bondsByIdx.push_back(b);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   while (curr < cnrs) {
    // RDKit✔️✔️:     fused.clear();
    // RDKit✔️✔️:     RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
    // RDKit✔️✔️:     const unsigned int maxFused = 6;
    // RDKit✔️✔️:     const unsigned int minRingSize = 6;
    // RDKit✔️✔️:     applyHuckelToFused(mol, cRings, brings, fused, edon, neighMap, narom,
    // RDKit✔️✔️:                        maxFused, bondsByIdx, minRingSize);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     int rix;
    // RDKit✔️✔️:     for (rix = 0; rix < cnrs; ++rix) {
    // RDKit✔️✔️:       if (!fusDone[rix]) {
    // RDKit✔️✔️:         curr = rix;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (rix == cnrs) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   mol.setProp(common_properties::numArom, narom, true);
    // RDKit✔️✔️:   return narom;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION mdlAromaticityHelper
    validate_inputs(topology, rings)?;
    let valence = crate::assign_valence_for_topology(topology, crate::ValenceModel::RdkitLike)?;
    validate_valence(topology, &valence)?;
    let mut seen = vec![false; topology.atoms.len()];
    let mut candidates = vec![false; topology.atoms.len()];
    let mut donors = vec![ElectronDonorType::None; topology.atoms.len()];
    let mut candidate_rings = Vec::new();
    let options = CandidateOptions {
        allow_third_row: false,
        allow_triple_bonds: false,
        allow_higher_exceptions: false,
        only_c_or_n: true,
        allow_exocyclic_multiple_bonds: false,
    };
    for source_ring in rings.atom_rings() {
        let mut all_aromatic = true;
        let mut all_dummy = true;
        for &atom_id in source_ring {
            let current_atom = atom(topology, atom_id)?;
            if all_dummy && current_atom.atomic_number() != 0 {
                all_dummy = false;
            }
            if seen[atom_id.index()] {
                if !candidates[atom_id.index()] {
                    all_aromatic = false;
                }
                continue;
            }
            seen[atom_id.index()] = true;
            donors[atom_id.index()] = atom_donor_type(topology, rings, &valence, atom_id, false)?;
            if donors[atom_id.index()] != ElectronDonorType::One {
                all_aromatic = false;
                continue;
            }
            candidates[atom_id.index()] = is_atom_candidate(
                topology,
                rings,
                &valence,
                atom_id,
                donors[atom_id.index()],
                options,
            )?;
            if !candidates[atom_id.index()] {
                all_aromatic = false;
            }
        }
        if all_aromatic && !all_dummy {
            candidate_rings.push(source_ring.clone());
        }
    }
    let bond_rings = convert_to_bonds(topology, &candidate_rings)?;
    let neighbor_map = make_ring_neighbor_map(&bond_rings, 24, 1);
    let mut working = topology.clone();
    let mut aromatic_ring_count = 0usize;
    let mut done = vec![false; candidate_rings.len()];
    while let Some(current) = done.iter().position(|is_done| !*is_done) {
        let fused = pick_fused_rings(current, &neighbor_map, &mut done)?;
        aromatic_ring_count = aromatic_ring_count
            .checked_add(apply_huckel_to_fused(
                &mut working,
                &candidate_rings,
                &bond_rings,
                &fused,
                &donors,
                &neighbor_map,
                6,
                6,
            )?)
            .ok_or(AromaticityError::IntegerOverflow {
                field: "aromatic ring count",
            })?;
    }
    working.validate()?;
    Ok(AromaticityAssignment {
        topology: working,
        aromatic_ring_count,
    })
}

fn set_mmff_aromaticity(
    topology: &mut TopologyBlock,
    rings: &RingInfo,
) -> Result<(), AromaticityError> {
    // BEGIN RDKIT CPP FUNCTION setMMFFAromaticity
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
    // END RDKIT CPP FUNCTION setMMFFAromaticity
    let atom_rings = rings.atom_rings();
    let mut perceived_atoms = vec![false; topology.atoms.len()];
    let mut aromatic_rings = vec![false; atom_rings.len()];
    let mut aromatic_rows_set = 0i32;
    let mut old_aromatic_rows_set = -1i32;
    let mut all_ring_atoms_perceived = false;

    while !all_ring_atoms_perceived
        && !atom_rings.is_empty()
        && aromatic_rows_set > old_aromatic_rows_set
    {
        let valence = crate::assign_valence_for_topology(topology, crate::ValenceModel::RdkitLike)?;
        validate_valence(topology, &valence)?;
        for (ring_index, ring) in atom_rings.iter().enumerate() {
            let mut move_to_next_ring = false;
            let mut has_n_o_or_divalent_s = false;
            let mut exocyclic_double_bond = false;
            let mut pi_electrons = 0usize;

            for (position, &atom_id) in ring.iter().enumerate() {
                if move_to_next_ring {
                    break;
                }
                let current_atom = atom(topology, atom_id)?;
                let atomic_number = current_atom.atomic_number();
                if atomic_number == 7
                    || atomic_number == 8
                    || (atomic_number == 16
                        && topology.adjacency.neighbors_of(atom_id.index()).len() == 2)
                {
                    has_n_o_or_divalent_s = true;
                }
                let next = ring[(position + 1) % ring.len()];
                let ring_bond_id = bond_between_atoms(topology, atom_id, next)?;
                if bond(topology, ring_bond_id)?.order() == BondOrder::Double {
                    pi_electrons =
                        pi_electrons
                            .checked_add(2)
                            .ok_or(AromaticityError::IntegerOverflow {
                                field: "MMFF94 pi electron count",
                            })?;
                    continue;
                }

                let (explicit_valence, implicit_hydrogens) = valence_row(&valence, atom_id)?;
                if atomic_number != 6
                    && !(atomic_number == 7
                        && explicit_valence.checked_add(implicit_hydrogens).ok_or(
                            AromaticityError::IntegerOverflow {
                                field: "MMFF94 nitrogen total valence",
                            },
                        )? == 4)
                {
                    continue;
                }

                for neighbor in topology.adjacency.neighbors_of(atom_id.index()) {
                    let neighbor_id = AtomId::new(neighbor.atom_index);
                    if ring.contains(&neighbor_id) {
                        continue;
                    }
                    let exocyclic_bond = bond(topology, neighbor.bond)?;
                    if exocyclic_bond.order() == BondOrder::Single {
                        continue;
                    }
                    if rings.num_atom_rings(neighbor_id) != 0
                        && !perceived_atoms[neighbor_id.index()]
                    {
                        move_to_next_ring = true;
                        break;
                    }
                    if exocyclic_bond.order() == BondOrder::Double {
                        if atom(topology, neighbor_id)?.is_aromatic() {
                            pi_electrons = pi_electrons.checked_add(1).ok_or(
                                AromaticityError::IntegerOverflow {
                                    field: "MMFF94 pi electron count",
                                },
                            )?;
                        } else {
                            exocyclic_double_bond = true;
                        }
                    }
                }
            }
            if move_to_next_ring {
                continue;
            }

            let mut can_be_aromatic = true;
            for &atom_id in ring {
                perceived_atoms[atom_id.index()] = true;
                let current_atom = atom(topology, atom_id)?;
                if matches!(current_atom.atomic_number(), 6 | 7)
                    && current_atom.hybridization() != Hybridization::Sp2
                {
                    can_be_aromatic = false;
                }
            }
            if !can_be_aromatic {
                continue;
            }
            if has_n_o_or_divalent_s && !exocyclic_double_bond && ring.len() % 2 == 1 {
                pi_electrons =
                    pi_electrons
                        .checked_add(2)
                        .ok_or(AromaticityError::IntegerOverflow {
                            field: "MMFF94 pi electron count",
                        })?;
            }
            if pi_electrons > 2 && (pi_electrons - 2) % 4 == 0 {
                aromatic_rings[ring_index] = true;
                for &atom_id in ring {
                    topology.atoms[atom_id.index()].set_aromatic(true);
                }
            }
        }

        old_aromatic_rows_set = aromatic_rows_set;
        aromatic_rows_set = 0;
        all_ring_atoms_perceived = true;
        for ring in atom_rings {
            for &atom_id in ring {
                if perceived_atoms[atom_id.index()] {
                    aromatic_rows_set = aromatic_rows_set.checked_add(1).ok_or(
                        AromaticityError::IntegerOverflow {
                            field: "MMFF94 perceived aromatic atom row count",
                        },
                    )?;
                } else {
                    all_ring_atoms_perceived = false;
                }
            }
        }
    }

    for (ring_index, ring) in atom_rings.iter().enumerate() {
        if !aromatic_rings[ring_index] {
            continue;
        }
        for (position, &atom_id) in ring.iter().enumerate() {
            let next = ring[(position + 1) % ring.len()];
            let bond_id = bond_between_atoms(topology, atom_id, next)?;
            let current_bond = &mut topology.bonds[bond_id.index()];
            current_bond.set_order(BondOrder::Aromatic);
            current_bond.set_aromatic(true);
        }
    }

    for (ring_index, ring) in atom_rings.iter().enumerate() {
        if !aromatic_rings[ring_index] {
            continue;
        }
        for &atom_id in ring {
            if topology.atoms[atom_id.index()].atomic_number() == 6 {
                continue;
            }
            let implicit = crate::calculate_implicit_valence_for_topology(
                topology, atom_id, -1, false, false,
            )?;
            let _ =
                crate::calculate_explicit_valence_for_topology(topology, atom_id, false, false)?;
            if implicit != 0 {
                let explicit_hydrogens =
                    u8::try_from(implicit).map_err(|_| AromaticityError::InvalidValenceRow {
                        atom: atom_id,
                        field: "MMFF94 implicit valence",
                        value: implicit,
                    })?;
                topology.atoms[atom_id.index()].set_explicit_hydrogens(explicit_hydrogens);
                let _ = crate::calculate_implicit_valence_for_topology(
                    topology, atom_id, -1, false, false,
                )?;
            }
        }
    }
    Ok(())
}

fn bond_between_atoms(
    topology: &TopologyBlock,
    begin: AtomId,
    end: AtomId,
) -> Result<BondId, AromaticityError> {
    topology
        .adjacency
        .neighbors_of(begin.index())
        .iter()
        .find(|neighbor| neighbor.atom_index == end.index())
        .map(|neighbor| neighbor.bond)
        .ok_or(AromaticityError::ExpectedRingBondNotFound { begin, end })
}

fn mmff94_aromaticity_helper(
    topology: &TopologyBlock,
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
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return narom;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION mmff94AromaticityHelper
    validate_inputs(topology, rings)?;
    let mut working = topology.clone();
    if working.atoms.iter().any(Atom::is_aromatic) {
        working = crate::kekulize(&working, &crate::KekulizeParams::default())?.topology;
    }
    // The source's private `_MMFFSanitized` cache flag has no detached-model
    // meaning. Each call owns a fresh working topology, so the source's guarded
    // initialization is completed locally without manufacturing a public prop.
    set_mmff_aromaticity(&mut working, rings)?;

    let mut aromatic_ring_count = 1usize;
    for ring in rings.atom_rings() {
        if ring
            .iter()
            .all(|atom_id| working.atoms[atom_id.index()].is_aromatic())
        {
            aromatic_ring_count =
                aromatic_ring_count
                    .checked_add(1)
                    .ok_or(AromaticityError::IntegerOverflow {
                        field: "MMFF94 aromatic ring count",
                    })?;
        }
    }
    working.validate()?;
    Ok(AromaticityAssignment {
        topology: working,
        aromatic_ring_count,
    })
}

pub fn assign_aromaticity(
    topology: &TopologyBlock,
    rings: &RingInfo,
    params: &AromaticityParams,
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
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int res;
    // RDKit✔️✔️:   switch (model) {
    // RDKit✔️✔️:     case AROMATICITY_DEFAULT:
    // RDKit✔️✔️:     case AROMATICITY_RDKIT:
    // RDKit✔️✔️:       res = aromaticityHelper(mol, srings, 0, 0, true);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case AROMATICITY_SIMPLE:
    // RDKit✔️✔️:       res = aromaticityHelper(mol, srings, 5, 6, false);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case AROMATICITY_MDL:
    // RDKit✔️✔️:       res = mdlAromaticityHelper(mol, srings);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case AROMATICITY_MMFF94:
    // RDKit✔️✔️:       res = mmff94AromaticityHelper(mol, srings);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case AROMATICITY_CUSTOM:
    // RDKit✔️✔️:       PRECONDITION(
    // RDKit✔️✔️:           func,
    // RDKit✔️✔️:           "function must be set when aromaticity model is AROMATICITY_CUSTOM");
    // RDKit✔️✔️:       res = func(mol);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       throw ValueErrorException("Bad AromaticityModel");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::setAromaticity
    validate_inputs(topology, rings)?;
    let assignment = match params.model {
        AromaticityModel::Rdkit => aromaticity_helper(topology, rings, 0, 0, true)?,
        AromaticityModel::Simple => aromaticity_helper(topology, rings, 5, 6, false)?,
        AromaticityModel::Mdl => mdl_aromaticity_helper(topology, rings)?,
        AromaticityModel::Mmff94 => mmff94_aromaticity_helper(topology, rings)?,
        AromaticityModel::Custom => {
            return Err(AromaticityError::UnsupportedModel {
                model: AromaticityModel::Custom,
                detail: "custom aromaticity requires a human-approved cross-language callback contract",
            });
        }
    };

    assignment.topology.validate()?;
    if assignment.topology.atoms.len() != topology.atoms.len()
        || assignment.topology.bonds.len() != topology.bonds.len()
    {
        return Err(AromaticityError::UnexpectedTopologyShape {
            input_atoms: topology.atoms.len(),
            input_bonds: topology.bonds.len(),
            output_atoms: assignment.topology.atoms.len(),
            output_bonds: assignment.topology.bonds.len(),
        });
    }
    for (row, (input, output)) in topology
        .atoms
        .iter()
        .zip(&assignment.topology.atoms)
        .enumerate()
    {
        if input.id() != output.id() {
            return Err(AromaticityError::AtomIdentityChanged {
                row,
                expected: input.id(),
                actual: output.id(),
            });
        }
    }
    for (row, (input, output)) in topology
        .bonds
        .iter()
        .zip(&assignment.topology.bonds)
        .enumerate()
    {
        if input.id() != output.id() {
            return Err(AromaticityError::BondIdentityChanged {
                row,
                expected: input.id(),
                actual: output.id(),
            });
        }
    }
    let maximum = if params.model == AromaticityModel::Mmff94 {
        rings
            .atom_rings()
            .len()
            .checked_add(1)
            .ok_or(AromaticityError::IntegerOverflow {
                field: "MMFF94 aromatic ring count upper bound",
            })?
    } else {
        rings.atom_rings().len()
    };
    if assignment.aromatic_ring_count > maximum {
        return Err(AromaticityError::AromaticRingCountOutOfRange {
            model: params.model,
            actual: assignment.aromatic_ring_count,
            maximum,
        });
    }
    Ok(assignment)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RingFindType;
    use cosmolkit_model::{AtomSpec, BondSpec};
    use cosmolkit_types::Element;

    fn topology(atom_specs: Vec<AtomSpec>, bonds: Vec<(usize, usize, BondOrder)>) -> TopologyBlock {
        let atoms = atom_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
            .collect();
        let bonds = bonds
            .into_iter()
            .enumerate()
            .map(|(index, (begin, end, order))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }

    fn valence(explicit: &[i32], implicit: &[i32]) -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: explicit.to_vec(),
            implicit_hydrogens: implicit.to_vec(),
        }
    }

    fn no_rings(topology: &TopologyBlock) -> RingInfo {
        RingInfo::new(
            RingFindType::Sssr,
            topology.atoms.len(),
            topology.bonds.len(),
        )
    }

    fn triangle_ring(topology: &TopologyBlock) -> RingInfo {
        let mut rings = no_rings(topology);
        rings.add_ring(&[0, 1, 2], &[0, 1, 2]).unwrap();
        rings
    }

    #[test]
    fn candidate_donor_classes_follow_the_source_acceptance_switch() {
        let graph = topology(
            vec![AtomSpec::new(Element::C).with_explicit_hydrogens(3)],
            vec![],
        );
        let rings = no_rings(&graph);
        let assigned = valence(&[3], &[0]);
        for donor in [
            ElectronDonorType::Vacant,
            ElectronDonorType::One,
            ElectronDonorType::Two,
            ElectronDonorType::OneOrTwo,
            ElectronDonorType::Any,
        ] {
            assert!(
                is_atom_candidate(
                    &graph,
                    &rings,
                    &assigned,
                    AtomId::new(0),
                    donor,
                    CandidateOptions::default(),
                )
                .unwrap(),
                "source donor {donor:?} should remain eligible"
            );
        }
        assert!(
            !is_atom_candidate(
                &graph,
                &rings,
                &assigned,
                AtomId::new(0),
                ElectronDonorType::None,
                CandidateOptions::default(),
            )
            .unwrap()
        );
    }

    #[test]
    fn candidate_dummy_donor_distinguishes_cyclic_multiple_bonds() {
        let acyclic = topology(
            vec![AtomSpec::new(Element::DUMMY), AtomSpec::new(Element::C)],
            vec![(0, 1, BondOrder::Double)],
        );
        assert_eq!(
            atom_donor_type(
                &acyclic,
                &no_rings(&acyclic),
                &valence(&[0, 2], &[0, 0]),
                AtomId::new(0),
                true,
            )
            .unwrap(),
            ElectronDonorType::Any
        );

        let cyclic = topology(
            vec![
                AtomSpec::new(Element::DUMMY),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
            ],
            vec![
                (0, 1, BondOrder::Double),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
            ],
        );
        assert_eq!(
            atom_donor_type(
                &cyclic,
                &triangle_ring(&cyclic),
                &valence(&[0, 3, 2], &[0, 0, 1]),
                AtomId::new(0),
                true,
            )
            .unwrap(),
            ElectronDonorType::One
        );
    }

    #[test]
    fn candidate_count_atom_electrons_covers_hydrogen_zero_and_dative_degree() {
        let carbon = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
            ],
            vec![(0, 1, BondOrder::Double), (0, 2, BondOrder::Single)],
        );
        assert_eq!(
            count_atom_electrons(&carbon, &valence(&[3, 2, 1], &[1, 0, 3]), AtomId::new(0))
                .unwrap(),
            1
        );

        let pyrrole_n = topology(
            vec![
                AtomSpec::new(Element::N).with_explicit_hydrogens(1),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
            ],
            vec![(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
        );
        assert_eq!(
            count_atom_electrons(&pyrrole_n, &valence(&[3, 1, 1], &[0, 3, 3]), AtomId::new(0),)
                .unwrap(),
            2
        );

        let zero_and_dative = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::N),
            ],
            vec![(0, 1, BondOrder::Zero), (0, 2, BondOrder::Dative)],
        );
        assert_eq!(
            count_atom_electrons(
                &zero_and_dative,
                &valence(&[0, 0, 1], &[0, 4, 2]),
                AtomId::new(0),
            )
            .unwrap(),
            4
        );
        assert!(!is_bond_order_query(&zero_and_dative.bonds[0]));
        assert!(!is_bond_order_query(&zero_and_dative.bonds[1]));

        let fluorine = topology(vec![AtomSpec::new(Element::F)], vec![]);
        assert_eq!(
            count_atom_electrons(&fluorine, &valence(&[0], &[0]), AtomId::new(0)).unwrap(),
            -1
        );
    }

    #[test]
    fn candidate_donor_charge_exocyclic_and_electronegativity_branches_are_exact() {
        let cation = topology(
            vec![
                AtomSpec::new(Element::C).with_formal_charge(1),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
            ],
            vec![
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );
        assert_eq!(
            atom_donor_type(
                &cation,
                &no_rings(&cation),
                &valence(&[3, 1, 1, 1], &[0, 3, 3, 3]),
                AtomId::new(0),
                true,
            )
            .unwrap(),
            ElectronDonorType::Vacant
        );

        let carbonyl = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::O),
                AtomSpec::new(Element::C),
            ],
            vec![(0, 1, BondOrder::Double), (0, 2, BondOrder::Single)],
        );
        let assigned = valence(&[3, 2, 1], &[1, 0, 3]);
        assert_eq!(
            atom_donor_type(
                &carbonyl,
                &no_rings(&carbonyl),
                &assigned,
                AtomId::new(0),
                true,
            )
            .unwrap(),
            ElectronDonorType::Vacant
        );
        assert_eq!(
            atom_donor_type(
                &carbonyl,
                &no_rings(&carbonyl),
                &assigned,
                AtomId::new(0),
                false,
            )
            .unwrap(),
            ElectronDonorType::One
        );
    }

    #[test]
    fn candidate_options_cover_rows_radicals_triples_and_exocyclic_bonds() {
        let isolated = |element: Element, charge: i8, radicals: u8| {
            topology(
                vec![
                    AtomSpec::new(element)
                        .with_formal_charge(charge)
                        .with_radical_electrons(radicals),
                ],
                vec![],
            )
        };
        for atomic_number in [34, 52] {
            let graph = isolated(Element::from_atomic_number(atomic_number).unwrap(), 0, 0);
            assert!(
                is_atom_candidate(
                    &graph,
                    &no_rings(&graph),
                    &valence(&[0], &[0]),
                    AtomId::new(0),
                    ElectronDonorType::One,
                    CandidateOptions::default(),
                )
                .unwrap()
            );
        }
        let gold = isolated(Element::AU, 0, 0);
        assert!(
            !is_atom_candidate(
                &gold,
                &no_rings(&gold),
                &valence(&[0], &[0]),
                AtomId::new(0),
                ElectronDonorType::One,
                CandidateOptions::default(),
            )
            .unwrap()
        );
        let sulfur = isolated(Element::S, 0, 0);
        assert!(
            !is_atom_candidate(
                &sulfur,
                &no_rings(&sulfur),
                &valence(&[0], &[0]),
                AtomId::new(0),
                ElectronDonorType::One,
                CandidateOptions {
                    allow_third_row: false,
                    ..CandidateOptions::default()
                },
            )
            .unwrap()
        );
        let oxygen = isolated(Element::O, 0, 0);
        assert!(
            !is_atom_candidate(
                &oxygen,
                &no_rings(&oxygen),
                &valence(&[0], &[0]),
                AtomId::new(0),
                ElectronDonorType::Two,
                CandidateOptions {
                    only_c_or_n: true,
                    ..CandidateOptions::default()
                },
            )
            .unwrap()
        );

        for (element, charge) in [(Element::N, 0), (Element::C, 1)] {
            let graph = isolated(element, charge, 1);
            assert!(
                !is_atom_candidate(
                    &graph,
                    &no_rings(&graph),
                    &valence(&[0], &[0]),
                    AtomId::new(0),
                    ElectronDonorType::One,
                    CandidateOptions::default(),
                )
                .unwrap()
            );
        }
        let neutral_radical_carbon = isolated(Element::C, 0, 1);
        assert!(
            is_atom_candidate(
                &neutral_radical_carbon,
                &no_rings(&neutral_radical_carbon),
                &valence(&[0], &[0]),
                AtomId::new(0),
                ElectronDonorType::One,
                CandidateOptions::default(),
            )
            .unwrap()
        );

        let triple = topology(
            vec![AtomSpec::new(Element::C), AtomSpec::new(Element::N)],
            vec![(0, 1, BondOrder::Triple)],
        );
        assert!(
            !is_atom_candidate(
                &triple,
                &no_rings(&triple),
                &valence(&[3, 3], &[0, 0]),
                AtomId::new(0),
                ElectronDonorType::One,
                CandidateOptions {
                    allow_triple_bonds: false,
                    ..CandidateOptions::default()
                },
            )
            .unwrap()
        );
        let double = topology(
            vec![AtomSpec::new(Element::C), AtomSpec::new(Element::O)],
            vec![(0, 1, BondOrder::Double)],
        );
        assert!(
            !is_atom_candidate(
                &double,
                &no_rings(&double),
                &valence(&[2, 2], &[2, 0]),
                AtomId::new(0),
                ElectronDonorType::One,
                CandidateOptions {
                    allow_exocyclic_multiple_bonds: false,
                    ..CandidateOptions::default()
                },
            )
            .unwrap()
        );
    }

    #[test]
    fn candidate_validation_rejects_malformed_detached_rows_before_indexing() {
        let graph = topology(vec![AtomSpec::new(Element::C)], vec![]);
        let wrong_rings = RingInfo::new(RingFindType::Sssr, 2, 0);
        assert!(matches!(
            validate_inputs(&graph, &wrong_rings),
            Err(AromaticityError::RingInfoDimensionMismatch {
                ring_atom_count: 2,
                topology_atom_count: 1,
                ..
            })
        ));
        assert!(matches!(
            validate_valence(&graph, &valence(&[], &[])),
            Err(AromaticityError::ValenceAssignmentLength {
                field: "explicit_valence",
                actual: 0,
                expected: 1,
            })
        ));
        assert!(matches!(
            incident_multiple_bond(&graph, &valence(&[], &[]), AtomId::new(0)),
            Err(AromaticityError::ValenceAssignmentLength { .. })
        ));
        assert!(matches!(
            atom_donor_type(
                &graph,
                &no_rings(&graph),
                &valence(&[0], &[0]),
                AtomId::new(3),
                true,
            ),
            Err(AromaticityError::AtomOutOfRange {
                atom,
                atom_count: 1,
            }) if atom == AtomId::new(3)
        ));
    }

    fn ring_info(topology: &TopologyBlock, atom_rings: &[Vec<usize>]) -> RingInfo {
        let mut rings = no_rings(topology);
        for atom_ring in atom_rings {
            let mut bond_ring = Vec::with_capacity(atom_ring.len());
            for index in 0..atom_ring.len() {
                let begin = atom_ring[index];
                let end = atom_ring[(index + 1) % atom_ring.len()];
                let bond = topology
                    .adjacency
                    .neighbors_of(begin)
                    .iter()
                    .find(|neighbor| neighbor.atom_index == end)
                    .expect("test ring edge must exist")
                    .bond;
                bond_ring.push(bond.index());
            }
            rings.add_ring(atom_ring, &bond_ring).unwrap();
        }
        rings
    }

    fn alternating_cycle(elements: &[Element]) -> (TopologyBlock, RingInfo) {
        let atoms = elements.iter().copied().map(AtomSpec::new).collect();
        let bonds = (0..elements.len())
            .map(|index| {
                (
                    index,
                    (index + 1) % elements.len(),
                    if index % 2 == 0 {
                        BondOrder::Double
                    } else {
                        BondOrder::Single
                    },
                )
            })
            .collect();
        let graph = topology(atoms, bonds);
        let ring = ring_info(&graph, &[(0..elements.len()).collect()]);
        (graph, ring)
    }

    fn alternating_sp2_cycle(elements: &[Element]) -> (TopologyBlock, RingInfo) {
        let atoms = elements
            .iter()
            .copied()
            .map(|element| AtomSpec::new(element).with_hybridization(Hybridization::Sp2))
            .collect();
        let bonds = (0..elements.len())
            .map(|index| {
                (
                    index,
                    (index + 1) % elements.len(),
                    if index % 2 == 0 {
                        BondOrder::Double
                    } else {
                        BondOrder::Single
                    },
                )
            })
            .collect();
        let graph = topology(atoms, bonds);
        let ring = ring_info(&graph, &[(0..elements.len()).collect()]);
        (graph, ring)
    }

    #[test]
    fn huckel_electron_ranges_any_limit_and_minimum_ring_size_are_exact() {
        let six: Vec<AtomId> = (0..6).map(AtomId::new).collect();
        assert!(apply_huckel(&six, &[ElectronDonorType::One; 6], 0).unwrap());
        assert!(!apply_huckel(&six, &[ElectronDonorType::One; 6], 7).unwrap());
        assert!(!apply_huckel(&six[..4], &[ElectronDonorType::One; 6], 0).unwrap());
        assert!(
            apply_huckel(
                &six[..1],
                &[
                    ElectronDonorType::OneOrTwo,
                    ElectronDonorType::None,
                    ElectronDonorType::None,
                    ElectronDonorType::None,
                    ElectronDonorType::None,
                    ElectronDonorType::None,
                ],
                0,
            )
            .unwrap()
        );
        assert!(
            !apply_huckel(
                &six[..2],
                &[
                    ElectronDonorType::Any,
                    ElectronDonorType::Any,
                    ElectronDonorType::None,
                    ElectronDonorType::None,
                    ElectronDonorType::None,
                    ElectronDonorType::None,
                ],
                0,
            )
            .unwrap()
        );
    }

    #[test]
    fn huckel_combination_order_and_terminal_state_follow_rdkit() {
        let mut combination = vec![0, 1];
        let mut observed = vec![combination.clone()];
        while next_combination(&mut combination, 4).is_some() {
            observed.push(combination.clone());
        }
        assert_eq!(
            observed,
            vec![
                vec![0, 1],
                vec![0, 2],
                vec![0, 3],
                vec![1, 2],
                vec![1, 3],
                vec![2, 3],
            ]
        );
        assert_eq!(next_combination(&mut [], 4), None);
    }

    #[test]
    fn huckel_neighbor_map_honors_size_overlap_and_subset_connectivity() {
        let bond_rings = vec![
            vec![BondId::new(0), BondId::new(1), BondId::new(2)],
            vec![BondId::new(2), BondId::new(3), BondId::new(4)],
            vec![BondId::new(3), BondId::new(4), BondId::new(5)],
            (10..35).map(BondId::new).collect(),
        ];
        let neighbors = make_ring_neighbor_map(&bond_rings, 24, 1);
        assert_eq!(neighbors, vec![vec![1], vec![0], vec![], vec![]]);
        assert!(check_fused(&[0, 1], &neighbors).unwrap());
        assert!(!check_fused(&[0, 2], &neighbors).unwrap());
        let mut done = vec![false; neighbors.len()];
        assert_eq!(
            pick_fused_rings(0, &neighbors, &mut done).unwrap(),
            vec![0, 1]
        );
        assert_eq!(done, vec![true, true, false, false]);
    }

    #[test]
    fn huckel_fused_marking_marks_the_complete_outer_envelope_only() {
        let graph = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
            ],
            vec![
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 0, BondOrder::Single),
                (1, 3, BondOrder::Single),
                (3, 0, BondOrder::Double),
            ],
        );
        let atom_rings = vec![
            vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)],
            vec![AtomId::new(0), AtomId::new(1), AtomId::new(3)],
        ];
        let bond_rings = convert_to_bonds(&graph, &atom_rings).unwrap();
        assert_eq!(bond_rings[0][0], bond_rings[1][0]);
        let shared = bond_rings[0][0].index();
        let mut working = graph.clone();
        let mut done = BTreeSet::new();
        mark_atoms_bonds_aromatic(&mut working, &bond_rings, &[0, 1], &mut done).unwrap();
        assert!(!working.bonds[shared].is_aromatic());
        assert_eq!(done.len(), 4);
        for (index, bond) in working.bonds.iter().enumerate() {
            assert_eq!(bond.is_aromatic(), index != shared);
            if index != shared {
                assert_eq!(bond.order(), BondOrder::Aromatic);
            }
        }
        assert!(working.atoms.iter().all(Atom::is_aromatic));
    }

    #[test]
    fn huckel_rdkit_fused_envelope_differs_from_simple_azulene_policy() {
        let graph = topology(
            vec![AtomSpec::new(Element::C); 10],
            vec![
                (0, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Double),
                (4, 1, BondOrder::Single),
                (1, 0, BondOrder::Single),
                (0, 5, BondOrder::Single),
                (5, 6, BondOrder::Double),
                (6, 7, BondOrder::Single),
                (7, 8, BondOrder::Double),
                (8, 9, BondOrder::Single),
                (9, 1, BondOrder::Double),
            ],
        );
        let rings = ring_info(&graph, &[vec![0, 2, 3, 4, 1], vec![0, 5, 6, 7, 8, 9, 1]]);
        let rdkit = aromaticity_helper(&graph, &rings, 0, 0, true).unwrap();
        let simple = aromaticity_helper(&graph, &rings, 5, 6, false).unwrap();
        assert_eq!(rdkit.aromatic_ring_count, 2);
        assert_eq!(simple.aromatic_ring_count, 0);
        assert_eq!(
            rdkit
                .topology
                .bonds
                .iter()
                .filter(|bond| bond.is_aromatic())
                .count(),
            10
        );
        assert!(!rdkit.topology.bonds[4].is_aromatic());
        assert!(simple.topology.bonds.iter().all(|bond| !bond.is_aromatic()));
        assert_eq!(
            graph.atoms.iter().filter(|atom| atom.is_aromatic()).count(),
            0
        );
    }

    #[test]
    fn huckel_mdl_accepts_six_carbon_one_electron_rows_and_rejects_hetero_mixed_rows() {
        let (benzene, benzene_rings) = alternating_cycle(&[Element::C; 6]);
        let positive = mdl_aromaticity_helper(&benzene, &benzene_rings).unwrap();
        assert_eq!(positive.aromatic_ring_count, 1);
        assert!(positive.topology.atoms.iter().all(Atom::is_aromatic));

        let hetero = topology(
            vec![
                AtomSpec::new(Element::N).with_explicit_hydrogens(1),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
            ],
            vec![
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Double),
                (4, 5, BondOrder::Single),
                (5, 0, BondOrder::Single),
            ],
        );
        let hetero_rings = ring_info(&hetero, &[vec![0, 1, 2, 3, 4, 5]]);
        let hetero_valence =
            crate::assign_valence_for_topology(&hetero, crate::ValenceModel::RdkitLike).unwrap();
        assert_eq!(
            atom_donor_type(
                &hetero,
                &hetero_rings,
                &hetero_valence,
                AtomId::new(0),
                false,
            )
            .unwrap(),
            ElectronDonorType::Two
        );
        let negative = mdl_aromaticity_helper(&hetero, &hetero_rings).unwrap();
        assert_eq!(negative.aromatic_ring_count, 0);
        assert!(
            negative
                .topology
                .bonds
                .iter()
                .all(|bond| !bond.is_aromatic())
        );

        let mixed = topology(
            vec![AtomSpec::new(Element::C); 8],
            vec![
                (0, 1, BondOrder::Double),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Double),
                (3, 4, BondOrder::Single),
                (4, 5, BondOrder::Double),
                (5, 0, BondOrder::Single),
                (2, 6, BondOrder::Single),
                (6, 7, BondOrder::Double),
                (7, 3, BondOrder::Single),
            ],
        );
        let mixed_rings = ring_info(&mixed, &[vec![0, 1, 2, 3, 4, 5], vec![2, 6, 7, 3]]);
        let mixed_result = mdl_aromaticity_helper(&mixed, &mixed_rings).unwrap();
        assert_eq!(mixed_result.aromatic_ring_count, 1);
    }

    #[test]
    fn mmff94_prearomatic_input_is_locally_kekulized_without_mutating_the_source() {
        let (mut graph, _) = alternating_sp2_cycle(&[Element::C; 6]);
        for atom in &mut graph.atoms {
            atom.set_aromatic(true);
        }
        for bond in &mut graph.bonds {
            bond.set_order(BondOrder::Aromatic);
            bond.set_aromatic(true);
        }
        let rings = ring_info(&graph, &[vec![0, 1, 2, 3, 4, 5]]);
        let source = graph.clone();

        let assignment = mmff94_aromaticity_helper(&graph, &rings).unwrap();

        assert_eq!(graph, source);
        assert_eq!(assignment.aromatic_ring_count, 2);
        assert!(assignment.topology.atoms.iter().all(Atom::is_aromatic));
        assert!(
            assignment
                .topology
                .bonds
                .iter()
                .all(|bond| { bond.is_aromatic() && bond.order() == BondOrder::Aromatic })
        );
        assert!(
            assignment
                .topology
                .atoms
                .iter()
                .all(|atom| atom.prop("_MMFFSanitized").is_none())
        );
        assert!(
            assignment
                .topology
                .bonds
                .iter()
                .all(|bond| bond.prop("_MMFFSanitized").is_none())
        );
    }

    #[test]
    fn mmff94_sp2_and_four_n_plus_two_boundaries_are_exact() {
        let (benzene, benzene_rings) = alternating_sp2_cycle(&[Element::C; 6]);
        let six_pi = mmff94_aromaticity_helper(&benzene, &benzene_rings).unwrap();
        assert_eq!(six_pi.aromatic_ring_count, 2);
        assert!(six_pi.topology.atoms.iter().all(Atom::is_aromatic));

        let (mut non_sp2, non_sp2_rings) = alternating_sp2_cycle(&[Element::C; 6]);
        non_sp2.atoms[3].set_hybridization(Hybridization::Sp3);
        let rejected = mmff94_aromaticity_helper(&non_sp2, &non_sp2_rings).unwrap();
        assert_eq!(rejected.aromatic_ring_count, 1);
        assert!(
            rejected
                .topology
                .atoms
                .iter()
                .all(|atom| !atom.is_aromatic())
        );

        let (four_pi, four_pi_rings) = alternating_sp2_cycle(&[Element::C; 4]);
        let antiaromatic = mmff94_aromaticity_helper(&four_pi, &four_pi_rings).unwrap();
        assert_eq!(antiaromatic.aromatic_ring_count, 1);
        assert!(
            antiaromatic
                .topology
                .bonds
                .iter()
                .all(|bond| !bond.is_aromatic())
        );

        let (ten_pi, ten_pi_rings) = alternating_sp2_cycle(&[Element::C; 10]);
        let larger = mmff94_aromaticity_helper(&ten_pi, &ten_pi_rings).unwrap();
        assert_eq!(larger.aromatic_ring_count, 2);
        assert!(larger.topology.bonds.iter().all(|bond| bond.is_aromatic()));

        let empty = topology(
            vec![AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2)],
            vec![],
        );
        assert_eq!(
            mmff94_aromaticity_helper(&empty, &no_rings(&empty))
                .unwrap()
                .aromatic_ring_count,
            1
        );
    }

    #[test]
    fn mmff94_n_o_and_divalent_s_odd_rings_receive_the_source_two_electron_adjustment() {
        for hetero in [Element::N, Element::O, Element::S] {
            let hetero_spec = AtomSpec::new(hetero)
                .with_hybridization(Hybridization::Sp2)
                .with_explicit_hydrogens(u8::from(hetero == Element::N));
            let graph = topology(
                vec![
                    hetero_spec,
                    AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
                    AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
                    AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
                    AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
                ],
                vec![
                    (0, 1, BondOrder::Single),
                    (1, 2, BondOrder::Double),
                    (2, 3, BondOrder::Single),
                    (3, 4, BondOrder::Double),
                    (4, 0, BondOrder::Single),
                ],
            );
            let rings = ring_info(&graph, &[vec![0, 1, 2, 3, 4]]);
            let assignment = mmff94_aromaticity_helper(&graph, &rings).unwrap();
            assert_eq!(assignment.aromatic_ring_count, 2, "hetero={hetero:?}");
            assert!(
                assignment
                    .topology
                    .bonds
                    .iter()
                    .all(|bond| bond.is_aromatic()),
                "hetero={hetero:?}"
            );
        }

        let trivalent_sulfur = topology(
            vec![
                AtomSpec::new(Element::S).with_hybridization(Hybridization::Sp2),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
                AtomSpec::new(Element::C),
            ],
            vec![
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Double),
                (4, 0, BondOrder::Single),
                (0, 5, BondOrder::Single),
            ],
        );
        let rings = ring_info(&trivalent_sulfur, &[vec![0, 1, 2, 3, 4]]);
        let assignment = mmff94_aromaticity_helper(&trivalent_sulfur, &rings).unwrap();
        assert_eq!(assignment.aromatic_ring_count, 1);
        assert!(
            assignment
                .topology
                .bonds
                .iter()
                .all(|bond| !bond.is_aromatic())
        );
    }

    #[test]
    fn mmff94_exocyclic_double_bonds_and_odd_electron_rows_follow_source_branches() {
        let ring_atoms = vec![
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::N).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
        ];
        let ring_bonds = vec![
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Single),
            (3, 4, BondOrder::Double),
            (4, 0, BondOrder::Single),
        ];
        let control = topology(ring_atoms.clone(), ring_bonds.clone());
        let control_rings = ring_info(&control, &[vec![0, 1, 2, 3, 4]]);
        assert_eq!(
            mmff94_aromaticity_helper(&control, &control_rings)
                .unwrap()
                .aromatic_ring_count,
            2
        );

        let mut exocyclic_atoms = ring_atoms;
        exocyclic_atoms.push(AtomSpec::new(Element::O));
        let mut exocyclic_bonds = ring_bonds;
        exocyclic_bonds.push((0, 5, BondOrder::Double));
        let exocyclic = topology(exocyclic_atoms, exocyclic_bonds);
        let exocyclic_rings = ring_info(&exocyclic, &[vec![0, 1, 2, 3, 4]]);
        let rejected = mmff94_aromaticity_helper(&exocyclic, &exocyclic_rings).unwrap();
        assert_eq!(rejected.aromatic_ring_count, 1);
        assert!(
            rejected.topology.atoms[..5]
                .iter()
                .all(|atom| !atom.is_aromatic())
        );

        let (mut radical, radical_rings) = alternating_sp2_cycle(&[Element::C; 6]);
        radical.atoms[0].set_radical_electrons(1);
        let radical_assignment = mmff94_aromaticity_helper(&radical, &radical_rings).unwrap();
        assert_eq!(radical_assignment.aromatic_ring_count, 2);
        assert!(
            radical_assignment
                .topology
                .atoms
                .iter()
                .all(Atom::is_aromatic)
        );
    }

    #[test]
    fn mmff94_ring_dependencies_are_deferred_until_both_exocyclic_neighbors_are_perceived() {
        let mut atoms = vec![AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2); 16];
        for &sulfur in &[6usize, 11] {
            atoms[sulfur] = AtomSpec::new(Element::S).with_hybridization(Hybridization::Sp2);
        }
        for &nitrogen in &[8usize, 13] {
            atoms[nitrogen] = AtomSpec::new(Element::N).with_hybridization(Hybridization::Sp2);
        }
        let graph = topology(
            atoms,
            vec![
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Single),
                (4, 5, BondOrder::Double),
                (5, 0, BondOrder::Single),
                (6, 7, BondOrder::Single),
                (7, 8, BondOrder::Double),
                (8, 9, BondOrder::Single),
                (9, 10, BondOrder::Double),
                (10, 6, BondOrder::Single),
                (11, 12, BondOrder::Single),
                (12, 13, BondOrder::Double),
                (13, 14, BondOrder::Single),
                (14, 15, BondOrder::Double),
                (15, 11, BondOrder::Single),
                (0, 6, BondOrder::Double),
                (3, 11, BondOrder::Double),
            ],
        );
        let rings = ring_info(
            &graph,
            &[
                vec![0, 1, 2, 3, 4, 5],
                vec![6, 7, 8, 9, 10],
                vec![11, 12, 13, 14, 15],
            ],
        );

        let assignment = mmff94_aromaticity_helper(&graph, &rings).unwrap();

        assert_eq!(assignment.aromatic_ring_count, 4);
        assert!(assignment.topology.atoms.iter().all(Atom::is_aromatic));
        assert!(!assignment.topology.bonds[16].is_aromatic());
        assert!(!assignment.topology.bonds[17].is_aromatic());
    }

    #[test]
    fn mmff94_noncarbon_hydrogen_adjustment_and_kekulization_failure_are_atomic() {
        let c = || AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2);
        let cationic_n = AtomSpec::new(Element::N)
            .with_formal_charge(1)
            .with_hybridization(Hybridization::Sp2);
        let graph = topology(
            vec![cationic_n, c(), c(), c(), c()],
            vec![
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Double),
                (4, 0, BondOrder::Single),
            ],
        );
        let rings = ring_info(&graph, &[vec![0, 1, 2, 3, 4]]);
        let assignment = mmff94_aromaticity_helper(&graph, &rings).unwrap();
        assert_eq!(assignment.aromatic_ring_count, 2);
        assert_eq!(assignment.topology.atoms[0].explicit_hydrogens(), 1);
        assert_eq!(graph.atoms[0].explicit_hydrogens(), 0);

        let mut impossible = topology(
            vec![c(), c(), c()],
            vec![
                (0, 1, BondOrder::Aromatic),
                (1, 2, BondOrder::Aromatic),
                (2, 0, BondOrder::Aromatic),
            ],
        );
        for atom in &mut impossible.atoms {
            atom.set_aromatic(true);
        }
        for bond in &mut impossible.bonds {
            bond.set_aromatic(true);
        }
        let impossible_rings = ring_info(&impossible, &[vec![0, 1, 2]]);
        let snapshot = impossible.clone();
        assert!(matches!(
            mmff94_aromaticity_helper(&impossible, &impossible_rings),
            Err(AromaticityError::Kekulize(
                KekulizeError::NotKekulizable { .. }
            ))
        ));
        assert_eq!(impossible, snapshot);
    }
}
