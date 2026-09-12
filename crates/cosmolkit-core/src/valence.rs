//! RDKit-aligned valence primitives over detached model values.
//!
//! This module is deliberately below the runtime boundary.  It computes a
//! value result from [`cosmolkit_model::TopologyBlock`] and never reads or
//! updates a live molecule property cache.

use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, Bond, BondId, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::BondOrder;

use crate::periodic_table;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValenceModel {
    RdkitLike,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ValenceParams {
    pub model: ValenceModel,
    pub strict: bool,
}

impl Default for ValenceParams {
    fn default() -> Self {
        Self {
            model: ValenceModel::RdkitLike,
            strict: true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValencePhase {
    EffectiveAtomicNumber,
    Explicit,
    Implicit,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValenceAssignment {
    pub explicit_valence: Vec<i32>,
    pub implicit_hydrogens: Vec<i32>,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum ValenceError {
    #[error("{message}")]
    InvalidValence {
        atom: AtomId,
        atomic_number: u8,
        formal_charge: i8,
        phase: ValencePhase,
        calculated: Option<i32>,
        reason: &'static str,
        message: String,
    },
    #[error("invalid topology: {source}")]
    InvalidTopology { source: TopologyValidationError },
    #[error("atom {atom} is out of range for {atom_count} atoms")]
    AtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error(
        "adjacency for atom {atom} references neighbor atom row {neighbor_atom}, out of range for {atom_count} atoms"
    )]
    AdjacencyAtomOutOfRange {
        atom: AtomId,
        neighbor_atom: usize,
        atom_count: usize,
    },
    #[error(
        "adjacency for atom {atom} references bond {bond}, out of range for {bond_count} bonds"
    )]
    AdjacencyBondOutOfRange {
        atom: AtomId,
        bond: BondId,
        bond_count: usize,
    },
    #[error(
        "adjacency for atom {atom} references neighbor row {neighbor_atom} through bond {bond}, but that bond has endpoints {begin}-{end}"
    )]
    AdjacencyEndpointMismatch {
        atom: AtomId,
        neighbor_atom: usize,
        bond: BondId,
        begin: AtomId,
        end: AtomId,
    },
    #[error("explicit valence input for atom {atom} must be nonnegative, got {value}")]
    InvalidExplicitValenceInput { atom: AtomId, value: i32 },
    #[error("periodic-table field {field} is unavailable for atomic number {atomic_number}")]
    PeriodicTableLookup {
        atomic_number: u8,
        field: &'static str,
    },
    #[error("explicit valence is not available for atom {atom}")]
    ExplicitValenceCacheNotInitialized { atom: AtomId },
    #[error("implicit valence is not available for atom {atom}")]
    ImplicitValenceCacheNotInitialized { atom: AtomId },
    #[error(
        "hydrogen count overflow at atom {atom}: explicit={explicit}, implicit={implicit}, neighbor_hydrogens={neighbor_hydrogens}"
    )]
    HydrogenCountOverflow {
        atom: AtomId,
        explicit: u32,
        implicit: u32,
        neighbor_hydrogens: usize,
    },
    #[error("Bad bond type")]
    BadBondType {
        bond: Option<BondId>,
        order: BondOrder,
    },
}

fn atom_from_parts(atoms: &[Atom], atom_id: AtomId) -> Result<&Atom, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION ROMol::getAtomWithIdx
    // RDKit✔️✔️: //! returns a pointer to a particular Atom
    // RDKit✔️✔️: Atom *getAtomWithIdx(unsigned int idx);
    // RDKit✔️✔️: //! \overload
    // RDKit✔️✔️: const Atom *getAtomWithIdx(unsigned int idx) const;
    // END RDKIT CPP FUNCTION ROMol::getAtomWithIdx
    atoms
        .get(atom_id.index())
        .ok_or(ValenceError::AtomOutOfRange {
            atom: atom_id,
            atom_count: atoms.len(),
        })
}

fn incident_bonds_from_parts<'a>(
    atom_count: usize,
    bonds: &'a [Bond],
    adjacency: &'a AdjacencyList,
    atom_id: AtomId,
) -> Result<impl Iterator<Item = &'a Bond> + 'a, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION ROMol::atomBonds
    // RDKit✔️✔️: CXXBondIterator<const MolGraph, Bond *const, MolGraph::out_edge_iterator>
    // RDKit✔️✔️: atomBonds(Atom const *at) const {
    // RDKit✔️✔️:   auto pr = getAtomBonds(at);
    // RDKit✔️✔️:   return {&d_graph, pr.first, pr.second};
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ROMol::atomBonds
    if atom_id.index() >= atom_count {
        return Err(ValenceError::AtomOutOfRange {
            atom: atom_id,
            atom_count,
        });
    }
    let neighbors = adjacency.neighbors_of(atom_id.index());
    for neighbor in neighbors {
        if neighbor.atom_index >= atom_count {
            return Err(ValenceError::AdjacencyAtomOutOfRange {
                atom: atom_id,
                neighbor_atom: neighbor.atom_index,
                atom_count,
            });
        }
        let Some(bond) = bonds.get(neighbor.bond.index()) else {
            return Err(ValenceError::AdjacencyBondOutOfRange {
                atom: atom_id,
                bond: neighbor.bond,
                bond_count: bonds.len(),
            });
        };
        let neighbor_id = AtomId::new(neighbor.atom_index);
        let matches_endpoints = (bond.begin() == atom_id && bond.end() == neighbor_id)
            || (bond.end() == atom_id && bond.begin() == neighbor_id);
        if !matches_endpoints {
            return Err(ValenceError::AdjacencyEndpointMismatch {
                atom: atom_id,
                neighbor_atom: neighbor.atom_index,
                bond: neighbor.bond,
                begin: bond.begin(),
                end: bond.end(),
            });
        }
    }
    Ok(neighbors
        .iter()
        .map(move |neighbor| &bonds[neighbor.bond.index()]))
}

fn atom(topology: &TopologyBlock, id: AtomId) -> Result<&Atom, ValenceError> {
    atom_from_parts(&topology.atoms, id)
}

pub(crate) fn incident<'a>(
    topology: &'a TopologyBlock,
    id: AtomId,
) -> Result<impl Iterator<Item = &'a Bond> + 'a, ValenceError> {
    incident_bonds_from_parts(
        topology.atoms.len(),
        &topology.bonds,
        &topology.adjacency,
        id,
    )
}

fn validate_topology(topology: &TopologyBlock) -> Result<(), ValenceError> {
    topology
        .validate()
        .map_err(|source| ValenceError::InvalidTopology { source })
}

fn validate_atom_id(topology: &TopologyBlock, atom_id: AtomId) -> Result<(), ValenceError> {
    if atom_id.index() >= topology.atoms.len() {
        return Err(ValenceError::AtomOutOfRange {
            atom: atom_id,
            atom_count: topology.atoms.len(),
        });
    }
    Ok(())
}

pub fn bond_type_as_double(order: BondOrder) -> Result<f64, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION Bond::getBondTypeAsDouble
    // RDKit✔️✔️: double Bond::getBondTypeAsDouble() const {
    // RDKit✔️✔️:   double res;
    // RDKit✔️✔️:   switch (getBondType()) {
    let value = match order {
        // RDKit✔️✔️:     case UNSPECIFIED:
        // RDKit✔️✔️:     case IONIC:
        // RDKit✔️✔️:     case ZERO:
        // RDKit✔️✔️:       res = 0;
        BondOrder::Unspecified | BondOrder::Ionic | BondOrder::Zero => 0.0,
        // RDKit✔️✔️:     case SINGLE:
        // RDKit✔️✔️:       res = 1;
        BondOrder::Single => 1.0,
        // RDKit✔️✔️:     case DOUBLE:
        // RDKit✔️✔️:       res = 2;
        BondOrder::Double => 2.0,
        // RDKit✔️✔️:     case TRIPLE:
        // RDKit✔️✔️:       res = 3;
        BondOrder::Triple => 3.0,
        // RDKit✔️✔️:     case QUADRUPLE:
        // RDKit✔️✔️:       res = 4;
        BondOrder::Quadruple => 4.0,
        // RDKit✔️✔️:     case QUINTUPLE:
        // RDKit✔️✔️:       res = 5;
        BondOrder::Quintuple => 5.0,
        // RDKit✔️✔️:     case HEXTUPLE:
        // RDKit✔️✔️:       res = 6;
        BondOrder::Hextuple => 6.0,
        // RDKit✔️✔️:     case ONEANDAHALF:
        // RDKit✔️✔️:       res = 1.5;
        BondOrder::OneAndHalf => 1.5,
        // RDKit✔️✔️:     case TWOANDAHALF:
        // RDKit✔️✔️:       res = 2.5;
        BondOrder::TwoAndHalf => 2.5,
        // RDKit✔️✔️:     case THREEANDAHALF:
        // RDKit✔️✔️:       res = 3.5;
        BondOrder::ThreeAndHalf => 3.5,
        // RDKit✔️✔️:     case FOURANDAHALF:
        // RDKit✔️✔️:       res = 4.5;
        BondOrder::FourAndHalf => 4.5,
        // RDKit✔️✔️:     case FIVEANDAHALF:
        // RDKit✔️✔️:       res = 5.5;
        BondOrder::FiveAndHalf => 5.5,
        // RDKit✔️✔️:     case AROMATIC:
        // RDKit✔️✔️:       res = 1.5;
        BondOrder::Aromatic => 1.5,
        // RDKit✔️✔️:     case DATIVEONE:
        // RDKit✔️✔️:       res = 1.0;
        // RDKit✔️✔️:       break;  // FIX: this should probably be different
        // RDKit✔️✔️:     case DATIVE:
        // RDKit✔️✔️:       res = 1.0;
        // RDKit✔️✔️:       break;  // FIX: again probably wrong
        BondOrder::Dative | BondOrder::DativeOne => 1.0,
        // RDKit✔️✔️:     case HYDROGEN:
        // RDKit✔️✔️:       res = 0.0;
        BondOrder::Hydrogen => 0.0,
        // RDKit✔️✔️:     default:
        // RDKit✔️✔️:       UNDER_CONSTRUCTION("Bad bond type");
        BondOrder::DativeLeft
        | BondOrder::DativeRight
        | BondOrder::ThreeCenter
        | BondOrder::Other => {
            return Err(ValenceError::BadBondType { bond: None, order });
        }
    };
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Bond::getBondTypeAsDouble
    Ok(value)
}

/// Return one bond's valence contribution at a selected endpoint.
pub fn bond_valence_contrib(bond: &Bond, atom: AtomId) -> Result<f64, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION Bond::getValenceContrib
    // RDKit✔️✔️: double Bond::getValenceContrib(const Atom *atom) const {
    // RDKit✔️✔️:   if (atom != getBeginAtom() && atom != getEndAtom()) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    if bond.begin() != atom && bond.end() != atom {
        return Ok(0.0);
    }
    // RDKit✔️✔️:   double res;
    // RDKit✔️✔️:   if ((getBondType() == DATIVE || getBondType() == DATIVEONE) &&
    // RDKit✔️✔️:       atom->getIdx() != getEndAtomIdx()) {
    // RDKit✔️✔️:     res = 0.0;
    // RDKit✔️✔️:   } else if (getBondType() == HYDROGEN) {
    // RDKit✔️✔️:     res = 0.0;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = getBondTypeAsDouble();
    // RDKit✔️✔️:   }
    let result =
        if matches!(bond.order(), BondOrder::Dative | BondOrder::DativeOne) && bond.end() != atom {
            0.0
        } else if bond.order() == BondOrder::Hydrogen {
            0.0
        } else {
            bond_type_as_double(bond.order()).map_err(|error| match error {
                ValenceError::BadBondType { order, .. } => ValenceError::BadBondType {
                    bond: Some(bond.id()),
                    order,
                },
                other => other,
            })?
        };
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Bond::getValenceContrib
    Ok(result)
}

pub fn get_effective_atomic_num(atom: &Atom, check_value: bool) -> Result<u8, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION getEffectiveAtomicNum
    // RDKit✔️✔️: unsigned int getEffectiveAtomicNum(const Atom &atom, bool checkValue) {
    // RDKit✔️✔️:   auto effectiveAtomicNum = atom.getAtomicNum() - atom.getFormalCharge();
    let value = i32::from(atom.atomic_number()) - i32::from(atom.formal_charge());
    // RDKit✔️✔️:   if (checkValue &&
    // RDKit✔️✔️:       (effectiveAtomicNum < 0 ||
    // RDKit✔️✔️:        effectiveAtomicNum >
    // RDKit✔️✔️:            static_cast<int>(PeriodicTable::getTable()->getMaxAtomicNumber()))) {
    // RDKit✔️✔️:     throw AtomValenceException("Effective atomic number out of range",
    // RDKit✔️✔️:                                atom.getIdx());
    // RDKit✔️✔️:   }
    if check_value && !(0..=118).contains(&value) {
        return Err(invalid_valence_with_message(
            atom,
            ValencePhase::EffectiveAtomicNumber,
            Some(value),
            "effective atomic number out of range",
            "Effective atomic number out of range".to_string(),
        ));
    }
    // RDKit✔️✔️:   effectiveAtomicNum = std::clamp(
    // RDKit✔️✔️:       effectiveAtomicNum, 0,
    // RDKit✔️✔️:       static_cast<int>(PeriodicTable::getTable()->getMaxAtomicNumber()));
    // RDKit✔️✔️:   return static_cast<unsigned int>(effectiveAtomicNum);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getEffectiveAtomicNum
    Ok(value.clamp(0, 118) as u8)
}

pub fn can_be_hypervalent(atom: &Atom, effective_atomic_num: u8) -> bool {
    // BEGIN RDKIT CPP FUNCTION canBeHypervalent
    // RDKit✔️✔️: bool canBeHypervalent(const Atom &atom, unsigned int effectiveAtomicNum) {
    // RDKit✔️✔️:   return (effectiveAtomicNum > 16 &&
    // RDKit✔️✔️:           (atom.getAtomicNum() == 15 || atom.getAtomicNum() == 16)) ||
    // RDKit✔️✔️:          (effectiveAtomicNum > 34 &&
    // RDKit✔️✔️:           (atom.getAtomicNum() == 33 || atom.getAtomicNum() == 34));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION canBeHypervalent
    (effective_atomic_num > 16 && matches!(atom.atomic_number(), 15 | 16))
        || (effective_atomic_num > 34 && matches!(atom.atomic_number(), 33 | 34))
}

fn is_aromatic_atom_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
) -> Result<bool, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION isAromaticAtom
    // RDKit✔️✔️: bool isAromaticAtom(const Atom &atom) {
    // RDKit✔️✔️:   if (atom.getIsAromatic()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    if atom_from_parts(atoms, atom_id)?.is_aromatic() {
        return Ok(true);
    }
    // RDKit✔️✔️:   if (atom.hasOwningMol()) {
    // RDKit✔️✔️:     for (const auto &bond : atom.getOwningMol().atomBonds(&atom)) {
    // RDKit✔️✔️:       if (bond->getIsAromatic() ||
    // RDKit✔️✔️:           bond->getBondType() == Bond::BondType::AROMATIC) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    for bond in incident_bonds_from_parts(atoms.len(), bonds, adjacency, atom_id)? {
        if bond.is_aromatic() || bond.order() == BondOrder::Aromatic {
            return Ok(true);
        }
    }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isAromaticAtom
    Ok(false)
}

pub fn required_valence_list(atomic_number: u8) -> Result<&'static [i32], ValenceError> {
    periodic_table::valences(atomic_number).ok_or(ValenceError::PeriodicTableLookup {
        atomic_number,
        field: "valences",
    })
}

pub fn rdkit_default_valence(atomic_number: u8) -> Result<i32, ValenceError> {
    Ok(required_valence_list(atomic_number)?[0])
}

pub fn periodic_table_row(atomic_number: u8) -> Option<u8> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getRow / atomicData::Row
    // RDKit✔️✔️: UINT getRow(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].Row();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int Row() const { return row; }
    // END RDKIT CPP FUNCTION PeriodicTable::getRow / atomicData::Row
    periodic_table::period(atomic_number)
}

pub fn periodic_table_outer_electrons(atomic_number: u8) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION atomicData::NumOuterShellElec / PeriodicTable::getNouterElecs
    // RDKit✔️✔️: int NumOuterShellElec() const { return nVal; }
    // RDKit✔️✔️: int getNouterElecs(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].NumOuterShellElec();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT extern const std::string periodicTableAtomData;
    // END RDKIT CPP FUNCTION atomicData::NumOuterShellElec / PeriodicTable::getNouterElecs
    periodic_table::outer_electrons(atomic_number).ok_or(ValenceError::PeriodicTableLookup {
        atomic_number,
        field: "outer_electrons",
    })
}

pub fn periodic_table_more_electronegative(
    atomic_number_1: u8,
    atomic_number_2: u8,
) -> Result<bool, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::moreElectroNegative
    // RDKit✔️✔️: bool moreElectroNegative(UINT anum1, UINT anum2) const {
    // RDKit✔️✔️:   PRECONDITION(anum1 < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   PRECONDITION(anum2 < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   // FIX: the atomic_data needs to have real electronegativity values
    // RDKit✔️✔️:   UINT ne1 = getNouterElecs(anum1);
    // RDKit✔️✔️:   UINT ne2 = getNouterElecs(anum2);
    // RDKit✔️✔️:   if (ne1 > ne2) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (ne1 == ne2) {
    // RDKit✔️✔️:     if (anum1 < anum2) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PeriodicTable::moreElectroNegative
    periodic_table::more_electronegative(atomic_number_1, atomic_number_2).ok_or(
        ValenceError::PeriodicTableLookup {
            atomic_number: atomic_number_1.max(atomic_number_2),
            field: "electronegativity",
        },
    )
}

#[must_use]
pub fn rdkit_rb0(atomic_number: u8) -> f64 {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getRb0 / atomicData::Rb0
    // RDKit✔️✔️: double getRb0(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].Rb0();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: double Rb0() const { return rB0; }
    // END RDKIT CPP FUNCTION PeriodicTable::getRb0 / atomicData::Rb0
    periodic_table::rb0(atomic_number)
}

#[must_use]
pub fn rdkit_atomic_number_from_symbol(symbol: &str) -> Option<u8> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getAtomicNumber
    // RDKit✔️🔝: int getAtomicNumber(const std::string &elementSymbol) const {
    // RDKit✔️🔝:   int anum = -1;
    // RDKit✔️🔝:   if (elementSymbol == "C") {
    // RDKit✔️🔝:     anum = 6;
    // RDKit✔️🔝:   } else if (elementSymbol == "N") {
    // RDKit✔️🔝:     anum = 7;
    // RDKit✔️🔝:   } else if (elementSymbol == "O") {
    // RDKit✔️🔝:     anum = 8;
    // RDKit✔️🔝:   } else {
    // RDKit✔️🔝:     STR_UINT_MAP::const_iterator iter = byname.find(elementSymbol);
    // RDKit✔️🔝:     if (iter != byname.end()) {
    // RDKit✔️🔝:       anum = iter->second;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   POSTCONDITION(anum > -1, "Element '" + elementSymbol + "' not found");
    // RDKit✔️🔝:   return anum;
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION PeriodicTable::getAtomicNumber
    // The model's generated static symbol match avoids RDKit's O(log n)
    // std::map lookup without allocation and preserves canonical and legacy
    // aliases for the modeled periodic-table rows.
    periodic_table::atomic_number_from_symbol(symbol)
}

fn invalid_valence_with_message(
    atom: &Atom,
    phase: ValencePhase,
    calculated: Option<i32>,
    reason: &'static str,
    message: String,
) -> ValenceError {
    ValenceError::InvalidValence {
        atom: atom.id(),
        atomic_number: atom.atomic_number(),
        formal_charge: atom.formal_charge(),
        phase,
        calculated,
        reason,
        message,
    }
}

fn calculate_explicit_valence(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
    check_it: bool,
) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION calculateExplicitValence
    // RDKit✔️✔️: int calculateExplicitValence(const Atom &atom, bool strict, bool checkIt) {
    // RDKit✔️✔️:   // FIX: contributions of bonds to valence are being done at best
    // RDKit✔️✔️:   // approximately
    // RDKit✔️✔️:   double accum = 0;
    // RDKit✔️✔️:   for (const auto bnd : atom.getOwningMol().atomBonds(&atom)) {
    // RDKit✔️✔️:     accum += bnd->getValenceContrib(&atom);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   accum += atom.getNumExplicitHs();
    let atom = atom_from_parts(atoms, atom_id)?;
    let mut accum = 0.0;
    for bond in incident_bonds_from_parts(atoms.len(), bonds, adjacency, atom_id)? {
        accum += bond_valence_contrib(bond, atom_id)?;
    }
    accum += f64::from(atom.explicit_hydrogens());

    // RDKit✔️✔️:   const auto &ovalens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(atom.getAtomicNum());
    // RDKit✔️✔️:   // if we start with an atom that doesn't have specified valences, we stick
    // RDKit✔️✔️:   // with that. otherwise we will use the effective valence
    // RDKit✔️✔️:   unsigned int effectiveAtomicNum = atom.getAtomicNum();
    // RDKit✔️✔️:   if (ovalens.size() > 1 || ovalens[0] != -1) {
    // RDKit✔️✔️:     effectiveAtomicNum = getEffectiveAtomicNum(atom, checkIt);
    // RDKit✔️✔️:   }
    let ovalens = required_valence_list(atom.atomic_number())?;
    let mut effective_atomic_num = atom.atomic_number();
    if ovalens.len() > 1 || ovalens[0] != -1 {
        effective_atomic_num = get_effective_atomic_num(atom, check_it)?;
    }

    // RDKit✔️✔️:   unsigned int dv =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getDefaultValence(effectiveAtomicNum);
    // RDKit✔️✔️:   const auto &valens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(effectiveAtomicNum);
    let default_valence = rdkit_default_valence(effective_atomic_num)?;
    let valens = required_valence_list(effective_atomic_num)?;

    // RDKit✔️✔️:   if (accum > dv && isAromaticAtom(atom)) {
    // RDKit✔️✔️:     // this needs some explanation : if the atom is aromatic and
    // RDKit✔️✔️:     // accum > dv we assume that no hydrogen can be added
    // RDKit✔️✔️:     // to this atom.  We set x = (v + chr) such that x is the
    // RDKit✔️✔️:     // closest possible integer to "accum" but less than
    // RDKit✔️✔️:     // "accum".
    // RDKit✔️✔️:     //
    // RDKit✔️✔️:     // "v" here is one of the allowed valences. For example:
    // RDKit✔️✔️:     //    sulfur here : O=c1ccs(=O)cc1
    // RDKit✔️✔️:     //    nitrogen here : c1cccn1C
    if accum > f64::from(default_valence)
        && is_aromatic_atom_from_parts(atoms, bonds, adjacency, atom_id)?
    {
        // RDKit✔️✔️:     int pval = dv;
        // RDKit✔️✔️:     for (auto val : valens) {
        // RDKit✔️✔️:       if (val == -1) {
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (val > accum) {
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         pval = val;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let mut pval = default_valence;
        for &valence in valens {
            if valence == -1 || f64::from(valence) > accum {
                break;
            }
            pval = valence;
        }
        // RDKit✔️✔️:     // if we're within 1.5 of the allowed valence, go ahead and take it.
        // RDKit✔️✔️:     // this reflects things like the N in c1cccn1C, which starts with
        // RDKit✔️✔️:     // accum of 4, but which can be kekulized to C1=CC=CN1C, where
        // RDKit✔️✔️:     // the valence is 3 or the bridging N in c1ccn2cncc2c1, which starts
        // RDKit✔️✔️:     // with a valence of 4.5, but can be happily kekulized down to a valence
        // RDKit✔️✔️:     // of 3
        // RDKit✔️✔️:     if (accum - pval <= 1.5) {
        // RDKit✔️✔️:       accum = pval;
        // RDKit✔️✔️:     }
        if accum - f64::from(pval) <= 1.5 {
            accum = f64::from(pval);
        }
    }

    // RDKit✔️✔️:   // despite promising to not to blame it on him - this a trick Greg
    // RDKit✔️✔️:   // came up with: if we have a bond order sum of x.5 (i.e. 1.5, 2.5
    // RDKit✔️✔️:   // etc) we would like it to round to the higher integer value --
    // RDKit✔️✔️:   // 2.5 to 3 instead of 2 -- so we will add 0.1 to accum.
    // RDKit✔️✔️:   // this plays a role in the number of hydrogen that are implicitly
    // RDKit✔️✔️:   // added. This will only happen when the accum is a non-integer
    // RDKit✔️✔️:   // value and less than the default valence (otherwise the above if
    // RDKit✔️✔️:   // statement should have caught it). An example of where this can
    // RDKit✔️✔️:   // happen is the following smiles:
    // RDKit✔️✔️:   //    C1ccccC1
    // RDKit✔️✔️:   // Daylight accepts this smiles and we should be able to Kekulize
    // RDKit✔️✔️:   // correctly.
    // RDKit✔️✔️:   accum += 0.1;
    // RDKit✔️✔️:   auto res = static_cast<int>(std::round(accum));
    accum += 0.1;
    let result = accum.round() as i32;

    // RDKit✔️✔️:   if (strict || checkIt) {
    // RDKit✔️✔️:     int maxValence = valens.back();
    // RDKit✔️✔️:     int offset = 0;
    if strict || check_it {
        let mut max_valence = *valens.last().expect("valence list is nonempty");
        let mut offset = 0;
        // RDKit✔️✔️:     // we have to include a special case here for negatively charged P, S, As,
        // RDKit✔️✔️:     // and Se, which all support "hypervalent" forms, but which can be
        // RDKit✔️✔️:     // isoelectronic to Cl/Ar or Br/Kr, which do not support hypervalent forms.
        // RDKit✔️✔️:     if (canBeHypervalent(atom, effectiveAtomicNum)) {
        // RDKit✔️✔️:       maxValence = ovalens.back();
        // RDKit✔️✔️:       offset -= atom.getFormalCharge();
        // RDKit✔️✔️:     }
        if can_be_hypervalent(atom, effective_atomic_num) {
            max_valence = *ovalens.last().expect("valence list is nonempty");
            offset -= i32::from(atom.formal_charge());
        }
        // RDKit✔️✔️:     // we have historically accepted two-coordinate [H-] as a valid atom. This
        // RDKit✔️✔️:     // is highly questionable, but changing it requires some thought. For now we
        // RDKit✔️✔️:     // will just keep accepting it
        // RDKit✔️✔️:     if (atom.getAtomicNum() == 1 && atom.getFormalCharge() == -1) {
        // RDKit✔️✔️:       maxValence = 2;
        // RDKit✔️✔️:     }
        if atom.atomic_number() == 1 && atom.formal_charge() == -1 {
            max_valence = 2;
        }
        // RDKit✔️✔️:     // maxValence == -1 signifies that we'll take anything at the high end
        // RDKit✔️✔️:     if (maxValence >= 0 && ovalens.back() >= 0 && (res + offset) > maxValence) {
        // RDKit✔️✔️:       // the explicit valence is greater than any
        // RDKit✔️✔️:       // allowed valence for the atoms
        // RDKit✔️✔️:       if (strict) {
        // RDKit✔️✔️:         // raise an error
        // RDKit✔️✔️:         std::ostringstream errout;
        // RDKit✔️✔️:         errout << "Explicit valence for atom # " << atom.getIdx() << " "
        // RDKit✔️✔️:                << PeriodicTable::getTable()->getElementSymbol(
        // RDKit✔️✔️:                       atom.getAtomicNum())
        // RDKit✔️✔️:                << ", " << res << ", is greater than permitted";
        // RDKit✔️✔️:         std::string msg = errout.str();
        // RDKit✔️✔️:         BOOST_LOG(rdErrorLog) << msg << std::endl;
        // RDKit✔️✔️:         throw AtomValenceException(msg, atom.getIdx());
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         return -1;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        if max_valence >= 0
            && *ovalens.last().expect("valence list is nonempty") >= 0
            && result + offset > max_valence
        {
            if strict {
                return Err(invalid_valence_with_message(
                    atom,
                    ValencePhase::Explicit,
                    Some(result),
                    "greater than permitted",
                    format!(
                        "Explicit valence for atom # {} {}, {}, is greater than permitted",
                        atom.id(),
                        rdkit_element_symbol(atom.atomic_number())?,
                        result
                    ),
                ));
            }
            return Ok(-1);
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION calculateExplicitValence
    Ok(result)
}

/// Calculate one atom's explicit valence from a detached topology.
///
/// `check_it` exposes RDKit's validation-only sentinel path: with
/// `strict == false`, an invalid valence returns `-1` instead of throwing.
pub fn calculate_explicit_valence_for_topology(
    topology: &TopologyBlock,
    atom_id: AtomId,
    strict: bool,
    check_it: bool,
) -> Result<i32, ValenceError> {
    calculate_explicit_valence(
        &topology.atoms,
        &topology.bonds,
        &topology.adjacency,
        atom_id,
        strict,
        check_it,
    )
}

/// Calculate one atom's explicit valence through the canonical detached API.
pub fn explicit_valence_for_atom(
    topology: &TopologyBlock,
    atom_id: AtomId,
    strict: bool,
) -> Result<i32, ValenceError> {
    validate_atom_id(topology, atom_id)?;
    validate_topology(topology)?;
    calculate_explicit_valence_for_topology(topology, atom_id, strict, false)
}

/// Calculate one atom's explicit valence from detached topology parts.
pub fn calculate_explicit_valence_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
    check_it: bool,
) -> Result<i32, ValenceError> {
    calculate_explicit_valence(atoms, bonds, adjacency, atom_id, strict, check_it)
}

fn calculate_implicit_valence(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    explicit_valence: i32,
    strict: bool,
    check_it: bool,
) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION calculateImplicitValence
    // RDKit✔️✔️: int calculateImplicitValence(const Atom &atom, bool strict, bool checkIt) {
    // RDKit✔️✔️:   if (atom.df_noImplicit) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    let atom = atom_from_parts(atoms, atom_id)?;
    if atom.no_implicit() {
        return Ok(0);
    }
    // RDKit✔️✔️:   auto explicitValence = atom.d_explicitValence;
    // RDKit✔️✔️:   if (explicitValence == -1) {
    // RDKit✔️✔️:     explicitValence = calculateExplicitValence(atom, strict, checkIt);
    // RDKit✔️✔️:   }
    let explicit_valence = if explicit_valence == -1 {
        calculate_explicit_valence(atoms, bonds, adjacency, atom_id, strict, check_it)?
    } else {
        explicit_valence
    };
    // RDKit✔️✔️:   // special cases
    // RDKit✔️✔️:   auto atomicNum = atom.d_atomicNum;
    // RDKit✔️✔️:   if (atomicNum == 0) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    let atomic_num = atom.atomic_number();
    if atomic_num == 0 {
        return Ok(0);
    }
    // RDKit✔️✔️:   for (const auto bnd : atom.getOwningMol().atomBonds(&atom)) {
    // RDKit✔️✔️:     if (QueryOps::hasComplexBondTypeQuery(*bnd)) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // `TopologyBlock` contains concrete bonds, so a query-bearing bond cannot
    // cross this detached API boundary. Query graphs remain in cosmolkit-model
    // as a distinct value type.
    // RDKit✔️✔️:   auto formalCharge = atom.d_formalCharge;
    // RDKit✔️✔️:   auto numRadicalElectrons = atom.d_numRadicalElectrons;
    // RDKit✔️✔️:   if (explicitValence == 0 && numRadicalElectrons == 0 && atomicNum == 1) {
    // RDKit✔️✔️:     if (formalCharge == 1 || formalCharge == -1) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     } else if (formalCharge == 0) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       if (strict) {
    // RDKit✔️✔️:         std::ostringstream errout;
    // RDKit✔️✔️:         errout << "Unreasonable formal charge on atom # " << atom.getIdx()
    // RDKit✔️✔️:                << ".";
    // RDKit✔️✔️:         std::string msg = errout.str();
    // RDKit✔️✔️:         BOOST_LOG(rdErrorLog) << msg << std::endl;
    // RDKit✔️✔️:         throw AtomValenceException(msg, atom.getIdx());
    // RDKit✔️✔️:       } else if (checkIt) {
    // RDKit✔️✔️:         return -1;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         return 0;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if explicit_valence == 0 && atom.radical_electrons() == 0 && atomic_num == 1 {
        return match atom.formal_charge() {
            1 | -1 => Ok(0),
            0 => Ok(1),
            _ if strict => Err(invalid_valence_with_message(
                atom,
                ValencePhase::Implicit,
                Some(explicit_valence),
                "unreasonable hydrogen formal charge",
                format!("Unreasonable formal charge on atom # {}.", atom.id()),
            )),
            _ if check_it => Ok(-1),
            _ => Ok(0),
        };
    }

    // RDKit✔️✔️:   int explicitPlusRadV = atom.d_explicitValence + atom.d_numRadicalElectrons;
    let mut explicit_plus_rad_v = explicit_valence + i32::from(atom.radical_electrons());

    // RDKit✔️✔️:   const auto &ovalens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(atom.d_atomicNum);
    // RDKit✔️✔️:   // if we start with an atom that doesn't have specified valences, we stick
    // RDKit✔️✔️:   // with that. otherwise we will use the effective valence for the rest of
    // RDKit✔️✔️:   // this.
    // RDKit✔️✔️:   unsigned int effectiveAtomicNum = atom.d_atomicNum;
    // RDKit✔️✔️:   if (ovalens.size() > 1 || ovalens[0] != -1) {
    // RDKit✔️✔️:     effectiveAtomicNum = getEffectiveAtomicNum(atom, checkIt);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (effectiveAtomicNum == 0) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    let ovalens = required_valence_list(atomic_num)?;
    let mut effective_atomic_num = atomic_num;
    if ovalens.len() > 1 || ovalens[0] != -1 {
        effective_atomic_num = get_effective_atomic_num(atom, check_it)?;
    }
    if effective_atomic_num == 0 {
        return Ok(0);
    }

    // RDKit✔️✔️:   // The d-block and f-block of the periodic table (i.e. transition metals,
    // RDKit✔️✔️:   // lanthanoids and actinoids) have no default valence.
    // RDKit✔️✔️:   int dv = PeriodicTable::getTable()->getDefaultValence(effectiveAtomicNum);
    // RDKit✔️✔️:   if (dv == -1) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    let default_valence = rdkit_default_valence(effective_atomic_num)?;
    if default_valence == -1 {
        return Ok(0);
    }

    // RDKit✔️✔️:   // here is how we are going to deal with the possibility of
    // RDKit✔️✔️:   // multiple valences
    // RDKit✔️✔️:   // - check the explicit valence "ev"
    // RDKit✔️✔️:   // - if it is already equal to one of the allowed valences for the
    // RDKit✔️✔️:   //    atom return 0
    // RDKit✔️✔️:   // - otherwise take return difference between next larger allowed
    // RDKit✔️✔️:   //   valence and "ev"
    // RDKit✔️✔️:   // if "ev" is greater than all allowed valences for the atom raise an
    // RDKit✔️✔️:   // exception
    // RDKit✔️✔️:   // finally aromatic cases are dealt with differently - these atoms are allowed
    // RDKit✔️✔️:   // only default valences
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // we have to include a special case here for negatively charged P, S, As,
    // RDKit✔️✔️:   // and Se, which all support "hypervalent" forms, but which can be
    // RDKit✔️✔️:   // isoelectronic to Cl/Ar or Br/Kr, which do not support hypervalent forms.
    // RDKit✔️✔️:   if (canBeHypervalent(atom, effectiveAtomicNum)) {
    // RDKit✔️✔️:     effectiveAtomicNum = atomicNum;
    // RDKit✔️✔️:     explicitPlusRadV -= atom.d_formalCharge;
    // RDKit✔️✔️:   }
    if can_be_hypervalent(atom, effective_atomic_num) {
        effective_atomic_num = atomic_num;
        explicit_plus_rad_v -= i32::from(atom.formal_charge());
    }
    // RDKit✔️✔️:   const auto &valens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(effectiveAtomicNum);
    let valens = required_valence_list(effective_atomic_num)?;

    // RDKit✔️✔️:   int res = 0;
    let result;
    // RDKit✔️✔️:   // if we have an aromatic case treat it differently
    // RDKit✔️✔️:   if (isAromaticAtom(atom)) {
    if is_aromatic_atom_from_parts(atoms, bonds, adjacency, atom_id)? {
        // RDKit✔️✔️:     if (explicitPlusRadV <= dv) {
        // RDKit✔️✔️:       res = dv - explicitPlusRadV;
        // RDKit✔️✔️:     } else {
        if explicit_plus_rad_v <= default_valence {
            result = default_valence - explicit_plus_rad_v;
        } else {
            // RDKit✔️✔️:       bool satis = false;
            // RDKit✔️✔️:       for (auto vi = valens.begin(); vi != valens.end() && *vi > 0; ++vi) {
            // RDKit✔️✔️:         if (explicitPlusRadV == *vi) {
            // RDKit✔️✔️:           satis = true;
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            let satisfied = valens
                .iter()
                .take_while(|&&valence| valence > 0)
                .any(|&valence| explicit_plus_rad_v == valence);
            // RDKit✔️✔️:       if (!satis && (strict || checkIt)) {
            // RDKit✔️✔️:         if (strict) {
            // RDKit✔️✔️:           std::ostringstream errout;
            // RDKit✔️✔️:           errout << "Explicit valence for aromatic atom # " << atom.getIdx()
            // RDKit✔️✔️:                  << " not equal to any accepted valence\n";
            // RDKit✔️✔️:           std::string msg = errout.str();
            // RDKit✔️✔️:           BOOST_LOG(rdErrorLog) << msg << std::endl;
            // RDKit✔️✔️:           throw AtomValenceException(msg, atom.getIdx());
            // RDKit✔️✔️:         } else {
            // RDKit✔️✔️:           return -1;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            if !satisfied && (strict || check_it) {
                if strict {
                    return Err(invalid_valence_with_message(
                        atom,
                        ValencePhase::Implicit,
                        Some(explicit_plus_rad_v),
                        "aromatic valence is not accepted",
                        format!(
                            "Explicit valence for aromatic atom # {} not equal to any accepted valence\n",
                            atom.id()
                        ),
                    ));
                }
                return Ok(-1);
            }
            // RDKit✔️✔️:       res = 0;
            result = 0;
        }
    } else {
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     // non-aromatic case we are allowed to have non default valences
        // RDKit✔️✔️:     // and be able to add hydrogens
        // RDKit✔️✔️:     res = -1;
        // RDKit✔️✔️:     for (auto vi = valens.begin(); vi != valens.end() && *vi >= 0; ++vi) {
        // RDKit✔️✔️:       int tot = *vi;
        // RDKit✔️✔️:       if (explicitPlusRadV <= tot) {
        // RDKit✔️✔️:         res = tot - explicitPlusRadV;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let mut candidate = -1;
        for &valence in valens.iter().take_while(|&&valence| valence >= 0) {
            if explicit_plus_rad_v <= valence {
                candidate = valence - explicit_plus_rad_v;
                break;
            }
        }
        // RDKit✔️✔️:     if (res < 0) {
        // RDKit✔️✔️:       if ((strict || checkIt) && valens.back() != -1 && ovalens.back() > 0) {
        // RDKit✔️✔️:         // this means that the explicit valence is greater than any
        // RDKit✔️✔️:         // allowed valence for the atoms
        // RDKit✔️✔️:         if (strict) {
        // RDKit✔️✔️:           std::ostringstream errout;
        // RDKit✔️✔️:           errout << "Explicit valence for atom # " << atom.getIdx() << " "
        // RDKit✔️✔️:                  << PeriodicTable::getTable()->getElementSymbol(atomicNum)
        // RDKit✔️✔️:                  << " greater than permitted";
        // RDKit✔️✔️:           std::string msg = errout.str();
        // RDKit✔️✔️:           BOOST_LOG(rdErrorLog) << msg << std::endl;
        // RDKit✔️✔️:           throw AtomValenceException(msg, atom.getIdx());
        // RDKit✔️✔️:         } else {
        // RDKit✔️✔️:           return -1;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         res = 0;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        if candidate < 0 {
            if (strict || check_it)
                && *valens.last().expect("valence list is nonempty") != -1
                && *ovalens.last().expect("valence list is nonempty") > 0
            {
                if strict {
                    return Err(invalid_valence_with_message(
                        atom,
                        ValencePhase::Implicit,
                        Some(explicit_plus_rad_v),
                        "greater than permitted",
                        format!(
                            "Explicit valence for atom # {} {} greater than permitted",
                            atom.id(),
                            rdkit_element_symbol(atomic_num)?
                        ),
                    ));
                }
                return Ok(-1);
            }
            candidate = 0;
        }
        result = candidate;
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION calculateImplicitValence
    Ok(result)
}

/// Calculate one atom's implicit valence from a detached topology.
pub fn calculate_implicit_valence_for_topology(
    topology: &TopologyBlock,
    atom_id: AtomId,
    explicit_valence: i32,
    strict: bool,
    check_it: bool,
) -> Result<i32, ValenceError> {
    calculate_implicit_valence(
        &topology.atoms,
        &topology.bonds,
        &topology.adjacency,
        atom_id,
        explicit_valence,
        strict,
        check_it,
    )
}

/// Calculate one atom's implicit valence through the canonical detached API.
pub fn implicit_valence_for_atom(
    topology: &TopologyBlock,
    atom_id: AtomId,
    explicit_valence: Option<i32>,
    strict: bool,
) -> Result<i32, ValenceError> {
    validate_atom_id(topology, atom_id)?;
    validate_topology(topology)?;
    if let Some(value) = explicit_valence
        && value < 0
    {
        return Err(ValenceError::InvalidExplicitValenceInput {
            atom: atom_id,
            value,
        });
    }
    calculate_implicit_valence_for_topology(
        topology,
        atom_id,
        explicit_valence.unwrap_or(-1),
        strict,
        false,
    )
}

/// Calculate one atom's implicit valence from detached topology parts.
pub fn calculate_implicit_valence_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    explicit_valence: i32,
    strict: bool,
    check_it: bool,
) -> Result<i32, ValenceError> {
    calculate_implicit_valence(
        atoms,
        bonds,
        adjacency,
        atom_id,
        explicit_valence,
        strict,
        check_it,
    )
}

/// Assign explicit valence and implicit hydrogens for detached topology.
pub fn assign_valence_for_topology(
    topology: &TopologyBlock,
    model: ValenceModel,
) -> Result<ValenceAssignment, ValenceError> {
    assign_valence(
        topology,
        &ValenceParams {
            model,
            strict: true,
        },
    )
}

pub fn assign_valence_with_options_for_topology(
    topology: &TopologyBlock,
    model: ValenceModel,
    strict: bool,
) -> Result<ValenceAssignment, ValenceError> {
    assign_valence(topology, &ValenceParams { model, strict })
}

/// Assign explicit valence and implicit hydrogens in stable atom-row order.
pub fn assign_valence(
    topology: &TopologyBlock,
    params: &ValenceParams,
) -> Result<ValenceAssignment, ValenceError> {
    validate_topology(topology)?;
    assign_valence_with_options_from_parts(
        &topology.atoms,
        &topology.bonds,
        &topology.adjacency,
        params.model,
        params.strict,
    )
}

/// Assign valence from borrowed detached topology parts without cloning them.
pub fn assign_valence_with_options_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    _model: ValenceModel,
    strict: bool,
) -> Result<ValenceAssignment, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION ROMol::updatePropertyCache / Atom::updatePropertyCache
    // RDKit✔️✔️: void ROMol::updatePropertyCache(bool strict) {
    // RDKit✔️✔️:   for (auto atom : atoms()) {
    // RDKit✔️✔️:     atom->updatePropertyCache(strict);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto bond : bonds()) {
    // RDKit✔️✔️:     bond->updatePropertyCache(strict);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: void Atom::updatePropertyCache(bool strict) {
    // RDKit✔️✔️:   calcExplicitValence(strict);
    // RDKit✔️✔️:   calcImplicitValence(strict);
    // RDKit✔️✔️: }
    let mut explicit_values = Vec::with_capacity(atoms.len());
    let mut implicit_hydrogens = Vec::with_capacity(atoms.len());
    for atom in atoms {
        let explicit =
            calculate_explicit_valence(atoms, bonds, adjacency, atom.id(), strict, false)?;
        explicit_values.push(explicit);
        implicit_hydrogens.push(calculate_implicit_valence(
            atoms,
            bonds,
            adjacency,
            atom.id(),
            explicit,
            strict,
            false,
        )?);
    }
    // RDKit✔️✔️: void updatePropertyCache(bool strict = true) { (void)strict; }
    // END RDKIT CPP FUNCTION ROMol::updatePropertyCache / Atom::updatePropertyCache
    Ok(ValenceAssignment {
        explicit_valence: explicit_values,
        implicit_hydrogens,
    })
}

/// Assign both valence fields for one atom from borrowed detached parts.
pub fn assign_valence_state_for_atom_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
) -> Result<(i32, i32), ValenceError> {
    let explicit = calculate_explicit_valence(atoms, bonds, adjacency, atom_id, strict, false)?;
    let implicit =
        calculate_implicit_valence(atoms, bonds, adjacency, atom_id, explicit, strict, false)?;
    Ok((explicit, implicit))
}

/// Assign one atom's explicit valence from borrowed detached parts.
pub fn assign_explicit_valence_for_atom_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
) -> Result<i32, ValenceError> {
    calculate_explicit_valence(atoms, bonds, adjacency, atom_id, strict, false)
}

/// Assign one atom's implicit valence using an explicit detached value.
pub fn assign_implicit_valence_for_atom_from_parts_with_explicit_valence(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    explicit_valence: i32,
    strict: bool,
) -> Result<i32, ValenceError> {
    calculate_implicit_valence(
        atoms,
        bonds,
        adjacency,
        atom_id,
        explicit_valence,
        strict,
        false,
    )
}

/// Return the source-backed preferred valence list for an atomic number.
pub fn rdkit_valence_list(atomic_number: u8) -> Result<Option<&'static [i32]>, ValenceError> {
    if atomic_number > 118 {
        return Err(ValenceError::PeriodicTableLookup {
            atomic_number,
            field: "valences",
        });
    }
    Ok(periodic_table::valences(atomic_number))
}

pub fn rdkit_element_symbol(atomic_number: u8) -> Result<&'static str, ValenceError> {
    periodic_table::symbol(atomic_number).ok_or(ValenceError::PeriodicTableLookup {
        atomic_number,
        field: "symbol",
    })
}

pub fn atom_has_valence_violation_for_topology(
    topology: &TopologyBlock,
    id: AtomId,
) -> Result<bool, ValenceError> {
    atom_has_valence_violation_from_parts(&topology.atoms, &topology.bonds, &topology.adjacency, id)
}

/// Classify one atom through the canonical detached violation API.
pub fn has_valence_violation(topology: &TopologyBlock, id: AtomId) -> Result<bool, ValenceError> {
    validate_atom_id(topology, id)?;
    validate_topology(topology)?;
    atom_has_valence_violation_for_topology(topology, id)
}

/// Check one atom for a valence violation using borrowed detached parts.
pub fn atom_has_valence_violation_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    id: AtomId,
) -> Result<bool, ValenceError> {
    let atom = atom_from_parts(atoms, id)?;
    // BEGIN RDKIT CPP FUNCTION Atom::hasValenceViolation
    // RDKit✔️✔️: bool Atom::hasValenceViolation() const {
    // RDKit✔️✔️: if (getAtomicNum() == 0 || hasQuery() ||
    // RDKit✔️✔️:     std::any_of(bonds.begin(), bonds.end(), is_query)) {
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // Query atoms and bonds cannot occur in a concrete `TopologyBlock`.
    if atom.atomic_number() == 0 {
        return Ok(false);
    }
    // RDKit✔️✔️:   unsigned int effectiveAtomicNum;
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     bool checkIt = true;
    // RDKit✔️✔️:     effectiveAtomicNum = getEffectiveAtomicNum(*this, checkIt);
    // RDKit✔️✔️:   } catch (const AtomValenceException &) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    let effective_atomic_num = match get_effective_atomic_num(atom, true) {
        Ok(value) => value,
        Err(ValenceError::InvalidValence { .. }) => return Ok(true),
        Err(error) => return Err(error),
    };
    // RDKit✔️✔️:   // special case for H:
    // RDKit✔️✔️:   if (getAtomicNum() == 1) {
    // RDKit✔️✔️:     if (getFormalCharge() > 1 || getFormalCharge() < -1) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     if (getFormalCharge() > getAtomicNum() ||
    // RDKit✔️✔️:         PeriodicTable::getTable()->getRow(d_atomicNum) !=
    // RDKit✔️✔️:             PeriodicTable::getTable()->getRow(effectiveAtomicNum)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if atom.atomic_number() == 1 {
        if !(-1..=1).contains(&atom.formal_charge()) {
            return Ok(true);
        }
    } else if atom.formal_charge() > atom.atomic_number() as i8
        || periodic_table_row(atom.atomic_number()) != periodic_table_row(effective_atomic_num)
    {
        return Ok(true);
    }
    // RDKit✔️✔️:   bool strict = false;
    // RDKit✔️✔️:   bool checkIt = true;
    // RDKit✔️✔️:   if (calculateExplicitValence(*this, strict, checkIt) == -1 ||
    // RDKit✔️✔️:       calculateImplicitValence(*this, strict, checkIt) == -1) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::hasValenceViolation
    let explicit = calculate_explicit_valence(atoms, bonds, adjacency, id, false, true)?;
    if explicit == -1 {
        return Ok(true);
    }
    Ok(calculate_implicit_valence(atoms, bonds, adjacency, id, explicit, false, true)? == -1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{AdjacencyList, AtomSpec, BondId, BondSpec};
    use cosmolkit_types::Element;

    #[test]
    fn detached_assignment_matches_basic_ethanol_valence() {
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::O)),
        ];
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
            ),
        ];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let assignment =
            assign_valence_for_topology(&topology, ValenceModel::RdkitLike).expect("CCO valence");
        assert_eq!(assignment.explicit_valence, vec![1, 2, 1]);
        assert_eq!(assignment.implicit_hydrogens, vec![3, 2, 1]);
    }

    #[test]
    fn valence_list_uses_shared_periodic_table_without_runtime_state() {
        assert_eq!(rdkit_valence_list(6).unwrap(), Some(&[4][..]));
        assert!(rdkit_valence_list(119).is_err());
    }

    #[test]
    fn valence_violation_is_scoped_to_the_requested_atom() {
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(
                AtomId::new(1),
                AtomSpec::new(Element::C).with_formal_charge(7),
            ),
        ];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(2, &[]),
            atoms,
            ..TopologyBlock::default()
        };
        assert!(!atom_has_valence_violation_for_topology(&topology, AtomId::new(0)).unwrap());
        assert!(atom_has_valence_violation_for_topology(&topology, AtomId::new(1)).unwrap());
    }
}
