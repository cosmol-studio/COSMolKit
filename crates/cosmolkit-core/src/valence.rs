//! RDKit-aligned valence primitives over detached model values.
//!
//! This module is deliberately below the runtime boundary.  It computes a
//! value result from [`cosmolkit_model::TopologyBlock`] and never reads or
//! updates a live molecule property cache.

use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, Bond, BondId, SourceAtomValenceFacts, TopologyBlock,
    TopologyValidationError,
};
use cosmolkit_types::BondOrder;

use crate::periodic_table;

pub use cosmolkit_model::{AtomMetadata, ValenceError, ValencePhase};

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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValenceAssignment {
    pub explicit_valence: Vec<i32>,
    pub implicit_hydrogens: Vec<i32>,
}

// Closed borrowed atom getter projection: both source Atom and QueryAtom
// inherit the same nonvirtual Atom cache algorithms. No query/identity coercion.
trait SourceValenceAtom {
    fn id(&self) -> AtomId;
    fn atomic_number(&self) -> u8;
    fn formal_charge(&self) -> i8;
    fn explicit_hydrogens(&self) -> u8;
    fn no_implicit(&self) -> bool;
    fn radical_electrons(&self) -> u8;
    fn is_aromatic(&self) -> bool;
}
impl SourceValenceAtom for Atom {
    fn id(&self) -> AtomId {
        Atom::id(self)
    }
    fn atomic_number(&self) -> u8 {
        Atom::atomic_number(self)
    }
    fn formal_charge(&self) -> i8 {
        Atom::formal_charge(self)
    }
    fn explicit_hydrogens(&self) -> u8 {
        Atom::explicit_hydrogens(self)
    }
    fn no_implicit(&self) -> bool {
        Atom::no_implicit(self)
    }
    fn radical_electrons(&self) -> u8 {
        Atom::radical_electrons(self)
    }
    fn is_aromatic(&self) -> bool {
        Atom::is_aromatic(self)
    }
}
impl SourceValenceAtom for cosmolkit_model::QueryAtom {
    fn id(&self) -> AtomId {
        cosmolkit_model::QueryAtom::id(self)
    }
    fn atomic_number(&self) -> u8 {
        cosmolkit_model::QueryAtom::atomic_number(self)
    }
    fn formal_charge(&self) -> i8 {
        cosmolkit_model::QueryAtom::formal_charge(self)
    }
    fn explicit_hydrogens(&self) -> u8 {
        cosmolkit_model::QueryAtom::explicit_hydrogens(self)
    }
    fn no_implicit(&self) -> bool {
        cosmolkit_model::QueryAtom::no_implicit(self)
    }
    fn radical_electrons(&self) -> u8 {
        cosmolkit_model::QueryAtom::radical_electrons(self)
    }
    fn is_aromatic(&self) -> bool {
        cosmolkit_model::QueryAtom::is_aromatic(self)
    }
}

fn atom_from_parts<A: SourceValenceAtom>(atoms: &[A], atom_id: AtomId) -> Result<&A, ValenceError> {
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
    get_effective_atomic_num_source(atom, check_value)
}

fn get_effective_atomic_num_source<A: SourceValenceAtom>(
    atom: &A,
    check_value: bool,
) -> Result<u8, ValenceError> {
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
    can_be_hypervalent_source(atom, effective_atomic_num)
}

fn can_be_hypervalent_source<A: SourceValenceAtom>(atom: &A, effective_atomic_num: u8) -> bool {
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

fn is_aromatic_atom_from_parts<A: SourceValenceAtom>(
    atoms: &[A],
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
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getValenceList complete source
    // RDKit✔️✔️:   const INT_VECT &getValenceList(UINT atomicNumber) const {
    // RDKit✔️✔️:     PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:     return byanum[atomicNumber].ValenceList();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: const INT_VECT &ValenceList() const { return valence; }
    // END RDKIT CPP FUNCTION PeriodicTable::getValenceList complete source
    // Borrow the complete canonical signed list in source insertion order,
    // including terminal -1. Range failure is a checked precondition error.
    // Native vector reference and Rust immutable slice both use one table
    // index without a per-call allocation, scan, sort, or copy.

    periodic_table::valences(atomic_number).ok_or(ValenceError::PeriodicTableLookup {
        atomic_number,
        field: "valences",
    })
}

pub fn rdkit_default_valence(atomic_number: u8) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getDefaultValence complete source
    // RDKit✔️✔️:   int getDefaultValence(UINT atomicNumber) const {
    // RDKit✔️✔️:     PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:     return byanum[atomicNumber].DefaultValence();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: int DefaultValence() const { return valence.front(); }
    // END RDKIT CPP FUNCTION PeriodicTable::getDefaultValence complete source
    // Canonical source rows 0..=118 all have a nonempty signed valence list.
    // Return its first entry verbatim, including -1 and 0; never infer a
    // preferred valence from period, charge, or the last allowed valence.
    // Range failure remains a structured checked precondition error.
    // One shared immutable-table lookup and one indexed i32 load; no clone,
    // scan, allocation, or alternate periodic data are introduced.

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

/// Source postcondition failure retaining the original symbol bytes.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("Element '{}' not found", String::from_utf8_lossy(.symbol.as_bytes()))]
pub struct AtomicSymbolLookupError {
    pub symbol: cosmolkit_model::PropertyText,
}

pub fn rdkit_atomic_number_from_symbol(
    symbol: impl AsRef<[u8]>,
) -> Result<u8, AtomicSymbolLookupError> {
    // RDKit❗✔️:   int getAtomicNumber(const std::string &elementSymbol) const {
    // RDKit❗✔️:     // this little optimization actually makes a measurable difference
    // RDKit❗✔️:     // in molecule-construction time
    // RDKit❗✔️:     int anum = -1;
    // RDKit❗✔️:     if (elementSymbol == "C") {
    // RDKit❗✔️:       anum = 6;
    // RDKit❗✔️:     } else if (elementSymbol == "N") {
    // RDKit❗✔️:       anum = 7;
    // RDKit❗✔️:     } else if (elementSymbol == "O") {
    // RDKit❗✔️:       anum = 8;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       STR_UINT_MAP::const_iterator iter = byname.find(elementSymbol);
    // RDKit❗✔️:       if (iter != byname.end()) {
    // RDKit❗✔️:         anum = iter->second;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     POSTCONDITION(anum > -1, "Element '" + elementSymbol + "' not found");
    // RDKit❗✔️:     return anum;
    // RDKit❗✔️:   }
    // Reuse the sole static vocabulary populated from all pinned byname rows,
    // including the dummy and Uut/Uup aliases. The source map keys are all ASCII:
    // invalid UTF8 cannot equal any key and follows the same failed postcondition.
    // Embedded NUL is part of std::string, never trimmed or truncated here.
    // Keep the native C/N/O short-circuit branches. Other lookups borrow input;
    // static vocabulary dispatch avoids table/map allocation, while error owns
    // only its exact symbol bytes as the native postcondition message does.
    let symbol = symbol.as_ref();
    let atomic_number = match symbol {
        b"C" => Some(6),
        b"N" => Some(7),
        b"O" => Some(8),
        _ => periodic_table::atomic_number_from_symbol(symbol),
    };
    atomic_number.ok_or_else(|| AtomicSymbolLookupError {
        symbol: cosmolkit_model::PropertyText::from_bytes(symbol),
    })
}

/// The source char-pointer overload over a valid NUL-terminated byte string.
pub fn rdkit_atomic_number_from_c_symbol(
    symbol: &std::ffi::CStr,
) -> Result<u8, AtomicSymbolLookupError> {
    // RDKit❗✔️:   int getAtomicNumber(const char *elementSymbol) const {
    // RDKit❗✔️:     std::string symb(elementSymbol);
    // RDKit❗✔️:
    // RDKit❗✔️:     return getAtomicNumber(symb);
    // RDKit❗✔️:   }
    // CStr supplies the same valid char-pointer input bytes before the first
    // NUL used by std::string(elementSymbol). Reuse the sole checked string
    // lookup. Bytes after the terminator are never included in lookup/errors.
    // Borrowing removes the native temporary string copy; no second symbol
    // matcher, trimming, decoding fallback or fabricated empty symbol exists.
    rdkit_atomic_number_from_symbol(symbol.to_bytes())
}

/// Native symbolic-valence precondition or delegated numeric lookup failure.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SymbolValenceLookupError {
    #[error("Element '{}' not found", String::from_utf8_lossy(.symbol.as_bytes()))]
    SymbolNotFound {
        symbol: cosmolkit_model::PropertyText,
    },
    #[error(transparent)]
    NumericLookup(#[from] ValenceError),
}

pub fn rdkit_default_valence_from_symbol(
    symbol: impl AsRef<[u8]>,
) -> Result<i32, SymbolValenceLookupError> {
    // RDKit❗✔️:   int getDefaultValence(const std::string &elementSymbol) const {
    // RDKit❗✔️:     PRECONDITION(byname.count(elementSymbol),
    // RDKit❗✔️:                  "Element '" + elementSymbol + "' not found");
    // RDKit❗✔️:     return getDefaultValence(byname.find(elementSymbol)->second);
    // RDKit❗✔️:   }
    // RDKit❗✔️: #define PRECONDITION(expr, mess)                                           \
    // RDKit❗✔️:   if (!(expr)) {                                                           \
    // RDKit❗✔️:     Invar::Invariant inv("Pre-condition Violation", mess, #expr, __FILE__, \
    // RDKit❗✔️:                          __LINE__);                                        \
    // RDKit❗✔️:     BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n" << std::endl; \
    // RDKit❗✔️:     throw inv;                                                             \
    // RDKit❗✔️:   }
    // PRECONDITION fails before numeric lookup. This error is distinct from
    // getAtomicNumber's POSTCONDITION; no guessed atomic number or valence.
    // Reuse the canonical static byname projection and sole numeric default
    // owner, retaining its first signed entry including -1 and the dummy row.
    // Immutable byname count/find can share one lookup without changing state;
    // no map/table clone, per-call table preparation or secondary valence rule.
    let symbol = symbol.as_ref();
    let number = periodic_table::atomic_number_from_symbol(symbol).ok_or_else(|| {
        SymbolValenceLookupError::SymbolNotFound {
            symbol: cosmolkit_model::PropertyText::from_bytes(symbol),
        }
    })?;
    Ok(rdkit_default_valence(number)?)
}

pub fn rdkit_default_valence_from_c_symbol(
    symbol: &std::ffi::CStr,
) -> Result<i32, SymbolValenceLookupError> {
    // RDKit❗✔️:   int getDefaultValence(const char *elementSymbol) const {
    // RDKit❗✔️:     return getDefaultValence(std::string(elementSymbol));
    // RDKit❗✔️:   }
    // Valid source char-pointer bytes stop at first NUL during std::string
    // construction. The complete string owner retains PRECONDITION order,
    // source signed first valence and structural lookup errors unchanged.
    // Borrowed CStr bytes avoid a temporary string copy and second name table.
    rdkit_default_valence_from_symbol(symbol.to_bytes())
}

pub fn rdkit_valence_list_from_symbol(
    symbol: impl AsRef<[u8]>,
) -> Result<&'static [i32], SymbolValenceLookupError> {
    // RDKit❗✔️:   const INT_VECT &getValenceList(const std::string &elementSymbol) const {
    // RDKit❗✔️:     PRECONDITION(byname.count(elementSymbol),
    // RDKit❗✔️:                  "Element '" + elementSymbol + "' not found");
    // RDKit❗✔️:     return getValenceList(byname.find(elementSymbol)->second);
    // RDKit❗✔️:   }
    // Native name presence PRECONDITION precedes the numeric getter. Preserve
    // the actual complete signed list, source order and terminal -1; a missing
    // name is an error, never a fabricated unrestricted or empty valence list.
    // Reuse both canonical owners: immutable byname projection and required
    // numeric slice. One shared lookup replaces source count/find; return the
    // same stable borrowed table storage without a per-call copy or allocation.
    let symbol = symbol.as_ref();
    let number = periodic_table::atomic_number_from_symbol(symbol).ok_or_else(|| {
        SymbolValenceLookupError::SymbolNotFound {
            symbol: cosmolkit_model::PropertyText::from_bytes(symbol),
        }
    })?;
    Ok(required_valence_list(number)?)
}

pub fn rdkit_valence_list_from_c_symbol(
    symbol: &std::ffi::CStr,
) -> Result<&'static [i32], SymbolValenceLookupError> {
    // RDKit❗✔️:   const INT_VECT &getValenceList(const char *elementSymbol) const {
    // RDKit❗✔️:     return getValenceList(std::string(elementSymbol));
    // RDKit❗✔️:   }
    // Source char-pointer construction ends at the first NUL. Delegate the
    // exact bytes to the unique complete string/list owner, retaining the
    // precondition, signed entries and borrowed table lifetime independently
    // of the temporary symbol. No copied list/string or secondary lookup.
    rdkit_valence_list_from_symbol(symbol.to_bytes())
}

fn invalid_valence_with_message<A: SourceValenceAtom>(
    atom: &A,
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

fn calculate_explicit_valence<A: SourceValenceAtom>(
    atoms: &[A],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
    check_it: bool,
) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION calculateExplicitValence complete source
    // RDKit✔️✔️: int calculateExplicitValence(const Atom &atom, bool strict, bool checkIt) {
    // RDKit✔️✔️:   // FIX: contributions of bonds to valence are being done at best
    // RDKit✔️✔️:   // approximately
    // RDKit✔️✔️:   double accum = 0;
    // RDKit✔️❌:   for (const auto bnd : atom.getOwningMol().atomBonds(&atom)) {
    // RDKit✔️✔️:     accum += bnd->getValenceContrib(&atom);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   accum += atom.getNumExplicitHs();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const auto &ovalens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(atom.getAtomicNum());
    // RDKit✔️✔️:   // if we start with an atom that doesn't have specified valences, we stick
    // RDKit✔️✔️:   // with that. otherwise we will use the effective valence
    // RDKit✔️✔️:   unsigned int effectiveAtomicNum = atom.getAtomicNum();
    // RDKit✔️✔️:   if (ovalens.size() > 1 || ovalens[0] != -1) {
    // RDKit✔️✔️:     effectiveAtomicNum = getEffectiveAtomicNum(atom, checkIt);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int dv =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getDefaultValence(effectiveAtomicNum);
    // RDKit✔️✔️:   const auto &valens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(effectiveAtomicNum);
    // RDKit✔️❌:   if (accum > dv && isAromaticAtom(atom)) {
    // RDKit✔️✔️:     // this needs some explanation : if the atom is aromatic and
    // RDKit✔️✔️:     // accum > dv we assume that no hydrogen can be added
    // RDKit✔️✔️:     // to this atom.  We set x = (v + chr) such that x is the
    // RDKit✔️✔️:     // closest possible integer to "accum" but less than
    // RDKit✔️✔️:     // "accum".
    // RDKit✔️✔️:     //
    // RDKit✔️✔️:     // "v" here is one of the allowed valences. For example:
    // RDKit✔️✔️:     //    sulfur here : O=c1ccs(=O)cc1
    // RDKit✔️✔️:     //    nitrogen here : c1cccn1C
    // RDKit✔️✔️:
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
    // RDKit✔️✔️:     // if we're within 1.5 of the allowed valence, go ahead and take it.
    // RDKit✔️✔️:     // this reflects things like the N in c1cccn1C, which starts with
    // RDKit✔️✔️:     // accum of 4, but which can be kekulized to C1=CC=CN1C, where
    // RDKit✔️✔️:     // the valence is 3 or the bridging N in c1ccn2cncc2c1, which starts
    // RDKit✔️✔️:     // with a valence of 4.5, but can be happily kekulized down to a valence
    // RDKit✔️✔️:     // of 3
    // RDKit✔️✔️:     if (accum - pval <= 1.5) {
    // RDKit✔️✔️:       accum = pval;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
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
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto res = static_cast<int>(std::round(accum));
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (strict || checkIt) {
    // RDKit✔️✔️:     int maxValence = valens.back();
    // RDKit✔️✔️:     int offset = 0;
    // RDKit✔️✔️:     // we have to include a special case here for negatively charged P, S, As,
    // RDKit✔️✔️:     // and Se, which all support "hypervalent" forms, but which can be
    // RDKit✔️✔️:     // isoelectronic to Cl/Ar or Br/Kr, which do not support hypervalent forms.
    // RDKit✔️✔️:     if (canBeHypervalent(atom, effectiveAtomicNum)) {
    // RDKit✔️✔️:       maxValence = ovalens.back();
    // RDKit✔️✔️:       offset -= atom.getFormalCharge();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // we have historically accepted two-coordinate [H-] as a valid atom. This
    // RDKit✔️✔️:     // is highly questionable, but changing it requires some thought. For now we
    // RDKit✔️✔️:     // will just keep accepting it
    // RDKit✔️✔️:     if (atom.getAtomicNum() == 1 && atom.getFormalCharge() == -1) {
    // RDKit✔️✔️:       maxValence = 2;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // maxValence == -1 signifies that we'll take anything at the high end
    // RDKit✔️✔️:     if (maxValence >= 0 && ovalens.back() >= 0 && (res + offset) > maxValence) {
    // RDKit✔️✔️:       // the explicit valence is greater than any
    // RDKit✔️✔️:       // allowed valence for the atoms
    // RDKit✔️✔️:
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
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION calculateExplicitValence complete source
    // RDKit declares dv unsigned: converting canonical -1 gives UINT32_MAX,
    // so a defined nonnegative int-range accum cannot enter aromatic lowering.
    // Preserve that conversion explicitly; pval is signed as in the source.
    // All other source branches, bond order, +0.1/round, original/max-valence
    // guards, hypervalent offset, hydride exception and strict/check error
    // paths remain unchanged; table data are borrowed from their sole owner.
    // O(degree + permitted valences), no successful-path heap allocation. The
    // detached incident-bond helper performs an additional validation prepass
    // compared with source-owned graph iteration; loop cost is marked worse.

    let atom = atom_from_parts(atoms, atom_id)?;
    let mut accum = 0.0;
    for bond in incident_bonds_from_parts(atoms.len(), bonds, adjacency, atom_id)? {
        accum += bond_valence_contrib(bond, atom_id)?;
    }
    accum += f64::from(atom.explicit_hydrogens());

    let ovalens = required_valence_list(atom.atomic_number())?;
    let mut effective_atomic_num = atom.atomic_number();
    if ovalens.len() > 1 || ovalens[0] != -1 {
        effective_atomic_num = get_effective_atomic_num_source(atom, check_it)?;
    }

    let default_valence = rdkit_default_valence(effective_atomic_num)? as u32;
    let valens = required_valence_list(effective_atomic_num)?;

    if accum > f64::from(default_valence)
        && is_aromatic_atom_from_parts(atoms, bonds, adjacency, atom_id)?
    {
        let mut pval = default_valence as i32;
        for &valence in valens {
            if valence == -1 || f64::from(valence) > accum {
                break;
            }
            pval = valence;
        }
        if accum - f64::from(pval) <= 1.5 {
            accum = f64::from(pval);
        }
    }

    accum += 0.1;
    let result = accum.round() as i32;

    if strict || check_it {
        let mut max_valence = *valens.last().expect("valence list is nonempty");
        let mut offset = 0;
        if can_be_hypervalent_source(atom, effective_atomic_num) {
            max_valence = *ovalens.last().expect("valence list is nonempty");
            offset -= i32::from(atom.formal_charge());
        }
        if atom.atomic_number() == 1 && atom.formal_charge() == -1 {
            max_valence = 2;
        }
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
    Ok(result)
}

fn calculate_cached_explicit_valence<A: SourceValenceAtom>(
    atoms: &[A],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION Atom::calcExplicitValence complete source
    // RDKit✔️✔️: int Atom::calcExplicitValence(bool strict) {
    // RDKit✔️✔️:   bool checkIt = false;
    // RDKit✔️✔️:   d_explicitValence = calculateExplicitValence(*this, strict, checkIt);
    // RDKit✔️✔️:   return d_explicitValence;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:   std::int8_t d_implicitValence, d_explicitValence;
    // END RDKIT CPP FUNCTION Atom::calcExplicitValence complete source
    // The pinned C++20 signed8 store precedes promotion of the returned cache.
    // Reuse this sole detached cache projection; scalar kernel returns int.
    // The fixed cast performs no allocation, graph copy, or extra traversal.
    Ok(i32::from(
        calculate_explicit_valence(atoms, bonds, adjacency, atom_id, strict, false)? as i8,
    ))
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

fn calculate_implicit_valence<A: SourceValenceAtom>(
    atoms: &[A],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    explicit_valence: i32,
    strict: bool,
    check_it: bool,
    source_complex_bonds: Option<&[bool]>,
) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION calculateImplicitValence complete source
    // RDKit✔️✔️: int calculateImplicitValence(const Atom &atom, bool strict, bool checkIt) {
    // RDKit✔️✔️:   if (atom.df_noImplicit) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto explicitValence = atom.d_explicitValence;
    // RDKit✔️✔️:   if (explicitValence == -1) {
    // RDKit✔️✔️:     explicitValence = calculateExplicitValence(atom, strict, checkIt);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // special cases
    // RDKit✔️✔️:   auto atomicNum = atom.d_atomicNum;
    // RDKit✔️✔️:   if (atomicNum == 0) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️❌:   for (const auto bnd : atom.getOwningMol().atomBonds(&atom)) {
    // RDKit✔️✔️:     if (QueryOps::hasComplexBondTypeQuery(*bnd)) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
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
    // RDKit✔️✔️:   int explicitPlusRadV = atom.d_explicitValence + atom.d_numRadicalElectrons;
    // RDKit✔️✔️:
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
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // this is basically the difference between the allowed valence of
    // RDKit✔️✔️:   // the atom and the explicit valence already specified - tells how
    // RDKit✔️✔️:   // many Hs to add
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // The d-block and f-block of the periodic table (i.e. transition metals,
    // RDKit✔️✔️:   // lanthanoids and actinoids) have no default valence.
    // RDKit✔️✔️:   int dv = PeriodicTable::getTable()->getDefaultValence(effectiveAtomicNum);
    // RDKit✔️✔️:   if (dv == -1) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
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
    // RDKit✔️✔️:   const auto &valens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(effectiveAtomicNum);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int res = 0;
    // RDKit✔️✔️:   // if we have an aromatic case treat it differently
    // RDKit✔️❌:   if (isAromaticAtom(atom)) {
    // RDKit✔️✔️:     if (explicitPlusRadV <= dv) {
    // RDKit✔️✔️:       res = dv - explicitPlusRadV;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       // As we assume when finding the explicitPlusRadValence if we are
    // RDKit✔️✔️:       // aromatic we should not be adding any hydrogen and already
    // RDKit✔️✔️:       // be at an accepted valence state,
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // FIX: this is just ERROR checking and probably moot - the
    // RDKit✔️✔️:       // explicitPlusRadValence function called above should assure us that
    // RDKit✔️✔️:       // we satisfy one of the accepted valence states for the
    // RDKit✔️✔️:       // atom. The only diff I can think of is in the way we handle
    // RDKit✔️✔️:       // formal charge here vs the explicit valence function.
    // RDKit✔️✔️:       bool satis = false;
    // RDKit✔️✔️:       for (auto vi = valens.begin(); vi != valens.end() && *vi > 0; ++vi) {
    // RDKit✔️✔️:         if (explicitPlusRadV == *vi) {
    // RDKit✔️✔️:           satis = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
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
    // RDKit✔️✔️:       res = 0;
    // RDKit✔️✔️:     }
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
    // RDKit✔️✔️:     if (res < 0) {
    // RDKit✔️✔️:       if ((strict || checkIt) && valens.back() != -1 && ovalens.back() > 0) {
    // RDKit✔️✔️:         // this means that the explicit valence is greater than any
    // RDKit✔️✔️:         // allowed valence for the atoms
    // RDKit✔️✔️:         if (strict) {
    // RDKit✔️✔️:           // raise an error
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
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION calculateImplicitValence complete source
    // The source local `auto` inherits int8 from d_explicitValence. Sentinel
    // recomputation narrows that local and does not update the original field;
    // explicitPlusRadV deliberately reads the original field plus radicals.
    // noImplicit is checked before any recomputation; the default valence here
    // is signed int, unlike calculateExplicitValence's unsigned default.
    // Existing Bond.query identity is inspected by the canonical source
    // complex-type-query classifier, including source sibling-order state.
    // O(degree + allowed valences), no success-path allocation/graph clone.
    // Query-bond iteration and the reached aromatic helper retain their extra
    // detached validation passes versus source-owned iteration, marked worse.

    let atom = atom_from_parts(atoms, atom_id)?;
    if atom.no_implicit() {
        return Ok(0);
    }
    let stored_explicit_valence = explicit_valence;
    let explicit_valence = if stored_explicit_valence == -1 {
        i32::from(
            calculate_explicit_valence(atoms, bonds, adjacency, atom_id, strict, check_it)? as i8,
        )
    } else {
        stored_explicit_valence
    };
    let atomic_num = atom.atomic_number();
    if atomic_num == 0 {
        return Ok(0);
    }
    for bond in incident_bonds_from_parts(atoms.len(), bonds, adjacency, atom_id)? {
        if source_complex_bonds.map_or_else(
            || crate::query_ops::bond_has_complex_type_query(bond),
            |mask| mask[bond.id().index()],
        ) {
            return Ok(0);
        }
    }
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

    let mut explicit_plus_rad_v = stored_explicit_valence + i32::from(atom.radical_electrons());

    let ovalens = required_valence_list(atomic_num)?;
    let mut effective_atomic_num = atomic_num;
    if ovalens.len() > 1 || ovalens[0] != -1 {
        effective_atomic_num = get_effective_atomic_num_source(atom, check_it)?;
    }
    if effective_atomic_num == 0 {
        return Ok(0);
    }

    let default_valence = rdkit_default_valence(effective_atomic_num)?;
    if default_valence == -1 {
        return Ok(0);
    }

    if can_be_hypervalent_source(atom, effective_atomic_num) {
        effective_atomic_num = atomic_num;
        explicit_plus_rad_v -= i32::from(atom.formal_charge());
    }
    let valens = required_valence_list(effective_atomic_num)?;

    let result;
    if is_aromatic_atom_from_parts(atoms, bonds, adjacency, atom_id)? {
        if explicit_plus_rad_v <= default_valence {
            result = default_valence - explicit_plus_rad_v;
        } else {
            let satisfied = valens
                .iter()
                .take_while(|&&valence| valence > 0)
                .any(|&valence| explicit_plus_rad_v == valence);
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
            result = 0;
        }
    } else {
        let mut candidate = -1;
        for &valence in valens.iter().take_while(|&&valence| valence >= 0) {
            if explicit_plus_rad_v <= valence {
                candidate = valence - explicit_plus_rad_v;
                break;
            }
        }
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
    Ok(result)
}

fn calculate_cached_implicit_valence<A: SourceValenceAtom>(
    atoms: &[A],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    explicit_valence: &mut i32,
    strict: bool,
    source_complex_bonds: Option<&[bool]>,
) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION Atom::calcImplicitValence complete source
    // RDKit✔️✔️: int Atom::calcImplicitValence(bool strict) {
    // RDKit✔️✔️:   if (d_explicitValence == -1) {
    // RDKit✔️✔️:     calcExplicitValence(strict);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool checkIt = false;
    // RDKit✔️✔️:   d_implicitValence = calculateImplicitValence(*this, strict, checkIt);
    // RDKit✔️✔️:   return d_implicitValence;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:   std::int8_t d_implicitValence, d_explicitValence;
    // END RDKIT CPP FUNCTION Atom::calcImplicitValence complete source
    // This mutable scalar is a detached projection of the explicit cache.
    // Recompute its sentinel before entering the kernel, including noImplicit.
    // The C++20 signed8 implicit cache store precedes return-value promotion.
    // Successful source state transitions remain ordered; errors propagate
    // before the implicit store. No graph clone or additional allocation.
    if *explicit_valence == -1 {
        *explicit_valence =
            calculate_cached_explicit_valence(atoms, bonds, adjacency, atom_id, strict)?;
    }
    Ok(i32::from(calculate_implicit_valence(
        atoms,
        bonds,
        adjacency,
        atom_id,
        *explicit_valence,
        strict,
        false,
        source_complex_bonds,
    )? as i8))
}

fn calculate_cached_valence_state<A: SourceValenceAtom>(
    atoms: &[A],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
    source_complex_bonds: Option<&[bool]>,
) -> Result<(i32, i32), ValenceError> {
    let mut facts = SourceAtomValenceFacts::UNINITIALIZED;
    calculate_cached_valence_state_with_facts(
        atoms,
        bonds,
        adjacency,
        atom_id,
        &mut facts,
        strict,
        source_complex_bonds,
    )
}

fn calculate_cached_valence_state_with_facts<A: SourceValenceAtom>(
    atoms: &[A],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    facts: &mut SourceAtomValenceFacts,
    strict: bool,
    source_complex_bonds: Option<&[bool]>,
) -> Result<(i32, i32), ValenceError> {
    // BEGIN RDKIT CPP FUNCTION Atom::updatePropertyCache complete source
    // RDKit✔️✔️: void Atom::updatePropertyCache(bool strict) {
    // RDKit✔️✔️:   calcExplicitValence(strict);
    // RDKit✔️✔️:   calcImplicitValence(strict);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::updatePropertyCache complete source
    // The detached facts preserve each completed native signed8 store, including
    // the explicit store if the following implicit calculation returns an error.
    // Both readonly assignment and mutable detached consumers use this one kernel.
    let mut explicit = calculate_cached_explicit_valence(atoms, bonds, adjacency, atom_id, strict)?;
    facts.explicit_valence = explicit as i8;
    let implicit = calculate_cached_implicit_valence(
        atoms,
        bonds,
        adjacency,
        atom_id,
        &mut explicit,
        strict,
        source_complex_bonds,
    );
    facts.explicit_valence = explicit as i8;
    let implicit = implicit?;
    facts.implicit_valence = implicit as i8;
    Ok((explicit, implicit))
}

pub(crate) fn update_source_atom_cache(
    topology: &mut TopologyBlock,
    atom_id: AtomId,
    strict: bool,
) -> Result<(), ValenceError> {
    let mut facts = atom_from_parts(&topology.atoms, atom_id)?.source_valence_facts();
    let result = calculate_cached_valence_state_with_facts(
        &topology.atoms,
        &topology.bonds,
        &topology.adjacency,
        atom_id,
        &mut facts,
        strict,
        None,
    );
    // Source effect values only: this helper has no live runtime cache authority.
    topology.atoms[atom_id.index()].set_source_valence_facts(facts);
    result.map(|_| ())
}

pub(crate) fn source_atom_needs_cache_update(atom: &Atom) -> bool {
    // BEGIN RDKIT CPP FUNCTION Atom::needsUpdatePropertyCache complete source
    // RDKit✔️✔️: bool Atom::needsUpdatePropertyCache() const {
    // RDKit✔️✔️:   return !(this->d_explicitValence >= 0 &&
    // RDKit✔️✔️:            (this->df_noImplicit || this->d_implicitValence >= 0));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::needsUpdatePropertyCache complete source
    source_atom_needs_cache_update_with_facts(atom, atom.source_valence_facts())
}

pub(crate) fn source_atom_needs_cache_update_with_facts(
    atom: &Atom,
    facts: SourceAtomValenceFacts,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION Atom::needsUpdatePropertyCache complete source
    // RDKit✔️✔️: bool Atom::needsUpdatePropertyCache() const {
    // RDKit✔️✔️:   return !(this->d_explicitValence >= 0 &&
    // RDKit✔️✔️:            (this->df_noImplicit || this->d_implicitValence >= 0));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::needsUpdatePropertyCache complete source
    // The same predicate serves attached facts and the detached scalar overlay.
    !(facts.explicit_valence >= 0 && (atom.no_implicit() || facts.implicit_valence >= 0))
}

pub(crate) fn source_atom_implicit_hydrogens(atom: &Atom) -> Result<u32, &'static str> {
    // BEGIN RDKIT CPP FUNCTION Atom::getNumImplicitHs complete source
    // RDKit✔️✔️: unsigned int Atom::getNumImplicitHs() const {
    // RDKit✔️✔️:   if (df_noImplicit) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   PRECONDITION(d_implicitValence > -1,
    // RDKit✔️✔️:                "getNumImplicitHs() called without preceding call to "
    // RDKit✔️✔️:                "calcImplicitValence()");
    // RDKit✔️✔️:   return getValence(ValenceType::IMPLICIT);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::getNumImplicitHs complete source
    // BEGIN RDKIT CPP FUNCTION Atom::getValence complete source
    // RDKit✔️✔️: unsigned int Atom::getValence(ValenceType which) const {
    // RDKit✔️✔️:   if (!dp_mol) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       (which == ValenceType::IMPLICIT || d_explicitValence > -1),
    // RDKit✔️✔️:       "getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()");
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       (which == ValenceType::EXPLICIT || df_noImplicit ||
    // RDKit✔️✔️:        d_implicitValence > -1),
    // RDKit✔️✔️:       "getValence(ValenceType::IMPLICIT) called without call to calcImplicitValence()");
    // RDKit✔️✔️:   if (which == ValenceType::EXPLICIT) {
    // RDKit✔️✔️:     return d_explicitValence;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return df_noImplicit ? 0 : d_implicitValence;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::getValence complete source
    // This input is an atom in an attached detached topology. The orphan
    // !dp_mol return branch is outside that structural input boundary.
    if atom.no_implicit() {
        return Ok(0);
    }
    let implicit = atom.source_valence_facts().implicit_valence;
    if implicit < 0 {
        return Err("getNumImplicitHs() called without preceding call to calcImplicitValence()");
    }
    Ok(implicit as u32)
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
        None,
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
    // The optional query supplies a freshly calculated source cache when
    // none is given; the lower-level kernel separately accepts original -1.
    let explicit_valence = match explicit_valence {
        Some(value) => value,
        None => calculate_cached_explicit_valence(
            &topology.atoms,
            &topology.bonds,
            &topology.adjacency,
            atom_id,
            strict,
        )?,
    };
    calculate_implicit_valence_for_topology(topology, atom_id, explicit_valence, strict, false)
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
        None,
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
    // BEGIN RDKIT CPP FUNCTION ROMol::updatePropertyCache complete source
    // RDKit✔️✔️: void ROMol::updatePropertyCache(bool strict) {
    // RDKit✔️❌:   for (auto atom : atoms()) {
    // RDKit✔️❌:     atom->updatePropertyCache(strict);
    // RDKit✔️✔️:   }
    // RDKit✔️🔝:   for (auto bond : bonds()) {
    // RDKit✔️🔝:     bond->updatePropertyCache(strict);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Reuse the complete source Atom update helper once per stable atom row.
    // The existing detached result needs two O(V) output allocations instead
    // of source in-place cache writes; this is a known cost difference.
    // Pinned Bond::updatePropertyCache below is a nonvirtual no-op. Omitting
    // its O(E) traversal improves that loop without changing any state/error.
    let mut explicit_values = Vec::with_capacity(atoms.len());
    let mut implicit_hydrogens = Vec::with_capacity(atoms.len());
    for atom in atoms {
        let (explicit, implicit) =
            calculate_cached_valence_state(atoms, bonds, adjacency, atom.id(), strict, None)?;
        explicit_values.push(explicit);
        implicit_hydrogens.push(implicit);
    }
    // RDKit✔️✔️: void updatePropertyCache(bool strict = true) { (void)strict; }
    // END RDKIT CPP FUNCTION ROMol::updatePropertyCache complete source
    Ok(ValenceAssignment {
        explicit_valence: explicit_values,
        implicit_hydrogens,
    })
}

pub(crate) fn source_calc_implicit_cache_row(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    explicit: &mut i32,
    implicit: &mut i32,
    strict: bool,
) -> Result<(), ValenceError> {
    // RDKit✔️✔️: Atom::calcImplicitValence full source is anchored in the
    // calculate_cached_implicit_valence kernel directly called here.
    // Detached ValenceAssignment is this owner's source cache projection.
    // Normalize signed8 E before the same existing scalar kernel; I commits
    // only on success, while E's completed write survives an implicit error.
    *explicit = i32::from(*explicit as i8);
    let next = calculate_cached_implicit_valence(
        atoms, bonds, adjacency, atom_id, explicit, strict, None,
    )?;
    *implicit = next;
    Ok(())
}

/// Assign both valence fields for one atom from borrowed detached parts.
pub fn assign_valence_state_for_atom_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
) -> Result<(i32, i32), ValenceError> {
    calculate_cached_valence_state(atoms, bonds, adjacency, atom_id, strict, None)
}

/// Assign one atom's explicit valence from borrowed detached parts.
pub fn assign_explicit_valence_for_atom_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
) -> Result<i32, ValenceError> {
    calculate_cached_explicit_valence(atoms, bonds, adjacency, atom_id, strict)
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
    let mut explicit_cache = explicit_valence;
    calculate_cached_implicit_valence(
        atoms,
        bonds,
        adjacency,
        atom_id,
        &mut explicit_cache,
        strict,
        None,
    )
    .map(|value| i32::from(value as i8))
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
    Ok(calculate_implicit_valence(atoms, bonds, adjacency, id, explicit, false, true, None)? == -1)
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

/// Source numPiElectrons over a detached topology and exact selected i8 cache.
/// Aromatic/Sp3 short circuits do not access or initialize the valence cache.
/// Canonical owned-topology boundary; unowned C++ Atom state is not modeled.
/// Fixed source cases pass; broader source-native parity and p1 audit remain pending.
/// Existing checked incident traversal uses two O(degree) passes, no allocations
/// or whole-topology cloning; source cache access stays lazy.
pub fn num_pi_electrons_for_topology(
    topology: &TopologyBlock,
    atom_id: AtomId,
    explicit_valence: Option<i8>,
) -> Result<u32, ValenceError> {
    // RDKit❗✔️: unsigned int numPiElectrons(const Atom &atom) {
    // RDKit❗✔️:   unsigned int res = 0;
    // RDKit❗✔️:   if (atom.getIsAromatic()) {
    // RDKit❗✔️:     res = 1;
    // RDKit❗✔️:   } else if (atom.getHybridization() != Atom::SP3) {
    // RDKit❗✔️:     auto val =
    // RDKit❗✔️:         static_cast<unsigned int>(atom.getValence(Atom::ValenceType::EXPLICIT));
    // RDKit❗✔️:     unsigned int physical_bonds = atom.getNumExplicitHs();
    // RDKit❗✔️:     const auto &mol = atom.getOwningMol();
    // RDKit❗✔️:     for (const auto bond : mol.atomBonds(&atom)) {
    // RDKit❗✔️:       if (bond->getValenceContrib(&atom) != 0.0) {
    // RDKit❗✔️:         ++physical_bonds;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     CHECK_INVARIANT(val >= physical_bonds,
    // RDKit❗✔️:                     "explicit valence exceeds atom degree");
    // RDKit❗✔️:     res = val - physical_bonds;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit❗✔️: }  // namespace RDKit
    // RDKit❗✔️: unsigned int Atom::getValence(ValenceType which) const {
    // RDKit❗✔️:   if (!dp_mol) {
    // RDKit❗✔️:     return 0;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   PRECONDITION(
    // RDKit❗✔️:       (which == ValenceType::IMPLICIT || d_explicitValence > -1),
    // RDKit❗✔️:       "getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()");
    // RDKit❗✔️:   PRECONDITION(
    // RDKit❗✔️:       (which == ValenceType::EXPLICIT || df_noImplicit ||
    // RDKit❗✔️:        d_implicitValence > -1),
    // RDKit❗✔️:       "getValence(ValenceType::IMPLICIT) called without call to calcImplicitValence()");
    // RDKit❗✔️:   if (which == ValenceType::EXPLICIT) {
    // RDKit❗✔️:     return d_explicitValence;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     return df_noImplicit ? 0 : d_implicitValence;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️:
    let atom = atom_from_parts(&topology.atoms, atom_id)?;
    if atom.is_aromatic() {
        return Ok(1);
    }
    if atom.hybridization() == cosmolkit_types::Hybridization::Sp3 {
        return Ok(0);
    }
    // Source cache is signed8. None and every negative cached value fail the
    // explicit getter precondition. Caller i32 assignments are not narrowed.
    let value = explicit_valence
        .filter(|value| *value >= 0)
        .ok_or(ValenceError::PiElectronExplicitValenceCacheNotInitialized { atom: atom_id })?;
    let valence = value as u32;
    let mut physical_bonds = u32::from(atom.explicit_hydrogens());
    for bond in incident_bonds_from_parts(
        topology.atoms.len(),
        &topology.bonds,
        &topology.adjacency,
        atom_id,
    )? {
        if bond_valence_contrib(bond, atom_id)? != 0.0 {
            physical_bonds = physical_bonds.wrapping_add(1);
        }
    }
    if valence < physical_bonds {
        return Err(ValenceError::PiElectronInvariant {
            atom: atom_id,
            explicit_valence: valence,
            physical_bonds,
        });
    }
    Ok(valence - physical_bonds)
}

/// Read all atom metadata with one existing valence-owner calculation.
/// This does not install a cache or grant callers any runtime authority.
pub fn atom_metadata(topology: &TopologyBlock) -> Result<Vec<AtomMetadata>, ValenceError> {
    let assignment = assign_valence(topology, &ValenceParams::default())?;
    atom_metadata_from_assignment(topology, Some(&assignment))
}

/// Project source atom getters from an existing cache without recalculating it.
/// The assignment is borrowed; absent or uninitialized entries remain typed errors.
pub fn atom_metadata_from_assignment(
    topology: &TopologyBlock,
    assignment: Option<&ValenceAssignment>,
) -> Result<Vec<AtomMetadata>, ValenceError> {
    // RDKit✔️❌: unsigned int Atom::getDegree() const {
    // RDKit✔️❌:   return dp_mol ? getOwningMol().getAtomDegree(this) : 0;
    // RDKit✔️❌: }
    // RDKit✔️❌:
    // RDKit✔️❌: unsigned int Atom::getValence(ValenceType which) const {
    // RDKit✔️❌:   if (!dp_mol) {
    // RDKit✔️❌:     return 0;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   PRECONDITION(
    // RDKit✔️❌:       (which == ValenceType::IMPLICIT || d_explicitValence > -1),
    // RDKit✔️❌:       "getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()");
    // RDKit✔️❌:   PRECONDITION(
    // RDKit✔️❌:       (which == ValenceType::EXPLICIT || df_noImplicit ||
    // RDKit✔️❌:        d_implicitValence > -1),
    // RDKit✔️❌:       "getValence(ValenceType::IMPLICIT) called without call to calcImplicitValence()");
    // RDKit✔️❌:   if (which == ValenceType::EXPLICIT) {
    // RDKit✔️❌:     return d_explicitValence;
    // RDKit✔️❌:   } else {
    // RDKit✔️❌:     return df_noImplicit ? 0 : d_implicitValence;
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // RDKit✔️❌:
    // RDKit✔️❌: unsigned int Atom::getTotalValence() const {
    // RDKit✔️❌:   return getValence(ValenceType::EXPLICIT) + getValence(ValenceType::IMPLICIT);
    // RDKit✔️❌: }
    // RDKit✔️❌: std::int8_t d_implicitValence, d_explicitValence;
    // RDKit✔️❌: int Atom::calcExplicitValence(bool strict) {
    // RDKit✔️❌:   bool checkIt = false;
    // RDKit✔️❌:   d_explicitValence = calculateExplicitValence(*this, strict, checkIt);
    // RDKit✔️❌:   return d_explicitValence;
    // RDKit✔️❌: }
    // Caller-supplied detached assignments can retain calculation-width i32 values. Reproduce
    // the source's int8 storage conversion at this cached-getter boundary,
    // before applying its > -1 precondition; never recalculate strict valence.
    // Detached validation adds O(V+E) work over O(V) source cached getters,
    // so the complexity marker records that extra boundary cost.
    // Owned rows require one O(V) allocation. Validating detached topology once
    // costs O(V+E); each getter then performs O(1) cache/adjacency access.
    validate_topology(topology)?;
    let mut rows = Vec::with_capacity(topology.atoms.len());
    for atom in &topology.atoms {
        let id = atom.id();
        let explicit_valence = cached_explicit_valence(atom, assignment)?;
        let total_hydrogens = crate::hcount::total_hydrogen_count_from_validated(
            topology,
            assignment.expect("explicit cache checked"),
            id,
            false,
        )? as i32;
        let implicit_hydrogens = total_hydrogens - i32::from(atom.explicit_hydrogens());
        rows.push(AtomMetadata {
            degree: topology.adjacency.neighbors_of(id.index()).len(),
            explicit_valence,
            implicit_hydrogens,
            total_hydrogens,
            total_valence: explicit_valence + implicit_hydrogens,
        });
    }
    Ok(rows)
}

/// Classify the source signed cache fields without preparing replacements.
#[doc(hidden)]
pub fn atom_valence_cache_needs_update(atom: &Atom, explicit: i32, implicit: i32) -> bool {
    // RDKit✔️✔️: bool Atom::needsUpdatePropertyCache() const {
    // RDKit✔️✔️:   return !(this->d_explicitValence >= 0 &&
    // RDKit✔️✔️:            (this->df_noImplicit || this->d_implicitValence >= 0));
    // RDKit✔️✔️: }
    // RDKit✔️✔️:   std::int8_t d_implicitValence, d_explicitValence;
    // Behavior: classify signed source storage, including NoImplicit's I bypass.
    // Complexity: O(1), no graph traversal, allocation or cached-field mutation.
    let explicit = explicit as i8;
    let implicit = implicit as i8;
    !(explicit >= 0 && (atom.no_implicit() || implicit >= 0))
}

pub(crate) fn cached_explicit_valence(
    atom: &Atom,
    assignment: Option<&ValenceAssignment>,
) -> Result<i32, ValenceError> {
    // RDKit✔️✔️: unsigned int Atom::getValence(ValenceType which) const {
    // RDKit✔️✔️:   if (!dp_mol) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       (which == ValenceType::IMPLICIT || d_explicitValence > -1),
    // RDKit✔️✔️:       "getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()");
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       (which == ValenceType::EXPLICIT || df_noImplicit ||
    // RDKit✔️✔️:        d_implicitValence > -1),
    // RDKit✔️✔️:       "getValence(ValenceType::IMPLICIT) called without call to calcImplicitValence()");
    // RDKit✔️✔️:   if (which == ValenceType::EXPLICIT) {
    // RDKit✔️✔️:     return d_explicitValence;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return df_noImplicit ? 0 : d_implicitValence;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: unsigned int Atom::getTotalValence() const {
    // RDKit✔️✔️:   return getValence(ValenceType::EXPLICIT) + getValence(ValenceType::IMPLICIT);
    // RDKit✔️✔️: }
    assignment
        .and_then(|a| a.explicit_valence.get(atom.id().index()))
        .copied()
        .map(|value| i32::from(value as i8))
        .filter(|value| *value >= 0)
        .ok_or(ValenceError::ExplicitValenceCacheNotInitialized { atom: atom.id() })
}

pub(crate) fn cached_total_valence(
    atom: &Atom,
    assignment: &ValenceAssignment,
) -> Result<i32, ValenceError> {
    // RDKit✔️✔️: unsigned int Atom::getTotalValence() const {
    // RDKit✔️✔️:   return getValence(ValenceType::EXPLICIT) + getValence(ValenceType::IMPLICIT);
    // RDKit✔️✔️: }
    // Both source getters return initialized nonnegative int8 values; their
    // sum is bounded by 254. Reuse the unique implicit getter, including its
    // noImplicit early return, with no allocation or graph traversal.
    let explicit = cached_explicit_valence(atom, Some(assignment))?;
    let implicit = crate::hcount::implicit_hydrogen_count(atom, assignment)? as i32;
    Ok(explicit + implicit)
}

#[cfg(test)]
mod source_cached_field_initialization_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, Element};
    #[test]
    fn source_initialization_observes_signed_width_and_no_implicit() {
        let rows = [
            (0, 0, false, false),
            (4, 0, false, false),
            (0, -1, false, true),
            (0, -1, true, false),
            (-1, 0, true, true),
            (128, 0, false, true),
            (255, 0, true, true),
            (256, 0, false, false),
            (-256, 0, false, false),
            (-129, 0, false, false),
            (0, 128, false, true),
            (0, 255, false, true),
            (0, 128, true, false),
            (0, 256, false, false),
            (0, -256, false, false),
            (127, 127, false, false),
        ];
        for (explicit, implicit, no_implicit, expected) in rows {
            let atom = Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::C).with_no_implicit(no_implicit),
            );
            assert_eq!(
                atom_valence_cache_needs_update(&atom, explicit, implicit),
                expected,
                "{explicit}/{implicit}/{no_implicit}"
            );
        }
    }
}

#[cfg(test)]
mod complete_atomic_symbol_string_source_tests {
    use super::*;

    #[test]
    fn every_pinned_byname_row_including_dummy_and_legacy_aliases_is_exact() {
        // Mechanically prepared from pinned atomic_data.cpp's 121 symbol rows;
        // duplicate atomic numbers retain both source names, without alias inference.
        for (symbol, number) in [
            (b"*".as_slice(), 0),
            (b"H".as_slice(), 1),
            (b"He".as_slice(), 2),
            (b"Li".as_slice(), 3),
            (b"Be".as_slice(), 4),
            (b"B".as_slice(), 5),
            (b"C".as_slice(), 6),
            (b"N".as_slice(), 7),
            (b"O".as_slice(), 8),
            (b"F".as_slice(), 9),
            (b"Ne".as_slice(), 10),
            (b"Na".as_slice(), 11),
            (b"Mg".as_slice(), 12),
            (b"Al".as_slice(), 13),
            (b"Si".as_slice(), 14),
            (b"P".as_slice(), 15),
            (b"S".as_slice(), 16),
            (b"Cl".as_slice(), 17),
            (b"Ar".as_slice(), 18),
            (b"K".as_slice(), 19),
            (b"Ca".as_slice(), 20),
            (b"Sc".as_slice(), 21),
            (b"Ti".as_slice(), 22),
            (b"V".as_slice(), 23),
            (b"Cr".as_slice(), 24),
            (b"Mn".as_slice(), 25),
            (b"Fe".as_slice(), 26),
            (b"Co".as_slice(), 27),
            (b"Ni".as_slice(), 28),
            (b"Cu".as_slice(), 29),
            (b"Zn".as_slice(), 30),
            (b"Ga".as_slice(), 31),
            (b"Ge".as_slice(), 32),
            (b"As".as_slice(), 33),
            (b"Se".as_slice(), 34),
            (b"Br".as_slice(), 35),
            (b"Kr".as_slice(), 36),
            (b"Rb".as_slice(), 37),
            (b"Sr".as_slice(), 38),
            (b"Y".as_slice(), 39),
            (b"Zr".as_slice(), 40),
            (b"Nb".as_slice(), 41),
            (b"Mo".as_slice(), 42),
            (b"Tc".as_slice(), 43),
            (b"Ru".as_slice(), 44),
            (b"Rh".as_slice(), 45),
            (b"Pd".as_slice(), 46),
            (b"Ag".as_slice(), 47),
            (b"Cd".as_slice(), 48),
            (b"In".as_slice(), 49),
            (b"Sn".as_slice(), 50),
            (b"Sb".as_slice(), 51),
            (b"Te".as_slice(), 52),
            (b"I".as_slice(), 53),
            (b"Xe".as_slice(), 54),
            (b"Cs".as_slice(), 55),
            (b"Ba".as_slice(), 56),
            (b"La".as_slice(), 57),
            (b"Ce".as_slice(), 58),
            (b"Pr".as_slice(), 59),
            (b"Nd".as_slice(), 60),
            (b"Pm".as_slice(), 61),
            (b"Sm".as_slice(), 62),
            (b"Eu".as_slice(), 63),
            (b"Gd".as_slice(), 64),
            (b"Tb".as_slice(), 65),
            (b"Dy".as_slice(), 66),
            (b"Ho".as_slice(), 67),
            (b"Er".as_slice(), 68),
            (b"Tm".as_slice(), 69),
            (b"Yb".as_slice(), 70),
            (b"Lu".as_slice(), 71),
            (b"Hf".as_slice(), 72),
            (b"Ta".as_slice(), 73),
            (b"W".as_slice(), 74),
            (b"Re".as_slice(), 75),
            (b"Os".as_slice(), 76),
            (b"Ir".as_slice(), 77),
            (b"Pt".as_slice(), 78),
            (b"Au".as_slice(), 79),
            (b"Hg".as_slice(), 80),
            (b"Tl".as_slice(), 81),
            (b"Pb".as_slice(), 82),
            (b"Bi".as_slice(), 83),
            (b"Po".as_slice(), 84),
            (b"At".as_slice(), 85),
            (b"Rn".as_slice(), 86),
            (b"Fr".as_slice(), 87),
            (b"Ra".as_slice(), 88),
            (b"Ac".as_slice(), 89),
            (b"Th".as_slice(), 90),
            (b"Pa".as_slice(), 91),
            (b"U".as_slice(), 92),
            (b"Np".as_slice(), 93),
            (b"Pu".as_slice(), 94),
            (b"Am".as_slice(), 95),
            (b"Cm".as_slice(), 96),
            (b"Bk".as_slice(), 97),
            (b"Cf".as_slice(), 98),
            (b"Es".as_slice(), 99),
            (b"Fm".as_slice(), 100),
            (b"Md".as_slice(), 101),
            (b"No".as_slice(), 102),
            (b"Lr".as_slice(), 103),
            (b"Rf".as_slice(), 104),
            (b"Db".as_slice(), 105),
            (b"Sg".as_slice(), 106),
            (b"Bh".as_slice(), 107),
            (b"Hs".as_slice(), 108),
            (b"Mt".as_slice(), 109),
            (b"Ds".as_slice(), 110),
            (b"Rg".as_slice(), 111),
            (b"Cn".as_slice(), 112),
            (b"Nh".as_slice(), 113),
            (b"Uut".as_slice(), 113),
            (b"Fl".as_slice(), 114),
            (b"Mc".as_slice(), 115),
            (b"Uup".as_slice(), 115),
            (b"Lv".as_slice(), 116),
            (b"Ts".as_slice(), 117),
            (b"Og".as_slice(), 118),
        ] {
            assert_eq!(
                rdkit_atomic_number_from_symbol(symbol),
                Ok(number),
                "{symbol:?}"
            );
        }
    }

    #[test]
    fn missing_symbols_are_errors_with_exact_original_bytes() {
        for bytes in [
            b"".as_slice(),
            b"c",
            b"cl",
            b"CL",
            b" C",
            b"C ",
            b"C\t",
            b"Carbon",
            b"Uuo",
            b"D",
            b"T",
            b"C\0",
            b"C\0garbage",
            b"\xff",
            b"C\xff",
            b"\xc3\xa9",
        ] {
            let error = rdkit_atomic_number_from_symbol(bytes).unwrap_err();
            assert_eq!(error.symbol.as_bytes(), bytes, "{bytes:?}");
        }
    }

    #[test]
    fn string_nul_is_not_a_c_string_terminator_and_known_symbols_borrow_input() {
        for (symbol, number) in [
            ("C", 6),
            ("N", 7),
            ("O", 8),
            ("Cl", 17),
            ("Og", 118),
            ("*", 0),
            ("Uut", 113),
            ("Uup", 115),
        ] {
            assert_eq!(rdkit_atomic_number_from_symbol(symbol), Ok(number));
            let mut raw = symbol.as_bytes().to_vec();
            raw.push(0);
            assert_eq!(
                rdkit_atomic_number_from_symbol(&raw)
                    .unwrap_err()
                    .symbol
                    .as_bytes(),
                raw
            );
        }
    }
}

#[cfg(test)]
mod complete_atomic_symbol_c_source_tests {
    use super::*;
    use std::ffi::CStr;

    #[test]
    fn every_real_table_name_delegates_through_actual_c_string() {
        for (symbol, number) in [
            (b"*".as_slice(), 0),
            (b"H".as_slice(), 1),
            (b"He".as_slice(), 2),
            (b"Li".as_slice(), 3),
            (b"Be".as_slice(), 4),
            (b"B".as_slice(), 5),
            (b"C".as_slice(), 6),
            (b"N".as_slice(), 7),
            (b"O".as_slice(), 8),
            (b"F".as_slice(), 9),
            (b"Ne".as_slice(), 10),
            (b"Na".as_slice(), 11),
            (b"Mg".as_slice(), 12),
            (b"Al".as_slice(), 13),
            (b"Si".as_slice(), 14),
            (b"P".as_slice(), 15),
            (b"S".as_slice(), 16),
            (b"Cl".as_slice(), 17),
            (b"Ar".as_slice(), 18),
            (b"K".as_slice(), 19),
            (b"Ca".as_slice(), 20),
            (b"Sc".as_slice(), 21),
            (b"Ti".as_slice(), 22),
            (b"V".as_slice(), 23),
            (b"Cr".as_slice(), 24),
            (b"Mn".as_slice(), 25),
            (b"Fe".as_slice(), 26),
            (b"Co".as_slice(), 27),
            (b"Ni".as_slice(), 28),
            (b"Cu".as_slice(), 29),
            (b"Zn".as_slice(), 30),
            (b"Ga".as_slice(), 31),
            (b"Ge".as_slice(), 32),
            (b"As".as_slice(), 33),
            (b"Se".as_slice(), 34),
            (b"Br".as_slice(), 35),
            (b"Kr".as_slice(), 36),
            (b"Rb".as_slice(), 37),
            (b"Sr".as_slice(), 38),
            (b"Y".as_slice(), 39),
            (b"Zr".as_slice(), 40),
            (b"Nb".as_slice(), 41),
            (b"Mo".as_slice(), 42),
            (b"Tc".as_slice(), 43),
            (b"Ru".as_slice(), 44),
            (b"Rh".as_slice(), 45),
            (b"Pd".as_slice(), 46),
            (b"Ag".as_slice(), 47),
            (b"Cd".as_slice(), 48),
            (b"In".as_slice(), 49),
            (b"Sn".as_slice(), 50),
            (b"Sb".as_slice(), 51),
            (b"Te".as_slice(), 52),
            (b"I".as_slice(), 53),
            (b"Xe".as_slice(), 54),
            (b"Cs".as_slice(), 55),
            (b"Ba".as_slice(), 56),
            (b"La".as_slice(), 57),
            (b"Ce".as_slice(), 58),
            (b"Pr".as_slice(), 59),
            (b"Nd".as_slice(), 60),
            (b"Pm".as_slice(), 61),
            (b"Sm".as_slice(), 62),
            (b"Eu".as_slice(), 63),
            (b"Gd".as_slice(), 64),
            (b"Tb".as_slice(), 65),
            (b"Dy".as_slice(), 66),
            (b"Ho".as_slice(), 67),
            (b"Er".as_slice(), 68),
            (b"Tm".as_slice(), 69),
            (b"Yb".as_slice(), 70),
            (b"Lu".as_slice(), 71),
            (b"Hf".as_slice(), 72),
            (b"Ta".as_slice(), 73),
            (b"W".as_slice(), 74),
            (b"Re".as_slice(), 75),
            (b"Os".as_slice(), 76),
            (b"Ir".as_slice(), 77),
            (b"Pt".as_slice(), 78),
            (b"Au".as_slice(), 79),
            (b"Hg".as_slice(), 80),
            (b"Tl".as_slice(), 81),
            (b"Pb".as_slice(), 82),
            (b"Bi".as_slice(), 83),
            (b"Po".as_slice(), 84),
            (b"At".as_slice(), 85),
            (b"Rn".as_slice(), 86),
            (b"Fr".as_slice(), 87),
            (b"Ra".as_slice(), 88),
            (b"Ac".as_slice(), 89),
            (b"Th".as_slice(), 90),
            (b"Pa".as_slice(), 91),
            (b"U".as_slice(), 92),
            (b"Np".as_slice(), 93),
            (b"Pu".as_slice(), 94),
            (b"Am".as_slice(), 95),
            (b"Cm".as_slice(), 96),
            (b"Bk".as_slice(), 97),
            (b"Cf".as_slice(), 98),
            (b"Es".as_slice(), 99),
            (b"Fm".as_slice(), 100),
            (b"Md".as_slice(), 101),
            (b"No".as_slice(), 102),
            (b"Lr".as_slice(), 103),
            (b"Rf".as_slice(), 104),
            (b"Db".as_slice(), 105),
            (b"Sg".as_slice(), 106),
            (b"Bh".as_slice(), 107),
            (b"Hs".as_slice(), 108),
            (b"Mt".as_slice(), 109),
            (b"Ds".as_slice(), 110),
            (b"Rg".as_slice(), 111),
            (b"Cn".as_slice(), 112),
            (b"Nh".as_slice(), 113),
            (b"Uut".as_slice(), 113),
            (b"Fl".as_slice(), 114),
            (b"Mc".as_slice(), 115),
            (b"Uup".as_slice(), 115),
            (b"Lv".as_slice(), 116),
            (b"Ts".as_slice(), 117),
            (b"Og".as_slice(), 118),
        ] {
            let mut bytes = symbol.to_vec();
            bytes.push(0);
            let c_symbol = CStr::from_bytes_with_nul(&bytes).unwrap();
            assert_eq!(rdkit_atomic_number_from_c_symbol(c_symbol), Ok(number));
        }
    }

    #[test]
    fn first_nul_terminates_the_source_constructor_before_lookup() {
        for (bytes, number) in [
            (b"C\0N".as_slice(), 6),
            (b"Cl\0invalid".as_slice(), 17),
            (b"Uut\0\xff".as_slice(), 113),
        ] {
            let c_symbol = CStr::from_bytes_until_nul(bytes).unwrap();
            assert_eq!(rdkit_atomic_number_from_c_symbol(c_symbol), Ok(number));
            assert_eq!(
                rdkit_atomic_number_from_symbol(bytes)
                    .unwrap_err()
                    .symbol
                    .as_bytes(),
                bytes
            );
        }
    }

    #[test]
    fn empty_case_whitespace_and_invalid_bytes_preserve_checked_error_payload() {
        for (raw, expected) in [
            (b"\0C".as_slice(), b"".as_slice()),
            (b"c\0C", b"c"),
            (b"C \0", b"C "),
            (b"\xff\0C", b"\xff"),
            (b"Uuo\0Og", b"Uuo"),
        ] {
            let error = rdkit_atomic_number_from_c_symbol(CStr::from_bytes_until_nul(raw).unwrap())
                .unwrap_err();
            assert_eq!(error.symbol.as_bytes(), expected);
        }
        assert_eq!(
            rdkit_atomic_number_from_c_symbol(c"Bogus")
                .unwrap_err()
                .to_string(),
            "Element 'Bogus' not found"
        );
    }
}

#[cfg(test)]
mod complete_default_valence_symbol_source_tests {
    use super::*;

    #[test]
    fn all_pinned_name_rows_use_actual_first_signed_numeric_valence() {
        // Source byanum retains the first numerical row; byname retains both
        // Uut/Uup names. Expected signed defaults come from pinned source data.
        for (symbol, expected) in [
            (b"*".as_slice(), -1),
            (b"H".as_slice(), 1),
            (b"He".as_slice(), 0),
            (b"Li".as_slice(), 1),
            (b"Be".as_slice(), 2),
            (b"B".as_slice(), 3),
            (b"C".as_slice(), 4),
            (b"N".as_slice(), 3),
            (b"O".as_slice(), 2),
            (b"F".as_slice(), 1),
            (b"Ne".as_slice(), 0),
            (b"Na".as_slice(), 1),
            (b"Mg".as_slice(), 2),
            (b"Al".as_slice(), 3),
            (b"Si".as_slice(), 4),
            (b"P".as_slice(), 3),
            (b"S".as_slice(), 2),
            (b"Cl".as_slice(), 1),
            (b"Ar".as_slice(), 0),
            (b"K".as_slice(), 1),
            (b"Ca".as_slice(), 2),
            (b"Sc".as_slice(), -1),
            (b"Ti".as_slice(), -1),
            (b"V".as_slice(), -1),
            (b"Cr".as_slice(), -1),
            (b"Mn".as_slice(), -1),
            (b"Fe".as_slice(), -1),
            (b"Co".as_slice(), -1),
            (b"Ni".as_slice(), -1),
            (b"Cu".as_slice(), -1),
            (b"Zn".as_slice(), -1),
            (b"Ga".as_slice(), 3),
            (b"Ge".as_slice(), 4),
            (b"As".as_slice(), 3),
            (b"Se".as_slice(), 2),
            (b"Br".as_slice(), 1),
            (b"Kr".as_slice(), 0),
            (b"Rb".as_slice(), 1),
            (b"Sr".as_slice(), 2),
            (b"Y".as_slice(), -1),
            (b"Zr".as_slice(), -1),
            (b"Nb".as_slice(), -1),
            (b"Mo".as_slice(), -1),
            (b"Tc".as_slice(), -1),
            (b"Ru".as_slice(), -1),
            (b"Rh".as_slice(), -1),
            (b"Pd".as_slice(), -1),
            (b"Ag".as_slice(), -1),
            (b"Cd".as_slice(), -1),
            (b"In".as_slice(), 3),
            (b"Sn".as_slice(), 2),
            (b"Sb".as_slice(), 3),
            (b"Te".as_slice(), 2),
            (b"I".as_slice(), 1),
            (b"Xe".as_slice(), 0),
            (b"Cs".as_slice(), 1),
            (b"Ba".as_slice(), 2),
            (b"La".as_slice(), -1),
            (b"Ce".as_slice(), -1),
            (b"Pr".as_slice(), -1),
            (b"Nd".as_slice(), -1),
            (b"Pm".as_slice(), -1),
            (b"Sm".as_slice(), -1),
            (b"Eu".as_slice(), -1),
            (b"Gd".as_slice(), -1),
            (b"Tb".as_slice(), -1),
            (b"Dy".as_slice(), -1),
            (b"Ho".as_slice(), -1),
            (b"Er".as_slice(), -1),
            (b"Tm".as_slice(), -1),
            (b"Yb".as_slice(), -1),
            (b"Lu".as_slice(), -1),
            (b"Hf".as_slice(), -1),
            (b"Ta".as_slice(), -1),
            (b"W".as_slice(), -1),
            (b"Re".as_slice(), -1),
            (b"Os".as_slice(), -1),
            (b"Ir".as_slice(), -1),
            (b"Pt".as_slice(), -1),
            (b"Au".as_slice(), -1),
            (b"Hg".as_slice(), -1),
            (b"Tl".as_slice(), -1),
            (b"Pb".as_slice(), 2),
            (b"Bi".as_slice(), 3),
            (b"Po".as_slice(), 2),
            (b"At".as_slice(), 1),
            (b"Rn".as_slice(), 0),
            (b"Fr".as_slice(), 1),
            (b"Ra".as_slice(), 2),
            (b"Ac".as_slice(), -1),
            (b"Th".as_slice(), -1),
            (b"Pa".as_slice(), -1),
            (b"U".as_slice(), -1),
            (b"Np".as_slice(), -1),
            (b"Pu".as_slice(), -1),
            (b"Am".as_slice(), -1),
            (b"Cm".as_slice(), -1),
            (b"Bk".as_slice(), -1),
            (b"Cf".as_slice(), -1),
            (b"Es".as_slice(), -1),
            (b"Fm".as_slice(), -1),
            (b"Md".as_slice(), -1),
            (b"No".as_slice(), -1),
            (b"Lr".as_slice(), -1),
            (b"Rf".as_slice(), -1),
            (b"Db".as_slice(), -1),
            (b"Sg".as_slice(), -1),
            (b"Bh".as_slice(), -1),
            (b"Hs".as_slice(), -1),
            (b"Mt".as_slice(), -1),
            (b"Ds".as_slice(), -1),
            (b"Rg".as_slice(), -1),
            (b"Cn".as_slice(), -1),
            (b"Nh".as_slice(), -1),
            (b"Uut".as_slice(), -1),
            (b"Fl".as_slice(), -1),
            (b"Mc".as_slice(), -1),
            (b"Uup".as_slice(), -1),
            (b"Lv".as_slice(), -1),
            (b"Ts".as_slice(), -1),
            (b"Og".as_slice(), -1),
        ] {
            assert_eq!(
                rdkit_default_valence_from_symbol(symbol),
                Ok(expected),
                "{symbol:?}"
            );
        }
    }

    #[test]
    fn unknown_byte_symbol_fails_precondition_before_numeric_default_lookup() {
        for raw in [
            b"".as_slice(),
            b"c",
            b"C ",
            b" C",
            b"C\0",
            b"C\0N",
            b"D",
            b"T",
            b"Uuo",
            b"\xff",
            b"\xc3\xa9",
        ] {
            let error = rdkit_default_valence_from_symbol(raw).unwrap_err();
            assert_eq!(
                error,
                SymbolValenceLookupError::SymbolNotFound {
                    symbol: cosmolkit_model::PropertyText::from_bytes(raw),
                }
            );
        }
        assert_eq!(
            rdkit_default_valence_from_symbol("Unknown")
                .unwrap_err()
                .to_string(),
            "Element 'Unknown' not found"
        );
    }
}

#[cfg(test)]
mod complete_default_valence_c_source_tests {
    use super::*;
    use std::ffi::CStr;

    #[test]
    fn all_source_names_preserve_actual_signed_default_after_c_string_construction() {
        for (symbol, expected) in [
            (b"*".as_slice(), -1),
            (b"H".as_slice(), 1),
            (b"He".as_slice(), 0),
            (b"Li".as_slice(), 1),
            (b"Be".as_slice(), 2),
            (b"B".as_slice(), 3),
            (b"C".as_slice(), 4),
            (b"N".as_slice(), 3),
            (b"O".as_slice(), 2),
            (b"F".as_slice(), 1),
            (b"Ne".as_slice(), 0),
            (b"Na".as_slice(), 1),
            (b"Mg".as_slice(), 2),
            (b"Al".as_slice(), 3),
            (b"Si".as_slice(), 4),
            (b"P".as_slice(), 3),
            (b"S".as_slice(), 2),
            (b"Cl".as_slice(), 1),
            (b"Ar".as_slice(), 0),
            (b"K".as_slice(), 1),
            (b"Ca".as_slice(), 2),
            (b"Sc".as_slice(), -1),
            (b"Ti".as_slice(), -1),
            (b"V".as_slice(), -1),
            (b"Cr".as_slice(), -1),
            (b"Mn".as_slice(), -1),
            (b"Fe".as_slice(), -1),
            (b"Co".as_slice(), -1),
            (b"Ni".as_slice(), -1),
            (b"Cu".as_slice(), -1),
            (b"Zn".as_slice(), -1),
            (b"Ga".as_slice(), 3),
            (b"Ge".as_slice(), 4),
            (b"As".as_slice(), 3),
            (b"Se".as_slice(), 2),
            (b"Br".as_slice(), 1),
            (b"Kr".as_slice(), 0),
            (b"Rb".as_slice(), 1),
            (b"Sr".as_slice(), 2),
            (b"Y".as_slice(), -1),
            (b"Zr".as_slice(), -1),
            (b"Nb".as_slice(), -1),
            (b"Mo".as_slice(), -1),
            (b"Tc".as_slice(), -1),
            (b"Ru".as_slice(), -1),
            (b"Rh".as_slice(), -1),
            (b"Pd".as_slice(), -1),
            (b"Ag".as_slice(), -1),
            (b"Cd".as_slice(), -1),
            (b"In".as_slice(), 3),
            (b"Sn".as_slice(), 2),
            (b"Sb".as_slice(), 3),
            (b"Te".as_slice(), 2),
            (b"I".as_slice(), 1),
            (b"Xe".as_slice(), 0),
            (b"Cs".as_slice(), 1),
            (b"Ba".as_slice(), 2),
            (b"La".as_slice(), -1),
            (b"Ce".as_slice(), -1),
            (b"Pr".as_slice(), -1),
            (b"Nd".as_slice(), -1),
            (b"Pm".as_slice(), -1),
            (b"Sm".as_slice(), -1),
            (b"Eu".as_slice(), -1),
            (b"Gd".as_slice(), -1),
            (b"Tb".as_slice(), -1),
            (b"Dy".as_slice(), -1),
            (b"Ho".as_slice(), -1),
            (b"Er".as_slice(), -1),
            (b"Tm".as_slice(), -1),
            (b"Yb".as_slice(), -1),
            (b"Lu".as_slice(), -1),
            (b"Hf".as_slice(), -1),
            (b"Ta".as_slice(), -1),
            (b"W".as_slice(), -1),
            (b"Re".as_slice(), -1),
            (b"Os".as_slice(), -1),
            (b"Ir".as_slice(), -1),
            (b"Pt".as_slice(), -1),
            (b"Au".as_slice(), -1),
            (b"Hg".as_slice(), -1),
            (b"Tl".as_slice(), -1),
            (b"Pb".as_slice(), 2),
            (b"Bi".as_slice(), 3),
            (b"Po".as_slice(), 2),
            (b"At".as_slice(), 1),
            (b"Rn".as_slice(), 0),
            (b"Fr".as_slice(), 1),
            (b"Ra".as_slice(), 2),
            (b"Ac".as_slice(), -1),
            (b"Th".as_slice(), -1),
            (b"Pa".as_slice(), -1),
            (b"U".as_slice(), -1),
            (b"Np".as_slice(), -1),
            (b"Pu".as_slice(), -1),
            (b"Am".as_slice(), -1),
            (b"Cm".as_slice(), -1),
            (b"Bk".as_slice(), -1),
            (b"Cf".as_slice(), -1),
            (b"Es".as_slice(), -1),
            (b"Fm".as_slice(), -1),
            (b"Md".as_slice(), -1),
            (b"No".as_slice(), -1),
            (b"Lr".as_slice(), -1),
            (b"Rf".as_slice(), -1),
            (b"Db".as_slice(), -1),
            (b"Sg".as_slice(), -1),
            (b"Bh".as_slice(), -1),
            (b"Hs".as_slice(), -1),
            (b"Mt".as_slice(), -1),
            (b"Ds".as_slice(), -1),
            (b"Rg".as_slice(), -1),
            (b"Cn".as_slice(), -1),
            (b"Nh".as_slice(), -1),
            (b"Uut".as_slice(), -1),
            (b"Fl".as_slice(), -1),
            (b"Mc".as_slice(), -1),
            (b"Uup".as_slice(), -1),
            (b"Lv".as_slice(), -1),
            (b"Ts".as_slice(), -1),
            (b"Og".as_slice(), -1),
        ] {
            let mut raw = symbol.to_vec();
            raw.push(0);
            assert_eq!(
                rdkit_default_valence_from_c_symbol(CStr::from_bytes_with_nul(&raw).unwrap()),
                Ok(expected)
            );
        }
    }

    #[test]
    fn terminator_precedes_lookup_and_error_payload_uses_only_constructed_symbol() {
        for (raw, expected) in [
            (b"O\0N".as_slice(), 2),
            (b"Uut\0N", -1),
            (b"*\0C", -1),
            (b"Xe\0S", 0),
        ] {
            assert_eq!(
                rdkit_default_valence_from_c_symbol(CStr::from_bytes_until_nul(raw).unwrap()),
                Ok(expected)
            );
            assert!(rdkit_default_valence_from_symbol(raw).is_err());
        }
        for (raw, expected) in [
            (b"\0C".as_slice(), b"".as_slice()),
            (b"c\0C", b"c"),
            (b"C \0N", b"C "),
            (b"\xff\0C", b"\xff"),
            (b"D\0C", b"D"),
        ] {
            assert_eq!(
                rdkit_default_valence_from_c_symbol(CStr::from_bytes_until_nul(raw).unwrap()),
                Err(SymbolValenceLookupError::SymbolNotFound {
                    symbol: cosmolkit_model::PropertyText::from_bytes(expected),
                })
            );
        }
    }
}

#[cfg(test)]
mod complete_valence_list_symbol_source_tests {
    use super::*;

    #[test]
    fn every_source_name_returns_entire_signed_list_in_pinned_order() {
        // Expected arrays mechanically prepared from pinned first byanum rows.
        for (symbol, expected) in [
            (b"*".as_slice(), &[-1][..]),
            (b"H".as_slice(), &[1][..]),
            (b"He".as_slice(), &[0][..]),
            (b"Li".as_slice(), &[1, -1][..]),
            (b"Be".as_slice(), &[2][..]),
            (b"B".as_slice(), &[3][..]),
            (b"C".as_slice(), &[4][..]),
            (b"N".as_slice(), &[3][..]),
            (b"O".as_slice(), &[2][..]),
            (b"F".as_slice(), &[1][..]),
            (b"Ne".as_slice(), &[0][..]),
            (b"Na".as_slice(), &[1, -1][..]),
            (b"Mg".as_slice(), &[2, -1][..]),
            (b"Al".as_slice(), &[3][..]),
            (b"Si".as_slice(), &[4][..]),
            (b"P".as_slice(), &[3, 5][..]),
            (b"S".as_slice(), &[2, 4, 6][..]),
            (b"Cl".as_slice(), &[1][..]),
            (b"Ar".as_slice(), &[0][..]),
            (b"K".as_slice(), &[1, -1][..]),
            (b"Ca".as_slice(), &[2, -1][..]),
            (b"Sc".as_slice(), &[-1][..]),
            (b"Ti".as_slice(), &[-1][..]),
            (b"V".as_slice(), &[-1][..]),
            (b"Cr".as_slice(), &[-1][..]),
            (b"Mn".as_slice(), &[-1][..]),
            (b"Fe".as_slice(), &[-1][..]),
            (b"Co".as_slice(), &[-1][..]),
            (b"Ni".as_slice(), &[-1][..]),
            (b"Cu".as_slice(), &[-1][..]),
            (b"Zn".as_slice(), &[-1][..]),
            (b"Ga".as_slice(), &[3][..]),
            (b"Ge".as_slice(), &[4][..]),
            (b"As".as_slice(), &[3, 5][..]),
            (b"Se".as_slice(), &[2, 4, 6][..]),
            (b"Br".as_slice(), &[1][..]),
            (b"Kr".as_slice(), &[0][..]),
            (b"Rb".as_slice(), &[1, -1][..]),
            (b"Sr".as_slice(), &[2, -1][..]),
            (b"Y".as_slice(), &[-1][..]),
            (b"Zr".as_slice(), &[-1][..]),
            (b"Nb".as_slice(), &[-1][..]),
            (b"Mo".as_slice(), &[-1][..]),
            (b"Tc".as_slice(), &[-1][..]),
            (b"Ru".as_slice(), &[-1][..]),
            (b"Rh".as_slice(), &[-1][..]),
            (b"Pd".as_slice(), &[-1][..]),
            (b"Ag".as_slice(), &[-1][..]),
            (b"Cd".as_slice(), &[-1][..]),
            (b"In".as_slice(), &[3][..]),
            (b"Sn".as_slice(), &[2, 4][..]),
            (b"Sb".as_slice(), &[3, 5][..]),
            (b"Te".as_slice(), &[2, 4, 6][..]),
            (b"I".as_slice(), &[1, 3, 5][..]),
            (b"Xe".as_slice(), &[0, 2, 4, 6][..]),
            (b"Cs".as_slice(), &[1][..]),
            (b"Ba".as_slice(), &[2, -1][..]),
            (b"La".as_slice(), &[-1][..]),
            (b"Ce".as_slice(), &[-1][..]),
            (b"Pr".as_slice(), &[-1][..]),
            (b"Nd".as_slice(), &[-1][..]),
            (b"Pm".as_slice(), &[-1][..]),
            (b"Sm".as_slice(), &[-1][..]),
            (b"Eu".as_slice(), &[-1][..]),
            (b"Gd".as_slice(), &[-1][..]),
            (b"Tb".as_slice(), &[-1][..]),
            (b"Dy".as_slice(), &[-1][..]),
            (b"Ho".as_slice(), &[-1][..]),
            (b"Er".as_slice(), &[-1][..]),
            (b"Tm".as_slice(), &[-1][..]),
            (b"Yb".as_slice(), &[-1][..]),
            (b"Lu".as_slice(), &[-1][..]),
            (b"Hf".as_slice(), &[-1][..]),
            (b"Ta".as_slice(), &[-1][..]),
            (b"W".as_slice(), &[-1][..]),
            (b"Re".as_slice(), &[-1][..]),
            (b"Os".as_slice(), &[-1][..]),
            (b"Ir".as_slice(), &[-1][..]),
            (b"Pt".as_slice(), &[-1][..]),
            (b"Au".as_slice(), &[-1][..]),
            (b"Hg".as_slice(), &[-1][..]),
            (b"Tl".as_slice(), &[-1][..]),
            (b"Pb".as_slice(), &[2, 4][..]),
            (b"Bi".as_slice(), &[3, 5][..]),
            (b"Po".as_slice(), &[2, 4, 6][..]),
            (b"At".as_slice(), &[1, 3, 5][..]),
            (b"Rn".as_slice(), &[0][..]),
            (b"Fr".as_slice(), &[1][..]),
            (b"Ra".as_slice(), &[2, -1][..]),
            (b"Ac".as_slice(), &[-1][..]),
            (b"Th".as_slice(), &[-1][..]),
            (b"Pa".as_slice(), &[-1][..]),
            (b"U".as_slice(), &[-1][..]),
            (b"Np".as_slice(), &[-1][..]),
            (b"Pu".as_slice(), &[-1][..]),
            (b"Am".as_slice(), &[-1][..]),
            (b"Cm".as_slice(), &[-1][..]),
            (b"Bk".as_slice(), &[-1][..]),
            (b"Cf".as_slice(), &[-1][..]),
            (b"Es".as_slice(), &[-1][..]),
            (b"Fm".as_slice(), &[-1][..]),
            (b"Md".as_slice(), &[-1][..]),
            (b"No".as_slice(), &[-1][..]),
            (b"Lr".as_slice(), &[-1][..]),
            (b"Rf".as_slice(), &[-1][..]),
            (b"Db".as_slice(), &[-1][..]),
            (b"Sg".as_slice(), &[-1][..]),
            (b"Bh".as_slice(), &[-1][..]),
            (b"Hs".as_slice(), &[-1][..]),
            (b"Mt".as_slice(), &[-1][..]),
            (b"Ds".as_slice(), &[-1][..]),
            (b"Rg".as_slice(), &[-1][..]),
            (b"Cn".as_slice(), &[-1][..]),
            (b"Nh".as_slice(), &[-1][..]),
            (b"Uut".as_slice(), &[-1][..]),
            (b"Fl".as_slice(), &[-1][..]),
            (b"Mc".as_slice(), &[-1][..]),
            (b"Uup".as_slice(), &[-1][..]),
            (b"Lv".as_slice(), &[-1][..]),
            (b"Ts".as_slice(), &[-1][..]),
            (b"Og".as_slice(), &[-1][..]),
        ] {
            assert_eq!(
                rdkit_valence_list_from_symbol(symbol),
                Ok(expected),
                "{symbol:?}"
            );
        }
    }

    #[test]
    fn missing_raw_name_is_precondition_error_without_empty_or_unlimited_list() {
        for raw in [
            b"".as_slice(),
            b"c",
            b"C ",
            b" C",
            b"C\0",
            b"C\0N",
            b"D",
            b"T",
            b"Uuo",
            b"\xff",
        ] {
            assert_eq!(
                rdkit_valence_list_from_symbol(raw),
                Err(SymbolValenceLookupError::SymbolNotFound {
                    symbol: cosmolkit_model::PropertyText::from_bytes(raw),
                })
            );
        }
    }

    #[test]
    fn aliases_and_temporary_input_share_the_stable_actual_numeric_slice() {
        for (canonical, alias, number) in [("Nh", "Uut", 113), ("Mc", "Uup", 115)] {
            let first = rdkit_valence_list_from_symbol(canonical).unwrap();
            let second = rdkit_valence_list_from_symbol(alias).unwrap();
            assert!(std::ptr::eq(first, second));
            assert!(std::ptr::eq(first, required_valence_list(number).unwrap()));
        }
        let list = {
            let owned_symbol = String::from("Xe");
            rdkit_valence_list_from_symbol(&owned_symbol).unwrap()
        };
        assert_eq!(list, &[0, 2, 4, 6]);
    }
}

#[cfg(test)]
mod complete_valence_list_c_source_tests {
    use super::*;
    use std::ffi::{CStr, CString};

    #[test]
    fn every_actual_c_name_returns_complete_pinned_signed_list() {
        for (symbol, expected) in [
            (b"*".as_slice(), &[-1][..]),
            (b"H".as_slice(), &[1][..]),
            (b"He".as_slice(), &[0][..]),
            (b"Li".as_slice(), &[1, -1][..]),
            (b"Be".as_slice(), &[2][..]),
            (b"B".as_slice(), &[3][..]),
            (b"C".as_slice(), &[4][..]),
            (b"N".as_slice(), &[3][..]),
            (b"O".as_slice(), &[2][..]),
            (b"F".as_slice(), &[1][..]),
            (b"Ne".as_slice(), &[0][..]),
            (b"Na".as_slice(), &[1, -1][..]),
            (b"Mg".as_slice(), &[2, -1][..]),
            (b"Al".as_slice(), &[3][..]),
            (b"Si".as_slice(), &[4][..]),
            (b"P".as_slice(), &[3, 5][..]),
            (b"S".as_slice(), &[2, 4, 6][..]),
            (b"Cl".as_slice(), &[1][..]),
            (b"Ar".as_slice(), &[0][..]),
            (b"K".as_slice(), &[1, -1][..]),
            (b"Ca".as_slice(), &[2, -1][..]),
            (b"Sc".as_slice(), &[-1][..]),
            (b"Ti".as_slice(), &[-1][..]),
            (b"V".as_slice(), &[-1][..]),
            (b"Cr".as_slice(), &[-1][..]),
            (b"Mn".as_slice(), &[-1][..]),
            (b"Fe".as_slice(), &[-1][..]),
            (b"Co".as_slice(), &[-1][..]),
            (b"Ni".as_slice(), &[-1][..]),
            (b"Cu".as_slice(), &[-1][..]),
            (b"Zn".as_slice(), &[-1][..]),
            (b"Ga".as_slice(), &[3][..]),
            (b"Ge".as_slice(), &[4][..]),
            (b"As".as_slice(), &[3, 5][..]),
            (b"Se".as_slice(), &[2, 4, 6][..]),
            (b"Br".as_slice(), &[1][..]),
            (b"Kr".as_slice(), &[0][..]),
            (b"Rb".as_slice(), &[1, -1][..]),
            (b"Sr".as_slice(), &[2, -1][..]),
            (b"Y".as_slice(), &[-1][..]),
            (b"Zr".as_slice(), &[-1][..]),
            (b"Nb".as_slice(), &[-1][..]),
            (b"Mo".as_slice(), &[-1][..]),
            (b"Tc".as_slice(), &[-1][..]),
            (b"Ru".as_slice(), &[-1][..]),
            (b"Rh".as_slice(), &[-1][..]),
            (b"Pd".as_slice(), &[-1][..]),
            (b"Ag".as_slice(), &[-1][..]),
            (b"Cd".as_slice(), &[-1][..]),
            (b"In".as_slice(), &[3][..]),
            (b"Sn".as_slice(), &[2, 4][..]),
            (b"Sb".as_slice(), &[3, 5][..]),
            (b"Te".as_slice(), &[2, 4, 6][..]),
            (b"I".as_slice(), &[1, 3, 5][..]),
            (b"Xe".as_slice(), &[0, 2, 4, 6][..]),
            (b"Cs".as_slice(), &[1][..]),
            (b"Ba".as_slice(), &[2, -1][..]),
            (b"La".as_slice(), &[-1][..]),
            (b"Ce".as_slice(), &[-1][..]),
            (b"Pr".as_slice(), &[-1][..]),
            (b"Nd".as_slice(), &[-1][..]),
            (b"Pm".as_slice(), &[-1][..]),
            (b"Sm".as_slice(), &[-1][..]),
            (b"Eu".as_slice(), &[-1][..]),
            (b"Gd".as_slice(), &[-1][..]),
            (b"Tb".as_slice(), &[-1][..]),
            (b"Dy".as_slice(), &[-1][..]),
            (b"Ho".as_slice(), &[-1][..]),
            (b"Er".as_slice(), &[-1][..]),
            (b"Tm".as_slice(), &[-1][..]),
            (b"Yb".as_slice(), &[-1][..]),
            (b"Lu".as_slice(), &[-1][..]),
            (b"Hf".as_slice(), &[-1][..]),
            (b"Ta".as_slice(), &[-1][..]),
            (b"W".as_slice(), &[-1][..]),
            (b"Re".as_slice(), &[-1][..]),
            (b"Os".as_slice(), &[-1][..]),
            (b"Ir".as_slice(), &[-1][..]),
            (b"Pt".as_slice(), &[-1][..]),
            (b"Au".as_slice(), &[-1][..]),
            (b"Hg".as_slice(), &[-1][..]),
            (b"Tl".as_slice(), &[-1][..]),
            (b"Pb".as_slice(), &[2, 4][..]),
            (b"Bi".as_slice(), &[3, 5][..]),
            (b"Po".as_slice(), &[2, 4, 6][..]),
            (b"At".as_slice(), &[1, 3, 5][..]),
            (b"Rn".as_slice(), &[0][..]),
            (b"Fr".as_slice(), &[1][..]),
            (b"Ra".as_slice(), &[2, -1][..]),
            (b"Ac".as_slice(), &[-1][..]),
            (b"Th".as_slice(), &[-1][..]),
            (b"Pa".as_slice(), &[-1][..]),
            (b"U".as_slice(), &[-1][..]),
            (b"Np".as_slice(), &[-1][..]),
            (b"Pu".as_slice(), &[-1][..]),
            (b"Am".as_slice(), &[-1][..]),
            (b"Cm".as_slice(), &[-1][..]),
            (b"Bk".as_slice(), &[-1][..]),
            (b"Cf".as_slice(), &[-1][..]),
            (b"Es".as_slice(), &[-1][..]),
            (b"Fm".as_slice(), &[-1][..]),
            (b"Md".as_slice(), &[-1][..]),
            (b"No".as_slice(), &[-1][..]),
            (b"Lr".as_slice(), &[-1][..]),
            (b"Rf".as_slice(), &[-1][..]),
            (b"Db".as_slice(), &[-1][..]),
            (b"Sg".as_slice(), &[-1][..]),
            (b"Bh".as_slice(), &[-1][..]),
            (b"Hs".as_slice(), &[-1][..]),
            (b"Mt".as_slice(), &[-1][..]),
            (b"Ds".as_slice(), &[-1][..]),
            (b"Rg".as_slice(), &[-1][..]),
            (b"Cn".as_slice(), &[-1][..]),
            (b"Nh".as_slice(), &[-1][..]),
            (b"Uut".as_slice(), &[-1][..]),
            (b"Fl".as_slice(), &[-1][..]),
            (b"Mc".as_slice(), &[-1][..]),
            (b"Uup".as_slice(), &[-1][..]),
            (b"Lv".as_slice(), &[-1][..]),
            (b"Ts".as_slice(), &[-1][..]),
            (b"Og".as_slice(), &[-1][..]),
        ] {
            let mut raw = symbol.to_vec();
            raw.push(0);
            assert_eq!(
                rdkit_valence_list_from_c_symbol(CStr::from_bytes_with_nul(&raw).unwrap()),
                Ok(expected),
                "{symbol:?}"
            );
        }
    }

    #[test]
    fn first_nul_controls_name_and_precondition_error_bytes() {
        for (raw, expected) in [
            (b"S\0C".as_slice(), &[2, 4, 6][..]),
            (b"Uut\0Xe", &[-1][..]),
            (b"Xe\0\xff", &[0, 2, 4, 6][..]),
            (b"*\0C", &[-1][..]),
        ] {
            assert_eq!(
                rdkit_valence_list_from_c_symbol(CStr::from_bytes_until_nul(raw).unwrap()),
                Ok(expected)
            );
            assert!(rdkit_valence_list_from_symbol(raw).is_err());
        }
        for (raw, expected) in [
            (b"\0C".as_slice(), b"".as_slice()),
            (b"c\0C", b"c"),
            (b"C \0", b"C "),
            (b"\xff\0C", b"\xff"),
            (b"D\0C", b"D"),
        ] {
            assert_eq!(
                rdkit_valence_list_from_c_symbol(CStr::from_bytes_until_nul(raw).unwrap()),
                Err(SymbolValenceLookupError::SymbolNotFound {
                    symbol: cosmolkit_model::PropertyText::from_bytes(expected),
                })
            );
        }
    }

    #[test]
    fn actual_numeric_storage_outlives_c_string_and_is_shared_by_alias() {
        let slice = {
            let owned = CString::new("Uup").unwrap();
            rdkit_valence_list_from_c_symbol(&owned).unwrap()
        };
        assert!(std::ptr::eq(
            slice,
            rdkit_valence_list_from_c_symbol(c"Mc").unwrap()
        ));
        assert!(std::ptr::eq(slice, required_valence_list(115).unwrap()));
        assert_eq!(slice, &[-1]);
    }
}

/// Apply the sole inherited Atom cache body to a query carrier, retaining
/// each successful signed8 store even if a later source calculation errors.
#[doc(hidden)]
pub fn update_query_atom_property_cache_source(
    atoms: &mut [cosmolkit_model::QueryAtom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    id: AtomId,
    strict: bool,
    source_complex_bonds: Option<&[bool]>,
) -> Result<(), ValenceError> {
    // RDKit❗✔️: void Atom::updatePropertyCache(bool strict) {
    // RDKit❗✔️:   calcExplicitValence(strict);
    // RDKit❗✔️:   calcImplicitValence(strict);
    // RDKit❗✔️: }
    // QueryAtom does not override these nonvirtual Atom functions. The closed
    // getter trait reaches the actual u8 identity, including invalid table IDs.
    // Shared helper preserves source store/error order without cloning atoms,
    // their property dictionaries, predicates, or temporary default carriers.
    let mut facts = atoms
        .get(id.index())
        .ok_or(ValenceError::AtomOutOfRange {
            atom: id,
            atom_count: atoms.len(),
        })?
        .source_valence_facts();
    let result = calculate_cached_valence_state_with_facts(
        atoms,
        bonds,
        adjacency,
        id,
        &mut facts,
        strict,
        source_complex_bonds,
    );
    atoms[id.index()].set_source_valence_facts(facts);
    result.map(|_| ())
}
