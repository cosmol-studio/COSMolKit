// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use cosmolkit_model::{
    Atom, AtomId, Bond, BondId, BondValueError, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Hybridization};

use crate::{
    DoubleBondControl, DoubleBondStereoDescriptor, DoubleBondStereoError,
    DoubleBondStereoSpecified, RingInfo, StereoOrderError, ValenceAssignment, ValenceError,
    bond_affects_atom_chirality, count_swaps_to_interconvert, double_bond_stereo_info,
    is_atom_bridgehead_from_topology, total_hydrogen_count,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PotentialStereoParams {
    pub clean: bool,
    pub flag_possible: bool,
    pub allow_nontetrahedral: bool,
}

impl Default for PotentialStereoParams {
    fn default() -> Self {
        Self {
            clean: false,
            flag_possible: true,
            allow_nontetrahedral: true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PotentialStereoType {
    AtomTetrahedral,
    AtomSquarePlanar,
    AtomTrigonalBipyramidal,
    AtomOctahedral,
    BondDouble,
    BondCumuleneEven,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PotentialStereoSpecified {
    Unspecified,
    Specified,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PotentialStereoDescriptor {
    None,
    TetrahedralClockwise,
    TetrahedralCounterclockwise,
    BondCis,
    BondTrans,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PotentialStereoCenter {
    Atom(AtomId),
    Bond(BondId),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PotentialStereoInfo {
    pub stereo_type: PotentialStereoType,
    pub specified: PotentialStereoSpecified,
    pub centered_on: PotentialStereoCenter,
    pub descriptor: PotentialStereoDescriptor,
    pub permutation: u32,
    pub controlling_atoms: Vec<Option<AtomId>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RingStereoRelation {
    pub atom: AtomId,
    pub other: AtomId,
    pub same_orientation: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PotentialStereoAssignment {
    pub stereo: Vec<PotentialStereoInfo>,
    pub atom_ranks: Vec<u32>,
    pub ring_relations: Vec<RingStereoRelation>,
    pub cleaned_topology: Option<TopologyBlock>,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PotentialStereoError {
    #[error("invalid topology: {0}")]
    InvalidTopology(#[from] TopologyValidationError),
    #[error("valence field {field} has {actual} rows, expected {atom_count}")]
    InvalidValence {
        field: &'static str,
        actual: usize,
        atom_count: usize,
    },
    #[error("valence field {field} has invalid value {value} for atom {atom}")]
    InvalidValenceValue {
        field: &'static str,
        atom: AtomId,
        value: i32,
    },
    #[error("invalid ring information: {reason} at row {row} (value {value}, limit {limit})")]
    InvalidRingInfo {
        reason: &'static str,
        row: usize,
        value: usize,
        limit: usize,
    },
    #[error("atom {atom} has invalid nonzero degree {degree}")]
    InvalidAtomDegree { atom: AtomId, degree: usize },
    #[error("bond {bond} {endpoint} endpoint has invalid degree {degree}")]
    InvalidBondDegree {
        bond: BondId,
        endpoint: &'static str,
        degree: usize,
    },
    #[error("bond {bond} has invalid stereo references: {reason}")]
    InvalidStereoReferences { bond: BondId, reason: &'static str },
    #[error("bond {bond} has unsupported order {order:?}")]
    UnsupportedBondOrder { bond: BondId, order: BondOrder },
    #[error("atom {atom} has invalid chiral permutation {permutation}")]
    InvalidChiralPermutation { atom: AtomId, permutation: u32 },
    #[error("atropisomer support is unavailable for bond {bond}")]
    AtropisomerDependencyUnavailable { bond: BondId },
    #[error("potential-stereo refinement did not converge after {iterations} iterations")]
    RefinementDidNotConverge { iterations: usize },
    #[error("valence calculation failed: {0}")]
    Valence(#[from] ValenceError),
    #[error("tetrahedral order calculation failed: {0}")]
    StereoOrder(#[from] StereoOrderError),
    #[error("double-bond stereo calculation failed: {0}")]
    DoubleStereo(#[from] DoubleBondStereoError),
    #[error("invalid detached bond value: {0}")]
    BondValue(#[from] BondValueError),
}

fn validate_inputs(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
) -> Result<(), PotentialStereoError> {
    topology.validate()?;
    let atom_count = topology.atoms.len();
    for (field, rows) in [
        ("explicit_valence", &valence.explicit_valence),
        ("implicit_hydrogens", &valence.implicit_hydrogens),
    ] {
        if rows.len() != atom_count {
            return Err(PotentialStereoError::InvalidValence {
                field,
                actual: rows.len(),
                atom_count,
            });
        }
        if let Some((index, value)) = rows
            .iter()
            .copied()
            .enumerate()
            .find(|(_, value)| *value < 0)
        {
            return Err(PotentialStereoError::InvalidValenceValue {
                field,
                atom: AtomId::new(index),
                value,
            });
        }
    }
    if !rings.is_symm_sssr() {
        return Err(PotentialStereoError::InvalidRingInfo {
            reason: "ring information is not symmetric SSSR",
            row: 0,
            value: rings.find_type() as usize,
            limit: 0,
        });
    }
    if rings.atom_row_count() != atom_count {
        return Err(PotentialStereoError::InvalidRingInfo {
            reason: "atom membership row count mismatch",
            row: 0,
            value: rings.atom_row_count(),
            limit: atom_count,
        });
    }
    if rings.bond_row_count() != topology.bonds.len() {
        return Err(PotentialStereoError::InvalidRingInfo {
            reason: "bond membership row count mismatch",
            row: 0,
            value: rings.bond_row_count(),
            limit: topology.bonds.len(),
        });
    }
    if rings.atom_rings().len() != rings.bond_rings().len() {
        return Err(PotentialStereoError::InvalidRingInfo {
            reason: "atom/bond ring table length mismatch",
            row: 0,
            value: rings.atom_rings().len(),
            limit: rings.bond_rings().len(),
        });
    }
    for (row, (atoms, bonds)) in rings
        .atom_rings()
        .iter()
        .zip(rings.bond_rings())
        .enumerate()
    {
        if atoms.len() != bonds.len() {
            return Err(PotentialStereoError::InvalidRingInfo {
                reason: "ring atom/bond size mismatch",
                row,
                value: atoms.len(),
                limit: bonds.len(),
            });
        }
        if let Some(atom) = atoms.iter().find(|atom| atom.index() >= atom_count) {
            return Err(PotentialStereoError::InvalidRingInfo {
                reason: "ring atom out of range",
                row,
                value: atom.index(),
                limit: atom_count,
            });
        }
        if let Some(bond) = bonds
            .iter()
            .find(|bond| bond.index() >= topology.bonds.len())
        {
            return Err(PotentialStereoError::InvalidRingInfo {
                reason: "ring bond out of range",
                row,
                value: bond.index(),
                limit: topology.bonds.len(),
            });
        }
    }
    for bond in &topology.bonds {
        if matches!(bond.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw) {
            return Err(PotentialStereoError::AtropisomerDependencyUnavailable { bond: bond.id() });
        }
    }
    Ok(())
}

fn graph_degree(topology: &TopologyBlock, atom: AtomId) -> usize {
    topology.adjacency.neighbors_of(atom.index()).len()
}

fn nonzero_degree(topology: &TopologyBlock, atom: AtomId) -> Result<usize, PotentialStereoError> {
    let mut degree = 0;
    for neighbor in topology.adjacency.neighbors_of(atom.index()) {
        let bond = &topology.bonds[neighbor.bond.index()];
        degree += usize::from(bond_affects_atom_chirality(bond, atom)?);
    }
    Ok(degree)
}

fn total_hydrogens(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom: AtomId,
    include_neighbors: bool,
) -> Result<usize, PotentialStereoError> {
    Ok(total_hydrogen_count(topology, valence, atom, include_neighbors)? as usize)
}

fn total_degree(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom: AtomId,
) -> Result<usize, PotentialStereoError> {
    let value = &topology.atoms[atom.index()];
    Ok(graph_degree(topology, atom)
        + usize::from(value.explicit_hydrogens())
        + usize::try_from(valence.implicit_hydrogens[atom.index()]).unwrap_or(0))
}

fn has_protium_neighbor(topology: &TopologyBlock, atom: AtomId) -> bool {
    topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .any(|neighbor| {
            let value = &topology.atoms[neighbor.atom_index];
            value.atomic_number() == 1 && value.isotope().is_none()
        })
}

fn has_conjugated_bond(topology: &TopologyBlock, atom: AtomId) -> bool {
    topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .any(|neighbor| topology.bonds[neighbor.bond.index()].is_conjugated())
}

fn is_potential_nontetrahedral(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom: AtomId,
) -> Result<bool, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION isAtomPotentialNontetrahedralCenter
    // RDKit✔️❌: auto nzdegree = Chirality::detail::getAtomNonzeroDegree(atom);
    // RDKit✔️❌: auto impHDegree = atom->getTotalNumHs();
    // RDKit✔️❌: auto tnzdegree = nzdegree + impHDegree;
    // RDKit✔️❌: if (tnzdegree > 6 || tnzdegree < 2 || (anum < 12 && anum != 4)) return false;
    // RDKit✔️❌: if (chiralType >= Atom::CHI_SQUAREPLANAR &&
    // RDKit✔️❌:     chiralType <= Atom::CHI_OCTAHEDRAL) return true;
    // RDKit✔️❌: if (chiralType == Atom::CHI_UNSPECIFIED && tnzdegree >= 4) return true;
    // END RDKIT CPP FUNCTION isAtomPotentialNontetrahedralCenter
    // Hydrogen composition uses the already source-backed detached helper,
    // whose full-topology validation makes this path slower than the source.
    let value = &topology.atoms[atom.index()];
    let degree = nonzero_degree(topology, atom)? + total_hydrogens(topology, valence, atom, false)?;
    if degree > 6 || degree < 2 || (value.atomic_number() < 12 && value.atomic_number() != 4) {
        return Ok(false);
    }
    Ok(matches!(
        value.chiral_tag(),
        ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral
    ) || (value.chiral_tag() == ChiralTag::Unspecified && degree >= 4))
}

fn is_potential_tetrahedral(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    atom: AtomId,
) -> Result<bool, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION isAtomPotentialTetrahedralCenter
    // RDKit✔️❌: auto nzDegree = getAtomNonzeroDegree(atom);
    // RDKit✔️❌: auto tnzDegree = nzDegree + atom->getTotalNumHs();
    // RDKit✔️❌: if (tnzDegree > 4) return false;
    // RDKit✔️❌: if (nzDegree == 4) return true;
    // RDKit✔️❌: if (nzDegree <= 1) return false;
    // RDKit✔️❌: if (nzDegree < 3 && anum != 15 && anum != 33) return false;
    // RDKit✔️❌: if (anum == 15 || anum == 33) return true;
    // RDKit✔️❌: if (nzDegree == 3 && atom->getTotalNumHs() == 1) {
    // RDKit✔️❌:   if (detail::has_protium_neighbor(mol, atom)) return false;
    // RDKit✔️❌:   return true;
    // RDKit✔️❌: }
    // RDKit✔️❌: if ((anum == 16 || anum == 34) &&
    // RDKit✔️❌:    (explicitValence == 4 || (explicitValence == 3 && charge == 1))) return true;
    // RDKit✔️❌: if (anum == 7 && hybridization == SP3 && !atomHasConjugatedBond(atom) &&
    // RDKit✔️❌:    (isAtomInRingOfSize(idx, 3) || queryIsAtomBridgehead(atom))) return true;
    // END RDKIT CPP FUNCTION isAtomPotentialTetrahedralCenter
    let value = &topology.atoms[atom.index()];
    let degree = nonzero_degree(topology, atom)?;
    let hydrogens = total_hydrogens(topology, valence, atom, false)?;
    if degree + hydrogens > 4 {
        return Ok(false);
    }
    if degree == 4 {
        return Ok(true);
    }
    if degree <= 1 {
        return Ok(false);
    }
    if degree < 3 && !matches!(value.atomic_number(), 15 | 33) {
        return Ok(false);
    }
    if matches!(value.atomic_number(), 15 | 33) {
        return Ok(true);
    }
    if degree == 3 {
        if hydrogens == 1 {
            return Ok(!has_protium_neighbor(topology, atom));
        }
        if matches!(value.atomic_number(), 16 | 34)
            && (valence.explicit_valence[atom.index()] == 4
                || (valence.explicit_valence[atom.index()] == 3 && value.formal_charge() == 1))
        {
            return Ok(true);
        }
        if value.atomic_number() == 7
            && value.hybridization() == Hybridization::Sp3
            && !has_conjugated_bond(topology, atom)
            && (rings.is_atom_in_ring_of_size(atom, 3)
                || is_atom_bridgehead_from_topology(topology, atom.index(), rings) != 0)
        {
            return Ok(true);
        }
    }
    Ok(false)
}

fn is_potential_atom(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    atom: AtomId,
    allow_nontetrahedral: bool,
) -> Result<bool, PotentialStereoError> {
    Ok(is_potential_tetrahedral(topology, valence, rings, atom)?
        || (allow_nontetrahedral && is_potential_nontetrahedral(topology, valence, atom)?))
}

fn atom_info(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom: AtomId,
    allow_nontetrahedral: bool,
) -> Result<PotentialStereoInfo, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION getStereoInfo(const Atom *)
    // RDKit✔️✔️: for (const auto &nbri : mol.getAtomBonds(atom)) {
    // RDKit✔️✔️:   if (bnd->getBondDir() == Bond::UNKNOWN) explicitUnknownStereo = 1;
    // RDKit✔️✔️:   sinfo.controllingAtoms.push_back(bnd->getOtherAtomIdx(atom->getIdx()));
    // RDKit✔️✔️: }
    // RDKit✔️✔️: std::vector<unsigned> origNbrOrder = sinfo.controllingAtoms;
    // RDKit✔️✔️: std::sort(sinfo.controllingAtoms.begin(), sinfo.controllingAtoms.end());
    // RDKit✔️✔️: if (explicitUnknownStereo) sinfo.specified = StereoSpecified::Unknown;
    // RDKit✔️✔️: else if (stereo == CHI_TETRAHEDRAL_CCW || stereo == CHI_TETRAHEDRAL_CW) {
    // RDKit✔️✔️:   unsigned nSwaps = countSwapsToInterconvert(origNbrOrder, controllingAtoms);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getStereoInfo(const Atom *)
    let value = &topology.atoms[atom.index()];
    let original = topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .map(|neighbor| AtomId::new(neighbor.atom_index))
        .collect::<Vec<_>>();
    let mut controls = original.clone();
    controls.sort_unstable();
    let explicit_unknown = value.unknown_stereo()
        || topology
            .adjacency
            .neighbors_of(atom.index())
            .iter()
            .any(|neighbor| {
                let bond = &topology.bonds[neighbor.bond.index()];
                bond.direction() == BondDirection::Unknown || bond.unknown_stereo()
            });
    let mut info = PotentialStereoInfo {
        stereo_type: PotentialStereoType::AtomTetrahedral,
        specified: PotentialStereoSpecified::Unspecified,
        centered_on: PotentialStereoCenter::Atom(atom),
        descriptor: PotentialStereoDescriptor::None,
        permutation: 0,
        controlling_atoms: controls.iter().copied().map(Some).collect(),
    };
    if explicit_unknown {
        info.specified = PotentialStereoSpecified::Unknown;
        return Ok(info);
    }
    if matches!(
        value.chiral_tag(),
        ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
    ) {
        info.specified = PotentialStereoSpecified::Specified;
        let odd = count_swaps_to_interconvert(&original, &controls)? % 2 == 1;
        let clockwise = (value.chiral_tag() == ChiralTag::TetrahedralCw) ^ odd;
        info.descriptor = if clockwise {
            PotentialStereoDescriptor::TetrahedralClockwise
        } else {
            PotentialStereoDescriptor::TetrahedralCounterclockwise
        };
        return Ok(info);
    }
    if allow_nontetrahedral && is_potential_nontetrahedral(topology, valence, atom)? {
        let degree = total_degree(topology, valence, atom)?;
        info.stereo_type = match value.chiral_tag() {
            ChiralTag::SquarePlanar => PotentialStereoType::AtomSquarePlanar,
            ChiralTag::TrigonalBipyramidal => PotentialStereoType::AtomTrigonalBipyramidal,
            ChiralTag::Octahedral => PotentialStereoType::AtomOctahedral,
            ChiralTag::Unspecified if degree == 5 => PotentialStereoType::AtomTrigonalBipyramidal,
            ChiralTag::Unspecified if degree == 6 => PotentialStereoType::AtomOctahedral,
            _ => PotentialStereoType::AtomTetrahedral,
        };
        if let Some(permutation) = value.chiral_permutation() {
            info.permutation = permutation;
            info.specified = if permutation == 0 {
                PotentialStereoSpecified::Unknown
            } else {
                PotentialStereoSpecified::Specified
            };
        }
    }
    Ok(info)
}

fn is_potential_bond(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    bond: &Bond,
) -> Result<bool, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION isBondPotentialStereoBond
    // RDKit✔️❌: if (bond->getBondType() != Bond::DOUBLE) return false;
    // RDKit✔️❌: if (begDegree > 1 && begDegree < 4 && endDegree > 1 && endDegree < 4 &&
    // RDKit✔️❌:     beginAtom->getTotalNumHs(true) < 2 && endAtom->getTotalNumHs(true) < 2) {
    // RDKit✔️❌:   for (const auto &bring : ri->bondRings()) {
    // RDKit✔️❌:     if (bring.size() < minRingSizeForDoubleBondStereo && contains(bring, bidx))
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return true;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION isBondPotentialStereoBond
    if bond.order() != BondOrder::Double {
        return Ok(false);
    }
    let begin_degree = total_degree(topology, valence, bond.begin())?;
    let end_degree = total_degree(topology, valence, bond.end())?;
    if !(2..4).contains(&begin_degree)
        || !(2..4).contains(&end_degree)
        || total_hydrogens(topology, valence, bond.begin(), true)? >= 2
        || total_hydrogens(topology, valence, bond.end(), true)? >= 2
    {
        return Ok(false);
    }
    Ok(!rings
        .bond_ring_sizes(bond.id())
        .into_iter()
        .any(|size| size < 8))
}

fn bond_info(
    topology: &TopologyBlock,
    bond: BondId,
) -> Result<PotentialStereoInfo, PotentialStereoError> {
    let source = double_bond_stereo_info(topology, bond)?;
    Ok(PotentialStereoInfo {
        stereo_type: PotentialStereoType::BondDouble,
        specified: match source.specified {
            DoubleBondStereoSpecified::Unspecified => PotentialStereoSpecified::Unspecified,
            DoubleBondStereoSpecified::Specified => PotentialStereoSpecified::Specified,
            DoubleBondStereoSpecified::Unknown => PotentialStereoSpecified::Unknown,
        },
        centered_on: PotentialStereoCenter::Bond(bond),
        descriptor: match source.descriptor {
            Some(DoubleBondStereoDescriptor::Cis) => PotentialStereoDescriptor::BondCis,
            Some(DoubleBondStereoDescriptor::Trans) => PotentialStereoDescriptor::BondTrans,
            None => PotentialStereoDescriptor::None,
        },
        permutation: 0,
        controlling_atoms: source
            .controlling_atoms
            .into_iter()
            .map(|control| match control {
                DoubleBondControl::Atom(atom) => Some(atom),
                DoubleBondControl::Implicit => None,
            })
            .collect(),
    })
}

fn atom_symbol(atom: &Atom) -> String {
    format!(
        "{}{}{}",
        atom.isotope().unwrap_or(0),
        atom.element().symbol(),
        atom.formal_charge()
    )
}

fn bond_symbol(bond: &Bond) -> &'static str {
    // BEGIN RDKIT CPP FUNCTION getBondSymbol
    // RDKit✔️✔️: if (bond->getIsAromatic()) res = ":";
    // RDKit✔️✔️: else switch (bond->getBondType()) {
    // RDKit✔️✔️: case SINGLE: res = "-"; break; case DOUBLE: res = "="; break;
    // RDKit✔️✔️: case TRIPLE: res = "#"; break; case AROMATIC: res = ":"; break;
    // RDKit✔️✔️: default: res = "?"; break; }
    // END RDKIT CPP FUNCTION getBondSymbol
    if bond.is_aromatic() {
        ":"
    } else {
        match bond.order() {
            BondOrder::Single => "-",
            BondOrder::Double => "=",
            BondOrder::Triple => "#",
            BondOrder::Aromatic => ":",
            _ => "?",
        }
    }
}

fn connectivity_ranks(
    topology: &TopologyBlock,
    atom_symbols: &[String],
    bond_symbols: &[String],
) -> Result<Vec<u32>, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION rankFragmentAtoms restricted call from runCleanup
    // RDKit✔️✔️: Canon::rankFragmentAtoms(mol, aranks, atomsInPlay, bondsInPlay,
    // RDKit✔️✔️:   &atomSymbols, &bondSymbols, false, false, false, false, false, false);
    // END RDKIT CPP FUNCTION rankFragmentAtoms restricted call from runCleanup
    // With every atom and bond in play and all six optional dimensions false,
    // the source comparison is the stable equitable refinement of the supplied
    // atom label and sorted incident (bond label, neighbor class) multiset.
    let count = topology.atoms.len();
    if count == 0 {
        return Ok(Vec::new());
    }
    let mut ranks = vec![0; count];
    for iteration in 0..=count {
        let keys = (0..count)
            .map(|index| {
                let mut neighbors = topology
                    .adjacency
                    .neighbors_of(index)
                    .iter()
                    .map(|neighbor| {
                        (
                            bond_symbols[neighbor.bond.index()].clone(),
                            ranks[neighbor.atom_index],
                        )
                    })
                    .collect::<Vec<_>>();
                neighbors.sort_unstable();
                (ranks[index], atom_symbols[index].clone(), neighbors)
            })
            .collect::<Vec<_>>();
        let unique = keys.iter().cloned().collect::<BTreeSet<_>>();
        let classes = unique
            .into_iter()
            .enumerate()
            .map(|(rank, key)| (key, rank as u32))
            .collect::<BTreeMap<_, _>>();
        let next = keys.iter().map(|key| classes[key]).collect::<Vec<_>>();
        if iteration > 0 && next == ranks {
            return Ok(next);
        }
        ranks = next;
    }
    Err(PotentialStereoError::RefinementDidNotConverge {
        iterations: count + 1,
    })
}

fn initialize_atoms(
    topology: &mut TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    params: &PotentialStereoParams,
    known: &mut [bool],
    possible: &mut [bool],
    symbols: &mut [String],
) -> Result<(), PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION initAtomInfo
    // RDKit✔️❌: atomSymbols[aidx] = getAtomCompareSymbol(*atom);
    // RDKit✔️❌: if (detail::isAtomPotentialStereoAtom(atom, allowNontetrahedralStereo)) {
    // RDKit✔️❌:   auto sinfo = detail::getStereoInfo(atom);
    // RDKit✔️❌:   switch (sinfo.specified) {
    // RDKit✔️❌:   case Unknown: knownAtoms.set(aidx); atomSymbols[aidx] += std::to_string(aidx); break;
    // RDKit✔️❌:   case Chirality::StereoSpecified::Specified:
    // RDKit✔️❌:     knownAtoms.set(aidx);
    // RDKit✔️❌:     if (sinfo.descriptor == StereoDescriptor::Tet_CCW) {
    // RDKit✔️❌:       atomSymbols[aidx] += "_CCW";
    // RDKit✔️❌:     } else if (sinfo.descriptor == StereoDescriptor::Tet_CW) {
    // RDKit✔️❌:       atomSymbols[aidx] += "_CW";
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       atomSymbols[aidx] += "_STEREO";
    // RDKit✔️❌:     }
    // RDKit✔️❌:     break;
    // RDKit✔️❌:   case Unspecified: if (flagPossible) possibleAtoms.set(aidx); break;
    // RDKit✔️❌:   }
    // RDKit✔️❌: } else if (cleanIt) atom->setChiralTag(CHI_UNSPECIFIED);
    // END RDKIT CPP FUNCTION initAtomInfo
    for index in 0..topology.atoms.len() {
        let atom = AtomId::new(index);
        symbols[index] = atom_symbol(&topology.atoms[index]);
        if is_potential_atom(topology, valence, rings, atom, params.allow_nontetrahedral)? {
            let info = atom_info(topology, valence, atom, params.allow_nontetrahedral)?;
            match info.specified {
                PotentialStereoSpecified::Unknown => {
                    known[index] = true;
                    symbols[index].push_str(&index.to_string());
                }
                PotentialStereoSpecified::Specified => {
                    known[index] = true;
                    symbols[index].push_str(match info.descriptor {
                        PotentialStereoDescriptor::TetrahedralCounterclockwise => "_CCW",
                        PotentialStereoDescriptor::TetrahedralClockwise => "_CW",
                        _ => "_STEREO",
                    });
                }
                PotentialStereoSpecified::Unspecified if params.flag_possible => {
                    possible[index] = true;
                    if !params.clean {
                        symbols[index].push('_');
                        symbols[index].push_str(&index.to_string());
                    }
                }
                PotentialStereoSpecified::Unspecified => {}
            }
        } else if params.clean {
            topology.atoms[index].set_chiral_tag(ChiralTag::Unspecified);
        }
    }
    Ok(())
}

fn clear_bond_stereo_value(bond: &mut Bond) -> Result<(), PotentialStereoError> {
    bond.set_stereo(BondStereo::None)?;
    bond.set_stereo_atoms(None);
    Ok(())
}

fn initialize_bonds(
    topology: &mut TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    params: &PotentialStereoParams,
    known: &mut [bool],
    possible: &mut [bool],
    symbols: &mut [String],
) -> Result<(), PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION initBondInfo
    // RDKit✔️❌: bondSymbols[bidx] = getBondSymbol(bond);
    // RDKit✔️❌: if (detail::isBondPotentialStereoBond(bond)) {
    // RDKit✔️❌:   auto sinfo = detail::getStereoInfo(bond);
    // RDKit✔️❌:   switch (sinfo.specified) {
    // RDKit✔️❌:   case Unknown: knownBonds.set(bidx); bondSymbols[bidx] += "_" + to_string(bidx); break;
    // RDKit✔️❌:   case Chirality::StereoSpecified::Specified:
    // RDKit✔️❌:     knownBonds.set(bidx);
    // RDKit✔️❌:     if (sinfo.descriptor == StereoDescriptor::Bond_Cis) {
    // RDKit✔️❌:       bondSymbols[bidx] += "_cis";
    // RDKit✔️❌:     } else if (sinfo.descriptor == StereoDescriptor::Bond_Trans) {
    // RDKit✔️❌:       bondSymbols[bidx] += "_trans";
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       bondSymbols[bidx] += "_STEREO";
    // RDKit✔️❌:     }
    // RDKit✔️❌:     break;
    // RDKit✔️❌:   case Unspecified: if (flagPossible) possibleBonds.set(bidx); break;
    // RDKit✔️❌:   }
    // RDKit✔️❌: } else if (cleanIt) bond->setStereo(STEREONONE);
    // END RDKIT CPP FUNCTION initBondInfo
    for index in 0..topology.bonds.len() {
        symbols[index] = bond_symbol(&topology.bonds[index]).to_owned();
        let candidate = {
            let bond = &topology.bonds[index];
            is_potential_bond(topology, valence, rings, bond)?
        };
        if candidate {
            let info = bond_info(topology, BondId::new(index))?;
            match info.specified {
                PotentialStereoSpecified::Unknown => {
                    known[index] = true;
                    symbols[index].push('_');
                    symbols[index].push_str(&index.to_string());
                }
                PotentialStereoSpecified::Specified => {
                    known[index] = true;
                    symbols[index].push_str(match info.descriptor {
                        PotentialStereoDescriptor::BondCis => "_cis",
                        PotentialStereoDescriptor::BondTrans => "_trans",
                        _ => "_STEREO",
                    });
                }
                PotentialStereoSpecified::Unspecified if params.flag_possible => {
                    possible[index] = true;
                    if !params.clean {
                        symbols[index].push('_');
                        symbols[index].push_str(&index.to_string());
                    }
                }
                PotentialStereoSpecified::Unspecified => {}
            }
        } else if params.clean {
            clear_bond_stereo_value(&mut topology.bonds[index])?;
        }
    }
    Ok(())
}

fn bond_between(topology: &TopologyBlock, left: AtomId, right: AtomId) -> Option<BondId> {
    topology
        .adjacency
        .neighbors_of(left.index())
        .iter()
        .find_map(|neighbor| (neighbor.atom_index == right.index()).then_some(neighbor.bond))
}

fn flag_ring_stereo(
    topology: &TopologyBlock,
    rings: &RingInfo,
    possible_ring_atoms: &mut [usize],
    possible_ring_bonds: &mut [usize],
    known_atoms: &[bool],
    possible_atoms: Option<&[bool]>,
    known_bonds: &[bool],
    possible_bonds: Option<&[bool]>,
) {
    // BEGIN RDKIT CPP FUNCTION flagRingStereo
    // RDKit✔️✔️: for (unsigned int ridx = 0; ridx < ringInfo->atomRings().size(); ++ridx) {
    // RDKit✔️✔️:   for (unsigned int ai = 0; ai < sz; ++ai) {
    // RDKit✔️✔️:     for (unsigned int ringDivisor : {2, 3}) {
    // RDKit✔️✔️:       bool ringIsMultipleOfDivisor = ((sz % ringDivisor) == 0);
    // RDKit✔️✔️:       auto incrementSize = sz / ringDivisor;
    // RDKit✔️✔️:       if (ringIsMultipleOfDivisor) {
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (ringInfo->numAtomRings(aidx) > 1) {
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (nHere > 1) {
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION flagRingStereo
    let active_atom =
        |index: usize| known_atoms[index] || possible_atoms.is_some_and(|rows| rows[index]);
    let active_bond =
        |index: usize| known_bonds[index] || possible_bonds.is_some_and(|rows| rows[index]);
    for (ring_index, ring) in rings.atom_rings().iter().enumerate() {
        let ring_bonds = &rings.bond_rings()[ring_index];
        let size = ring.len();
        if size == 0 {
            continue;
        }
        let half_size = size / 2 + usize::from(size % 2 != 0);
        let mut marked = vec![false; topology.atoms.len()];
        let mut count_here = 0usize;
        for position in 0..size {
            let atom = ring[position];
            if !active_atom(atom.index()) {
                continue;
            }
            for divisor in [2usize, 3] {
                if size % divisor != 0 {
                    continue;
                }
                let increment = size / divisor;
                let mut by_bond = 0usize;
                let mut by_atom = 0usize;
                let mut offset = increment;
                while offset < size {
                    let other = ring[(position + offset) % size];
                    let has_external =
                        topology
                            .adjacency
                            .neighbors_of(other.index())
                            .iter()
                            .any(|neighbor| {
                                active_bond(neighbor.bond.index())
                                    && !ring_bonds.contains(&neighbor.bond)
                            });
                    if has_external {
                        by_bond += 1;
                    } else if by_bond == 0 && active_atom(other.index()) {
                        by_atom += 1;
                    }
                    offset += increment;
                }
                if by_bond == divisor - 1 || by_atom == divisor - 1 {
                    count_here += 1 + by_bond;
                    let mut offset = 0;
                    while offset < size {
                        marked[ring[(position + offset) % size].index()] = true;
                        offset += increment;
                    }
                }
            }
            if rings.num_atom_rings(atom) > 1 {
                let mut previous = atom;
                for step in 1..=half_size {
                    let other = ring[(position + step) % size];
                    let Some(edge) = bond_between(topology, previous, other) else {
                        break;
                    };
                    if rings.num_bond_rings(edge) < 2 {
                        break;
                    }
                    if active_atom(other.index()) {
                        count_here += 2;
                        marked[atom.index()] = true;
                        marked[other.index()] = true;
                        break;
                    }
                    previous = other;
                }
            }
        }
        if count_here > 1 {
            for atom in ring {
                if marked[atom.index()] {
                    possible_ring_atoms[atom.index()] += 1;
                }
            }
            for bond in ring_bonds {
                possible_ring_bonds[bond.index()] += 1;
            }
        }
    }
}

fn controlling_atoms_are_duplicates(
    topology: &TopologyBlock,
    rings: &RingInfo,
    bond: BondId,
    left: AtomId,
    right: AtomId,
    ranks: &[u32],
    possible_atoms: &[bool],
    known_atoms: &[bool],
    possible_bonds: &[bool],
    known_bonds: &[bool],
) -> bool {
    // BEGIN RDKIT CPP FUNCTION areStereobondControllingAtomsDupes
    // RDKit✔️✔️: if (atomRanks[controllingAtom1] != atomRanks[controllingAtom2]) return false;
    // RDKit✔️✔️: if (atomRanks[controllingAtom1] != atomRanks[controllingAtom2]) {
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION areStereobondControllingAtomsDupes
    if ranks[left.index()] != ranks[right.index()] {
        return false;
    }
    let left_members = rings.atom_members(left);
    let right_members = rings.atom_members(right);
    let mut li = 0;
    let mut ri = 0;
    while li < left_members.len() && ri < right_members.len() {
        match left_members[li].cmp(&right_members[ri]) {
            std::cmp::Ordering::Less => li += 1,
            std::cmp::Ordering::Greater => ri += 1,
            std::cmp::Ordering::Equal => {
                let ring = &rings.atom_rings()[left_members[li]];
                li += 1;
                ri += 1;
                if ring.len() % 2 != 0 {
                    continue;
                }
                for endpoint in [
                    topology.bonds[bond.index()].begin(),
                    topology.bonds[bond.index()].end(),
                ] {
                    let Some(position) = ring.iter().position(|atom| *atom == endpoint) else {
                        continue;
                    };
                    let opposite = ring[(position + ring.len() / 2) % ring.len()];
                    if possible_atoms[opposite.index()] || known_atoms[opposite.index()] {
                        return false;
                    }
                    if graph_degree(topology, opposite) == 3 {
                        for neighbor in topology.adjacency.neighbors_of(opposite.index()) {
                            if !ring.contains(&AtomId::new(neighbor.atom_index))
                                && (possible_bonds[neighbor.bond.index()]
                                    || known_bonds[neighbor.bond.index()])
                            {
                                return false;
                            }
                        }
                    }
                }
            }
        }
    }
    true
}

#[allow(clippy::too_many_arguments)]
fn update_atoms(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    ranks: &[u32],
    symbols: &mut [String],
    possible: &mut [bool],
    known: &[bool],
    fixed: &mut [bool],
    possible_ring_atoms: &mut [usize],
    possible_ring_bonds: &mut [usize],
    rings: &RingInfo,
    allow_nontetrahedral: bool,
    output: &mut Vec<PotentialStereoInfo>,
) -> Result<bool, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION updateAtoms
    // RDKit✔️✔️: if (knownAtoms[aidx] || possibleAtoms[aidx]) {
    // RDKit✔️✔️:   auto sinfo = detail::getStereoInfo(atom);
    // RDKit✔️✔️:   if (fixedAtoms[aidx]) sinfos.push_back(std::move(sinfo));
    // RDKit✔️✔️:   if (fixedAtoms[aidx]) {
    // RDKit✔️✔️:     sinfos.push_back(std::move(sinfo));
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     std::vector<unsigned int> nbrs;
    // RDKit✔️✔️:     nbrs.reserve(sinfo.controllingAtoms.size());
    // RDKit✔️✔️:     bool haveADupe = false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION updateAtoms
    let mut another = false;
    for index in 0..topology.atoms.len() {
        if !known[index] && !possible[index] {
            continue;
        }
        let atom = AtomId::new(index);
        let mut info = atom_info(topology, valence, atom, allow_nontetrahedral)?;
        if fixed[index] {
            output.push(info);
            continue;
        }
        let mut neighbor_ranks = Vec::new();
        let mut duplicate = false;
        if info.stereo_type == PotentialStereoType::AtomTetrahedral {
            for neighbor in info.controlling_atoms.iter().flatten().copied() {
                let rank = ranks[neighbor.index()];
                if neighbor_ranks.contains(&rank) {
                    if possible_ring_atoms[index] > 0 {
                        let transmitting = bond_between(topology, atom, neighbor)
                            .is_some_and(|bond| possible_ring_bonds[bond.index()] > 0);
                        if !transmitting {
                            duplicate = true;
                            break;
                        }
                    } else {
                        duplicate = true;
                        break;
                    }
                } else {
                    neighbor_ranks.push(rank);
                }
            }
        }
        if !duplicate {
            let mut next_symbol = symbols[index].clone();
            if !possible[index] {
                let mut sorted = neighbor_ranks.clone();
                sorted.sort_unstable();
                if info.stereo_type == PotentialStereoType::AtomTetrahedral
                    && count_swaps_to_interconvert(&neighbor_ranks, &sorted)? % 2 == 1
                {
                    info.descriptor = match info.descriptor {
                        PotentialStereoDescriptor::TetrahedralClockwise => {
                            PotentialStereoDescriptor::TetrahedralCounterclockwise
                        }
                        PotentialStereoDescriptor::TetrahedralCounterclockwise => {
                            PotentialStereoDescriptor::TetrahedralClockwise
                        }
                        descriptor => descriptor,
                    };
                }
                next_symbol = atom_symbol(&topology.atoms[index]);
                next_symbol.push_str(match info.descriptor {
                    PotentialStereoDescriptor::TetrahedralClockwise => "_CW",
                    PotentialStereoDescriptor::TetrahedralCounterclockwise => "_CCW",
                    _ => "",
                });
                fixed[index] = true;
            }
            if symbols[index] != next_symbol {
                symbols[index] = next_symbol;
                another = true;
            }
            output.push(info);
        } else {
            another |= possible[index];
            possible[index] = false;
            symbols[index] = atom_symbol(&topology.atoms[index]);
            if possible_ring_atoms[index] > 0 {
                possible_ring_atoms[index] = 0;
                another = true;
                for (ring_index, ring) in rings.atom_rings().iter().enumerate() {
                    let mut remaining = 0usize;
                    for ring_atom in ring {
                        fixed[ring_atom.index()] = false;
                        remaining += usize::from(possible_ring_atoms[ring_atom.index()] > 0);
                    }
                    if remaining <= 1 {
                        if remaining == 1 {
                            if let Some(last) = ring
                                .iter()
                                .find(|ring_atom| possible_ring_atoms[ring_atom.index()] > 0)
                            {
                                possible_ring_atoms[last.index()] -= 1;
                            }
                        }
                        for ring_bond in &rings.bond_rings()[ring_index] {
                            if possible_ring_bonds[ring_bond.index()] > 0 {
                                possible_ring_bonds[ring_bond.index()] -= 1;
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(another)
}

fn swap_bond_descriptor(info: &mut PotentialStereoInfo) {
    info.descriptor = match info.descriptor {
        PotentialStereoDescriptor::BondCis => PotentialStereoDescriptor::BondTrans,
        PotentialStereoDescriptor::BondTrans => PotentialStereoDescriptor::BondCis,
        descriptor => descriptor,
    };
}

#[allow(clippy::too_many_arguments)]
fn update_bonds(
    topology: &TopologyBlock,
    ranks: &[u32],
    symbols: &mut [String],
    possible_atoms: &[bool],
    possible_bonds: &mut [bool],
    known_atoms: &[bool],
    known_bonds: &[bool],
    fixed_bonds: &mut [bool],
    rings: &RingInfo,
    output: &mut Vec<PotentialStereoInfo>,
) -> Result<bool, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION updateBonds
    // RDKit✔️✔️: if (knownBonds[bidx] || possibleBonds[bidx]) {
    // RDKit✔️✔️:   auto sinfo = detail::getStereoInfo(bond);
    // RDKit✔️✔️:   if (both controlling slots missing on either end) fixedBonds.set(bidx);
    // RDKit✔️✔️:   if (!fixedBonds[bidx]) {
    // RDKit✔️✔️:     if (sinfo.controllingAtoms[0] != Atom::NOATOM &&
    // RDKit✔️✔️:         sinfo.controllingAtoms[1] != Atom::NOATOM) {
    // RDKit✔️✔️:     if (haveADupe && possibleBonds[bidx]) possibleBonds[bidx] = 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION updateBonds
    let mut another = false;
    for index in 0..topology.bonds.len() {
        if !known_bonds[index] && !possible_bonds[index] {
            continue;
        }
        let bond_id = BondId::new(index);
        let mut info = bond_info(topology, bond_id)?;
        if info.controlling_atoms.len() != 4 {
            return Err(PotentialStereoError::InvalidStereoReferences {
                bond: bond_id,
                reason: "potential stereo bond must have four controlling slots",
            });
        }
        if (info.controlling_atoms[0].is_none() && info.controlling_atoms[1].is_none())
            || (info.controlling_atoms[2].is_none() && info.controlling_atoms[3].is_none())
        {
            if info.specified == PotentialStereoSpecified::Specified {
                return Err(PotentialStereoError::InvalidStereoReferences {
                    bond: bond_id,
                    reason: "specified bond has no controlling atom on one endpoint",
                });
            }
            fixed_bonds[index] = true;
        }
        if fixed_bonds[index] {
            output.push(info);
            continue;
        }
        let mut duplicate = false;
        let mut needs_swap = false;
        for offset in [0usize, 2] {
            if let (Some(left), Some(right)) = (
                info.controlling_atoms[offset],
                info.controlling_atoms[offset + 1],
            ) {
                if controlling_atoms_are_duplicates(
                    topology,
                    rings,
                    bond_id,
                    left,
                    right,
                    ranks,
                    possible_atoms,
                    known_atoms,
                    possible_bonds,
                    known_bonds,
                ) {
                    duplicate = true;
                } else if ranks[left.index()] < ranks[right.index()] {
                    info.controlling_atoms.swap(offset, offset + 1);
                    needs_swap = !needs_swap;
                }
            }
        }
        if !duplicate {
            if needs_swap {
                swap_bond_descriptor(&mut info);
            }
            let mut next_symbol = symbols[index].clone();
            match (info.specified, info.descriptor) {
                (PotentialStereoSpecified::Specified, PotentialStereoDescriptor::BondCis) => {
                    next_symbol.push_str("_cis")
                }
                (PotentialStereoSpecified::Specified, PotentialStereoDescriptor::BondTrans) => {
                    next_symbol.push_str("_trans")
                }
                (PotentialStereoSpecified::Unknown, _) => next_symbol.push_str("_unk"),
                _ => {}
            }
            if symbols[index] != next_symbol {
                symbols[index] = next_symbol;
                another = true;
            }
            if !possible_bonds[index] {
                fixed_bonds[index] = true;
            }
            output.push(info);
        } else if possible_bonds[index] {
            possible_bonds[index] = false;
            symbols[index] = bond_symbol(&topology.bonds[index]).to_owned();
            another = true;
        }
    }
    Ok(another)
}

fn clean_invalid_stereo(
    topology: &mut TopologyBlock,
    fixed_atoms: &[bool],
    known_atoms: &[bool],
    fixed_bonds: &[bool],
    known_bonds: &[bool],
) -> Result<(), PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION cleanMolStereo
    // RDKit✔️✔️: if (!fixedAtoms[i] && knownAtoms[i]) {
    // RDKit✔️✔️:   if (atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:       atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (!fixedBonds[i] && knownBonds[i]) {
    // RDKit✔️✔️:   bond->setStereo(STEREONONE); bond->setBondDir(NONE);
    // RDKit✔️✔️:   bond->getStereoAtoms().clear();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (removedStereo) {
    // END RDKIT CPP FUNCTION cleanMolStereo
    let mut wedge_centers = Vec::new();
    for index in 0..topology.atoms.len() {
        if fixed_atoms[index] || !known_atoms[index] {
            continue;
        }
        match topology.atoms[index].chiral_tag() {
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
                topology.atoms[index].set_chiral_tag(ChiralTag::Unspecified);
                wedge_centers.push(AtomId::new(index));
            }
            ChiralTag::Tetrahedral
            | ChiralTag::SquarePlanar
            | ChiralTag::TrigonalBipyramidal
            | ChiralTag::Octahedral => {
                topology.atoms[index].set_chiral_permutation(Some(0));
            }
            _ => {}
        }
    }
    for atom in wedge_centers {
        let bonds = topology
            .adjacency
            .neighbors_of(atom.index())
            .iter()
            .map(|neighbor| neighbor.bond)
            .collect::<Vec<_>>();
        for bond in bonds {
            if matches!(
                topology.bonds[bond.index()].direction(),
                BondDirection::BeginDash | BondDirection::BeginWedge
            ) {
                topology.bonds[bond.index()].set_direction(BondDirection::None);
            }
        }
    }
    let mut removed_bond_stereo = false;
    for index in 0..topology.bonds.len() {
        if !fixed_bonds[index] && known_bonds[index] {
            clear_bond_stereo_value(&mut topology.bonds[index])?;
            topology.bonds[index].set_direction(BondDirection::None);
            removed_bond_stereo = true;
        }
    }
    if removed_bond_stereo {
        for index in 0..topology.bonds.len() {
            if !matches!(
                topology.bonds[index].direction(),
                BondDirection::EndDownRight | BondDirection::EndUpRight
            ) {
                continue;
            }
            let endpoints = [topology.bonds[index].begin(), topology.bonds[index].end()];
            let mut direction_ok = false;
            for endpoint in endpoints {
                direction_ok = topology
                    .adjacency
                    .neighbors_of(endpoint.index())
                    .iter()
                    .any(|neighbor| {
                        neighbor.bond.index() != index
                            && topology.bonds[neighbor.bond.index()].stereo() != BondStereo::None
                    });
                if !direction_ok {
                    topology.bonds[index].set_direction(BondDirection::None);
                }
            }
        }
    }
    Ok(())
}

fn atom_is_candidate_for_ring_stereochemistry(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    ranks: &[u32],
    atom: AtomId,
) -> Result<bool, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION atomIsCandidateForRingStereochem
    // RDKit✔️❌: if (ringInfo->isInitialized() && ringInfo->numAtomRings(atom->getIdx())) {
    // RDKit✔️❌:   if (atom->getAtomicNum() == 7 && atom->getTotalDegree() == 3 &&
    // RDKit✔️❌:       !ringInfo->isAtomInRingOfSize(atom->getIdx(), 3) &&
    // RDKit✔️❌:       !queryIsAtomBridgehead(atom)) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   std::vector<const Atom *> nonRingNbrs;
    // RDKit✔️❌:   std::vector<const Atom *> ringNbrs;
    // RDKit✔️❌:   std::set<unsigned int> ringNbrRanks;
    // RDKit✔️❌:   for (const auto bond : mol.atomBonds(atom)) {
    // RDKit✔️❌:     if (!ringInfo->numBondRings(bond->getIdx())) {
    // RDKit✔️❌:       nonRingNbrs.push_back(bond->getOtherAtom(atom));
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       const Atom *nbr = bond->getOtherAtom(atom);
    // RDKit✔️❌:       ringNbrs.push_back(nbr);
    // RDKit✔️❌:       ringNbrRanks.insert(atomRanks[nbr->getIdx()]);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION atomIsCandidateForRingStereochem
    // The detached helper repeats validated total-H lookup for each candidate,
    // so behavior is reproduced while the performance axis remains lower.
    if !rings.is_initialized() || rings.num_atom_rings(atom) == 0 {
        return Ok(false);
    }
    let value = &topology.atoms[atom.index()];
    if value.atomic_number() == 7
        && total_degree(topology, valence, atom)? == 3
        && !rings.is_atom_in_ring_of_size(atom, 3)
        && is_atom_bridgehead_from_topology(topology, atom.index(), rings) == 0
    {
        return Ok(false);
    }
    let mut non_ring_neighbors = Vec::new();
    let mut ring_neighbor_count = 0usize;
    let mut ring_neighbor_ranks = BTreeSet::new();
    for neighbor in topology.adjacency.neighbors_of(atom.index()) {
        let neighbor_atom = AtomId::new(neighbor.atom_index);
        if rings.num_bond_rings(neighbor.bond) == 0 {
            non_ring_neighbors.push(neighbor_atom);
        } else {
            ring_neighbor_count += 1;
            ring_neighbor_ranks.insert(ranks[neighbor.atom_index]);
        }
    }
    Ok(match non_ring_neighbors.as_slice() {
        [left, right] => {
            ranks[left.index()] != ranks[right.index()]
                && ring_neighbor_count != ring_neighbor_ranks.len()
        }
        [_] => ring_neighbor_count > ring_neighbor_ranks.len(),
        [] => {
            (ring_neighbor_count == 4 && ring_neighbor_ranks.len() == 3)
                || (ring_neighbor_count == 3 && ring_neighbor_ranks.len() == 2)
        }
        _ => false,
    })
}

pub(crate) fn special_ring_relations(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    ranks: &[u32],
) -> Result<Vec<RingStereoRelation>, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION findChiralAtomSpecialCases
    // RDKit✔️✔️: for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:   if (atomsSeen[atom->getIdx()]) {
    // RDKit✔️✔️:     continue;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (atom->getChiralTag() == Atom::CHI_UNSPECIFIED ||
    // RDKit✔️✔️:       atom->hasProp(common_properties::_CIPCode) ||
    // RDKit✔️✔️:       !mol.getRingInfo()->numAtomRings(atom->getIdx()) ||
    // RDKit✔️✔️:       !atomIsCandidateForRingStereochem(mol, atom, atomRanks)) {
    // RDKit✔️✔️:     continue;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // do a BFS from this ring atom along ring bonds and find other
    // RDKit✔️✔️:   // stereochemistry candidates.
    // END RDKIT CPP FUNCTION findChiralAtomSpecialCases
    let mut candidates = vec![false; topology.atoms.len()];
    for index in 0..topology.atoms.len() {
        let atom = AtomId::new(index);
        let value = &topology.atoms[index];
        candidates[index] = value.chiral_tag() != ChiralTag::Unspecified
            && value.prop("_CIPCode").is_none()
            && rings.num_atom_rings(atom) > 0
            && atom_is_candidate_for_ring_stereochemistry(topology, valence, rings, ranks, atom)?;
    }

    let mut seen = vec![false; topology.atoms.len()];
    let mut relations = BTreeSet::new();
    for start_index in 0..topology.atoms.len() {
        if seen[start_index] || !candidates[start_index] {
            continue;
        }
        let mut queue = VecDeque::from([AtomId::new(start_index)]);
        let mut component = Vec::new();
        seen[start_index] = true;
        while let Some(atom) = queue.pop_front() {
            if candidates[atom.index()] {
                component.push(atom);
            }
            for neighbor in topology.adjacency.neighbors_of(atom.index()) {
                if rings.num_bond_rings(neighbor.bond) == 0 || seen[neighbor.atom_index] {
                    continue;
                }
                seen[neighbor.atom_index] = true;
                queue.push_back(AtomId::new(neighbor.atom_index));
            }
        }
        component.sort_unstable();
        for left_index in 0..component.len() {
            for right_index in (left_index + 1)..component.len() {
                let left = component[left_index];
                let right = component[right_index];
                let same = topology.atoms[left.index()].chiral_tag()
                    == topology.atoms[right.index()].chiral_tag();
                relations.insert(RingStereoRelation {
                    atom: left,
                    other: right,
                    same_orientation: same,
                });
                relations.insert(RingStereoRelation {
                    atom: right,
                    other: left,
                    same_orientation: same,
                });
            }
        }
    }
    Ok(relations.into_iter().collect())
}

#[allow(clippy::too_many_arguments)]
fn refinement_loop(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    allow_nontetrahedral: bool,
    atom_symbols: &mut [String],
    bond_symbols: &mut [String],
    possible_atoms: &mut [bool],
    known_atoms: &[bool],
    fixed_atoms: &mut [bool],
    possible_bonds: &mut [bool],
    known_bonds: &[bool],
    fixed_bonds: &mut [bool],
    possible_ring_atoms: &mut [usize],
    possible_ring_bonds: &mut [usize],
) -> Result<(Vec<PotentialStereoInfo>, Vec<u32>), PotentialStereoError> {
    let limit = topology.atoms.len() + topology.bonds.len() + 2;
    let mut ranks = vec![0; topology.atoms.len()];
    let mut result = Vec::new();
    for iteration in 0..limit {
        result.clear();
        ranks = connectivity_ranks(topology, atom_symbols, bond_symbols)?;
        let mut another = update_atoms(
            topology,
            valence,
            &ranks,
            atom_symbols,
            possible_atoms,
            known_atoms,
            fixed_atoms,
            possible_ring_atoms,
            possible_ring_bonds,
            rings,
            allow_nontetrahedral,
            &mut result,
        )?;
        another |= update_bonds(
            topology,
            &ranks,
            bond_symbols,
            possible_atoms,
            possible_bonds,
            known_atoms,
            known_bonds,
            fixed_bonds,
            rings,
            &mut result,
        )?;
        if !another {
            return Ok((result, ranks));
        }
        if iteration + 1 == limit {
            return Err(PotentialStereoError::RefinementDidNotConverge { iterations: limit });
        }
    }
    unreachable!("nonzero refinement limit always returns from the loop")
}

pub fn potential_stereo(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    params: &PotentialStereoParams,
) -> Result<PotentialStereoAssignment, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION findPotentialStereo
    // RDKit✔️❌: std::vector<StereoInfo> findPotentialStereo(ROMol &mol, bool cleanIt,
    // RDKit✔️❌:                                             bool findPossible) {
    // RDKit✔️❌:   if (!mol.getRingInfo()->isSymmSssr()) {
    // RDKit✔️❌:     MolOps::symmetrizeSSSR(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (mol.needsUpdatePropertyCache()) {
    // RDKit✔️❌:     mol.updatePropertyCache(false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   std::vector<StereoInfo> res = runCleanup(mol, findPossible, cleanIt);
    // RDKit✔️❌:   mol.setProp("_potentialStereo", res, true);
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION findPotentialStereo
    // The caller supplies the already materialized valence and symmetric-ring
    // assignments. Repeated detached total-H validation and the private rank
    // implementation make this port intentionally slower than the C++ path.
    validate_inputs(topology, valence, rings)?;
    let mut working = topology.clone();
    let atom_count = working.atoms.len();
    let bond_count = working.bonds.len();
    let mut known_atoms = vec![false; atom_count];
    let mut possible_atoms = vec![false; atom_count];
    let mut atom_symbols = vec![String::new(); atom_count];
    initialize_atoms(
        &mut working,
        valence,
        rings,
        params,
        &mut known_atoms,
        &mut possible_atoms,
        &mut atom_symbols,
    )?;
    let mut known_bonds = vec![false; bond_count];
    let mut possible_bonds = vec![false; bond_count];
    let mut bond_symbols = vec![String::new(); bond_count];
    initialize_bonds(
        &mut working,
        valence,
        rings,
        params,
        &mut known_bonds,
        &mut possible_bonds,
        &mut bond_symbols,
    )?;
    let original_possible_atoms = possible_atoms.clone();
    let original_possible_bonds = possible_bonds.clone();
    let mut possible_ring_atoms = vec![0usize; atom_count];
    let mut possible_ring_bonds = vec![0usize; bond_count];
    let possible_atom_view = if params.clean {
        None
    } else {
        Some(possible_atoms.as_slice())
    };
    let possible_bond_view = if params.clean {
        None
    } else {
        Some(possible_bonds.as_slice())
    };
    flag_ring_stereo(
        &working,
        rings,
        &mut possible_ring_atoms,
        &mut possible_ring_bonds,
        &known_atoms,
        possible_atom_view,
        &known_bonds,
        possible_bond_view,
    );
    let mut fixed_atoms = vec![false; atom_count];
    let mut fixed_bonds = vec![false; bond_count];
    let (mut stereo, mut ranks) = refinement_loop(
        &working,
        valence,
        rings,
        params.allow_nontetrahedral,
        &mut atom_symbols,
        &mut bond_symbols,
        &mut possible_atoms,
        &known_atoms,
        &mut fixed_atoms,
        &mut possible_bonds,
        &known_bonds,
        &mut fixed_bonds,
        &mut possible_ring_atoms,
        &mut possible_ring_bonds,
    )?;
    if params.clean {
        clean_invalid_stereo(
            &mut working,
            &fixed_atoms,
            &known_atoms,
            &fixed_bonds,
            &known_bonds,
        )?;
    }
    if params.flag_possible
        && (possible_atoms != original_possible_atoms || possible_bonds != original_possible_bonds)
    {
        possible_atoms = original_possible_atoms;
        for index in 0..atom_count {
            if !fixed_atoms[index] && known_atoms[index] {
                possible_atoms[index] = true;
                known_atoms[index] = false;
            }
            if possible_atoms[index] {
                atom_symbols[index].push('_');
                atom_symbols[index].push_str(&index.to_string());
            }
        }
        possible_bonds = original_possible_bonds;
        for index in 0..bond_count {
            if !fixed_bonds[index] && known_bonds[index] {
                possible_bonds[index] = true;
                known_bonds[index] = false;
            }
            if possible_bonds[index] {
                bond_symbols[index].push('_');
                bond_symbols[index].push_str(&index.to_string());
            }
        }
        flag_ring_stereo(
            &working,
            rings,
            &mut possible_ring_atoms,
            &mut possible_ring_bonds,
            &known_atoms,
            Some(&possible_atoms),
            &known_bonds,
            Some(&possible_bonds),
        );
        (stereo, ranks) = refinement_loop(
            &working,
            valence,
            rings,
            params.allow_nontetrahedral,
            &mut atom_symbols,
            &mut bond_symbols,
            &mut possible_atoms,
            &known_atoms,
            &mut fixed_atoms,
            &mut possible_bonds,
            &known_bonds,
            &mut fixed_bonds,
            &mut possible_ring_atoms,
            &mut possible_ring_bonds,
        )?;
    }
    let ring_relations = special_ring_relations(&working, valence, rings, &ranks)?;
    working.validate()?;
    Ok(PotentialStereoAssignment {
        stereo,
        atom_ranks: ranks,
        ring_relations,
        cleaned_topology: params.clean.then_some(working),
    })
}
