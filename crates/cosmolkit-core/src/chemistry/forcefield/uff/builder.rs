//! Source-backed RDKit UFF force-field builder helpers.

use crate::Molecule;
use crate::chemistry::forcefield::core::{ForceField, ForceFieldVec3};
use crate::chemistry::forcefield::torsion_query::{
    DEFAULT_TORSION_BOND_SMARTS, TorsionBondQueryError, match_torsion_bonds,
};
use crate::chemistry::rings::RingFindingError;
use crate::chemistry::valence::ValenceError;
use crate::notation::fragment::get_fragment_atom_mapping;
use crate::{AdjacencyList, AtomId, Hybridization};

use super::angle_bend::AngleBendContrib;
use super::atom_typer::{UffAtomTyperError, get_atom_types_for_uff};
use super::bond_stretch::BondStretchContrib;
use super::inversion::InversionContrib;
use super::nonbonded::VdwContrib;
use super::params::AtomicParams;
use super::torsion_angle::TorsionAngleContrib;
use super::utils::UffUtilsError;

pub(crate) const RELATION_1_2: u8 = 0;
pub(crate) const RELATION_1_3: u8 = 1;
pub(crate) const RELATION_1_4: u8 = 2;
pub(crate) const RELATION_1_X: u8 = 3;
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub(crate) enum UffBuilderError {
    #[error("RDKit UFF buildNeighborMatrix is unsupported for empty molecules")]
    EmptyMolecule,
    #[error("UFF constructForceField requires a 3D conformer when conf_id={conf_id}")]
    Missing3dConformer { conf_id: isize },
    #[error("UFF constructForceField conf_id must be >= -1, got {conf_id}")]
    InvalidConformerId { conf_id: isize },
    #[error("UFF addTrigonalBipyramidAngles requires central atom parameters for atom {atom_index}")]
    MissingTrigonalBipyramidCenterParams { atom_index: usize },
    #[error("UFF addBonds parameter length mismatch: atoms={atoms}, params={params}")]
    ParamsLengthMismatch { atoms: usize, params: usize },
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    TorsionBondQuery(#[from] TorsionBondQueryError),
    #[error(transparent)]
    AtomTyper(#[from] UffAtomTyperError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
    #[error(transparent)]
    UffUtils(#[from] UffUtilsError),
}

pub(crate) fn add_bonds(
    mol: &Molecule,
    params: &[Option<AtomicParams>],
    field: &mut ForceField,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addBonds (Builder.cpp:30-52)
    // RDKit✔️✔️: void addBonds(const ROMol &mol, const AtomicParamVect &params,
    // RDKit✔️✔️:               ForceFields::ForceField *field) {
    // RDKit✔️✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if mol.atoms().len() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: mol.atoms().len(),
            params: params.len(),
        });
    }
    // RDKit✔️✔️:   PRECONDITION(field, "bad forcefield");
    // Rust's mutable reference models RDKit's non-null forcefield precondition.

    // RDKit✔️✔️:   for (ROMol::ConstBondIterator bi = mol.beginBonds(); bi != mol.endBonds();
    // RDKit✔️✔️:        bi++) {
    for bond in mol.bonds() {
        // RDKit✔️✔️:     int idx1 = (*bi)->getBeginAtomIdx();
        let idx1 = bond.begin().index();
        // RDKit✔️✔️:     int idx2 = (*bi)->getEndAtomIdx();
        let idx2 = bond.end().index();

        // RDKit✔️✔️:     // FIX: recognize amide bonds here.

        // RDKit✔️✔️:     if (params[idx1] && params[idx2]) {
        if let (Some(end1_params), Some(end2_params)) = (&params[idx1], &params[idx2]) {
            // RDKit✔️✔️:       BondStretchContrib *contrib;
            // RDKit✔️✔️:       contrib = new BondStretchContrib(field, idx1, idx2,
            // RDKit✔️✔️:                                        (*bi)->getBondTypeAsDouble(),
            // RDKit✔️✔️:                                        params[idx1], params[idx2]);
            let contrib = BondStretchContrib::new(
                field,
                idx1,
                idx2,
                crate::chemistry::valence::bond_type_as_double(bond.order())?,
                end1_params,
                end2_params,
            )?;
            // RDKit✔️✔️:       field->contribs().push_back(ForceFields::ContribPtr(contrib));
            field.add_contrib(Box::new(contrib));
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::addBonds
    Ok(())
}

pub(crate) fn select_uff_conformer_index(mol: &Molecule, conf_id: isize) -> Result<usize, UffBuilderError> {
    let conformers = mol.conformers_3d();
    if conformers.is_empty() {
        return Err(UffBuilderError::Missing3dConformer { conf_id });
    }
    if conf_id == -1 {
        return Ok(0);
    }
    let requested = usize::try_from(conf_id).map_err(|_| UffBuilderError::InvalidConformerId { conf_id })?;
    conformers
        .iter()
        .position(|conformer| conformer.id() == requested)
        .ok_or(UffBuilderError::Missing3dConformer { conf_id })
}

fn direction_vector(from: ForceFieldVec3, to: ForceFieldVec3) -> ForceFieldVec3 {
    let mut direction = to - from;
    let length = direction.length();
    assert!(length > 0.0, "Cannot normalize a zero length vector");
    direction /= length;
    direction
}

fn bond_between_idx_simple(mol: &Molecule, a: usize, b: usize) -> Option<usize> {
    mol.bonds()
        .iter()
        .find(|bond| {
            (bond.begin().index() == a && bond.end().index() == b)
                || (bond.begin().index() == b && bond.end().index() == a)
        })
        .map(|bond| bond.id().index())
}

pub(crate) fn two_bit_cell_pos(n_atoms: usize, mut i: usize, mut j: usize) -> usize {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::twoBitCellPos (Builder.cpp:54-60)
    // RDKit✔️✔️: unsigned int twoBitCellPos(unsigned int nAtoms, int i, int j) {
    // RDKit✔️✔️:   if (j < i) {
    if j < i {
        // RDKit✔️✔️:     std::swap(i, j);
        std::mem::swap(&mut i, &mut j);
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return i * (nAtoms - 1) + i * (1 - i) / 2 + j;
    // Rust keeps the same triangular-index formula without signed temporary
    // underflow by using the algebraically equivalent subtraction form.
    i * (n_atoms - 1) - (i * i.saturating_sub(1)) / 2 + j
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::twoBitCellPos
}

pub(crate) fn set_two_bit_cell(res: &mut [u8], pos: usize, value: u8) {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::setTwoBitCell (Builder.cpp:62-68)
    // RDKit✔️✔️: void setTwoBitCell(boost::shared_array<std::uint8_t> &res, unsigned int pos,
    // RDKit✔️✔️:                    std::uint8_t value) {
    // RDKit✔️✔️:   unsigned int twoBitPos = pos / 4;
    let two_bit_pos = pos / 4;
    // RDKit✔️✔️:   unsigned int shift = 2 * (pos % 4);
    let shift = 2 * (pos % 4);
    // RDKit✔️✔️:   std::uint8_t twoBitMask = 3 << shift;
    let two_bit_mask = 3_u8 << shift;
    // RDKit✔️✔️:   res[twoBitPos] = ((res[twoBitPos] & (~twoBitMask)) | (value << shift));
    res[two_bit_pos] = (res[two_bit_pos] & !two_bit_mask) | (value << shift);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::setTwoBitCell
}

pub(crate) fn get_two_bit_cell(res: &[u8], pos: usize) -> u8 {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::getTwoBitCell (Builder.cpp:70-78)
    // RDKit✔️✔️: std::uint8_t getTwoBitCell(boost::shared_array<std::uint8_t> &res,
    // RDKit✔️✔️:                            unsigned int pos) {
    // RDKit✔️✔️:   unsigned int twoBitPos = pos / 4;
    let two_bit_pos = pos / 4;
    // RDKit✔️✔️:   unsigned int shift = 2 * (pos % 4);
    let shift = 2 * (pos % 4);
    // RDKit✔️✔️:   std::uint8_t twoBitMask = 3 << shift;
    let two_bit_mask = 3_u8 << shift;
    // RDKit✔️✔️:   return ((res[twoBitPos] & twoBitMask) >> shift);
    (res[two_bit_pos] & two_bit_mask) >> shift
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::getTwoBitCell
}

pub(crate) fn build_neighbor_matrix(mol: &Molecule) -> Result<Vec<u8>, UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::buildNeighborMatrix (Builder.cpp:91-130)
    // RDKit❗✔️: boost::shared_array<std::uint8_t> buildNeighborMatrix(const ROMol &mol) {
    // RDKit❗✔️:   const std::uint8_t RELATION_1_X_INIT = RELATION_1_X | (RELATION_1_X << 2) |
    // RDKit❗✔️:                                          (RELATION_1_X << 4) |
    // RDKit❗✔️:                                          (RELATION_1_X << 6);
    let relation_1_x_init = RELATION_1_X | (RELATION_1_X << 2) | (RELATION_1_X << 4) | (RELATION_1_X << 6);
    // RDKit❗✔️:   unsigned int nAtoms = mol.getNumAtoms();
    let n_atoms = mol.atoms().len();
    if n_atoms == 0 {
        // Rust models the empty-molecule edge explicitly instead of reproducing
        // RDKit's unsigned-arithmetic allocation path for nAtoms == 0.
        return Err(UffBuilderError::EmptyMolecule);
    }
    // RDKit❗✔️:   unsigned nTwoBitCells = (nAtoms * (nAtoms + 1) - 1) / 8 + 1;
    let n_two_bit_cells = (n_atoms * (n_atoms + 1) - 1) / 8 + 1;
    // RDKit❗✔️:   boost::shared_array<std::uint8_t> res(new std::uint8_t[nTwoBitCells]);
    // RDKit❗✔️:   std::memset(res.get(), RELATION_1_X_INIT, nTwoBitCells);
    let mut res = vec![relation_1_x_init; n_two_bit_cells];
    // RDKit❗✔️:   for (ROMol::ConstBondIterator bondi = mol.beginBonds();
    // RDKit❗✔️:        bondi != mol.endBonds(); ++bondi) {
    for (bond_i_index, bond_i) in mol.bonds().iter().enumerate() {
        // RDKit❗✔️:     setTwoBitCell(res,
        // RDKit❗✔️:                   twoBitCellPos(nAtoms, (*bondi)->getBeginAtomIdx(),
        // RDKit❗✔️:                                 (*bondi)->getEndAtomIdx()),
        // RDKit❗✔️:                   RELATION_1_2);
        set_two_bit_cell(
            &mut res,
            two_bit_cell_pos(n_atoms, bond_i.begin().index(), bond_i.end().index()),
            RELATION_1_2,
        );
        // RDKit❗✔️:     unsigned int bondiBeginAtomIdx = (*bondi)->getBeginAtomIdx();
        let bond_i_begin_atom_idx = bond_i.begin().index();
        // RDKit❗✔️:     unsigned int bondiEndAtomIdx = (*bondi)->getEndAtomIdx();
        let bond_i_end_atom_idx = bond_i.end().index();
        // RDKit❗✔️:     for (ROMol::ConstBondIterator bondj = bondi; ++bondj != mol.endBonds();) {
        for bond_j in mol.bonds().iter().skip(bond_i_index + 1) {
            // RDKit❗✔️:       int idx1 = -1;
            // RDKit❗✔️:       int idx3 = -1;
            let mut idx1 = None;
            let mut idx3 = None;
            // RDKit❗✔️:       unsigned int bondjBeginAtomIdx = (*bondj)->getBeginAtomIdx();
            let bond_j_begin_atom_idx = bond_j.begin().index();
            // RDKit❗✔️:       unsigned int bondjEndAtomIdx = (*bondj)->getEndAtomIdx();
            let bond_j_end_atom_idx = bond_j.end().index();
            // RDKit❗✔️:       if (bondiBeginAtomIdx == bondjBeginAtomIdx) {
            if bond_i_begin_atom_idx == bond_j_begin_atom_idx {
                // RDKit❗✔️:         idx1 = bondiEndAtomIdx;
                idx1 = Some(bond_i_end_atom_idx);
                // RDKit❗✔️:         idx3 = bondjEndAtomIdx;
                idx3 = Some(bond_j_end_atom_idx);
                // RDKit❗✔️:       } else if (bondiBeginAtomIdx == bondjEndAtomIdx) {
            } else if bond_i_begin_atom_idx == bond_j_end_atom_idx {
                // RDKit❗✔️:         idx1 = bondiEndAtomIdx;
                idx1 = Some(bond_i_end_atom_idx);
                // RDKit❗✔️:         idx3 = bondjBeginAtomIdx;
                idx3 = Some(bond_j_begin_atom_idx);
                // RDKit❗✔️:       } else if (bondiEndAtomIdx == bondjBeginAtomIdx) {
            } else if bond_i_end_atom_idx == bond_j_begin_atom_idx {
                // RDKit❗✔️:         idx1 = bondiBeginAtomIdx;
                idx1 = Some(bond_i_begin_atom_idx);
                // RDKit❗✔️:         idx3 = bondjEndAtomIdx;
                idx3 = Some(bond_j_end_atom_idx);
                // RDKit❗✔️:       } else if (bondiEndAtomIdx == bondjEndAtomIdx) {
            } else if bond_i_end_atom_idx == bond_j_end_atom_idx {
                // RDKit❗✔️:         idx1 = bondiBeginAtomIdx;
                idx1 = Some(bond_i_begin_atom_idx);
                // RDKit❗✔️:         idx3 = bondjBeginAtomIdx;
                idx3 = Some(bond_j_begin_atom_idx);
                // RDKit❗✔️:       }
            }
            // RDKit❗✔️:       if (idx1 > -1) {
            if let (Some(idx1), Some(idx3)) = (idx1, idx3) {
                // RDKit❗✔️:         setTwoBitCell(res, twoBitCellPos(nAtoms, idx1, idx3), RELATION_1_3);
                set_two_bit_cell(&mut res, two_bit_cell_pos(n_atoms, idx1, idx3), RELATION_1_3);
                // RDKit❗✔️:       }
            }
            // RDKit❗✔️:     }
        }
        // RDKit❗✔️:   }
    }
    // RDKit❗✔️:   return res;
    Ok(res)
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::buildNeighborMatrix
}

fn neighbor_bond_type_as_double(mol: &Molecule, center: usize, neighbor: usize) -> Result<f64, UffBuilderError> {
    let bond_id = mol
        .topology_block()
        .adjacency
        .neighbors_of(center)
        .iter()
        .find(|entry| entry.atom_index == neighbor)
        .expect("adjacency must contain bonded neighbor")
        .bond;
    crate::chemistry::valence::bond_type_as_double(mol.bonds()[bond_id.index()].order()).map_err(UffBuilderError::from)
}

pub(crate) fn add_angles(
    mol: &Molecule,
    params: &[Option<AtomicParams>],
    field: &mut ForceField,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addAngles (Builder.cpp:139-223)
    // RDKit✔️❌: void addAngles(const ROMol &mol, const AtomicParamVect &params,
    // RDKit✔️❌:                ForceFields::ForceField *field) {
    // RDKit✔️❌:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if mol.atoms().len() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: mol.atoms().len(),
            params: params.len(),
        });
    }
    // RDKit✔️❌:   PRECONDITION(field, "bad forcefield");
    // Rust's mutable reference models RDKit's non-null forcefield precondition.
    // RDKit✔️❌:   ROMol::ADJ_ITER nbr1Idx;
    // RDKit✔️❌:   ROMol::ADJ_ITER end1Nbrs;
    // RDKit✔️❌:   ROMol::ADJ_ITER nbr2Idx;
    // RDKit✔️❌:   ROMol::ADJ_ITER end2Nbrs;
    // RDKit✔️❌:   RingInfo *rings = mol.getRingInfo();
    // COSMolKit currently recomputes SSSR ring membership here because Molecule
    // does not yet expose RDKit-style cached RingInfo on the read path.
    let rings = crate::chemistry::rings::find_sssr(mol)?;

    // RDKit✔️❌:   unsigned int nAtoms = mol.getNumAtoms();
    let n_atoms = mol.atoms().len();
    // RDKit✔️❌:   for (unsigned int j = 0; j < nAtoms; j++) {
    for j in 0..n_atoms {
        // RDKit✔️❌:     if (!params[j]) {
        // RDKit✔️❌:       continue;
        // RDKit✔️❌:     }
        let Some(center_params) = &params[j] else {
            continue;
        };
        // RDKit✔️❌:     const Atom *atomJ = mol.getAtomWithIdx(j);
        let atom_j = &mol.atoms()[j];
        // RDKit✔️❌:     if (atomJ->getDegree() == 1) {
        // RDKit✔️❌:       continue;
        // RDKit✔️❌:     }
        let neighbors = mol.topology_block().adjacency.neighbors_of(j);
        if neighbors.len() == 1 {
            continue;
        }
        // RDKit✔️❌:     boost::tie(nbr1Idx, end1Nbrs) = mol.getAtomNeighbors(atomJ);
        // RDKit✔️❌:     for (; nbr1Idx != end1Nbrs; nbr1Idx++) {
        for (nbr1_pos, nbr1) in neighbors.iter().enumerate() {
            // RDKit✔️❌:       const Atom *atomI = mol[*nbr1Idx];
            // RDKit✔️❌:       unsigned int i = atomI->getIdx();
            let i = nbr1.atom_index;
            // RDKit✔️❌:       if (!params[i]) {
            // RDKit✔️❌:         continue;
            // RDKit✔️❌:       }
            let Some(atom_i_params) = &params[i] else {
                continue;
            };
            // RDKit✔️❌:       boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(atomJ);
            // RDKit✔️❌:       for (; nbr2Idx != end2Nbrs; nbr2Idx++) {
            for (nbr2_pos, nbr2) in neighbors.iter().enumerate() {
                // RDKit✔️❌:         if (nbr2Idx < (nbr1Idx + 1)) {
                // RDKit✔️❌:           continue;
                // RDKit✔️❌:         }
                if nbr2_pos < nbr1_pos + 1 {
                    continue;
                }
                // RDKit✔️❌:         const Atom *atomK = mol[*nbr2Idx];
                // RDKit✔️❌:         unsigned int k = atomK->getIdx();
                let k = nbr2.atom_index;
                // RDKit✔️❌:         if (!params[k]) {
                // RDKit✔️❌:           continue;
                // RDKit✔️❌:         }
                let Some(atom_k_params) = &params[k] else {
                    continue;
                };
                // RDKit✔️❌:         // skip special cases:
                // RDKit✔️❌:         if (!(atomJ->getHybridization() == Atom::SP3D &&
                // RDKit✔️❌:               atomJ->getDegree() == 5)) {
                if !(atom_j.hybridization() == Hybridization::Sp3d && neighbors.len() == 5) {
                    // RDKit✔️❌:           const Bond *b1 = mol.getBondBetweenAtoms(i, j);
                    // RDKit✔️❌:           const Bond *b2 = mol.getBondBetweenAtoms(k, j);
                    let bond_order_ij = neighbor_bond_type_as_double(mol, j, i)?;
                    let bond_order_kj = neighbor_bond_type_as_double(mol, j, k)?;
                    // RDKit✔️❌:           // FIX: recognize amide bonds here.
                    // RDKit✔️❌:           AngleBendContrib *contrib;
                    // RDKit✔️❌:           int order = 0;
                    let mut order = 0_u32;
                    // RDKit✔️❌:           switch (atomJ->getHybridization()) {
                    match atom_j.hybridization() {
                        // RDKit✔️❌:             case Atom::SP:
                        // RDKit✔️❌:               order = 1;
                        // RDKit✔️❌:               break;
                        Hybridization::Sp => order = 1,
                        // RDKit✔️❌:             case Atom::SP2:
                        // RDKit✔️❌:               order = 3;
                        Hybridization::Sp2 => {
                            order = 3;
                            // RDKit✔️❌:               // the following is a hack to get decent geometries
                            // RDKit✔️❌:               // with 3- and 4-membered rings incorporating sp2 atoms
                            // RDKit✔️❌:               // if the central atom is in a ring of size 3
                            if rings.is_atom_in_ring_of_size(AtomId::new(j), 3) {
                                // RDKit✔️❌:                 // if the central atom and one of the bonded atoms, but not the
                                // RDKit✔️❌:                 //  other one are inside a ring, then this angle is between a
                                // RDKit✔️❌:                 // ring substituent and a ring edge
                                if rings.is_atom_in_ring_of_size(AtomId::new(i), 3)
                                    != rings.is_atom_in_ring_of_size(AtomId::new(k), 3)
                                {
                                    // RDKit✔️❌:                   order = 30;
                                    order = 30;
                                // RDKit✔️❌:                 // if all atoms are inside the ring, then this is one of ring
                                // RDKit✔️❌:                 // angles
                                } else if rings.is_atom_in_ring_of_size(AtomId::new(i), 3)
                                    && rings.is_atom_in_ring_of_size(AtomId::new(k), 3)
                                {
                                    // RDKit✔️❌:                   order = 35;
                                    order = 35;
                                }
                            // RDKit✔️❌:               // if the central atom is in a ring of size 4
                            } else if rings.is_atom_in_ring_of_size(AtomId::new(j), 4) {
                                // RDKit✔️❌:                 // if the central atom and one of the bonded atoms, but not the
                                // RDKit✔️❌:                 //  other one are inside a ring, then this angle is between a
                                // RDKit✔️❌:                 // ring substituent and a ring edge
                                if rings.is_atom_in_ring_of_size(AtomId::new(i), 4)
                                    != rings.is_atom_in_ring_of_size(AtomId::new(k), 4)
                                {
                                    // RDKit✔️❌:                   order = 40;
                                    order = 40;
                                // RDKit✔️❌:                 // if all atoms are inside the ring, then this is one of ring
                                // RDKit✔️❌:                 // angles
                                } else if rings.is_atom_in_ring_of_size(AtomId::new(i), 4)
                                    && rings.is_atom_in_ring_of_size(AtomId::new(k), 4)
                                {
                                    // RDKit✔️❌:                   order = 45;
                                    order = 45;
                                }
                            }
                            // RDKit✔️❌:               // end of the hack
                            // RDKit✔️❌:               break;
                        }
                        // RDKit✔️❌:             case Atom::SP3D2:
                        // RDKit✔️❌:               order = 4;
                        // RDKit✔️❌:               break;
                        Hybridization::Sp3d2 => order = 4,
                        // RDKit✔️❌:             default:
                        // RDKit✔️❌:               order = 0;
                        // RDKit✔️❌:               break;
                        _ => order = 0,
                    }

                    // RDKit✔️❌:           contrib =
                    // RDKit✔️❌:               new AngleBendContrib(field, i, j, k, b1->getBondTypeAsDouble(),
                    // RDKit✔️❌:                                    b2->getBondTypeAsDouble(), params[i],
                    // RDKit✔️❌:                                    params[j], params[k], order);
                    let contrib = AngleBendContrib::new(
                        field,
                        i,
                        j,
                        k,
                        bond_order_ij,
                        bond_order_kj,
                        atom_i_params,
                        center_params,
                        atom_k_params,
                        order,
                    )?;
                    // RDKit✔️❌:           field->contribs().push_back(ForceFields::ContribPtr(contrib));
                    field.add_contrib(Box::new(contrib));
                    // RDKit✔️❌:         }
                }
                // RDKit✔️❌:       }
            }
            // RDKit✔️❌:     }
        }
        // RDKit✔️❌:   }
    }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION UFF::Tools::addAngles
    Ok(())
}

pub(crate) fn add_torsions(
    mol: &Molecule,
    params: &[Option<AtomicParams>],
    field: &mut ForceField,
    torsion_bond_smarts: &str,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addTorsions (Builder.cpp:499-585)
    // RDKit✔️❌: void addTorsions(const ROMol &mol, const AtomicParamVect &params,
    // RDKit✔️❌:                  ForceFields::ForceField *field,
    // RDKit✔️❌:                  const std::string &torsionBondSmarts) {
    // RDKit✔️❌:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if mol.atoms().len() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: mol.atoms().len(),
            params: params.len(),
        });
    }
    // RDKit✔️❌:   PRECONDITION(field, "bad forcefield");
    // Rust's mutable reference models RDKit's non-null forcefield precondition.

    // RDKit✔️❌:   // find all of the torsion bonds:
    // RDKit✔️❌:   std::vector<MatchVectType> matchVect;
    // RDKit✔️❌:   const ROMol *defaultQuery = DefaultTorsionBondSmarts::query();
    // RDKit✔️❌:   const ROMol *query = (torsionBondSmarts == DefaultTorsionBondSmarts::string())
    // RDKit✔️❌:                            ? defaultQuery
    // RDKit✔️❌:                            : SmartsToMol(torsionBondSmarts);
    // RDKit✔️❌:   TEST_ASSERT(query);
    let torsion_bond_matches = match_torsion_bonds(mol, torsion_bond_smarts)?;
    // RDKit✔️❌:   unsigned int nHits = SubstructMatch(mol, *query, matchVect);
    // RDKit✔️❌:   if (query != defaultQuery) {
    // RDKit✔️❌:     delete query;
    // RDKit✔️❌:   }

    // RDKit✔️❌:   for (unsigned int i = 0; i < nHits; i++) {
    for matched in &torsion_bond_matches {
        let idx1 = matched.begin_atom_index;
        let idx2 = matched.end_atom_index;
        let bond_idx = matched.bond_index;
        // RDKit✔️❌:     MatchVectType match = matchVect[i];
        // RDKit✔️❌:     TEST_ASSERT(match.size() == 2);
        // RDKit✔️❌:     int idx1 = match[0].second;
        // RDKit✔️❌:     int idx2 = match[1].second;
        // RDKit✔️❌:     if (!params[idx1] || !params[idx2]) {
        // RDKit✔️❌:       continue;
        // RDKit✔️❌:     }
        let (Some(atom1_params), Some(atom2_params)) = (&params[idx1], &params[idx2]) else {
            continue;
        };
        // RDKit✔️❌:     const Bond *bond = mol.getBondBetweenAtoms(idx1, idx2);
        let bond = &mol.bonds()[bond_idx];
        // RDKit✔️❌:     std::vector<TorsionAngleContrib *> contribsHere;
        let mut contribs_here = Vec::new();
        // RDKit✔️❌:     TEST_ASSERT(bond);
        // The query layer preserves source endpoint orientation and resolves
        // the corresponding topology bond before contribution construction.
        // RDKit✔️❌:     const Atom *atom1 = mol.getAtomWithIdx(idx1);
        let atom1 = &mol.atoms()[idx1];
        // RDKit✔️❌:     const Atom *atom2 = mol.getAtomWithIdx(idx2);
        let atom2 = &mol.atoms()[idx2];

        // RDKit✔️❌:     if ((atom1->getHybridization() == Atom::SP2 ||
        // RDKit✔️❌:          atom1->getHybridization() == Atom::SP3) &&
        // RDKit✔️❌:         (atom2->getHybridization() == Atom::SP2 ||
        // RDKit✔️❌:          atom2->getHybridization() == Atom::SP3)) {
        if (atom1.hybridization() == Hybridization::Sp2 || atom1.hybridization() == Hybridization::Sp3)
            && (atom2.hybridization() == Hybridization::Sp2 || atom2.hybridization() == Hybridization::Sp3)
        {
            // RDKit✔️❌:       ROMol::OEDGE_ITER beg1, end1;
            // RDKit✔️❌:       boost::tie(beg1, end1) = mol.getAtomBonds(atom1);
            // RDKit✔️❌:       while (beg1 != end1) {
            for first_neighbor in mol.topology_block().adjacency.neighbors_of(idx1) {
                let t_bond1 = &mol.bonds()[first_neighbor.bond.index()];
                // RDKit✔️❌:         const Bond *tBond1 = mol[*beg1];
                // RDKit✔️❌:         if (tBond1 != bond) {
                if t_bond1.id() != bond.id() {
                    // RDKit✔️❌:           int bIdx = tBond1->getOtherAtomIdx(idx1);
                    let b_idx = first_neighbor.atom_index;
                    // RDKit✔️❌:           ROMol::OEDGE_ITER beg2, end2;
                    // RDKit✔️❌:           boost::tie(beg2, end2) = mol.getAtomBonds(atom2);
                    // RDKit✔️❌:           while (beg2 != end2) {
                    for second_neighbor in mol.topology_block().adjacency.neighbors_of(idx2) {
                        let t_bond2 = &mol.bonds()[second_neighbor.bond.index()];
                        // RDKit✔️❌:             const Bond *tBond2 = mol[*beg2];
                        // RDKit✔️❌:             if (tBond2 != bond && tBond2 != tBond1) {
                        if t_bond2.id() != bond.id() && t_bond2.id() != t_bond1.id() {
                            // RDKit✔️❌:               int eIdx = tBond2->getOtherAtomIdx(idx2);
                            let e_idx = second_neighbor.atom_index;
                            // RDKit✔️❌:               // make sure this isn't a three-membered ring:
                            // RDKit✔️❌:               if (eIdx != bIdx) {
                            if e_idx != b_idx {
                                // RDKit✔️❌:                 // we now have a torsion involving atoms (bonds):
                                // RDKit✔️❌:                 //  bIdx - (tBond1) - idx1 - (bond) - idx2 - (tBond2) - eIdx
                                // RDKit✔️❌:                 TorsionAngleContrib *contrib;

                                // RDKit✔️❌:                 // if either of the end atoms is SP2 hybridized, set a flag
                                // RDKit✔️❌:                 // here.
                                // RDKit✔️❌:                 bool hasSP2 = false;
                                let has_sp2 = mol.atoms()[b_idx].hybridization() == Hybridization::Sp2
                                    || mol.atoms()[e_idx].hybridization() == Hybridization::Sp2;
                                // RDKit✔️❌:                 if (mol.getAtomWithIdx(bIdx)->getHybridization() == Atom::SP2 ||
                                // RDKit✔️❌:                     mol.getAtomWithIdx(eIdx)->getHybridization() == Atom::SP2) {
                                // RDKit✔️❌:                   hasSP2 = true;
                                // RDKit✔️❌:                 }
                                // RDKit✔️❌:                 // std::cout << "Torsion: " << bIdx << "-" << idx1 << "-" <<
                                // RDKit✔️❌:                 // idx2 << "-" << eIdx << std::endl;
                                // RDKit✔️❌:                 // if(okToIncludeTorsion(mol,bond,bIdx,idx1,idx2,eIdx)){
                                // RDKit✔️❌:                 // std::cout << "  INCLUDED" << std::endl;
                                let contrib = TorsionAngleContrib::new(
                                    field,
                                    b_idx,
                                    idx1,
                                    idx2,
                                    e_idx,
                                    crate::chemistry::valence::bond_type_as_double(bond.order())?,
                                    i32::from(atom1.atomic_number()),
                                    i32::from(atom2.atomic_number()),
                                    atom1.hybridization(),
                                    atom2.hybridization(),
                                    atom1_params,
                                    atom2_params,
                                    has_sp2,
                                );
                                // RDKit✔️❌:                 field->contribs().push_back(ForceFields::ContribPtr(contrib));
                                contribs_here.push(contrib);
                                // RDKit✔️❌:                 contribsHere.push_back(contrib);
                                // RDKit✔️❌:                 //}
                            }
                        }
                    }
                }
            }
        }
        // RDKit✔️❌:     }
        // RDKit✔️❌:     // now divide the force constant for each contribution to the torsion energy
        // RDKit✔️❌:     // about this bond by the number of contribs about this bond:
        // RDKit✔️❌:     for (auto chI = contribsHere.begin(); chI != contribsHere.end(); ++chI) {
        let contrib_count = contribs_here.len();
        for contrib in &mut contribs_here {
            // RDKit✔️❌:       (*chI)->scaleForceConstant(contribsHere.size());
            contrib.scale_force_constant(contrib_count);
            // RDKit✔️❌:     }
        }
        for contrib in contribs_here {
            field.add_contrib(Box::new(contrib));
        }
    }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION UFF::Tools::addTorsions
    Ok(())
}

pub(crate) fn add_inversions(
    mol: &Molecule,
    params: &[Option<AtomicParams>],
    field: &mut ForceField,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addInversions (Builder.cpp:593-657)
    // RDKit✔️✔️: void addInversions(const ROMol &mol, const AtomicParamVect &params,
    // RDKit✔️✔️:                    ForceFields::ForceField *field) {
    // RDKit✔️✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if mol.atoms().len() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: mol.atoms().len(),
            params: params.len(),
        });
    }
    // RDKit✔️✔️:   PRECONDITION(field, "bad forcefield");
    // Rust's mutable reference models RDKit's non-null forcefield precondition.

    // RDKit✔️✔️:   unsigned int idx[4];
    // RDKit✔️✔️:   unsigned int n[4];
    let mut idx = [0_usize; 4];
    let mut n = [0_usize; 4];
    // RDKit✔️✔️:   const Atom *atom[4];
    // Rust reads atoms directly from the molecule tables when needed.

    // RDKit✔️✔️:   for (idx[1] = 0; idx[1] < mol.getNumAtoms(); ++idx[1]) {
    for center in 0..mol.num_atoms() {
        idx[1] = center;
        // RDKit✔️✔️:     atom[1] = mol.getAtomWithIdx(idx[1]);
        let center_atom = &mol.atoms()[center];
        // RDKit✔️✔️:     int at2AtomicNum = atom[1]->getAtomicNum();
        let at2_atomic_num = i32::from(center_atom.atomic_number());
        let neighbors = mol.topology_block().adjacency.neighbors_of(center);
        // RDKit✔️✔️:     // if the central atom is not carbon, nitrogen, oxygen,
        // RDKit✔️✔️:     // phosphorous, arsenic, antimonium or bismuth, skip it
        // RDKit✔️✔️:     if (((at2AtomicNum != 6) && (at2AtomicNum != 7) && (at2AtomicNum != 8) &&
        // RDKit✔️✔️:          (at2AtomicNum != 15) && (at2AtomicNum != 33) && (at2AtomicNum != 51) &&
        // RDKit✔️✔️:          (at2AtomicNum != 83)) ||
        // RDKit✔️✔️:         (atom[1]->getDegree() != 3)) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if !matches!(at2_atomic_num, 6 | 7 | 8 | 15 | 33 | 51 | 83) || neighbors.len() != 3 {
            continue;
        }
        // RDKit✔️✔️:     // if the central atom is carbon, nitrogen or oxygen
        // RDKit✔️✔️:     // but hybridization is not sp2, skip it
        // RDKit✔️✔️:     if (((at2AtomicNum == 6) || (at2AtomicNum == 7) || (at2AtomicNum == 8)) &&
        // RDKit✔️✔️:         (atom[1]->getHybridization() != Atom::SP2)) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if matches!(at2_atomic_num, 6 | 7 | 8) && center_atom.hybridization() != Hybridization::Sp2 {
            continue;
        }
        // RDKit✔️✔️:     boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom[1]);
        // RDKit✔️✔️:     unsigned int i = 0;
        let mut i = 0_usize;
        // RDKit✔️✔️:     bool isBoundToSP2O = false;
        let mut is_bound_to_sp2_o = false;
        // RDKit✔️✔️:     for (; nbrIdx != endNbrs; ++nbrIdx) {
        for neighbor in neighbors {
            let neighbor_atom = &mol.atoms()[neighbor.atom_index];
            // RDKit✔️✔️:       atom[i] = mol[*nbrIdx];
            // RDKit✔️✔️:       idx[i] = atom[i]->getIdx();
            idx[i] = neighbor_atom.id().index();
            // RDKit✔️✔️:       // if the central atom is sp2 carbon and is
            // RDKit✔️✔️:       // bound to sp2 oxygen, set a flag
            // RDKit✔️✔️:       if (!isBoundToSP2O) {
            // RDKit✔️✔️:         isBoundToSP2O =
            // RDKit✔️✔️:             ((at2AtomicNum == 6) && (atom[i]->getAtomicNum() == 8) &&
            // RDKit✔️✔️:              (atom[i]->getHybridization() == Atom::SP2));
            // RDKit✔️✔️:       }
            if !is_bound_to_sp2_o {
                is_bound_to_sp2_o = at2_atomic_num == 6
                    && neighbor_atom.atomic_number() == 8
                    && neighbor_atom.hybridization() == Hybridization::Sp2;
            }
            // RDKit✔️✔️:       if (!i) {
            // RDKit✔️✔️:         ++i;
            // RDKit✔️✔️:       }
            if i == 0 {
                i += 1;
            }
            // RDKit✔️✔️:       ++i;
            i += 1;
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     for (unsigned int i = 0; i < 3; ++i) {
        for i in 0..3 {
            // RDKit✔️✔️:       n[1] = 1;
            n[1] = 1;
            // RDKit✔️✔️:       switch (i) {
            match i {
                // RDKit✔️✔️:         case 0:
                // RDKit✔️✔️:           n[0] = 0;
                // RDKit✔️✔️:           n[2] = 2;
                // RDKit✔️✔️:           n[3] = 3;
                // RDKit✔️✔️:           break;
                0 => {
                    n[0] = 0;
                    n[2] = 2;
                    n[3] = 3;
                }
                // RDKit✔️✔️:         case 1:
                // RDKit✔️✔️:           n[0] = 0;
                // RDKit✔️✔️:           n[2] = 3;
                // RDKit✔️✔️:           n[3] = 2;
                // RDKit✔️✔️:           break;
                1 => {
                    n[0] = 0;
                    n[2] = 3;
                    n[3] = 2;
                }
                // RDKit✔️✔️:         case 2:
                // RDKit✔️✔️:           n[0] = 2;
                // RDKit✔️✔️:           n[2] = 3;
                // RDKit✔️✔️:           n[3] = 0;
                // RDKit✔️✔️:           break;
                2 => {
                    n[0] = 2;
                    n[2] = 3;
                    n[3] = 0;
                }
                _ => unreachable!("loop bounds guarantee 0..3"),
            }
            // RDKit✔️✔️:       InversionContrib *contrib;
            // RDKit✔️✔️:       contrib = new InversionContrib(field, idx[n[0]], idx[n[1]], idx[n[2]],
            // RDKit✔️✔️:                                      idx[n[3]], at2AtomicNum, isBoundToSP2O);
            let contrib = InversionContrib::new(
                field,
                idx[n[0]],
                idx[n[1]],
                idx[n[2]],
                idx[n[3]],
                at2_atomic_num,
                is_bound_to_sp2_o,
                1.0,
            );
            // RDKit✔️✔️:       field->contribs().push_back(ForceFields::ContribPtr(contrib));
            field.add_contrib(Box::new(contrib));
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::addInversions
    Ok(())
}

pub(crate) fn add_nonbonded(
    mol: &Molecule,
    params: &[Option<AtomicParams>],
    field: &mut ForceField,
    neighbor_matrix: &[u8],
    vdw_thresh: f64,
    ignore_interfrag_interactions: bool,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addNonbonded (Builder.cpp:432-466)
    // RDKit✔️✔️: void addNonbonded(const ROMol &mol, int confId, const AtomicParamVect &params,
    // RDKit✔️✔️:                   ForceFields::ForceField *field,
    // RDKit✔️✔️:                   boost::shared_array<std::uint8_t> neighborMatrix,
    // RDKit✔️✔️:                   double vdwThresh, bool ignoreInterfragInteractions) {
    // RDKit✔️✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if mol.num_atoms() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: mol.num_atoms(),
            params: params.len(),
        });
    }
    // RDKit✔️✔️:   PRECONDITION(field, "bad forcefield");
    // Rust's mutable reference models RDKit's non-null forcefield precondition.

    // RDKit✔️✔️:   INT_VECT fragMapping;
    // RDKit✔️✔️:   if (ignoreInterfragInteractions) {
    // RDKit✔️✔️:     std::vector<ROMOL_SPTR> molFrags =
    // RDKit✔️✔️:         MolOps::getMolFrags(mol, true, &fragMapping);
    // RDKit✔️✔️:   }
    let frag_mapping = if ignore_interfrag_interactions {
        Some(get_fragment_atom_mapping(mol))
    } else {
        None
    };

    // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
    let n_atoms = mol.num_atoms();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; i++) {
    for i in 0..n_atoms {
        // RDKit✔️✔️:     if (!params[i]) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        let Some(at1_params) = &params[i] else {
            continue;
        };
        // RDKit✔️✔️:     for (unsigned int j = i + 1; j < nAtoms; j++) {
        for j in (i + 1)..n_atoms {
            // RDKit✔️✔️:       if (!params[j] ||
            // RDKit✔️✔️:           (ignoreInterfragInteractions && fragMapping[i] != fragMapping[j])) {
            // RDKit✔️✔️:         continue;
            // RDKit✔️✔️:       }
            let Some(at2_params) = &params[j] else {
                continue;
            };
            if ignore_interfrag_interactions && frag_mapping.as_ref().is_some_and(|mapping| mapping[i] != mapping[j]) {
                continue;
            }
            // RDKit✔️✔️:       if (getTwoBitCell(neighborMatrix, twoBitCellPos(nAtoms, i, j)) >=
            // RDKit✔️✔️:           RELATION_1_4) {
            if get_two_bit_cell(neighbor_matrix, two_bit_cell_pos(n_atoms, i, j)) >= RELATION_1_4 {
                // RDKit✔️✔️:         double dist = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
                // COSMolKit's builder path reads coordinates from the target
                // force field positions rather than from a Molecule conformer
                // handle.
                let pos_i = field.positions()[i];
                let pos_j = field.positions()[j];
                let dist = (pos_i - pos_j).length();
                // RDKit✔️✔️:         if (dist < vdwThresh *
                // RDKit✔️✔️:                        UFF::Utils::calcNonbondedMinimum(params[i], params[j])) {
                let cutoff = vdw_thresh * super::utils::calc_nonbonded_minimum(at1_params, at2_params);
                if dist < cutoff {
                    // RDKit✔️✔️:           vdWContrib *contrib;
                    // RDKit✔️✔️:           contrib = new vdWContrib(field, i, j, params[i], params[j]);
                    let contrib = VdwContrib::new(field, i, j, at1_params, at2_params, 10.0);
                    // RDKit✔️✔️:           field->contribs().push_back(ForceFields::ContribPtr(contrib));
                    field.add_contrib(Box::new(contrib));
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::addNonbonded
    Ok(())
}

fn add_trigonal_bipyramid_angles(
    atom_idx: usize,
    mol: &Molecule,
    adjacency: &AdjacencyList,
    params: &[Option<AtomicParams>],
    field: &mut ForceField,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles (Builder.cpp:226-410)
    // RDKit✔️❌: void addTrigonalBipyramidAngles(const Atom *atom, const ROMol &mol, int confId,
    // RDKit✔️❌:                                 const AtomicParamVect &params,
    // RDKit✔️❌:                                 ForceFields::ForceField *field) {
    let atom = &mol.atoms()[atom_idx];
    // RDKit✔️❌:   PRECONDITION(atom, "bad atom");
    // Rust's indexed atom lookup models RDKit's non-null atom precondition.
    // RDKit✔️❌:   PRECONDITION(atom->getHybridization() == Atom::SP3D, "bad hybridization");
    assert!(atom.hybridization() == Hybridization::Sp3d, "bad hybridization");
    // RDKit✔️❌:   PRECONDITION(atom->getDegree() == 5, "bad degree");
    assert!(adjacency.neighbors_of(atom_idx).len() == 5, "bad degree");
    // RDKit✔️❌:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if mol.num_atoms() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: mol.num_atoms(),
            params: params.len(),
        });
    }
    // RDKit✔️❌:   PRECONDITION(field, "bad forcefield");
    // Rust's mutable reference models RDKit's non-null forcefield precondition.
    let Some(center_params) = params[atom_idx].as_ref() else {
        return Err(UffBuilderError::MissingTrigonalBipyramidCenterParams { atom_index: atom_idx });
    };

    // RDKit✔️❌:   const Bond *ax1 = nullptr, *ax2 = nullptr;
    // RDKit✔️❌:   const Bond *eq1 = nullptr, *eq2 = nullptr, *eq3 = nullptr;
    let mut axial_pair = None;
    let mut most_neg = 100.0_f64;
    let neighbors = adjacency.neighbors_of(atom_idx);
    let atom_pos = field.positions()[atom_idx];

    // RDKit✔️❌:   const Conformer &conf = mol.getConformer(confId);
    // RDKit✔️❌:   //------------------------------------------------------------
    // RDKit✔️❌:   // identify the axial and equatorial bonds:
    // RDKit✔️❌:   double mostNeg = 100.0;
    // RDKit✔️❌:   ROMol::OEDGE_ITER beg1, end1;
    // RDKit✔️❌:   boost::tie(beg1, end1) = mol.getAtomBonds(atom);
    // RDKit✔️❌:   unsigned int aid = atom->getIdx();
    // RDKit✔️❌:   while (beg1 != end1) {
    // RDKit✔️❌:     const Bond *bond1 = mol[*beg1];
    // RDKit✔️❌:     unsigned int oaid = bond1->getOtherAtomIdx(aid);
    // RDKit✔️❌:     RDGeom::Point3D v1 =
    // RDKit✔️❌:         conf.getAtomPos(aid).directionVector(conf.getAtomPos(oaid));
    for (i, neighbor1) in neighbors.iter().enumerate() {
        let v1 = direction_vector(atom_pos, field.positions()[neighbor1.atom_index]);
        // RDKit✔️❌:     ROMol::OEDGE_ITER beg2, end2;
        // RDKit✔️❌:     boost::tie(beg2, end2) = mol.getAtomBonds(atom);
        // RDKit✔️❌:     while (beg2 != end2) {
        // RDKit✔️❌:       const Bond *bond2 = mol[*beg2];
        // RDKit✔️❌:       if (bond2->getIdx() > bond1->getIdx()) {
        // RDKit✔️❌:         unsigned int oaid2 = bond2->getOtherAtomIdx(aid);
        // RDKit✔️❌:         RDGeom::Point3D v2 =
        // RDKit✔️❌:             conf.getAtomPos(aid).directionVector(conf.getAtomPos(oaid2));
        // RDKit✔️❌:         double dot = v1.dotProduct(v2);
        // RDKit✔️❌:         if (dot < mostNeg) {
        // RDKit✔️❌:           mostNeg = dot;
        // RDKit✔️❌:           ax1 = bond1;
        // RDKit✔️❌:           ax2 = bond2;
        // RDKit✔️❌:         }
        // RDKit✔️❌:       }
        // RDKit✔️❌:       ++beg2;
        // RDKit✔️❌:     }
        // RDKit✔️❌:     ++beg1;
        // RDKit✔️❌:   }
        for neighbor2 in neighbors.iter().skip(i + 1) {
            let v2 = direction_vector(atom_pos, field.positions()[neighbor2.atom_index]);
            let dot = v1.dot_product(v2);
            if dot < most_neg {
                most_neg = dot;
                axial_pair = Some((neighbor1.bond.index(), neighbor2.bond.index()));
            }
        }
    }
    let (ax1_idx, ax2_idx) = axial_pair.expect("axial bond not found");

    // RDKit✔️❌:   CHECK_INVARIANT(ax1, "axial bond not found");
    // RDKit✔️❌:   CHECK_INVARIANT(ax2, "axial bond not found");
    let mut equatorial = Vec::with_capacity(3);
    // RDKit✔️❌:   boost::tie(beg1, end1) = mol.getAtomBonds(atom);
    // RDKit✔️❌:   while (beg1 != end1) {
    // RDKit✔️❌:     const Bond *bond = mol[*beg1];
    // RDKit✔️❌:     ++beg1;
    // RDKit✔️❌:     if (bond == ax1 || bond == ax2) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (!eq1) {
    // RDKit✔️❌:       eq1 = bond;
    // RDKit✔️❌:     } else if (!eq2) {
    // RDKit✔️❌:       eq2 = bond;
    // RDKit✔️❌:     } else if (!eq3) {
    // RDKit✔️❌:       eq3 = bond;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    for neighbor in neighbors {
        let bond_idx = neighbor.bond.index();
        if bond_idx != ax1_idx && bond_idx != ax2_idx {
            equatorial.push(bond_idx);
        }
    }

    // RDKit✔️❌:   CHECK_INVARIANT(eq1, "equatorial bond not found");
    // RDKit✔️❌:   CHECK_INVARIANT(eq2, "equatorial bond not found");
    // RDKit✔️❌:   CHECK_INVARIANT(eq3, "equatorial bond not found");
    assert!(equatorial.len() == 3, "equatorial bond not found");

    let add_angle = |field: &mut ForceField,
                     first_bond_idx: usize,
                     second_bond_idx: usize,
                     order: u32|
     -> Result<(), UffBuilderError> {
        let first_bond = &mol.bonds()[first_bond_idx];
        let second_bond = &mol.bonds()[second_bond_idx];
        let i = if first_bond.begin().index() == atom_idx {
            first_bond.end().index()
        } else {
            first_bond.begin().index()
        };
        let j = if second_bond.begin().index() == atom_idx {
            second_bond.end().index()
        } else {
            second_bond.begin().index()
        };
        if let (Some(first_params), Some(second_params)) = (&params[i], &params[j]) {
            let contrib = AngleBendContrib::new(
                field,
                i,
                atom_idx,
                j,
                crate::chemistry::valence::bond_type_as_double(first_bond.order())?,
                crate::chemistry::valence::bond_type_as_double(second_bond.order())?,
                first_params,
                center_params,
                second_params,
                order,
            )?;
            field.add_contrib(Box::new(contrib));
        }
        Ok(())
    };

    // RDKit✔️❌:   //------------------------------------------------------------
    // RDKit✔️❌:   // alright, add the angles:
    // RDKit✔️❌:   AngleBendContrib *contrib;
    // RDKit✔️❌:   int atomIdx = atom->getIdx();
    // RDKit✔️❌:   int i, j;
    // RDKit✔️❌:
    // RDKit✔️❌:   // Axial-Axial
    add_angle(field, ax1_idx, ax2_idx, 2)?;
    // RDKit✔️❌:   // Equatorial-Equatorial
    add_angle(field, equatorial[0], equatorial[1], 3)?;
    add_angle(field, equatorial[0], equatorial[2], 3)?;
    add_angle(field, equatorial[1], equatorial[2], 3)?;
    // RDKit✔️❌:   // Axial-Equatorial
    add_angle(field, ax1_idx, equatorial[0], 0)?;
    add_angle(field, ax1_idx, equatorial[1], 0)?;
    add_angle(field, ax1_idx, equatorial[2], 0)?;
    add_angle(field, ax2_idx, equatorial[0], 0)?;
    add_angle(field, ax2_idx, equatorial[1], 0)?;
    add_angle(field, ax2_idx, equatorial[2], 0)?;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles
    Ok(())
}

fn add_angle_special_cases(
    mol: &Molecule,
    params: &[Option<AtomicParams>],
    field: &mut ForceField,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addAngleSpecialCases (Builder.cpp:412-430)
    // RDKit✔️❌: void addAngleSpecialCases(const ROMol &mol, int confId,
    // RDKit✔️❌:                           const AtomicParamVect &params,
    // RDKit✔️❌:                           ForceFields::ForceField *field) {
    // RDKit✔️❌:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if mol.num_atoms() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: mol.num_atoms(),
            params: params.len(),
        });
    }
    // RDKit✔️❌:   PRECONDITION(field, "bad forcefield");
    // Rust's mutable reference models RDKit's non-null forcefield precondition.
    let adjacency = AdjacencyList::from_topology(mol.num_atoms(), mol.bonds());
    // RDKit✔️❌:   unsigned int nAtoms = mol.getNumAtoms();
    // RDKit✔️❌:   for (unsigned int i = 0; i < nAtoms; i++) {
    for atom_idx in 0..mol.num_atoms() {
        let atom = &mol.atoms()[atom_idx];
        // RDKit✔️❌:     const Atom *atom = mol.getAtomWithIdx(i);
        // RDKit✔️❌:     // trigonal bipyramidal:
        // RDKit✔️❌:     if ((atom->getHybridization() == Atom::SP3D && atom->getDegree() == 5)) {
        if atom.hybridization() == Hybridization::Sp3d && adjacency.neighbors_of(atom_idx).len() == 5 {
            // RDKit✔️❌:       addTrigonalBipyramidAngles(atom, mol, confId, params, field);
            add_trigonal_bipyramid_angles(atom_idx, mol, &adjacency, params, field)?;
            // RDKit✔️❌:     }
        }
    }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION UFF::Tools::addAngleSpecialCases
    Ok(())
}

pub(crate) fn construct_force_field_with_params(
    mol: &Molecule,
    params: &[Option<AtomicParams>],
    vdw_thresh: f64,
    conf_id: isize,
    ignore_interfrag_interactions: bool,
) -> Result<ForceField, UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::constructForceField(ROMol&, const AtomicParamVect&, double, int, bool) (Builder.cpp:674-709)
    // RDKit✔️✔️: ForceFields::ForceField *constructForceField(ROMol &mol,
    // RDKit✔️✔️:                                              const AtomicParamVect &params,
    // RDKit✔️✔️:                                              double vdwThresh, int confId,
    // RDKit✔️✔️:                                              bool ignoreInterfragInteractions) {
    // RDKit✔️✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if mol.num_atoms() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: mol.num_atoms(),
            params: params.len(),
        });
    }

    // RDKit✔️✔️:   if (MolOps::needsHs(mol)) {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "Molecule does not have explicit Hs. Consider calling AddHs()"
    // RDKit✔️✔️:         << std::endl;
    // RDKit✔️✔️:   }
    // COSMolKit currently omits RDKit's warning-only `needsHs()` logging side effect.

    // RDKit✔️✔️:   std::unique_ptr<ForceFields::ForceField> res(new ForceFields::ForceField());
    let mut res = ForceField::new(3);

    // RDKit✔️✔️:   // add the atomic positions:
    // RDKit✔️✔️:   Conformer &conf = mol.getConformer(confId);
    let conf_index = select_uff_conformer_index(mol, conf_id)?;
    let conf = &mol.conformers_3d()[conf_index];
    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    for coords in conf.coordinates() {
        // RDKit✔️✔️:     res->positions().push_back(&conf.getAtomPos(i));
        res.positions_mut()
            .push(ForceFieldVec3::new(coords[0], coords[1], coords[2]));
    }

    // RDKit✔️✔️:   Tools::addBonds(mol, params, res.get());
    add_bonds(mol, params, &mut res)?;
    // RDKit✔️✔️:   Tools::addAngles(mol, params, res.get());
    add_angles(mol, params, &mut res)?;
    // RDKit✔️✔️:   Tools::addAngleSpecialCases(mol, confId, params, res.get());
    add_angle_special_cases(mol, params, &mut res)?;
    // RDKit✔️✔️:   boost::shared_array<std::uint8_t> neighborMat =
    // RDKit✔️✔️:       Tools::buildNeighborMatrix(mol);
    let neighbor_matrix = build_neighbor_matrix(mol)?;
    // RDKit✔️✔️:   Tools::addNonbonded(mol, confId, params, res.get(), neighborMat, vdwThresh,
    // RDKit✔️✔️:                       ignoreInterfragInteractions);
    add_nonbonded(
        mol,
        params,
        &mut res,
        &neighbor_matrix,
        vdw_thresh,
        ignore_interfrag_interactions,
    )?;
    // RDKit✔️✔️:   Tools::addTorsions(mol, params, res.get());
    add_torsions(mol, params, &mut res, DEFAULT_TORSION_BOND_SMARTS)?;
    // RDKit✔️✔️:   Tools::addInversions(mol, params, res.get());
    add_inversions(mol, params, &mut res)?;

    // RDKit✔️✔️:   return res.release();
    Ok(res)
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::UFF::constructForceField(ROMol&, const AtomicParamVect&, double, int, bool)
}

pub(crate) fn construct_force_field(
    mol: &Molecule,
    total_valences: &[i32],
    hybridizations: &[Hybridization],
    atom_has_conjugated_bond: &[bool],
    vdw_thresh: f64,
    conf_id: isize,
    ignore_interfrag_interactions: bool,
) -> Result<ForceField, UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::constructForceField(ROMol&, double, int, bool) (Builder.cpp:712-718)
    // RDKit✔️❌: ForceFields::ForceField *constructForceField(ROMol &mol, double vdwThresh,
    // RDKit✔️❌:                                              int confId,
    // RDKit✔️❌:                                              bool ignoreInterfragInteractions) {
    // RDKit✔️❌:   bool foundAll;
    // RDKit✔️❌:   AtomicParamVect params;
    // RDKit✔️❌:   boost::tie(params, foundAll) = getAtomTypes(mol);
    let (params, _found_all) = get_atom_types_for_uff(mol, total_valences, hybridizations, atom_has_conjugated_bond)?;
    // RDKit✔️❌:   return constructForceField(mol, params, vdwThresh, confId,
    // RDKit✔️❌:                              ignoreInterfragInteractions);
    construct_force_field_with_params(mol, &params, vdw_thresh, conf_id, ignore_interfrag_interactions)
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION RDKit::UFF::constructForceField(ROMol&, double, int, bool)
}

#[cfg(test)]
mod tests {
    use super::{
        DEFAULT_TORSION_BOND_SMARTS, RELATION_1_2, RELATION_1_3, RELATION_1_4, RELATION_1_X, UffAtomTyperError,
        UffBuilderError, add_angles, add_bonds, add_inversions, add_nonbonded, add_torsions, build_neighbor_matrix,
        construct_force_field, construct_force_field_with_params, get_two_bit_cell, two_bit_cell_pos,
    };
    use crate::chemistry::forcefield::core::{ForceField, ForceFieldVec3};
    use crate::chemistry::forcefield::uff::angle_bend::AngleBendContrib;
    use crate::chemistry::forcefield::uff::inversion::InversionContrib;
    use crate::chemistry::forcefield::uff::nonbonded::VdwContrib;
    use crate::chemistry::forcefield::uff::params::AtomicParams;
    use crate::chemistry::forcefield::uff::torsion_angle::TorsionAngleContrib;
    use crate::chemistry::forcefield::uff::utils::{calc_bond_force_constant, calc_bond_rest_length};
    use crate::chemistry::valence::ValenceError;
    use crate::{AtomSpec, BondOrder, BondSpec, Conformer3D, Element, Hybridization, Molecule, MoleculeBuilder};

    fn carbon() -> AtomSpec {
        AtomSpec::new(Element::C)
    }

    fn carbon_with_hybridization(hybridization: Hybridization) -> AtomSpec {
        AtomSpec::new(Element::C).with_hybridization(hybridization)
    }

    fn atom_with_element_and_hybridization(element: Element, hybridization: Hybridization) -> AtomSpec {
        AtomSpec::new(element).with_hybridization(hybridization)
    }

    fn linear_chain(len: usize) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..len).map(|_| builder.add_atom(carbon())).collect();
        for pair in atoms.windows(2) {
            builder
                .add_bond(BondSpec::new(pair[0], pair[1], BondOrder::Single))
                .expect("linear chain bond should build");
        }
        builder.build().expect("linear chain should build")
    }

    fn single_bond_molecule(order: BondOrder) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(carbon());
        let a1 = builder.add_atom(carbon());
        builder
            .add_bond(BondSpec::new(a0, a1, order))
            .expect("single bond molecule should build");
        builder.build().expect("single bond molecule should build")
    }

    fn force_field_with_positions(n_atoms: usize) -> ForceField {
        let mut ff = ForceField::new(3);
        for i in 0..n_atoms {
            ff.positions_mut().push(ForceFieldVec3::new(i as f64, 0.0, 0.0));
        }
        ff
    }

    fn force_field_from_positions(positions: &[[f64; 3]]) -> ForceField {
        let mut ff = ForceField::new(3);
        for position in positions {
            ff.positions_mut()
                .push(ForceFieldVec3::new(position[0], position[1], position[2]));
        }
        ff
    }

    fn atomic_params(r1: f64, gmp_xi: f64, z1: f64) -> AtomicParams {
        AtomicParams {
            r1,
            theta0: 0.0,
            x1: 0.0,
            d1: 0.0,
            zeta: 0.0,
            z1,
            v1: 0.0,
            u1: 0.0,
            gmp_xi,
            gmp_hardness: 0.0,
            gmp_radius: 0.0,
        }
    }

    fn torsion_params(v1: f64, u1: f64) -> AtomicParams {
        AtomicParams {
            r1: 0.7,
            theta0: 0.0,
            x1: 0.0,
            d1: 0.0,
            zeta: 0.0,
            z1: 1.0,
            v1,
            u1,
            gmp_xi: 0.0,
            gmp_hardness: 0.0,
            gmp_radius: 0.0,
        }
    }

    fn torsion_chain(atoms: [(Element, Hybridization); 4], bond_orders: [BondOrder; 3]) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let atom_ids: Vec<_> = atoms
            .into_iter()
            .map(|(element, hybridization)| {
                builder.add_atom(atom_with_element_and_hybridization(element, hybridization))
            })
            .collect();
        for (pair, order) in atom_ids.windows(2).zip(bond_orders) {
            builder
                .add_bond(BondSpec::new(pair[0], pair[1], order))
                .expect("torsion chain bond should build");
        }
        builder.build().expect("torsion chain should build")
    }

    fn branched_torsion_molecule() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(carbon());
        let a1 = builder.add_atom(carbon_with_hybridization(Hybridization::Sp3));
        let a2 = builder.add_atom(carbon_with_hybridization(Hybridization::Sp3));
        let a3 = builder.add_atom(carbon());
        let a4 = builder.add_atom(carbon());
        let a5 = builder.add_atom(carbon());

        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .expect("bond 0-1 should build");
        builder
            .add_bond(BondSpec::new(a1, a2, BondOrder::Single))
            .expect("bond 1-2 should build");
        builder
            .add_bond(BondSpec::new(a2, a3, BondOrder::Single))
            .expect("bond 2-3 should build");
        builder
            .add_bond(BondSpec::new(a4, a1, BondOrder::Single))
            .expect("bond 4-1 should build");
        builder
            .add_bond(BondSpec::new(a2, a5, BondOrder::Single))
            .expect("bond 2-5 should build");

        builder.build().expect("branched torsion molecule should build")
    }

    fn trigonal_center_molecule(
        center: (Element, Hybridization),
        neighbors: [(Element, Hybridization); 3],
        bond_orders: [BondOrder; 3],
    ) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let center_id = builder.add_atom(atom_with_element_and_hybridization(center.0, center.1));
        let neighbor_ids: Vec<_> = neighbors
            .into_iter()
            .map(|(element, hybridization)| {
                builder.add_atom(atom_with_element_and_hybridization(element, hybridization))
            })
            .collect();
        for (neighbor_id, order) in neighbor_ids.into_iter().zip(bond_orders) {
            builder
                .add_bond(BondSpec::new(center_id, neighbor_id, order))
                .expect("trigonal center bond should build");
        }
        builder.build().expect("trigonal center molecule should build")
    }

    fn vdw_params(x1: f64, d1: f64) -> AtomicParams {
        AtomicParams {
            r1: 0.0,
            theta0: 0.0,
            x1,
            d1,
            zeta: 0.0,
            z1: 0.0,
            v1: 0.0,
            u1: 0.0,
            gmp_xi: 0.0,
            gmp_hardness: 0.0,
            gmp_radius: 0.0,
        }
    }

    fn disconnected_pair() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(carbon());
        builder.add_atom(carbon());
        builder.build().expect("disconnected pair should build")
    }

    fn disconnected_pair_with_3d(coords: [[f64; 3]; 2]) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(carbon());
        builder.add_atom(carbon());
        builder
            .add_3d_conformer(coords.to_vec())
            .expect("3d conformer should build");
        builder.build().expect("disconnected pair should build")
    }

    fn pair_with_elements_and_3d(first: Element, second: Element, coords: [[f64; 3]; 2]) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(first));
        builder.add_atom(AtomSpec::new(second));
        builder
            .add_3d_conformer(coords.to_vec())
            .expect("3d conformer should build");
        builder.build().expect("pair should build")
    }

    fn single_bond_molecule_with_3d(order: BondOrder, coords: [[f64; 3]; 2]) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(carbon());
        let a1 = builder.add_atom(carbon());
        builder
            .add_bond(BondSpec::new(a0, a1, order))
            .expect("single bond molecule should build");
        builder
            .add_3d_conformer(coords.to_vec())
            .expect("3d conformer should build");
        builder.build().expect("single bond molecule should build")
    }

    fn molecule_with_named_conformers(order: BondOrder) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(carbon());
        let a1 = builder.add_atom(carbon());
        builder
            .add_bond(BondSpec::new(a0, a1, order))
            .expect("bond should build");
        builder
            .add_conformer(Conformer3D::new(7, vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], true))
            .expect("first named conformer should build");
        builder
            .add_conformer(Conformer3D::new(9, vec![[10.0, 11.0, 12.0], [13.0, 14.0, 15.0]], true))
            .expect("second named conformer should build");
        builder.build().expect("named-conformer molecule should build")
    }

    fn trigonal_bipyramidal_molecule_with_3d() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::P).with_hybridization(Hybridization::Sp3d));
        let a1 = builder.add_atom(carbon());
        let a2 = builder.add_atom(carbon());
        let a3 = builder.add_atom(carbon());
        let a4 = builder.add_atom(carbon());
        let a5 = builder.add_atom(carbon());
        for neighbor in [a1, a2, a3, a4, a5] {
            builder
                .add_bond(BondSpec::new(center, neighbor, BondOrder::Single))
                .expect("trigonal bipyramidal bond should build");
        }
        builder
            .add_3d_conformer(vec![
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
                [1.0, 0.0, 0.0],
                [-0.5, 0.866_025_403_784_438_6, 0.0],
                [-0.5, -0.866_025_403_784_438_6, 0.0],
            ])
            .expect("3d conformer should build");
        builder.build().expect("trigonal bipyramidal molecule should build")
    }

    fn branch_coverage_molecule() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(carbon());
        let a1 = builder.add_atom(carbon());
        let a2 = builder.add_atom(carbon());
        let a3 = builder.add_atom(carbon());
        let a4 = builder.add_atom(carbon());
        let a5 = builder.add_atom(carbon());

        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .expect("bond 0-1 should build");
        builder
            .add_bond(BondSpec::new(a0, a2, BondOrder::Single))
            .expect("bond 0-2 should build");
        builder
            .add_bond(BondSpec::new(a3, a0, BondOrder::Single))
            .expect("bond 3-0 should build");
        builder
            .add_bond(BondSpec::new(a1, a4, BondOrder::Single))
            .expect("bond 1-4 should build");
        builder
            .add_bond(BondSpec::new(a5, a1, BondOrder::Single))
            .expect("bond 5-1 should build");

        builder.build().expect("branch coverage molecule should build")
    }

    fn three_atom_angle(
        center_hybridization: Hybridization,
        left_order: BondOrder,
        right_order: BondOrder,
    ) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let left = builder.add_atom(carbon());
        let center = builder.add_atom(carbon_with_hybridization(center_hybridization));
        let right = builder.add_atom(carbon());
        builder
            .add_bond(BondSpec::new(left, center, left_order))
            .expect("left bond should build");
        builder
            .add_bond(BondSpec::new(center, right, right_order))
            .expect("right bond should build");
        builder.build().expect("three-atom angle should build")
    }

    fn star_molecule(center_hybridization: Hybridization, degree: usize) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(carbon_with_hybridization(center_hybridization));
        for _ in 0..degree {
            let outer = builder.add_atom(carbon());
            builder
                .add_bond(BondSpec::new(center, outer, BondOrder::Single))
                .expect("star bond should build");
        }
        builder.build().expect("star molecule should build")
    }

    fn triangle_ring() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        let b = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        let c = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        builder
            .add_bond(BondSpec::new(a, b, BondOrder::Single))
            .expect("bond a-b should build");
        builder
            .add_bond(BondSpec::new(b, c, BondOrder::Single))
            .expect("bond b-c should build");
        builder
            .add_bond(BondSpec::new(c, a, BondOrder::Single))
            .expect("bond c-a should build");
        builder.build().expect("triangle ring should build")
    }

    fn triangle_ring_with_substituent() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        let b = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        let c = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        let d = builder.add_atom(carbon());
        builder
            .add_bond(BondSpec::new(a, b, BondOrder::Single))
            .expect("bond a-b should build");
        builder
            .add_bond(BondSpec::new(b, c, BondOrder::Single))
            .expect("bond b-c should build");
        builder
            .add_bond(BondSpec::new(c, a, BondOrder::Single))
            .expect("bond c-a should build");
        builder
            .add_bond(BondSpec::new(b, d, BondOrder::Single))
            .expect("bond b-d should build");
        builder.build().expect("triangle ring with substituent should build")
    }

    fn square_ring() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        let b = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        let c = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        let d = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        builder
            .add_bond(BondSpec::new(a, b, BondOrder::Single))
            .expect("bond a-b should build");
        builder
            .add_bond(BondSpec::new(b, c, BondOrder::Single))
            .expect("bond b-c should build");
        builder
            .add_bond(BondSpec::new(c, d, BondOrder::Single))
            .expect("bond c-d should build");
        builder
            .add_bond(BondSpec::new(d, a, BondOrder::Single))
            .expect("bond d-a should build");
        builder.build().expect("square ring should build")
    }

    fn square_ring_with_substituent() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        let b = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        let c = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        let d = builder.add_atom(carbon_with_hybridization(Hybridization::Sp2));
        let e = builder.add_atom(carbon());
        builder
            .add_bond(BondSpec::new(a, b, BondOrder::Single))
            .expect("bond a-b should build");
        builder
            .add_bond(BondSpec::new(b, c, BondOrder::Single))
            .expect("bond b-c should build");
        builder
            .add_bond(BondSpec::new(c, d, BondOrder::Single))
            .expect("bond c-d should build");
        builder
            .add_bond(BondSpec::new(d, a, BondOrder::Single))
            .expect("bond d-a should build");
        builder
            .add_bond(BondSpec::new(b, e, BondOrder::Single))
            .expect("bond b-e should build");
        builder.build().expect("square ring with substituent should build")
    }

    #[test]
    fn uff_builder_two_bit_cell_pos_matches_rdkit_upper_triangle_layout() {
        let n_atoms = 4;
        assert_eq!(two_bit_cell_pos(n_atoms, 1, 1), 4);
        assert_eq!(two_bit_cell_pos(n_atoms, 2, 1), 5);
        assert_eq!(two_bit_cell_pos(n_atoms, 1, 2), 5);
    }

    #[test]
    fn uff_builder_build_neighbor_matrix_marks_linear_chain_relations() {
        let mol = linear_chain(4);
        let nbr_mat = build_neighbor_matrix(&mol).expect("linear chain matrix should build");

        assert_eq!(nbr_mat.len(), 3);
        assert_eq!(get_two_bit_cell(&nbr_mat, 0), RELATION_1_X);
        assert_eq!(get_two_bit_cell(&nbr_mat, 1), RELATION_1_2);
        assert_eq!(get_two_bit_cell(&nbr_mat, 2), RELATION_1_3);
        assert_eq!(get_two_bit_cell(&nbr_mat, 3), RELATION_1_X);
        assert_eq!(get_two_bit_cell(&nbr_mat, 4), RELATION_1_X);
        assert_eq!(get_two_bit_cell(&nbr_mat, 5), RELATION_1_2);
        assert_eq!(get_two_bit_cell(&nbr_mat, 6), RELATION_1_3);
        assert_eq!(get_two_bit_cell(&nbr_mat, 7), RELATION_1_X);
        assert_eq!(get_two_bit_cell(&nbr_mat, 8), RELATION_1_2);
        assert_eq!(get_two_bit_cell(&nbr_mat, 9), RELATION_1_X);
        assert_eq!(
            get_two_bit_cell(&nbr_mat, two_bit_cell_pos(mol.num_atoms(), 0, 3)),
            RELATION_1_X
        );
        assert_ne!(
            get_two_bit_cell(&nbr_mat, two_bit_cell_pos(mol.num_atoms(), 0, 3)),
            RELATION_1_4
        );
    }

    #[test]
    fn uff_builder_build_neighbor_matrix_covers_all_shared_endpoint_branches() {
        let mol = branch_coverage_molecule();
        let nbr_mat = build_neighbor_matrix(&mol).expect("branch coverage matrix should build");

        assert_eq!(
            get_two_bit_cell(&nbr_mat, two_bit_cell_pos(mol.num_atoms(), 0, 1)),
            RELATION_1_2
        );
        assert_eq!(
            get_two_bit_cell(&nbr_mat, two_bit_cell_pos(mol.num_atoms(), 0, 2)),
            RELATION_1_2
        );
        assert_eq!(
            get_two_bit_cell(&nbr_mat, two_bit_cell_pos(mol.num_atoms(), 0, 3)),
            RELATION_1_2
        );
        assert_eq!(
            get_two_bit_cell(&nbr_mat, two_bit_cell_pos(mol.num_atoms(), 1, 4)),
            RELATION_1_2
        );
        assert_eq!(
            get_two_bit_cell(&nbr_mat, two_bit_cell_pos(mol.num_atoms(), 1, 5)),
            RELATION_1_2
        );

        assert_eq!(
            get_two_bit_cell(&nbr_mat, two_bit_cell_pos(mol.num_atoms(), 1, 2)),
            RELATION_1_3
        );
        assert_eq!(
            get_two_bit_cell(&nbr_mat, two_bit_cell_pos(mol.num_atoms(), 1, 3)),
            RELATION_1_3
        );
        assert_eq!(
            get_two_bit_cell(&nbr_mat, two_bit_cell_pos(mol.num_atoms(), 0, 4)),
            RELATION_1_3
        );
        assert_eq!(
            get_two_bit_cell(&nbr_mat, two_bit_cell_pos(mol.num_atoms(), 0, 5)),
            RELATION_1_3
        );
        assert_eq!(
            get_two_bit_cell(&nbr_mat, two_bit_cell_pos(mol.num_atoms(), 2, 4)),
            RELATION_1_X
        );
    }

    #[test]
    fn uff_builder_build_neighbor_matrix_rejects_empty_molecule() {
        let err =
            build_neighbor_matrix(&Molecule::default()).expect_err("empty molecule should stay explicitly unsupported");
        assert_eq!(err, UffBuilderError::EmptyMolecule);
    }

    #[test]
    fn uff_builder_add_bonds_adds_one_bond_stretch_per_parameterized_bond() {
        let mol = single_bond_molecule(BondOrder::Single);
        let mut ff = force_field_with_positions(mol.num_atoms());
        let end1 = atomic_params(0.757, 5.343, 1.912);
        let end2 = atomic_params(0.658, 8.741, 2.3);
        let params = vec![Some(end1.clone()), Some(end2.clone())];

        add_bonds(&mol, &params, &mut ff).expect("single bond should add one contrib");
        ff.initialize();

        assert_eq!(ff.contribs().len(), 1);
        let rest_len = calc_bond_rest_length(1.0, &end1, &end2).expect("valid rest length");
        let force_constant = calc_bond_force_constant(rest_len, &end1, &end2).expect("valid force constant");
        let pos = [0.0, 0.0, 0.0, rest_len + 0.25, 0.0, 0.0];
        let expected = 0.5 * force_constant * 0.25_f64 * 0.25_f64;
        assert!((ff.contribs()[0].get_energy(&pos) - expected).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_bonds_skips_bonds_with_missing_params() {
        let mol = linear_chain(3);
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![
            Some(atomic_params(0.757, 5.343, 1.912)),
            None,
            Some(atomic_params(0.658, 8.741, 2.3)),
        ];

        add_bonds(&mol, &params, &mut ff).expect("missing params should skip bond contribs");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_add_bonds_rejects_parameter_length_mismatch() {
        let mol = linear_chain(2);
        let mut ff = force_field_with_positions(mol.num_atoms());
        let err = add_bonds(&mol, &[Some(atomic_params(0.757, 5.343, 1.912))], &mut ff)
            .expect_err("length mismatch should error");

        assert_eq!(err, UffBuilderError::ParamsLengthMismatch { atoms: 2, params: 1 });
    }

    #[test]
    fn uff_builder_add_bonds_propagates_bad_bond_type() {
        let mol = single_bond_molecule(BondOrder::DativeLeft);
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![
            Some(atomic_params(0.757, 5.343, 1.912)),
            Some(atomic_params(0.658, 8.741, 2.3)),
        ];

        let err = add_bonds(&mol, &params, &mut ff).expect_err("unsupported bond type should fail");

        assert!(matches!(
            err,
            UffBuilderError::Valence(ValenceError::BadBondType { .. })
        ));
    }

    #[test]
    fn uff_builder_add_angles_adds_one_sp_angle_contrib() {
        let mol = three_atom_angle(Hybridization::Sp, BondOrder::Single, BondOrder::Single);
        let positions = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mut ff = force_field_from_positions(&positions);
        let params = vec![
            Some(atomic_params(0.757, 5.343, 1.912)),
            Some(atomic_params(0.658, 8.741, 2.3)),
            Some(atomic_params(0.701, 6.1, 1.8)),
        ];

        add_angles(&mol, &params, &mut ff).expect("sp angle should add one contrib");
        ff.initialize();

        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let expected = AngleBendContrib::new(
            &ff,
            0,
            1,
            2,
            1.0,
            1.0,
            params[0].as_ref().expect("left params should exist"),
            params[1].as_ref().expect("center params should exist"),
            params[2].as_ref().expect("right params should exist"),
            1,
        )
        .expect("direct sp angle contrib should build");

        assert_eq!(ff.contribs().len(), 1);
        assert!((ff.contribs()[0].get_energy(&pos) - expected.get_energy(&pos)).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_angles_uses_default_order_for_non_special_center() {
        let mol = three_atom_angle(Hybridization::Sp3, BondOrder::Single, BondOrder::Single);
        let positions = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mut ff = force_field_from_positions(&positions);
        let params = vec![
            Some(atomic_params(0.757, 5.343, 1.912)),
            Some(atomic_params(0.658, 8.741, 2.3)),
            Some(atomic_params(0.701, 6.1, 1.8)),
        ];

        add_angles(&mol, &params, &mut ff).expect("sp3 angle should add one contrib");
        ff.initialize();

        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let expected = AngleBendContrib::new(
            &ff,
            0,
            1,
            2,
            1.0,
            1.0,
            params[0].as_ref().expect("left params should exist"),
            params[1].as_ref().expect("center params should exist"),
            params[2].as_ref().expect("right params should exist"),
            0,
        )
        .expect("direct default-order angle contrib should build");

        assert_eq!(ff.contribs().len(), 1);
        assert!((ff.contribs()[0].get_energy(&pos) - expected.get_energy(&pos)).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_angles_uses_sp3d2_order_four() {
        let mol = three_atom_angle(Hybridization::Sp3d2, BondOrder::Single, BondOrder::Single);
        let positions = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mut ff = force_field_from_positions(&positions);
        let params = vec![
            Some(atomic_params(0.757, 5.343, 1.912)),
            Some(atomic_params(0.658, 8.741, 2.3)),
            Some(atomic_params(0.701, 6.1, 1.8)),
        ];

        add_angles(&mol, &params, &mut ff).expect("sp3d2 angle should add one contrib");
        ff.initialize();

        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let expected = AngleBendContrib::new(
            &ff,
            0,
            1,
            2,
            1.0,
            1.0,
            params[0].as_ref().expect("left params should exist"),
            params[1].as_ref().expect("center params should exist"),
            params[2].as_ref().expect("right params should exist"),
            4,
        )
        .expect("direct sp3d2 angle contrib should build");

        assert_eq!(ff.contribs().len(), 1);
        assert!((ff.contribs()[0].get_energy(&pos) - expected.get_energy(&pos)).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_angles_counts_unique_neighbor_pairs_once() {
        let mol = star_molecule(Hybridization::Sp3, 3);
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![Some(atomic_params(0.658, 8.741, 2.3)); mol.num_atoms()];

        add_angles(&mol, &params, &mut ff).expect("star center should add unique angle pairs");

        assert_eq!(ff.contribs().len(), 3);
    }

    #[test]
    fn uff_builder_add_angles_skips_degree_one_and_missing_params() {
        let degree_one = single_bond_molecule(BondOrder::Single);
        let mut degree_one_ff = force_field_with_positions(degree_one.num_atoms());
        let degree_one_params = vec![
            Some(atomic_params(0.757, 5.343, 1.912)),
            Some(atomic_params(0.658, 8.741, 2.3)),
        ];
        add_angles(&degree_one, &degree_one_params, &mut degree_one_ff).expect("degree-one atoms should simply skip");
        assert_eq!(degree_one_ff.contribs().len(), 0);

        let missing_center = three_atom_angle(Hybridization::Sp3, BondOrder::Single, BondOrder::Single);
        let mut missing_center_ff = force_field_with_positions(missing_center.num_atoms());
        let missing_center_params = vec![
            Some(atomic_params(0.757, 5.343, 1.912)),
            None,
            Some(atomic_params(0.701, 6.1, 1.8)),
        ];
        add_angles(&missing_center, &missing_center_params, &mut missing_center_ff)
            .expect("missing center params should skip");
        assert_eq!(missing_center_ff.contribs().len(), 0);

        let missing_end = three_atom_angle(Hybridization::Sp3, BondOrder::Single, BondOrder::Single);
        let mut missing_end_ff = force_field_with_positions(missing_end.num_atoms());
        let missing_end_params = vec![
            Some(atomic_params(0.757, 5.343, 1.912)),
            Some(atomic_params(0.658, 8.741, 2.3)),
            None,
        ];
        add_angles(&missing_end, &missing_end_params, &mut missing_end_ff)
            .expect("missing endpoint params should skip");
        assert_eq!(missing_end_ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_add_angles_skips_sp3d_degree_five_special_case() {
        let mol = star_molecule(Hybridization::Sp3d, 5);
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![Some(atomic_params(0.658, 8.741, 2.3)); mol.num_atoms()];

        add_angles(&mol, &params, &mut ff).expect("special sp3d degree-5 case should skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_add_angles_rejects_parameter_length_mismatch() {
        let mol = three_atom_angle(Hybridization::Sp3, BondOrder::Single, BondOrder::Single);
        let mut ff = force_field_with_positions(mol.num_atoms());

        let err = add_angles(&mol, &[Some(atomic_params(0.757, 5.343, 1.912))], &mut ff)
            .expect_err("length mismatch should error");

        assert_eq!(err, UffBuilderError::ParamsLengthMismatch { atoms: 3, params: 1 });
    }

    #[test]
    fn uff_builder_add_angles_propagates_bad_bond_type() {
        let mol = three_atom_angle(Hybridization::Sp3, BondOrder::DativeLeft, BondOrder::Single);
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![
            Some(atomic_params(0.757, 5.343, 1.912)),
            Some(atomic_params(0.658, 8.741, 2.3)),
            Some(atomic_params(0.701, 6.1, 1.8)),
        ];

        let err = add_angles(&mol, &params, &mut ff).expect_err("unsupported bond type should fail");

        assert!(matches!(
            err,
            UffBuilderError::Valence(ValenceError::BadBondType { .. })
        ));
    }

    #[test]
    fn uff_builder_add_angles_applies_ring3_substituent_hack_order30() {
        let mol = triangle_ring_with_substituent();
        let positions = [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [-0.5, 0.866_025_403_784_438_6, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let mut ff = force_field_from_positions(&positions);
        let params = vec![
            Some(atomic_params(0.757, 5.343, 1.912)),
            Some(atomic_params(0.658, 8.741, 2.3)),
            None,
            Some(atomic_params(0.701, 6.1, 1.8)),
        ];

        add_angles(&mol, &params, &mut ff).expect("ring-3 substituent angle should build");
        ff.initialize();

        let pos = [
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -0.5,
            0.866_025_403_784_438_6,
            0.0,
            0.0,
            1.0,
            0.0,
        ];
        let expected = AngleBendContrib::new(
            &ff,
            0,
            1,
            3,
            1.0,
            1.0,
            params[0].as_ref().expect("ring endpoint params should exist"),
            params[1].as_ref().expect("center params should exist"),
            params[3].as_ref().expect("substituent params should exist"),
            30,
        )
        .expect("direct order30 contrib should build");

        assert_eq!(ff.contribs().len(), 1);
        assert!((ff.contribs()[0].get_energy(&pos) - expected.get_energy(&pos)).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_angles_applies_ring3_internal_hack_order35() {
        let mol = triangle_ring();
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.866_025_403_784_438_6, 0.0]];
        let mut ff = force_field_from_positions(&positions);
        let params = vec![Some(atomic_params(0.658, 8.741, 2.3)); mol.num_atoms()];

        add_angles(&mol, &params, &mut ff).expect("triangle ring angles should build");
        ff.initialize();

        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.866_025_403_784_438_6, 0.0];
        let expected_total = AngleBendContrib::new(
            &ff,
            0,
            1,
            2,
            1.0,
            1.0,
            params[0].as_ref().expect("params should exist"),
            params[1].as_ref().expect("params should exist"),
            params[2].as_ref().expect("params should exist"),
            35,
        )
        .expect("order35 contrib should build")
        .get_energy(&pos)
            + AngleBendContrib::new(
                &ff,
                1,
                2,
                0,
                1.0,
                1.0,
                params[1].as_ref().expect("params should exist"),
                params[2].as_ref().expect("params should exist"),
                params[0].as_ref().expect("params should exist"),
                35,
            )
            .expect("order35 contrib should build")
            .get_energy(&pos)
            + AngleBendContrib::new(
                &ff,
                2,
                0,
                1,
                1.0,
                1.0,
                params[2].as_ref().expect("params should exist"),
                params[0].as_ref().expect("params should exist"),
                params[1].as_ref().expect("params should exist"),
                35,
            )
            .expect("order35 contrib should build")
            .get_energy(&pos);
        let actual_total = ff
            .contribs()
            .iter()
            .map(|contrib| contrib.get_energy(&pos))
            .sum::<f64>();

        assert_eq!(ff.contribs().len(), 3);
        assert!((actual_total - expected_total).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_angles_applies_ring4_substituent_hack_order40() {
        let mol = square_ring_with_substituent();
        let positions = [
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let mut ff = force_field_from_positions(&positions);
        let params = vec![
            Some(atomic_params(0.757, 5.343, 1.912)),
            Some(atomic_params(0.658, 8.741, 2.3)),
            None,
            None,
            Some(atomic_params(0.701, 6.1, 1.8)),
        ];

        add_angles(&mol, &params, &mut ff).expect("ring-4 substituent angle should build");
        ff.initialize();

        let pos = [
            -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
        ];
        let expected = AngleBendContrib::new(
            &ff,
            0,
            1,
            4,
            1.0,
            1.0,
            params[0].as_ref().expect("ring endpoint params should exist"),
            params[1].as_ref().expect("center params should exist"),
            params[4].as_ref().expect("substituent params should exist"),
            40,
        )
        .expect("direct order40 contrib should build");

        assert_eq!(ff.contribs().len(), 1);
        assert!((ff.contribs()[0].get_energy(&pos) - expected.get_energy(&pos)).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_angles_applies_ring4_internal_hack_order45() {
        let mol = square_ring();
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
        let mut ff = force_field_from_positions(&positions);
        let params = vec![Some(atomic_params(0.658, 8.741, 2.3)); mol.num_atoms()];

        add_angles(&mol, &params, &mut ff).expect("square ring angles should build");
        ff.initialize();

        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0];
        let expected_total = AngleBendContrib::new(
            &ff,
            0,
            1,
            2,
            1.0,
            1.0,
            params[0].as_ref().expect("params should exist"),
            params[1].as_ref().expect("params should exist"),
            params[2].as_ref().expect("params should exist"),
            45,
        )
        .expect("order45 contrib should build")
        .get_energy(&pos)
            + AngleBendContrib::new(
                &ff,
                1,
                2,
                3,
                1.0,
                1.0,
                params[1].as_ref().expect("params should exist"),
                params[2].as_ref().expect("params should exist"),
                params[3].as_ref().expect("params should exist"),
                45,
            )
            .expect("order45 contrib should build")
            .get_energy(&pos)
            + AngleBendContrib::new(
                &ff,
                2,
                3,
                0,
                1.0,
                1.0,
                params[2].as_ref().expect("params should exist"),
                params[3].as_ref().expect("params should exist"),
                params[0].as_ref().expect("params should exist"),
                45,
            )
            .expect("order45 contrib should build")
            .get_energy(&pos)
            + AngleBendContrib::new(
                &ff,
                3,
                0,
                1,
                1.0,
                1.0,
                params[3].as_ref().expect("params should exist"),
                params[0].as_ref().expect("params should exist"),
                params[1].as_ref().expect("params should exist"),
                45,
            )
            .expect("order45 contrib should build")
            .get_energy(&pos);
        let actual_total = ff
            .contribs()
            .iter()
            .map(|contrib| contrib.get_energy(&pos))
            .sum::<f64>();

        assert_eq!(ff.contribs().len(), 4);
        assert!((actual_total - expected_total).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_torsions_rejects_parameter_length_mismatch() {
        let mol = linear_chain(4);
        let mut ff = force_field_with_positions(mol.num_atoms());

        let err = add_torsions(
            &mol,
            &[Some(torsion_params(2.0, 1.5))],
            &mut ff,
            DEFAULT_TORSION_BOND_SMARTS,
        )
        .expect_err("length mismatch should error");

        assert_eq!(err, UffBuilderError::ParamsLengthMismatch { atoms: 4, params: 1 });
    }

    #[test]
    fn uff_builder_add_torsions_accepts_custom_bond_smarts() {
        let mol = torsion_chain(
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![Some(torsion_params(2.0, 1.5)); mol.num_atoms()];

        add_torsions(&mol, &params, &mut ff, "[*:1]-[*:2]")
            .expect("custom torsion SMARTS should select the three single bonds");

        assert_eq!(ff.contribs().len(), 1);
    }

    #[test]
    fn uff_builder_add_torsions_adds_one_sp3_sp3_contrib_for_default_query_hit() {
        let mol = torsion_chain(
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 1.0]];
        let mut ff = force_field_from_positions(&positions);
        let params = vec![
            Some(torsion_params(0.0, 0.0)),
            Some(torsion_params(4.0, 1.0)),
            Some(torsion_params(9.0, 1.0)),
            Some(torsion_params(0.0, 0.0)),
        ];

        add_torsions(&mol, &params, &mut ff, DEFAULT_TORSION_BOND_SMARTS)
            .expect("default torsion query hit should add a contrib");
        ff.initialize();

        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 2.0, 1.0, 1.0];
        let expected = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            3,
            1.0,
            i32::from(mol.atoms()[1].atomic_number()),
            i32::from(mol.atoms()[2].atomic_number()),
            mol.atoms()[1].hybridization(),
            mol.atoms()[2].hybridization(),
            params[1].as_ref().expect("central params should exist"),
            params[2].as_ref().expect("central params should exist"),
            false,
        );

        assert_eq!(ff.contribs().len(), 1);
        assert!((ff.contribs()[0].get_energy(&pos) - expected.get_energy(&pos)).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_torsions_skips_default_query_hits_with_triple_adjacent_atom() {
        let mol = torsion_chain(
            [
                (Element::C, Hybridization::Sp),
                (Element::C, Hybridization::Sp),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Triple, BondOrder::Single, BondOrder::Single],
        );
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![Some(torsion_params(4.0, 1.0)); mol.num_atoms()];

        add_torsions(&mol, &params, &mut ff, DEFAULT_TORSION_BOND_SMARTS)
            .expect("triple-adjacent default-query miss should skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_add_torsions_skips_bonds_with_missing_central_params() {
        let mol = torsion_chain(
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![
            Some(torsion_params(0.0, 0.0)),
            None,
            Some(torsion_params(9.0, 1.0)),
            Some(torsion_params(0.0, 0.0)),
        ];

        add_torsions(&mol, &params, &mut ff, DEFAULT_TORSION_BOND_SMARTS)
            .expect("missing central params should skip the torsion bond");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_add_torsions_skips_non_sp2_sp3_central_atoms() {
        let mol = torsion_chain(
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3d2),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![Some(torsion_params(4.0, 1.0)); mol.num_atoms()];

        add_torsions(&mol, &params, &mut ff, DEFAULT_TORSION_BOND_SMARTS)
            .expect("unsupported central hybridization should skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_add_torsions_skips_three_membered_ring_torsions() {
        let mol = triangle_ring();
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![Some(torsion_params(4.0, 1.5)); mol.num_atoms()];

        add_torsions(&mol, &params, &mut ff, DEFAULT_TORSION_BOND_SMARTS)
            .expect("three-membered ring torsions should skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_add_torsions_propagates_bad_central_bond_type() {
        let mol = torsion_chain(
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::DativeLeft, BondOrder::Single],
        );
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![Some(torsion_params(4.0, 1.0)); mol.num_atoms()];

        let err = add_torsions(&mol, &params, &mut ff, DEFAULT_TORSION_BOND_SMARTS)
            .expect_err("unsupported central bond type should fail");

        assert!(matches!(
            err,
            UffBuilderError::Valence(ValenceError::BadBondType { .. })
        ));
    }

    #[test]
    fn uff_builder_add_torsions_sets_end_atom_sp2_flag_for_sp2_sp3_case() {
        let mol = torsion_chain(
            [
                (Element::C, Hybridization::Sp2),
                (Element::C, Hybridization::Sp2),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 1.0]];
        let mut ff = force_field_from_positions(&positions);
        let params = vec![
            Some(torsion_params(0.0, 0.0)),
            Some(torsion_params(0.0, 2.25)),
            Some(torsion_params(0.0, 4.0)),
            Some(torsion_params(0.0, 0.0)),
        ];

        add_torsions(&mol, &params, &mut ff, DEFAULT_TORSION_BOND_SMARTS)
            .expect("sp2-sp3 torsion should add a contrib");
        ff.initialize();

        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 2.0, 1.0, 1.0];
        let expected = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            3,
            1.0,
            i32::from(mol.atoms()[1].atomic_number()),
            i32::from(mol.atoms()[2].atomic_number()),
            mol.atoms()[1].hybridization(),
            mol.atoms()[2].hybridization(),
            params[1].as_ref().expect("central params should exist"),
            params[2].as_ref().expect("central params should exist"),
            true,
        );

        assert_eq!(ff.contribs().len(), 1);
        assert!((ff.contribs()[0].get_energy(&pos) - expected.get_energy(&pos)).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_torsions_scales_force_constant_by_contribs_on_bond() {
        let mol = branched_torsion_molecule();
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [3.0, 1.0, 1.0],
            [1.0, -1.0, 1.0],
            [2.0, 2.0, 1.0],
        ];
        let mut ff = force_field_from_positions(&positions);
        let params = vec![
            Some(torsion_params(0.0, 0.0)),
            Some(torsion_params(4.0, 1.0)),
            Some(torsion_params(9.0, 1.0)),
            Some(torsion_params(0.0, 0.0)),
            Some(torsion_params(0.0, 0.0)),
            Some(torsion_params(0.0, 0.0)),
        ];

        add_torsions(&mol, &params, &mut ff, DEFAULT_TORSION_BOND_SMARTS)
            .expect("branched bond should add four torsion contribs");
        ff.initialize();

        let pos = [
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0, 3.0, 1.0, 1.0, 1.0, -1.0, 1.0, 2.0, 2.0, 1.0,
        ];
        let mut expected_total = 0.0;
        for (a, d) in [(0, 3), (0, 5), (4, 3), (4, 5)] {
            let mut contrib = TorsionAngleContrib::new(
                &ff,
                a,
                1,
                2,
                d,
                1.0,
                i32::from(mol.atoms()[1].atomic_number()),
                i32::from(mol.atoms()[2].atomic_number()),
                mol.atoms()[1].hybridization(),
                mol.atoms()[2].hybridization(),
                params[1].as_ref().expect("central params should exist"),
                params[2].as_ref().expect("central params should exist"),
                false,
            );
            contrib.scale_force_constant(4);
            expected_total += contrib.get_energy(&pos);
        }
        let actual_total = ff
            .contribs()
            .iter()
            .map(|contrib| contrib.get_energy(&pos))
            .sum::<f64>();

        assert_eq!(ff.contribs().len(), 4);
        assert!((actual_total - expected_total).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_inversions_rejects_parameter_length_mismatch() {
        let mol = trigonal_center_molecule(
            (Element::C, Hybridization::Sp2),
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let mut ff = force_field_with_positions(mol.num_atoms());

        let err = add_inversions(&mol, &[None], &mut ff).expect_err("length mismatch should fail");

        assert_eq!(err, UffBuilderError::ParamsLengthMismatch { atoms: 4, params: 1 });
    }

    #[test]
    fn uff_builder_add_inversions_adds_three_permuted_contribs_for_sp2_center() {
        let mol = trigonal_center_molecule(
            (Element::C, Hybridization::Sp2),
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let positions = [
            [0.0, 0.0, 0.1],
            [1.0, 0.0, 0.0],
            [-0.5, 0.866_025_403_784_438_6, 0.0],
            [-0.5, -0.866_025_403_784_438_6, 0.0],
        ];
        let mut ff = force_field_from_positions(&positions);
        let params = vec![None; mol.num_atoms()];

        add_inversions(&mol, &params, &mut ff).expect("sp2 trigonal center should add inversions");
        ff.initialize();

        let pos = [
            0.0,
            0.0,
            0.1,
            1.0,
            0.0,
            0.0,
            -0.5,
            0.866_025_403_784_438_6,
            0.0,
            -0.5,
            -0.866_025_403_784_438_6,
            0.0,
        ];
        let expected_total = InversionContrib::new(&ff, 1, 0, 2, 3, 6, false, 1.0).get_energy(&pos)
            + InversionContrib::new(&ff, 1, 0, 3, 2, 6, false, 1.0).get_energy(&pos)
            + InversionContrib::new(&ff, 2, 0, 3, 1, 6, false, 1.0).get_energy(&pos);
        let actual_total = ff
            .contribs()
            .iter()
            .map(|contrib| contrib.get_energy(&pos))
            .sum::<f64>();

        assert_eq!(ff.contribs().len(), 3);
        assert!((actual_total - expected_total).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_inversions_sets_carboxyl_sp2_oxygen_flag() {
        let mol = trigonal_center_molecule(
            (Element::C, Hybridization::Sp2),
            [
                (Element::O, Hybridization::Sp2),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Double, BondOrder::Single, BondOrder::Single],
        );
        let positions = [[0.0, 0.0, 0.15], [1.2, 0.0, 0.0], [-0.6, 0.9, 0.0], [-0.6, -0.9, 0.0]];
        let mut ff = force_field_from_positions(&positions);
        let params = vec![None; mol.num_atoms()];

        add_inversions(&mol, &params, &mut ff).expect("sp2 carbon bound to sp2 oxygen should add inversions");
        ff.initialize();

        let pos = [0.0, 0.0, 0.15, 1.2, 0.0, 0.0, -0.6, 0.9, 0.0, -0.6, -0.9, 0.0];
        let expected = InversionContrib::new(&ff, 1, 0, 2, 3, 6, true, 1.0);

        assert_eq!(ff.contribs().len(), 3);
        assert!((ff.contribs()[0].get_energy(&pos) - expected.get_energy(&pos)).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_inversions_skips_unsupported_central_atomic_number() {
        let mol = trigonal_center_molecule(
            (Element::S, Hybridization::Sp2),
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![None; mol.num_atoms()];

        add_inversions(&mol, &params, &mut ff).expect("unsupported central atom should simply skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_add_inversions_skips_degree_not_equal_three() {
        let mol = three_atom_angle(Hybridization::Sp2, BondOrder::Single, BondOrder::Single);
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![None; mol.num_atoms()];

        add_inversions(&mol, &params, &mut ff).expect("degree-two center should skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_add_inversions_skips_cnosp2_centers_without_sp2_hybridization() {
        let mol = trigonal_center_molecule(
            (Element::N, Hybridization::Sp3),
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![None; mol.num_atoms()];

        add_inversions(&mol, &params, &mut ff).expect("non-sp2 nitrogen center should skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_add_inversions_allows_group15_center_without_sp2_requirement() {
        let mol = trigonal_center_molecule(
            (Element::P, Hybridization::Sp3),
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let mut ff = force_field_with_positions(mol.num_atoms());
        let params = vec![None; mol.num_atoms()];

        add_inversions(&mol, &params, &mut ff).expect("phosphorus center should not require sp2 hybridization");

        assert_eq!(ff.contribs().len(), 3);
    }

    #[test]
    fn uff_builder_add_nonbonded_rejects_parameter_length_mismatch() {
        let mol = disconnected_pair();
        let mut ff = force_field_with_positions(mol.num_atoms());
        let neighbor_matrix = build_neighbor_matrix(&mol).expect("neighbor matrix should build");

        let err = add_nonbonded(
            &mol,
            &[Some(vdw_params(4.0, 0.25))],
            &mut ff,
            &neighbor_matrix,
            1.0,
            false,
        )
        .expect_err("length mismatch should fail");

        assert_eq!(err, UffBuilderError::ParamsLengthMismatch { atoms: 2, params: 1 });
    }

    #[test]
    fn uff_builder_add_nonbonded_adds_one_vdw_contrib_for_14_pair_within_threshold() {
        let mol = linear_chain(4);
        let neighbor_matrix = build_neighbor_matrix(&mol).expect("neighbor matrix should build");
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let mut ff = force_field_from_positions(&positions);
        let params = vec![
            Some(vdw_params(4.0, 0.25)),
            Some(vdw_params(4.0, 0.25)),
            Some(vdw_params(4.0, 0.25)),
            Some(vdw_params(4.0, 0.25)),
        ];

        add_nonbonded(&mol, &params, &mut ff, &neighbor_matrix, 1.0, false)
            .expect("1-4 pair inside threshold should add one contrib");
        ff.initialize();

        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0];
        let expected = VdwContrib::new(
            &ff,
            0,
            3,
            params[0].as_ref().expect("params should exist"),
            params[3].as_ref().expect("params should exist"),
            10.0,
        );

        assert_eq!(ff.contribs().len(), 1);
        assert!((ff.contribs()[0].get_energy(&pos) - expected.get_energy(&pos)).abs() < 1.0e-12);
    }

    #[test]
    fn uff_builder_add_nonbonded_skips_pairs_closer_than_14_relation() {
        let mol = linear_chain(3);
        let neighbor_matrix = build_neighbor_matrix(&mol).expect("neighbor matrix should build");
        let mut ff = force_field_from_positions(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let params = vec![Some(vdw_params(4.0, 0.25)); mol.num_atoms()];

        add_nonbonded(&mol, &params, &mut ff, &neighbor_matrix, 10.0, false)
            .expect("1-2 and 1-3 relations should skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_add_nonbonded_skips_pairs_at_threshold_boundary() {
        let mol = disconnected_pair();
        let neighbor_matrix = build_neighbor_matrix(&mol).expect("neighbor matrix should build");
        let mut ff = force_field_from_positions(&[[0.0, 0.0, 0.0], [6.0, 0.0, 0.0]]);
        let params = vec![Some(vdw_params(4.0, 0.25)), Some(vdw_params(9.0, 0.36))];

        add_nonbonded(&mol, &params, &mut ff, &neighbor_matrix, 1.0, false)
            .expect("distance exactly at threshold should skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_add_nonbonded_skips_pairs_with_missing_params() {
        let mol = disconnected_pair();
        let neighbor_matrix = build_neighbor_matrix(&mol).expect("neighbor matrix should build");
        let mut ff = force_field_from_positions(&[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);
        let params = vec![Some(vdw_params(4.0, 0.25)), None];

        add_nonbonded(&mol, &params, &mut ff, &neighbor_matrix, 1.0, false).expect("missing params should skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_add_nonbonded_ignores_interfragment_interactions_when_requested() {
        let mol = disconnected_pair();
        let neighbor_matrix = build_neighbor_matrix(&mol).expect("neighbor matrix should build");
        let params = vec![Some(vdw_params(4.0, 0.25)), Some(vdw_params(9.0, 0.36))];

        let mut ff_keep = force_field_from_positions(&[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);
        add_nonbonded(&mol, &params, &mut ff_keep, &neighbor_matrix, 1.0, false)
            .expect("cross-fragment interactions should be kept when not ignored");
        assert_eq!(ff_keep.contribs().len(), 1);

        let mut ff_skip = force_field_from_positions(&[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);
        add_nonbonded(&mol, &params, &mut ff_skip, &neighbor_matrix, 1.0, true)
            .expect("cross-fragment interactions should be skipped when requested");
        assert_eq!(ff_skip.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_construct_force_field_with_params_rejects_parameter_length_mismatch() {
        let mol = disconnected_pair_with_3d([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);

        let err = match construct_force_field_with_params(&mol, &[Some(vdw_params(4.0, 9.0))], 1.0, -1, false) {
            Ok(_) => panic!("parameter length mismatch should error"),
            Err(err) => err,
        };

        assert_eq!(err, UffBuilderError::ParamsLengthMismatch { atoms: 2, params: 1 });
    }

    #[test]
    fn uff_builder_construct_force_field_with_params_rejects_invalid_conf_id_below_minus_one() {
        let mol = disconnected_pair_with_3d([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);
        let params = vec![Some(vdw_params(4.0, 9.0)), Some(vdw_params(4.0, 9.0))];

        let err = match construct_force_field_with_params(&mol, &params, 1.0, -2, false) {
            Ok(_) => panic!("conf_id < -1 should error"),
            Err(err) => err,
        };

        assert_eq!(err, UffBuilderError::InvalidConformerId { conf_id: -2 });
    }

    #[test]
    fn uff_builder_construct_force_field_with_params_rejects_missing_3d_conformer() {
        let mol = disconnected_pair();
        let params = vec![Some(vdw_params(4.0, 9.0)), Some(vdw_params(4.0, 9.0))];

        let err = match construct_force_field_with_params(&mol, &params, 1.0, -1, false) {
            Ok(_) => panic!("missing 3d conformer should error"),
            Err(err) => err,
        };

        assert_eq!(err, UffBuilderError::Missing3dConformer { conf_id: -1 });
    }

    #[test]
    fn uff_builder_construct_force_field_with_params_uses_default_conformer_for_minus_one() {
        let mol = single_bond_molecule_with_3d(BondOrder::Single, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        let params = vec![Some(atomic_params(0.7, 1.2, 1.0)), Some(atomic_params(0.8, 1.4, 1.1))];

        let ff = construct_force_field_with_params(&mol, &params, 0.0, -1, false)
            .expect("default conformer should build force field");

        assert_eq!(
            ff.positions(),
            &[ForceFieldVec3::new(1.0, 2.0, 3.0), ForceFieldVec3::new(4.0, 5.0, 6.0)]
        );
        assert_eq!(ff.contribs().len(), 1);
    }

    #[test]
    fn uff_builder_construct_force_field_with_params_selects_requested_conformer_by_id() {
        let mol = molecule_with_named_conformers(BondOrder::Single);
        let params = vec![Some(atomic_params(0.7, 1.2, 1.0)), Some(atomic_params(0.8, 1.4, 1.1))];

        let ff = construct_force_field_with_params(&mol, &params, 0.0, 9, false)
            .expect("named conformer id should build force field");

        assert_eq!(
            ff.positions(),
            &[
                ForceFieldVec3::new(10.0, 11.0, 12.0),
                ForceFieldVec3::new(13.0, 14.0, 15.0)
            ]
        );
    }

    #[test]
    fn uff_builder_construct_force_field_with_params_skips_missing_params() {
        let mol = disconnected_pair_with_3d([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);
        let params = vec![Some(vdw_params(4.0, 9.0)), None];

        let ff = construct_force_field_with_params(&mol, &params, 1.0, -1, false)
            .expect("missing parameter slots should be skipped");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_construct_force_field_with_params_ignores_interfragment_interactions_when_requested() {
        let mol = disconnected_pair_with_3d([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);
        let params = vec![Some(vdw_params(4.0, 9.0)), Some(vdw_params(4.0, 9.0))];

        let ff_keep = construct_force_field_with_params(&mol, &params, 1.0, -1, false)
            .expect("interfragment interactions should be kept");
        let ff_skip = construct_force_field_with_params(&mol, &params, 1.0, -1, true)
            .expect("interfragment interactions should be skipped");

        assert_eq!(ff_keep.contribs().len(), 1);
        assert_eq!(ff_skip.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_construct_force_field_with_params_adds_trigonal_bipyramidal_special_case_angles() {
        let mol = trigonal_bipyramidal_molecule_with_3d();
        let params = vec![
            Some(atomic_params(0.7, 1.0, 1.0)),
            Some(atomic_params(0.7, 1.0, 1.0)),
            Some(atomic_params(0.7, 1.0, 1.0)),
            Some(atomic_params(0.7, 1.0, 1.0)),
            Some(atomic_params(0.7, 1.0, 1.0)),
            Some(atomic_params(0.7, 1.0, 1.0)),
        ];

        let ff = construct_force_field_with_params(&mol, &params, 0.0, -1, false)
            .expect("trigonal bipyramidal special case should build");

        assert_eq!(ff.contribs().len(), 15);
    }

    #[test]
    fn uff_builder_construct_force_field_rejects_context_length_mismatch() {
        let mol = disconnected_pair_with_3d([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);

        let err = match construct_force_field(
            &mol,
            &[4],
            &[Hybridization::Sp3, Hybridization::Sp3],
            &[false],
            1.0,
            -1,
            false,
        ) {
            Ok(_) => panic!("context mismatch should fail"),
            Err(err) => err,
        };

        assert!(matches!(
            err,
            UffBuilderError::AtomTyper(UffAtomTyperError::ContextLengthMismatch {
                atoms: 2,
                total_valences: 1,
                hybridizations: 2,
                conjugation_flags: 1
            })
        ));
    }

    #[test]
    fn uff_builder_construct_force_field_rejects_invalid_conf_id_below_minus_one() {
        let mol = disconnected_pair_with_3d([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);

        let err = match construct_force_field(
            &mol,
            &[4, 4],
            &[Hybridization::Sp3, Hybridization::Sp3],
            &[false, false],
            1.0,
            -2,
            false,
        ) {
            Ok(_) => panic!("conf_id < -1 should fail"),
            Err(err) => err,
        };

        assert_eq!(err, UffBuilderError::InvalidConformerId { conf_id: -2 });
    }

    #[test]
    fn uff_builder_construct_force_field_rejects_missing_3d_conformer() {
        let mol = disconnected_pair();

        let err = match construct_force_field(
            &mol,
            &[4, 4],
            &[Hybridization::Sp3, Hybridization::Sp3],
            &[false, false],
            1.0,
            -1,
            false,
        ) {
            Ok(_) => panic!("missing 3d conformer should fail"),
            Err(err) => err,
        };

        assert_eq!(err, UffBuilderError::Missing3dConformer { conf_id: -1 });
    }

    #[test]
    fn uff_builder_construct_force_field_uses_default_conformer_for_minus_one() {
        let mol = single_bond_molecule_with_3d(BondOrder::Single, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);

        let ff = construct_force_field(
            &mol,
            &[4, 4],
            &[Hybridization::Sp3, Hybridization::Sp3],
            &[false, false],
            0.0,
            -1,
            false,
        )
        .expect("default conformer should build force field");

        assert_eq!(
            ff.positions(),
            &[ForceFieldVec3::new(1.0, 2.0, 3.0), ForceFieldVec3::new(4.0, 5.0, 6.0)]
        );
        assert_eq!(ff.contribs().len(), 1);
    }

    #[test]
    fn uff_builder_construct_force_field_selects_requested_conformer_by_id() {
        let mol = molecule_with_named_conformers(BondOrder::Single);

        let ff = construct_force_field(
            &mol,
            &[4, 4],
            &[Hybridization::Sp3, Hybridization::Sp3],
            &[false, false],
            0.0,
            9,
            false,
        )
        .expect("named conformer id should build force field");

        assert_eq!(
            ff.positions(),
            &[
                ForceFieldVec3::new(10.0, 11.0, 12.0),
                ForceFieldVec3::new(13.0, 14.0, 15.0)
            ]
        );
    }

    #[test]
    fn uff_builder_construct_force_field_respects_vdw_threshold_for_nonbonded_terms() {
        let mol = disconnected_pair_with_3d([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);

        let ff_keep = construct_force_field(
            &mol,
            &[4, 4],
            &[Hybridization::Sp3, Hybridization::Sp3],
            &[false, false],
            1.0,
            -1,
            false,
        )
        .expect("nonbonded pair should be included");
        let ff_skip = construct_force_field(
            &mol,
            &[4, 4],
            &[Hybridization::Sp3, Hybridization::Sp3],
            &[false, false],
            0.0,
            -1,
            false,
        )
        .expect("zero threshold should skip nonbonded pair");

        assert_eq!(ff_keep.contribs().len(), 1);
        assert_eq!(ff_skip.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_construct_force_field_ignores_interfragment_interactions_when_requested() {
        let mol = disconnected_pair_with_3d([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);

        let ff_keep = construct_force_field(
            &mol,
            &[4, 4],
            &[Hybridization::Sp3, Hybridization::Sp3],
            &[false, false],
            1.0,
            -1,
            false,
        )
        .expect("interfragment interactions should be kept");
        let ff_skip = construct_force_field(
            &mol,
            &[4, 4],
            &[Hybridization::Sp3, Hybridization::Sp3],
            &[false, false],
            1.0,
            -1,
            true,
        )
        .expect("interfragment interactions should be skipped");

        assert_eq!(ff_keep.contribs().len(), 1);
        assert_eq!(ff_skip.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_construct_force_field_skips_missing_parameter_slots() {
        let mol = pair_with_elements_and_3d(Element::C, Element::F, [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);

        let ff = construct_force_field(
            &mol,
            &[0, 0],
            &[Hybridization::S, Hybridization::Sp3],
            &[false, false],
            1.0,
            -1,
            false,
        )
        .expect("missing parameter slots should be skipped");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn uff_builder_construct_force_field_adds_trigonal_bipyramidal_special_case_angles() {
        let mol = trigonal_bipyramidal_molecule_with_3d();

        let ff = construct_force_field(
            &mol,
            &[5, 4, 4, 4, 4, 4],
            &[
                Hybridization::Sp3d,
                Hybridization::Sp3,
                Hybridization::Sp3,
                Hybridization::Sp3,
                Hybridization::Sp3,
                Hybridization::Sp3,
            ],
            &[false, false, false, false, false, false],
            0.0,
            -1,
            false,
        )
        .expect("trigonal bipyramidal special case should build");

        assert_eq!(ff.contribs().len(), 15);
    }
}
