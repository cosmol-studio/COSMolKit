//! Source-backed RDKit UFF atom typing helpers.

use super::builder::{UffBuilderError, construct_force_field, select_uff_conformer_index};
use super::params::{
    AtomicParams, ParamCollection, RAD2DEG, UffAngle, UffBond, UffInv, UffParamError, UffTor,
    UffVdw,
};
use super::torsion_angle::calc_torsion_params;
use super::utils::{
    UffUtilsError, calc_angle_force_constant, calc_bond_force_constant, calc_bond_rest_length,
    calc_inversion_coefficients, calc_nonbonded_depth, calc_nonbonded_minimum,
};
use crate::chemistry::valence::{
    ValenceError, ValenceModel, assign_valence, bond_type_as_double,
    periodic_table_outer_electrons, rdkit_default_valence, rdkit_element_symbol,
};
use crate::rings::RingFindingError;
use crate::{Atom, Hybridization, Molecule};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub(crate) enum UffAtomTyperError {
    #[error(transparent)]
    Valence(#[from] ValenceError),
    #[error(transparent)]
    Params(#[from] UffParamError),
    #[error(
        "UFF atom typing context length mismatch: atoms={atoms}, total_valences={total_valences}, hybridizations={hybridizations}, conjugation_flags={conjugation_flags}"
    )]
    ContextLengthMismatch {
        atoms: usize,
        total_valences: usize,
        hybridizations: usize,
        conjugation_flags: usize,
    },
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum UffPublicApiError {
    #[error("UFF atom index out of range: atom_index={atom_index}, atoms={atoms}")]
    InvalidAtomIndex { atom_index: usize, atoms: usize },
    #[error("RDKit UFF buildNeighborMatrix is unsupported for empty molecules")]
    EmptyMolecule,
    #[error("UFF constructForceField requires a 3D conformer when conf_id={conf_id}")]
    Missing3dConformer { conf_id: isize },
    #[error("UFF constructForceField conf_id must be >= -1, got {conf_id}")]
    InvalidConformerId { conf_id: isize },
    #[error(
        "UFF addTrigonalBipyramidAngles requires central atom parameters for atom {atom_index}"
    )]
    MissingTrigonalBipyramidCenterParams { atom_index: usize },
    #[error("UFF parameter length mismatch: atoms={atoms}, params={params}")]
    ParamsLengthMismatch { atoms: usize, params: usize },
    #[error("unsupported torsion bond SMARTS for UFF addTorsions: {0}")]
    UnsupportedTorsionBondSmarts(String),
    #[error(
        "UFF atom typing context length mismatch: atoms={atoms}, total_valences={total_valences}, hybridizations={hybridizations}, conjugation_flags={conjugation_flags}"
    )]
    AtomTypingContextLengthMismatch {
        atoms: usize,
        total_valences: usize,
        hybridizations: usize,
        conjugation_flags: usize,
    },
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    Params(#[from] UffParamError),
    #[error(transparent)]
    UffUtils(#[from] UffUtilsError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
}

#[derive(Debug, Clone, PartialEq)]
pub struct UffOptimizeMoleculeResult {
    pub molecule: Molecule,
    pub needs_more: i32,
    pub energy: f64,
}

impl UffOptimizeMoleculeResult {
    #[must_use]
    pub fn needs_more(&self) -> bool {
        self.needs_more > 0
    }

    #[must_use]
    pub fn status_code(&self) -> i32 {
        self.needs_more
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UffOptimizeMoleculeConfResult {
    pub needs_more: i32,
    pub energy: f64,
}

impl UffOptimizeMoleculeConfResult {
    #[must_use]
    pub fn needs_more(&self) -> bool {
        self.needs_more > 0
    }

    #[must_use]
    pub fn status_code(&self) -> i32 {
        self.needs_more
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct UffOptimizeMoleculeConfsResult {
    pub molecule: Molecule,
    pub conformer_results: Vec<UffOptimizeMoleculeConfResult>,
}

struct UffTypingContext {
    total_valences: Vec<i32>,
    hybridizations: Vec<Hybridization>,
    atom_has_conjugated_bond: Vec<bool>,
}

fn public_error_from_atom_typer(err: UffAtomTyperError) -> UffPublicApiError {
    match err {
        UffAtomTyperError::Valence(err) => UffPublicApiError::Valence(err),
        UffAtomTyperError::Params(err) => UffPublicApiError::Params(err),
        UffAtomTyperError::ContextLengthMismatch {
            atoms,
            total_valences,
            hybridizations,
            conjugation_flags,
        } => UffPublicApiError::AtomTypingContextLengthMismatch {
            atoms,
            total_valences,
            hybridizations,
            conjugation_flags,
        },
    }
}

fn public_error_from_builder(err: UffBuilderError) -> UffPublicApiError {
    match err {
        UffBuilderError::EmptyMolecule => UffPublicApiError::EmptyMolecule,
        UffBuilderError::Missing3dConformer { conf_id } => {
            UffPublicApiError::Missing3dConformer { conf_id }
        }
        UffBuilderError::InvalidConformerId { conf_id } => {
            UffPublicApiError::InvalidConformerId { conf_id }
        }
        UffBuilderError::MissingTrigonalBipyramidCenterParams { atom_index } => {
            UffPublicApiError::MissingTrigonalBipyramidCenterParams { atom_index }
        }
        UffBuilderError::ParamsLengthMismatch { atoms, params } => {
            UffPublicApiError::ParamsLengthMismatch { atoms, params }
        }
        UffBuilderError::RingFinding(err) => UffPublicApiError::RingFinding(err),
        UffBuilderError::UnsupportedTorsionBondSmarts(smarts) => {
            UffPublicApiError::UnsupportedTorsionBondSmarts(smarts)
        }
        UffBuilderError::AtomTyper(err) => public_error_from_atom_typer(err),
        UffBuilderError::Valence(err) => UffPublicApiError::Valence(err),
        UffBuilderError::UffUtils(err) => UffPublicApiError::UffUtils(err),
    }
}

fn uff_typing_context(mol: &Molecule) -> Result<UffTypingContext, UffPublicApiError> {
    let adjacency = &mol.topology_block().adjacency;
    let assignment = assign_valence(mol, ValenceModel::RdkitLike)?;
    let hybridizations: Vec<_> = mol.atoms().iter().map(Atom::hybridization).collect();
    let atom_has_conjugated_bond: Vec<_> = mol
        .atoms()
        .iter()
        .enumerate()
        .map(|(idx, _atom)| {
            adjacency
                .neighbors_of(idx)
                .iter()
                .any(|neighbor| mol.bonds()[neighbor.bond.index()].is_conjugated())
        })
        .collect();
    let total_valences: Vec<_> = assignment
        .explicit_valence
        .iter()
        .zip(assignment.implicit_hydrogens.iter())
        .map(|(explicit, implicit)| explicit + implicit)
        .collect();
    Ok(UffTypingContext {
        total_valences,
        hybridizations,
        atom_has_conjugated_bond,
    })
}

pub fn get_uff_bond_stretch_params(
    mol: &Molecule,
    idx1: usize,
    idx2: usize,
) -> Result<Option<UffBond>, UffPublicApiError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFBondStretchParams (AtomTyper.cpp:535-557)
    // RDKit✔️❗: bool getUFFBondStretchParams(const ROMol &mol, unsigned int idx1,
    // RDKit✔️❗:                              unsigned int idx2, UFFBond &uffBondStretchParams) {
    // RDKit✔️❗:   auto params = ParamCollection::getParams();
    let params = ParamCollection::get_params("")?;
    // RDKit✔️❗:   unsigned int idx[2] = {idx1, idx2};
    let idx = [idx1, idx2];
    // RDKit✔️❗:   AtomicParamVect paramVect(2);
    let mut param_vect = [None; 2];
    // RDKit✔️❗:   unsigned int i;
    // RDKit✔️❗:   const Bond *bond = mol.getBondBetweenAtoms(idx1, idx2);
    let atom_count = mol.atoms().len();
    for &atom_index in &idx {
        if atom_index >= atom_count {
            return Err(UffPublicApiError::InvalidAtomIndex {
                atom_index,
                atoms: atom_count,
            });
        }
    }
    let adjacency = &mol.topology_block().adjacency;
    let bond = adjacency.neighbors_of(idx1).iter().find_map(|neighbor| {
        (neighbor.atom_index == idx2).then_some(&mol.bonds()[neighbor.bond.index()])
    });
    // RDKit✔️❗:   bool res = bond != nullptr;
    let Some(bond) = bond else {
        return Ok(None);
    };
    // RDKit✔️❗:   for (i = 0; res && (i < 2); ++i) {
    let assignment = assign_valence(mol, ValenceModel::RdkitLike)?;
    for i in 0..2 {
        // RDKit✔️❗:     const Atom *atom = mol.getAtomWithIdx(idx[i]);
        let atom = &mol.atoms()[idx[i]];
        // RDKit✔️❗:     std::string atomKey = Tools::getAtomLabel(atom);
        let atom_has_conjugated_bond = adjacency
            .neighbors_of(idx[i])
            .iter()
            .any(|neighbor| mol.bonds()[neighbor.bond.index()].is_conjugated());
        let atom_key = get_atom_label_for_uff(
            atom,
            idx[i],
            assignment.explicit_valence[idx[i]] + assignment.implicit_hydrogens[idx[i]],
            atom.hybridization(),
            atom_has_conjugated_bond,
            true,
        )?;
        // RDKit✔️❗:     paramVect[i] = (*params)(atomKey);
        param_vect[i] = params.get(&atom_key).copied();
        // RDKit✔️❗:     res = paramVect[i] != nullptr;
        if param_vect[i].is_none() {
            return Ok(None);
        }
        // RDKit✔️❗:   }
    }
    // RDKit✔️❗:   if (res) {
    let param0 = param_vect[0].expect("param existence checked");
    let param1 = param_vect[1].expect("param existence checked");
    // RDKit✔️❗:     double bondOrder = bond->getBondTypeAsDouble();
    let bond_order = bond_type_as_double(bond.order())?;
    // RDKit✔️❗:     uffBondStretchParams.r0 =
    // RDKit✔️❗:         UFF::Utils::calcBondRestLength(bondOrder, paramVect[0], paramVect[1]);
    let r0 = calc_bond_rest_length(bond_order, &param0, &param1)?;
    // RDKit✔️❗:     uffBondStretchParams.kb = UFF::Utils::calcBondForceConstant(
    // RDKit✔️❗:         uffBondStretchParams.r0, paramVect[0], paramVect[1]);
    let kb = calc_bond_force_constant(r0, &param0, &param1)?;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return res;
    Ok(Some(UffBond { kb, r0 }))
    // RDKit✔️❗: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFBondStretchParams
}

pub fn get_uff_angle_bend_params(
    mol: &Molecule,
    idx1: usize,
    idx2: usize,
    idx3: usize,
) -> Result<Option<UffAngle>, UffPublicApiError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFAngleBendParams (AtomTyper.cpp:559-589)
    // RDKit✔️❗: bool getUFFAngleBendParams(const ROMol &mol, unsigned int idx1,
    // RDKit✔️❗:                            unsigned int idx2, unsigned int idx3,
    // RDKit✔️❗:                            UFFAngle &uffAngleBendParams) {
    // RDKit✔️❗:   auto params = ParamCollection::getParams();
    let params = ParamCollection::get_params("")?;
    // RDKit✔️❗:   unsigned int idx[3] = {idx1, idx2, idx3};
    let idx = [idx1, idx2, idx3];
    // RDKit✔️❗:   AtomicParamVect paramVect(3);
    let mut param_vect = [None; 3];
    // RDKit✔️❗:   unsigned int i;
    // RDKit✔️❗:   const Bond *bond[2];
    let mut bonds = [None, None];
    // RDKit✔️❗:   bool res = true;
    let atom_count = mol.atoms().len();
    for &atom_index in &idx {
        if atom_index >= atom_count {
            return Err(UffPublicApiError::InvalidAtomIndex {
                atom_index,
                atoms: atom_count,
            });
        }
    }
    let adjacency = &mol.topology_block().adjacency;
    let assignment = assign_valence(mol, ValenceModel::RdkitLike)?;
    // RDKit✔️❗:   for (i = 0; res && (i < 3); ++i) {
    for i in 0..3 {
        // RDKit✔️❗:     if (i < 2) {
        if i < 2 {
            // RDKit✔️❗:       bond[i] = mol.getBondBetweenAtoms(idx[i], idx[i + 1]);
            bonds[i] = adjacency.neighbors_of(idx[i]).iter().find_map(|neighbor| {
                (neighbor.atom_index == idx[i + 1]).then_some(&mol.bonds()[neighbor.bond.index()])
            });
            // RDKit✔️❗:       res = bond[i] != nullptr;
            if bonds[i].is_none() {
                return Ok(None);
            }
            // RDKit✔️❗:     }
        }
        // RDKit✔️❗:     if (res) {
        // RDKit✔️❗:       const Atom *atom = mol.getAtomWithIdx(idx[i]);
        let atom = &mol.atoms()[idx[i]];
        // RDKit✔️❗:       std::string atomKey = Tools::getAtomLabel(atom);
        let atom_has_conjugated_bond = adjacency
            .neighbors_of(idx[i])
            .iter()
            .any(|neighbor| mol.bonds()[neighbor.bond.index()].is_conjugated());
        let atom_key = get_atom_label_for_uff(
            atom,
            idx[i],
            assignment.explicit_valence[idx[i]] + assignment.implicit_hydrogens[idx[i]],
            atom.hybridization(),
            atom_has_conjugated_bond,
            true,
        )?;
        // RDKit✔️❗:       paramVect[i] = (*params)(atomKey);
        param_vect[i] = params.get(&atom_key).copied();
        // RDKit✔️❗:       res = paramVect[i] != nullptr;
        if param_vect[i].is_none() {
            return Ok(None);
        }
        // RDKit✔️❗:     }
        // RDKit✔️❗:   }
    }
    let bond0 = bonds[0].expect("bond existence checked");
    let bond1 = bonds[1].expect("bond existence checked");
    let param0 = param_vect[0].expect("param existence checked");
    let param1 = param_vect[1].expect("param existence checked");
    let param2 = param_vect[2].expect("param existence checked");
    // RDKit✔️❗:   if (res) {
    // RDKit✔️❗:     double bondOrder12 = bond[0]->getBondTypeAsDouble();
    let bond_order12 = bond_type_as_double(bond0.order())?;
    // RDKit✔️❗:     double bondOrder23 = bond[1]->getBondTypeAsDouble();
    let bond_order23 = bond_type_as_double(bond1.order())?;
    // RDKit✔️❗:     uffAngleBendParams.theta0 = RAD2DEG * paramVect[1]->theta0;
    let theta0 = RAD2DEG * param1.theta0;
    // RDKit✔️❗:     uffAngleBendParams.ka = UFF::Utils::calcAngleForceConstant(
    // RDKit✔️❗:         paramVect[1]->theta0, bondOrder12, bondOrder23, paramVect[0],
    // RDKit✔️❗:         paramVect[1], paramVect[2]);
    let ka = calc_angle_force_constant(
        param1.theta0,
        bond_order12,
        bond_order23,
        &param0,
        &param1,
        &param2,
    )?;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return res;
    Ok(Some(UffAngle { ka, theta0 }))
    // RDKit✔️❗: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFAngleBendParams
}

pub fn get_uff_torsion_params(
    mol: &Molecule,
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
) -> Result<Option<UffTor>, UffPublicApiError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFTorsionParams (AtomTyper.cpp:591-665)
    // RDKit✔️❗: bool getUFFTorsionParams(const ROMol &mol, unsigned int idx1, unsigned int idx2,
    // RDKit✔️❗:                          unsigned int idx3, unsigned int idx4,
    // RDKit✔️❗:                          UFFTor &uffTorsionParams) {
    // RDKit✔️❗:   auto params = ParamCollection::getParams();
    let params = ParamCollection::get_params("")?;
    // RDKit✔️❗:   unsigned int idx[4] = {idx1, idx2, idx3, idx4};
    let idx = [idx1, idx2, idx3, idx4];
    // RDKit✔️❗:   AtomicParamVect paramVect(2);
    let mut param_vect = [None, None];
    // RDKit✔️❗:   unsigned int i;
    // RDKit✔️❗:   const Bond *bond = mol.getBondBetweenAtoms(idx2, idx3);
    let atom_count = mol.atoms().len();
    for &atom_index in &idx {
        if atom_index >= atom_count {
            return Err(UffPublicApiError::InvalidAtomIndex {
                atom_index,
                atoms: atom_count,
            });
        }
    }
    let adjacency = &mol.topology_block().adjacency;
    let center_bond = adjacency.neighbors_of(idx2).iter().find_map(|neighbor| {
        (neighbor.atom_index == idx3).then_some(&mol.bonds()[neighbor.bond.index()])
    });
    let Some(center_bond) = center_bond else {
        return Ok(None);
    };
    // RDKit✔️❗:   int atNum[2];
    let mut at_num = [0_i32; 2];
    // RDKit✔️❗:   Atom::HybridizationType hyb[2];
    let mut hyb = [Hybridization::Unspecified; 2];
    // RDKit✔️❗:   bool res = true;
    // RDKit✔️❗:   bool hasSP2 = false;
    let mut has_sp2 = false;
    let assignment = assign_valence(mol, ValenceModel::RdkitLike)?;
    // RDKit✔️❗:   for (i = 0; res && (i < 4); ++i) {
    for i in 0..4 {
        // RDKit✔️❗:     if (i < 3) {
        if i < 3 {
            // RDKit✔️❗:       res = mol.getBondBetweenAtoms(idx[i], idx[i + 1]) != nullptr;
            let has_bond = adjacency
                .neighbors_of(idx[i])
                .iter()
                .any(|neighbor| neighbor.atom_index == idx[i + 1]);
            if !has_bond {
                return Ok(None);
            }
            // RDKit✔️❗:     }
        }
        // RDKit✔️❗:     const Atom *atom = mol.getAtomWithIdx(idx[i]);
        let atom = &mol.atoms()[idx[i]];
        // RDKit✔️❗:     if ((i == 1) || (i == 2)) {
        if i == 1 || i == 2 {
            // RDKit✔️❗:       unsigned int j = i - 1;
            let j = i - 1;
            // RDKit✔️❗:       atNum[j] = atom->getAtomicNum();
            at_num[j] = i32::from(atom.atomic_number());
            // RDKit✔️❗:       hyb[j] = atom->getHybridization();
            hyb[j] = atom.hybridization();
            // RDKit✔️❗:       std::string atomKey = Tools::getAtomLabel(atom);
            let atom_has_conjugated_bond = adjacency
                .neighbors_of(idx[i])
                .iter()
                .any(|neighbor| mol.bonds()[neighbor.bond.index()].is_conjugated());
            let atom_key = get_atom_label_for_uff(
                atom,
                idx[i],
                assignment.explicit_valence[idx[i]] + assignment.implicit_hydrogens[idx[i]],
                atom.hybridization(),
                atom_has_conjugated_bond,
                true,
            )?;
            // RDKit✔️❗:       paramVect[j] = (*params)(atomKey);
            param_vect[j] = params.get(&atom_key).copied();
            // RDKit✔️❗:       res = paramVect[j] != nullptr;
            if param_vect[j].is_none() {
                return Ok(None);
            }
            // RDKit✔️❗:     } else if (atom->getHybridization() == Atom::SP2) {
        } else if atom.hybridization() == Hybridization::Sp2 {
            // RDKit✔️❗:       hasSP2 = true;
            has_sp2 = true;
            // RDKit✔️❗:     }
        }
        // RDKit✔️❗:   }
    }
    // RDKit✔️❗:   if (res) {
    // RDKit✔️❗:     res = (((hyb[0] == RDKit::Atom::SP2) || (hyb[0] == RDKit::Atom::SP3)) &&
    // RDKit✔️❗:            ((hyb[1] == RDKit::Atom::SP2) || (hyb[1] == RDKit::Atom::SP3)));
    if !matches!(hyb[0], Hybridization::Sp2 | Hybridization::Sp3)
        || !matches!(hyb[1], Hybridization::Sp2 | Hybridization::Sp3)
    {
        return Ok(None);
    }
    let param0 = param_vect[0].expect("param existence checked");
    let param1 = param_vect[1].expect("param existence checked");
    // RDKit✔️❗:   if (res) {
    // RDKit✔️❗:     double bondOrder = bond->getBondTypeAsDouble();
    let bond_order = bond_type_as_double(center_bond.order())?;
    // RDKit✔️❗:     if ((hyb[0] == RDKit::Atom::SP3) && (hyb[1] == RDKit::Atom::SP3)) {
    // RDKit✔️❗:       // general case:
    // RDKit✔️❗:       uffTorsionParams.V = sqrt(paramVect[0]->V1 * paramVect[1]->V1);
    // RDKit✔️❗:       // special case for single bonds between group 6 elements:
    // RDKit✔️❗:       if (((int)(bondOrder * 10) == 10) && UFF::Utils::isInGroup6(atNum[0]) &&
    // RDKit✔️❗:           UFF::Utils::isInGroup6(atNum[1])) {
    // RDKit✔️❗:         ...
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if ((hyb[0] == RDKit::Atom::SP2) && (hyb[1] == RDKit::Atom::SP2)) {
    // RDKit✔️❗:       uffTorsionParams.V =
    // RDKit✔️❗:           UFF::Utils::equation17(bondOrder, paramVect[0], paramVect[1]);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       // SP2 - SP3 ...
    // RDKit✔️❗:     }
    let (v, _order, _cos_term) = calc_torsion_params(
        bond_order, at_num[0], at_num[1], hyb[0], hyb[1], &param0, &param1, has_sp2,
    );
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return res;
    Ok(Some(UffTor { v }))
    // RDKit✔️❗: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFTorsionParams
}

pub fn get_uff_inversion_params(
    mol: &Molecule,
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
) -> Result<Option<UffInv>, UffPublicApiError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFInversionParams (AtomTyper.cpp:667-711)
    // RDKit✔️❗: bool getUFFInversionParams(const ROMol &mol, unsigned int idx1,
    // RDKit✔️❗:                            unsigned int idx2, unsigned int idx3,
    // RDKit✔️❗:                            unsigned int idx4, UFFInv &uffInversionParams) {
    // RDKit✔️❗:   unsigned int idx[4] = {idx1, idx2, idx3, idx4};
    let idx = [idx1, idx2, idx3, idx4];
    let atom_count = mol.atoms().len();
    for &atom_index in &idx {
        if atom_index >= atom_count {
            return Err(UffPublicApiError::InvalidAtomIndex {
                atom_index,
                atoms: atom_count,
            });
        }
    }
    let adjacency = &mol.topology_block().adjacency;
    let has_bond = |a: usize, b: usize| {
        adjacency
            .neighbors_of(a)
            .iter()
            .any(|neighbor| neighbor.atom_index == b)
    };
    // RDKit✔️❗:   bool res = (mol.getBondBetweenAtoms(idx1, idx2) &&
    // RDKit✔️❗:               mol.getBondBetweenAtoms(idx2, idx3) &&
    // RDKit✔️❗:               mol.getBondBetweenAtoms(idx2, idx4));
    if !(has_bond(idx1, idx2) && has_bond(idx2, idx3) && has_bond(idx2, idx4)) {
        return Ok(None);
    }
    // RDKit✔️❗:   unsigned int i;
    // RDKit✔️❗:   // bool isAtom2C = false;
    // RDKit✔️❗:   bool isBoundToSP2O = false;
    let mut is_bound_to_sp2_o = false;
    // RDKit✔️❗:   unsigned int at2AtomicNum = 0;
    let mut at2_atomic_num = 0_u8;
    // RDKit✔️❗:   for (i = 0; res && (i < 4); ++i) {
    for i in 0..4 {
        // RDKit✔️❗:     const Atom *atom = mol.getAtomWithIdx(idx[i]);
        let atom = &mol.atoms()[idx[i]];
        // RDKit✔️❗:     if (i == 1) {
        if i == 1 {
            // RDKit✔️❗:       at2AtomicNum = atom->getAtomicNum();
            at2_atomic_num = atom.atomic_number();
            // RDKit✔️❗:       if (res) {
            // RDKit✔️❗:         // if the central atom is not carbon, nitrogen, oxygen,
            // RDKit✔️❗:         // phosphorous, arsenic, antimonium or bismuth, skip it
            // RDKit✔️❗:         res = (!(((at2AtomicNum != 6) && (at2AtomicNum != 7) &&
            // RDKit✔️❗:                   (at2AtomicNum != 8) && (at2AtomicNum != 15) &&
            // RDKit✔️❗:                   (at2AtomicNum != 33) && (at2AtomicNum != 51) &&
            // RDKit✔️❗:                   (at2AtomicNum != 83)) ||
            // RDKit✔️❗:                  (atom->getDegree() != 3)));
            if !matches!(at2_atomic_num, 6 | 7 | 8 | 15 | 33 | 51 | 83)
                || adjacency.neighbors_of(idx[i]).len() != 3
            {
                return Ok(None);
            }
            // RDKit✔️❗:       }
            // RDKit✔️❗:       if (res) {
            // RDKit✔️❗:         // if the central atom is carbon, nitrogen or oxygen
            // RDKit✔️❗:         // but hybridization is not sp2, skip it
            // RDKit✔️❗:         res = (!(((at2AtomicNum == 6) || (at2AtomicNum == 7) ||
            // RDKit✔️❗:                   (at2AtomicNum == 8)) &&
            // RDKit✔️❗:                  (atom->getHybridization() != Atom::SP2)));
            if matches!(at2_atomic_num, 6 | 7 | 8) && atom.hybridization() != Hybridization::Sp2 {
                return Ok(None);
            }
            // RDKit✔️❗:       }
            // RDKit✔️❗:     } else if ((atom->getAtomicNum() == 8) &&
            // RDKit✔️❗:                (atom->getHybridization() == Atom::SP2)) {
        } else if atom.atomic_number() == 8 && atom.hybridization() == Hybridization::Sp2 {
            // RDKit✔️❗:       isBoundToSP2O = true;
            is_bound_to_sp2_o = true;
            // RDKit✔️❗:     }
        }
        // RDKit✔️❗:   }
    }
    // RDKit✔️❗:   if (res) {
    // RDKit✔️❗:     isBoundToSP2O = (isBoundToSP2O && (at2AtomicNum == 6));
    is_bound_to_sp2_o = is_bound_to_sp2_o && at2_atomic_num == 6;
    // RDKit✔️❗:     auto invCoeffForceCon =
    // RDKit✔️❗:         UFF::Utils::calcInversionCoefficientsAndForceConstant(at2AtomicNum,
    // RDKit✔️❗:                                                               isBoundToSP2O);
    let (k, _c0, _c1, _c2) =
        calc_inversion_coefficients(i32::from(at2_atomic_num), is_bound_to_sp2_o);
    // RDKit✔️❗:     uffInversionParams.K = std::get<0>(invCoeffForceCon);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return res;
    Ok(Some(UffInv { k }))
    // RDKit✔️❗: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFInversionParams
}

pub fn get_uff_vdw_params(
    mol: &Molecule,
    idx1: usize,
    idx2: usize,
) -> Result<Option<UffVdw>, UffPublicApiError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFVdWParams (AtomTyper.cpp:713-733)
    // RDKit✔️❗: bool getUFFVdWParams(const ROMol &mol, unsigned int idx1, unsigned int idx2,
    // RDKit✔️❗:                      UFFVdW &uffVdWParams) {
    // RDKit✔️❗:   bool res = true;
    // RDKit✔️❗:   auto params = ParamCollection::getParams();
    let params = ParamCollection::get_params("")?;
    // RDKit✔️❗:   unsigned int idx[2] = {idx1, idx2};
    let idx = [idx1, idx2];
    // RDKit✔️❗:   AtomicParamVect paramVect(2);
    let mut param_vect = [None; 2];
    // RDKit✔️❗:   unsigned int i;
    let atom_count = mol.atoms().len();
    for &atom_index in &idx {
        if atom_index >= atom_count {
            return Err(UffPublicApiError::InvalidAtomIndex {
                atom_index,
                atoms: atom_count,
            });
        }
    }
    let adjacency = &mol.topology_block().adjacency;
    let assignment = assign_valence(mol, ValenceModel::RdkitLike)?;
    // RDKit✔️❗:   for (i = 0; res && (i < 2); ++i) {
    for i in 0..2 {
        // RDKit✔️❗:     const Atom *atom = mol.getAtomWithIdx(idx[i]);
        let atom = &mol.atoms()[idx[i]];
        // RDKit✔️❗:     std::string atomKey = Tools::getAtomLabel(atom);
        let atom_has_conjugated_bond = adjacency
            .neighbors_of(idx[i])
            .iter()
            .any(|neighbor| mol.bonds()[neighbor.bond.index()].is_conjugated());
        let atom_key = get_atom_label_for_uff(
            atom,
            idx[i],
            assignment.explicit_valence[idx[i]] + assignment.implicit_hydrogens[idx[i]],
            atom.hybridization(),
            atom_has_conjugated_bond,
            true,
        )?;
        // RDKit✔️❗:     paramVect[i] = (*params)(atomKey);
        param_vect[i] = params.get(&atom_key).copied();
        // RDKit✔️❗:     res = paramVect[i] != nullptr;
        if param_vect[i].is_none() {
            return Ok(None);
        }
        // RDKit✔️❗:   }
    }
    // RDKit✔️❗:   if (res) {
    let param0 = param_vect[0].expect("param existence checked");
    let param1 = param_vect[1].expect("param existence checked");
    // RDKit✔️❗:     uffVdWParams.x_ij =
    // RDKit✔️❗:         UFF::Utils::calcNonbondedMinimum(paramVect[0], paramVect[1]);
    let x_ij = calc_nonbonded_minimum(&param0, &param1);
    // RDKit✔️❗:     uffVdWParams.D_ij =
    // RDKit✔️❗:         UFF::Utils::calcNonbondedDepth(paramVect[0], paramVect[1]);
    let d_ij = calc_nonbonded_depth(&param0, &param1);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return res;
    Ok(Some(UffVdw { x_ij, d_ij }))
    // RDKit✔️❗: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFVdWParams
}

pub fn uff_has_all_molecule_params(mol: &Molecule) -> Result<bool, UffPublicApiError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFFHasAllMoleculeParams (rdForceFields.cpp:107-112)
    // RDKit✔️❗: bool UFFHasAllMoleculeParams(const ROMol &mol) {
    // RDKit✔️❗:   UFF::AtomicParamVect types;
    // RDKit✔️❗:   bool foundAll;
    // RDKit✔️❗:   boost::tie(types, foundAll) = UFF::getAtomTypes(mol);
    let context = uff_typing_context(mol)?;
    let (_types, found_all) = get_atom_types_for_uff(
        mol,
        &context.total_valences,
        &context.hybridizations,
        &context.atom_has_conjugated_bond,
    )
    .map_err(public_error_from_atom_typer)?;
    // RDKit✔️❗:   return foundAll;
    Ok(found_all)
    // RDKit✔️❗: }
    // END RDKIT CPP FUNCTION RDKit::UFFHasAllMoleculeParams
}

#[doc(hidden)]
pub fn uff_initial_gradient_for_parity(
    mol: &Molecule,
    vdw_thresh: f64,
    conf_id: isize,
    ignore_interfrag_interactions: bool,
) -> Result<Vec<f64>, UffPublicApiError> {
    let context = uff_typing_context(mol)?;
    let mut ff = construct_force_field(
        mol,
        &context.total_valences,
        &context.hybridizations,
        &context.atom_has_conjugated_bond,
        vdw_thresh,
        conf_id,
        ignore_interfrag_interactions,
    )
    .map_err(public_error_from_builder)?;
    ff.initialize();
    let mut grad = vec![0.0; ff.dimension() * ff.positions().len()];
    ff.calc_grad_current(&mut grad);
    Ok(grad)
}

pub fn uff_optimize_molecule(
    mol: &Molecule,
    max_iters: usize,
    vdw_thresh: f64,
    conf_id: isize,
    ignore_interfrag_interactions: bool,
) -> Result<UffOptimizeMoleculeResult, UffPublicApiError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::UFFOptimizeMolecule (UFF.h:40-49)
    // RDKit❗✔️: inline std::pair<int, double> UFFOptimizeMolecule(
    // RDKit❗✔️:     ROMol &mol, int maxIters = 1000, double vdwThresh = 10.0, int confId = -1,
    // RDKit❗✔️:     bool ignoreInterfragInteractions = true) {
    // RDKit❗✔️:   std::unique_ptr<ForceFields::ForceField> ff(UFF::constructForceField(
    // RDKit❗✔️:       mol, vdwThresh, confId, ignoreInterfragInteractions));
    let context = uff_typing_context(mol)?;
    let conf_index = select_uff_conformer_index(mol, conf_id).map_err(public_error_from_builder)?;
    let mut ff = construct_force_field(
        mol,
        &context.total_valences,
        &context.hybridizations,
        &context.atom_has_conjugated_bond,
        vdw_thresh,
        conf_id,
        ignore_interfrag_interactions,
    )
    .map_err(public_error_from_builder)?;
    // RDKit❗✔️:   std::pair<int, double> res =
    // RDKit❗✔️:       ForceFieldsHelper::OptimizeMolecule(*ff, maxIters);
    //
    // BEGIN RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::OptimizeMolecule (FFConvenience.h:94-101)
    // RDKit✔️✔️: inline std::pair<int, double> OptimizeMolecule(ForceFields::ForceField &ff,
    // RDKit✔️✔️:                                                int maxIters = 1000) {
    // RDKit✔️✔️:   ff.initialize();
    ff.initialize();
    // RDKit✔️✔️:   int res = ff.minimize(maxIters);
    // RDKit default ForceField::minimize(maxIts) forwards forceTol=1e-4 and energyTol=1e-6.
    let needs_more = ff.minimize(max_iters, 1.0e-4, 1.0e-6);
    // RDKit✔️✔️:   double e = ff.calcEnergy();
    let energy = ff.calc_energy_current(None);
    // RDKit✔️✔️:   return std::make_pair(res, e);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::OptimizeMolecule
    // RDKit❗✔️:   return res;
    let mut molecule = mol.clone();
    {
        let coords = molecule.coordinate_block_mut().conformers_3d[conf_index].coordinates_mut();
        for (coord, position) in coords.iter_mut().zip(ff.positions()) {
            *coord = [position.x, position.y, position.z];
        }
    }
    Ok(UffOptimizeMoleculeResult {
        molecule,
        needs_more,
        energy,
    })
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::UFF::UFFOptimizeMolecule
}

pub fn uff_optimize_molecule_confs(
    mol: &Molecule,
    num_threads: i32,
    max_iters: usize,
    vdw_thresh: f64,
    ignore_interfrag_interactions: bool,
) -> Result<UffOptimizeMoleculeConfsResult, UffPublicApiError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::UFFOptimizeMoleculeConfs (UFF.h:67-76)
    // RDKit✔️❌: inline void UFFOptimizeMoleculeConfs(ROMol &mol,
    // RDKit✔️❌:                                      std::vector<std::pair<int, double>> &res,
    // RDKit✔️❌:                                      int numThreads = 1, int maxIters = 1000,
    // RDKit✔️❌:                                      double vdwThresh = 10.0,
    // RDKit✔️❌:                                      bool ignoreInterfragInteractions = true) {
    // RDKit✔️❌:   std::unique_ptr<ForceFields::ForceField> ff(UFF::constructForceField(
    // RDKit✔️❌:       mol, vdwThresh, -1, ignoreInterfragInteractions));
    let context = uff_typing_context(mol)?;
    let mut ff = construct_force_field(
        mol,
        &context.total_valences,
        &context.hybridizations,
        &context.atom_has_conjugated_bond,
        vdw_thresh,
        -1,
        ignore_interfrag_interactions,
    )
    .map_err(public_error_from_builder)?;
    // RDKit✔️❌:   ForceFieldsHelper::OptimizeMoleculeConfs(mol, *ff, res, numThreads, maxIters);
    //
    // BEGIN RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::OptimizeMoleculeConfs (FFConvenience.h:116-127)
    // RDKit✔️❌: inline void OptimizeMoleculeConfs(ROMol &mol, ForceFields::ForceField &ff,
    // RDKit✔️❌:                                   std::vector<std::pair<int, double>> &res,
    // RDKit✔️❌:                                   int numThreads = 1, int maxIters = 1000) {
    // RDKit✔️❌:   res.resize(mol.getNumConformers());
    let mut molecule = mol.clone();
    let conformer_count = molecule.conformers_3d().len();
    let mut conformer_results = vec![
        UffOptimizeMoleculeConfResult {
            needs_more: 0,
            energy: 0.0,
        };
        conformer_count
    ];
    // RDKit✔️❌:   numThreads = getNumThreadsToUse(numThreads);
    //
    // BEGIN RDKIT CPP FUNCTION RDKit::getNumThreadsToUse without RDK_BUILD_THREADSAFE_SSS (RDThreads.h:37-40)
    // RDKit✔️❌: inline unsigned int getNumThreadsToUse(int target) {
    // RDKit✔️❌:   RDUNUSED_PARAM(target);
    let _ = num_threads;
    // RDKit✔️❌:   return 1;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION RDKit::getNumThreadsToUse without RDK_BUILD_THREADSAFE_SSS
    // RDKit✔️❌:   if (numThreads == 1) {
    // RDKit✔️❌:     detail::OptimizeMoleculeConfsST(mol, ff, res, maxIters);
    //
    // BEGIN RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::detail::OptimizeMoleculeConfsST (FFConvenience.h:61-78)
    // RDKit✔️❌: inline void OptimizeMoleculeConfsST(ROMol &mol, ForceFields::ForceField &ff,
    // RDKit✔️❌:                                     std::vector<std::pair<int, double>> &res,
    // RDKit✔️❌:                                     int maxIters) {
    // RDKit✔️❌:   PRECONDITION(res.size() >= mol.getNumConformers(),
    // RDKit✔️❌:                "res.size() must be >= mol.getNumConformers()");
    assert!(conformer_results.len() >= conformer_count);
    // RDKit✔️❌:   unsigned int i = 0;
    // RDKit✔️❌:   for (ROMol::ConformerIterator cit = mol.beginConformers();
    // RDKit✔️❌:        cit != mol.endConformers(); ++cit, ++i) {
    for conf_index in 0..conformer_count {
        let start_coords = molecule.conformers_3d()[conf_index].coordinates().to_vec();
        // RDKit✔️❌:     for (unsigned int aidx = 0; aidx < mol.getNumAtoms(); ++aidx) {
        // RDKit✔️❌:       ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
        // RDKit stores pointers to the conformer rows; Rust owns ForceFieldVec3 values and
        // copies them back to the selected conformer after minimization.
        ff.positions_mut().clear();
        for coord in &start_coords {
            ff.positions_mut()
                .push(crate::chemistry::forcefield::core::ForceFieldVec3::new(
                    coord[0], coord[1], coord[2],
                ));
        }
        // RDKit✔️❌:     }
        // RDKit✔️❌:     ff.initialize();
        ff.initialize();
        // RDKit✔️❌:     int needsMore = ff.minimize(maxIters);
        let needs_more = ff.minimize(max_iters, 1.0e-4, 1.0e-6);
        // RDKit✔️❌:     double e = ff.calcEnergy();
        let energy = ff.calc_energy_current(None);
        // RDKit✔️❌:     res[i] = std::make_pair(needsMore, e);
        conformer_results[conf_index] = UffOptimizeMoleculeConfResult { needs_more, energy };
        {
            let coords =
                molecule.coordinate_block_mut().conformers_3d[conf_index].coordinates_mut();
            for (coord, position) in coords.iter_mut().zip(ff.positions()) {
                *coord = [position.x, position.y, position.z];
            }
        }
        // RDKit✔️❌:   }
    }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::detail::OptimizeMoleculeConfsST
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::OptimizeMoleculeConfs
    Ok(UffOptimizeMoleculeConfsResult {
        molecule,
        conformer_results,
    })
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION RDKit::UFF::UFFOptimizeMoleculeConfs
}

pub(crate) fn get_atom_label_for_uff(
    atom: &Atom,
    atom_index: usize,
    total_valence: i32,
    hybridization: Hybridization,
    atom_has_conjugated_bond: bool,
    tolerate_charge_mismatch: bool,
) -> Result<String, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::getAtomLabel (AtomTyper.cpp:415-498)
    // RDKit✔️✔️: std::string getAtomLabel(const Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // Rust references reproduce RDKit's non-null atom precondition.
    // RDKit✔️✔️:   int atNum = atom->getAtomicNum();
    let atomic_number = atom.atomic_number();
    // RDKit✔️✔️:   std::string atomKey = atom->getSymbol();
    let mut atom_key = rdkit_element_symbol(atomic_number)?.to_string();
    // RDKit✔️✔️:   if (atomKey.size() == 1) {
    if atom_key.len() == 1 {
        // RDKit✔️✔️:     atomKey += '_';
        atom_key.push('_');
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   PeriodicTable *table = PeriodicTable::getTable();
    // RDKit✔️✔️:   // FIX: handle main group/organometallic cases better:
    // RDKit✔️✔️:   if (atNum) {
    if atomic_number != 0 {
        // RDKit✔️✔️:     // do not do hybridization on alkali metals or halogens:
        // RDKit✔️✔️:     if (table->getDefaultValence(atNum) == -1 ||
        // RDKit✔️✔️:         (table->getNouterElecs(atNum) != 1 &&
        // RDKit✔️✔️:          table->getNouterElecs(atNum) != 7)) {
        let default_valence = rdkit_default_valence(atomic_number)?;
        let n_outer_electrons = periodic_table_outer_electrons(atomic_number)?;
        if default_valence == -1 || (n_outer_electrons != 1 && n_outer_electrons != 7) {
            // RDKit✔️✔️:       switch (atom->getAtomicNum()) {
            match atomic_number {
                // RDKit✔️✔️:         case 12:
                // RDKit✔️✔️:         case 13:
                // RDKit✔️✔️:         case 14:
                // RDKit✔️✔️:         case 15:
                // RDKit✔️✔️:         case 50:
                // RDKit✔️✔️:         case 51:
                // RDKit✔️✔️:         case 52:
                // RDKit✔️✔️:         case 81:
                // RDKit✔️✔️:         case 82:
                // RDKit✔️✔️:         case 83:
                // RDKit✔️✔️:         case 84:
                12 | 13 | 14 | 15 | 50 | 51 | 52 | 81 | 82 | 83 | 84 => {
                    // RDKit✔️✔️:           atomKey += '3';
                    atom_key.push('3');
                    // RDKit✔️✔️:           if (atom->getHybridization() != Atom::SP3) {
                    if hybridization != Hybridization::Sp3 {
                        // RDKit✔️✔️:             BOOST_LOG(rdWarningLog)
                        // RDKit✔️✔️:                 << "UFFTYPER: Warning: hybridization set to SP3 for atom "
                        // RDKit✔️✔️:                 << atom->getIdx() << std::endl;
                        let _ = atom_index;
                        // RDKit logs the warning; COSMolKit preserves atomKey behavior.
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:           break;
                }
                // RDKit✔️✔️:         case 80:
                80 => {
                    // RDKit✔️✔️:           atomKey += '1';
                    atom_key.push('1');
                    // RDKit✔️✔️:           if (atom->getHybridization() != Atom::SP) {
                    if hybridization != Hybridization::Sp {
                        // RDKit✔️✔️:             BOOST_LOG(rdWarningLog)
                        // RDKit✔️✔️:                 << "UFFTYPER: Warning: hybridization set to SP for atom "
                        // RDKit✔️✔️:                 << atom->getIdx() << std::endl;
                        let _ = atom_index;
                        // RDKit logs the warning; COSMolKit preserves atomKey behavior.
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:           break;
                }
                // RDKit✔️✔️:         default:
                _ => {
                    // RDKit✔️✔️:           switch (atom->getHybridization()) {
                    match hybridization {
                        // RDKit✔️✔️:             case Atom::S:
                        // RDKit✔️✔️:               // don't need to do anything here
                        // RDKit✔️✔️:               break;
                        Hybridization::S => {}
                        // RDKit✔️✔️:             case Atom::SP:
                        Hybridization::Sp => {
                            // RDKit✔️✔️:               atomKey += '1';
                            atom_key.push('1');
                            // RDKit✔️✔️:               break;
                        }
                        // RDKit✔️✔️:             case Atom::SP2:
                        Hybridization::Sp2 => {
                            // RDKit✔️✔️:               if ((atom->getIsAromatic() ||
                            // RDKit✔️✔️:                    MolOps::atomHasConjugatedBond(atom)) &&
                            // RDKit✔️✔️:                   (atNum == 6 || atNum == 7 || atNum == 8 || atNum == 16)) {
                            if (atom.is_aromatic() || atom_has_conjugated_bond)
                                && matches!(atomic_number, 6 | 7 | 8 | 16)
                            {
                                // RDKit✔️✔️:                 atomKey += 'R';
                                atom_key.push('R');
                                // RDKit✔️✔️:               } else {
                            } else {
                                // RDKit✔️✔️:                 atomKey += '2';
                                atom_key.push('2');
                                // RDKit✔️✔️:               }
                            }
                            // RDKit✔️✔️:               break;
                        }
                        // RDKit✔️✔️:             case Atom::SP3:
                        Hybridization::Sp3 => {
                            // RDKit✔️✔️:               atomKey += '3';
                            atom_key.push('3');
                            // RDKit✔️✔️:               break;
                        }
                        // RDKit✔️✔️:             case Atom::SP2D:
                        Hybridization::Sp2d => {
                            // RDKit✔️✔️:               atomKey += '4';
                            atom_key.push('4');
                            // RDKit✔️✔️:               break;
                        }
                        // RDKit✔️✔️:             case Atom::SP3D:
                        Hybridization::Sp3d => {
                            // RDKit✔️✔️:               atomKey += '5';
                            atom_key.push('5');
                            // RDKit✔️✔️:               break;
                        }
                        // RDKit✔️✔️:             case Atom::SP3D2:
                        Hybridization::Sp3d2 => {
                            // RDKit✔️✔️:               atomKey += '6';
                            atom_key.push('6');
                            // RDKit✔️✔️:               break;
                        }
                        // RDKit✔️✔️:             default:
                        Hybridization::Unspecified | Hybridization::Other => {
                            // RDKit✔️✔️:               BOOST_LOG(rdErrorLog)
                            // RDKit✔️✔️:                   << "UFFTYPER: Unrecognized hybridization for atom: "
                            // RDKit✔️✔️:                   << atom->getIdx() << std::endl;
                            let _ = atom_index;
                            // RDKit logs the error; COSMolKit preserves atomKey unchanged.
                            // RDKit✔️✔️:           }
                        }
                    }
                    // RDKit✔️✔️:       }
                }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   // special cases by element type:
    // RDKit✔️✔️:   addAtomChargeFlags(atom, atomKey);
    add_atom_charge_flags_for_uff(
        atom,
        atom_index,
        total_valence,
        &mut atom_key,
        hybridization,
        tolerate_charge_mismatch,
    );
    // RDKit✔️✔️:   return atomKey;
    Ok(atom_key)
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::getAtomLabel
}

pub(crate) fn add_atom_charge_flags_for_uff(
    atom: &Atom,
    atom_index: usize,
    total_valence: i32,
    atom_key: &mut String,
    hybridization: Hybridization,
    tolerate_charge_mismatch: bool,
) {
    // RDKit✔️✔️: void addAtomChargeFlags(const Atom *atom, std::string &atomKey,
    // RDKit✔️✔️:                         bool tolerateChargeMismatch) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // Rust references reproduce RDKit's non-null atom precondition.
    // RDKit✔️✔️:   int totalValence = atom->getTotalValence();
    // The caller supplies COSMolKit's RDKit-aligned total valence for the atom.
    // RDKit✔️✔️:   int fc = atom->getFormalCharge();
    let formal_charge = atom.formal_charge();
    // RDKit✔️✔️:   // FIX: come up with some way of handling metals here
    // RDKit✔️✔️:   switch (atom->getAtomicNum()) {
    match atom.atomic_number() {
        // RDKit✔️✔️:     // Atoms only +1 in default UFF params
        // RDKit✔️✔️:     case 29:  // Cu
        // RDKit✔️✔️:     case 47:  // Ag
        29 | 47 => {
            // RDKit✔️✔️:       if (totalValence == 1 || fc == 1 || tolerateChargeMismatch) {
            if total_valence == 1 || formal_charge == 1 || tolerate_charge_mismatch {
                // RDKit✔️✔️:         atomKey += "+1";
                atom_key.push_str("+1");
                // RDKit✔️✔️:       } else {
            } else {
                // RDKit✔️✔️:         BOOST_LOG(rdErrorLog)
                // RDKit✔️✔️:             << "UFFTYPER: Unrecognized charge state for atom: "
                // RDKit✔️✔️:             << atom->getIdx() << std::endl;
                let _ = atom_index;
                // RDKit logs the unrecognized state; COSMolKit preserves atomKey unchanged.
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:       break;
        }
        // RDKit✔️✔️:     // Atoms only +2 in default UFF params
        // RDKit✔️✔️:     case 4:   // Be
        // RDKit✔️✔️:     case 20:  // Ca
        // RDKit✔️✔️:     case 25:  // Mn
        // RDKit✔️✔️:     case 26:  // Fe
        // RDKit✔️✔️:     case 28:  // Ni
        // RDKit✔️✔️:               //  case 30:  // Zn
        // RDKit✔️✔️:     case 46:  // Pd
        // RDKit✔️✔️:               //  case 48:  // Cd
        // RDKit✔️✔️:     case 78:  // Pt
        4 | 20 | 25 | 26 | 28 | 46 | 78 => {
            // RDKit✔️✔️:       if (totalValence == 2 || fc == 2 || tolerateChargeMismatch) {
            if total_valence == 2 || formal_charge == 2 || tolerate_charge_mismatch {
                // RDKit✔️✔️:         atomKey += "+2";
                atom_key.push_str("+2");
                // RDKit✔️✔️:       } else {
            } else {
                // RDKit✔️✔️:         BOOST_LOG(rdErrorLog)
                // RDKit✔️✔️:             << "UFFTYPER: Unrecognized charge state for atom: "
                // RDKit✔️✔️:             << atom->getIdx() << std::endl;
                let _ = atom_index;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:       break;
        }
        // RDKit✔️✔️:     // Atoms only +3 in default UFF params
        // RDKit✔️✔️:     case 21:   // Sc
        // RDKit✔️✔️:     case 24:   // Cr
        // RDKit✔️✔️:     case 27:   // Co
        // RDKit✔️✔️:                //  case 49:  // In
        // RDKit✔️✔️:     case 79:   // Au
        // RDKit✔️✔️:     case 89:   // Ac
        // RDKit✔️✔️:     case 96:   // Cm
        // RDKit✔️✔️:     case 97:   // Bk
        // RDKit✔️✔️:     case 98:   // Cf
        // RDKit✔️✔️:     case 99:   // Es
        // RDKit✔️✔️:     case 100:  // Fm
        // RDKit✔️✔️:     case 101:  // Md
        // RDKit✔️✔️:     case 102:  // No
        // RDKit✔️✔️:     case 103:  // Lr/Lw
        21 | 24 | 27 | 79 | 89 | 96..=103 => {
            // RDKit✔️✔️:       if (totalValence == 3 || fc == 3 || tolerateChargeMismatch) {
            if total_valence == 3 || formal_charge == 3 || tolerate_charge_mismatch {
                // RDKit✔️✔️:         atomKey += "+3";
                atom_key.push_str("+3");
            }
            // RDKit✔️✔️:       break;
        }
        // RDKit✔️✔️:     // Atoms only +4 in default UFF params
        // RDKit✔️✔️:     case 2:   // He
        // RDKit✔️✔️:     case 18:  // Ar
        // RDKit✔️✔️:     case 22:  // Ti
        // RDKit✔️✔️:     case 36:  // Kr
        // RDKit✔️✔️:     case 54:  // Xe
        // RDKit✔️✔️:     case 90:  // Th
        // RDKit✔️✔️:     case 91:  // Pa
        // RDKit✔️✔️:     case 92:  // U
        // RDKit✔️✔️:     case 93:  // Np
        // RDKit✔️✔️:     case 94:  // Pu
        // RDKit✔️✔️:     case 95:  // Am
        2 | 18 | 22 | 36 | 54 | 90 | 91 | 92 | 93 | 94 | 95 => {
            // RDKit✔️✔️:       if (totalValence == 4 || fc == 4 || tolerateChargeMismatch) {
            if total_valence == 4 || formal_charge == 4 || tolerate_charge_mismatch {
                // RDKit✔️✔️:         atomKey += "+4";
                atom_key.push_str("+4");
            }
            // RDKit✔️✔️:       break;
        }
        // RDKit✔️✔️:     // Atoms only +5 in default UFF params
        // RDKit✔️✔️:     case 23:  // V
        // RDKit✔️✔️:     case 41:  // Nb
        // RDKit✔️✔️:     case 43:  // Tc
        // RDKit✔️✔️:     case 73:  // Ta
        23 | 41 | 43 | 73 => {
            // RDKit✔️✔️:       if (totalValence == 5 || fc == 5 || tolerateChargeMismatch) {
            if total_valence == 5 || formal_charge == 5 || tolerate_charge_mismatch {
                // RDKit✔️✔️:         atomKey += "+5";
                atom_key.push_str("+5");
            }
            // RDKit✔️✔️:       break;
        }
        // RDKit✔️✔️:     // Atoms only +6 in default UFF params
        // RDKit✔️✔️:     case 42:  // Mo
        42 => {
            // RDKit✔️✔️:       if (totalValence == 6 || fc == 6 || tolerateChargeMismatch) {
            if total_valence == 6 || formal_charge == 6 || tolerate_charge_mismatch {
                // RDKit✔️✔️:         atomKey += "+6";
                atom_key.push_str("+6");
            }
            // RDKit✔️✔️:       break;
        }
        // RDKit✔️✔️:     case 12:  // Mg
        // RDKit✔️✔️:       switch (totalValence) {
        12 => match total_valence {
            // RDKit✔️✔️:         case 2:
            2 => {
                // RDKit✔️✔️:           atomKey += "+2";
                atom_key.push_str("+2");
                // RDKit✔️✔️:           break;
            }
            // RDKit✔️✔️:         default:
            _ => {
                // RDKit✔️✔️:           if (tolerateChargeMismatch) {
                if tolerate_charge_mismatch {
                    // RDKit✔️✔️:             atomKey += "+2";
                    atom_key.push_str("+2");
                    // RDKit✔️✔️:           }
                }
                // RDKit✔️✔️:           BOOST_LOG(rdErrorLog)
                // RDKit✔️✔️:               << "UFFTYPER: Unrecognized charge state for atom: "
                // RDKit✔️✔️:               << atom->getIdx() << std::endl;
            } // RDKit✔️✔️:       }
              // RDKit✔️✔️:       break;
        },
        // RDKit✔️✔️:     case 13:  // Al
        13 => {
            // RDKit✔️✔️:       if (totalValence != 3) {
            // RDKit✔️✔️:         BOOST_LOG(rdErrorLog)
            // RDKit✔️✔️:             << "UFFTYPER: Unrecognized charge state for atom: "
            // RDKit✔️✔️:             << atom->getIdx() << std::endl;
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       break;
            let _ = total_valence != 3;
        }
        // RDKit✔️✔️:     case 14:  // Si
        14 => {
            // RDKit✔️✔️:       if (totalValence != 4) {
            // RDKit✔️✔️:         BOOST_LOG(rdErrorLog)
            // RDKit✔️✔️:             << "UFFTYPER: Unrecognized charge state for atom: "
            // RDKit✔️✔️:             << atom->getIdx() << std::endl;
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       break;
            let _ = total_valence != 4;
        }
        // RDKit✔️✔️:     case 15:  // P
        // RDKit✔️✔️:       switch (totalValence) {
        15 => match total_valence {
            // RDKit✔️✔️:         case 3:
            3 => {
                // RDKit✔️✔️:           atomKey += "+3";
                atom_key.push_str("+3");
            }
            // RDKit✔️✔️:         case 5:
            5 => {
                // RDKit✔️✔️:           atomKey += "+5";
                atom_key.push_str("+5");
            }
            // RDKit✔️✔️:         default:
            _ => {
                // RDKit✔️✔️:           if (tolerateChargeMismatch) {
                if tolerate_charge_mismatch {
                    // RDKit✔️✔️:             atomKey += "+5";
                    atom_key.push_str("+5");
                }
            }
        },
        // RDKit✔️✔️:     case 16:  // S
        16 => {
            // RDKit✔️✔️:       if (atom->getHybridization() != Atom::SP2) {
            if hybridization != Hybridization::Sp2 {
                // RDKit✔️✔️:         switch (totalValence) {
                match total_valence {
                    // RDKit✔️✔️:           case 2:
                    2 => atom_key.push_str("+2"),
                    // RDKit✔️✔️:           case 4:
                    4 => atom_key.push_str("+4"),
                    // RDKit✔️✔️:           case 6:
                    6 => atom_key.push_str("+6"),
                    // RDKit✔️✔️:           default:
                    _ => {
                        // RDKit✔️✔️:             if (tolerateChargeMismatch) {
                        if tolerate_charge_mismatch {
                            // RDKit✔️✔️:               atomKey += "+6";
                            atom_key.push_str("+6");
                        }
                    }
                }
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:       break;
        }
        // RDKit✔️✔️:     case 30:  // Zn
        // RDKit✔️✔️:     case 34:  // Se
        // RDKit✔️✔️:     case 48:  // Cd
        // RDKit✔️✔️:     case 52:  // Te
        // RDKit✔️✔️:     case 80:  // Hg
        // RDKit✔️✔️:     case 84:  // Po
        30 | 34 | 48 | 52 | 80 | 84 => match total_valence {
            // RDKit✔️✔️:         case 2:
            2 => atom_key.push_str("+2"),
            // RDKit✔️✔️:         default:
            _ if tolerate_charge_mismatch => atom_key.push_str("+2"),
            _ => {}
        },
        // RDKit✔️✔️:     case 31:  // Ga
        // RDKit✔️✔️:     case 33:  // As
        // RDKit✔️✔️:     case 49:  // In
        // RDKit✔️✔️:     case 51:  // Sb
        // RDKit✔️✔️:     case 81:  // Tl
        // RDKit✔️✔️:     case 82:  // Pb
        // RDKit✔️✔️:     case 83:  // Bi
        31 | 33 | 49 | 51 | 81 | 82 | 83 => match total_valence {
            // RDKit✔️✔️:         case 3:
            3 => atom_key.push_str("+3"),
            // RDKit✔️✔️:         default:
            _ if tolerate_charge_mismatch => atom_key.push_str("+3"),
            _ => {}
        },
        // RDKit✔️✔️:     case 75:  // Re
        75 => {
            // RDKit✔️✔️:       if (tolerateChargeMismatch) {
            if tolerate_charge_mismatch {
                // RDKit✔️✔️:         if (atomKey == "Re6") {
                if atom_key == "Re6" {
                    // RDKit✔️✔️:           atomKey = "Re6+5";
                    *atom_key = "Re6+5".to_string();
                    // RDKit✔️✔️:         } else if (atomKey == "Re3") {
                } else if atom_key == "Re3" {
                    // RDKit✔️✔️:           atomKey = "Re3+7";
                    *atom_key = "Re3+7".to_string();
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:       BOOST_LOG(rdErrorLog)
            // RDKit✔️✔️:           << "UFFTYPER: Unrecognized charge state for atom: " << atom->getIdx()
            // RDKit✔️✔️:           << std::endl;
            // RDKit✔️✔️:       break;
        }
        // RDKit✔️✔️:   }
        _ => {}
    }
    // RDKit✔️✔️:   // lanthanides
    // RDKit✔️✔️:   if (atom->getAtomicNum() >= 57 && atom->getAtomicNum() <= 71) {
    if (57..=71).contains(&atom.atomic_number()) {
        // RDKit✔️✔️:     switch (totalValence) {
        match total_valence {
            // RDKit✔️✔️:       case 6:
            6 => {
                // RDKit✔️✔️:         atomKey += "+3";
                atom_key.push_str("+3");
                // RDKit✔️✔️:         break;
            }
            // RDKit✔️✔️:       default:
            _ => {
                // RDKit✔️✔️:         if (tolerateChargeMismatch) {
                if tolerate_charge_mismatch {
                    // RDKit✔️✔️:           atomKey += "+3";
                    atom_key.push_str("+3");
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:         BOOST_LOG(rdErrorLog)
                // RDKit✔️✔️:             << "UFFTYPER: Unrecognized charge state for atom: "
                // RDKit✔️✔️:             << atom->getIdx() << std::endl;
            } // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
}

pub(crate) fn get_atom_types_for_uff(
    mol: &Molecule,
    total_valences: &[i32],
    hybridizations: &[Hybridization],
    atom_has_conjugated_bond: &[bool],
) -> Result<(Vec<Option<AtomicParams>>, bool), UffAtomTyperError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getAtomTypes (AtomTyper.cpp:502-528)
    // RDKit✔️✔️: std::pair<AtomicParamVect, bool> getAtomTypes(const ROMol &mol,
    // RDKit✔️✔️:                                               const std::string &) {
    let atom_count = mol.atoms().len();
    if total_valences.len() != atom_count
        || hybridizations.len() != atom_count
        || atom_has_conjugated_bond.len() != atom_count
    {
        return Err(UffAtomTyperError::ContextLengthMismatch {
            atoms: atom_count,
            total_valences: total_valences.len(),
            hybridizations: hybridizations.len(),
            conjugation_flags: atom_has_conjugated_bond.len(),
        });
    }
    // RDKit✔️✔️:   bool foundAll = true;
    let mut found_all = true;
    // RDKit✔️✔️:   auto params = ParamCollection::getParams();
    let params = ParamCollection::get_params("")?;
    // RDKit✔️✔️:   AtomicParamVect paramVect;
    // RDKit✔️✔️:   paramVect.resize(mol.getNumAtoms());
    let mut param_vect = Vec::with_capacity(atom_count);
    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    for atom_index in 0..atom_count {
        // RDKit✔️✔️:     const Atom *atom = mol.getAtomWithIdx(i);
        let atom = &mol.atoms()[atom_index];
        // RDKit✔️✔️:     // construct the atom key:
        // RDKit✔️✔️:     std::string atomKey = Tools::getAtomLabel(atom);
        let atom_key = get_atom_label_for_uff(
            atom,
            atom_index,
            total_valences[atom_index],
            hybridizations[atom_index],
            atom_has_conjugated_bond[atom_index],
            true,
        )?;
        // RDKit✔️✔️:     // ok, we've got the atom key, now get the parameters:
        // RDKit✔️✔️:     const AtomicParams *theParams = (*params)(atomKey);
        let the_params = params.get(&atom_key).copied();
        // RDKit✔️✔️:     if (!theParams) {
        if the_params.is_none() {
            // RDKit✔️✔️:       foundAll = false;
            found_all = false;
            // RDKit✔️✔️:       BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized atom type: " << atomKey
            // RDKit✔️✔️:                             << " (" << i << ")" << std::endl;
            // RDKit logs the missing type; COSMolKit preserves the null slot as None.
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     paramVect[i] = theParams;
        param_vect.push(the_params);
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return std::make_pair(paramVect, foundAll);
    Ok((param_vect, found_all))
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getAtomTypes
}

#[cfg(test)]
mod tests {
    use crate::builder::MoleculeBuilder;
    use crate::chemistry::forcefield::uff::utils::{
        calc_angle_force_constant, calc_bond_force_constant, calc_bond_rest_length,
        calc_inversion_coefficients, calc_nonbonded_depth, calc_nonbonded_minimum,
    };
    use crate::{AtomId, AtomSpec, BondOrder, BondSpec, Conformer3D, Element};

    use super::*;

    fn atom(atomic_num: u8, formal_charge: i8) -> Atom {
        atom_with_aromatic(atomic_num, formal_charge, false)
    }

    fn atom_with_aromatic(atomic_num: u8, formal_charge: i8, is_aromatic: bool) -> Atom {
        let element = Element::from_atomic_number(atomic_num).expect("test atomic number");
        Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(element)
                .with_formal_charge(formal_charge)
                .with_aromatic(is_aromatic),
        )
    }

    fn label(
        atomic_num: u8,
        total_valence: i32,
        formal_charge: i8,
        hybridization: Hybridization,
        atom_has_conjugated_bond: bool,
        tolerate_charge_mismatch: bool,
        is_aromatic: bool,
    ) -> String {
        let atom = atom_with_aromatic(atomic_num, formal_charge, is_aromatic);
        get_atom_label_for_uff(
            &atom,
            0,
            total_valence,
            hybridization,
            atom_has_conjugated_bond,
            tolerate_charge_mismatch,
        )
        .expect("UFF atom label")
    }

    fn molecule_with_atomic_numbers(atomic_numbers: &[u8]) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        for &atomic_number in atomic_numbers {
            let element = Element::from_atomic_number(atomic_number).expect("test atomic number");
            builder.add_atom(AtomSpec::new(element));
        }
        builder.build().expect("test molecule")
    }

    fn bonded_pair(
        first: AtomSpec,
        second: AtomSpec,
        order: BondOrder,
        is_conjugated: bool,
    ) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(first);
        let a1 = builder.add_atom(second);
        builder
            .add_bond(BondSpec::new(a0, a1, order).with_conjugated(is_conjugated))
            .expect("test bond");
        builder.build().expect("test bonded pair")
    }

    fn bonded_pair_with_3d_conformer(
        first: AtomSpec,
        second: AtomSpec,
        order: BondOrder,
        coords: Vec<[f64; 3]>,
    ) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(first);
        let a1 = builder.add_atom(second);
        builder
            .add_bond(BondSpec::new(a0, a1, order))
            .expect("test bond");
        builder.add_3d_conformer(coords).expect("test 3d conformer");
        builder.build().expect("test bonded pair with conformer")
    }

    fn bonded_pair_with_named_3d_conformers(
        first: AtomSpec,
        second: AtomSpec,
        order: BondOrder,
        first_coords: Vec<[f64; 3]>,
        second_coords: Vec<[f64; 3]>,
    ) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(first);
        let a1 = builder.add_atom(second);
        builder
            .add_bond(BondSpec::new(a0, a1, order))
            .expect("test bond");
        builder
            .add_conformer(Conformer3D::new(0, first_coords, true))
            .expect("first conformer");
        builder
            .add_conformer(Conformer3D::new(7, second_coords, true))
            .expect("second conformer");
        builder
            .build()
            .expect("test bonded pair with named conformers")
    }

    fn bonded_triple(
        first: AtomSpec,
        second: AtomSpec,
        third: AtomSpec,
        order12: BondOrder,
        order23: BondOrder,
        is_conjugated12: bool,
        is_conjugated23: bool,
    ) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(first);
        let a1 = builder.add_atom(second);
        let a2 = builder.add_atom(third);
        builder
            .add_bond(BondSpec::new(a0, a1, order12).with_conjugated(is_conjugated12))
            .expect("test bond 0-1");
        builder
            .add_bond(BondSpec::new(a1, a2, order23).with_conjugated(is_conjugated23))
            .expect("test bond 1-2");
        builder.build().expect("test bonded triple")
    }

    fn bonded_quad(
        first: AtomSpec,
        second: AtomSpec,
        third: AtomSpec,
        fourth: AtomSpec,
        order12: BondOrder,
        order23: BondOrder,
        order34: BondOrder,
        is_conjugated12: bool,
        is_conjugated23: bool,
        is_conjugated34: bool,
    ) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(first);
        let a1 = builder.add_atom(second);
        let a2 = builder.add_atom(third);
        let a3 = builder.add_atom(fourth);
        builder
            .add_bond(BondSpec::new(a0, a1, order12).with_conjugated(is_conjugated12))
            .expect("test bond 0-1");
        builder
            .add_bond(BondSpec::new(a1, a2, order23).with_conjugated(is_conjugated23))
            .expect("test bond 1-2");
        builder
            .add_bond(BondSpec::new(a2, a3, order34).with_conjugated(is_conjugated34))
            .expect("test bond 2-3");
        builder.build().expect("test bonded quad")
    }

    fn trigonal_center_molecule(
        center: AtomSpec,
        neighbors: [AtomSpec; 3],
        bond_orders: [BondOrder; 3],
    ) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let center_id = builder.add_atom(center);
        let neighbor_ids: Vec<_> = neighbors
            .into_iter()
            .map(|neighbor| builder.add_atom(neighbor))
            .collect();
        for (neighbor_id, bond_order) in neighbor_ids.into_iter().zip(bond_orders) {
            builder
                .add_bond(BondSpec::new(center_id, neighbor_id, bond_order))
                .expect("trigonal-center bond should build");
        }
        builder
            .build()
            .expect("trigonal-center molecule should build")
    }

    fn flagged(
        atomic_num: u8,
        total_valence: i32,
        formal_charge: i8,
        hybridization: Hybridization,
        tolerate_charge_mismatch: bool,
        start: &str,
    ) -> String {
        let atom = atom(atomic_num, formal_charge);
        let mut atom_key = start.to_string();
        add_atom_charge_flags_for_uff(
            &atom,
            0,
            total_valence,
            &mut atom_key,
            hybridization,
            tolerate_charge_mismatch,
        );
        atom_key
    }

    #[test]
    fn uff_atom_typer_add_charge_flags_handles_fixed_positive_charge_families() {
        let cases = [
            (29, 1, 0, "Cu+1"),
            (47, 0, 1, "Ag+1"),
            (4, 2, 0, "Be+2"),
            (78, 0, 2, "Pt+2"),
            (21, 3, 0, "Sc+3"),
            (103, 0, 3, "Lr+3"),
            (2, 4, 0, "He+4"),
            (95, 0, 4, "Am+4"),
            (23, 5, 0, "V+5"),
            (73, 0, 5, "Ta+5"),
            (42, 6, 0, "Mo+6"),
            (42, 0, 6, "Mo+6"),
        ];

        for (atomic_num, total_valence, formal_charge, expected) in cases {
            assert_eq!(
                flagged(
                    atomic_num,
                    total_valence,
                    formal_charge,
                    Hybridization::Unspecified,
                    false,
                    &expected[..expected.len() - 2],
                ),
                expected
            );
        }
    }

    #[test]
    fn uff_atom_typer_add_charge_flags_tolerates_fixed_family_mismatches() {
        assert_eq!(
            flagged(29, 0, 0, Hybridization::Unspecified, true, "Cu"),
            "Cu+1"
        );
        assert_eq!(
            flagged(4, 0, 0, Hybridization::Unspecified, true, "Be"),
            "Be+2"
        );
        assert_eq!(
            flagged(21, 0, 0, Hybridization::Unspecified, true, "Sc"),
            "Sc+3"
        );
        assert_eq!(
            flagged(2, 0, 0, Hybridization::Unspecified, true, "He"),
            "He+4"
        );
        assert_eq!(
            flagged(23, 0, 0, Hybridization::Unspecified, true, "V"),
            "V+5"
        );
        assert_eq!(
            flagged(42, 0, 0, Hybridization::Unspecified, true, "Mo"),
            "Mo+6"
        );
    }

    #[test]
    fn uff_atom_typer_add_charge_flags_leaves_unrecognized_fixed_family_mismatches_unchanged() {
        assert_eq!(
            flagged(29, 0, 0, Hybridization::Unspecified, false, "Cu"),
            "Cu"
        );
        assert_eq!(
            flagged(4, 0, 0, Hybridization::Unspecified, false, "Be"),
            "Be"
        );
        assert_eq!(
            flagged(21, 0, 0, Hybridization::Unspecified, false, "Sc"),
            "Sc"
        );
        assert_eq!(
            flagged(2, 0, 0, Hybridization::Unspecified, false, "He"),
            "He"
        );
        assert_eq!(
            flagged(23, 0, 0, Hybridization::Unspecified, false, "V"),
            "V"
        );
        assert_eq!(
            flagged(42, 0, 0, Hybridization::Unspecified, false, "Mo"),
            "Mo"
        );
    }

    #[test]
    fn uff_atom_typer_add_charge_flags_handles_mg_al_si_and_phosphorus() {
        assert_eq!(
            flagged(12, 2, 0, Hybridization::Unspecified, false, "Mg"),
            "Mg+2"
        );
        assert_eq!(
            flagged(12, 0, 0, Hybridization::Unspecified, true, "Mg"),
            "Mg+2"
        );
        assert_eq!(
            flagged(12, 0, 0, Hybridization::Unspecified, false, "Mg"),
            "Mg"
        );
        assert_eq!(
            flagged(13, 0, 0, Hybridization::Unspecified, true, "Al"),
            "Al"
        );
        assert_eq!(
            flagged(14, 0, 0, Hybridization::Unspecified, true, "Si"),
            "Si"
        );
        assert_eq!(
            flagged(15, 3, 0, Hybridization::Unspecified, false, "P_"),
            "P_+3"
        );
        assert_eq!(
            flagged(15, 5, 0, Hybridization::Unspecified, false, "P_"),
            "P_+5"
        );
        assert_eq!(
            flagged(15, 0, 0, Hybridization::Unspecified, true, "P_"),
            "P_+5"
        );
        assert_eq!(
            flagged(15, 0, 0, Hybridization::Unspecified, false, "P_"),
            "P_"
        );
    }

    #[test]
    fn uff_atom_typer_add_charge_flags_handles_sulfur_hybridization_and_valence() {
        assert_eq!(flagged(16, 2, 0, Hybridization::Sp3, false, "S_"), "S_+2");
        assert_eq!(flagged(16, 4, 0, Hybridization::Sp3, false, "S_"), "S_+4");
        assert_eq!(flagged(16, 6, 0, Hybridization::Sp3, false, "S_"), "S_+6");
        assert_eq!(flagged(16, 0, 0, Hybridization::Sp3, true, "S_"), "S_+6");
        assert_eq!(flagged(16, 0, 0, Hybridization::Sp3, false, "S_"), "S_");
        assert_eq!(flagged(16, 6, 0, Hybridization::Sp2, true, "S_"), "S_");
    }

    #[test]
    fn uff_atom_typer_add_charge_flags_handles_post_transition_special_groups() {
        for atomic_num in [30, 34, 48, 52, 80, 84] {
            assert_eq!(
                flagged(atomic_num, 2, 0, Hybridization::Unspecified, false, "X"),
                "X+2"
            );
            assert_eq!(
                flagged(atomic_num, 0, 0, Hybridization::Unspecified, true, "X"),
                "X+2"
            );
            assert_eq!(
                flagged(atomic_num, 0, 0, Hybridization::Unspecified, false, "X"),
                "X"
            );
        }
        for atomic_num in [31, 33, 49, 51, 81, 82, 83] {
            assert_eq!(
                flagged(atomic_num, 3, 0, Hybridization::Unspecified, false, "X"),
                "X+3"
            );
            assert_eq!(
                flagged(atomic_num, 0, 0, Hybridization::Unspecified, true, "X"),
                "X+3"
            );
            assert_eq!(
                flagged(atomic_num, 0, 0, Hybridization::Unspecified, false, "X"),
                "X"
            );
        }
    }

    #[test]
    fn uff_atom_typer_add_charge_flags_handles_rhenium_rewrites() {
        assert_eq!(
            flagged(75, 0, 0, Hybridization::Unspecified, true, "Re6"),
            "Re6+5"
        );
        assert_eq!(
            flagged(75, 0, 0, Hybridization::Unspecified, true, "Re3"),
            "Re3+7"
        );
        assert_eq!(
            flagged(75, 0, 0, Hybridization::Unspecified, true, "Re4"),
            "Re4"
        );
        assert_eq!(
            flagged(75, 0, 0, Hybridization::Unspecified, false, "Re6"),
            "Re6"
        );
    }

    #[test]
    fn uff_atom_typer_add_charge_flags_handles_lanthanide_boundaries() {
        assert_eq!(
            flagged(57, 6, 0, Hybridization::Unspecified, false, "La"),
            "La+3"
        );
        assert_eq!(
            flagged(71, 0, 0, Hybridization::Unspecified, true, "Lu"),
            "Lu+3"
        );
        assert_eq!(
            flagged(57, 0, 0, Hybridization::Unspecified, false, "La"),
            "La"
        );
        assert_eq!(
            flagged(56, 6, 0, Hybridization::Unspecified, true, "Ba"),
            "Ba"
        );
        assert_eq!(
            flagged(72, 6, 0, Hybridization::Unspecified, true, "Hf"),
            "Hf"
        );
    }

    #[test]
    fn uff_atom_typer_add_charge_flags_leaves_default_atoms_unchanged() {
        assert_eq!(flagged(6, 4, 0, Hybridization::Sp3, true, "C_3"), "C_3");
    }

    #[test]
    fn uff_atom_typer_get_atom_label_pads_one_character_symbols_and_skips_atomic_zero() {
        assert_eq!(label(0, 0, 0, Hybridization::Sp3, false, true, false), "*_");
        assert_eq!(label(6, 0, 0, Hybridization::S, false, true, false), "C_");
    }

    #[test]
    fn uff_atom_typer_get_atom_label_skips_hybridization_for_alkali_metals_and_halogens() {
        assert_eq!(label(1, 0, 0, Hybridization::Sp3, false, true, false), "H_");
        assert_eq!(label(3, 0, 0, Hybridization::Sp3, false, true, false), "Li");
        assert_eq!(
            label(17, 0, 0, Hybridization::Sp3, false, true, false),
            "Cl"
        );
    }

    #[test]
    fn uff_atom_typer_get_atom_label_handles_source_forced_hybridization_elements() {
        assert_eq!(
            label(12, 0, 0, Hybridization::Sp2, false, false, false),
            "Mg3"
        );
        assert_eq!(
            label(15, 0, 0, Hybridization::Sp2, false, false, false),
            "P_3"
        );
        assert_eq!(
            label(50, 0, 0, Hybridization::Sp, false, false, false),
            "Sn3"
        );
        assert_eq!(
            label(80, 0, 0, Hybridization::Sp3, false, false, false),
            "Hg1"
        );
    }

    #[test]
    fn uff_atom_typer_get_atom_label_maps_hybridization_suffixes() {
        let cases = [
            (Hybridization::S, "C_"),
            (Hybridization::Sp, "C_1"),
            (Hybridization::Sp2, "C_2"),
            (Hybridization::Sp3, "C_3"),
            (Hybridization::Sp2d, "C_4"),
            (Hybridization::Sp3d, "C_5"),
            (Hybridization::Sp3d2, "C_6"),
            (Hybridization::Unspecified, "C_"),
            (Hybridization::Other, "C_"),
        ];

        for (hybridization, expected) in cases {
            assert_eq!(label(6, 0, 0, hybridization, false, false, false), expected);
        }
    }

    #[test]
    fn uff_atom_typer_get_atom_label_uses_aromatic_or_conjugated_sp2_r_suffix() {
        assert_eq!(
            label(6, 0, 0, Hybridization::Sp2, false, false, true),
            "C_R"
        );
        assert_eq!(
            label(7, 0, 0, Hybridization::Sp2, true, false, false),
            "N_R"
        );
        assert_eq!(
            label(8, 0, 0, Hybridization::Sp2, false, false, true),
            "O_R"
        );
        assert_eq!(
            label(16, 0, 0, Hybridization::Sp2, true, false, false),
            "S_R"
        );
        assert_eq!(label(9, 0, 0, Hybridization::Sp2, true, false, true), "F_");
    }

    #[test]
    fn uff_atom_typer_get_atom_label_composes_charge_flags_after_label_suffixes() {
        assert_eq!(
            label(12, 2, 0, Hybridization::Sp3, false, false, false),
            "Mg3+2"
        );
        assert_eq!(
            label(29, 1, 0, Hybridization::Sp3, false, false, false),
            "Cu3+1"
        );
        assert_eq!(
            label(75, 0, 0, Hybridization::Sp3d2, false, true, false),
            "Re6+5"
        );
    }

    #[test]
    fn uff_atom_typer_get_atom_label_returns_error_for_out_of_range_atomic_number() {
        let atom = atom(119, 0);
        let err = get_atom_label_for_uff(&atom, 0, 0, Hybridization::Sp3, false, true)
            .expect_err("out-of-range atomic number should fail");
        assert!(err.to_string().contains("out of range"));
    }

    #[test]
    fn uff_atom_typer_get_atom_types_returns_params_and_found_all_for_known_labels() {
        let mol = molecule_with_atomic_numbers(&[6, 8]);
        let (params, found_all) = get_atom_types_for_uff(
            &mol,
            &[0, 0],
            &[Hybridization::Sp3, Hybridization::Sp2],
            &[false, false],
        )
        .expect("UFF atom types");

        assert!(found_all);
        assert_eq!(params.len(), 2);
        assert_eq!(params[0].expect("C_3 params").r1, 0.757);
        assert_eq!(params[1].expect("O_2 params").r1, 0.634);
    }

    #[test]
    fn uff_atom_typer_get_atom_types_preserves_missing_param_slots_and_found_all_false() {
        let mol = molecule_with_atomic_numbers(&[6, 9]);
        let (params, found_all) = get_atom_types_for_uff(
            &mol,
            &[0, 0],
            &[Hybridization::S, Hybridization::Sp3],
            &[false, false],
        )
        .expect("UFF atom types");

        assert!(!found_all);
        assert_eq!(params.len(), 2);
        assert!(params[0].is_none());
        assert_eq!(params[1].expect("F_ params").r1, 0.668);
    }

    #[test]
    fn uff_atom_typer_get_atom_types_rejects_context_length_mismatch() {
        let mol = molecule_with_atomic_numbers(&[6, 8]);
        let err = get_atom_types_for_uff(
            &mol,
            &[0],
            &[Hybridization::Sp3, Hybridization::Sp2],
            &[false],
        )
        .expect_err("context mismatch should fail");

        assert!(matches!(
            err,
            UffAtomTyperError::ContextLengthMismatch {
                atoms: 2,
                total_valences: 1,
                hybridizations: 2,
                conjugation_flags: 1
            }
        ));
    }

    #[test]
    fn uff_atom_typer_get_atom_types_propagates_atom_label_errors() {
        let mol = molecule_with_atomic_numbers(&[119]);
        let err = get_atom_types_for_uff(&mol, &[0], &[Hybridization::Sp3], &[false])
            .expect_err("out-of-range atom should fail");

        assert!(matches!(err, UffAtomTyperError::Valence(_)));
    }

    #[test]
    fn uff_public_api_get_uff_bond_stretch_params_returns_expected_values_for_known_bond() {
        let mol = bonded_pair(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::O).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            false,
        );

        let params = get_uff_bond_stretch_params(&mol, 0, 1)
            .expect("UFF bond stretch params should compute")
            .expect("bonded typed atoms should return params");

        let collection = ParamCollection::get_params("").expect("params");
        let c3 = collection.get("C_3").expect("C_3 params");
        let o3 = collection.get("O_3").expect("O_3 params");
        let expected_r0 = calc_bond_rest_length(1.0, c3, o3).expect("rest length");
        let expected_kb = calc_bond_force_constant(expected_r0, c3, o3).expect("force constant");

        assert_eq!(
            params,
            UffBond {
                kb: expected_kb,
                r0: expected_r0
            }
        );
    }

    #[test]
    fn uff_public_api_get_uff_bond_stretch_params_returns_none_when_bond_is_absent() {
        let mol = molecule_with_atomic_numbers(&[6, 8]);

        let params =
            get_uff_bond_stretch_params(&mol, 0, 1).expect("bond lookup without edge should work");

        assert!(params.is_none());
    }

    #[test]
    fn uff_public_api_get_uff_bond_stretch_params_returns_none_for_missing_atom_params() {
        let mol = bonded_pair(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::S),
            AtomSpec::new(Element::F).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            false,
        );

        let params = get_uff_bond_stretch_params(&mol, 0, 1)
            .expect("missing atom params should not error structurally");

        assert!(params.is_none());
    }

    #[test]
    fn uff_public_api_get_uff_bond_stretch_params_rejects_invalid_atom_index() {
        let mol = bonded_pair(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            false,
        );

        let err = get_uff_bond_stretch_params(&mol, 0, 2)
            .expect_err("out-of-range atom index should fail");

        assert_eq!(
            err,
            UffPublicApiError::InvalidAtomIndex {
                atom_index: 2,
                atoms: 2
            }
        );
    }

    #[test]
    fn uff_public_api_get_uff_angle_bend_params_returns_expected_values_for_known_angle() {
        let mol = bonded_triple(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::O).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            BondOrder::Single,
            false,
            false,
        );

        let params = get_uff_angle_bend_params(&mol, 0, 1, 2)
            .expect("UFF angle bend params should compute")
            .expect("connected typed atoms should return params");

        let collection = ParamCollection::get_params("").expect("params");
        let c3 = collection.get("C_3").expect("C_3 params");
        let o3 = collection.get("O_3").expect("O_3 params");
        let expected_theta0 = RAD2DEG * c3.theta0;
        let expected_ka =
            calc_angle_force_constant(c3.theta0, 1.0, 1.0, c3, c3, o3).expect("force constant");

        assert_eq!(
            params,
            UffAngle {
                ka: expected_ka,
                theta0: expected_theta0
            }
        );
    }

    #[test]
    fn uff_public_api_get_uff_angle_bend_params_returns_none_when_path_is_broken() {
        let mol = bonded_pair(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            false,
        );

        let params = get_uff_angle_bend_params(&mol, 0, 1, 2)
            .expect_err("out-of-range atom index should fail before missing path");

        assert_eq!(
            params,
            UffPublicApiError::InvalidAtomIndex {
                atom_index: 2,
                atoms: 2
            }
        );
    }

    #[test]
    fn uff_public_api_get_uff_angle_bend_params_returns_none_for_missing_bond() {
        let mol = molecule_with_atomic_numbers(&[6, 6, 8]);

        let params = get_uff_angle_bend_params(&mol, 0, 1, 2)
            .expect("disconnected angle lookup should work");

        assert!(params.is_none());
    }

    #[test]
    fn uff_public_api_get_uff_angle_bend_params_returns_none_for_missing_atom_params() {
        let mol = bonded_triple(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::S),
            AtomSpec::new(Element::O).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            BondOrder::Single,
            false,
            false,
        );

        let params = get_uff_angle_bend_params(&mol, 0, 1, 2)
            .expect("missing atom params should not error structurally");

        assert!(params.is_none());
    }

    #[test]
    fn uff_public_api_get_uff_angle_bend_params_rejects_invalid_atom_index() {
        let mol = bonded_triple(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::O).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            BondOrder::Single,
            false,
            false,
        );

        let err = get_uff_angle_bend_params(&mol, 0, 1, 3)
            .expect_err("out-of-range atom index should fail");

        assert_eq!(
            err,
            UffPublicApiError::InvalidAtomIndex {
                atom_index: 3,
                atoms: 3
            }
        );
    }

    #[test]
    fn uff_public_api_get_uff_torsion_params_returns_expected_value_for_sp3_sp3_case() {
        let mol = bonded_quad(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            false,
            false,
            false,
        );

        let params = get_uff_torsion_params(&mol, 0, 1, 2, 3)
            .expect("UFF torsion params should compute")
            .expect("connected typed atoms should return params");

        let collection = ParamCollection::get_params("").expect("params");
        let c3 = collection.get("C_3").expect("C_3 params");
        assert_eq!(
            params,
            UffTor {
                v: (c3.v1 * c3.v1).sqrt()
            }
        );
    }

    #[test]
    fn uff_public_api_get_uff_torsion_params_handles_group6_single_bond_special_case() {
        let mol = bonded_quad(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::O).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::S).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            false,
            false,
            false,
        );

        let params = get_uff_torsion_params(&mol, 0, 1, 2, 3)
            .expect("group6 torsion params should compute")
            .expect("group6 torsion should return params");

        assert_eq!(
            params,
            UffTor {
                v: (2.0_f64 * 6.8_f64).sqrt()
            }
        );
    }

    #[test]
    fn uff_public_api_get_uff_torsion_params_handles_sp2_sp2_and_sp2_sp3_cases() {
        let collection = ParamCollection::get_params("").expect("params");
        let c2 = collection.get("C_2").expect("C_2 params");
        let c3 = collection.get("C_3").expect("C_3 params");

        let sp2_sp2 = bonded_quad(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            false,
            true,
            false,
        );
        let sp2_sp2_params = get_uff_torsion_params(&sp2_sp2, 0, 1, 2, 3)
            .expect("sp2-sp2 torsion params should compute")
            .expect("sp2-sp2 torsion should return params");
        assert_eq!(
            sp2_sp2_params,
            UffTor {
                v: super::super::torsion_angle::equation17(2.0, c2, c2)
            }
        );

        let sp2_sp3_default = bonded_quad(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            false,
            false,
            false,
        );
        let sp2_sp3_default_params = get_uff_torsion_params(&sp2_sp3_default, 0, 1, 2, 3)
            .expect("sp2-sp3 torsion params should compute")
            .expect("sp2-sp3 torsion should return params");
        assert_eq!(sp2_sp3_default_params, UffTor { v: 1.0 });

        let sp2_sp3_end_sp2 = bonded_quad(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Single,
            true,
            false,
            false,
        );
        let sp2_sp3_end_sp2_params = get_uff_torsion_params(&sp2_sp3_end_sp2, 0, 1, 2, 3)
            .expect("sp2-sp3 end-sp2 torsion params should compute")
            .expect("sp2-sp3 end-sp2 torsion should return params");
        assert_eq!(sp2_sp3_end_sp2_params, UffTor { v: 2.0 });

        let sp3_group6_non_group6_sp2 = bonded_quad(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::O).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            false,
            false,
            false,
        );
        let sp3_group6_non_group6_sp2_params =
            get_uff_torsion_params(&sp3_group6_non_group6_sp2, 0, 1, 2, 3)
                .expect("group6 sp3/non-group6 sp2 torsion params should compute")
                .expect("group6 sp3/non-group6 sp2 torsion should return params");
        assert_eq!(
            sp3_group6_non_group6_sp2_params,
            UffTor {
                v: super::super::torsion_angle::equation17(1.0, collection.get("O_3").unwrap(), c2)
            }
        );

        let _ = c3;
    }

    #[test]
    fn uff_public_api_get_uff_torsion_params_returns_none_for_missing_bond_or_params_or_bad_hybridization()
     {
        let disconnected = molecule_with_atomic_numbers(&[6, 6, 6, 6]);
        assert!(
            get_uff_torsion_params(&disconnected, 0, 1, 2, 3)
                .expect("disconnected torsion lookup should work")
                .is_none()
        );

        let missing_params = bonded_quad(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::S),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            false,
            false,
            false,
        );
        assert!(
            get_uff_torsion_params(&missing_params, 0, 1, 2, 3)
                .expect("missing central params should not error structurally")
                .is_none()
        );

        let bad_hybridization = bonded_quad(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3d),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            false,
            false,
            false,
        );
        assert!(
            get_uff_torsion_params(&bad_hybridization, 0, 1, 2, 3)
                .expect("unsupported central hybridization should not error structurally")
                .is_none()
        );
    }

    #[test]
    fn uff_public_api_get_uff_torsion_params_rejects_invalid_atom_index() {
        let mol = bonded_quad(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            false,
            false,
            false,
        );

        let err = get_uff_torsion_params(&mol, 0, 1, 2, 4)
            .expect_err("out-of-range atom index should fail");

        assert_eq!(
            err,
            UffPublicApiError::InvalidAtomIndex {
                atom_index: 4,
                atoms: 4
            }
        );
    }

    #[test]
    fn uff_public_api_get_uff_inversion_params_returns_expected_force_constant_for_sp2_carbon() {
        let mol = trigonal_center_molecule(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            [
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );

        let params = get_uff_inversion_params(&mol, 1, 0, 2, 3)
            .expect("sp2 trigonal center should compute inversion params")
            .expect("supported inversion center should return params");

        let (k, _, _, _) = calc_inversion_coefficients(6, false);
        assert_eq!(params, UffInv { k });
    }

    #[test]
    fn uff_public_api_get_uff_inversion_params_detects_carbon_bound_to_sp2_oxygen() {
        let mol = trigonal_center_molecule(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            [
                AtomSpec::new(Element::O).with_hybridization(Hybridization::Sp2),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            ],
            [BondOrder::Double, BondOrder::Single, BondOrder::Single],
        );

        let params = get_uff_inversion_params(&mol, 1, 0, 2, 3)
            .expect("carboxyl-like center should compute inversion params")
            .expect("supported inversion center should return params");

        let (k, _, _, _) = calc_inversion_coefficients(6, true);
        assert_eq!(params, UffInv { k });
    }

    #[test]
    fn uff_public_api_get_uff_inversion_params_allows_group15_center_without_sp2_requirement() {
        let mol = trigonal_center_molecule(
            AtomSpec::new(Element::P).with_hybridization(Hybridization::Sp3),
            [
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );

        let params = get_uff_inversion_params(&mol, 1, 0, 2, 3)
            .expect("group15 center should compute inversion params")
            .expect("group15 center should return params");

        let (k, _, _, _) = calc_inversion_coefficients(15, false);
        assert_eq!(params, UffInv { k });
    }

    #[test]
    fn uff_public_api_get_uff_inversion_params_returns_none_for_unsupported_or_invalid_topology() {
        let unsupported_center = trigonal_center_molecule(
            AtomSpec::new(Element::S).with_hybridization(Hybridization::Sp2),
            [
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        assert!(
            get_uff_inversion_params(&unsupported_center, 1, 0, 2, 3)
                .expect("unsupported center should not error structurally")
                .is_none()
        );

        let degree_two_center = bonded_triple(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            BondOrder::Single,
            false,
            false,
        );
        assert!(
            get_uff_inversion_params(&degree_two_center, 0, 1, 2, 0)
                .expect("degree-two center should not error structurally")
                .is_none()
        );

        let non_sp2_nitrogen = trigonal_center_molecule(
            AtomSpec::new(Element::N).with_hybridization(Hybridization::Sp3),
            [
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        assert!(
            get_uff_inversion_params(&non_sp2_nitrogen, 1, 0, 2, 3)
                .expect("non-sp2 nitrogen should not error structurally")
                .is_none()
        );

        let missing_bond = molecule_with_atomic_numbers(&[6, 6, 6, 6]);
        assert!(
            get_uff_inversion_params(&missing_bond, 0, 1, 2, 3)
                .expect("missing-bond case should not error structurally")
                .is_none()
        );
    }

    #[test]
    fn uff_public_api_get_uff_inversion_params_rejects_invalid_atom_index() {
        let mol = trigonal_center_molecule(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            [
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
                AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );

        let err = get_uff_inversion_params(&mol, 1, 0, 2, 4)
            .expect_err("out-of-range atom index should fail");

        assert_eq!(
            err,
            UffPublicApiError::InvalidAtomIndex {
                atom_index: 4,
                atoms: 4
            }
        );
    }

    #[test]
    fn uff_public_api_get_uff_vdw_params_returns_expected_values_for_known_atoms() {
        let mol = bonded_pair(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::O).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            false,
        );

        let params = get_uff_vdw_params(&mol, 0, 1)
            .expect("UFF vdw params should compute")
            .expect("typed atoms should return params");

        let collection = ParamCollection::get_params("").expect("params");
        let c3 = collection.get("C_3").expect("C_3 params");
        let o3 = collection.get("O_3").expect("O_3 params");
        assert_eq!(
            params,
            UffVdw {
                x_ij: calc_nonbonded_minimum(c3, o3),
                d_ij: calc_nonbonded_depth(c3, o3)
            }
        );
    }

    #[test]
    fn uff_public_api_get_uff_vdw_params_does_not_require_bond_between_atoms() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3));
        builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3));
        let mol = builder
            .build()
            .expect("disconnected typed atoms should build");

        let params = get_uff_vdw_params(&mol, 0, 1)
            .expect("disconnected vdw lookup should work")
            .expect("typed disconnected atoms should return params");

        let collection = ParamCollection::get_params("").expect("params");
        let c3 = collection.get("C_3").expect("C_3 params");
        assert_eq!(
            params,
            UffVdw {
                x_ij: calc_nonbonded_minimum(c3, c3),
                d_ij: calc_nonbonded_depth(c3, c3)
            }
        );
    }

    #[test]
    fn uff_public_api_get_uff_vdw_params_returns_none_for_missing_atom_params() {
        let mol = bonded_pair(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::S),
            AtomSpec::new(Element::F).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            false,
        );

        let params = get_uff_vdw_params(&mol, 0, 1).expect("missing atom params should not error");

        assert!(params.is_none());
    }

    #[test]
    fn uff_public_api_get_uff_vdw_params_rejects_invalid_atom_index() {
        let mol = molecule_with_atomic_numbers(&[6, 8]);

        let err = get_uff_vdw_params(&mol, 0, 2).expect_err("out-of-range atom index should fail");

        assert_eq!(
            err,
            UffPublicApiError::InvalidAtomIndex {
                atom_index: 2,
                atoms: 2
            }
        );
    }

    #[test]
    fn uff_public_api_uff_has_all_molecule_params_returns_true_when_all_atoms_are_typed() {
        let mol = bonded_pair(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::O).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            false,
        );

        let found_all =
            uff_has_all_molecule_params(&mol).expect("UFF parameter coverage should compute");

        assert!(found_all);
    }

    #[test]
    fn uff_public_api_uff_has_all_molecule_params_returns_false_for_missing_atom_params() {
        let mol = bonded_pair(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::S),
            AtomSpec::new(Element::F).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            false,
        );

        let found_all =
            uff_has_all_molecule_params(&mol).expect("UFF parameter coverage should compute");

        assert!(!found_all);
    }

    #[test]
    fn uff_public_api_uff_has_all_molecule_params_preserves_empty_molecule_true_boundary() {
        let mol = MoleculeBuilder::new()
            .build()
            .expect("empty molecule should build");

        let found_all =
            uff_has_all_molecule_params(&mol).expect("empty UFF parameter coverage should compute");

        assert!(found_all);
    }

    #[test]
    fn uff_public_api_uff_optimize_molecule_returns_value_style_result_for_typed_molecule() {
        let mol = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        );
        let original_coords = mol.conformers_3d()[0].coordinates().to_vec();

        let result =
            uff_optimize_molecule(&mol, 25, 10.0, -1, true).expect("UFF optimize should run");

        assert!(matches!(result.needs_more, 0 | 1));
        assert!(result.energy.is_finite());
        assert_eq!(
            mol.conformers_3d()[0].coordinates(),
            original_coords.as_slice()
        );
        assert_ne!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_coords.as_slice()
        );
    }

    #[test]
    fn uff_public_api_uff_optimize_molecule_allows_missing_atom_parameters() {
        let mol = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::S),
            AtomSpec::new(Element::F).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]],
        );

        let result = uff_optimize_molecule(&mol, 5, 10.0, -1, true)
            .expect("RDKit constructForceField ignores foundAll for this wrapper");

        assert_eq!(result.needs_more, 0);
        assert!(result.energy.is_finite());
    }

    #[test]
    fn uff_public_api_uff_optimize_molecule_rejects_invalid_conformer_id() {
        let mol = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        );

        let err = uff_optimize_molecule(&mol, 5, 10.0, -2, true)
            .expect_err("conf_id below -1 should fail");

        assert_eq!(err, UffPublicApiError::InvalidConformerId { conf_id: -2 });
    }

    #[test]
    fn uff_public_api_uff_optimize_molecule_rejects_missing_conformer() {
        let mol = bonded_pair(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            false,
        );

        let err = uff_optimize_molecule(&mol, 5, 10.0, -1, true)
            .expect_err("missing conformer should fail");

        assert_eq!(err, UffPublicApiError::Missing3dConformer { conf_id: -1 });
    }

    #[test]
    fn uff_public_api_uff_optimize_molecule_reports_max_iteration_limit() {
        let mol = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
        );

        let result =
            uff_optimize_molecule(&mol, 0, 10.0, -1, true).expect("UFF optimize should run");

        assert_eq!(result.needs_more, 1);
        assert!(result.energy.is_finite());
    }

    #[test]
    fn uff_public_api_uff_optimize_molecule_updates_selected_named_conformer_only() {
        let mol = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
        );
        let first_coords = mol.conformers_3d()[0].coordinates().to_vec();
        let selected_coords = mol.conformers_3d()[1].coordinates().to_vec();

        let result =
            uff_optimize_molecule(&mol, 25, 10.0, 7, true).expect("UFF optimize should run");

        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            first_coords
        );
        assert_ne!(
            result.molecule.conformers_3d()[1].coordinates(),
            selected_coords.as_slice()
        );
        assert_eq!(
            mol.conformers_3d()[1].coordinates(),
            selected_coords.as_slice()
        );
    }

    #[test]
    fn uff_public_api_uff_optimize_molecule_confs_returns_value_style_results_for_all_conformers() {
        let mol = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.3, 0.0, 0.0]],
        );
        let original_first = mol.conformers_3d()[0].coordinates().to_vec();
        let original_second = mol.conformers_3d()[1].coordinates().to_vec();

        let result = uff_optimize_molecule_confs(&mol, 1, 25, 10.0, true)
            .expect("UFF conformer optimization should run");

        assert_eq!(result.conformer_results.len(), 2);
        assert!(
            result
                .conformer_results
                .iter()
                .all(|entry| matches!(entry.needs_more, 0 | 1) && entry.energy.is_finite())
        );
        assert_eq!(
            mol.conformers_3d()[0].coordinates(),
            original_first.as_slice()
        );
        assert_eq!(
            mol.conformers_3d()[1].coordinates(),
            original_second.as_slice()
        );
        assert_ne!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_first.as_slice()
        );
        assert_ne!(
            result.molecule.conformers_3d()[1].coordinates(),
            original_second.as_slice()
        );
    }

    #[test]
    fn uff_public_api_uff_optimize_molecule_confs_allows_missing_atom_parameters() {
        let mol = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::S),
            AtomSpec::new(Element::F).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [1.8, 0.0, 0.0]],
        );

        let result = uff_optimize_molecule_confs(&mol, 1, 5, 10.0, true)
            .expect("RDKit conformer wrapper ignores foundAll while building force field");

        assert_eq!(result.conformer_results.len(), 2);
        assert!(
            result
                .conformer_results
                .iter()
                .all(|entry| entry.needs_more == 0 && entry.energy.is_finite())
        );
    }

    #[test]
    fn uff_public_api_uff_optimize_molecule_confs_rejects_missing_conformer() {
        let mol = bonded_pair(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            false,
        );

        let err = uff_optimize_molecule_confs(&mol, 1, 5, 10.0, true)
            .expect_err("missing conformer should fail");

        assert_eq!(err, UffPublicApiError::Missing3dConformer { conf_id: -1 });
    }

    #[test]
    fn uff_public_api_uff_optimize_molecule_confs_reports_max_iteration_limit() {
        let mol = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
        );

        let result = uff_optimize_molecule_confs(&mol, 1, 0, 10.0, true)
            .expect("UFF conformer optimization should run");

        assert_eq!(result.conformer_results.len(), 1);
        assert_eq!(result.conformer_results[0].needs_more, 1);
        assert!(result.conformer_results[0].energy.is_finite());
    }

    #[test]
    fn uff_public_api_uff_optimize_molecule_confs_handles_non_positive_thread_request_like_non_threaded_rdkit_build()
     {
        let mol = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.2, 0.0, 0.0]],
        );

        let zero_threads = uff_optimize_molecule_confs(&mol, 0, 5, 10.0, true)
            .expect("zero-thread request should use non-threaded RDKit path");
        let negative_threads = uff_optimize_molecule_confs(&mol, -1, 5, 10.0, true)
            .expect("negative-thread request should use non-threaded RDKit path");

        assert_eq!(zero_threads.conformer_results.len(), 2);
        assert_eq!(negative_threads.conformer_results.len(), 2);
        for (left, right) in zero_threads
            .conformer_results
            .iter()
            .zip(negative_threads.conformer_results.iter())
        {
            assert_eq!(left.needs_more, right.needs_more);
            assert!((left.energy - right.energy).abs() < 1.0e-12);
        }
    }
}
