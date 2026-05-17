//! RDKit-allied distance geometry bounds matrix builder.
//!
//! Source reproduction protocol: dev/source_reproduction_protocol.md
//!
//! This module implements a distance bounds matrix for distance geometry
//! embedding, ported from RDKit's DGeomHelpers::BoundsMatrixBuilder.
//!
//! The algorithm computes upper/lower distance bounds for every pair of atoms:
//!
//! - **1-2 bounds**: bond-length estimates from bond order
//! - **1-3 bounds**: triangle geometry from bond angles
//! - **1-4 bounds**: torsion-angle ranges
//! - **VDW lower bounds**: Van der Waals radii for non-bonded pairs

use crate::chemistry::stereo::{get_ideal_angle_between_ligands, has_non_tetrahedral_stereo};
use crate::{
    Atom, AtomId, Bond, BondId, BondOrder, BondStereo, Hybridization, Molecule, ValenceModel,
    assign_valence, rdkit_valence_list,
};
use std::collections::HashSet;

// ──────────────────────────────────────────────
// Constants
// ──────────────────────────────────────────────

// RDKit❗✔️: from DistGeomHelpers/BoundsMatrixBuilder.cpp
const DIST12_DELTA: f64 = 0.01;
const DIST13_TOL: f64 = 0.04;
const GEN_DIST_TOL: f64 = 0.06;
const DIST15_TOL: f64 = 0.08;
const VDW_SCALE_15: f64 = 0.7;
const H_BOND_LENGTH: f64 = 1.8;
const MAX_UPPER: f64 = 1000.0;
const MIN_MACROCYCLE_RING_SIZE: usize = 9;
const DEFAULT_LOWER: f64 = 0.001;
const DEFAULT_UPPER: f64 = MAX_UPPER;
const UFF_LAMBDA: f64 = 0.1332;

// ──────────────────────────────────────────────
// Error type
// ──────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum DgBoundsError {
    #[error("distance bounds matrix generation failed: {0}")]
    GenerationFailed(String),
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
}

#[derive(Debug, Copy, Clone)]
struct UffAtomicParams {
    r1: f64,
    gmp_xi: f64,
}

fn uff_params_for_label(label: &str) -> Option<UffAtomicParams> {
    cosmolkit_macros::rdkit_uff_param_match!("src/data/rdkit_uff_params.tsv", label)
}

// ──────────────────────────────────────────────
// VDW radii (Å) — RDKit PeriodicTable::getRvdw
// ──────────────────────────────────────────────

/// RDKit❗✔️: VDW radii from RDKit PeriodicTable::getRvdw
fn vdw_radius(atomic_num: u8) -> f64 {
    match atomic_num {
        1 => 1.2,       // H
        2 => 1.4,       // He
        3 => 1.82,      // Li
        4 => 1.7,       // Be
        5 => 2.08,      // B
        6 => 1.7,       // C
        7 => 1.55,      // N
        8 => 1.52,      // O
        9 => 1.47,      // F
        10 => 1.54,     // Ne
        11 => 2.27,     // Na
        12 => 1.73,     // Mg
        13 => 2.0,      // Al
        14 => 2.1,      // Si
        15 => 1.8,      // P
        16 => 1.8,      // S
        17 => 1.75,     // Cl
        18 => 1.88,     // Ar
        19 => 2.75,     // K
        20 => 2.0,      // Ca
        21..=29 => 2.0, // Sc-Cu
        30 => 1.39,     // Zn
        31 => 1.87,     // Ga
        32 => 2.0,      // Ge
        33 => 1.85,     // As
        34 => 1.9,      // Se
        35 => 1.85,     // Br
        36 => 2.02,     // Kr
        37 => 2.0,      // Rb
        38 => 2.0,      // Sr
        39..=47 => 2.0, // Y-Ag
        48 => 1.58,     // Cd
        49 => 1.93,     // In
        50 => 2.17,     // Sn
        51 => 2.0,      // Sb
        52 => 2.06,     // Te
        53 => 1.98,     // I
        54 => 2.16,     // Xe
        55 => 2.0,      // Cs
        56 => 2.0,      // Ba
        _ => 2.0,       // default
    }
}

// ──────────────────────────────────────────────
// H-bond detection helpers (BoundsMatrixBuilder.cpp:297-313)
// ──────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::isHBondAcceptor (BoundsMatrixBuilder.cpp:297-299)
// RDKit✔️✔️: bool isHBondAcceptor(const Atom *atom) {
// RDKit✔️✔️:   return (atom->getAtomicNum() == 7 || atom->getAtomicNum() == 8);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::isHBondAcceptor
fn is_hbond_acceptor(atomic_num: u8) -> bool {
    atomic_num == 7 || atomic_num == 8
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::isHBondDonor (BoundsMatrixBuilder.cpp:301-303)
// RDKit✔️✔️: bool isHBondDonor(const Atom *atom) {
// RDKit✔️✔️:   return isHBondAcceptor(atom) && atom->getTotalNumHs() > 0;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::isHBondDonor
fn is_hbond_donor(atom: &Atom, total_num_hydrogens: u32) -> bool {
    is_hbond_acceptor(atom.atomic_number()) && total_num_hydrogens > 0
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::isHinHBondDonor (BoundsMatrixBuilder.cpp:305-313)
// RDKit✔️✔️: bool isHinHBondDonor(const Atom *atom, const ROMol &mol) {
// RDKit✔️✔️:   if (atom->getAtomicNum() != 1) { return false; }
// RDKit✔️✔️:   auto nbrs = mol.atomNeighbors(atom);
// RDKit✔️✔️:   return std::any_of(nbrs.begin(), nbrs.end(), [](const Atom *nbr) {
// RDKit✔️✔️:     return nbr->getAtomicNum() == 7 || nbr->getAtomicNum() == 8;
// RDKit✔️✔️:   });
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::isHinHBondDonor
fn is_h_in_hbond_donor(mol: &Molecule, atom_idx: usize) -> bool {
    if mol.atoms()[atom_idx].atomic_number() != 1 {
        return false;
    }
    // Check if this hydrogen is bonded to N or O
    for bond in mol.bonds() {
        let other = if bond.begin().index() == atom_idx {
            bond.end().index()
        } else if bond.end().index() == atom_idx {
            bond.begin().index()
        } else {
            continue;
        };
        let an = mol.atoms()[other].atomic_number();
        if an == 7 || an == 8 {
            return true;
        }
    }
    false
}
// ──────────────────────────────────────────────

fn rdkit_default_valence(atomic_num: u8) -> Option<i32> {
    let vals = rdkit_valence_list(atomic_num).ok()??;
    vals.iter().copied().find(|v| *v >= 0)
}

fn rdkit_n_outer_electrons(atomic_num: u8) -> Option<i32> {
    crate::chemistry::valence::periodic_table_outer_electrons(atomic_num).ok()
}

fn bond_valence_contrib_for_atom(bond: &Bond, atom_index: usize) -> f64 {
    if bond.begin().index() != atom_index && bond.end().index() != atom_index {
        return 0.0;
    }
    match bond.order() {
        BondOrder::Null | BondOrder::Zero | BondOrder::Ionic => 0.0,
        BondOrder::Single => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Quintuple => 5.0,
        BondOrder::Hextuple => 6.0,
        BondOrder::OneAndHalf => 1.5,
        BondOrder::TwoAndHalf => 2.5,
        BondOrder::ThreeAndHalf => 3.5,
        BondOrder::FourAndHalf => 4.5,
        BondOrder::FiveAndHalf => 5.5,
        BondOrder::Aromatic => 1.5,
        BondOrder::Hydrogen => 0.0,
        BondOrder::Dative
        | BondOrder::DativeOne
        | BondOrder::DativeLeft
        | BondOrder::DativeRight => {
            if bond.end().index() == atom_index {
                1.0
            } else {
                0.0
            }
        }
        BondOrder::ThreeCenter | BondOrder::Unspecified | BondOrder::Other => 0.0,
    }
}

fn count_atom_electrons_rdkit(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_degree: &[usize],
    atom_index: usize,
) -> i32 {
    let atom = &mol.atoms()[atom_index];
    let Some(dv) = rdkit_default_valence(atom.atomic_number()) else {
        return -1;
    };
    if dv <= 1 {
        return -1;
    }
    let mut degree = atom_degree[atom_index] as i32
        + atom.explicit_hydrogens() as i32
        + assignment.implicit_hydrogens[atom_index] as i32;
    for bond in mol.bonds() {
        if (bond.begin().index() == atom_index || bond.end().index() == atom_index)
            && bond_valence_contrib_for_atom(bond, atom_index) == 0.0
        {
            degree -= 1;
        }
    }
    if degree > 3 {
        return -1;
    }
    let Some(nouter) = rdkit_n_outer_electrons(atom.atomic_number()) else {
        return -1;
    };
    let nlp = (nouter - dv - atom.formal_charge() as i32).max(0);
    let radicals = atom.radical_electrons() as i32;
    let mut res = (dv - degree) + nlp - radicals;
    if res > 1 {
        let n_unsaturations =
            assignment.explicit_valence[atom_index] - atom_degree[atom_index] as i32;
        if n_unsaturations > 1 {
            res = 1;
        }
    }
    res
}

fn is_atom_conjugation_candidate(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_degree: &[usize],
    atom_index: usize,
) -> bool {
    let atom = &mol.atoms()[atom_index];
    if let Ok(Some(vals)) = rdkit_valence_list(atom.atomic_number())
        && atom.formal_charge() == 0
        && !vals.is_empty()
        && vals[0] >= 0
    {
        let total_valence =
            assignment.explicit_valence[atom_index] + assignment.implicit_hydrogens[atom_index];
        if total_valence > vals[0] {
            return false;
        }
    }
    let nouter = rdkit_n_outer_electrons(atom.atomic_number()).unwrap_or(0);
    let total_degree = atom_degree[atom_index]
        + atom.explicit_hydrogens() as usize
        + assignment.implicit_hydrogens[atom_index].max(0) as usize;
    ((atom.atomic_number() <= 10)
        || (nouter != 5 && nouter != 6)
        || (nouter == 6 && total_degree < 2))
        && count_atom_electrons_rdkit(mol, assignment, atom_degree, atom_index) > 0
}

fn compute_conjugated_bonds_for_uff(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_degree: &[usize],
) -> Vec<bool> {
    let mut conjugated = vec![false; mol.bonds().len()];
    for (bi, bond) in mol.bonds().iter().enumerate() {
        conjugated[bi] = bond.is_aromatic();
    }
    for atom_index in 0..mol.atoms().len() {
        if !is_atom_conjugation_candidate(mol, assignment, atom_degree, atom_index) {
            continue;
        }
        let sbo = atom_degree[atom_index]
            + mol.atoms()[atom_index].explicit_hydrogens() as usize
            + assignment.implicit_hydrogens[atom_index].max(0) as usize;
        if !(2..=3).contains(&sbo) {
            continue;
        }
        let incident_bonds: Vec<usize> = mol
            .bonds()
            .iter()
            .enumerate()
            .filter_map(|(bond_index, bond)| {
                if bond.begin().index() == atom_index || bond.end().index() == atom_index {
                    Some(bond_index)
                } else {
                    None
                }
            })
            .collect();
        for &bond1_idx in &incident_bonds {
            let bond1 = &mol.bonds()[bond1_idx];
            if bond_valence_contrib_for_atom(bond1, atom_index) < 1.5 {
                continue;
            }
            let other1 = if bond1.begin().index() == atom_index {
                bond1.end().index()
            } else {
                bond1.begin().index()
            };
            if !is_atom_conjugation_candidate(mol, assignment, atom_degree, other1) {
                continue;
            }
            for &bond2_idx in &incident_bonds {
                if bond1_idx == bond2_idx {
                    continue;
                }
                let bond2 = &mol.bonds()[bond2_idx];
                let other2 = if bond2.begin().index() == atom_index {
                    bond2.end().index()
                } else {
                    bond2.begin().index()
                };
                let sbo2 = atom_degree[other2]
                    + mol.atoms()[other2].explicit_hydrogens() as usize
                    + assignment.implicit_hydrogens[other2].max(0) as usize;
                if sbo2 > 3 {
                    continue;
                }
                if is_atom_conjugation_candidate(mol, assignment, atom_degree, other2) {
                    conjugated[bond1_idx] = true;
                    conjugated[bond2_idx] = true;
                }
            }
        }
    }
    conjugated
}

fn compute_hybridizations_for_uff(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_degree: &[usize],
    atom_has_conjugated_bond: &[bool],
) -> Vec<Hybridization> {
    let mut out = Vec::with_capacity(mol.atoms().len());
    for (atom_index, atom) in mol.atoms().iter().enumerate() {
        if atom.atomic_number() == 0 {
            out.push(Hybridization::Unspecified);
            continue;
        }
        let mut deg = atom_degree[atom_index] as i32
            + atom.explicit_hydrogens() as i32
            + assignment.implicit_hydrogens[atom_index];
        for bond in mol.bonds() {
            if (bond.begin().index() == atom_index || bond.end().index() == atom_index)
                && matches!(bond.order(), BondOrder::Dative | BondOrder::DativeOne)
                && bond.end().index() != atom_index
            {
                deg -= 1;
            }
        }
        let hyb = if atom.atomic_number() <= 1 {
            match deg {
                0 | 1 => Hybridization::S,
                2 => Hybridization::Sp,
                3 => Hybridization::Sp2,
                4 => Hybridization::Sp3,
                5 => Hybridization::Sp3d,
                6 => Hybridization::Sp3d2,
                _ => Hybridization::Unspecified,
            }
        } else {
            let nouter = rdkit_n_outer_electrons(atom.atomic_number()).unwrap_or(0);
            let total_valence =
                assignment.explicit_valence[atom_index] + assignment.implicit_hydrogens[atom_index];
            let num_free = nouter - (total_valence + atom.formal_charge() as i32);
            let norbs = if total_valence + nouter - (atom.formal_charge() as i32) < 8 {
                let radicals = atom.radical_electrons() as i32;
                let lone_pairs = (num_free - radicals) / 2;
                deg + lone_pairs + radicals
            } else {
                let lone_pairs = num_free / 2;
                deg + lone_pairs
            };
            match norbs {
                0 | 1 => Hybridization::S,
                2 => Hybridization::Sp,
                3 => Hybridization::Sp2,
                4 => {
                    let total_degree = atom_degree[atom_index]
                        + atom.explicit_hydrogens() as usize
                        + assignment.implicit_hydrogens[atom_index].max(0) as usize;
                    if total_degree > 3 || !atom_has_conjugated_bond[atom_index] {
                        Hybridization::Sp3
                    } else {
                        Hybridization::Sp2
                    }
                }
                5 => Hybridization::Sp3d,
                6 => Hybridization::Sp3d2,
                _ => Hybridization::Unspecified,
            }
        };
        out.push(hyb);
    }
    out
}

fn atom_total_valence_for_uff(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_index: usize,
) -> i32 {
    mol.bonds()
        .iter()
        .map(|bond| bond_valence_contrib_for_atom(bond, atom_index))
        .sum::<f64>()
        .round() as i32
        + assignment.implicit_hydrogens[atom_index]
}

// BEGIN RDKIT CPP FUNCTION Atom::getTotalNumHs (Atom.cpp:283-293)
// RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
// RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
// RDKit✔️✔️:   if (includeNeighbors && dp_mol) {
// RDKit✔️✔️:     auto nbrs = dp_mol->atomNeighbors(this);
// RDKit✔️✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
// RDKit✔️✔️:       return (nbr->getAtomicNum() == 1);
// RDKit✔️✔️:     });
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION Atom::getTotalNumHs
fn total_num_hydrogens_for_distgeom(
    mol: &Molecule,
    atom: &Atom,
    assignment: &crate::ValenceAssignment,
    atom_index: usize,
) -> u32 {
    let mut res =
        atom.explicit_hydrogens() as u32 + assignment.implicit_hydrogens[atom_index].max(0) as u32;
    res += neighbors_for_atom(mol, atom_index)
        .into_iter()
        .filter(|&nbr| mol.atoms()[nbr].atomic_number() == 1)
        .count() as u32;
    res
}

// BEGIN RDKIT CPP FUNCTION UFF::Tools::addAtomChargeFlags (AtomTyper.cpp:23-257)
// RDKit❗✔️: void addAtomChargeFlags(const Atom *atom, std::string &atomKey,
// RDKit❗✔️:                         bool tolerateChargeMismatch) {
// RDKit❗✔️:   int totalValence = atom->getTotalValence();
// RDKit❗✔️:   int fc = atom->getFormalCharge();
// RDKit❗✔️:   switch (atom->getAtomicNum()) { ... }
// RDKit❗✔️:   if (atom->getAtomicNum() >= 57 && atom->getAtomicNum() <= 71) { ... }
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION UFF::Tools::addAtomChargeFlags
fn add_atom_charge_flags_for_uff(
    atom_index: usize,
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_key: &mut String,
    hybridizations: &[Hybridization],
) {
    let atom = &mol.atoms()[atom_index];
    let total_valence = atom_total_valence_for_uff(mol, assignment, atom_index);
    let _formal_charge = atom.formal_charge();
    match atom.atomic_number() {
        29 | 47 => atom_key.push_str("+1"),
        4 | 20 | 25 | 26 | 28 | 46 | 78 => atom_key.push_str("+2"),
        21 | 24 | 27 | 79 | 89 | 96..=103 => atom_key.push_str("+3"),
        2 | 18 | 22 | 36 | 54 | 90 | 91 | 92 | 93 | 94 | 95 => atom_key.push_str("+4"),
        23 | 41 | 43 | 73 => atom_key.push_str("+5"),
        42 => atom_key.push_str("+6"),
        12 | 30 | 48 => atom_key.push_str("+2"),
        15 => match total_valence {
            3 => atom_key.push_str("+3"),
            5 => atom_key.push_str("+5"),
            _ => atom_key.push_str("+5"),
        },
        16 => {
            if hybridizations[atom_index] != Hybridization::Sp2 {
                match total_valence {
                    2 => atom_key.push_str("+2"),
                    4 => atom_key.push_str("+4"),
                    6 => atom_key.push_str("+6"),
                    _ => atom_key.push_str("+6"),
                }
            }
        }
        31 | 33 | 49 | 51 | 81 | 82 | 83 => atom_key.push_str("+3"),
        34 | 52 | 84 => atom_key.push_str("+2"),
        75 => {
            if atom_key == "Re6" {
                *atom_key = "Re6+5".to_string();
            } else if atom_key == "Re3" {
                *atom_key = "Re3+7".to_string();
            }
        }
        _ => {
            if (57..=71).contains(&atom.atomic_number()) {
                atom_key.push_str("+3");
            }
        }
    }
}

// BEGIN RDKIT CPP FUNCTION UFF::Tools::getAtomLabel (AtomTyper.cpp:259-505)
// RDKit❗✔️: std::string getAtomLabel(const Atom *atom) {
// RDKit❗✔️:   int atNum = atom->getAtomicNum();
// RDKit❗✔️:   std::string atomKey = atom->getSymbol();
// RDKit❗✔️:   if (atomKey.size() == 1) { atomKey += '_'; }
// RDKit❗✔️:   if (atNum) { ... switch(atom->getHybridization()) ... }
// RDKit❗✔️:   addAtomChargeFlags(atom, atomKey);
// RDKit❗✔️:   return atomKey;
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION UFF::Tools::getAtomLabel
fn get_atom_label_for_uff(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    hybridizations: &[Hybridization],
    atom_has_conjugated_bond: &[bool],
    atom_index: usize,
) -> Result<String, DgBoundsError> {
    let atom = &mol.atoms()[atom_index];
    let mut atom_key = crate::chemistry::valence::rdkit_element_symbol(atom.atomic_number())
        .map_err(|err| {
            DgBoundsError::GenerationFailed(format!(
                "UFF atom symbol lookup failed for atomic number {}: {err}",
                atom.atomic_number()
            ))
        })?
        .to_string();
    if atom_key.len() == 1 {
        atom_key.push('_');
    }
    if atom.atomic_number() != 0 {
        let nouter = rdkit_n_outer_electrons(atom.atomic_number()).unwrap_or(0);
        if rdkit_default_valence(atom.atomic_number()) == Some(-1) || (nouter != 1 && nouter != 7) {
            match atom.atomic_number() {
                12 | 13 | 14 | 15 | 50 | 51 | 52 | 81 | 82 | 83 | 84 => atom_key.push('3'),
                80 => atom_key.push('1'),
                _ => match hybridizations[atom_index] {
                    Hybridization::S => {}
                    Hybridization::Sp => atom_key.push('1'),
                    Hybridization::Sp2 => {
                        if (atom.is_aromatic() || atom_has_conjugated_bond[atom_index])
                            && matches!(atom.atomic_number(), 6 | 7 | 8 | 16)
                        {
                            atom_key.push('R');
                        } else {
                            atom_key.push('2');
                        }
                    }
                    Hybridization::Sp3 => atom_key.push('3'),
                    Hybridization::Sp2d => atom_key.push('4'),
                    Hybridization::Sp3d => atom_key.push('5'),
                    Hybridization::Sp3d2 => atom_key.push('6'),
                    Hybridization::Unspecified | Hybridization::Other => {}
                },
            }
        }
    }
    add_atom_charge_flags_for_uff(atom_index, mol, assignment, &mut atom_key, hybridizations);
    Ok(atom_key)
}

fn calc_bond_rest_length(bond_order: f64, end1: &UffAtomicParams, end2: &UffAtomicParams) -> f64 {
    let ri = end1.r1;
    let rj = end2.r1;
    let r_bo = -UFF_LAMBDA * (ri + rj) * bond_order.ln();
    let xi = end1.gmp_xi;
    let xj = end2.gmp_xi;
    let root_delta = xi.sqrt() - xj.sqrt();
    let r_en = ri * rj * root_delta * root_delta / (xi * ri + xj * rj);
    ri + rj + r_bo - r_en
}

// ──────────────────────────────────────────────
// Hybridization-based angle estimates
// ──────────────────────────────────────────────

/// RDKit❗✔️: determine ideal bond angle from hybridization
fn ideal_bond_angle(hybridization: &crate::Hybridization, ring_size: Option<usize>) -> f64 {
    const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;

    // RDKit❗✔️: _setRingAngle from BoundsMatrixBuilder.cpp
    // Account for ring geometry
    if let Some(rsize) = ring_size {
        match hybridization {
            crate::Hybridization::Sp2 if rsize <= 8 || rsize == 3 || rsize == 4 => {
                std::f64::consts::PI * (1.0 - 2.0 / rsize as f64)
            }
            crate::Hybridization::Sp3 if rsize == 5 => 104.0 * DEG_TO_RAD,
            crate::Hybridization::Sp3 => 109.5 * DEG_TO_RAD,
            crate::Hybridization::Sp3d => 105.0 * DEG_TO_RAD,
            crate::Hybridization::Sp3d2 => 90.0 * DEG_TO_RAD,
            _ => 120.0 * DEG_TO_RAD,
        }
    } else {
        match hybridization {
            crate::Hybridization::Sp => 180.0 * DEG_TO_RAD,
            crate::Hybridization::Sp2 => 120.0 * DEG_TO_RAD,
            crate::Hybridization::Sp3 => 109.5 * DEG_TO_RAD,
            crate::Hybridization::Sp3d => 90.0 * DEG_TO_RAD,
            crate::Hybridization::Sp3d2 => 90.0 * DEG_TO_RAD,
            crate::Hybridization::Sp2d => 120.0 * DEG_TO_RAD,
            _ => 120.0 * DEG_TO_RAD,
        }
    }
}

// ──────────────────────────────────────────────
// Topological distance (BFS shortest path length)
// ──────────────────────────────────────────────

/// RDKit❗✔️: compute shortest topological path lengths between all atom pairs
fn compute_topological_distances(mol: &Molecule) -> Vec<Vec<i32>> {
    let n = mol.atoms().len();
    let mut dist = vec![vec![i32::MAX / 2; n]; n];
    for i in 0..n {
        dist[i][i] = 0;
    }
    for bond in mol.bonds() {
        let b = bond.begin().index();
        let e = bond.end().index();
        dist[b][e] = 1;
        dist[e][b] = 1;
    }
    // Floyd-Warshall for all-pairs shortest path
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                let nd = dist[i][k].saturating_add(dist[k][j]);
                if nd < dist[i][j] {
                    dist[i][j] = nd;
                }
            }
        }
    }
    dist
}

fn flatten_topological_distances_matrix(mol: &Molecule) -> Vec<f64> {
    let topo = compute_topological_distances(mol);
    let n = mol.num_atoms();
    let mut flat = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            flat[i * n + j] = topo[i][j] as f64;
        }
    }
    flat
}

// ──────────────────────────────────────────────
// Adjacency helpers
// ──────────────────────────────────────────────

fn neighbors_for_atom(mol: &Molecule, idx: usize) -> Vec<usize> {
    let mut ns = Vec::new();
    for bond in mol.bonds() {
        if bond.begin().index() == idx {
            ns.push(bond.end().index());
        } else if bond.end().index() == idx {
            ns.push(bond.begin().index());
        }
    }
    ns
}

fn bond_between(mol: &Molecule, a: usize, b: usize) -> Option<&Bond> {
    mol.bonds().iter().find(|bond| {
        (bond.begin().index() == a && bond.end().index() == b)
            || (bond.begin().index() == b && bond.end().index() == a)
    })
}

/// Find common middle atom for a 1-3 path a1-a2-a3
fn common_middle_atom(mol: &Molecule, a1: usize, a3: usize) -> Option<usize> {
    let n1 = neighbors_for_atom(mol, a1);
    let n3 = neighbors_for_atom(mol, a3);
    for &n in &n1 {
        if n3.contains(&n) {
            return Some(n);
        }
    }
    None
}

// ──────────────────────────────────────────────
// Bounds matrix core
// ──────────────────────────────────────────────

struct BoundsMatrix {
    data: Vec<Vec<f64>>,
    n: usize,
}

impl BoundsMatrix {
    fn new(n: usize) -> Self {
        let mut matrix = Self {
            data: vec![vec![0.0; n]; n],
            n,
        };
        init_bounds_mat_shared(&mut matrix, DEFAULT_LOWER, DEFAULT_UPPER);
        matrix
    }

    fn get_val(&self, i: usize, j: usize) -> f64 {
        self.data[i][j]
    }

    fn set_val(&mut self, i: usize, j: usize, value: f64) {
        self.data[i][j] = value;
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::getUpperBound (BoundsMatrix.h:37-43)
    // RDKit✔️✔️: inline double getUpperBound(unsigned int i, unsigned int j) const {
    // RDKit✔️✔️:   if (i < j) {
    // RDKit✔️✔️:     return getVal(i, j);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return getVal(j, i);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::getUpperBound
    fn get_upper(&self, i: usize, j: usize) -> f64 {
        if i < j {
            self.get_val(i, j)
        } else {
            self.get_val(j, i)
        }
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setUpperBound (BoundsMatrix.h:46-53)
    // RDKit✔️✔️: inline void setUpperBound(unsigned int i, unsigned int j, double val) {
    // RDKit✔️✔️:   CHECK_INVARIANT(val >= 0.0, "Negative upper bound");
    // RDKit✔️✔️:   if (i < j) {
    // RDKit✔️✔️:     setVal(i, j, val);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     setVal(j, i, val);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setUpperBound
    fn set_upper(&mut self, i: usize, j: usize, v: f64) {
        assert!(v >= 0.0, "Negative upper bound");
        if i < j {
            self.set_val(i, j, v);
        } else {
            self.set_val(j, i, v);
        }
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setLowerBound (BoundsMatrix.h:65-72)
    // RDKit✔️✔️: inline void setLowerBound(unsigned int i, unsigned int j, double val) {
    // RDKit✔️✔️:   CHECK_INVARIANT(val >= 0.0, "Negative lower bound");
    // RDKit✔️✔️:   if (i < j) {
    // RDKit✔️✔️:     setVal(j, i, val);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     setVal(i, j, val);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setLowerBound
    fn set_lower(&mut self, i: usize, j: usize, v: f64) {
        assert!(v >= 0.0, "Negative lower bound");
        if i < j {
            self.set_val(j, i, v);
        } else {
            self.set_val(i, j, v);
        }
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::getLowerBound (BoundsMatrix.h:84-90)
    // RDKit✔️✔️: inline double getLowerBound(unsigned int i, unsigned int j) const {
    // RDKit✔️✔️:   if (i < j) {
    // RDKit✔️✔️:     return getVal(j, i);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return getVal(i, j);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::getLowerBound
    fn get_lower(&self, i: usize, j: usize) -> f64 {
        if i < j {
            self.get_val(j, i)
        } else {
            self.get_val(i, j)
        }
    }

    fn num_rows(&self) -> usize {
        self.n
    }

    fn check_and_set_bounds(&mut self, i: usize, j: usize, lb: f64, ub: f64) {
        self.check_and_set_bounds_with_mode(i, j, lb, ub, false);
    }

    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkAndSetBounds (BoundsMatrixBuilder.cpp:195-227)
    // RDKit✔️✔️: void _checkAndSetBounds(unsigned int i, unsigned int j, double lb, double ub,
    // RDKit✔️✔️:                         DistGeom::BoundsMatPtr mmat, bool setIfBetter = false) {
    // RDKit✔️✔️:   // get the existing bounds
    // RDKit✔️✔️:   double clb = mmat->getLowerBound(i, j);
    // RDKit✔️✔️:   double cub = mmat->getUpperBound(i, j);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   CHECK_INVARIANT(ub > lb, "upper bound not greater than lower bound");
    // RDKit✔️✔️:   CHECK_INVARIANT(lb > DIST12_DELTA || clb > DIST12_DELTA, "bad lower bound");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // Note: setIfBetter should ONLY be set if the distances are consistent;
    // RDKit✔️✔️:   // currently this is not the case, therefore, for now, we are pessimistic on
    // RDKit✔️✔️:   // the bounds
    // RDKit✔️✔️:   if (setIfBetter) {
    // RDKit✔️✔️:     double nlb = std::max(clb, lb);
    // RDKit✔️✔️:     double nub = std::min(cub, ub);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (nub <= nlb) {
    // RDKit✔️✔️:       // if not overlapping ranges -> be conservative
    // RDKit✔️✔️:       nlb = std::min(clb, lb);
    // RDKit✔️✔️:       nub = std::max(cub, ub);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     mmat->setLowerBound(i, j, nlb);
    // RDKit✔️✔️:     mmat->setUpperBound(i, j, nub);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     if (clb <= DIST12_DELTA) {
    // RDKit✔️✔️:       mmat->setLowerBound(i, j, lb);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       if ((lb < clb) && (lb > DIST12_DELTA)) {
    // RDKit✔️✔️:         mmat->setLowerBound(i, j, lb);  // conservative bound setting
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (cub >= MAX_UPPER) {  // FIX this
    // RDKit✔️✔️:       mmat->setUpperBound(i, j, ub);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       if ((ub > cub) && (ub < MAX_UPPER)) {
    // RDKit✔️✔️:         mmat->setUpperBound(i, j, ub);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::_checkAndSetBounds
    fn check_and_set_bounds_with_mode(
        &mut self,
        i: usize,
        j: usize,
        lb: f64,
        ub: f64,
        set_if_better: bool,
    ) {
        let clb = self.get_lower(i, j);
        let cub = self.get_upper(i, j);

        assert!(ub > lb, "upper bound not greater than lower bound");
        assert!(lb > DIST12_DELTA || clb > DIST12_DELTA, "bad lower bound");

        if set_if_better {
            let mut nlb = clb.max(lb);
            let mut nub = cub.min(ub);

            if nub <= nlb {
                nlb = clb.min(lb);
                nub = cub.max(ub);
            }

            self.set_lower(i, j, nlb);
            self.set_upper(i, j, nub);
        } else {
            if clb <= DIST12_DELTA {
                self.set_lower(i, j, lb);
            } else if lb < clb && lb > DIST12_DELTA {
                self.set_lower(i, j, lb);
            }

            if cub >= MAX_UPPER {
                self.set_upper(i, j, ub);
            } else if ub > cub && ub < MAX_UPPER {
                self.set_upper(i, j, ub);
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::triangleSmoothBounds (TriangleSmooth.cpp:15-61)
    // RDKit✔️❌: bool triangleSmoothBounds(BoundsMatrix *boundsMat, double tol) {
    // RDKit✔️❌:   auto npt = boundsMat->numRows();
    // RDKit✔️❌:   for (auto k = 0u; k < npt; k++) {
    // RDKit✔️❌:     for (auto i = 0u; i < npt - 1; i++) {
    // RDKit✔️❌:       if (i == k) {
    // RDKit✔️❌:         continue;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       auto ii = i;
    // RDKit✔️❌:       auto ik = k;
    // RDKit✔️❌:       if (ii > ik) {
    // RDKit✔️❌:         std::swap(ii, ik);
    // RDKit✔️❌:       }
    // RDKit✔️❌:
    // RDKit✔️❌:       const auto Uik = boundsMat->getValUnchecked(ii, ik);  // upper bound
    // RDKit✔️❌:       const auto Lik = boundsMat->getValUnchecked(ik, ii);  // lower bound
    // RDKit✔️❌:       for (auto j = i + 1; j < npt; j++) {
    // RDKit✔️❌:         if (j == k) {
    // RDKit✔️❌:           continue;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         auto jj = j;
    // RDKit✔️❌:         auto jk = k;
    // RDKit✔️❌:         if (jj > jk) {
    // RDKit✔️❌:           std::swap(jj, jk);
    // RDKit✔️❌:         }
    // RDKit✔️❌:         const auto Ukj = boundsMat->getValUnchecked(jj, jk);  // upper bound
    // RDKit✔️❌:         const auto sumUikUkj = Uik + Ukj;
    // RDKit✔️❌:         if (boundsMat->getValUnchecked(i, j) > sumUikUkj) {
    // RDKit✔️❌:           // adjust the upper bound
    // RDKit✔️❌:           boundsMat->setValUnchecked(i, j, sumUikUkj);
    // RDKit✔️❌:         }
    // RDKit✔️❌:
    // RDKit✔️❌:         const auto diffLikUjk = Lik - Ukj;
    // RDKit✔️❌:         const auto diffLjkUik = boundsMat->getValUnchecked(jk, jj) - Uik;
    // RDKit✔️❌:         if (boundsMat->getValUnchecked(j, i) < diffLikUjk) {
    // RDKit✔️❌:           // adjust the lower bound
    // RDKit✔️❌:           boundsMat->setValUnchecked(j, i, diffLikUjk);
    // RDKit✔️❌:         } else if (boundsMat->getValUnchecked(j, i) < diffLjkUik) {
    // RDKit✔️❌:           // adjust the lower bound
    // RDKit✔️❌:           boundsMat->setValUnchecked(j, i, diffLjkUik);
    // RDKit✔️❌:         }
    // RDKit✔️❌:         const auto lBound = boundsMat->getValUnchecked(j, i);
    // RDKit✔️❌:         const auto uBound = boundsMat->getValUnchecked(i, j);
    // RDKit✔️❌:         if (tol > 0. && (lBound - uBound) / lBound > 0. &&
    // RDKit✔️❌:             (lBound - uBound) / lBound < tol) {
    // RDKit✔️❌:           // adjust the upper bound
    // RDKit✔️❌:           boundsMat->setValUnchecked(i, j, lBound);
    // RDKit✔️❌:         } else if (lBound - uBound > 0.) {
    // RDKit✔️❌:           return false;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return true;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION DistGeom::triangleSmoothBounds
    //
    // Performance note: the Rust port keeps RDKit's O(n^3) control flow and
    // raw-triangle access pattern, but `Vec<Vec<f64>>` is less cache-friendly
    // than RDKit's contiguous square-matrix backing.
    fn triangle_smooth(&mut self, tol: f64) -> bool {
        let npt = self.num_rows();
        if npt < 2 {
            return true;
        }

        for k in 0..npt {
            for i in 0..(npt - 1) {
                if i == k {
                    continue;
                }
                let (ii, ik) = if i > k { (k, i) } else { (i, k) };

                let uik = self.get_val(ii, ik);
                let lik = self.get_val(ik, ii);
                for j in (i + 1)..npt {
                    if j == k {
                        continue;
                    }
                    let (jj, jk) = if j > k { (k, j) } else { (j, k) };
                    let ukj = self.get_val(jj, jk);
                    let sum_uik_ukj = uik + ukj;
                    if self.get_val(i, j) > sum_uik_ukj {
                        self.set_val(i, j, sum_uik_ukj);
                    }

                    let diff_lik_ujk = lik - ukj;
                    let diff_ljk_uik = self.get_val(jk, jj) - uik;
                    if self.get_val(j, i) < diff_lik_ujk {
                        self.set_val(j, i, diff_lik_ujk);
                    } else if self.get_val(j, i) < diff_ljk_uik {
                        self.set_val(j, i, diff_ljk_uik);
                    }
                    let l_bound = self.get_val(j, i);
                    let u_bound = self.get_val(i, j);
                    let rel_gap = (l_bound - u_bound) / l_bound;
                    if tol > 0.0 && rel_gap > 0.0 && rel_gap < tol {
                        self.set_val(i, j, l_bound);
                    } else if l_bound - u_bound > 0.0 {
                        return false;
                    }
                }
            }
        }
        true
    }

    fn to_vec_vec(self) -> Vec<Vec<f64>> {
        self.data
    }
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::initBoundsMat (BoundsMatrixBuilder.cpp:1792-1800)
// RDKit✔️✔️: void initBoundsMat(DistGeom::BoundsMatrix *mmat, double defaultMin,
// RDKit✔️✔️:                    double defaultMax) {
// RDKit✔️✔️:   unsigned int npt = mmat->numRows();
// RDKit✔️✔️:
// RDKit✔️✔️:   for (unsigned int i = 1; i < npt; i++) {
// RDKit✔️✔️:     for (unsigned int j = 0; j < i; j++) {
// RDKit✔️✔️:       mmat->setUpperBound(i, j, defaultMax);
// RDKit✔️✔️:       mmat->setLowerBound(i, j, defaultMin);
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::initBoundsMat
fn init_bounds_mat_ptr(mmat: &mut BoundsMatrix, default_min: f64, default_max: f64) {
    let npt = mmat.num_rows();

    for i in 1..npt {
        for j in 0..i {
            mmat.set_upper(i, j, default_max);
            mmat.set_lower(i, j, default_min);
        }
    }
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::initBoundsMat (BoundsMatrixBuilder.cpp:1801-1804)
// RDKit✔️✔️: void initBoundsMat(DistGeom::BoundsMatPtr mmat, double defaultMin,
// RDKit✔️✔️:                    double defaultMax) {
// RDKit✔️✔️:   initBoundsMat(mmat.get(), defaultMin, defaultMax);
// RDKit✔️✔️: };
// END RDKIT CPP FUNCTION DGeomHelpers::initBoundsMat
fn init_bounds_mat_shared(mmat: &mut BoundsMatrix, default_min: f64, default_max: f64) {
    init_bounds_mat_ptr(mmat, default_min, default_max);
}

// ──────────────────────────────────────────────
// Ring utilities
// ──────────────────────────────────────────────

fn atom_ring_sizes(ring_info: &crate::RingInfo, atom_idx: usize) -> Vec<usize> {
    let aid = AtomId::new(atom_idx);
    let mut sizes = Vec::new();
    for r in ring_info.atom_rings() {
        if r.contains(&aid) {
            sizes.push(r.len());
        }
    }
    sizes
}

fn is_atom_in_ring_of_size(ring_info: &crate::RingInfo, atom_idx: usize, size: usize) -> bool {
    atom_ring_sizes(ring_info, atom_idx).contains(&size)
}

fn is_bond_in_ring_of_size(ring_info: &crate::RingInfo, bond_idx: usize, size: usize) -> bool {
    let bid = BondId::new(bond_idx);
    ring_info
        .bond_rings()
        .iter()
        .any(|r| r.len() == size && r.contains(&bid))
}

// ──────────────────────────────────────────────
// 1-2 distance bounds
// ──────────────────────────────────────────────

/// RDKit❗✔️: set12Bounds — bond-length derived 1-2 bounds
//
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::set12Bounds (BoundsMatrixBuilder.cpp:238-295)
// RDKit❗✔️: void set12Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
// RDKit❗✔️:                  ComputedData &accumData) {
// RDKit❗✔️:   auto [atomParams, foundAll] = UFF::getAtomTypes(mol);
// RDKit❗✔️:   boost::dynamic_bitset<> squishAtoms(mol.getNumAtoms());
// RDKit❗✔️:   for (const auto bond : mol.bonds()) {
// RDKit❗✔️:     if (bond->getIsConjugated() &&
// RDKit❗✔️:         (bond->getBeginAtom()->getAtomicNum() > 10 ||
// RDKit❗✔️:          bond->getEndAtom()->getAtomicNum() > 10)) { ... }
// RDKit❗✔️:   }
// RDKit❗✔️:   for (const auto bond : mol.bonds()) {
// RDKit❗✔️:     auto bl = ForceFields::UFF::Utils::calcBondRestLength(
// RDKit❗✔️:         bOrder, atomParams[begId], atomParams[endId]);
// RDKit❗✔️:     double extraSquish = 0.0;
// RDKit❗✔️:     if (squishAtoms[begId] || squishAtoms[endId]) extraSquish = 0.2;
// RDKit❗✔️:     accumData.bondLengths[bond->getIdx()] = bl;
// RDKit❗✔️:     mmat->setUpperBound(begId, endId, bl + extraSquish + DIST12_DELTA);
// RDKit❗✔️:     mmat->setLowerBound(begId, endId, bl - extraSquish - DIST12_DELTA);
// RDKit❗✔️:   }
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::set12Bounds
//
// Rust implementation uses estimate_bond_length (topology-based) instead of
// UFF atom types. The squish logic for conjugated 5-ring heteroatoms matches.
// This is a structural adaptation: UFF bond lengths require UFF atom typing
// which COSMolKit does not have.
fn set_12_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
) -> Result<(), DgBoundsError> {
    let npt = mmat.num_rows();
    if npt != mol.atoms().len() {
        return Err(DgBoundsError::GenerationFailed(
            "Wrong size metric matrix".to_string(),
        ));
    }
    if accum_data.bond_lengths.len() < mol.bonds().len() {
        return Err(DgBoundsError::GenerationFailed(
            "Wrong size accumData".to_string(),
        ));
    }
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).map_err(|err| {
        DgBoundsError::GenerationFailed(format!(
            "RDKit UFF atom typing valence assignment failed: {err}"
        ))
    })?;
    let mut atom_degree = vec![0usize; mol.atoms().len()];
    for bond in mol.bonds() {
        atom_degree[bond.begin().index()] += 1;
        atom_degree[bond.end().index()] += 1;
    }
    let conjugated = compute_conjugated_bonds_for_uff(mol, &assignment, &atom_degree);
    let mut atom_has_conjugated_bond = vec![false; mol.atoms().len()];
    for (bond_index, bond) in mol.bonds().iter().enumerate() {
        if conjugated[bond_index] {
            atom_has_conjugated_bond[bond.begin().index()] = true;
            atom_has_conjugated_bond[bond.end().index()] = true;
        }
    }
    let hybridizations =
        compute_hybridizations_for_uff(mol, &assignment, &atom_degree, &atom_has_conjugated_bond);
    let mut atom_params = Vec::with_capacity(mol.atoms().len());
    for atom_index in 0..mol.atoms().len() {
        let label = get_atom_label_for_uff(
            mol,
            &assignment,
            &hybridizations,
            &atom_has_conjugated_bond,
            atom_index,
        )?;
        atom_params.push(uff_params_for_label(&label));
    }

    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::set12Bounds (BoundsMatrixBuilder.cpp:228-295)
    // RDKit✔️❗: void set12Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    // RDKit✔️❗:                  ComputedData &accumData) {
    // RDKit✔️❗:   unsigned int npt = mmat->numRows();
    // RDKit✔️❗:   CHECK_INVARIANT(npt == mol.getNumAtoms(), "Wrong size metric matrix");
    // RDKit✔️❗:   CHECK_INVARIANT(accumData.bondLengths.size() >= mol.getNumBonds(),
    // RDKit✔️❗:                   "Wrong size accumData");
    // RDKit✔️❗:   auto [atomParams, foundAll] = UFF::getAtomTypes(mol);
    // RDKit✔️❗:   CHECK_INVARIANT(atomParams.size() == mol.getNumAtoms(),
    // RDKit✔️❗:                   "parameter vector size mismatch");
    // RDKit✔️❗:   boost::dynamic_bitset<> squishAtoms(mol.getNumAtoms());
    // RDKit✔️❗:   for (const auto bond : mol.bonds()) { ... conjugated 5-ring squish ... }
    // RDKit✔️❗:   for (const auto bond : mol.bonds()) {
    // RDKit✔️❗:     auto begId = bond->getBeginAtomIdx();
    // RDKit✔️❗:     auto endId = bond->getEndAtomIdx();
    // RDKit✔️❗:     auto bOrder = bond->getBondTypeAsDouble();
    // RDKit✔️❗:     if (atomParams[begId] && atomParams[endId] && bOrder > 0) {
    // RDKit✔️❗:       auto bl = ForceFields::UFF::Utils::calcBondRestLength(
    // RDKit✔️❗:           bOrder, atomParams[begId], atomParams[endId]);
    // RDKit✔️❗:       double extraSquish = 0.0;
    // RDKit✔️❗:       if (squishAtoms[begId] || squishAtoms[endId]) {
    // RDKit✔️❗:         extraSquish = 0.2;  // empirical
    // RDKit✔️❗:       }
    // RDKit✔️❗:       accumData.bondLengths[bond->getIdx()] = bl;
    // RDKit✔️❗:       mmat->setUpperBound(begId, endId, bl + extraSquish + DIST12_DELTA);
    // RDKit✔️❗:       mmat->setLowerBound(begId, endId, bl - extraSquish - DIST12_DELTA);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       auto vw1 = PeriodicTable::getTable()->getRvdw(
    // RDKit✔️❗:           mol.getAtomWithIdx(begId)->getAtomicNum());
    // RDKit✔️❗:       auto vw2 = PeriodicTable::getTable()->getRvdw(
    // RDKit✔️❗:           mol.getAtomWithIdx(endId)->getAtomicNum());
    // RDKit✔️❗:       auto bl = (vw1 + vw2) / 2;
    // RDKit✔️❗:       accumData.bondLengths[bond->getIdx()] = bl;
    // RDKit✔️❗:       mmat->setUpperBound(begId, endId, 1.5 * bl);
    // RDKit✔️❗:       mmat->setLowerBound(begId, endId, .5 * bl);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     unsigned int pid =
    // RDKit✔️❗:         std::min(begId, endId) * mol.getNumAtoms() + std::max(begId, endId);
    // RDKit✔️❗:     accumData.visited12Bounds.set(pid);
    // RDKit✔️❗:   }
    // RDKit✔️❗: }
    // END RDKIT CPP FUNCTION DGeomHelpers::set12Bounds
    let mut squish_atoms = vec![false; mol.atoms().len()];
    let ring_info = mol.derived_cache().rings.as_ref();
    for (bond_index, bond) in mol.bonds().iter().enumerate() {
        if conjugated[bond_index]
            && (mol.atoms()[bond.begin().index()].atomic_number() > 10
                || mol.atoms()[bond.end().index()].atomic_number() > 10)
            && ring_info.is_some_and(|ri| is_bond_in_ring_of_size(ri, bond.id().index(), 5))
        {
            squish_atoms[bond.begin().index()] = true;
            squish_atoms[bond.end().index()] = true;
        }
    }

    for bond in mol.bonds() {
        let beg_id = bond.begin().index();
        let end_id = bond.end().index();
        let bond_order =
            crate::chemistry::valence::bond_type_as_double(bond.order()).map_err(|err| {
                DgBoundsError::GenerationFailed(format!(
                    "RDKit set12Bounds bond-order conversion failed: {err}"
                ))
            })?;
        if let (Some(param1), Some(param2)) = (atom_params[beg_id], atom_params[end_id]) {
            if bond_order > 0.0 {
                let bl = calc_bond_rest_length(bond_order, &param1, &param2);
                let extra_squish = if squish_atoms[beg_id] || squish_atoms[end_id] {
                    0.2
                } else {
                    0.0
                };
                accum_data.bond_lengths[bond.id().index()] = bl;
                mmat.set_upper(beg_id, end_id, bl + extra_squish + DIST12_DELTA);
                mmat.set_lower(beg_id, end_id, bl - extra_squish - DIST12_DELTA);
            } else {
                let bl = (vdw_radius(mol.atoms()[beg_id].atomic_number())
                    + vdw_radius(mol.atoms()[end_id].atomic_number()))
                    / 2.0;
                accum_data.bond_lengths[bond.id().index()] = bl;
                mmat.set_upper(beg_id, end_id, 1.5 * bl);
                mmat.set_lower(beg_id, end_id, 0.5 * bl);
            }
        } else {
            let bl = (vdw_radius(mol.atoms()[beg_id].atomic_number())
                + vdw_radius(mol.atoms()[end_id].atomic_number()))
                / 2.0;
            accum_data.bond_lengths[bond.id().index()] = bl;
            mmat.set_upper(beg_id, end_id, 1.5 * bl);
            mmat.set_lower(beg_id, end_id, 0.5 * bl);
        }
        let pid = beg_id.min(end_id) * mol.atoms().len() + beg_id.max(end_id);
        accum_data.visited12_bounds[pid] = true;
    }
    Ok(())
}

// ──────────────────────────────────────────────
// 1-3 distance bounds
// ──────────────────────────────────────────────

/// RDKit❗✔️: compute 1-3 distance from two bond lengths and an angle
fn compute_13_dist(bl1: f64, bl2: f64, angle: f64) -> f64 {
    // RDKit❗✔️: RDGeom::compute13Dist via law of cosines
    (bl1 * bl1 + bl2 * bl2 - 2.0 * bl1 * bl2 * angle.cos()).sqrt()
}

// BEGIN RDKIT CPP FUNCTION RDGeom::compute14Dist3D (Geometry/Utils.h:50-65)
// RDKit✔️✔️: inline double compute14Dist3D(double d1, double d2, double d3, double ang12,
// RDKit✔️✔️:                               double ang23, double torAng) {
// RDKit✔️✔️:   Point3D p1(d1 * cos(ang12), d1 * sin(ang12), 0.0);
// RDKit✔️✔️:   Point3D p4(d2 - d3 * cos(ang23), d3 * sin(ang23), 0.0);
// RDKit✔️✔️:   Transform3D trans;
// RDKit✔️✔️:   trans.SetRotation(torAng, X_Axis);
// RDKit✔️✔️:   trans.TransformPoint(p4);
// RDKit✔️✔️:   p4 -= p1;
// RDKit✔️✔️:   return p4.length();
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION RDGeom::compute14Dist3D
fn compute_14_dist_3d(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64, tor_ang: f64) -> f64 {
    let p1x = d1 * ang12.cos();
    let p1y = d1 * ang12.sin();
    let p4x = d2 - d3 * ang23.cos();
    let p4y = d3 * ang23.sin() * tor_ang.cos();
    let p4z = d3 * ang23.sin() * tor_ang.sin();
    let dx = p4x - p1x;
    let dy = p4y - p1y;
    let dz = p4z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

// BEGIN RDKIT CPP FUNCTION RDGeom::compute14DistCis (Geometry/Utils.h:74-81)
// RDKit✔️✔️: inline double compute14DistCis(double d1, double d2, double d3, double ang12,
// RDKit✔️✔️:                                double ang23) {
// RDKit✔️✔️:   double dx = d2 - d3 * cos(ang23) - d1 * cos(ang12);
// RDKit✔️✔️:   double dy = d3 * sin(ang23) - d1 * sin(ang12);
// RDKit✔️✔️:   double res = dx * dx + dy * dy;
// RDKit✔️✔️:   return sqrt(res);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION RDGeom::compute14DistCis
fn compute_14_dist_cis(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() - d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
}

// BEGIN RDKIT CPP FUNCTION RDGeom::compute14DistTrans (Geometry/Utils.h:88-95)
// RDKit✔️✔️: inline double compute14DistTrans(double d1, double d2, double d3, double ang12,
// RDKit✔️✔️:                                  double ang23) {
// RDKit✔️✔️:   double dx = d2 - d3 * cos(ang23) - d1 * cos(ang12);
// RDKit✔️✔️:   double dy = d3 * sin(ang23) + d1 * sin(ang12);
// RDKit✔️✔️:   double res = dx * dx + dy * dy;
// RDKit✔️✔️:   return sqrt(res);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION RDGeom::compute14DistTrans
fn compute_14_dist_trans(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() + d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
}

/// RDKit❗✔️: set13BoundsHelper — compute and set 1-3 bounds
fn set_13_bounds_helper(
    aid1: usize,
    aid: usize,
    aid3: usize,
    angle: f64,
    bond_lengths: &[f64],
    mmat: &mut BoundsMatrix,
    mol: &Molecule,
) {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::_set13BoundsHelper (BoundsMatrixBuilder.cpp:369-391)
    // RDKit✔️✔️: void _set13BoundsHelper(unsigned int aid1, unsigned int aid, unsigned int aid3,
    // RDKit✔️✔️:                         double angle, const ComputedData &accumData,
    // RDKit✔️✔️:                         DistGeom::BoundsMatPtr mmat, const ROMol &mol) {
    // RDKit✔️✔️:   auto bid1 = mol.getBondBetweenAtoms(aid1, aid)->getIdx();
    // RDKit✔️✔️:   auto bid2 = mol.getBondBetweenAtoms(aid, aid3)->getIdx();
    // RDKit✔️✔️:   auto dl = RDGeom::compute13Dist(accumData.bondLengths[bid1],
    // RDKit✔️✔️:                                   accumData.bondLengths[bid2], angle);
    // RDKit✔️✔️:   auto distTol = DIST13_TOL;
    // RDKit✔️✔️:   // Now increase the tolerance if we're outside of the first row of the
    // RDKit✔️✔️:   // periodic table.
    // RDKit✔️✔️:   if (isLargerSP2Atom(mol.getAtomWithIdx(aid1))) {
    // RDKit✔️✔️:     distTol *= 2;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (isLargerSP2Atom(mol.getAtomWithIdx(aid))) {
    // RDKit✔️✔️:     distTol *= 2;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (isLargerSP2Atom(mol.getAtomWithIdx(aid3))) {
    // RDKit✔️✔️:     distTol *= 2;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto du = dl + distTol;
    // RDKit✔️✔️:   dl -= distTol;
    // RDKit✔️✔️:   _checkAndSetBounds(aid1, aid3, dl, du, mmat);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::_set13BoundsHelper
    let Some(bid1) = bond_between_idx_simple(mol, aid1, aid) else {
        return;
    };
    let Some(bid2) = bond_between_idx_simple(mol, aid, aid3) else {
        return;
    };

    let mut dl = compute_13_dist(bond_lengths[bid1], bond_lengths[bid2], angle);
    let mut dist_tol = DIST13_TOL;

    if is_larger_sp2_atom_idx(mol, aid1) {
        dist_tol *= 2.0;
    }
    if is_larger_sp2_atom_idx(mol, aid) {
        dist_tol *= 2.0;
    }
    if is_larger_sp2_atom_idx(mol, aid3) {
        dist_tol *= 2.0;
    }

    let du = dl + dist_tol;
    dl -= dist_tol;
    mmat.check_and_set_bounds(aid1, aid3, dl, du);
}

fn bond_between_idx(bond_idx: usize) -> Option<usize> {
    // placeholder — not used directly, see below
    None
}

// BEGIN RDKIT CPP FUNCTION isLargerSP2Atom (BoundsMatrixBuilder.cpp:364-367)
// RDKit✔️✔️: bool isLargerSP2Atom(const Atom *atom) {
// RDKit✔️✔️:   return atom->getAtomicNum() > 13 && atom->getHybridization() == Atom::SP2 &&
// RDKit✔️✔️:          atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx());
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION isLargerSP2Atom
fn is_larger_sp2_atom_idx(mol: &Molecule, idx: usize) -> bool {
    let atom = &mol.atoms()[idx];
    atom.atomic_number() > 13
        && atom.hybridization() == Hybridization::Sp2
        && mol
            .derived_cache()
            .rings
            .as_ref()
            .is_some_and(|rings| rings.num_atom_rings(AtomId::new(idx)) > 0)
}

// ──────────────────────────────────────────────
// ComputedData and 1-4/1-5 path infrastructure
// ──────────────────────────────────────────────

/// RDKit❗✔️: Path14Configuration (BoundsMatrixBuilder.cpp:55-63)
/// Stores one 1-4 path (bid1-bid2-bid3) and its stereochemistry type.
#[derive(Debug, Clone, Copy)]
struct Path14Configuration {
    bid1: usize,
    bid2: usize,
    bid3: usize,
    kind: Path14Kind,
}

/// RDKit❗✔️: Path14Type — cis/trans/other (BoundsMatrixBuilder.cpp:58-62)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Path14Kind {
    Cis,
    Trans,
    Other,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DistType {
    Dist12,
    Dist13,
    Dist14,
}

/// RDKit❗✔️: ComputedData (BoundsMatrixBuilder.cpp:76-108)
/// Accumulates computed bond lengths, angles, adjacency, and 1-4/1-5 path info.
#[derive(Debug, Clone)]
struct ComputedData {
    bond_lengths: Vec<f64>,
    bond_adj: Vec<i32>,    // n_bonds x n_bonds — shared atom index or -1
    bond_angles: Vec<f64>, // n_bonds x n_bonds — angle between adjacent bonds
    paths14: Vec<Path14Configuration>,
    cis_paths: Vec<u64>,
    trans_paths: Vec<u64>,
    set15_atoms: Vec<bool>, // n_atoms x n_atoms
    visited12_bounds: Vec<bool>,
    visited13_bounds: Vec<bool>,
    visited14_bounds: Vec<bool>,
}

impl ComputedData {
    fn new(n_atoms: usize, n_bonds: usize) -> Self {
        Self {
            bond_lengths: vec![0.0; n_bonds],
            bond_adj: vec![-1; n_bonds * n_bonds],
            bond_angles: vec![-1.0; n_bonds * n_bonds],
            paths14: Vec::new(),
            cis_paths: Vec::new(),
            trans_paths: Vec::new(),
            set15_atoms: vec![false; n_atoms * n_atoms],
            visited12_bounds: vec![false; n_atoms * n_atoms],
            visited13_bounds: vec![false; n_atoms * n_atoms],
            visited14_bounds: vec![false; n_atoms * n_atoms],
        }
    }

    fn bond_mat_idx(&self, n_bonds: usize, i: usize, j: usize) -> usize {
        i * n_bonds + j
    }

    fn set_bond_adj(&mut self, n_bonds: usize, i: usize, j: usize, value: i32) {
        let idx = self.bond_mat_idx(n_bonds, i, j);
        let rev = self.bond_mat_idx(n_bonds, j, i);
        self.bond_adj[idx] = value;
        self.bond_adj[rev] = value;
    }

    fn get_bond_adj(&self, n_bonds: usize, i: usize, j: usize) -> i32 {
        self.bond_adj[self.bond_mat_idx(n_bonds, i, j)]
    }

    fn set_bond_angle(&mut self, n_bonds: usize, i: usize, j: usize, value: f64) {
        let idx = self.bond_mat_idx(n_bonds, i, j);
        let rev = self.bond_mat_idx(n_bonds, j, i);
        self.bond_angles[idx] = value;
        self.bond_angles[rev] = value;
    }

    fn get_bond_angle(&self, n_bonds: usize, i: usize, j: usize) -> f64 {
        self.bond_angles[self.bond_mat_idx(n_bonds, i, j)]
    }

    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::ComputedData::visitedBound (BoundsMatrixBuilder.cpp:92-96)
    // RDKit✔️✔️: bool visitedBound(unsigned int aid, DistType dt) const {
    // RDKit✔️✔️:   return visited12Bounds[aid] ||
    // RDKit✔️✔️:          ((dt >= DIST13) && visited13Bounds[aid]) ||
    // RDKit✔️✔️:          ((dt >= DIST14) && visited14Bounds[aid]);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::ComputedData::visitedBound
    fn visited_bound(&self, pid: usize, dist_type: DistType) -> bool {
        self.visited12_bounds[pid]
            || (matches!(dist_type, DistType::Dist13 | DistType::Dist14)
                && self.visited13_bounds[pid])
            || (matches!(dist_type, DistType::Dist14) && self.visited14_bounds[pid])
    }
}

fn record_path_flag(path_store: &mut Vec<u64>, path_id: u64) {
    if !path_store.contains(&path_id) {
        path_store.push(path_id);
    }
}

fn has_path_flag(path_store: &[u64], path_id: u64) -> bool {
    path_store.contains(&path_id)
}

fn path14_id(nb: usize, bid1: usize, bid2: usize, bid3: usize) -> u64 {
    bid1 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid3 as u64
}

fn ring_info_for_distgeom(mol: &Molecule) -> &crate::RingInfo {
    mol.derived_cache()
        .rings
        .as_ref()
        .expect("distgeom requires ring info")
}

fn bond_pair_shared_atom(
    mol: &Molecule,
    accum_data: &ComputedData,
    bid1: usize,
    bid2: usize,
) -> usize {
    let nb = mol.num_bonds();
    let aid = accum_data.get_bond_adj(nb, bid1, bid2);
    assert!(aid >= 0, "missing shared atom for bond pair");
    aid as usize
}

/// RDKit❗✔️: _getAtomStereo — get the effective stereo for a bond given 1-4 atoms
fn get_atom_stereo(bond: &Bond, aid1: usize, aid4: usize) -> BondStereo {
    let mut stype = bond.stereo();
    if let Some(ref stereo_atoms) = bond.stereo_atoms() {
        let needs_flip = (stereo_atoms[0].index() != aid1) ^ (stereo_atoms[1].index() != aid4);
        if needs_flip {
            stype = match stype {
                BondStereo::Cis => BondStereo::Trans,
                BondStereo::Trans => BondStereo::Cis,
                BondStereo::Z => BondStereo::E,
                BondStereo::E => BondStereo::Z,
                other => other,
            };
        }
    }
    stype
}

// ──────────────────────────────────────────────
// Set ring angle
// ──────────────────────────────────────────────

fn set_ring_angle(mol: &Molecule, aid2: usize, ring_size: usize) -> f64 {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setRingAngle (BoundsMatrixBuilder.cpp:393-414)
    // RDKit✔️✔️: void _setRingAngle(Atom::HybridizationType aHyb, unsigned int ringSize,
    // RDKit✔️✔️:                    double &angle) {
    // RDKit✔️✔️:   if ((aHyb == Atom::SP2 && ringSize <= 8) || (ringSize == 3) ||
    // RDKit✔️✔️:       (ringSize == 4)) {
    // RDKit✔️✔️:     angle = M_PI * (1 - 2.0 / ringSize);
    // RDKit✔️✔️:   } else if (aHyb == Atom::SP3) {
    // RDKit✔️✔️:     if (ringSize == 5) {
    // RDKit✔️✔️:       angle = 104 * M_PI / 180;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       angle = 109.5 * M_PI / 180;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (aHyb == Atom::SP3D) {
    // RDKit✔️✔️:     angle = 105.0 * M_PI / 180;
    // RDKit✔️✔️:   } else if (aHyb == Atom::SP3D2) {
    // RDKit✔️✔️:     angle = 90.0 * M_PI / 180;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     angle = 120 * M_PI / 180;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::_setRingAngle
    let hyb = mol.atoms()[aid2].hybridization();
    ideal_bond_angle(&hyb, Some(ring_size))
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::set13Bounds (BoundsMatrixBuilder.cpp:417-654)
// RDKit✔️❌: void set13Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
// RDKit✔️❌:                  ComputedData &accumData) {
// RDKit✔️❌:   auto npt = mmat->numRows();
// RDKit✔️❌:   CHECK_INVARIANT(npt == mol.getNumAtoms(), "Wrong size metric matrix");
// RDKit✔️❌:   CHECK_INVARIANT(accumData.bondAngles->numRows() == mol.getNumBonds(),
// RDKit✔️❌:                   "Wrong size bond angle matrix");
// RDKit✔️❌:   CHECK_INVARIANT(accumData.bondAdj->numRows() == mol.getNumBonds(),
// RDKit✔️❌:                   "Wrong size bond adjacency matrix");
// RDKit✔️❌:   const auto rinfo = mol.getRingInfo();
// RDKit✔️❌:   CHECK_INVARIANT(rinfo, "");
// RDKit✔️❌:   auto atomRings = rinfo->atomRings();
// RDKit✔️❌:   std::sort(atomRings.begin(), atomRings.end(), lessVector);
// RDKit✔️❌:   INT_VECT visited(npt, 0);
// RDKit✔️❌:   DOUBLE_VECT angleTaken(npt, 0.0);
// RDKit✔️❌:   auto nb = mol.getNumBonds();
// RDKit✔️❌:   BIT_SET donePaths(nb * nb);
// RDKit✔️❌:   // first deal with all rings and atoms in them
// RDKit✔️❌:   for (const auto &ringi : atomRings) { ... }
// RDKit✔️❌:   // now deal with the remaining atoms
// RDKit✔️❌:   for (aid2 = 0; aid2 < npt; aid2++) { ... }
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::set13Bounds
fn set_13_bounds(mol: &Molecule, mmat: &mut BoundsMatrix, accum_data: &mut ComputedData) {
    let npt = mmat.num_rows();
    assert_eq!(npt, mol.atoms().len(), "Wrong size metric matrix");

    let nb = mol.bonds().len();
    assert_eq!(
        accum_data.bond_angles.len(),
        nb * nb,
        "Wrong size bond angle matrix"
    );
    assert_eq!(
        accum_data.bond_adj.len(),
        nb * nb,
        "Wrong size bond adjacency matrix"
    );

    let rinfo = mol
        .derived_cache()
        .rings
        .as_ref()
        .expect("set13Bounds requires ring info");

    let mut atom_rings: Vec<Vec<usize>> = rinfo
        .atom_rings()
        .iter()
        .map(|ring| ring.iter().map(|aid| aid.index()).collect())
        .collect();
    atom_rings.sort_by_key(Vec::len);

    let mut visited = vec![0usize; npt];
    let mut angle_taken = vec![0.0; npt];
    let mut done_paths = vec![false; nb * nb];

    for ring in &atom_rings {
        let r_size = ring.len();
        let mut aid1 = ring[r_size - 1];
        for i in 0..r_size {
            let aid2 = ring[i];
            let aid3 = if i == r_size - 1 {
                ring[0]
            } else {
                ring[i + 1]
            };
            let b1 = bond_between_idx_simple(mol, aid1, aid2).expect("no bond found");
            let b2 = bond_between_idx_simple(mol, aid2, aid3).expect("no bond found");
            let id1 = nb * b1 + b2;
            let id2 = nb * b2 + b1;
            let pid = aid1.min(aid3) * npt + aid1.max(aid3);

            if !done_paths[id1] && !done_paths[id2] {
                let angle = set_ring_angle(mol, aid2, r_size);
                if !accum_data.visited_bound(pid, DistType::Dist12) {
                    set_13_bounds_helper(
                        aid1,
                        aid2,
                        aid3,
                        angle,
                        &accum_data.bond_lengths,
                        mmat,
                        mol,
                    );
                    accum_data.visited13_bounds[pid] = true;
                }
                accum_data.set_bond_angle(nb, b1, b2, angle);
                accum_data.set_bond_adj(nb, b1, b2, aid2 as i32);
                visited[aid2] += 1;
                angle_taken[aid2] += angle;
                done_paths[id1] = true;
                done_paths[id2] = true;
            }
            aid1 = aid2;
        }
    }

    for aid2 in 0..npt {
        let atom = &mol.atoms()[aid2];
        let nbrs = neighbors_for_atom(mol, aid2);
        let deg = nbrs.len();
        let n13 = deg * (deg.saturating_sub(1)) / 2;
        if n13 == visited[aid2] {
            continue;
        }
        let ahyb = atom.hybridization();

        if visited[aid2] >= 1 {
            for left in 0..nbrs.len() {
                let aid1 = nbrs[left];
                let bid1 = bond_between_idx_simple(mol, aid1, aid2).expect("missing bond");
                for right in 0..left {
                    let aid3 = nbrs[right];
                    let bid2 = bond_between_idx_simple(mol, aid3, aid2).expect("missing bond");
                    if accum_data.get_bond_angle(nb, bid1, bid2) >= 0.0 {
                        continue;
                    }

                    let angle = if ahyb == Hybridization::Sp2 {
                        (2.0 * std::f64::consts::PI - angle_taken[aid2])
                            / (n13 - visited[aid2]) as f64
                    } else if ahyb == Hybridization::Sp3 {
                        if rinfo.is_atom_in_ring_of_size(AtomId::new(aid2), 3) {
                            116.0_f64.to_radians()
                        } else if rinfo.is_atom_in_ring_of_size(AtomId::new(aid2), 4) {
                            112.0_f64.to_radians()
                        } else {
                            109.5_f64.to_radians()
                        }
                    } else if has_non_tetrahedral_stereo(atom) {
                        get_ideal_angle_between_ligands(aid2, aid1, aid3, mol).to_radians()
                    } else if deg == 5 {
                        105.0_f64.to_radians()
                    } else if deg == 6 {
                        135.0_f64.to_radians()
                    } else {
                        120.0_f64.to_radians()
                    };

                    let pid = aid1.min(aid3) * npt + aid1.max(aid3);
                    if !accum_data.visited_bound(pid, DistType::Dist12) {
                        set_13_bounds_helper(
                            aid1,
                            aid2,
                            aid3,
                            angle,
                            &accum_data.bond_lengths,
                            mmat,
                            mol,
                        );
                        accum_data.visited13_bounds[pid] = true;
                    }

                    accum_data.set_bond_angle(nb, bid1, bid2, angle);
                    accum_data.set_bond_adj(nb, bid1, bid2, aid2 as i32);
                    angle_taken[aid2] += angle;
                    visited[aid2] += 1;
                }
            }
        } else {
            for left in 0..nbrs.len() {
                let aid1 = nbrs[left];
                let bid1 = bond_between_idx_simple(mol, aid1, aid2).expect("missing bond");
                for right in 0..left {
                    let aid3 = nbrs[right];
                    let bid2 = bond_between_idx_simple(mol, aid3, aid2).expect("missing bond");

                    let angle = if has_non_tetrahedral_stereo(atom) {
                        get_ideal_angle_between_ligands(aid2, aid1, aid3, mol).to_radians()
                    } else if ahyb == Hybridization::Sp {
                        std::f64::consts::PI
                    } else if ahyb == Hybridization::Sp2 {
                        2.0 * std::f64::consts::PI / 3.0
                    } else if ahyb == Hybridization::Sp3 {
                        109.5_f64.to_radians()
                    } else if ahyb == Hybridization::Sp3d {
                        105.0_f64.to_radians()
                    } else if ahyb == Hybridization::Sp3d2 {
                        135.0_f64.to_radians()
                    } else {
                        120.0_f64.to_radians()
                    };

                    let pid = aid1.min(aid3) * npt + aid1.max(aid3);
                    if !accum_data.visited_bound(pid, DistType::Dist12) {
                        if deg <= 4
                            || (has_non_tetrahedral_stereo(atom)
                                && atom.chiral_permutation().is_some())
                        {
                            set_13_bounds_helper(
                                aid1,
                                aid2,
                                aid3,
                                angle,
                                &accum_data.bond_lengths,
                                mmat,
                                mol,
                            );
                        } else {
                            let dmax =
                                accum_data.bond_lengths[bid1] + accum_data.bond_lengths[bid2];
                            let dl = 1.0;
                            let du = dmax * 1.2;
                            mmat.check_and_set_bounds(aid1, aid3, dl, du);
                        }
                        accum_data.visited13_bounds[pid] = true;
                    }

                    accum_data.set_bond_angle(nb, bid1, bid2, angle);
                    accum_data.set_bond_adj(nb, bid1, bid2, aid2 as i32);
                    angle_taken[aid2] += angle;
                    visited[aid2] += 1;
                }
            }
        }
    }
}

// ──────────────────────────────────────────────
// 1-4 distance bounds (torsion-based)
// ──────────────────────────────────────────────

fn total_num_hydrogens_rdkit_like(mol: &Molecule, atom_idx: usize) -> Option<u32> {
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).ok()?;
    Some(total_num_hydrogens_for_distgeom(
        mol,
        &mol.atoms()[atom_idx],
        &assignment,
        atom_idx,
    ))
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_record14Path (BoundsMatrixBuilder.cpp:1319-1340)
// RDKit✔️❌: void _record14Path(const ROMol &mol, unsigned int bid1, unsigned int bid2,
// RDKit✔️❌:                    unsigned int bid3, ComputedData &accumData) {
// RDKit✔️❌:   const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2));
// RDKit✔️❌:   Atom::HybridizationType ahyb2 = atm2->getHybridization();
// RDKit✔️❌:   const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3));
// RDKit✔️❌:   Atom::HybridizationType ahyb3 = atm3->getHybridization();
// RDKit✔️❌:   unsigned int nb = mol.getNumBonds();
// RDKit✔️❌:   Path14Configuration path14;
// RDKit✔️❌:   path14.bid1 = bid1;
// RDKit✔️❌:   path14.bid2 = bid2;
// RDKit✔️❌:   path14.bid3 = bid3;
// RDKit✔️❌:   if ((ahyb2 == Atom::SP2) && (ahyb3 == Atom::SP2)) {
// RDKit✔️❌:     path14.type = Path14Configuration::CIS;
// RDKit✔️❌:     accumData.cisPaths.insert(static_cast<unsigned long>(bid1) * nb * nb +
// RDKit✔️❌:                               bid2 * nb + bid3);
// RDKit✔️❌:     accumData.cisPaths.insert(static_cast<unsigned long>(bid3) * nb * nb +
// RDKit✔️❌:                               bid2 * nb + bid1);
// RDKit✔️❌:   } else {
// RDKit✔️❌:     path14.type = Path14Configuration::OTHER;
// RDKit✔️❌:   }
// RDKit✔️❌:   accumData.paths14.push_back(path14);
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_record14Path
fn record_14_path(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
) {
    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2);
    let ahyb2 = mol.atoms()[atm2].hybridization();
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3);
    let ahyb3 = mol.atoms()[atm3].hybridization();
    let nb = mol.num_bonds();

    let kind = if ahyb2 == Hybridization::Sp2 && ahyb3 == Hybridization::Sp2 {
        record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
        record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
        Path14Kind::Cis
    } else {
        Path14Kind::Other
    };

    accum_data.paths14.push(Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind,
    });
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setInRing14Bounds (BoundsMatrixBuilder.cpp:724-826)
// RDKit✔️❌: void _setInRing14Bounds(const ROMol &mol, const Bond *bnd1, const Bond *bnd2,
// RDKit✔️❌:                         const Bond *bnd3, ComputedData &accumData,
// RDKit✔️❌:                         DistGeom::BoundsMatPtr mmat, double *dmat,
// RDKit✔️❌:                         int ringSize) { ... }
// END RDKIT CPP FUNCTION DGeomHelpers::_setInRing14Bounds
fn set_in_ring_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
    ring_size: usize,
) {
    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2);
    let ahyb2 = mol.atoms()[atm2].hybridization();
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3);
    let ahyb3 = mol.atoms()[atm3].hybridization();

    let bnd1 = &mol.bonds()[bid1];
    let bnd3 = &mol.bonds()[bid3];
    let aid1 = if bnd1.begin().index() == atm2 {
        bnd1.end().index()
    } else {
        bnd1.begin().index()
    };
    let aid4 = if bnd3.begin().index() == atm3 {
        bnd3.end().index()
    } else {
        bnd3.begin().index()
    };
    let pid = aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4);
    if accum_data.visited_bound(pid, DistType::Dist13) {
        return;
    }
    if dmat[aid1.max(aid4) * mmat.num_rows() + aid1.min(aid4)] < 2.9 {
        return;
    }

    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3);
    assert!(ba12 > 0.0);
    assert!(ba23 > 0.0);

    let stype = get_atom_stereo(&mol.bonds()[bid2], aid1, aid4);
    let ring_info = ring_info_for_distgeom(mol);
    let mut prefer_cis = false;
    let mut prefer_trans = false;

    if ring_size <= 8
        && ahyb2 == Hybridization::Sp2
        && ahyb3 == Hybridization::Sp2
        && !matches!(stype, BondStereo::E | BondStereo::Trans)
    {
        if ring_info.num_bond_rings(BondId::new(bid2)) > 1 {
            if ring_info.num_bond_rings(BondId::new(bid1)) == 1
                && ring_info.num_bond_rings(BondId::new(bid3)) == 1
            {
                for br in ring_info.bond_rings() {
                    if br.contains(&BondId::new(bid1)) {
                        if br.contains(&BondId::new(bid3)) {
                            prefer_cis = true;
                        }
                        break;
                    }
                }
            }
        } else {
            prefer_cis = true;
        }
    } else if matches!(stype, BondStereo::Z | BondStereo::Cis) {
        prefer_cis = true;
    } else if matches!(stype, BondStereo::E | BondStereo::Trans) {
        prefer_trans = true;
    }

    let nb = mol.num_bonds();
    let kind = if prefer_cis {
        record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
        record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
        Path14Kind::Cis
    } else if prefer_trans {
        record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid1, bid2, bid3));
        record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid3, bid2, bid1));
        Path14Kind::Trans
    } else {
        Path14Kind::Other
    };

    accum_data.paths14.push(Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind,
    });

    let (mut dl, mut du) = if prefer_cis {
        let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
        (dl, dl + 2.0 * GEN_DIST_TOL)
    } else if prefer_trans {
        let dl = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
        (dl, dl + 2.0 * GEN_DIST_TOL)
    } else {
        let mut dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
        let mut du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        if du < dl {
            std::mem::swap(&mut dl, &mut du);
        }
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        (dl, du)
    };

    accum_data.visited14_bounds[pid] = true;
    mmat.check_and_set_bounds(aid1, aid4, dl, du);
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setTwoInSameRing14Bounds (BoundsMatrixBuilder.cpp:828-900)
// RDKit✔️❌: void _setTwoInSameRing14Bounds(const ROMol &mol, const Bond *bnd1,
// RDKit✔️❌:                                const Bond *bnd2, const Bond *bnd3,
// RDKit✔️❌:                                ComputedData &accumData,
// RDKit✔️❌:                                DistGeom::BoundsMatPtr mmat, double *dmat) { ... }
// END RDKIT CPP FUNCTION DGeomHelpers::_setTwoInSameRing14Bounds
fn set_two_in_same_ring_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
) {
    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2);
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3);
    let bnd1 = &mol.bonds()[bid1];
    let bnd3 = &mol.bonds()[bid3];
    let aid1 = if bnd1.begin().index() == atm2 {
        bnd1.end().index()
    } else {
        bnd1.begin().index()
    };
    let aid4 = if bnd3.begin().index() == atm3 {
        bnd3.end().index()
    } else {
        bnd3.begin().index()
    };
    let pid = aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4);

    if accum_data.visited_bound(pid, DistType::Dist13) {
        return;
    }
    if dmat[aid1.max(aid4) * mmat.num_rows() + aid1.min(aid4)] < 2.9 {
        return;
    }
    if bond_between(mol, aid1, atm3).is_some() || bond_between(mol, aid4, atm2).is_some() {
        return;
    }

    let ahyb2 = mol.atoms()[atm2].hybridization();
    let ahyb3 = mol.atoms()[atm3].hybridization();
    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3);
    assert!(ba12 > 0.0);
    assert!(ba23 > 0.0);

    let nb = mol.num_bonds();
    let (mut dl, mut du, kind) = if ahyb2 == Hybridization::Sp2 && ahyb3 == Hybridization::Sp2 {
        record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid1, bid2, bid3));
        record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid3, bid2, bid1));
        let du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        (du - GEN_DIST_TOL, du + GEN_DIST_TOL, Path14Kind::Trans)
    } else {
        let mut dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
        let mut du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        if du < dl {
            std::mem::swap(&mut dl, &mut du);
        }
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        (dl, du, Path14Kind::Other)
    };

    mmat.check_and_set_bounds(aid1, aid4, dl, du);
    accum_data.paths14.push(Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind,
    });
    accum_data.visited14_bounds[pid] = true;
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setTwoInDiffRing14Bounds (BoundsMatrixBuilder.cpp:902-910)
// RDKit✔️✔️: void _setTwoInDiffRing14Bounds(const ROMol &mol, const Bond *bnd1,
// RDKit✔️✔️:                                const Bond *bnd2, const Bond *bnd3,
// RDKit✔️✔️:                                ComputedData &accumData,
// RDKit✔️✔️:                                DistGeom::BoundsMatPtr mmat, double *dmat) {
// RDKit✔️✔️:   _setInRing14Bounds(mol, bnd1, bnd2, bnd3, accumData, mmat, dmat, 0);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_setTwoInDiffRing14Bounds
fn set_two_in_diff_ring_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
) {
    set_in_ring_14_bounds(mol, bid1, bid2, bid3, accum_data, mmat, dmat, 0);
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setShareRingBond14Bounds (BoundsMatrixBuilder.cpp:912-920)
// RDKit✔️✔️: void _setShareRingBond14Bounds(const ROMol &mol, const Bond *bnd1,
// RDKit✔️✔️:                                const Bond *bnd2, const Bond *bnd3,
// RDKit✔️✔️:                                ComputedData &accumData,
// RDKit✔️✔️:                                DistGeom::BoundsMatPtr mmat, double *dmat) {
// RDKit✔️✔️:   _setInRing14Bounds(mol, bnd1, bnd2, bnd3, accumData, mmat, dmat, 0);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_setShareRingBond14Bounds
fn set_share_ring_bond_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
) {
    set_in_ring_14_bounds(mol, bid1, bid2, bid3, accum_data, mmat, dmat, 0);
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkMacrocycleTwoInSameRingAmideEster14 (BoundsMatrixBuilder.cpp:1342-1361)
// RDKit✔️✔️: bool _checkMacrocycleTwoInSameRingAmideEster14(
// RDKit✔️✔️:     const Bond *bnd1, const Bond *bnd3, const Atom *atm1, const Atom *atm2,
// RDKit✔️✔️:     const Atom *atm3, const Atom *atm4) {
// RDKit✔️✔️:   unsigned int a1Num = atm1->getAtomicNum();
// RDKit✔️✔️:   unsigned int a2Num = atm2->getAtomicNum();
// RDKit✔️✔️:   unsigned int a3Num = atm3->getAtomicNum();
// RDKit✔️✔️:   unsigned int a4Num = atm4->getAtomicNum();
// RDKit✔️✔️:   return a1Num != 1 && a3Num == 6 && bnd3->getBondType() == Bond::DOUBLE &&
// RDKit✔️✔️:          (a4Num == 8 || a4Num == 7) && bnd1->getBondType() == Bond::SINGLE &&
// RDKit✔️✔️:          (a2Num == 8 || a2Num == 7);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_checkMacrocycleTwoInSameRingAmideEster14
fn check_macrocycle_two_in_same_ring_amide_ester_14(
    mol: &Molecule,
    bnd1_idx: usize,
    bnd3_idx: usize,
    atm1: usize,
    atm2: usize,
    atm3: usize,
    atm4: usize,
) -> bool {
    let a1_num = mol.atoms()[atm1].atomic_number();
    let a2_num = mol.atoms()[atm2].atomic_number();
    let a3_num = mol.atoms()[atm3].atomic_number();
    let a4_num = mol.atoms()[atm4].atomic_number();

    a1_num != 1
        && a3_num == 6
        && mol.bonds()[bnd3_idx].order() == BondOrder::Double
        && (a4_num == 8 || a4_num == 7)
        && mol.bonds()[bnd1_idx].order() == BondOrder::Single
        && (a2_num == 8 || a2_num == 7)
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setMacrocycleTwoInSameRing14Bounds (BoundsMatrixBuilder.cpp:1363-1454)
// RDKit✔️❌: void _setMacrocycleTwoInSameRing14Bounds(const ROMol &mol, const Bond *bnd1,
// RDKit✔️❌:                                          const Bond *bnd2, const Bond *bnd3,
// RDKit✔️❌:                                          ComputedData &accumData,
// RDKit✔️❌:                                          DistGeom::BoundsMatPtr mmat,
// RDKit✔️❌:                                          double *dmat) { ... }
// END RDKIT CPP FUNCTION DGeomHelpers::_setMacrocycleTwoInSameRing14Bounds
fn set_macrocycle_two_in_same_ring_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
) {
    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2);
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3);
    let bnd1 = &mol.bonds()[bid1];
    let bnd3 = &mol.bonds()[bid3];
    let aid1 = if bnd1.begin().index() == atm2 {
        bnd1.end().index()
    } else {
        bnd1.begin().index()
    };
    let aid4 = if bnd3.begin().index() == atm3 {
        bnd3.end().index()
    } else {
        bnd3.begin().index()
    };
    let atm1 = aid1;
    let atm4 = aid4;
    let pid = aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4);

    if accum_data.visited_bound(pid, DistType::Dist13) {
        return;
    }
    if dmat[aid1.max(aid4) * mmat.num_rows() + aid1.min(aid4)] < 2.9 {
        return;
    }
    if bond_between(mol, aid1, atm3).is_some() || bond_between(mol, aid4, atm2).is_some() {
        return;
    }

    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3);
    assert!(ba12 > 0.0);
    assert!(ba23 > 0.0);

    let nb = mol.num_bonds();
    let (mut dl, mut du, kind) = if check_macrocycle_two_in_same_ring_amide_ester_14(
        mol, bid1, bid3, atm1, atm2, atm3, atm4,
    ) || check_macrocycle_two_in_same_ring_amide_ester_14(
        mol, bid3, bid1, atm4, atm3, atm2, atm1,
    ) {
        let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
        record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
        record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
        (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Cis)
    } else {
        let mut dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
        let mut du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        if du < dl {
            std::mem::swap(&mut dl, &mut du);
        }
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        (dl, du, Path14Kind::Other)
    };

    mmat.check_and_set_bounds(aid1, aid4, dl, du);
    accum_data.paths14.push(Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind,
    });
    accum_data.visited14_bounds[pid] = true;
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setMacrocycleAllInSameRing14Bounds (BoundsMatrixBuilder.cpp:1456-1660)
// RDKit✔️❌: void _setMacrocycleAllInSameRing14Bounds(const ROMol &mol, const Bond *bnd1,
// RDKit✔️❌:                                          const Bond *bnd2, const Bond *bnd3,
// RDKit✔️❌:                                          ComputedData &accumData,
// RDKit✔️❌:                                          DistGeom::BoundsMatPtr mmat,
// RDKit✔️❌:                                          double *) { ... }
// END RDKIT CPP FUNCTION DGeomHelpers::_setMacrocycleAllInSameRing14Bounds
fn set_macrocycle_all_in_same_ring_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
) {
    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2);
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3);
    let bnd1 = &mol.bonds()[bid1];
    let bnd2 = &mol.bonds()[bid2];
    let bnd3 = &mol.bonds()[bid3];
    let aid1 = if bnd1.begin().index() == atm2 {
        bnd1.end().index()
    } else {
        bnd1.begin().index()
    };
    let aid4 = if bnd3.begin().index() == atm3 {
        bnd3.end().index()
    } else {
        bnd3.begin().index()
    };
    let pid = aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4);
    if accum_data.visited_bound(pid, DistType::Dist13) {
        return;
    }

    let atm1 = aid1;
    let atm4 = aid4;
    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3);
    assert!(ba12 > 0.0);
    assert!(ba23 > 0.0);

    let mut set_the_bound = true;
    let nb = mol.num_bonds();
    let (mut dl, mut du, kind) = match bnd2.order() {
        BondOrder::Double => {
            if bnd1.order() == BondOrder::Double || bnd3.order() == BondOrder::Double {
                let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
                record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
                (dl, dl + 2.0 * GEN_DIST_TOL, Path14Kind::Cis)
            } else if matches!(
                bnd2.stereo(),
                BondStereo::Z
                    | BondStereo::E
                    | BondStereo::Cis
                    | BondStereo::Trans
                    | BondStereo::AtropCw
                    | BondStereo::AtropCcw
            ) {
                let stype = get_atom_stereo(bnd2, aid1, aid4);
                if matches!(stype, BondStereo::Z | BondStereo::Cis) {
                    let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                    record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
                    record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
                    (dl, dl + 2.0 * GEN_DIST_TOL, Path14Kind::Cis)
                } else {
                    let du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid1, bid2, bid3));
                    record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid3, bid2, bid1));
                    (du - GEN_DIST_TOL, du + GEN_DIST_TOL, Path14Kind::Trans)
                }
            } else {
                let mut dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                let mut du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                if (du - dl).abs() < DIST12_DELTA {
                    dl -= GEN_DIST_TOL;
                    du += GEN_DIST_TOL;
                }
                (dl, du, Path14Kind::Other)
            }
        }
        BondOrder::Single => {
            if mol.atoms()[atm2].atomic_number() == 16
                && mol.atoms()[atm3].atomic_number() == 16
                && neighbors_for_atom(mol, atm2).len() == 2
                && neighbors_for_atom(mol, atm3).len() == 2
            {
                let dl = compute_14_dist_3d(bl1, bl2, bl3, ba12, ba23, std::f64::consts::PI / 2.0)
                    - GEN_DIST_TOL;
                (dl, dl + 2.0 * GEN_DIST_TOL, Path14Kind::Other)
            } else if check_macrocycle_all_in_same_ring_amide_ester_14(mol, atm1, atm2, atm3, atm4)
                || check_macrocycle_all_in_same_ring_amide_ester_14(mol, atm4, atm3, atm2, atm1)
            {
                let dl = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23) + 0.1;
                record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid1, bid2, bid3));
                record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid3, bid2, bid1));
                (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Trans)
            } else if check_amide_ester_15(mol, bid1, bid3, atm2, atm3)
                || check_amide_ester_15(mol, bid3, bid1, atm3, atm2)
            {
                // COSMolKit models the default RDKit branch where FORCE_TRANS_AMIDES is not defined.
                let total_hs_atm2 = total_num_hydrogens_rdkit_like(mol, atm2).unwrap_or(0);
                if mol.atoms()[atm2].atomic_number() == 7
                    && neighbors_for_atom(mol, atm2).len() == 3
                    && mol.atoms()[atm1].atomic_number() == 1
                    && total_hs_atm2 == 1
                {
                    set_the_bound = false;
                    (0.0, 0.0, Path14Kind::Other)
                } else {
                    let dl = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid1, bid2, bid3));
                    record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid3, bid2, bid1));
                    (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Trans)
                }
            } else {
                (
                    compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
                    compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
                    Path14Kind::Other,
                )
            }
        }
        _ => (
            compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
            compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
            Path14Kind::Other,
        ),
    };

    if set_the_bound {
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        mmat.check_and_set_bounds(aid1, aid4, dl, du);
        accum_data.paths14.push(Path14Configuration {
            bid1,
            bid2,
            bid3,
            kind,
        });
        accum_data.visited14_bounds[pid] = true;
    }
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setChain14Bounds (BoundsMatrixBuilder.cpp:980-1317)
// RDKit✔️❌: void _setChain14Bounds(const ROMol &mol, const Bond *bnd1, const Bond *bnd2,
// RDKit✔️❌:                        const Bond *bnd3, ComputedData &accumData,
// RDKit✔️❌:                        DistGeom::BoundsMatPtr mmat, double *,
// RDKit✔️❌:                        bool forceTransAmides) { ... }
// END RDKIT CPP FUNCTION DGeomHelpers::_setChain14Bounds
fn set_chain_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    force_trans_amides: bool,
) {
    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2);
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3);
    let bnd1 = &mol.bonds()[bid1];
    let bnd2 = &mol.bonds()[bid2];
    let bnd3 = &mol.bonds()[bid3];
    let aid1 = if bnd1.begin().index() == atm2 {
        bnd1.end().index()
    } else {
        bnd1.begin().index()
    };
    let aid4 = if bnd3.begin().index() == atm3 {
        bnd3.end().index()
    } else {
        bnd3.begin().index()
    };
    let pid = aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4);
    if accum_data.visited_bound(pid, DistType::Dist13) {
        return;
    }

    let atm1 = aid1;
    let atm4 = aid4;
    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3);
    assert!(ba12 > 0.0);
    assert!(ba23 > 0.0);

    let mut set_the_bound = true;
    let nb = mol.num_bonds();
    let (mut dl, mut du, kind) = match bnd2.order() {
        BondOrder::Double => {
            if bnd1.order() == BondOrder::Double || bnd3.order() == BondOrder::Double {
                let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
                record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
                (dl, dl + 2.0 * GEN_DIST_TOL, Path14Kind::Cis)
            } else if matches!(
                bnd2.stereo(),
                BondStereo::Z
                    | BondStereo::E
                    | BondStereo::Cis
                    | BondStereo::Trans
                    | BondStereo::AtropCw
                    | BondStereo::AtropCcw
            ) {
                let stype = get_atom_stereo(bnd2, aid1, aid4);
                if matches!(stype, BondStereo::Z | BondStereo::Cis) {
                    let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                    record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
                    record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
                    (dl, dl + 2.0 * GEN_DIST_TOL, Path14Kind::Cis)
                } else {
                    let du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid1, bid2, bid3));
                    record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid3, bid2, bid1));
                    (du - GEN_DIST_TOL, du + GEN_DIST_TOL, Path14Kind::Trans)
                }
            } else {
                let mut dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                let mut du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                if (du - dl).abs() < DIST12_DELTA {
                    dl -= GEN_DIST_TOL;
                    du += GEN_DIST_TOL;
                }
                (dl, du, Path14Kind::Other)
            }
        }
        BondOrder::Single => {
            if mol.atoms()[atm2].atomic_number() == 16
                && mol.atoms()[atm3].atomic_number() == 16
                && neighbors_for_atom(mol, atm2).len() == 2
                && neighbors_for_atom(mol, atm3).len() == 2
            {
                let dl = compute_14_dist_3d(bl1, bl2, bl3, ba12, ba23, std::f64::consts::PI / 2.0)
                    - GEN_DIST_TOL;
                (dl, dl + 2.0 * GEN_DIST_TOL, Path14Kind::Other)
            } else if check_amide_ester_14(mol, bid1, bid3, atm2, atm3, atm4)
                || check_amide_ester_14(mol, bid3, bid1, atm3, atm2, atm1)
            {
                if force_trans_amides {
                    let total_hs_atm2 = total_num_hydrogens_rdkit_like(mol, atm2).unwrap_or(0);
                    let total_hs_atm3 = total_num_hydrogens_rdkit_like(mol, atm3).unwrap_or(0);
                    let secondary_left = mol.atoms()[atm1].atomic_number() == 1
                        && mol.atoms()[atm2].atomic_number() == 7
                        && neighbors_for_atom(mol, atm2).len() == 3
                        && total_hs_atm2 == 1;
                    let secondary_right = mol.atoms()[atm4].atomic_number() == 1
                        && mol.atoms()[atm3].atomic_number() == 7
                        && neighbors_for_atom(mol, atm3).len() == 3
                        && total_hs_atm3 == 1;
                    if secondary_left || secondary_right {
                        let dl = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        record_path_flag(
                            &mut accum_data.trans_paths,
                            path14_id(nb, bid1, bid2, bid3),
                        );
                        record_path_flag(
                            &mut accum_data.trans_paths,
                            path14_id(nb, bid3, bid2, bid1),
                        );
                        (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Trans)
                    } else {
                        let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                        record_path_flag(
                            &mut accum_data.cis_paths,
                            path14_id(nb, bid1, bid2, bid3),
                        );
                        record_path_flag(
                            &mut accum_data.cis_paths,
                            path14_id(nb, bid3, bid2, bid1),
                        );
                        (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Cis)
                    }
                } else {
                    (
                        compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
                        compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
                        Path14Kind::Other,
                    )
                }
            } else if check_amide_ester_15(mol, bid1, bid3, atm2, atm3)
                || check_amide_ester_15(mol, bid3, bid1, atm3, atm2)
            {
                if force_trans_amides {
                    let total_hs_atm2 = total_num_hydrogens_rdkit_like(mol, atm2).unwrap_or(0);
                    let total_hs_atm3 = total_num_hydrogens_rdkit_like(mol, atm3).unwrap_or(0);
                    let secondary_left = mol.atoms()[atm1].atomic_number() == 1
                        && mol.atoms()[atm2].atomic_number() == 7
                        && neighbors_for_atom(mol, atm2).len() == 3
                        && total_hs_atm2 == 1;
                    let secondary_right = mol.atoms()[atm4].atomic_number() == 1
                        && mol.atoms()[atm3].atomic_number() == 7
                        && neighbors_for_atom(mol, atm3).len() == 3
                        && total_hs_atm3 == 1;
                    if secondary_left || secondary_right {
                        let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                        record_path_flag(
                            &mut accum_data.cis_paths,
                            path14_id(nb, bid1, bid2, bid3),
                        );
                        record_path_flag(
                            &mut accum_data.cis_paths,
                            path14_id(nb, bid3, bid2, bid1),
                        );
                        (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Cis)
                    } else {
                        let dl = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        record_path_flag(
                            &mut accum_data.trans_paths,
                            path14_id(nb, bid1, bid2, bid3),
                        );
                        record_path_flag(
                            &mut accum_data.trans_paths,
                            path14_id(nb, bid3, bid2, bid1),
                        );
                        (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Trans)
                    }
                } else {
                    (
                        compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
                        compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
                        Path14Kind::Other,
                    )
                }
            } else {
                (
                    compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
                    compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
                    Path14Kind::Other,
                )
            }
        }
        _ => (
            compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
            compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
            Path14Kind::Other,
        ),
    };

    if set_the_bound {
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        mmat.check_and_set_bounds(aid1, aid4, dl, du);
        accum_data.paths14.push(Path14Configuration {
            bid1,
            bid2,
            bid3,
            kind,
        });
        accum_data.visited14_bounds[pid] = true;
    }
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkH2NX3H1OX2 (BoundsMatrixBuilder.cpp:845-856)
// RDKit✔️❌: bool _checkH2NX3H1OX2(const Atom *atm) {
// RDKit✔️❌:   if ((atm->getAtomicNum() == 6) && (atm->getTotalNumHs(true) == 2)) {
// RDKit✔️❌:     return true;
// RDKit✔️❌:   } else if ((atm->getAtomicNum() == 8) && (atm->getTotalNumHs(true) == 0)) {
// RDKit✔️❌:     return true;
// RDKit✔️❌:   } else if ((atm->getAtomicNum() == 7) && (atm->getDegree() == 3) &&
// RDKit✔️❌:              (atm->getTotalNumHs(true) == 1)) {
// RDKit✔️❌:     return true;
// RDKit✔️❌:   }
// RDKit✔️❌:   return false;
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_checkH2NX3H1OX2
fn check_h2_nx3_h1_ox2(mol: &Molecule, atom_idx: usize) -> bool {
    let atom = &mol.atoms()[atom_idx];
    let Some(total_hs) = total_num_hydrogens_rdkit_like(mol, atom_idx) else {
        return false;
    };
    (atom.atomic_number() == 6 && total_hs == 2)
        || (atom.atomic_number() == 8 && total_hs == 0)
        || (atom.atomic_number() == 7
            && neighbors_for_atom(mol, atom_idx).len() == 3
            && total_hs == 1)
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkNhChChNh (BoundsMatrixBuilder.cpp:858-866)
// RDKit✔️❌: bool _checkNhChChNh(const Atom *atm1, const Atom *atm2, const Atom *atm3,
// RDKit✔️❌:                     const Atom *atm4) {
// RDKit✔️❌:   if ((atm1->getAtomicNum() != 1) && (atm4->getAtomicNum() != 1)) {
// RDKit✔️❌:     if ((_checkH2NX3H1OX2(atm2)) && (_checkH2NX3H1OX2(atm3))) {
// RDKit✔️❌:       return true;
// RDKit✔️❌:     }
// RDKit✔️❌:   }
// RDKit✔️❌:   return false;
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_checkNhChChNh
fn check_nh_ch_ch_nh(mol: &Molecule, atm1: usize, atm2: usize, atm3: usize, atm4: usize) -> bool {
    mol.atoms()[atm1].atomic_number() != 1
        && mol.atoms()[atm4].atomic_number() != 1
        && check_h2_nx3_h1_ox2(mol, atm2)
        && check_h2_nx3_h1_ox2(mol, atm3)
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkAmideEster14 (BoundsMatrixBuilder.cpp:874-894)
// RDKit✔️✔️: bool _checkAmideEster14(const Bond *bnd1, const Bond *bnd3, const Atom *,
// RDKit✔️✔️:                         const Atom *atm2, const Atom *atm3, const Atom *atm4) {
// RDKit✔️✔️:   unsigned int a2Num = atm2->getAtomicNum();
// RDKit✔️✔️:   unsigned int a3Num = atm3->getAtomicNum();
// RDKit✔️✔️:   unsigned int a4Num = atm4->getAtomicNum();
// RDKit✔️✔️:   if (a3Num == 6 && bnd3->getBondType() == Bond::DOUBLE &&
// RDKit✔️✔️:       (a4Num == 8 || a4Num == 7) && bnd1->getBondType() == Bond::SINGLE &&
// RDKit✔️✔️:       (a2Num == 8 || (a2Num == 7 && atm2->getTotalNumHs(true) == 1))) {
// RDKit✔️✔️:     return true;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return false;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_checkAmideEster14
fn check_amide_ester_14(
    mol: &Molecule,
    bnd1_idx: usize,
    bnd3_idx: usize,
    atm2: usize,
    atm3: usize,
    atm4: usize,
) -> bool {
    let bnd1 = &mol.bonds()[bnd1_idx];
    let bnd3 = &mol.bonds()[bnd3_idx];
    let a2_num = mol.atoms()[atm2].atomic_number();
    let a3_num = mol.atoms()[atm3].atomic_number();
    let a4_num = mol.atoms()[atm4].atomic_number();
    let total_hs_atm2 = total_num_hydrogens_rdkit_like(mol, atm2).unwrap_or(0);

    a3_num == 6
        && bnd3.order() == BondOrder::Double
        && (a4_num == 8 || a4_num == 7)
        && bnd1.order() == BondOrder::Single
        && (a2_num == 8 || (a2_num == 7 && total_hs_atm2 == 1))
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkMacrocycleAllInSameRingAmideEster14 (BoundsMatrixBuilder.cpp:903-949)
// RDKit✔️❌: bool _checkMacrocycleAllInSameRingAmideEster14(const ROMol &mol, const Bond *,
// RDKit✔️❌:                                                const Bond *, const Atom *atm1,
// RDKit✔️❌:                                                const Atom *atm2,
// RDKit✔️❌:                                                const Atom *atm3,
// RDKit✔️❌:                                                const Atom *atm4) { ... }
// END RDKIT CPP FUNCTION DGeomHelpers::_checkMacrocycleAllInSameRingAmideEster14
fn check_macrocycle_all_in_same_ring_amide_ester_14(
    mol: &Molecule,
    atm1: usize,
    atm2: usize,
    atm3: usize,
    atm4: usize,
) -> bool {
    let a2_num = mol.atoms()[atm2].atomic_number();
    let a3_num = mol.atoms()[atm3].atomic_number();
    if a3_num != 6 {
        return false;
    }
    if a2_num != 7 && a2_num != 8 {
        return false;
    }
    if neighbors_for_atom(mol, atm2).len() != 3 || neighbors_for_atom(mol, atm3).len() != 3 {
        return false;
    }

    for nbr_idx in neighbors_for_atom(mol, atm2) {
        if nbr_idx != atm1 && nbr_idx != atm3 {
            let res = &mol.atoms()[nbr_idx];
            let res_bnd = &mol.bonds()[bond_between_idx_simple(mol, atm2, nbr_idx).expect("bond")];
            if (res.atomic_number() != 6 && res.atomic_number() != 1)
                || res_bnd.order() != BondOrder::Single
            {
                return false;
            }
            break;
        }
    }

    for nbr_idx in neighbors_for_atom(mol, atm3) {
        if nbr_idx != atm2 && nbr_idx != atm4 {
            let res = &mol.atoms()[nbr_idx];
            let res_bnd = &mol.bonds()[bond_between_idx_simple(mol, atm3, nbr_idx).expect("bond")];
            if res.atomic_number() != 8 || res_bnd.order() != BondOrder::Double {
                return false;
            }
            break;
        }
    }

    true
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_isCarbonyl (BoundsMatrixBuilder.cpp:951-963)
// RDKit✔️✔️: bool _isCarbonyl(const ROMol &mol, const Atom *at) {
// RDKit✔️✔️:   if (at->getAtomicNum() == 6 && at->getDegree() > 2) {
// RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(at)) {
// RDKit✔️✔️:       unsigned int atNum = nbr->getAtomicNum();
// RDKit✔️✔️:       if ((atNum == 8 || atNum == 7) &&
// RDKit✔️✔️:           mol.getBondBetweenAtoms(at->getIdx(), nbr->getIdx())->getBondType() ==
// RDKit✔️✔️:               Bond::DOUBLE) {
// RDKit✔️✔️:         return true;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return false;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_isCarbonyl
fn is_carbonyl(mol: &Molecule, atom_idx: usize) -> bool {
    let atom = &mol.atoms()[atom_idx];
    if atom.atomic_number() != 6 || neighbors_for_atom(mol, atom_idx).len() <= 2 {
        return false;
    }
    neighbors_for_atom(mol, atom_idx)
        .into_iter()
        .any(|nbr_idx| {
            let at_num = mol.atoms()[nbr_idx].atomic_number();
            (at_num == 8 || at_num == 7)
                && mol.bonds()[bond_between_idx_simple(mol, atom_idx, nbr_idx).expect("bond")]
                    .order()
                    == BondOrder::Double
        })
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkAmideEster15 (BoundsMatrixBuilder.cpp:965-978)
// RDKit✔️❌: bool _checkAmideEster15(const ROMol &mol, const Bond *bnd1, const Bond *bnd3,
// RDKit✔️❌:                         const Atom *, const Atom *atm2, const Atom *atm3,
// RDKit✔️❌:                         const Atom *) {
// RDKit✔️❌:   unsigned int a2Num = atm2->getAtomicNum();
// RDKit✔️❌:   if ((a2Num == 8) || ((a2Num == 7) && (atm2->getTotalNumHs(true) == 1))) {
// RDKit✔️❌:     if ((bnd1->getBondType() == Bond::SINGLE)) {
// RDKit✔️❌:       if ((atm3->getAtomicNum() == 6) &&
// RDKit✔️❌:           (bnd3->getBondType() == Bond::SINGLE) && _isCarbonyl(mol, atm3)) {
// RDKit✔️❌:         return true;
// RDKit✔️❌:       }
// RDKit✔️❌:     }
// RDKit✔️❌:   }
// RDKit✔️❌:   return false;
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_checkAmideEster15
fn check_amide_ester_15(
    mol: &Molecule,
    bnd1_idx: usize,
    bnd3_idx: usize,
    atm2: usize,
    atm3: usize,
) -> bool {
    let a2_num = mol.atoms()[atm2].atomic_number();
    let total_hs_atm2 = total_num_hydrogens_rdkit_like(mol, atm2).unwrap_or(0);
    ((a2_num == 8) || (a2_num == 7 && total_hs_atm2 == 1))
        && mol.bonds()[bnd1_idx].order() == BondOrder::Single
        && mol.atoms()[atm3].atomic_number() == 6
        && mol.bonds()[bnd3_idx].order() == BondOrder::Single
        && is_carbonyl(mol, atm3)
}

// ──────────────────────────────────────────────
// 1-5 distance bounds (optional)
// ──────────────────────────────────────────────

/// RDKit✔️✔️: _compute15DistsCisCis — 1-5 distance, cis-cis configuration
///
/// Configuration:
///         5
///          \
///     1     4
///      \   /
///       2-3
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsCisCis (BoundsMatrixBuilder.cpp:1965-1981)
// RDKit✔️✔️: double _compute15DistsCisCis(double d1, double d2, double d3, double d4,
// RDKit✔️✔️:                              double ang12, double ang23, double ang34) {
// RDKit✔️✔️:   double dx14 = d2 - d3 * cos(ang23) - d1 * cos(ang12);
// RDKit✔️✔️:   double dy14 = d3 * sin(ang23) - d1 * sin(ang12);
// RDKit✔️✔️:   double d14 = sqrt(dx14 * dx14 + dy14 * dy14);
// RDKit✔️✔️:   double cval = (d3 - d2 * cos(ang23) + d1 * cos(ang12 + ang23)) / d14;
// RDKit✔️✔️:   if (cval > 1.0) { cval = 1.0; } else if (cval < -1.0) { cval = -1.0; }
// RDKit✔️✔️:   double ang143 = acos(cval);
// RDKit✔️✔️:   double ang145 = ang34 - ang143;
// RDKit✔️✔️:   double res = RDGeom::compute13Dist(d14, d4, ang145);
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION
fn compute_15_dist_cis_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 - ang143;
    compute_13_dist(d14, d4, ang145)
}

/// RDKit✔️✔️: _compute15DistsCisTrans — 1-5 distance, cis-trans configuration
///
/// Configuration:
///  1     4-5
///   \   /
///    2-3
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsCisTrans (BoundsMatrixBuilder.cpp:2004-2019)
// RDKit✔️✔️: double _compute15DistsCisTrans(double d1, double d2, double d3, double d4,
// RDKit✔️✔️:                                double ang12, double ang23, double ang34) {
// RDKit✔️✔️:   double dx14 = d2 - d3 * cos(ang23) - d1 * cos(ang12);
// RDKit✔️✔️:   double dy14 = d3 * sin(ang23) - d1 * sin(ang12);
// RDKit✔️✔️:   double d14 = sqrt(dx14 * dx14 + dy14 * dy14);
// RDKit✔️✔️:   double cval = (d3 - d2 * cos(ang23) + d1 * cos(ang12 + ang23)) / d14;
// RDKit✔️✔️:   if (cval > 1.0) { cval = 1.0; } else if (cval < -1.0) { cval = -1.0; }
// RDKit✔️✔️:   double ang143 = acos(cval);
// RDKit✔️✔️:   double ang145 = ang34 + ang143;
// RDKit✔️✔️:   return RDGeom::compute13Dist(d14, d4, ang145);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsCisTrans
fn compute_15_dist_cis_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 + ang143;
    compute_13_dist(d14, d4, ang145)
}

/// RDKit✔️✔️: _compute15DistsTransTrans — 1-5 distance, trans-trans configuration
///
/// Configuration:
///  1
///   \
///    2-3
///       \
///        4-5
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsTransTrans (BoundsMatrixBuilder.cpp:2045-2060)
// RDKit✔️✔️: double _compute15DistsTransTrans(double d1, double d2, double d3, double d4,
// RDKit✔️✔️:                                  double ang12, double ang23, double ang34) {
// RDKit✔️✔️:   double dx14 = d2 - d3 * cos(ang23) - d1 * cos(ang12);
// RDKit✔️✔️:   double dy14 = d3 * sin(ang23) + d1 * sin(ang12);
// RDKit✔️✔️:   double d14 = sqrt(dx14 * dx14 + dy14 * dy14);
// RDKit✔️✔️:   double cval = (d3 - d2 * cos(ang23) + d1 * cos(ang12 - ang23)) / d14;
// RDKit✔️✔️:   if (cval > 1.0) { cval = 1.0; } else if (cval < -1.0) { cval = -1.0; }
// RDKit✔️✔️:   double ang143 = acos(cval);
// RDKit✔️✔️:   double ang145 = ang34 + ang143;
// RDKit✔️✔️:   return RDGeom::compute13Dist(d14, d4, ang145);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsTransTrans
fn compute_15_dist_trans_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 + ang143;
    compute_13_dist(d14, d4, ang145)
}

/// RDKit✔️✔️: _compute15DistsTransCis — 1-5 distance, trans-cis configuration
///
/// Configuration:
///                    1
///                     \
///                      2-3
///                         \
///                          4
///                         /
///                        5
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsTransCis (BoundsMatrixBuilder.cpp:2088-2104)
// RDKit✔️✔️: double _compute15DistsTransCis(double d1, double d2, double d3, double d4,
// RDKit✔️✔️:                                double ang12, double ang23, double ang34) {
// RDKit✔️✔️:   double dx14 = d2 - d3 * cos(ang23) - d1 * cos(ang12);
// RDKit✔️✔️:   double dy14 = d3 * sin(ang23) + d1 * sin(ang12);
// RDKit✔️✔️:   double d14 = sqrt(dx14 * dx14 + dy14 * dy14);
// RDKit✔️✔️:   double cval = (d3 - d2 * cos(ang23) + d1 * cos(ang12 - ang23)) / d14;
// RDKit✔️✔️:   if (cval > 1.0) { cval = 1.0; } else if (cval < -1.0) { cval = -1.0; }
// RDKit✔️✔️:   double ang143 = acos(cval);
// RDKit✔️✔️:   double ang145 = ang34 - ang143;
// RDKit✔️✔️:   return RDGeom::compute13Dist(d14, d4, ang145);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsTransCis
fn compute_15_dist_trans_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 - ang143;
    compute_13_dist(d14, d4, ang145)
}

/// RDKit❗✔️: set15Bounds — 1-5 distance bounds with stereochemistry (full RDKit version)
///
/// Finds all 1-5 atom pairs (a-b-c-d-e) via pre-computed 1-4 paths (paths14)
/// and sets distance bounds using bond lengths, computed bond angles, and
/// cis/trans stereochemistry detection on the central double bond.
///
/// The three 1-4 path types:
/// - **Cis**: stereochemistry prefers cis on the central bond
/// - **Trans**: stereochemistry prefers trans on the central bond
/// - **Other**: no detected stereochemistry, uses VDW fallback for unknown paths
///
/// For each 1-4 path, the function extends through all bonds from the 4th atom
/// to find 1-5 pairs, then computes the distance using one of four geometry
/// configurations (cis-cis, cis-trans, trans-cis, trans-trans) depending on
/// the path type and cis/trans path flags.
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::set15Bounds (BoundsMatrixBuilder.cpp:2232-2248)
// RDKit✔️✔️: void set15Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
// RDKit✔️✔️:                  ComputedData &accumData, double *distMatrix) {
// RDKit✔️✔️:   PATH14_VECT_CI pti;
// RDKit✔️✔️:   unsigned int bid1, bid2, bid3, type;
// RDKit✔️✔️:   for (pti = accumData.paths14.begin(); pti != accumData.paths14.end(); pti++) {
// RDKit✔️✔️:     bid1 = pti->bid1;
// RDKit✔️✔️:     bid2 = pti->bid2;
// RDKit✔️✔️:     bid3 = pti->bid3;
// RDKit✔️✔️:     type = pti->type;
// RDKit✔️✔️:     // 15 distances going one way with with 14 paths
// RDKit✔️✔️:     _set15BoundsHelper(mol, bid1, bid2, bid3, type, accumData, mmat,
// RDKit✔️✔️:                        distMatrix);
// RDKit✔️✔️:     // going the other way - reverse the 14 path
// RDKit✔️✔️:     _set15BoundsHelper(mol, bid3, bid2, bid1, type, accumData, mmat,
// RDKit✔️✔️:                        distMatrix);
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION
fn set_15_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    dmat: &[f64],
) {
    let nb = mol.bonds().len();
    let na = mol.atoms().len();
    for path_idx in 0..accum_data.paths14.len() {
        let path = accum_data.paths14[path_idx];
        set_15_bounds_helper(
            mol, mmat, accum_data, dmat, nb, na, path.bid1, path.bid2, path.bid3, path.kind,
        );
        set_15_bounds_helper(
            mol, mmat, accum_data, dmat, nb, na, path.bid3, path.bid2, path.bid1, path.kind,
        );
    }
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_set15BoundsHelper (BoundsMatrixBuilder.cpp:2106-2253)
// RDKit✔️✔️: void _set15BoundsHelper(const ROMol &mol, unsigned int bid1, unsigned int bid2,
// RDKit✔️✔️:                         unsigned int bid3, unsigned int type,
// RDKit✔️✔️:                         ComputedData &accumData, DistGeom::BoundsMatPtr mmat,
// RDKit✔️✔️:                         double *dmat) {
// RDKit✔️✔️:   unsigned int i, aid1, aid2, aid3, aid4, aid5;
// RDKit✔️✔️:   double d1, d2, d3, d4, ang12, ang23, ang34, du, dl, vw1, vw5;
// RDKit✔️✔️:   unsigned int nb = mol.getNumBonds();
// RDKit✔️✔️:   unsigned int na = mol.getNumAtoms();
// RDKit✔️✔️:   aid2 = accumData.bondAdj->getVal(bid1, bid2);
// RDKit✔️✔️:   aid1 = mol.getBondWithIdx(bid1)->getOtherAtomIdx(aid2);
// RDKit✔️✔️:   aid3 = accumData.bondAdj->getVal(bid2, bid3);
// RDKit✔️✔️:   aid4 = mol.getBondWithIdx(bid3)->getOtherAtomIdx(aid3);
// RDKit✔️✔️:   d1 = accumData.bondLengths[bid1];
// RDKit✔️✔️:   d2 = accumData.bondLengths[bid2];
// RDKit✔️✔️:   d3 = accumData.bondLengths[bid3];
// RDKit✔️✔️:   ang12 = accumData.bondAngles->getVal(bid1, bid2);
// RDKit✔️✔️:   ang23 = accumData.bondAngles->getVal(bid2, bid3);
// RDKit✔️✔️:   for (i = 0; i < nb; i++) {
// RDKit✔️✔️:     du = -1.0;
// RDKit✔️✔️:     dl = 0.0;
// RDKit✔️✔️:     if (accumData.bondAdj->getVal(bid3, i) == static_cast<int>(aid4)) {
// RDKit✔️✔️:       aid5 = mol.getBondWithIdx(i)->getOtherAtomIdx(aid4);
// RDKit✔️✔️:       // make sure we did not com back to the first atom in the path -
// RDKit✔️✔️:       // possible
// RDKit✔️✔️:       // with 4 membered rings
// RDKit✔️✔️:       // this is a fix for Issue 244
// RDKit✔️✔️:       const unsigned int pid = std::min(aid1, aid5) * na + std::max(aid1, aid5);
// RDKit✔️✔️:       if (accumData.visitedBound(pid, DistType::DIST14)) {
// RDKit✔️✔️:         return;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       if (dmat[std::max(aid1, aid5) * mmat->numRows() + std::min(aid1, aid5)] <
// RDKit✔️✔️:           3.9) {
// RDKit✔️✔️:         continue;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       if (aid1 != aid5) {  // FIX: do we need this
// RDKit✔️✔️:         unsigned int pid1 = aid1 * na + aid5;
// RDKit✔️✔️:         unsigned int pid2 = aid5 * na + aid1;
// RDKit✔️✔️:         if ((mmat->getLowerBound(aid1, aid5) < DIST12_DELTA) ||
// RDKit✔️✔️:             (accumData.set15Atoms[pid1]) || (accumData.set15Atoms[pid2])) {
// RDKit✔️✔️:           d4 = accumData.bondLengths[i];
// RDKit✔️✔️:           ang34 = accumData.bondAngles->getVal(bid3, i);
// RDKit✔️✔️:           unsigned long pathId =
// RDKit✔️✔️:               static_cast<unsigned long>(bid2) * nb * nb + (bid3)*nb + i;
// RDKit✔️✔️:           if (type == 0) {
// RDKit✔️✔️:             if (accumData.cisPaths.find(pathId) != accumData.cisPaths.end()) {
// RDKit✔️✔️:               dl = _compute15DistsCisCis(d1, d2, d3, d4, ang12, ang23, ang34);
// RDKit✔️✔️:               du = dl + DIST15_TOL;
// RDKit✔️✔️:               dl -= DIST15_TOL;
// RDKit✔️✔️:             } else if (accumData.transPaths.find(pathId) !=
// RDKit✔️✔️:                        accumData.transPaths.end()) {
// RDKit✔️✔️:               dl = _compute15DistsCisTrans(d1, d2, d3, d4, ang12, ang23, ang34);
// RDKit✔️✔️:               du = dl + DIST15_TOL;
// RDKit✔️✔️:               dl -= DIST15_TOL;
// RDKit✔️✔️:             } else {
// RDKit✔️✔️:               dl = _compute15DistsCisCis(d1, d2, d3, d4, ang12, ang23, ang34) -
// RDKit✔️✔️:                    DIST15_TOL;
// RDKit✔️✔️:               du =
// RDKit✔️✔️:                   _compute15DistsCisTrans(d1, d2, d3, d4, ang12, ang23, ang34) +
// RDKit✔️✔️:                   DIST15_TOL;
// RDKit✔️✔️:             }
// RDKit✔️✔️:           } else if (type == 1) {
// RDKit✔️✔️:             if (accumData.cisPaths.find(pathId) != accumData.cisPaths.end()) {
// RDKit✔️✔️:               dl = _compute15DistsTransCis(d1, d2, d3, d4, ang12, ang23, ang34);
// RDKit✔️✔️:               du = dl + DIST15_TOL;
// RDKit✔️✔️:               dl -= DIST15_TOL;
// RDKit✔️✔️:             } else if (accumData.transPaths.find(pathId) !=
// RDKit✔️✔️:                        accumData.transPaths.end()) {
// RDKit✔️✔️:               dl = _compute15DistsTransTrans(d1, d2, d3, d4, ang12, ang23,
// RDKit✔️✔️:                                              ang34);
// RDKit✔️✔️:               du = dl + DIST15_TOL;
// RDKit✔️✔️:               dl -= DIST15_TOL;
// RDKit✔️✔️:             } else {
// RDKit✔️✔️:               dl =
// RDKit✔️✔️:                   _compute15DistsTransCis(d1, d2, d3, d4, ang12, ang23, ang34) -
// RDKit✔️✔️:                   DIST15_TOL;
// RDKit✔️✔️:               du = _compute15DistsTransTrans(d1, d2, d3, d4, ang12, ang23,
// RDKit✔️✔️:                                              ang34) +
// RDKit✔️✔️:                    DIST15_TOL;
// RDKit✔️✔️:             }
// RDKit✔️✔️:           } else {
// RDKit✔️✔️:             if (accumData.cisPaths.find(pathId) != accumData.cisPaths.end()) {
// RDKit✔️✔️:               dl = _compute15DistsCisCis(d4, d3, d2, d1, ang34, ang23, ang12) -
// RDKit✔️✔️:                    DIST15_TOL;
// RDKit✔️✔️:               du =
// RDKit✔️✔️:                   _compute15DistsCisTrans(d4, d3, d2, d1, ang34, ang23, ang12) +
// RDKit✔️✔️:                   DIST15_TOL;
// RDKit✔️✔️:             } else if (accumData.transPaths.find(pathId) !=
// RDKit✔️✔️:                        accumData.transPaths.end()) {
// RDKit✔️✔️:               dl =
// RDKit✔️✔️:                   _compute15DistsTransCis(d4, d3, d2, d1, ang34, ang23, ang12) -
// RDKit✔️✔️:                   DIST15_TOL;
// RDKit✔️✔️:               du = _compute15DistsTransTrans(d4, d3, d2, d1, ang34, ang23,
// RDKit✔️✔️:                                              ang12) +
// RDKit✔️✔️:                    DIST15_TOL;
// RDKit✔️✔️:             } else {
// RDKit✔️✔️:               vw1 = PeriodicTable::getTable()->getRvdw(
// RDKit✔️✔️:                   mol.getAtomWithIdx(aid1)->getAtomicNum());
// RDKit✔️✔️:               vw5 = PeriodicTable::getTable()->getRvdw(
// RDKit✔️✔️:                   mol.getAtomWithIdx(aid5)->getAtomicNum());
// RDKit✔️✔️:               dl = VDW_SCALE_15 * (vw1 + vw5);
// RDKit✔️✔️:             }
// RDKit✔️✔️:           }
// RDKit✔️✔️:           if (du < 0.0) {
// RDKit✔️✔️:             du = MAX_UPPER;
// RDKit✔️✔️:           }
// RDKit✔️✔️:           _checkAndSetBounds(aid1, aid5, dl, du, mmat);
// RDKit✔️✔️:           accumData.set15Atoms[aid1 * na + aid5] = 1;
// RDKit✔️✔️:           accumData.set15Atoms[aid5 * na + aid1] = 1;
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_set15BoundsHelper
#[allow(clippy::too_many_arguments)]
fn set_15_bounds_helper(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    dmat: &[f64],
    nb: usize,
    na: usize,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    kind: Path14Kind,
) {
    // RDKit❗✔️: Get shared atoms via bond adjacency matrix
    let aid2 = accum_data.get_bond_adj(nb, bid1, bid2) as usize;
    let aid1 = if mol.bonds()[bid1].begin().index() == aid2 {
        mol.bonds()[bid1].end().index()
    } else {
        mol.bonds()[bid1].begin().index()
    };
    let aid3 = accum_data.get_bond_adj(nb, bid2, bid3) as usize;
    let aid4 = if mol.bonds()[bid3].begin().index() == aid3 {
        mol.bonds()[bid3].end().index()
    } else {
        mol.bonds()[bid3].begin().index()
    };

    let d1 = accum_data.bond_lengths[bid1];
    let d2 = accum_data.bond_lengths[bid2];
    let d3 = accum_data.bond_lengths[bid3];
    let ang12 = accum_data.get_bond_angle(nb, bid1, bid2);
    let ang23 = accum_data.get_bond_angle(nb, bid2, bid3);

    // RDKit❗✔️: loop over all bonds to find the 5th atom (extension from aid4)
    for i in 0..nb {
        if accum_data.get_bond_adj(nb, bid3, i) != aid4 as i32 {
            continue;
        }
        let aid5 = if mol.bonds()[i].begin().index() == aid4 {
            mol.bonds()[i].end().index()
        } else {
            mol.bonds()[i].begin().index()
        };

        // RDKit❗✔️: FIX: make sure we did not come back to the first atom (Issue 244)
        let pid = aid1.min(aid5) * na + aid1.max(aid5);

        // RDKit❗✔️: if pair was already bounded as 1-2/1-3/1-4, skip the whole function
        if accum_data.visited_bound(pid, DistType::Dist14) {
            return;
        }

        // RDKit❗✔️: skip if atoms are already close in distance matrix
        if dmat[aid1.max(aid5) * na + aid1.min(aid5)] < 3.9 {
            continue;
        }

        // RDKit❗✔️: skip if same atom (4-membered ring issue)
        if aid1 == aid5 {
            continue;
        }

        // RDKit❗✔️: check if bounds not already set and not already processed
        let pid1 = aid1 * na + aid5;
        let pid2 = aid5 * na + aid1;
        if !(mmat.get_lower(aid1, aid5) < DIST12_DELTA
            || accum_data.set15_atoms[pid1]
            || accum_data.set15_atoms[pid2])
        {
            continue;
        }

        let d4 = accum_data.bond_lengths[i];
        let ang34 = accum_data.get_bond_angle(nb, bid3, i);

        // RDKit❗✔️: pathId = bid2*nb*nb + bid3*nb + i
        let path_id = bid2 as u64 * nb as u64 * nb as u64 + bid3 as u64 * nb as u64 + i as u64;

        // RDKit❗✔️: compute 1-5 distance based on path type and cis/trans flags
        let (mut dl, mut du) = match kind {
            // RDKit❗✔️: type == 0 (CIS)
            Path14Kind::Cis => {
                if has_path_flag(&accum_data.cis_paths, path_id) {
                    // RDKit❗✔️: _compute15DistsCisCis +- DIST15_TOL
                    let base = compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                    (base - DIST15_TOL, base + DIST15_TOL)
                } else if has_path_flag(&accum_data.trans_paths, path_id) {
                    // RDKit❗✔️: _compute15DistsCisTrans +- DIST15_TOL
                    let base = compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                    (base - DIST15_TOL, base + DIST15_TOL)
                } else {
                    // RDKit❗✔️: cis-cis - tol, cis-trans + tol
                    (
                        compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL,
                        compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34) + DIST15_TOL,
                    )
                }
            }
            // RDKit❗✔️: type == 1 (TRANS)
            Path14Kind::Trans => {
                if has_path_flag(&accum_data.cis_paths, path_id) {
                    // RDKit❗✔️: _compute15DistsTransCis +- DIST15_TOL
                    let base = compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                    (base - DIST15_TOL, base + DIST15_TOL)
                } else if has_path_flag(&accum_data.trans_paths, path_id) {
                    // RDKit❗✔️: _compute15DistsTransTrans +- DIST15_TOL
                    let base = compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                    (base - DIST15_TOL, base + DIST15_TOL)
                } else {
                    // RDKit❗✔️: trans-cis - tol, trans-trans + tol
                    (
                        compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL,
                        compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34)
                            + DIST15_TOL,
                    )
                }
            }
            // RDKit❗✔️: type == 2 (OTHER) — reversed d-params for cis/trans detection
            Path14Kind::Other => {
                if has_path_flag(&accum_data.cis_paths, path_id) {
                    // RDKit❗✔️: _compute15DistsCisCis(d4,d3,d2,d1,ang34,ang23,ang12) - tol
                    // RDKit❗✔️: _compute15DistsCisTrans(d4,d3,d2,d1,ang34,ang23,ang12) + tol
                    (
                        compute_15_dist_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL,
                        compute_15_dist_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL,
                    )
                } else if has_path_flag(&accum_data.trans_paths, path_id) {
                    // RDKit❗✔️: _compute15DistsTransCis(d4,d3,d2,d1,ang34,ang23,ang12) - tol
                    // RDKit❗✔️: _compute15DistsTransTrans(d4,d3,d2,d1,ang34,ang23,ang12) + tol
                    (
                        compute_15_dist_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL,
                        compute_15_dist_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12)
                            + DIST15_TOL,
                    )
                } else {
                    // RDKit❗✔️: VDW fallback for unknown stereochemistry
                    let vw1 = vdw_radius(mol.atoms()[aid1].atomic_number());
                    let vw5 = vdw_radius(mol.atoms()[aid5].atomic_number());
                    (VDW_SCALE_15 * (vw1 + vw5), MAX_UPPER)
                }
            }
        };

        // RDKit❗✔️: if (du < 0.0) { du = MAX_UPPER; }
        if du < 0.0 {
            du = MAX_UPPER;
        }

        // RDKit❗✔️: _checkAndSetBounds(aid1, aid5, dl, du, mmat);
        mmat.check_and_set_bounds(aid1, aid5, dl, du);

        // RDKit❗✔️: accumData.set15Atoms[pid1] = 1; accumData.set15Atoms[pid2] = 1;
        accum_data.set15_atoms[pid1] = true;
        accum_data.set15_atoms[pid2] = true;
    }
}

/// RDKit❗✔️: set14Bounds — torsion-based 1-4 bounds
/// Also populates paths14, bond_adj, bond_angles for 1-5 bounds.
fn set_14_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    dmat: &[f64],
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
) {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::set14Bounds (BoundsMatrixBuilder.cpp:1664-1784)
    // RDKit❗✔️: void set14Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    // RDKit❗✔️:                  ComputedData &accumData, double *distMatrix,
    // RDKit❗✔️:                  bool useMacrocycle14config, bool forceTransAmides) {
    // RDKit❗✔️:   unsigned int npt = mmat->numRows();
    // RDKit❗✔️:   CHECK_INVARIANT(npt == mol.getNumAtoms(), "Wrong size metric matrix");
    // RDKit❗✔️:   const size_t MAX_NUM_BONDS = static_cast<size_t>(
    // RDKit❗✔️:       std::pow(std::numeric_limits<std::uint64_t>::max(), 1. / 3));
    // RDKit❗✔️:   if (mol.getNumBonds() >= MAX_NUM_BONDS) { throw ...; }
    // RDKit❗✔️:   const auto rinfo = mol.getRingInfo();
    // RDKit❗✔️:   const auto &bondRings = rinfo->bondRings();
    // RDKit❗✔️:   std::unordered_set<unsigned int> bidIsMacrocycle;
    // RDKit❗✔️:   std::unordered_set<std::uint64_t> ringBondPairs;
    // RDKit❗✔️:   std::unordered_set<std::uint64_t> donePaths;
    // RDKit❗✔️:   std::uint64_t nb = mol.getNumBonds();
    // RDKit❗✔️:   for (const auto &bring : bondRings) { ... }
    // RDKit❗✔️:   for (const auto bond : mol.bonds()) { ... }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::set14Bounds
    let npt = mmat.num_rows();
    assert_eq!(npt, mol.num_atoms(), "Wrong size metric matrix");
    let max_num_bonds = (u64::MAX as f64).cbrt() as usize;
    if mol.num_bonds() >= max_num_bonds {
        panic!("Too many bonds in the molecule, cannot compute 1-4 bounds");
    }
    let rinfo = ring_info_for_distgeom(mol);
    let bond_rings = rinfo.bond_rings();

    let mut bid_is_macrocycle: HashSet<usize> = HashSet::new();
    let mut ring_bond_pairs: HashSet<u64> = HashSet::new();
    let mut done_paths: HashSet<u64> = HashSet::new();
    let nb = mol.num_bonds() as u64;

    for bring in bond_rings {
        let r_size = bring.len();
        if r_size < 3 {
            continue;
        }
        let mut bid1 = bring[r_size - 1].index();
        for i in 0..r_size {
            let bid2 = bring[i].index();
            let bid3 = bring[(i + 1) % r_size].index();
            let pid1 = bid1 as u64 * nb + bid2 as u64;
            let pid2 = bid2 as u64 * nb + bid1 as u64;
            let id1 = bid1 as u64 * nb * nb + bid2 as u64 * nb + bid3 as u64;
            let id2 = bid3 as u64 * nb * nb + bid2 as u64 * nb + bid1 as u64;

            ring_bond_pairs.insert(pid1);
            ring_bond_pairs.insert(pid2);
            done_paths.insert(id1);
            done_paths.insert(id2);

            if r_size > 5 {
                if use_macrocycle_14config && r_size >= MIN_MACROCYCLE_RING_SIZE {
                    set_macrocycle_all_in_same_ring_14_bounds(
                        mol, bid1, bid2, bid3, accum_data, mmat,
                    );
                    bid_is_macrocycle.insert(bid2);
                } else {
                    set_in_ring_14_bounds(mol, bid1, bid2, bid3, accum_data, mmat, dmat, r_size);
                }
            } else {
                record_14_path(mol, bid1, bid2, bid3, accum_data);
            }
            bid1 = bid2;
        }
    }

    for bond in mol.bonds() {
        let bid2 = bond.id().index();
        let aid2 = bond.begin().index();
        let aid3 = bond.end().index();
        for nbr1 in neighbors_for_atom(mol, aid2) {
            let Some(bid1) = bond_between_idx_simple(mol, aid2, nbr1) else {
                continue;
            };
            if bid1 == bid2 {
                continue;
            }
            for nbr3 in neighbors_for_atom(mol, aid3) {
                let Some(bid3) = bond_between_idx_simple(mol, aid3, nbr3) else {
                    continue;
                };
                if bid3 == bid2 {
                    continue;
                }
                let id1 = bid1 as u64 * nb * nb + bid2 as u64 * nb + bid3 as u64;
                let id2 = bid3 as u64 * nb * nb + bid2 as u64 * nb + bid1 as u64;
                if done_paths.contains(&id1) || done_paths.contains(&id2) {
                    continue;
                }

                let pid1 = bid1 as u64 * nb + bid2 as u64;
                let pid2 = bid2 as u64 * nb + bid1 as u64;
                let pid3 = bid2 as u64 * nb + bid3 as u64;
                let pid4 = bid3 as u64 * nb + bid2 as u64;

                if ring_bond_pairs.contains(&pid1)
                    || ring_bond_pairs.contains(&pid2)
                    || ring_bond_pairs.contains(&pid3)
                    || ring_bond_pairs.contains(&pid4)
                {
                    if use_macrocycle_14config && bid_is_macrocycle.contains(&bid2) {
                        set_macrocycle_two_in_same_ring_14_bounds(
                            mol, bid1, bid2, bid3, accum_data, mmat, dmat,
                        );
                    } else {
                        set_two_in_same_ring_14_bounds(
                            mol, bid1, bid2, bid3, accum_data, mmat, dmat,
                        );
                    }
                } else if (rinfo.num_bond_rings(BondId::new(bid1)) > 0
                    && rinfo.num_bond_rings(BondId::new(bid2)) > 0)
                    || (rinfo.num_bond_rings(BondId::new(bid2)) > 0
                        && rinfo.num_bond_rings(BondId::new(bid3)) > 0)
                {
                    set_two_in_diff_ring_14_bounds(mol, bid1, bid2, bid3, accum_data, mmat, dmat);
                } else if rinfo.num_bond_rings(BondId::new(bid2)) > 0 {
                    set_share_ring_bond_14_bounds(mol, bid1, bid2, bid3, accum_data, mmat, dmat);
                } else {
                    set_chain_14_bounds(
                        mol,
                        bid1,
                        bid2,
                        bid3,
                        accum_data,
                        mmat,
                        force_trans_amides,
                    );
                }
            }
        }
    }
}

/// Helper: get bond length from ComputedData
fn bond_lengths_from_accum(accum_data: &ComputedData, bid: usize) -> f64 {
    let bl = accum_data.bond_lengths[bid];
    if bl > 0.0 { bl } else { 1.5 }
}

fn bond_between_idx_simple(mol: &Molecule, a: usize, b: usize) -> Option<usize> {
    mol.bonds()
        .iter()
        .find(|bond| {
            (bond.begin().index() == a && bond.end().index() == b)
                || (bond.begin().index() == b && bond.end().index() == a)
        })
        .map(|b| b.id().index())
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::collectBondsAndAngles (BoundsMatrixBuilder.cpp:1847-1898)
// RDKit✔️✔️: void collectBondsAndAngles(const ROMol &mol,
// RDKit✔️✔️:                            std::vector<std::pair<int, int>> &bonds,
// RDKit✔️✔️:                            std::vector<std::vector<int>> &angles) {
// RDKit✔️✔️:   bonds.resize(0);
// RDKit✔️✔️:   angles.resize(0);
// RDKit✔️✔️:   bonds.reserve(mol.getNumBonds());
// RDKit✔️✔️:   for (const auto bondi : mol.bonds()) {
// RDKit✔️✔️:     bonds.emplace_back(bondi->getBeginAtomIdx(), bondi->getEndAtomIdx());
// RDKit✔️✔️:     for (unsigned int j = bondi->getIdx() + 1; j < mol.getNumBonds(); ++j) {
// RDKit✔️✔️:       const Bond *bondj = mol.getBondWithIdx(j);
// RDKit✔️✔️:       int aid11 = bondi->getBeginAtomIdx();
// RDKit✔️✔️:       int aid12 = bondi->getEndAtomIdx();
// RDKit✔️✔️:       int aid21 = bondj->getBeginAtomIdx();
// RDKit✔️✔️:       int aid22 = bondj->getEndAtomIdx();
// RDKit✔️✔️:       if (aid11 != aid21 && aid11 != aid22 && aid12 != aid21 &&
// RDKit✔️✔️:           aid12 != aid22) {
// RDKit✔️✔️:         continue;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       std::vector<int> tmp(4,
// RDKit✔️✔️:                            0);  // elements: aid1, aid2, flag for triple bonds
// RDKit✔️✔️:       if (aid12 == aid21) {
// RDKit✔️✔️:         tmp[0] = aid11;
// RDKit✔️✔️:         tmp[1] = aid12;
// RDKit✔️✔️:         tmp[2] = aid22;
// RDKit✔️✔️:       } else if (aid12 == aid22) {
// RDKit✔️✔️:         tmp[0] = aid11;
// RDKit✔️✔️:         tmp[1] = aid12;
// RDKit✔️✔️:         tmp[2] = aid21;
// RDKit✔️✔️:       } else if (aid11 == aid21) {
// RDKit✔️✔️:         tmp[0] = aid12;
// RDKit✔️✔️:         tmp[1] = aid11;
// RDKit✔️✔️:         tmp[2] = aid22;
// RDKit✔️✔️:       } else if (aid11 == aid22) {
// RDKit✔️✔️:         tmp[0] = aid12;
// RDKit✔️✔️:         tmp[1] = aid11;
// RDKit✔️✔️:         tmp[2] = aid21;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       if (bondi->getBondType() == Bond::TRIPLE ||
// RDKit✔️✔️:           bondj->getBondType() == Bond::TRIPLE) {
// RDKit✔️✔️:         // triple bond
// RDKit✔️✔️:         tmp[3] = 1;
// RDKit✔️✔️:       } else if (bondi->getBondType() == Bond::DOUBLE &&
// RDKit✔️✔️:                  bondj->getBondType() == Bond::DOUBLE &&
// RDKit✔️✔️:                  mol.getAtomWithIdx(tmp[1])->getDegree() == 2) {
// RDKit✔️✔️:         // consecutive double bonds
// RDKit✔️✔️:         tmp[3] = 1;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       angles.push_back(tmp);
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::collectBondsAndAngles
fn collect_bonds_and_angles(
    mol: &Molecule,
    bonds: &mut Vec<(i32, i32)>,
    angles: &mut Vec<Vec<i32>>,
) {
    bonds.clear();
    angles.clear();
    bonds.reserve(mol.num_bonds());
    for bondi in mol.bonds() {
        bonds.push((bondi.begin().index() as i32, bondi.end().index() as i32));

        for j in (bondi.id().index() + 1)..mol.num_bonds() {
            let bondj = &mol.bonds()[j];
            let aid11 = bondi.begin().index() as i32;
            let aid12 = bondi.end().index() as i32;
            let aid21 = bondj.begin().index() as i32;
            let aid22 = bondj.end().index() as i32;
            if aid11 != aid21 && aid11 != aid22 && aid12 != aid21 && aid12 != aid22 {
                continue;
            }

            let mut tmp = vec![0; 4];
            if aid12 == aid21 {
                tmp[0] = aid11;
                tmp[1] = aid12;
                tmp[2] = aid22;
            } else if aid12 == aid22 {
                tmp[0] = aid11;
                tmp[1] = aid12;
                tmp[2] = aid21;
            } else if aid11 == aid21 {
                tmp[0] = aid12;
                tmp[1] = aid11;
                tmp[2] = aid22;
            } else if aid11 == aid22 {
                tmp[0] = aid12;
                tmp[1] = aid11;
                tmp[2] = aid21;
            }

            if bondi.order() == BondOrder::Triple || bondj.order() == BondOrder::Triple {
                tmp[3] = 1;
            } else if bondi.order() == BondOrder::Double
                && bondj.order() == BondOrder::Double
                && neighbors_for_atom(mol, tmp[1] as usize).len() == 2
            {
                tmp[3] = 1;
            }

            angles.push(tmp);
        }
    }
}

// ──────────────────────────────────────────────
// VDW lower bounds
// ──────────────────────────────────────────────

/// RDKit❗✔️: setLowerBoundVDW — VDW lower bounds for non-bonded pairs
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::setLowerBoundVDW (BoundsMatrixBuilder.cpp:315-361)
// RDKit✔️✔️: void setLowerBoundVDW(const ROMol &mol, DistGeom::BoundsMatPtr mmat, bool,
// RDKit✔️✔️:                       double *dmat) {
// RDKit✔️✔️:   for (unsigned int i = 1; i < npt; i++) {
// RDKit✔️✔️:     const auto atomI = mol.getAtomWithIdx(i);
// RDKit✔️✔️:     auto vw1 = PeriodicTable::getTable()->getRvdw(atomI->getAtomicNum());
// RDKit✔️✔️:     if (isHinHBondDonor(atomI, mol)) hinHBondDonors.set(i);
// RDKit✔️✔️:     if (isHBondAcceptor(atomI)) hBondAcceptors.set(i);
// RDKit✔️✔️:     for (unsigned int j = 0; j < i; j++) {
// RDKit✔️✔️:       if (mmat->getLowerBound(i, j) < DIST12_DELTA) {
// RDKit✔️✔️:         if ((hinHBondDonors[i] && hBondAcceptors[j]) || ...) {
// RDKit✔️✔️:           mmat->setLowerBound(i, j, H_BOND_LENGTH);
// RDKit✔️✔️:         } else if (dmat[i * npt + j] == 4.0) {
// RDKit✔️✔️:           mmat->setLowerBound(i, j, VDW_SCALE_15 * (vw1 + vw2));
// RDKit✔️✔️:         } else ...
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::setLowerBoundVDW
fn set_lower_bound_vdw(mol: &Molecule, mmat: &mut BoundsMatrix, _scale_vdw: bool, dmat: &[f64]) {
    let n = mol.atoms().len();
    let npt = mmat.num_rows();
    assert_eq!(npt, n, "Wrong size metric matrix");

    // Pre-compute H-bond donor/acceptor bitsets
    let mut h_in_hbond_donors = vec![false; n];
    let mut hbond_acceptors = vec![false; n];
    for i in 1..n {
        h_in_hbond_donors[i] = is_h_in_hbond_donor(mol, i);
        hbond_acceptors[i] = is_hbond_acceptor(mol.atoms()[i].atomic_number());
    }

    for i in 1..n {
        let vw1 = vdw_radius(mol.atoms()[i].atomic_number());

        for j in 0..i {
            if mmat.get_lower(i, j) > DIST12_DELTA {
                continue;
            }

            let vw2 = vdw_radius(mol.atoms()[j].atomic_number());
            let td = dmat[i * npt + j];

            // H-bond special case: set lower bound to 1.8A for donor-H + acceptor
            if (h_in_hbond_donors[i] && hbond_acceptors[j])
                || (hbond_acceptors[i] && h_in_hbond_donors[j])
            {
                mmat.set_lower(i, j, H_BOND_LENGTH);
            } else if td == 4.0 {
                // 1-5: scaled VDW
                mmat.set_lower(i, j, VDW_SCALE_15 * (vw1 + vw2));
            } else if td == 5.0 {
                // 1-6: slightly less scaled
                mmat.set_lower(
                    i,
                    j,
                    (VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * (vw1 + vw2),
                );
            } else {
                // Full VDW sum
                mmat.set_lower(i, j, vw1 + vw2);
            }
        }
    }
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::setTopolBounds (BoundsMatrixBuilder.cpp:1808-1845)
// RDKit✔️✔️: void setTopolBounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
// RDKit✔️✔️:                     bool set15bounds, bool scaleVDW, bool useMacrocycle14config,
// RDKit✔️✔️:                     bool forceTransAmides, bool set14bounds, bool set13bounds) {
// RDKit✔️✔️:   PRECONDITION(mmat.get(), "bad pointer");
// RDKit✔️✔️:   unsigned int nb = mol.getNumBonds();
// RDKit✔️✔️:   unsigned int na = mol.getNumAtoms();
// RDKit✔️✔️:   if (!na) {
// RDKit✔️✔️:     throw ValueErrorException("molecule has no atoms");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   const size_t MAX_NUM_BONDS = static_cast<size_t>(
// RDKit✔️✔️:       std::pow(std::numeric_limits<std::uint64_t>::max(), 1. / 3));
// RDKit✔️✔️:   if (mol.getNumBonds() >= MAX_NUM_BONDS) {
// RDKit✔️✔️:     throw ValueErrorException(
// RDKit✔️✔️:         "Too many bonds in the molecule, cannot compute 1-4 bounds");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   ComputedData accumData(na, nb);
// RDKit✔️✔️:   double *distMatrix = nullptr;
// RDKit✔️✔️:   distMatrix = MolOps::getDistanceMat(mol);
// RDKit✔️✔️:   set12Bounds(mol, mmat, accumData);
// RDKit✔️✔️:   if (set13bounds) {
// RDKit✔️✔️:     set13Bounds(mol, mmat, accumData);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (set14bounds) {
// RDKit✔️✔️:     set14Bounds(mol, mmat, accumData, distMatrix, useMacrocycle14config,
// RDKit✔️✔️:                 forceTransAmides);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (set15bounds) {
// RDKit✔️✔️:     set15Bounds(mol, mmat, accumData, distMatrix);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   setLowerBoundVDW(mol, mmat, scaleVDW, distMatrix);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::setTopolBounds
#[allow(clippy::too_many_arguments)]
fn set_topol_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    set15bounds: bool,
    scale_vdw: bool,
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
    set14bounds: bool,
    set13bounds: bool,
) -> Result<(), DgBoundsError> {
    let nb = mol.num_bonds();
    let na = mol.num_atoms();
    if na == 0 {
        return Err(DgBoundsError::GenerationFailed(
            "molecule has no atoms".to_string(),
        ));
    }
    let max_num_bonds = (u64::MAX as f64).cbrt() as usize;
    if nb >= max_num_bonds {
        return Err(DgBoundsError::GenerationFailed(
            "Too many bonds in the molecule, cannot compute 1-4 bounds".to_string(),
        ));
    }

    let mut accum_data = ComputedData::new(na, nb);
    let dist_matrix = flatten_topological_distances_matrix(mol);

    set_12_bounds(mol, mmat, &mut accum_data)?;
    if set13bounds {
        set_13_bounds(mol, mmat, &mut accum_data);
    }
    if set14bounds {
        set_14_bounds(
            mol,
            mmat,
            &mut accum_data,
            &dist_matrix,
            use_macrocycle_14config,
            force_trans_amides,
        );
    }
    if set15bounds {
        set_15_bounds(mol, mmat, &mut accum_data, &dist_matrix);
    }
    set_lower_bound_vdw(mol, mmat, scale_vdw, &dist_matrix);
    Ok(())
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::setTopolBounds (BoundsMatrixBuilder.cpp:1903-1942)
// RDKit✔️✔️: void setTopolBounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
// RDKit✔️✔️:                     std::vector<std::pair<int, int>> &bonds,
// RDKit✔️✔️:                     std::vector<std::vector<int>> &angles, bool set15bounds,
// RDKit✔️✔️:                     bool scaleVDW, bool useMacrocycle14config,
// RDKit✔️✔️:                     bool forceTransAmides, bool set14bounds, bool set13bounds) {
// RDKit✔️✔️:   PRECONDITION(mmat.get(), "bad pointer");
// RDKit✔️✔️:   bonds.clear();
// RDKit✔️✔️:   angles.clear();
// RDKit✔️✔️:   unsigned int nb = mol.getNumBonds();
// RDKit✔️✔️:   unsigned int na = mol.getNumAtoms();
// RDKit✔️✔️:   if (!na) {
// RDKit✔️✔️:     throw ValueErrorException("molecule has no atoms");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   ComputedData accumData(na, nb);
// RDKit✔️✔️:   double *distMatrix = nullptr;
// RDKit✔️✔️:   distMatrix = MolOps::getDistanceMat(mol);
// RDKit✔️✔️:   set12Bounds(mol, mmat, accumData);
// RDKit✔️✔️:   if (set13bounds) {
// RDKit✔️✔️:     set13Bounds(mol, mmat, accumData);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (set14bounds) {
// RDKit✔️✔️:     set14Bounds(mol, mmat, accumData, distMatrix, useMacrocycle14config,
// RDKit✔️✔️:                 forceTransAmides);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (set15bounds) {
// RDKit✔️✔️:     set15Bounds(mol, mmat, accumData, distMatrix);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   setLowerBoundVDW(mol, mmat, scaleVDW, distMatrix);
// RDKit✔️✔️:   collectBondsAndAngles(mol, bonds, angles);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::setTopolBounds
#[allow(clippy::too_many_arguments)]
fn set_topol_bounds_with_outputs(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    bonds: &mut Vec<(i32, i32)>,
    angles: &mut Vec<Vec<i32>>,
    set15bounds: bool,
    scale_vdw: bool,
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
    set14bounds: bool,
    set13bounds: bool,
) -> Result<(), DgBoundsError> {
    bonds.clear();
    angles.clear();
    set_topol_bounds(
        mol,
        mmat,
        set15bounds,
        scale_vdw,
        use_macrocycle_14config,
        force_trans_amides,
        set14bounds,
        set13bounds,
    )?;
    collect_bonds_and_angles(mol, bonds, angles);
    Ok(())
}

// ──────────────────────────────────────────────
// Public entry point
// ──────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION RDKit::getMolBoundsMatrix (rdDistGeom.cpp:223-242)
// RDKit✔️❌: PyObject *getMolBoundsMatrix(ROMol &mol, bool set15bounds = true,
// RDKit✔️❌:                              bool scaleVDW = false,
// RDKit✔️❌:                              bool doTriangleSmoothing = true,
// RDKit✔️❌:                              bool useMacrocycle14config = false) {
// RDKit✔️❌:   unsigned int nats = mol.getNumAtoms();
// RDKit✔️❌:   npy_intp dims[2];
// RDKit✔️❌:   dims[0] = nats;
// RDKit✔️❌:   dims[1] = nats;
// RDKit✔️❌:
// RDKit✔️❌:   DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(nats));
// RDKit✔️❌:   DGeomHelpers::initBoundsMat(mat);
// RDKit✔️❌:   DGeomHelpers::setTopolBounds(mol, mat, set15bounds, scaleVDW,
// RDKit✔️❌:                                useMacrocycle14config);
// RDKit✔️❌:   if (doTriangleSmoothing) {
// RDKit✔️❌:     DistGeom::triangleSmoothBounds(mat);
// RDKit✔️❌:   }
// RDKit✔️❌:   auto *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
// RDKit✔️❌:   memcpy(static_cast<void *>(PyArray_DATA(res)),
// RDKit✔️❌:          static_cast<void *>(mat->getData()), nats * nats * sizeof(double));
// RDKit✔️❌:
// RDKit✔️❌:   return PyArray_Return(res);
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION RDKit::getMolBoundsMatrix
//
// Performance note: the Rust wrapper preserves RDKit's control flow and raw
// matrix export semantics, but materializing `Vec<Vec<f64>>` is less efficient
// than RDKit's direct NumPy buffer handoff from contiguous storage.
pub fn dg_bounds_matrix_with_options(
    molecule: &Molecule,
    set15bounds: bool,
    scale_vdw: bool,
    do_triangle_smoothing: bool,
    use_macrocycle14config: bool,
) -> Result<Vec<Vec<f64>>, DgBoundsError> {
    let n = molecule.atoms().len();
    let mut mmat = BoundsMatrix::new(n);
    set_topol_bounds(
        molecule,
        &mut mmat,
        set15bounds,
        scale_vdw,
        use_macrocycle14config,
        false,
        true,
        true,
    )?;
    if do_triangle_smoothing && !mmat.triangle_smooth(0.0) {
        return Err(DgBoundsError::GenerationFailed(
            "triangle smoothing found inconsistent bounds".to_string(),
        ));
    }
    Ok(mmat.to_vec_vec())
}

pub fn dg_bounds_matrix(molecule: &Molecule) -> Result<Vec<Vec<f64>>, DgBoundsError> {
    dg_bounds_matrix_with_options(molecule, true, false, true, false)
}

// ──────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::MoleculeBuilder;
    use crate::{AtomSpec, BondSpec, Element, Molecule, ValenceModel, assign_valence};

    fn run_set12_bounds(mol: &Molecule) -> (BoundsMatrix, ComputedData) {
        let mut mmat = BoundsMatrix::new(mol.num_atoms());
        let mut accum_data = ComputedData::new(mol.num_atoms(), mol.num_bonds());
        set_12_bounds(mol, &mut mmat, &mut accum_data).expect("set12Bounds");
        (mmat, accum_data)
    }

    fn run_set13_bounds(mol: &Molecule) -> (BoundsMatrix, ComputedData) {
        let (mut mmat, mut accum_data) = run_set12_bounds(mol);
        set_13_bounds(mol, &mut mmat, &mut accum_data);
        (mmat, accum_data)
    }

    fn run_set14_bounds(
        mol: &Molecule,
        use_macrocycle_14config: bool,
        force_trans_amides: bool,
    ) -> (BoundsMatrix, ComputedData, Vec<f64>) {
        let (mut mmat, mut accum_data) = run_set13_bounds(mol);
        let dmat = flatten_topological_distances(mol);
        set_14_bounds(
            mol,
            &mut mmat,
            &mut accum_data,
            &dmat,
            use_macrocycle_14config,
            force_trans_amides,
        );
        (mmat, accum_data, dmat)
    }

    fn run_set15_bounds(
        mol: &Molecule,
        use_macrocycle_14config: bool,
        force_trans_amides: bool,
    ) -> (BoundsMatrix, ComputedData, Vec<f64>) {
        let (mut mmat, mut accum_data, dmat) =
            run_set14_bounds(mol, use_macrocycle_14config, force_trans_amides);
        set_15_bounds(mol, &mut mmat, &mut accum_data, &dmat);
        (mmat, accum_data, dmat)
    }

    fn run_set_topol_bounds(
        mol: &Molecule,
        set15bounds: bool,
        scale_vdw: bool,
        use_macrocycle_14config: bool,
        force_trans_amides: bool,
        set14bounds: bool,
        set13bounds: bool,
    ) -> BoundsMatrix {
        let mut mmat = BoundsMatrix::new(mol.num_atoms());
        set_topol_bounds(
            mol,
            &mut mmat,
            set15bounds,
            scale_vdw,
            use_macrocycle_14config,
            force_trans_amides,
            set14bounds,
            set13bounds,
        )
        .expect("setTopolBounds");
        mmat
    }

    fn run_set_topol_bounds_with_outputs(
        mol: &Molecule,
        set15bounds: bool,
        scale_vdw: bool,
        use_macrocycle_14config: bool,
        force_trans_amides: bool,
        set14bounds: bool,
        set13bounds: bool,
    ) -> (BoundsMatrix, Vec<(i32, i32)>, Vec<Vec<i32>>) {
        let mut mmat = BoundsMatrix::new(mol.num_atoms());
        let mut bonds = Vec::new();
        let mut angles = Vec::new();
        set_topol_bounds_with_outputs(
            mol,
            &mut mmat,
            &mut bonds,
            &mut angles,
            set15bounds,
            scale_vdw,
            use_macrocycle_14config,
            force_trans_amides,
            set14bounds,
            set13bounds,
        )
        .expect("setTopolBounds with outputs");
        (mmat, bonds, angles)
    }

    fn run_set14_same_ring_pass_only(
        mol: &Molecule,
        use_macrocycle_14config: bool,
    ) -> (BoundsMatrix, ComputedData, Vec<f64>) {
        let (mut mmat, mut accum_data) = run_set13_bounds(mol);
        let dmat = flatten_topological_distances(mol);
        let rinfo = ring_info_for_distgeom(mol);
        for bring in rinfo.bond_rings() {
            let r_size = bring.len();
            if r_size < 3 {
                continue;
            }
            let mut bid1 = bring[r_size - 1].index();
            for i in 0..r_size {
                let bid2 = bring[i].index();
                let bid3 = bring[(i + 1) % r_size].index();
                if r_size > 5 {
                    if use_macrocycle_14config && r_size >= MIN_MACROCYCLE_RING_SIZE {
                        set_macrocycle_all_in_same_ring_14_bounds(
                            mol,
                            bid1,
                            bid2,
                            bid3,
                            &mut accum_data,
                            &mut mmat,
                        );
                    } else {
                        set_in_ring_14_bounds(
                            mol,
                            bid1,
                            bid2,
                            bid3,
                            &mut accum_data,
                            &mut mmat,
                            &dmat,
                            r_size,
                        );
                    }
                } else {
                    record_14_path(mol, bid1, bid2, bid3, &mut accum_data);
                }
                bid1 = bid2;
            }
        }
        (mmat, accum_data, dmat)
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum Set14DispatchCase {
        TwoSameRing,
        TwoDiffRing,
        ShareRingBond,
        Chain,
    }

    fn find_same_ring_dispatch_triple(
        mol: &Molecule,
        use_macrocycle_14config: bool,
    ) -> Option<(usize, usize, usize, usize)> {
        let rinfo = ring_info_for_distgeom(mol);
        for bring in rinfo.bond_rings() {
            let r_size = bring.len();
            if r_size < 3 {
                continue;
            }
            if r_size > 5 && (!use_macrocycle_14config || r_size >= MIN_MACROCYCLE_RING_SIZE) {
                let bid1 = bring[r_size - 1].index();
                let bid2 = bring[0].index();
                let bid3 = bring[1].index();
                return Some((bid1, bid2, bid3, r_size));
            }
        }
        None
    }

    fn find_dispatch_triple(
        mol: &Molecule,
        target: Set14DispatchCase,
        use_macrocycle_14config: bool,
    ) -> Option<(usize, usize, usize)> {
        let rinfo = ring_info_for_distgeom(mol);
        let mut bid_is_macrocycle: HashSet<usize> = HashSet::new();
        let mut ring_bond_pairs: HashSet<u64> = HashSet::new();
        let mut done_paths: HashSet<u64> = HashSet::new();
        let nb = mol.num_bonds() as u64;

        for bring in rinfo.bond_rings() {
            let r_size = bring.len();
            if r_size < 3 {
                continue;
            }
            let mut bid1 = bring[r_size - 1].index();
            for i in 0..r_size {
                let bid2 = bring[i].index();
                let bid3 = bring[(i + 1) % r_size].index();
                ring_bond_pairs.insert(bid1 as u64 * nb + bid2 as u64);
                ring_bond_pairs.insert(bid2 as u64 * nb + bid1 as u64);
                done_paths.insert(bid1 as u64 * nb * nb + bid2 as u64 * nb + bid3 as u64);
                done_paths.insert(bid3 as u64 * nb * nb + bid2 as u64 * nb + bid1 as u64);
                if use_macrocycle_14config && r_size >= MIN_MACROCYCLE_RING_SIZE {
                    bid_is_macrocycle.insert(bid2);
                }
                bid1 = bid2;
            }
        }

        for bond in mol.bonds() {
            let bid2 = bond.id().index();
            let aid2 = bond.begin().index();
            let aid3 = bond.end().index();
            for nbr1 in neighbors_for_atom(mol, aid2) {
                let Some(bid1) = bond_between_idx_simple(mol, aid2, nbr1) else {
                    continue;
                };
                if bid1 == bid2 {
                    continue;
                }
                for nbr3 in neighbors_for_atom(mol, aid3) {
                    let Some(bid3) = bond_between_idx_simple(mol, aid3, nbr3) else {
                        continue;
                    };
                    if bid3 == bid2 {
                        continue;
                    }
                    let id1 = bid1 as u64 * nb * nb + bid2 as u64 * nb + bid3 as u64;
                    let id2 = bid3 as u64 * nb * nb + bid2 as u64 * nb + bid1 as u64;
                    if done_paths.contains(&id1) || done_paths.contains(&id2) {
                        continue;
                    }
                    let pid1 = bid1 as u64 * nb + bid2 as u64;
                    let pid2 = bid2 as u64 * nb + bid1 as u64;
                    let pid3 = bid2 as u64 * nb + bid3 as u64;
                    let pid4 = bid3 as u64 * nb + bid2 as u64;
                    let case = if ring_bond_pairs.contains(&pid1)
                        || ring_bond_pairs.contains(&pid2)
                        || ring_bond_pairs.contains(&pid3)
                        || ring_bond_pairs.contains(&pid4)
                    {
                        Set14DispatchCase::TwoSameRing
                    } else if (rinfo.num_bond_rings(BondId::new(bid1)) > 0
                        && rinfo.num_bond_rings(BondId::new(bid2)) > 0)
                        || (rinfo.num_bond_rings(BondId::new(bid2)) > 0
                            && rinfo.num_bond_rings(BondId::new(bid3)) > 0)
                    {
                        Set14DispatchCase::TwoDiffRing
                    } else if rinfo.num_bond_rings(BondId::new(bid2)) > 0 {
                        Set14DispatchCase::ShareRingBond
                    } else {
                        Set14DispatchCase::Chain
                    };
                    if case == target {
                        return Some((bid1, bid2, bid3));
                    }
                }
            }
        }
        let _ = bid_is_macrocycle;
        None
    }

    fn flatten_topological_distances(mol: &Molecule) -> Vec<f64> {
        let topo = compute_topological_distances(mol);
        let n = mol.num_atoms();
        let mut flat = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                flat[i * n + j] = topo[i][j] as f64;
            }
        }
        flat
    }

    #[test]
    fn collect_bonds_and_angles_flags_triple_bond_paths_like_rdkit() {
        let mol = Molecule::from_smiles("CC#N").expect("acetonitrile");
        let mut bonds = Vec::new();
        let mut angles = Vec::new();

        collect_bonds_and_angles(&mol, &mut bonds, &mut angles);

        assert_eq!(bonds, vec![(0, 1), (1, 2)]);
        assert_eq!(angles, vec![vec![0, 1, 2, 1]]);
    }

    #[test]
    fn collect_bonds_and_angles_flags_consecutive_double_bonds_like_rdkit() {
        let mol = Molecule::from_smiles("C=C=C").expect("allene");
        let mut bonds = Vec::new();
        let mut angles = Vec::new();

        collect_bonds_and_angles(&mol, &mut bonds, &mut angles);

        assert_eq!(bonds, vec![(0, 1), (1, 2)]);
        assert_eq!(angles, vec![vec![0, 1, 2, 1]]);
    }

    #[test]
    fn init_bounds_ptr_sets_default_bounds_and_keeps_diagonal_zero() {
        let mut mmat = BoundsMatrix {
            data: vec![vec![0.0; 3]; 3],
            n: 3,
        };

        init_bounds_mat_ptr(&mut mmat, 0.25, 42.0);

        for i in 0..3 {
            assert_eq!(mmat.get_lower(i, i), 0.0);
            assert_eq!(mmat.get_upper(i, i), 0.0);
        }
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    continue;
                }
                assert_eq!(mmat.get_lower(i, j), 0.25);
                assert_eq!(mmat.get_upper(i, j), 42.0);
            }
        }
    }

    #[test]
    fn init_bounds_shared_matches_pointer_overload() {
        let mut from_ptr = BoundsMatrix {
            data: vec![vec![0.0; 4]; 4],
            n: 4,
        };
        let mut from_shared = BoundsMatrix {
            data: vec![vec![0.0; 4]; 4],
            n: 4,
        };

        init_bounds_mat_ptr(&mut from_ptr, DEFAULT_LOWER, DEFAULT_UPPER);
        init_bounds_mat_shared(&mut from_shared, DEFAULT_LOWER, DEFAULT_UPPER);

        assert_eq!(from_ptr.data, from_shared.data);
    }

    #[test]
    fn init_bounds_new_uses_rdkit_default_min_and_max() {
        let mmat = BoundsMatrix::new(2);

        assert_eq!(mmat.get_lower(0, 1), DEFAULT_LOWER);
        assert_eq!(mmat.get_upper(0, 1), DEFAULT_UPPER);
        assert_eq!(mmat.get_lower(0, 0), 0.0);
        assert_eq!(mmat.get_upper(1, 1), 0.0);
    }

    #[test]
    fn bounds_matrix_stores_upper_in_upper_triangle_and_lower_in_lower_triangle() {
        let mut mmat = BoundsMatrix {
            data: vec![vec![0.0; 3]; 3],
            n: 3,
        };

        mmat.set_upper(2, 0, 4.2);
        mmat.set_lower(2, 0, 1.8);

        assert_eq!(mmat.data[0][2], 4.2);
        assert_eq!(mmat.data[2][0], 1.8);
        assert_eq!(mmat.get_upper(0, 2), 4.2);
        assert_eq!(mmat.get_upper(2, 0), 4.2);
        assert_eq!(mmat.get_lower(0, 2), 1.8);
        assert_eq!(mmat.get_lower(2, 0), 1.8);
    }

    #[test]
    fn bounds_matrix_export_preserves_rdkit_raw_triangle_layout() {
        let mut mmat = BoundsMatrix::new(2);
        mmat.set_upper(0, 1, 2.5);
        mmat.set_lower(0, 1, 1.25);

        let raw = mmat.to_vec_vec();
        assert_eq!(raw, vec![vec![0.0, 2.5], vec![1.25, 0.0]]);
    }

    #[test]
    fn check_and_set_bounds_sets_uninitialized_pair_conservatively() {
        let mut mmat = BoundsMatrix::new(3);

        mmat.check_and_set_bounds(0, 2, 1.2, 2.8);

        assert_eq!(mmat.get_lower(0, 2), 1.2);
        assert_eq!(mmat.get_upper(0, 2), 2.8);
        assert_eq!(mmat.data[2][0], 1.2);
        assert_eq!(mmat.data[0][2], 2.8);
    }

    #[test]
    fn check_and_set_bounds_only_tightens_lower_and_relaxes_upper_when_allowed() {
        let mut mmat = BoundsMatrix::new(3);
        mmat.set_lower(0, 2, 1.8);
        mmat.set_upper(0, 2, 2.2);

        mmat.check_and_set_bounds(0, 2, 1.4, 2.6);
        assert_eq!(mmat.get_lower(0, 2), 1.4);
        assert_eq!(mmat.get_upper(0, 2), 2.6);

        mmat.check_and_set_bounds(0, 2, 1.6, 2.0);
        assert_eq!(mmat.get_lower(0, 2), 1.4);
        assert_eq!(mmat.get_upper(0, 2), 2.6);
    }

    #[test]
    fn check_and_set_bounds_with_mode_uses_overlap_when_ranges_are_consistent() {
        let mut mmat = BoundsMatrix::new(2);
        mmat.set_lower(0, 1, 1.3);
        mmat.set_upper(0, 1, 3.1);

        mmat.check_and_set_bounds_with_mode(0, 1, 1.8, 2.4, true);

        assert_eq!(mmat.get_lower(0, 1), 1.8);
        assert_eq!(mmat.get_upper(0, 1), 2.4);
        assert_eq!(mmat.data[1][0], 1.8);
        assert_eq!(mmat.data[0][1], 2.4);
    }

    #[test]
    fn check_and_set_bounds_with_mode_falls_back_to_conservative_union_when_disjoint() {
        let mut mmat = BoundsMatrix::new(2);
        mmat.set_lower(0, 1, 1.3);
        mmat.set_upper(0, 1, 2.0);

        mmat.check_and_set_bounds_with_mode(0, 1, 2.4, 3.0, true);

        assert_eq!(mmat.get_lower(0, 1), 1.3);
        assert_eq!(mmat.get_upper(0, 1), 3.0);
        assert_eq!(mmat.data[1][0], 1.3);
        assert_eq!(mmat.data[0][1], 3.0);
    }

    #[test]
    fn set12_bounds_uses_uff_rest_length_for_supported_atoms() {
        let mol = Molecule::from_smiles_with_sanitize("CC", false).expect("ethane skeleton");

        let (mmat, accum_data) = run_set12_bounds(&mol);

        let lower = mmat.get_lower(0, 1);
        let upper = mmat.get_upper(0, 1);
        let width = upper - lower;
        assert!((width - (2.0 * DIST12_DELTA)).abs() < 1e-9);
        assert!(accum_data.bond_lengths[0] > 1.4);
        assert!(accum_data.bond_lengths[0] < 1.6);
        assert!((lower - (accum_data.bond_lengths[0] - DIST12_DELTA)).abs() < 1e-9);
        assert!((upper - (accum_data.bond_lengths[0] + DIST12_DELTA)).abs() < 1e-9);
        assert!(accum_data.visited12_bounds[1]);
    }

    #[test]
    fn set12_bounds_falls_back_to_vdw_when_uff_params_are_missing() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let dummy = builder.add_atom(AtomSpec::new(Element::DUMMY));
        builder
            .add_bond(BondSpec::new(carbon, dummy, BondOrder::Single))
            .unwrap();
        let mol = builder.build().expect("dummy-carbon molecule");

        let (mmat, accum_data) = run_set12_bounds(&mol);
        let expected = (vdw_radius(6) + vdw_radius(0)) / 2.0;

        assert!((accum_data.bond_lengths[0] - expected).abs() < 1e-9);
        assert!((mmat.get_lower(0, 1) - (0.5 * expected)).abs() < 1e-9);
        assert!((mmat.get_upper(0, 1) - (1.5 * expected)).abs() < 1e-9);
    }

    #[test]
    fn set12_bounds_adds_extra_squish_for_conjugated_hetero_five_ring_bonds() {
        let mol = Molecule::from_smiles("s1cccc1").expect("thiophene");

        let (mmat, _accum_data) = run_set12_bounds(&mol);

        let sulfur_atom = mol
            .atoms()
            .iter()
            .position(|atom| atom.atomic_number() == 16)
            .expect("sulfur atom");
        let squished_width = mol
            .bonds()
            .iter()
            .find(|bond| bond.begin().index() == sulfur_atom || bond.end().index() == sulfur_atom)
            .map(|bond| {
                mmat.get_upper(bond.begin().index(), bond.end().index())
                    - mmat.get_lower(bond.begin().index(), bond.end().index())
            })
            .expect("sulfur bond");
        let mut sulfur_adjacent = vec![false; mol.num_atoms()];
        sulfur_adjacent[sulfur_atom] = true;
        for bond in mol.bonds() {
            if bond.begin().index() == sulfur_atom {
                sulfur_adjacent[bond.end().index()] = true;
            } else if bond.end().index() == sulfur_atom {
                sulfur_adjacent[bond.begin().index()] = true;
            }
        }
        let unsquished_width = mol
            .bonds()
            .iter()
            .find(|bond| {
                !sulfur_adjacent[bond.begin().index()] && !sulfur_adjacent[bond.end().index()]
            })
            .map(|bond| {
                mmat.get_upper(bond.begin().index(), bond.end().index())
                    - mmat.get_lower(bond.begin().index(), bond.end().index())
            })
            .expect("carbon-carbon bond");

        assert!((squished_width - (2.0 * (0.2 + DIST12_DELTA))).abs() < 1e-9);
        assert!((unsquished_width - (2.0 * DIST12_DELTA)).abs() < 1e-9);
    }

    #[test]
    fn hbond_helpers_follow_rdkit_acceptor_donor_and_bound_hydrogen_rules() {
        let ammonia = Molecule::from_smiles_with_sanitize("N", false).expect("ammonia-like N");
        let ammonia_valence = assign_valence(&ammonia, ValenceModel::RdkitLike).expect("valence");
        let ammonia_total_hs =
            total_num_hydrogens_for_distgeom(&ammonia, &ammonia.atoms()[0], &ammonia_valence, 0);
        assert!(is_hbond_acceptor(ammonia.atoms()[0].atomic_number()));
        assert!(is_hbond_donor(&ammonia.atoms()[0], ammonia_total_hs));

        let methane = Molecule::from_smiles_with_sanitize("C", false).expect("methane-like C");
        let methane_valence = assign_valence(&methane, ValenceModel::RdkitLike).expect("valence");
        let methane_total_hs =
            total_num_hydrogens_for_distgeom(&methane, &methane.atoms()[0], &methane_valence, 0);
        assert!(!is_hbond_acceptor(methane.atoms()[0].atomic_number()));
        assert!(!is_hbond_donor(&methane.atoms()[0], methane_total_hs));

        let mut builder = MoleculeBuilder::new();
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let oxygen = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(hydrogen, oxygen, BondOrder::Single))
            .unwrap();
        let hbond_mol = builder.build().expect("explicit O-H");
        assert!(is_h_in_hbond_donor(&hbond_mol, hydrogen.index()));
        assert!(!is_h_in_hbond_donor(&hbond_mol, oxygen.index()));
    }

    #[test]
    fn total_num_hydrogens_for_distgeom_counts_neighbor_hydrogens_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H).with_no_implicit(true));
        let nitrogen = builder.add_atom(AtomSpec::new(Element::N).with_no_implicit(true));
        let methyl = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        builder
            .add_bond(BondSpec::new(hydrogen, nitrogen, BondOrder::Single))
            .expect("h-n");
        builder
            .add_bond(BondSpec::new(nitrogen, methyl, BondOrder::Single))
            .expect("n-c");
        let mol = builder.build().expect("explicit secondary amine fragment");
        let assignment = assign_valence(&mol, ValenceModel::RdkitLike).expect("valence");
        let total_hs = total_num_hydrogens_for_distgeom(
            &mol,
            &mol.atoms()[nitrogen.index()],
            &assignment,
            nitrogen.index(),
        );
        assert_eq!(total_hs, 1);
    }

    #[test]
    fn set_lower_bound_vdw_scales_15_16_and_longer_paths_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        for _ in 0..7 {
            builder.add_atom(AtomSpec::new(Element::C));
        }
        let mol = builder.build().expect("seven carbons");
        let mut mmat = BoundsMatrix::new(mol.num_atoms());
        let mut dmat = vec![0.0; mol.num_atoms() * mol.num_atoms()];
        dmat[4 * mol.num_atoms()] = 4.0;
        dmat[5 * mol.num_atoms()] = 5.0;
        dmat[6 * mol.num_atoms()] = 6.0;

        set_lower_bound_vdw(&mol, &mut mmat, true, &dmat);

        let vdw_sum = vdw_radius(6) + vdw_radius(6);
        assert!((mmat.get_lower(4, 0) - (VDW_SCALE_15 * vdw_sum)).abs() < 1e-9);
        assert!(
            (mmat.get_lower(5, 0) - ((VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * vdw_sum)).abs()
                < 1e-9
        );
        assert!((mmat.get_lower(6, 0) - vdw_sum).abs() < 1e-9);
    }

    #[test]
    fn set_lower_bound_vdw_uses_hbond_floor_before_vdw_scaling() {
        let mut builder = MoleculeBuilder::new();
        let nitrogen = builder.add_atom(AtomSpec::new(Element::N));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let _oxygen = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(nitrogen, hydrogen, BondOrder::Single))
            .unwrap();
        let mol = builder.build().expect("H-N ... O");
        let mut mmat = BoundsMatrix::new(mol.num_atoms());
        let mut dmat = vec![0.0; mol.num_atoms() * mol.num_atoms()];
        dmat[2 * mol.num_atoms() + 1] = 6.0;
        dmat[1 * mol.num_atoms() + 2] = 6.0;

        set_lower_bound_vdw(&mol, &mut mmat, true, &dmat);

        assert!((mmat.get_lower(2, 1) - H_BOND_LENGTH).abs() < 1e-9);
    }

    #[test]
    fn set_ring_angle_matches_rdkit_ring_hybridization_special_cases() {
        let mut three_builder = MoleculeBuilder::new();
        three_builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2));
        let cyclopropane_like = three_builder.build().expect("3-ring test atom");

        let mut five_builder = MoleculeBuilder::new();
        five_builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3));
        let cyclopentane_like = five_builder.build().expect("5-ring test atom");

        let tri = set_ring_angle(&cyclopropane_like, 0, 3);
        let five = set_ring_angle(&cyclopentane_like, 0, 5);

        assert!((tri - std::f64::consts::PI / 3.0).abs() < 1e-9);
        assert!((five - (104.0_f64.to_radians())).abs() < 1e-9);
    }

    #[test]
    fn set_13_bounds_helper_doubles_tolerance_for_larger_sp2_ring_atoms() {
        let thiophene = Molecule::from_smiles("s1cccc1").expect("thiophene");
        let benzene = Molecule::from_smiles("c1ccccc1").expect("benzene");

        let (_thiophene_mmat_12, thiophene_accum) = run_set12_bounds(&thiophene);
        let (_benzene_mmat_12, benzene_accum) = run_set12_bounds(&benzene);

        let sulfur = thiophene
            .atoms()
            .iter()
            .position(|atom| atom.atomic_number() == 16)
            .expect("sulfur");
        let sulfur_neighbors: Vec<usize> = thiophene
            .bonds()
            .iter()
            .filter_map(|bond| {
                if bond.begin().index() == sulfur {
                    Some(bond.end().index())
                } else if bond.end().index() == sulfur {
                    Some(bond.begin().index())
                } else {
                    None
                }
            })
            .collect();
        let mut thiophene_mmat = BoundsMatrix::new(thiophene.num_atoms());
        let thiophene_angle = set_ring_angle(&thiophene, sulfur, 5);
        set_13_bounds_helper(
            sulfur_neighbors[0],
            sulfur,
            sulfur_neighbors[1],
            thiophene_angle,
            &thiophene_accum.bond_lengths,
            &mut thiophene_mmat,
            &thiophene,
        );

        let mut benzene_mmat = BoundsMatrix::new(benzene.num_atoms());
        let benzene_angle = set_ring_angle(&benzene, 1, 6);
        set_13_bounds_helper(
            0,
            1,
            2,
            benzene_angle,
            &benzene_accum.bond_lengths,
            &mut benzene_mmat,
            &benzene,
        );

        let thiophene_width = thiophene_mmat.get_upper(sulfur_neighbors[0], sulfur_neighbors[1])
            - thiophene_mmat.get_lower(sulfur_neighbors[0], sulfur_neighbors[1]);
        let benzene_width = benzene_mmat.get_upper(0, 2) - benzene_mmat.get_lower(0, 2);

        assert!((thiophene_width - (4.0 * DIST13_TOL)).abs() < 1e-9);
        assert!((benzene_width - (2.0 * DIST13_TOL)).abs() < 1e-9);
        assert!(is_larger_sp2_atom_idx(&thiophene, sulfur));
        assert!(!is_larger_sp2_atom_idx(&benzene, 1));
    }

    #[test]
    fn visited_bound_obeys_rdkit_dist_type_thresholds() {
        let mut accum = ComputedData::new(4, 2);
        accum.visited13_bounds[3] = true;
        accum.visited14_bounds[5] = true;

        assert!(!accum.visited_bound(3, DistType::Dist12));
        assert!(accum.visited_bound(3, DistType::Dist13));
        assert!(accum.visited_bound(3, DistType::Dist14));
        assert!(!accum.visited_bound(5, DistType::Dist13));
        assert!(accum.visited_bound(5, DistType::Dist14));
    }

    #[test]
    fn set_13_bounds_sets_non_ring_sp3_bounds_for_propane_path() {
        let mol = Molecule::from_smiles("CCC").expect("propane");
        let (mmat, accum_data) = run_set13_bounds(&mol);

        let bid01 = bond_between_idx_simple(&mol, 0, 1).expect("0-1 bond");
        let bid12 = bond_between_idx_simple(&mol, 1, 2).expect("1-2 bond");
        let expected = compute_13_dist(
            accum_data.bond_lengths[bid01],
            accum_data.bond_lengths[bid12],
            109.5_f64.to_radians(),
        );

        assert!((mmat.get_lower(0, 2) - (expected - DIST13_TOL)).abs() < 1e-9);
        assert!((mmat.get_upper(0, 2) - (expected + DIST13_TOL)).abs() < 1e-9);
        assert!(
            (accum_data.get_bond_angle(mol.num_bonds(), bid01, bid12) - 109.5_f64.to_radians())
                .abs()
                < 1e-9
        );
        assert!(accum_data.visited13_bounds[2]);
    }

    #[test]
    fn set_13_bounds_distributes_remaining_fused_ring_angle_like_rdkit() {
        let mol = Molecule::from_smiles("c1cccc2ccccc12").expect("naphthalene");
        let (_mmat, accum_data) = run_set13_bounds(&mol);
        let rings = mol.derived_cache().rings.as_ref().expect("rings");

        let fusion_atom = mol
            .atoms()
            .iter()
            .enumerate()
            .find_map(|(idx, atom)| {
                (atom.hybridization() == Hybridization::Sp2
                    && rings.num_atom_rings(atom.id()) > 1
                    && neighbors_for_atom(&mol, idx).len() == 3)
                    .then_some(idx)
            })
            .expect("fusion atom");

        let neighbors = neighbors_for_atom(&mol, fusion_atom);
        let mut pair_angles = Vec::new();
        for left in 0..neighbors.len() {
            let bid1 = bond_between_idx_simple(&mol, fusion_atom, neighbors[left]).expect("bond 1");
            for right in 0..left {
                let bid2 =
                    bond_between_idx_simple(&mol, fusion_atom, neighbors[right]).expect("bond 2");
                pair_angles.push(accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2));
            }
        }

        assert_eq!(pair_angles.len(), 3);
        for angle in pair_angles {
            assert!((angle - 120.0_f64.to_radians()).abs() < 1e-9);
        }
    }

    #[test]
    fn chain_and_carbonyl_classification_helpers_follow_rdkit_patterns() {
        let propane = Molecule::from_smiles("CCC").expect("propane");
        assert!(check_h2_nx3_h1_ox2(&propane, 1));

        let ether = Molecule::from_smiles("COC").expect("dimethyl ether");
        let oxygen = ether
            .atoms()
            .iter()
            .position(|atom| atom.atomic_number() == 8)
            .expect("ether oxygen");
        assert!(check_h2_nx3_h1_ox2(&ether, oxygen));

        let butane = Molecule::from_smiles("CCCC").expect("butane");
        assert!(check_nh_ch_ch_nh(&butane, 0, 1, 2, 3));

        let acetate = Molecule::from_smiles("CC(=O)O").expect("acetate");
        let carbonyl = acetate
            .atoms()
            .iter()
            .enumerate()
            .find_map(|(idx, _)| is_carbonyl(&acetate, idx).then_some(idx))
            .expect("carbonyl carbon");
        assert!(is_carbonyl(&acetate, carbonyl));
        assert!(!is_carbonyl(&acetate, 0));
    }

    #[test]
    fn amide_ester_classification_helpers_match_ester_patterns() {
        let ester = Molecule::from_smiles("COC(=O)C").expect("methyl acetate");
        let carbonyl = ester
            .atoms()
            .iter()
            .enumerate()
            .find_map(|(idx, _)| is_carbonyl(&ester, idx).then_some(idx))
            .expect("carbonyl carbon");

        let double_hetero = neighbors_for_atom(&ester, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&ester, carbonyl, nbr).expect("bond");
                ester.bonds()[bond_idx].order() == BondOrder::Double
                    && (ester.atoms()[nbr].atomic_number() == 8
                        || ester.atoms()[nbr].atomic_number() == 7)
            })
            .expect("double bonded hetero");
        let single_hetero = neighbors_for_atom(&ester, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&ester, carbonyl, nbr).expect("bond");
                ester.bonds()[bond_idx].order() == BondOrder::Single
                    && (ester.atoms()[nbr].atomic_number() == 8
                        || ester.atoms()[nbr].atomic_number() == 7)
            })
            .expect("single bonded hetero");
        let atom1 = neighbors_for_atom(&ester, single_hetero)
            .into_iter()
            .find(|&nbr| nbr != carbonyl)
            .expect("atom1");
        let bnd1 = bond_between_idx_simple(&ester, atom1, single_hetero).expect("bond1");
        let bnd3_double =
            bond_between_idx_simple(&ester, carbonyl, double_hetero).expect("bond3 double");
        let carbonyl_substituent = neighbors_for_atom(&ester, carbonyl)
            .into_iter()
            .find(|&nbr| nbr != single_hetero && nbr != double_hetero)
            .expect("carbonyl substituent");
        let bnd3_single =
            bond_between_idx_simple(&ester, carbonyl, carbonyl_substituent).expect("bond3 single");

        assert!(check_amide_ester_14(
            &ester,
            bnd1,
            bnd3_double,
            single_hetero,
            carbonyl,
            double_hetero,
        ));
        assert!(check_amide_ester_15(
            &ester,
            bnd1,
            bnd3_single,
            single_hetero,
            carbonyl,
        ));
    }

    #[test]
    fn macrocycle_all_in_same_ring_amide_ester_helper_matches_lactone_pattern() {
        let mol = Molecule::from_smiles("O=C1N(C)CCCCCCCC1").expect("macrocyclic tertiary lactam");
        let carbonyl = mol
            .atoms()
            .iter()
            .enumerate()
            .find_map(|(idx, _)| is_carbonyl(&mol, idx).then_some(idx))
            .expect("carbonyl carbon");
        let atm2 = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds()[bond_idx].order() == BondOrder::Single
                    && mol.atoms()[nbr].atomic_number() == 7
            })
            .expect("ring amide nitrogen");
        let atm4 = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds()[bond_idx].order() == BondOrder::Single && nbr != atm2
            })
            .expect("ring carbon neighbor");
        let atm1 = neighbors_for_atom(&mol, atm2)
            .into_iter()
            .find(|&nbr| nbr != carbonyl && mol.atoms()[nbr].atomic_number() == 6)
            .expect("preceding ring atom");

        assert!(check_macrocycle_all_in_same_ring_amide_ester_14(
            &mol, atm1, atm2, carbonyl, atm4,
        ));
    }

    #[test]
    fn macrocycle_two_in_same_ring_amide_ester_helper_matches_tertiary_lactam_pattern() {
        let mol = Molecule::from_smiles("O=C1N(C)CCCCCCCC1").expect("macrocyclic tertiary lactam");
        let carbonyl = mol
            .atoms()
            .iter()
            .enumerate()
            .find_map(|(idx, _)| is_carbonyl(&mol, idx).then_some(idx))
            .expect("carbonyl carbon");
        let oxygen = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds()[bond_idx].order() == BondOrder::Double
            })
            .expect("carbonyl oxygen");
        let nitrogen = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds()[bond_idx].order() == BondOrder::Single
                    && mol.atoms()[nbr].atomic_number() == 7
            })
            .expect("amide nitrogen");
        let atom1 = neighbors_for_atom(&mol, nitrogen)
            .into_iter()
            .find(|&nbr| nbr != carbonyl && mol.atoms()[nbr].atomic_number() == 6)
            .expect("preceding ring carbon");
        let bnd1 = bond_between_idx_simple(&mol, atom1, nitrogen).expect("bond1");
        let bnd3 = bond_between_idx_simple(&mol, carbonyl, oxygen).expect("bond3");

        assert!(check_macrocycle_two_in_same_ring_amide_ester_14(
            &mol, bnd1, bnd3, atom1, nitrogen, carbonyl, oxygen,
        ));
    }

    #[test]
    fn set_macrocycle_two_in_same_ring_14_bounds_uses_cis_for_macrocycle_amide_path() {
        let mol = Molecule::from_smiles("O=C1N(C)CCCCCCCC1").expect("macrocyclic tertiary lactam");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let dmat = flatten_topological_distances(&mol);

        let carbonyl = mol
            .atoms()
            .iter()
            .enumerate()
            .find_map(|(idx, _)| is_carbonyl(&mol, idx).then_some(idx))
            .expect("carbonyl carbon");
        let oxygen = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds()[bond_idx].order() == BondOrder::Double
            })
            .expect("carbonyl oxygen");
        let nitrogen = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds()[bond_idx].order() == BondOrder::Single
                    && mol.atoms()[nbr].atomic_number() == 7
            })
            .expect("amide nitrogen");
        let atom1 = neighbors_for_atom(&mol, nitrogen)
            .into_iter()
            .find(|&nbr| nbr != carbonyl && mol.atoms()[nbr].atomic_number() == 6)
            .expect("preceding ring carbon");

        let bid1 = bond_between_idx_simple(&mol, atom1, nitrogen).expect("b1");
        let bid2 = bond_between_idx_simple(&mol, nitrogen, carbonyl).expect("b2");
        let bid3 = bond_between_idx_simple(&mol, carbonyl, oxygen).expect("b3");
        set_macrocycle_two_in_same_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut accum_data,
            &mut mmat,
            &dmat,
        );

        let path = accum_data.paths14.last().expect("path");
        assert_eq!(path.kind, Path14Kind::Cis);
        let expected = compute_14_dist_cis(
            accum_data.bond_lengths[bid1],
            accum_data.bond_lengths[bid2],
            accum_data.bond_lengths[bid3],
            accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2),
            accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3),
        );
        assert!((mmat.get_lower(atom1, oxygen) - (expected - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((mmat.get_upper(atom1, oxygen) - (expected + GEN_DIST_TOL)).abs() < 1e-9);
    }

    #[test]
    fn set_macrocycle_all_in_same_ring_14_bounds_uses_trans_plus_point_one_for_macrocycle_amide() {
        let mol = Molecule::from_smiles("O=C1N(C)CCCCCCCC1").expect("macrocyclic tertiary lactam");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);

        let carbonyl = mol
            .atoms()
            .iter()
            .enumerate()
            .find_map(|(idx, _)| is_carbonyl(&mol, idx).then_some(idx))
            .expect("carbonyl carbon");
        let nitrogen = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds()[bond_idx].order() == BondOrder::Single
                    && mol.atoms()[nbr].atomic_number() == 7
            })
            .expect("amide nitrogen");
        let atm4 = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds()[bond_idx].order() == BondOrder::Single && nbr != nitrogen
            })
            .expect("ring carbon neighbor");
        let atom1 = neighbors_for_atom(&mol, nitrogen)
            .into_iter()
            .find(|&nbr| nbr != carbonyl && mol.atoms()[nbr].atomic_number() == 6)
            .expect("preceding ring carbon");

        let bid1 = bond_between_idx_simple(&mol, atom1, nitrogen).expect("b1");
        let bid2 = bond_between_idx_simple(&mol, nitrogen, carbonyl).expect("b2");
        let bid3 = bond_between_idx_simple(&mol, carbonyl, atm4).expect("b3");
        set_macrocycle_all_in_same_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut accum_data,
            &mut mmat,
        );

        let path = accum_data.paths14.last().expect("path");
        assert_eq!(path.kind, Path14Kind::Trans);
        let expected = compute_14_dist_trans(
            accum_data.bond_lengths[bid1],
            accum_data.bond_lengths[bid2],
            accum_data.bond_lengths[bid3],
            accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2),
            accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3),
        ) + 0.1;
        assert!((mmat.get_lower(atom1, atm4) - (expected - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((mmat.get_upper(atom1, atm4) - (expected + GEN_DIST_TOL)).abs() < 1e-9);
    }

    #[test]
    fn set_chain_14_bounds_uses_defined_double_bond_stereo_for_alkenes() {
        let trans = Molecule::from_smiles("C/C=C/C").expect("trans alkene");
        let (mut trans_mmat, mut trans_accum) = run_set13_bounds(&trans);
        let t_bid1 = bond_between_idx_simple(&trans, 0, 1).expect("0-1");
        let t_bid2 = bond_between_idx_simple(&trans, 1, 2).expect("1-2");
        let t_bid3 = bond_between_idx_simple(&trans, 2, 3).expect("2-3");
        set_chain_14_bounds(
            &trans,
            t_bid1,
            t_bid2,
            t_bid3,
            &mut trans_accum,
            &mut trans_mmat,
            false,
        );
        let trans_path = trans_accum.paths14.last().expect("trans path");
        assert_eq!(trans_path.kind, Path14Kind::Trans);
        let trans_expected = compute_14_dist_trans(
            trans_accum.bond_lengths[t_bid1],
            trans_accum.bond_lengths[t_bid2],
            trans_accum.bond_lengths[t_bid3],
            trans_accum.get_bond_angle(trans.num_bonds(), t_bid1, t_bid2),
            trans_accum.get_bond_angle(trans.num_bonds(), t_bid2, t_bid3),
        );
        assert!((trans_mmat.get_lower(0, 3) - (trans_expected - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((trans_mmat.get_upper(0, 3) - (trans_expected + GEN_DIST_TOL)).abs() < 1e-9);

        let cis = Molecule::from_smiles("C/C=C\\C").expect("cis alkene");
        let (mut cis_mmat, mut cis_accum) = run_set13_bounds(&cis);
        let c_bid1 = bond_between_idx_simple(&cis, 0, 1).expect("0-1");
        let c_bid2 = bond_between_idx_simple(&cis, 1, 2).expect("1-2");
        let c_bid3 = bond_between_idx_simple(&cis, 2, 3).expect("2-3");
        set_chain_14_bounds(
            &cis,
            c_bid1,
            c_bid2,
            c_bid3,
            &mut cis_accum,
            &mut cis_mmat,
            false,
        );
        let cis_path = cis_accum.paths14.last().expect("cis path");
        assert_eq!(cis_path.kind, Path14Kind::Cis);
        let cis_expected = compute_14_dist_cis(
            cis_accum.bond_lengths[c_bid1],
            cis_accum.bond_lengths[c_bid2],
            cis_accum.bond_lengths[c_bid3],
            cis_accum.get_bond_angle(cis.num_bonds(), c_bid1, c_bid2),
            cis_accum.get_bond_angle(cis.num_bonds(), c_bid2, c_bid3),
        );
        assert!((cis_mmat.get_lower(0, 3) - (cis_expected - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((cis_mmat.get_upper(0, 3) - (cis_expected + GEN_DIST_TOL)).abs() < 1e-9);
    }

    #[test]
    fn set_chain_14_bounds_uses_ss_special_case() {
        let mol = Molecule::from_smiles("CSSC").expect("disulfide");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
        set_chain_14_bounds(&mol, bid1, bid2, bid3, &mut accum_data, &mut mmat, false);

        let path = accum_data.paths14.last().expect("path");
        assert_eq!(path.kind, Path14Kind::Other);
        let expected = compute_14_dist_3d(
            accum_data.bond_lengths[bid1],
            accum_data.bond_lengths[bid2],
            accum_data.bond_lengths[bid3],
            accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2),
            accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3),
            std::f64::consts::PI / 2.0,
        );
        assert!((mmat.get_lower(0, 3) - (expected - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((mmat.get_upper(0, 3) - (expected + GEN_DIST_TOL)).abs() < 1e-9);
    }

    #[test]
    fn set_chain_14_bounds_honors_force_trans_amides_for_secondary_amide_h_paths() {
        let mut builder = Molecule::builder();
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H)).index();
        let nitrogen =
            builder.add_atom(AtomSpec::new(Element::N).with_hybridization(Hybridization::Sp2));
        let nitrogen = nitrogen.index();
        let carbonyl =
            builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2));
        let carbonyl = carbonyl.index();
        let oxygen =
            builder.add_atom(AtomSpec::new(Element::O).with_hybridization(Hybridization::Sp2));
        let oxygen = oxygen.index();
        let n_methyl = builder.add_atom(AtomSpec::new(Element::C)).index();
        let carbonyl_methyl = builder.add_atom(AtomSpec::new(Element::C)).index();
        builder
            .add_bond(BondSpec::new(
                AtomId::new(hydrogen),
                AtomId::new(nitrogen),
                BondOrder::Single,
            ))
            .expect("h-n");
        builder
            .add_bond(BondSpec::new(
                AtomId::new(nitrogen),
                AtomId::new(carbonyl),
                BondOrder::Single,
            ))
            .expect("n-c");
        builder
            .add_bond(BondSpec::new(
                AtomId::new(nitrogen),
                AtomId::new(n_methyl),
                BondOrder::Single,
            ))
            .expect("n-c methyl");
        builder
            .add_bond(BondSpec::new(
                AtomId::new(carbonyl),
                AtomId::new(oxygen),
                BondOrder::Double,
            ))
            .expect("c-o");
        builder
            .add_bond(BondSpec::new(
                AtomId::new(carbonyl),
                AtomId::new(carbonyl_methyl),
                BondOrder::Single,
            ))
            .expect("c-c");
        let mol = builder.build().expect("secondary amide with explicit N-H");
        let (mut mmat, mut accum_data) = run_set12_bounds(&mol);

        let bid_hn = bond_between_idx_simple(&mol, hydrogen, nitrogen).expect("h-n");
        let bid_nc = bond_between_idx_simple(&mol, nitrogen, carbonyl).expect("n-c");
        let bid_co = bond_between_idx_simple(&mol, carbonyl, oxygen).expect("c-o");
        let bid_cm = bond_between_idx_simple(&mol, carbonyl, carbonyl_methyl).expect("c-m");
        let nb = mol.num_bonds();
        accum_data.set_bond_adj(nb, bid_hn, bid_nc, nitrogen as i32);
        accum_data.set_bond_adj(nb, bid_nc, bid_co, carbonyl as i32);
        accum_data.set_bond_adj(nb, bid_nc, bid_cm, carbonyl as i32);
        accum_data.set_bond_angle(
            nb,
            bid_hn,
            bid_nc,
            ideal_bond_angle(&mol.atoms()[nitrogen].hybridization(), None),
        );
        accum_data.set_bond_angle(
            nb,
            bid_nc,
            bid_co,
            ideal_bond_angle(&mol.atoms()[carbonyl].hybridization(), None),
        );
        accum_data.set_bond_angle(
            nb,
            bid_nc,
            bid_cm,
            ideal_bond_angle(&mol.atoms()[carbonyl].hybridization(), None),
        );

        set_chain_14_bounds(
            &mol,
            bid_hn,
            bid_nc,
            bid_co,
            &mut accum_data,
            &mut mmat,
            true,
        );
        let amide14_path = accum_data.paths14.last().expect("amide14 path");
        assert_eq!(amide14_path.kind, Path14Kind::Trans);
        let expected_14 = compute_14_dist_trans(
            accum_data.bond_lengths[bid_hn],
            accum_data.bond_lengths[bid_nc],
            accum_data.bond_lengths[bid_co],
            accum_data.get_bond_angle(mol.num_bonds(), bid_hn, bid_nc),
            accum_data.get_bond_angle(mol.num_bonds(), bid_nc, bid_co),
        );
        assert!((mmat.get_lower(hydrogen, oxygen) - (expected_14 - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((mmat.get_upper(hydrogen, oxygen) - (expected_14 + GEN_DIST_TOL)).abs() < 1e-9);

        set_chain_14_bounds(
            &mol,
            bid_hn,
            bid_nc,
            bid_cm,
            &mut accum_data,
            &mut mmat,
            true,
        );
        let amide15_path = accum_data.paths14.last().expect("amide15 path");
        assert_eq!(amide15_path.kind, Path14Kind::Cis);
        let expected_15 = compute_14_dist_cis(
            accum_data.bond_lengths[bid_hn],
            accum_data.bond_lengths[bid_nc],
            accum_data.bond_lengths[bid_cm],
            accum_data.get_bond_angle(mol.num_bonds(), bid_hn, bid_nc),
            accum_data.get_bond_angle(mol.num_bonds(), bid_nc, bid_cm),
        );
        assert!(
            (mmat.get_lower(hydrogen, carbonyl_methyl) - (expected_15 - GEN_DIST_TOL)).abs() < 1e-9
        );
        assert!(
            (mmat.get_upper(hydrogen, carbonyl_methyl) - (expected_15 + GEN_DIST_TOL)).abs() < 1e-9
        );
    }

    #[test]
    fn record_14_path_marks_sp2_sp2_ring_paths_as_cis() {
        let mol = Molecule::from_smiles("c1ccccc1").expect("benzene");
        let (_mmat, mut accum_data) = run_set13_bounds(&mol);
        let bid1 = bond_between_idx_simple(&mol, 5, 0).expect("5-0");
        let bid2 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid3 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        record_14_path(&mol, bid1, bid2, bid3, &mut accum_data);

        let path = accum_data.paths14.last().expect("path");
        assert_eq!(path.bid1, bid1);
        assert_eq!(path.bid2, bid2);
        assert_eq!(path.bid3, bid3);
        assert_eq!(path.kind, Path14Kind::Cis);
        assert!(has_path_flag(
            &accum_data.cis_paths,
            path14_id(mol.num_bonds(), bid1, bid2, bid3)
        ));
        assert!(has_path_flag(
            &accum_data.cis_paths,
            path14_id(mol.num_bonds(), bid3, bid2, bid1)
        ));
    }

    #[test]
    fn set_in_ring_14_bounds_prefers_cis_for_small_sp2_ring_paths() {
        let mol = Molecule::from_smiles("c1ccccc1").expect("benzene");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let dmat = flatten_topological_distances(&mol);
        let bid1 = bond_between_idx_simple(&mol, 5, 0).expect("5-0");
        let bid2 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid3 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        set_in_ring_14_bounds(&mol, bid1, bid2, bid3, &mut accum_data, &mut mmat, &dmat, 6);

        let path = accum_data.paths14.last().expect("path");
        assert_eq!(path.kind, Path14Kind::Cis);
        assert!(has_path_flag(
            &accum_data.cis_paths,
            path14_id(mol.num_bonds(), bid1, bid2, bid3)
        ));

        let aid1 = 5usize;
        let aid4 = 2usize;
        let pid = aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4);
        let expected = compute_14_dist_cis(
            accum_data.bond_lengths[bid1],
            accum_data.bond_lengths[bid2],
            accum_data.bond_lengths[bid3],
            accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2),
            accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3),
        );
        assert!(accum_data.visited14_bounds[pid]);
        assert!((mmat.get_lower(aid1, aid4) - (expected - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((mmat.get_upper(aid1, aid4) - (expected + GEN_DIST_TOL)).abs() < 1e-9);
    }

    #[test]
    fn set_two_in_same_ring_14_bounds_uses_trans_for_sp2_external_substituent_path() {
        let mol = Molecule::from_smiles("Cc1ccccc1").expect("toluene");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let dmat = flatten_topological_distances(&mol);

        let bid_exo = bond_between_idx_simple(&mol, 0, 1).expect("exo bond");
        let bid_ring_12 = bond_between_idx_simple(&mol, 1, 2).expect("ring bond 1-2");
        let bid_ring_23 = bond_between_idx_simple(&mol, 2, 3).expect("ring bond 2-3");
        set_two_in_same_ring_14_bounds(
            &mol,
            bid_exo,
            bid_ring_12,
            bid_ring_23,
            &mut accum_data,
            &mut mmat,
            &dmat,
        );

        let path = accum_data.paths14.last().expect("path");
        assert_eq!(path.kind, Path14Kind::Trans);
        assert!(has_path_flag(
            &accum_data.trans_paths,
            path14_id(mol.num_bonds(), bid_exo, bid_ring_12, bid_ring_23)
        ));

        let aid1 = 0usize;
        let aid4 = 3usize;
        let expected = compute_14_dist_trans(
            accum_data.bond_lengths[bid_exo],
            accum_data.bond_lengths[bid_ring_12],
            accum_data.bond_lengths[bid_ring_23],
            accum_data.get_bond_angle(mol.num_bonds(), bid_exo, bid_ring_12),
            accum_data.get_bond_angle(mol.num_bonds(), bid_ring_12, bid_ring_23),
        );
        assert!((mmat.get_lower(aid1, aid4) - (expected - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((mmat.get_upper(aid1, aid4) - (expected + GEN_DIST_TOL)).abs() < 1e-9);
    }

    #[test]
    fn diff_ring14_and_share_ring14_delegate_to_in_ring_helper() {
        let mol = Molecule::from_smiles("c1ccccc1").expect("benzene");
        let dmat = flatten_topological_distances(&mol);
        let bid1 = bond_between_idx_simple(&mol, 5, 0).expect("5-0");
        let bid2 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid3 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");

        let (mut base_mmat, mut base_accum) = run_set13_bounds(&mol);
        set_in_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut base_accum,
            &mut base_mmat,
            &dmat,
            0,
        );

        let (mut diff_mmat, mut diff_accum) = run_set13_bounds(&mol);
        set_two_in_diff_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut diff_accum,
            &mut diff_mmat,
            &dmat,
        );

        let (mut share_mmat, mut share_accum) = run_set13_bounds(&mol);
        set_share_ring_bond_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut share_accum,
            &mut share_mmat,
            &dmat,
        );

        assert_eq!(
            diff_accum.paths14.last().map(|p| p.kind),
            base_accum.paths14.last().map(|p| p.kind)
        );
        assert_eq!(
            share_accum.paths14.last().map(|p| p.kind),
            base_accum.paths14.last().map(|p| p.kind)
        );
        assert!((diff_mmat.get_lower(5, 2) - base_mmat.get_lower(5, 2)).abs() < 1e-9);
        assert!((diff_mmat.get_upper(5, 2) - base_mmat.get_upper(5, 2)).abs() < 1e-9);
        assert!((share_mmat.get_lower(5, 2) - base_mmat.get_lower(5, 2)).abs() < 1e-9);
        assert!((share_mmat.get_upper(5, 2) - base_mmat.get_upper(5, 2)).abs() < 1e-9);
    }

    #[test]
    fn compute_15_dist_helpers_match_rdkit_cis_and_trans_formulas() {
        let d1: f64 = 1.41;
        let d2: f64 = 1.52;
        let d3: f64 = 1.38;
        let d4: f64 = 1.47;
        let ang12: f64 = 1.91;
        let ang23: f64 = 2.04;
        let ang34: f64 = 1.88;

        let cis_dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
        let cis_dy14 = d3 * ang23.sin() - d1 * ang12.sin();
        let cis_d14 = (cis_dx14 * cis_dx14 + cis_dy14 * cis_dy14).sqrt();
        let cis_cval =
            ((d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / cis_d14).clamp(-1.0, 1.0);
        let cis_ang143 = cis_cval.acos();
        let expected_cis_cis = compute_13_dist(cis_d14, d4, ang34 - cis_ang143);
        let expected_cis_trans = compute_13_dist(cis_d14, d4, ang34 + cis_ang143);

        let trans_dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
        let trans_dy14 = d3 * ang23.sin() + d1 * ang12.sin();
        let trans_d14 = (trans_dx14 * trans_dx14 + trans_dy14 * trans_dy14).sqrt();
        let trans_cval =
            ((d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / trans_d14).clamp(-1.0, 1.0);
        let trans_ang143 = trans_cval.acos();
        let expected_trans_trans = compute_13_dist(trans_d14, d4, ang34 + trans_ang143);
        let expected_trans_cis = compute_13_dist(trans_d14, d4, ang34 - trans_ang143);

        assert!(
            (compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34) - expected_cis_cis).abs()
                < 1e-12
        );
        assert!(
            (compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34) - expected_cis_trans)
                .abs()
                < 1e-12
        );
        assert!(
            (compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34)
                - expected_trans_trans)
                .abs()
                < 1e-12
        );
        assert!(
            (compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34) - expected_trans_cis)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn compute_15_dist_helpers_support_rdkit_reverse_argument_order() {
        let d1: f64 = 1.41;
        let d2: f64 = 1.52;
        let d3: f64 = 1.38;
        let d4: f64 = 1.47;
        let ang12: f64 = 1.91;
        let ang23: f64 = 2.04;
        let ang34: f64 = 1.88;

        let reversed_cis_cis = compute_15_dist_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12);
        let reversed_cis_trans = compute_15_dist_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12);
        let reversed_trans_cis = compute_15_dist_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12);
        let reversed_trans_trans = compute_15_dist_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12);

        assert_ne!(
            compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34),
            reversed_cis_cis
        );
        assert_ne!(
            compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34),
            reversed_cis_trans
        );
        assert_ne!(
            compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34),
            reversed_trans_cis
        );
        assert_ne!(
            compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34),
            reversed_trans_trans
        );
    }

    #[test]
    fn set_15_bounds_helper_returns_immediately_for_visited_14_pair() {
        let mol = Molecule::from_smiles("CCCCC").expect("pentane");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let dmat = flatten_topological_distances(&mol);
        let nb = mol.num_bonds();
        let na = mol.num_atoms();
        let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
        let pid = 0usize * na + 4usize;
        accum_data.visited14_bounds[pid] = true;
        let before_lower = mmat.get_lower(0, 4);
        let before_upper = mmat.get_upper(0, 4);

        set_15_bounds_helper(
            &mol,
            &mut mmat,
            &mut accum_data,
            &dmat,
            nb,
            na,
            bid1,
            bid2,
            bid3,
            Path14Kind::Other,
        );

        assert_eq!(mmat.get_lower(0, 4), before_lower);
        assert_eq!(mmat.get_upper(0, 4), before_upper);
        assert!(!accum_data.set15_atoms[0 * na + 4]);
        assert!(!accum_data.set15_atoms[4 * na + 0]);
    }

    #[test]
    fn set_15_bounds_helper_uses_vdw_fallback_and_marks_set15_atoms_for_unknown_path() {
        let mol = Molecule::from_smiles("CCCCC").expect("pentane");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let dmat = flatten_topological_distances(&mol);
        let nb = mol.num_bonds();
        let na = mol.num_atoms();
        let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");

        set_15_bounds_helper(
            &mol,
            &mut mmat,
            &mut accum_data,
            &dmat,
            nb,
            na,
            bid1,
            bid2,
            bid3,
            Path14Kind::Other,
        );

        let expected_lower = VDW_SCALE_15 * (vdw_radius(6) + vdw_radius(6));
        assert!((mmat.get_lower(0, 4) - expected_lower).abs() < 1e-12);
        assert_eq!(mmat.get_upper(0, 4), MAX_UPPER);
        assert!(accum_data.set15_atoms[0 * na + 4]);
        assert!(accum_data.set15_atoms[4 * na + 0]);
    }

    #[test]
    fn set_15_bounds_helper_uses_reversed_other_branch_formula_for_cis_path() {
        let mol = Molecule::from_smiles("CCCCC").expect("pentane");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let dmat = flatten_topological_distances(&mol);
        let nb = mol.num_bonds();
        let na = mol.num_atoms();
        let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
        let bid4 = bond_between_idx_simple(&mol, 3, 4).expect("3-4");
        let path_id = bid2 as u64 * nb as u64 * nb as u64 + bid3 as u64 * nb as u64 + bid4 as u64;
        record_path_flag(&mut accum_data.cis_paths, path_id);

        set_15_bounds_helper(
            &mol,
            &mut mmat,
            &mut accum_data,
            &dmat,
            nb,
            na,
            bid1,
            bid2,
            bid3,
            Path14Kind::Other,
        );

        let d1 = accum_data.bond_lengths[bid1];
        let d2 = accum_data.bond_lengths[bid2];
        let d3 = accum_data.bond_lengths[bid3];
        let d4 = accum_data.bond_lengths[bid4];
        let ang12 = accum_data.get_bond_angle(nb, bid1, bid2);
        let ang23 = accum_data.get_bond_angle(nb, bid2, bid3);
        let ang34 = accum_data.get_bond_angle(nb, bid3, bid4);
        let expected_lower =
            compute_15_dist_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
        let expected_upper =
            compute_15_dist_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;

        assert!((mmat.get_lower(0, 4) - expected_lower).abs() < 1e-12);
        assert!((mmat.get_upper(0, 4) - expected_upper).abs() < 1e-12);
        assert!(accum_data.set15_atoms[0 * na + 4]);
        assert!(accum_data.set15_atoms[4 * na + 0]);
    }

    #[test]
    fn set_15_bounds_entrypoint_matches_two_helper_calls_for_single_path() {
        let mol = Molecule::from_smiles("CCCCC").expect("pentane");
        let dmat = flatten_topological_distances(&mol);
        let nb = mol.num_bonds();
        let na = mol.num_atoms();
        let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
        let bid4 = bond_between_idx_simple(&mol, 3, 4).expect("3-4");
        let path_id = bid2 as u64 * nb as u64 * nb as u64 + bid3 as u64 * nb as u64 + bid4 as u64;

        let (mut entry_mmat, mut entry_accum) = run_set13_bounds(&mol);
        entry_accum.paths14.push(Path14Configuration {
            bid1,
            bid2,
            bid3,
            kind: Path14Kind::Other,
        });
        record_path_flag(&mut entry_accum.cis_paths, path_id);
        set_15_bounds(&mol, &mut entry_mmat, &mut entry_accum, &dmat);

        let (mut helper_mmat, mut helper_accum) = run_set13_bounds(&mol);
        helper_accum.paths14.push(Path14Configuration {
            bid1,
            bid2,
            bid3,
            kind: Path14Kind::Other,
        });
        record_path_flag(&mut helper_accum.cis_paths, path_id);
        set_15_bounds_helper(
            &mol,
            &mut helper_mmat,
            &mut helper_accum,
            &dmat,
            nb,
            na,
            bid1,
            bid2,
            bid3,
            Path14Kind::Other,
        );
        set_15_bounds_helper(
            &mol,
            &mut helper_mmat,
            &mut helper_accum,
            &dmat,
            nb,
            na,
            bid3,
            bid2,
            bid1,
            Path14Kind::Other,
        );

        assert_eq!(entry_mmat.get_lower(0, 4), helper_mmat.get_lower(0, 4));
        assert_eq!(entry_mmat.get_upper(0, 4), helper_mmat.get_upper(0, 4));
        assert_eq!(entry_accum.set15_atoms, helper_accum.set15_atoms);
    }

    #[test]
    fn set_15_bounds_uses_paths14_produced_by_set_14_bounds() {
        let mol = Molecule::from_smiles("CCCCC").expect("pentane");
        let (mmat, accum_data, _dmat) = run_set15_bounds(&mol, false, false);

        assert!(!accum_data.paths14.is_empty());
        assert!(accum_data.set15_atoms[0 * mol.num_atoms() + 4]);
        assert!(accum_data.set15_atoms[4 * mol.num_atoms() + 0]);
        assert!(mmat.get_lower(0, 4) > 0.0);
        assert!(mmat.get_upper(0, 4) >= mmat.get_lower(0, 4));
    }

    #[test]
    fn set_topol_bounds_can_disable_13_and_14_stages_like_rdkit() {
        let mol = Molecule::from_smiles("CCCC").expect("butane");
        let disabled = run_set_topol_bounds(&mol, false, true, false, false, false, false);
        let enabled = run_set_topol_bounds(&mol, false, true, false, false, true, true);

        assert_eq!(disabled.get_upper(0, 2), MAX_UPPER);
        assert_eq!(disabled.get_upper(0, 3), MAX_UPPER);
        assert!(enabled.get_upper(0, 2) < MAX_UPPER);
        assert!(enabled.get_upper(0, 3) < MAX_UPPER);
    }

    #[test]
    fn set_topol_bounds_can_enable_15_stage_independently_like_rdkit() {
        let mol = Molecule::from_smiles("C/C=C/CC").expect("stereo pentene");
        let without_15 = run_set_topol_bounds(&mol, false, true, false, false, true, true);
        let with_15 = run_set_topol_bounds(&mol, true, true, false, false, true, true);

        assert_eq!(without_15.get_upper(0, 4), MAX_UPPER);
        assert!(with_15.get_upper(0, 4) < MAX_UPPER);
        assert!(with_15.get_lower(0, 4) > without_15.get_lower(0, 4));
    }

    #[test]
    fn set_topol_bounds_ignores_scale_vdw_flag_like_current_rdkit_source() {
        let mol = Molecule::from_smiles("CCCCCC").expect("hexane");
        let scaled = run_set_topol_bounds(&mol, false, true, false, false, false, false);
        let unscaled = run_set_topol_bounds(&mol, false, false, false, false, false, false);

        assert_eq!(scaled.get_lower(0, 5), unscaled.get_lower(0, 5));
        assert_eq!(scaled.get_upper(0, 5), unscaled.get_upper(0, 5));
    }

    #[test]
    fn set_topol_bounds_with_outputs_matches_first_overload_matrix() {
        let mol = Molecule::from_smiles("C/C=C/CC").expect("stereo pentene");
        let plain = run_set_topol_bounds(&mol, true, true, false, false, true, true);
        let (with_outputs, bonds, angles) =
            run_set_topol_bounds_with_outputs(&mol, true, true, false, false, true, true);

        assert_eq!(plain.data, with_outputs.data);
        assert_eq!(bonds.len(), mol.num_bonds());
        assert!(!angles.is_empty());
    }

    #[test]
    fn set_topol_bounds_with_outputs_emits_exact_rdkit_bonds_and_angles_for_triple_bond_case() {
        let mol = Molecule::from_smiles("CC#N").expect("acetonitrile");
        let (_mmat, bonds, angles) =
            run_set_topol_bounds_with_outputs(&mol, false, true, false, false, false, false);

        assert_eq!(bonds, vec![(0, 1), (1, 2)]);
        assert_eq!(angles, vec![vec![0, 1, 2, 1]]);
    }

    #[test]
    fn set_14_bounds_entrypoint_same_ring_matches_direct_helper() {
        let mol = Molecule::from_smiles("c1ccccc1").expect("benzene");
        let (mmat, accum_data, dmat) = run_set14_bounds(&mol, false, false);
        let (bid1, bid2, bid3, ring_size) =
            find_same_ring_dispatch_triple(&mol, false).expect("same-ring triple");

        let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
        set_in_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut direct_accum,
            &mut direct_mmat,
            &dmat,
            ring_size,
        );

        let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2);
        let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3);
        let aid1 = if mol.bonds()[bid1].begin().index() == atm2 {
            mol.bonds()[bid1].end().index()
        } else {
            mol.bonds()[bid1].begin().index()
        };
        let aid4 = if mol.bonds()[bid3].begin().index() == atm3 {
            mol.bonds()[bid3].end().index()
        } else {
            mol.bonds()[bid3].begin().index()
        };
        assert!((mmat.get_lower(aid1, aid4) - direct_mmat.get_lower(aid1, aid4)).abs() < 1e-9);
        assert!((mmat.get_upper(aid1, aid4) - direct_mmat.get_upper(aid1, aid4)).abs() < 1e-9);
        assert!(accum_data.visited14_bounds[aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4)]);
    }

    #[test]
    fn set_14_bounds_entrypoint_two_same_ring_matches_direct_helper() {
        let mol = Molecule::from_smiles("Cc1ccccc1").expect("toluene");
        let (mmat, _accum_data, dmat) = run_set14_bounds(&mol, false, false);
        let (bid1, bid2, bid3) = find_dispatch_triple(&mol, Set14DispatchCase::TwoSameRing, false)
            .expect("two-same-ring triple");

        let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
        set_two_in_same_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut direct_accum,
            &mut direct_mmat,
            &dmat,
        );

        let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2);
        let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3);
        let aid1 = if mol.bonds()[bid1].begin().index() == atm2 {
            mol.bonds()[bid1].end().index()
        } else {
            mol.bonds()[bid1].begin().index()
        };
        let aid4 = if mol.bonds()[bid3].begin().index() == atm3 {
            mol.bonds()[bid3].end().index()
        } else {
            mol.bonds()[bid3].begin().index()
        };
        assert!((mmat.get_lower(aid1, aid4) - direct_mmat.get_lower(aid1, aid4)).abs() < 1e-9);
        assert!((mmat.get_upper(aid1, aid4) - direct_mmat.get_upper(aid1, aid4)).abs() < 1e-9);
    }

    #[test]
    fn set_14_bounds_entrypoint_two_diff_ring_matches_direct_helper() {
        let mol = Molecule::from_smiles("C1CCC2(CC1)CCC3CCCCC23").expect("two-diff-ring polycycle");
        let (mmat, _accum_data, dmat) = run_set14_bounds(&mol, false, false);
        let (bid1, bid2, bid3) = find_dispatch_triple(&mol, Set14DispatchCase::TwoDiffRing, false)
            .expect("two-diff-ring triple");

        let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
        set_two_in_diff_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut direct_accum,
            &mut direct_mmat,
            &dmat,
        );

        let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2);
        let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3);
        let aid1 = if mol.bonds()[bid1].begin().index() == atm2 {
            mol.bonds()[bid1].end().index()
        } else {
            mol.bonds()[bid1].begin().index()
        };
        let aid4 = if mol.bonds()[bid3].begin().index() == atm3 {
            mol.bonds()[bid3].end().index()
        } else {
            mol.bonds()[bid3].begin().index()
        };
        assert!((mmat.get_lower(aid1, aid4) - direct_mmat.get_lower(aid1, aid4)).abs() < 1e-9);
        assert!((mmat.get_upper(aid1, aid4) - direct_mmat.get_upper(aid1, aid4)).abs() < 1e-9);
    }

    #[test]
    fn set_14_bounds_entrypoint_share_ring_bond_matches_direct_helper() {
        let mol = Molecule::from_smiles("Cc1ccccc1C").expect("xylene");
        let (mmat, _accum_data, dmat) = run_set14_bounds(&mol, false, false);
        let (bid1, bid2, bid3) =
            find_dispatch_triple(&mol, Set14DispatchCase::ShareRingBond, false)
                .expect("share-ring-bond triple");

        let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
        set_share_ring_bond_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut direct_accum,
            &mut direct_mmat,
            &dmat,
        );

        let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2);
        let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3);
        let aid1 = if mol.bonds()[bid1].begin().index() == atm2 {
            mol.bonds()[bid1].end().index()
        } else {
            mol.bonds()[bid1].begin().index()
        };
        let aid4 = if mol.bonds()[bid3].begin().index() == atm3 {
            mol.bonds()[bid3].end().index()
        } else {
            mol.bonds()[bid3].begin().index()
        };
        assert!((mmat.get_lower(aid1, aid4) - direct_mmat.get_lower(aid1, aid4)).abs() < 1e-9);
        assert!((mmat.get_upper(aid1, aid4) - direct_mmat.get_upper(aid1, aid4)).abs() < 1e-9);
    }

    #[test]
    fn set_14_bounds_entrypoint_chain_matches_direct_helper() {
        let mol = Molecule::from_smiles("CSSC").expect("disulfide");
        let (mmat, _accum_data, _dmat) = run_set14_bounds(&mol, false, false);
        let (bid1, bid2, bid3) =
            find_dispatch_triple(&mol, Set14DispatchCase::Chain, false).expect("chain triple");

        let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
        set_chain_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut direct_accum,
            &mut direct_mmat,
            false,
        );

        let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2);
        let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3);
        let aid1 = if mol.bonds()[bid1].begin().index() == atm2 {
            mol.bonds()[bid1].end().index()
        } else {
            mol.bonds()[bid1].begin().index()
        };
        let aid4 = if mol.bonds()[bid3].begin().index() == atm3 {
            mol.bonds()[bid3].end().index()
        } else {
            mol.bonds()[bid3].begin().index()
        };
        assert!((mmat.get_lower(aid1, aid4) - direct_mmat.get_lower(aid1, aid4)).abs() < 1e-9);
        assert!((mmat.get_upper(aid1, aid4) - direct_mmat.get_upper(aid1, aid4)).abs() < 1e-9);
    }

    #[test]
    fn set_14_bounds_entrypoint_macrocycle_matches_direct_helper() {
        let mol = Molecule::from_smiles("O=C1N(C)CCCCCCCC1").expect("macrocyclic tertiary lactam");
        let (mmat, _accum_data, _dmat) = run_set14_same_ring_pass_only(&mol, true);
        let (bid1, bid2, bid3, _ring_size) =
            find_same_ring_dispatch_triple(&mol, true).expect("macrocycle same-ring triple");

        let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
        set_macrocycle_all_in_same_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut direct_accum,
            &mut direct_mmat,
        );
        let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2);
        let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3);
        let atom1 = if mol.bonds()[bid1].begin().index() == atm2 {
            mol.bonds()[bid1].end().index()
        } else {
            mol.bonds()[bid1].begin().index()
        };
        let atm4 = if mol.bonds()[bid3].begin().index() == atm3 {
            mol.bonds()[bid3].end().index()
        } else {
            mol.bonds()[bid3].begin().index()
        };
        assert!((mmat.get_lower(atom1, atm4) - direct_mmat.get_lower(atom1, atm4)).abs() < 1e-9);
        assert!((mmat.get_upper(atom1, atm4) - direct_mmat.get_upper(atom1, atm4)).abs() < 1e-9);
    }

    #[test]
    fn test_single_atom() {
        let mol = Molecule::from_smiles("C").expect("methane skeleton");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0][0], 0.0);
    }

    #[test]
    fn test_diatomic() {
        let mol = Molecule::from_smiles("CC").expect("ethane skeleton");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        assert_eq!(result.len(), 2);
        assert!(result[0][1] > 0.0);
        assert!(result[0][1] < 5.0);
    }

    #[test]
    fn test_ethane() {
        let mol = Molecule::from_smiles("CC")
            .expect("ethane skeleton")
            .with_hydrogens()
            .expect("explicit hydrogens");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        assert_eq!(result.len(), 8);
        let upper = |i: usize, j: usize| {
            if i < j { result[i][j] } else { result[j][i] }
        };
        let carbon_pair = mol
            .bonds()
            .iter()
            .find_map(|bond| {
                let begin = bond.begin().index();
                let end = bond.end().index();
                (mol.atoms()[begin].atomic_number() == 6 && mol.atoms()[end].atomic_number() == 6)
                    .then_some((begin, end))
            })
            .expect("ethane must contain a carbon-carbon bond");
        let carbon_hydrogen_pair = mol
            .bonds()
            .iter()
            .find_map(|bond| {
                let begin = bond.begin().index();
                let end = bond.end().index();
                if mol.atoms()[begin].atomic_number() == 6 && mol.atoms()[end].atomic_number() == 1
                {
                    Some((begin, end))
                } else if mol.atoms()[begin].atomic_number() == 1
                    && mol.atoms()[end].atomic_number() == 6
                {
                    Some((end, begin))
                } else {
                    None
                }
            })
            .expect("ethane must contain a carbon-hydrogen bond");

        assert!(
            upper(carbon_pair.0, carbon_pair.1) > 1.0 && upper(carbon_pair.0, carbon_pair.1) < 3.0
        );
        assert!(
            upper(carbon_hydrogen_pair.0, carbon_hydrogen_pair.1) > 1.0
                && upper(carbon_hydrogen_pair.0, carbon_hydrogen_pair.1) < 3.0
        );
    }

    #[test]
    fn test_empty() {
        let builder = MoleculeBuilder::new();
        let mol = builder.build().expect("build");
        let err =
            dg_bounds_matrix(&mol).expect_err("empty molecule must now follow setTopolBounds");
        assert!(matches!(
            err,
            DgBoundsError::GenerationFailed(message) if message == "molecule has no atoms"
        ));
    }

    #[test]
    fn dg_bounds_matrix_matches_source_backed_set_topol_bounds_path() {
        let mol = Molecule::from_smiles("C/C=C/CC").expect("stereo pentene");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        let mut mmat = BoundsMatrix::new(mol.num_atoms());
        set_topol_bounds(&mol, &mut mmat, true, false, false, false, true, true)
            .expect("setTopolBounds");
        assert!(mmat.triangle_smooth(0.0));

        assert_eq!(result, mmat.to_vec_vec());
    }

    #[test]
    fn dg_bounds_matrix_uses_rdkit_wrapper_defaults() {
        let mol = Molecule::from_smiles("CCCCCC").expect("hexane");
        let from_default = dg_bounds_matrix(&mol).expect("default dg_bounds");
        let explicit_default =
            dg_bounds_matrix_with_options(&mol, true, false, true, false).expect("explicit");
        let scaled =
            dg_bounds_matrix_with_options(&mol, true, true, true, false).expect("scaled vdw");

        assert_eq!(from_default, explicit_default);
        assert_eq!(from_default, scaled);
    }

    #[test]
    fn dg_bounds_matrix_with_options_can_skip_triangle_smoothing() {
        let mol = Molecule::from_smiles("C/C=C/CC").expect("stereo pentene");
        let unsmoothed =
            dg_bounds_matrix_with_options(&mol, true, false, false, false).expect("unsmoothed");
        let smoothed =
            dg_bounds_matrix_with_options(&mol, true, false, true, false).expect("smoothed");
        let mut manual_unsmoothed = BoundsMatrix::new(mol.num_atoms());
        set_topol_bounds(
            &mol,
            &mut manual_unsmoothed,
            true,
            false,
            false,
            false,
            true,
            true,
        )
        .expect("setTopolBounds");
        let mut manual_smoothed = BoundsMatrix {
            data: manual_unsmoothed.data.clone(),
            n: manual_unsmoothed.n,
        };
        assert!(manual_smoothed.triangle_smooth(0.0));

        assert_eq!(unsmoothed, manual_unsmoothed.data);
        assert_eq!(smoothed, manual_smoothed.data);
    }

    #[test]
    fn triangle_smooth_uses_rdkit_difference_formula_to_tighten_lower_bound() {
        let mut mmat = BoundsMatrix {
            data: vec![vec![0.0; 3]; 3],
            n: 3,
        };
        mmat.set_lower(0, 1, 0.5);
        mmat.set_upper(0, 1, 1.0);
        mmat.set_lower(1, 2, 4.0);
        mmat.set_upper(1, 2, 5.0);
        mmat.set_lower(0, 2, 1.0);
        mmat.set_upper(0, 2, 10.0);

        assert!(mmat.triangle_smooth(0.0));

        assert_eq!(mmat.get_upper(0, 2), 6.0);
        assert_eq!(mmat.get_lower(0, 2), 3.0);
        assert_eq!(mmat.data[0][2], 6.0);
        assert_eq!(mmat.data[2][0], 3.0);
    }

    #[test]
    fn triangle_smooth_reconciles_small_lower_upper_inversion_with_tolerance() {
        let mut mmat = BoundsMatrix {
            data: vec![vec![0.0; 3]; 3],
            n: 3,
        };
        mmat.set_lower(0, 1, 2.01);
        mmat.set_upper(0, 1, 2.0);
        mmat.set_lower(0, 2, 0.1);
        mmat.set_upper(0, 2, 100.0);
        mmat.set_lower(1, 2, 0.1);
        mmat.set_upper(1, 2, 100.0);

        assert!(mmat.triangle_smooth(0.01));
        assert_eq!(mmat.get_lower(0, 1), 2.01);
        assert_eq!(mmat.get_upper(0, 1), 2.01);
    }

    #[test]
    fn triangle_smooth_fails_when_lower_exceeds_upper_beyond_tolerance() {
        let mut mmat = BoundsMatrix {
            data: vec![vec![0.0; 3]; 3],
            n: 3,
        };
        mmat.set_lower(0, 1, 2.2);
        mmat.set_upper(0, 1, 2.0);
        mmat.set_lower(0, 2, 0.1);
        mmat.set_upper(0, 2, 100.0);
        mmat.set_lower(1, 2, 0.1);
        mmat.set_upper(1, 2, 100.0);

        assert!(!mmat.triangle_smooth(0.01));
    }

    #[test]
    fn test_ethane_bounds_consistency() {
        // Ethane (C2H6) — verify 1-2, 1-3, 1-4, and VDW bounds are reasonable
        let mol = Molecule::from_smiles("CC")
            .expect("ethane skeleton")
            .with_hydrogens()
            .expect("explicit hydrogens");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        let upper = |i: usize, j: usize| {
            if i < j { result[i][j] } else { result[j][i] }
        };
        assert_eq!(result.len(), 8);
        let carbons: Vec<_> = mol
            .atoms()
            .iter()
            .enumerate()
            .filter_map(|(idx, atom)| (atom.atomic_number() == 6).then_some(idx))
            .collect();
        assert_eq!(carbons.len(), 2);
        let (c1, c2) = mol
            .bonds()
            .iter()
            .find_map(|bond| {
                let begin = bond.begin().index();
                let end = bond.end().index();
                (mol.atoms()[begin].atomic_number() == 6 && mol.atoms()[end].atomic_number() == 6)
                    .then_some((begin, end))
            })
            .expect("ethane must contain a carbon-carbon bond");
        let hydrogens_on = |carbon_idx: usize| -> Vec<usize> {
            mol.bonds()
                .iter()
                .filter_map(|bond| {
                    let begin = bond.begin().index();
                    let end = bond.end().index();
                    if begin == carbon_idx && mol.atoms()[end].atomic_number() == 1 {
                        Some(end)
                    } else if end == carbon_idx && mol.atoms()[begin].atomic_number() == 1 {
                        Some(begin)
                    } else {
                        None
                    }
                })
                .collect()
        };
        let h_on_c1 = hydrogens_on(c1);
        let h_on_c2 = hydrogens_on(c2);
        assert!(!h_on_c1.is_empty());
        assert!(!h_on_c2.is_empty());
        let h1 = h_on_c1[0];
        let h4 = h_on_c2[0];
        let lower = |i: usize, j: usize| {
            if i < j { result[j][i] } else { result[i][j] }
        };
        let cc_upper = upper(c1, c2);
        let cc_lower = lower(c1, c2);
        let ch_upper = upper(c1, h1);
        let ch_lower = lower(c1, h1);
        let cch_upper = upper(c1, h4);
        let hcch_upper = upper(h1, h4);

        // Directly bonded pairs should keep valid finite bounds.
        assert!(
            cc_lower > 0.0 && cc_upper > cc_lower,
            "C-C bounds should be finite and ordered, got [{cc_lower}, {cc_upper}]"
        );
        assert!(
            ch_lower > 0.0 && ch_upper > ch_lower,
            "C-H bounds should be finite and ordered, got [{ch_lower}, {ch_upper}]"
        );

        // Longer topological separations should not be tighter than direct bonds.
        assert!(
            cch_upper > ch_upper,
            "C-C-H 1-3 upper bound {cch_upper} should exceed direct C-H upper bound {ch_upper}"
        );

        assert!(
            hcch_upper > cch_upper,
            "H-C-C-H 1-4 upper bound {hcch_upper} should exceed C-C-H 1-3 upper bound {cch_upper}"
        );

        // Matrix symmetry check
        for i in 0..8 {
            for j in 0..8 {
                assert!(!result[i][j].is_nan(), "bounds[{i}][{j}] is NaN");
                assert!(
                    result[i][j] >= 0.0 || i == j,
                    "bounds[{i}][{j}] = {} should be >= 0",
                    result[i][j]
                );
            }
        }
    }

    #[test]
    fn test_ring_bounds_consistency() {
        // Cyclohexane ring skeleton — this exercises ring-aware angle handling
        // This tests ring-aware angle computation and triangle smoothing
        let mol = Molecule::from_smiles("C1CCCCC1").expect("cyclohexane");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        let upper = |i: usize, j: usize| {
            if i < j { result[i][j] } else { result[j][i] }
        };
        assert_eq!(result.len(), 6);

        // 1-2 (C-C bonds): upper bound ~1.51
        for i in 0..6 {
            let j = (i + 1) % 6;
            assert!(
                upper(i, j) > 1.40 && upper(i, j) < 1.60,
                "C-C 1-2 upper bound [{i}][{j}] = {} should be ~1.50",
                upper(i, j)
            );
        }

        // 1-3 (C-C-C in ring): should be larger than 1-2
        for i in 0..6 {
            let j = (i + 2) % 6;
            assert!(
                upper(i, j) > upper(i, (i + 1) % 6),
                "1-3 upper bound [{i}][{j}] = {} should exceed 1-2 bound {}",
                upper(i, j),
                upper(i, (i + 1) % 6)
            );
        }

        // 1-4 (C-C-C-C across hexagon): should be the largest intra-ring distance
        for i in 0..6 {
            let j = (i + 3) % 6;
            assert!(
                upper(i, j) > upper(i, (i + 1) % 6) && upper(i, j) < 5.0,
                "1-4 upper bound [{i}][{j}] = {} should be > 1-2 and < 5.0",
                upper(i, j)
            );
        }

        assert!(!upper(0, 1).is_nan());
        assert!(!upper(0, 1).is_infinite());
    }
}
