use crate::Molecule;

use super::{AtomMapping, BondMapping, IndexMapping, TopologyOpKind, TopologyOpSpec};

#[derive(Debug)]
pub(crate) struct TopologyEditResult {
    pub(crate) molecule: Molecule,
    pub(crate) atom_mapping: Option<AtomMapping>,
    pub(crate) bond_mapping: Option<BondMapping>,
}

#[inline]
pub(crate) fn apply_topology_op<E, F>(
    mol: &Molecule,
    spec: &'static TopologyOpSpec,
    f: F,
) -> Result<Molecule, E>
where
    F: FnOnce(&Molecule) -> Result<Molecule, E>,
{
    debug_assert_eq!(spec.kind, TopologyOpKind::Weak);
    let out = f(mol)?;
    debug_check_weak_topology_op(mol, &out, spec);
    Ok(out)
}

#[inline]
pub(crate) fn apply_strong_topology_op<E, F>(
    mol: &Molecule,
    spec: &'static TopologyOpSpec,
    f: F,
) -> Result<TopologyEditResult, E>
where
    F: FnOnce(&Molecule) -> Result<TopologyEditResult, E>,
{
    debug_assert_eq!(spec.kind, TopologyOpKind::Strong);
    let edit = f(mol)?;
    debug_check_strong_topology_op(mol, &edit, spec);
    Ok(edit)
}

pub(crate) fn identity_mapping(len: usize) -> IndexMapping {
    let old_to_new = (0..len).map(Some).collect();
    let new_to_old = (0..len).map(Some).collect();
    IndexMapping::new(old_to_new, new_to_old)
}

pub(crate) fn identity_bond_mapping(mol: &Molecule) -> BondMapping {
    identity_mapping(mol.bonds().len())
}

pub(crate) fn identity_then_created_mapping(old_len: usize, new_len: usize) -> IndexMapping {
    debug_assert!(new_len >= old_len);
    let old_to_new = (0..old_len).map(Some).collect();
    let mut new_to_old = Vec::with_capacity(new_len);
    for new_idx in 0..new_len {
        if new_idx < old_len {
            new_to_old.push(Some(new_idx));
        } else {
            new_to_old.push(None);
        }
    }
    IndexMapping::new(old_to_new, new_to_old)
}

#[cfg(debug_assertions)]
fn debug_check_weak_topology_op(before: &Molecule, after: &Molecule, spec: &TopologyOpSpec) {
    debug_assert_eq!(
        before.atoms().len(),
        after.atoms().len(),
        "weak topology operation {} changed atom count",
        spec.name
    );
    debug_assert_eq!(
        before.bonds().len(),
        after.bonds().len(),
        "weak topology operation {} changed bond count",
        spec.name
    );
    debug_assert_conformer_rows_match_atoms(after, spec.name);
    debug_assert_bond_endpoints_valid(after, spec.name);
}

#[cfg(not(debug_assertions))]
fn debug_check_weak_topology_op(_before: &Molecule, _after: &Molecule, _spec: &TopologyOpSpec) {}

#[cfg(debug_assertions)]
fn debug_check_strong_topology_op(
    _before: &Molecule,
    edit: &TopologyEditResult,
    spec: &TopologyOpSpec,
) {
    debug_assert_conformer_rows_match_atoms(&edit.molecule, spec.name);
    debug_assert_bond_endpoints_valid(&edit.molecule, spec.name);
}

#[cfg(not(debug_assertions))]
fn debug_check_strong_topology_op(
    _before: &Molecule,
    _edit: &TopologyEditResult,
    _spec: &TopologyOpSpec,
) {
}

#[cfg(debug_assertions)]
fn debug_assert_conformer_rows_match_atoms(mol: &Molecule, op_name: &str) {
    let atom_count = mol.atoms().len();
    if let Some(coords_2d) = mol.coords_2d() {
        debug_assert_eq!(
            coords_2d.len(),
            atom_count,
            "topology operation {op_name} left 2D coordinate rows misaligned"
        );
    }
    for (idx, coords) in mol.conformers_3d().iter().enumerate() {
        debug_assert_eq!(
            coords.len(),
            atom_count,
            "topology operation {op_name} left 3D conformer {idx} rows misaligned"
        );
    }
}

#[cfg(debug_assertions)]
fn debug_assert_bond_endpoints_valid(mol: &Molecule, op_name: &str) {
    let atom_count = mol.atoms().len();
    for bond in mol.bonds() {
        debug_assert!(
            bond.begin_atom < atom_count,
            "topology operation {op_name} left invalid bond begin endpoint"
        );
        debug_assert!(
            bond.end_atom < atom_count,
            "topology operation {op_name} left invalid bond end endpoint"
        );
        debug_assert_ne!(
            bond.begin_atom, bond.end_atom,
            "topology operation {op_name} left a self-loop bond"
        );
    }
}
