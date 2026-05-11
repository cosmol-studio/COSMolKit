use crate::valence::{assign_valence, valence_list};
use crate::{Atom, Bond, BondDirection, BondOrder, BondStereo, ChiralTag, Molecule, ValenceModel};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum AddHydrogensError {
    #[error("unsupported atom for AddHs at index {atom_index} with atomic number {atomic_num}")]
    UnsupportedAtom { atom_index: usize, atomic_num: u8 },
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum RemoveHydrogensError {
    #[error("unsupported hydrogen removal at atom index {atom_index}")]
    UnsupportedHydrogen { atom_index: usize },
    #[error("{0}")]
    Sanitize(#[from] crate::sanitize::SanitizeError),
}

#[derive(Debug, Copy, Clone)]
struct RemoveHydrogensParams {
    remove_degree_zero: bool,
    remove_higher_degrees: bool,
    remove_only_h_neighbors: bool,
    remove_isotopes: bool,
    remove_dummy_neighbors: bool,
    remove_defining_bond_stereo: bool,
    remove_mapped: bool,
    update_explicit_count: bool,
    remove_hydrides: bool,
}

impl Default for RemoveHydrogensParams {
    fn default() -> Self {
        Self {
            remove_degree_zero: false,
            remove_higher_degrees: false,
            remove_only_h_neighbors: false,
            remove_isotopes: false,
            remove_dummy_neighbors: false,
            remove_defining_bond_stereo: false,
            remove_mapped: true,
            update_explicit_count: false,
            remove_hydrides: false,
        }
    }
}

pub fn add_hydrogens_in_place(molecule: &mut Molecule) -> Result<(), AddHydrogensError> {
    let n = molecule.atoms().len();
    if n == 0 {
        return Ok(());
    }
    let assignment =
        assign_valence(molecule, ValenceModel::RdkitLike).map_err(|err| match err {
            crate::ValenceError::InvalidValence {
                atom_index,
                atomic_num,
                ..
            } => AddHydrogensError::UnsupportedAtom {
                atom_index,
                atomic_num,
            },
            crate::ValenceError::NotImplemented => AddHydrogensError::UnsupportedAtom {
                atom_index: 0,
                atomic_num: 0,
            },
        })?;

    let mut add_counts = vec![0usize; n];
    for i in 0..n {
        if molecule.atoms()[i].atomic_num == 1 {
            continue;
        }
        add_counts[i] = assignment.implicit_hydrogens[i] as usize
            + molecule.atoms()[i].explicit_hydrogens as usize;
    }

    for (i, cnt) in add_counts.into_iter().enumerate() {
        if cnt == 0 {
            continue;
        }
        molecule.atoms_mut()[i].explicit_hydrogens = 0;
        for _ in 0..cnt {
            let h_idx = molecule.add_atom(Atom {
                index: 0,
                atomic_num: 1,
                is_aromatic: false,
                formal_charge: 0,
                explicit_hydrogens: 0,
                no_implicit: false,
                num_radical_electrons: 0,
                chiral_tag: ChiralTag::Unspecified,
                isotope: None,
                atom_map_num: None,
                props: Default::default(),
                query: None,
                rdkit_cip_rank: None,
            });
            molecule.add_bond(Bond {
                index: 0,
                begin_atom: i,
                end_atom: h_idx,
                order: BondOrder::Single,
                is_aromatic: false,
                direction: BondDirection::None,
                stereo: BondStereo::None,
                stereo_atoms: Vec::new(),
                molfile_query_bond_code: None,
                props: Default::default(),
                query: None,
            });
        }
    }

    molecule.rebuild_adjacency();
    Ok(())
}

pub fn remove_hydrogens_in_place(molecule: &mut Molecule) -> Result<(), RemoveHydrogensError> {
    remove_hydrogens_with_sanitize_in_place(molecule, true)
}

pub fn remove_hydrogens_with_sanitize_in_place(
    molecule: &mut Molecule,
    sanitize: bool,
) -> Result<(), RemoveHydrogensError> {
    let removed =
        remove_hydrogens_with_params_in_place(molecule, RemoveHydrogensParams::default())?;
    if removed && sanitize {
        crate::sanitize::apply_sanitize_pipeline(
            molecule,
            crate::sanitize::SanitizeOps::SUPPORTED_ALL,
        )?;
    }
    Ok(())
}

pub(crate) fn remove_hydrogens_after_smiles_parse_in_place(
    molecule: &mut Molecule,
) -> Result<(), RemoveHydrogensError> {
    remove_hydrogens_with_params_in_place(
        molecule,
        RemoveHydrogensParams {
            update_explicit_count: true,
            ..Default::default()
        },
    )?;
    Ok(())
}

pub(crate) fn adjust_hydrogens_after_aromaticity_in_place(
    molecule: &mut Molecule,
    original_implicit_hydrogens: &[u8],
) {
    let Ok(assignment) = assign_valence(molecule, ValenceModel::RdkitLike) else {
        return;
    };
    for (idx, atom) in molecule.atoms_mut().iter_mut().enumerate() {
        let Some(&original_implicit) = original_implicit_hydrogens.get(idx) else {
            continue;
        };
        let new_implicit = assignment.implicit_hydrogens[idx];
        if new_implicit < original_implicit {
            atom.explicit_hydrogens = atom
                .explicit_hydrogens
                .saturating_add(original_implicit - new_implicit);
        }
    }
}

fn remove_hydrogens_with_params_in_place(
    molecule: &mut Molecule,
    params: RemoveHydrogensParams,
) -> Result<bool, RemoveHydrogensError> {
    let n = molecule.atoms().len();
    if n == 0 {
        return Ok(false);
    }

    let mut degree = vec![0usize; n];
    let mut h_neighbor = vec![None::<usize>; n];
    let mut h_bond = vec![None::<usize>; n];
    for bond in molecule.bonds() {
        degree[bond.begin_atom] += 1;
        degree[bond.end_atom] += 1;
        if molecule.atoms()[bond.begin_atom].atomic_num == 1 {
            h_neighbor[bond.begin_atom] = Some(bond.end_atom);
            h_bond[bond.begin_atom] = Some(bond.index);
        }
        if molecule.atoms()[bond.end_atom].atomic_num == 1 {
            h_neighbor[bond.end_atom] = Some(bond.begin_atom);
            h_bond[bond.end_atom] = Some(bond.index);
        }
    }

    let mut remove = vec![false; n];
    for (idx, atom) in molecule.atoms().iter().enumerate() {
        if should_remove_hydrogen(molecule, idx, atom, &degree, &h_neighbor, &h_bond, params)? {
            remove[idx] = true;
        }
    }

    if remove.iter().all(|x| !*x) {
        return Ok(false);
    }

    let cached_total_valence =
        assign_valence(molecule, ValenceModel::RdkitLike)
            .ok()
            .map(|assignment| {
                assignment
                    .explicit_valence
                    .iter()
                    .zip(assignment.implicit_hydrogens.iter())
                    .map(|(explicit, implicit)| *explicit as i32 + *implicit as i32)
                    .collect::<Vec<_>>()
            });

    for idx in (0..n).rev() {
        if remove[idx] {
            remove_hydrogen_atom_at(
                molecule,
                idx,
                params.update_explicit_count,
                cached_total_valence.as_deref(),
                h_neighbor[idx].zip(h_bond[idx]),
            )?;
        }
    }

    for atom in molecule.atoms_mut() {
        if !atom.no_implicit
            && !matches!(atom.chiral_tag, ChiralTag::Unspecified)
            && atom.explicit_hydrogens > 1
        {
            atom.explicit_hydrogens = 0;
        }
    }

    molecule.rebuild_adjacency();
    Ok(true)
}

fn should_remove_hydrogen(
    molecule: &Molecule,
    atom_index: usize,
    atom: &Atom,
    degree: &[usize],
    h_neighbor: &[Option<usize>],
    h_bond: &[Option<usize>],
    params: RemoveHydrogensParams,
) -> Result<bool, RemoveHydrogensError> {
    if atom.atomic_num != 1 {
        return Ok(false);
    }
    if !params.remove_degree_zero && degree[atom_index] == 0 {
        return Ok(false);
    }
    if !params.remove_higher_degrees && degree[atom_index] > 1 {
        return Ok(false);
    }
    if !params.remove_isotopes && atom.isotope.is_some() {
        return Ok(false);
    }
    if !params.remove_mapped && atom.atom_map_num.is_some() {
        return Ok(false);
    }
    if !params.remove_hydrides && atom.formal_charge == -1 {
        return Ok(false);
    }
    if degree[atom_index] == 0 {
        return Ok(true);
    }

    let Some(neighbor_idx) = h_neighbor[atom_index] else {
        return Err(RemoveHydrogensError::UnsupportedHydrogen { atom_index });
    };
    if !params.remove_dummy_neighbors && molecule.atoms()[neighbor_idx].atomic_num == 0 {
        return Ok(false);
    }

    if !params.remove_only_h_neighbors {
        let only_h_neighbors = neighbors_of(molecule, atom_index)
            .into_iter()
            .all(|idx| molecule.atoms()[idx].atomic_num == 1);
        if only_h_neighbors {
            return Ok(false);
        }
    }

    if !params.remove_defining_bond_stereo && degree[neighbor_idx] == 2 {
        let Some(h_bond_idx) = h_bond[atom_index] else {
            return Err(RemoveHydrogensError::UnsupportedHydrogen { atom_index });
        };
        let h_bond_direction = molecule.bonds()[h_bond_idx].direction;
        for bond in incident_bond_indices(molecule, neighbor_idx) {
            if matches!(molecule.bonds()[bond].order, BondOrder::Double)
                && (matches!(
                    molecule.bonds()[bond].stereo,
                    BondStereo::Cis | BondStereo::Trans
                ) || !matches!(h_bond_direction, BondDirection::None))
            {
                return Ok(false);
            }
        }
    }

    Ok(true)
}

fn remove_hydrogen_atom_at(
    molecule: &mut Molecule,
    atom_index: usize,
    update_explicit_count: bool,
    cached_total_valence: Option<&[i32]>,
    known_neighbor_and_bond: Option<(usize, usize)>,
) -> Result<(), RemoveHydrogensError> {
    let Some((bond_index, heavy_idx)) = known_neighbor_and_bond
        .and_then(|(neighbor_idx, bond_idx)| {
            valid_hydrogen_bond_and_neighbor(molecule, atom_index, neighbor_idx, bond_idx)
        })
        .or_else(|| hydrogen_bond_and_neighbor(molecule, atom_index))
    else {
        return Err(RemoveHydrogensError::UnsupportedHydrogen { atom_index });
    };

    update_neighbor_before_hydrogen_removal(
        molecule,
        atom_index,
        heavy_idx,
        bond_index,
        update_explicit_count,
        cached_total_valence,
    );
    remove_atom_at(molecule, atom_index, Some(bond_index));
    Ok(())
}

fn update_neighbor_before_hydrogen_removal(
    molecule: &mut Molecule,
    atom_index: usize,
    heavy_idx: usize,
    bond_index: usize,
    update_explicit_count: bool,
    cached_total_valence: Option<&[i32]>,
) {
    if should_increment_explicit_h_count(
        molecule,
        heavy_idx,
        update_explicit_count,
        cached_total_valence,
    ) {
        molecule.atoms_mut()[heavy_idx].explicit_hydrogens = molecule.atoms()[heavy_idx]
            .explicit_hydrogens
            .saturating_add(1);
    }

    if !matches!(
        molecule.atoms()[heavy_idx].chiral_tag,
        ChiralTag::Unspecified
    ) {
        let reference = incident_bond_indices(molecule, heavy_idx);
        let mut probe: Vec<usize> = reference
            .iter()
            .copied()
            .filter(|idx| *idx != bond_index)
            .collect();
        probe.push(bond_index);
        if count_swaps_to_interconvert(&probe, &reference).unwrap_or(0) % 2 == 1 {
            invert_chiral_tag(&mut molecule.atoms_mut()[heavy_idx]);
        }
    }

    if incident_bond_indices(molecule, heavy_idx).len() == 2 {
        for idx in incident_bond_indices(molecule, heavy_idx) {
            if idx == bond_index {
                continue;
            }
            if matches!(
                molecule.bonds()[idx].stereo,
                BondStereo::Cis | BondStereo::Trans
            ) {
                molecule.bonds_mut()[idx].stereo = BondStereo::None;
                molecule.bonds_mut()[idx].stereo_atoms.clear();
            }
            break;
        }
    }

    if !matches!(molecule.bonds()[bond_index].direction, BondDirection::None) {
        copy_hydrogen_bond_direction_to_neighbor(molecule, heavy_idx, bond_index);
    }
    adjust_stereo_atoms_if_required(molecule, atom_index, heavy_idx);
}

fn should_increment_explicit_h_count(
    molecule: &Molecule,
    heavy_idx: usize,
    update_explicit_count: bool,
    cached_total_valence: Option<&[i32]>,
) -> bool {
    let heavy = &molecule.atoms()[heavy_idx];
    if update_explicit_count
        || heavy.no_implicit
        || !matches!(heavy.chiral_tag, ChiralTag::Unspecified)
    {
        return true;
    }

    let total_valence = cached_total_valence
        .and_then(|values| values.get(heavy_idx))
        .copied()
        .unwrap_or(0);
    let non_default_valence = valence_list(heavy.atomic_num)
        .map(|values| values.iter().skip(1).any(|value| *value == total_valence))
        .unwrap_or(false);
    ((heavy.atomic_num == 7
        || heavy.atomic_num == 15
        || may_need_extra_h(molecule, heavy_idx, cached_total_valence))
        && heavy.is_aromatic)
        || non_default_valence
}

fn may_need_extra_h(
    molecule: &Molecule,
    atom_index: usize,
    cached_total_valence: Option<&[i32]>,
) -> bool {
    let mut single_bonds = 0usize;
    let mut aromatic_bonds = 0usize;
    for bond_idx in incident_bond_indices(molecule, atom_index) {
        match molecule.bonds()[bond_idx].order {
            BondOrder::Single => single_bonds += 1,
            BondOrder::Aromatic => aromatic_bonds += 1,
            _ => return false,
        }
    }
    let total_valence = cached_total_valence
        .and_then(|values| values.get(atom_index))
        .copied()
        .unwrap_or(0);
    single_bonds == 1 && aromatic_bonds == 2 && total_valence == 3
}

fn copy_hydrogen_bond_direction_to_neighbor(
    molecule: &mut Molecule,
    heavy_idx: usize,
    hydrogen_bond_idx: usize,
) {
    let mut found_direction = false;
    let mut other_single_bond = None::<usize>;
    for idx in incident_bond_indices(molecule, heavy_idx) {
        if idx == hydrogen_bond_idx || !matches!(molecule.bonds()[idx].order, BondOrder::Single) {
            continue;
        }
        if matches!(molecule.bonds()[idx].direction, BondDirection::None) {
            other_single_bond = Some(idx);
        } else {
            found_direction = true;
        }
    }
    if found_direction {
        return;
    }
    let Some(other_bond_idx) = other_single_bond else {
        return;
    };

    let mut direction = molecule.bonds()[hydrogen_bond_idx].direction;
    let flip = molecule.bonds()[other_bond_idx].begin_atom == heavy_idx
        && molecule.bonds()[hydrogen_bond_idx].begin_atom == heavy_idx;
    if flip {
        direction = opposite_bond_direction(direction);
    }
    molecule.bonds_mut()[other_bond_idx].direction = direction;
}

fn adjust_stereo_atoms_if_required(
    molecule: &mut Molecule,
    atom_index: usize,
    heavy_idx: usize,
) -> bool {
    if incident_bond_indices(molecule, heavy_idx).len() == 2 {
        return false;
    }

    let incident = incident_bond_indices(molecule, heavy_idx);
    for bond_idx in incident {
        if !matches!(molecule.bonds()[bond_idx].order, BondOrder::Double)
            || !matches!(
                molecule.bonds()[bond_idx].stereo,
                BondStereo::Cis | BondStereo::Trans
            )
        {
            continue;
        }
        let Some(stereo_pos) = molecule.bonds()[bond_idx]
            .stereo_atoms
            .iter()
            .position(|idx| *idx == atom_index)
        else {
            continue;
        };
        let double_neighbor = if molecule.bonds()[bond_idx].begin_atom == heavy_idx {
            molecule.bonds()[bond_idx].end_atom
        } else {
            molecule.bonds()[bond_idx].begin_atom
        };
        for neighbor in neighbors_of(molecule, heavy_idx) {
            if neighbor == double_neighbor || neighbor == atom_index {
                continue;
            }
            molecule.bonds_mut()[bond_idx].stereo_atoms[stereo_pos] = neighbor;
            molecule.bonds_mut()[bond_idx].stereo = match molecule.bonds()[bond_idx].stereo {
                BondStereo::Cis => BondStereo::Trans,
                BondStereo::Trans => BondStereo::Cis,
                other => other,
            };
            return true;
        }
    }
    false
}

fn remove_atom_at(molecule: &mut Molecule, atom_index: usize, known_bond_index: Option<usize>) {
    if remove_trailing_hydrogen_fast_path(molecule, atom_index, known_bond_index) {
        return;
    }

    molecule.atoms_mut().remove(atom_index);
    for (new_idx, atom) in molecule.atoms_mut().iter_mut().enumerate() {
        atom.index = new_idx;
    }

    let mut bonds = Vec::with_capacity(molecule.bonds().len());
    for bond in molecule.bonds() {
        if bond.begin_atom == atom_index || bond.end_atom == atom_index {
            continue;
        }
        let mut bond = bond.clone();
        bond.index = bonds.len();
        if bond.begin_atom > atom_index {
            bond.begin_atom -= 1;
        }
        if bond.end_atom > atom_index {
            bond.end_atom -= 1;
        }
        bond.stereo_atoms = bond
            .stereo_atoms
            .iter()
            .filter_map(|idx| {
                if *idx == atom_index {
                    None
                } else if *idx > atom_index {
                    Some(*idx - 1)
                } else {
                    Some(*idx)
                }
            })
            .collect();
        bonds.push(bond);
    }
    *molecule.bonds_mut() = bonds;
    remove_conformer_atom_at(molecule, atom_index);
    molecule.clear_adjacency_cache();
}

fn remove_trailing_hydrogen_fast_path(
    molecule: &mut Molecule,
    atom_index: usize,
    known_bond_index: Option<usize>,
) -> bool {
    if atom_index + 1 != molecule.atoms().len() {
        return false;
    }
    let Some(bond_index) = known_bond_index else {
        return false;
    };
    if bond_index + 1 != molecule.bonds().len() {
        return false;
    }
    let Some(bond) = molecule.bonds().last() else {
        return false;
    };
    if bond.begin_atom != atom_index && bond.end_atom != atom_index {
        return false;
    }

    molecule.atoms_mut().pop();
    molecule.bonds_mut().pop();
    remove_conformer_atom_at(molecule, atom_index);
    molecule.clear_adjacency_cache();
    true
}

fn remove_conformer_atom_at(molecule: &mut Molecule, atom_index: usize) {
    if molecule.coords_2d().is_some()
        && let Some(coords) = molecule.coords_2d_mut().as_mut()
    {
        if atom_index < coords.len() {
            coords.remove(atom_index);
        } else {
            molecule.set_coords_2d(None);
        }
    }

    let mut keep_source_dim = true;
    if !molecule.conformers_3d().is_empty() {
        let conformers = molecule.conformers_3d_mut();
        for coords in conformers.iter_mut() {
            if atom_index < coords.len() {
                coords.remove(atom_index);
            } else {
                keep_source_dim = false;
            }
        }
        if !keep_source_dim {
            conformers.clear();
        }
    }
    if !keep_source_dim && molecule.coords_2d().is_none() {
        molecule.set_source_coordinate_dim(None);
    }
}

fn valid_hydrogen_bond_and_neighbor(
    molecule: &Molecule,
    atom_index: usize,
    neighbor_idx: usize,
    bond_index: usize,
) -> Option<(usize, usize)> {
    let bond = molecule.bonds().get(bond_index)?;
    if bond.begin_atom == atom_index && bond.end_atom == neighbor_idx {
        Some((bond_index, neighbor_idx))
    } else if bond.end_atom == atom_index && bond.begin_atom == neighbor_idx {
        Some((bond_index, neighbor_idx))
    } else {
        None
    }
}

fn hydrogen_bond_and_neighbor(molecule: &Molecule, atom_index: usize) -> Option<(usize, usize)> {
    let mut out = None;
    for bond in molecule.bonds() {
        if bond.begin_atom == atom_index {
            if out.is_some() {
                return None;
            }
            out = Some((bond.index, bond.end_atom));
        } else if bond.end_atom == atom_index {
            if out.is_some() {
                return None;
            }
            out = Some((bond.index, bond.begin_atom));
        }
    }
    out
}

fn incident_bond_indices(molecule: &Molecule, atom_index: usize) -> Vec<usize> {
    molecule
        .bonds()
        .iter()
        .filter_map(|bond| {
            if bond.begin_atom == atom_index || bond.end_atom == atom_index {
                Some(bond.index)
            } else {
                None
            }
        })
        .collect()
}

fn neighbors_of(molecule: &Molecule, atom_index: usize) -> Vec<usize> {
    molecule
        .bonds()
        .iter()
        .filter_map(|bond| {
            if bond.begin_atom == atom_index {
                Some(bond.end_atom)
            } else if bond.end_atom == atom_index {
                Some(bond.begin_atom)
            } else {
                None
            }
        })
        .collect()
}

fn opposite_bond_direction(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        BondDirection::None => BondDirection::None,
        BondDirection::Unknown => BondDirection::Unknown,
    }
}

fn invert_chiral_tag(atom: &mut Atom) {
    atom.chiral_tag = match atom.chiral_tag {
        ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
        ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
        ChiralTag::TrigonalBipyramidal => ChiralTag::TrigonalBipyramidal,
        ChiralTag::Unspecified => ChiralTag::Unspecified,
    };
}

fn count_swaps_to_interconvert(probe: &[usize], reference: &[usize]) -> Option<usize> {
    if probe.len() != reference.len() {
        return None;
    }
    let mut work = probe.to_vec();
    let mut swaps = 0usize;
    for i in 0..work.len() {
        if work[i] == reference[i] {
            continue;
        }
        let mut found = None;
        for (j, value) in work.iter().enumerate().skip(i + 1) {
            if *value == reference[i] {
                found = Some(j);
                break;
            }
        }
        let j = found?;
        work.swap(i, j);
        swaps += 1;
    }
    Some(swaps)
}

#[cfg(test)]
mod tests {
    use crate::{
        Atom, Bond, BondDirection, BondOrder, BondStereo, ChiralTag, CoordinateDimension, Molecule,
    };
    use glam::{DVec2, DVec3};

    fn test_atom(atomic_num: u8) -> Atom {
        Atom {
            index: 0,
            atomic_num,
            is_aromatic: false,
            formal_charge: 0,
            explicit_hydrogens: 0,
            no_implicit: false,
            num_radical_electrons: 0,
            chiral_tag: ChiralTag::Unspecified,
            isotope: None,
            atom_map_num: None,
            props: Default::default(),
            query: None,
            rdkit_cip_rank: None,
        }
    }

    fn test_bond(begin_atom: usize, end_atom: usize) -> Bond {
        Bond {
            index: 0,
            begin_atom,
            end_atom,
            order: BondOrder::Single,
            is_aromatic: false,
            direction: BondDirection::None,
            stereo: BondStereo::None,
            stereo_atoms: Vec::new(),
            molfile_query_bond_code: None,
            props: Default::default(),
            query: None,
        }
    }

    fn interleaved_hydrogen_molecule() -> Molecule {
        let mut mol = Molecule::new();
        for atomic_num in [6, 1, 8, 1] {
            mol.add_atom(test_atom(atomic_num));
        }
        mol.add_bond(test_bond(0, 1));
        mol.add_bond(test_bond(0, 2));
        mol.add_bond(test_bond(2, 3));
        mol
    }

    #[test]
    fn remove_hydrogens_sanitize_flag_controls_rdkit_cleanup() {
        let raw = Molecule::from_smiles_with_sanitize("[H]CN(=O)=O", false)
            .expect("unsanitized SMILES should parse");

        let removed_only = raw
            .without_hydrogens_with_sanitize(false)
            .expect("hydrogen removal should succeed without sanitizing");
        assert_eq!(removed_only.atoms().len(), 4);
        assert_eq!(
            removed_only
                .atoms()
                .iter()
                .map(|atom| atom.formal_charge)
                .collect::<Vec<_>>(),
            vec![0, 0, 0, 0]
        );

        let sanitized = raw
            .without_hydrogens_with_sanitize(true)
            .expect("RDKit-style RemoveHs default should sanitize after removal");
        assert_eq!(sanitized.atoms().len(), 4);
        assert_eq!(
            sanitized
                .atoms()
                .iter()
                .map(|atom| atom.formal_charge)
                .collect::<Vec<_>>(),
            vec![0, 1, -1, 0]
        );
        assert_eq!(raw.atoms().len(), 5);
    }

    #[test]
    fn remove_hydrogens_filters_3d_conformer_rows_for_interleaved_hydrogens() {
        let mut mol = interleaved_hydrogen_molecule();
        mol.conformers_3d_mut().push(vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
            DVec3::new(3.0, 0.0, 0.0),
        ]);
        mol.conformers_3d_mut().push(vec![
            DVec3::new(10.0, 0.0, 0.0),
            DVec3::new(11.0, 0.0, 0.0),
            DVec3::new(12.0, 0.0, 0.0),
            DVec3::new(13.0, 0.0, 0.0),
        ]);
        mol.set_source_coordinate_dim(Some(CoordinateDimension::ThreeD));

        let removed = mol
            .without_hydrogens_with_sanitize(false)
            .expect("hydrogens should remove");

        assert_eq!(removed.atoms().len(), 2);
        assert_eq!(
            removed.coords_3d().expect("3D coords should remain"),
            &[DVec3::new(0.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0)]
        );
        assert_eq!(
            removed
                .conformer_3d(1)
                .expect("second conformer should remain"),
            &[DVec3::new(10.0, 0.0, 0.0), DVec3::new(12.0, 0.0, 0.0)]
        );
        assert_eq!(
            removed.source_coordinate_dim(),
            Some(CoordinateDimension::ThreeD)
        );
    }

    #[test]
    fn remove_hydrogens_filters_2d_rows_for_interleaved_hydrogens() {
        let mut mol = interleaved_hydrogen_molecule();
        mol.set_coords_2d(Some(vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 0.0),
            DVec2::new(2.0, 0.0),
            DVec2::new(3.0, 0.0),
        ]));
        mol.set_source_coordinate_dim(Some(CoordinateDimension::TwoD));

        let removed = mol
            .without_hydrogens_with_sanitize(false)
            .expect("hydrogens should remove");

        assert_eq!(removed.atoms().len(), 2);
        assert_eq!(
            removed.coords_2d().expect("2D coords should remain"),
            &[DVec2::new(0.0, 0.0), DVec2::new(2.0, 0.0)]
        );
        assert_eq!(
            removed.source_coordinate_dim(),
            Some(CoordinateDimension::TwoD)
        );
    }

    #[test]
    fn remove_hydrogens_filters_both_2d_and_3d_coordinate_stores() {
        let mut mol = interleaved_hydrogen_molecule();
        mol.set_coords_2d(Some(vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 0.0),
            DVec2::new(2.0, 0.0),
            DVec2::new(3.0, 0.0),
        ]));
        mol.conformers_3d_mut().push(vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
            DVec3::new(3.0, 0.0, 0.0),
        ]);
        mol.set_source_coordinate_dim(Some(CoordinateDimension::ThreeD));

        let removed = mol
            .without_hydrogens_with_sanitize(false)
            .expect("hydrogens should remove");

        assert_eq!(removed.atoms().len(), 2);
        assert_eq!(
            removed.coords_2d().expect("2D coords should remain").len(),
            2
        );
        assert_eq!(
            removed.coords_3d().expect("3D coords should remain").len(),
            2
        );
        assert_eq!(
            removed.coords_2d().expect("2D coords should remain"),
            &[DVec2::new(0.0, 0.0), DVec2::new(2.0, 0.0)]
        );
        assert_eq!(
            removed.coords_3d().expect("3D coords should remain"),
            &[DVec3::new(0.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0)]
        );
        assert_eq!(
            removed.source_coordinate_dim(),
            Some(CoordinateDimension::ThreeD)
        );
    }
}
