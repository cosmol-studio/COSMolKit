use crate::valence::assign_valence;
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
}

pub fn add_hydrogens_in_place(molecule: &mut Molecule) -> Result<(), AddHydrogensError> {
    let n = molecule.atoms.len();
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
        if molecule.atoms[i].atomic_num == 1 {
            continue;
        }
        add_counts[i] = assignment.implicit_hydrogens[i] as usize
            + molecule.atoms[i].explicit_hydrogens as usize;
    }

    for (i, cnt) in add_counts.into_iter().enumerate() {
        if cnt == 0 {
            continue;
        }
        molecule.atoms[i].explicit_hydrogens = 0;
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
            });
            molecule.add_bond(Bond {
                index: 0,
                begin_atom: i,
                end_atom: h_idx,
                order: BondOrder::Single,
                direction: BondDirection::None,
                stereo: BondStereo::None,
                stereo_atoms: Vec::new(),
            });
        }
    }

    molecule.rebuild_adjacency();
    Ok(())
}

pub fn remove_hydrogens_in_place(molecule: &mut Molecule) -> Result<(), RemoveHydrogensError> {
    let n = molecule.atoms.len();
    if n == 0 {
        return Ok(());
    }

    let mut degree = vec![0usize; n];
    let mut h_neighbor = vec![None::<usize>; n];
    for bond in &molecule.bonds {
        degree[bond.begin_atom] += 1;
        degree[bond.end_atom] += 1;
        if molecule.atoms[bond.begin_atom].atomic_num == 1 {
            h_neighbor[bond.begin_atom] = Some(bond.end_atom);
        }
        if molecule.atoms[bond.end_atom].atomic_num == 1 {
            h_neighbor[bond.end_atom] = Some(bond.begin_atom);
        }
    }

    let mut remove = vec![false; n];
    for (idx, atom) in molecule.atoms.iter().enumerate() {
        if atom.atomic_num != 1 {
            continue;
        }
        if degree[idx] == 0 || degree[idx] > 1 {
            continue;
        }
        if atom.isotope.is_some() || atom.formal_charge == -1 {
            continue;
        }
        let Some(neighbor_idx) = h_neighbor[idx] else {
            return Err(RemoveHydrogensError::UnsupportedHydrogen { atom_index: idx });
        };
        if molecule.atoms[neighbor_idx].atomic_num == 0 {
            continue;
        }
        remove[idx] = true;
    }

    if remove.iter().all(|x| !*x) {
        return Ok(());
    }

    for idx in (0..n).rev() {
        if !remove[idx] {
            continue;
        }
        let Some(heavy_idx) = h_neighbor[idx] else {
            return Err(RemoveHydrogensError::UnsupportedHydrogen { atom_index: idx });
        };
        let heavy = &mut molecule.atoms[heavy_idx];
        if heavy.no_implicit || !matches!(heavy.chiral_tag, ChiralTag::Unspecified) {
            heavy.explicit_hydrogens = heavy.explicit_hydrogens.saturating_add(1);
        }
    }

    let mut old_to_new = vec![None::<usize>; n];
    let mut atoms = Vec::with_capacity(n - remove.iter().filter(|x| **x).count());
    for (old_idx, atom) in molecule.atoms.iter().enumerate() {
        if remove[old_idx] {
            continue;
        }
        let mut atom = atom.clone();
        atom.index = atoms.len();
        old_to_new[old_idx] = Some(atom.index);
        atoms.push(atom);
    }

    let mut bonds = Vec::new();
    for bond in &molecule.bonds {
        if remove[bond.begin_atom] || remove[bond.end_atom] {
            continue;
        }
        let Some(begin_atom) = old_to_new[bond.begin_atom] else {
            continue;
        };
        let Some(end_atom) = old_to_new[bond.end_atom] else {
            continue;
        };
        let mut bond = bond.clone();
        bond.index = bonds.len();
        bond.begin_atom = begin_atom;
        bond.end_atom = end_atom;
        bond.stereo_atoms = bond
            .stereo_atoms
            .iter()
            .filter_map(|idx| old_to_new.get(*idx).and_then(|x| *x))
            .collect();
        bonds.push(bond);
    }

    molecule.atoms = atoms;
    molecule.bonds = bonds;
    molecule.rebuild_adjacency();
    Ok(())
}
