use crate::valence::assign_valence;
use crate::{Atom, Bond, BondOrder, Molecule, ValenceModel};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AddHydrogensError {
    UnsupportedAtom { atom_index: usize, atomic_num: u8 },
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
                isotope: None,
            });
            molecule.add_bond(Bond {
                index: 0,
                begin_atom: i,
                end_atom: h_idx,
                order: BondOrder::Single,
            });
        }
    }

    molecule.rebuild_adjacency();
    Ok(())
}
