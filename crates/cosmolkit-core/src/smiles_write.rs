use crate::{
    Atom, Bond, BondDirection, BondOrder, ChiralTag, Molecule, ValenceModel, assign_valence,
};
use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesWriteParams {
    pub do_isomeric_smiles: bool,
    pub do_kekule: bool,
    pub canonical: bool,
    pub clean_stereo: bool,
    pub all_bonds_explicit: bool,
    pub all_hs_explicit: bool,
    pub do_random: bool,
    pub rooted_at_atom: Option<usize>,
    pub include_dative_bonds: bool,
    pub ignore_atom_map_numbers: bool,
}

impl Default for SmilesWriteParams {
    fn default() -> Self {
        Self {
            do_isomeric_smiles: true,
            do_kekule: false,
            canonical: true,
            clean_stereo: true,
            all_bonds_explicit: false,
            all_hs_explicit: false,
            do_random: false,
            rooted_at_atom: None,
            include_dative_bonds: true,
            ignore_atom_map_numbers: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SmilesWriteError {
    #[error("SMILES writer rootedAtAtom {0} is out of range")]
    RootedAtAtomOutOfRange(usize),
    #[error("SMILES writer unsupported path: {0}")]
    UnsupportedPath(&'static str),
}

const ORGANIC_SUBSET_ATOMS: &[u8] = &[5, 6, 7, 8, 9, 15, 16, 17, 35, 53];

pub fn mol_to_smiles(
    mol: &Molecule,
    params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    if mol.atoms.is_empty() {
        return Ok(String::new());
    }
    if let Some(rooted_at_atom) = params.rooted_at_atom
        && rooted_at_atom >= mol.atoms.len()
    {
        return Err(SmilesWriteError::RootedAtAtomOutOfRange(rooted_at_atom));
    }

    if params.do_random {
        return Err(SmilesWriteError::UnsupportedPath(
            "doRandom path from RDKit SmilesWrite::detail::MolToSmiles is not ported yet",
        ));
    }
    if params.rooted_at_atom.is_some() {
        return Err(SmilesWriteError::UnsupportedPath(
            "rootedAtAtom path from RDKit SmilesWrite::detail::MolToSmiles is not ported yet",
        ));
    }
    if params.canonical {
        let atom_ranks = crate::canon_smiles::rank_mol_atoms(
            mol,
            true,
            params.do_isomeric_smiles,
            params.do_isomeric_smiles,
            true,
            false,
            params.do_isomeric_smiles,
            false,
            true,
        )?;
        let fragments = connected_components(mol);
        let mut pieces = Vec::with_capacity(fragments.len());
        let write_mol = if params.do_kekule {
            let mut owned = mol.clone();
            crate::kekulize::kekulize_in_place(&mut owned, true).map_err(|_| {
                SmilesWriteError::UnsupportedPath(
                    "RDKit FragmentSmilesConstruct doKekule path failed or is not fully ported",
                )
            })?;
            Some(owned)
        } else {
            None
        };
        let traversal_mol = write_mol.as_ref().unwrap_or(mol);
        for fragment in fragments {
            let start = fragment
                .iter()
                .copied()
                .min_by_key(|&atom_idx| atom_ranks[atom_idx])
                .ok_or(SmilesWriteError::UnsupportedPath(
                    "canonical fragment start selection failed",
                ))?;
            let traversal = crate::canon_smiles::canonicalize_fragment(
                traversal_mol,
                start,
                &atom_ranks,
                params.do_isomeric_smiles,
                params.do_random,
                true,
            )?;
            pieces.push(emit_fragment_smiles(traversal_mol, &traversal, params)?);
        }
        pieces.sort();
        return Ok(pieces.join("."));
    }
    let fragments = connected_components(mol);
    let mut pieces = Vec::with_capacity(fragments.len());
    let write_mol = if params.do_kekule {
        let mut owned = mol.clone();
        crate::kekulize::kekulize_in_place(&mut owned, true).map_err(|_| {
            SmilesWriteError::UnsupportedPath(
                "RDKit FragmentSmilesConstruct doKekule path failed or is not fully ported",
            )
        })?;
        Some(owned)
    } else {
        None
    };
    let traversal_mol = write_mol.as_ref().unwrap_or(mol);
    for fragment in fragments {
        let start = fragment[0];
        let atom_ranks: Vec<u32> = (0..mol.atoms.len()).map(|idx| idx as u32).collect();
        let traversal = crate::canon_smiles::canonicalize_fragment(
            traversal_mol,
            start,
            &atom_ranks,
            params.do_isomeric_smiles,
            false,
            true,
        )?;
        pieces.push(emit_fragment_smiles(traversal_mol, &traversal, params)?);
    }
    Ok(pieces.join("."))
}

fn connected_components(mol: &Molecule) -> Vec<Vec<usize>> {
    let adjacency = mol
        .adjacency
        .clone()
        .unwrap_or_else(|| crate::AdjacencyList::from_topology(mol.atoms.len(), &mol.bonds));
    let mut seen = vec![false; mol.atoms.len()];
    let mut out = Vec::new();
    for atom_idx in 0..mol.atoms.len() {
        if seen[atom_idx] {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![atom_idx];
        seen[atom_idx] = true;
        while let Some(curr) = stack.pop() {
            component.push(curr);
            for nbr in adjacency.neighbors_of(curr) {
                if !seen[nbr.atom_index] {
                    seen[nbr.atom_index] = true;
                    stack.push(nbr.atom_index);
                }
            }
        }
        component.sort_unstable();
        out.push(component);
    }
    out
}

fn emit_fragment_smiles(
    mol: &Molecule,
    traversal: &crate::canon_smiles::FragmentTraversal,
    params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    let mut res = String::new();
    let mut ring_closure_map = std::collections::BTreeMap::<usize, usize>::new();
    let mut ring_closures_to_erase = Vec::<usize>::new();
    for elem in &traversal.mol_stack {
        match elem {
            crate::canon_smiles::MolStackElem::Atom { atom_idx } => {
                for key in ring_closures_to_erase.drain(..) {
                    ring_closure_map.remove(&key);
                }
                let mut atom = mol.atoms[*atom_idx].clone();
                if let Some(tag) = traversal.chiral_tag_overrides[*atom_idx] {
                    atom.chiral_tag = tag;
                }
                res.push_str(&get_atom_smiles(mol, &atom, params)?);
            }
            crate::canon_smiles::MolStackElem::Bond {
                bond_idx,
                atom_to_left_idx,
                ..
            } => {
                if *atom_to_left_idx != usize::MAX {
                    let mut bond = mol.bonds[*bond_idx].clone();
                    if let Some(direction) = traversal.bond_direction_overrides[*bond_idx] {
                        bond.direction = direction;
                    }
                    res.push_str(&get_bond_smiles(
                        mol,
                        &bond,
                        params,
                        Some(*atom_to_left_idx),
                    ));
                }
            }
            crate::canon_smiles::MolStackElem::Ring { ring_idx } => {
                let closure_val = if let Some(existing) = ring_closure_map.get(ring_idx).copied() {
                    ring_closures_to_erase.push(*ring_idx);
                    existing
                } else {
                    let used: BTreeSet<usize> = ring_closure_map.values().copied().collect();
                    let mut candidate = 1usize;
                    while used.contains(&candidate) {
                        candidate += 1;
                    }
                    ring_closure_map.insert(*ring_idx, candidate);
                    candidate
                };
                if closure_val < 10 {
                    res.push(char::from(b'0' + closure_val as u8));
                } else if closure_val < 100 {
                    res.push('%');
                    res.push_str(&closure_val.to_string());
                } else {
                    res.push_str("%(");
                    res.push_str(&closure_val.to_string());
                    res.push(')');
                }
            }
            crate::canon_smiles::MolStackElem::BranchOpen { .. } => res.push('('),
            crate::canon_smiles::MolStackElem::BranchClose { .. } => res.push(')'),
        }
    }
    Ok(res)
}

pub fn in_organic_subset(atomic_number: u8) -> bool {
    ORGANIC_SUBSET_ATOMS.contains(&atomic_number)
}

fn atom_chirality_info(atom: &Atom) -> &'static str {
    match atom.chiral_tag {
        ChiralTag::TetrahedralCw => "@@",
        ChiralTag::TetrahedralCcw => "@",
        ChiralTag::Unspecified => "",
    }
}

fn element_symbol(atomic_num: u8) -> Result<&'static str, SmilesWriteError> {
    match atomic_num {
        0 => Ok("*"),
        1 => Ok("H"),
        3 => Ok("Li"),
        5 => Ok("B"),
        6 => Ok("C"),
        7 => Ok("N"),
        8 => Ok("O"),
        9 => Ok("F"),
        11 => Ok("Na"),
        14 => Ok("Si"),
        15 => Ok("P"),
        16 => Ok("S"),
        17 => Ok("Cl"),
        19 => Ok("K"),
        29 => Ok("Cu"),
        45 => Ok("Rh"),
        33 => Ok("As"),
        34 => Ok("Se"),
        35 => Ok("Br"),
        52 => Ok("Te"),
        53 => Ok("I"),
        _ => Err(SmilesWriteError::UnsupportedPath(
            "element symbol lookup outside current RDKit-aligned table is not ported yet",
        )),
    }
}

fn total_num_hs_rdkit_like(mol: &Molecule, atom_index: usize) -> Result<usize, SmilesWriteError> {
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).map_err(|_| {
        SmilesWriteError::UnsupportedPath(
            "RDKit-like valence assignment required by GetAtomSmiles failed",
        )
    })?;
    Ok(mol.atoms[atom_index].explicit_hydrogens as usize
        + assignment.implicit_hydrogens[atom_index] as usize)
}

fn total_valence_rdkit_like(mol: &Molecule, atom_index: usize) -> Result<i32, SmilesWriteError> {
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).map_err(|_| {
        SmilesWriteError::UnsupportedPath(
            "RDKit-like valence assignment required by GetAtomSmiles failed",
        )
    })?;
    Ok(i32::from(assignment.explicit_valence[atom_index])
        + i32::from(assignment.implicit_hydrogens[atom_index]))
}

fn bond_is_to_metal(mol: &Molecule, atom_index: usize) -> bool {
    mol.bonds.iter().any(|bond| {
        let other = if bond.begin_atom == atom_index {
            Some(bond.end_atom)
        } else if bond.end_atom == atom_index {
            Some(bond.begin_atom)
        } else {
            None
        };
        let Some(other) = other else {
            return false;
        };
        matches!(
            mol.atoms[other].atomic_num,
            3 | 4 | 11 | 12 | 13 | 19 | 20 | 21..=32 | 37..=51 | 55..=84 | 87..=116
        )
    })
}

fn atom_needs_bracket(
    mol: &Molecule,
    atom: &Atom,
    at_string: &str,
    params: &SmilesWriteParams,
) -> Result<bool, SmilesWriteError> {
    let num = atom.atomic_num;
    if !in_organic_subset(num) {
        return Ok(true);
    }
    if atom.formal_charge != 0 {
        return Ok(true);
    }
    if atom.atom_map_num.is_some() && !params.ignore_atom_map_numbers {
        return Ok(true);
    }
    if params.do_isomeric_smiles && (atom.isotope.is_some() || !at_string.is_empty()) {
        return Ok(true);
    }
    if atom.num_radical_electrons != 0 {
        return Ok(true);
    }
    if (num == 7 || num == 15) && atom.is_aromatic && atom.explicit_hydrogens > 0 {
        return Ok(true);
    }

    let total_valence = total_valence_rdkit_like(mol, atom.index)?;
    let total_num_hs = total_num_hs_rdkit_like(mol, atom.index)?;
    let default_valence = match num {
        5 => Some(3),
        6 => Some(4),
        7 => Some(3),
        8 => Some(2),
        9 => Some(1),
        15 => Some(3),
        16 => Some(2),
        17 => Some(1),
        35 => Some(1),
        53 => Some(1),
        _ => None,
    };
    if let Some(default_valence) = default_valence
        && total_valence != default_valence
        && total_num_hs > 0
    {
        return Ok(true);
    }
    if bond_is_to_metal(mol, atom.index) {
        return Ok(true);
    }
    Ok(false)
}

pub fn get_atom_smiles(
    mol: &Molecule,
    atom: &Atom,
    params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    let fc = atom.formal_charge;
    let isotope = atom.isotope;
    let mut symb = element_symbol(atom.atomic_num)?.to_string();
    let at_string = if params.do_isomeric_smiles {
        atom_chirality_info(atom)
    } else {
        ""
    };
    let needs_bracket = if !params.all_hs_explicit {
        atom_needs_bracket(mol, atom, at_string, params)?
    } else {
        true
    };

    let mut res = String::new();
    if needs_bracket {
        res.push('[');
    }
    if let Some(isotope) = isotope
        && params.do_isomeric_smiles
    {
        res.push_str(&isotope.to_string());
    }
    if !params.do_kekule && atom.is_aromatic && symb.as_bytes()[0].is_ascii_uppercase() {
        match atom.atomic_num {
            5 | 6 | 7 | 8 | 14 | 15 | 16 | 33 | 34 | 52 => {
                symb.replace_range(0..1, &symb[0..1].to_ascii_lowercase());
            }
            _ => {}
        }
    }
    res.push_str(&symb);
    res.push_str(at_string);

    if needs_bracket {
        let tot_num_hs = total_num_hs_rdkit_like(mol, atom.index)?;
        if tot_num_hs > 0 {
            res.push('H');
            if tot_num_hs > 1 {
                res.push_str(&tot_num_hs.to_string());
            }
        }
        if fc > 0 {
            res.push('+');
            if fc > 1 {
                res.push_str(&fc.to_string());
            }
        } else if fc < 0 {
            if fc < -1 {
                res.push_str(&fc.to_string());
            } else {
                res.push('-');
            }
        }
        if let Some(map_num) = atom.atom_map_num
            && !params.ignore_atom_map_numbers
        {
            res.push(':');
            res.push_str(&map_num.to_string());
        }
        res.push(']');
    }
    Ok(res)
}

fn aromatic_bond_smiles_context(
    mol: &Molecule,
    bond: &Bond,
    atom_to_left_idx: usize,
    params: &SmilesWriteParams,
) -> bool {
    if params.do_kekule {
        return false;
    }
    if !matches!(
        bond.order,
        BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
    ) {
        return false;
    }
    let a1 = &mol.atoms[atom_to_left_idx];
    let a2 = &mol.atoms[if bond.begin_atom == atom_to_left_idx {
        bond.end_atom
    } else {
        bond.begin_atom
    }];
    a1.is_aromatic && a2.is_aromatic && (a1.atomic_num != 0 || a2.atomic_num != 0)
}

pub fn get_bond_smiles(
    mol: &Molecule,
    bond: &Bond,
    params: &SmilesWriteParams,
    atom_to_left_idx: Option<usize>,
) -> String {
    let atom_to_left_idx = atom_to_left_idx.unwrap_or(bond.begin_atom);
    let aromatic = aromatic_bond_smiles_context(mol, bond, atom_to_left_idx, params);
    match bond.order {
        BondOrder::Single => match bond.direction {
            BondDirection::EndDownRight => {
                if params.all_bonds_explicit || params.do_isomeric_smiles {
                    "\\".to_string()
                } else {
                    String::new()
                }
            }
            BondDirection::EndUpRight => {
                if params.all_bonds_explicit || params.do_isomeric_smiles {
                    "/".to_string()
                } else {
                    String::new()
                }
            }
            BondDirection::None => {
                if params.all_bonds_explicit || (aromatic && !bond.is_aromatic) {
                    "-".to_string()
                } else {
                    String::new()
                }
            }
        },
        BondOrder::Double => {
            if !aromatic || !bond.is_aromatic || params.all_bonds_explicit {
                "=".to_string()
            } else {
                String::new()
            }
        }
        BondOrder::Triple => "#".to_string(),
        BondOrder::Quadruple => "$".to_string(),
        BondOrder::Aromatic => match bond.direction {
            BondDirection::EndDownRight => {
                if params.all_bonds_explicit || params.do_isomeric_smiles {
                    "\\".to_string()
                } else {
                    String::new()
                }
            }
            BondDirection::EndUpRight => {
                if params.all_bonds_explicit || params.do_isomeric_smiles {
                    "/".to_string()
                } else {
                    String::new()
                }
            }
            BondDirection::None => {
                if params.all_bonds_explicit || !aromatic {
                    ":".to_string()
                } else {
                    String::new()
                }
            }
        },
        BondOrder::Dative => {
            if bond.begin_atom == atom_to_left_idx {
                "->".to_string()
            } else {
                "<-".to_string()
            }
        }
        BondOrder::Null => "~".to_string(),
    }
}
