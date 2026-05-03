use crate::{Atom, Bond, BondDirection, BondOrder, ChiralTag, Molecule, ValenceAssignment};
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

struct SmilesWriteState<'a> {
    mol: &'a Molecule,
    params: &'a SmilesWriteParams,
    valence: ValenceAssignment,
}

impl<'a> SmilesWriteState<'a> {
    fn new(mol: &'a Molecule, params: &'a SmilesWriteParams) -> Result<Self, SmilesWriteError> {
        let valence = crate::assign_valence(mol, crate::ValenceModel::RdkitLike).map_err(|_| {
            SmilesWriteError::UnsupportedPath(
                "RDKit-like valence assignment required by GetAtomSmiles failed",
            )
        })?;
        Ok(Self {
            mol,
            params,
            valence,
        })
    }

    fn total_num_hs(&self, atom_index: usize) -> usize {
        self.mol.atoms()[atom_index].explicit_hydrogens as usize
            + self.valence.implicit_hydrogens[atom_index] as usize
    }

    fn total_valence(&self, atom_index: usize) -> i32 {
        i32::from(self.valence.explicit_valence[atom_index])
            + i32::from(self.valence.implicit_hydrogens[atom_index])
    }

    fn atom_needs_bracket(&self, atom: &Atom, at_string: &str) -> bool {
        let num = atom.atomic_num;
        if !in_organic_subset(num) {
            return true;
        }
        if atom.formal_charge != 0 {
            return true;
        }
        if atom.atom_map_num.is_some() && !self.params.ignore_atom_map_numbers {
            return true;
        }
        if self.params.do_isomeric_smiles && (atom.isotope.is_some() || !at_string.is_empty()) {
            return true;
        }
        if atom.num_radical_electrons != 0 {
            return true;
        }
        if (num == 7 || num == 15) && atom.is_aromatic && atom.explicit_hydrogens > 0 {
            return true;
        }

        let total_valence = self.total_valence(atom.index);
        let total_num_hs = self.total_num_hs(atom.index);
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
            return true;
        }
        if bond_is_to_metal(self.mol, atom.index) {
            return true;
        }
        false
    }

    fn atom_smiles(&self, atom: &Atom) -> Result<String, SmilesWriteError> {
        let fc = atom.formal_charge;
        let isotope = atom.isotope;
        let mut symb = crate::periodic_table::element_symbol(atom.atomic_num)
            .ok_or(SmilesWriteError::UnsupportedPath(
                "element symbol lookup outside the RDKit periodic table",
            ))?
            .to_string();
        let at_string = if self.params.do_isomeric_smiles {
            atom_chirality_info(atom)
        } else {
            String::new()
        };
        let needs_bracket = if !self.params.all_hs_explicit {
            self.atom_needs_bracket(atom, &at_string)
        } else {
            true
        };

        let mut res = String::new();
        if needs_bracket {
            res.push('[');
        }
        if let Some(isotope) = isotope
            && self.params.do_isomeric_smiles
        {
            res.push_str(&isotope.to_string());
        }
        if !self.params.do_kekule && atom.is_aromatic && symb.as_bytes()[0].is_ascii_uppercase() {
            match atom.atomic_num {
                5 | 6 | 7 | 8 | 14 | 15 | 16 | 33 | 34 | 52 => {
                    symb.replace_range(0..1, &symb[0..1].to_ascii_lowercase());
                }
                _ => {}
            }
        }
        res.push_str(&symb);
        res.push_str(&at_string);

        if needs_bracket {
            let tot_num_hs = self.total_num_hs(atom.index);
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
                && !self.params.ignore_atom_map_numbers
            {
                res.push(':');
                res.push_str(&map_num.to_string());
            }
            res.push(']');
        }
        Ok(res)
    }

    fn bond_smiles(&self, bond: &Bond, atom_to_left_idx: Option<usize>) -> String {
        get_bond_smiles(self.mol, bond, self.params, atom_to_left_idx)
    }
}

pub fn mol_to_smiles(
    mol: &Molecule,
    params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    if mol.atoms().is_empty() {
        return Ok(String::new());
    }
    if let Some(rooted_at_atom) = params.rooted_at_atom
        && rooted_at_atom >= mol.atoms().len()
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
        let rank_state = SmilesWriteState::new(mol, params)?;
        let rank_ring_stereo_atoms = if params.do_isomeric_smiles {
            Some(crate::canon_smiles::find_chiral_atom_special_cases(mol)?)
        } else {
            None
        };
        let atom_ranks = crate::canon_smiles::rank_mol_atoms_with_valence(
            mol,
            &rank_state.valence,
            true,
            params.do_isomeric_smiles,
            params.do_isomeric_smiles,
            true,
            false,
            false,
            rank_ring_stereo_atoms.as_deref(),
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
        let traversal_ring_stereo_atoms_storage;
        let traversal_state_storage;
        let (emit_state, traversal_ring_stereo_atoms) = if let Some(write_mol) = write_mol.as_ref()
        {
            traversal_state_storage = SmilesWriteState::new(write_mol, params)?;
            traversal_ring_stereo_atoms_storage = if params.do_isomeric_smiles {
                Some(crate::canon_smiles::find_chiral_atom_special_cases(
                    write_mol,
                )?)
            } else {
                None
            };
            (
                &traversal_state_storage,
                traversal_ring_stereo_atoms_storage.as_deref(),
            )
        } else {
            (&rank_state, rank_ring_stereo_atoms.as_deref())
        };
        for fragment in fragments {
            let start = fragment
                .iter()
                .copied()
                .min_by_key(|&atom_idx| atom_ranks[atom_idx])
                .ok_or(SmilesWriteError::UnsupportedPath(
                    "canonical fragment start selection failed",
                ))?;
            let traversal = crate::canon_smiles::canonicalize_fragment_with_valence(
                traversal_mol,
                start,
                &atom_ranks,
                params.do_isomeric_smiles,
                params.do_random,
                true,
                &emit_state.valence,
                traversal_ring_stereo_atoms,
            )?;
            pieces.push(emit_fragment_smiles(&emit_state, &traversal)?);
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
    let emit_state = SmilesWriteState::new(traversal_mol, params)?;
    for fragment in fragments {
        let start = fragment[0];
        let atom_ranks: Vec<u32> = (0..mol.atoms().len()).map(|idx| idx as u32).collect();
        let traversal = crate::canon_smiles::canonicalize_fragment_with_valence(
            traversal_mol,
            start,
            &atom_ranks,
            params.do_isomeric_smiles,
            false,
            true,
            &emit_state.valence,
            None,
        )?;
        pieces.push(emit_fragment_smiles(&emit_state, &traversal)?);
    }
    Ok(pieces.join("."))
}

fn connected_components(mol: &Molecule) -> Vec<Vec<usize>> {
    let adjacency = mol
        .adjacency()
        .cloned()
        .unwrap_or_else(|| crate::AdjacencyList::from_topology(mol.atoms().len(), mol.bonds()));
    let mut seen = vec![false; mol.atoms().len()];
    let mut out = Vec::new();
    for atom_idx in 0..mol.atoms().len() {
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
    state: &SmilesWriteState<'_>,
    traversal: &crate::canon_smiles::FragmentTraversal,
) -> Result<String, SmilesWriteError> {
    let mol = state.mol;
    let mut res = String::new();
    let mut ring_closure_map = std::collections::BTreeMap::<usize, usize>::new();
    let mut ring_closures_to_erase = Vec::<usize>::new();
    for elem in &traversal.mol_stack {
        match elem {
            crate::canon_smiles::MolStackElem::Atom { atom_idx } => {
                for key in ring_closures_to_erase.drain(..) {
                    ring_closure_map.remove(&key);
                }
                let mut atom = mol.atoms()[*atom_idx].clone();
                if let Some(tag) = traversal.chiral_tag_overrides[*atom_idx] {
                    atom.chiral_tag = tag;
                }
                if let Some(perm) = traversal.chiral_permutation_overrides[*atom_idx] {
                    atom.props
                        .insert("_chiralPermutation".to_owned(), perm.to_string());
                }
                res.push_str(&state.atom_smiles(&atom)?);
            }
            crate::canon_smiles::MolStackElem::Bond {
                bond_idx,
                atom_to_left_idx,
                ..
            } => {
                if *atom_to_left_idx != usize::MAX {
                    let mut bond = mol.bonds()[*bond_idx].clone();
                    if let Some(direction) = traversal.bond_direction_overrides[*bond_idx] {
                        bond.direction = direction;
                    }
                    res.push_str(&state.bond_smiles(&bond, Some(*atom_to_left_idx)));
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

fn atom_chirality_info(atom: &Atom) -> String {
    match atom.chiral_tag {
        ChiralTag::TetrahedralCw => "@@".to_owned(),
        ChiralTag::TetrahedralCcw => "@".to_owned(),
        ChiralTag::TrigonalBipyramidal => {
            let perm = atom
                .props
                .get("_chiralPermutation")
                .map(String::as_str)
                .unwrap_or("1");
            format!("@TB{perm}")
        }
        ChiralTag::Unspecified => String::new(),
    }
}

fn bond_is_to_metal(mol: &Molecule, atom_index: usize) -> bool {
    mol.bonds().iter().any(|bond| {
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
            mol.atoms()[other].atomic_num,
            3 | 4 | 11 | 12 | 13 | 19 | 20 | 21..=32 | 37..=51 | 55..=84 | 87..=116
        )
    })
}

pub fn get_atom_smiles(
    mol: &Molecule,
    atom: &Atom,
    params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    SmilesWriteState::new(mol, params)?.atom_smiles(atom)
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
    let a1 = &mol.atoms()[atom_to_left_idx];
    let a2 = &mol.atoms()[if bond.begin_atom == atom_to_left_idx {
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
