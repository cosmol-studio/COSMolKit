use crate::Molecule;

/// Valence model selector for future compatibility modes.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum ValenceModel {
    /// Intended RDKit-like valence behavior.
    RdkitLike,
}

/// Per-atom valence assignment result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValenceAssignment {
    pub explicit_valence: Vec<u8>,
    pub implicit_hydrogens: Vec<u8>,
}

/// Valence handling errors.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum ValenceError {
    #[error(
        "invalid valence at atom {atom_index}: atomic_num={atomic_num}, formal_charge={formal_charge}"
    )]
    InvalidValence {
        atom_index: usize,
        atomic_num: u8,
        formal_charge: i8,
    },
    #[error("valence assignment path is not implemented")]
    NotImplemented,
}

pub(crate) fn valence_list(atomic_num: u8) -> Option<&'static [i32]> {
    crate::periodic_table::valence_list(atomic_num)
}

/// Returns the RDKit 2026.03.1 valence list for an atomic number.
pub fn rdkit_valence_list(atomic_num: u8) -> Option<&'static [i32]> {
    valence_list(atomic_num)
}

struct ValenceContext {
    bond_valence_sums: Vec<f64>,
    aromatic_atoms: Vec<bool>,
    degrees: Vec<usize>,
}

impl ValenceContext {
    fn from_molecule(molecule: &Molecule) -> Self {
        let mut bond_valence_sums = vec![0.0; molecule.atoms().len()];
        let mut aromatic_atoms = molecule
            .atoms()
            .iter()
            .map(|atom| atom.is_aromatic)
            .collect::<Vec<_>>();
        let mut degrees = vec![0usize; molecule.atoms().len()];

        for bond in molecule.bonds() {
            bond_valence_sums[bond.begin_atom] +=
                bond_valence_contrib_for_atom(bond, bond.begin_atom);
            bond_valence_sums[bond.end_atom] += bond_valence_contrib_for_atom(bond, bond.end_atom);
            degrees[bond.begin_atom] += 1;
            degrees[bond.end_atom] += 1;
            if bond.is_aromatic {
                aromatic_atoms[bond.begin_atom] = true;
                aromatic_atoms[bond.end_atom] = true;
            }
        }

        Self {
            bond_valence_sums,
            aromatic_atoms,
            degrees,
        }
    }
}

fn is_aromatic_atom_with_context(context: &ValenceContext, atom_index: usize) -> bool {
    context.aromatic_atoms[atom_index]
}

fn get_effective_atomic_num(atomic_num: u8, formal_charge: i8) -> Option<u8> {
    let ea = atomic_num as i32 - formal_charge as i32;
    if !(0..=118).contains(&ea) {
        return None;
    }
    Some(ea as u8)
}

fn can_be_hypervalent(atomic_num: u8, effective_atomic_num: u8) -> bool {
    ((effective_atomic_num > 16) && (atomic_num == 15 || atomic_num == 16))
        || ((effective_atomic_num > 34) && (atomic_num == 33 || atomic_num == 34))
}

fn bond_type_as_double(order: crate::BondOrder) -> f64 {
    match order {
        crate::BondOrder::Null => 0.0,
        crate::BondOrder::Single => 1.0,
        crate::BondOrder::Double => 2.0,
        crate::BondOrder::Triple => 3.0,
        crate::BondOrder::Quadruple => 4.0,
        crate::BondOrder::Aromatic => 1.5,
        crate::BondOrder::Dative => 1.0,
        crate::BondOrder::Hydrogen => 0.0,
    }
}

fn bond_valence_contrib_for_atom(bond: &crate::Bond, atom_index: usize) -> f64 {
    if bond.begin_atom != atom_index && bond.end_atom != atom_index {
        return 0.0;
    }
    if matches!(bond.order, crate::BondOrder::Dative) {
        if bond.end_atom == atom_index {
            return 1.0;
        }
        return 0.0;
    }
    bond_type_as_double(bond.order)
}

pub(crate) fn calculate_explicit_valence(
    molecule: &Molecule,
    atom_index: usize,
    strict: bool,
) -> Result<i32, ValenceError> {
    let context = ValenceContext::from_molecule(molecule);
    calculate_explicit_valence_with_context(molecule, &context, atom_index, strict)
}

fn calculate_explicit_valence_with_context(
    molecule: &Molecule,
    context: &ValenceContext,
    atom_index: usize,
    strict: bool,
) -> Result<i32, ValenceError> {
    let atom = &molecule.atoms()[atom_index];
    let ovalens = valence_list(atom.atomic_num).ok_or(ValenceError::NotImplemented)?;
    let mut effective_atomic_num = atom.atomic_num;
    if ovalens.len() > 1 || ovalens[0] != -1 {
        effective_atomic_num = get_effective_atomic_num(atom.atomic_num, atom.formal_charge)
            .ok_or(ValenceError::InvalidValence {
                atom_index,
                atomic_num: atom.atomic_num,
                formal_charge: atom.formal_charge,
            })?;
    }
    let valens = valence_list(effective_atomic_num).ok_or(ValenceError::NotImplemented)?;
    let dv = valens[0];

    let mut accum = atom.explicit_hydrogens as f64 + context.bond_valence_sums[atom_index];

    if accum > dv as f64 && is_aromatic_atom_with_context(context, atom_index) {
        let mut pval = dv;
        for &v in valens {
            if v == -1 {
                break;
            }
            if (v as f64) > accum {
                break;
            }
            pval = v;
        }
        if accum - pval as f64 <= 1.5 {
            accum = pval as f64;
        }
    }
    accum += 0.1;
    let res = accum.round() as i32;

    if strict {
        let mut max_valence = *valens.last().unwrap_or(&-1);
        let mut offset = 0;
        if can_be_hypervalent(atom.atomic_num, effective_atomic_num) {
            max_valence = *ovalens.last().unwrap_or(&-1);
            offset -= atom.formal_charge as i32;
        }
        if atom.atomic_num == 1 && atom.formal_charge == -1 {
            max_valence = 2;
        }
        if max_valence >= 0 && *ovalens.last().unwrap_or(&-1) >= 0 && (res + offset) > max_valence {
            return Err(ValenceError::InvalidValence {
                atom_index,
                atomic_num: atom.atomic_num,
                formal_charge: atom.formal_charge,
            });
        }
    }
    Ok(res)
}

fn calculate_implicit_valence_with_context(
    molecule: &Molecule,
    context: &ValenceContext,
    atom_index: usize,
    explicit_valence: i32,
    strict: bool,
) -> Result<i32, ValenceError> {
    let atom = &molecule.atoms()[atom_index];
    if atom.atomic_num == 0 {
        return Ok(0);
    }
    if atom.no_implicit {
        return Ok(0);
    }
    if atom.atomic_num == 1 && explicit_valence == 0 {
        return match atom.formal_charge {
            -1 | 1 => Ok(0),
            0 => Ok(1),
            _ => {
                if strict {
                    Err(ValenceError::InvalidValence {
                        atom_index,
                        atomic_num: atom.atomic_num,
                        formal_charge: atom.formal_charge,
                    })
                } else {
                    Ok(0)
                }
            }
        };
    }

    let ovalens = valence_list(atom.atomic_num).ok_or(ValenceError::NotImplemented)?;
    let mut effective_atomic_num = atom.atomic_num;
    if ovalens.len() > 1 || ovalens[0] != -1 {
        effective_atomic_num = get_effective_atomic_num(atom.atomic_num, atom.formal_charge)
            .ok_or(ValenceError::InvalidValence {
                atom_index,
                atomic_num: atom.atomic_num,
                formal_charge: atom.formal_charge,
            })?;
    }
    if effective_atomic_num == 0 {
        return Ok(0);
    }

    let mut explicit_plus_rad_v = explicit_valence;
    let mut valens = valence_list(effective_atomic_num).ok_or(ValenceError::NotImplemented)?;
    let dv = valens[0];
    if dv == -1 {
        return Ok(0);
    }

    if can_be_hypervalent(atom.atomic_num, effective_atomic_num) {
        effective_atomic_num = atom.atomic_num;
        explicit_plus_rad_v -= atom.formal_charge as i32;
        valens = valence_list(effective_atomic_num).ok_or(ValenceError::NotImplemented)?;
    }

    if is_aromatic_atom_with_context(context, atom_index) {
        if explicit_plus_rad_v <= dv {
            return Ok(dv - explicit_plus_rad_v);
        }
        let mut satis = false;
        for &v in valens {
            if v <= 0 {
                break;
            }
            if explicit_plus_rad_v == v {
                satis = true;
                break;
            }
        }
        if !satis && strict {
            return Err(ValenceError::InvalidValence {
                atom_index,
                atomic_num: atom.atomic_num,
                formal_charge: atom.formal_charge,
            });
        }
        return Ok(0);
    }

    let mut res = -1;
    for &v in valens {
        if v < 0 {
            break;
        }
        if explicit_plus_rad_v <= v {
            res = v - explicit_plus_rad_v;
            break;
        }
    }
    if res < 0 {
        if strict && *valens.last().unwrap_or(&-1) != -1 && *ovalens.last().unwrap_or(&-1) > 0 {
            return Err(ValenceError::InvalidValence {
                atom_index,
                atomic_num: atom.atomic_num,
                formal_charge: atom.formal_charge,
            });
        }
        return Ok(0);
    }
    Ok(res)
}

fn any_unsupported_features(molecule: &Molecule) -> bool {
    // The current Bond type does not encode DATIVE direction tag variants
    // beyond begin/end orientation, and Atom does not yet carry noImplicit /
    // radical input state. Keep this hard check explicit for future extension.
    let _ = molecule;
    false
}

fn to_u8_checked(v: i32, atom_index: usize, atom: &crate::Atom) -> Result<u8, ValenceError> {
    if !(0..=u8::MAX as i32).contains(&v) {
        return Err(ValenceError::InvalidValence {
            atom_index,
            atomic_num: atom.atomic_num,
            formal_charge: atom.formal_charge,
        });
    }
    Ok(v as u8)
}

fn n_outer_electrons(atomic_num: u8) -> Option<i32> {
    crate::periodic_table::n_outer_electrons(atomic_num)
}

/// RDKit 2026.03.1 MolOps::assignRadicals parity helper.
pub fn assign_radicals_rdkit_2025(
    molecule: &Molecule,
    existing_explicit_valence: &[u8],
) -> Result<Vec<u8>, ValenceError> {
    if existing_explicit_valence.len() != molecule.atoms().len() {
        return Err(ValenceError::NotImplemented);
    }
    let mut radicals: Vec<u8> = molecule
        .atoms()
        .iter()
        .map(|a| a.num_radical_electrons)
        .collect();
    let context = ValenceContext::from_molecule(molecule);

    for (i, atom) in molecule.atoms().iter().enumerate() {
        if !atom.no_implicit || atom.atomic_num == 0 {
            continue;
        }
        let valens = valence_list(atom.atomic_num).ok_or(ValenceError::NotImplemented)?;
        let chg = atom.formal_charge as i32;
        let n_outer = n_outer_electrons(atom.atomic_num).ok_or(ValenceError::NotImplemented)?;
        let value = if valens.len() != 1 || valens[0] != -1 {
            let total_valence = if is_aromatic_atom_with_context(&context, i) {
                // RDKit runs assignRadicals after Kekulize(). In this codebase we
                // do not yet store a separate kekulized bond table, so we use the
                // already-RDKit-aligned explicit valence cache as the sanitized
                // total valence surrogate for aromatic atoms.
                existing_explicit_valence[i] as i32
            } else {
                (atom.explicit_hydrogens as f64 + context.bond_valence_sums[i] + 0.1) as i32
            };
            let base_count = if atom.atomic_num == 1 || atom.atomic_num == 2 {
                2
            } else {
                8
            };
            let mut num_radicals = base_count - n_outer - total_valence + chg;
            if num_radicals < 0 {
                num_radicals = 0;
                if valens.len() > 1 {
                    for &v in valens {
                        if v - total_valence + chg >= 0 {
                            num_radicals = v - total_valence + chg;
                            break;
                        }
                    }
                }
            }
            let num_radicals2 = n_outer - total_valence - chg;
            if num_radicals2 >= 0 {
                num_radicals = num_radicals.min(num_radicals2);
            }
            num_radicals
        } else {
            let degree = context.degrees[i];
            if degree > 0 {
                0
            } else {
                let mut n_valence = n_outer - chg;
                if n_valence < 0 {
                    n_valence = 0;
                }
                n_valence % 2
            }
        };
        radicals[i] = to_u8_checked(value, i, atom)?;
    }
    Ok(radicals)
}

/// Compute explicit valence and implicit hydrogen assignment.
pub fn assign_valence(
    molecule: &Molecule,
    model: ValenceModel,
) -> Result<ValenceAssignment, ValenceError> {
    match model {
        ValenceModel::RdkitLike => {}
    }
    if any_unsupported_features(molecule) {
        return Err(ValenceError::NotImplemented);
    }

    let context = ValenceContext::from_molecule(molecule);
    let mut explicit_valence = vec![0u8; molecule.atoms().len()];
    let mut implicit_hydrogens = vec![0u8; molecule.atoms().len()];
    for (i, atom) in molecule.atoms().iter().enumerate() {
        let ev = calculate_explicit_valence_with_context(molecule, &context, i, true)?;
        let ih = calculate_implicit_valence_with_context(molecule, &context, i, ev, true)?;
        explicit_valence[i] = to_u8_checked(ev, i, atom)?;
        implicit_hydrogens[i] = to_u8_checked(ih, i, atom)?;
    }

    Ok(ValenceAssignment {
        explicit_valence,
        implicit_hydrogens,
    })
}
