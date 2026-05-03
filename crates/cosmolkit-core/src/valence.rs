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
    // Complete RDKit valence cache (Release_2026_03_1), generated via:
    //   [list(Chem.GetPeriodicTable().GetValenceList(z)) for z in range(119)]
    // This intentionally mirrors RDKit's table entries and ordering.
    match atomic_num {
        0 => Some(&[-1]),
        1 => Some(&[1]),
        2 => Some(&[0]),
        3 => Some(&[1]),
        4 => Some(&[2]),
        5 => Some(&[3]),
        6 => Some(&[4]),
        7 => Some(&[3]),
        8 => Some(&[2]),
        9 => Some(&[1]),
        10 => Some(&[0]),
        11 => Some(&[1]),
        12 => Some(&[2, -1]),
        13 => Some(&[3, 6]),
        14 => Some(&[4, 6]),
        15 => Some(&[3, 5, 7]),
        16 => Some(&[2, 4, 6]),
        17 => Some(&[1]),
        18 => Some(&[0]),
        19 => Some(&[1]),
        20 => Some(&[2, -1]),
        21 => Some(&[-1]),
        22 => Some(&[-1]),
        23 => Some(&[-1]),
        24 => Some(&[-1]),
        25 => Some(&[-1]),
        26 => Some(&[-1]),
        27 => Some(&[-1]),
        28 => Some(&[-1]),
        29 => Some(&[-1]),
        30 => Some(&[-1]),
        31 => Some(&[3]),
        32 => Some(&[4]),
        33 => Some(&[3, 5, 7]),
        34 => Some(&[2, 4, 6]),
        35 => Some(&[1]),
        36 => Some(&[0]),
        37 => Some(&[1]),
        38 => Some(&[2, -1]),
        39 => Some(&[-1]),
        40 => Some(&[-1]),
        41 => Some(&[-1]),
        42 => Some(&[-1]),
        43 => Some(&[-1]),
        44 => Some(&[-1]),
        45 => Some(&[-1]),
        46 => Some(&[-1]),
        47 => Some(&[-1]),
        48 => Some(&[-1]),
        49 => Some(&[3]),
        50 => Some(&[2, 4]),
        51 => Some(&[3, 5, 7]),
        52 => Some(&[2, 4, 6]),
        53 => Some(&[1, 3, 5]),
        54 => Some(&[0, 2, 4, 6]),
        55 => Some(&[1]),
        56 => Some(&[2, -1]),
        57 => Some(&[-1]),
        58 => Some(&[-1]),
        59 => Some(&[-1]),
        60 => Some(&[-1]),
        61 => Some(&[-1]),
        62 => Some(&[-1]),
        63 => Some(&[-1]),
        64 => Some(&[-1]),
        65 => Some(&[-1]),
        66 => Some(&[-1]),
        67 => Some(&[-1]),
        68 => Some(&[-1]),
        69 => Some(&[-1]),
        70 => Some(&[-1]),
        71 => Some(&[-1]),
        72 => Some(&[-1]),
        73 => Some(&[-1]),
        74 => Some(&[-1]),
        75 => Some(&[-1]),
        76 => Some(&[-1]),
        77 => Some(&[-1]),
        78 => Some(&[-1]),
        79 => Some(&[-1]),
        80 => Some(&[-1]),
        81 => Some(&[3]),
        82 => Some(&[2, 4]),
        83 => Some(&[3, 5, 7]),
        84 => Some(&[2, 4, 6]),
        85 => Some(&[1, 3, 5]),
        86 => Some(&[0]),
        87 => Some(&[1]),
        88 => Some(&[2, -1]),
        89 => Some(&[-1]),
        90 => Some(&[-1]),
        91 => Some(&[-1]),
        92 => Some(&[-1]),
        93 => Some(&[-1]),
        94 => Some(&[-1]),
        95 => Some(&[-1]),
        96 => Some(&[-1]),
        97 => Some(&[-1]),
        98 => Some(&[-1]),
        99 => Some(&[-1]),
        100 => Some(&[-1]),
        101 => Some(&[-1]),
        102 => Some(&[-1]),
        103 => Some(&[-1]),
        104 => Some(&[-1]),
        105 => Some(&[-1]),
        106 => Some(&[-1]),
        107 => Some(&[-1]),
        108 => Some(&[-1]),
        109 => Some(&[-1]),
        110 => Some(&[-1]),
        111 => Some(&[-1]),
        112 => Some(&[-1]),
        113 => Some(&[-1]),
        114 => Some(&[-1]),
        115 => Some(&[-1]),
        116 => Some(&[-1]),
        117 => Some(&[-1]),
        118 => Some(&[-1]),
        _ => None,
    }
}

/// Returns the RDKit 2026.03.1 valence list for an atomic number.
pub fn rdkit_valence_list(atomic_num: u8) -> Option<&'static [i32]> {
    valence_list(atomic_num)
}

fn is_aromatic_atom(molecule: &Molecule, atom_index: usize) -> bool {
    if molecule.atoms()[atom_index].is_aromatic {
        return true;
    }
    molecule
        .bonds()
        .iter()
        .any(|b| (b.begin_atom == atom_index || b.end_atom == atom_index) && b.is_aromatic)
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

    let mut accum = atom.explicit_hydrogens as f64;
    for b in molecule.bonds() {
        accum += bond_valence_contrib_for_atom(b, atom_index);
    }

    if accum > dv as f64 && is_aromatic_atom(molecule, atom_index) {
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

fn calculate_implicit_valence(
    molecule: &Molecule,
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

    if is_aromatic_atom(molecule, atom_index) {
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
    match atomic_num {
        0 => Some(0),
        1 => Some(1),
        2 => Some(2),
        3 => Some(1),
        4 => Some(2),
        5 => Some(3),
        6 => Some(4),
        7 => Some(5),
        8 => Some(6),
        9 => Some(7),
        10 => Some(8),
        11 => Some(1),
        12 => Some(2),
        13 => Some(3),
        14 => Some(4),
        15 => Some(5),
        16 => Some(6),
        17 => Some(7),
        18 => Some(8),
        19 => Some(1),
        20 => Some(2),
        21 => Some(3),
        22 => Some(4),
        23 => Some(5),
        24 => Some(6),
        25 => Some(7),
        26 => Some(8),
        27 => Some(9),
        28 => Some(10),
        29 => Some(11),
        30 => Some(2),
        31 => Some(3),
        32 => Some(4),
        33 => Some(5),
        34 => Some(6),
        35 => Some(7),
        36 => Some(8),
        37 => Some(1),
        38 => Some(2),
        39 => Some(3),
        40 => Some(4),
        41 => Some(5),
        42 => Some(6),
        43 => Some(7),
        44 => Some(8),
        45 => Some(9),
        46 => Some(10),
        47 => Some(11),
        48 => Some(2),
        49 => Some(3),
        50 => Some(4),
        51 => Some(5),
        52 => Some(6),
        53 => Some(7),
        54 => Some(8),
        55 => Some(1),
        56 => Some(2),
        57 => Some(3),
        58 => Some(4),
        59 => Some(3),
        60 => Some(4),
        61 => Some(5),
        62 => Some(6),
        63 => Some(7),
        64 => Some(8),
        65 => Some(9),
        66 => Some(10),
        67 => Some(11),
        68 => Some(12),
        69 => Some(13),
        70 => Some(14),
        71 => Some(15),
        72 => Some(4),
        73 => Some(5),
        74 => Some(6),
        75 => Some(7),
        76 => Some(8),
        77 => Some(9),
        78 => Some(10),
        79 => Some(11),
        80 => Some(2),
        81 => Some(3),
        82 => Some(4),
        83 => Some(5),
        84 => Some(6),
        85 => Some(7),
        86 => Some(8),
        87 => Some(1),
        88 => Some(2),
        89 => Some(3),
        90 => Some(4),
        91 => Some(3),
        92 => Some(4),
        93 => Some(5),
        94 => Some(6),
        95 => Some(7),
        96 => Some(8),
        97 => Some(9),
        98 => Some(10),
        99 => Some(11),
        100 => Some(12),
        101 => Some(13),
        102 => Some(14),
        103 => Some(15),
        104 => Some(2),
        105 => Some(2),
        106 => Some(2),
        107 => Some(2),
        108 => Some(2),
        109 => Some(2),
        110 => Some(2),
        111 => Some(2),
        112 => Some(2),
        113 => Some(2),
        114 => Some(2),
        115 => Some(2),
        116 => Some(2),
        117 => Some(2),
        118 => Some(2),
        _ => None,
    }
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

    for (i, atom) in molecule.atoms().iter().enumerate() {
        if !atom.no_implicit || atom.atomic_num == 0 {
            continue;
        }
        let valens = valence_list(atom.atomic_num).ok_or(ValenceError::NotImplemented)?;
        let chg = atom.formal_charge as i32;
        let n_outer = n_outer_electrons(atom.atomic_num).ok_or(ValenceError::NotImplemented)?;
        let value = if valens.len() != 1 || valens[0] != -1 {
            let total_valence = if is_aromatic_atom(molecule, i) {
                // RDKit runs assignRadicals after Kekulize(). In this codebase we
                // do not yet store a separate kekulized bond table, so we use the
                // already-RDKit-aligned explicit valence cache as the sanitized
                // total valence surrogate for aromatic atoms.
                existing_explicit_valence[i] as i32
            } else {
                let mut accum = atom.explicit_hydrogens as f64;
                for b in molecule.bonds() {
                    accum += bond_valence_contrib_for_atom(b, i);
                }
                (accum + 0.1) as i32
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
            let degree = molecule
                .bonds()
                .iter()
                .filter(|b| b.begin_atom == i || b.end_atom == i)
                .count();
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

    let mut explicit_valence = vec![0u8; molecule.atoms().len()];
    let mut implicit_hydrogens = vec![0u8; molecule.atoms().len()];
    for (i, atom) in molecule.atoms().iter().enumerate() {
        let ev = calculate_explicit_valence(molecule, i, true)?;
        let ih = calculate_implicit_valence(molecule, i, ev, true)?;
        explicit_valence[i] = to_u8_checked(ev, i, atom)?;
        implicit_hydrogens[i] = to_u8_checked(ih, i, atom)?;
    }

    Ok(ValenceAssignment {
        explicit_valence,
        implicit_hydrogens,
    })
}
