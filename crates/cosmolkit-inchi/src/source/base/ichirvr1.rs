use crate::source_types::{
    BOND_TYPE_SINGLE, MAX_NUM_STEREO_BONDS, MAXVAL, NUM_H_ISOTOPES, RI_ERR_PROGR, RI_ERR_SYNTAX,
    SourceHeapError, inp_ATOM,
};

#[allow(non_snake_case)]
pub(crate) fn IncrZeroBonds(
    atoms: &mut [inp_ATOM],
    num_atoms: i32,
    component: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4172 IncrZeroBonds
    // INCHI✔️❌: void IncrZeroBonds( inp_ATOM *at, int num_at, int iComponent )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j;
    // INCHI✔️❌:     for (i = 0; i < num_at; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         at[i].component = iComponent;
    // INCHI✔️❌:         for (j = 0; j < at[i].valence; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!at[i].bond_type[j])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 at[i].bond_type[j] = BOND_TYPE_SINGLE;
    // INCHI✔️❌:                 at[i].chem_bonds_valence += BOND_TYPE_SINGLE;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: IncrZeroBonds
    // BEGIN ACTIVE INCHI MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/incomdef.h:67
    // INCHI✔️❌: #define BOND_TYPE_SINGLE    1
    // END ACTIVE INCHI MACRO

    let count =
        usize::try_from(num_atoms.max(0)).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    for index in 0..count {
        let atom = atoms
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        atom.component = component as u16;
        let valence = usize::try_from(i32::from(atom.valence).max(0))
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        for bond in 0..valence {
            let bond_type = atom
                .bond_type
                .get_mut(bond)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if *bond_type == 0 {
                *bond_type = BOND_TYPE_SINGLE as u8;
                atom.chem_bonds_valence =
                    atom.chem_bonds_valence.wrapping_add(BOND_TYPE_SINGLE as i8);
            }
        }
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn ClearEndpts(atoms: &mut [inp_ATOM], num_atoms: i32) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4191 ClearEndpts
    // INCHI✔️❌: void ClearEndpts( inp_ATOM *at, int num_at )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:     for (i = 0; i < num_at; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         at[i].endpoint = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ClearEndpts

    let count =
        usize::try_from(num_atoms.max(0)).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    for index in 0..count {
        atoms
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .endpoint = 0;
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn ConnectDisconnectedH(
    atoms: &mut [inp_ATOM],
    num_atoms: i32,
    num_deleted_h: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5480 ConnectDisconnectedH
    // INCHI✔️❌: int ConnectDisconnectedH( inp_ATOM *at, int num_atoms, int num_deleted_H )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, k, n, m, num_H;
    // INCHI✔️❌:     int tot_atoms = num_atoms + num_deleted_H;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = num_atoms; i < tot_atoms; i = j)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         k = at[i].neighbor[0]; /* a[k] is the atom connected to the explicit hydrogen at[i] */
    // INCHI✔️❌:
    // INCHI✔️❌:         for (j = i; j < tot_atoms && at[j].neighbor[0] == k; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         num_H = j - i; /* number of explicit H for at[k] */
    // INCHI✔️❌:         if (num_H > at[k].num_H)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return RI_ERR_PROGR;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (num_H + at[k].valence > MAXVAL)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return RI_ERR_SYNTAX;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* insert links to explicit H before all other links in the connection list */
    // INCHI✔️❌:         n = at[k].valence;
    // INCHI✔️❌:         memmove(at[k].neighbor + num_H, at[k].neighbor, sizeof(at[k].neighbor[0]) * n);
    // INCHI✔️❌:         memmove(at[k].bond_stereo + num_H, at[k].bond_stereo, sizeof(at[k].bond_stereo[0]) * n);
    // INCHI✔️❌:         memmove(at[k].bond_type + num_H, at[k].bond_type, sizeof(at[k].bond_type[0]) * n);
    // INCHI✔️❌:         for (n = 0; n < num_H; n++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             at[k].neighbor[n] = i + n;
    // INCHI✔️❌:             at[k].bond_stereo[n] = 0;
    // INCHI✔️❌:             at[k].bond_type[n] = BOND_TYPE_SINGLE;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         for (m = 0; m < MAX_NUM_STEREO_BONDS && at[k].sb_parity[m]; m++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             at[k].sb_ord[m] += num_H;
    // INCHI✔️❌:             if (at[k].sn_ord[m] < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (n = i; n < j; n++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (at[n].orig_at_number == at[k].sn_orig_at_num[m])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         at[k].sn_ord[m] = n - i;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (n == j)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return RI_ERR_PROGR;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 at[k].sn_ord[m] += num_H;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         at[k].valence += num_H;
    // INCHI✔️❌:         at[k].chem_bonds_valence += num_H;
    // INCHI✔️❌:         at[k].num_H -= num_H; /* cannot be negative */
    // INCHI✔️❌:
    // INCHI✔️❌:         /*memset( at[k].num_iso_H, 0, sizeof(at[0].num_iso_H) );*/ /* attached H must carry all isotopic shifts */
    // INCHI✔️❌:         for (n = i; n < j; n++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             at[n].chem_bonds_valence = BOND_TYPE_SINGLE;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* isotopic H */
    // INCHI✔️❌:         for (m = j - 1; i <= m && at[m].iso_atw_diff > 0; m--)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at[m].iso_atw_diff > NUM_H_ISOTOPES)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return RI_ERR_PROGR;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (0 >= at[k].num_iso_H[(int) at[m].iso_atw_diff - 1] --)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return RI_ERR_PROGR;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* subtract isotopic H */
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (m = 0; m < NUM_H_ISOTOPES; m++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             at[i].num_H -= at[i].num_iso_H[m];
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (0 > at[i].num_H)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return RI_ERR_PROGR;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return tot_atoms;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ConnectDisconnectedH
    // BEGIN ACTIVE INCHI MACROS
    // INCHI✔️❌: #define MAXVAL                20 /* max number of bonds per atom */
    // INCHI✔️❌: #define NUM_H_ISOTOPES        3  /* number of hydrogen isotopes: protium, deuterium, tritium */
    // INCHI✔️❌: #define BOND_TYPE_SINGLE       1
    // INCHI✔️❌: #define MAX_NUM_STEREO_BONDS   3
    // INCHI✔️❌: #define RI_ERR_SYNTAX          (-2)
    // INCHI✔️❌: #define RI_ERR_PROGR           (-3)
    // END ACTIVE INCHI MACROS

    let total_atoms = num_atoms
        .checked_add(num_deleted_h)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let atom_index = |value: i32| -> Result<usize, SourceHeapError> {
        usize::try_from(value).map_err(|_| SourceHeapError::PointerOutOfBounds)
    };

    let mut i = num_atoms;
    while i < total_atoms {
        let explicit_index = atom_index(i)?;
        let k = i32::from(
            atoms
                .get(explicit_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .neighbor[0],
        );
        let base_index = atom_index(k)?;

        let mut j = i;
        while j < total_atoms {
            let current = atoms
                .get(atom_index(j)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if i32::from(current.neighbor[0]) != k {
                break;
            }
            j += 1;
        }

        let num_h = j - i;
        let base = atoms
            .get(base_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if num_h > i32::from(base.num_H) {
            return Ok(RI_ERR_PROGR);
        }
        if num_h + i32::from(base.valence) > MAXVAL as i32 {
            return Ok(RI_ERR_SYNTAX);
        }

        let old_valence = usize::try_from(i32::from(base.valence))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let insert_count =
            usize::try_from(num_h).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let base = atoms
            .get_mut(base_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        base.neighbor.copy_within(0..old_valence, insert_count);
        base.bond_stereo.copy_within(0..old_valence, insert_count);
        base.bond_type.copy_within(0..old_valence, insert_count);
        for n in 0..insert_count {
            base.neighbor[n] = (i + n as i32) as u16;
            base.bond_stereo[n] = 0;
            base.bond_type[n] = BOND_TYPE_SINGLE as u8;
        }

        for m in 0..MAX_NUM_STEREO_BONDS as usize {
            if atoms[base_index].sb_parity[m] == 0 {
                break;
            }
            atoms[base_index].sb_ord[m] = atoms[base_index].sb_ord[m].wrapping_add(num_h as i8);
            if atoms[base_index].sn_ord[m] < 0 {
                let target = atoms[base_index].sn_orig_at_num[m];
                let mut n = i;
                while n < j {
                    if atoms
                        .get(atom_index(n)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .orig_at_number
                        == target
                    {
                        atoms[base_index].sn_ord[m] = (n - i) as i8;
                        break;
                    }
                    n += 1;
                }
                if n == j {
                    return Ok(RI_ERR_PROGR);
                }
            } else {
                atoms[base_index].sn_ord[m] = atoms[base_index].sn_ord[m].wrapping_add(num_h as i8);
            }
        }

        let base = atoms
            .get_mut(base_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        base.valence = base.valence.wrapping_add(num_h as i8);
        base.chem_bonds_valence = base.chem_bonds_valence.wrapping_add(num_h as i8);
        base.num_H = base.num_H.wrapping_sub(num_h as i8);

        let mut n = i;
        while n < j {
            atoms
                .get_mut(atom_index(n)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .chem_bonds_valence = BOND_TYPE_SINGLE as i8;
            n += 1;
        }

        let mut m = j - 1;
        while i <= m {
            let isotope_difference = atoms
                .get(atom_index(m)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .iso_atw_diff;
            if isotope_difference <= 0 {
                break;
            }
            if i32::from(isotope_difference) > NUM_H_ISOTOPES as i32 {
                return Ok(RI_ERR_PROGR);
            }
            let isotope_index = usize::try_from(i32::from(isotope_difference) - 1)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let count = &mut atoms
                .get_mut(base_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .num_iso_H[isotope_index];
            let old_count = *count;
            *count = count.wrapping_sub(1);
            if old_count <= 0 {
                return Ok(RI_ERR_PROGR);
            }
            m -= 1;
        }
        i = j;
    }

    let mut i = 0;
    while i < num_atoms {
        let atom = atoms
            .get_mut(atom_index(i)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for m in 0..NUM_H_ISOTOPES as usize {
            atom.num_H = atom.num_H.wrapping_sub(atom.num_iso_H[m]);
        }
        if atom.num_H < 0 {
            return Ok(RI_ERR_PROGR);
        }
        i += 1;
    }

    Ok(total_atoms)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_port__ichirvr1__incrzerobonds__line_4172() {
        let untouched = inp_ATOM {
            component: 9,
            valence: 1,
            bond_type: {
                let mut value = [0; 20];
                value[0] = 3;
                value
            },
            chem_bonds_valence: 3,
            ..inp_ATOM::default()
        };
        let mut atoms = vec![untouched.clone()];
        assert_eq!(IncrZeroBonds(&mut atoms, -1, 7), Ok(()));
        assert_eq!(atoms, vec![untouched.clone()]);
        assert_eq!(IncrZeroBonds(&mut atoms, 0, 7), Ok(()));
        assert_eq!(atoms, vec![untouched]);

        let mut first = inp_ATOM {
            component: 8,
            endpoint: 17,
            valence: 4,
            chem_bonds_valence: 126,
            ..inp_ATOM::default()
        };
        first.bond_type[..4].copy_from_slice(&[0, 1, 2, 0]);
        let mut second = inp_ATOM {
            component: 8,
            endpoint: 19,
            valence: -1,
            chem_bonds_valence: -7,
            ..inp_ATOM::default()
        };
        second.bond_type[0] = 0;
        let third = inp_ATOM {
            component: 21,
            valence: 1,
            ..inp_ATOM::default()
        };
        let mut atoms = vec![first, second, third.clone()];
        assert_eq!(IncrZeroBonds(&mut atoms, 2, -1), Ok(()));
        assert_eq!(atoms[0].component, u16::MAX);
        assert_eq!(atoms[0].endpoint, 17);
        assert_eq!(&atoms[0].bond_type[..4], &[1, 1, 2, 1]);
        assert_eq!(atoms[0].chem_bonds_valence, i8::MIN);
        assert_eq!(atoms[1].component, u16::MAX);
        assert_eq!(atoms[1].endpoint, 19);
        assert_eq!(atoms[1].bond_type[0], 0);
        assert_eq!(atoms[1].chem_bonds_valence, -7);
        assert_eq!(atoms[2], third);

        assert_eq!(IncrZeroBonds(&mut atoms, 1, 65_536), Ok(()));
        assert_eq!(atoms[0].component, 0);
        assert_eq!(&atoms[0].bond_type[..4], &[1, 1, 2, 1]);
        assert_eq!(atoms[0].chem_bonds_valence, i8::MIN);
    }

    #[test]
    fn source_port__ichirvr1__clearendpts__line_4191() {
        let original = vec![
            inp_ATOM {
                endpoint: 1,
                component: 11,
                charge: -3,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                endpoint: u16::MAX,
                component: 13,
                charge: 4,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                endpoint: 17,
                component: 19,
                charge: 5,
                ..inp_ATOM::default()
            },
        ];
        let mut atoms = original.clone();
        assert_eq!(ClearEndpts(&mut atoms, -1), Ok(()));
        assert_eq!(atoms, original);
        assert_eq!(ClearEndpts(&mut atoms, 0), Ok(()));
        assert_eq!(atoms, original);

        assert_eq!(ClearEndpts(&mut atoms, 2), Ok(()));
        assert_eq!(atoms[0].endpoint, 0);
        assert_eq!(atoms[1].endpoint, 0);
        assert_eq!(atoms[2].endpoint, 17);
        assert_eq!(atoms[0].component, 11);
        assert_eq!(atoms[1].component, 13);
        assert_eq!(atoms[2].component, 19);
        assert_eq!(atoms[0].charge, -3);
        assert_eq!(atoms[1].charge, 4);
        assert_eq!(atoms[2].charge, 5);

        assert_eq!(ClearEndpts(&mut atoms, 3), Ok(()));
        assert!(atoms.iter().all(|atom| atom.endpoint == 0));
    }

    #[test]
    fn source_port__ichirvr1__connectdisconnectedh__line_5480() {
        let mut base = inp_ATOM {
            valence: 1,
            chem_bonds_valence: 2,
            num_H: 5,
            num_iso_H: [1, 1, 1],
            ..inp_ATOM::default()
        };
        base.neighbor[0] = 9;
        base.bond_stereo[0] = -4;
        base.bond_type[0] = 2;
        base.sb_parity[..2].copy_from_slice(&[1, 2]);
        base.sb_ord[..2].copy_from_slice(&[0, 2]);
        base.sn_ord[..2].copy_from_slice(&[-1, 0]);
        base.sn_orig_at_num[0] = 102;
        let mut atoms = vec![
            base,
            inp_ATOM {
                neighbor: {
                    let mut value = [0; MAXVAL as usize];
                    value[0] = 0;
                    value
                },
                orig_at_number: 101,
                iso_atw_diff: 0,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                neighbor: {
                    let mut value = [0; MAXVAL as usize];
                    value[0] = 0;
                    value
                },
                orig_at_number: 102,
                iso_atw_diff: 2,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                neighbor: {
                    let mut value = [0; MAXVAL as usize];
                    value[0] = 0;
                    value
                },
                orig_at_number: 103,
                iso_atw_diff: 3,
                ..inp_ATOM::default()
            },
        ];
        assert_eq!(ConnectDisconnectedH(&mut atoms, 1, 3), Ok(4));
        assert_eq!(&atoms[0].neighbor[..4], &[1, 2, 3, 9]);
        assert_eq!(&atoms[0].bond_stereo[..4], &[0, 0, 0, -4]);
        assert_eq!(&atoms[0].bond_type[..4], &[1, 1, 1, 2]);
        assert_eq!(&atoms[0].sb_ord[..2], &[3, 5]);
        assert_eq!(&atoms[0].sn_ord[..2], &[1, 3]);
        assert_eq!(atoms[0].valence, 4);
        assert_eq!(atoms[0].chem_bonds_valence, 5);
        assert_eq!(atoms[0].num_iso_H, [1, 0, 0]);
        assert_eq!(atoms[0].num_H, 1);
        assert!(atoms[1..].iter().all(|atom| atom.chem_bonds_valence == 1));

        let mut isotope_only = vec![inp_ATOM {
            num_H: 6,
            num_iso_H: [1, 2, 3],
            ..inp_ATOM::default()
        }];
        assert_eq!(ConnectDisconnectedH(&mut isotope_only, 1, 0), Ok(1));
        assert_eq!(isotope_only[0].num_H, 0);

        let mut too_many = vec![
            inp_ATOM {
                num_H: 0,
                ..inp_ATOM::default()
            },
            inp_ATOM::default(),
        ];
        assert_eq!(ConnectDisconnectedH(&mut too_many, 1, 1), Ok(RI_ERR_PROGR));
        assert_eq!(too_many[0].valence, 0);

        let mut maxval = vec![
            inp_ATOM {
                valence: MAXVAL as i8,
                num_H: 1,
                ..inp_ATOM::default()
            },
            inp_ATOM::default(),
        ];
        assert_eq!(ConnectDisconnectedH(&mut maxval, 1, 1), Ok(RI_ERR_SYNTAX));
        assert_eq!(maxval[0].valence, MAXVAL as i8);

        let mut missing_stereo_h = vec![
            inp_ATOM {
                neighbor: {
                    let mut value = [0; MAXVAL as usize];
                    value[0] = 7;
                    value
                },
                bond_stereo: {
                    let mut value = [0; MAXVAL as usize];
                    value[0] = 6;
                    value
                },
                bond_type: {
                    let mut value = [0; MAXVAL as usize];
                    value[0] = 2;
                    value
                },
                valence: 1,
                num_H: 1,
                sb_ord: [4, 0, 0],
                sn_ord: [-1, 0, 0],
                sb_parity: [1, 0, 0],
                sn_orig_at_num: [999, 0, 0],
                ..inp_ATOM::default()
            },
            inp_ATOM {
                orig_at_number: 100,
                ..inp_ATOM::default()
            },
        ];
        assert_eq!(
            ConnectDisconnectedH(&mut missing_stereo_h, 1, 1),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(&missing_stereo_h[0].neighbor[..2], &[1, 7]);
        assert_eq!(&missing_stereo_h[0].bond_stereo[..2], &[0, 6]);
        assert_eq!(&missing_stereo_h[0].bond_type[..2], &[1, 2]);
        assert_eq!(missing_stereo_h[0].sb_ord[0], 5);
        assert_eq!(missing_stereo_h[0].valence, 1);

        let mut isotope_too_large = vec![
            inp_ATOM {
                num_H: 1,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                iso_atw_diff: (NUM_H_ISOTOPES + 1) as i8,
                ..inp_ATOM::default()
            },
        ];
        assert_eq!(
            ConnectDisconnectedH(&mut isotope_too_large, 1, 1),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(isotope_too_large[0].valence, 1);
        assert_eq!(isotope_too_large[1].chem_bonds_valence, 1);

        let mut isotope_underflow = vec![
            inp_ATOM {
                num_H: 1,
                num_iso_H: [0, 0, 0],
                ..inp_ATOM::default()
            },
            inp_ATOM {
                iso_atw_diff: 1,
                ..inp_ATOM::default()
            },
        ];
        assert_eq!(
            ConnectDisconnectedH(&mut isotope_underflow, 1, 1),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(isotope_underflow[0].num_iso_H[0], -1);

        let mut final_underflow = vec![inp_ATOM {
            num_H: 1,
            num_iso_H: [1, 1, 0],
            ..inp_ATOM::default()
        }];
        assert_eq!(
            ConnectDisconnectedH(&mut final_underflow, 1, 0),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(final_underflow[0].num_H, -1);
    }
}
