use crate::source_types::{
    AT_FLAG_ISO_H_POINT, AT_ISO_SORT_KEY, AT_ISO_SORT_KEY_MULT, SourceHeap, SourceHeapError, SourceMutPointer,
    T_GROUP_INFO, sp_ATOM,
};

#[allow(non_snake_case)]
pub(crate) fn make_iso_sort_key(iso_atw_diff: i32, num_1H: i32, num_2H: i32, num_3H: i32) -> AT_ISO_SORT_KEY {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiisot.c:47 make_iso_sort_key
    // INCHI✔️✔️: AT_ISO_SORT_KEY make_iso_sort_key( int iso_atw_diff, int num_1H, int num_2H, int num_3H )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     AT_ISO_SORT_KEY iso_sort_key = 0, mult = 1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     iso_sort_key += mult * num_1H;
    // INCHI✔️✔️:     mult *= AT_ISO_SORT_KEY_MULT;
    // INCHI✔️✔️:     iso_sort_key += mult * num_2H;
    // INCHI✔️✔️:     mult *= AT_ISO_SORT_KEY_MULT;
    // INCHI✔️✔️:     iso_sort_key += mult * num_3H;
    // INCHI✔️✔️:     mult *= AT_ISO_SORT_KEY_MULT;
    // INCHI✔️✔️:     iso_sort_key += mult * iso_atw_diff;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return iso_sort_key;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: make_iso_sort_key
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: make_iso_sort_key
    // INCHI✔️✔️: #define AT_ISO_SORT_KEY_MULT 32
    // INCHI✔️✔️: GCC/Linux LP64: typedef long AT_ISO_SORT_KEY; sizeof(long) == 8
    // END INCHI ACTIVE MACRO CONFIGURATION: make_iso_sort_key

    let mut iso_sort_key: AT_ISO_SORT_KEY = 0;
    let mut mult: AT_ISO_SORT_KEY = 1;
    iso_sort_key += mult * AT_ISO_SORT_KEY::from(num_1H);
    mult *= AT_ISO_SORT_KEY::from(AT_ISO_SORT_KEY_MULT);
    iso_sort_key += mult * AT_ISO_SORT_KEY::from(num_2H);
    mult *= AT_ISO_SORT_KEY::from(AT_ISO_SORT_KEY_MULT);
    iso_sort_key += mult * AT_ISO_SORT_KEY::from(num_3H);
    mult *= AT_ISO_SORT_KEY::from(AT_ISO_SORT_KEY_MULT);
    iso_sort_key += mult * AT_ISO_SORT_KEY::from(iso_atw_diff);
    iso_sort_key
}

#[allow(non_snake_case)]
pub(crate) fn set_atom_iso_sort_keys(
    heap: &mut SourceHeap,
    num_at: i32,
    at: SourceMutPointer<sp_ATOM>,
    t_group_info: Option<&T_GROUP_INFO>,
    mut bHasIsotopicInTautomerGroups: Option<&mut i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiisot.c:67 set_atom_iso_sort_keys
    // INCHI✔️✔️: int set_atom_iso_sort_keys( int num_at,
    // INCHI✔️✔️:                             sp_ATOM *at,
    // INCHI✔️✔️:                             T_GROUP_INFO* t_group_info,
    // INCHI✔️✔️:                             int *bHasIsotopicInTautomerGroups )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int             i, num_isotopic = 0, bMergedTgroup;
    // INCHI✔️✔️:     AT_ISO_SORT_KEY iso_sort_key;
    // INCHI✔️✔️:     T_GROUP        *t_group =
    // INCHI✔️✔️:         ( t_group_info &&
    // INCHI✔️✔️:          t_group_info->t_group &&
    // INCHI✔️✔️:          t_group_info->num_t_groups > 0 ) ? t_group_info->t_group : NULL;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (bHasIsotopicInTautomerGroups)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *bHasIsotopicInTautomerGroups = 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     for (i = 0; i < num_at; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         bMergedTgroup = ( t_group_info && t_group_info->nIsotopicEndpointAtomNumber && ( at[i].cFlags & AT_FLAG_ISO_H_POINT ) );
    // INCHI✔️✔️:         if (( !at[i].endpoint || !t_group ) && !bMergedTgroup)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             iso_sort_key = make_iso_sort_key( at[i].iso_atw_diff, at[i].num_iso_H[0], at[i].num_iso_H[1], at[i].num_iso_H[2] );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:              /*  H isotopes go to the tautomer part of the CT (name) */
    // INCHI✔️✔️:              /*  if (at[i].endpoint && t_group) ... */
    // INCHI✔️✔️:             iso_sort_key = make_iso_sort_key( at[i].iso_atw_diff, 0, 0, 0 );
    // INCHI✔️✔️:             if (bHasIsotopicInTautomerGroups)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 *bHasIsotopicInTautomerGroups += ( at[i].num_iso_H[0] || at[i].num_iso_H[1] || at[i].num_iso_H[2] || bMergedTgroup );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         at[i].iso_sort_key = iso_sort_key;
    // INCHI✔️✔️:         num_isotopic += ( iso_sort_key != 0 );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return num_isotopic;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: set_atom_iso_sort_keys
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: set_atom_iso_sort_keys
    // INCHI✔️✔️: #define AT_FLAG_ISO_H_POINT 0x01
    // INCHI✔️✔️: #define AT_ISO_SORT_KEY_MULT 32
    // INCHI✔️✔️: GCC/Linux LP64: typedef long AT_ISO_SORT_KEY; sizeof(long) == 8
    // END INCHI ACTIVE MACRO CONFIGURATION: set_atom_iso_sort_keys

    if let Some(output) = bHasIsotopicInTautomerGroups.as_deref_mut() {
        *output = 0;
    }

    let has_t_group = t_group_info.is_some_and(|info| !info.t_group.is_null() && info.num_t_groups > 0);
    let has_isotopic_endpoint_numbers = t_group_info.is_some_and(|info| !info.nIsotopicEndpointAtomNumber.is_null());
    let mut atoms = if num_at > 0 { Some(heap.slice_mut(at)?) } else { None };
    let mut num_isotopic = 0_i32;
    let mut i = 0_i32;
    while i < num_at {
        let atom = atoms
            .as_deref_mut()
            .and_then(|atoms| atoms.get_mut(i as usize))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let merged_t_group = has_isotopic_endpoint_numbers && (atom.cFlags & AT_FLAG_ISO_H_POINT as i8) != 0;
        let iso_sort_key = if (atom.endpoint == 0 || !has_t_group) && !merged_t_group {
            make_iso_sort_key(
                i32::from(atom.iso_atw_diff),
                i32::from(atom.num_iso_H[0]),
                i32::from(atom.num_iso_H[1]),
                i32::from(atom.num_iso_H[2]),
            )
        } else {
            if let Some(output) = bHasIsotopicInTautomerGroups.as_deref_mut() {
                *output = output.wrapping_add(i32::from(
                    atom.num_iso_H[0] != 0 || atom.num_iso_H[1] != 0 || atom.num_iso_H[2] != 0 || merged_t_group,
                ));
            }
            make_iso_sort_key(i32::from(atom.iso_atw_diff), 0, 0, 0)
        };
        atom.iso_sort_key = iso_sort_key;
        num_isotopic = num_isotopic.wrapping_add(i32::from(iso_sort_key != 0));
        i = i.wrapping_add(1);
    }

    Ok(num_isotopic)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{SourceMutPointer, T_GROUP};

    #[test]
    fn source_port__ichiisot__make_iso_sort_key__line_47() {
        assert_eq!(make_iso_sort_key(0, 0, 0, 0), 0);
        assert_eq!(make_iso_sort_key(0, 1, 0, 0), 1);
        assert_eq!(make_iso_sort_key(0, 0, 1, 0), 32);
        assert_eq!(make_iso_sort_key(0, 0, 0, 1), 1_024);
        assert_eq!(make_iso_sort_key(1, 0, 0, 0), 32_768);
        assert_eq!(make_iso_sort_key(-1, 3, 2, 1), -31_677);
        assert_eq!(make_iso_sort_key(7, 4, 3, 2), 231_524);

        assert_eq!(
            make_iso_sort_key(i32::MAX, i32::MAX, i32::MAX, i32::MAX),
            72_638_634_359_775
        );
        assert_eq!(
            make_iso_sort_key(i32::MIN, i32::MIN, i32::MIN, i32::MIN),
            -72_638_634_393_600
        );
        assert_eq!(
            make_iso_sort_key(i32::MIN, i32::MAX, i32::MIN, i32::MAX),
            -68_236_292_916_225
        );
    }

    #[test]
    fn source_port__ichiisot__set_atom_iso_sort_keys__line_67() {
        let mut empty_heap = SourceHeap::default();
        let mut empty_taut_count = 91;
        assert_eq!(
            set_atom_iso_sort_keys(
                &mut empty_heap,
                0,
                SourceMutPointer::null(),
                None,
                Some(&mut empty_taut_count),
            ),
            Ok(0)
        );
        assert_eq!(empty_taut_count, 0);
        empty_taut_count = 92;
        assert_eq!(
            set_atom_iso_sort_keys(
                &mut empty_heap,
                -1,
                SourceMutPointer::null(),
                None,
                Some(&mut empty_taut_count),
            ),
            Ok(0)
        );
        assert_eq!(empty_taut_count, 0);

        let mut heap = SourceHeap::default();
        let mut atoms = vec![sp_ATOM {
            iso_sort_key: 111,
            ..sp_ATOM::default()
        }];
        atoms.extend([
            sp_ATOM {
                iso_atw_diff: 0,
                num_iso_H: [1, 2, 3],
                iso_sort_key: 112,
                ..sp_ATOM::default()
            },
            sp_ATOM {
                endpoint: 1,
                iso_atw_diff: -2,
                num_iso_H: [4, 5, 6],
                iso_sort_key: 113,
                ..sp_ATOM::default()
            },
            sp_ATOM {
                cFlags: AT_FLAG_ISO_H_POINT as i8,
                iso_atw_diff: 7,
                num_iso_H: [0, 0, 0],
                iso_sort_key: 114,
                ..sp_ATOM::default()
            },
            sp_ATOM {
                endpoint: 2,
                iso_atw_diff: i8::MIN,
                num_iso_H: [i8::MAX, i8::MIN, -1],
                iso_sort_key: 115,
                ..sp_ATOM::default()
            },
        ]);
        let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
        let group_pointer = heap.allocate_model_storage(vec![T_GROUP::default()]).unwrap();
        let isotope_endpoints = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let active_info = T_GROUP_INFO {
            t_group: group_pointer,
            num_t_groups: 1,
            nIsotopicEndpointAtomNumber: isotope_endpoints,
            ..T_GROUP_INFO::default()
        };
        let mut taut_count = 99;
        assert_eq!(
            set_atom_iso_sort_keys(
                &mut heap,
                4,
                atom_pointer.offset(1).unwrap(),
                Some(&active_info),
                Some(&mut taut_count),
            ),
            Ok(4)
        );
        assert_eq!(taut_count, 3);
        let actual = heap.slice(atom_pointer.as_const()).unwrap();
        assert_eq!(actual[0].iso_sort_key, 111);
        assert_eq!(actual[1].iso_sort_key, make_iso_sort_key(0, 1, 2, 3));
        assert_eq!(actual[2].iso_sort_key, make_iso_sort_key(-2, 0, 0, 0));
        assert_eq!(actual[3].iso_sort_key, make_iso_sort_key(7, 0, 0, 0));
        assert_eq!(actual[4].iso_sort_key, make_iso_sort_key(-128, 0, 0, 0));

        let no_active_group_cases = [
            T_GROUP_INFO::default(),
            T_GROUP_INFO {
                t_group: group_pointer,
                num_t_groups: 0,
                ..T_GROUP_INFO::default()
            },
            T_GROUP_INFO {
                num_t_groups: 1,
                ..T_GROUP_INFO::default()
            },
        ];
        for info in &no_active_group_cases {
            let pointer = heap
                .allocate_model_storage(vec![sp_ATOM {
                    endpoint: 1,
                    iso_atw_diff: 1,
                    num_iso_H: [2, 3, 4],
                    iso_sort_key: -1,
                    ..sp_ATOM::default()
                }])
                .unwrap();
            let mut output = -1;
            assert_eq!(
                set_atom_iso_sort_keys(&mut heap, 1, pointer, Some(info), Some(&mut output)),
                Ok(1)
            );
            assert_eq!(output, 0);
            assert_eq!(
                heap.slice(pointer.as_const()).unwrap()[0].iso_sort_key,
                make_iso_sort_key(1, 2, 3, 4)
            );
        }

        let pointer = heap
            .allocate_model_storage(vec![sp_ATOM {
                cFlags: AT_FLAG_ISO_H_POINT as i8,
                num_iso_H: [1, 0, 0],
                iso_sort_key: -1,
                ..sp_ATOM::default()
            }])
            .unwrap();
        assert_eq!(
            set_atom_iso_sort_keys(&mut heap, 1, pointer, Some(&active_info), None),
            Ok(0)
        );
        assert_eq!(heap.slice(pointer.as_const()).unwrap()[0].iso_sort_key, 0);
    }
}
