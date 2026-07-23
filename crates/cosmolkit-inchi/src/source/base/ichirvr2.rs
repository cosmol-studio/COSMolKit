use crate::source_types::{SourceHeapError, inp_ATOM, inp_ATOM_STEREO};

#[allow(non_snake_case)]
pub(crate) fn CopySt2At(
    atoms: &mut [inp_ATOM],
    stereo: Option<&[inp_ATOM_STEREO]>,
    num_atoms: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:79 CopySt2At
    // INCHI✔️❌: void CopySt2At( inp_ATOM *at, inp_ATOM_STEREO * st, int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:     if (!st)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (st[i].p_parity)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memcpy(at[i].p_orig_at_num, st[i].p_orig_at_num, sizeof(at[0].p_orig_at_num));
    // INCHI✔️❌:             at[i].p_parity = st[i].p_parity;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (st[i].sb_parity[0])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memcpy(at[i].sb_ord, st[i].sb_ord, sizeof(st[0].sb_ord));
    // INCHI✔️❌:             memcpy(at[i].sb_parity, st[i].sb_parity, sizeof(at[0].sb_parity));
    // INCHI✔️❌:             memcpy(at[i].sn_ord, st[i].sn_ord, sizeof(at[0].sn_ord));
    // INCHI✔️❌:             memcpy(at[i].sn_orig_at_num, st[i].sn_orig_at_num, sizeof(at[0].sn_orig_at_num));
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CopySt2At

    let Some(stereo) = stereo else {
        return Ok(());
    };
    let count =
        usize::try_from(num_atoms.max(0)).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    for index in 0..count {
        let source = stereo
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let target = atoms
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if source.p_parity != 0 {
            target.p_orig_at_num = source.p_orig_at_num;
            target.p_parity = source.p_parity;
        }
        if source.sb_parity[0] != 0 {
            target.sb_ord = source.sb_ord;
            target.sb_parity = source.sb_parity;
            target.sn_ord = source.sn_ord;
            target.sn_orig_at_num = source.sn_orig_at_num;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_port__ichirvr2__copyst2at__line_79() {
        let original = inp_ATOM {
            p_parity: 9,
            p_orig_at_num: [91, 92, 93, 94],
            sb_ord: [11, 12, 13],
            sb_parity: [14, 15, 16],
            sn_ord: [17, 18, 19],
            sn_orig_at_num: [20, 21, 22],
            charge: -3,
            ..inp_ATOM::default()
        };
        let mut atoms = vec![original.clone()];
        assert_eq!(CopySt2At(&mut atoms, None, 1), Ok(()));
        assert_eq!(atoms, vec![original.clone()]);

        let ignored = inp_ATOM_STEREO {
            p_parity: 0,
            p_orig_at_num: [1, 2, 3, 4],
            sb_ord: [5, 6, 7],
            sb_parity: [0, 8, 9],
            sn_ord: [-1, -2, -3],
            sn_orig_at_num: [10, 11, 12],
            ..inp_ATOM_STEREO::default()
        };
        assert_eq!(CopySt2At(&mut atoms, Some(&[ignored]), -1), Ok(()));
        assert_eq!(atoms, vec![original.clone()]);
        assert_eq!(
            CopySt2At(
                &mut atoms,
                Some(&[inp_ATOM_STEREO {
                    p_parity: 0,
                    sb_parity: [0, 8, 9],
                    ..inp_ATOM_STEREO::default()
                }]),
                1,
            ),
            Ok(())
        );
        assert_eq!(atoms, vec![original.clone()]);

        let tetrahedral_only = inp_ATOM_STEREO {
            bUsed0DParity: 7,
            p_parity: -4,
            p_orig_at_num: [1, 2, 3, 4],
            sb_ord: [31, 32, 33],
            sb_parity: [0, 34, 35],
            sn_ord: [36, 37, 38],
            sn_orig_at_num: [39, 40, 41],
        };
        let bond_only = inp_ATOM_STEREO {
            bUsed0DParity: 8,
            p_parity: 0,
            p_orig_at_num: [51, 52, 53, 54],
            sb_ord: [-1, -2, -3],
            sb_parity: [1, 2, 3],
            sn_ord: [4, 5, 6],
            sn_orig_at_num: [61, 62, 63],
        };
        let mut atoms = vec![original.clone(), original.clone(), original.clone()];
        assert_eq!(
            CopySt2At(&mut atoms, Some(&[tetrahedral_only, bond_only]), 2),
            Ok(())
        );
        assert_eq!(atoms[0].p_parity, -4);
        assert_eq!(atoms[0].p_orig_at_num, [1, 2, 3, 4]);
        assert_eq!(atoms[0].sb_ord, original.sb_ord);
        assert_eq!(atoms[0].sb_parity, original.sb_parity);
        assert_eq!(atoms[0].sn_ord, original.sn_ord);
        assert_eq!(atoms[0].sn_orig_at_num, original.sn_orig_at_num);
        assert_eq!(atoms[1].p_parity, original.p_parity);
        assert_eq!(atoms[1].p_orig_at_num, original.p_orig_at_num);
        assert_eq!(atoms[1].sb_ord, [-1, -2, -3]);
        assert_eq!(atoms[1].sb_parity, [1, 2, 3]);
        assert_eq!(atoms[1].sn_ord, [4, 5, 6]);
        assert_eq!(atoms[1].sn_orig_at_num, [61, 62, 63]);
        assert_eq!(atoms[2], original);
        assert_eq!(atoms[0].charge, -3);
        assert_eq!(atoms[1].charge, -3);
    }
}
