use crate::source::base::mol2atom::FreeExtOrigAtData;
use crate::source::base::util::inchi_calloc;
use crate::source_types::{
    OAD_Polymer, OAD_V3000, SourceHeap, SourceHeapError, SourceMutPointer, inchi_Input_Polymer,
    inchi_Input_V3000,
};

const SOURCE_SIZEOF_OAD_POLYMER: u64 = 48;
const SOURCE_SIZEOF_OAD_POLYMER_UNIT: u64 = 240;
const SOURCE_SIZEOF_OAD_V3000: u64 = 104;
const SOURCE_SIZEOF_POINTER: u64 = 8;
const SOURCE_SIZEOF_INT: u64 = 4;

#[allow(non_snake_case)]
pub(crate) fn SetExtOrigAtDataByInChIExtInput(
    heap: &mut SourceHeap,
    pp_polymer: &mut SourceMutPointer<OAD_Polymer>,
    pp_v3000: &mut SourceMutPointer<OAD_V3000>,
    iep: SourceMutPointer<inchi_Input_Polymer>,
    iev: SourceMutPointer<inchi_Input_V3000>,
    nat: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3021 SetExtOrigAtDataByInChIExtInput
    // INCHI✔❌: int SetExtOrigAtDataByInChIExtInput( OAD_Polymer **ppPolymer,
    // INCHI✔❌:                                      OAD_V3000 **ppV3000,
    // INCHI✔❌:                                      inchi_Input_Polymer *iep,
    // INCHI✔❌:                                      inchi_Input_V3000 *iev,
    // INCHI✔❌:                                      int nat )
    // INCHI✔❌: {
    // INCHI✔❌:     int    k, m, err = 0;
    // INCHI✔❌:     OAD_V3000 *pv = NULL;
    // INCHI✔❌:
    // INCHI✔❌:     /* Polymers */
    // INCHI✔❌:     if (iep && iep->n)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* Prepare OAD_Polymer container */
    // INCHI✔❌:         *ppPolymer = (OAD_Polymer *) inchi_calloc( 1, sizeof( OAD_Polymer ) );
    // INCHI✔❌:         if (!*ppPolymer)
    // INCHI✔❌:         {
    // INCHI✔❌:             err = 9001;
    // INCHI✔❌:             goto exitf;
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         /* Convert Molfile's Sgroup's to OAD_PolymerUnit's */
    // INCHI✔❌:         ( *ppPolymer )->units = (OAD_PolymerUnit**) inchi_calloc( iep->n, sizeof( ( *ppPolymer )->units[0] ) );
    // INCHI✔❌:         if (!( *ppPolymer )->units)
    // INCHI✔❌:         {
    // INCHI✔❌:             err = 9001;
    // INCHI✔❌:             goto exitf;
    // INCHI✔❌:         }
    // INCHI✔❌:         memset( ( *ppPolymer )->units, 0, sizeof( *( *ppPolymer )->units ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:
    // INCHI✔❌:         ( *ppPolymer )->n = iep->n;
    // INCHI✔❌:         /*( *ppPolymer )->valid = -1;*/
    // INCHI✔❌:         ( *ppPolymer )->really_do_frame_shift = 0;
    // INCHI✔❌:
    // INCHI✔❌:         for (k = 0; k < iep->n; k++)
    // INCHI✔❌:         {
    // INCHI✔❌:             int q = 0;
    // INCHI✔❌:             OAD_PolymerUnit *unitk;
    // INCHI✔❌:
    // INCHI✔❌:             inchi_Input_PolymerUnit *groupk = iep->units[k];
    // INCHI✔❌:             ( *ppPolymer )->units[k] = (OAD_PolymerUnit*) inchi_calloc( 1, sizeof( OAD_PolymerUnit ) );
    // INCHI✔❌:             unitk = ( *ppPolymer )->units[k];
    // INCHI✔❌:             if (!unitk)
    // INCHI✔❌:             {
    // INCHI✔❌:                 err = 9001;
    // INCHI✔❌:                 goto exitf;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             memset( unitk, 0, sizeof( *unitk ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:             unitk->id = groupk->id;
    // INCHI✔❌:             unitk->type = groupk->type;
    // INCHI✔❌:             unitk->subtype = groupk->subtype;
    // INCHI✔❌:             unitk->conn = groupk->conn;
    // INCHI✔❌:             unitk->label = groupk->label;
    // INCHI✔❌:
    // INCHI✔❌:             for (q = 0; q < 4; q++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 unitk->xbr1[q] = groupk->xbr1[q];
    // INCHI✔❌:                 unitk->xbr2[q] = groupk->xbr2[q];
    // INCHI✔❌:             }
    // INCHI✔❌:             strcpy( unitk->smt, groupk->smt );
    // INCHI✔❌:             unitk->na = groupk->na;
    // INCHI✔❌:             unitk->alist = (int *) inchi_calloc( unitk->na, sizeof( int ) );
    // INCHI✔❌:             if (!unitk->alist )
    // INCHI✔❌:             {
    // INCHI✔❌:                 err = 9001;
    // INCHI✔❌:                 goto exitf;
    // INCHI✔❌:             }
    // INCHI✔❌:             for (m = 0; m < unitk->na; m++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 unitk->alist[m] = groupk->alist[m];
    // INCHI✔❌:             }
    // INCHI✔❌:             unitk->nb = groupk->nb;
    // INCHI✔❌:             if (unitk->nb > 0)
    // INCHI✔❌:             {
    // INCHI✔❌:                 unitk->blist = (int *) inchi_calloc( 2 * (long long)unitk->nb, sizeof( int ) ); /* djb-rwth: cast operator added */
    // INCHI✔❌:                 if (!unitk->blist )
    // INCHI✔❌:                 {
    // INCHI✔❌:                     err = 9001;
    // INCHI✔❌:                     goto exitf;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 for (m = 0; m < 2 * groupk->nb; m++)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     unitk->blist[m] = groupk->blist[m];
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 unitk->blist = NULL;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* V3000 Extensions */
    // INCHI✔❌:     if (iev)
    // INCHI✔❌:     {
    // INCHI✔❌:         int nn;
    // INCHI✔❌:         *ppV3000 = (OAD_V3000 *) inchi_calloc( 1, sizeof( OAD_V3000 ) );
    // INCHI✔❌:         pv = *ppV3000;
    // INCHI✔❌:         if (!pv)
    // INCHI✔❌:         {
    // INCHI✔❌:             err = 9001;
    // INCHI✔❌:             goto exitf;
    // INCHI✔❌:         }
    // INCHI✔❌:         memset( pv, 0, sizeof( *pv ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:
    // INCHI✔❌:         pv->n_collections = iev->n_collections;
    // INCHI✔❌:         pv->n_haptic_bonds = iev->n_haptic_bonds;
    // INCHI✔❌:         pv->n_non_haptic_bonds = iev->n_non_haptic_bonds;
    // INCHI✔❌:         pv->n_sgroups = iev->n_sgroups;
    // INCHI✔❌:         pv->n_non_star_atoms = iev->n_non_star_atoms;
    // INCHI✔❌:         pv->n_star_atoms = iev->n_star_atoms;
    // INCHI✔❌:         pv->n_steabs = iev->n_steabs;
    // INCHI✔❌:         pv->n_sterac = iev->n_sterac;
    // INCHI✔❌:         pv->n_sterel = iev->n_sterel;
    // INCHI✔❌:         pv->n_3d_constraints = iev->n_3d_constraints;
    // INCHI✔❌:
    // INCHI✔❌:         if (iev->atom_index_orig)
    // INCHI✔❌:         {
    // INCHI✔❌:             pv->atom_index_orig = (int *) inchi_calloc( nat, sizeof( int ) );
    // INCHI✔❌:             if (NULL == pv->atom_index_orig)
    // INCHI✔❌:             {
    // INCHI✔❌:                 err = 9001;
    // INCHI✔❌:                 goto exitf;
    // INCHI✔❌:             }
    // INCHI✔❌:             memcpy( pv->atom_index_orig, iev->atom_index_orig, nat );
    // INCHI✔❌:         }
    // INCHI✔❌:         if (iev->atom_index_fin)
    // INCHI✔❌:         {
    // INCHI✔❌:             pv->atom_index_fin = (int *) inchi_calloc( nat, sizeof( int ) );
    // INCHI✔❌:             if (NULL == pv->atom_index_fin)
    // INCHI✔❌:             {
    // INCHI✔❌:                 err = 9001;
    // INCHI✔❌:                 goto exitf;
    // INCHI✔❌:             }
    // INCHI✔❌:             memcpy( pv->atom_index_fin, iev->atom_index_fin, nat );
    // INCHI✔❌:         }
    // INCHI✔❌:         if (iev->n_haptic_bonds && iev->lists_haptic_bonds)
    // INCHI✔❌:         {
    // INCHI✔❌:             pv->lists_haptic_bonds = (int **) inchi_calloc( iev->n_haptic_bonds, sizeof( int* ) );
    // INCHI✔❌:             if (NULL == pv->lists_haptic_bonds)
    // INCHI✔❌:             {
    // INCHI✔❌:                 err = 9001;
    // INCHI✔❌:                 goto exitf;
    // INCHI✔❌:             }
    // INCHI✔❌:             for (m = 0; m < iev->n_haptic_bonds; m++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 int *lst = NULL;
    // INCHI✔❌:                 int *mol_lst = iev->lists_haptic_bonds[m];
    // INCHI✔❌:                 nn = mol_lst[2] + 3;
    // INCHI✔❌:                 lst = pv->lists_haptic_bonds[m] = (int *) inchi_calloc( nn, sizeof( int ) );
    // INCHI✔❌:                 if (NULL == lst)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     err = 9001;
    // INCHI✔❌:                     goto exitf;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 for (k = 0; k < nn; k++)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     lst[k] = mol_lst[k];
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         if (iev->n_steabs && iev->lists_steabs)
    // INCHI✔❌:         {
    // INCHI✔❌:             pv->lists_steabs = (int **) inchi_calloc( iev->n_steabs, sizeof( int* ) );
    // INCHI✔❌:             if (NULL == pv->lists_steabs) { err = 9001; goto exitf; }
    // INCHI✔❌:             for (m = 0; m < iev->n_steabs; m++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 int *lst = NULL;
    // INCHI✔❌:                 int *mol_lst = iev->lists_steabs[m];
    // INCHI✔❌:                 nn = mol_lst[1] + 2;
    // INCHI✔❌:                 lst = pv->lists_steabs[m] = (int *) inchi_calloc( nn, sizeof( int ) );
    // INCHI✔❌:                 if (NULL == lst)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     err = 9001;
    // INCHI✔❌:                     goto exitf;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 for (k = 0; k < nn; k++)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     lst[k] = mol_lst[k];
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         if (iev->n_sterac && iev->lists_sterac)
    // INCHI✔❌:         {
    // INCHI✔❌:             pv->lists_sterac = (int **) inchi_calloc( iev->n_sterac, sizeof( int* ) );
    // INCHI✔❌:             if (NULL == pv->lists_sterac) { err = 9001; goto exitf; }
    // INCHI✔❌:             for (m = 0; m < iev->n_sterac; m++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 int *lst = NULL;
    // INCHI✔❌:                 int *mol_lst = iev->lists_sterac[m];
    // INCHI✔❌:                 nn = mol_lst[1] + 2;
    // INCHI✔❌:                 lst = pv->lists_sterac[m] = (int *) inchi_calloc( nn, sizeof( int ) );
    // INCHI✔❌:                 if (NULL == lst)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     err = 9001;
    // INCHI✔❌:                     goto exitf;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 for (k = 0; k < nn; k++)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     lst[k] = mol_lst[k];
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         if (iev->n_sterel && iev->lists_sterel)
    // INCHI✔❌:         {
    // INCHI✔❌:             pv->lists_sterel = (int **) inchi_calloc( iev->n_sterel, sizeof( int* ) );
    // INCHI✔❌:             if (NULL == pv->lists_sterel) { err = 9001; goto exitf; }
    // INCHI✔❌:             for (m = 0; m < iev->n_sterel; m++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 int *lst = NULL;
    // INCHI✔❌:                 int *mol_lst = iev->lists_sterel[m];
    // INCHI✔❌:                 nn = mol_lst[1] + 2;
    // INCHI✔❌:                 lst = pv->lists_sterel[m] = (int *) inchi_calloc( nn, sizeof( int ) );
    // INCHI✔❌:                 if (NULL == lst)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     err = 9001;
    // INCHI✔❌:                     goto exitf;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 for (k = 0; k < nn; k++)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     lst[k] = mol_lst[k];
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌: exitf:
    // INCHI✔❌:     if (err)
    // INCHI✔❌:     {
    // INCHI✔❌:         FreeExtOrigAtData( *ppPolymer, pv );
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return err;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: SetExtOrigAtDataByInChIExtInput

    let mut err = 0;
    let mut pv = SourceMutPointer::null();

    'port: {
        if !iep.is_null() {
            let input_polymer = heap
                .slice(iep.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if input_polymer.n != 0 {
                *pp_polymer = match inchi_calloc(heap, 1, SOURCE_SIZEOF_OAD_POLYMER) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        SourceMutPointer::null()
                    }
                };
                if pp_polymer.is_null() {
                    break 'port;
                }
                let unit_count = u64::try_from(input_polymer.n).unwrap_or(u64::MAX);
                let units = match inchi_calloc(heap, unit_count, SOURCE_SIZEOF_POINTER) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                {
                    let polymer = &mut heap.slice_mut(*pp_polymer)?[0];
                    polymer.units = units;
                    polymer.n = input_polymer.n;
                    polymer.really_do_frame_shift = 0;
                }
                for k in 0..input_polymer.n {
                    let group_pointer = heap.slice(input_polymer.units.as_const())?[k as usize];
                    let group = heap
                        .slice(group_pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    let unit = match inchi_calloc(heap, 1, SOURCE_SIZEOF_OAD_POLYMER_UNIT) {
                        Ok(pointer) => pointer,
                        Err(_) => {
                            err = 9001;
                            break 'port;
                        }
                    };
                    heap.slice_mut(units)?[k as usize] = unit;
                    {
                        let output = &mut heap.slice_mut(unit)?[0];
                        output.id = group.id;
                        output.type_ = group.type_;
                        output.subtype = group.subtype;
                        output.conn = group.conn;
                        output.label = group.label;
                        output.xbr1 = group.xbr1;
                        output.xbr2 = group.xbr2;
                        let nul = group
                            .smt
                            .iter()
                            .position(|byte| *byte == 0)
                            .ok_or(SourceHeapError::MissingNulTerminator)?;
                        output.smt[..=nul].copy_from_slice(&group.smt[..=nul]);
                        output.na = group.na;
                    }
                    let atom_count = u64::try_from(group.na).unwrap_or(u64::MAX);
                    let alist = match inchi_calloc(heap, atom_count, SOURCE_SIZEOF_INT) {
                        Ok(pointer) => pointer,
                        Err(_) => {
                            err = 9001;
                            break 'port;
                        }
                    };
                    heap.slice_mut(unit)?[0].alist = alist;
                    for m in 0..group.na {
                        heap.slice_mut(alist)?[m as usize] =
                            heap.slice(group.alist.as_const())?[m as usize];
                    }
                    heap.slice_mut(unit)?[0].nb = group.nb;
                    if group.nb > 0 {
                        let bond_count = u64::try_from(i64::from(group.nb) * 2).unwrap_or(u64::MAX);
                        let blist = match inchi_calloc(heap, bond_count, SOURCE_SIZEOF_INT) {
                            Ok(pointer) => pointer,
                            Err(_) => {
                                err = 9001;
                                break 'port;
                            }
                        };
                        heap.slice_mut(unit)?[0].blist = blist;
                        for m in 0..group.nb * 2 {
                            heap.slice_mut(blist)?[m as usize] =
                                heap.slice(group.blist.as_const())?[m as usize];
                        }
                    } else {
                        heap.slice_mut(unit)?[0].blist = SourceMutPointer::null();
                    }
                }
            }
        }

        if !iev.is_null() {
            let input = heap
                .slice(iev.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            *pp_v3000 = match inchi_calloc(heap, 1, SOURCE_SIZEOF_OAD_V3000) {
                Ok(pointer) => pointer,
                Err(_) => {
                    err = 9001;
                    SourceMutPointer::null()
                }
            };
            pv = *pp_v3000;
            if pv.is_null() {
                break 'port;
            }
            {
                let output = &mut heap.slice_mut(pv)?[0];
                output.n_collections = input.n_collections;
                output.n_haptic_bonds = input.n_haptic_bonds;
                output.n_non_haptic_bonds = input.n_non_haptic_bonds;
                output.n_sgroups = input.n_sgroups;
                output.n_non_star_atoms = input.n_non_star_atoms;
                output.n_star_atoms = input.n_star_atoms;
                output.n_steabs = input.n_steabs;
                output.n_sterac = input.n_sterac;
                output.n_sterel = input.n_sterel;
                output.n_3d_constraints = input.n_3d_constraints;
            }
            let nat_count = u64::try_from(nat).unwrap_or(u64::MAX);
            if !input.atom_index_orig.is_null() {
                let destination = match inchi_calloc(heap, nat_count, SOURCE_SIZEOF_INT) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                heap.slice_mut(pv)?[0].atom_index_orig = destination;
                for byte_index in 0..usize::try_from(nat).unwrap_or(0) {
                    let word = byte_index / 4;
                    let shift = (byte_index % 4) * 8;
                    let byte = ((heap.slice(input.atom_index_orig.as_const())?[word] as u32
                        >> shift)
                        & 0xff) as u32;
                    let target = &mut heap.slice_mut(destination)?[word];
                    *target = (((*target as u32) & !(0xff << shift)) | (byte << shift)) as i32;
                }
            }
            if !input.atom_index_fin.is_null() {
                let destination = match inchi_calloc(heap, nat_count, SOURCE_SIZEOF_INT) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                heap.slice_mut(pv)?[0].atom_index_fin = destination;
                for byte_index in 0..usize::try_from(nat).unwrap_or(0) {
                    let word = byte_index / 4;
                    let shift = (byte_index % 4) * 8;
                    let byte = ((heap.slice(input.atom_index_fin.as_const())?[word] as u32
                        >> shift)
                        & 0xff) as u32;
                    let target = &mut heap.slice_mut(destination)?[word];
                    *target = (((*target as u32) & !(0xff << shift)) | (byte << shift)) as i32;
                }
            }

            if input.n_haptic_bonds != 0 && !input.lists_haptic_bonds.is_null() {
                let count = u64::try_from(input.n_haptic_bonds).unwrap_or(u64::MAX);
                let outer = match inchi_calloc(heap, count, SOURCE_SIZEOF_POINTER) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                heap.slice_mut(pv)?[0].lists_haptic_bonds = outer;
                for m in 0..input.n_haptic_bonds {
                    let source = heap.slice(input.lists_haptic_bonds.as_const())?[m as usize];
                    let nn = heap.slice(source.as_const())?[2].wrapping_add(3);
                    let row = match inchi_calloc(
                        heap,
                        u64::try_from(nn).unwrap_or(u64::MAX),
                        SOURCE_SIZEOF_INT,
                    ) {
                        Ok(pointer) => pointer,
                        Err(_) => {
                            err = 9001;
                            break 'port;
                        }
                    };
                    heap.slice_mut(outer)?[m as usize] = row;
                    for k in 0..nn {
                        heap.slice_mut(row)?[k as usize] =
                            heap.slice(source.as_const())?[k as usize];
                    }
                }
            }
            if input.n_steabs != 0 && !input.lists_steabs.is_null() {
                let outer = match inchi_calloc(
                    heap,
                    u64::try_from(input.n_steabs).unwrap_or(u64::MAX),
                    SOURCE_SIZEOF_POINTER,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                heap.slice_mut(pv)?[0].lists_steabs = outer;
                for m in 0..input.n_steabs {
                    let source = heap.slice(input.lists_steabs.as_const())?[m as usize];
                    let nn = heap.slice(source.as_const())?[1].wrapping_add(2);
                    let row = match inchi_calloc(
                        heap,
                        u64::try_from(nn).unwrap_or(u64::MAX),
                        SOURCE_SIZEOF_INT,
                    ) {
                        Ok(pointer) => pointer,
                        Err(_) => {
                            err = 9001;
                            break 'port;
                        }
                    };
                    heap.slice_mut(outer)?[m as usize] = row;
                    for k in 0..nn {
                        heap.slice_mut(row)?[k as usize] =
                            heap.slice(source.as_const())?[k as usize];
                    }
                }
            }
            if input.n_sterac != 0 && !input.lists_sterac.is_null() {
                let outer = match inchi_calloc(
                    heap,
                    u64::try_from(input.n_sterac).unwrap_or(u64::MAX),
                    SOURCE_SIZEOF_POINTER,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                heap.slice_mut(pv)?[0].lists_sterac = outer;
                for m in 0..input.n_sterac {
                    let source = heap.slice(input.lists_sterac.as_const())?[m as usize];
                    let nn = heap.slice(source.as_const())?[1].wrapping_add(2);
                    let row = match inchi_calloc(
                        heap,
                        u64::try_from(nn).unwrap_or(u64::MAX),
                        SOURCE_SIZEOF_INT,
                    ) {
                        Ok(pointer) => pointer,
                        Err(_) => {
                            err = 9001;
                            break 'port;
                        }
                    };
                    heap.slice_mut(outer)?[m as usize] = row;
                    for k in 0..nn {
                        heap.slice_mut(row)?[k as usize] =
                            heap.slice(source.as_const())?[k as usize];
                    }
                }
            }
            if input.n_sterel != 0 && !input.lists_sterel.is_null() {
                let outer = match inchi_calloc(
                    heap,
                    u64::try_from(input.n_sterel).unwrap_or(u64::MAX),
                    SOURCE_SIZEOF_POINTER,
                ) {
                    Ok(pointer) => pointer,
                    Err(_) => {
                        err = 9001;
                        break 'port;
                    }
                };
                heap.slice_mut(pv)?[0].lists_sterel = outer;
                for m in 0..input.n_sterel {
                    let source = heap.slice(input.lists_sterel.as_const())?[m as usize];
                    let nn = heap.slice(source.as_const())?[1].wrapping_add(2);
                    let row = match inchi_calloc(
                        heap,
                        u64::try_from(nn).unwrap_or(u64::MAX),
                        SOURCE_SIZEOF_INT,
                    ) {
                        Ok(pointer) => pointer,
                        Err(_) => {
                            err = 9001;
                            break 'port;
                        }
                    };
                    heap.slice_mut(outer)?[m as usize] = row;
                    for k in 0..nn {
                        heap.slice_mut(row)?[k as usize] =
                            heap.slice(source.as_const())?[k as usize];
                    }
                }
            }
        }
    }

    if err != 0 {
        FreeExtOrigAtData(heap, *pp_polymer, pv)?;
    }
    Ok(err)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::inchi_Input_PolymerUnit;

    struct Fixture {
        polymer: SourceMutPointer<inchi_Input_Polymer>,
        v3000: SourceMutPointer<inchi_Input_V3000>,
    }

    fn allocate_fixture(heap: &mut SourceHeap, bonds: i32) -> Fixture {
        let alist = heap.allocate(vec![3_i32, 4]).unwrap();
        let blist = if bonds > 0 {
            heap.allocate(vec![3_i32, 4]).unwrap()
        } else {
            SourceMutPointer::null()
        };
        let mut smt = [0_i8; 80];
        smt[..7].copy_from_slice(&[b'A' as i8, b'B' as i8, 0, 9, 8, 7, 6]);
        let unit = heap
            .allocate(vec![inchi_Input_PolymerUnit {
                id: 11,
                type_: 12,
                subtype: 13,
                conn: 14,
                label: 15,
                na: 2,
                nb: bonds,
                xbr1: [1.0, 2.0, 3.0, 4.0],
                xbr2: [5.0, 6.0, 7.0, 8.0],
                smt,
                alist,
                blist,
            }])
            .unwrap();
        let units = heap.allocate(vec![unit]).unwrap();
        let polymer = heap
            .allocate(vec![inchi_Input_Polymer { units, n: 1 }])
            .unwrap();

        let atom_index_orig = heap
            .allocate(vec![0x1122_3344_i32, 0x5566_7788_i32])
            .unwrap();
        let atom_index_fin = heap
            .allocate(vec![0x1020_3040_i32, 0x5060_7080_i32])
            .unwrap();
        let haptic_row = heap.allocate(vec![7_i32, 8, 2, 11, 12]).unwrap();
        let steabs_row = heap.allocate(vec![0_i32, 2, 21, 22]).unwrap();
        let sterac_row = heap.allocate(vec![31_i32, 2, 32, 33]).unwrap();
        let sterel_row = heap.allocate(vec![41_i32, 2, 42, 43]).unwrap();
        let lists_haptic_bonds = heap.allocate(vec![haptic_row]).unwrap();
        let lists_steabs = heap.allocate(vec![steabs_row]).unwrap();
        let lists_sterac = heap.allocate(vec![sterac_row]).unwrap();
        let lists_sterel = heap.allocate(vec![sterel_row]).unwrap();
        let v3000 = heap
            .allocate(vec![inchi_Input_V3000 {
                n_non_star_atoms: 51,
                n_star_atoms: 52,
                atom_index_orig,
                atom_index_fin,
                n_sgroups: 53,
                n_3d_constraints: 54,
                n_collections: 55,
                n_non_haptic_bonds: 56,
                n_haptic_bonds: 1,
                lists_haptic_bonds,
                n_steabs: 1,
                lists_steabs,
                n_sterel: 1,
                lists_sterel,
                n_sterac: 1,
                lists_sterac,
            }])
            .unwrap();
        Fixture { polymer, v3000 }
    }

    #[test]
    fn source_port__inchi_dll__setextorigatdatabyinchiextinput__line_3021() {
        let mut heap = SourceHeap::default();
        let old_polymer = heap.allocate(vec![OAD_Polymer::default()]).unwrap();
        let old_v3000 = heap.allocate(vec![OAD_V3000::default()]).unwrap();
        let mut polymer_slot = old_polymer;
        let mut v3000_slot = old_v3000;
        assert_eq!(
            SetExtOrigAtDataByInChIExtInput(
                &mut heap,
                &mut polymer_slot,
                &mut v3000_slot,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                i32::MIN,
            ),
            Ok(0)
        );
        assert_eq!(polymer_slot, old_polymer);
        assert_eq!(v3000_slot, old_v3000);
        FreeExtOrigAtData(&mut heap, old_polymer, old_v3000).unwrap();

        let zero_input = heap
            .allocate(vec![inchi_Input_Polymer {
                units: SourceMutPointer::null(),
                n: 0,
            }])
            .unwrap();
        let untouched = heap.allocate(vec![OAD_Polymer::default()]).unwrap();
        let mut polymer_slot = untouched;
        let mut v3000_slot = SourceMutPointer::null();
        assert_eq!(
            SetExtOrigAtDataByInChIExtInput(
                &mut heap,
                &mut polymer_slot,
                &mut v3000_slot,
                zero_input,
                SourceMutPointer::null(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(polymer_slot, untouched);
        FreeExtOrigAtData(&mut heap, untouched, SourceMutPointer::null()).unwrap();

        let fixture = allocate_fixture(&mut heap, 1);
        let mut polymer_slot = SourceMutPointer::null();
        let mut v3000_slot = SourceMutPointer::null();
        assert_eq!(
            SetExtOrigAtDataByInChIExtInput(
                &mut heap,
                &mut polymer_slot,
                &mut v3000_slot,
                fixture.polymer,
                fixture.v3000,
                5,
            ),
            Ok(0)
        );
        let polymer = &heap.slice(polymer_slot.as_const()).unwrap()[0];
        assert_eq!(polymer.n, 1);
        assert_eq!(polymer.really_do_frame_shift, 0);
        let unit_pointer = heap.slice(polymer.units.as_const()).unwrap()[0];
        let unit = &heap.slice(unit_pointer.as_const()).unwrap()[0];
        assert_eq!(
            (unit.id, unit.type_, unit.subtype, unit.conn, unit.label),
            (11, 12, 13, 14, 15)
        );
        assert_eq!(unit.xbr1, [1.0, 2.0, 3.0, 4.0]);
        assert_eq!(unit.xbr2, [5.0, 6.0, 7.0, 8.0]);
        assert_eq!(&unit.smt[..7], &[b'A' as i8, b'B' as i8, 0, 0, 0, 0, 0]);
        assert_eq!(heap.slice(unit.alist.as_const()).unwrap(), &[3, 4]);
        assert_eq!(heap.slice(unit.blist.as_const()).unwrap(), &[3, 4]);

        let output = &heap.slice(v3000_slot.as_const()).unwrap()[0];
        assert_eq!(
            (
                output.n_non_star_atoms,
                output.n_star_atoms,
                output.n_sgroups,
                output.n_3d_constraints,
                output.n_collections,
                output.n_non_haptic_bonds,
            ),
            (51, 52, 53, 54, 55, 56)
        );
        assert_eq!(
            heap.slice(output.atom_index_orig.as_const()).unwrap(),
            &[0x1122_3344, 0x88, 0, 0, 0]
        );
        assert_eq!(
            heap.slice(output.atom_index_fin.as_const()).unwrap(),
            &[0x1020_3040, 0x80, 0, 0, 0]
        );
        let haptic = heap.slice(output.lists_haptic_bonds.as_const()).unwrap()[0];
        let steabs = heap.slice(output.lists_steabs.as_const()).unwrap()[0];
        let sterac = heap.slice(output.lists_sterac.as_const()).unwrap()[0];
        let sterel = heap.slice(output.lists_sterel.as_const()).unwrap()[0];
        assert_eq!(heap.slice(haptic.as_const()).unwrap(), &[7, 8, 2, 11, 12]);
        assert_eq!(heap.slice(steabs.as_const()).unwrap(), &[0, 2, 21, 22]);
        assert_eq!(heap.slice(sterac.as_const()).unwrap(), &[31, 2, 32, 33]);
        assert_eq!(heap.slice(sterel.as_const()).unwrap(), &[41, 2, 42, 43]);
        FreeExtOrigAtData(&mut heap, polymer_slot, v3000_slot).unwrap();

        let no_bond_fixture = allocate_fixture(&mut heap, 0);
        let mut polymer_slot = SourceMutPointer::null();
        let mut v3000_slot = SourceMutPointer::null();
        assert_eq!(
            SetExtOrigAtDataByInChIExtInput(
                &mut heap,
                &mut polymer_slot,
                &mut v3000_slot,
                no_bond_fixture.polymer,
                SourceMutPointer::null(),
                0,
            ),
            Ok(0)
        );
        let polymer = &heap.slice(polymer_slot.as_const()).unwrap()[0];
        let unit_pointer = heap.slice(polymer.units.as_const()).unwrap()[0];
        assert!(
            heap.slice(unit_pointer.as_const()).unwrap()[0]
                .blist
                .is_null()
        );
        FreeExtOrigAtData(&mut heap, polymer_slot, SourceMutPointer::null()).unwrap();

        for successful_allocations in 0..16_u64 {
            let mut heap = SourceHeap::default();
            let fixture = allocate_fixture(&mut heap, 1);
            heap.fail_after_allocations(successful_allocations);
            let mut polymer_slot = SourceMutPointer::null();
            let mut v3000_slot = SourceMutPointer::null();
            assert_eq!(
                SetExtOrigAtDataByInChIExtInput(
                    &mut heap,
                    &mut polymer_slot,
                    &mut v3000_slot,
                    fixture.polymer,
                    fixture.v3000,
                    5,
                ),
                Ok(9001),
                "allocation ordinal {successful_allocations}"
            );
            assert_eq!(
                heap.source_allocation_calls(),
                successful_allocations + 1,
                "allocation ordinal {successful_allocations}"
            );
            if !polymer_slot.is_null() {
                assert_eq!(
                    heap.slice(polymer_slot.as_const()),
                    Err(SourceHeapError::MissingAllocation),
                    "allocation ordinal {successful_allocations}"
                );
            }
            if !v3000_slot.is_null() {
                assert_eq!(
                    heap.slice(v3000_slot.as_const()),
                    Err(SourceHeapError::MissingAllocation),
                    "allocation ordinal {successful_allocations}"
                );
            }
        }
    }
}
