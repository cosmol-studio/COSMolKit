use crate::source::base::ichinorm::MarkRingSystemsInp;
use crate::source::base::ichisort::iisort;
use crate::source::base::strutil::{imat_free, imat_new};
use crate::source::base::util::{inchi_calloc, inchi_free, is_in_the_ilist};
use crate::source_types::{
    CLOSING_SRU_DIRADICAL, CLOSING_SRU_HIGHER_ORDER_BOND, CLOSING_SRU_RING, OAD_AtProps,
    OAD_Polymer, OAD_PolymerUnit, SourceConstPointer, SourceHeap, SourceHeapError,
    SourceMutPointer, inp_ATOM, tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE,
    tagINCHIBondType_INCHI_BOND_TYPE_SINGLE, tagINCHIBondType_INCHI_BOND_TYPE_TRIPLE,
};

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_Free(
    heap: &mut SourceHeap,
    unit: SourceMutPointer<OAD_PolymerUnit>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1342 OAD_PolymerUnit_Free
    // INCHI✔❌: void OAD_PolymerUnit_Free( OAD_PolymerUnit *unit )
    // INCHI✔❌: {
    // INCHI✔❌:
    // INCHI✔❌:     ITRACE_( "\n************** About to free OAD_PolymerUnit @ %-p\n", unit );
    // INCHI✔❌:     OAD_PolymerUnit_DebugTrace( unit );
    // INCHI✔❌:
    // INCHI✔❌:     if (unit)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (unit->alist)
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_free( unit->alist );
    // INCHI✔❌:             unit->alist = NULL;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (unit->blist)
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_free( unit->blist );
    // INCHI✔❌:             unit->blist = NULL;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (unit->bkbonds)
    // INCHI✔❌:         {
    // INCHI✔❌:             imat_free( unit->maxbkbonds, unit->bkbonds );
    // INCHI✔❌:             unit->bkbonds = NULL;
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     inchi_free( unit );
    // INCHI✔❌:
    // INCHI✔❌:     return;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_Free
    // BEGIN INCHI ACTIVE HEADER MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.h:133 ITRACE_
    // INCHI✔❌: #define ITRACE_ 0 && _inchi_trace
    // END INCHI ACTIVE HEADER MACRO: ITRACE_

    if unit.is_null() {
        OAD_PolymerUnit_DebugTrace(None);
    } else {
        let unit_ref = heap
            .slice(unit.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        OAD_PolymerUnit_DebugTrace(Some(unit_ref));
    }

    if !unit.is_null() {
        let (alist, blist, maxbkbonds, bkbonds) = {
            let unit_ref = heap
                .slice(unit.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (
                unit_ref.alist,
                unit_ref.blist,
                unit_ref.maxbkbonds,
                unit_ref.bkbonds,
            )
        };
        if !alist.is_null() {
            inchi_free(heap, alist)?;
            heap.slice_mut(unit)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .alist = SourceMutPointer::null();
        }
        if !blist.is_null() {
            inchi_free(heap, blist)?;
            heap.slice_mut(unit)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .blist = SourceMutPointer::null();
        }
        if !bkbonds.is_null() {
            imat_free(heap, maxbkbonds, bkbonds)?;
            heap.slice_mut(unit)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .bkbonds = SourceMutPointer::null();
        }
    }

    inchi_free(heap, unit)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_CreateCopy(
    heap: &mut SourceHeap,
    unit: SourceMutPointer<OAD_PolymerUnit>,
) -> Result<SourceMutPointer<OAD_PolymerUnit>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1260 OAD_PolymerUnit_CreateCopy
    // INCHI✔️❌: OAD_PolymerUnit* OAD_PolymerUnit_CreateCopy( OAD_PolymerUnit *u )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int k, err = 0;
    // INCHI✔️❌:     OAD_PolymerUnit *u2 = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     u2 = (OAD_PolymerUnit*) inchi_calloc( 1, sizeof( OAD_PolymerUnit ) );
    // INCHI✔️❌:     if (NULL == u2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         err = 1;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     u2->id = u->id;
    // INCHI✔️❌:     u2->type = u->type;
    // INCHI✔️❌:     u2->subtype = u->subtype;
    // INCHI✔️❌:     u2->conn = u->conn;
    // INCHI✔️❌:     u2->label = u->label;
    // INCHI✔️❌:     u2->na = u->na;
    // INCHI✔️❌:     u2->nb = u->nb;
    // INCHI✔️❌:     u2->cyclizable = u->cyclizable;
    // INCHI✔️❌:     u2->cyclized = u->cyclized;
    // INCHI✔️❌:     u2->cap1_is_undef = u->cap1_is_undef;
    // INCHI✔️❌:     u2->cap2_is_undef = u->cap2_is_undef;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (k = 0; k < 4; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u2->xbr1[k] = u->xbr1[k];
    // INCHI✔️❌:         u2->xbr2[k] = u->xbr2[k];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     strcpy(u2->smt, u->smt);
    // INCHI✔️❌:
    // INCHI✔️❌:     u2->cap1 = u->cap1;
    // INCHI✔️❌:     u2->end_atom1 = u->end_atom1;
    // INCHI✔️❌:     u2->cap2 = u->cap2;
    // INCHI✔️❌:     u2->end_atom2 = u->end_atom2;
    // INCHI✔️❌:     u2->nbkbonds = u->nbkbonds;
    // INCHI✔️❌:     u2->maxbkbonds = inchi_max( u->maxbkbonds, u->nbkbonds );
    // INCHI✔️❌:
    // INCHI✔️❌:     u2->alist = (int *) inchi_calloc( u2->na, sizeof( int ) );
    // INCHI✔️❌:     if (!u2->alist)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         err = 2;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (k = 0; k < u2->na; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u2->alist[k] = u->alist[k];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     u2->blist = (int *) inchi_calloc( 2 * (long long)u2->nb, sizeof( int ) );
    // INCHI✔️❌:     if (!u2->blist)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         err = 2;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (k = 0; k < 2 * u2->nb; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u2->blist[k] = u->blist[k];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     err = imat_new( u2->maxbkbonds, 2, &( u2->bkbonds ) );
    // INCHI✔️❌:     if (!err)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (k = 0; k < u2->nbkbonds; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             u2->bkbonds[k][0] = u->bkbonds[k][0];
    // INCHI✔️❌:             u2->bkbonds[k][1] = u->bkbonds[k][1];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     if (err)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         OAD_PolymerUnit_Free( u2 );
    // INCHI✔️❌:         return NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return u2;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_CreateCopy

    let source = heap
        .slice(unit.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let copy = match inchi_calloc::<OAD_PolymerUnit>(
        heap,
        1,
        std::mem::size_of::<OAD_PolymerUnit>() as u64,
    ) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(SourceMutPointer::null()),
        Err(error) => return Err(error),
    };
    let mut copied = source.clone();
    copied.alist = SourceMutPointer::null();
    copied.blist = SourceMutPointer::null();
    copied.bkbonds = SourceMutPointer::null();
    copied.representation = 0;
    copied.maxbkbonds = source.maxbkbonds.max(source.nbkbonds);
    copied.smt = [0; 80];
    let smt_length = source
        .smt
        .iter()
        .position(|byte| *byte == 0)
        .unwrap_or(source.smt.len());
    copied.smt[..smt_length].copy_from_slice(&source.smt[..smt_length]);
    if smt_length < copied.smt.len() {
        copied.smt[smt_length] = 0;
    }
    heap.slice_mut(copy)?[0] = copied;

    let fail_copy = |heap: &mut SourceHeap,
                     copy: SourceMutPointer<OAD_PolymerUnit>|
     -> Result<SourceMutPointer<OAD_PolymerUnit>, SourceHeapError> {
        OAD_PolymerUnit_Free(heap, copy)?;
        Ok(SourceMutPointer::null())
    };

    let alist_count = match u64::try_from(source.na) {
        Ok(value) => value,
        Err(_) => return fail_copy(heap, copy),
    };
    let alist = match inchi_calloc::<i32>(heap, alist_count, 4) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return fail_copy(heap, copy),
        Err(error) => return Err(error),
    };
    if !source.alist.is_null() {
        let source_values = heap.slice(source.alist.as_const())?.to_vec();
        let target = heap.slice_mut(alist)?;
        for (index, value) in target.iter_mut().enumerate() {
            *value = *source_values
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
        }
    }
    heap.slice_mut(copy)?[0].alist = alist;

    let blist_count = i64::from(source.nb)
        .checked_mul(2)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let blist_count = match u64::try_from(blist_count) {
        Ok(value) => value,
        Err(_) => return fail_copy(heap, copy),
    };
    let blist = match inchi_calloc::<i32>(heap, blist_count, 4) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return fail_copy(heap, copy),
        Err(error) => return Err(error),
    };
    if !source.blist.is_null() {
        let source_values = heap.slice(source.blist.as_const())?.to_vec();
        let target = heap.slice_mut(blist)?;
        for (index, value) in target.iter_mut().enumerate() {
            *value = *source_values
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
        }
    }
    heap.slice_mut(copy)?[0].blist = blist;

    let max_bonds = heap.slice(copy.as_const())?[0].maxbkbonds;
    let mut backbone_bonds = heap.slice(copy.as_const())?[0].bkbonds;
    let bond_error = imat_new(heap, max_bonds, 2, &mut backbone_bonds)?;
    heap.slice_mut(copy)?[0].bkbonds = backbone_bonds;
    if bond_error != 0 {
        return fail_copy(heap, copy);
    }
    if !source.bkbonds.is_null() {
        for index in 0..source.nbkbonds {
            let source_row = heap
                .slice(source.bkbonds.as_const())?
                .get(index as usize)
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let target_rows = heap.slice(copy.as_const())?[0].bkbonds;
            let target_row = heap
                .slice(target_rows.as_const())?
                .get(index as usize)
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let values = heap
                .slice(source_row.as_const())?
                .get(..2)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let values = values.to_vec();
            heap.slice_mut(target_row)?
                .get_mut(..2)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .copy_from_slice(&values);
        }
    }
    Ok(copy)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_CompareAtomListsMod(
    heap: &SourceHeap,
    first: SourceConstPointer<OAD_PolymerUnit>,
    second: SourceConstPointer<OAD_PolymerUnit>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1377 OAD_PolymerUnit_CompareAtomListsMod
    // INCHI❌❌: int  OAD_PolymerUnit_CompareAtomListsMod( OAD_PolymerUnit *u1,
    // INCHI❌❌:                                           OAD_PolymerUnit *u2 )
    // INCHI❌❌: {
    // INCHI❌❌:     int i;
    // INCHI❌❌:     int n1 = u1->na;
    // INCHI❌❌:     int n2 = u2->na;
    // INCHI❌❌:     int n = n1;
    // INCHI❌❌:     if (n1 < n2)    return -1;
    // INCHI❌❌:     if (n1 > n2)    return 1;
    // INCHI❌❌:     /* n1 == n2 == n */
    // INCHI❌❌:     for (i = 0; i < n; i++)
    // INCHI❌❌:     {
    // INCHI❌❌:         if (u1->alist[i] < u2->alist[i])    return -1;
    // INCHI❌❌:         if (u1->alist[i] > u2->alist[i])    return    1;
    // INCHI❌❌:     }
    // INCHI❌❌:
    // INCHI❌❌:     return 0;
    // INCHI❌❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_CompareAtomListsMod

    let first = heap
        .slice(first)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second = heap
        .slice(second)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if first.na < second.na {
        return Ok(-1);
    }
    if first.na > second.na {
        return Ok(1);
    }
    if first.na <= 0 {
        return Ok(0);
    }
    let count = usize::try_from(first.na).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let first_values = heap.slice(first.alist.as_const())?;
    let second_values = heap.slice(second.alist.as_const())?;
    for index in 0..count {
        if first_values[index] < second_values[index] {
            return Ok(-1);
        }
        if first_values[index] > second_values[index] {
            return Ok(1);
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves(
    heap: &mut SourceHeap,
    unit: &mut OAD_PolymerUnit,
    number_of_star_atoms: i32,
    star_atoms: SourceConstPointer<i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1437 OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves
    // INCHI❌❌: int  OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves( OAD_PolymerUnit  *u,
    // INCHI❌❌:                                                        int n_star_atoms,
    // INCHI❌❌:                                                        int *star_atoms )
    // INCHI❌❌: {
    // INCHI❌❌:     int k;
    // INCHI❌❌:
    // INCHI❌❌:     /* Sort bond atoms */
    // INCHI❌❌:     for (k = 0; k < u->nb; k++)
    // INCHI❌❌:     {
    // INCHI❌❌:         /* Place not-in-unit bond end to first place */
    // INCHI❌❌:         int a1 = u->blist[2 * k];
    // INCHI❌❌:         int a2 = u->blist[2 * k + 1];
    // INCHI❌❌:         int a1_is_not_in_alist = 0;
    // INCHI❌❌:         int a1_is_star_atom = 0;
    // INCHI❌❌:         int a2_is_not_in_alist = 0;
    // INCHI❌❌:         int a2_is_star_atom = 0;
    // INCHI❌❌:
    // INCHI❌❌:         if (!is_in_the_ilist( u->alist, a1, u->na ))
    // INCHI❌❌:         {
    // INCHI❌❌:             a1_is_not_in_alist = 1;
    // INCHI❌❌:         }
    // INCHI❌❌:         if (is_in_the_ilist( star_atoms, a1, n_star_atoms ))
    // INCHI❌❌:         {
    // INCHI❌❌:             a1_is_star_atom = 1;
    // INCHI❌❌:         }
    // INCHI❌❌:
    // INCHI❌❌:         if (!is_in_the_ilist( u->alist, a2, u->na ))
    // INCHI❌❌:         {
    // INCHI❌❌:             a2_is_not_in_alist = 1;
    // INCHI❌❌:         }
    // INCHI❌❌:         if (is_in_the_ilist( star_atoms, a2, n_star_atoms ))
    // INCHI❌❌:         {
    // INCHI❌❌:             a2_is_star_atom = 1;
    // INCHI❌❌:         }
    // INCHI❌❌:
    // INCHI❌❌:         if (( a1_is_not_in_alist || a1_is_star_atom ) &&
    // INCHI❌❌:             ( a2_is_not_in_alist || a2_is_star_atom ))
    // INCHI❌❌:         {
    // INCHI❌❌:             /* Both the ends are out of unit: the crossing bond is invalid */
    // INCHI❌❌:             return 1;
    // INCHI❌❌:         }
    // INCHI❌❌:         /* If a2 is star atom or non-star external to the current unit, swap(a2,a1) */
    // INCHI❌❌:         if (a2_is_star_atom || a2_is_not_in_alist)
    // INCHI❌❌:         {
    // INCHI❌❌:             u->blist[2 * k] = a2;
    // INCHI❌❌:             u->blist[2 * k + 1] = a1;
    // INCHI❌❌:         }
    // INCHI❌❌:     }
    // INCHI❌❌:
    // INCHI❌❌:     /* Sort bond themselves
    // INCHI❌❌:         for now, consider only the simplest cases of 2 bonds
    // INCHI❌❌:     */
    // INCHI❌❌:     if (u->nb == 2)            /* two bonds in SBL */
    // INCHI❌❌:     {
    // INCHI❌❌:         int b1a1 = u->blist[0];
    // INCHI❌❌:         int b1a2 = u->blist[1];
    // INCHI❌❌:         int b2a1 = u->blist[2];
    // INCHI❌❌:         int b2a2 = u->blist[3];
    // INCHI❌❌:         if (b1a1 > b2a1)
    // INCHI❌❌:         {
    // INCHI❌❌:             /* swap */
    // INCHI❌❌:             u->blist[0] = b2a1; u->blist[1] = b2a2;
    // INCHI❌❌:             u->blist[2] = b1a1; u->blist[3] = b1a2;
    // INCHI❌❌:         }
    // INCHI❌❌:     }
    // INCHI❌❌:
    // INCHI❌❌:     /* for single or no bonds, do nothing
    // INCHI❌❌:     else
    // INCHI❌❌:         ;
    // INCHI❌❌:     */
    // INCHI❌❌:
    // INCHI❌❌:     return 0;
    // INCHI❌❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves

    if unit.nb <= 0 {
        return Ok(0);
    }
    let alist = if unit.na == 0 {
        None
    } else {
        Some(heap.slice(unit.alist.as_const())?.to_vec())
    };
    let stars = if number_of_star_atoms == 0 {
        None
    } else {
        Some(heap.slice(star_atoms)?.to_vec())
    };
    for bond_index in 0..unit.nb {
        let offset = usize::try_from(i64::from(bond_index) * 2)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let (first, second) = {
            let bonds = heap.slice(unit.blist.as_const())?;
            (
                *bonds
                    .get(offset)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                *bonds
                    .get(offset + 1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        };
        let first_external = is_in_the_ilist(alist.as_deref(), first, unit.na)?.is_none();
        let first_star = is_in_the_ilist(stars.as_deref(), first, number_of_star_atoms)?.is_some();
        let second_external = is_in_the_ilist(alist.as_deref(), second, unit.na)?.is_none();
        let second_star =
            is_in_the_ilist(stars.as_deref(), second, number_of_star_atoms)?.is_some();
        if (first_external || first_star) && (second_external || second_star) {
            return Ok(1);
        }
        if second_star || second_external {
            let bonds = heap.slice_mut(unit.blist)?;
            bonds[offset] = second;
            bonds[offset + 1] = first;
        }
    }
    if unit.nb == 2 {
        let bonds = heap.slice_mut(unit.blist)?;
        let [first_first, first_second, second_first, second_second] =
            [bonds[0], bonds[1], bonds[2], bonds[3]];
        if first_first > second_first {
            bonds[0] = second_first;
            bonds[1] = second_second;
            bonds[2] = first_first;
            bonds[3] = first_second;
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_PrepareWorkingSet(
    heap: &mut SourceHeap,
    polymer: &mut OAD_Polymer,
    canonical_numbers: SourceConstPointer<i32>,
    _component_numbers: SourceConstPointer<i32>,
    units2: SourceMutPointer<SourceMutPointer<OAD_PolymerUnit>>,
    unit_numbers: SourceMutPointer<i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2342 OAD_Polymer_PrepareWorkingSet
    // INCHI✔️❌: int OAD_Polymer_PrepareWorkingSet( OAD_Polymer     *p,
    // INCHI✔️❌:                                    int             *cano_nums,
    // INCHI✔️❌:                                    int             *compnt_nums,
    // INCHI✔️❌:                                    OAD_PolymerUnit **units2,
    // INCHI✔️❌:                                    int             *unum )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, k, err = 0, cano_num1 = -1, cano_num2 = -1;
    // INCHI✔️❌:     OAD_PolymerUnit *u;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*OAD_Polymer_DebugTrace( p );*/
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  Replace original atom numbers in polymer data with canonical plus 1.                    */
    // INCHI✔️❌:     /*  Note that we use 'cano1 nums', that is, 1-based (InChI internal 'cano nums' are 0-based)*/
    // INCHI✔️❌:     /*  Also remove from the list atoms who mapped to cano number 0  ( == -1 + 1_offset ),      */
    // INCHI✔️❌:     /*  they are explicit H's which have already been deleted.                                  */
    // INCHI✔️❌:     for (k = 0; k < p->n_pzz; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         cano_num1 = cano_nums[p->pzz[k]] + 1;
    // INCHI✔️❌:         if (cano_num1 == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* we shouldn't arrive here */
    // INCHI✔️❌:             err = 10;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         p->pzz[k] = cano_num1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < p->n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int na_new = -1;
    // INCHI✔️❌:         u = units2[i];
    // INCHI✔️❌:
    // INCHI✔️❌:         for (k = 0; k < u->na; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cano_num1 = cano_nums[u->alist[k]] + 1;
    // INCHI✔️❌:             if (cano_num1 == 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             u->alist[++na_new] = cano_num1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         u->na = na_new + 1;
    // INCHI✔️❌:         for (k = 0; k < 2 * u->nb; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cano_num1 = cano_nums[u->blist[k]] + 1;
    // INCHI✔️❌:             if (cano_num1 == 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* Can not proceed further as one of PU crossing bond ends
    // INCHI✔️❌:                    leads to explicit H which has been removed already       */
    // INCHI✔️❌:                 err = 11;
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             u->blist[k] = cano_num1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         cano_num1 = cano_nums[u->cap1] + 1;
    // INCHI✔️❌:         if (cano_num1 == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             err = 11;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         u->cap1 = cano_num1;
    // INCHI✔️❌:
    // INCHI✔️❌:         cano_num1 = cano_nums[u->cap2] + 1;
    // INCHI✔️❌:         if (cano_num1 == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             err = 11;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         u->cap2 = cano_num1;
    // INCHI✔️❌:
    // INCHI✔️❌:         cano_num1 = cano_nums[u->end_atom1] + 1;
    // INCHI✔️❌:         if (cano_num1 == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             err = 11;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         u->end_atom1 = cano_num1;
    // INCHI✔️❌:
    // INCHI✔️❌:         cano_num1 = cano_nums[u->end_atom2] + 1;
    // INCHI✔️❌:         if (cano_num1 == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             err = 11;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         u->end_atom2 = cano_num1;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (k = 0; k < u->nbkbonds; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cano_num1 = cano_nums[u->bkbonds[k][0]] + 1;
    // INCHI✔️❌:             if (cano_num1 == 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             cano_num2 = cano_nums[u->bkbonds[k][1]] + 1;
    // INCHI✔️❌:             if (cano_num2 == 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             u->bkbonds[k][0] = inchi_min( cano_num1, cano_num2 );
    // INCHI✔️❌:             u->bkbonds[k][1] = inchi_max( cano_num1, cano_num2 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Sort the atoms and the bonds in all units */
    // INCHI✔️❌:     for (i = 0; i < p->n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u = units2[i];
    // INCHI✔️❌:
    // INCHI✔️❌:         /* sort atoms (alist) */
    // INCHI✔️❌:         iisort( u->alist, u->na );
    // INCHI✔️❌:
    // INCHI✔️❌:         /*ITRACE_( "\n*** Polymer unit %-d : ( ", i );
    // INCHI✔️❌:         for (k = 0; k < u->na - 1; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ITRACE_( "%-d-", u->alist[k] );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         ITRACE_( "%-d )\n", u->alist[u->na - 1] );*/
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Sort bonds (blist) */
    // INCHI✔️❌:         err = OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves( u, p->n_pzz, p->pzz );
    // INCHI✔️❌:         if (err)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* crossing bonds in blist are invalid */
    // INCHI✔️❌:             err = 12;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Check each unit for >1 connected components */
    // INCHI✔️❌: #if 0
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int icompnt;
    // INCHI✔️❌:             icompnt = compnt_nums[u->alist[0] - 1];
    // INCHI✔️❌:             for (k = 1; k < u->na; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (compnt_nums[u->alist[k] - 1] != icompnt)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     u->disjoint = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Sort all units in modified alist's lexicographic order
    // INCHI✔️❌:     (modification is: longer list always go first )\t\t\t*/
    // INCHI✔️❌:     for (i = 0; i < p->n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         unum[i] = i;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 1; i < p->n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int tmp = unum[i];
    // INCHI✔️❌:         int j = i - 1;
    // INCHI✔️❌:         while (j >= 0 && OAD_PolymerUnit_CompareAtomListsMod( units2[unum[j]], units2[tmp] ) > 0)
    // INCHI✔️❌:         /*while ( j >= 0 &&    OAD_PolymerUnit_CompareAtomLists( units2[ unum[j] ], units2[ tmp ] ) > 0  )*/
    // INCHI✔️❌:         {
    // INCHI✔️❌:             unum[j + 1] = unum[j];
    // INCHI✔️❌:             j--;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         unum[j + 1] = tmp;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     return err;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OAD_Polymer_PrepareWorkingSet

    let canonical_values = if polymer.n_pzz > 0 || polymer.n > 0 {
        heap.slice(canonical_numbers)?.to_vec()
    } else {
        Vec::new()
    };
    let map = |number: i32| -> Result<i32, SourceHeapError> {
        let index = usize::try_from(number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        Ok(canonical_values
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .wrapping_add(1))
    };
    for index in 0..polymer.n_pzz {
        let pointer = polymer.pzz.offset(i64::from(index))?;
        let original = *heap
            .slice(pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mapped = map(original)?;
        if mapped == 0 {
            return Ok(10);
        }
        heap.slice_mut(pointer)?[0] = mapped;
    }

    for index in 0..polymer.n {
        let unit_pointer = *heap
            .slice(units2.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let unit = heap
            .slice(unit_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let mut kept = 0_i32;
        for atom_index in 0..unit.na {
            let atom_pointer = unit.alist.offset(i64::from(atom_index))?;
            let original = *heap
                .slice(atom_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mapped = map(original)?;
            if mapped != 0 {
                heap.slice_mut(unit.alist.offset(i64::from(kept))?)?[0] = mapped;
                kept = kept.wrapping_add(1);
            }
        }
        heap.slice_mut(unit_pointer)?[0].na = kept;
        for bond_index in 0..unit.nb.wrapping_mul(2) {
            let pointer = unit.blist.offset(i64::from(bond_index))?;
            let original = *heap
                .slice(pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mapped = map(original)?;
            if mapped == 0 {
                return Ok(11);
            }
            heap.slice_mut(pointer)?[0] = mapped;
        }
        for field_index in 0..4 {
            let original = {
                let current = &heap.slice(unit_pointer.as_const())?[0];
                match field_index {
                    0 => current.cap1,
                    1 => current.cap2,
                    2 => current.end_atom1,
                    _ => current.end_atom2,
                }
            };
            let mapped = map(original)?;
            if mapped == 0 {
                return Ok(11);
            }
            let current = &mut heap.slice_mut(unit_pointer)?[0];
            match field_index {
                0 => current.cap1 = mapped,
                1 => current.cap2 = mapped,
                2 => current.end_atom1 = mapped,
                _ => current.end_atom2 = mapped,
            }
        }
        for bond_index in 0..unit.nbkbonds {
            let row_pointer = *heap
                .slice(unit.bkbonds.as_const().offset(i64::from(bond_index))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let first = *heap
                .slice(row_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let first_mapped = map(first)?;
            if first_mapped == 0 {
                continue;
            }
            let second = *heap
                .slice(row_pointer.as_const().offset(1)?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let second_mapped = map(second)?;
            if second_mapped == 0 {
                continue;
            }
            let row = heap.slice_mut(row_pointer)?;
            row[0] = first_mapped.min(second_mapped);
            row[1] = first_mapped.max(second_mapped);
        }
    }

    for index in 0..polymer.n {
        let unit_pointer = *heap
            .slice(units2.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let unit = heap
            .slice(unit_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        iisort(heap, unit.alist, unit.na)?;
        let mut sorted_unit = unit;
        let order_error = OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves(
            heap,
            &mut sorted_unit,
            polymer.n_pzz,
            polymer.pzz.as_const(),
        )?;
        if order_error != 0 {
            return Ok(12);
        }
        heap.slice_mut(unit_pointer)?[0] = sorted_unit;
    }

    for index in 0..polymer.n {
        heap.slice_mut(unit_numbers.offset(i64::from(index))?)?[0] = index;
    }
    for index in 1..polymer.n {
        let temporary = *heap
            .slice(unit_numbers.offset(i64::from(index))?.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut previous = index - 1;
        while previous >= 0 {
            let left_number = *heap
                .slice(unit_numbers.offset(i64::from(previous))?.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let left_pointer = *heap
                .slice(units2.offset(i64::from(left_number))?.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let right_pointer = *heap
                .slice(units2.offset(i64::from(temporary))?.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if OAD_PolymerUnit_CompareAtomListsMod(
                heap,
                left_pointer.as_const(),
                right_pointer.as_const(),
            )? <= 0
            {
                break;
            }
            heap.slice_mut(unit_numbers.offset(i64::from(previous + 1))?)?[0] = left_number;
            previous -= 1;
        }
        heap.slice_mut(unit_numbers.offset(i64::from(previous + 1))?)?[0] = temporary;
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OrigAtData_RemoveHalfBond(
    heap: &mut SourceHeap,
    this_atom: i32,
    other_atom: i32,
    atoms: SourceMutPointer<inp_ATOM>,
    bond_type: &mut i32,
    bond_stereo: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2517 OrigAtData_RemoveHalfBond
    // INCHI✔️✔️: int  OrigAtData_RemoveHalfBond( int      this_atom,
    // INCHI✔️✔️:                                 int      other_atom,
    // INCHI✔️✔️:                                 inp_ATOM *at,
    // INCHI✔️✔️:                                 int      *bond_type,
    // INCHI✔️✔️:                                 int      *bond_stereo )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int k, kk;
    // INCHI✔️✔️:     /* djb-rwth: fixing oss-fuzz issues #68286, #30342 */
    // INCHI✔️✔️:     if (at && (this_atom >= 0) && (other_atom >= 0))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inp_ATOM* a = &(at[this_atom]);
    // INCHI✔️✔️:         if (a)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             for (k = 0; k < a->valence; k++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (a->neighbor[k] != other_atom)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     continue;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 *bond_type = a->bond_type[k];
    // INCHI✔️✔️:                 *bond_stereo = a->bond_stereo[k];
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 a->neighbor[k] = a->bond_type[k] = a->bond_stereo[k] = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 for (kk = k + 1; kk < a->valence; kk++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     a->neighbor[kk - 1] = a->neighbor[kk];
    // INCHI✔️✔️:                     a->bond_type[kk - 1] = a->bond_type[kk];
    // INCHI✔️✔️:                     a->bond_stereo[kk - 1] = a->bond_stereo[kk];
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 for (kk = a->valence - 1; kk < MAXVAL; kk++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     a->neighbor[kk] = 0;
    // INCHI✔️✔️:                     a->bond_type[kk] = (U_CHAR)0;
    // INCHI✔️✔️:                     a->bond_stereo[kk] = (S_CHAR)0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 return 1;
    // INCHI✔️✔️:             } /* k */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OrigAtData_RemoveHalfBond

    if atoms.is_null() || this_atom < 0 || other_atom < 0 {
        return Ok(0);
    }
    let atom_pointer = atoms.offset(i64::from(this_atom))?;
    let atom = heap
        .slice_mut(atom_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let valence = i32::from(atom.valence);
    for k in 0..valence {
        let index = usize::try_from(k).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if i32::from(
            *atom
                .neighbor
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        ) != other_atom
        {
            continue;
        }
        *bond_type = i32::from(atom.bond_type[index]);
        *bond_stereo = i32::from(atom.bond_stereo[index]);
        atom.neighbor[index] = 0;
        atom.bond_type[index] = 0;
        atom.bond_stereo[index] = 0;
        for kk in (k + 1)..valence {
            let source = usize::try_from(kk).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let target = source - 1;
            atom.neighbor[target] = atom.neighbor[source];
            atom.bond_type[target] = atom.bond_type[source];
            atom.bond_stereo[target] = atom.bond_stereo[source];
        }
        for kk in (valence - 1)..20 {
            let index = usize::try_from(kk).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            atom.neighbor[index] = 0;
            atom.bond_type[index] = 0;
            atom.bond_stereo[index] = 0;
        }
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OrigAtData_RemoveBond(
    heap: &mut SourceHeap,
    this_atom: i32,
    other_atom: i32,
    atoms: SourceMutPointer<inp_ATOM>,
    bond_type: &mut i32,
    bond_stereo: &mut i32,
    num_inp_bonds: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2577 OrigAtData_RemoveBond
    // INCHI✔️✔️: int  OrigAtData_RemoveBond( int      this_atom,
    // INCHI✔️✔️:                             int      other_atom,
    // INCHI✔️✔️:                             inp_ATOM *at,
    // INCHI✔️✔️:                             int      *bond_type,
    // INCHI✔️✔️:                             int      *bond_stereo,
    // INCHI✔️✔️:                             int      *num_inp_bonds )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int del = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (at && (this_atom >= 0) && (other_atom >= 0)) /* djb-rwth: fixing oss-fuzz issue #68329, #68286 */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         del = OrigAtData_RemoveHalfBond(this_atom, other_atom, at, bond_type, bond_stereo);
    // INCHI✔️✔️:         del += OrigAtData_RemoveHalfBond(other_atom, this_atom, at, bond_type, bond_stereo);
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (del == 2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             (*num_inp_bonds)--;
    // INCHI✔️✔️:             at[this_atom].valence--;
    // INCHI✔️✔️:             at[this_atom].chem_bonds_valence -= *bond_type;
    // INCHI✔️✔️:             at[other_atom].valence--;
    // INCHI✔️✔️:             at[other_atom].chem_bonds_valence -= *bond_type;
    // INCHI✔️✔️:             return 1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OrigAtData_RemoveBond

    if atoms.is_null() || this_atom < 0 || other_atom < 0 {
        return Ok(0);
    }
    let first =
        OrigAtData_RemoveHalfBond(heap, this_atom, other_atom, atoms, bond_type, bond_stereo)?;
    let second =
        OrigAtData_RemoveHalfBond(heap, other_atom, this_atom, atoms, bond_type, bond_stereo)?;
    if first + second == 2 {
        *num_inp_bonds = num_inp_bonds.wrapping_sub(1);
        let first_pointer = atoms.offset(i64::from(this_atom))?;
        let second_pointer = atoms.offset(i64::from(other_atom))?;
        {
            let first_atom = heap
                .slice_mut(first_pointer)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            first_atom.valence = first_atom.valence.wrapping_sub(1);
            first_atom.chem_bonds_valence =
                first_atom.chem_bonds_valence.wrapping_sub(*bond_type as i8);
        }
        {
            let second_atom = heap
                .slice_mut(second_pointer)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            second_atom.valence = second_atom.valence.wrapping_sub(1);
            second_atom.chem_bonds_valence = second_atom
                .chem_bonds_valence
                .wrapping_sub(*bond_type as i8);
        }
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OrigAtData_AddBond(
    heap: &mut SourceHeap,
    this_atom: i32,
    other_atom: i32,
    atoms: SourceMutPointer<inp_ATOM>,
    mut bond_type: i32,
    bond_stereo: i32,
    num_bonds: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2607 OrigAtData_AddBond
    // INCHI✔️✔️: int  OrigAtData_AddBond( int        this_atom,
    // INCHI✔️✔️:                          int        other_atom,
    // INCHI✔️✔️:                          inp_ATOM   *at,
    // INCHI✔️✔️:                          int        bond_type,
    // INCHI✔️✔️:                          int        bond_stereo,
    // INCHI✔️✔️:                          int        *num_bonds )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (at)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* djb-rwth: fixing oss-fuzz issue #68286 */
    // INCHI✔️✔️:         int i, k, already_here;
    // INCHI✔️✔️:         inp_ATOM* a = &(at[this_atom]);
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (at[this_atom].valence >= MAXVAL ||
    // INCHI✔️✔️:             at[other_atom].valence >= MAXVAL)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (bond_type != INCHI_BOND_TYPE_DOUBLE && bond_type != INCHI_BOND_TYPE_TRIPLE)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         k = a->valence;
    // INCHI✔️✔️:         already_here = 0;
    // INCHI✔️✔️:         for (i = 0; i < k; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (a->neighbor[i] == other_atom)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 already_here = 1; break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (!already_here)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a->neighbor[k] = other_atom;
    // INCHI✔️✔️:             a->bond_type[k] = (U_CHAR)bond_type;
    // INCHI✔️✔️:             a->bond_stereo[k] = (S_CHAR)bond_stereo;
    // INCHI✔️✔️:             a->chem_bonds_valence += bond_type;
    // INCHI✔️✔️:             a->valence++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         a = &(at[other_atom]);
    // INCHI✔️✔️:         k = a->valence;
    // INCHI✔️✔️:         already_here = 0;
    // INCHI✔️✔️:         for (i = 0; i < k; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (a->neighbor[i] == this_atom)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 already_here = 1; break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (!already_here && (k < MAXVAL)) /* djb-rwth: condition added to prevent buffer overrun */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a->neighbor[k] = this_atom;
    // INCHI✔️✔️:             a->bond_type[k] = (U_CHAR)bond_type;
    // INCHI✔️✔️:             a->bond_stereo[k] = (S_CHAR)bond_stereo;
    // INCHI✔️✔️:             a->chem_bonds_valence += bond_type;
    // INCHI✔️✔️:             a->valence++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         (*num_bonds)++;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OrigAtData_AddBond

    if atoms.is_null() {
        return Ok(0);
    }
    let first_pointer = atoms.offset(i64::from(this_atom))?;
    let second_pointer = atoms.offset(i64::from(other_atom))?;
    let first_valence = heap
        .slice(first_pointer.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .valence;
    let second_valence = heap
        .slice(second_pointer.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .valence;
    if first_valence >= 20 || second_valence >= 20 {
        return Ok(0);
    }
    if bond_type != tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE as i32
        && bond_type != tagINCHIBondType_INCHI_BOND_TYPE_TRIPLE as i32
    {
        bond_type = tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i32;
    }
    {
        let atom = heap
            .slice_mut(first_pointer)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let k =
            usize::try_from(atom.valence).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let already_here = atom.neighbor[..k]
            .iter()
            .any(|neighbor| i32::from(*neighbor) == other_atom);
        if !already_here {
            atom.neighbor[k] = other_atom as u16;
            atom.bond_type[k] = bond_type as u8;
            atom.bond_stereo[k] = bond_stereo as i8;
            atom.chem_bonds_valence = atom.chem_bonds_valence.wrapping_add(bond_type as i8);
            atom.valence = atom.valence.wrapping_add(1);
        }
    }
    {
        let atom = heap
            .slice_mut(second_pointer)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let k =
            usize::try_from(atom.valence).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let already_here = atom.neighbor[..k]
            .iter()
            .any(|neighbor| i32::from(*neighbor) == this_atom);
        if !already_here && k < 20 {
            atom.neighbor[k] = this_atom as u16;
            atom.bond_type[k] = bond_type as u8;
            atom.bond_stereo[k] = bond_stereo as i8;
            atom.chem_bonds_valence = atom.chem_bonds_valence.wrapping_add(bond_type as i8);
            atom.valence = atom.valence.wrapping_add(1);
        }
    }
    *num_bonds = num_bonds.wrapping_add(1);
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn UnMarkRingSystemsInp(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1961 UnMarkRingSystemsInp
    // INCHI✔️✔️: int UnMarkRingSystemsInp( inp_ATOM *at, int num_atoms )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         at[i].bCutVertex = 0;
    // INCHI✔️✔️:         at[i].nRingSystem = 0;
    // INCHI✔️✔️:         at[i].nNumAtInRingSystem = 0;
    // INCHI✔️✔️:         at[i].nBlockSystem = 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: UnMarkRingSystemsInp

    for index in 0..num_atoms {
        let atom = heap
            .slice_mut(atoms.offset(i64::from(index))?)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        atom.bCutVertex = 0;
        atom.nRingSystem = 0;
        atom.nNumAtInRingSystem = 0;
        atom.nBlockSystem = 0;
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OrigAtData_AddSingleStereolessBond(
    heap: &mut SourceHeap,
    this_atom: i32,
    other_atom: i32,
    atoms: SourceMutPointer<inp_ATOM>,
    num_bonds: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2682 OrigAtData_AddSingleStereolessBond
    // INCHI✔️✔️: int  OrigAtData_AddSingleStereolessBond( int      this_atom,
    // INCHI✔️✔️:                                          int      other_atom,
    // INCHI✔️✔️:                                          inp_ATOM *at,
    // INCHI✔️✔️:                                          int      *num_bonds )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     return OrigAtData_AddBond( this_atom, other_atom, at, INCHI_BOND_TYPE_SINGLE, 0, num_bonds );
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OrigAtData_AddSingleStereolessBond

    OrigAtData_AddBond(
        heap,
        this_atom,
        other_atom,
        atoms,
        tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i32,
        0,
        num_bonds,
    )
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_FindRingSystems(
    heap: &mut SourceHeap,
    polymer: SourceConstPointer<OAD_Polymer>,
    atoms: SourceMutPointer<inp_ATOM>,
    nat: i32,
    num_inp_bonds: &mut i32,
    num_ring_sys: SourceMutPointer<i32>,
    size_ring_sys: SourceMutPointer<i32>,
    start: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3236 OAD_Polymer_FindRingSystems
    // INCHI✔️❌: int  OAD_Polymer_FindRingSystems( OAD_Polymer *pd,
    // INCHI✔️❌:                                   inp_ATOM    *at,
    // INCHI✔️❌:                                   int         nat,
    // INCHI✔️❌:                                   int         *num_inp_bonds,
    // INCHI✔️❌:                                   int         *num_ring_sys,
    // INCHI✔️❌:                                   int         *size_ring_sys,
    // INCHI✔️❌:                                   int         start )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, nrings = 0, bond_type, bond_stereo;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (NULL == num_ring_sys)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Remove polymer SRU 'cyclizing' bonds if any */
    // INCHI✔️❌:     for (j = 0; j < pd->n; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pd->units[j]->cyclized)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             OrigAtData_RemoveBond( pd->units[j]->end_atom1 - 1,
    // INCHI✔️❌:                                    pd->units[j]->end_atom2 - 1,
    // INCHI✔️❌:                                    at, &bond_type, &bond_stereo,
    // INCHI✔️❌:                                    num_inp_bonds );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     MarkRingSystemsInp( at, nat, start ); /*0 );*/
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i <= nat; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_ring_sys[i] = -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 0; i < nat; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (at[i].nNumAtInRingSystem > 2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int atnum = at[i].orig_at_number;
    // INCHI✔️❌:             num_ring_sys[atnum] = at[i].nRingSystem;
    // INCHI✔️❌:             if (NULL != size_ring_sys)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 size_ring_sys[atnum] = at[i].nNumAtInRingSystem;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     UnMarkRingSystemsInp( at, nat );
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < nat; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (num_ring_sys[i] > -1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nrings++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Restore polymer SRU 'cyclizing' bonds if applicable */
    // INCHI✔️❌:     for (j = 0; j < pd->n; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pd->units[j]->cyclized)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             OrigAtData_AddSingleStereolessBond( pd->units[j]->end_atom1 - 1,
    // INCHI✔️❌:                                                 pd->units[j]->end_atom2 - 1,
    // INCHI✔️❌:                                                 at,
    // INCHI✔️❌:                                                 num_inp_bonds );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return nrings;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OAD_Polymer_FindRingSystems

    if num_ring_sys.is_null() {
        return Ok(0);
    }
    let polymer_value = heap
        .slice(polymer)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let mut bond_type = 0;
    let mut bond_stereo = 0;
    for index in 0..polymer_value.n {
        let unit_pointer = *heap
            .slice(polymer_value.units.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let unit = heap
            .slice(unit_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        if unit.cyclized != 0 {
            let _ = OrigAtData_RemoveBond(
                heap,
                unit.end_atom1.wrapping_sub(1),
                unit.end_atom2.wrapping_sub(1),
                atoms,
                &mut bond_type,
                &mut bond_stereo,
                num_inp_bonds,
            )?;
        }
    }
    let _ = MarkRingSystemsInp(heap, atoms, nat, start)?;
    for index in 0..=nat {
        heap.slice_mut(num_ring_sys.offset(i64::from(index))?)?[0] = -1;
    }
    for index in 0..nat {
        let atom = heap
            .slice(atoms.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        if atom.nNumAtInRingSystem > 2 {
            let atom_index = i32::from(atom.orig_at_number);
            heap.slice_mut(num_ring_sys.offset(i64::from(atom_index))?)?[0] =
                i32::from(atom.nRingSystem);
            if !size_ring_sys.is_null() {
                heap.slice_mut(size_ring_sys.offset(i64::from(atom_index))?)?[0] =
                    i32::from(atom.nNumAtInRingSystem);
            }
        }
    }
    UnMarkRingSystemsInp(heap, atoms, nat)?;
    let mut rings = 0;
    for index in 0..nat {
        if heap
            .slice(num_ring_sys.as_const().offset(i64::from(index))?)?
            .first()
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            > -1
        {
            rings += 1;
        }
    }
    for index in 0..polymer_value.n {
        let unit_pointer = *heap
            .slice(polymer_value.units.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let unit = heap
            .slice(unit_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        if unit.cyclized != 0 {
            let _ = OrigAtData_AddSingleStereolessBond(
                heap,
                unit.end_atom1.wrapping_sub(1),
                unit.end_atom2.wrapping_sub(1),
                atoms,
                num_inp_bonds,
            )?;
        }
    }
    Ok(rings)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_SetAtProps(
    heap: &mut SourceHeap,
    polymer: SourceConstPointer<OAD_Polymer>,
    atoms: SourceMutPointer<inp_ATOM>,
    nat: i32,
    num_inp_bonds: &mut i32,
    atom_properties: SourceMutPointer<OAD_AtProps>,
    canonical_numbers: SourceConstPointer<i32>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3312 OAD_Polymer_SetAtProps
    // INCHI✔️❌: void OAD_Polymer_SetAtProps( OAD_Polymer *pd,
    // INCHI✔️❌:                              inp_ATOM    *at,
    // INCHI✔️❌:                              int         nat,
    // INCHI✔️❌:                              int         *num_inp_bonds,
    // INCHI✔️❌:                              OAD_AtProps *aprops,
    // INCHI✔️❌:                              int         *cano_nums )
    // INCHI✔️❌: {
    // INCHI✔️❌: /*  Max rank for in-ring atom is 216 which is achieved for N (element number 7 in Periodic system & erank_rule2[] ),*/
    // INCHI✔️❌: /*    then goes O with rank 215 (element number 8), and so on... lowest rank is 1 for H .                           */
    // INCHI✔️❌: /*                                                                                                                  */
    // INCHI✔️❌: /*  This follows to IUPAC rule 2 [Pure Appl. Chem., Vol. 74, No. 10, 2002, p. 1926] which states:                   */
    // INCHI✔️❌: /*  a. a ring or ring system containing nitrogen;                                                                   */
    // INCHI✔️❌: /*  b. a ring or ring system containing the heteroatom occurring earliest in the order given in Rule 4;             */
    // INCHI✔️❌: /*  ( which is     O > S > Se > Te > N > P > As > Sb > Bi > Si > Ge > Sn > Pb > B > Hg )                            */
    // INCHI✔️❌: /*  ...                                                                                                             */
    // INCHI✔️❌:
    // INCHI✔️❌:     int erank_rule2[] = { 0,1,198,197,196,202,2,216,215,191,190,189,188,187,206,210,214,183,182,181,180,179,178,177,176,
    // INCHI✔️❌:                           175,174,173,172,171,170,169,205,209,213,165,164,163,162,161,160,159,158,157,156,155,154,153,152,
    // INCHI✔️❌:                           151,204,208,212,147,146,145,144,143,142,141,140,139,138,137,136,135,134,133,132,131,130,129,128,
    // INCHI✔️❌:                           127,126,125,124,123,122,121,201,119,203,207,116,115,114,113,112,111,110,109,108,107,106,105,104,
    // INCHI✔️❌:                           103,102,101,100,99,98,97,96,95,94,93,92,91,90,89,88,87,86,85,84,83,82,81 };
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  Max rank for chain atom is 215 which is achieved for O (element number 8 in Periodic system & erank_rule4[] ),  */
    // INCHI✔️❌:     /*  then goes N with rank 212 (element number 8), and so on... lowest rank is 1 for H .                             */
    // INCHI✔️❌:     /*                                                                                                                  */
    // INCHI✔️❌:     /*  This follows to IUPAC rule 4 [Pure Appl. Chem., Vol. 74, No. 10, 2002, p. 1927] which states:                   */
    // INCHI✔️❌:     /*  O > S > Se > Te > N > P > As > Sb > Bi > Si > Ge > Sn > Pb > B > Hg                                             */
    // INCHI✔️❌:     /*  Note: Other heteroatoms may be placed within this order as indicated by their positions in the                  */
    // INCHI✔️❌:     /*  periodic table [5].                                                                                             */
    // INCHI✔️❌:
    // INCHI✔️❌:     int erank_rule4[] = { 0,1,198,197,196,202,2,211,215,191,190,189,188,187,206,210,214,183,182,181,180,179,178,177,176,
    // INCHI✔️❌:                           175,174,173,172,171,170,169,205,209,213,165,164,163,162,161,160,159,158,157,156,155,154,153,152,
    // INCHI✔️❌:                           151,204,208,212,147,146,145,144,143,142,141,140,139,138,137,136,135,134,133,132,131,130,129,128,
    // INCHI✔️❌:                           127,126,125,124,123,122,121,201,119,203,207,116,115,114,113,112,111,110,109,108,107,106,105,104,
    // INCHI✔️❌:                           103,102,101,100,99,98,97,96,95,94,93,92,91,90,89,88,87,86,85,84,83,82,81 };
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     int i, j, k, nrings = 0;
    // INCHI✔️❌:     int a1, a2, dummy = 0, bond_type = -1, bond_stereo = -1;
    // INCHI✔️❌:     int *num_ring_sys = NULL, *size_ring_sys = NULL;
    // INCHI✔️❌:     /* djb-rwth: fixing oss-fuzz issue #68112 */
    // INCHI✔️❌:     int err2_len = sizeof(erank_rule2) / sizeof(erank_rule2[0]);
    // INCHI✔️❌:     int err4_len = sizeof(erank_rule4) / sizeof(erank_rule4[0]);
    // INCHI✔️❌:
    // INCHI✔️❌:     if ((NULL == aprops) || !at || !pd) /* djb-rwth: fixing oss-fuzz issue #68329, #68286 */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Establish element ranks for atoms */
    // INCHI✔️❌:     for (k = 0; k < nat; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int atnum = at[k].orig_at_number, index = k;
    // INCHI✔️❌:         U_CHAR err4_ind = at[k].el_number;
    // INCHI✔️❌:         if (cano_nums)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             index = cano_nums[atnum];
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (index >= 0 && err4_ind < err4_len)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             aprops[index].erank = erank_rule4[err4_ind];
    // INCHI✔️❌:             aprops[index].ring_erank = 0;
    // INCHI✔️❌:             aprops[index].ring_size = 0;
    // INCHI✔️❌:             aprops[index].ring_num = -1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* deleted H's go here */
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Establish ring systems assignments for atoms */
    // INCHI✔️❌:     num_ring_sys = (int *) inchi_calloc( (long long)nat + 1, sizeof( int ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     if (NULL == num_ring_sys)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     size_ring_sys = (int *) inchi_calloc( (long long)nat + 1, sizeof( int ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     if (NULL == size_ring_sys)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Note that we get here on the way of InChI2Struct conversion.            */
    // INCHI✔️❌:     /* Break temporarily any of (actually, the first) SRU 'cyclizing' bonds    */
    // INCHI✔️❌:     for (j = 0; j < pd->n; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pd->units[j]->na > 2 && pd->units[j]->nbkbonds > 0 &&
    // INCHI✔️❌:              pd->units[j]->cyclized == 0 &&
    // INCHI✔️❌:              pd->units[j]->cyclizable == CLOSING_SRU_RING)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             a1 = pd->units[j]->bkbonds[0][0] - 1;
    // INCHI✔️❌:             a2 = pd->units[j]->bkbonds[0][1] - 1;
    // INCHI✔️❌:             OrigAtData_RemoveBond( a1, a2, at, &bond_type, &bond_stereo, &dummy );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nrings = OAD_Polymer_FindRingSystems( pd, at, nat, num_inp_bonds, num_ring_sys, size_ring_sys, 0 );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Immediately restore just broken bond(s) */
    // INCHI✔️❌:     for (j = 0; j < pd->n; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pd->units[j]->na > 2 &&
    // INCHI✔️❌:              pd->units[j]->nbkbonds > 0 &&
    // INCHI✔️❌:              pd->units[j]->cyclized == 0 &&
    // INCHI✔️❌:              pd->units[j]->cyclizable == CLOSING_SRU_RING)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             a1 = pd->units[j]->bkbonds[0][0] - 1;
    // INCHI✔️❌:             a2 = pd->units[j]->bkbonds[0][1] - 1;
    // INCHI✔️❌:             /* OrigAtData_AddSingleStereolessBond( a1, a2, at, &dummy ); */
    // INCHI✔️❌:             OrigAtData_AddBond( a1, a2, at, bond_type, bond_stereo, &dummy );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nrings)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int max_ring_num = 0;
    // INCHI✔️❌:         /* SRU contains ring[s], proceed with them following (not totally) the IUPAC guidelines */
    // INCHI✔️❌:         for (k = 0; k < nat; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* Browse 0-based original atoms, go to 1-based cano nums domain if cano_nums mapping is suppied */
    // INCHI✔️❌:             int atnum = at[k].orig_at_number, index = k;
    // INCHI✔️❌:             if (cano_nums)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 index = cano_nums[atnum] + 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (num_ring_sys[atnum] >= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 aprops[index].ring_num = num_ring_sys[atnum];  /* temporarily */
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (max_ring_num < aprops[index].ring_num)
    // INCHI✔️❌:                     max_ring_num = aprops[index].ring_num;          /* NB: OAD_Polymer_FindRingSystems may return num_ring_sys[]  */
    // INCHI✔️❌:                                                                     /* which is not a list of consecutive numbers                 */
    // INCHI✔️❌:
    // INCHI✔️❌:                 aprops[index].ring_size = size_ring_sys[atnum];     /* Size of ring system which includes the atom k .            */
    // INCHI✔️❌:                                                                     /* It is used as an additional score for in-ring              */
    // INCHI✔️❌:                                                                     /* atoms' prioritizing (instead of criteria in                */
    // INCHI✔️❌:                                                                     /* 2c-2h of IUPAC rule 2 which deal with ring sizes).         */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0; i <= max_ring_num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int erank, max_erank = 0;
    // INCHI✔️❌:             for (k = 0; k < nat; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int atnum = at[k].orig_at_number, index = k;
    // INCHI✔️❌:                 if (cano_nums)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     index = cano_nums[atnum] + 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (aprops[index].ring_num == i)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     erank = erank_rule2[at[k].el_number];
    // INCHI✔️❌:                     if (erank > max_erank)
    // INCHI✔️❌:                         max_erank = erank;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (k = 0; k < nat; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int atnum = at[k].orig_at_number, index = k;
    // INCHI✔️❌:                 if (cano_nums)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     index = cano_nums[atnum] + 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (aprops[index].ring_num == i)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (aprops[index].ring_size > 2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         aprops[index].ring_erank = max_erank;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     if (num_ring_sys)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( num_ring_sys );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (size_ring_sys)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( size_ring_sys );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OAD_Polymer_SetAtProps

    const ERANK_RULE2: &[i32] = &[
        0, 1, 198, 197, 196, 202, 2, 216, 215, 191, 190, 189, 188, 187, 206, 210, 214, 183, 182,
        181, 180, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 169, 205, 209, 213, 165, 164,
        163, 162, 161, 160, 159, 158, 157, 156, 155, 154, 153, 152, 151, 204, 208, 212, 147, 146,
        145, 144, 143, 142, 141, 140, 139, 138, 137, 136, 135, 134, 133, 132, 131, 130, 129, 128,
        127, 126, 125, 124, 123, 122, 121, 201, 119, 203, 207, 116, 115, 114, 113, 112, 111, 110,
        109, 108, 107, 106, 105, 104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90,
        89, 88, 87, 86, 85, 84, 83, 82, 81,
    ];
    const ERANK_RULE4: &[i32] = &[
        0, 1, 198, 197, 196, 202, 2, 211, 215, 191, 190, 189, 188, 187, 206, 210, 214, 183, 182,
        181, 180, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 169, 205, 209, 213, 165, 164,
        163, 162, 161, 160, 159, 158, 157, 156, 155, 154, 153, 152, 151, 204, 208, 212, 147, 146,
        145, 144, 143, 142, 141, 140, 139, 138, 137, 136, 135, 134, 133, 132, 131, 130, 129, 128,
        127, 126, 125, 124, 123, 122, 121, 201, 119, 203, 207, 116, 115, 114, 113, 112, 111, 110,
        109, 108, 107, 106, 105, 104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90,
        89, 88, 87, 86, 85, 84, 83, 82, 81,
    ];

    if atom_properties.is_null() || atoms.is_null() || polymer.is_null() {
        return Ok(());
    }
    let canonical_index = |heap: &SourceHeap, atom_number: i32| -> Result<i32, SourceHeapError> {
        if canonical_numbers.is_null() {
            return Ok(atom_number);
        }
        heap.slice(canonical_numbers.offset(i64::from(atom_number))?)?
            .first()
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    };
    for k in 0..nat {
        let atom = heap
            .slice(atoms.as_const().offset(i64::from(k))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let atom_number = i32::from(atom.orig_at_number);
        let index = canonical_index(heap, atom_number)?;
        let element = usize::from(atom.el_number);
        if index >= 0 && element < ERANK_RULE4.len() {
            let properties = heap
                .slice_mut(atom_properties.offset(i64::from(index))?)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            properties.erank = ERANK_RULE4[element];
            properties.ring_erank = 0;
            properties.ring_size = 0;
            properties.ring_num = -1;
        }
    }

    let allocation_count = match usize::try_from(i64::from(nat) + 1) {
        Ok(count) => count,
        Err(_) => return Ok(()),
    };
    let ring_numbers = match heap.allocate(vec![0_i32; allocation_count]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(()),
        Err(error) => return Err(error),
    };
    let ring_sizes = match heap.allocate(vec![0_i32; allocation_count]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => {
            inchi_free(heap, ring_numbers)?;
            return Ok(());
        }
        Err(error) => {
            inchi_free(heap, ring_numbers)?;
            return Err(error);
        }
    };

    let result = (|| -> Result<(), SourceHeapError> {
        let polymer_value = heap
            .slice(polymer)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let mut dummy = 0_i32;
        let mut bond_type = -1_i32;
        let mut bond_stereo = -1_i32;
        for j in 0..polymer_value.n {
            let unit_pointer = *heap
                .slice(polymer_value.units.as_const().offset(i64::from(j))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let unit = heap
                .slice(unit_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if unit.na > 2
                && unit.nbkbonds > 0
                && unit.cyclized == 0
                && unit.cyclizable == CLOSING_SRU_RING as i32
            {
                let row = *heap
                    .slice(unit.bkbonds.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let bond = heap.slice(row.as_const())?;
                let a1 = bond[0].wrapping_sub(1);
                let a2 = bond[1].wrapping_sub(1);
                let _ = OrigAtData_RemoveBond(
                    heap,
                    a1,
                    a2,
                    atoms,
                    &mut bond_type,
                    &mut bond_stereo,
                    &mut dummy,
                )?;
            }
        }

        let rings = OAD_Polymer_FindRingSystems(
            heap,
            polymer,
            atoms,
            nat,
            num_inp_bonds,
            ring_numbers,
            ring_sizes,
            0,
        )?;

        for j in 0..polymer_value.n {
            let unit_pointer = *heap
                .slice(polymer_value.units.as_const().offset(i64::from(j))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let unit = heap
                .slice(unit_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if unit.na > 2
                && unit.nbkbonds > 0
                && unit.cyclized == 0
                && unit.cyclizable == CLOSING_SRU_RING as i32
            {
                let row = *heap
                    .slice(unit.bkbonds.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let bond = heap.slice(row.as_const())?;
                let a1 = bond[0].wrapping_sub(1);
                let a2 = bond[1].wrapping_sub(1);
                let _ =
                    OrigAtData_AddBond(heap, a1, a2, atoms, bond_type, bond_stereo, &mut dummy)?;
            }
        }

        if rings != 0 {
            let mut max_ring_number = 0_i32;
            for k in 0..nat {
                let atom = heap
                    .slice(atoms.as_const().offset(i64::from(k))?)?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                let atom_number = i32::from(atom.orig_at_number);
                let index = if canonical_numbers.is_null() {
                    k
                } else {
                    canonical_index(heap, atom_number)?.wrapping_add(1)
                };
                let ring_number =
                    heap.slice(ring_numbers.as_const())?[usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                if ring_number >= 0 {
                    let ring_size =
                        heap.slice(ring_sizes.as_const())?[usize::try_from(atom_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                    let properties = heap
                        .slice_mut(atom_properties.offset(i64::from(index))?)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    properties.ring_num = ring_number;
                    if max_ring_number < properties.ring_num {
                        max_ring_number = properties.ring_num;
                    }
                    properties.ring_size = ring_size;
                }
            }
            for ring_number in 0..=max_ring_number {
                let mut max_erank = 0_i32;
                for k in 0..nat {
                    let atom = heap
                        .slice(atoms.as_const().offset(i64::from(k))?)?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    let atom_number = i32::from(atom.orig_at_number);
                    let index = if canonical_numbers.is_null() {
                        k
                    } else {
                        canonical_index(heap, atom_number)?.wrapping_add(1)
                    };
                    let properties = heap
                        .slice(atom_properties.as_const().offset(i64::from(index))?)?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if properties.ring_num == ring_number {
                        max_erank = max_erank.max(ERANK_RULE2[usize::from(atom.el_number)]);
                    }
                }
                for k in 0..nat {
                    let atom = heap
                        .slice(atoms.as_const().offset(i64::from(k))?)?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    let atom_number = i32::from(atom.orig_at_number);
                    let index = if canonical_numbers.is_null() {
                        k
                    } else {
                        canonical_index(heap, atom_number)?.wrapping_add(1)
                    };
                    let properties = heap
                        .slice_mut(atom_properties.offset(i64::from(index))?)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if properties.ring_num == ring_number && properties.ring_size > 2 {
                        properties.ring_erank = max_erank;
                    }
                }
            }
        }
        Ok(())
    })();

    inchi_free(heap, ring_numbers)?;
    inchi_free(heap, ring_sizes)?;
    result
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_Free(
    heap: &mut SourceHeap,
    mut pd: SourceMutPointer<OAD_Polymer>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2071 OAD_Polymer_Free
    // INCHI✔❌: void OAD_Polymer_Free( OAD_Polymer *pd )
    // INCHI✔❌: {
    // INCHI✔❌:     if (pd)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (pd->pzz)
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_free( pd->pzz );
    // INCHI✔❌:             pd->pzz = NULL;
    // INCHI✔❌:             pd->n_pzz = 0;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (pd->n && pd->units)
    // INCHI✔❌:         {
    // INCHI✔❌:             int k;
    // INCHI✔❌:             for (k = 0; k < pd->n; k++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 OAD_PolymerUnit_Free( pd->units[k] );
    // INCHI✔❌:             }
    // INCHI✔❌:             inchi_free( pd->units );
    // INCHI✔❌:             pd->units = NULL;
    // INCHI✔❌:             pd->n = 0;
    // INCHI✔❌:         }
    // INCHI✔❌:         inchi_free( pd );
    // INCHI✔❌:         pd = NULL;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: OAD_Polymer_Free

    if !pd.is_null() {
        let pzz = heap
            .slice(pd.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .pzz;
        if !pzz.is_null() {
            inchi_free(heap, pzz)?;
            let polymer = heap
                .slice_mut(pd)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            polymer.pzz = SourceMutPointer::null();
            polymer.n_pzz = 0;
        }

        let (n, units) = {
            let polymer = heap
                .slice(pd.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (polymer.n, polymer.units)
        };
        if n != 0 && !units.is_null() {
            for k in 0..n {
                let unit = heap.slice(units.as_const())?[k as usize];
                OAD_PolymerUnit_Free(heap, unit)?;
            }
            inchi_free(heap, units)?;
            let polymer = heap
                .slice_mut(pd)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            polymer.units = SourceMutPointer::null();
            polymer.n = 0;
        }
        inchi_free(heap, pd)?;
        pd = SourceMutPointer::null();
    }
    let _ = pd;
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_SetReopeningDetails(
    heap: &SourceHeap,
    unit: &mut OAD_PolymerUnit,
    atoms: &[inp_ATOM],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4007 OAD_PolymerUnit_SetReopeningDetails
    // INCHI✔️✔️: int OAD_PolymerUnit_SetReopeningDetails( OAD_PolymerUnit *u, inp_ATOM *at )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int k;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Check reopening  type */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Caps are separated by one atom - that's not error but do nothing */
    // INCHI✔️✔️:     if (u->nbkbonds == 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (u->nbkbonds == 1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         u->end_atom1 = u->bkbonds[0][0];
    // INCHI✔️✔️:         u->end_atom2 = u->bkbonds[0][1];
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (u->end_atom1 == u->end_atom2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️: #ifdef ALLOW_CLOSING_SRU_VIA_DIRADICAL
    // INCHI✔️✔️:             u->cyclizable = CLOSING_SRU_DIRADICAL;
    // INCHI✔️✔️: #else
    // INCHI✔️✔️:             u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* If caps are separated by two atoms - that's not error but do nothing */
    // INCHI✔️✔️:             for (k = 0; k < at[u->end_atom1 - 1].valence; k++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (at[u->end_atom1 - 1].neighbor[k] == u->end_atom2 - 1)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     if (at[u->end_atom1 - 1].bond_type[k] > 1)
    // INCHI✔️✔️: #ifdef ALLOW_CLOSING_SRU_VIA_HIGHER_ORDER_BOND
    // INCHI✔️✔️:                         u->cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND;
    // INCHI✔️✔️: #else
    // INCHI✔️✔️: /*                  u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;*/
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️: break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     } /*    */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return u->nbkbonds;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_SetReopeningDetails

    if unit.nbkbonds == 0 {
        return Ok(0);
    }
    if unit.nbkbonds == 1 {
        let first_row = *heap
            .slice(unit.bkbonds.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let row = heap.slice(first_row.as_const())?;
        unit.end_atom1 = *row.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
        unit.end_atom2 = *row.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if unit.end_atom1 == unit.end_atom2 {
            unit.cyclizable = CLOSING_SRU_DIRADICAL as i32;
        } else {
            let atom_index = usize::try_from(
                unit.end_atom1
                    .checked_sub(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?,
            )
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atom = atoms
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let target = unit
                .end_atom2
                .checked_sub(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let valence = usize::try_from(i32::from(atom.valence).max(0))
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            for index in 0..valence {
                if i32::from(atom.neighbor[index]) == target {
                    if atom.bond_type[index] > 1 {
                        unit.cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND as i32;
                    }
                    break;
                }
            }
        }
    }
    Ok(unit.nbkbonds)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(
    heap: &mut SourceHeap,
    unit: &mut OAD_PolymerUnit,
    _atoms: SourceMutPointer<inp_ATOM>,
    atom_properties: &[OAD_AtProps],
    senior_bond: &mut i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4056 OAD_PolymerUnit_SortBackboneBondsAndSetSeniors
    // INCHI✔️❌: void OAD_PolymerUnit_SortBackboneBondsAndSetSeniors( OAD_PolymerUnit *u,
    // INCHI✔️❌:                                                      inp_ATOM        *at,
    // INCHI✔️❌:                                                      OAD_AtProps     *aprops,
    // INCHI✔️❌:                                                      int             *senior_bond )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int j, *bnum = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     *senior_bond = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Sort backbone (== frame shiftable) bonds if necessary */
    // INCHI✔️❌:     if (u->nbkbonds > 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bnum = (int *) inchi_calloc( u->nbkbonds, sizeof( int ) );
    // INCHI✔️❌:         if (bnum)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (j = 0; j < u->nbkbonds; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bnum[j] = j;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             OAD_PolymerUnit_SortBackboneBonds( u, aprops, bnum );
    // INCHI✔️❌:             *senior_bond = bnum[0];
    // INCHI✔️❌:             inchi_free( bnum );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* v. 1.05+ : place senior atom the first ("left") in the senior bond */
    // INCHI✔️❌:     if (OAD_Polymer_IsFirstAtomRankLower( u->bkbonds[*senior_bond][0], u->bkbonds[*senior_bond][1], aprops ) == 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int tmp = u->bkbonds[*senior_bond][0];
    // INCHI✔️❌:         u->bkbonds[*senior_bond][0] = u->bkbonds[*senior_bond][1];
    // INCHI✔️❌:         u->bkbonds[*senior_bond][1] = tmp;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     u->end_atom1 = u->bkbonds[*senior_bond][0];
    // INCHI✔️❌:     u->end_atom2 = u->bkbonds[*senior_bond][1];
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_SortBackboneBondsAndSetSeniors

    *senior_bond = 0;
    if unit.nbkbonds > 1 {
        match inchi_calloc::<i32>(heap, unit.nbkbonds as u64, 4) {
            Ok(bond_number_pointer) => {
                for (index, value) in heap.slice_mut(bond_number_pointer)?.iter_mut().enumerate() {
                    *value =
                        i32::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                }
                let sort_result = OAD_PolymerUnit_SortBackboneBonds(
                    heap,
                    unit,
                    atom_properties,
                    bond_number_pointer,
                );
                if sort_result.is_ok() {
                    *senior_bond = *heap
                        .slice(bond_number_pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                }
                inchi_free(heap, bond_number_pointer)?;
                sort_result?;
            }
            Err(SourceHeapError::AllocationFailed) => {}
            Err(error) => return Err(error),
        }
    }

    let senior_index =
        usize::try_from(*senior_bond).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let row_pointer = *heap
        .slice(unit.bkbonds.as_const())?
        .get(senior_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let (first, second) = {
        let row = heap.slice(row_pointer.as_const())?;
        (
            *row.first().ok_or(SourceHeapError::PointerOutOfBounds)?,
            *row.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    if OAD_Polymer_IsFirstAtomRankLower(first, second, atom_properties)? == 1 {
        let row = heap.slice_mut(row_pointer)?;
        row[0] = second;
        row[1] = first;
    }
    let row = heap.slice(row_pointer.as_const())?;
    unit.end_atom1 = row[0];
    unit.end_atom2 = row[1];
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_SortBackboneBonds(
    heap: &mut SourceHeap,
    unit: &OAD_PolymerUnit,
    atom_properties: &[OAD_AtProps],
    bond_numbers: SourceMutPointer<i32>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4097 OAD_PolymerUnit_SortBackboneBonds
    // INCHI✔️❌: void OAD_PolymerUnit_SortBackboneBonds( OAD_PolymerUnit *u,
    // INCHI✔️❌:                                         OAD_AtProps     *aprops,
    // INCHI✔️❌:                                         int             *bnum )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, tmp;
    // INCHI✔️❌:     int n = u->nbkbonds;
    // INCHI✔️❌:     if (NULL == bnum)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 1; i < n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         tmp = bnum[i];
    // INCHI✔️❌:         j = i - 1;
    // INCHI✔️❌:         while (j >= 0 && OAD_Polymer_CompareBackboneBondsSeniority( u->bkbonds[bnum[j]], u->bkbonds[tmp], aprops ) > 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             bnum[j + 1] = bnum[j];
    // INCHI✔️❌:             j--;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         bnum[j + 1] = tmp;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_SortBackboneBonds

    if bond_numbers.is_null() {
        return Ok(());
    }
    let bond_at = |heap: &SourceHeap, number: i32| -> Result<[i32; 2], SourceHeapError> {
        let index = usize::try_from(number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let row_pointer = *heap
            .slice(unit.bkbonds.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let row = heap.slice(row_pointer.as_const())?;
        Ok([
            *row.first().ok_or(SourceHeapError::PointerOutOfBounds)?,
            *row.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?,
        ])
    };

    for i in 1..unit.nbkbonds {
        let i = usize::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let temporary = *heap
            .slice(bond_numbers.as_const())?
            .get(i)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut j = i32::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)? - 1;
        while j >= 0 {
            let j_index = usize::try_from(j).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let current = *heap
                .slice(bond_numbers.as_const())?
                .get(j_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if OAD_Polymer_CompareBackboneBondsSeniority(
                &bond_at(heap, current)?,
                &bond_at(heap, temporary)?,
                atom_properties,
            )? <= 0
            {
                break;
            }
            heap.slice_mut(bond_numbers)?[j_index + 1] = current;
            j -= 1;
        }
        let destination =
            usize::try_from(j + 1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        heap.slice_mut(bond_numbers)?[destination] = temporary;
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_CompareBackboneBondsSeniority(
    bond1: &[i32; 2],
    bond2: &[i32; 2],
    atom_properties: &[OAD_AtProps],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4128 OAD_Polymer_CompareBackboneBondsSeniority
    // INCHI✔️✔️: int  OAD_Polymer_CompareBackboneBondsSeniority( int* b1, int* b2, OAD_AtProps *aprops )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int b1min, b1max, b2min, b2max, tmp, cmp = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Find min and max ext-ranked ends of the both bonds */
    // INCHI✔️✔️:     b1max = b1[0]; b1min = b1[1];
    // INCHI✔️✔️:     b2max = b2[0]; b2min = b2[1];
    // INCHI✔️✔️:     if (OAD_Polymer_IsFirstAtomRankLower( b1min, b1max, aprops ) == -1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         tmp = b1max;
    // INCHI✔️✔️:         b1max = b1min;
    // INCHI✔️✔️:         b1min = tmp;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (OAD_Polymer_IsFirstAtomRankLower( b2min, b2max, aprops ) == -1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         tmp = b2max;
    // INCHI✔️✔️:         b2max = b2min;
    // INCHI✔️✔️:         b2min = tmp;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Compare bonds' seniority */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* First, favor the bond which has greater ext-rank end
    // INCHI✔️✔️:        NB: the result may be 0, that is, equal max ext. ranks */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     cmp = OAD_Polymer_CompareRanksOfTwoAtoms( b1max, b2max, aprops );
    // INCHI✔️✔️:     if (cmp == 1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return   1;        /* rank(b1max) < rank(b2max), so bond2 is senior */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (cmp == -1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  -1;        /* rank(b1max) > rank(b2max), so bond1 is senior */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Max ends are of the same rank, so favor the bond with lesser min-rank end
    // INCHI✔️✔️:        NB: the result may NOT be 0, that is, the case is always resolved    */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     cmp = OAD_Polymer_CompareRanksOfTwoAtoms( b1min, b2min, aprops ); /*OAD_Polymer_IsFirstAtomRankLower( b1min, b2min, aprops );*/
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (cmp == 1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  -1;         /* rank(b1min) < rank(b2min), so bond1 is senior */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (cmp == -1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return   1;         /* rank(b1min) > rank(b2min), so bond2 is senior */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Min ends are of the same rank. Here is the time to compare directly
    // INCHI✔️✔️:        which canonical number is larger of max-ends ... */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (b1max < b2max)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (b1max > b2max)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* ... they are the same, so compare which canonical number is larger for min-ends ... */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (b1min < b2min)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;          /* b1min < b2min, so bond1 is senior */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (b1min > b2min)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  1;          /* b1min > b2min, so bond2 is senior */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;    /* we should not reach there */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OAD_Polymer_CompareBackboneBondsSeniority

    let (mut bond1_max, mut bond1_min) = (bond1[0], bond1[1]);
    let (mut bond2_max, mut bond2_min) = (bond2[0], bond2[1]);
    if OAD_Polymer_IsFirstAtomRankLower(bond1_min, bond1_max, atom_properties)? == -1 {
        std::mem::swap(&mut bond1_max, &mut bond1_min);
    }
    if OAD_Polymer_IsFirstAtomRankLower(bond2_min, bond2_max, atom_properties)? == -1 {
        std::mem::swap(&mut bond2_max, &mut bond2_min);
    }

    let comparison = OAD_Polymer_CompareRanksOfTwoAtoms(bond1_max, bond2_max, atom_properties)?;
    if comparison == 1 {
        return Ok(1);
    } else if comparison == -1 {
        return Ok(-1);
    }

    let comparison = OAD_Polymer_CompareRanksOfTwoAtoms(bond1_min, bond2_min, atom_properties)?;
    if comparison == 1 {
        return Ok(-1);
    } else if comparison == -1 {
        return Ok(1);
    }

    if bond1_max < bond2_max {
        return Ok(1);
    }
    if bond1_max > bond2_max {
        return Ok(-1);
    }
    if bond1_min < bond2_min {
        return Ok(-1);
    }
    if bond1_min > bond2_min {
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_CompareRanksOfTwoAtoms(
    atom1: i32,
    atom2: i32,
    atom_properties: &[OAD_AtProps],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4209 OAD_Polymer_CompareRanksOfTwoAtoms
    // INCHI✔️✔️: int OAD_Polymer_CompareRanksOfTwoAtoms( int atom1, int atom2, OAD_AtProps *aprops )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     const int HETEROCYC = 3, HETEROAT = 2, CARBOCYC = 1, CARBOAT = 0;
    // INCHI✔️✔️:         /* NB: Carbon's rank is always 2, next to the lowest */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     int a1 = atom1 - 1;
    // INCHI✔️✔️:     int a2 = atom2 - 1;
    // INCHI✔️✔️:     int a1typ = CARBOAT;
    // INCHI✔️✔️:     int a2typ = CARBOAT;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* djb-rwth: fixing oss-fuzz issue #69501, #68277 */
    // INCHI✔️✔️:     if ((a1 < 0) || (a2 < 0))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (aprops[a1].ring_size > 2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (aprops[a1].ring_erank <= 2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a1typ = CARBOCYC;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a1typ = HETEROCYC;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (aprops[a1].erank == 2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a1typ = CARBOAT;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a1typ = HETEROAT;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (aprops[a2].ring_size > 2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (aprops[a2].ring_erank <= 2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a2typ = CARBOCYC;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a2typ = HETEROCYC;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (aprops[a2].erank == 2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a2typ = CARBOAT;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a2typ = HETEROAT;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Compare */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:         Follow IUPAC Rule 1
    // INCHI✔️✔️:             'The basic order of seniority of subunits is:
    // INCHI✔️✔️:                 heterocyclic rings and ring systems > heteroatom chains >
    // INCHI✔️✔️:                     > carbocyclic rings and ring systems > acyclic carbon chains'
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (a1typ == HETEROCYC && a2typ == HETEROCYC)   /* a1 and a2 are HETEROCYC */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* Try resolving by senior-heteroatom ring */
    // INCHI✔️✔️:         if (aprops[a1].ring_erank < aprops[a2].ring_erank)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return  1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (aprops[a1].ring_erank > aprops[a2].ring_erank)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return -1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* Same senior-heteroatom rings, try resolving by total ring size */
    // INCHI✔️✔️:         if (aprops[a1].ring_size < aprops[a2].ring_size)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return  1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (aprops[a1].ring_size > aprops[a2].ring_size)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return -1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* Could not resolve... */
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (a1typ == HETEROCYC)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;  /* a1 is HETEROCYC, a2 is any other (==junior) */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (a2typ == HETEROCYC)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  1;  /* a2 is HETEROCYC, a1 is any other (==junior) */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* HETEROCYC left out */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (a1typ == HETEROAT && a2typ == HETEROAT)  /* a1 and a2 are HETEROAT */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (aprops[a1].erank < aprops[a2].erank)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return  1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (aprops[a1].erank > aprops[a2].erank)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return -1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* Could not resolve... */
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (a1typ == HETEROAT)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;  /* a1 is HETEROAT, a2 is any other (==junior) */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (a2typ == HETEROAT)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  1;  /* a2 is HETEROAT, a1 is any other (==junior) */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* HETEROAT left out */
    // INCHI✔️✔️:     if (a1typ == CARBOCYC && a2typ == CARBOCYC) /* a1 and a2 are CARBOCYC */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* Same senior-atom (C) ring, try resolving by total ring size */
    // INCHI✔️✔️:         if (aprops[a1].ring_size < aprops[a2].ring_size)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return  1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (aprops[a1].ring_size > aprops[a2].ring_size)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return -1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* Could not resolve... */
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (a1typ == CARBOCYC)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (a2typ == CARBOCYC)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;        /* 0 means unresolved. It is legal here */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OAD_Polymer_CompareRanksOfTwoAtoms

    const HETEROCYC: i32 = 3;
    const HETEROAT: i32 = 2;
    const CARBOCYC: i32 = 1;
    const CARBOAT: i32 = 0;

    let index1 = atom1
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let index2 = atom2
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    if index1 < 0 || index2 < 0 {
        return Ok(0);
    }
    let property1 = atom_properties
        .get(usize::try_from(index1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let property2 = atom_properties
        .get(usize::try_from(index2).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let classify = |property: &OAD_AtProps| {
        if property.ring_size > 2 {
            if property.ring_erank <= 2 {
                CARBOCYC
            } else {
                HETEROCYC
            }
        } else if property.erank == 2 {
            CARBOAT
        } else {
            HETEROAT
        }
    };
    let type1 = classify(property1);
    let type2 = classify(property2);
    if type1 == HETEROCYC && type2 == HETEROCYC {
        if property1.ring_erank < property2.ring_erank {
            return Ok(1);
        }
        if property1.ring_erank > property2.ring_erank {
            return Ok(-1);
        }
        if property1.ring_size < property2.ring_size {
            return Ok(1);
        }
        if property1.ring_size > property2.ring_size {
            return Ok(-1);
        }
        return Ok(0);
    } else if type1 == HETEROCYC {
        return Ok(-1);
    } else if type2 == HETEROCYC {
        return Ok(1);
    }
    if type1 == HETEROAT && type2 == HETEROAT {
        if property1.erank < property2.erank {
            return Ok(1);
        }
        if property1.erank > property2.erank {
            return Ok(-1);
        }
        return Ok(0);
    } else if type1 == HETEROAT {
        return Ok(-1);
    } else if type2 == HETEROAT {
        return Ok(1);
    }
    if type1 == CARBOCYC && type2 == CARBOCYC {
        if property1.ring_size < property2.ring_size {
            return Ok(1);
        }
        if property1.ring_size > property2.ring_size {
            return Ok(-1);
        }
        return Ok(0);
    } else if type1 == CARBOCYC {
        return Ok(-1);
    } else if type2 == CARBOCYC {
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_IsFirstAtomRankLower(
    atom1: i32,
    atom2: i32,
    atom_properties: &[OAD_AtProps],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4369 OAD_Polymer_IsFirstAtomRankLower
    // INCHI✔️✔️: int OAD_Polymer_IsFirstAtomRankLower( int atom1, int atom2, OAD_AtProps *aprops )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /* Compare ext-ranks */
    // INCHI✔️✔️:     int result = OAD_Polymer_CompareRanksOfTwoAtoms( atom1, atom2, aprops );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (result)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return result;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Could not resolve who is junior by extended-ranks...             */
    // INCHI✔️✔️:     /* As a last resort, simply check which canonical number is lesser  */
    // INCHI✔️✔️:     if (atom1 < atom2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (atom1 > atom2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* should not reach there */
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OAD_Polymer_IsFirstAtomRankLower

    let result = OAD_Polymer_CompareRanksOfTwoAtoms(atom1, atom2, atom_properties)?;
    if result != 0 {
        return Ok(result);
    }
    if atom1 < atom2 {
        return Ok(1);
    }
    if atom1 > atom2 {
        return Ok(-1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_DebugTrace(unit: Option<&OAD_PolymerUnit>) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3642 OAD_PolymerUnit_DebugTrace
    // INCHI✔❌: void OAD_PolymerUnit_DebugTrace( OAD_PolymerUnit *u )
    // INCHI✔❌: {
    // INCHI✔❌:     char *conn = "ABSENT", *typ = "ABSENT", *styp = "ABSENT";
    // INCHI✔❌:
    // INCHI✔❌:     if (!u)
    // INCHI✔❌:     {
    // INCHI✔❌:         return;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (u->conn == 1)
    // INCHI✔❌:     {
    // INCHI✔❌:         conn = "HT"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->conn == 2)
    // INCHI✔❌:     {
    // INCHI✔❌:         conn = "HH"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->conn == 3)
    // INCHI✔❌:     {
    // INCHI✔❌:         conn = "EU"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (u->type == 0)
    // INCHI✔❌:     {
    // INCHI✔❌:         typ = "NONE"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->type == 1)
    // INCHI✔❌:     {
    // INCHI✔❌:         typ = "SRU"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->type == 2)
    // INCHI✔❌:     {
    // INCHI✔❌:         typ = "MON"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->type == 3)
    // INCHI✔❌:     {
    // INCHI✔❌:         typ = "COP"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->type == 4)
    // INCHI✔❌:     {
    // INCHI✔❌:         typ = "MOD"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->type == 5)
    // INCHI✔❌:     {
    // INCHI✔❌:         typ = "MER"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (u->subtype == 1)
    // INCHI✔❌:     {
    // INCHI✔❌:         styp = "ALT"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->subtype == 2)
    // INCHI✔❌:     {
    // INCHI✔❌:         styp = "RAN"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->subtype == 3)
    // INCHI✔❌:     {
    // INCHI✔❌:         styp = "BLK"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     {
    // INCHI✔❌:         int i, k;
    // INCHI✔❌:         int na, nb;
    // INCHI✔❌:
    // INCHI✔❌:         ITRACE_("\n\nPOLYMER UNIT @ %-p", u);
    // INCHI✔❌:
    // INCHI✔❌:         ITRACE_( "\n\tid=%-d   label=%-d   type=%-s   subtype=%-s   conn=%-s   subscr='%-s'\n",
    // INCHI✔❌:                 u->id, u->label, typ, styp, conn, u->smt );
    // INCHI✔❌:
    // INCHI✔❌:         ITRACE_( "\tBracket1 coords: %-f, %-f, %-f, %-f\n", u->xbr1[0], u->xbr1[1], u->xbr1[2], u->xbr1[3] );
    // INCHI✔❌:         ITRACE_( "\tBracket2 coords: %-f, %-f, %-f, %-f\n", u->xbr2[0], u->xbr2[1], u->xbr2[2], u->xbr2[3] );
    // INCHI✔❌:
    // INCHI✔❌:         na = u->na;
    // INCHI✔❌:         ITRACE_( "\t%-d atoms { ", na );
    // INCHI✔❌:         for (k = 0; k < na - 1; k++)
    // INCHI✔❌:         {
    // INCHI✔❌:             ITRACE_( " %-d, ", u->alist[k] );
    // INCHI✔❌:         }
    // INCHI✔❌:         ITRACE_( " %-d }\n", u->alist[na - 1] );
    // INCHI✔❌:
    // INCHI✔❌:         nb = u->nb;
    // INCHI✔❌:         ITRACE_( "\t%-d bonds crossing unit borders { ", nb );
    // INCHI✔❌:
    // INCHI✔❌:         for (k = 0; k < nb; k++)
    // INCHI✔❌:         {
    // INCHI✔❌:             ITRACE_( " %-d-%-d ", u->blist[2 * k], u->blist[2 * k + 1] );
    // INCHI✔❌:         }
    // INCHI✔❌:         ITRACE_( "}\n" );
    // INCHI✔❌:
    // INCHI✔❌:         ITRACE_( "\tCRU caps and end atoms { " );
    // INCHI✔❌:
    // INCHI✔❌:         ITRACE_( "*%-d-[-%-d(end1) ... ", u->cap1, u->end_atom1 );
    // INCHI✔❌:         ITRACE_( "%-d(end2)-]-*%-d", u->end_atom2, u->cap2 );
    // INCHI✔❌:         ITRACE_( " }\n" );
    // INCHI✔❌:
    // INCHI✔❌:         ITRACE_( "\tBackbone bonds (may include 'artificially cyclizing' one) : %-d bonds ", u->nbkbonds );
    // INCHI✔❌:         if (u->nbkbonds)
    // INCHI✔❌:         {
    // INCHI✔❌:             ITRACE_(" { ");
    // INCHI✔❌:             for (i = 0; i < u->nbkbonds; i++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ITRACE_( "(%-d, %-d)  ", u->bkbonds[i][0], u->bkbonds[i][1] );
    // INCHI✔❌:             }
    // INCHI✔❌:             ITRACE_("}\n");
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:     return;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_DebugTrace
    // BEGIN INCHI ACTIVE HEADER MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.h:133 ITRACE_
    // INCHI✔❌: #define ITRACE_ 0 && _inchi_trace
    // END INCHI ACTIVE HEADER MACRO: ITRACE_

    let Some(unit) = unit else {
        return;
    };
    let connection = match unit.conn {
        1 => "HT",
        2 => "HH",
        3 => "EU",
        _ => "ABSENT",
    };
    let unit_type = match unit.type_ {
        0 => "NONE",
        1 => "SRU",
        2 => "MON",
        3 => "COP",
        4 => "MOD",
        5 => "MER",
        _ => "ABSENT",
    };
    let subtype = match unit.subtype {
        1 => "ALT",
        2 => "RAN",
        3 => "BLK",
        _ => "ABSENT",
    };
    let _ = (connection, unit_type, subtype);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_port__runichi3__oad_polymerunit_orderbondatomsandbondsthemselves__line_1437() {
        let mut heap = SourceHeap::default();
        let alist = heap.allocate_model_storage(vec![1_i32, 2]).unwrap();
        let stars = heap.allocate_model_storage(vec![9_i32, 8, 7]).unwrap();
        let bonds = heap.allocate_model_storage(vec![1_i32, 9, 2, 8]).unwrap();
        let mut unit = OAD_PolymerUnit {
            na: 2,
            alist,
            nb: 2,
            blist: bonds,
            ..OAD_PolymerUnit::default()
        };
        assert_eq!(
            OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves(
                &mut heap,
                &mut unit,
                3,
                stars.as_const(),
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(bonds.as_const()).unwrap(), &[8, 2, 9, 1]);

        let partial_bonds = heap.allocate_model_storage(vec![1_i32, 9, 8, 7]).unwrap();
        let mut partial = OAD_PolymerUnit {
            na: 2,
            alist,
            nb: 2,
            blist: partial_bonds,
            ..OAD_PolymerUnit::default()
        };
        assert_eq!(
            OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves(
                &mut heap,
                &mut partial,
                3,
                stars.as_const(),
            ),
            Ok(1)
        );
        assert_eq!(heap.slice(partial_bonds.as_const()).unwrap(), &[9, 1, 8, 7]);

        let external = heap.allocate_model_storage(vec![6_i32, 1]).unwrap();
        let mut single = OAD_PolymerUnit {
            na: 2,
            alist,
            nb: 1,
            blist: external,
            ..OAD_PolymerUnit::default()
        };
        assert_eq!(
            OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves(
                &mut heap,
                &mut single,
                0,
                SourceConstPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(external.as_const()).unwrap(), &[6, 1]);

        let mut negative = OAD_PolymerUnit {
            na: i32::MIN,
            nb: i32::MIN,
            alist: SourceMutPointer::null(),
            blist: SourceMutPointer::null(),
            ..OAD_PolymerUnit::default()
        };
        assert_eq!(
            OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves(
                &mut heap,
                &mut negative,
                i32::MIN,
                SourceConstPointer::null(),
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymer_prepareworkingset__line_2342() {
        fn allocate_unit(
            heap: &mut SourceHeap,
            atom_list: Vec<i32>,
            bond_list: Vec<i32>,
            fields: [i32; 4],
            backbone_bonds: &[[i32; 2]],
        ) -> (
            SourceMutPointer<OAD_PolymerUnit>,
            Vec<SourceMutPointer<i32>>,
        ) {
            let na = i32::try_from(atom_list.len()).unwrap();
            let nb = i32::try_from(bond_list.len() / 2).unwrap();
            let alist = if atom_list.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(atom_list).unwrap()
            };
            let blist = if bond_list.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(bond_list).unwrap()
            };
            let rows = backbone_bonds
                .iter()
                .map(|row| heap.allocate_model_storage(row.to_vec()).unwrap())
                .collect::<Vec<_>>();
            let bkbonds = if rows.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(rows.clone()).unwrap()
            };
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    na,
                    nb,
                    cap1: fields[0],
                    cap2: fields[1],
                    end_atom1: fields[2],
                    end_atom2: fields[3],
                    alist,
                    blist,
                    nbkbonds: i32::try_from(rows.len()).unwrap(),
                    maxbkbonds: i32::try_from(rows.len()).unwrap(),
                    bkbonds,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            (unit, rows)
        }

        let mut heap = SourceHeap::default();
        let canonical = heap
            .allocate_model_storage(vec![4_i32, -1, 2, 0, 5, 1, 3, 6])
            .unwrap();
        let pzz = heap.allocate_model_storage(vec![6_i32]).unwrap();
        let (first, first_rows) = allocate_unit(
            &mut heap,
            vec![3, 1, 0, 2],
            vec![2, 4, 3, 0],
            [0, 2, 3, 5],
            &[[4, 0], [1, 2], [0, 1]],
        );
        let (second, _) = allocate_unit(&mut heap, vec![7, 5], vec![], [0, 0, 0, 0], &[]);
        let units = heap.allocate_model_storage(vec![first, second]).unwrap();
        let unit_numbers = heap.allocate_model_storage(vec![-9_i32, -9]).unwrap();
        let mut polymer = OAD_Polymer {
            n: 2,
            n_pzz: 1,
            pzz,
            ..OAD_Polymer::default()
        };
        assert_eq!(
            OAD_Polymer_PrepareWorkingSet(
                &mut heap,
                &mut polymer,
                canonical.as_const(),
                SourceConstPointer::null(),
                units,
                unit_numbers,
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(pzz.as_const()).unwrap(), &[4]);
        assert_eq!(heap.slice(unit_numbers.as_const()).unwrap(), &[1, 0]);
        let first_value = &heap.slice(first.as_const()).unwrap()[0];
        assert_eq!(first_value.na, 3);
        assert_eq!(
            (
                first_value.cap1,
                first_value.cap2,
                first_value.end_atom1,
                first_value.end_atom2
            ),
            (5, 3, 1, 2)
        );
        assert_eq!(
            heap.slice(first_value.alist.as_const()).unwrap(),
            &[1, 3, 5, 2]
        );
        assert_eq!(
            heap.slice(first_value.blist.as_const()).unwrap(),
            &[1, 5, 6, 3]
        );
        assert_eq!(heap.slice(first_rows[0].as_const()).unwrap(), &[5, 6]);
        assert_eq!(heap.slice(first_rows[1].as_const()).unwrap(), &[1, 2]);
        assert_eq!(heap.slice(first_rows[2].as_const()).unwrap(), &[0, 1]);
        let second_value = &heap.slice(second.as_const()).unwrap()[0];
        assert_eq!(heap.slice(second_value.alist.as_const()).unwrap(), &[2, 7]);

        let mut error10_heap = SourceHeap::default();
        let error10_canonical = error10_heap
            .allocate_model_storage(vec![0_i32, -1])
            .unwrap();
        let error10_pzz = error10_heap.allocate_model_storage(vec![0_i32, 1]).unwrap();
        let mut error10_polymer = OAD_Polymer {
            n_pzz: 2,
            pzz: error10_pzz,
            ..OAD_Polymer::default()
        };
        assert_eq!(
            OAD_Polymer_PrepareWorkingSet(
                &mut error10_heap,
                &mut error10_polymer,
                error10_canonical.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(10)
        );
        assert_eq!(error10_heap.slice(error10_pzz.as_const()).unwrap(), &[1, 1]);

        let mut bond_error_heap = SourceHeap::default();
        let bond_error_canonical = bond_error_heap
            .allocate_model_storage(vec![4_i32, 3, -1])
            .unwrap();
        let (bond_error_unit, _) =
            allocate_unit(&mut bond_error_heap, vec![0], vec![1, 2], [0, 0, 0, 0], &[]);
        let bond_error_units = bond_error_heap
            .allocate_model_storage(vec![bond_error_unit])
            .unwrap();
        let mut bond_error_polymer = OAD_Polymer {
            n: 1,
            ..OAD_Polymer::default()
        };
        assert_eq!(
            OAD_Polymer_PrepareWorkingSet(
                &mut bond_error_heap,
                &mut bond_error_polymer,
                bond_error_canonical.as_const(),
                SourceConstPointer::null(),
                bond_error_units,
                SourceMutPointer::null(),
            ),
            Ok(11)
        );
        let bond_error_value = &bond_error_heap.slice(bond_error_unit.as_const()).unwrap()[0];
        assert_eq!(bond_error_value.na, 1);
        assert_eq!(
            bond_error_heap
                .slice(bond_error_value.blist.as_const())
                .unwrap(),
            &[4, 2]
        );
        assert_eq!(bond_error_value.cap1, 0);

        for failed_field in 0..4 {
            let mut error_heap = SourceHeap::default();
            let mut canonical_values = vec![9_i32, 8, 7, 6, 5, 4];
            canonical_values[failed_field + 1] = -1;
            let canonical = error_heap.allocate_model_storage(canonical_values).unwrap();
            let (first_unit, _) =
                allocate_unit(&mut error_heap, vec![5, 0], vec![], [0, 0, 0, 0], &[]);
            let (failing_unit, _) =
                allocate_unit(&mut error_heap, vec![0], vec![], [1, 2, 3, 4], &[]);
            let units = error_heap
                .allocate_model_storage(vec![first_unit, failing_unit])
                .unwrap();
            let mut polymer = OAD_Polymer {
                n: 2,
                ..OAD_Polymer::default()
            };
            assert_eq!(
                OAD_Polymer_PrepareWorkingSet(
                    &mut error_heap,
                    &mut polymer,
                    canonical.as_const(),
                    SourceConstPointer::null(),
                    units,
                    SourceMutPointer::null(),
                ),
                Ok(11)
            );
            let first_value = &error_heap.slice(first_unit.as_const()).unwrap()[0];
            assert_eq!(
                error_heap.slice(first_value.alist.as_const()).unwrap(),
                &[5, 10]
            );
            let value = &error_heap.slice(failing_unit.as_const()).unwrap()[0];
            let actual = [value.cap1, value.cap2, value.end_atom1, value.end_atom2];
            for index in 0..4 {
                let expected = if index < failed_field {
                    9 - i32::try_from(index).unwrap()
                } else {
                    i32::try_from(index + 1).unwrap()
                };
                assert_eq!(
                    actual[index], expected,
                    "failed field {failed_field}, field {index}"
                );
            }
        }

        let mut error12_heap = SourceHeap::default();
        let error12_canonical = error12_heap
            .allocate_model_storage(vec![0_i32, 1, 2])
            .unwrap();
        let (error12_unit, _) =
            allocate_unit(&mut error12_heap, vec![0], vec![1, 2], [0, 0, 0, 0], &[]);
        let error12_units = error12_heap
            .allocate_model_storage(vec![error12_unit])
            .unwrap();
        let error12_numbers = error12_heap.allocate_model_storage(vec![-7_i32]).unwrap();
        let mut error12_polymer = OAD_Polymer {
            n: 1,
            ..OAD_Polymer::default()
        };
        assert_eq!(
            OAD_Polymer_PrepareWorkingSet(
                &mut error12_heap,
                &mut error12_polymer,
                error12_canonical.as_const(),
                SourceConstPointer::null(),
                error12_units,
                error12_numbers,
            ),
            Ok(12)
        );
        assert_eq!(
            error12_heap.slice(error12_numbers.as_const()).unwrap(),
            &[-7]
        );

        let mut empty_heap = SourceHeap::default();
        let mut empty_polymer = OAD_Polymer {
            n: i32::MIN,
            n_pzz: i32::MIN,
            pzz: SourceMutPointer::null(),
            ..OAD_Polymer::default()
        };
        assert_eq!(
            OAD_Polymer_PrepareWorkingSet(
                &mut empty_heap,
                &mut empty_polymer,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__runichi3__origatdata_removehalfbond__line_2517() {
        let mut atom = inp_ATOM {
            valence: 3,
            neighbor: {
                let mut values = [0; 20];
                values[..3].copy_from_slice(&[4, 7, 9]);
                values
            },
            bond_type: {
                let mut values = [0; 20];
                values[..3].copy_from_slice(&[1, 2, 3]);
                values
            },
            bond_stereo: {
                let mut values = [0; 20];
                values[..3].copy_from_slice(&[-1, 0, 1]);
                values
            },
            ..inp_ATOM::default()
        };
        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(vec![atom.clone()]).unwrap();
        let mut bond_type = -8;
        let mut bond_stereo = -9;
        assert_eq!(
            OrigAtData_RemoveHalfBond(&mut heap, 0, 7, atoms, &mut bond_type, &mut bond_stereo,),
            Ok(1)
        );
        assert_eq!((bond_type, bond_stereo), (2, 0));
        atom = heap.slice(atoms.as_const()).unwrap()[0].clone();
        assert_eq!(&atom.neighbor[..4], &[4, 9, 0, 0]);
        assert_eq!(&atom.bond_type[..4], &[1, 3, 0, 0]);
        assert_eq!(&atom.bond_stereo[..4], &[-1, 1, 0, 0]);

        let before = atom.clone();
        assert_eq!(
            OrigAtData_RemoveHalfBond(&mut heap, 0, 99, atoms, &mut bond_type, &mut bond_stereo,),
            Ok(0)
        );
        assert_eq!(heap.slice(atoms.as_const()).unwrap()[0], before);
        assert_eq!((bond_type, bond_stereo), (2, 0));
        assert_eq!(
            OrigAtData_RemoveHalfBond(&mut heap, -1, 4, atoms, &mut bond_type, &mut bond_stereo,),
            Ok(0)
        );
        assert_eq!(
            OrigAtData_RemoveHalfBond(
                &mut heap,
                0,
                4,
                SourceMutPointer::null(),
                &mut bond_type,
                &mut bond_stereo,
            ),
            Ok(0)
        );

        let mut full = inp_ATOM {
            valence: 20,
            ..inp_ATOM::default()
        };
        for index in 0..20 {
            full.neighbor[index] = u16::try_from(index + 1).unwrap();
            full.bond_type[index] = u8::try_from(index + 1).unwrap();
            full.bond_stereo[index] = i8::try_from(index as i32 - 10).unwrap();
        }
        let full_atoms = heap.allocate_model_storage(vec![full]).unwrap();
        assert_eq!(
            OrigAtData_RemoveHalfBond(
                &mut heap,
                0,
                1,
                full_atoms,
                &mut bond_type,
                &mut bond_stereo,
            ),
            Ok(1)
        );
        assert_eq!((bond_type, bond_stereo), (1, -10));
        let full = &heap.slice(full_atoms.as_const()).unwrap()[0];
        assert_eq!(full.neighbor[0], 2);
        assert_eq!(full.bond_type[0], 2);
        assert_eq!(full.bond_stereo[0], -9);
        assert_eq!(full.neighbor[19], 0);
        assert_eq!(full.bond_type[19], 0);
        assert_eq!(full.bond_stereo[19], 0);
    }

    #[test]
    fn source_port__runichi3__origatdata_removebond__line_2577() {
        fn atom(neighbors: &[u16], types: &[u8], stereos: &[i8], chem: i8) -> inp_ATOM {
            let mut atom = inp_ATOM {
                valence: i8::try_from(neighbors.len()).unwrap(),
                chem_bonds_valence: chem,
                ..inp_ATOM::default()
            };
            atom.neighbor[..neighbors.len()].copy_from_slice(neighbors);
            atom.bond_type[..types.len()].copy_from_slice(types);
            atom.bond_stereo[..stereos.len()].copy_from_slice(stereos);
            atom
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![
                atom(&[1, 2], &[2, 1], &[-1, 0], 3),
                atom(&[0, 2], &[2, 3], &[1, 0], 5),
            ])
            .unwrap();
        let mut bond_type = -1;
        let mut bond_stereo = -2;
        let mut num_bonds = 2;
        assert_eq!(
            OrigAtData_RemoveBond(
                &mut heap,
                0,
                1,
                atoms,
                &mut bond_type,
                &mut bond_stereo,
                &mut num_bonds,
            ),
            Ok(1)
        );
        assert_eq!((bond_type, bond_stereo, num_bonds), (2, 1, 1));
        let values = heap.slice(atoms.as_const()).unwrap();
        assert_eq!((values[0].valence, values[0].chem_bonds_valence), (1, 1));
        assert_eq!((values[1].valence, values[1].chem_bonds_valence), (1, 3));
        assert_eq!(&values[0].neighbor[..3], &[2, 0, 0]);
        assert_eq!(&values[1].neighbor[..3], &[2, 0, 0]);

        let partial = heap
            .allocate_model_storage(vec![atom(&[1], &[3], &[-1], 3), atom(&[2], &[4], &[1], 4)])
            .unwrap();
        bond_type = -1;
        bond_stereo = -2;
        num_bonds = 7;
        assert_eq!(
            OrigAtData_RemoveBond(
                &mut heap,
                0,
                1,
                partial,
                &mut bond_type,
                &mut bond_stereo,
                &mut num_bonds,
            ),
            Ok(0)
        );
        assert_eq!((bond_type, bond_stereo, num_bonds), (3, -1, 7));
        let values = heap.slice(partial.as_const()).unwrap();
        assert_eq!(values[0].neighbor[0], 0);
        assert_eq!((values[0].valence, values[0].chem_bonds_valence), (1, 3));
        assert_eq!(values[1].neighbor[0], 2);

        let before = heap.slice(partial.as_const()).unwrap().to_vec();
        assert_eq!(
            OrigAtData_RemoveBond(
                &mut heap,
                -1,
                1,
                partial,
                &mut bond_type,
                &mut bond_stereo,
                &mut num_bonds,
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(partial.as_const()).unwrap(), before.as_slice());
        assert_eq!(
            OrigAtData_RemoveBond(
                &mut heap,
                0,
                1,
                SourceMutPointer::null(),
                &mut bond_type,
                &mut bond_stereo,
                &mut num_bonds,
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__runichi3__origatdata_addbond__line_2607() {
        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default(), inp_ATOM::default()])
            .unwrap();
        let mut num_bonds = 4;
        assert_eq!(
            OrigAtData_AddBond(&mut heap, 0, 1, atoms, 99, 200, &mut num_bonds),
            Ok(1)
        );
        assert_eq!(num_bonds, 5);
        let values = heap.slice(atoms.as_const()).unwrap();
        for (atom, neighbor) in [(&values[0], 1_u16), (&values[1], 0_u16)] {
            assert_eq!(atom.valence, 1);
            assert_eq!(atom.chem_bonds_valence, 1);
            assert_eq!(atom.neighbor[0], neighbor);
            assert_eq!(atom.bond_type[0], 1);
            assert_eq!(atom.bond_stereo[0], -56);
        }

        let before = heap.slice(atoms.as_const()).unwrap().to_vec();
        assert_eq!(
            OrigAtData_AddBond(&mut heap, 0, 1, atoms, 3, -7, &mut num_bonds),
            Ok(1)
        );
        assert_eq!(num_bonds, 6);
        assert_eq!(heap.slice(atoms.as_const()).unwrap(), before.as_slice());

        let asymmetric = heap
            .allocate_model_storage(vec![
                {
                    let mut atom = inp_ATOM {
                        valence: 1,
                        chem_bonds_valence: 2,
                        ..inp_ATOM::default()
                    };
                    atom.neighbor[0] = 1;
                    atom.bond_type[0] = 2;
                    atom
                },
                inp_ATOM::default(),
            ])
            .unwrap();
        num_bonds = 0;
        assert_eq!(
            OrigAtData_AddBond(&mut heap, 0, 1, asymmetric, 2, -3, &mut num_bonds),
            Ok(1)
        );
        let values = heap.slice(asymmetric.as_const()).unwrap();
        assert_eq!(values[0].valence, 1);
        assert_eq!(values[1].valence, 1);
        assert_eq!(values[1].bond_type[0], 2);
        assert_eq!(values[1].bond_stereo[0], -3);

        let mut full = inp_ATOM {
            valence: 20,
            chem_bonds_valence: 77,
            ..inp_ATOM::default()
        };
        full.neighbor = [8; 20];
        let full_atoms = heap
            .allocate_model_storage(vec![full, inp_ATOM::default()])
            .unwrap();
        let before = heap.slice(full_atoms.as_const()).unwrap().to_vec();
        num_bonds = 10;
        assert_eq!(
            OrigAtData_AddBond(&mut heap, 0, 1, full_atoms, 2, 3, &mut num_bonds),
            Ok(0)
        );
        assert_eq!(num_bonds, 10);
        assert_eq!(
            heap.slice(full_atoms.as_const()).unwrap(),
            before.as_slice()
        );

        let self_atom = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        num_bonds = -1;
        assert_eq!(
            OrigAtData_AddBond(&mut heap, 0, 0, self_atom, 3, 1, &mut num_bonds),
            Ok(1)
        );
        let value = &heap.slice(self_atom.as_const()).unwrap()[0];
        assert_eq!((value.valence, value.chem_bonds_valence), (1, 3));
        assert_eq!((value.neighbor[0], value.bond_type[0]), (0, 3));
        assert_eq!(num_bonds, 0);

        assert_eq!(
            OrigAtData_AddBond(
                &mut heap,
                0,
                1,
                SourceMutPointer::null(),
                2,
                0,
                &mut num_bonds,
            ),
            Ok(0)
        );
        assert_eq!(num_bonds, 0);
    }

    #[test]
    fn source_port__runichi3__unmarkringsystemsinp__line_1961() {
        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![
                inp_ATOM {
                    bCutVertex: -7,
                    nRingSystem: 11,
                    nNumAtInRingSystem: 12,
                    nBlockSystem: 13,
                    charge: -2,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    bCutVertex: 8,
                    nRingSystem: u16::MAX,
                    nNumAtInRingSystem: u16::MAX,
                    nBlockSystem: u16::MAX,
                    charge: 3,
                    ..inp_ATOM::default()
                },
            ])
            .unwrap();
        assert_eq!(UnMarkRingSystemsInp(&mut heap, atoms, 2), Ok(0));
        let values = heap.slice(atoms.as_const()).unwrap();
        for atom in values {
            assert_eq!(
                (
                    atom.bCutVertex,
                    atom.nRingSystem,
                    atom.nNumAtInRingSystem,
                    atom.nBlockSystem,
                ),
                (0, 0, 0, 0)
            );
        }
        assert_eq!((values[0].charge, values[1].charge), (-2, 3));
        assert_eq!(
            UnMarkRingSystemsInp(&mut heap, SourceMutPointer::null(), 0),
            Ok(0)
        );
        assert_eq!(
            UnMarkRingSystemsInp(&mut heap, SourceMutPointer::null(), i32::MIN),
            Ok(0)
        );
    }

    #[test]
    fn source_port__runichi3__origatdata_addsinglestereolessbond__line_2682() {
        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default(), inp_ATOM::default()])
            .unwrap();
        let mut num_bonds = 0;
        assert_eq!(
            OrigAtData_AddSingleStereolessBond(&mut heap, 0, 1, atoms, &mut num_bonds,),
            Ok(1)
        );
        assert_eq!(num_bonds, 1);
        let values = heap.slice(atoms.as_const()).unwrap();
        assert_eq!((values[0].bond_type[0], values[0].bond_stereo[0]), (1, 0));
        assert_eq!((values[1].bond_type[0], values[1].bond_stereo[0]), (1, 0));
        assert_eq!(
            OrigAtData_AddSingleStereolessBond(
                &mut heap,
                0,
                1,
                SourceMutPointer::null(),
                &mut num_bonds,
            ),
            Ok(0)
        );
        assert_eq!(num_bonds, 1);
    }

    #[test]
    fn source_port__runichi3__oad_polymer_findringsystems__line_3236() {
        let mut heap = SourceHeap::default();
        let mut triangle_atoms = vec![inp_ATOM::default(); 3];
        triangle_atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
        triangle_atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
        triangle_atoms[2].neighbor[..2].copy_from_slice(&[0, 1]);
        for (index, atom) in triangle_atoms.iter_mut().enumerate() {
            atom.valence = 2;
            atom.orig_at_number = u16::try_from(index).unwrap();
        }
        let atoms = heap.allocate_model_storage(triangle_atoms).unwrap();
        let polymer = heap
            .allocate_model_storage(vec![OAD_Polymer::default()])
            .unwrap();
        let ring_numbers = heap.allocate_model_storage(vec![-8_i32; 4]).unwrap();
        let ring_sizes = heap.allocate_model_storage(vec![-9_i32; 4]).unwrap();
        let mut num_bonds = 3;
        assert_eq!(
            OAD_Polymer_FindRingSystems(
                &mut heap,
                polymer.as_const(),
                atoms,
                3,
                &mut num_bonds,
                ring_numbers,
                ring_sizes,
                0,
            ),
            Ok(3)
        );
        assert_eq!(num_bonds, 3);
        assert_eq!(heap.slice(ring_numbers.as_const()).unwrap(), &[1, 1, 1, -1]);
        assert_eq!(heap.slice(ring_sizes.as_const()).unwrap(), &[3, 3, 3, -9]);
        assert!(
            heap.slice(atoms.as_const())
                .unwrap()
                .iter()
                .all(|atom| atom.nRingSystem == 0 && atom.nBlockSystem == 0)
        );

        assert_eq!(
            OAD_Polymer_FindRingSystems(
                &mut heap,
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                3,
                &mut num_bonds,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                0,
            ),
            Ok(0)
        );

        let mut path_atoms = vec![inp_ATOM::default(); 3];
        path_atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
        path_atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
        path_atoms[2].neighbor[..2].copy_from_slice(&[0, 1]);
        for atom in &mut path_atoms {
            atom.valence = 2;
            atom.chem_bonds_valence = 2;
            atom.bond_type[..2].copy_from_slice(&[1, 1]);
        }
        let path_atoms = heap.allocate_model_storage(path_atoms).unwrap();
        let unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                cyclized: 1,
                end_atom1: 1,
                end_atom2: 3,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let units = heap.allocate_model_storage(vec![unit]).unwrap();
        let polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                n: 1,
                units,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let numbers = heap.allocate_model_storage(vec![-4_i32; 4]).unwrap();
        let mut num_bonds = 3;
        assert_eq!(
            OAD_Polymer_FindRingSystems(
                &mut heap,
                polymer.as_const(),
                path_atoms,
                3,
                &mut num_bonds,
                numbers,
                SourceMutPointer::null(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(num_bonds, 3);
        assert_eq!(heap.slice(numbers.as_const()).unwrap(), &[-1, -1, -1, -1]);
    }

    #[test]
    fn source_port__runichi3__oad_polymer_setatprops__line_3312() {
        let mut heap = SourceHeap::default();
        let mut atoms_values = vec![inp_ATOM::default(); 3];
        atoms_values[0].neighbor[..2].copy_from_slice(&[1, 2]);
        atoms_values[1].neighbor[..2].copy_from_slice(&[0, 2]);
        atoms_values[2].neighbor[..2].copy_from_slice(&[0, 1]);
        for (index, atom) in atoms_values.iter_mut().enumerate() {
            atom.valence = 2;
            atom.orig_at_number = u16::try_from(index).unwrap();
            atom.el_number = [6, 7, 8][index];
        }
        let atoms = heap.allocate_model_storage(atoms_values).unwrap();
        let polymer = heap
            .allocate_model_storage(vec![OAD_Polymer::default()])
            .unwrap();
        let properties = heap
            .allocate_model_storage(vec![OAD_AtProps::default(); 3])
            .unwrap();
        let mut num_bonds = 3;
        assert_eq!(
            OAD_Polymer_SetAtProps(
                &mut heap,
                polymer.as_const(),
                atoms,
                3,
                &mut num_bonds,
                properties,
                SourceConstPointer::null(),
            ),
            Ok(())
        );
        let values = heap.slice(properties.as_const()).unwrap();
        assert_eq!(
            values,
            &[
                OAD_AtProps {
                    erank: 2,
                    ring_erank: 216,
                    ring_num: 1,
                    ring_size: 3,
                },
                OAD_AtProps {
                    erank: 211,
                    ring_erank: 216,
                    ring_num: 1,
                    ring_size: 3,
                },
                OAD_AtProps {
                    erank: 215,
                    ring_erank: 216,
                    ring_num: 1,
                    ring_size: 3,
                },
            ]
        );
        assert_eq!(num_bonds, 3);

        let canonical = heap.allocate_model_storage(vec![2_i32, 0, 1]).unwrap();
        let mapped_properties = heap
            .allocate_model_storage(vec![OAD_AtProps::default(); 4])
            .unwrap();
        assert_eq!(
            OAD_Polymer_SetAtProps(
                &mut heap,
                polymer.as_const(),
                atoms,
                3,
                &mut num_bonds,
                mapped_properties,
                canonical.as_const(),
            ),
            Ok(())
        );
        let values = heap.slice(mapped_properties.as_const()).unwrap();
        assert_eq!(
            (values[1].ring_num, values[2].ring_num, values[3].ring_num),
            (1, 1, 1)
        );
        assert_eq!(
            (
                values[1].ring_size,
                values[2].ring_size,
                values[3].ring_size
            ),
            (3, 3, 3)
        );

        let failure_properties = heap
            .allocate_model_storage(vec![OAD_AtProps::default(); 3])
            .unwrap();
        heap.trace_source_allocations();
        heap.fail_after_allocations(0);
        assert_eq!(
            OAD_Polymer_SetAtProps(
                &mut heap,
                polymer.as_const(),
                atoms,
                3,
                &mut num_bonds,
                failure_properties,
                SourceConstPointer::null(),
            ),
            Ok(())
        );
        assert_eq!(
            heap.slice(failure_properties.as_const()).unwrap()[0].erank,
            2
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_compareatomlistsmod__line_1377() {
        let mut heap = SourceHeap::default();
        let first_values = heap.allocate_model_storage(vec![1_i32, 4]).unwrap();
        let second_values = heap.allocate_model_storage(vec![1_i32, 5]).unwrap();
        let first = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: 2,
                alist: first_values,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let second = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: 2,
                alist: second_values,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        assert_eq!(
            OAD_PolymerUnit_CompareAtomListsMod(&heap, first.as_const(), second.as_const()),
            Ok(-1)
        );
        assert_eq!(
            OAD_PolymerUnit_CompareAtomListsMod(&heap, second.as_const(), first.as_const()),
            Ok(1)
        );

        let same_values = heap.allocate_model_storage(vec![1_i32, 4]).unwrap();
        let same = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: 2,
                alist: same_values,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        assert_eq!(
            OAD_PolymerUnit_CompareAtomListsMod(&heap, first.as_const(), same.as_const()),
            Ok(0)
        );

        let short = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: 1,
                alist: first_values,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        assert_eq!(
            OAD_PolymerUnit_CompareAtomListsMod(&heap, short.as_const(), first.as_const()),
            Ok(-1)
        );
        assert_eq!(
            OAD_PolymerUnit_CompareAtomListsMod(&heap, first.as_const(), short.as_const()),
            Ok(1)
        );

        let negative = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: -1,
                alist: SourceMutPointer::null(),
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let negative_other = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: -1,
                alist: SourceMutPointer::null(),
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        assert_eq!(
            OAD_PolymerUnit_CompareAtomListsMod(
                &heap,
                negative.as_const(),
                negative_other.as_const()
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_createcopy__line_1260() {
        fn source_unit(heap: &mut SourceHeap) -> SourceMutPointer<OAD_PolymerUnit> {
            let alist = heap.allocate_model_storage(vec![4_i32, 5]).unwrap();
            let blist = heap.allocate_model_storage(vec![1_i32, 2, 3, 4]).unwrap();
            let first_row = heap.allocate_model_storage(vec![6_i32, 7]).unwrap();
            let second_row = heap.allocate_model_storage(vec![8_i32, 9]).unwrap();
            let bkbonds = heap
                .allocate_model_storage(vec![first_row, second_row])
                .unwrap();
            let mut smt = [99_i8; 80];
            smt[..4].copy_from_slice(&[b'A' as i8, b'B' as i8, b'C' as i8, 0]);
            heap.allocate_model_storage(vec![OAD_PolymerUnit {
                id: 11,
                type_: 12,
                subtype: 13,
                conn: 14,
                label: 15,
                na: 2,
                nb: 2,
                cyclizable: 16,
                cyclized: 17,
                xbr1: [1.0, -0.0, f64::INFINITY, f64::NAN],
                xbr2: [-1.0, 0.0, f64::NEG_INFINITY, f64::MIN],
                smt,
                representation: 91,
                cap1: 18,
                end_atom1: 19,
                cap2: 20,
                end_atom2: 21,
                cap1_is_undef: 22,
                cap2_is_undef: 23,
                alist,
                blist,
                maxbkbonds: 1,
                nbkbonds: 2,
                bkbonds,
            }])
            .unwrap()
        }

        let mut heap = SourceHeap::default();
        let source = source_unit(&mut heap);
        let copy = OAD_PolymerUnit_CreateCopy(&mut heap, source).unwrap();
        assert!(!copy.is_null());
        let source_value = heap.slice(source.as_const()).unwrap()[0].clone();
        let copy_value = heap.slice(copy.as_const()).unwrap()[0].clone();
        assert_eq!(copy_value.id, source_value.id);
        assert_eq!(copy_value.type_, source_value.type_);
        assert_eq!(copy_value.subtype, source_value.subtype);
        assert_eq!(copy_value.conn, source_value.conn);
        assert_eq!(copy_value.label, source_value.label);
        assert_eq!(copy_value.na, source_value.na);
        assert_eq!(copy_value.nb, source_value.nb);
        assert_eq!(copy_value.cyclizable, source_value.cyclizable);
        assert_eq!(copy_value.cyclized, source_value.cyclized);
        for (actual, expected) in copy_value.xbr1.iter().zip(source_value.xbr1.iter()) {
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
        for (actual, expected) in copy_value.xbr2.iter().zip(source_value.xbr2.iter()) {
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
        assert_eq!(
            &copy_value.smt[..4],
            &[b'A' as i8, b'B' as i8, b'C' as i8, 0]
        );
        assert!(copy_value.smt[4..].iter().all(|byte| *byte == 0));
        assert_eq!(copy_value.representation, 0);
        assert_eq!(copy_value.maxbkbonds, 2);
        assert_ne!(copy_value.alist, source_value.alist);
        assert_ne!(copy_value.blist, source_value.blist);
        assert_ne!(copy_value.bkbonds, source_value.bkbonds);
        assert_eq!(heap.slice(copy_value.alist.as_const()).unwrap(), &[4, 5]);
        assert_eq!(
            heap.slice(copy_value.blist.as_const()).unwrap(),
            &[1, 2, 3, 4]
        );
        let copy_rows = heap.slice(copy_value.bkbonds.as_const()).unwrap();
        let source_rows = heap.slice(source_value.bkbonds.as_const()).unwrap();
        assert_ne!(copy_rows[0], source_rows[0]);
        assert_ne!(copy_rows[1], source_rows[1]);
        assert_eq!(heap.slice(copy_rows[0].as_const()).unwrap(), &[6, 7]);
        assert_eq!(heap.slice(copy_rows[1].as_const()).unwrap(), &[8, 9]);
        OAD_PolymerUnit_Free(&mut heap, copy).unwrap();
        assert!(heap.slice(source.as_const()).is_ok());
        OAD_PolymerUnit_Free(&mut heap, source).unwrap();

        for successful_allocations in 0..=5 {
            let mut failure_heap = SourceHeap::default();
            let source = source_unit(&mut failure_heap);
            let source_value = failure_heap.slice(source.as_const()).unwrap()[0].clone();
            failure_heap.fail_after_allocations(successful_allocations);
            assert!(
                OAD_PolymerUnit_CreateCopy(&mut failure_heap, source)
                    .unwrap()
                    .is_null(),
                "successful_allocations={successful_allocations}"
            );
            assert!(failure_heap.slice(source.as_const()).is_ok());
            assert_eq!(
                failure_heap.slice(source_value.alist.as_const()).unwrap(),
                &[4, 5]
            );
            OAD_PolymerUnit_Free(&mut failure_heap, source).unwrap();
        }
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_free__line_1342() {
        let mut heap = SourceHeap::default();

        OAD_PolymerUnit_Free(&mut heap, SourceMutPointer::null()).unwrap();

        let empty_unit = heap.allocate(vec![OAD_PolymerUnit::default()]).unwrap();
        OAD_PolymerUnit_Free(&mut heap, empty_unit).unwrap();
        assert_eq!(
            heap.slice(empty_unit.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let alist = heap.allocate(vec![1_i32, 2, 3]).unwrap();
        let blist = heap.allocate(vec![1_i32, 2, 2, 3]).unwrap();
        let row_zero = heap.allocate(vec![1_i32, 2]).unwrap();
        let row_one = heap.allocate(vec![2_i32, 3]).unwrap();
        let bkbonds = heap.allocate(vec![row_zero, row_one]).unwrap();
        let complete_unit = heap
            .allocate(vec![OAD_PolymerUnit {
                alist,
                blist,
                maxbkbonds: 2,
                bkbonds,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        OAD_PolymerUnit_Free(&mut heap, complete_unit).unwrap();
        for result in [
            heap.slice(alist.as_const()),
            heap.slice(blist.as_const()),
            heap.slice(row_zero.as_const()),
            heap.slice(row_one.as_const()),
        ] {
            assert_eq!(result, Err(SourceHeapError::MissingAllocation));
        }
        assert_eq!(
            heap.slice(bkbonds.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(complete_unit.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let only_blist = heap.allocate(vec![4_i32, 5]).unwrap();
        let mixed_unit = heap
            .allocate(vec![OAD_PolymerUnit {
                na: i32::MAX,
                nb: i32::MIN,
                nbkbonds: i32::MAX,
                alist: SourceMutPointer::null(),
                blist: only_blist,
                bkbonds: SourceMutPointer::null(),
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        OAD_PolymerUnit_Free(&mut heap, mixed_unit).unwrap();
        assert_eq!(
            heap.slice(only_blist.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let unvisited_row = heap.allocate(vec![6_i32, 7]).unwrap();
        let negative_bkbonds = heap.allocate(vec![unvisited_row]).unwrap();
        let negative_count_unit = heap
            .allocate(vec![OAD_PolymerUnit {
                maxbkbonds: i32::MIN,
                bkbonds: negative_bkbonds,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        OAD_PolymerUnit_Free(&mut heap, negative_count_unit).unwrap();
        assert_eq!(heap.slice(unvisited_row.as_const()).unwrap(), &[6, 7]);
        assert_eq!(
            heap.slice(negative_bkbonds.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        heap.free(unvisited_row).unwrap();
    }

    #[test]
    fn source_port__runichi3__oad_polymer_free__line_2071() {
        let mut heap = SourceHeap::default();

        OAD_Polymer_Free(&mut heap, SourceMutPointer::null()).unwrap();

        let empty_polymer = heap.allocate(vec![OAD_Polymer::default()]).unwrap();
        OAD_Polymer_Free(&mut heap, empty_polymer).unwrap();
        assert_eq!(
            heap.slice(empty_polymer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let first_alist = heap.allocate(vec![1_i32]).unwrap();
        let first_unit = heap
            .allocate(vec![OAD_PolymerUnit {
                alist: first_alist,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let second_blist = heap.allocate(vec![1_i32, 2]).unwrap();
        let second_unit = heap
            .allocate(vec![OAD_PolymerUnit {
                blist: second_blist,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let units = heap.allocate(vec![first_unit, second_unit]).unwrap();
        let pzz = heap.allocate(vec![7_i32, 8]).unwrap();
        let complete_polymer = heap
            .allocate(vec![OAD_Polymer {
                units,
                n: 2,
                n_pzz: 2,
                pzz,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        OAD_Polymer_Free(&mut heap, complete_polymer).unwrap();
        assert_eq!(
            heap.slice(first_alist.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(second_blist.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(first_unit.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(second_unit.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(units.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(pzz.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(complete_polymer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let retained_zero_unit = heap.allocate(vec![OAD_PolymerUnit::default()]).unwrap();
        let retained_zero_units = heap.allocate(vec![retained_zero_unit]).unwrap();
        let zero_count_polymer = heap
            .allocate(vec![OAD_Polymer {
                units: retained_zero_units,
                n: 0,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        OAD_Polymer_Free(&mut heap, zero_count_polymer).unwrap();
        assert_eq!(
            heap.slice(retained_zero_units.as_const()).unwrap(),
            &[retained_zero_unit]
        );
        assert!(heap.slice(retained_zero_unit.as_const()).is_ok());
        OAD_PolymerUnit_Free(&mut heap, retained_zero_unit).unwrap();
        heap.free(retained_zero_units).unwrap();

        let retained_negative_unit = heap.allocate(vec![OAD_PolymerUnit::default()]).unwrap();
        let negative_units = heap.allocate(vec![retained_negative_unit]).unwrap();
        let negative_count_polymer = heap
            .allocate(vec![OAD_Polymer {
                units: negative_units,
                n: i32::MIN,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        OAD_Polymer_Free(&mut heap, negative_count_polymer).unwrap();
        assert_eq!(
            heap.slice(negative_units.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert!(heap.slice(retained_negative_unit.as_const()).is_ok());
        OAD_PolymerUnit_Free(&mut heap, retained_negative_unit).unwrap();

        let no_units_pzz = heap.allocate(vec![9_i32]).unwrap();
        let no_units_polymer = heap
            .allocate(vec![OAD_Polymer {
                units: SourceMutPointer::null(),
                n: i32::MAX,
                n_pzz: i32::MIN,
                pzz: no_units_pzz,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        OAD_Polymer_Free(&mut heap, no_units_polymer).unwrap();
        assert_eq!(
            heap.slice(no_units_pzz.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_debugtrace__line_3642() {
        OAD_PolymerUnit_DebugTrace(None);

        for connection in [i32::MIN, 0, 1, 2, 3, 4, i32::MAX] {
            for unit_type in [i32::MIN, -1, 0, 1, 2, 3, 4, 5, 6, i32::MAX] {
                for subtype in [i32::MIN, 0, 1, 2, 3, 4, i32::MAX] {
                    let unit = OAD_PolymerUnit {
                        conn: connection,
                        type_: unit_type,
                        subtype,
                        na: i32::MIN,
                        nb: i32::MAX,
                        nbkbonds: i32::MAX,
                        alist: SourceMutPointer::null(),
                        blist: SourceMutPointer::null(),
                        bkbonds: SourceMutPointer::null(),
                        ..OAD_PolymerUnit::default()
                    };
                    let before = unit.clone();
                    OAD_PolymerUnit_DebugTrace(Some(&unit));
                    assert_eq!(unit, before);
                }
            }
        }
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_setreopeningdetails__line_4007() {
        fn run(
            nbkbonds: i32,
            endpoints: [i32; 2],
            atom: inp_ATOM,
            initial_cyclizable: i32,
        ) -> (i32, OAD_PolymerUnit) {
            let mut heap = SourceHeap::default();
            let row = heap.allocate_model_storage(endpoints.to_vec()).unwrap();
            let rows = heap.allocate_model_storage(vec![row]).unwrap();
            let mut unit = OAD_PolymerUnit {
                nbkbonds,
                bkbonds: rows,
                end_atom1: 91,
                end_atom2: 92,
                cyclizable: initial_cyclizable,
                ..OAD_PolymerUnit::default()
            };
            let mut atoms = vec![inp_ATOM::default(); 4];
            atoms[0] = atom;
            let status = OAD_PolymerUnit_SetReopeningDetails(&heap, &mut unit, &atoms).unwrap();
            (status, unit)
        }

        for count in [i32::MIN, -1, 0, 2, i32::MAX] {
            let (status, unit) = run(count, [1, 2], inp_ATOM::default(), 7);
            assert_eq!(status, count);
            assert_eq!(
                (unit.end_atom1, unit.end_atom2, unit.cyclizable),
                (91, 92, 7)
            );
        }

        let (_, diradical) = run(1, [1, 1], inp_ATOM::default(), 7);
        assert_eq!((diradical.end_atom1, diradical.end_atom2), (1, 1));
        assert_eq!(diradical.cyclizable, CLOSING_SRU_DIRADICAL as i32);

        let mut higher_order = inp_ATOM {
            valence: 2,
            ..inp_ATOM::default()
        };
        higher_order.neighbor[..2].copy_from_slice(&[3, 1]);
        higher_order.bond_type[..2].copy_from_slice(&[1, 2]);
        let (_, higher) = run(1, [1, 2], higher_order, 7);
        assert_eq!(higher.cyclizable, CLOSING_SRU_HIGHER_ORDER_BOND as i32);

        let mut first_match_is_single = inp_ATOM {
            valence: 3,
            ..inp_ATOM::default()
        };
        first_match_is_single.neighbor[..3].copy_from_slice(&[1, 1, 2]);
        first_match_is_single.bond_type[..3].copy_from_slice(&[1, 3, 3]);
        let (_, single) = run(1, [1, 2], first_match_is_single, 7);
        assert_eq!(single.cyclizable, 7);

        let mut no_match = inp_ATOM {
            valence: -1,
            ..inp_ATOM::default()
        };
        no_match.neighbor[0] = 1;
        no_match.bond_type[0] = 3;
        let (_, unchanged) = run(1, [1, 2], no_match, 7);
        assert_eq!(unchanged.cyclizable, 7);
    }

    #[test]
    fn source_port__runichi3__oad_polymer_compareranksoftwoatoms__line_4209() {
        let categories = [
            OAD_AtProps {
                erank: 2,
                ring_size: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 2,
                ring_erank: 2,
                ring_size: 3,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 3,
                ring_size: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 3,
                ring_erank: 3,
                ring_size: 3,
                ..OAD_AtProps::default()
            },
        ];
        for first in 0..categories.len() {
            for second in 0..categories.len() {
                let expected = if first > second {
                    -1
                } else if first < second {
                    1
                } else {
                    0
                };
                assert_eq!(
                    OAD_Polymer_CompareRanksOfTwoAtoms(
                        i32::try_from(first + 1).unwrap(),
                        i32::try_from(second + 1).unwrap(),
                        &categories,
                    )
                    .unwrap(),
                    expected
                );
            }
        }

        for (first, second, expected) in [
            (
                OAD_AtProps {
                    ring_erank: 3,
                    ring_size: 4,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    ring_erank: 4,
                    ring_size: 3,
                    ..OAD_AtProps::default()
                },
                1,
            ),
            (
                OAD_AtProps {
                    ring_erank: 4,
                    ring_size: 3,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    ring_erank: 3,
                    ring_size: 4,
                    ..OAD_AtProps::default()
                },
                -1,
            ),
            (
                OAD_AtProps {
                    ring_erank: 3,
                    ring_size: 5,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    ring_erank: 3,
                    ring_size: 4,
                    ..OAD_AtProps::default()
                },
                -1,
            ),
            (
                OAD_AtProps {
                    erank: 8,
                    ring_size: 2,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    erank: 7,
                    ring_size: 2,
                    ..OAD_AtProps::default()
                },
                -1,
            ),
            (
                OAD_AtProps {
                    erank: 2,
                    ring_erank: 2,
                    ring_size: 8,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    erank: 2,
                    ring_erank: 2,
                    ring_size: 7,
                    ..OAD_AtProps::default()
                },
                -1,
            ),
        ] {
            let values = [first, second];
            assert_eq!(
                OAD_Polymer_CompareRanksOfTwoAtoms(1, 2, &values).unwrap(),
                expected
            );
            assert_eq!(
                OAD_Polymer_CompareRanksOfTwoAtoms(2, 1, &values).unwrap(),
                -expected
            );
        }

        assert_eq!(OAD_Polymer_CompareRanksOfTwoAtoms(0, 1, &[]).unwrap(), 0);
        assert_eq!(OAD_Polymer_CompareRanksOfTwoAtoms(-1, 1, &[]).unwrap(), 0);
        assert!(matches!(
            OAD_Polymer_CompareRanksOfTwoAtoms(i32::MIN, 1, &[]),
            Err(SourceHeapError::SourceIntegerOverflow)
        ));
    }

    #[test]
    fn source_port__runichi3__oad_polymer_isfirstatomranklower__line_4369() {
        let properties = [
            OAD_AtProps {
                erank: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 8,
                ..OAD_AtProps::default()
            },
        ];
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(1, 2, &properties).unwrap(),
            1
        );
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(2, 1, &properties).unwrap(),
            -1
        );
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(1, 1, &properties).unwrap(),
            0
        );
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(1, 3, &properties).unwrap(),
            1
        );
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(3, 1, &properties).unwrap(),
            -1
        );
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(0, 1, &properties).unwrap(),
            1
        );
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(-1, -2, &properties).unwrap(),
            -1
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymer_comparebackbonebondsseniority__line_4128() {
        let carbon = OAD_AtProps {
            erank: 2,
            ..OAD_AtProps::default()
        };
        let heteroatom = OAD_AtProps {
            erank: 8,
            ..OAD_AtProps::default()
        };
        let properties = [
            carbon.clone(),
            carbon.clone(),
            carbon.clone(),
            carbon,
            heteroatom,
        ];

        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[1, 5], &[2, 3], &properties).unwrap(),
            -1
        );
        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[2, 3], &[5, 1], &properties).unwrap(),
            1
        );

        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[1, 4], &[2, 5], &properties).unwrap(),
            1
        );
        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[2, 5], &[4, 1], &properties).unwrap(),
            -1
        );

        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[1, 3], &[2, 4], &properties).unwrap(),
            1
        );
        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[2, 4], &[3, 1], &properties).unwrap(),
            -1
        );
        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[4, 1], &[4, 2], &properties).unwrap(),
            -1
        );
        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[4, 2], &[1, 4], &properties).unwrap(),
            1
        );
        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[3, 1], &[1, 3], &properties).unwrap(),
            0
        );

        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[i32::MIN, 1], &[2, 3], &properties),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_sortbackbonebonds__line_4097() {
        let mut heap = SourceHeap::default();
        let row0 = heap.allocate(vec![1_i32, 2]).unwrap();
        let row1 = heap.allocate(vec![3_i32, 5]).unwrap();
        let row2 = heap.allocate(vec![1_i32, 4]).unwrap();
        let row3 = heap.allocate(vec![2_i32, 1]).unwrap();
        let rows = heap.allocate(vec![row0, row1, row2, row3]).unwrap();
        let unit = OAD_PolymerUnit {
            nbkbonds: 4,
            bkbonds: rows,
            ..OAD_PolymerUnit::default()
        };
        let properties = [
            OAD_AtProps {
                erank: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 8,
                ..OAD_AtProps::default()
            },
        ];

        let bond_numbers = heap.allocate(vec![0, 1, 2, 3]).unwrap();
        OAD_PolymerUnit_SortBackboneBonds(&mut heap, &unit, &properties, bond_numbers).unwrap();
        assert_eq!(heap.slice(bond_numbers.as_const()).unwrap(), &[1, 2, 0, 3]);

        let unchanged = heap.allocate(vec![3, 2, 1, 0]).unwrap();
        let no_bonds = OAD_PolymerUnit {
            nbkbonds: 0,
            ..unit.clone()
        };
        OAD_PolymerUnit_SortBackboneBonds(&mut heap, &no_bonds, &properties, unchanged).unwrap();
        assert_eq!(heap.slice(unchanged.as_const()).unwrap(), &[3, 2, 1, 0]);
        OAD_PolymerUnit_SortBackboneBonds(&mut heap, &unit, &properties, SourceMutPointer::null())
            .unwrap();

        let invalid = heap.allocate(vec![0, -1, 2, 3]).unwrap();
        assert_eq!(
            OAD_PolymerUnit_SortBackboneBonds(&mut heap, &unit, &properties, invalid,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_sortbackbonebondsandsetseniors__line_4056() {
        fn properties() -> Vec<OAD_AtProps> {
            vec![
                OAD_AtProps {
                    erank: 2,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    erank: 2,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    erank: 2,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    erank: 2,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    erank: 8,
                    ..OAD_AtProps::default()
                },
            ]
        }

        let mut heap = SourceHeap::default();
        let row0 = heap.allocate(vec![1_i32, 2]).unwrap();
        let row1 = heap.allocate(vec![3_i32, 5]).unwrap();
        let rows = heap.allocate(vec![row0, row1]).unwrap();
        let mut unit = OAD_PolymerUnit {
            nbkbonds: 2,
            bkbonds: rows,
            ..OAD_PolymerUnit::default()
        };
        let mut senior = -1;
        OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(
            &mut heap,
            &mut unit,
            SourceMutPointer::null(),
            &properties(),
            &mut senior,
        )
        .unwrap();
        assert_eq!(senior, 1);
        assert_eq!(&heap.slice(row1.as_const()).unwrap()[..2], &[5, 3]);
        assert_eq!((unit.end_atom1, unit.end_atom2), (5, 3));

        let mut failure_heap = SourceHeap::default();
        let failure_row0 = failure_heap.allocate(vec![1_i32, 2]).unwrap();
        let failure_row1 = failure_heap.allocate(vec![3_i32, 5]).unwrap();
        let failure_rows = failure_heap
            .allocate(vec![failure_row0, failure_row1])
            .unwrap();
        let mut failure_unit = OAD_PolymerUnit {
            nbkbonds: 2,
            bkbonds: failure_rows,
            ..OAD_PolymerUnit::default()
        };
        failure_heap.fail_after_allocations(0);
        senior = -1;
        OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(
            &mut failure_heap,
            &mut failure_unit,
            SourceMutPointer::null(),
            &properties(),
            &mut senior,
        )
        .unwrap();
        assert_eq!(senior, 0);
        assert_eq!(
            &failure_heap.slice(failure_row0.as_const()).unwrap()[..2],
            &[2, 1]
        );
        assert_eq!((failure_unit.end_atom1, failure_unit.end_atom2), (2, 1));

        let mut single_heap = SourceHeap::default();
        let single_row = single_heap.allocate(vec![2_i32, 1]).unwrap();
        let single_rows = single_heap.allocate(vec![single_row]).unwrap();
        let mut single_unit = OAD_PolymerUnit {
            nbkbonds: 1,
            bkbonds: single_rows,
            ..OAD_PolymerUnit::default()
        };
        senior = -1;
        OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(
            &mut single_heap,
            &mut single_unit,
            SourceMutPointer::null(),
            &properties(),
            &mut senior,
        )
        .unwrap();
        assert_eq!(senior, 0);
        assert_eq!(
            &single_heap.slice(single_row.as_const()).unwrap()[..2],
            &[2, 1]
        );
        assert_eq!((single_unit.end_atom1, single_unit.end_atom2), (2, 1));
    }
}
