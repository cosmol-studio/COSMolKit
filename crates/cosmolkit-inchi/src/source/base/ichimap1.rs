use crate::source_types::{
    AB_MAX_KNOWN_PARITY, AB_MIN_KNOWN_PARITY, AB_PARITY_UNDF, AB_PARITY_UNKN, AT_NUMB, AT_RANK,
    AT_STEREO_CARB, AT_STEREO_DBLE, BEST_PARITY, BITS_PARITY, CT_STEREOCOUNT_ERR, CUR_TREE,
    MAX_NUM_STEREO_BONDS, SB_PARITY_MASK, SB_PARITY_SHFT, STEREO_AT_MARK, SourceConstPointer,
    SourceHeap, SourceHeapError, SourceMutPointer, WORSE_PARITY, ppAT_RANK, sp_ATOM,
};

use super::util::{inchi_calloc, inchi_free};

#[allow(non_snake_case)]
pub(crate) fn CurTreeReAlloc(
    heap: &mut SourceHeap,
    cur_tree: Option<&mut CUR_TREE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:850 CurTreeReAlloc
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int CurTreeReAlloc( CUR_TREE *cur_tree )
    {
        if (cur_tree)
        {
            if (cur_tree->tree && cur_tree->max_len > 0 && cur_tree->incr_len > 0)
            {
                void *p = cur_tree->tree;
                if ((cur_tree->tree = (AT_NUMB *) inchi_calloc( (long long)cur_tree->max_len + (long long)cur_tree->incr_len, sizeof( cur_tree->tree[0] ) ))) /* djb-rwth: cast operators added; addressing LLVM warning */
                {
                    memcpy(cur_tree->tree, p, cur_tree->cur_len * sizeof(cur_tree->tree[0]));
                    inchi_free( p );
                    cur_tree->max_len += cur_tree->incr_len;
                    return 0; /*  ok */
                }
            }
        }

        return -1; /*  error */ /*   <BRKPT> */
    }
    */
    // END INCHI C FUNCTION: CurTreeReAlloc

    let Some(cur_tree) = cur_tree else {
        return Ok(-1);
    };
    if cur_tree.tree.is_null() || cur_tree.max_len <= 0 || cur_tree.incr_len <= 0 {
        return Ok(-1);
    }
    let old = cur_tree.tree;
    let new_count = i64::from(cur_tree.max_len) + i64::from(cur_tree.incr_len);
    let replacement = match inchi_calloc::<AT_NUMB>(heap, new_count as u64, 2) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed)
        | Err(SourceHeapError::AllocationElementCountOutOfRange)
        | Err(SourceHeapError::AllocationSizeOverflow) => {
            cur_tree.tree = SourceMutPointer::null();
            return Ok(-1);
        }
        Err(error) => return Err(error),
    };
    cur_tree.tree = replacement;
    let cur_len =
        usize::try_from(cur_tree.cur_len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    heap.with_slice_mut_and_heap(replacement, |new_values, heap| {
        let old_values = heap.slice(old.as_const())?;
        let source = old_values
            .get(..cur_len)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let target = new_values
            .get_mut(..cur_len)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        target.copy_from_slice(source);
        Ok(())
    })?;
    inchi_free(heap, old)?;
    cur_tree.max_len = cur_tree.max_len.wrapping_add(cur_tree.incr_len);
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CurTreeAddRank(
    heap: &mut SourceHeap,
    cur_tree: Option<&mut CUR_TREE>,
    rank: AT_NUMB,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:881 CurTreeAddRank
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int CurTreeAddRank( CUR_TREE *cur_tree, AT_NUMB rank )
    {
        if (cur_tree)
        {
            if (cur_tree->cur_len + 2 > cur_tree->max_len)
            {
                if (CurTreeReAlloc( cur_tree ))
                {
                    return -1; /*  error */ /*   <BRKPT> */
                }
            }
            cur_tree->tree[cur_tree->cur_len++] = rank;
            cur_tree->tree[cur_tree->cur_len++] = 1;
            return 0;
        }

        return -1;  /*  error  */ /*   <BRKPT> */
    }
    */
    // END INCHI C FUNCTION: CurTreeAddRank

    let Some(cur_tree) = cur_tree else {
        return Ok(-1);
    };
    if cur_tree.cur_len.wrapping_add(2) > cur_tree.max_len
        && CurTreeReAlloc(heap, Some(cur_tree))? != 0
    {
        return Ok(-1);
    }
    let first =
        usize::try_from(cur_tree.cur_len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(cur_tree.tree)?
        .get_mut(first)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = rank;
    cur_tree.cur_len = cur_tree.cur_len.wrapping_add(1);
    let second =
        usize::try_from(cur_tree.cur_len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(cur_tree.tree)?
        .get_mut(second)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 1;
    cur_tree.cur_len = cur_tree.cur_len.wrapping_add(1);
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CurTreeIsLastRank(
    heap: &SourceHeap,
    cur_tree: Option<&CUR_TREE>,
    rank: AT_NUMB,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:902 CurTreeIsLastRank
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int CurTreeIsLastRank( CUR_TREE *cur_tree, AT_NUMB rank )
    {
        if (cur_tree && cur_tree->cur_len > 0)
        {
            int rank_pos;
            rank_pos = cur_tree->cur_len - 1;
            rank_pos -= cur_tree->tree[rank_pos];
            if (rank_pos >= 0)
            {
                return ( rank == cur_tree->tree[rank_pos] );
            }
        }

        return 0;  /*  not found */
    }
    */
    // END INCHI C FUNCTION: CurTreeIsLastRank

    let Some(cur_tree) = cur_tree else {
        return Ok(0);
    };
    if cur_tree.cur_len <= 0 {
        return Ok(0);
    }
    let end =
        usize::try_from(cur_tree.cur_len - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let tree = heap.slice(cur_tree.tree.as_const())?;
    let rank_pos = (cur_tree.cur_len - 1).wrapping_sub(i32::from(
        *tree.get(end).ok_or(SourceHeapError::PointerOutOfBounds)?,
    ));
    if rank_pos < 0 {
        return Ok(0);
    }
    let rank_pos = usize::try_from(rank_pos).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    Ok(i32::from(
        rank == *tree
            .get(rank_pos)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    ))
}

#[allow(non_snake_case)]
pub(crate) fn CurTreeRemoveLastRankIfNoAtoms(
    heap: &SourceHeap,
    cur_tree: Option<&mut CUR_TREE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:920 CurTreeRemoveLastRankIfNoAtoms
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int CurTreeRemoveLastRankIfNoAtoms( CUR_TREE *cur_tree )
    {
        if (cur_tree && cur_tree->tree && cur_tree->cur_len >= 2)
        {
            if (1 == cur_tree->tree[cur_tree->cur_len - 1])
            {
                return CurTreeRemoveLastRank( cur_tree ); /*  0=> success, -1=>failed */
            }
            return 1; /*  cannot remove */
        }
        return -1; /*  error */ /*   <BRKPT> */
    }
    */
    // END INCHI C FUNCTION: CurTreeRemoveLastRankIfNoAtoms

    let Some(cur_tree) = cur_tree else {
        return Ok(-1);
    };
    if cur_tree.tree.is_null() || cur_tree.cur_len < 2 {
        return Ok(-1);
    }
    let end =
        usize::try_from(cur_tree.cur_len - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if *heap
        .slice(cur_tree.tree.as_const())?
        .get(end)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        == 1
    {
        CurTreeRemoveLastRank(heap, Some(cur_tree))
    } else {
        Ok(1)
    }
}

#[allow(non_snake_case)]
pub(crate) fn CurTreeAddAtom(
    heap: &mut SourceHeap,
    cur_tree: Option<&mut CUR_TREE>,
    at_no: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:935 CurTreeAddAtom
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int CurTreeAddAtom( CUR_TREE *cur_tree, int at_no )
    {
        if (cur_tree)
        {
            if (cur_tree->cur_len + 1 > cur_tree->max_len)
            {
                if (CurTreeReAlloc( cur_tree ))
                {
                    return -1; /*  error */ /*   <BRKPT> */
                }
            }
            if (cur_tree->cur_len > 0)
            {
                AT_NUMB new_len = cur_tree->tree[--cur_tree->cur_len] + 1;
                cur_tree->tree[cur_tree->cur_len++] = (AT_NUMB) at_no;
                cur_tree->tree[cur_tree->cur_len++] = new_len;
                return 0;
            }
        }

        return -1;
    }
    */
    // END INCHI C FUNCTION: CurTreeAddAtom

    let Some(cur_tree) = cur_tree else {
        return Ok(-1);
    };
    if cur_tree.cur_len.wrapping_add(1) > cur_tree.max_len
        && CurTreeReAlloc(heap, Some(cur_tree))? != 0
    {
        return Ok(-1);
    }
    if cur_tree.cur_len <= 0 {
        return Ok(-1);
    }
    cur_tree.cur_len = cur_tree.cur_len.wrapping_sub(1);
    let marker_index =
        usize::try_from(cur_tree.cur_len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let new_len = heap
        .slice(cur_tree.tree.as_const())?
        .get(marker_index)
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .wrapping_add(1);
    *heap
        .slice_mut(cur_tree.tree)?
        .get_mut(marker_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = at_no as AT_NUMB;
    cur_tree.cur_len = cur_tree.cur_len.wrapping_add(1);
    let new_marker_index =
        usize::try_from(cur_tree.cur_len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(cur_tree.tree)?
        .get_mut(new_marker_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = new_len;
    cur_tree.cur_len = cur_tree.cur_len.wrapping_add(1);
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CurTreeKeepLastAtomsOnly(
    heap: &mut SourceHeap,
    cur_tree: Option<&mut CUR_TREE>,
    tpos: i32,
    shift: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:960 CurTreeKeepLastAtomsOnly
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    void CurTreeKeepLastAtomsOnly( CUR_TREE *cur_tree, int tpos, int shift )
    {   /*  on first entry: shift = 1; other values may occur in subsequent recursion */
        /*  cur_tree[cur_tree->cur_len - shift] is the length of a segment */
        /*  action: remove all atoms except the last from all segments
                    that have length value positon to the right from tpos */
        int cur_length_pos;
        if (cur_tree && cur_tree->tree && ( cur_length_pos = cur_tree->cur_len - shift ) > tpos)
        {
            if (cur_tree->tree[cur_length_pos] > 2)
            {
                /*  current segment contains more than 1 atom. Leave in the segment: rank, the last atom, length value */
                /*  subtract (old segment length)-(new segment length) from the tree length  */
                /*  actual segment length including segment length value = (cur_tree->tree[cur_length_pos]+1) */
                cur_tree->cur_len -= (int) cur_tree->tree[cur_length_pos] - 2;
                memmove(cur_tree->tree + cur_length_pos - cur_tree->tree[cur_length_pos] + 1, /*  1st atom pos */
                    cur_tree->tree + cur_length_pos - 1,  /*  last atom in the current segment position */
                    ((long long)shift + 1) * sizeof(cur_tree->tree[0])); /* djb-rwth: cast operator added */
                /*  (current segment length) distance from the last tree element has not changed */
                cur_tree->tree[cur_tree->cur_len - shift] = 2;
                /*  add 3 to move to the previous segment length position */
                shift += 3; /*  lenghth = 3 accounts for 3 currently present. segment items:
                                (1) the last atom, (2) rank, (3) length value */
            }
            else
            {
                shift += (int) cur_tree->tree[cur_length_pos] + 1; /*  cur_tree->cur_len - (previous segment length position) */
            }
            CurTreeKeepLastAtomsOnly( cur_tree, tpos, shift );
        }
    }
    */
    // END INCHI C FUNCTION: CurTreeKeepLastAtomsOnly

    let Some(cur_tree) = cur_tree else {
        return Ok(());
    };
    if cur_tree.tree.is_null() {
        return Ok(());
    }
    let cur_length_pos = cur_tree.cur_len.wrapping_sub(shift);
    if cur_length_pos <= tpos {
        return Ok(());
    }
    let position =
        usize::try_from(cur_length_pos).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let segment_length = *heap
        .slice(cur_tree.tree.as_const())?
        .get(position)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let next_shift;
    if segment_length > 2 {
        cur_tree.cur_len = cur_tree
            .cur_len
            .wrapping_sub(i32::from(segment_length).wrapping_sub(2));
        let destination = cur_length_pos
            .wrapping_sub(i32::from(segment_length))
            .wrapping_add(1);
        let source = cur_length_pos.wrapping_sub(1);
        let count = i64::from(shift) + 1;
        let destination =
            usize::try_from(destination).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let source = usize::try_from(source).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let count = usize::try_from(count).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let source_end = source
            .checked_add(count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let destination_end = destination
            .checked_add(count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let tree = heap.slice_mut(cur_tree.tree)?;
        if source_end > tree.len() || destination_end > tree.len() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        tree.copy_within(source..source_end, destination);
        let marker = usize::try_from(cur_tree.cur_len.wrapping_sub(shift))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *tree
            .get_mut(marker)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = 2;
        next_shift = shift.wrapping_add(3);
    } else {
        next_shift = shift
            .wrapping_add(i32::from(segment_length))
            .wrapping_add(1);
    }
    CurTreeKeepLastAtomsOnly(heap, Some(cur_tree), tpos, next_shift)
}

#[allow(non_snake_case)]
pub(crate) fn CurTreeRemoveIfLastAtom(
    heap: &mut SourceHeap,
    cur_tree: Option<&mut CUR_TREE>,
    at_no: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:993 CurTreeRemoveIfLastAtom
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int CurTreeRemoveIfLastAtom( CUR_TREE *cur_tree, int at_no )
    {
        if (cur_tree && cur_tree->tree && cur_tree->cur_len > 2)
        {
            AT_NUMB len = cur_tree->tree[cur_tree->cur_len - 1];
            if (len >= 2 && (int) cur_tree->tree[cur_tree->cur_len - 2] == at_no)
            {
                cur_tree->tree[--cur_tree->cur_len - 1] = len - 1;
                return 0;
            }
            return 1; /*  not found */
        }

        return -1; /*  error */ /*   <BRKPT> */
    }
    */
    // END INCHI C FUNCTION: CurTreeRemoveIfLastAtom

    let Some(cur_tree) = cur_tree else {
        return Ok(-1);
    };
    if cur_tree.tree.is_null() || cur_tree.cur_len <= 2 {
        return Ok(-1);
    }
    let end =
        usize::try_from(cur_tree.cur_len - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom_position =
        usize::try_from(cur_tree.cur_len - 2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let tree = heap.slice(cur_tree.tree.as_const())?;
    let len = *tree.get(end).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let atom = *tree
        .get(atom_position)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if len >= 2 && i32::from(atom) == at_no {
        cur_tree.cur_len = cur_tree.cur_len.wrapping_sub(1);
        let marker_position = usize::try_from(cur_tree.cur_len.wrapping_sub(1))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(cur_tree.tree)?
            .get_mut(marker_position)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = len - 1;
        Ok(0)
    } else {
        Ok(1)
    }
}

#[allow(non_snake_case)]
pub(crate) fn SetUseAtomForStereo(
    heap: &mut SourceHeap,
    bAtomUsedForStereo: SourceMutPointer<i8>,
    at: SourceConstPointer<sp_ATOM>,
    num_atoms: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:804 SetUseAtomForStereo
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    void SetUseAtomForStereo( S_CHAR *bAtomUsedForStereo, sp_ATOM *at, int num_atoms )
    {
        int i, k;
        memset( bAtomUsedForStereo, 0, sizeof( bAtomUsedForStereo[0] )*num_atoms ); /* djb-rwth: memset_s C11/Annex K variant? */
        for (i = 0; i < num_atoms; i++)
        {
            if (at[i].parity)
            {
                for (k = 0; k < MAX_NUM_STEREO_BONDS && at[i].stereo_bond_neighbor[k]; k++)
                {
                    ;
                }
                bAtomUsedForStereo[i] = k ? k : STEREO_AT_MARK;
            }
        }
    }
    */
    // END INCHI C FUNCTION: SetUseAtomForStereo

    let count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    heap.with_slice_mut_and_heap(bAtomUsedForStereo, |used, heap| {
        let used = used
            .get_mut(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atoms = heap.slice(at)?;
        let atoms = atoms
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        used.fill(0);
        for (index, atom) in atoms.iter().enumerate() {
            if atom.parity != 0 {
                let mut stereo_count = 0_usize;
                while stereo_count < MAX_NUM_STEREO_BONDS as usize
                    && atom.stereo_bond_neighbor[stereo_count] != 0
                {
                    stereo_count += 1;
                }
                used[index] = if stereo_count != 0 {
                    stereo_count as i8
                } else {
                    STEREO_AT_MARK as i8
                };
            }
        }
        Ok(())
    })
}

#[allow(non_snake_case)]
pub(crate) fn CompareLinCtStereoDoubleToValues(
    heap: &SourceHeap,
    LinearCTStereoDble: SourceConstPointer<AT_STEREO_DBLE>,
    at_rank_canon1: AT_RANK,
    at_rank_canon2: AT_RANK,
    bond_parity: u8,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:764 CompareLinCtStereoDoubleToValues
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int CompareLinCtStereoDoubleToValues( AT_STEREO_DBLE *LinearCTStereoDble,
                                          AT_RANK at_rank_canon1,
                                          AT_RANK at_rank_canon2,
                                          U_CHAR bond_parity )
    {
        if (LinearCTStereoDble->at_num1 CT_GREATER_THAN at_rank_canon1)
        {
            return 1;
        }
        if (LinearCTStereoDble->at_num1 != at_rank_canon1)
        {
            return -1;
        }
        if (LinearCTStereoDble->at_num2 CT_GREATER_THAN at_rank_canon2)
        {
            return 1;
        }
        if (LinearCTStereoDble->at_num2 != at_rank_canon2)
        {
            return -1;
        }
        if (LinearCTStereoDble->parity CT_GREATER_THAN bond_parity)
        {
            return 1;
        }
        if (LinearCTStereoDble->parity != bond_parity)
        {
            return -1;
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: CompareLinCtStereoDoubleToValues

    let value = heap
        .slice(LinearCTStereoDble)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if value.at_num1 > at_rank_canon1 {
        return Ok(1);
    }
    if value.at_num1 != at_rank_canon1 {
        return Ok(-1);
    }
    if value.at_num2 > at_rank_canon2 {
        return Ok(1);
    }
    if value.at_num2 != at_rank_canon2 {
        return Ok(-1);
    }
    if value.parity > bond_parity {
        return Ok(1);
    }
    if value.parity != bond_parity {
        return Ok(-1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CompareLinCtStereoAtomToValues(
    heap: &SourceHeap,
    LinearCTStereoCarb: SourceConstPointer<AT_STEREO_CARB>,
    at_rank_canon1: AT_RANK,
    parity: u8,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:287 CompareLinCtStereoAtomToValues
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int CompareLinCtStereoAtomToValues( AT_STEREO_CARB *LinearCTStereoCarb,
                                        AT_RANK at_rank_canon1,
                                        U_CHAR parity )
    {
        if (LinearCTStereoCarb->at_num CT_GREATER_THAN at_rank_canon1)
        {
            return 1;
        }
        if (LinearCTStereoCarb->at_num != at_rank_canon1)
        {
            return -1;
        }
        if (LinearCTStereoCarb->parity CT_GREATER_THAN parity)
        {
            return 1;
        }
        if (LinearCTStereoCarb->parity != parity)
        {
            return -1;
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: CompareLinCtStereoAtomToValues

    let value = heap
        .slice(LinearCTStereoCarb)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if value.at_num > at_rank_canon1 {
        return Ok(1);
    }
    if value.at_num != at_rank_canon1 {
        return Ok(-1);
    }
    if value.parity > parity {
        return Ok(1);
    }
    if value.parity != parity {
        return Ok(-1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn bUniqueAtNbrFromMappingRank(
    heap: &SourceHeap,
    pRankStack: ppAT_RANK,
    nAtRank: AT_RANK,
    nAtNumber: &mut AT_NUMB,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:316 bUniqueAtNbrFromMappingRank
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int bUniqueAtNbrFromMappingRank( AT_RANK **pRankStack, AT_RANK nAtRank, AT_NUMB *nAtNumber )
    {
        int       r = (int) nAtRank - 1;
        AT_NUMB   i = pRankStack[1][r];
        if (nAtRank == pRankStack[0][(int) i] &&
            ( !r || nAtRank != pRankStack[0][pRankStack[1][r - 1]] )
           )
        {
            *nAtNumber = i;
            return 1;
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: bUniqueAtNbrFromMappingRank

    let stack = heap.slice(pRankStack.as_const())?;
    let ranks = *stack.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
    let order = *stack.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let ranks = heap.slice(ranks.as_const())?;
    let order = heap.slice(order.as_const())?;
    let r = i32::from(nAtRank).wrapping_sub(1);
    let r_index = usize::try_from(r).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = *order
        .get(r_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let atom_rank = *ranks
        .get(usize::from(atom))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let unique_from_previous = if r == 0 {
        true
    } else {
        let previous_atom = *order
            .get(r_index - 1)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        nAtRank
            != *ranks
                .get(usize::from(previous_atom))
                .ok_or(SourceHeapError::PointerOutOfBounds)?
    };
    if nAtRank == atom_rank && unique_from_previous {
        *nAtNumber = atom;
        Ok(1)
    } else {
        Ok(0)
    }
}

#[allow(non_snake_case)]
pub(crate) fn nGetMcr(
    heap: &mut SourceHeap,
    nEqArray: SourceMutPointer<AT_RANK>,
    n: AT_RANK,
) -> Result<AT_RANK, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:336 nGetMcr
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    AT_RANK nGetMcr( AT_RANK *nEqArray, AT_RANK n )
    {
        AT_RANK n1, n2, mcr; /*  recursive version is much shorter. */

        n1 = nEqArray[(int) n];
        if (n == n1)
        {
            return n;
        }
        /*  1st pass: find mcr */
        while (n1 != ( n2 = nEqArray[(int) n1] ))
        {
            n1 = n2;
        }
        /*  2nd pass: copy mcr to each element of the set starting from nEqArray[n] */
        mcr = n1;
        n1 = n;
        while ( /*n1*/ mcr != ( n2 = nEqArray[(int) n1] ))
        {
            nEqArray[(int) n1] = mcr;
            n1 = n2;
        }

        return ( mcr );
    }
    */
    // END INCHI C FUNCTION: nGetMcr

    let get = |heap: &SourceHeap, index: AT_RANK| -> Result<AT_RANK, SourceHeapError> {
        heap.slice(nEqArray.as_const())?
            .get(usize::from(index))
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    };
    let mut n1 = get(heap, n)?;
    if n == n1 {
        return Ok(n);
    }
    loop {
        let n2 = get(heap, n1)?;
        if n1 == n2 {
            break;
        }
        n1 = n2;
    }
    let mcr = n1;
    n1 = n;
    loop {
        let n2 = get(heap, n1)?;
        if mcr == n2 {
            break;
        }
        *heap
            .slice_mut(nEqArray)?
            .get_mut(usize::from(n1))
            .ok_or(SourceHeapError::PointerOutOfBounds)? = mcr;
        n1 = n2;
    }
    Ok(mcr)
}

#[allow(non_snake_case)]
pub(crate) fn nJoin2Mcrs(
    heap: &mut SourceHeap,
    nEqArray: SourceMutPointer<AT_RANK>,
    mut n1: AT_RANK,
    mut n2: AT_RANK,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:366 nJoin2Mcrs
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int nJoin2Mcrs( AT_RANK *nEqArray, AT_RANK n1, AT_RANK n2 )
    {
        n1 = nGetMcr( nEqArray, n1 );
        n2 = nGetMcr( nEqArray, n2 );
        if (n1 < n2)
        {
            nEqArray[n2] = n1;
            return 1; /*  a change has been made */
        }
        if (n2 < n1)
        {
            nEqArray[n1] = n2;
            return 1; /*  a change has been made */
        }

        return 0; /*  no changes */
    }
    */
    // END INCHI C FUNCTION: nJoin2Mcrs

    n1 = nGetMcr(heap, nEqArray, n1)?;
    n2 = nGetMcr(heap, nEqArray, n2)?;
    if n1 < n2 {
        *heap
            .slice_mut(nEqArray)?
            .get_mut(usize::from(n2))
            .ok_or(SourceHeapError::PointerOutOfBounds)? = n1;
        return Ok(1);
    }
    if n2 < n1 {
        *heap
            .slice_mut(nEqArray)?
            .get_mut(usize::from(n1))
            .ok_or(SourceHeapError::PointerOutOfBounds)? = n2;
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn All_SC_Same(
    heap: &SourceHeap,
    canon_rank1: AT_RANK,
    pRankStack1: ppAT_RANK,
    pRankStack2: ppAT_RANK,
    nAtomNumberCanonFrom: SourceConstPointer<AT_RANK>,
    at: SourceConstPointer<sp_ATOM>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:53 All_SC_Same
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int All_SC_Same( AT_RANK canon_rank1, /*  canonical number */
                      const ppAT_RANK pRankStack1,
                      const ppAT_RANK pRankStack2,
                      const AT_RANK *nAtomNumberCanonFrom,
                      const sp_ATOM *at )
    {
        int     n1 = (int) nAtomNumberCanonFrom[(int) canon_rank1 - 1];
        AT_RANK r1 = pRankStack1[0][n1];
        int     iMax1 = (int) r1;
        int     i1, s1;
        int     bFound = 0, stereo_atom_parity = -1;

        /*  find one stereo atom such that canon_rank1 can be mapped on it */
        for (i1 = 1; i1 <= iMax1 && r1 == pRankStack2[0][s1 = (int) pRankStack2[1][iMax1 - i1]]; i1++)
        {
            if (at[s1].stereo_bond_neighbor[0])
            {
                bFound = 0; /* at[s1] is not sp3-stereogenic: it belongs to a stereobond */
                break;
            }
            else
                if (i1 == 1)
                {
                    stereo_atom_parity = PARITY_VAL( at[s1].stereo_atom_parity );
                    if (!ATOM_PARITY_KNOWN( stereo_atom_parity ))
                    {
                        bFound = 0;  /* at[s1] does not have a KNOWN parity */
                        break;
                    }
                }
                else
                    if (stereo_atom_parity != PARITY_VAL( at[s1].stereo_atom_parity ))
                    {
                        bFound = 0; /* two equivalent atoms have different parities */
                        break;
                    }
            bFound++;
        }

        return bFound;
    }
    */
    // END INCHI C FUNCTION: All_SC_Same

    let canon_order = heap.slice(nAtomNumberCanonFrom)?;
    let canon_index = i32::from(canon_rank1).wrapping_sub(1);
    let n1 = *canon_order
        .get(usize::try_from(canon_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let stack1 = heap.slice(pRankStack1.as_const())?;
    let rank1_pointer = *stack1.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
    let rank1 = heap.slice(rank1_pointer.as_const())?;
    let r1 = *rank1
        .get(usize::from(n1))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let stack2 = heap.slice(pRankStack2.as_const())?;
    let rank2_pointer = *stack2.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
    let order2_pointer = *stack2.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let rank2 = heap.slice(rank2_pointer.as_const())?;
    let order2 = heap.slice(order2_pointer.as_const())?;
    let atoms = heap.slice(at)?;

    let iMax1 = i32::from(r1);
    let mut bFound = 0_i32;
    let mut stereo_atom_parity = -1_i32;
    let mut i1 = 1_i32;
    while i1 <= iMax1 {
        let order_index = iMax1.wrapping_sub(i1);
        let s1 = *order2
            .get(usize::try_from(order_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if r1
            != *rank2
                .get(usize::from(s1))
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            break;
        }
        let atom = atoms
            .get(usize::from(s1))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atom.stereo_bond_neighbor[0] != 0 {
            bFound = 0;
            break;
        } else if i1 == 1 {
            stereo_atom_parity = i32::from(atom.stereo_atom_parity) & BITS_PARITY as i32;
            if !(AB_MIN_KNOWN_PARITY as i32..=AB_MAX_KNOWN_PARITY as i32)
                .contains(&stereo_atom_parity)
            {
                bFound = 0;
                break;
            }
        } else if stereo_atom_parity != (i32::from(atom.stereo_atom_parity) & BITS_PARITY as i32) {
            bFound = 0;
            break;
        }
        bFound = bFound.wrapping_add(1);
        i1 = i1.wrapping_add(1);
    }
    Ok(bFound)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn Next_SC_At_CanonRank2(
    heap: &SourceHeap,
    canon_rank1: &mut AT_RANK,
    canon_rank1_min: &mut AT_RANK,
    bFirstTime: &mut i32,
    bAtomUsedForStereo: SourceConstPointer<i8>,
    pRankStack1: ppAT_RANK,
    pRankStack2: ppAT_RANK,
    nAtomNumberCanonFrom: SourceConstPointer<AT_RANK>,
    num_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:99 Next_SC_At_CanonRank2
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int Next_SC_At_CanonRank2( AT_RANK *canon_rank1,        /*  1st call input: largest canon number mapped so far or 0 */
                                                            /*  output: suggested canon. rank > than input if success */
                               AT_RANK *canon_rank1_min,    /*  1st call:0 next calls: first tried canon. number */
                               int *bFirstTime,             /*  1 at the time of the 1st call  */
                               S_CHAR *bAtomUsedForStereo,  /*  STEREO_AT_MARK if the atom has not been mapped yet */
                               const ppAT_RANK pRankStack1, /*  mapping ranks/sort order of atoms with canon. numbers (from) */
                               const ppAT_RANK pRankStack2, /*  mapping ranks/sort order of atoms with stereo (to) */
                               const AT_RANK *nAtomNumberCanonFrom, /*  sorted order of the canon. numbers */
                               int num_atoms )
    {
        AT_RANK canon_rank1_inp = *canon_rank1;
        AT_RANK cr1;  /*  canonical rank (canonical number) */
        AT_RANK r1;   /*  mapping rank */
        int     n1;   /*  ord. number of an atom with the canon. number */
        int     s1;   /*  ord. number of an atom with stereo */
        int     i1, bFound = 0;
        int     iMax1;

        if (canon_rank1_inp < *canon_rank1_min)
        {
            canon_rank1_inp = *canon_rank1_min;
        }
        else
        {
            if (canon_rank1_inp < 1)
            {
                canon_rank1_inp = 1;
            }
            else
            {
                canon_rank1_inp++; /*  next canonical rank */
            }
        }
        cr1 = canon_rank1_inp;

        while ((int) cr1 <= num_atoms)
        {
            n1 = (int) nAtomNumberCanonFrom[(int) cr1 - 1]; /*  atom1 (which has canon. rank cr1) ord. number */
            iMax1 = (int) ( r1 = pRankStack1[0][n1] ); /*  mapping rank of atom1 */
            /*  find atoms "to" to which the canon. number can be mapped; they have mapping rank r1, number s1 */
            for (i1 = 1; i1 <= iMax1 && r1 == pRankStack2[0][s1 = (int) pRankStack2[1][iMax1 - i1]]; i1++)
            {
                /*  looking for a stereo center atom that has mapping rank r1 */
                if (bAtomUsedForStereo[s1] == STEREO_AT_MARK)
                {
                    /*  found a sterogenic atom that has not been mapped yet */
                    bFound = 1;
                    break;
                }
            }
            if (bFound)
            {
                /*  one sterogenic not mapped yet atom "to" has been found */
                if (*bFirstTime)
                {
                    *canon_rank1_min = cr1;
                    *bFirstTime = 0;
                }
                break;
            }
            else
            {
                 /*  a not mapped yet stereogenic atom has not found */
                 /*  for the mapping rank r1 defined by the canonical rank cr1; try next cr1 */
                cr1++;
            }
        }
        if (bFound)
        {
            /*  success */
            *canon_rank1 = cr1;
            return 1;
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: Next_SC_At_CanonRank2

    let mut canon_rank1_inp = *canon_rank1;
    if canon_rank1_inp < *canon_rank1_min {
        canon_rank1_inp = *canon_rank1_min;
    } else if canon_rank1_inp < 1 {
        canon_rank1_inp = 1;
    } else {
        canon_rank1_inp = canon_rank1_inp.wrapping_add(1);
    }
    let mut cr1 = canon_rank1_inp;
    let mut bFound = 0_i32;

    if i32::from(cr1) <= num_atoms {
        let canon_order = heap.slice(nAtomNumberCanonFrom)?;
        let stack1 = heap.slice(pRankStack1.as_const())?;
        let rank1_pointer = *stack1.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
        let rank1 = heap.slice(rank1_pointer.as_const())?;
        let stack2 = heap.slice(pRankStack2.as_const())?;
        let rank2_pointer = *stack2.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
        let order2_pointer = *stack2.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let rank2 = heap.slice(rank2_pointer.as_const())?;
        let order2 = heap.slice(order2_pointer.as_const())?;
        let used = heap.slice(bAtomUsedForStereo)?;

        while i32::from(cr1) <= num_atoms {
            let canon_index = i32::from(cr1).wrapping_sub(1);
            let n1 = *canon_order
                .get(
                    usize::try_from(canon_index)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let r1 = *rank1
                .get(usize::from(n1))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let iMax1 = i32::from(r1);
            let mut i1 = 1_i32;
            while i1 <= iMax1 {
                let order_index = iMax1.wrapping_sub(i1);
                let s1 = *order2
                    .get(
                        usize::try_from(order_index)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if r1
                    != *rank2
                        .get(usize::from(s1))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                {
                    break;
                }
                if *used
                    .get(usize::from(s1))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    == STEREO_AT_MARK as i8
                {
                    bFound = 1;
                    break;
                }
                i1 = i1.wrapping_add(1);
            }
            if bFound != 0 {
                if *bFirstTime != 0 {
                    *canon_rank1_min = cr1;
                    *bFirstTime = 0;
                }
                break;
            }
            cr1 = cr1.wrapping_add(1);
        }
    }
    if bFound != 0 {
        *canon_rank1 = cr1;
        Ok(1)
    } else {
        Ok(0)
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn NextStereoParity2Test(
    stereo_bond_parity: &mut i32,
    sb_parity_calc: &mut i32,
    nNumBest: i32,
    nNumWorse: i32,
    nNumUnkn: i32,
    nNumUndf: i32,
    nNumCalc: i32,
    vABParityUnknown: i32,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:659 NextStereoParity2Test
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int NextStereoParity2Test( int *stereo_bond_parity,
                               int *sb_parity_calc,
                               int nNumBest,
                               int nNumWorse,
                               int nNumUnkn,
                               int nNumUndf,
                               int nNumCalc,
                               int vABParityUnknown )
    {
        /* sequence of (stereo_bond_parity, sb_parity_calc) pairs:

              (BEST_PARITY, BEST_PARITY)  <calc>
                          |
              (BEST_PARITY, WORSE_PARITY) <known>
                          |
              (WORSE_PARITY, WORSE_PARITY) <calc>                (BEST_PARITY, 0) <known>
                           \___________________________________________/
                                                  |
                                           (WORSE_PARITY, 0)   <known>
                                                  |
                                           (AB_PARITY_UNKN, 0) <known>
                                                  |
                                           (AB_PARITY_UNDF, 0) <known>
                                                  |
                                           <next pair of ranks>
          Meaning:
          stereo_bond_parity is the parity we are looking for
          stereo_bond_parity==sb_parity_calc  => parity to be calculated from canonical numbers
          stereo_bond_parity!=sb_parity_calc  => parity is already known
         */
    get_next_parity:
        switch (*stereo_bond_parity)
        {
            case BEST_PARITY:
                switch (*sb_parity_calc)
                {
                    case 0:                                 /*  BEST_PARITY(known) : (BEST_PARITY, 0) -> */
                        *stereo_bond_parity = WORSE_PARITY;  /*  WORSE_PARITY(known): (WORSE_PARITY, 0) */
                        if (!nNumWorse)
                        {
                            goto get_next_parity;
                        }
                        break;
                    case BEST_PARITY:                       /*  BEST_PARITY(calc) : (BEST_PARITY, BEST_PARITY) -> */
                        *sb_parity_calc = WORSE_PARITY;      /*  BEST_PARITY(known): (BEST_PARITY, WORSE_PARITY) */
                        if (!nNumBest)
                        {
                            goto get_next_parity;
                        }
                        break;
                    case WORSE_PARITY:                      /*  BEST_PARITY(known): (BEST_PARITY, WORSE_PARITY)-> */
                        *stereo_bond_parity = WORSE_PARITY;  /*  WORSE_PARITY(calc): (WORSE_PARITY,WORSE_PARITY) */
                        if (!nNumCalc)
                        { /* added 12-17-2003 */
                            goto get_next_parity;
                        }
                        break;
                }
                break;
            case WORSE_PARITY:
                switch (*sb_parity_calc)
                {
                    case 0:                                 /*  WORSE_PARITY(known)  : (WORSE_PARITY, 0) -> */
                        *stereo_bond_parity = vABParityUnknown /* AB_PARITY_UNKN */;/*  AB_PARITY_UNKN(known): (AB_PARITY_UNKN, 0) */
                        if (!nNumUnkn)
                        {
                            goto get_next_parity;
                        }
                        break;
                    case BEST_PARITY:                       /*  error */
                        return CT_STEREOCOUNT_ERR;          /*   <BRKPT> */
                    case WORSE_PARITY:                      /*  WORSE_PARITY(calc) : (WORSE_PARITY,WORSE_PARITY)-> */
                        *sb_parity_calc = 0;                 /*  WORSE_PARITY(known): (WORSE_PARITY, 0) */
                        if (!nNumWorse)
                        {
                            goto get_next_parity;
                        }
                        break;
                }
                break;

            case AB_PARITY_UNKN:                        /* AB_PARITY_UNKN(known): (AB_PARITY_UNKN, 0) -> */
                if (*sb_parity_calc)                 /*  error */
                {
                    return CT_STEREOCOUNT_ERR;          /*   <BRKPT> */
                }
                *stereo_bond_parity = AB_PARITY_UNDF;    /* AB_PARITY_UNDF(known): (AB_PARITY_UNDF, 0) */
                if (!nNumUndf)
                {
                    return 1; /*goto next_canon_ranks;*/
                }
                break;

            case AB_PARITY_UNDF:                        /*  AB_PARITY_UNDF(known): (AB_PARITY_UNDF, 0) -> */
                if (*sb_parity_calc)
                {                /*  error */
                    return CT_STEREOCOUNT_ERR;          /*   <BRKPT> */
                }
                return 1; /*goto next_canon_ranks;*/     /*  next canon ranks */
        }
        return 0;
    }
    */
    // END INCHI C FUNCTION: NextStereoParity2Test

    loop {
        match *stereo_bond_parity {
            value if value == BEST_PARITY as i32 => match *sb_parity_calc {
                0 => {
                    *stereo_bond_parity = WORSE_PARITY as i32;
                    if nNumWorse == 0 {
                        continue;
                    }
                }
                value if value == BEST_PARITY as i32 => {
                    *sb_parity_calc = WORSE_PARITY as i32;
                    if nNumBest == 0 {
                        continue;
                    }
                }
                value if value == WORSE_PARITY as i32 => {
                    *stereo_bond_parity = WORSE_PARITY as i32;
                    if nNumCalc == 0 {
                        continue;
                    }
                }
                _ => {}
            },
            value if value == WORSE_PARITY as i32 => match *sb_parity_calc {
                0 => {
                    *stereo_bond_parity = vABParityUnknown;
                    if nNumUnkn == 0 {
                        continue;
                    }
                }
                value if value == BEST_PARITY as i32 => return CT_STEREOCOUNT_ERR,
                value if value == WORSE_PARITY as i32 => {
                    *sb_parity_calc = 0;
                    if nNumWorse == 0 {
                        continue;
                    }
                }
                _ => {}
            },
            value if value == AB_PARITY_UNKN as i32 => {
                if *sb_parity_calc != 0 {
                    return CT_STEREOCOUNT_ERR;
                }
                *stereo_bond_parity = AB_PARITY_UNDF as i32;
                if nNumUndf == 0 {
                    return 1;
                }
            }
            value if value == AB_PARITY_UNDF as i32 => {
                if *sb_parity_calc != 0 {
                    return CT_STEREOCOUNT_ERR;
                }
                return 1;
            }
            _ => {}
        }
        return 0;
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn Next_SB_At_CanonRanks2(
    heap: &SourceHeap,
    canon_rank1: &mut AT_RANK,
    canon_rank2: &mut AT_RANK,
    canon_rank1_min: &mut AT_RANK,
    canon_rank2_min: &mut AT_RANK,
    bFirstTime: &mut i32,
    bAtomUsedForStereo: SourceConstPointer<i8>,
    pRankStack1: ppAT_RANK,
    pRankStack2: ppAT_RANK,
    nCanonRankFrom: SourceConstPointer<AT_RANK>,
    nAtomNumberCanonFrom: SourceConstPointer<AT_RANK>,
    at: SourceConstPointer<sp_ATOM>,
    num_atoms: i32,
    bAllene: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:518 Next_SB_At_CanonRanks2
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int Next_SB_At_CanonRanks2( AT_RANK *canon_rank1,
                                AT_RANK *canon_rank2, /*  canonical numbers */
                                AT_RANK *canon_rank1_min,
                                AT_RANK *canon_rank2_min,
                                int *bFirstTime,
                                S_CHAR *bAtomUsedForStereo,
                                const ppAT_RANK pRankStack1,
                                const ppAT_RANK pRankStack2,
                                const AT_RANK *nCanonRankFrom,
                                const AT_RANK *nAtomNumberCanonFrom,
                                const sp_ATOM *at,
                                int num_atoms,
                                int bAllene )
    {
        AT_RANK canon_rank1_inp = *canon_rank1;
        AT_RANK canon_rank2_inp = *canon_rank2;
        AT_RANK cr1, cr2; /*  canonical ranks (canonical numbers) */
        AT_RANK r1, r2;   /*  mapping ranks */
        int     n1, n2;   /*  ord. numbers of atoms with stereo */
        int     s1, s2;   /*  ord. numbers of atoms with canon. numbers */
        int     i1, i2, k, m;
        int     iMax1, iMax2;

        if (canon_rank1_inp < *canon_rank1_min ||
             (canon_rank1_inp == *canon_rank1_min &&
             canon_rank2_inp < *canon_rank2_min)) /* djb-rwth: addressing LLVM warning */
        {

            canon_rank1_inp = *canon_rank1_min;
            canon_rank2_inp = *canon_rank2_min;
        }
        else
            if (canon_rank1_inp < 2)
            {
                canon_rank1_inp = 2;
                canon_rank2_inp = 0;
            }
        cr1 = canon_rank1_inp;
        cr2 = num_atoms; /* initialize. 1/8/2002 */
        while ((int) cr1 <= num_atoms)
        {
            cr2 = cr1;
            n1 = (int) nAtomNumberCanonFrom[(int) cr1 - 1]; /*  atom1=at[n1] (which has canon. rank) ord. number */
            iMax1 = (int) ( r1 = pRankStack1[0][n1] ); /*  mapping rank of atom1 */
            for (i1 = 1; i1 <= iMax1 && r1 == pRankStack2[0][s1 = (int) pRankStack2[1][iMax1 - i1]]; i1++)
            {
                /*  looking for a stereo bond atom that has mapping rank r1 */
                /*  found at[s1] such that rank cr1 can be mapped on at[s1] because cr1 and s1 have equal */
                /*  mapping rank = r1. Check at[s1] stereo bonds */
                if (bAtomUsedForStereo[s1] && bAtomUsedForStereo[s1] < STEREO_AT_MARK)
                {
                    for (k = 0; k < MAX_NUM_STEREO_BONDS && ( s2 = (int) at[s1].stereo_bond_neighbor[k] ); k++) /* djb-rwth: removing redundant code */
                    {
                        /*  stereo bond at[s1]-at[s2] has been found */
                        if (bAtomUsedForStereo[--s2])
                        {
                            /*  stereo bonds have not been mapped. however, this check is not needed */
                            int cumulene_len = BOND_CHAIN_LEN( at[s1].stereo_bond_parity[k] );
                            if ((cumulene_len % 2 && !bAllene) || /* 09-26-2003 */
                                 (!( cumulene_len % 2 ) && bAllene)) /* djb-rwth: addressing LLVM warning */
                            { /* 08-17-2003 Fix05 */
                                continue;
                            }
                            iMax2 = (int) ( r2 = pRankStack2[0][s2] ); /*  mapping rank of atom2 */
                            /*  Go back to canonical ranks and find an atom that has mapping rank r2 */
                            /*  and is connected to the atom with canonical rank cr1 (possibly by cumulene chain) */
                            /*  These cr1-cr2 canon. ranks possibly can be mapped on at[s1]-at[s2] stereo bond */
                            for (i2 = 1; i2 <= iMax2 && r2 == pRankStack1[0][n2 = (int) pRankStack1[1][iMax2 - i2]]; i2++)
                            {
                                if (cumulene_len)
                                {
                                    int prev, next, len, j;
                                    for (m = 0; m < at[n1].valence; m++)
                                    {
                                        for (prev = n1, len = 0, next = (int) at[n1].neighbor[m]; len < cumulene_len; len++)
                                        {
                                            if (at[next].valence == 2 && !at[next].num_H)
                                            {
                                                j = ( (int) at[next].neighbor[0] == prev );
                                                prev = next;
                                                next = at[next].neighbor[j];
                                            }
                                            else
                                            {
                                                break; /*  cannot continue */
                                            }
                                        }
                                        if (len == cumulene_len && n2 == next)
                                        {
                                            break;
                                        }
                                    }
                                }
                                else
                                {
                                    for (m = 0; m < at[n1].valence && n2 != (int) at[n1].neighbor[m]; m++)
                                        ;
                                }
                                if (m < at[n1].valence &&
                                     nCanonRankFrom[n2] < cr2 &&
                                     nCanonRankFrom[n2] > canon_rank2_inp)
                                {

                                    cr2 = nCanonRankFrom[n2]; /*  found a candidate for cr2 */
                                }
                            }
                        }
                    }
                }
            }
            if (cr2 >= cr1)
            {
                /*  not found for this r1 */
                cr1++;
                canon_rank2_inp = 0;
            }
            else
            {
                 /* found cr2 < cr1 */
                if (*bFirstTime)
                {
                    *canon_rank1_min = cr1;
                    *canon_rank2_min = cr2;
                    *bFirstTime = 0;
                }
                break;
            }
        }
        if (cr1 > cr2 && cr1 <= num_atoms)
        {
            /*  success */
            *canon_rank1 = cr1;
            *canon_rank2 = cr2;
            return 1;
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: Next_SB_At_CanonRanks2

    let mut canon_rank1_inp = *canon_rank1;
    let mut canon_rank2_inp = *canon_rank2;
    if canon_rank1_inp < *canon_rank1_min
        || (canon_rank1_inp == *canon_rank1_min && canon_rank2_inp < *canon_rank2_min)
    {
        canon_rank1_inp = *canon_rank1_min;
        canon_rank2_inp = *canon_rank2_min;
    } else if canon_rank1_inp < 2 {
        canon_rank1_inp = 2;
        canon_rank2_inp = 0;
    }

    let stack1 = heap.slice(pRankStack1.as_const())?;
    let stack2 = heap.slice(pRankStack2.as_const())?;
    let rank1 = heap.slice(
        stack1
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .as_const(),
    )?;
    let order1 = heap.slice(
        stack1
            .get(1)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .as_const(),
    )?;
    let rank2 = heap.slice(
        stack2
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .as_const(),
    )?;
    let order2 = heap.slice(
        stack2
            .get(1)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .as_const(),
    )?;
    let canon_from = heap.slice(nCanonRankFrom)?;
    let atom_from_canon = heap.slice(nAtomNumberCanonFrom)?;
    let atoms = heap.slice(at)?;
    let used = heap.slice(bAtomUsedForStereo)?;

    let mut cr1 = canon_rank1_inp;
    let mut cr2 = num_atoms as AT_RANK;
    while i32::from(cr1) <= num_atoms {
        cr2 = cr1;
        let n1 = usize::from(
            *atom_from_canon
                .get(
                    usize::from(cr1)
                        .checked_sub(1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let r1 = *rank1.get(n1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let i_max1 = usize::from(r1);
        for i1_count in 1..=i_max1 {
            let s1 = usize::from(
                *order2
                    .get(i_max1 - i1_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            if *rank2.get(s1).ok_or(SourceHeapError::PointerOutOfBounds)? != r1 {
                break;
            }
            let used_s1 = *used.get(s1).ok_or(SourceHeapError::PointerOutOfBounds)?;
            if used_s1 == 0 || i32::from(used_s1) >= STEREO_AT_MARK as i32 {
                continue;
            }
            for k in 0..MAX_NUM_STEREO_BONDS as usize {
                let neighbor = atoms
                    .get(s1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .stereo_bond_neighbor[k];
                if neighbor == 0 {
                    break;
                }
                let s2 = usize::from(neighbor - 1);
                if *used.get(s2).ok_or(SourceHeapError::PointerOutOfBounds)? == 0 {
                    continue;
                }
                let parity = i32::from(
                    atoms
                        .get(s1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .stereo_bond_parity[k],
                );
                let cumulene_len = (parity & 0x38) / (1_i32 << SB_PARITY_SHFT);
                if (cumulene_len % 2 != 0 && bAllene == 0)
                    || (cumulene_len % 2 == 0 && bAllene != 0)
                {
                    continue;
                }
                let r2 = *rank2.get(s2).ok_or(SourceHeapError::PointerOutOfBounds)?;
                let i_max2 = usize::from(r2);
                for i2_count in 1..=i_max2 {
                    let n2 = usize::from(
                        *order1
                            .get(i_max2 - i2_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    if *rank1.get(n2).ok_or(SourceHeapError::PointerOutOfBounds)? != r2 {
                        break;
                    }
                    let atom1 = atoms.get(n1).ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let valence = i32::from(atom1.valence);
                    let mut m = 0_i32;
                    if cumulene_len != 0 {
                        while m < valence {
                            let mut prev = n1;
                            let mut next = usize::from(
                                *atom1
                                    .neighbor
                                    .get(
                                        usize::try_from(m)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            );
                            let mut len = 0_i32;
                            while len < cumulene_len {
                                let next_atom =
                                    atoms.get(next).ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if next_atom.valence == 2 && next_atom.num_H == 0 {
                                    let j = usize::from(usize::from(next_atom.neighbor[0]) == prev);
                                    prev = next;
                                    next = usize::from(next_atom.neighbor[j]);
                                } else {
                                    break;
                                }
                                len += 1;
                            }
                            if len == cumulene_len && n2 == next {
                                break;
                            }
                            m += 1;
                        }
                    } else {
                        while m < valence {
                            if n2
                                == usize::from(
                                    *atom1
                                        .neighbor
                                        .get(
                                            usize::try_from(m)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                                )
                            {
                                break;
                            }
                            m += 1;
                        }
                    }
                    if m < valence {
                        let candidate = *canon_from
                            .get(n2)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if candidate < cr2 && candidate > canon_rank2_inp {
                            cr2 = candidate;
                        }
                    }
                }
            }
        }
        if cr2 >= cr1 {
            cr1 = cr1.wrapping_add(1);
            canon_rank2_inp = 0;
        } else {
            if *bFirstTime != 0 {
                *canon_rank1_min = cr1;
                *canon_rank2_min = cr2;
                *bFirstTime = 0;
            }
            break;
        }
    }
    if cr1 > cr2 && i32::from(cr1) <= num_atoms {
        *canon_rank1 = cr1;
        *canon_rank2 = cr2;
        Ok(1)
    } else {
        Ok(0)
    }
}

#[allow(non_snake_case)]
pub(crate) fn All_SB_Same(
    heap: &SourceHeap,
    canon_rank1: AT_RANK,
    canon_rank2: AT_RANK,
    pRankStack1: ppAT_RANK,
    pRankStack2: ppAT_RANK,
    nAtomNumberCanonFrom: SourceConstPointer<AT_RANK>,
    at: SourceConstPointer<sp_ATOM>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:393 All_SB_Same
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int All_SB_Same( AT_RANK canon_rank1,
                     AT_RANK canon_rank2, /*  canonical numbers */
                     const ppAT_RANK pRankStack1,
                     const ppAT_RANK pRankStack2,
                     const AT_RANK *nAtomNumberCanonFrom,
                     sp_ATOM *at )
    {
        int     n1 = (int) nAtomNumberCanonFrom[(int) canon_rank1 - 1]; /* at1 has canon_rank1 */
        int     n2 = (int) nAtomNumberCanonFrom[(int) canon_rank2 - 1]; /* at2 has canon_rank2 */
        AT_RANK r1 = pRankStack1[0][n1]; /* at1 mapping rank */
        AT_RANK r2 = pRankStack1[0][n2]; /* at2 mapping rank */
        AT_RANK rNeigh1, rNeigh2;
        int     iMax1 = (int) r1;
        /* int     iMax2 = (int)r2; */
        int     i1, i2, s1 = 0, s2 = 0, k1 = 0, k2, m, k, num_equal;
        int     bNotFound = 1, cumulene_len, stereo_bond_parity;

        /*  at the first atom that possibly may have canon_rank1 find one stereo bond such that */
        /*  canon_rank1-canon_rank2 possibly may be mapped on it */
        for (i1 = 1; i1 <= iMax1 && r1 == pRankStack2[0][s1 = (int) pRankStack2[1][iMax1 - i1]]; i1++)
        {
            /* at[n1] may be possible to map on at[s1] */
            for (k1 = 0, s2 = 0, bNotFound = 1;
                  k1 < MAX_NUM_STEREO_BONDS && ( s2 = (int) at[s1].stereo_bond_neighbor[k1] ) &&
                  ( bNotFound = ( r2 != pRankStack2[0][--s2] ) ); k1++)
                ; /* continue until the 1st at[s2] (to which at[n2] may be mapped) have been found */
            if (!bNotFound)
            {
                break; /* stop at 1st found */
            }
        }
        if (bNotFound)
        {
            return -1; /*  error: no mapping exists */
        }
        for (k2 = 0, m = 0; k2 < MAX_NUM_STEREO_BONDS && ( m = (int) at[s2].stereo_bond_neighbor[k2] ) && m - 1 != s1; k2++)
            ;
        if (m - 1 != s1)
        {
            return -1; /*  program error: stereo bond in opposite direction not found */
        }
        stereo_bond_parity = at[s1].stereo_bond_parity[k1];
        if (!PARITY_KNOWN( stereo_bond_parity ))
        {
            return 0;
        }
        cumulene_len = BOND_CHAIN_LEN( stereo_bond_parity );
        rNeigh1 = pRankStack2[0][(int) at[s1].neighbor[(int) at[s1].stereo_bond_ord[k1]]];
        rNeigh2 = pRankStack2[0][(int) at[s2].neighbor[(int) at[s2].stereo_bond_ord[k2]]];

        num_equal = 0;
        /*  Search among ALL neighbors because sometimes a stereo bond may be mapped on a non-stereo bond. */
        /*  If is so then return 0: not all mappings are stereo-equivalent */
        for (s1 = 1; s1 <= iMax1 && r1 == pRankStack2[0][i1 = (int) pRankStack2[1][iMax1 - s1]]; s1++)
        {
            for (k = 0; k < at[i1].valence; k++)
            {
                n1 = at[i1].neighbor[k];
                if (rNeigh1 != pRankStack2[0][n1])
                {
                    continue; /*  wrong neighbor */
                }
                if (cumulene_len)
                {
                    int prev, next, len, j;
                    for (prev = i1, len = 0, next = n1; len < cumulene_len; len++)
                    {
                        if (at[next].valence == 2 && !at[next].num_H)
                        {
                            j = ( (int) at[next].neighbor[0] == prev );
                            prev = next;
                            next = at[next].neighbor[j];
                        }
                        else
                        {
                            break; /*  cannot continue */
                        }
                    }
                    if (len != cumulene_len ||
                         r2 != pRankStack2[0][next] ||
                         rNeigh2 != pRankStack2[0][prev])
                    {
                        /*  cumulene chain not found */
                        continue;
                    }
                    i2 = next;
                }
                else
                {
                    i2 = n1;
                }
                /*  find if a stereogenic bond between at[i1]-at[i2] exists */
                for (k1 = 0; k1 < MAX_NUM_STEREO_BONDS &&
                    ( m = (int) at[i1].stereo_bond_neighbor[k1] ) && m - 1 != i2; k1++)
                    ;
                if (m - 1 != i2)
                {
                    return 0;
                }
                for (k2 = 0; k2 < MAX_NUM_STEREO_BONDS &&
                    ( m = (int) at[i2].stereo_bond_neighbor[k2] ) && m - 1 != i1; k2++)
                    ;
                if (m - 1 != i1)
                {
                    return 0;
                }
                if (at[i1].stereo_bond_parity[k1] != at[i2].stereo_bond_parity[k2])
                {
                    return -1; /*  program error */
                }
                if (stereo_bond_parity != at[i1].stereo_bond_parity[k1])
                {
                    return 0;
                }
                num_equal++;
            }
        }

        return num_equal;
    }
    */
    // END INCHI C FUNCTION: All_SB_Same

    let stack1 = heap.slice(pRankStack1.as_const())?;
    let stack2 = heap.slice(pRankStack2.as_const())?;
    let rank1 = heap.slice(
        stack1
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .as_const(),
    )?;
    let rank2 = heap.slice(
        stack2
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .as_const(),
    )?;
    let order2 = heap.slice(
        stack2
            .get(1)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .as_const(),
    )?;
    let canon_from = heap.slice(nAtomNumberCanonFrom)?;
    let atoms = heap.slice(at)?;

    let canon1 = usize::from(canon_rank1)
        .checked_sub(1)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let canon2 = usize::from(canon_rank2)
        .checked_sub(1)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut n1 = usize::from(
        *canon_from
            .get(canon1)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let n2 = usize::from(
        *canon_from
            .get(canon2)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let r1 = *rank1.get(n1).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let r2 = *rank1.get(n2).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let i_max1 = usize::from(r1);

    let mut mapped_s1 = 0_usize;
    let mut mapped_s2 = 0_usize;
    let mut mapped_k1 = 0_usize;
    let mut not_found = true;
    for i1 in 1..=i_max1 {
        let candidate = usize::from(
            *order2
                .get(i_max1 - i1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        if *rank2
            .get(candidate)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != r1
        {
            break;
        }
        mapped_s1 = candidate;
        mapped_s2 = 0;
        mapped_k1 = 0;
        not_found = true;
        while mapped_k1 < MAX_NUM_STEREO_BONDS as usize {
            let neighbor = atoms
                .get(mapped_s1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .stereo_bond_neighbor[mapped_k1];
            if neighbor == 0 {
                break;
            }
            mapped_s2 = usize::from(neighbor - 1);
            not_found = *rank2
                .get(mapped_s2)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                != r2;
            if !not_found {
                break;
            }
            mapped_k1 += 1;
        }
        if !not_found {
            break;
        }
    }
    if not_found {
        return Ok(-1);
    }

    let mut k2 = 0_usize;
    let mut m = 0_usize;
    while k2 < MAX_NUM_STEREO_BONDS as usize {
        m = usize::from(
            atoms
                .get(mapped_s2)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .stereo_bond_neighbor[k2],
        );
        if m == 0 || m - 1 == mapped_s1 {
            break;
        }
        k2 += 1;
    }
    if m == 0 || m - 1 != mapped_s1 {
        return Ok(-1);
    }

    let stereo_bond_parity = i32::from(
        atoms
            .get(mapped_s1)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .stereo_bond_parity[mapped_k1],
    );
    let parity_value = stereo_bond_parity & SB_PARITY_MASK as i32;
    if !(AB_MIN_KNOWN_PARITY as i32..=AB_MAX_KNOWN_PARITY as i32).contains(&parity_value) {
        return Ok(0);
    }
    let cumulene_len = (stereo_bond_parity & 0x38) / (1_i32 << SB_PARITY_SHFT);
    let first_atom = atoms
        .get(mapped_s1)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let first_order = usize::try_from(first_atom.stereo_bond_ord[mapped_k1])
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let r_neigh1 = *rank2
        .get(usize::from(
            *first_atom
                .neighbor
                .get(first_order)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        ))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second_atom = atoms
        .get(mapped_s2)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second_order = usize::try_from(second_atom.stereo_bond_ord[k2])
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let r_neigh2 = *rank2
        .get(usize::from(
            *second_atom
                .neighbor
                .get(second_order)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        ))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;

    let mut num_equal = 0_i32;
    for s1_count in 1..=i_max1 {
        let i1 = usize::from(
            *order2
                .get(i_max1 - s1_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        if *rank2.get(i1).ok_or(SourceHeapError::PointerOutOfBounds)? != r1 {
            break;
        }
        let atom1 = atoms.get(i1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let valence =
            usize::try_from(atom1.valence).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        for k in 0..valence {
            n1 = usize::from(
                *atom1
                    .neighbor
                    .get(k)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            if *rank2.get(n1).ok_or(SourceHeapError::PointerOutOfBounds)? != r_neigh1 {
                continue;
            }
            let i2 = if cumulene_len != 0 {
                let mut prev = i1;
                let mut next = n1;
                let mut len = 0_i32;
                while len < cumulene_len {
                    let next_atom = atoms.get(next).ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if next_atom.valence == 2 && next_atom.num_H == 0 {
                        let j = usize::from(usize::from(next_atom.neighbor[0]) == prev);
                        prev = next;
                        next = usize::from(next_atom.neighbor[j]);
                    } else {
                        break;
                    }
                    len += 1;
                }
                if len != cumulene_len
                    || *rank2.get(next).ok_or(SourceHeapError::PointerOutOfBounds)? != r2
                    || *rank2.get(prev).ok_or(SourceHeapError::PointerOutOfBounds)? != r_neigh2
                {
                    continue;
                }
                next
            } else {
                n1
            };

            let atom1 = atoms.get(i1).ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut k1 = 0_usize;
            let mut m = 0_usize;
            while k1 < MAX_NUM_STEREO_BONDS as usize {
                m = usize::from(atom1.stereo_bond_neighbor[k1]);
                if m == 0 || m - 1 == i2 {
                    break;
                }
                k1 += 1;
            }
            if m == 0 || m - 1 != i2 {
                return Ok(0);
            }
            let atom2 = atoms.get(i2).ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut k2 = 0_usize;
            m = 0;
            while k2 < MAX_NUM_STEREO_BONDS as usize {
                m = usize::from(atom2.stereo_bond_neighbor[k2]);
                if m == 0 || m - 1 == i1 {
                    break;
                }
                k2 += 1;
            }
            if m == 0 || m - 1 != i1 {
                return Ok(0);
            }
            if atom1.stereo_bond_parity[k1] != atom2.stereo_bond_parity[k2] {
                return Ok(-1);
            }
            if stereo_bond_parity != i32::from(atom1.stereo_bond_parity[k1]) {
                return Ok(0);
            }
            num_equal = num_equal.wrapping_add(1);
        }
    }
    Ok(num_equal)
}

#[allow(non_snake_case)]
pub(crate) fn CurTreeGetPos(cur_tree: Option<&CUR_TREE>) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1011 CurTreeGetPos
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int CurTreeGetPos( CUR_TREE *cur_tree )
    {
        if (cur_tree)
        {
            return cur_tree->cur_len;
        }

        return -1;
    }
    */
    // END INCHI C FUNCTION: CurTreeGetPos

    cur_tree.map_or(-1, |cur_tree| cur_tree.cur_len)
}

#[allow(non_snake_case)]
pub(crate) fn CurTreeRemoveLastRank(
    heap: &SourceHeap,
    cur_tree: Option<&mut CUR_TREE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1036 CurTreeRemoveLastRank
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int CurTreeRemoveLastRank( CUR_TREE *cur_tree )
    {
        if (cur_tree && cur_tree->cur_len > 0)
        {
            cur_tree->cur_len -= cur_tree->tree[cur_tree->cur_len - 1] + 1;
            if (cur_tree->cur_len >= 0)
            {
                return 0;
            }
        }

        return -1;
    }
    */
    // END INCHI C FUNCTION: CurTreeRemoveLastRank

    let Some(cur_tree) = cur_tree else {
        return Ok(-1);
    };
    if cur_tree.cur_len <= 0 {
        return Ok(-1);
    }
    let end =
        usize::try_from(cur_tree.cur_len - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let count = *heap
        .slice(cur_tree.tree.as_const())?
        .get(end)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    cur_tree.cur_len = cur_tree
        .cur_len
        .wrapping_sub(i32::from(count).wrapping_add(1));
    Ok(if cur_tree.cur_len >= 0 { 0 } else { -1 })
}

#[allow(non_snake_case)]
pub(crate) fn CurTreeIsLastAtomEqu(
    heap: &SourceHeap,
    cur_tree: Option<&CUR_TREE>,
    at_no: i32,
    nSymmStereo: SourceConstPointer<AT_NUMB>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1054 CurTreeIsLastAtomEqu
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int CurTreeIsLastAtomEqu( CUR_TREE *cur_tree, int at_no, AT_NUMB *nSymmStereo )
    {
        if (cur_tree && cur_tree->tree && nSymmStereo && cur_tree->cur_len > 1)
        {
            AT_NUMB nEq = nSymmStereo[at_no];
            int end = cur_tree->cur_len - 1;
            int len = cur_tree->tree[end] - 1;
            for (; len > 0; len--)
            {
                if (nSymmStereo[(int) cur_tree->tree[end - len]] == nEq)
                    return 1;
            }
            return 0;
        }

        return -1; /*  error */ /*   <BRKPT> */
    }
    */
    // END INCHI C FUNCTION: CurTreeIsLastAtomEqu

    let Some(cur_tree) = cur_tree else {
        return Ok(-1);
    };
    if cur_tree.tree.is_null() || nSymmStereo.is_null() || cur_tree.cur_len <= 1 {
        return Ok(-1);
    }
    let end =
        usize::try_from(cur_tree.cur_len - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let tree = heap.slice(cur_tree.tree.as_const())?;
    let len = i32::from(*tree.get(end).ok_or(SourceHeapError::PointerOutOfBounds)?) - 1;
    let at_no = usize::try_from(at_no).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let symm = heap.slice(nSymmStereo)?;
    let n_eq = *symm.get(at_no).ok_or(SourceHeapError::PointerOutOfBounds)?;
    for offset in (1..=len).rev() {
        let atom_index = end
            .checked_sub(usize::try_from(offset).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atom = usize::from(
            *tree
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        if *symm.get(atom).ok_or(SourceHeapError::PointerOutOfBounds)? == n_eq {
            return Ok(1);
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CurTreeAlloc(
    heap: &mut SourceHeap,
    cur_tree: Option<&mut CUR_TREE>,
    num_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:823 CurTreeAlloc
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int CurTreeAlloc( CUR_TREE *cur_tree, int num_atoms )
    {
        if (cur_tree)
        {
            if (cur_tree->tree && cur_tree->max_len > 0 && !( cur_tree->max_len % num_atoms ))
            {
                /*  do not reallocate */
                cur_tree->cur_len = 0;
                cur_tree->incr_len = num_atoms;
                memset( cur_tree->tree, 0, cur_tree->max_len * sizeof( cur_tree->tree[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                return 0; /*  ok */
            }
            inchi_free( cur_tree->tree );
            memset( cur_tree, 0, sizeof( *cur_tree ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            if ((cur_tree->tree = (AT_NUMB *) inchi_calloc( num_atoms, sizeof( cur_tree->tree[0] ) ))) /* djb-rwth: addressing LLVM warning */
            {
                cur_tree->incr_len =
                    cur_tree->max_len = num_atoms;
                return 0; /*  ok */
            }
        }

        return -1; /*  error */ /*   <BRKPT> */
    }
    */
    // END INCHI C FUNCTION: CurTreeAlloc

    let Some(cur_tree) = cur_tree else {
        return Ok(-1);
    };

    if !cur_tree.tree.is_null() && cur_tree.max_len > 0 {
        if num_atoms == 0 {
            return Err(SourceHeapError::UnsupportedSourceBehavior);
        }
        if cur_tree.max_len % num_atoms == 0 {
            cur_tree.cur_len = 0;
            cur_tree.incr_len = num_atoms;
            let max_len = usize::try_from(cur_tree.max_len)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let tree = heap.slice_mut(cur_tree.tree)?;
            let tree = tree
                .get_mut(..max_len)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            tree.fill(0);
            return Ok(0);
        }
    }

    inchi_free(heap, cur_tree.tree)?;
    *cur_tree = CUR_TREE::default();
    match inchi_calloc::<u16>(heap, num_atoms as u64, 2) {
        Ok(tree) => {
            cur_tree.tree = tree;
            cur_tree.max_len = num_atoms;
            cur_tree.incr_len = num_atoms;
            Ok(0)
        }
        Err(SourceHeapError::AllocationFailed)
        | Err(SourceHeapError::AllocationElementCountOutOfRange)
        | Err(SourceHeapError::AllocationSizeOverflow) => Ok(-1),
        Err(error) => Err(error),
    }
}

#[allow(non_snake_case)]
pub(crate) fn CurTreeFree(
    heap: &mut SourceHeap,
    cur_tree: Option<&mut CUR_TREE>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:871 CurTreeFree
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    void CurTreeFree( CUR_TREE *cur_tree )
    {
        if (cur_tree)
        {
            inchi_free( cur_tree->tree );
            memset( cur_tree, 0, sizeof( *cur_tree ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        }
    }
    */
    // END INCHI C FUNCTION: CurTreeFree

    if let Some(cur_tree) = cur_tree {
        inchi_free(heap, cur_tree.tree)?;
        *cur_tree = CUR_TREE::default();
    }
    Ok(())
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CompareLinCtStereo(
    heap: &SourceHeap,
    LinearCTStereoDble1: SourceConstPointer<AT_STEREO_DBLE>,
    nLenLinearCTStereoDble1: i32,
    LinearCTStereoCarb1: SourceConstPointer<AT_STEREO_CARB>,
    nLenLinearCTStereoCarb1: i32,
    LinearCTStereoDble2: SourceConstPointer<AT_STEREO_DBLE>,
    nLenLinearCTStereoDble2: i32,
    LinearCTStereoCarb2: SourceConstPointer<AT_STEREO_CARB>,
    nLenLinearCTStereoCarb2: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:262 CompareLinCtStereo
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int CompareLinCtStereo( AT_STEREO_DBLE *LinearCTStereoDble1,
                            int nLenLinearCTStereoDble1,
                            AT_STEREO_CARB *LinearCTStereoCarb1,
                            int nLenLinearCTStereoCarb1,
                            AT_STEREO_DBLE *LinearCTStereoDble2,
                            int nLenLinearCTStereoDble2,
                            AT_STEREO_CARB *LinearCTStereoCarb2,
                            int nLenLinearCTStereoCarb2 )
    {
        int ret;

        /* compare double bonds */
        ret = CompareLinCtStereoDble( LinearCTStereoDble1, nLenLinearCTStereoDble1,
                                       LinearCTStereoDble2, nLenLinearCTStereoDble2 );
        if (!ret)
        {
            ret = CompareLinCtStereoCarb( LinearCTStereoCarb1, nLenLinearCTStereoCarb1,

                                          LinearCTStereoCarb2, nLenLinearCTStereoCarb2 );
        }
        return ret;
    }
    */
    // END INCHI C FUNCTION: CompareLinCtStereo

    let ret = CompareLinCtStereoDble(
        heap,
        LinearCTStereoDble1,
        nLenLinearCTStereoDble1,
        LinearCTStereoDble2,
        nLenLinearCTStereoDble2,
    )?;
    if ret != 0 {
        Ok(ret)
    } else {
        CompareLinCtStereoCarb(
            heap,
            LinearCTStereoCarb1,
            nLenLinearCTStereoCarb1,
            LinearCTStereoCarb2,
            nLenLinearCTStereoCarb2,
        )
    }
}

#[allow(non_snake_case)]
pub(crate) fn CompareLinCtStereoCarb(
    heap: &SourceHeap,
    LinearCTStereoCarb1: SourceConstPointer<AT_STEREO_CARB>,
    nLenLinearCTStereoCarb1: i32,
    LinearCTStereoCarb2: SourceConstPointer<AT_STEREO_CARB>,
    nLenLinearCTStereoCarb2: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:223 CompareLinCtStereoCarb
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int CompareLinCtStereoCarb( AT_STEREO_CARB *LinearCTStereoCarb1,
                                int nLenLinearCTStereoCarb1,
                                AT_STEREO_CARB *LinearCTStereoCarb2,
                                int nLenLinearCTStereoCarb2 )
    {
        int i, num, ret = 0;

        /* compare stereocenters */
        if (LinearCTStereoCarb1 && LinearCTStereoCarb2)
        {
            num = inchi_min( nLenLinearCTStereoCarb1, nLenLinearCTStereoCarb2 );
            for (i = 0; i < num; i++)
            {
                if ((ret = (int) LinearCTStereoCarb1[i].at_num - (int) LinearCTStereoCarb2[i].at_num)) /* djb-rwth: addressing LLVM warning */
                    break;
                if ((ret = (int) LinearCTStereoCarb1[i].parity - (int) LinearCTStereoCarb2[i].parity)) /* djb-rwth: addressing LLVM warning */
                    break;
            }
            if (!ret)
            {
                ret = nLenLinearCTStereoCarb1 - nLenLinearCTStereoCarb2;
            }
        }
        else
            if (LinearCTStereoCarb1 && nLenLinearCTStereoCarb1 > 0)
            {
                ret = 1;
            }
            else
                if (LinearCTStereoCarb2 && nLenLinearCTStereoCarb2 > 0)
                {
                    ret = -1;
                }

        return ret;
    }
    */
    // END INCHI C FUNCTION: CompareLinCtStereoCarb

    if !LinearCTStereoCarb1.is_null() && !LinearCTStereoCarb2.is_null() {
        let num = nLenLinearCTStereoCarb1.min(nLenLinearCTStereoCarb2);
        if num > 0 {
            let num = usize::try_from(num).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let first = heap.slice(LinearCTStereoCarb1)?;
            let second = heap.slice(LinearCTStereoCarb2)?;
            if first.len() < num || second.len() < num {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            for i in 0..num {
                let mut ret = i32::from(first[i].at_num) - i32::from(second[i].at_num);
                if ret != 0 {
                    return Ok(ret);
                }
                ret = i32::from(first[i].parity) - i32::from(second[i].parity);
                if ret != 0 {
                    return Ok(ret);
                }
            }
        }
        Ok(nLenLinearCTStereoCarb1.wrapping_sub(nLenLinearCTStereoCarb2))
    } else if !LinearCTStereoCarb1.is_null() && nLenLinearCTStereoCarb1 > 0 {
        Ok(1)
    } else if !LinearCTStereoCarb2.is_null() && nLenLinearCTStereoCarb2 > 0 {
        Ok(-1)
    } else {
        Ok(0)
    }
}

#[allow(non_snake_case)]
pub(crate) fn CompareLinCtStereoDble(
    heap: &SourceHeap,
    LinearCTStereoDble1: SourceConstPointer<AT_STEREO_DBLE>,
    nLenLinearCTStereoDble1: i32,
    LinearCTStereoDble2: SourceConstPointer<AT_STEREO_DBLE>,
    nLenLinearCTStereoDble2: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:178 CompareLinCtStereoDble
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int CompareLinCtStereoDble( AT_STEREO_DBLE *LinearCTStereoDble1,
                                int nLenLinearCTStereoDble1,
                                AT_STEREO_DBLE *LinearCTStereoDble2,
                                int nLenLinearCTStereoDble2 )
    {
        int i, num, ret = 0;

        /* compare double bonds */
        if (LinearCTStereoDble1 && LinearCTStereoDble2)
        {
            num = inchi_min( nLenLinearCTStereoDble1, nLenLinearCTStereoDble2 );
            for (i = 0; i < num; i++)
            {
                if ((ret = (int) LinearCTStereoDble1[i].at_num1 - (int) LinearCTStereoDble2[i].at_num1)) /* djb-rwth: addressing LLVM warning */
                    break;
                if ((ret = (int) LinearCTStereoDble1[i].at_num2 - (int) LinearCTStereoDble2[i].at_num2)) /* djb-rwth: addressing LLVM warning */
                    break;
                if ((ret = (int) LinearCTStereoDble1[i].parity - (int) LinearCTStereoDble2[i].parity)) /* djb-rwth: addressing LLVM warning */
                    break;
            }
            if (!ret)
            {
                ret = nLenLinearCTStereoDble1 - nLenLinearCTStereoDble2;
            }
        }
        else
        {
            if (LinearCTStereoDble1 && nLenLinearCTStereoDble1 > 0)
            {
                ret = 1;
            }
            else
            {
                if (LinearCTStereoDble2 && nLenLinearCTStereoDble2 > 0)
                {
                    ret = -1;
                }
            }
        }

        return ret;
    }
    */
    // END INCHI C FUNCTION: CompareLinCtStereoDble

    if !LinearCTStereoDble1.is_null() && !LinearCTStereoDble2.is_null() {
        let num = nLenLinearCTStereoDble1.min(nLenLinearCTStereoDble2);
        if num > 0 {
            let num = usize::try_from(num).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let first = heap.slice(LinearCTStereoDble1)?;
            let second = heap.slice(LinearCTStereoDble2)?;
            if first.len() < num || second.len() < num {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            for i in 0..num {
                let mut ret = i32::from(first[i].at_num1) - i32::from(second[i].at_num1);
                if ret != 0 {
                    return Ok(ret);
                }
                ret = i32::from(first[i].at_num2) - i32::from(second[i].at_num2);
                if ret != 0 {
                    return Ok(ret);
                }
                ret = i32::from(first[i].parity) - i32::from(second[i].parity);
                if ret != 0 {
                    return Ok(ret);
                }
            }
        }
        Ok(nLenLinearCTStereoDble1.wrapping_sub(nLenLinearCTStereoDble2))
    } else if !LinearCTStereoDble1.is_null() && nLenLinearCTStereoDble1 > 0 {
        Ok(1)
    } else if !LinearCTStereoDble2.is_null() && nLenLinearCTStereoDble2 > 0 {
        Ok(-1)
    } else {
        Ok(0)
    }
}

#[allow(non_snake_case)]
pub(crate) fn CurTreeSetPos(cur_tree: Option<&mut CUR_TREE>, len: i32) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1023 CurTreeSetPos
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int CurTreeSetPos( CUR_TREE *cur_tree, int len )
    {
        if (cur_tree)
        {
            cur_tree->cur_len = len;
            return 0;
        }

        return -1;
    }
    */
    // END INCHI C FUNCTION: CurTreeSetPos

    if let Some(cur_tree) = cur_tree {
        cur_tree.cur_len = len;
        0
    } else {
        -1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_port__ichimap1__curtreeremoveiflastatom__line_993() {
        assert_eq!(
            CurTreeRemoveIfLastAtom(&mut SourceHeap::default(), None, 2),
            Ok(-1)
        );
        let mut heap = SourceHeap::default();
        let pointer = heap.allocate_model_storage(vec![7_u16, 1, 2, 3]).unwrap();
        let mut tree = CUR_TREE {
            tree: pointer,
            max_len: 4,
            cur_len: 4,
            incr_len: 4,
        };
        assert_eq!(
            CurTreeRemoveIfLastAtom(&mut heap, Some(&mut tree), 1),
            Ok(1)
        );
        assert_eq!(tree.cur_len, 4);
        assert_eq!(
            CurTreeRemoveIfLastAtom(&mut heap, Some(&mut tree), 2),
            Ok(0)
        );
        assert_eq!(tree.cur_len, 3);
        assert_eq!(heap.slice(tree.tree.as_const()).unwrap(), &[7, 1, 2, 3]);
        assert_eq!(
            CurTreeRemoveIfLastAtom(&mut heap, Some(&mut tree), 1),
            Ok(0)
        );
        assert_eq!(tree.cur_len, 2);
        assert_eq!(heap.slice(tree.tree.as_const()).unwrap(), &[7, 1, 2, 3]);
        assert_eq!(
            CurTreeRemoveIfLastAtom(&mut heap, Some(&mut tree), 7),
            Ok(-1)
        );

        let short_len = heap.allocate_model_storage(vec![9_u16, 4, 1]).unwrap();
        let mut short_tree = CUR_TREE {
            tree: short_len,
            max_len: 3,
            cur_len: 3,
            incr_len: 3,
        };
        assert_eq!(
            CurTreeRemoveIfLastAtom(&mut heap, Some(&mut short_tree), 4),
            Ok(1)
        );
        assert_eq!(short_tree.cur_len, 3);

        let mut null_tree = CUR_TREE {
            cur_len: 3,
            ..CUR_TREE::default()
        };
        assert_eq!(
            CurTreeRemoveIfLastAtom(&mut heap, Some(&mut null_tree), 0),
            Ok(-1)
        );
        let mut too_long = CUR_TREE {
            tree: pointer,
            max_len: 4,
            cur_len: 5,
            incr_len: 4,
        };
        assert_eq!(
            CurTreeRemoveIfLastAtom(&mut heap, Some(&mut too_long), 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(too_long.cur_len, 5);
    }

    #[test]
    fn source_port__ichimap1__curtreekeeplastatomsonly__line_960() {
        assert_eq!(
            CurTreeKeepLastAtomsOnly(&mut SourceHeap::default(), None, -1, 1),
            Ok(())
        );

        let mut heap = SourceHeap::default();
        let single_pointer = heap
            .allocate_model_storage(vec![10_u16, 1, 2, 3, 4])
            .unwrap();
        let mut single = CUR_TREE {
            tree: single_pointer,
            max_len: 5,
            cur_len: 5,
            incr_len: 5,
        };
        assert_eq!(
            CurTreeKeepLastAtomsOnly(&mut heap, Some(&mut single), -1, 1),
            Ok(())
        );
        assert_eq!(single.cur_len, 3);
        assert_eq!(
            heap.slice(single.tree.as_const()).unwrap(),
            &[10, 3, 2, 3, 4]
        );

        let multiple_pointer = heap
            .allocate_model_storage(vec![10_u16, 1, 2, 3, 4, 20, 8, 2])
            .unwrap();
        let mut multiple = CUR_TREE {
            tree: multiple_pointer,
            max_len: 8,
            cur_len: 8,
            incr_len: 8,
        };
        assert_eq!(
            CurTreeKeepLastAtomsOnly(&mut heap, Some(&mut multiple), -1, 1),
            Ok(())
        );
        assert_eq!(multiple.cur_len, 6);
        assert_eq!(
            heap.slice(multiple.tree.as_const()).unwrap(),
            &[10, 3, 2, 20, 8, 2, 8, 2]
        );

        let unchanged = heap.allocate_model_storage(vec![10_u16, 1, 2]).unwrap();
        let mut short = CUR_TREE {
            tree: unchanged,
            max_len: 3,
            cur_len: 3,
            incr_len: 3,
        };
        assert_eq!(
            CurTreeKeepLastAtomsOnly(&mut heap, Some(&mut short), -1, 1),
            Ok(())
        );
        assert_eq!(short.cur_len, 3);
        assert_eq!(heap.slice(short.tree.as_const()).unwrap(), &[10, 1, 2]);
        assert_eq!(
            CurTreeKeepLastAtomsOnly(&mut heap, Some(&mut short), 2, 1),
            Ok(())
        );

        let malformed_pointer = heap.allocate_model_storage(vec![7_u16, 9]).unwrap();
        let mut malformed = CUR_TREE {
            tree: malformed_pointer,
            max_len: 2,
            cur_len: 2,
            incr_len: 2,
        };
        assert_eq!(
            CurTreeKeepLastAtomsOnly(&mut heap, Some(&mut malformed), -1, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(malformed.cur_len, -5);

        let mut null_tree = CUR_TREE {
            cur_len: i32::MAX,
            ..CUR_TREE::default()
        };
        assert_eq!(
            CurTreeKeepLastAtomsOnly(&mut heap, Some(&mut null_tree), -1, 1),
            Ok(())
        );
        assert_eq!(null_tree.cur_len, i32::MAX);
    }

    #[test]
    fn source_port__ichimap1__curtreeaddatom__line_935() {
        assert_eq!(CurTreeAddAtom(&mut SourceHeap::default(), None, 3), Ok(-1));
        let mut heap = SourceHeap::default();
        let pointer = heap.allocate_model_storage(vec![7_u16, 1, 0, 0]).unwrap();
        let mut tree = CUR_TREE {
            tree: pointer,
            max_len: 4,
            cur_len: 2,
            incr_len: 4,
        };
        assert_eq!(CurTreeAddAtom(&mut heap, Some(&mut tree), -1), Ok(0));
        assert_eq!(tree.cur_len, 3);
        assert_eq!(
            heap.slice(tree.tree.as_const()).unwrap(),
            &[7, u16::MAX, 2, 0]
        );
        assert_eq!(CurTreeAddAtom(&mut heap, Some(&mut tree), 65_536), Ok(0));
        assert_eq!(tree.cur_len, 4);
        assert_eq!(
            heap.slice(tree.tree.as_const()).unwrap(),
            &[7, u16::MAX, 0, 3]
        );

        let old = tree.tree;
        assert_eq!(CurTreeAddAtom(&mut heap, Some(&mut tree), i32::MAX), Ok(0));
        assert_ne!(tree.tree, old);
        assert_eq!((tree.max_len, tree.cur_len, tree.incr_len), (8, 5, 4));
        assert_eq!(
            heap.slice(tree.tree.as_const()).unwrap(),
            &[7, u16::MAX, 0, u16::MAX, 4, 0, 0, 0]
        );

        let mut empty = CUR_TREE {
            tree: tree.tree,
            max_len: 8,
            cur_len: 0,
            incr_len: 4,
        };
        assert_eq!(CurTreeAddAtom(&mut heap, Some(&mut empty), 1), Ok(-1));
        assert_eq!(empty.cur_len, 0);

        let overflow_marker = heap
            .allocate_model_storage(vec![9_u16, u16::MAX, 0])
            .unwrap();
        let mut overflow_tree = CUR_TREE {
            tree: overflow_marker,
            max_len: 3,
            cur_len: 2,
            incr_len: 3,
        };
        assert_eq!(
            CurTreeAddAtom(&mut heap, Some(&mut overflow_tree), 4),
            Ok(0)
        );
        assert_eq!(
            heap.slice(overflow_tree.tree.as_const()).unwrap(),
            &[9, 4, 0]
        );

        let mut failure_heap = SourceHeap::default();
        let leaked = failure_heap.allocate_model_storage(vec![5_u16, 1]).unwrap();
        let mut failure = CUR_TREE {
            tree: leaked,
            max_len: 2,
            cur_len: 2,
            incr_len: 2,
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            CurTreeAddAtom(&mut failure_heap, Some(&mut failure), 3),
            Ok(-1)
        );
        assert!(failure.tree.is_null());
        assert_eq!(failure.cur_len, 2);
        assert_eq!(failure_heap.slice(leaked.as_const()).unwrap(), &[5, 1]);
    }

    #[test]
    fn source_port__ichimap1__curtreeremovelastrankifnoatoms__line_920() {
        assert_eq!(
            CurTreeRemoveLastRankIfNoAtoms(&SourceHeap::default(), None),
            Ok(-1)
        );
        let mut heap = SourceHeap::default();
        let pointer = heap.allocate_model_storage(vec![7_u16, 1, 8, 2]).unwrap();
        let mut tree = CUR_TREE {
            tree: pointer,
            max_len: 4,
            cur_len: 2,
            incr_len: 4,
        };
        assert_eq!(
            CurTreeRemoveLastRankIfNoAtoms(&heap, Some(&mut tree)),
            Ok(0)
        );
        assert_eq!(tree.cur_len, 0);

        tree.cur_len = 4;
        assert_eq!(
            CurTreeRemoveLastRankIfNoAtoms(&heap, Some(&mut tree)),
            Ok(1)
        );
        assert_eq!(tree.cur_len, 4);
        tree.cur_len = 1;
        assert_eq!(
            CurTreeRemoveLastRankIfNoAtoms(&heap, Some(&mut tree)),
            Ok(-1)
        );
        assert_eq!(tree.cur_len, 1);

        let mut null_tree = CUR_TREE {
            cur_len: 2,
            ..CUR_TREE::default()
        };
        assert_eq!(
            CurTreeRemoveLastRankIfNoAtoms(&heap, Some(&mut null_tree)),
            Ok(-1)
        );
        let mut too_long = CUR_TREE { cur_len: 5, ..tree };
        assert_eq!(
            CurTreeRemoveLastRankIfNoAtoms(&heap, Some(&mut too_long)),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(too_long.cur_len, 5);
    }

    #[test]
    fn source_port__ichimap1__curtreeislastrank__line_902() {
        assert_eq!(CurTreeIsLastRank(&SourceHeap::default(), None, 5), Ok(0));
        let mut heap = SourceHeap::default();
        let pointer = heap.allocate_model_storage(vec![5_u16, 20, 21, 3]).unwrap();
        let tree = CUR_TREE {
            tree: pointer,
            max_len: 4,
            cur_len: 4,
            incr_len: 4,
        };
        assert_eq!(CurTreeIsLastRank(&heap, Some(&tree), 5), Ok(1));
        assert_eq!(CurTreeIsLastRank(&heap, Some(&tree), 6), Ok(0));

        let empty = CUR_TREE {
            cur_len: 0,
            ..tree.clone()
        };
        assert_eq!(CurTreeIsLastRank(&heap, Some(&empty), 5), Ok(0));
        let negative = CUR_TREE {
            cur_len: -1,
            ..tree.clone()
        };
        assert_eq!(CurTreeIsLastRank(&heap, Some(&negative), 5), Ok(0));

        heap.slice_mut(pointer).unwrap()[3] = 5;
        assert_eq!(CurTreeIsLastRank(&heap, Some(&tree), 5), Ok(0));
        let too_long = CUR_TREE {
            cur_len: 5,
            ..tree.clone()
        };
        assert_eq!(
            CurTreeIsLastRank(&heap, Some(&too_long), 5),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let null_tree = CUR_TREE {
            cur_len: 1,
            ..CUR_TREE::default()
        };
        assert_eq!(
            CurTreeIsLastRank(&heap, Some(&null_tree), 5),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichimap1__curtreeaddrank__line_881() {
        assert_eq!(CurTreeAddRank(&mut SourceHeap::default(), None, 7), Ok(-1));

        let mut heap = SourceHeap::default();
        let pointer = heap.allocate_model_storage(vec![0_u16; 4]).unwrap();
        let mut tree = CUR_TREE {
            tree: pointer,
            max_len: 4,
            cur_len: 0,
            incr_len: 4,
        };
        assert_eq!(CurTreeAddRank(&mut heap, Some(&mut tree), u16::MAX), Ok(0));
        assert_eq!(tree.cur_len, 2);
        assert_eq!(
            heap.slice(tree.tree.as_const()).unwrap(),
            &[u16::MAX, 1, 0, 0]
        );
        assert_eq!(CurTreeAddRank(&mut heap, Some(&mut tree), 42), Ok(0));
        assert_eq!(tree.cur_len, 4);
        assert_eq!(
            heap.slice(tree.tree.as_const()).unwrap(),
            &[u16::MAX, 1, 42, 1]
        );

        let old = tree.tree;
        assert_eq!(CurTreeAddRank(&mut heap, Some(&mut tree), 77), Ok(0));
        assert_ne!(tree.tree, old);
        assert_eq!((tree.max_len, tree.cur_len, tree.incr_len), (8, 6, 4));
        assert_eq!(
            heap.slice(tree.tree.as_const()).unwrap(),
            &[u16::MAX, 1, 42, 1, 77, 1, 0, 0]
        );
        assert_eq!(
            heap.slice(old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let mut failure_heap = SourceHeap::default();
        let leaked = failure_heap.allocate_model_storage(vec![5_u16, 1]).unwrap();
        let mut failure = CUR_TREE {
            tree: leaked,
            max_len: 2,
            cur_len: 2,
            incr_len: 2,
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            CurTreeAddRank(&mut failure_heap, Some(&mut failure), 9),
            Ok(-1)
        );
        assert!(failure.tree.is_null());
        assert_eq!(
            (failure.max_len, failure.cur_len, failure.incr_len),
            (2, 2, 2)
        );
        assert_eq!(failure_heap.slice(leaked.as_const()).unwrap(), &[5, 1]);

        let mut malformed = CUR_TREE {
            tree: tree.tree,
            max_len: i32::MAX,
            cur_len: -1,
            incr_len: 1,
        };
        assert_eq!(
            CurTreeAddRank(&mut heap, Some(&mut malformed), 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(malformed.cur_len, -1);
    }

    #[test]
    fn source_port__ichimap1__curtreerealloc__line_850() {
        assert_eq!(CurTreeReAlloc(&mut SourceHeap::default(), None), Ok(-1));

        let mut heap = SourceHeap::default();
        let old = heap.allocate_model_storage(vec![4_u16, 5, 6]).unwrap();
        let mut tree = CUR_TREE {
            tree: old,
            max_len: 3,
            cur_len: 2,
            incr_len: 2,
        };
        assert_eq!(CurTreeReAlloc(&mut heap, Some(&mut tree)), Ok(0));
        assert_ne!(tree.tree, old);
        assert_eq!((tree.max_len, tree.cur_len, tree.incr_len), (5, 2, 2));
        assert_eq!(heap.slice(tree.tree.as_const()).unwrap(), &[4, 5, 0, 0, 0]);
        assert_eq!(
            heap.slice(old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        for invalid in [
            CUR_TREE::default(),
            CUR_TREE {
                tree: tree.tree,
                max_len: 0,
                cur_len: 0,
                incr_len: 1,
            },
            CUR_TREE {
                tree: tree.tree,
                max_len: 5,
                cur_len: 0,
                incr_len: 0,
            },
        ] {
            let mut invalid = invalid;
            let original = invalid.clone();
            assert_eq!(CurTreeReAlloc(&mut heap, Some(&mut invalid)), Ok(-1));
            assert_eq!(invalid, original);
        }

        let mut failure_heap = SourceHeap::default();
        let leaked = failure_heap.allocate_model_storage(vec![9_u16, 8]).unwrap();
        let mut failure = CUR_TREE {
            tree: leaked,
            max_len: 2,
            cur_len: 1,
            incr_len: 2,
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            CurTreeReAlloc(&mut failure_heap, Some(&mut failure)),
            Ok(-1)
        );
        assert!(failure.tree.is_null());
        assert_eq!(
            (failure.max_len, failure.cur_len, failure.incr_len),
            (2, 1, 2)
        );
        assert_eq!(failure_heap.slice(leaked.as_const()).unwrap(), &[9, 8]);

        let mut malformed_heap = SourceHeap::default();
        let malformed_old = malformed_heap.allocate_model_storage(vec![1_u16]).unwrap();
        let mut malformed = CUR_TREE {
            tree: malformed_old,
            max_len: 1,
            cur_len: 2,
            incr_len: 1,
        };
        assert_eq!(
            CurTreeReAlloc(&mut malformed_heap, Some(&mut malformed)),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_ne!(malformed.tree, malformed_old);
        assert_eq!(
            malformed_heap.slice(malformed_old.as_const()).unwrap(),
            &[1]
        );
        assert_eq!(
            malformed_heap.slice(malformed.tree.as_const()).unwrap(),
            &[0, 0]
        );
        assert_eq!(malformed.max_len, 1);
    }

    #[test]
    fn source_port__ichimap1__setuseatomforstereo__line_804() {
        let mut no_parity = sp_ATOM::default();
        no_parity.stereo_bond_neighbor = [1, 2, 3];
        let mut atom_mark = sp_ATOM::default();
        atom_mark.parity = -1;
        let mut two_bonds = sp_ATOM::default();
        two_bonds.parity = 1;
        two_bonds.stereo_bond_neighbor = [2, 3, 0];
        let mut three_bonds = sp_ATOM::default();
        three_bonds.parity = 2;
        three_bonds.stereo_bond_neighbor = [2, 3, 4];

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![no_parity, atom_mark, two_bonds, three_bonds])
            .unwrap();
        let used = heap.allocate_model_storage(vec![99_i8; 6]).unwrap();
        assert_eq!(
            SetUseAtomForStereo(&mut heap, used, atoms.as_const(), 4),
            Ok(())
        );
        assert_eq!(
            heap.slice(used.as_const()).unwrap(),
            &[0, STEREO_AT_MARK as i8, 2, 3, 99, 99]
        );

        assert_eq!(
            SetUseAtomForStereo(&mut heap, used, atoms.as_const(), 0),
            Ok(())
        );
        assert_eq!(
            heap.slice(used.as_const()).unwrap(),
            &[0, STEREO_AT_MARK as i8, 2, 3, 99, 99]
        );
        assert_eq!(
            SetUseAtomForStereo(&mut heap, used, atoms.as_const(), -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            SetUseAtomForStereo(&mut heap, used, atoms.as_const(), 5),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let short_used = heap.allocate_model_storage(vec![7_i8; 2]).unwrap();
        assert_eq!(
            SetUseAtomForStereo(&mut heap, short_used, atoms.as_const(), 4),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(short_used.as_const()).unwrap(), &[7, 7]);
        assert_eq!(
            SetUseAtomForStereo(&mut heap, SourceMutPointer::null(), atoms.as_const(), 0),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichimap1__comparelinctstereodoubletovalues__line_764() {
        let mut heap = SourceHeap::default();
        let value = heap
            .allocate_model_storage(vec![stereo_dble(10, 20, 3)])
            .unwrap();
        assert_eq!(
            CompareLinCtStereoDoubleToValues(&heap, value.as_const(), 10, 20, 3),
            Ok(0)
        );
        assert_eq!(
            CompareLinCtStereoDoubleToValues(&heap, value.as_const(), 9, u16::MAX, u8::MAX),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereoDoubleToValues(&heap, value.as_const(), 11, 0, 0),
            Ok(-1)
        );
        assert_eq!(
            CompareLinCtStereoDoubleToValues(&heap, value.as_const(), 10, 19, u8::MAX),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereoDoubleToValues(&heap, value.as_const(), 10, 21, 0),
            Ok(-1)
        );
        assert_eq!(
            CompareLinCtStereoDoubleToValues(&heap, value.as_const(), 10, 20, 2),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereoDoubleToValues(&heap, value.as_const(), 10, 20, 4),
            Ok(-1)
        );

        let extrema = heap
            .allocate_model_storage(vec![stereo_dble(u16::MAX, u16::MAX, u8::MAX)])
            .unwrap();
        assert_eq!(
            CompareLinCtStereoDoubleToValues(
                &heap,
                extrema.as_const(),
                u16::MAX,
                u16::MAX,
                u8::MAX
            ),
            Ok(0)
        );
        assert_eq!(
            CompareLinCtStereoDoubleToValues(&heap, SourceConstPointer::null(), 0, 0, 0),
            Err(SourceHeapError::NullPointer)
        );
        let empty = heap
            .allocate_model_storage(Vec::<AT_STEREO_DBLE>::new())
            .unwrap();
        assert_eq!(
            CompareLinCtStereoDoubleToValues(&heap, empty.as_const(), 0, 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichimap1__comparelinctstereoatomtovalues__line_287() {
        let mut heap = SourceHeap::default();
        let value = heap
            .allocate_model_storage(vec![stereo_carb(10, 3)])
            .unwrap();
        assert_eq!(
            CompareLinCtStereoAtomToValues(&heap, value.as_const(), 10, 3),
            Ok(0)
        );
        assert_eq!(
            CompareLinCtStereoAtomToValues(&heap, value.as_const(), 9, u8::MAX),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereoAtomToValues(&heap, value.as_const(), 11, 0),
            Ok(-1)
        );
        assert_eq!(
            CompareLinCtStereoAtomToValues(&heap, value.as_const(), 10, 2),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereoAtomToValues(&heap, value.as_const(), 10, 4),
            Ok(-1)
        );

        let extrema = heap
            .allocate_model_storage(vec![stereo_carb(u16::MAX, u8::MAX)])
            .unwrap();
        assert_eq!(
            CompareLinCtStereoAtomToValues(&heap, extrema.as_const(), u16::MAX, u8::MAX),
            Ok(0)
        );
        assert_eq!(
            CompareLinCtStereoAtomToValues(&heap, SourceConstPointer::null(), 0, 0),
            Err(SourceHeapError::NullPointer)
        );
        let empty = heap
            .allocate_model_storage(Vec::<AT_STEREO_CARB>::new())
            .unwrap();
        assert_eq!(
            CompareLinCtStereoAtomToValues(&heap, empty.as_const(), 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichimap1__buniqueatnbrfrommappingrank__line_316() {
        let mut heap = SourceHeap::default();
        let ranks = heap.allocate_model_storage(vec![1_u16, 3, 3, 4]).unwrap();
        let order = heap.allocate_model_storage(vec![0_u16, 1, 2, 3]).unwrap();
        let stack = heap.allocate_model_storage(vec![ranks, order]).unwrap();
        let mut atom = 99_u16;
        assert_eq!(
            bUniqueAtNbrFromMappingRank(&heap, stack, 1, &mut atom),
            Ok(1)
        );
        assert_eq!(atom, 0);
        atom = 99;
        assert_eq!(
            bUniqueAtNbrFromMappingRank(&heap, stack, 2, &mut atom),
            Ok(0)
        );
        assert_eq!(atom, 99);
        assert_eq!(
            bUniqueAtNbrFromMappingRank(&heap, stack, 3, &mut atom),
            Ok(0)
        );
        assert_eq!(atom, 99);
        assert_eq!(
            bUniqueAtNbrFromMappingRank(&heap, stack, 4, &mut atom),
            Ok(1)
        );
        assert_eq!(atom, 3);

        atom = 77;
        assert_eq!(
            bUniqueAtNbrFromMappingRank(&heap, stack, 0, &mut atom),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(atom, 77);
        assert_eq!(
            bUniqueAtNbrFromMappingRank(&heap, stack, 5, &mut atom),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(atom, 77);
        assert_eq!(
            bUniqueAtNbrFromMappingRank(&heap, SourceMutPointer::null(), 1, &mut atom),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(atom, 77);

        let short_stack = heap.allocate_model_storage(vec![ranks]).unwrap();
        assert_eq!(
            bUniqueAtNbrFromMappingRank(&heap, short_stack, 1, &mut atom),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(atom, 77);
    }

    #[test]
    fn source_port__ichimap1__ngetmcr__line_336() {
        let mut heap = SourceHeap::default();
        let classes = heap
            .allocate_model_storage(vec![0_u16, 0, 1, 2, 4, 4])
            .unwrap();
        assert_eq!(nGetMcr(&mut heap, classes, 0), Ok(0));
        assert_eq!(heap.slice(classes.as_const()).unwrap(), &[0, 0, 1, 2, 4, 4]);
        assert_eq!(nGetMcr(&mut heap, classes, 3), Ok(0));
        assert_eq!(heap.slice(classes.as_const()).unwrap(), &[0, 0, 0, 0, 4, 4]);
        assert_eq!(nGetMcr(&mut heap, classes, 5), Ok(4));
        assert_eq!(heap.slice(classes.as_const()).unwrap(), &[0, 0, 0, 0, 4, 4]);
        assert_eq!(
            nGetMcr(&mut heap, classes, 6),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            nGetMcr(&mut heap, SourceMutPointer::null(), 0),
            Err(SourceHeapError::NullPointer)
        );

        let dangling = heap.allocate_model_storage(vec![2_u16]).unwrap();
        assert_eq!(
            nGetMcr(&mut heap, dangling, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(dangling.as_const()).unwrap(), &[2]);
    }

    #[test]
    fn source_port__ichimap1__njoin2mcrs__line_366() {
        let mut heap = SourceHeap::default();
        let classes = heap.allocate_model_storage(vec![0_u16, 0, 2, 2]).unwrap();
        assert_eq!(nJoin2Mcrs(&mut heap, classes, 1, 3), Ok(1));
        assert_eq!(heap.slice(classes.as_const()).unwrap(), &[0, 0, 0, 2]);
        assert_eq!(nJoin2Mcrs(&mut heap, classes, 3, 1), Ok(0));
        assert_eq!(heap.slice(classes.as_const()).unwrap(), &[0, 0, 0, 0]);

        let reverse = heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
        assert_eq!(nJoin2Mcrs(&mut heap, reverse, 1, 0), Ok(1));
        assert_eq!(heap.slice(reverse.as_const()).unwrap(), &[0, 0]);

        let partial = heap.allocate_model_storage(vec![0_u16, 0, 1, 2]).unwrap();
        assert_eq!(
            nJoin2Mcrs(&mut heap, partial, 3, 9),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(partial.as_const()).unwrap(), &[0, 0, 0, 0]);
        assert_eq!(
            nJoin2Mcrs(&mut heap, SourceMutPointer::null(), 0, 0),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichimap1__all_sc_same__line_53() {
        fn run(atoms: Vec<sp_ATOM>, ranks_to: Vec<AT_RANK>) -> Result<i32, SourceHeapError> {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let ranks_from = heap.allocate_model_storage(vec![3_u16]).unwrap();
            let stack_from = heap.allocate_model_storage(vec![ranks_from]).unwrap();
            let ranks_to = heap.allocate_model_storage(ranks_to).unwrap();
            let order_to = heap.allocate_model_storage(vec![0_u16, 1, 2]).unwrap();
            let stack_to = heap
                .allocate_model_storage(vec![ranks_to, order_to])
                .unwrap();
            let canon_order = heap.allocate_model_storage(vec![0_u16]).unwrap();
            All_SC_Same(
                &heap,
                1,
                stack_from,
                stack_to,
                canon_order.as_const(),
                atoms.as_const(),
            )
        }

        let mut atoms = vec![sp_ATOM::default(); 3];
        for atom in &mut atoms {
            atom.stereo_atom_parity = 1;
        }
        assert_eq!(run(atoms.clone(), vec![3, 3, 3]), Ok(3));

        let mut masked = atoms.clone();
        for atom in &mut masked {
            atom.stereo_atom_parity = 9;
        }
        assert_eq!(run(masked, vec![3, 3, 3]), Ok(3));

        let mut different = atoms.clone();
        different[1].stereo_atom_parity = 2;
        assert_eq!(run(different, vec![3, 3, 3]), Ok(0));

        let mut stereo_bond = atoms.clone();
        stereo_bond[1].stereo_bond_neighbor[0] = 1;
        assert_eq!(run(stereo_bond, vec![3, 3, 3]), Ok(0));

        let mut unknown_first = atoms.clone();
        unknown_first[2].stereo_atom_parity = 5;
        assert_eq!(run(unknown_first, vec![3, 3, 3]), Ok(0));
        assert_eq!(run(atoms, vec![3, 2, 3]), Ok(1));

        let mut invalid_heap = SourceHeap::default();
        let ranks = invalid_heap.allocate_model_storage(vec![1_u16]).unwrap();
        let from_stack = invalid_heap.allocate_model_storage(vec![ranks]).unwrap();
        let order = invalid_heap.allocate_model_storage(vec![0_u16]).unwrap();
        let to_stack = invalid_heap
            .allocate_model_storage(vec![ranks, order])
            .unwrap();
        let canon = invalid_heap.allocate_model_storage(vec![0_u16]).unwrap();
        let atom = invalid_heap
            .allocate_model_storage(vec![sp_ATOM::default()])
            .unwrap();
        assert_eq!(
            All_SC_Same(
                &invalid_heap,
                0,
                from_stack,
                to_stack,
                canon.as_const(),
                atom.as_const(),
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            All_SC_Same(
                &invalid_heap,
                1,
                SourceMutPointer::null(),
                to_stack,
                canon.as_const(),
                atom.as_const(),
            ),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichimap1__next_sc_at_canonrank2__line_99() {
        let mut heap = SourceHeap::default();
        let ranks_from = heap.allocate_model_storage(vec![1_u16, 2, 3, 4]).unwrap();
        let stack_from = heap.allocate_model_storage(vec![ranks_from]).unwrap();
        let ranks_to = heap.allocate_model_storage(vec![1_u16, 2, 3, 4]).unwrap();
        let order_to = heap.allocate_model_storage(vec![0_u16, 1, 2, 3]).unwrap();
        let stack_to = heap
            .allocate_model_storage(vec![ranks_to, order_to])
            .unwrap();
        let canon_order = heap.allocate_model_storage(vec![0_u16, 1, 2, 3]).unwrap();
        let used = heap
            .allocate_model_storage(vec![0_i8, 0, STEREO_AT_MARK as i8, 0])
            .unwrap();

        let mut canon_rank = 0_u16;
        let mut canon_min = 0_u16;
        let mut first = 1_i32;
        assert_eq!(
            Next_SC_At_CanonRank2(
                &heap,
                &mut canon_rank,
                &mut canon_min,
                &mut first,
                used.as_const(),
                stack_from,
                stack_to,
                canon_order.as_const(),
                4,
            ),
            Ok(1)
        );
        assert_eq!((canon_rank, canon_min, first), (3, 3, 0));
        assert_eq!(
            Next_SC_At_CanonRank2(
                &heap,
                &mut canon_rank,
                &mut canon_min,
                &mut first,
                used.as_const(),
                stack_from,
                stack_to,
                canon_order.as_const(),
                4,
            ),
            Ok(0)
        );
        assert_eq!((canon_rank, canon_min, first), (3, 3, 0));

        canon_rank = 1;
        first = 9;
        assert_eq!(
            Next_SC_At_CanonRank2(
                &heap,
                &mut canon_rank,
                &mut canon_min,
                &mut first,
                used.as_const(),
                stack_from,
                stack_to,
                canon_order.as_const(),
                4,
            ),
            Ok(1)
        );
        assert_eq!((canon_rank, canon_min, first), (3, 3, 0));

        let tie_from_rank = heap.allocate_model_storage(vec![2_u16, 2]).unwrap();
        let tie_from_stack = heap.allocate_model_storage(vec![tie_from_rank]).unwrap();
        let tie_to_rank = heap.allocate_model_storage(vec![2_u16, 2]).unwrap();
        let tie_to_order = heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
        let tie_to_stack = heap
            .allocate_model_storage(vec![tie_to_rank, tie_to_order])
            .unwrap();
        let tie_canon = heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
        let tie_used = heap
            .allocate_model_storage(vec![STEREO_AT_MARK as i8, 0])
            .unwrap();
        canon_rank = 0;
        canon_min = 0;
        first = 1;
        assert_eq!(
            Next_SC_At_CanonRank2(
                &heap,
                &mut canon_rank,
                &mut canon_min,
                &mut first,
                tie_used.as_const(),
                tie_from_stack,
                tie_to_stack,
                tie_canon.as_const(),
                2,
            ),
            Ok(1)
        );
        assert_eq!((canon_rank, canon_min, first), (1, 1, 0));

        canon_rank = 0;
        canon_min = 0;
        first = 1;
        assert_eq!(
            Next_SC_At_CanonRank2(
                &heap,
                &mut canon_rank,
                &mut canon_min,
                &mut first,
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                -1,
            ),
            Ok(0)
        );
        assert_eq!((canon_rank, canon_min, first), (0, 0, 1));

        canon_rank = u16::MAX;
        assert_eq!(
            Next_SC_At_CanonRank2(
                &heap,
                &mut canon_rank,
                &mut canon_min,
                &mut first,
                used.as_const(),
                stack_from,
                stack_to,
                canon_order.as_const(),
                4,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(canon_rank, u16::MAX);
    }

    fn next_parity(
        mut parity: i32,
        mut calc: i32,
        counts: [i32; 5],
        unknown: i32,
    ) -> (i32, i32, i32) {
        let status = NextStereoParity2Test(
            &mut parity,
            &mut calc,
            counts[0],
            counts[1],
            counts[2],
            counts[3],
            counts[4],
            unknown,
        );
        (status, parity, calc)
    }

    #[test]
    fn source_port__ichimap1__nextstereoparity2test__line_659() {
        let best = BEST_PARITY as i32;
        let worse = WORSE_PARITY as i32;
        let unknown = AB_PARITY_UNKN as i32;
        let undefined = AB_PARITY_UNDF as i32;

        assert_eq!(next_parity(best, 0, [1; 5], unknown), (0, worse, 0));
        assert_eq!(next_parity(best, best, [1; 5], unknown), (0, best, worse));
        assert_eq!(next_parity(best, worse, [1; 5], unknown), (0, worse, worse));
        assert_eq!(next_parity(best, 99, [1; 5], unknown), (0, best, 99));

        assert_eq!(next_parity(worse, 0, [1; 5], unknown), (0, unknown, 0));
        assert_eq!(
            next_parity(worse, best, [1; 5], unknown),
            (CT_STEREOCOUNT_ERR, worse, best)
        );
        assert_eq!(next_parity(worse, worse, [1; 5], unknown), (0, worse, 0));
        assert_eq!(next_parity(worse, 99, [1; 5], unknown), (0, worse, 99));

        assert_eq!(next_parity(unknown, 0, [1; 5], unknown), (0, undefined, 0));
        assert_eq!(
            next_parity(unknown, best, [1; 5], unknown),
            (CT_STEREOCOUNT_ERR, unknown, best)
        );
        assert_eq!(
            next_parity(undefined, 0, [1; 5], unknown),
            (1, undefined, 0)
        );
        assert_eq!(
            next_parity(undefined, worse, [1; 5], unknown),
            (CT_STEREOCOUNT_ERR, undefined, worse)
        );
        assert_eq!(next_parity(99, 88, [0; 5], unknown), (0, 99, 88));

        assert_eq!(next_parity(best, best, [0; 5], unknown), (1, undefined, 0));
        assert_eq!(
            next_parity(best, 0, [1, 0, 0, 0, 1], unknown),
            (1, undefined, 0)
        );
        assert_eq!(
            next_parity(best, worse, [1, 0, 1, 1, 0], unknown),
            (0, unknown, 0)
        );
        assert_eq!(
            next_parity(unknown, 0, [1, 1, 1, 0, 1], unknown),
            (1, undefined, 0)
        );
        assert_eq!(next_parity(worse, 0, [0; 5], 77), (0, 77, 0));
    }

    #[allow(clippy::too_many_arguments)]
    fn run_next_sb(
        atoms: Vec<sp_ATOM>,
        used: Vec<i8>,
        ranks: Vec<AT_RANK>,
        order: Vec<AT_RANK>,
        canon_from: Vec<AT_RANK>,
        atom_from_canon: Vec<AT_RANK>,
        mut cr1: AT_RANK,
        mut cr2: AT_RANK,
        mut cr1_min: AT_RANK,
        mut cr2_min: AT_RANK,
        mut first_time: i32,
        allene: i32,
    ) -> Result<(i32, AT_RANK, AT_RANK, AT_RANK, AT_RANK, i32), SourceHeapError> {
        let num_atoms = i32::try_from(atoms.len()).unwrap();
        let mut heap = SourceHeap::default();
        let rank1 = heap.allocate_model_storage(ranks.clone()).unwrap();
        let order1 = heap.allocate_model_storage(order.clone()).unwrap();
        let rank2 = heap.allocate_model_storage(ranks).unwrap();
        let order2 = heap.allocate_model_storage(order).unwrap();
        let stack1 = heap.allocate_model_storage(vec![rank1, order1]).unwrap();
        let stack2 = heap.allocate_model_storage(vec![rank2, order2]).unwrap();
        let canon_from = heap.allocate_model_storage(canon_from).unwrap();
        let atom_from_canon = heap.allocate_model_storage(atom_from_canon).unwrap();
        let atoms = heap.allocate_model_storage(atoms).unwrap();
        let used = heap.allocate_model_storage(used).unwrap();
        let status = Next_SB_At_CanonRanks2(
            &heap,
            &mut cr1,
            &mut cr2,
            &mut cr1_min,
            &mut cr2_min,
            &mut first_time,
            used.as_const(),
            stack1,
            stack2,
            canon_from.as_const(),
            atom_from_canon.as_const(),
            atoms.as_const(),
            num_atoms,
            allene,
        )?;
        Ok((status, cr1, cr2, cr1_min, cr2_min, first_time))
    }

    #[test]
    fn source_port__ichimap1__next_sb_at_canonranks2__line_518() {
        let pair = direct_stereo_pair(1, 1);
        assert_eq!(
            run_next_sb(
                pair.clone(),
                vec![1, 1],
                vec![1, 2],
                vec![0, 1],
                vec![1, 2],
                vec![0, 1],
                0,
                0,
                0,
                0,
                1,
                0,
            ),
            Ok((1, 2, 1, 2, 1, 0))
        );
        assert_eq!(
            run_next_sb(
                pair.clone(),
                vec![1, 1],
                vec![1, 2],
                vec![0, 1],
                vec![1, 2],
                vec![0, 1],
                2,
                1,
                2,
                1,
                0,
                0,
            ),
            Ok((0, 2, 1, 2, 1, 0))
        );
        assert_eq!(
            run_next_sb(
                pair.clone(),
                vec![1, 1],
                vec![1, 2],
                vec![0, 1],
                vec![1, 2],
                vec![0, 1],
                0,
                0,
                2,
                1,
                1,
                0,
            ),
            Ok((0, 0, 0, 2, 1, 1))
        );
        assert_eq!(
            run_next_sb(
                pair.clone(),
                vec![1, 1],
                vec![1, 2],
                vec![0, 1],
                vec![1, 2],
                vec![0, 1],
                0,
                0,
                0,
                0,
                1,
                1,
            ),
            Ok((0, 0, 0, 0, 0, 1))
        );
        assert_eq!(
            run_next_sb(
                pair,
                vec![STEREO_AT_MARK as i8, STEREO_AT_MARK as i8],
                vec![1, 2],
                vec![0, 1],
                vec![1, 2],
                vec![0, 1],
                0,
                0,
                0,
                0,
                1,
                0,
            ),
            Ok((0, 0, 0, 0, 0, 1))
        );

        let mut start = sp_ATOM::default();
        start.valence = 1;
        start.neighbor[0] = 1;
        start.stereo_bond_neighbor[0] = 3;
        start.stereo_bond_parity[0] = 9;
        let mut middle = sp_ATOM::default();
        middle.valence = 2;
        middle.neighbor[0] = 0;
        middle.neighbor[1] = 2;
        let mut end = sp_ATOM::default();
        end.valence = 1;
        end.neighbor[0] = 1;
        end.stereo_bond_neighbor[0] = 1;
        end.stereo_bond_parity[0] = 9;
        assert_eq!(
            run_next_sb(
                vec![start, middle, end],
                vec![1, 1, 1],
                vec![1, 2, 3],
                vec![0, 1, 2],
                vec![1, 2, 3],
                vec![0, 1, 2],
                0,
                0,
                0,
                0,
                1,
                1,
            ),
            Ok((1, 3, 1, 3, 1, 0))
        );
    }

    fn run_all_sb_same(
        atoms: Vec<sp_ATOM>,
        ranks1: Vec<AT_RANK>,
        ranks2: Vec<AT_RANK>,
        order2: Vec<AT_RANK>,
        canon_from: Vec<AT_RANK>,
        canon_rank1: AT_RANK,
        canon_rank2: AT_RANK,
    ) -> Result<i32, SourceHeapError> {
        let mut heap = SourceHeap::default();
        let ranks1 = heap.allocate_model_storage(ranks1).unwrap();
        let ranks2 = heap.allocate_model_storage(ranks2).unwrap();
        let order2 = heap.allocate_model_storage(order2).unwrap();
        let unused_order1 = heap.allocate_model_storage(Vec::<AT_RANK>::new()).unwrap();
        let stack1 = heap
            .allocate_model_storage(vec![ranks1, unused_order1])
            .unwrap();
        let stack2 = heap.allocate_model_storage(vec![ranks2, order2]).unwrap();
        let canon_from = heap.allocate_model_storage(canon_from).unwrap();
        let atoms = heap.allocate_model_storage(atoms).unwrap();
        All_SB_Same(
            &heap,
            canon_rank1,
            canon_rank2,
            stack1,
            stack2,
            canon_from.as_const(),
            atoms.as_const(),
        )
    }

    fn direct_stereo_pair(parity1: i8, parity2: i8) -> Vec<sp_ATOM> {
        let mut first = sp_ATOM::default();
        first.valence = 1;
        first.neighbor[0] = 1;
        first.stereo_bond_neighbor[0] = 2;
        first.stereo_bond_ord[0] = 0;
        first.stereo_bond_parity[0] = parity1;

        let mut second = sp_ATOM::default();
        second.valence = 1;
        second.neighbor[0] = 0;
        second.stereo_bond_neighbor[0] = 1;
        second.stereo_bond_ord[0] = 0;
        second.stereo_bond_parity[0] = parity2;
        vec![first, second]
    }

    #[test]
    fn source_port__ichimap1__all_sb_same__line_393() {
        assert_eq!(
            run_all_sb_same(
                direct_stereo_pair(1, 1),
                vec![1, 2],
                vec![1, 2],
                vec![0, 1],
                vec![0, 1],
                1,
                2,
            ),
            Ok(1)
        );

        let mut no_mapping = direct_stereo_pair(1, 1);
        no_mapping[0].stereo_bond_neighbor[0] = 0;
        assert_eq!(
            run_all_sb_same(
                no_mapping,
                vec![1, 2],
                vec![1, 2],
                vec![0, 1],
                vec![0, 1],
                1,
                2,
            ),
            Ok(-1)
        );

        let mut no_reverse = direct_stereo_pair(1, 1);
        no_reverse[1].stereo_bond_neighbor[0] = 0;
        assert_eq!(
            run_all_sb_same(
                no_reverse,
                vec![1, 2],
                vec![1, 2],
                vec![0, 1],
                vec![0, 1],
                1,
                2,
            ),
            Ok(-1)
        );

        assert_eq!(
            run_all_sb_same(
                direct_stereo_pair(0, 0),
                vec![1, 2],
                vec![1, 2],
                vec![0, 1],
                vec![0, 1],
                1,
                2,
            ),
            Ok(0)
        );
        assert_eq!(
            run_all_sb_same(
                direct_stereo_pair(1, 2),
                vec![1, 2],
                vec![1, 2],
                vec![0, 1],
                vec![0, 1],
                1,
                2,
            ),
            Ok(-1)
        );

        let mut non_stereo_mapping = direct_stereo_pair(1, 1);
        let mut equivalent = sp_ATOM::default();
        equivalent.valence = 1;
        equivalent.neighbor[0] = 1;
        non_stereo_mapping.push(equivalent);
        assert_eq!(
            run_all_sb_same(
                non_stereo_mapping,
                vec![2, 3, 2],
                vec![2, 3, 2],
                vec![0, 2, 1],
                vec![0, 1, 2],
                1,
                2,
            ),
            Ok(0)
        );

        let mut allene_start = sp_ATOM::default();
        allene_start.valence = 1;
        allene_start.neighbor[0] = 1;
        allene_start.stereo_bond_neighbor[0] = 3;
        allene_start.stereo_bond_ord[0] = 0;
        allene_start.stereo_bond_parity[0] = 9;
        let mut allene_middle = sp_ATOM::default();
        allene_middle.valence = 2;
        allene_middle.neighbor[0] = 0;
        allene_middle.neighbor[1] = 2;
        let mut allene_end = sp_ATOM::default();
        allene_end.valence = 1;
        allene_end.neighbor[0] = 1;
        allene_end.stereo_bond_neighbor[0] = 1;
        allene_end.stereo_bond_ord[0] = 0;
        allene_end.stereo_bond_parity[0] = 9;
        let allene = vec![allene_start, allene_middle, allene_end];
        assert_eq!(
            run_all_sb_same(
                allene.clone(),
                vec![1, 2, 3],
                vec![1, 2, 3],
                vec![0, 1, 2],
                vec![0, 2, 1],
                1,
                2,
            ),
            Ok(1)
        );
        let mut broken_allene = allene;
        broken_allene[1].num_H = 1;
        assert_eq!(
            run_all_sb_same(
                broken_allene,
                vec![1, 2, 3],
                vec![1, 2, 3],
                vec![0, 1, 2],
                vec![0, 2, 1],
                1,
                2,
            ),
            Ok(0)
        );

        assert_eq!(
            All_SB_Same(
                &SourceHeap::default(),
                1,
                2,
                ppAT_RANK::null(),
                ppAT_RANK::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
            ),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichimap1__curtreegetpos__line_1011() {
        assert_eq!(CurTreeGetPos(None), -1);
        for value in [i32::MIN, -1, 0, 1, i32::MAX] {
            let tree = CUR_TREE {
                max_len: 91,
                cur_len: value,
                incr_len: 92,
                ..CUR_TREE::default()
            };
            assert_eq!(CurTreeGetPos(Some(&tree)), value);
        }
    }

    #[test]
    fn source_port__ichimap1__curtreeremovelastrank__line_1036() {
        assert_eq!(CurTreeRemoveLastRank(&SourceHeap::default(), None), Ok(-1));

        let mut heap = SourceHeap::default();
        let pointer = heap.allocate_model_storage(vec![10_u16, 11, 2]).unwrap();
        let mut tree = CUR_TREE {
            tree: pointer,
            max_len: 3,
            cur_len: 3,
            incr_len: 3,
        };
        assert_eq!(CurTreeRemoveLastRank(&heap, Some(&mut tree)), Ok(0));
        assert_eq!(tree.cur_len, 0);

        assert_eq!(CurTreeRemoveLastRank(&heap, Some(&mut tree)), Ok(-1));
        assert_eq!(tree.cur_len, 0);
        tree.cur_len = -1;
        assert_eq!(CurTreeRemoveLastRank(&heap, Some(&mut tree)), Ok(-1));
        assert_eq!(tree.cur_len, -1);

        heap.slice_mut(pointer).unwrap()[0] = 1;
        tree.cur_len = 1;
        assert_eq!(CurTreeRemoveLastRank(&heap, Some(&mut tree)), Ok(-1));
        assert_eq!(tree.cur_len, -1);

        tree.cur_len = 4;
        assert_eq!(
            CurTreeRemoveLastRank(&heap, Some(&mut tree)),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(tree.cur_len, 4);

        let mut null_tree = CUR_TREE {
            cur_len: 1,
            ..CUR_TREE::default()
        };
        assert_eq!(
            CurTreeRemoveLastRank(&heap, Some(&mut null_tree)),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(null_tree.cur_len, 1);
    }

    #[test]
    fn source_port__ichimap1__curtreeislastatomequ__line_1054() {
        let mut heap = SourceHeap::default();
        let tree_pointer = heap
            .allocate_model_storage(vec![99_u16, 0, 1, 2, 3])
            .unwrap();
        let symm = heap
            .allocate_model_storage(vec![5_u16, 7, 8, 7, 9])
            .unwrap();
        let tree = CUR_TREE {
            tree: tree_pointer,
            max_len: 5,
            cur_len: 5,
            incr_len: 5,
        };
        assert_eq!(
            CurTreeIsLastAtomEqu(&heap, Some(&tree), 3, symm.as_const()),
            Ok(1)
        );
        assert_eq!(
            CurTreeIsLastAtomEqu(&heap, Some(&tree), 2, symm.as_const()),
            Ok(1)
        );
        assert_eq!(
            CurTreeIsLastAtomEqu(&heap, Some(&tree), 4, symm.as_const()),
            Ok(0)
        );

        let no_tree = CUR_TREE {
            cur_len: 2,
            ..CUR_TREE::default()
        };
        assert_eq!(
            CurTreeIsLastAtomEqu(&heap, Some(&no_tree), 0, symm.as_const()),
            Ok(-1)
        );
        assert_eq!(
            CurTreeIsLastAtomEqu(&heap, None, 0, symm.as_const()),
            Ok(-1)
        );
        assert_eq!(
            CurTreeIsLastAtomEqu(&heap, Some(&tree), 0, SourceConstPointer::null()),
            Ok(-1)
        );

        let short = CUR_TREE {
            tree: tree_pointer,
            cur_len: 1,
            ..tree.clone()
        };
        assert_eq!(
            CurTreeIsLastAtomEqu(&heap, Some(&short), 0, symm.as_const()),
            Ok(-1)
        );

        let malformed = heap.allocate_model_storage(vec![0_u16, 0, 0]).unwrap();
        let malformed_tree = CUR_TREE {
            tree: malformed,
            max_len: 3,
            cur_len: 3,
            incr_len: 3,
        };
        assert_eq!(
            CurTreeIsLastAtomEqu(&heap, Some(&malformed_tree), 0, symm.as_const()),
            Ok(0)
        );

        assert_eq!(
            CurTreeIsLastAtomEqu(&heap, Some(&tree), -1, symm.as_const()),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let out_of_bounds = CUR_TREE { cur_len: 6, ..tree };
        assert_eq!(
            CurTreeIsLastAtomEqu(&heap, Some(&out_of_bounds), 0, symm.as_const()),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichimap1__curtreealloc__line_823() {
        assert_eq!(CurTreeAlloc(&mut SourceHeap::default(), None, 4), Ok(-1));

        let mut heap = SourceHeap::default();
        let mut tree = CUR_TREE::default();
        assert_eq!(CurTreeAlloc(&mut heap, Some(&mut tree), 4), Ok(0));
        assert!(!tree.tree.is_null());
        assert_eq!((tree.max_len, tree.cur_len, tree.incr_len), (4, 0, 4));
        assert_eq!(heap.slice(tree.tree.as_const()).unwrap(), &[0, 0, 0, 0]);

        heap.slice_mut(tree.tree)
            .unwrap()
            .copy_from_slice(&[1, 2, 3, 4]);
        tree.cur_len = 3;
        let reused = tree.tree;
        assert_eq!(CurTreeAlloc(&mut heap, Some(&mut tree), 2), Ok(0));
        assert_eq!(tree.tree, reused);
        assert_eq!((tree.max_len, tree.cur_len, tree.incr_len), (4, 0, 2));
        assert_eq!(heap.slice(tree.tree.as_const()).unwrap(), &[0, 0, 0, 0]);

        heap.slice_mut(tree.tree)
            .unwrap()
            .copy_from_slice(&[5, 6, 7, 8]);
        tree.cur_len = 2;
        assert_eq!(CurTreeAlloc(&mut heap, Some(&mut tree), -2), Ok(0));
        assert_eq!(tree.tree, reused);
        assert_eq!((tree.max_len, tree.cur_len, tree.incr_len), (4, 0, -2));
        assert_eq!(heap.slice(tree.tree.as_const()).unwrap(), &[0, 0, 0, 0]);

        assert_eq!(CurTreeAlloc(&mut heap, Some(&mut tree), 3), Ok(0));
        assert_ne!(tree.tree, reused);
        assert_eq!((tree.max_len, tree.cur_len, tree.incr_len), (3, 0, 3));
        assert_eq!(heap.slice(tree.tree.as_const()).unwrap(), &[0, 0, 0]);
        assert_eq!(
            heap.slice(reused.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        assert_eq!(
            CurTreeAlloc(&mut heap, Some(&mut tree), 0),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
        assert_eq!((tree.max_len, tree.cur_len, tree.incr_len), (3, 0, 3));

        let mut zero_heap = SourceHeap::default();
        let mut zero_tree = CUR_TREE::default();
        assert_eq!(CurTreeAlloc(&mut zero_heap, Some(&mut zero_tree), 0), Ok(0));
        assert!(!zero_tree.tree.is_null());
        assert_eq!(
            (zero_tree.max_len, zero_tree.cur_len, zero_tree.incr_len),
            (0, 0, 0)
        );
        assert_eq!(
            zero_heap.slice(zero_tree.tree.as_const()).unwrap(),
            &[] as &[u16]
        );

        let mut failure_heap = SourceHeap::default();
        let old = failure_heap.allocate_model_storage(vec![9_u16, 8]).unwrap();
        let mut failure_tree = CUR_TREE {
            tree: old,
            max_len: 2,
            cur_len: 1,
            incr_len: 2,
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            CurTreeAlloc(&mut failure_heap, Some(&mut failure_tree), 3),
            Ok(-1)
        );
        assert_eq!(failure_tree, CUR_TREE::default());
        assert_eq!(
            failure_heap.slice(old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__ichimap1__curtreefree__line_871() {
        let mut heap = SourceHeap::default();
        let tree_pointer = heap.allocate_model_storage(vec![1_u16, 2, 3]).unwrap();
        let mut tree = CUR_TREE {
            tree: tree_pointer,
            max_len: 3,
            cur_len: 2,
            incr_len: 3,
        };
        assert_eq!(CurTreeFree(&mut heap, Some(&mut tree)), Ok(()));
        assert_eq!(tree, CUR_TREE::default());
        assert_eq!(
            heap.slice(tree_pointer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        assert_eq!(CurTreeFree(&mut heap, None), Ok(()));
        let mut no_tree = CUR_TREE::default();
        assert_eq!(CurTreeFree(&mut heap, Some(&mut no_tree)), Ok(()));
        assert_eq!(no_tree, CUR_TREE::default());

        let mut broken_heap = SourceHeap::default();
        let mut broken = CUR_TREE {
            tree: SourceConstPointer::<u16>::null().as_mut(),
            max_len: 9,
            cur_len: 8,
            incr_len: 7,
        };
        assert_eq!(CurTreeFree(&mut broken_heap, Some(&mut broken)), Ok(()));
        assert_eq!(broken, CUR_TREE::default());
    }

    fn stereo_dble(at_num1: u16, at_num2: u16, parity: u8) -> AT_STEREO_DBLE {
        AT_STEREO_DBLE {
            at_num1,
            at_num2,
            parity,
        }
    }

    fn stereo_carb(at_num: u16, parity: u8) -> AT_STEREO_CARB {
        AT_STEREO_CARB { at_num, parity }
    }

    #[test]
    fn source_port__ichimap1__comparelinctstereocarb__line_223() {
        let mut heap = SourceHeap::default();
        let first = heap
            .allocate_model_storage(vec![stereo_carb(1, 2), stereo_carb(3, 4)])
            .unwrap();
        let equal = heap
            .allocate_model_storage(vec![stereo_carb(1, 2), stereo_carb(3, 4)])
            .unwrap();
        let at_num = heap
            .allocate_model_storage(vec![stereo_carb(1, 2), stereo_carb(8, 0)])
            .unwrap();
        let parity = heap
            .allocate_model_storage(vec![stereo_carb(1, 2), stereo_carb(3, 9)])
            .unwrap();
        let precedence = heap
            .allocate_model_storage(vec![stereo_carb(0, u8::MAX)])
            .unwrap();
        let longer = heap
            .allocate_model_storage(vec![
                stereo_carb(1, 2),
                stereo_carb(3, 4),
                stereo_carb(5, 6),
            ])
            .unwrap();

        assert_eq!(
            CompareLinCtStereoCarb(&heap, first.as_const(), 2, equal.as_const(), 2),
            Ok(0)
        );
        assert_eq!(
            CompareLinCtStereoCarb(&heap, first.as_const(), 2, at_num.as_const(), 2),
            Ok(-5)
        );
        assert_eq!(
            CompareLinCtStereoCarb(&heap, at_num.as_const(), 2, first.as_const(), 2),
            Ok(5)
        );
        assert_eq!(
            CompareLinCtStereoCarb(&heap, first.as_const(), 2, parity.as_const(), 2),
            Ok(-5)
        );
        assert_eq!(
            CompareLinCtStereoCarb(&heap, parity.as_const(), 2, first.as_const(), 2),
            Ok(5)
        );
        assert_eq!(
            CompareLinCtStereoCarb(&heap, first.as_const(), 1, precedence.as_const(), 1),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereoCarb(&heap, first.as_const(), 2, longer.as_const(), 3),
            Ok(-1)
        );
        assert_eq!(
            CompareLinCtStereoCarb(&heap, longer.as_const(), 3, first.as_const(), 2),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereoCarb(
                &heap,
                first.as_const(),
                i32::MIN,
                equal.as_const(),
                i32::MAX
            ),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereoCarb(
                &heap,
                first.as_const(),
                i32::MAX,
                equal.as_const(),
                i32::MIN
            ),
            Ok(-1)
        );

        let null = SourceConstPointer::null();
        assert_eq!(CompareLinCtStereoCarb(&heap, null, 1, null, 1), Ok(0));
        for length in [0, -1, i32::MIN] {
            assert_eq!(
                CompareLinCtStereoCarb(&heap, first.as_const(), length, null, 1),
                Ok(0)
            );
            assert_eq!(
                CompareLinCtStereoCarb(&heap, null, 1, first.as_const(), length),
                Ok(0)
            );
        }
        assert_eq!(
            CompareLinCtStereoCarb(&heap, first.as_const(), 1, null, 0),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereoCarb(&heap, null, 0, first.as_const(), 1),
            Ok(-1)
        );

        let empty = heap
            .allocate_model_storage(Vec::<AT_STEREO_CARB>::new())
            .unwrap();
        assert_eq!(
            CompareLinCtStereoCarb(&heap, empty.as_const(), 0, empty.as_const(), 0),
            Ok(0)
        );
        assert_eq!(
            CompareLinCtStereoCarb(&heap, empty.as_const(), 1, first.as_const(), 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichimap1__comparelinctstereo__line_262() {
        let mut heap = SourceHeap::default();
        let dble_low = heap
            .allocate_model_storage(vec![stereo_dble(1, 2, 3)])
            .unwrap();
        let dble_high = heap
            .allocate_model_storage(vec![stereo_dble(2, 2, 3)])
            .unwrap();
        let dble_equal = heap
            .allocate_model_storage(vec![stereo_dble(1, 2, 3)])
            .unwrap();
        let carb_low = heap
            .allocate_model_storage(vec![stereo_carb(4, 5)])
            .unwrap();
        let carb_high = heap
            .allocate_model_storage(vec![stereo_carb(4, 7)])
            .unwrap();
        let empty_carb = heap
            .allocate_model_storage(Vec::<AT_STEREO_CARB>::new())
            .unwrap();

        assert_eq!(
            CompareLinCtStereo(
                &heap,
                dble_low.as_const(),
                1,
                empty_carb.as_const(),
                1,
                dble_high.as_const(),
                1,
                empty_carb.as_const(),
                1,
            ),
            Ok(-1)
        );
        assert_eq!(
            CompareLinCtStereo(
                &heap,
                dble_high.as_const(),
                1,
                empty_carb.as_const(),
                1,
                dble_low.as_const(),
                1,
                empty_carb.as_const(),
                1,
            ),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereo(
                &heap,
                dble_low.as_const(),
                1,
                carb_low.as_const(),
                1,
                dble_equal.as_const(),
                1,
                carb_high.as_const(),
                1,
            ),
            Ok(-2)
        );
        assert_eq!(
            CompareLinCtStereo(
                &heap,
                dble_low.as_const(),
                1,
                carb_low.as_const(),
                1,
                dble_equal.as_const(),
                1,
                carb_low.as_const(),
                1,
            ),
            Ok(0)
        );
        assert_eq!(
            CompareLinCtStereo(
                &heap,
                SourceConstPointer::null(),
                0,
                empty_carb.as_const(),
                1,
                SourceConstPointer::null(),
                0,
                empty_carb.as_const(),
                1,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichimap1__comparelinctstereodble__line_178() {
        let mut heap = SourceHeap::default();
        let first = heap
            .allocate_model_storage(vec![stereo_dble(1, 2, 3), stereo_dble(4, 5, 6)])
            .unwrap();
        let equal = heap
            .allocate_model_storage(vec![stereo_dble(1, 2, 3), stereo_dble(4, 5, 6)])
            .unwrap();
        let at_num1 = heap
            .allocate_model_storage(vec![stereo_dble(1, 2, 3), stereo_dble(7, 0, 0)])
            .unwrap();
        let at_num2 = heap
            .allocate_model_storage(vec![stereo_dble(1, 2, 3), stereo_dble(4, 8, 0)])
            .unwrap();
        let parity = heap
            .allocate_model_storage(vec![stereo_dble(1, 2, 3), stereo_dble(4, 5, 9)])
            .unwrap();
        let precedence = heap
            .allocate_model_storage(vec![stereo_dble(0, u16::MAX, u8::MAX)])
            .unwrap();
        let longer = heap
            .allocate_model_storage(vec![
                stereo_dble(1, 2, 3),
                stereo_dble(4, 5, 6),
                stereo_dble(7, 8, 9),
            ])
            .unwrap();

        assert_eq!(
            CompareLinCtStereoDble(&heap, first.as_const(), 2, equal.as_const(), 2),
            Ok(0)
        );
        assert_eq!(
            CompareLinCtStereoDble(&heap, first.as_const(), 2, at_num1.as_const(), 2),
            Ok(-3)
        );
        assert_eq!(
            CompareLinCtStereoDble(&heap, at_num1.as_const(), 2, first.as_const(), 2),
            Ok(3)
        );
        assert_eq!(
            CompareLinCtStereoDble(&heap, first.as_const(), 2, at_num2.as_const(), 2),
            Ok(-3)
        );
        assert_eq!(
            CompareLinCtStereoDble(&heap, at_num2.as_const(), 2, first.as_const(), 2),
            Ok(3)
        );
        assert_eq!(
            CompareLinCtStereoDble(&heap, first.as_const(), 2, parity.as_const(), 2),
            Ok(-3)
        );
        assert_eq!(
            CompareLinCtStereoDble(&heap, parity.as_const(), 2, first.as_const(), 2),
            Ok(3)
        );
        assert_eq!(
            CompareLinCtStereoDble(&heap, first.as_const(), 1, precedence.as_const(), 1),
            Ok(1)
        );

        assert_eq!(
            CompareLinCtStereoDble(&heap, first.as_const(), 2, longer.as_const(), 3),
            Ok(-1)
        );
        assert_eq!(
            CompareLinCtStereoDble(&heap, longer.as_const(), 3, first.as_const(), 2),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereoDble(&heap, first.as_const(), 0, equal.as_const(), 0),
            Ok(0)
        );
        assert_eq!(
            CompareLinCtStereoDble(
                &heap,
                first.as_const(),
                i32::MIN,
                equal.as_const(),
                i32::MAX
            ),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereoDble(
                &heap,
                first.as_const(),
                i32::MAX,
                equal.as_const(),
                i32::MIN
            ),
            Ok(-1)
        );

        let null = SourceConstPointer::null();
        assert_eq!(CompareLinCtStereoDble(&heap, null, 1, null, 1), Ok(0));
        for length in [0, -1, i32::MIN] {
            assert_eq!(
                CompareLinCtStereoDble(&heap, first.as_const(), length, null, 1),
                Ok(0)
            );
            assert_eq!(
                CompareLinCtStereoDble(&heap, null, 1, first.as_const(), length),
                Ok(0)
            );
        }
        assert_eq!(
            CompareLinCtStereoDble(&heap, first.as_const(), 1, null, 0),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereoDble(&heap, null, 0, first.as_const(), 1),
            Ok(-1)
        );
        assert_eq!(
            CompareLinCtStereoDble(&heap, first.as_const(), i32::MAX, null, 0),
            Ok(1)
        );
        assert_eq!(
            CompareLinCtStereoDble(&heap, null, 0, first.as_const(), i32::MAX),
            Ok(-1)
        );

        let empty = heap
            .allocate_model_storage(Vec::<AT_STEREO_DBLE>::new())
            .unwrap();
        assert_eq!(
            CompareLinCtStereoDble(&heap, empty.as_const(), 0, empty.as_const(), 0),
            Ok(0)
        );
        assert_eq!(
            CompareLinCtStereoDble(&heap, empty.as_const(), 1, first.as_const(), 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichimap1__curtreesetpos__line_1023() {
        let mut tree = CUR_TREE {
            max_len: 17,
            cur_len: 18,
            incr_len: 19,
            ..CUR_TREE::default()
        };
        assert_eq!(CurTreeSetPos(Some(&mut tree), 0), 0);
        assert_eq!((tree.max_len, tree.cur_len, tree.incr_len), (17, 0, 19));

        assert_eq!(CurTreeSetPos(Some(&mut tree), -1), 0);
        assert_eq!(tree.cur_len, -1);
        assert_eq!(CurTreeSetPos(Some(&mut tree), i32::MIN), 0);
        assert_eq!(tree.cur_len, i32::MIN);
        assert_eq!(CurTreeSetPos(Some(&mut tree), i32::MAX), 0);
        assert_eq!(tree.cur_len, i32::MAX);

        assert_eq!(CurTreeSetPos(None, 123), -1);
    }
}
